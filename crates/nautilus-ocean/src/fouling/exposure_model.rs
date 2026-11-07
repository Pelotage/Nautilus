//! Exposure model: one hull exposure clock, coating roughness drift and propeller
//! efficiency drift.
//!
//! The hull clock `E` advances by `T(sst)·S(salinity)·gate(segment)·Δ` per
//! segment; the hull fouling resistance fraction is `min(hull_cap, slope·E)`.
//! This model has a single clock, so there is no freshwater cull: freshwater only
//! pauses the clock (salinity factor 0). The age-clock model handles the cull.

use nautilus_units::{DAYS_PER_YEAR, Days, Micrometres};

use super::environment::Segment;
use super::events::MaintenanceEvent;
use super::params::ExposureParams;
use super::params::{check_fraction, check_non_negative, check_responses};
use super::power::{PowerPenalty, coating_resistance_fraction};
use super::{FoulingError, FoulingModel, OperatingPoint};

/// Exposure model fouling state. Starts freshly docked.
#[derive(Debug, Clone, PartialEq)]
pub struct ExposureModel {
    params: ExposureParams,
    exposure: f64,
    ahr: f64,
    prop: f64,
    /// True between a dry-docking and the next segment with time in it.
    just_docked: bool,
}

impl ExposureModel {
    /// Builds a freshly docked model. Rejects non-finite or negative parameters,
    /// fractions outside [0, 1] and `prop_cap ≥ 1` with
    /// [`FoulingError::InvalidParameter`].
    pub fn new(params: ExposureParams) -> Result<Self, FoulingError> {
        check_responses(&[params.temperature], &[params.salinity], &params.gate)?;
        for (v, what) in [
            (params.slope_per_day, "slope_per_day"),
            (params.hull_cap, "hull_cap"),
            (params.ahr_out_docking.value(), "ahr_out_docking"),
            (params.ahr_drift_per_year.value(), "ahr_drift_per_year"),
            (params.prop_drift_per_year, "prop_drift_per_year"),
            (params.prop_cap, "prop_cap"),
            (params.residual_exposure.value(), "residual_exposure"),
            (params.spot_repair_recovery, "spot_repair_recovery"),
        ] {
            check_non_negative(v, what)?;
        }
        check_fraction(params.spot_repair_recovery, "spot_repair_recovery")?;
        if params.prop_cap >= 1.0 {
            return Err(FoulingError::InvalidParameter("prop_cap"));
        }
        Ok(ExposureModel {
            params,
            exposure: 0.0,
            ahr: params.ahr_out_docking.value(),
            prop: 0.0,
            just_docked: true,
        })
    }

    /// Hull exposure clock (equivalent idle days at 25 °C, full salinity).
    pub fn hull_exposure(&self) -> Days {
        Days(self.exposure)
    }
    /// Current coating average hull roughness.
    pub fn ahr(&self) -> Micrometres {
        Micrometres(self.ahr)
    }
    /// Current propeller efficiency loss.
    pub fn prop_eta_loss(&self) -> f64 {
        self.prop
    }
    /// Hull fouling resistance fraction `min(hull_cap, slope·E)`.
    pub fn hull_fouling_fraction(&self) -> f64 {
        let f = self.params.slope_per_day * self.exposure;
        // Not `f64::min`: that would turn a NaN state into the cap and hide it.
        if f > self.params.hull_cap {
            self.params.hull_cap
        } else {
            f
        }
    }
}

impl FoulingModel for ExposureModel {
    fn advance(&mut self, seg: &Segment) {
        let p = &self.params;
        let dt = seg.duration().value();
        if dt > 0.0 {
            self.just_docked = false;
        }
        self.exposure += p.temperature.factor(seg.sst())
            * p.salinity.factor(seg.salinity())
            * p.gate.gate(seg)
            * dt;
        self.ahr += p.ahr_drift_per_year.value() * dt / DAYS_PER_YEAR;
        self.prop = (self.prop + p.prop_drift_per_year * dt / DAYS_PER_YEAR).min(p.prop_cap);
    }

    fn apply(&mut self, event: MaintenanceEvent) {
        let p = &self.params;
        let out = p.ahr_out_docking.value();
        let e_res = p.residual_exposure.value();
        match event {
            MaintenanceEvent::DryDockFullBlast => {
                self.exposure = 0.0;
                self.ahr = out;
                self.prop = 0.0;
                self.just_docked = true;
            }
            MaintenanceEvent::DryDockSpotRepair => {
                self.exposure = 0.0;
                // A repeated docking record with no time in between is the same docking.
                if !self.just_docked {
                    self.ahr -= p.spot_repair_recovery * (self.ahr - out);
                }
                self.prop = 0.0;
                self.just_docked = true;
            }
            MaintenanceEvent::InWaterHullClean { calcareous_remnant } => {
                let rho = if calcareous_remnant.is_nan() {
                    0.0
                } else {
                    calcareous_remnant.clamp(0.0, 1.0)
                };
                self.exposure = self.exposure.min(e_res + rho * (self.exposure - e_res));
            }
            MaintenanceEvent::HullGrooming => self.exposure = self.exposure.min(e_res),
            MaintenanceEvent::PropellerPolish => self.prop = 0.0,
        }
    }

    fn power_penalty(&self, op: &OperatingPoint) -> PowerPenalty {
        PowerPenalty {
            hull_fouling_fraction: self.hull_fouling_fraction(),
            coating_fraction: coating_resistance_fraction(
                op,
                self.params.ahr_out_docking,
                self.ahr(),
            ),
            propeller_eta_loss: self.prop,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fouling::environment::{GateParams, SalinityResponse, TemperatureResponse};
    use crate::fouling::environment::{PortStatus, Segment};
    use crate::fouling::events::MaintenanceEvent;
    use crate::fouling::granville::GranvilleTable;
    use crate::fouling::power::ShipParticulars;
    use crate::fouling::{FoulingError, FoulingModel, OperatingPoint};
    use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }
    fn berth(days: f64, t: f64) -> Segment {
        Segment::new(
            Days(days),
            Knots(0.0),
            Celsius(t),
            Psu(35.0),
            PortStatus::Berth,
            Kilometres(0.0),
        )
        .unwrap()
    }
    fn sail(days: f64) -> Segment {
        Segment::new(
            Days(days),
            Knots(14.0),
            Celsius(25.0),
            Psu(35.0),
            PortStatus::AtSea,
            Kilometres(300.0),
        )
        .unwrap()
    }
    fn m() -> ExposureModel {
        ExposureModel::new(ExposureParams::default()).unwrap()
    }

    /// Exposure model core idea: fouling clock runs while the ship sits in port, practically stops
    /// at 14 kn offshore, and runs at half speed in 15 °C water (Q10 = 2).
    #[test]
    fn idle_accumulates_sailing_does_not() {
        let mut s = m();
        s.advance(&berth(10.0, 25.0));
        assert!(close(s.hull_exposure().value(), 10.0, 1e-9));
        s.advance(&sail(100.0));
        assert!(close(s.hull_exposure().value(), 10.0, 1e-6));
        s.advance(&berth(10.0, 15.0));
        assert!(close(s.hull_exposure().value(), 15.0, 1e-9)); // Q10 = 2 ⇒ half rate at 15 °C
    }
    /// Exposure model penalty grows linearly with exposure (0.0012 per idle day → 19.98% at 166.5 d)
    /// and saturates at the 86% heavy-calcareous cap (Schultz).
    #[test]
    fn observation_linear_then_capped() {
        let mut s = m();
        s.advance(&berth(166.5, 25.0));
        assert!(close(s.hull_fouling_fraction(), 0.1998, 1e-6));
        s.advance(&berth(1.0e4, 25.0));
        assert_eq!(s.hull_fouling_fraction(), 0.86);
    }
    /// Paint roughness (+20 µm/yr) and propeller loss (+2.6%/yr) grow with calendar time
    /// regardless of activity; the propeller loss caps at 32%.
    #[test]
    fn calendar_drifts() {
        let mut s = m();
        s.advance(&sail(365.25));
        assert!(close(s.ahr().value(), 140.0, 1e-9));
        assert!(close(s.prop_eta_loss(), 0.026, 1e-12));
        s.advance(&sail(20.0 * 365.25));
        assert_eq!(s.prop_eta_loss(), 0.32);
    }
    /// What each maintenance event does in exposure model: in-water clean leaves a 13-day residual
    /// plus 20% of the excess, grooming leaves only the residual, dry-dock resets everything.
    #[test]
    fn reset_maps() {
        let mut s = m();
        s.advance(&berth(500.0, 25.0));
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.2,
        });
        assert!(close(s.hull_exposure().value(), 13.0 + 0.2 * 487.0, 1e-9));
        s.apply(MaintenanceEvent::HullGrooming);
        assert!(close(s.hull_exposure().value(), 13.0, 1e-9));
        s.apply(MaintenanceEvent::DryDockFullBlast);
        assert_eq!(
            (
                s.hull_exposure().value(),
                s.ahr().value(),
                s.prop_eta_loss()
            ),
            (0.0, 120.0, 0.0)
        );
    }
    /// Cleaning an already-clean hull changes nothing; docking twice equals docking once;
    /// polishing zeroes the propeller loss. (Review focus 3.)
    #[test]
    fn events_never_increase_fouling_and_are_idempotent() {
        let mut s = m();
        s.advance(&berth(5.0, 25.0));
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.0,
        });
        assert!(close(s.hull_exposure().value(), 5.0, 1e-9)); // already below residual
        s.apply(MaintenanceEvent::DryDockFullBlast);
        let once = (s.hull_exposure(), s.ahr(), s.prop_eta_loss());
        s.apply(MaintenanceEvent::DryDockFullBlast);
        assert_eq!(once, (s.hull_exposure(), s.ahr(), s.prop_eta_loss()));
        s.apply(MaintenanceEvent::PropellerPolish);
        assert_eq!(s.prop_eta_loss(), 0.0);
    }
    /// A spot-repair docking only recovers half of the paint roughness gained since docking.
    #[test]
    fn spot_repair_partial_coating_reset() {
        let mut s = m();
        s.advance(&sail(2.0 * 365.25)); // AHR 160
        s.apply(MaintenanceEvent::DryDockSpotRepair);
        assert!(close(s.ahr().value(), 140.0, 1e-9));
    }
    /// 30 years at 30 °C with no maintenance: penalties stay finite and within their caps.
    /// (Review focus 2.)
    #[test]
    fn thirty_years_finite() {
        let mut s = m();
        for _ in 0..(30 * 365) {
            s.advance(&berth(1.0, 30.0));
        }
        let (ship, table) = (
            ShipParticulars::new(Metres(176.0), 0.5).unwrap(),
            GranvilleTable::new(),
        );
        let p = s.power_penalty(&OperatingPoint::new(&ship, Knots(14.0), &table).unwrap());
        assert!(
            p.power_factor().is_finite()
                && p.hull_fouling_fraction <= 0.86
                && p.propeller_eta_loss <= 0.32
        );
    }
    /// A ship straight out of dry-dock has a power factor of exactly 1 (penalties are measured
    /// from the as-applied hull).
    #[test]
    fn freshly_docked_has_no_penalty() {
        let (ship, table) = (
            ShipParticulars::new(Metres(176.0), 0.5).unwrap(),
            GranvilleTable::new(),
        );
        let p = m().power_penalty(&OperatingPoint::new(&ship, Knots(14.0), &table).unwrap());
        assert_eq!(p.power_factor(), 1.0);
    }
    /// Parameter values that would make the power factor NaN, infinite or negative
    /// (e.g. a 100% propeller-loss cap, zero settlement width, zero Q10) are rejected up front.
    #[test]
    fn rejects_parameters_that_break_the_power_factor() {
        let bad = |p: ExposureParams| {
            matches!(
                ExposureModel::new(p),
                Err(FoulingError::InvalidParameter(_))
            )
        };
        let d = ExposureParams::default();
        assert!(bad(ExposureParams { prop_cap: 1.0, ..d }));
        assert!(bad(ExposureParams {
            gate: GateParams {
                settle_width: Knots(0.0),
                ..d.gate
            },
            ..d
        }));
        assert!(bad(ExposureParams {
            gate: GateParams {
                sailing_fraction: 1.5,
                ..d.gate
            },
            ..d
        }));
        assert!(bad(ExposureParams {
            temperature: TemperatureResponse {
                q10: 0.0,
                ..d.temperature
            },
            ..d
        }));
        assert!(bad(ExposureParams {
            spot_repair_recovery: 2.0,
            ..d
        }));
        let inverted = SalinityResponse {
            zero_below: Psu(20.0),
            half_at: Psu(5.0),
            full_at: Psu(2.0),
        };
        assert!(bad(ExposureParams {
            salinity: inverted,
            ..d
        }));
        assert!(ExposureModel::new(d).is_ok());
    }
    /// A spot-repair docking recorded twice (duplicate record, no time at sea in between)
    /// counts once: AHR 160 → 140, not 140 → 130. A later docking applies again.
    #[test]
    fn duplicate_spot_repair_counts_once() {
        let mut s = m();
        s.advance(&sail(2.0 * 365.25)); // AHR 160
        s.apply(MaintenanceEvent::DryDockSpotRepair);
        s.apply(MaintenanceEvent::DryDockSpotRepair);
        assert!(close(s.ahr().value(), 140.0, 1e-9));
        s.advance(&sail(365.25)); // AHR 160 again
        s.apply(MaintenanceEvent::DryDockSpotRepair);
        assert!(close(s.ahr().value(), 140.0, 1e-9));
    }
    /// Exposure model with growth while sailing (δ = 0.2): 100 days at 14 kn far offshore adds
    /// ≈ 20 exposure days instead of ≈ 0.
    #[test]
    fn sailing_fraction_drives_exposure_clock() {
        let p = ExposureParams {
            gate: GateParams {
                sailing_fraction: 0.2,
                ..GateParams::default()
            },
            ..ExposureParams::default()
        };
        let mut s = ExposureModel::new(p).unwrap();
        s.advance(&sail(100.0));
        assert!(close(s.hull_exposure().value(), 20.0, 1e-6));
    }
    /// Every numeric exposure model parameter rejects NaN and negative values.
    #[test]
    fn every_parameter_rejects_nan_and_negative() {
        let d = ExposureParams::default();
        type Set = fn(&mut ExposureParams, f64);
        let setters: [(&str, Set); 11] = [
            ("slope", |p, v| p.slope_per_day = v),
            ("hull_cap", |p, v| p.hull_cap = v),
            ("ahr_out", |p, v| {
                p.ahr_out_docking = nautilus_units::Micrometres(v)
            }),
            ("ahr_drift", |p, v| {
                p.ahr_drift_per_year = nautilus_units::Micrometres(v)
            }),
            ("prop_drift", |p, v| p.prop_drift_per_year = v),
            ("prop_cap", |p, v| p.prop_cap = v),
            ("residual", |p, v| p.residual_exposure = Days(v)),
            ("recovery", |p, v| p.spot_repair_recovery = v),
            ("floor", |p, v| p.temperature.floor = v),
            ("t_min", |p, v| {
                p.temperature.t_min = Celsius(if v < 0.0 { f64::NAN } else { v })
            }),
            ("near_coast", |p, v| p.gate.near_coast = Kilometres(v)),
        ];
        for (name, set) in setters {
            for v in [f64::NAN, -1.0] {
                let mut p = d;
                set(&mut p, v);
                assert!(
                    matches!(
                        ExposureModel::new(p),
                        Err(FoulingError::InvalidParameter(_))
                    ),
                    "{name} = {v}"
                );
            }
        }
        let far_inside = ExposureParams {
            gate: GateParams {
                near_coast: Kilometres(200.0),
                ..d.gate
            },
            ..d
        };
        assert!(
            ExposureModel::new(far_inside).is_err(),
            "near_coast beyond far_coast"
        );
    }
    /// A corrupted (NaN) exposure state must surface as NaN, not be hidden behind the 86%
    /// cap (`f64::min` silently returns the cap for NaN).
    #[test]
    fn nan_state_is_not_masked_by_cap() {
        let mut s = m();
        s.exposure = f64::NAN;
        assert!(s.hull_fouling_fraction().is_nan());
    }
}
