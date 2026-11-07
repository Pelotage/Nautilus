//! Stage 1: one hull exposure clock, coating roughness drift and propeller
//! efficiency drift.
//!
//! The hull clock `E` advances by `T(sst)·S(salinity)·gate(segment)·Δ` per
//! segment; the hull fouling resistance fraction is `min(hull_cap, slope·E)`.
//! Stage 1 has a single clock, so there is no freshwater cull: freshwater only
//! pauses the clock (salinity factor 0). Stage 2 models the cull.

use nautilus_units::{DAYS_PER_YEAR, Days, Micrometres};

use super::environment::{GateParams, SalinityResponse, Segment, TemperatureResponse};
use super::events::MaintenanceEvent;
use super::power::{PowerPenalty, coating_resistance_fraction};
use super::{FoulingError, FoulingModel, OperatingPoint};

/// Stage 1 parameters. See [`Stage1Params::default`] for sources.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stage1Params {
    /// Temperature response of the hull clock (Q10 2).
    pub temperature: TemperatureResponse,
    /// Salinity response of the hull clock (euryhaline).
    pub salinity: SalinityResponse,
    /// Settlement / propagule gate.
    pub gate: GateParams,
    /// Hull fouling resistance fraction per exposure day. 0.0012 (±50%), derived
    /// from GloFouling's equatorial profile. Assumption — calibrate.
    pub slope_per_day: f64,
    /// Cap on the hull fouling resistance fraction (0.86, Schultz heavy
    /// calcareous shaft-power increase).
    pub hull_cap: f64,
    /// Coating average hull roughness out of dry-dock (120 µm).
    pub ahr_out_docking: Micrometres,
    /// Coating AHR growth per year (20 µm; range 3–40). Assumption — calibrate.
    pub ahr_drift_per_year: Micrometres,
    /// Propeller efficiency loss per year (0.026, van Ballegooijen via GloFouling).
    pub prop_drift_per_year: f64,
    /// Cap on propeller efficiency loss (0.32).
    pub prop_cap: f64,
    /// Exposure left after an in-water clean or grooming (13 days ≈ +1.6% at the
    /// 40 µm GloFouling post-cleaning floor, tanker). Assumption — calibrate.
    pub residual_exposure: Days,
    /// Fraction of excess coating roughness recovered by a spot repair (0.5).
    /// Assumption — calibrate.
    pub spot_repair_recovery: f64,
}

impl Default for Stage1Params {
    fn default() -> Self {
        Stage1Params {
            temperature: TemperatureResponse::default(),
            salinity: SalinityResponse::EURYHALINE,
            gate: GateParams::default(),
            slope_per_day: 0.0012,
            hull_cap: 0.86,
            ahr_out_docking: Micrometres(120.0),
            ahr_drift_per_year: Micrometres(20.0),
            prop_drift_per_year: 0.026,
            prop_cap: 0.32,
            residual_exposure: Days(13.0),
            spot_repair_recovery: 0.5,
        }
    }
}

/// Rejects non-finite or negative values.
pub(crate) fn check_non_negative(v: f64, what: &'static str) -> Result<(), FoulingError> {
    if v.is_finite() && v >= 0.0 {
        Ok(())
    } else {
        Err(FoulingError::InvalidParameter(what))
    }
}

/// Rejects values outside [0, 1].
pub(crate) fn check_fraction(v: f64, what: &'static str) -> Result<(), FoulingError> {
    if (0.0..=1.0).contains(&v) {
        Ok(())
    } else {
        Err(FoulingError::InvalidParameter(what))
    }
}

/// Rejects non-finite or non-positive values.
fn check_positive(v: f64, what: &'static str) -> Result<(), FoulingError> {
    if v.is_finite() && v > 0.0 {
        Ok(())
    } else {
        Err(FoulingError::InvalidParameter(what))
    }
}

/// Validates the shared response parameters: all finite; Q10 and settlement
/// width positive; salinity breakpoints and coast distances non-negative and
/// ordered; gate fractions in [0, 1].
pub(crate) fn check_responses(
    temps: &[TemperatureResponse],
    salinities: &[SalinityResponse],
    gate: &GateParams,
) -> Result<(), FoulingError> {
    for t in temps {
        check_positive(t.q10, "temperature.q10")?;
        check_non_negative(t.floor, "temperature.floor")?;
        for v in [t.t_ref.value(), t.t_min.value(), t.t_max.value()] {
            if !v.is_finite() {
                return Err(FoulingError::InvalidParameter("temperature limits"));
            }
        }
    }
    for s in salinities {
        let (z, h, f) = (s.zero_below.value(), s.half_at.value(), s.full_at.value());
        for v in [z, h, f] {
            check_non_negative(v, "salinity response")?;
        }
        if !(z <= h && h <= f) {
            return Err(FoulingError::InvalidParameter("salinity response order"));
        }
    }
    check_non_negative(gate.settle_v50.value(), "gate.settle_v50")?;
    check_positive(gate.settle_width.value(), "gate.settle_width")?;
    check_non_negative(gate.near_coast.value(), "gate.near_coast")?;
    check_non_negative(gate.far_coast.value(), "gate.far_coast")?;
    if gate.near_coast.value() > gate.far_coast.value() {
        return Err(FoulingError::InvalidParameter("gate coast distances order"));
    }
    check_fraction(gate.offshore_floor, "gate.offshore_floor")?;
    check_fraction(gate.sailing_fraction, "gate.sailing_fraction")?;
    Ok(())
}

/// Stage 1 fouling state. Starts freshly docked.
#[derive(Debug, Clone, PartialEq)]
pub struct Stage1Model {
    params: Stage1Params,
    exposure: f64,
    ahr: f64,
    prop: f64,
}

impl Stage1Model {
    /// Builds a freshly docked model. Rejects non-finite or negative parameters,
    /// fractions outside [0, 1] and `prop_cap ≥ 1` with
    /// [`FoulingError::InvalidParameter`].
    pub fn new(params: Stage1Params) -> Result<Self, FoulingError> {
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
        Ok(Stage1Model {
            params,
            exposure: 0.0,
            ahr: params.ahr_out_docking.value(),
            prop: 0.0,
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
        (self.params.slope_per_day * self.exposure).min(self.params.hull_cap)
    }
}

impl FoulingModel for Stage1Model {
    fn advance(&mut self, seg: &Segment) {
        let p = &self.params;
        let dt = seg.duration().value();
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
            }
            MaintenanceEvent::DryDockSpotRepair => {
                self.exposure = 0.0;
                self.ahr -= p.spot_repair_recovery * (self.ahr - out);
                self.prop = 0.0;
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
    fn m() -> Stage1Model {
        Stage1Model::new(Stage1Params::default()).unwrap()
    }

    /// Stage 1 core idea: fouling clock runs while the ship sits in port, practically stops
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
    /// Stage 1 penalty grows linearly with exposure (0.0012 per idle day → 19.98% at 166.5 d)
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
    /// What each maintenance event does in Stage 1: in-water clean leaves a 13-day residual
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
        let p = s.power_penalty(&OperatingPoint {
            ship: &ship,
            speed: Knots(14.0),
            table: &table,
        });
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
        let p = m().power_penalty(&OperatingPoint {
            ship: &ship,
            speed: Knots(14.0),
            table: &table,
        });
        assert_eq!(p.power_factor(), 1.0);
    }
    /// Parameter values that would make the power factor NaN, infinite or negative
    /// (e.g. a 100% propeller-loss cap, zero settlement width, zero Q10) are rejected up front.
    #[test]
    fn rejects_parameters_that_break_the_power_factor() {
        let bad =
            |p: Stage1Params| matches!(Stage1Model::new(p), Err(FoulingError::InvalidParameter(_)));
        let d = Stage1Params::default();
        assert!(bad(Stage1Params { prop_cap: 1.0, ..d }));
        assert!(bad(Stage1Params {
            gate: GateParams {
                settle_width: Knots(0.0),
                ..d.gate
            },
            ..d
        }));
        assert!(bad(Stage1Params {
            gate: GateParams {
                sailing_fraction: 1.5,
                ..d.gate
            },
            ..d
        }));
        assert!(bad(Stage1Params {
            temperature: TemperatureResponse {
                q10: 0.0,
                ..d.temperature
            },
            ..d
        }));
        assert!(bad(Stage1Params {
            spot_repair_recovery: 2.0,
            ..d
        }));
        let inverted = SalinityResponse {
            zero_below: Psu(20.0),
            half_at: Psu(5.0),
            full_at: Psu(2.0),
        };
        assert!(bad(Stage1Params {
            salinity: inverted,
            ..d
        }));
        assert!(Stage1Model::new(d).is_ok());
    }
}
