//! Stage 2: Uzun et al. (2019) three fouling classes (slime, soft, calcareous)
//! driven by temperature-scaled *age clocks*, with an Uzun → k_s → Granville →
//! ΔC_F observation chain.
//!
//! Each class keeps an age clock in equivalent days at the 25 °C reference
//! (Uzun's equatorial site). Per segment the clock advances by
//! `T(sst)·S(salinity)·gate·Δ`, so at 25 °C, full salinity and idle the model
//! reproduces the published equatorial curves exactly. The calcareous coverage
//! asymptotes depend on the exposure-weighted mean SST, interpolated between the
//! Mediterranean (17 °C) and equatorial (25 °C) sites.
//!
//! Q10 values are fitted from the two Uzun sites (17 °C, 25 °C). The calcareous
//! Q10 (2.0417) is fitted to the **logistic midpoints** of coverage
//! (393.12 d equatorial → 695.86 d Mediterranean), not the spec's 2.3, which came
//! from the calcareous rating `t0` that this model does not use.
//!
//! Parameters are for self-polishing (SPC) coatings only; foul-release coatings
//! are refused with [`FoulingError::UnsupportedCoating`].

use nautilus_units::{Celsius, DAYS_PER_YEAR, Days, Micrometres, Psu};

use super::environment::{GateParams, SalinityResponse, Segment, TemperatureResponse};
use super::events::MaintenanceEvent;
use super::power::{
    PowerPenalty, coating_resistance_fraction, fouling_resistance_fraction, propeller_eta_loss,
};
use super::stage1::{check_fraction, check_non_negative, check_responses};
use super::uzun::{self, EQ_CALCAREOUS, EQ_FAST_ASYM, EQ_MAIN_ASYM, EQ_SLIME, EQ_SOFT};
use super::{FoulingError, FoulingModel, OperatingPoint};

/// Slime rating after an in-water clean or grooming (0.75 ≡ k_s 40 µm,
/// GloFouling's post-cleaning floor).
const CLEANED_SLIME_RATING: f64 = 0.75;
/// Temperature cap used in the exposure-weighted mean SST.
const MEAN_SST_CAP: f64 = 32.0;

/// Hull coating type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoatingType {
    /// Self-polishing copolymer antifouling (the coating Uzun's panels used).
    SelfPolishing,
    /// Foul-release (silicone/fluoropolymer) coating. Not supported yet.
    FoulRelease,
}

/// Stage 2 parameters. See [`Stage2Params::default`] for values and sources.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stage2Params {
    /// Coating type; only [`CoatingType::SelfPolishing`] is supported.
    pub coating: CoatingType,
    /// Settlement / propagule gate.
    pub gate: GateParams,
    /// Slime clock temperature response (Q10 4.1554, fitted to Uzun's two sites).
    pub slime_temp: TemperatureResponse,
    /// Soft-fouling clock temperature response (Q10 1.5406, fitted to Uzun's two sites).
    pub soft_temp: TemperatureResponse,
    /// Calcareous clock temperature response (Q10 2.0417, fitted to the coverage
    /// logistic midpoints; the spec's 2.3 is not used).
    pub calc_temp: TemperatureResponse,
    /// Salinity response of slime and soft fouling (stenohaline). Assumption — calibrate.
    pub soft_salinity: SalinityResponse,
    /// Salinity response of calcareous fouling (euryhaline). Assumption — calibrate.
    pub hard_salinity: SalinityResponse,
    /// Mean SST at which the cold (Mediterranean) asymptotes apply (17 °C).
    pub cold_asym_temp: Celsius,
    /// Mean SST at which the warm (equatorial) asymptotes apply (25 °C).
    pub warm_asym_temp: Celsius,
    /// Cold-site main coverage asymptote (50%). Placeholder from two sites —
    /// assumption — calibrate.
    pub cold_main_asym: f64,
    /// Cold-site fast coverage asymptote (0%). Placeholder from two sites —
    /// assumption — calibrate.
    pub cold_fast_asym: f64,
    /// Salinity below which a long enough freshwater stay culls soft fouling (2 PSU).
    pub fresh_threshold: Psu,
    /// Minimum duration of a freshwater stay (consecutive segments below
    /// `fresh_threshold`) for its one cull (3 h).
    pub cull_min_duration: Days,
    /// Rating multiplier applied by a freshwater cull (0.5; range 0–0.5).
    /// Assumption — calibrate.
    pub cull_factor: f64,
    /// Coating AHR out of dry-dock (120 µm).
    pub ahr_out_docking: Micrometres,
    /// Coating AHR growth per year (20 µm; range 3–40). Assumption — calibrate.
    pub ahr_drift_per_year: Micrometres,
    /// Fraction of excess coating roughness recovered by a spot repair (0.5).
    /// Assumption — calibrate.
    pub spot_repair_recovery: f64,
    /// Propeller roughness when new or after dry-dock (30 µm). Assumption — calibrate.
    pub prop_clean: Micrometres,
    /// Propeller roughness after polishing (40 µm). Assumption — calibrate.
    pub prop_polished: Micrometres,
    /// Propeller roughness growth per year (100 µm, Rubert scale). Assumption — calibrate.
    pub prop_drift_per_year: Micrometres,
    /// Cap on drift-driven propeller roughness (300 µm). Assumption — calibrate.
    pub prop_rough_cap: Micrometres,
}

impl Default for Stage2Params {
    fn default() -> Self {
        let q = |q10| TemperatureResponse {
            q10,
            ..TemperatureResponse::default()
        };
        Stage2Params {
            coating: CoatingType::SelfPolishing,
            gate: GateParams::default(),
            slime_temp: q(4.1554),
            soft_temp: q(1.5406),
            calc_temp: q(2.0417),
            soft_salinity: SalinityResponse::STENOHALINE,
            hard_salinity: SalinityResponse::EURYHALINE,
            cold_asym_temp: Celsius(17.0),
            warm_asym_temp: Celsius(25.0),
            cold_main_asym: 50.0,
            cold_fast_asym: 0.0,
            fresh_threshold: Psu(2.0),
            cull_min_duration: Days(0.125),
            cull_factor: 0.5,
            ahr_out_docking: Micrometres(120.0),
            ahr_drift_per_year: Micrometres(20.0),
            spot_repair_recovery: 0.5,
            prop_clean: Micrometres(30.0),
            prop_polished: Micrometres(40.0),
            prop_drift_per_year: Micrometres(100.0),
            prop_rough_cap: Micrometres(300.0),
        }
    }
}

/// Stage 2 fouling state. Starts freshly docked.
#[derive(Debug, Clone, PartialEq)]
pub struct Stage2Model {
    params: Stage2Params,
    a_slime: f64,
    a_soft: f64,
    a_calc: f64,
    a_prop_calc: f64,
    calc_w_sum: f64,
    calc_wt_sum: f64,
    ahr: f64,
    prop_rough: f64,
    fresh_stay: f64,
}

impl Stage2Model {
    /// Builds a freshly docked model. Returns
    /// [`FoulingError::UnsupportedCoating`] for foul-release coatings and
    /// [`FoulingError::InvalidParameter`] for non-finite or negative values.
    pub fn new(params: Stage2Params) -> Result<Self, FoulingError> {
        if params.coating != CoatingType::SelfPolishing {
            return Err(FoulingError::UnsupportedCoating);
        }
        check_responses(
            &[params.slime_temp, params.soft_temp, params.calc_temp],
            &[params.soft_salinity, params.hard_salinity],
            &params.gate,
        )?;
        for (v, what) in [
            (params.cold_main_asym, "cold_main_asym"),
            (params.cold_fast_asym, "cold_fast_asym"),
            (params.fresh_threshold.value(), "fresh_threshold"),
            (params.cull_min_duration.value(), "cull_min_duration"),
            (params.cull_factor, "cull_factor"),
            (params.ahr_out_docking.value(), "ahr_out_docking"),
            (params.ahr_drift_per_year.value(), "ahr_drift_per_year"),
            (params.spot_repair_recovery, "spot_repair_recovery"),
            (params.prop_clean.value(), "prop_clean"),
            (params.prop_polished.value(), "prop_polished"),
            (params.prop_drift_per_year.value(), "prop_drift_per_year"),
            (params.prop_rough_cap.value(), "prop_rough_cap"),
        ] {
            check_non_negative(v, what)?;
        }
        check_fraction(params.cull_factor, "cull_factor")?;
        check_fraction(params.spot_repair_recovery, "spot_repair_recovery")?;
        let (cold, warm) = (params.cold_asym_temp.value(), params.warm_asym_temp.value());
        if !(cold.is_finite() && warm.is_finite() && warm > cold) {
            return Err(FoulingError::InvalidParameter("asymptote temperatures"));
        }
        Ok(Stage2Model {
            params,
            a_slime: 0.0,
            a_soft: 0.0,
            a_calc: 0.0,
            a_prop_calc: 0.0,
            calc_w_sum: 0.0,
            calc_wt_sum: 0.0,
            ahr: params.ahr_out_docking.value(),
            prop_rough: params.prop_clean.value(),
            fresh_stay: 0.0,
        })
    }

    /// Slime fouling rating.
    pub fn slime_rating(&self) -> f64 {
        EQ_SLIME.rating(self.a_slime)
    }
    /// Soft-fouling rating.
    pub fn soft_rating(&self) -> f64 {
        EQ_SOFT.rating(self.a_soft)
    }
    /// Calcareous surface coverage on the hull (%).
    pub fn calcareous_coverage(&self) -> f64 {
        self.coverage_at(self.a_calc)
    }
    /// Hull roughness from total rating and calcareous coverage (Uzun scale; a
    /// freshly docked hull reads slightly above 30 µm because of the curve seeds,
    /// and the power penalty is taken relative to that docked value).
    pub fn hull_ks(&self) -> Micrometres {
        uzun::hull_ks(
            self.slime_rating() + self.soft_rating(),
            self.calcareous_coverage(),
        )
    }
    /// Propeller roughness: the larger of drift roughness and the roughness of
    /// calcareous coverage on the blades.
    pub fn prop_ks(&self) -> Micrometres {
        let from_coverage =
            uzun::ks_from_coverage(self.coverage_at(self.a_prop_calc)).map_or(0.0, |k| k.value());
        Micrometres(self.prop_rough.max(from_coverage))
    }
    /// Current coating average hull roughness.
    pub fn ahr(&self) -> Micrometres {
        Micrometres(self.ahr)
    }

    /// Coverage asymptotes `(main, fast)` from the exposure-weighted mean SST,
    /// linear between the cold and warm sites and clamped.
    fn asymptotes(&self) -> (f64, f64) {
        let p = &self.params;
        let mean_t = if self.calc_w_sum > 0.0 {
            self.calc_wt_sum / self.calc_w_sum
        } else {
            p.warm_asym_temp.value()
        };
        let (cold, warm) = (p.cold_asym_temp.value(), p.warm_asym_temp.value());
        let f = ((mean_t - cold) / (warm - cold)).clamp(0.0, 1.0);
        (
            p.cold_main_asym + f * (EQ_MAIN_ASYM - p.cold_main_asym),
            p.cold_fast_asym + f * (EQ_FAST_ASYM - p.cold_fast_asym),
        )
    }

    fn coverage_at(&self, age: f64) -> f64 {
        let (main, fast) = self.asymptotes();
        EQ_CALCAREOUS.coverage(age, main, fast)
    }

    /// Clock age at which the hull coverage would equal `sc`, at or below `age`.
    fn age_for_coverage(&self, sc: f64, age: f64) -> f64 {
        let (main, fast) = self.asymptotes();
        EQ_CALCAREOUS.age_for_coverage(sc, main, fast, age).min(age)
    }

    fn reset_soft_classes(&mut self) {
        self.a_slime = self
            .a_slime
            .min(EQ_SLIME.age_for_rating(CLEANED_SLIME_RATING));
        self.a_soft = 0.0;
    }

    fn dry_dock(&mut self) {
        self.a_slime = 0.0;
        self.a_soft = 0.0;
        self.a_calc = 0.0;
        self.calc_w_sum = 0.0;
        self.calc_wt_sum = 0.0;
        self.prop_rough = self.params.prop_clean.value();
        self.a_prop_calc = 0.0;
    }
}

impl FoulingModel for Stage2Model {
    fn advance(&mut self, seg: &Segment) {
        let p = self.params;
        let dt = seg.duration().value();
        let (t, s) = (seg.sst(), seg.salinity());
        let g = p.gate.gate(seg);
        let soft_s = p.soft_salinity.factor(s);
        self.a_slime += p.slime_temp.factor(t) * soft_s * g * dt;
        self.a_soft += p.soft_temp.factor(t) * soft_s * g * dt;
        let w_c = p.calc_temp.factor(t) * p.hard_salinity.factor(s) * g;
        self.a_calc += w_c * dt;
        self.a_prop_calc += w_c * dt;
        self.calc_w_sum += w_c * dt;
        self.calc_wt_sum += w_c * dt * t.value().min(MEAN_SST_CAP);
        self.ahr += p.ahr_drift_per_year.value() * dt / DAYS_PER_YEAR;
        self.prop_rough = (self.prop_rough + p.prop_drift_per_year.value() * dt / DAYS_PER_YEAR)
            .min(p.prop_rough_cap.value());

        // Freshwater cull: kills soft fouling; dead calcareous shells stay. It fires
        // once per freshwater stay (consecutive segments below the threshold),
        // when the stay first reaches `cull_min_duration`, so the result does not
        // depend on how the stay is split into segments.
        if s.value() < p.fresh_threshold.value() {
            let before = self.fresh_stay;
            self.fresh_stay += dt;
            if before < p.cull_min_duration.value()
                && self.fresh_stay >= p.cull_min_duration.value()
            {
                self.a_slime = EQ_SLIME
                    .age_for_rating(p.cull_factor * self.slime_rating())
                    .min(self.a_slime);
                self.a_soft = EQ_SOFT
                    .age_for_rating(p.cull_factor * self.soft_rating())
                    .min(self.a_soft);
            }
        } else {
            self.fresh_stay = 0.0;
        }
    }

    fn apply(&mut self, event: MaintenanceEvent) {
        let out = self.params.ahr_out_docking.value();
        match event {
            MaintenanceEvent::DryDockFullBlast => {
                self.dry_dock();
                self.ahr = out;
            }
            MaintenanceEvent::DryDockSpotRepair => {
                self.dry_dock();
                self.ahr -= self.params.spot_repair_recovery * (self.ahr - out);
            }
            MaintenanceEvent::InWaterHullClean { calcareous_remnant } => {
                let rho = if calcareous_remnant.is_nan() {
                    0.0
                } else {
                    calcareous_remnant.clamp(0.0, 1.0)
                };
                self.reset_soft_classes();
                self.a_calc = self.age_for_coverage(rho * self.calcareous_coverage(), self.a_calc);
            }
            MaintenanceEvent::HullGrooming => {
                self.reset_soft_classes();
                self.a_calc = 0.0;
            }
            MaintenanceEvent::PropellerPolish => {
                self.prop_rough = self.prop_rough.min(self.params.prop_polished.value());
                self.a_prop_calc = 0.0;
            }
        }
    }

    fn power_penalty(&self, op: &OperatingPoint) -> PowerPenalty {
        PowerPenalty {
            // Relative to the freshly docked hull: Uzun's Gaussians have non-zero
            // seeds, so a docked hull maps to slightly above the 30 µm reference.
            hull_fouling_fraction: (fouling_resistance_fraction(op, self.hull_ks())
                - fouling_resistance_fraction(
                    op,
                    uzun::ks_from_rating(EQ_SLIME.seed() + EQ_SOFT.seed()),
                ))
            .max(0.0),
            coating_fraction: coating_resistance_fraction(
                op,
                self.params.ahr_out_docking,
                self.ahr(),
            ),
            propeller_eta_loss: propeller_eta_loss(self.prop_ks()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fouling::environment::{PortStatus, Segment};
    use crate::fouling::events::MaintenanceEvent;
    use crate::fouling::granville::GranvilleTable;
    use crate::fouling::power::ShipParticulars;
    use crate::fouling::uzun::{EQ_SOFT, ks_from_rating};
    use crate::fouling::{FoulingError, FoulingModel, OperatingPoint};
    use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }
    fn idle(days: f64, t: f64, s: f64) -> Segment {
        Segment::new(
            Days(days),
            Knots(0.0),
            Celsius(t),
            Psu(s),
            PortStatus::Berth,
            Kilometres(0.0),
        )
        .unwrap()
    }
    fn m() -> Stage2Model {
        Stage2Model::new(Stage2Params::default()).unwrap()
    }

    /// Stage 2 is only calibrated for self-polishing paint; foul-release is refused.
    #[test]
    fn foul_release_refused() {
        let p = Stage2Params {
            coating: CoatingType::FoulRelease,
            ..Stage2Params::default()
        };
        assert!(matches!(
            Stage2Model::new(p),
            Err(FoulingError::UnsupportedCoating)
        ));
    }
    /// At 25 °C in port the Stage 2 model reproduces Uzun's published equatorial curves
    /// exactly (slime at 50 d, soft at 200 d, coverage at 393 d).
    #[test]
    fn reproduces_equatorial_curves_at_25c() {
        let mut s = m();
        s.advance(&idle(50.0, 25.0, 35.0));
        assert!(close(s.slime_rating(), 7.38937, 1e-4));
        s.advance(&idle(150.0, 25.0, 35.0));
        assert!(close(s.soft_rating(), 19.2643, 1e-3));
        s.advance(&idle(193.1204, 25.0, 35.0));
        assert!(close(s.calcareous_coverage(), 51.9495, 1e-2));
    }
    /// At 17 °C the Q10 clocks reproduce Uzun's Mediterranean timing: slime plateaus near day
    /// 271.9 and coverage reaches half its 50% plateau at day ~696.
    #[test]
    fn reproduces_mediterranean_delays_at_17c() {
        let mut s = m();
        s.advance(&idle(268.0, 17.0, 35.0));
        assert!(s.slime_rating() < 20.0);
        s.advance(&idle(8.0, 17.0, 35.0)); // 276 d > Med t0 271.9 (±1.5%)
        assert_eq!(s.slime_rating(), 20.0);
        let mut c = m();
        c.advance(&idle(695.864, 17.0, 35.0));
        assert!(close(c.calcareous_coverage(), 25.0, 0.1)); // half of the 50% Med plateau
    }
    /// A 12 h freshwater transit (e.g. Panama's Gatun Lake) halves slime and soft fouling but
    /// leaves barnacles; a 2.4 h transit is too short to kill anything.
    #[test]
    fn freshwater_culls_soft_not_hard() {
        let mut s = m();
        s.advance(&idle(500.0, 25.0, 35.0));
        let (sl, so, sc) = (s.slime_rating(), s.soft_rating(), s.calcareous_coverage());
        s.advance(&idle(0.5, 25.0, 0.3)); // 12 h Gatun-like transit
        assert!(close(s.slime_rating(), 0.5 * sl, 1e-6));
        assert!(close(s.soft_rating(), 0.5 * so, 1e-6));
        assert!(close(s.calcareous_coverage(), sc, 1e-9));
        let mut short = m();
        short.advance(&idle(500.0, 25.0, 35.0));
        short.advance(&idle(0.1, 25.0, 0.3)); // 2.4 h < 3 h ⇒ pause only
        assert!(close(short.soft_rating(), so, 1e-9));
    }
    /// In 10 PSU brackish water soft fouling cannot grow while barnacles grow at the reduced
    /// euryhaline rate.
    #[test]
    fn brackish_freezes_soft_slows_hard() {
        let mut s = m();
        s.advance(&idle(100.0, 25.0, 10.0));
        assert!(close(s.soft_rating(), EQ_SOFT.seed(), 1e-12)); // stenohaline: 0 at 10 PSU
        let mut r = m();
        r.advance(&idle(100.0 * (0.5 + 0.5 * 5.0 / 15.0), 25.0, 35.0)); // euryhaline factor at 10 PSU
        assert!(close(
            s.calcareous_coverage(),
            r.calcareous_coverage(),
            1e-9
        ));
    }
    /// An in-water clean leaves the 40 µm post-cleaning roughness floor (GloFouling), and a
    /// second clean changes nothing.
    #[test]
    fn in_water_clean_resets_to_40um_floor() {
        let mut s = m();
        s.advance(&idle(300.0, 25.0, 35.0)); // SC ≈ 5.6% at 300 d; ρ = 0 removes it
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.0,
        });
        assert!(close(
            s.hull_ks().value(),
            ks_from_rating(0.75 + EQ_SOFT.seed()).value(),
            1e-6
        ));
        let before = s.hull_ks();
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.0,
        });
        assert_eq!(before, s.hull_ks()); // idempotent
    }
    /// An in-water clean that leaves 30% of barnacles leaves 30% of the coverage.
    #[test]
    fn calcareous_remnant_scales_coverage() {
        let mut s = m();
        s.advance(&idle(600.0, 25.0, 35.0));
        let sc = s.calcareous_coverage();
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.3,
        });
        assert!(close(s.calcareous_coverage(), 0.3 * sc, 1e-6));
    }
    /// Propeller roughness grows 100 µm/yr, polishing brings it to 40 µm, docking to 30 µm.
    #[test]
    fn propeller_drift_polish_dock() {
        let mut s = m();
        s.advance(
            &Segment::new(
                Days(365.25),
                Knots(14.0),
                Celsius(25.0),
                Psu(35.0),
                PortStatus::AtSea,
                Kilometres(300.0),
            )
            .unwrap(),
        );
        assert!(close(s.prop_ks().value(), 130.0, 1e-6));
        s.apply(MaintenanceEvent::PropellerPolish);
        assert!(close(s.prop_ks().value(), 40.0, 1e-9));
        s.apply(MaintenanceEvent::DryDockFullBlast);
        assert!(close(s.prop_ks().value(), 30.0, 1e-9));
    }
    /// 30 years at 31 °C: ratings ≤ 70, coverage ≤ 100%, finite power factor. (Review focus 2.)
    #[test]
    fn thirty_years_bounded() {
        let mut s = m();
        for _ in 0..(30 * 365) {
            s.advance(&idle(1.0, 31.0, 35.0));
        }
        assert!(s.slime_rating() + s.soft_rating() <= 70.0);
        assert!(s.calcareous_coverage() <= 100.0);
        let (ship, table) = (
            ShipParticulars::new(Metres(176.0), 0.5).unwrap(),
            GranvilleTable::new(),
        );
        assert!(
            s.power_penalty(&OperatingPoint {
                ship: &ship,
                speed: Knots(14.0),
                table: &table
            })
            .power_factor()
            .is_finite()
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
        let op = OperatingPoint {
            ship: &ship,
            speed: Knots(14.0),
            table: &table,
        };
        assert_eq!(m().power_penalty(&op).power_factor(), 1.0);
        let mut s = m();
        s.advance(&idle(400.0, 25.0, 35.0));
        s.apply(MaintenanceEvent::DryDockFullBlast);
        assert_eq!(s.power_penalty(&op).power_factor(), 1.0);
    }
    /// The freshwater kill depends on how long the stay is, not on how the AIS track happens
    /// to be cut into segments (one 24 h piece, two 12 h, or twelve 2 h give the same result).
    #[test]
    fn freshwater_cull_independent_of_split() {
        let mut base = m();
        base.advance(&idle(500.0, 25.0, 35.0));
        let soft_after = |pieces: usize| {
            let mut s = base.clone();
            for _ in 0..pieces {
                s.advance(&idle(1.0 / pieces as f64, 25.0, 0.3));
            }
            s.soft_rating()
        };
        let one = soft_after(1);
        assert!(close(one, 25.0, 1e-6));
        for k in [2, 12] {
            assert!(close(soft_after(k), one, 1e-9), "k={k}: {}", soft_after(k));
        }
        // A marine segment ends the stay; the next freshwater stay culls again.
        let mut s = base.clone();
        s.advance(&idle(0.5, 25.0, 0.3));
        s.advance(&idle(0.0, 25.0, 35.0));
        s.advance(&idle(0.5, 25.0, 0.3));
        assert!(close(s.soft_rating(), 12.5, 1e-6));
    }
    /// Stage 2 rejects fraction parameters outside [0, 1] (e.g. a cull that would add fouling).
    #[test]
    fn rejects_out_of_range_fractions() {
        let p = Stage2Params {
            cull_factor: 1.5,
            ..Stage2Params::default()
        };
        assert!(matches!(
            Stage2Model::new(p),
            Err(FoulingError::InvalidParameter(_))
        ));
    }
}
