//! Parameters of the fouling models, their defaults (with sources) and the shared
//! validation helpers.

use nautilus_units::{Celsius, Days, Knots, Micrometres, MilligramsPerCubicMetre, Psu};

use super::FoulingError;
use super::environment::{GateParams, SalinityResponse, TemperatureResponse};

/// Exposure model parameters. See [`ExposureParams::default`] for sources.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExposureParams {
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

impl Default for ExposureParams {
    fn default() -> Self {
        ExposureParams {
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

/// Hull coating type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoatingType {
    /// Self-polishing copolymer antifouling (the coating Uzun's panels used).
    SelfPolishing,
    /// Foul-release (silicone/fluoropolymer) coating. Not supported yet.
    FoulRelease,
}

/// Age-clock model parameters. See [`AgeClockParams::default`] for values and sources.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AgeClockParams {
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

impl Default for AgeClockParams {
    fn default() -> Self {
        let q = |q10| TemperatureResponse {
            q10,
            ..TemperatureResponse::default()
        };
        AgeClockParams {
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

/// Thermal window of one guild (a community adapted to warm or to cold water).
/// All values: assumption — calibrate (literature priors, research spec Section on
/// temperature tolerance).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GuildWindow {
    /// Temperature at which [`ClassRates`] `*_ref` values apply (warm 25 °C; cold 9 °C so
    /// both guilds have the same peak rate at their own optimum).
    pub rate_ref_temp: Celsius,
    /// Optimum temperature: rates rise by Q10 up to here.
    pub t_opt: Celsius,
    /// Rates fall linearly from `t_opt` to zero here.
    pub t_max: Celsius,
    /// Below this temperature the guild dies off.
    pub lethal_low: Celsius,
    /// Above this temperature the guild dies off.
    pub lethal_high: Celsius,
    /// Lower end of the breeding (larval settlement) range.
    pub breed_low: Celsius,
    /// Upper end of the breeding range.
    pub breed_high: Celsius,
}

/// Settlement and growth rates of one fouling class at the guild's reference temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassRates {
    /// Settlement rate at the reference temperature (rating or % per day).
    pub settle_ref: f64,
    /// Growth rate at the reference temperature (per day).
    pub growth_ref: f64,
    /// Q10 of settlement.
    pub settle_q10: f64,
    /// Q10 of growth.
    pub growth_q10: f64,
}

/// Speed-driven detachment of a soft class: `rate_per_knot · max(0, v − threshold)` per day.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Detachment {
    /// Detachment rate per day per knot above the threshold.
    pub rate_per_knot: f64,
    /// Speed through water below which nothing detaches.
    pub threshold: Knots,
}

/// Index of the warm guild in [`GuildParams::guilds`].
pub const WARM: usize = 0;
/// Index of the cold guild in [`GuildParams::guilds`].
pub const COLD: usize = 1;
/// Index of slime in [`GuildParams::classes`] (and of [`GuildParams::detachment`]).
pub const SLIME: usize = 0;
/// Index of soft fouling in [`GuildParams::classes`] (and of [`GuildParams::detachment`]).
pub const SOFT: usize = 1;
/// Index of calcareous fouling in [`GuildParams::classes`].
pub const CALCAREOUS: usize = 2;

/// Parameters of the guild model (`guild_model.rs`). See [`GuildParams::default`] for
/// values and their status (prior or fitted).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GuildParams {
    /// Coating type; only self-polishing is supported.
    pub coating: CoatingType,
    /// Settlement / propagule / growth gate.
    pub gate: GateParams,
    /// Salinity response of slime and soft fouling (stenohaline). Assumption — calibrate.
    pub soft_salinity: SalinityResponse,
    /// Salinity response of calcareous fouling (euryhaline). Assumption — calibrate.
    pub hard_salinity: SalinityResponse,
    /// Thermal windows, indexed by [`WARM`] and [`COLD`].
    pub guilds: [GuildWindow; 2],
    /// Class rates, indexed by [`SLIME`], [`SOFT`], [`CALCAREOUS`] (warm-guild reference;
    /// the cold guild uses the same values at its own reference temperature). Fitted
    /// (Layer 1). Growth rates and their Q10s are physically plausible; the *settlement*
    /// Q10s are not (calcareous ≈ 3.6e9, soft ≈ 0.012). Settlement only sets when each curve
    /// takes off, and Uzun's Mediterranean barnacle curve starts from an almost-zero seed.
    /// In practice barnacles barely settle below about 22 °C, and settlement is capped at the
    /// guild optimum (29 °C).
    pub classes: [ClassRates; 3],
    /// Maximum thermal and brackish mortality, per day (0.1: a stranded community largely
    /// dies within a month). Assumption — calibrate.
    pub mortality_max: f64,
    /// Width of the ramp from zero to maximum mortality outside the lethal window (3 °C).
    /// Assumption — calibrate.
    pub mortality_ramp: Celsius,
    /// Width of the ramp at the edges of the breeding range (1 °C). Assumption — calibrate.
    pub breeding_ramp: Celsius,
    /// Salinity below which freshwater mortality applies (2 PSU).
    pub fresh_threshold: Psu,
    /// Freshwater mortality of slime and soft fouling, per day (ln 2 / 0.125: half killed per
    /// 3 h, from 3–4 h tunicate kill times).
    pub fresh_mortality_soft: f64,
    /// Freshwater mortality of live barnacles, per day (ln 2 / 3). Assumption — calibrate.
    pub fresh_mortality_hard: f64,
    /// Food multiplier strength α in `1 + α·chl/(chl + K)` (0: no effect until calibrated).
    pub food_alpha: f64,
    /// Food half-saturation K (1 mg/m³). Assumption — calibrate.
    pub food_half_saturation: MilligramsPerCubicMetre,
    /// Detachment of slime and soft fouling, indexed by [`SLIME`] and [`SOFT`] (rate 0 for
    /// self-polishing paint; thresholds 30 kn and 10 kn from foul-release release speeds).
    pub detachment: [Detachment; 2],
    /// Calcareous capacity at 17 °C, % coverage. Fitted (Layer 1): 49.9 (Uzun plateau 50).
    pub capacity_17: f64,
    /// Calcareous capacity at 25 °C, % coverage. Fitted (Layer 1): 96.8 (Uzun main term 96.9).
    pub capacity_25: f64,
    /// Floor of the calcareous capacity line, % coverage (10). Assumption — calibrate.
    pub capacity_min: f64,
    /// Weathering of dead shells, per day (0.001: half-life about 2 years). Neither Uzun site
    /// has any die-off, so Layer 1 cannot identify it (a free fit ran to 200/day).
    /// Assumption — calibrate.
    pub weathering_per_day: f64,
    /// Slope of the coverage-to-roughness mapping above 50 % coverage, µm per % (80, Uzun
    /// Eq. 18). Kept at Uzun's value: Layer 2 showed the GloFouling gap is not mainly above
    /// 50 % coverage.
    pub high_coverage_slope: f64,
    /// Scale λ on all Uzun-derived fouling roughness above the clean 30 µm reference:
    /// k_s = 30 + λ·(k_s,Uzun − 30). Fitted (Layer 2).
    pub roughness_scale: f64,
    /// Paint roughness added per 100 % calcareous coverage scraped off by an in-water clean
    /// (10 µm). Layer 3 showed realistic damage barely moves the dry-dock / in-water ratio, so
    /// it is a fixed prior. Assumption — calibrate.
    pub clean_damage: Micrometres,
    /// Share of the hull's fouling an in-water clean actually treats (niche areas, parts of the
    /// flat bottom and sea chests are missed). Fitted (Layer 3).
    pub clean_efficiency: f64,
    /// Coating AHR out of dry-dock (120 µm).
    pub ahr_out_docking: Micrometres,
    /// Coating AHR growth per year (20 µm). Assumption — calibrate.
    pub ahr_drift_per_year: Micrometres,
    /// Fraction of excess coating roughness recovered by a spot repair (0.5). Assumption —
    /// calibrate.
    pub spot_repair_recovery: f64,
    /// Propeller roughness when new or after dry-dock (30 µm). Assumption — calibrate.
    pub prop_clean: Micrometres,
    /// Propeller roughness after polishing (40 µm). Assumption — calibrate.
    pub prop_polished: Micrometres,
    /// Propeller roughness growth per year (100 µm). Assumption — calibrate.
    pub prop_drift_per_year: Micrometres,
    /// Cap on drift-driven propeller roughness (300 µm). Assumption — calibrate.
    pub prop_rough_cap: Micrometres,
}

impl Default for GuildParams {
    fn default() -> Self {
        let window = |r, o, m, ll, lh, bl, bh| GuildWindow {
            rate_ref_temp: Celsius(r),
            t_opt: Celsius(o),
            t_max: Celsius(m),
            lethal_low: Celsius(ll),
            lethal_high: Celsius(lh),
            breed_low: Celsius(bl),
            breed_high: Celsius(bh),
        };
        let rates = |s, g, qs, qg| ClassRates {
            settle_ref: s,
            growth_ref: g,
            settle_q10: qs,
            growth_q10: qg,
        };
        GuildParams {
            coating: CoatingType::SelfPolishing,
            gate: GateParams::default(),
            soft_salinity: SalinityResponse::STENOHALINE,
            hard_salinity: SalinityResponse::EURYHALINE,
            guilds: [
                window(25.0, 29.0, 35.0, 10.0, 35.0, 15.0, 35.0),
                window(9.0, 13.0, 24.0, -2.0, 24.0, -2.0, 12.0),
            ],
            // Layer 1 fit: examples/calibrate_fouling.rs layer1, 2026-11-07. Worst gaps to Uzun
            // 2019 at 25 °C / 17 °C: slime 0.081 / 0.065, soft 0.005 / 0.005 rating points,
            // coverage 0.074 / 0.074 percentage points (limits 1, 2, 3).
            classes: [
                rates(
                    1.0025880039562401e-2,
                    5.252675505260153e-2,
                    3.966155775502118e1,
                    3.3036727918417013e0,
                ),
                rates(
                    4.926792468332563e-6,
                    2.7350951847529014e-2,
                    1.2200085392013421e-2,
                    1.944356796762606e0,
                ),
                rates(
                    4.4801957879398244e-7,
                    4.067578624297247e-2,
                    3.550678779933242e9,
                    8.229349560531856e-1,
                ),
            ],
            mortality_max: 0.1,
            mortality_ramp: Celsius(3.0),
            breeding_ramp: Celsius(1.0),
            fresh_threshold: Psu(2.0),
            fresh_mortality_soft: std::f64::consts::LN_2 / 0.125,
            fresh_mortality_hard: std::f64::consts::LN_2 / 3.0,
            food_alpha: 0.0,
            food_half_saturation: MilligramsPerCubicMetre(1.0),
            detachment: [
                Detachment {
                    rate_per_knot: 0.0,
                    threshold: Knots(30.0),
                },
                Detachment {
                    rate_per_knot: 0.0,
                    threshold: Knots(10.0),
                },
            ],
            capacity_17: 4.9930829219449535e1,
            capacity_25: 9.682533388237766e1,
            capacity_min: 10.0,
            weathering_per_day: 0.001,
            high_coverage_slope: 80.0,
            // Layer 2 fit: examples/calibrate_fouling.rs layer2, 2026-11-07. 5-year mean hull
            // penalty 0.330 equatorial (GloFouling 0.30) and 0.172 Mediterranean (0.22); limits ± 0.05.
            roughness_scale: 5.190796405078837e-1,
            clean_damage: Micrometres(10.0),
            // Layer 3 fit: examples/calibrate_fouling.rs layer3, 2026-11-07. Dry-dock / in-water
            // recovery ratio 2.000 (Adland ~2; limits 1.6–2.4).
            clean_efficiency: 6.300333712297965e-1,
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

impl GuildParams {
    /// Checks every parameter: all finite; rates and ramps non-negative (Q10s, ramps and the
    /// food half-saturation positive); fractions in [0, 1]; each guild window ordered
    /// (`t_opt < t_max`, `lethal_low < lethal_high`, `breed_low < breed_high`); capacities
    /// in (0, 100] with `capacity_min ≤ min(capacity_17, capacity_25)`. The coating type is
    /// checked by the model, not here.
    pub fn validate(&self) -> Result<(), FoulingError> {
        check_responses(&[], &[self.soft_salinity, self.hard_salinity], &self.gate)?;
        for w in &self.guilds {
            let t = [
                w.rate_ref_temp,
                w.t_opt,
                w.t_max,
                w.lethal_low,
                w.lethal_high,
                w.breed_low,
                w.breed_high,
            ];
            if t.iter().any(|c| !c.value().is_finite()) {
                return Err(FoulingError::InvalidParameter("guild window temperature"));
            }
            if !(w.t_opt.value() < w.t_max.value()
                && w.lethal_low.value() < w.lethal_high.value()
                && w.breed_low.value() < w.breed_high.value())
            {
                return Err(FoulingError::InvalidParameter("guild window order"));
            }
        }
        for c in &self.classes {
            check_non_negative(c.settle_ref, "settle_ref")?;
            check_non_negative(c.growth_ref, "growth_ref")?;
            check_positive(c.settle_q10, "settle_q10")?;
            check_positive(c.growth_q10, "growth_q10")?;
        }
        for d in &self.detachment {
            check_non_negative(d.rate_per_knot, "detachment rate")?;
            check_non_negative(d.threshold.value(), "detachment threshold")?;
        }
        for (v, what) in [
            (self.mortality_max, "mortality_max"),
            (self.fresh_threshold.value(), "fresh_threshold"),
            (self.fresh_mortality_soft, "fresh_mortality_soft"),
            (self.fresh_mortality_hard, "fresh_mortality_hard"),
            (self.food_alpha, "food_alpha"),
            (self.weathering_per_day, "weathering_per_day"),
            (self.high_coverage_slope, "high_coverage_slope"),
            (self.clean_damage.value(), "clean_damage"),
            (self.ahr_out_docking.value(), "ahr_out_docking"),
            (self.ahr_drift_per_year.value(), "ahr_drift_per_year"),
            (self.prop_clean.value(), "prop_clean"),
            (self.prop_polished.value(), "prop_polished"),
            (self.prop_drift_per_year.value(), "prop_drift_per_year"),
            (self.prop_rough_cap.value(), "prop_rough_cap"),
        ] {
            check_non_negative(v, what)?;
        }
        check_positive(self.mortality_ramp.value(), "mortality_ramp")?;
        check_positive(self.breeding_ramp.value(), "breeding_ramp")?;
        check_positive(self.food_half_saturation.value(), "food_half_saturation")?;
        check_fraction(self.spot_repair_recovery, "spot_repair_recovery")?;
        check_positive(self.roughness_scale, "roughness_scale")?;
        if !(self.clean_efficiency > 0.0 && self.clean_efficiency <= 1.0) {
            return Err(FoulingError::InvalidParameter("clean_efficiency"));
        }
        for (v, what) in [
            (self.capacity_17, "capacity_17"),
            (self.capacity_25, "capacity_25"),
            (self.capacity_min, "capacity_min"),
        ] {
            if !(v > 0.0 && v <= 100.0) {
                return Err(FoulingError::InvalidParameter(what));
            }
        }
        if self.capacity_min > self.capacity_17 || self.capacity_min > self.capacity_25 {
            return Err(FoulingError::InvalidParameter(
                "capacity_min above capacity",
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bad(p: GuildParams) -> bool {
        matches!(p.validate(), Err(FoulingError::InvalidParameter(_)))
    }

    /// The default guild parameters (literature priors plus fitted values) pass validation.
    #[test]
    fn defaults_are_valid() {
        assert_eq!(GuildParams::default().validate(), Ok(()));
    }

    /// Each numeric guild parameter rejects NaN and negative values, so no invalid number
    /// can reach the solver.
    #[test]
    fn every_guild_rule_rejects_bad_values() {
        type Set = fn(&mut GuildParams, f64);
        let setters: [(&str, Set); 20] = [
            ("settle_ref", |p, v| p.classes[0].settle_ref = v),
            ("growth_ref", |p, v| p.classes[1].growth_ref = v),
            ("settle_q10", |p, v| p.classes[2].settle_q10 = v),
            ("growth_q10", |p, v| p.classes[0].growth_q10 = v),
            ("mortality_max", |p, v| p.mortality_max = v),
            ("mortality_ramp", |p, v| p.mortality_ramp = Celsius(v)),
            ("breeding_ramp", |p, v| p.breeding_ramp = Celsius(v)),
            ("fresh_threshold", |p, v| p.fresh_threshold = Psu(v)),
            ("fresh_soft", |p, v| p.fresh_mortality_soft = v),
            ("fresh_hard", |p, v| p.fresh_mortality_hard = v),
            ("food_alpha", |p, v| p.food_alpha = v),
            ("food_half", |p, v| {
                p.food_half_saturation = MilligramsPerCubicMetre(v)
            }),
            ("detach_rate", |p, v| p.detachment[1].rate_per_knot = v),
            ("detach_threshold", |p, v| {
                p.detachment[0].threshold = Knots(v)
            }),
            ("weathering", |p, v| p.weathering_per_day = v),
            ("slope", |p, v| p.high_coverage_slope = v),
            ("clean_damage", |p, v| p.clean_damage = Micrometres(v)),
            ("ahr_out", |p, v| p.ahr_out_docking = Micrometres(v)),
            ("prop_drift", |p, v| p.prop_drift_per_year = Micrometres(v)),
            ("recovery", |p, v| {
                p.spot_repair_recovery = if v < 0.0 { 2.0 } else { v }
            }),
        ];
        for (name, set) in setters {
            for v in [f64::NAN, -1.0] {
                let mut p = GuildParams::default();
                set(&mut p, v);
                assert!(bad(p), "{name} = {v}");
            }
        }
        let mut p = GuildParams::default();
        p.guilds[1].t_opt = Celsius(f64::NAN);
        assert!(bad(p), "NaN window temperature");
    }

    /// Windows whose edges coincide or are reversed are rejected (a zero-width window would
    /// divide by zero in the rate functions).
    #[test]
    fn window_edges_rejected() {
        type Set = fn(&mut GuildWindow);
        let cases: [(&str, Set); 4] = [
            ("t_opt == t_max", |w| w.t_max = w.t_opt),
            ("lethal_low == lethal_high", |w| {
                w.lethal_high = w.lethal_low
            }),
            ("breed_low == breed_high", |w| w.breed_high = w.breed_low),
            ("t_opt > t_max", |w| {
                w.t_opt = Celsius(w.t_max.value() + 1.0)
            }),
        ];
        for (name, set) in cases {
            for g in 0..2 {
                let mut p = GuildParams::default();
                set(&mut p.guilds[g]);
                assert!(bad(p), "guild {g}: {name}");
            }
        }
    }

    /// Calcareous capacities must lie in (0, 100] % and the floor must not exceed either
    /// fitted capacity.
    #[test]
    fn capacity_bounds() {
        let d = GuildParams::default();
        assert!(bad(GuildParams {
            capacity_17: 0.0,
            ..d
        }));
        assert!(bad(GuildParams {
            capacity_25: 101.0,
            ..d
        }));
        assert!(bad(GuildParams {
            capacity_min: d.capacity_17 + 1.0,
            ..d
        }));
        assert!(bad(GuildParams {
            capacity_min: 0.0,
            ..d
        }));
        assert!(bad(GuildParams {
            roughness_scale: 0.0,
            ..d
        }));
        assert!(bad(GuildParams {
            roughness_scale: f64::NAN,
            ..d
        }));
        for eta in [0.0, 1.5, f64::NAN] {
            assert!(
                bad(GuildParams {
                    clean_efficiency: eta,
                    ..d
                }),
                "clean_efficiency {eta}"
            );
        }
    }
}
