//! Parameters of the fouling models, their defaults (with sources) and the shared
//! validation helpers.

use nautilus_units::{Celsius, Days, Micrometres, Psu};

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
