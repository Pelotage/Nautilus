//! Hull and propeller biofouling dynamics.
//!
//! The model is a hybrid (impulsive) dynamical system. Between maintenance
//! events the fouling state grows continuously on *exposure clocks*: clocks that
//! advance at a rate set by sea temperature, salinity, ship speed and distance to
//! coast rather than by calendar time. Maintenance events (dry-docking, in-water
//! cleaning, grooming, propeller polishing) are jumps that reset parts of the state,
//! each leaving a residual.
//!
//! Forcing is piecewise constant per [`environment::Segment`], so every update is
//! closed form; no ODE solver is involved and runs are fully deterministic.
//!
//! Spec: `research/reports/Biofouling growth dynamical model.md`.
//!
//! # Example
//!
//! A tanker idles 200 days at a 25 °C berth; the exposure model model turns that into a power
//! penalty, and a dry-docking resets it to the clean hull:
//!
//! ```
//! use nautilus_ocean::fouling::environment::{PortStatus, Segment};
//! use nautilus_ocean::fouling::events::MaintenanceEvent;
//! use nautilus_ocean::fouling::granville::GranvilleTable;
//! use nautilus_ocean::fouling::power::ShipParticulars;
//! use nautilus_ocean::fouling::{exposure_model::{ExposureModel}, params::{ExposureParams}};
//! use nautilus_ocean::fouling::{FoulingModel, OperatingPoint};
//! use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};
//!
//! let berth = Segment::new(Days(200.0), Knots(0.0), Celsius(25.0), Psu(35.0), PortStatus::Berth, Kilometres(0.0))?;
//! let mut model = ExposureModel::new(ExposureParams::default())?;
//! model.advance(&berth);
//! assert!((model.hull_fouling_fraction() - 0.24).abs() < 1e-12); // 0.0012 per idle day
//!
//! let (ship, table) = (ShipParticulars::new(Metres(176.0), 0.5)?, GranvilleTable::new());
//! let op = OperatingPoint::new(&ship, Knots(14.0), &table)?;
//! assert!(model.power_penalty(&op).power_factor() > 1.24);
//!
//! model.apply(MaintenanceEvent::DryDockFullBlast);
//! assert_eq!(model.power_penalty(&op).power_factor(), 1.0);
//! # Ok::<(), nautilus_ocean::fouling::FoulingError>(())
//! ```

use std::fmt;

use nautilus_units::Knots;

pub mod age_clock_model;
pub mod environment;
pub mod events;
pub mod exposure_model;
pub mod granville;
pub(crate) mod integrate;
pub mod params;
pub mod power;
pub mod rates;
pub mod simulate;
pub mod uzun;

/// A fouling model: continuous growth along segments, jumps at events, and an
/// observation that turns the state into a power penalty.
pub trait FoulingModel {
    /// Advances the state through one segment of constant forcing. A
    /// zero-duration segment is a no-op.
    fn advance(&mut self, seg: &environment::Segment);
    /// Applies a maintenance event. Events never increase fouling.
    fn apply(&mut self, event: events::MaintenanceEvent);
    /// Power penalty of the current state at an operating point, relative to
    /// the clean as-applied reference.
    fn power_penalty(&self, op: &OperatingPoint) -> power::PowerPenalty;
}

/// Where and how fast the ship is operating when a penalty is evaluated.
#[derive(Debug, Clone, Copy)]
pub struct OperatingPoint<'a> {
    ship: &'a power::ShipParticulars,
    speed: Knots,
    table: &'a granville::GranvilleTable,
}

impl<'a> OperatingPoint<'a> {
    /// Builds an operating point. `speed` (through water) must be finite and
    /// positive: friction coefficients are undefined for a ship that is not moving.
    pub fn new(
        ship: &'a power::ShipParticulars,
        speed: Knots,
        table: &'a granville::GranvilleTable,
    ) -> Result<Self, FoulingError> {
        if !(speed.value().is_finite() && speed.value() > 0.0) {
            return Err(FoulingError::InvalidParameter("operating speed"));
        }
        Ok(OperatingPoint { ship, speed, table })
    }
    /// Ship particulars.
    pub fn ship(&self) -> &'a power::ShipParticulars {
        self.ship
    }
    /// Speed through water.
    pub fn speed(&self) -> Knots {
        self.speed
    }
    /// Granville lookup table.
    pub fn table(&self) -> &'a granville::GranvilleTable {
        self.table
    }
}

/// Errors raised by the fouling model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FoulingError {
    /// A segment field is non-finite or out of range; the message names it.
    InvalidSegment(&'static str),
    /// A model parameter is non-finite or out of range; the message names it.
    InvalidParameter(&'static str),
    /// The requested coating type has no calibrated parameters.
    UnsupportedCoating,
}

impl fmt::Display for FoulingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FoulingError::InvalidSegment(what) => write!(f, "invalid segment: {what}"),
            FoulingError::InvalidParameter(what) => write!(f, "invalid parameter: {what}"),
            FoulingError::UnsupportedCoating => {
                write!(f, "unsupported coating type (no calibrated parameters)")
            }
        }
    }
}

impl std::error::Error for FoulingError {}
#[cfg(test)]
mod tests {
    use super::*;

    /// Error messages name what was wrong, so a bad input can be traced.
    #[test]
    fn error_messages_name_the_problem() {
        assert_eq!(
            FoulingError::InvalidSegment("salinity").to_string(),
            "invalid segment: salinity"
        );
        assert_eq!(
            FoulingError::InvalidParameter("q10").to_string(),
            "invalid parameter: q10"
        );
        assert!(
            FoulingError::UnsupportedCoating
                .to_string()
                .contains("unsupported coating")
        );
        let e: Box<dyn std::error::Error> = Box::new(FoulingError::UnsupportedCoating);
        assert!(e.source().is_none());
    }
}
