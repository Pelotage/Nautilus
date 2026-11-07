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

use std::fmt;

use nautilus_units::Knots;

pub mod environment;
pub mod events;
pub mod granville;
pub mod power;
pub mod stage1;
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
    /// Ship particulars.
    pub ship: &'a power::ShipParticulars,
    /// Speed through water.
    pub speed: Knots,
    /// Granville lookup table.
    pub table: &'a granville::GranvilleTable,
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
