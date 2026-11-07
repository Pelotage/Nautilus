//! Maintenance events: the jumps of the hybrid system.

/// A maintenance event applied between segments. Every event leaves a
/// residual; none ever increases fouling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaintenanceEvent {
    /// Dry-docking with full blast and recoat: hull, coating and propeller back
    /// to the as-applied state.
    DryDockFullBlast,
    /// Dry-docking with spot repair: fouling removed, coating roughness only
    /// partly recovered.
    DryDockSpotRepair,
    /// In-water hull cleaning. `calcareous_remnant` ∈ [0, 1] is the fraction of
    /// hard (calcareous) fouling left behind; values outside are clamped.
    InWaterHullClean {
        /// Fraction of calcareous fouling remaining after the clean.
        calcareous_remnant: f64,
    },
    /// Light, frequent proactive grooming: removes slime and soft fouling down to
    /// the post-cleaning floor and early calcareous settlers.
    HullGrooming,
    /// Propeller polishing.
    PropellerPolish,
}
