//! Simulation driver: runs a model along a timeline of segments and events,
//! with optional state-triggered cleaning policies.

use nautilus_units::Days;

use super::environment::Segment;
use super::events::MaintenanceEvent;
use super::power::PowerPenalty;
use super::{FoulingError, FoulingModel, OperatingPoint};

/// One item of a ship's timeline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TimelineItem {
    /// A stretch of constant forcing.
    Segment(Segment),
    /// A scheduled maintenance event (e.g. a known dry-docking).
    Event(MaintenanceEvent),
}

/// In-water cleaning policy applied by [`simulate`] on top of the timeline: either no
/// policy cleaning, or clean when the hull fouling fraction first reaches a trigger and
/// then every `interval` until the next dry-docking.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CleaningPolicy {
    trigger: Option<(f64, Days)>,
}

impl CleaningPolicy {
    /// No policy cleaning.
    pub const NONE: CleaningPolicy = CleaningPolicy { trigger: None };
    /// GloFouling proactive policy: first clean at 20%, then every 6 months.
    pub const PROACTIVE: CleaningPolicy = CleaningPolicy {
        trigger: Some((0.20, Days(182.5))),
    };
    /// GloFouling reactive policy: first clean at 40%, then every 12 months.
    pub const REACTIVE: CleaningPolicy = CleaningPolicy {
        trigger: Some((0.40, Days(365.0))),
    };

    /// A triggered policy: first clean when the hull fouling resistance fraction reaches
    /// `first_trigger` (finite, ≥ 0), then every `interval` (finite, > 0).
    pub fn triggered(first_trigger: f64, interval: Days) -> Result<Self, FoulingError> {
        if !(first_trigger.is_finite() && first_trigger >= 0.0) {
            return Err(FoulingError::InvalidParameter(
                "cleaning policy first_trigger",
            ));
        }
        if !(interval.value().is_finite() && interval.value() > 0.0) {
            return Err(FoulingError::InvalidParameter("cleaning policy interval"));
        }
        Ok(CleaningPolicy {
            trigger: Some((first_trigger, interval)),
        })
    }
    /// Hull fouling fraction that triggers the first clean, if this policy cleans.
    pub fn first_trigger(&self) -> Option<f64> {
        self.trigger.map(|t| t.0)
    }
    /// Interval between subsequent cleans, if this policy cleans.
    pub fn interval(&self) -> Option<Days> {
        self.trigger.map(|t| t.1)
    }
}

/// The state of the run after one timeline item or policy event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrajectoryPoint {
    /// Cumulative segment time at this point.
    pub elapsed: Days,
    /// Duration of the segment this point closes (0 for event points).
    pub segment_duration: Days,
    /// Power penalty after the segment or event.
    pub penalty: PowerPenalty,
    /// The event applied at this point, if it is an event point.
    pub event: Option<MaintenanceEvent>,
}

/// Runs `model` along `timeline`, applying `policy`, and returns one point per
/// segment (after advancing) and one per event, timeline or policy.
///
/// A policy clean is an `InWaterHullClean { calcareous_remnant: 0.0 }` followed
/// by a `PropellerPolish`, each its own zero-duration point. Any timeline
/// dry-docking resets the policy to wait for its first trigger again.
///
/// ```
/// use nautilus_ocean::fouling::environment::{PortStatus, Segment};
/// use nautilus_ocean::fouling::events::MaintenanceEvent;
/// use nautilus_ocean::fouling::granville::GranvilleTable;
/// use nautilus_ocean::fouling::power::ShipParticulars;
/// use nautilus_ocean::fouling::simulate::{simulate, CleaningPolicy, TimelineItem};
/// use nautilus_ocean::fouling::{exposure_model::{ExposureModel}, params::{ExposureParams}};
/// use nautilus_ocean::fouling::OperatingPoint;
/// use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};
///
/// let berth = Segment::new(Days(200.0), Knots(0.0), Celsius(25.0), Psu(35.0), PortStatus::Berth, Kilometres(0.0))?;
/// let (ship, table) = (ShipParticulars::new(Metres(176.0), 0.5)?, GranvilleTable::new());
/// let op = OperatingPoint::new(&ship, Knots(14.0), &table)?;
/// let mut model = ExposureModel::new(ExposureParams::default())?;
/// let traj = simulate(&mut model, &[TimelineItem::Segment(berth)], CleaningPolicy::PROACTIVE, &op);
/// // 200 idle days → 24% ≥ the 20% trigger, so the policy cleans hull and propeller.
/// let events: Vec<_> = traj.iter().filter_map(|p| p.event).collect();
/// assert_eq!(events, [MaintenanceEvent::InWaterHullClean { calcareous_remnant: 0.0 }, MaintenanceEvent::PropellerPolish]);
/// # Ok::<(), nautilus_ocean::fouling::FoulingError>(())
/// ```
pub fn simulate<M: FoulingModel>(
    model: &mut M,
    timeline: &[TimelineItem],
    policy: CleaningPolicy,
    op: &OperatingPoint,
) -> Vec<TrajectoryPoint> {
    let mut points = Vec::with_capacity(timeline.len());
    let mut elapsed = 0.0;
    let mut last_clean: Option<f64> = None;
    let event_point = |model: &mut M, points: &mut Vec<TrajectoryPoint>, elapsed: f64, event| {
        model.apply(event);
        points.push(TrajectoryPoint {
            elapsed: Days(elapsed),
            segment_duration: Days(0.0),
            penalty: model.power_penalty(op),
            event: Some(event),
        });
    };
    for item in timeline {
        match *item {
            TimelineItem::Event(event) => {
                if matches!(
                    event,
                    MaintenanceEvent::DryDockFullBlast | MaintenanceEvent::DryDockSpotRepair
                ) {
                    last_clean = None;
                }
                event_point(model, &mut points, elapsed, event);
            }
            TimelineItem::Segment(seg) => {
                model.advance(&seg);
                elapsed += seg.duration().value();
                let penalty = model.power_penalty(op);
                points.push(TrajectoryPoint {
                    elapsed: Days(elapsed),
                    segment_duration: seg.duration(),
                    penalty,
                    event: None,
                });
                if let Some((first_trigger, interval)) = policy.trigger {
                    let due = match last_clean {
                        None => penalty.hull_fouling_fraction >= first_trigger,
                        Some(t) => elapsed - t >= interval.value(),
                    };
                    if due {
                        event_point(
                            model,
                            &mut points,
                            elapsed,
                            MaintenanceEvent::InWaterHullClean {
                                calcareous_remnant: 0.0,
                            },
                        );
                        event_point(
                            model,
                            &mut points,
                            elapsed,
                            MaintenanceEvent::PropellerPolish,
                        );
                        last_clean = Some(elapsed);
                    }
                }
            }
        }
    }
    points
}

/// Duration-weighted mean power factor over segment points (event points have
/// zero weight). 1.0 for an empty or zero-duration trajectory.
///
/// Each segment is weighted by the penalty at its *end*, and [`simulate`]
/// checks policy triggers only at segment ends. Both are accurate only when
/// segments are short relative to how fast fouling changes (days, not weeks):
/// split long segments before simulating.
pub fn mean_power_factor(trajectory: &[TrajectoryPoint]) -> f64 {
    let (mut weighted, mut total) = (0.0, 0.0);
    for p in trajectory {
        let w = p.segment_duration.value();
        weighted += w * p.penalty.power_factor();
        total += w;
    }
    if total > 0.0 { weighted / total } else { 1.0 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fouling::OperatingPoint;
    use crate::fouling::environment::{PortStatus, Segment};
    use crate::fouling::events::MaintenanceEvent;
    use crate::fouling::granville::GranvilleTable;
    use crate::fouling::power::ShipParticulars;
    use crate::fouling::{exposure_model::ExposureModel, params::ExposureParams};
    use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};

    fn berth(days: f64) -> TimelineItem {
        TimelineItem::Segment(
            Segment::new(
                Days(days),
                Knots(0.0),
                Celsius(25.0),
                Psu(35.0),
                PortStatus::Berth,
                Kilometres(0.0),
            )
            .unwrap(),
        )
    }
    fn sail(days: f64) -> TimelineItem {
        TimelineItem::Segment(
            Segment::new(
                Days(days),
                Knots(14.0),
                Celsius(25.0),
                Psu(35.0),
                PortStatus::AtSea,
                Kilometres(300.0),
            )
            .unwrap(),
        )
    }
    fn run(timeline: &[TimelineItem], policy: CleaningPolicy) -> Vec<TrajectoryPoint> {
        let (ship, table) = (
            ShipParticulars::new(Metres(176.0), 0.5).unwrap(),
            GranvilleTable::new(),
        );
        let op = OperatingPoint::new(&ship, Knots(14.0), &table).unwrap();
        simulate(
            &mut ExposureModel::new(ExposureParams::default()).unwrap(),
            timeline,
            policy,
            &op,
        )
    }

    /// An empty trajectory has no penalty (mean power factor 1).
    #[test]
    fn mean_of_empty_is_one() {
        assert_eq!(mean_power_factor(&[]), 1.0);
    }
    /// The simulator emits one point per segment and per event, in order, with elapsed time
    /// adding up correctly.
    #[test]
    fn points_follow_timeline_with_cumulative_elapsed() {
        let traj = run(
            &[
                berth(10.0),
                TimelineItem::Event(MaintenanceEvent::PropellerPolish),
                sail(5.0),
            ],
            CleaningPolicy::NONE,
        );
        assert_eq!(traj.len(), 3);
        assert_eq!(
            traj.iter()
                .map(|p| (p.elapsed.value(), p.segment_duration.value()))
                .collect::<Vec<_>>(),
            vec![(10.0, 10.0), (10.0, 0.0), (15.0, 5.0)]
        );
        assert_eq!(traj[1].event, Some(MaintenanceEvent::PropellerPolish));
        assert!(traj[0].event.is_none() && traj[2].event.is_none());
    }
    /// When a cleaning policy triggers it records a hull clean and a propeller polish.
    #[test]
    fn policy_clean_emits_clean_and_polish_points() {
        let traj = run(&[berth(200.0)], CleaningPolicy::PROACTIVE); // 0.24 ≥ 0.20
        let events: Vec<_> = traj.iter().filter_map(|p| p.event).collect();
        assert_eq!(
            events,
            vec![
                MaintenanceEvent::InWaterHullClean {
                    calcareous_remnant: 0.0
                },
                MaintenanceEvent::PropellerPolish
            ]
        );
        assert!(traj.iter().all(|p| p.elapsed.value() == 200.0));
    }
    /// After a dry-dock the policy waits for the first fouling trigger again rather than
    /// cleaning on the old schedule.
    #[test]
    fn dry_dock_resets_policy_to_first_trigger() {
        // Clean at day 200; dock; then 190 sailing days (≥ 182.5 since the clean, but
        // ~no fouling). A reset policy waits for the 20% trigger, so no further clean.
        let traj = run(
            &[
                berth(200.0),
                TimelineItem::Event(MaintenanceEvent::DryDockFullBlast),
                sail(190.0),
            ],
            CleaningPolicy::PROACTIVE,
        );
        let after_dock: Vec<_> = traj
            .iter()
            .skip_while(|p| p.event != Some(MaintenanceEvent::DryDockFullBlast))
            .skip(1)
            .collect();
        assert_eq!(after_dock.len(), 1);
        assert!(after_dock[0].event.is_none());
    }
    /// The mean power factor is weighted by how long each segment lasts.
    #[test]
    fn mean_is_duration_weighted() {
        let traj = run(&[berth(100.0), sail(100.0)], CleaningPolicy::NONE);
        let expected = (traj[0].penalty.power_factor() * 100.0
            + traj[1].penalty.power_factor() * 100.0)
            / 200.0;
        assert!((mean_power_factor(&traj) - expected).abs() < 1e-12);
    }
    /// Boundary of the cleaning interval: a clean is due when *exactly* `interval` days
    /// have passed (≥, not >). Clean at day 200, then 182.5 days of sailing → due at 382.5.
    #[test]
    fn interval_boundary_is_inclusive() {
        let traj = run(&[berth(200.0), sail(182.5)], CleaningPolicy::PROACTIVE);
        let cleans: Vec<f64> = traj
            .iter()
            .filter(|p| matches!(p.event, Some(MaintenanceEvent::InWaterHullClean { .. })))
            .map(|p| p.elapsed.value())
            .collect();
        assert_eq!(cleans, vec![200.0, 382.5]);
    }
    /// Cleaning-policy settings are checked: the interval must be positive and finite and
    /// the trigger finite and non-negative; the GloFouling presets equal their constructor.
    #[test]
    fn policy_validation() {
        assert!(CleaningPolicy::triggered(0.2, Days(0.0)).is_err());
        assert!(CleaningPolicy::triggered(0.2, Days(-5.0)).is_err());
        assert!(CleaningPolicy::triggered(f64::NAN, Days(10.0)).is_err());
        assert!(CleaningPolicy::triggered(-0.1, Days(10.0)).is_err());
        assert!(CleaningPolicy::triggered(0.2, Days(f64::INFINITY)).is_err());
        assert_eq!(
            CleaningPolicy::triggered(0.20, Days(182.5)).unwrap(),
            CleaningPolicy::PROACTIVE
        );
        assert_eq!(
            CleaningPolicy::triggered(0.40, Days(365.0)).unwrap(),
            CleaningPolicy::REACTIVE
        );
    }
}
