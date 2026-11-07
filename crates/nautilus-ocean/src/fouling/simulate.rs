//! Simulation driver: runs a model along a timeline of segments and events,
//! with optional state-triggered cleaning policies.

use nautilus_units::Days;

use super::environment::Segment;
use super::events::MaintenanceEvent;
use super::power::PowerPenalty;
use super::{FoulingModel, OperatingPoint};

/// One item of a ship's timeline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TimelineItem {
    /// A stretch of constant forcing.
    Segment(Segment),
    /// A scheduled maintenance event (e.g. a known dry-docking).
    Event(MaintenanceEvent),
}

/// In-water cleaning policy applied by [`simulate`] on top of the timeline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CleaningPolicy {
    /// No policy cleaning.
    None,
    /// Clean when the hull fouling fraction first reaches `first_trigger`, then
    /// every `interval` until the next dry-docking.
    Triggered {
        /// Hull fouling resistance fraction that triggers the first clean.
        first_trigger: f64,
        /// Interval between subsequent cleans.
        interval: Days,
    },
}

impl CleaningPolicy {
    /// GloFouling proactive policy: first clean at 20%, then every 6 months.
    pub const PROACTIVE: CleaningPolicy = CleaningPolicy::Triggered {
        first_trigger: 0.20,
        interval: Days(182.5),
    };
    /// GloFouling reactive policy: first clean at 40%, then every 12 months.
    pub const REACTIVE: CleaningPolicy = CleaningPolicy::Triggered {
        first_trigger: 0.40,
        interval: Days(365.0),
    };
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
                if let CleaningPolicy::Triggered {
                    first_trigger,
                    interval,
                } = policy
                {
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
    use crate::fouling::stage1::{Stage1Model, Stage1Params};
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
        let op = OperatingPoint {
            ship: &ship,
            speed: Knots(14.0),
            table: &table,
        };
        simulate(
            &mut Stage1Model::new(Stage1Params::default()).unwrap(),
            timeline,
            policy,
            &op,
        )
    }

    #[test]
    fn mean_of_empty_is_one() {
        assert_eq!(mean_power_factor(&[]), 1.0);
    }
    #[test]
    fn points_follow_timeline_with_cumulative_elapsed() {
        let traj = run(
            &[
                berth(10.0),
                TimelineItem::Event(MaintenanceEvent::PropellerPolish),
                sail(5.0),
            ],
            CleaningPolicy::None,
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
    #[test]
    fn mean_is_duration_weighted() {
        let traj = run(&[berth(100.0), sail(100.0)], CleaningPolicy::None);
        let expected = (traj[0].penalty.power_factor() * 100.0
            + traj[1].penalty.power_factor() * 100.0)
            / 200.0;
        assert!((mean_power_factor(&traj) - expected).abs() < 1e-12);
    }
}
