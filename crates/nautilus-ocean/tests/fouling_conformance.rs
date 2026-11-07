//! Scenario-level conformance tests for the fouling model.

use nautilus_ocean::fouling::environment::{PortStatus, Segment};
use nautilus_ocean::fouling::granville::GranvilleTable;
use nautilus_ocean::fouling::power::ShipParticulars;
use nautilus_ocean::fouling::simulate::{
    CleaningPolicy, TimelineItem, mean_power_factor, simulate,
};
use nautilus_ocean::fouling::stage1::{Stage1Model, Stage1Params};
use nautilus_ocean::fouling::stage2::{Stage2Model, Stage2Params};
use nautilus_ocean::fouling::{FoulingModel, OperatingPoint};
use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};

/// 10-day cycles: 3.0411 d idle at berth + 6.9589 d sailing at 14 kn, 300 km offshore, 25 °C, 35 PSU.
/// (GloFouling equatorial virtual profile: 556 idle days in 5 years, spread evenly.)
fn equatorial_profile(years: f64) -> Vec<TimelineItem> {
    let cycles = (years * 36.525).ceil() as usize;
    let idle = Segment::new(
        Days(3.0411),
        Knots(0.0),
        Celsius(25.0),
        Psu(35.0),
        PortStatus::Berth,
        Kilometres(0.0),
    )
    .unwrap();
    let sail = Segment::new(
        Days(6.9589),
        Knots(14.0),
        Celsius(25.0),
        Psu(35.0),
        PortStatus::AtSea,
        Kilometres(300.0),
    )
    .unwrap();
    (0..cycles)
        .flat_map(|_| [TimelineItem::Segment(idle), TimelineItem::Segment(sail)])
        .collect()
}

fn with_op<R>(f: impl FnOnce(&OperatingPoint) -> R) -> R {
    let (ship, table) = (
        ShipParticulars::new(Metres(176.0), 0.5).unwrap(),
        GranvilleTable::new(),
    );
    f(&OperatingPoint::new(&ship, Knots(14.0), &table).unwrap())
}

fn stage1() -> Stage1Model {
    Stage1Model::new(Stage1Params::default()).unwrap()
}
fn stage2() -> Stage2Model {
    Stage2Model::new(Stage2Params::default()).unwrap()
}

/// Stage 1 on GloFouling's equatorial ship profile (556 idle days in 5 years) reaches the
/// published 20% penalty at 1.5 years and 40% at 3 years.
/// Oracle: GloFouling 2022 trigger points (the source Stage 1's slope was derived from).
#[test]
fn stage1_hits_glofouling_trigger_points() {
    let traj = with_op(|op| {
        simulate(
            &mut stage1(),
            &equatorial_profile(3.0),
            CleaningPolicy::NONE,
            op,
        )
    });
    let at = |d: f64| {
        traj.iter()
            .rfind(|p| p.elapsed.value() <= d)
            .unwrap()
            .penalty
            .hull_fouling_fraction
    };
    assert!((at(547.5) - 0.20).abs() < 0.005, "{}", at(547.5));
    assert!((at(1095.0) - 0.40).abs() < 0.01, "{}", at(1095.0));
}

fn policy_means<M: FoulingModel>(make: impl Fn() -> M) -> [f64; 3] {
    let timeline = equatorial_profile(5.0);
    with_op(|op| {
        [
            CleaningPolicy::PROACTIVE,
            CleaningPolicy::REACTIVE,
            CleaningPolicy::NONE,
        ]
        .map(|policy| mean_power_factor(&simulate(&mut make(), &timeline, policy, op)))
    })
}

/// Over 5 years, cleaning proactively beats cleaning reactively, which beats never
/// cleaning — for both model stages. Oracle: the qualitative GloFouling result.
#[test]
fn policies_order_mean_penalty_both_stages() {
    let s1 = policy_means(stage1);
    let s2 = policy_means(stage2);
    eprintln!("stage1 proactive/reactive/none = {s1:?}; stage2 = {s2:?}");
    assert!(s1[0] < s1[1] && s1[1] < s1[2], "stage 1: {s1:?}");
    assert!(s2[0] < s2[1] && s2[1] < s2[2], "stage 2: {s2:?}");
}

/// Running the same 5-year scenario twice gives bit-identical output (no hidden randomness
/// or iteration-order effects).
#[test]
fn deterministic_replay() {
    let timeline = equatorial_profile(5.0);
    let (a, b) = with_op(|op| {
        (
            simulate(&mut stage2(), &timeline, CleaningPolicy::REACTIVE, op),
            simulate(&mut stage2(), &timeline, CleaningPolicy::REACTIVE, op),
        )
    });
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    assert_eq!(a, b);
}
