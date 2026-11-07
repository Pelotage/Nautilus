//! Helpers shared by the integration tests (deterministic RNG, scenarios, tolerances).
#![allow(dead_code)]

use nautilus_ocean::fouling::OperatingPoint;
use nautilus_ocean::fouling::environment::{PortStatus, Segment};
use nautilus_ocean::fouling::granville::GranvilleTable;
use nautilus_ocean::fouling::power::ShipParticulars;
use nautilus_ocean::fouling::simulate::TimelineItem;
use nautilus_ocean::fouling::stage1::{Stage1Model, Stage1Params};
use nautilus_ocean::fouling::stage2::{Stage2Model, Stage2Params};
use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};

/// SplitMix64 — deterministic, seedable random numbers (no external crates).
pub struct Rng(pub u64);
impl Rng {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn index(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// |a − b| within `tol` relative to the larger magnitude (plus a 1e-12 absolute floor).
pub fn rel_close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()) + 1e-12
}

pub fn stage1() -> Stage1Model {
    Stage1Model::new(Stage1Params::default()).unwrap()
}
pub fn stage2() -> Stage2Model {
    Stage2Model::new(Stage2Params::default()).unwrap()
}

pub fn idle(days: f64, t: f64, s: f64) -> Segment {
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

/// Any legal segment: 0.05–20 d, 0–20 kn, −2–34 °C, 0–38 PSU, any port status, 0–400 km.
pub fn random_segment(rng: &mut Rng) -> Segment {
    let (d, v, t, s) = (
        rng.uniform(0.05, 20.0),
        rng.uniform(0.0, 20.0),
        rng.uniform(-2.0, 34.0),
        rng.uniform(0.0, 38.0),
    );
    let p = [PortStatus::AtSea, PortStatus::Anchorage, PortStatus::Berth][rng.index(3)];
    Segment::new(
        Days(d),
        Knots(v),
        Celsius(t),
        Psu(s),
        p,
        Kilometres(rng.uniform(0.0, 400.0)),
    )
    .unwrap()
}

pub fn with_duration(seg: &Segment, d: f64) -> Segment {
    Segment::new(
        Days(d),
        seg.speed_through_water(),
        seg.sst(),
        seg.salinity(),
        seg.port(),
        seg.distance_to_coast(),
    )
    .unwrap()
}

/// GloFouling-style virtual profile: 10-day cycles of 3.0411 d idle at berth and 6.9589 d
/// sailing at 14 kn 300 km offshore (556 idle days in 5 years), at a constant SST, 35 PSU.
pub fn glofouling_profile(years: f64, sst: f64) -> Vec<TimelineItem> {
    let cycles = (years * 36.525).ceil() as usize;
    let idle = idle(3.0411, sst, 35.0);
    let sail = Segment::new(
        Days(6.9589),
        Knots(14.0),
        Celsius(sst),
        Psu(35.0),
        PortStatus::AtSea,
        Kilometres(300.0),
    )
    .unwrap();
    (0..cycles)
        .flat_map(|_| [TimelineItem::Segment(idle), TimelineItem::Segment(sail)])
        .collect()
}

/// 176 m tanker, 50% frictional resistance share, 14 kn.
pub fn with_tanker<R>(f: impl FnOnce(&OperatingPoint) -> R) -> R {
    let (ship, table) = (
        ShipParticulars::new(Metres(176.0), 0.5).unwrap(),
        GranvilleTable::new(),
    );
    f(&OperatingPoint::new(&ship, Knots(14.0), &table).unwrap())
}
