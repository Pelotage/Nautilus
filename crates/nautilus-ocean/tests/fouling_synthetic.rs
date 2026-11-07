//! Synthetic-data accuracy and invariance tests: the model checked against
//! independently computed values on random (deterministic) forcing.

use nautilus_ocean::fouling::FoulingModel;
use nautilus_ocean::fouling::OperatingPoint;
use nautilus_ocean::fouling::environment::{
    GateParams, PortStatus, SalinityResponse, Segment, TemperatureResponse,
};
use nautilus_ocean::fouling::events::MaintenanceEvent;
use nautilus_ocean::fouling::granville::{GranvilleTable, delta_cf_ratio_direct};
use nautilus_ocean::fouling::power::{ShipParticulars, reynolds};
use nautilus_ocean::fouling::stage1::{Stage1Model, Stage1Params};
use nautilus_ocean::fouling::stage2::{Stage2Model, Stage2Params};
use nautilus_ocean::fouling::uzun::{EQ_CALCAREOUS, EQ_FAST_ASYM, EQ_MAIN_ASYM, EQ_SLIME, EQ_SOFT};
use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Micrometres, Psu};

/// SplitMix64 — deterministic, seedable. Each test fixes its own seed.
struct Rng(u64);
impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn index(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

const N: usize = 200;

fn port(rng: &mut Rng) -> PortStatus {
    [PortStatus::AtSea, PortStatus::Anchorage, PortStatus::Berth][rng.index(3)]
}

fn random_segment(rng: &mut Rng) -> Segment {
    let (d, v, t, s) = (
        rng.uniform(0.05, 20.0),
        rng.uniform(0.0, 20.0),
        rng.uniform(-2.0, 34.0),
        rng.uniform(0.0, 38.0),
    );
    let p = port(rng);
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

/// Random segment with marine salinity (≥ 20 PSU, so no freshwater cull) and,
/// if given, a fixed SST.
fn marine_segment(rng: &mut Rng, sst: Option<f64>) -> Segment {
    let (d, v) = (rng.uniform(0.05, 20.0), rng.uniform(0.0, 20.0));
    let t = sst.unwrap_or_else(|| rng.uniform(-2.0, 34.0));
    let s = rng.uniform(20.0, 38.0);
    let p = port(rng);
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

fn with_duration(seg: &Segment, d: f64) -> Segment {
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

fn rel_close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()) + 1e-12
}

fn stage1() -> Stage1Model {
    Stage1Model::new(Stage1Params::default()).unwrap()
}
fn stage2() -> Stage2Model {
    Stage2Model::new(Stage2Params::default()).unwrap()
}

fn with_op<R>(f: impl FnOnce(&OperatingPoint) -> R) -> R {
    let (ship, table) = (
        ShipParticulars::new(Metres(176.0), 0.5).unwrap(),
        GranvilleTable::new(),
    );
    f(&OperatingPoint {
        ship: &ship,
        speed: Knots(14.0),
        table: &table,
    })
}

fn stage1_rate(seg: &Segment) -> f64 {
    let p = Stage1Params::default();
    p.temperature.factor(seg.sst()) * p.salinity.factor(seg.salinity()) * p.gate.gate(seg)
}

/// Stage 1's exposure after 200 random segments equals the sum of rate × duration computed
/// in the test from the public response functions. Oracle: independent re-computation.
#[test]
fn stage1_exposure_equals_hand_sum() {
    let mut rng = Rng(1);
    let segs: Vec<Segment> = (0..N).map(|_| random_segment(&mut rng)).collect();
    let mut m = stage1();
    segs.iter().for_each(|s| m.advance(s));
    let expected: f64 = segs
        .iter()
        .map(|s| stage1_rate(s) * s.duration().value())
        .sum();
    assert!(
        rel_close(m.hull_exposure().value(), expected, 1e-9),
        "{} vs {expected}",
        m.hull_exposure().value()
    );
}

/// Shuffling the order of segments does not change Stage 1's exposure (it is a sum).
#[test]
fn stage1_exposure_is_order_independent() {
    let mut rng = Rng(2);
    let segs: Vec<Segment> = (0..N).map(|_| random_segment(&mut rng)).collect();
    let mut shuffled = segs.clone();
    for i in (1..shuffled.len()).rev() {
        let j = rng.index(i + 1);
        shuffled.swap(i, j);
    }
    let (mut a, mut b) = (stage1(), stage1());
    segs.iter().for_each(|s| a.advance(s));
    shuffled.iter().for_each(|s| b.advance(s));
    assert!(rel_close(
        a.hull_exposure().value(),
        b.hull_exposure().value(),
        1e-9
    ));
}

/// Cutting a segment into 2, 5 or 17 equal pieces gives the same state as the whole
/// segment, for both stages — proves the closed-form updates are exact, not approximations.
#[test]
fn segment_splitting_invariance_both_stages() {
    let mut rng = Rng(3);
    for _ in 0..N {
        let seg = marine_segment(&mut rng, None);
        let (mut one1, mut one2) = (stage1(), stage2());
        one1.advance(&seg);
        one2.advance(&seg);
        for k in [2usize, 5, 17] {
            let part = with_duration(&seg, seg.duration().value() / k as f64);
            let (mut s1, mut s2) = (stage1(), stage2());
            for _ in 0..k {
                s1.advance(&part);
                s2.advance(&part);
            }
            assert!(
                rel_close(
                    one1.hull_exposure().value(),
                    s1.hull_exposure().value(),
                    1e-9
                ),
                "k={k} {seg:?}"
            );
            for (a, b, what) in [
                (one2.slime_rating(), s2.slime_rating(), "slime"),
                (one2.soft_rating(), s2.soft_rating(), "soft"),
                (
                    one2.calcareous_coverage(),
                    s2.calcareous_coverage(),
                    "coverage",
                ),
                (one2.hull_ks().value(), s2.hull_ks().value(), "hull_ks"),
                (one2.prop_ks().value(), s2.prop_ks().value(), "prop_ks"),
            ] {
                assert!(
                    rel_close(a, b, 1e-9),
                    "{what} k={k}: {a} vs {b} for {seg:?}"
                );
            }
        }
    }
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

/// At any temperature 12–32 °C, d days of Stage 2 growth equals d × Q10^((T−25)/10) days
/// on Uzun's equatorial curve. Oracle: the Q10 law applied to the published curve.
#[test]
fn stage2_time_temperature_equivalence() {
    let mut rng = Rng(4);
    for _ in 0..N {
        let (t, d) = (rng.uniform(12.0, 32.0), rng.uniform(1.0, 800.0));
        let mut m = stage2();
        m.advance(&idle(d, t, 35.0));
        let slime = EQ_SLIME.rating(d * 4.1554_f64.powf((t - 25.0) / 10.0));
        let soft = EQ_SOFT.rating(d * 1.5406_f64.powf((t - 25.0) / 10.0));
        assert!(rel_close(m.slime_rating(), slime, 1e-9), "T={t} d={d}");
        assert!(rel_close(m.soft_rating(), soft, 1e-9), "T={t} d={d}");
    }
}

/// Without maintenance and at constant temperature, the power factor never goes down.
#[test]
fn penalties_monotone_without_events() {
    let mut rng = Rng(5);
    with_op(|op| {
        for _ in 0..20 {
            let sst = rng.uniform(-2.0, 34.0);
            let (mut s1, mut s2) = (stage1(), stage2());
            let (mut prev1, mut prev2) = (1.0_f64, 1.0_f64);
            for i in 0..N {
                let seg = marine_segment(&mut rng, Some(sst));
                s1.advance(&seg);
                s2.advance(&seg);
                let (p1, p2) = (
                    s1.power_penalty(op).power_factor(),
                    s2.power_penalty(op).power_factor(),
                );
                assert!(
                    p1 >= prev1 - 1e-12,
                    "stage1 SST={sst} step {i}: {p1} < {prev1}"
                );
                assert!(
                    p2 >= prev2 - 1e-12,
                    "stage2 SST={sst} step {i}: {p2} < {prev2}"
                );
                (prev1, prev2) = (p1, p2);
            }
        }
    });
}

/// A ship in water 2 °C warmer never ends up with a smoother hull.
#[test]
fn warmer_never_fouls_less() {
    let mut rng = Rng(6);
    for _ in 0..N {
        let t = rng.uniform(12.0, 30.0);
        let (mut cool, mut warm) = (stage2(), stage2());
        for _ in 0..20 {
            let d = rng.uniform(0.05, 60.0);
            cool.advance(&idle(d, t, 35.0));
            warm.advance(&idle(d, t + 2.0, 35.0));
            assert!(warm.hull_ks().value() >= cool.hull_ks().value(), "T={t}");
        }
    }
}

/// After 200 random voyages, every maintenance event (any remnant fraction) leaves the hull
/// and propeller no rougher than before. (Review focus 3.)
#[test]
fn events_never_increase_fouling_random() {
    let mut rng = Rng(7);
    for _ in 0..N {
        let (mut s1, mut s2) = (stage1(), stage2());
        for _ in 0..rng.index(60) + 1 {
            let seg = random_segment(&mut rng);
            s1.advance(&seg);
            s2.advance(&seg);
        }
        let rho = rng.uniform(0.0, 1.0);
        for event in [
            MaintenanceEvent::DryDockFullBlast,
            MaintenanceEvent::DryDockSpotRepair,
            MaintenanceEvent::InWaterHullClean {
                calcareous_remnant: rho,
            },
            MaintenanceEvent::HullGrooming,
            MaintenanceEvent::PropellerPolish,
        ] {
            let (mut a, mut b) = (s1.clone(), s2.clone());
            a.apply(event);
            b.apply(event);
            assert!(
                a.hull_fouling_fraction() <= s1.hull_fouling_fraction(),
                "{event:?}"
            );
            assert!(b.hull_ks().value() <= s2.hull_ks().value(), "{event:?}");
            assert!(b.prop_ks().value() <= s2.prop_ks().value(), "{event:?}");
        }
    }
}

/// 200 random ships (80–400 m, 8–25 kn, 30 µm–10 mm): lookup table within 0.6 percentage
/// points of the exact solver, and the exact solver never says rougher is faster.
#[test]
fn granville_table_accuracy_random() {
    let mut rng = Rng(8);
    let table = GranvilleTable::new();
    for _ in 0..N {
        let (l, v, ks) = (
            rng.uniform(80.0, 400.0),
            rng.uniform(8.0, 25.0),
            rng.uniform(30.0, 10000.0),
        );
        let re = reynolds(Metres(l), Knots(v));
        let rel = Micrometres(ks).to_metres() / l;
        let direct = delta_cf_ratio_direct(re, rel);
        let tab = table.delta_cf_ratio(re, rel);
        assert!(
            (100.0 * (tab - direct)).abs() <= 0.6,
            "L={l} V={v} ks={ks}: table {tab} direct {direct}"
        );
        let rougher = delta_cf_ratio_direct(re, rel * 1.1);
        assert!(
            rougher >= direct - 1e-12,
            "not monotone at L={l} V={v} ks={ks}"
        );
    }
}

/// A four-leg voyage (hot berth → Panama freshwater → ocean crossing → anchorage) whose
/// expected ratings and coverage are worked out step by step in the test.
#[test]
fn hand_computed_voyage_stage2() {
    // Voyage: 120 d berth 28 °C 35 PSU → 0.5 d Panama-like freshwater transit
    // (4 kn, 0.3 PSU, 28 °C) → 30 d sailing 16 kn 300 km offshore 26 °C 35 PSU →
    // 60 d anchorage 22 °C 33 PSU.
    let s1 = idle(120.0, 28.0, 35.0);
    let s2 = Segment::new(
        Days(0.5),
        Knots(4.0),
        Celsius(28.0),
        Psu(0.3),
        PortStatus::AtSea,
        Kilometres(0.0),
    )
    .unwrap();
    let s3 = Segment::new(
        Days(30.0),
        Knots(16.0),
        Celsius(26.0),
        Psu(35.0),
        PortStatus::AtSea,
        Kilometres(300.0),
    )
    .unwrap();
    let s4 = Segment::new(
        Days(60.0),
        Knots(0.0),
        Celsius(22.0),
        Psu(33.0),
        PortStatus::Anchorage,
        Kilometres(5.0),
    )
    .unwrap();
    let mut m = stage2();
    for s in [&s1, &s2, &s3, &s4] {
        m.advance(s);
    }

    // Clock arithmetic, written out. Q10 factors relative to 25 °C (all T within [12, 32]).
    let q = |q10: f64, t: f64| q10.powf((t - 25.0) / 10.0);
    let (q_sl, q_so, q_c) = (4.1554, 1.5406, 2.0417);
    let gate = GateParams::default();
    // Segment 1: berth ⇒ gate 1; 35 PSU ⇒ salinity factors 1.
    let (mut a_sl, mut a_so, mut a_c) = (
        120.0 * q(q_sl, 28.0),
        120.0 * q(q_so, 28.0),
        120.0 * q(q_c, 28.0),
    );
    let (mut w, mut wt) = (a_c, a_c * 28.0);
    // Segment 2: 0.3 PSU < 2 ⇒ both salinity factors 0 (no growth); 12 h ≥ 3 h ⇒ cull
    // halves slime and soft ratings; calcareous untouched.
    assert_eq!(SalinityResponse::STENOHALINE.factor(Psu(0.3)), 0.0);
    assert_eq!(SalinityResponse::EURYHALINE.factor(Psu(0.3)), 0.0);
    a_sl = EQ_SLIME
        .age_for_rating(0.5 * EQ_SLIME.rating(a_sl))
        .min(a_sl);
    a_so = EQ_SOFT.age_for_rating(0.5 * EQ_SOFT.rating(a_so)).min(a_so);
    // Segment 3: offshore sailing, tiny gate; 35 PSU ⇒ salinity factors 1.
    let g3 = gate.gate(&s3);
    a_sl += 30.0 * q(q_sl, 26.0) * g3;
    a_so += 30.0 * q(q_so, 26.0) * g3;
    let wc3 = 30.0 * q(q_c, 26.0) * g3;
    a_c += wc3;
    w += wc3;
    wt += wc3 * 26.0;
    // Segment 4: anchorage, 0 kn ⇒ gate 1; 33 PSU ⇒ salinity factors 1.
    a_sl += 60.0 * q(q_sl, 22.0);
    a_so += 60.0 * q(q_so, 22.0);
    let wc4 = 60.0 * q(q_c, 22.0);
    a_c += wc4;
    w += wc4;
    wt += wc4 * 22.0;
    // Asymptotes from mean SST, linear between 17 °C (50, 0) and 25 °C (equatorial).
    let f = ((wt / w - 17.0) / 8.0).clamp(0.0, 1.0);
    let (main, fast) = (50.0 + f * (EQ_MAIN_ASYM - 50.0), f * EQ_FAST_ASYM);

    assert!(
        rel_close(m.slime_rating(), EQ_SLIME.rating(a_sl), 1e-9),
        "{} vs {}",
        m.slime_rating(),
        EQ_SLIME.rating(a_sl)
    );
    assert!(
        rel_close(m.soft_rating(), EQ_SOFT.rating(a_so), 1e-9),
        "{} vs {}",
        m.soft_rating(),
        EQ_SOFT.rating(a_so)
    );
    let sc = EQ_CALCAREOUS.coverage(a_c, main, fast);
    assert!(
        rel_close(m.calcareous_coverage(), sc, 1e-9),
        "{} vs {sc}",
        m.calcareous_coverage()
    );
    // Sanity: the default temperature response matches the raw Q10 law in range.
    let t = TemperatureResponse {
        q10: q_sl,
        ..TemperatureResponse::default()
    };
    assert!(rel_close(t.factor(Celsius(22.0)), q(q_sl, 22.0), 1e-12));
}

/// 200 segments of extreme but legal input (−2/34 °C, 0/40 PSU, up to 10 years long): no
/// penalty is ever NaN or infinite. (Review focus 1 and 2.)
#[test]
fn no_nan_under_extremes() {
    let mut rng = Rng(9);
    with_op(|op| {
        let (mut s1, mut s2) = (stage1(), stage2());
        for _ in 0..N {
            let t = if rng.index(2) == 0 { -2.0 } else { 34.0 };
            let s = if rng.index(2) == 0 { 0.0 } else { 40.0 };
            let d = rng.uniform(0.0, 3650.0);
            let v = rng.uniform(0.0, 20.0);
            let p = port(&mut rng);
            let seg = Segment::new(
                Days(d),
                Knots(v),
                Celsius(t),
                Psu(s),
                p,
                Kilometres(rng.uniform(0.0, 400.0)),
            )
            .unwrap();
            s1.advance(&seg);
            s2.advance(&seg);
            for pen in [s1.power_penalty(op), s2.power_penalty(op)] {
                assert!(
                    pen.hull_fouling_fraction.is_finite()
                        && pen.coating_fraction.is_finite()
                        && pen.propeller_eta_loss.is_finite()
                        && pen.power_factor().is_finite(),
                    "{pen:?}"
                );
            }
        }
    });
}
