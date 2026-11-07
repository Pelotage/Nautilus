//! Mathematical properties and adversarial inputs.
//!
//! These tests do not compare against measured data. They check that the equations are
//! solved correctly (residuals, identities, inverses), that invariants hold over very
//! long random histories, and that hand-derived schedules come out exactly. Tolerances
//! are near machine precision wherever the maths allows it.

mod common;

use common::*;
use nautilus_ocean::fouling::FoulingModel;
use nautilus_ocean::fouling::OperatingPoint;
use nautilus_ocean::fouling::events::MaintenanceEvent;
use nautilus_ocean::fouling::granville::{GranvilleTable, delta_cf_ratio_direct, schoenherr_cf};
use nautilus_ocean::fouling::power::{ShipParticulars, fouling_resistance_fraction, reynolds};
use nautilus_ocean::fouling::simulate::{CleaningPolicy, simulate};
use nautilus_ocean::fouling::uzun::{EQ_CALCAREOUS, EQ_SLIME, EQ_SOFT};
use nautilus_units::{Knots, Metres, Micrometres};

/// The Schoenherr solver's answer satisfies its own equation, 0.242/√C_F = log10(Re·C_F),
/// to 1e-12 over 200 random Reynolds numbers from model scale (1e5) to the largest
/// ships (1e10).
#[test]
fn schoenherr_solves_its_equation() {
    let mut rng = Rng(101);
    for _ in 0..200 {
        let re = 10f64.powf(rng.uniform(5.0, 10.0));
        let cf = schoenherr_cf(re);
        let residual = 0.242 / cf.sqrt() - (re * cf).log10();
        assert!(residual.abs() < 1e-12, "Re={re:e}: residual {residual:e}");
    }
}

/// Identity of the Granville method: with no roughness the roughness function is zero,
/// so the rough-plate friction must reduce *exactly* to the smooth Schoenherr line.
#[test]
fn granville_reduces_to_smooth_line_without_roughness() {
    let mut rng = Rng(102);
    for _ in 0..50 {
        let re = 10f64.powf(rng.uniform(7.0, 10.0));
        let d = delta_cf_ratio_direct(re, 1e-14);
        assert!(d.abs() < 1e-12, "Re={re:e}: {d:e}");
    }
}

/// Physics sanity: for any ship, a rougher hull or a faster ship never gives a *smaller*
/// relative friction increase (rough-wall effects grow with Reynolds number).
#[test]
fn granville_monotone_in_roughness_and_speed() {
    for l in [60.0, 120.0, 176.0, 300.0, 400.0] {
        for ks in [30.0, 60.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0] {
            let mut prev = f64::NEG_INFINITY;
            for v in [4.0, 8.0, 12.0, 16.0, 20.0, 25.0] {
                let d = delta_cf_ratio_direct(reynolds(Metres(l), Knots(v)), ks * 1e-6 / l);
                assert!(
                    d >= prev - 1e-12,
                    "speed: L={l} ks={ks} v={v}: {d} < {prev}"
                );
                prev = d;
            }
        }
        let re = reynolds(Metres(l), Knots(14.0));
        let mut prev = f64::NEG_INFINITY;
        for i in 0..60 {
            let ks = 10f64.powf(-6.0 + 4.0 * i as f64 / 59.0); // 1 µm .. 10 mm
            let d = delta_cf_ratio_direct(re, ks / l);
            assert!(
                d >= prev - 1e-12 && d >= -1e-12,
                "roughness: L={l} ks={ks}: {d} < {prev}"
            );
            prev = d;
        }
    }
}

/// The age-clock rewrite rests on Uzun's Gaussian curve solving the ODE
/// dFR/dA = (2/τ)·FR·√ln(max/FR) (see `uzun.rs`). Checked by numerical differentiation
/// at 200 ages on both curves.
#[test]
fn uzun_gaussian_solves_square_root_gompertz_ode() {
    let mut rng = Rng(103);
    for c in [EQ_SLIME, EQ_SOFT] {
        for _ in 0..200 {
            let a = rng.uniform(1.0, c.t0 - 1.0);
            let h = 1e-4;
            let numeric = (c.rating(a + h) - c.rating(a - h)) / (2.0 * h);
            let fr = c.rating(a);
            let ode = 2.0 / c.tau * fr * (c.max / fr).ln().sqrt();
            assert!(rel_close(numeric, ode, 1e-6), "A={a}: {numeric} vs {ode}");
        }
    }
}

/// rating → age → rating returns the same rating to 1e-12 at 1000 points spanning the
/// whole curve; reset maps depend on this inverse.
#[test]
fn rating_inverse_round_trip() {
    for c in [EQ_SLIME, EQ_SOFT] {
        for i in 1..1000 {
            let fr = c.seed() + (c.max - c.seed()) * i as f64 / 1000.0;
            assert!(
                rel_close(c.rating(c.age_for_rating(fr)), fr, 1e-12),
                "fr={fr}"
            );
        }
    }
}

/// coverage → age → coverage round trip with random (temperature-dependent) asymptotes.
#[test]
fn coverage_inverse_round_trip_random_asymptotes() {
    let mut rng = Rng(104);
    for _ in 0..500 {
        let (main, fast) = (rng.uniform(0.0, 100.0), rng.uniform(0.0, 10.0));
        let (lo, hi) = (
            EQ_CALCAREOUS.coverage(0.0, main, fast),
            EQ_CALCAREOUS.coverage(3000.0, main, fast),
        );
        let sc = rng.uniform(lo, hi);
        let age = EQ_CALCAREOUS.age_for_coverage(sc, main, fast, 3000.0);
        assert!(
            rel_close(EQ_CALCAREOUS.coverage(age, main, fast), sc, 1e-9),
            "sc={sc} main={main} fast={fast}"
        );
    }
}

/// The hull penalty scales exactly with the frictional share of resistance (a ship whose
/// resistance is 25% friction suffers a quarter of the full-friction penalty).
#[test]
fn hull_penalty_linear_in_frictional_share() {
    let table = GranvilleTable::new();
    let full = ShipParticulars::new(Metres(200.0), 1.0).unwrap();
    for share in [0.1, 0.25, 0.5, 0.8] {
        let ship = ShipParticulars::new(Metres(200.0), share).unwrap();
        for ks in [50.0, 500.0, 5000.0] {
            let f = |s: &ShipParticulars| {
                fouling_resistance_fraction(
                    &OperatingPoint::new(s, Knots(15.0), &table).unwrap(),
                    Micrometres(ks),
                )
            };
            assert!(rel_close(f(&ship), share * f(&full), 1e-14));
        }
    }
}

fn clean_days(policy: CleaningPolicy) -> Vec<f64> {
    let traj =
        with_tanker(|op| simulate(&mut exposure(), &glofouling_profile(5.0, 25.0), policy, op));
    traj.iter()
        .filter(|p| matches!(p.event, Some(MaintenanceEvent::InWaterHullClean { .. })))
        .map(|p| p.elapsed.value())
        .collect()
}

fn assert_days(got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len(), "{got:?}");
    for (g, w) in got.iter().zip(want) {
        assert!((g - w).abs() < 1e-6, "{got:?} vs {want:?}");
    }
}

/// Proactive policy on the 5-year equatorial profile, worked out by hand:
/// 20% needs 166.67 idle days ÷ 3.0411 per cycle → the 55th berth stay, ending day
/// 543.0411. Later cleans happen at the first segment end ≥ 182.5 days after the last:
/// 730, 913.0411, 1100, 1283.0411, 1470, 1653.0411 (the next, 1835.5, is past day 1830).
#[test]
fn proactive_schedule_matches_hand_calculation() {
    assert_days(
        &clean_days(CleaningPolicy::PROACTIVE),
        &[
            543.0411, 730.0, 913.0411, 1100.0, 1283.0411, 1470.0, 1653.0411,
        ],
    );
}

/// Reactive policy, by hand: 40% needs 333.33 idle days → the 110th berth stay, ending
/// day 1093.0411; then every ≥ 365 days at segment ends: 1460 and 1830.
#[test]
fn reactive_schedule_matches_hand_calculation() {
    assert_days(
        &clean_days(CleaningPolicy::REACTIVE),
        &[1093.0411, 1460.0, 1830.0],
    );
}

/// Harder splitting test than the equal-pieces one in `fouling_synthetic`: random,
/// *uneven* cut points, and the full salinity range including freshwater (which now
/// culls per stay). One segment and its pieces must give the same state in both stages.
#[test]
fn uneven_random_splits_give_identical_state() {
    let mut rng = Rng(105);
    for _ in 0..300 {
        let seg = random_segment(&mut rng);
        let total = seg.duration().value();
        let mut cuts: Vec<f64> = (0..rng.index(8) + 1)
            .map(|_| rng.uniform(0.0, total))
            .collect();
        cuts.push(0.0);
        cuts.push(total);
        cuts.sort_by(f64::total_cmp);
        let (mut w1, mut w2, mut p1, mut p2) = (exposure(), age_clock(), exposure(), age_clock());
        w1.advance(&seg);
        w2.advance(&seg);
        for w in cuts.windows(2) {
            let piece = with_duration(&seg, w[1] - w[0]);
            p1.advance(&piece);
            p2.advance(&piece);
        }
        assert!(
            rel_close(w1.hull_exposure().value(), p1.hull_exposure().value(), 1e-9),
            "{seg:?}"
        );
        for (a, b, what) in [
            (w2.slime_rating(), p2.slime_rating(), "slime"),
            (w2.soft_rating(), p2.soft_rating(), "soft"),
            (
                w2.calcareous_coverage(),
                p2.calcareous_coverage(),
                "coverage",
            ),
            (w2.prop_ks().value(), p2.prop_ks().value(), "prop_ks"),
        ] {
            assert!(
                rel_close(a, b, 1e-9),
                "{what}: {a} vs {b} for {seg:?} cuts {cuts:?}"
            );
        }
    }
}

/// 20,000 random segments with maintenance events mixed in (about 1 in 20), both stages.
/// After every step: every penalty finite and non-negative, power factor ≥ 1, exposure model
/// within its caps, age-clock model ratings within class maxima and coverage within [0, 100],
/// and no event ever makes the hull or propeller rougher.
#[test]
fn long_random_history_keeps_invariants() {
    let mut rng = Rng(106);
    with_tanker(|op| {
        let (mut s1, mut s2) = (exposure(), age_clock());
        for step in 0..20_000 {
            if rng.index(20) == 0 {
                let event = match rng.index(5) {
                    0 => MaintenanceEvent::DryDockFullBlast,
                    1 => MaintenanceEvent::DryDockSpotRepair,
                    2 => MaintenanceEvent::InWaterHullClean {
                        calcareous_remnant: rng.uniform(0.0, 1.0),
                    },
                    3 => MaintenanceEvent::HullGrooming,
                    _ => MaintenanceEvent::PropellerPolish,
                };
                let (h1, h2, p2) = (
                    s1.hull_fouling_fraction(),
                    s2.hull_ks().value(),
                    s2.prop_ks().value(),
                );
                s1.apply(event);
                s2.apply(event);
                assert!(s1.hull_fouling_fraction() <= h1, "step {step} {event:?}");
                assert!(
                    s2.hull_ks().value() <= h2 && s2.prop_ks().value() <= p2,
                    "step {step} {event:?}"
                );
            } else {
                let seg = random_segment(&mut rng);
                s1.advance(&seg);
                s2.advance(&seg);
            }
            assert!(
                (0.0..=0.86).contains(&s1.hull_fouling_fraction())
                    && (0.0..=0.32).contains(&s1.prop_eta_loss())
            );
            assert!(s2.slime_rating() <= 20.0 && s2.soft_rating() <= 50.0);
            assert!((0.0..=100.0).contains(&s2.calcareous_coverage()));
            for p in [s1.power_penalty(op), s2.power_penalty(op)] {
                for v in [
                    p.hull_fouling_fraction,
                    p.coating_fraction,
                    p.propeller_eta_loss,
                ] {
                    assert!(v.is_finite() && v >= 0.0, "step {step}: {p:?}");
                }
                assert!(
                    p.power_factor() >= 1.0 && p.power_factor().is_finite(),
                    "step {step}: {p:?}"
                );
            }
        }
    });
}

/// 20,000 random world segments (−2–34 °C, 0–38 PSU, 0–20 kn, any port, chlorophyll 0–20 mg/m³
/// on half of them) with maintenance events mixed in (about 1 in 20), for the guild model.
/// After every step: penalties finite and non-negative, power factor ≥ 1, slime ≤ 20,
/// soft ≤ 50, live + dead coverage ≤ 100 %, and no event makes the hull or propeller rougher.
#[test]
fn long_random_history_guild_model() {
    use nautilus_ocean::fouling::guild_model::GuildModel;
    use nautilus_ocean::fouling::params::GuildParams;
    use nautilus_units::MilligramsPerCubicMetre;
    let mut rng = Rng(107);
    with_tanker(|op| {
        for p in [
            GuildParams::default(),
            GuildParams {
                food_alpha: 0.5,
                ..GuildParams::default()
            },
        ] {
            let mut m = GuildModel::new(p).unwrap();
            for step in 0..20_000 {
                if rng.index(20) == 0 {
                    let event = match rng.index(5) {
                        0 => MaintenanceEvent::DryDockFullBlast,
                        1 => MaintenanceEvent::DryDockSpotRepair,
                        2 => MaintenanceEvent::InWaterHullClean {
                            calcareous_remnant: rng.uniform(0.0, 1.0),
                        },
                        3 => MaintenanceEvent::HullGrooming,
                        _ => MaintenanceEvent::PropellerPolish,
                    };
                    let (h, pk) = (m.hull_ks(), m.prop_ks());
                    m.apply(event);
                    assert!(
                        m.hull_ks() <= h && m.prop_ks() <= pk,
                        "step {step} {event:?}"
                    );
                } else {
                    let mut seg = random_segment(&mut rng);
                    if rng.index(2) == 0 {
                        seg = seg
                            .with_chlorophyll(MilligramsPerCubicMetre(rng.uniform(0.0, 20.0)))
                            .unwrap();
                    }
                    m.advance(&seg);
                }
                assert!(
                    m.slime_rating() <= 20.0 && m.soft_rating() <= 50.0,
                    "step {step}"
                );
                assert!(
                    m.calcareous_coverage() <= 100.0 && m.live_coverage() >= 0.0,
                    "step {step}"
                );
                let pen = m.power_penalty(op);
                for v in [
                    pen.hull_fouling_fraction,
                    pen.coating_fraction,
                    pen.propeller_eta_loss,
                ] {
                    assert!(v.is_finite() && v >= 0.0, "step {step}: {pen:?}");
                }
                assert!(pen.power_factor() >= 1.0, "step {step}: {pen:?}");
            }
        }
    });
}
