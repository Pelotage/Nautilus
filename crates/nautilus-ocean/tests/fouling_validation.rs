//! Validation against **published data that the model was not fitted to**, or was
//! fitted to only in part.
//!
//! Passing tests here show where the model agrees with the literature. Tests marked
//! `#[ignore = "calibration target: …"]` are known gaps: they fail today and document
//! how far the model is from a published result. Run them with
//! `cargo test -p nautilus-ocean --test fouling_validation -- --ignored`.
//!
//! Sources (all quoted in `research/reports/Biofouling growth dynamical model.md`):
//! Uzun et al. 2019 Tables 2–3 (Mediterranean panel site); GloFouling 2022 (5-year
//! cumulative penalties); Adland et al. (dry-dock vs in-water recovery); Hundley via
//! GloFouling (destroyer cleaning sequence).

mod common;

use common::*;
use nautilus_ocean::fouling::FoulingModel;
use nautilus_ocean::fouling::calibration::{layer1_gaps, layer2_penalties, layer3_ratio};
use nautilus_ocean::fouling::events::MaintenanceEvent;
use nautilus_ocean::fouling::guild_model::GuildModel;
use nautilus_ocean::fouling::params::GuildParams;
use nautilus_ocean::fouling::simulate::{CleaningPolicy, TimelineItem, simulate};
use nautilus_ocean::fouling::uzun::{GaussianRating, hull_ks};

/// Uzun 2019 Table 2, Mediterranean site (17 °C): slime and soft-fouling Gaussians.
const MED_SLIME: GaussianRating = GaussianRating {
    max: 20.0,
    t0: 271.9,
    tau: 99.31,
};
const MED_SOFT: GaussianRating = GaussianRating {
    max: 50.0,
    t0: 383.5,
    tau: 124.4,
};
/// Uzun 2019 Table 3, Mediterranean calcareous coverage (two logistics; 50% plateau).
fn med_coverage(age: f64) -> f64 {
    0.00517 / (1.0 + (10.0 - 40.0 * age).exp()) + 50.0 / (1.0 + (32.81 - 0.04715 * age).exp())
}

/// Runs age-clock model day by day at the Mediterranean site (17 °C, berth, 35 PSU) and returns,
/// for each day, (slime, soft, coverage, hull k_s µm) from the model and from Uzun's curves.
fn med_comparison(days: usize) -> Vec<([f64; 4], [f64; 4])> {
    let mut m = age_clock();
    (1..=days)
        .map(|d| {
            m.advance(&idle(1.0, 17.0, 35.0));
            let a = d as f64;
            let (sl, so, sc) = (MED_SLIME.rating(a), MED_SOFT.rating(a), med_coverage(a));
            (
                [
                    m.slime_rating(),
                    m.soft_rating(),
                    m.calcareous_coverage(),
                    m.hull_ks().value(),
                ],
                [sl, so, sc, hull_ks(sl + so, sc).value()],
            )
        })
        .collect()
}

fn first_day(series: &[([f64; 4], [f64; 4])], k: usize, level: f64) -> usize {
    series
        .iter()
        .position(|(m, _)| m[k] >= level)
        .expect("never reached")
        + 1
}

/// Timing at the Mediterranean site: the model's slime plateaus on day 272 and soft
/// fouling on day 384, against Uzun's published 271.9 and 383.5 days (±1 day), and
/// coverage reaches half its 50% plateau on day 696 (Uzun midpoint 695.9).
/// These are the quantities the Q10 values were fitted to, so this checks the fit and
/// the age-clock machinery, not an independent prediction.
#[test]
fn med_site_timing_matches_uzun() {
    let s = med_comparison(1500);
    assert!(
        (first_day(&s, 0, 20.0) as f64 - 271.9).abs() <= 1.0,
        "slime plateau day {}",
        first_day(&s, 0, 20.0)
    );
    assert!(
        (first_day(&s, 1, 50.0) as f64 - 383.5).abs() <= 1.0,
        "soft plateau day {}",
        first_day(&s, 1, 50.0)
    );
    assert!(
        (first_day(&s, 2, 25.0) as f64 - 695.9).abs() <= 1.0,
        "coverage half day {}",
        first_day(&s, 2, 25.0)
    );
    let plateau = s.last().unwrap();
    assert!(
        (plateau.0[2] - plateau.1[2]).abs() < 0.05,
        "plateau {} vs {}",
        plateau.0[2],
        plateau.1[2]
    );
}

/// Shape at the Mediterranean site, as a regression guard. Scaling time by a Q10 cannot
/// change the *shape* of a curve, and Uzun's Mediterranean curves have different widths
/// (slime τ 99 d vs the model's 116 d; soft 124 vs 103; coverage rate 0.047/d vs 0.023/d).
/// Measured worst-case gaps over 1500 days: slime 2.3, soft 6.8 rating points, coverage
/// 7.8 percentage points, hull roughness 1.4 mm (around day 750). This test fails if a
/// change makes any of them worse; `med_site_shape_calibration_target` is the goal.
#[test]
fn med_site_shape_gap_does_not_grow() {
    let s = med_comparison(1500);
    let worst = |k: usize| {
        s.iter()
            .map(|(m, u)| (m[k] - u[k]).abs())
            .fold(0.0, f64::max)
    };
    for (k, bound, what) in [
        (0, 2.4, "slime"),
        (1, 6.9, "soft"),
        (2, 7.9, "coverage %"),
        (3, 1450.0, "hull k_s µm"),
    ] {
        assert!(
            worst(k) <= bound,
            "{what}: worst gap {} > {bound}",
            worst(k)
        );
    }
}

/// Layer 1 (growth) of the guild model: at both Uzun panel sites (25 °C equatorial and
/// 17 °C Mediterranean, idle at berth in sea water, 1,500 days) the worst gap to Uzun's
/// published curves is at most 1 rating point for slime, 2 for soft fouling and 3 percentage
/// points for calcareous coverage. Oracle: Uzun et al. 2019 Tables 2–3. This replaces the
/// age-clock model's calibration target, which could match timing but not shape.
#[test]
fn guild_layer1_growth_matches_uzun_both_sites() {
    for (site, g) in ["equatorial", "mediterranean"]
        .iter()
        .zip(layer1_gaps(&GuildParams::default()))
    {
        assert!(
            g.slime <= 1.0 && g.soft <= 2.0 && g.coverage <= 3.0,
            "{site}: slime {:.3} soft {:.3} coverage {:.3}",
            g.slime,
            g.soft,
            g.coverage
        );
    }
}

/// Mean hull fouling penalty over 5 years without cleaning on the GloFouling profile.
fn five_year_mean_hull<M: FoulingModel>(mut model: M, sst: f64) -> f64 {
    let traj = with_tanker(|op| {
        simulate(
            &mut model,
            &glofouling_profile(5.0, sst),
            CleaningPolicy::NONE,
            op,
        )
    });
    let total: f64 = traj.iter().map(|p| p.segment_duration.value()).sum();
    traj.iter()
        .map(|p| p.segment_duration.value() * p.penalty.hull_fouling_fraction)
        .sum::<f64>()
        / total
}

/// GloFouling 2022 reports cumulative 5-year penalties without cleaning of +30%
/// (equatorial) and +22% (Mediterranean). Exposure model gives 0.336 and 0.193 — within 5
/// percentage points. Assumes GloFouling's Mediterranean ship has the same idle pattern
/// as the equatorial one, at 17 °C.
#[test]
fn glofouling_five_year_penalty_exposure() {
    let (eq, med) = (
        five_year_mean_hull(exposure(), 25.0),
        five_year_mean_hull(exposure(), 17.0),
    );
    assert!((eq - 0.30).abs() <= 0.05, "equatorial {eq}");
    assert!((med - 0.22).abs() <= 0.05, "mediterranean {med}");
    assert!(eq > med, "warm water must foul more");
}

/// Layer 2 (roughness → power) of the guild model: the 5-year mean hull fouling penalty
/// without cleaning on the GloFouling profile is within 5 percentage points of GloFouling's
/// +30 % (equatorial, 25 °C) and +22 % (Mediterranean, 17 °C), and warm water fouls more.
/// Oracle: GloFouling 2022. Replaces the age-clock model's calibration target (0.43 vs 0.30).
#[test]
fn guild_layer2_glofouling_five_year_penalties() {
    let (eq, med) = layer2_penalties(&GuildParams::default());
    assert!((eq - 0.30).abs() <= 0.05, "equatorial {eq}");
    assert!((med - 0.22).abs() <= 0.05, "mediterranean {med}");
    assert!(eq > med);
}

/// Layer 3 (cleaning and docking) of the guild model: after 1,100 days on the equatorial
/// profile, a full-blast dry-dock recovers 1.6–2.4 times as much power as an in-water clean
/// leaving 20 % base plates (each with a propeller polish). Oracle: Adland et al. (about 2:1
/// on tanker-like profiles). Replaces the exposure / age-clock calibration target (1.1–1.3).
#[test]
fn guild_layer3_dry_dock_recovers_twice_in_water_clean() {
    let ratio = layer3_ratio(&GuildParams::default());
    assert!((1.6..=2.4).contains(&ratio), "ratio {ratio}");
}

/// Hundley's destroyer (via GloFouling): +115% → +70% after propeller cleaning → +9%
/// after hull and propeller cleaning. The model is not calibrated to that ship, so this
/// checks only the pattern: polishing alone helps, cleaning the hull as well helps much
/// more, and a small non-zero penalty remains (paint roughness and the cleaning floor).
#[test]
fn hundley_cleaning_sequence_pattern() {
    hundley_pattern(age_clock());
}

/// The Hundley pattern for the guild model: propeller polish helps, the hull clean helps more,
/// and a small non-zero penalty remains. Oracle: Hundley via GloFouling (+70 % → +9 %).
///
/// Settled decision (user, 2026-11-07): Hundley's destroyer had a thorough naval clean, so this
/// runs with the whole hull cleanable (η = 1). The default η = 0.39 represents routine
/// commercial cleans (Adland, Layer 3). The two cannot share one default; do not re-litigate.
#[test]
fn guild_hundley_cleaning_sequence_pattern() {
    let thorough = GuildParams {
        clean_efficiency: 1.0,
        ..GuildParams::default()
    };
    hundley_pattern(GuildModel::new(thorough).unwrap());
}

fn hundley_pattern<M: FoulingModel>(mut m: M) {
    with_tanker(|op| {
        for item in &glofouling_profile(5.0, 25.0) {
            if let TimelineItem::Segment(s) = item {
                m.advance(s);
            }
        }
        let before = m.power_penalty(op).power_factor();
        m.apply(MaintenanceEvent::PropellerPolish);
        let after_prop = m.power_penalty(op).power_factor();
        m.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.0,
        });
        let after_both = m.power_penalty(op).power_factor();
        assert!(after_prop < before, "{after_prop} !< {before}");
        assert!(
            before - after_prop < after_prop - after_both,
            "hull clean should matter most"
        );
        assert!(
            after_both > 1.0 && after_both < 1.10,
            "residual {after_both}"
        );
    });
}
