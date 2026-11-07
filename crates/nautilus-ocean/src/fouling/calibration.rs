//! Layered calibration of the guild model (spec Section 6): reference data, the scoring
//! function of each layer, and a deterministic Nelder–Mead optimizer.
//!
//! The fitting itself runs in `examples/calibrate_fouling.rs`; its results are hard-coded
//! in [`GuildParams::default`]. The validation tests call the same scoring functions to check
//! that the defaults meet every pass criterion.

use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};

use super::environment::{PortStatus, Segment};
use super::events::MaintenanceEvent;
use super::granville::GranvilleTable;
use super::guild_model::GuildModel;
use super::params::GuildParams;
use super::power::ShipParticulars;
use super::simulate::{CleaningPolicy, TimelineItem, simulate};
use super::uzun::{EQ_CALCAREOUS, EQ_MAIN_ASYM, EQ_SLIME, EQ_SOFT, GaussianRating};
use super::{FoulingModel, OperatingPoint};

/// Uzun 2019 Table 2, Mediterranean (17 °C) slime curve.
pub const MED_SLIME: GaussianRating = GaussianRating {
    max: 20.0,
    t0: 271.9,
    tau: 99.31,
};
/// Uzun 2019 Table 2, Mediterranean (17 °C) soft-fouling curve.
pub const MED_SOFT: GaussianRating = GaussianRating {
    max: 50.0,
    t0: 383.5,
    tau: 124.4,
};

/// Uzun 2019 Table 3, Mediterranean calcareous coverage (%) at `age_days` (two logistics,
/// 50 % plateau).
pub fn med_coverage(age_days: f64) -> f64 {
    0.00517 / (1.0 + (10.0 - 40.0 * age_days).exp())
        + 50.0 / (1.0 + (32.81 - 0.04715 * age_days).exp())
}

/// Uzun 2019 Table 3, equatorial calcareous coverage (%) at `age_days`: the **main**
/// logistic only. The small fast term (3.5 % within about two weeks) is excluded: the research
/// spec lists it as unconfirmed (open item 1), and it stays below the 5 % coverage at which
/// Uzun's roughness mapping starts, so it never affects roughness or power. Decided with the
/// user on 2026-11-07.
pub fn eq_coverage(age_days: f64) -> f64 {
    EQ_CALCAREOUS.coverage(age_days, EQ_MAIN_ASYM, 0.0)
}

/// Worst absolute gap between the model and Uzun's curve at one site.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiteGaps {
    /// Slime rating points.
    pub slime: f64,
    /// Soft-fouling rating points.
    pub soft: f64,
    /// Calcareous coverage percentage points.
    pub coverage: f64,
}

/// Layer 1 score: worst gaps at the equatorial (25 °C) and Mediterranean (17 °C) sites, over
/// 1,500 daily idle segments at a berth in 35 PSU water. Invalid parameters score infinity.
pub fn layer1_gaps(p: &GuildParams) -> [SiteGaps; 2] {
    let site = |t: f64, slime: GaussianRating, soft: GaussianRating, cov: fn(f64) -> f64| {
        let inf = SiteGaps {
            slime: f64::INFINITY,
            soft: f64::INFINITY,
            coverage: f64::INFINITY,
        };
        let Ok(mut m) = GuildModel::new(*p) else {
            return inf;
        };
        let Ok(day) = Segment::new(
            Days(1.0),
            Knots(0.0),
            Celsius(t),
            Psu(35.0),
            PortStatus::Berth,
            Kilometres(0.0),
        ) else {
            return inf;
        };
        let mut g = SiteGaps {
            slime: 0.0,
            soft: 0.0,
            coverage: 0.0,
        };
        for d in 1..=1500 {
            m.advance(&day);
            let a = d as f64;
            // NaN-safe running maxima: a NaN gap becomes infinity.
            let worse = |cur: f64, gap: f64| {
                if gap.is_nan() {
                    f64::INFINITY
                } else {
                    cur.max(gap)
                }
            };
            g.slime = worse(g.slime, (m.slime_rating() - slime.rating(a)).abs());
            g.soft = worse(g.soft, (m.soft_rating() - soft.rating(a)).abs());
            g.coverage = worse(g.coverage, (m.calcareous_coverage() - cov(a)).abs());
        }
        g
    };
    [
        site(25.0, EQ_SLIME, EQ_SOFT, eq_coverage),
        site(17.0, MED_SLIME, MED_SOFT, med_coverage),
    ]
}

/// GloFouling virtual profile: 10-day cycles of 3.0411 d idle at a berth and 6.9589 d sailing
/// at 14 kn 300 km offshore (556 idle days in 5 years), at a constant SST, 35 PSU.
pub fn glofouling_profile(years: f64, sst: f64) -> Vec<TimelineItem> {
    let cycles = (years * 36.525).ceil() as usize;
    let seg = |d, v, port, km| {
        Segment::new(
            Days(d),
            Knots(v),
            Celsius(sst),
            Psu(35.0),
            port,
            Kilometres(km),
        )
    };
    match (
        seg(3.0411, 0.0, PortStatus::Berth, 0.0),
        seg(6.9589, 14.0, PortStatus::AtSea, 300.0),
    ) {
        (Ok(idle), Ok(sail)) => (0..cycles)
            .flat_map(|_| [TimelineItem::Segment(idle), TimelineItem::Segment(sail)])
            .collect(),
        _ => Vec::new(),
    }
}

/// Runs `f` at the calibration operating point: 176 m tanker, frictional share 0.5, 14 kn.
fn with_tanker<R>(f: impl FnOnce(&OperatingPoint) -> R) -> Option<R> {
    let ship = ShipParticulars::new(Metres(176.0), 0.5).ok()?;
    let table = GranvilleTable::new();
    let op = OperatingPoint::new(&ship, Knots(14.0), &table).ok()?;
    Some(f(&op))
}

/// Layer 2 score: 5-year mean hull fouling fraction without cleaning on the GloFouling
/// profile, `(equatorial 25 °C, Mediterranean 17 °C)`. Invalid parameters give infinity.
pub fn layer2_penalties(p: &GuildParams) -> (f64, f64) {
    let run = |sst: f64| {
        let Ok(mut m) = GuildModel::new(*p) else {
            return f64::INFINITY;
        };
        with_tanker(|op| {
            let traj = simulate(
                &mut m,
                &glofouling_profile(5.0, sst),
                CleaningPolicy::NONE,
                op,
            );
            let total: f64 = traj.iter().map(|t| t.segment_duration.value()).sum();
            traj.iter()
                .map(|t| t.segment_duration.value() * t.penalty.hull_fouling_fraction)
                .sum::<f64>()
                / total
        })
        .unwrap_or(f64::INFINITY)
    };
    (run(25.0), run(17.0))
}

/// Layer 3 score: after 1,100 days (220 profile items) on the equatorial profile, the power
/// recovered by a full-blast dry-dock divided by that recovered by an in-water clean leaving
/// 20 % base plates (each followed by a propeller polish). Adland et al. find about 2.
pub fn layer3_ratio(p: &GuildParams) -> f64 {
    let Ok(mut m) = GuildModel::new(*p) else {
        return f64::INFINITY;
    };
    for item in glofouling_profile(3.0, 25.0).iter().take(220) {
        if let TimelineItem::Segment(s) = item {
            m.advance(s);
        }
    }
    with_tanker(|op| {
        let before = m.power_penalty(op).power_factor();
        let after = |e: MaintenanceEvent| {
            let mut c = m.clone();
            c.apply(e);
            c.apply(MaintenanceEvent::PropellerPolish);
            c.power_penalty(op).power_factor()
        };
        (before - after(MaintenanceEvent::DryDockFullBlast))
            / (before
                - after(MaintenanceEvent::InWaterHullClean {
                    calcareous_remnant: 0.2,
                }))
    })
    .unwrap_or(f64::INFINITY)
}

/// Minimises `f` by Nelder–Mead (reflection 1, expansion 2, contraction 0.5, shrink 0.5)
/// from `x0` with initial simplex steps `step`, for exactly `iterations` iterations.
/// Non-finite objective values count as +∞. Deterministic.
pub fn nelder_mead(
    f: impl Fn(&[f64]) -> f64,
    x0: &[f64],
    step: &[f64],
    iterations: usize,
) -> Vec<f64> {
    let n = x0.len();
    let eval = |x: &[f64]| {
        let v = f(x);
        if v.is_finite() { v } else { f64::INFINITY }
    };
    let mut simplex: Vec<(Vec<f64>, f64)> = (0..=n)
        .map(|i| {
            let mut x = x0.to_vec();
            if i > 0 {
                x[i - 1] += step[i - 1];
            }
            let v = eval(&x);
            (x, v)
        })
        .collect();
    let lerp = |a: &[f64], b: &[f64], t: f64| {
        a.iter()
            .zip(b)
            .map(|(a, b)| a + t * (b - a))
            .collect::<Vec<_>>()
    };
    for _ in 0..iterations {
        simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
        let centroid: Vec<f64> = (0..n)
            .map(|j| simplex[..n].iter().map(|p| p.0[j]).sum::<f64>() / n as f64)
            .collect();
        let worst = simplex[n].clone();
        let reflected = lerp(&centroid, &worst.0, -1.0);
        let fr = eval(&reflected);
        if fr < simplex[0].1 {
            let expanded = lerp(&centroid, &worst.0, -2.0);
            let fe = eval(&expanded);
            simplex[n] = if fe < fr {
                (expanded, fe)
            } else {
                (reflected, fr)
            };
        } else if fr < simplex[n - 1].1 {
            simplex[n] = (reflected, fr);
        } else {
            let contracted = if fr < worst.1 {
                lerp(&centroid, &reflected, 0.5)
            } else {
                lerp(&centroid, &worst.0, 0.5)
            };
            let fc = eval(&contracted);
            if fc < worst.1.min(fr) {
                simplex[n] = (contracted, fc);
            } else {
                let best = simplex[0].0.clone();
                for p in simplex.iter_mut().skip(1) {
                    p.0 = lerp(&best, &p.0, 0.5);
                    p.1 = eval(&p.0);
                }
            }
        }
    }
    simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
    simplex.swap_remove(0).0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The optimizer finds the Rosenbrock minimum (1, 1) from the classic start (−1.2, 1).
    #[test]
    fn nelder_mead_finds_rosenbrock_minimum() {
        let f = |x: &[f64]| (1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0] * x[0]).powi(2);
        let x = nelder_mead(f, &[-1.2, 1.0], &[0.5, 0.5], 2000);
        assert!(
            (x[0] - 1.0).abs() < 1e-3 && (x[1] - 1.0).abs() < 1e-3,
            "{x:?}"
        );
    }

    /// The optimizer is deterministic: the same call gives bit-identical results.
    #[test]
    fn nelder_mead_is_deterministic() {
        let f = |x: &[f64]| (x[0] - 3.0).powi(2) + (x[1] + 1.0).abs();
        assert_eq!(
            nelder_mead(f, &[0.0, 0.0], &[1.0, 1.0], 300),
            nelder_mead(f, &[0.0, 0.0], &[1.0, 1.0], 300)
        );
    }

    /// Reference curves: Uzun's Mediterranean coverage is half its 50 % plateau at its
    /// 695.86-day midpoint, and the equatorial main term is half its 96.899 % plateau at its
    /// 393.12-day midpoint (Uzun 2019 Table 3; fast term excluded, see `eq_coverage`).
    #[test]
    fn reference_curves() {
        assert!((med_coverage(695.86) - 25.0).abs() < 0.01);
        assert!((eq_coverage(393.1204) - 96.899 / 2.0).abs() < 1e-3);
        assert_eq!((MED_SLIME.t0, MED_SOFT.t0), (271.9, 383.5));
    }
}
