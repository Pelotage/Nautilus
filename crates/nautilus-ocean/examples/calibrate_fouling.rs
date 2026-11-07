//! Layered calibration of the guild model (spec Section 6).
//!
//! ```text
//! cargo run --release -p nautilus-ocean --example calibrate_fouling -- layer1
//! cargo run --release -p nautilus-ocean --example calibrate_fouling -- layer2
//! cargo run --release -p nautilus-ocean --example calibrate_fouling -- layer3
//! ```
//!
//! Each layer starts from the current `GuildParams::default()` (so earlier layers' fitted
//! values are fixed), fits only its own parameters with a deterministic Nelder–Mead (fixed
//! iterations, fixed starts, log-space where a parameter must stay positive), and prints the
//! fitted constants and the score achieved. Copy the printed values into
//! `GuildParams::default()`.

use nautilus_ocean::fouling::calibration::{
    layer1_gaps, layer2_penalties, layer3_ratio, nelder_mead,
};
use nautilus_ocean::fouling::params::{CALCAREOUS, GuildParams, SLIME, SOFT};

/// Nelder–Mead with `rounds` restarts from the best point so far (fresh simplex each time).
fn fit(
    f: impl Fn(&[f64]) -> f64,
    x0: &[f64],
    step: &[f64],
    rounds: usize,
    iterations: usize,
) -> Vec<f64> {
    let mut x = x0.to_vec();
    for _ in 0..rounds {
        x = nelder_mead(&f, &x, step, iterations);
    }
    x
}

fn with_class(base: &GuildParams, c: usize, v: &[f64]) -> GuildParams {
    let mut p = *base;
    let r = &mut p.classes[c];
    (r.settle_ref, r.growth_ref, r.settle_q10, r.growth_q10) =
        (v[0].exp(), v[1].exp(), v[2].exp(), v[3].exp());
    if c == CALCAREOUS {
        (p.capacity_17, p.capacity_25) = (v[4], v[5]);
    }
    p
}

fn layer1(base: GuildParams) {
    // Analytic starts: growth 2/τ (Uzun Gaussian widths) and logistic rates; Q10s from the
    // ratio of the two sites' rates.
    // Calcareous settlement starts are the analytic logistic-with-settlement values that put
    // the coverage midpoints at Uzun's 393 d (25 °C) and 696 d (17 °C): s ≈ r·K·e^(−r·t_mid).
    let starts: [(usize, &str, f64, Vec<f64>); 3] = [
        (
            SLIME,
            "slime",
            1.0,
            vec![
                0.05f64.ln(),
                (2.0f64 / 37.08).ln(),
                3.0f64.ln(),
                3.44f64.ln(),
            ],
        ),
        (
            SOFT,
            "soft",
            2.0,
            vec![
                0.002f64.ln(),
                (2.0f64 / 73.11).ln(),
                3.0f64.ln(),
                1.94f64.ln(),
            ],
        ),
        (
            CALCAREOUS,
            "calcareous",
            3.0,
            vec![
                4.6e-7f64.ln(),
                0.0407f64.ln(),
                2.8e9f64.ln(),
                0.83f64.ln(),
                50.0,
                98.6,
            ],
        ),
    ];
    let mut p = base;
    for (c, name, threshold, x0) in starts {
        let current = p;
        let score = |v: &[f64]| {
            let g = layer1_gaps(&with_class(&current, c, v));
            let pick = |s: &nautilus_ocean::fouling::calibration::SiteGaps| match c {
                SLIME => s.slime,
                SOFT => s.soft,
                _ => s.coverage,
            };
            pick(&g[0]).max(pick(&g[1])) / threshold
        };
        let step: Vec<f64> = x0
            .iter()
            .enumerate()
            .map(|(i, _)| {
                if c == CALCAREOUS && (i == 4 || i == 5) {
                    10.0
                } else {
                    0.5
                }
            })
            .collect();
        let v = fit(score, &x0, &step, 4, 600);
        let achieved = score(&v);
        p = with_class(&p, c, &v);
        let r = p.classes[c];
        println!(
            "{name}: settle_ref {:e}, growth_ref {:e}, settle_q10 {:e}, growth_q10 {:e}; score {:.4} (≤ 1 passes)",
            r.settle_ref, r.growth_ref, r.settle_q10, r.growth_q10, achieved
        );
    }
    println!(
        "capacity_17 {:e}, capacity_25 {:e} (weathering fixed at {:e}/day)",
        p.capacity_17, p.capacity_25, p.weathering_per_day
    );
    let g = layer1_gaps(&p);
    for (site, s) in ["equatorial", "mediterranean"].iter().zip(g) {
        println!(
            "{site}: slime {:.3} (≤ 1), soft {:.3} (≤ 2), coverage {:.3} (≤ 3)",
            s.slime, s.soft, s.coverage
        );
    }
}

fn layer2(base: GuildParams) {
    // One roughness scale λ fitted to both GloFouling sites at once.
    let set = |v: &[f64]| GuildParams {
        roughness_scale: v[0],
        ..base
    };
    let score = |v: &[f64]| {
        let (eq, med) = layer2_penalties(&set(v));
        (eq - 0.30).powi(2) + (med - 0.22).powi(2)
    };
    let v = fit(score, &[base.roughness_scale], &[0.2], 2, 200);
    let (eq, med) = layer2_penalties(&set(&v));
    println!(
        "roughness_scale {:e}; equatorial {eq:.4} (0.30 ± 0.05), mediterranean {med:.4} (0.22 ± 0.05)",
        v[0]
    );
}

fn layer3(base: GuildParams) {
    // Cleaning efficiency η fitted through a logistic transform (keeps 0 < η < 1); paint
    // damage stays at its 10 µm prior.
    let eta = |u: f64| 1.0 / (1.0 + (-u).exp());
    let set = |v: &[f64]| GuildParams {
        clean_efficiency: eta(v[0]),
        ..base
    };
    let v = fit(
        |v| (layer3_ratio(&set(v)) - 2.0).powi(2),
        &[0.0],
        &[1.0],
        2,
        200,
    );
    println!(
        "clean_efficiency {:e}; dry-dock / in-water ratio {:.4} (1.6–2.4)",
        eta(v[0]),
        layer3_ratio(&set(&v))
    );
}

fn main() {
    let base = GuildParams::default();
    match std::env::args().nth(1).as_deref() {
        Some("layer1") => layer1(base),
        Some("layer2") => layer2(base),
        Some("layer3") => layer3(base),
        _ => eprintln!("usage: calibrate_fouling layer1|layer2|layer3"),
    }
}
