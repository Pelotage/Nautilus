//! Granville similarity-law scaling from sand-grain roughness `k_s` to the
//! frictional resistance coefficient of a full-scale flat plate.
//!
//! The roughness function is the Nikuradse-type (sand-grain) function. The
//! Schultz / Uzun `k_s` scale used elsewhere in this module tree is tied to this
//! function: pairing that `k_s` with a Colebrook-type function instead gives a
//! clean-hull penalty about +27% at 30 µm (spec, §observation chain), so do not
//! swap the function without re-deriving the `k_s` mappings.
//!
//! All solvers are fixed-iteration bisections (64 steps), so results are
//! deterministic. Ported from `granville_prototype.py` in the research notes.

use std::sync::OnceLock;

/// Von Kármán constant used by the roughness function and Granville law.
pub const KAPPA: f64 = 0.42;

const BISECTION_STEPS: usize = 64;

/// Schoenherr smooth-plate friction line: solves `0.242/√C_F = log10(Re·C_F)`
/// by 64-step bisection on `C_F ∈ [1e-4, 0.02]`.
pub fn schoenherr_cf(re: f64) -> f64 {
    let (mut lo, mut hi) = (1e-4_f64, 0.02_f64);
    let mut mid = 0.5 * (lo + hi);
    for _ in 0..BISECTION_STEPS {
        mid = 0.5 * (lo + hi);
        let f = 0.242 / mid.sqrt() - (re * mid).log10();
        if f > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    mid
}

/// Nikuradse-type roughness function ΔU⁺(k⁺): 0 for k⁺ < 3 (hydraulically
/// smooth), `(1/κ)·ln(0.26·k⁺)` for k⁺ ≥ 15 (fully rough), and in the transitional
/// range the same multiplied by `sin((π/2)·log10(k⁺/3)/log10(5))`. Floored at 0.
pub fn nikuradse_roughness_function(k_plus: f64) -> f64 {
    if k_plus < 3.0 {
        return 0.0;
    }
    let mut v = (1.0 / KAPPA) * (0.26 * k_plus).ln();
    if k_plus < 15.0 {
        v *= (std::f64::consts::FRAC_PI_2 * (k_plus / 3.0).log10() / 5f64.log10()).sin();
    }
    v.max(0.0)
}

/// For a given plate length in wall units `L⁺`, solve for the rough-plate `C_F`
/// and the Reynolds number it corresponds to.
fn granville_inner(l_plus: f64, ks_over_l: f64) -> (f64, f64) {
    let d = nikuradse_roughness_function(ks_over_l * l_plus);
    let re_of = |cf: f64| {
        let u = (cf / 2.0).sqrt();
        l_plus / (u * (1.0 - u / KAPPA))
    };
    let (mut lo, mut hi) = (5e-4_f64, 0.03_f64);
    let mut mid = 0.5 * (lo + hi);
    for _ in 0..BISECTION_STEPS {
        mid = 0.5 * (lo + hi);
        let g = mid - schoenherr_cf(re_of(mid) * (-KAPPA * d).exp());
        if g > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (mid, re_of(mid))
}

/// Rough-plate frictional coefficient at Reynolds number `re` (based on plate
/// length) and relative roughness `ks_over_l`, by the Granville similarity law
/// with the Nikuradse-type roughness function.
pub fn rough_cf(re: f64, ks_over_l: f64) -> f64 {
    let (mut lo, mut hi) = (1e3_f64, 1e9_f64);
    let mut cf = 0.0;
    for _ in 0..BISECTION_STEPS {
        let l_plus = (lo * hi).sqrt();
        let (c, r) = granville_inner(l_plus, ks_over_l);
        cf = c;
        if r < re {
            lo = l_plus;
        } else {
            hi = l_plus;
        }
    }
    cf
}

/// Relative frictional increase `rough_cf/schoenherr_cf − 1`, solved directly.
pub fn delta_cf_ratio_direct(re: f64, ks_over_l: f64) -> f64 {
    rough_cf(re, ks_over_l) / schoenherr_cf(re) - 1.0
}

/// ITTC-78 (Townsin) roughness allowance
/// `0.044·((k_s/L)^(1/3) − 10·Re^(−1/3)) + 0.000125`. Used here for the
/// coating-roughness (AHR) term only.
pub fn ittc_roughness_allowance(ks_over_l: f64, re: f64) -> f64 {
    0.044 * (ks_over_l.cbrt() - 10.0 * re.powf(-1.0 / 3.0)) + 0.000125
}

const LOG_RE_MIN: f64 = 7.0;
const LOG_RE_MAX: f64 = 10.0;
const N_RE: usize = 32;
const LOG_KS_MIN: f64 = -8.0;
const LOG_KS_MAX: f64 = -3.0;
const N_KS: usize = 48;

/// Precomputed grid of [`delta_cf_ratio_direct`] over log10(Re) ∈ [7, 10]
/// (32 points) and log10(k_s/L) ∈ [−8, −3] (48 points), bilinearly
/// interpolated. Inputs outside the grid are clamped to its edges, never
/// extrapolated.
#[derive(Debug, Clone, PartialEq)]
pub struct GranvilleTable {
    values: Vec<f64>,
}

impl Default for GranvilleTable {
    fn default() -> Self {
        Self::new()
    }
}

impl GranvilleTable {
    /// Returns the table. The grid is solved once per process (the Granville
    /// law at every grid point) and shared; later calls copy the cached grid.
    pub fn new() -> GranvilleTable {
        static GRID: OnceLock<Vec<f64>> = OnceLock::new();
        let values = GRID.get_or_init(|| {
            let mut values = Vec::with_capacity(N_RE * N_KS);
            for i in 0..N_RE {
                let re = 10f64.powf(axis(LOG_RE_MIN, LOG_RE_MAX, N_RE, i));
                for j in 0..N_KS {
                    let ks = 10f64.powf(axis(LOG_KS_MIN, LOG_KS_MAX, N_KS, j));
                    values.push(delta_cf_ratio_direct(re, ks));
                }
            }
            values
        });
        GranvilleTable {
            values: values.clone(),
        }
    }

    /// Interpolated relative frictional increase `ΔC_F/C_F` at (`re`, `ks_over_l`).
    pub fn delta_cf_ratio(&self, re: f64, ks_over_l: f64) -> f64 {
        let (i, fi) = locate(re.log10(), LOG_RE_MIN, LOG_RE_MAX, N_RE);
        let (j, fj) = locate(ks_over_l.log10(), LOG_KS_MIN, LOG_KS_MAX, N_KS);
        let v = |i: usize, j: usize| self.values[i * N_KS + j];
        let a = v(i, j) * (1.0 - fj) + v(i, j + 1) * fj;
        let b = v(i + 1, j) * (1.0 - fj) + v(i + 1, j + 1) * fj;
        a * (1.0 - fi) + b * fi
    }
}

fn axis(min: f64, max: f64, n: usize, i: usize) -> f64 {
    min + (max - min) * i as f64 / (n - 1) as f64
}

/// Cell index and fractional position for `x` on an `n`-point axis, clamped.
/// Non-finite `x` (e.g. log of 0) clamps to the lower edge.
fn locate(x: f64, min: f64, max: f64, n: usize) -> (usize, f64) {
    let x = if x.is_nan() { min } else { x.clamp(min, max) };
    let pos = (x - min) / (max - min) * (n - 1) as f64;
    let i = (pos.floor() as usize).min(n - 2);
    (i, pos - i as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }
    fn re(l: f64, kn: f64) -> f64 {
        kn * 1852.0 / 3600.0 * l / 1.19e-6
    }
    #[test]
    fn schoenherr_matches_prototype() {
        assert!(close(schoenherr_cf(re(176.0, 14.0)) * 1e3, 1.519, 0.001));
        assert!(close(schoenherr_cf(re(230.0, 24.0)) * 1e3, 1.379, 0.001));
    }
    #[test]
    fn granville_direct_vectors_176m_14kn() {
        let r = re(176.0, 14.0);
        for (ks_um, pct) in [
            (30.0, 1.0),
            (100.0, 21.0),
            (300.0, 42.2),
            (1000.0, 72.4),
            (3000.0, 108.7),
            (10000.0, 162.7),
        ] {
            let got = 100.0 * delta_cf_ratio_direct(r, ks_um * 1e-6 / 176.0);
            assert!(close(got, pct, 0.15), "ks={ks_um}: {got}");
        }
    }
    #[test]
    fn granville_direct_vectors_other_ships() {
        assert!(close(
            100.0 * delta_cf_ratio_direct(re(124.0, 15.0), 10000e-6 / 124.0),
            172.6,
            0.15
        ));
        assert!(close(
            100.0 * delta_cf_ratio_direct(re(320.0, 15.0), 1000e-6 / 320.0),
            69.4,
            0.15
        ));
        assert!(close(
            100.0 * delta_cf_ratio_direct(re(100.0, 12.0), 300e-6 / 100.0),
            42.0,
            0.15
        ));
    }
    #[test]
    fn table_matches_direct_off_grid() {
        let t = GranvilleTable::new();
        let r = re(230.0, 16.0);
        assert!(close(
            100.0 * t.delta_cf_ratio(r, 150e-6 / 230.0),
            29.5,
            0.4
        ));
        assert!(close(
            100.0 * t.delta_cf_ratio(r, 3000e-6 / 230.0),
            108.8,
            0.6
        ));
    }
    #[test]
    fn table_clamps_outside_grid() {
        let t = GranvilleTable::new();
        let a = t.delta_cf_ratio(1e12, 1e-12);
        let b = t.delta_cf_ratio(1e10, 1e-8);
        assert!(a.is_finite() && close(a, b, 1e-12));
        assert!(t.delta_cf_ratio(1e9, 1.0).is_finite());
    }
    #[test]
    fn ittc_allowance_vectors() {
        let r = re(176.0, 14.0);
        assert!(close(
            ittc_roughness_allowance(150e-6 / 176.0, r) * 1e3,
            0.111,
            0.002
        ));
        assert!(close(
            ittc_roughness_allowance(30e-6 / 176.0, r) * 1e3,
            -0.062,
            0.002
        ));
    }
}
