//! Uzun et al. (2019) fouling-growth curves (equatorial site, self-polishing
//! coating) and their roughness mappings.
//!
//! **Age-clock rewrite.** Uzun gives fouling rating (FR) for slime and soft
//! fouling as a Gaussian rise to a plateau in calendar age, and calcareous
//! surface coverage (SC) as a sum of two logistics. A Gaussian in age is the
//! solution of a square-root-Gompertz ODE, `dFR/dA = (2/τ)·FR·√ln(max/FR)`,
//! so the curves can be driven by an *age clock* that runs faster or slower than
//! calendar time (temperature, salinity, settlement gate) while reproducing the
//! published curve exactly at the reference conditions. The inverses
//! ([`GaussianRating::age_for_rating`], [`CoverageLogistic::age_for_coverage`])
//! let reset maps act on ratings and coverage and be written back as clock values.
//!
//! **Open item 1 (spec).** The text-extracted Uzun Eq. 4 calcareous coefficient,
//! the evaluation point of the derivatives in Eqs. 6 and 8, and whether the fast
//! equatorial logistic term is applied are unconfirmed against the typeset
//! article. The constants here follow the spec's reading.
//!
//! The calcareous *rating* (Uzun class A) is deliberately not implemented: the
//! spec excludes it from the total rating; calcareous fouling enters through
//! coverage only.

use nautilus_units::Micrometres;

const BISECTION_STEPS: usize = 64;

/// Gaussian rise of a fouling rating to a plateau `max` at age `t0` (days).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussianRating {
    /// Plateau rating.
    pub max: f64,
    /// Age (days) at which the plateau is reached.
    pub t0: f64,
    /// Width (days) of the Gaussian rise.
    pub tau: f64,
}

impl GaussianRating {
    /// Rating at age `age_days`: `max·exp(−((A − t0)/τ)²)` for A < t0, else `max`.
    pub fn rating(&self, age_days: f64) -> f64 {
        if age_days < self.t0 {
            let z = (age_days - self.t0) / self.tau;
            self.max * (-z * z).exp()
        } else {
            self.max
        }
    }

    /// Rating at age 0 (the seed value of the Gaussian).
    pub fn seed(&self) -> f64 {
        self.rating(0.0)
    }

    /// Age at which the rating equals `fr`: `t0 − τ·√ln(max/fr)`; 0 if
    /// `fr ≤ seed`, `t0` if `fr ≥ max`.
    ///
    /// ```
    /// use nautilus_ocean::fouling::uzun::EQ_SLIME;
    /// let age = EQ_SLIME.age_for_rating(10.0);
    /// assert!((EQ_SLIME.rating(age) - 10.0).abs() < 1e-12);
    /// ```
    pub fn age_for_rating(&self, fr: f64) -> f64 {
        if fr.is_nan() || fr <= self.seed() {
            0.0
        } else if fr >= self.max {
            self.t0
        } else {
            self.t0 - self.tau * (self.max / fr).ln().sqrt()
        }
    }
}

/// Uzun equatorial slime rating curve.
pub const EQ_SLIME: GaussianRating = GaussianRating {
    max: 20.0,
    t0: 87.0,
    tau: 37.08,
};
/// Uzun equatorial soft-fouling rating curve.
pub const EQ_SOFT: GaussianRating = GaussianRating {
    max: 50.0,
    t0: 271.4,
    tau: 73.11,
};

/// Calcareous surface coverage (%) as the sum of a main and a fast logistic in
/// age; the asymptotes are passed in so they can depend on mean temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoverageLogistic {
    /// Offset of the main logistic.
    pub main_offset: f64,
    /// Rate (1/day) of the main logistic.
    pub main_rate: f64,
    /// Offset of the fast logistic.
    pub fast_offset: f64,
    /// Rate (1/day) of the fast logistic.
    pub fast_rate: f64,
}

impl CoverageLogistic {
    /// Coverage (%) at age `age_days`:
    /// `main_asym/(1+exp(main_offset − main_rate·A)) + fast_asym/(1+exp(fast_offset − fast_rate·A))`,
    /// clamped to [0, 100].
    pub fn coverage(&self, age_days: f64, main_asym: f64, fast_asym: f64) -> f64 {
        let main = main_asym / (1.0 + (self.main_offset - self.main_rate * age_days).exp());
        let fast = fast_asym / (1.0 + (self.fast_offset - self.fast_rate * age_days).exp());
        (main + fast).clamp(0.0, 100.0)
    }

    /// Age at which coverage equals `sc`, by 64-step bisection on
    /// `[0, max_age]` (coverage is monotone in age). Returns the bracket edge
    /// when `sc` lies outside the range reached on it.
    pub fn age_for_coverage(&self, sc: f64, main_asym: f64, fast_asym: f64, max_age: f64) -> f64 {
        let (mut lo, mut hi) = (0.0_f64, max_age);
        if sc.is_nan() || self.coverage(lo, main_asym, fast_asym) >= sc {
            return lo;
        }
        if self.coverage(hi, main_asym, fast_asym) <= sc {
            return hi;
        }
        for _ in 0..BISECTION_STEPS {
            let mid = 0.5 * (lo + hi);
            if self.coverage(mid, main_asym, fast_asym) < sc {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

/// Uzun equatorial calcareous coverage curve.
pub const EQ_CALCAREOUS: CoverageLogistic = CoverageLogistic {
    main_offset: 16.0,
    main_rate: 0.0407,
    fast_offset: 10.32,
    fast_rate: 0.7759,
};
/// Equatorial main-logistic asymptote (% coverage).
pub const EQ_MAIN_ASYM: f64 = 96.899;
/// Equatorial fast-logistic asymptote (% coverage).
pub const EQ_FAST_ASYM: f64 = 3.5;

/// Roughness from total fouling rating (Uzun Eq. 17):
/// `0.007143·FR² + 13.36·FR + 30` µm, FR clamped to [0, 70].
pub fn ks_from_rating(fr: f64) -> Micrometres {
    let fr = if fr.is_nan() {
        0.0
    } else {
        fr.clamp(0.0, 70.0)
    };
    Micrometres(0.007143 * fr * fr + 13.36 * fr + 30.0)
}

/// Roughness from calcareous coverage (Uzun Eq. 18): `None` for SC ≤ 5%;
/// `2.4669·SC² − 24.84·SC + 1065.7` µm for 5 < SC ≤ 50; `80·SC + 2000` above.
pub fn ks_from_coverage(sc: f64) -> Option<Micrometres> {
    if sc.is_nan() || sc <= 5.0 {
        None
    } else if sc <= 50.0 {
        Some(Micrometres(2.4669 * sc * sc - 24.84 * sc + 1065.7))
    } else {
        Some(Micrometres(80.0 * sc + 2000.0))
    }
}

/// Hull roughness: the larger of the rating and coverage roughness (coverage
/// overrides when present).
pub fn hull_ks(fr_total: f64, sc: f64) -> Micrometres {
    let from_rating = ks_from_rating(fr_total);
    match ks_from_coverage(sc) {
        Some(k) if k.value() > from_rating.value() => k,
        _ => from_rating,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }
    /// Uzun 2019 equatorial slime and soft-fouling curves at sample ages.
    /// Expected values: Uzun's Table 2 parameters evaluated independently.
    #[test]
    fn gaussian_ratings() {
        assert!(close(EQ_SLIME.seed(), 0.081326, 1e-5));
        assert!(close(EQ_SLIME.rating(50.0), 7.38937, 1e-4));
        assert_eq!(EQ_SLIME.rating(400.0), 20.0);
        assert!(close(EQ_SOFT.rating(200.0), 19.2643, 1e-3));
        assert!(close(EQ_SOFT.seed(), 5.178e-5, 1e-7));
    }
    /// Inverting the curves (rating → age), which reset maps rely on; edge cases below the
    /// seed and above the plateau.
    #[test]
    fn gaussian_inverse() {
        assert!(close(EQ_SLIME.age_for_rating(0.75), 19.8103, 1e-3));
        assert_eq!(EQ_SLIME.age_for_rating(1e-9), 0.0);
        assert_eq!(EQ_SLIME.age_for_rating(25.0), 87.0);
        let a = EQ_SOFT.age_for_rating(10.0);
        assert!(close(EQ_SOFT.rating(a), 10.0, 1e-9));
    }
    /// Uzun calcareous coverage curve: 3.5% at 100 days, ~52% at the 393-day midpoint, capped
    /// at 100%; inverse returns the right age. Expected values: Uzun Table 3.
    #[test]
    fn coverage_curve_and_inverse() {
        let c = |a| EQ_CALCAREOUS.coverage(a, EQ_MAIN_ASYM, EQ_FAST_ASYM);
        assert!(close(c(100.0), 3.50064, 1e-4));
        assert!(close(c(393.1204), 51.9495, 1e-3));
        assert_eq!(c(5000.0), 100.0); // 96.899 + 3.5 > 100 ⇒ clamped
        let a = EQ_CALCAREOUS.age_for_coverage(30.0, EQ_MAIN_ASYM, EQ_FAST_ASYM, 5000.0);
        assert!(close(c(a), 30.0, 1e-6));
    }
    /// Uzun Eqs. 17–18, converting fouling rating and coverage into an equivalent sand-grain
    /// roughness; coverage overrides rating when present. Expected values: the equations by hand.
    #[test]
    fn ks_mappings() {
        assert!(close(ks_from_rating(0.75).value(), 40.024, 1e-3));
        assert!(close(ks_from_rating(20.0).value(), 300.057, 1e-3));
        assert!(close(ks_from_rating(500.0).value(), 1000.2007, 1e-3)); // clamped at 70
        assert_eq!(ks_from_coverage(4.9), None);
        assert!(close(
            ks_from_coverage(25.0).unwrap().value(),
            1986.5125,
            1e-3
        ));
        assert!(close(ks_from_coverage(60.0).unwrap().value(), 6800.0, 1e-9));
        assert!(close(hull_ks(10.0, 25.0).value(), 1986.5125, 1e-3));
        assert!(close(
            hull_ks(10.0, 2.0).value(),
            ks_from_rating(10.0).value(),
            1e-12
        ));
    }
}
