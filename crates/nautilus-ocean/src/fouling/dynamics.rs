//! Right-hand side of the guild model (spec Section 5.2): gains (settlement and growth)
//! and linear losses (mortality, detachment, weathering) for the state vector.
// TEMP(Task 7): remove once GuildModel uses the dynamics.
#![allow(dead_code)]

use super::integrate::LinearLoss;
pub(crate) use super::rates::SegmentRates;

/// Length of the state vector.
pub(crate) const N: usize = 10;

/// Positions in the state vector (`[warm, cold]` pairs for each class).
pub(crate) mod ix {
    /// Slime rating per guild.
    pub const SLIME: [usize; 2] = [0, 1];
    /// Soft-fouling rating per guild.
    pub const SOFT: [usize; 2] = [2, 3];
    /// Live calcareous coverage on the hull per guild (%).
    pub const CALC: [usize; 2] = [4, 5];
    /// Dead shell coverage on the hull (%).
    pub const DEAD: usize = 6;
    /// Live calcareous coverage on the propeller per guild (%).
    pub const PROP_CALC: [usize; 2] = [7, 8];
    /// Dead shell coverage on the propeller (%).
    pub const PROP_DEAD: usize = 9;
}

/// Rating capacities of slime and soft fouling (Uzun maxima).
pub(crate) const RATING_CAPACITY: [f64; 2] = [20.0, 50.0];

/// Gain terms dX/dt (settlement + growth). Never negative for a valid state.
pub(crate) fn gains(x: &[f64; N], r: &SegmentRates) -> [f64; N] {
    let mut d = [0.0; N];
    for (c, idx) in [ix::SLIME, ix::SOFT].into_iter().enumerate() {
        let k = RATING_CAPACITY[c];
        let total = x[idx[0]] + x[idx[1]];
        let ln = (k / total).ln();
        for g in 0..2 {
            let settle = if total < k {
                r.settle[c][g] * (1.0 - total / k)
            } else {
                0.0
            };
            // Growth is 0 at F = 0 (its limit) and whenever ln(K/F) ≤ 0 (round-off above K).
            let grow = if total > 0.0 && ln > 0.0 {
                r.grow[c][g] * x[idx[g]] * ln.sqrt()
            } else {
                0.0
            };
            d[idx[g]] = settle + grow;
        }
    }
    for (live, dead) in [(ix::CALC, ix::DEAD), (ix::PROP_CALC, ix::PROP_DEAD)] {
        let free = r.capacity - x[live[0]] - x[live[1]] - x[dead];
        if free > 0.0 {
            for g in 0..2 {
                d[live[g]] = (r.settle[2][g] + r.grow[2][g] * x[live[g]]) * free / r.capacity;
            }
        }
    }
    d
}

/// Linear losses: per-guild mortality and detachment; live shells go to the matching
/// dead-shell pool; dead shells weather away.
pub(crate) fn losses(r: &SegmentRates) -> LinearLoss<N> {
    let mut l = LinearLoss {
        rate: [0.0; N],
        into: [None; N],
    };
    for g in 0..2 {
        l.rate[ix::SLIME[g]] = r.loss[0][g];
        l.rate[ix::SOFT[g]] = r.loss[1][g];
        l.rate[ix::CALC[g]] = r.loss[2][g];
        l.into[ix::CALC[g]] = Some(ix::DEAD);
        l.rate[ix::PROP_CALC[g]] = r.loss[2][g];
        l.into[ix::PROP_CALC[g]] = Some(ix::PROP_DEAD);
    }
    l.rate[ix::DEAD] = r.weathering;
    l.rate[ix::PROP_DEAD] = r.weathering;
    l
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rates(settle: f64, grow: f64, capacity: f64) -> SegmentRates {
        SegmentRates {
            settle: [[settle; 2]; 3],
            grow: [[grow; 2]; 3],
            loss: [[0.02, 0.03]; 3],
            capacity,
            weathering: 0.001,
        }
    }

    /// A clean hull gains nothing when no larvae settle: growth needs existing fouling.
    #[test]
    fn clean_hull_without_settlement_stays_clean() {
        assert_eq!(gains(&[0.0; N], &rates(0.0, 0.1, 100.0)), [0.0; N]);
    }

    /// Slime and soft fouling stop gaining at their rating capacity (20 and 50): no
    /// settlement into a full class and no growth (√ln(K/F) = 0).
    #[test]
    fn growth_and_settlement_stop_at_capacity() {
        let mut x = [0.0; N];
        x[ix::SOFT[0]] = 30.0;
        x[ix::SOFT[1]] = 20.0;
        x[ix::SLIME[0]] = 20.0;
        let g = gains(&x, &rates(0.5, 0.1, 100.0));
        for i in ix::SOFT.iter().chain(&ix::SLIME) {
            assert_eq!(g[*i], 0.0, "index {i}");
        }
    }

    /// A rating pushed just above capacity by round-off must not make √ln(K/F) NaN: all gains
    /// stay finite and non-negative. (Review focus 1.)
    #[test]
    fn round_off_above_capacity_is_finite() {
        let mut x = [0.0; N];
        x[ix::SOFT[0]] = 50.0 * (1.0 + 1e-12);
        let g = gains(&x, &rates(0.5, 0.1, 100.0));
        assert!(g.iter().all(|v| v.is_finite() && *v >= 0.0), "{g:?}");
    }

    /// When the temperature-dependent capacity falls below the coverage already present
    /// (warm hull in cold water), barnacles gain nothing and nothing goes negative.
    /// (Review focus 2.)
    #[test]
    fn no_gain_when_capacity_below_coverage() {
        let mut x = [0.0; N];
        x[ix::CALC[0]] = 30.0;
        x[ix::PROP_CALC[1]] = 15.0;
        let g = gains(&x, &rates(0.5, 0.1, 10.0));
        for i in ix::CALC.iter().chain(&ix::PROP_CALC) {
            assert_eq!(g[*i], 0.0, "index {i}");
        }
    }

    /// Losses: each guild uses its own rate; dead live shells (hull and propeller) go to the
    /// matching dead-shell pool, and dead shells weather away.
    #[test]
    fn losses_route_live_shells_to_dead() {
        let l = losses(&rates(0.0, 0.0, 100.0));
        for g in 0..2 {
            assert_eq!(l.rate[ix::SOFT[g]], [0.02, 0.03][g]);
            assert_eq!(l.into[ix::SOFT[g]], None);
            assert_eq!(l.into[ix::CALC[g]], Some(ix::DEAD));
            assert_eq!(l.into[ix::PROP_CALC[g]], Some(ix::PROP_DEAD));
        }
        assert_eq!((l.rate[ix::DEAD], l.rate[ix::PROP_DEAD]), (0.001, 0.001));
        assert_eq!((l.into[ix::DEAD], l.into[ix::PROP_DEAD]), (None, None));
    }

    /// Gains are never negative for any valid state and rates (1,000 random cases), so
    /// fouling only decreases through the explicit loss terms.
    #[test]
    fn gains_never_negative() {
        let mut s = 61u64;
        let mut u = |lo: f64, hi: f64| {
            s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = s;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            lo + (hi - lo) * ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..1000 {
            let mut x = [0.0; N];
            for (i, v) in x.iter_mut().enumerate() {
                *v = match i {
                    0 | 1 => u(0.0, 10.0),
                    2 | 3 => u(0.0, 25.0),
                    _ => u(0.0, 33.0),
                };
            }
            let r = rates(u(0.0, 1.0), u(0.0, 0.2), u(10.0, 100.0));
            assert!(gains(&x, &r).iter().all(|v| *v >= 0.0 && v.is_finite()));
        }
    }
}
