//! Fixed-step solver for the guild model: integrating-factor (Lawson) RK4.
//!
//! The state obeys `dx/dt = −L·x + G(x)`. `L` holds the linear losses (mortality,
//! detachment, weathering), with mass optionally routed from one component into another
//! (live shells → dead shells). The flow of `L` is solved **exactly**. RK4 handles the
//! nonlinear gains `G` (settlement and growth) in the integrating-factor frame, so the scheme
//! is fourth order overall and exact when `G = 0`. A 3-hour freshwater half-life therefore
//! stays exact even with quarter-day steps.
//!
//! Steps are fixed: a segment of length Δ takes `max(1, ⌈Δ/0.25⌉)` equal steps. There are no
//! tolerance loops, so results are deterministic.
// TEMP(Task 7): remove once GuildModel drives the solver.
#![allow(dead_code)]

/// Largest step, in days.
pub(crate) const MAX_STEP_DAYS: f64 = 0.25;

/// Linear losses: component `j` decays at `rate[j]` per day, and the lost mass goes to
/// component `into[j]` (or leaves the system when `None`). A component that receives mass
/// must not itself route mass onward (one level of transfer, live → dead).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LinearLoss<const N: usize> {
    /// Loss rate per day of each component.
    pub rate: [f64; N],
    /// Where each component's lost mass goes.
    pub into: [Option<usize>; N],
}

/// Number of equal steps for a segment of `duration` days (at least 1).
pub(crate) fn step_count(duration: f64, max_step: f64) -> usize {
    ((duration / max_step).ceil() as usize).max(1)
}

/// Advances `x` over `duration` days with the default step bound.
pub(crate) fn advance<const N: usize>(
    x: &mut [f64; N],
    duration: f64,
    loss: &LinearLoss<N>,
    gain: impl Fn(&[f64; N]) -> [f64; N],
) {
    advance_with_step(x, duration, MAX_STEP_DAYS, loss, gain);
}

/// Advances `x` over `duration` days with steps of at most `max_step` days.
pub(crate) fn advance_with_step<const N: usize>(
    x: &mut [f64; N],
    duration: f64,
    max_step: f64,
    loss: &LinearLoss<N>,
    gain: impl Fn(&[f64; N]) -> [f64; N],
) {
    if duration == 0.0 {
        return;
    }
    let n = step_count(duration, max_step);
    let h = duration / n as f64;
    for _ in 0..n {
        let (xa, xb) = (flow(x, h / 2.0, loss), flow(x, h, loss));
        let k1 = gain(x);
        let k2 = gain(&axpy(&xa, h / 2.0, &flow(&k1, h / 2.0, loss)));
        let k3 = gain(&axpy(&xa, h / 2.0, &k2));
        let k4 = gain(&axpy(&xb, h, &flow(&k3, h / 2.0, loss)));
        let (f1, f23) = (flow(&k1, h, loss), flow(&add(&k2, &k3), h / 2.0, loss));
        for j in 0..N {
            let v = xb[j] + h / 6.0 * (f1[j] + 2.0 * f23[j] + k4[j]);
            // Clamp round-off only; a NaN passes through unchanged.
            x[j] = if v < 0.0 { 0.0 } else { v };
        }
    }
}

/// Exact solution of `dx/dt = −L·x` over `tau` days (a linear map, so it also applies to
/// derivative vectors).
fn flow<const N: usize>(x: &[f64; N], tau: f64, loss: &LinearLoss<N>) -> [f64; N] {
    let mut y = [0.0; N];
    for j in 0..N {
        y[j] += x[j] * (-loss.rate[j] * tau).exp();
        if let Some(k) = loss.into[j] {
            y[k] += x[j] * loss.rate[j] * transfer_kernel(loss.rate[j], loss.rate[k], tau);
        }
    }
    y
}

/// `(e^(−bτ) − e^(−aτ))/(a − b)`: the fraction of a source decaying at `a` that is found
/// in a target decaying at `b` after `tau`, per unit source rate. Limit `τ·e^(−bτ)` when a ≈ b.
fn transfer_kernel(a: f64, b: f64, tau: f64) -> f64 {
    let d = a - b;
    if (d * tau).abs() < 1e-8 {
        tau * (-b * tau).exp() * (1.0 - d * tau / 2.0)
    } else {
        -(-b * tau).exp() * (-d * tau).exp_m1() / d
    }
}

fn axpy<const N: usize>(x: &[f64; N], a: f64, y: &[f64; N]) -> [f64; N] {
    let mut z = *x;
    for j in 0..N {
        z[j] += a * y[j];
    }
    z
}

fn add<const N: usize>(x: &[f64; N], y: &[f64; N]) -> [f64; N] {
    axpy(x, 1.0, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rel(a: f64, b: f64) -> f64 {
        (a - b).abs() / a.abs().max(b.abs())
    }
    fn logistic(x: &[f64; 1]) -> [f64; 1] {
        [0.1 * x[0] * (1.0 - x[0])]
    }
    const NO_LOSS: LinearLoss<1> = LinearLoss {
        rate: [0.0],
        into: [None],
    };

    /// Linear losses are solved exactly: pure decay at 0.3/day over 10 days gives x0·e^(−3)
    /// (to round-off), whatever the step size.
    #[test]
    fn pure_decay_is_exact() {
        let mut x = [2.0];
        advance(
            &mut x,
            10.0,
            &LinearLoss {
                rate: [0.3],
                into: [None],
            },
            |_| [0.0],
        );
        assert!(rel(x[0], 2.0 * (-3.0f64).exp()) < 1e-13, "{}", x[0]);
    }

    /// Mass lost by a decaying component and routed into another (live → dead shells) is
    /// conserved exactly, and the target's own decay follows the closed form
    /// x1(T) = x1₀e^(−ωT) + x0₀·λ(e^(−ωT) − e^(−λT))/(λ − ω).
    #[test]
    fn transfer_conserves_mass_and_matches_closed_form() {
        let mut x = [3.0, 1.0];
        advance(
            &mut x,
            7.3,
            &LinearLoss {
                rate: [0.2, 0.0],
                into: [Some(1), None],
            },
            |_| [0.0; 2],
        );
        assert!((x[0] + x[1] - 4.0).abs() < 1e-12);
        assert!(rel(x[0], 3.0 * (-0.2f64 * 7.3).exp()) < 1e-13);
        let (l, w, t) = (0.2f64, 0.01f64, 30.0f64);
        let mut y = [3.0, 1.0];
        advance(
            &mut y,
            t,
            &LinearLoss {
                rate: [l, w],
                into: [Some(1), None],
            },
            |_| [0.0; 2],
        );
        let want = 1.0 * (-w * t).exp() + 3.0 * l * ((-w * t).exp() - (-l * t).exp()) / (l - w);
        assert!(rel(y[1], want) < 1e-12, "{} vs {want}", y[1]);
        let mut z = [3.0, 1.0]; // equal rates use the limit form t·e^(−λt)
        advance(
            &mut z,
            t,
            &LinearLoss {
                rate: [l, l],
                into: [Some(1), None],
            },
            |_| [0.0; 2],
        );
        assert!(rel(z[1], 1.0 * (-l * t).exp() + 3.0 * l * t * (-l * t).exp()) < 1e-12);
    }

    /// Nonlinear gains: logistic growth dx/dt = 0.1·x(1 − x) from 0.01 over 100 days matches
    /// the closed form 1/(1 + 99·e^(−10)) to 1e-8 with the default quarter-day step.
    #[test]
    fn logistic_matches_closed_form() {
        let mut x = [0.01];
        advance(&mut x, 100.0, &NO_LOSS, logistic);
        assert!((x[0] - 1.0 / (1.0 + 99.0 * (-10.0f64).exp())).abs() < 1e-8);
    }

    /// The scheme is fourth order when gains and losses act together: halving the step cuts
    /// the error by 12–20× (2⁴ = 16).
    #[test]
    fn fourth_order_convergence() {
        let loss = LinearLoss {
            rate: [0.05],
            into: [None],
        };
        let run = |h: f64| {
            let mut x = [0.01];
            advance_with_step(&mut x, 100.0, h, &loss, logistic);
            x[0]
        };
        let reference = run(1e-3);
        let ratio = (run(0.5) - reference).abs() / (run(0.25) - reference).abs();
        assert!((12.0..=20.0).contains(&ratio), "ratio {ratio}");
    }

    /// A ten-year stay as one segment equals the same time as 3,650 daily segments (to 1e-6),
    /// and a microsecond-scale segment still takes exactly one step. (Review focus 5.)
    #[test]
    fn long_and_tiny_segments() {
        let loss = LinearLoss {
            rate: [0.05],
            into: [None],
        };
        let mut one = [0.01];
        advance(&mut one, 3650.0, &loss, logistic);
        let mut daily = [0.01];
        for _ in 0..3650 {
            advance(&mut daily, 1.0, &loss, logistic);
        }
        assert!(rel(one[0], daily[0]) < 1e-6, "{} vs {}", one[0], daily[0]);
        assert_eq!(step_count(1e-6, MAX_STEP_DAYS), 1);
        assert_eq!(step_count(1.0, MAX_STEP_DAYS), 4);
        assert_eq!(step_count(1.01, MAX_STEP_DAYS), 5);
    }

    /// A zero-length segment changes nothing; the same input always gives bit-identical output.
    #[test]
    fn zero_duration_and_replay() {
        let mut x = [0.3];
        advance(
            &mut x,
            0.0,
            &LinearLoss {
                rate: [5.0],
                into: [None],
            },
            |_| [1.0],
        );
        assert_eq!(x, [0.3]);
        let run = || {
            let mut x = [0.01];
            advance(
                &mut x,
                123.456,
                &LinearLoss {
                    rate: [0.07],
                    into: [None],
                },
                logistic,
            );
            x
        };
        assert_eq!(run(), run());
    }

    /// Round-off below zero is clamped to zero, but a NaN is never hidden.
    #[test]
    fn clamps_round_off_but_not_nan() {
        let mut x = [1e-300];
        advance(&mut x, 1.0, &NO_LOSS, |_| [-1.0]);
        assert_eq!(x, [0.0]);
        let mut y = [f64::NAN];
        advance(&mut y, 1.0, &NO_LOSS, |_| [0.0]);
        assert!(y[0].is_nan());
    }
}
