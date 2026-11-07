//! Rate functions of the guild model: pure functions of the environment and the
//! parameters (spec Section 5.3), and [`SegmentRates`], the constant rates of one segment.

use nautilus_units::{Celsius, Knots, MilligramsPerCubicMetre, Psu};

use super::environment::{SalinityResponse, Segment};
use super::params::{CALCAREOUS, GuildParams};
pub use super::params::{Detachment, GuildWindow};

/// Thermal performance curve: `ref_rate·q10^((T − rate_ref_temp)/10)` up to `t_opt`, then
/// linear from that value to 0 at `t_max`, and 0 above.
pub fn thermal_rate(ref_rate: f64, q10: f64, w: &GuildWindow, t: Celsius) -> f64 {
    let q = |t: f64| ref_rate * q10.powf((t - w.rate_ref_temp.value()) / 10.0);
    let (t, opt, max) = (t.value(), w.t_opt.value(), w.t_max.value());
    if t <= opt {
        q(t)
    } else if t < max {
        q(opt) * (max - t) / (max - opt)
    } else {
        0.0
    }
}

/// Breeding factor: 1 on `[breed_low, breed_high]`, falling linearly to 0 over `ramp`
/// outside it.
pub fn breeding(w: &GuildWindow, ramp: Celsius, t: Celsius) -> f64 {
    let outside = (w.breed_low.value() - t.value()).max(t.value() - w.breed_high.value());
    if outside <= 0.0 {
        1.0
    } else if outside < ramp.value() {
        1.0 - outside / ramp.value()
    } else {
        0.0
    }
}

/// Thermal mortality per day: 0 on `[lethal_low, lethal_high]`, rising linearly to `max`
/// over `ramp` outside it.
pub fn thermal_mortality(w: &GuildWindow, max: f64, ramp: Celsius, t: Celsius) -> f64 {
    let outside = (w.lethal_low.value() - t.value()).max(t.value() - w.lethal_high.value());
    if outside <= 0.0 {
        0.0
    } else {
        max * (outside / ramp.value()).min(1.0)
    }
}

/// Freshwater mortality per day: `rate` below `threshold`, else 0.
pub fn fresh_mortality(rate: f64, threshold: Psu, s: Psu) -> f64 {
    if s.value() < threshold.value() {
        rate
    } else {
        0.0
    }
}

/// Brackish mortality of stenohaline (slime and soft) fouling per day: `max` below the
/// salinity at which their growth stops, else 0.
pub fn brackish_mortality(max: f64, stenohaline: &SalinityResponse, s: Psu) -> f64 {
    if s.value() < stenohaline.zero_below.value() {
        max
    } else {
        0.0
    }
}

/// Food multiplier `1 + α·chl/(chl + K)`; exactly 1 when chlorophyll is unknown.
pub fn food(
    alpha: f64,
    half_saturation: MilligramsPerCubicMetre,
    chl: Option<MilligramsPerCubicMetre>,
) -> f64 {
    match chl {
        None => 1.0,
        Some(c) => 1.0 + alpha * c.value() / (c.value() + half_saturation.value()),
    }
}

/// Detachment per day: `rate_per_knot·max(0, v − threshold)`.
pub fn detachment(d: &Detachment, v: Knots) -> f64 {
    d.rate_per_knot * (v.value() - d.threshold.value()).max(0.0)
}

/// Calcareous capacity (% coverage): the line through (17 °C, `capacity_17`) and
/// (25 °C, `capacity_25`), clamped to `[capacity_min, 100]`.
pub fn calcareous_capacity(p: &GuildParams, t: Celsius) -> f64 {
    let k = p.capacity_17 + (p.capacity_25 - p.capacity_17) * (t.value() - 17.0) / 8.0;
    k.clamp(p.capacity_min, 100.0)
}

/// Constant rates of one segment, indexed `[class][guild]` (class: slime, soft,
/// calcareous; guild: warm, cold).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SegmentRates {
    /// Settlement: gate × thermal settlement rate × breeding × salinity × food.
    pub settle: [[f64; 2]; 3],
    /// Growth: growth gate × thermal growth rate × salinity × food.
    pub grow: [[f64; 2]; 3],
    /// Linear losses per day: thermal + freshwater (+ brackish and detachment for slime/soft).
    pub loss: [[f64; 2]; 3],
    /// Calcareous capacity K_h(T), % coverage.
    pub capacity: f64,
    /// Dead-shell weathering per day.
    pub weathering: f64,
}

impl SegmentRates {
    /// Evaluates every rate for the segment's environment.
    pub(crate) fn new(seg: &Segment, p: &GuildParams) -> SegmentRates {
        let (t, s, v) = (seg.sst(), seg.salinity(), seg.speed_through_water());
        let x = p.gate.settlement(v) * p.gate.propagule(seg.port(), seg.distance_to_coast());
        let gamma = p.gate.gate(seg);
        let phi = food(p.food_alpha, p.food_half_saturation, seg.chlorophyll());
        let mut r = SegmentRates {
            settle: [[0.0; 2]; 3],
            grow: [[0.0; 2]; 3],
            loss: [[0.0; 2]; 3],
            capacity: calcareous_capacity(p, t),
            weathering: p.weathering_per_day,
        };
        for c in 0..3 {
            let rates = &p.classes[c];
            let hard = c == CALCAREOUS;
            let sal = if hard {
                p.hard_salinity.factor(s)
            } else {
                p.soft_salinity.factor(s)
            };
            let fresh_rate = if hard {
                p.fresh_mortality_hard
            } else {
                p.fresh_mortality_soft
            };
            let extra = if hard {
                0.0
            } else {
                brackish_mortality(p.mortality_max, &p.soft_salinity, s)
                    + detachment(&p.detachment[c], v)
            };
            for (g, w) in p.guilds.iter().enumerate() {
                r.settle[c][g] = x
                    * thermal_rate(rates.settle_ref, rates.settle_q10, w, t)
                    * breeding(w, p.breeding_ramp, t)
                    * sal
                    * phi;
                r.grow[c][g] =
                    gamma * thermal_rate(rates.growth_ref, rates.growth_q10, w, t) * sal * phi;
                r.loss[c][g] = thermal_mortality(w, p.mortality_max, p.mortality_ramp, t)
                    + fresh_mortality(fresh_rate, p.fresh_threshold, s)
                    + extra;
            }
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fouling::environment::{PortStatus, Segment};
    use crate::fouling::params::{COLD, SOFT, WARM};
    use nautilus_units::{Celsius, Days, Kilometres, Knots, MilligramsPerCubicMetre, Psu};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }
    fn warm() -> GuildWindow {
        GuildParams::default().guilds[WARM]
    }
    fn cold() -> GuildWindow {
        GuildParams::default().guilds[COLD]
    }

    /// Thermal performance curve: Q10 rise to the optimum, linear fall to zero at T_max,
    /// zero above; the cold guild peaks at the same rate at its own optimum.
    /// Expected values: the formula by hand (0.05·2^((15−25)/10) = 0.025, etc.).
    #[test]
    fn thermal_rate_curve() {
        let w = warm();
        let r = |t: f64| thermal_rate(0.05, 2.0, &w, Celsius(t));
        assert!(close(r(15.0), 0.025, 1e-15));
        assert!(close(r(29.0), 0.05 * 2f64.powf(0.4), 1e-15));
        assert!(close(r(32.0), 0.5 * r(29.0), 1e-15));
        assert_eq!((r(35.0), r(40.0)), (0.0, 0.0));
        assert!(close(
            thermal_rate(0.05, 2.0, &cold(), Celsius(13.0)),
            r(29.0),
            1e-15
        ));
    }

    /// Breeding range with 1 °C ramps: warm guild breeds at ≥ 15 °C, half at 14.5 °C, not at
    /// 13.9 °C, and half again at 35.5 °C.
    #[test]
    fn breeding_edges() {
        let b = |t: f64| breeding(&warm(), Celsius(1.0), Celsius(t));
        assert!(close(b(14.5), 0.5, 1e-12));
        assert_eq!((b(15.0), b(13.9), b(25.0)), (1.0, 0.0, 1.0));
        assert!(close(b(35.5), 0.5, 1e-12));
    }

    /// Thermal mortality: zero inside the warm guild's 10–35 °C window, rising over 3 °C to
    /// 0.1/day outside it; finite exactly on the window edges.
    #[test]
    fn thermal_mortality_ramps() {
        let m = |t: f64| thermal_mortality(&warm(), 0.1, Celsius(3.0), Celsius(t));
        assert_eq!(m(20.0), 0.0);
        assert!(close(m(8.5), 0.05, 1e-15));
        assert_eq!((m(6.0), m(-2.0)), (0.1, 0.1));
        assert!(close(m(36.5), 0.05, 1e-15));
        assert_eq!((m(10.0), m(35.0)), (0.0, 0.0));
    }

    /// Freshwater kills half of soft fouling per 3 h (1/16 left after 12 h) below 2 PSU; brackish
    /// water below the stenohaline 10 PSU limit kills it at 0.1/day.
    #[test]
    fn fresh_and_brackish() {
        let rate = std::f64::consts::LN_2 / 0.125;
        let f = fresh_mortality(rate, Psu(2.0), Psu(0.3));
        assert!(close((-f * 0.125).exp(), 0.5, 1e-12));
        assert!(close((-f * 0.5).exp(), 1.0 / 16.0, 1e-12));
        assert_eq!(fresh_mortality(rate, Psu(2.0), Psu(2.0)), 0.0);
        let steno = GuildParams::default().soft_salinity;
        assert_eq!(brackish_mortality(0.1, &steno, Psu(6.0)), 0.1);
        assert_eq!(brackish_mortality(0.1, &steno, Psu(12.0)), 0.0);
    }

    /// Food multiplier 1 + α·chl/(chl + K): exactly 1 when α = 0 or chlorophyll is unknown,
    /// 1.25 at chl = K with α = 0.5, and bounded below 1 + α for absurd concentrations.
    #[test]
    fn food_factor() {
        let k = MilligramsPerCubicMetre(1.0);
        assert_eq!(food(0.0, k, Some(MilligramsPerCubicMetre(5.0))), 1.0);
        assert_eq!(food(0.5, k, None), 1.0);
        assert!(close(food(0.5, k, Some(k)), 1.25, 1e-15));
        let extreme = food(0.5, k, Some(MilligramsPerCubicMetre(1e6)));
        assert!(extreme > 1.4999 && extreme < 1.5);
    }

    /// Detachment acts only above its speed threshold: 0.05/day/kn above 10 kn gives 0 at 9 kn
    /// and 0.2/day at 14 kn.
    #[test]
    fn detachment_threshold() {
        let d = Detachment {
            rate_per_knot: 0.05,
            threshold: Knots(10.0),
        };
        assert_eq!(detachment(&d, Knots(9.0)), 0.0);
        assert!(close(detachment(&d, Knots(14.0)), 0.2, 1e-15));
    }

    /// Calcareous capacity: the line through (17 °C, c17) and (25 °C, c25), clamped to
    /// [capacity_min, 100].
    #[test]
    fn capacity_line() {
        let p = GuildParams {
            capacity_17: 50.0,
            capacity_25: 90.0,
            capacity_min: 10.0,
            ..GuildParams::default()
        };
        let k = |t: f64| calcareous_capacity(&p, Celsius(t));
        assert_eq!((k(17.0), k(25.0)), (50.0, 90.0));
        assert!(close(k(21.0), 70.0, 1e-12));
        assert_eq!((k(5.0), k(35.0)), (10.0, 100.0));
    }

    /// Per-segment rates combine gate, temperature, breeding, salinity and food: at a 25 °C
    /// berth in sea water the warm soft-fouling settlement rate equals its reference value,
    /// and the cold guild (1 °C above its 24 °C limit) dies at 0.1/3 per day.
    #[test]
    fn segment_rates_assemble() {
        let p = GuildParams::default();
        let seg = Segment::new(
            Days(1.0),
            Knots(0.0),
            Celsius(25.0),
            Psu(35.0),
            PortStatus::Berth,
            Kilometres(0.0),
        )
        .unwrap();
        let r = SegmentRates::new(&seg, &p);
        assert!(close(
            r.settle[SOFT][WARM],
            p.classes[SOFT].settle_ref,
            1e-15
        ));
        assert!(close(r.grow[SOFT][WARM], p.classes[SOFT].growth_ref, 1e-15));
        assert!(close(r.loss[SOFT][COLD], 0.1 / 3.0, 1e-15));
        assert_eq!(r.settle[SOFT][COLD], 0.0); // cold guild does not breed at 25 °C
        assert_eq!(r.capacity, p.capacity_25);
        assert_eq!(r.weathering, p.weathering_per_day);
    }
}
