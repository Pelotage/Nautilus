//! Segment forcing and the environmental response functions that set how fast
//! the exposure clocks run.

use nautilus_units::{Celsius, Days, Kilometres, Knots, Psu};

use super::FoulingError;

/// Where the ship is during a segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortStatus {
    /// Under way or drifting at sea.
    AtSea,
    /// At anchor (waiting area, roads).
    Anchorage,
    /// Alongside a berth.
    Berth,
}

/// A stretch of time with constant forcing (speed, sea conditions, location).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    duration: Days,
    speed_through_water: Knots,
    sst: Celsius,
    salinity: Psu,
    port: PortStatus,
    distance_to_coast: Kilometres,
}

impl Segment {
    /// Builds a segment. Rejects non-finite values and negative duration,
    /// speed, salinity or distance with [`FoulingError::InvalidSegment`].
    /// A zero-duration segment is valid and is a no-op for every model.
    ///
    /// ```
    /// use nautilus_ocean::fouling::environment::{PortStatus, Segment};
    /// use nautilus_ocean::fouling::FoulingError;
    /// use nautilus_units::{Celsius, Days, Kilometres, Knots, Psu};
    ///
    /// let ok = Segment::new(Days(1.5), Knots(12.0), Celsius(18.0), Psu(34.0), PortStatus::AtSea, Kilometres(80.0));
    /// assert!(ok.is_ok());
    /// let bad = Segment::new(Days(-1.0), Knots(12.0), Celsius(18.0), Psu(34.0), PortStatus::AtSea, Kilometres(80.0));
    /// assert_eq!(bad, Err(FoulingError::InvalidSegment("duration")));
    /// ```
    pub fn new(
        duration: Days,
        speed_through_water: Knots,
        sst: Celsius,
        salinity: Psu,
        port: PortStatus,
        distance_to_coast: Kilometres,
    ) -> Result<Segment, FoulingError> {
        fn non_negative(v: f64, what: &'static str) -> Result<(), FoulingError> {
            if v.is_finite() && v >= 0.0 {
                Ok(())
            } else {
                Err(FoulingError::InvalidSegment(what))
            }
        }
        non_negative(duration.value(), "duration")?;
        non_negative(speed_through_water.value(), "speed_through_water")?;
        if !sst.value().is_finite() {
            return Err(FoulingError::InvalidSegment("sst"));
        }
        non_negative(salinity.value(), "salinity")?;
        non_negative(distance_to_coast.value(), "distance_to_coast")?;
        Ok(Segment {
            duration,
            speed_through_water,
            sst,
            salinity,
            port,
            distance_to_coast,
        })
    }

    /// Segment duration.
    pub fn duration(&self) -> Days {
        self.duration
    }
    /// Speed through water.
    pub fn speed_through_water(&self) -> Knots {
        self.speed_through_water
    }
    /// Sea-surface temperature.
    pub fn sst(&self) -> Celsius {
        self.sst
    }
    /// Sea-surface salinity.
    pub fn salinity(&self) -> Psu {
        self.salinity
    }
    /// Port status.
    pub fn port(&self) -> PortStatus {
        self.port
    }
    /// Distance to the nearest coast.
    pub fn distance_to_coast(&self) -> Kilometres {
        self.distance_to_coast
    }
}

/// Q10 temperature law for a growth-clock rate, relative to `t_ref`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TemperatureResponse {
    /// Rate multiplier per +10 °C.
    pub q10: f64,
    /// Reference temperature where the factor is 1 (25 °C, Uzun equatorial site).
    pub t_ref: Celsius,
    /// Below this temperature the factor drops to `floor`. Assumption — calibrate.
    pub t_min: Celsius,
    /// Minimum factor (cold dormancy). Assumption — calibrate.
    pub floor: f64,
    /// Temperatures above this are clamped to it (thermal optimum). Assumption — calibrate.
    pub t_max: Celsius,
}

impl Default for TemperatureResponse {
    fn default() -> Self {
        TemperatureResponse {
            q10: 2.0,
            t_ref: Celsius(25.0),
            t_min: Celsius(12.0),
            floor: 0.1,
            t_max: Celsius(32.0),
        }
    }
}

impl TemperatureResponse {
    /// Rate factor at temperature `t`: `max(floor, q10^((min(t, t_max) − t_ref)/10))`,
    /// or `floor` below `t_min`.
    pub fn factor(&self, t: Celsius) -> f64 {
        let t = t.value().min(self.t_max.value());
        if t < self.t_min.value() {
            return self.floor;
        }
        self.q10
            .powf((t - self.t_ref.value()) / 10.0)
            .max(self.floor)
    }
}

/// Piecewise-linear salinity response: 0 below `zero_below`, 0.5 at `half_at`,
/// 1 at and above `full_at`. Assumption — calibrate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SalinityResponse {
    /// Factor is 0 at and below this salinity.
    pub zero_below: Psu,
    /// Factor is 0.5 at this salinity.
    pub half_at: Psu,
    /// Factor is 1 at and above this salinity.
    pub full_at: Psu,
}

impl SalinityResponse {
    /// Tolerant of brackish water (hard fouling: barnacles, tubeworms). Assumption — calibrate.
    pub const EURYHALINE: SalinityResponse = SalinityResponse {
        zero_below: Psu(2.0),
        half_at: Psu(5.0),
        full_at: Psu(20.0),
    };
    /// Drops out below 10–25 PSU (soft fouling and slime). Assumption — calibrate.
    pub const STENOHALINE: SalinityResponse = SalinityResponse {
        zero_below: Psu(10.0),
        half_at: Psu(17.5),
        full_at: Psu(25.0),
    };

    /// Rate factor at salinity `s`, in [0, 1].
    pub fn factor(&self, s: Psu) -> f64 {
        let (s, z, h, f) = (
            s.value(),
            self.zero_below.value(),
            self.half_at.value(),
            self.full_at.value(),
        );
        if s <= z {
            0.0
        } else if s < h {
            0.5 * (s - z) / (h - z)
        } else if s < f {
            0.5 + 0.5 * (s - h) / (f - h)
        } else {
            1.0
        }
    }
}

/// Settlement and propagule-supply gate: how much of a segment counts as
/// fouling exposure, from speed through water and distance to coast.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GateParams {
    /// Speed at which settlement is halved (logistic midpoint).
    pub settle_v50: Knots,
    /// Logistic width of the settlement drop-off.
    pub settle_width: Knots,
    /// At sea within this distance of the coast, propagule supply is full. Assumption — calibrate.
    pub near_coast: Kilometres,
    /// Beyond this distance, propagule supply is `offshore_floor`. Assumption — calibrate.
    pub far_coast: Kilometres,
    /// Offshore propagule supply. Assumption — calibrate.
    pub offshore_floor: f64,
    /// Fraction of non-gated time that still counts as exposure (δ). 0 is
    /// Uzun-compatible (growth only when idle). Assumption — calibrate.
    pub sailing_fraction: f64,
}

impl Default for GateParams {
    fn default() -> Self {
        GateParams {
            settle_v50: Knots(1.25),
            settle_width: Knots(0.2547),
            near_coast: Kilometres(15.0),
            far_coast: Kilometres(100.0),
            offshore_floor: 0.02,
            sailing_fraction: 0.0,
        }
    }
}

impl GateParams {
    /// Normalised settlement factor `A(v) = L(v)/L(0)` with
    /// `L(v) = 1/(1 + exp((v − v50)/width))`, so `A(0) = 1`.
    pub fn settlement(&self, v: Knots) -> f64 {
        let l = |v: f64| {
            1.0 / (1.0 + ((v - self.settle_v50.value()) / self.settle_width.value()).exp())
        };
        l(v.value()) / l(0.0)
    }

    /// Propagule supply: 1 at berth or anchorage; at sea 1 within `near_coast`,
    /// linear to `offshore_floor` at `far_coast`, `offshore_floor` beyond.
    pub fn propagule(&self, port: PortStatus, d: Kilometres) -> f64 {
        match port {
            PortStatus::Berth | PortStatus::Anchorage => 1.0,
            PortStatus::AtSea => {
                let (d, near, far) = (d.value(), self.near_coast.value(), self.far_coast.value());
                if d <= near {
                    1.0
                } else if d >= far {
                    self.offshore_floor
                } else {
                    1.0 + (self.offshore_floor - 1.0) * (d - near) / (far - near)
                }
            }
        }
    }

    /// Exposure gate for a segment: `x + sailing_fraction·(1 − x)` with
    /// `x = settlement(v)·propagule(port, d)`.
    pub fn gate(&self, seg: &Segment) -> f64 {
        let x = self.settlement(seg.speed_through_water())
            * self.propagule(seg.port(), seg.distance_to_coast());
        x + self.sailing_fraction * (1.0 - x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fouling::FoulingError;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    /// Temperature response: doubles per +10 °C (Q10 = 2), equals 1 at the 25 °C reference,
    /// stops rising above 32 °C, and drops to the 0.1 floor in polar water.
    /// Expected values: the Q10 formula worked by hand (2^(-1) = 0.5 at 15 °C).
    #[test]
    fn temperature_q10_and_clamps() {
        let t = TemperatureResponse::default();
        assert!(close(t.factor(Celsius(25.0)), 1.0, 1e-12));
        assert!(close(t.factor(Celsius(15.0)), 0.5, 1e-12));
        assert!(close(
            t.factor(Celsius(35.0)),
            t.factor(Celsius(32.0)),
            1e-12
        )); // warm clamp
        assert_eq!(t.factor(Celsius(-2.0)), 0.1); // cold floor
        assert!(t.factor(Celsius(12.0)) >= 0.1);
    }
    /// Salinity response: full growth in sea water, half at 5 PSU, none in fresh water, and
    /// soft fouling (stenohaline) already zero at 10 PSU.
    /// Expected values: the piecewise-linear ramp worked by hand.
    #[test]
    fn salinity_ramps() {
        let s = SalinityResponse::EURYHALINE;
        assert_eq!(s.factor(Psu(35.0)), 1.0);
        assert!(close(s.factor(Psu(10.0)), 0.5 + 0.5 * 5.0 / 15.0, 1e-12));
        assert!(close(s.factor(Psu(5.0)), 0.5, 1e-12));
        assert_eq!(s.factor(Psu(0.3)), 0.0);
        assert_eq!(SalinityResponse::STENOHALINE.factor(Psu(10.0)), 0.0);
    }
    /// Settlement vs speed: larvae settle freely when the ship is still, ~96% at 0.5 kn,
    /// ~5% at 2 kn, essentially never at 14 kn service speed.
    /// Expected values: the logistic formula evaluated independently (spec §settlement gate).
    #[test]
    fn settlement_gate_values() {
        let g = GateParams::default();
        assert!(close(g.settlement(Knots(0.0)), 1.0, 1e-12));
        assert!(close(g.settlement(Knots(0.5)), 0.95703, 1e-4));
        assert!(close(g.settlement(Knots(2.0)), 0.05036, 1e-4));
        assert!(g.settlement(Knots(14.0)) < 1e-20);
    }
    /// Larval supply vs location: full in port and near the coast, falling linearly to a 2%
    /// floor beyond 100 km offshore; a berthed, stationary ship gets the full gate.
    /// Expected values: linear interpolation worked by hand (57.5 km is halfway → 0.51).
    #[test]
    fn propagule_and_gate() {
        let g = GateParams::default();
        assert_eq!(g.propagule(PortStatus::Berth, Kilometres(500.0)), 1.0);
        assert_eq!(g.propagule(PortStatus::AtSea, Kilometres(10.0)), 1.0);
        assert!(close(
            g.propagule(PortStatus::AtSea, Kilometres(57.5)),
            0.51,
            1e-12
        ));
        assert_eq!(g.propagule(PortStatus::AtSea, Kilometres(300.0)), 0.02);
        let berth = Segment::new(
            Days(1.0),
            Knots(0.0),
            Celsius(25.0),
            Psu(35.0),
            PortStatus::Berth,
            Kilometres(0.0),
        )
        .unwrap();
        assert!(close(g.gate(&berth), 1.0, 1e-12));
    }
    /// Input checking: zero-length segments are allowed, but negative, NaN or infinite
    /// duration, speed, temperature or salinity are rejected instead of silently corrupting state.
    #[test]
    fn segment_validation() {
        let ok = |d: f64| {
            Segment::new(
                Days(d),
                Knots(0.0),
                Celsius(20.0),
                Psu(35.0),
                PortStatus::Berth,
                Kilometres(0.0),
            )
        };
        assert!(ok(0.0).is_ok());
        assert!(matches!(ok(-1.0), Err(FoulingError::InvalidSegment(_))));
        assert!(matches!(ok(f64::NAN), Err(FoulingError::InvalidSegment(_))));
        assert!(
            Segment::new(
                Days(1.0),
                Knots(-1.0),
                Celsius(20.0),
                Psu(35.0),
                PortStatus::AtSea,
                Kilometres(0.0)
            )
            .is_err()
        );
        assert!(
            Segment::new(
                Days(1.0),
                Knots(0.0),
                Celsius(f64::INFINITY),
                Psu(35.0),
                PortStatus::AtSea,
                Kilometres(0.0)
            )
            .is_err()
        );
        assert!(
            Segment::new(
                Days(1.0),
                Knots(0.0),
                Celsius(20.0),
                Psu(-0.1),
                PortStatus::AtSea,
                Kilometres(0.0)
            )
            .is_err()
        );
    }
    /// Growth while sailing (δ, spec open item 8): with δ = 0.2 a ship at speed far offshore
    /// still gets 20% exposure, a berthed ship stays at 100%, and a stationary ship 57.5 km
    /// out (x = 0.51) gets 0.51 + 0.2 × 0.49. Expected values: gate = x + δ·(1 − x) by hand.
    #[test]
    fn sailing_fraction_blends_toward_full_exposure() {
        let g = GateParams {
            sailing_fraction: 0.2,
            ..GateParams::default()
        };
        let seg = |v: f64, p: PortStatus, d: f64| {
            Segment::new(
                Days(1.0),
                Knots(v),
                Celsius(25.0),
                Psu(35.0),
                p,
                Kilometres(d),
            )
            .unwrap()
        };
        let x_sea = g.settlement(Knots(14.0)) * 0.02;
        assert!(close(
            g.gate(&seg(14.0, PortStatus::AtSea, 300.0)),
            x_sea + 0.2 * (1.0 - x_sea),
            1e-15
        ));
        assert!(close(g.gate(&seg(0.0, PortStatus::Berth, 0.0)), 1.0, 1e-15));
        assert!(close(
            g.gate(&seg(0.0, PortStatus::AtSea, 57.5)),
            0.51 + 0.2 * 0.49,
            1e-12
        ));
        assert!(close(
            g.propagule(PortStatus::Anchorage, Kilometres(250.0)),
            1.0,
            0.0
        ));
    }
    /// Segment getters return exactly what was passed in.
    #[test]
    fn segment_getters_round_trip() {
        let s = Segment::new(
            Days(2.5),
            Knots(3.0),
            Celsius(18.0),
            Psu(30.0),
            PortStatus::Anchorage,
            Kilometres(12.0),
        )
        .unwrap();
        assert_eq!(
            (
                s.duration(),
                s.speed_through_water(),
                s.sst(),
                s.salinity(),
                s.port(),
                s.distance_to_coast()
            ),
            (
                Days(2.5),
                Knots(3.0),
                Celsius(18.0),
                Psu(30.0),
                PortStatus::Anchorage,
                Kilometres(12.0)
            )
        );
        assert!(matches!(
            Segment::new(
                Days(1.0),
                Knots(0.0),
                Celsius(20.0),
                Psu(35.0),
                PortStatus::AtSea,
                Kilometres(-1.0)
            ),
            Err(FoulingError::InvalidSegment("distance_to_coast"))
        ));
    }
}
