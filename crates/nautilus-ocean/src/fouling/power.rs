//! Power chain: converts roughness and coating state into a shaft-power penalty,
//! always relative to the clean as-applied reference (a freshly docked ship has
//! zero penalty).
//!
//! The penalty describes the modelled effect of fouling on a modelled clean-hull
//! power. It must never be added on top of measured fuel consumption, which
//! already contains the fouling effect.

use nautilus_units::{Knots, Metres, Micrometres};

use super::granville::schoenherr_cf;
use super::{FoulingError, OperatingPoint};

/// Seawater kinematic viscosity, m²/s.
const NU_SEAWATER: f64 = 1.19e-6;

/// Ship properties the power chain needs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShipParticulars {
    length: Metres,
    frictional_share: f64,
}

impl ShipParticulars {
    /// Builds particulars. `length` (waterline length) must be finite and > 0;
    /// `frictional_share` (frictional share of total calm-water resistance) must
    /// be in (0, 1].
    pub fn new(length: Metres, frictional_share: f64) -> Result<Self, FoulingError> {
        if !(length.value().is_finite() && length.value() > 0.0) {
            return Err(FoulingError::InvalidParameter("length"));
        }
        if !(frictional_share > 0.0 && frictional_share <= 1.0) {
            return Err(FoulingError::InvalidParameter("frictional_share"));
        }
        Ok(ShipParticulars {
            length,
            frictional_share,
        })
    }
    /// Waterline length.
    pub fn length(&self) -> Metres {
        self.length
    }
    /// Frictional share of total calm-water resistance, in (0, 1].
    pub fn frictional_share(&self) -> f64 {
        self.frictional_share
    }
}

/// Length-based Reynolds number `V·L/ν` with ν = 1.19e-6 m²/s.
pub fn reynolds(length: Metres, speed: Knots) -> f64 {
    speed.to_metres_per_second() * length.value() / NU_SEAWATER
}

/// Clean as-applied hull roughness `k_s` (the zero-penalty reference).
pub const REFERENCE_KS: Micrometres = Micrometres(30.0);

/// Fractional increase in total resistance from fouling roughness `ks`:
/// `share · max(0, G(ks) − G(30 µm))`, `G` = Granville `ΔC_F/C_F` from the table.
/// `ks` below the 30 µm reference counts as the reference.
pub fn fouling_resistance_fraction(op: &OperatingPoint, ks: Micrometres) -> f64 {
    let length = op.ship.length();
    let re = reynolds(length, op.speed);
    let g = |k: Micrometres| op.table.delta_cf_ratio(re, k.to_metres() / length.value());
    let ks = if ks.value() > REFERENCE_KS.value() {
        ks
    } else {
        REFERENCE_KS
    };
    op.ship.frictional_share() * (g(ks) - g(REFERENCE_KS)).max(0.0)
}

/// Fractional increase in total resistance from coating roughness (AHR) above
/// its out-of-docking value:
/// `share · 0.044·((now/L)^(1/3) − (out/L)^(1/3)) / C_F,Schoenherr(Re)`, floored at 0.
/// This is the difference of two ITTC roughness allowances; the Re and constant
/// terms cancel.
pub fn coating_resistance_fraction(
    op: &OperatingPoint,
    ahr_out: Micrometres,
    ahr_now: Micrometres,
) -> f64 {
    let length = op.ship.length().value();
    let re = reynolds(op.ship.length(), op.speed);
    let d = 0.044 * ((ahr_now.to_metres() / length).cbrt() - (ahr_out.to_metres() / length).cbrt());
    (op.ship.frictional_share() * d / schoenherr_cf(re)).max(0.0)
}

/// Propeller open-water efficiency loss by blade roughness `(k_p µm, Δη₀/η₀)`.
/// ITTC strip-method values for an illustrative 7.9 m 5-blade propeller
/// (spec open item 6). Assumption — calibrate.
const PROPELLER_TABLE: [(f64, f64); 6] = [
    (30.0, 0.0),
    (100.0, 0.024),
    (300.0, 0.052),
    (1000.0, 0.094),
    (3000.0, 0.145),
    (10000.0, 0.222),
];

/// Fractional loss of propeller open-water efficiency at blade roughness `k_p`:
/// linear in log10(k_p) through the strip-method table, 0 at and below 30 µm,
/// 0.222 at and above 10 mm.
pub fn propeller_eta_loss(k_p: Micrometres) -> f64 {
    let k = k_p.value();
    let (first, last) = (
        PROPELLER_TABLE[0],
        PROPELLER_TABLE[PROPELLER_TABLE.len() - 1],
    );
    if k.is_nan() || k <= first.0 {
        return first.1;
    }
    if k >= last.0 {
        return last.1;
    }
    for w in PROPELLER_TABLE.windows(2) {
        let ((k0, e0), (k1, e1)) = (w[0], w[1]);
        if k <= k1 {
            let f = (k.log10() - k0.log10()) / (k1.log10() - k0.log10());
            return e0 + f * (e1 - e0);
        }
    }
    last.1
}

/// Power penalty components, all relative to the clean as-applied reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PowerPenalty {
    /// Fractional total-resistance increase from biofouling.
    pub hull_fouling_fraction: f64,
    /// Fractional total-resistance increase from coating roughness (AHR) growth.
    pub coating_fraction: f64,
    /// Fractional loss of propeller open-water efficiency.
    pub propeller_eta_loss: f64,
}

impl PowerPenalty {
    /// Total hull resistance fraction `r` (fouling + coating).
    pub fn hull_resistance_fraction(&self) -> f64 {
        self.hull_fouling_fraction + self.coating_fraction
    }

    /// Shaft-power multiplier `(1 + r)·(1 + 0.1·r) / (1 − η-loss)`.
    ///
    /// The `(1 + 0.1·r)` propeller-loading term is back-calculated from Schultz
    /// (86% power vs 80% resistance), not a published coefficient.
    pub fn power_factor(&self) -> f64 {
        let r = self.hull_resistance_fraction();
        (1.0 + r) * (1.0 + 0.1 * r) / (1.0 - self.propeller_eta_loss)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fouling::OperatingPoint;
    use crate::fouling::granville::GranvilleTable;
    use nautilus_units::{Knots, Metres, Micrometres};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }
    fn tanker() -> ShipParticulars {
        ShipParticulars::new(Metres(176.0), 0.5).unwrap()
    }
    #[test]
    fn power_factor_arithmetic() {
        let p = PowerPenalty {
            hull_fouling_fraction: 0.2,
            coating_fraction: 0.0,
            propeller_eta_loss: 0.1,
        };
        assert!(close(p.power_factor(), 1.2 * 1.02 / 0.9, 1e-12));
        let zero = PowerPenalty {
            hull_fouling_fraction: 0.0,
            coating_fraction: 0.0,
            propeller_eta_loss: 0.0,
        };
        assert_eq!(zero.power_factor(), 1.0);
    }
    #[test]
    fn fouling_fraction_relative_to_reference() {
        let (ship, table) = (tanker(), GranvilleTable::new());
        let op = OperatingPoint {
            ship: &ship,
            speed: Knots(14.0),
            table: &table,
        };
        assert!(close(
            fouling_resistance_fraction(&op, Micrometres(30.0)),
            0.0,
            1e-12
        ));
        assert_eq!(fouling_resistance_fraction(&op, Micrometres(5.0)), 0.0);
        // (21.0% − 1.0%) × 0.5 share at 100 µm; (162.7 − 1.0) × 0.5 at 10 mm
        assert!(close(
            fouling_resistance_fraction(&op, Micrometres(100.0)),
            0.100,
            0.004
        ));
        assert!(close(
            fouling_resistance_fraction(&op, Micrometres(10000.0)),
            0.8085,
            0.006
        ));
    }
    #[test]
    fn coating_fraction_120_to_150() {
        let (ship, table) = (
            ShipParticulars::new(Metres(176.0), 1.0).unwrap(),
            GranvilleTable::new(),
        );
        let op = OperatingPoint {
            ship: &ship,
            speed: Knots(14.0),
            table: &table,
        };
        assert!(close(
            coating_resistance_fraction(&op, Micrometres(120.0), Micrometres(150.0)),
            0.0197,
            0.0005
        ));
        assert_eq!(
            coating_resistance_fraction(&op, Micrometres(120.0), Micrometres(100.0)),
            0.0
        );
    }
    #[test]
    fn propeller_table() {
        assert_eq!(propeller_eta_loss(Micrometres(30.0)), 0.0);
        assert!(close(propeller_eta_loss(Micrometres(1000.0)), 0.094, 1e-12));
        assert_eq!(propeller_eta_loss(Micrometres(50000.0)), 0.222);
        let mid = propeller_eta_loss(Micrometres(173.205)); // log-midpoint of 100 and 300
        assert!(close(mid, 0.038, 1e-3));
    }
    #[test]
    fn ship_validation() {
        assert!(ShipParticulars::new(Metres(0.0), 0.5).is_err());
        assert!(ShipParticulars::new(Metres(100.0), 1.5).is_err());
    }
}
