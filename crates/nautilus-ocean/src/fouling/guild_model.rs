//! Guild model: a state-based fouling model for ships that move around the world
//! (spec `docs/superpowers/specs/2026-11-07-fouling-guild-model-design.md`).
//!
//! The state is the fouling actually present: slime and soft-fouling ratings and live
//! calcareous coverage for a warm-water and a cold-water guild, plus dead shells, on the
//! hull and the propeller. Settlement, growth, thermal and freshwater death, brackish death,
//! detachment and shell weathering all depend on the current environment *and* on the
//! state, so curve shapes, plateaus, die-offs and breeding seasons follow from the route.
//! Integration: fixed-step integrating-factor RK4 (`integrate.rs`).
//!
//! # Example
//!
//! ```
//! use nautilus_ocean::fouling::environment::{PortStatus, Segment};
//! use nautilus_ocean::fouling::events::MaintenanceEvent;
//! use nautilus_ocean::fouling::granville::GranvilleTable;
//! use nautilus_ocean::fouling::guild_model::{Guild, GuildModel};
//! use nautilus_ocean::fouling::params::GuildParams;
//! use nautilus_ocean::fouling::power::ShipParticulars;
//! use nautilus_ocean::fouling::{FoulingModel, OperatingPoint};
//! use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};
//!
//! let mut model = GuildModel::new(GuildParams::default())?;
//! model.advance(&Segment::new(Days(200.0), Knots(0.0), Celsius(28.0), Psu(35.0), PortStatus::Berth, Kilometres(0.0))?);
//! assert!(model.guild(Guild::Warm).slime > 0.0); // tropical berth: the warm guild settles
//!
//! let (ship, table) = (ShipParticulars::new(Metres(176.0), 0.5)?, GranvilleTable::new());
//! let op = OperatingPoint::new(&ship, Knots(14.0), &table)?;
//! assert!(model.power_penalty(&op).power_factor() > 1.0);
//! model.apply(MaintenanceEvent::DryDockFullBlast);
//! assert_eq!(model.power_penalty(&op).power_factor(), 1.0);
//! # Ok::<(), nautilus_ocean::fouling::FoulingError>(())
//! ```

use nautilus_units::{DAYS_PER_YEAR, Micrometres};

use super::dynamics::{N, RATING_CAPACITY, SegmentRates, gains, ix, losses};
use super::environment::Segment;
use super::events::MaintenanceEvent;
use super::integrate;
use super::params::{COLD, CoatingType, GuildParams, WARM};
use super::power::{
    PowerPenalty, REFERENCE_KS, coating_resistance_fraction, fouling_resistance_fraction,
    propeller_eta_loss,
};
use super::uzun::{ks_from_coverage, ks_from_rating};
use super::{FoulingError, FoulingModel, OperatingPoint};

/// Slime rating left by an in-water clean or grooming (≡ k_s 40 µm, GloFouling's floor).
const CLEANED_SLIME_RATING: f64 = 0.75;

/// A thermal guild.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Guild {
    /// Warm-water community (tropical and subtropical).
    Warm,
    /// Cold-water community (boreal and temperate).
    Cold,
}

/// The fouling of one guild on the hull.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GuildState {
    /// Slime rating.
    pub slime: f64,
    /// Soft-fouling rating.
    pub soft: f64,
    /// Live calcareous coverage, %.
    pub live_coverage: f64,
}

/// Guild-model state. Starts freshly docked (no fouling).
#[derive(Debug, Clone, PartialEq)]
pub struct GuildModel {
    params: GuildParams,
    x: [f64; N],
    /// Hull base plates left by the last in-water clean (part of the dead-shell pool).
    plates: f64,
    ahr: f64,
    prop_rough: f64,
    just_docked: bool,
    /// True after an in-water clean that removed fouling, until the next segment with time
    /// in it (a repeated clean misses the same areas, so it is a no-op).
    just_cleaned: bool,
}

impl GuildModel {
    /// Builds a freshly docked model. Rejects invalid parameters with
    /// [`FoulingError::InvalidParameter`] and foul-release paint with
    /// [`FoulingError::UnsupportedCoating`].
    pub fn new(params: GuildParams) -> Result<Self, FoulingError> {
        if params.coating != CoatingType::SelfPolishing {
            return Err(FoulingError::UnsupportedCoating);
        }
        params.validate()?;
        Ok(GuildModel {
            params,
            x: [0.0; N],
            plates: 0.0,
            ahr: params.ahr_out_docking.value(),
            prop_rough: params.prop_clean.value(),
            just_docked: true,
            just_cleaned: false,
        })
    }

    /// The model's parameters.
    pub fn params(&self) -> &GuildParams {
        &self.params
    }
    /// Slime rating, both guilds.
    pub fn slime_rating(&self) -> f64 {
        self.x[ix::SLIME[0]] + self.x[ix::SLIME[1]]
    }
    /// Soft-fouling rating, both guilds.
    pub fn soft_rating(&self) -> f64 {
        self.x[ix::SOFT[0]] + self.x[ix::SOFT[1]]
    }
    /// Live calcareous coverage on the hull, %.
    pub fn live_coverage(&self) -> f64 {
        self.x[ix::CALC[0]] + self.x[ix::CALC[1]]
    }
    /// Dead shell coverage on the hull, %.
    pub fn dead_coverage(&self) -> f64 {
        self.x[ix::DEAD]
    }
    /// All calcareous coverage on the hull (live + dead), %: what roughness sees.
    pub fn calcareous_coverage(&self) -> f64 {
        self.live_coverage() + self.dead_coverage()
    }
    /// Hull roughness: the larger of the rating and coverage roughness.
    pub fn hull_ks(&self) -> Micrometres {
        let rating = ks_from_rating(self.slime_rating() + self.soft_rating());
        let uzun = match self.coverage_ks(self.calcareous_coverage()) {
            Some(k) if k.value() > rating.value() => k,
            _ => rating,
        };
        self.scaled(uzun)
    }
    /// Propeller roughness: the larger of drift roughness and blade-shell roughness.
    pub fn prop_ks(&self) -> Micrometres {
        let shells = self.x[ix::PROP_CALC[0]] + self.x[ix::PROP_CALC[1]] + self.x[ix::PROP_DEAD];
        match self.coverage_ks(shells).map(|k| self.scaled(k)) {
            Some(k) if k.value() > self.prop_rough => k,
            _ => Micrometres(self.prop_rough),
        }
    }
    /// Coating average hull roughness.
    pub fn ahr(&self) -> Micrometres {
        Micrometres(self.ahr)
    }
    /// The fouling of one guild on the hull.
    pub fn guild(&self, g: Guild) -> GuildState {
        let i = match g {
            Guild::Warm => WARM,
            Guild::Cold => COLD,
        };
        GuildState {
            slime: self.x[ix::SLIME[i]],
            soft: self.x[ix::SOFT[i]],
            live_coverage: self.x[ix::CALC[i]],
        }
    }

    /// Applies the Layer 2 roughness scale: 30 µm + λ·(k_s − 30 µm).
    fn scaled(&self, ks: Micrometres) -> Micrometres {
        Micrometres(
            REFERENCE_KS.value()
                + self.params.roughness_scale * (ks.value() - REFERENCE_KS.value()),
        )
    }

    /// Coverage roughness: Uzun Eq. 18 up to 50 % coverage, then continued linearly with the
    /// calibrated slope (`None` at or below 5 %).
    fn coverage_ks(&self, sc: f64) -> Option<Micrometres> {
        if sc > 50.0 {
            let k50 = 2.4669 * 2500.0 - 24.84 * 50.0 + 1065.7;
            Some(Micrometres(
                k50 + self.params.high_coverage_slope * (sc - 50.0),
            ))
        } else {
            ks_from_coverage(sc)
        }
    }

    /// Trims slime to the cleaning floor (both guilds scaled together) and removes soft fouling.
    fn clean_soft_classes(&mut self) {
        let slime = self.slime_rating();
        if slime > CLEANED_SLIME_RATING {
            let f = CLEANED_SLIME_RATING / slime;
            for i in ix::SLIME {
                self.x[i] *= f;
            }
        }
        for i in ix::SOFT {
            self.x[i] = 0.0;
        }
    }

    /// Scales any class total that a solver step pushed above its hard limit back onto it
    /// (slime 20, soft 50, hull and propeller coverage 100 %). This absorbs RK4 overshoot near
    /// saturation. A NaN total fails the `>` test, so it is never hidden.
    fn project_onto_limits(&mut self) {
        let groups: [(&[usize], f64); 4] = [
            (&ix::SLIME, RATING_CAPACITY[0]),
            (&ix::SOFT, RATING_CAPACITY[1]),
            (&[ix::CALC[0], ix::CALC[1], ix::DEAD], 100.0),
            (&[ix::PROP_CALC[0], ix::PROP_CALC[1], ix::PROP_DEAD], 100.0),
        ];
        for (idx, limit) in groups {
            let total: f64 = idx.iter().map(|&i| self.x[i]).sum();
            if total > limit {
                for &i in idx {
                    self.x[i] *= limit / total;
                }
            }
        }
    }

    fn dry_dock(&mut self) {
        self.x = [0.0; N];
        self.just_cleaned = false;
        self.plates = 0.0;
        self.prop_rough = self.params.prop_clean.value();
        self.just_docked = true;
    }
}

impl FoulingModel for GuildModel {
    fn advance(&mut self, seg: &Segment) {
        let dt = seg.duration().value();
        if dt > 0.0 {
            self.just_docked = false;
            self.just_cleaned = false;
        }
        let r = SegmentRates::new(seg, &self.params);
        integrate::advance(&mut self.x, dt, &losses(&r), |s| gains(s, &r));
        self.project_onto_limits();
        let p = &self.params;
        self.ahr += p.ahr_drift_per_year.value() * dt / DAYS_PER_YEAR;
        let rough = self.prop_rough + p.prop_drift_per_year.value() * dt / DAYS_PER_YEAR;
        self.prop_rough = if rough > p.prop_rough_cap.value() {
            p.prop_rough_cap.value()
        } else {
            rough
        };
    }

    fn apply(&mut self, event: MaintenanceEvent) {
        let out = self.params.ahr_out_docking.value();
        match event {
            MaintenanceEvent::DryDockFullBlast => {
                self.dry_dock();
                self.ahr = out;
            }
            MaintenanceEvent::DryDockSpotRepair => {
                // A repeated docking record with no time in between is the same docking.
                if !self.just_docked {
                    self.ahr -= self.params.spot_repair_recovery * (self.ahr - out);
                }
                self.dry_dock();
            }
            MaintenanceEvent::InWaterHullClean { calcareous_remnant } => {
                if self.just_cleaned {
                    return;
                }
                let rho = if calcareous_remnant.is_nan() {
                    0.0
                } else {
                    calcareous_remnant.clamp(0.0, 1.0)
                };
                let eta = self.params.clean_efficiency;
                let before = self.clone();
                // Treated share η: cleaned as on a fully cleaned hull. Untreated share 1 − η:
                // unchanged. The state is the area-weighted mix of the two.
                let slime = self.slime_rating();
                let trim = if slime > CLEANED_SLIME_RATING {
                    CLEANED_SLIME_RATING / slime
                } else {
                    1.0
                };
                for i in ix::SLIME {
                    self.x[i] *= 1.0 - eta + eta * trim;
                }
                for i in ix::SOFT.into_iter().chain(ix::CALC) {
                    self.x[i] *= 1.0 - eta;
                }
                // Base plates from earlier cleans stay; of the shells above them, a fraction ρ
                // stays as new base plates in the treated area.
                let shells = before.calcareous_coverage();
                let dead = before.x[ix::DEAD];
                let old_plates = if self.plates < dead {
                    self.plates
                } else {
                    dead
                };
                let plates = old_plates + rho * (shells - old_plates);
                self.x[ix::DEAD] = eta * plates + (1.0 - eta) * dead;
                self.plates = eta * plates + (1.0 - eta) * old_plates;
                self.ahr += self.params.clean_damage.value() * eta * (shells - plates) / 100.0;
                if self.x != before.x {
                    self.just_cleaned = true;
                }
            }
            MaintenanceEvent::HullGrooming => {
                self.clean_soft_classes();
                for i in ix::CALC {
                    self.x[i] = 0.0;
                }
            }
            MaintenanceEvent::PropellerPolish => {
                let polished = self.params.prop_polished.value();
                if self.prop_rough > polished {
                    self.prop_rough = polished;
                }
                for i in ix::PROP_CALC.into_iter().chain([ix::PROP_DEAD]) {
                    self.x[i] = 0.0;
                }
            }
        }
    }

    fn power_penalty(&self, op: &OperatingPoint) -> PowerPenalty {
        PowerPenalty {
            hull_fouling_fraction: fouling_resistance_fraction(op, self.hull_ks()),
            coating_fraction: coating_resistance_fraction(
                op,
                self.params.ahr_out_docking,
                self.ahr(),
            ),
            propeller_eta_loss: propeller_eta_loss(self.prop_ks()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fouling::environment::{PortStatus, Segment};
    use crate::fouling::granville::GranvilleTable;
    use crate::fouling::power::ShipParticulars;
    use crate::fouling::uzun;
    use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};

    fn m() -> GuildModel {
        GuildModel::new(GuildParams::default()).unwrap()
    }
    fn idle(days: f64, t: f64) -> Segment {
        Segment::new(
            Days(days),
            Knots(0.0),
            Celsius(t),
            Psu(35.0),
            PortStatus::Berth,
            Kilometres(0.0),
        )
        .unwrap()
    }
    fn fouled() -> GuildModel {
        let mut s = m();
        s.advance(&idle(400.0, 28.0));
        s
    }

    /// The guild model is calibrated for self-polishing paint only; foul-release is refused.
    #[test]
    fn foul_release_refused() {
        let p = GuildParams {
            coating: CoatingType::FoulRelease,
            ..GuildParams::default()
        };
        assert!(matches!(
            GuildModel::new(p),
            Err(FoulingError::UnsupportedCoating)
        ));
        let bad = GuildParams {
            mortality_max: f64::NAN,
            ..GuildParams::default()
        };
        assert!(matches!(
            GuildModel::new(bad),
            Err(FoulingError::InvalidParameter(_))
        ));
    }

    /// A fresh hull has no fouling at all, so k_s is exactly the 30 µm reference and the
    /// power factor is exactly 1 (no seed offset, unlike the age-clock model).
    #[test]
    fn fresh_hull_has_unit_power_factor() {
        let (ship, table) = (
            ShipParticulars::new(Metres(176.0), 0.5).unwrap(),
            GranvilleTable::new(),
        );
        let op = OperatingPoint::new(&ship, Knots(14.0), &table).unwrap();
        let s = m();
        assert_eq!(
            (s.hull_ks(), s.prop_ks()),
            (Micrometres(30.0), Micrometres(30.0))
        );
        assert_eq!(s.power_penalty(&op).power_factor(), 1.0);
    }

    /// A zero-length segment changes nothing.
    #[test]
    fn zero_duration_is_noop() {
        let mut s = fouled();
        let before = s.clone();
        s.advance(&idle(0.0, 30.0));
        assert_eq!(s, before);
    }

    /// Coverage roughness equals Uzun Eq. 18 up to 50 % and continues linearly with the
    /// calibrated slope above it, with no jump at 50 %.
    #[test]
    fn coverage_roughness_continuous_at_50() {
        let s = m();
        let k = |sc: f64| s.coverage_ks(sc).unwrap().value();
        assert_eq!(k(25.0), uzun::ks_from_coverage(25.0).unwrap().value());
        assert_eq!(k(50.0), uzun::ks_from_coverage(50.0).unwrap().value());
        assert!((k(60.0) - (k(50.0) + 10.0 * s.params().high_coverage_slope)).abs() < 1e-9);
        assert_eq!(s.coverage_ks(4.0), None);
    }

    /// A full-blast dry-dock returns everything to the as-applied state, and docking twice
    /// equals docking once.
    #[test]
    fn dry_dock_resets_and_is_idempotent() {
        let mut s = fouled();
        s.apply(MaintenanceEvent::DryDockFullBlast);
        let once = s.clone();
        s.apply(MaintenanceEvent::DryDockFullBlast);
        assert_eq!(s, once);
        assert_eq!(
            (s.calcareous_coverage(), s.slime_rating(), s.soft_rating()),
            (0.0, 0.0, 0.0)
        );
        assert_eq!(
            (s.ahr(), s.prop_ks()),
            (Micrometres(120.0), Micrometres(30.0))
        );
    }

    /// A spot repair recovers half the paint roughness gained since docking (160 → 140 µm), and
    /// a duplicated record with no time in between counts once.
    #[test]
    fn spot_repair_counts_once() {
        let mut s = m();
        s.advance(&idle(2.0 * 365.25, 25.0));
        assert!((s.ahr().value() - 160.0).abs() < 1e-9);
        s.apply(MaintenanceEvent::DryDockSpotRepair);
        s.apply(MaintenanceEvent::DryDockSpotRepair);
        assert!((s.ahr().value() - 140.0).abs() < 1e-9);
    }

    /// In-water clean of the whole hull (η = 1, spec 5.4): soft fouling removed, slime trimmed
    /// to the 0.75 floor, live barnacles removed, base plates ρ·(live + dead) left as dead
    /// shell, and paint damage κ × coverage removed / 100.
    #[test]
    fn in_water_clean_effects() {
        let p = GuildParams {
            clean_damage: Micrometres(10.0),
            clean_efficiency: 1.0,
            ..GuildParams::default()
        };
        let mut s = GuildModel::new(p).unwrap();
        s.advance(&idle(400.0, 28.0));
        let (shells, ahr) = (s.calcareous_coverage(), s.ahr().value());
        assert!(shells > 1.0, "need some barnacles to test with: {shells}");
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.2,
        });
        assert_eq!((s.soft_rating(), s.live_coverage()), (0.0, 0.0));
        assert!(s.slime_rating() <= 0.75 + 1e-12);
        assert!((s.dead_coverage() - 0.2 * shells).abs() < 1e-9);
        assert!((s.ahr().value() - (ahr + 10.0 * 0.8 * shells / 100.0)).abs() < 1e-9);
    }

    /// Cleaning is idempotent: a fresh hull is unchanged by a clean, and a second clean right
    /// after the first (no time in between, so the same areas are missed) removes nothing more
    /// and does no more paint damage.
    #[test]
    fn cleaning_a_clean_hull_is_noop() {
        let p = GuildParams {
            clean_damage: Micrometres(10.0),
            ..GuildParams::default()
        };
        let mut fresh = GuildModel::new(p).unwrap();
        let before = fresh.clone();
        fresh.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.2,
        });
        assert_eq!(fresh, before);
        let mut s = GuildModel::new(p).unwrap();
        s.advance(&idle(400.0, 28.0));
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.2,
        });
        let once = s.clone();
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.2,
        });
        assert_eq!(s, once);
    }

    /// Grooming removes slime (to the floor), soft fouling and live barnacles but leaves dead
    /// shells, and does no paint damage.
    #[test]
    fn grooming_keeps_dead_shells() {
        let mut s = fouled();
        s.advance(&idle(60.0, 5.0)); // warm guild dies → dead shells
        let (dead, ahr) = (s.dead_coverage(), s.ahr());
        assert!(dead > 0.0);
        s.apply(MaintenanceEvent::HullGrooming);
        assert_eq!(
            (
                s.live_coverage(),
                s.soft_rating(),
                s.dead_coverage(),
                s.ahr()
            ),
            (0.0, 0.0, dead, ahr)
        );
        assert!(s.slime_rating() <= 0.75 + 1e-12);
    }

    /// Polishing the propeller removes blade barnacles and brings roughness to at most 40 µm.
    #[test]
    fn polish_resets_propeller() {
        let mut s = fouled();
        s.apply(MaintenanceEvent::PropellerPolish);
        assert_eq!(s.prop_ks(), Micrometres(40.0));
        let mut fresh = m();
        fresh.apply(MaintenanceEvent::PropellerPolish);
        assert_eq!(fresh.prop_ks(), Micrometres(30.0)); // never rougher than before
    }

    /// The per-guild breakdown adds up to the totals.
    #[test]
    fn guild_breakdown_sums() {
        let mut s = fouled();
        s.advance(&idle(90.0, 9.0)); // cold guild settles in cold water
        let (w, c) = (s.guild(Guild::Warm), s.guild(Guild::Cold));
        assert!((w.slime + c.slime - s.slime_rating()).abs() < 1e-12);
        assert!((w.soft + c.soft - s.soft_rating()).abs() < 1e-12);
        assert!((w.live_coverage + c.live_coverage - s.live_coverage()).abs() < 1e-12);
    }

    /// No maintenance event makes the hull or propeller rougher, after 200 random histories.
    #[test]
    fn events_never_increase_fouling_roughness() {
        let mut seed = 77u64;
        let mut u = |lo: f64, hi: f64| {
            seed = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = seed;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            lo + (hi - lo) * ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..200 {
            let mut s = m();
            for _ in 0..10 {
                let seg = Segment::new(
                    Days(u(1.0, 60.0)),
                    Knots(u(0.0, 3.0)),
                    Celsius(u(-2.0, 34.0)),
                    Psu(u(0.0, 38.0)),
                    PortStatus::Berth,
                    Kilometres(0.0),
                )
                .unwrap();
                s.advance(&seg);
            }
            for e in [
                MaintenanceEvent::DryDockFullBlast,
                MaintenanceEvent::DryDockSpotRepair,
                MaintenanceEvent::InWaterHullClean {
                    calcareous_remnant: u(0.0, 1.0),
                },
                MaintenanceEvent::HullGrooming,
                MaintenanceEvent::PropellerPolish,
            ] {
                let mut t = s.clone();
                t.apply(e);
                assert!(
                    t.hull_ks() <= s.hull_ks() && t.prop_ks() <= s.prop_ks(),
                    "{e:?}"
                );
            }
        }
    }
    /// Layer 2's roughness scale λ shrinks every Uzun-derived roughness toward the clean 30 µm
    /// reference: k_s = 30 + λ·(k_s,Uzun − 30) on the hull and the blade shells. Propeller drift
    /// roughness is not scaled, and a clean hull stays at exactly 30 µm.
    #[test]
    fn roughness_scale_applies_above_reference() {
        let full = GuildModel::new(GuildParams {
            roughness_scale: 1.0,
            ..GuildParams::default()
        })
        .unwrap();
        let half = GuildModel::new(GuildParams {
            roughness_scale: 0.5,
            ..GuildParams::default()
        })
        .unwrap();
        assert_eq!(half.hull_ks(), Micrometres(30.0));
        let (mut a, mut b) = (full, half);
        for s in [&mut a, &mut b] {
            s.advance(&idle(500.0, 28.0));
        }
        assert!((b.hull_ks().value() - (30.0 + 0.5 * (a.hull_ks().value() - 30.0))).abs() < 1e-9);
        a.apply(MaintenanceEvent::PropellerPolish);
        b.apply(MaintenanceEvent::PropellerPolish);
        assert_eq!(a.prop_ks(), b.prop_ks()); // drift roughness only: unscaled
    }
    /// Partial cleaning (efficiency η = 0.6): the treated 60 % of the hull is cleaned as
    /// usual, the untreated 40 % keeps its fouling, and paint damage scales with what was
    /// actually removed. Expected values: the area-weighted mix worked by hand.
    #[test]
    fn partial_clean_mixes_cleaned_and_untouched_areas() {
        let p = GuildParams {
            clean_damage: Micrometres(10.0),
            clean_efficiency: 0.6,
            ..GuildParams::default()
        };
        let mut s = GuildModel::new(p).unwrap();
        s.advance(&idle(500.0, 28.0));
        let (soft, live, dead, slime, ahr) = (
            s.soft_rating(),
            s.live_coverage(),
            s.dead_coverage(),
            s.slime_rating(),
            s.ahr().value(),
        );
        let shells = live + dead;
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.2,
        });
        assert!((s.soft_rating() - 0.4 * soft).abs() < 1e-9);
        assert!((s.live_coverage() - 0.4 * live).abs() < 1e-9);
        assert!((s.dead_coverage() - (0.6 * 0.2 * shells + 0.4 * dead)).abs() < 1e-9);
        assert!((s.slime_rating() - (0.6 * slime.min(0.75) + 0.4 * slime)).abs() < 1e-9);
        assert!((s.ahr().value() - (ahr + 10.0 * 0.6 * 0.8 * shells / 100.0)).abs() < 1e-9);
        let once = s.clone();
        s.apply(MaintenanceEvent::InWaterHullClean {
            calcareous_remnant: 0.2,
        });
        assert_eq!(s, once);
    }
    /// Ratings and coverage never exceed their hard limits (slime 20, soft 50, coverage 100 %),
    /// even when fast, well-fed growth makes an RK4 step overshoot near capacity.
    #[test]
    fn saturation_never_overshoots_limits() {
        use nautilus_units::MilligramsPerCubicMetre;
        let mut s = GuildModel::new(GuildParams {
            food_alpha: 0.5,
            ..GuildParams::default()
        })
        .unwrap();
        let day = idle(1.0, 29.0)
            .with_chlorophyll(MilligramsPerCubicMetre(3.0))
            .unwrap();
        for _ in 0..1500 {
            s.advance(&day);
            assert!(
                s.slime_rating() <= 20.0 && s.soft_rating() <= 50.0,
                "{} {}",
                s.slime_rating(),
                s.soft_rating()
            );
            assert!(s.calcareous_coverage() <= 100.0);
        }
    }
}
