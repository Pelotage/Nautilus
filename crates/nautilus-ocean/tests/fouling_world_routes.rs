//! The guild model on synthetic world routes: the behaviour a ship crossing climates should
//! show (spec Section 8). Each route is split into daily segments, as an AIS track would be.

mod common;

use common::*;
use nautilus_ocean::fouling::FoulingModel;
use nautilus_ocean::fouling::environment::{PortStatus, Segment};
use nautilus_ocean::fouling::events::MaintenanceEvent;
use nautilus_ocean::fouling::guild_model::{Guild, GuildModel};
use nautilus_ocean::fouling::params::{GuildParams, SOFT};
use nautilus_units::MilligramsPerCubicMetre;

fn sail(m: &mut GuildModel, segs: &[Segment]) {
    segs.iter().for_each(|s| m.advance(s));
}
fn live(m: &GuildModel, g: Guild) -> f64 {
    let s = m.guild(g);
    s.slime + s.soft + s.live_coverage
}

/// Singapore → Rotterdam in winter: a hull fouled for 300 days at 29 °C then berthed 45 days
/// at 8 °C loses over 90 % of its warm community (8 °C is 2 °C below its 10 °C lethal limit),
/// but the dead shells keep the hull at least 90 % as rough, and over the next 60 days the
/// cold-water guild starts settling.
#[test]
fn singapore_to_rotterdam_winter_kills_warm_guild() {
    let mut m = guild();
    sail(&mut m, &berth_leg(300.0, 29.0, 35.0));
    let (warm, ks) = (live(&m, Guild::Warm), m.hull_ks().value());
    sail(&mut m, &berth_leg(45.0, 8.0, 35.0));
    assert!(
        live(&m, Guild::Warm) < 0.1 * warm,
        "warm {} of {warm}",
        live(&m, Guild::Warm)
    );
    assert!(
        m.hull_ks().value() >= 0.9 * ks,
        "k_s {} of {ks}",
        m.hull_ks().value()
    );
    sail(&mut m, &berth_leg(60.0, 8.0, 35.0));
    assert!(m.guild(Guild::Cold).slime > 0.0);
}

/// Rotterdam → Persian Gulf: a cold-water community grown over 200 days at 9 °C dies off in
/// 60 days at 33 °C (9 °C above its 24 °C limit), while the warm guild settles.
#[test]
fn rotterdam_to_gulf_kills_cold_guild() {
    let mut m = guild();
    sail(&mut m, &berth_leg(200.0, 9.0, 35.0));
    let cold = live(&m, Guild::Cold);
    assert!(cold > 0.0);
    sail(&mut m, &berth_leg(60.0, 33.0, 35.0));
    assert!(live(&m, Guild::Cold) < 0.1 * cold);
    assert!(m.guild(Guild::Warm).slime > 0.0);
}

/// A 60-day Baltic stay at 6 PSU: marine soft fouling (stenohaline, nothing below 10 PSU) dies
/// off, while barnacles (euryhaline) survive and keep at least half their coverage.
#[test]
fn baltic_stay_removes_soft_keeps_barnacles() {
    let mut m = guild();
    sail(&mut m, &berth_leg(300.0, 25.0, 35.0));
    let (soft, cover) = (m.soft_rating(), m.live_coverage());
    sail(&mut m, &berth_leg(60.0, 15.0, 6.0));
    assert!(
        m.soft_rating() < 0.01 * soft,
        "soft {} of {soft}",
        m.soft_rating()
    );
    assert!(m.live_coverage() > 0.5 * cover);
}

/// Panama transit, 12 h at 0.3 PSU: soft fouling halves every 3 h (1/16 left) and also suffers
/// brackish death (0.1/day below 10 PSU), so (1/16)·e^(−0.05) remains. Both are exact decays.
#[test]
fn panama_transit_freshwater_kill() {
    let mut m = guild();
    sail(&mut m, &berth_leg(300.0, 25.0, 35.0));
    let soft = m.soft_rating();
    sail(&mut m, &berth_leg(0.5, 25.0, 0.3));
    let want = soft / 16.0 * (-0.1f64 * 0.5).exp();
    assert!(
        ((m.soft_rating() - want) / want).abs() < 1e-9,
        "{} vs {want}",
        m.soft_rating()
    );
}

/// Seasons without a calendar: on a clean hull, a cold month (9 °C) lets only the cold guild
/// settle, and a warm month (20 °C) lets only the warm guild settle. That is how a Cape Town –
/// Durban route moves through the year, with no calendar input.
#[test]
fn seasons_switch_settlement_between_guilds_without_calendar() {
    let mut winter = guild();
    sail(&mut winter, &berth_leg(30.0, 9.0, 35.0));
    assert!(winter.guild(Guild::Cold).slime > 0.0 && winter.guild(Guild::Warm).slime == 0.0);
    let mut summer = guild();
    sail(&mut summer, &berth_leg(30.0, 20.0, 35.0));
    assert!(summer.guild(Guild::Warm).slime > 0.0 && summer.guild(Guild::Cold).slime == 0.0);
}

/// No dead band: at every temperature between the two guilds' breeding ranges (12.5–15.5 °C;
/// San Francisco, Cape Town, the Channel), a hull idle for a year fouls. (Review finding C1:
/// with non-overlapping ranges a ship laid up at 13.5 °C stayed perfectly clean forever.)
#[test]
fn no_dead_band_between_guild_breeding_ranges() {
    for t in [12.5, 13.0, 13.5, 14.0, 14.5, 15.5] {
        let mut m = guild();
        sail(&mut m, &berth_leg(365.0, t, 35.0));
        assert!(m.slime_rating() > 1.0, "{t} °C: slime {}", m.slime_rating());
    }
}

fn year_idle_hull_fraction(t: f64) -> f64 {
    let mut m = guild();
    sail(&mut m, &berth_leg(365.0, t, 35.0));
    with_tanker(|op| m.power_penalty(op).hull_fouling_fraction)
}

/// Cold water fouls less than temperate water: a year idle at −2 °C or 5 °C gives a smaller
/// hull penalty than at 17 °C. (Review finding I1: soft fouling saturated in near-freezing water.)
#[test]
fn colder_water_fouls_less() {
    let (polar, cold, med) = (
        year_idle_hull_fraction(-2.0),
        year_idle_hull_fraction(5.0),
        year_idle_hull_fraction(17.0),
    );
    assert!(
        polar <= cold && cold < med,
        "−2 °C {polar}, 5 °C {cold}, 17 °C {med}"
    );
}

/// Hot water is not a refuge: over 5 years on the GloFouling profile, a Persian Gulf summer
/// temperature (33 °C) fouls at least as much as the 25 °C equatorial site. (Review finding I2:
/// fouling collapsed above 32 °C.)
#[test]
fn gulf_heat_fouls_at_least_as_much_as_equator() {
    use nautilus_ocean::fouling::calibration::glofouling_profile;
    use nautilus_ocean::fouling::simulate::{CleaningPolicy, simulate};
    let mean = |sst: f64| {
        let mut m = guild();
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
    };
    let (gulf, eq) = (mean(33.0), mean(25.0));
    assert!(gulf >= eq, "33 °C {gulf} < 25 °C {eq}");
}

/// Light grooming never beats a full in-water clean on the same heavily fouled hull (it keeps
/// dead shells and treats the same reachable area). (Review finding I4.)
#[test]
fn grooming_never_beats_full_clean() {
    let mut m = guild();
    sail(&mut m, &berth_leg(900.0, 28.0, 35.0));
    let after = |e: MaintenanceEvent| {
        let mut c = m.clone();
        c.apply(e);
        with_tanker(|op| c.power_penalty(op).power_factor())
    };
    let groom = after(MaintenanceEvent::HullGrooming);
    let clean = after(MaintenanceEvent::InWaterHullClean {
        calcareous_remnant: 0.2,
    });
    assert!(groom >= clean, "groom {groom} < clean {clean}");
}

fn mixed_route(chl: Option<f64>) -> Vec<Segment> {
    let mut route = berth_leg(120.0, 28.0, 35.0);
    route.extend(leg(30.0, 26.0, 35.0, 14.0, PortStatus::AtSea, 300.0));
    route.extend(berth_leg(200.0, 18.0, 34.0));
    route.extend(leg(20.0, 22.0, 30.0, 3.0, PortStatus::Anchorage, 10.0));
    route.extend(berth_leg(360.0, 27.0, 35.0));
    match chl {
        None => route,
        Some(c) => route
            .into_iter()
            .map(|s| s.with_chlorophyll(MilligramsPerCubicMetre(c)).unwrap())
            .collect(),
    }
}

/// Food supply: with α = 0, adding chlorophyll changes nothing (bit-identical state); with
/// α = 0.5 the hull never ends up with less fouling in total than without food.
#[test]
fn chlorophyll_alpha_zero_is_identical_and_positive_never_less() {
    let run = |p: GuildParams, chl| {
        let mut m = GuildModel::new(p).unwrap();
        sail(&mut m, &mixed_route(chl));
        m
    };
    let d = GuildParams::default();
    assert_eq!(run(d, None), run(d, Some(3.0)));
    let fed = GuildParams {
        food_alpha: 0.5,
        ..d
    };
    let (base, more) = (run(fed, None), run(fed, Some(3.0)));
    assert!(more.slime_rating() >= base.slime_rating() && more.soft_rating() >= base.soft_rating());
    assert!(
        more.calcareous_coverage() >= base.calcareous_coverage()
            && more.hull_ks() >= base.hull_ks()
    );
}

/// Fast liner: with detachment enabled (0.05/day per knot above 10 kn) soft fouling falls during
/// 20 days at 18 kn offshore; with detachment off (self-polishing default) it is unchanged.
#[test]
fn fast_liner_detachment() {
    let run = |rate: f64| {
        let mut p = GuildParams::default();
        p.detachment[SOFT].rate_per_knot = rate;
        let mut m = GuildModel::new(p).unwrap();
        sail(&mut m, &berth_leg(200.0, 25.0, 35.0));
        let before = m.soft_rating();
        sail(
            &mut m,
            &leg(20.0, 25.0, 35.0, 18.0, PortStatus::AtSea, 300.0),
        );
        (before, m.soft_rating())
    };
    let (before, after) = run(0.05);
    assert!(after < 0.1 * before, "{after} of {before}");
    let (before, after) = run(0.0);
    assert!(((after - before) / before).abs() < 1e-6);
}

/// Moving a heavily fouled hull to 5 °C drops the barnacle capacity to its 10 % floor, but that
/// alone removes nothing: with mortality and weathering switched off, coverage is unchanged.
/// (Review focus 2.)
#[test]
fn capacity_drop_does_not_remove_coverage() {
    let p = GuildParams {
        mortality_max: 0.0,
        weathering_per_day: 0.0,
        ..GuildParams::default()
    };
    let mut m = GuildModel::new(p).unwrap();
    sail(&mut m, &berth_leg(500.0, 28.0, 35.0));
    let cover = m.calcareous_coverage();
    assert!(cover > 10.0);
    sail(&mut m, &berth_leg(30.0, 5.0, 35.0));
    assert!((m.calcareous_coverage() - cover).abs() < 1e-12);
}

const GOLDEN: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/golden/fouling_guild_world_route.csv"
);

/// A fixed 3-year world route (Singapore → Red Sea → Mediterranean → Rotterdam → Baltic →
/// Gulf, dry-dock after 2 years) gives exactly the frozen trajectory in the golden file (to
/// the 6 printed decimals). Regression guard only; correctness rests on the other tests.
/// Regenerate after an intended model change with `UPDATE_GOLDEN=1`.
#[test]
fn guild_world_route_matches_golden() {
    let legs: Vec<Vec<Segment>> = vec![
        berth_leg(60.0, 29.0, 33.0),
        leg(12.0, 30.0, 40.0, 14.0, PortStatus::AtSea, 50.0),
        berth_leg(30.0, 19.0, 38.0),
        leg(10.0, 16.0, 36.0, 14.0, PortStatus::AtSea, 200.0),
        berth_leg(120.0, 8.0, 32.0),
        berth_leg(60.0, 10.0, 6.0),
        leg(25.0, 18.0, 35.0, 13.0, PortStatus::AtSea, 300.0),
        berth_leg(150.0, 33.0, 39.0),
        leg(25.0, 28.0, 35.0, 13.0, PortStatus::AtSea, 300.0),
        berth_leg(238.0, 29.0, 33.0),
        berth_leg(365.0, 27.0, 35.0),
    ];
    let mut out = String::from("day,hull_ks_um,prop_ks_um,power_factor\n");
    let mut m = guild();
    let mut day = 0.0;
    with_tanker(|op| {
        for (i, segs) in legs.iter().enumerate() {
            if i == 10 {
                m.apply(MaintenanceEvent::DryDockFullBlast);
            }
            for s in segs {
                m.advance(s);
                day += s.duration().value();
                out += &format!(
                    "{day:.4},{:.6},{:.6},{:.6}\n",
                    m.hull_ks().value(),
                    m.prop_ks().value(),
                    m.power_penalty(op).power_factor()
                );
            }
        }
    });
    if std::env::var_os("UPDATE_GOLDEN").is_some() {
        std::fs::write(GOLDEN, &out).unwrap();
        return;
    }
    let golden =
        std::fs::read_to_string(GOLDEN).expect("golden file missing: run with UPDATE_GOLDEN=1");
    let (a, b): (Vec<&str>, Vec<&str>) = (out.lines().collect(), golden.lines().collect());
    assert_eq!(a.len(), b.len(), "line count changed");
    for (i, (x, y)) in a.iter().zip(&b).enumerate().skip(1) {
        for (u, v) in x.split(',').zip(y.split(',')) {
            let (p, q): (f64, f64) = (u.parse().unwrap(), v.parse().unwrap());
            assert!(
                (p - q).abs() <= 1.5e-6 * q.abs().max(1.0),
                "line {}: {x}\n   golden: {y}",
                i + 1
            );
        }
    }
}
