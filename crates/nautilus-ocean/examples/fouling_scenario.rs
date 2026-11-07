//! Prints a CSV fouling trajectory for the GloFouling equatorial profile
//! (5 years), for Stage 1 and Stage 2 under each cleaning policy, followed by one
//! summary line per (stage, policy) with the mean power factor.
//!
//! Run: `cargo run -p nautilus-ocean --example fouling_scenario`

use nautilus_ocean::fouling::environment::{PortStatus, Segment};
use nautilus_ocean::fouling::granville::GranvilleTable;
use nautilus_ocean::fouling::power::ShipParticulars;
use nautilus_ocean::fouling::simulate::{
    CleaningPolicy, TimelineItem, TrajectoryPoint, mean_power_factor, simulate,
};
use nautilus_ocean::fouling::stage1::{Stage1Model, Stage1Params};
use nautilus_ocean::fouling::stage2::{Stage2Model, Stage2Params};
use nautilus_ocean::fouling::{FoulingError, OperatingPoint};
use std::error::Error;
use std::io::{self, Write};

use nautilus_units::{Celsius, Days, Kilometres, Knots, Metres, Psu};

/// 10-day cycles: 3.0411 d idle at berth + 6.9589 d sailing at 14 kn, 300 km
/// offshore, 25 °C, 35 PSU (556 idle days in 5 years).
fn equatorial_profile(years: f64) -> Result<Vec<TimelineItem>, FoulingError> {
    let cycles = (years * 36.525).ceil() as usize;
    let idle = Segment::new(
        Days(3.0411),
        Knots(0.0),
        Celsius(25.0),
        Psu(35.0),
        PortStatus::Berth,
        Kilometres(0.0),
    )?;
    let sail = Segment::new(
        Days(6.9589),
        Knots(14.0),
        Celsius(25.0),
        Psu(35.0),
        PortStatus::AtSea,
        Kilometres(300.0),
    )?;
    Ok((0..cycles)
        .flat_map(|_| [TimelineItem::Segment(idle), TimelineItem::Segment(sail)])
        .collect())
}

fn main() -> Result<(), Box<dyn Error>> {
    match run() {
        // Output piped into `head` and closed early: not an error.
        Err(e)
            if e.downcast_ref::<io::Error>()
                .is_some_and(|e| e.kind() == io::ErrorKind::BrokenPipe) =>
        {
            Ok(())
        }
        other => other,
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    write_scenario(&mut io::stdout().lock())
}

/// Writes the CSV trajectories and summary lines to `out` (used by `main` and by the
/// `fouling_example` test).
pub fn write_scenario(out: &mut impl Write) -> Result<(), Box<dyn Error>> {
    let timeline = equatorial_profile(5.0)?;
    let (ship, table) = (
        ShipParticulars::new(Metres(176.0), 0.5)?,
        GranvilleTable::new(),
    );
    let op = OperatingPoint::new(&ship, Knots(14.0), &table)?;
    let policies = [
        ("none", CleaningPolicy::NONE),
        ("reactive", CleaningPolicy::REACTIVE),
        ("proactive", CleaningPolicy::PROACTIVE),
    ];

    let mut runs: Vec<(&str, &str, Vec<TrajectoryPoint>)> = Vec::new();
    for (name, policy) in policies {
        runs.push((
            "stage1",
            name,
            simulate(
                &mut Stage1Model::new(Stage1Params::default())?,
                &timeline,
                policy,
                &op,
            ),
        ));
        runs.push((
            "stage2",
            name,
            simulate(
                &mut Stage2Model::new(Stage2Params::default())?,
                &timeline,
                policy,
                &op,
            ),
        ));
    }

    writeln!(
        out,
        "day,stage,policy,hull_fouling_fraction,coating_fraction,propeller_eta_loss,power_factor,event"
    )?;
    for (stage, policy, traj) in &runs {
        for p in traj {
            let event = p.event.map(|e| format!("{e:?}")).unwrap_or_default();
            writeln!(
                out,
                "{:.4},{stage},{policy},{:.6},{:.6},{:.6},{:.6},\"{event}\"",
                p.elapsed.value(),
                p.penalty.hull_fouling_fraction,
                p.penalty.coating_fraction,
                p.penalty.propeller_eta_loss,
                p.penalty.power_factor()
            )?;
        }
    }
    for (stage, policy, traj) in &runs {
        writeln!(
            out,
            "# summary {stage} {policy} mean_power_factor={:.6}",
            mean_power_factor(traj)
        )?;
    }
    Ok(())
}
