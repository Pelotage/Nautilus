//! The example program and a frozen ("golden") copy of its full output.
//!
//! `tests/golden/fouling_scenario.csv` is the example's output when this test was last
//! accepted. Any change to any number in any of the six 5-year trajectories fails the
//! test, so numerical drift from refactoring cannot slip in unnoticed. The golden file is
//! only a regression guard: its correctness rests on the other test suites.
//!
//! After an *intended* model change, regenerate it with
//! `UPDATE_GOLDEN=1 cargo test -p nautilus-ocean --test fouling_example`
//! and review the diff before committing.

#[path = "../examples/fouling_scenario.rs"]
#[allow(dead_code)]
mod scenario;

const GOLDEN: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/golden/fouling_scenario.csv"
);

fn output() -> String {
    let mut buf = Vec::new();
    scenario::write_scenario(&mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

/// The example writes the documented CSV header, rows for all six (stage, policy) runs
/// with finite values and power factor ≥ 1, and six summary lines in which proactive <
/// reactive < none for each stage.
#[test]
fn example_output_is_well_formed() {
    let out = output();
    let mut lines = out.lines();
    assert_eq!(
        lines.next().unwrap(),
        "day,stage,policy,hull_fouling_fraction,coating_fraction,propeller_eta_loss,power_factor,event"
    );
    let (mut rows, mut summaries) = (0, std::collections::BTreeMap::new());
    for line in lines {
        if let Some(rest) = line.strip_prefix("# summary ") {
            let parts: Vec<&str> = rest.split(' ').collect();
            let value: f64 = parts[2]
                .strip_prefix("mean_power_factor=")
                .unwrap()
                .parse()
                .unwrap();
            summaries.insert((parts[0].to_string(), parts[1].to_string()), value);
            continue;
        }
        let f: Vec<&str> = line.split(',').collect();
        assert_eq!(f.len(), 8, "{line}");
        assert!(
            ["stage1", "stage2"].contains(&f[1])
                && ["none", "reactive", "proactive"].contains(&f[2]),
            "{line}"
        );
        let nums: Vec<f64> = f[3..7].iter().map(|x| x.parse().unwrap()).collect();
        assert!(
            nums.iter().all(|x| x.is_finite() && *x >= 0.0) && nums[3] >= 1.0,
            "{line}"
        );
        rows += 1;
    }
    assert!(rows > 6 * 366, "only {rows} rows");
    assert_eq!(summaries.len(), 6);
    for stage in ["stage1", "stage2"] {
        let get = |p: &str| summaries[&(stage.to_string(), p.to_string())];
        assert!(
            get("proactive") < get("reactive") && get("reactive") < get("none"),
            "{stage}: {summaries:?}"
        );
    }
}

/// Every number in the example output matches the golden file (to the 6 printed
/// decimals, ±1 in the last place for cross-platform rounding); every text field matches
/// exactly.
#[test]
fn example_output_matches_golden() {
    let out = output();
    if std::env::var_os("UPDATE_GOLDEN").is_some() {
        std::fs::write(GOLDEN, &out).unwrap();
        return;
    }
    let golden =
        std::fs::read_to_string(GOLDEN).expect("golden file missing: run with UPDATE_GOLDEN=1");
    let (a, b): (Vec<&str>, Vec<&str>) = (out.lines().collect(), golden.lines().collect());
    assert_eq!(a.len(), b.len(), "line count changed");
    for (i, (x, y)) in a.iter().zip(&b).enumerate() {
        let (xs, ys): (Vec<&str>, Vec<&str>) =
            (x.split([',', '=']).collect(), y.split([',', '=']).collect());
        assert_eq!(xs.len(), ys.len(), "line {}: {x} vs {y}", i + 1);
        for (u, v) in xs.iter().zip(&ys) {
            match (u.parse::<f64>(), v.parse::<f64>()) {
                (Ok(p), Ok(q)) => assert!(
                    (p - q).abs() <= 1.5e-6,
                    "line {}: {x}\n   golden: {y}",
                    i + 1
                ),
                _ => assert_eq!(u, v, "line {}", i + 1),
            }
        }
    }
}
