# Global Guild Fouling Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `GuildModel` to `nautilus-ocean::fouling`. It is a state-based fouling model with
warm and cold thermal guilds, mortality, breeding seasons, food and detachment, calibrated in
three layers, so that hull and propeller fouling is sensible on any world route. The existing
"stage" models are renamed by how they work.

**Architecture:** Fouling amounts (ratings, live coverage, dead shells) are the state. Each
segment's environment is turned into a `SegmentRates` struct of constant rates. A Strang-split
solver advances the state, with exact exponential decay for linear losses and RK4 for gains,
in fixed steps of at most 0.25 day. The model implements the existing `FoulingModel` trait, so
`simulate`, the cleaning policies and the power chain are reused unchanged. Calibration scoring
lives in the library (`calibration.rs`); a deterministic Nelder–Mead example fits the
constants, which are then hard-coded as defaults.

**Tech Stack:** Rust 1.99 (edition 2024), std only.

**Spec:** `docs/superpowers/specs/2026-11-07-fouling-guild-model-design.md` (read it; section
numbers below refer to it). Background: `research/reports/Biofouling growth dynamical model.md`.

## Global Constraints

- `f64` only; std only; no new dependencies; no `unsafe`.
- Every `pub` item has a doc comment. `cargo clippy -p nautilus-ocean --all-targets -- -D warnings`
  must be clean.
- Physical quantities use `nautilus-units` newtypes. Ratios, rates (per day), ratings and
  coverage % are plain `f64`.
- Deterministic: fixed iteration and step counts; no tolerance loops; no `HashMap` in outputs.
- Never use `f64::min` / `f64::max` / `clamp` on *state* values where a NaN could be hidden.
  Use `if x > cap { cap } else { x }` (NaN passes through).
- Every test has a `///` comment saying what it checks and where its expected values come from.
- Every default that is a prior and not fitted carries "Assumption — calibrate." in its doc
  comment. Fitted defaults name `examples/calibrate_fouling.rs`, the date and the objective value.
- Commit messages have **no** AI co-author or attribution trailers.
- Toolchain: run `export PATH="$HOME/.cargo/bin:$PATH"` before cargo. Use `set -o pipefail`
  when piping test output.
- Existing behaviour must not change: after Task 1 all 97 existing tests pass (3 ignored), and
  `tests/golden/fouling_scenario.csv` numbers are identical.

## Review Focus

1. **Rating slightly above capacity from round-off** (F_i = K_i·(1 + 1e-12)): `√ln(K/F)` must not
   become NaN. The growth term is 0 whenever ln(K/F) ≤ 0. Test in Task 6.
2. **Capacity dropping below present coverage** (warm hull moved to 5 °C, so K_h = 10%): free
   space Φ is 0, no gain term goes negative, and coverage does not fall from the capacity drop
   alone. Test in Task 6, and the route test in Task 11.
3. **Temperature exactly on a window edge** (T = T_opt, T_max, L_lo, B_lo): rates are finite and
   continuous where the spec defines a ramp. T_opt = T_max is rejected at construction. Tests in
   Tasks 3 and 4.
4. **Extreme chlorophyll** (1e6 mg/m³ with α = 0.5): the food factor stays below 1 + α and the
   state stays finite. Test in Task 4, and the random history in Task 11.
5. **Very long or very short segments** (a 3,650-day berth stay; 1e-6 d): the step count is
   ⌈Δ/0.25⌉ with at least 1 step. A long stay equals the same time split into daily segments to
   within 1e-6. Test in Task 5, and the random history in Task 11.

---

## File map

| File | Responsibility |
|---|---|
| `crates/nautilus-units/src/lib.rs` | add `MilligramsPerCubicMetre` |
| `fouling/environment.rs` | `Segment::with_chlorophyll`, `Segment::chlorophyll` |
| `fouling/params.rs` | `ExposureParams`, `AgeClockParams`, `CoatingType`, `GuildParams` (+ sub-structs), validation helpers |
| `fouling/exposure_model.rs` | renamed `stage1.rs` (`ExposureModel`) |
| `fouling/age_clock_model.rs` | renamed `stage2.rs` (`AgeClockModel`) |
| `fouling/rates.rs` | pure rate functions; `SegmentRates` |
| `fouling/integrate.rs` | `LinearLoss`, `advance`, `advance_with_step` |
| `fouling/dynamics.rs` | state layout `ix`, `gains`, `losses` |
| `fouling/guild_model.rs` | `Guild`, `GuildState`, `GuildModel` |
| `fouling/calibration.rs` | Uzun Mediterranean reference curves, layer scoring, `nelder_mead` |
| `examples/calibrate_fouling.rs` | runs the three fits and prints constants |
| `tests/fouling_validation.rs` | required calibration-target tests for `GuildModel` |
| `tests/fouling_world_routes.rs` | world-route scenarios and the guild golden file |
| `tests/fouling_properties.rs` | long random history extended to `GuildModel` |
| `docs/methodology/ocean-fouling.md` | equations, parameters, calibration results, names |

All `fouling/` paths are under `crates/nautilus-ocean/src/`.

---

### Task 1: Rename the "stage" models and centralise parameters

**Files:**
- Rename: `fouling/stage1.rs` → `fouling/exposure_model.rs`; `fouling/stage2.rs` → `fouling/age_clock_model.rs` (`git mv`)
- Create: `fouling/params.rs`
- Modify: `fouling/mod.rs`, `fouling/simulate.rs` (tests), every file in `crates/nautilus-ocean/tests/` and `examples/`, `tests/golden/fouling_scenario.csv`, `docs/methodology/ocean-fouling.md`

**Interfaces:**
- Produces:
  - `fouling::exposure_model::ExposureModel` (was `Stage1Model`);
  - `fouling::age_clock_model::AgeClockModel` (was `Stage2Model`);
  - `fouling::params::{ExposureParams, AgeClockParams, CoatingType}`;
  - crate-internal helpers in `params.rs`: `check_non_negative`, `check_positive`,
    `check_fraction`, `check_responses` (moved from `stage1.rs`, same signatures).
  - The models' `new(params)` signatures are otherwise unchanged.

- [ ] **Step 1: Capture the baseline.**
  `cp crates/nautilus-ocean/tests/golden/fouling_scenario.csv $CLAUDE_JOB_DIR/tmp/golden_before.csv`
  and record the per-binary pass counts from `cargo test --workspace`.
- [ ] **Step 2: Move and rename.**
  - `git mv` the two files.
  - Move the parameter structs, their `Default` impls, `CoatingType` and the check helpers into
    `params.rs`.
  - Rename every identifier across the crate, tests and example: `Stage1Model`→`ExposureModel`,
    `Stage1Params`→`ExposureParams`, `Stage2Model`→`AgeClockModel`,
    `Stage2Params`→`AgeClockParams`, `stage1::`→`exposure_model::`,
    `stage2::`→`age_clock_model::`; test function names containing `stage1`/`stage2` →
    `exposure`/`age_clock`.
  - The example's `stage` column values become `exposure` and `age_clock`. Update module docs and
    the methodology doc ("Stage 1" → "exposure model", "Stage 2" → "age-clock model").
- [ ] **Step 3: Regenerate the golden file and prove the numbers did not move.**
  `UPDATE_GOLDEN=1 cargo test -p nautilus-ocean --test fouling_example`, then
  `diff <(cut -d, -f1,4- $CLAUDE_JOB_DIR/tmp/golden_before.csv) <(cut -d, -f1,4- crates/nautilus-ocean/tests/golden/fouling_scenario.csv)`.
  Expected: no output (only the `stage` column and the summary-line labels change).
- [ ] **Step 4: Verify.**
  `cargo test --workspace && cargo clippy --workspace --all-targets -- -D warnings`.
  Expected: the same pass counts as Step 1 (3 ignored), clippy clean.
  `grep -rn "Stage1\|Stage2\|stage1\|stage2" crates docs/methodology` → no matches.
- [ ] **Step 5: Commit** `refactor(ocean/fouling): name models by mechanism (exposure, age-clock); centralise params`

---

### Task 2: Chlorophyll unit and segment field

**Files:** Modify `crates/nautilus-units/src/lib.rs`, `fouling/environment.rs`

**Interfaces:**
- Produces:
  - `nautilus_units::MilligramsPerCubicMetre(pub f64)` with `value()`; derives as the other units.
  - `Segment::with_chlorophyll(self, chl: MilligramsPerCubicMetre) -> Result<Segment, FoulingError>`.
  - `Segment::chlorophyll(&self) -> Option<MilligramsPerCubicMetre>`; `Segment::new` sets it to `None`.

- [ ] **Step 1: Write the failing tests.**
  - In units, `chlorophyll_unit`: `assert_eq!(MilligramsPerCubicMetre(2.5).value(), 2.5)`.
  - In environment, `chlorophyll_is_optional_and_validated`:
    - a new segment's `chlorophyll()` is `None`;
    - `with_chlorophyll(MilligramsPerCubicMetre(0.8))` gives `Some(MilligramsPerCubicMetre(0.8))`
      and leaves the other getters unchanged;
    - `-0.1`, `NaN` and `INFINITY` give `Err(FoulingError::InvalidSegment("chlorophyll"))`.
- [ ] **Step 2: Run** `cargo test -p nautilus-units -p nautilus-ocean --lib chlorophyll`.
  Expected: compile FAIL (missing items).
- [ ] **Step 3: Implement.** Add a private field `chlorophyll: Option<MilligramsPerCubicMetre>` to
  `Segment`.
- [ ] **Step 4: Run** the same command. Expected: PASS, and the full suite is unchanged.
- [ ] **Step 5: Commit** `feat(ocean/fouling): optional chlorophyll on segments`

---

### Task 3: `GuildParams` and validation

**Files:** Modify `fouling/params.rs` (tests in its `mod tests`)

**Interfaces:**
- Produces, in `fouling::params`. All are `#[derive(Debug, Clone, Copy, PartialEq)]` with pub
  fields:

```rust
pub struct GuildWindow {          // one thermal guild
    pub rate_ref_temp: Celsius,   // temperature at which ClassRates' *_ref apply (warm 25, cold 9)
    pub t_opt: Celsius, pub t_max: Celsius,
    pub lethal_low: Celsius, pub lethal_high: Celsius,
    pub breed_low: Celsius, pub breed_high: Celsius,
}
pub struct ClassRates { pub settle_ref: f64, pub growth_ref: f64, pub settle_q10: f64, pub growth_q10: f64 }
pub struct Detachment { pub rate_per_knot: f64, pub threshold: Knots }
pub struct GuildParams {
    pub coating: CoatingType,
    pub gate: GateParams,
    pub soft_salinity: SalinityResponse, pub hard_salinity: SalinityResponse,
    pub guilds: [GuildWindow; 2],           // index 0 = warm, 1 = cold
    pub classes: [ClassRates; 3],           // index 0 = slime, 1 = soft, 2 = calcareous
    pub mortality_max: f64,                 // per day
    pub mortality_ramp: Celsius,            // width of the lethal ramp (3 °C)
    pub breeding_ramp: Celsius,             // width of the breeding ramp (1 °C)
    pub fresh_threshold: Psu,               // 2 PSU
    pub fresh_mortality_soft: f64,          // per day, ln 2 / 0.125
    pub fresh_mortality_hard: f64,          // per day, ln 2 / 3
    pub food_alpha: f64, pub food_half_saturation: MilligramsPerCubicMetre,
    pub detachment: [Detachment; 2],        // index 0 = slime, 1 = soft
    pub capacity_17: f64, pub capacity_25: f64, pub capacity_min: f64,  // % coverage
    pub weathering_per_day: f64,
    pub high_coverage_slope: f64,           // µm per % above 50 % coverage (Uzun: 80)
    pub clean_damage: Micrometres,          // κ: AHR added per 100 % coverage removed
    pub ahr_out_docking: Micrometres, pub ahr_drift_per_year: Micrometres, pub spot_repair_recovery: f64,
    pub prop_clean: Micrometres, pub prop_polished: Micrometres,
    pub prop_drift_per_year: Micrometres, pub prop_rough_cap: Micrometres,
}
impl GuildParams { pub fn validate(&self) -> Result<(), FoulingError>; }
impl Default for GuildParams
```

Prior defaults (spec Section 5.3):
- **Warm guild:** `{25, 29, 35, 10, 35, 15, 35}`. **Cold guild:** `{9, 13, 24, −2, 24, −2, 12}`.
- `mortality_max 0.1`, ramps `3` and `1`, `fresh_threshold 2`, fresh rates `ln2/0.125` and `ln2/3`.
- Food: `alpha 0`, `half_saturation 1`.
- Detachment: slime `{0, 30 kn}`, soft `{0, 10 kn}`.
- `capacity_min 10`, `spot_repair_recovery 0.5`.
- Coating and propeller values as `AgeClockParams`.
- Provisional values for fitted fields (replaced in Tasks 8–10):
  - slime `{0.01, 2/37.08, 2, 4.16}`;
  - soft `{0.001, 2/73.11, 2, 1.54}`;
  - calcareous `{0.01, 0.0407, 2, 2.04}`;
  - `capacity_17 50`, `capacity_25 100`, `weathering 0.001`, `high_coverage_slope 80`,
    `clean_damage 0`.

Validation rules: every number finite; rates, ramps, capacities and Q10s > 0 where they are
divisors or exponent bases; fractions (`food_alpha` excepted) in [0, 1]; per guild `t_opt < t_max`,
`lethal_low < lethal_high`, `breed_low < breed_high`; capacities in (0, 100];
`capacity_min ≤ min(capacity_17, capacity_25)`; `prop_cap` style limits as `AgeClockParams`;
salinity and gate rules via `check_responses`. A foul-release coating is rejected only by
`GuildModel::new` (`UnsupportedCoating`), not by `validate`.

- [ ] **Step 1: Write the failing tests.**
  - `defaults_are_valid`.
  - `every_guild_rule_rejects_bad_values`: a table of setters, each applying NaN, a negative
    value or an out-of-order window, asserting `Err(InvalidParameter(_))`.
  - `window_edges_rejected` (Review Focus 3): `t_opt == t_max`, `lethal_low == lethal_high` and
    `breed_low == breed_high` are each rejected.
  - `capacity_bounds`: 0, 101, and `capacity_min > capacity_17` are each rejected.
- [ ] **Step 2: Run** `cargo test -p nautilus-ocean --lib params`. Expected: compile FAIL.
- [ ] **Step 3: Implement** the structs, `Default` and `validate`.
- [ ] **Step 4: Run** the same command. Expected: PASS.
- [ ] **Step 5: Commit** `feat(ocean/fouling): GuildParams with priors and validation`

---

### Task 4: Rate functions (`rates.rs`)

**Files:** Create `fouling/rates.rs`; add `pub mod rates;` in `mod.rs`

**Interfaces:**
- Consumes: `GuildParams`, `GuildWindow`, `ClassRates`, `Detachment` (Task 3); `Segment` (Task 2).
- Produces:

```rust
pub fn thermal_rate(ref_rate: f64, q10: f64, w: &GuildWindow, t: Celsius) -> f64
    // ref_rate·q10^((T − rate_ref_temp)/10) for T ≤ t_opt; linear from that value at t_opt to 0 at t_max; 0 above
pub fn breeding(w: &GuildWindow, ramp: Celsius, t: Celsius) -> f64          // 1 on [breed_low, breed_high], linear to 0 over `ramp`
pub fn thermal_mortality(w: &GuildWindow, max: f64, ramp: Celsius, t: Celsius) -> f64
pub fn fresh_mortality(rate: f64, threshold: Psu, s: Psu) -> f64              // rate if s < threshold else 0
pub fn brackish_mortality(max: f64, stenohaline: &SalinityResponse, s: Psu) -> f64  // max if s < zero_below else 0
pub fn food(alpha: f64, half_saturation: MilligramsPerCubicMetre, chl: Option<MilligramsPerCubicMetre>) -> f64
pub fn detachment(d: &Detachment, v: Knots) -> f64                          // rate_per_knot·max(0, v − threshold)
pub fn calcareous_capacity(p: &GuildParams, t: Celsius) -> f64              // line through (17, c17), (25, c25), clamped [capacity_min, 100]
pub(crate) struct SegmentRates {
    pub settle: [[f64; 2]; 3],  // [class][guild] = x·s·b_g·σ·φ      (x = settlement·propagule gate)
    pub grow:   [[f64; 2]; 3],  // [class][guild] = γ·r·σ·φ          (γ = GateParams::gate)
    pub loss:   [[f64; 2]; 3],  // [class][guild] = μ_g + μ_fresh(class) + μ_brack(class; slime/soft only) + δ(class; slime/soft only)
    pub capacity: f64,          // K_h(T), %
    pub weathering: f64,        // ω
}
impl SegmentRates { pub(crate) fn new(seg: &Segment, p: &GuildParams) -> SegmentRates }
```

- [ ] **Step 1: Write the failing tests.** The expected values are worked by hand from the formulas.
  - `thermal_rate_curve`, using the warm window and `ref_rate 0.05`, `q10 2`:
    - `T 15` → `0.025`;
    - `T 29` → `0.05·2^0.4`;
    - `T 32` → half of the `T 29` value;
    - `T 35` and `40` → `0`;
    - the cold window at `T 13` equals the warm window at `T 29` (same peak).
  - `breeding_edges`, warm, ramp 1: `14.5` → `0.5`; `15` → `1`; `13.9` → `0`; `35.5` → `0.5`.
  - `thermal_mortality_ramps`, warm, max 0.1, ramp 3:
    - `20` → `0`; `8.5` → `0.05`; `6.0` → `0.1`; `−2` → `0.1`; `36.5` → `0.05`;
    - finite at exactly `10` and `35` (Review Focus 3).
  - `fresh_and_brackish`:
    - `exp(−fresh_mortality(ln2/0.125, 2, 0.3)·0.125) = 0.5` (1e-12);
    - the same rate over 12 h leaves `1/16`;
    - `S = 2` → `0`;
    - brackish at `6 PSU` → `0.1`, at `12 PSU` → `0`.
  - `food_factor` (Review Focus 4): `alpha 0` → `1`; `None` → `1`; `alpha 0.5, chl = K` → `1.25`;
    `chl 1e6` → in `(1.4999, 1.5)`.
  - `detachment_threshold`: `{0.05, 10 kn}` at `9 kn` → `0`, at `14 kn` → `0.2`.
  - `capacity_line`: `17` → `c17`; `25` → `c25`; `21` → midpoint; `5` → `capacity_min`; `35` → `100`.
  - `segment_rates_assemble`:
    - at a berth (x = 1, γ = 1) at 25 °C and 35 PSU, `settle[1][0]` equals
      `thermal_rate(soft.settle_ref, soft.settle_q10, warm, 25)·breeding(warm, 25)` (= `settle_ref`);
    - the cold-guild `loss` at 25 °C is `0.1`.
- [ ] **Step 2: Run** `cargo test -p nautilus-ocean --lib rates`. Expected: compile FAIL.
- [ ] **Step 3: Implement.** Settlement uses `x = gate.settlement(v)·gate.propagule(port, d)`;
  growth uses `gate.gate(seg)`.
- [ ] **Step 4: Run.** Expected: PASS.
- [ ] **Step 5: Commit** `feat(ocean/fouling): guild rate functions`

---

### Task 5: Fixed-step Strang-split solver (`integrate.rs`)

**Files:** Create `fouling/integrate.rs`; `pub(crate) mod integrate;`

**Interfaces:**
- Produces:

```rust
pub(crate) const MAX_STEP_DAYS: f64 = 0.25;
pub(crate) struct LinearLoss<const N: usize> { pub rate: [f64; N], pub into: [Option<usize>; N] }
pub(crate) fn advance<const N: usize>(x: &mut [f64; N], duration: f64, loss: &LinearLoss<N>, gain: impl Fn(&[f64; N]) -> [f64; N])
pub(crate) fn advance_with_step<const N: usize>(x: &mut [f64; N], duration: f64, max_step: f64, loss: &LinearLoss<N>, gain: impl Fn(&[f64; N]) -> [f64; N])
```

Algorithm (`advance` = `advance_with_step(.., MAX_STEP_DAYS, ..)`):
- If `duration == 0` → return.
- `n = max(1, ceil(duration/max_step))`, `h = duration/n`.
- Each step:
  1. **decay(h/2):** compute every `lost_j = x_j·(1 − exp(−rate_j·h/2))` from the pre-decay
     values; then `x_j −= lost_j`; then `x_k += lost_j` for each `into_j = Some(k)`.
  2. **RK4** on `gain` with step `h`.
  3. **decay(h/2).**
  4. Clamp: `if x_j < 0.0 { x_j = 0.0 }` (NaN passes through).

- [ ] **Step 1: Write the failing tests.**
  - `pure_decay_is_exact`: `N = 1`, rate 0.3, zero gain, 10 days → `x0·exp(−3)` (rel 1e-14).
  - `transfer_conserves_mass`: `N = 2`, `rate [0.2, 0]`, `into [Some(1), None]`, zero gain →
    `x0 + x1` constant (1e-12); `x0 = exp(−0.2·T)` exact.
  - `logistic_matches_closed_form`: gain `0.1·x(1 − x)`, no loss, `x0 = 0.01`, 100 days →
    within 1e-8 of `1/(1 + 99·exp(−10))`.
  - `fourth_order_convergence`: the same logistic plus loss 0.05.
    `e(h) = |x_h − x_ref|` with `x_ref` from `max_step 1e-3`; `e(0.5)/e(0.25)` ∈ [12, 20].
  - `splitting_matches` (Review Focus 5):
    - one 3,650-day segment equals 3,650 one-day segments to within 1e-6 (rel);
    - a 1e-6-day segment takes exactly 1 step.
  - `zero_duration_noop` and `bit_identical_replay`.
- [ ] **Step 2: Run** `cargo test -p nautilus-ocean --lib integrate`. Expected: compile FAIL.
- [ ] **Step 3: Implement** per the algorithm.
- [ ] **Step 4: Run.** Expected: PASS.
- [ ] **Step 5: Commit** `feat(ocean/fouling): fixed-step Strang-split integrator`

---

### Task 6: Guild dynamics (`dynamics.rs`)

**Files:** Create `fouling/dynamics.rs`; `pub(crate) mod dynamics;`

**Interfaces:**
- Consumes: `SegmentRates` (Task 4), `LinearLoss` (Task 5).
- Produces:

```rust
pub(crate) const N: usize = 10;
pub(crate) mod ix { pub const SLIME: [usize; 2] = [0, 1]; pub const SOFT: [usize; 2] = [2, 3];
    pub const CALC: [usize; 2] = [4, 5]; pub const DEAD: usize = 6;
    pub const PROP_CALC: [usize; 2] = [7, 8]; pub const PROP_DEAD: usize = 9; }
pub(crate) const RATING_CAPACITY: [f64; 2] = [20.0, 50.0];   // slime, soft
pub(crate) fn gains(x: &[f64; N], r: &SegmentRates) -> [f64; N]
pub(crate) fn losses(r: &SegmentRates) -> LinearLoss<N>
```

Gains (spec 5.2):
- **Slime and soft:** `settle·(1 − F/K) + grow·F_g·√ln(K/F)`, with the growth term `0` when
  `F == 0` or `ln(K/F) ≤ 0`, and the settlement term `0` when `F ≥ K`.
- **Hull and propeller calcareous:** `Φ = K_h − C_w − C_c − D` (or 0 if negative);
  gain `= (settle + grow·C_g)·Φ/K_h`.

Losses:
- slime and soft at `loss[class][g]`, with `into = None`;
- live calcareous at `loss[2][g]`, `into = Some(DEAD)` (propeller: `PROP_DEAD`);
- `DEAD` and `PROP_DEAD` at `weathering`, `into = None`.

- [ ] **Step 1: Write the failing tests.**
  - `clean_hull_without_settlement_stays_clean`.
  - `growth_and_settlement_stop_at_capacity`: soft `F = 50` → soft gains 0.
  - `round_off_above_capacity_is_finite` (Review Focus 1): `F = 50·(1 + 1e-12)` → all gains finite
    and ≥ 0.
  - `no_gain_when_capacity_below_coverage` (Review Focus 2): `capacity 10`, `C_w 30` → calcareous
    gains exactly 0, none negative.
  - `losses_route_live_shells_to_dead`.
  - `gains_never_negative`: 1,000 random valid states and rates (SplitMix64, seed 61) → every gain ≥ 0.
- [ ] **Step 2: Run** `cargo test -p nautilus-ocean --lib dynamics`. Expected: compile FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run.** Expected: PASS.
- [ ] **Step 5: Commit** `feat(ocean/fouling): guild dynamics right-hand side`

---

### Task 7: `GuildModel`

**Files:** Create `fouling/guild_model.rs`; `pub mod guild_model;`; extend the `mod.rs` docs to
describe all three models.

**Interfaces:**
- Consumes: Tasks 3–6; `FoulingModel`, `MaintenanceEvent`, `OperatingPoint`, `PowerPenalty`,
  `fouling_resistance_fraction`, `coating_resistance_fraction`, `propeller_eta_loss`,
  `uzun::{ks_from_rating, ks_from_coverage}`.
- Produces:

```rust
pub enum Guild { Warm, Cold }                                  // derive Debug, Clone, Copy, PartialEq, Eq
pub struct GuildState { pub slime: f64, pub soft: f64, pub live_coverage: f64 }   // derive Debug, Clone, Copy, PartialEq
pub struct GuildModel { /* params, x: [f64; N], ahr: f64, prop_rough: f64, just_docked: bool */ }  // derive Debug, Clone, PartialEq
impl GuildModel {
    pub fn new(params: GuildParams) -> Result<Self, FoulingError>;   // validate(); FoulRelease → UnsupportedCoating
    pub fn params(&self) -> &GuildParams;
    pub fn slime_rating(&self) -> f64; pub fn soft_rating(&self) -> f64;
    pub fn live_coverage(&self) -> f64; pub fn dead_coverage(&self) -> f64; pub fn calcareous_coverage(&self) -> f64;
    pub fn hull_ks(&self) -> Micrometres; pub fn prop_ks(&self) -> Micrometres; pub fn ahr(&self) -> Micrometres;
    pub fn guild(&self, g: Guild) -> GuildState;
}
impl FoulingModel for GuildModel
```

- **`advance`:** if `dt > 0`, clear `just_docked`; `SegmentRates::new`, then
  `integrate::advance(&mut x, dt, &losses(&r), |s| gains(s, &r))`; paint and propeller drift as in
  `AgeClockModel`.
- **Coverage roughness:** `uzun::ks_from_coverage(sc)` for `sc ≤ 50`; above that,
  `k(50) + high_coverage_slope·(sc − 50)`, where `k(50) = 2.4669·2500 − 24.84·50 + 1065.7`.
- **Events:** spec Section 5.4 table exactly.
  - In-water clean: removed = `(1 − ρ)·(ΣC + D)` and `ahr += clean_damage·removed/100`.
  - Slime is scaled by `min(1, 0.75/ΣF_sl)`, applied to both guilds.
- **Observation:** spec Section 5.6.

- [ ] **Step 1: Write the failing tests.**
  - `foul_release_refused`.
  - `fresh_hull_has_unit_power_factor`: `hull_ks == 30 µm`, `power_factor == 1.0`.
  - `zero_duration_is_noop`.
  - `coverage_roughness_continuous_at_50`: equal to Uzun at 25 and 50; at 60 equal to
    `k(50) + 10·slope`.
  - `dry_dock_resets_and_is_idempotent` (state equality after docking twice).
  - `spot_repair_counts_once`: AHR 160 → 140 → 140.
  - `in_water_clean_effects`: after 400 idle days at 28 °C with `clean_damage = 10 µm`, ρ = 0.2 →
    soft 0; slime ≤ 0.75; live 0; dead = `0.2·(live + dead)_before`;
    AHR `+= 10·0.8·(live + dead)_before/100`.
  - `cleaning_a_clean_hull_is_noop` (state equality, AHR unchanged).
  - `grooming_keeps_dead_shells`.
  - `polish_resets_propeller`.
  - `guild_breakdown_sums`: `guild(Warm) + guild(Cold)` slime/soft/live equal the totals.
  - `events_never_increase_fouling_roughness`: random histories, every event → hull and propeller
    k_s ≤ before.
  - A doc example on `GuildModel` (200 idle days at 28 °C, then a dry-dock gives power factor 1).
- [ ] **Step 2: Run** `cargo test -p nautilus-ocean --lib guild_model && cargo test -p nautilus-ocean --doc`.
  Expected: compile FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** both. Expected: PASS. Then run the full `cargo test -p nautilus-ocean` and clippy.
- [ ] **Step 5: Commit** `feat(ocean/fouling): GuildModel with events and observation`

---

### Task 8: Calibration module and Layer 1 (growth)

**Files:** Create `fouling/calibration.rs` (`pub mod calibration;`) and
`examples/calibrate_fouling.rs`; modify `GuildParams::default` and `tests/fouling_validation.rs`.

**Interfaces:**
- Produces, in `fouling::calibration`:

```rust
pub const MED_SLIME: GaussianRating;  pub const MED_SOFT: GaussianRating;      // Uzun Table 2, values in tests/fouling_validation.rs today
pub fn med_coverage(age_days: f64) -> f64;                                       // Uzun Table 3 Mediterranean logistic pair
pub fn eq_coverage(age_days: f64) -> f64;                                        // EQ_CALCAREOUS with EQ_MAIN_ASYM, EQ_FAST_ASYM
#[derive(Debug, Clone, Copy, PartialEq)] pub struct SiteGaps { pub slime: f64, pub soft: f64, pub coverage: f64 }   // worst |model − Uzun|
pub fn layer1_gaps(p: &GuildParams) -> [SiteGaps; 2];   // [25 °C, 17 °C]; 1,500 daily idle berth segments, 35 PSU
pub fn glofouling_profile(years: f64, sst: f64) -> Vec<TimelineItem>;            // same profile as tests/common
pub fn layer2_penalties(p: &GuildParams) -> (f64, f64);  // 5-yr mean hull_fouling_fraction (eq 25 °C, med 17 °C), no cleaning, 176 m / 0.5 / 14 kn
pub fn layer3_ratio(p: &GuildParams) -> f64;             // (pf_before − pf_after DryDockFullBlast)/(pf_before − pf_after InWaterHullClean{0.2}), both + PropellerPolish, after 220 profile items at 25 °C
pub fn nelder_mead(f: impl Fn(&[f64]) -> f64, x0: &[f64], step: &[f64], iterations: usize) -> Vec<f64>;   // standard coefficients 1, 2, 0.5, 0.5; fixed iterations
```

- **Layer 1 fit, per class.**
  - Parameters, in log space: `ln settle_ref`, `ln growth_ref`, `ln settle_q10`, `ln growth_q10`;
    for calcareous also `capacity_17`, `capacity_25` and `ln weathering`.
  - Objective: `max(gap_eq, gap_med) / threshold` for that class.
  - Thresholds: slime 1, soft 2, coverage 3.
  - Nelder–Mead with 800 iterations, starting from the Task 3 provisional values with steps of 0.5
    in log space.

- [ ] **Step 1: Write the failing tests.**
  - In `calibration.rs`:
    - `nelder_mead_finds_rosenbrock_minimum`: (−1.2, 1), 2,000 iterations → within 1e-3 of (1, 1);
    - `nelder_mead_is_deterministic`;
    - `reference_curves`: `med_coverage(695.86)` within 0.01 of 25.0; `eq_coverage(393.1204)` within
      1e-3 of 51.9495.
  - In `tests/fouling_validation.rs`, replace `med_site_shape_calibration_target` (and remove its
    `#[ignore]`) with `guild_layer1_growth_matches_uzun_both_sites`: for both sites,
    `slime ≤ 1.0`, `soft ≤ 2.0`, `coverage ≤ 3.0`.
- [ ] **Step 2: Run** `cargo test -p nautilus-ocean --lib calibration && cargo test -p nautilus-ocean --test fouling_validation guild_layer1`.
  Expected: the unit tests FAIL on compile; once they compile, the Layer 1 test FAILS on the
  provisional defaults.
- [ ] **Step 3: Implement** `calibration.rs` and the `layer1` mode of the example
  (`cargo run --release -p nautilus-ocean --example calibrate_fouling -- layer1` prints the fitted
  `ClassRates`, capacities, weathering and the gaps achieved).
- [ ] **Step 4: Run the fit and copy the printed constants** into `GuildParams::default()`. The doc
  comments must name the tool, the date and the gaps. **If any gap exceeds its threshold:** stop
  and report the best gaps to the user. Do not loosen a threshold.
- [ ] **Step 5: Run** the Step 2 commands again. Expected: PASS. Then run the full suite.
- [ ] **Step 6: Commit** `feat(ocean/fouling): calibration module; Layer 1 growth fit to Uzun 2019`

---

### Task 9: Layer 2 (roughness → power)

**Files:** Modify `examples/calibrate_fouling.rs` (`layer2` mode), `GuildParams::default`
(`high_coverage_slope`), `tests/fouling_validation.rs`.

**Interfaces:** Consumes `layer2_penalties` (Task 8).

- [ ] **Step 1: Write the failing test.** Replace the ignored `glofouling_five_year_penalty_stage2`
  with `guild_layer2_glofouling_five_year_penalties`: `(eq − 0.30).abs() ≤ 0.05`,
  `(med − 0.22).abs() ≤ 0.05`, and `eq > med`.
- [ ] **Step 2: Run** it. Expected: FAIL (slope 80 gives an equatorial value that is too high).
- [ ] **Step 3: Fit.** One-dimensional Nelder–Mead on `high_coverage_slope` (start 80, step 20,
  200 iterations), minimising `(eq − 0.30)²` with Layer 1 fixed. Copy the result into the
  defaults. Same stop rule as Task 8.
- [ ] **Step 4: Run** the test and the full suite. Expected: PASS.
- [ ] **Step 5: Commit** `feat(ocean/fouling): Layer 2 coverage-roughness fit to GloFouling`

---

### Task 10: Layer 3 (cleaning and docking)

**Files:** Modify the example (`layer3` mode), `GuildParams::default` (`clean_damage`), and
`tests/fouling_validation.rs`.

**Interfaces:** Consumes `layer3_ratio` (Task 8).

- [ ] **Step 1: Write the failing tests.**
  - Replace the ignored Adland test with `guild_layer3_dry_dock_recovers_twice_in_water_clean`:
    ratio in [1.6, 2.4].
  - Add `guild_hundley_cleaning_sequence_pattern`: the same assertions as the existing Hundley test,
    run on `GuildModel`.
- [ ] **Step 2: Run** them. Expected: the Adland test FAILS (κ = 0).
- [ ] **Step 3: Fit** `clean_damage` (start 10 µm, step 10, 200 iterations), minimising
  `(ratio − 2)²` with Layers 1–2 fixed. Copy the result into the defaults. Same stop rule.
- [ ] **Step 4: Run.** Expected: PASS. Then run the full suite, including the Task 7 tests.
  `in_water_clean_effects` sets its own `clean_damage`, so it is unaffected.
- [ ] **Step 5: Commit** `feat(ocean/fouling): Layer 3 in-water cleaning damage fit to Adland`

---

### Task 11: World routes, random histories, guild golden file

**Files:** Create `tests/fouling_world_routes.rs` and `tests/golden/fouling_guild_world_route.csv`;
modify `tests/fouling_properties.rs` and `tests/common/mod.rs`.

**Interfaces:**
- Consumes: the public `GuildModel` API.
- Adds `tests/common::guild()` → `GuildModel::new(GuildParams::default()).unwrap()`.
- Adds a voyage helper:
  `fn leg(days: f64, sst: f64, psu: f64, knots: f64, port: PortStatus, km: f64) -> Segment`, which
  splits each leg into daily segments.

Tests (each with its route described in its `///` comment). "Warm live" means the warm guild's
`slime + soft + live_coverage`; the same for "cold live".

| Test | Route | Assertions |
|---|---|---|
| `singapore_to_rotterdam_winter_kills_warm_guild` | 300 d idle at 29 °C, then 45 d at berth at 8 °C | warm live after < 0.1 × before; `hull_ks` after ≥ 0.9 × before; then 60 more days at 8 °C → `guild(Cold)` slime > 0 |
| `rotterdam_to_gulf_kills_cold_guild` | 200 d at 9 °C, then 60 d at 33 °C | cold live < 0.1 × before; `guild(Warm)` slime > 0 |
| `baltic_stay_removes_soft_keeps_barnacles` | 300 d at 25 °C marine, then 60 d at 6 PSU 15 °C | soft < 0.01 × before; live coverage > 0.5 × before |
| `panama_transit_leaves_one_sixteenth_of_soft` | 300 d at 25 °C, then 0.5 d at 0.3 PSU | soft after / before ∈ 1/16 ± 2% (relative) |
| `southern_seasons_switch_settlement_without_calendar` | clean hull, 30 d at 14 °C | stays clean (all guild amounts 0); then 30 d at 18 °C → warm slime > 0 |
| `chlorophyll_alpha_zero_is_identical_and_positive_never_less` | a 2-year mixed route | with α = 0, chlorophyll segments give identical state; with α = 0.5 every state amount ≥ the α = 0 run |
| `fast_liner_detachment` | 200 d idle at 25 °C, then 20 d at 18 kn 300 km offshore | `d_soft = 0.05` → soft falls; `d = 0` → soft unchanged within 1e-6 (rel) |
| `capacity_drop_does_not_remove_coverage` (Review Focus 2) | 500 d at 28 °C, then 30 d at 5 °C | `calcareous_coverage` after ≥ before − loss from mortality: compare with a run where `mortality_max` = 0 and `weathering` = 0, whose coverage must be exactly unchanged |
| `guild_world_route_matches_golden` | a fixed 3-year route (Singapore → Suez → Rotterdam → Baltic → Gulf, dry-dock at 2 years) | the CSV of day, hull_ks, prop_ks, power factor matches the golden file to 1.5e-6 (same `UPDATE_GOLDEN` mechanism as `fouling_example.rs`) |

In `fouling_properties.rs`, extend `long_random_history_keeps_invariants` to `GuildModel`:
- random chlorophyll 0–20 mg/m³ on half the segments;
- all penalties finite and ≥ 0; power factor ≥ 1;
- `slime ≤ 20`, `soft ≤ 50`, `live + dead ≤ 100`;
- no event increases hull or propeller k_s.

- [ ] **Step 1: Write the tests.** Generate the golden file with `UPDATE_GOLDEN=1`.
- [ ] **Step 2: Run** `cargo test -p nautilus-ocean --test fouling_world_routes --test fouling_properties`.
  Expected: PASS. Any failure is a real finding: fix the production code (TDD) and re-run the full
  suite.
- [ ] **Step 3: Commit** `test(ocean/fouling): world-route scenarios, guild random history, guild golden`

---

### Task 12: Documentation and finishing checks

**Files:** Modify `docs/methodology/ocean-fouling.md` and `fouling/mod.rs` docs.

- [ ] **Step 1: Update the methodology note.**
  - Add the guild model section: the equations (spec 5.2), the rate table, the events, the solver,
    and a parameter table with each value, its source and its "fitted" or "assumption — calibrate"
    flag.
  - Add the calibration results per layer, with the gaps achieved.
  - Add a "which model to use" paragraph: exposure for quick estimates, age-clock as the Uzun
    baseline, guild for routes that cross climates.
- [ ] **Step 2: Deliberate-bug sweep.** Plant at least 15 single-line bugs in `rates.rs`,
  `dynamics.rs`, `integrate.rs`, `guild_model.rs` and `calibration.rs`, covering:
  - each rate function;
  - the settlement and growth terms;
  - loss routing to dead;
  - the decay half-step;
  - the RK4 weights;
  - each event row;
  - the coverage-roughness branch.

  Run `cargo test -p nautilus-ocean --no-fail-fast` for each and restore the file afterwards.
  Expected: every planted bug is caught. Add a test for any that survives, confirm it catches the
  bug, and re-run.
- [ ] **Step 3: Final verification.**
  `cargo fmt --all --check && cargo clippy --workspace --all-targets -- -D warnings && cargo test --workspace`.
  Expected: all pass, with 0 ignored tests in `fouling_validation.rs` for the guild model. The
  remaining baseline characterisation tests are not ignored.
- [ ] **Step 4: Commit** `docs(ocean/fouling): guild model methodology and calibration results`
