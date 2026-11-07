# Global guild fouling model — design

**Date:** 2026-11-07
**Status:** approved in conversation; awaiting written-spec review
**Crate:** `nautilus-ocean`, module `fouling`
**Research spec:** `research/reports/Biofouling growth dynamical model.md`
**Predecessor:** `plans/0002-biofouling-dynamics-core.md` (exposure and age-clock models)

## 1. Purpose

Ships move around the world. Between ports they pass through water from −2 °C to 34 °C, from
fresh to fully marine, rich or poor in food, at speeds from 0 to 25 kn. The fouling model must
give a physically sensible hull and propeller state for **any** such route, and it must agree
with the published data it can be checked against.

The current models (the exposure model and the Uzun age-clock model) cannot do this. Every
growth rate in them depends only on the environment, never on the fouling present. Each state
variable is therefore a time integral of the forcing, and a site's curve is always the
equatorial curve stretched in time. That causes four problems:

1. Curve **shape** and **plateau** cannot change between sites, only timing. This is the
   measured Mediterranean gap: soft fouling up to 6.8 rating points off, coverage up to 7.8
   points off.
2. Nothing dies except in fresh water. A tropical community sailing into the North Sea simply
   grows slowly.
3. The barnacle plateau depends on a running mean temperature, a stand-in that lowers
   coverage with no event when a ship moves to colder water.
4. There are no breeding seasons and no food supply.

The other two calibration gaps follow from the same structure. The equatorial 5-year
penalty is 0.43 against GloFouling's 0.30, and in-water cleaning recovers too much (the
dry-dock to in-water ratio is 1.1–1.3 against Adland's about 2).

**Success means:**
- the three `#[ignore]`d calibration targets in `tests/fouling_validation.rs` pass as required
  tests for the new model;
- synthetic world routes show the expected qualitative behaviour (Section 8);
- all existing tests still pass.

## 2. Scope

**In scope:**
- a new state-based model, `GuildModel`;
- warm-water and cold-water thermal guilds per fouling class;
- thermal and freshwater mortality, with dead shells persisting;
- temperature-driven breeding seasons;
- an optional food (chlorophyll) multiplier;
- speed-driven detachment of soft fouling;
- coating damage from in-water cleaning;
- layered calibration;
- renaming the existing "stage" models.

**Out of scope:**
- sampling real ocean data along AIS tracks (recorded in `docs/future-implementations.md`);
- foul-release coatings (still refused);
- barnacle size, coating chemistry and polishing (the research spec's Stage 3);
- calibration against fleet MRV data.

## 3. Naming and file layout

All files sit directly in `crates/nautilus-ocean/src/fouling/`. There are no subfolders, and
models are named for how they work, not by "stage".

| Now | Becomes | Contents |
|---|---|---|
| `stage1.rs`, `Stage1Model`, `Stage1Params` | `exposure_model.rs`, `ExposureModel`, `ExposureParams` | single exposure clock |
| `stage2.rs`, `Stage2Model`, `Stage2Params` | `age_clock_model.rs`, `AgeClockModel`, `AgeClockParams` | Uzun age clocks, frozen as the baseline |
| — | `params.rs` | `ExposureParams`, `AgeClockParams`, `GuildParams` and the shared validation helpers |
| — | `rates.rs` | pure rate functions of the environment (Section 5.3) |
| — | `dynamics.rs` | the right-hand side dX/dt of the guild model |
| — | `integrate.rs` | the fixed-step solver (Section 5.5) |
| — | `guild_model.rs` | `GuildModel`: state, `FoulingModel` implementation, events, observation |
| `calibration.rs` | the three layers' scoring functions and the Nelder–Mead optimizer, shared by the calibration example and the validation tests |

Unchanged: `environment.rs`, `events.rs`, `granville.rs`, `power.rs`, `simulate.rs`, `uzun.rs`.

The rename also covers:
- the tests;
- the example's `stage` column, which becomes `exposure` / `age_clock` / `guild`;
- the golden file, whose text is regenerated with numbers that must be identical;
- `docs/methodology/ocean-fouling.md`.

It is done first, as a pure rename with no behaviour change.

New in `nautilus-units`: `MilligramsPerCubicMetre(pub f64)`, for chlorophyll-a.

New in `environment.rs`:
- `Segment::with_chlorophyll(self, chl: MilligramsPerCubicMetre) -> Result<Segment, FoulingError>`,
  which rejects negative or non-finite values;
- `Segment::chlorophyll(&self) -> Option<MilligramsPerCubicMetre>`.

Existing `Segment::new` calls are unchanged.

## 4. The dynamical system currently in use (for reference)

**Forcing.** It is piecewise constant per segment k of length Δₖ:
uₖ = (T, S, v, port, d).

**Age-clock model state:**
x = (A_sl, A_so, A_c, A_pc, W, Θ, R, k_p, φ), where:

| Symbol | Meaning |
|---|---|
| A_sl, A_so, A_c | hull age clocks for slime, soft and calcareous fouling |
| A_pc | propeller calcareous clock |
| W, Θ | exposure-weight and weighted-temperature accumulators |
| R | paint roughness (AHR) |
| k_p | propeller roughness |
| φ | freshwater-stay timer |

**Flow:**

```
dA_sl/dt = q_sl(T) σ_soft(S) γ          dA_so/dt = q_so(T) σ_soft(S) γ
dA_c/dt  = dA_pc/dt = dW/dt = q_c(T) σ_hard(S) γ
dΘ/dt    = q_c(T) σ_hard(S) γ · min(T, 32)
dR/dt    = 20 µm/yr                     dk_p/dt = 100 µm/yr (capped at 300 µm)
dφ/dt    = 1 if S < 2 PSU, else φ ← 0
q_i(T)   = Q10_i^((min(T,32) − 25)/10) for T ≥ 12 °C, else 0.1
γ        = x + δ(1 − x),   x = A(v) · P(port, d)
```

**Jumps.**
- A freshwater stay reaching 3 h halves the slime and soft ratings: A ← G⁻¹(½·G(A)).
- Maintenance events follow the reset tables in `docs/methodology/ocean-fouling.md`.

**Observation:**
- FR = G_sl(A_sl) + G_so(A_so), with the G being Uzun's equatorial Gaussians;
- SC = logistic(A_c), with asymptotes from Θ/W;
- k_s = max(k_s(FR), k_s(SC)) from Uzun Eqs. 17–18;
- Granville friction from the docked hull; paint and propeller terms;
- P/P₀ = (1 + r)(1 + 0.1 r)/(1 − Δη).

The exposure model is the one-clock version of the same system: dE/dt = q(T) σ(S) γ, with
hull fraction = min(0.86, 0.0012 E).

## 5. The guild model

### 5.1 State

The hull and the propeller each have their own state. Guild g ∈ {warm, cold}.

| Variable | Meaning | Range |
|---|---|---|
| F_sl,g | slime rating of guild g | Σ_g F_sl,g ≤ K_sl = 20 |
| F_so,g | soft-fouling rating of guild g | Σ_g F_so,g ≤ K_so = 50 |
| C_g | live calcareous coverage of guild g (%) | |
| D | dead shell coverage (%) | C_w + C_c + D ≤ 100 |
| R | paint roughness (AHR, µm) | ≥ out-docking value |
| k_p | propeller blade roughness (µm) | [prop_clean, prop_rough_cap] |
| C_p,g, D_p | calcareous coverage on the propeller | as on the hull |
| `just_docked` | duplicate-docking guard (as now) | bool |

A fresh hull starts at all fouling = 0, R = R_out and k_p = prop_clean.

### 5.2 Equations

Gates: x = A(v)·P(port, d) is settlement and larval supply, and γ = x + δ(1 − x) is growth.
Both are unchanged from today.

**Slime and soft fouling**, for i ∈ {sl, so}, with F_i = Σ_g F_i,g:

```
dF_i,g/dt = σ_i(S)·φ(chl)·[ x·s_i,g(T)·b_g(T)·(1 − F_i/K_i)        settlement into free capacity
                          + γ·r_i,g(T)·F_i,g·√ln(K_i/F_i) ]        growth (Uzun's curve family)
            − [ μ_g(T) + μ_fresh,i(S) + μ_brack,i(S) + δ_i(v) ]·F_i,g   thermal, fresh/brackish, detachment
```

The growth term is defined as 0 at F_i = 0 (its limit). Then F_i,g√ln(K_i/F_i) → 0.

**Calcareous fouling**, with temperature-dependent capacity K_h(T):

```
Φ        = max(0, K_h(T) − C_w − C_c − D)                         free space now
dC_g/dt  = σ_h(S)·φ(chl)·[ x·s_c,g(T)·b_g(T)·Φ/K_h + γ·r_c,g(T)·C_g·Φ/K_h ]
           − [ μ_g(T) + μ_fresh,h(S) ]·C_g
dD/dt    = Σ_g [ μ_g(T) + μ_fresh,h(S) ]·C_g − ω·D                shells persist, weather slowly
```

When a ship moves to colder water and K_h(T) drops below the present coverage, Φ = 0. Growth
stops, but nothing is removed. Coverage never falls without a mortality, weathering or
maintenance cause.

The **propeller** uses the calcareous equations with its own C_p,g and D_p (blades in port
foul like the hull). Its drift roughness k_p grows as today.

**Paint:** dR/dt = ahr_drift (20 µm/yr). **Propeller drift:**
dk_p/dt = 100 µm/yr, capped at 300 µm.

### 5.3 Rate functions (`rates.rs`)

All are pure functions of the environment and parameters. Each guild has its own window
parameters.

| Function | Definition |
|---|---|
| Growth r_i,g(T) | r_ref·Qr^((T − 25)/10) for T ≤ T_opt,g; linear to 0 at T_max,g; 0 above |
| Settlement s_i,g(T) | s_ref·Qs^((T − 25)/10) for T ≤ T_opt,g; linear to 0 at T_max,g; 0 above |
| Breeding b_g(T) | 1 on [B_lo,g, B_hi,g]; linear to 0 over 1 °C outside; 0 beyond |
| Mortality μ_g(T) | 0 on [L_lo,g, L_hi,g]; μ_max·min(1, distance outside / 3 °C) |
| Freshwater μ_fresh,i(S) | μ_f,i for S < 2 PSU, else 0. Soft and slime: μ_f = ln 2 / 0.125 d (half after 3 h, 1/16 left after 12 h, >99.6% killed after 24 h). Live barnacles: ln 2 / 3 d |
| Brackish μ_brack,i(S) | slime and soft only: μ_max for S below the stenohaline zero (10 PSU), else 0. Marine soft fouling cannot survive a long brackish stay |
| Food φ(chl) | 1 + α·chl/(chl + K_chl); 1 when chlorophyll is absent |
| Detachment δ_i(v) | d_i·max(0, v − v_d,i) |
| Salinity σ_i(S) | unchanged: stenohaline for slime and soft, euryhaline for calcareous |
| Capacity K_h(T) | linear through (17 °C, K17) and (25 °C, K25), clamped to [K_min, 100] |

Starting values (priors, all "assumption — calibrate" unless fitted in Section 6):

| Parameter | Warm guild | Cold guild | Basis |
|---|---|---|---|
| T_opt / T_max | 29 / 35 °C | 13 / 24 °C | spec: warm barnacles fail below 12 °C; boreal species |
| Lethal window [L_lo, L_hi] | [10, 35] °C | [−2, 24] °C | Nasrolahi 2016, Anil 1995; MarLIN |
| Breeding [B_lo, B_hi] | [15, 35] °C | [−2, 12] °C | boreal breeding below about 7 °C monthly mean |
| μ_max | 0.1 /day | 0.1 /day | a stranded community largely dies within a month |
| s_ref, r_ref, Qs, Qr per class | fitted (Layer 1) | same peak rates as warm, about its own optimum | two sites cannot identify both guilds |

Shared defaults:

| Parameter | Value | Note |
|---|---|---|
| Sailing fraction δ | 0 | |
| Food α / K_chl | 0 / 1 mg/m³ | |
| Detachment, self-polishing paint | d = 0, v_d = 10 kn (soft), 30 kn (slime) | foul-release release speeds from the spec |
| ω, K17, K25 | fitted (Layer 1) | K_min = 10% |
| ρ (base plates left by in-water clean) | 0.2 | spec range 0–0.3; used by the Adland calibration |
| κ (paint damage per unit coverage removed) | fitted (Layer 3) | |

### 5.4 Events (`guild_model.rs`)

| Event | Slime / soft | Calcareous (hull) | Paint R | Propeller |
|---|---|---|---|---|
| DryDockFullBlast | 0 | C_g = D = 0 | R_out | k_p = clean; C_p = D_p = 0 |
| DryDockSpotRepair | 0 | C_g = D = 0 | R − recovery·(R − R_out), skipped if `just_docked` | as full blast |
| InWaterHullClean{ρ} | slime scaled so Σ F_sl = min(Σ F_sl, 0.75); soft 0 | D ← ρ·(ΣC_g + D); C_g = 0 | R + κ·(coverage removed)/100 | — |
| HullGrooming | as in-water clean | C_g = 0; D kept | unchanged | — |
| PropellerPolish | — | — | — | k_p = min(k_p, polished); C_p = D_p = 0 |

Properties:
- No event increases fouling roughness.
- Docking twice equals docking once.
- Cleaning an already-clean hull does nothing, because the damage scales with what was removed.
- The cleaning policies in `simulate` keep passing ρ = 0, so existing outputs do not change.

### 5.5 Solver (`integrate.rs`)

Each segment is cut into n = ⌈Δ/h_max⌉ equal steps, h = Δ/n, h_max = 0.25 day. A zero-length
segment is a no-op. The environment, and so every rate, is constant within a segment.

Each step uses integrating-factor (Lawson) RK4. The linear losses, including the transfer of
dead live-shell mass into D, are solved exactly by their closed-form flow, and classical RK4
handles the gain terms (settlement and growth) in that frame. The scheme is fourth order and
exact when there are no gains, so a 3 h freshwater half-life is exact even with quarter-day
steps. After each step, negative round-off is clamped to 0. The clamp never hides a NaN, and a
test checks it is never needed by more than 1e-12.

(An earlier draft used a Strang split: half-step decay, RK4, half-step decay. That is only
second order whenever gains and losses interact, which contradicts the fourth-order test in
Section 8.)

Paint and propeller drift are integrated exactly, as now. Results are deterministic: fixed
step count and no tolerance loops. Splitting a segment reproduces the unsplit result to
within 1e-6 (tested).

### 5.6 Observation

- **Hull k_s** = max(ks_from_rating(Σ F_sl + Σ F_so), ks_cov(C_w + C_c + D)).
- **ks_cov** is Uzun Eq. 18 with the branch above 50% replaced by a·(SC − 50) + k(50), so it is
  continuous at 50%. The slope a is fitted in Layer 2; Uzun's value is 80.
- **Propeller k_p,eff** = max(k_p, ks_cov(C_p,w + C_p,c + D_p)).
- **Penalty** = { fouling_resistance_fraction(op, k_s), coating_resistance_fraction(op, R_out, R),
  propeller_eta_loss(k_p,eff) }. A fresh or just-docked hull has k_s = 30 µm exactly and a
  power factor of exactly 1.

### 5.7 Interface

- `GuildModel::new(GuildParams) -> Result<GuildModel, FoulingError>` rejects:
  - non-finite or negative parameters;
  - guild windows out of order (T_opt < T_max, L_lo < L_hi, B_lo < B_hi);
  - K17, K25 or K_min outside (0, 100], or K_min > min(K17, K25);
  - a foul-release coating (`UnsupportedCoating`).
- It implements `FoulingModel`, so it works with `simulate`, `CleaningPolicy`,
  `mean_power_factor` and the power chain unchanged.
- Getters:
  - `slime_rating()`, `soft_rating()` (summed over guilds);
  - `live_coverage()`, `dead_coverage()`, `calcareous_coverage()`;
  - `hull_ks()`, `prop_ks()`, `ahr()`;
  - `guild(Guild) -> GuildState { slime, soft, live_coverage }`.
- `pub enum Guild { Warm, Cold }`.

## 6. Calibration (layered)

Each layer fits only its own parameters, in order. A later layer never changes an earlier
layer's parameters.

| Layer | Data | Fitted | Fixed | Pass criterion |
|---|---|---|---|---|
| 1 Growth | Uzun 2019 Tables 2–3, equatorial (25 °C) and Mediterranean (17 °C), idle at berth, 35 PSU, simulated day by day for 1,500 days | per class: s_ref, r_ref, Qs, Qr (warm guild); K17, K25, ω | guild windows, cold-guild peak rates, μ_max | worst-case gap at both sites: slime ≤ 1 point, soft ≤ 2 points, coverage ≤ 3 points |
| 2 Roughness → power | GloFouling 2022 5-year penalties without cleaning (+30% equatorial, +22% Mediterranean); GloFouling profile; 176 m tanker, frictional share 0.5, 14 kn | roughness scale λ: k_s = 30 + λ·(k_s,Uzun − 30) (amended 2026-11-07, see note) | Layer 1 | mean hull fouling fraction 0.30 ± 0.05 (equatorial) and 0.22 ± 0.05 (Mediterranean) |
| 3 Events | Adland dry-dock / in-water recovery ratio (about 2); Hundley destroyer sequence | κ | ρ = 0.2, Layers 1–2 | ratio in [1.6, 2.4] on the 3-year equatorial profile; Hundley pattern holds |

**Layer 2 amendment (2026-11-07, agreed with the user).** The planned knob, the coverage slope
above 50%, could not reach the equatorial target: 0.39 even at slope 0. The gap comes mainly
from Uzun's rating-to-roughness mapping (rating 70 → 1 mm → 36% penalty), which applies before
coverage passes 50%. Layer 2 therefore fits one scale λ on all Uzun-derived roughness, and the
slope stays at Uzun's 80. The result is λ = 0.519, giving 0.330 equatorial and 0.172
Mediterranean.

**Tool.** `crates/nautilus-ocean/examples/calibrate_fouling.rs`:
- runs the three layers with Nelder–Mead, using a fixed iteration count and fixed starting
  points, so it is deterministic and adds no dependency;
- fits Q10s and rates in log space;
- prints the fitted constants and the gaps achieved.

The constants are copied into `GuildParams::default()` with a comment naming the tool, the
date and the objective value. Tests check that the defaults meet every pass criterion; they do
not run the optimizer.

**If a criterion cannot be met**, implementation stops and the best fit is reported to the
user. Thresholds are never loosened silently.

## 7. Error handling

- Invalid segments and parameters fail at construction with `FoulingError`, as today.
- `Segment::with_chlorophyll` rejects negative or non-finite values with `InvalidSegment`.
- With valid inputs the solver cannot produce NaN; a property test over random world histories
  checks this.
- No `f64::min` or `f64::max` on state values in a way that could hide a NaN (the lesson from
  the exposure model's cap).

## 8. Testing

Every test carries a `///` comment saying what it checks and what the expected values come
from.

**Unit tests:**

| File | Checks |
|---|---|
| `params.rs` | every validation rule accepts the defaults and rejects each bad value; window ordering; capacity bounds |
| `rates.rs` | growth and settlement curves (Q10 rise, optimum, zero point); mortality zero inside the window and ramping outside; breeding edges; food exactly 1 at α = 0 or with no chlorophyll, and saturating; freshwater 50% kill at 3 h and > 99% at 24 h; detachment only above v_d; capacity line and clamps |
| `dynamics.rs` | derivative is zero for a clean hull with no settlement; growth stops at capacity; no negative derivative at zero state; live mass lost equals dead mass gained |
| `integrate.rs` | pure decay exact against exp; logistic growth against its closed form; 4th-order convergence (halving h cuts the error by 12–20×); segment splitting matches to 1e-6; bit-identical replay; clamping never exceeds 1e-12 |
| `guild_model.rs` | every event's effect; double docking equals single docking; cleaning a clean hull is a no-op; coating damage scales with coverage removed; fresh hull power factor exactly 1; per-guild breakdown; a runnable usage example in the docs |

**Integration tests (`crates/nautilus-ocean/tests/`):**
- **`fouling_validation.rs`:** the Layer 1, 2 and 3 pass criteria as required tests for
  `GuildModel`. The baseline characterisation tests stay.
- **`fouling_world_routes.rs`** (new), with synthetic world voyages:

  | Route | Expected behaviour |
  |---|---|
  | Singapore → Rotterdam (winter, 8 °C) | warm guild live amount falls > 90% within 45 days; dead coverage keeps hull k_s within 10% of its arrival value; the cold guild starts settling |
  | Rotterdam → Persian Gulf (33 °C) | cold guild dies off; warm guild settles |
  | Baltic stay (6 PSU, 60 days) | soft fouling → near 0; live calcareous survives (euryhaline) |
  | Panama transit (0.3 PSU, 12 h) | soft rating falls to 1/16 of its value (± 2% relative), from the 3 h half-life |
  | Southern-hemisphere seasonal route (Cape Town–Durban SST cycle) | settlement switches on and off with temperature, with no calendar input |
  | Same route with and without chlorophyll | α = 0 output identical to no-chlorophyll output; α > 0 never gives less fouling |
  | Fast liner, detachment enabled | soft fouling falls while sailing above v_d; unchanged with d = 0 |

- **Baseline comparison:** at constant 25 °C idle, `GuildModel` stays within the Layer 1
  tolerances of Uzun's equatorial curves.
- **Long random histories:** 20,000 random world segments (the ranges of the existing
  `random_segment` test helper, plus chlorophyll 0–20 mg/m³),
  with random events, for all three models. All invariants hold and every value stays finite.
- **Golden file:** `tests/golden/fouling_guild_world_route.csv`, a frozen trajectory on a fixed
  world route (`UPDATE_GOLDEN=1` regenerates it).
- **Existing suites:** all pass after the rename, with the existing golden numbers unchanged.

**Finishing checks:**
- a deliberate-bug sweep over the new files, where every planted bug must be caught;
- `cargo test --workspace`, `cargo clippy --workspace --all-targets -- -D warnings` and
  `cargo fmt --check`.

## 9. Documentation

- `docs/methodology/ocean-fouling.md`:
  - the guild model's equations (Section 5);
  - the parameter table with sources and "assumption — calibrate" flags;
  - the calibration results per layer;
  - updated model names.
- The module docs in `fouling/mod.rs` describe the three models and when to use each.

## 10. Constraints carried over

- `f64` only, std only, no new dependencies, no `unsafe`.
- Doc comments on every public item; clippy with `-D warnings`.
- Penalties are relative to the clean as-applied hull and are never added on top of measured
  fuel.
- Deterministic: fixed iteration counts, no `HashMap` iteration in outputs.
- No AI co-author trailers in commits.
