# Ocean: hull and propeller biofouling model

Code: `crates/nautilus-ocean/src/fouling/`. Spec: `research/reports/Biofouling growth dynamical model.md`.
Plan: `plans/0002-biofouling-dynamics-core.md`.

## What it computes

Given a ship's timeline (segments of constant forcing, plus maintenance events), the model
tracks a fouling state and turns it into a **shaft-power multiplier relative to the clean,
as-applied hull**. A freshly docked ship has a multiplier of exactly 1.

## Formulation: a hybrid (impulsive) dynamical system

- **Flow.** Between events, the state grows on *exposure clocks*. A clock does not count
  calendar time. It counts time weighted by how favourable conditions are for fouling:
  `dE = T(sst) · S(salinity) · gate(speed, port, distance to coast) · Δt`.
  Forcing is piecewise constant per `Segment`, so each update is closed form. There is no ODE
  solver, and a segment split into equal parts gives the same state as the whole segment
  (tested).
- **Jumps.** `MaintenanceEvent`s map the state to a new state. Every reset leaves a residual,
  no event ever increases fouling, and dry-docking twice in a row has the same effect as once
  (tested).
- **Observation.** The state maps to a `PowerPenalty { hull_fouling_fraction, coating_fraction,
  propeller_eta_loss }`, and then to
  `power_factor = (1 + r)·(1 + 0.1·r) / (1 − Δη₀/η₀)` with `r` = hull fouling + coating
  resistance fraction. The `0.1·r` propeller-loading term is back-calculated from Schultz (86%
  power vs 80% resistance). It is not a published coefficient.

### Shared response functions (`environment.rs`)

| Function | Form |
|---|---|
| Temperature `T` | `max(floor, Q10^((min(t, t_max) − 25)/10))`; `floor` below `t_min` |
| Salinity `S` | piecewise linear: 0 below `zero_below`, 0.5 at `half_at`, 1 at and above `full_at` |
| Settlement `A(v)` | `L(v)/L(0)`, `L(v) = 1/(1 + exp((v − v50)/width))` |
| Propagule `P` | 1 at berth/anchorage; at sea 1 within `near`, linear to `offshore_floor` at `far` |
| Gate | `x + δ·(1 − x)`, `x = A·P`; δ = 0 matches Uzun (growth only while idle) |

## Exposure model (`exposure_model.rs`)

One hull clock `E`, plus calendar drifts:

- `E += T·S·gate·Δ`, `AHR += drift·Δ/365.25`, `prop = min(cap, prop + drift·Δ/365.25)`
- `hull_fouling_fraction = min(0.86, 0.0012·E)`; coating fraction from AHR (below); propeller
  loss = `prop`.

Exposure model has a single clock, so it has no freshwater cull. Freshwater only pauses the clock
(S = 0).

| Event | E | AHR | prop |
|---|---|---|---|
| DryDockFullBlast | 0 | out-docking | 0 |
| DryDockSpotRepair | 0 | `AHR − recovery·(AHR − out)` | 0 |
| InWaterHullClean{ρ} | `min(E, E_res + ρ·(E − E_res))` | — | — |
| HullGrooming | `min(E, E_res)` | — | — |
| PropellerPolish | — | — | 0 |

## Age-clock model (`age_clock_model.rs`, `uzun.rs`)

This model uses the Uzun et al. (2019) equatorial curves for a self-polishing (SPC) coating,
rewritten as **age clocks**. A Gaussian rise in age, `FR = max·exp(−((A − t0)/τ)²)`, solves
the square-root-Gompertz ODE `dFR/dA = (2/τ)·FR·√ln(max/FR)`. The curve can therefore be
driven by a clock that runs faster or slower than calendar time, and it reproduces the
published curve exactly at 25 °C, full salinity, idle.

- Clocks (days at 25 °C): `a_slime += T_slime·S_soft·g·Δ`, `a_soft += T_soft·S_soft·g·Δ`,
  `a_calc += T_calc·S_hard·g·Δ` (same for the propeller calcareous clock).
- Calcareous coverage `SC` = main + fast logistic in `a_calc`. Its asymptotes depend on the
  exposure-weighted mean SST, linear between 17 °C (main 50%, fast 0%) and 25 °C (96.899%,
  3.5%), clamped.
- **Freshwater cull:** a freshwater *stay* (consecutive segments below 2 PSU) multiplies the
  slime and soft ratings by 0.5 once, when the stay first reaches 3 h, written back as clock
  values via the closed-form inverse. Because it acts per stay rather than per segment, the
  result does not depend on how the stay is split. Dead calcareous shells stay.
- Observation chain: `FR_tot = FR_slime + FR_soft` → `k_s` (Uzun Eq. 17); `SC` → `k_s`
  (Eq. 18, applies above 5%); hull `k_s` = the larger of the two → Granville `ΔC_F/C_F`.
  Propeller `k_p = max(drift roughness, k_s(SC_prop))` → Δη₀ table.
- The hull fouling fraction is taken relative to the freshly docked hull. Uzun's Gaussians
  have non-zero seeds, so a docked hull reads k_s ≈ 31 µm, not 30 µm; the penalty subtracts
  that docked value so a docked ship has a power factor of exactly 1.
- Foul-release coatings are refused (`UnsupportedCoating`). The Uzun parameters are SPC-only.

| Event | slime / soft | a_calc | coating | propeller |
|---|---|---|---|---|
| DryDockFullBlast | 0 / 0 | 0 (and SST accumulators) | AHR ← out | rough ← `prop_clean` (30 µm), a_prop_calc ← 0 |
| DryDockSpotRepair | 0 / 0 | 0 (and SST accumulators) | partial recovery | as full blast |
| InWaterHullClean{ρ} | slime ← min(a, age(FR 0.75)); soft ← 0 | age at `ρ·SC` | — | — |
| HullGrooming | as in-water clean | 0 | — | — |
| PropellerPolish | — | — | — | rough ← min(rough, `prop_polished` (40 µm)), a_prop_calc ← 0 |

FR 0.75 corresponds to k_s ≈ 40 µm, GloFouling's post-cleaning floor.

## Guild model (`guild_model.rs`, `rates.rs`, `dynamics.rs`, `integrate.rs`)

The exposure and age-clock models make every growth rate depend on the environment only.
Each site's curve is then a time-stretched equatorial curve, nothing dies except in fresh
water, and the barnacle plateau needs a running mean temperature. The guild model is
built for ships that move between climates. Its state is the fouling actually present,
and the rates depend on that state.

**State**, with g ∈ {warm, cold}. The hull is **two patches** that follow the same equations: the cleanable share η, which in-water cleaning and grooming reach, and the niche share 1 − η (sea chests, recesses, parts of the flat bottom), which only a dry-dock resets. Each patch carries:

- slime ratings F_sl,g (sum ≤ 20);
- soft ratings F_so,g (sum ≤ 50);
- live calcareous coverage C_g, plus dead shells D (live + dead ≤ 100 %);
- paint roughness R and propeller roughness k_p.

A fresh hull has no fouling, so its power factor is exactly 1.

**Equations.** Here x = settlement·supply gate and γ = growth gate, as for the other models:

```
dF_i,g/dt = σ_i(S)·φ(chl)·[ x·s_i,g(T)·b_g(T)·(1 − F_i/K_i) + γ·r_i,g(T)·F_i,g·√ln(K_i/F_i) ]
            − [ μ_g(T) + μ_fresh,i(S) + μ_brack,i(S) + δ_i(v) ]·F_i,g          (slime, soft)
Φ         = max(0, K_h(T) − C_w − C_c − D)
dC_g/dt   = σ_h(S)·φ(chl)·[ x·s_c,g(T)·b_g(T) + γ·r_c,g(T)·C_g ]·Φ/K_h − [ μ_g(T) + μ_fresh,h(S) ]·C_g
dD/dt     = Σ_g [ μ_g + μ_fresh,h ]·C_g − ω·D
```

| Rate | Form |
|---|---|
| Settlement s, growth r | Q10 rise to the guild optimum T_opt, linear fall to 0 at T_max |
| Breeding b_g | 1 inside the guild's breeding range, 1 °C ramps at the edges |
| Thermal death μ_g | 0 inside the lethal window, rising to 0.1/day within 3 °C outside it |
| Freshwater death μ_fresh | soft and slime: ln 2 per 3 h below 2 PSU; live barnacles: ln 2 per 3 days |
| Brackish death μ_brack | slime and soft: 0.1/day below 10 PSU |
| Food φ | 1 + α·chl/(chl + K); α = 0 by default |
| Detachment δ | d·max(0, v − v_d); d = 0 for self-polishing paint |
| Capacity K_h(T) | line through (17 °C, 49.9 %) and (25 °C, 96.8 %), clamped to [10, 100] % |

Guild windows (priors):

| Guild | Optimum / zero | Lethal window | Breeding range |
|---|---|---|---|
| Warm | 31 / 38 °C | 10–38 °C | ≥ 15 °C |
| Cold | 13 / 24 °C | −2–24 °C | ≤ 15 °C |

Both guilds share one thermal performance curve, referenced to 25 °C. The cold guild survives
and breeds in cold water; it does not grow faster there. The breeding ranges overlap, so there
is no temperature at which nothing settles. These priors were revised after the final review
(see Settled decisions).

**Observation:**

- per patch, k_s = 30 µm + λ·(k_s,Uzun − 30 µm), where k_s,Uzun is the larger of Uzun Eq. 17 (total
  rating) and Eq. 18 (live + dead coverage), and λ is the Layer 2 roughness scale;
- the hull fouling penalty is the **area-weighted mean of the two patches' resistance**. The
  `hull_ks()` getter reports the area-weighted mean roughness;
- the propeller uses the larger of its drift roughness and its scaled blade-shell roughness;
- the power chain is shared with the other models.

**Events:**

- **Dry-dock:** resets both patches. A duplicate record counts once.
- **In-water clean:** acts on the cleanable patch only. Slime goes to the 0.75 floor, soft
  fouling and live barnacles are removed, a fraction ρ of the shells above earlier base plates
  stays as new plates, and the paint gains 10 µm per 100 % coverage removed (area-weighted).
  The niche patch is untouched, so cleaning twice in a row is a no-op.
- **Grooming:** light brushing of the cleanable patch. Slime goes to the floor and soft fouling
  is removed; calcified barnacles (live and dead) are left, so grooming never beats a clean.
- **Propeller polish:** resets the blades.

**Solver.** Fixed-step integrating-factor RK4, with steps of at most 0.25 day. The linear
losses, including the transfer of dead barnacles into D, are solved exactly. The scheme is
fourth order and deterministic. After each segment, class totals are projected onto their hard
limits, so RK4 overshoot near saturation cannot leave a rating above its maximum.

### Calibration (layered, `examples/calibrate_fouling.rs`)

| Layer | Data | Fitted | Result |
|---|---|---|---|
| 1 Growth | Uzun 2019, 25 °C and 17 °C panels, 1,500 days | settlement (Q10 ≥ 1) and growth rates and their Q10s per class; capacities | worst gaps: slime 0.08, soft 1.42 rating points; coverage 0.07 percentage points (limits 1 / 2 / 3) |
| 2 Roughness → power | GloFouling 5-year penalties (+30 % equatorial, +22 % Mediterranean) | roughness scale λ = 0.519 | 0.3305 and 0.1707 (limits ± 0.05) |
| 3 Events | Adland dry-dock / in-water recovery ratio of about 2 | cleanable share η = 0.392 | 2.000 (limits 1.6–2.4) |

### Settled decisions (do not re-litigate)

All of these were decided with the user on 2026-11-07.

1. **Fast barnacle term:** equatorial coverage is scored against Uzun's main logistic only. The
   3.5 % fast term is unconfirmed (research spec open item 1) and stays below the 5 % threshold
   of Eq. 18.
2. **Layer 2 knob:** one roughness scale λ on all Uzun-derived roughness. The coverage slope
   above 50 % could not close the gap; Uzun's rating mapping (a rating of 70 maps to 1 mm) is
   the cause.
3. **Partial cleaning:** two hull patches with area-weighted resistance, instead of
   hull-average coverage.
4. **Cleanable share η = 0.39 is the default.** It represents routine commercial cleans
   (Adland's 2:1). A thorough clean like Hundley's destroyer (+70 % → +9 %) needs about 0.95.
   One default cannot satisfy both: the Hundley test runs with η = 1, and users set η per fleet.
5. **Priors after the final review:**
   - overlapping breeding ranges (no 13–14 °C dead band);
   - one shared thermal curve for both guilds, so cold water fouls less;
   - warm guild optimum 31 °C and limit 38 °C, so Gulf and Red Sea heat still fouls;
   - settlement Q10 ≥ 1;
   - grooming leaves barnacles.

**Caveats:**

- The fitted calcareous *settlement* Q10 is a fit artefact (≈ 3.6e9). In practice barnacles
  barely settle below about 22 °C, and settlement is capped at the 31 °C optimum.
- Shell weathering (0.001/day) is a prior: it cannot be identified from the panel sites.
- Layer 2's Mediterranean result (0.1707) sits only 0.0007 inside its band.

### Which model to use

- **Exposure model:** quick estimates and GloFouling-style trigger points.
- **Age-clock model:** the published Uzun baseline at a fixed climate.
- **Guild model:** routes that cross climates or pass through fresh or brackish water.

## Roughness → resistance (`granville.rs`, `power.rs`)

- Granville similarity law with the **Nikuradse-type** roughness function, κ = 0.42, against the
  Schoenherr line. All solvers are 64-step bisections, so results are deterministic. The
  Schultz/Uzun `k_s` scale is tied to this function. A Colebrook-type function would give
  +27% at 30 µm.
- `GranvilleTable`: 32 × 48 grid over log10 Re ∈ [7, 10] and log10(k_s/L) ∈ [−8, −3], bilinear,
  clamped at the edges (never extrapolated). Random-point accuracy vs the direct solver is
  within 0.6 percentage points (tested). The grid is solved once per process.
- Hull fouling fraction = `frictional_share · max(0, G(k_s) − G(30 µm))`. The penalty is
  measured from the as-applied reference.
  Roughness above 20 mm (beyond every source used) is treated as 20 mm.
- The operating speed must be finite and positive (`OperatingPoint::new` rejects anything
  else); friction coefficients are undefined for a ship that is not moving.
- Coating fraction = `share · 0.044·((AHR/L)^(1/3) − (AHR_out/L)^(1/3)) / C_F`. This is the
  difference of two ITTC roughness allowances, so it is relative to the out-docking AHR.

## Input checks and duplicate records

- Segments with negative or non-finite duration, speed, salinity or distance, or a
  non-finite temperature, are rejected.
- Model parameters are rejected if they would make the power factor NaN, infinite or
  negative: non-positive Q10 or settlement width, a propeller-loss cap of 1 or more,
  fractions outside [0, 1], or salinity breakpoints and coast distances out of order.
- Cleaning policies are built with `CleaningPolicy::triggered`, which requires a finite,
  non-negative trigger and a finite, positive interval.
- A spot-repair docking recorded twice with no time in between counts once, so a
  duplicated docking record cannot compound the paint recovery.

## Segmentation

`mean_power_factor` weights each segment by the penalty at its end, and `simulate` checks
cleaning triggers only at segment ends. Results are accurate only when segments are short
relative to how fast fouling changes (days, not weeks). For example, the exposure model at berth for 1000
days gives a mean power factor of 2.20 as one segment and 1.73 as 1000 one-day segments. Split
long segments (for example, multi-week AIS berth stays) before simulating.

## Parameters

| Parameter | Value | Source | Assumption — calibrate? |
|---|---|---|---|
| Reference temperature | 25 °C | Uzun equatorial site (per GloFouling) | |
| Q10 the exposure model | 2.0 | spec | yes |
| T floor / t_min / t_max | 0.1 / 12 °C / 32 °C | spec | yes — the factor jumps at `t_min` (0.1 → 0.41 for Q10 2) by design |
| Euryhaline salinity (0 / 0.5 / 1) | 2 / 5 / 20 PSU | spec | yes |
| Stenohaline salinity | 10 / 17.5 / 25 PSU | spec (soft drop-out 10–25 PSU) | yes |
| Settlement v50 / width | 1.25 kn / 0.2547 kn | spec | |
| Propagule near / far / floor | 15 km / 100 km / 0.02 | spec | yes |
| Sailing fraction δ | 0 | Uzun-compatible | yes (open item 8) |
| Exposure model slope | 0.0012 per exposure day (±50%) | derived from GloFouling | yes |
| Exposure model hull cap | 0.86 | Schultz heavy calcareous | |
| AHR out-docking | 120 µm | spec | |
| AHR drift | 20 µm/yr (range 3–40) | spec | yes |
| Exposure model propeller drift / cap | 0.026/yr / 0.32 | van Ballegooijen via GloFouling | |
| Residual exposure | 13 d (≈ +1.6% at 40 µm, tanker) | derived | yes |
| Spot-repair recovery | 0.5 | spec | yes |
| Slime / soft / calcareous Q10 | 4.1554 / 1.5406 / 2.0417 | fitted to Uzun's two sites | site SSTs: open item 2 |
| Cold asymptotes (main / fast) | 50% / 0% | two-site placeholder | yes |
| Fresh threshold / min duration / cull | 2 PSU / 3 h / 0.5 (range 0–0.5) | spec | yes |
| Propeller clean / polished / drift / cap | 30 / 40 µm / 100 µm/yr / 300 µm | Rubert scale | yes |
| Propeller Δη₀ table | (30, 0) … (10 000 µm, 0.222) | ITTC strip, 7.9 m 5-blade | yes (open item 6) |
| Cleaning policies | proactive 20% then 6 mo; reactive 40% then 12 mo | GloFouling | |

### Deviation: calcareous Q10

The spec gives 2.3, derived from the calcareous **rating** `t0`. The calcareous rating is not
used here: the spec excludes it from `FR_tot`, and calcareous fouling enters only through
coverage. The model therefore fits Q10 to the **coverage logistic midpoints**, 393.12 d
(equatorial, 25 °C) → 695.86 d (Mediterranean, 17 °C), which gives
`Q10 = (695.86/393.12)^(10/8) = 2.0417`.

## Double counting

The penalty describes the modelled effect of fouling on a *modelled* clean-hull power. **It is
never added on top of measured fuel**, which already contains the fouling. When the model feeds
the ocean pipeline's modelled tier, the fleet-mean penalty must stay near the IMO estimate of
7–9%. A materially higher fleet mean points to a calibration error, not a finding.

## Reference outputs

`cargo run -p nautilus-ocean --example fouling_scenario` gives these 5-year mean power factors
for the GloFouling equatorial profile (176 m tanker, frictional share 0.5, 14 kn):

| | none | reactive | proactive |
|---|---|---|---|
| Exposure model | 1.513 | 1.229 | 1.101 |
| Age-clock model | 1.736 | 1.304 | 1.152 |

This profile is the GloFouling worst case: 25 °C all year and 30% idle time. It is not a fleet
average.

## Open items carried from the spec

1. Uzun Eq. 4 calcareous coefficient; the evaluation point of the derivatives in Eqs. 6 and 8;
   whether the fast equatorial logistic term is applied. *Resolve from the typeset article or
   the Uzun PhD thesis.*
2. Uzun test-site SSTs (17/25 °C are GloFouling's values). All the Q10 values derived here
   depend on them.
3. The roughness function and constant B used by Schultz 2007. These would confirm the
   Nikuradse pairing.
6. Farkas 2021 relation between k_s and Δη₀. It would give a better propeller observation than
   the ITTC strip table.
8. Growth of attached fouling while sailing (δ), and the barnacle growth rate k_B(T). This is
   the largest structural uncertainty.
9. Residual calcareous fraction ρ after in-water cleaning, and the size of coating damage.
10. The rule for combining coating ΔC_F (ITTC) with fouling ΔC_F (Granville). The current sum
    may double-count the as-applied roughness.
11. Frictional share by ship type and Froude number.
