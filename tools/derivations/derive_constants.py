"""Recomputes the hand-derived constants and test expectations written into the Rust code,
and checks each against the value hard-coded there. Standard library only.

    python3 -I tools/derivations/derive_constants.py      # exits non-zero on any mismatch
"""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from granville_prototype import cf_rough, cf_s, dU_nik

checks = []
def check(name, rust_location, got, want, tol):
    ok = abs(got - want) <= tol
    checks.append(ok)
    print(f"{'OK ' if ok else 'BAD'} {name:<46} {got:<22.10g} rust {want:<14.10g} [{rust_location}]")

# Age-clock Q10s: ratio of Uzun's Mediterranean to equatorial timing over 8 °C (17 vs 25 °C).
q = lambda med, eq: (med / eq) ** (10 / 8)
check("slime Q10 (t0 271.9 / 87)", "params.rs AgeClockParams", q(271.9, 87.0), 4.1554, 5e-4)
check("soft Q10 (t0 383.5 / 271.4)", "params.rs AgeClockParams", q(383.5, 271.4), 1.5406, 5e-4)
check("calcareous Q10 (midpoints 695.86 / 393.12)", "params.rs AgeClockParams", q(695.86, 393.12), 2.0417, 5e-4)

# Exposure-model slope: GloFouling 20 % at 1.5 years with 556 idle days in 5 years.
idle_per_day = 556 / (5 * 365.25)
check("exposure slope per idle day", "params.rs ExposureParams", 0.20 / (547.5 * idle_per_day), 0.0012, 2e-5)

# Settlement gate: logistic midpoint 1.25 kn, width such that the logistic is 0.95 at 0.5 kn
# and 0.05 at 2 kn (spec: "1 below 0.5 kn → 0 by 2 kn"): (2 − 1.25)/w = ln 19.
check("settlement width 0.75 / ln 19", "environment.rs GateParams", 0.75 / math.log(19), 0.2547, 5e-5)

# Residual exposure: 13 idle days ≈ the +1.6 % hull penalty of GloFouling's 40 µm floor (tanker).
re = 14 * 1852 / 3600 * 176 / 1.19e-6
g = lambda ks: cf_rough(re, ks, 176, dU_nik) / cf_s(re) - 1
floor_penalty = 0.5 * (g(40e-6) - g(30e-6))
check("residual exposure (days)", "params.rs ExposureParams", floor_penalty / 0.0012, 13.0, 1.0)

# Uzun 2019 curve test vectors (uzun.rs tests).
gauss = lambda a, mx, t0, tau: mx * math.exp(-(((a - t0) / tau) ** 2)) if a < t0 else mx
check("EQ slime seed", "uzun.rs tests", gauss(0, 20, 87, 37.08), 0.081326, 1e-5)
check("EQ slime rating at 50 d", "uzun.rs tests", gauss(50, 20, 87, 37.08), 7.38937, 1e-4)
check("EQ soft rating at 200 d", "uzun.rs tests", gauss(200, 50, 271.4, 73.11), 19.2643, 1e-3)
check("EQ slime age for rating 0.75", "uzun.rs tests", 87 - 37.08 * math.sqrt(math.log(20 / 0.75)), 19.8103, 1e-3)
cov = lambda a: 96.899 / (1 + math.exp(16 - 0.0407 * a)) + 3.5 / (1 + math.exp(10.32 - 0.7759 * a))
check("EQ coverage at 100 d", "uzun.rs tests", cov(100), 3.50064, 1e-4)
check("k_s from coverage 25 % (Eq. 18)", "uzun.rs tests", 2.4669 * 625 - 24.84 * 25 + 1065.7, 1986.5125, 1e-3)

# ITTC-78 allowance to full precision (granville.rs ittc_allowance_full_precision).
ittc = lambda ks: 0.044 * ((ks / 176) ** (1 / 3) - 10 * re ** (-1 / 3)) + 0.000125
check("ITTC allowance at 150 µm", "granville.rs tests", ittc(150e-6), 0.00011133652775339888, 1e-15)
check("ITTC allowance at 1000 µm", "granville.rs tests", ittc(1000e-6), 0.0004793097972998913, 1e-15)

# Coating penalty 120 → 150 µm on a 176 m ship at 14 kn, full frictional share (power.rs tests).
coat = 0.044 * ((150e-6 / 176) ** (1 / 3) - (120e-6 / 176) ** (1 / 3)) / cf_s(re)
check("coating fraction 120 → 150 µm", "power.rs tests", coat, 0.0197, 5e-4)

# Guild model: analytic starts for the calcareous Layer 1 fit (examples/calibrate_fouling.rs).
s25 = 0.0407 * 98.6 * math.exp(-0.0407 * 393.12)
s17 = 0.04715 * 50 * math.exp(-0.04715 * 695.86)
check("calcareous settle start s(25 °C)", "calibrate_fouling.rs", s25, 4.6e-7, 1e-8)
check("calcareous growth Q10 start", "calibrate_fouling.rs", (0.0407 / 0.04715) ** (10 / 8), 0.83, 5e-3)
check("slime growth start 2/τ", "calibrate_fouling.rs", 2 / 37.08, 0.05394, 1e-4)

print(f"\n{sum(checks)}/{len(checks)} constants reproduced")
sys.exit(0 if all(checks) else 1)
