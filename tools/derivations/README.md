# Derivations

The scripts that produced the constants and test expectations hard-coded in the Rust fouling
model (`crates/nautilus-ocean/src/fouling/`), kept so every number can be reproduced.

| Script | What it produces | Where the values live |
|---|---|---|
| `derive_constants.py` | Recomputes and **checks** the hand-derived values: the age-clock Q10s, the exposure-model slope and residual, the settlement width, the Uzun curve test vectors, the ITTC allowance, the coating penalty, and the analytic starts for the guild calibration | `params.rs`, `environment.rs`, the tests in `uzun.rs`, `granville.rs` and `power.rs`, and `examples/calibrate_fouling.rs` |
| `granville_prototype.py` | The friction-increase table behind the Granville test vectors | `granville.rs` and `power.rs` tests |
| `pdf_text/` | Text extractors used to read the source tables (Uzun 2019, GloFouling 2022) out of the papers' PDFs. `ext.py` needs `pypdf`; the other two use the standard library only | the inputs to all of the above |

Run `python3 -I tools/derivations/derive_constants.py`. It prints each value next to the Rust
value and exits non-zero on any mismatch.

**Fitted constants are not derived in Python.** The guild model's Layer 1–3 values in
`GuildParams::default()` come from the Rust calibration tool:

```
cargo run --release -p nautilus-ocean --example calibrate_fouling -- layer1   # then layer2, layer3
```

The doc comment on each fitted field records the date and the score achieved.
