# Nautilus

Freight emissions estimation engine (ISO 14083): traceable, reproducible, point-to-point
estimates across air, ocean, road and parcel legs.

Architecture and rationale: `plans/0001-engine-architecture.md` (local, gitignored).

## Quick start (Ubuntu)

```bash
./scripts/setup-ubuntu.sh        # one-time: installs everything needed
cargo build --workspace
cargo run -p nautilus-cli
cargo test --workspace
```

See `scripts/setup-ubuntu.sh --help` for options.
