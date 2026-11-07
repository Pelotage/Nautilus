#!/usr/bin/env bash
# One-shot, re-runnable setup for building and running Nautilus on Ubuntu.
#
#   ./scripts/setup-ubuntu.sh              core toolchain + dev tools + Python sidecar
#   ./scripts/setup-ubuntu.sh --with-fuzz  also install nightly + cargo-fuzz (fuzzing only)
#   ./scripts/setup-ubuntu.sh --no-python  skip uv / Python sidecar setup
#
# Uses sudo only for apt. Everything else installs into your home directory.
set -euo pipefail

WITH_FUZZ=0
WITH_PYTHON=1
for arg in "$@"; do
  case "$arg" in
    --with-fuzz) WITH_FUZZ=1 ;;
    --no-python) WITH_PYTHON=0 ;;
    -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

step() { printf '\n==> %s\n' "$*"; }

step "System packages (apt)"
sudo apt-get update
sudo apt-get install -y \
  build-essential pkg-config libssl-dev curl git ca-certificates \
  python3 python3-venv python3-dev

step "Rust toolchain (rustup)"
if ! command -v rustup >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain none
fi
# shellcheck disable=SC1091
[ -f "$HOME/.cargo/env" ] && . "$HOME/.cargo/env"

# rust-toolchain.toml pins the exact version and components; this installs them.
rustup show active-toolchain >/dev/null 2>&1 || rustup toolchain install
rustc --version
cargo --version

step "Cargo dev tools (pinned by --locked)"
for tool in cargo-deny cargo-audit cargo-insta; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    cargo install --locked "$tool"
  else
    echo "$tool already installed"
  fi
done

if [ "$WITH_FUZZ" -eq 1 ]; then
  step "Fuzzing (needs nightly; used only for 'cargo +nightly fuzz')"
  rustup toolchain install nightly --profile minimal
  command -v cargo-fuzz >/dev/null 2>&1 || cargo install --locked cargo-fuzz
fi

if [ "$WITH_PYTHON" -eq 1 ]; then
  step "Python sidecar (uv)"
  if ! command -v uv >/dev/null 2>&1; then
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.local/bin:$PATH"
  fi
  uv --version
  (cd research && uv sync)
fi

step "Verifying the workspace builds"
cargo build --workspace
cargo test --workspace

cat <<'EOF'

Setup complete. If this is a fresh shell, load Rust and uv into PATH with:
  . "$HOME/.cargo/env"; export PATH="$HOME/.local/bin:$PATH"
Then try:  cargo run -p nautilus-cli
EOF
