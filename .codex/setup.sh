#!/usr/bin/env bash
# Setup script for the poasta repository.
# Installs Rust if not present and builds the project in release mode.
set -euo pipefail

if ! command -v cargo >/dev/null; then
    echo "Installing Rust toolchain" >&2
    curl https://sh.rustup.rs -sSf | sh -s -- -y
    source "$HOME/.cargo/env"
fi

cargo build --release
