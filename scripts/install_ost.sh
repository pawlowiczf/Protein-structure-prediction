#!/usr/bin/env bash
# Usage: ./scripts/install_ost.sh
# Installs OpenStructure into tools/ost-env using a local micromamba binary.
# No root, no conda init, nothing outside the repo. Safe to re-run.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS="$ROOT/tools"
MAMBA="$TOOLS/micromamba"
ENV_DIR="$TOOLS/ost-env"

mkdir -p "$TOOLS"

if [[ -x "$ENV_DIR/bin/ost" ]]; then
    echo "OpenStructure already installed: $ENV_DIR/bin/ost"
    exit 0
fi

if [[ ! -x "$MAMBA" ]]; then
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64)  PLATFORM=linux-64 ;;
        Linux-aarch64) PLATFORM=linux-aarch64 ;;
        Darwin-arm64)  PLATFORM=osx-arm64 ;;
        Darwin-x86_64) PLATFORM=osx-64 ;;
        *) echo "Unsupported platform: $(uname -s)-$(uname -m)" >&2; exit 1 ;;
    esac
    echo "Downloading micromamba ($PLATFORM)..."
    curl -fsSL "https://micro.mamba.pm/api/micromamba/$PLATFORM/latest" \
        | tar -xj -C "$TOOLS" --strip-components=1 bin/micromamba
fi

export MAMBA_ROOT_PREFIX="$TOOLS/mamba-root"
"$MAMBA" create -y -p "$ENV_DIR" -f "$ROOT/environment-ost.yml"
echo "Installed: $ENV_DIR/bin/ost"
"$ENV_DIR/bin/ost" --help | head -3
