#!/usr/bin/env bash
# Usage: ./scripts/install_tmscore.sh
# Downloads TMscore.cpp and builds tools/TMscore. Safe to re-run.

set -euo pipefail

TOOLS_DIR="$(cd "$(dirname "$0")/.." && pwd)/tools"
mkdir -p "$TOOLS_DIR"
if [[ -x "$TOOLS_DIR/TMscore" ]]; then
    echo "TMscore already built: $TOOLS_DIR/TMscore"
    exit 0
fi

curl -fsSL -o "$TOOLS_DIR/TMscore.cpp" https://zhanggroup.org/TM-score/TMscore.cpp
g++ -O3 -ffast-math -o "$TOOLS_DIR/TMscore" "$TOOLS_DIR/TMscore.cpp" -lm
echo "Built $TOOLS_DIR/TMscore"
