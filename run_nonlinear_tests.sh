#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$ROOT/build"
mkdir -p "$LOG_DIR"

export CI=1
export GTEST_COLOR=0
export LD_LIBRARY_PATH="$ROOT/external/flint-install/lib:${LD_LIBRARY_PATH:-}"
export RUR_CRT_DIAG=${RUR_CRT_DIAG:-1}

cd "$ROOT"

declare -a TESTS=(
  "NonlinearSystemTests.ParabolaLineIntersection"
  "NonlinearSystemTests.TwoCirclesIntersection"
  "NonlinearSystemTests.QuadraticSystem"
)

for t in "${TESTS[@]}"; do
  LOG_FILE="$LOG_DIR/${t//[^A-Za-z0-9_]/_}.log"
  echo "=== Running $t ===" | tee "$LOG_FILE"
  timeout 120 "$ROOT/build/rur_tests" --gtest_brief=1 --gtest_filter="$t" \
    >> "$LOG_FILE" 2>&1 || true
  echo "Exit: $?" | tee -a "$LOG_FILE"
  echo "--- summary (last 60 lines) ---" | tee -a "$LOG_FILE"
  tail -n 60 "$LOG_FILE" || true
  echo "Log: $LOG_FILE"
done

echo "--- CRT diag tail ---"
tail -n 80 "$LOG_DIR/crt_diag.log" 2>/dev/null || true


