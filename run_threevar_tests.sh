#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$ROOT/build"
LOG_FILE="$LOG_DIR/threevar_tests.log"

mkdir -p "$LOG_DIR"

export CI=1
export GTEST_COLOR=0
export LD_LIBRARY_PATH="$ROOT/external/flint-install/lib:${LD_LIBRARY_PATH:-}"
export RUR_CRT_DIAG=${RUR_CRT_DIAG:-1}

cd "$ROOT"

timeout 120 "$ROOT/build/rur_tests" --gtest_brief=1 --gtest_filter="*ThreeVariable*" \
  > "$LOG_FILE" 2>&1 || true

echo "Exit: $?"
echo "--- summary (last 80 lines) ---"
tail -n 80 "$LOG_FILE" || true
echo "--- grep summary ---"
grep -E 'FAILED|PASSED|\[ *[A-Z]+ *\]| tests from | Failure' "$LOG_FILE" | tail -n 40 || true
echo "Full log: $LOG_FILE"


