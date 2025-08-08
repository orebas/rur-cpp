Developer notes: non-interactive builds/tests in Cursor

Background
Some Cursor shells may appear to hang after a long-running command with large output. Use one-shot, non-interactive patterns and log to files to avoid blocking and to keep the UI responsive.

Patterns
- Build (one-shot):
  cd /home/orebas/code/rur-cpp/build && stdbuf -oL -eL ninja -j4 | cat; exit

- Run a focused test with compact output (log to file, print summary tail):
  cd /home/orebas/code/rur-cpp && \
  CI=1 GTEST_COLOR=0 \
  LD_LIBRARY_PATH=/home/orebas/code/rur-cpp/external/flint-install/lib:$LD_LIBRARY_PATH \
  timeout 120 ./build/rur_tests --gtest_brief=1 --gtest_filter="*ThreeVariable*" \
  > build/threevar_tests.log 2>&1; \
  status=$?; echo "Exit: $status"; \
  echo "--- summary (last 80 lines) ---"; tail -n 80 build/threevar_tests.log; \
  echo "--- grep summary ---"; \
  grep -E 'FAILED|PASSED|\\[ *[A-Z]+ *\\]| tests from | Failure' build/threevar_tests.log | tail -n 40; \
  echo "Full log: build/threevar_tests.log"; \
  exit

Tips
- Lower tail size for extra-quiet summaries (e.g., tail -n 40).
- Use script -qec 'CMD' /dev/null to force a non-interactive PTY if needed.
- Replace the shell with the command using exec to ensure auto-exit: exec CMD


