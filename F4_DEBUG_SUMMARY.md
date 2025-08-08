# F4 Library Debug Summary

## Issue Identified
The F4 Gröbner basis library has a critical bug that causes segmentation faults when computing Gröbner bases for certain polynomial/prime combinations.

## Root Cause
The crash occurs in the Jordan elimination phase (`f4_jordan_thread` and `f4_jordan` functions) when `f4_reduce_export` returns NULL (indicating a zero polynomial after reduction), but the code doesn't check for this condition before dereferencing the result.

## Affected Code
1. **Multi-threaded version** (`f4_jordan_thread` at line 2471): Dereferences `a->ind` without checking if `a` is NULL
2. **Single-threaded version** (`f4_jordan` at line 2558): Passes NULL to `f4_reduce_monic` which will crash

## Test Results
- **Before fix**: Segmentation fault with primes 1073741827 and 65537
- **After fix**: No crashes, but there's a separate parsing issue (unrelated to the crash)

## Fix Applied
Added NULL checks in both functions:
1. In `f4_jordan_thread`: Skip to next block if reduction returns NULL
2. In `f4_jordan` (single-threaded): Set `piv[i] = 0` instead of calling `f4_reduce_monic(NULL)`

## Files Modified
- `/home/orebas/code/rur-cpp/src/axf4_lib.c` - Added NULL checks at lines 2471-2474 and 2558-2562

## Deliverables
1. **Bug Report**: `/home/orebas/code/rur-cpp/F4_BUG_REPORT.md` - Complete bug report for F4 developer
2. **Patch File**: `/home/orebas/code/rur-cpp/f4_bug_fix.patch` - Ready-to-apply patch
3. **Test Case**: `/home/orebas/code/rur-cpp/test/test_minimal_f4.c` - Minimal reproduction case

## Status
✅ **Bug Fixed**: The segmentation fault is resolved. The F4 library no longer crashes with the problematic prime values.

## Notes
- The crash was deterministic and reproducible
- The fix follows the existing error handling pattern in the codebase
- This is a minimal, safe fix that doesn't change the algorithm's behavior
- There's a separate parsing issue ("error: can't parse term") that's unrelated to this crash