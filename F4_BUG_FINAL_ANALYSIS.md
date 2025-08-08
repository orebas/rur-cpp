# F4 Bug - Final Analysis and Fix

## Executive Summary

The F4 Gröbner basis library had a classic C programming bug: an uninitialized variable. This bug caused inconsistent polynomial parsing and led to incorrect Gröbner bases for most prime moduli. The fix is a simple one-line change.

## The Real Bug

In the `getmon` function in `axf4_lib.c` (line ~2860), the variable `t` was declared but not initialized:

```c
INT getmon(char *s, INT *c, int *l)
{
    INT z[MAXVARS] = {0};
    INT b, e, m, n, p, t;  // 't' is not initialized!
    int i, j, k;
    char *v;

    n = f4_nvars;
    p = f4_prime; 
    *c = 1; i = 0;  // Should be: *c = 1; i = 0; t = +1;
```

The variable `t` stores the sign of a polynomial term (+1 or -1). It's only set when a term explicitly starts with '+' or '-'. For polynomials like `"1*x^2-1"`, the first term `1*x^2` has no leading sign, leaving `t` uninitialized with stack garbage.

## Why This Caused Inconsistent Results

When parsing coefficients, the code uses `t` to determine the sign:
```c
if (j > 0 && t==+1) {
    *c = b % p;
    NORMAL(*c, p);
}
if (j > 0 && t==-1) *c = (p-b) % p;
```

With uninitialized `t`, the behavior was undefined and depended on whatever value happened to be on the stack. This explains why:
- Different primes gave different results (different call patterns → different stack states)
- Some primes "worked" by chance (stack happened to contain a good value)
- The results were inconsistent between runs

## The Fix

Change line 2866 in `axf4_lib.c`:

**From:**
```c
*c = 1; i = 0;
```

**To:**
```c
*c = 1; i = 0; t = +1;
```

This ensures `t` defaults to +1 (positive) for terms without an explicit sign.

## Verification

After applying the fix, all tested primes now produce correct and consistent results:
- Prime 1073741827: ✓ Correct (y^2 + 1073741826)
- Prime 1048573: ✓ Correct (y^2 + 1048572)  
- Prime 268435399: ✓ Correct (y^2 + 268435398)
- Prime 65537: ✓ Correct (y^2 + 65536)

All coefficients now correctly represent -1 as (p-1) in modular arithmetic.

## Lessons Learned

1. **The obvious explanation isn't always right**: We initially suspected complex mathematical issues with modular arithmetic and the NORMAL macro, but it was just an uninitialized variable.

2. **Undefined behavior is unpredictable**: The bug manifested differently for different primes due to varying stack states, making it look like a deep algorithmic issue.

3. **Classic C bugs still happen**: Even in well-established libraries, simple bugs like uninitialized variables can lurk for years.

4. **Collaborative debugging works**: Multiple perspectives (yours, mine, Gemini's, O3's, and the other developer's) were crucial in finding the real issue.

## Credit

The other developer deserves full credit for identifying this bug. Their insight to check for initialization issues rather than mathematical problems was spot-on.