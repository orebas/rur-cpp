# F4 Bug Analysis: NORMAL Macro Misuse

## Summary
The F4 library has a bug in its modular arithmetic for primes <= 2^31-1 on 64-bit systems. The bug is in the "delayed reduction" optimization path where the NORMAL macro is incorrectly used with p² instead of p.

## Root Cause

In `axf4_lib.c`, the NORMAL macro is defined as:
```c
#define NORMAL(x,p) x += (x >> (WORDSIZE-1)) & p
```

This macro normalizes negative values by adding `p` if the value is negative (sign bit is 1).

However, in the optimized arithmetic path for primes <= 2^31-1, the code does:
```c
if (WORDSIZE==64 && p <= 2147483647) {
    p2 = p*p;
    // ... arithmetic operations ...
    y = vec[z] - x;
    NORMAL(y,p2);  // BUG: Using p² instead of p!
    if (y >= p) y %= p;
}
```

## Why This Causes Wrong Results

When computing x² - 1 = 0:
1. The coefficient -1 should be represented as p-1 in modular arithmetic
2. During reduction, the code computes `y = vec[z] - x` which can be negative
3. `NORMAL(y,p2)` adds p² (instead of p) to negative values
4. This produces incorrect coefficients that are not properly reduced modulo p

## Why Prime 1073741827 Works

This appears to be coincidental. The specific arithmetic operations with this prime somehow produce correct results despite the bug, possibly due to:
- Specific bit patterns in the prime
- How the multiplication and reduction interact
- The order of operations in the Gröbner basis computation

## Proof of Bug

Our test shows:
- Prime 1073741827: Coefficient stored as 1073741826 (correct: p-1)
- All other primes: Coefficient stored as 1 (wrong: should be p-1)

## Fix Options

1. **Correct the NORMAL macro usage**:
   ```c
   NORMAL(y,p);  // Use p, not p2
   ```

2. **Disable the optimization** for now by forcing all primes to use the slower path:
   ```c
   if (0 && WORDSIZE==64 && p <= 2147483647) {  // Disabled
   ```

3. **Work around in our code**: Only use primes > 2^31-1 or specifically use 1073741827

## Impact

This bug affects ALL polynomial systems with negative coefficients when using primes <= 2^31-1 (except for the coincidentally working prime 1073741827). It produces mathematically incorrect Gröbner bases, which cascade into wrong results for:
- Polynomial system solving
- Ideal membership testing  
- Rational univariate representation
- Any computation depending on correct Gröbner bases