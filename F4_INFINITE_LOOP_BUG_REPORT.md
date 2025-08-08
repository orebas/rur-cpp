# F4 Gröbner Basis Algorithm Infinite Loop Bug Report

## Summary
The F4 algorithm implementation in `axf4_lib.c` enters an infinite loop for certain prime moduli when computing Gröbner bases. The algorithm gets stuck repeatedly processing the same pair without making progress, outputting "0 x 0 with 0 non-zero, -nan per row" indefinitely.

## Environment
- File: `/home/orebas/code/rur-cpp/src/axf4_lib.c`
- Wrapper: `/home/orebas/code/rur-cpp/src/axf4_wrapper.c`
- Platform: Linux 6.6.87.2-microsoft-standard-WSL2

## Bug Details

### Symptoms
1. **Infinite Loop**: The algorithm enters an infinite loop with the step counter incrementing indefinitely
2. **Zero Matrix**: Output shows "0 x 0 with 0 non-zero, -nan per row" indicating no polynomials in the matrix
3. **No Progress**: `new=0 basis=1 extra=1` shows no new polynomials are being generated
4. **NaN Output**: Division by zero when computing average non-zeros per row (0/0)

### Affected Primes
- **Working**: 1073741827 (and likely other primes < 2^31)
- **Failing**: 2147483629, 2147483587, 2147483579, 2147483563, 2147483549 (all near 2^31)

### Test Case
```
System: x^2 + y^2 - 1 = 0 (circle)
        x - y = 0 (line)
Variables: x, y
```

### Root Cause Analysis

The issue appears to be related to integer overflow or precision problems when the prime modulus approaches or exceeds 2^31. The key observations:

1. **Prime Size Threshold**: All failing primes are close to 2^31 (2147483648)
2. **Arithmetic Operations**: The code has special handling for primes ≤ 2147483647:
   ```c
   if (WORDSIZE==64 && p <= 2147483647) {
       // Fast path using p2 = p*p
   } else {
       // Slower path using Montgomery multiplication
   }
   ```

3. **Potential Overflow**: When p > 2^31-1, the fast path is disabled, and the algorithm may have issues with:
   - Modular arithmetic computations
   - Polynomial reduction producing all-zero results
   - Pair selection criteria failing

### Detailed Analysis

The infinite loop occurs in the main algorithm loop in `f4gb_mod()`:

```c
while (f4_pload > 0) {
    // ... algorithm steps ...
    for (i=0; i < f4_aload; i++) {
        b = f4_array[i];
        b->sug = d;
        f4_update(b);  // Should generate new pairs
    }
}
```

When `f4_aload` is 0 (no new polynomials), but `f4_pload` remains 1 (one pair to process), the loop cannot terminate. The pair keeps being selected but produces no new polynomials after reduction.

### Critical Code Sections

1. **Division by Zero** (line in f4_symbol):
   ```c
   printf("%lld x %lld with %lld non-zero, %.1f per row\n", 
       (long long int)f4_aload,
       (long long int)f4_uload,
       (long long int)l,
       (double)l/f4_aload  // Division by zero when f4_aload = 0
   );
   ```

2. **Prime-dependent arithmetic** (f4_reduce_vector and related functions):
   - Different code paths for p ≤ 2147483647 vs larger primes
   - Potential issues with Montgomery multiplication for large primes

## Recommendations

### Immediate Workaround
Restrict prime modulus to values < 2^31 (2147483647) until the bug is fixed.

### Suggested Fixes

1. **Add Zero-Check Guard**:
   ```c
   // In f4_symbol() before the printf
   if (f4_aload == 0) {
       printf("0 x %lld with 0 non-zero, 0.0 per row\n", 
              (long long int)f4_uload);
   } else {
       printf("%lld x %lld with %lld non-zero, %.1f per row\n", 
              (long long int)f4_aload,
              (long long int)f4_uload,
              (long long int)l,
              (double)l/f4_aload);
   }
   ```

2. **Add Loop Termination Check**:
   ```c
   // In f4gb_mod() after f4_update loop
   if (f4_aload == 0 && f4_pload > 0) {
       printf("WARNING: No progress made, %lld pairs remaining\n", 
              (long long int)f4_pload);
       // Either break or handle the degenerate case
       break;
   }
   ```

3. **Review Large Prime Arithmetic**:
   - Verify Montgomery multiplication implementation for primes ≥ 2^31
   - Check for integer overflow in intermediate calculations
   - Consider using INT = int64_t consistently for large prime support

## Minimal Reproducible Example

```c
#include <stdio.h>
#include "axf4_wrapper.h"

int main() {
    const char* vars[] = {"x", "y"};
    
    // This works
    axf4_session_t session1 = axf4_create_session(1073741827, vars, 2);
    axf4_add_polynomial(session1, "x^2+y^2-1");
    axf4_add_polynomial(session1, "x-y");
    axf4_result_t result1 = axf4_compute_groebner_basis(session1);
    printf("Prime 1073741827: %s\n", 
           result1.status == 0 ? "SUCCESS" : "FAILED");
    axf4_free_result(&result1);
    axf4_destroy_session(session1);
    
    // This hangs
    axf4_session_t session2 = axf4_create_session(2147483629, vars, 2);
    axf4_add_polynomial(session2, "x^2+y^2-1");
    axf4_add_polynomial(session2, "x-y");
    // WARNING: This will hang indefinitely!
    axf4_result_t result2 = axf4_compute_groebner_basis(session2);
    
    return 0;
}
```

## Impact
This bug prevents the use of F4 algorithm with prime moduli ≥ 2^31, which may be required for certain applications needing larger coefficient spaces to avoid modular arithmetic artifacts.

## Additional Notes
- The bug is deterministic and reproducible
- No data corruption occurs; the algorithm simply doesn't terminate
- The issue appears to be in the core F4 algorithm, not the wrapper