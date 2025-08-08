# CRT Reconstruction Failure - Root Cause Analysis

## Executive Summary

The CRT rational reconstruction failure for the system `x²-1, y²-2, z²-3` is caused by an **empty canonicalization function** that fails to ensure consistent coefficient vector ordering across different prime moduli. This causes position `[0][6]` to refer to different monomial terms for different primes, leading to CRT attempting to reconstruct coefficients from incompatible values.

## The Problem

### System Under Test
- Polynomials: `x² - 1`, `y² - 2`, `z² - 3`
- Expected: 8 solutions (2³ combinations)
- Failure: CRT reconstruction fails at coefficient `[0][6]` of the minimal polynomial

### Symptom
```
CRT DEBUG: Coefficient [0][6] remainders (first 10):
  Prime 2147483629: remainder 2147480101  
  Prime 2147483587: remainder 2147479867
  Prime 2147483579: remainder 2147480183
  ...
```
These remainders vary by thousands, causing CRT to produce a 2000+ digit result instead of a small integer.

## Root Cause Identified

### The Empty Function
In `/home/orebas/code/rur-cpp/src/julia_rur/rur_main_algorithm.cpp:29`:

```cpp
static void canonicalize_rur_permutation(
    std::vector<std::vector<ModularCoeff>> &current_table,
    const std::vector<std::vector<ModularCoeff>> &reference_table,
    ModularCoeff prime,
    size_t num_variables) {}  // <-- COMPLETELY EMPTY!
```

### What This Function Should Do

Based on analysis with O3 and Gemini, and comparison with Julia's RUR implementation:

1. **Ensure Consistent Vector Lengths**
   - Pad coefficient vectors to the same length across all primes
   - Handle trailing zeros consistently

2. **Maintain Fixed Monomial Ordering**
   - Store coefficients in ascending degree order (T⁰, T¹, T², ...)
   - Ensure position `[i][j]` refers to the same monomial term across all primes

3. **Normalize Scalars**
   - Make minimal polynomial monic (leading coefficient = 1)
   - Scale all parameterizations consistently

## How Julia Handles This

Julia's RationalUnivariateRepresentation.jl maintains consistency through:

1. **Deterministic Separating Element Selection**
   - Uses strategies like `:current`, `:l0_norm`, `:deterministic`
   - Same separating element is found for all primes

2. **Fixed Monomial Basis Ordering**
   - Quotient basis computed with consistent ordering
   - For our system: always `{1, x, y, z, xy, xz, yz, xyz}`

3. **Automatic Alignment**
   - Julia's `crt_and_ratrec!` function assumes `tables_zp[k][i][j]` refers to the same monomial across all primes `k`

## The Fix

### Minimal Implementation

```cpp
static void canonicalize_rur_permutation(
    std::vector<std::vector<ModularCoeff>> &current_table,
    const std::vector<std::vector<ModularCoeff>> &reference_table,
    ModularCoeff prime,
    size_t num_variables) {
    
    // 1. Ensure consistent vector lengths
    for (size_t i = 0; i < current_table.size(); ++i) {
        if (i < reference_table.size()) {
            size_t ref_size = reference_table[i].size();
            current_table[i].resize(ref_size, 0);
        }
    }
    
    // 2. Make minimal polynomial monic
    if (!current_table[0].empty()) {
        ModularCoeff lead = current_table[0].back();
        if (lead != 0 && lead != 1) {
            ModularCoeff inv = mod_inverse(lead, prime);
            for (auto& poly : current_table) {
                for (auto& coeff : poly) {
                    coeff = (static_cast<AccModularCoeff>(coeff) * inv) % prime;
                }
            }
        }
    }
}
```

## Evidence Trail

### O3's Analysis
"The wildly-different remainders shown for coefficient [0][6] do not come from CRT itself – they are already different before CRT is applied. The same term of the minimal polynomial is being read from different positions in the per-prime tables."

### Gemini's Insight
"Without canonicalization, position [0][6] in the minimal polynomial refers to different monomial terms across different primes... The function name `canonicalize_rur_permutation` suggests its purpose is to re-order the computed coefficients to match the canonical ordering."

### Debug Output Confirmation
Our debug programs showed:
- Coefficient vectors have different lengths across primes
- Position `[0][6]` contains wildly different values
- These values appear to be small negative numbers in modular arithmetic

## Impact

Without proper canonicalization:
1. CRT attempts to combine coefficients from different monomials
2. The reconstruction produces massive integers (2000+ digits)
3. Rational reconstruction fails completely
4. The entire RUR computation fails

## Next Steps

1. **Implement the canonicalization function** with at least:
   - Vector length alignment
   - Monic normalization
   
2. **Consider more robust fixes**:
   - Ensure consistent monomial basis ordering throughout
   - Add validation checks for coefficient consistency
   - Implement permutation detection and correction

3. **Test thoroughly**:
   - Verify CRT reconstruction succeeds for the 3-variable system
   - Test with other polynomial systems
   - Compare results with Julia's implementation

## Lessons Learned

1. **Consistency is Critical**: In modular algorithms, maintaining consistent data representation across different moduli is essential.

2. **Empty Functions are Dangerous**: An unimplemented canonicalization function can cause subtle, hard-to-debug failures.

3. **Reference Implementations Help**: Comparing with Julia's working implementation was invaluable for understanding the correct behavior.

4. **AI Tools are Powerful**: O3 and Gemini quickly identified the root cause by analyzing the symptoms and code structure.