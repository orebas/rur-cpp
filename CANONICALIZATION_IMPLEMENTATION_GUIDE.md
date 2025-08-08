# Canonicalization Implementation Guide

## Executive Summary

The `canonicalize_rur_permutation` function ensures coefficient vectors have consistent ordering across different prime moduli, enabling successful CRT reconstruction. Without it, position `[i][j]` refers to different monomials for different primes, causing CRT to fail catastrophically.

## Minimal vs Complete Implementation

### Minimal Implementation (INSUFFICIENT)

A minimal implementation might only:
```cpp
// Just pad vectors and make monic
for (size_t i = 0; i < current_table.size(); ++i) {
    current_table[i].resize(reference_table[i].size(), 0);
}
make_monic(current_table[0], prime);
```

**Why This Fails:**

1. **Variable Permutation Ignored**
   - Prime 1: `[minpoly, param_x, param_y, param_z]`
   - Prime 2: `[minpoly, param_z, param_y, param_x]`
   - CRT would combine x's coefficients with z's coefficients → wrong results

2. **Coefficient Reversal Ignored**
   - Prime 1: `[c₀, c₁, c₂, c₃]` (ascending degree)
   - Prime 2: `[c₃, c₂, c₁, c₀]` (descending degree)
   - CRT would combine degree-0 with degree-3 terms → nonsense polynomial

3. **Global Scaling Ignored**
   - Prime 1: Minimal polynomial `T⁸ - 108T⁶ + ...`
   - Prime 2: Minimal polynomial `-3T⁸ + 324T⁶ + ...` (scaled by -3)
   - Different roots, different solutions → wrong system

### Complete Implementation (CORRECT)

Our complete implementation handles all cases:

#### 1. **Detects and Fixes Coefficient Ordering**
```cpp
// Check both normal and reversed ordering
if (are_proportional(cur_minpoly, ref_minpoly, prime, false, scale)) {
    needs_reversal = false;
} else if (are_proportional(cur_minpoly, ref_minpoly, prime, true, scale)) {
    needs_reversal = true;
}
```

#### 2. **Handles Variable Permutations**
```cpp
// Match each parameterization to the correct variable
for (size_t cur_idx = 1; cur_idx < current_table.size(); ++cur_idx) {
    for (size_t ref_idx = 1; ref_idx <= num_variables; ++ref_idx) {
        if (are_proportional(current[cur_idx], reference[ref_idx], ...)) {
            aligned_table[ref_idx] = current[cur_idx];
        }
    }
}
```

#### 3. **Normalizes Global Scaling**
```cpp
// Apply global scale to all polynomials
for (auto& poly : current_table) {
    scale_and_reverse(poly, global_scale, needs_reversal, prime);
}
make_monic(current_table[0], prime);
```

#### 4. **Validates Structural Consistency**
```cpp
// Detect degree mismatches that indicate bad primes
if (current_table[i][j] != 0 && j >= reference_size) {
    throw std::runtime_error("Degree mismatch - bad prime");
}
```

## Key Differences from Julia

### Julia's Approach (Proactive)
- Uses **same separating element** for all primes
- Maintains **fixed monomial ordering** throughout
- No post-hoc canonicalization needed

### C++ Approach (Reactive)
- Computes RUR **independently** for each prime
- Requires **post-hoc alignment** to reference
- More complex but allows parallel computation

## Critical Implementation Details

### 1. **Proportionality Testing**
Two polynomials are proportional if one is a scalar multiple of the other:
- Must handle trailing zeros correctly
- Must check both forward and reverse ordering
- Must compute scale factor using modular inverse

### 2. **Greedy vs Optimal Matching**
- **Greedy** (our approach): Fast, works for non-degenerate systems
- **Optimal** (bipartite matching): Handles all cases but more complex

### 3. **Error Handling Strategy**
```cpp
try {
    canonicalize_rur_permutation(prime_table, reference_table, prime, num_vars);
} catch (const std::exception& e) {
    // Mark prime as bad and continue
    bad_primes_other++;
    continue;
}
```

## Testing Requirements

### Essential Test Cases

1. **Identity Test**
   - Same computation twice should yield identical tables

2. **Permutation Test**
   - Swap variable order, verify correct matching

3. **Reversal Test**
   - Reverse coefficient order, verify detection

4. **Scaling Test**
   - Scale by constant, verify normalization

5. **Bad Prime Test**
   - Structural mismatch should throw exception

## Integration Checklist

- [ ] Replace empty stub in `rur_main_algorithm.cpp:29`
- [ ] Add `normalize_first_table` for first prime
- [ ] Add try-catch around canonicalization
- [ ] Update bad prime counters
- [ ] Add verbose logging for debugging
- [ ] Test with 3-variable system `x²-1, y²-2, z²-3`

## Performance Considerations

- **Time Complexity**: O(n²m) where n = num_variables, m = polynomial degree
- **Space Complexity**: O(nm) for temporary copies
- **Optimization**: Could cache proportionality results if matching many primes

## Conclusion

The complete implementation is **necessary for correctness**. A minimal implementation that only pads vectors will produce wrong results whenever:
- Variables are permuted between primes
- Coefficient ordering differs (ascending vs descending)
- Global scaling differs between computations

The complete implementation detects and fixes all these issues, ensuring CRT reconstruction succeeds.