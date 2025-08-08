# Julia RUR Algorithm Documentation

This document provides a comprehensive analysis of the RationalUnivariateRepresentation.jl implementation, focusing on the algorithmic aspects that need to be ported to C++.

## Table of Contents
1. [Overview](#overview)
2. [Key Algorithmic Insights](#key-algorithmic-insights)
3. [Data Structures](#data-structures)
4. [Algorithm Flow](#algorithm-flow)
5. [Multiplication Tables](#multiplication-tables)
6. [Bivariate Lexicographic Algorithm](#bivariate-lexicographic-algorithm)
7. [Modular Arithmetic and Multi-Modular Computation](#modular-arithmetic-and-multi-modular-computation)
8. [Rational Reconstruction](#rational-reconstruction)
9. [Implementation Strategy for C++](#implementation-strategy-for-c)

## Overview

The Julia RUR implementation computes a Rational Univariate Representation of a zero-dimensional polynomial system. The key innovation is that it **avoids re-entrancy issues** by:

1. Computing the Gröbner basis **only once** using a learning phase
2. Building multiplication tables directly from the GB without additional GB calls
3. Using a bivariate lexicographic algorithm to compute the parameterization
4. Lifting results from modular arithmetic to Q using CRT and rational reconstruction

## Key Algorithmic Insights

### 1. Avoiding Re-entrancy
- The algorithm computes GB once with `_gb_4_rur_learn` and captures a "reduction graph"
- For subsequent primes, it uses `_gb_4_rur_apply!` to replay the graph
- All further computations use only linear algebra on pre-computed tables

### 2. Bivariate Approach
- Instead of computing trace formulas (which would require multiple GB calls)
- Uses `first_variable` to compute univariate minimal polynomial
- Uses `biv_lex!` to compute bivariate relations for parameterization
- This completely avoids needing to call GB solver multiple times

## Data Structures

### Core Types
```julia
PP = Vector{Deg}                    # Power product (monomial exponent vector)
Deg = UInt32                        # Degree type
ModularCoeff = UInt32               # Modular coefficient
AccModularCoeff = UInt64            # Accumulator for modular arithmetic
```

### Key Structures

#### StackVect (Lines 57-62)
```julia
mutable struct StackVect
    pos::Int32      # Position index
    mon::PP         # Monomial (power product)
    prev::Int32     # Previous element in chain
    var::Int32      # Variable index
end
```
Used to explore the border of the quotient ring during multiplication table construction.

### Main Data Arrays

1. **ltg**: Leading terms of Gröbner basis
2. **quo**: Quotient basis monomials (sorted by DRL order)
3. **t_xw**: Vector of StackVect - the border structure
4. **i_xw**: For each variable xi, maps basis elements to their products with xi
5. **t_v**: Coefficient vectors - expansion of each border element in quotient basis

## Algorithm Flow

### Main Pipeline (zdim_parameterization → _zdim_multi_modular_RUR!)

1. **Choose first prime and compute GB once**
   ```
   _zdim_modular_RUR_LV (or variant) 
   ├─ _gb_4_rur_learn: Compute GB and learn reduction graph
   ├─ compute_quotient_basis: Find monomial basis of quotient ring
   ├─ prepare_table_mxi: Build border structure and indices
   ├─ compute_fill_quo_gb!: Initialize coefficient vectors
   ├─ learn_compute_table!: Fill multiplication tables
   └─ _zdim_parameterization: Compute RUR using bivariate algorithm
   ```

2. **Multi-modular loop**
   ```
   For each block of primes:
   ├─ _gb_4_rur_apply!: Replay GB computation
   ├─ _zdim_modular_RUR_LV_apply!: Replay table construction
   └─ crt_and_ratrec!: Lift to rational numbers
   ```

## Multiplication Tables

### Construction Process

1. **compute_quotient_basis** (Lines 181-214)
   - BFS traversal of monomial lattice
   - Stop when monomial is divisible by GB leading term
   - Returns sorted list of quotient basis monomials

2. **prepare_table_mxi** (Lines 341-392)
   - For each basis monomial m and variable xi:
     - Compute xi·m
     - Classify as: quotient element, GB leading term, or new border element
   - Build fast index arrays i_xw[i][j] = index of xi·mj

3. **learn_compute_table!** (Lines 439-463)
   - Iteratively compute coefficient vectors for border elements
   - Uses chain rule: if we know expansion of m, we can compute xi·m

### Multiplication Operation (_mul_var_quo!, Lines 398-431)
- Computes xi·v where v is a vector in quotient basis
- Uses pre-computed indices for O(1) lookup
- Aggressive loop unrolling and packed modular reduction

## Bivariate Lexicographic Algorithm

This is the heart of avoiding re-entrancy!

### first_variable (Lines 544-633)
Computes minimal polynomial of separating element T = xn:

1. Start with v = xn in quotient basis
2. Iteratively compute T^k = xn·T^(k-1)
3. Use Gaussian reduction to detect linear dependence
4. When T^d lies in span of {1, T, ..., T^(d-1)}, extract minimal polynomial

### biv_lex! (Lines 635-684)
Extends to two variables (T and xi):

1. Maintain bivariate monomials as pairs (deg_in_T, deg_in_xi)
2. Start with univariate basis from first_variable
3. Multiply by xi and reduce
4. Track new relations as they appear
5. Output: bivariate basis and generating relations

### Parameterization Formula
For each variable xi:
- Run biv_lex! to get bivariate relations
- Extract rational function xi = fi(T)/f'(T)
- Where f'(T) is derivative of square-free part of minimal polynomial

## Modular Arithmetic and Multi-Modular Computation

### Arithmetic Abstractions
- **ModularArithZp**: Single prime arithmetic
- **ModularArithMZp**: Composite arithmetic (multiple primes packed)

### Multi-Modular Strategy

1. **Learning Phase** (first prime)
   - Full computation including GB
   - Store reduction graph and all intermediate data

2. **Application Phase** (subsequent primes)
   - Replay GB using stored graph
   - Replay table construction
   - No Buchberger algorithm re-entry!

3. **Parallelization**
   - Blocks of primes processed in parallel
   - Each thread has private copy of reduction graph
   - Composite arithmetic for SIMD efficiency

## Rational Reconstruction

### CRT and Rational Recovery (crt_and_ratrec!, Lines 1169-1246)

For each coefficient position:
- Maintain: mod p images, CRT integer, rational result
- After each prime block:
  1. Update CRT using fast multi-modular CRT
  2. Try multiple rational reconstruction strategies:
     - Unbalanced with known denominator
     - Very unbalanced (CRT-like)
     - Standard balanced
     - Using previous coefficient's denominator
  3. Track common denominator across coefficients

### Key Optimization
Process coefficients right-to-left so later coefficients provide denominator hints for earlier ones.

## Implementation Strategy for C++

### Phase 1: Core Data Structures
```cpp
using PP = std::vector<uint32_t>;       // Power product
using Coeff = uint32_t;                 // Modular coefficient
using AccumCoeff = uint64_t;            // Accumulator

struct StackVect {
    int32_t pos;
    PP mon;
    int32_t prev;
    int32_t var;
};

struct QuotientRing {
    std::vector<PP> basis;              // Monomial basis
    std::vector<StackVect> border;      // t_xw
    std::vector<std::vector<int32_t>> var_indices;  // i_xw
    std::vector<std::vector<Coeff>> coeff_vectors;  // t_v
};
```

### Phase 2: Single Prime Implementation
1. Integrate with F4 wrapper to get GB
2. Implement compute_quotient_basis
3. Build multiplication tables
4. Implement first_variable algorithm
5. Implement biv_lex!
6. Extract parameterization

### Phase 3: Multi-Modular Lift
1. Implement modular arithmetic traits
2. Add CRT accumulation
3. Implement rational reconstruction
4. Add parallel prime processing

### Critical Implementation Notes

1. **Memory Layout**: The Julia code uses column-major arrays. C++ should use row-major for cache efficiency.

2. **Modular Reduction**: Pack reduction operations (see pack_value function) to minimize overhead.

3. **Graph Storage**: The GB reduction graph needs careful serialization to enable replay.

4. **FLINT Integration**: Use FLINT's fmpq_t for rational arithmetic and nmod_t for modular.

5. **Parallelization**: Use thread-local storage for reduction graphs in parallel mode.

## Separating Element Strategies

The implementation provides multiple strategies for finding a separating linear form:

1. **current** (default): Try last variable, then swap variables, then add linear forms
2. **random**: Random coefficients in bounded range
3. **l0_norm**: Prefer sparse linear forms
4. **deterministic**: Systematic search using power sequences

When the last variable is not separating, the system is extended with a new variable T and polynomial T - (a1·x1 + ... + an·xn).

## Summary

The Julia RUR implementation is sophisticated and well-optimized. The key insight for C++ porting is that the bivariate lexicographic algorithm completely avoids the need for multiple GB computations, making it compatible with the non-reentrant F4 implementation. The multi-modular lift with rational reconstruction provides the path from Fp to Q.