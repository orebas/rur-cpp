# Modern C++ Multivariate Polynomial System Solver - Implementation Plan

## Executive Summary

We are building a modern C++ library for solving multivariate polynomial systems using:
1. The new F4 parallel Gröbner basis algorithm (axf4.c)
2. Rational Univariate Representation (RUR) for converting to univariate form
3. Root isolation libraries for finding solutions

## Critical Update

**The Julia implementation avoids re-entrancy issues by using a different algorithm!** Instead of computing multiple Gröbner bases, it:
- Builds multiplication tables once
- Uses a bivariate lexicographic conversion algorithm locally
- Avoids additional Gröbner basis calls

This is the approach described in the Demin, Rouillier et al. paper (2402.07141v3.pdf).

## Architecture Overview

```
Input (Multivariate System)
    ↓
F4 Wrapper (single call to axf4.c)
    ↓
Gröbner Basis
    ↓
Multiplication Matrices (MX₁, ..., MXₙ)
    ↓
Trace-based RUR Computation
    ↓
Univariate Representation
    ↓
Root Isolation
    ↓
Solutions
```

## Phase 1: Foundation (2-3 weeks)

### Goal
Create C++ wrapper for axf4.c and implement core data structures.

### Tasks
1. **Data Structures**
   - `Monomial` class using `std::vector<int>` for exponents
   - `MultivariatePolynomial<CoeffT>` template class
   - Use FLINT's `fmpq_t` for rational coefficients
   - Implement monomial orderings (grevlex, lex)

2. **F4 C Wrapper**
   - Create `axf4_wrapper.c/h` exposing clean C API
   - Handle all memory management and thread safety
   - Input: string representation of polynomials
   - Output: string representation of Gröbner basis

3. **F4 C++ Wrapper**
   - `F4Solver` class wrapping the C interface
   - Convert between C++ polynomials and string format
   - Parse results back to C++ data structures

### Deliverable
Ability to compute Gröbner basis from C++ for simple systems.

## Phase 2: RUR via Bivariate Lex Algorithm (3-4 weeks)

### Goal
Implement RUR computation using the Julia/paper algorithm (no additional GB calls).

### Algorithm (Based on Demin, Rouillier et al. 2024)
1. Compute quotient ring basis B = {b₁, ..., bD}
2. Build multiplication tables/matrices for each variable
3. For separating element t = Σ cᵢXᵢ:
   - Compute minimal polynomial of t using FGLM-like iteration
   - For each variable Xᵢ:
     - Use bivariate lexicographic algorithm (biv_lex!)
     - Extract parameterization from bivariate representation
     - No new Gröbner basis computation needed!

### Tasks
1. **Quotient Ring Operations**
   - Implement `NormalForm(poly, basis)` function
   - Compute quotient basis from Gröbner basis

2. **Multiplication Tables**
   - Build multiplication tables (sparse representation)
   - Store as vectors of coefficient vectors
   - Implement table lookup and multiplication

3. **FGLM-like Iteration**
   - Compute minimal polynomial of separating element
   - Use Gaussian elimination on expanding matrix
   - Track free/dependent monomials

4. **Bivariate Lex Algorithm**
   - Port `biv_lex!` function from Julia
   - Build bivariate representation incrementally
   - Extract coefficients for parameterization

5. **Parameterization Extraction**
   - Convert bivariate lex basis to RUR format
   - Handle special cases (shape position, etc.)

### Risk Mitigation
- Study the Julia implementation in detail for algorithm correctness
- The bivariate approach avoids re-entrancy issues completely
- Validate against Julia reference implementation on same test cases

## Phase 3: Integration & Optimization (2-3 weeks)

### Tasks
1. **Root Finding**
   - Integrate MPSolve or similar library
   - Find roots of h(t)
   - Back-substitute to get full solutions

2. **API Design**
   ```cpp
   class SystemSolver {
   public:
       SystemSolver(const std::vector<MultivariatePolynomial>& system);
       
       void compute_groebner_basis();
       int get_dimension(); // Check if zero-dimensional
       void compute_rur();
       void isolate_roots();
       
       const std::vector<Solution>& get_solutions() const;
   };
   ```

3. **Command-Line Interface**
   - Input: polynomial system file
   - Output: solutions in standard format

## Phase 4: Performance & Hardening (2-3 weeks)

### Tasks
1. **Modular Arithmetic**
   - Template specialization for `nmod_t` coefficients
   - Compute modulo multiple primes
   - Chinese Remainder Theorem reconstruction

2. **Optimization**
   - Profile and optimize bottlenecks
   - Implement fast Hankel solver if needed
   - Parallel trace computations

3. **Testing & Documentation**
   - Comprehensive test suite
   - Benchmark against existing tools
   - User documentation

## Dependencies

- **FLINT**: For number theory and polynomial arithmetic
- **Eigen**: For linear algebra operations
- **MPSolve**: For univariate root finding
- **C++ Standard**: C++17 or later

## Risk Summary

1. **High Risk**: Coefficient growth in rational arithmetic
   - Mitigation: Modular computation strategy

2. **Medium Risk**: Performance of trace computations
   - Mitigation: Start simple, optimize based on profiling

3. **Low Risk**: Integration complexity
   - Mitigation: Clean interfaces between components

## Success Criteria

1. Correctly solve benchmark polynomial systems
2. Performance competitive with existing tools (msolve)
3. Clean, maintainable C++ codebase
4. Both library and standalone usage supported