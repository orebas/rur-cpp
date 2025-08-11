# Revised RUR Implementation Plan: Mirror Julia Structure

## 🎯 **Overall Strategy**

**Primary Goal**: Mirror the Julia RationalUnivariateRepresentation.jl implementation as closely as possible while adapting to C++ idioms. Focus on correctness first, optimization last.

**Key Principles**:
1. **Mirror Julia Exactly**: Keep variable names, algorithm structure, and data layouts identical where possible
2. **Incremental Testing**: Unit test each component individually (better than Julia package testing)
3. **Defer Optimization**: Make notes for later but implement straightforward versions first
4. **Property Verification**: Test mathematical properties, not just integration

---

## 📋 **Detailed Phase Breakdown**

### **Phase 1: Core Data Structures (Mirror Julia exactly)**
**Goal**: Implement the exact data structures from Julia with proper C++ types

#### Task 1.1: StackVect Structure
```cpp
// Mirror Julia's StackVect exactly
struct StackVect {
    int32_t pos;     // Position index  
    PP mon;          // Monomial (power product)
    int32_t prev;    // Previous element in chain
    int32_t var;     // Variable index
};
```
**Test**: Create StackVect instances and verify field access

#### Task 1.2: Type Definitions (Match Julia)
```cpp
using PP = std::vector<uint32_t>;           // Power product (Julia: Vector{Deg})
using ModularCoeff = uint32_t;              // Modular coefficient (Julia: UInt32)
using AccModularCoeff = uint64_t;           // Accumulator (Julia: UInt64)
using Deg = uint32_t;                       // Degree type (Julia: UInt32)
```
**Test**: Verify size requirements and arithmetic operations

#### Task 1.3: Core Arrays (Match Julia variable names)
```cpp
class JuliaStyleQuotientRing {
    std::vector<PP> quo;                           // Quotient basis monomials (Julia: quo)
    std::vector<StackVect> t_xw;                  // Border structure (Julia: t_xw)  
    std::vector<std::vector<int32_t>> i_xw;       // Variable multiplication indices (Julia: i_xw)
    std::vector<std::vector<ModularCoeff>> t_v;   // Coefficient vectors (Julia: t_v)
};
```
**Test**: Memory layout, resize operations, index access patterns

---

### **Phase 2: Multiplication Tables (Julia algorithms)**
**Goal**: Implement the exact Julia table construction process

#### Task 2.1: `prepare_table_mxi` Implementation
**Julia Reference**: Lines 341-392 in RationalUnivariateRepresentation.jl
```cpp
void prepare_table_mxi(
    const std::vector<PP>& quo,                    // Quotient basis
    const std::vector<PP>& ltg,                    // GB leading terms  
    std::vector<StackVect>& t_xw,                  // Border structure (output)
    std::vector<std::vector<int32_t>>& i_xw        // Variable indices (output)
);
```
**Algorithm**:
- For each basis monomial m and variable xi:
  - Compute xi·m  
  - Classify as: quotient element, GB leading term, or new border element
- Build fast index arrays i_xw[i][j] = index of xi·mj

**Test**: 
- Verify border classification for x²-1, y²-1 system
- Check i_xw indices point to correct elements
- Validate all quotient elements are reachable

#### Task 2.2: `learn_compute_table!` Implementation  
**Julia Reference**: Lines 439-463
```cpp
void learn_compute_table(
    std::vector<std::vector<ModularCoeff>>& t_v,   // Coefficient vectors (output)
    const std::vector<StackVect>& t_xw,            // Border structure
    const std::vector<std::vector<int32_t>>& i_xw, // Variable indices
    const ModularArithmetic& arithm                // Modular arithmetic context
);
```
**Algorithm**:
- Iteratively compute coefficient vectors for border elements
- Use chain rule: if we know expansion of m, we can compute xi·m
- Fill t_v with proper modular reductions

**Test**:
- Verify coefficient vectors sum correctly  
- Check modular arithmetic consistency
- Test against known small examples (manual calculation)

#### Task 2.3: `_mul_var_quo!` Implementation
**Julia Reference**: Lines 398-431  
```cpp
void mul_var_quo(
    std::vector<ModularCoeff>& result,             // Output vector
    const std::vector<ModularCoeff>& input,       // Input vector in quotient basis
    int32_t var_index,                             // Variable to multiply by
    const std::vector<std::vector<int32_t>>& i_xw, // Pre-computed indices
    const std::vector<std::vector<ModularCoeff>>& t_v, // Coefficient vectors
    const ModularArithmetic& arithm                // Modular arithmetic
);
```
**Algorithm**:
- Use pre-computed indices for O(1) lookup
- Apply coefficient vectors from t_v
- Perform packed modular reduction

**Test**:
- Compare with direct polynomial multiplication
- Verify distributivity: (a+b)*x = a*x + b*x  
- Check associativity: (a*x)*y = a*(x*y)
- Performance test: should be significantly faster than direct method

---

### **Phase 3: Bivariate Algorithm (Critical missing piece)**
**Goal**: Implement the exact Julia bivariate lexicographic approach

#### Task 3.1: `first_variable` Refinement
**Current Status**: ✅ Implemented and working
**Improvement**: Match Julia coefficient handling exactly
**Test**: Verify minimal polynomial matches Julia output exactly

#### Task 3.2: `biv_lex!` Implementation
**Julia Reference**: Lines 635-684
```cpp
struct BivariateResult {
    std::vector<std::tuple<int32_t, int32_t>> monomial_basis;    // (deg_in_T, deg_in_xi)
    std::vector<std::tuple<int32_t, int32_t>> leading_monomials; // Leading terms
    std::vector<std::vector<int32_t>> generators;                // Relation generators
};

BivariateResult biv_lex(
    const std::vector<std::vector<ModularCoeff>>& t_v,    // Coefficient vectors
    const std::vector<std::vector<int32_t>>& i_xw,        // Variable indices  
    const std::vector<ModularCoeff>& initial_basis,       // From first_variable
    int32_t variable_index,                                // Variable to extend with
    const ModularArithmetic& arithm                       // Modular arithmetic
);
```
**Algorithm**:
1. Maintain bivariate monomials as pairs (deg_in_T, deg_in_xi)
2. Start with univariate basis from first_variable  
3. Multiply by xi and reduce using quotient operations
4. Track new relations as they appear
5. Output: bivariate basis and generating relations

**Test**:
- Verify bivariate monomial generation is systematic
- Check that relations are linearly independent
- Test with x²-1, y²-1 system: should get correct x,y parameterizations

#### Task 3.3: `convert_biv_lex_2_biv_pol` Implementation
**Julia Reference**: Lines 686-720
```cpp
std::vector<std::vector<std::vector<ModularCoeff>>> convert_biv_lex_2_biv_pol(
    const std::vector<std::vector<int32_t>>& generators,           // From biv_lex
    const std::vector<std::tuple<int32_t, int32_t>>& monomial_basis, // Bivariate basis
    const std::vector<std::tuple<int32_t, int32_t>>& leading_terms   // Leading terms  
);
```
**Algorithm**:
- Convert bivariate relations to polynomial matrices
- Each matrix represents lt_g[i] - sum(n_g[i][j] * m_b[j])
- Organize by degree in second variable

**Test**:
- Verify polynomial reconstruction from bivariate data
- Check that converted polynomials evaluate correctly
- Test rational function extraction works

---

### **Phase 4: Integration and End-to-End Testing**
**Goal**: Connect all pieces and verify complete RUR computation

#### Task 4.1: Complete RUR Pipeline Integration
```cpp
class JuliaStyleRURSolver {
    // Phase 1: Data structures
    JuliaStyleQuotientRing quotient_ring_;
    
    // Phase 2: Multiplication tables  
    void build_multiplication_tables();
    
    // Phase 3: Bivariate computation
    void compute_variable_parameterizations();
    
    // Complete pipeline
    void solve_single_prime();
};
```

#### Task 4.2: Property-Based Testing
**Tests to implement**:
- **Correctness**: RUR solutions satisfy original polynomial system
- **Completeness**: Number of solutions matches expected (Bézout bound, etc.)
- **Consistency**: Same results across different prime choices
- **Regression**: Match Julia results on benchmark systems

#### Task 4.3: Comprehensive Unit Testing  
**Test each component**:
- Data structure operations (creation, access, modification)
- Table construction (border classification, index building, coefficient computation)
- Bivariate algorithm (monomial generation, relation detection, polynomial conversion)
- Integration (pipeline execution, error handling, edge cases)

---

### **Phase 5: Multi-Modular Computation**
**Goal**: Lift from single prime to rational numbers

#### Task 5.1: Multi-Modular Framework
```cpp
class MultiModularRUR {
    std::vector<uint32_t> primes_;                    // Prime list
    std::vector<JuliaStyleRURSolver> solvers_;        // One per prime
    
    void compute_multi_modular();                     // Parallel execution
    void combine_results();                           // CRT combination
};
```

#### Task 5.2: CRT and Rational Reconstruction
**Julia Reference**: Lines 1169-1246 (`crt_and_ratrec!`)
```cpp
class RationalReconstructor {
    void chinese_remainder_theorem(const std::vector<ModularCoeff>& residues,
                                  const std::vector<uint32_t>& primes);
    void rational_reconstruction(const BigInteger& crt_value,
                                const BigInteger& modulus);
};
```

#### Task 5.3: F4 Learning/Replay Mechanism  
**Goal**: Avoid F4 re-entrancy issues
```cpp
class F4LearnReplay {
    void learn_phase(const PolynomialSystem& system);    // Store reduction graph
    void apply_phase(uint32_t prime);                    // Replay with new prime
};
```

---

### **Phase 6: Production Hardening**
**Goal**: Make the implementation robust and reliable

#### Task 6.1: Error Handling and Validation
- Input validation (polynomial system requirements)
- Numerical stability checks (prime choice, condition numbers)
- Graceful fallback strategies (separating element failure, etc.)

#### Task 6.2: Memory Management and RAII
- Ensure all FLINT resources are properly managed
- Memory leak detection and prevention
- Exception safety guarantees

#### Task 6.3: Comprehensive Documentation
- API documentation with examples
- Algorithm explanation with references
- Performance characteristics and limitations

---

### **Phase 7: Optimization (LAST)**
**Goal**: Performance improvements after correctness is established

#### Task 7.1: Performance Profiling
- Identify bottlenecks in table construction
- Memory access pattern analysis
- Arithmetic operation optimization opportunities

#### Task 7.2: Algorithmic Optimizations
```cpp
// TODO OPTIMIZATION: Use packed arithmetic like Julia's pack_value()
// TODO OPTIMIZATION: SIMD for modular operations  
// TODO OPTIMIZATION: Cache-friendly memory layout
// TODO OPTIMIZATION: Vectorized coefficient operations
```

#### Task 7.3: Parallelization
```cpp
// TODO OPTIMIZATION: Parallel prime processing
// TODO OPTIMIZATION: Thread-local contexts and reduction graphs
// TODO OPTIMIZATION: Lock-free data structures for shared results
```

---

## 🧪 **Testing Strategy**

### **Unit Testing Philosophy**
**Better than Julia package**: Test each component individually rather than just integration testing

1. **Data Structure Tests**
   - Creation, modification, access patterns
   - Memory layout and alignment
   - Iterator behavior and bounds checking

2. **Algorithm Tests**  
   - Mathematical properties (commutativity, associativity, etc.)
   - Edge cases (empty inputs, single elements, maximum sizes)
   - Error conditions and recovery

3. **Integration Tests**
   - End-to-end pipeline execution
   - Cross-component data flow
   - Performance benchmarks

### **Property-Based Testing**
Use mathematical properties to verify correctness:

```cpp
// Property: RUR solutions satisfy original system
PROPERTY_TEST(rur_solutions_satisfy_system) {
    auto system = generate_polynomial_system();
    auto rur = solve_rur(system);
    auto solutions = evaluate_rur(rur);
    
    for (const auto& solution : solutions) {
        REQUIRE(system.evaluate(solution) == zero_vector);
    }
}

// Property: Solution count matches expected
PROPERTY_TEST(solution_count_correct) {
    auto system = generate_zero_dimensional_system();
    auto rur = solve_rur(system);
    auto expected_count = compute_bezout_bound(system);
    
    REQUIRE(rur.solution_count() <= expected_count);
}
```

### **Regression Testing**
Test against known polynomial systems with verified results:

```cpp
namespace RegressionTests {
    void test_cyclic_systems();      // Cyclic-n polynomial systems
    void test_katsura_systems();     // Katsura-n systems  
    void test_eco_systems();         // Eco-n systems
    void test_random_systems();      // Randomly generated systems
}
```

---

## 🔍 **Code Review Process**

**CRITICAL: After each function implementation, we MUST:**

1. **Gemini Pro Review**: Use `mcp__zen__codereview` with Gemini Pro model
2. **Provide Julia Reference**: Include the corresponding Julia code for comparison
3. **Cross-validation**: Ensure C++ implementation faithfully mirrors Julia algorithm
4. **Document Findings**: Note any discrepancies or improvements suggested

**Review Template:**
```
Files to review: [C++ implementation files]
Focus on: Faithfulness to Julia implementation, correctness, potential bugs
Include Julia reference code for comparison
Model: gemini-2.5-pro-preview-06-05
```

**Review Checklist:**
- [ ] Algorithm logic matches Julia exactly
- [ ] Variable names and data structures align
- [ ] Edge cases handled correctly  
- [ ] No obvious bugs or memory issues
- [ ] Performance considerations noted
- [ ] Areas for future optimization identified

## 📝 **Implementation Guidelines**

### **C++ Adaptations from Julia**
1. **Memory Management**: Use RAII instead of manual allocation/deallocation
2. **Error Handling**: Use exceptions instead of error return codes
3. **Type Safety**: Use strong typing instead of dynamic typing where beneficial
4. **Const Correctness**: Mark immutable data structures and methods as const
5. **Iterator Support**: Provide standard C++ iterators for container-like structures

### **Julia Compatibility**
1. **Variable Names**: Keep Julia names (`t_xw`, `i_xw`, `t_v`, etc.) for easy cross-reference
2. **Algorithm Structure**: Follow Julia control flow and iteration patterns
3. **Data Layouts**: Use similar memory layouts (row-major adaptations where needed)
4. **Modular Arithmetic**: Match Julia's modular reduction patterns exactly

### **Code Organization**
```
src/julia_rur/
├── data_structures.hpp/cpp     # StackVect, core arrays
├── multiplication_tables.hpp/cpp # prepare_table_mxi, learn_compute_table
├── bivariate_algorithm.hpp/cpp   # biv_lex, convert_biv_lex_2_biv_pol  
├── multi_modular.hpp/cpp         # CRT, rational reconstruction
├── julia_rur_solver.hpp/cpp      # Main integration class
└── tests/
    ├── test_data_structures.cpp
    ├── test_multiplication_tables.cpp
    ├── test_bivariate_algorithm.cpp
    ├── test_integration.cpp
    └── test_regression.cpp
```

---

## 🎯 **Success Criteria**

### **Phase 1-3 Success**: 
- Single prime RUR computation works correctly
- Results match Julia output for test systems
- All unit tests pass

### **Phase 4-5 Success**:
- Multi-modular computation produces rational results  
- Performance is reasonable (within 10x of Julia initially)
- Comprehensive test suite passes

### **Phase 6-7 Success**:
- Production-ready implementation  
- Performance competitive with Julia
- Full documentation and examples

---

## 🚀 **Immediate Next Steps**

1. **Start with Task 1.1**: Implement StackVect structure and basic types
2. **Create test framework**: Set up unit testing infrastructure  
3. **Implement incrementally**: Each task should be fully tested before moving to next
4. **Document as we go**: Maintain clear mapping between C++ and Julia code

This plan provides a clear roadmap for implementing a Julia-faithful RUR solver in C++ with better testing and proper optimization deferral.