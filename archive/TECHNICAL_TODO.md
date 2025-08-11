# Technical TODO for RUR Implementation

## Immediate Next Steps (High Priority)

### 1. Quotient Ring from Gröbner Basis
```cpp
// File: src/julia_rur/quotient_ring_gb.hpp/cpp
namespace julia_rur {

// Extract quotient basis from Gröbner basis
std::vector<PP> extract_quotient_basis(
    const std::vector<PP>& gb_leading_terms,
    int32_t num_variables
);

// Fill quotient ring structure from GB
void compute_fill_quo_gb(
    std::vector<StackVect>& t_xw,          // Output: border structure
    const std::vector<PP>& gb_exponents,   // Gröbner basis exponents
    const std::vector<std::vector<ModularCoeff>>& gb_coeffs,  // GB coefficients
    const std::vector<PP>& quotient_basis, // Quotient basis
    ModularCoeff prime
);

}
```

### 2. F4 Learning/Replay Integration
```cpp
// File: src/julia_rur/f4_integration.hpp/cpp
namespace julia_rur {

struct F4Graph {
    // Opaque structure from F4 learning
    void* internal_graph;
};

// Learn F4 computation graph
F4Graph gb_learn(
    const std::vector<std::vector<PP>>& exponents,
    const std::vector<std::vector<ModularCoeff>>& coefficients,
    ModularCoeff prime
);

// Apply learned graph to new prime
std::pair<std::vector<PP>, std::vector<std::vector<ModularCoeff>>>
gb_apply(
    const F4Graph& graph,
    const std::vector<std::vector<PP>>& exponents,
    const std::vector<std::vector<ModularCoeff>>& coefficients,
    ModularCoeff prime
);

}
```

### 3. Univariate Parameterization
```cpp
// File: src/julia_rur/univariate_parameterization.hpp/cpp
namespace julia_rur {

struct UnivariateParmResult {
    bool success;
    std::vector<std::vector<ModularCoeff>> rur_polynomials;
    int32_t degree;
};

// Main univariate parameterization
UnivariateParmResult zdim_parameterization(
    const std::vector<std::vector<ModularCoeff>>& multiplication_tables,
    const std::vector<std::vector<int32_t>>& variable_indices,
    int32_t expected_degree,
    int32_t start_column,
    ModularCoeff prime
);

}
```

### 4. Separating Element Selection
```cpp
// File: src/julia_rur/separating_element.hpp/cpp
namespace julia_rur {

struct SeparatingElementResult {
    bool found;
    std::vector<int32_t> linear_form;  // Coefficients of separating element
    std::vector<std::vector<ModularCoeff>> parameterization;
};

// Different strategies
SeparatingElementResult find_separating_element_current(
    const std::vector<std::vector<PP>>& gb_exponents,
    const std::vector<std::vector<ModularCoeff>>& gb_coeffs,
    const std::vector<PP>& quotient_basis,
    ModularCoeff prime
);

SeparatingElementResult find_separating_element_random(
    // ... parameters ...
    uint32_t random_seed = 42
);

}
```

### 5. Main RUR Algorithm
```cpp
// File: src/julia_rur/rur_algorithm.hpp/cpp
namespace julia_rur {

struct ModularRURResult {
    bool success;
    std::vector<int32_t> separating_form;
    std::vector<std::vector<ModularCoeff>> parameterization;
    std::vector<PP> leading_terms;
    std::vector<PP> quotient_basis;
    // ... other fields ...
};

// Main modular RUR computation
ModularRURResult zdim_modular_RUR(
    const std::vector<std::vector<PP>>& input_exponents,
    const std::vector<std::vector<ModularCoeff>>& input_coeffs,
    ModularCoeff prime,
    const std::string& strategy = "current",
    bool learn = false
);

}
```

### 6. Multi-modular Main Loop
```cpp
// File: src/julia_rur/multi_modular_rur.hpp/cpp
namespace julia_rur {

struct RURResult {
    std::vector<std::vector<mpq_class>> rur_polynomials;
    std::vector<int32_t> separating_element;
    bool is_radical;
};

// Main multi-modular RUR
RURResult zdim_multi_modular_RUR(
    const std::vector<std::vector<PP>>& input_exponents,
    const std::vector<std::vector<mpq_class>>& input_coeffs,
    int32_t prime_bits = 28,
    const std::string& strategy = "current",
    int32_t num_threads = 1
);

}
```

## Integration Points

### Connect F4 to Multiplication Tables
1. Parse F4 output format (currently using axf4 C wrapper)
2. Extract leading terms for quotient basis computation
3. Build multiplication tables from reduced polynomials

### Connect Everything in Main Entry Point
```cpp
// File: src/julia_rur/rur_solver.hpp/cpp
namespace julia_rur {

class RURSolver {
public:
    // Main entry point matching Julia
    RURResult solve(
        const std::vector<Polynomial>& system,
        const SolverOptions& options = SolverOptions()
    );
    
private:
    // Internal methods for each phase
    std::pair<std::vector<PP>, std::vector<mpq_class>> 
    parse_polynomial_system(const std::vector<Polynomial>& system);
    
    // ... other internal methods ...
};

}
```

## Testing Requirements

1. **Unit Tests for Each Component**
   - Quotient basis extraction
   - Univariate parameterization 
   - Separating element strategies
   - Full pipeline tests

2. **Integration Tests**
   - Known systems from Julia test suite
   - Cyclic polynomials
   - Katsura benchmarks

3. **Validation Tests**
   - Compare results with Julia implementation
   - Verify RUR properties (evaluation at roots)

## Performance Considerations

1. **Memory Management**
   - Pre-allocate buffers for modular computations
   - Reuse F4 graph structures across primes

2. **Parallelization Points**
   - Multiple primes in parallel
   - Matrix operations in parameterization
   - Independent polynomial evaluations

3. **Optimization Opportunities**
   - Cache multiplication table results
   - Vectorize modular arithmetic operations
   - Early termination when enough primes

## Documentation Needs

1. **API Documentation**
   - Public interface for RURSolver
   - Options and parameters
   - Result interpretation

2. **Algorithm Documentation**
   - Mathematical background
   - Implementation choices
   - Differences from Julia version

3. **Usage Examples**
   - Simple polynomial systems
   - Benchmarks
   - Integration with other tools