# RUR C++ API Guide

This guide provides detailed information about using the RUR solver library in your C++ applications.

## Table of Contents

1. [Basic Usage](#basic-usage)
2. [Core API Functions](#core-api-functions)
3. [Data Structures](#data-structures)
4. [Configuration Options](#configuration-options)
5. [Error Handling](#error-handling)
6. [Integration Examples](#integration-examples)
7. [Performance Tips](#performance-tips)

## Basic Usage

### Minimal Example

```cpp
#include <julia_rur/rur_main_algorithm.hpp>
#include <iostream>

int main() {
    using namespace julia_rur;
    
    // Define polynomial system
    std::vector<std::string> polys = {
        "1*x^2+1*y^2+100002",  // x² + y² - 1
        "1*x+100002*y"         // x - y
    };
    std::vector<std::string> vars = {"x", "y"};
    
    // Compute RUR
    auto result = compute_rational_rur(polys, vars);
    
    if (result.success) {
        std::cout << format_rur_result(result, vars) << std::endl;
    }
    
    return 0;
}
```

### Compilation

```bash
g++ -std=c++17 myapp.cpp -lrur -lgmpxx -lgmp -leigen3 -o myapp
```

## Core API Functions

### Main Computation Functions

#### `compute_rational_rur`

Computes RUR with exact rational coefficients using multi-modular arithmetic.

```cpp
RationalRURResult compute_rational_rur(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const RURConfig& config = RURConfig()
);
```

**Parameters:**
- `polynomials`: System of polynomial equations in string format
- `variables`: Ordered list of variable names
- `config`: Optional configuration settings

**Returns:** `RationalRURResult` with exact rational coefficients

**Example:**
```cpp
RURConfig config;
config.verbose = true;
config.initial_prime_bits = 20;  // Use smaller primes

auto result = compute_rational_rur(polys, vars, config);
```

#### `compute_modular_rur`

Computes RUR modulo a single prime (faster for when exact rationals aren't needed).

```cpp
ModularRURResult compute_modular_rur(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    ModularCoeff prime,
    const RURConfig& config = RURConfig()
);
```

**Use when:**
- You only need solutions modulo a prime
- You're doing a feasibility check
- You want to test if a system is solvable

**Example:**
```cpp
auto mod_result = compute_modular_rur(polys, vars, 1048583);
if (mod_result.success) {
    std::cout << "System is solvable modulo 1048583\n";
}
```

### Utility Functions

#### `is_zero_dimensional_system`

Quick check if system has finitely many solutions.

```cpp
bool is_zero_dimensional_system(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    ModularCoeff prime = 100003
);
```

**Example:**
```cpp
if (!is_zero_dimensional_system(polys, vars)) {
    std::cerr << "System has infinitely many solutions!\n";
    return;
}
```

#### `compute_quotient_dimension`

Returns the number of solutions (counting multiplicities).

```cpp
int compute_quotient_dimension(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    ModularCoeff prime = 100003
);
```

**Example:**
```cpp
int num_solutions = compute_quotient_dimension(polys, vars);
std::cout << "System has " << num_solutions << " solutions\n";
```

#### `format_rur_result`

Pretty-prints the RUR result.

```cpp
std::string format_rur_result(
    const RationalRURResult& result,
    const std::vector<std::string>& variables
);
```

## Data Structures

### `RationalRURResult`

Contains the complete rational univariate representation.

```cpp
struct RationalRURResult {
    // Minimal polynomial f(T) coefficients (low to high degree)
    std::vector<mpq_class> minimal_polynomial;
    
    // Numerator polynomials g_i(T) for each variable
    std::vector<std::vector<mpq_class>> numerators;
    
    // Common denominator f'(T)
    mpq_class denominator_derivative;
    
    // Monomial basis of quotient ring
    std::vector<PP> quotient_basis;
    
    // Success flag
    bool success;
    
    // Error message if failed
    std::string error_message;
};
```

**Usage Example:**
```cpp
if (result.success) {
    // Access minimal polynomial
    std::cout << "Degree: " << result.minimal_polynomial.size() - 1 << "\n";
    
    // Access leading coefficient
    auto leading = result.minimal_polynomial.back();
    std::cout << "Leading coeff: " << leading << "\n";
    
    // Check number of variables
    std::cout << "Variables: " << result.numerators.size() << "\n";
}
```

### `ModularRURResult`

Result from single prime computation.

```cpp
struct ModularRURResult {
    ModularCoeff prime;
    MinimalPolynomialResult minimal_polynomial;
    std::vector<BivariateResult> parameterizations;
    std::vector<PP> quotient_basis;
    bool success;
};
```

### `RURConfig`

Configuration options for RUR computation.

```cpp
struct RURConfig {
    size_t initial_prime_bits = 31;     // Starting prime size
    size_t max_prime_bits = 63;         // Maximum prime size
    size_t num_threads = 1;             // Parallelization (future)
    SeparatingStrategy separating_strategy = SeparatingStrategy::CURRENT;
    bool verbose = false;               // Progress output
};
```

### `PP` (Power Product)

Represents a monomial as exponent vector.

```cpp
using PP = std::vector<uint32_t>;

// Example: x²y³ in variables [x,y,z]
PP monomial = {2, 3, 0};  // x^2 * y^3 * z^0
```

## Configuration Options

### Separating Element Strategies

```cpp
enum class SeparatingStrategy {
    CURRENT,        // Default: smart heuristic
    RANDOM,         // Random linear combinations
    L0_NORM,        // Prefer sparse combinations
    DETERMINISTIC,  // Systematic search
    MRON_0L        // Reverse L0 norm
};

// Usage
RURConfig config;
config.separating_strategy = SeparatingStrategy::L0_NORM;
```

### Prime Size Selection

```cpp
// For systems with small coefficients
config.initial_prime_bits = 20;  // ~1 million

// For systems with large coefficients
config.initial_prime_bits = 31;  // ~2 billion (default)

// For very large systems
config.max_prime_bits = 63;      // Up to 64-bit primes
```

### Verbosity Levels

```cpp
config.verbose = false;  // Silent (default)
config.verbose = true;   // Progress messages
```

## Error Handling

### Common Errors and Solutions

```cpp
RationalRURResult result = compute_rational_rur(polys, vars);

if (!result.success) {
    // Check specific error
    if (result.error_message.find("zero-dimensional") != std::string::npos) {
        std::cerr << "System has infinitely many solutions\n";
    } else if (result.error_message.find("no solutions") != std::string::npos) {
        std::cerr << "System is inconsistent\n";
    } else if (result.error_message.find("prime") != std::string::npos) {
        std::cerr << "Prime selection failed\n";
    } else {
        std::cerr << "Error: " << result.error_message << "\n";
    }
}
```

### Exception Safety

```cpp
try {
    auto result = compute_rational_rur(polys, vars);
    // Process result
} catch (const std::bad_alloc& e) {
    std::cerr << "Out of memory\n";
} catch (const std::exception& e) {
    std::cerr << "Unexpected error: " << e.what() << "\n";
}
```

## Integration Examples

### With Eigen for Numerical Solutions

```cpp
#include <Eigen/Dense>

// Extract double coefficients
std::vector<double> get_double_coeffs(const std::vector<mpq_class>& rationals) {
    std::vector<double> doubles;
    for (const auto& r : rationals) {
        doubles.push_back(r.get_d());
    }
    return doubles;
}

// Find numerical roots
std::vector<std::complex<double>> find_roots(const RationalRURResult& rur) {
    auto coeffs = get_double_coeffs(rur.minimal_polynomial);
    int degree = coeffs.size() - 1;
    
    // Build companion matrix
    Eigen::MatrixXd C(degree, degree);
    C.setZero();
    
    // Fill companion matrix
    for (int i = 1; i < degree; ++i) {
        C(i, i-1) = 1.0;
    }
    
    double lead = coeffs[degree];
    for (int i = 0; i < degree; ++i) {
        C(i, degree-1) = -coeffs[i] / lead;
    }
    
    // Compute eigenvalues
    Eigen::EigenSolver<Eigen::MatrixXd> solver(C);
    
    std::vector<std::complex<double>> roots;
    for (int i = 0; i < solver.eigenvalues().size(); ++i) {
        roots.push_back(solver.eigenvalues()[i]);
    }
    
    return roots;
}
```

### With GMP for Exact Arithmetic

```cpp
#include <gmpxx.h>

// Work with exact rationals
void process_exact_solution(const RationalRURResult& rur) {
    // Get exact leading coefficient
    mpq_class leading = rur.minimal_polynomial.back();
    
    // Normalize to monic polynomial
    std::vector<mpq_class> monic_poly;
    for (const auto& coeff : rur.minimal_polynomial) {
        monic_poly.push_back(coeff / leading);
    }
    
    // Now monic_poly has leading coefficient 1
    // Can be passed to algebraic number libraries
}
```

### Export to Other Systems

```cpp
// Export to Maple format
std::string to_maple(const RationalRURResult& rur, 
                    const std::vector<std::string>& vars) {
    std::stringstream maple;
    
    // Minimal polynomial
    maple << "f := ";
    for (int i = rur.minimal_polynomial.size()-1; i >= 0; --i) {
        if (i < rur.minimal_polynomial.size()-1) maple << " + ";
        maple << "(" << rur.minimal_polynomial[i] << ")*T^" << i;
    }
    maple << ";\n";
    
    // Solve
    maple << "sols := [solve(f=0, T)];\n";
    
    // Parameterizations
    maple << "fprime := diff(f, T);\n";
    for (size_t i = 0; i < vars.size(); ++i) {
        maple << vars[i] << "_vals := map(t -> ";
        // Output numerator evaluation
        maple << "subs(T=t, ...) / subs(T=t, fprime), sols);\n";
    }
    
    return maple.str();
}

// Export to Python/SymPy
std::string to_python(const RationalRURResult& rur) {
    std::stringstream py;
    
    py << "import sympy as sp\n";
    py << "T = sp.Symbol('T')\n\n";
    
    // Minimal polynomial
    py << "f = ";
    for (int i = rur.minimal_polynomial.size()-1; i >= 0; --i) {
        if (i > 0) {
            py << rur.minimal_polynomial[i] << "*T**" << i;
            if (i > 1) py << " + ";
        } else {
            py << " + " << rur.minimal_polynomial[0];
        }
    }
    py << "\n\n";
    
    py << "# Find roots\n";
    py << "roots = sp.solve(f, T)\n";
    py << "print(f'Found {len(roots)} roots')\n";
    
    return py.str();
}
```

## Performance Tips

### 1. Pre-check Dimensionality

```cpp
// Always check before full computation
if (!is_zero_dimensional_system(polys, vars)) {
    // Skip RUR computation
    return handle_infinite_solutions();
}
```

### 2. Use Modular Computation First

```cpp
// Quick feasibility check
auto mod_result = compute_modular_rur(polys, vars, 100003);
if (!mod_result.success) {
    // System likely has no solutions
    return;
}

// Now do full computation if needed
auto full_result = compute_rational_rur(polys, vars);
```

### 3. Choose Appropriate Prime Sizes

```cpp
// For small systems (< 10 solutions)
config.initial_prime_bits = 20;

// For medium systems (10-100 solutions)
config.initial_prime_bits = 25;

// For large systems (> 100 solutions)
config.initial_prime_bits = 31;  // default
```

### 4. Reuse Sessions for Multiple Systems

```cpp
// If solving many similar systems
class RURSolver {
    RURConfig config;
    
public:
    RURSolver() {
        config.verbose = false;
        config.initial_prime_bits = 25;
    }
    
    RationalRURResult solve(const std::vector<std::string>& polys,
                           const std::vector<std::string>& vars) {
        return compute_rational_rur(polys, vars, config);
    }
};
```

### 5. Memory Management

```cpp
// For large systems, process results immediately
{
    auto result = compute_rational_rur(large_system, vars);
    // Extract only what you need
    std::vector<double> roots = find_roots(result);
    // result goes out of scope, memory freed
}
```

## Advanced Topics

### Custom Modular Arithmetic

```cpp
// Convert between rational and modular representations
ModularCoeff to_modular(const mpq_class& rational, ModularCoeff prime) {
    mpz_class num = rational.get_num();
    mpz_class den = rational.get_den();
    
    // Compute num * den^(-1) mod prime
    // ... (implementation using modular inverse)
}
```

### Parallel Computation (Future)

```cpp
// Future API for parallel multi-modular computation
RURConfig config;
config.num_threads = 4;  // Use 4 threads for different primes
```

### Custom Separating Elements

```cpp
// Future API for custom separating element
config.custom_separating_element = {1, 2, 3};  // x + 2y + 3z
```