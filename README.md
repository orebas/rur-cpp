# RUR-CPP: Rational Univariate Representation Solver

A high-performance C++ implementation of the Rational Univariate Representation (RUR) algorithm for solving systems of polynomial equations. This implementation is based on the Julia package `RationalUnivariateRepresentation.jl` and provides exact solutions to zero-dimensional polynomial systems.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Format](#input-format)
- [Output Format](#output-format)
- [API Reference](#api-reference)
- [Examples](#examples)
- [Integration with Univariate Solvers](#integration-with-univariate-solvers)
- [Algorithm Details](#algorithm-details)

## Overview

The RUR algorithm transforms a system of multivariate polynomial equations into a univariate representation, making it possible to find all solutions by solving a single univariate polynomial. This implementation features:

- **Exact Arithmetic**: Uses multi-modular computation with Chinese Remainder Theorem
- **High Performance**: Leverages the F4 algorithm for Gröbner basis computation
- **No Re-entrancy**: Innovative bivariate algorithm avoids multiple Gröbner basis computations
- **Flexible Input**: Supports various polynomial formats and variable orderings

## Installation

### Prerequisites

- C++17 compatible compiler (GCC 7+, Clang 5+)
- CMake 3.14+
- GMP library for arbitrary precision arithmetic
- Eigen3 for linear algebra operations

### Building from Source

```bash
git clone https://github.com/yourusername/rur-cpp.git
cd rur-cpp
mkdir build && cd build
cmake ..
make -j4
```

## Quick Start

```cpp
#include "julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    // Define polynomial system: x^2 + y^2 - 1 = 0, x - y = 0
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2+100002",  // In modular form (100002 ≡ -1 mod 100003)
        "1*x+100002*y"         // x - y
    };
    std::vector<std::string> variables = {"x", "y"};
    
    // Compute RUR
    RURConfig config;
    config.verbose = true;
    
    RationalRURResult result = compute_rational_rur(polynomials, variables, config);
    
    if (result.success) {
        std::cout << format_rur_result(result, variables) << std::endl;
    }
    
    return 0;
}
```

## Input Format

### Polynomial String Format

Polynomials are specified as strings with the following syntax:

```
coefficient*variable^power + coefficient*variable^power + ...
```

Rules:
- Coefficients must be explicit (use `1*x` not just `x`)
- Powers use `^` notation: `x^2`, `y^3`
- Multiple variables: `2*x^2*y^3`
- Constants are written as coefficients without variables
- For modular arithmetic, negative numbers are represented as `prime + negative_value`

Examples:
```cpp
"1*x^2+1*y^2+100002"     // x² + y² - 1 (where 100002 ≡ -1 mod 100003)
"1*x*y+5*x+100000*y"     // xy + 5x - 3y (where 100000 ≡ -3 mod 100003)
"1*x^3+1*y^3+1*z^3+100002"  // x³ + y³ + z³ - 1
```

### Variable Ordering

Variables should be provided in the order you want them to appear in the output:

```cpp
std::vector<std::string> variables = {"x", "y", "z"};
```

## Output Format

### RationalRURResult Structure

```cpp
struct RationalRURResult {
    std::vector<mpq_class> minimal_polynomial;     // f(T) coefficients
    std::vector<std::vector<mpq_class>> numerators; // gi(T) for each xi
    mpq_class denominator_derivative;              // f'(T)
    std::vector<PP> quotient_basis;                // Monomial basis
    bool success;
    std::string error_message;
};
```

The RUR provides a parametric representation where:
- Each variable xi = gi(T)/f'(T)
- T is the separating element (usually a variable or linear combination)
- f(T) is the minimal polynomial whose roots give all solutions

### Example Output

For the system {x² + y² - 1 = 0, x - y = 0}:

```
Rational Univariate Representation:
===================================

Minimal polynomial f(T):
f(T) = T^2 - 1/2

Quotient ring dimension: 2

Parameterizations:
x = (T) / f'(T)
y = (T) / f'(T)

Solutions: T = ±√(1/2), giving (x,y) = (±√(1/2), ±√(1/2))
```

## API Reference

### Main Functions

#### `compute_rational_rur`

Computes the complete RUR with exact rational coefficients.

```cpp
RationalRURResult compute_rational_rur(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const RURConfig& config = RURConfig()
);
```

#### `compute_modular_rur`

Computes RUR modulo a single prime (faster, for when modular arithmetic suffices).

```cpp
ModularRURResult compute_modular_rur(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    ModularCoeff prime,
    const RURConfig& config = RURConfig()
);
```

#### `is_zero_dimensional_system`

Quick check if the system has finitely many solutions.

```cpp
bool is_zero_dimensional_system(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    ModularCoeff prime = 100003
);
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

### Configuration Options

```cpp
struct RURConfig {
    size_t initial_prime_bits = 31;     // Starting prime size
    size_t max_prime_bits = 63;         // Maximum prime size
    size_t num_threads = 1;             // For future parallelization
    SeparatingStrategy separating_strategy = SeparatingStrategy::CURRENT;
    bool verbose = false;               // Enable progress output
};
```

Separating strategies:
- `CURRENT`: Try last variable first, then systematic search
- `RANDOM`: Random linear combinations
- `L0_NORM`: Prefer sparse linear forms
- `DETERMINISTIC`: Systematic power-based search

## Examples

### Example 1: Circle-Line Intersection

```cpp
// System: x² + y² = 1, x = y
// Solutions: (1/√2, 1/√2) and (-1/√2, -1/√2)

std::vector<std::string> polynomials = {
    "1*x^2+1*y^2+100002",  // x² + y² - 1
    "1*x+100002*y"         // x - y
};
std::vector<std::string> variables = {"x", "y"};

auto result = compute_rational_rur(polynomials, variables);
// Result: f(T) = T² - 1/2, x = T/f'(T), y = T/f'(T)
```

### Example 2: Cyclic-3 System

```cpp
// Famous benchmark: x + y + z = 0, xy + yz + zx = 0, xyz = 1
std::vector<std::string> polynomials = {
    "1*x+1*y+1*z",           // x + y + z
    "1*x*y+1*y*z+1*z*x",     // xy + yz + zx
    "1*x*y*z+100002"         // xyz - 1
};
std::vector<std::string> variables = {"x", "y", "z"};

auto result = compute_rational_rur(polynomials, variables);
// Result: 6 solutions corresponding to permutations of cube roots
```

### Example 3: Robot Kinematics

```cpp
// 2D robot arm with constraints
std::vector<std::string> polynomials = {
    "1*c1^2+1*s1^2+100002",     // cos²(θ₁) + sin²(θ₁) = 1
    "1*c2^2+1*s2^2+100002",     // cos²(θ₂) + sin²(θ₂) = 1
    "2*c1+1*c1*c2+100002*s1*s2+100001",  // Forward kinematics x = 2
    "2*s1+1*s1*c2+1*c1*s2"      // Forward kinematics y = 0
};
std::vector<std::string> variables = {"c1", "s1", "c2", "s2"};

auto result = compute_rational_rur(polynomials, variables);
```

## Integration with Univariate Solvers

Once you have the RUR, solving the system reduces to finding roots of the univariate polynomial f(T). Here's how to integrate with various solvers:

### Using Eigen for Numerical Roots

```cpp
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

std::vector<std::complex<double>> solve_rur_numerically(
    const RationalRURResult& rur
) {
    size_t degree = rur.minimal_polynomial.size() - 1;
    
    // Build companion matrix
    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();
    
    // Fill subdiagonal with 1s
    for (size_t i = 1; i < degree; ++i) {
        companion(i, i-1) = 1.0;
    }
    
    // Last column contains -coefficients/leading_coeff
    double leading = rur.minimal_polynomial[degree].get_d();
    for (size_t i = 0; i < degree; ++i) {
        companion(i, degree-1) = -rur.minimal_polynomial[i].get_d() / leading;
    }
    
    // Compute eigenvalues
    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    
    std::vector<std::complex<double>> roots;
    for (int i = 0; i < solver.eigenvalues().size(); ++i) {
        roots.push_back(solver.eigenvalues()[i]);
    }
    
    return roots;
}

// Evaluate parameterizations at roots
std::vector<std::vector<std::complex<double>>> evaluate_solutions(
    const RationalRURResult& rur,
    const std::vector<std::complex<double>>& roots
) {
    std::vector<std::vector<std::complex<double>>> solutions;
    
    for (const auto& root : roots) {
        std::vector<std::complex<double>> point;
        
        // Evaluate f'(T) at root
        std::complex<double> fprime = evaluate_derivative(
            rur.minimal_polynomial, root
        );
        
        // For each variable, evaluate gi(T)/f'(T)
        for (const auto& numerator : rur.numerators) {
            std::complex<double> gi = evaluate_polynomial(numerator, root);
            point.push_back(gi / fprime);
        }
        
        solutions.push_back(point);
    }
    
    return solutions;
}
```

### Using FLINT for Exact Algebraic Roots

```cpp
#include <flint/fmpq_poly.h>

void solve_rur_exactly(const RationalRURResult& rur) {
    // Convert to FLINT polynomial
    fmpq_poly_t f;
    fmpq_poly_init(f);
    
    for (size_t i = 0; i < rur.minimal_polynomial.size(); ++i) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_mpq(coeff, rur.minimal_polynomial[i].get_mpq_t());
        fmpq_poly_set_coeff_fmpq(f, i, coeff);
        fmpq_clear(coeff);
    }
    
    // Factor or find roots using FLINT's algebraic number theory
    // ... (depends on your exact needs)
    
    fmpq_poly_clear(f);
}
```

### Symbolic Output for Computer Algebra Systems

```cpp
std::string export_to_maple(const RationalRURResult& rur,
                           const std::vector<std::string>& vars) {
    std::stringstream maple;
    
    // Define minimal polynomial
    maple << "f := ";
    for (size_t i = rur.minimal_polynomial.size(); i > 0; --i) {
        if (i > 1) maple << rur.minimal_polynomial[i-1] << "*T^" << (i-1);
        else maple << rur.minimal_polynomial[0];
        if (i > 1) maple << " + ";
    }
    maple << ";\n\n";
    
    // Define parameterizations
    maple << "fprime := diff(f, T);\n";
    for (size_t i = 0; i < vars.size(); ++i) {
        maple << vars[i] << "_param := (";
        // Output numerator polynomial
        maple << ") / fprime;\n";
    }
    
    // Solve command
    maple << "solutions := [solve(f = 0, T)];\n";
    
    return maple.str();
}
```

## Algorithm Details

The RUR algorithm consists of several key phases:

1. **Gröbner Basis Computation**: Uses the F4 algorithm for efficiency
2. **Quotient Basis Extraction**: Finds monomial basis of the quotient ring
3. **Multiplication Table Construction**: Precomputes multiplication operations
4. **Separating Element Selection**: Finds a linear form that separates all solutions
5. **Minimal Polynomial Computation**: Uses linear algebra in the quotient ring
6. **Bivariate Parameterization**: Extends to find rational parameterizations
7. **Multi-Modular Lifting**: Uses CRT for exact rational reconstruction

The key innovation is the bivariate algorithm that avoids computing multiple Gröbner bases, making the method significantly more efficient than traditional approaches.

## Performance Considerations

- **Prime Selection**: Larger primes give more accurate results but slower computation
- **Early Termination**: The multi-modular loop can terminate early if reconstruction succeeds
- **Memory Usage**: Scales with the dimension of the quotient ring (number of solutions)
- **Parallelization**: Different primes can be computed in parallel (future enhancement)

## Troubleshooting

### Common Issues

1. **"System has no solutions"**: The Gröbner basis is {1}
2. **"Not zero-dimensional"**: The system has infinitely many solutions
3. **"Prime too small"**: F4 requires primes > 2^16
4. **Memory exhaustion**: Try smaller primes or fewer variables

### Debug Output

Enable verbose mode to see detailed progress:

```cpp
RURConfig config;
config.verbose = true;
```

## License

This implementation is based on the Julia package RationalUnivariateRepresentation.jl and follows its algorithmic approach. See LICENSE file for details.

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.

## References

1. Rouillier, F. (1999). "Solving Zero-Dimensional Systems Through the Rational Univariate Representation"
2. Faugère, J.-C. (1999). "A new efficient algorithm for computing Gröbner bases (F4)"
3. RationalUnivariateRepresentation.jl - Julia implementation by [authors]