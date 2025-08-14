# RUR-CPP: Rational Univariate Representation Solver

A high-performance C++ implementation of the Rational Univariate Representation (RUR) algorithm for solving systems of polynomial equations. This implementation is based on the Julia package `RationalUnivariateRepresentation.jl` and provides exact solutions to zero-dimensional polynomial systems.

## Key Features

- **Exact Arithmetic**: Uses multi-modular computation with Chinese Remainder Theorem for exact rational solutions
- **High Performance**: Leverages the F4 algorithm for Gröbner basis computation
- **Robust**: Handles systems with up to thousands of solutions (tested on Cyclic-7 with 924 solutions)
- **Flexible Input**: Supports standard polynomial notation

## Installation

### Prerequisites

- C++23 compatible compiler (GCC 14+, Clang 16+)
- CMake 3.16+
- GMP, MPFR libraries
- Eigen3 for numerical root finding
- FLINT 3.3.1 with ARB support (automatically built by setup script)

### Quick Start (Recommended)

The easiest way to build RUR-CPP:

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/rur-cpp.git
cd rur-cpp

# 2. Install system dependencies (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install build-essential cmake git autoconf libtool
sudo apt-get install libgmp-dev libmpfr-dev libeigen3-dev libgtest-dev

# 3. Run the automated setup script (builds FLINT 3.3.1)
./setup-deps.sh

# 4. Build the project
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# 5. Run tests to verify everything works
./rur_tests
```

That's it! The setup script handles all the complex dependency management automatically.

### Platform-Specific Installation

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install build-essential cmake git autoconf libtool
sudo apt-get install libgmp-dev libmpfr-dev libeigen3-dev libgtest-dev
```

#### Fedora/RHEL
```bash
sudo dnf install gcc gcc-c++ make cmake git autoconf libtool
sudo dnf install gmp-devel mpfr-devel eigen3-devel gtest-devel
```

#### macOS (with Homebrew)
```bash
brew install cmake git autoconf automake libtool
brew install gmp mpfr eigen googletest
```

### Manual Build (Advanced)

If you prefer to manage dependencies manually:

```bash
# Build FLINT 3.3.1 manually
mkdir -p external && cd external
wget https://www.flintlib.org/flint-3.3.1.tar.gz
tar xzf flint-3.3.1.tar.gz
cd flint-3.3.1
./configure --prefix=$(pwd)/../flint-install --with-gmp --with-mpfr
make -j$(nproc)
make install
cd ../..

# Build the project
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

## Quick Start

### Example 1: Simple Circle-Line Intersection

```cpp
#include "julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>

int main() {
    // System: x^2 + y^2 = 1 (unit circle)
    //         x - y = 0     (diagonal line)
    // Expected: 2 solutions at (±√2/2, ±√2/2)
    
    std::vector<std::string> polynomials = {
        "x^2 + y^2 - 1",
        "x - y"
    };
    std::vector<std::string> variables = {"x", "y"};
    
    // Solve the system
    auto solution = julia_rur::solve_polynomial_system_enhanced(
        polynomials, variables
    );
    
    if (solution.success) {
        std::cout << "Found " << solution.solutions.size() << " solutions:\n";
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            std::cout << "Solution " << (i+1) << ": (";
            for (size_t j = 0; j < variables.size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << variables[j] << "=" << solution.solutions[i][j].real();
            }
            std::cout << ")\n";
        }
    }
    
    return 0;
}
```

### Example 2: The Cyclic-3 System

```cpp
// Famous benchmark system with 6 solutions
std::vector<std::string> polynomials = {
    "x + y + z",           // Sum = 0
    "x*y + y*z + z*x",     // Symmetric sum of products = 0  
    "x*y*z - 1"            // Product = 1
};
std::vector<std::string> variables = {"x", "y", "z"};

auto solution = julia_rur::solve_polynomial_system_enhanced(
    polynomials, variables
);
// Result: 6 complex solutions (permutations of cube roots of unity)
```

### Example 3: Using the Lower-Level RUR API

```cpp
#include "julia_rur/rur_main_algorithm.hpp"

// For when you need more control over the computation
julia_rur::RURConfig config;
config.verbose = true;  // Show progress

// Compute RUR with rational coefficients
auto result = julia_rur::compute_rational_rur(
    polynomials, variables, config
);

if (result.success) {
    std::cout << "Minimal polynomial degree: " 
              << result.minimal_polynomial.size() - 1 << "\n";
    std::cout << "Number of solutions: " 
              << result.quotient_basis.size() << "\n";
    
    // The RUR gives parametric representation:
    // x_i = g_i(T) / f'(T) where f(T) is the minimal polynomial
}
```

## Supported Systems

The solver handles:
- **Zero-dimensional systems**: Systems with finitely many solutions
- **Systems with rational coefficients**: Exact arithmetic throughout
- **Large systems**: Tested on benchmarks like Cyclic-7 (924 solutions)
- **Systems with multiplicities**: Proper handling of repeated roots

### Benchmark Results

| System | Variables | Solutions | Time (s) |
|--------|-----------|-----------|----------|
| Circle-Line | 2 | 2 | < 0.01 |
| Cyclic-3 | 3 | 6 | < 0.01 |
| Cyclic-4 | 4 | 16 | 0.02 |
| Cyclic-5 | 5 | 70 | 0.1 |
| Cyclic-6 | 6 | 156 | 4 |
| Cyclic-7 | 7 | 924 | ~120 |
| Katsura-4 | 5 | 16 | 0.05 |

## API Reference

### High-Level Solver

```cpp
namespace julia_rur {

// Complete polynomial system solver with numerical roots
EnhancedSolution solve_polynomial_system_enhanced(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const EnhancedSolverConfig& config = {}
);

struct EnhancedSolution {
    bool success;
    std::vector<std::vector<std::complex<double>>> solutions;
    std::vector<int> multiplicities;  // If tracking enabled
    std::string error_message;
};
}
```

### RUR Computation

```cpp
// Compute RUR with exact rational coefficients
RationalRURResult compute_rational_rur(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const RURConfig& config = {}
);

// Single prime computation (faster, modular arithmetic)
std::pair<ModularRURResult, std::vector<int>> compute_modular_rur(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    ModularCoeff prime,
    const RURConfig& config = {},
    const std::vector<int>& separating_element = {}
);
```

### System Analysis

```cpp
// Check if system has finitely many solutions
bool is_zero_dimensional(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables
);

// Get dimension of quotient ring (number of solutions with multiplicity)
int get_quotient_dimension(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables
);
```

## Input Format

Polynomials can be specified in natural notation:
- Standard operations: `+`, `-`, `*`, `^`
- Variables: any alphanumeric string
- Constants: integer or rational numbers
- Parentheses for grouping

Examples:
```
"x^2 + y^2 - 1"           // Circle
"x*y - 1"                 // Hyperbola
"x^3 - 2*x + 1"           // Cubic
"(x + y)^2 - 2*x*y"       // Expanded binomial
```

## Implementation Details

### Algorithm Phases

1. **Polynomial Parsing**: Convert string input to internal representation
2. **Gröbner Basis**: Compute via F4 algorithm using axf4 library
3. **Quotient Basis**: Extract monomial basis of quotient ring
4. **Multiplication Tables**: Build multiplication matrices
5. **Separating Element**: Find linear form that separates all roots
6. **Minimal Polynomial**: Compute characteristic polynomial
7. **Parameterization**: Express each variable as rational function
8. **Multi-Modular CRT**: Reconstruct exact rational coefficients
9. **Root Finding**: Numerical approximation via eigenvalue methods

### Performance Optimizations

- **Separating Element Caching**: Reuses separating element across primes
- **Smart Prime Selection**: Avoids "bad" primes that change system dimension
- **Early CRT Termination**: Stops when coefficients stabilize
- **Efficient F4**: Uses specialized implementation with optimized linear algebra

## Testing

Run the comprehensive test suite:

```bash
cd build-release
./rur_tests                    # Google Test suite
./test_benchmark_systems       # Benchmark systems
./demo_rur examples/cyclic3.txt  # Demo program
```

## Directory Structure

```
rur-cpp/
├── src/
│   ├── julia_rur/         # Main RUR implementation
│   ├── axf4_wrapper.c     # F4 algorithm interface
│   └── ...
├── test/
│   ├── gtest/            # Google Test suite
│   └── ...
├── examples/             # Example polynomial systems
├── external/             # External dependencies (FLINT)
├── vcpkg/               # Package manager
└── build-*/             # Build directories
```

## Contributing

Contributions are welcome! Areas of interest:
- Parallelization of multi-modular computations
- GPU acceleration for linear algebra
- Support for positive-dimensional systems
- Certified/interval arithmetic for numerical roots

## License

This project is open source. See LICENSE file for details.

## Acknowledgments

- Based on algorithms from `RationalUnivariateRepresentation.jl`
- F4 implementation from the axf4 library
- FLINT library for arbitrary precision arithmetic

## References

1. Rouillier, F. (1999). "Solving Zero-Dimensional Systems Through the Rational Univariate Representation"
2. Faugère, J.-C. (1999). "A new efficient algorithm for computing Gröbner bases (F4)"
3. Noro, M., Takeshima, T. (1992). "Risa/Asir - A Computer Algebra System"

## Contact

For questions, issues, or contributions, please open an issue on GitHub.