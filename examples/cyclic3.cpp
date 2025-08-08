/**
 * @file cyclic3.cpp
 * @brief Example: Solving the Cyclic-3 polynomial system using RUR
 * 
 * The Cyclic-3 system is a famous benchmark in polynomial system solving:
 *   x + y + z = 0
 *   xy + yz + zx = 0
 *   xyz - 1 = 0
 * 
 * This system has 6 complex solutions corresponding to the 3! = 6
 * permutations of the three cube roots of unity scaled appropriately.
 */

#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <complex>
#include <vector>

using namespace julia_rur;

// Helper function to evaluate a polynomial at a point
std::complex<double> evaluate_poly(const std::vector<mpq_class>& coeffs, 
                                  std::complex<double> x) {
    std::complex<double> result = 0.0;
    std::complex<double> x_power = 1.0;
    
    for (const auto& coeff : coeffs) {
        result += coeff.get_d() * x_power;
        x_power *= x;
    }
    
    return result;
}

int main() {
    std::cout << "Cyclic-3 System Example\n";
    std::cout << "======================\n\n";
    
    // Define the Cyclic-3 system
    // Using modular arithmetic where 100002 ≡ -1 (mod 100003)
    std::vector<std::string> polynomials = {
        "1*x+1*y+1*z",           // x + y + z = 0
        "1*x*y+1*y*z+1*z*x",     // xy + yz + zx = 0
        "1*x*y*z+100002"         // xyz - 1 = 0
    };
    
    std::vector<std::string> variables = {"x", "y", "z"};
    
    std::cout << "System of equations:\n";
    std::cout << "  x + y + z = 0\n";
    std::cout << "  xy + yz + zx = 0\n";
    std::cout << "  xyz = 1\n\n";
    
    // Check dimensionality
    int dim = compute_quotient_dimension(polynomials, variables);
    std::cout << "Quotient ring dimension: " << dim << "\n";
    std::cout << "This means the system has " << dim << " solutions (counting multiplicity).\n\n";
    
    // Compute RUR
    std::cout << "Computing Rational Univariate Representation...\n";
    
    RURConfig config;
    config.verbose = true;  // Show progress
    
    RationalRURResult result = compute_rational_rur(polynomials, variables, config);
    
    if (!result.success) {
        std::cerr << "Error: " << result.error_message << "\n";
        return 1;
    }
    
    std::cout << "\n" << format_rur_result(result, variables) << "\n";
    
    // Display the minimal polynomial
    std::cout << "Minimal polynomial degree: " << result.minimal_polynomial.size() - 1 << "\n";
    std::cout << "Coefficients (from constant to highest degree):\n";
    for (size_t i = 0; i < result.minimal_polynomial.size(); ++i) {
        std::cout << "  f[" << i << "] = " << result.minimal_polynomial[i] << "\n";
    }
    
    // Mathematical interpretation
    std::cout << "\nMathematical Background:\n";
    std::cout << "========================\n";
    std::cout << "The Cyclic-3 system is related to the roots of unity.\n";
    std::cout << "If ω = e^(2πi/3) is a primitive cube root of unity,\n";
    std::cout << "then the solutions involve permutations of {ω^a, ω^b, ω^c}\n";
    std::cout << "scaled to satisfy xyz = 1.\n\n";
    
    std::cout << "The RUR gives us:\n";
    std::cout << "- A univariate polynomial f(T) of degree " << dim << "\n";
    std::cout << "- Rational parameterizations x(T), y(T), z(T)\n";
    std::cout << "- Each root T_i of f(T) gives a solution (x_i, y_i, z_i)\n\n";
    
    // Numerical approximation (if desired)
    std::cout << "To find the actual solutions:\n";
    std::cout << "1. Find all roots of the minimal polynomial f(T)\n";
    std::cout << "2. For each root T_i, evaluate:\n";
    std::cout << "   x_i = g_x(T_i) / f'(T_i)\n";
    std::cout << "   y_i = g_y(T_i) / f'(T_i)\n";
    std::cout << "   z_i = g_z(T_i) / f'(T_i)\n";
    std::cout << "where g_x, g_y, g_z are the numerator polynomials.\n";
    
    return 0;
}