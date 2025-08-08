#pragma once

#include "rur_main_algorithm.hpp"
#include "numerical_roots_eigen.hpp"
#include <vector>
#include <string>
#include <complex>

namespace julia_rur {

/**
 * @brief Complete solution for a polynomial system
 */
struct PolynomialSystemSolution {
    bool success;
    std::string error_message;
    
    // Solutions: each inner vector contains values for all variables
    std::vector<std::vector<std::complex<double>>> solutions;
    
    // Whether each solution is real (all components have zero imaginary part)
    std::vector<bool> is_real_solution;
    
    // Variable names for reference
    std::vector<std::string> variable_names;
    
    // Additional information
    int quotient_dimension;  // Dimension of the quotient ring
    std::vector<mpq_class> minimal_polynomial;  // Coefficients of minimal polynomial
};

/**
 * @brief Main entry point for solving polynomial systems
 * 
 * This function combines:
 * 1. RUR computation for algebraic solving
 * 2. Numerical root finding for the minimal polynomial
 * 3. Back-substitution to find all variable values
 * 
 * @param polynomials Vector of polynomial strings (e.g., "x^2 + y^2 - 1")
 * @param variables Vector of variable names (e.g., ["x", "y"])
 * @param config Configuration options
 * @return Complete solution with all roots
 */
PolynomialSystemSolution solve_polynomial_system_complete(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const RURConfig& config = RURConfig()
);

/**
 * @brief Pretty-print the solution
 */
void print_solution(const PolynomialSystemSolution& solution, std::ostream& out = std::cout);

/**
 * @brief Convert solution to string format
 */
std::string solution_to_string(const PolynomialSystemSolution& solution);

} // namespace julia_rur