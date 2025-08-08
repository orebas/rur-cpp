#pragma once

#include "rur_main_algorithm.hpp"
#include <unsupported/Eigen/Polynomials>
#include <vector>
#include <complex>

namespace julia_rur {

/**
 * @brief Result of numerical root finding
 */
struct NumericalSolution {
    std::vector<std::vector<std::complex<double>>> solutions;  // Each solution is a vector of values for each variable
    std::vector<bool> is_real;                                 // Whether each solution is real
    bool success;
    std::string error_message;
};

/**
 * @brief Find roots of polynomial using Eigen
 * 
 * @param polynomial_coeffs Coefficients of univariate polynomial (low to high degree)
 * @return Complex roots
 */
std::vector<std::complex<double>> find_polynomial_roots(
    const std::vector<mpq_class>& polynomial_coeffs
);

/**
 * @brief Back-substitute roots to find all variable values
 * 
 * Given roots of the minimal polynomial and parameterizations,
 * compute values for all variables
 */
NumericalSolution back_substitute_roots(
    const RationalRURResult& rur_result,
    const std::vector<std::string>& variables
);

/**
 * @brief Complete polynomial system solver
 * 
 * Combines RUR computation with numerical root finding
 */
NumericalSolution solve_polynomial_system(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const RURConfig& config = RURConfig()
);

} // namespace julia_rur