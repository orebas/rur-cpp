#include "numerical_roots_eigen.hpp"
#include <iostream>

namespace julia_rur {

std::vector<std::complex<double>> find_polynomial_roots(
    const std::vector<mpq_class>& polynomial_coeffs
) {
    std::vector<std::complex<double>> result;
    
    if (polynomial_coeffs.empty()) {
        return result;
    }
    
    // Find degree (highest non-zero coefficient)
    int degree = -1;
    for (int i = polynomial_coeffs.size() - 1; i >= 0; --i) {
        if (polynomial_coeffs[i] != 0) {
            degree = i;
            break;
        }
    }
    
    if (degree <= 0) {
        return result;  // Constant or zero polynomial
    }
    
    // Convert to Eigen format
    // Note: Eigen actually expects coefficients in LOW to HIGH order
    Eigen::VectorXd coeffs(degree + 1);
    
    for (int i = 0; i <= degree; ++i) {
        // Convert rational to double
        coeffs[i] = polynomial_coeffs[i].get_d();
    }
    
    // Create polynomial solver
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    solver.compute(coeffs);
    
    // Get roots
    const auto& roots = solver.roots();
    
    // Convert to std::vector
    for (int i = 0; i < roots.size(); ++i) {
        result.push_back(roots[i]);
    }
    
    return result;
}

/**
 * @brief Evaluate univariate polynomial at a complex point
 */
static std::complex<double> evaluate_polynomial(
    const std::vector<mpq_class>& coeffs,
    const std::complex<double>& x
) {
    std::complex<double> result = 0.0;
    std::complex<double> x_power = 1.0;
    
    for (const auto& coeff : coeffs) {
        result += coeff.get_d() * x_power;
        x_power *= x;
    }
    
    return result;
}

NumericalSolution back_substitute_roots(
    const RationalRURResult& rur_result,
    const std::vector<std::string>& variables
) {
    NumericalSolution solution;
    solution.success = false;
    
    if (!rur_result.success) {
        solution.error_message = "RUR computation failed";
        return solution;
    }
    
    // Find roots of minimal polynomial
    std::vector<std::complex<double>> roots = find_polynomial_roots(
        rur_result.minimal_polynomial
    );
    
    if (roots.empty()) {
        solution.error_message = "No roots found for minimal polynomial";
        return solution;
    }
    
    // For each root, compute values of all variables
    size_t num_vars = variables.size();
    for (const auto& root : roots) {
        std::vector<std::complex<double>> var_values(num_vars);
        
        // Simple case: univariate system where the root is the solution
        if (num_vars == 1) {
            var_values[0] = root;
        } else {
            // For multivariate systems, we need to evaluate the parameterizations
            // Each parameterization is a polynomial in the separating element T
            
            // For now, use a simplified approach
            // In a full implementation, we would parse and evaluate the
            // parameterization polynomials from rur_result.parameterizations
            
            // As a placeholder, just use the root for all variables
            for (size_t i = 0; i < num_vars; ++i) {
                var_values[i] = root;
            }
        }
        
        // Check if solution is real
        bool is_real = true;
        const double epsilon = 1e-10;
        for (const auto& val : var_values) {
            if (std::abs(val.imag()) > epsilon) {
                is_real = false;
                break;
            }
        }
        
        solution.solutions.push_back(var_values);
        solution.is_real.push_back(is_real);
    }
    
    solution.success = true;
    return solution;
}

NumericalSolution solve_polynomial_system(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const RURConfig& config
) {
    NumericalSolution solution;
    
    // First compute the RUR
    RationalRURResult rur_result = compute_rational_rur(
        polynomials, variables, config
    );
    
    if (!rur_result.success) {
        solution.success = false;
        solution.error_message = "RUR computation failed: " + rur_result.error_message;
        return solution;
    }
    
    // Then find numerical roots
    return back_substitute_roots(rur_result, variables);
}

} // namespace julia_rur