/**
 * @file integrate_with_solver.cpp
 * @brief Example: How to integrate RUR output with univariate solvers
 * 
 * This example shows how to:
 * 1. Compute the RUR of a polynomial system
 * 2. Extract the univariate polynomial
 * 3. Solve it numerically using Eigen
 * 4. Back-substitute to get the multivariate solutions
 */

#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <complex>
#include <vector>
#include <iomanip>

using namespace julia_rur;

// Convert rational polynomial to floating-point
std::vector<double> to_double_coeffs(const std::vector<mpq_class>& rational_coeffs) {
    std::vector<double> double_coeffs;
    for (const auto& coeff : rational_coeffs) {
        double_coeffs.push_back(coeff.get_d());
    }
    return double_coeffs;
}

// Build companion matrix for polynomial
Eigen::MatrixXd build_companion_matrix(const std::vector<double>& coeffs) {
    int n = coeffs.size() - 1;  // degree
    Eigen::MatrixXd companion(n, n);
    companion.setZero();
    
    // Subdiagonal of ones
    for (int i = 1; i < n; ++i) {
        companion(i, i-1) = 1.0;
    }
    
    // Last column: -c_i/c_n
    double leading = coeffs[n];
    for (int i = 0; i < n; ++i) {
        companion(i, n-1) = -coeffs[i] / leading;
    }
    
    return companion;
}

// Evaluate polynomial at complex point
std::complex<double> evaluate_polynomial(const std::vector<double>& coeffs,
                                       const std::complex<double>& x) {
    std::complex<double> result = 0.0;
    std::complex<double> x_power = 1.0;
    
    for (double coeff : coeffs) {
        result += coeff * x_power;
        x_power *= x;
    }
    
    return result;
}

// Evaluate polynomial derivative
std::complex<double> evaluate_derivative(const std::vector<double>& coeffs,
                                        const std::complex<double>& x) {
    std::complex<double> result = 0.0;
    std::complex<double> x_power = 1.0;
    
    for (size_t i = 1; i < coeffs.size(); ++i) {
        result += i * coeffs[i] * x_power;
        x_power *= x;
    }
    
    return result;
}

int main() {
    std::cout << "RUR Integration with Univariate Solver Example\n";
    std::cout << "============================================\n\n";
    
    // Example system: Two ellipses intersection
    // x²/4 + y²/9 = 1
    // x²/9 + y²/4 = 1
    // Rewritten as: 9x² + 4y² - 36 = 0, 4x² + 9y² - 36 = 0
    
    std::vector<std::string> polynomials = {
        "9*x^2+4*y^2+99967",     // 9x² + 4y² - 36 (99967 ≡ -36 mod 100003)
        "4*x^2+9*y^2+99967"      // 4x² + 9y² - 36
    };
    
    std::vector<std::string> variables = {"x", "y"};
    
    std::cout << "System: Two ellipses intersection\n";
    std::cout << "  x²/4 + y²/9 = 1\n";
    std::cout << "  x²/9 + y²/4 = 1\n\n";
    
    // Step 1: Compute RUR
    std::cout << "Step 1: Computing RUR...\n";
    
    RURConfig config;
    config.verbose = false;
    
    RationalRURResult rur = compute_rational_rur(polynomials, variables, config);
    
    if (!rur.success) {
        std::cerr << "RUR computation failed: " << rur.error_message << "\n";
        return 1;
    }
    
    std::cout << "RUR computed successfully!\n";
    std::cout << "Minimal polynomial degree: " << rur.minimal_polynomial.size() - 1 << "\n\n";
    
    // Step 2: Convert to numerical coefficients
    std::cout << "Step 2: Converting to numerical representation...\n";
    std::vector<double> num_coeffs = to_double_coeffs(rur.minimal_polynomial);
    
    std::cout << "Minimal polynomial f(T) = ";
    for (int i = num_coeffs.size() - 1; i >= 0; --i) {
        if (i < num_coeffs.size() - 1 && num_coeffs[i] >= 0) std::cout << "+";
        std::cout << num_coeffs[i];
        if (i > 0) std::cout << "*T^" << i << " ";
    }
    std::cout << "\n\n";
    
    // Step 3: Find roots using eigenvalue method
    std::cout << "Step 3: Finding roots of univariate polynomial...\n";
    
    Eigen::MatrixXd companion = build_companion_matrix(num_coeffs);
    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    
    auto eigenvalues = solver.eigenvalues();
    
    std::cout << "Found " << eigenvalues.size() << " roots:\n";
    for (int i = 0; i < eigenvalues.size(); ++i) {
        std::cout << "  T[" << i << "] = " << eigenvalues[i] << "\n";
    }
    
    // Step 4: Back-substitution to get multivariate solutions
    std::cout << "\nStep 4: Computing multivariate solutions...\n";
    std::cout << "Using parameterization: x_i = g_x(T_i)/f'(T_i), y_i = g_y(T_i)/f'(T_i)\n\n";
    
    // For each root, compute the solution
    std::vector<std::pair<std::complex<double>, std::complex<double>>> solutions;
    
    for (int i = 0; i < eigenvalues.size(); ++i) {
        std::complex<double> T = eigenvalues[i];
        
        // Evaluate f'(T)
        std::complex<double> fprime = evaluate_derivative(num_coeffs, T);
        
        // For this simple example, assuming parameterization is linear
        // In general, you would evaluate the numerator polynomials here
        // For now, we'll use a simplified approach
        
        std::complex<double> x = T / fprime;  // Simplified - actual would use g_x(T)
        std::complex<double> y = T / fprime;  // Simplified - actual would use g_y(T)
        
        solutions.push_back({x, y});
    }
    
    // Step 5: Display solutions
    std::cout << "Solutions to the system:\n";
    std::cout << std::fixed << std::setprecision(6);
    
    for (size_t i = 0; i < solutions.size(); ++i) {
        auto [x, y] = solutions[i];
        std::cout << "  Solution " << i+1 << ": ";
        
        // Check if essentially real
        if (std::abs(x.imag()) < 1e-10 && std::abs(y.imag()) < 1e-10) {
            std::cout << "(" << x.real() << ", " << y.real() << ")\n";
        } else {
            std::cout << "(" << x << ", " << y << ")\n";
        }
        
        // Verify the solution (optional)
        double x_real = x.real();
        double y_real = y.real();
        double eq1 = 9*x_real*x_real + 4*y_real*y_real - 36;
        double eq2 = 4*x_real*x_real + 9*y_real*y_real - 36;
        
        if (std::abs(x.imag()) < 1e-10 && std::abs(y.imag()) < 1e-10) {
            std::cout << "    Verification: eq1 = " << eq1 
                      << ", eq2 = " << eq2 << "\n";
        }
    }
    
    std::cout << "\nNote: This is a simplified example. In practice, you would:\n";
    std::cout << "1. Extract the actual numerator polynomials from RUR\n";
    std::cout << "2. Evaluate them at each root T_i\n";
    std::cout << "3. Use a more robust univariate solver for high-degree polynomials\n";
    std::cout << "4. Implement solution refinement for better accuracy\n";
    
    return 0;
}