#include <iostream>
#include <cmath>
#include "julia_rur/rur_main_algorithm.hpp"
#include "julia_rur/polynomial_solver.hpp"
#include "julia_rur/numerical_roots_eigen.hpp"

using namespace julia_rur;

int main() {
    std::cout << "Debug: x^2 - 2 = 0\n";
    std::cout << "Expected roots: ±√2 ≈ ±1.41421356\n\n";
    
    // Test the polynomial system
    std::vector<std::string> polynomials = {"x^2 - 2"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    std::cout << "=== Computing RUR ===\n";
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (!rur_result.success) {
        std::cerr << "RUR computation failed\n";
        return 1;
    }
    
    std::cout << "Minimal polynomial coefficients: ";
    for (const auto& c : rur_result.minimal_polynomial) {
        std::cout << c << " ";
    }
    std::cout << "\n";
    
    // Interpret the minimal polynomial
    std::cout << "\nMinimal polynomial: ";
    bool first = true;
    for (int i = rur_result.minimal_polynomial.size() - 1; i >= 0; --i) {
        const auto& coeff = rur_result.minimal_polynomial[i];
        if (coeff != 0) {
            if (!first && coeff > 0) std::cout << " + ";
            if (coeff < 0) std::cout << " - ";
            
            mpq_class abs_coeff = abs(coeff);
            if (abs_coeff != 1 || i == 0) {
                std::cout << abs_coeff;
            }
            
            if (i > 0) {
                if (abs_coeff != 1 || i == 0) std::cout << "*";
                std::cout << "T";
                if (i > 1) std::cout << "^" << i;
            }
            first = false;
        }
    }
    std::cout << " = 0\n";
    
    // Now find roots using our polynomial root finder
    std::cout << "\n=== Finding roots of minimal polynomial ===\n";
    std::vector<std::complex<double>> roots = find_polynomial_roots(rur_result.minimal_polynomial);
    
    std::cout << "Found " << roots.size() << " roots:\n";
    for (size_t i = 0; i < roots.size(); ++i) {
        std::cout << "  Root " << i << ": " << roots[i].real();
        if (std::abs(roots[i].imag()) > 1e-10) {
            std::cout << " + " << roots[i].imag() << "i";
        }
        std::cout << "\n";
    }
    
    // Let's also manually check what Eigen is doing
    std::cout << "\n=== Manual check with Eigen ===\n";
    std::vector<mpq_class> test_poly = rur_result.minimal_polynomial;
    
    // Convert to double coefficients
    std::vector<double> double_coeffs;
    std::cout << "Double coefficients: ";
    for (const auto& c : test_poly) {
        double val = c.get_d();
        double_coeffs.push_back(val);
        std::cout << val << " ";
    }
    std::cout << "\n";
    
    // Create companion matrix manually to debug
    int n = double_coeffs.size() - 1;
    std::cout << "\nDegree of polynomial: " << n << "\n";
    
    if (n == 2) {
        // For quadratic aT^2 + bT + c = 0, roots are (-b ± √(b²-4ac)) / 2a
        double a = double_coeffs[2];
        double b = double_coeffs[1];
        double c = double_coeffs[0];
        
        std::cout << "Quadratic formula: a=" << a << ", b=" << b << ", c=" << c << "\n";
        double discriminant = b*b - 4*a*c;
        std::cout << "Discriminant: b²-4ac = " << discriminant << "\n";
        
        if (discriminant >= 0) {
            double root1 = (-b + std::sqrt(discriminant)) / (2*a);
            double root2 = (-b - std::sqrt(discriminant)) / (2*a);
            std::cout << "Real roots: " << root1 << ", " << root2 << "\n";
        } else {
            std::cout << "Complex roots with discriminant " << discriminant << "\n";
        }
    }
    
    // Now test the full solver
    std::cout << "\n=== Full polynomial solver ===\n";
    PolynomialSystemSolution solution = solve_polynomial_system_complete(
        polynomials, variables, config
    );
    
    if (solution.success) {
        std::cout << "Solutions:\n";
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            std::cout << "  " << i << ". x = " << solution.solutions[i][0].real();
            if (std::abs(solution.solutions[i][0].imag()) > 1e-10) {
                std::cout << " + " << solution.solutions[i][0].imag() << "i";
            }
            std::cout << "\n";
        }
    }
    
    return 0;
}