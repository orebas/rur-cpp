#include "../src/julia_rur/numerical_roots_eigen.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace julia_rur;

void test_polynomial_roots() {
    std::cout << "Testing polynomial root finding with x^2 - 2" << std::endl;
    
    // Coefficients of x^2 - 2: [-2, 0, 1]
    std::vector<mpq_class> coeffs = {
        mpq_class(-2),
        mpq_class(0),
        mpq_class(1)
    };
    
    std::cout << "Polynomial coefficients (low to high degree): ";
    for (const auto& c : coeffs) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    auto roots = find_polynomial_roots(coeffs);
    
    std::cout << "Found " << roots.size() << " roots:" << std::endl;
    for (const auto& root : roots) {
        std::cout << "  " << std::setprecision(10) 
                  << root.real();
        if (std::abs(root.imag()) > 1e-10) {
            std::cout << " + " << root.imag() << "i";
        }
        std::cout << std::endl;
    }
    
    // Check accuracy
    double expected = std::sqrt(2.0);
    bool found_positive = false;
    bool found_negative = false;
    
    for (const auto& root : roots) {
        if (std::abs(root.real() - expected) < 1e-10 && std::abs(root.imag()) < 1e-10) {
            found_positive = true;
        }
        if (std::abs(root.real() + expected) < 1e-10 && std::abs(root.imag()) < 1e-10) {
            found_negative = true;
        }
    }
    
    if (found_positive && found_negative) {
        std::cout << "✓ Correctly found ±√2" << std::endl;
    } else {
        std::cout << "✗ Failed to find correct roots" << std::endl;
    }
}

void test_cubic_roots() {
    std::cout << "\nTesting cubic polynomial: x^3 - 2x - 5" << std::endl;
    
    // Coefficients: [-5, -2, 0, 1]
    std::vector<mpq_class> coeffs = {
        mpq_class(-5),
        mpq_class(-2),
        mpq_class(0),
        mpq_class(1)
    };
    
    // Test with known root: x ≈ 2.094551482
    double known_root = 2.094551482;
    std::cout << "Known real root: " << known_root << std::endl;
    std::cout << "Verification: " << known_root << "^3 - 2*" << known_root << " - 5 = " 
              << (known_root * known_root * known_root - 2 * known_root - 5) << std::endl;
    
    auto roots = find_polynomial_roots(coeffs);
    
    std::cout << "Found " << roots.size() << " roots:" << std::endl;
    for (const auto& root : roots) {
        std::cout << "  " << std::setprecision(10) 
                  << root.real();
        if (std::abs(root.imag()) > 1e-10) {
            std::cout << " + " << root.imag() << "i";
        }
        
        // Verify by substitution
        std::complex<double> value = root * root * root - 2.0 * root - 5.0;
        std::cout << " (residual: " << std::abs(value) << ")";
        std::cout << std::endl;
    }
}

void test_full_system() {
    std::cout << "\nTesting full system solver for x^2 - 2" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x^2-2"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    // First compute RUR to see what polynomial we get
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);
    if (rur_result.success) {
        std::cout << "RUR minimal polynomial coefficients: ";
        for (const auto& c : rur_result.minimal_polynomial) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }
    
    auto solution = solve_polynomial_system(polynomials, variables, config);
    
    if (solution.success) {
        std::cout << "Found " << solution.solutions.size() << " solutions:" << std::endl;
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            std::cout << "  Solution " << (i+1) << ": ";
            for (size_t j = 0; j < variables.size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << variables[j] << " = " << solution.solutions[i][j].real();
                if (std::abs(solution.solutions[i][j].imag()) > 1e-10) {
                    std::cout << " + " << solution.solutions[i][j].imag() << "i";
                }
            }
            if (solution.is_real[i]) {
                std::cout << " (real)";
            }
            std::cout << std::endl;
        }
    } else {
        std::cout << "Failed: " << solution.error_message << std::endl;
    }
}

int main() {
    test_polynomial_roots();
    test_cubic_roots();
    test_full_system();
    return 0;
}