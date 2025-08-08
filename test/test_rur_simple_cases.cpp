#include "../src/julia_rur/rur_main_algorithm.hpp"
#include "../src/julia_rur/numerical_roots_eigen.hpp"
#include <iostream>
#include <iomanip>

using namespace julia_rur;

void test_polynomial(const std::string& name, 
                    const std::vector<std::string>& polynomials,
                    const std::vector<std::string>& variables) {
    std::cout << "\n=== " << name << " ===" << std::endl;
    
    RURConfig config;
    config.verbose = false;
    
    auto rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (!rur_result.success) {
        std::cout << "RUR failed: " << rur_result.error_message << std::endl;
        return;
    }
    
    // Print minimal polynomial
    std::cout << "Minimal polynomial: ";
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
                std::cout << "T";
                if (i > 1) std::cout << "^" << i;
            }
            first = false;
        }
    }
    std::cout << " = 0" << std::endl;
    
    // Print parameterizations
    std::cout << "Parameterizations:" << std::endl;
    for (size_t i = 0; i < variables.size(); ++i) {
        std::cout << "  " << variables[i] << " = ";
        if (i < rur_result.numerators.size() && !rur_result.numerators[i].empty()) {
            // Print numerator polynomial
            bool first_term = true;
            for (int j = rur_result.numerators[i].size() - 1; j >= 0; --j) {
                const auto& coeff = rur_result.numerators[i][j];
                if (coeff != 0) {
                    if (!first_term && coeff > 0) std::cout << " + ";
                    if (coeff < 0) std::cout << " - ";
                    
                    mpq_class abs_coeff = abs(coeff);
                    if (abs_coeff != 1 || j == 0) {
                        std::cout << abs_coeff;
                    }
                    
                    if (j > 0) {
                        std::cout << "T";
                        if (j > 1) std::cout << "^" << j;
                    }
                    first_term = false;
                }
            }
            std::cout << " / f'(T)";
        } else {
            std::cout << "T";  // Direct parameterization
        }
        std::cout << std::endl;
    }
    
    // Find and print roots
    auto roots = find_polynomial_roots(rur_result.minimal_polynomial);
    std::cout << "Roots of minimal polynomial:" << std::endl;
    for (const auto& root : roots) {
        std::cout << "  T = " << std::setprecision(10) << root << std::endl;
    }
}

int main() {
    // Test univariate cases
    test_polynomial("x^2 - 2", {"1*x^2-2"}, {"x"});
    test_polynomial("x - 3", {"1*x-3"}, {"x"});
    test_polynomial("x^3 - 2x - 5", {"1*x^3-2*x-5"}, {"x"});
    
    // Test bivariate case
    test_polynomial("Circle and line", 
                   {"1*x^2+1*y^2-1", "1*x-1*y"}, 
                   {"x", "y"});
    
    return 0;
}