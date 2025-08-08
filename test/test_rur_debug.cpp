#include "../src/julia_rur/rur_main_algorithm.hpp"
#include "../src/julia_rur/numerical_roots_eigen.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    // Test 1: x^2 - 2
    {
        std::cout << "=== Test: x^2 - 2 ===" << std::endl;
        std::vector<std::string> polynomials = {"1*x^2-2"};
        std::vector<std::string> variables = {"x"};
        
        RURConfig config;
        config.verbose = true;
        
        auto rur_result = compute_rational_rur(polynomials, variables, config);
        
        if (rur_result.success) {
            std::cout << "\nRUR Result:" << std::endl;
            std::cout << "Minimal polynomial coefficients: ";
            for (const auto& c : rur_result.minimal_polynomial) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
            
            std::cout << "Numerators size: " << rur_result.numerators.size() << std::endl;
            if (!rur_result.numerators.empty() && !rur_result.numerators[0].empty()) {
                std::cout << "Numerator[0] coefficients: ";
                for (const auto& c : rur_result.numerators[0]) {
                    std::cout << c << " ";
                }
                std::cout << std::endl;
            }
            
            // Find roots
            auto roots = find_polynomial_roots(rur_result.minimal_polynomial);
            std::cout << "\nRoots of minimal polynomial:" << std::endl;
            for (const auto& r : roots) {
                std::cout << "  " << r << std::endl;
            }
        }
    }
    
    // Test 2: x - 3
    {
        std::cout << "\n=== Test: x - 3 ===" << std::endl;
        std::vector<std::string> polynomials = {"1*x-3"};
        std::vector<std::string> variables = {"x"};
        
        RURConfig config;
        config.verbose = false;
        
        auto rur_result = compute_rational_rur(polynomials, variables, config);
        
        if (rur_result.success) {
            std::cout << "\nRUR Result:" << std::endl;
            std::cout << "Minimal polynomial coefficients: ";
            for (const auto& c : rur_result.minimal_polynomial) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
            
            // Find roots
            auto roots = find_polynomial_roots(rur_result.minimal_polynomial);
            std::cout << "\nRoots of minimal polynomial:" << std::endl;
            for (const auto& r : roots) {
                std::cout << "  " << r << std::endl;
            }
        }
    }
    
    return 0;
}