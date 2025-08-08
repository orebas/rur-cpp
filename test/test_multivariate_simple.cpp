#include "../src/julia_rur/polynomial_solver.hpp"
#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <cmath>

using namespace julia_rur;

int main() {
    std::cout << "=== Testing simple multivariate system ===" << std::endl;
    std::cout << "System: x^2 + y^2 - 1 = 0 (circle), x - y = 0 (line)" << std::endl;
    std::cout << "Expected solutions: (1/√2, 1/√2) and (-1/√2, -1/√2)\n" << std::endl;
    
    std::vector<std::string> polynomials = {
        "x^2 + y^2 - 1",  // Circle
        "x - y"           // Line x = y
    };
    std::vector<std::string> variables = {"x", "y"};
    
    RURConfig config;
    config.verbose = true;  // Enable verbose to see what's happening
    
    std::cout << "Computing RUR..." << std::endl;
    RationalRURResult rur = compute_rational_rur(polynomials, variables, config);
    
    if (!rur.success) {
        std::cout << "RUR computation failed: " << rur.error_message << std::endl;
        return 1;
    }
    
    std::cout << "\nRUR Results:" << std::endl;
    std::cout << "- Minimal polynomial coefficients: ";
    for (const auto& c : rur.minimal_polynomial) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    std::cout << "- Number of parameterizations: " << rur.numerators.size() << std::endl;
    for (size_t i = 0; i < rur.numerators.size(); ++i) {
        std::cout << "- Parameterization " << i << " size: " << rur.numerators[i].size() << std::endl;
        if (!rur.numerators[i].empty()) {
            std::cout << "  Coefficients: ";
            for (const auto& c : rur.numerators[i]) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }
    }
    
    return 0;
}