#include "../src/julia_rur/rur_main_algorithm.hpp"
#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    std::cout << "=== Debugging RUR for x^3 - 2x - 5 ===\n" << std::endl;
    
    std::vector<std::string> polynomials = {"x^3 - 2*x - 5"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    RationalRURResult rur = compute_rational_rur(polynomials, variables, config);
    
    std::cout << "RUR Results:" << std::endl;
    std::cout << "- Success: " << (rur.success ? "yes" : "no") << std::endl;
    std::cout << "- Minimal polynomial coefficients: ";
    for (const auto& c : rur.minimal_polynomial) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    std::cout << "- Number of numerators: " << rur.numerators.size() << std::endl;
    if (!rur.numerators.empty()) {
        std::cout << "- First numerator size: " << rur.numerators[0].size() << std::endl;
        if (!rur.numerators[0].empty()) {
            std::cout << "- First numerator coefficients: ";
            for (const auto& c : rur.numerators[0]) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "\nFor x^3 - 2x - 5 = 0 where x is the separating element:" << std::endl;
    std::cout << "- Expected: numerators should be empty or represent x = T" << std::endl;
    std::cout << "- Actual: " << (rur.numerators.empty() ? "empty" : "not empty") << std::endl;
    
    return 0;
}