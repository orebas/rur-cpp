#include "../src/julia_rur/rur_main_algorithm.hpp"
#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    std::cout << "Testing RUR for x - 3 = 0" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x-3"};
    std::vector<std::string> variables = {"x"};
    RURConfig config;
    config.verbose = false;
    
    auto rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (rur_result.success) {
        std::cout << "\nRUR Result:" << std::endl;
        std::cout << "  Quotient basis size: " << rur_result.quotient_basis.size() << std::endl;
        std::cout << "  Minimal polynomial: ";
        for (auto c : rur_result.minimal_polynomial) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        
        std::cout << "  Number of numerators: " << rur_result.numerators.size() << std::endl;
        for (size_t i = 0; i < rur_result.numerators.size(); ++i) {
            std::cout << "    Numerator[" << i << "] size: " << rur_result.numerators[i].size() << std::endl;
            if (!rur_result.numerators[i].empty()) {
                std::cout << "      Coefficients: ";
                for (auto c : rur_result.numerators[i]) {
                    std::cout << c << " ";
                }
                std::cout << std::endl;
            }
        }
        
        // Now solve
        auto solution = solve_polynomial_system_complete(polynomials, variables, config);
        if (solution.success) {
            std::cout << "\nSolution:" << std::endl;
            for (size_t i = 0; i < solution.solutions.size(); ++i) {
                std::cout << "  x = " << solution.solutions[i][0] << std::endl;
            }
        }
    }
    
    return 0;
}