#include <iostream>
#include "julia_rur/rur_main_algorithm.hpp"
#include "julia_rur/polynomial_solver.hpp"

using namespace julia_rur;

int main() {
    std::cout << "Debug: Linear system x - 3 = 0 (full solver)\n\n";
    
    // Test system: x - 3 = 0
    std::vector<std::string> polynomials = {"x - 3"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    config.initial_prime_bits = 31;
    
    std::cout << "=== Computing RUR ===\n";
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (!rur_result.success) {
        std::cerr << "RUR computation failed: " << rur_result.error_message << std::endl;
        return 1;
    }
    
    std::cout << "\n=== RUR Result ===\n";
    std::cout << "Minimal polynomial: ";
    for (const auto& coeff : rur_result.minimal_polynomial) {
        std::cout << coeff << " ";
    }
    std::cout << "\n";
    
    std::cout << "Number of numerators: " << rur_result.numerators.size() << "\n";
    for (size_t i = 0; i < rur_result.numerators.size(); ++i) {
        std::cout << "Numerator " << i << ": ";
        if (rur_result.numerators[i].empty()) {
            std::cout << "[empty]";
        } else {
            for (const auto& coeff : rur_result.numerators[i]) {
                std::cout << coeff << " ";
            }
        }
        std::cout << "\n";
    }
    
    std::cout << "\n=== Solving system ===\n";
    PolynomialSystemSolution solution = solve_polynomial_system_complete(
        polynomials, variables, config
    );
    
    print_solution(solution, std::cout);
    
    return 0;
}