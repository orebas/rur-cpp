#include <iostream>
#include "julia_rur/rur_main_algorithm.hpp"
#include "julia_rur/polynomial_solver.hpp"

using namespace julia_rur;

int main() {
    std::cout << "Testing simple linear system x - 3 = 0\n\n";
    
    // Test system: x - 3 = 0
    std::vector<std::string> polynomials = {"x - 3"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    config.initial_prime_bits = 20;  // Use smaller primes to avoid F4 issues
    
    std::cout << "=== Computing RUR with smaller primes ===\n";
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (!rur_result.success) {
        std::cerr << "RUR computation failed: " << rur_result.error_message << std::endl;
        return 1;
    }
    
    std::cout << "RUR succeeded!\n";
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
    
    // Try to solve the system
    std::cout << "\n=== Solving polynomial system ===\n";
    PolynomialSystemSolution solution = solve_polynomial_system_complete(
        polynomials, variables, config
    );
    
    if (!solution.success) {
        std::cerr << "Failed to solve: " << solution.error_message << std::endl;
        return 1;
    }
    
    std::cout << "Solutions:\n";
    for (size_t i = 0; i < solution.solutions.size(); ++i) {
        std::cout << "  " << (i+1) << ". x = " << solution.solutions[i][0].real();
        if (std::abs(solution.solutions[i][0].imag()) > 1e-10) {
            std::cout << " + " << solution.solutions[i][0].imag() << "i";
        }
        std::cout << "\n";
    }
    
    return 0;
}