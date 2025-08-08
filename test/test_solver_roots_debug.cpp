#include "../src/julia_rur/polynomial_solver.hpp"
#include "../src/julia_rur/numerical_roots_eigen.hpp"
#include <iostream>
#include <cmath>

using namespace julia_rur;

int main() {
    std::cout << "=== Debugging polynomial solver roots for x^2 - 2 ===\n" << std::endl;
    
    std::vector<std::string> polynomials = {"x^2 - 2"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    // Get RUR result
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (!rur_result.success) {
        std::cout << "RUR failed!" << std::endl;
        return 1;
    }
    
    std::cout << "RUR Results:" << std::endl;
    std::cout << "- Minimal polynomial coefficients: ";
    for (const auto& c : rur_result.minimal_polynomial) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    std::cout << "- Quotient basis size: " << rur_result.quotient_basis.size() << std::endl;
    std::cout << "- Number of numerators: " << rur_result.numerators.size() << std::endl;
    
    if (!rur_result.numerators.empty() && !rur_result.numerators[0].empty()) {
        std::cout << "- First numerator: ";
        for (const auto& c : rur_result.numerators[0]) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }
    
    // Find roots of minimal polynomial directly
    std::cout << "\nDirect root finding of minimal polynomial:" << std::endl;
    auto roots = find_polynomial_roots(rur_result.minimal_polynomial);
    for (const auto& root : roots) {
        std::cout << "  T = " << root.real();
        if (std::abs(root.imag()) > 1e-10) {
            std::cout << " + " << root.imag() << "i";
        }
        std::cout << std::endl;
    }
    
    // Now test the complete solver
    std::cout << "\nComplete polynomial solver:" << std::endl;
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    
    if (solution.success) {
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            std::cout << "Solution " << i+1 << ": x = " << solution.solutions[i][0].real();
            if (std::abs(solution.solutions[i][0].imag()) > 1e-10) {
                std::cout << " + " << solution.solutions[i][0].imag() << "i";
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "\nExpected: x = ±√2 ≈ ±" << sqrt(2.0) << std::endl;
    
    return 0;
}