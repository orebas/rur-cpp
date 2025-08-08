#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    std::cout << "Debugging polynomial solver for x^2 - 2..." << std::endl;
    
    std::vector<std::string> polynomials = {"x^2 - 2"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    // First get the RUR result
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (rur_result.success) {
        std::cout << "\nRUR result:" << std::endl;
        std::cout << "Minimal polynomial coefficients: ";
        for (const auto& c : rur_result.minimal_polynomial) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        
        // Now solve
        auto solution = solve_polynomial_system_complete(polynomials, variables, config);
        
        std::cout << "\nSolution object:" << std::endl;
        std::cout << "Minimal polynomial in solution: ";
        for (const auto& c : solution.minimal_polynomial) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        
        // Print using print_solution
        std::cout << "\nOutput from print_solution:" << std::endl;
        print_solution(solution, std::cout);
    }
    
    return 0;
}