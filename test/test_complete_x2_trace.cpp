#include "../src/julia_rur/polynomial_solver.hpp"
#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <iomanip>

using namespace julia_rur;

int main() {
    std::cout << "=== Complete trace for x^2 - 2 ===\n" << std::endl;
    
    std::vector<std::string> polynomials = {"x^2 - 2"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    // Step 1: Compute RUR
    std::cout << "1. Computing RUR..." << std::endl;
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (!rur_result.success) {
        std::cout << "RUR failed: " << rur_result.error_message << std::endl;
        return 1;
    }
    
    std::cout << "   - Minimal polynomial coefficients: ";
    for (const auto& c : rur_result.minimal_polynomial) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    // Print the polynomial in human-readable form
    std::cout << "   - Minimal polynomial: ";
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
    
    // Step 2: Call polynomial solver
    std::cout << "\n2. Calling polynomial solver..." << std::endl;
    PolynomialSystemSolution solution = solve_polynomial_system_complete(polynomials, variables, config);
    
    std::cout << "   - Solution minimal polynomial coefficients: ";
    for (const auto& c : solution.minimal_polynomial) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    // Step 3: Check what print_solution outputs
    std::cout << "\n3. Output from print_solution:" << std::endl;
    print_solution(solution, std::cout);
    
    // Step 4: Direct coefficient check
    std::cout << "\n4. Direct coefficient analysis:" << std::endl;
    if (solution.minimal_polynomial.size() >= 3) {
        std::cout << "   - Coefficient at index 0 (constant): " << solution.minimal_polynomial[0] << std::endl;
        std::cout << "   - Coefficient at index 1 (linear): " << solution.minimal_polynomial[1] << std::endl;
        std::cout << "   - Coefficient at index 2 (quadratic): " << solution.minimal_polynomial[2] << std::endl;
    }
    
    // Step 5: Manual polynomial printing
    std::cout << "\n5. Manual polynomial construction:" << std::endl;
    std::cout << "   - From coefficients [-2, 0, 1]: ";
    if (solution.minimal_polynomial.size() >= 3 && 
        solution.minimal_polynomial[0] == -2 && 
        solution.minimal_polynomial[1] == 0 && 
        solution.minimal_polynomial[2] == 1) {
        std::cout << "T^2 - 2 = 0 (CORRECT)" << std::endl;
    } else {
        std::cout << "Coefficients don't match expected [-2, 0, 1]" << std::endl;
    }
    
    return 0;
}