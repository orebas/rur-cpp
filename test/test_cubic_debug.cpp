#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>
#include <cmath>

using namespace julia_rur;

int main() {
    std::cout << "=== Debugging x^3 - 2x - 5 ===\n" << std::endl;
    
    std::vector<std::string> polynomials = {"x^3 - 2*x - 5"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    
    if (solution.success) {
        std::cout << "Minimal polynomial: ";
        for (const auto& c : solution.minimal_polynomial) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        
        std::cout << "\nRoots:" << std::endl;
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            double x = solution.solutions[i][0].real();
            double imag = solution.solutions[i][0].imag();
            std::cout << "  x = " << x;
            if (std::abs(imag) > 1e-10) {
                std::cout << " + " << imag << "i";
            }
            
            // Verify by substituting back
            double value = x*x*x - 2*x - 5;
            std::cout << "  (check: x^3 - 2x - 5 = " << value << ")";
            std::cout << std::endl;
        }
        
        std::cout << "\nExpected real root: 2.094551482" << std::endl;
        std::cout << "Check: 2.094551482^3 - 2*2.094551482 - 5 = " 
                  << (2.094551482*2.094551482*2.094551482 - 2*2.094551482 - 5) << std::endl;
    }
    
    return 0;
}