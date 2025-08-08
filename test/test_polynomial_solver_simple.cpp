#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    std::cout << "Testing simple polynomial solver" << std::endl;
    
    // Test 1: x^2 - 2
    {
        std::cout << "\nTest: x^2 - 2" << std::endl;
        std::vector<std::string> polynomials = {"1*x^2-2"};
        std::vector<std::string> variables = {"x"};
        
        RURConfig config;
        config.verbose = false;
        
        auto solution = solve_polynomial_system_complete(polynomials, variables, config);
        
        if (solution.success) {
            std::cout << "Success! Found " << solution.solutions.size() << " solutions" << std::endl;
            for (size_t i = 0; i < solution.solutions.size(); ++i) {
                std::cout << "  x = " << solution.solutions[i][0] << std::endl;
            }
        } else {
            std::cout << "Failed: " << solution.error_message << std::endl;
        }
    }
    
    std::cout << "\nTest completed without errors" << std::endl;
    return 0;
}