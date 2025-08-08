#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    // Test the elem_sym_3 system that's failing
    std::vector<std::string> polys = {
        "x + y + z - 3",
        "x*y + y*z + z*x - 3", 
        "x*y*z - 1"
    };
    std::vector<std::string> vars = {"x", "y", "z"};
    
    std::cout << "Testing elementary symmetric polynomial system of degree 3\n";
    std::cout << "Expected: 6 solutions\n";
    std::cout << "================================================\n\n";
    
    EnhancedSolverConfig config;
    config.verbose = true;  // Enable verbose output to see what's happening
    
    auto result = solve_polynomial_system_enhanced(polys, vars, config);
    
    std::cout << "\n================================================\n";
    if (result.success) {
        std::cout << "SUCCESS: Found " << result.solutions.size() << " solutions\n";
        for (size_t i = 0; i < result.solutions.size(); ++i) {
            std::cout << "Solution " << (i+1) << ": ";
            for (size_t j = 0; j < vars.size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << vars[j] << " = " << result.solutions[i][j];
            }
            std::cout << "\n";
        }
    } else {
        std::cout << "FAILED: " << result.error_message << "\n";
        
        // Let's try a simpler strategy - increase the number of attempts
        std::cout << "\nTrying with more attempts...\n";
        // We'll need to modify the code to increase attempts
    }
    
    return result.success ? 0 : 1;
}