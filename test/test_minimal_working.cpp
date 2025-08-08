#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include "../src/julia_rur/polynomial_evaluator.hpp"
#include <iostream>
#include <iomanip>

using namespace julia_rur;

int main() {
    std::cout << "Testing minimal polynomial systems...\n\n";
    
    // Test 1: Simple univariate
    {
        std::cout << "Test 1: x^2 - 4 = 0\n";
        std::vector<std::string> polys = {"x^2 - 4"};
        std::vector<std::string> vars = {"x"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        
        auto result = solve_polynomial_system_enhanced(polys, vars, config);
        
        if (result.success) {
            std::cout << "  SUCCESS: Found " << result.solutions.size() << " solutions\n";
            for (size_t i = 0; i < result.solutions.size(); ++i) {
                std::cout << "    x = " << result.solutions[i][0] << "\n";
            }
        } else {
            std::cout << "  FAILED: " << result.error_message << "\n";
        }
    }
    
    // Test 2: Simple 2x2 linear
    {
        std::cout << "\nTest 2: Linear system\n";
        std::vector<std::string> polys = {"x + y - 3", "x - y - 1"};
        std::vector<std::string> vars = {"x", "y"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        
        auto result = solve_polynomial_system_enhanced(polys, vars, config);
        
        if (result.success) {
            std::cout << "  SUCCESS: Found " << result.solutions.size() << " solutions\n";
            for (size_t i = 0; i < result.solutions.size(); ++i) {
                std::cout << "    (x,y) = (" << result.solutions[i][0] 
                          << ", " << result.solutions[i][1] << ")\n";
            }
        } else {
            std::cout << "  FAILED: " << result.error_message << "\n";
        }
    }
    
    // Test 3: Circle and line
    {
        std::cout << "\nTest 3: Circle and line\n";
        std::vector<std::string> polys = {"x^2 + y^2 - 1", "y - x"};
        std::vector<std::string> vars = {"x", "y"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        
        auto result = solve_polynomial_system_enhanced(polys, vars, config);
        
        if (result.success) {
            std::cout << "  SUCCESS: Found " << result.solutions.size() << " solutions\n";
            for (size_t i = 0; i < result.solutions.size(); ++i) {
                auto x = result.solutions[i][0];
                auto y = result.solutions[i][1];
                std::cout << "    (x,y) = (" << x << ", " << y << ")\n";
                
                // Check residual
                double res = PolynomialEvaluator::compute_residual(polys, vars, result.solutions[i]);
                std::cout << "    Residual: " << std::scientific << res << std::fixed << "\n";
            }
        } else {
            std::cout << "  FAILED: " << result.error_message << "\n";
        }
    }
    
    return 0;
}