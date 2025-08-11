#include <iostream>
#include <vector>
#include <string>
#include "julia_rur/polynomial_solver_enhanced.hpp"

int main() {
    // Circle equation (obvious positive-dimensional)
    std::vector<std::string> polynomials = { "x^2 + y^2 - 1" };
    std::vector<std::string> variables = { "x", "y" };
    
    julia_rur::EnhancedSolverConfig config;
    config.fail_on_positive_dimensional = true;
    config.auto_hyperplane_sections = false;
    config.verbose = true;
    config.enable_dimension_precheck = false;  // As per test fixture
    
    auto solution = julia_rur::solve_polynomial_system_enhanced(polynomials, variables, config);
    
    std::cout << "\nFinal Result:\n";
    std::cout << "Success: " << (solution.success ? "true" : "false") << "\n";
    std::cout << "Computed dimension: " << solution.computed_dimension << "\n";
    std::cout << "Error message: " << solution.error_message << "\n";
    
    return 0;
}