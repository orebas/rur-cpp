#include <iostream>
#include <vector>
#include <string>
#include "julia_rur/polynomial_solver_enhanced.hpp"

int main() {
    // Cyclic-5 system
    std::vector<std::string> polynomials = {
        "x0+x1+x2+x3+x4",
        "x0*x1+x1*x2+x2*x3+x3*x4+x4*x0",
        "x0*x1*x2+x1*x2*x3+x2*x3*x4+x3*x4*x0+x4*x0*x1",
        "x0*x1*x2*x3+x1*x2*x3*x4+x2*x3*x4*x0+x3*x4*x0*x1+x4*x0*x1*x2",
        "x0*x1*x2*x3*x4-1"
    };
    std::vector<std::string> variables = {"x0", "x1", "x2", "x3", "x4"};
    
    std::cout << "Testing Cyclic-5 with dimension precheck disabled...\n";
    
    // Configure to skip dimension precheck
    julia_rur::EnhancedSolverConfig config;
    config.enable_dimension_precheck = false;  // Skip the problematic dimension check
    config.auto_hyperplane_sections = false;   // Don't try hyperplane sections
    config.fail_on_positive_dimensional = false;  // Don't fail even if it seems positive-dim
    config.verbose = false;
    
    auto solution = julia_rur::solve_polynomial_system_enhanced(polynomials, variables, config);
    
    if (solution.success) {
        std::cout << "SUCCESS! Found " << solution.solutions.size() << " solutions\n";
        std::cout << "Expected: 70 solutions for Cyclic-5\n";
    } else {
        std::cout << "FAILED: " << solution.error_message << "\n";
    }
    
    return solution.success ? 0 : 1;
}