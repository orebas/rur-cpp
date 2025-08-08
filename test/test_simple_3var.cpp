#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    // Try a really simple 3-variable system that should definitely work
    std::vector<std::string> polys = {
        "x^2 - 1",
        "y^2 - 2", 
        "z^2 - 3"
    };
    std::vector<std::string> vars = {"x", "y", "z"};
    
    std::cout << "Testing simple 3-variable system:\n";
    for (const auto& p : polys) {
        std::cout << "  " << p << "\n";
    }
    std::cout << "\nThis should have 8 solutions (2^3).\n\n";
    
    EnhancedSolverConfig config;
    config.verbose = true;
    
    auto result = solve_polynomial_system_enhanced(polys, vars, config);
    
    if (result.success) {
        std::cout << "\n✓ SUCCESS: Found " << result.solutions.size() << " solutions\n";
        for (size_t i = 0; i < result.solutions.size(); ++i) {
            std::cout << "Solution " << (i+1) << ": (";
            for (size_t j = 0; j < vars.size(); ++j) {
                if (j > 0) std::cout << ", ";
                auto val = result.solutions[i][j];
                if (std::abs(val.imag()) < 1e-10) {
                    std::cout << val.real();
                } else {
                    std::cout << val.real() << "+" << val.imag() << "i";
                }
            }
            std::cout << ")\n";
        }
    } else {
        std::cout << "\n✗ FAILED: " << result.error_message << "\n";
    }
    
    return result.success ? 0 : 1;
}