#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    // Test a simple 3x3 linear system
    std::vector<std::string> polys = {
        "x + 2*y + 3*z - 6",
        "2*x - y + z - 1", 
        "3*x + y - z - 2"
    };
    std::vector<std::string> vars = {"x", "y", "z"};
    
    std::cout << "Testing 3x3 linear system:\n";
    for (const auto& p : polys) {
        std::cout << "  " << p << "\n";
    }
    std::cout << "\nThis should have exactly 1 solution.\n\n";
    
    EnhancedSolverConfig config;
    config.verbose = false;
    
    auto result = solve_polynomial_system_enhanced(polys, vars, config);
    
    if (result.success) {
        std::cout << "✓ SUCCESS: Found " << result.solutions.size() << " solutions\n";
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
        std::cout << "✗ FAILED: " << result.error_message << "\n";
    }
    
    return result.success ? 0 : 1;
}