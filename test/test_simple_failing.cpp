#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>
#include <unistd.h>

using namespace julia_rur;

int main() {
    // Suppress stderr for cleaner output
    int saved_stderr = dup(2);
    
    EnhancedSolverConfig config;
    config.verbose = false;
    
    struct TestCase {
        const char* name;
        std::vector<std::string> polys;
        std::vector<std::string> vars;
    };
    
    // Test cases from simplest to most complex
    TestCase tests[] = {
        {"x-1", {"x - 1"}, {"x"}},
        {"x^2-4", {"x^2 - 4"}, {"x"}},
        {"2x2_linear", {"x + y - 3", "x - y - 1"}, {"x", "y"}},
        {"elem_sym_2", {"x + y - 2", "x*y - 1"}, {"x", "y"}},
        {"circle_line", {"x^2 + y^2 - 1", "y"}, {"x", "y"}},
        {"elem_sym_3", {"x + y + z - 3", "x*y + y*z + z*x - 3", "x*y*z - 1"}, {"x", "y", "z"}},
        {"cyclic_3", {"x + y + z", "x*y + y*z + z*x", "x*y*z - 1"}, {"x", "y", "z"}}
    };
    
    for (const auto& test : tests) {
        std::cout << "\nTesting: " << test.name << "\n";
        std::cout << "Polynomials: ";
        for (const auto& p : test.polys) std::cout << p << "; ";
        std::cout << "\nVariables: ";
        for (const auto& v : test.vars) std::cout << v << " ";
        std::cout << "\n";
        
        // Redirect stderr to suppress debug output
        freopen("/dev/null", "w", stderr);
        
        auto result = solve_polynomial_system_enhanced(test.polys, test.vars, config);
        
        // Restore stderr
        dup2(saved_stderr, 2);
        
        if (result.success) {
            std::cout << "✓ SUCCESS: " << result.solutions.size() << " solutions found\n";
            // Print first solution if exists
            if (!result.solutions.empty()) {
                std::cout << "  First solution: ";
                for (size_t i = 0; i < result.solutions[0].size(); ++i) {
                    auto val = result.solutions[0][i];
                    std::cout << test.vars[i] << "=" << val.real();
                    if (std::abs(val.imag()) > 1e-10) {
                        std::cout << "+" << val.imag() << "i";
                    }
                    std::cout << " ";
                }
                std::cout << "\n";
            }
        } else {
            std::cout << "✗ FAILED: " << result.error_message << "\n";
            // This is our first failure - let's debug it
            if (test.name == std::string("elem_sym_2")) {
                std::cout << "\n=== DEBUGGING FIRST FAILURE ===\n";
                config.verbose = true;
                auto debug_result = solve_polynomial_system_enhanced(test.polys, test.vars, config);
                break;
            }
        }
    }
    
    close(saved_stderr);
    return 0;
}