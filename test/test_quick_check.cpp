#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>
#include <unistd.h>

using namespace julia_rur;

int main() {
    // Suppress all stderr output
    int saved_stderr = dup(2);
    freopen("/dev/null", "w", stderr);
    
    EnhancedSolverConfig config;
    config.verbose = false;
    
    struct Test {
        const char* name;
        std::vector<std::string> polys;
        std::vector<std::string> vars;
    };
    
    Test tests[] = {
        {"x-1", {"x - 1"}, {"x"}},
        {"x^2-4", {"x^2 - 4"}, {"x"}},
        {"2x2 linear", {"x + y - 3", "x - y - 1"}, {"x", "y"}},
        {"circle-line", {"x^2 + y^2 - 1", "y"}, {"x", "y"}},
        {"elem-sym-2", {"x + y - 2", "x*y - 1"}, {"x", "y"}},
        {"elem-sym-3", {"x + y + z - 3", "x*y + y*z + z*x - 3", "x*y*z - 1"}, {"x", "y", "z"}},
        {"cyclic-3", {"x + y + z", "x*y + y*z + z*x", "x*y*z - 1"}, {"x", "y", "z"}},
        {"katsura-3", {"x + 2*y + 2*z - 1", "x^2 + 2*y^2 + 2*z^2 - x", "2*x*y + 2*y*z - y"}, {"x", "y", "z"}},
    };
    
    dup2(saved_stderr, 2);
    close(saved_stderr);
    
    std::cout << "Quick Test Results:\n";
    std::cout << "==================\n";
    
    for (const auto& test : tests) {
        freopen("/dev/null", "w", stderr);
        auto result = solve_polynomial_system_enhanced(test.polys, test.vars, config);
        
        std::cout << test.name << ": ";
        if (result.success) {
            std::cout << "PASS (" << result.solutions.size() << " solutions)\n";
        } else {
            std::cout << "FAIL - " << result.error_message << "\n";
        }
    }
    
    return 0;
}