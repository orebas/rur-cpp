#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>
#include <iomanip>
#include <fcntl.h>
#include <unistd.h>

using namespace julia_rur;

int main() {
    // Save original stdout/stderr
    int saved_stdout = dup(1);
    int saved_stderr = dup(2);
    
    EnhancedSolverConfig config;
    config.verbose = false;
    
    struct TestCase {
        const char* name;
        std::vector<std::string> polys;
        std::vector<std::string> vars;
        int expected_solutions;
    };
    
    // Test cases from simplest to most complex
    TestCase tests[] = {
        {"x-1", {"x - 1"}, {"x"}, 1},
        {"x^2-4", {"x^2 - 4"}, {"x"}, 2},
        {"linear_2x2", {"x + y - 3", "x - y - 1"}, {"x", "y"}, 1},
        {"elem_sym_2", {"x + y - 2", "x*y - 1"}, {"x", "y"}, 2},
        {"circle_line", {"x^2 + y^2 - 1", "y"}, {"x", "y"}, 2},
        {"elem_sym_3", {"x + y + z - 3", "x*y + y*z + z*x - 3", "x*y*z - 1"}, {"x", "y", "z"}, 6},
        {"cyclic_3", {"x + y + z", "x*y + y*z + z*x", "x*y*z - 1"}, {"x", "y", "z"}, 6}
    };
    
    std::cout << "Testing polynomial systems (suppressing debug output)...\n";
    std::cout << "=" << std::string(60, '=') << "\n";
    
    const TestCase* first_failure = nullptr;
    
    for (const auto& test : tests) {
        // Redirect all output to /dev/null
        int devnull = open("/dev/null", O_WRONLY);
        dup2(devnull, 1);  // stdout
        dup2(devnull, 2);  // stderr
        close(devnull);
        
        auto result = solve_polynomial_system_enhanced(test.polys, test.vars, config);
        
        // Restore output
        dup2(saved_stdout, 1);
        dup2(saved_stderr, 2);
        
        std::cout << std::left << std::setw(15) << test.name << ": ";
        
        if (result.success) {
            std::cout << "✓ SUCCESS (" << result.solutions.size() << " solutions";
            if (test.expected_solutions > 0) {
                if (result.solutions.size() == static_cast<size_t>(test.expected_solutions)) {
                    std::cout << ", as expected";
                } else {
                    std::cout << ", expected " << test.expected_solutions;
                }
            }
            std::cout << ")\n";
        } else {
            std::cout << "✗ FAILED: " << result.error_message << "\n";
            if (!first_failure) {
                first_failure = &test;
            }
        }
    }
    
    std::cout << "=" << std::string(60, '=') << "\n";
    
    if (first_failure) {
        std::cout << "\nFirst failing test: " << first_failure->name << "\n";
        std::cout << "Polynomials: ";
        for (const auto& p : first_failure->polys) {
            std::cout << p << "; ";
        }
        std::cout << "\nVariables: ";
        for (const auto& v : first_failure->vars) {
            std::cout << v << " ";
        }
        std::cout << "\n\nRerunning with verbose output for debugging...\n";
        std::cout << "-" << std::string(60, '-') << "\n";
        
        // Rerun with verbose output
        config.verbose = true;
        auto debug_result = solve_polynomial_system_enhanced(
            first_failure->polys, 
            first_failure->vars, 
            config
        );
        
        if (!debug_result.success) {
            std::cout << "\nFinal error: " << debug_result.error_message << "\n";
        }
    } else {
        std::cout << "\nAll tests passed!\n";
    }
    
    close(saved_stdout);
    close(saved_stderr);
    
    return first_failure ? 1 : 0;
}