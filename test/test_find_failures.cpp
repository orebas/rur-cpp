/**
 * @file test_find_failures.cpp
 * @brief Find the simplest failing test case
 */

#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include "../src/julia_rur/polynomial_evaluator.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <unistd.h>  // For dup, dup2, close

using namespace julia_rur;

struct TestCase {
    std::string name;
    std::vector<std::string> polynomials;
    std::vector<std::string> variables;
    int expected_solutions;
};

void run_test(const TestCase& test) {
    std::cout << "\n" << test.name << ": ";
    std::cout.flush();
    
    EnhancedSolverConfig config;
    config.verbose = false;
    
    // Redirect stderr to suppress debug output
    int saved_stderr = dup(2);
    freopen("/dev/null", "w", stderr);
    
    auto result = solve_polynomial_system_enhanced(test.polynomials, test.variables, config);
    
    // Restore stderr
    dup2(saved_stderr, 2);
    close(saved_stderr);
    
    if (result.success) {
        std::cout << "SUCCESS (" << result.solutions.size() << " solutions";
        if (test.expected_solutions >= 0) {
            if (result.solutions.size() == static_cast<size_t>(test.expected_solutions)) {
                std::cout << " ✓";
            } else {
                std::cout << ", expected " << test.expected_solutions << " ✗";
            }
        }
        std::cout << ")";
        
        // Check residuals
        double max_residual = 0.0;
        for (const auto& sol : result.solutions) {
            double res = PolynomialEvaluator::compute_residual(test.polynomials, test.variables, sol);
            max_residual = std::max(max_residual, res);
        }
        if (max_residual > 1e-6) {
            std::cout << " [residual: " << std::scientific << max_residual << "]";
        }
    } else {
        std::cout << "FAILED: " << result.error_message;
    }
}

int main() {
    std::cout << "===========================================\n";
    std::cout << "Finding Simplest Failing Test Cases\n";
    std::cout << "===========================================\n";
    
    std::vector<TestCase> tests = {
        // Start with the absolute simplest
        {"Linear x = 1", {"x - 1"}, {"x"}, 1},
        {"Quadratic x^2 = 4", {"x^2 - 4"}, {"x"}, 2},
        {"Cubic x^3 = 8", {"x^3 - 8"}, {"x"}, 3},
        
        // Two variable linear
        {"Linear 2x2 trivial", {"x", "y"}, {"x", "y"}, 1},
        {"Linear 2x2 simple", {"x + y - 2", "x - y"}, {"x", "y"}, 1},
        {"Linear 2x2 standard", {"2*x + y - 3", "x - y - 1"}, {"x", "y"}, 1},
        
        // Simple nonlinear 2x2
        {"Circle-Line", {"x^2 + y^2 - 1", "y"}, {"x", "y"}, 2},
        {"Parabola-Line", {"y - x^2", "y - 1"}, {"x", "y"}, 2},
        {"Two circles", {"x^2 + y^2 - 1", "(x-1)^2 + y^2 - 1"}, {"x", "y"}, 2},
        
        // Three variables - simplest
        {"Linear 3x3 trivial", {"x", "y", "z"}, {"x", "y", "z"}, 1},
        {"Linear 3x3 simple", {"x + y + z - 3", "x - y", "y - z"}, {"x", "y", "z"}, 1},
        
        // Classic small systems that should work
        {"Elementary symmetric 2", {"x + y - 2", "x*y - 1"}, {"x", "y"}, 2},
        {"Elementary symmetric 3", {"x + y + z - 3", "x*y + y*z + z*x - 3", "x*y*z - 1"}, {"x", "y", "z"}, 6},
        
        // The problematic ones
        {"Cyclic-3", {"x + y + z", "x*y + y*z + z*x", "x*y*z - 1"}, {"x", "y", "z"}, 6},
        {"Katsura-3", {"x + 2*y + 2*z - 1", "x^2 + 2*y^2 + 2*z^2 - x", "2*x*y + 2*y*z - y"}, {"x", "y", "z"}, 4},
        
        // Edge cases
        {"Empty system", {}, {}, -1},
        {"Overdetermined", {"x - 1", "x - 2"}, {"x"}, 0},
        {"Underdetermined", {"x + y - 1"}, {"x", "y"}, -1},
        {"Constant true", {"0"}, {"x"}, -1},
        {"Constant false", {"1"}, {"x"}, 0},
    };
    
    int passed = 0;
    int failed = 0;
    std::vector<std::string> failures;
    
    for (const auto& test : tests) {
        run_test(test);
        
        // Check if last test failed
        EnhancedSolverConfig config;
        config.verbose = false;
        freopen("/dev/null", "w", stderr);
        auto result = solve_polynomial_system_enhanced(test.polynomials, test.variables, config);
        if (!result.success) {
            failures.push_back(test.name);
            failed++;
        } else {
            passed++;
        }
    }
    
    std::cout << "\n\n===========================================\n";
    std::cout << "Summary: " << passed << " passed, " << failed << " failed\n";
    
    if (!failures.empty()) {
        std::cout << "\nFailed tests (simplest first):\n";
        for (const auto& name : failures) {
            std::cout << "  - " << name << "\n";
        }
        
        std::cout << "\n🎯 Simplest failure: " << failures[0] << "\n";
    }
    
    std::cout << "===========================================\n";
    
    return failed > 0 ? 1 : 0;
}