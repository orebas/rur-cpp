#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include "../src/julia_rur/polynomial_evaluator.hpp"
#include <iostream>
#include <iomanip>
#include <unistd.h>  // For dup, dup2, close

using namespace julia_rur;

void test_simple_case(const std::string& name, 
                     const std::vector<std::string>& polys,
                     const std::vector<std::string>& vars) {
    std::cout << "\n" << name << ":" << std::endl;
    std::cout << "  System: ";
    for (size_t i = 0; i < polys.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << polys[i];
    }
    std::cout << std::endl;
    
    EnhancedSolverConfig config;
    config.verbose = false;
    
    // Redirect stderr to suppress debug output
    int old_stderr = dup(2);
    freopen("/dev/null", "w", stderr);
    
    auto result = solve_polynomial_system_enhanced(polys, vars, config);
    
    // Restore stderr
    dup2(old_stderr, 2);
    close(old_stderr);
    
    if (!result.success) {
        std::cout << "  FAILED: " << result.error_message << std::endl;
        return;
    }
    
    std::cout << "  Found " << result.solutions.size() << " solutions:" << std::endl;
    
    for (size_t i = 0; i < result.solutions.size(); ++i) {
        std::cout << "    Solution " << (i+1) << ": ";
        for (size_t j = 0; j < vars.size(); ++j) {
            if (j > 0) std::cout << ", ";
            auto val = result.solutions[i][j];
            
            // Pretty print
            if (std::abs(val.imag()) < 1e-10) {
                std::cout << vars[j] << " = " << val.real();
            } else {
                std::cout << vars[j] << " = " << val;
            }
        }
        
        // Compute residual
        double residual = PolynomialEvaluator::compute_residual(polys, vars, result.solutions[i]);
        std::cout << " (residual: " << std::scientific << std::setprecision(2) << residual << ")" << std::endl;
    }
}

int main() {
    std::cout << "\n==================================\n";
    std::cout << "Simple Multiplicity Tests\n";
    std::cout << "==================================\n";
    
    // Test 1: Simple double root
    test_simple_case("Double root (x-1)^2", 
                    {"(x-1)^2"}, 
                    {"x"});
    
    // Test 2: Triple root
    test_simple_case("Triple root (x-2)^3", 
                    {"(x-2)^3"}, 
                    {"x"});
    
    // Test 3: Product of powers
    test_simple_case("Product (x-1)^2 * (x+1)^2",
                    {"(x-1)^2 * (x+1)^2"},
                    {"x"});
    
    // Test 4: Tangent line to circle
    test_simple_case("Tangent: x^2+y^2-1, y-1",
                    {"x^2 + y^2 - 1", "y - 1"},
                    {"x", "y"});
    
    // Test 5: Standard polynomial
    test_simple_case("Standard: x^2 - 4",
                    {"x^2 - 4"},
                    {"x"});
    
    std::cout << "\n==================================\n";
    return 0;
}