#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>
#include <cmath>

using namespace julia_rur;

void test_univariate() {
    std::cout << "=== Test 1: Univariate polynomial x^2 - 2 ===" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x^2-2"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    // Verify solution
    if (solution.success) {
        double expected = std::sqrt(2.0);
        bool found_positive = false, found_negative = false;
        
        for (const auto& sol : solution.solutions) {
            double x = sol[0].real();
            if (std::abs(x - expected) < 1e-9) found_positive = true;
            if (std::abs(x + expected) < 1e-9) found_negative = true;
        }
        
        if (found_positive && found_negative) {
            std::cout << "\n✓ Correctly found ±√2" << std::endl;
        } else {
            std::cout << "\n✗ Failed to find correct roots" << std::endl;
        }
    }
    std::cout << std::endl;
}

void test_linear() {
    std::cout << "=== Test 2: Linear system x - 3 ===" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x-3"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    if (solution.success && solution.solutions.size() == 1) {
        double x = solution.solutions[0][0].real();
        if (std::abs(x - 3.0) < 1e-9) {
            std::cout << "\n✓ Correctly found x = 3" << std::endl;
        } else {
            std::cout << "\n✗ Wrong solution: expected 3, got " << x << std::endl;
        }
    }
    std::cout << std::endl;
}

void test_cubic() {
    std::cout << "=== Test 3: Cubic polynomial x^3 - 2x - 5 ===" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x^3-2*x-5"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    // Known real root: approximately 2.094551482
    if (solution.success) {
        bool found_real_root = false;
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            if (solution.is_real_solution[i]) {
                double x = solution.solutions[i][0].real();
                if (std::abs(x - 2.094551482) < 1e-8) {
                    found_real_root = true;
                    std::cout << "\n✓ Found expected real root" << std::endl;
                }
            }
        }
        if (!found_real_root) {
            std::cout << "\n✗ Failed to find expected real root" << std::endl;
        }
    }
    std::cout << std::endl;
}

void test_multivariate() {
    std::cout << "=== Test 4: Multivariate system (circle and line) ===" << std::endl;
    std::cout << "System: x^2 + y^2 - 1 = 0, x - y = 0" << std::endl;
    
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2-1",  // Circle
        "1*x-1*y"         // Line x = y
    };
    std::vector<std::string> variables = {"x", "y"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    // Expected solutions: (1/√2, 1/√2) and (-1/√2, -1/√2)
    if (solution.success) {
        double expected = 1.0 / std::sqrt(2.0);
        bool found_positive = false, found_negative = false;
        
        for (const auto& sol : solution.solutions) {
            if (sol.size() >= 2) {
                double x = sol[0].real();
                double y = sol[1].real();
                
                if (std::abs(x - expected) < 1e-9 && std::abs(y - expected) < 1e-9) {
                    found_positive = true;
                }
                if (std::abs(x + expected) < 1e-9 && std::abs(y + expected) < 1e-9) {
                    found_negative = true;
                }
            }
        }
        
        if (found_positive && found_negative) {
            std::cout << "\n✓ Correctly found intersection points" << std::endl;
        } else {
            std::cout << "\n✗ Failed to find correct intersection points" << std::endl;
        }
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "Testing Complete Polynomial System Solver" << std::endl;
    std::cout << "=========================================" << std::endl << std::endl;
    
    test_univariate();
    test_linear();
    test_cubic();
    test_multivariate(); // F4 issue has been fixed
    
    return 0;
}