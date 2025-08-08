#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>
#include <cmath>

using namespace julia_rur;

void test_simple_quadratic() {
    std::cout << "=== Test: Simple quadratic system ===\n";
    std::cout << "System: y - 1 = 0, x^2 - 2 = 0\n\n";
    std::cout << "Expected: y = 1, x = ±√2\n\n";
    
    std::vector<std::string> polynomials = {
        "1*y-1",      // y = 1
        "1*x^2-2"     // x^2 = 2
    };
    std::vector<std::string> variables = {"x", "y"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    if (solution.success) {
        std::cout << "\nVerifying solutions:\n";
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            if (solution.is_real_solution[i]) {
                double x = solution.solutions[i][0].real();
                double y = solution.solutions[i][1].real();
                
                std::cout << "Solution " << (i+1) << ": x = " << x << ", y = " << y << "\n";
                std::cout << "  Error in y - 1 = " << std::abs(y - 1) << "\n";
                std::cout << "  Error in x^2 - 2 = " << std::abs(x*x - 2) << "\n";
            }
        }
    }
}

void test_coupled_quadratic() {
    std::cout << "\n=== Test: Coupled quadratic system ===\n";
    std::cout << "System: y^2 - 1 = 0, x^2 - y - 1 = 0\n\n";
    std::cout << "Expected: y = ±1, x^2 = y + 1\n";
    std::cout << "For y = 1: x = ±√2\n";
    std::cout << "For y = -1: x = 0\n\n";
    
    std::vector<std::string> polynomials = {
        "1*y^2-1",        // y^2 = 1
        "1*x^2-1*y-1"     // x^2 = y + 1
    };
    std::vector<std::string> variables = {"x", "y"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    if (solution.success) {
        std::cout << "\nVerifying solutions:\n";
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            if (solution.is_real_solution[i]) {
                double x = solution.solutions[i][0].real();
                double y = solution.solutions[i][1].real();
                
                std::cout << "Solution " << (i+1) << ": x = " << x << ", y = " << y << "\n";
                std::cout << "  Error in y^2 - 1 = " << std::abs(y*y - 1) << "\n";
                std::cout << "  Error in x^2 - y - 1 = " << std::abs(x*x - y - 1) << "\n";
            }
        }
    }
}

int main() {
    test_simple_quadratic();
    test_coupled_quadratic();
    return 0;
}