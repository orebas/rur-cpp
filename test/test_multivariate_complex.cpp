#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>
#include <cmath>

using namespace julia_rur;

void test_three_variable_system() {
    std::cout << "=== Test: Three variable system ===\n";
    std::cout << "System: x + y + z - 3 = 0, x^2 + y^2 + z^2 - 5 = 0, x*y + y*z + z*x - 4 = 0\n\n";
    
    std::vector<std::string> polynomials = {
        "1*x+1*y+1*z-3",          // x + y + z - 3 = 0  
        "1*x^2+1*y^2+1*z^2-5",    // x^2 + y^2 + z^2 - 5 = 0
        "1*x*y+1*y*z+1*z*x-4"     // xy + yz + zx - 4 = 0
    };
    std::vector<std::string> variables = {"x", "y", "z"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    // Verify symmetry - expect solutions with permutations of coordinates
    if (solution.success) {
        std::cout << "\nChecking solution properties:\n";
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            if (solution.is_real_solution[i]) {
                double x = solution.solutions[i][0].real();
                double y = solution.solutions[i][1].real();
                double z = solution.solutions[i][2].real();
                
                double sum = x + y + z;
                double sum_sq = x*x + y*y + z*z;
                double sum_prod = x*y + y*z + z*x;
                
                std::cout << "Solution " << (i+1) << ": ";
                std::cout << "sum=" << sum << " (expect 3), ";
                std::cout << "sum_sq=" << sum_sq << " (expect 5), ";
                std::cout << "sum_prod=" << sum_prod << " (expect 4)\n";
            }
        }
    }
}

void test_quadratic_curve_intersection() {
    std::cout << "\n=== Test: Two quadratic curves intersection ===\n";
    std::cout << "System: x^2 + y^2 - 4 = 0, x^2 - y - 1 = 0\n\n";
    
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2-4",    // Circle: x^2 + y^2 = 4
        "1*x^2-1*y-1"       // Parabola: x^2 = y + 1
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
                
                double circle_val = x*x + y*y - 4;
                double parabola_val = x*x - y - 1;
                
                std::cout << "Solution " << (i+1) << " at (" << x << ", " << y << "): ";
                std::cout << "circle error=" << std::abs(circle_val) << ", ";
                std::cout << "parabola error=" << std::abs(parabola_val) << "\n";
            }
        }
    }
}

void test_ellipse_hyperbola() {
    std::cout << "\n=== Test: Ellipse and hyperbola intersection ===\n";
    std::cout << "System: x^2/4 + y^2 - 1 = 0, x*y - 1 = 0\n";
    std::cout << "Rewritten as: x^2 + 4*y^2 - 4 = 0, x*y - 1 = 0\n\n";
    
    std::vector<std::string> polynomials = {
        "1*x^2+4*y^2-4",    // Ellipse: x^2/4 + y^2 = 1
        "1*x*y-1"           // Hyperbola: xy = 1
    };
    std::vector<std::string> variables = {"x", "y"};
    
    RURConfig config;
    config.verbose = false;
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    if (solution.success) {
        std::cout << "\nVerifying solutions lie on both curves:\n";
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            if (solution.is_real_solution[i]) {
                double x = solution.solutions[i][0].real();
                double y = solution.solutions[i][1].real();
                
                double ellipse_val = x*x/4.0 + y*y - 1;
                double hyperbola_val = x*y - 1;
                
                std::cout << "Solution " << (i+1) << " at (" << x << ", " << y << "): ";
                std::cout << "ellipse error=" << std::abs(ellipse_val) << ", ";
                std::cout << "hyperbola error=" << std::abs(hyperbola_val) << "\n";
            }
        }
    }
}

int main() {
    std::cout << "Testing Complex Multivariate Polynomial Systems\n";
    std::cout << "==============================================\n\n";
    
    test_quadratic_curve_intersection();
    test_ellipse_hyperbola();
    test_three_variable_system();
    
    return 0;
}