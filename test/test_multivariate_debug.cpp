#include "../src/julia_rur/polynomial_solver.hpp"
#include <iostream>
#include <cmath>

using namespace julia_rur;

void test_circle_parabola() {
    std::cout << "=== Debug: Circle and parabola intersection ===\n";
    std::cout << "System: x^2 + y^2 - 4 = 0, x^2 - y - 1 = 0\n\n";
    
    // From first equation: x^2 + y^2 = 4
    // From second equation: x^2 = y + 1
    // Substituting: (y + 1) + y^2 = 4
    // y^2 + y - 3 = 0
    // y = (-1 ± √13) / 2 ≈ 1.303, -2.303
    // x^2 = y + 1, so x = ±√(y + 1)
    
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2-4",    // Circle: x^2 + y^2 = 4
        "1*x^2-1*y-1"       // Parabola: x^2 = y + 1
    };
    std::vector<std::string> variables = {"x", "y"};
    
    RURConfig config;
    config.verbose = true;  // Enable verbose output
    
    auto solution = solve_polynomial_system_complete(polynomials, variables, config);
    print_solution(solution);
    
    std::cout << "\nExpected solutions:\n";
    std::cout << "y = (-1 + √13) / 2 ≈ 1.303, x = ±√(y + 1) ≈ ±1.549\n";
    std::cout << "y = (-1 - √13) / 2 ≈ -2.303, x = ±√(y + 1) (complex)\n";
    
    if (solution.success) {
        std::cout << "\nVerifying solutions:\n";
        for (size_t i = 0; i < solution.solutions.size(); ++i) {
            double x = solution.solutions[i][0].real();
            double y = solution.solutions[i][1].real();
            double x_imag = solution.solutions[i][0].imag();
            double y_imag = solution.solutions[i][1].imag();
            
            std::cout << "Solution " << (i+1) << ": ";
            if (solution.is_real_solution[i]) {
                std::cout << "x = " << x << ", y = " << y << " (real)\n";
            } else {
                std::cout << "x = " << x << " + " << x_imag << "i, ";
                std::cout << "y = " << y << " + " << y_imag << "i (complex)\n";
            }
            
            // Check if solution satisfies equations
            std::complex<double> xc = solution.solutions[i][0];
            std::complex<double> yc = solution.solutions[i][1];
            
            std::complex<double> circle_val = xc*xc + yc*yc - 4.0;
            std::complex<double> parabola_val = xc*xc - yc - 1.0;
            
            std::cout << "  Circle residual: |" << std::abs(circle_val) << "|\n";
            std::cout << "  Parabola residual: |" << std::abs(parabola_val) << "|\n";
        }
    }
}

int main() {
    test_circle_parabola();
    return 0;
}