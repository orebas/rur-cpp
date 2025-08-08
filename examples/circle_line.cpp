/**
 * @file circle_line.cpp
 * @brief Example: Finding intersection of a circle and line using RUR
 * 
 * This example demonstrates solving the system:
 *   x² + y² = 1  (unit circle)
 *   x = y        (diagonal line)
 * 
 * The solutions are the two points where the line intersects the circle:
 *   (1/√2, 1/√2) and (-1/√2, -1/√2)
 */

#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace julia_rur;

int main() {
    std::cout << "Circle-Line Intersection Example\n";
    std::cout << "================================\n\n";
    
    // Define the polynomial system
    // Note: We use modular arithmetic, so -1 is represented as 100002 (mod 100003)
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2+100002",  // x² + y² - 1 = 0
        "1*x+100002*y"         // x - y = 0
    };
    
    std::vector<std::string> variables = {"x", "y"};
    
    std::cout << "System of equations:\n";
    std::cout << "  x² + y² = 1\n";
    std::cout << "  x = y\n\n";
    
    // Check if the system is zero-dimensional (finite number of solutions)
    if (!is_zero_dimensional_system(polynomials, variables)) {
        std::cerr << "Error: System is not zero-dimensional!\n";
        return 1;
    }
    
    // Get the number of solutions
    int dim = compute_quotient_dimension(polynomials, variables);
    std::cout << "Number of solutions (with multiplicity): " << dim << "\n\n";
    
    // Compute the RUR with default configuration
    RURConfig config;
    config.verbose = false;  // Set to true for detailed output
    
    std::cout << "Computing Rational Univariate Representation...\n";
    RationalRURResult result = compute_rational_rur(polynomials, variables, config);
    
    if (!result.success) {
        std::cerr << "Error: " << result.error_message << "\n";
        return 1;
    }
    
    // Display the results
    std::cout << "\n" << format_rur_result(result, variables) << "\n";
    
    // Show the minimal polynomial explicitly
    std::cout << "Minimal polynomial coefficients (low to high degree):\n";
    for (size_t i = 0; i < result.minimal_polynomial.size(); ++i) {
        std::cout << "  c[" << i << "] = " << result.minimal_polynomial[i] << "\n";
    }
    
    // Interpret the results
    std::cout << "\nInterpretation:\n";
    std::cout << "The minimal polynomial f(T) = ";
    
    // Print polynomial in readable form
    bool first = true;
    for (int i = result.minimal_polynomial.size() - 1; i >= 0; --i) {
        auto coeff = result.minimal_polynomial[i];
        if (coeff != 0) {
            if (!first && coeff > 0) std::cout << " + ";
            std::cout << coeff;
            if (i > 0) {
                std::cout << "*T";
                if (i > 1) std::cout << "^" << i;
            }
            first = false;
        }
    }
    std::cout << "\n\n";
    
    // For this example, we expect f(T) = T² - 1/2
    // The roots are T = ±1/√2
    std::cout << "The roots of f(T) give the x-coordinates (and y-coordinates) of the intersection points.\n";
    std::cout << "Since x = y on the line, and we have a separating element T,\n";
    std::cout << "the solutions are approximately:\n";
    std::cout << "  T = ±" << std::setprecision(6) << 1.0/std::sqrt(2) << "\n";
    std::cout << "  giving (x,y) = (±0.707107, ±0.707107)\n";
    
    return 0;
}