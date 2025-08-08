#include <iostream>
#include <gmpxx.h>
#include "../src/julia_rur/polynomial_solver.hpp"

int main() {
    std::cout << "Testing polynomial printing..." << std::endl;
    
    // Create a minimal polynomial T^2 - 2 = 0
    // Coefficients are [-2, 0, 1] in low-to-high order
    std::vector<mpq_class> min_poly;
    min_poly.push_back(mpq_class(-2));
    min_poly.push_back(mpq_class(0));
    min_poly.push_back(mpq_class(1));
    
    std::cout << "\nCoefficients in order: ";
    for (size_t i = 0; i < min_poly.size(); ++i) {
        std::cout << "[" << i << "]=" << min_poly[i] << " ";
    }
    std::cout << std::endl;
    
    // Test the printing logic
    std::cout << "\nPrinting polynomial: ";
    bool first = true;
    for (int i = min_poly.size() - 1; i >= 0; --i) {
        const auto& coeff = min_poly[i];
        std::cout << "\n  i=" << i << ", coeff=" << coeff << ", first=" << first;
        
        if (coeff != 0) {
            if (!first && coeff > 0) {
                std::cout << " -> print ' + '";
            }
            if (coeff < 0) {
                std::cout << " -> print ' - '";
            }
            
            mpq_class abs_coeff = abs(coeff);
            std::cout << " -> print '" << abs_coeff << "'";
            
            if (i > 0) {
                std::cout << " -> print 'T'";
                if (i > 1) {
                    std::cout << " -> print '^" << i << "'";
                }
            }
            first = false;
        }
    }
    std::cout << std::endl;
    
    // Now test with actual solution
    julia_rur::PolynomialSystemSolution solution;
    solution.success = true;
    solution.variable_names = {"x"};
    solution.minimal_polynomial = min_poly;
    solution.quotient_dimension = 2;
    solution.solutions.push_back({std::complex<double>(0, sqrt(2))});
    solution.solutions.push_back({std::complex<double>(0, -sqrt(2))});
    solution.is_real_solution = {false, false};
    
    std::cout << "\nActual output from print_solution:" << std::endl;
    julia_rur::print_solution(solution, std::cout);
    
    return 0;
}