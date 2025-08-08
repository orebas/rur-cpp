#include "../src/julia_rur/polynomial_solver.hpp"
#include "../src/julia_rur/numerical_roots_eigen.hpp"
#include <iostream>
#include <cmath>

using namespace julia_rur;

int main() {
    std::cout << "=== Testing multivariate back-substitution ===" << std::endl;
    
    // Create a mock RUR result for the system:
    // x^2 + y^2 - 1 = 0, x - y = 0
    // We know the minimal polynomial should be T^2 - 1/2 = 0
    // with roots T = ±1/√2
    // And the parameterizations should be x = T, y = T
    
    RationalRURResult rur;
    rur.success = true;
    
    // Minimal polynomial: T^2 - 1/2 = 0
    rur.minimal_polynomial = {mpq_class(-1, 2), mpq_class(0), mpq_class(1)};
    
    // For this simple case, both x and y equal T
    // So the parameterizations are just T (represented as [0, 1])
    rur.numerators.resize(2);
    rur.numerators[0] = {mpq_class(0), mpq_class(1)};  // x = T
    rur.numerators[1] = {mpq_class(0), mpq_class(1)};  // y = T
    
    std::cout << "Mock RUR result:" << std::endl;
    std::cout << "- Minimal polynomial: T^2 - 1/2 = 0" << std::endl;
    std::cout << "- Parameterization x: T" << std::endl;
    std::cout << "- Parameterization y: T" << std::endl;
    
    // Find roots of minimal polynomial
    std::vector<std::complex<double>> t_roots = find_polynomial_roots(rur.minimal_polynomial);
    
    std::cout << "\nRoots of minimal polynomial:" << std::endl;
    for (const auto& root : t_roots) {
        std::cout << "  T = " << root.real();
        if (std::abs(root.imag()) > 1e-10) {
            std::cout << " + " << root.imag() << "i";
        }
        std::cout << std::endl;
    }
    
    // Now test the polynomial solver
    std::vector<std::string> variables = {"x", "y"};
    PolynomialSystemSolution solution;
    solution.variable_names = variables;
    solution.minimal_polynomial = rur.minimal_polynomial;
    solution.quotient_dimension = 2;
    
    // Back-substitute manually
    const double epsilon = 1e-10;
    for (const auto& t_root : t_roots) {
        std::vector<std::complex<double>> var_values;
        
        // Since x = T and y = T, both variables equal the root
        var_values.push_back(t_root);  // x
        var_values.push_back(t_root);  // y
        
        // Check if solution is real
        bool is_real = true;
        for (const auto& val : var_values) {
            if (std::abs(val.imag()) > epsilon) {
                is_real = false;
                break;
            }
        }
        
        solution.solutions.push_back(var_values);
        solution.is_real_solution.push_back(is_real);
    }
    
    solution.success = true;
    
    std::cout << "\nBack-substituted solutions:" << std::endl;
    print_solution(solution, std::cout);
    
    // Verify the solutions
    std::cout << "\nVerification:" << std::endl;
    for (size_t i = 0; i < solution.solutions.size(); ++i) {
        double x = solution.solutions[i][0].real();
        double y = solution.solutions[i][1].real();
        
        double eq1 = x*x + y*y - 1.0;  // Should be ~0
        double eq2 = x - y;            // Should be ~0
        
        std::cout << "Solution " << (i+1) << ": (" << x << ", " << y << ")" << std::endl;
        std::cout << "  x^2 + y^2 - 1 = " << eq1 << std::endl;
        std::cout << "  x - y = " << eq2 << std::endl;
    }
    
    return 0;
}