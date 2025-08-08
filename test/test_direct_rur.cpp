#include "direct_rur_solver.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include <iostream>
#include <vector>
#include <cassert>

void test_simple_known_gb() {
    std::cout << "Testing with known simple Gröbner basis..." << std::endl;
    
    try {
        std::vector<std::string> vars = {"x", "y"};
        DirectRURSolver<int> solver(100003, vars);
        
        // Create known GB for x^2 - 1, y^2 - 1
        // The GB should be {x^2 - 1, y^2 - 1} (already reduced)
        std::vector<MultivariatePolynomial<int>> gb;
        
        // x^2 - 1
        MultivariatePolynomial<int> poly1;
        Monomial x2(2);
        x2[0] = 2;  // x^2
        poly1.set_coefficient(x2, 1);
        poly1.set_coefficient(Monomial(2), -1);  // constant term -1
        gb.push_back(poly1);
        
        // y^2 - 1
        MultivariatePolynomial<int> poly2;
        Monomial y2(2);
        y2[1] = 2;  // y^2
        poly2.set_coefficient(y2, 1);
        poly2.set_coefficient(Monomial(2), -1);  // constant term -1
        gb.push_back(poly2);
        
        std::cout << "Setting Gröbner basis:" << std::endl;
        for (size_t i = 0; i < gb.size(); ++i) {
            std::cout << "  GB[" << i << "]: " << gb[i].to_string(vars) << std::endl;
        }
        
        solver.set_groebner_basis(gb);
        
        std::cout << "Is zero-dimensional: " << (solver.is_zero_dimensional() ? "yes" : "no") << std::endl;
        
        // Add assertion for zero-dimensionality
        assert(solver.is_zero_dimensional());  // System x^2-1, y^2-1 should be zero-dimensional
        
        if (solver.is_zero_dimensional()) {
            std::cout << "Computing RUR..." << std::endl;
            solver.compute_rur();
            std::cout << "Solution count: " << solver.get_solution_count() << std::endl;
            
            // Add assertion for solution count
            assert(solver.get_solution_count() == 4);  // Should have 4 solutions: (±1, ±1)
        } else {
            std::cout << "System is not zero-dimensional, cannot compute RUR" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

void test_univariate_system() {
    std::cout << "\nTesting simple univariate system..." << std::endl;
    
    try {
        std::vector<std::string> vars = {"x"};
        DirectRURSolver<int> solver(100003, vars);
        
        // Create GB: {x^3 - 1} (has 3 solutions: cube roots of unity)
        std::vector<MultivariatePolynomial<int>> gb;
        
        MultivariatePolynomial<int> poly;
        Monomial x3(1);
        x3[0] = 3;  // x^3
        poly.set_coefficient(x3, 1);
        poly.set_coefficient(Monomial(1), -1);  // constant term -1
        gb.push_back(poly);
        
        std::cout << "Setting Gröbner basis: " << poly.to_string(vars) << std::endl;
        
        solver.set_groebner_basis(gb);
        
        std::cout << "Is zero-dimensional: " << (solver.is_zero_dimensional() ? "yes" : "no") << std::endl;
        
        // Add assertion for zero-dimensionality
        assert(solver.is_zero_dimensional());  // System x^3-1 should be zero-dimensional
        
        if (solver.is_zero_dimensional()) {
            std::cout << "Computing RUR..." << std::endl;
            solver.compute_rur();
            std::cout << "Solution count: " << solver.get_solution_count() << std::endl;
            
            // Add assertion for solution count
            assert(solver.get_solution_count() == 3);  // Should have 3 solutions: cube roots of unity
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

int main() {
    std::cout << "=== Direct RUR Solver Tests ===" << std::endl;
    
    test_univariate_system();
    test_simple_known_gb();
    
    std::cout << "\n=== Tests Complete ===" << std::endl;
    return 0;
}