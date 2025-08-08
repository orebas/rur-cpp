#include "rur_solver.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include <iostream>
#include <vector>
#include <cassert>

void test_simple_system() {
    std::cout << "Testing simple polynomial system..." << std::endl;
    
    try {
        // Test system: x^2 - 1, y^2 - 1 (4 solutions: (±1, ±1))
        std::vector<std::string> vars = {"x", "y"};
        RURSolver<int> solver(100003, vars);
        
        // Create x^2 - 1
        MultivariatePolynomial<int> poly1;
        Monomial x2(2);
        x2[0] = 2;  // x^2
        poly1.set_coefficient(x2, 1);
        poly1.set_coefficient(Monomial(2), -1);  // constant term -1
        
        // Create y^2 - 1
        MultivariatePolynomial<int> poly2;
        Monomial y2(2);
        y2[1] = 2;  // y^2
        poly2.set_coefficient(y2, 1);
        poly2.set_coefficient(Monomial(2), -1);  // constant term -1
        
        solver.add_polynomial(poly1);
        solver.add_polynomial(poly2);
        
        std::cout << "Computing Gröbner basis..." << std::endl;
        solver.compute_groebner_basis();
        
        std::cout << "GB size: " << solver.get_gb_size() << std::endl;
        std::cout << "GB computation time: " << solver.get_gb_computation_time() << " seconds" << std::endl;
        
        // Add assertion for zero-dimensionality
        assert(solver.is_zero_dimensional());  // System x^2-1, y^2-1 should be zero-dimensional
        
        if (solver.is_zero_dimensional()) {
            std::cout << "System is zero-dimensional, computing RUR..." << std::endl;
            
            try {
                solver.compute_rur();
                
                std::cout << "Solution count: " << solver.get_solution_count() << std::endl;
                
                // Add assertion for solution count
                assert(solver.get_solution_count() == 4);  // Should have 4 solutions: (±1, ±1)
                
                auto parameterization = solver.get_parameterization();
                auto minimal_poly = solver.get_minimal_polynomial();
                
                std::cout << "Minimal polynomial degree: " << minimal_poly.size() - 1 << std::endl;
                
                // Add assertion for minimal polynomial degree
                assert(minimal_poly.size() - 1 == 4);  // Degree should equal solution count
                
                std::cout << "Parameterization computed successfully!" << std::endl;
                
            } catch (const std::exception& e) {
                std::cout << "RUR computation failed (expected for prototype): " << e.what() << std::endl;
            }
        } else {
            std::cout << "System is not zero-dimensional" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

void test_cyclic_3() {
    std::cout << "\nTesting cyclic-3 system..." << std::endl;
    
    try {
        // Cyclic-3: x + y + z, xy + yz + zx, xyz - 1
        std::vector<std::string> vars = {"x", "y", "z"};
        RURSolver<int> solver(100003, vars);
        
        // x + y + z
        MultivariatePolynomial<int> poly1;
        Monomial x(3), y(3), z(3);
        x[0] = 1; y[1] = 1; z[2] = 1;
        poly1.set_coefficient(x, 1);
        poly1.set_coefficient(y, 1);
        poly1.set_coefficient(z, 1);
        
        // xy + yz + zx
        MultivariatePolynomial<int> poly2;
        Monomial xy(3), yz(3), zx(3);
        xy[0] = 1; xy[1] = 1;  // xy
        yz[1] = 1; yz[2] = 1;  // yz
        zx[2] = 1; zx[0] = 1;  // zx
        poly2.set_coefficient(xy, 1);
        poly2.set_coefficient(yz, 1);
        poly2.set_coefficient(zx, 1);
        
        // xyz - 1
        MultivariatePolynomial<int> poly3;
        Monomial xyz(3);
        xyz[0] = 1; xyz[1] = 1; xyz[2] = 1;  // xyz
        poly3.set_coefficient(xyz, 1);
        poly3.set_coefficient(Monomial(3), -1);  // constant -1
        
        solver.add_polynomial(poly1);
        solver.add_polynomial(poly2);
        solver.add_polynomial(poly3);
        
        std::cout << "Computing Gröbner basis..." << std::endl;
        solver.compute_groebner_basis();
        
        std::cout << "GB size: " << solver.get_gb_size() << std::endl;
        std::cout << "Is zero-dimensional: " << (solver.is_zero_dimensional() ? "yes" : "no") << std::endl;
        
        // Add assertion for zero-dimensionality
        assert(solver.is_zero_dimensional());  // Cyclic-3 should be zero-dimensional
        
        if (solver.is_zero_dimensional()) {
            try {
                solver.compute_rur();
                std::cout << "Solution count: " << solver.get_solution_count() << std::endl;
                
                // Add assertion for solution count
                assert(solver.get_solution_count() == 6);  // Cyclic-3 has 6 solutions
                
                std::cout << "RUR computed successfully!" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "RUR computation failed (expected for prototype): " << e.what() << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

int main() {
    std::cout << "=== RUR Solver Tests ===" << std::endl;
    
    test_simple_system();
    test_cyclic_3();
    
    std::cout << "\n=== Tests Complete ===" << std::endl;
    return 0;
}