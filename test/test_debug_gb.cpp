#include "f4_solver.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include <iostream>
#include <vector>

void debug_simple_system() {
    std::cout << "Debugging simple system: x^2 - 1, y^2 - 1" << std::endl;
    
    try {
        std::vector<std::string> vars = {"x", "y"};
        F4Solver<int> solver(100003, vars);
        
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
        
        std::cout << "Input polynomial 1: " << poly1.to_string(vars) << std::endl;
        std::cout << "Input polynomial 2: " << poly2.to_string(vars) << std::endl;
        
        solver.add_polynomial(poly1);
        solver.add_polynomial(poly2);
        
        solver.compute_groebner_basis();
        
        auto gb = solver.get_groebner_basis();
        std::cout << "\nGröbner basis (" << gb.size() << " polynomials):" << std::endl;
        
        for (size_t i = 0; i < gb.size(); ++i) {
            std::cout << "GB[" << i << "]: " << gb[i].to_string(vars) << std::endl;
            std::cout << "  Leading monomial: " << gb[i].leading_monomial().to_string(vars) << std::endl;
            std::cout << "  Degree: " << gb[i].degree() << std::endl;
            
            // Count non-zero variables in leading monomial
            auto leading_mon = gb[i].leading_monomial();
            size_t non_zero_vars = 0;
            for (size_t j = 0; j < vars.size(); ++j) {
                if (leading_mon[j] > 0) {
                    non_zero_vars++;
                    std::cout << "  Variable " << vars[j] << " has exponent " << leading_mon[j] << std::endl;
                }
            }
            std::cout << "  Non-zero variables: " << non_zero_vars << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

void debug_cyclic_3() {
    std::cout << "\nDebugging cyclic-3 system" << std::endl;
    
    try {
        std::vector<std::string> vars = {"x", "y", "z"};
        F4Solver<int> solver(100003, vars);
        
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
        
        std::cout << "Input polynomial 1: " << poly1.to_string(vars) << std::endl;
        std::cout << "Input polynomial 2: " << poly2.to_string(vars) << std::endl;
        std::cout << "Input polynomial 3: " << poly3.to_string(vars) << std::endl;
        
        solver.add_polynomial(poly1);
        solver.add_polynomial(poly2);
        solver.add_polynomial(poly3);
        
        solver.compute_groebner_basis();
        
        auto gb = solver.get_groebner_basis();
        std::cout << "\nGröbner basis (" << gb.size() << " polynomials):" << std::endl;
        
        for (size_t i = 0; i < gb.size(); ++i) {
            std::cout << "GB[" << i << "]: " << gb[i].to_string(vars) << std::endl;
            std::cout << "  Leading monomial: " << gb[i].leading_monomial().to_string(vars) << std::endl;
            std::cout << "  Degree: " << gb[i].degree() << std::endl;
            
            // Count non-zero variables in leading monomial
            auto leading_mon = gb[i].leading_monomial();
            size_t non_zero_vars = 0;
            for (size_t j = 0; j < vars.size(); ++j) {
                if (leading_mon[j] > 0) {
                    non_zero_vars++;
                    std::cout << "  Variable " << vars[j] << " has exponent " << leading_mon[j] << std::endl;
                }
            }
            std::cout << "  Non-zero variables: " << non_zero_vars << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}

int main() {
    debug_simple_system();
    debug_cyclic_3();
    return 0;
}