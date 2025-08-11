#include "src/flint_rur_solver.hpp"
#include <iostream>

int main() {
    try {
        mp_limb_t prime = 100003;
        std::vector<std::string> vars = {"x", "y"};
        
        // Create simple system: x^2 - 1, y^2 - 1
        MultivariatePolynomial<int> poly1;
        poly1.set_coefficient(Monomial({2, 0}), 1);   // x^2
        poly1.set_coefficient(Monomial({0, 0}), -1);  // -1
        
        MultivariatePolynomial<int> poly2;
        poly2.set_coefficient(Monomial({0, 2}), 1);   // y^2  
        poly2.set_coefficient(Monomial({0, 0}), -1);  // -1
        
        std::vector<MultivariatePolynomial<int>> gb = {poly1, poly2};
        
        std::cout << "Polynomials:" << std::endl;
        std::cout << "p1: " << poly1.to_string(vars) << std::endl;
        std::cout << "p2: " << poly2.to_string(vars) << std::endl;
        
        // Create polynomial basis
        auto basis = PolynomialBasis<int>::create(gb, 2, prime, ORD_DEGREVLEX);
        
        std::cout << "\nQuotient basis monomials:" << std::endl;
        for (size_t i = 0; i < basis->quotient_basis_monomials.size(); i++) {
            std::cout << i << ": " << basis->quotient_basis_monomials[i].to_string(vars) << std::endl;
        }
        
        // Create reducer and test some reductions
        FLINTReducer<int> reducer(basis);
        
        // Test reducing x*1, y*1, x*y*1, etc.
        std::vector<Monomial> test_monomials = {
            Monomial({1, 0}),  // x
            Monomial({0, 1}),  // y  
            Monomial({1, 1}),  // xy
            Monomial({2, 0}),  // x^2 (should reduce)
            Monomial({0, 2}),  // y^2 (should reduce)
        };
        
        for (const auto& mon : test_monomials) {
            MultivariatePolynomial<int> test_poly;
            test_poly.set_coefficient(mon, 1);
            auto reduced = reducer.reduce(test_poly);
            
            std::cout << "\nReducing " << mon.to_string(vars) << ":" << std::endl;
            std::cout << "  Coefficients: ";
            for (size_t i = 0; i < reduced.size(); i++) {
                std::cout << reduced[i] << " ";
            }
            std::cout << std::endl;
        }
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}