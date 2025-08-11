#include "src/flint_reducer_v2.hpp"
#include "src/polynomial_basis.hpp"
#include "src/polynomial.hpp"
#include <iostream>

int main() {
    try {
        // Create a simple polynomial system x^2 - 1, y^2 - 1
        std::vector<MultivariatePolynomial<int>> gb;
        
        MultivariatePolynomial<int> poly1;
        poly1.set_coefficient(Monomial({2, 0}), 1);   // x^2
        poly1.set_coefficient(Monomial({0, 0}), -1);  // -1
        gb.push_back(poly1);
        
        MultivariatePolynomial<int> poly2;
        poly2.set_coefficient(Monomial({0, 2}), 1);   // y^2
        poly2.set_coefficient(Monomial({0, 0}), -1);  // -1
        gb.push_back(poly2);
        
        auto basis = PolynomialBasis<int>::create(gb, 2, 100003, ORD_DEGREVLEX);
        
        // Create reducer - this should trigger the linking error
        FLINTReducer<int> reducer(basis);
        
        std::cout << "FLINTReducer created successfully!" << std::endl;
        
        // Now try to actually use it to trigger the divrem_ideal call
        MultivariatePolynomial<int> test_poly;
        test_poly.set_coefficient(Monomial({3, 0}), 1);  // x^3
        
        auto coeffs = reducer.reduce(test_poly);
        
        std::cout << "Reduction completed successfully! Coeffs size: " << coeffs.size() << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}