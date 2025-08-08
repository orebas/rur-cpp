#include "monomial.hpp"
#include "polynomial.hpp"
#include <iostream>
#include <vector>
#include <string>

int main() {
    std::cout << "Testing basic polynomial operations...\n";
    
    // Test monomials
    Monomial m1({1, 2, 0});  // x0*x1^2
    Monomial m2({0, 1, 1});  // x1*x2
    
    std::cout << "m1 = " << m1 << std::endl;
    std::cout << "m2 = " << m2 << std::endl;
    std::cout << "m1 * m2 = " << (m1 * m2) << std::endl;
    std::cout << "m1 degree = " << m1.degree() << std::endl;
    
    // Test polynomials
    MultivariatePolynomial<int> p1;
    p1.set_coefficient(Monomial({2, 0, 0}), 3);  // 3*x0^2
    p1.set_coefficient(Monomial({1, 1, 0}), -2); // -2*x0*x1
    p1.set_coefficient(Monomial({0, 0, 0}), 1);  // +1
    
    MultivariatePolynomial<int> p2;
    p2.set_coefficient(Monomial({1, 0, 0}), 1);  // x0
    p2.set_coefficient(Monomial({0, 1, 0}), 1);  // x1
    
    std::vector<std::string> varnames = {"x", "y", "z"};
    
    std::cout << "\nPolynomial tests:\n";
    std::cout << "p1 = " << p1.to_string(varnames) << std::endl;
    std::cout << "p2 = " << p2.to_string(varnames) << std::endl;
    std::cout << "p1 + p2 = " << (p1 + p2).to_string(varnames) << std::endl;
    std::cout << "p1 * p2 = " << (p1 * p2).to_string(varnames) << std::endl;
    
    std::cout << "\nLeading terms:\n";
    std::cout << "p1 leading monomial = " << p1.leading_monomial().to_string(varnames) << std::endl;
    std::cout << "p1 leading coefficient = " << p1.leading_coefficient() << std::endl;
    
    std::cout << "\nAll tests passed!\n";
    return 0;
}