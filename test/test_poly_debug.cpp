#include "polynomial.hpp"
#include <iostream>

int main() {
    std::cout << "Testing polynomial creation and display...\n";
    
    // Create x - 1
    MultivariatePolynomial<int> p;
    p.set_coefficient(Monomial({1}), 1);   // x
    p.set_coefficient(Monomial({0}), -1);  // -1
    
    std::vector<std::string> vars = {"x"};
    
    std::cout << "Polynomial has " << p.nterms() << " terms\n";
    std::cout << "To string: " << p.to_string(vars) << std::endl;
    
    // Test monomial creation
    Monomial m1({1});  // Should be x
    Monomial m2({0});  // Should be constant
    
    std::cout << "\nTesting monomials:\n";
    std::cout << "m1 ({1}): degree=" << m1.degree() << ", exponents=[" << m1[0] << "], string=" << m1.to_string(vars) << std::endl;
    std::cout << "m2 ({0}): degree=" << m2.degree() << ", exponents=[" << m2[0] << "], string=" << m2.to_string(vars) << std::endl;
    
    // Check monomials directly
    std::cout << "\nMonomials in polynomial:\n";
    for (const auto& [mon, coeff] : p) {
        std::cout << "Monomial degree: " << mon.degree() << ", exponents=[" << mon[0] << "], coeff: " << coeff << std::endl;
        std::cout << "Monomial to_string: " << mon.to_string(vars) << std::endl;
    }
    
    return 0;
}