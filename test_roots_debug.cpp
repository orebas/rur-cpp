#include <iostream>
#include <vector>
#include <gmpxx.h>
#include "julia_rur/numerical_roots_eigen.hpp"

int main() {
    // Test polynomial T (coefficients 0, 1), which is T = 0 + 1*T
    std::vector<mpq_class> coeffs = {mpq_class(0), mpq_class(1)};
    
    std::cout << "Testing polynomial with coefficients: ";
    for (const auto& c : coeffs) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    auto roots = julia_rur::find_polynomial_roots(coeffs);
    
    std::cout << "Found " << roots.size() << " roots:" << std::endl;
    for (const auto& root : roots) {
        std::cout << "  " << root << std::endl;
    }
    
    // Also test x - 4 (coefficients -4, 1)
    std::vector<mpq_class> coeffs2 = {mpq_class(-4), mpq_class(1)};
    
    std::cout << "\nTesting polynomial with coefficients: ";
    for (const auto& c : coeffs2) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    auto roots2 = julia_rur::find_polynomial_roots(coeffs2);
    
    std::cout << "Found " << roots2.size() << " roots:" << std::endl;
    for (const auto& root : roots2) {
        std::cout << "  " << root << std::endl;
    }
    
    return 0;
}