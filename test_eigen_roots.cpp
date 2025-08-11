#include <iostream>
#include <vector>
#include <complex>
#include <gmpxx.h>
#include "../src/julia_rur/numerical_roots_eigen.hpp"

using namespace julia_rur;

int main() {
    // Test: T - 3 = 0
    std::vector<mpq_class> coeffs = {mpq_class(-3), mpq_class(1)};
    
    std::cout << "Testing polynomial: ";
    for (int i = coeffs.size() - 1; i >= 0; --i) {
        if (i < coeffs.size() - 1 && coeffs[i] > 0) std::cout << " + ";
        std::cout << coeffs[i];
        if (i > 0) std::cout << "*T^" << i;
    }
    std::cout << " = 0" << std::endl;
    
    auto roots = find_polynomial_roots(coeffs);
    
    std::cout << "Found " << roots.size() << " roots:" << std::endl;
    for (const auto& root : roots) {
        std::cout << "  " << root << std::endl;
    }
    
    return 0;
}