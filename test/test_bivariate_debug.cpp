#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    // Test x^2 - 2 with single prime
    std::cout << "Testing x^2 - 2 with single prime" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x^2-2"};
    std::vector<std::string> variables = {"x"};
    ModularCoeff prime = 1073741827;
    
    RURConfig config;
    config.verbose = true;
    
    auto result = compute_modular_rur(polynomials, variables, prime, config);
    
    if (result.success) {
        std::cout << "\nSuccess!" << std::endl;
        std::cout << "Minimal polynomial degree: " << result.minimal_polynomial.degree << std::endl;
        std::cout << "Minimal polynomial coefficients: ";
        for (const auto& c : result.minimal_polynomial.coefficients) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "\nFailed!" << std::endl;
    }
    
    return 0;
}