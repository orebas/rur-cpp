#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    std::cout << "Testing modular minimal polynomial computation" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x^2-2"};
    std::vector<std::string> variables = {"x"};
    
    ModularCoeff prime = 131063;
    RURConfig config;
    config.verbose = true;
    
    auto result = compute_modular_rur(polynomials, variables, prime, config);
    
    if (result.success) {
        std::cout << "\nMinimal polynomial coefficients:";
        for (auto c : result.minimal_polynomial.coefficients) {
            std::cout << " " << c;
            // Also show what it is modulo prime
            if (c > prime/2) {
                std::cout << "(" << static_cast<int64_t>(c) - prime << ")";
            }
        }
        std::cout << std::endl;
        std::cout << "Degree: " << result.minimal_polynomial.degree << std::endl;
    } else {
        std::cout << "Failed!" << std::endl;
    }
    
    return 0;
}