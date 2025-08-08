#include <iostream>
#include <gmpxx.h>
#include "../src/julia_rur/rur_main_algorithm.hpp"

int main() {
    std::cout << "Testing rational coefficient representation..." << std::endl;
    
    std::vector<std::string> polynomials = {"x^2 - 2"};
    std::vector<std::string> variables = {"x"};
    
    julia_rur::RURConfig config;
    config.verbose = false;
    
    julia_rur::RationalRURResult result = julia_rur::compute_rational_rur(
        polynomials, variables, config
    );
    
    if (result.success && !result.minimal_polynomial.empty()) {
        std::cout << "\nMinimal polynomial has " << result.minimal_polynomial.size() << " coefficients" << std::endl;
        
        for (size_t i = 0; i < result.minimal_polynomial.size(); ++i) {
            const mpq_class& coeff = result.minimal_polynomial[i];
            
            std::cout << "\nCoeff[" << i << "]:" << std::endl;
            std::cout << "  Direct print: " << coeff << std::endl;
            std::cout << "  Numerator: " << coeff.get_num() << std::endl;
            std::cout << "  Denominator: " << coeff.get_den() << std::endl;
            std::cout << "  get_d(): " << coeff.get_d() << std::endl;
            
            // Check sign
            if (coeff < 0) {
                std::cout << "  Sign: negative" << std::endl;
            } else if (coeff > 0) {
                std::cout << "  Sign: positive" << std::endl;
            } else {
                std::cout << "  Sign: zero" << std::endl;
            }
        }
    }
    
    return 0;
}