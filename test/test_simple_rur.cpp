#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    std::cout << "Testing RUR for x^2 - 2" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x^2-2"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = true;
    config.initial_prime_bits = 17;
    
    RationalRURResult result = compute_rational_rur(polynomials, variables, config);
    
    if (result.success) {
        std::cout << "\nSUCCESS! Minimal polynomial coefficients:";
        for (const auto& c : result.minimal_polynomial) {
            std::cout << " " << c;
        }
        std::cout << std::endl;
    } else {
        std::cout << "\nFAILED: " << result.error_message << std::endl;
    }
    
    return result.success ? 0 : 1;
}