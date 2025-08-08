#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    // Test x^2 - 2 with verbose output
    std::cout << "Testing x^2 - 2 to debug sign issue" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x^2-2"};
    std::vector<std::string> variables = {"x"};
    
    // First test with single prime
    ModularCoeff prime = 1073741827;
    RURConfig config;
    config.verbose = true;
    
    std::cout << "\n=== Testing modulo " << prime << " ===" << std::endl;
    auto mod_result = compute_modular_rur(polynomials, variables, prime, config);
    
    if (mod_result.success) {
        std::cout << "\nMinimal polynomial coefficients mod " << prime << ": ";
        for (const auto& c : mod_result.minimal_polynomial.coefficients) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        
        // Convert to symmetric range to see actual values
        std::cout << "In symmetric range: ";
        for (const auto& c : mod_result.minimal_polynomial.coefficients) {
            ModularCoeff sym = c;
            if (sym > prime / 2) sym = sym - prime;
            std::cout << (int32_t)sym << " ";
        }
        std::cout << std::endl;
    }
    
    // Now test rational reconstruction with just a few primes
    std::cout << "\n=== Testing rational reconstruction ===" << std::endl;
    config.verbose = false;
    config.initial_prime_bits = 31;
    
    auto rat_result = compute_rational_rur(polynomials, variables, config);
    
    if (rat_result.success) {
        std::cout << "\nRational minimal polynomial coefficients: ";
        for (const auto& c : rat_result.minimal_polynomial) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
    }
    
    return 0;
}