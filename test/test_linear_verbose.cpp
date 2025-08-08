#include <iostream>
#include "julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    std::cout << "Testing RUR for x - 3 = 0 with verbose output\n\n";
    
    std::vector<std::string> polynomials = {"x - 3"};
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = true;
    config.initial_prime_bits = 17;  // Use smaller primes like test_simple_rur
    
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);
    
    if (!rur_result.success) {
        std::cerr << "RUR computation failed: " << rur_result.error_message << std::endl;
        return 1;
    }
    
    std::cout << "\n=== Final RUR Result ===\n";
    std::cout << "Minimal polynomial: ";
    for (const auto& coeff : rur_result.minimal_polynomial) {
        std::cout << coeff << " ";
    }
    std::cout << "\n";
    
    std::cout << "Number of numerators: " << rur_result.numerators.size() << "\n";
    for (size_t i = 0; i < rur_result.numerators.size(); ++i) {
        std::cout << "Numerator " << i << ": ";
        if (rur_result.numerators[i].empty()) {
            std::cout << "[empty]";
        } else {
            for (const auto& coeff : rur_result.numerators[i]) {
                std::cout << coeff << " ";
            }
        }
        std::cout << "\n";
    }
    
    return 0;
}