#include <iostream>
#include <vector>
#include "julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    std::cout << "=== Debug Coefficient Ordering Across Primes ===" << std::endl;
    
    // Test the specific failing system
    std::vector<std::string> polynomials = {"x^2 - 1", "y^2 - 2", "z^2 - 3"};
    std::vector<std::string> variables = {"x", "y", "z"};
    
    // Test with a few specific primes to see coefficient ordering differences
    std::vector<ModularCoeff> test_primes = {1073741827, 1073741783, 1073741681};
    
    for (size_t p_idx = 0; p_idx < test_primes.size(); ++p_idx) {
        ModularCoeff prime = test_primes[p_idx];
        std::cout << "\n========================================" << std::endl;
        std::cout << "Testing prime " << prime << " (index " << p_idx << ")" << std::endl;
        std::cout << "========================================" << std::endl;
        
        RURConfig config;
        config.verbose = false;  // Reduce noise
        
        try {
            auto mod_result = compute_modular_rur(polynomials, variables, prime, config);
            
            if (mod_result.success) {
                std::cout << "✓ SUCCESS for prime " << prime << std::endl;
                
                // Print minimal polynomial coefficients in detail
                std::cout << "Minimal polynomial coefficients (size=" << mod_result.minimal_polynomial.coefficients.size() << "):" << std::endl;
                for (size_t i = 0; i < mod_result.minimal_polynomial.coefficients.size(); ++i) {
                    std::cout << "  coeff[" << i << "] = " << mod_result.minimal_polynomial.coefficients[i];
                    if (i == 6) {
                        std::cout << " <-- THIS IS POSITION [0][6] CAUSING CRT FAILURE";
                    }
                    std::cout << std::endl;
                }
                
                // Print parameterizations
                for (size_t var_idx = 0; var_idx < mod_result.parameterizations.size(); ++var_idx) {
                    if (!mod_result.parameterizations[var_idx].generators.empty()) {
                        const auto& gen = mod_result.parameterizations[var_idx].generators[0];
                        std::cout << "Variable " << variables[var_idx] << " parameterization (size=" << gen.size() << "):" << std::endl;
                        for (size_t i = 0; i < std::min(size_t(10), gen.size()); ++i) {
                            std::cout << "  param[" << i << "] = " << gen[i] << std::endl;
                        }
                    }
                }
                
            } else {
                std::cout << "✗ FAILED for prime " << prime << std::endl;
            }
            
        } catch (const std::exception& e) {
            std::cout << "Exception for prime " << prime << ": " << e.what() << std::endl;
        }
    }
    
    std::cout << "\n=== Analysis ===" << std::endl;
    std::cout << "Look for differences in:" << std::endl;
    std::cout << "1. Size of minimal polynomial coefficient vectors" << std::endl;
    std::cout << "2. Value at position [0][6] across different primes" << std::endl;
    std::cout << "3. Ordering of coefficients (constant term first vs last)" << std::endl;
    std::cout << "\nIf coefficient [0][6] varies significantly across primes," << std::endl;
    std::cout << "then the canonicalization function is not working properly." << std::endl;
    
    return 0;
}