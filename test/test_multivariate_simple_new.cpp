#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <chrono>

using namespace julia_rur;

int main() {
    std::cout << "Testing rational RUR for simple multivariate system" << std::endl;
    
    // Simple system with rational solutions: x^2 - 1 = 0, y - x = 0
    // Solutions: (x,y) = (1,1) and (-1,-1)
    std::vector<std::string> polynomials = {
        "1*x^2-1",      // x^2 - 1 = 0
        "1*y-1*x"       // y - x = 0
    };
    std::vector<std::string> variables = {"x", "y"};
    
    RURConfig config;
    config.verbose = true;
    config.initial_prime_bits = 20;  // Use smaller primes
    
    std::cout << "\n=== Computing rational RUR ===\n" << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    auto rat_result = compute_rational_rur(polynomials, variables, config);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "\nTime taken: " << duration.count() << " ms" << std::endl;
    
    if (rat_result.success) {
        std::cout << "\n✓ Rational RUR computation succeeded" << std::endl;
        std::cout << "Quotient dimension: " << rat_result.quotient_basis.size() << std::endl;
        std::cout << "Minimal polynomial degree: " << rat_result.minimal_polynomial.size() - 1 << std::endl;
        
        std::cout << "\nMinimal polynomial: T(t) = ";
        for (size_t i = 0; i < rat_result.minimal_polynomial.size(); ++i) {
            if (i > 0 && rat_result.minimal_polynomial[i] >= 0) std::cout << " + ";
            std::cout << rat_result.minimal_polynomial[i];
            if (i > 0) std::cout << "*t";
            if (i > 1) std::cout << "^" << i;
        }
        std::cout << std::endl;
        
        std::cout << "\nParameterizations:" << std::endl;
        for (size_t i = 0; i < variables.size(); ++i) {
            std::cout << "  " << variables[i] << " = ";
            if (i < rat_result.numerators.size()) {
                std::cout << "(size=" << rat_result.numerators[i].size() << ") ";
                if (rat_result.numerators[i].empty()) {
                    std::cout << "[empty]";
                } else {
                    bool first = true;
                    for (size_t j = 0; j < rat_result.numerators[i].size(); ++j) {
                        if (rat_result.numerators[i][j] != 0) {
                            if (!first && rat_result.numerators[i][j] > 0) std::cout << " + ";
                            std::cout << rat_result.numerators[i][j];
                            if (j > 0) std::cout << "*t";
                            if (j > 1) std::cout << "^" << j;
                            first = false;
                        }
                    }
                    if (first) {
                        // All coefficients were zero
                        std::cout << "0";
                    }
                }
                std::cout << " / (dT/dt)";
            } else {
                std::cout << "[missing]";
            }
            std::cout << std::endl;
        }
        
        // Verify the minimal polynomial has roots ±1
        std::cout << "\nVerification:" << std::endl;
        std::cout << "Expected minimal polynomial: (t-1)(t+1) = t^2 - 1" << std::endl;
        std::cout << "Expected x parameterization: x = t" << std::endl;
        std::cout << "Expected y parameterization: y = t" << std::endl;
    } else {
        std::cout << "\n✗ Rational RUR computation failed" << std::endl;
        std::cout << "Error: " << rat_result.error_message << std::endl;
    }
    
    return 0;
}