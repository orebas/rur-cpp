#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <chrono>
#include <signal.h>

using namespace julia_rur;

volatile bool timeout_occurred = false;

void timeout_handler(int sig) {
    timeout_occurred = true;
    std::cerr << "\nTIMEOUT: Rational RUR computation taking too long!" << std::endl;
    exit(1);
}

int main() {
    std::cout << "Testing rational RUR for multivariate system" << std::endl;
    
    // Simple system: x^2 + y^2 - 1 = 0, x - y = 0
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2-1",  // Circle
        "1*x-1*y"         // Line x = y
    };
    std::vector<std::string> variables = {"x", "y"};
    
    RURConfig config;
    config.verbose = true;
    config.initial_prime_bits = 28;  // Use 28-bit primes like Julia does by default
    // config.max_primes = 5;  // Field doesn't exist
    
    // Set timeout
    signal(SIGALRM, timeout_handler);
    alarm(30);  // 30 second timeout
    
    std::cout << "\n=== Computing rational RUR ===" << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    try {
        auto rat_result = compute_rational_rur(polynomials, variables, config);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "\nTime taken: " << duration.count() << " ms" << std::endl;
        
        if (rat_result.success) {
            std::cout << "\n✓ Rational RUR computation succeeded" << std::endl;
            std::cout << "Quotient dimension: " << rat_result.quotient_basis.size() << std::endl;
            std::cout << "Minimal polynomial degree: " << rat_result.minimal_polynomial.size() - 1 << std::endl;
            
            std::cout << "\nMinimal polynomial coefficients: ";
            for (const auto& c : rat_result.minimal_polynomial) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
            
            std::cout << "\nParameterizations:" << std::endl;
            for (size_t i = 0; i < variables.size(); ++i) {
                std::cout << "  Variable " << variables[i] << ": ";
                if (i < rat_result.numerators.size()) {
                    for (const auto& c : rat_result.numerators[i]) {
                        std::cout << c << " ";
                    }
                } else {
                    std::cout << "(missing)";
                }
                std::cout << std::endl;
            }
        } else {
            std::cout << "\n✗ Rational RUR computation failed" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "\nException: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}