#include <iostream>
#include <vector>
#include "src/julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    std::cout << "Testing Katsura-4 with RUR implementation\n";
    std::cout << "==========================================\n\n";
    
    // Katsura-4 polynomials
    std::vector<std::string> polynomials = {
        "x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1",
        "x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1",
        "x0*x2 + 2*x1*x3 + 2*x2*x4 - x2", 
        "x0*x3 + 2*x1*x4 - x3",
        "x0*x4 - x4"
    };
    
    std::vector<std::string> variables = {"x0", "x1", "x2", "x3", "x4"};
    
    // Configure for verbose output
    RURConfig config;
    config.verbose = true;
    config.timing = true;
    
    // Try computing with a single prime first
    ModularCoeff prime = 1073741827; // 30-bit prime
    std::cout << "Computing modular RUR with prime " << prime << "\n";
    
    auto [result, sep_coeffs] = compute_modular_rur(polynomials, variables, prime, config, {});
    
    if (result.success) {
        std::cout << "\nSUCCESS!\n";
        std::cout << "Quotient basis dimension: " << result.quotient_basis.size() << "\n";
        std::cout << "Minimal polynomial degree: " << result.minimal_polynomial.degree << "\n";
        std::cout << "Separating element coefficients: [";
        for (size_t i = 0; i < sep_coeffs.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << sep_coeffs[i];
        }
        std::cout << "]\n";
    } else {
        std::cout << "\nFAILED!\n";
        std::cout << "Quotient basis dimension: " << result.quotient_basis.size() << "\n";
        if (result.minimal_polynomial.success) {
            std::cout << "Minimal polynomial was found (degree " << result.minimal_polynomial.degree << ")\n";
        } else {
            std::cout << "Failed to find minimal polynomial\n";
        }
        std::cout << "Number of parameterizations: " << result.parameterizations.size() << "\n";
    }
    
    return result.success ? 0 : 1;
}