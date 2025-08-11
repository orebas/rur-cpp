#include <iostream>
#include <vector>
#include "src/julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    std::cout << "Testing if x4 is a separating element for Katsura-4\n";
    std::cout << "====================================================\n\n";
    
    // Katsura-4 polynomials
    std::vector<std::string> polynomials = {
        "x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1",
        "x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1",
        "x0*x2 + 2*x1*x3 + 2*x2*x4 - x2", 
        "x0*x3 + 2*x1*x4 - x3",
        "x0*x4 - x4"
    };
    
    std::vector<std::string> variables = {"x0", "x1", "x2", "x3", "x4"};
    
    // Try with a specific prime
    ModularCoeff prime = 1073741827; // 30-bit prime
    
    RURConfig config;
    config.verbose = true;
    
    std::cout << "Computing modular RUR with prime " << prime << "\n";
    std::cout << "Using x4 (variable 5) as separating element\n\n";
    
    // Force x4 as separating element: coefficients [0,0,0,0,1] means 0*x0 + 0*x1 + 0*x2 + 0*x3 + 1*x4
    std::vector<int> x4_coeffs = {0, 0, 0, 0, 1};
    
    auto [result, sep_coeffs] = compute_modular_rur(polynomials, variables, prime, config, x4_coeffs);
    
    if (result.success) {
        std::cout << "\n=== RESULT ===\n";
        std::cout << "SUCCESS: RUR computed\n";
        std::cout << "Quotient basis dimension: " << result.quotient_basis.size() << "\n";
        std::cout << "Minimal polynomial degree: " << result.minimal_polynomial.degree << "\n";
        
        if (result.minimal_polynomial.degree == 16) {
            std::cout << "\n✓ x4 IS a valid separating element (degree = 16)\n";
        } else {
            std::cout << "\n✗ x4 is NOT a separating element (degree = " 
                      << result.minimal_polynomial.degree << " ≠ 16)\n";
        }
    } else {
        std::cout << "\nFAILED to compute RUR\n";
        std::cout << "This could mean x4 is not a separating element\n";
    }
    
    return 0;
}