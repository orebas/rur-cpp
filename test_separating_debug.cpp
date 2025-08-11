#include <iostream>
#include <vector>
#include "src/julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    std::cout << "Debugging separating element search for Katsura-4\n";
    std::cout << "==================================================\n\n";
    
    // Katsura-4 polynomials
    std::vector<std::string> polynomials = {
        "x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1",
        "x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1",
        "x0*x2 + 2*x1*x3 + 2*x2*x4 - x2", 
        "x0*x3 + 2*x1*x4 - x3",
        "x0*x4 - x4"
    };
    
    std::vector<std::string> variables = {"x0", "x1", "x2", "x3", "x4"};
    
    // Use a known good prime
    ModularCoeff prime = 1073741827; // 30-bit prime
    
    RURConfig config;
    config.verbose = false;  // Less noise for this test
    
    std::cout << "Testing various separating elements with prime " << prime << "\n\n";
    
    // Test 1: Try x4 (last variable - what our code tries first)
    std::cout << "1. Testing x4 as separating element:\n";
    std::vector<int> x4_coeffs = {0, 0, 0, 0, 1};
    auto [result1, sep1] = compute_modular_rur(polynomials, variables, prime, config, x4_coeffs);
    if (result1.success) {
        std::cout << "   Quotient dimension: " << result1.quotient_basis.size() << "\n";
        std::cout << "   Minimal poly degree: " << result1.minimal_polynomial.degree << "\n";
        std::cout << "   Is separating? " << (result1.minimal_polynomial.degree == 16 ? "YES" : "NO") << "\n";
    } else {
        std::cout << "   Failed to compute RUR\n";
    }
    
    // Test 2: Try x0 (first variable)
    std::cout << "\n2. Testing x0 as separating element:\n";
    std::vector<int> x0_coeffs = {1, 0, 0, 0, 0};
    auto [result2, sep2] = compute_modular_rur(polynomials, variables, prime, config, x0_coeffs);
    if (result2.success) {
        std::cout << "   Quotient dimension: " << result2.quotient_basis.size() << "\n";
        std::cout << "   Minimal poly degree: " << result2.minimal_polynomial.degree << "\n";
        std::cout << "   Is separating? " << (result2.minimal_polynomial.degree == 16 ? "YES" : "NO") << "\n";
    }
    
    // Test 3: Try x4 - x3 (what Julia tries after single variables fail)
    std::cout << "\n3. Testing x4 - x3 as separating element:\n";
    std::vector<int> x4_minus_x3 = {0, 0, 0, -1, 1};
    auto [result3, sep3] = compute_modular_rur(polynomials, variables, prime, config, x4_minus_x3);
    if (result3.success) {
        std::cout << "   Quotient dimension: " << result3.quotient_basis.size() << "\n";
        std::cout << "   Minimal poly degree: " << result3.minimal_polynomial.degree << "\n";
        std::cout << "   Is separating? " << (result3.minimal_polynomial.degree == 16 ? "YES" : "NO") << "\n";
    }
    
    // Test 4: Try a simple linear combination
    std::cout << "\n4. Testing x0 + 2*x1 + 3*x2 + 5*x3 + 7*x4:\n";
    std::vector<int> linear_combo = {1, 2, 3, 5, 7};
    auto [result4, sep4] = compute_modular_rur(polynomials, variables, prime, config, linear_combo);
    if (result4.success) {
        std::cout << "   Quotient dimension: " << result4.quotient_basis.size() << "\n";
        std::cout << "   Minimal poly degree: " << result4.minimal_polynomial.degree << "\n";
        std::cout << "   Is separating? " << (result4.minimal_polynomial.degree == 16 ? "YES" : "NO") << "\n";
    }
    
    // Test 5: Let the algorithm find its own
    std::cout << "\n5. Letting algorithm choose separating element:\n";
    auto [result5, sep5] = compute_modular_rur(polynomials, variables, prime, config, {});
    if (result5.success) {
        std::cout << "   Quotient dimension: " << result5.quotient_basis.size() << "\n";
        std::cout << "   Minimal poly degree: " << result5.minimal_polynomial.degree << "\n";
        std::cout << "   Is separating? " << (result5.minimal_polynomial.degree == 16 ? "YES" : "NO") << "\n";
        if (!sep5.empty()) {
            std::cout << "   Found coefficients: [";
            for (size_t i = 0; i < sep5.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << sep5[i];
            }
            std::cout << "]\n";
        }
    } else {
        std::cout << "   FAILED - algorithm could not find separating element\n";
    }
    
    return 0;
}