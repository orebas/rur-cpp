#include "../src/julia_rur/data_structures.hpp"
#include "../src/julia_rur/f4_integration.hpp"
#include "../src/julia_rur/f4_polynomial_formatter.hpp"
#include "../src/julia_rur/quotient_basis.hpp"
#include <iostream>

using namespace julia_rur;

void test_gb_for_prime(ModularCoeff prime) {
    std::cout << "\n=== Testing GB computation for prime " << prime << " ===" << std::endl;
    
    // System: x^2 - 1 = 0, y - x = 0
    std::vector<std::string> polynomials = {
        "1*x^2-1",
        "1*y-1*x"
    };
    std::vector<std::string> variables = {"x", "y"};
    
    // Format polynomials for F4
    auto f4_input = format_polynomials_for_f4(polynomials, variables);
    if (!f4_input.success) {
        std::cerr << "Failed to format polynomials: " << f4_input.error_message << std::endl;
        return;
    }
    
    // Compute Gröbner basis
    F4Result result = compute_groebner_basis_f4(
        f4_input.polynomial_strings,
        f4_input.variable_names,
        prime,
        true  // verbose
    );
    
    if (!result.success) {
        std::cerr << "F4 failed: " << result.error_message << std::endl;
        return;
    }
    
    std::cout << "\nGröbner basis has " << result.polynomial_strings.size() << " polynomials:" << std::endl;
    for (size_t i = 0; i < result.polynomial_strings.size(); ++i) {
        std::cout << "  " << result.polynomial_strings[i] << std::endl;
        
        // Show the coefficient details
        if (i < result.polynomial_coefficients.size()) {
            std::cout << "    Coefficients: [";
            for (size_t j = 0; j < result.polynomial_coefficients[i].size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << result.polynomial_coefficients[i][j];
                
                // Show modular interpretation
                if (result.polynomial_coefficients[i][j] > prime/2) {
                    int64_t signed_val = static_cast<int64_t>(result.polynomial_coefficients[i][j]) - prime;
                    std::cout << " (" << signed_val << ")";
                }
            }
            std::cout << "]" << std::endl;
        }
    }
}

int main() {
    std::cout << "Testing F4 Gröbner basis consistency across primes" << std::endl;
    
    // Test with different primes
    test_gb_for_prime(1073741827);  // Large prime where it works
    test_gb_for_prime(1048573);      // Smaller prime where we get wrong result
    test_gb_for_prime(524287);       // Another small prime
    test_gb_for_prime(268435399);    // 28-bit prime
    
    return 0;
}