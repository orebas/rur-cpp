#include "../src/julia_rur/data_structures.hpp"
#include "../src/julia_rur/f4_polynomial_formatter.hpp"
extern "C" {
#include "../src/axf4_wrapper.h"
}
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
    
    // Create F4 session
    std::vector<const char*> var_ptrs;
    for (const auto& var : variables) {
        var_ptrs.push_back(var.c_str());
    }
    
    axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), variables.size());
    if (!session) {
        std::cerr << "Failed to create F4 session" << std::endl;
        return;
    }
    
    // Add polynomials
    for (const auto& poly : polynomials) {
        std::string formatted = format_polynomial_for_f4(poly);
        std::cout << "Adding polynomial: " << formatted << std::endl;
        if (axf4_add_polynomial(session, formatted.c_str()) != 0) {
            std::cerr << "Failed to add polynomial: " << poly << std::endl;
            axf4_destroy_session(session);
            return;
        }
    }
    
    // Compute Gröbner basis
    axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
    if (gb_result.status != 0) {
        std::cerr << "F4 failed with status " << gb_result.status << std::endl;
        axf4_free_result(&gb_result);
        axf4_destroy_session(session);
        return;
    }
    
    std::cout << "\nGröbner basis computed, size: " << gb_result.basis_size << std::endl;
    std::cout << "Gröbner basis string: " << gb_result.groebner_basis << std::endl;
    
    // Clean up
    axf4_free_result(&gb_result);
    axf4_destroy_session(session);
}

int main() {
    std::cout << "Testing F4 Gröbner basis consistency across primes" << std::endl;
    std::cout << "System: x^2 - 1 = 0, y - x = 0" << std::endl;
    std::cout << "Expected GB: should contain y^2 - 1 (or equivalent)" << std::endl;
    
    // Test with different primes
    test_gb_for_prime(1073741827);  // Large prime where it works
    test_gb_for_prime(1048573);      // Smaller prime where we get wrong result
    test_gb_for_prime(524287);       // Another small prime
    test_gb_for_prime(268435399);    // 28-bit prime where 2 is QR
    
    // Test the primes we just saw fail
    test_gb_for_prime(1048447);      // Should be OK but fails
    test_gb_for_prime(1048433);      // Should be OK but fails
    
    return 0;
}