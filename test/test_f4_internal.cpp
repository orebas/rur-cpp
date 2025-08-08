#include <iostream>
#include <cstdio>
extern "C" {
#include "../src/axf4_wrapper.h"
}

void test_system_internal(unsigned int prime) {
    std::cout << "\n=== Testing with prime " << prime << " ===" << std::endl;
    
    // Create F4 session with variables x, y
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(prime, vars, 2);
    
    // Add polynomial x^2 - 1
    std::cout << "Adding: x^2 - 1" << std::endl;
    axf4_add_polynomial(session, "1*x^2-1");
    
    // Add polynomial y - x  
    std::cout << "Adding: y - x" << std::endl;
    axf4_add_polynomial(session, "1*y-1*x");
    
    // Compute GB but keep data
    axf4_result_t result = axf4_compute_groebner_basis_keep_data(session);
    
    std::cout << "GB string output:\n" << result.groebner_basis << std::endl;
    
    // Now examine the internal data
    int basis_size = axf4_get_basis_size();
    std::cout << "\nInternal representation:" << std::endl;
    std::cout << "Basis size: " << basis_size << std::endl;
    
    for (int i = 0; i < basis_size; i++) {
        int term_count = axf4_get_poly_term_count(i);
        std::cout << "\nPolynomial " << i << " has " << term_count << " terms:" << std::endl;
        
        if (term_count > 0) {
            unsigned int* coeffs = new unsigned int[term_count];
            unsigned int* monomials = new unsigned int[term_count];
            
            axf4_get_poly_data(i, coeffs, monomials);
            
            for (int j = 0; j < term_count; j++) {
                std::cout << "  Term " << j << ": coeff=" << coeffs[j];
                
                // Check if this looks like a negative value
                if (coeffs[j] > prime / 2) {
                    std::cout << " (could be -" << (prime - coeffs[j]) << ")";
                }
                
                // Check if coefficient is exactly 1
                if (coeffs[j] == 1) {
                    std::cout << " <-- This is literally 1, not -1!";
                }
                
                std::cout << ", monomial=" << monomials[j] << std::endl;
            }
            
            delete[] coeffs;
            delete[] monomials;
        }
    }
    
    axf4_cleanup_basis_data();
    axf4_free_result(&result);
    axf4_destroy_session(session);
}

int main() {
    std::cout << "F4 Internal representation test" << std::endl;
    
    // Test with different primes
    test_system_internal(1073741827);  // Working prime
    test_system_internal(1048573);      // Failing prime
    
    return 0;
}