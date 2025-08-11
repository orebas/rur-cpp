#include <stdio.h>
#include <stdlib.h>
#include "axf4.h"

int main() {
    // Test x^2 + y^2 - 1 with F4 directly
    const char* vars[] = {"x", "y"};
    
    // Test with various primes
    ModularCoeff primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 1073741827};
    
    for (int i = 0; i < 11; i++) {
        ModularCoeff prime = primes[i];
        printf("\nTesting with prime %u:\n", prime);
        
        axf4_session_t session = axf4_create_session(prime, vars, 2);
        if (!session) {
            printf("  Failed to create session\n");
            continue;
        }
        
        // Add polynomial x^2 + y^2 - 1
        // Need to format it correctly for F4
        if (axf4_add_polynomial(session, "x^2 + y^2 - 1") != 0) {
            printf("  Failed to add polynomial\n");
            axf4_destroy_session(session);
            continue;
        }
        
        // Compute Groebner basis
        axf4_result_t result = axf4_compute_groebner_basis_keep_data(session);
        if (result.status != 0) {
            printf("  Failed to compute GB\n");
        } else {
            int basis_size = axf4_get_basis_size();
            printf("  GB has %d polynomial(s)\n", basis_size);
            
            // Get leading monomials
            if (basis_size > 0) {
                unsigned int* leading_monomials = malloc(basis_size * sizeof(unsigned int));
                int num_leading = axf4_get_all_leading_monomials(leading_monomials);
                printf("  Got %d leading monomial(s)\n", num_leading);
                
                // Decode and print first leading monomial
                if (num_leading > 0) {
                    unsigned int encoded = leading_monomials[0];
                    printf("  First leading monomial (encoded): %u\n", encoded);
                    // The encoding is: bits 0-14 for first var, 15-29 for second var
                    int exp_x = encoded & 0x7FFF;
                    int exp_y = (encoded >> 15) & 0x7FFF;
                    printf("  Decoded as: x^%d * y^%d\n", exp_x, exp_y);
                }
                free(leading_monomials);
            }
            
            // Try to get the actual polynomial
            char* poly_str = axf4_get_basis_polynomial(0);
            if (poly_str) {
                printf("  First polynomial: %s\n", poly_str);
                free(poly_str);
            }
        }
        
        axf4_free_result(&result);
        axf4_cleanup_basis_data();
        axf4_destroy_session(session);
    }
    
    return 0;
}