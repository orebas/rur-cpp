#include "../src/axf4_wrapper.h"
#include <iostream>

int main() {
    int primes[] = {131063, 131059, 100003};
    
    for (auto prime : primes) {
        std::cout << "\nPrime " << prime << ":" << std::endl;
        
        const char* vars[] = {"x"};
        axf4_session_t session = axf4_create_session(prime, vars, 1);
        
        // Test different ways of representing x^2 - 2
        const char* polys[] = {
            "1*x^2-2",
            "1*x^2+-2",
            "1*x^2+(-2)"
        };
        
        for (auto poly : polys) {
            axf4_session_t test_session = axf4_create_session(prime, vars, 1);
            int status = axf4_add_polynomial(test_session, poly);
            
            if (status == 0) {
                axf4_result_t result = axf4_compute_groebner_basis(test_session);
                std::cout << "  Input: " << poly << " -> GB: " << result.groebner_basis << std::endl;
                axf4_free_result(&result);
            } else {
                std::cout << "  Input: " << poly << " -> FAILED TO PARSE" << std::endl;
            }
            
            axf4_destroy_session(test_session);
        }
    }
    
    return 0;
}