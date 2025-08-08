#include <iostream>
#include <string>
#include <vector>
#include "axf4_wrapper.h"

int main() {
    std::cout << "Testing F4 parsing of negative constants..." << std::endl;
    
    // Test prime
    unsigned long prime = 131071;
    
    // Create session
    const char* vars[] = {"x"};
    axf4_session_t session = axf4_create_session(prime, vars, 1);
    if (!session) {
        std::cerr << "Failed to create session" << std::endl;
        return 1;
    }
    
    // Test various formats
    std::vector<std::string> test_cases = {
        "1*x^2-2",       // Direct negative
        "1*x^2+-2",      // Plus-minus (current formatter output)
        "1*x^2+131069",  // Modular representation of -2 mod 131071
        "x^2+131069",    // Without coefficient
        "1*x^2+(-2)",    // Parentheses
        "x-3",           // Simple case
        "1*x-3",         // With coefficient
        "1*x+131068"     // Modular representation of -3 mod 131071
    };
    
    for (const auto& poly : test_cases) {
        std::cout << "\nTesting: \"" << poly << "\"" << std::endl;
        
        int result = axf4_add_polynomial(session, poly.c_str());
        if (result == 0) {
            std::cout << "  Success! Polynomial added." << std::endl;
            
            // Compute GB to see result
            axf4_result_t gb = axf4_compute_groebner_basis(session);
            if (gb.status == 0) {
                std::cout << "  GB: " << gb.groebner_basis << std::endl;
            }
            axf4_free_result(&gb);
            
            // Clear for next test
            axf4_cleanup_basis_data();
            axf4_destroy_session(session);
            session = axf4_create_session(prime, vars, 1);
        } else {
            const char* error = axf4_get_last_error(session);
            std::cout << "  Failed: " << (error ? error : "unknown error") << std::endl;
        }
    }
    
    axf4_destroy_session(session);
    return 0;
}