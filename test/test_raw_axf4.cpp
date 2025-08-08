#include "axf4_wrapper.h"
#include <iostream>
#include <cstring>

int main() {
    // Test simple system: x^2 - 1, y^2 - 1
    const char* vars[] = {"x", "y"};
    
    axf4_session_t session = axf4_create_session(100003, vars, 2);
    if (!session) {
        std::cout << "Failed to create session" << std::endl;
        return 1;
    }
    
    // Add polynomials in axf4 format
    if (axf4_add_polynomial(session, "1*x^2+100002") != 0) {  // x^2 - 1 = x^2 + (p-1)
        std::cout << "Failed to add polynomial 1" << std::endl;
    }
    
    if (axf4_add_polynomial(session, "1*y^2+100002") != 0) {  // y^2 - 1 = y^2 + (p-1)
        std::cout << "Failed to add polynomial 2" << std::endl;
    }
    
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    if (result.status == 0) {
        std::cout << "Raw axf4 output:" << std::endl;
        std::cout << "----------------" << std::endl;
        std::cout << result.groebner_basis << std::endl;
        std::cout << "----------------" << std::endl;
        std::cout << "Basis size: " << result.basis_size << std::endl;
    } else {
        std::cout << "Error: " << (result.error_message ? result.error_message : "unknown") << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    return 0;
}