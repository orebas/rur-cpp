#include "../src/axf4_wrapper.h"
#include <iostream>
#include <vector>
#include <cassert>

int main() {
    std::cout << "Debug: Testing structured API by copying working wrapper approach..." << std::endl;
    
    // Use the exact same approach as test_f4_wrapper.cpp but add structured calls
    // This mimics the working F4Solver but calls structured API
    
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(100003, vars, 2);
    if (!session) {
        std::cerr << "Failed to create session" << std::endl;
        return 1;
    }
    
    std::cout << "Debug: Adding polynomials like working test..." << std::endl;
    // Add the same polynomials as the working test: x^2 + y^2 - 1, x - y
    if (axf4_add_polynomial(session, "x^2 + y^2 - 1") != 0) {
        std::cerr << "Failed to add first polynomial" << std::endl;
        axf4_destroy_session(session);
        return 1;
    }
    
    if (axf4_add_polynomial(session, "x - y") != 0) {
        std::cerr << "Failed to add second polynomial" << std::endl;
        axf4_destroy_session(session);
        return 1;
    }
    
    std::cout << "Debug: Computing basis..." << std::endl;
    // Use the standard function that we know works
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    if (result.status != 0) {
        std::cerr << "Computation failed: " << result.error_message << std::endl;
        axf4_free_result(&result);
        axf4_destroy_session(session);
        return 1;
    }
    
    std::cout << "Debug: Computation succeeded!" << std::endl;
    std::cout << "Basis size: " << result.basis_size << std::endl;
    
    // Now test if structured API works after the standard computation
    // (This should fail since the data was cleaned up)
    std::cout << "Debug: Testing structured API after cleanup..." << std::endl;
    int api_basis_size = axf4_get_basis_size();
    std::cout << "Structured API basis size: " << api_basis_size << std::endl;
    
    if (api_basis_size != result.basis_size) {
        std::cout << "As expected, structured API doesn't work after cleanup" << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    std::cout << "✓ Debug test completed - confirmed structured API needs data kept alive" << std::endl;
    return 0;
}