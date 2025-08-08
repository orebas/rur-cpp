#include <iostream>
#include "../src/axf4_wrapper.h"

// From axf4_lib.c - this controls the number of threads
extern "C" {
    extern int num_thread;
}

int main() {
    std::cout << "Testing F4 with single thread to avoid crashes\n\n";
    
    // Force single-threaded execution
    num_thread = 1;
    std::cout << "Set num_thread = 1\n";
    
    // Test the problematic case
    const char* vars[] = {"x"};
    axf4_session_t session = axf4_create_session(1073741827, vars, 1);
    
    if (!session) {
        std::cerr << "Failed to create session" << std::endl;
        return 1;
    }
    
    if (axf4_add_polynomial(session, "x - 3") != 0) {
        std::cerr << "Failed to add polynomial" << std::endl;
        axf4_destroy_session(session);
        return 1;
    }
    
    std::cout << "Computing GB with prime 1073741827 for x - 3..." << std::endl;
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    if (result.status == 0) {
        std::cout << "Success! GB: " << result.groebner_basis << std::endl;
    } else {
        std::cout << "Failed with status: " << result.status << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    // Test with different primes
    std::cout << "\nTesting with smaller prime 131063..." << std::endl;
    session = axf4_create_session(131063, vars, 1);
    
    if (!session) {
        std::cerr << "Failed to create session" << std::endl;
        return 1;
    }
    
    if (axf4_add_polynomial(session, "x - 3") != 0) {
        std::cerr << "Failed to add polynomial" << std::endl;
        axf4_destroy_session(session);
        return 1;
    }
    
    result = axf4_compute_groebner_basis(session);
    
    if (result.status == 0) {
        std::cout << "Success! GB: " << result.groebner_basis << std::endl;
    } else {
        std::cout << "Failed with status: " << result.status << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    // Reset to use all threads
    num_thread = 0;
    std::cout << "\nReset num_thread = 0 (use all cores)" << std::endl;
    
    return 0;
}