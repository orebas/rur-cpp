#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <chrono>
#include <thread>
#include <atomic>
#include <signal.h>

extern "C" {
#include "src/axf4_wrapper.h"
}

// Timeout mechanism to detect infinite loops
std::atomic<bool> computation_complete(false);
std::atomic<bool> timeout_occurred(false);

void timeout_handler(int signum) {
    if (!computation_complete.load()) {
        timeout_occurred.store(true);
        std::cout << "\nTIMEOUT: Computation exceeded time limit!" << std::endl;
        exit(1);
    }
}

void test_prime(int prime, int timeout_seconds = 10) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Testing F4 with prime = " << prime << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Reset atomic flags
    computation_complete.store(false);
    timeout_occurred.store(false);
    
    // Setup timeout
    signal(SIGALRM, timeout_handler);
    alarm(timeout_seconds);
    
    // Create F4 session
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(prime, vars, 2);
    
    if (!session) {
        std::cout << "ERROR: Failed to create F4 session" << std::endl;
        return;
    }
    
    // Add polynomials: circle (x^2 + y^2 - 1) and line (x - y)
    std::cout << "Input polynomials:" << std::endl;
    std::cout << "  f1 = x^2 + y^2 - 1" << std::endl;
    std::cout << "  f2 = x - y" << std::endl;
    
    // F4 expects polynomials in a specific format
    if (axf4_add_polynomial(session, "x^2+y^2-1") != 0) {
        std::cout << "ERROR: Failed to add polynomial 1" << std::endl;
        axf4_destroy_session(session);
        return;
    }
    
    if (axf4_add_polynomial(session, "x-y") != 0) {
        std::cout << "ERROR: Failed to add polynomial 2" << std::endl;
        axf4_destroy_session(session);
        return;
    }
    
    // Compute Gröbner basis
    std::cout << "\nComputing Gröbner basis..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    // Cancel alarm - computation completed
    computation_complete.store(true);
    alarm(0);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    if (result.status == 0) {
        std::cout << "\nSUCCESS! Gröbner basis computed in " 
                  << duration.count() << " ms" << std::endl;
        std::cout << "Computation time (F4 internal): " 
                  << result.computation_time << " seconds" << std::endl;
        std::cout << "Basis size: " << result.basis_size << std::endl;
        std::cout << "\nGröbner basis:" << std::endl;
        std::cout << result.groebner_basis << std::endl;
    } else {
        std::cout << "\nERROR: Gröbner basis computation failed!" << std::endl;
        if (result.error_message) {
            std::cout << "Error message: " << result.error_message << std::endl;
        }
    }
    
    // Cleanup
    axf4_free_result(&result);
    axf4_destroy_session(session);
}

int main(int argc, char* argv[]) {
    std::cout << "F4 Infinite Loop Test Case" << std::endl;
    std::cout << "Testing Gröbner basis computation for system:" << std::endl;
    std::cout << "  x^2 + y^2 - 1 = 0  (circle)" << std::endl;
    std::cout << "  x - y = 0          (line)" << std::endl;
    
    // Test with known working prime
    test_prime(1073741827, 30);  // Works correctly
    
    // Test with problematic primes that cause infinite loops
    std::vector<int> problematic_primes = {
        2147483629,  // First reported problematic prime
        2147483587,
        2147483579,
        2147483563,
        2147483549
    };
    
    for (int prime : problematic_primes) {
        test_prime(prime, 10);  // 10 second timeout for infinite loop detection
    }
    
    return 0;
}