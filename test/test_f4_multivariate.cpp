#include "../src/axf4_wrapper.h"
#include "../src/julia_rur/data_structures.hpp"
#include <iostream>
#include <chrono>
#include <signal.h>
#include <vector>
#include <string>

using namespace julia_rur;

volatile bool timeout_occurred = false;

void timeout_handler(int sig) {
    timeout_occurred = true;
    std::cerr << "\nTIMEOUT: F4 computation taking too long!" << std::endl;
    exit(1);
}

int main() {
    std::cout << "Testing F4 on multivariate system" << std::endl;
    
    // Simple system: x^2 + y^2 - 1 = 0, x - y = 0
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2-1",  // Circle
        "1*x-1*y"         // Line x = y
    };
    std::vector<std::string> variables = {"x", "y"};
    
    ModularCoeff prime = 65537;
    
    // Set timeout
    signal(SIGALRM, timeout_handler);
    alarm(10);  // 10 second timeout
    
    std::cout << "\n=== Testing F4 modulo " << prime << " ===" << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Create F4 session
    std::vector<const char*> var_ptrs;
    for (const auto& var : variables) {
        var_ptrs.push_back(var.c_str());
    }
    
    axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), variables.size());
    if (!session) {
        std::cerr << "Failed to create F4 session" << std::endl;
        return 1;
    }
    
    // Note: axf4_set_verbosity doesn't exist in the C interface
    
    // Add polynomials
    for (const auto& poly : polynomials) {
        if (axf4_add_polynomial(session, poly.c_str()) != 0) {
            std::cerr << "Failed to add polynomial: " << poly << std::endl;
            axf4_destroy_session(session);
            return 1;
        }
    }
    
    std::cout << "Computing Gröbner basis..." << std::endl;
    
    // Compute GB
    axf4_result_t gb_result = axf4_compute_groebner_basis(session);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "\nTime taken: " << duration.count() << " ms" << std::endl;
    
    if (gb_result.status == 0) {
        std::cout << "\n✓ F4 computation succeeded" << std::endl;
        std::cout << "GB size: " << gb_result.basis_size << std::endl;
        std::cout << "\nGröbner basis:" << std::endl;
        std::cout << gb_result.groebner_basis << std::endl;
    } else {
        std::cout << "\n✗ F4 computation failed with status: " << gb_result.status << std::endl;
    }
    
    axf4_free_result(&gb_result);
    axf4_destroy_session(session);
    
    // Cancel alarm
    alarm(0);
    
    return 0;
}