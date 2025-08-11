#include <iostream>
#include <cstdio>
#include <cstring>

extern "C" {
#include "src/axf4_wrapper.h"
}

void test_system_with_prime(int prime) {
    std::cout << "\n=== Testing with prime " << prime << " ===" << std::endl;
    
    // Redirect stdout to capture F4 output
    char buffer[65536];
    memset(buffer, 0, sizeof(buffer));
    FILE* memstream = fmemopen(buffer, sizeof(buffer)-1, "w");
    FILE* old_stdout = stdout;
    stdout = memstream;
    
    // Create F4 session
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(prime, vars, 2);
    
    if (!session) {
        stdout = old_stdout;
        std::cout << "Failed to create session" << std::endl;
        return;
    }
    
    // Add simple system: x^2 + y^2 - 1, x - y
    axf4_add_polynomial(session, "x^2+y^2-1");
    axf4_add_polynomial(session, "x-y");
    
    // Compute with timeout using alarm
    alarm(5);  // 5 second timeout
    
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    alarm(0);  // Cancel alarm
    
    // Restore stdout and close memstream
    fflush(memstream);
    stdout = old_stdout;
    fclose(memstream);
    
    // Print captured output
    std::cout << "F4 Output:" << std::endl;
    std::cout << buffer << std::endl;
    
    if (result.status == 0) {
        std::cout << "SUCCESS: Basis size = " << result.basis_size << std::endl;
        std::cout << "Result:\n" << result.groebner_basis << std::endl;
    } else {
        std::cout << "FAILED: " << (result.error_message ? result.error_message : "Unknown error") << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
}

int main() {
    // Test a working prime
    test_system_with_prime(1073741827);
    
    // Test a problematic prime
    test_system_with_prime(2147483629);
    
    return 0;
}