#include <iostream>
#include <cstdio>
extern "C" {
#include "../src/axf4_wrapper.h"
}

void test_system(unsigned int prime) {
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
    
    // Compute GB
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    std::cout << "GB status: " << result.status << std::endl;
    std::cout << "GB size: " << result.basis_size << std::endl;
    std::cout << "GB string:\n" << result.groebner_basis << std::endl;
    
    // Expected: The GB should contain y^2 - 1
    // This would be represented as either:
    // - "y^2 + (p-1)" where p-1 represents -1 mod p
    // - "y^2 - 1" in some form
    
    std::cout << "\nAnalysis:" << std::endl;
    std::cout << "For x^2 - 1 = 0 and y - x = 0:" << std::endl;
    std::cout << "- We have y = x" << std::endl;
    std::cout << "- Substituting: y^2 = x^2 = 1" << std::endl;
    std::cout << "- So y^2 - 1 = 0 should be in the GB" << std::endl;
    std::cout << "- In F4 format with prime " << prime << ":" << std::endl;
    std::cout << "  - Correct: y^2 + " << (prime - 1) << " (representing y^2 - 1)" << std::endl;
    std::cout << "  - Wrong: y^2 + 1 (representing y^2 + 1)" << std::endl;
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
}

int main() {
    std::cout << "Direct F4 test - checking polynomial representation" << std::endl;
    
    // Test with a few primes
    test_system(1073741827);  // 30-bit prime
    test_system(1048573);      // 20-bit prime  
    test_system(524287);       // 19-bit Mersenne prime
    test_system(65537);        // 17-bit Fermat prime
    
    return 0;
}