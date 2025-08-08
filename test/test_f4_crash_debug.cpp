#include <iostream>
#include "../src/axf4_wrapper.h"
#include <vector>
#include <cstdint>

void test_prime(uint32_t prime, const std::string& poly_str) {
    std::cout << "Testing prime " << prime << " with polynomial: " << poly_str << std::endl;
    
    const char* vars[] = {"x"};
    axf4_session_t session = axf4_create_session(prime, vars, 1);
    
    if (!session) {
        std::cerr << "  Failed to create session" << std::endl;
        return;
    }
    
    if (axf4_add_polynomial(session, poly_str.c_str()) != 0) {
        std::cerr << "  Failed to add polynomial" << std::endl;
        axf4_destroy_session(session);
        return;
    }
    
    std::cout << "  Computing GB..." << std::endl;
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    if (result.status == 0) {
        std::cout << "  Success! GB size: " << result.basis_size << std::endl;
    } else {
        std::cout << "  Failed with status: " << result.status << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
}

int main() {
    std::cout << "Testing F4 with different primes and polynomials\n\n";
    
    // Test different polynomial systems
    std::vector<std::string> polys = {
        "x - 3",
        "x^2 - 2",
        "x^3 - 2*x - 5",
        "x^2 + x + 1"
    };
    
    // Test different prime sizes
    std::vector<uint32_t> primes = {
        65537,        // 2^16 + 1
        131063,       // Small prime that works
        1048583,      // ~2^20
        1073741827,   // The problematic prime
        2147483629    // Large 31-bit prime
    };
    
    for (const auto& poly : polys) {
        std::cout << "\n=== Testing polynomial: " << poly << " ===\n";
        for (uint32_t prime : primes) {
            test_prime(prime, poly);
        }
    }
    
    // Try to isolate the crash
    std::cout << "\n=== Specific crash test ===\n";
    std::cout << "Testing the known problematic case: x - 3 with prime 1073741827\n";
    test_prime(1073741827, "x - 3");
    
    std::cout << "\nTesting with x^2 - 2 and prime 1073741827\n";
    test_prime(1073741827, "x^2 - 2");
    
    return 0;
}