#include "../src/axf4_wrapper.h"
#include <iostream>
#include <vector>
#include <string>

int main() {
    std::cout << "Testing simple multivariate system with specific primes" << std::endl;
    
    // Test system: x^2 + y^2 - 1 = 0, x - y = 0
    // Solutions: x = y = ±1/√2
    // For this to work, we need 2 to be a quadratic residue mod p
    
    std::vector<std::string> variables = {"x", "y"};
    std::vector<const char*> var_ptrs = {variables[0].c_str(), variables[1].c_str()};
    
    // Prime 65537: 2 is a QR since 2^((p-1)/2) = 2^32768 ≡ 1 (mod 65537)
    // Let's verify with a few specific primes
    std::vector<int> primes = {65537, 131071, 262147};
    
    for (int prime : primes) {
        std::cout << "\n=== Testing prime " << prime << " ===" << std::endl;
        
        // Check if 2 is a quadratic residue mod p
        // Using Euler's criterion: a^((p-1)/2) ≡ 1 (mod p) iff a is QR
        long long half = (prime - 1) / 2;
        long long result = 1;
        long long base = 2;
        long long p = prime;
        
        // Fast modular exponentiation
        while (half > 0) {
            if (half & 1) {
                result = (result * base) % p;
            }
            base = (base * base) % p;
            half >>= 1;
        }
        
        std::cout << "2^((p-1)/2) mod p = " << result << std::endl;
        std::cout << "2 is " << (result == 1 ? "a" : "NOT a") << " quadratic residue mod " << prime << std::endl;
        
        // Create F4 session
        axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), 2);
        if (!session) {
            std::cerr << "Failed to create session" << std::endl;
            continue;
        }
        
        // Add polynomials
        if (axf4_add_polynomial(session, "x^2 + y^2 - 1") != 0) {
            std::cerr << "Failed to add first polynomial" << std::endl;
            axf4_destroy_session(session);
            continue;
        }
        
        if (axf4_add_polynomial(session, "x - y") != 0) {
            std::cerr << "Failed to add second polynomial" << std::endl;
            axf4_destroy_session(session);
            continue;
        }
        
        // Compute GB
        axf4_result_t gb_result = axf4_compute_groebner_basis(session);
        
        std::cout << "\nGröbner basis computation:" << std::endl;
        std::cout << "Status: " << gb_result.status << std::endl;
        std::cout << "Basis size: " << gb_result.basis_size << std::endl;
        
        if (gb_result.status == 0) {
            std::cout << "Gröbner basis:" << std::endl;
            std::cout << gb_result.groebner_basis << std::endl;
        }
        
        axf4_free_result(&gb_result);
        axf4_destroy_session(session);
    }
    
    return 0;
}