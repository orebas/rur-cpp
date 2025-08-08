#include "../src/julia_rur/multi_modular.hpp"
#include <iostream>
#include <gmpxx.h>

using namespace julia_rur;

int main() {
    // Test rational reconstruction with a negative value
    
    // Example: -2 mod 1073741827 = 1073741825
    mpz_class modulus = 1073741827;
    mpz_class value = 1073741825;  // This represents -2
    
    std::cout << "Testing rational reconstruction:" << std::endl;
    std::cout << "Value: " << value << " (mod " << modulus << ")" << std::endl;
    
    // Convert to symmetric range
    mpz_class symmetric = value;
    if (symmetric > modulus / 2) {
        symmetric -= modulus;
    }
    std::cout << "In symmetric range: " << symmetric << std::endl;
    
    // Test 1: Direct reconstruction of the positive value
    std::cout << "\nTest 1: Reconstruct positive value" << std::endl;
    auto result1 = rational_reconstruction(value, modulus, 0, 0);
    if (result1.success) {
        std::cout << "Result: " << result1.rational << std::endl;
    } else {
        std::cout << "Failed" << std::endl;
    }
    
    // Test 2: Reconstruct after converting back to positive
    std::cout << "\nTest 2: Reconstruct symmetric value converted to positive" << std::endl;
    mpz_class positive = symmetric;
    if (positive < 0) {
        positive += modulus;
    }
    std::cout << "Positive value: " << positive << std::endl;
    auto result2 = rational_reconstruction(positive, modulus, 0, 0);
    if (result2.success) {
        std::cout << "Result: " << result2.rational << std::endl;
    } else {
        std::cout << "Failed" << std::endl;
    }
    
    // Test 3: Check CRT with multiple primes
    std::cout << "\nTest 3: CRT reconstruction" << std::endl;
    std::vector<ModularCoeff> primes = {1073741827, 2147483629};
    std::vector<ModularCoeff> remainders = {1073741825, 2147483627}; // Both represent -2
    
    std::cout << "Primes: ";
    for (auto p : primes) std::cout << p << " ";
    std::cout << std::endl;
    
    std::cout << "Remainders: ";
    for (auto r : remainders) std::cout << r << " ";
    std::cout << std::endl;
    
    // Compute modulus
    mpz_class full_modulus = 1;
    for (auto p : primes) {
        full_modulus *= p;
    }
    
    // Apply CRT manually for 2 primes
    mpz_class crt_result;
    mpz_class p1 = primes[0];
    mpz_class p2 = primes[1];
    mpz_class r1 = remainders[0];
    mpz_class r2 = remainders[1];
    
    // CRT: find x such that x ≡ r1 (mod p1) and x ≡ r2 (mod p2)
    mpz_class m1_inv;
    mpz_invert(m1_inv.get_mpz_t(), p2.get_mpz_t(), p1.get_mpz_t());
    mpz_class m2_inv;
    mpz_invert(m2_inv.get_mpz_t(), p1.get_mpz_t(), p2.get_mpz_t());
    
    crt_result = (r1 * p2 * m1_inv + r2 * p1 * m2_inv) % full_modulus;
    
    std::cout << "CRT result: " << crt_result << std::endl;
    
    // Normalize to symmetric range
    if (crt_result > full_modulus / 2) {
        crt_result -= full_modulus;
    }
    std::cout << "In symmetric range: " << crt_result << std::endl;
    
    // Convert to positive for reconstruction
    mpz_class crt_positive = crt_result;
    if (crt_positive < 0) {
        crt_positive += full_modulus;
    }
    
    auto result3 = rational_reconstruction(crt_positive, full_modulus, 0, 0);
    if (result3.success) {
        std::cout << "Rational reconstruction result: " << result3.rational << std::endl;
    } else {
        std::cout << "Rational reconstruction failed" << std::endl;
    }
    
    return 0;
}