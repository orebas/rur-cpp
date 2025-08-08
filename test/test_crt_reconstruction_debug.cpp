#include <iostream>
#include <vector>
#include <gmpxx.h>
#include "julia_rur/multi_modular.hpp"

int main() {
    std::cout << "Testing CRT and rational reconstruction for -2..." << std::endl;
    
    // Test primes
    std::vector<julia_rur::ModularCoeff> primes = {131071, 131063, 131059};
    
    // For each prime p, we have the value (p - 2) which represents -2 mod p
    std::vector<julia_rur::ModularCoeff> remainders;
    for (auto p : primes) {
        remainders.push_back(p - 2);
        std::cout << "Prime " << p << ": " << (p - 2) << " represents -2" << std::endl;
    }
    
    // Apply CRT
    mpz_class result;
    std::vector<mpz_class> crt_multipliers;
    julia_rur::chinese_remainder_theorem(result, remainders, primes, crt_multipliers);
    
    std::cout << "\nCRT result: " << result << std::endl;
    
    // Compute modulus
    mpz_class modulus = 1;
    for (auto p : primes) {
        modulus *= p;
    }
    std::cout << "Modulus: " << modulus << std::endl;
    std::cout << "Modulus/2: " << (modulus/2) << std::endl;
    
    // Check if result is in upper half
    if (result > modulus / 2) {
        std::cout << "\nResult is in upper half of range" << std::endl;
        mpz_class signed_result = result - modulus;
        std::cout << "Signed interpretation: " << signed_result << std::endl;
    }
    
    // Test rational reconstruction directly
    std::cout << "\nTesting rational reconstruction..." << std::endl;
    
    // Test with the unsigned value
    auto bounds = julia_rur::compute_balanced_bounds(modulus);
    auto rat1 = julia_rur::rational_reconstruction(result, modulus, bounds.N, bounds.D);
    std::cout << "Reconstruction of " << result << ": ";
    if (rat1.success) {
        std::cout << rat1.rational << std::endl;
    } else {
        std::cout << "failed" << std::endl;
    }
    
    // Test with signed value
    if (result > modulus / 2) {
        mpz_class signed_val = result - modulus;
        // Note: rational_reconstruction expects non-negative input
        // So we need to convert back
        mpz_class positive_val = signed_val + modulus;
        auto rat2 = julia_rur::rational_reconstruction(positive_val, modulus, bounds.N, bounds.D);
        std::cout << "Reconstruction of signed->positive " << positive_val << ": ";
        if (rat2.success) {
            std::cout << rat2.rational << std::endl;
        } else {
            std::cout << "failed" << std::endl;
        }
    }
    
    // Test the full CRT and reconstruction pipeline
    std::cout << "\nTesting full pipeline..." << std::endl;
    std::vector<std::vector<std::vector<julia_rur::ModularCoeff>>> modular_tables;
    
    // Create a simple table with one polynomial having one coefficient
    for (size_t i = 0; i < primes.size(); ++i) {
        std::vector<std::vector<julia_rur::ModularCoeff>> prime_table;
        std::vector<julia_rur::ModularCoeff> poly_coeffs = {primes[i] - 2}; // -2 mod p
        prime_table.push_back(poly_coeffs);
        modular_tables.push_back(prime_table);
    }
    
    // Reference for verification (not used in this test)
    std::vector<std::vector<julia_rur::ModularCoeff>> reference_mod_p = {{0}};
    julia_rur::ModularCoeff verification_prime = 1073741827;
    
    std::vector<std::vector<mpq_class>> qq_result;
    std::vector<std::vector<mpz_class>> zz_temp;
    mpz_class known_denominator = 1;
    
    auto [success, new_denominator] = julia_rur::crt_and_rational_reconstruction(
        qq_result, zz_temp, modular_tables, primes,
        reference_mod_p, verification_prime, known_denominator
    );
    
    if (success && !qq_result.empty() && !qq_result[0].empty()) {
        std::cout << "Full pipeline result: " << qq_result[0][0] << std::endl;
        std::cout << "Expected: -2" << std::endl;
    } else {
        std::cout << "Full pipeline failed" << std::endl;
    }
    
    return 0;
}