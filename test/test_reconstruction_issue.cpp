#include <iostream>
#include <vector>
#include <gmpxx.h>
#include "julia_rur/multi_modular.hpp"

int main() {
    std::cout << "Testing specific reconstruction issue..." << std::endl;
    
    // The actual values from our computation
    std::vector<julia_rur::ModularCoeff> primes = {
        131063, 131059, 131041, 131023, 131011
    };
    
    // For -2, we have p-2 for each prime
    std::vector<julia_rur::ModularCoeff> remainders;
    for (auto p : primes) {
        julia_rur::ModularCoeff val = p - 2;
        remainders.push_back(val);
        std::cout << "Prime " << p << ": remainder = " << val << std::endl;
    }
    
    // Compute modulus
    mpz_class modulus = 1;
    for (auto p : primes) {
        modulus *= p;
    }
    std::cout << "\nModulus = " << modulus << std::endl;
    
    // Apply CRT
    mpz_class crt_result;
    std::vector<mpz_class> crt_multipliers;
    julia_rur::chinese_remainder_theorem(crt_result, remainders, primes, crt_multipliers);
    
    std::cout << "CRT result = " << crt_result << std::endl;
    std::cout << "Modulus/2 = " << (modulus/2) << std::endl;
    
    // Check symmetric range
    mpz_class signed_value = crt_result;
    if (crt_result > modulus / 2) {
        signed_value = crt_result - modulus;
        std::cout << "In upper half, signed value = " << signed_value << std::endl;
    }
    
    // Now test rational reconstruction with known denominator = 1
    mpz_class known_denom = 1;
    
    // Test different reconstruction strategies
    auto bounds_9to1 = julia_rur::compute_unbalanced_bounds_9to1(modulus);
    auto bounds_const = julia_rur::compute_constant_bounds(modulus);
    auto bounds_balanced = julia_rur::compute_balanced_bounds(modulus);
    
    std::cout << "\nTesting reconstruction strategies:" << std::endl;
    
    // Strategy 1: 9:1 with known denominator
    auto rat1 = julia_rur::rational_reconstruction_with_denominator_no_verify(
        signed_value, known_denom, modulus, bounds_9to1.N, bounds_9to1.D
    );
    std::cout << "Strategy 1 (9:1): ";
    if (rat1.success) {
        std::cout << rat1.rational << std::endl;
    } else {
        std::cout << "failed" << std::endl;
    }
    
    // Try with unsigned value
    auto rat1b = julia_rur::rational_reconstruction_with_denominator_no_verify(
        crt_result, known_denom, modulus, bounds_9to1.N, bounds_9to1.D
    );
    std::cout << "Strategy 1b (9:1 unsigned): ";
    if (rat1b.success) {
        std::cout << rat1b.rational << std::endl;
    } else {
        std::cout << "failed" << std::endl;
    }
    
    // Strategy 2: constant bounds
    auto rat2 = julia_rur::rational_reconstruction_with_denominator_no_verify(
        signed_value, known_denom, modulus, bounds_const.N, bounds_const.D
    );
    std::cout << "Strategy 2 (constant): ";
    if (rat2.success) {
        std::cout << rat2.rational << std::endl;
    } else {
        std::cout << "failed" << std::endl;
    }
    
    // Direct reconstruction test
    std::cout << "\nDirect reconstruction of CRT result:" << std::endl;
    auto direct = julia_rur::rational_reconstruction(
        crt_result, modulus, bounds_balanced.N, bounds_balanced.D
    );
    if (direct.success) {
        std::cout << "Direct: " << direct.rational << std::endl;
    } else {
        std::cout << "Direct: failed" << std::endl;
    }
    
    return 0;
}