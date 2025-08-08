#include "../src/julia_rur/multi_modular.hpp"
#include <iostream>
#include <cassert>

using namespace julia_rur;

void test_simple_reconstruction() {
    std::cout << "Testing simple rational reconstruction" << std::endl;
    
    // Test reconstructing -2/1
    ModularCoeff p1 = 131063;  // First prime
    ModularCoeff p2 = 131059;  // Second prime
    
    // -2 mod p1 = p1 - 2
    ModularCoeff val_p1 = p1 - 2;
    // -2 mod p2 = p2 - 2  
    ModularCoeff val_p2 = p2 - 2;
    
    std::cout << "-2 mod " << p1 << " = " << val_p1 << std::endl;
    std::cout << "-2 mod " << p2 << " = " << val_p2 << std::endl;
    
    // Apply CRT
    std::vector<ModularCoeff> remainders = {val_p1, val_p2};
    std::vector<ModularCoeff> primes = {p1, p2};
    std::vector<mpz_class> crt_multipliers;
    
    mpz_class result;
    chinese_remainder_theorem(result, remainders, primes, crt_multipliers);
    
    mpz_class modulus = static_cast<mpz_class>(p1) * p2;
    std::cout << "\nCRT result: " << result << std::endl;
    std::cout << "Modulus: " << modulus << std::endl;
    
    // Compute bounds
    auto bounds = compute_balanced_bounds(modulus);
    std::cout << "\nBalanced bounds: N=" << bounds.N << ", D=" << bounds.D << std::endl;
    
    // Try standard rational reconstruction
    auto rat_result = rational_reconstruction(result, modulus, bounds.N, bounds.D);
    
    if (rat_result.success) {
        std::cout << "Reconstruction succeeded: " << rat_result.rational << std::endl;
        assert(rat_result.rational == mpq_class(-2, 1));
    } else {
        std::cout << "Reconstruction FAILED!" << std::endl;
        
        // Try with larger bounds
        mpz_class large_N = modulus / 4;
        mpz_class large_D = 2;
        std::cout << "\nTrying with larger bounds: N=" << large_N << ", D=" << large_D << std::endl;
        
        auto rat2 = rational_reconstruction(result, modulus, large_N, large_D);
        if (rat2.success) {
            std::cout << "Reconstruction succeeded: " << rat2.rational << std::endl;
        } else {
            std::cout << "Still failed!" << std::endl;
        }
    }
}

int main() {
    test_simple_reconstruction();
    return 0;
}