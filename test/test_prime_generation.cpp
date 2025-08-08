#include "../src/julia_rur/multi_modular.hpp"
#include <iostream>

using namespace julia_rur;

int main() {
    std::cout << "Testing prime generation at different bit sizes:\n" << std::endl;
    
    // Test 17-bit primes
    std::cout << "17-bit primes:" << std::endl;
    auto primes17 = generate_primes(17, 5);
    for (auto p : primes17) {
        std::cout << "  " << p << " (bits: " << (32 - __builtin_clz(p)) << ")" << std::endl;
    }
    
    // Test 20-bit primes
    std::cout << "\n20-bit primes:" << std::endl;
    auto primes20 = generate_primes(20, 5);
    for (auto p : primes20) {
        std::cout << "  " << p << " (bits: " << (32 - __builtin_clz(p)) << ")" << std::endl;
    }
    
    // Test 28-bit primes (Julia's default)
    std::cout << "\n28-bit primes (Julia's default):" << std::endl;
    auto primes28 = generate_primes(28, 10);
    for (auto p : primes28) {
        std::cout << "  " << p << " (bits: " << (32 - __builtin_clz(p)) << ")" << std::endl;
    }
    
    // Test 30-bit primes
    std::cout << "\n30-bit primes:" << std::endl;
    auto primes30 = generate_primes(30, 5);
    for (auto p : primes30) {
        std::cout << "  " << p << " (bits: " << (32 - __builtin_clz(p)) << ")" << std::endl;
    }
    
    // Test 31-bit primes
    std::cout << "\n31-bit primes (F4 bug threshold):" << std::endl;
    auto primes31 = generate_primes(31, 5);
    for (auto p : primes31) {
        std::cout << "  " << p << " (bits: " << (32 - __builtin_clz(p)) << ")" << std::endl;
        if (p > 2147483647) {
            std::cout << "    WARNING: This prime exceeds F4's safe limit!" << std::endl;
        }
    }
    
    // Check all 28-bit primes are below F4 threshold
    std::cout << "\nVerifying all 28-bit primes are safe for F4:" << std::endl;
    bool all_safe = true;
    for (auto p : primes28) {
        if (p > 2147483647) {
            std::cout << "  ERROR: Prime " << p << " exceeds F4 limit!" << std::endl;
            all_safe = false;
        }
    }
    if (all_safe) {
        std::cout << "  ✓ All 28-bit primes are safe for F4" << std::endl;
    }
    
    return 0;
}