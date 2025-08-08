#include <iostream>
#include <vector>
#include <cstdint>

using ModularCoeff = uint32_t;

// Simple implementation for testing
std::vector<ModularCoeff> generate_primes(size_t bits, size_t count) {
    std::vector<ModularCoeff> primes;
    
    if (bits <= 28) {
        // Primes around 2^28 - matching Julia's default pr_max_bitsize = 28
        const ModularCoeff primes_28bit[] = {
            268435459, 268435463, 268435493, 268435537, 268435561,
            268435577, 268435579, 268435597, 268435607, 268435631
        };
        
        for (size_t i = 0; i < count && i < 10; ++i) {
            primes.push_back(primes_28bit[i]);
        }
    }
    
    return primes;
}

int main() {
    std::cout << "Testing prime generation for 28-bit primes (Julia's default):\n" << std::endl;
    
    auto primes = generate_primes(28, 10);
    
    for (auto p : primes) {
        std::cout << "Prime: " << p << std::endl;
        
        // Check if it's safe for F4
        if (p <= 2147483647) {
            std::cout << "  ✓ Safe for F4 (below 2^31-1)" << std::endl;
        } else {
            std::cout << "  ✗ UNSAFE for F4 (exceeds 2^31-1)" << std::endl;
        }
        
        // Show bit size
        int bits = 32 - __builtin_clz(p);
        std::cout << "  Bit size: " << bits << std::endl;
    }
    
    std::cout << "\nConclusion: All 28-bit primes are well below the F4 bug threshold of 2^31-1 = 2147483647" << std::endl;
    
    return 0;
}