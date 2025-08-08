#include <iostream>
#include <vector>

using ModularCoeff = unsigned int;

// Test how F4 handles negative coefficients
void test_negative_representation(ModularCoeff prime) {
    std::cout << "\nPrime " << prime << ":" << std::endl;
    
    // Test -1 mod p
    ModularCoeff neg_one = prime - 1;
    std::cout << "  -1 mod p = " << neg_one << std::endl;
    
    // Test if it's stored differently
    if (neg_one > prime / 2) {
        std::cout << "  Interpreted as negative: -" << (prime - neg_one) << std::endl;
    }
    
    // Test 0 - 1 mod p
    ModularCoeff zero = 0;
    ModularCoeff one = 1;
    ModularCoeff result;
    
    // Standard modular subtraction
    if (zero >= one) {
        result = zero - one;
    } else {
        result = zero + prime - one;
    }
    std::cout << "  0 - 1 mod p = " << result << std::endl;
    
    // Check if result equals neg_one
    if (result == neg_one) {
        std::cout << "  ✓ Correct: 0 - 1 = -1 mod p" << std::endl;
    } else {
        std::cout << "  ✗ ERROR: 0 - 1 ≠ -1 mod p!" << std::endl;
    }
}

int main() {
    std::cout << "Testing modular arithmetic representation:" << std::endl;
    
    // Working prime
    test_negative_representation(1073741827);
    
    // Failing primes
    test_negative_representation(1048573);
    test_negative_representation(524287);
    test_negative_representation(268435399);
    
    return 0;
}