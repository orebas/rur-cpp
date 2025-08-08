#include <iostream>
#include <cmath>

void analyze_prime(unsigned int p) {
    std::cout << "\nPrime " << p << ":" << std::endl;
    
    // Check if it's a Mersenne prime (2^k - 1)
    unsigned int val = p + 1;
    int k = 0;
    while ((val & 1) == 0) {
        val >>= 1;
        k++;
    }
    if (val == 1 && k > 0) {
        std::cout << "  Mersenne prime: 2^" << k << " - 1" << std::endl;
    }
    
    // Check if it's 2^k - c for small c
    for (int k = 2; k <= 32; k++) {
        unsigned long long two_k = 1ULL << k;
        if (two_k > p && two_k - p < 100) {
            std::cout << "  Form: 2^" << k << " - " << (two_k - p) << std::endl;
        }
    }
    
    // Check mod 8 (for quadratic residue properties)
    std::cout << "  mod 8 = " << (p % 8) << std::endl;
    
    // Check bit size
    int bits = 0;
    unsigned int temp = p;
    while (temp > 0) {
        bits++;
        temp >>= 1;
    }
    std::cout << "  bit size: " << bits << std::endl;
}

int main() {
    std::cout << "Analysis of primes used in F4 testing:" << std::endl;
    
    // Working prime
    std::cout << "\n=== WORKING PRIME ===" << std::endl;
    analyze_prime(1073741827);
    
    // Failing primes
    std::cout << "\n=== FAILING PRIMES ===" << std::endl;
    analyze_prime(1048573);
    analyze_prime(524287);
    analyze_prime(268435399);
    
    return 0;
}