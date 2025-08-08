#include <iostream>

bool is_problematic_prime_for_f4(unsigned int prime) {
    // Check if prime is a Mersenne prime (2^k - 1)
    unsigned int val = prime + 1;
    int k = 0;
    while ((val & 1) == 0) {
        val >>= 1;
        k++;
    }
    if (val == 1 && k > 0) {
        std::cout << "  Mersenne: 2^" << k << " - 1" << std::endl;
        return true;
    }
    
    // Check if prime is 2^k - c for small c (c < 100)
    for (int bits = 2; bits <= 31; bits++) {
        unsigned long long two_k = 1ULL << bits;
        if (two_k > prime && two_k - prime < 100) {
            std::cout << "  Form: 2^" << bits << " - " << (two_k - prime) << std::endl;
            return true;
        }
    }
    
    return false;
}

int main() {
    unsigned int test_primes[] = {
        1048447, 1048433, 1048423, 1048417, 1048391,
        524287, 1048573, 268435399
    };
    
    for (auto p : test_primes) {
        std::cout << "\nPrime " << p << ":";
        if (is_problematic_prime_for_f4(p)) {
            std::cout << " PROBLEMATIC";
        } else {
            std::cout << " OK";
        }
        std::cout << std::endl;
    }
    
    return 0;
}