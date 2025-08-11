#include <iostream>
#include <vector>
#include <string>
#include "julia_rur/prime_utils.hpp"

int main() {
    std::cout << "Testing prime generation for Cyclic-5...\n\n";
    
    // Generate a bunch of random primes to see what we're getting
    for (int i = 0; i < 10; i++) {
        auto prime = julia_rur::generate_random_prime(28, 30);
        std::cout << "Random prime " << i+1 << ": " << prime << "\n";
    }
    
    std::cout << "\nTesting prime sequence generation:\n";
    auto primes = julia_rur::get_random_prime_sequence(5, 28, 30);
    for (size_t i = 0; i < primes.size(); i++) {
        std::cout << "Prime " << i+1 << ": " << primes[i] << "\n";
    }
    
    // Check the default prime that was previously hardcoded
    std::cout << "\nPrevious default prime: 1073741827\n";
    std::cout << "Is it prime? " << (julia_rur::is_probable_prime(1073741827) ? "yes" : "no") << "\n";
    
    return 0;
}