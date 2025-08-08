#include <iostream>
#include <vector>

using ModularCoeff = unsigned int;

// Simulated prime generation
std::vector<ModularCoeff> generate_primes(size_t bits, size_t count) {
    if (bits == 28) {
        // Return our good primes where 2 is QR
        const ModularCoeff primes[] = {
            268435399, 268435367, 268435361, 268435337, 268435313
        };
        std::vector<ModularCoeff> result;
        for (size_t i = 0; i < count && i < 5; ++i) {
            result.push_back(primes[i]);
        }
        return result;
    }
    return {};
}

// Simulated prev_primes - goes backwards from start
std::vector<ModularCoeff> prev_primes(ModularCoeff start, size_t count) {
    std::cout << "prev_primes called with start=" << start << ", count=" << count << std::endl;
    
    // In the real code, this would generate primes going backwards
    // This shows the issue - it doesn't use our curated list!
    std::vector<ModularCoeff> result;
    ModularCoeff p = start;
    
    // Simulate finding primes backwards
    if (p == 268435312) {  // This would be 268435313 - 1
        // Next primes going backwards would include bad ones
        const ModularCoeff bad_primes[] = {
            268435291, // mod 8 = 3, BAD
            268435243, // mod 8 = 3, BAD  
            268435171, // mod 8 = 3, BAD
        };
        for (size_t i = 0; i < count && i < 3; ++i) {
            result.push_back(bad_primes[i]);
        }
    }
    
    return result;
}

int main() {
    std::cout << "Simulating prime generation flow:\n" << std::endl;
    
    // Initial batch - uses generate_primes
    size_t batch_size = 5;
    auto initial_primes = generate_primes(28, batch_size);
    
    std::cout << "Initial batch (from generate_primes):" << std::endl;
    for (auto p : initial_primes) {
        std::cout << "  " << p << " (mod 8 = " << (p % 8) << ")" << std::endl;
    }
    
    // Subsequent batch - uses prev_primes
    ModularCoeff current_prime = initial_primes.back() - 1;
    std::cout << "\nNext batch starting from " << current_prime << ":" << std::endl;
    
    auto next_primes = prev_primes(current_prime, batch_size);
    for (auto p : next_primes) {
        std::cout << "  " << p << " (mod 8 = " << (p % 8) << ") - BAD for sqrt(2)!" << std::endl;
    }
    
    std::cout << "\nThis shows the issue: prev_primes doesn't use our curated list!" << std::endl;
    
    return 0;
}