#include "prime_utils.hpp"
#include <chrono>
#include <algorithm>
#include <set>

namespace julia_rur {

// Simple Miller-Rabin primality test
bool is_probable_prime(ModularCoeff n, int k) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;
    
    // Write n-1 as 2^r * d
    ModularCoeff d = n - 1;
    int r = 0;
    while (d % 2 == 0) {
        d /= 2;
        r++;
    }
    
    // Miller-Rabin test with k rounds
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<ModularCoeff> dis(2, n - 2);
    
    for (int i = 0; i < k; i++) {
        ModularCoeff a = dis(gen);
        ModularCoeff x = 1;
        ModularCoeff base = a;
        ModularCoeff exp = d;
        
        // Compute a^d mod n using binary exponentiation
        while (exp > 0) {
            if (exp % 2 == 1) {
                x = (static_cast<uint64_t>(x) * base) % n;
            }
            base = (static_cast<uint64_t>(base) * base) % n;
            exp /= 2;
        }
        
        if (x == 1 || x == n - 1) continue;
        
        bool composite = true;
        for (int j = 0; j < r - 1; j++) {
            x = (static_cast<uint64_t>(x) * x) % n;
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        
        if (composite) return false;
    }
    
    return true;
}

ModularCoeff generate_random_prime(int min_bits, int max_bits, std::mt19937* rng) {
    // Use provided RNG or create a new one
    std::mt19937 local_gen;
    std::mt19937* gen_ptr = rng;
    if (!gen_ptr) {
        auto seed = std::chrono::steady_clock::now().time_since_epoch().count();
        local_gen = std::mt19937(seed);
        gen_ptr = &local_gen;
    }
    
    // Generate random odd numbers in the range and test for primality
    ModularCoeff min_val = (1u << min_bits);
    ModularCoeff max_val = (1u << max_bits) - 1;
    
    // Generate truly random primes instead of using biased seed primes
    // This avoids systematically choosing unlucky primes
    std::uniform_int_distribution<ModularCoeff> dis(min_val, max_val);
    
    // Try random odd numbers until we find a prime
    for (int attempts = 0; attempts < 10000; attempts++) {
        ModularCoeff candidate = dis(*gen_ptr);
        // Make it odd
        candidate |= 1;
        
        // Ensure it's in bounds
        if (candidate > max_val) candidate = max_val;
        if (candidate < min_val) candidate = min_val | 1;
        
        if (is_probable_prime(candidate, 10)) {
            return candidate;
        }
    }
    
    // Fallback to some known good primes that work well for polynomial systems
    // These are chosen to avoid the problematic 1073741789
    const std::vector<ModularCoeff> fallback_primes = {
        // 30-bit primes that are known to work well
        1073741827, 1073741831, 1073741833, 1073741839, 1073741843,
        1073741857, 1073741873, 1073741909, 1073741939, 1073741941,
        1073741953, 1073741969, 1073741987, 1073742013, 1073742043,
        // 29-bit primes
        536870923, 536870939, 536870951, 536870969, 536871013,
        // 28-bit primes
        268435459, 268435463, 268435501, 268435537, 268435561
    };
    
    // Select randomly from fallback list
    std::uniform_int_distribution<size_t> fallback_dis(0, fallback_primes.size() - 1);
    return fallback_primes[fallback_dis(*gen_ptr)];
}

std::vector<ModularCoeff> get_random_prime_sequence(int num_primes, int min_bits, int max_bits) {
    std::vector<ModularCoeff> result;
    result.reserve(num_primes);
    
    std::mt19937 gen(std::chrono::steady_clock::now().time_since_epoch().count());
    std::set<ModularCoeff> used_primes;
    
    for (int i = 0; i < num_primes; i++) {
        ModularCoeff prime;
        int attempts = 0;
        do {
            prime = generate_random_prime(min_bits, max_bits, &gen);
            attempts++;
        } while (used_primes.count(prime) > 0 && attempts < 100);
        
        used_primes.insert(prime);
        result.push_back(prime);
    }
    
    return result;
}

} // namespace julia_rur