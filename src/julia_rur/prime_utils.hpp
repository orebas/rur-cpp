#ifndef JULIA_RUR_PRIME_UTILS_HPP
#define JULIA_RUR_PRIME_UTILS_HPP

#include "data_structures.hpp"
#include <random>
#include <vector>

namespace julia_rur {

/**
 * @brief Generate a random prime in a specified bit range
 * 
 * @param min_bits Minimum number of bits (e.g., 28)
 * @param max_bits Maximum number of bits (e.g., 30)
 * @param rng Random number generator (optional)
 * @return A random prime in the specified range
 */
ModularCoeff generate_random_prime(int min_bits = 28, int max_bits = 30, 
                                   std::mt19937* rng = nullptr);

/**
 * @brief Get a list of known good primes for multi-modular computation
 * 
 * These are primes that are known to work well with F4 and are
 * unlikely to be systematically unlucky for common polynomial systems.
 * 
 * @param num_primes Number of primes to return
 * @param min_bits Minimum bit size
 * @param max_bits Maximum bit size  
 * @return Vector of distinct primes
 */
std::vector<ModularCoeff> get_random_prime_sequence(int num_primes, 
                                                    int min_bits = 28, 
                                                    int max_bits = 30);

/**
 * @brief Check if a number is likely prime using Miller-Rabin test
 * 
 * @param n Number to test
 * @param k Number of rounds (default 5 gives error probability < 1/1000)
 * @return true if probably prime, false if definitely composite
 */
bool is_probable_prime(ModularCoeff n, int k = 5);

} // namespace julia_rur

#endif // JULIA_RUR_PRIME_UTILS_HPP