#pragma once

#include "bivariate_algorithm.hpp" // For modular_inverse function
#include "data_structures.hpp"
#include <gmpxx.h> // For big integer arithmetic
#include <numeric>
#include <vector>


/**
 * Multi-modular computation framework
 * Mirrors RationalUnivariateRepresentation.jl multi-modular algorithms
 *
 * This module provides:
 * - Chinese Remainder Theorem (CRT) for combining modular results
 * - Rational reconstruction to recover Q coefficients from Z_p
 * - Multi-threaded modular computation management
 * - Prime generation and management
 */

namespace julia_rur {

/**
 * Result of rational reconstruction
 */
struct RationalReconstructionResult {
    bool success;
    mpq_class rational; // GMP rational number

    RationalReconstructionResult()
      : success(false)
      , rational(0) {}
    RationalReconstructionResult(bool s, const mpq_class &r)
      : success(s)
      , rational(r) {}
};

/**
 * Multi-modular computation context
 */
struct MultiModularContext {
    std::vector<ModularCoeff> primes;       // List of primes used
    std::vector<mpz_class> crt_multipliers; // Precomputed CRT multipliers
    mpz_class modulus;                      // Product of all primes

    MultiModularContext() = default;
};

/**
 * Generate the next prime less than the given value
 * Mirrors Julia's PrevPrime function
 */
ModularCoeff
prev_prime(ModularCoeff n);

/**
 * Generate a list of n primes less than start
 * Mirrors Julia's PrevPrimes function
 */
std::vector<ModularCoeff>
prev_primes(ModularCoeff start, size_t count);

/**
 * Check if a number is prime
 * Helper function for prime generation
 */
bool
is_prime(ModularCoeff n);

/**
 * Chinese Remainder Theorem - combine modular results
 * Mirrors Groebner.crt! function
 *
 * @param result Output: combined result as big integer
 * @param remainders Remainders modulo each prime
 * @param primes List of primes used
 * @param multipliers Precomputed CRT multipliers (can be empty on first call)
 * @return Updated modulus (product of all primes)
 */
mpz_class
chinese_remainder_theorem(mpz_class &result,
                          const std::vector<ModularCoeff> &remainders,
                          const std::vector<ModularCoeff> &primes,
                          std::vector<mpz_class> &multipliers);

/**
 * Rational reconstruction using extended GCD
 * Mirrors Groebner.ratrec_nemo function
 *
 * Given an integer a and modulus m, find rational p/q such that:
 * - a ≡ p/q (mod m)
 * - |p| ≤ N and 0 < q ≤ D
 *
 * @param a Integer to reconstruct
 * @param m Modulus
 * @param N Numerator bound (default: sqrt(m/2))
 * @param D Denominator bound (default: sqrt(m/2))
 * @return RationalReconstructionResult with success flag and rational
 */
RationalReconstructionResult
rational_reconstruction(const mpz_class &a,
                        const mpz_class &m,
                        const mpz_class &N = mpz_class(0),
                        const mpz_class &D = mpz_class(0));

/**
 * Try rational reconstruction with known denominator
 * Mirrors ratrec_try! function
 *
 * @param zz Integer value from CRT
 * @param den Known denominator to try
 * @param modulus Current modulus
 * @param N Numerator bound
 * @param D Denominator bound
 * @param zp Original value mod p (for verification)
 * @param p Prime used for verification
 * @return RationalReconstructionResult
 */
RationalReconstructionResult
rational_reconstruction_with_denominator(const mpz_class &zz,
                                         const mpz_class &den,
                                         const mpz_class &modulus,
                                         const mpz_class &N,
                                         const mpz_class &D,
                                         ModularCoeff zp,
                                         ModularCoeff p);

/**
 * Compute bounds for rational reconstruction
 * Different strategies as in Julia code
 */
struct RationalReconstructionBounds {
    mpz_class N; // Numerator bound
    mpz_class D; // Denominator bound
};

RationalReconstructionBounds
compute_balanced_bounds(const mpz_class &modulus);
RationalReconstructionBounds
compute_unbalanced_bounds_9to1(const mpz_class &modulus);
RationalReconstructionBounds
compute_unbalanced_bounds_99to1(const mpz_class &modulus);
RationalReconstructionBounds
compute_constant_bounds(const mpz_class &modulus);

/**
 * Combine CRT and rational reconstruction for a matrix of values
 * Mirrors crt_and_ratrec! function
 *
 * @param qq_result Output: rational matrix
 * @param zz_temp Working space: integer matrix for CRT results
 * @param modular_tables Input: tables from different primes
 * @param primes List of primes used
 * @param reference_mod_p Reference values for verification
 * @param verification_prime Prime for verification
 * @param known_denominator Known common denominator (starts at 1)
 * @return Pair of (success, lcm_of_denominators)
 */
std::pair<bool, mpz_class>
crt_and_rational_reconstruction(std::vector<std::vector<mpq_class>> &qq_result,
                                std::vector<std::vector<mpz_class>> &zz_temp,
                                const std::vector<std::vector<std::vector<ModularCoeff>>> &modular_tables,
                                const std::vector<ModularCoeff> &primes,
                                const std::vector<std::vector<ModularCoeff>> &reference_mod_p,
                                ModularCoeff verification_prime,
                                mpz_class &known_denominator);

/**
 * Align a number to the nearest multiple of n (round up)
 * Mirrors align_to function
 */
/**
 * @brief Align a number to the nearest multiple of n (round up)
 *
 * @param x Number to align
 * @param n Alignment size
 * @return Aligned number (next multiple of n)
 */
int32_t inline align_to(int32_t x, int32_t n) { return ((x + n - 1) / n) * n; }

/**
 * Modular coefficient vector conversion
 * Convert rational coefficients to modular arithmetic
 */
std::vector<ModularCoeff>
modular_coeffs_vect(const std::vector<mpq_class> &rational_coeffs, ModularCoeff prime);

/**
 * Extended GCD for rational reconstruction
 * Helper function that computes gcd(a,b) and Bezout coefficients
 */
struct ExtendedGCDResult {
    mpz_class gcd;
    mpz_class x; // Coefficient for a
    mpz_class y; // Coefficient for b
};

ExtendedGCDResult
extended_gcd(const mpz_class &a, const mpz_class &b);

/**
 * Rational reconstruction with denominator but without verification
 * Internal use for debugging
 */
RationalReconstructionResult
rational_reconstruction_with_denominator_no_verify(const mpz_class &zz,
                                                   const mpz_class &den,
                                                   const mpz_class &modulus,
                                                   const mpz_class &N,
                                                   const mpz_class &D);

/**
 * @brief Generate primes for multi-modular computation
 *
 * Generates n primes of approximately the specified bit size
 * suitable for modular arithmetic.
 *
 * @param bits Approximate bit size of primes
 * @param count Number of primes to generate
 * @return Vector of primes
 */
std::vector<ModularCoeff>
generate_primes(size_t bits, size_t count);

void
crt_precompute(const mpz_class &modulus,
               mpz_class &n1,
               mpz_class &n2,
               std::vector<mpz_class> &mults,
               const std::vector<ModularCoeff> &primes);


} // namespace julia_rur