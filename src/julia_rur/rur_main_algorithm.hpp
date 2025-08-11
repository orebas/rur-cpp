#ifndef JULIA_RUR_MAIN_ALGORITHM_HPP
#define JULIA_RUR_MAIN_ALGORITHM_HPP

#include "../axf4_wrapper.h"
#include "data_structures.hpp"
#include "f4_integration.hpp"
#include "multi_modular.hpp"
#include "multiplication_tables.hpp"
#include "quotient_basis.hpp"
#include "univariate_parameterization.hpp"
#include <memory>
#include <vector>


namespace julia_rur {

/**
 * @brief Result of a single modular RUR computation
 *
 * Contains all data from one prime computation for CRT reconstruction
 */
struct ModularRURResult {
    ModularCoeff prime;
    MinimalPolynomialResult minimal_polynomial;
    std::vector<BivariateResult> parameterizations;
    std::vector<PP> quotient_basis;

    // Variable position mapping: var_positions[i] = position of variable i+1 in quotient basis
    // -1 if variable is not in basis (shouldn't happen for proper systems)
    std::vector<int32_t> var_positions;

    bool success;
};

/**
 * @brief Final RUR result with rational coefficients
 *
 * Contains the complete rational univariate representation
 */
struct RationalRURResult {
    std::vector<mpq_class> minimal_polynomial;      // Rational coefficients of f(T)
    std::vector<std::vector<mpq_class>> numerators; // Numerators of xi = gi(T)/f'(T)
    mpq_class denominator_derivative;               // f'(T) as common denominator
    std::vector<PP> quotient_basis;                 // For reference
    bool success;
    std::string error_message;
};

/**
 * @brief Configuration for RUR computation
 */
struct RURConfig {
    size_t initial_prime_bits = 31; // Starting prime size
    size_t max_prime_bits = 63;     // Maximum prime size
    size_t num_threads = 1;         // Parallelization (future)
    SeparatingStrategy separating_strategy = SeparatingStrategy::CURRENT;
    bool verbose = false;
    bool timing = false; // enable coarse timing logs
    bool track_multiplicities = true; // Track and report multiplicities for non-radical ideals
    bool apply_square_free = true;    // Apply square-free reduction (required for non-radical ideals)
};

/**
 * @brief Compute RUR for a single prime
 *
 * This is the core single-prime RUR algorithm that:
 * 1. Computes Gröbner basis using F4
 * 2. Extracts quotient basis
 * 3. Builds multiplication tables
 * 4. Finds separating element
 * 5. Computes univariate parameterization
 *
 * @param polynomials Input polynomial system
 * @param variables Variable names
 * @param prime Prime for modular arithmetic
 * @param config Configuration options
 * @return ModularRURResult
 */
std::pair<ModularRURResult, std::vector<int>>
compute_modular_rur(const std::vector<std::string> &polynomials,
                    const std::vector<std::string> &variables,
                    ModularCoeff prime,
                    const RURConfig &config,
                    const std::vector<int> &separating_element = std::vector<int>());

/**
 * @brief Multi-modular RUR computation with rational reconstruction
 *
 * Main algorithm that:
 * 1. Computes RUR modulo multiple primes
 * 2. Uses Chinese Remainder Theorem for lifting
 * 3. Performs rational reconstruction
 * 4. Returns exact rational result
 *
 * @param polynomials Input polynomial system
 * @param variables Variable names
 * @param config Configuration options
 * @return RationalRURResult
 */
RationalRURResult
compute_rational_rur(const std::vector<std::string> &polynomials,
                     const std::vector<std::string> &variables,
                     const RURConfig &config = RURConfig());

/**
 * @brief Check if system is zero-dimensional
 *
 * Quick check before full RUR computation
 *
 * @param polynomials Input polynomial system
 * @param variables Variable names
 * @param prime Prime for test computation
 * @return true if zero-dimensional, false otherwise
 */
bool
is_zero_dimensional_system(const std::vector<std::string> &polynomials,
                           const std::vector<std::string> &variables,
                           ModularCoeff prime = 100003);

/**
 * @brief Compute only the quotient ring dimension
 *
 * Useful for checking system properties without full RUR
 *
 * @param polynomials Input polynomial system
 * @param variables Variable names
 * @param prime Prime for computation
 * @return Dimension of quotient ring, or -1 on error
 */
int
compute_quotient_dimension(const std::vector<std::string> &polynomials,
                           const std::vector<std::string> &variables,
                           ModularCoeff prime = 100003);

/**
 * @brief Pretty-print RUR result
 *
 * Formats the rational RUR result for display
 *
 * @param result The RUR result to print
 * @param variables Variable names for display
 * @return Formatted string representation
 */
std::string
format_rur_result(const RationalRURResult &result, const std::vector<std::string> &variables);

} // namespace julia_rur

#endif // JULIA_RUR_MAIN_ALGORITHM_HPP