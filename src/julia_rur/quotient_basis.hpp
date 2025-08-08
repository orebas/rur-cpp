#ifndef JULIA_RUR_QUOTIENT_BASIS_HPP
#define JULIA_RUR_QUOTIENT_BASIS_HPP

#include "data_structures.hpp"
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace julia_rur {

/**
 * @brief Extract quotient basis from Gröbner basis leading terms
 * 
 * Computes the standard monomials (those not divisible by any leading monomial).
 * These form a basis for the quotient ring k[x1,...,xn]/I.
 * 
 * Algorithm: Starting from monomial 1, explore all monomials by increasing
 * each variable degree, keeping only those not divisible by any leading term.
 * 
 * @param leading_terms Leading terms from Gröbner basis (degrevlex order)
 * @return Vector of standard monomials forming quotient basis (sorted degrevlex)
 * @throws std::domain_error if ideal is not zero-dimensional
 * @throws std::runtime_error if system has no solutions (GB = {1})
 */
std::vector<PP> compute_quotient_basis(const std::vector<PP>& leading_terms);

/**
 * @brief Find if a monomial is divisible by any leading term
 * 
 * @param monomial The monomial to test
 * @param leading_terms Leading terms from Gröbner basis
 * @return Index (1-based) of first divisor, or 0 if not divisible
 */
int find_divisor(const PP& monomial, const std::vector<PP>& leading_terms);


/**
 * @brief Check if ideal is zero-dimensional
 * 
 * An ideal is zero-dimensional if for each variable xi, there exists
 * a leading term that is a pure power of xi (has degree > 0 in xi only).
 * 
 * @param leading_terms Leading terms from Gröbner basis
 * @param num_variables Number of variables
 * @return true if zero-dimensional, false otherwise
 */
bool is_zero_dimensional(const std::vector<PP>& leading_terms, int num_variables);

} // namespace julia_rur

#endif // JULIA_RUR_QUOTIENT_BASIS_HPP