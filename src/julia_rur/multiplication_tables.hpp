#pragma once

#include "data_structures.hpp"
#include <algorithm>
#include <tuple>

/**
 * Julia-style multiplication table construction
 * Mirrors RationalUnivariateRepresentation.jl multiplication table algorithms exactly
 */

namespace julia_rur {

/**
 * Result structure for find_in_border function
 */
struct BorderSearchResult {
    int32_t flag;        // 0=not found, 1=found exact match, 2=found predecessor
    int32_t var_index;   // Variable index (directly from Julia's dd value)
    int32_t pos;         // Position in border structure
    
    BorderSearchResult(int32_t f = 0, int32_t v_idx = 0, int32_t p = 0) 
        : flag(f), var_index(v_idx), pos(p) {}
};

/**
 * Result structure for Julia's divides function that returns (bool, var_index)
 */
struct DividesResult {
    bool divides;
    int32_t var_index;
    
    DividesResult(bool d = false, int32_t v = 0) : divides(d), var_index(v) {}
};

/**
 * Julia-style divides function: returns (bool, var_index)
 * Returns true if monomial a divides monomial b, and the variable index where division fails
 */
DividesResult divides_with_var(const PP& a, const PP& b);

/**
 * Multiply power product by variable: m * x_i
 * Mirrors Julia's mul_pp_by_var! function
 */
PP mul_pp_by_var(const PP& m, int32_t var_index);

/**
 * Find element in border structure
 * Mirrors Julia's find_in_border function exactly
 * 
 * @param m Target monomial to find
 * @param t_xw Border structure to search
 * @return BorderSearchResult with flag, variable/predecessor, and position
 */
BorderSearchResult find_in_border(const PP& m, const std::vector<StackVect>& t_xw);

/**
 * Prepare multiplication table indices
 * Mirrors Julia's prepare_table_mxi function exactly
 * 
 * @param ltg Leading terms of Gröbner basis
 * @param kb Quotient basis (quo)
 * @param t_xw Output: border structure 
 * @param i_xw Output: variable multiplication indices
 */
void prepare_table_mxi(
    const std::vector<PP>& ltg,                           // GB leading terms
    const std::vector<PP>& kb,                            // Quotient basis
    std::vector<StackVect>& t_xw,                         // Border structure (output)
    std::vector<std::vector<int32_t>>& i_xw               // Variable indices (output)
);

/**
 * Learn and compute multiplication table coefficients
 * Mirrors Julia's learn_compute_table! function
 * 
 * @param t_v Output: coefficient vectors
 * @param t_xw Border structure
 * @param i_xw Variable indices
 * @param prime Modular arithmetic prime
 */
void learn_compute_table(
    std::vector<std::vector<ModularCoeff>>& t_v,          // Coefficient vectors (output)
    const std::vector<StackVect>& t_xw,                   // Border structure
    const std::vector<std::vector<int32_t>>& i_xw,        // Variable indices
    ModularCoeff prime                                    // Modular arithmetic prime
);

/**
 * Multiply quotient ring element by variable using precomputed tables
 * Mirrors Julia's _mul_var_quo! function
 * 
 * @param result Output vector
 * @param input Input vector in quotient basis
 * @param var_index Variable to multiply by
 * @param i_xw Pre-computed indices
 * @param t_v Coefficient vectors
 * @param prime Modular arithmetic prime
 */
void mul_var_quo(
    std::vector<ModularCoeff>& result,                    // Output vector
    const std::vector<ModularCoeff>& input,              // Input vector in quotient basis
    int32_t var_index,                                    // Variable to multiply by
    const std::vector<std::vector<int32_t>>& i_xw,       // Pre-computed indices
    const std::vector<std::vector<ModularCoeff>>& t_v,   // Coefficient vectors
    ModularCoeff prime                                    // Modular arithmetic prime
);

/**
 * Internal helper function for multiplication with accumulator buffer
 * Used by learn_compute_table for efficient computation
 */
bool mul_var_quo_internal(
    std::vector<ModularCoeff>& result,                    // Output vector
    const std::vector<ModularCoeff>& input,              // Input vector in quotient basis
    int32_t var_index,                                    // Variable to multiply by (1-based)
    const std::vector<std::vector<int32_t>>& i_xw,       // Pre-computed indices
    const std::vector<std::vector<ModularCoeff>>& t_v,   // Coefficient vectors
    ModularCoeff prime,                                   // Modular arithmetic prime
    std::vector<AccModularCoeff>& buf                     // Accumulator buffer
);

/**
 * Compute and fill quotient basis multiplication tables from Gröbner basis
 * Mirrors Julia's compute_fill_quo_gb! function exactly
 * 
 * This is the critical bridge function that converts F4 Gröbner basis output
 * into multiplication table coefficients. For each element in the border structure:
 * - If it's a Gröbner basis element: vectorize the polynomial in quotient basis
 * - If it's a quotient element: store direct reference
 * 
 * @param t_v Output coefficient vectors (one per border element)
 * @param t_xw Border structure from prepare_table_mxi
 * @param groebner_exponents Exponent vectors from F4 (indexed by GB position)
 * @param groebner_coefficients Coefficient vectors from F4 (indexed by GB position)
 * @param quotient_basis The quotient basis monomials
 * @param prime Modular arithmetic prime
 */
void compute_fill_quo_gb(
    std::vector<std::vector<ModularCoeff>>& t_v,
    const std::vector<StackVect>& t_xw,
    const std::vector<std::vector<PP>>& groebner_exponents,
    const std::vector<std::vector<ModularCoeff>>& groebner_coefficients,
    const std::vector<PP>& quotient_basis,
    ModularCoeff prime
);

/**
 * Vectorize a Gröbner basis polynomial in the quotient basis
 * Helper function for compute_fill_quo_gb that converts a polynomial
 * from monomial-coefficient representation to quotient basis vector
 * 
 * @param exponents Polynomial exponents (power products)
 * @param coefficients Polynomial coefficients
 * @param quotient_basis The quotient basis monomials
 * @param result Output vector in quotient basis (already allocated)
 * @param prime Modular arithmetic prime
 */
void vectorize_polynomial_in_quotient_basis(
    const std::vector<PP>& exponents,
    const std::vector<ModularCoeff>& coefficients,
    const std::vector<PP>& quotient_basis,
    std::vector<ModularCoeff>& result,
    ModularCoeff prime
);

/**
 * Initialize coefficient vectors for basic elements
 * Simplified version for testing basic functionality
 * 
 * @param t_v Output coefficient vectors
 * @param t_xw Border structure
 * @param quotient_basis_size Size of quotient basis
 */
void initialize_coefficient_vectors(
    std::vector<std::vector<ModularCoeff>>& t_v,
    const std::vector<StackVect>& t_xw,
    size_t quotient_basis_size
);

} // namespace julia_rur