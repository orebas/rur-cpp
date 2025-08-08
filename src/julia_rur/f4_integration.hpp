#ifndef JULIA_RUR_F4_INTEGRATION_HPP
#define JULIA_RUR_F4_INTEGRATION_HPP

#include "data_structures.hpp"
#include "quotient_basis.hpp"
#include "multiplication_tables.hpp"
#include "f4_monomial_decoder.hpp"
#include "../axf4_wrapper.h"
#include <vector>

namespace julia_rur {

/**
 * Extract structured Gröbner basis data from F4
 * Converts F4's internal representation to our polynomial format
 * 
 * @param groebner_exponents Output: exponent vectors for each GB polynomial
 * @param groebner_coefficients Output: coefficient vectors for each GB polynomial  
 * @param num_variables Number of variables
 * @return true if successful, false on error
 */
bool extract_f4_groebner_basis(
    std::vector<std::vector<PP>>& groebner_exponents,
    std::vector<std::vector<ModularCoeff>>& groebner_coefficients,
    int num_variables
);

/**
 * Complete F4 to multiplication tables pipeline
 * Integrates all steps: F4 → quotient basis → multiplication tables
 * 
 * @param session F4 session (must have computed GB with keep_data)
 * @param t_v Output: coefficient vectors for multiplication tables
 * @param t_xw Output: border structure
 * @param i_xw Output: variable indices
 * @param quotient_basis Output: quotient basis
 * @param prime Modular arithmetic prime
 * @return true if successful, false on error
 */
bool f4_to_multiplication_tables(
    axf4_session_t session,
    std::vector<std::vector<ModularCoeff>>& t_v,
    std::vector<StackVect>& t_xw,
    std::vector<std::vector<int32_t>>& i_xw,
    std::vector<PP>& quotient_basis,
    ModularCoeff prime
);

} // namespace julia_rur

#endif // JULIA_RUR_F4_INTEGRATION_HPP