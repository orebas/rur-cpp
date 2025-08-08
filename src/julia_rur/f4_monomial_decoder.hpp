#ifndef JULIA_RUR_F4_MONOMIAL_DECODER_HPP
#define JULIA_RUR_F4_MONOMIAL_DECODER_HPP

#include "data_structures.hpp"
#include <vector>
#include <iostream>
#include <cstdint>

// Forward declare F4 internal structures we need access to
extern "C" {
    extern int f4_nvars;
    // F4 uses platform-specific INT type (int32 or int64)
    #if UINTPTR_MAX==0xFFFFFFFF
        typedef int INT;
    #else
        typedef long long int INT;
    #endif
    extern INT* f4_monom;
    // EXPON macro from F4: returns pointer to exponent vector
    #define F4_EXPON(m) (f4_monom + (m)*(f4_nvars+1))
}

namespace julia_rur {

/**
 * @brief Decode F4 monomial index to our power product (PP) format
 * 
 * F4 stores monomials as indices into a hash table. The actual exponent
 * vectors are stored in the f4_monom array, accessed via EXPON(m).
 * 
 * @param monomial_index F4 monomial index (as returned by structured API)
 * @param num_variables Number of variables in the polynomial ring
 * @return PP Power product with exponents for each variable
 */
inline PP decode_f4_monomial(unsigned int monomial_index, int num_variables) {
    PP result(num_variables, 0);
    
    // Get pointer to exponent vector from F4's internal storage
    INT* exponents = F4_EXPON(monomial_index);
    
    // Debug output disabled - decoder is working correctly now
    // std::cout << "decode_f4_monomial(" << monomial_index << "): f4_nvars=" << f4_nvars 
    //           << ", exponents at " << (void*)exponents << " = ";
    // for (int i = 0; i < num_variables; i++) {
    //     std::cout << exponents[i] << " ";
    // }
    // std::cout << std::endl;
    
    // Copy exponents, converting from F4's int to our uint32_t
    for (int i = 0; i < num_variables; i++) {
        if (exponents[i] < 0) {
            throw std::runtime_error("Negative exponent in F4 monomial");
        }
        result[i] = static_cast<uint32_t>(exponents[i]);
    }
    
    return result;
}

/**
 * @brief Extract leading monomials from F4 Gröbner basis
 * 
 * Uses F4 structured API to get leading terms and decode them to PP format.
 * Assumes axf4_compute_groebner_basis_keep_data() has been called.
 * 
 * @param basis_size Number of polynomials in the basis (from axf4_get_basis_size())
 * @param num_variables Number of variables in the polynomial ring
 * @return Vector of leading monomials as PPs
 */
std::vector<PP> extract_f4_leading_monomials(int basis_size, int num_variables);

} // namespace julia_rur

#endif // JULIA_RUR_F4_MONOMIAL_DECODER_HPP