#include "f4_monomial_decoder.hpp"
#include "../axf4_wrapper.h"
#include <stdexcept>

namespace julia_rur {

std::vector<PP> extract_f4_leading_monomials(int basis_size, int num_variables) {
    std::vector<PP> leading_monomials;
    leading_monomials.reserve(basis_size);
    
    for (int i = 0; i < basis_size; i++) {
        unsigned int lead_coeff, lead_monomial;
        
        if (axf4_get_leading_term(i, &lead_coeff, &lead_monomial) != 0) {
            throw std::runtime_error("Failed to get leading term for polynomial " + std::to_string(i));
        }
        
        // Decode the F4 monomial index to our PP format
        PP lt = decode_f4_monomial(lead_monomial, num_variables);
        leading_monomials.push_back(lt);
    }
    
    return leading_monomials;
}

} // namespace julia_rur