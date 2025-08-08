#include "../src/julia_rur/f4_integration.hpp"
#include "../src/axf4_wrapper.h"
#include <iostream>

using namespace julia_rur;

void debug_f4_output(int prime) {
    std::cout << "\n=== Debugging F4 output for prime " << prime << " ===" << std::endl;
    
    const char* vars[] = {"x"};
    axf4_session_t session = axf4_create_session(prime, vars, 1);
    axf4_add_polynomial(session, "1*x^2-2");
    
    axf4_result_t result = axf4_compute_groebner_basis_keep_data(session);
    std::cout << "GB string: " << result.groebner_basis << std::endl;
    
    // Extract polynomial data
    int basis_size = axf4_get_basis_size();
    std::cout << "Basis size: " << basis_size << std::endl;
    
    for (int poly_idx = 0; poly_idx < basis_size; ++poly_idx) {
        int term_count = axf4_get_poly_term_count(poly_idx);
        std::cout << "\nPolynomial " << poly_idx << " has " << term_count << " terms:" << std::endl;
        
        std::vector<unsigned int> coeffs(term_count);
        std::vector<unsigned int> monomials(term_count);
        axf4_get_poly_data(poly_idx, coeffs.data(), monomials.data());
        
        for (int i = 0; i < term_count; ++i) {
            PP pp = decode_f4_monomial(monomials[i], 1);
            std::cout << "  Term " << i << ": coeff=" << coeffs[i];
            if (coeffs[i] > prime/2) {
                std::cout << " (" << static_cast<int>(coeffs[i]) - prime << ")";
            }
            std::cout << ", monomial=" << monomials[i] << " -> PP=[" << pp[0] << "]" << std::endl;
        }
    }
    
    // Now test the full multiplication table building
    std::vector<std::vector<ModularCoeff>> t_v;
    std::vector<StackVect> t_xw;
    std::vector<std::vector<int32_t>> i_xw;
    std::vector<PP> quotient_basis;
    
    bool success = f4_to_multiplication_tables(session, t_v, t_xw, i_xw, quotient_basis, prime);
    
    std::cout << "\nMultiplication table t_v:" << std::endl;
    for (size_t i = 0; i < t_v.size(); ++i) {
        std::cout << "  Row " << i << ": [";
        for (size_t j = 0; j < t_v[i].size(); ++j) {
            if (j > 0) std::cout << ", ";
            std::cout << t_v[i][j];
            if (t_v[i][j] > prime/2) {
                std::cout << "(" << static_cast<int>(t_v[i][j]) - prime << ")";
            }
        }
        std::cout << "]" << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);
}

int main() {
    debug_f4_output(131063);
    debug_f4_output(131059);
    return 0;
}