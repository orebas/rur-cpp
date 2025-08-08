#include "../src/axf4_wrapper.h"
#include <iostream>
#include <vector>
#include <cassert>

int main() {
    std::cout << "Simple F4 structured API test..." << std::endl;
    
    // Test with no basis computed first
    std::cout << "Testing with no basis computed:" << std::endl;
    int size = axf4_get_basis_size();
    std::cout << "Basis size (should be 0): " << size << std::endl;
    
    // Create session with a trivial polynomial
    std::cout << "Creating session..." << std::endl;
    const char* vars[] = {"x"};
    axf4_session_t session = axf4_create_session(100003, vars, 1);
    if (!session) {
        std::cerr << "Failed to create session" << std::endl;
        return 1;
    }
    
    // Add a very simple polynomial: x - 1 (in F4 format: 1*x + 100002)
    std::cout << "Adding polynomial: 1*x+100002 (x-1 in modular form)" << std::endl;
    int ret = axf4_add_polynomial(session, "1*x+100002");
    if (ret != 0) {
        std::cerr << "Failed to add polynomial: " << axf4_get_last_error(session) << std::endl;
        axf4_destroy_session(session);
        return 1;
    }
    
    std::cout << "Computing Gröbner basis..." << std::endl;
    axf4_result_t result = axf4_compute_groebner_basis_keep_data(session);
    
    if (result.status != 0) {
        std::cerr << "Gröbner basis computation failed: " << result.error_message << std::endl;
        axf4_free_result(&result);
        axf4_destroy_session(session);
        return 1;
    }
    
    std::cout << "Basis computed successfully!" << std::endl;
    std::cout << "Result basis size: " << result.basis_size << std::endl;
    std::cout << "String result: " << (result.groebner_basis ? result.groebner_basis : "NULL") << std::endl;
    
    // Test structured API
    size = axf4_get_basis_size();
    std::cout << "Structured API basis size: " << size << std::endl;
    
    if (size > 0) {
        std::cout << "Testing first polynomial:" << std::endl;
        
        int term_count = axf4_get_poly_term_count(0);
        std::cout << "  Term count: " << term_count << std::endl;
        
        if (term_count > 0) {
            std::vector<unsigned int> coeffs(term_count);
            std::vector<unsigned int> monomials(term_count);
            
            ret = axf4_get_poly_data(0, coeffs.data(), monomials.data());
            std::cout << "  Get data result: " << ret << std::endl;
            
            if (ret == 0) {
                std::cout << "  Polynomial data:" << std::endl;
                for (int i = 0; i < term_count; i++) {
                    std::cout << "    Term " << i << ": coeff=" << coeffs[i] << " monomial=" << monomials[i] << std::endl;
                }
                
                unsigned int lead_coeff, lead_monomial;
                ret = axf4_get_leading_term(0, &lead_coeff, &lead_monomial);
                std::cout << "  Leading term result: " << ret << std::endl;
                if (ret == 0) {
                    std::cout << "  Leading term: coeff=" << lead_coeff << " monomial=" << lead_monomial << std::endl;
                }
            }
        }
    }
    
    // Cleanup structured data first, then result and session
    axf4_cleanup_basis_data();
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    std::cout << "✓ Simple structured API test completed!" << std::endl;
    return 0;
}