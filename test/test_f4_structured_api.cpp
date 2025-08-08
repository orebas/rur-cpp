#include "../src/axf4_wrapper.h"
#include <iostream>
#include <vector>
#include <cassert>

void test_structured_api() {
    std::cout << "Testing F4 structured data extraction API..." << std::endl;
    
    // Create a simple session
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(100003, vars, 2);
    assert(session != nullptr);
    
    // Add simple polynomials in F4 format: x^2 + y^2 - 1, x - y
    // x^2 + y^2 - 1 -> 1*x^2+1*y^2+100002
    // x - y -> 1*x+100002*y (since -1 = 100002 mod 100003)
    assert(axf4_add_polynomial(session, "1*x^2+1*y^2+100002") == 0);
    assert(axf4_add_polynomial(session, "1*x+100002*y") == 0);
    
    std::cout << "  Computing Gröbner basis..." << std::endl;
    axf4_result_t result = axf4_compute_groebner_basis_keep_data(session);
    assert(result.status == 0);
    std::cout << "  Basis computed with " << result.basis_size << " polynomials" << std::endl;
    
    // Test structured API
    std::cout << "  Testing structured API functions..." << std::endl;
    
    // Get basis size
    int basis_size = axf4_get_basis_size();
    std::cout << "    Basis size: " << basis_size << std::endl;
    assert(basis_size > 0);
    assert(basis_size == result.basis_size);
    
    // Test each polynomial
    for (int i = 0; i < basis_size; i++) {
        std::cout << "    Polynomial " << i << ":" << std::endl;
        
        // Get term count
        int term_count = axf4_get_poly_term_count(i);
        std::cout << "      Terms: " << term_count << std::endl;
        assert(term_count > 0);
        
        // Get polynomial data
        std::vector<unsigned int> coeffs(term_count);
        std::vector<unsigned int> monomials(term_count);
        
        int ret = axf4_get_poly_data(i, coeffs.data(), monomials.data());
        assert(ret == 0);
        
        std::cout << "      Data: ";
        for (int j = 0; j < term_count; j++) {
            std::cout << "[c:" << coeffs[j] << " m:" << monomials[j] << "] ";
        }
        std::cout << std::endl;
        
        // Get leading term
        unsigned int lead_coeff, lead_monomial;
        ret = axf4_get_leading_term(i, &lead_coeff, &lead_monomial);
        assert(ret == 0);
        assert(lead_coeff == coeffs[term_count - 1]);  // Should match last coefficient (F4 sorts ascending)
        assert(lead_monomial == monomials[term_count - 1]);  // Should match last monomial
        
        std::cout << "      Leading term: [c:" << lead_coeff << " m:" << lead_monomial << "]" << std::endl;
    }
    
    // Get all leading monomials
    std::vector<unsigned int> leading_monomials(basis_size);
    int lead_count = axf4_get_all_leading_monomials(leading_monomials.data());
    assert(lead_count == basis_size);
    
    std::cout << "    All leading monomials: ";
    for (int i = 0; i < lead_count; i++) {
        std::cout << leading_monomials[i] << " ";
    }
    std::cout << std::endl;
    
    // Test error cases
    std::cout << "  Testing error cases..." << std::endl;
    
    // Out of bounds access
    assert(axf4_get_poly_term_count(-1) == -1);
    assert(axf4_get_poly_term_count(basis_size) == -1);
    
    unsigned int dummy_coeff, dummy_monomial;
    assert(axf4_get_leading_term(-1, &dummy_coeff, &dummy_monomial) == -1);
    assert(axf4_get_leading_term(basis_size, &dummy_coeff, &dummy_monomial) == -1);
    
    // Null pointer checks  
    std::vector<unsigned int> test_coeffs(10);
    std::vector<unsigned int> test_monomials(10);
    assert(axf4_get_poly_data(0, nullptr, test_monomials.data()) == -1);
    assert(axf4_get_poly_data(0, test_coeffs.data(), nullptr) == -1);
    assert(axf4_get_leading_term(0, nullptr, &dummy_monomial) == -1);
    assert(axf4_get_leading_term(0, &dummy_coeff, nullptr) == -1);
    assert(axf4_get_all_leading_monomials(nullptr) == -1);
    
    std::cout << "  ✓ Error handling works correctly" << std::endl;
    
    // Cleanup structured data first, then result and session
    axf4_cleanup_basis_data();
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    std::cout << "✓ Structured API test passed!" << std::endl;
}

void test_empty_basis() {
    std::cout << "Testing API with no computed basis..." << std::endl;
    
    // Test API functions when no basis is computed
    assert(axf4_get_basis_size() == 0);  // Should be 0, not -1, when no computation done
    assert(axf4_get_poly_term_count(0) == -1);  // Should fail
    
    unsigned int dummy;
    assert(axf4_get_all_leading_monomials(&dummy) == -1);  // Should fail
    
    std::cout << "✓ Empty basis test passed!" << std::endl;
}

int main() {
    try {
        test_empty_basis();
        test_structured_api();
        
        std::cout << "\n🎉 All structured API tests passed!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with error: " << e.what() << std::endl;
        return 1;
    }
}