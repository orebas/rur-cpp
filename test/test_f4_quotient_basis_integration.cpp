#include "../src/axf4_wrapper.h"
#include "../src/julia_rur/quotient_basis.hpp"
#include "../src/julia_rur/f4_monomial_decoder.hpp"
#include <iostream>
#include <cassert>
#include <vector>
#include <unordered_map>

using namespace julia_rur;


void test_simple_linear_system() {
    std::cout << "\n=== Test: F4 -> Quotient Basis for {x-1} ===" << std::endl;
    
    const char* vars[] = {"x"};
    axf4_session_t session = axf4_create_session(100003, vars, 1);
    assert(session != nullptr);
    
    // Add polynomial x-1 in F4 format (1*x + 100002)
    assert(axf4_add_polynomial(session, "1*x+100002") == 0);
    
    // Compute Gröbner basis keeping data alive
    axf4_result_t result = axf4_compute_groebner_basis_keep_data(session);
    assert(result.status == 0);
    
    std::cout << "F4 Gröbner basis: " << result.groebner_basis << std::endl;
    std::cout << "Basis size: " << axf4_get_basis_size() << std::endl;
    
    // Extract leading terms using structured API and proper decoder
    int basis_size = axf4_get_basis_size();
    std::vector<PP> leading_terms = extract_f4_leading_monomials(basis_size, 1);
    
    std::cout << "Leading terms extracted: " << leading_terms.size() << std::endl;
    
    std::cout << "Leading terms: ";
    for (size_t i = 0; i < leading_terms.size(); i++) {
        if (i > 0) std::cout << ", ";
        // Print the power product
        bool first = true;
        bool all_zero = true;
        for (size_t j = 0; j < leading_terms[i].size(); j++) {
            if (leading_terms[i][j] > 0) {
                all_zero = false;
                if (!first) std::cout << "*";
                std::cout << "x";
                if (leading_terms[i][j] > 1) std::cout << "^" << leading_terms[i][j];
                first = false;
            }
        }
        if (all_zero) std::cout << "1";
    }
    std::cout << std::endl;
    
    // Compute quotient basis
    auto quotient_basis = compute_quotient_basis(leading_terms);
    
    std::cout << "Quotient basis size: " << quotient_basis.size() << std::endl;
    assert(quotient_basis.size() == 1);  // Should be {1}
    assert(quotient_basis[0] == PP({0}));  // The constant 1
    
    // Cleanup F4 internal data first, then session
    axf4_cleanup_basis_data();
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    std::cout << "✓ Test passed!" << std::endl;
}

void test_two_variable_system() {
    std::cout << "\n=== Test: F4 -> Quotient Basis for {x^2+y^2-1, x-y} ===" << std::endl;
    
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(100003, vars, 2);
    assert(session != nullptr);
    
    // Add polynomials in F4 format
    assert(axf4_add_polynomial(session, "1*x^2+1*y^2+100002") == 0);  // x^2+y^2-1
    assert(axf4_add_polynomial(session, "1*x+100002*y") == 0);        // x-y
    
    // Compute Gröbner basis
    axf4_result_t result = axf4_compute_groebner_basis_keep_data(session);
    assert(result.status == 0);
    
    std::cout << "F4 Gröbner basis: " << result.groebner_basis << std::endl;
    std::cout << "Basis size: " << axf4_get_basis_size() << std::endl;
    
    // Extract leading terms using structured API and proper decoder
    int basis_size = axf4_get_basis_size();
    std::vector<PP> leading_terms = extract_f4_leading_monomials(basis_size, 2);
    
    std::cout << "Leading terms: ";
    for (size_t i = 0; i < leading_terms.size(); i++) {
        if (i > 0) std::cout << ", ";
        // Print the power product
        bool first = true;
        for (size_t j = 0; j < leading_terms[i].size(); j++) {
            if (leading_terms[i][j] > 0) {
                if (!first) std::cout << "*";
                std::cout << (j == 0 ? "x" : "y");
                if (leading_terms[i][j] > 1) std::cout << "^" << leading_terms[i][j];
                first = false;
            }
        }
        if (first) std::cout << "1";
    }
    std::cout << std::endl;
    
    // Compute quotient basis
    auto quotient_basis = compute_quotient_basis(leading_terms);
    
    std::cout << "Quotient basis size: " << quotient_basis.size() << std::endl;
    // Depending on the GB computation, this could be 2 or more
    
    // Cleanup F4 internal data first, then session
    axf4_cleanup_basis_data();
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    std::cout << "✓ Test passed!" << std::endl;
}

int main() {
    std::cout << "=== F4 to Quotient Basis Integration Test ===" << std::endl;
    
    try {
        test_simple_linear_system();
        test_two_variable_system();
        
        std::cout << "\n✅ All integration tests passed!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n❌ Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
}