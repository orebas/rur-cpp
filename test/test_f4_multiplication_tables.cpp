#include "../src/julia_rur/f4_integration.hpp"
#include "../src/axf4_wrapper.h"
#include <iostream>
#include <cassert>

using namespace julia_rur;

void test_f4_to_multiplication_tables() {
    std::cout << "\n=== Testing F4 to Multiplication Tables Pipeline ===" << std::endl;
    
    // Create F4 session
    const char* variables[] = {"x"};
    axf4_session_t session = axf4_create_session(100003, variables, 1);
    assert(session != nullptr);
    
    // Add polynomial x^2 (should give quotient basis {1, x})
    int result = axf4_add_polynomial(session, "1*x^2");
    assert(result == 0);
    
    // Compute Gröbner basis with keep_data
    axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
    assert(gb_result.status == 0);
    assert(gb_result.basis_size == 1);
    
    std::cout << "F4 computation successful, basis size: " << gb_result.basis_size << std::endl;
    
    // Test the complete pipeline
    std::vector<std::vector<ModularCoeff>> t_v;
    std::vector<StackVect> t_xw;
    std::vector<std::vector<int32_t>> i_xw;
    std::vector<PP> quotient_basis;
    
    bool success = f4_to_multiplication_tables(session, t_v, t_xw, i_xw, quotient_basis, 100003);
    assert(success);
    
    std::cout << "Pipeline successful!" << std::endl;
    std::cout << "Quotient basis size: " << quotient_basis.size() << std::endl;
    std::cout << "Border structure size: " << t_xw.size() << std::endl;
    std::cout << "Coefficient vectors size: " << t_v.size() << std::endl;
    
    // Validate quotient basis
    assert(quotient_basis.size() == 2);  // For x^2, quotient basis should be {1, x}
    assert(quotient_basis[0] == PP({0})); // The constant monomial
    assert(quotient_basis[1] == PP({1})); // The x monomial
    
    // Print some details
    std::cout << "\nQuotient basis: ";
    for (const auto& pp : quotient_basis) {
        std::cout << "deg=" << pp[0] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "\nBorder structure:" << std::endl;
    for (size_t i = 0; i < t_xw.size(); ++i) {
        std::cout << "  [" << i << "] pos=" << t_xw[i].pos 
                  << " var=" << t_xw[i].var 
                  << " prev=" << t_xw[i].prev 
                  << " mon=deg" << (t_xw[i].mon.empty() ? -1 : (int)t_xw[i].mon[0])
                  << " t_v_size=" << t_v[i].size() << std::endl;
    }
    
    // Cleanup
    axf4_free_result(&gb_result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);
    
    std::cout << "✓ F4 to multiplication tables test passed!" << std::endl;
}

void test_two_variable_system() {
    std::cout << "\n=== Testing Two Variable System ===" << std::endl;
    
    // Create F4 session for system {x^2, y^2} 
    const char* variables[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(100003, variables, 2);
    assert(session != nullptr);
    
    // Add polynomials x^2 and y^2
    assert(axf4_add_polynomial(session, "1*x^2") == 0);
    assert(axf4_add_polynomial(session, "1*y^2") == 0);
    
    // Compute Gröbner basis
    axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
    assert(gb_result.status == 0);
    
    std::cout << "F4 computation successful, basis size: " << gb_result.basis_size << std::endl;
    
    // Test the pipeline
    std::vector<std::vector<ModularCoeff>> t_v;
    std::vector<StackVect> t_xw;
    std::vector<std::vector<int32_t>> i_xw;
    std::vector<PP> quotient_basis;
    
    bool success = f4_to_multiplication_tables(session, t_v, t_xw, i_xw, quotient_basis, 100003);
    assert(success);
    
    std::cout << "Pipeline successful!" << std::endl;
    std::cout << "Quotient basis size: " << quotient_basis.size() << std::endl;
    
    // For {x^2, y^2}, quotient basis should be {1, x, y, xy}
    assert(quotient_basis.size() == 4);
    
    std::cout << "\nQuotient basis:" << std::endl;
    for (size_t i = 0; i < quotient_basis.size(); ++i) {
        const auto& pp = quotient_basis[i];
        std::cout << "  [" << i << "] x^" << pp[0] << "*y^" << pp[1] << std::endl;
    }
    
    // Cleanup
    axf4_free_result(&gb_result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);
    
    std::cout << "✓ Two variable system test passed!" << std::endl;
}

int main() {
    std::cout << "=== Testing F4 to Multiplication Tables Integration ===" << std::endl;
    
    try {
        test_f4_to_multiplication_tables();
        test_two_variable_system();
        
        std::cout << "\n✅ All F4 multiplication table tests passed!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Test failed with unknown exception" << std::endl;
        return 1;
    }
}