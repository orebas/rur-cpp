#include "../src/julia_rur/rur_main_algorithm.hpp"
#include "../src/julia_rur/multiplication_tables.hpp"
#include <iostream>

using namespace julia_rur;

void test_mult_table_for_prime(ModularCoeff prime) {
    std::cout << "\n=== Testing multiplication table for prime " << prime << " ===" << std::endl;
    
    std::vector<std::string> polynomials = {"1*x^2-2"};
    std::vector<std::string> variables = {"x"};
    
    // Create F4 session
    std::vector<const char*> var_ptrs;
    for (const auto& var : variables) {
        var_ptrs.push_back(var.c_str());
    }
    
    axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), variables.size());
    
    // Add polynomial
    axf4_add_polynomial(session, polynomials[0].c_str());
    
    // Compute GB
    axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
    std::cout << "GB: " << gb_result.groebner_basis << std::endl;
    
    // Extract multiplication tables
    std::vector<std::vector<ModularCoeff>> t_v;
    std::vector<StackVect> t_xw;
    std::vector<std::vector<int32_t>> i_xw;
    std::vector<PP> quotient_basis;
    
    bool success = f4_to_multiplication_tables(
        session, t_v, t_xw, i_xw, quotient_basis, prime
    );
    
    std::cout << "Quotient basis size: " << quotient_basis.size() << std::endl;
    
    if (success && quotient_basis.size() == 2) {
        std::cout << "Multiplication table t_v:" << std::endl;
        for (size_t i = 0; i < t_v.size(); ++i) {
            std::cout << "  Row " << i << ": [";
            for (size_t j = 0; j < t_v[i].size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << t_v[i][j];
                if (t_v[i][j] > prime/2) {
                    std::cout << "(" << static_cast<int64_t>(t_v[i][j]) - prime << ")";
                }
            }
            std::cout << "]" << std::endl;
        }
        
        // Test multiplication: x * x
        std::cout << "\nComputing x * x in quotient ring:" << std::endl;
        std::vector<ModularCoeff> x_vec(2, 0);
        x_vec[1] = 1;  // x = [0, 1]
        
        std::vector<ModularCoeff> result(2, 0);
        mul_var_quo(result, x_vec, 1, i_xw, t_v, prime);
        
        std::cout << "x * x = [" << result[0];
        if (result[0] > prime/2) {
            std::cout << "(" << static_cast<int64_t>(result[0]) - prime << ")";
        }
        std::cout << ", " << result[1] << "]" << std::endl;
    }
    
    axf4_free_result(&gb_result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);
}

int main() {
    test_mult_table_for_prime(131063);
    test_mult_table_for_prime(131059);
    
    return 0;
}