#include <iostream>
#include <vector>
#include "julia_rur/univariate_parameterization.hpp"
#include "julia_rur/multiplication_tables.hpp"
#include "julia_rur/f4_integration.hpp"
#include "julia_rur/quotient_basis.hpp"

using namespace julia_rur;

int main() {
    std::cout << "=== Testing Minimal Polynomial Implementations ===" << std::endl;
    
    // Test with the 3-variable symmetric system
    // We'll use a simple test case with known multiplication tables
    
    ModularCoeff prime = 1073741827;
    std::cout << "Prime: " << prime << std::endl;
    
    // Quotient basis: {1, z, y, z^2, yz, yz^2}
    std::vector<PP> quotient_basis = {
        {0, 0, 0},  // 1
        {0, 0, 1},  // z
        {0, 1, 0},  // y
        {0, 0, 2},  // z^2
        {0, 1, 1},  // yz
        {0, 1, 2}   // yz^2
    };
    
    // Separating element: -2x - y + z
    // In quotient basis representation (from our debug output):
    // [1073741815, 3, 1, 0, 0, 0] mod 1073741827
    std::vector<ModularCoeff> separating_element = {1073741815, 3, 1, 0, 0, 0};
    
    // For this test, we need actual multiplication tables
    // These would normally come from F4, but for testing we can use dummy values
    // or load from a previous computation
    
    std::vector<std::vector<int32_t>> i_xw;
    std::vector<std::vector<ModularCoeff>> t_v;
    
    // Initialize dummy multiplication tables for testing
    // In a real test, we'd load these from F4 output
    // For now, just create minimal structure to avoid crashes
    i_xw.resize(3);  // 3 variables
    for (int i = 0; i < 3; ++i) {
        i_xw[i].resize(6);  // 6 basis elements
        for (int j = 0; j < 6; ++j) {
            i_xw[i][j] = j + 1;  // Dummy indices
        }
    }
    
    t_v.resize(18);  // 3 variables * 6 basis elements
    for (size_t i = 0; i < t_v.size(); ++i) {
        t_v[i].resize(6, 0);  // Initialize with zeros
        t_v[i][0] = 1;  // Some non-zero values
    }
    
    std::cout << "\nSeparating element: [";
    for (size_t i = 0; i < separating_element.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << separating_element[i];
    }
    std::cout << "]" << std::endl;
    
    // Test OLD implementation
    std::cout << "\n--- Testing OLD implementation ---" << std::endl;
    auto result_old = compute_minimal_polynomial(
        separating_element, i_xw, t_v, quotient_basis, prime
    );
    
    if (result_old.success) {
        std::cout << "OLD: Success! Degree = " << result_old.degree << std::endl;
        std::cout << "OLD: Coefficients = [";
        for (size_t i = 0; i < result_old.coefficients.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << result_old.coefficients[i];
        }
        std::cout << "]" << std::endl;
    } else {
        std::cout << "OLD: FAILED" << std::endl;
    }
    
    // Test NEW implementation
    std::cout << "\n--- Testing NEW implementation ---" << std::endl;
    auto result_new = compute_minimal_polynomial_flint(
        separating_element, i_xw, t_v, quotient_basis, prime
    );
    
    if (result_new.success) {
        std::cout << "NEW: Success! Degree = " << result_new.degree << std::endl;
        std::cout << "NEW: Coefficients = [";
        for (size_t i = 0; i < result_new.coefficients.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << result_new.coefficients[i];
        }
        std::cout << "]" << std::endl;
    } else {
        std::cout << "NEW: FAILED" << std::endl;
    }
    
    // Compare results
    std::cout << "\n--- Comparison ---" << std::endl;
    if (result_old.success && result_new.success) {
        if (result_old.degree == result_new.degree) {
            std::cout << "Degrees match: " << result_old.degree << std::endl;
        } else {
            std::cout << "MISMATCH: Degrees differ! OLD=" << result_old.degree 
                      << " NEW=" << result_new.degree << std::endl;
        }
        
        bool coeffs_match = true;
        if (result_old.coefficients.size() == result_new.coefficients.size()) {
            for (size_t i = 0; i < result_old.coefficients.size(); ++i) {
                if (result_old.coefficients[i] != result_new.coefficients[i]) {
                    coeffs_match = false;
                    std::cout << "  Coeff[" << i << "] differs: OLD=" 
                              << result_old.coefficients[i] 
                              << " NEW=" << result_new.coefficients[i] << std::endl;
                }
            }
        } else {
            coeffs_match = false;
            std::cout << "Coefficient array sizes differ!" << std::endl;
        }
        
        if (coeffs_match) {
            std::cout << "Coefficients match!" << std::endl;
        }
    }
    
    // Expected Julia coefficients for reference
    std::cout << "\n--- Expected (Julia) ---" << std::endl;
    std::cout << "Julia coefficients: [1260, 2952, 2545, 1056, 226, 24, 1]" << std::endl;
    
    return 0;
}