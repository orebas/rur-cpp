#include "../src/julia_rur/univariate_parameterization.hpp"
#include "../src/julia_rur/multiplication_tables.hpp"
#include <iostream>
#include <vector>

using namespace julia_rur;

int main() {
    std::cout << "=== Testing biv_lex for Parabola-Line System ===" << std::endl;
    std::cout << "System: x^2 - y = 0, y - x - 1 = 0" << std::endl;
    std::cout << "This gives x^2 - x - 1 = 0, so x should be a separating element" << std::endl;
    std::cout << "Expected: y = x + 1 should give LINEAR parameterization" << std::endl;
    std::cout << "\n";

    // Prime
    ModularCoeff prime = 1073741827;
    
    // Quotient basis: {1, y}
    std::vector<PP> quotient_basis = {
        {0, 0},  // 1
        {0, 1}   // y
    };
    
    // From the actual test output, we have:
    // x = [1073741826, 1] = -1 + y (mod p)
    // y = [0, 1]
    // y^2 = [1073741826, 3] = -1 + 3y (mod p)
    
    // Multiplication tables (from debug output)
    // i_xw[0] corresponds to x, i_xw[1] corresponds to y
    std::vector<std::vector<int32_t>> i_xw(2);
    i_xw[0] = {2, 3};  // x * {1, y}
    i_xw[1] = {1, 4};  // y * {1, y}
    
    // t_v contains the actual multiplication results
    std::vector<std::vector<ModularCoeff>> t_v(4);
    t_v[0] = {0, 1};           // y * 1 = y = [0, 1]
    t_v[1] = {1073741826, 1};  // x * 1 = x = [-1, 1] mod p
    t_v[2] = {1073741826, 2};  // x * y = xy = [-1, 2] mod p  
    t_v[3] = {1073741826, 3};  // y * y = y^2 = [-1, 3] mod p
    
    std::cout << "Testing with y as separating element:" << std::endl;
    std::cout << "================================" << std::endl;
    
    // Use y = [0, 1] as separating element
    std::vector<ModularCoeff> separating_element_y = {0, 1};
    
    // Compute minimal polynomial for y
    MinimalPolynomialResult min_poly_y = compute_minimal_polynomial(
        separating_element_y, i_xw, t_v, quotient_basis, prime
    );
    
    std::cout << "Minimal polynomial degree: " << min_poly_y.degree << " (expected: 2)" << std::endl;
    
    if (min_poly_y.degree == 2) {
        std::cout << "✓ Minimal polynomial has correct degree" << std::endl;
        
        // Now test biv_lex for x (variable index 1)
        std::cout << "\nTesting biv_lex for variable x (index 1):" << std::endl;
        BivariateResult result_x = biv_lex(1, min_poly_y, i_xw, t_v, quotient_basis.size(), prime);
        
        if (result_x.success) {
            std::cout << "biv_lex returned success" << std::endl;
            std::cout << "Number of generators: " << result_x.generators.size() << std::endl;
            if (!result_x.generators.empty()) {
                std::cout << "First generator size: " << result_x.generators[0].size() << std::endl;
            }
            
            // Check the basis to see if it's linear
            std::cout << "Basis monomials:" << std::endl;
            for (size_t i = 0; i < result_x.basis.size(); ++i) {
                const auto& pair = result_x.basis[i];
                std::cout << "  Monomial " << i << ": T^" << pair.first[0] << " * x^" << pair.second[0] << std::endl;
                
                // Check if deg_xi > 1 (non-linear)
                if (pair.second[0] > 1) {
                    std::cout << "ERROR: Found non-linear term with deg_xi = " << pair.second[0] << std::endl;
                    std::cout << "This means y is incorrectly classified as non-separating!" << std::endl;
                }
            }
        } else {
            std::cout << "ERROR: biv_lex failed!" << std::endl;
        }
    }
    
    std::cout << "\n\nTesting with x as separating element:" << std::endl;
    std::cout << "================================" << std::endl;
    
    // Use x = [1073741826, 1] as separating element  
    std::vector<ModularCoeff> separating_element_x = {1073741826, 1};
    
    // Compute minimal polynomial for x
    MinimalPolynomialResult min_poly_x = compute_minimal_polynomial(
        separating_element_x, i_xw, t_v, quotient_basis, prime
    );
    
    std::cout << "Minimal polynomial degree: " << min_poly_x.degree << " (expected: 2)" << std::endl;
    
    if (min_poly_x.degree == 2) {
        std::cout << "✓ Minimal polynomial has correct degree" << std::endl;
        
        // Test biv_lex for y (variable index 2)
        std::cout << "\nTesting biv_lex for variable y (index 2):" << std::endl;
        BivariateResult result_y = biv_lex(2, min_poly_x, i_xw, t_v, quotient_basis.size(), prime);
        
        if (result_y.success) {
            std::cout << "biv_lex returned success" << std::endl;
            
            // Check the basis to see if it's linear
            std::cout << "Basis monomials:" << std::endl;
            for (size_t i = 0; i < result_y.basis.size(); ++i) {
                const auto& pair = result_y.basis[i];
                std::cout << "  Monomial " << i << ": T^" << pair.first[0] << " * y^" << pair.second[0] << std::endl;
                
                // Check if deg_yi > 1 (non-linear)
                if (pair.second[0] > 1) {
                    std::cout << "ERROR: Found non-linear term with deg_yi = " << pair.second[0] << std::endl;
                    std::cout << "This means x is incorrectly classified as non-separating!" << std::endl;
                }
            }
        } else {
            std::cout << "ERROR: biv_lex failed!" << std::endl;
        }
    }
    
    std::cout << "\n\nTesting with x+y as separating element:" << std::endl;
    std::cout << "=====================================" << std::endl;
    
    // x + y = [1073741826, 1] + [0, 1] = [1073741826, 2]
    std::vector<ModularCoeff> separating_element_xy = {1073741826, 2};
    
    MinimalPolynomialResult min_poly_xy = compute_minimal_polynomial(
        separating_element_xy, i_xw, t_v, quotient_basis, prime
    );
    
    std::cout << "Minimal polynomial degree: " << min_poly_xy.degree << " (expected: 2)" << std::endl;
    
    if (min_poly_xy.degree == 2) {
        std::cout << "✓ Minimal polynomial has correct degree" << std::endl;
        
        // Test both variables
        for (int var = 1; var <= 2; ++var) {
            std::cout << "\nTesting biv_lex for variable " << (var == 1 ? "x" : "y") << ":" << std::endl;
            BivariateResult result = biv_lex(var, min_poly_xy, i_xw, t_v, quotient_basis.size(), prime);
            
            if (result.success) {
                bool is_linear = true;
                for (const auto& pair : result.basis) {
                    if (pair.second[0] > 1) {
                        is_linear = false;
                        std::cout << "ERROR: Non-linear term with deg = " << pair.second[0] << std::endl;
                    }
                }
                if (is_linear) {
                    std::cout << "✓ LINEAR parameterization found!" << std::endl;
                }
            } else {
                std::cout << "ERROR: biv_lex failed!" << std::endl;
            }
        }
    }
    
    return 0;
}