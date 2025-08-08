#include <iostream>
#include <vector>
#include "../src/julia_rur/bivariate_algorithm.hpp"
#include "../src/julia_rur/multiplication_tables.hpp"
#include "../src/julia_rur/polynomial_operations.hpp"

using namespace julia_rur;

int main() {
    std::cout << "=== Testing Univariate Fix ===" << std::endl;
    
    // For x^3 - 8 = 0, the quotient ring basis is {1, x, x^2}
    // and x^3 = 8 in the quotient ring
    
    ModularCoeff prime = 1073741827;
    
    // Multiplication table for the quotient ring Q[x]/(x^3 - 8)
    // We have basis {1, x, x^2}
    // x * 1 = x (basis element 1)
    // x * x = x^2 (basis element 2) 
    // x * x^2 = x^3 = 8 (reduces to 8*1 = basis element 0 with coeff 8)
    
    std::vector<std::vector<ModularCoeff>> t_v;
    t_v.push_back({0, 1, 0});  // x * 1 = x
    t_v.push_back({0, 0, 1});  // x * x = x^2
    t_v.push_back({8, 0, 0});  // x * x^2 = 8
    
    std::vector<std::vector<int32_t>> i_xw = {{1, 2, 3}}; // Indices for x multiplication
    
    // The minimal polynomial for T = x is T^3 - 8
    // So we have powers T^0 = 1, T^1 = x, T^2 = x^2
    std::vector<std::vector<ModularCoeff>> free_set;
    free_set.push_back({1, 0, 0});  // T^0 = 1
    free_set.push_back({0, 1, 0});  // T^1 = x
    free_set.push_back({0, 0, 1});  // T^2 = x^2
    
    // Run biv_lex to get the relation for x in terms of T
    GaussianReductionContext context;
    BivLexResult result = biv_lex(t_v, i_xw, context, free_set, 1, prime);
    
    std::cout << "\nbiv_lex results:" << std::endl;
    std::cout << "Number of generators: " << result.generators.size() << std::endl;
    std::cout << "Number of leading monomials: " << result.leading_monomials.size() << std::endl;
    
    if (!result.generators.empty()) {
        std::cout << "\nFirst generator coefficients: ";
        for (auto c : result.generators[0]) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        
        if (!result.leading_monomials.empty()) {
            auto &lm = result.leading_monomials[0];
            std::cout << "Leading monomial: T^" << lm.deg_T << " * x^" << lm.deg_xi << std::endl;
            
            // For univariate case with T = x, we expect:
            // The relation x - T = 0, which means:
            // - Leading monomial should be x^1 (deg_T=0, deg_xi=1)
            // - Generator should encode the coefficients [T, -1] meaning T - x = 0
            // But in the quotient ring representation, T is represented as [0,1,0]
            
            std::cout << "\nExpected for x = T:" << std::endl;
            std::cout << "  Leading monomial: T^0 * x^1" << std::endl;
            std::cout << "  Generator should represent x - T = 0" << std::endl;
            std::cout << "  In basis {1, x, x^2}, T is [0,1,0]" << std::endl;
            std::cout << "  So generator should be [0,1,0] for the T term" << std::endl;
        }
    }
    
    return 0;
}