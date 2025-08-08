#pragma once

#include <string>
#include <vector>
#include <map>

namespace julia_rur {

/**
 * @brief Expand polynomial expressions to standard form
 * 
 * Handles:
 * - Powers of binomials: (x-1)^2 -> x^2 - 2*x + 1
 * - Products: (x-1)*(x+1) -> x^2 - 1
 * - Nested expressions
 */
class PolynomialExpander {
public:
    /**
     * @brief Expand a polynomial expression to standard form
     * 
     * @param expression Input expression like "(x-1)^2"
     * @return Expanded form like "x^2 - 2*x + 1"
     */
    static std::string expand(const std::string& expression);
    
private:
    // Simple polynomial representation for expansion
    struct Monomial {
        double coefficient;
        std::map<std::string, int> variables; // variable -> power
    };
    
    using Polynomial = std::vector<Monomial>;
    
    // Convert expanded polynomial back to string
    static std::string polynomial_to_string(const Polynomial& poly);
    
    // Multiply two polynomials
    static Polynomial multiply(const Polynomial& p1, const Polynomial& p2);
    
    // Raise polynomial to a power
    static Polynomial power(const Polynomial& p, int n);
    
    // Parse simple expressions for expansion
    static Polynomial parse_for_expansion(const std::string& expr);
};

} // namespace julia_rur