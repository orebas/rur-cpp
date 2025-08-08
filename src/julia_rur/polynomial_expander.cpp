#include "polynomial_expander.hpp"
#include <regex>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <cmath>

namespace julia_rur {

std::string PolynomialExpander::expand(const std::string& expression) {
    // For now, implement a simple expansion for common patterns
    // This is a basic implementation that handles (x±a)^n patterns
    
    std::string expr = expression;
    // Remove spaces
    expr.erase(std::remove_if(expr.begin(), expr.end(), ::isspace), expr.end());
    
    // Pattern for (var ± number)^power
    std::regex binomial_power(R"(\(([a-zA-Z][0-9]*)([-+])([0-9]+(?:\.[0-9]+)?)\)\^([0-9]+))");
    std::smatch match;
    
    if (std::regex_match(expr, match, binomial_power)) {
        std::string var = match[1];
        std::string op = match[2];
        double constant = std::stod(match[3]);
        int power = std::stoi(match[4]);
        
        if (op == "-") constant = -constant;
        
        // Use binomial expansion
        // (x + a)^n = sum(C(n,k) * x^(n-k) * a^k) for k=0 to n
        std::ostringstream result;
        bool first = true;
        
        for (int k = 0; k <= power; ++k) {
            // Calculate binomial coefficient C(n,k)
            int coeff = 1;
            for (int i = 1; i <= k; ++i) {
                coeff = coeff * (power - i + 1) / i;
            }
            
            // Calculate a^k
            double a_power = std::pow(constant, k);
            double term_coeff = coeff * a_power;
            
            // Skip if coefficient is effectively zero
            if (std::abs(term_coeff) < 1e-10) continue;
            
            // Build the term
            if (!first && term_coeff > 0) result << " + ";
            else if (!first && term_coeff < 0) {
                result << " - ";
                term_coeff = -term_coeff;
            }
            first = false;
            
            int x_power = power - k;
            
            if (x_power == 0) {
                // Just the constant
                result << term_coeff;
            } else {
                // Coefficient (if not 1)
                if (std::abs(term_coeff - 1.0) > 1e-10) {
                    result << term_coeff << "*";
                }
                
                // Variable
                result << var;
                if (x_power > 1) {
                    result << "^" << x_power;
                }
            }
        }
        
        return result.str();
    }
    
    // Pattern for simple (var)^power
    std::regex simple_power(R"(\(([a-zA-Z][0-9]*)\)\^([0-9]+))");
    if (std::regex_match(expr, match, simple_power)) {
        std::string var = match[1];
        int power = std::stoi(match[2]);
        
        if (power == 1) return var;
        return var + "^" + std::to_string(power);
    }
    
    // Pattern for (number)^power
    std::regex number_power(R"(\(([0-9]+(?:\.[0-9]+)?)\)\^([0-9]+))");
    if (std::regex_match(expr, match, number_power)) {
        double base = std::stod(match[1]);
        int power = std::stoi(match[2]);
        double result = std::pow(base, power);
        return std::to_string(result);
    }
    
    // If no patterns match, return original
    return expression;
}

std::string PolynomialExpander::polynomial_to_string(const Polynomial& poly) {
    std::ostringstream oss;
    
    // Sort monomials by total degree and lexicographic order
    auto sorted_poly = poly;
    std::sort(sorted_poly.begin(), sorted_poly.end(), 
        [](const Monomial& a, const Monomial& b) {
            // Calculate total degree
            int deg_a = 0, deg_b = 0;
            for (const auto& [var, pow] : a.variables) deg_a += pow;
            for (const auto& [var, pow] : b.variables) deg_b += pow;
            
            if (deg_a != deg_b) return deg_a > deg_b;
            
            // Same degree, compare lexicographically
            return a.variables < b.variables;
        });
    
    bool first = true;
    for (const auto& mono : sorted_poly) {
        if (std::abs(mono.coefficient) < 1e-10) continue;
        
        if (!first) {
            if (mono.coefficient > 0) oss << " + ";
            else oss << " - ";
        } else if (mono.coefficient < 0) {
            oss << "-";
        }
        first = false;
        
        double abs_coeff = std::abs(mono.coefficient);
        
        // Write coefficient if not 1 or if it's a constant term
        if (mono.variables.empty() || std::abs(abs_coeff - 1.0) > 1e-10) {
            oss << abs_coeff;
            if (!mono.variables.empty()) oss << "*";
        }
        
        // Write variables
        bool first_var = true;
        for (const auto& [var, power] : mono.variables) {
            if (!first_var) oss << "*";
            first_var = false;
            
            oss << var;
            if (power > 1) oss << "^" << power;
        }
    }
    
    return oss.str();
}

} // namespace julia_rur