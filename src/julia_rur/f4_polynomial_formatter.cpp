#include "f4_polynomial_formatter.hpp"
#include <string>
#include <vector>
#include <stdexcept>
#include <cctype>
#include <map>
#include <numeric>
#include <algorithm>
#include <gmpxx.h>
#include <sstream>

namespace julia_rur {

// Represents a monomial as a map from variable name to its power
using Monomial = std::map<std::string, int>;

// Represents a polynomial as a map from a monomial to its rational coefficient
using Polynomial = std::map<Monomial, mpq_class>;

// Forward declarations for the recursive parser
Polynomial parse_expression(std::string::const_iterator& it, const std::string::const_iterator& end);

// --- Helper Functions ---
long long gcd(long long a, long long b) { return b == 0 ? a : gcd(b, a % b); }
long long lcm(long long a, long long b) {
    if (a == 0 || b == 0) return 0;
    if (a == 1) return b;
    if (b == 1) return a;
    return std::abs(a * b) / gcd(a, b);
}

void trim_whitespace(std::string::const_iterator& it, const std::string::const_iterator& end) {
    while (it != end && std::isspace(*it)) {
        ++it;
    }
}

// --- Polynomial Arithmetic ---
Polynomial poly_add(const Polynomial& a, const Polynomial& b) {
    Polynomial result = a;
    for (const auto& [mon, coeff] : b) {
        result[mon] += coeff;
    }
    return result;
}

Polynomial poly_subtract(const Polynomial& a, const Polynomial& b) {
    Polynomial result = a;
    for (const auto& [mon, coeff] : b) {
        result[mon] -= coeff;
    }
    return result;
}

Polynomial poly_multiply(const Polynomial& a, const Polynomial& b) {
    Polynomial result;
    for (const auto& [mon_a, coeff_a] : a) {
        for (const auto& [mon_b, coeff_b] : b) {
            Monomial new_mon = mon_a;
            for (const auto& [var, exp] : mon_b) {
                new_mon[var] += exp;
            }
            result[new_mon] += coeff_a * coeff_b;
        }
    }
    return result;
}

Polynomial poly_pow(Polynomial base, int exp) {
    Polynomial result = {{{}, 1}}; // Start with identity polynomial "1"
    if (exp < 0) return {}; // Or throw error for negative exponent
    while (exp > 0) {
        if (exp % 2 == 1) result = poly_multiply(result, base);
        base = poly_multiply(base, base);
        exp /= 2;
    }
    return result;
}


// --- Recursive-Descent Parser ---

// factor ::= number | variable | '(' expression ')'
Polynomial parse_factor(std::string::const_iterator& it, const std::string::const_iterator& end) {
    trim_whitespace(it, end);
    if (it == end) throw std::runtime_error("Unexpected end of expression in parse_factor");

    if (*it == '(') {
        ++it; // Consume '('
        Polynomial result = parse_expression(it, end);
        trim_whitespace(it, end);
        if (it == end || *it != ')') throw std::runtime_error("Mismatched parentheses");
        ++it; // Consume ')'
        return result;
    }

    if (std::isdigit(*it) || *it == '.') {
        std::string num_str;
        while (it != end && (std::isdigit(*it) || *it == '.')) {
            num_str += *it;
            ++it;
        }
        return {{{}, mpq_class(num_str)}};
    }
    
    if (std::isalpha(*it)) {
        std::string var_name;
        while (it != end && std::isalnum(*it)) {
            var_name += *it;
            ++it;
        }
        return {{{{var_name, 1}}, 1}}; // Represents "1 * var_name^1"
    }

    throw std::runtime_error("Unexpected character in expression");
}

// term ::= factor ('^' factor)*
Polynomial parse_power(std::string::const_iterator& it, const std::string::const_iterator& end) {
    Polynomial result = parse_factor(it, end);
    trim_whitespace(it, end);
    while (it != end && *it == '^') {
        ++it; // Consume '^'
        Polynomial exponent_poly = parse_factor(it, end);
        // Exponent must be a simple integer constant
        if (exponent_poly.size() != 1 || !exponent_poly.begin()->first.empty() || exponent_poly.begin()->second.get_den() != 1) {
            throw std::runtime_error("Exponent must be an integer");
        }
        int exp = exponent_poly.begin()->second.get_num().get_si();
        result = poly_pow(result, exp);
        trim_whitespace(it, end);
    }
    return result;
}

// product ::= power ('*' power)*
Polynomial parse_product(std::string::const_iterator& it, const std::string::const_iterator& end) {
    Polynomial result = parse_power(it, end);
    trim_whitespace(it, end);
    while (it != end && *it == '*') {
        ++it; // Consume '*'
        result = poly_multiply(result, parse_power(it, end));
        trim_whitespace(it, end);
    }
    return result;
}

// expression ::= product (('+' | '-') product)*
Polynomial parse_expression(std::string::const_iterator& it, const std::string::const_iterator& end) {
    trim_whitespace(it, end);
    
    // Handle leading sign
    bool negative_first = false;
    if (it != end && *it == '-') {
        negative_first = true;
        ++it;
    } else if (it != end && *it == '+') {
        ++it;
    }
    
    Polynomial result = parse_product(it, end);
    if (negative_first) {
        // Negate the first term
        Polynomial negated;
        for (const auto& [mon, coeff] : result) {
            negated[mon] = -coeff;
        }
        result = negated;
    }
    
    trim_whitespace(it, end);
    while (it != end && (*it == '+' || *it == '-')) {
        char op = *it;
        ++it;
        if (op == '+') {
            result = poly_add(result, parse_product(it, end));
        } else {
            result = poly_subtract(result, parse_product(it, end));
        }
        trim_whitespace(it, end);
    }
    return result;
}

// Main formatting function
std::string format_polynomial_for_f4(const std::string& polynomial_str) {
    auto it = polynomial_str.begin();
    Polynomial p = parse_expression(it, polynomial_str.end());

    // Find LCM of all denominators
    long long global_lcm = 1;
    for (const auto& [mon, coeff] : p) {
        if (coeff.get_den() != 1) {
            global_lcm = lcm(global_lcm, coeff.get_den().get_si());
        }
    }

    // Rebuild string with integer coefficients
    std::stringstream result;
    bool first_term = true;
    for (const auto& [mon, coeff] : p) {
        mpq_class scaled = coeff * static_cast<long>(global_lcm);
        mpz_class num = scaled.get_num();
        if (num == 0) continue;

        if (num > 0 && !first_term) {
            result << "+";
        }
        
        bool is_constant = mon.empty() || std::all_of(mon.begin(), mon.end(), [](auto const& p){ return p.second == 0; });
        
        if (num == 1 && !is_constant) {
            result << "1*";
        } else if (num == -1 && !is_constant) {
            result << "-1*";
        } else {
            result << num.get_str();
            if(!is_constant) result << "*";
        }

        bool first_var = true;
        for (const auto& [var, exp] : mon) {
            if (exp > 0) {
                if (!first_var) result << "*";
                result << var;
                if (exp > 1) result << "^" << exp;
                first_var = false;
            }
        }
        first_term = false;
    }
    
    std::string final_str = result.str();
    if (final_str.empty()) return "0";
    if (final_str[0] == '+') return final_str.substr(1);
    return final_str;
}

}  // namespace julia_rur