#pragma once

#include "data_structures.hpp"
#include <vector>
#include <stdexcept>

namespace julia_rur {

/**
 * Polynomial operations modulo a prime p
 * Polynomials are represented as coefficient vectors (dense representation)
 * poly[i] is the coefficient of t^i
 */

// Compute modular inverse of a modulo p using Fermat's little theorem
inline ModularCoeff modular_inverse(ModularCoeff a, ModularCoeff p) {
    ModularCoeff result = 1;
    ModularCoeff base = a;
    ModularCoeff exp = p - 2;
    while (exp > 0) {
        if (exp & 1) {
            result = (static_cast<AccModularCoeff>(result) * base) % p;
        }
        base = (static_cast<AccModularCoeff>(base) * base) % p;
        exp >>= 1;
    }
    return result;
}

// Remove leading zeros from polynomial
inline void normalize_polynomial(std::vector<ModularCoeff>& poly) {
    while (poly.size() > 1 && poly.back() == 0) {
        poly.pop_back();
    }
}

// Compute derivative of polynomial
inline std::vector<ModularCoeff> polynomial_derivative(
    const std::vector<ModularCoeff>& poly,
    ModularCoeff prime
) {
    if (poly.size() <= 1) {
        return {0};  // Derivative of constant is 0
    }
    
    std::vector<ModularCoeff> deriv(poly.size() - 1);
    for (size_t i = 1; i < poly.size(); ++i) {
        // deriv[i-1] = i * poly[i] mod p
        deriv[i-1] = static_cast<ModularCoeff>(
            (static_cast<AccModularCoeff>(i) * poly[i]) % prime
        );
    }
    
    normalize_polynomial(deriv);
    return deriv;
}

// Add two polynomials modulo prime
inline std::vector<ModularCoeff> polynomial_add(
    const std::vector<ModularCoeff>& a,
    const std::vector<ModularCoeff>& b,
    ModularCoeff prime
) {
    size_t max_size = std::max(a.size(), b.size());
    std::vector<ModularCoeff> result(max_size, 0);
    
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i];
    }
    
    for (size_t i = 0; i < b.size(); ++i) {
        result[i] = (result[i] + b[i]) % prime;
    }
    
    normalize_polynomial(result);
    return result;
}

// Subtract two polynomials modulo prime
inline std::vector<ModularCoeff> polynomial_subtract(
    const std::vector<ModularCoeff>& a,
    const std::vector<ModularCoeff>& b,
    ModularCoeff prime
) {
    size_t max_size = std::max(a.size(), b.size());
    std::vector<ModularCoeff> result(max_size, 0);
    
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i];
    }
    
    for (size_t i = 0; i < b.size(); ++i) {
        if (result[i] >= b[i]) {
            result[i] = (result[i] - b[i]) % prime;
        } else {
            result[i] = (result[i] + prime - b[i]) % prime;
        }
    }
    
    normalize_polynomial(result);
    return result;
}

// Multiply two polynomials modulo prime (without reduction by minimal polynomial)
inline std::vector<ModularCoeff> polynomial_multiply(
    const std::vector<ModularCoeff>& a,
    const std::vector<ModularCoeff>& b,
    ModularCoeff prime
) {
    if (a.empty() || b.empty()) {
        return {0};
    }
    
    std::vector<ModularCoeff> result(a.size() + b.size() - 1, 0);
    
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            AccModularCoeff prod = static_cast<AccModularCoeff>(a[i]) * b[j];
            result[i + j] = (result[i + j] + prod) % prime;
        }
    }
    
    normalize_polynomial(result);
    return result;
}

// Divide polynomial a by polynomial b, returning (quotient, remainder)
inline std::pair<std::vector<ModularCoeff>, std::vector<ModularCoeff>> polynomial_divmod(
    std::vector<ModularCoeff> a,  // Copy since we modify it
    const std::vector<ModularCoeff>& b,
    ModularCoeff prime
) {
    normalize_polynomial(a);
    
    if (b.empty() || (b.size() == 1 && b[0] == 0)) {
        throw std::runtime_error("Division by zero polynomial");
    }
    
    // Handle constant divisor
    if (b.size() == 1) {
        ModularCoeff b_inv = modular_inverse(b[0], prime);
        std::vector<ModularCoeff> quotient(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            quotient[i] = static_cast<ModularCoeff>(
                (static_cast<AccModularCoeff>(a[i]) * b_inv) % prime
            );
        }
        return {quotient, {0}};
    }
    
    std::vector<ModularCoeff> quotient;
    
    // Get leading coefficient of b and its inverse
    ModularCoeff b_lead = b.back();
    ModularCoeff b_lead_inv = modular_inverse(b_lead, prime);
    
    while (a.size() >= b.size() && !(a.size() == 1 && a[0] == 0)) {
        // Compute leading coefficient of quotient term
        ModularCoeff q_coeff = static_cast<ModularCoeff>(
            (static_cast<AccModularCoeff>(a.back()) * b_lead_inv) % prime
        );
        quotient.push_back(q_coeff);
        
        // Subtract q_coeff * t^(deg(a)-deg(b)) * b from a
        size_t degree_diff = a.size() - b.size();
        for (size_t i = 0; i < b.size(); ++i) {
            AccModularCoeff sub = static_cast<AccModularCoeff>(q_coeff) * b[i];
            sub %= prime;
            
            if (a[i + degree_diff] >= sub) {
                a[i + degree_diff] = (a[i + degree_diff] - sub) % prime;
            } else {
                a[i + degree_diff] = (a[i + degree_diff] + prime - sub) % prime;
            }
        }
        
        // Remove leading zero
        a.pop_back();
    }
    
    // Reverse quotient (we built it backwards)
    std::reverse(quotient.begin(), quotient.end());
    
    if (quotient.empty()) {
        quotient = {0};
    }
    
    normalize_polynomial(a);  // a is now the remainder
    return {quotient, a};
}

// Reduce polynomial modulo another polynomial
inline std::vector<ModularCoeff> polynomial_mod(
    const std::vector<ModularCoeff>& a,
    const std::vector<ModularCoeff>& b,
    ModularCoeff prime
) {
    return polynomial_divmod(a, b, prime).second;
}

// Extended Euclidean algorithm for polynomials
// Returns (gcd, x, y) such that a*x + b*y = gcd
struct PolynomialGCDResult {
    std::vector<ModularCoeff> gcd;
    std::vector<ModularCoeff> x;
    std::vector<ModularCoeff> y;
};

inline PolynomialGCDResult polynomial_extended_gcd(
    const std::vector<ModularCoeff>& a,
    const std::vector<ModularCoeff>& b,
    ModularCoeff prime
) {
    // Initialize
    std::vector<ModularCoeff> old_r = a;
    std::vector<ModularCoeff> r = b;
    std::vector<ModularCoeff> old_s = {1};
    std::vector<ModularCoeff> s = {0};
    std::vector<ModularCoeff> old_t = {0};
    std::vector<ModularCoeff> t = {1};
    
    while (r.size() > 1 || (r.size() == 1 && r[0] != 0)) {
        auto [quotient, remainder] = polynomial_divmod(old_r, r, prime);
        
        old_r = r;
        r = remainder;
        
        // Update s
        auto temp = s;
        s = polynomial_subtract(old_s, polynomial_multiply(quotient, s, prime), prime);
        old_s = temp;
        
        // Update t
        temp = t;
        t = polynomial_subtract(old_t, polynomial_multiply(quotient, t, prime), prime);
        old_t = temp;
    }
    
    // Normalize gcd to be monic
    if (!old_r.empty() && old_r.back() != 0 && old_r.back() != 1) {
        ModularCoeff lead_inv = modular_inverse(old_r.back(), prime);
        for (auto& coeff : old_r) {
            coeff = static_cast<ModularCoeff>(
                (static_cast<AccModularCoeff>(coeff) * lead_inv) % prime
            );
        }
        for (auto& coeff : old_s) {
            coeff = static_cast<ModularCoeff>(
                (static_cast<AccModularCoeff>(coeff) * lead_inv) % prime
            );
        }
        for (auto& coeff : old_t) {
            coeff = static_cast<ModularCoeff>(
                (static_cast<AccModularCoeff>(coeff) * lead_inv) % prime
            );
        }
    }
    
    return {old_r, old_s, old_t};
}

// Compute modular inverse of polynomial a modulo polynomial m
inline std::vector<ModularCoeff> polynomial_inverse_mod(
    const std::vector<ModularCoeff>& a,
    const std::vector<ModularCoeff>& m,
    ModularCoeff prime
) {
    auto [gcd, x, y] = polynomial_extended_gcd(a, m, prime);
    
    // Check if gcd is 1
    if (gcd.size() != 1 || gcd[0] != 1) {
        throw std::runtime_error("Polynomial is not invertible modulo the given polynomial");
    }
    
    // x is the inverse of a modulo m
    return x;
}

// Multiply polynomial by scalar modulo prime
inline std::vector<ModularCoeff> polynomial_scalar_multiply(
    const std::vector<ModularCoeff>& poly,
    ModularCoeff scalar,
    ModularCoeff prime
) {
    std::vector<ModularCoeff> result(poly.size());
    for (size_t i = 0; i < poly.size(); ++i) {
        result[i] = static_cast<ModularCoeff>(
            (static_cast<AccModularCoeff>(poly[i]) * scalar) % prime
        );
    }
    normalize_polynomial(result);
    return result;
}

// Compute the RUR numerator from biv_lex output
// N(t) = (f(t) * h(t)^(-1) * T'(t)) mod T(t)
inline std::vector<ModularCoeff> compute_rur_numerator(
    const std::vector<ModularCoeff>& f_t,  // f(t) from relation h(t)*x - f(t) = 0
    const std::vector<ModularCoeff>& h_t,  // h(t) from relation
    const std::vector<ModularCoeff>& minimal_poly,  // T(t)
    ModularCoeff prime
) {
    // Compute T'(t)
    auto t_prime = polynomial_derivative(minimal_poly, prime);
    
    // Compute h(t)^(-1) mod T(t)
    auto h_inv = polynomial_inverse_mod(h_t, minimal_poly, prime);
    
    // Compute f(t) * h(t)^(-1) mod T(t)
    auto f_h_inv = polynomial_multiply(f_t, h_inv, prime);
    f_h_inv = polynomial_mod(f_h_inv, minimal_poly, prime);
    
    // Compute (f(t) * h(t)^(-1)) * T'(t) mod T(t)
    auto result = polynomial_multiply(f_h_inv, t_prime, prime);
    result = polynomial_mod(result, minimal_poly, prime);
    
    return result;
}

} // namespace julia_rur