#pragma once

#include <cstdint>
#include <iostream>
#include <vector>


/**
 * Julia-style RUR data structures
 * Mirrors RationalUnivariateRepresentation.jl exactly
 */

namespace julia_rur {

enum class SeparatingStrategy { CURRENT, RANDOM, DETERMINISTIC };

// Type definitions matching Julia exactly
using PP = std::vector<uint32_t>; // Power product (Julia: Vector{Deg})
using ModularCoeff = uint32_t;    // Modular coefficient (Julia: UInt32)
using AccModularCoeff = uint64_t; // Accumulator (Julia: UInt64)
using Deg = uint32_t;             // Degree type (Julia: UInt32)

/**
 * StackVect structure - exact mirror of Julia implementation
 * Used to explore the border of the quotient ring during multiplication table construction
 */
struct StackVect {
    int32_t pos;  // Position index
    PP mon;       // Monomial (power product)
    int32_t prev; // Previous element in chain
    int32_t var;  // Variable index

    // Constructor with parameters (no default constructor for safety)
    StackVect(int32_t position, const PP &monomial, int32_t previous, int32_t variable)
      : pos(position)
      , mon(monomial)
      , prev(previous)
      , var(variable) {}

    // Copy constructor (default is fine)
    StackVect(const StackVect &) = default;
    StackVect &operator=(const StackVect &) = default;

    // Move constructor
    StackVect(StackVect &&) = default;
    StackVect &operator=(StackVect &&) = default;

    // Equality comparison for testing
    bool operator==(const StackVect &other) const {
        return pos == other.pos && mon == other.mon && prev == other.prev && var == other.var;
    }

    bool operator!=(const StackVect &other) const { return !(*this == other); }

    // Debug output
    void print(std::ostream &os = std::cout) const {
        os << "StackVect{pos=" << pos << ", mon=[";
        for (size_t i = 0; i < mon.size(); ++i) {
            if (i > 0) os << ",";
            os << mon[i];
        }
        os << "], prev=" << prev << ", var=" << var << "}";
    }
};

/**
 * Core data arrays matching Julia variable names exactly
 * This represents the quotient ring with multiplication tables
 */
class QuotientRingData {
  public:
    // Core arrays - matching Julia names exactly
    std::vector<PP> quo;                        // Quotient basis monomials (Julia: quo)
    std::vector<StackVect> t_xw;                // Border structure (Julia: t_xw)
    std::vector<std::vector<int32_t>> i_xw;     // Variable multiplication indices (Julia: i_xw)
    std::vector<std::vector<ModularCoeff>> t_v; // Coefficient vectors (Julia: t_v)

    // Metadata
    size_t nvars;       // Number of variables
    ModularCoeff prime; // Modular arithmetic prime

    // Constructor
    explicit QuotientRingData(size_t num_variables, ModularCoeff modular_prime)
      : nvars(num_variables)
      , prime(modular_prime) {
        i_xw.resize(nvars); // One vector per variable (others default-constructed empty)
    }

    // Size accessors
    size_t quotient_basis_size() const { return quo.size(); }
    size_t border_size() const { return t_xw.size(); }
    size_t total_elements() const { return quotient_basis_size() + border_size(); }

    // Validation
    bool is_valid() const {
        // Check basic consistency
        if (i_xw.size() != nvars) return false;
        if (quo.empty()) return false;

        // Check that all monomials have correct dimension
        for (const auto &monomial : quo) {
            if (monomial.size() != nvars) return false;
        }

        for (const auto &stack_elem : t_xw) {
            if (stack_elem.mon.size() != nvars) return false;
        }

        return true;
    }

    // Debug output
    void print_summary(std::ostream &os = std::cout) const {
        os << "QuotientRingData Summary:\n";
        os << "  nvars: " << nvars << "\n";
        os << "  prime: " << prime << "\n";
        os << "  quotient_basis_size: " << quotient_basis_size() << "\n";
        os << "  border_size: " << border_size() << "\n";
        os << "  total_elements: " << total_elements() << "\n";

        if (!quo.empty()) {
            os << "  Quotient basis (first few):\n";
            for (size_t i = 0; i < std::min(quo.size(), size_t(5)); ++i) {
                os << "    " << i << ": [";
                for (size_t j = 0; j < quo[i].size(); ++j) {
                    if (j > 0) os << ",";
                    os << quo[i][j];
                }
                os << "]\n";
            }
        }

        if (!t_xw.empty()) {
            os << "  Border elements (first few):\n";
            for (size_t i = 0; i < std::min(t_xw.size(), size_t(3)); ++i) {
                os << "    " << i << ": ";
                t_xw[i].print(os);
                os << "\n";
            }
        }
    }

    // Clear all data
    void clear() {
        quo.clear();
        t_xw.clear();
        i_xw = std::vector<std::vector<int32_t>>(nvars); // More idiomatic
        t_v.clear();
    }

    // Reserve space for efficiency
    void reserve_quotient_basis(size_t size) { quo.reserve(size); }

    void reserve_border(size_t size) {
        t_xw.reserve(size);
        t_v.reserve(size);
    }
};

/**
 * Helper functions for power product (monomial) operations
 * Matching Julia's monomial arithmetic
 */
namespace power_product {

// Create monomial with all zeros
inline PP
zero_monomial(size_t nvars) {
    return PP(nvars, 0);
}

// Create unit monomial for variable i (e_i)
inline PP
unit_monomial(size_t nvars, size_t var_index) {
    PP result(nvars, 0);
    if (var_index < nvars) { result[var_index] = 1; }
    return result;
}

// Multiply two monomials: result[i] = a[i] + b[i]
inline PP
multiply(const PP &a, const PP &b) {
    if (a.size() != b.size()) { throw std::invalid_argument("Monomial dimensions must match"); }

    PP result(a.size());
    for (size_t i = 0; i < a.size(); ++i) { result[i] = a[i] + b[i]; }
    return result;
}

// Check if monomial a divides monomial b
inline bool
divides(const PP &a, const PP &b) {
    if (a.size() != b.size()) { throw std::invalid_argument("Monomial dimensions must match for division check"); }

    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) return false;
    }
    return true;
}

// Compute total degree of monomial
inline uint32_t
total_degree(const PP &monomial) {
    uint32_t deg = 0;
    for (uint32_t exp : monomial) { deg += exp; }
    return deg;
}

// Compare monomials for degree reverse lexicographic order
inline int
compare_degrevlex(const PP &a, const PP &b) {
    if (a.size() != b.size()) { throw std::invalid_argument("Monomial dimensions must match"); }

    // First compare total degrees
    uint32_t deg_a = total_degree(a);
    uint32_t deg_b = total_degree(b);

    if (deg_a < deg_b) return -1;
    if (deg_a > deg_b) return 1;

    // Same degree - reverse lexicographic comparison
    // For degrevlex: compare from right to left (reverse of lexicographic)
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
        if (a[i] < b[i]) return 1;  // SMALLER exponent in rightmost variables comes FIRST
        if (a[i] > b[i]) return -1; // (this is the "reverse" part)
    }

    return 0; // Equal
}

// Print monomial for debugging
inline void
print_monomial(const PP &monomial, std::ostream &os = std::cout) {
    os << "[";
    for (size_t i = 0; i < monomial.size(); ++i) {
        if (i > 0) os << ",";
        os << monomial[i];
    }
    os << "]";
}
}

} // namespace julia_rur