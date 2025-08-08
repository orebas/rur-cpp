#pragma once

#include "monomial.hpp"
#include "polynomial.hpp"
#include <vector>
#include <cstdint>
#include <map>
#include <memory>

// Type aliases matching Julia implementation
using PP = std::vector<uint32_t>;        // Power product (monomial exponent vector)
using Deg = uint32_t;                    // Degree type
using ModularCoeff = uint32_t;           // Modular coefficient
using AccumCoeff = uint64_t;             // Accumulator for modular arithmetic

/**
 * StackVect structure for exploring the border of the quotient ring
 * Corresponds to Julia's StackVect (Lines 57-62 in JULIA_RUR_ALGORITHM.md)
 */
struct StackVect {
    int32_t pos;        // Position index
    PP mon;             // Monomial (power product)
    int32_t prev;       // Previous element in chain
    int32_t var;        // Variable index
    
    StackVect() : pos(-1), prev(-1), var(-1) {}
    StackVect(int32_t p, const PP& m, int32_t pr, int32_t v) 
        : pos(p), mon(m), prev(pr), var(v) {}
};

/**
 * Quotient ring representation with multiplication tables
 * Implements the core data structures from Julia RUR algorithm
 */
template<typename CoeffT>
class QuotientRing {
private:
    int prime_;
    size_t nvars_;
    
    // Core data arrays (matching Julia implementation)
    std::vector<PP> ltg_;                                    // Leading terms of Gröbner basis
    std::vector<PP> quo_;                                    // Quotient basis monomials (sorted by DRL)
    std::vector<StackVect> t_xw_;                           // Border structure
    std::vector<std::vector<int32_t>> i_xw_;                // Variable indices: i_xw[var][basis_idx] = index
    std::vector<std::vector<ModularCoeff>> t_v_;            // Coefficient vectors
    
    // Gröbner basis polynomials
    std::vector<MultivariatePolynomial<CoeffT>> groebner_basis_;
    
    bool tables_computed_;
    
public:
    QuotientRing(int prime, size_t nvars) 
        : prime_(prime), nvars_(nvars), tables_computed_(false) {
        i_xw_.resize(nvars);
    }
    
    // Initialize with Gröbner basis
    void set_groebner_basis(const std::vector<MultivariatePolynomial<CoeffT>>& basis);
    
    // Compute quotient basis using BFS (Julia: compute_quotient_basis, Lines 181-214)
    void compute_quotient_basis();
    
    // Build multiplication table structure (Julia: prepare_table_mxi, Lines 341-392)
    void prepare_multiplication_tables();
    
    // Fill multiplication tables (Julia: learn_compute_table!, Lines 439-463)
    void compute_multiplication_tables();
    
    // Compute normal form of a monomial in quotient ring
    std::vector<CoeffT> reduce_monomial(const PP& monomial) const;
    
    // Multiply vector by variable xi (Julia: _mul_var_quo!, Lines 398-431)
    std::vector<CoeffT> multiply_by_variable(const std::vector<CoeffT>& vec, size_t var_idx) const;
    
    // Check if monomial is divisible by any GB leading term
    bool is_divisible_by_leading_terms(const PP& monomial) const;
    
    // Get quotient basis
    const std::vector<PP>& get_quotient_basis() const { return quo_; }
    
    // Get basis size
    size_t dimension() const { return quo_.size(); }
    
    // Get variable count
    size_t num_variables() const { return nvars_; }
    
    // Check if tables are computed
    bool are_tables_computed() const { return tables_computed_; }
    
private:
    // Helper: compare monomials in DRL (degree reverse lexicographic) order
    bool drl_compare(const PP& a, const PP& b) const;
    
    // Helper: compute total degree of monomial
    uint32_t total_degree(const PP& monomial) const;
    
    // Helper: monomial division
    bool divides(const PP& divisor, const PP& dividend) const;
    
    // Helper: monomial multiplication
    PP multiply_monomials(const PP& a, const PP& b) const;
};

// Explicit template instantiation declarations
extern template class QuotientRing<int>;
extern template class QuotientRing<long>;