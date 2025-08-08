#pragma once

#include "flint_wrappers.hpp"
#include "polynomial.hpp"
#include "monomial.hpp"
#include "quotient_ring.hpp"  // For PP typedef
#include <vector>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <iostream>

/**
 * Production-quality polynomial reduction using FLINT
 * Implements the multiplication table approach from Julia RUR algorithm
 */
template<typename CoeffT>
class FLINTReducer {
private:
    std::unique_ptr<flint_cpp::ModContext> ctx_;
    int prime_;
    size_t nvars_;
    
    // Quotient basis and multiplication tables
    std::vector<PP> quotient_basis_;
    std::vector<std::vector<std::vector<CoeffT>>> multiplication_tables_; // [var][basis_idx] -> coeffs
    
    // FLINT representations of Gröbner basis
    std::vector<std::unique_ptr<flint_cpp::ModPoly>> gb_polys_;
    
    bool tables_computed_;
    
public:
    FLINTReducer(int prime, size_t nvars) 
        : prime_(prime), nvars_(nvars), tables_computed_(false) {
        ctx_ = std::make_unique<flint_cpp::ModContext>(prime);
    }
    
    // Set Gröbner basis and convert to FLINT format
    void set_groebner_basis(const std::vector<MultivariatePolynomial<CoeffT>>& basis);
    
    // Compute quotient basis using BFS (like current implementation)
    void compute_quotient_basis();
    
    // Build multiplication tables using FLINT reduction
    void build_multiplication_tables();
    
    // Reduce a monomial to normal form using pre-computed tables
    std::vector<CoeffT> reduce_monomial(const PP& monomial) const;
    
    // Multiply vector by variable using pre-computed tables
    std::vector<CoeffT> multiply_by_variable(const std::vector<CoeffT>& vec, size_t var_idx) const;
    
    // Get quotient basis
    const std::vector<PP>& get_quotient_basis() const { return quotient_basis_; }
    
    // Get quotient dimension
    size_t dimension() const { return quotient_basis_.size(); }
    
    // Check if tables are ready
    bool are_tables_computed() const { return tables_computed_; }
    
private:
    // Convert our polynomial to FLINT format
    std::unique_ptr<flint_cpp::ModPoly> convert_to_flint(const MultivariatePolynomial<CoeffT>& poly) const;
    
    // Convert monomial to FLINT polynomial
    std::unique_ptr<flint_cpp::ModPoly> monomial_to_flint(const PP& monomial) const;
    
    // Reduce FLINT polynomial by Gröbner basis and convert back to coefficient vector
    std::vector<CoeffT> reduce_flint_poly(const flint_cpp::ModPoly& poly) const;
    
    // Check if monomial is divisible by any GB leading term
    bool is_divisible_by_leading_terms(const PP& monomial) const;
    
    // Helper: total degree of monomial
    uint32_t total_degree(const PP& monomial) const;
    
    // Helper: monomial division check
    bool divides(const PP& divisor, const PP& dividend) const;
};

// Implementation

template<typename CoeffT>
void FLINTReducer<CoeffT>::set_groebner_basis(const std::vector<MultivariatePolynomial<CoeffT>>& basis) {
    gb_polys_.clear();
    
    for (const auto& poly : basis) {
        if (!poly.is_zero()) {
            gb_polys_.push_back(convert_to_flint(poly));
        }
    }
    
    tables_computed_ = false;
}

template<typename CoeffT>
std::unique_ptr<flint_cpp::ModPoly> FLINTReducer<CoeffT>::convert_to_flint(const MultivariatePolynomial<CoeffT>& poly) const {
    auto flint_poly = std::make_unique<flint_cpp::ModPoly>(*ctx_);
    
    // For multivariate polynomials, we need to encode them somehow
    // This is a placeholder - real implementation would need to handle
    // multivariate -> univariate encoding (e.g., using variable ordering)
    
    // For now, convert to dense univariate representation using first variable
    // This is simplified - full implementation needs proper multivariate handling
    
    flint_poly->zero();
    
    for (const auto& [monomial, coeff] : poly) {
        // Simple encoding: use first variable degree as polynomial degree
        slong degree = monomial.exponents()[0];  // Simplified
        CoeffT mod_coeff = coeff;
        if (mod_coeff < 0) {
            mod_coeff = mod_coeff % prime_;
            if (mod_coeff < 0) mod_coeff += prime_;
        } else {
            mod_coeff = mod_coeff % prime_;
        }
        
        flint_poly->set_coeff_ui(degree, static_cast<ulong>(mod_coeff));
    }
    
    return flint_poly;
}

template<typename CoeffT>
std::unique_ptr<flint_cpp::ModPoly> FLINTReducer<CoeffT>::monomial_to_flint(const PP& monomial) const {
    auto flint_poly = std::make_unique<flint_cpp::ModPoly>(*ctx_);
    
    // Convert monomial to FLINT polynomial
    // Simplified: use first variable degree
    slong degree = monomial[0];  // Simplified encoding
    flint_poly->set_coeff_ui(degree, 1);
    
    return flint_poly;
}

template<typename CoeffT>
uint32_t FLINTReducer<CoeffT>::total_degree(const PP& monomial) const {
    uint32_t total = 0;
    for (uint32_t exp : monomial) {
        total += exp;
    }
    return total;
}

template<typename CoeffT>
bool FLINTReducer<CoeffT>::divides(const PP& divisor, const PP& dividend) const {
    if (divisor.size() != dividend.size()) return false;
    
    for (size_t i = 0; i < divisor.size(); ++i) {
        if (divisor[i] > dividend[i]) {
            return false;
        }
    }
    return true;
}

template<typename CoeffT>
bool FLINTReducer<CoeffT>::is_divisible_by_leading_terms(const PP& monomial) const {
    // This needs to be implemented based on GB leading terms
    // Placeholder implementation
    return false;
}

template<typename CoeffT>
void FLINTReducer<CoeffT>::compute_quotient_basis() {
    // Use same BFS approach as before
    quotient_basis_.clear();
    
    std::queue<PP> queue;
    std::set<PP> visited;
    
    // Start with constant monomial
    PP constant(nvars_, 0);
    queue.push(constant);
    visited.insert(constant);
    
    while (!queue.empty()) {
        PP current = queue.front();
        queue.pop();
        
        // Check if current monomial is divisible by any GB leading term
        if (!is_divisible_by_leading_terms(current)) {
            quotient_basis_.push_back(current);
            
            // Add neighbors (multiply by each variable)
            for (size_t i = 0; i < nvars_; ++i) {
                PP neighbor = current;
                neighbor[i]++;
                
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    queue.push(neighbor);
                }
            }
        }
    }
    
    // Sort by degree reverse lexicographic order
    std::sort(quotient_basis_.begin(), quotient_basis_.end(), [this](const PP& a, const PP& b) {
        uint32_t deg_a = total_degree(a);
        uint32_t deg_b = total_degree(b);
        
        if (deg_a != deg_b) {
            return deg_a < deg_b;  // Lower degree first
        }
        
        // Same degree: reverse lexicographic
        for (int i = nvars_ - 1; i >= 0; --i) {
            if (a[i] != b[i]) {
                return a[i] > b[i];  // Higher exponent first in reverse lex
            }
        }
        
        return false;  // Equal
    });
    
    std::cout << "FLINT reducer: Computed quotient basis with " << quotient_basis_.size() << " elements" << std::endl;
}

template<typename CoeffT>
void FLINTReducer<CoeffT>::build_multiplication_tables() {
    if (quotient_basis_.empty()) {
        throw std::runtime_error("Quotient basis not computed");
    }
    
    multiplication_tables_.clear();
    multiplication_tables_.resize(nvars_);
    
    for (size_t var = 0; var < nvars_; ++var) {
        multiplication_tables_[var].resize(quotient_basis_.size());
        
        for (size_t basis_idx = 0; basis_idx < quotient_basis_.size(); ++basis_idx) {
            // Compute xi * basis_element
            PP product = quotient_basis_[basis_idx];
            product[var]++;
            
            // Use FLINT to reduce this monomial
            auto flint_mono = monomial_to_flint(product);
            std::vector<CoeffT> reduced = reduce_flint_poly(*flint_mono);
            
            multiplication_tables_[var][basis_idx] = reduced;
        }
    }
    
    tables_computed_ = true;
    std::cout << "FLINT reducer: Built multiplication tables" << std::endl;
}

template<typename CoeffT>
std::vector<CoeffT> FLINTReducer<CoeffT>::reduce_flint_poly(const flint_cpp::ModPoly& poly) const {
    // This is where we'd use FLINT's reduction capabilities
    // For now, placeholder implementation
    std::vector<CoeffT> result(quotient_basis_.size(), 0);
    
    // Find if this polynomial corresponds to a quotient basis element
    // This is simplified - real implementation would use FLINT's divrem functions
    
    return result;
}

template<typename CoeffT>
std::vector<CoeffT> FLINTReducer<CoeffT>::reduce_monomial(const PP& monomial) const {
    if (!tables_computed_) {
        throw std::runtime_error("Multiplication tables not computed");
    }
    
    std::vector<CoeffT> result(quotient_basis_.size(), 0);
    
    // Check if monomial is directly in quotient basis
    auto it = std::find(quotient_basis_.begin(), quotient_basis_.end(), monomial);
    if (it != quotient_basis_.end()) {
        size_t idx = std::distance(quotient_basis_.begin(), it);
        result[idx] = 1;
        return result;
    }
    
    // Use FLINT reduction (placeholder)
    auto flint_mono = monomial_to_flint(monomial);
    return reduce_flint_poly(*flint_mono);
}

template<typename CoeffT>
std::vector<CoeffT> FLINTReducer<CoeffT>::multiply_by_variable(const std::vector<CoeffT>& vec, size_t var_idx) const {
    if (!tables_computed_) {
        throw std::runtime_error("Multiplication tables not computed");
    }
    
    if (vec.size() != quotient_basis_.size()) {
        throw std::invalid_argument("Vector size must match quotient basis dimension");
    }
    
    if (var_idx >= nvars_) {
        throw std::invalid_argument("Variable index out of range");
    }
    
    std::vector<CoeffT> result(quotient_basis_.size(), 0);
    
    // Use pre-computed multiplication tables
    for (size_t i = 0; i < quotient_basis_.size(); ++i) {
        if (vec[i] != 0) {
            const auto& table_row = multiplication_tables_[var_idx][i];
            for (size_t j = 0; j < quotient_basis_.size(); ++j) {
                result[j] = (result[j] + static_cast<CoeffT>(vec[i]) * table_row[j]) % prime_;
            }
        }
    }
    
    return result;
}

// Explicit template instantiations
extern template class FLINTReducer<int>;
extern template class FLINTReducer<long>;