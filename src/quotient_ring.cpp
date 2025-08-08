#include "quotient_ring.hpp"
#include <algorithm>
#include <queue>
#include <set>
#include <numeric>
#include <stdexcept>
#include <iostream>

template<typename CoeffT>
void QuotientRing<CoeffT>::set_groebner_basis(const std::vector<MultivariatePolynomial<CoeffT>>& basis) {
    groebner_basis_ = basis;
    
    // Extract leading terms
    ltg_.clear();
    for (const auto& poly : basis) {
        if (!poly.is_zero()) {
            auto exps = poly.leading_monomial().exponents();
            PP leading_term(exps.begin(), exps.end());  // Convert int to uint32_t
            ltg_.push_back(leading_term);
        }
    }
    
    tables_computed_ = false;
}

template<typename CoeffT>
uint32_t QuotientRing<CoeffT>::total_degree(const PP& monomial) const {
    return std::accumulate(monomial.begin(), monomial.end(), 0u);
}

template<typename CoeffT>
bool QuotientRing<CoeffT>::drl_compare(const PP& a, const PP& b) const {
    // Degree reverse lexicographic order
    uint32_t deg_a = total_degree(a);
    uint32_t deg_b = total_degree(b);
    
    if (deg_a != deg_b) {
        return deg_a < deg_b;  // Lower degree first
    }
    
    // Same degree: reverse lexicographic (rightmost variable first)
    for (int i = nvars_ - 1; i >= 0; --i) {
        if (a[i] != b[i]) {
            return a[i] > b[i];  // Higher exponent first in reverse lex
        }
    }
    
    return false;  // Equal
}

template<typename CoeffT>
bool QuotientRing<CoeffT>::divides(const PP& divisor, const PP& dividend) const {
    if (divisor.size() != dividend.size()) return false;
    
    for (size_t i = 0; i < divisor.size(); ++i) {
        if (divisor[i] > dividend[i]) {
            return false;
        }
    }
    return true;
}

template<typename CoeffT>
PP QuotientRing<CoeffT>::multiply_monomials(const PP& a, const PP& b) const {
    PP result(nvars_);
    for (size_t i = 0; i < nvars_; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template<typename CoeffT>
bool QuotientRing<CoeffT>::is_divisible_by_leading_terms(const PP& monomial) const {
    for (const auto& lt : ltg_) {
        if (divides(lt, monomial)) {
            return true;
        }
    }
    return false;
}

template<typename CoeffT>
void QuotientRing<CoeffT>::compute_quotient_basis() {
    // BFS traversal to find quotient basis (Julia: compute_quotient_basis, Lines 181-214)
    quo_.clear();
    
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
            quo_.push_back(current);
            
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
    
    // Sort by DRL order
    std::sort(quo_.begin(), quo_.end(), [this](const PP& a, const PP& b) {
        return drl_compare(a, b);
    });
    
    std::cout << "Computed quotient basis with " << quo_.size() << " elements" << std::endl;
}

template<typename CoeffT>
void QuotientRing<CoeffT>::prepare_multiplication_tables() {
    if (quo_.empty()) {
        throw std::runtime_error("Quotient basis not computed");
    }
    
    // Initialize structures (Julia: prepare_table_mxi, Lines 341-392)
    t_xw_.clear();
    for (auto& var_indices : i_xw_) {
        var_indices.clear();
        var_indices.resize(quo_.size(), -1);
    }
    
    // For each basis monomial and each variable
    for (size_t var = 0; var < nvars_; ++var) {
        for (size_t basis_idx = 0; basis_idx < quo_.size(); ++basis_idx) {
            PP product = quo_[basis_idx];
            product[var]++;
            
            // Check if product is in quotient basis
            auto it = std::find(quo_.begin(), quo_.end(), product);
            if (it != quo_.end()) {
                // Product is in quotient basis
                i_xw_[var][basis_idx] = std::distance(quo_.begin(), it);
            } else if (!is_divisible_by_leading_terms(product)) {
                // Product is a new border element
                int32_t border_idx = t_xw_.size();
                t_xw_.emplace_back(border_idx, product, -1, var);
                i_xw_[var][basis_idx] = -(border_idx + 1);  // Negative index for border elements
            } else {
                // Product reduces to zero (divisible by GB leading term)
                i_xw_[var][basis_idx] = -1;
            }
        }
    }
    
    std::cout << "Prepared multiplication tables with " << t_xw_.size() << " border elements" << std::endl;
}

template<typename CoeffT>
void QuotientRing<CoeffT>::compute_multiplication_tables() {
    // For simple cases where all multiplications stay in quotient basis,
    // we may have no border elements (t_xw_ is empty)
    // In this case, the multiplication tables are already complete
    
    if (t_xw_.empty()) {
        std::cout << "No border elements - all multiplications stay in quotient basis" << std::endl;
        tables_computed_ = true;
        return;
    }
    
    // Initialize coefficient vectors (Julia: compute_fill_quo_gb!, Lines 464-484)
    t_v_.clear();
    t_v_.resize(t_xw_.size(), std::vector<ModularCoeff>(quo_.size(), 0));
    
    // For each border element, compute its expansion in quotient basis
    // This is a simplified version - in practice, we need to reduce using GB
    for (size_t border_idx = 0; border_idx < t_xw_.size(); ++border_idx) {
        const auto& border_elem = t_xw_[border_idx];
        
        // Find the monomial that reduces to this border element
        PP monomial = border_elem.mon;
        
        // Reduce monomial using Gröbner basis (simplified)
        std::vector<CoeffT> reduced = reduce_monomial(monomial);
        
        // Convert to modular coefficients
        for (size_t i = 0; i < reduced.size() && i < quo_.size(); ++i) {
            CoeffT coeff = reduced[i];
            if (coeff < 0) {
                coeff = coeff % prime_;
                if (coeff < 0) coeff += prime_;
            } else {
                coeff = coeff % prime_;
            }
            t_v_[border_idx][i] = static_cast<ModularCoeff>(coeff);
        }
    }
    
    tables_computed_ = true;
    std::cout << "Computed multiplication tables" << std::endl;
}

template<typename CoeffT>
std::vector<CoeffT> QuotientRing<CoeffT>::reduce_monomial(const PP& monomial) const {
    // This is a placeholder implementation
    // In the full implementation, we would reduce the monomial using the Gröbner basis
    // For now, return zero vector except for monomials in the quotient basis
    
    std::vector<CoeffT> result(quo_.size(), 0);
    
    auto it = std::find(quo_.begin(), quo_.end(), monomial);
    if (it != quo_.end()) {
        size_t idx = std::distance(quo_.begin(), it);
        result[idx] = 1;
    }
    
    return result;
}

template<typename CoeffT>
std::vector<CoeffT> QuotientRing<CoeffT>::multiply_by_variable(const std::vector<CoeffT>& vec, size_t var_idx) const {
    if (!tables_computed_) {
        throw std::runtime_error("Multiplication tables not computed");
    }
    
    if (vec.size() != quo_.size()) {
        throw std::invalid_argument("Vector size must match quotient basis dimension");
    }
    
    if (var_idx >= nvars_) {
        throw std::invalid_argument("Variable index out of range");
    }
    
    std::vector<CoeffT> result(quo_.size(), 0);
    
    // Multiply each basis element by the variable (Julia: _mul_var_quo!, Lines 398-431)
    for (size_t i = 0; i < quo_.size(); ++i) {
        if (vec[i] == 0) continue;
        
        int32_t target_idx = i_xw_[var_idx][i];
        
        if (target_idx >= 0) {
            // Maps to quotient basis element
            result[target_idx] = (result[target_idx] + vec[i]) % prime_;
        } else if (target_idx < -1) {
            // Maps to border element
            int32_t border_idx = -(target_idx + 1);
            
            // Add contribution from border element's expansion
            for (size_t j = 0; j < quo_.size(); ++j) {
                CoeffT contribution = (static_cast<CoeffT>(vec[i]) * static_cast<CoeffT>(t_v_[border_idx][j])) % prime_;
                result[j] = (result[j] + contribution) % prime_;
            }
        }
        // target_idx == -1 means multiplication gives zero (no contribution)
    }
    
    return result;
}

// Explicit template instantiations
template class QuotientRing<int>;
template class QuotientRing<long>;