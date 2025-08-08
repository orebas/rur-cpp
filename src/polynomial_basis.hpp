#pragma once

#include "flint_mpoly_wrappers.hpp" 
#include "polynomial.hpp"
#include "monomial.hpp"
#include <vector>
#include <unordered_map>
#include <memory>
#include <queue>
#include <set>

/**
 * Immutable representation of pre-computed polynomial basis information
 * Separates the expensive basis computation from the fast reduction operations
 */
template<typename CoeffT>
struct PolynomialBasis {
    // FLINT representation of Gröbner basis  
    std::vector<flint_mpoly::Mpoly> grobner_basis_flint;
    
    // Quotient basis (standard monomials)
    std::vector<Monomial> quotient_basis_monomials;
    
    // Fast lookup: monomial -> index in quotient basis
    std::unordered_map<Monomial, size_t> monomial_to_index;
    
    // Metadata
    size_t nvars;
    ulong prime;
    ordering_t ordering;
    
    /**
     * Factory method to create polynomial basis from Gröbner basis
     */
    static std::shared_ptr<const PolynomialBasis> create(
        const std::vector<MultivariatePolynomial<CoeffT>>& groebner_basis,
        size_t nvars,
        ulong prime,
        ordering_t ordering = ORD_DEGREVLEX
    );
    
private:
    // Private constructor - use factory method
    PolynomialBasis() = default;
    
    // Convert our polynomial to FLINT format with proper modular reduction
    static flint_mpoly::Mpoly to_flint_mpoly(
        const MultivariatePolynomial<CoeffT>& poly,
        std::shared_ptr<const flint_mpoly::MpolyContext> ctx,
        ulong prime
    );
    
    // Compute quotient basis using BFS on monomials not divisible by GB leading terms
    void compute_quotient_basis(
        const std::vector<MultivariatePolynomial<CoeffT>>& groebner_basis,
        std::shared_ptr<const flint_mpoly::MpolyContext> ctx
    );
    
    // Check if monomial is divisible by any GB leading term
    bool is_divisible_by_leading_terms(
        const Monomial& monomial,
        const std::vector<MultivariatePolynomial<CoeffT>>& groebner_basis
    ) const;
    
    // Convert FLINT monomial ordering to our MonomialOrder enum
    static MonomialOrder flint_to_monomial_order(ordering_t flint_ord);
};

// Implementation

template<typename CoeffT>
std::shared_ptr<const PolynomialBasis<CoeffT>> PolynomialBasis<CoeffT>::create(
    const std::vector<MultivariatePolynomial<CoeffT>>& groebner_basis,
    size_t nvars,
    ulong prime,
    ordering_t ordering
) {
    if (groebner_basis.empty()) {
        throw std::invalid_argument("Gröbner basis cannot be empty");
    }
    
    if (nvars == 0) {
        throw std::invalid_argument("Number of variables must be positive");
    }
    
    // Create the basis object
    auto basis = std::shared_ptr<PolynomialBasis>(new PolynomialBasis());
    basis->nvars = nvars;
    basis->prime = prime;
    basis->ordering = ordering;
    
    // Create FLINT context
    auto ctx = std::make_shared<const flint_mpoly::MpolyContext>(nvars, prime, ordering);
    
    // Convert Gröbner basis to FLINT format
    basis->grobner_basis_flint.reserve(groebner_basis.size());
    for (const auto& poly : groebner_basis) {
        if (!poly.is_zero()) {
            basis->grobner_basis_flint.push_back(to_flint_mpoly(poly, ctx, prime));
        }
    }
    
    if (basis->grobner_basis_flint.empty()) {
        throw std::invalid_argument("Gröbner basis contains only zero polynomials");
    }
    
    // Compute quotient basis
    basis->compute_quotient_basis(groebner_basis, ctx);
    
    return basis;
}

template<typename CoeffT>
flint_mpoly::Mpoly PolynomialBasis<CoeffT>::to_flint_mpoly(
    const MultivariatePolynomial<CoeffT>& poly,
    std::shared_ptr<const flint_mpoly::MpolyContext> ctx,
    ulong prime
) {
    flint_mpoly::Mpoly result(ctx);
    
    if (poly.is_zero()) {
        return result; // Already initialized to zero
    }
    
    for (const auto& [monomial, coeff] : poly) {
        // Handle negative coefficients properly for modular arithmetic
        long long c = static_cast<long long>(coeff);
        long long remainder = c % static_cast<long long>(prime);
        ulong mod_coeff = (remainder < 0) ? 
            static_cast<ulong>(remainder + prime) : 
            static_cast<ulong>(remainder);
        
        if (mod_coeff != 0) {
            // Convert monomial exponents to FLINT format
            const auto& exps = monomial.exponents();
            if (exps.size() != ctx->nvars()) {
                throw std::runtime_error("Monomial dimension mismatch");
            }
            
            // Convert to ulong array for FLINT
            std::vector<ulong> flint_exps(exps.size());
            for (size_t i = 0; i < exps.size(); ++i) {
                if (exps[i] < 0) {
                    throw std::runtime_error("Negative exponents not allowed");
                }
                flint_exps[i] = static_cast<ulong>(exps[i]);
            }
            
            result.set_coeff_ui_monomial(mod_coeff, flint_exps.data());
        }
    }
    
    return result;
}

template<typename CoeffT>
void PolynomialBasis<CoeffT>::compute_quotient_basis(
    const std::vector<MultivariatePolynomial<CoeffT>>& groebner_basis,
    std::shared_ptr<const flint_mpoly::MpolyContext> ctx
) {
    quotient_basis_monomials.clear();
    monomial_to_index.clear();
    
    // BFS to find standard monomials
    std::queue<Monomial> queue;
    std::set<Monomial> visited;
    
    // Start with constant monomial (all exponents zero)
    Monomial constant(nvars);
    queue.push(constant);
    visited.insert(constant);
    
    while (!queue.empty()) {
        Monomial current = queue.front();
        queue.pop();
        
        // Check if current monomial is divisible by any GB leading term
        if (!is_divisible_by_leading_terms(current, groebner_basis)) {
            // This is a standard monomial
            size_t index = quotient_basis_monomials.size();
            quotient_basis_monomials.push_back(current);
            monomial_to_index[current] = index;
            
            // Add neighbors (multiply by each variable)
            for (size_t var = 0; var < nvars; ++var) {
                Monomial neighbor = current;
                neighbor[var]++;
                
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    queue.push(neighbor);
                }
            }
        }
    }
    
    // Sort quotient basis by degree and then by ordering
    MonomialOrder our_order = flint_to_monomial_order(ordering);
    MonomialComparator comp(our_order);
    
    // Create temporary vector with indices for stable sorting
    std::vector<std::pair<Monomial, size_t>> indexed_basis;
    for (size_t i = 0; i < quotient_basis_monomials.size(); ++i) {
        indexed_basis.emplace_back(quotient_basis_monomials[i], i);
    }
    
    std::sort(indexed_basis.begin(), indexed_basis.end(),
        [&comp](const auto& a, const auto& b) {
            return comp(a.first, b.first);
        });
    
    // Rebuild the structures with new ordering
    std::vector<Monomial> new_basis;
    std::unordered_map<Monomial, size_t> new_index_map;
    
    for (size_t i = 0; i < indexed_basis.size(); ++i) {
        new_basis.push_back(indexed_basis[i].first);
        new_index_map[indexed_basis[i].first] = i;
    }
    
    quotient_basis_monomials = std::move(new_basis);
    monomial_to_index = std::move(new_index_map);
}

template<typename CoeffT>
bool PolynomialBasis<CoeffT>::is_divisible_by_leading_terms(
    const Monomial& monomial,
    const std::vector<MultivariatePolynomial<CoeffT>>& groebner_basis
) const {
    for (const auto& poly : groebner_basis) {
        if (!poly.is_zero()) {
            const Monomial& leading_term = poly.leading_monomial();
            if (leading_term.divides(monomial)) {
                return true;
            }
        }
    }
    return false;
}

template<typename CoeffT>
MonomialOrder PolynomialBasis<CoeffT>::flint_to_monomial_order(ordering_t flint_ord) {
    switch (flint_ord) {
        case ORD_LEX:
            return MonomialOrder::LEX;
        case ORD_DEGLEX:
            return MonomialOrder::DEGLEX;
        case ORD_DEGREVLEX:
            return MonomialOrder::GREVLEX;
        default:
            return MonomialOrder::GREVLEX; // Safe default
    }
}

// Explicit template instantiations
extern template struct PolynomialBasis<int>;
extern template struct PolynomialBasis<long>;