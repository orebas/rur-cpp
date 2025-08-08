#pragma once

#include "flint_mpoly_wrappers.hpp"
#include "polynomial_basis.hpp"
#include "polynomial.hpp"
#include <memory>
#include <vector>

/**
 * Production-quality polynomial reducer using FLINT multivariate polynomials
 * 
 * This class provides fast, mathematically correct polynomial reduction
 * using FLINT's nmod_mpoly_divrem_ideal for simultaneous division by 
 * multiple polynomials (Gröbner basis).
 * 
 * Thread Safety: Each FLINTReducer instance should be used by only one thread.
 * Create separate instances for multi-threaded applications.
 */
template<typename CoeffT>
class FLINTReducer {
private:
    std::shared_ptr<const flint_mpoly::MpolyContext> ctx_;
    std::shared_ptr<const PolynomialBasis<CoeffT>> basis_;
    ulong prime_;
    
public:
    /**
     * Constructor takes pre-computed polynomial basis
     * 
     * @param basis Pre-computed PolynomialBasis containing GB and quotient basis
     * @param prime Prime modulus for finite field arithmetic  
     */
    explicit FLINTReducer(std::shared_ptr<const PolynomialBasis<CoeffT>> basis);
    
    /**
     * Reduce polynomial modulo the Gröbner basis
     * 
     * Computes the unique normal form of the input polynomial with respect
     * to the Gröbner basis. The result is expressed as coefficients with
     * respect to the quotient basis.
     * 
     * @param poly_in Input polynomial to reduce
     * @param result_coeffs Output array (must be pre-allocated to quotient basis size)
     */
    void reduce(const MultivariatePolynomial<CoeffT>& poly_in, CoeffT* result_coeffs);
    
    /**
     * Reduce polynomial and return coefficient vector
     * 
     * Convenience method that returns std::vector instead of writing to array
     * 
     * @param poly_in Input polynomial to reduce
     * @return Coefficient vector with respect to quotient basis
     */
    std::vector<CoeffT> reduce(const MultivariatePolynomial<CoeffT>& poly_in);
    
    /**
     * Get quotient basis information
     */
    const std::vector<Monomial>& get_quotient_basis() const {
        return basis_->quotient_basis_monomials;
    }
    
    size_t quotient_dimension() const {
        return basis_->quotient_basis_monomials.size();
    }
    
    ulong prime() const { return prime_; }
    size_t nvars() const { return basis_->nvars; }
    
private:
    /**
     * Convert our polynomial to FLINT format with modular reduction
     */
    flint_mpoly::Mpoly to_flint(const MultivariatePolynomial<CoeffT>& poly) const;
    
    /**
     * Extract coefficients from reduced FLINT polynomial
     * 
     * The input polynomial should only contain monomials from the quotient basis.
     * If any other monomial is found, this indicates a bug in the reduction.
     */
    void extract_coefficients(const flint_mpoly::Mpoly& reduced_poly, CoeffT* result) const;
    
    /**
     * Validate coefficient for conversion back to CoeffT
     */
    CoeffT validate_and_convert_coeff(ulong flint_coeff) const;
};

// Implementation

template<typename CoeffT>
FLINTReducer<CoeffT>::FLINTReducer(std::shared_ptr<const PolynomialBasis<CoeffT>> basis) 
    : basis_(basis), prime_(basis->prime) {
    
    if (!basis_) {
        throw std::invalid_argument("PolynomialBasis cannot be null");
    }
    
    if (basis_->grobner_basis_flint.empty()) {
        throw std::invalid_argument("Gröbner basis cannot be empty");
    }
    
    if (basis_->quotient_basis_monomials.empty()) {
        throw std::invalid_argument("Quotient basis cannot be empty");
    }
    
    // Get context from the first polynomial in the basis
    ctx_ = basis_->grobner_basis_flint[0].context();
    
    // Verify all polynomials use the same context
    for (const auto& poly : basis_->grobner_basis_flint) {
        if (poly.context() != ctx_) {
            throw std::runtime_error("All Gröbner basis polynomials must use same context");
        }
    }
}

template<typename CoeffT>
void FLINTReducer<CoeffT>::reduce(const MultivariatePolynomial<CoeffT>& poly_in, CoeffT* result_coeffs) {
    if (!result_coeffs) {
        throw std::invalid_argument("Result coefficients array cannot be null");
    }
    
    // Initialize result to zero
    for (size_t i = 0; i < quotient_dimension(); ++i) {
        result_coeffs[i] = CoeffT(0);
    }
    
    // Convert input to FLINT format
    auto flint_input = to_flint(poly_in);
    
    // Prepare quotient polynomials (we need them for the API even though we discard them)
    std::vector<flint_mpoly::Mpoly> quotients;
    quotients.reserve(basis_->grobner_basis_flint.size());
    for (size_t i = 0; i < basis_->grobner_basis_flint.size(); ++i) {
        quotients.emplace_back(ctx_);
    }
    
    // Prepare remainder polynomial
    flint_mpoly::Mpoly remainder(ctx_);
    
    // Perform the ideal division - this is the core mathematical operation
    flint_mpoly::divrem_ideal(quotients, remainder, flint_input, basis_->grobner_basis_flint);
    
    // Extract coefficients from remainder
    extract_coefficients(remainder, result_coeffs);
}

template<typename CoeffT>
std::vector<CoeffT> FLINTReducer<CoeffT>::reduce(const MultivariatePolynomial<CoeffT>& poly_in) {
    std::vector<CoeffT> result(quotient_dimension());
    reduce(poly_in, result.data());
    return result;
}

template<typename CoeffT>
flint_mpoly::Mpoly FLINTReducer<CoeffT>::to_flint(const MultivariatePolynomial<CoeffT>& poly) const {
    flint_mpoly::Mpoly result(ctx_);
    
    if (poly.is_zero()) {
        return result; // Already initialized to zero
    }
    
    for (const auto& [monomial, coeff] : poly) {
        // Handle negative coefficients properly for modular arithmetic
        long long c = static_cast<long long>(coeff);
        long long remainder = c % static_cast<long long>(prime_);
        ulong mod_coeff = (remainder < 0) ? 
            static_cast<ulong>(remainder + prime_) : 
            static_cast<ulong>(remainder);
        
        if (mod_coeff != 0) {
            // Convert monomial exponents to FLINT format
            const auto& exps = monomial.exponents();
            if (exps.size() != nvars()) {
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
void FLINTReducer<CoeffT>::extract_coefficients(const flint_mpoly::Mpoly& reduced_poly, CoeffT* result) const {
    // The reduced polynomial should only contain monomials from the quotient basis
    for (const auto& term : reduced_poly) {
        // Convert FLINT exponent vector to our Monomial
        std::vector<int> exps(term.exps.begin(), term.exps.end());
        Monomial monomial(exps);
        
        // Look up this monomial in the quotient basis
        auto it = basis_->monomial_to_index.find(monomial);
        if (it == basis_->monomial_to_index.end()) {
            // This should never happen if reduction is correct
            throw std::runtime_error("Reduced polynomial contains monomial not in quotient basis: " + 
                                   monomial.to_string());
        }
        
        size_t index = it->second;
        result[index] = validate_and_convert_coeff(term.coeff);
    }
}

template<typename CoeffT>
CoeffT FLINTReducer<CoeffT>::validate_and_convert_coeff(ulong flint_coeff) const {
    // Convert back from finite field to our coefficient type
    // For now, just do a direct cast (this assumes CoeffT can represent the value)
    if (flint_coeff >= prime_) {
        throw std::runtime_error("FLINT coefficient outside expected range");
    }
    
    return static_cast<CoeffT>(flint_coeff);
}

// Explicit template instantiations
extern template class FLINTReducer<int>;
extern template class FLINTReducer<long>;