#pragma once

#include "quotient_ring.hpp"
#include "polynomial.hpp"
#include <vector>
#include <utility>

/**
 * Bivariate relation for parameterization
 * Represents a relation of the form: xi^a * T^b = sum(coeffs[j] * T^j)
 */
template<typename CoeffT>
struct BivariateRelation {
    size_t var_index;           // Which variable xi
    uint32_t var_degree;        // Degree in xi (usually 1)
    uint32_t sep_degree;        // Degree in separating element T
    std::vector<CoeffT> coeffs; // Coefficients of T^j terms
    
    BivariateRelation(size_t var_idx, uint32_t var_deg, uint32_t sep_deg, const std::vector<CoeffT>& c)
        : var_index(var_idx), var_degree(var_deg), sep_degree(sep_deg), coeffs(c) {}
};

/**
 * Bivariate lexicographic algorithm implementation
 * Corresponds to Julia's first_variable and biv_lex! functions
 */
template<typename CoeffT>
class BivariateAlgorithm {
private:
    QuotientRing<CoeffT>& quotient_ring_;
    int prime_;
    size_t nvars_;
    
    // Separating element data
    size_t separating_var_;
    std::vector<CoeffT> minimal_polynomial_;
    uint32_t separating_degree_;
    
    // Bivariate basis and relations
    std::vector<std::pair<uint32_t, uint32_t>> bivariate_basis_; // (T_degree, xi_degree) pairs
    std::vector<BivariateRelation<CoeffT>> relations_;
    
public:
    BivariateAlgorithm(QuotientRing<CoeffT>& qr, int prime)
        : quotient_ring_(qr), prime_(prime), nvars_(qr.num_variables()),
          separating_var_(nvars_ - 1), separating_degree_(0) {}
    
    // Compute minimal polynomial of separating element (Julia: first_variable, Lines 544-633)
    void compute_separating_element();
    
    // Compute bivariate relations for given variable (Julia: biv_lex!, Lines 635-684)
    std::vector<BivariateRelation<CoeffT>> compute_bivariate_relations(size_t var_index);
    
    // Extract parameterization for all variables
    std::vector<MultivariatePolynomial<CoeffT>> extract_parameterization();
    
    // Get minimal polynomial of separating element
    const std::vector<CoeffT>& get_minimal_polynomial() const { return minimal_polynomial_; }
    
    // Get degree of separating element
    uint32_t get_separating_degree() const { return separating_degree_; }
    
    // Set separating variable (default is last variable)
    void set_separating_variable(size_t var_idx) {
        if (var_idx >= nvars_) {
            throw std::invalid_argument("Separating variable index out of range");
        }
        separating_var_ = var_idx;
    }
    
private:
    // Helper: Gaussian elimination to detect linear dependence
    std::vector<CoeffT> solve_linear_system(const std::vector<std::vector<CoeffT>>& matrix);
    
    // Helper: compute powers of separating element
    std::vector<std::vector<CoeffT>> compute_separating_powers(uint32_t max_degree);
    
    // Helper: convert rational function to polynomial representation
    MultivariatePolynomial<CoeffT> rational_to_polynomial(const std::vector<CoeffT>& numerator,
                                                          const std::vector<CoeffT>& denominator);
    
    // Helper: compute derivative of polynomial
    std::vector<CoeffT> polynomial_derivative(const std::vector<CoeffT>& poly);
    
    // Helper: modular inverse
    CoeffT modular_inverse(CoeffT a) const;
    
    // Helper: extended Euclidean algorithm
    std::tuple<CoeffT, CoeffT, CoeffT> extended_gcd(CoeffT a, CoeffT b) const;
};

// Explicit template instantiation declarations
extern template class BivariateAlgorithm<int>;
extern template class BivariateAlgorithm<long>;