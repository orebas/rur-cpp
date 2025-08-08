#include "bivariate_lex.hpp"
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>

template<typename CoeffT>
std::tuple<CoeffT, CoeffT, CoeffT> BivariateAlgorithm<CoeffT>::extended_gcd(CoeffT a, CoeffT b) const {
    if (b == 0) {
        return std::make_tuple(a, 1, 0);
    }
    
    auto [gcd, x1, y1] = extended_gcd(b, a % b);
    CoeffT x = y1;
    CoeffT y = x1 - (a / b) * y1;
    
    return std::make_tuple(gcd, x, y);
}

template<typename CoeffT>
CoeffT BivariateAlgorithm<CoeffT>::modular_inverse(CoeffT a) const {
    auto [gcd, x, y] = extended_gcd(a, prime_);
    
    if (gcd != 1) {
        throw std::runtime_error("Modular inverse does not exist");
    }
    
    return (x % prime_ + prime_) % prime_;
}

template<typename CoeffT>
void BivariateAlgorithm<CoeffT>::compute_separating_element() {
    if (!quotient_ring_.are_tables_computed()) {
        throw std::runtime_error("Quotient ring multiplication tables not computed");
    }
    
    size_t dim = quotient_ring_.dimension();
    
    // Start with separating variable as basis vector
    std::vector<CoeffT> current(dim, 0);
    
    // Find index of separating variable in quotient basis
    PP sep_monomial(nvars_, 0);
    sep_monomial[separating_var_] = 1;
    
    const auto& basis = quotient_ring_.get_quotient_basis();
    auto it = std::find(basis.begin(), basis.end(), sep_monomial);
    
    if (it == basis.end()) {
        throw std::runtime_error("Separating variable not found in quotient basis");
    }
    
    size_t sep_index = std::distance(basis.begin(), it);
    current[sep_index] = 1;
    
    // Compute powers of T until linear dependence is found
    std::vector<std::vector<CoeffT>> powers;
    powers.push_back(std::vector<CoeffT>(dim, 0));  // T^0 = 1
    powers[0][0] = 1;  // Assuming constant monomial is at index 0
    
    powers.push_back(current);  // T^1
    
    // Iteratively compute T^k = T * T^(k-1)
    for (uint32_t degree = 2; degree <= dim; ++degree) {
        std::vector<CoeffT> next_power = quotient_ring_.multiply_by_variable(powers.back(), separating_var_);
        powers.push_back(next_power);
        
        // Check for linear dependence using Gaussian elimination
        std::vector<std::vector<CoeffT>> matrix;
        for (const auto& power : powers) {
            matrix.push_back(power);
        }
        
        // Try to solve the system: c0*T^0 + c1*T^1 + ... + cd*T^d = 0
        std::vector<CoeffT> coeffs = solve_linear_system(matrix);
        
        if (!coeffs.empty()) {
            // Found linear dependence - extract minimal polynomial
            minimal_polynomial_ = coeffs;
            separating_degree_ = degree;
            
            std::cout << "Found separating element of degree " << degree << std::endl;
            return;
        }
    }
    
    throw std::runtime_error("Failed to find minimal polynomial for separating element");
}

template<typename CoeffT>
std::vector<CoeffT> BivariateAlgorithm<CoeffT>::solve_linear_system(const std::vector<std::vector<CoeffT>>& matrix) {
    if (matrix.empty()) return {};
    
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    
    // Create augmented matrix for homogeneous system
    std::vector<std::vector<CoeffT>> aug_matrix(rows, std::vector<CoeffT>(cols + 1, 0));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            aug_matrix[i][j] = matrix[i][j];
        }
    }
    
    // Gaussian elimination
    size_t pivot_row = 0;
    for (size_t col = 0; col < cols && pivot_row < rows; ++col) {
        // Find pivot
        size_t best_row = pivot_row;
        for (size_t row = pivot_row + 1; row < rows; ++row) {
            if (abs(aug_matrix[row][col]) > abs(aug_matrix[best_row][col])) {
                best_row = row;
            }
        }
        
        if (aug_matrix[best_row][col] == 0) continue;
        
        // Swap rows
        if (best_row != pivot_row) {
            std::swap(aug_matrix[pivot_row], aug_matrix[best_row]);
        }
        
        // Eliminate
        CoeffT pivot = aug_matrix[pivot_row][col];
        CoeffT inv_pivot = modular_inverse(pivot);
        
        for (size_t row = pivot_row + 1; row < rows; ++row) {
            if (aug_matrix[row][col] != 0) {
                CoeffT factor = (aug_matrix[row][col] * inv_pivot) % prime_;
                for (size_t j = col; j < cols + 1; ++j) {
                    aug_matrix[row][j] = (aug_matrix[row][j] - factor * aug_matrix[pivot_row][j] % prime_ + prime_) % prime_;
                }
            }
        }
        
        pivot_row++;
    }
    
    // Check if we have a free variable (indicates linear dependence)
    if (pivot_row < rows) {
        // Found linear dependence - construct solution
        std::vector<CoeffT> solution(rows, 0);
        solution[rows - 1] = 1;  // Set last coefficient to 1
        
        // Back substitution to find other coefficients
        for (int i = pivot_row - 1; i >= 0; --i) {
            CoeffT sum = 0;
            for (size_t j = i + 1; j < rows; ++j) {
                sum = (sum + aug_matrix[i][j] * solution[j]) % prime_;
            }
            
            if (aug_matrix[i][i] != 0) {
                solution[i] = (prime_ - sum * modular_inverse(aug_matrix[i][i]) % prime_) % prime_;
            }
        }
        
        return solution;
    }
    
    return {};  // No linear dependence found
}

template<typename CoeffT>
std::vector<BivariateRelation<CoeffT>> BivariateAlgorithm<CoeffT>::compute_bivariate_relations(size_t var_index) {
    if (minimal_polynomial_.empty()) {
        throw std::runtime_error("Separating element not computed");
    }
    
    if (var_index >= nvars_) {
        throw std::invalid_argument("Variable index out of range");
    }
    
    std::vector<BivariateRelation<CoeffT>> relations;
    
    // Compute powers of separating element
    std::vector<std::vector<CoeffT>> T_powers = compute_separating_powers(separating_degree_);
    
    // Start with xi as vector in quotient basis
    size_t dim = quotient_ring_.dimension();
    std::vector<CoeffT> xi_vec(dim, 0);
    
    PP xi_monomial(nvars_, 0);
    xi_monomial[var_index] = 1;
    
    const auto& basis = quotient_ring_.get_quotient_basis();
    auto it = std::find(basis.begin(), basis.end(), xi_monomial);
    
    if (it != basis.end()) {
        size_t xi_index = std::distance(basis.begin(), it);
        xi_vec[xi_index] = 1;
    }
    
    // Maintain bivariate basis: start with univariate T basis
    bivariate_basis_.clear();
    for (uint32_t t_deg = 0; t_deg < separating_degree_; ++t_deg) {
        bivariate_basis_.emplace_back(t_deg, 0);  // (T^t_deg, xi^0)
    }
    
    // Iteratively add xi * (previous elements) and reduce
    std::vector<std::vector<CoeffT>> current_basis = T_powers;
    
    for (uint32_t xi_deg = 1; xi_deg <= separating_degree_; ++xi_deg) {
        std::vector<std::vector<CoeffT>> next_basis;
        
        for (const auto& vec : current_basis) {
            // Multiply by xi
            std::vector<CoeffT> xi_product = quotient_ring_.multiply_by_variable(vec, var_index);
            
            // Try to reduce using existing bivariate basis
            bool found_relation = false;
            
            // Check if xi_product can be expressed in terms of T powers
            for (uint32_t t_deg = 0; t_deg < separating_degree_; ++t_deg) {
                // This is a simplified check - full implementation would use Gaussian elimination
                bool linearly_dependent = true;  // Placeholder
                
                if (linearly_dependent) {
                    // Found a relation: xi^xi_deg * T^some_deg = linear_combination_of_T_powers
                    std::vector<CoeffT> relation_coeffs = xi_product;  // Simplified
                    relations.emplace_back(var_index, xi_deg, t_deg, relation_coeffs);
                    found_relation = true;
                    break;
                }
            }
            
            if (!found_relation) {
                next_basis.push_back(xi_product);
                bivariate_basis_.emplace_back(0, xi_deg);  // Placeholder T degree
            }
        }
        
        current_basis = next_basis;
        
        if (current_basis.empty()) {
            break;  // No more elements to add
        }
    }
    
    return relations;
}

template<typename CoeffT>
std::vector<std::vector<CoeffT>> BivariateAlgorithm<CoeffT>::compute_separating_powers(uint32_t max_degree) {
    std::vector<std::vector<CoeffT>> powers;
    size_t dim = quotient_ring_.dimension();
    
    // T^0 = 1
    std::vector<CoeffT> t0(dim, 0);
    t0[0] = 1;  // Assuming constant is at index 0
    powers.push_back(t0);
    
    // T^1 = separating variable
    std::vector<CoeffT> t1(dim, 0);
    PP sep_monomial(nvars_, 0);
    sep_monomial[separating_var_] = 1;
    
    const auto& basis = quotient_ring_.get_quotient_basis();
    auto it = std::find(basis.begin(), basis.end(), sep_monomial);
    if (it != basis.end()) {
        size_t sep_index = std::distance(basis.begin(), it);
        t1[sep_index] = 1;
    }
    powers.push_back(t1);
    
    // T^k = T * T^(k-1)
    for (uint32_t k = 2; k < max_degree; ++k) {
        std::vector<CoeffT> tk = quotient_ring_.multiply_by_variable(powers.back(), separating_var_);
        powers.push_back(tk);
    }
    
    return powers;
}

template<typename CoeffT>
std::vector<MultivariatePolynomial<CoeffT>> BivariateAlgorithm<CoeffT>::extract_parameterization() {
    std::vector<MultivariatePolynomial<CoeffT>> parameterization(nvars_);
    
    // Compute derivative of minimal polynomial for denominator
    std::vector<CoeffT> derivative = polynomial_derivative(minimal_polynomial_);
    
    for (size_t var = 0; var < nvars_; ++var) {
        if (var == separating_var_) {
            // Separating variable parameterizes as T itself
            MultivariatePolynomial<CoeffT> T_poly;
            Monomial T_mon(nvars_);
            T_mon[separating_var_] = 1;
            T_poly.set_coefficient(T_mon, 1);
            parameterization[var] = T_poly;
        } else {
            // Compute bivariate relations for this variable
            auto relations = compute_bivariate_relations(var);
            
            if (!relations.empty()) {
                // Extract rational function from first relation
                // xi = numerator(T) / derivative(T)
                parameterization[var] = rational_to_polynomial(relations[0].coeffs, derivative);
            }
        }
    }
    
    return parameterization;
}

template<typename CoeffT>
std::vector<CoeffT> BivariateAlgorithm<CoeffT>::polynomial_derivative(const std::vector<CoeffT>& poly) {
    if (poly.size() <= 1) {
        return {0};
    }
    
    std::vector<CoeffT> derivative(poly.size() - 1);
    for (size_t i = 1; i < poly.size(); ++i) {
        derivative[i - 1] = (static_cast<CoeffT>(i) * poly[i]) % prime_;
    }
    
    return derivative;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT> BivariateAlgorithm<CoeffT>::rational_to_polynomial(
    const std::vector<CoeffT>& numerator, const std::vector<CoeffT>& denominator) {
    
    // This is a placeholder implementation
    // Full implementation would handle rational function arithmetic properly
    MultivariatePolynomial<CoeffT> result;
    
    // For now, just return the numerator as a polynomial in the separating variable
    for (size_t i = 0; i < numerator.size(); ++i) {
        if (numerator[i] != 0) {
            Monomial mon(nvars_);
            mon[separating_var_] = i;
            result.set_coefficient(mon, numerator[i]);
        }
    }
    
    return result;
}

// Explicit template instantiations
template class BivariateAlgorithm<int>;
template class BivariateAlgorithm<long>;