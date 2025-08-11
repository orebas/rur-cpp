#include "bivariate_algorithm.hpp"
#include "multiplication_tables.hpp" // For mul_var_quo
#include "polynomial_operations.hpp" // For modular_inverse
#include <iostream>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/flint.h>

namespace julia_rur {

// Helper to perform Gaussian reduction on a row using FLINT nmod_vec
bool gaussian_reduce_row(std::vector<ModularCoeff>& row, GaussianReductionContext& context, ModularCoeff prime) {
    if (row.empty()) return true;
    
    // Set up FLINT nmod context
    nmod_t mod;
    nmod_init(&mod, prime);
    
    size_t n = row.size();
    
    // Allocate FLINT vector for the row (will be freed automatically)
    nn_ptr row_vec = _nmod_vec_init(n);
    
    // Copy input row to FLINT format
    for (size_t j = 0; j < n; ++j) {
        row_vec[j] = row[j];
    }
    
    // Reduce against each row in the context matrix
    for (size_t i = 0; i < context.matrix.size(); ++i) {
        int32_t pivot_col = context.pivot_columns[i];
        if (pivot_col < 0 || pivot_col >= static_cast<int32_t>(n)) continue;

        ModularCoeff factor = row_vec[pivot_col];
        if (factor == 0) continue;
        
        // Allocate FLINT vector for the matrix row
        nn_ptr matrix_row_vec = _nmod_vec_init(n);
        
        // Copy matrix row to FLINT format
        for (size_t j = 0; j < n; ++j) {
            matrix_row_vec[j] = context.matrix[i][j];
        }

        // row_vec -= factor * matrix_row_vec using FLINT optimized operation
        // This is equivalent to: row = row - factor * matrix[i]
        _nmod_vec_scalar_addmul_nmod(row_vec, matrix_row_vec, n, nmod_neg(factor, mod), mod);
        
        _nmod_vec_clear(matrix_row_vec);
    }
    
    // Check if row is zero using FLINT
    bool is_zero = _nmod_vec_is_zero(row_vec, n);
    
    // Copy result back to std::vector
    for (size_t j = 0; j < n; ++j) {
        row[j] = row_vec[j];
    }
    
    _nmod_vec_clear(row_vec);
    return is_zero;
}


BivLexResult biv_lex(const std::vector<std::vector<ModularCoeff>>& t_v,
                     const std::vector<std::vector<int32_t>>& i_xw,
                     GaussianReductionContext& context,
                     const std::vector<std::vector<ModularCoeff>>& initial_free_set,
                     int32_t var_index,
                     ModularCoeff prime) {

    BivLexResult result;
    context.matrix = initial_free_set;

    // Initialize pivots for the initial matrix (assumed to be in row-echelon form)
    context.pivot_columns.assign(context.matrix.size(), -1);
    for(size_t i = 0; i < context.matrix.size(); ++i) {
        for(size_t j = 0; j < context.matrix[i].size(); ++j) {
            if (context.matrix[i][j] != 0) {
                context.pivot_columns[i] = j;
                break;
            }
        }
    }
    
    std::vector<std::vector<ModularCoeff>> current_free_set = initial_free_set;
    std::vector<BivariateMonomial> current_monomials;
    for (size_t i = 0; i < initial_free_set.size(); ++i) {
        current_monomials.push_back({(int32_t)i, 0});
    }
    result.monomial_basis = current_monomials;

    while (!current_free_set.empty()) {
        std::vector<std::vector<ModularCoeff>> next_free_set;
        std::vector<BivariateMonomial> next_monomials;

        for (size_t i = 0; i < current_free_set.size(); ++i) {
            auto& vec = current_free_set[i];
            auto mon = current_monomials[i];

            std::vector<ModularCoeff> new_vec(vec.size());
            mul_var_quo(new_vec, vec, var_index, i_xw, t_v, prime);

            BivariateMonomial new_mon = {mon.deg_T, mon.deg_xi + 1};

            if (!gaussian_reduce_row(new_vec, context, prime)) {
                // Not linearly dependent, add to basis
                
                // Normalize the new row to have a pivot of 1
                int32_t pivot = -1;
                for(size_t j = 0; j < new_vec.size(); ++j) {
                    if (new_vec[j] != 0) {
                        pivot = j;
                        break;
                    }
                }

                if (pivot != -1) {
                    // Set up FLINT context for normalization and reduction
                    nmod_t mod;
                    nmod_init(&mod, prime);
                    size_t n = new_vec.size();
                    
                    // Normalize the new row using FLINT
                    ModularCoeff inv = modular_inverse(new_vec[pivot], prime);
                    nn_ptr new_vec_ptr = _nmod_vec_init(n);
                    for (size_t j = 0; j < n; ++j) {
                        new_vec_ptr[j] = new_vec[j];
                    }
                    _nmod_vec_scalar_mul_nmod(new_vec_ptr, new_vec_ptr, n, inv, mod);
                    
                    // Copy normalized result back
                    for (size_t j = 0; j < n; ++j) {
                        new_vec[j] = new_vec_ptr[j];
                    }
                    
                    // Reduce other rows in the matrix using FLINT
                    for(auto& row : context.matrix) {
                        ModularCoeff factor = row[pivot];
                        if (factor == 0) continue;
                        
                        // Convert row to FLINT format
                        nn_ptr row_ptr = _nmod_vec_init(n);
                        for (size_t j = 0; j < n; ++j) {
                            row_ptr[j] = row[j];
                        }
                        
                        // row -= factor * new_vec using FLINT
                        _nmod_vec_scalar_addmul_nmod(row_ptr, new_vec_ptr, n, nmod_neg(factor, mod), mod);
                        
                        // Copy result back
                        for (size_t j = 0; j < n; ++j) {
                            row[j] = row_ptr[j];
                        }
                        
                        _nmod_vec_clear(row_ptr);
                    }
                    
                    _nmod_vec_clear(new_vec_ptr);

                    context.matrix.push_back(new_vec);
                    context.pivot_columns.push_back(pivot);
                    next_free_set.push_back(new_vec);
                    next_monomials.push_back(new_mon);
                    result.monomial_basis.push_back(new_mon);
                }

            } else {
                // Linearly dependent, we found a generator
                result.leading_monomials.push_back(new_mon);
                
                // To find the generator coefficients, we need to express the dependency
                // The original vector before reduction gives us the coefficients
                // Recompute the original vector
                std::vector<ModularCoeff> original_vec(vec.size());
                mul_var_quo(original_vec, vec, var_index, i_xw, t_v, prime);
                
                // The generator is the original vector since it represents
                // the linear combination that equals the new monomial
                result.generators.push_back(original_vec);
            }
        }
        current_free_set = next_free_set;
        current_monomials = next_monomials;
    }

    return result;
}

} // namespace julia_rur