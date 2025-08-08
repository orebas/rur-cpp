#include "bivariate_algorithm.hpp"
#include "multiplication_tables.hpp" // For mul_var_quo
#include "polynomial_operations.hpp" // For modular_inverse
#include <iostream>

namespace julia_rur {

// Helper to perform Gaussian reduction on a row
bool gaussian_reduce_row(std::vector<ModularCoeff>& row, GaussianReductionContext& context, ModularCoeff prime) {
    for (size_t i = 0; i < context.matrix.size(); ++i) {
        int32_t pivot_col = context.pivot_columns[i];
        if (pivot_col < 0) continue;

        ModularCoeff factor = row[pivot_col];
        if (factor == 0) continue;

        // row -= factor * matrix[i]
        for (size_t j = 0; j < row.size(); ++j) {
            AccModularCoeff prod = static_cast<AccModularCoeff>(factor) * context.matrix[i][j];
            row[j] = (row[j] + prime - (prod % prime)) % prime;
        }
    }

    // Check if row is zero
    for (ModularCoeff c : row) {
        if (c != 0) return false;
    }
    return true;
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
                    ModularCoeff inv = modular_inverse(new_vec[pivot], prime);
                    for (size_t j = 0; j < new_vec.size(); ++j) {
                        new_vec[j] = (static_cast<AccModularCoeff>(new_vec[j]) * inv) % prime;
                    }
                    
                    // Reduce other rows in the matrix
                    for(auto& row : context.matrix) {
                        ModularCoeff factor = row[pivot];
                        if (factor == 0) continue;
                        for(size_t j = 0; j < row.size(); ++j) {
                            AccModularCoeff prod = static_cast<AccModularCoeff>(factor) * new_vec[j];
                            row[j] = (row[j] + prime - (prod % prime)) % prime;
                        }
                    }

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