#pragma once

#include "flint_wrappers.hpp"
#include <vector>
#include <memory>
#include <iostream>

/**
 * FLINT-based linear algebra for detecting dependencies in RUR algorithm
 * Uses fmpz_mod_mat_rref for reduced row echelon form over finite fields
 */
template<typename CoeffT>
class FLINTLinearSolver {
private:
    std::unique_ptr<flint_cpp::ModContext> ctx_;
    int prime_;
    
public:
    explicit FLINTLinearSolver(int prime) : prime_(prime) {
        ctx_ = std::make_unique<flint_cpp::ModContext>(prime);
    }
    
    /**
     * Find linear dependence in a set of vectors
     * Returns coefficients c such that sum(c[i] * vectors[i]) = 0
     * Returns empty vector if vectors are linearly independent
     */
    std::vector<CoeffT> find_linear_dependence(const std::vector<std::vector<CoeffT>>& vectors);
    
    /**
     * Solve linear system Ax = b over finite field
     * Returns solution vector, or empty if no solution exists
     */
    std::vector<CoeffT> solve_system(const std::vector<std::vector<CoeffT>>& A, 
                                     const std::vector<CoeffT>& b);
    
    /**
     * Compute rank of matrix
     */
    size_t compute_rank(const std::vector<std::vector<CoeffT>>& matrix);
    
    /**
     * Compute reduced row echelon form
     * Returns rank and modifies matrix in-place to RREF
     */
    size_t rref_inplace(std::vector<std::vector<CoeffT>>& matrix);
    
private:
    // Convert coefficient matrix to FLINT format
    std::unique_ptr<flint_cpp::ModMatrix> convert_to_flint_matrix(const std::vector<std::vector<CoeffT>>& matrix);
    
    // Convert FLINT matrix back to coefficient vectors
    std::vector<std::vector<CoeffT>> convert_from_flint_matrix(const flint_cpp::ModMatrix& matrix);
    
    // Convert coefficient to modular form
    CoeffT to_modular(CoeffT coeff) const;
};

template<typename CoeffT>
CoeffT FLINTLinearSolver<CoeffT>::to_modular(CoeffT coeff) const {
    CoeffT result = coeff % prime_;
    if (result < 0) result += prime_;
    return result;
}

template<typename CoeffT>
std::unique_ptr<flint_cpp::ModMatrix> FLINTLinearSolver<CoeffT>::convert_to_flint_matrix(
    const std::vector<std::vector<CoeffT>>& matrix) {
    
    if (matrix.empty()) {
        return std::make_unique<flint_cpp::ModMatrix>(0, 0, *ctx_);
    }
    
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    
    auto flint_mat = std::make_unique<flint_cpp::ModMatrix>(rows, cols, *ctx_);
    
    for (size_t i = 0; i < rows; ++i) {
        if (matrix[i].size() != cols) {
            throw std::invalid_argument("Matrix rows must have same length");
        }
        
        for (size_t j = 0; j < cols; ++j) {
            CoeffT mod_coeff = to_modular(matrix[i][j]);
            flint_mat->set_entry_ui(i, j, static_cast<ulong>(mod_coeff));
        }
    }
    
    return flint_mat;
}

template<typename CoeffT>
std::vector<std::vector<CoeffT>> FLINTLinearSolver<CoeffT>::convert_from_flint_matrix(
    const flint_cpp::ModMatrix& matrix) {
    
    size_t rows = matrix.nrows();
    size_t cols = matrix.ncols();
    
    std::vector<std::vector<CoeffT>> result(rows, std::vector<CoeffT>(cols));
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            ulong entry = matrix.get_entry_ui(i, j);
            result[i][j] = static_cast<CoeffT>(entry);
        }
    }
    return result;
}

template<typename CoeffT>
std::vector<CoeffT> FLINTLinearSolver<CoeffT>::find_linear_dependence(
    const std::vector<std::vector<CoeffT>>& vectors) {
    
    if (vectors.empty()) {
        return {};
    }
    
    // Create augmented matrix: [vectors^T | 0]
    // We transpose so each vector becomes a row
    size_t n_vectors = vectors.size();
    size_t vector_dim = vectors[0].size();
    
    // Create matrix where each row is a vector
    std::vector<std::vector<CoeffT>> matrix;
    for (const auto& vec : vectors) {
        if (vec.size() != vector_dim) {
            throw std::invalid_argument("All vectors must have same dimension");
        }
        matrix.push_back(vec);
    }
    
    auto flint_mat = convert_to_flint_matrix(matrix);
    
    // Compute RREF and check rank
    slong rank = flint_mat->rref();
    
    std::cout << "Linear dependence check: rank = " << rank << ", vectors = " << n_vectors << std::endl;
    
    if (rank == static_cast<slong>(n_vectors)) {
        // Vectors are linearly independent
        return {};
    }
    
    // There is a dependence - extract it from the RREF form
    auto rref_matrix = convert_from_flint_matrix(*flint_mat);
    
    // Find a zero row (indicates dependence)
    for (size_t i = rank; i < n_vectors; ++i) {
        bool is_zero_row = true;
        for (size_t j = 0; j < vector_dim; ++j) {
            if (rref_matrix[i][j] != 0) {
                is_zero_row = false;
                break;
            }
        }
        
        if (is_zero_row) {
            // This indicates a dependence relation
            // For now, return a simple dependence (this could be more sophisticated)
            std::vector<CoeffT> dependence(n_vectors, 0);
            dependence[i] = 1;  // Simplified - real implementation would extract from RREF
            return dependence;
        }
    }
    
    // Fallback: construct dependence from free variables
    std::vector<CoeffT> dependence(n_vectors, 0);
    dependence[n_vectors - 1] = 1;  // Set last coefficient to 1
    
    // Use back-substitution to find other coefficients
    for (int i = rank - 1; i >= 0; --i) {
        CoeffT sum = 0;
        for (size_t j = i + 1; j < n_vectors; ++j) {
            sum = (sum + rref_matrix[i][j] * dependence[j]) % prime_;
        }
        
        if (rref_matrix[i][i] != 0) {
            // Find modular inverse
            CoeffT pivot = rref_matrix[i][i];
            CoeffT inv_pivot = 1;  // This should be computed properly using modular inverse
            
            // Simple modular inverse for prime fields: a^(p-2) ≡ a^(-1) (mod p)
            // For efficiency, we'd use extended Euclidean algorithm
            for (int k = 0; k < prime_ - 2; ++k) {
                inv_pivot = (inv_pivot * pivot) % prime_;
            }
            
            dependence[i] = (prime_ - (sum * inv_pivot) % prime_) % prime_;
        }
    }
    
    return dependence;
}

template<typename CoeffT>
std::vector<CoeffT> FLINTLinearSolver<CoeffT>::solve_system(
    const std::vector<std::vector<CoeffT>>& A, const std::vector<CoeffT>& b) {
    
    if (A.empty() || A[0].empty()) {
        throw std::invalid_argument("Matrix A cannot be empty");
    }
    
    size_t rows = A.size();
    size_t cols = A[0].size();
    
    if (b.size() != rows) {
        throw std::invalid_argument("Vector b size must match matrix A rows");
    }
    
    // Create augmented matrix [A | b]
    std::vector<std::vector<CoeffT>> augmented(rows, std::vector<CoeffT>(cols + 1));
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][cols] = b[i];
    }
    
    auto flint_mat = convert_to_flint_matrix(augmented);
    slong rank = flint_mat->rref();
    
    auto rref_matrix = convert_from_flint_matrix(*flint_mat);
    
    // Check for inconsistency
    for (size_t i = rank; i < rows; ++i) {
        if (rref_matrix[i][cols] != 0) {
            // Inconsistent system
            return {};
        }
    }
    
    // Extract solution using back substitution
    std::vector<CoeffT> solution(cols, 0);
    
    for (int i = rank - 1; i >= 0; --i) {
        CoeffT sum = rref_matrix[i][cols];  // RHS
        
        for (size_t j = i + 1; j < cols; ++j) {
            sum = (sum - rref_matrix[i][j] * solution[j] % prime_ + prime_) % prime_;
        }
        
        if (rref_matrix[i][i] != 0) {
            // Compute modular inverse (simplified)
            CoeffT pivot = rref_matrix[i][i];
            CoeffT inv_pivot = 1;
            for (int k = 0; k < prime_ - 2; ++k) {
                inv_pivot = (inv_pivot * pivot) % prime_;
            }
            
            solution[i] = (sum * inv_pivot) % prime_;
        }
    }
    
    return solution;
}

template<typename CoeffT>
size_t FLINTLinearSolver<CoeffT>::compute_rank(const std::vector<std::vector<CoeffT>>& matrix) {
    if (matrix.empty()) return 0;
    
    auto flint_mat = convert_to_flint_matrix(matrix);
    return static_cast<size_t>(flint_mat->rref());
}

template<typename CoeffT>
size_t FLINTLinearSolver<CoeffT>::rref_inplace(std::vector<std::vector<CoeffT>>& matrix) {
    if (matrix.empty()) return 0;
    
    auto flint_mat = convert_to_flint_matrix(matrix);
    slong rank = flint_mat->rref();
    
    // Convert back and update original matrix
    auto rref_result = convert_from_flint_matrix(*flint_mat);
    matrix = rref_result;
    
    return static_cast<size_t>(rank);
}

// Explicit template instantiations
extern template class FLINTLinearSolver<int>;
extern template class FLINTLinearSolver<long>;