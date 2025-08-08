#pragma once

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <memory>
#include <stdexcept>
#include <vector>

/**
 * RAII wrappers for FLINT dense linear algebra over finite fields
 * Provides safe memory management and clean C++ integration for nmod_mat operations
 */

namespace flint_linalg {

/**
 * RAII wrapper for nmod_mat_t (dense matrix over Z/pZ)
 */
class NModMat {
private:
    nmod_mat_t mat_;
    ulong prime_;
    
public:
    NModMat(slong rows, slong cols, ulong prime) : prime_(prime) {
        if (rows < 0 || cols < 0) {
            throw std::invalid_argument("Matrix dimensions must be non-negative");
        }
        if (prime < 2) {
            throw std::invalid_argument("Prime must be >= 2");
        }
        nmod_mat_init(mat_, rows, cols, prime);
    }
    
    ~NModMat() {
        nmod_mat_clear(mat_);
    }
    
    // Non-copyable for now (can add copy constructor if needed)
    NModMat(const NModMat&) = delete;
    NModMat& operator=(const NModMat&) = delete;
    
    // Movable
    NModMat(NModMat&& other) noexcept : prime_(other.prime_) {
        *mat_ = *other.mat_;
        // Re-initialize moved-from object to valid empty state
        nmod_mat_init(other.mat_, 0, 0, prime_);
    }
    
    NModMat& operator=(NModMat&& other) {
        if (this != &other) {
            if (prime_ != other.prime_) {
                throw std::runtime_error("Cannot move matrices with different moduli");
            }
            
            nmod_mat_clear(mat_);
            *mat_ = *other.mat_;
            nmod_mat_init(other.mat_, 0, 0, prime_);
        }
        return *this;
    }
    
    // Access
    nmod_mat_t& get() { return mat_; }
    const nmod_mat_t& get() const { return mat_; }
    
    ulong prime() const { return prime_; }
    slong rows() const { return nmod_mat_nrows(mat_); }
    slong cols() const { return nmod_mat_ncols(mat_); }
    
    // Element access with bounds checking
    ulong get_entry(slong i, slong j) const {
        if (i < 0 || i >= rows() || j < 0 || j >= cols()) {
            throw std::out_of_range("Matrix index out of range");
        }
        return nmod_mat_get_entry(mat_, i, j);
    }
    
    void set_entry(slong i, slong j, ulong val) {
        if (i < 0 || i >= rows() || j < 0 || j >= cols()) {
            throw std::out_of_range("Matrix index out of range");
        }
        nmod_mat_set_entry(mat_, i, j, val);
    }
    
    // Matrix operations
    void zero() {
        nmod_mat_zero(mat_);
    }
    
    void one() {
        if (rows() != cols()) {
            throw std::runtime_error("Identity matrix requires square matrix");
        }
        nmod_mat_one(mat_);
    }
    
    bool is_zero() const {
        return nmod_mat_is_zero(mat_);
    }
    
    // Create window view (no memory allocation)
    NModMat window(slong r1, slong c1, slong r2, slong c2) const {
        if (r1 < 0 || r2 > rows() || c1 < 0 || c2 > cols() || r1 >= r2 || c1 >= c2) {
            throw std::out_of_range("Invalid window bounds");
        }
        
        NModMat result(0, 0, prime_);  // Empty initialization
        nmod_mat_clear(result.mat_);   // Clear the empty matrix
        nmod_mat_window_init(result.mat_, mat_, r1, c1, r2, c2);
        return result;
    }
};

/**
 * RAII wrapper for nmod_poly_t (polynomial over Z/pZ)
 */
class NModPoly {
private:
    nmod_poly_t poly_;
    ulong prime_;
    
public:
    explicit NModPoly(ulong prime) : prime_(prime) {
        if (prime < 2) {
            throw std::invalid_argument("Prime must be >= 2");
        }
        nmod_poly_init(poly_, prime);
    }
    
    ~NModPoly() {
        nmod_poly_clear(poly_);
    }
    
    // Non-copyable for now
    NModPoly(const NModPoly&) = delete;
    NModPoly& operator=(const NModPoly&) = delete;
    
    // Movable
    NModPoly(NModPoly&& other) noexcept : prime_(other.prime_) {
        *poly_ = *other.poly_;
        nmod_poly_init(other.poly_, prime_);
    }
    
    NModPoly& operator=(NModPoly&& other) {
        if (this != &other) {
            if (prime_ != other.prime_) {
                throw std::runtime_error("Cannot move polynomials with different moduli");
            }
            
            nmod_poly_clear(poly_);
            *poly_ = *other.poly_;
            nmod_poly_init(other.poly_, prime_);
        }
        return *this;
    }
    
    // Access
    nmod_poly_t& get() { return poly_; }
    const nmod_poly_t& get() const { return poly_; }
    
    ulong prime() const { return prime_; }
    slong degree() const { return nmod_poly_degree(poly_); }
    slong length() const { return nmod_poly_length(poly_); }
    
    // Coefficient access
    ulong get_coeff(slong n) const {
        return nmod_poly_get_coeff_ui(poly_, n);
    }
    
    void set_coeff(slong n, ulong val) {
        nmod_poly_set_coeff_ui(poly_, n, val);
    }
    
    // Polynomial operations
    void zero() {
        nmod_poly_zero(poly_);
    }
    
    void one() {
        nmod_poly_one(poly_);
    }
    
    bool is_zero() const {
        return nmod_poly_is_zero(poly_);
    }
};

/**
 * Core linear algebra operations
 */
class LinearSolver {
private:
    ulong prime_;
    
public:
    explicit LinearSolver(ulong prime) : prime_(prime) {
        if (prime < 2) {
            throw std::invalid_argument("Prime must be >= 2");
        }
    }
    
    /**
     * Compute reduced row echelon form
     * Returns the rank of the matrix
     */
    slong rref(NModMat& matrix) const {
        if (matrix.prime() != prime_) {
            throw std::runtime_error("Matrix modulus must match solver modulus");
        }
        return nmod_mat_rref(matrix.get());
    }
    
    /**
     * Compute nullspace of matrix
     * Returns basis vectors as columns of the result matrix
     */
    slong nullspace(NModMat& nullspace_basis, const NModMat& matrix) const {
        if (matrix.prime() != prime_ || nullspace_basis.prime() != prime_) {
            throw std::runtime_error("Matrix modulus must match solver modulus");
        }
        return nmod_mat_nullspace(nullspace_basis.get(), matrix.get());
    }
    
    /**
     * Solve linear system AX = B
     * Returns true if solution exists, false otherwise
     */
    bool solve(NModMat& solution, const NModMat& matrix, const NModMat& rhs) const {
        if (matrix.prime() != prime_ || rhs.prime() != prime_ || solution.prime() != prime_) {
            throw std::runtime_error("Matrix modulus must match solver modulus");
        }
        return nmod_mat_solve(solution.get(), matrix.get(), rhs.get()) != 0;
    }
    
    /**
     * Compute minimal polynomial of a matrix using nullspace method
     * For a square matrix A, finds the monic polynomial p(x) of minimal degree
     * such that p(A) = 0
     */
    void minimal_polynomial(NModPoly& min_poly, const NModMat& matrix) const {
        if (matrix.prime() != prime_) {
            throw std::runtime_error("Matrix modulus must match solver modulus");
        }
        if (matrix.rows() != matrix.cols()) {
            throw std::runtime_error("Minimal polynomial requires square matrix");
        }
        
        slong n = matrix.rows();
        
        // Build companion-style matrix for powers of A
        // [I | A | A^2 | ... | A^n] and find nullspace
        NModMat powers(n, n * (n + 1), prime_);
        powers.zero();
        
        // Set first block to identity
        for (slong i = 0; i < n; i++) {
            powers.set_entry(i, i, 1);
        }
        
        // Compute powers A^k iteratively
        NModMat current_power(n, n, prime_);
        current_power.one();  // A^0 = I
        
        for (slong k = 1; k <= n; k++) {
            // Multiply current_power by matrix to get A^k
            nmod_mat_mul(current_power.get(), current_power.get(), matrix.get());
            
            // Copy A^k to appropriate block in powers matrix
            for (slong i = 0; i < n; i++) {
                for (slong j = 0; j < n; j++) {
                    powers.set_entry(i, k * n + j, current_power.get_entry(i, j));
                }
            }
        }
        
        // Find nullspace - should be 1-dimensional for minimal polynomial
        NModMat nullspace_basis(n * (n + 1), n * (n + 1), prime_);
        slong nullity = nullspace(nullspace_basis, powers);
        
        if (nullity == 0) {
            throw std::runtime_error("No minimal polynomial found - matrix may be degenerate");
        }
        
        // Extract coefficients from first nullspace vector
        min_poly.zero();
        for (slong i = 0; i <= n; i++) {
            ulong coeff = nullspace_basis.get_entry(i, 0);
            if (coeff != 0) {
                min_poly.set_coeff(i, coeff);
            }
        }
        
        // Make monic if possible
        slong deg = min_poly.degree();
        if (deg >= 0) {
            ulong leading = min_poly.get_coeff(deg);
            if (leading != 0) {
                ulong inv_leading = n_invmod(leading, prime_);
                for (slong i = 0; i <= deg; i++) {
                    ulong old_coeff = min_poly.get_coeff(i);
                    ulong new_coeff = n_mulmod2_preinv(old_coeff, inv_leading, prime_, 0);
                    min_poly.set_coeff(i, new_coeff);
                }
            }
        }
    }
    
    ulong prime() const { return prime_; }
};

} // namespace flint_linalg