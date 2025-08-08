#pragma once

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <stdexcept>
#include <vector>

/**
 * RAII wrappers for FLINT types to ensure safe memory management
 * and clean integration with C++ code
 */

namespace flint_cpp {

/**
 * RAII wrapper for nmod_t (modular arithmetic context)
 */
class ModContext {
private:
    nmod_t mod_;
    
public:
    explicit ModContext(ulong prime) {
        nmod_init(&mod_, prime);
    }
    
    ~ModContext() {
        // nmod_t doesn't need explicit cleanup
    }
    
    // Non-copyable
    ModContext(const ModContext&) = delete;
    ModContext& operator=(const ModContext&) = delete;
    
    // Movable
    ModContext(ModContext&& other) noexcept {
        mod_ = other.mod_;
        // Other doesn't need cleanup for nmod_t
    }
    
    ModContext& operator=(ModContext&& other) noexcept {
        if (this != &other) {
            mod_ = other.mod_;
        }
        return *this;
    }
    
    const nmod_t& get() const { return mod_; }
    nmod_t& get() { return mod_; }
    ulong prime() const { return mod_.n; }
};

/**
 * RAII wrapper for nmod_poly_t (polynomial mod prime)
 */
class ModPoly {
private:
    nmod_poly_t poly_;
    const ModContext* ctx_;
    
public:
    explicit ModPoly(const ModContext& ctx) : ctx_(&ctx) {
        nmod_poly_init(poly_, ctx.prime());
    }
    
    ~ModPoly() {
        nmod_poly_clear(poly_);
    }
    
    // Non-copyable for now (can add if needed)
    ModPoly(const ModPoly&) = delete;
    ModPoly& operator=(const ModPoly&) = delete;
    
    // Movable
    ModPoly(ModPoly&& other) noexcept : ctx_(other.ctx_) {
        *poly_ = *other.poly_;
        nmod_poly_init(other.poly_, ctx_->prime()); // Reset other
    }
    
    ModPoly& operator=(ModPoly&& other) noexcept {
        if (this != &other) {
            nmod_poly_clear(poly_);
            ctx_ = other.ctx_;
            *poly_ = *other.poly_;
            nmod_poly_init(other.poly_, ctx_->prime()); // Reset other
        }
        return *this;
    }
    
    const nmod_poly_t& get() const { return poly_; }
    nmod_poly_t& get() { return poly_; }
    
    // Set coefficient at given index
    void set_coeff_ui(slong i, ulong c) {
        nmod_poly_set_coeff_ui(poly_, i, c);
    }
    
    // Get coefficient at given index
    ulong get_coeff_ui(slong i) const {
        return nmod_poly_get_coeff_ui(poly_, i);
    }
    
    // Get degree
    slong degree() const {
        return nmod_poly_degree(poly_);
    }
    
    // Set to zero
    void zero() {
        nmod_poly_zero(poly_);
    }
    
    // Polynomial division: this = this % divisor
    void rem(const ModPoly& divisor) {
        nmod_poly_rem(poly_, poly_, divisor.get());
    }
    
    // Polynomial division: compute quotient and remainder
    void divrem(ModPoly& quotient, const ModPoly& divisor) {
        nmod_poly_divrem(quotient.get(), poly_, poly_, divisor.get());
    }
};

/**
 * RAII wrapper for nmod_mat_t (matrix mod prime)
 */
class ModMatrix {
private:
    nmod_mat_t mat_;
    const ModContext* ctx_;
    
public:
    ModMatrix(slong rows, slong cols, const ModContext& ctx) : ctx_(&ctx) {
        nmod_mat_init(mat_, rows, cols, ctx.prime());
    }
    
    ~ModMatrix() {
        nmod_mat_clear(mat_);
    }
    
    // Non-copyable for now
    ModMatrix(const ModMatrix&) = delete;
    ModMatrix& operator=(const ModMatrix&) = delete;
    
    // Movable
    ModMatrix(ModMatrix&& other) noexcept : ctx_(other.ctx_) {
        *mat_ = *other.mat_;
        nmod_mat_init(other.mat_, 0, 0, ctx_->prime()); // Reset other
    }
    
    ModMatrix& operator=(ModMatrix&& other) noexcept {
        if (this != &other) {
            nmod_mat_clear(mat_);
            ctx_ = other.ctx_;
            *mat_ = *other.mat_;
            nmod_mat_init(other.mat_, 0, 0, ctx_->prime()); // Reset other
        }
        return *this;
    }
    
    const nmod_mat_t& get() const { return mat_; }
    nmod_mat_t& get() { return mat_; }
    
    // Get dimensions
    slong nrows() const { return nmod_mat_nrows(mat_); }
    slong ncols() const { return nmod_mat_ncols(mat_); }
    
    // Set/get entry
    void set_entry_ui(slong i, slong j, ulong c) {
        nmod_mat_set_entry(mat_, i, j, c);
    }
    
    ulong get_entry_ui(slong i, slong j) const {
        return nmod_mat_get_entry(mat_, i, j);
    }
    
    // Zero the matrix
    void zero() {
        nmod_mat_zero(mat_);
    }
    
    // Reduced row echelon form - returns rank
    slong rref() {
        return nmod_mat_rref(mat_);
    }
    
    // Create copy in RREF form
    slong rref_copy(ModMatrix& result) const {
        if (result.nrows() != nrows() || result.ncols() != ncols()) {
            throw std::invalid_argument("Matrix dimensions must match for RREF copy");
        }
        nmod_mat_set(result.get(), mat_);
        return nmod_mat_rref(result.get());
    }
};

/**
 * RAII wrapper for fmpq_t (rational number)
 */
class Rational {
private:
    fmpq_t rat_;
    
public:
    Rational() {
        fmpq_init(rat_);
    }
    
    explicit Rational(slong num, slong den = 1) {
        fmpq_init(rat_);
        fmpq_set_si(rat_, num, den);
    }
    
    ~Rational() {
        fmpq_clear(rat_);
    }
    
    // Non-copyable for now
    Rational(const Rational&) = delete;
    Rational& operator=(const Rational&) = delete;
    
    // Movable
    Rational(Rational&& other) noexcept {
        *rat_ = *other.rat_;
        fmpq_init(other.rat_); // Reset other
    }
    
    Rational& operator=(Rational&& other) noexcept {
        if (this != &other) {
            fmpq_clear(rat_);
            *rat_ = *other.rat_;
            fmpq_init(other.rat_); // Reset other
        }
        return *this;
    }
    
    const fmpq_t& get() const { return rat_; }
    fmpq_t& get() { return rat_; }
    
    // Set from integers
    void set_si(slong num, slong den = 1) {
        fmpq_set_si(rat_, num, den);
    }
    
    // Canonicalize (reduce to lowest terms)
    void canonicalise() {
        fmpq_canonicalise(rat_);
    }
};

/**
 * RAII wrapper for fmpq_poly_t (polynomial with rational coefficients)
 */
class RationalPoly {
private:
    fmpq_poly_t poly_;
    
public:
    RationalPoly() {
        fmpq_poly_init(poly_);
    }
    
    ~RationalPoly() {
        fmpq_poly_clear(poly_);
    }
    
    // Non-copyable for now
    RationalPoly(const RationalPoly&) = delete;
    RationalPoly& operator=(const RationalPoly&) = delete;
    
    // Movable
    RationalPoly(RationalPoly&& other) noexcept {
        *poly_ = *other.poly_;
        fmpq_poly_init(other.poly_); // Reset other
    }
    
    RationalPoly& operator=(RationalPoly&& other) noexcept {
        if (this != &other) {
            fmpq_poly_clear(poly_);
            *poly_ = *other.poly_;
            fmpq_poly_init(other.poly_); // Reset other
        }
        return *this;
    }
    
    const fmpq_poly_t& get() const { return poly_; }
    fmpq_poly_t& get() { return poly_; }
    
    // Set coefficient
    void set_coeff_si(slong i, slong num, slong den = 1) {
        fmpq_poly_set_coeff_si(poly_, i, num);
        if (den != 1) {
            fmpq_t coeff;
            fmpq_init(coeff);
            fmpq_set_si(coeff, num, den);
            fmpq_poly_set_coeff_fmpq(poly_, i, coeff);
            fmpq_clear(coeff);
        }
    }
    
    // Get degree
    slong degree() const {
        return fmpq_poly_degree(poly_);
    }
    
    // Set to zero
    void zero() {
        fmpq_poly_zero(poly_);
    }
    
    // Addition
    void add(const RationalPoly& other) {
        fmpq_poly_add(poly_, poly_, other.get());
    }
    
    // Multiplication
    void mul(const RationalPoly& other) {
        fmpq_poly_mul(poly_, poly_, other.get());
    }
};

} // namespace flint_cpp