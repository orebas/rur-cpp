#pragma once

#include <flint/flint.h>
#include <flint/nmod_mpoly.h>
#include <memory>
#include <stdexcept>
#include <vector>

/**
 * RAII wrappers for FLINT multivariate polynomial types
 * Provides safe memory management and clean C++ integration
 */

namespace flint_mpoly {

/**
 * RAII wrapper for nmod_mpoly_ctx_t (multivariate polynomial context)
 */
class MpolyContext {
private:
    nmod_mpoly_ctx_t ctx_;
    ulong prime_;
    
public:
    explicit MpolyContext(slong nvars, ulong prime, ordering_t ordering = ORD_DEGREVLEX) 
        : prime_(prime) {
        if (nvars <= 0) {
            throw std::invalid_argument("Number of variables must be positive");
        }
        if (prime < 2) {
            throw std::invalid_argument("Prime must be >= 2");
        }
        
        nmod_mpoly_ctx_init(ctx_, nvars, ordering, prime);
    }
    
    ~MpolyContext() {
        nmod_mpoly_ctx_clear(ctx_);
    }
    
    // Non-copyable
    MpolyContext(const MpolyContext&) = delete;
    MpolyContext& operator=(const MpolyContext&) = delete;
    
    // Non-movable (context should be created once and shared)
    MpolyContext(MpolyContext&&) = delete;
    MpolyContext& operator=(MpolyContext&&) = delete;
    
    const nmod_mpoly_ctx_t& get() const { return ctx_; }
    nmod_mpoly_ctx_t& get() { return ctx_; }
    
    ulong prime() const { return prime_; }
    slong nvars() const { return ctx_->minfo->nvars; }
    ordering_t ordering() const { return ctx_->minfo->ord; }
};

/**
 * RAII wrapper for nmod_mpoly_t (multivariate polynomial)
 */
class Mpoly {
private:
    nmod_mpoly_t poly_;
    std::shared_ptr<const MpolyContext> ctx_;
    
public:
    explicit Mpoly(std::shared_ptr<const MpolyContext> ctx) : ctx_(ctx) {
        if (!ctx) {
            throw std::invalid_argument("Context cannot be null");
        }
        nmod_mpoly_init(poly_, ctx_->get());
    }
    
    ~Mpoly() {
        nmod_mpoly_clear(poly_, ctx_->get());
    }
    
    // Non-copyable for now (can add copy constructor if needed)
    Mpoly(const Mpoly&) = delete;
    Mpoly& operator=(const Mpoly&) = delete;
    
    // Movable
    Mpoly(Mpoly&& other) noexcept : ctx_(other.ctx_) {
        // Steal the polynomial data
        *poly_ = *other.poly_;
        // Re-initialize the moved-from object to valid empty state
        nmod_mpoly_init(other.poly_, ctx_->get());
    }
    
    Mpoly& operator=(Mpoly&& other) {
        if (this != &other) {
            if (ctx_ != other.ctx_) {
                throw std::runtime_error("Cannot move polynomials with different contexts");
            }
            
            // Clear current state
            nmod_mpoly_clear(poly_, ctx_->get());
            
            // Steal the data
            *poly_ = *other.poly_;
            
            // Re-initialize moved-from object
            nmod_mpoly_init(other.poly_, ctx_->get());
        }
        return *this;
    }
    
    nmod_mpoly_t& get() { return poly_; }
    const nmod_mpoly_t& get() const { return poly_; }
    
    std::shared_ptr<const MpolyContext> context() const { return ctx_; }
    
    // Convenience methods
    bool is_zero() const {
        return nmod_mpoly_is_zero(poly_, ctx_->get());
    }
    
    void zero() {
        nmod_mpoly_zero(poly_, ctx_->get());
    }
    
    slong length() const {
        return nmod_mpoly_length(poly_, ctx_->get());
    }
    
    // Set coefficient for a monomial (exponents must be array of size nvars)
    void set_coeff_ui_monomial(ulong coeff, const ulong* exps) {
        nmod_mpoly_set_coeff_ui_ui(poly_, coeff, exps, ctx_->get());
    }
    
    // Get coefficient for a monomial
    ulong get_coeff_ui_monomial(const ulong* exps) const {
        return nmod_mpoly_get_coeff_ui_ui(poly_, exps, ctx_->get());
    }
    
    // Iterator support for accessing terms
    class TermIterator {
    private:
        const Mpoly* poly_;
        slong index_;
        
    public:
        TermIterator(const Mpoly* poly, slong index) 
            : poly_(poly), index_(index) {}
        
        bool operator!=(const TermIterator& other) const {
            return index_ != other.index_;
        }
        
        TermIterator& operator++() {
            ++index_;
            return *this;
        }
        
        struct Term {
            ulong coeff;
            std::vector<ulong> exps;
        };
        
        Term operator*() const {
            Term term;
            term.coeff = nmod_mpoly_get_term_coeff_ui(poly_->poly_, index_, poly_->ctx_->get());
            
            term.exps.resize(poly_->ctx_->nvars());
            ulong* exp_vec = term.exps.data();
            nmod_mpoly_get_term_exp_ui(exp_vec, poly_->poly_, index_, poly_->ctx_->get());
            
            return term;
        }
    };
    
    TermIterator begin() const {
        return TermIterator(this, 0);
    }
    
    TermIterator end() const {
        return TermIterator(this, length());
    }
};

/**
 * Helper functions for polynomial operations
 */

// Perform ideal division: compute Q and R such that A = sum(Q[i] * B[i]) + R
void divrem_ideal(std::vector<Mpoly>& quotients, Mpoly& remainder,
                  const Mpoly& dividend, const std::vector<Mpoly>& divisors);

} // namespace flint_mpoly