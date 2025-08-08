#include "flint_mpoly_wrappers.hpp"

namespace flint_mpoly {

// Perform ideal division: compute Q and R such that A = sum(Q[i] * B[i]) + R
void divrem_ideal(std::vector<Mpoly>& quotients, Mpoly& remainder,
                  const Mpoly& dividend, const std::vector<Mpoly>& divisors) {
    if (quotients.size() != divisors.size()) {
        throw std::invalid_argument("Quotients and divisors must have same size");
    }
    
    auto ctx = dividend.context();
    
    // Prepare C-style arrays for FLINT
    std::vector<nmod_mpoly_struct*> q_ptrs;
    std::vector<nmod_mpoly_struct*> b_ptrs;
    
    q_ptrs.reserve(quotients.size());
    b_ptrs.reserve(divisors.size());
    
    for (auto& q : quotients) {
        if (q.context() != ctx) {
            throw std::runtime_error("All polynomials must have same context");
        }
        q_ptrs.push_back(q.get());
    }
    
    for (const auto& b : divisors) {
        if (b.context() != ctx) {
            throw std::runtime_error("All polynomials must have same context");
        }
        b_ptrs.push_back(const_cast<nmod_mpoly_struct*>(b.get()));
    }
    
    if (remainder.context() != ctx) {
        throw std::runtime_error("All polynomials must have same context");
    }
    
    // Call FLINT
    nmod_mpoly_divrem_ideal(q_ptrs.data(), remainder.get(),
                           dividend.get(), 
                           const_cast<nmod_mpoly_struct* const*>(b_ptrs.data()),
                           divisors.size(), ctx->get());
}

} // namespace flint_mpoly