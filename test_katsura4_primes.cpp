#include <iostream>
#include <vector>
#include "src/axf4.h"
#include "src/julia_rur/f4_monomial_decoder.hpp"
#include "src/julia_rur/quotient_basis.hpp"

using namespace julia_rur;

// Test Katsura-4 with a specific prime
void test_katsura4_prime(uint32_t prime) {
    std::cout << "\n=== Testing Katsura-4 with prime " << prime 
              << " (approx 2^" << (int)(log2(prime)) << ") ===" << std::endl;
    
    // Initialize F4 session
    axf4_session_t session = axf4_create_session(5, prime);
    
    // Katsura-4 polynomials
    // x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1
    std::vector<uint32_t> coeff1 = {1, 2, 2, 2, 2, prime - 1};
    std::vector<uint32_t> mono1 = {16, 8, 4, 2, 1, 0};
    axf4_add_polynomial(session, 6, coeff1.data(), mono1.data());
    
    // x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1
    std::vector<uint32_t> coeff2 = {1, 2, 2, 2, prime - 1};
    std::vector<uint32_t> mono2 = {24, 12, 6, 3, 8};
    axf4_add_polynomial(session, 5, coeff2.data(), mono2.data());
    
    // x0*x2 + 2*x1*x3 + 2*x2*x4 - x2
    std::vector<uint32_t> coeff3 = {1, 2, 2, prime - 1};
    std::vector<uint32_t> mono3 = {20, 10, 5, 4};
    axf4_add_polynomial(session, 4, coeff3.data(), mono3.data());
    
    // x0*x3 + 2*x1*x4 - x3
    std::vector<uint32_t> coeff4 = {1, 2, prime - 1};
    std::vector<uint32_t> mono4 = {18, 9, 2};
    axf4_add_polynomial(session, 3, coeff4.data(), mono4.data());
    
    // x0*x4 - x4
    std::vector<uint32_t> coeff5 = {1, prime - 1};
    std::vector<uint32_t> mono5 = {17, 1};
    axf4_add_polynomial(session, 2, coeff5.data(), mono5.data());
    
    // Compute Gröbner basis
    axf4_result_t result = axf4_compute_groebner_basis_keep_data(session);
    
    if (result.status != 0) {
        std::cout << "  F4 failed with status " << result.status << std::endl;
        axf4_free_result(&result);
        axf4_destroy_session(session);
        return;
    }
    
    // Get basis size
    int basis_size = axf4_get_basis_size(session);
    std::cout << "  GB size: " << basis_size << " polynomials" << std::endl;
    
    // Extract leading terms
    std::vector<PP> leading_terms;
    for (int i = 0; i < basis_size; ++i) {
        int term_count = axf4_get_poly_term_count(i);
        if (term_count > 0) {
            std::vector<uint32_t> coeffs(term_count);
            std::vector<uint32_t> monomials(term_count);
            axf4_get_poly_data(i, coeffs.data(), monomials.data());
            
            // Find leading term (highest in degrevlex)
            PP lt = decode_f4_monomial(monomials[0], 5);
            for (int j = 1; j < term_count; ++j) {
                PP pp = decode_f4_monomial(monomials[j], 5);
                if (degrevlex_greater(pp, lt)) {
                    lt = pp;
                }
            }
            leading_terms.push_back(lt);
            
            // Check for single-term polynomials
            if (term_count == 1) {
                int degree = 0;
                for (int d : lt) degree += d;
                if (degree >= 4) {
                    std::cout << "  WARNING: Polynomial " << i << " has only 1 term with degree " 
                              << degree << " : [";
                    for (size_t k = 0; k < lt.size(); ++k) {
                        if (k > 0) std::cout << ",";
                        std::cout << lt[k];
                    }
                    std::cout << "]" << std::endl;
                }
            }
        }
    }
    
    // Compute quotient basis
    try {
        auto qb = compute_quotient_basis(leading_terms);
        std::cout << "  Quotient basis dimension: " << qb.size() << std::endl;
        
        // Expected dimension for Katsura-4 is 12
        if (qb.size() == 12) {
            std::cout << "  ✓ GOOD PRIME: Correct dimension!" << std::endl;
        } else {
            std::cout << "  ✗ BAD PRIME: Expected dimension 12" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  Error computing quotient basis: " << e.what() << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
}

// Previous prime function
uint32_t prev_prime(uint32_t n) {
    if (n <= 2) return 2;
    if (n == 3) return 2;
    n = (n % 2 == 0) ? n - 1 : n - 2;
    
    while (n > 2) {
        bool is_prime = true;
        uint32_t sqrt_n = (uint32_t)sqrt(n);
        for (uint32_t i = 3; i <= sqrt_n; i += 2) {
            if (n % i == 0) {
                is_prime = false;
                break;
            }
        }
        if (is_prime) return n;
        n -= 2;
    }
    return 2;
}

int main() {
    std::cout << "Testing Katsura-4 with different primes\n";
    std::cout << "========================================\n";
    
    // Test with Julia's default range (2^28)
    std::cout << "\nJulia's default range (2^28):" << std::endl;
    uint32_t p28 = prev_prime((1U << 28) - 1);
    for (int i = 0; i < 3; ++i) {
        test_katsura4_prime(p28);
        p28 = prev_prime(p28 - 1);
    }
    
    // Test with our range (2^30)
    std::cout << "\n\nOur current range (2^30):" << std::endl;
    uint32_t p30 = prev_prime((1U << 30) - 1);
    for (int i = 0; i < 3; ++i) {
        test_katsura4_prime(p30);
        p30 = prev_prime(p30 - 1);
    }
    
    // Test with smaller primes
    std::cout << "\n\nSmaller range (2^20):" << std::endl;
    uint32_t p20 = prev_prime((1U << 20) - 1);
    for (int i = 0; i < 3; ++i) {
        test_katsura4_prime(p20);
        p20 = prev_prime(p20 - 1);
    }
    
    // Test with even smaller primes
    std::cout << "\n\nEven smaller range (2^16):" << std::endl;
    uint32_t p16 = prev_prime((1U << 16) - 1);
    for (int i = 0; i < 3; ++i) {
        test_katsura4_prime(p16);
        p16 = prev_prime(p16 - 1);
    }
    
    return 0;
}