#include "multi_modular.hpp"
#include "polynomial_operations.hpp" // For modular_inverse
// #include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <stdexcept>


namespace julia_rur {

// diagnostics removed

bool
is_prime(ModularCoeff n) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;

    // Check divisibility up to sqrt(n)
    for (ModularCoeff i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) { return false; }
    }
    return true;
}

// No longer needed - F4 bug has been fixed

ModularCoeff
prev_prime(ModularCoeff n) {
    // Find the largest prime less than n
    if (n <= 2) { throw std::invalid_argument("No prime less than 2"); }

    n--; // Start checking from n-1
    while (n > 1) {
        if (is_prime(n)) { return n; }
        n--;
    }

    if (n <= 1) { throw std::invalid_argument("No prime found"); }

    return n;
}

std::vector<ModularCoeff>
prev_primes(ModularCoeff start, size_t count) {
    std::vector<ModularCoeff> primes;
    primes.reserve(count);

    // Debug print
    if (start > 1000000000u) {
        std::cerr << "WARNING: prev_primes called with very large start value: " << start << std::endl;
    }

    ModularCoeff current = start;
    for (size_t i = 0; i < count; ++i) {
        try {
            current = prev_prime(current);
            primes.push_back(current);
        } catch (const std::exception &e) {
            // No more primes available
            break;
        }
    }

    return primes;
}

ExtendedGCDResult
extended_gcd(const mpz_class &a, const mpz_class &b) {
    ExtendedGCDResult result;

    if (b == 0) {
        result.gcd = a;
        result.x = 1;
        result.y = 0;
        return result;
    }

    mpz_class old_r = a, r = b;
    mpz_class old_s = 1, s = 0;
    mpz_class old_t = 0, t = 1;

    while (r != 0) {
        mpz_class quotient = old_r / r;

        mpz_class temp = r;
        r = old_r - quotient * r;
        old_r = temp;

        temp = s;
        s = old_s - quotient * s;
        old_s = temp;

        temp = t;
        t = old_t - quotient * t;
        old_t = temp;
    }

    result.gcd = old_r;
    result.x = old_s;
    result.y = old_t;

    return result;
}

mpz_class
chinese_remainder_theorem(mpz_class &result,
                          const std::vector<ModularCoeff> &remainders,
                          const std::vector<ModularCoeff> &primes,
                          std::vector<mpz_class> &multipliers) {
    if (remainders.size() != primes.size()) {
        throw std::invalid_argument("Remainders and primes must have same size");
    }

    size_t n = primes.size();

    // Compute total modulus M = product of all primes
    mpz_class M = 1;
    for (ModularCoeff p : primes) { M *= p; }

    // Precompute multipliers if not already done
    if (multipliers.size() != n) {
        multipliers.clear();
        multipliers.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            mpz_class Mi = M / primes[i];

            // Find Mi_inv such that Mi * Mi_inv ≡ 1 (mod primes[i])
            ExtendedGCDResult egcd = extended_gcd(Mi, mpz_class(primes[i]));
            mpz_class Mi_inv = egcd.x;

            // Ensure Mi_inv is positive
            while (Mi_inv < 0) { Mi_inv += primes[i]; }

            multipliers.push_back(Mi * Mi_inv);
        }
    }

    // Apply CRT formula: result = sum(remainders[i] * multipliers[i]) mod M
    result = 0;
    for (size_t i = 0; i < n; ++i) { result += mpz_class(remainders[i]) * multipliers[i]; }
    result %= M;

    // Ensure result is positive
    if (result < 0) { result += M; }

    return M;
}

RationalReconstructionResult
rational_reconstruction(const mpz_class &a, const mpz_class &m, const mpz_class &N_input, const mpz_class &D_input) {
    // Default bounds if not provided
    mpz_class N = N_input;
    mpz_class D = D_input;

    if (N == 0 || D == 0) {
        // Use balanced bounds: N = D = floor(sqrt(m/2))
        mpz_class temp, m_half = m / 2;
        mpz_sqrt(temp.get_mpz_t(), m_half.get_mpz_t());
        if (N == 0) N = temp;
        if (D == 0) D = temp;
    }

    // Special case: a = 0
    if (a == 0) {
        // Reconstruction for 0 is 0/1. Check if it fits the bounds.
        // The numerator 0 always satisfies |0| <= N.
        // The denominator 1 must satisfy |1| <= D.
        if (D >= 1) { return RationalReconstructionResult(true, mpq_class(0, 1)); }
        return RationalReconstructionResult(false, mpq_class(0));
    }

    // Extended Euclidean algorithm with early termination
    mpz_class r0 = m, r1 = a % m;
    mpz_class s0 = 1, s1 = 0;
    mpz_class t0 = 0, t1 = 1;

    while (r1 != 0 && r0 > N) {
        mpz_class q = r0 / r1;

        mpz_class temp = r1;
        r1 = r0 - q * r1;
        r0 = temp;

        temp = s1;
        s1 = s0 - q * s1;
        s0 = temp;

        temp = t1;
        t1 = t0 - q * t1;
        t0 = temp;
    }

    // Check if we found a valid rational reconstruction
    if (abs(t0) <= D && abs(r0) <= N && t0 != 0) {
        // Normalize to positive denominator
        if (t0 < 0) {
            r0 = -r0;
            t0 = -t0;
        }

        // Verify: a * t0 ≡ r0 (mod m)
        mpz_class verification = (a * t0 - r0) % m;
        if (verification < 0) verification += m; // Normalize to [0, m-1]
        if (verification == 0) {
            mpq_class result(r0, t0);
            result.canonicalize();
            return RationalReconstructionResult(true, result);
        }
    }

    return RationalReconstructionResult(false, mpq_class(0));
}

RationalReconstructionResult
rational_reconstruction_with_denominator(const mpz_class &zz,
                                         const mpz_class &den,
                                         const mpz_class &modulus,
                                         const mpz_class &N,
                                         const mpz_class &D,
                                         ModularCoeff zp,
                                         ModularCoeff p) {
    // Compute rem_nemo = (den * zz) mod modulus
    mpz_class rem_nemo = (den * zz) % modulus;
    if (rem_nemo < 0) rem_nemo += modulus;

    // Try rational reconstruction
    RationalReconstructionResult rat_result = rational_reconstruction(rem_nemo, modulus, N, D);

    if (!rat_result.success) { return RationalReconstructionResult(false, mpq_class(0)); }

    // Compute PQ = pq / den
    mpq_class PQ = rat_result.rational / den;
    PQ.canonicalize();

    // Verify: zp == (numerator(PQ) * invmod(denominator(PQ), p)) mod p
    mpz_class num = PQ.get_num();
    mpz_class denom = PQ.get_den();

    // Compute modular inverse of denominator
    mpz_class denom_temp = denom % p;
    ModularCoeff denom_mod = static_cast<ModularCoeff>(denom_temp.get_ui());
    ModularCoeff inv_denom = modular_inverse(denom_mod, p);

    // Compute verification value
    mpz_class num_temp = num % p;
    ModularCoeff num_mod = static_cast<ModularCoeff>(num_temp.get_ui());
    if (num < 0) {
        mpz_class neg_num_temp = (-num) % p;
        num_mod = p - static_cast<ModularCoeff>(neg_num_temp.get_ui());
    }

    ModularCoeff verification = static_cast<ModularCoeff>((static_cast<AccModularCoeff>(num_mod) * inv_denom) % p);

    if (verification != zp) { return RationalReconstructionResult(false, mpq_class(0)); }

    return RationalReconstructionResult(true, PQ);
}

RationalReconstructionResult
rational_reconstruction_with_denominator_no_verify(const mpz_class &zz,
                                                   const mpz_class &den,
                                                   const mpz_class &modulus,
                                                   const mpz_class &N,
                                                   const mpz_class &D) {
    // Compute rem_nemo = (den * zz) mod modulus
    mpz_class rem_nemo = (den * zz) % modulus;
    if (rem_nemo < 0) rem_nemo += modulus;

    // Try rational reconstruction
    RationalReconstructionResult rat_result = rational_reconstruction(rem_nemo, modulus, N, D);

    if (!rat_result.success) { return RationalReconstructionResult(false, mpq_class(0)); }

    // Compute PQ = pq / den
    mpq_class PQ = rat_result.rational / den;
    PQ.canonicalize();

    return RationalReconstructionResult(true, PQ);
}

RationalReconstructionBounds
compute_balanced_bounds(const mpz_class &modulus) {
    RationalReconstructionBounds bounds;
    mpz_class modulus_half = modulus / 2;
    mpz_sqrt(bounds.N.get_mpz_t(), modulus_half.get_mpz_t());
    bounds.D = bounds.N;
    return bounds;
}

RationalReconstructionBounds
compute_unbalanced_bounds_9to1(const mpz_class &modulus) {
    RationalReconstructionBounds bounds;

    // N gets 90% of the bits, D gets 10%
    // Total bits available: log2(modulus/2)
    mpz_class temp = modulus / 2;
    size_t total_bits = mpz_sizeinbase(temp.get_mpz_t(), 2);

    size_t n_bits = (total_bits * 9) / 10;
    size_t d_bits = total_bits - n_bits;

    bounds.N = mpz_class(1) << n_bits;
    bounds.D = mpz_class(1) << d_bits;

    return bounds;
}

RationalReconstructionBounds
compute_unbalanced_bounds_99to1(const mpz_class &modulus) {
    RationalReconstructionBounds bounds;

    // N gets 99% of the bits, D gets 1%
    mpz_class temp = modulus / 2;
    size_t total_bits = mpz_sizeinbase(temp.get_mpz_t(), 2);

    size_t n_bits = (total_bits * 99) / 100;
    size_t d_bits = total_bits - n_bits;

    bounds.N = mpz_class(1) << n_bits;
    bounds.D = mpz_class(1) << d_bits;

    return bounds;
}

RationalReconstructionBounds
compute_constant_bounds(const mpz_class &modulus) {
    RationalReconstructionBounds bounds;

    // Very unbalanced: D = 1 (denominators must be 1)
    bounds.N = modulus / 2;
    bounds.D = 1;

    return bounds;
}

std::pair<bool, mpz_class>
crt_and_rational_reconstruction(std::vector<std::vector<mpq_class>> &qq_result,
                                std::vector<std::vector<mpz_class>> &zz_temp,
                                const std::vector<std::vector<std::vector<ModularCoeff>>> &modular_tables,
                                const std::vector<ModularCoeff> &used_primes,
                                const std::vector<std::vector<ModularCoeff>> &reference_mod_p,
                                ModularCoeff verification_prime,
                                mpz_class &known_denominator) {
    if (used_primes.empty()) { return { false, 1 }; }

    const bool verbose_crt = false;
    if (verbose_crt) {
        std::cout << "CRT DEBUG: Starting reconstruction with " << used_primes.size() << " primes" << std::endl;
    }
    if (std::getenv("RUR_VERBOSE_CRT")) {
        std::cout << "  Reference dimensions: " << reference_mod_p.size() << " polynomials" << std::endl;
    }
    if (std::getenv("RUR_VERBOSE_CRT")) {
        for (size_t i = 0; i < std::min(reference_mod_p.size(), size_t(3)); ++i) {
            std::cout << "    Poly " << i << ": " << reference_mod_p[i].size() << " coefficients" << std::endl;
        }
    }

    if (qq_result.empty()) {
        qq_result.resize(reference_mod_p.size());
        zz_temp.resize(reference_mod_p.size());
        for (size_t i = 0; i < reference_mod_p.size(); ++i) {
            qq_result[i].resize(reference_mod_p[i].size());
            zz_temp[i].resize(reference_mod_p[i].size());
        }
    }

    mpz_class modulus = 1;
    for (ModularCoeff p : used_primes) { modulus *= p; }

    auto bounds_balanced = compute_balanced_bounds(modulus);
    auto bounds_9to1 = compute_unbalanced_bounds_9to1(modulus);
    auto bounds_99to1 = compute_unbalanced_bounds_99to1(modulus);
    auto bounds_const = compute_constant_bounds(modulus);

    bool all_success = true;
    // diagnostics disabled
    for (size_t i = 0; i < qq_result.size(); ++i) {
        for (int j = qq_result[i].size() - 1; j >= 0; --j) {
            std::vector<ModularCoeff> rems;
            for (size_t k = 0; k < used_primes.size(); ++k) {
                if (i < modular_tables[k].size() && j < modular_tables[k][i].size()) {
                    rems.push_back(modular_tables[k][i][j]);
                } else {
                    rems.push_back(0);
                }
            }

            // Special debugging for the failing coefficient
            if (i == 0 && j == 6) {
                if (verbose_crt) { std::cout << "CRT DEBUG: Coefficient [0][6] remainders (first 10):" << std::endl; }
                for (size_t k = 0; k < std::min(used_primes.size(), size_t(10)); ++k) {
                    if (std::getenv("RUR_VERBOSE_CRT")) {
                        std::cout << "  Prime " << used_primes[k] << ": remainder " << rems[k] << std::endl;
                    }
                }
                if (std::getenv("RUR_VERBOSE_CRT")) {
                    std::cout << "  ... (showing first 10 of " << used_primes.size() << " primes)" << std::endl;
                }
            }

            mpz_class crt_result;
            std::vector<mpz_class> mults; // compute fresh to match original behavior
            chinese_remainder_theorem(crt_result, rems, used_primes, mults);
            zz_temp[i][j] = crt_result;

            ModularCoeff zp_verify = reference_mod_p[i][j];
            // diagnostics disabled

            // Original order: try 9:1 first, then const, then balanced, then reuse next-den with 99:1
            auto res1 = rational_reconstruction_with_denominator(
              crt_result, known_denominator, modulus, bounds_9to1.N, bounds_9to1.D, zp_verify, verification_prime);
            if (res1.success) {
                qq_result[i][j] = res1.rational;
            } else {
                auto res2 = rational_reconstruction_with_denominator(crt_result,
                                                                     known_denominator,
                                                                     modulus,
                                                                     bounds_const.N,
                                                                     bounds_const.D,
                                                                     zp_verify,
                                                                     verification_prime);
                if (res2.success) {
                    qq_result[i][j] = res2.rational;
                } else {
                    auto res3 = rational_reconstruction_with_denominator(
                      crt_result, 1, modulus, bounds_balanced.N, bounds_balanced.D, zp_verify, verification_prime);
                    if (res3.success) {
                        qq_result[i][j] = res3.rational;
                    } else if (j + 1 < static_cast<int>(qq_result[i].size())) {
                        auto res4 = rational_reconstruction_with_denominator(crt_result,
                                                                             qq_result[i][j + 1].get_den(),
                                                                             modulus,
                                                                             bounds_99to1.N,
                                                                             bounds_99to1.D,
                                                                             zp_verify,
                                                                             verification_prime);
                        if (res4.success) {
                            qq_result[i][j] = res4.rational;
                        } else {
                            // Final fallback even when reuse is available: try no-verify balanced bounds
                            auto res_nv = rational_reconstruction_with_denominator_no_verify(
                              crt_result, known_denominator, modulus, bounds_balanced.N, bounds_balanced.D);
                            if (res_nv.success) {
                                qq_result[i][j] = res_nv.rational;
                            } else {
                                std::cout << "CRT FAILURE: Rational reconstruction failed for coefficient [" << i
                                          << "][" << j << "]" << std::endl;
                                std::cout << "  CRT result: " << crt_result << std::endl;
                                std::cout << "  Modulus: " << modulus << std::endl;
                                std::cout << "  Verification: expected " << zp_verify << " mod " << verification_prime
                                          << std::endl;
                                std::cout << "  All 4 reconstruction methods failed" << std::endl;
                                all_success = false;
                                break;
                            }
                        }
                    } else {
                        // Final fallback: attempt no-verify rational reconstruction with balanced bounds
                        auto res_nv = rational_reconstruction_with_denominator_no_verify(
                          crt_result, known_denominator, modulus, bounds_balanced.N, bounds_balanced.D);
                        if (res_nv.success) {
                            qq_result[i][j] = res_nv.rational;
                        } else {
                            std::cout << "CRT FAILURE: Rational reconstruction failed for coefficient [" << i << "]["
                                      << j << "] (no previous coefficient available for reuse strategy)" << std::endl;
                            std::cout << "  CRT result: " << crt_result << std::endl;
                            std::cout << "  Modulus: " << modulus << std::endl;
                            std::cout << "  Verification: expected " << zp_verify << " mod " << verification_prime
                                      << std::endl;
                            std::cout << "  First 3 reconstruction methods failed" << std::endl;
                            all_success = false;
                            break;
                        }
                    }
                }
            }
            known_denominator = lcm(known_denominator, qq_result[i][j].get_den());
        }
        if (!all_success) break;
    }

    return { all_success, known_denominator };
}
std::vector<ModularCoeff>
modular_coeffs_vect(const std::vector<mpq_class> &rational_coeffs, ModularCoeff prime) {
    std::vector<ModularCoeff> result;
    result.reserve(rational_coeffs.size());

    for (const auto &coeff : rational_coeffs) {
        // Compute (numerator * inverse(denominator)) mod prime
        mpz_class num = coeff.get_num();
        mpz_class denom = coeff.get_den();

        // Handle negative numerators
        bool negative = (num < 0);
        if (negative) num = -num;

        // Reduce modulo prime
        mpz_class num_temp = num % prime;
        mpz_class denom_temp = denom % prime;
        ModularCoeff num_mod = static_cast<ModularCoeff>(num_temp.get_ui());
        ModularCoeff denom_mod = static_cast<ModularCoeff>(denom_temp.get_ui());

        if (denom_mod == 0) { throw std::runtime_error("Denominator is zero modulo prime"); }

        // Compute modular inverse of denominator
        ModularCoeff inv_denom = modular_inverse(denom_mod, prime);

        // Compute result
        ModularCoeff value = static_cast<ModularCoeff>((static_cast<AccModularCoeff>(num_mod) * inv_denom) % prime);

        // Handle negative case
        if (negative && value != 0) { value = prime - value; }

        result.push_back(value);
    }

    return result;
}

std::vector<ModularCoeff>
generate_primes(size_t bits, size_t count) {
    std::vector<ModularCoeff> primes;

    // F4 requires primes > 65536 (2^16) per axf4_wrapper.c
    if (bits < 17) { bits = 17; }

    // Simple prime generation for modular arithmetic
    // Start with some known primes suitable for F4
    if (bits == 17) {
        // Primes just above 2^16 (F4's minimum requirement)
        const ModularCoeff small_primes[] = { 65537, 65539, 65543, 65551, 65557, 65563, 65579, 65581, 65587, 65599 };

        for (size_t i = 0; i < count && i < 10; ++i) { primes.push_back(small_primes[i]); }
    } else if (bits <= 20) {
        // Primes around 2^20
        const ModularCoeff medium_primes[] = { 1048583, 1048589, 1048601, 1048607, 1048609,
                                               1048613, 1048633, 1048661, 1048673, 1048681 };

        for (size_t i = 0; i < count && i < 10; ++i) { primes.push_back(medium_primes[i]); }
    } else if (bits <= 28) {
        // Primes around 2^28 - matching Julia's default pr_max_bitsize = 28
        // First list primes where 2 is a quadratic residue (p ≡ 1,7 mod 8)
        // These are good for systems involving sqrt(2)
        const ModularCoeff primes_28bit_2qr[] = {
            268435399, // mod 8 = 7
            268435367, // mod 8 = 7
            268435361, // mod 8 = 1
            268435337, // mod 8 = 1
            268435313, // mod 8 = 1
            268435273, // mod 8 = 1
            268435183, // mod 8 = 7
            268435129, // mod 8 = 1
            268435121, // mod 8 = 1
            268435039, // mod 8 = 7
            268434979, // mod 8 = 3 (backup)
            268434961, // mod 8 = 1
            268434949, // mod 8 = 5 (backup)
            268434937, // mod 8 = 1
            268434873, // mod 8 = 1
            268434837, // mod 8 = 5 (backup)
            268434817, // mod 8 = 1
            268434793, // mod 8 = 1
            268434721, // mod 8 = 1
            268434697  // mod 8 = 1
        };

        for (size_t i = 0; i < count && i < 20; ++i) { primes.push_back(primes_28bit_2qr[i]); }
    } else if (bits <= 30) {
        // Primes around 2^30 - still safe for F4
        const ModularCoeff primes_30bit[] = { 1073741827, 1073741831, 1073741833, 1073741839, 1073741843,
                                              1073741857, 1073741873, 1073741909, 1073741939, 1073741941,
                                              1073741953, 1073741969, 1073741971, 1073741981, 1073741987,
                                              1073742013, 1073742047, 1073742053, 1073742071, 1073742091 };

        for (size_t i = 0; i < count && i < 20; ++i) { primes.push_back(primes_30bit[i]); }
    } else {
        // 31-bit primes - use with caution, only the first one is safe for F4
        const ModularCoeff large_primes[] = {
            2147483647, // Largest safe prime (2^31 - 1) for F4
            2147483629, // WARNING: These exceed F4's safe limit
            2147483587, // WARNING: These exceed F4's safe limit
            2147483579, // WARNING: These exceed F4's safe limit
            2147483563  // WARNING: These exceed F4's safe limit
        };

        // Only include the first prime which is safe for F4
        if (count > 0) { primes.push_back(large_primes[0]); }

        // Fill remaining with 30-bit primes which are safe
        const ModularCoeff fallback_primes[] = { 1073741783, 1073741741, 1073741723, 1073741719, 1073741717 };

        for (size_t i = 1; i < count && i < 6; ++i) { primes.push_back(fallback_primes[i - 1]); }
    }

    return primes;
}
void
crt_precompute(const mpz_class &modulus,
               mpz_class &n1,
               mpz_class &n2,
               std::vector<mpz_class> &mults,
               const std::vector<ModularCoeff> &primes) {
    mults.resize(primes.size());
    for (size_t i = 0; i < primes.size(); ++i) {
        mpz_class mi = primes[i];
        mpz_class ni = modulus / mi;
        mpz_class inv;
        mpz_invert(inv.get_mpz_t(), ni.get_mpz_t(), mi.get_mpz_t());
        mults[i] = inv * ni;
    }
}
} // namespace julia_rur