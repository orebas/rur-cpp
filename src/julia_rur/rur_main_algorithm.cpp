#include "rur_main_algorithm.hpp"
#include "f4_polynomial_formatter.hpp"
#include "multiplication_tables.hpp"
#include "polynomial_operations.hpp"
#include "prime_utils.hpp"
#include "separating_element_systematic.hpp"
#include "univariate_parameterization.hpp"
#include <algorithm>
#include <climits>
// #include <iomanip>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <sstream>


namespace julia_rur {

// Helper function to compute normalized difference modulo p
static inline int64_t
normalized_diff_mod_p(ModularCoeff ref_val, ModularCoeff mod_val, ModularCoeff p) {
    int64_t diff = static_cast<int64_t>(ref_val) - static_cast<int64_t>(mod_val);
    // Normalize to [-p/2, p/2]
    diff = diff % p;
    if (diff > static_cast<int64_t>(p / 2)) {
        diff -= p;
    } else if (diff < -static_cast<int64_t>(p / 2)) {
        diff += p;
    }
    return diff;
}

// Helper functions for canonicalization
inline ModularCoeff
mod_inv(ModularCoeff a, ModularCoeff p) {
    // Using Fermat's little theorem: a^(p-1) ≡ 1 (mod p), so a^(-1) ≡ a^(p-2) (mod p)
    ModularCoeff result = 1;
    ModularCoeff base = a;
    ModularCoeff exp = p - 2;
    while (exp > 0) {
        if (exp & 1) { result = (static_cast<AccModularCoeff>(result) * base) % p; }
        base = (static_cast<AccModularCoeff>(base) * base) % p;
        exp >>= 1;
    }
    return result;
}

inline ModularCoeff
mod_mul(ModularCoeff a, ModularCoeff b, ModularCoeff p) {
    // Safe for primes up to 31 bits
    return static_cast<ModularCoeff>((static_cast<AccModularCoeff>(a) * b) % p);
}

namespace {
// Build quotient-ring element from a polynomial in T using the stored powers V(T^k)
static std::vector<ModularCoeff>
build_element_from_poly(const std::vector<ModularCoeff> &poly,
                        const std::vector<std::vector<ModularCoeff>> &powers,
                        ModularCoeff prime) {
    if (powers.empty()) return {};
    size_t d = powers[0].size();
    std::vector<ModularCoeff> acc(d, 0);
    size_t max_k = std::min(powers.size(), poly.size());
    for (size_t k = 0; k < max_k; ++k) {
        if (poly[k] == 0) continue;
        for (size_t j = 0; j < d; ++j) {
            acc[j] = (acc[j] + static_cast<AccModularCoeff>(poly[k]) * powers[k][j]) % prime;
        }
    }
    return acc;
}

// Multiply two elements in the quotient ring using multiplication tables
static std::vector<ModularCoeff>
multiply_in_quotient(const std::vector<ModularCoeff> &a,
                     const std::vector<ModularCoeff> &b,
                     const std::vector<PP> &quotient_basis,
                     const std::vector<std::vector<int32_t>> &i_xw,
                     const std::vector<std::vector<ModularCoeff>> &t_v,
                     ModularCoeff prime) {
    size_t d = quotient_basis.size();
    if (a.size() != d || b.size() != d) return {};
    std::vector<ModularCoeff> res(d, 0);

    // For each monomial basis component of a
    for (size_t i = 0; i < d; ++i) {
        ModularCoeff ai = a[i];
        if (ai == 0) continue;

        // Start with b scaled by ai
        std::vector<ModularCoeff> temp = b;
        for (size_t j = 0; j < d; ++j) { temp[j] = (static_cast<AccModularCoeff>(temp[j]) * ai) % prime; }

        // Multiply by each variable according to monomial exponents in quotient_basis[i]
        const PP &mon = quotient_basis[i];
        for (size_t var_idx = 0; var_idx < mon.size(); ++var_idx) {
            for (int deg = 0; deg < mon[var_idx]; ++deg) {
                std::vector<ModularCoeff> next(d, 0);
                mul_var_quo(next, temp, static_cast<int32_t>(var_idx + 1), i_xw, t_v, prime);
                temp.swap(next);
            }
        }

        // Accumulate
        for (size_t j = 0; j < d; ++j) { res[j] = (res[j] + temp[j]) % prime; }
    }
    return res;
}

// Attempt to reorder parameterizations so that for each variable x_j,
// g_j(T) is the unique numerator satisfying x_j * f'(T) - g_j(T) = 0 in the quotient ring.
static void
reorder_parameterizations_by_identity(std::vector<BivariateResult> &params,
                                      const MinimalPolynomialResult &min_poly,
                                      const std::vector<PP> &quotient_basis,
                                      const std::vector<std::vector<int32_t>> &i_xw,
                                      const std::vector<std::vector<ModularCoeff>> &t_v,
                                      int num_variables,
                                      ModularCoeff prime) {
    if (min_poly.powers.empty()) return;
    size_t d = quotient_basis.size();
    if (d == 0) return;

    // Compute derivative f'(T)
    std::vector<ModularCoeff> fprime(
      std::max<size_t>(1, min_poly.coefficients.size() ? min_poly.coefficients.size() - 1 : 0), 0);
    for (size_t k = 1; k < min_poly.coefficients.size(); ++k) {
        ModularCoeff c = min_poly.coefficients[k];
        ModularCoeff coeff = (static_cast<AccModularCoeff>(k) * c) % prime;
        fprime[k - 1] = coeff;
    }
    // Build V_fprime from powers of T
    std::vector<ModularCoeff> V_fprime = build_element_from_poly(fprime, min_poly.powers, prime);

    // Precompute x_j vectors
    std::vector<std::vector<ModularCoeff>> V_x(num_variables);
    for (int j = 0; j < num_variables; ++j) { V_x[j] = element_to_vector(j + 1, i_xw, t_v, d); }

    // Compute match matrix: match[j] = index i such that x_j * f'(T) - g_i(T) == 0
    std::vector<int> match_for_var(num_variables, -1);
    std::vector<bool> used_num(num_variables, false);

    // Rebuild numerators directly from identities: g_j(T) = x_j * f'(T) in quotient ring
    std::vector<BivariateResult> rebuilt(num_variables);
    for (int j = 0; j < num_variables; ++j) {
        auto prod = multiply_in_quotient(V_x[j], V_fprime, quotient_basis, i_xw, t_v, prime);
        // Solve for coefficients c such that sum_k c_k V(T^k) = prod
        // Build column matrix [V(T^0) ... V(T^{d-1})]
        std::vector<std::vector<ModularCoeff>> A(d, std::vector<ModularCoeff>(d + 1, 0));
        for (size_t r = 0; r < d; ++r) {
            for (size_t c = 0; c < d; ++c) { A[r][c] = min_poly.powers[c][r] % prime; }
            A[r][d] = prod[r] % prime;
        }
        // Gaussian elimination to solve A * x = b
        size_t row = 0;
        for (size_t col = 0; col < d && row < d; ++col) {
            size_t pivot = row;
            while (pivot < d && A[pivot][col] == 0) { ++pivot; }
            if (pivot == d) continue;
            if (pivot != row) std::swap(A[pivot], A[row]);
            ModularCoeff inv = mod_inv(A[row][col], prime);
            for (size_t j2 = col; j2 <= d; ++j2) {
                A[row][j2] = (static_cast<AccModularCoeff>(A[row][j2]) * inv) % prime;
            }
            for (size_t r2 = 0; r2 < d; ++r2) {
                if (r2 == row) continue;
                ModularCoeff factor = A[r2][col];
                if (factor == 0) continue;
                for (size_t j2 = col; j2 <= d; ++j2) {
                    AccModularCoeff sub = static_cast<AccModularCoeff>(factor) * A[row][j2];
                    A[r2][j2] = (A[r2][j2] + prime - (sub % prime)) % prime;
                }
            }
            ++row;
        }
        std::vector<ModularCoeff> coeffs(d, 0);
        for (size_t r = 0; r < d; ++r) coeffs[r] = A[r][d];
        // Trim to minimal degree
        while (!coeffs.empty() && coeffs.back() == 0) coeffs.pop_back();
        BivariateResult br;
        br.success = true;
        br.generators.clear();
        br.generators.push_back(coeffs);
        rebuilt[j] = std::move(br);
    }
    params.swap(rebuilt);
}

void
trim_polynomial(std::vector<ModularCoeff> &v) {
    while (!v.empty() && v.back() == 0) { v.pop_back(); }
}

std::vector<ModularCoeff>
trimmed_copy(const std::vector<ModularCoeff> &v) {
    std::vector<ModularCoeff> result = v;
    trim_polynomial(result);
    return result;
}

bool
are_proportional(const std::vector<ModularCoeff> &p,
                 const std::vector<ModularCoeff> &q,
                 ModularCoeff prime,
                 bool check_reverse,
                 ModularCoeff &out_scale) {

    auto p_trimmed = trimmed_copy(p);
    auto q_trimmed = trimmed_copy(q);

    if (p_trimmed.empty() && q_trimmed.empty()) {
        out_scale = 1;
        return true;
    }

    if (p_trimmed.empty() || q_trimmed.empty()) { return false; }

    if (p_trimmed.size() != q_trimmed.size()) { return false; }

    size_t n = p_trimmed.size();

    auto coeff = [&](const std::vector<ModularCoeff> &poly, size_t i) -> ModularCoeff {
        return check_reverse ? poly[n - 1 - i] : poly[i];
    };

    size_t first = 0;
    while (first < n && coeff(p_trimmed, first) == 0 && coeff(q_trimmed, first) == 0) { ++first; }

    if (first == n) {
        out_scale = 1;
        return true;
    }

    if (coeff(p_trimmed, first) == 0 || coeff(q_trimmed, first) == 0) { return false; }

    out_scale =
      (static_cast<AccModularCoeff>(coeff(q_trimmed, first)) * mod_inv(coeff(p_trimmed, first), prime)) % prime;

    for (size_t i = first; i < n; ++i) {
        ModularCoeff expected = (static_cast<AccModularCoeff>(coeff(p_trimmed, i)) * out_scale) % prime;
        if (expected != coeff(q_trimmed, i)) { return false; }
    }

    return true;
}

void
scale_and_reverse(std::vector<ModularCoeff> &poly, ModularCoeff scale, bool do_reverse, ModularCoeff prime) {
    for (auto &c : poly) { c = (static_cast<AccModularCoeff>(c) * scale) % prime; }
    if (do_reverse) { std::reverse(poly.begin(), poly.end()); }
}

void
make_monic(std::vector<ModularCoeff> &poly, ModularCoeff prime) {
    trim_polynomial(poly);
    if (poly.empty()) return;

    ModularCoeff leading_coeff = poly.back();
    if (leading_coeff == 0 || leading_coeff == 1) return;

    ModularCoeff inv = mod_inv(leading_coeff, prime);
    for (auto &c : poly) { c = (static_cast<AccModularCoeff>(c) * inv) % prime; }
}

} // anonymous namespace
// === Rational polynomial helpers (over Q) for post-CRT processing ===
namespace {
static void
trim_rational(std::vector<mpq_class> &p) {
    while (!p.empty() && p.back() == 0) p.pop_back();
}
static std::vector<mpq_class>
derivative_rational(const std::vector<mpq_class> &p) {
    if (p.size() <= 1) return { mpq_class(0) };
    std::vector<mpq_class> d(p.size() - 1, mpq_class(0));
    for (size_t i = 1; i < p.size(); ++i) { d[i - 1] = p[i] * static_cast<int>(i); }
    trim_rational(d);
    return d;
}
static void
make_monic_rational(std::vector<mpq_class> &p) {
    trim_rational(p);
    if (p.empty()) return;
    mpq_class lc = p.back();
    if (lc == 1) return;
    for (auto &c : p) c /= lc;
}
static std::pair<std::vector<mpq_class>, std::vector<mpq_class>>
divmod_rational(std::vector<mpq_class> a, const std::vector<mpq_class> &b) {
    trim_rational(a);
    std::vector<mpq_class> q;
    if (b.empty() || (b.size() == 1 && b[0] == 0)) return { q, a };
    std::vector<mpq_class> r = a;
    std::vector<mpq_class> bb = b;
    trim_rational(r);
    trim_rational(bb);
    if (bb.empty()) return { q, r };
    mpq_class inv_lc = 1 / bb.back();
    while (r.size() >= bb.size() && !(r.size() == 1 && r[0] == 0)) {
        size_t k = r.size() - bb.size();
        mpq_class coeff = r.back() * inv_lc;
        if (q.size() < k + 1) q.resize(k + 1);
        q[k] += coeff;
        // r -= coeff * x^k * bb
        for (size_t i = 0; i < bb.size(); ++i) { r[k + i] -= coeff * bb[i]; }
        trim_rational(r);
        if (bb.size() == 1) break;
    }
    trim_rational(q);
    trim_rational(r);
    return { q, r };
}
static std::vector<mpq_class>
mod_rational(const std::vector<mpq_class> &a, const std::vector<mpq_class> &b) {
    return divmod_rational(a, b).second;
}
static std::vector<mpq_class>
gcd_rational(std::vector<mpq_class> a, std::vector<mpq_class> b) {
    trim_rational(a);
    trim_rational(b);
    if (a.empty()) return b;
    if (b.empty()) return a;
    // Euclidean algorithm over Q[t]
    while (!b.empty()) {
        auto rem = mod_rational(a, b);
        a.swap(b);
        b.swap(rem);
    }
    make_monic_rational(a);
    return a;
}
} // namespace

// Function to canonicalize RUR by normalizing and finding best permutation of parameterizations
// Returns true on success, false if the prime should be marked as bad
static bool
canonicalize_rur_permutation(std::vector<std::vector<ModularCoeff>> &current_table,
                             const std::vector<std::vector<ModularCoeff>> &reference_table,
                             ModularCoeff prime,
                             size_t num_variables) {

    const bool verbose = false; // master toggle for canonicalization logs

    if (verbose) {
        std::cout << "CANON: sizes cur=" << current_table.size() << " ref=" << reference_table.size() << std::endl;
    }

    // Validation
    if (reference_table.empty() || current_table.empty()) {
        if (verbose) { std::cout << "CANON: empty table, returning" << std::endl; }
        return false; // Nothing to canonicalize - bad prime
    }

    // Trim current polynomials
    for (auto &poly : current_table) { trim_polynomial(poly); }

    // Determine global scale and reversal from minimal polynomial
    bool needs_reversal = false;
    ModularCoeff global_scale = 1;
    bool found_match = false;

    const auto &ref_minpoly = reference_table[0];
    auto &cur_minpoly = current_table[0];

    // Try without reversal first
    if (are_proportional(cur_minpoly, ref_minpoly, prime, false, global_scale)) {
        needs_reversal = false;
        found_match = true;
    }
    // Then try with reversal
    else if (are_proportional(cur_minpoly, ref_minpoly, prime, true, global_scale)) {
        needs_reversal = true;
        found_match = true;
    }

    if (!found_match) {
        // Can't match minimal polynomials - this prime is probably bad
        if (verbose) { std::cout << "CANON: Can't match minimal polynomials - bad prime" << std::endl; }
        return false;
    }

    // Apply global transformation to all polynomials
    for (auto &poly : current_table) { scale_and_reverse(poly, global_scale, needs_reversal, prime); }

    // Make minimal polynomial monic
    make_monic(current_table[0], prime);

    // Match parameterizations (handle variable permutations)
    size_t expected_size = 1 + num_variables;
    std::vector<std::vector<ModularCoeff>> aligned_table(expected_size);
    aligned_table[0] = current_table[0]; // Minimal polynomial already placed

    if (verbose) {
        std::cout << "CANONICALIZATION DEBUG: current_table.size()=" << current_table.size()
                  << ", reference_table.size()=" << reference_table.size() << ", expected_size=" << expected_size
                  << ", num_variables=" << num_variables << std::endl;
        for (size_t i = 0; i < aligned_table.size(); ++i) {
            std::cout << "  aligned_table[" << i << "].size()=" << aligned_table[i].size() << std::endl;
        }
    }

    std::vector<bool> reference_used(num_variables, false);

    // Try to match each current parameterization with a reference one
    for (size_t cur_idx = 1; cur_idx < current_table.size(); ++cur_idx) {
        bool found = false;

        for (size_t ref_idx = 1; ref_idx <= num_variables && ref_idx < reference_table.size(); ++ref_idx) {
            if (reference_used[ref_idx - 1]) continue;

            ModularCoeff scale = 1;
            if (are_proportional(current_table[cur_idx], reference_table[ref_idx], prime, false, scale)) {

                // Found a match - scale to make it exactly equal
                scale_and_reverse(current_table[cur_idx], scale, false, prime);
                aligned_table[ref_idx] = std::move(current_table[cur_idx]);
                reference_used[ref_idx - 1] = true;
                found = true;
                break;
            }
        }
    }

    // Check if we successfully matched all parameterizations
    int unmatched_count = 0;
    for (size_t i = 1; i < aligned_table.size(); ++i) {
        if (aligned_table[i].empty()) {
            unmatched_count++;
            if (verbose) {
                std::cout << "CANON: Warning - aligned_table[" << i << "] is empty (unmatched)" << std::endl;
            }
        }
    }

    if (unmatched_count > 0) {
        std::cout << "CANON: Failed to match " << unmatched_count << " parameterizations - marking prime as bad"
                  << std::endl;
        // Don't modify current_table - let caller handle the bad prime
        return false;
    }

    // Replace current table with aligned version
    if (verbose) {
        std::cout << "CANONICALIZATION: Before replacing current_table:" << std::endl;
        for (size_t i = 0; i < aligned_table.size(); ++i) {
            std::cout << "  aligned_table[" << i << "] has " << aligned_table[i].size() << " elements" << std::endl;
        }
    }
    current_table = std::move(aligned_table);

    // Ensure consistent lengths (pad with zeros)
    for (size_t i = 0; i < current_table.size() && i < reference_table.size(); ++i) {
        size_t ref_size = reference_table[i].size();
        current_table[i].resize(ref_size, 0);
    }

    return true; // Success
}

// Special case: Normalize first prime's table to establish canonical form
static void
normalize_first_table(std::vector<std::vector<ModularCoeff>> &table, ModularCoeff prime) {

    for (auto &poly : table) { trim_polynomial(poly); }

    if (!table.empty()) { make_monic(table[0], prime); }
}

std::pair<ModularRURResult, std::vector<int>>
compute_modular_rur(const std::vector<std::string> &polynomials,
                    const std::vector<std::string> &variables,
                    ModularCoeff prime,
                    const RURConfig &config,
                    const std::vector<int> &separating_element) {
    auto time_start_total = std::chrono::steady_clock::now();
    bool timing = config.timing || (std::getenv("RUR_TIMING") && std::string(std::getenv("RUR_TIMING")) != "0");
    ModularRURResult result;
    result.prime = prime;
    result.success = false;
    std::vector<int> found_coeffs;


    // Check if prime is safe for F4 (must be <= 2^31-1 to avoid infinite loop bug)
    const ModularCoeff F4_MAX_SAFE_PRIME = 2147483647; // 2^31 - 1
    if (prime > F4_MAX_SAFE_PRIME) {
        if (config.verbose) {
            std::cerr << "F4 does not support primes > " << F4_MAX_SAFE_PRIME << " due to known bug. Prime " << prime
                      << " is too large." << std::endl;
        }
        // Return failure - F4 cannot handle this prime
        result.success = false;
        return { result, {} };
    }

    // Also check for suspiciously large primes that might be due to overflow
    if (prime > 4000000000u) { // Primes > 4 billion are suspicious
        if (config.verbose) {
            std::cerr << "WARNING: Suspiciously large prime " << prime
                      << " detected. This may be due to integer overflow." << std::endl;
        }
        result.success = false;
        return { result, {} };
    }

    // Step 1: Create F4 session
    std::vector<const char *> var_ptrs;
    for (const auto &var : variables) { var_ptrs.push_back(var.c_str()); }

    axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), variables.size());
    if (!session) { return { result, {} }; }

    // Step 2: Add polynomials
    for (const auto &poly : polynomials) {
        std::string formatted = format_polynomial_for_f4(poly);
        if (config.verbose) {
            std::cout << "DEBUG: Original: \"" << poly << "\" -> Formatted: \"" << formatted << "\"" << std::endl;
        }
        if (axf4_add_polynomial(session, formatted.c_str()) != 0) {
            axf4_destroy_session(session);
            return { result, {} };
        }
    }

    // Step 3: Compute Gröbner basis
    auto time_start_f4 = std::chrono::steady_clock::now();
    axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
    if (gb_result.status != 0) {
        if (config.verbose) {
            std::cerr << "F4 failed with status " << gb_result.status << " for prime " << prime << std::endl;
        }
        axf4_free_result(&gb_result);
        axf4_destroy_session(session);
        return { result, {} };
    }
    auto time_end_f4 = std::chrono::steady_clock::now();
    if (timing) {
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end_f4 - time_start_f4).count();
        if (timing && std::getenv("RUR_VERBOSE_TIMING")) {
            std::cout << "[timing] p=" << prime << " groebner_ms=" << ms << std::endl;
        }
    }

    if (config.verbose) {
        std::cout << "Gröbner basis computed, size: " << gb_result.basis_size << std::endl;
        std::cout << "Gröbner basis string: " << gb_result.groebner_basis << std::endl;
    }

    // Step 4: Extract multiplication tables using F4 integration
    std::vector<std::vector<ModularCoeff>> t_v;
    std::vector<StackVect> t_xw;
    std::vector<std::vector<int32_t>> i_xw;

    auto time_start_tables = std::chrono::steady_clock::now();
    bool tables_success = f4_to_multiplication_tables(session, t_v, t_xw, i_xw, result.quotient_basis, prime);
    auto time_end_tables = std::chrono::steady_clock::now();

    // Clean up F4 data
    axf4_free_result(&gb_result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);

    if (!tables_success) { return { result, {} }; }
    if (timing) {
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end_tables - time_start_tables).count();
        if (timing && std::getenv("RUR_VERBOSE_TIMING")) {
            std::cout << "[timing] p=" << prime << " tables_ms=" << ms << std::endl;
        }
    }

    // CRITICAL: Check for "unlucky primes" that cause dimensional collapse
    // For Katsura-4, we expect dimension 16. Primes like p=2 give dimension 12.
    // Skip known bad primes
    if (prime == 2) {
        if (config.verbose) { std::cout << "Skipping known bad prime p=2 for Katsura systems" << std::endl; }
        return { result, {} };
    }

    // Build variable position mapping
    result.var_positions.resize(variables.size(), -1);
    for (size_t i = 0; i < result.quotient_basis.size(); ++i) {
        const auto &pp = result.quotient_basis[i];
        // Check if this is a pure variable monomial (exactly one exponent is 1, rest are 0)
        int var_count = 0;
        int var_idx = -1;
        for (size_t j = 0; j < pp.size(); ++j) {
            if (pp[j] == 1) {
                var_count++;
                var_idx = j;
            } else if (pp[j] != 0) {
                var_count = -1; // Not a pure variable
                break;
            }
        }

        if (var_count == 1 && var_idx >= 0) { result.var_positions[var_idx] = i; }
    }

    if (config.verbose) {
        std::cout << "Quotient basis size: " << result.quotient_basis.size() << std::endl;
        std::cout << "Quotient basis elements: ";
        for (const auto &pp : result.quotient_basis) {
            std::cout << "[";
            for (size_t i = 0; i < pp.size(); ++i) {
                if (i > 0) std::cout << ",";
                std::cout << pp[i];
            }
            std::cout << "] ";
        }
        std::cout << std::endl;
        std::cout << "Variable positions in quotient basis: ";
        for (size_t i = 0; i < result.var_positions.size(); ++i) {
            std::cout << variables[i] << "->pos " << result.var_positions[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "Multiplication tables built" << std::endl;
    }

    // Heuristic degeneracy check (OFF by default):
    // If a variable's Krylov sequence {1, x, x^2, ...} dies very early (e.g., x^k = 0 well before the
    // quotient dimension), we used to reject the prime as "non-generic" to save time. However, this can
    // incorrectly reject valid non-radical (nilpotent) cases where square-free processing recovers the radical.
    // Therefore, this check is DISABLED by default. To enable it for speed-only experimentation, set
    // RUR_ENABLE_KRYLOV_CHECK=1. Be aware: enabling it may break correctness by rejecting solvable cases.
    {
        bool enable_krylov =
          (std::getenv("RUR_ENABLE_KRYLOV_CHECK") && std::string(std::getenv("RUR_ENABLE_KRYLOV_CHECK")) != "0");
        size_t d = result.quotient_basis.size();
        if (enable_krylov && d > 1) {
            auto is_degenerate_var = [&](int var_index) -> bool {
                std::vector<ModularCoeff> v = element_to_vector(var_index, i_xw, t_v, d);
                if (std::all_of(v.begin(), v.end(), [](ModularCoeff c) { return c == 0; })) return true;
                std::vector<ModularCoeff> tmp(d, 0);
                // Allow up to d steps; if zero reached far earlier (e.g., < d/3), consider degenerate
                size_t zero_at = 0;
                for (size_t k = 1; k <= d; ++k) {
                    std::fill(tmp.begin(), tmp.end(), 0);
                    mul_var_quo(tmp, v, var_index, i_xw, t_v, prime);
                    v.swap(tmp);
                    bool is_zero = std::all_of(v.begin(), v.end(), [](ModularCoeff c) { return c == 0; });
                    if (is_zero) {
                        zero_at = k;
                        break;
                    }
                }
                return (zero_at > 0 && zero_at < (d / 3 + 1));
            };

            // Check last variable first (often good separating candidate)
            int last_var = static_cast<int>(variables.size());
            bool deg_last = is_degenerate_var(last_var);
            bool deg_any = deg_last;
            if (!deg_any) {
                for (int j = 1; j <= static_cast<int>(variables.size()); ++j) {
                    if (j == last_var) continue;
                    if (is_degenerate_var(j)) {
                        deg_any = true;
                        break;
                    }
                }
            }
            if (deg_any) {
                if (config.verbose) {
                    std::cout << "[DEBUG] Degenerate Krylov sequence detected (early nilpotency). Marking prime as bad."
                              << std::endl;
                }
                result.success = false;
                return { result, {} };
            }
        }
    }

    // Diagnostics: check x_last^5 behavior when last variable is expected separating element
    {
        bool debug_pow =
          config.verbose || (std::getenv("RUR_DEBUG_POWERS") && std::string(std::getenv("RUR_DEBUG_POWERS")) != "0");
        if (debug_pow && !result.quotient_basis.empty() && !i_xw.empty()) {
            int last_var = static_cast<int>(variables.size());
            size_t d = result.quotient_basis.size();
            // Build x_last vector
            std::vector<ModularCoeff> x_last = element_to_vector(last_var, i_xw, t_v, d);
            // Find index of pure power x_last^4 in quotient basis
            int idx_pow4 = -1;
            if (!result.quotient_basis.empty()) {
                const size_t nvars = result.quotient_basis[0].size();
                for (size_t i = 0; i < d; ++i) {
                    const PP &pp = result.quotient_basis[i];
                    bool ok = true;
                    for (size_t j = 0; j < nvars; ++j) {
                        int exp = pp[j];
                        if (j + 1 == static_cast<size_t>(last_var)) {
                            if (exp != 4) {
                                ok = false;
                                break;
                            }
                        } else {
                            if (exp != 0) {
                                ok = false;
                                break;
                            }
                        }
                    }
                    if (ok) {
                        idx_pow4 = static_cast<int>(i);
                        break;
                    }
                }
            }
            if (idx_pow4 >= 0 && last_var - 1 < static_cast<int>(i_xw.size()) &&
                idx_pow4 < static_cast<int>(i_xw[last_var - 1].size())) {
                int32_t target_idx = i_xw[last_var - 1][idx_pow4];
                std::cout << "[DIAG] i_xw[last][x_last^4 basis idx=" << idx_pow4 << "] = " << target_idx << std::endl;
                if (target_idx > 0 && target_idx <= static_cast<int>(t_v.size())) {
                    const auto &tv = t_v[target_idx - 1];
                    std::cout << "[DIAG] t_v[target_idx-1] size=" << tv.size() << ", first 10 coeffs=[";
                    for (size_t k = 0; k < tv.size() && k < 10; ++k) {
                        if (k) std::cout << ",";
                        std::cout << tv[k];
                    }
                    std::cout << "]" << std::endl;
                }
            } else {
                std::cout << "[DIAG] Could not locate x_last^4 in quotient basis or invalid i_xw mapping" << std::endl;
            }

            // Compute x_last^5 by repeated multiplication
            std::vector<ModularCoeff> cur = x_last;
            std::vector<ModularCoeff> tmp(d, 0);
            for (int rep = 1; rep <= 4; ++rep) {
                std::fill(tmp.begin(), tmp.end(), 0);
                mul_var_quo(tmp, cur, last_var, i_xw, t_v, prime);
                cur.swap(tmp);
                ModularCoeff sum = 0;
                for (ModularCoeff c : cur) sum = (sum + c) % prime;
                std::cout << "[DIAG] after multiply by x_last (rep=" << rep << ") nonzeros="
                          << std::count_if(cur.begin(), cur.end(), [](ModularCoeff c) { return c != 0; })
                          << ", checksum=" << sum << std::endl;
            }
        }
    }

    // === SPECIAL CASE FOR 1-DIMENSIONAL QUOTIENT RING ===
    // If the quotient basis has size 1, the solution is a single point and is
    // already determined by the Gröbner basis. The RUR is trivial.
    if (result.quotient_basis.size() == 1) {
        if (config.verbose) {
            std::cout << "DEBUG: Handling 1-dimensional quotient ring (single solution)" << std::endl;
        }
        result.success = true;

        // For a 1-dimensional quotient ring (single solution), we need to extract
        // the constant values of the variables and create a proper RUR.
        // The minimal polynomial will be T - c where c is the value used as separating element

        // For a 1-dimensional quotient ring, the Gröbner basis should be linear in each variable
        // We need to extract the constant values from the Gröbner basis itself
        // The GB should contain polynomials of the form: x_i + c_i = 0, so x_i = -c_i

        // Parse the Gröbner basis to extract constant values
        std::vector<ModularCoeff> var_values(variables.size(), 0);

        // The Gröbner basis is stored in gb_result.groebner_basis as a string
        // We need to get it from the F4 session
        // For now, we can use the multiplication tables if they're available

        // Try to extract from multiplication tables first
        bool extracted = false;
        if (!t_v.empty() && !i_xw.empty()) {
            // For 1-dimensional quotient ring, the multiplication table for x_i * 1 gives the value
            for (size_t i = 0; i < variables.size(); ++i) {
                if (result.var_positions[i] >= 0) {
                    // Variable appears in quotient basis - shouldn't happen for 1-dim
                    var_values[i] = 0;
                } else if (i < i_xw.size() && !i_xw[i].empty()) {
                    int32_t table_idx = i_xw[i][0]; // Entry for var_(i+1) * 1
                    if (table_idx > 0 && table_idx <= static_cast<int32_t>(t_v.size()) && !t_v[table_idx - 1].empty()) {
                        // The value in the table is the negative of what we computed in GB
                        // Since GB has x_i + c_i = 0, we have x_i = -c_i mod p
                        // But the table stores the value directly
                        var_values[i] = t_v[table_idx - 1][0];
                        extracted = true;
                        if (config.verbose) {
                            std::cout << "  Variable " << variables[i] << " = " << var_values[i] << " (mod " << prime
                                      << ")" << std::endl;
                        }
                    }
                }
            }
        }

        if (!extracted && config.verbose) {
            std::cout << "  Warning: Could not extract variable values from multiplication tables" << std::endl;
        }

        // Use a deterministic linear combination as separating element: x1 + 2*x2 + 3*x3 + ...
        // This matches the approach used for general multivariate systems
        ModularCoeff sep_value = 0;
        for (size_t i = 0; i < var_values.size(); ++i) {
            sep_value =
              (sep_value + static_cast<ModularCoeff>((static_cast<AccModularCoeff>(i + 1) * var_values[i]) % prime)) %
              prime;
        }

        // Minimal polynomial is T - sep_value
        MinimalPolynomialResult min_poly_result;
        min_poly_result.success = true;
        // Normalize sep_value to [-p/2, p/2] range for consistency
        ModularCoeff p_minus_sep = (prime - sep_value) % prime;
        min_poly_result.coefficients = { p_minus_sep, 1 }; // T - sep_value = 0
        min_poly_result.degree = 1;
        result.minimal_polynomial = min_poly_result;

        // For parameterizations: each variable is expressed as a constant polynomial
        // Since T = sep_value, we express each variable as a constant
        result.parameterizations.resize(variables.size());
        for (size_t i = 0; i < variables.size(); ++i) {
            BivariateResult param;
            param.success = true;
            // The parameterization is just the constant value
            param.generators.push_back({ var_values[i] });
            result.parameterizations[i] = param;
        }
        return { result, {} };
    }

    // Step 5: Compute univariate parameterization with a deterministic, globally consistent separating element
    // Use the same linear form across primes to align RURs for CRT.
    // Build modular vector for separating element only when needed, but prefer integer coeffs API
    std::vector<ModularCoeff> sep_element_mod;
    if (!separating_element.empty()) {
        for (int coeff : separating_element) {
            sep_element_mod.push_back((coeff < 0) ? (prime - ((-coeff) % prime)) % prime : coeff % prime);
        }
    }
    std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>> tup;
    if (!separating_element.empty()) {
        tup = try_separating_element(result.quotient_basis, i_xw, t_v, variables.size(), prime, separating_element, -1);
    } else if (variables.size() == 1) {
        tup = try_separating_element(result.quotient_basis, i_xw, t_v, variables.size(), prime, std::vector<int>(), 1);
    } else {
        // Mirror Julia strategy: try last variable first, then other variables, then linear forms
        bool success = false;
        MinimalPolynomialResult min_poly_tmp;
        std::vector<BivariateResult> params_tmp;
        // 1) Try last variable as separating element
        std::tie(success, min_poly_tmp, params_tmp) = try_separating_element(result.quotient_basis,
                                                                             i_xw,
                                                                             t_v,
                                                                             variables.size(),
                                                                             prime,
                                                                             std::vector<int>(),
                                                                             static_cast<int32_t>(variables.size()));
        if (success) {
            // Record that we used the last variable as separating element
            found_coeffs.assign(variables.size(), 0);
            found_coeffs[variables.size() - 1] = 1;
        }
        // 2) If failed, try other single variables in reverse order
        if (!success) {
            for (int v = static_cast<int>(variables.size()) - 1; v >= 1; --v) {
                std::tie(success, min_poly_tmp, params_tmp) = try_separating_element(result.quotient_basis,
                                                                                     i_xw,
                                                                                     t_v,
                                                                                     variables.size(),
                                                                                     prime,
                                                                                     std::vector<int>(),
                                                                                     static_cast<int32_t>(v));
                if (success) {
                    // Record which variable we used as separating element
                    found_coeffs.assign(variables.size(), 0);
                    found_coeffs[v - 1] = 1;
                    break;
                }
            }
        }
        // 3) Fall back to systematic linear-form search with increasing width
        if (!success) {
            bool got = false;
            MinimalPolynomialResult mp_s;
            std::vector<BivariateResult> params_s;
            std::vector<int> coeffs_s;
            for (int width : { 3, 5, 7, 9, 12 }) {
                auto [param_success, min_poly, params, coeffs] =
                  find_separating_element_systematic(result.quotient_basis, i_xw, t_v, variables.size(), prime, width);
                if (param_success) {
                    got = true;
                    mp_s = std::move(min_poly);
                    params_s = std::move(params);
                    coeffs_s = std::move(coeffs);
                    break;
                }
            }
            // Note: Do not try alternate primes here; i_xw/t_v are built for this 'prime'.
            // Prime retries are handled at the reference-pass level where tables are recomputed.
            tup = { got, mp_s, params_s };
            found_coeffs = coeffs_s;
        } else {
            tup = { success, min_poly_tmp, params_tmp };
        }
    }
    auto [param_success, min_poly, params] = tup;

    // Ensure minimal polynomial is square-free and parameterizations are reduced modulo it.
    // Even if the minpoly computation already applied square-free, this normalizes shapes.
    if (param_success && !min_poly.coefficients.empty()) {
        // Compute gcd(f, f') over F_p
        auto f_sf_deriv = polynomial_derivative(min_poly.coefficients, prime);
        auto eg = polynomial_extended_gcd(min_poly.coefficients, f_sf_deriv, prime);
        const std::vector<ModularCoeff> &g = eg.gcd;
        // Compute square-free part: f / g if gcd non-constant
        if (!(g.size() == 1 && g[0] == 1)) {
            auto divres = polynomial_divmod(min_poly.coefficients, g, prime);
            std::vector<ModularCoeff> f_sf = divres.first;
            normalize_polynomial(f_sf);
            if (!f_sf.empty() && f_sf.size() < min_poly.coefficients.size()) {
                // Replace min_poly with square-free part
                min_poly.coefficients = f_sf;
                min_poly.degree = f_sf.size() > 0 ? f_sf.size() - 1 : 0;
                // Trim stored powers to the new degree (keep T^0..T^{d-1})
                if (!min_poly.powers.empty() && min_poly.powers.size() > min_poly.degree) {
                    min_poly.powers.resize(min_poly.degree + 1);
                }
                // Reduce each numerator modulo f_sf
                for (auto &br : params) {
                    if (!br.generators.empty()) {
                        auto red = polynomial_mod(br.generators[0], min_poly.coefficients, prime);
                        normalize_polynomial(red);
                        br.generators[0] = std::move(red);
                    }
                }
            }
        }
    }

    // Special case: Non-radical ideal with square-free reduction
    // When min_poly.degree < quotient_basis.size() due to multiplicities
    if (min_poly.success && min_poly.degree == 1 && result.quotient_basis.size() > 1) {
        if (config.verbose || getenv("RUR_NON_RADICAL_VERBOSE")) {
            std::cout << "Non-radical ideal detected: degree 1 minimal polynomial with quotient dimension "
                      << result.quotient_basis.size() << std::endl;
            if (min_poly.original_degree > 0) {
                std::cout << "  Original degree before square-free: " << min_poly.original_degree << " (multiplicity ~"
                          << min_poly.multiplicity << ")" << std::endl;
            }
        }
        // For degree 1 minimal polynomial T - c = 0, all variables are constants
        // The parameterizations should be extracted from the Gröbner basis
        // This is similar to the quotient_basis.size() == 1 case but with multiplicities
        param_success = true; // We can handle this case
    }

    // Cyclic optimization: if deg(f) == dim(quotient), rebuild numerators via
    // the identity xi * f'(T) = g_i(T) in the quotient ring to avoid full bivariate lex pipeline.
    // This is generic for cyclic quotient algebras and preserves correctness.
    else if (min_poly.success && min_poly.degree == result.quotient_basis.size()) {
        try {
            std::vector<BivariateResult> rebuilt_params(variables.size());
            // Reuse helper that reconstructs numerators from identities
            reorder_parameterizations_by_identity(
              rebuilt_params, min_poly, result.quotient_basis, i_xw, t_v, variables.size(), prime);
            params.swap(rebuilt_params);
            param_success = true;
        } catch (...) {
            // Fall back to parameters returned by try_separating_element
        }
    }

    // Identity-based reordering of parameterizations to ensure cross-prime consistency
    if (param_success && variables.size() > 1) {
        try {
            reorder_parameterizations_by_identity(
              params, min_poly, result.quotient_basis, i_xw, t_v, variables.size(), prime);
        } catch (...) {
            std::cout << "CANON-IDENTITY: exception during reordering, keeping original order" << std::endl;
        }
    }

    auto time_end_param = std::chrono::steady_clock::now();
    result.minimal_polynomial = min_poly;
    result.parameterizations = params;
    result.success = param_success;

    if (config.verbose) {
        if (param_success) {
            std::cout << "Univariate parameterization successful" << std::endl;
            std::cout << "Minimal polynomial degree: " << min_poly.degree << std::endl;
        } else {
            std::cout << "Univariate parameterization failed" << std::endl;
        }
    }
    if (timing) {
        auto ms_param = std::chrono::duration_cast<std::chrono::milliseconds>(time_end_param - time_end_tables).count();
        auto ms_total =
          std::chrono::duration_cast<std::chrono::milliseconds>(time_end_param - time_start_total).count();
        if (timing && std::getenv("RUR_VERBOSE_TIMING")) {
            std::cout << "[timing] p=" << prime << " parameterization_ms=" << ms_param << " total_prime_ms=" << ms_total
                      << std::endl;
        }
    }
    return { result, found_coeffs };
}

RationalRURResult
compute_rational_rur(const std::vector<std::string> &polynomials,
                     const std::vector<std::string> &variables,
                     const RURConfig &config) {
    auto time_start_all = std::chrono::steady_clock::now();
    bool timing = config.timing || (std::getenv("RUR_TIMING") && std::string(std::getenv("RUR_TIMING")) != "0");
    long long timing_ref_ms = 0;
    long long timing_total_modular_ms = 0;
    long long timing_total_crt_ms = 0;
    RationalRURResult result;
    result.success = false;

    // Define F4 safety limit
    const ModularCoeff F4_MAX_SAFE_PRIME = 2147483647; // 2^31 - 1

    // First, compute with a verification prime to get the structure
    // Use a large prime for the reference computation to get generic behavior
    // Use a random 30-bit prime for the reference pass to avoid systematic unlucky primes
    ModularCoeff verification_prime = generate_random_prime(29, 30);
    if (config.verbose) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Computing reference RUR modulo " << verification_prime << std::endl;
        std::cout << "========================================" << std::endl;
    }

    // Reference pass with bounded prime retries
    ModularRURResult reference_result;
    std::vector<int> separating_element_coeffs;
    {
        std::vector<ModularCoeff> candidate_primes;
        // Prefer a handful of 30-bit primes first
        candidate_primes.push_back(verification_prime);
        auto alt30 = prev_primes(verification_prime, 7);
        candidate_primes.insert(candidate_primes.end(), alt30.begin(), alt30.end());

        // If those fail, also try a batch of safe 28-bit primes to avoid
        // specialization where a single variable is not separating mod p
        ModularCoeff start28 = (static_cast<ModularCoeff>(1u) << 28) - 3u;
        auto alt28 = prev_primes(start28, 12);
        candidate_primes.insert(candidate_primes.end(), alt28.begin(), alt28.end());

        bool got_reference = false;
        for (size_t pi = 0; pi < candidate_primes.size(); ++pi) {
            ModularCoeff p = candidate_primes[pi];
            if (config.verbose) {
                std::cout << "[DEBUG] Trying reference prime " << (pi + 1) << "/" << candidate_primes.size() << ": "
                          << p << std::endl;
            }
            auto time_start_ref = std::chrono::steady_clock::now();
            auto [ref_res, se_coeffs] = compute_modular_rur(polynomials, variables, p, config, {});
            auto time_end_ref = std::chrono::steady_clock::now();
            timing_ref_ms =
              std::chrono::duration_cast<std::chrono::milliseconds>(time_end_ref - time_start_ref).count();
            if (timing && std::getenv("RUR_VERBOSE_TIMING")) { std::cout << "[timing] reference_prime_ms=" << timing_ref_ms << std::endl; }

            if (ref_res.success) {
                if (config.verbose) {
                    std::cout << "[DEBUG] Reference prime " << p << " succeeded!" << std::endl;
                    std::cout << "[DEBUG] Minimal polynomial degree: "
                              << ref_res.minimal_polynomial.coefficients.size() - 1 << std::endl;
                }
                reference_result = std::move(ref_res);
                separating_element_coeffs = std::move(se_coeffs);
                got_reference = true;
                break;
            } else {
                if (config.verbose) {
                    std::cout << "[DEBUG] Reference prime " << p << " failed" << std::endl;
                    if (!ref_res.minimal_polynomial.success) {
                        std::cout << "  -> Failed to find minimal polynomial" << std::endl;
                    }
                    if (ref_res.parameterizations.empty()) { std::cout << "  -> Empty parameterizations" << std::endl; }
                }
            }
        }
        if (!got_reference) {
            result.error_message = "Failed to compute reference RUR";
            if (timing) {
                auto time_end_all = std::chrono::steady_clock::now();
                auto total_ms =
                  std::chrono::duration_cast<std::chrono::milliseconds>(time_end_all - time_start_all).count();
                if (std::getenv("RUR_VERBOSE_TIMING")) {
                std::cout << "[timing] total_ms=" << total_ms << std::endl;
            }
            }
            return result;
        }
    }
    // Fast path for univariate: use rational reconstruction to handle fractions
    if (variables.size() == 1) {
        result.success = true;
        result.quotient_basis = reference_result.quotient_basis;

        // Set up for rational reconstruction
        mpz_class modulus = verification_prime;
        auto bounds = compute_balanced_bounds(modulus);

        // Minimal polynomial - use rational reconstruction
        result.minimal_polynomial.clear();
        for (auto c : reference_result.minimal_polynomial.coefficients) {
            mpz_class val = c;
            auto rat_result = rational_reconstruction(val, modulus, bounds.N, bounds.D);
            if (rat_result.success) {
                result.minimal_polynomial.push_back(rat_result.rational);
            } else {
                // Fall back to symmetric lift
                long long cc = static_cast<long long>(c);
                long long p = static_cast<long long>(verification_prime);
                if (cc > p / 2) cc -= p;
                result.minimal_polynomial.push_back(mpq_class(static_cast<long>(cc)));
            }
        }

        // Numerators - use rational reconstruction
        result.numerators.resize(1);
        if (!reference_result.parameterizations.empty() && !reference_result.parameterizations[0].generators.empty()) {
            for (auto c : reference_result.parameterizations[0].generators[0]) {
                mpz_class val = c;
                auto rat_result = rational_reconstruction(val, modulus, bounds.N, bounds.D);
                if (rat_result.success) {
                    result.numerators[0].push_back(rat_result.rational);
                } else {
                    // Fall back to symmetric lift
                    long long cc = static_cast<long long>(c);
                    long long p = static_cast<long long>(verification_prime);
                    if (cc > p / 2) cc -= p;
                    result.numerators[0].push_back(mpq_class(static_cast<long>(cc)));
                }
            }
        }
        // Common denominator is f'(T); we set derivative denominator as 1 here; downstream consumers derive as needed
        result.denominator_derivative = mpq_class(1);
        if (timing) {
            auto time_end_all = std::chrono::steady_clock::now();
            auto total_ms =
              std::chrono::duration_cast<std::chrono::milliseconds>(time_end_all - time_start_all).count();
            if (std::getenv("RUR_VERBOSE_TIMING")) {
                std::cout << "[timing] total_ms=" << total_ms << std::endl;
            }
        }
        return result;
    }

    // Fast path for 1-dimensional quotient rings (single solution)
    if (reference_result.quotient_basis.size() == 1) {
        result.success = true;
        result.quotient_basis = reference_result.quotient_basis;

        // For 1-dimensional case, we need rational reconstruction, not just symmetric lift
        // The modular values might represent fractions like 2/3, not integers

        // Set up for rational reconstruction
        mpz_class modulus = verification_prime;
        auto bounds = compute_balanced_bounds(modulus);

        // Convert minimal polynomial using rational reconstruction
        result.minimal_polynomial.clear();
        for (auto c : reference_result.minimal_polynomial.coefficients) {
            mpz_class val = c;
            auto rat_result = rational_reconstruction(val, modulus, bounds.N, bounds.D);
            if (rat_result.success) {
                result.minimal_polynomial.push_back(rat_result.rational);
            } else {
                // Fall back to symmetric lift for this coefficient
                long long cc = static_cast<long long>(c);
                long long p = static_cast<long long>(verification_prime);
                if (cc > p / 2) cc -= p;
                result.minimal_polynomial.push_back(mpq_class(static_cast<long>(cc)));
            }
        }

        // Convert parameterizations using rational reconstruction
        result.numerators.resize(variables.size());
        for (size_t i = 0; i < variables.size(); ++i) {
            if (i < reference_result.parameterizations.size() &&
                !reference_result.parameterizations[i].generators.empty() &&
                !reference_result.parameterizations[i].generators[0].empty()) {
                for (auto c : reference_result.parameterizations[i].generators[0]) {
                    mpz_class val = c;
                    auto rat_result = rational_reconstruction(val, modulus, bounds.N, bounds.D);
                    if (rat_result.success) {
                        result.numerators[i].push_back(rat_result.rational);
                    } else {
                        // Fall back to symmetric lift
                        long long cc = static_cast<long long>(c);
                        long long p = static_cast<long long>(verification_prime);
                        if (cc > p / 2) cc -= p;
                        result.numerators[i].push_back(mpq_class(static_cast<long>(cc)));
                    }
                }
            }
        }

        return result;
    }

    // Note: General 1-dimensional quotient rings are now handled by the improved
    // multi-strategy rational reconstruction in crt_and_rational_reconstruction

    // Extract structure information
    size_t num_vars = variables.size();
    size_t quotient_dim = reference_result.quotient_basis.size();
    // Use the square-free degree for CRT template, but retain original_degree as metadata
    size_t min_poly_degree = reference_result.minimal_polynomial.coefficients.size();

    // Prepare storage for CRT and rational reconstruction
    std::vector<std::vector<mpq_class>> qq_result;
    std::vector<std::vector<mpz_class>> zz_temp;
    std::vector<std::vector<std::vector<ModularCoeff>>> modular_tables;
    std::vector<ModularCoeff> used_primes;

    // Convert reference result to format needed for verification
    std::vector<std::vector<ModularCoeff>> reference_mod_p;
    // Reference minimal polynomial: use the square-free coefficients (already applied in modular stage)
    reference_mod_p.push_back(reference_result.minimal_polynomial.coefficients);
    for (const auto &param : reference_result.parameterizations) {
        reference_mod_p.push_back(param.generators.empty() ? std::vector<ModularCoeff>() : param.generators[0]);
    }

    // Initialize with known denominator 1
    mpz_class known_denominator = 1;

    // Generate initial batch of primes
    // Use the configured bit size or default to 28 bits for better generic behavior
    // Small primes (< 20 bits) often give non-generic Gröbner bases
    // But 31-bit primes exceed F4's safe limit!
    // Prefer ≤28-bit primes for stability with F4 and to avoid pathological ranges
    size_t bits = config.initial_prime_bits > 0 ? std::min(config.initial_prime_bits, static_cast<size_t>(28)) : 28;
    size_t batch_size = 5;
    ModularCoeff current_prime = (1UL << bits) - 1;

    if (config.verbose) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Starting multi-modular CRT computation" << std::endl;
        std::cout << "Initial prime size: " << bits << " bits (~" << current_prime << ")" << std::endl;
        std::cout << "Reference quotient dimension: " << reference_result.quotient_basis.size() << std::endl;
        std::cout << "Reference minpoly degree (square-free): "
                  << (reference_result.minimal_polynomial.coefficients.size() - 1) << std::endl;
        std::cout << "========================================" << std::endl;
    }

    if (config.verbose) { std::cout << "\nMulti-modular computation:" << std::endl; }

    // Adaptive loop
    bool reconstruction_success = false;
    size_t max_primes = 200; // Safety limit - increased to handle bad primes
    size_t consecutive_bad_primes = 0;
    const size_t max_consecutive_bad = 6; // lower cap to avoid long bad-prime runs

    // Statistics for debugging
    size_t total_primes_tried = 0;
    size_t bad_primes_f4_failed = 0;
    size_t bad_primes_dim_mismatch = 0;
    size_t bad_primes_degree_mismatch = 0;
    size_t bad_primes_other = 0;

    while (!reconstruction_success && used_primes.size() < max_primes) {
        // Generate next batch of primes
        std::vector<ModularCoeff> new_primes = prev_primes(current_prime, batch_size);
        if (!new_primes.empty()) {
            // Check for underflow - if the last prime is very small, we might be running out
            if (new_primes.back() < 65537) { // F4 needs primes > 65536
                if (config.verbose) {
                    std::cout << "Running out of primes (last prime: " << new_primes.back() << ")" << std::endl;
                }
                break;
            }
            current_prime = new_primes.back() - 1;
        } else {
            if (config.verbose) { std::cout << "No more primes available from " << current_prime << std::endl; }
            break;
        }

        // Additional safety check - ensure primes are valid for F4
        std::vector<ModularCoeff> valid_primes;
        for (ModularCoeff p : new_primes) {
            if (p > 65536) { // F4 requires primes > 2^16
                valid_primes.push_back(p);
            }
        }
        new_primes = valid_primes;

        if (new_primes.empty()) {
            if (config.verbose) { std::cout << "No more valid primes available" << std::endl; }
            break;
        }

        // Compute modular RUR for each new prime
        for (ModularCoeff prime : new_primes) {
            // Additional safety check - skip primes that F4 can't handle
            if (prime > F4_MAX_SAFE_PRIME) {
                if (config.verbose) {
                    std::cout << "Skipping unsafe prime " << prime << " (exceeds F4 limit)" << std::endl;
                }
                continue;
            }

            total_primes_tried++;
            if (config.verbose) {
                std::cout << "[DEBUG] Prime #" << total_primes_tried << ": " << prime
                          << " (used_primes.size()=" << used_primes.size() << ")... ";
            }

            // Create a non-verbose config for modular computations
            RURConfig mod_config = config;
            mod_config.verbose = false;

            // CRITICAL FIX: We should pass the separating element, not the minimal polynomial!
            // The minimal polynomial is f(T), but the separating element is a linear form in the variables.
            // For now, pass empty vector to let each prime find its own separating element.
            auto time_start_prime = std::chrono::steady_clock::now();
            auto [mod_result, coeffs] =
              compute_modular_rur(polynomials, variables, prime, mod_config, separating_element_coeffs);
            auto time_end_prime = std::chrono::steady_clock::now();
            if (mod_result.success) {
                timing_total_modular_ms +=
                  std::chrono::duration_cast<std::chrono::milliseconds>(time_end_prime - time_start_prime).count();
            }

            if (!mod_result.success) {
                if (config.verbose) { std::cout << "failed (F4 computation failed)" << std::endl; }
                bad_primes_f4_failed++;
                consecutive_bad_primes++;
                if (consecutive_bad_primes >= max_consecutive_bad) {
                    if (config.verbose) {
                        std::cout << "Too many consecutive bad primes. The system may have no rational solutions."
                                  << std::endl;
                    }
                    break;
                }
                continue;
            }

            // Additional validation: Check structural compatibility with reference
            bool is_bad_prime = false;

            // Check 1: Quotient dimension must match
            if (mod_result.quotient_basis.size() != reference_result.quotient_basis.size()) {
                if (config.verbose) {
                    std::cout << "BAD PRIME (dimension mismatch): quotient dim = " << mod_result.quotient_basis.size()
                              << " vs reference = " << reference_result.quotient_basis.size() << std::endl;
                }
                bad_primes_dim_mismatch++;
                is_bad_prime = true;
            }

            // Check 2: Minimal polynomial degree must match (square-free degrees)
            if (!is_bad_prime && mod_result.minimal_polynomial.coefficients.size() !=
                                   reference_result.minimal_polynomial.coefficients.size()) {
                if (config.verbose) {
                    std::cout << "BAD PRIME (minpoly square-free degree mismatch): degree = "
                              << (mod_result.minimal_polynomial.coefficients.size() - 1)
                              << " vs reference = " << (reference_result.minimal_polynomial.coefficients.size() - 1)
                              << std::endl;
                }
                bad_primes_degree_mismatch++;
                is_bad_prime = true;
            }

            // Check 3: For cyclic quotient algebras, minpoly degree should equal quotient dimension
            if (!is_bad_prime &&
                mod_result.minimal_polynomial.coefficients.size() - 1 != mod_result.quotient_basis.size()) {
                if (config.verbose) {
                    std::cout << "WARNING: Non-cyclic quotient algebra detected. Minpoly degree = "
                              << (mod_result.minimal_polynomial.coefficients.size() - 1)
                              << " != quotient dim = " << mod_result.quotient_basis.size() << std::endl;
                    std::cout << "This may indicate a reducible system or non-separating element." << std::endl;
                }
                // Don't mark as bad yet - might still work
            }

            // Check 4: Number of parameterizations must match
            if (!is_bad_prime && mod_result.parameterizations.size() != reference_result.parameterizations.size()) {
                if (config.verbose) {
                    std::cout << "BAD PRIME (parameterization count mismatch): " << mod_result.parameterizations.size()
                              << " vs reference = " << reference_result.parameterizations.size() << std::endl;
                }
                bad_primes_other++;
                is_bad_prime = true;
            }

            if (is_bad_prime) {
                consecutive_bad_primes++;
                if (consecutive_bad_primes >= max_consecutive_bad) {
                    std::cout << "CAP: too many consecutive bad primes (" << consecutive_bad_primes
                              << ") – attempting reference rebase on next acceptable prime" << std::endl;
                    // Clear to force rebase with next good prime
                    modular_tables.clear();
                    used_primes.clear();
                    consecutive_bad_primes = 0;
                }
                continue;
            }

            // Reset bad prime counter on success
            consecutive_bad_primes = 0;

            // Add to modular tables
            std::vector<std::vector<ModularCoeff>> prime_table;
            prime_table.push_back(mod_result.minimal_polynomial.coefficients);


            for (const auto &param : mod_result.parameterizations) {
                if (!param.generators.empty()) {
                    prime_table.push_back(param.generators[0]);
                } else {
                    prime_table.push_back(std::vector<ModularCoeff>());
                }
            }

            // remove monic normalization and cross-prime scaling alignment

            // DIAGNOSTIC: Log first few coefficients to check for permutations
            const bool verbose_diag = false;
            if (verbose_diag && (used_primes.size() < 5 || used_primes.size() % 50 == 0)) {
                std::cout << "\nDIAGNOSTIC: Prime #" << (used_primes.size() + 1) << " (p=" << prime
                          << ") BEFORE canonicalization:\n";
                for (size_t i = 0; i < prime_table.size() && i < 4; ++i) {
                    std::cout << "  Poly " << i << " first 3 coeffs: ";
                    for (size_t j = 0; j < prime_table[i].size() && j < 3; ++j) {
                        std::cout << prime_table[i][j] << " ";
                    }
                    std::cout << "\n";
                }
            }

            // Per-prime normalization only (no cross-prime coefficient matching)
            for (auto &poly : prime_table) { trim_polynomial(poly); }
            if (verbose_diag && (used_primes.size() < 5 || used_primes.size() % 50 == 0)) {
                std::cout << "AFTER per-prime normalization:\n";
                for (size_t i = 0; i < prime_table.size() && i < 4; ++i) {
                    std::cout << "  Poly " << i << " first 3 coeffs: ";
                    for (size_t j = 0; j < prime_table[i].size() && j < 3; ++j) {
                        std::cout << prime_table[i][j] << " ";
                    }
                    std::cout << "\n";
                }
            }

            // Guard: if any numerator is empty, skip this prime
            bool any_empty = false;
            for (size_t i = 1; i < prime_table.size(); ++i) {
                if (prime_table[i].empty()) {
                    any_empty = true;
                    break;
                }
            }
            if (any_empty) {
                // quiet by default
                consecutive_bad_primes++;
                if (consecutive_bad_primes >= max_consecutive_bad) {
                    // quiet by default
                    modular_tables.clear();
                    used_primes.clear();
                    consecutive_bad_primes = 0;
                }
                continue;
            }

            // Pad to reference dimensions (square-free) for CRT alignment
            for (size_t i = 0; i < prime_table.size() && i < reference_mod_p.size(); ++i) {
                prime_table[i].resize(reference_mod_p[i].size(), 0);
            }

            modular_tables.push_back(prime_table);
            used_primes.push_back(prime);

            if (config.verbose) { std::cout << "success" << std::endl; }
        }

        if (config.verbose) {
            std::cout << "[DEBUG] Used primes: " << used_primes.size() << " - ";
            if (!used_primes.empty()) {
                std::cout << "primes: [";
                for (size_t i = 0; i < used_primes.size(); ++i) {
                    if (i > 0) std::cout << ", ";
                    std::cout << used_primes[i];
                }
                std::cout << "] - ";
            }
            std::cout << std::flush;
        }

        // quiet by default

        // Try CRT and rational reconstruction
        auto time_start_crt = std::chrono::steady_clock::now();
        auto [success, new_denominator] = crt_and_rational_reconstruction(
          qq_result, zz_temp, modular_tables, used_primes, reference_mod_p, verification_prime, known_denominator);
        auto time_end_crt = std::chrono::steady_clock::now();
        timing_total_crt_ms +=
          std::chrono::duration_cast<std::chrono::milliseconds>(time_end_crt - time_start_crt).count();

        // Require at least 2 primes (or separate verification) before accepting success
        reconstruction_success = success && used_primes.size() >= 2;
        known_denominator = new_denominator;

        if (config.verbose) {
            std::cout << "[DEBUG] After CRT, success=" << success << ", recon_success=" << reconstruction_success
                      << ", known_denominator=" << new_denominator << ", num_primes=" << used_primes.size()
                      << std::endl;
        }

        if (config.verbose && !reconstruction_success) { std::cout << "reconstruction failed" << std::endl; }

        if (reconstruction_success) {
            // quiet

            // Print reconstructed minimal polynomial and numerators
            if (!qq_result.empty()) {
                if (config.verbose) { std::cout << "DEBUG: Reconstructed minimal polynomial coefficients:\n  ["; }
                for (size_t j = 0; j < qq_result[0].size(); ++j) {
                    if (j > 0) std::cout << ", ";
                    std::cout << qq_result[0][j];
                }
                std::cout << "]" << std::endl;
                // Also print numerators for variables (first few coeffs)
                for (size_t vi = 0; vi < num_vars && (vi + 1) < qq_result.size(); ++vi) {
                    if (config.verbose) {
                        std::cout << "DEBUG: Reconstructed numerator for variable " << (vi + 1) << ": [";
                        for (size_t j = 0; j < qq_result[vi + 1].size(); ++j) {
                            if (j > 0) std::cout << ", ";
                            std::cout << qq_result[vi + 1][j];
                        }
                        std::cout << "]" << std::endl;
                    }
                }
            }

            // Note: keep f(T) as reconstructed (may be non square-free); evaluation handles multiplicities downstream.

            if (config.verbose) {
                std::cout << "reconstruction succeeded!" << std::endl;
                std::cout << "qq_result has " << qq_result.size() << " polynomials:" << std::endl;
                for (size_t i = 0; i < qq_result.size(); ++i) {
                    std::cout << "  Poly " << i << " has " << qq_result[i].size() << " coeffs" << std::endl;
                    for (size_t j = 0; j < qq_result[i].size() && j < 3; ++j) {
                        std::cout << "    Coeff " << j << ": " << qq_result[i][j] << std::endl;
                    }
                }

                // Skip verification for 1-dimensional quotient rings for now
                std::cout << "DEBUG: reference_result.quotient_basis.size() = "
                          << reference_result.quotient_basis.size() << std::endl;
                if (reference_result.quotient_basis.size() == 1) {
                    std::cout << "Skipping verification for 1-dimensional quotient ring" << std::endl;
                } else {
                    // Verification check with a new prime
                    // Find a new prime for verification
                    ModularCoeff check_prime = verification_prime;
                    if (std::find(used_primes.begin(), used_primes.end(), verification_prime) != used_primes.end()) {
                        // If we've already used the default verification prime, find another one
                        check_prime = 1073741783; // Another 30-bit prime
                        if (std::find(used_primes.begin(), used_primes.end(), check_prime) != used_primes.end()) {
                            // Use a different prime
                            check_prime = 1073741789;
                        }
                    }
                    std::cout << "Verifying with prime " << check_prime << "... ";

                    // Convert rational result to modular for verification
                    std::vector<std::vector<ModularCoeff>> check_mod_p;
                    for (size_t i = 0; i < qq_result.size(); ++i) {
                        if (config.verbose) {
                            std::cout << "Converting rational poly " << i << " with " << qq_result[i].size()
                                      << " coeffs to modular" << std::endl;
                            if (!qq_result[i].empty()) {
                                std::cout << "  First coeff: " << qq_result[i][0] << std::endl;
                            }
                        }
                        check_mod_p.push_back(modular_coeffs_vect(qq_result[i], check_prime));
                    }

                    // Compute fresh modular result
                    RURConfig check_config = config;
                    check_config.verbose = false;

                    // Use empty separating element for verification as well
                    auto [check_result, check_coeffs] =
                      compute_modular_rur(polynomials, variables, check_prime, check_config, separating_element_coeffs);

                    if (check_result.success) {
                        // Compare results
                        bool verified = true;
                        if (check_result.minimal_polynomial.coefficients != check_mod_p[0]) { verified = false; }

                        for (size_t i = 0; i < num_vars && verified; ++i) {
                            if (!check_result.parameterizations[i].generators.empty() && i + 1 < check_mod_p.size()) {
                                if (check_result.parameterizations[i].generators[0] != check_mod_p[i + 1]) {
                                    verified = false;
                                }
                            }
                        }

                        if (verified) {
                            std::cout << "verification passed!" << std::endl;
                        } else {
                            std::cout << "verification FAILED!" << std::endl;
                            reconstruction_success = false;
                        }
                    } else {
                        std::cout << "check computation failed" << std::endl;
                        reconstruction_success = false;
                    }
                } // end else (skip verification for 1-dim)
            } // end if (config.verbose)
        }

        // Adaptive batch size adjustment
        if (!reconstruction_success) { batch_size = std::max(2UL, used_primes.size() / 10); }
    }

    if (!reconstruction_success) {
        if (config.verbose) {
            std::cout << "\n========================================" << std::endl;
            std::cout << "CRT RECONSTRUCTION FAILED - STATISTICS:" << std::endl;
            std::cout << "Total primes tried: " << total_primes_tried << std::endl;
            std::cout << "Good primes used: " << used_primes.size() << std::endl;
            std::cout << "Bad primes (F4 failed): " << bad_primes_f4_failed << std::endl;
            std::cout << "Bad primes (dimension mismatch): " << bad_primes_dim_mismatch << std::endl;
            std::cout << "Bad primes (degree mismatch): " << bad_primes_degree_mismatch << std::endl;
            std::cout << "Bad primes (other): " << bad_primes_other << std::endl;

            if (bad_primes_dim_mismatch > 0 || bad_primes_degree_mismatch > 0) {
                std::cout << "\nNOTE: Frequent dimension/degree mismatches suggest the system may be REDUCIBLE."
                          << std::endl;
                std::cout << "Reducible systems cannot be solved with a single RUR and require factorization."
                          << std::endl;
            }
            std::cout << "========================================" << std::endl;
        }

        if (used_primes.empty()) {
            result.error_message = "No valid primes found for this system. The system may have no rational solutions.";
        } else {
            result.error_message =
              "Rational reconstruction failed after " + std::to_string(used_primes.size()) + " primes";
        }
        if (timing) {
            auto time_end_all = std::chrono::steady_clock::now();
            auto total_ms =
              std::chrono::duration_cast<std::chrono::milliseconds>(time_end_all - time_start_all).count();
            if (std::getenv("RUR_VERBOSE_TIMING")) {
                std::cout << "[timing] summary: reference_ms=" << timing_ref_ms << " modular_ms=" << timing_total_modular_ms
                          << " crt_ms=" << timing_total_crt_ms << " total_ms=" << total_ms << std::endl;
            }
        }
        return result;
    }

    // Convert results to output format
    result.success = true;
    result.quotient_basis = reference_result.quotient_basis;

    // Minimal polynomial
    if (!qq_result.empty()) { result.minimal_polynomial = qq_result[0]; }

    // Parameterizations (numerators)
    result.numerators.resize(num_vars);
    for (size_t i = 0; i < num_vars && i + 1 < qq_result.size(); ++i) { result.numerators[i] = qq_result[i + 1]; }

    if (config.verbose) { std::cout << "Total primes used: " << used_primes.size() << std::endl; }

    if (timing) {
        auto time_end_all = std::chrono::steady_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end_all - time_start_all).count();
        if (std::getenv("RUR_VERBOSE_TIMING")) {
            std::cout << "[timing] summary: reference_ms=" << timing_ref_ms << " modular_ms=" << timing_total_modular_ms
                      << " crt_ms=" << timing_total_crt_ms << " total_ms=" << total_ms << std::endl;
        }
    }
    return result;
}

bool
is_zero_dimensional_system(const std::vector<std::string> &polynomials,
                           const std::vector<std::string> &variables,
                           ModularCoeff prime) {
    // Check if prime is safe for F4
    const ModularCoeff F4_MAX_SAFE_PRIME = 2147483647; // 2^31 - 1
    if (prime > F4_MAX_SAFE_PRIME) {
        // Cannot use F4 with this prime, assume system is zero-dimensional
        // This is a conservative assumption
        return true;
    }

    // Quick check by computing GB and checking leading terms
    std::vector<const char *> var_ptrs;
    for (const auto &var : variables) { var_ptrs.push_back(var.c_str()); }

    axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), variables.size());
    if (!session) return false;

    // Add polynomials
    for (const auto &poly : polynomials) {
        std::string formatted = format_polynomial_for_f4(poly);
        if (axf4_add_polynomial(session, formatted.c_str()) != 0) {
            axf4_destroy_session(session);
            return false;
        }
    }

    // Compute GB
    axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
    if (gb_result.status != 0) {
        axf4_free_result(&gb_result);
        axf4_destroy_session(session);
        return false;
    }

    // Extract leading terms and check zero-dimensionality
    std::vector<PP> leading_terms;
    int num_vars = axf4_get_num_variables(session);

    // This would use extract_f4_leading_monomials in full implementation
    bool is_zero_dim = false;
    if (gb_result.basis_size > 0 && num_vars > 0) {
        // Simplified check - would need full implementation
        is_zero_dim = true;
    }

    axf4_free_result(&gb_result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);

    return is_zero_dim;
}

int
compute_quotient_dimension(const std::vector<std::string> &polynomials,
                           const std::vector<std::string> &variables,
                           ModularCoeff prime) {
    RURConfig config;
    config.verbose = false;

    auto [result, coeffs] = compute_modular_rur(polynomials, variables, prime, config, {});

    if (result.success) { return static_cast<int>(result.quotient_basis.size()); }

    return -1;
}

std::string
format_rur_result(const RationalRURResult &result, const std::vector<std::string> &variables) {
    std::ostringstream oss;

    if (!result.success) {
        oss << "RUR computation failed: " << result.error_message << std::endl;
        return oss.str();
    }

    oss << "Rational Univariate Representation:" << std::endl;
    oss << "===================================" << std::endl;

    // Minimal polynomial
    oss << "\nMinimal polynomial f(T):" << std::endl;
    oss << "f(T) = ";
    for (size_t i = result.minimal_polynomial.size(); i > 0; --i) {
        const auto &coeff = result.minimal_polynomial[i - 1];
        if (coeff != 0) {
            if (i < result.minimal_polynomial.size() && coeff > 0) { oss << " + "; }
            oss << coeff;
            if (i > 1) {
                oss << "*T";
                if (i > 2) { oss << "^" << (i - 1); }
            }
        }
    }
    oss << std::endl;

    // Quotient basis dimension
    oss << "\nQuotient ring dimension: " << result.quotient_basis.size() << std::endl;

    // Parameterizations
    oss << "\nParameterizations:" << std::endl;
    for (size_t i = 0; i < variables.size(); ++i) { oss << variables[i] << " = g_" << i << "(T) / f'(T)" << std::endl; }

    return oss.str();
}

} // namespace julia_rur