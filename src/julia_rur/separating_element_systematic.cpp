#include "separating_element_systematic.hpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>

namespace julia_rur {

static std::vector<int>
get_next_coeffs(const std::vector<int> &current_coeffs, int max_c) {
    std::vector<int> next_coeffs = current_coeffs;
    int n = next_coeffs.size();
    int carry = 1;

    for (int i = 0; i < n; ++i) {
        if (carry == 0) break;
        int val = next_coeffs[i] + carry;
        if (val > max_c) {
            next_coeffs[i] = -max_c;
            carry = 1;
        } else {
            next_coeffs[i] = val;
            carry = 0;
        }
    }

    if (carry) { // overflow
        return {};
    }
    return next_coeffs;
}

std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>, std::vector<int>>
find_separating_element_systematic(std::vector<PP> quotient_basis,
                                   const std::vector<std::vector<int32_t>> &i_xw,
                                   std::vector<std::vector<ModularCoeff>> &t_v,
                                   int num_variables,
                                   ModularCoeff prime,
                                   int max_coeffs) {
    bool timing = (std::getenv("RUR_TIMING") && std::string(std::getenv("RUR_TIMING")) != "0");
    auto t_start_search = std::chrono::steady_clock::now();
    long long candidates_tested = 0;
    long long ms_minpoly_accum = 0;
    long long ms_try_accum = 0;

    std::cout << "[DEBUG] Starting systematic separating element search:" << std::endl;
    std::cout << "  quotient_basis.size() = " << quotient_basis.size() << std::endl;
    std::cout << "  num_variables = " << num_variables << std::endl;
    std::cout << "  prime = " << prime << std::endl;
    std::cout << "  max_coeffs = " << max_coeffs << std::endl;

    // Cap total candidates per call to keep runtime predictable
    long long max_candidates = 1000;
    if (const char *env = std::getenv("RUR_SEP_MAX_ATTEMPTS")) {
        try {
            long long v = std::stoll(env);
            if (v > 0) max_candidates = v;
        } catch (...) {}
    }

    std::vector<std::vector<ModularCoeff>> var_vectors(num_variables);
    for (int i = 0; i < num_variables; ++i) {
        if (i < static_cast<int>(i_xw.size()) && !i_xw[i].empty()) {
            int32_t table_idx = i_xw[i][0];
            if (table_idx > 0 && table_idx <= static_cast<int32_t>(t_v.size())) {
                var_vectors[i] = t_v[table_idx - 1];
            } else {
                var_vectors[i].assign(quotient_basis.size(), 0);
            }
        } else {
            var_vectors[i].assign(quotient_basis.size(), 0);
        }
    }

    std::vector<int> coeffs(num_variables, -max_coeffs);
    if (!coeffs.empty()) {
        coeffs[0] = 1; // Start with x1
    }

    // Deterministic randomized pre-pass: sample a few random linear forms
    // to quickly escape hostile structures; seed by prime for reproducibility
    {
        std::mt19937 gen(static_cast<uint32_t>(prime));
        std::uniform_int_distribution<int> dis(-max_coeffs, max_coeffs);
        int rand_trials = static_cast<int>(std::min<long long>(max_candidates / 4, 200));
        for (int t = 0; t < rand_trials && candidates_tested < max_candidates; ++t) {
            std::vector<int> rc(num_variables, 0);
            bool all_zero = true;
            for (int i = 0; i < num_variables; ++i) {
                rc[i] = dis(gen);
                if (rc[i] != 0) all_zero = false;
            }
            if (all_zero) { rc[0] = 1; }

            std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);
            for (int i = 0; i < num_variables; ++i) {
                ModularCoeff c_mod = (rc[i] < 0) ? (prime - ((-rc[i]) % prime)) % prime : rc[i] % prime;
                for (size_t j = 0; j < linear_form.size(); ++j) {
                    linear_form[j] = (linear_form[j] + static_cast<AccModularCoeff>(c_mod) * var_vectors[i][j]) % prime;
                }
            }

            auto t_start_mp = std::chrono::steady_clock::now();
            MinimalPolynomialResult test_mp =
              compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);
            auto t_end_mp = std::chrono::steady_clock::now();
            candidates_tested++;
            ms_minpoly_accum += std::chrono::duration_cast<std::chrono::milliseconds>(t_end_mp - t_start_mp).count();

            // Accept if degree matches OR if we have a non-radical ideal after square-free
            bool is_separating = (test_mp.degree == quotient_basis.size()) ||
                                 (test_mp.original_degree > 0 && test_mp.original_degree == quotient_basis.size());

            if (test_mp.success && is_separating) {
                auto t_start_try = std::chrono::steady_clock::now();
                auto [success, mp, params] =
                  try_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, rc, -1);
                auto t_end_try = std::chrono::steady_clock::now();
                ms_try_accum += std::chrono::duration_cast<std::chrono::milliseconds>(t_end_try - t_start_try).count();
                if (success) {
                    if (timing) {
                        auto t_end_search = std::chrono::steady_clock::now();
                        auto ms_total =
                          std::chrono::duration_cast<std::chrono::milliseconds>(t_end_search - t_start_search).count();
                        std::cout << "[timing.search] candidates=" << candidates_tested
                                  << " minpoly_ms=" << ms_minpoly_accum << " try_ms=" << ms_try_accum
                                  << " total_ms=" << ms_total << std::endl;
                    }
                    return { true, mp, params, rc };
                }
            }
        }
    }

    while (true) {
        std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);
        for (int i = 0; i < num_variables; ++i) {
            ModularCoeff c_mod = (coeffs[i] < 0) ? (prime - ((-coeffs[i]) % prime)) % prime : coeffs[i] % prime;
            for (size_t j = 0; j < linear_form.size(); ++j) {
                linear_form[j] = (linear_form[j] + static_cast<AccModularCoeff>(c_mod) * var_vectors[i][j]) % prime;
            }
        }

        auto t_start_mp = std::chrono::steady_clock::now();
        MinimalPolynomialResult test_mp =
          compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);
        auto t_end_mp = std::chrono::steady_clock::now();
        candidates_tested++;
        ms_minpoly_accum += std::chrono::duration_cast<std::chrono::milliseconds>(t_end_mp - t_start_mp).count();

        // Accept if degree matches OR if we have a non-radical ideal after square-free
        bool is_separating = (test_mp.degree == quotient_basis.size()) ||
                             (test_mp.original_degree > 0 && test_mp.original_degree == quotient_basis.size());

        if (test_mp.success && is_separating) {
            // Pass integer coefficients instead of quotient-basis vector
            std::vector<int> coeffs_int(num_variables, 0);
            for (int i = 0; i < num_variables && i < static_cast<int>(coeffs.size()); ++i) {
                coeffs_int[i] = coeffs[i];
            }
            auto t_start_try = std::chrono::steady_clock::now();
            auto [success, mp, params] =
              try_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, coeffs_int, -1);
            auto t_end_try = std::chrono::steady_clock::now();
            ms_try_accum += std::chrono::duration_cast<std::chrono::milliseconds>(t_end_try - t_start_try).count();
            if (success) {
                if (timing) {
                    auto t_end_search = std::chrono::steady_clock::now();
                    auto ms_total =
                      std::chrono::duration_cast<std::chrono::milliseconds>(t_end_search - t_start_search).count();
                    std::cout << "[timing.search] candidates=" << candidates_tested
                              << " minpoly_ms=" << ms_minpoly_accum << " try_ms=" << ms_try_accum
                              << " total_ms=" << ms_total << std::endl;
                }
                return { true, mp, params, coeffs };
            }
        }

        if (coeffs.empty()) { break; }
        coeffs = get_next_coeffs(coeffs, max_coeffs);
        if (coeffs.empty()) { break; }
        if (candidates_tested >= max_candidates) { break; }
    }

    if (timing) {
        auto t_end_search = std::chrono::steady_clock::now();
        auto ms_total = std::chrono::duration_cast<std::chrono::milliseconds>(t_end_search - t_start_search).count();
        std::cout << "[timing.search] candidates=" << candidates_tested << " minpoly_ms=" << ms_minpoly_accum
                  << " try_ms=" << ms_try_accum << " total_ms=" << ms_total << std::endl;
    }
    std::cout << "[DEBUG] Systematic search FAILED after testing " << candidates_tested << " candidates" << std::endl;
    return { false, {}, {}, {} };
}

} // namespace julia_rur
