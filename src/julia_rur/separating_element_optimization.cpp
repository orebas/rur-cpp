/**
 * Separating Element Optimization
 *
 * This file implements the optimization for reusing successful separating elements
 * across different primes in the multi-modular RUR algorithm.
 *
 * When we find a separating element that works for the reference prime,
 * we should try those same coefficients (modulo the new prime) as our first
 * attempt for subsequent primes.
 */

#include "multiplication_tables.hpp"
#include "univariate_parameterization.hpp"
#include <iostream>
#include <vector>

namespace julia_rur {

// Structure to store successful separating element information
struct SeparatingElementInfo {
    bool is_variable;                    // true if a single variable, false if linear form
    int variable_index;                  // if is_variable, which variable (1-based)
    std::vector<int> linear_form_coeffs; // if !is_variable, coefficients for linear form

    SeparatingElementInfo()
      : is_variable(false)
      , variable_index(-1) {}
};

/**
 * Try a separating element with given coefficients (modulo prime)
 * Returns true if it successfully separates
 */
bool
try_separating_coefficients(const std::vector<int> &coeffs,
                            const std::vector<PP> &quotient_basis,
                            const std::vector<std::vector<int32_t>> &i_xw,
                            const std::vector<std::vector<ModularCoeff>> &t_v,
                            int num_variables,
                            ModularCoeff prime) {
    std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);

    // Build the linear form from coefficients
    for (int i = 0; i < num_variables && i < static_cast<int>(coeffs.size()); ++i) {
        if (coeffs[i] == 0) continue;

        // Get the vector representation of variable i+1
        std::vector<ModularCoeff> var_vec = element_to_vector(i + 1, i_xw, t_v, quotient_basis.size());

        // Convert coefficient to modular form
        ModularCoeff lambda_mod = (coeffs[i] < 0) ? (prime - ((-coeffs[i]) % prime)) % prime : coeffs[i] % prime;

        // Add to linear form
        for (size_t j = 0; j < linear_form.size(); ++j) {
            AccModularCoeff prod = static_cast<AccModularCoeff>(lambda_mod) * var_vec[j];
            linear_form[j] = (linear_form[j] + prod % prime) % prime;
        }
    }

    // Check if it's separating
    MinimalPolynomialResult test_mp = compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);

    return test_mp.success && (test_mp.degree == quotient_basis.size() ||
                               (test_mp.original_degree > 0 && test_mp.original_degree == quotient_basis.size()));
}

/**
 * Enhanced compute_univariate_parameterization that tries a hint first
 */
std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>>
compute_univariate_parameterization_with_hint(std::vector<PP> quotient_basis, // Pass by value to allow modification
                                              const std::vector<std::vector<int32_t>> &i_xw,
                                              std::vector<std::vector<ModularCoeff>> &t_v,
                                              int num_variables,
                                              ModularCoeff prime,
                                              const SeparatingElementInfo &hint,
                                              int max_retries) {
    // First ensure quotient basis ordering (same as before)
    auto it = std::find_if(quotient_basis.begin(), quotient_basis.end(), [](const PP &m) {
        return std::all_of(m.begin(), m.end(), [](int e) { return e == 0; });
    });

    if (it != quotient_basis.end() && it != quotient_basis.begin()) {
        std::cout << "INFO: Reordering quotient basis to place '1' at the front." << std::endl;
        std::iter_swap(quotient_basis.begin(), it);
    }

    MinimalPolynomialResult min_poly;
    std::vector<BivariateResult> parameterizations(num_variables);

    // Try the hint first if provided
    if (hint.is_variable && hint.variable_index > 0) {
        std::cout << "\\nTrying hint: variable x" << hint.variable_index << " as separating element..." << std::endl;

        auto [success, mp, params] = try_separating_element(
          quotient_basis, i_xw, t_v, num_variables, prime, std::vector<ModularCoeff>(), hint.variable_index);

        if (success) {
            std::cout << "SUCCESS: Hint worked! Variable x" << hint.variable_index << " is separating for prime "
                      << prime << std::endl;
            return { true, mp, params };
        }
    } else if (!hint.is_variable && !hint.linear_form_coeffs.empty()) {
        std::cout << "\\nTrying hint: linear form [";
        for (size_t i = 0; i < hint.linear_form_coeffs.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << hint.linear_form_coeffs[i];
        }
        std::cout << "] as separating element..." << std::endl;

        // Build linear form from hint coefficients
        std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);
        for (int i = 0; i < num_variables && i < static_cast<int>(hint.linear_form_coeffs.size()); ++i) {
            if (hint.linear_form_coeffs[i] == 0) continue;

            std::vector<ModularCoeff> var_vec = element_to_vector(i + 1, i_xw, t_v, quotient_basis.size());
            ModularCoeff lambda_mod = (hint.linear_form_coeffs[i] < 0)
                                        ? (prime - ((-hint.linear_form_coeffs[i]) % prime)) % prime
                                        : hint.linear_form_coeffs[i] % prime;

            for (size_t j = 0; j < linear_form.size(); ++j) {
                AccModularCoeff prod = static_cast<AccModularCoeff>(lambda_mod) * var_vec[j];
                linear_form[j] = (linear_form[j] + prod % prime) % prime;
            }
        }

        // Check if it's separating
        MinimalPolynomialResult test_mp =
          compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);

        if (test_mp.success && test_mp.degree == quotient_basis.size()) {
            std::cout << "SUCCESS: Hint worked! Linear form is separating for prime " << prime << std::endl;

            auto [success, mp, params] =
              try_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, linear_form, -1);

            if (success) { return { true, mp, params }; }
        } else {
            std::cout << "Hint didn't work (degree " << test_mp.degree << " != " << quotient_basis.size()
                      << "), trying other methods..." << std::endl;
        }
    }

    // If hint didn't work or wasn't provided, fall back to the original retry logic
    return compute_univariate_parameterization_with_retry(quotient_basis, i_xw, t_v, num_variables, prime, max_retries);
}

/**
 * Extract separating element info from a successful computation
 */
SeparatingElementInfo
extract_separating_element_info(const MinimalPolynomialResult &min_poly,
                                const std::vector<BivariateResult> &params,
                                int num_variables,
                                int sep_var,
                                const std::vector<int> &sep_coeffs) {
    SeparatingElementInfo info;

    if (sep_var > 0) {
        info.is_variable = true;
        info.variable_index = sep_var;
    } else if (!sep_coeffs.empty()) {
        info.is_variable = false;
        info.linear_form_coeffs = sep_coeffs;
    }

    return info;
}

} // namespace julia_rur