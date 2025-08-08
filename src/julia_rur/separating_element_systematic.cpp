#include "separating_element_systematic.hpp"
#include <iostream>

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

    while (true) {
        std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);
        for (int i = 0; i < num_variables; ++i) {
            ModularCoeff c_mod = (coeffs[i] < 0) ? (prime - ((-coeffs[i]) % prime)) % prime : coeffs[i] % prime;
            for (size_t j = 0; j < linear_form.size(); ++j) {
                linear_form[j] = (linear_form[j] + static_cast<AccModularCoeff>(c_mod) * var_vectors[i][j]) % prime;
            }
        }

        MinimalPolynomialResult test_mp =
          compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);

        if (test_mp.success && test_mp.degree == quotient_basis.size()) {
            // Pass integer coefficients instead of quotient-basis vector
            std::vector<int> coeffs_int(num_variables, 0);
            for (int i = 0; i < num_variables && i < static_cast<int>(coeffs.size()); ++i) { coeffs_int[i] = coeffs[i]; }
            auto [success, mp, params] =
              try_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, coeffs_int, -1);
            if (success) { return { true, mp, params, coeffs }; }
        }

        if (coeffs.empty()) { break; }
        coeffs = get_next_coeffs(coeffs, max_coeffs);
        if (coeffs.empty()) { break; }
    }

    return { false, {}, {}, {} };
}

} // namespace julia_rur
