#include "multiplication_tables.hpp"
#include "univariate_parameterization.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <set>


namespace julia_rur {

// Try to find a separating element and compute univariate parameterization

// Wrapper that retries with different separating elements until success
std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>>
compute_univariate_parameterization_with_retry(
  std::vector<PP> quotient_basis, // Changed to pass by value to allow modification
  const std::vector<std::vector<int32_t>> &i_xw,
  std::vector<std::vector<ModularCoeff>> &t_v,
  int num_variables,
  ModularCoeff prime,
  int max_retries) {
    // CRITICAL FIX: Ensure '1' (a monomial of all zeros) is the first element of the quotient basis.
    // The entire julia_rur library assumes quotient_basis[0] is the constant monomial.

    // Debug: Print current quotient basis
    std::cout << "DEBUG: Checking quotient basis ordering (size=" << quotient_basis.size() << "):" << std::endl;
    for (size_t i = 0; i < quotient_basis.size() && i < 3; ++i) {
        std::cout << "  quotient_basis[" << i << "] = [";
        for (size_t j = 0; j < quotient_basis[i].size(); ++j) {
            if (j > 0) std::cout << ",";
            std::cout << quotient_basis[i][j];
        }
        std::cout << "]" << std::endl;
    }

    auto it = std::find_if(quotient_basis.begin(), quotient_basis.end(), [](const PP &m) {
        return std::all_of(m.begin(), m.end(), [](int e) { return e == 0; });
    });

    if (it != quotient_basis.end() && it != quotient_basis.begin()) {
        std::cout << "INFO: Reordering quotient basis to place '1' at the front." << std::endl;
        std::cout << "  Before: quotient_basis[0] = [";
        for (size_t i = 0; i < quotient_basis[0].size(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << quotient_basis[0][i];
        }
        std::cout << "]" << std::endl;

        std::iter_swap(quotient_basis.begin(), it);

        std::cout << "  After: quotient_basis[0] = [";
        for (size_t i = 0; i < quotient_basis[0].size(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << quotient_basis[0][i];
        }
        std::cout << "]" << std::endl;
    } else if (it == quotient_basis.end()) {
        std::cout << "ERROR: No constant monomial '1' found in quotient basis!" << std::endl;
    }

    MinimalPolynomialResult min_poly;
    std::vector<BivariateResult> parameterizations(num_variables);

    // For single variable systems, just use the variable itself
    if (num_variables == 1) {
        std::cout << "Single variable system: using x1 as separating element" << std::endl;
        auto [success, mp, params] = try_separating_element(quotient_basis,
                                                            i_xw,
                                                            t_v,
                                                            num_variables,
                                                            prime,
                                                            std::vector<int>(),
                                                            1 // Use variable 1 as separating element
        );

        if (success) {
            std::cout << "SUCCESS: Single variable parameterization!" << std::endl;
            return { true, mp, params };
        } else {
            std::cerr << "Failed to parameterize single variable system!" << std::endl;
            return { false, min_poly, parameterizations };
        }
    }

    // For multivariate systems, first try individual variables
    for (int var = num_variables; var >= 1; --var) {
        std::cout << "\nTrying variable x" << var << " as separating element..." << std::endl;
        auto [success, mp, params] = try_separating_element(quotient_basis,
                                                            i_xw,
                                                            t_v,
                                                            num_variables,
                                                            prime,
                                                            std::vector<int>(),
                                                            var // Use variable var as separating element
        );

        if (success) {
            std::cout << "SUCCESS: Variable x" << var << " is a separating element!" << std::endl;
            return { true, mp, params };
        }
    }

    // Pre-calculate the vector representation for each variable from the multiplication tables.
    // This is the correct way to get these vectors, as they may not be basis elements themselves.
    std::vector<std::vector<ModularCoeff>> var_vectors(num_variables);
    for (int i = 0; i < num_variables; ++i) {
        if (i < static_cast<int>(i_xw.size()) && !i_xw[i].empty()) {
            int32_t table_idx = i_xw[i][0]; // Entry for var_(i+1) * 1
            if (table_idx > 0 && table_idx <= static_cast<int32_t>(t_v.size())) {
                var_vectors[i] = t_v[table_idx - 1];
            } else {
                // This case should ideally not happen in a valid setup
                var_vectors[i].assign(quotient_basis.size(), 0);
            }
        } else {
            var_vectors[i].assign(quotient_basis.size(), 0);
        }
    }

    // Try structured linear combinations first
    std::cout << "\nTrying structured linear combinations...\n";

    // For 2 variables, try systematic small coefficients
    if (num_variables == 2) {
        // First try simple combinations that often work
        std::vector<std::pair<int, int>> structured_coeffs = { { 1, 1 },  { 1, 2 },  { 2, 1 }, { 1, -1 },
                                                               { 2, -1 }, { -1, 2 }, { 3, 1 }, { 1, 3 },
                                                               { 3, -1 }, { -1, 3 }, { 1, 0 }, { 0, 1 } };

        // For tensor product systems (like x³=1, y²=1), we need to mix both factors
        // Try more combinations systematically
        for (int a = -5; a <= 5; ++a) {
            if (a == 0) continue; // Skip a=0
            for (int b = -5; b <= 5; ++b) {
                // Skip if both are 0 or if we already tried this in structured_coeffs
                if (b == 0 && std::abs(a) <= 1) continue;

                structured_coeffs.push_back({ a, b });
            }
        }

        for (const auto &[c1, c2] : structured_coeffs) {
            std::cout << "Trying " << c1 << "*x1 + " << c2 << "*x2..." << std::endl;
            std::vector<int> coeffs = { c1, c2 };

            // Check if it's separating
            MinimalPolynomialResult test_mp =
              compute_minimal_polynomial_flint(coeffs, i_xw, t_v, quotient_basis, prime);

            if (test_mp.success) {
                if (test_mp.degree == quotient_basis.size()) {
                    std::cout << "  -> Minimal polynomial degree " << test_mp.degree << " matches quotient size!"
                              << std::endl;
                    auto [success, mp, params] =
                      try_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, coeffs, -1);

                    if (success) {
                        std::cout << "SUCCESS: Found separating linear combination " << c1 << "*x1 + " << c2 << "*x2!"
                                  << std::endl;
                        return { true, mp, params };
                    }
                } else if ((c1 == 1 && c2 == 1) || (std::abs(c1) <= 3 && std::abs(c2) <= 3)) {
                    // Only print degree mismatch for simple combinations to avoid spam
                    std::cout << "  -> Minimal polynomial degree " << test_mp.degree << " != quotient size "
                              << quotient_basis.size() << std::endl;
                }
            }
        }
    }

    // Try systematic approach similar to Julia's :current strategy
    std::cout << "\nTrying systematic separating element search (Julia-style)...\n";

    // Julia's strategy: Start with [0, 0, ..., -1, 1] and systematically modify
    std::vector<int> sep_coeffs(num_variables, 0);
    if (num_variables >= 2) {
        sep_coeffs[num_variables - 1] = 1;  // Last variable
        sep_coeffs[num_variables - 2] = -1; // Second to last

        int vv = num_variables - 2;
        int max_attempts = 100; // Safety limit

        for (int attempt = 0; attempt < max_attempts; ++attempt) {
            std::cout << "  Trying separating form: [";
            for (int i = 0; i < num_variables; ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << sep_coeffs[i];
            }
            std::cout << "]" << std::endl;

            // Use ONLY FLINT implementation (OLD method disabled - produces wrong coefficients)
            MinimalPolynomialResult test_mp =
              compute_minimal_polynomial_flint(sep_coeffs, i_xw, t_v, quotient_basis, prime);

            if (test_mp.success && test_mp.degree == quotient_basis.size()) {
                std::cout << "  SUCCESS: Found separating element with coeffs [";
                for (int i = 0; i < num_variables; ++i) {
                    if (i > 0) std::cout << ", ";
                    std::cout << sep_coeffs[i];
                }
                std::cout << "]" << std::endl;

                auto [success, mp, params] =
                  try_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, sep_coeffs, -1);
                if (success) { return { true, mp, params }; }
            }

            // Julia's strategy for updating coefficients
            // If the minimal polynomial has the wrong degree, find which variable caused it
            // and modify that coefficient
            int uu = test_mp.degree - 1; // This is a simplified heuristic
            if (uu < 0) uu = 0;
            if (uu >= num_variables) uu = num_variables - 1;

            // Update the coefficient
            if (sep_coeffs[uu] < 0) {
                sep_coeffs[uu] = -sep_coeffs[uu];
            } else {
                sep_coeffs[uu] = -sep_coeffs[uu] - 1;
            }
            vv = uu;

            if (vv < 0) {
                std::cout << "  Systematic search failed, moving to random search" << std::endl;
                break;
            }
        }
    }

    // Then try random linear forms
    std::cout << "\nTrying random linear combinations...\n";
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int attempt = 0; attempt < max_retries; ++attempt) {
        std::cout << "\nAttempt " << (attempt + 1) << ": Trying random linear form as separating element..."
                  << std::endl;

        auto linear_form_coeffs = compute_random_linear_form(quotient_basis, i_xw, t_v, num_variables, prime, gen);

        // Print the linear form coefficients
        std::cout << "  Linear form coefficients: [";
        for (size_t i = 0; i < linear_form_coeffs.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << linear_form_coeffs[i];
        }
        std::cout << "] mod " << prime << std::endl;

        // Check if it's potentially separating
        MinimalPolynomialResult test_mp =
          compute_minimal_polynomial_flint(linear_form_coeffs, i_xw, t_v, quotient_basis, prime);
        if (!test_mp.success || test_mp.degree != quotient_basis.size()) {
            std::cout << "  Random linear form is not separating (minimal poly degree " << test_mp.degree
                      << " != quotient size " << quotient_basis.size() << ")" << std::endl;
            std::cout << "  Minpoly coefficients: [";
            for (size_t i = 0; i < test_mp.coefficients.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << test_mp.coefficients[i];
            }
            std::cout << "]" << std::endl;
            continue;
        }

        auto [success, mp, params] =
          try_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, linear_form_coeffs, -1);

        if (success) {
            std::cout << "SUCCESS: Found valid separating linear form with all linear parameterizations!" << std::endl;
            return { true, mp, params };
        }

        std::cout << "Random linear form failed (non-linear parameterization detected)" << std::endl;
    }

    std::cerr << "Failed to find separating element after " << max_retries << " attempts" << std::endl;

    // Final attempt: deterministic search
    std::cout << "\nTrying deterministic linear combinations as a last resort...\n";
    auto [sep_coeffs_det, sep_var_det] =
      find_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, SeparatingStrategy::DETERMINISTIC);

    if (!sep_coeffs_det.empty()) {
        return try_separating_element(quotient_basis, i_xw, t_v, num_variables, prime, sep_coeffs_det, sep_var_det);
    }


    return { false, min_poly, parameterizations };
}

} // namespace julia_rur