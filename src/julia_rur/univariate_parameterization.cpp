#include "univariate_parameterization.hpp"
#include "bivariate_algorithm.hpp"
#include "flint_linear_algebra.hpp" // For FLINT-based implementation
#include "multiplication_tables.hpp"
#include "polynomial_operations.hpp"
#include <algorithm>
#include <iostream>
#include <random>

extern "C" {
#include <flint/nmod_poly.h>
}

namespace julia_rur {

// Helper function to multiply two elements in the quotient ring
// This uses the multiplication tables to compute element1 * element2
std::vector<ModularCoeff>
multiply_elements_in_quotient(const std::vector<ModularCoeff> &element1,
                              const std::vector<ModularCoeff> &element2,
                              const std::vector<PP> &quotient_basis,
                              const std::vector<std::vector<int32_t>> &i_xw,
                              const std::vector<std::vector<ModularCoeff>> &t_v,
                              ModularCoeff prime) {
    size_t quotient_basis_size = quotient_basis.size();
    if (quotient_basis_size == 0) { return std::vector<ModularCoeff>(); }
    if (element1.size() != quotient_basis_size || element2.size() != quotient_basis_size) {
        // Size mismatch - elements must have same size as quotient basis
        return std::vector<ModularCoeff>();
    }
    size_t n_vars = quotient_basis[0].size();
    std::vector<ModularCoeff> result(quotient_basis_size, 0);

    // Pre-allocate buffers to reduce allocations in inner loops
    std::vector<ModularCoeff> temp;
    temp.reserve(quotient_basis_size);
    std::vector<ModularCoeff> next(quotient_basis_size, 0);

    // For each monomial in element1
    for (size_t i = 0; i < quotient_basis_size; ++i) {
        if (element1[i] == 0) continue;

        // Get the monomial from quotient basis
        const PP &monomial = quotient_basis[i];

        // Start with element2 scaled by coefficient
        if (element2.size() != quotient_basis_size) {
            // Size mismatch - can't multiply
            return std::vector<ModularCoeff>();
        }
        temp.assign(element2.begin(), element2.end());
        for (size_t j = 0; j < quotient_basis_size; ++j) {
            temp[j] = (static_cast<AccModularCoeff>(temp[j]) * element1[i]) % prime;
        }

        // Multiply by each variable in the monomial
        for (size_t var_idx = 0; var_idx < n_vars; ++var_idx) {
            for (int deg = 0; deg < monomial[var_idx]; ++deg) {
                std::fill(next.begin(), next.end(), 0);
                mul_var_quo(next, temp, var_idx + 1, i_xw, t_v, prime);
                temp.swap(next);
            }
        }

        // Add to result
        for (size_t j = 0; j < quotient_basis_size; ++j) { result[j] = (result[j] + temp[j]) % prime; }
    }

    return result;
}

// Helper to compute random linear form
std::vector<int>
compute_random_linear_form(const std::vector<PP> &quotient_basis,
                           const std::vector<std::vector<int32_t>> &i_xw,
                           const std::vector<std::vector<ModularCoeff>> &t_v,
                           int num_variables,
                           ModularCoeff prime,
                           std::mt19937 &gen) {
    std::uniform_int_distribution<> dis(-20, 20);
    std::vector<int> coefficients(num_variables, 0);

    // Keep trying until we get a non-zero linear form
    bool all_zero = true;
    while (all_zero) {
        for (int i = 0; i < num_variables; ++i) {
            coefficients[i] = dis(gen);
            if (coefficients[i] != 0) { all_zero = false; }
        }
    }
    return coefficients;
}

// Helper function to compute element representation in quotient basis
std::vector<ModularCoeff>
element_to_vector(int32_t var_index,
                  const std::vector<std::vector<int32_t>> &i_xw,
                  const std::vector<std::vector<ModularCoeff>> &t_v,
                  size_t quotient_basis_size) {
    // Get the vector for variable 'var_index' from the multiplication table.
    // This is the correct way to get the representation of a variable,
    // as it may not be a basis element itself.

    // Ensure var_index is 1-based and valid
    if (var_index < 1 || var_index > static_cast<int32_t>(i_xw.size()) || i_xw[var_index - 1].empty()) {
        return std::vector<ModularCoeff>(quotient_basis_size, 0);
    }

    // Get the index for `var * 1` from the multiplication table map.
    // The `[0]` corresponds to multiplication by the first basis element, which is `1`.
    int32_t table_idx = i_xw[var_index - 1][0];

    if (table_idx <= 0 || table_idx > static_cast<int32_t>(t_v.size())) {
        return std::vector<ModularCoeff>(quotient_basis_size, 0);
    }

    // The multiplication table is 1-indexed, so subtract 1.
    return t_v[table_idx - 1];
}

// Helper function for Gaussian reduction to detect linear dependence
static bool
gaussian_reduce(std::vector<std::vector<ModularCoeff>> &matrix,
                std::vector<ModularCoeff> &new_row,
                ModularCoeff prime,
                std::vector<ModularCoeff> &reduction_coeffs) {
    size_t num_cols = new_row.size();
    reduction_coeffs.assign(matrix.size(), 0);

    for (size_t j = 0; j < num_cols; ++j) {
        if (new_row[j] == 0) continue;

        bool found_pivot = false;
        for (size_t i = 0; i < matrix.size(); ++i) {
            if (matrix[i][j] != 0) {
                ModularCoeff factor = new_row[j];
                ModularCoeff pivot = matrix[i][j];

                ModularCoeff pivot_inv = 1;
                ModularCoeff temp = pivot;
                ModularCoeff exp = prime - 2;
                while (exp > 0) {
                    if (exp & 1) pivot_inv = (static_cast<AccModularCoeff>(pivot_inv) * temp) % prime;
                    temp = (static_cast<AccModularCoeff>(temp) * temp) % prime;
                    exp >>= 1;
                }

                ModularCoeff scale = (static_cast<AccModularCoeff>(factor) * pivot_inv) % prime;
                reduction_coeffs[i] = (reduction_coeffs[i] + scale) % prime;

                for (size_t c = 0; c < num_cols; ++c) {
                    AccModularCoeff prod = static_cast<AccModularCoeff>(scale) * matrix[i][c];
                    new_row[c] = (new_row[c] + prime - (prod % prime)) % prime;
                }
            }
        }
    }

    for (ModularCoeff val : new_row) {
        if (val != 0) return false;
    }
    return true;
}

/*MinimalPolynomialResult
compute_minimal_polynomial(const std::vector<ModularCoeff> &element,
                           const std::vector<std::vector<int32_t>> &i_xw,
                           const std::vector<std::vector<ModularCoeff>> &t_v,
                           const std::vector<PP> &quotient_basis,
                           ModularCoeff prime) {
    MinimalPolynomialResult result;
    result.success = false;

    size_t quotient_basis_size = quotient_basis.size();

    // Matrix to store powers of element
    std::vector<std::vector<ModularCoeff>> power_matrix;
    std::vector<ModularCoeff> current_power = element;

    // Add T^0 = 1
    std::vector<ModularCoeff> identity(quotient_basis_size, 0);
    identity[0] = 1; // Assuming first element in quotient basis is 1
    power_matrix.push_back(identity);
    result.powers.push_back(identity);

    // Maximum degree is the dimension of quotient ring
    size_t max_degree = quotient_basis_size;

    for (size_t deg = 1; deg <= max_degree; ++deg) {
        // Check if current_power is linearly dependent on previous powers
        std::vector<ModularCoeff> test_row = current_power;

        std::vector<ModularCoeff> reduction_coeffs;
        if (gaussian_reduce(power_matrix, test_row, prime, reduction_coeffs)) {
            // Found linear dependence! Extract minimal polynomial
            result.degree = deg;
            result.coefficients.resize(deg + 1, 0);

            // The last power is a linear combination of previous powers
            // T^d = c_{d-1}T^{d-1} + ... + c_0
            // Minimal polynomial: T^d - c_{d-1}T^{d-1} - ... - c_0 = 0
            result.coefficients[deg] = 1;
            for (size_t i = 0; i < deg; ++i) { result.coefficients[i] = (prime - reduction_coeffs[i]) % prime; }

            result.success = true;
            return result;
        }

        // Add current power to matrix and powers list
        power_matrix.push_back(current_power);
        result.powers.push_back(current_power);

        if (deg < max_degree) {
            // Compute next power: multiply current by element
            // Use proper multiplication in the quotient ring
            std::vector<ModularCoeff> next_power =
              multiply_elements_in_quotient(current_power, element, quotient_basis, i_xw, t_v, prime);

            current_power = next_power;
        }
    }

    // Should not reach here for zero-dimensional ideals
    return result;
}*/ // OLD CODE FROM CLAUDE

MinimalPolynomialResult
compute_minimal_polynomial_OLD_DISABLED(const std::vector<ModularCoeff> &element,
                                        const std::vector<std::vector<int32_t>> &i_xw,
                                        const std::vector<std::vector<ModularCoeff>> &t_v,
                                        const std::vector<PP> &quotient_basis,
                                        ModularCoeff prime) {
    MinimalPolynomialResult result;
    result.success = false;
    size_t d = quotient_basis.size();

    std::cout << "\n=== compute_minimal_polynomial DEBUG ===" << std::endl;
    std::cout << "Quotient basis size d = " << d << std::endl;
    std::cout << "Element vector: [";
    for (size_t i = 0; i < element.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << element[i];
    }
    std::cout << "]" << std::endl;

    // Special case for 1-dimensional quotient ring (e.g., simple linear equation)
    if (d == 1) {
        // Any element in this ring is just a constant. The 'element' vector
        // should have size 1 and contain this constant value 'c'.
        // The minimal polynomial is T - c.
        std::cout << "DEBUG: 1-dimensional quotient ring case, element size = " << element.size() << std::endl;
        if (element.size() != 1) {
            // This indicates a critical logic error elsewhere if it ever happens.
            std::cout << "ERROR: Element size mismatch for 1-dim case!" << std::endl;
            return result;
        }
        ModularCoeff c = element[0];
        result.coefficients = { (prime - c) % prime, 1 }; // Represents T - c
        result.degree = 1;
        result.powers.push_back(element); // T^0 is the element itself
        result.success = true;
        std::cout << "DEBUG: Returning minimal polynomial of degree " << result.degree << " for 1-dim case"
                  << std::endl;
        return result;
    }

    if (d == 0) return result;

    std::vector<std::vector<ModularCoeff>> matrix;
    std::vector<ModularCoeff> current_power(d, 0);
    current_power[0] = 1; // T^0 = 1

    for (size_t k = 0; k <= d; ++k) {
        std::cout << "\nIteration k=" << k << ", computing T^" << k << std::endl;
        std::cout << "  current_power = [";
        for (size_t i = 0; i < current_power.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << current_power[i];
        }
        std::cout << "]" << std::endl;

        std::vector<ModularCoeff> temp_row = current_power;

        // Reduce temp_row by the current basis in matrix
        std::cout << "  Matrix has " << matrix.size() << " rows, reducing temp_row..." << std::endl;
        for (const auto &row : matrix) {
            if (temp_row.empty()) break;

            size_t pivot_idx = 0;
            while (pivot_idx < d && row[pivot_idx] == 0) { pivot_idx++; }
            if (pivot_idx >= d) continue;

            ModularCoeff factor = temp_row[pivot_idx];
            if (factor == 0) continue;

            ModularCoeff inv_pivot = modular_inverse(row[pivot_idx], prime);
            ModularCoeff scale = (static_cast<AccModularCoeff>(factor) * inv_pivot) % prime;

            for (size_t i = 0; i < d; ++i) {
                AccModularCoeff prod = static_cast<AccModularCoeff>(scale) * row[i];
                temp_row[i] = (temp_row[i] + prime - (prod % prime)) % prime;
            }
        }

        bool is_zero = true;
        for (const auto &val : temp_row) {
            if (val != 0) {
                is_zero = false;
                break;
            }
        }

        if (is_zero) { // Linearly dependent
            std::cout << "  temp_row reduced to zero! Found linear dependence at degree " << k << std::endl;
            result.success = true;
            result.degree = k;

            // To find coefficients, solve matrix * c = -current_power
            // This is complex, so for now, we just indicate success
            result.coefficients.resize(k + 1);
            result.coefficients[k] = 1; // Monic

            // Back-substitution to find other coeffs
            std::vector<ModularCoeff> rhs = current_power;
            for (size_t i = 0; i < rhs.size(); ++i) { rhs[i] = (prime - rhs[i]) % prime; }

            for (int i = matrix.size() - 1; i >= 0; --i) {
                size_t pivot = 0;
                while (pivot < d && matrix[i][pivot] == 0) { pivot++; }

                // Check if pivot is valid before accessing rhs
                if (pivot >= d || pivot >= rhs.size()) {
                    continue; // Skip this row if no valid pivot found
                }

                ModularCoeff val = rhs[pivot];
                for (size_t j = i + 1; j < matrix.size(); ++j) {
                    if (pivot < result.coefficients.size() && j < d) {
                        AccModularCoeff prod = static_cast<AccModularCoeff>(matrix[i][j]) * result.coefficients[j];
                        val = (val + prime - (prod % prime)) % prime;
                    }
                }
                if (pivot < result.coefficients.size()) { result.coefficients[i] = val; }
            }

            return result;
        }

        std::cout << "  temp_row after reduction = [";
        for (size_t i = 0; i < temp_row.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << temp_row[i];
        }
        std::cout << "]" << std::endl;

        // Add to basis and perform Gaussian elimination
        size_t pivot_col = 0;
        while (pivot_col < d && temp_row[pivot_col] == 0) { pivot_col++; }

        ModularCoeff inv = modular_inverse(temp_row[pivot_col], prime);
        for (size_t i = 0; i < d; ++i) { temp_row[i] = (static_cast<AccModularCoeff>(temp_row[i]) * inv) % prime; }

        for (auto &row : matrix) {
            ModularCoeff factor = row[pivot_col];
            if (factor == 0) continue;
            for (size_t i = 0; i < d; ++i) {
                AccModularCoeff prod = static_cast<AccModularCoeff>(factor) * temp_row[i];
                row[i] = (row[i] + prime - (prod % prime)) % prime;
            }
        }

        matrix.push_back(temp_row);
        result.powers.push_back(current_power);

        if (k < d) {
            std::cout << "  Computing next power: [";
            for (size_t i = 0; i < current_power.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << current_power[i];
            }
            std::cout << "] * [";
            for (size_t i = 0; i < element.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << element[i];
            }
            std::cout << "]" << std::endl;

            current_power = multiply_elements_in_quotient(current_power, element, quotient_basis, i_xw, t_v, prime);

            std::cout << "    Result = [";
            for (size_t i = 0; i < current_power.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << current_power[i];
            }
            std::cout << "]" << std::endl;
        }
    }

    return result;
}


MinimalPolynomialResult
first_variable(int32_t separating_var_index,
               const std::vector<std::vector<int32_t>> &i_xw,
               std::vector<std::vector<ModularCoeff>> &t_v,
               const std::vector<PP> &quotient_basis,
               ModularCoeff prime) {
    // This function finds the vector representation of the separating variable
    // and then calls compute_minimal_polynomial to do the actual work.

    // Find the vector for the separating variable from the multiplication table.
    // This is the representation of the variable in the quotient basis.
    if (separating_var_index < 1 || separating_var_index > static_cast<int32_t>(i_xw.size()) ||
        i_xw[separating_var_index - 1].empty()) {
        std::cerr << "Error in first_variable: Invalid separating variable index " << separating_var_index << std::endl;
        MinimalPolynomialResult res;
        res.success = false;
        return res;
    }

    int32_t table_idx = i_xw[separating_var_index - 1][0]; // Entry for var * 1
    const bool verbose_uv = false;  // Disabled for clean output
    if (verbose_uv) {
        std::cout << "[DEBUG] first_variable: separating_var_index=" << separating_var_index
                  << ", i_xw size=" << i_xw.size() << ", table_idx=" << table_idx << ", t_v size=" << t_v.size()
                  << std::endl;
    }

    if (table_idx <= 0 || table_idx > static_cast<int32_t>(t_v.size())) {
        std::cerr << "Error in first_variable: Invalid table index " << table_idx << " from i_xw." << std::endl;
        MinimalPolynomialResult res;
        res.success = false;
        return res;
    }

    const std::vector<ModularCoeff> &element = t_v[table_idx - 1];
    if (verbose_uv) {
        std::cout << "[DEBUG] first_variable: element size=" << element.size() << ", values=[";
        for (size_t i = 0; i < element.size() && i < 20; ++i) {  // Limit output
            if (i > 0) std::cout << ",";
            std::cout << element[i];
        }
        if (element.size() > 20) std::cout << ",...";
        std::cout << "]" << std::endl;
        std::cout << "[DEBUG] Quotient basis size=" << quotient_basis.size() << std::endl;
    }

    // This centralized function handles all cases, including 1-D.
    return compute_minimal_polynomial_flint(element, i_xw, t_v, quotient_basis, prime);
}

BivariateResult
biv_lex(int32_t var_index,
        const MinimalPolynomialResult &minimal_poly_result,
        const std::vector<std::vector<int32_t>> &i_xw,
        const std::vector<std::vector<ModularCoeff>> &t_v,
        size_t quotient_basis_size,
        ModularCoeff prime) {
    BivariateResult result;
    result.success = false;

    const bool verbose_biv = false;
    if (verbose_biv) {
        std::cout << "biv_lex wrapper: var_index = " << var_index
                  << ", minimal_poly degree = " << minimal_poly_result.degree << std::endl;
    }

    // Initialize Gaussian reduction context
    GaussianReductionContext context;

    // The free set starts with the powers of T from the minimal polynomial computation
    // We have T^0, T^1, ..., T^(d-1) where d is the degree of the minimal polynomial
    std::vector<std::vector<ModularCoeff>> free_set = minimal_poly_result.powers;

    if (verbose_biv) {
        std::cout << "Free set size: " << free_set.size() << std::endl;
        for (size_t i = 0; i < free_set.size(); ++i) {
            std::cout << "  T^" << i << " = [";
            for (size_t j = 0; j < free_set[i].size(); ++j) {
                if (j > 0) std::cout << ",";
                std::cout << free_set[i][j];
            }
            std::cout << "]" << std::endl;
        }
    }

    // Special handling for linear case (quotient basis = {1})
    if (quotient_basis_size == 1 && minimal_poly_result.degree == 1) {
        // For a linear polynomial like x - c = 0, we have x = c in the quotient ring
        // The minimal polynomial is T - c = 0
        // The parameterization should be x = c (constant)

        // Extract the constant from the multiplication table
        // For variable var_index, check what it reduces to
        if (var_index > 0 && var_index <= static_cast<int32_t>(i_xw.size())) {
            int32_t var_idx = var_index - 1; // Convert to 0-based
            if (!i_xw[var_idx].empty()) {
                int32_t table_idx = i_xw[var_idx][0] - 1; // i_xw uses 1-based indices
                if (table_idx >= 0 && table_idx < static_cast<int32_t>(t_v.size()) && !t_v[table_idx].empty()) {
                    // The variable reduces to a constant
                    ModularCoeff stored_value = t_v[table_idx][0];

                    // For the linear case x - c = 0, we have x = c
                    // The constant c should be positive
                    // But due to modular representation, it might be stored as prime - c
                    ModularCoeff constant_value = stored_value;
                    if (stored_value > prime / 2) {
                        // Convert from negative representation
                        constant_value = prime - stored_value;
                    }

                    // Create a constant generator polynomial
                    result.generators.push_back({ constant_value });
                    result.basis.push_back({ { 0 }, { 0 } }); // Dummy basis for compatibility
                    result.success = true;

                    return result;
                }
            }
        }
    }

    // Call the actual bivariate lexicographic algorithm
    try {
        BivLexResult biv_result = ::julia_rur::biv_lex(t_v,
                                                       i_xw,
                                                       context,
                                                       free_set,
                                                       var_index, // Variable index (1-based)
                                                       prime);

        // Convert BivLexResult to BivariateResult
        result.success = !biv_result.generators.empty();

        if (result.success) {
            // Convert monomial basis to PP pairs
            result.basis.clear();
            for (const auto &mon : biv_result.monomial_basis) {
                PP deg_T(quotient_basis_size, 0);
                PP deg_xi(quotient_basis_size, 0);

                // Simple encoding: deg_T[0] = mon.deg_T, deg_xi[0] = mon.deg_xi
                if (quotient_basis_size > 0) {
                    deg_T[0] = mon.deg_T;
                    deg_xi[0] = mon.deg_xi;
                }

                result.basis.push_back({ deg_T, deg_xi });
            }

            // Extract h(t) and f(t) from the biv_lex result
            // The relation is h(t)*x - f(t) = 0
            // biv_lex returns the coefficients of this bivariate polynomial
            if (!biv_result.generators.empty() && !biv_result.leading_monomials.empty()) {
                // Get the coefficients from the first generator
                std::vector<ModularCoeff> coeffs = biv_result.generators[0];

                // Determine h(t) and f(t) based on the leading monomial
                auto leading_mon = biv_result.leading_monomials[0];

                // For the relation h(t)*x - f(t) = 0:
                // - Terms with x have coefficient from h(t)
                // - Terms without x have coefficient from -f(t)

                // Build h(t) and f(t) polynomials
                std::vector<ModularCoeff> h_t, f_t;

                // Simple case: if leading monomial is t^a * x^1, then we have relation
                if (leading_mon.deg_xi == 1) {
                    // The relation is of form h(t)*x - f(t) = 0
                    // The generator coefficients directly give us f(t) when h(t) = t^(deg_T)

                    if (coeffs.size() >= 1) {
                        // The generator gives us the representation of what x equals
                        // For leading monomial T^a * x^1, we have:
                        // - If deg_T = 0: relation is x = f(t), so h(t) = 1, f(t) = coeffs
                        // - If deg_T > 0: relation is t^a * x = f(t), so h(t) = t^a, f(t) = coeffs

                        if (leading_mon.deg_T == 0) {
                            // Simple case: x = f(t), so h(t) = 1
                            h_t = { 1 };
                            f_t = coeffs; // Direct assignment, no negation
                        } else {
                            // General case: t^a * x = f(t)
                            h_t.resize(leading_mon.deg_T + 1, 0);
                            h_t[leading_mon.deg_T] = 1; // h(t) = t^(deg_T)
                            f_t = coeffs;               // Direct assignment
                        }

                        // Debug output
                        if (verbose_biv) {
                            std::cout << "DEBUG: f(t) = ";
                            for (size_t i = 0; i < f_t.size(); ++i) {
                                if (i > 0) std::cout << ", ";
                                std::cout << f_t[i];
                            }
                            std::cout << std::endl;
                            std::cout << "DEBUG: h(t) = ";
                            for (size_t i = 0; i < h_t.size(); ++i) {
                                if (i > 0) std::cout << ", ";
                                std::cout << h_t[i];
                            }
                            std::cout << std::endl;
                            std::cout << "DEBUG: T(t) = ";
                            for (size_t i = 0; i < minimal_poly_result.coefficients.size(); ++i) {
                                if (i > 0) std::cout << ", ";
                                std::cout << minimal_poly_result.coefficients[i];
                            }
                            std::cout << std::endl;
                        }

                        // Check for special case: if h(t) = 1 and f(t) = t, then x = t
                        // This means the numerator should be t * T'(t) / T'(t) = t
                        bool is_direct_parameterization = false;
                        if (h_t.size() == 1 && h_t[0] == 1 &&                // h(t) = 1
                            f_t.size() == 2 && f_t[0] == 0 && f_t[1] == 1) { // f(t) = t
                            is_direct_parameterization = true;
                        }

                        std::vector<ModularCoeff> numerator;
                        if (is_direct_parameterization) {
                            // Direct parameterization: x = t
                            // The numerator is t * T'(t)
                            auto t_prime = polynomial_derivative(minimal_poly_result.coefficients, prime);
                            numerator = polynomial_multiply({ 0, 1 }, t_prime, prime); // t * T'(t)
                            numerator = polynomial_mod(numerator, minimal_poly_result.coefficients, prime);
                            std::cout << "Direct parameterization detected: x = t" << std::endl;
                        } else {
                            // General case: compute the RUR numerator using the formula
                            numerator = compute_rur_numerator(f_t, h_t, minimal_poly_result.coefficients, prime);
                        }

                        // Clear existing generators and add the computed numerator
                        result.generators.clear();
                        result.generators.push_back(numerator);

                        if (verbose_biv) {
                            std::cout << "Computed RUR numerator from h(t)=t^" << leading_mon.deg_T << ", f(t) has "
                                      << f_t.size() << " coeffs" << std::endl;
                            std::cout << "DEBUG: Numerator = ";
                            for (size_t i = 0; i < numerator.size(); ++i) {
                                if (i > 0) std::cout << ", ";
                                std::cout << numerator[i];
                            }
                            std::cout << std::endl;
                        }
                    }
                } else {
                    // Non-linear case: deg_xi > 1 means the separating element is not valid
                    std::cerr << "ERROR: biv_lex returned non-linear generator with deg_xi = " << leading_mon.deg_xi
                              << std::endl;
                    std::cerr << "This indicates the chosen element is not separating!" << std::endl;
                    result.success = false;
                    return result;
                }
            } else {
                // No generators or leading monomials
                std::cerr << "ERROR: biv_lex returned empty generators or leading monomials" << std::endl;
                result.success = false;
                return result;
            }

            if (verbose_biv) {
                std::cout << "biv_lex successful: " << biv_result.generators.size() << " generators, "
                          << biv_result.monomial_basis.size() << " basis monomials" << std::endl;
            }
        } else {
            std::cerr << "biv_lex failed: no generators found" << std::endl;
        }
    } catch (const std::exception &e) {
        std::cerr << "biv_lex exception: " << e.what() << std::endl;
        result.success = false;
    }

    return result;
}


bool
is_separating_element(const std::vector<ModularCoeff> &linear_form,
                      const std::vector<PP> &quotient_basis,
                      const std::vector<std::vector<int32_t>> &i_xw,
                      const std::vector<std::vector<ModularCoeff>> &t_v,
                      ModularCoeff prime) {
    // Compute minimal polynomial of the linear form
    MinimalPolynomialResult min_poly = compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);

    // Element is separating if minimal polynomial degree equals quotient ring dimension
    return min_poly.success && min_poly.degree == quotient_basis.size();
}

/*std::pair<std::vector<ModularCoeff>, int32_t>
find_separating_element(const std::vector<PP> &quotient_basis,
                        const std::vector<std::vector<int32_t>> &i_xw,
                        const std::vector<std::vector<ModularCoeff>> &t_v,
                        int num_variables,
                        ModularCoeff prime,
                        SeparatingStrategy strategy) {
    std::vector<ModularCoeff> coefficients;
    int32_t var_index = -1;

    // Strategy: CURRENT (default)
    // First try the last variable
    std::vector<ModularCoeff> last_var = element_to_vector(num_variables, i_xw, t_v, quotient_basis.size());

    if (is_separating_element(last_var, quotient_basis, i_xw, t_v, prime)) { return { last_var, num_variables }; }

    // Try other single variables
    for (int i = num_variables - 1; i >= 1; --i) {
        std::vector<ModularCoeff> var_vec = element_to_vector(i, i_xw, t_v, quotient_basis.size());

        if (is_separating_element(var_vec, quotient_basis, i_xw, t_v, prime)) { return { var_vec, i }; }
    }

    // Try deterministic linear forms if requested
    if (strategy == SeparatingStrategy::DETERMINISTIC) {
        std::cout << "\nTrying deterministic linear combinations...\n";
        // Iterate up to a reasonable bound based on quotient basis size
        int deterministic_trials = 2 * num_variables * static_cast<int>(quotient_basis.size());
        for (int i = 1; i <= deterministic_trials; ++i) {
            std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);
            ModularCoeff current_power = 1;

            for (int j = 1; j <= num_variables; ++j) {
                std::vector<ModularCoeff> var_vec = element_to_vector(j, i_xw, t_v, quotient_basis.size());

                for (size_t k = 0; k < linear_form.size(); ++k) {
                    linear_form[k] =
                      (linear_form[k] + (static_cast<AccModularCoeff>(current_power) * var_vec[k]) % prime) % prime;
                }
                current_power = (static_cast<AccModularCoeff>(current_power) * i) % prime;
            }

            if (is_separating_element(linear_form, quotient_basis, i_xw, t_v, prime)) {
                std::cout << "Found separating element with deterministic search (i=" << i << ")" << std::endl;
                return { linear_form, -1 };
            }
        }
    }

    // Try random linear forms
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(-10, 10); // Small coefficients

    const int max_trials = 20;
    for (int trial = 0; trial < max_trials; ++trial) {
        // Generate random linear form: λ₁x₁ + λ₂x₂ + ... + λₙxₙ
        std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);

        // Compute the linear combination in the quotient ring
        for (int i = 1; i <= num_variables; ++i) {
            int lambda = dis(gen);
            if (lambda == 0) continue; // Skip zero coefficients

            // Get representation of variable i in quotient basis
            std::vector<ModularCoeff> var_vec = element_to_vector(i, i_xw, t_v, quotient_basis.size());

            // Add λᵢ * xᵢ to the linear form
            ModularCoeff lambda_mod = (lambda < 0) ? (prime - ((-lambda) % prime)) % prime : lambda % prime;

            for (size_t j = 0; j < linear_form.size(); ++j) {
                linear_form[j] =
                  (linear_form[j] + (static_cast<AccModularCoeff>(lambda_mod) * var_vec[j]) % prime) % prime;
            }
        }

        // Check if this linear form is separating
        if (is_separating_element(linear_form, quotient_basis, i_xw, t_v, prime)) {
            std::cout << "Found separating linear form after " << (trial + 1) << " trials" << std::endl;
            return { linear_form, -1 }; // -1 indicates linear form, not single variable
        }
    }

    // Failed to find separating element
    std::cerr << "Failed to find separating element after " << max_trials << " trials" << std::endl;
    return { coefficients, -1 };
} */ //OLD CODE


std::pair<std::vector<int>, int32_t>
find_separating_element(const std::vector<PP> &quotient_basis,
                        const std::vector<std::vector<int32_t>> &i_xw,
                        const std::vector<std::vector<ModularCoeff>> &t_v,
                        int num_variables,
                        ModularCoeff prime,
                        SeparatingStrategy strategy) {
    std::vector<int> coefficients;
    int32_t var_index = -1;

    // Try single variables first
    std::cout << "DEBUG find_separating_element: num_variables=" << num_variables << ", i_xw.size()=" << i_xw.size()
              << std::endl;
    for (int i = num_variables; i >= 1; --i) {
        // Get the vector for variable 'i' from the multiplication table.
        // This is the correct way to get the representation of a variable,
        // as it may not be a basis element itself.
        std::cout << "  Checking variable " << i << ": i_xw[" << (i - 1)
                  << "].size()=" << (i <= static_cast<int32_t>(i_xw.size()) ? i_xw[i - 1].size() : 0) << std::endl;
        if (i > static_cast<int32_t>(i_xw.size()) || i_xw[i - 1].empty()) continue;
        int32_t table_idx = i_xw[i - 1][0]; // Index for var_i * 1
        if (table_idx <= 0 || table_idx > static_cast<int32_t>(t_v.size())) continue;

        const std::vector<ModularCoeff> &var_vec = t_v[table_idx - 1];

        if (is_separating_element(var_vec, quotient_basis, i_xw, t_v, prime)) {
            std::vector<int> sep_coeffs(num_variables, 0);
            sep_coeffs[i - 1] = 1;
            return { sep_coeffs, i };
        }
    }

    // Try deterministic linear forms if requested
    if (strategy == SeparatingStrategy::DETERMINISTIC) {
        std::cout << "\nTrying deterministic linear combinations...\n";
        int deterministic_trials = 2 * num_variables * static_cast<int>(quotient_basis.size());
        for (int i = 1; i <= deterministic_trials; ++i) {
            std::vector<int> sep_coeffs_int(num_variables, 0);
            ModularCoeff current_power = 1;

            for (int j = 1; j <= num_variables; ++j) {
                sep_coeffs_int[j - 1] = static_cast<int>(current_power);
                current_power = (static_cast<AccModularCoeff>(current_power) * i) % prime;
            }

            MinimalPolynomialResult test_mp =
              compute_minimal_polynomial_flint(sep_coeffs_int, i_xw, t_v, quotient_basis, prime);
            if (test_mp.success && test_mp.degree == quotient_basis.size()) {
                std::cout << "Found separating element with deterministic search (i=" << i << ")" << std::endl;
                return { sep_coeffs_int, -1 };
            }
        }
    }

    // Try random linear forms
    std::random_device rd;
    std::mt19937 gen(rd());
    const int max_trials = 20;
    for (int trial = 0; trial < max_trials; ++trial) {
        std::vector<int> linear_form_coeffs =
          compute_random_linear_form(quotient_basis, i_xw, t_v, num_variables, prime, gen);

        MinimalPolynomialResult test_mp =
          compute_minimal_polynomial_flint(linear_form_coeffs, i_xw, t_v, quotient_basis, prime);
        if (test_mp.success && test_mp.degree == quotient_basis.size()) {
            std::cout << "Found separating linear form after " << (trial + 1) << " trials" << std::endl;
            return { linear_form_coeffs, -1 };
        }
    }

    std::cerr << "Failed to find separating element after " << max_trials << " trials" << std::endl;
    return { coefficients, -1 };
}

std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>>
compute_univariate_parameterization(const std::vector<PP> &quotient_basis,
                                    const std::vector<std::vector<int32_t>> &i_xw,
                                    std::vector<std::vector<ModularCoeff>> &t_v,
                                    int num_variables,
                                    ModularCoeff prime) {
    MinimalPolynomialResult min_poly;
    std::vector<BivariateResult> parameterizations(num_variables);

    std::cout << "compute_univariate_parameterization: quotient_basis size = " << quotient_basis.size()
              << ", num_variables = " << num_variables << std::endl;
    std::cout << "i_xw size = " << i_xw.size() << ", t_v size = " << t_v.size() << std::endl;

    // Step 1: Find separating element
    auto [sep_coeffs, sep_var] = find_separating_element(quotient_basis, i_xw, t_v, num_variables, prime);

    if (sep_coeffs.empty()) {
        std::cerr << "Failed to find separating element" << std::endl;
        return { false, min_poly, parameterizations };
    }

    std::cout << "Found separating element, sep_var = " << sep_var << std::endl;

    // Step 2: Compute minimal polynomial
    if (sep_var > 0) {
        // Single variable is separating
        min_poly = first_variable(sep_var, i_xw, t_v, quotient_basis, prime);
    } else {
        // Linear form is separating
        min_poly = compute_minimal_polynomial_flint(sep_coeffs, i_xw, t_v, quotient_basis, prime);
    }

    if (!min_poly.success) {
        std::cerr << "Failed to compute minimal polynomial" << std::endl;
        return { false, min_poly, parameterizations };
    }

    std::cout << "Computed minimal polynomial of degree " << min_poly.degree << std::endl;

    // Step 3: Compute parameterizations for each variable
    for (int i = 0; i < num_variables; ++i) {
        if (std::getenv("RUR_VERBOSE_PROGRESS")) {
            std::cout << "Computing parameterization for variable " << i << std::endl;
        }
        parameterizations[i] = biv_lex(i + 1, min_poly, i_xw, t_v, quotient_basis.size(), prime);
        if (!parameterizations[i].success) {
            std::cerr << "Failed to compute parameterization for variable " << i << std::endl;
        }
    }

    // Check if all succeeded
    bool all_success = true;
    for (const auto &param : parameterizations) {
        if (!param.success) {
            all_success = false;
            break;
        }
    }

    return { all_success, min_poly, parameterizations };
}

// NEW IMPLEMENTATION USING FLINT
MinimalPolynomialResult
compute_minimal_polynomial_flint(const std::vector<int> &element_coeffs,
                                 const std::vector<std::vector<int32_t>> &i_xw,
                                 const std::vector<std::vector<ModularCoeff>> &t_v,
                                 const std::vector<PP> &quotient_basis,
                                 ModularCoeff prime) {
    // Build the element vector from integer coefficients of variables
    std::vector<ModularCoeff> element(quotient_basis.size(), 0);
    for (size_t i = 0; i < element_coeffs.size(); ++i) {
        ModularCoeff c_mod =
          (element_coeffs[i] < 0) ? (prime - ((-element_coeffs[i]) % prime)) % prime : element_coeffs[i] % prime;
        std::vector<ModularCoeff> var_vec =
          element_to_vector(static_cast<int32_t>(i + 1), i_xw, t_v, quotient_basis.size());
        for (size_t j = 0; j < element.size(); ++j) {
            element[j] = (element[j] + static_cast<AccModularCoeff>(c_mod) * var_vec[j]) % prime;
        }
    }
    // Delegate to overload that works with element vector directly
    return compute_minimal_polynomial_flint(element, i_xw, t_v, quotient_basis, prime);
}

// Overload: element is already represented in the quotient basis
MinimalPolynomialResult
compute_minimal_polynomial_flint(const std::vector<ModularCoeff> &element,
                                 const std::vector<std::vector<int32_t>> &i_xw,
                                 const std::vector<std::vector<ModularCoeff>> &t_v,
                                 const std::vector<PP> &quotient_basis,
                                 ModularCoeff prime) {
    MinimalPolynomialResult result;
    result.success = false;
    size_t d = quotient_basis.size();

    const bool verbose_minpoly = false;  // Disabled for clean output
    if (verbose_minpoly) {
        std::cout << "\n=== compute_minimal_polynomial_flint (NEW) ===" << std::endl;
        std::cout << "Quotient basis size d = " << d << std::endl;
        std::cout << "Element vector: [";
        for (size_t i = 0; i < element.size() && i < 20; ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << element[i];
        }
        if (element.size() > 20) std::cout << ",...";
        std::cout << "]" << std::endl;
    }

    // Special case for 1-dimensional quotient ring
    if (d == 1) {
        if (element.size() != 1) {
            std::cout << "ERROR: Element size mismatch for 1-dim case!" << std::endl;
            return result;
        }
        ModularCoeff c = element[0];
        result.coefficients = { (prime - c) % prime, 1 }; // T - c
        result.degree = 1;
        result.powers.push_back(element);
        result.success = true;
        return result;
    }

    // Precompute multiplication-by-element linear map M (d x d)
    // Column j is the product (basis_j * element)
    std::vector<std::vector<ModularCoeff>> M_rows(d, std::vector<ModularCoeff>(d, 0));
    {
        std::vector<ModularCoeff> basis_vec(d, 0);
        for (size_t j = 0; j < d; ++j) {
            std::fill(basis_vec.begin(), basis_vec.end(), 0);
            basis_vec[j] = 1;
            std::vector<ModularCoeff> col =
              multiply_elements_in_quotient(basis_vec, element, quotient_basis, i_xw, t_v, prime);
            for (size_t i = 0; i < d; ++i) { M_rows[i][j] = col[i]; }
        }
    }

    // Use M to build matrix of powers and find linear dependence
    flint_linalg::NModMat power_matrix(d + 1, d, prime);

    // Current power of the element
    std::vector<ModularCoeff> current_power(d, 0);
    current_power[0] = 1; // T^0 = 1

    // Store powers for later use
    result.powers.clear();

    auto matvec = [&](const std::vector<ModularCoeff> &x) {
        std::vector<ModularCoeff> y(d, 0);
        for (size_t i = 0; i < d; ++i) {
            AccModularCoeff acc = 0;
            const auto &row = M_rows[i];
            for (size_t j = 0; j < d; ++j) { acc += static_cast<AccModularCoeff>(row[j]) * x[j]; }
            y[i] = static_cast<ModularCoeff>(acc % prime);
        }
        return y;
    };

    for (size_t k = 0; k <= d; ++k) {
        if (verbose_minpoly) { std::cout << "  Power T^" << k << " = ["; }
        for (size_t j = 0; j < d; ++j) {
            if (verbose_minpoly) {
                if (j > 0) std::cout << ", ";
                std::cout << current_power[j];
            }
            power_matrix.set_entry(k, j, current_power[j]);
        }
        if (verbose_minpoly) { std::cout << "]" << std::endl; }

        // Store this power
        result.powers.push_back(current_power);

        // Compute next power if we're not at the last iteration
        if (k < d) { current_power = matvec(current_power); }
    }

    // One-shot nullspace on full transposed power matrix (d x (d+1))
    nmod_mat_t transposed;
    nmod_mat_init(transposed, d, d + 1, prime);
    for (size_t i = 0; i <= d; ++i) {
        for (size_t j = 0; j < d; ++j) { nmod_mat_set_entry(transposed, j, i, power_matrix.get_entry(i, j)); }
    }

    nmod_mat_t nullsp;
    nmod_mat_init(nullsp, d + 1, d, prime);
    slong nullity = nmod_mat_nullspace(nullsp, transposed);
    if (verbose_minpoly) { std::cout << "  Nullspace dimension: " << nullity << std::endl; }

    if (nullity <= 0) {
        if (verbose_minpoly) { std::cout << "ERROR: No linear dependence found up to degree " << d << std::endl; }
        nmod_mat_clear(nullsp);
        nmod_mat_clear(transposed);
        return result;
    }

    // Choose the null vector with smallest degree (highest non-zero index minimal)
    size_t best_deg = d + 1;
    std::vector<ModularCoeff> best_coeffs;
    for (slong col = 0; col < nullity; ++col) {
        std::vector<ModularCoeff> coeffs(d + 1);
        size_t deg = 0;
        for (size_t i = 0; i <= d; ++i) { coeffs[i] = nmod_mat_get_entry(nullsp, i, col); }
        for (size_t i = coeffs.size(); i-- > 0;) {
            if (coeffs[i] != 0) {
                deg = i;
                break;
            }
        }
        if (deg < best_deg) {
            best_deg = deg;
            best_coeffs.swap(coeffs);
        }
    }

    if (best_deg == d + 1) {
        nmod_mat_clear(nullsp);
        nmod_mat_clear(transposed);
        return result;
    }

    // Make monic
    if (best_deg > 0) {
        ModularCoeff leading = best_coeffs[best_deg];
        if (leading != 1 && leading != 0) {
            ModularCoeff inv = modular_inverse(leading, prime);
            for (auto &c : best_coeffs) { c = (static_cast<AccModularCoeff>(c) * inv) % prime; }
        }
    }

    // Apply square-free reduction to handle non-radical ideals
    size_t original_deg = best_deg;
    auto square_free_coeffs = compute_square_free_part(best_coeffs, prime, &original_deg);
    
    // Update result with square-free polynomial
    result.coefficients = std::move(square_free_coeffs);
    result.degree = result.coefficients.size() > 0 ? result.coefficients.size() - 1 : 0;
    result.original_degree = original_deg;  // Store original degree for multiplicity info
    
    if (original_deg > result.degree && result.degree > 0) {
        result.multiplicity = original_deg / result.degree;
        if (verbose_minpoly || getenv("RUR_MULTIPLICITY_VERBOSE")) {
            std::cout << "  [Square-free applied] Original degree: " << original_deg 
                     << " -> Square-free degree: " << result.degree 
                     << " (multiplicity ~" << result.multiplicity << ")" << std::endl;
        }
    }
    
    result.success = true;

    nmod_mat_clear(nullsp);
    nmod_mat_clear(transposed);
    return result;
}

/**
 * @brief Compute square-free part of polynomial using FLINT
 * 
 * Given polynomial f, computes f / gcd(f, f') to remove repeated factors.
 * This is essential for handling non-radical ideals where the minimal polynomial
 * may have high multiplicities.
 */
std::vector<ModularCoeff>
compute_square_free_part(const std::vector<ModularCoeff> &poly, 
                         ModularCoeff prime,
                         size_t *original_degree) {
    
    // Handle trivial cases
    if (poly.empty()) return poly;
    
    // Find actual degree (highest non-zero coefficient)
    size_t deg = 0;
    for (size_t i = poly.size(); i-- > 0;) {
        if (poly[i] != 0) {
            deg = i;
            break;
        }
    }
    
    if (original_degree) *original_degree = deg;
    
    // Degree 0 or 1 polynomials are already square-free
    if (deg <= 1) return poly;
    
    // Initialize FLINT polynomial
    nmod_poly_t f, f_prime, g, result_poly;
    nmod_poly_init(f, prime);
    nmod_poly_init(f_prime, prime);
    nmod_poly_init(g, prime);
    nmod_poly_init(result_poly, prime);
    
    // Set coefficients of f
    for (size_t i = 0; i <= deg; ++i) {
        nmod_poly_set_coeff_ui(f, i, poly[i]);
    }
    
    // Compute derivative f'
    nmod_poly_derivative(f_prime, f);
    
    // Compute g = gcd(f, f')
    nmod_poly_gcd(g, f, f_prime);
    
    // Compute square-free part: f / g
    if (nmod_poly_degree(g) > 0) {
        // Non-trivial GCD means there are repeated factors
        nmod_poly_div(result_poly, f, g);
        
        if (getenv("RUR_SQUARE_FREE_VERBOSE")) {
            std::cout << "[Square-free reduction]" << std::endl;
            std::cout << "  Original degree: " << deg << std::endl;
            std::cout << "  GCD degree: " << nmod_poly_degree(g) << std::endl;
            std::cout << "  Square-free degree: " << nmod_poly_degree(result_poly) << std::endl;
            
            // Estimate multiplicity
            if (nmod_poly_degree(result_poly) > 0) {
                size_t multiplicity = deg / nmod_poly_degree(result_poly);
                std::cout << "  Estimated multiplicity: " << multiplicity << std::endl;
            }
        }
    } else {
        // Already square-free
        nmod_poly_set(result_poly, f);
    }
    
    // Extract coefficients
    size_t result_deg = nmod_poly_degree(result_poly);
    std::vector<ModularCoeff> result(result_deg + 1);
    for (size_t i = 0; i <= result_deg; ++i) {
        result[i] = nmod_poly_get_coeff_ui(result_poly, i);
    }
    
    // Clean up
    nmod_poly_clear(f);
    nmod_poly_clear(f_prime);
    nmod_poly_clear(g);
    nmod_poly_clear(result_poly);
    
    return result;
}

} // namespace julia_rur