#ifndef JULIA_RUR_UNIVARIATE_PARAMETERIZATION_HPP
#define JULIA_RUR_UNIVARIATE_PARAMETERIZATION_HPP

#include "data_structures.hpp"
// #include "multiplication_tables.hpp" // not needed in header
#include <algorithm>
#include <random>
#include <tuple>
#include <vector>


namespace julia_rur {

// Helper function to compute element representation in quotient basis
std::vector<ModularCoeff>
element_to_vector(int32_t var_index,
                  const std::vector<std::vector<int32_t>> &i_xw,
                  const std::vector<std::vector<ModularCoeff>> &t_v,
                  size_t quotient_basis_size);

/**
 * @brief Result of the minimal polynomial computation for the separating element
 *
 * Contains the minimal polynomial f(T) such that f(sep_element) = 0 in the quotient ring
 */
struct MinimalPolynomialResult {
    std::vector<ModularCoeff> coefficients;        // Coefficients of minimal polynomial (low to high degree)
    std::vector<std::vector<ModularCoeff>> powers; // Powers of T in quotient basis (T^0, T^1, ..., T^(d-1))
    size_t degree;                                 // Degree of minimal polynomial
    bool success;                                  // Whether computation succeeded
};

/**
 * @brief Result of bivariate lexicographic computation
 *
 * Contains the bivariate ideal generators that allow expressing xi in terms of T
 */
struct BivariateResult {
    std::vector<std::pair<PP, PP>> basis;              // Bivariate basis: pairs (deg_in_T, deg_in_xi)
    std::vector<std::vector<ModularCoeff>> generators; // Generator polynomials
    bool success;                                      // Whether computation succeeded
};

/**
 * @brief Separating element strategy enumeration
 *
 * Different strategies for finding or constructing a separating linear form
 */
// enum class SeparatingStrategy {
//     CURRENT,        // Default: try last variable, then systematic search
//     RANDOM,         // Random coefficients in bounded range
//     L0_NORM,        // Prefer sparse linear forms (powerset approach)
//     DETERMINISTIC,  // Systematic search using power sequences
//     MRON_0L         // Reverse order l0 norm strategy
// };

/**
 * @brief Compute minimal polynomial of an element in the quotient ring
 *
 * Uses multiplication tables to find the minimal polynomial of the given element
 * by iteratively computing powers until linear dependence is detected.
 *
 * @param element Initial element vector in quotient basis
 * @param i_xw Variable multiplication indices
 * @param t_v Coefficient vectors from multiplication tables
 * @param quotient_basis_size Size of quotient basis
 * @param prime Modular arithmetic prime
 * @return MinimalPolynomialResult containing polynomial and powers
 */
// DISABLED: OLD IMPLEMENTATION - PRODUCES WRONG COEFFICIENTS
// Use compute_minimal_polynomial_flint instead
MinimalPolynomialResult
compute_minimal_polynomial_OLD_DISABLED(const std::vector<ModularCoeff> &element,
                                        const std::vector<std::vector<int32_t>> &i_xw,
                                        const std::vector<std::vector<ModularCoeff>> &t_v,
                                        const std::vector<PP> &quotient_basis,
                                        ModularCoeff prime);

/**
 * @brief NEW: Compute minimal polynomial using FLINT linear algebra
 *
 * This is a new implementation using FLINT's robust finite field linear algebra.
 * Kept separate from the original to allow testing and comparison.
 *
 * @param element Initial element vector in quotient basis
 * @param i_xw Variable multiplication indices
 * @param t_v Coefficient vectors from multiplication tables
 * @param quotient_basis_size Size of quotient basis
 * @param prime Modular arithmetic prime
 * @return MinimalPolynomialResult containing polynomial and powers
 */
MinimalPolynomialResult
compute_minimal_polynomial_flint(const std::vector<int> &element,
                                 const std::vector<std::vector<int32_t>> &i_xw,
                                 const std::vector<std::vector<ModularCoeff>> &t_v,
                                 const std::vector<PP> &quotient_basis,
                                 ModularCoeff prime);

// Overload: compute minimal polynomial when the element is already given in the
// quotient basis as a ModularCoeff vector (e.g., a single variable's image)
MinimalPolynomialResult
compute_minimal_polynomial_flint(const std::vector<ModularCoeff> &element,
                                 const std::vector<std::vector<int32_t>> &i_xw,
                                 const std::vector<std::vector<ModularCoeff>> &t_v,
                                 const std::vector<PP> &quotient_basis,
                                 ModularCoeff prime);

// API COMPATIBILITY: Redirect old function name to FLINT implementation
inline MinimalPolynomialResult
compute_minimal_polynomial(const std::vector<int> &element,
                           const std::vector<std::vector<int32_t>> &i_xw,
                           const std::vector<std::vector<ModularCoeff>> &t_v,
                           const std::vector<PP> &quotient_basis,
                           ModularCoeff prime) {
    return compute_minimal_polynomial_flint(element, i_xw, t_v, quotient_basis, prime);
}

// Overload: when the element is already represented in the quotient basis
inline MinimalPolynomialResult
compute_minimal_polynomial(const std::vector<ModularCoeff> &element,
                           const std::vector<std::vector<int32_t>> &i_xw,
                           const std::vector<std::vector<ModularCoeff>> &t_v,
                           const std::vector<PP> &quotient_basis,
                           ModularCoeff prime) {
    return compute_minimal_polynomial_flint(element, i_xw, t_v, quotient_basis, prime);
}

/**
 * @brief First variable algorithm - computes minimal polynomial of separating element
 *
 * This is the Julia first_variable function that computes the minimal polynomial
 * of the separating element (typically the last variable or a linear form).
 *
 * @param separating_var_index Index of separating variable (1-based)
 * @param i_xw Variable multiplication indices
 * @param t_v Coefficient vectors
 * @param quotient_basis_size Size of quotient basis
 * @param prime Modular arithmetic prime
 * @return MinimalPolynomialResult
 */
MinimalPolynomialResult
first_variable(int32_t separating_var_index,
               const std::vector<std::vector<int32_t>> &i_xw,
               std::vector<std::vector<ModularCoeff>> &t_v,
               const std::vector<PP> &quotient_basis,
               ModularCoeff prime);

/**
 * @brief Bivariate lexicographic algorithm
 *
 * Extends the univariate ideal to bivariate by introducing another variable.
 * Computes generators that express xi in terms of the separating element T.
 *
 * @param var_index Variable index to parameterize (1-based)
 * @param minimal_poly_result Result from first_variable
 * @param i_xw Variable multiplication indices
 * @param t_v Coefficient vectors
 * @param quotient_basis_size Size of quotient basis
 * @param prime Modular arithmetic prime
 * @return BivariateResult containing generators
 */
BivariateResult
biv_lex(int32_t var_index,
        const MinimalPolynomialResult &minimal_poly_result,
        const std::vector<std::vector<int32_t>> &i_xw,
        const std::vector<std::vector<ModularCoeff>> &t_v,
        size_t quotient_basis_size,
        ModularCoeff prime);

/**
 * @brief Check if a linear form is separating
 *
 * A linear form is separating if its minimal polynomial has degree equal
 * to the dimension of the quotient ring.
 *
 * @param linear_form Coefficients of the linear form
 * @param quotient_basis The quotient basis
 * @param i_xw Variable multiplication indices
 * @param t_v Coefficient vectors
 * @param prime Modular arithmetic prime
 * @return true if separating, false otherwise
 */
bool
is_separating_element(const std::vector<ModularCoeff> &linear_form,
                      const std::vector<PP> &quotient_basis,
                      const std::vector<std::vector<int32_t>> &i_xw,
                      const std::vector<std::vector<ModularCoeff>> &t_v,
                      ModularCoeff prime);

/**
 * @brief Find or construct a separating element
 *
 * Tries various strategies to find a linear form that separates all points.
 *
 * @param quotient_basis The quotient basis
 * @param i_xw Variable multiplication indices
 * @param t_v Coefficient vectors
 * @param num_variables Number of variables
 * @param prime Modular arithmetic prime
 * @param strategy Strategy to use for finding separating element
 * @return Pair of (coefficients, variable_index) where variable_index is -1 for linear form
 */
std::pair<std::vector<int>, int32_t>
find_separating_element(const std::vector<PP> &quotient_basis,
                        const std::vector<std::vector<int32_t>> &i_xw,
                        const std::vector<std::vector<ModularCoeff>> &t_v,
                        int num_variables,
                        ModularCoeff prime,
                        SeparatingStrategy strategy = SeparatingStrategy::CURRENT);

/**
 * @brief Complete univariate parameterization
 *
 * Main function that orchestrates the entire parameterization:
 * 1. Finds separating element
 * 2. Computes minimal polynomial
 * 3. Computes parameterizations for each variable
 *
 * @param quotient_basis The quotient basis
 * @param i_xw Variable multiplication indices
 * @param t_v Coefficient vectors
 * @param num_variables Number of variables
 * @param prime Modular arithmetic prime
 * @return Tuple of (success, minimal_poly, parameterizations) where parameterizations[i] = fi(T)/f'(T)
 */
std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>>
compute_univariate_parameterization(const std::vector<PP> &quotient_basis,
                                    const std::vector<std::vector<int32_t>> &i_xw,
                                    std::vector<std::vector<ModularCoeff>> &t_v,
                                    int num_variables,
                                    ModularCoeff prime);


// Forward declaration - function moved to .cpp file
std::vector<int>
compute_random_linear_form(const std::vector<PP> &quotient_basis,
                           const std::vector<std::vector<int32_t>> &i_xw,
                           const std::vector<std::vector<ModularCoeff>> &t_v,
                           int num_variables,
                           ModularCoeff prime,
                           std::mt19937 &gen);


// (removed) deterministic builder for linear form


inline std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>>
try_separating_element(const std::vector<PP> &quotient_basis,
                       const std::vector<std::vector<int32_t>> &i_xw,
                       std::vector<std::vector<ModularCoeff>> &t_v,
                       int num_variables,
                       ModularCoeff prime,
                       const std::vector<int> &sep_coeffs,
                       int32_t sep_var) {
    // CRITICAL CHECK: Ensure quotient_basis[0] is the constant monomial '1'
    if (!quotient_basis.empty()) {
        bool is_constant =
          std::all_of(quotient_basis[0].begin(), quotient_basis[0].end(), [](int e) { return e == 0; });
        if (!is_constant) {
            std::cout << "ERROR in try_separating_element: quotient_basis[0] is not the constant '1'!" << std::endl;
            std::cout << "  quotient_basis[0] = [";
            for (size_t i = 0; i < quotient_basis[0].size(); ++i) {
                if (i > 0) std::cout << ",";
                std::cout << quotient_basis[0][i];
            }
            std::cout << "]" << std::endl;
            std::cout << "  This will cause incorrect multiplication table lookups!" << std::endl;

            // We can't fix it here since quotient_basis is const
            // But at least we can warn about the problem
        }
    }

    MinimalPolynomialResult min_poly;
    std::vector<BivariateResult> parameterizations(num_variables);

    // Compute minimal polynomial
    if (sep_var > 0) {
        // Single variable is separating
        min_poly = first_variable(sep_var, i_xw, t_v, quotient_basis, prime);
    } else {
        // Linear form is separating (specified by integer coefficients)
        min_poly = compute_minimal_polynomial_flint(sep_coeffs, i_xw, t_v, quotient_basis, prime);
    }
    if (!min_poly.success || min_poly.degree != quotient_basis.size()) {
        if (min_poly.degree != quotient_basis.size()) {
            std::cout << "Random linear form is not separating (minimal poly degree " << min_poly.degree
                      << " != quotient size " << quotient_basis.size() << ")" << std::endl;
        } else {
            std::cerr << "Failed to compute minimal polynomial" << std::endl;
        }
        return { false, min_poly, parameterizations };
    }
    std::cout << "Computed minimal polynomial of degree " << min_poly.degree << std::endl;

    // Compute parameterizations for each variable
    bool all_success = true;
    for (int i = 0; i < num_variables; ++i) {
        std::cout << "Computing parameterization for variable " << i << std::endl;
        parameterizations[i] = biv_lex(i + 1, min_poly, i_xw, t_v, quotient_basis.size(), prime);
        if (!parameterizations[i].success) {
            std::cerr << "Failed to compute parameterization for variable " << i << std::endl;
            all_success = false;
            break; // No point continuing if one failed
        }
    }

    return { all_success, min_poly, parameterizations };
}


} // namespace julia_rur

#endif // JULIA_RUR_UNIVARIATE_PARAMETERIZATION_HPP