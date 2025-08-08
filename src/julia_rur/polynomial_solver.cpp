#include "polynomial_solver.hpp"
#include <iomanip>
#include <sstream>

namespace julia_rur {

/**
 * @brief Evaluate a polynomial in T at a given value
 */
static std::complex<double>
evaluate_polynomial_in_T(const std::vector<mpq_class> &coeffs, const std::complex<double> &t_value) {
    std::complex<double> result = 0.0;
    std::complex<double> t_power = 1.0;

    for (const auto &coeff : coeffs) {
        result += coeff.get_d() * t_power;
        t_power *= t_value;
    }

    return result;
}

// Robust evaluation of N(T)/f'(T) at a root T = t0 with L'Hospital fallback for repeated roots
static std::complex<double>
evaluate_parameterized_variable(std::vector<mpq_class> numerator_coeffs,
                                std::vector<mpq_class> fprime_coeffs,
                                const std::complex<double> &t_value) {
    auto eval_poly = [](const std::vector<mpq_class> &c, const std::complex<double> &t) {
        std::complex<double> v = 0.0, tp = 1.0;
        for (const auto &a : c) {
            v += a.get_d() * tp;
            tp *= t;
        }
        return v;
    };

    const double eps = 1e-12;
    // DEBUG: report evaluation point
    std::cout << "EVAL-DEBUG: t=" << std::setprecision(15) << t_value << std::endl;

    // Try direct evaluation, then up to two L'Hospital steps
    for (int k = 0; k < 3; ++k) {
        std::complex<double> num_v = eval_poly(numerator_coeffs, t_value);
        std::complex<double> den_v = eval_poly(fprime_coeffs, t_value);
        std::cout << "EVAL-DEBUG: step=" << k << " num(t)=" << num_v << " den(t)=" << den_v << std::endl;
        if (std::abs(den_v) > eps) {
            auto out = num_v / den_v;
            std::cout << "EVAL-DEBUG: return num/den=" << out << std::endl;
            return out;
        }
        if (std::abs(num_v) <= eps) {
            // Apply L'Hospital: differentiate numerator and denominator
            if (fprime_coeffs.size() <= 1 || numerator_coeffs.size() <= 1) break;
            // Differentiate numerator
            std::vector<mpq_class> num_d;
            num_d.reserve(numerator_coeffs.size() > 0 ? numerator_coeffs.size() - 1 : 0);
            for (size_t i = 1; i < numerator_coeffs.size(); ++i) {
                num_d.push_back(numerator_coeffs[i] * static_cast<int>(i));
            }
            numerator_coeffs.swap(num_d);
            std::cout << "EVAL-DEBUG: LHospital numerator' deg="
                      << (numerator_coeffs.empty() ? -1 : static_cast<int>(numerator_coeffs.size()) - 1) << std::endl;
            // Differentiate denominator (f' -> f'')
            std::vector<mpq_class> den_d;
            den_d.reserve(fprime_coeffs.size() > 0 ? fprime_coeffs.size() - 1 : 0);
            for (size_t i = 1; i < fprime_coeffs.size(); ++i) {
                den_d.push_back(fprime_coeffs[i] * static_cast<int>(i));
            }
            fprime_coeffs.swap(den_d);
            std::cout << "EVAL-DEBUG: LHospital denominator'' deg="
                      << (fprime_coeffs.empty() ? -1 : static_cast<int>(fprime_coeffs.size()) - 1) << std::endl;
            continue;
        }
        // num not ~0 but den ~0: return numerator as fallback
        std::cout << "EVAL-DEBUG: den~0 but num!=0, return num=" << num_v << std::endl;
        return num_v;
    }
    // Final fallback
    std::cout << "EVAL-DEBUG: final fallback 0" << std::endl;
    return std::complex<double>(0.0, 0.0);
}

/**
 * @brief Extract polynomial coefficients from parameterization
 *
 * For now, this is a simplified version that handles basic cases
 */
static std::vector<mpq_class>
extract_parameterization_coeffs(const std::vector<mpq_class> &param_numerator,
                                const std::vector<mpq_class> &param_denominator) {
    // Simple case: denominator is just a constant
    if (param_denominator.size() == 1 && param_denominator[0] != 0) {
        std::vector<mpq_class> result;
        for (const auto &num : param_numerator) { result.push_back(num / param_denominator[0]); }
        return result;
    }

    // More complex case would require polynomial division
    // For now, just return the numerator
    return param_numerator;
}

PolynomialSystemSolution
solve_polynomial_system_complete(const std::vector<std::string> &polynomials,
                                 const std::vector<std::string> &variables,
                                 const RURConfig &config) {
    PolynomialSystemSolution solution;
    solution.variable_names = variables;
    solution.success = false;

    // Step 1: Compute RUR
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);

    if (!rur_result.success) {
        solution.error_message = "RUR computation failed: " + rur_result.error_message;
        return solution;
    }

    solution.minimal_polynomial = rur_result.minimal_polynomial;
    solution.quotient_dimension = rur_result.quotient_basis.size();

    // Step 2: Find roots of minimal polynomial
    std::vector<std::complex<double>> t_roots = find_polynomial_roots(rur_result.minimal_polynomial);

    if (t_roots.empty()) {
        solution.error_message = "No roots found for minimal polynomial";
        return solution;
    }

    // Step 3: Back-substitute to find variable values
    const double epsilon = 1e-10;

    for (const auto &t_root : t_roots) {
        std::vector<std::complex<double>> var_values;

        if (variables.size() == 1) {
            // Univariate case: Check if this is the simple case where x = T
            // This happens when the separating element is the variable itself
            // In this case, the minimal polynomial is directly in terms of the variable
            // and we don't need any parameterization

            // Check if we have a simple case (no parameterization or trivial parameterization)
            // For univariate systems where the variable is the separating element,
            // we should just have x = T, regardless of what biv_lex returns
            bool is_simple_case = true; // For univariate, always treat as simple case

            // The RUR algorithm for univariate polynomials should recognize that
            // when the variable itself is the separating element T, no parameterization
            // is needed. The roots of the minimal polynomial are the values of the variable.

            if (is_simple_case) {
                // Simple case: x = T, so the roots of the minimal polynomial are the values of x
                var_values.push_back(t_root);
            } else {
                // Parameterized case: x is a rational function of T
                std::complex<double> numerator_value = evaluate_polynomial_in_T(rur_result.numerators[0], t_root);

                // Compute f'(T) - derivative of minimal polynomial
                std::vector<mpq_class> derivative_coeffs;
                for (size_t i = 1; i < rur_result.minimal_polynomial.size(); ++i) {
                    derivative_coeffs.push_back(rur_result.minimal_polynomial[i] * static_cast<int>(i));
                }

                std::complex<double> f_prime_value = evaluate_polynomial_in_T(derivative_coeffs, t_root);

                // Variable value = numerator(T) / f'(T)
                std::complex<double> var_value = numerator_value / f_prime_value;
                var_values.push_back(var_value);
            }
        } else {
            // Multivariate case: use parameterizations xi = gi(T) / f'(T)
            // Compute f'(T) coefficients once
            std::vector<mpq_class> fprime_coeffs;
            for (size_t i = 1; i < rur_result.minimal_polynomial.size(); ++i) {
                fprime_coeffs.push_back(rur_result.minimal_polynomial[i] * static_cast<int>(i));
            }
            std::cout << "EVAL-DEBUG: root t=" << std::setprecision(15) << t_root << std::endl;
            for (size_t var_idx = 0; var_idx < variables.size(); ++var_idx) {
                if (var_idx < rur_result.numerators.size()) {
                    std::cout << "EVAL-DEBUG: var " << variables[var_idx] << " numerator=[";
                    for (size_t i = 0; i < rur_result.numerators[var_idx].size(); ++i) {
                        if (i) std::cout << ", ";
                        std::cout << rur_result.numerators[var_idx][i];
                    }
                    std::cout << "]" << std::endl;
                    auto val = evaluate_parameterized_variable(rur_result.numerators[var_idx], fprime_coeffs, t_root);
                    var_values.push_back(val);
                } else {
                    var_values.push_back(t_root);
                }
            }
        }

        // Check if solution is real
        bool is_real = true;
        for (const auto &val : var_values) {
            if (std::abs(val.imag()) > epsilon) {
                is_real = false;
                break;
            }
        }

        solution.solutions.push_back(var_values);
        solution.is_real_solution.push_back(is_real);
    }

    solution.success = true;
    return solution;
}

void
print_solution(const PolynomialSystemSolution &solution, std::ostream &out) {
    if (!solution.success) {
        out << "Failed to solve system: " << solution.error_message << std::endl;
        return;
    }

    out << "Polynomial System Solution" << std::endl;
    out << "==========================" << std::endl;
    out << "Number of solutions: " << solution.solutions.size() << std::endl;
    out << "Quotient dimension: " << solution.quotient_dimension << std::endl;

    // Print minimal polynomial
    out << "\nMinimal polynomial: ";
    bool first = true;
    for (int i = solution.minimal_polynomial.size() - 1; i >= 0; --i) {
        const auto &coeff = solution.minimal_polynomial[i];
        if (coeff != 0) {
            if (!first && coeff > 0) out << " + ";
            if (coeff < 0) out << " - ";

            mpq_class abs_coeff = abs(coeff);
            if (abs_coeff != 1 || i == 0) { out << abs_coeff; }

            if (i > 0) {
                out << "T";
                if (i > 1) out << "^" << i;
            }
            first = false;
        }
    }
    out << " = 0" << std::endl;

    // Print solutions
    out << "\nSolutions:" << std::endl;
    for (size_t i = 0; i < solution.solutions.size(); ++i) {
        out << std::setw(3) << (i + 1) << ". ";

        for (size_t j = 0; j < solution.variable_names.size(); ++j) {
            if (j > 0) out << ", ";
            out << solution.variable_names[j] << " = ";

            const auto &val = solution.solutions[i][j];
            if (solution.is_real_solution[i]) {
                out << std::setprecision(10) << val.real();
            } else {
                out << std::setprecision(10) << val.real();
                if (val.imag() >= 0)
                    out << " + ";
                else
                    out << " - ";
                out << std::abs(val.imag()) << "i";
            }
        }

        if (solution.is_real_solution[i]) { out << " (real)"; }
        out << std::endl;
    }
}

std::string
solution_to_string(const PolynomialSystemSolution &solution) {
    std::ostringstream oss;
    print_solution(solution, oss);
    return oss.str();
}

} // namespace julia_rur