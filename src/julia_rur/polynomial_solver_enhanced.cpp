#include "polynomial_solver_enhanced.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "hyperplane_sections.hpp"
#include <map>


namespace julia_rur {

RootFindingMethod
analyze_polynomial_for_method(const std::vector<mpq_class> &polynomial_coeffs, const EnhancedSolverConfig &config) {
    if (config.root_method != RootFindingMethod::AUTO) { return config.root_method; }

    // Find degree
    int degree = -1;
    for (int i = polynomial_coeffs.size() - 1; i >= 0; --i) {
        if (polynomial_coeffs[i] != 0) {
            degree = i;
            break;
        }
    }

    // High degree polynomials benefit from FLINT's robustness
    if (degree > config.degree_threshold_for_flint) { return RootFindingMethod::FLINT; }

    // Check coefficient sizes - large coefficients need arbitrary precision
    bool has_large_coeffs = false;
    for (const auto &coeff : polynomial_coeffs) {
        if (coeff != 0) {
            // Check if numerator or denominator is large
            if (coeff.get_num().get_str().length() > 15 || coeff.get_den().get_str().length() > 15) {
                has_large_coeffs = true;
                break;
            }
        }
    }

    if (has_large_coeffs || config.prefer_certified_roots) { return RootFindingMethod::FLINT; }

    // For moderate degree with reasonable coefficients, Eigen is faster
    return RootFindingMethod::EIGEN;
}

std::vector<std::complex<double>>
find_roots_auto(const std::vector<mpq_class> &polynomial_coeffs,
                const EnhancedSolverConfig &config,
                std::optional<std::vector<double>> &error_bounds) {
    RootFindingMethod method = analyze_polynomial_for_method(polynomial_coeffs, config);

    if (method == RootFindingMethod::FLINT) {
        std::vector<double> bounds;
        auto roots = find_polynomial_roots_certified(polynomial_coeffs, bounds, config.flint_config);
        error_bounds = bounds;
        return roots;
    } else {
        // Use Eigen
        error_bounds = std::nullopt;
        return find_polynomial_roots(polynomial_coeffs);
    }
}

void
analyze_root_structure(const std::vector<std::complex<double>> &roots,
                       EnhancedPolynomialSolution &solution,
                       const EnhancedSolverConfig &config) {
    const double eps = config.imaginary_tolerance;
    const double cluster_eps = config.clustering_threshold;

    // Clear previous analysis
    solution.root_multiplicities.clear();
    solution.conjugate_pairs.clear();
    solution.root_multiplicities.resize(roots.size(), 1);

    // Find conjugate pairs
    std::vector<bool> paired(roots.size(), false);

    for (size_t i = 0; i < roots.size(); ++i) {
        if (paired[i]) continue;

        // Check if this root has non-zero imaginary part
        if (std::abs(roots[i].imag()) > eps) {
            // Look for its conjugate
            for (size_t j = i + 1; j < roots.size(); ++j) {
                if (paired[j]) continue;

                // Check if j is conjugate of i
                double real_diff = std::abs(roots[i].real() - roots[j].real());
                double imag_sum = std::abs(roots[i].imag() + roots[j].imag());

                if (real_diff < cluster_eps && imag_sum < cluster_eps) {
                    solution.conjugate_pairs.push_back({ static_cast<int>(i), static_cast<int>(j) });
                    paired[i] = true;
                    paired[j] = true;
                    break;
                }
            }
        }
    }

    // Estimate multiplicities by clustering nearby roots
    std::vector<bool> counted(roots.size(), false);

    for (size_t i = 0; i < roots.size(); ++i) {
        if (counted[i]) continue;

        int multiplicity = 1;
        counted[i] = true;

        // Find all roots close to this one
        for (size_t j = i + 1; j < roots.size(); ++j) {
            if (counted[j]) continue;

            double dist = std::abs(roots[i] - roots[j]);
            if (dist < cluster_eps) {
                multiplicity++;
                counted[j] = true;
                solution.root_multiplicities[j] = 0; // Mark as part of cluster
            }
        }

        solution.root_multiplicities[i] = multiplicity;
    }
}

/**
 * @brief Evaluate polynomial parameterization at root value
 */
static std::complex<double>
evaluate_parameterization(const std::vector<mpq_class> &numerator,
                          const std::complex<double> &t_value,
                          const std::vector<mpq_class> &f_prime_coeffs,
                          bool verbose) {
    if (numerator.empty()) { return t_value; }

    auto eval_poly = [](const std::vector<mpq_class> &coeffs, const std::complex<double> &t) {
        std::complex<double> v = 0.0, tp = 1.0;
        for (const auto &c : coeffs) {
            v += c.get_d() * tp;
            tp *= t;
        }
        return v;
    };

    const double eps = 1e-12;

    if (verbose) {
        std::cout << "EVAL-DEBUG-ENH: t=" << std::setprecision(15) << t_value << std::endl;
        std::cout << "EVAL-DEBUG-ENH: numerator=[";
        for (size_t i = 0; i < numerator.size(); ++i) {
            if (i) std::cout << ", ";
            std::cout << numerator[i];
        }
        std::cout << "]" << std::endl;
    }

    // Evaluate f'(t)
    std::vector<mpq_class> f_deriv = f_prime_coeffs; // start from f'(T)
    std::complex<double> den = eval_poly(f_deriv, t_value);
    if (verbose) { std::cout << "EVAL-DEBUG-ENH: f'(t)=" << den << std::endl; }

    // Evaluate numerator at t
    std::vector<mpq_class> num = numerator;

    // Bound retries by combined degrees to avoid loops
    int max_steps = static_cast<int>(f_prime_coeffs.size() + numerator.size() + 2);
    for (int k = 0; k < max_steps; ++k) {
        std::complex<double> num_v = eval_poly(num, t_value);
        if (verbose) {
            std::cout << "EVAL-DEBUG-ENH: step=" << k << " num(t)=" << num_v << " den(t)=" << den << std::endl;
        }
        if (std::abs(den) > eps) {
            auto out = num_v / den;
            if (verbose) { std::cout << "EVAL-DEBUG-ENH: return num/den=" << out << std::endl; }
            return out;
        }

        // If both numerator and denominator vanish, apply one step of L'Hospital
        bool num_zeroish = std::abs(num_v) <= eps;
        if (num_zeroish) {
            if (num.size() > 1) {
                // Differentiate numerator symbolically: num <- num'
                std::vector<mpq_class> num_d;
                num_d.reserve(num.size() - 1);
                for (size_t i = 1; i < num.size(); ++i) { num_d.push_back(num[i] * static_cast<int>(i)); }
                num.swap(num_d);
                if (verbose) {
                    std::cout << "EVAL-DEBUG-ENH: LHospital numerator' deg=" << ((int)num.size() - 1) << std::endl;
                }
            }

            // Advance denominator derivative order: f^{(m)}(T) -> f^{(m+1)}(T)
            std::vector<mpq_class> f_next;
            if (f_deriv.size() > 1) {
                f_next.reserve(f_deriv.size() - 1);
                for (size_t i = 1; i < f_deriv.size(); ++i) { f_next.push_back(f_deriv[i] * static_cast<int>(i)); }
            }
            f_deriv.swap(f_next);
            den = f_deriv.empty() ? std::complex<double>(0.0, 0.0) : eval_poly(f_deriv, t_value);
            if (verbose) { std::cout << "EVAL-DEBUG-ENH: advanced denominator derivative, value=" << den << std::endl; }
            continue;
        }

        // num not ~0 but den ~0: ambiguous, return numerator value as last resort
        if (verbose) { std::cout << "EVAL-DEBUG-ENH: den~0 but num!=0, return num=" << num_v << std::endl; }
        return num_v;
    }
    // Fallback
    if (verbose) { std::cout << "EVAL-DEBUG-ENH: final fallback 0" << std::endl; }
    return std::complex<double>(0.0, 0.0);
}

EnhancedPolynomialSolution
solve_polynomial_system_enhanced(const std::vector<std::string> &polynomials,
                                 const std::vector<std::string> &variables,
                                 const EnhancedSolverConfig &config) {
    auto start_time = std::chrono::high_resolution_clock::now();

    EnhancedPolynomialSolution solution;
    solution.variable_names = variables;
    solution.success = false;

    // Step 0: Quick dimensionality pre-check
    solution.computed_dimension = -1;
    {
        DimensionAnalysis dim = analyze_system_dimension(polynomials, variables);
        if (dim.dimension >= 0) {
            solution.computed_dimension = dim.dimension;
            if (!dim.is_zero_dimensional) {
                solution.error_message = "Positive-dimensional ideal (dimension " +
                                         std::to_string(dim.dimension) +
                                         "). Use adaptive solver for hyperplane sections.";
                return solution;
            }
        }
    }

    // Step 1: Compute RUR
    RationalRURResult rur_result = compute_rational_rur(polynomials, variables, config);

    if (!rur_result.success) {
        solution.error_message = "RUR computation failed: " + rur_result.error_message;
        return solution;
    }

    solution.minimal_polynomial = rur_result.minimal_polynomial;
    solution.quotient_dimension = rur_result.quotient_basis.size();

    // Step 2: Find roots with enhanced method
    auto root_start = std::chrono::high_resolution_clock::now();

    std::optional<std::vector<double>> error_bounds;
    std::vector<std::complex<double>> t_roots = find_roots_auto(rur_result.minimal_polynomial, config, error_bounds);

    solution.method_used = analyze_polynomial_for_method(rur_result.minimal_polynomial, config);

    if (error_bounds.has_value()) {
        solution.root_error_bounds = error_bounds.value();
        solution.roots_certified = true;
    }

    auto root_end = std::chrono::high_resolution_clock::now();
    solution.root_finding_time_ms = std::chrono::duration<double, std::milli>(root_end - root_start).count();

    if (t_roots.empty()) {
        solution.error_message = "No roots found for minimal polynomial";
        return solution;
    }

    // Step 3: Compute f'(T) for parameterization denominators
    std::vector<mpq_class> derivative_coeffs;
    for (size_t i = 1; i < rur_result.minimal_polynomial.size(); ++i) {
        derivative_coeffs.push_back(rur_result.minimal_polynomial[i] * static_cast<int>(i));
    }

    // Step 4: Back-substitute to find variable values
    for (const auto &t_root : t_roots) {
        std::vector<std::complex<double>> var_values;

        if (config.verbose) {
            std::cout << "EVAL-DEBUG-ENH: solving for t_root=" << std::setprecision(15) << t_root << std::endl;
        }

        if (variables.size() == 1) {
            // Univariate case - check if simple
            if (rur_result.numerators.empty() || rur_result.numerators[0].empty() ||
                (rur_result.numerators[0].size() == 2 && rur_result.numerators[0][0] == 0 &&
                 rur_result.numerators[0][1] != 0)) {
                // Simple case: x = T
                var_values.push_back(t_root);
            } else {
                // Parameterized case
                var_values.push_back(
                  evaluate_parameterization(rur_result.numerators[0], t_root, derivative_coeffs, config.verbose));
            }
        } else {
            // Multivariate case
            for (size_t var_idx = 0; var_idx < variables.size(); ++var_idx) {
                if (var_idx < rur_result.numerators.size()) {
                    if (config.verbose) {
                        std::cout << "EVAL-DEBUG-ENH: var " << variables[var_idx] << " numerator=[";
                        for (size_t i = 0; i < rur_result.numerators[var_idx].size(); ++i) {
                            if (i) std::cout << ", ";
                            std::cout << rur_result.numerators[var_idx][i];
                        }
                        std::cout << "]" << std::endl;
                    }
                    var_values.push_back(evaluate_parameterization(
                      rur_result.numerators[var_idx], t_root, derivative_coeffs, config.verbose));
                } else {
                    // Fallback
                    var_values.push_back(t_root);
                }
            }
        }

        // Check if solution is real
        bool is_real = true;
        for (const auto &val : var_values) {
            if (std::abs(val.imag()) > config.imaginary_tolerance) {
                is_real = false;
                break;
            }
        }

        // Apply filtering if requested
        if (!config.real_roots_only || is_real) {
            solution.solutions.push_back(var_values);
            solution.is_real_solution.push_back(is_real);
        }
    }

    // Step 5: Analyze root structure
    analyze_root_structure(t_roots, solution, config);

    solution.success = true;

    auto end_time = std::chrono::high_resolution_clock::now();
    solution.total_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    return solution;
}

void
print_enhanced_solution(const EnhancedPolynomialSolution &solution, std::ostream &out) {
    if (!solution.success) {
        out << "Failed to solve system: " << solution.error_message << std::endl;
        return;
    }

    out << "Enhanced Polynomial System Solution" << std::endl;
    out << "===================================" << std::endl;
    out << "Method used: ";
    switch (solution.method_used) {
        case RootFindingMethod::EIGEN:
            out << "Eigen";
            break;
        case RootFindingMethod::FLINT:
            out << "FLINT (certified)";
            break;
        case RootFindingMethod::HYBRID:
            out << "Hybrid";
            break;
        default:
            out << "Auto";
            break;
    }
    out << std::endl;

    out << "Number of solutions: " << solution.solutions.size() << std::endl;
    out << "Quotient dimension: " << solution.quotient_dimension << std::endl;
    out << "Root finding time: " << std::fixed << std::setprecision(2) << solution.root_finding_time_ms << " ms"
        << std::endl;
    out << "Total time: " << solution.total_time_ms << " ms" << std::endl;

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

    // Print conjugate pairs if any
    if (!solution.conjugate_pairs.empty()) {
        out << "\nComplex conjugate pairs: ";
        for (const auto &pair : solution.conjugate_pairs) { out << "(" << pair[0] + 1 << "," << pair[1] + 1 << ") "; }
        out << std::endl;
    }

    // Print solutions with error bounds if available
    out << "\nSolutions:" << std::endl;
    for (size_t i = 0; i < solution.solutions.size(); ++i) {
        out << std::setw(3) << (i + 1) << ". ";

        for (size_t j = 0; j < solution.variable_names.size(); ++j) {
            if (j > 0) out << ", ";
            out << solution.variable_names[j] << " = ";

            const auto &val = solution.solutions[i][j];
            if (solution.is_real_solution[i]) {
                out << std::setprecision(15) << val.real();
            } else {
                out << std::setprecision(15) << val.real();
                if (val.imag() >= 0)
                    out << " + ";
                else
                    out << " - ";
                out << std::abs(val.imag()) << "i";
            }
        }

        if (solution.roots_certified && i < solution.root_error_bounds.size()) {
            out << " (±" << std::scientific << std::setprecision(2) << solution.root_error_bounds[i] << ")";
        }

        if (i < solution.root_multiplicities.size() && solution.root_multiplicities[i] > 1) {
            out << " [mult: " << solution.root_multiplicities[i] << "]";
        }

        if (solution.is_real_solution[i]) { out << " (real)"; }
        out << std::endl;
    }
}

} // namespace julia_rur