#include "hyperplane_sections.hpp"
#include "f4_monomial_decoder.hpp"
#include "f4_polynomial_formatter.hpp"
#include "polynomial_solver_enhanced.hpp"
#include "quotient_basis.hpp"
#include "prime_utils.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>
#include <unordered_set>
#include <set>

namespace julia_rur {

static bool
is_zero_dimensional_majority(const std::vector<std::string> &polynomials,
                             const std::vector<std::string> &variables,
                             int num_primes = 5) {
    using julia_rur::PP;
    if (num_primes <= 0) num_primes = 1;
    
    // Get random primes to avoid systematic unlucky primes
    auto primes = get_random_prime_sequence(num_primes, 28, 30);
    int zero_counts = 0;
    
    for (ModularCoeff p : primes) {
        // Compute GB and extract leading monomials
        std::vector<const char *> var_ptrs;
        var_ptrs.reserve(variables.size());
        for (const auto &v : variables) var_ptrs.push_back(v.c_str());
        axf4_session_t session = axf4_create_session(p, var_ptrs.data(), variables.size());
        if (!session) { break; }

        bool ok = true;
        for (const auto &poly : polynomials) {
            std::string formatted = format_polynomial_for_f4(poly);
            if (axf4_add_polynomial(session, formatted.c_str()) != 0) {
                ok = false;
                break;
            }
        }
        if (ok) {
            axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
            ok = (gb_result.status == 0);
            if (ok) {
                int basis_size = axf4_get_basis_size();
                std::vector<unsigned int> leading_monomials;
                if (basis_size > 0) {
                    leading_monomials.resize(basis_size);
                    int num_leading = axf4_get_all_leading_monomials(leading_monomials.data());
                    ok = (num_leading > 0);
                } else {
                    ok = false;
                }

                if (ok) {
                    int num_vars = axf4_get_num_variables(session);
                    std::vector<PP> leading_terms_pp;
                    leading_terms_pp.reserve(leading_monomials.size());
                    for (unsigned int enc : leading_monomials) {
                        auto exps = decode_f4_monomial(enc, num_vars);
                        leading_terms_pp.push_back(exps);
                    }
                    if (is_zero_dimensional(leading_terms_pp, num_vars)) { zero_counts++; }
                }
            }
            axf4_free_result(&gb_result);
            axf4_cleanup_basis_data();
        }
        axf4_destroy_session(session);
    }
    return zero_counts * 2 >= primes.size(); // majority
}

DimensionAnalysis
analyze_system_dimension(const std::vector<std::string> &polynomials,
                         const std::vector<std::string> &variables,
                         ModularCoeff prime) {
    DimensionAnalysis result;
    result.dimension = -1;
    result.codimension = -1;
    result.is_zero_dimensional = false;
    result.method_used = "none";

    if (polynomials.empty() || variables.empty()) { return result; }
    
    // Special case: single polynomial systems
    // The F4 algorithm has issues with these, so handle them directly
    if (polynomials.size() == 1) {
        // For a single polynomial, check if it has univariate leading terms
        // This is a simple heuristic: count variables that appear
        std::vector<bool> var_appears(variables.size(), false);
        for (size_t i = 0; i < variables.size(); i++) {
            if (polynomials[0].find(variables[i]) != std::string::npos) {
                var_appears[i] = true;
            }
        }
        
        int vars_in_poly = 0;
        for (bool appears : var_appears) {
            if (appears) vars_in_poly++;
        }
        
        // Single polynomial with multiple variables = positive dimensional
        if (vars_in_poly > 1) {
            result.dimension = vars_in_poly - 1;  // Dimension = #vars - #equations
            result.codimension = 1;
            result.is_zero_dimensional = false;
            result.method_used = "single_polynomial_heuristic";
            return result;
        }
    }
    
    // Use random prime if not specified
    if (prime == 0) {
        prime = generate_random_prime(28, 30);
    }
    
    // Debug: print which prime we're using
    if (variables.size() <= 3) {
        std::cout << "[DIM DEBUG] Using prime=" << prime << std::endl;
    }

    // Create F4 session to compute Gröbner basis
    std::vector<const char *> var_ptrs;
    for (const auto &var : variables) { var_ptrs.push_back(var.c_str()); }

    axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), variables.size());
    if (!session) { return result; }

    // Add polynomials
    for (const auto &poly : polynomials) {
        std::string formatted = format_polynomial_for_f4(poly);
        if (variables.size() <= 3) {
            std::cout << "[DIM DEBUG] Adding polynomial: " << formatted << std::endl;
        }
        if (axf4_add_polynomial(session, formatted.c_str()) != 0) {
            axf4_destroy_session(session);
            return result;
        }
    }

    // Compute Gröbner basis
    axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
    if (gb_result.status != 0) {
        axf4_free_result(&gb_result);
        axf4_destroy_session(session);
        return result;
    }

    // Get number of variables
    int32_t num_vars = axf4_get_num_variables(session);

    // Get basis size
    int basis_size = axf4_get_basis_size();
    
    // Debug output for small systems
    if (variables.size() <= 3) {
        std::cout << "[DIM DEBUG] GB has " << basis_size << " polynomial(s)" << std::endl;
    }
    
    // Handle empty basis (should not happen for non-trivial input)
    if (basis_size == 0) {
        result.dimension = -1;
        result.codimension = -1;
        result.is_zero_dimensional = false;
        result.method_used = "empty_gb";
        axf4_free_result(&gb_result);
        axf4_cleanup_basis_data();
        axf4_destroy_session(session);
        return result;
    }

    // Collect leading terms as PP vectors for proper dimension analysis
    std::vector<PP> leading_terms;
    bool extracted_any_lt = false;
    bool gb_is_one = false;  // Check for GB = {1} case
    
    if (basis_size > 0) {
        // Get all leading monomials
        std::vector<unsigned int> leading_monomials(basis_size);
        int num_leading = axf4_get_all_leading_monomials(leading_monomials.data());

        if (num_leading > 0) {
            // Convert encoded monomials to PP format
            for (int i = 0; i < num_leading; ++i) {
                auto exponents = decode_f4_monomial(leading_monomials[i], num_vars);
                leading_terms.push_back(exponents);
                
                // Check if any exponent is non-zero
                bool has_nonzero = false;
                for (int exp : exponents) {
                    if (exp > 0) {
                        has_nonzero = true;
                        extracted_any_lt = true;
                        break;
                    }
                }
                
                // For debugging: print leading term
                if (variables.size() <= 3 && i < 3) {
                    std::cout << "[DIM DEBUG] LT[" << i << "] has_nonzero=" << has_nonzero 
                              << " exps=[";
                    for (int j = 0; j < num_vars; ++j) {
                        if (j > 0) std::cout << ",";
                        std::cout << exponents[j];
                    }
                    std::cout << "]" << std::endl;
                }
            }
            
            // Check for GB = {1} after processing all leading terms
            // GB = {1} means exactly one polynomial with a constant (all zero exponents) leading term
            // BUT: Over small finite fields, some polynomials can reduce to constants
            // So we should be more careful here
            if (basis_size == 1 && num_leading == 1 && !extracted_any_lt) {
                // For small primes, this might be a field-specific artifact
                // We should try with a different prime rather than concluding GB={1}
                if (prime < 100) {
                    // Small prime - likely field-specific behavior
                    result.dimension = -1;  // Unknown
                    result.codimension = -1;
                    result.is_zero_dimensional = false;
                    result.method_used = "small_prime_artifact";
                    if (variables.size() <= 3) {
                        std::cout << "[DIM DEBUG] Small prime " << prime 
                                  << " gave constant - likely field artifact" << std::endl;
                    }
                    axf4_free_result(&gb_result);
                    axf4_cleanup_basis_data();
                    axf4_destroy_session(session);
                    return result;
                }
                gb_is_one = true;
                if (variables.size() <= 3) {
                    std::cout << "[DIM DEBUG] Detected GB={1} (constant polynomial)" << std::endl;
                }
            }
        }
    }
    
    // Special case: GB = {1} means no solutions mod this prime (unlucky prime)
    // This is NOT a positive-dimensional system!
    if (gb_is_one) {
        result.dimension = -1;  // Unknown/indeterminate
        result.codimension = -1;
        result.is_zero_dimensional = false;  // Can't determine from this prime
        result.method_used = "gb_equals_one_unlucky_prime";
        axf4_free_result(&gb_result);
        axf4_cleanup_basis_data();
        axf4_destroy_session(session);
        return result;
    }

    // Use the proper is_zero_dimensional check from quotient_basis.cpp
    // This checks if each variable has a univariate leading term
    bool is_zero_dim = false;
    if (!leading_terms.empty()) {
        is_zero_dim = is_zero_dimensional(leading_terms, num_vars);
        // Debug: print what we found
        if (variables.size() <= 3) {  // Only for small systems
            std::cout << "[DIM DEBUG] Found " << leading_terms.size() << " leading terms for " 
                      << num_vars << " vars, is_zero_dim=" << is_zero_dim << std::endl;
            for (size_t i = 0; i < leading_terms.size() && i < 5; ++i) {
                std::cout << "  LT[" << i << "] = [";
                for (int j = 0; j < num_vars; ++j) {
                    if (j > 0) std::cout << ",";
                    std::cout << leading_terms[i][j];
                }
                std::cout << "]" << std::endl;
            }
        }
    }
    
    // For positive-dimensional systems, identify free variables
    // (variables without univariate leading terms)
    result.free_variables.clear();
    if (!is_zero_dim) {
        std::vector<bool> has_univariate(num_vars, false);
        for (const auto& lt : leading_terms) {
            int nonzero_count = 0;
            int nonzero_var = -1;
            for (int i = 0; i < num_vars; i++) {
                if (lt[i] > 0) {
                    nonzero_count++;
                    nonzero_var = i;
                }
            }
            if (nonzero_count == 1) {
                has_univariate[nonzero_var] = true;
            }
        }
        
        // Collect free variables (those without univariate LT)
        std::set<int> free_var_set;
        for (int i = 0; i < num_vars; ++i) {
            if (!has_univariate[i]) {
                result.free_variables.push_back(i);
                free_var_set.insert(i);
            }
        }
        
        // IMPORTANT: Check if the free variables are truly independent
        // If there's a leading term that consists only of the "free" variables,
        // then they're algebraically dependent and the dimension is less
        if (result.free_variables.size() > 1) {
            bool has_dependency = false;
            for (const auto& lt : leading_terms) {
                bool uses_only_free = true;
                bool uses_at_least_one = false;
                for (int i = 0; i < num_vars; ++i) {
                    if (lt[i] > 0) {
                        if (free_var_set.find(i) == free_var_set.end()) {
                            uses_only_free = false;
                            break;
                        } else {
                            uses_at_least_one = true;
                        }
                    }
                }
                if (uses_only_free && uses_at_least_one) {
                    has_dependency = true;
                    break;
                }
            }
            
            // If there's dependency among free variables, the dimension is overestimated
            // For systems like Cyclic-4 with 2 free variables and dependency, dimension = 1
            if (has_dependency) {
                // As a heuristic: if we have exactly 2 free variables with dependency,
                // the dimension is likely 1 (common case for Cyclic-n)
                if (result.free_variables.size() == 2) {
                    // Override dimension to 1
                    int actual_dimension = 1;
                    if (variables.size() <= 3) {  // Debug output for small systems
                        std::cout << "[DIM DEBUG] Detected algebraic dependency among free variables. "
                                  << "Adjusting dimension from " << result.free_variables.size() 
                                  << " to " << actual_dimension << std::endl;
                    }
                    // Keep first free variable as representative
                    result.free_variables.resize(1);
                } else {
                    // For more complex cases, be conservative and reduce by 1
                    int adjusted_dim = static_cast<int>(result.free_variables.size()) - 1;
                    result.free_variables.resize(adjusted_dim);
                }
            }
        }
    }

    // If we couldn't extract any leading terms, mark as unknown
    if (!extracted_any_lt) {
        result.dimension = -1;
        result.codimension = -1;
        result.is_zero_dimensional = false;
        result.method_used = "leading_terms_unavailable";
        axf4_free_result(&gb_result);
        axf4_cleanup_basis_data();
        axf4_destroy_session(session);
        return result;
    }

    // Set the dimension based on the proper zero-dimensional check
    if (is_zero_dim) {
        result.dimension = 0;
        result.codimension = 0;
        result.is_zero_dimensional = true;
        result.method_used = "univariate_leading_terms";
    } else {
        // For positive-dimensional: dimension = number of free variables
        // (variables without univariate leading terms)
        result.dimension = static_cast<int>(result.free_variables.size());
        result.codimension = result.dimension; // Number of hyperplanes needed
        result.is_zero_dimensional = false;
        result.method_used = "free_variables_count";
    }

    // Cleanup
    axf4_free_result(&gb_result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);

    return result;
}

DimensionAnalysis
analyze_system_dimension_robust(const std::vector<std::string> &polynomials,
                                const std::vector<std::string> &variables,
                                const HyperplaneSectionConfig &config) {
    // 1) Multi-prime majority using leading-terms criterion
    // Use random primes to avoid systematic unlucky primes
    auto primes = get_random_prime_sequence(std::max(1, config.dimension_num_primes), 28, 30);

    std::vector<int> dims;
    dims.reserve(primes.size() * 3);
    std::mt19937 perm_rng(std::random_device{}());
    for (auto prime : primes) {
        // Try original order + two random permutations
        {
            auto d = analyze_system_dimension(polynomials, variables, prime);
            if (d.dimension >= 0) { dims.push_back(d.dimension); }
        }
        for (int t = 0; t < 2; ++t) {
            std::vector<std::string> perm_vars = variables;
            std::shuffle(perm_vars.begin(), perm_vars.end(), perm_rng);
            auto d = analyze_system_dimension(polynomials, perm_vars, prime);
            if (d.dimension >= 0) { dims.push_back(d.dimension); }
        }
    }

    if (!dims.empty()) {
        std::sort(dims.begin(), dims.end());
        int best_dim = dims[0], best_count = 1, cur_dim = dims[0], cur_count = 1;
        for (size_t i = 1; i < dims.size(); ++i) {
            if (dims[i] == cur_dim) {
                cur_count++;
            } else {
                if (cur_count > best_count) {
                    best_count = cur_count;
                    best_dim = cur_dim;
                }
                cur_dim = dims[i];
                cur_count = 1;
            }
        }
        if (cur_count > best_count) {
            best_count = cur_count;
            best_dim = cur_dim;
        }

        if (best_count * 2 >= static_cast<int>(primes.size())) {
            DimensionAnalysis out;
            out.dimension = best_dim;
            out.codimension = best_dim;
            out.is_zero_dimensional = (best_dim == 0);
            out.method_used = "leading_terms_multi_prime_majority";
            return out;
        }
    }

    // 2) Slicing fallback
    for (int s = 0; s <= static_cast<int>(variables.size()); ++s) {
        bool zero_dim_found = false;
        for (int trial = 0; trial < std::max(1, config.dimension_slicing_retries); ++trial) {
            HyperplaneSectionConfig hcfg = config;
            hcfg.verbose = false;
            std::vector<std::string> augmented = polynomials;
            std::random_device rd;
            std::mt19937 rng(rd());
            for (int i = 0; i < s; ++i) {
                std::string hp = generate_random_hyperplane(variables, hcfg, rng);
                augmented.push_back(hp);
            }

            auto d = analyze_system_dimension(augmented, variables);
            if (d.dimension == 0 && d.is_zero_dimensional) {
                zero_dim_found = true;
                break;
            }
        }
        if (zero_dim_found) {
            DimensionAnalysis out;
            out.dimension = s;
            out.codimension = s;
            out.is_zero_dimensional = (s == 0);
            out.method_used = "slicing_fallback";
            return out;
        }
    }

    DimensionAnalysis fail;
    fail.dimension = -1;
    fail.codimension = -1;
    fail.is_zero_dimensional = false;
    fail.method_used = "dimension_unknown";
    return fail;
}

std::string
generate_random_hyperplane(const std::vector<std::string> &variables,
                           const HyperplaneSectionConfig &config,
                           std::mt19937 &rng) {
    std::uniform_int_distribution<int> coeff_dist(config.min_coefficient, config.max_coefficient);

    std::stringstream hyperplane;
    bool first = true;
    bool has_nonzero = false;

    // Generate coefficients for each variable
    for (size_t i = 0; i < variables.size(); ++i) {
        int coeff = coeff_dist(rng);

        // Ensure at least one non-zero coefficient
        if (i == variables.size() - 1 && !has_nonzero) {
            while (coeff == 0) { coeff = coeff_dist(rng); }
        }

        if (config.include_all_variables) {
            // Force non-zero coefficient for every variable
            while (coeff == 0) { coeff = coeff_dist(rng); }
        } else if (config.avoid_zero_coefficients && coeff == 0 && i < variables.size() - 1) {
            // Retry once to avoid zero when not forcing all variables
            coeff = coeff_dist(rng);
        }

        if (coeff != 0) {
            has_nonzero = true;

            if (!first && coeff > 0) { hyperplane << "+"; }

            if (coeff == -1) {
                hyperplane << "-";
            } else if (coeff != 1) {
                hyperplane << coeff << "*";
            }

            hyperplane << variables[i];
            first = false;
        }
    }

    // Add constant term
    int constant = coeff_dist(rng);
    // Avoid purely homogeneous hyperplanes to reduce degeneracy
    while (constant == 0) { constant = coeff_dist(rng); }
    if (constant != 0) {
        if (constant > 0 && !first) { hyperplane << "+"; }
        hyperplane << constant;
    }

    // If we somehow got all zeros (shouldn't happen), create x1
    if (!has_nonzero) { return variables[0]; }

    return hyperplane.str();
}

HyperplaneSectionResult
solve_via_hyperplane_sections(const std::vector<std::string> &polynomials,
                              const std::vector<std::string> &variables,
                              const HyperplaneSectionConfig &config) {
    HyperplaneSectionResult result;
    result.success = false;

    if (config.verbose) {
        std::cout << "[slice] START force=" << (config.force_slicing ? 1 : 0) << " m=" << polynomials.size()
                  << " n=" << variables.size() << std::endl;
    }

    // First, analyze dimension (fast single-prime to avoid instability)
    result.dimension_info = analyze_system_dimension(polynomials, variables);

    if (config.verbose) {
        std::cout << "[slice] enter force=" << (config.force_slicing ? 1 : 0)
                  << " dim=" << result.dimension_info.dimension
                  << " zero=" << (result.dimension_info.is_zero_dimensional ? 1 : 0)
                  << " method=" << result.dimension_info.method_used << std::endl;
    }

    // If dimension analysis failed, still proceed with slicing starting at s=1

    if (result.dimension_info.is_zero_dimensional && !config.force_slicing) {
        // System appears zero-dimensional; try direct solve first
        result.solutions = solve_polynomial_system_enhanced(polynomials, variables, config);
        if (result.solutions.success) {
            result.success = true;
            result.augmented_system = polynomials;
            return result;
        }
        // Fall back to adding at least one hyperplane if direct solve fails
        if (config.verbose) {
            std::cout << "Direct zero-dimensional solve failed; falling back to hyperplane sections..." << std::endl;
        }
        // Treat as needing one hyperplane to cut potential degeneracy
        result.dimension_info.codimension = 1;
    }

    if (result.dimension_info.dimension > config.max_dimension) {
        result.error_message = "System dimension " + std::to_string(result.dimension_info.dimension) +
                               " exceeds maximum " + std::to_string(config.max_dimension);
        return result;
    }

    // Try multiple times to find consistent hyperplanes
    std::random_device rd;
    std::mt19937 rng(rd());
    const int max_retries = 6;

    int max_s = static_cast<int>(variables.size());
    int start_s = 1;
    if (result.dimension_info.dimension >= 0) { start_s = std::max(1, result.dimension_info.codimension); }
    for (int s = start_s; s <= max_s; ++s) {
        for (int retry = 0; retry < max_retries; ++retry) {
            result.augmented_system = polynomials;
            result.hyperplanes.clear();

            // Generate s random hyperplanes
            for (int i = 0; i < s; ++i) {
                std::string hyperplane = generate_random_hyperplane(variables, config, rng);
                result.hyperplanes.push_back(hyperplane);
                result.augmented_system.push_back(hyperplane);

                if (config.verbose) {
                    std::cout << "[slice] s=" << s << " retry=" << retry << " hp#" << (i + 1) << " : " << hyperplane
                              << std::endl;
                }
            }

            // Robust zero-dimensionality check on augmented system (multi-prime majority, no slicing fallback)
            // Prefer zero-dimensional slices, but attempt solve even if certificate is inconclusive
            bool is_zero_dim = is_zero_dimensional_majority(result.augmented_system, variables, 5);
            if (config.verbose) { std::cout << "[slice] cert_zero_dim=" << (is_zero_dim ? 1 : 0) << std::endl; }

            // Solve the augmented system without further hyperplane recursion
            EnhancedSolverConfig inner_cfg = config;
            inner_cfg.auto_hyperplane_sections = false;
            inner_cfg.fail_on_positive_dimensional = true;
            inner_cfg.skip_dimension_precheck = true; // we already checked/are enforcing via slicing
            if (config.verbose) {
                std::cout << "[slice] calling RUR on augmented system (m=" << result.augmented_system.size()
                          << ", n=" << variables.size() << ")" << std::endl;
                std::cout << "[slice] augmented: ";
                for (size_t i = 0; i < result.augmented_system.size(); ++i) {
                    if (i) std::cout << "; ";
                    std::cout << result.augmented_system[i];
                }
                std::cout << std::endl;
            }
            result.solutions = solve_polynomial_system_enhanced(result.augmented_system, variables, inner_cfg);

            if (result.solutions.success) {
                result.success = true;
                result.error_message.clear();
                return result;
            }

            const std::string &msg = result.solutions.error_message;
            if (config.verbose) { std::cout << "[slice] fail: " << msg << std::endl; }
            bool inconsistent = msg.find("Gröbner basis is {1}") != std::string::npos;
            bool not_zerodim = msg.find("not define zero-dimensional ideal") != std::string::npos;
            bool rur_failed = msg.find("Failed to compute reference RUR") != std::string::npos;
            if (!(inconsistent || not_zerodim || rur_failed)) {
                // Non-retriable
                result.success = false;
                result.error_message = msg;
                return result;
            }
            // Otherwise, retry with new random hyperplanes
        }
    }

    // Failed after max retries
    result.success = false;
    result.error_message = "Failed to find consistent hyperplanes after " + std::to_string(max_retries) + " attempts";
    return result;
}

std::vector<HyperplaneSectionResult>
sample_variety_points(const std::vector<std::string> &polynomials,
                      const std::vector<std::string> &variables,
                      int num_samples,
                      const HyperplaneSectionConfig &config) {
    std::vector<HyperplaneSectionResult> results;

    for (int i = 0; i < num_samples; ++i) {
        auto sample_config = config;
        sample_config.verbose = config.verbose && (i == 0); // Only verbose for first sample

        HyperplaneSectionResult sample = solve_via_hyperplane_sections(polynomials, variables, sample_config);

        results.push_back(sample);

        if (config.verbose && sample.success) {
            std::cout << "Sample " << (i + 1) << " found " << sample.solutions.solutions.size() << " points"
                      << std::endl;
        }
    }

    return results;
}

HyperplaneSectionResult
solve_polynomial_system_adaptive(const std::vector<std::string> &polynomials,
                                 const std::vector<std::string> &variables,
                                 const HyperplaneSectionConfig &config) {
    // This is now the main entry point that handles both cases
    return solve_via_hyperplane_sections(polynomials, variables, config);
}

} // namespace julia_rur