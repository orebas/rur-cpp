#include "hyperplane_sections.hpp"
#include "polynomial_solver_enhanced.hpp"
#include "f4_polynomial_formatter.hpp"
#include "f4_monomial_decoder.hpp"
#include "f4_integration.hpp"
#include <algorithm>
#include <unordered_set>
#include <sstream>
#include <iostream>

namespace julia_rur {

DimensionAnalysis analyze_system_dimension(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    ModularCoeff prime
) {
    DimensionAnalysis result;
    result.dimension = -1;
    result.codimension = -1;
    result.is_zero_dimensional = false;
    result.method_used = "none";
    
    if (polynomials.empty() || variables.empty()) {
        return result;
    }
    
    // Create F4 session to compute Gröbner basis
    std::vector<const char*> var_ptrs;
    for (const auto& var : variables) {
        var_ptrs.push_back(var.c_str());
    }
    
    axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), variables.size());
    if (!session) {
        return result;
    }
    
    // Add polynomials
    for (const auto& poly : polynomials) {
        std::string formatted = format_polynomial_for_f4(poly);
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
    
    // Extract leading monomials to determine which variables appear
    std::unordered_set<int> leading_variables;
    
    // Get number of variables
    int32_t num_vars = axf4_get_num_variables(session);
    
    // Get basis size
    int basis_size = axf4_get_basis_size();
    
    if (basis_size > 0) {
        // Get all leading monomials
        std::vector<unsigned int> leading_monomials(basis_size);
        int num_leading = axf4_get_all_leading_monomials(leading_monomials.data());
        
        if (num_leading > 0) {
            // Analyze which variables appear in leading terms
            // Decode the monomials to find variable indices
            for (int i = 0; i < num_leading; ++i) {
                // Use the decode_f4_monomial function directly
                auto exponents = decode_f4_monomial(leading_monomials[i], num_vars);
                for (int var_idx = 0; var_idx < num_vars; ++var_idx) {
                    if (exponents[var_idx] > 0) {
                        leading_variables.insert(var_idx);
                    }
                }
            }
        }
    }
    
    // Variables that don't appear in any leading term
    result.free_variables.clear();
    for (int i = 0; i < num_vars; ++i) {
        if (leading_variables.find(i) == leading_variables.end()) {
            result.free_variables.push_back(i);
        }
    }
    
    // The dimension of the variety is:
    // dim(V) = #variables - #equations_in_generic_position
    // For an ideal in generic position, this equals #free_variables
    // However, we need to be more careful...
    
    // For now, use a simple heuristic:
    // If we have n variables and k polynomials with k < n,
    // the expected dimension is n - k (assuming generic position)
    int expected_dim = std::max(0, static_cast<int>(variables.size()) - static_cast<int>(polynomials.size()));
    
    // The actual dimension based on leading terms
    int leading_term_dim = result.free_variables.size();
    
    // Use the minimum of the two as a conservative estimate
    result.dimension = std::min(expected_dim, leading_term_dim);
    result.codimension = result.dimension;  // Number of hyperplanes needed
    result.is_zero_dimensional = (result.dimension == 0);
    result.method_used = "leading_terms";
    
    // Cleanup
    axf4_free_result(&gb_result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);
    
    return result;
}

std::string generate_random_hyperplane(
    const std::vector<std::string>& variables,
    const HyperplaneSectionConfig& config,
    std::mt19937& rng
) {
    std::uniform_int_distribution<int> coeff_dist(config.min_coefficient, config.max_coefficient);
    
    std::stringstream hyperplane;
    bool first = true;
    bool has_nonzero = false;
    
    // Generate coefficients for each variable
    for (size_t i = 0; i < variables.size(); ++i) {
        int coeff = coeff_dist(rng);
        
        // Ensure at least one non-zero coefficient
        if (i == variables.size() - 1 && !has_nonzero) {
            while (coeff == 0) {
                coeff = coeff_dist(rng);
            }
        }
        
        if (config.avoid_zero_coefficients && coeff == 0 && i < variables.size() - 1) {
            // Retry once to avoid zero
            coeff = coeff_dist(rng);
        }
        
        if (coeff != 0) {
            has_nonzero = true;
            
            if (!first && coeff > 0) {
                hyperplane << "+";
            }
            
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
    if (constant != 0) {
        if (constant > 0 && !first) {
            hyperplane << "+";
        }
        hyperplane << constant;
    }
    
    // If we somehow got all zeros (shouldn't happen), create x1
    if (!has_nonzero) {
        return variables[0];
    }
    
    return hyperplane.str();
}

HyperplaneSectionResult solve_via_hyperplane_sections(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const HyperplaneSectionConfig& config
) {
    HyperplaneSectionResult result;
    result.success = false;
    
    // First, analyze dimension
    result.dimension_info = analyze_system_dimension(polynomials, variables);
    
    if (result.dimension_info.dimension < 0) {
        result.error_message = "Failed to analyze system dimension";
        return result;
    }
    
    if (result.dimension_info.is_zero_dimensional) {
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
        result.error_message = "System dimension " + std::to_string(result.dimension_info.dimension) 
                             + " exceeds maximum " + std::to_string(config.max_dimension);
        return result;
    }
    
    // Try multiple times to find consistent hyperplanes
    std::random_device rd;
    std::mt19937 rng(rd());
    const int max_retries = 10;
    
    for (int retry = 0; retry < max_retries; ++retry) {
        result.augmented_system = polynomials;
        result.hyperplanes.clear();
        
        // Generate random hyperplanes
        for (int i = 0; i < result.dimension_info.codimension; ++i) {
            std::string hyperplane = generate_random_hyperplane(variables, config, rng);
            result.hyperplanes.push_back(hyperplane);
            result.augmented_system.push_back(hyperplane);
            
            if (config.verbose && retry == 0) {
                std::cout << "Added hyperplane " << (i+1) << ": " << hyperplane << std::endl;
            }
        }
        
        // Solve the augmented system
        result.solutions = solve_polynomial_system_enhanced(result.augmented_system, variables, config);

        // Retry conditions on failure
        if (!result.solutions.success) {
            const std::string &msg = result.solutions.error_message;
            bool inconsistent = msg.find("Gröbner basis is {1}") != std::string::npos;
            bool not_zerodim = msg.find("not define zero-dimensional ideal") != std::string::npos;
            bool rur_failed = msg.find("Failed to compute reference RUR") != std::string::npos;
            if (inconsistent || not_zerodim || rur_failed) {
                if (config.verbose) {
                    std::cout << "Retry " << (retry + 1)
                              << ": Augmented system failed (" << msg
                              << "), trying different hyperplanes..." << std::endl;
                }
                continue;
            }
        }

        // Either success or non-retriable error - accept the result
        result.success = result.solutions.success;
        result.error_message = result.solutions.error_message;
        return result;
    }
    
    // Failed after max retries
    result.success = false;
    result.error_message = "Failed to find consistent hyperplanes after " + std::to_string(max_retries) + " attempts";
    return result;
}

std::vector<HyperplaneSectionResult> sample_variety_points(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    int num_samples,
    const HyperplaneSectionConfig& config
) {
    std::vector<HyperplaneSectionResult> results;
    
    for (int i = 0; i < num_samples; ++i) {
        auto sample_config = config;
        sample_config.verbose = config.verbose && (i == 0);  // Only verbose for first sample
        
        HyperplaneSectionResult sample = solve_via_hyperplane_sections(
            polynomials, variables, sample_config
        );
        
        results.push_back(sample);
        
        if (config.verbose && sample.success) {
            std::cout << "Sample " << (i+1) << " found " 
                      << sample.solutions.solutions.size() << " points" << std::endl;
        }
    }
    
    return results;
}

HyperplaneSectionResult solve_polynomial_system_adaptive(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const HyperplaneSectionConfig& config
) {
    // This is now the main entry point that handles both cases
    return solve_via_hyperplane_sections(polynomials, variables, config);
}

} // namespace julia_rur