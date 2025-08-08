#pragma once

#include "rur_main_algorithm.hpp"
#include "polynomial_solver_enhanced.hpp"
#include <random>
#include <optional>

namespace julia_rur {

/**
 * @brief Result of dimension analysis
 */
struct DimensionAnalysis {
    int dimension;                    // -1 if analysis failed
    int codimension;                  // Number of hyperplanes needed
    std::vector<int> free_variables;  // Indices of variables not in leading terms
    bool is_zero_dimensional;
    std::string method_used;          // "leading_terms", "hilbert", etc.
};

/**
 * @brief Hyperplane section result
 */
struct HyperplaneSectionResult {
    bool success;
    
    // The augmented polynomial system (original + hyperplanes)
    std::vector<std::string> augmented_system;
    
    // The hyperplanes added
    std::vector<std::string> hyperplanes;
    
    // Solutions found
    EnhancedPolynomialSolution solutions;
    
    // Original dimension info
    DimensionAnalysis dimension_info;
    
    std::string error_message;
};

/**
 * @brief Configuration for hyperplane sections
 */
struct HyperplaneSectionConfig : public EnhancedSolverConfig {
    // Random hyperplane generation
    int min_coefficient = -100;
    int max_coefficient = 100;
    bool avoid_zero_coefficients = true;
    
    // Sampling strategy
    int num_sample_sets = 1;  // Number of different hyperplane sets to try
    bool include_all_variables = true;  // Ensure all variables appear in hyperplanes
    
    // Dimension detection
    bool use_gb_for_dimension = true;  // Use Gröbner basis for dimension
    int max_dimension = 10;  // Maximum dimension to handle
};

/**
 * @brief Analyze dimension of polynomial system
 * 
 * Uses leading terms of Gröbner basis to quickly determine dimension.
 * Variables that don't appear as leading terms indicate positive dimension.
 * 
 * @param polynomials System of polynomials
 * @param variables Variable names
 * @param prime Prime for modular computation
 * @return Dimension analysis result
 */
DimensionAnalysis analyze_system_dimension(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    ModularCoeff prime = 1073741827
);

/**
 * @brief Generate random hyperplane
 * 
 * Creates a hyperplane of the form: c₁x₁ + c₂x₂ + ... + cₙxₙ + c₀ = 0
 * with random integer coefficients from the specified range.
 * 
 * @param variables Variable names
 * @param config Configuration
 * @param rng Random number generator
 * @return Hyperplane as polynomial string
 */
std::string generate_random_hyperplane(
    const std::vector<std::string>& variables,
    const HyperplaneSectionConfig& config,
    std::mt19937& rng
);

/**
 * @brief Solve positive-dimensional system via hyperplane sections
 * 
 * Main entry point for handling underdetermined systems:
 * 1. Detects dimension
 * 2. Adds random hyperplanes to make system zero-dimensional
 * 3. Solves augmented system
 * 4. Returns intersection points and hyperplane info
 * 
 * @param polynomials Original polynomial system
 * @param variables Variable names
 * @param config Configuration options
 * @return Solutions on random hyperplane sections
 */
HyperplaneSectionResult solve_via_hyperplane_sections(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const HyperplaneSectionConfig& config = HyperplaneSectionConfig()
);

/**
 * @brief Sample multiple hyperplane sections
 * 
 * Generates multiple sets of random hyperplanes to get better
 * coverage of the variety, especially for disconnected components.
 * 
 * @param polynomials Original polynomial system
 * @param variables Variable names
 * @param num_samples Number of sample sets
 * @param config Configuration options
 * @return Vector of results from different hyperplane sections
 */
std::vector<HyperplaneSectionResult> sample_variety_points(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    int num_samples,
    const HyperplaneSectionConfig& config = HyperplaneSectionConfig()
);

/**
 * @brief Enhanced solver that handles both zero and positive dimensional systems
 * 
 * Automatically detects dimension and applies hyperplane sections if needed.
 * For zero-dimensional systems, behaves like solve_polynomial_system_enhanced.
 * For positive-dimensional, returns sample points via hyperplane sections.
 */
HyperplaneSectionResult solve_polynomial_system_adaptive(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const HyperplaneSectionConfig& config = HyperplaneSectionConfig()
);

} // namespace julia_rur