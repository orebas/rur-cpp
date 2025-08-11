#pragma once

#include "numerical_roots_flint.hpp"
#include "polynomial_solver.hpp"
#include <optional>


namespace julia_rur {

/**
 * @brief Root finding method selection
 */
enum class RootFindingMethod {
    AUTO,  // Automatically choose best method
    EIGEN, // Use Eigen's polynomial solver
    FLINT, // Use FLINT's arbitrary precision solver
    HYBRID // Use FLINT for certification, Eigen for speed
};

/**
 * @brief Enhanced configuration for polynomial solving
 */
struct EnhancedSolverConfig : public RURConfig {
    RootFindingMethod root_method = RootFindingMethod::AUTO;
    FlintRootConfig flint_config;

    // Options for automatic method selection
    int degree_threshold_for_flint = 20; // Use FLINT for degree > this
    bool prefer_certified_roots = false; // Prefer FLINT when accuracy matters

    // Root clustering and separation
    double clustering_threshold = 1e-10;  // Distance to consider roots clustered
    bool separate_conjugate_pairs = true; // Group complex conjugate pairs

    // Solution filtering
    bool real_roots_only = false;       // Only return real solutions
    double imaginary_tolerance = 1e-12; // Tolerance for considering root real
    
    // Positive-dimensional handling
    bool fail_on_positive_dimensional = false;  // If true, return error instead of using hyperplanes
    bool auto_hyperplane_sections = true;       // If true, automatically add hyperplanes for positive-dim

    // Internal: allow callers (e.g., hyperplane slicing) to bypass the dimension pre-check
    // and go straight to RUR for augmented systems known (or expected) to be zero-dimensional.
    bool skip_dimension_precheck = false;

    // Verbose diagnostics for evaluation steps
    bool verbose = false;

    // Performance knobs for precheck
    bool enable_dimension_precheck = true; // Run GB-based precheck before RUR
    int max_precheck_variables = 8;        // Skip precheck when variable count exceeds this
};

/**
 * @brief Enhanced solution with additional metadata
 */
struct EnhancedPolynomialSolution : public PolynomialSystemSolution {
    // Root certification (if FLINT was used)
    std::vector<double> root_error_bounds;
    bool roots_certified = false;

    // Root clustering information
    std::vector<int> root_multiplicities;          // Estimated multiplicities
    std::vector<std::vector<int>> conjugate_pairs; // Indices of conjugate pairs

    // Method used
    RootFindingMethod method_used;

    // Performance metrics
    double root_finding_time_ms = 0.0;
    double total_time_ms = 0.0;

    // Dimension and augmentation metadata
    int computed_dimension = -1;               // -1 unknown, >=0 if analyzed
    std::vector<std::string> hyperplanes_used; // Non-empty when adaptive/hyperplane mode is used
};

/**
 * @brief Enhanced polynomial system solver with multiple backends
 *
 * Provides high-accuracy root finding using either Eigen or FLINT,
 * with automatic method selection based on problem characteristics.
 */
EnhancedPolynomialSolution
solve_polynomial_system_enhanced(const std::vector<std::string> &polynomials,
                                 const std::vector<std::string> &variables,
                                 const EnhancedSolverConfig &config = EnhancedSolverConfig());

/**
 * @brief Find roots using the best available method
 *
 * Automatically selects between Eigen and FLINT based on:
 * - Polynomial degree
 * - Coefficient size
 * - Accuracy requirements
 */
std::vector<std::complex<double>>
find_roots_auto(const std::vector<mpq_class> &polynomial_coeffs,
                const EnhancedSolverConfig &config,
                std::optional<std::vector<double>> &error_bounds);

/**
 * @brief Analyze polynomial to determine best root finding method
 */
RootFindingMethod
analyze_polynomial_for_method(const std::vector<mpq_class> &polynomial_coeffs, const EnhancedSolverConfig &config);

/**
 * @brief Post-process roots to identify special structures
 *
 * - Identifies real roots
 * - Finds conjugate pairs
 * - Estimates multiplicities
 * - Clusters nearby roots
 */
void
analyze_root_structure(const std::vector<std::complex<double>> &roots,
                       EnhancedPolynomialSolution &solution,
                       const EnhancedSolverConfig &config);

/**
 * @brief Pretty-print enhanced solution with all metadata
 */
void
print_enhanced_solution(const EnhancedPolynomialSolution &solution, std::ostream &out = std::cout);

} // namespace julia_rur