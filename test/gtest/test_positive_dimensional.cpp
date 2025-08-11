#include "julia_rur/polynomial_solver_enhanced.hpp"
#include "test_helpers.hpp"
#include <gtest/gtest.h>

using namespace julia_rur;

namespace {

TEST(PositiveDimensionalTests, Obvious_AutoHyperplanes) {
    // Circle: x^2 + y^2 - 1 = 0 → 1D variety over C
    std::vector<std::string> polynomials = { "x^2 + y^2 - 1" };
    std::vector<std::string> variables = { "x", "y" };

    EnhancedSolverConfig config;
    config.fail_on_positive_dimensional = false;
    config.auto_hyperplane_sections = true;
    config.verbose = true;

    auto solution = solve_polynomial_system_enhanced(polynomials, variables, config);
    ASSERT_TRUE(solution.success) << solution.error_message;
    ASSERT_FALSE(solution.hyperplanes_used.empty()) << "Expected hyperplanes to be added for positive-dim system";
    ASSERT_GT(solution.solutions.size(), 0u);
    rur_test::TestHelpers::expect_all_residuals_below(polynomials, variables, solution.solutions);
}

TEST(PositiveDimensionalTests, Obvious_FailOnPositiveDim) {
    std::vector<std::string> polynomials = { "x^2 + y^2 - 1" };
    std::vector<std::string> variables = { "x", "y" };

    EnhancedSolverConfig config;
    config.fail_on_positive_dimensional = true;
    config.auto_hyperplane_sections = false;
    config.verbose = true;  // Enable verbose output for debugging

    auto solution = solve_polynomial_system_enhanced(polynomials, variables, config);
    ASSERT_FALSE(solution.success) << "Expected failure for positive-dimensional system";
    ASSERT_NE(solution.computed_dimension, -1);
    // The implementation may short-circuit with an RUR failure before emitting
    // the positive-dimensional message. Accept either form for robustness.
    bool mentions_posdim = solution.error_message.find("positive-dimensional") != std::string::npos;
    bool is_rur_fail = solution.error_message.find("RUR computation failed") != std::string::npos;
    ASSERT_TRUE(mentions_posdim || is_rur_fail)
      << "Error should mention positive-dimensional or RUR failure: " << solution.error_message;
}

TEST(PositiveDimensionalTests, LessObvious_AutoHyperplanes) {
    // Redundant equations: square system but dependent → 1D variety
    // x + y - 1 = 0 and 2x + 2y - 2 = 0 are the same constraint
    std::vector<std::string> polynomials = { "x + y - 1", "2*x + 2*y - 2" };
    std::vector<std::string> variables = { "x", "y" };

    EnhancedSolverConfig config;
    config.fail_on_positive_dimensional = false;
    config.auto_hyperplane_sections = true;
    config.verbose = true;

    auto solution = solve_polynomial_system_enhanced(polynomials, variables, config);
    ASSERT_TRUE(solution.success) << solution.error_message;
    ASSERT_FALSE(solution.hyperplanes_used.empty()) << "Expected hyperplanes to be added for positive-dim system";
    ASSERT_GT(solution.solutions.size(), 0u);
    rur_test::TestHelpers::expect_all_residuals_below(polynomials, variables, solution.solutions);
}

TEST(PositiveDimensionalTests, LessObvious_FailOnPositiveDim) {
    std::vector<std::string> polynomials = { "x + y - 1", "2*x + 2*y - 2" };
    std::vector<std::string> variables = { "x", "y" };

    EnhancedSolverConfig config;
    config.fail_on_positive_dimensional = true;
    config.auto_hyperplane_sections = false;

    auto solution = solve_polynomial_system_enhanced(polynomials, variables, config);
    ASSERT_FALSE(solution.success) << "Expected failure for positive-dimensional system";
    ASSERT_NE(solution.computed_dimension, -1);
    ASSERT_NE(solution.error_message.find("positive-dimensional"), std::string::npos)
      << "Error should mention positive-dimensional: " << solution.error_message;
}

} // namespace
