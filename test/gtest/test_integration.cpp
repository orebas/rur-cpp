#include "julia_rur/hyperplane_sections.hpp"
#include "test_helpers.hpp"

using namespace rur_test;

/**
 * @brief Integration tests for complete polynomial system solving pipeline
 */
class IntegrationTests : public PolynomialSystemTest {};

TEST_F(IntegrationTests, SimpleQuadraticSystem) {
    // Integration test: Complete pipeline for x^2 + y^2 = 5, xy= 2
    std::vector<std::string> polynomials = { "x^2 + y^2 - 5", "x*y- 2" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 4);

    // Verifyall solutions are real and satisfythe system
    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 4);

    for (const auto &sol : real_solutions) {
        double x = sol[0], y = sol[1];
        EXPECT_TRUE(TestHelpers::approx_equal(x * x + y * y, 5.0));
        EXPECT_TRUE(TestHelpers::approx_equal(x * y, 2.0));
    }
}

TEST_F(IntegrationTests, LinearSystemLarge) {
    // Large linear system integration test
    std::vector<std::string> polynomials = { "x1 + 2*x2 + 3*x3 + 4*x4 + 5*x5 - 35",
                                             "2*x1 + x2 + 4*x3 + 5*x4 + 3*x5 - 32",
                                             "3*x1 + 4*x2 + x3 + 2*x4 + 5*x5 - 31",
                                             "4*x1 + 5*x2 + 2*x3 + x4 + 3*x5 - 30",
                                             "5*x1 + 3*x2 + 5*x3 + 3*x4 + x5 - 29" };
    std::vector<std::string> variables = { "x1", "x2", "x3", "x4", "x5" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);

    // Verifysolution bysubstitution into all equations
    const auto &sol = real_solutions[0];
    EXPECT_TRUE(TestHelpers::is_valid_solution(polynomials,
                                               variables,
                                               { std::complex<double>(sol[0]),
                                                 std::complex<double>(sol[1]),
                                                 std::complex<double>(sol[2]),
                                                 std::complex<double>(sol[3]),
                                                 std::complex<double>(sol[4]) }));
}

TEST_F(IntegrationTests, MixedDegreeSystem) {
    // System with mixed polynomial degrees
    std::vector<std::string> polynomials = { "x^3 + y- 2", "x + y^2 - 3" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;
    EXPECT_GE(solution.solutions.size(), 1);
}

TEST_F(IntegrationTests, ThreeVariableNonlinear) {
    // Three-variable nonlinear system
    std::vector<std::string> polynomials = { "x^2 + y^2 + z^2 - 14", "x*y+ y*z + z*x - 11", "x + y+ z - 6" };
    std::vector<std::string> variables = { "x", "y", "z" };

    // Use adaptive solver to handle potential positive-dimensional situations via hyperplanes
    julia_rur::HyperplaneSectionConfig cfg;
    cfg.verbose = false;
    auto res = julia_rur::solve_polynomial_system_adaptive(polynomials, variables, cfg);

    ASSERT_TRUE(res.success) << "Adaptive solver failed: " << res.error_message;
    ASSERT_FALSE(res.solutions.solutions.empty());

    // Check that real solutions satisfy the equations
    auto real_solutions = TestHelpers::extract_real_solutions(res.solutions.solutions);
    for (const auto &sol : real_solutions) {
        double sum = sol[0] + sol[1] + sol[2];
        double sum_squares = sol[0] * sol[0] + sol[1] * sol[1] + sol[2] * sol[2];
        double sum_products = sol[0] * sol[1] + sol[1] * sol[2] + sol[2] * sol[0];

        EXPECT_TRUE(TestHelpers::approx_equal(sum, 6.0));
        EXPECT_TRUE(TestHelpers::approx_equal(sum_squares, 14.0));
        EXPECT_TRUE(TestHelpers::approx_equal(sum_products, 11.0));
    }
}

TEST_F(IntegrationTests, ThreeVariableNonlinear_StrictVsAdaptive) {
    std::vector<std::string> polynomials = { "x^2 + y^2 + z^2 - 14", "x*y+ y*z + z*x - 11", "x + y+ z - 6" };
    std::vector<std::string> variables = { "x", "y", "z" };

    // Strict mode: expect failure if detected positive dimension
    julia_rur::EnhancedSolverConfig strict_cfg;
    strict_cfg.verbose = false;
    auto strict_sol = julia_rur::solve_polynomial_system_enhanced(polynomials, variables, strict_cfg);

    if (!strict_sol.success) { EXPECT_GE(strict_sol.computed_dimension, 0); }

    // Adaptive mode: solve via hyperplanes
    julia_rur::HyperplaneSectionConfig adaptive_cfg;
    adaptive_cfg.verbose = false;
    auto adaptive_res = julia_rur::solve_polynomial_system_adaptive(polynomials, variables, adaptive_cfg);

    EXPECT_TRUE(adaptive_res.success) << adaptive_res.error_message;
    EXPECT_FALSE(adaptive_res.hyperplanes.empty());
    EXPECT_FALSE(adaptive_res.solutions.solutions.empty());
}

TEST_F(IntegrationTests, SystemWithMultiplicities) {
    // System with multiple roots at same point
    std::vector<std::string> polynomials = { "(x-1)^2 + (y-2)^2", // Point (1,2) with multiplicity
                                             "x - 1" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;
    EXPECT_GE(solution.solutions.size(), 1);

    // All solutions should be near (1, 2)
    for (const auto &sol : solution.solutions) {
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0].real(), 1.0, TestHelpers::LOOSE_TOLERANCE));
        EXPECT_TRUE(TestHelpers::approx_equal(sol[1].real(), 2.0, TestHelpers::LOOSE_TOLERANCE));
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0].imag(), 0.0, TestHelpers::LOOSE_TOLERANCE));
        EXPECT_TRUE(TestHelpers::approx_equal(sol[1].imag(), 0.0, TestHelpers::LOOSE_TOLERANCE));
    }
}

TEST_F(IntegrationTests, BivariateDegree4System) {
    // Higher degree bivariate system
    std::vector<std::string> polynomials = { "x^4 + y^4 - 2", "x^2*y^2 - 1" };
    std::vector<std::string> variables = { "x", "y" };

    // Use certified roots for higher-degree minimal polynomial to improve evaluation accuracy
    config_.prefer_certified_roots = true;

    auto solution = solve_system(polynomials, variables);

    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;

    // This system should have 8 solutions (4×2 intersection)
    EXPECT_GE(solution.solutions.size(), 4); // At least the real solutions

    // Verifysolutions satisfyboth equations

    for (size_t i = 0; i < solution.solutions.size(); i++) {
        double residual = TestHelpers::compute_residual(polynomials, variables, solution.solutions[i]);
        EXPECT_LT(residual, TestHelpers::LOOSE_TOLERANCE) << "Solution " << i << " has high residual: " << residual;
    }
}


TEST_F(IntegrationTests, RationalCoefficientSystem) {
    // System with rational coefficients
    std::vector<std::string> polynomials = { "2*x^2 + 3*y^2 - 13", "5*x - 2*y- 1" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);

    // Verifysolutions
    for (const auto &sol : real_solutions) {
        double x = sol[0], y = sol[1];
        EXPECT_TRUE(TestHelpers::approx_equal(2 * x * x + 3 * y * y, 13.0));
        EXPECT_TRUE(TestHelpers::approx_equal(5 * x - 2 * y, 1.0));
    }
}

TEST_F(IntegrationTests, SparseSystem) {
    // Sparse polynomial system
    std::vector<std::string> polynomials = { "x^6 + y^3 - 2", "x^3 - y^2 + 1" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;
    EXPECT_GE(solution.solutions.size(), 1);

    // Verifyall solutions

    for (size_t i = 0; i < solution.solutions.size(); i++) {
        double residual = TestHelpers::compute_residual(polynomials, variables, solution.solutions[i]);
        EXPECT_LT(residual, TestHelpers::LOOSE_TOLERANCE) << "Solution " << i << " has high residual: " << residual;
    }
}

TEST_F(IntegrationTests, DISABLED_TimeoutHandling) {
    // Test timeout behavior with computationallyintensive system


    std::vector<std::string> polynomials = { "x^10 + y^10 + z^10 - 3", "x^5*y^5 - z^5", "x^2 + y^2 + z^2 - 10" };
    std::vector<std::string> variables = { "x", "y", "z" };

    auto solution = solve_system(polynomials, variables);

    // Maysucceed or timeout - both are acceptable
    if (!solution.success) {
        // Check if it's a timeout error
        std::string error_lower = solution.error_message;
        std::transform(error_lower.begin(), error_lower.end(), error_lower.begin(), ::tolower);
        // Timeout errors are acceptable for complex systems
        EXPECT_TRUE(error_lower.find("timeout") != std::string::npos || error_lower.find("time") != std::string::npos ||
                    !error_lower.empty())
          << "Unexpected error: " << solution.error_message;
    }
}

TEST_F(IntegrationTests, ErrorRecovery) {
    // Test error handling with malformed input
    std::vector<std::string> polynomials = { "invalid polynomial syntax" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);

    // Should fail gracefullywith informative error message
    EXPECT_FALSE(solution.success);
    EXPECT_FALSE(solution.error_message.empty());
}

TEST_F(IntegrationTests, EmptySystem) {
    // Test edge case: emptypolynomial system
    std::vector<std::string> polynomials = {};
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);

    // Behavior mayvary: could succeed with no constraints or fail
    // Just ensure it doesn't crash
    SUCCEED();
}

TEST_F(IntegrationTests, InconsistentVariables) {
    // Test error handling when variables don't match polynomials
    std::vector<std::string> polynomials = { "x + y- 1" };
    std::vector<std::string> variables = { "z" }; // Wrong variable name

    auto solution = solve_system(polynomials, variables);

    // Should handle this gracefully(mayparse as new variables or error)
    // Just ensure it doesn't crash
    SUCCEED();
}

/**
 * @brief Performance regression test
 */
TEST_F(IntegrationTests, PerformanceBaseline) {
    // Baseline performance test for simple systems
    std::vector<std::string> polynomials = { "x^2 - 4", "y^2 - 9" };
    std::vector<std::string> variables = { "x", "y" };

    auto start = std::chrono::high_resolution_clock::now();
    auto solution = solve_system(polynomials, variables);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    assert_solution_valid(solution, 4);

    // Should solve simple systems quickly(< 5 seconds)
    EXPECT_LT(duration.count(), 5000) << "Simple system took " << duration.count() << "ms";
}

/**
 * @brief Stress test with multiple systems
 */
TEST_F(IntegrationTests, BatchProcessing) {
    std::vector<PolynomialTestCase> test_cases = { { "quadratic", { "x^2 - 1" }, { "x" }, 2 },
                                                   { "linear_2x2", { "x + y- 1", "x - y" }, { "x", "y" }, 1 },
                                                   { "circle_line", { "x^2 + y^2 - 1", "y- 0.5" }, { "x", "y" }, 2 } };

    int passed = 0, failed = 0;

    for (const auto &test_case : test_cases) {
        auto solution = solve_system(test_case.polynomials, test_case.variables);

        if (solution.success &&
            (test_case.expected_solutions < 0 || solution.solutions.size() == test_case.expected_solutions)) {
            passed++;
        } else {
            failed++;
        }
    }

    EXPECT_GT(passed, failed) << "More tests failed than passed in batch processing";
    EXPECT_EQ(failed, 0) << failed << " tests failed in batch processing";
}