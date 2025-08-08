#include "test_helpers.hpp"

using namespace rur_test;

/**
 * @brief Test nonlinear polynomial systems
 */
class NonlinearSystemTests : public PolynomialSystemTest {};

TEST_F(NonlinearSystemTests, CircleLineIntersection) {
    // Test: x^2 + y^2 - 1 = 0, y - x = 0
    // Circle x^2 + y^2 = 1 intersected with line y = x
    // Solutions: (±√2/2, ±√2/2)
    std::vector<std::string> polynomials = { "x^2 + y^2 - 1", "y - x" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);

    // Sort solutions by x-coordinate
    TestHelpers::sort_solutions(solution.solutions);
    real_solutions = TestHelpers::extract_real_solutions(solution.solutions);

    double sqrt_half = std::sqrt(0.5);
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], -sqrt_half));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], -sqrt_half));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][0], sqrt_half));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][1], sqrt_half));
}

TEST_F(NonlinearSystemTests, ParabolaLineIntersection) {
    // Test: y - x^2 = 0, y - x - 2 = 0
    // Parabola y = x^2 intersected with line y = x + 2
    // Solutions: x^2 = x + 2 → x^2 - x - 2 = 0 → x = -1, 2
    std::vector<std::string> polynomials = { "y - x^2", "y - x - 2" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);

    // Sort by x-coordinate
    std::sort(real_solutions.begin(),
              real_solutions.end(),
              [](const std::vector<double> &a, const std::vector<double> &b) { return a[0] < b[0]; });

    // Check solutions: (-1, 1) and (2, 4)
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], -1.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 1.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][0], 2.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][1], 4.0));
}

TEST_F(NonlinearSystemTests, TwoCirclesIntersection) {
    // Test: x^2 + y^2 - 1 = 0, (x-1)^2 + y^2 - 1 = 0
    // Unit circle at origin and unit circle at (1,0)
    // Solutions: (0.5, ±√3/2)
    std::vector<std::string> polynomials = { "x^2 + y^2 - 1", "(x-1)^2 + y^2 - 1" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);

    // Both solutions should have x = 0.5
    for (const auto &sol : real_solutions) {
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0], 0.5));
        EXPECT_TRUE(TestHelpers::approx_equal(std::abs(sol[1]), std::sqrt(3.0) / 2.0));
    }
}

TEST_F(NonlinearSystemTests, QuadraticSystem) {
    // Test: x^2 + y^2 = 5, x*y = 2
    // Solutions can be found by substitution
    std::vector<std::string> polynomials = { "x^2 + y^2 - 5", "x*y - 2" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 4);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 4);

    // Verify all solutions satisfy both equations
    for (const auto &sol : real_solutions) {
        double x = sol[0], y = sol[1];
        EXPECT_TRUE(TestHelpers::approx_equal(x * x + y * y, 5.0));
        EXPECT_TRUE(TestHelpers::approx_equal(x * y, 2.0));
    }
}

TEST_F(NonlinearSystemTests, CubicSystem) {
    // Test: x^3 - y = 0, x^2 + y^2 - 2 = 0
    // Cubic curve y = x^3 intersected with circle x^2 + y^2 = 2
    std::vector<std::string> polynomials = { "x^3 - y", "x^2 + y^2 - 2" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    // This system typically has multiple solutions (some may be complex)
    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;
    EXPECT_GE(solution.solutions.size(), 1);
}

TEST_F(NonlinearSystemTests, SymmetricSystem) {
    // Elementary symmetric polynomials: x + y = 3, x*y = 2
    // Solutions: (1, 2) and (2, 1)
    std::vector<std::string> polynomials = { "x + y - 3", "x*y - 2" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);

    // Sort solutions
    std::sort(real_solutions.begin(),
              real_solutions.end(),
              [](const std::vector<double> &a, const std::vector<double> &b) { return a[0] < b[0]; });

    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 1.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 2.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][0], 2.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][1], 1.0));
}

TEST_F(NonlinearSystemTests, ThreeVariableSystem) {
    // Test: x^2 - 1, y^2 - 2, z^2 - 3
    // This system has 8 solutions (2^3 combinations)
    std::vector<std::string> polynomials = { "x^2 - 1", "y^2 - 2", "z^2 - 3" };
    std::vector<std::string> variables = { "x", "y", "z" };

    config_.verbose = true;
    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 8); // 2^3 = 8 solutions

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 8);

    // Each solution should be a combination of (±1, ±√2, ±√3)
    double sqrt2 = std::sqrt(2.0);
    double sqrt3 = std::sqrt(3.0);
    
    for (size_t i = 0; i < real_solutions.size(); ++i) {
        const auto &sol = real_solutions[i];
        
        std::cout << "Solution " << i << ": ["
                  << sol[0] << ", " << sol[1] << ", " << sol[2] << "]"
                  << std::endl;

        // Check that x = ±1
        EXPECT_TRUE(std::abs(std::abs(sol[0]) - 1.0) < 1e-6);
        // Check that y = ±√2
        EXPECT_TRUE(std::abs(std::abs(sol[1]) - sqrt2) < 1e-6);
        // Check that z = ±√3
        EXPECT_TRUE(std::abs(std::abs(sol[2]) - sqrt3) < 1e-6);
    }
}

TEST_F(NonlinearSystemTests, BilinearSystem) {
    // Test: x*y - 1 = 0, x + y - 3 = 0
    // From first equation: y = 1/x
    // Substituting into second: x + 1/x = 3 → x^2 - 3x + 1 = 0
    // Solutions: x = (3 ± √5)/2
    std::vector<std::string> polynomials = { "x*y - 1", "x + y - 3" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);

    // Verify that x*y = 1 and x + y = 3 for both solutions
    for (const auto &sol : real_solutions) {
        double x = sol[0], y = sol[1];
        EXPECT_TRUE(TestHelpers::approx_equal(x * y, 1.0));
        EXPECT_TRUE(TestHelpers::approx_equal(x + y, 3.0));
    }
}

TEST_F(NonlinearSystemTests, SquareSystem) {
    // Test: x^2 - y = 0, y^2 - x = 0
    // Solutions satisfy x^4 - x = 0 → x(x^3 - 1) = 0
    // So x ∈ {0, 1, ω, ω^2} where ω = cube root of unity
    std::vector<std::string> polynomials = { "x^2 - y", "y^2 - x" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;
    EXPECT_GE(solution.solutions.size(), 2); // At least (0,0) and (1,1)

    // Check that real solutions include (0,0) and (1,1)
    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    bool found_origin = false, found_unit = false;

    for (const auto &sol : real_solutions) {
        if (TestHelpers::approx_equal(sol[0], 0.0) && TestHelpers::approx_equal(sol[1], 0.0)) { found_origin = true; }
        if (TestHelpers::approx_equal(sol[0], 1.0) && TestHelpers::approx_equal(sol[1], 1.0)) { found_unit = true; }
    }

    EXPECT_TRUE(found_origin) << "Solution (0,0) not found";
    EXPECT_TRUE(found_unit) << "Solution (1,1) not found";
}

TEST_F(NonlinearSystemTests, HigherOrderSystem) {
    // Test: x^3 + y^3 = 2, x + y = 2
    // From second equation: y = 2 - x
    // Substituting: x^3 + (2-x)^3 = 2
    // Expanding and simplifying leads to solutions
    std::vector<std::string> polynomials = { "x^3 + y^3 - 2", "x + y - 2" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;
    EXPECT_GE(solution.solutions.size(), 1);

    // Verify all solutions satisfy both equations
    for (size_t i = 0; i < solution.solutions.size(); i++) {
        double residual = TestHelpers::compute_residual(polynomials, variables, solution.solutions[i]);
        EXPECT_LT(residual, TestHelpers::LOOSE_TOLERANCE) << "Solution " << i << " has high residual: " << residual;
    }

    // Check that real solutions satisfy the constraint x + y = 2
    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    for (const auto &sol : real_solutions) {
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0] + sol[1], 2.0))
          << "Solution (" << sol[0] << ", " << sol[1] << ") doesn't satisfy x + y = 2";
    }
}