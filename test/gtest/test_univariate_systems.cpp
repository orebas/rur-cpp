#include "test_helpers.hpp"

using namespace rur_test;

/**
 * @brief Test univariate polynomial systems
 */
class UnivariateSystemTests : public PolynomialSystemTest {};

TEST_F(UnivariateSystemTests, SimpleLinearPolynomials) {
    // Test x - 4 = 0, solution: x = 4
    {
        std::vector<std::string> polynomials = { "x - 4" };
        std::vector<std::string> variables = { "x" };
        
        auto solution = solve_system(polynomials, variables);
        assert_solution_valid(solution, 1);
        
        auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
        EXPECT_EQ(real_solutions.size(), 1);
        EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 4.0));
    }
    
    // Test 5*x - 4 = 0, solution: x = 4/5 = 0.8
    {
        std::vector<std::string> polynomials = { "5*x - 4" };
        std::vector<std::string> variables = { "x" };
        
        auto solution = solve_system(polynomials, variables);
        assert_solution_valid(solution, 1);
        
        auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
        EXPECT_EQ(real_solutions.size(), 1);
        EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 0.8));
    }
    
    // Test x - 0 = 0 (i.e., x = 0), solution: x = 0
    {
        std::vector<std::string> polynomials = { "x" };
        std::vector<std::string> variables = { "x" };
        
        auto solution = solve_system(polynomials, variables);
        assert_solution_valid(solution, 1);
        
        auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
        EXPECT_EQ(real_solutions.size(), 1);
        EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 0.0));
    }
    
    // Test 5*x - 0 = 0 (i.e., 5*x = 0), solution: x = 0
    {
        std::vector<std::string> polynomials = { "5*x" };
        std::vector<std::string> variables = { "x" };
        
        auto solution = solve_system(polynomials, variables);
        assert_solution_valid(solution, 1);
        
        auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
        EXPECT_EQ(real_solutions.size(), 1);
        EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 0.0));
    }
}

TEST_F(UnivariateSystemTests, QuadraticRoots) {
    // Test x^2 - 4 = 0, solutions: x = ±2
    std::vector<std::string> polynomials = { "x^2 - 4" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    // Extract solutions and check they are ±2
    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);

    // Sort solutions for consistent comparison
    std::sort(real_solutions.begin(),
              real_solutions.end(),
              [](const std::vector<double> &a, const std::vector<double> &b) { return a[0] < b[0]; });

    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], -2.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][0], 2.0));
}

TEST_F(UnivariateSystemTests, CubicRoots) {
    // Test x^3 - 8 = 0, solutions: x = 2, 2*ω, 2*ω^2 where ω = exp(2πi/3)
    std::vector<std::string> polynomials = { "x^3 - 8" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 3);

    // Should have exactly one real root (x = 2)
    int real_count = TestHelpers::count_real_solutions(solution.solutions);
    EXPECT_EQ(real_count, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 2.0));
}

TEST_F(UnivariateSystemTests, RootsOfUnity) {
    // Test x^4 - 1 = 0, solutions: ±1, ±i
    std::vector<std::string> polynomials = { "x^4 - 1" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 4);

    // Should have exactly 2 real roots (±1)
    int real_count = TestHelpers::count_real_solutions(solution.solutions);
    EXPECT_EQ(real_count, 2);
}

TEST_F(UnivariateSystemTests, HighDegreePolynomial) {
    // Test x^10 - 1 = 0 (10th roots of unity)
    std::vector<std::string> polynomials = { "x^10 - 1" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 10);

    // Should have exactly 2 real roots (±1)
    int real_count = TestHelpers::count_real_solutions(solution.solutions);
    EXPECT_EQ(real_count, 2);

    // Check that all solutions have magnitude 1 (lie on unit circle)
    for (const auto &sol : solution.solutions) {
        double magnitude = std::abs(sol[0]);
        EXPECT_TRUE(TestHelpers::approx_equal(magnitude, 1.0, TestHelpers::LOOSE_TOLERANCE))
          << "Solution " << sol[0] << " has magnitude " << magnitude;
    }
}

TEST_F(UnivariateSystemTests, RepeatedRoots) {
    // Test (x-1)^3 = 0, solution: x = 1 with multiplicity 3
    std::vector<std::string> polynomials = { "(x-1)^3" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);

    // The solver should find the root (may report it once or with multiplicity)
    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;
    EXPECT_GE(solution.solutions.size(), 1);

    // All solutions should be x = 1
    for (const auto &sol : solution.solutions) {
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0].real(), 1.0));
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0].imag(), 0.0));
    }
}

TEST_F(UnivariateSystemTests, SparsePolynomial) {
    // Test x^8 + x^4 + x^2 + 1 = 0
    std::vector<std::string> polynomials = { "x^8 + x^4 + x^2 + 1" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);

    // This polynomial has 8 complex roots
    ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;
    EXPECT_EQ(solution.solutions.size(), 8);
}

/**
 * @brief Parametrized tests for various polynomial degrees
 */

// Custom name generator for degree parameter
struct DegreeNameGenerator {
    std::string operator()(const ::testing::TestParamInfo<int>& info) const {
        return "Degree" + std::to_string(info.param);
    }
};

class UnivariateParametrized : public ParametrizedPolynomialTest<int> {};

TEST_P(UnivariateParametrized, SimpleRoots) {
    int degree = GetParam();

    // Test x^n - 1 = 0 (nth roots of unity)
    std::vector<std::string> polynomials = { "x^" + std::to_string(degree) + " - 1" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, degree);

    // All roots should lie on the unit circle
    for (const auto &sol : solution.solutions) {
        double magnitude = std::abs(sol[0]);
        EXPECT_TRUE(TestHelpers::approx_equal(magnitude, 1.0, TestHelpers::LOOSE_TOLERANCE))
          << "Root " << sol[0] << " not on unit circle";
    }
}

INSTANTIATE_TEST_SUITE_P(VariousDegrees, UnivariateParametrized, 
                         ::testing::Values(2, 3, 4, 5, 6, 8, 12),
                         DegreeNameGenerator());

TEST_F(UnivariateSystemTests, PolynomialWithRationalCoefficients) {
    // Test 2*x^2 - 3*x + 1 = 0, solutions: x = 1/2, 1
    std::vector<std::string> polynomials = { "2*x^2 - 3*x + 1" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);

    // Sort solutions
    std::sort(real_solutions.begin(),
              real_solutions.end(),
              [](const std::vector<double> &a, const std::vector<double> &b) { return a[0] < b[0]; });

    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 0.5));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][0], 1.0));
}

TEST_F(UnivariateSystemTests, ComplexCoefficients) {
    // Test x^2 + 1 = 0, solutions: ±i
    std::vector<std::string> polynomials = { "x^2 + 1" };
    std::vector<std::string> variables = { "x" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);

    // Should have no real roots
    int real_count = TestHelpers::count_real_solutions(solution.solutions);
    EXPECT_EQ(real_count, 0);

    // All roots should have real part = 0, imaginary part = ±1
    for (const auto &sol : solution.solutions) {
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0].real(), 0.0));
        EXPECT_TRUE(TestHelpers::approx_equal(std::abs(sol[0].imag()), 1.0));
    }
}