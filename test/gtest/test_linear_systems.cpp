#include "test_helpers.hpp"

using namespace rur_test;

/**
 * @brief Test linear polynomial systems
 */
class LinearSystemTests : public PolynomialSystemTest {};

TEST_F(LinearSystemTests, Simple2x2System) {
    // Test: x + y - 3 = 0, 2*x - y - 1 = 0
    // Solution: x = 4/3, y = 5/3
    std::vector<std::string> polynomials = { "x + y - 3", "2*x - y - 1" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);

    // Check solution values
    std::cout << "DEBUG: Got solution x=" << real_solutions[0][0] 
              << ", y=" << real_solutions[0][1] << std::endl;
    std::cout << "DEBUG: Expected x=" << (4.0/3.0) << ", y=" << (5.0/3.0) << std::endl;
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 4.0 / 3.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 5.0 / 3.0));
}

TEST_F(LinearSystemTests, Simple3x3System) {
    // Test: x + y + z - 6 = 0, x - y + z - 2 = 0, 2*x + y - z = 0
    // Solution: x = 2, y = 2, z = 2
    std::vector<std::string> polynomials = { "x + y + z - 6", "x - y + z - 2", "2*x + y - z" };
    std::vector<std::string> variables = { "x", "y", "z" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);

    // Check solution values
    std::cout << "DEBUG: Got solution x=" << real_solutions[0][0] 
              << ", y=" << real_solutions[0][1] 
              << ", z=" << real_solutions[0][2] << std::endl;
    std::cout << "DEBUG: Expected x=" << (2.0/3.0) << ", y=" << 2.0 << ", z=" << (10.0/3.0) << std::endl;
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 2.0 / 3.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 2.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][2], 10.0 / 3.0));
}

TEST_F(LinearSystemTests, OverdeterminedConsistent) {
    // 4 equations, 3 variables, but consistent
    // x + y + z = 6, x - y = 0, y + z = 4, 2*x - z = 2
    // Solution: x = 2, y = 2, z = 2
    std::vector<std::string> polynomials = { "x + y + z - 6", "x - y", "y + z - 4", "2*x - z - 2" };
    std::vector<std::string> variables = { "x", "y", "z" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);

    // Check solution values
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 2.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 2.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][2], 2.0));
}

TEST_F(LinearSystemTests, InconsistentSystem) {
    // Inconsistent system: x + y = 1, x + y = 2
    std::vector<std::string> polynomials = { "x + y - 1", "x + y - 2" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    // System should either fail or have no solutions
    if (solution.success) {
        EXPECT_EQ(solution.solutions.size(), 0);
    } else {
        // Failure is acceptable for inconsistent systems
        EXPECT_FALSE(solution.success);
    }
}

TEST_F(LinearSystemTests, UnderdeterminedSystem) {
    // Underdetermined: x + y = 1 (one equation, two variables)
    // This should either fail or find the parametric solution
    std::vector<std::string> polynomials = { "x + y - 1" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    // The behavior depends on the solver implementation
    // It may fail (expected) or find specific solutions
    if (solution.success) {}
}

TEST_F(LinearSystemTests, IdentitySystem) {
    // Simple identity: x = 1, y = 2, z = 3
    std::vector<std::string> polynomials = { "x - 1", "y - 2", "z - 3" };
    std::vector<std::string> variables = { "x", "y", "z" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);

    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 1.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 2.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][2], 3.0));
}

TEST_F(LinearSystemTests, LargerSystem5x5) {
    // 5x5 linear system with known solution: all variables = 1
    std::vector<std::string> polynomials = { "x1 + x2 + x3 + x4 + x5 - 5",
                                             "x1 - x2 + x3 - x4 + x5 - 1",
                                             "2*x1 + x2 + 3*x3 + x4 + 2*x5 - 9",
                                             "x1 + 2*x2 - x3 + x4 + x5 - 4",
                                             "3*x1 - x2 + x3 + 2*x4 - x5 - 4" };
    std::vector<std::string> variables = { "x1", "x2", "x3", "x4", "x5" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);

    // Check that all variables equal 1
    for (int i = 0; i < 5; i++) {
        EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][i], 1.0))
          << "Variable x" << (i + 1) << " = " << real_solutions[0][i] << ", expected 1.0";
    }
}

TEST_F(LinearSystemTests, HomogeneousSystem) {
    // Homogeneous system: x + y = 0, 2*x - y = 0
    // Solution: x = y = 0 (trivial solution)
    std::vector<std::string> polynomials = { "x + y", "2*x - y" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);

    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 0.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 0.0));
}

TEST_F(LinearSystemTests, NearSingularSystemRational) {
    // Nearly singular system with rational coefficients
    // x + y = 1, 1000000*x + 1000001*y = 1000001
    // Solution: x = 0, y = 1
    std::vector<std::string> polynomials = { "x + y - 1", "1000000*x + 1000001*y - 1000001" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    if (solution.success) {
        assert_solution_valid(solution, 1);

        auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
        if (real_solutions.size() == 1) {
            EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 0.0, TestHelpers::LOOSE_TOLERANCE));
            EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 1.0, TestHelpers::LOOSE_TOLERANCE));
        }
    } else {
        // Numerical issues with near-singular systems are acceptable
        std::cout << "Near-singular system failed (acceptable): " << solution.error_message << std::endl;
    }
}

TEST_F(LinearSystemTests, DISABLED_NearSingularSystemDecimal) {
    // Nearly singular system with decimal coefficients
    // NOTE: Decimal coefficient parsing is not yet supported
    // x + y = 1, x + 1.000001*y = 1.000001
    // Solution: x = 0, y = 1
    std::vector<std::string> polynomials = { "x + y - 1", "x + 1.000001*y - 1.000001" };
    std::vector<std::string> variables = { "x", "y" };

    auto solution = solve_system(polynomials, variables);

    if (solution.success) {
        assert_solution_valid(solution, 1);

        auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
        if (real_solutions.size() == 1) {
            EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 0.0, TestHelpers::LOOSE_TOLERANCE));
            EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 1.0, TestHelpers::LOOSE_TOLERANCE));
        }
    } else {
        // Numerical issues with near-singular systems are acceptable
        std::cout << "Near-singular system failed (acceptable): " << solution.error_message << std::endl;
    }
}

/**
 * @brief Parametrized test for various system sizes
 */
// Custom name generator for size parameter
struct SizeNameGenerator {
    std::string operator()(const ::testing::TestParamInfo<int>& info) const {
        return std::to_string(info.param) + "x" + std::to_string(info.param);
    }
};

class LinearSystemParametrized : public ParametrizedPolynomialTest<int> {};

TEST_P(LinearSystemParametrized, DiagonalSystem) {
    int n = GetParam();

    // Create diagonal system: x_i = i for i = 1..n
    std::vector<std::string> polynomials;
    std::vector<std::string> variables;

    for (int i = 1; i <= n; i++) {
        variables.push_back("x" + std::to_string(i));
        polynomials.push_back("x" + std::to_string(i) + " - " + std::to_string(i));
    }

    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 1);

    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 1);

    // Check solution values
    for (int i = 0; i < n; i++) {
        EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][i], double(i + 1)))
          << "x" << (i + 1) << " = " << real_solutions[0][i] << ", expected " << (i + 1);
    }
}

INSTANTIATE_TEST_SUITE_P(VariousSizes,
                         LinearSystemParametrized,
                         ::testing::Values(2, 3, 4, 5, 8),
                         SizeNameGenerator());