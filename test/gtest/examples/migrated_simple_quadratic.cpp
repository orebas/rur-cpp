#include "../test_helpers.hpp"

using namespace rur_test;

/**
 * @brief Example of migrating existing test to modern gtest framework
 * 
 * This demonstrates converting test/test_simple_quadratic.cpp to the new framework.
 * Original test had manual verification - now uses automated assertions.
 */
class MigratedQuadratic : public PolynomialSystemTest {
};

TEST_F(MigratedQuadratic, SimpleQuadraticSystem) {
    // Original: y - 1 = 0, x^2 - 2 = 0
    // Expected: y = 1, x = ±√2
    std::vector<std::string> polynomials = {
        "y - 1",      // y = 1
        "x^2 - 2"     // x^2 = 2
    };
    std::vector<std::string> variables = {"x", "y"};
    
    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);
    
    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);
    
    // Sort solutions by x value
    std::sort(real_solutions.begin(), real_solutions.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
                  return a[0] < b[0];
              });
    
    // Verify solutions: (-√2, 1) and (√2, 1)
    double sqrt2 = std::sqrt(2.0);
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], -sqrt2));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], 1.0));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][0], sqrt2));
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][1], 1.0));
    
    // Verify by substitution (original test's verification logic)
    for (const auto& sol : real_solutions) {
        double x = sol[0], y = sol[1];
        EXPECT_TRUE(TestHelpers::approx_equal(y, 1.0)) << "Error in y - 1";
        EXPECT_TRUE(TestHelpers::approx_equal(x*x, 2.0)) << "Error in x^2 - 2";
    }
}

TEST_F(MigratedQuadratic, CoupledQuadraticSystem) {
    // Original: y^2 - 1 = 0, x^2 - y - 1 = 0
    // Expected: y = ±1, x^2 = y + 1
    // For y = 1: x = ±√2, For y = -1: x = 0
    std::vector<std::string> polynomials = {
        "y^2 - 1",        // y^2 = 1
        "x^2 - y - 1"     // x^2 = y + 1
    };
    std::vector<std::string> variables = {"x", "y"};
    
    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 3);  // Should have 3 real solutions
    
    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 3);
    
    // Sort solutions by y value, then by x value
    std::sort(real_solutions.begin(), real_solutions.end(),
              [](const std::vector<double>& a, const std::vector<double>& b) {
                  if (std::abs(a[1] - b[1]) > 1e-10) return a[1] < b[1];
                  return a[0] < b[0];
              });
    
    // Expected solutions: (0, -1), (-√2, 1), (√2, 1)
    double sqrt2 = std::sqrt(2.0);
    
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][0], 0.0));    // x = 0
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[0][1], -1.0));   // y = -1
    
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][0], -sqrt2)); // x = -√2
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[1][1], 1.0));    // y = 1
    
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[2][0], sqrt2));  // x = √2
    EXPECT_TRUE(TestHelpers::approx_equal(real_solutions[2][1], 1.0));    // y = 1
    
    // Verify system equations for all solutions
    for (const auto& sol : real_solutions) {
        double x = sol[0], y = sol[1];
        EXPECT_TRUE(TestHelpers::approx_equal(y*y, 1.0)) << "Solution (" << x << ", " << y << ") fails y^2 = 1";
        EXPECT_TRUE(TestHelpers::approx_equal(x*x, y + 1.0)) << "Solution (" << x << ", " << y << ") fails x^2 = y + 1";
    }
}

/**
 * @brief Migration Notes:
 * 
 * Improvements over original test:
 * 1. Automatic pass/fail determination instead of manual inspection
 * 2. Structured assertions with descriptive error messages
 * 3. Consistent tolerance handling via TestHelpers
 * 4. Integration with test framework for reporting and filtering
 * 5. Clear separation of test setup, execution, and verification
 * 6. Better error reporting when tests fail
 * 
 * To migrate other tests:
 * 1. Convert main() or test functions to TEST_F methods
 * 2. Replace manual prints with EXPECT/ASSERT statements
 * 3. Use TestHelpers for common operations (residual computation, solution sorting)
 * 4. Structure verification logic as assertions rather than prints
 * 5. Use the PolynomialSystemTest fixture for common setup
 */