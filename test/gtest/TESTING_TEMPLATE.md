# RUR C++ Testing Framework Guide

This guide shows how to add new polynomial system tests to the unified Google Test framework.

## Quick Start

1. **Build with Google Test support:**
   ```bash
   # Install Google Test if not available
   ./vcpkg/vcpkg install gtest
   
   # Configure and build
   mkdir build && cd build
   cmake ..
   make rur_tests
   ```

2. **Run all tests:**
   ```bash
   ./rur_tests
   ```

3. **Run specific test categories:**
   ```bash
   ./rur_tests --gtest_filter="*Linear*"           # Only linear systems
   ./rur_tests --gtest_filter="*Univariate*"       # Only univariate
   ./rur_tests --gtest_filter="BasicOperations.*"  # Only basic ops
   ```

## Test Categories

The framework organizes tests into categories:

- **BasicOperations.***: Monomial and polynomial operations
- **UnivariateSystem.***: Single-variable polynomial systems  
- **LinearSystem.***: Linear equation systems
- **NonlinearSystem.***: Nonlinear polynomial systems
- **NumericalRoots.***: Root-finding algorithms (Eigen/FLINT)
- **Integration.***: End-to-end system tests

## Adding New Tests

### 1. Simple Single Test

Add to the appropriate category file (e.g., `test_nonlinear_systems.cpp`):

```cpp
TEST_F(NonlinearSystem, MyNewSystem) {
    std::vector<std::string> polynomials = {"x^2 + y^2 - 1", "x - y"};
    std::vector<std::string> variables = {"x", "y"};
    
    auto solution = solve_system(polynomials, variables);
    assert_solution_valid(solution, 2);  // Expect 2 solutions
    
    // Additional verification
    auto real_solutions = TestHelpers::extract_real_solutions(solution.solutions);
    EXPECT_EQ(real_solutions.size(), 2);
    
    // Check specific solution properties
    for (const auto& sol : real_solutions) {
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0] + sol[1], 0.0));  // x = y
        EXPECT_TRUE(TestHelpers::approx_equal(sol[0]*sol[0] + sol[1]*sol[1], 1.0));  // On unit circle
    }
}
```

### 2. Parameterized Tests

For testing similar cases with different parameters:

```cpp
struct MyTestCase {
    std::string name;
    std::vector<std::string> polynomials;
    int expected_solutions;
};

class MyParameterizedTest : public ParametrizedPolynomialTest<MyTestCase> {};

TEST_P(MyParameterizedTest, TestVariousCases) {
    MyTestCase test_case = GetParam();
    
    std::vector<std::string> variables = {"x", "y"};
    auto solution = solve_system(test_case.polynomials, variables);
    assert_solution_valid(solution, test_case.expected_solutions);
}

INSTANTIATE_TEST_SUITE_P(
    VariousSystems,
    MyParameterizedTest,
    ::testing::Values(
        MyTestCase{"circle", {"x^2 + y^2 - 1"}, -1},
        MyTestCase{"line", {"x - y"}, -1},
        MyTestCase{"intersection", {"x^2 + y^2 - 1", "x - y"}, 2}
    )
);
```

### 3. Custom Test Fixtures

For tests requiring special setup:

```cpp
class MySpecialTest : public PolynomialSystemTest {
protected:
    void SetUp() override {
        PolynomialSystemTest::SetUp();
        
        // Custom configuration
        config_.timeout_seconds = 60.0;
        config_.use_rational_reconstruction = false;
    }
    
    // Helper methods
    void verify_my_property(const julia_rur::EnhancedNumericalSolution& solution) {
        // Custom verification logic
    }
};

TEST_F(MySpecialTest, ComplexSystem) {
    // Test implementation
}
```

## Helper Functions

The `TestHelpers` class provides utilities:

```cpp
// Tolerance checking
TestHelpers::approx_equal(a, b, tolerance);
TestHelpers::approx_equal(complex1, complex2, tolerance);

// Solution analysis  
TestHelpers::compute_residual(polynomials, variables, solution);
TestHelpers::is_valid_solution(polynomials, variables, solution, tolerance);
TestHelpers::sort_solutions(solutions);
TestHelpers::count_real_solutions(solutions);
TestHelpers::extract_real_solutions(solutions);
```

## Test Assertions

Use these patterns for common checks:

```cpp
// Basic success check
ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;

// Solution count check  
EXPECT_EQ(solution.solutions.size(), expected_count);

// Use the helper for comprehensive validation
assert_solution_valid(solution, expected_count);

// Residual checks
for (size_t i = 0; i < solution.solutions.size(); i++) {
    EXPECT_LT(solution.residuals[i], TestHelpers::DEFAULT_TOLERANCE)
        << "Solution " << i << " has high residual";
}

// Real/complex solution checks
int real_count = TestHelpers::count_real_solutions(solution.solutions);
EXPECT_EQ(real_count, expected_real_count);
```

## Benchmark/Performance Tests

For performance-sensitive tests:

```cpp
TEST_F(Integration, PerformanceTest) {
    auto start = std::chrono::high_resolution_clock::now();
    
    auto solution = solve_system(polynomials, variables);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    assert_solution_valid(solution);
    EXPECT_LT(duration.count(), 1000) << "Test took " << duration.count() << "ms";
}
```

## Error Testing

For testing error conditions:

```cpp
TEST_F(NonlinearSystem, InconsistentSystem) {
    std::vector<std::string> polynomials = {"x^2 - 1", "x^2 - 4"};  // No solutions
    std::vector<std::string> variables = {"x"};
    
    auto solution = solve_system(polynomials, variables);
    
    // Should either fail or find no solutions
    if (solution.success) {
        EXPECT_EQ(solution.solutions.size(), 0);
    } else {
        EXPECT_FALSE(solution.error_message.empty());
    }
}
```

## Adding New Test Files

1. Create new file in `/test/gtest/` directory:
   - Name it `test_[category]_[description].cpp`
   - Include `test_helpers.hpp`
   - Use appropriate test fixture class

2. Add to CMakeLists.txt:
   ```cmake
   add_executable(rur_tests
       # ... existing files ...
       test/gtest/test_my_new_category.cpp
   )
   ```

3. Rebuild and test:
   ```bash
   make rur_tests
   ./rur_tests --gtest_filter="MyNewCategory.*"
   ```

## Best Practices

1. **Naming**: Use descriptive test names that indicate what is being tested
2. **Independence**: Each test should be independent and not rely on other tests
3. **Tolerance**: Use appropriate tolerance levels for numerical comparisons
4. **Documentation**: Add comments explaining complex test cases
5. **Coverage**: Test both success and failure cases
6. **Performance**: Include timing assertions for performance-critical code
7. **Parameterization**: Use parameterized tests for similar test cases

## Running Specific Tests

```bash
# List all available tests
./rur_tests --gtest_list_tests

# Run tests matching pattern
./rur_tests --gtest_filter="*Circle*"

# Run with verbose output
./rur_tests --gtest_filter="BasicOperations.MonomialCreation" --verbose

# Run tests and show detailed failures
./rur_tests --gtest_output=xml:test_results.xml
```

## Debugging Tests

```bash
# Run single test with debugging info
gdb --args ./rur_tests --gtest_filter="NonlinearSystem.CircleLineIntersection"

# Run with AddressSanitizer (if compiled with it)
export ASAN_OPTIONS=abort_on_error=1
./rur_tests
```

This framework makes it easy to add comprehensive tests while maintaining good organization and clear reporting.