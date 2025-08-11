#include <gtest/gtest.h>
#include "julia_rur/polynomial_solver_enhanced.hpp"
#include "julia_rur/hyperplane_sections.hpp"
#include "test_helpers.hpp"
#include <iostream>

using namespace julia_rur;

class UnderdeterminedSystemTests : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up configuration for hyperplane slicing
        config.auto_hyperplane_sections = true;
        config.verbose = true;
        config.force_slicing = false;
        config.max_dimension = 10;
        config.num_sample_sets = 2;
    }

    HyperplaneSectionConfig config;
};

// Test 1: Trivial underdetermined system - single linear equation in 3 variables
// Dimension = 2, needs 2 hyperplanes
TEST_F(UnderdeterminedSystemTests, SingleLinearEquationIn3D) {
    std::cout << "\n=== Testing: Single linear equation in 3 variables ===" << std::endl;
    
    std::vector<std::string> polynomials = {
        "x + y + z - 1"
    };
    std::vector<std::string> variables = {"x", "y", "z"};
    
    // Analyze dimension first
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Dimension analysis: dim=" << dim_info.dimension 
              << ", method=" << dim_info.method_used << std::endl;
    
    // Should be dimension 2 (plane in 3D)
    EXPECT_EQ(dim_info.dimension, 2);
    
    // Solve with hyperplane sections
    auto result = solve_via_hyperplane_sections(polynomials, variables, config);
    
    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_EQ(result.hyperplanes.size(), 2) << "Should need 2 hyperplanes for dimension 2";
    EXPECT_GT(result.solutions.solutions.size(), 0) << "Should find at least one solution";
    
    // Verify the solution satisfies the original equation
    if (!result.solutions.solutions.empty()) {
        auto& sol = result.solutions.solutions[0];
        double sum = sol[0].real() + sol[1].real() + sol[2].real();
        EXPECT_NEAR(sum, 1.0, 1e-6) << "Solution should satisfy x + y + z = 1";
    }
}

// Test 2: Two linear equations in 5 variables
// Dimension = 3, needs 3 hyperplanes
TEST_F(UnderdeterminedSystemTests, TwoLinearEquationsIn5D) {
    std::cout << "\n=== Testing: Two linear equations in 5 variables ===" << std::endl;
    
    std::vector<std::string> polynomials = {
        "a + b + c + d + e - 5",
        "a - b + c - d + e - 1"
    };
    std::vector<std::string> variables = {"a", "b", "c", "d", "e"};
    
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Dimension analysis: dim=" << dim_info.dimension 
              << ", method=" << dim_info.method_used << std::endl;
    
    // Should be dimension 3 (3-dimensional subspace in 5D)
    EXPECT_EQ(dim_info.dimension, 3);
    
    auto result = solve_via_hyperplane_sections(polynomials, variables, config);
    
    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_EQ(result.hyperplanes.size(), 3) << "Should need 3 hyperplanes for dimension 3";
    EXPECT_GT(result.solutions.solutions.size(), 0) << "Should find at least one solution";
}

// Test 3: Single quadratic equation in 4 variables (sphere in 4D)
// Dimension = 3, needs 3 hyperplanes
TEST_F(UnderdeterminedSystemTests, SphereIn4D) {
    std::cout << "\n=== Testing: Sphere in 4D ===" << std::endl;
    
    std::vector<std::string> polynomials = {
        "x^2 + y^2 + z^2 + w^2 - 4"
    };
    std::vector<std::string> variables = {"x", "y", "z", "w"};
    
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Dimension analysis: dim=" << dim_info.dimension 
              << ", method=" << dim_info.method_used << std::endl;
    
    // Should be dimension 3 (3-sphere in 4D)
    EXPECT_EQ(dim_info.dimension, 3);
    
    auto result = solve_via_hyperplane_sections(polynomials, variables, config);
    
    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_EQ(result.hyperplanes.size(), 3) << "Should need 3 hyperplanes for dimension 3";
    EXPECT_GT(result.solutions.solutions.size(), 0) << "Should find at least one solution";
    
    // Verify the solution satisfies the sphere equation
    if (!result.solutions.solutions.empty()) {
        auto& sol = result.solutions.solutions[0];
        double sum_squares = 0;
        for (int i = 0; i < 4; ++i) {
            double val = sol[i].real();
            sum_squares += val * val;
        }
        EXPECT_NEAR(sum_squares, 4.0, 1e-6) << "Solution should be on the sphere";
    }
}

// Test 4: Extremely underdetermined - single equation in 6 variables
// Dimension = 5, needs 5 hyperplanes
TEST_F(UnderdeterminedSystemTests, ExtremelyUnderdetermined6D) {
    std::cout << "\n=== Testing: Single equation in 6 variables ===" << std::endl;
    
    std::vector<std::string> polynomials = {
        "x1 + x2^2 + x3 + x4^2 + x5 + x6 - 10"
    };
    std::vector<std::string> variables = {"x1", "x2", "x3", "x4", "x5", "x6"};
    
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Dimension analysis: dim=" << dim_info.dimension 
              << ", method=" << dim_info.method_used << std::endl;
    
    // Should be dimension 5 (hypersurface in 6D)
    EXPECT_EQ(dim_info.dimension, 5);
    
    // For extremely high dimension, might need to adjust config
    config.max_dimension = 6;
    auto result = solve_via_hyperplane_sections(polynomials, variables, config);
    
    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_EQ(result.hyperplanes.size(), 5) << "Should need 5 hyperplanes for dimension 5";
    EXPECT_GT(result.solutions.solutions.size(), 0) << "Should find at least one solution";
}

// Test 5: Mixed system - one quadratic, one linear in 4 variables
// Dimension = 2, needs 2 hyperplanes
TEST_F(UnderdeterminedSystemTests, MixedQuadraticLinear4D) {
    std::cout << "\n=== Testing: Mixed quadratic/linear system in 4D ===" << std::endl;
    
    std::vector<std::string> polynomials = {
        "x^2 + y^2 - z^2 - 1",  // Hyperboloid
        "x + y + z + w - 2"      // Hyperplane
    };
    std::vector<std::string> variables = {"x", "y", "z", "w"};
    
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Dimension analysis: dim=" << dim_info.dimension 
              << ", method=" << dim_info.method_used << std::endl;
    
    // Should be dimension 2 (intersection of hyperboloid and hyperplane)
    EXPECT_EQ(dim_info.dimension, 2);
    
    auto result = solve_via_hyperplane_sections(polynomials, variables, config);
    
    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_EQ(result.hyperplanes.size(), 2) << "Should need 2 hyperplanes for dimension 2";
    EXPECT_GT(result.solutions.solutions.size(), 0) << "Should find at least one solution";
}

// Test 6: Wildly underdetermined with dependencies
// This is subtle - variables appear algebraically dependent
TEST_F(UnderdeterminedSystemTests, AlgebraicallyDependentVariables) {
    std::cout << "\n=== Testing: System with algebraically dependent variables ===" << std::endl;
    
    std::vector<std::string> polynomials = {
        "x*y - z",           // z depends on x,y
        "x^2 + y^2 - 4"      // Circle constraint on x,y
    };
    std::vector<std::string> variables = {"x", "y", "z"};
    
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Dimension analysis: dim=" << dim_info.dimension 
              << ", method=" << dim_info.method_used << std::endl;
    
    // Should be dimension 1 (parametric curve)
    // Even though we have 3 vars and 2 equations, z is determined by x,y
    EXPECT_EQ(dim_info.dimension, 1);
    
    auto result = solve_via_hyperplane_sections(polynomials, variables, config);
    
    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_EQ(result.hyperplanes.size(), 1) << "Should need 1 hyperplane for dimension 1";
    EXPECT_GT(result.solutions.solutions.size(), 0) << "Should find at least one solution";
}

// Test 7: Incremental hyperplane addition - test that we can handle many hyperplanes
TEST_F(UnderdeterminedSystemTests, IncrementalHyperplaneAddition) {
    std::cout << "\n=== Testing: Incremental hyperplane addition ===" << std::endl;
    
    // Start with empty system in 4 variables - dimension 4
    std::vector<std::string> polynomials = {};
    std::vector<std::string> variables = {"a", "b", "c", "d"};
    
    // Add constraints one by one
    for (int i = 0; i < 4; ++i) {
        if (i > 0) {
            // Add a polynomial constraint
            if (i == 1) {
                polynomials.push_back("a + b + c + d - 1");
            } else if (i == 2) {
                polynomials.push_back("a^2 + b^2 - 1");
            } else if (i == 3) {
                polynomials.push_back("c*d - 1");
            }
        }
        
        auto dim_info = analyze_system_dimension(polynomials, variables);
        std::cout << "After " << polynomials.size() << " constraint(s): dim=" 
                  << dim_info.dimension << std::endl;
        
        if (dim_info.dimension > 0) {
            // System is still positive dimensional
            auto result = solve_via_hyperplane_sections(polynomials, variables, config);
            
            if (result.success) {
                std::cout << "  Added " << result.hyperplanes.size() 
                          << " hyperplanes, found " << result.solutions.solutions.size() 
                          << " solutions" << std::endl;
            } else {
                std::cout << "  Failed: " << result.error_message << std::endl;
            }
        }
    }
}

// Test 8: Sample variety points - get multiple samples from positive-dimensional variety
TEST_F(UnderdeterminedSystemTests, SampleVarietyPoints) {
    std::cout << "\n=== Testing: Sampling points from variety ===" << std::endl;
    
    // Ellipse in 3D (dimension 1)
    std::vector<std::string> polynomials = {
        "x^2 + 2*y^2 - 4",
        "x + y + z - 1"
    };
    std::vector<std::string> variables = {"x", "y", "z"};
    
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Variety dimension: " << dim_info.dimension << std::endl;
    EXPECT_EQ(dim_info.dimension, 1);
    
    // Sample multiple points from the variety
    int num_samples = 5;
    auto samples = sample_variety_points(polynomials, variables, num_samples, config);
    
    int successful_samples = 0;
    for (const auto& sample : samples) {
        if (sample.success) {
            successful_samples++;
            std::cout << "Sample " << successful_samples << " found " 
                      << sample.solutions.solutions.size() << " points with hyperplanes: ";
            for (const auto& hp : sample.hyperplanes) {
                std::cout << "[" << hp << "] ";
            }
            std::cout << std::endl;
        }
    }
    
    EXPECT_GT(successful_samples, 0) << "Should get at least one successful sample";
    
    // Each sample should use different hyperplanes (randomized)
    if (samples.size() >= 2 && samples[0].success && samples[1].success) {
        EXPECT_NE(samples[0].hyperplanes, samples[1].hyperplanes) 
            << "Different samples should use different random hyperplanes";
    }
}

// Test 9: High-dimensional linear system
TEST_F(UnderdeterminedSystemTests, HighDimensionalLinearSystem) {
    std::cout << "\n=== Testing: High-dimensional linear system ===" << std::endl;
    
    // 2 equations in 7 variables - dimension 5
    std::vector<std::string> polynomials = {
        "x1 + x2 + x3 + x4 + x5 + x6 + x7 - 7",
        "x1 - x2 + x3 - x4 + x5 - x6 + x7 - 1"
    };
    std::vector<std::string> variables = {"x1", "x2", "x3", "x4", "x5", "x6", "x7"};
    
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Dimension analysis: dim=" << dim_info.dimension 
              << ", method=" << dim_info.method_used << std::endl;
    
    // Should be dimension 5
    EXPECT_EQ(dim_info.dimension, 5);
    
    config.max_dimension = 7;
    auto result = solve_via_hyperplane_sections(polynomials, variables, config);
    
    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_EQ(result.hyperplanes.size(), 5) << "Should need 5 hyperplanes for dimension 5";
    EXPECT_GT(result.solutions.solutions.size(), 0) << "Should find at least one solution";
}

// Test 10: Degenerate case - redundant constraints
TEST_F(UnderdeterminedSystemTests, RedundantConstraints) {
    std::cout << "\n=== Testing: System with redundant constraints ===" << std::endl;
    
    // These three equations define the same plane (scaled versions)
    std::vector<std::string> polynomials = {
        "x + y + z - 3",
        "2*x + 2*y + 2*z - 6",
        "3*x + 3*y + 3*z - 9"
    };
    std::vector<std::string> variables = {"x", "y", "z"};
    
    auto dim_info = analyze_system_dimension(polynomials, variables);
    std::cout << "Dimension analysis: dim=" << dim_info.dimension 
              << ", method=" << dim_info.method_used << std::endl;
    
    // Should still be dimension 2 (single plane)
    EXPECT_EQ(dim_info.dimension, 2);
    
    auto result = solve_via_hyperplane_sections(polynomials, variables, config);
    
    ASSERT_TRUE(result.success) << result.error_message;
    EXPECT_EQ(result.hyperplanes.size(), 2) << "Should need 2 hyperplanes for dimension 2";
    EXPECT_GT(result.solutions.solutions.size(), 0) << "Should find at least one solution";
}