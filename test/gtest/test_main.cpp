#include <gtest/gtest.h>

/**
 * @brief Main entry point for Google Test suite
 * 
 * This is the unified test runner for the RUR C++ polynomial solver project.
 * It automatically discovers and runs all tests defined across the test files.
 * 
 * Usage:
 *   ./rur_tests                           # Run all tests
 *   ./rur_tests --gtest_filter="*Linear*" # Run only linear system tests
 *   ./rur_tests --gtest_list_tests        # List all available tests
 *   ./rur_tests --gtest_help              # Show all options
 * 
 * Test Categories:
 *   - BasicOperations.*    : Basic polynomial and monomial operations
 *   - UnivariateSystem.*   : Single-variable polynomial systems
 *   - LinearSystem.*       : Linear equation systems
 *   - NonlinearSystem.*    : Nonlinear polynomial systems
 *   - NumericalRoots.*     : Root-finding algorithms (Eigen/FLINT)
 *   - Integration.*        : End-to-end system tests
 */

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    
    // Enable colored output if supported
    ::testing::GTEST_FLAG(color) = "auto";
    
    return RUN_ALL_TESTS();
}