#pragma once

#include "julia_rur/numerical_roots_eigen.hpp"
#include "julia_rur/polynomial_evaluator.hpp"
#include "julia_rur/polynomial_solver_enhanced.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace rur_test {

/**
 * @brief Helper functions and utilities for polynomial system testing
 */
class TestHelpers {
  public:
    // Tolerance constants
    static constexpr double DEFAULT_TOLERANCE = 1e-8;
    static constexpr double LOOSE_TOLERANCE = 1e-6;
    static constexpr double TIGHT_TOLERANCE = 1e-10;

    /**
     * @brief Check if two complex numbers are approximately equal
     */
    static bool approx_equal(const std::complex<double> &a,
                             const std::complex<double> &b,
                             double tolerance = DEFAULT_TOLERANCE) {
        return std::abs(a - b) < tolerance;
    }

    /**
     * @brief Check if two real numbers are approximately equal
     */
    static bool approx_equal(double a, double b, double tolerance = DEFAULT_TOLERANCE) {
        return std::abs(a - b) < tolerance;
    }

    /**
     * @brief Compute residual for a solution of polynomial system
     */
    static double compute_residual(const std::vector<std::string> &polynomials,
                                   const std::vector<std::string> &variables,
                                   const std::vector<std::complex<double>> &solution) {
        return julia_rur::PolynomialEvaluator::compute_residual(polynomials, variables, solution);
    }

    /**
     * @brief Check if a solution satisfies the polynomial system within tolerance
     */
    static bool is_valid_solution(const std::vector<std::string> &polynomials,
                                  const std::vector<std::string> &variables,
                                  const std::vector<std::complex<double>> &solution,
                                  double tolerance = DEFAULT_TOLERANCE) {
        return compute_residual(polynomials, variables, solution) < tolerance;
    }

    /**
     * @brief Sort solutions for consistent comparison in tests
     */
    static void sort_solutions(std::vector<std::vector<std::complex<double>>> &solutions) {
        std::sort(solutions.begin(),
                  solutions.end(),
                  [](const std::vector<std::complex<double>> &a, const std::vector<std::complex<double>> &b) {
                      // Sort by first component, then by second, etc.
                      for (size_t i = 0; i < std::min(a.size(), b.size()); i++) {
                          if (std::abs(a[i].real() - b[i].real()) > 1e-10) { return a[i].real() < b[i].real(); }
                          if (std::abs(a[i].imag() - b[i].imag()) > 1e-10) { return a[i].imag() < b[i].imag(); }
                      }
                      return a.size() < b.size();
                  });
    }

    /**
     * @brief Count real solutions (imaginary part near zero)
     */
    static int count_real_solutions(const std::vector<std::vector<std::complex<double>>> &solutions,
                                    double tolerance = DEFAULT_TOLERANCE) {
        int count = 0;
        for (const auto &solution : solutions) {
            bool is_real = true;
            for (const auto &val : solution) {
                if (std::abs(val.imag()) > tolerance) {
                    is_real = false;
                    break;
                }
            }
            if (is_real) count++;
        }
        return count;
    }

    /**
     * @brief Extract real parts from solutions (assuming they are real)
     */
    static std::vector<std::vector<double>> extract_real_solutions(
      const std::vector<std::vector<std::complex<double>>> &solutions,
      double tolerance = DEFAULT_TOLERANCE) {

        std::vector<std::vector<double>> real_solutions;
        for (const auto &solution : solutions) {
            bool is_real = true;
            std::vector<double> real_solution;

            for (const auto &val : solution) {
                if (std::abs(val.imag()) > tolerance) {
                    is_real = false;
                    break;
                }
                real_solution.push_back(val.real());
            }

            if (is_real) { real_solutions.push_back(real_solution); }
        }
        return real_solutions;
    }

    // Verify residuals for all solutions are below tolerance
    static void expect_all_residuals_below(const std::vector<std::string> &polynomials,
                                           const std::vector<std::string> &variables,
                                           const std::vector<std::vector<std::complex<double>>> &solutions,
                                           double tolerance = DEFAULT_TOLERANCE) {
        for (size_t i = 0; i < solutions.size(); ++i) {
            double residual = compute_residual(polynomials, variables, solutions[i]);
            EXPECT_LT(residual, tolerance) << "Solution " << i << " residual=" << residual;
        }
    }

    // Verify pairwise distinctness of solutions using max component-wise separation
    static void expect_pairwise_distinct(const std::vector<std::vector<std::complex<double>>> &solutions,
                                         double separation_tolerance = 1e-8) {
        for (size_t i = 0; i < solutions.size(); ++i) {
            for (size_t j = i + 1; j < solutions.size(); ++j) {
                double max_comp_dist = 0.0;
                size_t dim = std::min(solutions[i].size(), solutions[j].size());
                for (size_t k = 0; k < dim; ++k) {
                    max_comp_dist = std::max(max_comp_dist, std::abs(solutions[i][k] - solutions[j][k]));
                }
                EXPECT_GT(max_comp_dist, separation_tolerance)
                  << "Solutions " << i << " and " << j << " too close (potential duplicate): " << max_comp_dist;
            }
        }
    }
};

/**
 * @brief Base test fixture for polynomial system tests
 */
class PolynomialSystemTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Default configuration for tests
        config_.verbose = false; // Reduce output during tests
    }

    /**
     * @brief Solve polynomial system with standard configuration
     */
    julia_rur::EnhancedPolynomialSolution solve_system(const std::vector<std::string> &polynomials,
                                                       const std::vector<std::string> &variables) {

        return julia_rur::solve_polynomial_system_enhanced(polynomials, variables, config_);
    }

    /**
     * @brief Assert that solution is successful and has expected number of solutions
     */
    void assert_solution_valid(const julia_rur::EnhancedPolynomialSolution &solution, int expected_count = -1) {
        ASSERT_TRUE(solution.success) << "Solver failed: " << solution.error_message;

        if (expected_count >= 0) {
            EXPECT_EQ(solution.solutions.size(), expected_count)
              << "Expected " << expected_count << " solutions, got " << solution.solutions.size();
        }
    }

    julia_rur::EnhancedSolverConfig config_;
};

/**
 * @brief Parametrized test fixture for testing multiple similar cases
 */
template<typename T>
class ParametrizedPolynomialTest
  : public PolynomialSystemTest
  , public ::testing::WithParamInterface<T> {};

/**
 * @brief Test case structure for polynomial systems
 */
struct PolynomialTestCase {
    std::string name;
    std::vector<std::string> polynomials;
    std::vector<std::string> variables;
    int expected_solutions;
    double tolerance;

    PolynomialTestCase(const std::string &n,
                       const std::vector<std::string> &p,
                       const std::vector<std::string> &v,
                       int expected = -1,
                       double tol = TestHelpers::DEFAULT_TOLERANCE)
      : name(n)
      , polynomials(p)
      , variables(v)
      , expected_solutions(expected)
      , tolerance(tol) {}
};

} // namespace rur_test