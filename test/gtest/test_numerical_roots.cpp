#include "julia_rur/numerical_roots_eigen.hpp"
#include "julia_rur/numerical_roots_flint.hpp"
#include "test_helpers.hpp"

using namespace rur_test;

/**
 * @brief Test numerical root finding algorithms
 */
class NumericalRoots : public ::testing::Test {
  protected:
    void SetUp() override {
        // Set up common test polynomials
    }
};

TEST_F(NumericalRoots, EigenQuadraticRoots) {
    // Test Eigen root finder with x^2 - 4 (roots: ±2)
    std::vector<mpq_class> coeffs = {
        mpq_class(-4), // constant term
        mpq_class(0),  // x term
        mpq_class(1)   // x^2 term
    };

    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    EXPECT_EQ(roots.size(), 2);

    // Sort roots by real part
    std::sort(roots.begin(), roots.end(), [](const std::complex<double> &a, const std::complex<double> &b) {
        return a.real() < b.real();
    });

    EXPECT_TRUE(TestHelpers::approx_equal(roots[0], std::complex<double>(-2.0, 0.0)));
    EXPECT_TRUE(TestHelpers::approx_equal(roots[1], std::complex<double>(2.0, 0.0)));
}

TEST_F(NumericalRoots, EigenComplexRoots) {
    // Test x^2 + 1 (roots: ±i)
    std::vector<mpq_class> coeffs = {
        mpq_class(1), // constant term
        mpq_class(0), // x term
        mpq_class(1)  // x^2 term
    };

    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    EXPECT_EQ(roots.size(), 2);

    // Check that roots have real part ≈ 0 and imaginary part ≈ ±1
    for (const auto &root : roots) {
        EXPECT_TRUE(TestHelpers::approx_equal(root.real(), 0.0));
        EXPECT_TRUE(TestHelpers::approx_equal(std::abs(root.imag()), 1.0));
    }
}

TEST_F(NumericalRoots, EigenCubicRoots) {
    // Test x^3 - 8 (roots: 2, 2ω, 2ω^2 where ω = cube root of unity)
    std::vector<mpq_class> coeffs = {
        mpq_class(-8), // constant term
        mpq_class(0),  // x term
        mpq_class(0),  // x^2 term
        mpq_class(1)   // x^3 term
    };

    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    EXPECT_EQ(roots.size(), 3);

    // One root should be real (≈ 2), two should be complex
    int real_count = 0;
    for (const auto &root : roots) {
        if (std::abs(root.imag()) < TestHelpers::DEFAULT_TOLERANCE) {
            real_count++;
            EXPECT_TRUE(TestHelpers::approx_equal(root.real(), 2.0));
        } else {
            // Complex roots should have magnitude 2
            EXPECT_TRUE(TestHelpers::approx_equal(std::abs(root), 2.0, TestHelpers::LOOSE_TOLERANCE));
        }
    }

    EXPECT_EQ(real_count, 1);
}

TEST_F(NumericalRoots, EigenHighDegreeRoots) {
    // Test x^8 - 1 (8th roots of unity)
    std::vector<mpq_class> coeffs(9); // degree 8 + 1
    for (int i = 0; i < 8; i++) { coeffs[i] = mpq_class(0); }
    coeffs[0] = mpq_class(-1); // constant term
    coeffs[8] = mpq_class(1);  // x^8 term

    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    EXPECT_EQ(roots.size(), 8);

    // All roots should lie on the unit circle
    for (const auto &root : roots) {
        double magnitude = std::abs(root);
        EXPECT_TRUE(TestHelpers::approx_equal(magnitude, 1.0, TestHelpers::LOOSE_TOLERANCE))
          << "Root " << root << " has magnitude " << magnitude;
    }

    // Should have exactly 2 real roots (±1)
    int real_count = 0;
    for (const auto &root : roots) {
        if (std::abs(root.imag()) < TestHelpers::DEFAULT_TOLERANCE) {
            real_count++;
            EXPECT_TRUE(TestHelpers::approx_equal(std::abs(root.real()), 1.0));
        }
    }

    EXPECT_EQ(real_count, 2);
}

#ifdef HAVE_FLINT_ARB
TEST_F(NumericalRoots, FlintQuadraticRoots) {
    // Test FLINT root finder with x^2 - 4
    std::vector<mpq_class> coeffs = {
        mpq_class(-4), // constant term
        mpq_class(0),  // x term
        mpq_class(1)   // x^2 term
    };

    auto roots = julia_rur::find_polynomial_roots_flint(coeffs);

    EXPECT_EQ(roots.size(), 2);

    // Sort roots by real part
    std::sort(roots.begin(), roots.end(), [](const std::complex<double> &a, const std::complex<double> &b) {
        return a.real() < b.real();
    });

    EXPECT_TRUE(TestHelpers::approx_equal(roots[0], std::complex<double>(-2.0, 0.0)));
    EXPECT_TRUE(TestHelpers::approx_equal(roots[1], std::complex<double>(2.0, 0.0)));
}

TEST_F(NumericalRoots, FlintCertifiedRoots) {
    // Test FLINT certified root finding
    std::vector<mpq_class> coeffs = {
        mpq_class(-4), // constant term
        mpq_class(0),  // x term
        mpq_class(1)   // x^2 term
    };

    std::vector<double> error_bounds;
    auto roots = julia_rur::find_polynomial_roots_certified(coeffs, error_bounds);

    EXPECT_EQ(roots.size(), 2);
    EXPECT_EQ(error_bounds.size(), 2);

    // Error bounds should be small
    for (const auto &bound : error_bounds) { EXPECT_LT(bound, 1e-10); }
}

TEST_F(NumericalRoots, FlintHighPrecisionRoots) {
    // Test FLINT with custom precision settings
    julia_rur::FlintRootConfig config;
    config.initial_prec = 256;
    config.max_prec = 1024;
    config.epsilon = 1e-20;

    std::vector<mpq_class> coeffs = {
        mpq_class(-2), // constant term
        mpq_class(0),  // x term
        mpq_class(1)   // x^2 term
    };

    auto roots = julia_rur::find_polynomial_roots_flint(coeffs, config);

    EXPECT_EQ(roots.size(), 2);

    // Roots should be ±√2 with high precision
    for (const auto &root : roots) {
        double expected_magnitude = std::sqrt(2.0);
        EXPECT_TRUE(TestHelpers::approx_equal(std::abs(root.real()), expected_magnitude, 1e-12));
        EXPECT_TRUE(TestHelpers::approx_equal(root.imag(), 0.0, 1e-12));
    }
}
#endif

TEST_F(NumericalRoots, FallbackBehavior) {
    // Test that fallback to Eigen works when FLINT is not available
    std::vector<mpq_class> coeffs = {
        mpq_class(1), // constant term (x^2 + 1)
        mpq_class(0), // x term
        mpq_class(1)  // x^2 term
    };

    // This should work regardless of FLINT availability
    auto roots = julia_rur::find_polynomial_roots_flint(coeffs);

    EXPECT_EQ(roots.size(), 2);

    // Roots should be ±i
    for (const auto &root : roots) {
        EXPECT_TRUE(TestHelpers::approx_equal(root.real(), 0.0));
        EXPECT_TRUE(TestHelpers::approx_equal(std::abs(root.imag()), 1.0));
    }
}

TEST_F(NumericalRoots, RationalCoefficients) {
    // Test with rational coefficients: (2x-1)(3x+2) = 6x^2 + 4x - 3x - 2 = 6x^2 + x - 2
    std::vector<mpq_class> coeffs = {
        mpq_class(-2), // constant term
        mpq_class(1),  // x term
        mpq_class(6)   // x^2 term
    };

    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    EXPECT_EQ(roots.size(), 2);

    // Sort roots by real part
    std::sort(roots.begin(), roots.end(), [](const std::complex<double> &a, const std::complex<double> &b) {
        return a.real() < b.real();
    });

    // Roots should be -2/3 and 1/2
    EXPECT_TRUE(TestHelpers::approx_equal(roots[0].real(), -2.0 / 3.0));
    EXPECT_TRUE(TestHelpers::approx_equal(roots[0].imag(), 0.0));
    EXPECT_TRUE(TestHelpers::approx_equal(roots[1].real(), 0.5));
    EXPECT_TRUE(TestHelpers::approx_equal(roots[1].imag(), 0.0));
}

TEST_F(NumericalRoots, ZeroPolynomial) {
    // Test edge case: zero polynomial
    std::vector<mpq_class> coeffs = { mpq_class(0) };

    // This should either return empty or handle gracefully
    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    // The behavior may vary - either empty result or exception is acceptable
    // Just ensure it doesn't crash
    SUCCEED();
}

TEST_F(NumericalRoots, ConstantPolynomial) {
    // Test edge case: constant nonzero polynomial (no roots)
    std::vector<mpq_class> coeffs = { mpq_class(1) };

    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    // Should return empty vector (no roots)
    EXPECT_EQ(roots.size(), 0);
}

TEST_F(NumericalRoots, LinearPolynomial) {
    // Test x - 3 = 0 (root: x = 3)
    std::vector<mpq_class> coeffs = {
        mpq_class(-3), // constant term
        mpq_class(1)   // x term
    };

    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    EXPECT_EQ(roots.size(), 1);
    EXPECT_TRUE(TestHelpers::approx_equal(roots[0], std::complex<double>(3.0, 0.0)));
}

/**
 * @brief Comparative test between Eigen and FLINT (when available)
 */
TEST_F(NumericalRoots, CompareEigenAndFlint) {
    std::vector<mpq_class> coeffs = {
        mpq_class(-6), // constant term
        mpq_class(11), // x term
        mpq_class(-6), // x^2 term
        mpq_class(1)   // x^3 term
    };

    // Get Eigen results
    julia_rur::EnhancedSolverConfig config;
    std::optional<std::vector<double>> error_bounds;
    auto eigen_roots = julia_rur::find_roots_auto(coeffs, config, error_bounds);

    // Get FLINT results (may fallback to Eigen)
    auto flint_roots = julia_rur::find_polynomial_roots_flint(coeffs);

    EXPECT_EQ(eigen_roots.size(), flint_roots.size());

    // Sort both results for comparison
    std::sort(eigen_roots.begin(), eigen_roots.end(), [](const std::complex<double> &a, const std::complex<double> &b) {
        if (std::abs(a.real() - b.real()) > 1e-10) return a.real() < b.real();
        return a.imag() < b.imag();
    });

    std::sort(flint_roots.begin(), flint_roots.end(), [](const std::complex<double> &a, const std::complex<double> &b) {
        if (std::abs(a.real() - b.real()) > 1e-10) return a.real() < b.real();
        return a.imag() < b.imag();
    });

    // Results should be approximately equal
    for (size_t i = 0; i < eigen_roots.size(); i++) {
        EXPECT_TRUE(TestHelpers::approx_equal(eigen_roots[i], flint_roots[i], TestHelpers::LOOSE_TOLERANCE))
          << "Eigen root " << eigen_roots[i] << " != FLINT root " << flint_roots[i];
    }
}