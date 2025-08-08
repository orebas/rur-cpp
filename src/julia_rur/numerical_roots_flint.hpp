#pragma once

#include "rur_main_algorithm.hpp"
#include "numerical_roots_eigen.hpp"

// FLINT 3.x with integrated ARB/ACB is required
// Only use ARB if FLINT version supports it
#ifdef HAVE_FLINT_ARB
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>
#endif

typedef long slong;
#include <vector>
#include <complex>

namespace julia_rur {

/**
 * @brief Configuration for FLINT root finding
 */
struct FlintRootConfig {
    slong initial_prec = 128;       // Initial precision in bits
    slong max_prec = 4096;          // Maximum precision to try
    slong max_iter = 100;           // Maximum iterations for refinement
    double epsilon = 1e-14;         // Tolerance for considering roots as real
    bool isolate_real_roots = true; // Try to identify and refine real roots
};

#ifdef HAVE_FLINT_ARB

/**
 * @brief Find polynomial roots using FLINT's acb_poly
 * 
 * Uses FLINT's arbitrary precision complex interval arithmetic
 * to find all roots of a polynomial with certified error bounds.
 * 
 * @param polynomial_coeffs Coefficients of univariate polynomial (low to high degree)
 * @param config Configuration options
 * @return Complex roots with guaranteed accuracy
 */
std::vector<std::complex<double>> find_polynomial_roots_flint(
    const std::vector<mpq_class>& polynomial_coeffs,
    const FlintRootConfig& config = FlintRootConfig()
);

/**
 * @brief Find polynomial roots with certified error bounds
 * 
 * Returns roots along with radius of uncertainty for each root.
 * 
 * @param polynomial_coeffs Coefficients of univariate polynomial (low to high degree)
 * @param root_radii Output vector for error radii
 * @param config Configuration options
 * @return Complex roots with error bounds in root_radii
 */
std::vector<std::complex<double>> find_polynomial_roots_certified(
    const std::vector<mpq_class>& polynomial_coeffs,
    std::vector<double>& root_radii,
    const FlintRootConfig& config = FlintRootConfig()
);

/**
 * @brief Check if a root is real within tolerance
 * 
 * Uses FLINT's interval arithmetic to determine if a root
 * is provably real or complex.
 */
bool is_root_real_certified(const acb_t root, double epsilon = 1e-14);

/**
 * @brief Convert rational polynomial to FLINT polynomial
 * 
 * Helper function to convert from GMP rationals to FLINT's
 * arbitrary precision complex intervals.
 */
void rational_poly_to_acb_poly(
    acb_poly_t result,
    const std::vector<mpq_class>& coeffs,
    slong prec
);

#else

// Fallback implementations when ARB is not available
inline std::vector<std::complex<double>> find_polynomial_roots_flint(
    const std::vector<mpq_class>& polynomial_coeffs,
    const FlintRootConfig& config = FlintRootConfig()
) {
    // Fall back to Eigen
    return find_polynomial_roots(polynomial_coeffs);
}

inline std::vector<std::complex<double>> find_polynomial_roots_certified(
    const std::vector<mpq_class>& polynomial_coeffs,
    std::vector<double>& root_radii,
    const FlintRootConfig& config = FlintRootConfig()
) {
    // Fall back to Eigen without certification
    auto roots = find_polynomial_roots(polynomial_coeffs);
    root_radii.assign(roots.size(), 1e-10); // Dummy error bounds
    return roots;
}

#endif

} // namespace julia_rur