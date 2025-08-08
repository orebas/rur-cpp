#include "julia_rur/numerical_roots_eigen.hpp"
#include "julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <vector>


int
main() {
    std::cout << "Testing minimal polynomial computation for x^2 - 2..." << std::endl;

    std::vector<std::string> polynomials = { "x^2 - 2" };
    std::vector<std::string> variables = { "x" };

    julia_rur::RURConfig config;
    config.verbose = false;

    // Test with a single prime
    julia_rur::ModularCoeff prime = 131071;
    std::cout << "\nComputing modular RUR with prime " << prime << std::endl;

    auto [mod_result, coeffs] = julia_rur::compute_modular_rur(polynomials, variables, prime, config);

    if (mod_result.success) {
        std::cout << "Success! Minimal polynomial coefficients (modular):" << std::endl;
        for (size_t i = 0; i < mod_result.minimal_polynomial.coefficients.size(); ++i) {
            julia_rur::ModularCoeff coeff = mod_result.minimal_polynomial.coefficients[i];
            std::cout << "  coeff[" << i << "] = " << coeff;

            // Show signed interpretation
            if (coeff > prime / 2) {
                int64_t signed_val = static_cast<int64_t>(coeff) - static_cast<int64_t>(prime);
                std::cout << " (represents " << signed_val << ")";
            }
            std::cout << std::endl;
        }

        std::cout << "\nExpected: coeff[0] = " << (prime - 2) << " (represents -2)" << std::endl;
        std::cout << "          coeff[1] = 0" << std::endl;
        std::cout << "          coeff[2] = 1" << std::endl;
    }

    // Now test rational reconstruction
    std::cout << "\n\nTesting rational RUR computation..." << std::endl;
    julia_rur::RationalRURResult rat_result = julia_rur::compute_rational_rur(polynomials, variables, config);

    if (rat_result.success) {
        std::cout << "Success! Minimal polynomial coefficients (rational):" << std::endl;
        for (size_t i = 0; i < rat_result.minimal_polynomial.size(); ++i) {
            std::cout << "  coeff[" << i << "] = " << rat_result.minimal_polynomial[i] << std::endl;
        }

        std::cout << "\nExpected: coeff[0] = -2, coeff[1] = 0, coeff[2] = 1" << std::endl;
        std::cout << "This represents the polynomial: T^2 - 2 = 0" << std::endl;

        // Test root finding
        std::cout << "\nFinding roots of minimal polynomial..." << std::endl;
        auto roots = julia_rur::find_polynomial_roots(rat_result.minimal_polynomial);

        std::cout << "Found " << roots.size() << " roots:" << std::endl;
        for (size_t i = 0; i < roots.size(); ++i) { std::cout << "  root[" << i << "] = " << roots[i] << std::endl; }

        std::cout << "\nExpected: ±√2 ≈ ±1.41421356..." << std::endl;
    }

    return 0;
}