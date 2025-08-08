#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

int
main() {
    std::cout << "Debug: Linear system x - 3 = 0" << std::endl;

    std::vector<std::string> polynomials = { "1*x-3" };
    std::vector<std::string> variables = { "x" };

    // Test with a single prime first
    ModularCoeff prime = 65537;
    RURConfig config;
    config.verbose = true;

    std::cout << "\n=== Modular computation mod " << prime << " ===" << std::endl;
    auto [mod_result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (mod_result.success) {
        std::cout << "\nModular result:" << std::endl;
        std::cout << "Quotient basis size: " << mod_result.quotient_basis.size() << std::endl;
        std::cout << "Minimal polynomial coefficients: ";
        for (auto c : mod_result.minimal_polynomial.coefficients) { std::cout << c << " "; }
        std::cout << std::endl;

        // The polynomial x - 3 = 0 should give:
        // - Quotient basis: {1}
        // - GB: {x - 3}
        // - The variable x reduces to the constant 3 in the quotient ring
        // - So the minimal polynomial should be T - 3 = 0

        // Let's check what value x takes in the quotient ring
        // In the quotient ring R[x]/<x-3>, we have x = 3

        std::cout << "\nExpected: T - 3 = 0 (coefficients: [-3, 1] or [65534, 1] mod 65537)" << std::endl;
        std::cout << "Note: -3 mod 65537 = " << (65537 - 3) << std::endl;
    }

    // Now test rational reconstruction
    std::cout << "\n=== Rational computation ===" << std::endl;
    config.verbose = false;
    auto rat_result = compute_rational_rur(polynomials, variables, config);

    if (rat_result.success) {
        std::cout << "\nRational minimal polynomial: ";
        for (auto c : rat_result.minimal_polynomial) { std::cout << c << " "; }
        std::cout << std::endl;
        std::cout << "Expected: [-3, 1] giving T - 3 = 0" << std::endl;
    }

    return 0;
}