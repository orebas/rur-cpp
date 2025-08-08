#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

int
main() {
    std::cout << "=== Testing Parabola-Line System ===" << std::endl;

    // System: y = x^2, y = x + 1
    // Equivalent to: x^2 - y = 0, y - x - 1 = 0
    std::vector<std::string> polynomials = { "x^2 - y", "y - x - 1" };
    std::vector<std::string> variables = { "x", "y" };

    RURConfig config;
    config.verbose = true;

    std::cout << "\n--- Computing with first prime ---" << std::endl;
    ModularCoeff prime = 1073741827;
    auto [result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (result.success) {
        std::cout << "\nSUCCESS with prime " << prime << std::endl;
        std::cout << "Quotient basis size: " << result.quotient_basis.size() << std::endl;
        std::cout << "Minimal polynomial degree: " << result.minimal_polynomial.degree << std::endl;
    } else {
        std::cout << "\nFAILED with prime " << prime << std::endl;
    }

    return 0;
}