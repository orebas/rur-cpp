#include "julia_rur/f4_polynomial_formatter.hpp"
#include "julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <string>
#include <vector>


int
main() {
    std::cout << "Testing polynomial formatter..." << std::endl;

    // Test the formatter directly
    std::string poly = "x^2 - 2";
    std::string formatted = julia_rur::format_polynomial_for_f4(poly);
    std::cout << "Direct formatter test:" << std::endl;
    std::cout << "  Original: \"" << poly << "\"" << std::endl;
    std::cout << "  Formatted: \"" << formatted << "\"" << std::endl;

    // Also test with variations
    std::vector<std::string> test_cases = { "x^2 - 2",    "x^2-2", "1*x^2 - 2", "1*x^2-2",
                                            "x^2 + (-2)", "x - 3", "1*x - 3",   "x + (-3)" };

    std::cout << "\nTesting various formats:" << std::endl;
    for (const auto &test : test_cases) {
        std::string fmt = julia_rur::format_polynomial_for_f4(test);
        std::cout << "  \"" << test << "\" -> \"" << fmt << "\"" << std::endl;
    }

    // Now test with RUR algorithm
    std::cout << "\nTesting full RUR with verbose output:" << std::endl;
    std::vector<std::string> polynomials = { "x^2 - 2" };
    std::vector<std::string> variables = { "x" };

    julia_rur::RURConfig config;
    config.verbose = true;

    auto [mod_result, coeffs] = julia_rur::compute_modular_rur(polynomials, variables, 131071, config);

    if (mod_result.success) {
        std::cout << "RUR computation succeeded" << std::endl;
        std::cout << "Minimal polynomial has " << mod_result.minimal_polynomial.coefficients.size() << " coefficients"
                  << std::endl;
    } else {
        std::cout << "RUR computation failed" << std::endl;
    }

    return 0;
}