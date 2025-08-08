#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>


int
main() {
    std::cout << "Testing x^2 - 2 with small prime 131063..." << std::endl;

    std::vector<std::string> polynomials = { "x^2 - 2" };
    std::vector<std::string> variables = { "x" };

    julia_rur::RURConfig config;
    config.verbose = true;

    julia_rur::ModularCoeff prime = 131063;
    auto [result, coeffs] = julia_rur::compute_modular_rur(polynomials, variables, prime, config);

    if (result.success) {
        std::cout << "\nMinimal polynomial coefficients:" << std::endl;
        for (size_t i = 0; i < result.minimal_polynomial.coefficients.size(); ++i) {
            julia_rur::ModularCoeff coeff = result.minimal_polynomial.coefficients[i];
            std::cout << "  coeff[" << i << "] = " << coeff;

            if (coeff > prime / 2) {
                int64_t signed_val = static_cast<int64_t>(coeff) - static_cast<int64_t>(prime);
                std::cout << " (represents " << signed_val << ")";
            }
            std::cout << std::endl;
        }
    }

    return 0;
}