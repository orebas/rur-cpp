#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

void
test_tensor_product_system() {
    std::cout << "=== Testing Tensor Product System x^3=1, y^2=1 ===\n\n";

    std::vector<std::string> polynomials = { "1*x^3-1", "1*y^2-1" };
    std::vector<std::string> variables = { "x", "y" };

    // Test with the standard prime
    ModularCoeff prime = 1073741827;

    std::cout << "Testing with prime " << prime << ":\n";

    RURConfig config;
    config.verbose = true;

    auto [result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (result.success) {
        std::cout << "\nSUCCESS! Found separating element\n";
        std::cout << "Minimal polynomial degree: " << result.minimal_polynomial.degree << "\n";

        // Print the minimal polynomial coefficients
        std::cout << "Minimal polynomial coefficients: ";
        for (size_t i = 0; i < result.minimal_polynomial.coefficients.size(); ++i) {
            std::cout << result.minimal_polynomial.coefficients[i] << " ";
        }
        std::cout << "\n";
    } else {
        std::cout << "\nFailed to find separating element\n";
    }

    // Also test with a smaller prime to see behavior
    std::cout << "\n\n=== Testing with smaller prime 31 ===\n";
    ModularCoeff small_prime = 31;

    auto [result2, coeffs2] = compute_modular_rur(polynomials, variables, small_prime, config);

    if (result2.success) {
        std::cout << "\nSUCCESS with prime " << small_prime << "!\n";
        std::cout << "Minimal polynomial degree: " << result2.minimal_polynomial.degree << "\n";
    } else {
        std::cout << "\nFailed with prime " << small_prime << "\n";
    }
}

int
main() {
    test_tensor_product_system();
    return 0;
}