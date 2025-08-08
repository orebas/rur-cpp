#include "../src/julia_rur/multi_modular.hpp"
#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <chrono>
#include <iostream>


using namespace julia_rur;

int
main() {
    std::cout << "Testing multivariate system with different primes" << std::endl;

    // Simple system: x^2 + y^2 - 1 = 0, x - y = 0
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2-1", // Circle
        "1*x-1*y"        // Line x = y
    };
    std::vector<std::string> variables = { "x", "y" };

    // Test with various primes
    std::vector<ModularCoeff> test_primes = {
        65537,      // 2^16 + 1
        131071,     // 2^17 - 1
        262147,     // 2^18 + 3
        524287,     // 2^19 - 1
        1048583,    // First prime > 2^20
        2097169,    // First prime > 2^21
        1073741827, // Large prime
        1073741783, // Another large prime
        1073741789  // Yet another large prime
    };

    RURConfig config;
    config.verbose = false;

    for (ModularCoeff prime : test_primes) {
        std::cout << "\n=== Testing prime " << prime << " ===" << std::endl;

        auto start = std::chrono::high_resolution_clock::now();

        try {
            auto [mod_result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

            if (mod_result.success) {
                std::cout << "✓ Success in " << duration.count() << " ms" << std::endl;
                std::cout << "  Quotient dimension: " << mod_result.quotient_basis.size() << std::endl;
                std::cout << "  Minimal polynomial degree: " << mod_result.minimal_polynomial.degree << std::endl;
            } else {
                std::cout << "✗ Failed after " << duration.count() << " ms" << std::endl;
            }
        } catch (const std::exception &e) { std::cerr << "✗ Exception: " << e.what() << std::endl; }
    }

    // Now test if running multiple primes in sequence causes issues
    std::cout << "\n=== Testing multiple primes in sequence ===" << std::endl;

    for (size_t i = 0; i < 3; ++i) {
        std::cout << "\nRound " << (i + 1) << ":" << std::endl;
        for (size_t j = 0; j < 3 && j < test_primes.size(); ++j) {
            ModularCoeff prime = test_primes[j];
            std::cout << "  Prime " << prime << ": ";

            auto start = std::chrono::high_resolution_clock::now();

            try {
                auto [mod_result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

                if (mod_result.success) {
                    std::cout << "✓ " << duration.count() << " ms" << std::endl;
                } else {
                    std::cout << "✗ failed" << std::endl;
                }
            } catch (const std::exception &e) { std::cerr << "✗ exception" << std::endl; }
        }
    }

    return 0;
}