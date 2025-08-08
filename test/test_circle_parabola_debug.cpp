#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

void
test_circle_parabola() {
    std::cout << "=== Testing circle-parabola system ===\n";

    // Circle-parabola: x^2 + y^2 - 4 = 0, y - x^2 = 0
    std::vector<std::string> polynomials = { "1*x^2+1*y^2-4", "1*y-1*x^2" };
    std::vector<std::string> variables = { "x", "y" };

    RURConfig config;
    config.verbose = true;

    // Test with a single prime first
    ModularCoeff prime = 1073741827;

    auto [result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (result.success) {
        std::cout << "\nModular RUR successful!\n";
        std::cout << "Minimal polynomial degree: " << result.minimal_polynomial.degree << "\n";

        // Print parameterizations
        for (size_t i = 0; i < variables.size(); ++i) {
            std::cout << "\nParameterization for " << variables[i] << ":\n";
            if (i < result.parameterizations.size() && result.parameterizations[i].success) {
                std::cout << "  Success: " << result.parameterizations[i].generators.size() << " generators\n";
            } else {
                std::cout << "  Failed!\n";
            }
        }
    } else {
        std::cout << "\nModular RUR failed!\n";

        // Let's manually check what separating elements might work
        std::cout << "\nTrying to find a working separating element...\n";

        // The system should have 4 solutions corresponding to:
        // x = ±sqrt((-1 + sqrt(17))/2), y = (-1 + sqrt(17))/2
        // x = ±sqrt((-1 - sqrt(17))/2), y = (-1 - sqrt(17))/2

        // A good separating element might be x + y or x + 2y
        std::cout << "\nNote: For this system, we might need x+y or similar linear combinations\n";
    }
}

int
main() {
    test_circle_parabola();
    return 0;
}