#include "../src/julia_rur/polynomial_solver.hpp"
#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>


using namespace julia_rur;

void
test_single_variable() {
    std::cout << "=== Testing single variable system: x^2 - 2 ===\n";

    // Simple system: x^2 - 2 = 0
    std::vector<std::string> polynomials = { "1*x^2-2" };
    std::vector<std::string> variables = { "x" };

    // Solve
    PolynomialSystemSolution solutions = solve_polynomial_system_complete(polynomials, variables);

    if (solutions.success) {
        std::cout << "System solved successfully!\n";
        std::cout << "Number of solutions: " << solutions.solutions.size() << "\n";
        for (size_t i = 0; i < solutions.solutions.size(); ++i) {
            const auto &sol = solutions.solutions[i];
            std::cout << "Solution " << (i + 1) << ": x = " << sol[0] << "\n";
        }
    } else {
        std::cerr << "Failed to solve system: " << solutions.error_message << "\n";
    }
}

void
test_multivariate() {
    std::cout << "\n=== Testing multivariate system: circle-parabola ===\n";

    // Circle-parabola: x^2 + y^2 - 4 = 0, y - x^2 = 0
    std::vector<std::string> polynomials = { "1*x^2+1*y^2-4", "1*y-1*x^2" };
    std::vector<std::string> variables = { "x", "y" };

    RURConfig config;
    config.verbose = true;

    // Test modular RUR first
    ModularCoeff prime = 1073741827;
    std::cout << "\nTesting modular RUR with prime " << prime << ":\n";

    auto [mod_result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (mod_result.success) {
        std::cout << "Modular RUR successful!\n";
        std::cout << "Minimal polynomial degree: " << mod_result.minimal_polynomial.degree << "\n";
        std::cout << "Quotient basis size: " << mod_result.quotient_basis.size() << "\n";
    } else {
        std::cerr << "Modular RUR failed!\n";
    }
}

int
main() {
    test_single_variable();
    test_multivariate();
    return 0;
}