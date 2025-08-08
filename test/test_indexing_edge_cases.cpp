#include "../src/julia_rur/multiplication_tables.hpp"
#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>


using namespace julia_rur;

void
test_three_variable_system() {
    std::cout << "=== Testing 3-Variable System Indexing ===\n\n";

    // System: x^2 - 1 = 0, y^2 - 1 = 0, z^2 - 1 = 0
    // Quotient basis should be {1, x, y, z, xy, xz, yz, xyz}
    std::vector<std::string> polynomials = { "1*x^2-1", "1*y^2-1", "1*z^2-1" };
    std::vector<std::string> variables = { "x", "y", "z" };

    RURConfig config;
    config.verbose = false;

    ModularCoeff prime = 1073741827;
    auto [result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (!result.success) {
        std::cerr << "Failed to compute modular RUR\n";
        return;
    }

    std::cout << "Quotient basis size: " << result.quotient_basis.size() << "\n";
    std::cout << "Expected: 8 (2^3 = 8)\n\n";

    std::cout << "Quotient basis elements:\n";
    for (size_t i = 0; i < result.quotient_basis.size(); ++i) {
        const auto &pp = result.quotient_basis[i];
        std::cout << "  [" << i << "] = (";
        for (size_t j = 0; j < pp.size(); ++j) {
            if (j > 0) std::cout << ",";
            std::cout << pp[j];
        }
        std::cout << ") = ";

        // Interpret the monomial
        if (std::all_of(pp.begin(), pp.end(), [](int x) { return x == 0; })) {
            std::cout << "1";
        } else {
            bool first = true;
            std::vector<std::string> var_names = { "x", "y", "z" };
            for (size_t j = 0; j < pp.size() && j < var_names.size(); ++j) {
                if (pp[j] > 0) {
                    if (!first) std::cout << "*";
                    std::cout << var_names[j];
                    if (pp[j] > 1) std::cout << "^" << pp[j];
                    first = false;
                }
            }
        }
        std::cout << "\n";
    }

    // Test element_to_vector for all variables
    std::cout << "\nVariable representations: [DISABLED - needs i_xw and t_v parameters]\n";
}

void
test_variable_ordering() {
    std::cout << "\n=== Testing Variable Ordering ===\n\n";

    // System with specific variable ordering: z > y > x lexicographically
    // This tests if our indexing respects the variable ordering
    std::vector<std::string> polynomials = {
        "1*x-1", // x = 1
        "1*y-2", // y = 2
        "1*z-3"  // z = 3
    };
    std::vector<std::string> variables = { "x", "y", "z" };

    RURConfig config;
    config.verbose = false;

    ModularCoeff prime = 1073741827;
    auto [result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (!result.success) {
        std::cerr << "Failed to compute modular RUR\n";
        return;
    }

    std::cout << "For system where x=1, y=2, z=3:\n";
    std::cout << "Quotient basis should be {1} only\n";
    std::cout << "Actual quotient basis size: " << result.quotient_basis.size() << "\n";

    if (result.quotient_basis.size() == 1) {
        std::cout << "✓ Correct: System reduced to constants as expected\n";
    } else {
        std::cout << "✗ ERROR: Unexpected quotient basis size\n";
    }
}

void
test_mixed_degrees() {
    std::cout << "\n=== Testing Mixed Degree System ===\n\n";

    // System with different degrees: x^3 - 1 = 0, y^2 - 1 = 0
    // This creates a more complex quotient basis structure
    std::vector<std::string> polynomials = {
        "1*x^3-1", // x^3 = 1
        "1*y^2-1"  // y^2 = 1
    };
    std::vector<std::string> variables = { "x", "y" };

    RURConfig config;
    config.verbose = false;

    ModularCoeff prime = 1073741827;
    auto [result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (!result.success) {
        std::cerr << "Failed to compute modular RUR\n";
        return;
    }

    std::cout << "For x^3=1, y^2=1:\n";
    std::cout << "Expected quotient basis size: 6 (= 3*2)\n";
    std::cout << "Actual quotient basis size: " << result.quotient_basis.size() << "\n";

    std::cout << "\nQuotient basis:\n";
    for (size_t i = 0; i < result.quotient_basis.size(); ++i) {
        const auto &pp = result.quotient_basis[i];
        std::cout << "  [" << i << "] degree in x: " << pp[0] << ", degree in y: " << pp[1] << "\n";
    }

    // Check that highest degrees are x^2 and y^1
    bool correct_degrees = true;
    for (const auto &pp : result.quotient_basis) {
        if (pp[0] > 2) {
            std::cout << "ERROR: Found x degree > 2 (should be reduced by x^3=1)\n";
            correct_degrees = false;
        }
        if (pp[1] > 1) {
            std::cout << "ERROR: Found y degree > 1 (should be reduced by y^2=1)\n";
            correct_degrees = false;
        }
    }

    if (correct_degrees) { std::cout << "✓ All degrees are within expected bounds\n"; }
}

int
main() {
    test_three_variable_system();
    test_variable_ordering();
    test_mixed_degrees();
    return 0;
}