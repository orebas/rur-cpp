#include "../src/julia_rur/f4_integration.hpp"
#include "../src/julia_rur/f4_polynomial_formatter.hpp"
#include "../src/julia_rur/multiplication_tables.hpp"
#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>


using namespace julia_rur;

void
test_variable_representation() {
    std::cout << "=== Testing Variable Representation in Quotient Basis ===\n\n";

    // Test system: x^2 - y = 0, y^2 - x = 0
    // This system has quotient basis {1, x, y, xy}
    std::vector<std::string> polynomials = { "1*x^2-1*y", "1*y^2-1*x" };
    std::vector<std::string> variables = { "x", "y" };

    RURConfig config;
    config.verbose = true;

    ModularCoeff prime = 1073741827;
    auto [result, coeffs] = compute_modular_rur(polynomials, variables, prime, config);

    if (!result.success) {
        std::cerr << "Failed to compute modular RUR\n";
        return;
    }

    std::cout << "\nQuotient basis analysis:\n";
    std::cout << "Quotient basis size: " << result.quotient_basis.size() << "\n";
    std::cout << "Elements:\n";
    for (size_t i = 0; i < result.quotient_basis.size(); ++i) {
        std::cout << "  [" << i << "] = ";
        const auto &pp = result.quotient_basis[i];
        std::cout << "(";
        for (size_t j = 0; j < pp.size(); ++j) {
            if (j > 0) std::cout << ",";
            std::cout << pp[j];
        }
        std::cout << ")";

        // Interpret the monomial
        std::cout << " = ";
        if (pp[0] == 0 && pp[1] == 0)
            std::cout << "1";
        else {
            bool first = true;
            if (pp[0] > 0) {
                std::cout << "x";
                if (pp[0] > 1) std::cout << "^" << pp[0];
                first = false;
            }
            if (pp[1] > 0) {
                if (!first) std::cout << "*";
                std::cout << "y";
                if (pp[1] > 1) std::cout << "^" << pp[1];
            }
        }
        std::cout << "\n";
    }

    // Test element_to_vector for each variable
    std::cout << "\nTesting element_to_vector: [DISABLED - needs i_xw and t_v parameters]\n";
}

void
test_multiplication_indices() {
    std::cout << "\n=== Testing Multiplication Table Indices ===\n\n";

    // Simple system: x^2 - 1 = 0, y^2 - 1 = 0
    // Quotient basis should be {1, x, y, xy}
    std::vector<std::string> polynomials = { "1*x^2-1", "1*y^2-1" };
    std::vector<std::string> variables = { "x", "y" };

    RURConfig config;
    config.verbose = false;

    ModularCoeff prime = 1073741827;

    // Get the multiplication tables directly
    auto session = axf4_create_session(prime, nullptr, 0);

    // Add variables
    std::vector<const char *> var_ptrs;
    for (const auto &var : variables) { var_ptrs.push_back(var.c_str()); }

    session = axf4_create_session(prime, var_ptrs.data(), variables.size());

    // Add polynomials
    for (const auto &poly : polynomials) {
        std::string formatted = format_polynomial_for_f4(poly);
        axf4_add_polynomial(session, formatted.c_str());
    }

    // Compute GB
    auto gb_result = axf4_compute_groebner_basis_keep_data(session);

    if (gb_result.status == 0) {
        std::vector<std::vector<ModularCoeff>> t_v;
        std::vector<StackVect> t_xw;
        std::vector<std::vector<int32_t>> i_xw;
        std::vector<PP> quotient_basis;

        bool success = f4_to_multiplication_tables(session, t_v, t_xw, i_xw, quotient_basis, prime);

        if (success) {
            std::cout << "Multiplication table i_xw:\n";
            for (size_t var = 0; var < i_xw.size(); ++var) {
                std::cout << "Variable " << (var + 1) << " multiplication indices: [";
                for (size_t i = 0; i < i_xw[var].size(); ++i) {
                    if (i > 0) std::cout << ", ";
                    std::cout << i_xw[var][i];
                }
                std::cout << "]\n";
            }

            std::cout << "\nChecking multiplication results:\n";
            // For quotient basis {1, x, y, xy}, test multiplications
            std::cout << "1 * x = x (should go to position 1)\n";
            std::cout << "x * x = 1 (should reduce via x^2 - 1 = 0)\n";
            std::cout << "1 * y = y (should go to position 2)\n";
            std::cout << "x * y = xy (should go to position 3)\n";
        }
    }

    axf4_free_result(&gb_result);
    axf4_cleanup_basis_data();
    axf4_destroy_session(session);
}

int
main() {
    test_variable_representation();
    test_multiplication_indices();
    return 0;
}