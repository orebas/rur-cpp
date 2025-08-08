#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>

using namespace julia_rur;

// Override to disable verification during reconstruction
namespace julia_rur {

std::pair<bool, mpz_class>
crt_and_rational_reconstruction_no_verify(std::vector<std::vector<mpq_class>> &qq_result,
                                          std::vector<std::vector<mpz_class>> &zz_temp,
                                          const std::vector<std::vector<std::vector<ModularCoeff>>> &modular_tables,
                                          const std::vector<ModularCoeff> &primes) {
    size_t num_primes = primes.size();
    size_t num_polys = modular_tables[0].size();

    // Initialize if needed
    if (qq_result.size() != num_polys) {
        qq_result.resize(num_polys);
        zz_temp.resize(num_polys);
        for (size_t i = 0; i < num_polys; ++i) {
            size_t poly_size = modular_tables[0][i].size();
            qq_result[i].resize(poly_size, mpq_class(0));
            zz_temp[i].resize(poly_size, mpz_class(0));
        }
    }

    // Compute current modulus
    mpz_class modulus = 1;
    for (ModularCoeff p : primes) { modulus *= p; }

    auto bounds_balanced = compute_balanced_bounds(modulus);
    std::vector<mpz_class> crt_multipliers;

    bool success = true;
    mpz_class lcm_denominator = 1;

    // Process each polynomial
    for (size_t i = 0; i < num_polys && success; ++i) {
        size_t poly_size = modular_tables[0][i].size();

        // Process each coefficient
        for (size_t j = 0; j < poly_size && success; ++j) {
            // Collect remainders for this coefficient across all primes
            std::vector<ModularCoeff> remainders(num_primes);
            for (size_t k = 0; k < num_primes; ++k) { remainders[k] = modular_tables[k][i][j]; }

            // Apply CRT
            chinese_remainder_theorem(zz_temp[i][j], remainders, primes, crt_multipliers);

            // Normalize to symmetric range (-m/2, m/2]
            if (zz_temp[i][j] > modulus / 2) { zz_temp[i][j] -= modulus; }

            // Debug print
            if (i == 0) { // First polynomial (minimal polynomial)
                std::cout << "Coeff " << j << ": zz=" << zz_temp[i][j] << ", modulus=" << modulus
                          << ", bounds: N=" << bounds_balanced.N << " D=" << bounds_balanced.D << std::endl;

                // Try reconstruction to see what's happening
                // First try with the actual bounds
                auto test_rat = rational_reconstruction(zz_temp[i][j], modulus, bounds_balanced.N, bounds_balanced.D);
                std::cout << "  Direct reconstruction: success=" << test_rat.success;
                if (test_rat.success) { std::cout << ", value=" << test_rat.rational; }

                // The issue is our value is -131061 but bounds are only 92674
                // The correct value after CRT should be -2
                // Let's check the remainders
                std::cout << "\n  Remainders: [" << remainders[0] << ", " << remainders[1] << "]" << std::endl;
            }

            // Try simple balanced reconstruction (no verification)
            auto rat_result = rational_reconstruction_with_denominator_no_verify(
              zz_temp[i][j], mpz_class(1), modulus, bounds_balanced.N, bounds_balanced.D);

            if (!rat_result.success) {
                std::cout << "Failed to reconstruct poly " << i << " coeff " << j << std::endl;
                success = false;
                break;
            }

            qq_result[i][j] = rat_result.rational;

            // Update LCM of denominators
            mpz_lcm(
              lcm_denominator.get_mpz_t(), lcm_denominator.get_mpz_t(), rat_result.rational.get_den().get_mpz_t());
        }
    }

    return { success, lcm_denominator };
}

} // namespace julia_rur

int
main() {
    std::cout << "Testing rational reconstruction without verification" << std::endl;

    std::vector<std::string> polynomials = { "1*x^2-2" };
    std::vector<std::string> variables = { "x" };

    // First compute modular results
    std::vector<ModularCoeff> primes = { 131063, 131059 };
    std::vector<ModularRURResult> mod_results;

    RURConfig config;
    config.verbose = false;

    for (ModularCoeff p : primes) {
        auto [result, coeffs] = compute_modular_rur(polynomials, variables, p, config, {});
        if (result.success) { mod_results.push_back(result); }
    }

    std::cout << "Computed " << mod_results.size() << " modular results" << std::endl;

    // Prepare data for CRT
    std::vector<std::vector<std::vector<ModularCoeff>>> modular_tables;
    for (const auto &mod_result : mod_results) {
        std::vector<std::vector<ModularCoeff>> prime_table;

        // Debug: print minimal polynomial coefficients
        std::cout << "Prime " << mod_result.prime << " minimal poly coeffs:";
        for (auto c : mod_result.minimal_polynomial.coefficients) { std::cout << " " << c; }
        std::cout << std::endl;

        prime_table.push_back(mod_result.minimal_polynomial.coefficients);
        for (const auto &param : mod_result.parameterizations) {
            if (!param.generators.empty()) { prime_table.push_back(param.generators[0]); }
        }
        modular_tables.push_back(prime_table);
    }

    // Try reconstruction
    std::vector<std::vector<mpq_class>> qq_result;
    std::vector<std::vector<mpz_class>> zz_temp;

    auto [success, denom] = crt_and_rational_reconstruction_no_verify(qq_result, zz_temp, modular_tables, primes);

    if (success) {
        std::cout << "\nReconstruction succeeded!" << std::endl;
        std::cout << "Minimal polynomial:";
        for (const auto &c : qq_result[0]) { std::cout << " " << c; }
        std::cout << std::endl;
    } else {
        std::cout << "\nReconstruction failed" << std::endl;
    }

    return 0;
}