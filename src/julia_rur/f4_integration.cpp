#include "f4_integration.hpp"
#include "polynomial_operations.hpp" // For modular_inverse
#include <iostream>
#include <numeric>

namespace julia_rur {

// DO NOT normalize to symmetric range - keep coefficients as F4 provides them
// F4 already gives us coefficients in the correct modular representation
inline ModularCoeff
normalize_symmetric(ModularCoeff coeff, ModularCoeff prime) {
    // Just return the coefficient as-is
    // For x^2 - 2, F4 gives us [1, prime-2] which is correct
    return coeff;
}

// Degree reverse lexicographic ordering
bool
degrevlex_greater(const PP &a, const PP &b) {
    // First compare total degree
    uint32_t deg_a = std::accumulate(a.begin(), a.end(), 0u);
    uint32_t deg_b = std::accumulate(b.begin(), b.end(), 0u);

    if (deg_a != deg_b) { return deg_a > deg_b; }

    // Same degree - compare lexicographically from the end (reverse)
    for (int i = a.size() - 1; i >= 0; --i) {
        if (a[i] != b[i]) {
            return a[i] < b[i]; // Note: reversed comparison for degrevlex
        }
    }

    return false; // Equal
}

bool
extract_f4_groebner_basis(std::vector<std::vector<PP>> &groebner_exponents,
                          std::vector<std::vector<ModularCoeff>> &groebner_coefficients,
                          int num_variables,
                          ModularCoeff prime) {
    int basis_size = axf4_get_basis_size();
    // std::cout << "extract_f4_groebner_basis: basis_size = " << basis_size << std::endl;

    if (basis_size <= 0) {
        std::cerr << "extract_f4_groebner_basis: invalid basis size" << std::endl;
        return false;
    }

    // Initialize output vectors
    groebner_exponents.clear();
    groebner_exponents.resize(basis_size);
    groebner_coefficients.clear();
    groebner_coefficients.resize(basis_size);

    // Extract each polynomial from F4
    for (int poly_idx = 0; poly_idx < basis_size; ++poly_idx) {
        // Get polynomial data from F4
        std::vector<unsigned int> coeffs;
        std::vector<unsigned int> monomials;

        // Get the term count
        int term_count = axf4_get_poly_term_count(poly_idx);
        // std::cout << "Poly " << poly_idx << " has " << term_count << " terms" << std::endl;

        if (term_count <= 0) {
            std::cerr << "Skipping empty polynomial " << poly_idx << std::endl;
            continue; // Skip empty polynomials
        }

        // Allocate space and get actual data
        coeffs.resize(term_count);
        monomials.resize(term_count);
        int ret = axf4_get_poly_data(poly_idx, coeffs.data(), monomials.data());

        if (ret != 0) {
            std::cerr << "Failed to get poly data for polynomial " << poly_idx << std::endl;
            return false; // Error getting data
        }

        // Convert to our format
        groebner_exponents[poly_idx].resize(term_count);
        groebner_coefficients[poly_idx].resize(term_count);

        for (int term_idx = 0; term_idx < term_count; ++term_idx) {
            // Decode F4 monomial to power product
            PP pp = decode_f4_monomial(monomials[term_idx], num_variables);
            groebner_exponents[poly_idx][term_idx] = pp;

            // Normalize coefficient to symmetric range
            ModularCoeff coeff = static_cast<ModularCoeff>(coeffs[term_idx]);
            groebner_coefficients[poly_idx][term_idx] = normalize_symmetric(coeff, prime);

            // Debug output disabled
            // std::cout << "  Poly " << poly_idx << " Term " << term_idx
            //           << ": monomial=" << monomials[term_idx]
            //           << " -> PP=[";
            // for (size_t j = 0; j < pp.size(); ++j) {
            //     if (j > 0) std::cout << ",";
            //     std::cout << pp[j];
            // }
            // std::cout << "], coeff=" << coeffs[term_idx] << std::endl;
        }
    }

    return true;
}

bool
f4_to_multiplication_tables(axf4_session_t session,
                            std::vector<std::vector<ModularCoeff>> &t_v,
                            std::vector<StackVect> &t_xw,
                            std::vector<std::vector<int32_t>> &i_xw,
                            std::vector<PP> &quotient_basis,
                            ModularCoeff prime) {
    if (!session) {
        std::cerr << "f4_to_multiplication_tables: null session" << std::endl;
        return false;
    }

    // Step 1: Extract F4 Gröbner basis data
    std::vector<std::vector<PP>> groebner_exponents;
    std::vector<std::vector<ModularCoeff>> groebner_coefficients;

    int num_vars = axf4_get_num_variables(session);
    // std::cout << "f4_to_multiplication_tables: num_vars = " << num_vars << std::endl;

    if (num_vars <= 0 || !extract_f4_groebner_basis(groebner_exponents, groebner_coefficients, num_vars, prime)) {
        std::cerr << "f4_to_multiplication_tables: failed to extract groebner basis" << std::endl;
        return false;
    }

    // Step 2: Extract leading terms from Gröbner basis
    std::vector<PP> leading_terms;
    leading_terms.reserve(groebner_exponents.size());

    for (size_t i = 0; i < groebner_exponents.size(); ++i) {
        const auto &poly_exps = groebner_exponents[i];
        if (!poly_exps.empty()) {
            // F4 should sort in ascending order, but for linear polynomials
            // we need to find the actual leading term in degrevlex order
            PP lt = poly_exps[0];
            for (size_t j = 1; j < poly_exps.size(); ++j) {
                if (degrevlex_greater(poly_exps[j], lt)) { lt = poly_exps[j]; }
            }
            leading_terms.push_back(lt);

            // Debug output disabled
            // std::cout << "Poly " << i << " leading term (idx=" << leading_idx << "): [";
        }
    }

    // Step 3: Compute quotient basis
    const bool verbose_f4 = false;
    if (verbose_f4) {
        std::cout << "DEBUG: Computing quotient basis from " << leading_terms.size() << " leading terms:" << std::endl;
    }
    for (size_t i = 0; i < leading_terms.size(); ++i) {
        std::cout << "  LT[" << i << "] = [";
        for (size_t j = 0; j < leading_terms[i].size(); ++j) {
            if (j > 0) std::cout << ",";
            std::cout << leading_terms[i][j];
        }
        std::cout << "]" << std::endl;
    }

    try {
        quotient_basis = compute_quotient_basis(leading_terms);
    } catch (const std::runtime_error &e) {
        // Check if this is the "GB = {1}" case (system has no solutions mod p)
        std::string error_msg = e.what();
        if (error_msg.find("Groebner basis = {1}") != std::string::npos) {
            // This is a bad prime - system has no solutions modulo this prime
            // This is expected for some primes, so we don't print an error
            return false;
        }
        std::cerr << "Error computing quotient basis: " << e.what() << std::endl;
        return false;
    } catch (const std::exception &e) {
        std::cerr << "Error computing quotient basis: " << e.what() << std::endl;
        return false;
    }

    if (verbose_f4) { std::cout << "DEBUG: Computed quotient basis size = " << quotient_basis.size() << std::endl; }
    for (size_t i = 0; i < quotient_basis.size(); ++i) {
        std::cout << "  QB[" << i << "] = [";
        for (size_t j = 0; j < quotient_basis[i].size(); ++j) {
            if (j > 0) std::cout << ",";
            std::cout << quotient_basis[i][j];
        }
        std::cout << "]" << std::endl;
    }

    // Step 4: Prepare multiplication table structure
    prepare_table_mxi(leading_terms, quotient_basis, t_xw, i_xw);

    // Debug: print GB elements (disabled for now)
    /*
    std::cout << "DEBUG f4_integration: GB has " << groebner_exponents.size() << " polynomials" << std::endl;
    for (size_t i = 0; i < groebner_exponents.size() && i < 3; ++i) {
        std::cout << "  GB[" << i << "]: " << groebner_exponents[i].size() << " terms" << std::endl;
        for (size_t j = 0; j < groebner_exponents[i].size() && j < 5; ++j) {
            std::cout << "    Term " << j << ": exp=[";
            for (size_t k = 0; k < groebner_exponents[i][j].size(); ++k) {
                if (k > 0) std::cout << ",";
                std::cout << groebner_exponents[i][j][k];
            }
            std::cout << "], coeff=" << groebner_coefficients[i][j] << std::endl;
        }
    }
    */

    // Step 5: Fill multiplication tables from Gröbner basis
    compute_fill_quo_gb(t_v, t_xw, groebner_exponents, groebner_coefficients, quotient_basis, prime);

    // Step 6: Complete multiplication tables
    learn_compute_table(t_v, t_xw, i_xw, prime);

    // Special handling for 1-dimensional quotient ring
    // In this case, the GB should contain elements like x-c, y-d, z-e
    // and we need to extract the constant values directly
    if (quotient_basis.size() == 1 && num_vars > 0) {
        // std::cout << "DEBUG: 1-dim quotient ring, extracting values from GB ("
        //           << groebner_exponents.size() << " elements)" << std::endl;
        // For each GB element, check if it's a linear polynomial of the form variable - constant
        for (size_t gb_idx = 0; gb_idx < groebner_exponents.size(); ++gb_idx) {
            const auto &exps = groebner_exponents[gb_idx];
            const auto &coeffs = groebner_coefficients[gb_idx];

            // Check if this is a linear polynomial with exactly 2 terms
            if (exps.size() == 2) {
                // Find which term is the variable and which is the constant
                int var_term = -1;
                int const_term = -1;

                for (size_t t = 0; t < 2; ++t) {
                    int var_count = 0;
                    int var_idx = -1;
                    for (size_t v = 0; v < exps[t].size(); ++v) {
                        if (exps[t][v] == 1) {
                            var_count++;
                            var_idx = v;
                        } else if (exps[t][v] != 0) {
                            var_count = -1; // Not a linear term
                            break;
                        }
                    }

                    if (var_count == 1) {
                        var_term = t;
                    } else if (var_count == 0) {
                        const_term = t;
                    }
                }

                // If we found a polynomial of the form variable - constant
                if (var_term >= 0 && const_term >= 0) {
                    // Find which variable this is
                    int var_idx = -1;
                    for (size_t v = 0; v < exps[var_term].size(); ++v) {
                        if (exps[var_term][v] == 1) {
                            var_idx = v;
                            break;
                        }
                    }

                    if (var_idx >= 0 && var_idx < num_vars) {
                        // Get the constant value (negated because polynomial is var - const = 0)
                        ModularCoeff const_val = coeffs[const_term];
                        ModularCoeff var_coeff = coeffs[var_term];

                        // std::cout << "  Found linear GB element: var[" << var_idx << "] * "
                        //           << var_coeff << " + const " << const_val << " = 0" << std::endl;

                        // If the coefficient of the variable term is not 1, we need to adjust
                        if (var_coeff != 1) {
                            // Solve var_coeff * x - const_val = 0 for x
                            // x = const_val / var_coeff (mod prime)
                            ModularCoeff inv_coeff = modular_inverse(var_coeff, prime);
                            const_val = (static_cast<AccModularCoeff>(const_val) * inv_coeff) % prime;
                        }

                        // The value is -const_val/var_coeff
                        // Since polynomial is var_coeff*x + const_val = 0, x = -const_val/var_coeff
                        ModularCoeff value = (prime - const_val) % prime;

                        // std::cout << "    Variable " << var_idx << " = " << value
                        //           << " (mod " << prime << ")" << std::endl;

                        // Store this value in the appropriate position in t_v
                        // The index should correspond to variable multiplication by 1
                        if (var_idx < static_cast<int>(i_xw.size()) && !i_xw[var_idx].empty()) {
                            int32_t table_idx = i_xw[var_idx][0];
                            if (table_idx > 0 && table_idx <= static_cast<int32_t>(t_v.size())) {
                                t_v[table_idx - 1].resize(1);
                                t_v[table_idx - 1][0] = value;
                                // std::cout << "    Stored in t_v[" << (table_idx-1) << "]" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }

    return true;
}

} // namespace julia_rur