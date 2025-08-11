#include "multiplication_tables.hpp"
#include "polynomial_operations.hpp" // for modular_inverse
#include <stdexcept>

namespace julia_rur {

DividesResult
divides_with_var(const PP &a, const PP &b) {
    if (a.size() != b.size()) { throw std::invalid_argument("Monomial dimensions must match for division check"); }

    int32_t var = 0;
    for (size_t j = 0; j < a.size(); ++j) {
        if (a[j] > b[j]) {
            var = static_cast<int32_t>(j + 1); // Julia uses 1-based indexing
        }
        if (a[j] < b[j]) {
            // Division fails: return meaningful var index where failure occurred
            return DividesResult(false, static_cast<int32_t>(j + 1));
        }
    }
    return DividesResult(true, var);
}

PP
mul_pp_by_var(const PP &m, int32_t var_index) {
    PP result = m; // Copy the monomial
    if (var_index >= 1 && var_index <= static_cast<int32_t>(m.size())) {
        result[var_index - 1] += 1; // Convert to 0-based indexing
    }
    return result;
}

BorderSearchResult
find_in_border(const PP &m, const std::vector<StackVect> &t_xw) {
    using namespace power_product;

    // TODO OPTIMIZATION: This function has O(b) complexity where b is border size
    // TODO OPTIMIZATION: Could use degree-indexed data structures for O(log b) lookup
    // TODO OPTIMIZATION: Or pre-sort border by degree and use binary search

    int32_t pos = static_cast<int32_t>(t_xw.size());
    int32_t res_flag = 0;
    int32_t res_dd = 0; // Changed from PP to int32_t
    int32_t res_pos = 0;

    // Get total degree of target monomial
    uint32_t tm = total_degree(m);

    // Search backwards through border (by decreasing degrees to limit tests)
    while (pos > 0) {
        const PP &mm = t_xw[pos - 1].mon; // Convert to 0-based indexing
        uint32_t tmm = total_degree(mm);

        if (tmm == tm) {
            // Same degree - check for exact match
            if (m == mm) {
                return BorderSearchResult(1, 0, pos); // For exact match, var_index not used
            }
        } else if (tmm == tm - 1) {
            // One degree lower - check if this could be a predecessor
            // Test if mm is not in the quotient basis
            // (quotient elements have prev > 0 and var == 0)
            if (!((t_xw[pos - 1].prev > 0) && (t_xw[pos - 1].var == 0))) {
                // Fix: Call divides_with_var(m, mm) to match Julia's divides(m, mm)
                DividesResult div_result = divides_with_var(m, mm);
                if (div_result.divides) {
                    res_flag = 2;
                    res_dd = div_result.var_index; // Now directly an int32_t
                    res_pos = pos;
                }
            }
        } else {
            // Total degree is too high to find a predecessor
            if (res_flag > 0) {
                return BorderSearchResult(res_flag, res_dd, res_pos);
            } else {
                throw std::runtime_error("Error find in border: " + std::to_string(tm) + " " + std::to_string(tmm));
            }
        }
        pos = pos - 1;
    }

    // End of search
    if (res_flag > 0) {
        return BorderSearchResult(res_flag, res_dd, res_pos);
    } else {
        throw std::runtime_error("Error in finding a predecessor");
    }
}

void
prepare_table_mxi(const std::vector<PP> &ltg,             // GB leading terms
                  const std::vector<PP> &kb,              // Quotient basis
                  std::vector<StackVect> &t_xw,           // Border structure (output)
                  std::vector<std::vector<int32_t>> &i_xw // Variable indices (output)
) {
    if (ltg.empty() || kb.empty()) { throw std::invalid_argument("Leading terms and quotient basis cannot be empty"); }

    // TODO OPTIMIZATION: Replace linear searches with hash maps for O(1) lookup
    // TODO OPTIMIZATION: Use std::unordered_map<PP, int32_t> for kb and ltg lookups
    // TODO OPTIMIZATION: Use std::unordered_map<PP, int32_t> for t_xw monomial lookups
    // TODO OPTIMIZATION: This would reduce complexity from O(n²·m) to O(n·m)

    // Get number of variables from first leading term
    int32_t nbv = static_cast<int32_t>(ltg[0].size());

    // Initialize tablex (i_xw) - one vector per variable, one entry per quotient basis element
    i_xw.clear();
    i_xw.resize(nbv);
    for (int32_t i = 0; i < nbv; ++i) { i_xw[i].resize(kb.size(), 0); }

    // Initialize border structure
    t_xw.clear();
    int32_t nb_stack = 0;

    // For each quotient basis element
    for (size_t j = 0; j < kb.size(); ++j) {
        PP m = kb[j]; // Copy the monomial

        // For each variable (Julia iterates in reverse order: nbv down to 1)
        for (int32_t ii = 1; ii <= nbv; ++ii) {
            int32_t i = nbv - ii + 1; // Convert to Julia's reverse iteration

            // Multiply m by variable i: nm = m * x_i
            PP nm = mul_pp_by_var(m, i);

            // Check if nm is in quotient basis (kb)
            // TODO OPTIMIZATION: Replace with hash map lookup for O(1) instead of O(k)
            auto pos_iter = std::find(kb.begin(), kb.end(), nm);
            if (pos_iter == kb.end()) {
                // Not in quotient basis - check if it's in Gröbner basis leading terms
                // TODO OPTIMIZATION: Replace with hash map lookup for O(1) instead of O(g)
                auto gb_pos_iter = std::find(ltg.begin(), ltg.end(), nm);
                if (gb_pos_iter == ltg.end()) {
                    // Not in GB either - need to search border
                    BorderSearchResult border_result = find_in_border(nm, t_xw);

                    if (border_result.flag == 1) {
                        // Found exact match in border
                        i_xw[i - 1][j] = border_result.pos;
                    } else if (border_result.flag == 2) {
                        // Found predecessor - add new border element
                        nb_stack = nb_stack + 1;
                        i_xw[i - 1][j] = nb_stack;
                        // For border elements from find_in_border, var field contains the variable index
                        t_xw.emplace_back(nb_stack, nm, border_result.pos, border_result.var_index);
                    } else {
                        throw std::runtime_error("Error search table");
                    }
                } else {
                    // nm is a leading monomial of an element of the GB
                    // Insert with flags prev=0 and var=pos in gb
                    int32_t gb_pos = static_cast<int32_t>(gb_pos_iter - ltg.begin() + 1); // 1-based

                    // Check if already exists in border (using reverse iterators for findlast)
                    // TODO OPTIMIZATION: Replace with hash map lookup for O(1) instead of O(b)
                    auto existing =
                      std::find_if(t_xw.rbegin(), t_xw.rend(), [&nm](const StackVect &sv) { return sv.mon == nm; });

                    if (existing == t_xw.rend()) {
                        // Not found - add new border element
                        nb_stack = nb_stack + 1;
                        i_xw[i - 1][j] = nb_stack;
                        t_xw.emplace_back(nb_stack, ltg[gb_pos - 1], 0, gb_pos);
                    } else {
                        // Found existing - use its position
                        i_xw[i - 1][j] = existing->pos;
                    }
                }
            } else {
                // nm is an element of the quotient basis
                // Insert with flags prev=pos and var=0
                int32_t quo_pos = static_cast<int32_t>(pos_iter - kb.begin() + 1); // 1-based

                // Check if already exists in border (using reverse iterators for findlast)
                // TODO OPTIMIZATION: Replace with hash map lookup for O(1) instead of O(b)
                auto existing =
                  std::find_if(t_xw.rbegin(), t_xw.rend(), [&nm](const StackVect &sv) { return sv.mon == nm; });

                if (existing == t_xw.rend()) {
                    // Not found - add new border element
                    nb_stack = nb_stack + 1;
                    i_xw[i - 1][j] = nb_stack;
                    t_xw.emplace_back(nb_stack, kb[quo_pos - 1], quo_pos, 0);
                } else {
                    // Found existing - use its position
                    i_xw[i - 1][j] = existing->pos;
                }
            }
        }
    }
}

void
learn_compute_table(std::vector<std::vector<ModularCoeff>> &t_v,   // Coefficient vectors (output)
                    const std::vector<StackVect> &t_xw,            // Border structure
                    const std::vector<std::vector<int32_t>> &i_xw, // Variable indices
                    ModularCoeff prime                             // Modular arithmetic prime
) {
    // t_v should already be initialized - don't clear it!
    // Just ensure it has the right size
    if (t_v.size() != t_xw.size()) { t_v.resize(t_xw.size()); }

    // Julia algorithm: iteratively fill coefficient vectors
    int32_t nb = 1;               // Number of elements computed in this iteration
    std::vector<int32_t> t_learn; // Track which elements were computed

    // Buffer for accumulator arithmetic (Julia uses AccModularCoeff)
    std::vector<AccModularCoeff> buf;

    while (nb > 0) {
        nb = 0;

        // For each border element
        for (size_t i = 0; i < t_xw.size(); ++i) {
            // Check if already computed
            if (t_v[i].empty()) {
                bool continuer = false;

                // Check if the ancestor/predecessor is computed
                int32_t prev_idx = t_xw[i].prev;
                if (prev_idx > 0 && prev_idx <= static_cast<int32_t>(t_v.size()) &&
                    !t_v[prev_idx - 1].empty()) { // Convert to 0-based indexing

                    // Get the variable to multiply by
                    int32_t var_idx = t_xw[i].var;

                    // Multiply predecessor by variable to get this element
                    std::vector<ModularCoeff> result;
                    bool success = mul_var_quo_internal(result,
                                                        t_v[prev_idx - 1], // Predecessor coefficient vector
                                                        var_idx,           // Variable to multiply by
                                                        i_xw,              // Variable indices
                                                        t_v,               // Coefficient vectors
                                                        prime,             // Modular prime
                                                        buf                // Accumulator buffer
                    );

                    if (success) {
                        t_v[i] = std::move(result);
                        t_learn.push_back(static_cast<int32_t>(i + 1)); // 1-based for Julia compatibility
                        nb = nb + 1;
                    }
                }
            }
        }
    }
}

// Internal helper function for multiplication with accumulator buffer
bool
mul_var_quo_internal(std::vector<ModularCoeff> &result,                 // Output vector
                     const std::vector<ModularCoeff> &input,            // Input vector in quotient basis
                     int32_t var_index,                                 // Variable to multiply by (1-based)
                     const std::vector<std::vector<int32_t>> &i_xw,     // Pre-computed indices
                     const std::vector<std::vector<ModularCoeff>> &t_v, // Coefficient vectors
                     ModularCoeff prime,                                // Modular arithmetic prime
                     std::vector<AccModularCoeff> &buf                  // Accumulator buffer
) {
    size_t dim = input.size();
    result.resize(dim);

    // Optimized reduction frequency based on prime size
    // Larger primes can handle more operations before reduction
    const int32_t pack = (prime > (1UL << 30)) ? 8 :   // 31-bit primes: reduce more often
                         (prime > (1UL << 28)) ? 16 :   // 29-30 bit primes: standard
                         (prime > (1UL << 20)) ? 32 :   // 21-28 bit primes: less frequent
                         64;                            // Small primes: least frequent

    // Initialize accumulator buffer
    buf.resize(dim);
    for (size_t i = 0; i < dim; ++i) { buf[i] = 0; }

    bool continuer = true;

    // Convert var_index to 0-based for array access
    int32_t var_idx_0based = var_index - 1;

    // Check bounds
    if (var_idx_0based < 0 || var_idx_0based >= static_cast<int32_t>(i_xw.size())) { return false; }

    // For each element in input vector
    for (size_t j = 0; j < dim; ++j) {
        if (input[j] == 0) continue; // Skip zero coefficients

        // Get the index where xi * mj maps to
        if (j >= i_xw[var_idx_0based].size()) {
            continuer = false;
            break;
        }

        int32_t target_idx = i_xw[var_idx_0based][j];
        if (target_idx <= 0 || target_idx > static_cast<int32_t>(t_v.size())) {
            continuer = false;
            break;
        }

        // Convert to 0-based indexing
        target_idx -= 1;

        const auto &target_coeffs = t_v[target_idx];

        if (target_coeffs.size() > 1) {
            // Multiple coefficients: add_mul! operation
            // buf += input[j] * target_coeffs
            for (size_t k = 0; k < target_coeffs.size() && k < dim; ++k) {
                buf[k] += static_cast<AccModularCoeff>(input[j]) * static_cast<AccModularCoeff>(target_coeffs[k]);
            }

            // Periodic reduction to prevent overflow
            if ((j + 1) % pack == 0) {
                for (size_t k = 0; k < dim; ++k) { buf[k] %= prime; }
            }
        } else if (target_coeffs.size() == 1) {
            // Single coefficient: direct addition
            uint32_t kk = target_coeffs[0]; // kk is a 1-based index
            if (kk > 0 && kk <= dim) {
                // Convert 1-based kk to 0-based index for buffer access
                buf[kk - 1] = (buf[kk - 1] + static_cast<AccModularCoeff>(input[j])) % prime;
            }
        } else {
            // Empty coefficient vector: computation not ready
            continuer = false;
            break;
        }
    }

    // Final reduction: convert accumulator to modular coefficients
    for (size_t i = 0; i < dim; ++i) { result[i] = static_cast<ModularCoeff>(buf[i] % prime); }

    return continuer;
}

void
mul_var_quo(std::vector<ModularCoeff> &result,                 // Output vector
            const std::vector<ModularCoeff> &input,            // Input vector in quotient basis
            int32_t var_index,                                 // Variable to multiply by
            const std::vector<std::vector<int32_t>> &i_xw,     // Pre-computed indices
            const std::vector<std::vector<ModularCoeff>> &t_v, // Coefficient vectors
            ModularCoeff prime                                 // Modular arithmetic prime
) {
    std::vector<AccModularCoeff> buf;
    mul_var_quo_internal(result, input, var_index, i_xw, t_v, prime, buf);
}

void
vectorize_polynomial_in_quotient_basis(const std::vector<PP> &exponents,
                                       const std::vector<ModularCoeff> &coefficients,
                                       const std::vector<PP> &quotient_basis,
                                       std::vector<ModularCoeff> &result,
                                       ModularCoeff prime,
                                       bool skip_leading_term) {
    // Initialize result vector with zeros
    std::fill(result.begin(), result.end(), 0);

    // Debug: print input polynomial
    bool debug = false;
    if (debug && exponents.size() <= 3) {
        std::cout << "Vectorizing polynomial with " << exponents.size() << " terms:" << std::endl;
        for (size_t i = 0; i < exponents.size(); ++i) {
            std::cout << "  Term " << i << ": coeff=" << coefficients[i] << ", exp=[";
            for (auto e : exponents[i]) std::cout << e << " ";
            std::cout << "]" << std::endl;
        }
    }

    // Julia's vectorize_pol_gro! starts from index 2, skipping the leading term
    // This is crucial for GB elements
    size_t start_idx = skip_leading_term ? 1 : 0;

    // For each term in the polynomial (skipping leading if requested)
    for (size_t term_idx = start_idx; term_idx < exponents.size(); ++term_idx) {
        const PP &monomial = exponents[term_idx];
        ModularCoeff coeff = coefficients[term_idx];

        // Find this monomial in the quotient basis
        auto it = std::find(quotient_basis.begin(), quotient_basis.end(), monomial);
        if (it != quotient_basis.end()) {
            size_t basis_idx = std::distance(quotient_basis.begin(), it);
            // Julia does: res[pos] = (ModularPrime(arithm) % Coeff) - p_coeffs[i]
            // We negate the coefficient when vectorizing GB elements
            if (skip_leading_term) {
                // Negate for GB element (lt = -tail)
                result[basis_idx] = (result[basis_idx] + (prime - coeff)) % prime;
            } else {
                // Normal addition
                result[basis_idx] = (result[basis_idx] + coeff) % prime;
            }
        } else if (debug) {
            std::cout << "  Warning: monomial not in quotient basis" << std::endl;
        }
        // Note: If monomial not in quotient basis, it should be reducible by GB
        // This would indicate an error in the GB computation or quotient basis extraction
    }
}

// Overload for backward compatibility
void
vectorize_polynomial_in_quotient_basis(const std::vector<PP> &exponents,
                                       const std::vector<ModularCoeff> &coefficients,
                                       const std::vector<PP> &quotient_basis,
                                       std::vector<ModularCoeff> &result,
                                       ModularCoeff prime) {
    vectorize_polynomial_in_quotient_basis(exponents, coefficients, quotient_basis, result, prime, false);
}

void
compute_fill_quo_gb(std::vector<std::vector<ModularCoeff>> &t_v,
                    const std::vector<StackVect> &t_xw,
                    const std::vector<std::vector<PP>> &groebner_exponents,
                    const std::vector<std::vector<ModularCoeff>> &groebner_coefficients,
                    const std::vector<PP> &quotient_basis,
                    ModularCoeff prime) {
    // Initialize output vectors
    t_v.clear();
    t_v.resize(t_xw.size());

    // Process each element in the border structure
    for (size_t i = 0; i < t_xw.size(); ++i) {
        const StackVect &element = t_xw[i];

        if (element.var > 0 && element.prev == 0) {
            // This is a Gröbner basis element
            // var contains the 1-based index into the Gröbner basis
            size_t gb_index = static_cast<size_t>(element.var - 1); // Convert to 0-based

            // std::cout << "DEBUG: Processing GB element " << gb_index << " for border element " << i << std::endl;

            if (gb_index < groebner_exponents.size()) {
                // For a GB element g(x) = 0, we need to express the leading monomial
                // in terms of the other monomials. If g(x) = lt(g) + tail(g), then
                // lt(g) = -tail(g) in the quotient ring.

                // Debug: print the GB element (disabled)
                /*
                std::cout << "  GB element has " << groebner_exponents[gb_index].size() << " terms" << std::endl;
                for (size_t j = 0; j < groebner_exponents[gb_index].size(); ++j) {
                    std::cout << "    Term " << j << ": exp=[";
                    for (size_t k = 0; k < groebner_exponents[gb_index][j].size(); ++k) {
                        if (k > 0) std::cout << ",";
                        std::cout << groebner_exponents[gb_index][j][k];
                    }
                    std::cout << "], coeff=" << groebner_coefficients[gb_index][j] << std::endl;
                }
                */

                // CRITICAL FIX: Don't assume leading term is at index 0
                // Find the actual position of the leading monomial
                size_t leading_idx = 0;

                // The leading monomial should match element.mon (which comes from t_xw[i].mon)
                // t_xw[i].mon should be the leading monomial of this GB element
                const PP &expected_leading = element.mon;

                // Find which term in the GB polynomial matches the expected leading monomial
                bool found_leading = false;
                for (size_t j = 0; j < groebner_exponents[gb_index].size(); ++j) {
                    if (groebner_exponents[gb_index][j] == expected_leading) {
                        leading_idx = j;
                        found_leading = true;
                        break;
                    }
                }

                if (!found_leading) {
                    // This is a serious error - the border structure doesn't match the GB
                    std::cerr << "ERROR: Could not find expected leading monomial in GB element " << gb_index
                              << std::endl;
                    std::cerr << "  Expected leading: [";
                    for (size_t k = 0; k < expected_leading.size(); ++k) {
                        if (k > 0) std::cerr << ",";
                        std::cerr << expected_leading[k];
                    }
                    std::cerr << "]" << std::endl;
                    std::cerr << "  GB terms:" << std::endl;
                    for (size_t j = 0; j < groebner_exponents[gb_index].size(); ++j) {
                        std::cerr << "    Term " << j << ": [";
                        for (size_t k = 0; k < groebner_exponents[gb_index][j].size(); ++k) {
                            if (k > 0) std::cerr << ",";
                            std::cerr << groebner_exponents[gb_index][j][k];
                        }
                        std::cerr << "]" << std::endl;
                    }
                    // Continue with index 0 as fallback, but this will likely cause incorrect results
                    leading_idx = 0;
                }

                std::vector<PP> tail_exponents;
                std::vector<ModularCoeff> tail_coefficients;

                // CRITICAL: Julia skips the leading term when vectorizing GB polynomials
                // Their loop starts from index 2: "for i in 2:length(p_coeffs)"
                // This means for GB elements with only one term, nothing is added to the vector
                
                if (groebner_exponents[gb_index].size() == 1) {
                    // GB element with only leading term: maps to zero in quotient ring
                    // This is actually normal for some polynomial systems
                    // Leave tail empty - the vector will remain zeros
                    if (std::getenv("RUR_DEBUG_GB") && std::string(std::getenv("RUR_DEBUG_GB")) != "0") {
                        std::cout << "INFO: GB element " << gb_index << " has only leading term [";
                        for (size_t k = 0; k < expected_leading.size(); ++k) {
                            if (k > 0) std::cout << ",";
                            std::cout << expected_leading[k];
                        }
                        std::cout << "] = 0 in quotient ring" << std::endl;
                    }
                } else {
                    // Normal case: GB element has leading term and tail
                    // Skip the leading term (Julia starts from index 2)
                    // IMPORTANT: Scale tail by inverse of leading coefficient
                    ModularCoeff lt_coeff = groebner_coefficients[gb_index][leading_idx];
                    // Guard: if lt_coeff == 0 (should not happen), treat as 1
                    if (lt_coeff == 0) lt_coeff = 1;
                    ModularCoeff inv_lt = modular_inverse(lt_coeff, prime);
                    for (size_t j = 0; j < groebner_exponents[gb_index].size(); ++j) {
                        if (j != leading_idx) {
                            tail_exponents.push_back(groebner_exponents[gb_index][j]);
                            // Get the coefficient from the tail
                            ModularCoeff raw_coeff = groebner_coefficients[gb_index][j];

                            // For g(x) = c_lt*lt + sum c_j*m_j = 0, we need
                            // lt = - (sum (c_j / c_lt) * m_j) in the quotient ring.
                            // So add scaled, negated tail coefficients.
                            ModularCoeff neg = (prime - raw_coeff) % prime;
                            ModularCoeff neg_scaled = static_cast<ModularCoeff>(
                                (static_cast<AccModularCoeff>(neg) * inv_lt) % prime
                            );
                            tail_coefficients.push_back(neg_scaled);
                        }
                    }
                }

                // Reduce tail to standard monomials in quotient basis.
                // Some tail monomials may not be in the standard monomial set (non-reduced GB);
                // reduce them using GB relations lt(g_k) = -tail(g_k).
                t_v[i].assign(quotient_basis.size(), 0);

                auto debug_reduce =
                  (std::getenv("RUR_DEBUG_REDUCTION") && std::string(std::getenv("RUR_DEBUG_REDUCTION")) != "0");

                // Helper: find GB index for a leading monomial exactly equal to 'm'; returns -1 if none
                auto find_gb_by_leading_eq = [&](const PP &m) -> int {
                    // leading monomials are present in t_xw entries with var>0 && prev==0
                    // but we also can scan groebner_exponents to be safe
                    for (size_t gi = 0; gi < groebner_exponents.size(); ++gi) {
                        // Find actual leading term of GB[gi]
                        size_t lt_idx = 0;
                        for (size_t jj = 1; jj < groebner_exponents[gi].size(); ++jj) {
                            // We don't have monomial order here; rely on equality test via t_xw mapping where possible
                        }
                        // Fallback: try equality against any term tagged in border as leading for that GB index
                        // We can scan t_xw once for matching (var==gi+1 && prev==0)
                    }
                    // Cheaper exact check: scan t_xw for leading entries
                    for (const auto &sv : t_xw) {
                        if (sv.var > 0 && sv.prev == 0) {
                            if (sv.mon == m) return static_cast<int>(sv.var - 1);
                        }
                    }
                    return -1;
                };

                // Reduce a single (monomial, coeff) pair into quotient basis and accumulate into vec
                // Local recursive reducer (lambda with explicit type via auto)
                auto reduce_and_accum = [&](const PP &mon,
                                            ModularCoeff coeff,
                                            std::vector<ModularCoeff> &vec,
                                            int depth,
                                            const auto &self) -> void {
                    // Safety cap to avoid infinite loops on malformed GB
                    if (depth > 64) return;
                    auto it = std::find(quotient_basis.begin(), quotient_basis.end(), mon);
                    if (it != quotient_basis.end()) {
                        size_t idx = static_cast<size_t>(std::distance(quotient_basis.begin(), it));
                        vec[idx] = (vec[idx] + coeff) % prime;
                        return;
                    }
                    // Try exact match with a leading term
                    int gb_idx = find_gb_by_leading_eq(mon);
                    if (gb_idx >= 0 && gb_idx < static_cast<int>(groebner_exponents.size())) {
                        // Expand by the tail of GB[gb_idx]: mon = -sum c_j * tail_j
                        const auto &exps_k = groebner_exponents[gb_idx];
                        const auto &coeffs_k = groebner_coefficients[gb_idx];
                        if (exps_k.size() == 1) {
                            // Truly monomial GB element: treat as zero in this algebra
                            if (debug_reduce) {
                                std::cout << "[REDUCE] monomial ";
                                for (size_t z = 0; z < mon.size(); ++z) {
                                    if (z) std::cout << ",";
                                    std::cout << mon[z];
                                }
                                std::cout << " maps to 0 (monomial GB)" << std::endl;
                            }
                            return;
                        }
                        // Find the actual leading term position for this GB polynomial
                        size_t lt_idx_k = 0;
                        bool found_lt_k = false;
                        for (size_t jj = 0; jj < exps_k.size(); ++jj) {
                            if (exps_k[jj] == mon) { lt_idx_k = jj; found_lt_k = true; break; }
                        }
                        ModularCoeff lt_coeff_k = coeffs_k[lt_idx_k];
                        if (lt_coeff_k == 0) lt_coeff_k = 1;
                        ModularCoeff inv_lt_k = modular_inverse(lt_coeff_k, prime);

                        for (size_t tj = 0; tj < exps_k.size(); ++tj) {
                            // Skip the leading term; others form the tail
                            if (tj == lt_idx_k) continue;
                            ModularCoeff negc = (prime - coeffs_k[tj]) % prime;
                            // Scale by inverse leading coefficient and by current coeff
                            AccModularCoeff scaled = static_cast<AccModularCoeff>(negc) * inv_lt_k;
                            scaled %= prime;
                            scaled = (scaled * coeff) % prime;
                            ModularCoeff newc = static_cast<ModularCoeff>(scaled);
                            if (debug_reduce && depth < 3) {
                                std::cout << "[REDUCE] ";
                                for (size_t z = 0; z < mon.size(); ++z) {
                                    if (z) std::cout << ",";
                                    std::cout << mon[z];
                                }
                                std::cout << " -> add coeff " << (int64_t)newc << " * [";
                                for (size_t z = 0; z < exps_k[tj].size(); ++z) {
                                    if (z) std::cout << ",";
                                    std::cout << exps_k[tj][z];
                                }
                                std::cout << "]" << std::endl;
                            }
                            self(exps_k[tj], newc, vec, depth + 1, self);
                        }
                        return;
                    }
                    // If still not reduced, we cannot express this monomial with available data
                    if (debug_reduce) {
                        std::cout << "[REDUCE] WARNING: monomial not reducible and not in basis: [";
                        for (size_t z = 0; z < mon.size(); ++z) {
                            if (z) std::cout << ",";
                            std::cout << mon[z];
                        }
                        std::cout << "] (coeff=" << (int64_t)coeff << ")" << std::endl;
                    }
                };

                // Accumulate reductions for all tail terms
                for (size_t tj = 0; tj < tail_exponents.size(); ++tj) {
                    reduce_and_accum(tail_exponents[tj], tail_coefficients[tj], t_v[i], 0, reduce_and_accum);
                }

                // Debug (disabled)
                /*
                if (groebner_exponents[gb_index].size() <= 3) {
                    std::cout << "    Result t_v[" << i << "] = [";
                    for (size_t k = 0; k < t_v[i].size(); ++k) {
                        if (k > 0) std::cout << ",";
                        std::cout << t_v[i][k];
                    }
                    std::cout << "]" << std::endl;
                }
                */
            } else {
                // Error: invalid GB index
                t_v[i].clear();
            }
        } else if (element.var == 0 && element.prev > 0) {
            // This is a quotient element
            // prev contains the 1-based index into the quotient basis
            // Create a one-hot vector for this basis element
            // std::cout << "DEBUG: Processing quotient element for border element " << i
            //           << ", prev=" << element.prev << std::endl;
            t_v[i].resize(quotient_basis.size(), 0);
            if (element.prev > 0 && element.prev <= static_cast<int32_t>(quotient_basis.size())) {
                t_v[i][element.prev - 1] = 1; // 1-based to 0-based
                // std::cout << "    Created one-hot vector at position " << (element.prev - 1) << std::endl;
            }
        } else {
            // This is a border element (to be computed later by learn_compute_table)
            t_v[i].clear();
        }
    }
}

void
initialize_coefficient_vectors(std::vector<std::vector<ModularCoeff>> &t_v,
                               const std::vector<StackVect> &t_xw,
                               size_t quotient_basis_size) {
    t_v.clear();
    t_v.resize(t_xw.size());

    for (size_t i = 0; i < t_xw.size(); ++i) {
        if (t_xw[i].var > 0 && t_xw[i].prev == 0) {
            // GB element: For testing, just initialize with zeros
            // In full implementation, this would be vectorized GB polynomial
            t_v[i].resize(quotient_basis_size, 0);
            // For testing: set a simple pattern to make it identifiable
            if (t_xw[i].var <= static_cast<int32_t>(quotient_basis_size)) {
                t_v[i][t_xw[i].var - 1] = 1; // Simple unit vector pattern
            }
        } else if (t_xw[i].var == 0 && t_xw[i].prev > 0) {
            // Quotient element: single element vector [prev]
            t_v[i] = { static_cast<ModularCoeff>(t_xw[i].prev) };
        } else {
            // Border element: empty vector (to be computed by learn_compute_table)
            t_v[i].clear();
        }
    }
}

} // namespace julia_rur