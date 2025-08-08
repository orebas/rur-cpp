#include "multiplication_tables.hpp"
#include <stdexcept>

namespace julia_rur {

DividesResult divides_with_var(const PP& a, const PP& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Monomial dimensions must match for division check");
    }
    
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

PP mul_pp_by_var(const PP& m, int32_t var_index) {
    PP result = m;  // Copy the monomial
    if (var_index >= 1 && var_index <= static_cast<int32_t>(m.size())) {
        result[var_index - 1] += 1;  // Convert to 0-based indexing
    }
    return result;
}

BorderSearchResult find_in_border(const PP& m, const std::vector<StackVect>& t_xw) {
    using namespace power_product;
    
    // TODO OPTIMIZATION: This function has O(b) complexity where b is border size
    // TODO OPTIMIZATION: Could use degree-indexed data structures for O(log b) lookup
    // TODO OPTIMIZATION: Or pre-sort border by degree and use binary search
    
    int32_t pos = static_cast<int32_t>(t_xw.size());
    int32_t res_flag = 0;
    int32_t res_dd = 0;  // Changed from PP to int32_t
    int32_t res_pos = 0;
    
    // Get total degree of target monomial
    uint32_t tm = total_degree(m);
    
    // Search backwards through border (by decreasing degrees to limit tests)
    while (pos > 0) {
        const PP& mm = t_xw[pos - 1].mon;  // Convert to 0-based indexing
        uint32_t tmm = total_degree(mm);
        
        if (tmm == tm) {
            // Same degree - check for exact match
            if (m == mm) {
                return BorderSearchResult(1, 0, pos);  // For exact match, var_index not used
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
                    res_dd = div_result.var_index;  // Now directly an int32_t
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

void prepare_table_mxi(
    const std::vector<PP>& ltg,                           // GB leading terms
    const std::vector<PP>& kb,                            // Quotient basis
    std::vector<StackVect>& t_xw,                         // Border structure (output)
    std::vector<std::vector<int32_t>>& i_xw               // Variable indices (output)
) {
    if (ltg.empty() || kb.empty()) {
        throw std::invalid_argument("Leading terms and quotient basis cannot be empty");
    }
    
    // TODO OPTIMIZATION: Replace linear searches with hash maps for O(1) lookup
    // TODO OPTIMIZATION: Use std::unordered_map<PP, int32_t> for kb and ltg lookups  
    // TODO OPTIMIZATION: Use std::unordered_map<PP, int32_t> for t_xw monomial lookups
    // TODO OPTIMIZATION: This would reduce complexity from O(n²·m) to O(n·m)
    
    // Get number of variables from first leading term
    int32_t nbv = static_cast<int32_t>(ltg[0].size());
    
    // Initialize tablex (i_xw) - one vector per variable, one entry per quotient basis element
    i_xw.clear();
    i_xw.resize(nbv);
    for (int32_t i = 0; i < nbv; ++i) {
        i_xw[i].resize(kb.size(), 0);
    }
    
    // Initialize border structure
    t_xw.clear();
    int32_t nb_stack = 0;
    
    // For each quotient basis element
    for (size_t j = 0; j < kb.size(); ++j) {
        PP m = kb[j];  // Copy the monomial
        
        // For each variable (Julia iterates in reverse order: nbv down to 1)
        for (int32_t ii = 1; ii <= nbv; ++ii) {
            int32_t i = nbv - ii + 1;  // Convert to Julia's reverse iteration
            
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
                    auto existing = std::find_if(t_xw.rbegin(), t_xw.rend(),
                        [&nm](const StackVect& sv) { return sv.mon == nm; });
                    
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
                auto existing = std::find_if(t_xw.rbegin(), t_xw.rend(),
                    [&nm](const StackVect& sv) { return sv.mon == nm; });
                
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

void learn_compute_table(
    std::vector<std::vector<ModularCoeff>>& t_v,          // Coefficient vectors (output)
    const std::vector<StackVect>& t_xw,                   // Border structure
    const std::vector<std::vector<int32_t>>& i_xw,        // Variable indices
    ModularCoeff prime                                    // Modular arithmetic prime
) {
    // t_v should already be initialized - don't clear it!
    // Just ensure it has the right size
    if (t_v.size() != t_xw.size()) {
        t_v.resize(t_xw.size());
    }
    
    // Julia algorithm: iteratively fill coefficient vectors
    int32_t nb = 1;  // Number of elements computed in this iteration
    std::vector<int32_t> t_learn;  // Track which elements were computed
    
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
                    !t_v[prev_idx - 1].empty()) {  // Convert to 0-based indexing
                    
                    // Get the variable to multiply by
                    int32_t var_idx = t_xw[i].var;
                    
                    // Multiply predecessor by variable to get this element
                    std::vector<ModularCoeff> result;
                    bool success = mul_var_quo_internal(
                        result,
                        t_v[prev_idx - 1],     // Predecessor coefficient vector
                        var_idx,               // Variable to multiply by
                        i_xw,                  // Variable indices
                        t_v,                   // Coefficient vectors
                        prime,                 // Modular prime
                        buf                    // Accumulator buffer
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
bool mul_var_quo_internal(
    std::vector<ModularCoeff>& result,                    // Output vector
    const std::vector<ModularCoeff>& input,              // Input vector in quotient basis
    int32_t var_index,                                    // Variable to multiply by (1-based)
    const std::vector<std::vector<int32_t>>& i_xw,       // Pre-computed indices
    const std::vector<std::vector<ModularCoeff>>& t_v,   // Coefficient vectors
    ModularCoeff prime,                                   // Modular arithmetic prime
    std::vector<AccModularCoeff>& buf                     // Accumulator buffer
) {
    size_t dim = input.size();
    result.resize(dim);
    
    // For optimization, we could pack reductions (Julia's pack_value)
    // For now, use simple approach: reduce every 16 operations
    const int32_t pack = 16;
    
    // Initialize accumulator buffer
    buf.resize(dim);
    for (size_t i = 0; i < dim; ++i) {
        buf[i] = 0;
    }
    
    bool continuer = true;
    
    // Convert var_index to 0-based for array access
    int32_t var_idx_0based = var_index - 1;
    
    // Check bounds
    if (var_idx_0based < 0 || var_idx_0based >= static_cast<int32_t>(i_xw.size())) {
        return false;
    }
    
    // For each element in input vector
    for (size_t j = 0; j < dim; ++j) {
        if (input[j] == 0) continue;  // Skip zero coefficients
        
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
        
        const auto& target_coeffs = t_v[target_idx];
        
        if (target_coeffs.size() > 1) {
            // Multiple coefficients: add_mul! operation
            // buf += input[j] * target_coeffs
            for (size_t k = 0; k < target_coeffs.size() && k < dim; ++k) {
                buf[k] += static_cast<AccModularCoeff>(input[j]) * 
                          static_cast<AccModularCoeff>(target_coeffs[k]);
            }
            
            // Periodic reduction to prevent overflow
            if ((j + 1) % pack == 0) {
                for (size_t k = 0; k < dim; ++k) {
                    buf[k] %= prime;
                }
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
    for (size_t i = 0; i < dim; ++i) {
        result[i] = static_cast<ModularCoeff>(buf[i] % prime);
    }
    
    return continuer;
}

void mul_var_quo(
    std::vector<ModularCoeff>& result,                    // Output vector
    const std::vector<ModularCoeff>& input,              // Input vector in quotient basis
    int32_t var_index,                                    // Variable to multiply by
    const std::vector<std::vector<int32_t>>& i_xw,       // Pre-computed indices
    const std::vector<std::vector<ModularCoeff>>& t_v,   // Coefficient vectors
    ModularCoeff prime                                    // Modular arithmetic prime
) {
    std::vector<AccModularCoeff> buf;
    mul_var_quo_internal(result, input, var_index, i_xw, t_v, prime, buf);
}

void vectorize_polynomial_in_quotient_basis(
    const std::vector<PP>& exponents,
    const std::vector<ModularCoeff>& coefficients,
    const std::vector<PP>& quotient_basis,
    std::vector<ModularCoeff>& result,
    ModularCoeff prime
) {
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
    
    // For each term in the polynomial
    for (size_t term_idx = 0; term_idx < exponents.size(); ++term_idx) {
        const PP& monomial = exponents[term_idx];
        ModularCoeff coeff = coefficients[term_idx];
        
        // Find this monomial in the quotient basis
        auto it = std::find(quotient_basis.begin(), quotient_basis.end(), monomial);
        if (it != quotient_basis.end()) {
            size_t basis_idx = std::distance(quotient_basis.begin(), it);
            // Add coefficient (modular arithmetic)
            result[basis_idx] = (result[basis_idx] + coeff) % prime;
        } else if (debug) {
            std::cout << "  Warning: monomial not in quotient basis" << std::endl;
        }
        // Note: If monomial not in quotient basis, it should be reducible by GB
        // This would indicate an error in the GB computation or quotient basis extraction
    }
}

void compute_fill_quo_gb(
    std::vector<std::vector<ModularCoeff>>& t_v,
    const std::vector<StackVect>& t_xw,
    const std::vector<std::vector<PP>>& groebner_exponents,
    const std::vector<std::vector<ModularCoeff>>& groebner_coefficients,
    const std::vector<PP>& quotient_basis,
    ModularCoeff prime
) {
    // Initialize output vectors
    t_v.clear();
    t_v.resize(t_xw.size());
    
    // Process each element in the border structure
    for (size_t i = 0; i < t_xw.size(); ++i) {
        const StackVect& element = t_xw[i];
        
        if (element.var > 0 && element.prev == 0) {
            // This is a Gröbner basis element
            // var contains the 1-based index into the Gröbner basis
            size_t gb_index = static_cast<size_t>(element.var - 1);  // Convert to 0-based
            
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
                
                size_t leading_idx = 0;
                std::vector<PP> tail_exponents;
                std::vector<ModularCoeff> tail_coefficients;
                
                for (size_t j = 0; j < groebner_exponents[gb_index].size(); ++j) {
                    if (j != leading_idx) {
                        tail_exponents.push_back(groebner_exponents[gb_index][j]);
                        // Get the coefficient from the GB
                        // The F4 library returns coefficients as ModularCoeff (uint32_t)
                        // but they might actually represent signed values
                        ModularCoeff raw_coeff = groebner_coefficients[gb_index][j];
                        
                        // If the raw coefficient is very large (> 2^31), it's likely a negative
                        // value stored as unsigned due to the F4 library's representation
                        int64_t signed_coeff;
                        if (raw_coeff > 0x80000000U) {
                            // Interpret as negative 32-bit signed integer
                            signed_coeff = static_cast<int64_t>(static_cast<int32_t>(raw_coeff));
                        } else if (raw_coeff > prime / 2) {
                            // This is a negative value mod prime
                            signed_coeff = static_cast<int64_t>(raw_coeff) - static_cast<int64_t>(prime);
                        } else {
                            signed_coeff = static_cast<int64_t>(raw_coeff);
                        }
                        
                        
                        // IMPORTANT: For a GB element g(x) = lt(g) + tail(g) = 0,
                        // we need lt(g) = -tail(g) in the quotient ring.
                        // So we must negate the tail coefficients.
                        // For example: y^2 - 1 = 0 means y^2 = 1, so we negate -1 to get 1
                        ModularCoeff negated = (prime - raw_coeff) % prime;
                        tail_coefficients.push_back(negated);
                        
                        // Debug (disabled)
                        /*
                        if (groebner_exponents[gb_index].size() <= 3) {
                            std::cout << "    Tail coeff[" << j << "] = " << raw_coeff;
                            if (raw_coeff > prime / 2) {
                                std::cout << " (represents " << (static_cast<int64_t>(raw_coeff) - prime) << ")";
                            }
                            std::cout << " -> negated to " << negated;
                            std::cout << std::endl;
                        }
                        */
                    }
                }
                
                // Vectorize the tail in the quotient basis
                t_v[i].resize(quotient_basis.size(), 0);
                // std::cout << "    Quotient basis size: " << quotient_basis.size() << ", t_v[" << i << "] size after resize: " << t_v[i].size() << std::endl;
                vectorize_polynomial_in_quotient_basis(
                    tail_exponents,
                    tail_coefficients,
                    quotient_basis,
                    t_v[i],
                    prime
                );
                
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

void initialize_coefficient_vectors(
    std::vector<std::vector<ModularCoeff>>& t_v,
    const std::vector<StackVect>& t_xw,
    size_t quotient_basis_size
) {
    t_v.clear();
    t_v.resize(t_xw.size());
    
    for (size_t i = 0; i < t_xw.size(); ++i) {
        if (t_xw[i].var > 0 && t_xw[i].prev == 0) {
            // GB element: For testing, just initialize with zeros
            // In full implementation, this would be vectorized GB polynomial
            t_v[i].resize(quotient_basis_size, 0);
            // For testing: set a simple pattern to make it identifiable
            if (t_xw[i].var <= static_cast<int32_t>(quotient_basis_size)) {
                t_v[i][t_xw[i].var - 1] = 1;  // Simple unit vector pattern
            }
        } else if (t_xw[i].var == 0 && t_xw[i].prev > 0) {
            // Quotient element: single element vector [prev]
            t_v[i] = {static_cast<ModularCoeff>(t_xw[i].prev)};
        } else {
            // Border element: empty vector (to be computed by learn_compute_table)
            t_v[i].clear();
        }
    }
}

} // namespace julia_rur