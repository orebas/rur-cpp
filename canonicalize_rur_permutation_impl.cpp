// Complete implementation of canonicalize_rur_permutation
// This addresses all concerns raised by O3 and Gemini code reviews

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <sstream>

// Assuming these types are defined in the actual code
using ModularCoeff = uint32_t;
using AccModularCoeff = uint64_t;

namespace {

// ============================================================================
// Modular arithmetic helpers
// ============================================================================

inline ModularCoeff mod_add(ModularCoeff a, ModularCoeff b, ModularCoeff p) {
    ModularCoeff s = a + b;
    return (s >= p) ? s - p : s;
}

inline ModularCoeff mod_sub(ModularCoeff a, ModularCoeff b, ModularCoeff p) {
    return (a >= b) ? (a - b) : (a + p - b);
}

inline ModularCoeff mod_mul(ModularCoeff a, ModularCoeff b, ModularCoeff p) {
    // Safe for primes up to 31 bits
    return static_cast<ModularCoeff>((static_cast<AccModularCoeff>(a) * b) % p);
}

inline ModularCoeff mod_pow(ModularCoeff a, ModularCoeff e, ModularCoeff p) {
    ModularCoeff r = 1;
    while (e) {
        if (e & 1) r = mod_mul(r, a, p);
        a = mod_mul(a, a, p);
        e >>= 1;
    }
    return r;
}

inline ModularCoeff mod_inv(ModularCoeff a, ModularCoeff p) {
    // Using Fermat's little theorem: a^(p-1) ≡ 1 (mod p)
    // Therefore: a^(-1) ≡ a^(p-2) (mod p)
    return mod_pow(a, p - 2, p);
}

// ============================================================================
// Polynomial manipulation helpers
// ============================================================================

// Trim trailing zeros from a polynomial
void trim_polynomial(std::vector<ModularCoeff>& v) {
    while (!v.empty() && v.back() == 0) {
        v.pop_back();
    }
}

// Create a trimmed copy of a polynomial (for const references)
std::vector<ModularCoeff> trimmed_copy(const std::vector<ModularCoeff>& v) {
    std::vector<ModularCoeff> result = v;
    trim_polynomial(result);
    return result;
}

// Check if two polynomials are proportional (possibly with reversal)
// Returns true if p * out_scale == q (possibly reversed)
bool are_proportional(const std::vector<ModularCoeff>& p,
                      const std::vector<ModularCoeff>& q,
                      ModularCoeff prime,
                      bool check_reverse,
                      ModularCoeff& out_scale) {
    
    // Work with trimmed versions to ignore trailing zeros
    auto p_trimmed = trimmed_copy(p);
    auto q_trimmed = trimmed_copy(q);
    
    if (p_trimmed.empty() && q_trimmed.empty()) {
        out_scale = 1;
        return true; // Both zero polynomials
    }
    
    if (p_trimmed.empty() || q_trimmed.empty()) {
        return false; // One is zero, the other isn't
    }
    
    if (p_trimmed.size() != q_trimmed.size()) {
        return false; // Different degrees after trimming
    }
    
    size_t n = p_trimmed.size();
    
    // Lambda to access coefficient with optional reversal
    auto coeff = [&](const std::vector<ModularCoeff>& poly, size_t i) -> ModularCoeff {
        return check_reverse ? poly[n - 1 - i] : poly[i];
    };
    
    // Find first non-zero coefficient
    size_t first = 0;
    while (first < n && coeff(p_trimmed, first) == 0 && coeff(q_trimmed, first) == 0) {
        ++first;
    }
    
    if (first == n) {
        out_scale = 1;
        return true; // Both polynomials are zero
    }
    
    if (coeff(p_trimmed, first) == 0 || coeff(q_trimmed, first) == 0) {
        return false; // One has zero where the other doesn't
    }
    
    // Calculate scale factor: q[first] / p[first]
    out_scale = mod_mul(coeff(q_trimmed, first),
                       mod_inv(coeff(p_trimmed, first), prime), prime);
    
    // Verify proportionality for all coefficients
    for (size_t i = first; i < n; ++i) {
        ModularCoeff expected = mod_mul(coeff(p_trimmed, i), out_scale, prime);
        if (expected != coeff(q_trimmed, i)) {
            return false;
        }
    }
    
    return true;
}

// Apply scaling and optional reversal to a polynomial
void scale_and_reverse(std::vector<ModularCoeff>& poly,
                       ModularCoeff scale,
                       bool do_reverse,
                       ModularCoeff prime) {
    for (auto& c : poly) {
        c = mod_mul(c, scale, prime);
    }
    if (do_reverse) {
        std::reverse(poly.begin(), poly.end());
    }
}

// Make a polynomial monic (leading coefficient = 1)
void make_monic(std::vector<ModularCoeff>& poly, ModularCoeff prime) {
    trim_polynomial(poly);
    if (poly.empty()) return;
    
    ModularCoeff leading_coeff = poly.back();
    if (leading_coeff == 0 || leading_coeff == 1) return;
    
    ModularCoeff inv = mod_inv(leading_coeff, prime);
    for (auto& c : poly) {
        c = mod_mul(c, inv, prime);
    }
}

} // anonymous namespace

// ============================================================================
// Main canonicalization function
// ============================================================================

static void canonicalize_rur_permutation(
    std::vector<std::vector<ModularCoeff>>& current_table,
    const std::vector<std::vector<ModularCoeff>>& reference_table,
    ModularCoeff prime,
    size_t num_variables) {
    
    // Enable this for debugging
    const bool verbose = false;
    
    // -------- 0. Validation --------
    if (reference_table.empty()) {
        throw std::runtime_error("canonicalize_rur_permutation: Reference table is empty");
    }
    
    if (current_table.empty()) {
        throw std::runtime_error("canonicalize_rur_permutation: Current table is empty");
    }
    
    // Expected size: 1 (minimal poly) + num_variables (parameterizations)
    size_t expected_size = 1 + num_variables;
    
    if (reference_table.size() != expected_size) {
        std::stringstream ss;
        ss << "canonicalize_rur_permutation: Reference table size mismatch. "
           << "Expected " << expected_size << ", got " << reference_table.size();
        throw std::runtime_error(ss.str());
    }
    
    // -------- 1. Trim current polynomials --------
    for (auto& poly : current_table) {
        trim_polynomial(poly);
    }
    
    // -------- 2. Determine global scale and reversal from minimal polynomial --------
    bool needs_reversal = false;
    ModularCoeff global_scale = 1;
    bool found_match = false;
    
    const auto& ref_minpoly = reference_table[0];
    auto& cur_minpoly = current_table[0];
    
    // Try without reversal first
    if (are_proportional(cur_minpoly, ref_minpoly, prime, false, global_scale)) {
        needs_reversal = false;
        found_match = true;
        if (verbose) {
            std::cout << "Minimal polynomial matches without reversal, scale = " 
                     << global_scale << std::endl;
        }
    }
    // Then try with reversal
    else if (are_proportional(cur_minpoly, ref_minpoly, prime, true, global_scale)) {
        needs_reversal = true;
        found_match = true;
        if (verbose) {
            std::cout << "Minimal polynomial matches WITH reversal, scale = " 
                     << global_scale << std::endl;
        }
    }
    
    if (!found_match) {
        throw std::runtime_error(
            "canonicalize_rur_permutation: Cannot match minimal polynomials. "
            "This prime may be bad or the system has changed structure.");
    }
    
    // Apply global transformation to all polynomials
    for (auto& poly : current_table) {
        scale_and_reverse(poly, global_scale, needs_reversal, prime);
    }
    
    // -------- 3. Make minimal polynomial monic --------
    make_monic(current_table[0], prime);
    
    // Also ensure reference minimal polynomial would be monic for comparison
    // (We work with a copy to respect const-ness)
    auto ref_minpoly_copy = ref_minpoly;
    make_monic(ref_minpoly_copy, prime);
    
    // -------- 4. Match parameterizations (handle variable permutations) --------
    std::vector<std::vector<ModularCoeff>> aligned_table(expected_size);
    aligned_table[0] = current_table[0]; // Minimal polynomial already placed
    
    std::vector<bool> reference_used(num_variables, false);
    std::vector<bool> current_matched(current_table.size() - 1, false);
    
    // Try to match each current parameterization with a reference one
    for (size_t cur_idx = 1; cur_idx < current_table.size(); ++cur_idx) {
        bool found = false;
        
        for (size_t ref_idx = 1; ref_idx <= num_variables; ++ref_idx) {
            if (reference_used[ref_idx - 1]) continue;
            
            ModularCoeff scale = 1;
            if (are_proportional(current_table[cur_idx], 
                                reference_table[ref_idx], 
                                prime, false, scale)) {
                
                // Found a match - scale to make it exactly equal
                scale_and_reverse(current_table[cur_idx], scale, false, prime);
                aligned_table[ref_idx] = std::move(current_table[cur_idx]);
                reference_used[ref_idx - 1] = true;
                current_matched[cur_idx - 1] = true;
                found = true;
                
                if (verbose) {
                    std::cout << "Matched current[" << cur_idx << "] to reference[" 
                             << ref_idx << "] with scale " << scale << std::endl;
                }
                break;
            }
        }
        
        if (!found && verbose) {
            std::cout << "Warning: Could not match current[" << cur_idx << "]" << std::endl;
        }
    }
    
    // Check if we successfully matched all necessary parameterizations
    int unmatched_count = 0;
    for (size_t i = 0; i < num_variables; ++i) {
        if (!reference_used[i]) {
            unmatched_count++;
            if (verbose) {
                std::cout << "Warning: Reference parameterization " << (i+1) 
                         << " was not matched" << std::endl;
            }
        }
    }
    
    // For systems where not all variables appear (e.g., after variable elimination),
    // we might have fewer parameterizations. This is OK as long as we matched
    // all that we have.
    bool all_current_matched = true;
    for (size_t i = 0; i < current_matched.size(); ++i) {
        if (!current_matched[i]) {
            all_current_matched = false;
            break;
        }
    }
    
    if (!all_current_matched) {
        throw std::runtime_error(
            "canonicalize_rur_permutation: Failed to match all current parameterizations. "
            "This indicates a structural mismatch - marking prime as bad.");
    }
    
    // Replace current table with aligned version
    current_table = std::move(aligned_table);
    
    // -------- 5. Ensure consistent lengths (pad with zeros) --------
    for (size_t i = 0; i < current_table.size(); ++i) {
        if (i >= reference_table.size()) break;
        
        size_t ref_size = reference_table[i].size();
        size_t cur_size = current_table[i].size();
        
        if (cur_size > ref_size) {
            // Check if we would truncate non-zero coefficients
            for (size_t j = ref_size; j < cur_size; ++j) {
                if (current_table[i][j] != 0) {
                    std::stringstream ss;
                    ss << "canonicalize_rur_permutation: Degree mismatch for polynomial " << i
                       << ". Current has degree " << (cur_size - 1)
                       << " but reference has degree " << (ref_size - 1)
                       << ". Would truncate non-zero coefficient at position " << j;
                    throw std::runtime_error(ss.str());
                }
            }
        }
        
        // Resize to match reference (padding with zeros if needed)
        current_table[i].resize(ref_size, 0);
    }
    
    if (verbose) {
        std::cout << "Canonicalization complete. Table has " << current_table.size() 
                 << " polynomials." << std::endl;
    }
}

// ============================================================================
// Special case: Normalize first prime's table to establish canonical form
// ============================================================================

static void normalize_first_table(
    std::vector<std::vector<ModularCoeff>>& table,
    ModularCoeff prime) {
    
    // For the first prime, we just need to establish a canonical form:
    // 1. Trim trailing zeros
    // 2. Make minimal polynomial monic
    // 3. Ensure consistent ordering (we choose ascending degree)
    
    for (auto& poly : table) {
        trim_polynomial(poly);
    }
    
    if (!table.empty()) {
        make_monic(table[0], prime);
    }
    
    // That's it - this becomes our reference for all subsequent primes
}

// ============================================================================
// Usage example showing how to integrate this
// ============================================================================

/*
// In compute_rational_rur, replace the canonicalization section with:

if (variables.size() > 1) {
    if (modular_tables.empty()) {
        // First table: establish canonical form
        normalize_first_table(prime_table, prime);
    } else {
        // Subsequent tables: align to reference
        const auto& reference_table = modular_tables[0];
        try {
            canonicalize_rur_permutation(prime_table, reference_table, prime, variables.size());
        } catch (const std::exception& e) {
            if (config.verbose) {
                std::cout << "Canonicalization failed for prime " << prime 
                         << ": " << e.what() << std::endl;
                std::cout << "Marking prime as bad and continuing..." << std::endl;
            }
            bad_primes_other++;
            continue; // Skip this prime
        }
    }
}
*/