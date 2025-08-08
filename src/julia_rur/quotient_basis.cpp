#include "quotient_basis.hpp"
#include <queue>
#include <unordered_set>
#include <algorithm>

namespace julia_rur {

// Hash function for PP to enable unordered_set usage
struct PPHash {
    std::size_t operator()(const PP& pp) const {
        std::size_t seed = 0;
        for (auto deg : pp) {
            seed ^= std::hash<int>{}(deg) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

int find_divisor(const PP& monomial, const std::vector<PP>& leading_terms) {
    // Check from end to beginning (Julia iterates backwards)
    for (int j = leading_terms.size() - 1; j >= 0; j--) {
        if (power_product::divides(leading_terms[j], monomial)) {
            return j + 1;  // Return 1-based index
        }
    }
    return 0;  // Not divisible by any leading term
}


bool is_zero_dimensional(const std::vector<PP>& leading_terms, int num_variables) {
    // For each variable, check if there's a leading term that's a pure power of it
    std::vector<bool> has_univariate(num_variables, false);
    
    for (const auto& lt : leading_terms) {
        int nonzero_count = 0;
        int nonzero_var = -1;
        
        for (int i = 0; i < num_variables; i++) {
            if (lt[i] > 0) {
                nonzero_count++;
                nonzero_var = i;
            }
        }
        
        // If this is a univariate term, mark that variable as covered
        if (nonzero_count == 1) {
            has_univariate[nonzero_var] = true;
        }
    }
    
    // Zero-dimensional if all variables have univariate terms
    for (bool covered : has_univariate) {
        if (!covered) {
            return false;
        }
    }
    return true;
}

std::vector<PP> compute_quotient_basis(const std::vector<PP>& leading_terms) {
    if (leading_terms.empty()) {
        return {};
    }
    
    // Check for GB = {1} (system has no solutions)
    if (leading_terms.size() == 1) {
        bool all_zero = true;
        for (auto deg : leading_terms[0]) {
            if (deg != 0) {
                all_zero = false;
                break;
            }
        }
        if (all_zero) {
            throw std::runtime_error("System has no solutions (Groebner basis = {1})");
        }
    }
    
    int num_variables = leading_terms[0].size();
    
    // Check if ideal is zero-dimensional
    if (!is_zero_dimensional(leading_terms, num_variables)) {
        throw std::domain_error("Input does not define zero-dimensional ideal");
    }
    
    // Quotient basis computation using BFS exploration
    std::vector<PP> quotient_basis;
    std::queue<PP> todo;
    std::unordered_set<PP, PPHash> inspected;
    
    // Start with monomial 1 (all degrees zero)
    PP one(num_variables);
    todo.push(one);
    inspected.insert(one);
    
    while (!todo.empty()) {
        PP current = todo.front();
        todo.pop();
        
        // Check if current monomial is divisible by any leading term
        if (find_divisor(current, leading_terms) == 0) {
            // Not divisible - add to quotient basis
            quotient_basis.push_back(current);
            
            // Generate all neighbors by incrementing each variable
            for (int i = 0; i < num_variables; i++) {
                PP neighbor = current;
                neighbor[i] += 1;
                
                // Only add if we haven't seen it before (O(1) lookup)
                if (inspected.find(neighbor) == inspected.end()) {
                    todo.push(neighbor);
                    inspected.insert(neighbor);
                }
            }
        }
    }
    
    // Sort quotient basis in degrevlex order (matching Julia)
    auto pp_isless_drl = [](const PP& a, const PP& b) {
        return power_product::compare_degrevlex(a, b) < 0;
    };
    std::sort(quotient_basis.begin(), quotient_basis.end(), pp_isless_drl);
    
    // CRITICAL VERIFICATION: Ensure the constant monomial '1' is first
    if (!quotient_basis.empty()) {
        bool first_is_constant = std::all_of(quotient_basis[0].begin(), quotient_basis[0].end(),
                                            [](int e){ return e == 0; });
        if (!first_is_constant) {
            // This should never happen with correct degrevlex ordering
            std::cerr << "CRITICAL ERROR in compute_quotient_basis: First element is not the constant '1'!" << std::endl;
            std::cerr << "  quotient_basis[0] = [";
            for (size_t i = 0; i < quotient_basis[0].size(); ++i) {
                if (i > 0) std::cerr << ",";
                std::cerr << quotient_basis[0][i];
            }
            std::cerr << "]" << std::endl;
            
            // Force fix it - find the constant and move it to the front
            auto it = std::find_if(quotient_basis.begin(), quotient_basis.end(),
                [](const PP& m) {
                    return std::all_of(m.begin(), m.end(), [](int e){ return e == 0; });
                });
            
            if (it != quotient_basis.end() && it != quotient_basis.begin()) {
                std::cerr << "  Found constant at position " << std::distance(quotient_basis.begin(), it) << std::endl;
                std::cerr << "  Moving it to the front..." << std::endl;
                std::iter_swap(quotient_basis.begin(), it);
            } else if (it == quotient_basis.end()) {
                std::cerr << "  ERROR: No constant monomial found in quotient basis at all!" << std::endl;
                // This is a serious error - the quotient basis should always contain '1'
            }
        }
    }
    
    return quotient_basis;
}

} // namespace julia_rur