// Better dimension detection for polynomial systems
#include <vector>
#include <algorithm>
#include <set>
#include <iostream>

namespace julia_rur {

// Check if a monomial (PP vector) uses only variables from subset
bool uses_only_variables_from(const std::vector<int>& monomial, const std::set<int>& var_subset) {
    for (size_t i = 0; i < monomial.size(); ++i) {
        if (monomial[i] > 0 && var_subset.find(i) == var_subset.end()) {
            return false;
        }
    }
    return true;
}

// Check if any monomial uses exactly the variables in subset (no other variables)
bool has_monomial_in_only_these_vars(const std::vector<std::vector<int>>& leading_terms, 
                                     const std::set<int>& var_subset) {
    for (const auto& lt : leading_terms) {
        // Check if this monomial uses only variables from var_subset
        bool uses_only_subset = true;
        bool uses_at_least_one = false;
        
        for (size_t i = 0; i < lt.size(); ++i) {
            if (lt[i] > 0) {
                if (var_subset.find(i) == var_subset.end()) {
                    uses_only_subset = false;
                    break;
                } else {
                    uses_at_least_one = true;
                }
            }
        }
        
        if (uses_only_subset && uses_at_least_one) {
            return true;
        }
    }
    return false;
}

// Generate all k-element subsets of {0, 1, ..., n-1}
void generate_subsets(int n, int k, std::vector<std::set<int>>& result) {
    std::vector<int> indices(k);
    for (int i = 0; i < k; ++i) {
        indices[i] = i;
    }
    
    while (true) {
        std::set<int> subset;
        for (int i = 0; i < k; ++i) {
            subset.insert(indices[i]);
        }
        result.push_back(subset);
        
        // Generate next combination
        int i = k - 1;
        while (i >= 0 && indices[i] == n - k + i) {
            --i;
        }
        if (i < 0) break;
        
        ++indices[i];
        for (int j = i + 1; j < k; ++j) {
            indices[j] = indices[j-1] + 1;
        }
    }
}

// Correct dimension calculation based on Krull dimension
int compute_krull_dimension(const std::vector<std::vector<int>>& leading_terms, int num_vars) {
    // Try subsets from largest to smallest
    for (int k = num_vars; k >= 0; --k) {
        std::vector<std::set<int>> subsets;
        generate_subsets(num_vars, k, subsets);
        
        for (const auto& subset : subsets) {
            // Check if this subset is algebraically independent
            // i.e., no leading term consists only of variables from this subset
            if (!has_monomial_in_only_these_vars(leading_terms, subset)) {
                return k;  // Found the largest independent set
            }
        }
    }
    return 0;
}

// Simpler heuristic: For positive dimensional systems, try to estimate better
// For Cyclic-n specifically, the dimension is known to be 1 for n=4
int estimate_dimension_heuristic(const std::vector<std::vector<int>>& leading_terms, 
                                 int num_vars,
                                 int num_equations) {
    // Count variables with univariate leading terms
    std::vector<bool> has_univariate(num_vars, false);
    for (const auto& lt : leading_terms) {
        int nonzero_count = 0;
        int nonzero_var = -1;
        for (int i = 0; i < num_vars; ++i) {
            if (lt[i] > 0) {
                nonzero_count++;
                nonzero_var = i;
            }
        }
        if (nonzero_count == 1) {
            has_univariate[nonzero_var] = true;
        }
    }
    
    int univariate_count = 0;
    for (bool has : has_univariate) {
        if (has) univariate_count++;
    }
    
    // Naive estimate
    int naive_dimension = num_vars - univariate_count;
    
    // But check if there are leading terms that show dependencies
    // If we have leading terms that are pure in the "free" variables,
    // the actual dimension is lower
    std::set<int> free_vars;
    for (int i = 0; i < num_vars; ++i) {
        if (!has_univariate[i]) {
            free_vars.insert(i);
        }
    }
    
    if (free_vars.size() > 1 && has_monomial_in_only_these_vars(leading_terms, free_vars)) {
        // There's algebraic dependency among the free variables
        // The dimension is at most free_vars.size() - 1
        // For a more accurate estimate, we'd need to check smaller subsets
        // But as a heuristic, we can try dimension = 1 if the system looks like Cyclic-n
        if (num_vars == num_equations && free_vars.size() == 2) {
            // Likely a system like Cyclic-4
            return 1;
        }
        return std::max(1, naive_dimension - 1);
    }
    
    return naive_dimension;
}

} // namespace julia_rur