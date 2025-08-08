/**
 * Enhanced separating element search with coefficient tracking
 * 
 * This version tracks the successful separating element coefficients
 * so they can be reused for subsequent primes.
 */

#include "multiplication_tables.hpp"
#include "univariate_parameterization.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <set>
#include <vector>

namespace julia_rur {

// Structure to track successful separating element
struct SeparatingElementData {
    bool found = false;
    bool is_single_variable = false;
    int variable_index = -1;  // 1-based if is_single_variable
    std::vector<int> coefficients;  // coefficients for linear form if !is_single_variable
};

/**
 * Enhanced version that returns the separating element data
 */
std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>, SeparatingElementData>
compute_univariate_parameterization_with_tracking(
    std::vector<PP> quotient_basis,  // Pass by value to allow modification
    const std::vector<std::vector<int32_t>>& i_xw,
    std::vector<std::vector<ModularCoeff>>& t_v,
    int num_variables,
    ModularCoeff prime,
    const SeparatingElementData* hint,  // Optional hint from previous computation
    int max_retries
) {
    // Ensure '1' is first in quotient basis
    auto it = std::find_if(quotient_basis.begin(), quotient_basis.end(), 
        [](const PP& m) {
            return std::all_of(m.begin(), m.end(), [](int e){ return e == 0; });
        });

    if (it != quotient_basis.end() && it != quotient_basis.begin()) {
        std::cout << "INFO: Reordering quotient basis to place '1' at the front." << std::endl;
        std::iter_swap(quotient_basis.begin(), it);
    }
    
    MinimalPolynomialResult min_poly;
    std::vector<BivariateResult> parameterizations(num_variables);
    SeparatingElementData result_data;
    
    // Try the hint first if provided
    if (hint && hint->found) {
        if (hint->is_single_variable && hint->variable_index > 0) {
            std::cout << "\\nTrying hint: variable x" << hint->variable_index 
                      << " as separating element for prime " << prime << "..." << std::endl;
            
            auto [success, mp, params] = try_separating_element(
                quotient_basis, i_xw, t_v, num_variables, prime,
                std::vector<ModularCoeff>(), hint->variable_index
            );
            
            if (success) {
                std::cout << "SUCCESS: Hint worked! Variable x" << hint->variable_index 
                          << " is separating for prime " << prime << std::endl;
                result_data.found = true;
                result_data.is_single_variable = true;
                result_data.variable_index = hint->variable_index;
                return { true, mp, params, result_data };
            }
        } else if (!hint->is_single_variable && !hint->coefficients.empty()) {
            std::cout << "\\nTrying hint: linear form [";
            for (size_t i = 0; i < hint->coefficients.size(); ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << hint->coefficients[i];
            }
            std::cout << "] as separating element for prime " << prime << "..." << std::endl;
            
            // Build linear form from hint coefficients
            std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);
            for (int i = 0; i < num_variables && i < static_cast<int>(hint->coefficients.size()); ++i) {
                if (hint->coefficients[i] == 0) continue;
                
                std::vector<ModularCoeff> var_vec = element_to_vector(i + 1, i_xw, t_v, quotient_basis.size());
                ModularCoeff lambda_mod = (hint->coefficients[i] < 0) 
                    ? (prime - ((-hint->coefficients[i]) % prime)) % prime 
                    : hint->coefficients[i] % prime;
                
                for (size_t j = 0; j < linear_form.size(); ++j) {
                    AccModularCoeff prod = static_cast<AccModularCoeff>(lambda_mod) * var_vec[j];
                    linear_form[j] = (linear_form[j] + prod % prime) % prime;
                }
            }
            
            // Check if it's separating
            MinimalPolynomialResult test_mp = compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);
            
            if (test_mp.success && test_mp.degree == quotient_basis.size()) {
                std::cout << "SUCCESS: Hint worked! Linear form is separating for prime " << prime << std::endl;
                
                auto [success, mp, params] = try_separating_element(
                    quotient_basis, i_xw, t_v, num_variables, prime, linear_form, -1
                );
                
                if (success) {
                    result_data.found = true;
                    result_data.is_single_variable = false;
                    result_data.coefficients = hint->coefficients;
                    return { true, mp, params, result_data };
                }
            } else {
                std::cout << "Hint didn't work for prime " << prime << " (degree " << test_mp.degree 
                          << " != " << quotient_basis.size() << "), trying other methods..." << std::endl;
            }
        }
    }
    
    // For single variable systems, just use the variable itself
    if (num_variables == 1) {
        std::cout << "Single variable system: using x1 as separating element" << std::endl;
        auto [success, mp, params] = try_separating_element(
            quotient_basis, i_xw, t_v, num_variables, prime,
            std::vector<ModularCoeff>(), 1
        );

        if (success) {
            std::cout << "SUCCESS: Single variable parameterization!" << std::endl;
            result_data.found = true;
            result_data.is_single_variable = true;
            result_data.variable_index = 1;
            return { true, mp, params, result_data };
        } else {
            std::cerr << "Failed to parameterize single variable system!" << std::endl;
            return { false, min_poly, parameterizations, result_data };
        }
    }
    
    // Try individual variables (in reverse order as per Julia)
    for (int var = num_variables; var >= 1; --var) {
        std::cout << "\\nTrying variable x" << var << " as separating element..." << std::endl;
        auto [success, mp, params] = try_separating_element(
            quotient_basis, i_xw, t_v, num_variables, prime,
            std::vector<ModularCoeff>(), var
        );

        if (success) {
            std::cout << "SUCCESS: Variable x" << var << " is a separating element!" << std::endl;
            result_data.found = true;
            result_data.is_single_variable = true;
            result_data.variable_index = var;
            return { true, mp, params, result_data };
        }
    }
    
    // Pre-calculate variable vectors
    std::vector<std::vector<ModularCoeff>> var_vectors(num_variables);
    for (int i = 0; i < num_variables; ++i) {
        if (i < static_cast<int>(i_xw.size()) && !i_xw[i].empty()) {
            int32_t table_idx = i_xw[i][0];
            if (table_idx > 0 && table_idx <= static_cast<int32_t>(t_v.size())) {
                var_vectors[i] = t_v[table_idx - 1];
            } else {
                var_vectors[i].assign(quotient_basis.size(), 0);
            }
        } else {
            var_vectors[i].assign(quotient_basis.size(), 0);
        }
    }
    
    // Try structured linear combinations
    std::cout << "\\nTrying structured linear combinations...\\n";
    
    if (num_variables == 2) {
        // Try systematic small coefficients for 2 variables
        std::vector<std::pair<int, int>> structured_coeffs = {
            {1, 1}, {1, 2}, {2, 1}, {1, -1}, {2, -1}, {-1, 2}, 
            {3, 1}, {1, 3}, {3, -1}, {-1, 3}, {1, 0}, {0, 1}
        };
        
        // Add more systematic combinations
        for (int a = -5; a <= 5; ++a) {
            if (a == 0) continue;
            for (int b = -5; b <= 5; ++b) {
                if (b == 0 && std::abs(a) <= 1) continue;
                structured_coeffs.push_back({a, b});
            }
        }
        
        for (const auto& [c1, c2] : structured_coeffs) {
            std::cout << "Trying " << c1 << "*x1 + " << c2 << "*x2..." << std::endl;
            
            std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);
            
            // Build linear form
            const auto& vec1 = var_vectors[0];
            ModularCoeff c1_mod = (c1 < 0) ? (prime - ((-c1) % prime)) % prime : c1 % prime;
            for (size_t j = 0; j < linear_form.size(); ++j) {
                linear_form[j] = (linear_form[j] + static_cast<AccModularCoeff>(c1_mod) * vec1[j]) % prime;
            }
            
            const auto& vec2 = var_vectors[1];
            ModularCoeff c2_mod = (c2 < 0) ? (prime - ((-c2) % prime)) % prime : c2 % prime;
            for (size_t j = 0; j < linear_form.size(); ++j) {
                linear_form[j] = (linear_form[j] + static_cast<AccModularCoeff>(c2_mod) * vec2[j]) % prime;
            }
            
            // Check if it's separating
            MinimalPolynomialResult test_mp = compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);
            
            if (test_mp.success && test_mp.degree == quotient_basis.size()) {
                std::cout << "  -> Minimal polynomial degree " << test_mp.degree 
                          << " matches quotient size!" << std::endl;
                
                auto [success, mp, params] = try_separating_element(
                    quotient_basis, i_xw, t_v, num_variables, prime, linear_form, -1
                );
                
                if (success) {
                    std::cout << "SUCCESS: Found separating linear combination " 
                              << c1 << "*x1 + " << c2 << "*x2!" << std::endl;
                    result_data.found = true;
                    result_data.is_single_variable = false;
                    result_data.coefficients = {c1, c2};
                    return { true, mp, params, result_data };
                }
            }
        }
    }
    
    // Try Julia-style systematic search
    std::cout << "\\nTrying systematic separating element search (Julia-style)...\\n";
    
    std::vector<int> sep_coeffs(num_variables, 0);
    if (num_variables >= 2) {
        sep_coeffs[num_variables - 1] = 1;   // Last variable
        sep_coeffs[num_variables - 2] = -1;  // Second to last
        
        int max_attempts = 100;
        
        for (int attempt = 0; attempt < max_attempts; ++attempt) {
            std::cout << "  Trying separating form: [";
            for (int i = 0; i < num_variables; ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << sep_coeffs[i];
            }
            std::cout << "]" << std::endl;
            
            // Build linear form
            std::vector<ModularCoeff> linear_form(quotient_basis.size(), 0);
            for (int i = 0; i < num_variables; ++i) {
                if (sep_coeffs[i] == 0) continue;
                
                std::vector<ModularCoeff> var_vec = element_to_vector(i + 1, i_xw, t_v, quotient_basis.size());
                ModularCoeff lambda_mod = (sep_coeffs[i] < 0) 
                    ? (prime - ((-sep_coeffs[i]) % prime)) % prime 
                    : sep_coeffs[i] % prime;
                
                for (size_t j = 0; j < linear_form.size(); ++j) {
                    AccModularCoeff prod = static_cast<AccModularCoeff>(lambda_mod) * var_vec[j];
                    linear_form[j] = (linear_form[j] + prod % prime) % prime;
                }
            }
            
            MinimalPolynomialResult test_mp = compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);
            
            if (test_mp.success && test_mp.degree == quotient_basis.size()) {
                std::cout << "  SUCCESS: Found separating element with coeffs [";
                for (int i = 0; i < num_variables; ++i) {
                    if (i > 0) std::cout << ", ";
                    std::cout << sep_coeffs[i];
                }
                std::cout << "]" << std::endl;
                
                auto [success, mp, params] = try_separating_element(
                    quotient_basis, i_xw, t_v, num_variables, prime, linear_form, -1
                );
                
                if (success) {
                    result_data.found = true;
                    result_data.is_single_variable = false;
                    result_data.coefficients = sep_coeffs;
                    return { true, mp, params, result_data };
                }
            }
            
            // Update coefficients (Julia-style heuristic)
            int uu = test_mp.degree - 1;
            if (uu < 0) uu = 0;
            if (uu >= num_variables) uu = num_variables - 1;
            
            if (sep_coeffs[uu] < 0) {
                sep_coeffs[uu] = -sep_coeffs[uu];
            } else {
                sep_coeffs[uu] = -sep_coeffs[uu] - 1;
            }
        }
    }
    
    // Fall back to random search
    std::cout << "\\nTrying random linear combinations...\\n";
    std::random_device rd;
    std::mt19937 gen(rd());
    
    for (int attempt = 0; attempt < max_retries; ++attempt) {
        std::cout << "\\nAttempt " << (attempt + 1) << ": Trying random linear form..." << std::endl;
        
        std::vector<ModularCoeff> linear_form = compute_random_linear_form(
            quotient_basis, i_xw, t_v, num_variables, prime, gen
        );
        
        MinimalPolynomialResult test_mp = compute_minimal_polynomial_flint(linear_form, i_xw, t_v, quotient_basis, prime);
        
        if (test_mp.success && test_mp.degree == quotient_basis.size()) {
            auto [success, mp, params] = try_separating_element(
                quotient_basis, i_xw, t_v, num_variables, prime, linear_form, -1
            );
            
            if (success) {
                std::cout << "SUCCESS: Found valid separating linear form!" << std::endl;
                // For random, we don't store coefficients as they may be too complex
                result_data.found = true;
                result_data.is_single_variable = false;
                // Leave coefficients empty for random
                return { true, mp, params, result_data };
            }
        }
    }
    
    std::cerr << "Failed to find separating element after " << max_retries << " attempts" << std::endl;
    return { false, min_poly, parameterizations, result_data };
}

} // namespace julia_rur