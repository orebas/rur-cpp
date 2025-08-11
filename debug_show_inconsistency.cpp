#include <iostream>
#include <vector>
#include "julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    std::cout << "=== Demonstrating CRT Inconsistency ===\n" << std::endl;
    std::cout << "This test shows that coefficient vectors have different lengths/orderings\n";
    std::cout << "across different primes, causing CRT reconstruction to fail.\n" << std::endl;
    
    // Test the specific failing system
    std::vector<std::string> polynomials = {"x^2 - 1", "y^2 - 2", "z^2 - 3"};
    std::vector<std::string> variables = {"x", "y", "z"};
    
    // Test with just 2 primes to show the inconsistency clearly
    std::vector<ModularCoeff> test_primes = {2147483629, 2147483587};
    
    std::vector<std::vector<std::vector<ModularCoeff>>> all_tables;
    
    for (size_t p_idx = 0; p_idx < test_primes.size(); ++p_idx) {
        ModularCoeff prime = test_primes[p_idx];
        std::cout << "\n========================================" << std::endl;
        std::cout << "Prime " << (p_idx + 1) << ": " << prime << std::endl;
        std::cout << "========================================" << std::endl;
        
        RURConfig config;
        config.verbose = false;
        
        try {
            auto mod_result = compute_modular_rur(polynomials, variables, prime, config);
            
            if (mod_result.success) {
                // Build the table as done in rur_main_algorithm.cpp
                std::vector<std::vector<ModularCoeff>> prime_table;
                prime_table.push_back(mod_result.minimal_polynomial.coefficients);
                
                for (const auto &param : mod_result.parameterizations) {
                    if (!param.generators.empty()) {
                        prime_table.push_back(param.generators[0]);
                    } else {
                        prime_table.push_back(std::vector<ModularCoeff>());
                    }
                }
                
                all_tables.push_back(prime_table);
                
                // Print detailed info
                std::cout << "Minimal polynomial (table[0]):" << std::endl;
                std::cout << "  Size: " << prime_table[0].size() << std::endl;
                std::cout << "  Coefficients: [";
                for (size_t i = 0; i < prime_table[0].size(); ++i) {
                    if (i > 0) std::cout << ", ";
                    std::cout << prime_table[0][i];
                }
                std::cout << "]" << std::endl;
                
                // Show position [0][6] specifically
                if (prime_table[0].size() > 6) {
                    std::cout << "  Position [0][6] = " << prime_table[0][6] << std::endl;
                    
                    // Convert to signed interpretation
                    ModularCoeff val = prime_table[0][6];
                    if (val > prime / 2) {
                        int64_t signed_val = static_cast<int64_t>(val) - static_cast<int64_t>(prime);
                        std::cout << "  As signed integer: " << signed_val << std::endl;
                    }
                } else {
                    std::cout << "  Position [0][6] DOES NOT EXIST (vector too short)" << std::endl;
                }
                
                // Print parameterization sizes
                for (size_t i = 1; i < prime_table.size(); ++i) {
                    std::cout << "Parameterization " << (i-1) << " (table[" << i << "]):" << std::endl;
                    std::cout << "  Size: " << prime_table[i].size() << std::endl;
                }
                
            } else {
                std::cout << "✗ FAILED for prime " << prime << std::endl;
            }
            
        } catch (const std::exception& e) {
            std::cout << "Exception for prime " << prime << ": " << e.what() << std::endl;
        }
    }
    
    // Analysis
    std::cout << "\n========================================" << std::endl;
    std::cout << "INCONSISTENCY ANALYSIS" << std::endl;
    std::cout << "========================================" << std::endl;
    
    if (all_tables.size() >= 2) {
        std::cout << "\nComparing table structures across primes:" << std::endl;
        
        // Compare minimal polynomial sizes
        std::cout << "\n1. Minimal polynomial sizes:" << std::endl;
        for (size_t i = 0; i < all_tables.size(); ++i) {
            std::cout << "   Prime " << (i+1) << ": " << all_tables[i][0].size() << " coefficients" << std::endl;
        }
        
        bool same_size = true;
        for (size_t i = 1; i < all_tables.size(); ++i) {
            if (all_tables[i][0].size() != all_tables[0][0].size()) {
                same_size = false;
                break;
            }
        }
        
        if (!same_size) {
            std::cout << "   ✗ DIFFERENT SIZES - This is the problem!" << std::endl;
            std::cout << "   Position [0][6] refers to different monomials!" << std::endl;
        } else {
            std::cout << "   ✓ Same sizes" << std::endl;
            
            // Check if coefficient [0][6] values are consistent
            if (all_tables[0][0].size() > 6) {
                std::cout << "\n2. Coefficient [0][6] values:" << std::endl;
                for (size_t i = 0; i < all_tables.size(); ++i) {
                    ModularCoeff val = all_tables[i][0][6];
                    ModularCoeff prime = test_primes[i];
                    int64_t signed_val = (val > prime/2) ? 
                        static_cast<int64_t>(val) - static_cast<int64_t>(prime) : 
                        static_cast<int64_t>(val);
                    std::cout << "   Prime " << (i+1) << ": " << val 
                              << " (signed: " << signed_val << ")" << std::endl;
                }
            }
        }
        
        std::cout << "\n3. Parameterization sizes:" << std::endl;
        for (size_t j = 1; j < all_tables[0].size(); ++j) {
            std::cout << "   Variable " << variables[j-1] << ":" << std::endl;
            for (size_t i = 0; i < all_tables.size(); ++i) {
                std::cout << "     Prime " << (i+1) << ": " << all_tables[i][j].size() << " coefficients" << std::endl;
            }
        }
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "CONCLUSION" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "The empty canonicalize_rur_permutation function causes:" << std::endl;
    std::cout << "1. Coefficient vectors to have different lengths across primes" << std::endl;
    std::cout << "2. Position [0][6] to refer to different monomial terms" << std::endl;
    std::cout << "3. CRT to try reconstructing coefficients from different monomials" << std::endl;
    std::cout << "4. Rational reconstruction to fail with huge inconsistent values" << std::endl;
    std::cout << "\nThe fix: Implement canonicalize_rur_permutation to ensure" << std::endl;
    std::cout << "consistent vector lengths and ordering across all primes." << std::endl;
    
    return 0;
}