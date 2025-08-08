#include "../src/julia_rur/data_structures.hpp"
#include "../src/julia_rur/f4_polynomial_formatter.hpp"
extern "C" {
#include "../src/axf4_wrapper.h"
}
#include <iostream>
#include <vector>

using namespace julia_rur;

bool test_gb_for_prime(ModularCoeff prime) {
    // System: x^2 - 1 = 0, y - x = 0
    std::vector<std::string> polynomials = {
        "1*x^2-1",
        "1*y-1*x"
    };
    std::vector<std::string> variables = {"x", "y"};
    
    // Create F4 session
    std::vector<const char*> var_ptrs;
    for (const auto& var : variables) {
        var_ptrs.push_back(var.c_str());
    }
    
    axf4_session_t session = axf4_create_session(prime, var_ptrs.data(), variables.size());
    if (!session) {
        return false;
    }
    
    // Add polynomials
    for (const auto& poly : polynomials) {
        std::string formatted = format_polynomial_for_f4(poly);
        if (axf4_add_polynomial(session, formatted.c_str()) != 0) {
            axf4_destroy_session(session);
            return false;
        }
    }
    
    // Compute Gröbner basis
    axf4_result_t gb_result = axf4_compute_groebner_basis_keep_data(session);
    if (gb_result.status != 0) {
        axf4_free_result(&gb_result);
        axf4_destroy_session(session);
        return false;
    }
    
    // Check if the GB is correct (should contain y^2 - 1, represented as coefficient p-1)
    std::string gb_str(gb_result.groebner_basis);
    bool correct = gb_str.find(std::to_string(prime - 1)) != std::string::npos;
    
    // Clean up
    axf4_free_result(&gb_result);
    axf4_destroy_session(session);
    
    return correct;
}

int main() {
    std::cout << "Finding threshold where F4 bug occurs" << std::endl;
    std::cout << "Testing primes to see which produce correct GB" << std::endl;
    
    // Test various primes
    std::vector<ModularCoeff> test_primes = {
        // 30-bit primes
        1073741827, 1073741783, 1073741741, 1073741723,
        // 29-bit primes  
        536870909, 536870879, 536870869,
        // 28-bit primes
        268435459, 268435399, 268435367,
        // 27-bit primes
        134217727, 134217721, 134217689,
        // 26-bit primes
        67108859, 67108837, 67108819,
        // 25-bit primes
        33554431, 33554393, 33554383,
        // 24-bit primes
        16777213, 16777199, 16777183,
        // 23-bit primes
        8388607, 8388593, 8388587,
        // 22-bit primes
        4194301, 4194287, 4194277,
        // 21-bit primes
        2097143, 2097133, 2097131,
        // 20-bit primes
        1048573, 1048571, 1048559
    };
    
    ModularCoeff last_correct = 0;
    ModularCoeff first_wrong = 0;
    
    for (auto p : test_primes) {
        bool correct = test_gb_for_prime(p);
        std::cout << "Prime " << p << " (";
        
        // Show bit size
        int bits = 0;
        ModularCoeff temp = p;
        while (temp > 0) {
            bits++;
            temp >>= 1;
        }
        std::cout << bits << " bits): ";
        
        if (correct) {
            std::cout << "CORRECT (y^2 - 1)" << std::endl;
            if (p > last_correct) last_correct = p;
        } else {
            std::cout << "WRONG (y^2 + 1)" << std::endl;
            if (first_wrong == 0 || p < first_wrong) first_wrong = p;
        }
    }
    
    std::cout << "\nSummary:" << std::endl;
    if (last_correct > 0) {
        std::cout << "Largest prime with correct result: " << last_correct << std::endl;
    }
    if (first_wrong > 0) {
        std::cout << "Smallest prime with wrong result: " << first_wrong << std::endl;
    }
    
    return 0;
}