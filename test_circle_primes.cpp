#include <iostream>
#include <vector>
#include <string>
#include "julia_rur/hyperplane_sections.hpp"
#include "julia_rur/prime_utils.hpp"

int main() {
    std::vector<std::string> polynomials = { "x^2 + y^2 - 1" };
    std::vector<std::string> variables = { "x", "y" };
    
    std::cout << "Testing x^2 + y^2 - 1 over various primes:\n\n";
    
    // Test specific small primes
    std::vector<julia_rur::ModularCoeff> test_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
        101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199
    };
    
    for (auto prime : test_primes) {
        auto result = julia_rur::analyze_system_dimension(polynomials, variables, prime);
        std::cout << "Prime " << prime << ": ";
        std::cout << "dim=" << result.dimension;
        std::cout << ", method=" << result.method_used;
        std::cout << ", is_zero_dim=" << result.is_zero_dimensional;
        std::cout << std::endl;
    }
    
    std::cout << "\n\nNow testing with large random primes:\n";
    
    // Test with random large primes
    for (int i = 0; i < 10; i++) {
        auto prime = julia_rur::generate_random_prime(28, 30);
        auto result = julia_rur::analyze_system_dimension(polynomials, variables, prime);
        std::cout << "Prime " << prime << ": ";
        std::cout << "dim=" << result.dimension;
        std::cout << ", method=" << result.method_used;
        std::cout << ", is_zero_dim=" << result.is_zero_dimensional;
        std::cout << std::endl;
    }
    
    return 0;
}