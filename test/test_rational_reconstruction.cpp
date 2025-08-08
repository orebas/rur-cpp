#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <cassert>

using namespace julia_rur;

void test_simple_univariate() {
    std::cout << "Testing simple univariate system: x^2 - 2" << std::endl;
    
    std::vector<std::string> polynomials = {
        "1*x^2-2"  // x^2 - 2
    };
    
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = true;
    config.initial_prime_bits = 17;  // Use smaller primes for faster testing
    
    RationalRURResult result = compute_rational_rur(polynomials, variables, config);
    
    if (result.success) {
        std::cout << "\nSuccess! Minimal polynomial: ";
        for (size_t i = 0; i < result.minimal_polynomial.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << result.minimal_polynomial[i];
        }
        std::cout << std::endl;
        
        // Should be [-2, 0, 1] representing x^2 - 2
        assert(result.minimal_polynomial.size() == 3);
        assert(result.minimal_polynomial[0] == -2);
        assert(result.minimal_polynomial[1] == 0);
        assert(result.minimal_polynomial[2] == 1);
    } else {
        std::cerr << "Failed: " << result.error_message << std::endl;
    }
}

void test_linear_system() {
    std::cout << "\nTesting linear system: x - 3" << std::endl;
    
    std::vector<std::string> polynomials = {
        "1*x-3"  // x - 3
    };
    
    std::vector<std::string> variables = {"x"};
    
    RURConfig config;
    config.verbose = false;
    config.initial_prime_bits = 17;
    
    RationalRURResult result = compute_rational_rur(polynomials, variables, config);
    
    if (result.success) {
        std::cout << "Success! Minimal polynomial: ";
        for (size_t i = 0; i < result.minimal_polynomial.size(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << result.minimal_polynomial[i];
        }
        std::cout << std::endl;
        
        // Should be [-3, 1] representing x - 3
        assert(result.minimal_polynomial.size() == 2);
        assert(result.minimal_polynomial[0] == -3);
        assert(result.minimal_polynomial[1] == 1);
    } else {
        std::cerr << "Failed: " << result.error_message << std::endl;
    }
}

int main() {
    try {
        test_linear_system();
        test_simple_univariate();
        
        std::cout << "\nAll tests passed!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}