#include <iostream>
#include "julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    std::cout << "=== Testing Minimal Polynomial Consistency Across Primes ===" << std::endl;
    
    // Test the specific failing system
    std::vector<std::string> polynomials = {"x^2 - 1", "y^2 - 2", "z^2 - 3"};
    std::vector<std::string> variables = {"x", "y", "z"};
    
    RURConfig config;
    config.verbose = true;
    config.initial_prime_bits = 28;  // Use moderate size primes
    config.max_prime_bits = 28;      // Keep them small for debugging
    config.num_threads = 1;          // Single threaded for clarity
    
    std::cout << "\nTesting system: x^2-1, y^2-2, z^2-3" << std::endl;
    std::cout << "Will show minimal polynomial computation details for first few primes..." << std::endl;
    
    try {
        auto result = compute_rational_rur(polynomials, variables, config);
        
        if (result.success) {
            std::cout << "\n✓ SUCCESS: Found RUR solution" << std::endl;
            std::cout << "Minimal polynomial found" << std::endl;
        } else {
            std::cout << "\n✗ EXPECTED: CRT reconstruction failed" << std::endl;
            std::cout << "Error: " << result.error_message << std::endl;
            std::cout << "\nThis should show debugging information about coefficient inconsistency" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
    }
    
    std::cout << "\n=== Test Complete ===" << std::endl;
    std::cout << "Check the output above for CRT DEBUG messages showing coefficient values" << std::endl;
    std::cout << "Look for inconsistent remainders across different primes for the same coefficient" << std::endl;
    
    return 0;
}