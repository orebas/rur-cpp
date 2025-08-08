#include "../src/julia_rur/bivariate_algorithm.hpp"
#include "../src/julia_rur/multiplication_tables.hpp"
#include "../src/julia_rur/polynomial_operations.hpp"
#include <iostream>
#include <cassert>
#include <map>

using namespace julia_rur;

void test_modular_inverse() {
    std::cout << "Testing modular_inverse function..." << std::endl;
    
    ModularCoeff prime = 100003;
    
    // Test basic cases
    ModularCoeff inv2 = modular_inverse(2, prime);
    assert((2 * static_cast<uint64_t>(inv2)) % prime == 1);
    std::cout << "  ✓ inv(2) mod " << prime << " = " << inv2 << std::endl;
    
    ModularCoeff inv7 = modular_inverse(7, prime);
    assert((7 * static_cast<uint64_t>(inv7)) % prime == 1);
    std::cout << "  ✓ inv(7) mod " << prime << " = " << inv7 << std::endl;
    
    // Test edge case
    ModularCoeff inv_max = modular_inverse(prime - 1, prime);
    assert(((static_cast<uint64_t>(prime - 1) * inv_max) % prime) == 1);
    std::cout << "  ✓ inv(" << (prime - 1) << ") mod " << prime << " = " << inv_max << std::endl;
    
    std::cout << "✓ modular_inverse function tests passed" << std::endl;
}

int main() {
    try {
        std::cout << "Note: Most bivariate algorithm tests are disabled due to refactoring." << std::endl;
        std::cout << "The new biv_lex implementation is tested via integration tests:" << std::endl;
        std::cout << "  - test_parabola_line_debug" << std::endl;
        std::cout << "  - test_circle_line_fix" << std::endl;
        std::cout << "  - test_biv_lex_bug" << std::endl;
        std::cout << std::endl;
        
        // Only test the functions that still exist
        test_modular_inverse();
        
        std::cout << "\n✓ Limited bivariate algorithm tests passed!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with error: " << e.what() << std::endl;
        return 1;
    }
}