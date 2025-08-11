#include <iostream>
#include <vector>
#include "src/julia_rur/rur_main_algorithm.hpp"

using namespace julia_rur;

int main() {
    // Katsura-4 polynomials
    std::vector<std::string> polynomials = {
        "x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1",
        "x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1",
        "x0*x2 + 2*x1*x3 + 2*x2*x4 - x2", 
        "x0*x3 + 2*x1*x4 - x3",
        "x0*x4 - x4"
    };
    
    std::vector<std::string> variables = {"x0", "x1", "x2", "x3", "x4"};
    
    // Try multiple primes
    std::vector<ModularCoeff> primes = {
        1073741827,  // 30-bit prime
        1073741783,  // Another 30-bit
        268435459,   // 28-bit prime
        536870923,   // 29-bit prime
        2,           // The "bad" prime
        3, 5, 7, 11, 13  // Small primes
    };
    
    std::cout << "Testing Katsura-4 quotient dimension with different primes:\n";
    std::cout << "Expected: 16 (2^4 for Katsura-4)\n\n";
    
    for (ModularCoeff p : primes) {
        int dim = compute_quotient_dimension(polynomials, variables, p);
        std::cout << "Prime " << p << ": dimension = " << dim;
        if (dim != 16) {
            std::cout << " <- BAD PRIME (expected 16)";
        }
        std::cout << "\n";
    }
    
    return 0;
}