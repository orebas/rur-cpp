#include <iostream>
#include <vector>
#include <string>
#include "../src/julia_rur/f4_polynomial_formatter.hpp"

int main() {
    std::cout << "Testing new polynomial parser:" << std::endl;
    
    // Test cases
    std::vector<std::string> test_cases = {
        "(x-1)^2",
        "(x-1)^2 + y^2 - 1",
        "(x+2)^3",
        "x^2 - 2*x + 1",  // Expanded form of (x-1)^2
        "(x-1)*(x+1)",
        "x*(x-1)^2"
    };
    
    for (const auto& test : test_cases) {
        std::cout << "\nInput: " << test << std::endl;
        std::string formatted = julia_rur::format_polynomial_for_f4(test);
        std::cout << "Output: " << formatted << std::endl;
    }
    
    return 0;
}