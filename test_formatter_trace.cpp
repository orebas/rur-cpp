#include <iostream>
#include <string>
#include <vector>
#include "src/julia_rur/f4_polynomial_formatter.hpp"

int main() {
    std::vector<std::string> test_inputs = {
        "x^2 + y^2 - 1",
        "x^2+y^2-1", 
        "x + y - 1",
        "2*x + 2*y - 2"
    };
    
    for (const auto& input : test_inputs) {
        std::string formatted = julia_rur::format_polynomial_for_f4(input);
        std::cout << "Input:  \"" << input << "\"" << std::endl;
        std::cout << "Output: \"" << formatted << "\"" << std::endl << std::endl;
    }
    
    return 0;
}
