#include <iostream>
#include <vector>
#include <string>

int main() {
    std::cout << "==============================================\n";
    std::cout << "Katsura-4 System as used in C++ test\n";
    std::cout << "==============================================\n\n";
    
    std::vector<std::string> polynomials = {
        "x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1",
        "x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1",
        "x0*x2 + 2*x1*x3 + 2*x2*x4 - x2", 
        "x0*x3 + 2*x1*x4 - x3",
        "x0*x4 - x4"
    };
    
    std::cout << "C++ code uses these polynomials:\n";
    std::cout << "---------------------------------\n";
    for (size_t i = 0; i < polynomials.size(); ++i) {
        std::cout << "f" << (i+1) << " = " << polynomials[i] << std::endl;
    }
    
    std::cout << "\nVariables: x0, x1, x2, x3, x4\n";
    std::cout << "\nThis is a system of " << polynomials.size() 
              << " polynomial equations in 5 variables.\n";
    
    return 0;
}