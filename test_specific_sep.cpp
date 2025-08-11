#include <iostream>
#include <vector>
#include "../src/julia_rur/rur_main_algorithm.hpp"
#include "../src/polynomial_parser.hpp"

int main() {
    // Test the 3-variable symmetric system with Julia's separating element
    std::vector<std::string> poly_strs = {
        "x + y + z - 6",
        "x*y + y*z + z*x - 11",
        "x*y*z - 6"
    };
    std::vector<std::string> var_names = {"x", "y", "z"};
    
    std::cout << "System: elementary symmetric polynomials" << std::endl;
    std::cout << "Solutions: permutations of (1,2,3)" << std::endl;
    std::cout << "\nJulia found separating element: -2x - y + z" << std::endl;
    std::cout << "This gives a degree-6 minimal polynomial" << std::endl;
    
    // Test with specific prime
    ModularCoeff prime = 1073741827;  // 30-bit prime
    std::cout << "\nUsing prime: " << prime << std::endl;
    
    // The separating element coefficients from Julia: [-2, -1, 1]
    // In modular arithmetic: 
    // -2 mod p = p - 2
    // -1 mod p = p - 1
    //  1 mod p = 1
    ModularCoeff c1 = prime - 2;
    ModularCoeff c2 = prime - 1;
    ModularCoeff c3 = 1;
    
    std::cout << "Separating element coefficients (mod p):" << std::endl;
    std::cout << "  x: " << c1 << " (= -2 mod p)" << std::endl;
    std::cout << "  y: " << c2 << " (= -1 mod p)" << std::endl;
    std::cout << "  z: " << c3 << " (= 1 mod p)" << std::endl;
    
    return 0;
}
