#include <iostream>
#include <vector>
#include <string>

int main() {
    std::cout << "===========================================\n";
    std::cout << "The Katsura-4 System (5 equations, 5 variables)\n";
    std::cout << "===========================================\n\n";
    
    std::cout << "Variables: x0, x1, x2, x3, x4\n\n";
    
    std::cout << "Equations:\n";
    std::cout << "----------\n";
    std::cout << "f1: x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1 = 0\n";
    std::cout << "f2: x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1 = 0\n";
    std::cout << "f3: x0*x2 + 2*x1*x3 + 2*x2*x4 - x2 = 0\n";
    std::cout << "f4: x0*x3 + 2*x1*x4 - x3 = 0\n";
    std::cout << "f5: x0*x4 - x4 = 0\n\n";
    
    std::cout << "Note about f5:\n";
    std::cout << "  x0*x4 - x4 = 0\n";
    std::cout << "  x4*(x0 - 1) = 0\n";
    std::cout << "  So either x4 = 0 OR x0 = 1\n\n";
    
    std::cout << "If x4 = 0, then from f4:\n";
    std::cout << "  x0*x3 - x3 = 0 → x3*(x0 - 1) = 0\n";
    std::cout << "  So either x3 = 0 OR x0 = 1\n\n";
    
    std::cout << "This cascade structure is what gives Katsura-4 its special properties.\n\n";
    
    std::cout << "Known facts:\n";
    std::cout << "  - This system has exactly 16 complex solutions\n";
    std::cout << "  - The quotient ring Q[x]/I has dimension 16\n";
    std::cout << "  - Some solutions have x4 = 0, others have x0 = 1\n";
    std::cout << "  - The system exhibits symmetries that make single variables poor separating elements\n";
    
    return 0;
}