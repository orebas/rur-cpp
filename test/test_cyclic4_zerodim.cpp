#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <cmath>

// Test the overdetermined cyclic-4 system
// This creates a file that axf4 can process

int main() {
    std::cout << "Creating overdetermined cyclic-4 system file for axf4\n";
    std::cout << "=====================================================\n\n";
    
    // Write the polynomials to a file
    std::ofstream outfile("cyclic4_zerodim.txt");
    
    // Original cyclic-4 system
    outfile << "x0+x1+x2+x3\n";
    outfile << "x0*x1+x1*x2+x2*x3+x3*x0\n";
    outfile << "x0*x1*x2+x1*x2*x3+x2*x3*x0+x3*x0*x1\n";
    outfile << "x0*x1*x2*x3-1\n";
    // Additional linear constraint to make it zero-dimensional
    outfile << "x0+2*x1+3*x2+5*x3-1\n";
    
    outfile.close();
    
    std::cout << "System of equations (5 equations, 4 variables):\n";
    std::cout << "  x0 + x1 + x2 + x3 = 0\n";
    std::cout << "  x0*x1 + x1*x2 + x2*x3 + x3*x0 = 0\n";
    std::cout << "  x0*x1*x2 + x1*x2*x3 + x2*x3*x0 + x3*x0*x1 = 0\n";
    std::cout << "  x0*x1*x2*x3 - 1 = 0\n";
    std::cout << "  x0 + 2*x1 + 3*x2 + 5*x3 - 1 = 0  (additional constraint)\n";
    std::cout << "\n";
    
    std::cout << "File 'cyclic4_zerodim.txt' created for axf4 processing.\n";
    std::cout << "\n";
    std::cout << "This overdetermined system (5 equations in 4 variables) should have\n";
    std::cout << "fewer solutions than the original cyclic-4 system (which has 24).\n";
    std::cout << "\n";
    std::cout << "The additional linear constraint creates a hyperplane that intersects\n";
    std::cout << "the original solution variety, reducing the number of solutions.\n";
    
    return 0;
}