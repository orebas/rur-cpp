#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>
#include <gmpxx.h>
#include <vector>
#include <cmath>

int main() {
    std::cout << "=== Testing Eigen roots for T^3 - 2T - 5 ===\n" << std::endl;
    
    // Coefficients in low-to-high order: -5, -2, 0, 1
    Eigen::VectorXd coeffs(4);
    coeffs[0] = -5.0;  // constant
    coeffs[1] = -2.0;  // T^1
    coeffs[2] = 0.0;   // T^2
    coeffs[3] = 1.0;   // T^3
    
    std::cout << "Polynomial: T^3 - 2T - 5 = 0" << std::endl;
    std::cout << "Coefficients (low to high): " << coeffs.transpose() << std::endl;
    
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    solver.compute(coeffs);
    
    const auto& roots = solver.roots();
    std::cout << "\nRoots:" << std::endl;
    for (int i = 0; i < roots.size(); ++i) {
        std::cout << "  T = " << roots[i].real();
        if (std::abs(roots[i].imag()) > 1e-10) {
            std::cout << " + " << roots[i].imag() << "i";
        }
        
        // Verify
        auto t = roots[i];
        auto value = t*t*t - 2.0*t - 5.0;
        std::cout << "  (check: " << std::abs(value) << ")";
        std::cout << std::endl;
    }
    
    std::cout << "\nExpected real root: 2.094551482" << std::endl;
    
    return 0;
}