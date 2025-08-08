#include <unsupported/Eigen/Polynomials>
#include <iostream>

int main() {
    // Test x^2 - 2
    // According to Eigen docs, coefficients should be high to low degree
    
    // Try different orders
    std::cout << "Testing x^2 - 2 with Eigen PolynomialSolver" << std::endl;
    
    // Method 1: High to low [1, 0, -2]
    {
        Eigen::Vector3d coeffs;
        coeffs << 1, 0, -2;  // x^2 + 0*x - 2
        
        Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
        solver.compute(coeffs);
        
        std::cout << "\nMethod 1 (high to low): coeffs = [1, 0, -2]" << std::endl;
        auto roots = solver.roots();
        for (int i = 0; i < roots.size(); ++i) {
            std::cout << "  Root " << i << ": " << roots[i] << std::endl;
        }
    }
    
    // Method 2: Low to high [-2, 0, 1]
    {
        Eigen::Vector3d coeffs;
        coeffs << -2, 0, 1;  // -2 + 0*x + x^2
        
        Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
        solver.compute(coeffs);
        
        std::cout << "\nMethod 2 (low to high): coeffs = [-2, 0, 1]" << std::endl;
        auto roots = solver.roots();
        for (int i = 0; i < roots.size(); ++i) {
            std::cout << "  Root " << i << ": " << roots[i] << std::endl;
        }
    }
    
    // Test simple polynomial x - 2
    {
        Eigen::Vector2d coeffs;
        coeffs << 1, -2;  // x - 2
        
        Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
        solver.compute(coeffs);
        
        std::cout << "\nSimple test x - 2: coeffs = [1, -2]" << std::endl;
        auto roots = solver.roots();
        for (int i = 0; i < roots.size(); ++i) {
            std::cout << "  Root " << i << ": " << roots[i] << std::endl;
        }
    }
    
    return 0;
}