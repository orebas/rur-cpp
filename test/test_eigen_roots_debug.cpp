#include "../src/julia_rur/numerical_roots_eigen.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace julia_rur;

int main() {
    std::cout << "=== Testing Eigen polynomial root finding ===" << std::endl;
    
    // Test polynomial T^2 - 2 = 0
    std::vector<mpq_class> coeffs = {-2, 0, 1};
    
    std::cout << "\nPolynomial coefficients (low to high): ";
    for (const auto& c : coeffs) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    
    // Find roots using our function
    auto roots = find_polynomial_roots(coeffs);
    
    std::cout << "\nRoots found by find_polynomial_roots:" << std::endl;
    for (const auto& root : roots) {
        std::cout << "  " << root.real();
        if (std::abs(root.imag()) > 1e-10) {
            std::cout << " + " << root.imag() << "i";
        }
        std::cout << std::endl;
    }
    
    // Test directly with Eigen
    std::cout << "\nDirect Eigen test:" << std::endl;
    Eigen::VectorXd eigen_coeffs(3);
    eigen_coeffs[0] = -2.0;  // constant term
    eigen_coeffs[1] = 0.0;   // linear term
    eigen_coeffs[2] = 1.0;   // quadratic term
    
    std::cout << "Eigen coefficient vector: " << eigen_coeffs.transpose() << std::endl;
    
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    solver.compute(eigen_coeffs);
    
    const auto& eigen_roots = solver.roots();
    std::cout << "Eigen roots:" << std::endl;
    for (int i = 0; i < eigen_roots.size(); ++i) {
        std::cout << "  " << eigen_roots[i].real();
        if (std::abs(eigen_roots[i].imag()) > 1e-10) {
            std::cout << " + " << eigen_roots[i].imag() << "i";
        }
        std::cout << std::endl;
    }
    
    // Verify roots
    std::cout << "\nVerification (should be ±√2 ≈ ±1.414):" << std::endl;
    std::cout << "Expected: " << sqrt(2.0) << " and " << -sqrt(2.0) << std::endl;
    
    return 0;
}