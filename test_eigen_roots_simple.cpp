#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Dense>

std::vector<std::complex<double>> find_roots_simple(const std::vector<double>& coeffs) {
    int n = coeffs.size() - 1;  // degree
    if (n == 0) return {};
    
    // For linear case: ax + b = 0 => x = -b/a
    if (n == 1) {
        if (coeffs[1] != 0) {
            return {std::complex<double>(-coeffs[0] / coeffs[1], 0)};
        }
        return {};
    }
    
    // Build companion matrix
    Eigen::MatrixXd companion(n, n);
    companion.setZero();
    
    // Last column contains -coeffs[0..n-1] / coeffs[n]
    for (int i = 0; i < n; ++i) {
        companion(i, n-1) = -coeffs[i] / coeffs[n];
    }
    
    // Sub-diagonal contains 1s
    for (int i = 1; i < n; ++i) {
        companion(i, i-1) = 1.0;
    }
    
    // Find eigenvalues
    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    std::vector<std::complex<double>> roots;
    for (int i = 0; i < n; ++i) {
        roots.push_back(solver.eigenvalues()[i]);
    }
    
    return roots;
}

int main() {
    // Test 1: T - 3 = 0
    std::vector<double> coeffs1 = {-3, 1};
    std::cout << "Polynomial: T - 3 = 0" << std::endl;
    auto roots1 = find_roots_simple(coeffs1);
    for (const auto& root : roots1) {
        std::cout << "  Root: " << root << std::endl;
    }
    
    // Test 2: T + 3 = 0
    std::vector<double> coeffs2 = {3, 1};
    std::cout << "\nPolynomial: T + 3 = 0" << std::endl;
    auto roots2 = find_roots_simple(coeffs2);
    for (const auto& root : roots2) {
        std::cout << "  Root: " << root << std::endl;
    }
    
    // Test 3: T^2 - 2 = 0
    std::vector<double> coeffs3 = {-2, 0, 1};
    std::cout << "\nPolynomial: T^2 - 2 = 0" << std::endl;
    auto roots3 = find_roots_simple(coeffs3);
    for (const auto& root : roots3) {
        std::cout << "  Root: " << root << std::endl;
    }
    
    return 0;
}