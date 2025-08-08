#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <cmath>
#include "../src/julia_rur/numerical_roots_flint.hpp"

using namespace julia_rur;

void test_basic_polynomial_root_finding() {
    std::cout << "Testing FLINT/ARB polynomial root finding..." << std::endl;
    
    // Test simple quadratic: x^2 - 1 = 0 (roots: 1, -1)
    std::vector<mpq_class> quadratic = {-1, 0, 1}; // -1 + 0*x + 1*x^2
    
    FlintRootConfig config;
    auto roots = find_polynomial_roots_flint(quadratic, config);
    
    std::cout << "Found " << roots.size() << " roots for x^2 - 1 = 0:" << std::endl;
    for (size_t i = 0; i < roots.size(); ++i) {
        std::cout << "  Root " << i << ": " << roots[i].real() << " + " << roots[i].imag() << "i" << std::endl;
        
        // Verify this is actually a root
        std::complex<double> z = roots[i];
        std::complex<double> value = z * z - 1.0;
        double error = std::abs(value);
        std::cout << "    Verification: |f(" << z << ")| = " << error << std::endl;
        assert(error < 1e-10);
    }
    
    std::cout << "✓ Basic root finding test passed" << std::endl;
}

void test_certified_root_finding() {
    std::cout << "\nTesting certified root finding with error bounds..." << std::endl;
    
    // Test cubic: x^3 - 2*x + 1 = 0
    std::vector<mpq_class> cubic = {1, -2, 0, 1}; // 1 - 2*x + 0*x^2 + 1*x^3
    
    FlintRootConfig config;
    std::vector<double> radii;
    auto roots = find_polynomial_roots_certified(cubic, radii, config);
    
    std::cout << "Found " << roots.size() << " roots for x^3 - 2x + 1 = 0:" << std::endl;
    for (size_t i = 0; i < roots.size(); ++i) {
        std::cout << "  Root " << i << ": " << roots[i].real() << " + " << roots[i].imag() 
                  << "i (radius: " << radii[i] << ")" << std::endl;
    }
    
    assert(roots.size() > 0);
    assert(radii.size() == roots.size());
    
    std::cout << "✓ Certified root finding test passed" << std::endl;
}

void test_configuration_options() {
    std::cout << "\nTesting configuration options..." << std::endl;
    
    FlintRootConfig config;
    config.initial_prec = 256;
    config.max_prec = 1024;
    config.epsilon = 1e-15;
    config.isolate_real_roots = true;
    
    // Simple polynomial with known real and complex roots
    std::vector<mpq_class> poly = {-1, 0, 1}; // x^2 - 1 = 0
    auto roots = find_polynomial_roots_flint(poly, config);
    
    std::cout << "Configuration test with high precision settings completed" << std::endl;
    std::cout << "✓ Configuration options test passed" << std::endl;
}

int main() {
    std::cout << "=== FLINT/ARB Support Test Suite ===" << std::endl;
    
#ifdef HAVE_FLINT_ARB
    std::cout << "ARB support: ENABLED (using libarb)" << std::endl;
#else
    std::cout << "ARB support: DISABLED (using Eigen fallback)" << std::endl;
#endif
    
    try {
        test_basic_polynomial_root_finding();
        test_certified_root_finding();
        test_configuration_options();
        
        std::cout << "\n=== All tests passed! ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
}