#include <iostream>
#include <vector>
#include <string>
#include "julia_rur/hyperplane_sections.hpp"
#include "julia_rur/polynomial_solver_enhanced.hpp"

int main() {
    // Cyclic-5 system
    std::vector<std::string> polynomials = {
        "x0+x1+x2+x3+x4",
        "x0*x1+x1*x2+x2*x3+x3*x4+x4*x0",
        "x0*x1*x2+x1*x2*x3+x2*x3*x4+x3*x4*x0+x4*x0*x1",
        "x0*x1*x2*x3+x1*x2*x3*x4+x2*x3*x4*x0+x3*x4*x0*x1+x4*x0*x1*x2",
        "x0*x1*x2*x3*x4-1"
    };
    std::vector<std::string> variables = {"x0", "x1", "x2", "x3", "x4"};
    
    std::cout << "Testing dimension detection for Cyclic-5 system...\n";
    
    // Test with a few different primes
    for (int i = 0; i < 5; i++) {
        auto dim_analysis = julia_rur::analyze_system_dimension(polynomials, variables);
        std::cout << "Attempt " << i+1 << ":\n";
        std::cout << "  Dimension: " << dim_analysis.dimension << "\n";
        std::cout << "  Is zero-dimensional: " << (dim_analysis.is_zero_dimensional ? "yes" : "no") << "\n";
        std::cout << "  Method used: " << dim_analysis.method_used << "\n";
        std::cout << "  Free variables: ";
        for (int v : dim_analysis.free_variables) {
            std::cout << v << " ";
        }
        std::cout << "\n";
    }
    
    // Test with robust analysis
    julia_rur::HyperplaneSectionConfig config;
    auto robust_dim = julia_rur::analyze_system_dimension_robust(polynomials, variables, config);
    std::cout << "\nRobust analysis:\n";
    std::cout << "  Dimension: " << robust_dim.dimension << "\n";
    std::cout << "  Is zero-dimensional: " << (robust_dim.is_zero_dimensional ? "yes" : "no") << "\n";
    std::cout << "  Method used: " << robust_dim.method_used << "\n";
    
    return 0;
}