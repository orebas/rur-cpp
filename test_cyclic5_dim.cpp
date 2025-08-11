#include <iostream>
#include <vector>
#include <string>
#include "julia_rur/hyperplane_sections.hpp"

int main() {
    std::vector<std::string> polynomials = {
        "x0+x1+x2+x3+x4",
        "x0*x1+x1*x2+x2*x3+x3*x4+x4*x0",
        "x0*x1*x2+x1*x2*x3+x2*x3*x4+x3*x4*x0+x4*x0*x1",
        "x0*x1*x2*x3+x1*x2*x3*x4+x2*x3*x4*x0+x3*x4*x0*x1+x4*x0*x1*x2",
        "x0*x1*x2*x3*x4-1"
    };
    std::vector<std::string> variables = {"x0", "x1", "x2", "x3", "x4"};
    
    std::cout << "Testing Cyclic5 dimension detection:\n";
    
    julia_rur::DimensionAnalysis dim = julia_rur::analyze_system_dimension(polynomials, variables);
    
    std::cout << "Dimension: " << dim.dimension << std::endl;
    std::cout << "Is zero-dimensional: " << (dim.is_zero_dimensional ? "yes" : "no") << std::endl;
    std::cout << "Free variables: ";
    for (int v : dim.free_variables) {
        std::cout << "x" << v << " ";
    }
    std::cout << std::endl;
    std::cout << "Method used: " << dim.method_used << std::endl;
    
    return 0;
}