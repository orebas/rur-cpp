#include <iostream>
#include "julia_rur/rur_main_algorithm.hpp"

int main() {
    std::vector<std::string> polynomials = {"x - 4"};
    std::vector<std::string> variables = {"x"};
    
    julia_rur::RURConfig config;
    config.verbose = true;
    
    // Compute rational RUR
    auto rational_result = julia_rur::compute_rational_rur(
        polynomials, variables, config
    );
    
    std::cout << "\n=== Rational RUR Result ===" << std::endl;
    std::cout << "Success: " << rational_result.success << std::endl;
    if (!rational_result.success) {
        std::cout << "Error: " << rational_result.error_message << std::endl;
    } else {
        std::cout << "Quotient basis size: " << rational_result.quotient_basis.size() << std::endl;
        std::cout << "Minimal polynomial size: " << rational_result.minimal_polynomial.size() << std::endl;
        std::cout << "Minimal polynomial: ";
        for (const auto& c : rational_result.minimal_polynomial) {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        
        std::cout << "Number of numerators: " << rational_result.numerators.size() << std::endl;
        for (size_t i = 0; i < rational_result.numerators.size(); ++i) {
            std::cout << "Numerator " << i << ": ";
            for (const auto& c : rational_result.numerators[i]) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }
    }
    
    return 0;
}