#include <iostream>
#include "julia_rur/rur_main_algorithm.hpp"
#include "julia_rur/polynomial_solver_enhanced.hpp"

int main() {
    std::vector<std::string> polynomials = {"x - 4"};
    std::vector<std::string> variables = {"x"};
    
    julia_rur::RURConfig config;
    config.verbose = true;
    
    // Compute modular RUR
    julia_rur::ModularCoeff prime = 1073741827;
    auto modular_result = julia_rur::compute_modular_rur(
        polynomials, variables, prime, config, {}
    );
    
    std::cout << "\n=== Modular RUR Result ===" << std::endl;
    std::cout << "Success: " << modular_result.success << std::endl;
    std::cout << "Quotient basis size: " << modular_result.quotient_basis.size() << std::endl;
    std::cout << "Minimal polynomial coefficients: ";
    for (auto c : modular_result.minimal_polynomial.coefficients) {
        std::cout << c << " ";
    }
    std::cout << std::endl;
    std::cout << "Minimal polynomial degree: " << modular_result.minimal_polynomial.degree << std::endl;
    
    std::cout << "Number of parameterizations: " << modular_result.parameterizations.size() << std::endl;
    for (size_t i = 0; i < modular_result.parameterizations.size(); ++i) {
        const auto& param = modular_result.parameterizations[i];
        std::cout << "Param " << i << " success: " << param.success << std::endl;
        std::cout << "Param " << i << " generators size: " << param.generators.size() << std::endl;
        if (!param.generators.empty() && !param.generators[0].empty()) {
            std::cout << "Param " << i << " generator[0]: ";
            for (auto c : param.generators[0]) {
                std::cout << c << " ";
            }
            std::cout << std::endl;
        }
    }
    
    return 0;
}