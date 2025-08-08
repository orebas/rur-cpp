#include <iostream>
#include <vector>
#include <string>
#include "../src/julia_rur/numerical_roots_eigen.hpp"

using namespace julia_rur;

int main() {
    std::cout << "Testing circle-line intersection with fractional coefficients\n" << std::endl;
    
    // Test case 1: Circle and horizontal line with decimal
    {
        std::cout << "Test 1: Circle x^2 + y^2 - 1 = 0 and line y - 0.5 = 0" << std::endl;
        std::vector<std::string> polynomials = {"x^2 + y^2 - 1", "y - 0.5"};
        std::vector<std::string> variables = {"x", "y"};
        
        try {
            auto result = solve_polynomial_system(polynomials, variables);
            
            if (result.success) {
                std::cout << "  SUCCESS: System solved!" << std::endl;
                std::cout << "  Number of solutions: " << result.solutions.size() << std::endl;
                std::cout << "  Expected: 2 solutions at (±√3/2, 1/2)" << std::endl;
                
                for (size_t i = 0; i < result.solutions.size(); i++) {
                    std::cout << "  Solution " << (i+1) << ": (";
                    for (size_t j = 0; j < result.solutions[i].size(); j++) {
                        if (j > 0) std::cout << ", ";
                        auto& val = result.solutions[i][j];
                        std::cout << val.real();
                        if (std::abs(val.imag()) > 1e-10) {
                            std::cout << " + " << val.imag() << "i";
                        }
                    }
                    std::cout << ")" << std::endl;
                }
            } else {
                std::cout << "  FAILED: " << result.error_message << std::endl;
            }
        } catch (const std::exception& e) {
            std::cout << "  EXCEPTION: " << e.what() << std::endl;
        }
    }
    
    std::cout << std::endl;
    
    // Test case 2: More complex fractional coefficients
    {
        std::cout << "Test 2: System with multiple fractional coefficients" << std::endl;
        std::vector<std::string> polynomials = {"0.5*x^2 + 0.25*y^2 - 1", "x + 0.5*y - 1"};
        std::vector<std::string> variables = {"x", "y"};
        
        try {
            auto result = solve_polynomial_system(polynomials, variables);
            
            if (result.success) {
                std::cout << "  SUCCESS: System solved!" << std::endl;
                std::cout << "  Number of solutions: " << result.solutions.size() << std::endl;
                
                for (size_t i = 0; i < result.solutions.size(); i++) {
                    std::cout << "  Solution " << (i+1) << ": (";
                    for (size_t j = 0; j < result.solutions[i].size(); j++) {
                        if (j > 0) std::cout << ", ";
                        auto& val = result.solutions[i][j];
                        std::cout << val.real();
                        if (std::abs(val.imag()) > 1e-10) {
                            std::cout << " + " << val.imag() << "i";
                        }
                    }
                    std::cout << ")" << std::endl;
                }
            } else {
                std::cout << "  FAILED: " << result.error_message << std::endl;
            }
        } catch (const std::exception& e) {
            std::cout << "  EXCEPTION: " << e.what() << std::endl;
        }
    }
    
    std::cout << std::endl;
    
    // Test case 3: The original problematic case
    {
        std::cout << "Test 3: The y - 0.5 case (without spaces)" << std::endl;
        std::vector<std::string> polynomials = {"x^2 + y^2 - 1", "y- 0.5"};
        std::vector<std::string> variables = {"x", "y"};
        
        try {
            auto result = solve_polynomial_system(polynomials, variables);
            
            if (result.success) {
                std::cout << "  SUCCESS: System solved!" << std::endl;
                std::cout << "  Number of solutions: " << result.solutions.size() << std::endl;
            } else {
                std::cout << "  FAILED: " << result.error_message << std::endl;
            }
        } catch (const std::exception& e) {
            std::cout << "  EXCEPTION: " << e.what() << std::endl;
        }
    }
    
    return 0;
}