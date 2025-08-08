#include "../src/julia_rur/rur_main_algorithm.hpp"
#include <iostream>
#include <vector>

using namespace julia_rur;

void test_simple_system() {
    std::cout << "\n=== Testing Simple System: Circle and Line ===" << std::endl;
    std::cout << "System: {x^2 + y^2 - 1, x - y}" << std::endl;
    
    std::vector<std::string> polynomials = {
        "1*x^2+1*y^2-1",  // x^2 + y^2 - 1
        "1*x-1*y"         // x - y
    };
    
    std::vector<std::string> variables = {"x", "y"};
    
    // Check if zero-dimensional
    if (is_zero_dimensional_system(polynomials, variables)) {
        std::cout << "System is zero-dimensional ✓" << std::endl;
        
        // Compute dimension
        int dim = compute_quotient_dimension(polynomials, variables);
        std::cout << "Quotient ring dimension: " << dim << std::endl;
        
        // Compute rational RUR
        RURConfig config;
        config.verbose = false;
        
        std::cout << "\nComputing rational RUR..." << std::endl;
        RationalRURResult result = compute_rational_rur(
            polynomials, variables, config
        );
        
        if (result.success) {
            std::cout << "\n" << format_rur_result(result, variables) << std::endl;
        } else {
            std::cout << "RUR computation failed: " << result.error_message << std::endl;
        }
        
    } else {
        std::cout << "System is not zero-dimensional" << std::endl;
    }
}

void test_cyclic3_system() {
    std::cout << "\n=== Testing Cyclic-3 System ===" << std::endl;
    std::cout << "System: {x+y+z, xy+yz+zx, xyz-1}" << std::endl;
    
    std::vector<std::string> polynomials = {
        "1*x+1*y+1*z",           // x + y + z
        "1*x*y+1*y*z+1*z*x",     // xy + yz + zx
        "1*x*y*z-1"              // xyz - 1
    };
    
    std::vector<std::string> variables = {"x", "y", "z"};
    
    // Check dimension
    int dim = compute_quotient_dimension(polynomials, variables);
    std::cout << "Quotient ring dimension: " << dim << std::endl;
    
    if (dim > 0) {
        // Compute rational RUR
        RURConfig config;
        config.verbose = false;
        
        std::cout << "\nComputing rational RUR..." << std::endl;
        RationalRURResult result = compute_rational_rur(
            polynomials, variables, config
        );
        
        if (result.success) {
            std::cout << "\n" << format_rur_result(result, variables) << std::endl;
        } else {
            std::cout << "RUR computation failed: " << result.error_message << std::endl;
        }
    }
}

void demo_user_input() {
    std::cout << "\n=== Interactive RUR Demo ===" << std::endl;
    std::cout << "Enter polynomial system (press Ctrl+D when done):" << std::endl;
    std::cout << "Example format: 1*x^2+1*y^2+100002" << std::endl;
    
    // This is a placeholder for interactive input
    // In a full implementation, would parse user input
}

int main(int argc, char* argv[]) {
    std::cout << "Julia-style RUR Implementation Demo" << std::endl;
    std::cout << "===================================" << std::endl;
    
    try {
        test_simple_system();
        test_cyclic3_system();
        
        if (argc > 1 && std::string(argv[1]) == "--interactive") {
            demo_user_input();
        }
        
        std::cout << "\nDemo completed successfully!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}