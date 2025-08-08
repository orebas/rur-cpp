#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>
#include <iomanip>

using namespace julia_rur;

int main() {
    // Create a generic dense system with 3 variables, total degree up to 4
    // This should have a finite number of isolated solutions
    // Using dense polynomials with "random" coefficients (actually chosen to be generic)
    
    std::vector<std::string> test_systems[][2] = {
        // Simple test first - should have 2 solutions
        {
            {"x^2 - 1", "y - x"},
            {"x", "y"}
        },
        
        // Slightly harder - should have 4 solutions  
        {
            {"x^2 + y^2 - 5", "x*y - 2"},
            {"x", "y"}
        },
        
        // Generic dense system degree 2 - should have 8 solutions typically
        {
            {"x^2 + 2*x*y + 3*y^2 + 4*x + 5*y - 6",
             "2*x^2 - x*y + y^2 + 3*x - 2*y + 1"},
            {"x", "y"}
        },
        
        // Generic dense system degree 3 with 3 variables - this is what you asked for
        // Using integer coefficients that should give a generic system
        {
            {"2*x^3 + 3*x^2*y - x^2*z + 4*x*y^2 + 2*x*y*z - x*z^2 + 5*y^3 - y^2*z + 3*y*z^2 + z^3 + x^2 - 2*x*y + 3*x*z + y^2 - y*z + 2*z^2 + x - y + z - 1",
             "x^3 - 2*x^2*y + x^2*z + 3*x*y^2 - x*y*z + 2*x*z^2 - y^3 + 4*y^2*z - y*z^2 + 2*z^3 - x^2 + x*y - x*z + 2*y^2 + y*z - z^2 + 2*x + y - z + 3",
             "3*x^3 + x^2*y + 2*x^2*z - x*y^2 + 3*x*y*z - 2*x*z^2 + 2*y^3 + y^2*z - 2*y*z^2 - z^3 + 2*x^2 + x*y + x*z - y^2 + 2*y*z + z^2 - x + 2*y - 3*z - 2"},
            {"x", "y", "z"}
        },
        
        // Even more generic - degree 4 terms
        {
            {"x^4 + 2*x^3*y - x^3*z + x^2*y^2 + 3*x^2*y*z - 2*x^2*z^2 + x*y^3 - 2*x*y^2*z + x*y*z^2 + x*z^3 + y^4 - y^3*z + 2*y^2*z^2 - y*z^3 + z^4 + x^2 + y^2 - z^2 + x - y + 1",
             "2*x^4 - x^3*y + x^3*z - 2*x^2*y^2 + x^2*y*z + x^2*z^2 - x*y^3 + 3*x*y^2*z - x*y*z^2 - 2*x*z^3 - y^4 + 2*y^3*z - y^2*z^2 + y*z^3 - 2*z^4 - x^2 + 2*y^2 + z^2 - x + y - z - 2",  
             "x^4 + x^3*y + x^3*z + x^2*y^2 - x^2*y*z - x^2*z^2 + 2*x*y^3 - x*y^2*z - 2*x*y*z^2 + x*z^3 + 2*y^4 + y^3*z - y^2*z^2 - 2*y*z^3 + z^4 + 2*x^2 - y^2 - 2*z^2 + 2*x + y + z + 3"},
            {"x", "y", "z"}
        }
    };
    
    EnhancedSolverConfig config;
    config.verbose = false;
    
    for (size_t i = 0; i < sizeof(test_systems)/sizeof(test_systems[0]); ++i) {
        auto& polys = test_systems[i][0];
        auto& vars = test_systems[i][1];
        
        std::cout << "\n============================================\n";
        std::cout << "Test " << (i+1) << ": ";
        if (vars.size() == 2) {
            std::cout << "2 variables, ";
        } else {
            std::cout << "3 variables, ";
        }
        
        // Determine max degree
        int max_deg = 2;
        for (const auto& p : polys) {
            if (p.find("^4") != std::string::npos) max_deg = 4;
            else if (p.find("^3") != std::string::npos && max_deg < 3) max_deg = 3;
        }
        std::cout << "degree " << max_deg << " system\n";
        std::cout << "--------------------------------------------\n";
        
        // Show the system
        std::cout << "System:\n";
        for (size_t j = 0; j < polys.size(); ++j) {
            std::cout << "  f" << (j+1) << " = " << polys[j] << "\n";
        }
        std::cout << "Variables: ";
        for (const auto& v : vars) std::cout << v << " ";
        std::cout << "\n\n";
        
        // Solve it
        auto result = solve_polynomial_system_enhanced(polys, vars, config);
        
        if (result.success) {
            std::cout << "✓ SUCCESS: Found " << result.solutions.size() << " solutions\n";
            
            // For small systems, show the solutions
            if (result.solutions.size() <= 8) {
                for (size_t j = 0; j < result.solutions.size(); ++j) {
                    std::cout << "  Solution " << (j+1) << ": (";
                    for (size_t k = 0; k < vars.size(); ++k) {
                        if (k > 0) std::cout << ", ";
                        auto val = result.solutions[j][k];
                        std::cout << std::fixed << std::setprecision(4);
                        if (std::abs(val.imag()) < 1e-10) {
                            std::cout << val.real();
                        } else {
                            std::cout << val.real() << "+" << val.imag() << "i";
                        }
                    }
                    std::cout << ")\n";
                }
            }
        } else {
            std::cout << "✗ FAILED: " << result.error_message << "\n";
            
            // Try with verbose to see what's happening
            if (i == 3) {  // The 3-variable degree-3 system you asked for
                std::cout << "\nRetrying with verbose output for debugging...\n";
                config.verbose = true;
                auto debug_result = solve_polynomial_system_enhanced(polys, vars, config);
                config.verbose = false;
            }
        }
    }
    
    std::cout << "\n============================================\n";
    std::cout << "Testing complete.\n";
    
    return 0;
}