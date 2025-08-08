#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>

struct QuickTest {
    std::string name;
    std::vector<std::string> polynomials;
    std::vector<std::string> variables;
    int expected_solutions;
};

bool run_quick_test(const QuickTest& test) {
    julia_rur::EnhancedSolverConfig config;
    config.verbose = false;
    
    // Redirect stdout and stderr to suppress debug output
    std::streambuf* orig_cout = std::cout.rdbuf();
    std::streambuf* orig_cerr = std::cerr.rdbuf();
    std::ostringstream null_stream;
    std::cout.rdbuf(null_stream.rdbuf());
    std::cerr.rdbuf(null_stream.rdbuf());
    
    auto start = std::chrono::high_resolution_clock::now();
    auto result = julia_rur::solve_polynomial_system_enhanced(
        test.polynomials, test.variables, config
    );
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    
    // Restore stdout and stderr
    std::cout.rdbuf(orig_cout);
    std::cerr.rdbuf(orig_cerr);
    
    bool passed = result.success;
    if (test.expected_solutions >= 0 && result.solutions.size() != test.expected_solutions) {
        passed = false;
    }
    
    std::cout << std::setw(40) << std::left << test.name << " ";
    if (passed) {
        std::cout << "PASS";
    } else {
        std::cout << "FAIL";
    }
    std::cout << " (" << result.solutions.size() << " solutions, " 
              << std::fixed << std::setprecision(2) << elapsed << "s)" << std::endl;
    
    return passed;
}

int main() {
    std::vector<QuickTest> tests = {
        // Basic tests
        {"Univariate x^2 - 4", {"x^2 - 4"}, {"x"}, 2},
        {"Linear 2x2", {"x + y - 3", "2*x - y"}, {"x", "y"}, 1},
        {"Circle and line", {"x^2 + y^2 - 1", "y - x"}, {"x", "y"}, 2},
        
        // Julia package tests
        {"Julia example", 
         {"x0^2 + 2*x1^2 + 2*x2^2 - x0", "2*x0*x1 + 2*x1*x2 - x1", "x0 + 2*x1 + 2*x2 - 1"},
         {"x0", "x1", "x2"}, -1},
        
        // Univariate with multiplicity
        {"Double root (x-1)^2", {"(x-1)^2"}, {"x"}, 1},
        {"Triple root (x-2)^3", {"(x-2)^3"}, {"x"}, 1},
        
        // Higher degree
        {"Degree 10: x^10 - 1", {"x^10 - 1"}, {"x"}, 10},
        
        // Linear systems
        {"Linear 3x3", 
         {"x + y + z - 6", "2*x - y + z - 3", "x + 2*y - z - 1"},
         {"x", "y", "z"}, 1},
        
        {"Linear 5x5",
         {"x1 + x2 + x3 + x4 + x5 - 5",
          "x1 - x2 + x3 - x4 + x5 - 1",
          "x1 + 2*x2 + 3*x3 + 4*x4 + 5*x5 - 15",
          "2*x1 + 3*x2 + x3 + x4 + x5 - 8",
          "x1 + x2 + x3 + x4 - x5 - 2"},
         {"x1", "x2", "x3", "x4", "x5"}, 1},
        
        // Symmetric polynomial
        {"Elementary symmetric 3",
         {"x + y + z - 3", "x*y + y*z + z*x - 3", "x*y*z - 1"},
         {"x", "y", "z"}, 6},
         
        // Sparse high degree
        {"Sparse degree 20", {"x^20 + x^10 - x^5 + x - 1"}, {"x"}, -1},
    };
    
    std::cout << "\n===========================================" << std::endl;
    std::cout << "RUR C++ Test Suite - Quick Validation" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << std::endl;
    
    int passed = 0;
    int total = tests.size();
    
    for (const auto& test : tests) {
        if (run_quick_test(test)) {
            passed++;
        }
    }
    
    std::cout << "\n===========================================" << std::endl;
    std::cout << "Summary: " << passed << "/" << total << " tests passed" << std::endl;
    std::cout << "===========================================" << std::endl;
    
    return (passed == total) ? 0 : 1;
}