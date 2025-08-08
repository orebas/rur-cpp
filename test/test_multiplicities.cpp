#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include "../src/julia_rur/polynomial_evaluator.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace julia_rur;

struct MultiplicityTest {
    std::string name;
    std::vector<std::string> polynomials;
    std::vector<std::string> variables;
    int expected_distinct_solutions;
    std::string description;
};

void run_multiplicity_test(const MultiplicityTest& test) {
    std::cout << "\n" << test.name << ":" << std::endl;
    std::cout << "  Description: " << test.description << std::endl;
    std::cout << "  System: ";
    for (size_t i = 0; i < test.polynomials.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << test.polynomials[i] << " = 0";
    }
    std::cout << std::endl;
    
    // Solve the system
    EnhancedSolverConfig config;
    config.verbose = false;
    
    auto start = std::chrono::high_resolution_clock::now();
    auto result = solve_polynomial_system_enhanced(test.polynomials, test.variables, config);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    
    if (!result.success) {
        std::cout << "  FAILED to solve: " << result.error_message << std::endl;
        return;
    }
    
    std::cout << "  Found " << result.solutions.size() << " solutions in " 
              << std::fixed << std::setprecision(3) << elapsed << " seconds" << std::endl;
    
    // Verify each solution
    double max_residual = 0.0;
    for (size_t i = 0; i < result.solutions.size(); ++i) {
        double residual = PolynomialEvaluator::compute_residual(
            test.polynomials, test.variables, result.solutions[i]
        );
        max_residual = std::max(max_residual, residual);
        
        std::cout << "  Solution " << (i+1) << ": ";
        for (size_t j = 0; j < test.variables.size(); ++j) {
            if (j > 0) std::cout << ", ";
            std::cout << test.variables[j] << " = " << result.solutions[i][j];
        }
        std::cout << " (residual: " << std::scientific << std::setprecision(2) << residual << ")" << std::endl;
    }
    
    std::cout << "  Max residual: " << std::scientific << std::setprecision(2) << max_residual << std::endl;
    
    // Check expected count
    if (test.expected_distinct_solutions >= 0 && 
        result.solutions.size() != static_cast<size_t>(test.expected_distinct_solutions)) {
        std::cout << "  WARNING: Expected " << test.expected_distinct_solutions 
                  << " distinct solutions, got " << result.solutions.size() << std::endl;
    }
}

int main() {
    std::cout << "\n=========================================\n";
    std::cout << "Testing Systems with Multiplicities\n";
    std::cout << "=========================================\n";
    
    std::vector<MultiplicityTest> tests = {
        // Double roots
        {
            "Double root univariate",
            {"(x-1)^2"},
            {"x"},
            1,
            "Single root x=1 with multiplicity 2"
        },
        
        {
            "Double root in 2D",
            {"(x-1)^2", "y"},
            {"x", "y"},
            1,
            "Point (1,0) with multiplicity 2"
        },
        
        // Triple roots
        {
            "Triple root univariate",
            {"(x-2)^3"},
            {"x"},
            1,
            "Single root x=2 with multiplicity 3"
        },
        
        // Tangent multiplicities
        {
            "Tangent circle and line",
            {"x^2 + y^2 - 1", "y - 1"},
            {"x", "y"},
            1,
            "Line y=1 is tangent to unit circle at (0,1)"
        },
        
        {
            "Tangent parabola and line",
            {"y - x^2", "y"},
            {"x", "y"},
            1,
            "Line y=0 is tangent to parabola at origin"
        },
        
        // Self-intersection
        {
            "Figure-eight curve",
            {"x^2 - x^4 - y^2", "y"},
            {"x", "y"},
            1,
            "Self-intersection at origin"
        },
        
        // Multiple double roots
        {
            "Product of squares",
            {"(x-1)^2 * (x+1)^2"},
            {"x"},
            2,
            "Two double roots at x=1 and x=-1"
        },
        
        // Degenerate systems
        {
            "Repeated equation",
            {"x^2 - 1", "x^2 - 1"},
            {"x"},
            2,
            "Same equation twice, roots at x=±1"
        },
        
        // High multiplicity
        {
            "Fourth power",
            {"(x-3)^4"},
            {"x"},
            1,
            "Single root x=3 with multiplicity 4"
        },
        
        // Mixed multiplicities
        {
            "Mixed multiplicities",
            {"x^2 * (x-1)^3"},
            {"x"},
            2,
            "Double root at x=0, triple root at x=1"
        },
        
        // Tangent surfaces in 3D
        {
            "Tangent sphere and plane",
            {"x^2 + y^2 + z^2 - 1", "z - 1", "x"},
            {"x", "y", "z"},
            1,
            "Plane z=1 tangent to unit sphere at (0,0,1)"
        },
        
        // Cusp
        {
            "Cusp curve",
            {"y^2 - x^3", "y"},
            {"x", "y"},
            1,
            "Cusp at origin"
        },
        
        // Near-multiplicities (clustered roots)
        {
            "Nearly double root",
            {"(x-1)*(x-1.001)"},
            {"x"},
            2,
            "Two very close roots that should be distinguished"
        },
        
        // Symbolic multiplicities
        {
            "Wilkinson-like polynomial",
            {"(x-1)*(x-2)*(x-3)*(x-4)*(x-5)"},
            {"x"},
            5,
            "Five distinct roots, no multiplicities"
        }
    };
    
    // Run all tests
    for (const auto& test : tests) {
        run_multiplicity_test(test);
    }
    
    std::cout << "\n=========================================\n";
    std::cout << "Testing Complete\n";
    std::cout << "=========================================\n";
    
    return 0;
}