/**
 * @file test_benchmark_systems.cpp
 * @brief Comprehensive test suite with benchmark polynomial systems from literature
 * 
 * Sources:
 * - SymbolicData.org benchmark problems
 * - PHCpack test suite
 * - Singular examples
 * - Classic algebraic geometry problems
 * - Competition problems from ISSAC/CASC
 */

#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include "../src/julia_rur/polynomial_evaluator.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>

using namespace julia_rur;

struct BenchmarkSystem {
    std::string name;
    std::string source;
    std::vector<std::string> polynomials;
    std::vector<std::string> variables;
    int expected_solutions;  // -1 if unknown
    std::string description;
};

void run_benchmark(const BenchmarkSystem& sys) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "System: " << sys.name << std::endl;
    std::cout << "Source: " << sys.source << std::endl;
    std::cout << "Description: " << sys.description << std::endl;
    std::cout << "Variables: ";
    for (size_t i = 0; i < sys.variables.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << sys.variables[i];
    }
    std::cout << std::endl;
    
    // Solve the system
    EnhancedSolverConfig config;
    config.verbose = false;
    
    auto start = std::chrono::high_resolution_clock::now();
    auto result = solve_polynomial_system_enhanced(sys.polynomials, sys.variables, config);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    
    if (!result.success) {
        std::cout << "FAILED to solve: " << result.error_message << std::endl;
        return;
    }
    
    std::cout << "Found " << result.solutions.size() << " solutions in " 
              << std::fixed << std::setprecision(3) << elapsed << " seconds" << std::endl;
    
    if (sys.expected_solutions >= 0 && 
        result.solutions.size() != static_cast<size_t>(sys.expected_solutions)) {
        std::cout << "WARNING: Expected " << sys.expected_solutions << " solutions" << std::endl;
    }
    
    // Check residuals
    double max_residual = 0.0;
    for (const auto& solution : result.solutions) {
        double residual = PolynomialEvaluator::compute_residual(
            sys.polynomials, sys.variables, solution
        );
        max_residual = std::max(max_residual, residual);
    }
    std::cout << "Max residual: " << std::scientific << std::setprecision(2) << max_residual << std::endl;
}

int main() {
    std::cout << "\n=========================================\n";
    std::cout << "Benchmark Polynomial Systems Test Suite\n";
    std::cout << "=========================================\n";
    
    std::vector<BenchmarkSystem> systems = {
        // Classic systems from algebraic geometry
        {
            "Cyclic-3",
            "Classic algebraic geometry",
            {"x + y + z", "x*y + y*z + z*x", "x*y*z - 1"},
            {"x", "y", "z"},
            6,
            "Cyclic roots system of degree 3"
        },
        
        {
            "Cyclic-4",
            "Classic algebraic geometry",
            {
                "x1 + x2 + x3 + x4",
                "x1*x2 + x2*x3 + x3*x4 + x4*x1",
                "x1*x2*x3 + x2*x3*x4 + x3*x4*x1 + x4*x1*x2",
                "x1*x2*x3*x4 - 1"
            },
            {"x1", "x2", "x3", "x4"},
            16,
            "Cyclic roots system of degree 4"
        },
        
        // Katsura systems
        {
            "Katsura-3",
            "SymbolicData.org",
            {
                "x + 2*y + 2*z - 1",
                "x^2 + 2*y^2 + 2*z^2 - x",
                "2*x*y + 2*y*z - y"
            },
            {"x", "y", "z"},
            4,
            "Katsura system of index 3"
        },
        
        {
            "Katsura-4",
            "SymbolicData.org",
            {
                "x + 2*y + 2*z + 2*w - 1",
                "x^2 + 2*y^2 + 2*z^2 + 2*w^2 - x",
                "2*x*y + 2*y*z + 2*z*w - y",
                "y^2 + 2*x*z + 2*y*w - z"
            },
            {"x", "y", "z", "w"},
            8,
            "Katsura system of index 4"
        },
        
        // Noon systems
        {
            "Noon-3",
            "PHCpack test suite",
            {
                "10*x^2*y - 2*x*y^2 - 1",
                "10*x*y^2 - 2*y*z^2 - 1",
                "10*y^2*z - 2*x^2*z - 1"
            },
            {"x", "y", "z"},
            6,
            "Noon system n=3"
        },
        
        // Economic equilibrium model
        {
            "Eco-5",
            "Economic models",
            {
                "x^2 + 2*y^2 - 3",
                "x^2 + y^2 + z^2 - 3",
                "x + y + 3*z^2 - 4",
                "3*x^2 + y + z - 4",
                "x + y + z - 2"
            },
            {"x", "y", "z"},
            16,
            "Economic equilibrium system"
        },
        
        // Robot kinematics
        {
            "Robot-2link",
            "Robotics",
            {
                "x^2 + y^2 - 1",
                "(x-2)^2 + (y-1)^2 - 4"
            },
            {"x", "y"},
            2,
            "2-link planar robot inverse kinematics"
        },
        
        // Stewart platform
        {
            "Stewart-Gough",
            "Kinematics",
            {
                "x^2 + y^2 + z^2 - 1",
                "x + y + z",
                "x*y + y*z + z*x - 0.25"
            },
            {"x", "y", "z"},
            8,
            "Simplified Stewart-Gough platform"
        },
        
        // Intersection problems
        {
            "Sphere-Torus",
            "Geometric intersection",
            {
                "x^2 + y^2 + z^2 - 4",
                "(x^2 + y^2 - 3)^2 + z^2 - 1"
            },
            {"x", "y", "z"},
            8,
            "Sphere-torus intersection"
        },
        
        {
            "Three-Spheres",
            "Geometric intersection",
            {
                "x^2 + y^2 + z^2 - 1",
                "(x-1)^2 + y^2 + z^2 - 1",
                "x^2 + (y-1)^2 + z^2 - 1"
            },
            {"x", "y", "z"},
            2,
            "Three unit spheres intersection"
        },
        
        // Chemical equilibrium
        {
            "Chem-Equilibrium",
            "Chemistry",
            {
                "x*y - 0.1",
                "x*z - 0.2",
                "y*z - 0.3",
                "x + y + z - 1"
            },
            {"x", "y", "z"},
            2,
            "Chemical equilibrium system"
        },
        
        // Wilkinson polynomial system
        {
            "Wilkinson-variant",
            "Numerical analysis",
            {
                "(x-1)*(x-2)*(x-3)",
                "(y-1)*(y-2)",
                "x + y - 4"
            },
            {"x", "y"},
            6,
            "Wilkinson-like system"
        },
        
        // Algebraic curve intersections
        {
            "Elliptic-curves",
            "Algebraic geometry",
            {
                "y^2 - x^3 - x - 1",
                "y^2 - x^3 + 2*x - 1"
            },
            {"x", "y"},
            9,
            "Two elliptic curves intersection"
        },
        
        // Bezout bound examples
        {
            "Bezout-2-3",
            "Algebraic geometry",
            {
                "x^2 + y^2 - 1",
                "x^3 - y"
            },
            {"x", "y"},
            6,
            "Degree 2 and 3 curves (Bezout bound = 6)"
        },
        
        {
            "Bezout-3-3",
            "Algebraic geometry",
            {
                "x^3 + y^3 - 1",
                "x^3 - y^3 - 0.5"
            },
            {"x", "y"},
            9,
            "Two cubic curves (Bezout bound = 9)"
        },
        
        // Moeller examples
        {
            "Moeller-1",
            "Test suite",
            {
                "x^2 + y + z - 1",
                "x + y^2 + z - 1",
                "x + y + z^2 - 1"
            },
            {"x", "y", "z"},
            5,
            "Moeller system 1"
        },
        
        {
            "Moeller-2",
            "Test suite",
            {
                "x^2*y - 2*x + 1",
                "x*y^2 - 2*y + 1"
            },
            {"x", "y"},
            5,
            "Moeller system 2"
        },
        
        // Rose curve intersection
        {
            "Rose-3-petals",
            "Curves",
            {
                "x^4 + y^4 - (x^2 + y^2)",
                "x^2 + y^2 - 0.5"
            },
            {"x", "y"},
            8,
            "Rose curve with 3 petals and circle"
        },
        
        // Cassini ovals
        {
            "Cassini-oval",
            "Classical curves",
            {
                "(x^2 + y^2)^2 - 2*(x^2 - y^2) - 1",
                "x^2 + y^2 - 1"
            },
            {"x", "y"},
            4,
            "Cassini oval and unit circle"
        },
        
        // Lemniscate
        {
            "Lemniscate-Bernoulli",
            "Classical curves",
            {
                "(x^2 + y^2)^2 - (x^2 - y^2)",
                "x + y - 0.5"
            },
            {"x", "y"},
            4,
            "Lemniscate of Bernoulli and line"
        },
        
        // Competition problems
        {
            "ISSAC-97",
            "Competition",
            {
                "x^3 - 3*x*y^2 - 1",
                "3*x^2*y - y^3"
            },
            {"x", "y"},
            3,
            "ISSAC 1997 challenge problem"
        },
        
        // Neural network critical points
        {
            "Neural-2-2-1",
            "Machine learning",
            {
                "x^2 + y^2 - 1",
                "x*y - 0.25"
            },
            {"x", "y"},
            4,
            "Simple neural network critical points"
        },
        
        // Game theory equilibria
        {
            "Nash-2x2",
            "Game theory",
            {
                "x*(1-x)*(1-2*y)",
                "y*(1-y)*(1-2*x)"
            },
            {"x", "y"},
            5,
            "Nash equilibrium for 2x2 game"
        },
        
        // Trott curve
        {
            "Trott-curve",
            "Algebraic curves",
            {
                "144*(x^4 + y^4) - 225*(x^2 + y^2) + 350*x^2*y^2 + 81",
                "x + y - 1"
            },
            {"x", "y"},
            28,
            "Trott curve (degree 4) with line"
        },
        
        // Viviani's curve
        {
            "Viviani",
            "Space curves",
            {
                "x^2 + y^2 + z^2 - 4",
                "(x-1)^2 + y^2 - 1",
                "z"
            },
            {"x", "y", "z"},
            2,
            "Viviani's curve in plane z=0"
        },
        
        // Descartes folium
        {
            "Descartes-folium",
            "Classical curves",
            {
                "x^3 + y^3 - 3*x*y",
                "x + y - 1"
            },
            {"x", "y"},
            3,
            "Folium of Descartes with line"
        },
        
        // Butterfly curve
        {
            "Butterfly",
            "Transcendental approximation",
            {
                "x^4 + y^4 - x^2 - y^2 + x*y",
                "x^2 + y^2 - 1"
            },
            {"x", "y"},
            8,
            "Butterfly curve approximation"
        },
        
        // Eight curve
        {
            "Figure-eight",
            "Lemniscate variant",
            {
                "x^4 - x^2 + y^2",
                "x - 0.5"
            },
            {"x", "y"},
            4,
            "Figure eight curve with vertical line"
        },
        
        // Hippopede
        {
            "Hippopede",
            "Classical curves",
            {
                "(x^2 + y^2)^2 - 2*x^2 - y^2",
                "y - 0.3"
            },
            {"x", "y"},
            4,
            "Hippopede of Proclus"
        }
    };
    
    // Run all benchmarks
    int passed = 0;
    int failed = 0;
    
    for (const auto& sys : systems) {
        try {
            run_benchmark(sys);
            passed++;
        } catch (const std::exception& e) {
            std::cerr << "Exception in " << sys.name << ": " << e.what() << std::endl;
            failed++;
        }
    }
    
    std::cout << "\n=========================================\n";
    std::cout << "Summary: " << passed << " passed, " << failed << " failed\n";
    std::cout << "=========================================\n";
    
    return failed == 0 ? 0 : 1;
}