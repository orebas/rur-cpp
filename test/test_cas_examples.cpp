/**
 * @file test_cas_examples.cpp
 * @brief Test systems from various Computer Algebra Systems
 * 
 * Sources:
 * - Maple examples
 * - Mathematica examples  
 * - Macaulay2 examples
 * - SageMath examples
 * - Singular examples
 * - CoCoA examples
 */

#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include "../src/julia_rur/polynomial_evaluator.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace julia_rur;

struct CASExample {
    std::string name;
    std::string source;
    std::vector<std::string> polynomials;
    std::vector<std::string> variables;
    int expected_solutions;
    std::string notes;
};

void test_cas_example(const CASExample& example) {
    std::cout << "\nTesting: " << example.name << " (" << example.source << ")" << std::endl;
    if (!example.notes.empty()) {
        std::cout << "  Notes: " << example.notes << std::endl;
    }
    
    EnhancedSolverConfig config;
    config.verbose = false;
    
    auto start = std::chrono::high_resolution_clock::now();
    auto result = solve_polynomial_system_enhanced(example.polynomials, example.variables, config);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    
    if (!result.success) {
        std::cout << "  FAILED: " << result.error_message << std::endl;
        return;
    }
    
    std::cout << "  Solutions: " << result.solutions.size();
    if (example.expected_solutions >= 0) {
        std::cout << " (expected: " << example.expected_solutions << ")";
        if (result.solutions.size() != static_cast<size_t>(example.expected_solutions)) {
            std::cout << " MISMATCH!";
        }
    }
    std::cout << std::endl;
    std::cout << "  Time: " << std::fixed << std::setprecision(3) << elapsed << "s" << std::endl;
    
    // Compute max residual
    double max_residual = 0.0;
    for (const auto& sol : result.solutions) {
        double res = PolynomialEvaluator::compute_residual(example.polynomials, example.variables, sol);
        max_residual = std::max(max_residual, res);
    }
    std::cout << "  Max residual: " << std::scientific << std::setprecision(2) << max_residual << std::endl;
}

int main() {
    std::cout << "=============================================\n";
    std::cout << "Computer Algebra Systems Examples Test Suite\n";
    std::cout << "=============================================\n";
    
    std::vector<CASExample> examples = {
        // Maple examples
        {
            "Maple-RootFinding1",
            "Maple documentation",
            {"x^2 + y^2 - 4", "x*y - 1"},
            {"x", "y"},
            4,
            "Basic root finding example from Maple"
        },
        
        {
            "Maple-Groebner1",
            "Maple Groebner package",
            {"x^2 + x*y + 1", "y^2 - x + 2"},
            {"x", "y"},
            4,
            "Groebner basis example"
        },
        
        // Mathematica examples
        {
            "Mathematica-Solve1",
            "Mathematica documentation",
            {"x^2 + y^2 + z^2 - 1", "x + y + z", "x*y - z^2"},
            {"x", "y", "z"},
            4,
            "NSolve example from Mathematica"
        },
        
        {
            "Mathematica-Reduce",
            "Mathematica Reduce",
            {"x^3 - 2*x - 5", "y^2 - x"},
            {"x", "y"},
            6,
            "Cylindrical algebraic decomposition example"
        },
        
        // Macaulay2 examples
        {
            "Macaulay2-Tutorial",
            "Macaulay2 tutorial",
            {"x^3 - y^2", "x^2 - y - 1"},
            {"x", "y"},
            3,
            "Intersection of curves"
        },
        
        {
            "Macaulay2-Ideal",
            "Macaulay2 ideal theory",
            {"x^2*y - z^3", "x*z - y^2", "y*z - x^3"},
            {"x", "y", "z"},
            7,
            "Twisted cubic ideal"
        },
        
        // SageMath examples
        {
            "Sage-Variety1",
            "SageMath algebraic geometry",
            {"x^2 - y*z", "x*y - z^2", "y^2 - x*z"},
            {"x", "y", "z"},
            3,
            "Veronese surface equations"
        },
        
        {
            "Sage-Elliptic",
            "SageMath elliptic curves",
            {"y^2 - x^3 - x", "2*x + y - 3"},
            {"x", "y"},
            3,
            "Elliptic curve with line"
        },
        
        // Singular examples
        {
            "Singular-Example1",
            "Singular manual",
            {"x^2 + y^2 + z^2 - 1", "x*y + y*z + x*z", "x*y*z"},
            {"x", "y", "z"},
            8,
            "Standard Singular example"
        },
        
        {
            "Singular-Invariants",
            "Singular invariant theory",
            {"x^4 + y^4 - 1", "x^2*y^2 - 0.25"},
            {"x", "y"},
            8,
            "Invariant polynomial system"
        },
        
        // CoCoA examples
        {
            "CoCoA-Introduction",
            "CoCoA introduction",
            {"x^2 - 2", "y^2 - 3", "x + y - z"},
            {"x", "y", "z"},
            4,
            "Basic CoCoA example"
        },
        
        {
            "CoCoA-Buchberger",
            "CoCoA Buchberger",
            {"x^2*y + x*y^2 - 2*y", "x^3 - y^2 + 1"},
            {"x", "y"},
            5,
            "Buchberger algorithm example"
        },
        
        // GAP (computational algebra)
        {
            "GAP-Polynomial",
            "GAP system",
            {"x^3 + x + 1", "y^2 - x^2 + 1"},
            {"x", "y"},
            6,
            "GAP polynomial example"
        },
        
        // Magma examples
        {
            "Magma-Scheme",
            "Magma handbook",
            {"x*y - z^2", "x^2 - y*z", "y^2 - x*z"},
            {"x", "y", "z"},
            3,
            "Scheme intersection in Magma"
        },
        
        // Reduce examples
        {
            "Reduce-Solve",
            "Reduce CAS",
            {"x^2 + y^2 - 5", "x*y - 2"},
            {"x", "y"},
            4,
            "Reduce solve example"
        },
        
        // Axiom/FriCAS examples
        {
            "Axiom-Algebraic",
            "Axiom/FriCAS",
            {"x^3 - 3*x + 1", "y^2 - x - 1"},
            {"x", "y"},
            6,
            "Axiom algebraic example"
        },
        
        // SymPy examples
        {
            "SymPy-Nonlinear",
            "SymPy Python",
            {"x^2 + y^2 - 1", "x^2 - y^2 - 0.5"},
            {"x", "y"},
            4,
            "SymPy nonlinear system"
        },
        
        {
            "SymPy-Polynomial",
            "SymPy polys module",
            {"x^3 + 2*x^2 - x - 2", "y - x^2"},
            {"x", "y"},
            3,
            "SymPy polynomial module"
        },
        
        // Maxima examples
        {
            "Maxima-Algsys",
            "Maxima algsys",
            {"x^2 + y^2 - 2", "x - y^2"},
            {"x", "y"},
            3,
            "Maxima algebraic system solver"
        },
        
        // Pari/GP examples
        {
            "PariGP-Polroots",
            "Pari/GP",
            {"x^4 - 10*x^2 + 1", "y - x^2 + 2"},
            {"x", "y"},
            4,
            "Pari/GP polynomial roots"
        },
        
        // OSCAR examples
        {
            "OSCAR-Ideal",
            "OSCAR system",
            {"x^2 + y^2 + z^2 - 3", "x + y + z - 3", "x*y + y*z + x*z - 3"},
            {"x", "y", "z"},
            8,
            "OSCAR ideal computation"
        },
        
        // Classic textbook examples
        {
            "CLO-Example",
            "Cox-Little-O'Shea",
            {"x^2 + y^2 - 1", "x^2 + y^2 + (z-2)^2 - 4", "z"},
            {"x", "y", "z"},
            2,
            "Ideals, Varieties, and Algorithms example"
        },
        
        {
            "Eisenbud-Example",
            "Eisenbud textbook",
            {"x^2 - y*z", "x*y - z^2", "x*z - y^2"},
            {"x", "y", "z"},
            3,
            "Commutative Algebra example"
        },
        
        // Robotics applications
        {
            "Robot-3R",
            "Robotics",
            {"x^2 + y^2 - 2", "(x-1)^2 + y^2 - 1", "x + y - 1.5"},
            {"x", "y"},
            2,
            "3R robot workspace boundary"
        },
        
        // Control theory
        {
            "Control-Stabilization",
            "Control theory",
            {"x^2 - y", "y^2 - z", "z^2 - x"},
            {"x", "y", "z"},
            7,
            "Stabilization polynomial"
        },
        
        // Optimization critical points
        {
            "Optimization-Lagrange",
            "Optimization",
            {"2*x + y - 3", "x + 2*y - 3", "x^2 + y^2 - 2"},
            {"x", "y"},
            2,
            "Lagrange multiplier system"
        },
        
        // Computer graphics
        {
            "Graphics-Bezier",
            "Computer graphics",
            {"x^3 - 3*x + 2", "y - x^2 + x"},
            {"x", "y"},
            3,
            "Bezier curve intersection"
        },
        
        // Cryptography
        {
            "Crypto-EC",
            "Cryptography",
            {"y^2 - x^3 - 7", "2*x + 3*y - 12"},
            {"x", "y"},
            3,
            "Elliptic curve cryptography"
        },
        
        // Physics applications
        {
            "Physics-Pendulum",
            "Physics",
            {"x^2 + y^2 - 1", "x - 2*y^2 + 0.5"},
            {"x", "y"},
            3,
            "Double pendulum equilibrium"
        },
        
        // Differential equations
        {
            "ODE-Equilibrium",
            "Differential equations",
            {"x*(1 - x - y)", "y*(1 - 1.5*x - 0.5*y)"},
            {"x", "y"},
            4,
            "Lotka-Volterra equilibrium"
        }
    };
    
    int passed = 0;
    int failed = 0;
    
    for (const auto& example : examples) {
        try {
            test_cas_example(example);
            passed++;
        } catch (const std::exception& e) {
            std::cerr << "  Exception: " << e.what() << std::endl;
            failed++;
        }
    }
    
    std::cout << "\n=============================================\n";
    std::cout << "Results: " << passed << " passed, " << failed << " failed\n";
    std::cout << "Total systems tested: " << (passed + failed) << "\n";
    std::cout << "=============================================\n";
    
    return failed == 0 ? 0 : 1;
}