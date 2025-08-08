#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include <iostream>
#include <iomanip>

using namespace julia_rur;

void test_complex_roots() {
    std::cout << "=== Testing Complex Roots with Enhanced Solver ===" << std::endl;
    
    // Test 1: Circle intersected with line (zero-dimensional)
    // x^2 + y^2 = 1, y = x gives 2 real points
    {
        std::cout << "\n1. Circle with line: x^2 + y^2 - 1 = 0, y - x = 0" << std::endl;
        std::vector<std::string> polys = {"x^2 + y^2 - 1", "y - x"};
        std::vector<std::string> vars = {"x", "y"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        config.root_method = RootFindingMethod::AUTO;
        
        auto solution = solve_polynomial_system_enhanced(polys, vars, config);
        print_enhanced_solution(solution);
    }
    
    // Test 2: Complex univariate (x^4 + 1)
    // Has 4 complex roots on unit circle at 45, 135, 225, 315 degrees
    {
        std::cout << "\n2. Complex polynomial: x^4 + 1 = 0" << std::endl;
        std::vector<std::string> polys = {"x^4 + 1"};
        std::vector<std::string> vars = {"x"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        config.root_method = RootFindingMethod::FLINT;  // Force FLINT for certified roots
        config.separate_conjugate_pairs = true;
        
        auto solution = solve_polynomial_system_enhanced(polys, vars, config);
        print_enhanced_solution(solution);
    }
    
    // Test 3: Mixed real/complex system
    // x^2 + y^2 = 2, x*y = 1
    // Solutions: (1,1), (-1,-1) [real], and two complex conjugate pairs
    {
        std::cout << "\n3. Mixed system: x^2 + y^2 - 2 = 0, x*y - 1 = 0" << std::endl;
        std::vector<std::string> polys = {"x^2 + y^2 - 2", "x*y - 1"};
        std::vector<std::string> vars = {"x", "y"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        config.prefer_certified_roots = true;
        
        auto solution = solve_polynomial_system_enhanced(polys, vars, config);
        print_enhanced_solution(solution);
    }
    
    // Test 4: High degree polynomial with clustering
    // Wilkinson polynomial - product of (x-k) for k=1 to 6
    {
        std::cout << "\n4. Wilkinson-like: (x-1)(x-2)(x-3)(x-4)(x-5)(x-6)" << std::endl;
        // Expanded: x^6 - 21x^5 + 175x^4 - 735x^3 + 1624x^2 - 1764x + 720
        std::vector<std::string> polys = {
            "x^6 - 21*x^5 + 175*x^4 - 735*x^3 + 1624*x^2 - 1764*x + 720"
        };
        std::vector<std::string> vars = {"x"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        config.degree_threshold_for_flint = 5;  // Use FLINT for better accuracy
        
        auto solution = solve_polynomial_system_enhanced(polys, vars, config);
        print_enhanced_solution(solution);
    }
    
    // Test 5: Filter for real roots only
    {
        std::cout << "\n5. Same as test 2 but real roots only: x^4 + 1 = 0" << std::endl;
        std::vector<std::string> polys = {"x^4 + 1"};
        std::vector<std::string> vars = {"x"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        config.real_roots_only = true;  // Filter out complex solutions
        
        auto solution = solve_polynomial_system_enhanced(polys, vars, config);
        print_enhanced_solution(solution);
    }
    
    // Test 6: Nearly degenerate system to test root clustering
    {
        std::cout << "\n6. Nearly degenerate: (x-1)^3 * (x-1.0001)" << std::endl;
        // This has a triple root at x=1 and a simple root very close by
        std::vector<std::string> polys = {
            "x^4 - 4.0003*x^3 + 6.0009*x^2 - 4.0009*x + 1.0003"
        };
        std::vector<std::string> vars = {"x"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        config.clustering_threshold = 1e-3;  // Cluster roots within 0.001
        config.root_method = RootFindingMethod::FLINT;
        
        auto solution = solve_polynomial_system_enhanced(polys, vars, config);
        print_enhanced_solution(solution);
    }
}

void test_method_selection() {
    std::cout << "\n\n=== Testing Automatic Method Selection ===" << std::endl;
    
    // Small degree - should use Eigen
    {
        std::cout << "\n1. Small polynomial (degree 3) - expect Eigen" << std::endl;
        std::vector<std::string> polys = {"x^3 - 2*x + 1"};
        std::vector<std::string> vars = {"x"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        config.root_method = RootFindingMethod::AUTO;
        
        auto solution = solve_polynomial_system_enhanced(polys, vars, config);
        std::cout << "Method used: ";
        switch (solution.method_used) {
            case RootFindingMethod::EIGEN: std::cout << "Eigen"; break;
            case RootFindingMethod::FLINT: std::cout << "FLINT"; break;
            default: std::cout << "Unknown"; break;
        }
        std::cout << " (degree = 3)" << std::endl;
    }
    
    // High degree - should use FLINT
    {
        std::cout << "\n2. High degree polynomial (degree 25) - expect FLINT" << std::endl;
        // Create a simple high-degree polynomial: x^25 - 1
        std::string poly = "x^25 - 1";
        std::vector<std::string> polys = {poly};
        std::vector<std::string> vars = {"x"};
        
        EnhancedSolverConfig config;
        config.verbose = false;
        config.root_method = RootFindingMethod::AUTO;
        config.degree_threshold_for_flint = 20;
        
        auto solution = solve_polynomial_system_enhanced(polys, vars, config);
        std::cout << "Method used: ";
        switch (solution.method_used) {
            case RootFindingMethod::EIGEN: std::cout << "Eigen"; break;
            case RootFindingMethod::FLINT: std::cout << "FLINT"; break;
            default: std::cout << "Unknown"; break;
        }
        std::cout << " (degree = 25)" << std::endl;
        std::cout << "Found " << solution.solutions.size() << " roots" << std::endl;
    }
}

int main() {
    test_complex_roots();
    test_method_selection();
    
    std::cout << "\n=== All tests completed ===" << std::endl;
    return 0;
}