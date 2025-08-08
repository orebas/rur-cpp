#include "../src/julia_rur/hyperplane_sections.hpp"
#include <iostream>
#include <iomanip>

using namespace julia_rur;

void test_dimension_detection() {
    std::cout << "=== Testing Dimension Detection ===" << std::endl;
    
    // Test 1: Zero-dimensional system (2 equations, 2 variables, finite solutions)
    {
        std::cout << "\n1. Zero-dimensional: x^2 + y^2 - 1 = 0, y - x = 0" << std::endl;
        std::vector<std::string> polys = {"x^2 + y^2 - 1", "y - x"};
        std::vector<std::string> vars = {"x", "y"};
        
        auto dim_info = analyze_system_dimension(polys, vars);
        std::cout << "   Dimension: " << dim_info.dimension << " (expected: 0)" << std::endl;
        std::cout << "   Is zero-dimensional: " << (dim_info.is_zero_dimensional ? "yes" : "no") << std::endl;
        std::cout << "   Free variables: " << dim_info.free_variables.size() << std::endl;
    }
    
    // Test 2: One-dimensional system (circle)
    {
        std::cout << "\n2. One-dimensional: x^2 + y^2 - 1 = 0" << std::endl;
        std::vector<std::string> polys = {"x^2 + y^2 - 1"};
        std::vector<std::string> vars = {"x", "y"};
        
        auto dim_info = analyze_system_dimension(polys, vars);
        std::cout << "   Dimension: " << dim_info.dimension << " (expected: 1)" << std::endl;
        std::cout << "   Codimension: " << dim_info.codimension << " (hyperplanes needed)" << std::endl;
        std::cout << "   Free variables: ";
        for (int idx : dim_info.free_variables) {
            std::cout << vars[idx] << " ";
        }
        std::cout << std::endl;
    }
    
    // Test 3: Two-dimensional system (single equation in 3 variables)
    {
        std::cout << "\n3. Two-dimensional: x + y + z - 1 = 0" << std::endl;
        std::vector<std::string> polys = {"x + y + z - 1"};
        std::vector<std::string> vars = {"x", "y", "z"};
        
        auto dim_info = analyze_system_dimension(polys, vars);
        std::cout << "   Dimension: " << dim_info.dimension << " (expected: 2)" << std::endl;
        std::cout << "   Codimension: " << dim_info.codimension << " (hyperplanes needed)" << std::endl;
    }
}

void test_hyperplane_generation() {
    std::cout << "\n\n=== Testing Hyperplane Generation ===" << std::endl;
    
    std::vector<std::string> vars = {"x", "y", "z"};
    HyperplaneSectionConfig config;
    config.min_coefficient = -10;
    config.max_coefficient = 10;
    
    std::random_device rd;
    std::mt19937 rng(rd());
    
    std::cout << "\nGenerating 5 random hyperplanes:" << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::string hyperplane = generate_random_hyperplane(vars, config, rng);
        std::cout << "  " << (i+1) << ". " << hyperplane << " = 0" << std::endl;
    }
}

void test_circle_with_hyperplane() {
    std::cout << "\n\n=== Testing Circle with Hyperplane Section ===" << std::endl;
    
    // Circle: x^2 + y^2 = 1 (1-dimensional)
    std::vector<std::string> polys = {"x^2 + y^2 - 1"};
    std::vector<std::string> vars = {"x", "y"};
    
    HyperplaneSectionConfig config;
    config.verbose = true;
    config.min_coefficient = -5;
    config.max_coefficient = 5;
    
    auto result = solve_via_hyperplane_sections(polys, vars, config);
    
    if (result.success) {
        std::cout << "\nSuccess! Found " << result.solutions.solutions.size() << " intersection points" << std::endl;
        std::cout << "Hyperplane added: " << result.hyperplanes[0] << " = 0" << std::endl;
        
        std::cout << "\nIntersection points:" << std::endl;
        for (size_t i = 0; i < result.solutions.solutions.size(); ++i) {
            std::cout << "  " << (i+1) << ". ";
            for (size_t j = 0; j < vars.size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << vars[j] << " = " << std::setprecision(6) 
                          << result.solutions.solutions[i][j].real();
                if (std::abs(result.solutions.solutions[i][j].imag()) > 1e-10) {
                    std::cout << " + " << result.solutions.solutions[i][j].imag() << "i";
                }
            }
            std::cout << std::endl;
        }
        
        // Verify points are on the circle
        std::cout << "\nVerification (x^2 + y^2 should equal 1):" << std::endl;
        for (size_t i = 0; i < result.solutions.solutions.size(); ++i) {
            auto x = result.solutions.solutions[i][0].real();
            auto y = result.solutions.solutions[i][1].real();
            double sum = x*x + y*y;
            std::cout << "  Point " << (i+1) << ": " << sum << std::endl;
        }
    } else {
        std::cout << "\nFailed: " << result.error_message << std::endl;
    }
}

void test_surface_sampling() {
    std::cout << "\n\n=== Testing Surface Sampling (Paraboloid) ===" << std::endl;
    
    // Paraboloid: z = x^2 + y^2 (2-dimensional surface in 3D)
    std::vector<std::string> polys = {"z - x^2 - y^2"};
    std::vector<std::string> vars = {"x", "y", "z"};
    
    HyperplaneSectionConfig config;
    config.verbose = false;
    config.min_coefficient = -3;
    config.max_coefficient = 3;
    
    std::cout << "\nSampling 3 different hyperplane sections:" << std::endl;
    auto samples = sample_variety_points(polys, vars, 3, config);
    
    for (size_t s = 0; s < samples.size(); ++s) {
        if (samples[s].success) {
            std::cout << "\nSample " << (s+1) << ":" << std::endl;
            std::cout << "  Hyperplanes: ";
            for (const auto& h : samples[s].hyperplanes) {
                std::cout << h << " = 0; ";
            }
            std::cout << std::endl;
            std::cout << "  Points found: " << samples[s].solutions.solutions.size() << std::endl;
            
            // Show first 2 points
            for (size_t i = 0; i < std::min(size_t(2), samples[s].solutions.solutions.size()); ++i) {
                std::cout << "    (";
                for (size_t j = 0; j < vars.size(); ++j) {
                    if (j > 0) std::cout << ", ";
                    std::cout << std::setprecision(3) << samples[s].solutions.solutions[i][j].real();
                }
                std::cout << ")" << std::endl;
            }
        }
    }
}

void test_adaptive_solver() {
    std::cout << "\n\n=== Testing Adaptive Solver ===" << std::endl;
    
    // Test both zero and positive dimensional systems
    std::cout << "\n1. Adaptive solver on zero-dimensional system:" << std::endl;
    {
        std::vector<std::string> polys = {"x^2 - 2", "y^2 - 3"};
        std::vector<std::string> vars = {"x", "y"};
        
        auto result = solve_polynomial_system_adaptive(polys, vars);
        std::cout << "   Dimension detected: " << result.dimension_info.dimension << std::endl;
        std::cout << "   Hyperplanes added: " << result.hyperplanes.size() << std::endl;
        std::cout << "   Solutions found: " << result.solutions.solutions.size() << std::endl;
    }
    
    std::cout << "\n2. Adaptive solver on positive-dimensional system:" << std::endl;
    {
        std::vector<std::string> polys = {"x^2 + y^2 + z^2 - 1"};  // Sphere
        std::vector<std::string> vars = {"x", "y", "z"};
        
        HyperplaneSectionConfig config;
        config.verbose = true;
        
        auto result = solve_polynomial_system_adaptive(polys, vars, config);
        std::cout << "   Dimension detected: " << result.dimension_info.dimension << std::endl;
        std::cout << "   Hyperplanes added: " << result.hyperplanes.size() << std::endl;
        std::cout << "   Solutions found: " << result.solutions.solutions.size() << std::endl;
    }
}

int main() {
    test_dimension_detection();
    test_hyperplane_generation();
    test_circle_with_hyperplane();
    test_surface_sampling();
    test_adaptive_solver();
    
    std::cout << "\n=== All tests completed ===" << std::endl;
    return 0;
}