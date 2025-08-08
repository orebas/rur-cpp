#include "test_framework.hpp"
#include <filesystem>

using namespace test_framework;

int main() {
    // Create a small but comprehensive test suite
    std::vector<TestCase> tests;
    
    // 1. Basic univariate tests
    tests.push_back({
        "Univariate quadratic",
        {"x^2 - 4"},
        {"x"},
        2,  // ±2
        1e-8
    });
    
    tests.push_back({
        "Univariate degree 10",
        {"x^10 - 1"},
        {"x"},
        10,  // 10th roots of unity
        1e-8
    });
    
    // 2. Linear systems
    tests.push_back({
        "Linear 2x2",
        {"x + y - 3", "2*x - y"},
        {"x", "y"},
        1,
        1e-8
    });
    
    tests.push_back({
        "Linear 5x5",
        {
            "x1 + x2 + x3 + x4 + x5 - 5",
            "x1 - x2 + x3 - x4 + x5 - 1",
            "x1 + 2*x2 + 3*x3 + 4*x4 + 5*x5 - 15",
            "2*x1 + 3*x2 + x3 + x4 + x5 - 8",
            "x1 + x2 + x3 + x4 - x5 - 2"
        },
        {"x1", "x2", "x3", "x4", "x5"},
        1,
        1e-8
    });
    
    // 3. Nonlinear systems
    tests.push_back({
        "Circle-line intersection",
        {"x^2 + y^2 - 1", "y - x"},
        {"x", "y"},
        2,
        1e-8
    });
    
    tests.push_back({
        "Squared linear 5x5 (first 3 equations only)",
        {
            "x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 5",
            "x1^2 - x2^2 + x3^2 - x4^2 + x5^2 - 1",
            "x1^2 + 2*x2^2 + 3*x3^2 + 4*x4^2 + 5*x5^2 - 15"
        },
        {"x1", "x2", "x3", "x4", "x5"},
        -1,  // Unknown - positive dimensional
        1e-8
    });
    
    // 4. System with multiplicity
    tests.push_back({
        "Double root system",
        {"(x-1)^2", "y"},
        {"x", "y"},
        1,  // One solution with multiplicity 2
        1e-8
    });
    
    tests.push_back({
        "Triple root univariate",
        {"(x-2)^3"},
        {"x"},
        1,  // One solution with multiplicity 3
        1e-8
    });
    
    // 5. AxF4 benchmark (small ones only)
    std::string axf4_dir = "../axf4";
    if (std::filesystem::exists(axf4_dir + "/cyclic7.txt")) {
        try {
            auto cyclic7 = PolynomialSystemTester::load_axf4_file(
                axf4_dir + "/cyclic7.txt", "Cyclic-7"
            );
            cyclic7.expected_num_solutions = 924;
            cyclic7.residual_tolerance = 1e-6;
            tests.push_back(cyclic7);
        } catch (...) {
            std::cout << "Warning: Could not load Cyclic-7" << std::endl;
        }
    }
    
    if (std::filesystem::exists(axf4_dir + "/katsura7.txt")) {
        try {
            auto katsura7 = PolynomialSystemTester::load_axf4_file(
                axf4_dir + "/katsura7.txt", "Katsura-7"
            );
            katsura7.expected_num_solutions = 256;
            katsura7.residual_tolerance = 1e-6;
            tests.push_back(katsura7);
        } catch (...) {
            std::cout << "Warning: Could not load Katsura-7" << std::endl;
        }
    }
    
    // 6. Sparse high degree
    tests.push_back({
        "Sparse degree 20",
        {"x^20 + x^10 - x^5 + x - 1"},
        {"x"},
        -1,  // Up to 20 roots
        1e-6
    });
    
    // 7. Symmetric system
    tests.push_back({
        "Elementary symmetric",
        {"x + y + z - 3", "x*y + y*z + z*x - 3", "x*y*z - 1"},
        {"x", "y", "z"},
        6,  // 3! permutations
        1e-6
    });
    
    // Configure and run tests
    PolynomialSystemTester tester(true);  // verbose
    
    std::cout << "\n===================================================" << std::endl;
    std::cout << "RUR C++ Implementation - Comprehensive Test Suite" << std::endl;
    std::cout << "===================================================" << std::endl;
    
    tester.run_all_tests(tests);
    
    return 0;
}