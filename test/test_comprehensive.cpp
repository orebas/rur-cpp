#include "test_framework.hpp"
#include <filesystem>
#include <algorithm>

using namespace test_framework;

// Test cases from Julia package
std::vector<TestCase> get_julia_tests() {
    std::vector<TestCase> tests;
    
    // Example from Julia tests
    tests.push_back({
        "Julia example 1",
        {"x0^2 + 2*x1^2 + 2*x2^2 - x0", "2*x0*x1 + 2*x1*x2 - x1", "x0 + 2*x1 + 2*x2 - 1"},
        {"x0", "x1", "x2"},
        -1,  // Unknown number of solutions
        1e-8
    });
    
    // Trivial system
    tests.push_back({
        "Julia trivial",
        {"x0", "x1", "x2"},
        {"x0", "x1", "x2"},
        1,  // Origin
        1e-8
    });
    
    // Large coefficient system
    tests.push_back({
        "Julia large coefficients",
        {"x0 - 1267650600228229401496703205376", "x1 - 1", "x2 + 2535301200456458802993406410752"},
        {"x0", "x1", "x2"},
        1,
        1e-8
    });
    
    return tests;
}

// Load all axf4 test files
std::vector<TestCase> load_axf4_tests(const std::string& axf4_dir) {
    std::vector<TestCase> tests;
    
    // Systems to load with expected characteristics
    struct AxF4Test {
        std::string filename;
        std::string name;
        int expected_solutions;  // -1 if unknown
        double tolerance;
    };
    
    std::vector<AxF4Test> axf4_tests = {
        {"cyclic7.txt", "Cyclic-7", 924, 1e-6},      // Known to have 924 solutions
        {"cyclic8.txt", "Cyclic-8", -1, 1e-6},       // Too large for quick test
        {"katsura7.txt", "Katsura-7", 256, 1e-6},    // 2^7 solutions
        {"katsura8.txt", "Katsura-8", -1, 1e-6},     // Skip - too large
        {"noon9.txt", "Noon-9", -1, 1e-6},           // Unknown count
        {"gametwo.txt", "Game Two", -1, 1e-6},       // Unknown count
        {"eco12.txt", "Eco-12", -1, 1e-6},           // Skip if too large
        {"jason210.txt", "Jason 2-10", -1, 1e-6}     // Unknown count
    };
    
    for (const auto& test_info : axf4_tests) {
        try {
            std::string filepath = axf4_dir + "/" + test_info.filename;
            if (std::filesystem::exists(filepath)) {
                auto test_case = PolynomialSystemTester::load_axf4_file(filepath, test_info.name);
                test_case.expected_num_solutions = test_info.expected_solutions;
                test_case.residual_tolerance = test_info.tolerance;
                
                // Skip tests that are known to be too large for 20 second timeout
                if (test_info.name == "Cyclic-8" || test_info.name == "Katsura-8" || 
                    test_info.name == "Eco-12") {
                    std::cout << "Skipping " << test_info.name << " (too large for quick test)" << std::endl;
                    continue;
                }
                
                tests.push_back(test_case);
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not load " << test_info.filename << ": " << e.what() << std::endl;
        }
    }
    
    return tests;
}

// Additional stress tests
std::vector<TestCase> get_stress_tests() {
    std::vector<TestCase> tests;
    
    // Random dense quadratic system
    tests.push_back({
        "Dense quadratic 3x3",
        {
            "x^2 + 2*y^2 + 3*z^2 + x*y + x*z + y*z - 1",
            "2*x^2 + y^2 + z^2 + 2*x*y + 3*x*z + y*z - 2",
            "x^2 + y^2 + 2*z^2 + 3*x*y + x*z + 2*y*z - 3"
        },
        {"x", "y", "z"},
        -1,
        1e-6
    });
    
    // Sparse high degree univariate
    tests.push_back({
        "Sparse degree 20",
        {"x^20 + x^10 + x^5 + x - 1"},
        {"x"},
        -1,  // Up to 20 roots
        1e-6
    });
    
    // System with symmetry
    tests.push_back({
        "Symmetric system",
        {"x^2 + y^2 + z^2 - 3", "x*y + y*z + z*x - 1", "x*y*z - 1"},
        {"x", "y", "z"},
        -1,
        1e-6
    });
    
    return tests;
}

int main(int argc, char* argv[]) {
    bool verbose = true;
    std::string axf4_dir = "../axf4";  // Default location
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--quiet" || arg == "-q") {
            verbose = false;
        } else if (arg == "--axf4-dir" && i + 1 < argc) {
            axf4_dir = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  --quiet, -q          Reduce output verbosity" << std::endl;
            std::cout << "  --axf4-dir <path>    Path to axf4 test files (default: ../axf4)" << std::endl;
            std::cout << "  --help, -h           Show this help message" << std::endl;
            return 0;
        }
    }
    
    PolynomialSystemTester tester(verbose);
    std::vector<TestCase> all_tests;
    
    // 1. Standard generated tests
    std::cout << "Loading standard tests..." << std::endl;
    auto standard_tests = PolynomialSystemTester::generate_standard_tests();
    all_tests.insert(all_tests.end(), standard_tests.begin(), standard_tests.end());
    
    // 2. Julia package tests
    std::cout << "Loading Julia package tests..." << std::endl;
    auto julia_tests = get_julia_tests();
    all_tests.insert(all_tests.end(), julia_tests.begin(), julia_tests.end());
    
    // 3. AxF4 benchmark tests
    std::cout << "Loading AxF4 tests from " << axf4_dir << "..." << std::endl;
    auto axf4_tests = load_axf4_tests(axf4_dir);
    all_tests.insert(all_tests.end(), axf4_tests.begin(), axf4_tests.end());
    
    // 4. Additional stress tests
    std::cout << "Loading stress tests..." << std::endl;
    auto stress_tests = get_stress_tests();
    all_tests.insert(all_tests.end(), stress_tests.begin(), stress_tests.end());
    
    // Run all tests
    tester.run_all_tests(all_tests);
    
    return 0;
}