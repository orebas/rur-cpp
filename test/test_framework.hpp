#pragma once

#include "../src/julia_rur/polynomial_solver_enhanced.hpp"
#include "../src/julia_rur/polynomial_evaluator.hpp"
#include <complex>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <set>
#include <map>
#include <algorithm>

namespace test_framework {

struct TestResult {
    std::string test_name;
    bool passed;
    std::string error_message;
    double elapsed_seconds;
    int num_solutions;
    double max_residual;
};

struct TestCase {
    std::string name;
    std::vector<std::string> polynomials;
    std::vector<std::string> variables;
    int expected_num_solutions;  // -1 if unknown
    double residual_tolerance;
};

class PolynomialSystemTester {
public:
    PolynomialSystemTester(bool verbose = true) : verbose_(verbose) {}
    
    // Run a single test case
    TestResult run_test(const TestCase& test_case, double max_seconds = 20.0);
    
    // Run all test cases and print summary
    void run_all_tests(const std::vector<TestCase>& test_cases);
    
    // Load polynomial system from axf4 format file
    static TestCase load_axf4_file(const std::string& filename, const std::string& test_name);
    
    // Generate test cases
    static std::vector<TestCase> generate_standard_tests();
    
    // Compute residual for a solution
    static double compute_residual(
        const std::vector<std::string>& polynomials,
        const std::vector<std::string>& variables,
        const std::vector<std::complex<double>>& solution
    );

private:
    bool verbose_;
    std::vector<TestResult> results_;
    
    void print_test_result(const TestResult& result);
    void print_summary();
};

// Implementation

TestResult PolynomialSystemTester::run_test(const TestCase& test_case, double max_seconds) {
    TestResult result;
    result.test_name = test_case.name;
    result.passed = false;
    result.num_solutions = 0;
    result.max_residual = 0.0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Configure solver
        julia_rur::EnhancedSolverConfig config;
        config.verbose = false;  // Reduce output during tests
        // Note: timeout not implemented yet in solver
        
        // Solve the system
        auto solution = julia_rur::solve_polynomial_system_enhanced(
            test_case.polynomials, 
            test_case.variables, 
            config
        );
        
        auto end_time = std::chrono::high_resolution_clock::now();
        result.elapsed_seconds = std::chrono::duration<double>(end_time - start_time).count();
        
        if (!solution.success) {
            result.error_message = "Solver failed: " + solution.error_message;
            return result;
        }
        
        result.num_solutions = solution.solutions.size();
        
        // Check residuals
        result.max_residual = 0.0;
        for (const auto& sol : solution.solutions) {
            double residual = compute_residual(test_case.polynomials, test_case.variables, sol);
            result.max_residual = std::max(result.max_residual, residual);
        }
        
        // Check if test passed
        bool residual_ok = result.max_residual < test_case.residual_tolerance;
        bool count_ok = (test_case.expected_num_solutions < 0) || 
                       (result.num_solutions == test_case.expected_num_solutions);
        
        if (!residual_ok) {
            result.error_message = "Max residual " + std::to_string(result.max_residual) + 
                                 " exceeds tolerance " + std::to_string(test_case.residual_tolerance);
        } else if (!count_ok) {
            result.error_message = "Expected " + std::to_string(test_case.expected_num_solutions) + 
                                 " solutions, got " + std::to_string(result.num_solutions);
        } else {
            result.passed = true;
        }
        
    } catch (const std::exception& e) {
        auto end_time = std::chrono::high_resolution_clock::now();
        result.elapsed_seconds = std::chrono::duration<double>(end_time - start_time).count();
        result.error_message = "Exception: " + std::string(e.what());
    }
    
    return result;
}

void PolynomialSystemTester::run_all_tests(const std::vector<TestCase>& test_cases) {
    results_.clear();
    
    std::cout << "\n=== Running Polynomial System Tests ===\n" << std::endl;
    
    for (const auto& test_case : test_cases) {
        if (verbose_) {
            std::cout << "Running: " << test_case.name << "... " << std::flush;
        }
        
        auto result = run_test(test_case);
        results_.push_back(result);
        
        if (verbose_) {
            print_test_result(result);
        }
    }
    
    print_summary();
}

void PolynomialSystemTester::print_test_result(const TestResult& result) {
    if (result.passed) {
        std::cout << "PASSED";
    } else {
        std::cout << "FAILED";
    }
    
    std::cout << " (" << std::fixed << std::setprecision(2) << result.elapsed_seconds << "s, "
              << result.num_solutions << " solutions, "
              << "max residual: " << std::scientific << std::setprecision(2) << result.max_residual << ")";
    
    if (!result.passed && !result.error_message.empty()) {
        std::cout << "\n  Error: " << result.error_message;
    }
    
    std::cout << std::endl;
}

void PolynomialSystemTester::print_summary() {
    int passed = 0;
    int failed = 0;
    double total_time = 0.0;
    
    for (const auto& result : results_) {
        if (result.passed) passed++;
        else failed++;
        total_time += result.elapsed_seconds;
    }
    
    std::cout << "\n=== Test Summary ===" << std::endl;
    std::cout << "Total tests: " << (passed + failed) << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total time: " << std::fixed << std::setprecision(2) << total_time << " seconds" << std::endl;
    
    if (failed > 0) {
        std::cout << "\nFailed tests:" << std::endl;
        for (const auto& result : results_) {
            if (!result.passed) {
                std::cout << "  - " << result.test_name << ": " << result.error_message << std::endl;
            }
        }
    }
}

TestCase PolynomialSystemTester::load_axf4_file(const std::string& filename, const std::string& test_name) {
    TestCase test_case;
    test_case.name = test_name;
    test_case.expected_num_solutions = -1;  // Unknown
    test_case.residual_tolerance = 1e-8;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
            continue;
        }
        test_case.polynomials.push_back(line);
    }
    
    // Extract variables (assume x0, x1, ..., xn format)
    std::set<std::string> var_set;
    for (const auto& poly : test_case.polynomials) {
        size_t pos = 0;
        while ((pos = poly.find('x', pos)) != std::string::npos) {
            if (pos + 1 < poly.size() && std::isdigit(poly[pos + 1])) {
                size_t end = pos + 1;
                while (end < poly.size() && std::isdigit(poly[end])) {
                    end++;
                }
                var_set.insert(poly.substr(pos, end - pos));
            }
            pos++;
        }
    }
    
    // Sort variables by index
    std::vector<std::pair<int, std::string>> indexed_vars;
    for (const auto& var : var_set) {
        int idx = std::stoi(var.substr(1));  // Skip 'x'
        indexed_vars.push_back(std::make_pair(idx, var));
    }
    std::sort(indexed_vars.begin(), indexed_vars.end());
    
    for (const auto& pair : indexed_vars) {
        test_case.variables.push_back(pair.second);
    }
    
    return test_case;
}


double PolynomialSystemTester::compute_residual(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const std::vector<std::complex<double>>& solution
) {
    // Use the PolynomialEvaluator to compute residuals
    return julia_rur::PolynomialEvaluator::compute_residual(polynomials, variables, solution);
}

std::vector<TestCase> PolynomialSystemTester::generate_standard_tests() {
    std::vector<TestCase> tests;
    
    // 1. Trivial univariate
    tests.push_back({
        "Univariate x^2 - 5",
        {"x^2 - 5"},
        {"x"},
        2,  // ±sqrt(5)
        1e-8
    });
    
    // 2. Linear 2x2 system
    tests.push_back({
        "Linear 2x2",
        {"x + y - 3", "2*x - y - 1"},
        {"x", "y"},
        1,  // One solution
        1e-8
    });
    
    // 3. Circle and line
    tests.push_back({
        "Circle and line",
        {"x^2 + y^2 - 1", "y - x"},
        {"x", "y"},
        2,  // Two intersection points
        1e-8
    });
    
    // 4. Univariate degree 10
    tests.push_back({
        "Univariate degree 10",
        {"x^10 - 1"},
        {"x"},
        10,  // 10 roots (complex roots of unity)
        1e-8
    });
    
    // 5. System with multiplicity
    tests.push_back({
        "Double root",
        {"(x-1)^2", "y"},
        {"x", "y"},
        1,  // One solution with multiplicity 2
        1e-8
    });
    
    // 6. Linear 5x5 system
    tests.push_back({
        "Linear 5x5",
        {
            "x1 + 2*x2 + 3*x3 + 4*x4 + 5*x5 - 15",
            "2*x1 + 3*x2 + 4*x3 + 5*x4 + x5 - 14",
            "3*x1 + 4*x2 + 5*x3 + x4 + 2*x5 - 13",
            "4*x1 + 5*x2 + x3 + 2*x4 + 3*x5 - 12",
            "5*x1 + x2 + 2*x3 + 3*x4 + 4*x5 - 11"
        },
        {"x1", "x2", "x3", "x4", "x5"},
        1,  // One solution
        1e-8
    });
    
    // 7. Squared linear 5x5 system
    tests.push_back({
        "Squared linear 5x5",
        {
            "x1^2 + 2*x2^2 + 3*x3^2 + 4*x4^2 + 5*x5^2 - 15",
            "2*x1^2 + 3*x2^2 + 4*x3^2 + 5*x4^2 + x5^2 - 14",
            "3*x1^2 + 4*x2^2 + 5*x3^2 + x4^2 + 2*x5^2 - 13",
            "4*x1^2 + 5*x2^2 + x3^2 + 2*x4^2 + 3*x5^2 - 12",
            "5*x1^2 + x2^2 + 2*x3^2 + 3*x4^2 + 4*x5^2 - 11"
        },
        {"x1", "x2", "x3", "x4", "x5"},
        32,  // 2^5 solutions (each variable can be ±)
        1e-8
    });
    
    return tests;
}

} // namespace test_framework