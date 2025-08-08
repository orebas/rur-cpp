#include "../src/julia_rur/polynomial_evaluator.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace julia_rur;

void test_basic_evaluation() {
    std::cout << "Testing basic polynomial evaluation..." << std::endl;
    
    // Test 1: Simple univariate polynomial x^2 - 4
    {
        std::string poly = "x^2 - 4";
        std::vector<std::string> vars = {"x"};
        std::vector<std::complex<double>> values = {{2.0, 0.0}};
        
        auto result = PolynomialEvaluator::evaluate(poly, vars, values);
        std::cout << "  x^2 - 4 at x=2: " << result << " (expected: 0)" << std::endl;
        assert(std::abs(result) < 1e-10);
    }
    
    // Test 2: Multivariate polynomial x^2 + y^2 - 1
    {
        std::string poly = "x^2 + y^2 - 1";
        std::vector<std::string> vars = {"x", "y"};
        std::vector<std::complex<double>> values = {{0.6, 0.0}, {0.8, 0.0}};
        
        auto result = PolynomialEvaluator::evaluate(poly, vars, values);
        std::cout << "  x^2 + y^2 - 1 at (0.6, 0.8): " << result << " (expected: 0)" << std::endl;
        assert(std::abs(result) < 1e-10);
    }
    
    // Test 3: Complex values
    {
        std::string poly = "x^2 + 1";
        std::vector<std::string> vars = {"x"};
        std::vector<std::complex<double>> values = {{0.0, 1.0}};  // x = i
        
        auto result = PolynomialEvaluator::evaluate(poly, vars, values);
        std::cout << "  x^2 + 1 at x=i: " << result << " (expected: 0)" << std::endl;
        assert(std::abs(result) < 1e-10);
    }
    
    // Test 4: Polynomial with parentheses
    {
        std::string poly = "(x-1)^2";
        std::vector<std::string> vars = {"x"};
        std::vector<std::complex<double>> values = {{1.0, 0.0}};
        
        auto result = PolynomialEvaluator::evaluate(poly, vars, values);
        std::cout << "  (x-1)^2 at x=1: " << result << " (expected: 0)" << std::endl;
        assert(std::abs(result) < 1e-10);
    }
    
    // Test 5: Coefficients and multiplication
    {
        std::string poly = "2*x*y + 3*x - 5*y + 7";
        std::vector<std::string> vars = {"x", "y"};
        std::vector<std::complex<double>> values = {{2.0, 0.0}, {3.0, 0.0}};
        
        auto result = PolynomialEvaluator::evaluate(poly, vars, values);
        double expected = 2*2*3 + 3*2 - 5*3 + 7;  // = 12 + 6 - 15 + 7 = 10
        std::cout << "  2*x*y + 3*x - 5*y + 7 at (2,3): " << result << " (expected: " << expected << ")" << std::endl;
        assert(std::abs(result - std::complex<double>(expected, 0.0)) < 1e-10);
    }
    
    std::cout << "Basic evaluation tests PASSED!" << std::endl << std::endl;
}

void test_residual_computation() {
    std::cout << "Testing residual computation..." << std::endl;
    
    // Test system: x^2 - 4 = 0, solution x = 2
    {
        std::vector<std::string> polys = {"x^2 - 4"};
        std::vector<std::string> vars = {"x"};
        std::vector<std::complex<double>> solution = {{2.0, 0.0}};
        
        double residual = PolynomialEvaluator::compute_residual(polys, vars, solution);
        std::cout << "  System {x^2 - 4} at x=2, residual: " << residual << std::endl;
        assert(residual < 1e-10);
    }
    
    // Test system: x + y - 3 = 0, 2*x - y = 0, solution (1, 2)
    {
        std::vector<std::string> polys = {"x + y - 3", "2*x - y"};
        std::vector<std::string> vars = {"x", "y"};
        std::vector<std::complex<double>> solution = {{1.0, 0.0}, {2.0, 0.0}};
        
        double residual = PolynomialEvaluator::compute_residual(polys, vars, solution);
        std::cout << "  Linear system at (1,2), residual: " << residual << std::endl;
        assert(residual < 1e-10);
    }
    
    // Test with approximate solution
    {
        std::vector<std::string> polys = {"x^2 + y^2 - 1"};
        std::vector<std::string> vars = {"x", "y"};
        std::vector<std::complex<double>> solution = {{0.7071, 0.0}, {0.7071, 0.0}};  // Approximate sqrt(2)/2
        
        double residual = PolynomialEvaluator::compute_residual(polys, vars, solution);
        std::cout << "  Circle equation at approximate (√2/2, √2/2), residual: " << residual << std::endl;
        assert(residual < 1e-4);  // Should be small but not zero
    }
    
    std::cout << "Residual computation tests PASSED!" << std::endl << std::endl;
}

void test_edge_cases() {
    std::cout << "Testing edge cases..." << std::endl;
    
    // Test with indexed variables
    {
        std::string poly = "x1^2 + x2^2 - 5";
        std::vector<std::string> vars = {"x1", "x2"};
        std::vector<std::complex<double>> values = {{1.0, 0.0}, {2.0, 0.0}};
        
        auto result = PolynomialEvaluator::evaluate(poly, vars, values);
        std::cout << "  x1^2 + x2^2 - 5 at (1,2): " << result << " (expected: 0)" << std::endl;
        assert(std::abs(result) < 1e-10);
    }
    
    // Test with negative coefficients
    {
        std::string poly = "-3*x + 6";
        std::vector<std::string> vars = {"x"};
        std::vector<std::complex<double>> values = {{2.0, 0.0}};
        
        auto result = PolynomialEvaluator::evaluate(poly, vars, values);
        std::cout << "  -3*x + 6 at x=2: " << result << " (expected: 0)" << std::endl;
        assert(std::abs(result) < 1e-10);
    }
    
    // Test with decimal coefficients
    {
        std::string poly = "0.5*x - 1.5";
        std::vector<std::string> vars = {"x"};
        std::vector<std::complex<double>> values = {{3.0, 0.0}};
        
        auto result = PolynomialEvaluator::evaluate(poly, vars, values);
        std::cout << "  0.5*x - 1.5 at x=3: " << result << " (expected: 0)" << std::endl;
        assert(std::abs(result) < 1e-10);
    }
    
    std::cout << "Edge case tests PASSED!" << std::endl << std::endl;
}

int main() {
    std::cout << "\n=====================================\n";
    std::cout << "Polynomial Evaluator Test Suite\n";
    std::cout << "=====================================\n\n";
    
    try {
        test_basic_evaluation();
        test_residual_computation();
        test_edge_cases();
        
        std::cout << "All tests PASSED!\n";
        std::cout << "=====================================\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test FAILED with exception: " << e.what() << std::endl;
        return 1;
    }
}