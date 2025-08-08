#include <iostream>
#include <string>
#include <vector>
#include "../src/julia_rur/f4_polynomial_formatter.hpp"

using namespace julia_rur;

struct FormatterTest {
    std::string input;
    std::string expected;
    std::string description;
};

bool run_test(const FormatterTest& test) {
    std::string result = format_polynomial_for_f4(test.input);
    bool passed = (result == test.expected);
    
    std::cout << "\nTest: " << test.description << std::endl;
    std::cout << "  Input:    \"" << test.input << "\"" << std::endl;
    std::cout << "  Expected: \"" << test.expected << "\"" << std::endl;
    std::cout << "  Got:      \"" << result << "\"" << std::endl;
    std::cout << "  Status:   " << (passed ? "PASSED" : "FAILED") << std::endl;
    
    return passed;
}

int main() {
    std::vector<FormatterTest> tests = {
        // Basic decimal coefficients
        {"0.5*x", "1*x", "Simple decimal 0.5 -> clear to integer"},
        {"y - 0.5", "2*y-1", "Circle-line test case: y - 0.5"},
        {"y- 0.5", "2*y-1", "Circle-line test case without space"},
        {"x^2 + y^2 - 1", "1*x^2+1*y^2-1", "Circle equation"},
        
        // Various decimal forms
        {"0.25*x + 0.5*y - 0.75", "1*x+2*y-3", "Multiple decimals with LCM 4"},
        {"0.333333333*x", "333333333*x", "Repeating decimal approximation"},
        {"0.1*x + 0.2*y", "1*x+2*y", "Decimals 0.1 and 0.2"},
        {"1.5*x^2 - 2.5*y", "3*x^2-5*y", "Mixed integers and decimals"},
        
        // Fractions that are already cleared
        {"2*x + 3*y - 5", "2*x+3*y-5", "Already integer coefficients"},
        {"x - y", "1*x-1*y", "Implicit coefficients of 1"},
        {"-x + y", "-1*x+1*y", "Implicit coefficient of -1"},
        
        // Edge cases
        {"0", "0", "Zero polynomial"},
        {"1", "1", "Constant polynomial"},
        {"0.5", "1", "Decimal constant"},
        {"-0.5", "-1", "Negative decimal constant"},
        
        // Complex expressions with multiple variables
        {"0.5*x*y + 0.25*z", "2*x*y+1*z", "Product terms with decimals"},
        {"0.3*x^2*y - 0.6*x*y^2", "1*x^2*y-2*x*y^2", "Multiple variable products"},
        
        // Scientific notation (if supported)
        {"1e-1*x", "1*x", "Scientific notation 1e-1 = 0.1"},
        {"2.5e-1*x", "1*x", "Scientific notation 2.5e-1 = 0.25"},
        
        // Mixed positive and negative
        {"0.5*x - 0.25*y + 0.125*z", "4*x-2*y+1*z", "Mixed signs with LCM 8"},
        {"-0.5*x - 0.5*y", "-1*x-1*y", "All negative decimals"},
        
        // Stress test with many terms
        {"0.1*a + 0.2*b + 0.3*c + 0.4*d + 0.5*e", 
         "1*a+2*b+3*c+4*d+5*e", 
         "Many terms with decimals"},
         
        // Parenthetical expressions (if expanded correctly)
        {"(x - 0.5)^2", "1*x^2-2*x+1", "Expanded and cleared: (x - 0.5)^2"},
        
        // Very small decimals
        {"0.001*x", "1*x", "Very small decimal 0.001"},
        {"0.0001*x + 0.0002*y", "1*x+2*y", "Very small decimals"}
    };
    
    int passed = 0;
    int failed = 0;
    
    std::cout << "========================================" << std::endl;
    std::cout << "Polynomial Formatter Fraction Tests" << std::endl;
    std::cout << "========================================" << std::endl;
    
    for (const auto& test : tests) {
        if (run_test(test)) {
            passed++;
        } else {
            failed++;
        }
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Summary:" << std::endl;
    std::cout << "  Passed: " << passed << std::endl;
    std::cout << "  Failed: " << failed << std::endl;
    std::cout << "  Total:  " << (passed + failed) << std::endl;
    std::cout << "========================================" << std::endl;
    
    return failed > 0 ? 1 : 0;
}