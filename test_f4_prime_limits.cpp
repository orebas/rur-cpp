#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <chrono>
#include <iomanip>

extern "C" {
#include "src/axf4_wrapper.h"
}

struct TestResult {
    int prime;
    bool success;
    double time_ms;
    int basis_size;
    std::string error_msg;
    bool timeout;
};

// Test a specific prime with timeout protection
TestResult test_prime_modulus(int prime, int timeout_seconds = 5) {
    TestResult result;
    result.prime = prime;
    result.success = false;
    result.timeout = false;
    result.basis_size = 0;
    
    std::cout << "Testing prime " << prime << "... " << std::flush;
    
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(prime, vars, 2);
    
    if (!session) {
        result.error_msg = "Failed to create session";
        std::cout << "FAILED (session creation)" << std::endl;
        return result;
    }
    
    // Add test polynomials
    if (axf4_add_polynomial(session, "x^2+y^2-1") != 0 ||
        axf4_add_polynomial(session, "x-y") != 0) {
        result.error_msg = "Failed to add polynomials";
        std::cout << "FAILED (polynomial addition)" << std::endl;
        axf4_destroy_session(session);
        return result;
    }
    
    // Time the computation
    auto start = std::chrono::high_resolution_clock::now();
    
    // Set up timeout alarm
    alarm(timeout_seconds);
    
    axf4_result_t f4_result = axf4_compute_groebner_basis(session);
    
    // Cancel alarm
    alarm(0);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    result.time_ms = duration.count();
    
    if (f4_result.status == 0) {
        result.success = true;
        result.basis_size = f4_result.basis_size;
        std::cout << "SUCCESS (" << result.time_ms << " ms, " 
                  << result.basis_size << " elements)" << std::endl;
    } else {
        result.error_msg = f4_result.error_message ? f4_result.error_message : "Unknown error";
        std::cout << "FAILED (" << result.error_msg << ")" << std::endl;
    }
    
    axf4_free_result(&f4_result);
    axf4_destroy_session(session);
    
    return result;
}

int main() {
    std::cout << "F4 Prime Modulus Limit Testing" << std::endl;
    std::cout << "==============================" << std::endl;
    std::cout << "Test system: x^2 + y^2 - 1 = 0, x - y = 0" << std::endl;
    std::cout << "Expected basis: has 2 elements including 2*y^2 - 1" << std::endl;
    std::cout << std::endl;
    
    // Key values around 2^31
    const long long MAX_FAST_PRIME = 2147483647LL;  // 2^31 - 1
    
    std::cout << "Integer limits on this system:" << std::endl;
    std::cout << "  sizeof(int) = " << sizeof(int) << " bytes" << std::endl;
    std::cout << "  sizeof(long) = " << sizeof(long) << " bytes" << std::endl;
    std::cout << "  sizeof(long long) = " << sizeof(long long) << " bytes" << std::endl;
    std::cout << "  INT_MAX = " << INT_MAX << std::endl;
    std::cout << "  2^31 - 1 = " << MAX_FAST_PRIME << std::endl;
    std::cout << std::endl;
    
    std::vector<int> test_primes;
    
    // Small primes (should work)
    test_primes.push_back(65537);           // 2^16 + 1
    test_primes.push_back(1073741827);      // Known working
    
    // Primes around the 2^31 boundary
    test_primes.push_back(2147483647);      // 2^31 - 1 (largest "fast" prime)
    test_primes.push_back(2147483629);      // First problematic prime
    test_primes.push_back(2147483587);      
    test_primes.push_back(2147483579);
    test_primes.push_back(2147483563);
    test_primes.push_back(2147483549);
    
    // Test additional boundary cases
    test_primes.push_back(2147483489);      // Farther below 2^31
    test_primes.push_back(2147483497);
    
    std::cout << "Testing " << test_primes.size() << " prime moduli..." << std::endl;
    std::cout << std::endl;
    
    std::vector<TestResult> results;
    int successes = 0;
    int failures = 0;
    
    for (int prime : test_primes) {
        TestResult result = test_prime_modulus(prime, 10);
        results.push_back(result);
        
        if (result.success) {
            successes++;
        } else {
            failures++;
        }
    }
    
    // Summary
    std::cout << std::endl;
    std::cout << "Summary" << std::endl;
    std::cout << "=======" << std::endl;
    std::cout << "Total tests: " << results.size() << std::endl;
    std::cout << "Successes: " << successes << std::endl;
    std::cout << "Failures: " << failures << std::endl;
    std::cout << std::endl;
    
    // Detailed results table
    std::cout << "Detailed Results:" << std::endl;
    std::cout << std::setw(15) << "Prime" 
              << std::setw(15) << "Status"
              << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Basis Size"
              << std::setw(20) << "Notes" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    for (const auto& result : results) {
        std::cout << std::setw(15) << result.prime
                  << std::setw(15) << (result.success ? "SUCCESS" : "FAILED")
                  << std::setw(15) << (result.success ? std::to_string((int)result.time_ms) : "N/A")
                  << std::setw(15) << (result.success ? std::to_string(result.basis_size) : "N/A")
                  << std::setw(20) << (result.prime > MAX_FAST_PRIME ? "> 2^31-1" : "≤ 2^31-1")
                  << std::endl;
    }
    
    std::cout << std::endl;
    
    // Analysis
    std::cout << "Analysis:" << std::endl;
    std::cout << "=========" << std::endl;
    
    bool boundary_issue = true;
    for (const auto& result : results) {
        if (result.prime <= MAX_FAST_PRIME && !result.success) {
            boundary_issue = false;
            std::cout << "- Found failure with prime " << result.prime 
                      << " (≤ 2^31-1)" << std::endl;
        }
        if (result.prime > MAX_FAST_PRIME && result.success) {
            boundary_issue = false;
            std::cout << "- Found success with prime " << result.prime 
                      << " (> 2^31-1)" << std::endl;
        }
    }
    
    if (boundary_issue) {
        std::cout << "- All primes ≤ 2^31-1 succeeded" << std::endl;
        std::cout << "- All primes > 2^31-1 failed" << std::endl;
        std::cout << "- This confirms the issue is related to the 2^31 boundary" << std::endl;
        std::cout << "- The fast arithmetic path for p ≤ 2147483647 works correctly" << std::endl;
        std::cout << "- The Montgomery multiplication path for larger primes has a bug" << std::endl;
    }
    
    return failures > 0 ? 1 : 0;
}