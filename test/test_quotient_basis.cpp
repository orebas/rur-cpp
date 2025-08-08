#include "../src/julia_rur/quotient_basis.hpp"
#include <iostream>
#include <cassert>
#include <algorithm>

using namespace julia_rur;

void print_pp(const PP& pp, const std::vector<std::string>& vars) {
    bool first = true;
    bool all_zero = true;
    
    for (size_t i = 0; i < pp.size(); i++) {
        if (pp[i] > 0) {
            all_zero = false;
            if (!first) std::cout << "*";
            std::cout << vars[i];
            if (pp[i] > 1) std::cout << "^" << pp[i];
            first = false;
        }
    }
    
    if (all_zero) std::cout << "1";
}

void print_quotient_basis(const std::vector<PP>& basis, const std::vector<std::string>& vars) {
    std::cout << "Quotient basis has " << basis.size() << " elements: {";
    for (size_t i = 0; i < basis.size(); i++) {
        if (i > 0) std::cout << ", ";
        print_pp(basis[i], vars);
    }
    std::cout << "}" << std::endl;
}

void test_simple_linear() {
    std::cout << "\n=== Test 1: Simple linear system {x-1} ===" << std::endl;
    
    // Leading term of x-1 is x
    std::vector<PP> leading_terms = {PP({1})};  // x^1
    
    auto quotient_basis = compute_quotient_basis(leading_terms);
    print_quotient_basis(quotient_basis, {"x"});
    
    // Should be {1} since x divides all higher powers
    assert(quotient_basis.size() == 1);
    assert(quotient_basis[0] == PP({0}));  // The constant 1
    
    std::cout << "✓ Test passed!" << std::endl;
}

void test_two_variable_system() {
    std::cout << "\n=== Test 2: Two variable system {x^2+y^2-1, x-y} ===" << std::endl;
    
    // In degrevlex order, the Gröbner basis would typically have leading terms
    // like {x^2, y^2} or similar (depends on the specific computation)
    // For this test, let's assume a typical result: {y^2, x-y} -> LT = {y^2, x}
    std::vector<PP> leading_terms = {
        PP({0, 2}),  // y^2
        PP({1, 0})   // x
    };
    
    auto quotient_basis = compute_quotient_basis(leading_terms);
    print_quotient_basis(quotient_basis, {"x", "y"});
    
    // Quotient basis should be {1, y} since:
    // - 1 is not divisible by y^2 or x
    // - y is not divisible by y^2 or x
    // - x is divisible by x (excluded)
    // - y^2 is divisible by y^2 (excluded)
    // - xy is divisible by x (excluded)
    assert(quotient_basis.size() == 2);
    assert(quotient_basis[0] == PP({0, 0}));  // 1
    assert(quotient_basis[1] == PP({0, 1}));  // y
    
    std::cout << "✓ Test passed!" << std::endl;
}

void test_cyclic3_leading_terms() {
    std::cout << "\n=== Test 3: Cyclic-3 system leading terms ===" << std::endl;
    
    // For a zero-dimensional ideal, we need univariate terms for each variable
    // Real cyclic-3 example might have these leading terms:
    std::vector<PP> leading_terms;
    leading_terms.push_back(PP({0, 0, 3}));  // z^3 (univariate in z)
    leading_terms.push_back(PP({0, 2, 0}));  // y^2 (univariate in y) 
    leading_terms.push_back(PP({1, 0, 1}));  // x*z
    leading_terms.push_back(PP({2, 0, 0}));  // x^2 (univariate in x)
    
    
    auto quotient_basis = compute_quotient_basis(leading_terms);
    print_quotient_basis(quotient_basis, {"x", "y", "z"});
    
    // Expected quotient basis elements:
    // 1, z, z^2, y (but not x, yz, y*z^2, z^3, etc.)
    std::cout << "Quotient basis size: " << quotient_basis.size() << std::endl;
    
    // Verify expected properties
    // With LT = {z^3, y^2, x*z, x^2}, quotient basis should be:
    // All monomials not divisible by any leading term
    assert(quotient_basis.size() == 8);  // Expected size
    // Helper lambda to check if monomial is in basis
    auto is_in_basis = [&](const PP& monomial) {
        return std::find(quotient_basis.begin(), quotient_basis.end(), monomial) != quotient_basis.end();
    };
    
    assert(is_in_basis(PP({0, 0, 0})));  // 1 is in basis
    assert(is_in_basis(PP({0, 0, 1})));  // z is in basis
    assert(is_in_basis(PP({0, 0, 2})));  // z^2 is in basis
    assert(is_in_basis(PP({0, 1, 0})));  // y is in basis
    assert(is_in_basis(PP({1, 0, 0})));  // x IS in basis (not divisible by x^2 or x*z)
    assert(is_in_basis(PP({0, 1, 1})));  // y*z IS in basis
    assert(!is_in_basis(PP({2, 0, 0}))); // x^2 is NOT in basis (divisible by x^2)
    assert(!is_in_basis(PP({0, 0, 3}))); // z^3 is NOT in basis (divisible by z^3)
    
    std::cout << "✓ Test passed!" << std::endl;
}

void test_non_zero_dimensional() {
    std::cout << "\n=== Test 4: Non-zero-dimensional ideal ===" << std::endl;
    
    // System like {x*y} - missing univariate terms for both x and y
    std::vector<PP> leading_terms = {
        PP({1, 1})  // xy
    };
    
    try {
        auto quotient_basis = compute_quotient_basis(leading_terms);
        assert(false && "Should have thrown exception for non-zero-dimensional ideal");
    } catch (const std::domain_error& e) {
        std::cout << "✓ Correctly detected non-zero-dimensional ideal: " << e.what() << std::endl;
    }
}

void test_no_solutions() {
    std::cout << "\n=== Test 5: System with no solutions (GB = {1}) ===" << std::endl;
    
    // Gröbner basis = {1}
    std::vector<PP> leading_terms = {
        PP({0, 0})  // The constant 1
    };
    
    try {
        auto quotient_basis = compute_quotient_basis(leading_terms);
        assert(false && "Should have thrown exception for system with no solutions");
    } catch (const std::runtime_error& e) {
        std::cout << "✓ Correctly detected no solutions: " << e.what() << std::endl;
    }
}

void test_larger_quotient_basis() {
    std::cout << "\n=== Test 6: System with larger quotient basis ===" << std::endl;
    
    // A system where quotient basis has dimension > 2
    // Example: leading terms that create a 4-dimensional quotient space
    std::vector<PP> leading_terms = {
        PP({2, 0}),   // x^2
        PP({0, 2})    // y^2
    };
    
    auto quotient_basis = compute_quotient_basis(leading_terms);
    print_quotient_basis(quotient_basis, {"x", "y"});
    
    // Quotient basis should be {1, x, y, xy}
    assert(quotient_basis.size() == 4);
    assert(quotient_basis[0] == PP({0, 0}));  // 1
    assert(quotient_basis[1] == PP({0, 1}));  // y (comes before x in degrevlex)
    assert(quotient_basis[2] == PP({1, 0}));  // x
    assert(quotient_basis[3] == PP({1, 1}));  // xy
    
    std::cout << "✓ Test passed!" << std::endl;
}

int main() {
    std::cout << "=== Testing Quotient Basis Extraction ===" << std::endl;
    
    try {
        test_simple_linear();
        test_two_variable_system();
        test_cyclic3_leading_terms();
        test_non_zero_dimensional();
        test_no_solutions();
        test_larger_quotient_basis();
        
        std::cout << "\n✅ All quotient basis tests passed!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n❌ Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
}