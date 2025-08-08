#include "f4_solver.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

void test_simple_system() {
    std::cout << "Testing simple polynomial system...\n";
    
    // Test system: x^2 + y^2 - 1, x - y
    // Should have Gröbner basis with 2y^2 - 1, x - y (or similar)
    
    std::vector<std::string> vars = {"x", "y"};
    F4Solver<int> solver(100003, vars);  // Use prime > 65536
    
    // Create polynomials: x^2 + y^2 - 1
    MultivariatePolynomial<int> p1;
    p1.set_coefficient(Monomial(std::vector<int>{2, 0}), 1);   // x^2
    p1.set_coefficient(Monomial(std::vector<int>{0, 2}), 1);   // y^2  
    p1.set_coefficient(Monomial(std::vector<int>{0, 0}), -1);  // -1
    
    // x - y
    MultivariatePolynomial<int> p2;
    p2.set_coefficient(Monomial(std::vector<int>{1, 0}), 1);   // x
    p2.set_coefficient(Monomial(std::vector<int>{0, 1}), -1);  // -y
    
    solver.add_polynomial(p1);
    solver.add_polynomial(p2);
    
    std::cout << "Input polynomials:\n";
    std::cout << "p1 = " << p1.to_string(vars) << std::endl;
    std::cout << "p2 = " << p2.to_string(vars) << std::endl;
    
    solver.compute_groebner_basis();
    
    const auto& basis = solver.get_groebner_basis();
    
    std::cout << "Gröbner basis (" << basis.size() << " polynomials):\n";
    for (size_t i = 0; i < basis.size(); i++) {
        std::cout << "g" << i+1 << " = " << basis[i].to_string(vars) << std::endl;
    }
    
    std::cout << "Computation time: " << solver.get_computation_time() << " seconds\n";
    
    assert(basis.size() >= 2);
    std::cout << "✓ Simple system test passed!\n\n";
}

void test_cyclic3_system() {
    std::cout << "Testing cyclic-3 system...\n";
    
    // Cyclic-3: x + y + z, xy + yz + zx, xyz - 1
    std::vector<std::string> vars = {"x", "y", "z"};
    F4Solver<int> solver(100003, vars);  // Use prime > 65536
    
    // x + y + z
    MultivariatePolynomial<int> p1;
    p1.set_coefficient(Monomial(std::vector<int>{1, 0, 0}), 1);  // x
    p1.set_coefficient(Monomial(std::vector<int>{0, 1, 0}), 1);  // y
    p1.set_coefficient(Monomial(std::vector<int>{0, 0, 1}), 1);  // z
    
    // xy + yz + zx  
    MultivariatePolynomial<int> p2;
    p2.set_coefficient(Monomial(std::vector<int>{1, 1, 0}), 1);  // xy
    p2.set_coefficient(Monomial(std::vector<int>{0, 1, 1}), 1);  // yz
    p2.set_coefficient(Monomial(std::vector<int>{1, 0, 1}), 1);  // zx
    
    // xyz - 1
    MultivariatePolynomial<int> p3;
    p3.set_coefficient(Monomial(std::vector<int>{1, 1, 1}), 1);  // xyz
    p3.set_coefficient(Monomial(std::vector<int>{0, 0, 0}), -1); // -1
    
    solver.add_polynomials({p1, p2, p3});
    
    std::cout << "Input polynomials:\n";
    std::cout << "p1 = " << p1.to_string(vars) << std::endl;
    std::cout << "p2 = " << p2.to_string(vars) << std::endl;
    std::cout << "p3 = " << p3.to_string(vars) << std::endl;
    
    solver.compute_groebner_basis();
    
    const auto& basis = solver.get_groebner_basis();
    
    std::cout << "Gröbner basis (" << basis.size() << " polynomials):\n";
    for (size_t i = 0; i < basis.size(); i++) {
        std::cout << "g" << i+1 << " = " << basis[i].to_string(vars) << std::endl;
    }
    
    std::cout << "Computation time: " << solver.get_computation_time() << " seconds\n";
    
    // Cyclic-3 should have a finite number of solutions, so basis should be non-trivial
    assert(basis.size() >= 3);
    std::cout << "✓ Cyclic-3 system test passed!\n\n";
}

void test_error_handling() {
    std::cout << "Testing error handling...\n";
    
    // Test invalid prime
    try {
        F4Solver<int> solver(1000, {"x", "y"});  // Prime too small
        assert(false && "Should have thrown exception");
    } catch (const std::invalid_argument& e) {
        std::cout << "✓ Caught expected error for small prime: " << e.what() << std::endl;
    }
    
    // Test empty variables
    try {
        F4Solver<int> solver(32003, {});
        assert(false && "Should have thrown exception");  
    } catch (const std::invalid_argument& e) {
        std::cout << "✓ Caught expected error for empty variables: " << e.what() << std::endl;
    }
    
    // Test computing basis without polynomials
    F4Solver<int> solver(100003, {"x"});
    try {
        solver.compute_groebner_basis();
        assert(false && "Should have thrown exception");
    } catch (const std::runtime_error& e) {
        std::cout << "✓ Caught expected error for empty system: " << e.what() << std::endl;
    }
    
    std::cout << "✓ Error handling tests passed!\n\n";
}

void test_solver_properties() {
    std::cout << "Testing solver properties...\n";
    
    std::vector<std::string> vars = {"a", "b"};
    F4Solver<int> solver(100003, vars);
    
    assert(solver.prime() == 100003);
    assert(solver.variable_names() == vars);
    assert(!solver.is_basis_computed());
    assert(solver.basis_size() == 0);
    
    // Add a simple polynomial
    MultivariatePolynomial<int> p;
    p.set_coefficient(Monomial(std::vector<int>{1, 0}), 1);  // a
    p.set_coefficient(Monomial(std::vector<int>{0, 0}), -1); // -1
    
    solver.add_polynomial(p);
    solver.compute_groebner_basis();
    
    assert(solver.is_basis_computed());
    assert(solver.basis_size() > 0);
    assert(solver.get_computation_time() >= 0);
    
    std::cout << "✓ Solver properties tests passed!\n\n";
}

int main() {
    try {
        std::cout << "=== F4 Solver C++ Wrapper Tests ===\n\n";
        
        test_error_handling();
        std::cout << "About to test solver properties...\n";
        test_solver_properties();
        std::cout << "About to test simple system...\n";
        test_simple_system();
        std::cout << "About to test cyclic3 system...\n";
        test_cyclic3_system();
        
        std::cout << "=== All tests passed! ===\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
}