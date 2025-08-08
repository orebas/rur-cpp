#include "../src/flint_linear_algebra.hpp"
#include <iostream>
#include <cassert>
#include <vector>

using namespace flint_linalg;

void test_nmod_mat_basic() {
    std::cout << "Testing NModMat basic operations..." << std::endl;
    
    ulong prime = 7;
    NModMat mat(3, 3, prime);
    
    // Test initialization
    assert(mat.rows() == 3);
    assert(mat.cols() == 3);
    assert(mat.prime() == prime);
    
    // Test zero initialization
    for (slong i = 0; i < 3; i++) {
        for (slong j = 0; j < 3; j++) {
            assert(mat.get_entry(i, j) == 0);
        }
    }
    
    // Test element setting/getting
    mat.set_entry(0, 0, 3);
    mat.set_entry(1, 1, 5);
    mat.set_entry(2, 2, 2);
    
    assert(mat.get_entry(0, 0) == 3);
    assert(mat.get_entry(1, 1) == 5);
    assert(mat.get_entry(2, 2) == 2);
    
    // Test identity matrix
    NModMat id(3, 3, prime);
    id.one();
    for (slong i = 0; i < 3; i++) {
        for (slong j = 0; j < 3; j++) {
            if (i == j) {
                assert(id.get_entry(i, j) == 1);
            } else {
                assert(id.get_entry(i, j) == 0);
            }
        }
    }
    
    std::cout << "✓ NModMat basic operations passed" << std::endl;
}

void test_linear_solver_rref() {
    std::cout << "Testing LinearSolver RREF..." << std::endl;
    
    ulong prime = 7;
    LinearSolver solver(prime);
    
    // Test simple 2x2 matrix [[1, 1], [1, -1]] mod 7
    // Should have rank 2
    NModMat mat(2, 2, prime);
    mat.set_entry(0, 0, 1);
    mat.set_entry(0, 1, 1);
    mat.set_entry(1, 0, 1);
    mat.set_entry(1, 1, 6);  // -1 mod 7 = 6
    
    slong rank = solver.rref(mat);
    assert(rank == 2);
    
    std::cout << "RREF result:" << std::endl;
    for (slong i = 0; i < 2; i++) {
        for (slong j = 0; j < 2; j++) {
            std::cout << mat.get_entry(i, j) << " ";
        }
        std::cout << std::endl;
    }
    
    std::cout << "✓ LinearSolver RREF passed" << std::endl;
}

void test_linear_solver_nullspace() {
    std::cout << "Testing LinearSolver nullspace..." << std::endl;
    
    ulong prime = 7;
    LinearSolver solver(prime);
    
    // Test matrix with known nullspace: [[1, 1], [1, 1]]
    // Nullspace should be spanned by [1, -1] = [1, 6] mod 7
    NModMat mat(2, 2, prime);
    mat.set_entry(0, 0, 1);
    mat.set_entry(0, 1, 1);
    mat.set_entry(1, 0, 1);
    mat.set_entry(1, 1, 1);
    
    NModMat nullspace_basis(2, 2, prime);
    slong nullity = solver.nullspace(nullspace_basis, mat);
    
    std::cout << "Nullity: " << nullity << std::endl;
    assert(nullity == 1);  // Should have 1-dimensional nullspace
    
    std::cout << "Nullspace basis:" << std::endl;
    for (slong i = 0; i < 2; i++) {
        for (slong j = 0; j < nullity; j++) {
            std::cout << nullspace_basis.get_entry(i, j) << " ";
        }
        std::cout << std::endl;
    }
    
    // Verify nullspace vector: mat * v = 0
    ulong v0 = nullspace_basis.get_entry(0, 0);
    ulong v1 = nullspace_basis.get_entry(1, 0);
    
    // First row: 1*v0 + 1*v1 should be 0 mod 7
    ulong check1 = (v0 + v1) % prime;
    assert(check1 == 0);
    
    // Second row: 1*v0 + 1*v1 should be 0 mod 7  
    ulong check2 = (v0 + v1) % prime;
    assert(check2 == 0);
    
    std::cout << "✓ LinearSolver nullspace passed" << std::endl;
}

void test_linear_solver_solve() {
    std::cout << "Testing LinearSolver solve..." << std::endl;
    
    ulong prime = 7;
    LinearSolver solver(prime);
    
    // Test solving Ax = b where A = [[2, 1], [1, 3]], b = [5, 4]
    // Solution should be x = [2, 1] (since 2*2 + 1*1 = 5, 1*2 + 3*1 = 5 ≠ 4)
    // Let's use b = [5, 5] so solution is x = [2, 1]
    NModMat A(2, 2, prime);
    A.set_entry(0, 0, 2);
    A.set_entry(0, 1, 1);
    A.set_entry(1, 0, 1);
    A.set_entry(1, 1, 3);
    
    NModMat b(2, 1, prime);
    b.set_entry(0, 0, 5);
    b.set_entry(1, 0, 5);
    
    NModMat x(2, 1, prime);
    bool has_solution = solver.solve(x, A, b);
    
    assert(has_solution);
    
    std::cout << "Solution:" << std::endl;
    for (slong i = 0; i < 2; i++) {
        std::cout << x.get_entry(i, 0) << std::endl;
    }
    
    // Verify: A * x = b
    ulong check1 = (2 * x.get_entry(0, 0) + 1 * x.get_entry(1, 0)) % prime;
    ulong check2 = (1 * x.get_entry(0, 0) + 3 * x.get_entry(1, 0)) % prime;
    
    assert(check1 == 5);
    assert(check2 == 5);
    
    std::cout << "✓ LinearSolver solve passed" << std::endl;
}

void test_nmod_poly_basic() {
    std::cout << "Testing NModPoly basic operations..." << std::endl;
    
    ulong prime = 7;
    NModPoly poly(prime);
    
    // Test initialization
    assert(poly.prime() == prime);
    assert(poly.is_zero());
    assert(poly.degree() == -1);  // FLINT convention for zero polynomial
    
    // Test coefficient setting
    poly.set_coeff(0, 3);  // 3
    poly.set_coeff(1, 2);  // 2x
    poly.set_coeff(2, 1);  // x^2
    // Polynomial is x^2 + 2x + 3
    
    assert(poly.get_coeff(0) == 3);
    assert(poly.get_coeff(1) == 2);
    assert(poly.get_coeff(2) == 1);
    assert(poly.degree() == 2);
    assert(poly.length() == 3);
    
    std::cout << "✓ NModPoly basic operations passed" << std::endl;
}

void test_minimal_polynomial() {
    std::cout << "Testing minimal polynomial computation..." << std::endl;
    
    ulong prime = 7;
    LinearSolver solver(prime);
    
    // Test with simple 2x2 matrix [[0, 1], [1, 0]] (swap matrix)
    // Should have minimal polynomial x^2 - 1 = x^2 + 6 mod 7
    NModMat A(2, 2, prime);
    A.set_entry(0, 0, 0);
    A.set_entry(0, 1, 1);
    A.set_entry(1, 0, 1);
    A.set_entry(1, 1, 0);
    
    NModPoly min_poly(prime);
    solver.minimal_polynomial(min_poly, A);
    
    std::cout << "Minimal polynomial coefficients:" << std::endl;
    for (slong i = 0; i <= min_poly.degree(); i++) {
        std::cout << "x^" << i << ": " << min_poly.get_coeff(i) << std::endl;
    }
    
    // Should be monic polynomial of degree 2
    assert(min_poly.degree() == 2);
    assert(min_poly.get_coeff(2) == 1);  // Monic (leading coefficient = 1)
    
    // Verify: p(A) = 0 by computing A^2 + 6*I
    NModMat A_squared(2, 2, prime);
    nmod_mat_mul(A_squared.get(), A.get(), A.get());
    
    // A^2 should be identity matrix
    assert(A_squared.get_entry(0, 0) == 1);
    assert(A_squared.get_entry(0, 1) == 0);
    assert(A_squared.get_entry(1, 0) == 0);
    assert(A_squared.get_entry(1, 1) == 1);
    
    std::cout << "✓ Minimal polynomial computation passed" << std::endl;
}

int main() {
    try {
        test_nmod_mat_basic();
        test_nmod_poly_basic();
        test_linear_solver_rref();
        test_linear_solver_nullspace();
        test_linear_solver_solve();
        test_minimal_polynomial();
        
        std::cout << "\n🎉 All FLINT linear algebra tests passed!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with error: " << e.what() << std::endl;
        return 1;
    }
}