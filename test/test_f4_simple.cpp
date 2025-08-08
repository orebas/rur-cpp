#include "f4_solver.hpp"
#include <iostream>
#include <cassert>

int main() {
    try {
        std::cout << "Creating simple F4 solver...\n";
        
        std::vector<std::string> vars = {"x"};
        F4Solver<int> solver(100003, vars);
        
        std::cout << "Solver created successfully\n";
        std::cout << "Prime: " << solver.prime() << std::endl;
        std::cout << "Variables: ";
        for (const auto& var : solver.variable_names()) {
            std::cout << var << " ";
        }
        std::cout << std::endl;
        
        // Try a very simple polynomial: x - 1
        std::cout << "Creating polynomial...\n" << std::flush;
        MultivariatePolynomial<int> p;
        std::cout << "Setting x coefficient...\n" << std::flush;
        Monomial x_mon(1);  // Create monomial with 1 variable
        x_mon[0] = 1;       // Set x^1
        p.set_coefficient(x_mon, 1);
        std::cout << "Setting constant coefficient...\n" << std::flush;
        Monomial const_mon(1);  // Constant term (all exponents 0)
        p.set_coefficient(const_mon, -1);
        std::cout << "Polynomial created\n" << std::flush;
        
        std::cout << "Polynomial has " << p.nterms() << " terms\n";
        for (const auto& [mon, coeff] : p) {
            std::cout << "  Monomial exponents: [" << mon[0] << "]"
                      << " (degree " << mon.degree() << "), coeff: " << coeff << std::endl;
        }
        std::cout << "Polynomial: " << p.to_string(vars) << std::endl;
        
        // Test the axf4 conversion directly
        std::cout << "Testing axf4 format conversion...\n";
        
        std::cout << "Adding polynomial to solver...\n";
        solver.add_polynomial(p);
        
        std::cout << "Computing Gröbner basis...\n";
        solver.compute_groebner_basis();
        
        const auto& basis = solver.get_groebner_basis();
        std::cout << "Basis computed with " << basis.size() << " polynomials\n";
        
        // Add assertions to verify the result
        assert(basis.size() == 1);  // Should have exactly one polynomial in the basis
        assert(basis[0].nterms() == 2);  // Should have 2 terms (x and constant)
        
        for (size_t i = 0; i < basis.size(); i++) {
            std::cout << "g" << i+1 << " = " << basis[i].to_string(vars) << std::endl;
        }
        
        std::cout << "Test completed successfully!\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}