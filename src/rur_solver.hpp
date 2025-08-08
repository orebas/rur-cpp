#pragma once

#include "f4_solver.hpp"
#include "quotient_ring.hpp"
#include "bivariate_lex.hpp"
#include <vector>
#include <string>
#include <memory>

/**
 * Rational Univariate Representation solver
 * Implements the Julia RUR algorithm in C++
 */
template<typename CoeffT>
class RURSolver {
private:
    int prime_;
    std::vector<std::string> variable_names_;
    size_t nvars_;
    
    // Core components
    std::unique_ptr<F4Solver<CoeffT>> f4_solver_;
    std::unique_ptr<QuotientRing<CoeffT>> quotient_ring_;
    std::unique_ptr<BivariateAlgorithm<CoeffT>> bivariate_algo_;
    
    // Input system
    std::vector<MultivariatePolynomial<CoeffT>> input_system_;
    
    // Results
    std::vector<MultivariatePolynomial<CoeffT>> groebner_basis_;
    std::vector<MultivariatePolynomial<CoeffT>> parameterization_;
    std::vector<CoeffT> minimal_polynomial_;
    
    bool gb_computed_;
    bool rur_computed_;
    
public:
    RURSolver(int prime, const std::vector<std::string>& variable_names)
        : prime_(prime), variable_names_(variable_names), nvars_(variable_names.size()),
          gb_computed_(false), rur_computed_(false) {
        
        if (prime <= 65536) {
            throw std::invalid_argument("Prime must be > 65536");
        }
        
        if (nvars_ == 0 || nvars_ > 256) {
            throw std::invalid_argument("Number of variables must be between 1 and 256");
        }
        
        // Initialize components
        f4_solver_ = std::make_unique<F4Solver<CoeffT>>(prime, variable_names);
        quotient_ring_ = std::make_unique<QuotientRing<CoeffT>>(prime, nvars_);
        bivariate_algo_ = std::make_unique<BivariateAlgorithm<CoeffT>>(*quotient_ring_, prime);
    }
    
    // Add polynomial to the system
    void add_polynomial(const MultivariatePolynomial<CoeffT>& poly);
    
    // Add multiple polynomials
    void add_polynomials(const std::vector<MultivariatePolynomial<CoeffT>>& polys);
    
    // Compute Gröbner basis
    void compute_groebner_basis();
    
    // Compute RUR (requires GB to be computed first)
    void compute_rur();
    
    // Complete solve: GB + RUR in one call
    void solve();
    
    // Get results
    const std::vector<MultivariatePolynomial<CoeffT>>& get_groebner_basis() const;
    const std::vector<MultivariatePolynomial<CoeffT>>& get_parameterization() const;
    const std::vector<CoeffT>& get_minimal_polynomial() const;
    
    // Get quotient ring dimension (number of solutions)
    size_t get_solution_count() const;
    
    // Check if system is zero-dimensional
    bool is_zero_dimensional() const;
    
    // Get computation statistics
    double get_gb_computation_time() const;
    size_t get_gb_size() const;
    
    // Set separating variable (default is last variable)
    void set_separating_variable(size_t var_index);
    
    // Clear all data and reset
    void clear();
    
    // Validate that the system is suitable for RUR
    void validate_system() const;
    
private:
    // Helper: convert F4 output to internal format
    void process_groebner_basis();
    
    // Helper: check if system is zero-dimensional using GB
    bool check_zero_dimensional() const;
};

// Explicit template instantiation declarations
extern template class RURSolver<int>;
extern template class RURSolver<long>;