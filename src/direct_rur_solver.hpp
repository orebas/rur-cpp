#pragma once

#include "quotient_ring.hpp"
#include "bivariate_lex.hpp"
#include "polynomial.hpp"
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>

/**
 * Direct RUR solver that works with polynomial objects directly,
 * bypassing F4 string parsing issues. For algorithm development and testing.
 */
template<typename CoeffT>
class DirectRURSolver {
private:
    int prime_;
    std::vector<std::string> variable_names_;
    size_t nvars_;
    
    // Core components
    std::unique_ptr<QuotientRing<CoeffT>> quotient_ring_;
    std::unique_ptr<BivariateAlgorithm<CoeffT>> bivariate_algo_;
    
    // Input system and results
    std::vector<MultivariatePolynomial<CoeffT>> input_system_;
    std::vector<MultivariatePolynomial<CoeffT>> groebner_basis_;
    std::vector<MultivariatePolynomial<CoeffT>> parameterization_;
    std::vector<CoeffT> minimal_polynomial_;
    
    bool gb_computed_;
    bool rur_computed_;
    
public:
    DirectRURSolver(int prime, const std::vector<std::string>& variable_names)
        : prime_(prime), variable_names_(variable_names), nvars_(variable_names.size()),
          gb_computed_(false), rur_computed_(false) {
        
        if (prime <= 65536) {
            throw std::invalid_argument("Prime must be > 65536");
        }
        
        if (nvars_ == 0 || nvars_ > 256) {
            throw std::invalid_argument("Number of variables must be between 1 and 256");
        }
        
        // Initialize components
        quotient_ring_ = std::make_unique<QuotientRing<CoeffT>>(prime, nvars_);
        bivariate_algo_ = std::make_unique<BivariateAlgorithm<CoeffT>>(*quotient_ring_, prime);
    }
    
    // Add polynomial directly
    void add_polynomial(const MultivariatePolynomial<CoeffT>& poly) {
        if (gb_computed_) {
            throw std::runtime_error("Cannot add polynomials after Gröbner basis is set");
        }
        
        if (!poly.is_zero()) {
            input_system_.push_back(poly);
        }
    }
    
    // Set precomputed Gröbner basis directly (for testing)
    void set_groebner_basis(const std::vector<MultivariatePolynomial<CoeffT>>& basis) {
        groebner_basis_ = basis;
        quotient_ring_->set_groebner_basis(basis);
        gb_computed_ = true;
    }
    
    // Compute RUR from current GB
    void compute_rur() {
        if (!gb_computed_) {
            throw std::runtime_error("Gröbner basis must be set first");
        }
        
        if (rur_computed_) {
            return; // Already computed
        }
        
        std::cout << "Computing RUR directly..." << std::endl;
        
        // Step 1: Compute quotient basis
        std::cout << "Computing quotient basis..." << std::endl;
        quotient_ring_->compute_quotient_basis();
        
        // Step 2: Prepare multiplication tables
        std::cout << "Preparing multiplication tables..." << std::endl;
        quotient_ring_->prepare_multiplication_tables();
        
        // Step 3: Compute multiplication tables
        std::cout << "Computing multiplication tables..." << std::endl;
        quotient_ring_->compute_multiplication_tables();
        
        // Step 4: Compute separating element
        std::cout << "Computing separating element..." << std::endl;
        bivariate_algo_->compute_separating_element();
        minimal_polynomial_ = bivariate_algo_->get_minimal_polynomial();
        
        // Step 5: Extract parameterization
        std::cout << "Extracting parameterization..." << std::endl;
        parameterization_ = bivariate_algo_->extract_parameterization();
        
        rur_computed_ = true;
        
        std::cout << "RUR computation complete!" << std::endl;
        std::cout << "Quotient dimension: " << quotient_ring_->dimension() << std::endl;
        std::cout << "Separating element degree: " << bivariate_algo_->get_separating_degree() << std::endl;
    }
    
    // Get results
    const std::vector<MultivariatePolynomial<CoeffT>>& get_groebner_basis() const {
        if (!gb_computed_) {
            throw std::runtime_error("Gröbner basis not set yet");
        }
        return groebner_basis_;
    }
    
    const std::vector<MultivariatePolynomial<CoeffT>>& get_parameterization() const {
        if (!rur_computed_) {
            throw std::runtime_error("RUR not computed yet");
        }
        return parameterization_;
    }
    
    const std::vector<CoeffT>& get_minimal_polynomial() const {
        if (!rur_computed_) {
            throw std::runtime_error("RUR not computed yet");
        }
        return minimal_polynomial_;
    }
    
    size_t get_solution_count() const {
        if (!rur_computed_) {
            throw std::runtime_error("RUR not computed yet");
        }
        return quotient_ring_->dimension();
    }
    
    bool is_zero_dimensional() const {
        if (!gb_computed_) {
            return false;
        }
        return check_zero_dimensional();
    }
    
    void set_separating_variable(size_t var_index) {
        if (var_index >= nvars_) {
            throw std::invalid_argument("Variable index out of range");
        }
        
        if (rur_computed_) {
            throw std::runtime_error("Cannot change separating variable after RUR is computed");
        }
        
        bivariate_algo_->set_separating_variable(var_index);
    }
    
private:
    bool check_zero_dimensional() const {
        // Check if for each variable, there exists a univariate polynomial in the GB
        std::vector<bool> has_univariate(nvars_, false);
        
        for (const auto& poly : groebner_basis_) {
            // Count non-zero variables in leading monomial
            auto leading_mon = poly.leading_monomial();
            size_t non_zero_vars = 0;
            size_t last_var = 0;
            
            for (size_t i = 0; i < nvars_; ++i) {
                if (leading_mon[i] > 0) {
                    non_zero_vars++;
                    last_var = i;
                }
            }
            
            if (non_zero_vars == 1) {
                has_univariate[last_var] = true;
            }
        }
        
        // System is zero-dimensional if we have univariate polynomials for all variables
        return std::all_of(has_univariate.begin(), has_univariate.end(), [](bool b) { return b; });
    }
};

// Explicit template instantiation declarations
extern template class DirectRURSolver<int>;
extern template class DirectRURSolver<long>;