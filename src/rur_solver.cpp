#include "rur_solver.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>

template<typename CoeffT>
void RURSolver<CoeffT>::add_polynomial(const MultivariatePolynomial<CoeffT>& poly) {
    if (gb_computed_) {
        throw std::runtime_error("Cannot add polynomials after Gröbner basis is computed");
    }
    
    if (poly.is_zero()) {
        return; // Skip zero polynomials
    }
    
    input_system_.push_back(poly);
    f4_solver_->add_polynomial(poly);
}

template<typename CoeffT>
void RURSolver<CoeffT>::add_polynomials(const std::vector<MultivariatePolynomial<CoeffT>>& polys) {
    for (const auto& poly : polys) {
        add_polynomial(poly);
    }
}

template<typename CoeffT>
void RURSolver<CoeffT>::compute_groebner_basis() {
    if (input_system_.empty()) {
        throw std::runtime_error("No polynomials in the system");
    }
    
    if (gb_computed_) {
        return; // Already computed
    }
    
    std::cout << "Computing Gröbner basis for " << input_system_.size() << " polynomials..." << std::endl;
    
    // Compute GB using F4
    f4_solver_->compute_groebner_basis();
    
    // Process results
    process_groebner_basis();
    
    gb_computed_ = true;
    std::cout << "Gröbner basis computed with " << groebner_basis_.size() << " elements" << std::endl;
}

template<typename CoeffT>
void RURSolver<CoeffT>::process_groebner_basis() {
    groebner_basis_ = f4_solver_->get_groebner_basis();
    
    // Set GB in quotient ring
    quotient_ring_->set_groebner_basis(groebner_basis_);
}

template<typename CoeffT>
bool RURSolver<CoeffT>::check_zero_dimensional() const {
    if (!gb_computed_) {
        return false;
    }
    
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

template<typename CoeffT>
void RURSolver<CoeffT>::validate_system() const {
    if (!gb_computed_) {
        throw std::runtime_error("Gröbner basis must be computed first");
    }
    
    if (!is_zero_dimensional()) {
        throw std::runtime_error("System is not zero-dimensional - RUR not applicable");
    }
    
    if (groebner_basis_.empty()) {
        throw std::runtime_error("Empty Gröbner basis");
    }
    
    // Check if GB contains the constant polynomial 1 (inconsistent system)
    for (const auto& poly : groebner_basis_) {
        if (poly.is_constant() && poly.leading_coefficient() != 0) {
            throw std::runtime_error("System is inconsistent (contains non-zero constant)");
        }
    }
}

template<typename CoeffT>
void RURSolver<CoeffT>::compute_rur() {
    if (!gb_computed_) {
        compute_groebner_basis();
    }
    
    if (rur_computed_) {
        return; // Already computed
    }
    
    std::cout << "Computing RUR..." << std::endl;
    
    // Validate system
    validate_system();
    
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

template<typename CoeffT>
void RURSolver<CoeffT>::solve() {
    compute_groebner_basis();
    compute_rur();
}

template<typename CoeffT>
const std::vector<MultivariatePolynomial<CoeffT>>& RURSolver<CoeffT>::get_groebner_basis() const {
    if (!gb_computed_) {
        throw std::runtime_error("Gröbner basis not computed yet");
    }
    return groebner_basis_;
}

template<typename CoeffT>
const std::vector<MultivariatePolynomial<CoeffT>>& RURSolver<CoeffT>::get_parameterization() const {
    if (!rur_computed_) {
        throw std::runtime_error("RUR not computed yet");
    }
    return parameterization_;
}

template<typename CoeffT>
const std::vector<CoeffT>& RURSolver<CoeffT>::get_minimal_polynomial() const {
    if (!rur_computed_) {
        throw std::runtime_error("RUR not computed yet");
    }
    return minimal_polynomial_;
}

template<typename CoeffT>
size_t RURSolver<CoeffT>::get_solution_count() const {
    if (!rur_computed_) {
        throw std::runtime_error("RUR not computed yet");
    }
    return quotient_ring_->dimension();
}

template<typename CoeffT>
bool RURSolver<CoeffT>::is_zero_dimensional() const {
    if (!gb_computed_) {
        return false;
    }
    return check_zero_dimensional();
}

template<typename CoeffT>
double RURSolver<CoeffT>::get_gb_computation_time() const {
    if (!gb_computed_) {
        return 0.0;
    }
    return f4_solver_->get_computation_time();
}

template<typename CoeffT>
size_t RURSolver<CoeffT>::get_gb_size() const {
    if (!gb_computed_) {
        return 0;
    }
    return groebner_basis_.size();
}

template<typename CoeffT>
void RURSolver<CoeffT>::set_separating_variable(size_t var_index) {
    if (var_index >= nvars_) {
        throw std::invalid_argument("Variable index out of range");
    }
    
    if (rur_computed_) {
        throw std::runtime_error("Cannot change separating variable after RUR is computed");
    }
    
    bivariate_algo_->set_separating_variable(var_index);
}

template<typename CoeffT>
void RURSolver<CoeffT>::clear() {
    // Clear F4 solver
    f4_solver_->clear();
    
    // Recreate components
    quotient_ring_ = std::make_unique<QuotientRing<CoeffT>>(prime_, nvars_);
    bivariate_algo_ = std::make_unique<BivariateAlgorithm<CoeffT>>(*quotient_ring_, prime_);
    
    // Clear data
    input_system_.clear();
    groebner_basis_.clear();
    parameterization_.clear();
    minimal_polynomial_.clear();
    
    gb_computed_ = false;
    rur_computed_ = false;
}

// Explicit template instantiations
template class RURSolver<int>;
template class RURSolver<long>;