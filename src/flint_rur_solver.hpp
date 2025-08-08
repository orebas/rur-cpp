#pragma once

#include "flint_linear_algebra.hpp"
#include "flint_mpoly_wrappers.hpp"
#include "flint_reducer_v2.hpp"
#include "polynomial_basis.hpp"
#include "polynomial.hpp"
#include "monomial.hpp"
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <iostream>

using namespace flint_linalg;
using namespace flint_mpoly;

/**
 * FLINT-based Rational Univariate Representation solver
 * Uses FLINT's finite field arithmetic for production-quality RUR computation
 */
template<typename CoeffT>
class FLINTRURSolver {
private:
    ulong prime_;
    std::vector<std::string> variable_names_;
    slong nvars_;
    
    // FLINT context and linear algebra solver
    std::shared_ptr<MpolyContext> poly_ctx_;
    std::unique_ptr<LinearSolver> linear_solver_;
    
    // Polynomial basis and reducer
    std::shared_ptr<const PolynomialBasis<CoeffT>> basis_;
    std::unique_ptr<FLINTReducer<CoeffT>> reducer_;
    
    // Input system and results
    std::vector<MultivariatePolynomial<CoeffT>> input_system_;
    std::vector<MultivariatePolynomial<CoeffT>> groebner_basis_;
    
    // RUR results
    NModPoly minimal_polynomial_;
    std::vector<NModPoly> variable_numerators_;
    NModPoly common_denominator_;
    
    // Computation state
    bool gb_computed_;
    bool rur_computed_;
    slong separating_var_index_;
    
public:
    FLINTRURSolver(ulong prime, const std::vector<std::string>& variable_names)
        : prime_(prime), variable_names_(variable_names), nvars_(variable_names.size()),
          minimal_polynomial_(prime), common_denominator_(prime),
          gb_computed_(false), rur_computed_(false), separating_var_index_(nvars_ - 1) {
        
        if (prime < 2) {
            throw std::invalid_argument("Prime must be >= 2");
        }
        
        if (nvars_ <= 0 || nvars_ > 256) {
            throw std::invalid_argument("Number of variables must be between 1 and 256");
        }
        
        // Initialize FLINT context (degree reverse lexicographic ordering)
        poly_ctx_ = std::make_shared<MpolyContext>(nvars_, prime, ORD_DEGREVLEX);
        linear_solver_ = std::make_unique<LinearSolver>(prime);
        
        // Initialize variable numerators
        variable_numerators_.reserve(nvars_);
        for (slong i = 0; i < nvars_; i++) {
            variable_numerators_.emplace_back(prime);
        }
    }
    
    // Add polynomial to the system
    void add_polynomial(const MultivariatePolynomial<CoeffT>& poly) {
        if (gb_computed_) {
            throw std::runtime_error("Cannot add polynomials after Gröbner basis is computed");
        }
        
        if (!poly.is_zero()) {
            input_system_.push_back(poly);
        }
    }
    
    // Add multiple polynomials
    void add_polynomials(const std::vector<MultivariatePolynomial<CoeffT>>& polys) {
        for (const auto& poly : polys) {
            add_polynomial(poly);
        }
    }
    
    // Set Gröbner basis directly (when computed externally)
    void set_groebner_basis(const std::vector<MultivariatePolynomial<CoeffT>>& gb) {
        groebner_basis_ = gb;
        
        // Create polynomial basis
        basis_ = PolynomialBasis<CoeffT>::create(gb, nvars_, prime_, ORD_DEGREVLEX);
        reducer_ = std::make_unique<FLINTReducer<CoeffT>>(basis_);
        
        gb_computed_ = true;
    }
    
    // Compute RUR representation
    void compute_rur() {
        if (!gb_computed_) {
            throw std::runtime_error("Gröbner basis must be computed first");
        }
        
        if (rur_computed_) {
            return; // Already computed
        }
        
        std::cout << "Computing RUR with " << basis_->quotient_basis_monomials.size() 
                  << " basis elements..." << std::endl;
        
        // Step 1: Compute minimal polynomial of separating element
        compute_minimal_polynomial();
        
        // Step 2: Compute rational representations for each variable
        compute_variable_representations();
        
        rur_computed_ = true;
        std::cout << "RUR computation completed" << std::endl;
    }
    
    // Get results
    const std::vector<MultivariatePolynomial<CoeffT>>& get_groebner_basis() const {
        return groebner_basis_;
    }
    
    const NModPoly& get_minimal_polynomial() const {
        if (!rur_computed_) {
            throw std::runtime_error("RUR not computed yet");
        }
        return minimal_polynomial_;
    }
    
    const std::vector<NModPoly>& get_variable_numerators() const {
        if (!rur_computed_) {
            throw std::runtime_error("RUR not computed yet");
        }
        return variable_numerators_;
    }
    
    const NModPoly& get_common_denominator() const {
        if (!rur_computed_) {
            throw std::runtime_error("RUR not computed yet");
        }
        return common_denominator_;
    }
    
    // Get quotient ring dimension
    size_t get_solution_count() const {
        if (!gb_computed_) {
            return 0;
        }
        return basis_->quotient_basis_monomials.size();
    }
    
    // Set separating variable (default is last variable)
    void set_separating_variable(slong var_index) {
        if (var_index < 0 || var_index >= nvars_) {
            throw std::out_of_range("Invalid variable index");
        }
        separating_var_index_ = var_index;
    }
    
    // Get computation info
    ulong prime() const { return prime_; }
    slong nvars() const { return nvars_; }
    const std::vector<std::string>& variable_names() const { return variable_names_; }
    
private:
    /**
     * Compute minimal polynomial following Julia first_variable algorithm
     * Iteratively compute powers of separating element and detect linear dependence
     */
    void compute_minimal_polynomial() {
        size_t dim = basis_->quotient_basis_monomials.size();
        
        std::cout << "Computing minimal polynomial for dimension " << dim << " using Julia first_variable algorithm..." << std::endl;
        
        // Initialize with basis vector for separating variable
        // Start with e_sep_var = [0, 0, ..., 1, ..., 0] where 1 is at position of x_sep_var
        std::vector<ulong> current_power(dim, 0);
        
        // Find which basis element corresponds to x_sep_var
        Monomial sep_var_monomial(nvars_);
        sep_var_monomial[separating_var_index_] = 1;
        
        size_t sep_var_basis_idx = 0;
        for (size_t i = 0; i < dim; i++) {
            if (basis_->quotient_basis_monomials[i] == sep_var_monomial) {
                sep_var_basis_idx = i;
                break;
            }
        }
        current_power[sep_var_basis_idx] = 1;
        
        // Matrix to store powers: [1, t, t^2, t^3, ...]
        std::vector<std::vector<ulong>> power_matrix;
        
        // Add constant 1 = [1, 0, 0, ..., 0]
        std::vector<ulong> constant_vector(dim, 0);
        constant_vector[0] = 1; // Assume constant term is first basis element
        power_matrix.push_back(constant_vector);
        
        // Add t = current_power
        power_matrix.push_back(current_power);
        
        // Iteratively compute t^k and check for linear dependence
        for (size_t k = 2; k <= dim; k++) {
            // Compute t * t^(k-1) by applying multiplication
            MultivariatePolynomial<CoeffT> t_times_prev;
            for (size_t i = 0; i < dim; i++) {
                if (current_power[i] != 0) {
                    const Monomial& basis_mon = basis_->quotient_basis_monomials[i];
                    MultivariatePolynomial<CoeffT> scaled_basis;
                    scaled_basis.set_coefficient(basis_mon, static_cast<CoeffT>(current_power[i]));
                    
                    MultivariatePolynomial<CoeffT> t_times_scaled = compute_separating_element_times_basis(i);
                    for (const auto& [mon, coeff] : t_times_scaled) {
                        CoeffT existing = t_times_prev.get_coefficient(mon);
                        t_times_prev.set_coefficient(mon, existing + coeff * static_cast<CoeffT>(current_power[i]));
                    }
                }
            }
            
            // Reduce modulo ideal
            auto reduced_coeffs = reducer_->reduce(t_times_prev);
            
            // Convert to modular coefficients
            std::vector<ulong> new_power(dim);
            for (size_t i = 0; i < dim; i++) {
                if (i < reduced_coeffs.size()) {
                    ulong coeff_mod = static_cast<ulong>(reduced_coeffs[i]) % prime_;
                    if (reduced_coeffs[i] < 0) {
                        coeff_mod = (coeff_mod + prime_) % prime_;
                    }
                    new_power[i] = coeff_mod;
                } else {
                    new_power[i] = 0;
                }
            }
            
            // Check for linear dependence using Gaussian elimination
            NModMat dep_matrix(dim, k, prime_);
            for (size_t col = 0; col < k; col++) {
                if (col < power_matrix.size()) {
                    for (size_t row = 0; row < dim; row++) {
                        dep_matrix.set_entry(row, col, power_matrix[col][row]);
                    }
                } else {
                    for (size_t row = 0; row < dim; row++) {
                        dep_matrix.set_entry(row, col, new_power[row]);
                    }
                }
            }
            
            slong rank = linear_solver_->rref(dep_matrix);
            
            if (rank < k) {
                // Linear dependence detected - extract minimal polynomial
                std::cout << "Linear dependence detected at power " << k-1 << ", rank = " << rank << std::endl;
                
                // The minimal polynomial coefficients are the null space
                NModMat nullspace_matrix(k, k, prime_);
                slong nullity = linear_solver_->nullspace(nullspace_matrix, dep_matrix);
                
                if (nullity > 0) {
                    minimal_polynomial_.zero();
                    for (slong deg = 0; deg < k; deg++) {
                        ulong coeff = nullspace_matrix.get_entry(deg, 0);
                        minimal_polynomial_.set_coeff(deg, coeff);
                    }
                }
                break;
            }
            
            // Add new power to matrix
            power_matrix.push_back(new_power);
            current_power = new_power;
        }
        
        std::cout << "Minimal polynomial degree: " << minimal_polynomial_.degree() << std::endl;
    }
    
    /**
     * Compute rational representations x_i = p_i(t) / q(t)
     * Solves linear systems to find numerator polynomials
     */
    void compute_variable_representations() {
        size_t dim = basis_->quotient_basis_monomials.size();
        slong min_poly_deg = minimal_polynomial_.degree();
        
        std::cout << "Computing variable representations..." << std::endl;
        
        // Use derivative of minimal polynomial as common denominator (standard in RUR)
        common_denominator_ = std::move(NModPoly(prime_));
        for (slong i = 1; i <= min_poly_deg; i++) {
            ulong deriv_coeff = (static_cast<ulong>(i) * minimal_polynomial_.get_coeff(i)) % prime_;
            common_denominator_.set_coeff(i - 1, deriv_coeff);
        }
        
        // For each variable x_i, solve for numerator polynomial p_i(t)
        for (slong var_idx = 0; var_idx < nvars_; var_idx++) {
            compute_variable_numerator(var_idx, dim, min_poly_deg);
        }
    }
    
    /**
     * Compute numerator polynomial for variable x_i
     */
    void compute_variable_numerator(slong var_idx, size_t dim, slong min_poly_deg) {
        std::cout << "Computing numerator for variable " << var_idx << std::endl;
        
        // Create linear system to solve for p_i(t) such that
        // x_i = p_i(t) / q'(t) where q'(t) is derivative of minimal polynomial
        
        NModMat system_matrix(dim, min_poly_deg, prime_);
        NModMat rhs_vector(dim, 1, prime_);
        system_matrix.zero();
        rhs_vector.zero();
        
        // Build system: p_i(M_t) = x_i applied to quotient basis
        for (size_t basis_idx = 0; basis_idx < dim; basis_idx++) {
            // Compute x_i * basis_element
            MultivariatePolynomial<CoeffT> xi_times_basis = compute_variable_times_basis(var_idx, basis_idx);
            
            // Reduce modulo ideal to get coefficients in quotient basis
            auto xi_reduced = reducer_->reduce(xi_times_basis);
            
            // Set RHS: coefficient of this basis element after reduction
            if (basis_idx < xi_reduced.size()) {
                ulong coeff_mod = static_cast<ulong>(xi_reduced[basis_idx]) % prime_;
                if (xi_reduced[basis_idx] < 0) {
                    coeff_mod = (coeff_mod + prime_) % prime_;
                }
                rhs_vector.set_entry(basis_idx, 0, coeff_mod);
            } else {
                rhs_vector.set_entry(basis_idx, 0, 0);
            }
            
            // Fill system matrix columns: t^deg applied to this basis element
            for (slong deg = 0; deg < min_poly_deg; deg++) {
                MultivariatePolynomial<CoeffT> t_power_basis = compute_separating_power_times_basis(deg, basis_idx);
                auto t_power_reduced = reducer_->reduce(t_power_basis);
                
                // Coefficient of this basis element in the reduction
                ulong matrix_entry = 0;
                if (basis_idx < t_power_reduced.size()) {
                    matrix_entry = static_cast<ulong>(t_power_reduced[basis_idx]) % prime_;
                    if (t_power_reduced[basis_idx] < 0) {
                        matrix_entry = (matrix_entry + prime_) % prime_;
                    }
                }
                system_matrix.set_entry(basis_idx, deg, matrix_entry);
            }
        }
        
        std::cout << "System matrix " << dim << "x" << min_poly_deg << " created" << std::endl;
        
        // Debug: Print system matrix rank
        NModMat test_matrix(dim, min_poly_deg, prime_);
        for (size_t i = 0; i < dim; i++) {
            for (slong j = 0; j < min_poly_deg; j++) {
                test_matrix.set_entry(i, j, system_matrix.get_entry(i, j));
            }
        }
        slong rank = linear_solver_->rref(test_matrix);
        std::cout << "System matrix rank: " << rank << " (should be " << min_poly_deg << ")" << std::endl;
        
        // Solve the linear system
        NModMat solution(min_poly_deg, 1, prime_);
        bool solved = linear_solver_->solve(solution, system_matrix, rhs_vector);
        
        if (!solved) {
            std::cout << "Linear system failed for variable " << var_idx << std::endl;
            // Try least squares or pseudoinverse approach
            throw std::runtime_error("Failed to solve for variable representation");
        }
        
        // Store solution as polynomial coefficients
        variable_numerators_[var_idx].zero();
        for (slong deg = 0; deg < min_poly_deg; deg++) {
            ulong coeff = solution.get_entry(deg, 0);
            variable_numerators_[var_idx].set_coeff(deg, coeff);
        }
        
        std::cout << "Variable " << var_idx << " numerator computed successfully" << std::endl;
    }
    
    /**
     * Helper: compute separating element times basis element
     * Following Julia implementation: use last variable as separating element
     */
    MultivariatePolynomial<CoeffT> compute_separating_element_times_basis(size_t basis_idx) {
        const Monomial& basis_monomial = basis_->quotient_basis_monomials[basis_idx];
        MultivariatePolynomial<CoeffT> result;
        
        // Use last variable as separating element (Julia default)
        slong sep_var = separating_var_index_;
        
        // Create monomial by multiplying basis monomial by x_sep_var
        std::vector<int> new_exps = basis_monomial.exponents();
        if (sep_var < new_exps.size()) {
            new_exps[sep_var]++;
        } else {
            new_exps.resize(sep_var + 1, 0);
            new_exps[sep_var] = 1;
        }
        
        Monomial new_monomial(new_exps);
        result.set_coefficient(new_monomial, static_cast<CoeffT>(1));
        
        return result;
    }
    
    /**
     * Helper: compute variable x_i times basis element
     */
    MultivariatePolynomial<CoeffT> compute_variable_times_basis(slong var_idx, size_t basis_idx) {
        const Monomial& basis_monomial = basis_->quotient_basis_monomials[basis_idx];
        MultivariatePolynomial<CoeffT> result;
        
        // Create monomial by multiplying basis monomial by x_var_idx
        std::vector<int> new_exps = basis_monomial.exponents();
        if (var_idx < new_exps.size()) {
            new_exps[var_idx]++;
        } else {
            new_exps.resize(var_idx + 1, 0);
            new_exps[var_idx] = 1;
        }
        
        Monomial new_monomial(new_exps);
        result.set_coefficient(new_monomial, static_cast<CoeffT>(1));
        
        return result;
    }
    
    /**
     * Helper: compute t^deg times basis element  
     */
    MultivariatePolynomial<CoeffT> compute_separating_power_times_basis(slong deg, size_t basis_idx) {
        MultivariatePolynomial<CoeffT> result;
        result.set_coefficient(basis_->quotient_basis_monomials[basis_idx], static_cast<CoeffT>(1));
        
        // Apply t multiplication deg times
        for (slong i = 0; i < deg; i++) {
            MultivariatePolynomial<CoeffT> temp_result;
            
            for (const auto& [monomial, coeff] : result) {
                MultivariatePolynomial<CoeffT> t_times_monomial;
                
                // Apply separating element
                for (slong var = 0; var < nvars_; var++) {
                    CoeffT sep_coeff = static_cast<CoeffT>(var + 1);
                    
                    std::vector<int> new_exps = monomial.exponents();
                    if (var < new_exps.size()) {
                        new_exps[var]++;
                    } else {
                        new_exps.resize(var + 1, 0);
                        new_exps[var] = 1;
                    }
                    
                    Monomial new_monomial(new_exps);
                    CoeffT existing_coeff = temp_result.get_coefficient(new_monomial);
                    temp_result.set_coefficient(new_monomial, existing_coeff + coeff * sep_coeff);
                }
            }
            
            result = temp_result;
        }
        
        return result;
    }
};

// Explicit template instantiation declarations
extern template class FLINTRURSolver<int>;
extern template class FLINTRURSolver<long>;