#pragma once

#include "polynomial.hpp"
#include "axf4_wrapper.h"
#include <vector>
#include <string>
#include <chrono>
#include <stdexcept>

/**
 * @brief C++ wrapper for the F4 Gröbner basis algorithm
 */
template<typename CoeffT>
class F4Solver {
private:
    axf4_session_t session_;
    std::vector<std::string> variable_names_;
    int prime_;
    std::vector<MultivariatePolynomial<CoeffT>> input_polynomials_;
    std::vector<MultivariatePolynomial<CoeffT>> groebner_basis_;
    bool basis_computed_;
    double computation_time_;
    
    std::string polynomial_to_axf4_string(const MultivariatePolynomial<CoeffT>& poly) const;
    MultivariatePolynomial<CoeffT> axf4_string_to_polynomial(const std::string& str) const;
    
public:
    /**
     * @brief Construct F4 solver with given prime and variables
     * @param prime The prime modulus (must be > 65536)
     * @param variable_names Names of variables (e.g., {"x", "y", "z"})
     */
    F4Solver(int prime, const std::vector<std::string>& variable_names);
    
    /**
     * @brief Destructor - cleans up session
     */
    ~F4Solver();
    
    // Non-copyable but movable
    F4Solver(const F4Solver&) = delete;
    F4Solver& operator=(const F4Solver&) = delete;
    F4Solver(F4Solver&& other) noexcept;
    F4Solver& operator=(F4Solver&& other) noexcept;
    
    /**
     * @brief Add polynomial to the ideal
     * @param poly The polynomial to add
     */
    void add_polynomial(const MultivariatePolynomial<CoeffT>& poly);
    
    /**
     * @brief Add multiple polynomials
     * @param polys Vector of polynomials to add
     */
    void add_polynomials(const std::vector<MultivariatePolynomial<CoeffT>>& polys);
    
    /**
     * @brief Compute the Gröbner basis
     * @throws std::runtime_error if computation fails
     */
    void compute_groebner_basis();
    
    /**
     * @brief Get the computed Gröbner basis
     * @return Vector of polynomials forming the basis
     * @throws std::runtime_error if basis not computed
     */
    const std::vector<MultivariatePolynomial<CoeffT>>& get_groebner_basis() const;
    
    /**
     * @brief Check if Gröbner basis has been computed
     */
    bool is_basis_computed() const { return basis_computed_; }
    
    /**
     * @brief Get computation time in seconds
     */
    double get_computation_time() const { return computation_time_; }
    
    /**
     * @brief Get number of polynomials in basis
     */
    size_t basis_size() const { return groebner_basis_.size(); }
    
    /**
     * @brief Get variable names
     */
    const std::vector<std::string>& variable_names() const { return variable_names_; }
    
    /**
     * @brief Get prime modulus
     */
    int prime() const { return prime_; }
    
    /**
     * @brief Clear all polynomials and reset
     */
    void clear();
};

// Template instantiations
extern template class F4Solver<int>;
extern template class F4Solver<long>;