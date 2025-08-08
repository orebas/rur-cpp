#include "f4_solver.hpp"
#include <sstream>
#include <regex>
#include <algorithm>
#include <iostream>

template<typename CoeffT>
F4Solver<CoeffT>::F4Solver(int prime, const std::vector<std::string>& variable_names)
    : session_(nullptr), variable_names_(variable_names), prime_(prime),
      basis_computed_(false), computation_time_(0.0) {
    
    if (variable_names.empty() || variable_names.size() > 256) {
        throw std::invalid_argument("Number of variables must be between 1 and 256");
    }
    
    if (prime <= 65536) {
        throw std::invalid_argument("Prime must be > 65536");
    }
    
    // Convert to C-style strings for axf4
    std::vector<const char*> var_ptrs;
    for (const auto& name : variable_names_) {
        var_ptrs.push_back(name.c_str());
    }
    
    session_ = axf4_create_session(prime, var_ptrs.data(), variable_names_.size());
    if (!session_) {
        throw std::runtime_error("Failed to create F4 session");
    }
}

template<typename CoeffT>
F4Solver<CoeffT>::~F4Solver() {
    if (session_) {
        axf4_destroy_session(session_);
    }
}

template<typename CoeffT>
F4Solver<CoeffT>::F4Solver(F4Solver&& other) noexcept
    : session_(other.session_), variable_names_(std::move(other.variable_names_)),
      prime_(other.prime_), input_polynomials_(std::move(other.input_polynomials_)),
      groebner_basis_(std::move(other.groebner_basis_)),
      basis_computed_(other.basis_computed_), computation_time_(other.computation_time_) {
    other.session_ = nullptr;
}

template<typename CoeffT>
F4Solver<CoeffT>& F4Solver<CoeffT>::operator=(F4Solver&& other) noexcept {
    if (this != &other) {
        if (session_) {
            axf4_destroy_session(session_);
        }
        
        session_ = other.session_;
        variable_names_ = std::move(other.variable_names_);
        prime_ = other.prime_;
        input_polynomials_ = std::move(other.input_polynomials_);
        groebner_basis_ = std::move(other.groebner_basis_);
        basis_computed_ = other.basis_computed_;
        computation_time_ = other.computation_time_;
        
        other.session_ = nullptr;
    }
    return *this;
}

template<typename CoeffT>
std::string F4Solver<CoeffT>::polynomial_to_axf4_string(const MultivariatePolynomial<CoeffT>& poly) const {
    if (poly.is_zero()) {
        return "0";
    }
    
    std::ostringstream oss;
    bool first = true;
    
    // Need to iterate in reverse order (highest degree first) for axf4 format
    for (auto it = poly.rbegin(); it != poly.rend(); ++it) {
        const auto& [mon, coeff] = *it;
        
        // Convert coefficient to positive modular representation first
        CoeffT mod_coeff = coeff;
        if (mod_coeff < 0) {
            mod_coeff = mod_coeff % prime_;
            if (mod_coeff < 0) mod_coeff += prime_;
        } else {
            mod_coeff = mod_coeff % prime_;
        }
        
        // Always add + for non-first terms (axf4 expects this format)
        if (!first) {
            oss << "+";
        }
        first = false;
        
        // Always output coefficient with *
        oss << mod_coeff;
        
        // Add monomial part
        bool has_vars = false;
        for (size_t i = 0; i < variable_names_.size(); i++) {
            int exp = mon[i];
            if (exp > 0) {
                oss << "*" << variable_names_[i];
                if (exp > 1) {
                    oss << "^" << exp;
                }
                has_vars = true;
            }
        }
    }
    
    return oss.str();
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT> F4Solver<CoeffT>::axf4_string_to_polynomial(const std::string& str) const {
    MultivariatePolynomial<CoeffT> result;
    
    if (str.empty() || str == "0") {
        return result;
    }
    
    // Parse axf4 format: "+3*x0^2+5*x1*x2-1"
    std::regex term_regex(R"([+-]?\d+(?:\*[a-zA-Z_][a-zA-Z0-9_]*(?:\^?\d+)?)*?)");
    std::sregex_iterator iter(str.begin(), str.end(), term_regex);
    std::sregex_iterator end;
    
    for (; iter != end; ++iter) {
        std::string term = iter->str();
        if (term.empty()) continue;
        
        // Parse coefficient
        CoeffT coeff = 1;
        size_t star_pos = term.find('*');
        if (star_pos != std::string::npos) {
            std::string coeff_str = term.substr(0, star_pos);
            if (coeff_str == "+") coeff = 1;
            else if (coeff_str == "-") coeff = -1;
            else coeff = static_cast<CoeffT>(std::stoll(coeff_str));
        } else {
            // Pure coefficient term
            coeff = static_cast<CoeffT>(std::stoll(term));
        }
        
        // Create monomial
        Monomial mon(variable_names_.size());
        
        if (star_pos != std::string::npos) {
            std::string var_part = term.substr(star_pos + 1);
            
            // Parse variables and exponents
            std::regex var_regex(R"([a-zA-Z_][a-zA-Z0-9_]*(?:\^?\d+)?)");
            std::sregex_iterator var_iter(var_part.begin(), var_part.end(), var_regex);
            std::sregex_iterator var_end;
            
            for (; var_iter != var_end; ++var_iter) {
                std::string var_exp = var_iter->str();
                
                size_t caret_pos = var_exp.find('^');
                std::string var_name;
                int exp = 1;
                
                if (caret_pos != std::string::npos) {
                    var_name = var_exp.substr(0, caret_pos);
                    exp = std::stoi(var_exp.substr(caret_pos + 1));
                } else {
                    var_name = var_exp;
                }
                
                // Find variable index
                auto var_it = std::find(variable_names_.begin(), variable_names_.end(), var_name);
                if (var_it != variable_names_.end()) {
                    size_t var_idx = std::distance(variable_names_.begin(), var_it);
                    mon[var_idx] += exp;
                }
            }
        }
        
        result.set_coefficient(mon, result.get_coefficient(mon) + coeff);
    }
    
    return result;
}

template<typename CoeffT>
void F4Solver<CoeffT>::add_polynomial(const MultivariatePolynomial<CoeffT>& poly) {
    if (basis_computed_) {
        throw std::runtime_error("Cannot add polynomials after basis is computed");
    }
    
    input_polynomials_.push_back(poly);
    
    std::string axf4_str = polynomial_to_axf4_string(poly);
    int result = axf4_add_polynomial(session_, axf4_str.c_str());
    
    if (result != 0) {
        const char* error = axf4_get_last_error(session_);
        throw std::runtime_error(std::string("Failed to add polynomial: ") + 
                                (error ? error : "unknown error"));
    }
}

template<typename CoeffT>
void F4Solver<CoeffT>::add_polynomials(const std::vector<MultivariatePolynomial<CoeffT>>& polys) {
    for (const auto& poly : polys) {
        add_polynomial(poly);
    }
}

template<typename CoeffT>
void F4Solver<CoeffT>::compute_groebner_basis() {
    if (input_polynomials_.empty()) {
        throw std::runtime_error("No polynomials to compute basis for");
    }
    
    axf4_result_t result = axf4_compute_groebner_basis(session_);
    
    if (result.status != 0) {
        std::string error_msg = "F4 computation failed";
        if (result.error_message) {
            error_msg += ": ";
            error_msg += result.error_message;
        }
        axf4_free_result(&result);
        throw std::runtime_error(error_msg);
    }
    
    computation_time_ = result.computation_time;
    
    // Parse result into C++ polynomials
    groebner_basis_.clear();
    if (result.groebner_basis) {
        std::istringstream iss(result.groebner_basis);
        std::string line;
        
        while (std::getline(iss, line)) {
            if (!line.empty()) {
                try {
                    MultivariatePolynomial<CoeffT> poly = axf4_string_to_polynomial(line);
                    if (!poly.is_zero()) {
                        groebner_basis_.push_back(poly);
                    }
                } catch (const std::exception& e) {
                    axf4_free_result(&result);
                    throw std::runtime_error(std::string("Failed to parse basis polynomial: ") + e.what());
                }
            }
        }
    }
    
    basis_computed_ = true;
    axf4_free_result(&result);
}

template<typename CoeffT>
const std::vector<MultivariatePolynomial<CoeffT>>& F4Solver<CoeffT>::get_groebner_basis() const {
    if (!basis_computed_) {
        throw std::runtime_error("Gröbner basis not computed yet");
    }
    return groebner_basis_;
}

template<typename CoeffT>
void F4Solver<CoeffT>::clear() {
    if (session_) {
        axf4_destroy_session(session_);
        
        // Create new session
        std::vector<const char*> var_ptrs;
        for (const auto& name : variable_names_) {
            var_ptrs.push_back(name.c_str());
        }
        
        session_ = axf4_create_session(prime_, var_ptrs.data(), variable_names_.size());
        if (!session_) {
            throw std::runtime_error("Failed to recreate F4 session");
        }
    }
    
    input_polynomials_.clear();
    groebner_basis_.clear();
    basis_computed_ = false;
    computation_time_ = 0.0;
}

// Explicit template instantiations
template class F4Solver<int>;
template class F4Solver<long>;