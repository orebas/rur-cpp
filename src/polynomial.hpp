#pragma once

#include "monomial.hpp"
#include <map>
#include <string>
#include <ostream>
#include <vector>
#include <nlohmann/json.hpp>

template<typename CoeffT>
class MultivariatePolynomial {
private:
    std::map<Monomial, CoeffT, MonomialComparator> terms_;
    MonomialComparator comparator_;
    
public:
    explicit MultivariatePolynomial(MonomialOrder order = MonomialOrder::GREVLEX)
        : comparator_(order) {}
    
    explicit MultivariatePolynomial(const CoeffT& constant, 
                                  size_t nvars = 0,
                                  MonomialOrder order = MonomialOrder::GREVLEX)
        : comparator_(order) {
        if (constant != CoeffT(0)) {
            terms_[Monomial(nvars)] = constant;
        }
    }
    
    void set_coefficient(const Monomial& mon, const CoeffT& coeff) {
        if (coeff == CoeffT(0)) {
            terms_.erase(mon);
        } else {
            terms_[mon] = coeff;
        }
    }
    
    CoeffT get_coefficient(const Monomial& mon) const {
        auto it = terms_.find(mon);
        return (it != terms_.end()) ? it->second : CoeffT(0);
    }
    
    bool is_zero() const { return terms_.empty(); }
    
    size_t nterms() const { return terms_.size(); }
    
    Monomial leading_monomial() const {
        if (is_zero()) {
            throw std::runtime_error("Leading monomial of zero polynomial");
        }
        return terms_.rbegin()->first;
    }
    
    CoeffT leading_coefficient() const {
        if (is_zero()) {
            return CoeffT(0);
        }
        return terms_.rbegin()->second;
    }
    
    int degree() const {
        if (is_zero()) return -1;
        return leading_monomial().degree();
    }
    
    bool is_constant() const {
        if (is_zero()) return true;
        if (terms_.size() != 1) return false;
        return terms_.begin()->first.degree() == 0;
    }
    
    MultivariatePolynomial operator+(const MultivariatePolynomial& other) const;
    MultivariatePolynomial operator-(const MultivariatePolynomial& other) const;
    MultivariatePolynomial operator*(const MultivariatePolynomial& other) const;
    MultivariatePolynomial operator*(const CoeffT& scalar) const;
    
    MultivariatePolynomial& operator+=(const MultivariatePolynomial& other);
    MultivariatePolynomial& operator-=(const MultivariatePolynomial& other);
    MultivariatePolynomial& operator*=(const MultivariatePolynomial& other);
    MultivariatePolynomial& operator*=(const CoeffT& scalar);
    
    bool operator==(const MultivariatePolynomial& other) const {
        return terms_ == other.terms_;
    }
    
    bool operator!=(const MultivariatePolynomial& other) const {
        return !(*this == other);
    }
    
    std::string to_string(const std::vector<std::string>& varnames = {}) const;
    
    static MultivariatePolynomial parse(const std::string& str, 
                                      const std::vector<std::string>& varnames,
                                      MonomialOrder order = MonomialOrder::GREVLEX);
    
    // JSON serialization
    nlohmann::json to_json() const;
    static MultivariatePolynomial from_json(const nlohmann::json& j);
    
    auto begin() const { return terms_.begin(); }
    auto end() const { return terms_.end(); }
    auto rbegin() const { return terms_.rbegin(); }
    auto rend() const { return terms_.rend(); }
    
    template<typename T>
    friend std::ostream& operator<<(std::ostream& os, const MultivariatePolynomial<T>& p);
};

template<typename CoeffT>
std::ostream& operator<<(std::ostream& os, const MultivariatePolynomial<CoeffT>& p) {
    return os << p.to_string();
}