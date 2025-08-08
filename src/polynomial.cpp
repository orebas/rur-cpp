#include "polynomial.hpp"
#include <sstream>
#include <algorithm>

template<typename CoeffT>
MultivariatePolynomial<CoeffT> MultivariatePolynomial<CoeffT>::operator+(const MultivariatePolynomial& other) const {
    MultivariatePolynomial result(comparator_.order);
    
    for (const auto& [mon, coeff] : terms_) {
        result.terms_[mon] = coeff;
    }
    
    for (const auto& [mon, coeff] : other.terms_) {
        CoeffT new_coeff = result.get_coefficient(mon) + coeff;
        result.set_coefficient(mon, new_coeff);
    }
    
    return result;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT> MultivariatePolynomial<CoeffT>::operator-(const MultivariatePolynomial& other) const {
    MultivariatePolynomial result(comparator_.order);
    
    for (const auto& [mon, coeff] : terms_) {
        result.terms_[mon] = coeff;
    }
    
    for (const auto& [mon, coeff] : other.terms_) {
        CoeffT new_coeff = result.get_coefficient(mon) - coeff;
        result.set_coefficient(mon, new_coeff);
    }
    
    return result;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT> MultivariatePolynomial<CoeffT>::operator*(const MultivariatePolynomial& other) const {
    MultivariatePolynomial result(comparator_.order);
    
    for (const auto& [mon1, coeff1] : terms_) {
        for (const auto& [mon2, coeff2] : other.terms_) {
            Monomial product_mon = mon1 * mon2;
            CoeffT product_coeff = coeff1 * coeff2;
            CoeffT new_coeff = result.get_coefficient(product_mon) + product_coeff;
            result.set_coefficient(product_mon, new_coeff);
        }
    }
    
    return result;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT> MultivariatePolynomial<CoeffT>::operator*(const CoeffT& scalar) const {
    if (scalar == CoeffT(0)) {
        return MultivariatePolynomial(comparator_.order);
    }
    
    MultivariatePolynomial result(comparator_.order);
    for (const auto& [mon, coeff] : terms_) {
        result.terms_[mon] = coeff * scalar;
    }
    
    return result;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT>& MultivariatePolynomial<CoeffT>::operator+=(const MultivariatePolynomial& other) {
    for (const auto& [mon, coeff] : other.terms_) {
        CoeffT new_coeff = get_coefficient(mon) + coeff;
        set_coefficient(mon, new_coeff);
    }
    return *this;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT>& MultivariatePolynomial<CoeffT>::operator-=(const MultivariatePolynomial& other) {
    for (const auto& [mon, coeff] : other.terms_) {
        CoeffT new_coeff = get_coefficient(mon) - coeff;
        set_coefficient(mon, new_coeff);
    }
    return *this;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT>& MultivariatePolynomial<CoeffT>::operator*=(const MultivariatePolynomial& other) {
    *this = *this * other;
    return *this;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT>& MultivariatePolynomial<CoeffT>::operator*=(const CoeffT& scalar) {
    if (scalar == CoeffT(0)) {
        terms_.clear();
    } else {
        for (auto& [mon, coeff] : terms_) {
            coeff *= scalar;
        }
    }
    return *this;
}

template<typename CoeffT>
std::string MultivariatePolynomial<CoeffT>::to_string(const std::vector<std::string>& varnames) const {
    if (is_zero()) {
        return "0";
    }
    
    std::ostringstream oss;
    bool first = true;
    
    for (auto it = terms_.rbegin(); it != terms_.rend(); ++it) {
        const auto& [mon, coeff] = *it;
        
        if (!first) {
            if (coeff > CoeffT(0)) {
                oss << "+";
            }
        }
        first = false;
        
        if (coeff == CoeffT(-1) && mon.degree() > 0) {
            oss << "-";
        } else if (coeff != CoeffT(1) || mon.degree() == 0) {
            oss << coeff;
            if (mon.degree() > 0) {
                oss << "*";
            }
        }
        
        if (mon.degree() > 0) {
            oss << mon.to_string(varnames);
        }
    }
    
    return oss.str();
}

template<typename CoeffT>
nlohmann::json MultivariatePolynomial<CoeffT>::to_json() const {
    nlohmann::json j = nlohmann::json::array();
    
    for (const auto& [monomial, coeff] : terms_) {
        nlohmann::json term;
        term["coefficient"] = coeff;
        term["exponents"] = monomial.exponents();
        j.push_back(term);
    }
    
    return j;
}

template<typename CoeffT>
MultivariatePolynomial<CoeffT> MultivariatePolynomial<CoeffT>::from_json(const nlohmann::json& j) {
    MultivariatePolynomial<CoeffT> result;
    
    for (const auto& term : j) {
        CoeffT coeff = term["coefficient"].get<CoeffT>();
        std::vector<int> exponents = term["exponents"].get<std::vector<int>>();
        
        Monomial mon(exponents);
        result.set_coefficient(mon, coeff);
    }
    
    return result;
}

template class MultivariatePolynomial<int>;
template class MultivariatePolynomial<long>;
template class MultivariatePolynomial<double>;