#include "monomial.hpp"
#include <numeric>
#include <sstream>
#include <algorithm>

int Monomial::degree() const {
    return std::accumulate(exponents_.begin(), exponents_.end(), 0);
}

Monomial Monomial::operator*(const Monomial& other) const {
    if (nvars() != other.nvars()) {
        throw std::invalid_argument("Monomials must have same number of variables");
    }
    
    std::vector<int> result_exps(nvars());
    for (size_t i = 0; i < nvars(); ++i) {
        result_exps[i] = exponents_[i] + other.exponents_[i];
    }
    
    return Monomial(result_exps);
}

bool Monomial::operator<(const Monomial& other) const {
    return MonomialComparator(MonomialOrder::GREVLEX)(*this, other);
}

bool Monomial::divides(const Monomial& other) const {
    if (nvars() != other.nvars()) {
        return false;
    }
    
    for (size_t i = 0; i < nvars(); ++i) {
        if (exponents_[i] > other.exponents_[i]) {
            return false;
        }
    }
    return true;
}

Monomial Monomial::quotient(const Monomial& divisor) const {
    if (!divisor.divides(*this)) {
        throw std::invalid_argument("Divisor does not divide this monomial");
    }
    
    std::vector<int> result_exps(nvars());
    for (size_t i = 0; i < nvars(); ++i) {
        result_exps[i] = exponents_[i] - divisor.exponents_[i];
    }
    
    return Monomial(result_exps);
}

std::string Monomial::to_string(const std::vector<std::string>& varnames) const {
    std::ostringstream oss;
    
    bool first = true;
    for (size_t i = 0; i < nvars(); ++i) {
        if (exponents_[i] > 0) {
            if (!first) oss << "*";
            first = false;
            
            if (varnames.empty() || i >= varnames.size()) {
                oss << "x" << i;
            } else {
                oss << varnames[i];
            }
            
            if (exponents_[i] > 1) {
                oss << "^" << exponents_[i];
            }
        }
    }
    
    if (first) {
        oss << "1";
    }
    
    return oss.str();
}

std::ostream& operator<<(std::ostream& os, const Monomial& m) {
    return os << m.to_string();
}

bool MonomialComparator::operator()(const Monomial& a, const Monomial& b) const {
    switch (order) {
        case MonomialOrder::LEX:
            return std::lexicographical_compare(
                a.exponents_.begin(), a.exponents_.end(),
                b.exponents_.begin(), b.exponents_.end(),
                std::greater<int>()
            );
            
        case MonomialOrder::GREVLEX: {
            int deg_a = a.degree();
            int deg_b = b.degree();
            if (deg_a != deg_b) {
                return deg_a < deg_b;
            }
            return std::lexicographical_compare(
                a.exponents_.rbegin(), a.exponents_.rend(),
                b.exponents_.rbegin(), b.exponents_.rend(),
                std::greater<int>()
            );
        }
        
        case MonomialOrder::DEGLEX: {
            int deg_a = a.degree();
            int deg_b = b.degree();
            if (deg_a != deg_b) {
                return deg_a < deg_b;
            }
            return std::lexicographical_compare(
                a.exponents_.begin(), a.exponents_.end(),
                b.exponents_.begin(), b.exponents_.end(),
                std::greater<int>()
            );
        }
    }
    return false;
}