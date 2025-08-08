#pragma once

#include <vector>
#include <string>
#include <ostream>
#include <functional>

class Monomial {
private:
    std::vector<int> exponents_;
    
public:
    explicit Monomial(size_t nvars = 0) : exponents_(nvars, 0) {}
    
    explicit Monomial(const std::vector<int>& exponents) : exponents_(exponents) {}
    
    size_t nvars() const { return exponents_.size(); }
    
    int degree() const;
    
    int operator[](size_t i) const { return exponents_[i]; }
    int& operator[](size_t i) { return exponents_[i]; }
    
    Monomial operator*(const Monomial& other) const;
    
    bool operator==(const Monomial& other) const { return exponents_ == other.exponents_; }
    bool operator!=(const Monomial& other) const { return !(*this == other); }
    
    bool operator<(const Monomial& other) const;
    
    bool divides(const Monomial& other) const;
    Monomial quotient(const Monomial& divisor) const;
    
    std::string to_string(const std::vector<std::string>& varnames = {}) const;
    
    // Get exponents vector (for RUR implementation)
    const std::vector<int>& exponents() const { return exponents_; }
    
    friend std::ostream& operator<<(std::ostream& os, const Monomial& m);
    friend struct MonomialComparator;
};

enum class MonomialOrder {
    LEX,
    GREVLEX, 
    DEGLEX
};

struct MonomialComparator {
    MonomialOrder order;
    
    explicit MonomialComparator(MonomialOrder ord = MonomialOrder::GREVLEX) : order(ord) {}
    
    bool operator()(const Monomial& a, const Monomial& b) const;
};

// Hash specialization for std::unordered_map support
namespace std {
    template<>
    struct hash<Monomial> {
        std::size_t operator()(const Monomial& monomial) const noexcept {
            std::size_t result = 0;
            const auto& exps = monomial.exponents();
            
            // Use boost-style hash_combine algorithm
            for (const auto& exp : exps) {
                result ^= std::hash<int>{}(exp) + 0x9e3779b9 + (result << 6) + (result >> 2);
            }
            
            return result;
        }
    };
}