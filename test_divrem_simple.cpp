#include "src/flint_mpoly_wrappers.hpp"
#include <iostream>

int main() {
    try {
        // Create context with 2 variables, prime 100003
        auto ctx = std::make_shared<const flint_mpoly::MpolyContext>(2, 100003, ORD_DEGREVLEX);
        
        // Create dividend: x^2 + y^2
        flint_mpoly::Mpoly dividend(ctx);
        ulong exp_x2[] = {2, 0}; // x^2
        ulong exp_y2[] = {0, 2}; // y^2
        dividend.set_coeff_ui_monomial(1, exp_x2);
        dividend.set_coeff_ui_monomial(1, exp_y2);
        
        // Create divisors: [x, y]
        std::vector<flint_mpoly::Mpoly> divisors;
        
        flint_mpoly::Mpoly div_x(ctx);
        ulong exp_x[] = {1, 0}; // x
        div_x.set_coeff_ui_monomial(1, exp_x);
        divisors.push_back(std::move(div_x));
        
        flint_mpoly::Mpoly div_y(ctx);
        ulong exp_y[] = {0, 1}; // y
        div_y.set_coeff_ui_monomial(1, exp_y);
        divisors.push_back(std::move(div_y));
        
        // Create quotients vector
        std::vector<flint_mpoly::Mpoly> quotients;
        quotients.reserve(divisors.size());
        for (size_t i = 0; i < divisors.size(); ++i) {
            quotients.emplace_back(ctx);
        }
        
        // Create remainder
        flint_mpoly::Mpoly remainder(ctx);
        
        // Perform division - this should trigger the linking error
        flint_mpoly::divrem_ideal(quotients, remainder, dividend, divisors);
        
        std::cout << "Division completed successfully!" << std::endl;
        std::cout << "Remainder has " << remainder.length() << " terms" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}