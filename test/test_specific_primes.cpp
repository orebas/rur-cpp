#include "../src/julia_rur/data_structures.hpp"
#include "../src/julia_rur/f4_polynomial_formatter.hpp"
extern "C" {
#include "../src/axf4_wrapper.h"
}
#include <iostream>

using namespace julia_rur;

bool test_prime(ModularCoeff prime) {
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(prime, vars, 2);
    
    axf4_add_polynomial(session, "1*x^2-1");
    axf4_add_polynomial(session, "1*y-1*x");
    
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    std::string gb_str(result.groebner_basis);
    bool correct = gb_str.find(std::to_string(prime - 1)) != std::string::npos;
    
    std::cout << "Prime " << prime << ": ";
    if (correct) {
        std::cout << "CORRECT (contains " << (prime-1) << ")" << std::endl;
    } else {
        std::cout << "WRONG (contains +1)" << std::endl;
    }
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
    return correct;
}

int main() {
    std::cout << "Testing specific 30-bit primes:" << std::endl;
    
    // Test the 30-bit primes from our list
    ModularCoeff primes_30bit[] = {
        1073741827, 1073741831, 1073741833, 1073741839, 1073741843,
        1073741857, 1073741873, 1073741909, 1073741939, 1073741941,
        1073741783, 1073741741, 1073741723, 1073741719, 1073741717
    };
    
    int correct_count = 0;
    for (auto p : primes_30bit) {
        if (test_prime(p)) correct_count++;
    }
    
    std::cout << "\nSummary: " << correct_count << " out of " 
              << (sizeof(primes_30bit)/sizeof(primes_30bit[0])) 
              << " primes gave correct results" << std::endl;
    
    return 0;
}