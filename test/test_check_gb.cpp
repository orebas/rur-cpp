#include "../src/axf4_wrapper.h"
#include <iostream>

int main() {
    const char* vars[] = {"x"};
    
    axf4_session_t session = axf4_create_session(131063, vars, 1);
    axf4_add_polynomial(session, "1*x^2-2");
    
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    std::cout << "Gröbner basis: " << result.groebner_basis << std::endl;
    std::cout << "Status: " << result.status << std::endl;
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
    
    return 0;
}