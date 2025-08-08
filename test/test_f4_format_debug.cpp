#include "../src/axf4_wrapper.h"
#include <iostream>

int main() {
    std::cout << "Testing polynomial format issue..." << std::endl;
    
    const char* vars[] = {"x"};
    axf4_session_t session = axf4_create_session(100003, vars, 1);
    if (!session) {
        std::cerr << "Failed to create session" << std::endl;
        return 1;
    }
    
    // Try different polynomial formats:
    std::cout << "Trying format 1: '1*x+100002' (F4Solver format for x-1)" << std::endl;
    int ret = axf4_add_polynomial(session, "1*x+100002");
    if (ret != 0) {
        std::cerr << "Failed to add polynomial: " << axf4_get_last_error(session) << std::endl;
        axf4_destroy_session(session);
        return 1;
    }
    
    std::cout << "Computing basis..." << std::endl;
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    if (result.status != 0) {
        std::cerr << "Computation failed: " << result.error_message << std::endl;
        axf4_free_result(&result);
        axf4_destroy_session(session);
        return 1;
    }
    
    std::cout << "Success! Format worked." << std::endl;
    std::cout << "Result: " << (result.groebner_basis ? result.groebner_basis : "NULL") << std::endl;
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
    return 0;
}