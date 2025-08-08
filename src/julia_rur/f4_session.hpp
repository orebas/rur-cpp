#ifndef F4_SESSION_HPP
#define F4_SESSION_HPP

#include "axf4_wrapper.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>

namespace julia_rur {

// RAII wrapper for the axf4_session_t
class F4Session {
public:
    // Constructor creates the session, throws on failure
    F4Session(int prime, const std::vector<std::string>& variables) {
        std::vector<const char*> var_ptrs;
        var_ptrs.reserve(variables.size());
        for (const auto& var : variables) {
            var_ptrs.push_back(var.c_str());
        }

        session_ = axf4_create_session(prime, var_ptrs.data(), variables.size());
        if (!session_) {
            throw std::runtime_error("Failed to create axf4 session.");
        }
    }

    // Destructor guarantees that the session is destroyed
    ~F4Session() {
        if (session_) {
            axf4_destroy_session(session_);
        }
    }

    // Delete copy constructor and assignment operator to prevent double-free
    F4Session(const F4Session&) = delete;
    F4Session& operator=(const F4Session&) = delete;

    // Move constructor and assignment for safe transfers of ownership
    F4Session(F4Session&& other) noexcept : session_(other.session_) {
        other.session_ = nullptr;
    }

    F4Session& operator=(F4Session&& other) noexcept {
        if (this != &other) {
            if (session_) {
                axf4_destroy_session(session_);
            }
            session_ = other.session_;
            other.session_ = nullptr;
        }
        return *this;
    }

    void add_polynomial(const std::string& polynomial) {
        if (axf4_add_polynomial(session_, polynomial.c_str()) != 0) {
            throw std::runtime_error(std::string("Failed to add polynomial: ") + polynomial);
        }
    }

    axf4_result_t compute_groebner_basis() {
        return axf4_compute_groebner_basis(session_);
    }

    axf4_result_t compute_groebner_basis_keep_data() {
        return axf4_compute_groebner_basis_keep_data(session_);
    }
    
    // This function can now be a static member or a free function,
    // as it operates on global state, not a specific session instance.
    static void cleanup_basis_data() {
        axf4_cleanup_basis_data();
    }

    // Expose the raw handle if needed for other C API calls
    axf4_session_t get() const {
        return session_;
    }

    // Helper to get number of variables
    int get_num_variables() const {
        return axf4_get_num_variables(session_);
    }

private:
    axf4_session_t session_ = nullptr;
};

} // namespace julia_rur

#endif // F4_SESSION_HPP