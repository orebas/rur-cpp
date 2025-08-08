#ifndef AXF4_WRAPPER_H
#define AXF4_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Handle for F4 computation session
 */
typedef struct axf4_session* axf4_session_t;

/**
 * @brief Result structure for F4 computation
 */
typedef struct {
    char* groebner_basis;    /**< String representation of Gröbner basis */
    int status;              /**< 0 for success, non-zero for error */
    char* error_message;     /**< Error description if status != 0 */
    double computation_time; /**< Time taken in seconds */
    int basis_size;          /**< Number of polynomials in basis */
} axf4_result_t;

/**
 * @brief Create a new F4 computation session
 * @param prime The prime modulus to use (must be > 2^16)
 * @param variables Array of variable names in the form ["x0", "x1", ...]
 * @param nvars Number of variables
 * @return New session handle or NULL on error
 */
axf4_session_t axf4_create_session(int prime, const char** variables, int nvars);

/**
 * @brief Add polynomial to the session
 * @param session The session handle
 * @param polynomial String representation of polynomial (e.g., "x0^2 + 2*x0*x1 - 1")
 * @return 0 on success, non-zero on error
 */
int axf4_add_polynomial(axf4_session_t session, const char* polynomial);

/**
 * @brief Compute Gröbner basis
 * @param session The session handle
 * @return Result structure containing the basis
 */
axf4_result_t axf4_compute_groebner_basis(axf4_session_t session);

/**
 * @brief Free result structure
 * @param result The result to free
 */
void axf4_free_result(axf4_result_t* result);

/**
 * @brief Destroy session and free all resources
 * @param session The session to destroy
 */
void axf4_destroy_session(axf4_session_t session);

/**
 * @brief Get last error message for a session
 * @param session The session handle
 * @return Error message string or NULL
 */
const char* axf4_get_last_error(axf4_session_t session);

/* ========================================================================
 * STRUCTURED DATA EXTRACTION API
 * These functions provide direct access to computed Groebner basis data
 * without string parsing, for efficient integration with RUR algorithms.
 * ======================================================================== */

/**
 * @brief Get the number of polynomials in the computed Groebner basis
 * @return Number of basis polynomials, or -1 if no basis computed
 */
int axf4_get_basis_size(void);

/**
 * @brief Get the number of terms in a specific basis polynomial
 * @param poly_index Index of polynomial (0 to axf4_get_basis_size()-1)
 * @return Number of terms, or -1 if index out of bounds
 */
int axf4_get_poly_term_count(int poly_index);

/**
 * @brief Extract coefficient and monomial data for a polynomial
 * @param poly_index Index of polynomial
 * @param out_coeffs Buffer for coefficients (must have space for term_count elements)
 * @param out_monomials Buffer for monomials (must have space for term_count elements)
 * @return 0 on success, -1 on error
 */
int axf4_get_poly_data(int poly_index, unsigned int* out_coeffs, unsigned int* out_monomials);

/**
 * @brief Get the leading term of a polynomial
 * @param poly_index Index of polynomial
 * @param out_lead_coeff Pointer to store leading coefficient
 * @param out_lead_monomial Pointer to store leading monomial
 * @return 0 on success, -1 on error (out of bounds or empty polynomial)
 */
int axf4_get_leading_term(int poly_index, unsigned int* out_lead_coeff, unsigned int* out_lead_monomial);

/**
 * @brief Get all leading monomials for quotient basis extraction
 * @param out_leading_monomials Buffer for leading monomials (must have space for basis_size elements)
 * @return Number of leading monomials written, or -1 on error
 */
int axf4_get_all_leading_monomials(unsigned int* out_leading_monomials);

/**
 * @brief Compute Gröbner basis but keep internal data alive for structured access
 * @param session The session handle
 * @return Result structure containing the basis (internal data remains valid)
 */
axf4_result_t axf4_compute_groebner_basis_keep_data(axf4_session_t session);

/**
 * @brief Clean up internal F4 data after structured access is complete
 * Call this after you're done using the structured API functions
 */
void axf4_cleanup_basis_data(void);

/**
 * @brief Get number of variables in session
 * @param session The session handle
 * @return Number of variables, or -1 on error
 */
int axf4_get_num_variables(axf4_session_t session);

#ifdef __cplusplus
}
#endif

#endif /* AXF4_WRAPPER_H */