#include "axf4_wrapper.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

/* Forward declarations of axf4.c functions */
extern void f4mod_init(long long nvars, long long elim, long long prime);
extern void f4mod_free(void);
extern void f4gb_mod(void);
extern void f4_addrow(void* row);
extern void* getpol(const char* str, int* length);
extern void putpol(void* row, FILE* file);
extern long long f4_aload;  /* number of polynomials in basis */
extern void** f4_array; /* array of basis polynomials */
extern char* vars[1024]; /* variable names - MAXVARS in axf4_lib.c */

/* Global session tracking - F4 only supports one active session at a time */
static axf4_session_t axf4_global_session = NULL;

struct axf4_session {
    int prime;
    char** variables;
    int nvars;
    int initialized;
    char* error_message;
    char** polynomials;
    int num_polynomials;
    int poly_capacity;
};

static void set_error(axf4_session_t session, const char* message) {
    if (session->error_message) {
        free(session->error_message);
    }
    session->error_message = strdup(message);
}

axf4_session_t axf4_create_session(int prime, const char** variables, int nvars) {
    if (prime <= 65536 || nvars <= 0 || nvars > 256) {
        return NULL;
    }
    
    /* Ensure any previous F4 session is properly cleaned up
     * F4 library has global state that persists between sessions
     */
    if (axf4_global_session && axf4_global_session->initialized) {
        f4mod_free();
        axf4_global_session->initialized = 0;
    }
    
    axf4_session_t session = malloc(sizeof(struct axf4_session));
    if (!session) return NULL;
    
    session->prime = prime;
    session->nvars = nvars;
    session->initialized = 0;
    session->error_message = NULL;
    session->num_polynomials = 0;
    session->poly_capacity = 16;
    
    /* Copy variable names */
    session->variables = malloc(nvars * sizeof(char*));
    if (!session->variables) {
        free(session);
        return NULL;
    }
    
    for (int i = 0; i < nvars; i++) {
        session->variables[i] = strdup(variables[i]);
        if (!session->variables[i]) {
            /* Cleanup on error */
            for (int j = 0; j < i; j++) {
                free(session->variables[j]);
            }
            free(session->variables);
            free(session);
            return NULL;
        }
    }
    
    /* Initialize polynomial array */
    session->polynomials = malloc(session->poly_capacity * sizeof(char*));
    if (!session->polynomials) {
        for (int i = 0; i < nvars; i++) {
            free(session->variables[i]);
        }
        free(session->variables);
        free(session);
        return NULL;
    }
    
    return session;
}

int axf4_add_polynomial(axf4_session_t session, const char* polynomial) {
    if (!session || !polynomial) {
        return -1;
    }
    
    /* Expand polynomial array if needed */
    if (session->num_polynomials >= session->poly_capacity) {
        session->poly_capacity *= 2;
        char** new_polys = realloc(session->polynomials, 
                                  session->poly_capacity * sizeof(char*));
        if (!new_polys) {
            set_error(session, "Memory allocation failed");
            return -1;
        }
        session->polynomials = new_polys;
    }
    
    /* Copy polynomial string and add newline if not present */
    size_t len = strlen(polynomial);
    session->polynomials[session->num_polynomials] = malloc(len + 2);
    if (!session->polynomials[session->num_polynomials]) {
        set_error(session, "Memory allocation failed");
        return -1;
    }
    
    strcpy(session->polynomials[session->num_polynomials], polynomial);
    if (len == 0 || polynomial[len-1] != '\n') {
        session->polynomials[session->num_polynomials][len] = '\n';
        session->polynomials[session->num_polynomials][len+1] = '\0';
    }
    
    session->num_polynomials++;
    return 0;
}

/* Common implementation for Gröbner basis computation */
static axf4_result_t axf4_compute_groebner_basis_impl(axf4_session_t session, int keep_data) {
    axf4_result_t result = {0};
    
    if (!session) {
        result.status = -1;
        result.error_message = strdup("Invalid session");
        return result;
    }
    
    if (session->num_polynomials == 0) {
        result.status = -1; 
        result.error_message = strdup("No polynomials added");
        return result;
    }
    
    clock_t start_time = clock();
    
    /* Force cleanup of any previous F4 state before initializing */
    if (axf4_global_session && axf4_global_session->initialized) {
        f4mod_free();
        axf4_global_session->initialized = 0;
        axf4_global_session = NULL;
    }
    
    /* Initialize F4 module */
    f4mod_init(session->nvars, 0, session->prime);
    session->initialized = 1;
    
    /* Track this as the global session */
    axf4_global_session = session;
    
    /* Set variable names in global vars array */
    for (int i = 0; i < session->nvars; i++) {
        vars[i] = session->variables[i];
    }
    
    /* Parse and add all polynomials */
    for (int i = 0; i < session->num_polynomials; i++) {
        int len;
        void* poly = getpol(session->polynomials[i], &len);
        if (!poly) {
            result.status = -1;
            result.error_message = malloc(256);
            snprintf(result.error_message, 256, 
                    "Failed to parse polynomial %d: %s", 
                    i+1, session->polynomials[i]);
            f4mod_free();
            session->initialized = 0;
            return result;
        }
        f4_addrow(poly);
    }
    
    /* Compute Gröbner basis */
    f4gb_mod();
    
    clock_t end_time = clock();
    result.computation_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    result.basis_size = (int)f4_aload;  /* Cast to int for API compatibility */
    
    /* Create output string using memory stream (secure alternative to tmpnam) */
    char* buffer = NULL;
    size_t buffer_size = 0;
    FILE* mem_stream = open_memstream(&buffer, &buffer_size);
    if (!mem_stream) {
        result.status = -1;
        result.error_message = strdup("Failed to create memory stream");
        f4mod_free();
        session->initialized = 0;
        return result;
    }
    
    /* Write all basis polynomials to memory stream */
    for (long long i = 0; i < f4_aload; i++) {
        putpol(f4_array[i], mem_stream);
        fputc('\n', mem_stream);
    }
    fclose(mem_stream);
    
    /* Transfer ownership of buffer to result */
    result.groebner_basis = buffer;
    
    /* Cleanup based on keep_data flag */
    if (!keep_data) {
        f4mod_free();
        session->initialized = 0;
    }
    
    result.status = 0;
    return result;
}

axf4_result_t axf4_compute_groebner_basis(axf4_session_t session) {
    return axf4_compute_groebner_basis_impl(session, 0);
}

void axf4_free_result(axf4_result_t* result) {
    if (result) {
        if (result->groebner_basis) {
            free(result->groebner_basis);
            result->groebner_basis = NULL;
        }
        if (result->error_message) {
            free(result->error_message);
            result->error_message = NULL;
        }
    }
}

void axf4_destroy_session(axf4_session_t session) {
    if (!session) return;
    
    /* Cleanup F4 if still initialized */
    if (session->initialized) {
        f4mod_free();
        session->initialized = 0;  /* Mark as cleaned up */
    }
    
    /* Clear global session if this is it */
    if (axf4_global_session == session) {
        axf4_global_session = NULL;
    }
    
    /* Free variable names */
    if (session->variables) {
        for (int i = 0; i < session->nvars; i++) {
            if (session->variables[i]) {
                free(session->variables[i]);
            }
        }
        free(session->variables);
    }
    
    /* Free polynomial strings */
    if (session->polynomials) {
        for (int i = 0; i < session->num_polynomials; i++) {
            if (session->polynomials[i]) {
                free(session->polynomials[i]);
            }
        }
        free(session->polynomials);
    }
    
    if (session->error_message) {
        free(session->error_message);
    }
    
    free(session);
}

const char* axf4_get_last_error(axf4_session_t session) {
    return session ? session->error_message : NULL;
}

/* ========================================================================
 * STRUCTURED DATA EXTRACTION API IMPLEMENTATION
 * Direct access to F4 internal data structures for RUR integration
 * ======================================================================== */

/* Import F4 internal row structure definition */
typedef struct f4row {
    long long len;    /* number of terms in the matrix row */
    long long sug;    /* phantom sugar degree for the poly */
    long long fac;    /* monomial cofactor from the syzygy */
    long long *cof;   /* array of coefficients modulo p */
    long long *mon;   /* array of monomials from the basis */
    unsigned char *ind; /* array of column index differences */
    long long siz;    /* bytes of encoded sparsity pattern */
} f4row;

int axf4_get_basis_size(void) {
    return (int)f4_aload;  /* Cast to int for API compatibility */
}

int axf4_get_poly_term_count(int poly_index) {
    if (poly_index < 0 || poly_index >= f4_aload) {
        return -1;  /* Error: index out of bounds */
    }
    if (!f4_array || !f4_array[poly_index]) {
        return -1;  /* Error: invalid data */
    }
    
    f4row* poly = (f4row*)f4_array[poly_index];
    return (int)poly->len;
}

int axf4_get_poly_data(int poly_index, unsigned int* out_coeffs, unsigned int* out_monomials) {
    if (poly_index < 0 || poly_index >= f4_aload) {
        return -1;  /* Error: index out of bounds */
    }
    if (!f4_array || !f4_array[poly_index]) {
        return -1;  /* Error: invalid data */
    }
    if (!out_coeffs || !out_monomials) {
        return -1;  /* Error: null output buffers */
    }
    
    f4row* poly = (f4row*)f4_array[poly_index];
    
    /* Copy coefficient and monomial data with safe conversion */
    for (long long i = 0; i < poly->len; i++) {
        out_coeffs[i] = (unsigned int)poly->cof[i];
        out_monomials[i] = (unsigned int)poly->mon[i];
    }
    
    return 0;  /* Success */
}

int axf4_get_leading_term(int poly_index, unsigned int* out_lead_coeff, unsigned int* out_lead_monomial) {
    if (poly_index < 0 || poly_index >= f4_aload) {
        return -1;  /* Error: index out of bounds */
    }
    if (!f4_array || !f4_array[poly_index]) {
        return -1;  /* Error: invalid data */
    }
    if (!out_lead_coeff || !out_lead_monomial) {
        return -1;  /* Error: null output pointers */
    }
    
    f4row* poly = (f4row*)f4_array[poly_index];
    if (poly->len == 0) {
        return -1;  /* Error: empty polynomial has no leading term */
    }
    
    /* The leading term is the first term (F4 stores highest degree first) */
    *out_lead_coeff = (unsigned int)poly->cof[0];
    *out_lead_monomial = (unsigned int)poly->mon[0];
    
    return 0;  /* Success */
}

int axf4_get_all_leading_monomials(unsigned int* out_leading_monomials) {
    if (!out_leading_monomials) {
        return -1;  /* Error: null output buffer */
    }
    if (f4_aload <= 0 || !f4_array) {
        return -1;  /* Error: no basis computed */
    }
    
    int count = 0;
    for (long long i = 0; i < f4_aload; i++) {
        if (!f4_array[i]) {
            continue;  /* Skip null entries */
        }
        
        f4row* poly = (f4row*)f4_array[i];
        if (poly->len > 0) {
            /* Leading term should be the first term after sorting */
            out_leading_monomials[count] = (unsigned int)poly->mon[0];
            count++;
        }
    }
    
    return count;  /* Number of leading monomials written */
}

axf4_result_t axf4_compute_groebner_basis_keep_data(axf4_session_t session) {
    /* Use common implementation with keep_data=1 */
    return axf4_compute_groebner_basis_impl(session, 1);
}

void axf4_cleanup_basis_data(void) {
    // Note: This should only be called after compute_groebner_basis_keep_data()
    // and before destroy_session()
    f4mod_free();
    
    // Mark all sessions as not initialized to prevent double free
    // This is a temporary fix - ideally we'd track which session owns the F4 state
    // For now, we assume only one session is active at a time
    if (axf4_global_session && axf4_global_session->initialized) {
        axf4_global_session->initialized = 0;
    }
}

int axf4_get_num_variables(axf4_session_t session) {
    if (!session) {
        return -1;
    }
    return session->nvars;
}