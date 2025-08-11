/*
 * Direct test of F4 implementation to isolate bug with x^2 + y^2 - 1
 * This bypasses our wrapper entirely
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Direct F4 function declarations from axf4_lib.c
extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void f4gb_mod(void);
extern void f4_addrow(void* row);
extern void* getpol(const char* str, int* length);
extern void putpol(void* row, FILE* file);
extern int f4_aload;  // number of polynomials in basis
extern void** f4_array; // array of basis polynomials
extern char* vars[256]; // variable names

// F4 row structure
typedef struct f4row {
    long long len;    // number of terms
    long long sug;    // sugar degree
    long long fac;    // monomial cofactor
    long long *cof;   // coefficients
    long long *mon;   // monomials
    unsigned char *ind; // column indices
    long long siz;    // size
} f4row;

void test_circle_with_prime(int prime) {
    printf("\n==========================================\n");
    printf("Testing x^2 + y^2 - 1 with prime p = %d\n", prime);
    printf("Prime mod 4 = %d\n", prime % 4);
    printf("==========================================\n");
    
    // Set up variable names
    vars[0] = "x";
    vars[1] = "y";
    
    // Initialize F4 with 2 variables, no elimination, and the given prime
    printf("Initializing F4 with nvars=2, elim=0, prime=%d\n", prime);
    f4mod_init(2, 0, prime);
    
    // Parse the polynomial x^2 + y^2 - 1
    const char* poly_str = "x^2 + y^2 - 1\n";
    printf("Parsing polynomial: %s", poly_str);
    
    int len;
    void* poly = getpol(poly_str, &len);
    if (!poly) {
        printf("ERROR: Failed to parse polynomial\n");
        f4mod_free();
        return;
    }
    
    printf("Successfully parsed polynomial (length=%d)\n", len);
    
    // Examine the parsed polynomial structure
    f4row* row = (f4row*)poly;
    printf("Parsed polynomial details:\n");
    printf("  Number of terms: %lld\n", row->len);
    printf("  Sugar degree: %lld\n", row->sug);
    
    if (row->len > 0 && row->len <= 10) {
        printf("  Terms:\n");
        for (int i = 0; i < row->len; i++) {
            printf("    Term %d: coeff=%lld, monomial=%lld (encoded)\n", 
                   i, row->cof[i], row->mon[i]);
            // Decode monomial (assuming encoding is: low 15 bits = x exp, high 15 bits = y exp)
            int x_exp = row->mon[i] & 0x7FFF;
            int y_exp = (row->mon[i] >> 15) & 0x7FFF;
            printf("           decoded as x^%d * y^%d\n", x_exp, y_exp);
        }
    }
    
    // Add polynomial to the basis
    printf("\nAdding polynomial to basis\n");
    f4_addrow(poly);
    
    // Compute Groebner basis
    printf("Computing Groebner basis with F4...\n");
    f4gb_mod();
    
    // Check the result
    printf("\nGroebner basis computed:\n");
    printf("Number of polynomials in basis: %d\n", f4_aload);
    
    // Print each polynomial in the basis
    for (int i = 0; i < f4_aload && i < 5; i++) {
        printf("\nBasis polynomial %d:\n", i);
        
        // Print to stdout
        putpol(f4_array[i], stdout);
        printf("\n");
        
        // Examine the structure
        f4row* basis_row = (f4row*)f4_array[i];
        printf("  Structure: %lld terms, sugar=%lld\n", 
               basis_row->len, basis_row->sug);
        
        // Show leading term
        if (basis_row->len > 0) {
            printf("  Leading term: coeff=%lld, monomial=%lld\n",
                   basis_row->cof[0], basis_row->mon[0]);
            int x_exp = basis_row->mon[0] & 0x7FFF;
            int y_exp = (basis_row->mon[0] >> 15) & 0x7FFF;
            printf("  Leading monomial decoded: x^%d * y^%d\n", x_exp, y_exp);
            
            if (x_exp == 0 && y_exp == 0) {
                printf("  WARNING: Leading term is a CONSTANT!\n");
                if (basis_row->len == 1) {
                    printf("  This polynomial is just the constant %lld\n", basis_row->cof[0]);
                }
            }
        }
    }
    
    // Clean up
    printf("\nCleaning up F4\n");
    f4mod_free();
}

int main() {
    printf("Direct F4 Test for x^2 + y^2 - 1\n");
    printf("=================================\n");
    
    // Test with various primes
    int test_primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
                         101, 1009, 10007, 100003,
                         793493497, 1073741827};
    int num_primes = sizeof(test_primes) / sizeof(test_primes[0]);
    
    for (int i = 0; i < num_primes; i++) {
        test_circle_with_prime(test_primes[i]);
    }
    
    return 0;
}