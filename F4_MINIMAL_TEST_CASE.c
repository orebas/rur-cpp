/*
 * Minimal test case showing inconsistent F4 behavior across different primes
 * 
 * For the polynomial system {x^2 - 1 = 0, y - x = 0}, the Gröbner basis
 * should contain y^2 - 1, which is represented as y^2 + (p-1) in modular arithmetic.
 * 
 * However, F4 produces different results for different primes:
 * - Some primes: Correctly produces y^2 + (p-1)
 * - Most primes: Incorrectly produces y^2 + 1
 */

#include <stdio.h>
#include <string.h>

// F4 function declarations (from axf4_lib.c)
extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void f4gb_mod(void);
extern void f4_addrow(void* row);
extern void* getpol(const char* str, int* length);
extern void putpol(void* row, FILE* file);
extern int f4_aload;
extern void** f4_array;
extern char* vars[256];

void test_prime(long long prime) {
    printf("\nTesting prime %lld:\n", prime);
    
    // Initialize F4
    vars[0] = "x";
    vars[1] = "y";
    f4mod_init(2, 0, prime);
    
    // Add polynomials: x^2 - 1 and y - x
    int len;
    void* poly1 = getpol("1*x^2-1\n", &len);
    void* poly2 = getpol("1*y-1*x\n", &len);
    
    if (!poly1 || !poly2) {
        printf("Failed to parse polynomials\n");
        f4mod_free();
        return;
    }
    
    f4_addrow(poly1);
    f4_addrow(poly2);
    
    // Compute Gröbner basis
    f4gb_mod();
    
    // Print results
    printf("Gröbner basis:\n");
    for (int i = 0; i < f4_aload; i++) {
        putpol(f4_array[i], stdout);
        printf("\n");
    }
    
    // Check for y^2 term in output
    char buffer[1024];
    FILE* memstream = fmemopen(buffer, sizeof(buffer), "w");
    for (int i = 0; i < f4_aload; i++) {
        putpol(f4_array[i], memstream);
    }
    fclose(memstream);
    
    // Analyze result
    char expected[32];
    sprintf(expected, "y^2+%lld", prime - 1);
    
    if (strstr(buffer, expected)) {
        printf("✓ CORRECT: Found %s (represents y^2 - 1)\n", expected);
    } else if (strstr(buffer, "y^2+1")) {
        printf("✗ WRONG: Found y^2+1 (should be y^2+%lld)\n", prime - 1);
    } else {
        printf("? UNEXPECTED: y^2 term not found or has unexpected coefficient\n");
    }
    
    f4mod_free();
}

int main() {
    printf("F4 Inconsistency Test Case\n");
    printf("==========================\n");
    printf("Testing polynomial system: {x^2 - 1 = 0, y - x = 0}\n");
    printf("Expected: Gröbner basis contains y^2 + (p-1)\n");
    
    // Test various primes
    test_prime(65537);      // 2^16 + 1 (works correctly)
    test_prime(1073741827); // 2^30 + 3 (works correctly)
    test_prime(1048573);    // Fails
    test_prime(268435399);  // Fails
    
    printf("\nSummary: F4 produces different Gröbner bases for the same polynomial system\n");
    printf("when using different primes. Is this expected behavior?\n");
    
    return 0;
}