/*
 * Test different polynomial formats to see what F4 parser accepts
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// F4 function declarations
extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void* getpol(const char* str, int* length);
extern char* vars[256];

void test_format(const char* poly_str, int prime) {
    printf("Testing format: '%s'\n", poly_str);
    
    vars[0] = "x";
    vars[1] = "y";
    
    f4mod_init(2, 0, prime);
    
    int len;
    void* poly = getpol(poly_str, &len);
    
    if (!poly || len == 0) {
        printf("  FAILED to parse (returned %s, len=%d)\n", 
               poly ? "non-null" : "null", len);
    } else {
        printf("  SUCCESS: parsed with len=%d\n", len);
    }
    
    f4mod_free();
}

int main() {
    printf("F4 Parser Format Test\n");
    printf("=====================\n\n");
    
    int prime = 1073741827;
    
    // Test various formats
    printf("Test 1: Standard mathematical notation\n");
    test_format("x^2 + y^2 - 1\n", prime);
    
    printf("\nTest 2: Without spaces\n");
    test_format("x^2+y^2-1\n", prime);
    
    printf("\nTest 3: With explicit multiplication\n");
    test_format("1*x^2 + 1*y^2 - 1\n", prime);
    
    printf("\nTest 4: Leading negative\n");
    test_format("-1 + x^2 + y^2\n", prime);
    
    printf("\nTest 5: Our formatter output\n");
    test_format("-1+1*x^2+1*y^2\n", prime);
    
    printf("\nTest 6: With modular arithmetic\n");
    test_format("1073741826 + x^2 + y^2\n", prime);
    
    printf("\nTest 7: Simple monomial\n");
    test_format("x^2\n", prime);
    
    printf("\nTest 8: Constant\n");
    test_format("1\n", prime);
    
    printf("\nTest 9: Negative constant\n");
    test_format("-1\n", prime);
    
    printf("\nTest 10: F4 internal format guess\n");
    test_format("x^2 + y^2 + 1073741826\n", prime);
    
    return 0;
}