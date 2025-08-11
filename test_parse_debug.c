#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void* getpol(const char* str, int* length);
extern char* vars[256];
extern int f4mon_new(int* exponents);

// Reimplement getmon to debug
extern int getint(char* s, int* l);
extern char* getvar(char* s, int* l);

void debug_parse(const char* poly_str) {
    printf("\n=== Parsing: %s", poly_str);
    printf("=== Character by character:\n");
    for (int i = 0; poly_str[i] && poly_str[i] != '\n'; i++) {
        printf("  [%d] = '%c' (ASCII %d)\n", i, poly_str[i], (int)poly_str[i]);
    }
    
    vars[0] = "x";
    vars[1] = "y";
    
    f4mod_init(2, 0, 1073741827);
    
    int len;
    void* poly = getpol(poly_str, &len);
    
    printf("Result: getpol returned %s, consumed %d chars\n", 
           poly ? "non-null" : "null", len);
    
    f4mod_free();
}

int main() {
    printf("Testing F4 parser behavior\n");
    printf("===========================\n");
    
    // Test different formats
    debug_parse("-1+1*x^2+1*y^2\n");
    debug_parse("-1+x^2+y^2\n");
    debug_parse("x^2+y^2-1\n");
    debug_parse("1073741826+x^2+y^2\n");
    
    return 0;
}