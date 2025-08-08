#include <stdio.h>
#include <stdlib.h>

/* External declarations from axf4_lib.c */
extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void f4gb_mod(void);
extern void f4_addrow(void* row);
extern void* getpol(const char* str, int* length);
extern void putpol(void* row, FILE* file);
extern int f4_aload;
extern void** f4_array;
extern char* vars[256];
extern int num_thread;

int main() {
    printf("Testing F4 with minimal setup\n");
    
    // Force single thread
    num_thread = 1;
    
    // Test case that crashes
    long long prime = 1073741827;
    const char* poly = "x - 3\n";
    
    printf("Initializing F4 with prime=%lld\n", prime);
    vars[0] = "x";
    f4mod_init(1, 0, prime);
    
    printf("Parsing polynomial: %s", poly);
    int len;
    void* parsed_poly = getpol(poly, &len);
    if (!parsed_poly) {
        printf("Failed to parse polynomial\n");
        return 1;
    }
    
    printf("Adding polynomial to F4\n");
    f4_addrow(parsed_poly);
    
    printf("Computing Gröbner basis...\n");
    f4gb_mod();
    
    printf("Result: %d polynomials in basis\n", f4_aload);
    for (int i = 0; i < f4_aload; i++) {
        printf("Polynomial %d: ", i);
        putpol(f4_array[i], stdout);
        printf("\n");
    }
    
    f4mod_free();
    printf("Done!\n");
    
    return 0;
}