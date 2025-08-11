#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void f4gb_mod(void);
extern void f4_addrow(void* row);
extern void* getpol(const char* str, int* length);
extern void putpol(void* row, FILE* file);
extern int f4_aload;
extern void** f4_array;
extern char* vars[256];

// Import decoder function
extern int f4_nvars;
extern int* f4_monom;
#define F4_EXPON(m) (f4_monom + (m)*(f4_nvars+1))

typedef struct f4row {
    long long len, sug, fac;
    long long *cof, *mon;
    unsigned char *ind;
    long long siz;
} f4row;

int main() {
    vars[0] = "x";
    vars[1] = "y";
    
    f4mod_init(2, 0, 1073741827);
    
    // Parse x^2 + y^2 - 1
    const char* poly_str = "-1+1*x^2+1*y^2\n";
    printf("Parsing: %s", poly_str);
    
    int len;
    void* poly = getpol(poly_str, &len);
    f4_addrow(poly);
    f4gb_mod();
    
    printf("GB computed, %d polynomials\n\n", f4_aload);
    
    // Get the first GB polynomial
    f4row* gb = (f4row*)f4_array[0];
    printf("GB poly 0 has %lld terms\n", gb->len);
    printf("String form: ");
    putpol(f4_array[0], stdout);
    printf("\n\n");
    
    // Now let's decode each term
    printf("Decoding terms:\n");
    for (int i = 0; i < gb->len && i < 5; i++) {
        unsigned int mon_idx = gb->mon[i];
        long long coeff = gb->cof[i];
        
        printf("Term[%d]: coeff=%lld, mon_idx=%u\n", i, coeff, mon_idx);
        
        // Try to decode it
        int* exps = F4_EXPON(mon_idx);
        printf("  Exponents at %p: ", (void*)exps);
        for (int j = 0; j < f4_nvars; j++) {
            printf("var[%d]^%d ", j, exps[j]);
        }
        printf("\n");
        
        // Also show raw memory
        printf("  Raw bytes: ");
        unsigned char* bytes = (unsigned char*)exps;
        for (int j = 0; j < (f4_nvars+1)*sizeof(int); j++) {
            printf("%02x ", bytes[j]);
        }
        printf("\n");
    }
    
    // Also check our monomial table
    printf("\nMonomial table entries:\n");
    for (int i = 0; i < 10; i++) {
        int* exps = F4_EXPON(i);
        printf("  [%d]: ", i);
        for (int j = 0; j < f4_nvars; j++) {
            printf("v%d^%d ", j, exps[j]);
        }
        printf("\n");
    }
    
    f4mod_free();
    return 0;
}
