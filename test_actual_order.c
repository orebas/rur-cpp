#include <stdio.h>
#include <stdlib.h>

// Build a standalone F4 test without our wrapper
extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void f4gb_mod(void);
extern void f4_addrow(void* row);
extern void* getpol(const char* str, int* length);
extern void putpol(void* row, FILE* file);
extern int f4_aload;
extern void** f4_array;
extern char* vars[256];

typedef struct f4row {
    long long len, sug, fac;
    long long *cof, *mon;
    unsigned char *ind;
    long long siz;
} f4row;

// Try to use F4's actual monomial decoding
extern int f4mon_deg(int m);  // Get degree of monomial

int main() {
    vars[0] = "x";
    vars[1] = "y";
    
    f4mod_init(2, 0, 1073741827);
    
    // Simple test: x^3 + x
    const char* poly = "x^3+x\n";
    printf("Input: %s", poly);
    
    int len;
    void* p = getpol(poly, &len);
    f4_addrow(p);
    f4gb_mod();
    
    printf("GB has %d poly(s)\n", f4_aload);
    
    f4row* gb = (f4row*)f4_array[0];
    printf("Poly 0: ");
    putpol(f4_array[0], stdout);
    printf("\n");
    printf("  Has %lld terms\n", gb->len);
    
    for (int i = 0; i < gb->len; i++) {
        int mon_idx = gb->mon[i];
        long long coeff = gb->cof[i];
        
        // Try to get degree
        int deg = f4mon_deg(mon_idx);
        printf("  Term[%d]: coeff=%lld, mon_idx=%d, degree=%d\n", 
               i, coeff, mon_idx, deg);
    }
    
    printf("\nConclusion: ");
    if (gb->len > 1) {
        int first_deg = f4mon_deg(gb->mon[0]);
        int last_deg = f4mon_deg(gb->mon[gb->len-1]);
        printf("First term deg=%d, Last term deg=%d\n", first_deg, last_deg);
        printf("Leading term (highest degree) is at: %s\n",
               last_deg > first_deg ? "LAST" : "FIRST");
    }
    
    f4mod_free();
    return 0;
}
