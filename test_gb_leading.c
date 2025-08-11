#include <stdio.h>
#include <stdlib.h>

extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void f4gb_mod(void);
extern void f4_addrow(void* row);
extern void* getpol(const char* str, int* length);
extern void putpol(void* row, FILE* file);
extern int f4_aload;
extern void** f4_array;
extern char* vars[256];
extern int f4_nvars;
extern int* f4_monom;
extern int f4_mload;

typedef struct f4row {
    long long len, sug, fac;
    long long *cof, *mon;
    unsigned char *ind;
    long long siz;
} f4row;

#define F4_EXPON(m) (f4_monom + (m)*(f4_nvars+1))

void test_polynomial(const char* name, const char* poly_str) {
    printf("\n========== Testing: %s ==========\n", name);
    printf("Input: %s", poly_str);
    
    int len;
    void* poly = getpol(poly_str, &len);
    f4_addrow(poly);
    f4gb_mod();
    
    printf("GB computed. %d polynomial(s)\n", f4_aload);
    
    for (int p = 0; p < f4_aload && p < 3; p++) {
        f4row* gb = (f4row*)f4_array[p];
        printf("\nPoly %d: ", p);
        putpol(f4_array[p], stdout);
        printf("\n");
        printf("  %lld terms\n", gb->len);
        
        if (gb->len > 0) {
            // Check first term
            int* exp0 = F4_EXPON(gb->mon[0]);
            int deg0 = exp0[0] + exp0[1];
            printf("  FIRST [0]: mon_idx=%lld, x^%d*y^%d, deg=%d\n", 
                   gb->mon[0], exp0[0], exp0[1], deg0);
            
            // Check last term
            int last = gb->len - 1;
            int* expL = F4_EXPON(gb->mon[last]);
            int degL = expL[0] + expL[1];
            printf("  LAST [%d]: mon_idx=%lld, x^%d*y^%d, deg=%d\n", 
                   last, gb->mon[last], expL[0], expL[1], degL);
            
            printf("  => Leading term is at: %s\n", 
                   degL > deg0 ? "LAST" : deg0 > degL ? "FIRST" : "SAME DEGREE");
        }
    }
    
    // Clear for next test
    f4_aload = 0;
}

int main() {
    vars[0] = "x";
    vars[1] = "y";
    f4mod_init(2, 0, 1073741827);
    
    // Test various polynomials
    test_polynomial("Circle", "x^2+y^2-1\n");
    test_polynomial("Simple univariate", "x^3+2*x+1\n");
    test_polynomial("Mixed degrees", "1+x+y+x^2+x*y+y^2\n");
    
    f4mod_free();
    return 0;
}
