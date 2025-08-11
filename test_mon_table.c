#include <stdio.h>
#include <stdlib.h>

extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4mod_free(void);
extern void* getpol(const char* str, int* length);
extern char* vars[256];
extern int f4_nvars;
extern int* f4_monom;
extern int f4_msize;  // monomial table size

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
    
    // Parse a polynomial to populate monomial table
    const char* ps = "1+y+x+x*y+x^2+y^2+x^2*y+x*y^2+x^3+y^3\n";
    int len;
    void* poly = getpol(ps, &len);
    
    printf("Monomial table after parsing:\n");
    printf("f4_msize = %d\n", f4_msize);
    printf("f4_nvars = %d\n\n", f4_nvars);
    
    // Print first few entries in monomial table
    printf("Monomial table entries (index: x_exp y_exp):\n");
    for (int i = 0; i < 20 && i < f4_msize; i++) {
        int* exp = f4_monom + i * (f4_nvars + 1);
        printf("  [%2d]: x^%d y^%d", i, exp[0], exp[1]);
        if (exp[0] == 0 && exp[1] == 0) printf(" (constant)");
        printf("\n");
    }
    
    f4mod_free();
    return 0;
}
