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

typedef struct f4row {
    long long len, sug, fac;
    long long *cof, *mon;
    unsigned char *ind;
    long long siz;
} f4row;

#define F4_EXPON(m) (f4_monom + (m)*(f4_nvars+1))

void print_mon(unsigned int idx) {
    int* exps = F4_EXPON(idx);
    int deg = 0;
    printf("    idx %u: ", idx);
    for (int i = 0; i < f4_nvars; i++) {
        if (exps[i] > 0) {
            printf("%c^%d ", 120+i, exps[i]);
            deg += exps[i];
        }
    }
    if (deg == 0) printf("1");
    printf(" (deg=%d)\n", deg);
}

int main() {
    printf("F4 term ordering test\n\n");
    vars[0] = "x";
    vars[1] = "y";
    f4mod_init(2, 0, 1073741827);
    
    const char* ps = "y+x^2+x*y^2+x^3*y\n";
    printf("Input: %s", ps);
    
    int len;
    void* poly = getpol(ps, &len);
    if (!poly) return 1;
    
    f4row* row = (f4row*)poly;
    printf("Terms: %lld\n\n", row->len);
    
    for (int i = 0; i < row->len; i++) {
        printf("  [%d]: coeff=%lld, ", i, row->cof[i]);
        print_mon(row->mon[i]);
    }
    
    f4_addrow(poly);
    f4gb_mod();
    
    printf("\nGB (%d polys):\n", f4_aload);
    f4row* gb = (f4row*)f4_array[0];
    printf("  ");
    putpol(f4_array[0], stdout);
    printf("\n");
    
    if (gb->len > 0) {
        printf("  FIRST: ");
        print_mon(gb->mon[0]);
        if (gb->len > 1) {
            printf("  LAST:  ");
            print_mon(gb->mon[gb->len-1]);
        }
    }
    
    f4mod_free();
    return 0;
}