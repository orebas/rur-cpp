#include <stdio.h>
#include <stdlib.h>

// F4 declarations
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
    long long len;
    long long sug;
    long long fac;
    long long *cof;
    long long *mon;
    unsigned char *ind;
    long long siz;
} f4row;

void decode_monomial(unsigned int mon, int* x_exp, int* y_exp) {
    *x_exp = mon & 0x7FFF;
    *y_exp = (mon >> 15) & 0x7FFF;
}

int main() {
    printf("Testing F4 term ordering to verify leading term position\n");
    printf("=========================================================\n\n");
    
    vars[0] = "x";
    vars[1] = "y";
    
    f4mod_init(2, 0, 1073741827);
    
    // Test 1: x^2 - should have leading term x^2
    printf("Test 1: x^2\n");
    const char* poly1 = "x^2\n";
    int len;
    void* p1 = getpol(poly1, &len);
    if (p1) {
        f4row* row = (f4row*)p1;
        printf("  Number of terms: %lld\n", row->len);
        for (int i = 0; i < row->len; i++) {
            int x_exp, y_exp;
            decode_monomial(row->mon[i], &x_exp, &y_exp);
            printf("  Term[%d]: coeff=%lld, x^%d*y^%d (monomial=%u)\n", 
                   i, row->cof[i], x_exp, y_exp, (unsigned int)row->mon[i]);
        }
        f4_addrow(p1);
    }
    
    // Test 2: 1 + x + x^2 - should have terms in ascending order
    printf("\nTest 2: 1 + x + x^2\n");
    const char* poly2 = "1+x+x^2\n";
    void* p2 = getpol(poly2, &len);
    if (p2) {
        f4row* row = (f4row*)p2;
        printf("  Number of terms: %lld\n", row->len);
        for (int i = 0; i < row->len; i++) {
            int x_exp, y_exp;
            decode_monomial(row->mon[i], &x_exp, &y_exp);
            printf("  Term[%d]: coeff=%lld, x^%d*y^%d (monomial=%u)\n", 
                   i, row->cof[i], x_exp, y_exp, (unsigned int)row->mon[i]);
        }
    }
    
    // Compute GB and check leading terms
    f4gb_mod();
    
    printf("\nGröbner basis (%d polynomials):\n", f4_aload);
    for (int i = 0; i < f4_aload && i < 3; i++) {
        f4row* gb_row = (f4row*)f4_array[i];
        printf("\nPolynomial %d:\n", i);
        printf("  Full: ");
        putpol(f4_array[i], stdout);
        printf("\n");
        
        if (gb_row->len > 0) {
            printf("  Number of terms: %lld\n", gb_row->len);
            
            // Show first term
            int x0, y0;
            decode_monomial(gb_row->mon[0], &x0, &y0);
            printf("  FIRST term [0]: coeff=%lld, x^%d*y^%d\n", 
                   gb_row->cof[0], x0, y0);
            
            // Show last term
            int xL, yL;
            decode_monomial(gb_row->mon[gb_row->len-1], &xL, &yL);
            printf("  LAST term [%lld]: coeff=%lld, x^%d*y^%d\n", 
                   gb_row->len-1, gb_row->cof[gb_row->len-1], xL, yL);
            
            // Determine which has higher degree
            int deg_first = x0 + y0;
            int deg_last = xL + yL;
            printf("  First term total degree: %d\n", deg_first);
            printf("  Last term total degree: %d\n", deg_last);
            printf("  => Leading term is at position: %s\n", 
                   deg_first >= deg_last ? "FIRST [0]" : "LAST");
        }
    }
    
    f4mod_free();
    return 0;
}
