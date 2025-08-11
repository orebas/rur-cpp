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

typedef struct f4row {
    long long len;
    long long sug;
    long long fac;
    long long *cof;
    long long *mon;
    unsigned char *ind;
    long long siz;
} f4row;

int main() {
    printf("Testing formatted polynomial with F4\n\n");
    
    vars[0] = "x";
    vars[1] = "y";
    
    int prime = 1073741827;
    f4mod_init(2, 0, prime);
    
    // Test our formatter output
    const char* poly_str = "-1+1*x^2+1*y^2\n";
    printf("Parsing: %s", poly_str);
    
    int len;
    void* poly = getpol(poly_str, &len);
    
    if (poly == NULL || len == 0) {
        printf("FAILED to parse\n");
        f4mod_free();
        return 1;
    }
    
    printf("SUCCESS: parsed with %d terms\n", len);
    
    f4row* row = (f4row*)poly;
    printf("Details: %lld terms\n", row->len);
    for (int i = 0; i < row->len && i < 5; i++) {
        int x_exp = row->mon[i] & 0x7FFF;
        int y_exp = (row->mon[i] >> 15) & 0x7FFF;
        printf("  Term %d: %lld * x^%d * y^%d\n", i, row->cof[i], x_exp, y_exp);
    }
    
    // Add and compute GB
    f4_addrow(poly);
    f4gb_mod();
    
    printf("\nGB has %d polynomial(s)\n", f4_aload);
    for (int i = 0; i < f4_aload && i < 3; i++) {
        f4row* gb_row = (f4row*)f4_array[i];
        printf("Poly %d: ", i);
        putpol(f4_array[i], stdout);
        printf("\n");
        
        if (gb_row->len > 0) {
            int x_exp = gb_row->mon[0] & 0x7FFF;
            int y_exp = (gb_row->mon[0] >> 15) & 0x7FFF;
            printf("  Leading term: %lld * x^%d * y^%d\n", gb_row->cof[0], x_exp, y_exp);
        }
    }
    
    f4mod_free();
    return 0;
}