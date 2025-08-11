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
extern int f4_nvars;
extern int* f4_monom;

typedef struct f4row {
    long long len;
    long long sug;
    long long fac;
    long long *cof;
    long long *mon;
    unsigned char *ind;
    long long siz;
} f4row;

#define F4_EXPON(m) (f4_monom + (m)*(f4_nvars+1))

void print_monomial_details(unsigned int mon_idx) {
    int* exps = F4_EXPON(mon_idx);
    int total_deg = 0;
    printf("    Monomial index %u: ", mon_idx);
    for (int i = 0; i < f4_nvars; i++) {
        if (exps[i] > 0) {
            printf("%c^%d ", 'x' + i, exps[i]);
            total_deg += exps[i];
        }
    }
    if (total_deg == 0) printf("1");
    printf(" (total degree = %d)\n", total_deg);
}

int main() {
    printf("Testing F4 polynomial term ordering\n");
    printf("====================================\n\n");
    
    vars[0] = "x";
    vars[1] = "y";
    
    f4mod_init(2, 0, 1073741827);
    
    // Add a polynomial with known structure: y + x^2 + x*y^2 + x^3*y
    // This has degrees: 1, 2, 3, 4
    const char* poly_str = "y+x^2+x*y^2+x^3*y\n";
    printf("Input polynomial: %s", poly_str);
    
    int len;
    void* poly = getpol(poly_str, &len);
    if (\!poly) {
        printf("Failed to parse polynomial\n");
        f4mod_free();
        return 1;
    }
    
    f4row* row = (f4row*)poly;
    printf("Parsed with %lld terms\n\n", row->len);
    
    printf("Term ordering in f4row structure:\n");
    for (int i = 0; i < row->len; i++) {
        printf("  Term[%d]: coeff=%lld, monomial_index=%lld\n", 
               i, row->cof[i], row->mon[i]);
        print_monomial_details(row->mon[i]);
    }
    
    // Add polynomial and compute GB
    f4_addrow(poly);
    f4gb_mod();
    
    printf("\n\nGroebner basis has %d polynomial(s)\n", f4_aload);
    for (int i = 0; i < f4_aload && i < 3; i++) {
        f4row* gb_row = (f4row*)f4_array[i];
        printf("\nGB Polynomial %d (%lld terms):\n", i, gb_row->len);
        printf("  String form: ");
        putpol(f4_array[i], stdout);
        printf("\n");
        
        if (gb_row->len > 0) {
            printf("  FIRST term (index 0):\n");
            printf("    coeff=%lld, monomial_index=%lld\n", 
                   gb_row->cof[0], gb_row->mon[0]);
            print_monomial_details(gb_row->mon[0]);
            
            if (gb_row->len > 1) {
                printf("  LAST term (index %lld):\n", gb_row->len - 1);
                printf("    coeff=%lld, monomial_index=%lld\n", 
                       gb_row->cof[gb_row->len-1], gb_row->mon[gb_row->len-1]);
                print_monomial_details(gb_row->mon[gb_row->len-1]);
            }
        }
    }
    
    printf("\n\nCONCLUSION:\n");
    printf("The leading term (highest degree) is at position: ");
    if (row->len > 1) {
        int first_idx = row->mon[0];
        int last_idx = row->mon[row->len-1];
        int* first_exp = F4_EXPON(first_idx);
        int* last_exp = F4_EXPON(last_idx);
        int first_deg = 0, last_deg = 0;
        for (int i = 0; i < f4_nvars; i++) {
            first_deg += first_exp[i];
            last_deg += last_exp[i];
        }
        printf("%s\n", last_deg > first_deg ? "LAST" : "FIRST");
    }
    
    f4mod_free();
    return 0;
}
