#include <stdio.h>

// Test how F4 encodes monomials
// Based on looking at the axf4 code, it seems to pack exponents

void test_encoding(unsigned int mon) {
    // Try different interpretations
    printf("Monomial value: %u (0x%x)\n", mon, mon);
    
    // Interpretation 1: 15 bits per variable
    int x_exp = mon & 0x7FFF;
    int y_exp = (mon >> 15) & 0x7FFF;
    printf("  15-bit encoding: x^%d * y^%d\n", x_exp, y_exp);
    
    // Interpretation 2: Maybe it's packed differently?
    // Check if it's actually the exponent vector index
    printf("  Direct value: %u\n", mon);
}

int main() {
    printf("Testing monomial encodings:\n");
    printf("===========================\n\n");
    
    // From our test: x^2 gave monomial=1
    printf("x^2 gave monomial=1:\n");
    test_encoding(1);
    
    printf("\nx gave monomial=3:\n");
    test_encoding(3);
    
    printf("\n1 (constant) gave monomial=2:\n");
    test_encoding(2);
    
    // These don't look like exponent encodings - they look like indices\!
    printf("\n\nThese look like monomial BASIS INDICES, not exponent vectors\!\n");
    
    return 0;
}
