#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Minimal reproduction of the overflow issue
// This simulates the problematic part of f4_reduce_export

#define INT64 long long int
typedef INT64 INT;

void test_overflow_reproduction(int n) {
    unsigned char *buf;
    INT *vec;
    int i, j, k, l;
    
    // Create a worst-case sparse vector
    vec = calloc(n + 1, sizeof(INT));
    
    // Set up sparse non-zero values with large gaps
    // This triggers the worst-case memory usage
    for (i = 0; i <= n; i += 300) {
        vec[i] = i + 1;  // Non-zero value
    }
    
    printf("Testing with n=%d\n", n);
    printf("Processing %d elements (0 to %d)\n", n+1, n);
    
    // Original buggy allocation (as in line 1961)
    printf("Buggy allocation: %zu bytes\n", n*sizeof(INT)+n);
    buf = malloc(n*sizeof(INT)+n);
    
    // Encoding loop (simplified from lines 1963-1979)
    j = k = 0; 
    l = -1;
    
    for (i = n; i >= 0; i--) {
        if (vec[i] == 0) continue;
        
        if (l == -1) {
            // First non-zero element
            if (k + sizeof(INT) > n*sizeof(INT)+n) {
                printf("ERROR: Would overflow at first element write (k=%d)\n", k);
                break;
            }
            *(INT *)(buf+k) = i;
            k += sizeof(INT);
        }
        else if (l - i <= 255) {
            // Small gap - single byte
            if (k + 1 > n*sizeof(INT)+n) {
                printf("ERROR: Would overflow at small gap write (k=%d)\n", k);
                break;
            }
            buf[k++] = (unsigned char)(l-i);
        }
        else {
            // Large gap - null byte + INT
            if (k + 1 + sizeof(INT) > n*sizeof(INT)+n) {
                printf("ERROR: Would overflow at large gap write (k=%d)\n", k);
                break;
            }
            buf[k++] = 0;
            *(INT *)(buf+k) = i;
            k += sizeof(INT);
        }
        l = i;
        j++;
    }
    
    printf("After encoding loop: k=%d, allocated=%zu\n", k, n*sizeof(INT)+n);
    
    // The problematic line 1982 - extra null byte
    if (k + 1 > n*sizeof(INT)+n) {
        printf("ERROR: OVERFLOW at line 1982! k=%d, trying to write at position %d\n", k, k);
        printf("Allocated size: %zu, Required size: %d\n", n*sizeof(INT)+n, k+1);
    } else {
        buf[k++] = 0;  // This would cause the overflow
        printf("Success: No overflow. Final k=%d\n", k);
    }
    
    free(buf);
    
    // Now test with the fix
    printf("\nFixed allocation: %zu bytes\n", (n+1)*sizeof(INT)+(n+1));
    buf = malloc((n+1)*sizeof(INT)+(n+1));
    
    // Re-run the same encoding
    j = k = 0; 
    l = -1;
    
    for (i = n; i >= 0; i--) {
        if (vec[i] == 0) continue;
        
        if (l == -1) {
            *(INT *)(buf+k) = i;
            k += sizeof(INT);
        }
        else if (l - i <= 255) {
            buf[k++] = (unsigned char)(l-i);
        }
        else {
            buf[k++] = 0;
            *(INT *)(buf+k) = i;
            k += sizeof(INT);
        }
        l = i;
        j++;
    }
    
    printf("After encoding loop with fix: k=%d, allocated=%zu\n", k, (n+1)*sizeof(INT)+(n+1));
    
    // The line that previously caused overflow
    if (k + 1 > (n+1)*sizeof(INT)+(n+1)) {
        printf("ERROR: Still overflows with fix!\n");
    } else {
        buf[k++] = 0;
        printf("Success with fix: No overflow. Final k=%d\n", k);
    }
    
    free(buf);
    free(vec);
}

int main() {
    printf("Testing heap buffer overflow in f4_reduce_export\n");
    printf("sizeof(INT) = %zu\n\n", sizeof(INT));
    
    // Test with various sizes that trigger the overflow
    test_overflow_reproduction(1000);
    printf("\n");
    test_overflow_reproduction(500);
    
    return 0;
}