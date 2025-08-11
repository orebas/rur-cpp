#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Test case that triggers the actual overflow

#define INT64 long long int
typedef INT64 INT;

void test_dense_overflow(int n) {
    unsigned char *buf;
    INT *vec;
    int i, j, k, l;
    
    // Create a dense vector with all non-zero values
    vec = calloc(n + 1, sizeof(INT));
    
    // Fill with non-zero values - this maximizes memory usage
    for (i = 0; i <= n; i++) {
        vec[i] = i + 1;  // All non-zero
    }
    
    printf("Testing DENSE vector with n=%d\n", n);
    printf("Processing %d elements (0 to %d), ALL NON-ZERO\n", n+1, n);
    
    // Original buggy allocation
    size_t buggy_size = n*sizeof(INT)+n;
    printf("Buggy allocation: %zu bytes\n", buggy_size);
    buf = malloc(buggy_size);
    
    // Encoding loop
    j = k = 0; 
    l = -1;
    
    for (i = n; i >= 0; i--) {
        if (vec[i] == 0) continue;
        
        if (l == -1) {
            // First non-zero element
            printf("  First element at i=%d, writing %zu bytes at k=%d\n", i, sizeof(INT), k);
            if (k + sizeof(INT) > buggy_size) {
                printf("  ERROR: Overflow at first element!\n");
                goto cleanup1;
            }
            *(INT *)(buf+k) = i;
            k += sizeof(INT);
        }
        else if (l - i <= 255) {
            // Small gap - single byte
            if (k >= buggy_size) {
                printf("  ERROR: Overflow at i=%d, k=%d (gap=%d)\n", i, k, l-i);
                goto cleanup1;
            }
            buf[k++] = (unsigned char)(l-i);
        }
        else {
            // Large gap - should not happen with dense data
            printf("  Large gap at i=%d\n", i);
            if (k + 1 + sizeof(INT) > buggy_size) {
                printf("  ERROR: Overflow at large gap!\n");
                goto cleanup1;
            }
            buf[k++] = 0;
            *(INT *)(buf+k) = i;
            k += sizeof(INT);
        }
        l = i;
        j++;
    }
    
    printf("After encoding %d non-zeros: k=%d, allocated=%zu\n", j, k, buggy_size);
    printf("Remaining space: %zu bytes\n", buggy_size - k);
    
    // The problematic line 1982 - extra null byte
    if (k >= buggy_size) {
        printf("*** OVERFLOW DETECTED! ***\n");
        printf("  Trying to write at position k=%d\n", k);
        printf("  Allocated size: %zu\n", buggy_size);
        printf("  This is the heap buffer overflow!\n");
    } else {
        buf[k++] = 0;
        printf("No overflow with this test case. Final k=%d\n", k);
    }

cleanup1:
    free(buf);
    
    // Now test with the fix
    size_t fixed_size = (n+1)*sizeof(INT)+(n+1);
    printf("\nFixed allocation: %zu bytes (increase of %zu bytes)\n", 
           fixed_size, fixed_size - buggy_size);
    buf = malloc(fixed_size);
    
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
    
    printf("After encoding with fix: k=%d, allocated=%zu\n", k, fixed_size);
    printf("Remaining space: %zu bytes\n", fixed_size - k);
    
    if (k >= fixed_size) {
        printf("ERROR: Still overflows with fix!\n");
    } else {
        buf[k++] = 0;
        printf("Success with fix: No overflow. Final k=%d\n", k);
        printf("Final remaining space: %zu bytes\n", fixed_size - k);
    }
    
    free(buf);
    free(vec);
}

int main() {
    printf("Testing heap buffer overflow in f4_reduce_export\n");
    printf("sizeof(INT) = %zu\n\n", sizeof(INT));
    
    // Test with dense data that should trigger overflow
    test_dense_overflow(100);
    printf("\n============================================================\n\n");
    test_dense_overflow(1000);
    
    return 0;
}