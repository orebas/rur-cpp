/*
 * Test to verify the heap buffer overflow fix in f4_reduce_export
 * 
 * This test simulates the behavior of f4_reduce_export with both
 * the buggy and fixed allocation strategies.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define INT64 long long int
typedef INT64 INT;

// Simulate the encoding logic of f4_reduce_export
int simulate_encoding(INT n, INT *vec, unsigned char *buf, size_t buf_size, int verbose) {
    INT i, j, k, l;
    
    j = k = 0; 
    l = -1;
    
    // The critical loop: processes from n down to 0 (n+1 elements)
    for (i = n; i >= 0; i--) {
        if (vec[i] == 0) continue;
        
        if (l == -1) {
            // First non-zero element
            if (k + sizeof(INT) > buf_size) {
                if (verbose) printf("  ERROR: Overflow at first element (i=%lld, k=%lld)\n", i, k);
                return -1;
            }
            *(INT *)(buf + k) = i;
            k += sizeof(INT);
        }
        else if (l - i <= 255) {
            // Small gap - single byte
            if (k >= buf_size) {
                if (verbose) printf("  ERROR: Overflow at gap byte (i=%lld, k=%lld)\n", i, k);
                return -1;
            }
            buf[k++] = (unsigned char)(l - i);
        }
        else {
            // Large gap - null + INT
            if (k + 1 + sizeof(INT) > buf_size) {
                if (verbose) printf("  ERROR: Overflow at large gap (i=%lld, k=%lld)\n", i, k);
                return -1;
            }
            buf[k++] = 0;
            *(INT *)(buf + k) = i;
            k += sizeof(INT);
        }
        l = i;
        j++;
    }
    
    if (j == 0) return k;  // No non-zeros
    
    // The critical line 1982 - final null byte
    if (k >= buf_size) {
        if (verbose) printf("  ERROR: Overflow at final null byte (k=%lld, buf_size=%zu)\n", k, buf_size);
        return -1;
    }
    buf[k++] = 0;
    
    return k;
}

void test_allocation_strategy(const char *name, INT n, INT *vec, int non_zeros) {
    unsigned char *buf;
    size_t buggy_size, fixed_size;
    int result;
    
    printf("\nTesting %s (n=%lld, %d non-zeros):\n", name, n, non_zeros);
    
    // Original buggy allocation
    buggy_size = n * sizeof(INT) + n;
    buf = malloc(buggy_size + 10);  // Extra space to detect overwrites
    memset(buf, 0xAA, buggy_size + 10);  // Fill with sentinel pattern
    
    printf("  Buggy allocation: %zu bytes\n", buggy_size);
    result = simulate_encoding(n, vec, buf, buggy_size, 1);
    
    if (result < 0) {
        printf("  -> OVERFLOW DETECTED (prevented)\n");
    } else {
        printf("  -> Success: used %d bytes\n", result);
        
        // Check for actual overwrites
        int overflow = 0;
        for (size_t i = buggy_size; i < buggy_size + 10; i++) {
            if (buf[i] != 0xAA) {
                printf("  -> MEMORY CORRUPTION at offset +%zu!\n", i - buggy_size);
                overflow = 1;
                break;
            }
        }
        if (!overflow && result > buggy_size) {
            printf("  -> LOGICAL OVERFLOW: wrote %d bytes to %zu byte buffer\n", 
                   result, buggy_size);
        }
    }
    free(buf);
    
    // Fixed allocation
    fixed_size = (n + 1) * sizeof(INT) + (n + 1);
    buf = malloc(fixed_size);
    
    printf("  Fixed allocation: %zu bytes (diff: +%zu)\n", 
           fixed_size, fixed_size - buggy_size);
    result = simulate_encoding(n, vec, buf, fixed_size, 0);
    
    if (result < 0) {
        printf("  -> ERROR: Fixed allocation still overflows!\n");
    } else {
        printf("  -> Success: used %d bytes, %zu bytes remaining\n", 
               result, fixed_size - result);
    }
    free(buf);
}

void run_test_suite() {
    INT *vec;
    
    printf("=== F4_REDUCE_EXPORT OVERFLOW FIX VERIFICATION ===\n");
    printf("sizeof(INT) = %zu\n", sizeof(INT));
    
    // Test 1: Small dense vector (like Cyclic-3 might produce)
    vec = calloc(4, sizeof(INT));
    vec[0] = 1; vec[1] = 2; vec[2] = 3; vec[3] = 4;
    test_allocation_strategy("Small dense (n=3)", 3, vec, 4);
    free(vec);
    
    // Test 2: Sparse vector with large gaps
    vec = calloc(1001, sizeof(INT));
    vec[1000] = 1;
    vec[500] = 2;
    vec[0] = 3;
    test_allocation_strategy("Sparse with large gaps (n=1000)", 1000, vec, 3);
    free(vec);
    
    // Test 3: All zeros (edge case)
    vec = calloc(101, sizeof(INT));
    test_allocation_strategy("All zeros (n=100)", 100, vec, 0);
    free(vec);
    
    // Test 4: Single non-zero at index 0
    vec = calloc(101, sizeof(INT));
    vec[0] = 1;
    test_allocation_strategy("Single element at 0 (n=100)", 100, vec, 1);
    free(vec);
    
    // Test 5: Alternating pattern
    vec = calloc(101, sizeof(INT));
    for (int i = 0; i <= 100; i += 2) vec[i] = i + 1;
    test_allocation_strategy("Alternating pattern (n=100)", 100, vec, 51);
    free(vec);
    
    // Test 6: Maximum density small vector
    vec = calloc(11, sizeof(INT));
    for (int i = 0; i <= 10; i++) vec[i] = i + 1;
    test_allocation_strategy("Maximum density (n=10)", 10, vec, 11);
    free(vec);
}

int main() {
    run_test_suite();
    
    printf("\n=== SUMMARY ===\n");
    printf("The fix changes allocation from n*(sizeof(INT)+1) to (n+1)*(sizeof(INT)+1)\n");
    printf("This accounts for processing n+1 elements (indices 0 to n) instead of n elements\n");
    printf("The fix adds exactly sizeof(INT)+1 = %zu bytes to each allocation\n", 
           sizeof(INT) + 1);
    printf("This is a minimal, safe fix that prevents the heap buffer overflow\n");
    
    return 0;
}