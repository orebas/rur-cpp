#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// This test specifically triggers the heap buffer overflow

#define INT64 long long int
typedef INT64 INT;

void test_specific_overflow_pattern(int n) {
    unsigned char *buf;
    INT *vec;
    int i, j, k, l;
    int non_zeros = 0;
    
    // Create vector
    vec = calloc(n + 1, sizeof(INT));
    
    // Create a pattern that maximizes encoding size:
    // We need enough non-zeros with large gaps to exhaust the buffer
    // The key is having n+1 potential positions but allocation for only n
    
    // Strategy: Create non-zeros at positions with gaps > 255
    // This forces each to use 1 + sizeof(INT) bytes
    
    // Place non-zeros strategically
    if (n >= 1000) {
        // For large n, place non-zeros at specific intervals
        for (i = n; i >= 0; i -= 300) {
            vec[i] = i + 1;
            non_zeros++;
        }
    } else {
        // For smaller n, we need a different pattern
        // Place at edges and middle with large gaps
        vec[n] = n + 1;
        vec[0] = 1;
        non_zeros = 2;
        
        if (n >= 512) {
            vec[256] = 257;
            non_zeros++;
        }
    }
    
    printf("Testing with n=%d, %d non-zero elements\n", n, non_zeros);
    
    // Calculate actual required size
    size_t actual_required = 0;
    l = -1;
    for (i = n; i >= 0; i--) {
        if (vec[i] == 0) continue;
        
        if (l == -1) {
            actual_required += sizeof(INT);  // First element
        } else if (l - i <= 255) {
            actual_required += 1;  // Small gap
        } else {
            actual_required += 1 + sizeof(INT);  // Large gap
        }
        l = i;
    }
    actual_required += 1;  // Final null byte
    
    // Original buggy allocation
    size_t buggy_size = n * sizeof(INT) + n;
    
    printf("Buggy allocation: %zu bytes\n", buggy_size);
    printf("Actually required: %zu bytes\n", actual_required);
    
    if (actual_required > buggy_size) {
        printf("*** OVERFLOW WILL OCCUR! ***\n");
        printf("Overflow by: %zu bytes\n", actual_required - buggy_size);
    } else {
        printf("No overflow with this pattern (need %zu more bytes to trigger)\n", 
               buggy_size - actual_required + 1);
    }
    
    // Now test the fixed allocation
    size_t fixed_size = (n + 1) * sizeof(INT) + (n + 1);
    printf("\nFixed allocation: %zu bytes\n", fixed_size);
    printf("Safety margin: %zu bytes\n", fixed_size - actual_required);
    
    free(vec);
}

// Test with exact pattern that should maximize encoding
void test_worst_case_encoding() {
    printf("\n=== WORST CASE ENCODING TEST ===\n");
    
    // For the worst case, we need all n+1 elements to be non-zero
    // AND have them arranged so gaps are > 255
    // This is only possible if n is small enough
    
    int n = 3;  // Small n to demonstrate
    unsigned char *buf;
    INT *vec = calloc(n + 1, sizeof(INT));
    
    // Set specific pattern - this might match Cyclic-3
    vec[3] = 1;  // Highest index
    vec[2] = 2;  
    vec[1] = 3;
    vec[0] = 4;  // All 4 elements non-zero
    
    printf("Testing n=%d with all %d elements non-zero\n", n, n+1);
    
    // Calculate encoding size
    size_t encoding_size = sizeof(INT);  // First element at index 3
    encoding_size += 1;  // Gap from 3 to 2 is 1
    encoding_size += 1;  // Gap from 2 to 1 is 1  
    encoding_size += 1;  // Gap from 1 to 0 is 1
    encoding_size += 1;  // Final null
    
    size_t buggy_alloc = n * sizeof(INT) + n;  // 3*8 + 3 = 27
    size_t fixed_alloc = (n+1) * sizeof(INT) + (n+1);  // 4*8 + 4 = 36
    
    printf("Encoding requires: %zu bytes\n", encoding_size);
    printf("Buggy allocation: %zu bytes\n", buggy_alloc);
    printf("Fixed allocation: %zu bytes\n", fixed_alloc);
    
    if (encoding_size > buggy_alloc) {
        printf("*** This triggers overflow! ***\n");
    }
    
    // Actually try the encoding with bounds checking
    buf = malloc(buggy_alloc + 100);  // Extra space to avoid actual crash
    memset(buf, 0xFF, buggy_alloc + 100);  // Fill with sentinel
    
    int k = 0, l = -1;
    for (int i = n; i >= 0; i--) {
        if (vec[i] == 0) continue;
        
        if (l == -1) {
            if (k + sizeof(INT) > buggy_alloc) {
                printf("Would overflow at first element (k=%d)\n", k);
                break;
            }
            *(INT *)(buf + k) = i;
            k += sizeof(INT);
        } else if (l - i <= 255) {
            if (k + 1 > buggy_alloc) {
                printf("Would overflow at gap write (k=%d, i=%d)\n", k, i);
                break;
            }
            buf[k++] = (unsigned char)(l - i);
        } else {
            if (k + 1 + sizeof(INT) > buggy_alloc) {
                printf("Would overflow at large gap (k=%d, i=%d)\n", k, i);
                break;
            }
            buf[k++] = 0;
            *(INT *)(buf + k) = i;
            k += sizeof(INT);
        }
        l = i;
    }
    
    printf("After encoding: k=%d\n", k);
    
    // The critical line that causes overflow
    if (k >= buggy_alloc) {
        printf("*** OVERFLOW at final null byte! k=%d >= allocated=%zu ***\n", 
               k, buggy_alloc);
    } else {
        buf[k++] = 0;
        printf("No overflow. Final k=%d\n", k);
    }
    
    // Check sentinel bytes
    int overflow_detected = 0;
    for (size_t i = buggy_alloc; i < buggy_alloc + 10; i++) {
        if (buf[i] != 0xFF) {
            printf("*** MEMORY CORRUPTION DETECTED at offset %zu ***\n", 
                   i - buggy_alloc);
            overflow_detected = 1;
            break;
        }
    }
    
    if (!overflow_detected && k > buggy_alloc) {
        printf("*** LOGICAL OVERFLOW: Wrote %d bytes into %zu byte buffer ***\n",
               k, buggy_alloc);
    }
    
    free(buf);
    free(vec);
}

int main() {
    printf("Testing heap buffer overflow conditions in f4_reduce_export\n");
    printf("sizeof(INT) = %zu\n\n", sizeof(INT));
    
    // Test various sizes
    test_specific_overflow_pattern(3);     // Like Cyclic-3
    printf("\n");
    test_specific_overflow_pattern(100);
    printf("\n");
    test_specific_overflow_pattern(1000);
    
    // Test worst case
    test_worst_case_encoding();
    
    return 0;
}