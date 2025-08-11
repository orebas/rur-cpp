/*
 * Minimal test case for F4 bug with x^2 + y^2 - 1
 * This uses the axf4 library directly to isolate whether the bug is in F4 or our wrapper
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Forward declarations for axf4 functions
// These are defined in axf4_lib.c

extern int axf4_main(int argc, char** argv);

int main() {
    printf("Testing x^2 + y^2 - 1 with F4 directly\n");
    printf("=====================================\n\n");
    
    // Test with different primes
    int primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 
                    101, 1009, 10007, 100003, 1000003,
                    793493497, 1073741827};
    int num_primes = sizeof(primes) / sizeof(primes[0]);
    
    for (int i = 0; i < num_primes; i++) {
        int prime = primes[i];
        printf("Testing with prime p = %d (p mod 4 = %d):\n", prime, prime % 4);
        
        // Create input file for F4
        char filename[100];
        sprintf(filename, "/tmp/f4_test_%d.txt", prime);
        FILE* f = fopen(filename, "w");
        if (!f) {
            printf("  ERROR: Could not create input file\n");
            continue;
        }
        
        // Write the polynomial system in F4 format
        fprintf(f, "# Testing x^2 + y^2 - 1 over F_%d\n", prime);
        fprintf(f, "p = %d\n", prime);
        fprintf(f, "vars = [x, y]\n");
        fprintf(f, "polys = [\n");
        fprintf(f, "  x^2 + y^2 - 1\n");
        fprintf(f, "]\n");
        fclose(f);
        
        // Now we need to call F4 directly
        // The F4 implementation expects command line arguments
        char* argv[] = {
            "axf4",
            filename,
            NULL
        };
        
        printf("  Input file created: %s\n", filename);
        printf("  Would call: axf4_main(2, argv) if it were exported\n");
        
        // Note: The actual F4 implementation needs to be called here
        // but we need to check how it's exported/linked
        
        printf("\n");
    }
    
    return 0;
}