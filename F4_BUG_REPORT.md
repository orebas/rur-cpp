# Bug Report: Segmentation Fault in F4 Gröbner Basis Library

## Summary
The F4 library (`axf4_lib.c`) crashes with a segmentation fault when computing Gröbner bases for certain polynomial/prime combinations. The crash occurs in the `f4_jordan_thread` function due to a NULL pointer dereference.

## Environment
- Platform: Linux x86_64
- Compiler: GCC
- F4 Library Version: As included in RUR-CPP project

## Bug Details

### Description
The F4 library crashes during the Jordan elimination phase when `f4_reduce_export` returns NULL (indicating a zero polynomial after reduction), but the code doesn't check for this condition before dereferencing the result.

### Reproduction Steps
1. Initialize F4 with prime = 1073741827 (or 65537)
2. Add polynomial "x - 3" 
3. Call `f4gb_mod()` to compute Gröbner basis
4. Program crashes with segmentation fault

### Test Case
```c
#include <stdio.h>
extern void f4mod_init(int nvars, int elim, long long prime);
extern void f4gb_mod(void);
extern void f4_addrow(void* row);
extern void* getpol(const char* str, int* length);
extern char* vars[256];

int main() {
    vars[0] = "x";
    f4mod_init(1, 0, 1073741827);
    
    int len;
    void* poly = getpol("x - 3\n", &len);
    f4_addrow(poly);
    
    f4gb_mod(); // CRASHES HERE
    return 0;
}
```

### Stack Trace
```
Thread received signal SIGSEGV, Segmentation fault.
0x000055555555b958 in f4_jordan_thread (par=0x7fffffffd150) at src/axf4_lib.c:2471
2471    x = *(INT *)(a->ind);

#0  0x000055555555b958 in f4_jordan_thread at src/axf4_lib.c:2471
#1  0x00007ffff7c9caa4 in start_thread
#2  0x00007ffff7d29c3c in clone3

Register rdi = 0xffffffffffffffff (-1)
```

### Root Cause Analysis
In `f4_jordan_thread` (line 2470), the code calls:
```c
a = f4_reduce_export(s,vec);
x = *(INT *)(a->ind);  // CRASH: a is NULL
```

The function `f4_reduce_export` can return NULL when the reduced vector has no non-zero elements. The code doesn't check for this condition before dereferencing `a->ind`.

### Why Specific Primes?
The bug manifests with certain primes (1073741827, 65537) but not others (131063) due to the arithmetic properties of the reduction modulo different primes. With certain primes, the reduction of "x - 3" results in a zero polynomial, triggering the NULL return.

## Fix

### Patch
The bug affects both the multi-threaded and single-threaded implementations of the Jordan elimination:

```diff
--- src/axf4_lib.c.orig
+++ src/axf4_lib.c
@@ -2468,6 +2468,11 @@ void f4_jordan_thread(void *par)
 	s = f4_reduce_import(vec, f4_array[st], 1);
 	f4_reduce_vector(s-1,vec,piv);
 	a = f4_reduce_export(s,vec);
+	if (!a) {
+		/* Handle case where reduction results in zero polynomial */
+		goto block;
+	}
+	
 	x = *(INT *)(a->ind);
 	red[x] = a;
 
@@ -2555,7 +2560,11 @@ void f4_jordan()
 	piv[i] = 0;
 	f4_reduce_vector(i,vec,piv);
 	a = f4_reduce_export(i,vec);
-	piv[i] = f4_reduce_monic(a);
+	if (a) {
+		piv[i] = f4_reduce_monic(a);
+	} else {
+		piv[i] = 0;
+	}
```

### Verification
After applying this fix, the test case runs without crashing for all tested prime values.

## Additional Notes

1. The same NULL check pattern is already implemented correctly in other parts of the code (e.g., `f4_reduce_thread` at line 2255).

2. The crash occurs even with `num_thread = 1` because the multi-threaded code path is still used (it just spawns 1 thread).

3. This is a critical bug that can cause unpredictable crashes depending on the input polynomials and prime modulus.

## Recommendation
Apply the provided patch to add proper NULL checking in `f4_jordan_thread`. This is a minimal, safe fix that follows the existing error handling pattern in the codebase.