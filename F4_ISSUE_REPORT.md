# F4 Library Issue: Inconsistent Gröbner Basis Computation

## Summary
The F4 library appears to compute different Gröbner bases for the same polynomial system when using different prime moduli. This is mathematically incorrect - the structure of a Gröbner basis should be invariant across different prime fields (only the coefficients should differ by modular reduction).

## Minimal Test Case

### System Being Tested
```
Polynomials: x^2 - 1 = 0, y - x = 0
Variables: x, y
```

### Expected Mathematical Result
The Gröbner basis should contain:
- x + y (or some ordering)
- y^2 - 1 (since y = x and x^2 = 1 implies y^2 = 1)

In modular arithmetic with prime p, "y^2 - 1" is represented as "y^2 + (p-1)".

### Actual F4 Output

Using our wrapper code, F4 outputs:

1. **Prime 1073741827**: 
   ```
   +1*x+1*y
   +1*y^2+1073741826
   ```
   This is CORRECT: 1073741826 = 1073741827 - 1, representing y^2 - 1

2. **Prime 1048573** (and most other primes):
   ```
   +1*x+1*y
   +1*y^2+1
   ```
   This is WRONG: This represents y^2 + 1, which is mathematically incorrect

## Minimal C Test Program

```c
#include <stdio.h>
#include "axf4_wrapper.h"

void test_prime(unsigned int prime) {
    const char* vars[] = {"x", "y"};
    axf4_session_t session = axf4_create_session(prime, vars, 2);
    
    // Add x^2 - 1
    axf4_add_polynomial(session, "1*x^2-1");
    
    // Add y - x  
    axf4_add_polynomial(session, "1*y-1*x");
    
    // Compute GB
    axf4_result_t result = axf4_compute_groebner_basis(session);
    
    printf("Prime %u: %s\n", prime, result.groebner_basis);
    
    axf4_free_result(&result);
    axf4_destroy_session(session);
}

int main() {
    test_prime(1073741827);  // Shows correct: y^2+1073741826
    test_prime(1048573);     // Shows wrong: y^2+1
    return 0;
}
```

## Key Observations

1. Only ONE prime (1073741827) produces correct results in our testing
2. The F4 output shows different "pairs limit" values for different primes
3. We verified that F4's internal coefficient storage actually contains the value 1, not prime-1

## Questions for F4 Developers

1. **Is there a specific initialization or configuration we're missing?**
2. **Are there restrictions on which primes can be used with F4?**
3. **Is the polynomial input format "1*x^2-1" correct for representing x^2 - 1?**
4. **Does F4 have different internal algorithms based on prime size/form?**

## Our Implementation Details

We're using F4 through a C wrapper (axf4_wrapper.c) that:
- Calls f4mod_init() with the prime
- Parses polynomials with getpol() 
- Adds them with f4_addrow()
- Computes GB with f4gb_mod()
- Outputs results with putpol()

## Hypothesis

Given that:
- F4 is a well-tested library used in major CAS systems
- The bug affects nearly ALL primes except one specific value
- The polynomial parsing/output seems correct

We suspect we may be:
1. Using F4 incorrectly (missing initialization, wrong API usage)
2. Misunderstanding F4's polynomial representation
3. Having a bug in our wrapper code

## Next Steps

1. Compare our F4 usage with reference implementations
2. Check if we need to normalize/reduce input polynomials first
3. Verify our understanding of F4's internal representation

## Request for Help

Can you help us understand:
- The correct way to use F4 for modular Gröbner basis computation
- Any restrictions or requirements on prime selection
- Whether our test case reveals a real bug or user error