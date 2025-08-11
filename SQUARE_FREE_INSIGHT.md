# The Square-Free Part: Key to Understanding Katsura4Degenerate

## The Critical Discovery

Julia's RUR implementation computes the minimal polynomial and then **takes its square-free part** before using it for the rest of the RUR computation. This is the line:

```julia
f = f / Nemo.gcd(f, Nemo.derivative(f))  # Square-free part
```

## What This Means for Katsura4Degenerate

1. **Initial computation**: The minimal polynomial of degree 16 is computed (with multiplicities)
2. **Square-free reduction**: `f / gcd(f, f')` removes repeated factors
3. **Result**: A degree 1 polynomial (just T or T-c)

This explains why Julia returns degree 1 for our Katsura4Degenerate system!

## Mathematical Background

For a polynomial `f(T)`:
- If `f(T) = (T - a)^16` (our case with multiplicity 16)
- Then `f'(T) = 16(T - a)^15`
- So `gcd(f, f') = (T - a)^15`
- Therefore `f / gcd(f, f') = (T - a)` (degree 1!)

## Why Take the Square-Free Part?

The square-free part is used because:
1. **Algebraic correctness**: Working modulo the square-free part gives the correct field structure
2. **Numerical stability**: Avoids dealing with high multiplicities
3. **Efficiency**: Smaller degree polynomials for subsequent computations

## The RUR Components

Julia then computes the other RUR components (the parametrizations g_i) **modulo this square-free polynomial**. This is why we see:
- Minimal polynomial: [0, 1] (which is just T)
- Parametrizations: mostly empty or trivial

## For Our C++ Implementation

We have two choices:

### Option 1: Follow Julia's Approach
```cpp
// After computing minimal polynomial of degree 16
auto f_squared_free = compute_square_free_part(minimal_poly);
// Use f_squared_free for the rest of the computation
```

### Option 2: Keep Full Multiplicity Information
```cpp
// Use the full minimal polynomial with multiplicities
// This preserves more information but is computationally harder
```

## Testing Implications

For Katsura4Degenerate:
- The minimal polynomial computed is `(T - c)^16` for some constant c
- The square-free part is `(T - c)` (degree 1)
- This correctly identifies the single geometric solution
- The multiplicity 16 is "lost" in the square-free reduction

## The Deeper Question

Should RUR algorithms:
1. **Always use square-free parts** (Julia's approach)
   - Pros: Simpler, more stable, focuses on distinct solutions
   - Cons: Loses multiplicity information

2. **Preserve multiplicities** (theoretical approach)
   - Pros: Complete algebraic information
   - Cons: Computationally expensive, numerically unstable

3. **Offer both options** (flexible approach)
   - Let users choose based on their needs

## Conclusion

Julia's RUR is working correctly! It:
1. Computes the full minimal polynomial (degree 16)
2. Takes the square-free part (degree 1)
3. Uses this for the RUR parametrization

This is a deliberate design choice to handle non-radical ideals gracefully by focusing on the geometric solutions rather than the algebraic multiplicities.