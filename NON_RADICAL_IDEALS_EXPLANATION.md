# Understanding Non-Radical Ideals in RUR Algorithms

## The Katsura4Degenerate Case Study

We discovered an interesting test case where our "Katsura-4" system was actually a degenerate variant with a non-radical ideal. This document explains what this means and how RUR algorithms should handle such cases.

## What is a Non-Radical Ideal?

An ideal I is **radical** if I = √I (its radical). For zero-dimensional ideals:
- The **quotient ring dimension** = dim(R/I) counts solutions WITH multiplicities
- The **radical degree** = dim(R/√I) counts DISTINCT solutions

In our Katsura4Degenerate system:
- Quotient ring dimension = 16
- Radical degree = 1
- This means: 1 distinct solution with total multiplicity 16!

## The Katsura4Degenerate System

```
f1: x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1
f2: x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1
f3: x0*x2 + 2*x1*x3 + 2*x2*x4 - x2
f4: x0*x3 + 2*x1*x4 - x3
f5: x0*x4 - x4  // <-- The problematic equation
```

The last equation factors as: x4*(x0 - 1) = 0

This creates a cascade effect:
- If x4 = 0, then f4 gives x3*(x0-1) = 0
- If x3 = 0, then f3 gives x2*(x0-1) = 0
- And so on...

The unique solution is (1, 0, 0, 0, 0) with multiplicity 16.

## How RUR Should Handle Non-Radical Ideals

### 1. **Julia's Approach (Current)**
Julia's RUR returns:
- Minimal polynomial: [0, 1] (just T)
- This is technically correct but not useful
- It indicates the separating element doesn't separate the single point

### 2. **Ideal Approach for RUR**

When detecting a non-radical ideal (quotient dim ≠ radical degree):

**Option A: Work with the Radical**
1. Compute the radical of the ideal
2. Apply RUR to the radical ideal
3. Report multiplicities separately

**Option B: Use a Generic Separating Element**
1. Use a random linear combination of variables
2. This should separate the "infinitesimally close" points
3. The minimal polynomial will have degree = quotient dimension

**Option C: Report Failure with Explanation**
1. Detect non-radical case early
2. Return a special status indicating non-radical ideal
3. Provide quotient dimension and radical degree

### 3. **Our C++ Implementation's Behavior**

Currently, our implementation:
1. Correctly computes quotient dimension = 16
2. Tries x4 as separating element
3. Gets stuck because x4 only achieves degree 5 (not fully separating)
4. Should eventually try random linear combinations

The issue: Even with multiplicities, we need a separating element that distinguishes the "tangent directions" at the single solution point.

## Mathematical Background

For a non-radical ideal with a single solution p with multiplicity m:
- The quotient ring R/I has dimension m
- It represents the "jet space" at p (derivatives up to order m-1)
- A generic linear form will separate these "infinitesimal" directions
- The minimal polynomial will have degree m

## Recommended Implementation Strategy

```cpp
// Pseudocode for handling non-radical ideals
if (quotient_dimension != radical_degree) {
    // Non-radical ideal detected
    if (config.handle_multiplicities) {
        // Try increasingly generic separating elements
        for (int attempts = 0; attempts < max_attempts; ++attempts) {
            auto sep = generate_random_linear_combination(attempts);
            auto [min_poly, param] = compute_rur_with_separator(sep);
            if (degree(min_poly) == quotient_dimension) {
                // Success! We've separated all infinitesimal directions
                return {min_poly, param, multiplicities_info};
            }
        }
    } else {
        // User doesn't want multiplicities
        // Compute and return RUR of the radical
        auto radical_ideal = compute_radical(ideal);
        return compute_rur(radical_ideal);
    }
}
```

## Testing Non-Radical Ideals

Good test cases for non-radical ideals:
1. **Katsura4Degenerate**: Our discovered case
2. **Double roots**: (x^2, y^2) - single point (0,0) with multiplicity 4
3. **Tangent intersections**: Curves meeting tangentially
4. **Fat points**: Ideals of the form (x-a)^k, (y-b)^k

## Conclusion

Non-radical ideals are mathematically interesting and computationally challenging. They arise naturally in:
- Singular solutions of polynomial systems
- Tangent intersections in geometry
- Degenerate configurations in kinematics
- Critical points with multiplicity

A robust RUR implementation should:
1. Detect non-radical cases
2. Provide options for handling them
3. Use sufficiently generic separating elements
4. Report multiplicities when relevant