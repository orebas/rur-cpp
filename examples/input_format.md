# RUR Input Format Examples

This document describes the various input formats supported by the RUR solver.

## Basic Polynomial Format

### String Format (Recommended)

Polynomials are input as strings with explicit coefficients:

```cpp
// Format: coefficient*variable^power
"1*x^2+2*x+3"           // x² + 2x + 3
"1*x^2*y+3*x*y^2"       // x²y + 3xy²
"5*x^3+100002*x"        // 5x³ - x (using 100002 ≡ -1 mod 100003)
```

Rules:
- Always include coefficient (use `1*x` not just `x`)
- Powers use caret: `x^2`, `y^3`
- Multiple variables: `x^2*y^3`
- Space around + and - is optional
- Negative coefficients in modular form: `-1` → `100002` (mod 100003)

### JSON Input Format

For batch processing or programmatic use:

```json
{
  "system": {
    "polynomials": [
      "1*x^2+1*y^2+100002",
      "1*x+100002*y"
    ],
    "variables": ["x", "y"],
    "prime": 100003,
    "options": {
      "verbose": true,
      "separating_strategy": "CURRENT"
    }
  }
}
```

### File Input Format

Save systems in `.rur` files:

```
# Circle and line intersection
# This is a comment
VARIABLES: x, y
PRIME: 100003
POLYNOMIALS:
1*x^2+1*y^2+100002
1*x+100002*y
END
```

## Common Examples

### 1. Linear System

```cpp
// System: 2x + 3y = 5, x - y = 1
// Rewrite as: 2x + 3y - 5 = 0, x - y - 1 = 0
std::vector<std::string> polynomials = {
    "2*x+3*y+99998",    // 2x + 3y - 5 (99998 ≡ -5 mod 100003)
    "1*x+100002*y+100002" // x - y - 1
};
```

### 2. Quadratic System

```cpp
// System: x² + y² = 25, x² - y² = 7
std::vector<std::string> polynomials = {
    "1*x^2+1*y^2+99978",    // x² + y² - 25
    "1*x^2+100002*y^2+99996" // x² - y² - 7
};
```

### 3. Mixed Degree System

```cpp
// System: x³ + y = 0, x + y³ = 0
std::vector<std::string> polynomials = {
    "1*x^3+1*y",
    "1*x+1*y^3"
};
```

### 4. Three Variable System

```cpp
// Cyclic-3: x+y+z=0, xy+yz+zx=0, xyz=1
std::vector<std::string> polynomials = {
    "1*x+1*y+1*z",
    "1*x*y+1*y*z+1*z*x",
    "1*x*y*z+100002"    // xyz - 1
};
```

### 5. Parametric System

```cpp
// With parameter a: x² + y² = a², x + y = a
// Substitute specific value, e.g., a = 2
std::vector<std::string> polynomials = {
    "1*x^2+1*y^2+99999",    // x² + y² - 4
    "1*x+1*y+100001"        // x + y - 2
};
```

## Modular Arithmetic Reference

For prime p = 100003:
- -1 → 100002
- -2 → 100001
- -3 → 100000
- -4 → 99999
- -5 → 99998
- ...
- -36 → 99967

Formula: negative value n → p + n

## Variable Ordering

The order of variables in the input affects the output:

```cpp
std::vector<std::string> variables = {"x", "y", "z"};
// Output will parameterize in order: x(T), y(T), z(T)

std::vector<std::string> variables = {"z", "y", "x"};
// Output will parameterize in order: z(T), y(T), x(T)
```

## Special Cases

### Empty System
```cpp
std::vector<std::string> polynomials = {};
// Error: No polynomials
```

### Inconsistent System
```cpp
std::vector<std::string> polynomials = {"1"};  // Just the constant 1
// Result: No solutions (Gröbner basis = {1})
```

### Under-determined System
```cpp
std::vector<std::string> polynomials = {"1*x+1*y"};  // Single equation
std::vector<std::string> variables = {"x", "y"};
// Error: Not zero-dimensional (infinitely many solutions)
```

### Over-determined System
```cpp
// Three equations, two unknowns (may still have finite solutions)
std::vector<std::string> polynomials = {
    "1*x^2+1*y^2+100002",
    "1*x+100002*y",
    "1*x*y+100002"
};
```

## Programmatic Input Generation

```cpp
// Helper function to convert integer coefficient to modular form
std::string mod_coeff(int c, int prime = 100003) {
    if (c < 0) {
        c = prime + c;  // Convert negative to positive modular form
    }
    return std::to_string(c);
}

// Build polynomial programmatically
std::string build_poly(const std::vector<std::tuple<int, std::string>>& terms,
                      int prime = 100003) {
    std::string poly;
    bool first = true;
    
    for (const auto& [coeff, monomial] : terms) {
        if (!first) poly += "+";
        poly += mod_coeff(coeff, prime) + "*" + monomial;
        first = false;
    }
    
    return poly;
}

// Example usage
auto p1 = build_poly({{1, "x^2"}, {1, "y^2"}, {-1, "1"}});  // x² + y² - 1
```

## Error Handling

Common input errors and their messages:

1. **Invalid polynomial syntax**
   - Missing `*`: "2x" → Error, use "2*x"
   - Invalid character: "x²" → Error, use "x^2"

2. **Unknown variable**
   - Polynomial uses 'z' but variables only lists ["x", "y"]

3. **Prime too small**
   - Prime must be > 65536 for F4 algorithm

4. **Coefficient overflow**
   - Coefficients must fit in 32-bit unsigned integer