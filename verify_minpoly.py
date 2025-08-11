#!/usr/bin/env python3
"""
Verification script for minimal polynomial computation.
Tests whether computed polynomials actually evaluate to zero at the separating element.
"""

def mod_inverse(a, m):
    """Compute modular inverse of a modulo m"""
    def extended_gcd(a, b):
        if a == 0:
            return b, 0, 1
        gcd, x1, y1 = extended_gcd(b % a, a)
        x = y1 - (b // a) * x1
        y = x1
        return gcd, x, y
    
    gcd, x, _ = extended_gcd(a % m, m)
    if gcd != 1:
        raise ValueError(f"Modular inverse of {a} mod {m} does not exist")
    return (x % m + m) % m

def polynomial_evaluation_finite_field(poly_coeffs, element_powers, prime):
    """
    Evaluate polynomial p(T) = sum(coeffs[i] * T^i) at T = separating_element
    where T^i is represented as element_powers[i] in the quotient basis.
    """
    result = [0] * len(element_powers[0])
    
    for i, coeff in enumerate(poly_coeffs):
        if coeff != 0 and i < len(element_powers):
            for j in range(len(result)):
                result[j] = (result[j] + coeff * element_powers[i][j]) % prime
    
    return result

def test_minimal_polynomial():
    """Test the minimal polynomials from our C++ implementations"""
    
    prime = 1073741827
    
    # Separating element and its powers from test output
    separating_element = [1073741815, 3, 1, 0, 0, 0]
    
    # Powers T^0, T^1, T^2 from test output
    powers = [
        [1, 0, 0, 0, 0, 0],                    # T^0
        [1073741815, 3, 1, 0, 0, 0],          # T^1  
        [112, 1073741791, 1073741815, 0, 0, 0] # T^2
    ]
    
    print("=== Minimal Polynomial Verification ===")
    print(f"Prime: {prime}")
    print(f"Separating element: {separating_element}")
    print()
    
    # Test OLD implementation polynomial: [1073741715, 36, 1]  
    old_coeffs = [1073741715, 36, 1]
    print("OLD polynomial coefficients:", old_coeffs)
    
    result_old = polynomial_evaluation_finite_field(old_coeffs, powers, prime)
    print("OLD polynomial evaluation result:", result_old)
    
    is_zero_old = all(x == 0 for x in result_old)
    print("OLD polynomial evaluates to zero:", is_zero_old)
    print()
    
    # Test NEW implementation polynomial: [32, 12, 1]
    new_coeffs = [32, 12, 1]
    print("NEW polynomial coefficients:", new_coeffs)
    
    result_new = polynomial_evaluation_finite_field(new_coeffs, powers, prime)
    print("NEW polynomial evaluation result:", result_new)
    
    is_zero_new = all(x == 0 for x in result_new)
    print("NEW polynomial evaluates to zero:", is_zero_new)
    print()
    
    # Manual verification: compute T^2 manually
    print("=== Manual verification of T^2 ===")
    
    # T^1 * T^1 should equal our T^2
    # This uses the multiplication tables from the quotient ring
    t1 = powers[1]
    
    # For manual computation, we'd need the actual multiplication tables
    # For now, just verify our computed T^2 matches
    print("T^1:", t1)
    print("T^2 (computed):", powers[2])
    
    # Test the polynomial T^2 + 12*T + 32 = 0
    manual_result = []
    for i in range(len(powers[0])):
        val = (powers[2][i] + 12 * powers[1][i] + 32 * powers[0][i]) % prime
        manual_result.append(val)
    
    print("Manual T^2 + 12*T + 32:", manual_result)
    print("Manual result is zero:", all(x == 0 for x in manual_result))
    
    return is_zero_new

if __name__ == "__main__":
    success = test_minimal_polynomial()
    if success:
        print("\n✓ FLINT implementation appears correct!")
    else:
        print("\n✗ FLINT implementation has issues")