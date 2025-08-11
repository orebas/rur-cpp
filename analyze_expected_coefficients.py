#!/usr/bin/env python3
"""
Analyze the expected size of minimal polynomial coefficients for the system:
x² - 1 = 0, y² - 2 = 0, z² - 3 = 0

This system has 8 solutions: (±1, ±√2, ±√3)
"""

import sympy as sp
from sympy import symbols, expand, Poly, sqrt, I
import itertools

def analyze_system():
    """Analyze the 3-variable system and expected coefficient sizes"""
    print("=== Analysis of x² - 1, y² - 2, z² - 3 ===")
    
    # The 8 solutions are all combinations of signs
    solutions = []
    for sx, sy, sz in itertools.product([1, -1], repeat=3):
        x_val = sx * 1
        y_val = sy * sqrt(2) 
        z_val = sz * sqrt(3)
        solutions.append((x_val, y_val, z_val))
    
    print(f"Solutions ({len(solutions)} total):")
    for i, sol in enumerate(solutions):
        print(f"  {i+1}: {sol}")
    
    # Test potential separating elements
    print("\n=== Testing separating elements ===")
    
    separating_elements = [
        ("x", lambda x,y,z: x),
        ("y", lambda x,y,z: y), 
        ("z", lambda x,y,z: z),
        ("x+y", lambda x,y,z: x+y),
        ("x+z", lambda x,y,z: x+z),
        ("y+z", lambda x,y,z: y+z),
        ("x+y+z", lambda x,y,z: x+y+z),
        ("2x+y+z", lambda x,y,z: 2*x+y+z)
    ]
    
    for name, func in separating_elements:
        values = []
        for sol in solutions:
            val = func(*sol)
            values.append(val)
        
        unique_values = set(values)
        print(f"  {name}: {len(unique_values)} unique values")
        if len(unique_values) == 8:  # Separating
            print(f"    ✓ SEPARATING: {sorted(values, key=str)}")
            
            # Try to compute minimal polynomial
            try:
                t = symbols('t')
                # Build polynomial with these roots
                poly = 1
                for val in unique_values:
                    poly = expand(poly * (t - val))
                
                poly_coeffs = Poly(poly, t).all_coeffs()
                print(f"    Minimal polynomial degree: {len(poly_coeffs)-1}")
                print(f"    Coefficients: {poly_coeffs}")
                
                # Check coefficient magnitudes
                max_coeff = max(abs(c) for c in poly_coeffs if c != 0)
                print(f"    Max coefficient magnitude: {max_coeff}")
                print(f"    Max coefficient (float): {float(max_coeff)}")
                
            except Exception as e:
                print(f"    Error computing minimal polynomial: {e}")
            print()

def test_specific_separating_element():
    """Test a specific separating element that might match our C++ output"""
    print("\n=== Testing specific separating elements ===")
    
    # Based on our C++ debug output, we found a separating element
    # Let's see if we can identify it
    solutions = []
    for sx, sy, sz in itertools.product([1, -1], repeat=3):
        x_val = sx * 1
        y_val = sy * sqrt(2) 
        z_val = sz * sqrt(3)
        solutions.append((x_val, y_val, z_val))
    
    # Try various linear combinations
    test_combinations = [
        (0, 1, 1),   # y + z 
        (0, -1, 1),  # -y + z
        (1, 1, 1),   # x + y + z
        (2, 1, 1),   # 2x + y + z
        (1, 2, 3),   # x + 2y + 3z
    ]
    
    for a, b, c in test_combinations:
        print(f"\nTesting separating element {a}*x + {b}*y + {c}*z:")
        values = []
        for x_val, y_val, z_val in solutions:
            val = a*x_val + b*y_val + c*z_val
            values.append(val)
        
        unique_values = list(set(values))
        print(f"  Values: {sorted(values, key=str)}")
        print(f"  Unique: {len(unique_values)} ({'separating' if len(unique_values) == 8 else 'NOT separating'})")
        
        if len(unique_values) == 8:
            try:
                t = symbols('t')
                poly = 1
                for val in unique_values:
                    poly = expand(poly * (t - val))
                
                coeffs = Poly(poly, t).all_coeffs()
                print(f"  Minimal polynomial: {poly}")
                print(f"  Coefficients: {coeffs}")
                
                # Convert to numerical approximations
                numerical_coeffs = [complex(c.evalf()) for c in coeffs]
                print(f"  Numerical coefficients: {numerical_coeffs}")
                
            except Exception as e:
                print(f"  Error: {e}")

if __name__ == "__main__":
    analyze_system()
    test_specific_separating_element()