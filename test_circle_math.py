#!/usr/bin/env python3
"""
Test which primes p make x^2 + y^2 - 1 behave specially.

Over F_p, x^2 + y^2 - 1 = 0 has solutions iff -1 is a sum of two squares mod p.
This is related to whether -1 is a quadratic residue mod p.
"""

def has_solutions_mod_p(p):
    """Check if x^2 + y^2 = 1 has solutions mod p."""
    solutions = []
    for x in range(p):
        for y in range(p):
            if (x*x + y*y) % p == 1:
                solutions.append((x, y))
    return len(solutions) > 0, len(solutions)

def check_minus_one_qr(p):
    """Check if -1 is a quadratic residue mod p."""
    # -1 is QR mod p iff p = 1 mod 4 (for odd primes)
    if p == 2:
        return True
    return p % 4 == 1

print("Testing x^2 + y^2 - 1 over various primes:")
print("=" * 60)

# Test small primes
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

for p in primes:
    has_sol, num_sol = has_solutions_mod_p(p)
    minus_one_qr = check_minus_one_qr(p)
    
    # Special analysis
    special = ""
    if not has_sol:
        special = " [NO SOLUTIONS!]"
    elif num_sol == p:
        special = " [TRIVIAL - all points satisfy]"
    
    print(f"p={p:3d}: {num_sol:3d} solutions, -1 QR: {str(minus_one_qr):5s}, p mod 4 = {p%4}{special}")

print("\n" + "=" * 60)
print("Pattern analysis:")
print("- Primes ≡ 3 (mod 4): -1 is NOT a quadratic residue")
print("- These primes have special behavior for x^2 + y^2")
print("\nProblematic primes (no solutions):")
problematic = []
for p in primes:
    has_sol, _ = has_solutions_mod_p(p)
    if not has_sol:
        problematic.append(p)
        
if problematic:
    print(f"  {problematic}")
else:
    print("  None found")

# Check if certain primes make the polynomial reduce to a constant
print("\n" + "=" * 60)
print("Checking if polynomial becomes constant:")
for p in [2, 3, 5, 7, 11]:
    print(f"\nPrime {p}:")
    # In F_p, check if x^2 + y^2 - 1 factors completely
    # or has some special structure
    if p == 2:
        # F_2: x^2 = x, y^2 = y, so x^2 + y^2 - 1 = x + y + 1
        print("  Over F_2: x^2 + y^2 - 1 = x + y + 1 (since x^2=x, y^2=y, -1=1)")
        print("  This is LINEAR, not constant!")
    elif p == 3:
        # Check if it has no solutions
        has_sol, num = has_solutions_mod_p(3)
        print(f"  Has {num} solutions mod 3")
        if not has_sol:
            print("  No solutions means ideal contains 1 (inconsistent)")