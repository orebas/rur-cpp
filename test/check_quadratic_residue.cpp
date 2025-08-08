#include <iostream>
#include <vector>

// Check if a is a quadratic residue mod p using Euler's criterion
bool is_quadratic_residue(unsigned int a, unsigned int p) {
    if (a % p == 0) return true;
    
    // Compute a^((p-1)/2) mod p
    unsigned long long exp = (p - 1) / 2;
    unsigned long long result = 1;
    unsigned long long base = a % p;
    unsigned long long mod = p;
    
    while (exp > 0) {
        if (exp & 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp >>= 1;
    }
    
    return result == 1;
}

int main() {
    std::cout << "Checking if 2 is a quadratic residue for primes around 2^28:\n" << std::endl;
    
    // Check some primes around 2^28
    std::vector<unsigned int> primes = {
        268435459, 268435463, 268435493, 268435537, 268435561,
        268435577, 268435579, 268435597, 268435607, 268435631,
        268435649, 268435657, 268435661, 268435667, 268435669,
        268435693, 268435697, 268435721, 268435741, 268435757
    };
    
    int good_count = 0;
    int bad_count = 0;
    
    for (auto p : primes) {
        bool is_qr = is_quadratic_residue(2, p);
        if (is_qr) {
            std::cout << p << ": 2 is QR ✓" << std::endl;
            good_count++;
        } else {
            std::cout << p << ": 2 is NOT QR ✗" << std::endl;
            bad_count++;
        }
    }
    
    std::cout << "\nSummary: " << good_count << " good primes, " << bad_count << " bad primes" << std::endl;
    std::cout << "Ratio: " << (100.0 * good_count / primes.size()) << "% good primes" << std::endl;
    
    // For primes p ≡ 1 or 7 (mod 8), 2 is a QR
    // For primes p ≡ 3 or 5 (mod 8), 2 is NOT a QR
    std::cout << "\nChecking mod 8 pattern:" << std::endl;
    for (auto p : primes) {
        int mod8 = p % 8;
        bool should_be_qr = (mod8 == 1 || mod8 == 7);
        bool is_qr = is_quadratic_residue(2, p);
        std::cout << p << " ≡ " << mod8 << " (mod 8), ";
        std::cout << "should be QR: " << (should_be_qr ? "yes" : "no") << ", ";
        std::cout << "is QR: " << (is_qr ? "yes" : "no");
        if (should_be_qr != is_qr) {
            std::cout << " ERROR!";
        }
        std::cout << std::endl;
    }
    
    return 0;
}