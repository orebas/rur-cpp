#include <iostream>

int main() {
    unsigned int p = 1073741827;
    
    std::cout << "Analysis of prime 1073741827:" << std::endl;
    std::cout << "Binary: ";
    for (int i = 31; i >= 0; i--) {
        std::cout << ((p >> i) & 1);
        if (i % 4 == 0) std::cout << " ";
    }
    std::cout << std::endl;
    
    std::cout << "Hex: 0x" << std::hex << p << std::dec << std::endl;
    
    // Check if it's 2^k + c
    for (int k = 2; k <= 31; k++) {
        unsigned long long two_k = 1ULL << k;
        if (two_k > p) {
            std::cout << "Form: 2^" << k << " - " << (two_k - p) << std::endl;
            break;
        }
    }
    
    std::cout << "mod 8 = " << (p % 8) << std::endl;
    
    // Check the "pairs limit" progression
    std::cout << "\nNotice the 'pairs limit' progression in F4 output:" << std::endl;
    std::cout << "1073741827: pairs limit 2048" << std::endl;
    std::cout << "1073741831: pairs limit 4096" << std::endl;
    std::cout << "1073741833: pairs limit 8192" << std::endl;
    std::cout << "etc. (doubling each time)" << std::endl;
    std::cout << "\nThis suggests F4 uses different code paths based on prime!" << std::endl;
    
    return 0;
}