#include "../src/julia_rur/data_structures.hpp"
#include <iostream>
#include <cassert>
#include <sstream>

using namespace julia_rur;

void test_stackvect_basic() {
    std::cout << "Testing StackVect basic operations..." << std::endl;
    
    // Test parameterized constructor with empty monomial
    PP empty_monomial;
    StackVect sv1(0, empty_monomial, -1, -1);
    assert(sv1.pos == 0);
    assert(sv1.prev == -1);
    assert(sv1.var == -1);
    assert(sv1.mon.empty());
    
    // Test parameterized constructor
    PP monomial = {1, 2, 0};
    StackVect sv2(5, monomial, 3, 1);
    assert(sv2.pos == 5);
    assert(sv2.mon == monomial);
    assert(sv2.prev == 3);
    assert(sv2.var == 1);
    
    // Test equality
    StackVect sv3(5, monomial, 3, 1);
    assert(sv2 == sv3);
    assert(!(sv1 == sv2));
    
    // Test copy constructor
    StackVect sv4(sv2);
    assert(sv4 == sv2);
    
    // Test assignment
    StackVect sv5(0, empty_monomial, -1, -1);
    sv5 = sv2;
    assert(sv5 == sv2);
    
    std::cout << "✓ StackVect basic operations passed" << std::endl;
}

void test_stackvect_print() {
    std::cout << "Testing StackVect print functionality..." << std::endl;
    
    PP monomial = {1, 0, 2};
    StackVect sv(10, monomial, 7, 2);
    
    std::ostringstream oss;
    sv.print(oss);
    
    std::string expected = "StackVect{pos=10, mon=[1,0,2], prev=7, var=2}";
    assert(oss.str() == expected);
    
    std::cout << "✓ StackVect print functionality passed" << std::endl;
}

void test_power_product_operations() {
    std::cout << "Testing power product operations..." << std::endl;
    
    using namespace power_product;
    
    // Test zero monomial
    PP zero = zero_monomial(3);
    assert(zero.size() == 3);
    assert(zero[0] == 0 && zero[1] == 0 && zero[2] == 0);
    
    // Test unit monomial
    PP unit1 = unit_monomial(3, 1);
    assert(unit1.size() == 3);
    assert(unit1[0] == 0 && unit1[1] == 1 && unit1[2] == 0);
    
    // Test multiplication
    PP a = {1, 2, 0};
    PP b = {0, 1, 3};
    PP product = multiply(a, b);
    PP expected_product = {1, 3, 3};
    assert(product == expected_product);
    
    // Test division
    PP divisor = {1, 1, 0};
    PP dividend = {2, 3, 1};
    assert(divides(divisor, dividend));
    assert(!divides(dividend, divisor));
    
    // Test total degree
    PP monomial = {2, 1, 3};
    assert(total_degree(monomial) == 6);
    assert(total_degree(zero) == 0);
    
    // Test degree reverse lexicographic comparison
    PP m1 = {1, 1};  // degree 2
    PP m2 = {2, 0};  // degree 2 
    PP m3 = {0, 1};  // degree 1
    
    // m3 < m1, m2 (lower degree)
    assert(compare_degrevlex(m3, m1) < 0);
    assert(compare_degrevlex(m3, m2) < 0);
    
    // Debug: print comparison
    std::cout << "    Comparing m1=[1,1] vs m2=[2,0]: " << compare_degrevlex(m1, m2) << std::endl;
    
    // For same degree, reverse lexicographic: [2,0] > [1,1] in degrevlex
    // (smaller exponent in last variable comes first)
    assert(compare_degrevlex(m2, m1) > 0);
    
    // Self-comparison
    assert(compare_degrevlex(m1, m1) == 0);
    
    std::cout << "✓ Power product operations passed" << std::endl;
}

void test_quotient_ring_data_basic() {
    std::cout << "Testing QuotientRingData basic operations..." << std::endl;
    
    // Test constructor
    QuotientRingData qr(3, 100003);
    assert(qr.nvars == 3);
    assert(qr.prime == 100003);
    assert(qr.quotient_basis_size() == 0);
    assert(qr.border_size() == 0);
    assert(qr.i_xw.size() == 3);
    
    // Test adding quotient basis elements
    qr.quo.push_back({0, 0, 0});  // 1
    qr.quo.push_back({1, 0, 0});  // x
    qr.quo.push_back({0, 1, 0});  // y
    qr.quo.push_back({0, 0, 1});  // z
    
    assert(qr.quotient_basis_size() == 4);
    assert(qr.total_elements() == 4);
    
    // Test adding border elements
    PP border_mon1 = {2, 0, 0};  // x^2
    PP border_mon2 = {0, 2, 0};  // y^2
    
    qr.t_xw.emplace_back(0, border_mon1, -1, 0);
    qr.t_xw.emplace_back(1, border_mon2, -1, 1);
    
    assert(qr.border_size() == 2);
    assert(qr.total_elements() == 6);
    
    // Test validation
    assert(qr.is_valid());
    
    std::cout << "✓ QuotientRingData basic operations passed" << std::endl;
}

void test_quotient_ring_data_validation() {
    std::cout << "Testing QuotientRingData validation..." << std::endl;
    
    QuotientRingData qr(2, 100003);
    
    // Empty quotient basis should be invalid
    assert(!qr.is_valid());
    
    // Add valid quotient basis
    qr.quo.push_back({0, 0});  // 1
    qr.quo.push_back({1, 0});  // x
    assert(qr.is_valid());
    
    // Add monomial with wrong dimension - should be invalid
    qr.quo.push_back({1, 0, 1});  // Wrong dimension
    assert(!qr.is_valid());
    
    // Fix it
    qr.quo.pop_back();
    qr.quo.push_back({0, 1});  // y
    assert(qr.is_valid());
    
    // Add border element with wrong dimension
    PP wrong_dim = {1, 0, 1};  // 3 variables instead of 2
    qr.t_xw.emplace_back(0, wrong_dim, -1, 0);
    assert(!qr.is_valid());
    
    std::cout << "✓ QuotientRingData validation passed" << std::endl;
}

void test_quotient_ring_data_print() {
    std::cout << "Testing QuotientRingData print functionality..." << std::endl;
    
    QuotientRingData qr(2, 100003);
    qr.quo.push_back({0, 0});  // 1
    qr.quo.push_back({1, 0});  // x
    qr.quo.push_back({0, 1});  // y
    qr.quo.push_back({1, 1});  // xy
    
    // Add a border element
    PP border_mon = {2, 0};  // x^2
    qr.t_xw.emplace_back(0, border_mon, -1, 0);
    
    // Test print (just verify it doesn't crash)
    std::ostringstream oss;
    qr.print_summary(oss);
    
    std::string output = oss.str();
    assert(output.find("nvars: 2") != std::string::npos);
    assert(output.find("prime: 100003") != std::string::npos);
    assert(output.find("quotient_basis_size: 4") != std::string::npos);
    assert(output.find("border_size: 1") != std::string::npos);
    
    std::cout << "✓ QuotientRingData print functionality passed" << std::endl;
}

void test_quotient_ring_data_memory_management() {
    std::cout << "Testing QuotientRingData memory management..." << std::endl;
    
    QuotientRingData qr(2, 100003);
    
    // Test reserve functionality
    qr.reserve_quotient_basis(100);
    qr.reserve_border(50);
    
    // Add elements and verify capacity
    for (int i = 0; i < 10; ++i) {
        qr.quo.push_back({static_cast<uint32_t>(i), 0});
        PP border_mon = {static_cast<uint32_t>(i + 1), 1};
        qr.t_xw.emplace_back(i, border_mon, -1, 0);
    }
    
    assert(qr.quotient_basis_size() == 10);
    assert(qr.border_size() == 10);
    
    // Test clear functionality
    qr.clear();
    assert(qr.quotient_basis_size() == 0);
    assert(qr.border_size() == 0);
    assert(qr.i_xw.size() == 2);  // Should still have nvars vectors
    
    for (const auto& var_indices : qr.i_xw) {
        assert(var_indices.empty());
    }
    
    std::cout << "✓ QuotientRingData memory management passed" << std::endl;
}

void test_power_product_print() {
    std::cout << "Testing power product print functionality..." << std::endl;
    
    using namespace power_product;
    
    PP monomial = {1, 0, 2};
    std::ostringstream oss;
    print_monomial(monomial, oss);
    
    assert(oss.str() == "[1,0,2]");
    
    // Test empty monomial
    PP empty_monomial;
    std::ostringstream oss2;
    print_monomial(empty_monomial, oss2);
    assert(oss2.str() == "[]");
    
    std::cout << "✓ Power product print functionality passed" << std::endl;
}

int main() {
    try {
        test_stackvect_basic();
        test_stackvect_print();
        test_power_product_operations();
        test_quotient_ring_data_basic();
        test_quotient_ring_data_validation();
        test_quotient_ring_data_print();
        test_quotient_ring_data_memory_management();
        test_power_product_print();
        
        std::cout << "\n🎉 All Julia data structure tests passed!" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with error: " << e.what() << std::endl;
        return 1;
    }
}