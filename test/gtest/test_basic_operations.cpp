#include "monomial.hpp"
#include "polynomial.hpp"
#include "test_helpers.hpp"

using namespace rur_test;

/**
 * @brief Test basic monomial operations
 */
class BasicOperations : public ::testing::Test {
  protected:
    void SetUp() override {
        // Common test data
        varnames = { "x", "y", "z" };
    }

    std::vector<std::string> varnames;
};

TEST_F(BasicOperations, MonomialCreation) {
    Monomial m1({ 1, 2, 0 }); // x*y^2
    Monomial m2({ 0, 1, 1 }); // y*z

    EXPECT_EQ(m1.degree(), 3);
    EXPECT_EQ(m2.degree(), 2);

    // Test string representation
    std::string m1_str = m1.to_string(varnames);
    EXPECT_TRUE(m1_str.find("x") != std::string::npos);
    EXPECT_TRUE(m1_str.find("y") != std::string::npos);
}

TEST_F(BasicOperations, MonomialMultiplication) {
    Monomial m1({ 1, 2, 0 }); // x*y^2
    Monomial m2({ 0, 1, 1 }); // y*z

    Monomial product = m1 * m2; // x*y^3*z

    EXPECT_EQ(product.degree(), 5);
    EXPECT_EQ(product.exponents()[0], 1); // x^1
    EXPECT_EQ(product.exponents()[1], 3); // y^3
    EXPECT_EQ(product.exponents()[2], 1); // z^1
}

TEST_F(BasicOperations, MonomialComparison) {
    Monomial m1({ 1, 2, 0 });
    Monomial m2({ 1, 2, 0 });
    Monomial m3({ 2, 1, 0 });

    EXPECT_TRUE(m1 == m2);
    EXPECT_FALSE(m1 == m3);
    EXPECT_TRUE(m1 != m3);
}

TEST_F(BasicOperations, PolynomialConstruction) {
    MultivariatePolynomial<int> p;

    // Build polynomial: 3*x^2 - 2*x*y + 1
    p.set_coefficient(Monomial({ 2, 0, 0 }), 3);  // 3*x^2
    p.set_coefficient(Monomial({ 1, 1, 0 }), -2); // -2*x*y
    p.set_coefficient(Monomial({ 0, 0, 0 }), 1);  // +1

    EXPECT_FALSE(p.is_zero());
    EXPECT_EQ(p.nterms(), 3);

    // Check get_coefficients
    EXPECT_EQ(p.get_coefficient(Monomial({ 2, 0, 0 })), 3);
    EXPECT_EQ(p.get_coefficient(Monomial({ 1, 1, 0 })), -2);
    EXPECT_EQ(p.get_coefficient(Monomial({ 0, 0, 0 })), 1);

    // Check leading term
    EXPECT_EQ(p.leading_coefficient(), 3);
}

TEST_F(BasicOperations, PolynomialAddition) {
    MultivariatePolynomial<int> p1, p2;

    // p1 = x^2 + y
    p1.set_coefficient(Monomial({ 2, 0, 0 }), 1); // x^2
    p1.set_coefficient(Monomial({ 0, 1, 0 }), 1); // y

    // p2 = x^2 + z
    p2.set_coefficient(Monomial({ 2, 0, 0 }), 1); // x^2
    p2.set_coefficient(Monomial({ 0, 0, 1 }), 1); // z

    MultivariatePolynomial<int> sum = p1 + p2;

    // sum = 2*x^2 + y + z
    EXPECT_EQ(sum.get_coefficient(Monomial({ 2, 0, 0 })), 2); // 2*x^2
    EXPECT_EQ(sum.get_coefficient(Monomial({ 0, 1, 0 })), 1); // y
    EXPECT_EQ(sum.get_coefficient(Monomial({ 0, 0, 1 })), 1); // z
}

TEST_F(BasicOperations, PolynomialMultiplication) {
    MultivariatePolynomial<int> p1, p2;

    // p1 = x + 1
    p1.set_coefficient(Monomial({ 1, 0, 0 }), 1); // x
    p1.set_coefficient(Monomial({ 0, 0, 0 }), 1); // 1

    // p2 = x + y
    p2.set_coefficient(Monomial({ 1, 0, 0 }), 1); // x
    p2.set_coefficient(Monomial({ 0, 1, 0 }), 1); // y

    MultivariatePolynomial<int> product = p1 * p2;

    // product = x^2 + x*y + x + y
    EXPECT_EQ(product.get_coefficient(Monomial({ 2, 0, 0 })), 1); // x^2
    EXPECT_EQ(product.get_coefficient(Monomial({ 1, 1, 0 })), 1); // x*y
    EXPECT_EQ(product.get_coefficient(Monomial({ 1, 0, 0 })), 1); // x
    EXPECT_EQ(product.get_coefficient(Monomial({ 0, 1, 0 })), 1); // y
}

TEST_F(BasicOperations, PolynomialStringRepresentation) {
    MultivariatePolynomial<int> p;

    // Build: 3*x^2 - 2*x*y + 1
    p.set_coefficient(Monomial({ 2, 0, 0 }), 3);  // 3*x^2
    p.set_coefficient(Monomial({ 1, 1, 0 }), -2); // -2*x*y
    p.set_coefficient(Monomial({ 0, 0, 0 }), 1);  // +1

    std::string p_str = p.to_string(varnames);

    // Should contain all terms (order may vary)
    EXPECT_TRUE(p_str.find("x") != std::string::npos);
    EXPECT_TRUE(p_str.find("y") != std::string::npos);
    EXPECT_TRUE(p_str.find("3") != std::string::npos);
    EXPECT_TRUE(p_str.find("-2") != std::string::npos || p_str.find("- 2") != std::string::npos);
}

TEST_F(BasicOperations, ZeroPolynomial) {
    MultivariatePolynomial<int> zero_poly;

    EXPECT_TRUE(zero_poly.is_zero());
    EXPECT_EQ(zero_poly.nterms(), 0);
    EXPECT_EQ(zero_poly.leading_coefficient(), 0);
}

TEST_F(BasicOperations, PolynomialDegree) {
    MultivariatePolynomial<int> p;

    // Build: x^3*y^2 + x*z^4 + 1
    p.set_coefficient(Monomial({ 3, 2, 0 }), 1); // x^3*y^2 (degree 5)
    p.set_coefficient(Monomial({ 1, 0, 4 }), 1); // x*z^4 (degree 5)
    p.set_coefficient(Monomial({ 0, 0, 0 }), 1); // 1 (degree 0)

    EXPECT_EQ(p.degree(), 5);
}