#pragma once

#include <complex>
#include <vector>
#include <string>
#include <map>
#include <stdexcept>

namespace julia_rur {

/**
 * @brief Token types for polynomial parsing
 */
enum class TokenType {
    NUMBER,
    VARIABLE,
    PLUS,
    MINUS,
    MULTIPLY,
    POWER,
    LPAREN,
    RPAREN,
    END
};

/**
 * @brief Token structure for lexical analysis
 */
struct Token {
    TokenType type;
    std::string value;
    double number_value;
};

/**
 * @brief Simple polynomial parser and evaluator
 * 
 * Handles basic polynomial expressions with:
 * - Variables (x, y, z, x1, x2, etc.)
 * - Integer and decimal coefficients
 * - Basic operations: +, -, *, ^
 * - Parentheses for grouping
 */
class PolynomialEvaluator {
public:
    /**
     * @brief Evaluate a polynomial expression at given variable values
     * 
     * @param expression Polynomial string like "x^2 + 2*x*y - 3"
     * @param variables Variable names in order
     * @param values Complex values for each variable
     * @return Complex result of evaluation
     */
    static std::complex<double> evaluate(
        const std::string& expression,
        const std::vector<std::string>& variables,
        const std::vector<std::complex<double>>& values
    );
    
    /**
     * @brief Compute residual norm for a polynomial system at a solution
     * 
     * @param polynomials System of polynomial equations
     * @param variables Variable names
     * @param solution Values for each variable
     * @return Maximum absolute value among all polynomial evaluations
     */
    static double compute_residual(
        const std::vector<std::string>& polynomials,
        const std::vector<std::string>& variables,
        const std::vector<std::complex<double>>& solution
    );

private:
    // Lexical analysis
    static std::vector<Token> tokenize(const std::string& expression);
    static Token next_token(const std::string& expr, size_t& pos);
    static void skip_whitespace(const std::string& expr, size_t& pos);
    static Token parse_number(const std::string& expr, size_t& pos);
    static Token parse_variable(const std::string& expr, size_t& pos);
    
    // Parsing and evaluation
    static std::complex<double> parse_expression(
        const std::vector<Token>& tokens,
        size_t& pos,
        const std::map<std::string, std::complex<double>>& var_map
    );
    
    static std::complex<double> parse_term(
        const std::vector<Token>& tokens,
        size_t& pos,
        const std::map<std::string, std::complex<double>>& var_map
    );
    
    static std::complex<double> parse_factor(
        const std::vector<Token>& tokens,
        size_t& pos,
        const std::map<std::string, std::complex<double>>& var_map
    );
    
    static std::complex<double> parse_power(
        const std::vector<Token>& tokens,
        size_t& pos,
        const std::map<std::string, std::complex<double>>& var_map
    );
    
    static std::complex<double> parse_atom(
        const std::vector<Token>& tokens,
        size_t& pos,
        const std::map<std::string, std::complex<double>>& var_map
    );
};

} // namespace julia_rur