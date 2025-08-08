#include "polynomial_evaluator.hpp"
#include <cctype>
#include <cmath>
#include <sstream>
#include <algorithm>

namespace julia_rur {

std::complex<double> PolynomialEvaluator::evaluate(
    const std::string& expression,
    const std::vector<std::string>& variables,
    const std::vector<std::complex<double>>& values
) {
    if (variables.size() != values.size()) {
        throw std::runtime_error("Number of variables must match number of values");
    }
    
    // Create variable map
    std::map<std::string, std::complex<double>> var_map;
    for (size_t i = 0; i < variables.size(); ++i) {
        var_map[variables[i]] = values[i];
    }
    
    // Tokenize expression
    auto tokens = tokenize(expression);
    
    // Parse and evaluate
    size_t pos = 0;
    auto result = parse_expression(tokens, pos, var_map);
    
    if (pos < tokens.size() && tokens[pos].type != TokenType::END) {
        throw std::runtime_error("Unexpected token after expression");
    }
    
    return result;
}

double PolynomialEvaluator::compute_residual(
    const std::vector<std::string>& polynomials,
    const std::vector<std::string>& variables,
    const std::vector<std::complex<double>>& solution
) {
    double max_residual = 0.0;
    
    for (const auto& poly : polynomials) {
        try {
            auto value = evaluate(poly, variables, solution);
            double residual = std::abs(value);
            max_residual = std::max(max_residual, residual);
        } catch (const std::exception& e) {
            // If evaluation fails, return a large residual
            return 1e10;
        }
    }
    
    return max_residual;
}

std::vector<Token> PolynomialEvaluator::tokenize(const std::string& expression) {
    std::vector<Token> tokens;
    size_t pos = 0;
    
    while (pos < expression.length()) {
        Token token = next_token(expression, pos);
        if (token.type == TokenType::END) break;
        tokens.push_back(token);
    }
    
    tokens.push_back({TokenType::END, "", 0.0});
    return tokens;
}

Token PolynomialEvaluator::next_token(const std::string& expr, size_t& pos) {
    skip_whitespace(expr, pos);
    
    if (pos >= expr.length()) {
        return {TokenType::END, "", 0.0};
    }
    
    char ch = expr[pos];
    
    // Single character tokens
    if (ch == '+') {
        pos++;
        return {TokenType::PLUS, "+", 0.0};
    } else if (ch == '-') {
        pos++;
        return {TokenType::MINUS, "-", 0.0};
    } else if (ch == '*') {
        pos++;
        return {TokenType::MULTIPLY, "*", 0.0};
    } else if (ch == '^') {
        pos++;
        return {TokenType::POWER, "^", 0.0};
    } else if (ch == '(') {
        pos++;
        return {TokenType::LPAREN, "(", 0.0};
    } else if (ch == ')') {
        pos++;
        return {TokenType::RPAREN, ")", 0.0};
    }
    
    // Numbers
    if (std::isdigit(ch) || ch == '.') {
        return parse_number(expr, pos);
    }
    
    // Variables
    if (std::isalpha(ch)) {
        return parse_variable(expr, pos);
    }
    
    throw std::runtime_error("Unexpected character: " + std::string(1, ch));
}

void PolynomialEvaluator::skip_whitespace(const std::string& expr, size_t& pos) {
    while (pos < expr.length() && std::isspace(expr[pos])) {
        pos++;
    }
}

Token PolynomialEvaluator::parse_number(const std::string& expr, size_t& pos) {
    size_t start = pos;
    bool has_dot = false;
    
    while (pos < expr.length() && (std::isdigit(expr[pos]) || expr[pos] == '.')) {
        if (expr[pos] == '.') {
            if (has_dot) break;  // Second dot
            has_dot = true;
        }
        pos++;
    }
    
    std::string num_str = expr.substr(start, pos - start);
    double value = std::stod(num_str);
    
    return {TokenType::NUMBER, num_str, value};
}

Token PolynomialEvaluator::parse_variable(const std::string& expr, size_t& pos) {
    size_t start = pos;
    
    // First character is letter
    pos++;
    
    // Continue with letters or digits
    while (pos < expr.length() && (std::isalnum(expr[pos]) || expr[pos] == '_')) {
        pos++;
    }
    
    std::string var_name = expr.substr(start, pos - start);
    return {TokenType::VARIABLE, var_name, 0.0};
}

std::complex<double> PolynomialEvaluator::parse_expression(
    const std::vector<Token>& tokens,
    size_t& pos,
    const std::map<std::string, std::complex<double>>& var_map
) {
    auto result = parse_term(tokens, pos, var_map);
    
    while (pos < tokens.size()) {
        if (tokens[pos].type == TokenType::PLUS) {
            pos++;
            result += parse_term(tokens, pos, var_map);
        } else if (tokens[pos].type == TokenType::MINUS) {
            pos++;
            result -= parse_term(tokens, pos, var_map);
        } else {
            break;
        }
    }
    
    return result;
}

std::complex<double> PolynomialEvaluator::parse_term(
    const std::vector<Token>& tokens,
    size_t& pos,
    const std::map<std::string, std::complex<double>>& var_map
) {
    auto result = parse_factor(tokens, pos, var_map);
    
    while (pos < tokens.size() && tokens[pos].type == TokenType::MULTIPLY) {
        pos++;
        result *= parse_factor(tokens, pos, var_map);
    }
    
    return result;
}

std::complex<double> PolynomialEvaluator::parse_factor(
    const std::vector<Token>& tokens,
    size_t& pos,
    const std::map<std::string, std::complex<double>>& var_map
) {
    // Handle unary minus
    bool negative = false;
    if (pos < tokens.size() && tokens[pos].type == TokenType::MINUS) {
        negative = true;
        pos++;
    } else if (pos < tokens.size() && tokens[pos].type == TokenType::PLUS) {
        pos++;  // Skip unary plus
    }
    
    auto result = parse_power(tokens, pos, var_map);
    
    if (negative) {
        result = -result;
    }
    
    return result;
}

std::complex<double> PolynomialEvaluator::parse_power(
    const std::vector<Token>& tokens,
    size_t& pos,
    const std::map<std::string, std::complex<double>>& var_map
) {
    auto base = parse_atom(tokens, pos, var_map);
    
    if (pos < tokens.size() && tokens[pos].type == TokenType::POWER) {
        pos++;
        auto exponent = parse_atom(tokens, pos, var_map);
        
        // Handle integer exponents specially for accuracy
        if (exponent.imag() == 0.0 && std::floor(exponent.real()) == exponent.real()) {
            int exp = static_cast<int>(exponent.real());
            if (exp >= 0) {
                std::complex<double> result = 1.0;
                for (int i = 0; i < exp; ++i) {
                    result *= base;
                }
                return result;
            }
        }
        
        // General case
        return std::pow(base, exponent);
    }
    
    return base;
}

std::complex<double> PolynomialEvaluator::parse_atom(
    const std::vector<Token>& tokens,
    size_t& pos,
    const std::map<std::string, std::complex<double>>& var_map
) {
    if (pos >= tokens.size()) {
        throw std::runtime_error("Unexpected end of expression");
    }
    
    const Token& token = tokens[pos];
    
    if (token.type == TokenType::NUMBER) {
        pos++;
        return std::complex<double>(token.number_value, 0.0);
    }
    
    if (token.type == TokenType::VARIABLE) {
        pos++;
        auto it = var_map.find(token.value);
        if (it == var_map.end()) {
            throw std::runtime_error("Unknown variable: " + token.value);
        }
        return it->second;
    }
    
    if (token.type == TokenType::LPAREN) {
        pos++;
        auto result = parse_expression(tokens, pos, var_map);
        if (pos >= tokens.size() || tokens[pos].type != TokenType::RPAREN) {
            throw std::runtime_error("Missing closing parenthesis");
        }
        pos++;
        return result;
    }
    
    throw std::runtime_error("Unexpected token: " + token.value);
}

} // namespace julia_rur