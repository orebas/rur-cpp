#ifndef BIVARIATE_ALGORITHM_HPP
#define BIVARIATE_ALGORITHM_HPP

#include "data_structures.hpp"
#include "univariate_parameterization.hpp" // For MinimalPolynomialResult

namespace julia_rur {

// Represents a monomial in two variables (T, xi)
struct BivariateMonomial {
    int32_t deg_T;
    int32_t deg_xi;
};

// Result of the bivariate lexicographic algorithm
struct BivLexResult {
    std::vector<BivariateMonomial> monomial_basis;
    std::vector<BivariateMonomial> leading_monomials;
    std::vector<std::vector<ModularCoeff>> generators;
};

// Context for Gaussian elimination
struct GaussianReductionContext {
    std::vector<std::vector<ModularCoeff>> matrix;
    std::vector<int32_t> pivot_columns;
};

// The main biv_lex function
BivLexResult
biv_lex(const std::vector<std::vector<ModularCoeff>> &t_v,
        const std::vector<std::vector<int32_t>> &i_xw,
        GaussianReductionContext &context,
        const std::vector<std::vector<ModularCoeff>> &free_set,
        int32_t var_index,
        ModularCoeff prime);

} // namespace julia_rur

#endif // BIVARIATE_ALGORITHM_HPP