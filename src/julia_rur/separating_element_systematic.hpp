#pragma once

#include "data_structures.hpp"
#include "univariate_parameterization.hpp"
#include <tuple>
#include <vector>

namespace julia_rur {

// Systematic search for a separating element
std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>, std::vector<int>>
find_separating_element_systematic(std::vector<PP> quotient_basis,
                                   const std::vector<std::vector<int32_t>> &i_xw,
                                   std::vector<std::vector<ModularCoeff>> &t_v,
                                   int num_variables,
                                   ModularCoeff prime,
                                   int max_coeffs = 5);

} // namespace julia_rur
