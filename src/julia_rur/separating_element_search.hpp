#pragma once

#include "data_structures.hpp"
#include "univariate_parameterization.hpp"
#include <tuple>
#include <vector>

namespace julia_rur {


// Wrapper that retries with different separating elements until success
// Note: quotient_basis is now passed by value to allow reordering
std::tuple<bool, MinimalPolynomialResult, std::vector<BivariateResult>>
compute_univariate_parameterization_with_retry(std::vector<PP> quotient_basis,
                                               const std::vector<std::vector<int32_t>> &i_xw,
                                               std::vector<std::vector<ModularCoeff>> &t_v,
                                               int num_variables,
                                               ModularCoeff prime,
                                               int max_retries = 10);

} // namespace julia_rur