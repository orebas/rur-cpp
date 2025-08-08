#pragma once

#include <string>

namespace julia_rur {

/**
 * @brief Convert user-friendly polynomial string to F4 library format
 * 
 * The F4 library expects polynomials in a specific format:
 * - All terms must have explicit coefficients (e.g., "1*x" not just "x")
 * - Negative coefficients should be formatted carefully
 * - No spaces allowed
 * 
 * @param input User-provided polynomial string like "x^2 - 2"
 * @return F4-compatible string like "1*x^2-2"
 */
std::string format_polynomial_for_f4(const std::string& input);

} // namespace julia_rur