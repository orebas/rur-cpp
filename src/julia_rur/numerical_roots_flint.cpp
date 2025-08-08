#include "numerical_roots_flint.hpp"
#include <iomanip>
#include <iostream>


namespace julia_rur {

#ifdef HAVE_FLINT_ARB
// FLINT 3.x with ARB implementation

// --- Local helpers over Q[t] to check and compute square-free part ---
static void
trim_q(std::vector<mpq_class> &p) {
    while (!p.empty() && p.back() == 0) p.pop_back();
}
static std::vector<mpq_class>
deriv_q(const std::vector<mpq_class> &p) {
    if (p.size() <= 1) return { mpq_class(0) };
    std::vector<mpq_class> d(p.size() - 1, mpq_class(0));
    for (size_t i = 1; i < p.size(); ++i) d[i - 1] = p[i] * static_cast<int>(i);
    trim_q(d);
    return d;
}
static std::pair<std::vector<mpq_class>, std::vector<mpq_class>>
divmod_q(std::vector<mpq_class> a, const std::vector<mpq_class> &b) {
    trim_q(a);
    std::vector<mpq_class> q;
    if (b.empty() || (b.size() == 1 && b[0] == 0)) return { q, a };
    std::vector<mpq_class> r = a;
    std::vector<mpq_class> bb = b;
    trim_q(r);
    trim_q(bb);
    if (bb.empty()) return { q, r };
    mpq_class inv_lc = 1 / bb.back();
    while (r.size() >= bb.size() && !(r.size() == 1 && r[0] == 0)) {
        size_t k = r.size() - bb.size();
        mpq_class coeff = r.back() * inv_lc;
        if (q.size() < k + 1) q.resize(k + 1);
        q[k] += coeff;
        for (size_t i = 0; i < bb.size(); ++i) r[k + i] -= coeff * bb[i];
        trim_q(r);
        if (bb.size() == 1) break;
    }
    trim_q(q);
    trim_q(r);
    return { q, r };
}
static std::vector<mpq_class>
mod_q(const std::vector<mpq_class> &a, const std::vector<mpq_class> &b) {
    return divmod_q(a, b).second;
}
static std::vector<mpq_class>
gcd_q(std::vector<mpq_class> a, std::vector<mpq_class> b) {
    trim_q(a);
    trim_q(b);
    if (a.empty()) return b;
    if (b.empty()) return a;
    while (!b.empty()) {
        auto rem = mod_q(a, b);
        a.swap(b);
        b.swap(rem);
    }
    // make monic
    if (!a.empty()) {
        mpq_class lc = a.back();
        if (lc != 1) {
            for (auto &c : a) c /= lc;
        }
    }
    return a;
}
static bool
is_square_free_q(const std::vector<mpq_class> &f) {
    auto g = gcd_q(f, deriv_q(f));
    return g.size() <= 1; // gcd == 1 or 0
}
static std::vector<mpq_class>
square_free_part_q(const std::vector<mpq_class> &f) {
    auto g = gcd_q(f, deriv_q(f));
    if (g.size() <= 1) return f; // already square-free
    auto qr = divmod_q(f, g);
    return qr.first;
}

void
rational_poly_to_acb_poly(acb_poly_t result, const std::vector<mpq_class> &coeffs, slong prec) {
    acb_poly_zero(result);

    // Find degree (highest non-zero coefficient)
    int degree = -1;
    for (int i = coeffs.size() - 1; i >= 0; --i) {
        if (coeffs[i] != 0) {
            degree = i;
            break;
        }
    }

    if (degree < 0) {
        return; // Zero polynomial
    }

    // Set coefficients
    for (int i = 0; i <= degree; ++i) {
        if (coeffs[i] != 0) {
            acb_t coeff;
            acb_init(coeff);

            // Convert rational to acb
            // First set the real part
            arb_t real_part;
            arb_init(real_part);

            // Convert numerator and denominator
            fmpz_t num, den;
            fmpz_init(num);
            fmpz_init(den);

            fmpz_set_mpz(num, coeffs[i].get_num().get_mpz_t());
            fmpz_set_mpz(den, coeffs[i].get_den().get_mpz_t());

            // Set as exact rational
            arb_fmpz_div_fmpz(real_part, num, den, prec);

            // Set complex number with zero imaginary part
            acb_set_arb(coeff, real_part);

            // Set coefficient in polynomial
            acb_poly_set_coeff_acb(result, i, coeff);

            // Cleanup
            fmpz_clear(num);
            fmpz_clear(den);
            arb_clear(real_part);
            acb_clear(coeff);
        }
    }
}

bool
is_root_real_certified(const acb_t root, double epsilon) {
    // Check if the imaginary part contains zero
    if (arb_contains_zero(acb_imagref(root))) {
        // The root might be real - check the radius (not the magnitude)
        mag_t bound;
        mag_init(bound);
        mag_set(bound, arb_radref(acb_imagref(root)));
        double radius = mag_get_d(bound);
        mag_clear(bound);

        return radius < epsilon;
    }

    return false;
}

std::vector<std::complex<double>>
find_polynomial_roots_flint(const std::vector<mpq_class> &polynomial_coeffs, const FlintRootConfig &config) {
    std::vector<double> radii;
    return find_polynomial_roots_certified(polynomial_coeffs, radii, config);
}

std::vector<std::complex<double>>
find_polynomial_roots_certified(const std::vector<mpq_class> &polynomial_coeffs,
                                std::vector<double> &root_radii,
                                const FlintRootConfig &config) {
    std::vector<std::complex<double>> result;
    root_radii.clear();

    if (polynomial_coeffs.empty()) { return result; }

    // Find degree
    int degree = -1;
    for (int i = polynomial_coeffs.size() - 1; i >= 0; --i) {
        if (polynomial_coeffs[i] != 0) {
            degree = i;
            break;
        }
    }

    if (degree <= 0) {
        return result; // Constant or zero polynomial
    }

    // If requested by FLINT/ARB behavior, deflate to square-free part for robust root isolation
    const std::vector<mpq_class> &poly_in = polynomial_coeffs;
    std::vector<mpq_class> sqfree = square_free_part_q(poly_in);

    // Initialize FLINT polynomial
    acb_poly_t poly;
    acb_poly_init(poly);

    // Convert to FLINT format
    rational_poly_to_acb_poly(poly, sqfree, config.initial_prec);

    // Allocate space for roots
    acb_ptr roots = _acb_vec_init(degree);

    // Find roots with increasing precision if needed
    slong prec = config.initial_prec;
    slong found_roots = 0;

    while (prec <= config.max_prec) {
        // Clear previous roots
        _acb_vec_zero(roots, degree);

        // Find roots
        found_roots = acb_poly_find_roots(roots, poly, NULL, 0, prec);

        if (found_roots == degree) {
            // All roots found successfully
            break;
        }

        // Increase precision and try again
        prec *= 2;
        if (prec > config.max_prec) {
            std::cerr << "Warning: Only found " << found_roots << " out of " << degree << " roots" << std::endl;
            break;
        }

        // Recompute polynomial with higher precision
        rational_poly_to_acb_poly(poly, sqfree, prec);
    }

    // Convert roots to std::complex<double>
    for (slong i = 0; i < found_roots; ++i) {
        // Get midpoint
        arb_t real_mid, imag_mid;
        arb_init(real_mid);
        arb_init(imag_mid);

        arb_get_mid_arb(real_mid, acb_realref(roots + i));
        arb_get_mid_arb(imag_mid, acb_imagref(roots + i));

        double re = arf_get_d(arb_midref(real_mid), ARF_RND_NEAR);
        double im = arf_get_d(arb_midref(imag_mid), ARF_RND_NEAR);

        // Check if root should be considered real
        if (config.isolate_real_roots && is_root_real_certified(roots + i, config.epsilon)) {
            im = 0.0; // Force to real axis
        }

        result.push_back(std::complex<double>(re, im));

        // Get certified error radius: max(radius(real), radius(imag))
        mag_t err_bound;
        mag_t temp;
        mag_init(err_bound);
        mag_init(temp);

        mag_set(err_bound, arb_radref(acb_realref(roots + i)));
        mag_set(temp, arb_radref(acb_imagref(roots + i)));
        mag_max(err_bound, err_bound, temp);

        root_radii.push_back(mag_get_d(err_bound));

        mag_clear(err_bound);
        mag_clear(temp);
        arb_clear(real_mid);
        arb_clear(imag_mid);
    }

    // Cleanup
    _acb_vec_clear(roots, degree);
    acb_poly_clear(poly);

    return result;
}

#endif // HAVE_FLINT_ARB

} // namespace julia_rur