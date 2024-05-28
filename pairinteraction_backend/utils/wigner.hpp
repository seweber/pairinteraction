#pragma once

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include "utils/maths.hpp"
#include "utils/traits.hpp"

namespace wigner {

namespace {

template <typename Real>
inline constexpr Real PI = 3.141592653589793238462643383279502884;

template <typename Real>
inline Real wigner_uppercase_d_matrix_pi_half(Real f, Real m_initial, Real m_final) {
    static_assert(std::is_floating_point_v<Real>);

    Real result = 0;
    for (Real k = std::max(0, static_cast<int>(m_final - m_initial));
         k <= f + std::min(m_final, -m_initial); ++k) {
        result += std::pow(-1, k) * maths::binomial_coefficient(f + m_final, k) *
            maths::binomial_coefficient(f - m_final, k + m_initial - m_final);
    }
    result *= std::pow(-1, m_initial - m_final) / std::pow(2, f) *
        std::sqrt(maths::factorial(f + m_initial) * maths::factorial(f - m_initial) /
                  (maths::factorial(f + m_final) * maths::factorial(f - m_final)));
    return result;
}

} // namespace

template <typename Scalar>
inline Scalar wigner_uppercase_d_matrix(typename traits::NumTraits<Scalar>::real_t f,
                                        typename traits::NumTraits<Scalar>::real_t m_initial,
                                        typename traits::NumTraits<Scalar>::real_t m_final,
                                        typename traits::NumTraits<Scalar>::real_t alpha,
                                        typename traits::NumTraits<Scalar>::real_t beta,
                                        typename traits::NumTraits<Scalar>::real_t gamma) {
    static_assert(traits::is_complex_or_floating_point_v<Scalar>);

    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        return Scalar(std::cos(-m_initial * alpha), std::sin(-m_initial * alpha)) *
            wigner_uppercase_d_matrix<typename traits::NumTraits<Scalar>::real_t>(
                   f, m_initial, m_final, 0, beta, 0) *
            Scalar(std::cos(-m_final * gamma), std::sin(-m_final * gamma));
    } else {
        if (std::abs(
                std::remainder(m_initial * alpha, PI<typename traits::NumTraits<Scalar>::real_t>)) >
            10 * std::numeric_limits<typename traits::NumTraits<Scalar>::real_t>::epsilon()) {
            throw std::invalid_argument(
                "The scalar type must be complex if m_initial*alpha is not a multiple of pi");
        }
        if (std::abs(
                std::remainder(m_final * gamma, PI<typename traits::NumTraits<Scalar>::real_t>)) >
            10 * std::numeric_limits<typename traits::NumTraits<Scalar>::real_t>::epsilon()) {
            throw std::invalid_argument(
                "The scalar type must be complex if m_final*gamma is not a multiple of pi");
        }
        std::complex<Scalar> result = 0;
        for (typename traits::NumTraits<Scalar>::real_t m = -f; m <= f; ++m) {
            result += wigner_uppercase_d_matrix_pi_half(f, m_initial, m) *
                std::complex<typename traits::NumTraits<Scalar>::real_t>(std::cos(-m * beta),
                                                                         std::sin(-m * beta)) *
                wigner_uppercase_d_matrix_pi_half(f, m, -m_final);
        }
        result *= std::pow(std::complex<typename traits::NumTraits<Scalar>::real_t>(0, 1),
                           2 * f - m_initial - m_final) *
            static_cast<typename traits::NumTraits<Scalar>::real_t>(std::pow(-1, 2 * m_initial));
        return result.real();
    }
}

} // namespace wigner

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

DOCTEST_TEST_CASE("construction of wigner D matrix") {
    constexpr double PI = 3.141592653589793238462643383279502884;
    constexpr double EPSILON = 100 * std::numeric_limits<double>::epsilon();

    auto wigner_real_entry =
        wigner::wigner_uppercase_d_matrix<double>(0.5, 0.5, -0.5, 4 * PI, PI / 3, 2 * PI);
    auto wigner_real_entry_reference = -0.5;
    DOCTEST_CHECK((wigner_real_entry - wigner_real_entry_reference) <= EPSILON);

    std::string error_msg =
        "The scalar type must be complex if m_initial*alpha is not a multiple of pi";
    DOCTEST_CHECK_THROWS_WITH_AS(
        wigner::wigner_uppercase_d_matrix<double>(0.5, 0.5, -0.5, 0.1 * PI, 0, 0);
        , error_msg.c_str(), std::invalid_argument);

    error_msg = "The scalar type must be complex if "
                "m_final*gamma is not a multiple of pi";
    DOCTEST_CHECK_THROWS_WITH_AS(
        wigner::wigner_uppercase_d_matrix<double>(0.5, 0.5, -0.5, 0, 0, 0.1 * PI);
        , error_msg.c_str(), std::invalid_argument);

    auto wigner_complex_entry = wigner::wigner_uppercase_d_matrix<std::complex<double>>(
        0.5, 0.5, -0.5, 0.5 * PI, PI, -0.5 * PI);
    auto wigner_complex_entry_reference = std::complex<double>(0, 1);
    DOCTEST_CHECK(std::abs(wigner_complex_entry - wigner_complex_entry_reference) <= EPSILON);
}
