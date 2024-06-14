#pragma once

#include "utils/eigen_assertion.hpp"
#include "utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <array>
#include <type_traits>

namespace euler {

/**
 * @function get_rotation_matrix
 *
 * @brief Build a matrix that rotates the coordinate system to the new z-axis and y-axis.
 *
 * @tparam Real The type of a real number.
 *
 * @param to_z_axis The new z-axis.
 *
 * @param to_y_axis The new y-axis.
 *
 * @return The passive rotation matrix.
 */

template <typename Real>
inline Eigen::Matrix3<Real> get_rotation_matrix(std::array<Real, 3> to_z_axis,
                                                std::array<Real, 3> to_y_axis) {
    static_assert(std::is_floating_point_v<Real>);

    auto to_z_axis_mapped = Eigen::Map<Eigen::Vector3<Real>>(to_z_axis.data()).normalized();
    auto to_y_axis_mapped = Eigen::Map<Eigen::Vector3<Real>>(to_y_axis.data()).normalized();

    if (std::abs(to_z_axis_mapped.dot(to_y_axis_mapped)) > std::numeric_limits<Real>::epsilon()) {
        throw std::runtime_error("The z-axis and the y-axis are not orhogonal.");
    }

    Eigen::Matrix3<Real> rotator;
    rotator << to_y_axis_mapped.cross(to_z_axis_mapped), to_y_axis_mapped, to_z_axis_mapped;

    return rotator;
}

/**
 * @function get_euler_angles
 *
 * @brief Extract the Euler angles alpha, beta, gamma
 *
 * @tparam Real The type of a real number.
 *
 * @param to_z_axis The new z-axis.
 *
 * @param to_y_axis The new y-axis.
 *
 * @return The Euler angles.
 */

template <typename Real>
inline std::array<Real, 3> get_euler_angles(std::array<Real, 3> to_z_axis,
                                            std::array<Real, 3> to_y_axis) {
    static_assert(std::is_floating_point_v<Real>);

    auto rotator = get_rotation_matrix(to_z_axis, to_y_axis);
    std::array<Real, 3> euler_zyz;
    Eigen::Map<Eigen::Vector3<Real>>(euler_zyz.data()) = rotator.eulerAngles(2, 1, 2);
    return euler_zyz;
}

} // namespace euler

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

DOCTEST_TEST_CASE("construction of rotation matrixes") {
    auto rotator = euler::get_rotation_matrix<double>({0, 0, 1}, {0, 1, 0});
    auto rotator_reference = Eigen::Matrix3<double>::Identity();
    DOCTEST_CHECK((rotator - rotator_reference).norm() == 0);

    rotator = euler::get_rotation_matrix<double>({0, 0, 1}, {1, 1, 0});
    auto y_axis = Eigen::Vector3<double>{0, 1, 0};
    auto rotated_y_axis = rotator * y_axis;
    auto rotated_y_axis_reference = Eigen::Vector3<double>{1, 1, 0}.normalized();
    DOCTEST_CHECK((rotated_y_axis - rotated_y_axis_reference).norm() == 0);

    rotator = euler::get_rotation_matrix<double>({1, 0, 0}, {0, 1, 0});
    auto z_axis = Eigen::Vector3<double>{0, 0, 1};
    auto rotated_z_axis = rotator * z_axis;
    auto rotated_z_axis_reference = Eigen::Vector3<double>{1, 0, 0};
    DOCTEST_CHECK((rotated_z_axis - rotated_z_axis_reference).norm() == 0);

    std::string error_msg = "The z-axis and the y-axis are not orhogonal.";
    DOCTEST_CHECK_THROWS_WITH_AS(euler::get_rotation_matrix<double>({0, 0, 1}, {0, 1, 1});
                                 , error_msg.c_str(), std::runtime_error);
}

DOCTEST_TEST_CASE("construction of zyz euler angles") {
    constexpr double PI = 3.141592653589793238462643383279502884;

    auto euler_angles = euler::get_euler_angles<double>({1, 1, 0}, {0, 0, 1});
    std::array<double, 3> euler_angles_reference{0.25 * PI, 0.5 * PI, 0.5 * PI};
    DOCTEST_CHECK(std::abs(euler_angles[0] - euler_angles_reference[0]) +
                      std::abs(euler_angles[1] - euler_angles_reference[1]) +
                      std::abs(euler_angles[2] - euler_angles_reference[2]) <=
                  1e-6);
}
