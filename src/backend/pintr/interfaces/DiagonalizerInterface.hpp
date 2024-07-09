#pragma once

#include "pintr/utils/eigen_assertion.hpp"
#include "pintr/utils/eigen_compat.hpp"
#include "pintr/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pintr {
template <typename Scalar>
struct EigenSystemH {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> eigenvectors;
    Eigen::VectorX<real_t> eigenvalues;
};

template <typename Scalar>
class DiagonalizerInterface {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    virtual ~DiagonalizerInterface() = default;
    virtual EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                      real_t min_eigenvalue, real_t max_eigenvalue,
                                      int precision) const = 0;
};
} // namespace pintr