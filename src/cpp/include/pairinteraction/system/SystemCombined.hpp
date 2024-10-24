#pragma once

#include "pairinteraction/system/System.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <array>
#include <memory>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class OperatorCombined;

template <typename Scalar>
class BasisCombined;

template <typename Real>
class KetCombined;

template <typename T>
class SystemCombined;

template <typename Scalar>
class SystemAtom;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
struct traits::CrtpTraits<SystemCombined<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetCombined<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
    using basis_t = BasisCombined<scalar_t>;
    using operator_t = OperatorCombined<scalar_t>;
};

template <typename Scalar>
class SystemCombined : public System<SystemCombined<Scalar>> {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using Type = SystemCombined<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using basis_t = typename traits::CrtpTraits<Type>::basis_t;
    using ket_t = typename traits::CrtpTraits<Type>::ket_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;

    SystemCombined(std::shared_ptr<const basis_t> basis);

    Type &set_distance(real_t distance);
    Type &set_distance_vector(const std::array<real_t, 3> &vector);

private:
    void construct_hamiltonian() const override;
    static Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    calculate_tensor_product(const std::shared_ptr<const basis_t> &basis,
                             const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix1,
                             const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix2);

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> green_function_dipole_dipole{3, 3};
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> green_function_dipole_quadrupole{3, 5};
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> green_function_quadrupole_dipole{5, 3};
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> green_function_quadrupole_quadrupole{5, 5};
};

extern template class SystemCombined<float>;
extern template class SystemCombined<double>;
extern template class SystemCombined<std::complex<float>>;
extern template class SystemCombined<std::complex<double>>;
} // namespace pairinteraction
