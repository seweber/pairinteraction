#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class BasisCombinedCreator;

template <typename Real>
class KetCombined;

template <typename Scalar>
class BasisCombined;

template <typename Scalar>
class BasisAtom;

template <typename Real>
class KetAtom;

template <typename Scalar>
struct traits::CrtpTraits<BasisCombined<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetCombined<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

template <typename Scalar>
class BasisCombined : public Basis<BasisCombined<Scalar>> {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    friend class BasisCombinedCreator<Scalar>;
    struct Private {};

public:
    using Type = BasisCombined<Scalar>;
    using real_t = typename traits::CrtpTraits<Type>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;
    using range_t = Range<size_t>;
    using map_size_t = std::unordered_map<size_t, size_t>;
    using map_range_t = std::unordered_map<size_t, range_t>;

    BasisCombined(Private /*unused*/, ketvec_t &&kets, std::string &&id_of_kets,
                  map_range_t &&map_range_of_state_index2,
                  std::shared_ptr<const BasisAtom<Scalar>> basis1,
                  std::shared_ptr<const BasisAtom<Scalar>> basis2);
    const range_t &get_index_range(size_t state_index1) const;
    std::shared_ptr<const BasisAtom<Scalar>> get_basis1() const;
    std::shared_ptr<const BasisAtom<Scalar>> get_basis2() const;

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_amplitudes(std::shared_ptr<const KetAtom<real_t>> ket1,
                   std::shared_ptr<const KetAtom<real_t>> ket2) const;
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
    get_amplitudes(std::shared_ptr<const BasisAtom<Scalar>> other1,
                   std::shared_ptr<const BasisAtom<Scalar>> other2) const;
    Eigen::SparseMatrix<real_t, Eigen::RowMajor>
    get_overlaps(std::shared_ptr<const KetAtom<real_t>> ket1,
                 std::shared_ptr<const KetAtom<real_t>> ket2) const;
    Eigen::SparseMatrix<real_t, Eigen::RowMajor>
    get_overlaps(std::shared_ptr<const BasisAtom<Scalar>> other1,
                 std::shared_ptr<const BasisAtom<Scalar>> other2) const;

private:
    map_range_t map_range_of_state_index2;
    std::shared_ptr<const BasisAtom<Scalar>> basis1;
    std::shared_ptr<const BasisAtom<Scalar>> basis2;
};

extern template class BasisCombined<float>;
extern template class BasisCombined<double>;
extern template class BasisCombined<std::complex<float>>;
extern template class BasisCombined<std::complex<double>>;
} // namespace pairinteraction
