#pragma once

#include "pairinteraction/basis/Basis.hpp"
#include "pairinteraction/ket/KetClassicalLight.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>

namespace pairinteraction {
template <typename Scalar>
class BasisClassicalLightCreator;

template <typename Scalar>
class BasisClassicalLight;

template <typename Scalar>
struct traits::CrtpTraits<BasisClassicalLight<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetClassicalLight<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

/**
 * @class BasisAtom
 *
 * @brief Class for creating a basis of atomic kets.
 *
 * @tparam T Complex number type.
 */
template <typename Scalar>
class BasisClassicalLight : public Basis<BasisClassicalLight<Scalar>> {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    friend class BasisClassicalLightCreator<Scalar>;
    struct Private {};

public:
    using Type = BasisClassicalLight<Scalar>;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;

    BasisClassicalLight(Private /*unused*/, ketvec_t &&kets, std::string &&id_of_kets);
};

extern template class BasisClassicalLight<float>;
extern template class BasisClassicalLight<double>;
extern template class BasisClassicalLight<std::complex<float>>;
extern template class BasisClassicalLight<std::complex<double>>;
} // namespace pairinteraction
