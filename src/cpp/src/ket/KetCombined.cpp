#include "pairinteraction/ket/KetCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <limits>
#include <string>

namespace pairinteraction {
template <typename Scalar>
KetCombined<Scalar>::KetCombined(
    Private /*unused*/, std::initializer_list<size_t> atomic_indices,
    std::initializer_list<std::string> atomic_labels,
    std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases, real_t energy)
    : Ket<real_t>(energy, calculate_quantum_number_f(atomic_indices, atomic_bases),
                  calculate_quantum_number_m(atomic_indices, atomic_bases),
                  calculate_parity(atomic_indices, atomic_bases)),
      atomic_indices(atomic_indices), atomic_labels(atomic_labels), atomic_bases(atomic_bases) {
    if (atomic_indices.size() != atomic_labels.size() ||
        atomic_indices.size() != atomic_bases.size()) {
        throw std::invalid_argument(
            "The number of atomic indices, atomic labels, and atomic bases must be the same.");
    }
}

template <typename Scalar>
std::string KetCombined<Scalar>::get_label() const {
    std::string label;
    std::string separator;
    for (const auto &atomic_label : atomic_labels) {
        label += separator + atomic_label;
        separator = "; ";
    }
    return label;
}

template <typename Scalar>
std::shared_ptr<KetCombined<Scalar>>
KetCombined<Scalar>::get_ket_for_different_quantum_number_m(real_t /*new_quantum_number_m*/) const {
    // If we use symmetrized states so that the quantum_number_f is the total
    // angular quantum number, the quantum_number_m is the magnetic quantum number
    // corresponding to the total angular quantum number and we can implement this
    // method.
    throw std::runtime_error("Not implemented.");
}

template <typename Scalar>
typename KetCombined<Scalar>::real_t KetCombined<Scalar>::calculate_quantum_number_f(
    const std::vector<size_t> & /*indices*/,
    const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> & /*bases*/) {
    // Because this ket state is not symmetrized, the quantum_number_f is not well-defined.
    return std::numeric_limits<real_t>::max();
}

template <typename Scalar>
typename KetCombined<Scalar>::real_t KetCombined<Scalar>::calculate_quantum_number_m(
    const std::vector<size_t> &indices,
    const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> &bases) {
    for (const auto &basis : bases) {
        if (!basis->has_quantum_number_m()) {
            return std::numeric_limits<real_t>::max();
        }
    }
    real_t total_quantum_number_m = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
        total_quantum_number_m += bases[i]->get_quantum_number_m(indices[i]);
    }
    return total_quantum_number_m;
}

template <typename Scalar>
Parity KetCombined<Scalar>::calculate_parity(
    const std::vector<size_t> & /*indices*/,
    const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> & /*bases*/) {
    // Because this ket state is not symmetrized, the parity is not well-defined.
    return Parity::UNKNOWN;
}

template <typename Scalar>
bool KetCombined<Scalar>::operator==(const KetCombined<Scalar> &other) const {
    return Ket<real_t>::operator==(other) && atomic_indices == other.atomic_indices &&
        atomic_bases == other.atomic_bases;
}

template <typename Scalar>
size_t KetCombined<Scalar>::hash::operator()(const KetCombined<Scalar> &k) const {
    size_t seed = typename Ket<real_t>::hash()(k);
    for (const auto &index : k.atomic_indices) {
        utils::hash_combine(seed, index);
    }
    for (const auto &basis : k.atomic_bases) {
        utils::hash_combine(seed, reinterpret_cast<std::uintptr_t>(basis.get()));
    }
    return seed;
}

// Explicit instantiations
template class KetCombined<float>;
template class KetCombined<double>;
template class KetCombined<std::complex<float>>;
template class KetCombined<std::complex<double>>;
} // namespace pairinteraction
