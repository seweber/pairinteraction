#include "pairinteraction/convenience/matrix_elements.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>

namespace pairinteraction {

template <typename Scalar>
typename traits::NumTraits<Scalar>::real_t
calculate_energy(std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> ket,
                 const SystemAtom<Scalar> &system) {
    if (!system.is_diagonal()) {
        throw std::invalid_argument("The system must be diagonalized.");
    }

    size_t state_index = system.get_basis()->get_state_index_with_largest_overlap(ket);

    if (system.get_basis()->get_overlaps(ket)[state_index] < 0.5) {
        throw std::invalid_argument("There is no eigenstate that corresponds clearly to the ket.");
    }

    return std::real(system.get_matrix().coeff(state_index, state_index));
}

template <typename Real>
Real calculate_energy(std::shared_ptr<const KetAtom<Real>> ket) {
    return ket->get_energy();
}

template <typename Scalar>
Scalar calculate_electric_dipole_matrix_element(
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> initial_ket,
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> final_ket,
    const SystemAtom<Scalar> &system, int q) {
    if (!system.is_diagonal()) {
        throw std::invalid_argument("The system must be diagonalized.");
    }

    auto initial_state = system.get_basis()->get_state_with_largest_overlap(initial_ket);
    auto final_state = system.get_basis()->get_state_with_largest_overlap(final_ket);

    if (initial_state->get_overlaps(initial_ket)[0] < 0.5) {
        throw std::invalid_argument(
            "There is no eigenstate that corresponds clearly to the initial ket.");
    }
    if (final_state->get_overlaps(final_ket)[0] < 0.5) {
        throw std::invalid_argument(
            "There is no eigenstate that corresponds clearly to the final ket.");
    }

    auto matrix_element = initial_ket->get_database().get_matrix_elements(
        initial_state, final_state, OperatorType::ELECTRIC_DIPOLE, q);

    assert(matrix_element.rows() == 1);
    assert(matrix_element.cols() == 1);

    return matrix_element.coeff(0, 0);
}

template <typename Real>
Real calculate_electric_dipole_matrix_element(std::shared_ptr<const KetAtom<Real>> initial_ket,
                                              std::shared_ptr<const KetAtom<Real>> final_ket,
                                              int q) {
    auto basis = BasisAtomCreator<Real>()
                     .append_ket(initial_ket)
                     .append_ket(final_ket)
                     .create(initial_ket->get_database());
    SystemAtom<Real> system(basis);

    return calculate_electric_dipole_matrix_element(initial_ket, final_ket, system, q);
}

// Explicit instantiations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define INSTANTIATE_CALCULATE_FUNCTIONS_REAL(REAL)                                                 \
    template REAL calculate_energy<REAL>(std::shared_ptr<const KetAtom<REAL>>);                    \
    template REAL calculate_electric_dipole_matrix_element<REAL>(                                  \
        std::shared_ptr<const KetAtom<REAL>>, std::shared_ptr<const KetAtom<REAL>>, int);
#define INSTANTIATE_CALCULATE_FUNCTIONS(SCALAR, REAL)                                              \
    template REAL calculate_energy<SCALAR>(std::shared_ptr<const KetAtom<REAL>>,                   \
                                           const SystemAtom<SCALAR> &);                            \
    template SCALAR calculate_electric_dipole_matrix_element<SCALAR>(                              \
        std::shared_ptr<const KetAtom<REAL>>, std::shared_ptr<const KetAtom<REAL>>,                \
        const SystemAtom<SCALAR> &, int);
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

INSTANTIATE_CALCULATE_FUNCTIONS_REAL(float)
INSTANTIATE_CALCULATE_FUNCTIONS_REAL(double)

INSTANTIATE_CALCULATE_FUNCTIONS(float, float)
INSTANTIATE_CALCULATE_FUNCTIONS(double, double)
INSTANTIATE_CALCULATE_FUNCTIONS(std::complex<float>, float)
INSTANTIATE_CALCULATE_FUNCTIONS(std::complex<double>, double)

#undef INSTANTIATE_CALCULATE_FUNCTIONS_REAL
#undef INSTANTIATE_CALCULATE_FUNCTIONS
} // namespace pairinteraction
