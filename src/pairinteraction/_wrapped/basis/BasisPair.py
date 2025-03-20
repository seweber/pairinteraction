from collections.abc import Collection, Sequence
from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union, overload

import numpy as np
from typing_extensions import TypeGuard

from pairinteraction import _backend
from pairinteraction._wrapped.basis.Basis import BasisBase
from pairinteraction._wrapped.basis.BasisAtom import BasisAtomComplex, BasisAtomReal
from pairinteraction._wrapped.cpp_types import OperatorType, Parity, get_cpp_operator_type, get_cpp_parity
from pairinteraction._wrapped.ket.KetPair import (
    KetPair,
    KetPairComplex,
    KetPairReal,
    is_ket_atom_tuple,
    is_ket_pair_like,
)
from pairinteraction._wrapped.system.SystemAtom import SystemAtomComplex, SystemAtomReal
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix

    from pairinteraction._wrapped.basis.BasisAtom import BasisAtom
    from pairinteraction._wrapped.ket.KetAtom import KetAtom  # noqa: F401  # required fo sphinx for KetPairLike
    from pairinteraction._wrapped.ket.KetPair import KetPairLike
    from pairinteraction._wrapped.system.SystemAtom import SystemAtom
    from pairinteraction.units import NDArray, PintArray, PintFloat, PintSparse

KetPairType = TypeVar("KetPairType", bound=KetPair, covariant=True)
BasisPairLike = Union[
    "BasisPairReal",
    "BasisPairComplex",
    tuple["BasisAtomReal", "BasisAtomReal"],
    Sequence["BasisAtomReal"],
    tuple["BasisAtomComplex", "BasisAtomComplex"],
    Sequence["BasisAtomComplex"],
]

UnionCPPBasisPair = Union[_backend.BasisPairReal, _backend.BasisPairComplex]
UnionTypeCPPBasisPairCreator = Union[type[_backend.BasisPairCreatorReal], type[_backend.BasisPairCreatorComplex]]


def is_basis_pair_like(obj: Any) -> TypeGuard[BasisPairLike]:
    return isinstance(obj, BasisPair) or (len(obj) == 2 and all(isinstance(x, BasisAtom) for x in obj))


def is_basis_atom_real_tuple(obj: Any) -> TypeGuard[tuple["BasisAtomReal", "BasisAtomReal"]]:
    return len(obj) == 2 and all(isinstance(x, BasisAtomReal) for x in obj)


def is_basis_atom_complex_tuple(obj: Any) -> TypeGuard[tuple["BasisAtomComplex", "BasisAtomComplex"]]:
    return len(obj) == 2 and all(isinstance(x, BasisAtomComplex) for x in obj)


class BasisPair(BasisBase[KetPairType]):
    """Basis for a pair of atoms.

    Add all product states of the eigenstates of two given SystemAtom objects to the basis,
    which pair energy is within the given energy range.
    You can also specify which total magnetic quantum number m the pair should have (if it is conserved)
    and the product of the parities of the two atoms.
    Due to the possible restrictions of the basis states, the BasisPair coefficients matrix will in general
    not be square but (n x d),
    where n is the number of all involved kets (typically basis1.number_of_kets * basis2.number_of_kets)
    and d is the number of basis states (after applying the restrictions).

    Examples:
        >>> import pairinteraction.real as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
        >>> system = pi.SystemAtom(basis).set_magnetic_field([0, 0, 1], unit="G").diagonalize()
        >>> pair_energy = 2 * system.get_corresponding_energy(ket, unit="GHz")
        >>> pair_basis = pi.BasisPair(
        ...     [system, system],
        ...     energy=(pair_energy - 3, pair_energy + 3),
        ...     energy_unit="GHz",
        ... )
        >>> print(pair_basis)
        BasisPair(|Rb:59,S_1/2,-1/2; Rb:61,S_1/2,-1/2⟩ ... |Rb:58,F_7/2,7/2; Rb:59,S_1/2,1/2⟩)

    """

    _cpp: UnionCPPBasisPair
    _cpp_creator: ClassVar[UnionTypeCPPBasisPairCreator]

    def __init__(
        self,
        systems: Collection["SystemAtom[Any]"],
        m: Optional[tuple[float, float]] = None,
        product_of_parities: Optional[Parity] = None,
        energy: Union[tuple[float, float], tuple["PintFloat", "PintFloat"], None] = None,
        energy_unit: Optional[str] = None,
    ) -> None:
        """Create a basis for a pair of atoms.

        Args:
            systems: tuple of two SystemAtom objects, which define the two atoms, from which the BasisPair is build.
                Both systems have to be diagonalized before creating the BasisPair.
            m: tuple of (min, max) values for the total magnetic quantum number m of the pair state.
                Default None, i.e. no restriction.
            product_of_parities: The product parity of the states to consider.
                Default None, i.e. add all available states.
            energy: tuple of (min, max) value for the pair energy. Default None, i.e. add all available states.
            energy_unit: In which unit the energy values are given, e.g. "GHz".
                Default None, i.e. energy is provided as pint object.

        """
        creator = self._cpp_creator()
        for system in systems:
            if isinstance(creator, _backend.BasisPairCreatorReal) and isinstance(system, SystemAtomReal):
                creator.add(system._cpp)
            elif isinstance(creator, _backend.BasisPairCreatorComplex) and isinstance(system, SystemAtomComplex):
                creator.add(system._cpp)
            else:
                raise ValueError(f"Incompatible types: {type(system)=}; {type(creator)=}")
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if product_of_parities is not None:
            creator.restrict_product_of_parities(get_cpp_parity(product_of_parities))
        if energy is not None:
            min_energy_au = QuantityScalar.from_pint_or_unit(energy[0], energy_unit, "energy").to_base_unit()
            max_energy_au = QuantityScalar.from_pint_or_unit(energy[1], energy_unit, "energy").to_base_unit()
            creator.restrict_energy(min_energy_au, max_energy_au)
        self._cpp = creator.create()

    @overload
    def get_amplitudes(self, ket_or_basis: "KetPairLike") -> "NDArray": ...

    @overload
    def get_amplitudes(self, ket_or_basis: BasisPairLike) -> "csr_matrix": ...

    def get_amplitudes(self, ket_or_basis: Union["KetPairLike", BasisPairLike]) -> Union["NDArray", "csr_matrix"]:
        if is_ket_pair_like(ket_or_basis):
            if is_ket_atom_tuple(ket_or_basis):
                ket_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                return np.array(self._cpp.get_amplitudes(*ket_cpp))
            elif isinstance(ket_or_basis, KetPairReal) and isinstance(self._cpp, _backend.BasisPairReal):
                return np.array(self._cpp.get_amplitudes(ket_or_basis._cpp))  # type: ignore [call-overload] # _backend.pyi is incorrect
            elif isinstance(ket_or_basis, KetPairComplex) and isinstance(self._cpp, _backend.BasisPairComplex):
                return np.array(self._cpp.get_amplitudes(ket_or_basis._cpp))  # type: ignore [call-overload] # _backend.pyi is incorrect
            else:
                raise ValueError(f"Incompatible types: {type(ket_or_basis)=}; {type(self)=}")
        elif is_basis_pair_like(ket_or_basis):
            if is_basis_atom_real_tuple(ket_or_basis) and isinstance(self._cpp, _backend.BasisPairReal):
                basis_real_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                return self._cpp.get_amplitudes(*basis_real_cpp)
            elif is_basis_atom_complex_tuple(ket_or_basis) and isinstance(self._cpp, _backend.BasisPairComplex):
                basis_complex_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                return self._cpp.get_amplitudes(*basis_complex_cpp)
            elif isinstance(ket_or_basis, BasisPairReal) and isinstance(self._cpp, _backend.BasisPairReal):
                return self._cpp.get_amplitudes(ket_or_basis._cpp)  # type: ignore [call-overload,no-any-return] # _backend.pyi is incorrect
            elif isinstance(ket_or_basis, BasisPairComplex) and isinstance(self._cpp, _backend.BasisPairComplex):
                return self._cpp.get_amplitudes(ket_or_basis._cpp)  # type: ignore [call-overload,no-any-return] # _backend.pyi is incorrect
            else:
                raise ValueError(f"Incompatible types: {type(ket_or_basis)=}; {type(self)=}")
        else:
            raise ValueError(f"Unknown type: {type(ket_or_basis)=}")

    @overload
    def get_overlaps(self, ket_or_basis: "KetPairLike") -> "NDArray": ...

    @overload
    def get_overlaps(self, ket_or_basis: BasisPairLike) -> "csr_matrix": ...

    def get_overlaps(self, ket_or_basis: Union["KetPairLike", BasisPairLike]) -> Union["NDArray", "csr_matrix"]:
        if is_ket_pair_like(ket_or_basis):
            if is_ket_atom_tuple(ket_or_basis):
                ket_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                return np.array(self._cpp.get_overlaps(*ket_cpp))
            elif isinstance(ket_or_basis, KetPairReal) and isinstance(self._cpp, _backend.BasisPairReal):
                return np.array(self._cpp.get_overlaps(ket_or_basis._cpp))  # type: ignore [call-overload] # _backend.pyi is incorrect
            elif isinstance(ket_or_basis, KetPairComplex) and isinstance(self._cpp, _backend.BasisPairComplex):
                return np.array(self._cpp.get_overlaps(ket_or_basis._cpp))  # type: ignore [call-overload] # _backend.pyi is incorrect
            else:
                raise ValueError(f"Incompatible types: {type(ket_or_basis)=}; {type(self)=}")
        elif is_basis_pair_like(ket_or_basis):
            if is_basis_atom_real_tuple(ket_or_basis) and isinstance(self._cpp, _backend.BasisPairReal):
                basis_real_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                return self._cpp.get_overlaps(*basis_real_cpp)
            elif is_basis_atom_complex_tuple(ket_or_basis) and isinstance(self._cpp, _backend.BasisPairComplex):
                basis_complex_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                return self._cpp.get_overlaps(*basis_complex_cpp)
            elif isinstance(ket_or_basis, BasisPairReal) and isinstance(self._cpp, _backend.BasisPairReal):
                return self._cpp.get_overlaps(ket_or_basis._cpp)  # type: ignore [call-overload,no-any-return] # _backend.pyi is incorrect
            elif isinstance(ket_or_basis, BasisPairComplex) and isinstance(self._cpp, _backend.BasisPairComplex):
                return self._cpp.get_overlaps(ket_or_basis._cpp)  # type: ignore [call-overload,no-any-return] # _backend.pyi is incorrect
            else:
                raise ValueError(f"Incompatible types: {type(ket_or_basis)=}; {type(self)=}")
        else:
            raise ValueError(f"Unknown type: {type(ket_or_basis)=}")

    @overload
    def get_matrix_elements(
        self,
        ket_or_basis: "KetPairLike",
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: None = None,
    ) -> "PintArray": ...

    @overload
    def get_matrix_elements(
        self, ket_or_basis: "KetPairLike", operators: tuple[OperatorType, OperatorType], qs: tuple[int, int], unit: str
    ) -> "NDArray": ...

    @overload
    def get_matrix_elements(
        self,
        ket_or_basis: BasisPairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: None = None,
    ) -> "PintSparse": ...  # type: ignore [type-var] # see PintSparse

    @overload
    def get_matrix_elements(
        self,
        ket_or_basis: BasisPairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: str,
    ) -> "csr_matrix": ...

    def get_matrix_elements(
        self,
        ket_or_basis: Union["KetPairLike", BasisPairLike],
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: Optional[str] = None,
    ) -> Union["NDArray", "PintArray", "csr_matrix", "PintSparse"]:
        operators_cpp = (get_cpp_operator_type(operators[0]), get_cpp_operator_type(operators[1]))

        matrix_elements_au: Union[NDArray, csr_matrix]
        matrix_elements: Union[QuantityArray, QuantitySparse]
        if is_ket_pair_like(ket_or_basis):
            if is_ket_atom_tuple(ket_or_basis):
                ket_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                matrix_elements_au = np.array(self._cpp.get_matrix_elements(*ket_cpp, *operators_cpp, *qs))
            elif isinstance(ket_or_basis, KetPairReal) and isinstance(self._cpp, _backend.BasisPairReal):
                matrix_elements_au = np.array(self._cpp.get_matrix_elements(ket_or_basis._cpp, *operators_cpp, *qs))
            elif isinstance(ket_or_basis, KetPairComplex) and isinstance(self._cpp, _backend.BasisPairComplex):
                matrix_elements_au = np.array(self._cpp.get_matrix_elements(ket_or_basis._cpp, *operators_cpp, *qs))
            else:
                raise ValueError(f"Incompatible types: {type(ket_or_basis)=}; {type(self)=}")
            matrix_elements = QuantityArray.from_base_unit(matrix_elements_au, operators)
        elif is_basis_pair_like(ket_or_basis):
            if is_basis_atom_real_tuple(ket_or_basis) and isinstance(self._cpp, _backend.BasisPairReal):
                basis_real_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                matrix_elements_au = self._cpp.get_matrix_elements(*basis_real_cpp, *operators_cpp, *qs)
            elif is_basis_atom_complex_tuple(ket_or_basis) and isinstance(self._cpp, _backend.BasisPairComplex):
                basis_complex_cpp = (ket_or_basis[0]._cpp, ket_or_basis[1]._cpp)
                matrix_elements_au = self._cpp.get_matrix_elements(*basis_complex_cpp, *operators_cpp, *qs)
            elif isinstance(ket_or_basis, BasisPairReal) and isinstance(self._cpp, _backend.BasisPairReal):
                matrix_elements_au = self._cpp.get_matrix_elements(ket_or_basis._cpp, *operators_cpp, *qs)
            elif isinstance(ket_or_basis, BasisPairComplex) and isinstance(self._cpp, _backend.BasisPairComplex):
                matrix_elements_au = self._cpp.get_matrix_elements(ket_or_basis._cpp, *operators_cpp, *qs)
            else:
                raise ValueError(f"Incompatible types: {type(ket_or_basis)=}; {type(self)=}")
            matrix_elements = QuantitySparse.from_base_unit(matrix_elements_au, operators)
        else:
            raise ValueError(f"Unknown type: {type(ket_or_basis)=}")

        return matrix_elements.to_pint_or_unit(unit)


class BasisPairReal(BasisPair[KetPairReal]):
    _cpp: _backend.BasisPairReal
    _cpp_creator = _backend.BasisPairCreatorReal
    _TypeKet = KetPairReal


class BasisPairComplex(BasisPair[KetPairComplex]):
    _cpp: _backend.BasisPairComplex
    _cpp_creator = _backend.BasisPairCreatorComplex
    _TypeKet = KetPairComplex
