from typing import TYPE_CHECKING, ClassVar, Optional, TypeVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.Basis import BasisBase
from pairinteraction.backend._wrapped.Database import Database
from pairinteraction.backend._wrapped.ket.KetAtom import KetAtomBase, KetAtomDouble, KetAtomFloat
from pairinteraction.backend._wrapped.Parity import Parity, get_cpp_parity
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

Ket_t = TypeVar("Ket_t", bound=KetAtomBase)
UnionCPPBasisAtom = Union[
    _backend.BasisAtomFloat, _backend.BasisAtomComplexFloat, _backend.BasisAtomDouble, _backend.BasisAtomComplexDouble
]
UnionTypeCPPBasisAtomCreator = Union[
    type[_backend.BasisAtomCreatorFloat],
    type[_backend.BasisAtomCreatorComplexFloat],
    type[_backend.BasisAtomCreatorDouble],
    type[_backend.BasisAtomCreatorComplexDouble],
]


class BasisAtomBase(BasisBase[Ket_t]):
    _cpp: UnionCPPBasisAtom
    _cpp_creator: ClassVar[UnionTypeCPPBasisAtomCreator]

    def __init__(
        self,
        species: str,
        n: Optional[tuple[int, int]] = None,
        nu: Optional[tuple[float, float]] = None,
        l: Optional[tuple[float, float]] = None,
        s: Optional[tuple[float, float]] = None,
        j: Optional[tuple[float, float]] = None,
        f: Optional[tuple[float, float]] = None,
        m: Optional[tuple[float, float]] = None,
        energy: Union[tuple[float, float], tuple["PlainQuantity[float]", "PlainQuantity[float]"], None] = None,
        energy_unit: str = "pint",
        parity: Optional[Parity] = None,
        database: Optional[Database] = None,
        additional_kets: Optional[list[Ket_t]] = None,
    ) -> None:
        creator = self._cpp_creator()
        creator.set_species(species)
        if f is not None:
            creator.restrict_quantum_number_f(*f)
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if parity is not None:
            creator.restrict_parity(get_cpp_parity(parity))
        if n is not None:
            creator.restrict_quantum_number_n(*n)
        if nu is not None:
            creator.restrict_quantum_number_nu(*nu)
        if l is not None:
            creator.restrict_quantum_number_l(*l)
        if s is not None:
            creator.restrict_quantum_number_s(*s)
        if j is not None:
            creator.restrict_quantum_number_j(*j)
        if energy is not None:
            min_energy_au = QuantityScalar(energy[0], energy_unit).to_base("ENERGY")
            max_energy_au = QuantityScalar(energy[1], energy_unit).to_base("ENERGY")
            creator.restrict_energy(min_energy_au, max_energy_au)
        if database is None:
            database = Database.get_global_instance()
        if additional_kets is not None:
            for ket in additional_kets:
                creator.append_ket(ket._cpp)  # type: ignore [reportPrivateUsage]
        self._cpp = creator.create(database._cpp)  # type: ignore [reportPrivateUsage]


class BasisAtomFloat(BasisAtomBase[KetAtomFloat]):
    _cpp: _backend.BasisAtomFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorFloat
    _TypeKet = KetAtomFloat


class BasisAtomComplexFloat(BasisAtomBase[KetAtomFloat]):
    _cpp: _backend.BasisAtomComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorComplexFloat
    _TypeKet = KetAtomFloat


class BasisAtomDouble(BasisAtomBase[KetAtomDouble]):
    _cpp: _backend.BasisAtomDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorDouble
    _TypeKet = KetAtomDouble


class BasisAtomComplexDouble(BasisAtomBase[KetAtomDouble]):
    _cpp: _backend.BasisAtomComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorComplexDouble
    _TypeKet = KetAtomDouble
