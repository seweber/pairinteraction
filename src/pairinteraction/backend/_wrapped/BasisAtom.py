from functools import cached_property
from typing import Optional, Union

import scipy.sparse

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.KetAtom import KetAtomDouble, KetAtomFloat


class BasisAtomBase:
    _BasisAtomCreator: Union[
        type[_backend.BasisAtomCreatorFloat],
        type[_backend.BasisAtomCreatorComplexFloat],
        type[_backend.BasisAtomCreatorDouble],
        type[_backend.BasisAtomCreatorComplexDouble],
    ]
    _KetAtom: Union[type[KetAtomFloat], type[KetAtomDouble]]

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
        energy: Optional[tuple[float, float]] = None,
        parity: Optional[_backend.Parity] = None,
        database: Optional[_backend.Database] = None,
    ) -> None:
        self._cpp = self._create_obj(species, n, nu, l, s, j, f, m, energy, parity, database)

    @classmethod
    def _create_obj(
        cls,
        species: str,
        n: Optional[tuple[int, int]] = None,
        nu: Optional[tuple[float, float]] = None,
        l: Optional[tuple[float, float]] = None,
        s: Optional[tuple[float, float]] = None,
        j: Optional[tuple[float, float]] = None,
        f: Optional[tuple[float, float]] = None,
        m: Optional[tuple[float, float]] = None,
        energy: Optional[tuple[float, float]] = None,
        parity: Optional[_backend.Parity] = None,
        database: Optional[_backend.Database] = None,
    ):
        creator = cls._BasisAtomCreator()
        creator.set_species(species)
        if energy is not None:
            creator.restrict_energy(*energy)
        if f is not None:
            creator.restrict_quantum_number_f(*f)
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if parity is not None:
            creator.restrict_parity(parity)
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
        if database is None:
            database = _backend.Database.get_global_instance(
                download_missing=False, wigner_in_memory=True, database_dir=""
            )  # Use the default database
        cpp_obj = creator.create(database)
        return cpp_obj

    @cached_property
    def kets(self) -> list[Union[KetAtomFloat, KetAtomDouble]]:
        kets = [self._KetAtom._from_cpp_object(ket) for ket in self._cpp.get_kets()]
        return kets

    @property
    def number_of_states(self) -> int:
        return self._cpp.get_number_of_states()

    @property
    def number_of_kets(self) -> int:
        return self._cpp.get_number_of_kets()

    @property
    def coefficients(self) -> scipy.sparse.csr_matrix:
        return self._cpp.get_coefficients()


class BasisAtomFloat(BasisAtomBase):
    _cpp: _backend.BasisAtomFloat
    _BasisAtomCreator = _backend.BasisAtomCreatorFloat
    _KetAtom = KetAtomFloat


class BasisAtomComplexFloat(BasisAtomBase):
    _cpp: _backend.BasisAtomComplexFloat
    _BasisAtomCreator = _backend.BasisAtomCreatorComplexFloat
    _KetAtom = KetAtomFloat


class BasisAtomDouble(BasisAtomBase):
    _cpp: _backend.BasisAtomDouble
    _BasisAtomCreator = _backend.BasisAtomCreatorDouble
    _KetAtom = KetAtomDouble


class BasisAtomComplexDouble(BasisAtomBase):
    _cpp: _backend.BasisAtomComplexDouble
    _BasisAtomCreator = _backend.BasisAtomCreatorComplexDouble
    _KetAtom = KetAtomDouble
