"""Test the conversion from BaseUnits (=atomic units) used in the backend and the units input and output to the user."""

import numpy as np

import pairinteraction.backend.real as pi
from pairinteraction.units import BaseUnits, QuantityScalar, ureg


def test_magnetic() -> None:
    """Test magnetic units."""
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)

    mu = ket.get_matrix_element(ket, "MAGNETIC_DIPOLE", q=0)
    mu = mu.to("bohr_magneton")
    lande_factor = 2.002319304363
    assert np.isclose(mu.magnitude, -1 / 2 * lande_factor)

    # check magnetic field conversion is correct
    B_z = QuantityScalar.from_unit(1, "gauss", "MAGNETIC_FIELD")
    B_z_pint = ureg.Quantity(1, "gauss").to("T", "Gaussian")
    assert np.isclose(B_z.to_base_unit(), B_z_pint.to_base_units().magnitude)

    # such that mu * B_z is of dimension energy
    zeeman_energy = -mu * B_z_pint
    assert zeeman_energy.dimensionality == BaseUnits["ENERGY"].dimensionality


def test_electric_dipole() -> None:
    """Test electric dipole units."""
    ket_a = pi.KetAtom("Rb", n=60, l=0, m=0.5)
    ket_b = pi.KetAtom("Rb", n=61, l=0, m=0.5)
    ket_c = pi.KetAtom("Rb", n=60, l=1, j=3 / 2, m=0.5)

    dipole_a_c = ket_a.get_matrix_element(ket_c, "ELECTRIC_DIPOLE", q=0)
    dipole_b_c = ket_b.get_matrix_element(ket_c, "ELECTRIC_DIPOLE", q=0)

    kappa = ureg.Quantity(1 / (4 * np.pi), "1 / epsilon_0")
    C3 = kappa * dipole_a_c * dipole_b_c

    GHz = ureg.Quantity(1, "GHz")
    C3 = C3 * GHz / GHz.to("J", "spectroscopy")
    C3 = C3.to("GHz micrometer^3")

    distance = ureg.Quantity(10, "micrometer")
    basis = pi.BasisAtom("Rb", additional_kets=[ket_a, ket_b, ket_c])
    system = pi.SystemAtom(basis)
    basis_pair = pi.BasisPair([system, system])
    system_pair = pi.SystemPair(basis_pair)
    system_pair.set_order(3)
    system_pair.set_distance(distance)
    system.get_hamiltonian()

    ket_ab_idx = np.argmax(basis_pair.get_overlaps([ket_a, ket_b]))
    ket_cc_idx = np.argmax(basis_pair.get_overlaps([ket_c, ket_c]))
    H = system_pair.get_hamiltonian("GHz") * distance.to("micrometer").magnitude ** 3  # GHz * micrometer^3

    assert np.isclose(-2 * C3.magnitude, H[ket_ab_idx, ket_cc_idx])
    assert np.isclose(-2 * C3.magnitude, 5.73507543166919)
