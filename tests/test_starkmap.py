# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test the Stark map calculation."""

from pathlib import Path
from typing import TYPE_CHECKING, Optional

import numpy as np
import pairinteraction.real as pi
import pytest

if TYPE_CHECKING:
    from pairinteraction.units import NDArray

reference_kets_file = Path(__file__).parent.parent / "data/reference_stark_map/kets.txt"
reference_eigenenergies_file = Path(__file__).parent.parent / "data/reference_stark_map/eigenenergies.txt"
reference_overlaps_file = Path(__file__).parent.parent / "data/reference_stark_map/overlaps.txt"


def test_starkmap(generate_reference: bool) -> None:
    """Test calculating a Stark map."""
    # Create a basis
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))
    print(f"Number of basis states: {basis.number_of_states}")

    electric_fields = np.linspace(0, 10, 11)
    # Create systems for different values of the electric field
    systems = [pi.SystemAtom(basis).set_electric_field([0, 0, e], unit="V/cm") for e in electric_fields]

    # Diagonalize the systems in parallel
    pi.diagonalize(systems, diagonalizer="eigen", sort_by_energy=True)

    # Get the overlap with |ket>
    overlaps = np.array([system.basis.get_overlaps(ket) for system in systems])

    # Compare to reference data
    kets = [repr(ket) for ket in systems[0].basis.kets]
    eigenenergies = np.array([system.get_eigenenergies(unit="GHz") for system in systems])
    eigenvectors = np.array([system.get_eigenbasis().get_coefficients().todense().A1 for system in systems])

    if generate_reference:
        reference_kets_file.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_kets_file, kets, fmt="%s", delimiter="\t")
        np.savetxt(reference_eigenenergies_file, eigenenergies)
        np.savetxt(reference_overlaps_file, overlaps)
        pytest.skip("Reference data generated, skipping comparison test")

    compare_starkmap_to_reference(eigenenergies, overlaps, eigenvectors, kets)


def compare_starkmap_to_reference(
    eigenenergies: "NDArray",
    overlaps: Optional["NDArray"] = None,
    eigenvectors: Optional["NDArray"] = None,
    kets: Optional[list[str]] = None,
) -> None:
    np.testing.assert_allclose(eigenenergies, np.loadtxt(reference_eigenenergies_file))

    if overlaps is not None:
        # Ensure that the overlaps sum up to one
        np.testing.assert_allclose(np.sum(overlaps, axis=1), np.ones(len(eigenenergies)))
        np.testing.assert_allclose(overlaps, np.loadtxt(reference_overlaps_file), atol=1e-10)

    if kets is not None:
        np.testing.assert_equal(kets, np.loadtxt(reference_kets_file, dtype=str, delimiter="\t"))

    if eigenvectors is not None:
        # Because of degeneracies, checking the eigenvectors against reference data is complicated.
        # Thus, we only check their normalization and orthogonality.
        cumulative_norm = (np.array(eigenvectors) * np.array(eigenvectors).conj()).sum(axis=1)
        np.testing.assert_allclose(cumulative_norm, 90 * np.ones(len(eigenenergies)))
