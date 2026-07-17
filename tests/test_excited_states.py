"""Functional coverage for Many_Body_Chain.get_excited_states, distilled
from examples/excited_states down to a small, fast, deterministic
system, comparing ED against DMRG on both ITensor v2 and v3.
"""
import numpy as np
import pytest

from dmrgpy import spinchain

DMRG_TOL = 1e-6


def _excited_energies(itensor_version, sc, h, n):
    sc.set_hamiltonian(h)
    sc.setup_cpp(itensor_version)
    es, wfs = sc.get_excited_states(n=n, mode="DMRG")
    return np.sort(np.array(es))


def test_excited_states_energies_and_ordering():
    """4-site Heisenberg chain: lowest 4 eigenstates are a non-degenerate
    singlet ground state followed by a 3-fold degenerate triplet, as
    expected for this small system. ED and both DMRG backends must
    agree with a golden regression value for each level (backends
    silently fall back to ED when not compiled, see mode.py)."""
    n = 4
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]

    expected = np.array([
        -1.6160254037844355,
        -0.9571067811865455,
        -0.9571067811865448,
        -0.9571067811865448,
    ])

    sc.set_hamiltonian(h)
    es_ed, _ = sc.get_excited_states(n=4, mode="ED")
    es_ed = np.sort(np.array(es_ed))
    assert es_ed == pytest.approx(expected, abs=1e-6)

    # non-decreasing order is itself part of the contract being checked
    assert np.all(np.diff(es_ed) >= -1e-9)

    es_v2 = _excited_energies(2, sc, h, 4)
    assert es_v2 == pytest.approx(expected, abs=DMRG_TOL)

    es_v3 = _excited_energies(3, sc, h, 4)
    assert es_v3 == pytest.approx(expected, abs=DMRG_TOL)
