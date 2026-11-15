"""State-averaged two-site DMRG on itensor_version="python"
(get_excited_states(method="state_averaged"), pyitensor/sadmrg.py),
cross-checked against ED.

The chains are small enough that maxm=64 holds every state exactly (a
10-site spin-1/2 chain needs at most 2^5 = 32 at its middle bond, and a
bond next to the state index no more than the Hilbert space on its own
side), so the only error left is that of the local solves and the
tolerance on the energies is 1e-8: loose against the 1e-10 or better a
converged run reaches, tight against a wrong level (the singlet-triplet
gap of the 10-site chain is about 0.4). The triplet is asserted degenerate
to the same 1e-8, and the returned states orthonormal to 1e-10, since the
final Rayleigh-Ritz step orthonormalizes them exactly up to rounding.
"""

import numpy as np
import pytest

from dmrgpy import cppext, spinchain
from dmrgpy.pyitensor import sadmrg

TOL_E = 1e-8
TOL_ORTHO = 1e-10


def _heisenberg(n=10, field=0.0, dm=0.0, iv="python"):
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["1/2"] * n, itensor_version=iv)
    h = 0
    for i in range(n - 1):
        h = h + sc.SS(i, i + 1)
        if dm:   # Dzyaloshinskii-Moriya: makes the states complex
            h = h + dm * (sc.Sx[i] * sc.Sy[i + 1] - sc.Sy[i] * sc.Sx[i + 1])
    if field:
        for i in range(n):
            h = h + field * sc.Sz[i] + 0.13 * (-1) ** i * sc.Sx[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 64, 12, 1e-12
    return sc


def _ed(sc, n):
    return np.sort(np.real(sc.get_excited(mode="ED", n=n)))


def test_one_state_is_the_ground_state():
    sc = _heisenberg()
    e = np.real(sc.get_excited(n=1, method="state_averaged"))
    assert len(e) == 1
    assert abs(e[0] - _ed(sc, 1)[0]) < TOL_E


@pytest.mark.parametrize("batched", [True, False])
def test_four_states_match_ed_with_the_triplet(batched, monkeypatch):
    # batched=False applies H_eff column by column with dmrg.py's own
    # matvec, so a failure of the block plan alone shows up as one case
    monkeypatch.setattr(sadmrg, "SA_BATCHED_MATVEC", batched)
    sc = _heisenberg()
    es = np.sort(np.real(sc.get_excited(n=4, method="state_averaged")))
    ref = _ed(sc, 4)
    assert np.max(np.abs(es - ref)) < TOL_E
    assert es[1] - es[0] > 1e-2              # singlet below the triplet
    assert es[3] - es[1] < TOL_E             # the triplet, all three members


def test_guard_roots_complete_the_triplet():
    # n=2 cuts the triplet; two guard roots carry all of it
    sc = _heisenberg()
    es = np.sort(np.real(sc.get_excited(n=2, method="state_averaged", nguard=2)))
    assert len(es) == 2
    assert np.max(np.abs(es - _ed(sc, 2))) < TOL_E


@pytest.mark.parametrize("dm", [0.0, 0.3])
def test_returned_states_are_orthonormal(dm):
    sc = _heisenberg(n=8, field=0.3, dm=dm)
    es, ws = sc.get_excited_states(n=4, method="state_averaged")
    O = np.array([[complex(sc.overlap(a, b)) for b in ws] for a in ws])
    assert np.max(np.abs(O - np.eye(len(ws)))) < TOL_ORTHO
    assert np.max(np.abs(np.sort(np.real(es)) - _ed(sc, 4))) < TOL_E


def test_unknown_method_raises():
    sc = _heisenberg(n=6)
    with pytest.raises(ValueError):
        sc.get_excited(n=2, method="no_such_method")


@pytest.mark.skipif(not cppext.available(3), reason="v3 extension not compiled")
def test_other_backends_refuse():
    sc = _heisenberg(n=6, iv=3)
    with pytest.raises(NotImplementedError, match="python"):
        sc.get_excited(n=2, method="state_averaged")
