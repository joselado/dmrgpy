"""The "python" Davidson local solve (pyitensor/dmrg.py, since 2026-10-25)
stopped short on a nearly diagonal H: the preconditioned correction
(D - theta)^-1 r was almost the current vector, vanished after
orthogonalization, and the solve stopped with a residual of 5e-4 on
energies of 1e-2, so the excited states of a field-polarized S=1 chain were
1e-5 to 6e-5 above ED (tests/test_audit_2026_09_24b_groundstate.py's n_gs
and EX tests).

And groundstate.best_gs kept the lowest of n solves that were the same
solve: the session warm-started every solve after the first from the
previous answer. A "python" solve from a random start on a 3-site flux ring
ends in a one-particle state instead of the vacuum in about 1 of 40 starts
(left open, docs/development.md 4.5), so best_gs(n=3) is how
tests/test_audit_2026_09_24b_realtime.py::test_flux_ring_circulates_forward
gets its ground state, and that needs n fresh starts.
"""

import warnings

import numpy as np
import pytest

from dmrgpy import fermionchain, groundstate, spinchain


def test_excited_states_of_a_nearly_diagonal_chain_match_ed():
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version="python")
    sc.maxm, sc.nsweeps = 30, 15
    D = 1e-3
    h = D*sc.Sz[0]*sc.Sz[0] + (D + 1e-5)*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    np.random.seed(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        es = np.array(sc.get_excited(n=3, mode="DMRG"))
    es_ed = np.array(sc.get_excited(n=3, mode="ED"))
    assert np.max(np.abs(es - es_ed)) < 1e-12


def test_best_gs_starts_every_solve_fresh(monkeypatch):
    """On "python", where the session's methods can be wrapped; v2 and v3
    take the same call (clear_wavefunction, which
    tests/test_restart_drops_session_state.py exercises there)."""
    fc = fermionchain.Fermionic_Chain(4, itensor_version="python")
    h = 0
    for j in range(3):
        h = h - (fc.Cdag[j]*fc.C[j + 1] + fc.Cdag[j + 1]*fc.C[j])
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 10, 4
    fc.gs_energy()
    calls = []
    session = fc._session
    orig = session.clear_wavefunction
    monkeypatch.setattr(session, "clear_wavefunction",
                        lambda: (calls.append(1), orig())[1])
    groundstate.best_gs(fc, n=3)
    assert len(calls) == 3
    assert fc.gs_energy() == pytest.approx(fc.gs_energy(mode="ED"), abs=1e-10)
