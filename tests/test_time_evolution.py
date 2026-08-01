"""Functional coverage for dmrgpy.timeevolution.evolve_WF, distilled
from examples/time_evolution/time_evolution_wavefunction down to a
small, fast, deterministic system.

Unlike the other test files here, this doesn't do an ED/v2/v3
cross-backend comparison: the DMRG-side time evolution
(mpsalgebra.exponential_dmrg) is a Trotterized/fitted approximation,
while the ED side is an exact matrix exponential, so the two are not
expected to agree to DMRG_TOL-style precision without also tuning the
Trotter step count -- exercising exact ED evolution here is enough to
cover the feature.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain
from dmrgpy import spinchain
from dmrgpy import timedependent
from dmrgpy.timeevolution import evolve_WF


def test_time_evolution_conserves_norm_and_matches_reference_trajectory():
    """Two-site chain: prepare a Neel-like initial state as the ground
    state of a staggered field, then time-evolve under the Heisenberg
    Hamiltonian. The evolution must be norm-conserving (it's a unitary
    exponential), and <Sz_0>(t) must match a golden regression
    trajectory."""
    n = 2
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)

    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf0 = sc.get_gs(mode="ED")

    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    h = h + h.get_dagger()

    ts = np.linspace(0., 4., 5)
    wfs = evolve_WF(h, wf0, ts=ts)

    norms = np.array([w.dot(w) for w in wfs])
    assert norms.real == pytest.approx(np.ones(len(ts)), abs=1e-8)
    assert norms.imag == pytest.approx(np.zeros(len(ts)), abs=1e-8)

    sz0 = np.array([w.dot(sc.Sz[0] * w).real for w in wfs])
    expected = np.array([
        -0.5,
        0.20807341827357126,
        0.32682181043180597,
        -0.48008514332518304,
        0.07275001690430655,
    ])
    assert sz0 == pytest.approx(expected, abs=1e-6)


def test_tebd_matches_ed_spin_chain():
    """pyitensor/tebd.py's TEBDEvolver (self.tevol_method="TEBD",
    itensor_version="python" only) against exact diagonalization: prepare
    the GS of a Neel-favoring staggered field, quench to the (strictly
    nearest-neighbor) Heisenberg Hamiltonian, and require <Sz_0>(t) to
    match ED closely -- modeled on
    examples/time_evolution/tdvp_VS_ED_time_evolution, the same pattern
    used to validate TDVP."""
    n = 4
    spins = [2 for _ in range(n)]  # S=1/2
    sc = spinchain.Spin_Chain(spins)
    sc.setup_python()
    sc.tevol_method = "TEBD"

    h0 = 0
    for i in range(n):
        h0 = h0 + (-1) ** i * sc.Sz[i]
    h1 = 0
    for i in range(n - 1):
        h1 = h1 + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]

    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    wfED = sc.get_gs(mode="ED")
    sc.set_hamiltonian(h1)

    nt, dt = 50, 0.05
    (_ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wf)
    (_tsED, szED) = timedependent.evolve_and_measure(
            sc, operator=sc.Sz[0], nt=nt, dt=dt, wf=wfED, mode="ED")

    assert np.array(sz).real == pytest.approx(np.array(szED).real, abs=1e-4)


def test_tebd_matches_ed_fermion_chain():
    """Same cross-check as test_tebd_matches_ed_spin_chain, but for a
    spinless-fermion nearest-neighbor hopping Hamiltonian, to exercise
    TEBD's Jordan-Wigner handling (bond_hamiltonians() extracts each
    term's two adjacent-site matrices directly from HTerm.resolve(),
    relying on JW strings never needing to extend past an already-
    adjacent bond -- this is the case that would break if that reasoning
    were wrong)."""
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    fc.tevol_method = "TEBD"

    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)

    nt, dt = 40, 0.05
    (_x, y) = timedependent.evolution_ABA(fc, nt=nt, dt=dt, mode="ED", A=fc.Cdag[0], B=fc.N[0])
    (_x1, y1) = timedependent.evolution_ABA(fc, nt=nt, dt=dt, mode="DMRG", A=fc.Cdag[0], B=fc.N[0])

    assert np.array(y1) == pytest.approx(np.array(y), abs=1e-4)


def test_tebd_rejects_non_nearest_neighbor_hamiltonian():
    """bond_hamiltonians() must fail loudly, not silently drop the
    long-range piece, when the Hamiltonian isn't strictly nearest-
    neighbor -- TEBD has no representation for a term spanning 3+ sites."""
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()
    fc.tevol_method = "TEBD"

    h = fc.Cdag[0] * fc.C[2]
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)

    with pytest.raises(NotImplementedError):
        timedependent.evolution_ABA(fc, nt=5, dt=0.05, mode="DMRG", A=fc.Cdag[0], B=fc.N[0])
