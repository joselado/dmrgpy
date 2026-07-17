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

from dmrgpy import spinchain
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
