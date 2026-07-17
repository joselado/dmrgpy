"""Functional coverage for Many_Body_Chain.get_dynamical_correlator,
distilled from examples/hubbard_chain and examples/dynamical_correlator/
dynamical_correlator_minimal down to a small, fast, deterministic
system. Uses the ED "INV" submode (the same one examples/hubbard_chain
compares its DMRG result against), which is exact and fast for small
systems -- no KPM polynomial-order tuning needed."""
import numpy as np
import pytest

from dmrgpy import spinchain

DELTA = 0.05


def test_dynamical_correlator_peaks_at_excitation_gap():
    """4-site Heisenberg chain: the ground state (singlet) has a 3-fold
    degenerate triplet excited state at gap 0.658919 above it (see
    test_excited_states.py). The <Sz_0(t) Sz_0(0)> dynamical correlator
    must show a resonance there and be small away from it -- checked
    both qualitatively (peak >> off-resonance values) and against golden
    regression values at each frequency."""
    n = 4
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)

    name = (sc.Sz[0], sc.Sz[0])
    es = np.array([0.0, 0.658919, 1.5])
    x, y = sc.get_dynamical_correlator(mode="ED", submode="INV", name=name,
                                        es=es, delta=DELTA)

    expected = np.array([0.0066913407688090855, 1.045617999495007, 0.06866722208702186])
    assert y.real == pytest.approx(expected, abs=1e-6)

    # the resonance is the dominant feature by a wide margin
    assert y.real[1] > 10 * y.real[0]
    assert y.real[1] > 10 * y.real[2]
