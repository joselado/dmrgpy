"""Benchmark the current MultiOperator implementation's DMRG/ED energies
and correlators against golden values computed with the
pre-optimization implementation (see reference_data.py).

These are the main "did the optimization change any physics" checks:
they build genuinely non-trivial Hamiltonians (via the "H = H + term"
loop idiom, set_exchange, get_dagger, ...) and run real DMRG/ED
calculations through the C++ backend and the ED fallback.
"""
import pytest

from dmrgpy import spinchain, fermionchain

from reference_data import (
    HEISENBERG_DMRG_ENERGY,
    HEISENBERG_ED_ENERGY,
    SET_EXCHANGE_NTERMS,
    SET_EXCHANGE_DMRG_ENERGY,
    SET_EXCHANGE_ED_ENERGY,
    FERMION_CORRELATOR_U0,
    FERMION_CORRELATOR_U2,
    FERMION_CORRELATOR_U10,
)

DMRG_TOL = 1e-6  # DMRG is iterative/variational, not bit-exact
ED_TOL = 1e-8


def test_heisenberg_chain_energy_matches_reference():
    """10-site Heisenberg chain built with the "h = h + term" idiom
    used throughout this codebase; ground-state energy via DMRG and ED
    must match the original implementation."""
    n = 10
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1]
        h = h + sc.Sy[i] * sc.Sy[i + 1]
        h = h + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    assert sc.gs_energy(mode="DMRG") == pytest.approx(HEISENBERG_DMRG_ENERGY, abs=DMRG_TOL)
    assert sc.gs_energy(mode="ED") == pytest.approx(HEISENBERG_ED_ENERGY, abs=ED_TOL)


def test_set_exchange_energy_and_term_count_matches_reference():
    """spinchain.set_exchange (the O(9*L^2) exchange-coupling builder,
    now assembled via multioperator.msum) must produce the same term
    count (after an explicit clean()) and ground-state energy as the
    original implementation."""
    n = 6
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)

    def fj(i, j):
        return 1.0 if abs(i - j) == 1 else 0.0

    sc.set_exchange(fj)
    h = sc.exchange.copy()
    h.clean()
    assert len(h.op) == SET_EXCHANGE_NTERMS
    assert sc.gs_energy(mode="DMRG") == pytest.approx(SET_EXCHANGE_DMRG_ENERGY, abs=DMRG_TOL)
    assert sc.gs_energy(mode="ED") == pytest.approx(SET_EXCHANGE_ED_ENERGY, abs=ED_TOL)


@pytest.mark.parametrize("U,expected", [
    (0., FERMION_CORRELATOR_U0),
    (2., FERMION_CORRELATOR_U2),
    (10., FERMION_CORRELATOR_U10),
])
def test_spinful_fermion_correlator_matches_reference(U, expected):
    """Spinful fermion chain (nearest-neighbor hopping + Hubbard U),
    built with the "H = H + term" idiom and H.get_dagger(); the
    <Cdag_0up C_iup> ED correlator must match the original
    implementation for several interaction strengths."""
    L = 6
    fc = fermionchain.Spinful_Fermionic_Chain(L)
    H = 0
    for i in range(L - 1):
        H = H + fc.Cdagup[i] * fc.Cup[i + 1]
        H = H + fc.Cdagdn[i] * fc.Cdn[i + 1]
    H = H + H.get_dagger()
    for i in range(L):
        H = H + U * (fc.Nup[i] - 0.5) * (fc.Ndn[i] - 0.5)
    fc.set_hamiltonian(H)
    wf0 = fc.get_gs(mode="ED")
    correlator = [wf0.dot((fc.Cdagup[0] * fc.Cup[i]) * wf0).real for i in range(L)]
    assert correlator == pytest.approx(expected, abs=ED_TOL)
