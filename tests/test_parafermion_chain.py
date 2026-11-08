"""Functional coverage for parafermionchain.Parafermionic_Chain,
distilled from examples/parafermion_energy down to a small, fast, and
deterministic system, comparing ED against DMRG on both ITensor v2 and
v3 (see _helpers.py).

Note: the Hamiltonian below is deliberately built without ever adding
a plain Python int/float "0" as the initial accumulator (unlike the
"h = 0; h = h + term" idiom used everywhere else in this codebase). A
bare 0 promotes to a MultiOperator term tagged with a literal "Id"
operator (multioperator.py's scalar-promotion path), and the ED
backend for this chain (pyparafermion/parafermion.py) never registers
an "Id" operator for any site -- so any Hamiltonian that ends up
carrying that placeholder term hits a KeyError in ED. Starting the
accumulation from the first real term avoids ever creating it.
"""
import numpy as np
import pytest

from dmrgpy import parafermionchain

from _helpers import energy_ed_v2_v3
from _helpers import julia_live_param
from _julia_site_helpers import check_julia_live_against_ed

DMRG_TOL = 1e-6


def test_parafermion_chain_energy():
    """3-site Z3 parafermion chain with a deterministic (seeded) mix of
    number-operator and clock-operator couplings, checked against a
    golden regression value for the ground-state energy."""
    np.random.seed(7)
    n = 3
    pc = parafermionchain.Parafermionic_Chain(n)
    h = -1.0 * pc.N[0]
    for i in range(n):
        if i > 0:
            h = h - pc.N[i]
        for j in range(n):
            if i < j:
                c = 0.4
                h = h + c * (pc.Tau[i] * pc.Taud[j])
                h = h + c * (pc.Tau[j] * pc.Taud[i])

    e_ed, e_v2, e_v3 = energy_ed_v2_v3(pc, h)
    assert e_ed == pytest.approx(-5.0, abs=1e-8)
    assert e_v2 == pytest.approx(-5.0, abs=DMRG_TOL)
    assert e_v3 == pytest.approx(-5.0, abs=DMRG_TOL)


@pytest.mark.parametrize("Z", [2, 3, 4])
@pytest.mark.parametrize("v", [julia_live_param()])
def test_parafermion_chain_julia_live_matches_ed(Z, v):
    """julia_live builds the Z2/Z3/Z4 clock sites of sites_custom.jl:
    4-site chain with a transverse term Tau + Tau^dagger and the
    parafermion hopping Chi_i^dagger Chi_(i+1) + h.c., whose Tau strings
    are explicit in the terms, against ED"""
    n = 4
    pc = parafermionchain.Parafermionic_Chain(n, Z=Z)
    h = -1.2*(pc.Tau[0] + pc.Taud[0])
    for i in range(1, n):
        h = h - 1.2*(pc.Tau[i] + pc.Taud[i])
    for i in range(n - 1):
        hop = pc.Chid[i]*pc.Chi[i + 1]
        h = h - (hop + hop.get_dagger())
    h = h + 0.3*pc.N[0]
    obs = [pc.Tau[1], pc.N[0], pc.Sig[0]*pc.Sigd[1], pc.Chid[0]*pc.Chi[2]]
    check_julia_live_against_ed(pc, h, obs, maxm=40, nsweeps=20)
