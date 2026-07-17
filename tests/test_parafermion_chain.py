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
