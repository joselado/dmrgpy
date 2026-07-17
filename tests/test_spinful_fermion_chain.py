"""Functional coverage for fermionchain.Spinful_Fermionic_Chain (the
Hubbard-model chain), distilled from examples/hubbard_chain and
examples/fermionic_static_correlator down to small, fast systems,
comparing ED against DMRG on both ITensor v2 and v3 (see _helpers.py)."""
import pytest

from dmrgpy import fermionchain

from _helpers import energy_ed_v2_v3, vev_ed_v2_v3

DMRG_TOL = 1e-6


def test_hubbard_chain_ground_state_energy():
    """3-site Hubbard chain (hopping + particle-hole symmetric on-site
    U), as in examples/hubbard_chain shrunk from n=4: ground-state
    energy checked against a golden regression value."""
    n = 3
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + fc.Cdagdn[i] * fc.Cdn[i + 1]
    U = 2.0
    for i in range(n):
        h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
    h = h + h.get_dagger()

    e_ed, e_v2, e_v3 = energy_ed_v2_v3(fc, h)
    assert e_ed == pytest.approx(-4.236067977499792, abs=1e-8)
    assert e_v2 == pytest.approx(e_ed, abs=DMRG_TOL)
    assert e_v3 == pytest.approx(e_ed, abs=DMRG_TOL)


def test_static_hopping_correlator():
    """Static correlator <Cdagup_0 Cup_j> on the same 3-site Hubbard
    chain (examples/fermionic_static_correlator), checked against golden
    regression values for each site j."""
    n = 3
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + fc.Cdagdn[i] * fc.Cdn[i + 1]
    U = 2.0
    for i in range(n):
        h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
    h = h + h.get_dagger()

    expected = [0.1881966011250103, -0.26180339887498927, 0.11180339887498948]
    for j in range(n):
        op = fc.Cdagup[0] * fc.Cup[j]
        v_ed, v_v2, v_v3 = vev_ed_v2_v3(fc, h, op)
        assert v_ed.real == pytest.approx(expected[j], abs=1e-6)
        assert v_v2.real == pytest.approx(expected[j], abs=DMRG_TOL)
        assert v_v3.real == pytest.approx(expected[j], abs=DMRG_TOL)
