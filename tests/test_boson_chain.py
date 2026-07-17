"""Functional coverage for bosonchain.Bosonic_Chain, distilled from
examples/bosonic_chain down to a small, fast system, comparing ED
against DMRG on both ITensor v2 and v3 (see _helpers.py)."""
import pytest

from dmrgpy import bosonchain

from _helpers import energy_ed_v2_v3, vev_ed_v2_v3

DMRG_TOL = 1e-6


def test_bosonic_chain_energy_and_occupation():
    """3-site bosonic chain with hopping and an on-site interaction that
    favors single occupation (U*(N-1)^2), as in examples/bosonic_chain
    shrunk from n=4: ground-state energy and per-site occupation checked
    against golden regression values."""
    n = 3
    bc = bosonchain.Bosonic_Chain(n)
    h = 0
    t = 1.0
    U = 1.0
    for i in range(n - 1):
        h = h + t * bc.Adag[i] * bc.A[i + 1]
    for i in range(n):
        h = h + U * (bc.N[i] - 1.0) * (bc.N[i] - 1.0)
    h = h + h.get_dagger()

    e_ed, e_v2, e_v3 = energy_ed_v2_v3(bc, h)
    assert e_ed == pytest.approx(-1.7649845186790567, abs=1e-8)
    assert e_v2 == pytest.approx(e_ed, abs=DMRG_TOL)
    assert e_v3 == pytest.approx(e_ed, abs=DMRG_TOL)

    expected_occ = [1.2709856540948459, 1.4580286918103107, 1.2709856540948437]
    for i in range(n):
        occ_ed, occ_v2, occ_v3 = vev_ed_v2_v3(bc, h, bc.N[i])
        assert occ_ed.real == pytest.approx(expected_occ[i], abs=1e-6)
        assert occ_v2.real == pytest.approx(expected_occ[i], abs=DMRG_TOL)
        assert occ_v3.real == pytest.approx(expected_occ[i], abs=DMRG_TOL)
