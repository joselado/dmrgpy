"""Functional coverage for spinchain.Spin_Chain, distilled from the
examples/heisenberg_12, examples/transverse_ising_model and
examples/total_spin scripts down to small, fast systems. Each test
compares ED against DMRG on both the ITensor v2 and v3 C++ backends
(see _helpers.energy_ed_v2_v3/vev_ed_v2_v3): when a backend isn't
compiled, mode.py transparently falls back to ED for that call, so
these still run (as ED-only smoke tests) without a compiled extension,
and become real cross-backend regression checks wherever both are
built.
"""
import numpy as np
import pytest

from dmrgpy import spinchain

from _helpers import energy_ed_v2_v3, vev_ed_v2_v3

DMRG_TOL = 1e-6


def test_heisenberg_dimer_exact_singlet_energy():
    """Two-site Heisenberg dimer: exact ground state is the singlet,
    with energy -3/4 (in units where the exchange coupling is 1).

    v3 only (not both v2 and v3): itensor_version=3 crashes hard
    ("LocalOp is default constructed", an ITensor v3 internal check)
    for any exactly-2-physical-site chain, independent of physics type
    -- a genuine bug in mpscpp3, not something this test can route
    around via mode.py's ED fallback. See _helpers.energy_ed_v2_v3's
    docstring and test_fermion_chain.py's dimer test, which hits the
    same thing."""
    spins = ["S=1/2", "S=1/2"]
    sc = spinchain.Spin_Chain(spins)
    h = sc.Sx[0] * sc.Sx[1] + sc.Sy[0] * sc.Sy[1] + sc.Sz[0] * sc.Sz[1]

    e_ed, e_v2 = energy_ed_v2_v3(sc, h, versions=(2,))
    assert e_ed == pytest.approx(-0.75, abs=1e-8)
    assert e_v2 == pytest.approx(-0.75, abs=DMRG_TOL)


def test_transverse_field_ising_energy_and_zero_magnetization():
    """4-site transverse-field Ising chain (examples/transverse_ising_model,
    shrunk from n=40): the Hamiltonian is symmetric under a global spin
    flip along z, so <Sz_total> must vanish by symmetry regardless of the
    field strength; the ground-state energy is checked against a golden
    regression value."""
    n = 4
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h - sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + 1.0 * sc.Sx[i]

    Mz = 0
    for i in range(n):
        Mz = Mz + sc.Sz[i]

    e_ed, e_v2, e_v3 = energy_ed_v2_v3(sc, h)
    assert e_ed == pytest.approx(-2.0941996592125895, abs=1e-8)
    assert e_v2 == pytest.approx(e_ed, abs=DMRG_TOL)
    assert e_v3 == pytest.approx(e_ed, abs=DMRG_TOL)

    mz_ed, mz_v2, mz_v3 = vev_ed_v2_v3(sc, h, Mz)
    assert mz_ed.real == pytest.approx(0.0, abs=1e-6)
    assert mz_v2.real == pytest.approx(0.0, abs=DMRG_TOL)
    assert mz_v3.real == pytest.approx(0.0, abs=DMRG_TOL)


def test_odd_chain_ground_state_is_spin_doublet():
    """3-site Heisenberg chain (examples/total_spin, shrunk and made
    uniform spin-1/2): by the Lieb-Mattis theorem an antiferromagnetic
    Heisenberg chain with an odd number of spin-1/2 sites has a
    ground-state total spin S=1/2, so <S_total^2> = S(S+1) = 3/4
    exactly."""
    n = 3
    spins = ["S=1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]

    sxt = sum(sc.Sx)
    syt = sum(sc.Sy)
    szt = sum(sc.Sz)
    S2 = sxt * sxt + syt * syt + szt * szt

    e_ed, e_v2, e_v3 = energy_ed_v2_v3(sc, h)
    assert e_ed == pytest.approx(-1.0, abs=1e-8)
    assert e_v2 == pytest.approx(-1.0, abs=DMRG_TOL)
    assert e_v3 == pytest.approx(-1.0, abs=DMRG_TOL)

    s2_ed, s2_v2, s2_v3 = vev_ed_v2_v3(sc, h, S2)
    assert s2_ed.real == pytest.approx(0.75, abs=1e-6)
    assert s2_v2.real == pytest.approx(0.75, abs=DMRG_TOL)
    assert s2_v3.real == pytest.approx(0.75, abs=DMRG_TOL)
