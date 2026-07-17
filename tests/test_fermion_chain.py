"""Functional coverage for fermionchain.Fermionic_Chain (spinless
fermions), distilled from examples/spinless_fermions and
examples/total_particle_number down to small, fast systems, comparing
ED against DMRG on both ITensor v2 and v3 (see _helpers.py)."""
import pytest

from dmrgpy import fermionchain

from _helpers import energy_ed_v2_v3, vev_ed_v2_v3

DMRG_TOL = 1e-6


def test_spinless_fermion_dimer_hopping_energy():
    """Two-site spinless fermion hopping Hamiltonian H = t(Cdag_0 C_1 +
    h.c.): the single-particle sector diagonalizes to bonding/antibonding
    orbitals at energies -t/+t, and the many-body ground state (vacuum
    has E=0, doubly-occupied has E=0 too since there's no interaction)
    is the singly-occupied bonding orbital at E=-t.

    v2 only (not both v2 and v3): itensor_version=3 crashes hard for any
    exactly-2-physical-site chain regardless of physics type -- see
    _helpers.energy_ed_v2_v3's docstring and test_spin_chain.py's dimer
    test, which hits the same mpscpp3 bug."""
    n = 2
    fc = fermionchain.Fermionic_Chain(n)
    h = fc.Cdag[0] * fc.C[1]
    h = h + h.get_dagger()

    e_ed, e_v2 = energy_ed_v2_v3(fc, h, versions=(2,))
    assert e_ed == pytest.approx(-1.0, abs=1e-8)
    assert e_v2 == pytest.approx(-1.0, abs=DMRG_TOL)


def test_total_particle_number_conservation():
    """4-site spinless fermion chain with hopping, nearest-neighbor
    repulsion and a chemical potential (examples/total_particle_number,
    shrunk from L=4 with a scan over mu to one fixed mu): checks the
    total-N ground-state expectation value against a golden regression
    value (integer, as expected for a symmetry-conserving Hamiltonian)."""
    n = 4
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1]
    h = h + h.get_dagger()
    V = 1.0
    for i in range(n - 1):
        h = h + V * fc.N[i] * fc.N[i + 1]
    mu = -0.5
    for i in range(n):
        h = h + mu * fc.N[i]

    Nop = 0
    for i in range(n):
        Nop = Nop + fc.N[i]

    e_ed, e_v2, e_v3 = energy_ed_v2_v3(fc, h)
    assert e_ed == pytest.approx(-3.0, abs=1e-8)
    assert e_v2 == pytest.approx(-3.0, abs=DMRG_TOL)
    assert e_v3 == pytest.approx(-3.0, abs=DMRG_TOL)

    n_ed, n_v2, n_v3 = vev_ed_v2_v3(fc, h, Nop)
    assert n_ed.real == pytest.approx(2.0, abs=1e-6)
    assert n_v2.real == pytest.approx(2.0, abs=DMRG_TOL)
    assert n_v3.real == pytest.approx(2.0, abs=DMRG_TOL)
