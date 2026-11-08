"""Functional coverage for fermionchain.Spinful_Fermionic_Chain_Native (the
native-spinful-site Hubbard-model chain: one dimension-4 Electron/Hubbard
tensor-network site per orbital, mpscpp3/get_sites.h's site-type code 1,
instead of Spinful_Fermionic_Chain's two interleaved dimension-2 sites per
orbital). Same physical Hamiltonians/expected numbers as
test_spinful_fermion_chain.py -- the two classes describe the exact same
physics, so the golden values are identical -- restricted to
itensor_version=3, the only DMRG backend this class wires up (see the
class docstring)."""
import pytest

from dmrgpy import fermionchain

from _helpers import energy_ed_v2_v3, vev_ed_v2_v3
from _helpers import julia_live_param
from _julia_site_helpers import check_julia_live_against_ed

DMRG_TOL = 1e-6


def test_hubbard_chain_ground_state_energy():
    """3-site Hubbard chain (hopping + particle-hole symmetric on-site
    U): ground-state energy checked against the same golden regression
    value as the interleaved-site test_spinful_fermion_chain.py."""
    n = 3
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + fc.Cdagdn[i] * fc.Cdn[i + 1]
    U = 2.0
    for i in range(n):
        h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
    h = h + h.get_dagger()

    e_ed, e_v3 = energy_ed_v2_v3(fc, h, versions=(3,))
    assert e_ed == pytest.approx(-4.236067977499792, abs=1e-8)
    assert e_v3 == pytest.approx(e_ed, abs=DMRG_TOL)


def test_static_hopping_correlator():
    """Static correlator <Cdagup_0 Cup_j> on the same 3-site Hubbard
    chain, checked against the same golden regression values as
    test_spinful_fermion_chain.py's version of this test (see that file's
    docstring for why the explicit Zeeman-like term is needed to lift the
    up/down ground-state degeneracy)."""
    n = 3
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdagup[i] * fc.Cup[i + 1]
        h = h + fc.Cdagdn[i] * fc.Cdn[i + 1]
    U = 2.0
    for i in range(n):
        h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
    h = h + h.get_dagger()
    eps = 0.3  # lifts the up/down degeneracy
    h = h + eps * (fc.Nup[0] - fc.Ndn[0])

    expected = [0.09204667896906626, -0.20764359692352258, 0.09821625446016191]
    for j in range(n):
        op = fc.Cdagup[0] * fc.Cup[j]
        v_ed, v_v3 = vev_ed_v2_v3(fc, h, op, versions=(3,))
        assert v_ed.real == pytest.approx(expected[j], abs=1e-6)
        assert v_v3.real == pytest.approx(expected[j], abs=DMRG_TOL)


def test_matches_interleaved_chain_ground_state():
    """Spinful_Fermionic_Chain_Native and Spinful_Fermionic_Chain describe
    the exact same physical Hamiltonian (same Jordan-Wigner sign
    convention, see jordanwigner_spinful.py), just packaged into
    different tensor-network sites -- their ED ground-state energies must
    agree exactly (not merely to DMRG tolerance) for the same couplings.
    """
    n = 4
    U = 1.7
    fc_native = fermionchain.Spinful_Fermionic_Chain_Native(n)
    fc_double = fermionchain.Spinful_Fermionic_Chain(n)

    def build(fc):
        h = 0
        for i in range(n - 1):
            h = h + fc.Cdagup[i] * fc.Cup[i + 1]
            h = h + fc.Cdagdn[i] * fc.Cdn[i + 1]
        for i in range(n):
            h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
        return h + h.get_dagger()

    fc_native.set_hamiltonian(build(fc_native))
    fc_double.set_hamiltonian(build(fc_double))
    e_native = fc_native.gs_energy(mode="ED")
    e_double = fc_double.gs_energy(mode="ED")
    assert e_native == pytest.approx(e_double, abs=1e-10)


def test_kpm_dynamical_correlator_matches_interleaved_chain():
    """The KPM dynamical correlator <Cup_0(t) Cdagup_0(0)>-type spectral
    function must agree between the native and interleaved chains (same
    physics, same Hamiltonian, cross-checked via DMRG itensor_version=3
    on both -- KPM is stochastic-free but still an approximate moment
    expansion, so a loose tolerance is used)."""
    import numpy as np

    n = 4
    U = 2.0

    def build(chain_cls):
        fc = chain_cls(n)
        h = 0
        for i in range(n - 1):
            h = h + fc.Cdagup[i] * fc.Cup[i + 1]
            h = h + fc.Cdagdn[i] * fc.Cdn[i + 1]
        for i in range(n):
            h = h + U * (fc.Nup[i] - .5) * (fc.Ndn[i] - .5)
        h = h + h.get_dagger()
        fc.set_hamiltonian(h)
        fc.setup_cpp(3)
        return fc

    fc_native = build(fermionchain.Spinful_Fermionic_Chain_Native)
    fc_double = build(fermionchain.Spinful_Fermionic_Chain)

    es = np.linspace(-4, 4, 60)
    _, y_native = fc_native.get_dynamical_correlator(
            name=(fc_native.Cup[0], fc_native.Cdagup[0]), es=es, delta=0.3)
    _, y_double = fc_double.get_dynamical_correlator(
            name=(fc_double.Cup[0], fc_double.Cdagup[0]), es=es, delta=0.3)

    assert y_native == pytest.approx(y_double, abs=5e-3)


@pytest.mark.parametrize("v", [julia_live_param()])
def test_native_hubbard_chain_julia_live_matches_ed(v):
    """julia_live builds ITensors.jl's Electron site (Emp, Up, Dn, UpDn)
    and threads the Jordan-Wigner strings itself: 4-site Hubbard chain,
    t=1, U=4, a Zeeman field, an onsite Sx (a same-site Cdagup Cdn) and
    site potentials, against ED, with hoppings that cross a site"""
    n = 4
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdagup[i]*fc.Cup[i + 1] + fc.Cdagdn[i]*fc.Cdn[i + 1]
    h = h + h.get_dagger()
    for i in range(n):
        h = h + 4.0*(fc.Nup[i] - .5)*(fc.Ndn[i] - .5)
        h = h + 0.2*fc.Sz[i] + 0.05*(i + 1)*fc.Ntot[i]
    h = h + 0.15*fc.Sx[0]
    obs = [fc.Nup[1], fc.Sz[2], fc.Sx[0], fc.Cdagup[0]*fc.Cup[2],
           fc.Cdagdn[0]*fc.Cdn[3], fc.Cdagup[1]*fc.Cdn[2]]
    check_julia_live_against_ed(fc, h, obs, maxm=60, nsweeps=20)
