"""The Hamiltonian of a v3 session with the names of its two-fermion terms
kept, so that ITensor's AutoMPO threads their Jordan-Wigner strings
(groundstate.native_strings, MultiOperator.to_terms(native_strings=)).

The route is taken only above groundstate.NATIVE_STRING_FACTORS; the
tests that check its physics lower that threshold to -1, so that a small
chain takes it and ED can be the reference. A complex hopping, a
pairing, a parity-odd field and an interaction are all in one
Hamiltonian, since a sign error of the string shows in one of them."""

import numpy as np
import pytest

from dmrgpy import fermionchain, groundstate


def _cpp3():
    from dmrgpy import cppext
    return pytest.mark.skipif(not cppext.available(3),
                              reason="the v3 extension is not compiled")


pytestmark = _cpp3()


@pytest.fixture
def forced(monkeypatch):
    monkeypatch.setattr(groundstate, "NATIVE_STRING_FACTORS", -1)


def _long_range(fc, n, odd=0.0):
    rng = np.random.default_rng(3)
    H = 0
    for i in range(n):
        H = H + rng.uniform(-0.5, 0.5)*fc.N[i]
        for j in range(i+1, n):
            t = rng.normal()*np.exp(1j*rng.uniform(0, 2*np.pi))/(1+j-i)
            H = H + t*fc.Cdag[i]*fc.C[j] + np.conj(t)*fc.Cdag[j]*fc.C[i]
            d = 0.3*rng.normal()/(1+j-i)  # pairing, real and antisymmetric
            H = H + d*(fc.Cdag[i]*fc.Cdag[j] + fc.C[j]*fc.C[i])
        if i < n-1:
            H = H + 0.7*fc.N[i]*fc.N[i+1]
    if odd:
        H = H + odd*(fc.C[3] + fc.Cdag[3])  # a parity-odd field
    return H


@pytest.mark.parametrize("odd", [0.0, 0.2])
def test_native_route_against_ed(forced, odd):
    n = 8
    fc = fermionchain.Fermionic_Chain(n, itensor_version=3)
    fc.maxm, fc.nsweeps = 80, 30
    fc.set_hamiltonian(_long_range(fc, n, odd))
    assert groundstate.native_strings(fc)
    key, terms = groundstate._send_key(fc)
    names = {name for (c, f) in terms for (name, s) in f}
    assert "Cdag" in names and "C" in names  # the two-fermion terms
    if odd:
        assert "A" in names  # the odd field keeps its explicit string
    e_ed = fc.gs_energy(mode="ED")
    assert fc.gs_energy() == pytest.approx(e_ed, abs=1e-8)
    for op in (fc.Cdag[1]*fc.C[6], fc.N[2]*fc.N[5], fc.C[0]*fc.C[7]):
        assert fc.vev(op) == pytest.approx(fc.vev(op, mode="ED"), abs=1e-6)


def test_odd_operator_on_the_state_of_the_native_route(forced):
    """C applied to the state solved on the native route: <C_i^dag C_j>
    as the overlap of C_i|GS> and C_j|GS> (the applied operators keep
    their explicit strings) against ED."""
    n = 8
    fc = fermionchain.Fermionic_Chain(n, itensor_version=3)
    fc.maxm, fc.nsweeps = 80, 30
    fc.set_hamiltonian(_long_range(fc, n))
    wf = fc.get_gs()
    for (i, j) in [(1, 6), (0, 7), (3, 4)]:
        ci = fc.applyoperator(fc.C[i], wf)
        cj = fc.applyoperator(fc.C[j], wf)
        g = fc.overlap(ci, cj)
        assert g == pytest.approx(fc.vev(fc.Cdag[i]*fc.C[j], mode="ED"),
                                  abs=1e-6)


def test_native_route_in_a_sector(forced):
    n = 5
    fc = fermionchain.Spinful_Fermionic_Chain(n, itensor_version=3)
    fc.maxm, fc.nsweeps = 80, 80  # at 30 sweeps both routes stop 2.4e-5 up
    rng = np.random.default_rng(1)
    H = 0
    for i in range(n):
        for j in range(i+1, n):
            t = rng.normal()/(1+j-i)
            H = H + t*(fc.Cdagup[i]*fc.Cup[j] + fc.Cdagdn[i]*fc.Cdn[j])
            H = H + t*(fc.Cdagup[j]*fc.Cup[i] + fc.Cdagdn[j]*fc.Cdn[i])
    H = H + 2.0*sum(fc.Nup[i]*fc.Ndn[i] for i in range(n))
    fc.set_hamiltonian(H)
    fc.set_conserved_sector(Nf=4)
    assert groundstate.native_strings(fc)
    assert fc.gs_energy() == pytest.approx(fc.gs_energy(mode="ED"),
                                           abs=1e-8)


def test_short_range_hamiltonian_is_sent_as_before():
    """Below the threshold the terms are exactly the explicit ones, so the
    MPO of every short-range Hamiltonian is built as before."""
    n = 10
    fc = fermionchain.Spinful_Fermionic_Chain(n, itensor_version=3)
    H = 0
    for i in range(n-1):
        H = H - (fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
        H = H - (fc.Cdagup[i+1]*fc.Cup[i] + fc.Cdagdn[i+1]*fc.Cdn[i])
    H = H + 4.0*sum(fc.Nup[i]*fc.Ndn[i] for i in range(n))
    fc.set_hamiltonian(H)
    assert not groundstate.native_strings(fc)
    assert groundstate._send_key(fc)[1] == H.to_terms()


def test_gaussian_basis_hamiltonian_takes_the_native_route():
    n = 32
    fc = fermionchain.Spinful_Fermionic_Chain(n, itensor_version=3)
    h = np.zeros((n, n))
    for i in range(n-1):
        h[i, i+1] = h[i+1, i] = -1.0
    H = sum(h[i, j]*(fc.Cdagup[i]*fc.Cup[j]+fc.Cdagdn[i]*fc.Cdn[j])
            for i in range(n) for j in range(n) if h[i, j] != 0)
    H = H + 4.0*(fc.Nup[0]-0.5)*(fc.Ndn[0]-0.5)
    b = fc.get_gaussian_basis(h=h, nf=n, fixed=[0])
    fc.set_hamiltonian(b.rotate(H))
    assert groundstate.native_strings(fc)
    fc2 = fermionchain.Spinful_Fermionic_Chain(n, itensor_version="python")
    fc2.set_hamiltonian(b.rotate(H))
    assert not groundstate.native_strings(fc2)  # v3 only
