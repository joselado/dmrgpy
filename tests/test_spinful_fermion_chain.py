"""Functional coverage for fermionchain.Spinful_Fermionic_Chain (the
Hubbard-model chain), distilled from examples/hubbard_chain and
examples/fermionic_static_correlator down to small, fast systems,
comparing ED against DMRG on both ITensor v2 and v3 (see _helpers.py)."""
import numpy as np
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
    regression values for each site j.

    The plain Hamiltonian is symmetric under swapping up and down spins,
    so its ground state is exactly 2-fold degenerate (confirmed via ED:
    the two lowest eigenvalues coincide to machine precision) -- an
    unconstrained DMRG search (no fixed Sz/N sector, see CLAUDE.md) can
    converge to either member of that degenerate subspace, and <Nup_0>
    is *not* symmetric under the up/down swap, so it differs between
    them. This was only invisible while every "DMRG" call in this repo's
    dev environment silently fell back to ED (see mode.py): running the
    literal same test against a real compiled DMRG backend converged to
    the right energy but a correlator wildly off from the ED value here.
    A small explicit Zeeman-like term breaks the degeneracy so the
    ground state -- and this observable -- are uniquely defined."""
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
    eps = 0.3  # lifts the up/down degeneracy
    h = h + eps * (fc.Nup[0] - fc.Ndn[0])

    expected = [0.09204667896906626, -0.20764359692352258, 0.09821625446016191]
    for j in range(n):
        op = fc.Cdagup[0] * fc.Cup[j]
        v_ed, v_v2, v_v3 = vev_ed_v2_v3(fc, h, op)
        assert v_ed.real == pytest.approx(expected[j], abs=1e-6)
        assert v_v2.real == pytest.approx(expected[j], abs=DMRG_TOL)
        assert v_v3.real == pytest.approx(expected[j], abs=DMRG_TOL)


def test_long_range_hopping_correlator_matches_free_fermion_exact():
    """<Cdagup_i0 Cup_j> across a long separation, against the exact
    free-fermion one-body density matrix rather than a golden value.

    Two things the golden-value test above cannot check. First, the oracle
    here is independent of dmrgpy entirely: at U=0 the ground state is a
    Slater determinant, so <c^dag_i c_j> = sum over occupied single-particle
    states of conj(psi_n(i)) psi_n(j), computed here with plain numpy.
    Second, and the actual point: `Spinful_Fermionic_Chain` represents each
    orbital as two interleaved spinless-fermion sites, so <Cdagup_i0
    Cup_{i0+5}> spans NINE intervening tensor-network sites, every one of
    which must carry a Jordan-Wigner F factor. A missing or mis-signed
    string is invisible at the short separations a 3-site chain can reach,
    and would show up not as noise but as a sign flip or a systematic drift
    growing with distance -- i.e. exactly as a wrong exponential decay rate,
    the quantity such a correlator is usually measured for.

    The chain is dimerized (t1 != t2) so it is gapped and the occupied set
    is unambiguous: there is no chemical-potential term, so DMRG fills
    precisely the negative-energy single-particle levels.
    """
    n, t1, t2, i0 = 8, 1.0, 0.35, 2

    h1 = np.zeros((n, n))
    for i in range(n - 1):
        h1[i, i + 1] = h1[i + 1, i] = t1 if i % 2 == 0 else t2
    w, v = np.linalg.eigh(h1)
    occ = v[:, w < 0]
    P = occ.conj() @ occ.T          # P[i,j] = <c^dag_i c_j>, per spin flavour

    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        t = t1 if i % 2 == 0 else t2
        h = h + t * (fc.Cdagup[i] * fc.Cup[i + 1] + fc.Cdagup[i + 1] * fc.Cup[i])
        h = h + t * (fc.Cdagdn[i] * fc.Cdn[i + 1] + fc.Cdagdn[i + 1] * fc.Cdn[i])
    fc.set_hamiltonian(h)
    fc.maxm = 60

    assert fc.gs_energy() == pytest.approx(2 * w[w < 0].sum(), abs=1e-5)
    for d in range(n - i0):
        got = complex(fc.vev(fc.Cdagup[i0] * fc.Cup[i0 + d]))
        assert got.real == pytest.approx(P[i0, i0 + d], abs=1e-5), \
            "separation d={} disagrees with the exact one-body density matrix".format(d)
