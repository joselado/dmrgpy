"""Tests of nrgmps: the shells of an NRG read off the left bases of a DMRG
state (N2) and the hand-over from a DMRG head to the NRG tail (N3) of
docs/nrg_mps_kondo_plan.md."""

import numpy as np
import pytest

from dmrgpy import nrgmps, nrgspectral, wilsonchain, ysr
from dmrgpy import cppext

ZS = np.array([0.3+0.05j, -0.2+0.01j, 0.01+0.003j, 1.5+0.2j, -0.05+0.02j,
               0.002+0.001j])
VERSIONS = ["python"] + ([3] if cppext.available(3) else [])
JK, JS = 0.4, 0.3


def resolvent(w, W):
    return np.array([np.sum(W/(z - w)) for z in ZS])


def _wilson(n):
    return wilsonchain.wilson_chain(n, Lambda=4.0, band=(-1.0, 1.0), z=0.25)


def _head(ns, n, version):
    eps, t = _wilson(n)
    fc = ysr.spin_chain_on_metal(ns, n, JK=JK, Js=JS, Lambda=4.0,
                                 wilson=(eps, t), itensor_version=version)
    fc.maxm, fc.nsweeps = 64, 20
    return nrgmps.head_state(fc)


@pytest.mark.parametrize("version", VERSIONS)
def test_flow_diagram_is_exact_where_the_left_basis_is_complete(version):
    """One spin and 5 Wilson sites: the left basis of the Wilson sites 0
    and 1 holds the whole physical block (8 and 32 states), so the block
    spectrum in it is that of the exact NRG shell (4e-15); at site 2 only
    the 16 states of the right part fit, and the levels are off by 0.34
    omega_2."""
    c, e0 = _head(1, 5, version)
    eps, t = _wilson(5)
    start = nrgspectral.spin_cluster_shell(ns=1, JK=JK, Js=JS, eps0=eps[0],
                                           omega=2*abs(t[0]))
    shells = nrgspectral.iterate(start, eps[1:], t, keep=None)
    fd = nrgmps.flow_diagram(c, Lambda=4.0)
    for k in (0, 1):
        ex = np.sort(np.concatenate(list(shells[k].e.values())))
        ex = ex/shells[k].omega
        assert fd[k]["chi"] == 8*4**k
        assert np.max(np.abs(fd[k]["levels"][:20] - ex[:20])) < 1e-12
        assert fd[k]["E0"] == pytest.approx(shells[k].e0, abs=1e-12)
    ex2 = np.sort(np.concatenate(list(shells[2].e.values())))/shells[2].omega
    assert fd[2]["chi"] == 16
    assert np.max(np.abs(fd[2]["levels"] - ex2[:16])) > 0.1


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("k0,T", [(0, 0.0), (1, 0.0), (1, 0.05)])
def test_handover_with_a_complete_basis_equals_the_full_nrg(version, k0, T):
    """pad_bonds="full" completes the left basis at the hand-over to the
    whole space of the block, so the DMRG head plus the NRG tail is the
    exact NRG of the same chain, whatever the accuracy of the solve
    (measured 5e-14 on the resolvent, at T = 0 and T > 0)."""
    c, e0 = _head(1, 5, version)
    eps, t = _wilson(5)
    sh, ex = nrgmps.handover_shell(c, k0, U=0.2, pad_bonds="full",
                                   Lambda=4.0)
    r = nrgspectral.spectral_poles(sh, eps[k0+1:], t[k0:], keep=None, T=T,
                                   bins_per_decade=None)
    rr = nrgspectral.chain_poles(ns=1, n=5, JK=JK, Js=JS, U=0.2,
                                 wilson=(eps, t), keep=None, T=T,
                                 bins_per_decade=None, Lambda=4.0)
    assert np.max(np.abs(resolvent(r.omega, r.weight)
                         - resolvent(rr.omega, rr.weight))) < 1e-11
    assert r.sum_rule == pytest.approx(rr.total, abs=1e-12)
    assert r.total == pytest.approx(r.sum_rule, abs=1e-12)


@pytest.mark.parametrize("version", VERSIONS)
def test_two_spins_handover_with_a_complete_basis(version):
    """The same with two spins and the bath on the first one, a
    long-range Kondo term: the hand-over at the Wilson site 0."""
    eps, t = _wilson(4)
    fc = ysr.spin_chain_on_metal(2, 4, site=0, JK=JK, Js=JS, Lambda=4.0,
                                 wilson=(eps, t), itensor_version=version)
    fc.maxm, fc.nsweeps = 64, 20
    c, e0 = nrgmps.head_state(fc)
    sh, ex = nrgmps.handover_shell(c, 0, U=0.1, pad_bonds="full", Lambda=4.0)
    r = nrgspectral.spectral_poles(sh, eps[1:], t, keep=None,
                                   bins_per_decade=None)
    rr = nrgspectral.chain_poles(ns=2, n=4, site=0, JK=JK, Js=JS, U=0.1,
                                 wilson=(eps, t), keep=None,
                                 bins_per_decade=None, Lambda=4.0)
    assert np.max(np.abs(resolvent(r.omega, r.weight)
                         - resolvent(rr.omega, rr.weight))) < 1e-11


def test_real_gauge_of_a_complex_python_state():
    """The "python" backend returns a state that is real up to a phase in
    a complex gauge; the export puts it in a real one, and the residual
    of the conjugation test is at roundoff."""
    c, e0 = _head(1, 5, "python")
    ex = nrgmps.Export(c)
    assert all(not np.iscomplexobj(A) for A in ex.As)
    assert ex.real_residual < 1e-10


def test_graded_basis_refuses_a_mixed_basis():
    N = np.diag([1.0, 2.0])
    S = np.diag([1.0, 0.0])
    V, lab = nrgmps.graded_basis(N, S)
    assert sorted(lab) == [(1, 2), (2, 0)]
    R = np.array([[np.cos(0.3), -np.sin(0.3)], [np.sin(0.3), np.cos(0.3)]])
    P = R[:, :1] # one vector mixing the two charges
    with pytest.raises(RuntimeError):
        nrgmps.graded_basis(P.T @ N @ P, P.T @ S @ P)


def test_block_matrix_refuses_an_operator_outside_the_block():
    c, e0 = _head(1, 5, "python")
    ex = nrgmps.Export(c)
    with pytest.raises(ValueError):
        ex.block_matrix(c.Nup[4], 2)
