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


def test_submode_nrg_exact_head_is_the_exact_start_shell():
    """get_stm_spectrum(submode="NRG", head="exact") on a kondo_chain is
    chain_poles with the chain's own parameters, continued past its end."""
    J = 0.4*ysr.lambda_factor(4.0)
    fc = ysr.kondo_chain(8, J=J, Lambda=4.0, z=0.3, itensor_version="python")
    es = np.array([0.0, 1e-4, 1e-2, -0.3])
    _, A, rs = ysr.get_stm_spectrum(fc, es, submode="NRG", head="exact",
                                    keep=200, return_poles=True)
    eps, t = nrgmps.tail_chain(fc)
    r = nrgspectral.chain_poles(ns=1, n=len(eps), JK=J, Lambda=4.0,
                                wilson=(eps, t), keep=200, stop_tol=1e-4,
                                stop_extra=6)
    assert rs[0].nshells == r.nshells
    ref = nrgspectral.curve(es, r.omega, r.weight, omega_N=r.omega_n[-1],
                            b=0.8)
    assert np.max(np.abs(A - ref)) < 1e-12*np.max(np.abs(ref))


def test_submode_nrg_refusals():
    fc = ysr.kondo_chain(6, J=0.4, Lambda=4.0, itensor_version="python")
    es = np.array([0.1])
    for kw in (dict(imp=0), dict(site=1), dict(delta=0.1), dict(mode="ED"),
               dict(head="ED")):
        with pytest.raises(ValueError):
            ysr.get_stm_spectrum(fc, es, submode="NRG", **kw)
    sc = ysr.kondo_chain(6, J=0.4, Delta=0.1, Lambda=4.0,
                         itensor_version="python")
    with pytest.raises(ValueError):
        ysr.get_stm_spectrum(sc, es, submode="NRG", head="exact")


def test_twin_chain_is_the_same_model():
    """The twist of fc itself rebuilt by twin_chain has the ED ground
    energy of fc, for a spin chain on a metal and for a kondo_chain."""
    fc = ysr.spin_chain_on_metal(2, 3, site=0, Js=0.3, JK=0.5, B=0.01,
                                 Lambda=4.0, z=0.25, itensor_version="python")
    c = nrgmps.twin_chain(fc, 0.25)
    assert c.gs_energy(mode="ED") == pytest.approx(fc.gs_energy(mode="ED"),
                                                   abs=1e-12)
    kc = ysr.kondo_chain(3, J=0.4, B=0.02, Lambda=4.0, z=0.5,
                         itensor_version="python")
    c = nrgmps.twin_chain(kc, 0.5)
    assert c.gs_energy(mode="ED") == pytest.approx(kc.gs_energy(mode="ED"),
                                                   abs=1e-12)


def test_head_tail_join_carries_the_sum_rule():
    """With a complete hand-over basis the tail misses nothing (roundoff).
    With the basis of the DMRG state it misses the part of X^dag|GS>
    outside it, at high energy (1.4% here, half of the weight on a head of
    10 sites with 59 states), and the tail poles below w_join with the TASK
    poles of the head above it add up to the sum rule of the head state to
    2% (one spin, a head of 7 Wilson sites, measured 1.1% at nkry = 80):
    nothing conserves the weight across the cut."""
    J = 0.4*ysr.lambda_factor(4.0)
    fc = ysr.spin_chain_on_metal(1, 7, JK=J, Lambda=4.0,
                                 itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 20
    r = nrgmps.head_tail_poles(fc, m=4, keep=200, nkry=80, pad_bonds="full",
                               task=False)
    assert r.k0 == 2 and r.chi == 4**4
    assert abs(r.missing_weight) < 1e-12
    r = nrgmps.head_tail_poles(fc, m=4, keep=200, nkry=80)
    assert r.tail.total == pytest.approx(r.tail.sum_rule - r.missing_weight,
                                         abs=1e-12)
    assert r.missing_weight > 1e-3*r.sum_rule
    assert abs(np.sum(r.weight) - r.sum_rule) < 0.02*r.sum_rule


def test_real_gauge_with_tiny_schmidt_values():
    """A 6-site "python" state whose smallest Schmidt values are 1e-6: the
    weighted real gauge gives a residual at the accuracy of the solve, where
    the unweighted one found 3.7e-4 (independent review), and the hand-over
    works."""
    c, e0 = _head(1, 6, "python")
    ex = nrgmps.Export(c)
    assert ex.real_residual < 1e-5
    sh, _ = nrgmps.handover_shell(c, 0, export=ex, Lambda=4.0)
    assert sh.info["chi"] == 8 # the spin and the Wilson site 0


@pytest.mark.parametrize("version", VERSIONS)
def test_head_state_takes_the_lower_zeeman_member(version):
    """With B > 0 on the spins and an odd number of electrons the head is
    the 2Sz = -1 member, the ground state, and not the 2Sz = +1 one."""
    eps, t = _wilson(6)
    fc = ysr.spin_chain_on_metal(1, 6, JK=JK, Js=JS, B=0.05, Lambda=4.0,
                                 wilson=(eps, t), itensor_version=version)
    fc.maxm, fc.nsweeps = 64, 20
    c, e0 = nrgmps.head_state(fc)
    assert c.head_sector == dict(Nf=7, Sz=-1)
    assert e0 == pytest.approx(fc.gs_energy(mode="ED"), abs=1e-8)


def test_tail_chain_stays_above_the_broken_end():
    """The tail chain ends before its scale falls below 1e-30, where
    wilson_chain loses its precision (99 sites at Lambda = 4)."""
    assert nrgmps.max_sites(4.0) == 99
    fc = ysr.kondo_chain(6, J=0.4, Lambda=4.0, itensor_version="python")
    eps, t = nrgmps.tail_chain(fc)
    assert len(eps) == 99
    assert np.max(np.abs(eps)/np.abs(np.append(t, t[-1]))) < 1.0
    assert abs(t[-1]) < 1e-28
