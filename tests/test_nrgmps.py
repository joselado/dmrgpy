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
    10 sites with 59 states), and the tail poles joined with the TASK
    poles of the head (by default the fade of join_share around w_join)
    add up to the sum rule of the head state to 2% (one spin, a head of 7
    Wilson sites, measured 1.1% at nkry = 80 with the cut): nothing
    conserves the weight across the join."""
    J = 0.4*ysr.lambda_factor(4.0)
    fc = ysr.spin_chain_on_metal(1, 7, JK=J, Lambda=4.0,
                                 itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 20
    r = nrgmps.head_tail_poles(fc, m=4, keep=200, nkry=80, pad_bonds="full",
                               task=False)
    assert r.k0 == 2 and r.chi == 4**4
    assert abs(r.missing_weight) < 1e-12
    r = nrgmps.head_tail_poles(fc, m=4, keep=200, nkry=80)
    assert r.join == "fade"
    assert r.tail.total == pytest.approx(r.tail.sum_rule - r.missing_weight,
                                         abs=1e-12)
    assert r.missing_weight > 1e-3*r.sum_rule
    assert abs(np.sum(r.weight) - r.sum_rule) < 0.02*r.sum_rule


def test_join_share_is_a_step_or_a_cosine_in_log_w():
    """The share of a tail pole: a step at w_join for "cut", and for
    "fade" 1 below w_join/window, 0 above w_join*window, 1/2 at w_join,
    monotone, the same on both signs of w; a TASK pole takes the rest."""
    w = np.concatenate([-np.logspace(-6, 1, 141), np.logspace(-6, 1, 141)])
    cut = nrgmps.join_share(w, 1e-2, 4.0, "cut")
    assert np.array_equal(cut, (np.abs(w) < 1e-2).astype(float))
    f = nrgmps.join_share(w, 1e-2, 4.0, "fade")
    a = np.abs(w)
    assert np.all(f[a <= 1e-2/4] == 1.0) and np.all(f[a >= 4e-2] == 0.0)
    assert nrgmps.join_share([1e-2, -1e-2], 1e-2, 4.0, "fade") == \
        pytest.approx([0.5, 0.5], abs=1e-15)
    pos = w > 0
    assert np.all(np.diff(f[pos]) <= 0.0)
    assert np.array_equal(f[pos], f[~pos]) # -w runs over the same |w|
    with pytest.raises(ValueError):
        nrgmps.join_share(w, 1e-2, 4.0, "smooth")
    with pytest.raises(ValueError):
        nrgmps.join_share(w, 1e-2, 1.0, "fade")


def test_cut_and_fade_joins_of_tail_and_task():
    """On one head state: join="cut" is the tail poles below w_join and
    the TASK poles from it on, exactly as before join existed; "fade" is
    every pole of both routes times its share, the tail alone below
    w_join/window and TASK alone above w_join*window, and it carries the
    sum rule as the cut does (measured on 3 spins: the cut 0.02 to 0.3%
    off, the fade 0.2 to 0.4%, README section 6)."""
    J = 0.4*ysr.lambda_factor(4.0)
    fc = ysr.spin_chain_on_metal(1, 7, JK=J, Lambda=4.0,
                                 itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 20
    kw = dict(m=4, keep=200, nkry=80)
    rc = nrgmps.head_tail_poles(fc, join="cut", **kw)
    c = rc.chain
    rf = nrgmps.head_tail_poles(c, solved=True, join="fade", **kw)
    t, p, w, s = rc.tail, rc.task, rc.w_join, rc.window
    assert rc.join == "cut" and rf.join == "fade" and s[0] == w/4.0
    below, above = np.abs(t.omega) < w, np.abs(p["omega"]) >= w
    assert np.array_equal(rc.omega, np.concatenate([t.omega[below],
                                                    p["omega"][above]]))
    assert np.array_equal(rc.weight, np.concatenate([t.weight[below],
                                                     p["weight"][above]]))
    np.testing.assert_allclose(rf.task["weight"], p["weight"], atol=1e-14)
    f = nrgmps.join_share(rf.tail.omega, w, 4.0, "fade")
    g = 1.0 - nrgmps.join_share(rf.task["omega"], w, 4.0, "fade")
    total = np.sum(rf.tail.weight*f) + np.sum(rf.task["weight"]*g)
    assert np.sum(rf.weight) == pytest.approx(total, rel=1e-12)
    a = np.abs(rf.omega)
    low = np.sort(rf.weight[a < w/4.0])
    assert np.array_equal(low, np.sort(rf.tail.weight[np.abs(rf.tail.omega)
                                                      < w/4.0]))
    high = np.sort(rf.weight[a > 4.0*w])
    assert np.array_equal(high, np.sort(rf.task["weight"][
        np.abs(rf.task["omega"]) > 4.0*w]))
    assert abs(np.sum(rf.weight) - rf.sum_rule) < 0.02*rf.sum_rule
    with pytest.raises(ValueError):
        nrgmps.head_tail_poles(c, solved=True, join="smooth", **kw)


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


def test_a_small_handover_basis_warns():
    """Two sites after the hand-over hold at most 16 states, far below the
    about 100 the tail needs, and head_tail_poles says so."""
    J = 0.4*ysr.lambda_factor(4.0)
    fc = ysr.spin_chain_on_metal(1, 6, JK=J, Lambda=4.0,
                                 itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 10
    with pytest.warns(RuntimeWarning, match="hand-over basis"):
        nrgmps.head_tail_poles(fc, m=2, keep=100, task=False, n_tail=20)


def _doublet_head(version):
    """3 spins and 4 Wilson sites, the metal under the middle spin: an odd
    number of electrons, so the head is one member (2Sz = 1) of a
    doublet."""
    eps, t = _wilson(4)
    fc = ysr.spin_chain_on_metal(3, 4, site=1, JK=JK, Js=JS, Lambda=4.0,
                                 wilson=(eps, t), itensor_version=version)
    fc.maxm, fc.nsweeps = 64, 20
    c, e0 = nrgmps.head_state(fc)
    return c, eps, t


@pytest.mark.parametrize("version", VERSIONS)
def test_the_spin_flip_gives_the_partner_member(version):
    """R = exp(-i pi S^y_tot), site by site (spin_flip_matrix), maps the
    head onto the other member of its doublet: the same energy, the
    opposite S^z, and no overlap with the state."""
    c, eps, t = _doublet_head(version)
    ex = nrgmps.Export(c)
    exr = ex.flipped(nrgmps.spin_flip_matrix(c, ex))
    last = ex.L - 1
    Sz = sum(c.Sz[i] for i in range(ex.L))
    H = c.hamiltonian
    assert exr.block_matrix(H, last)[0, 0] == \
        pytest.approx(ex.block_matrix(H, last)[0, 0], abs=1e-12)
    assert ex.block_matrix(Sz, last)[0, 0] == pytest.approx(0.5, abs=1e-10)
    assert exr.block_matrix(Sz, last)[0, 0] == pytest.approx(-0.5, abs=1e-10)
    assert abs(ex.overlap(exr, last)[0, 0]) < 1e-12


@pytest.mark.parametrize("version", VERSIONS)
def test_a_symmetric_handover_is_closed_under_the_flip(version):
    """The left basis of one member of a doublet is not closed under the
    flip (partner blocks of different sizes); the joined basis of both
    members is, with the spectra of (Q, 2Sz) and (Q, -2Sz) equal. With a
    complete basis it is still the exact NRG."""
    c, eps, t = _doublet_head(version)
    def asymmetric(sh):
        return [(q, s) for (q, s), e in sh.E.items() if (q, -s) not in sh.E
                or len(sh.E[(q, -s)]) != len(e)
                or np.max(np.abs(sh.E[(q, -s)] - e)) > 1e-9]
    plain, _ = nrgmps.handover_shell(c, 1, U=0.1, Lambda=4.0)
    sym, _ = nrgmps.handover_shell(c, 1, U=0.1, Lambda=4.0, symmetric=True)
    assert asymmetric(plain)
    assert asymmetric(sym) == []
    assert sym.info["chi"] > plain.info["chi"]
    full, _ = nrgmps.handover_shell(c, 0, U=0.1, Lambda=4.0, symmetric=True,
                                    pad_bonds="full")
    r = nrgspectral.spectral_poles(full, eps[1:], t, keep=None,
                                   bins_per_decade=None)
    rr = nrgspectral.chain_poles(ns=3, n=4, site=1, JK=JK, Js=JS, U=0.1,
                                 wilson=(eps, t), keep=None,
                                 bins_per_decade=None, Lambda=4.0)
    assert np.max(np.abs(resolvent(r.omega, r.weight)
                         - resolvent(rr.omega, rr.weight))) < 1e-11


@pytest.mark.parametrize("version", VERSIONS)
def test_auto_joins_the_flipped_basis_only_for_a_head_with_a_moment(version):
    """symmetric="auto" (the default) takes the joined basis for a head
    with an odd number of electrons (3 spins and 4 Wilson sites, 2Sz = 1)
    and the plain one for an even head (5 Wilson sites, 2Sz = 0)."""
    kw = dict(k0=1, keep=100, task=False, n_tail=20, solved=True)
    for n, moment in ((4, True), (5, False)):
        fc = ysr.spin_chain_on_metal(3, n, site=1, JK=JK, Js=JS, Lambda=4.0,
                                     itensor_version=version)
        fc.maxm, fc.nsweeps = 64, 20
        c, e0 = nrgmps.head_state(fc)
        assert (c.head_sector["Sz"] != 0) == moment
        chi = {s: nrgmps.head_tail_poles(c, symmetric=s, **kw).chi
               for s in ("auto", True, False)}
        assert chi["auto"] == (chi[True] if moment else chi[False])
        if moment:
            assert chi[True] > chi[False]
    with pytest.raises(ValueError, match="symmetric"):
        nrgmps.head_tail_poles(c, symmetric="yes", **kw)


def test_a_symmetric_handover_needs_zero_field():
    eps, t = _wilson(4)
    fc = ysr.spin_chain_on_metal(1, 4, JK=JK, Lambda=4.0, B=1e-3,
                                 wilson=(eps, t), itensor_version="python")
    fc.maxm, fc.nsweeps = 32, 10
    c, e0 = nrgmps.head_state(fc)
    with pytest.raises(ValueError, match="B = 0"):
        nrgmps.handover_shell(c, 1, Lambda=4.0, symmetric=True)


@pytest.mark.parametrize("version", VERSIONS)
def test_the_mpo_keeps_the_deep_wilson_bonds(version):
    """The depth floor of README section 9 is gone: on 30 Wilson sites at
    Lambda = 4 the MPO used to drop the channels of the hoppings from
    t_21 = 3.0e-7 on (v3's svdMPO cut a squared singular value below an
    absolute 1e-13) and from t_20 on "python"; since the cutoff of the
    build is one per bond (mo_terms.h bond_cutoff_args, mpobuilder.py
    _bond_cutoff_factor) every bond keeps its channels. What is left is
    the Hamiltonian's own clean threshold, 1e-12 of its largest
    coefficient, which removes the hoppings from H itself (past the Wilson
    site 37 on 50 sites)."""
    J = 0.3*ysr.lambda_factor(4.0)
    deep = ysr.spin_chain_on_metal(1, 30, Js=0.01, JK=J, Lambda=4.0,
                                   itensor_version=version)
    assert nrgmps.lost_wilson_bonds(deep) == []
    from dmrgpy.pyitensor import mpobuilder
    ampo = deep._session._ampo(deep.hamiltonian.to_terms()) \
        if version == "python" else None
    if ampo is not None:
        f = mpobuilder._bond_cutoff_factor(ampo)
        assert 0.0 < f < 1e-12


def test_a_python_sector_solve_keeps_the_deep_wilson_hoppings():
    """On "python" a sector solve minimizes H plus a charge penalty, summed
    by mps_sum. Cut at a relative 1e-14, that sum lost the deep Wilson
    hoppings next to the penalty's constant lambda*q^2 (about 1.3e4 on 20
    Wilson sites): <H + P> on the solved state was 2.4e-4 above
    <H> + <P>. Chain._penalized_cutoff lowers the cutoff by the per-bond
    factor of H and by (cmax/scale)^2; what is left, 5.6e-7, is the same
    with no cutoff at all, the rounding of the sum next to that constant."""
    from dmrgpy.pyitensor.mpsalgebra import inner
    J = 0.3*ysr.lambda_factor(4.0)
    n = 20
    fc = ysr.spin_chain_on_metal(1, n, Js=0.01, JK=J, Lambda=4.0,
                                 itensor_version="python")
    fc.maxm, fc.nsweeps = 16, 4
    fc.set_conserved_sector(Nf=n + 1, Sz=1)
    fc.gs_energy()
    s = fc._session
    psi = s.wf0
    E = lambda m: complex(inner(psi, m, psi)).real
    assert abs(E(s._solve_H()) - E(s.H) - E(s._penalty_mpo())) < 1e-5


def test_bond_cutoff_factor_reads_the_weakest_bond():
    """m^2 of the per-bond rule: 1 for a uniform chain, s^2 for a weak link
    of strength s (merged over duplicate terms), and a one-site term sets
    the largest coefficient but crosses no bond."""
    from dmrgpy import spinchain
    from dmrgpy.pyitensor import mpobuilder
    sc = spinchain.Spin_Chain(["S=1/2"]*5, itensor_version="python")
    S = lambda i: sc.Sz[i]*sc.Sz[i+1]
    h = sum(S(i) for i in range(4))
    f = lambda op: mpobuilder._bond_cutoff_factor(
        sc._session._ampo(op.to_terms()))
    assert f(h) == 1.0
    assert f(h + 3.0*sc.Sz[0]) == pytest.approx(1/9.0, rel=1e-12)
    weak = h - S(2) + 1e-6*S(2)
    assert f(weak) == pytest.approx(1e-12, rel=1e-6)


def test_a_handover_past_a_lost_bond_warns(monkeypatch):
    """A hand-over whose block holds a Wilson bond the MPO dropped is not
    the block of the tail, and head_tail_poles says so; one at or above
    the lost bond does not."""
    J = 0.4*ysr.lambda_factor(4.0)
    fc = ysr.spin_chain_on_metal(1, 6, JK=J, Lambda=4.0,
                                 itensor_version="python")
    fc.maxm, fc.nsweeps = 64, 10
    monkeypatch.setattr(nrgmps, "lost_wilson_bonds", lambda c: [1, 2, 3, 4])
    with pytest.warns(RuntimeWarning, match="dropped channels of the Wilson bond 1"):
        nrgmps.head_tail_poles(fc, m=2, keep=100, task=False, n_tail=20)
    monkeypatch.setattr(nrgmps, "lost_wilson_bonds", lambda c: [3, 4])
    with pytest.warns(RuntimeWarning) as record:
        nrgmps.head_tail_poles(fc, m=2, keep=100, task=False, n_tail=20)
    assert not any("dropped channels" in str(w.message) for w in record)
