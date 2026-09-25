"""Regression tests for the `julia` cluster of the 2026-09-25b hole hunt
(docs/audit_2026_09_25b_hole_hunt.md), all on itensor_version="julia_live".

  #7  Every Julia DMRG-family solve started from a bond-dimension-1
      product state (ITensorMPS random_mps(sites)), and the Hermitian ones
      got no noise, so a Hamiltonian coupling two sites across a site that
      carries no term stopped at a product state: gs_energy() -0.5 against
      -1.0 on the Heisenberg chain on the even sites of 6, -1.5 against
      -3.232051 on an 8-site J2-only chain, Thermal_Spin_Chain at T=0 <H>
      -0.5 against -1.0, and the generalized and non-Hermitian solves with
      them. Each solve now starts from a random MPS of link dimension
      min(maxm, bond_ramp_start), and the Hermitian ones get the chain's
      noise on the first half of the schedule.
  #11 densitymatrix.jl divided the state by <psi|psi> instead of its
      norm, so after set_gs(c*s) get_rdm was rho/c^2.
  #12 get_rdm(i=ns-1) raised BoundsError (the Julia half of 2026-08 #16).
  #23 kpm.jl's same_mps took ||vi-vj|| < 1e-10 absolutely, so a cross pair
      of small images took the auto-correlator recursion: C[eps*Sz0,
      eps*Sz3] came back as C[Sz3,Sz3].
  lead session-julia-ed-guard-carveout-rdm-bond-entropy: get_rdm and the
      bond entropy skipped the "no ED implementation" guard on julia_live.

Anchors are ED on the same Hamiltonian, the generalized eigenproblem from
scipy, and exact identities (the reduced density matrix of a ray, the
bilinearity of the correlator). julia_live pays its JIT once per process,
so the chains are small and the checks share their Julia signatures.
"""

import io
import contextlib
import warnings

import numpy as np
import pytest
import scipy.linalg as sla

from dmrgpy import cppext, spinchain, thermal

from _helpers import julia_available


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


JULIA = [_backend("julia_live")]


def _require_julia(version):
    """Deferred to the test body, so that collecting this file (and
    -k "not julia_live") never boots a Julia session"""
    if version == "julia_live":
        ok, reason = julia_available()
        if not ok:
            pytest.skip("requires a working juliacall/Julia toolchain: %s"
                        % reason)


def _quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()


def _heis_pairs(sc, pairs):
    h = 0
    for (i, j) in pairs:
        h = h + sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j]
    return h


# two chains whose couplings all skip a site that carries no term: the
# 3-site Heisenberg chain on the even sites of 6 (three free spins beside
# it), and two decoupled 4-site chains interleaved on 8 sites
EVEN6 = ([(0, 2), (2, 4)], 6)
J2_8 = ([(i, i + 2) for i in range(6)], 8)


def _ed_matrix(pairs, n, extra=None):
    ed = spinchain.Spin_Chain(["S=1/2"]*n)
    h = _heis_pairs(ed, pairs)
    if extra is not None:
        h = h + extra(ed)
    ed.set_hamiltonian(h)
    return ed, np.asarray(ed.get_ED_obj().get_hamiltonian().todense())


def _chain(version, pairs, n, extra=None):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    h = _heis_pairs(sc, pairs)
    if extra is not None:
        h = h + extra(sc)
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 20
    return sc


# ----------------------------------------------------------------- #7

@pytest.mark.parametrize("version", JULIA + [_backend("python")])
def test_decoupled_sublattices_reach_the_ground_state(version):
    """gs_energy() on both decoupled-sublattice chains against ED, at the
    chain's own noise and at noise=0, where only the start can rescue the
    solve. julia_live gave -0.5 and -1.5 at every noise, from a state of
    link dimension 1; "python" is the backend that always got them."""
    _require_julia(version)
    for pairs, n in (EVEN6, J2_8):
        _, H = _ed_matrix(pairs, n)
        e_ed = np.linalg.eigvalsh(H)[0]
        for noise in (None, 0.0):
            sc = _chain(version, pairs, n)
            if noise is not None:
                sc.noise = noise
            assert np.real(_quiet(sc.gs_energy)) == pytest.approx(e_ed, abs=1e-7)
            if version == "julia_live":
                from dmrgpy.mpsjulialive.juliasession import Main as Mainjl
                assert int(Mainjl.maxlinkdim(sc.get_gs().jlmps)) > 1


@pytest.mark.parametrize("version", JULIA)
def test_the_chain_noise_reaches_the_julia_solve(version):
    """A sweep from the caller's product-state guess has no bond to grow
    through, so on the even-site chain only the noise gets it off the
    product state; julia_live never forwarded self.noise, and the guess
    stayed at -0.5 whatever it was set to. The noise=0 control is the
    discriminant: the same guess, swept without noise, does not reach the
    ground state."""
    _require_julia(version)
    pairs, n = EVEN6
    results = {}
    for noise in (1e-7, 0.0):
        sc = _chain(version, pairs, n)
        sc.noise = noise
        guess = sc.random_state()  # a product state (mpsalgebra.jl's random_state)
        sc.set_initial_wf_guess(guess)
        results[noise] = np.real(_quiet(sc.gs_energy))
    assert results[1e-7] == pytest.approx(-1.0, abs=1e-7)
    assert results[0.0] > -0.9


@pytest.mark.parametrize("version", JULIA)
def test_make_sweeps_tapers_the_noise_like_the_session_backends(version):
    """First half of the schedule only (mpscpp2/mpscpp3's make_sweeps,
    pyitensor's _make_sweeps), and taper=false for the one-sweep schedules
    generalized.jl builds per outer iteration."""
    _require_julia(version)
    from dmrgpy.mpsjulialive.juliasession import Main as Mainjl
    noise_of = Mainjl.seval("s -> collect(s.noise)")
    for ns, want in ((1, [0.0]), (4, [1e-3, 1e-3, 0.0, 0.0]),
                     (5, [1e-3, 1e-3, 0.0, 0.0, 0.0])):
        got = list(noise_of(Mainjl.make_sweeps(ns, 30, 1e-12, noise=1e-3)))
        assert got == pytest.approx(want, abs=0.0)
    got = list(noise_of(Mainjl.make_sweeps(1, 30, 1e-12, noise=1e-3, taper=False)))
    assert got == pytest.approx([1e-3], abs=0.0)
    # the non-Hermitian callers pass no noise and keep a noise-free schedule
    got = list(noise_of(Mainjl.make_sweeps(4, 30, 1e-12)))
    assert got == pytest.approx([0.0]*4, abs=0.0)


@pytest.mark.parametrize("version", JULIA)
def test_thermal_chain_at_zero_temperature(version):
    """Thermal_Spin_Chain at T=0 is the ground manifold of the physical H,
    which lives on the even sites of the doubled chain: <H> = -1 exactly on
    3 sites, and <Sz0 Sz1> = -1/6 on every state of that manifold (the
    operator conserves the physical Sz and is -1/6 on both members of the
    doublet). julia_live gave <H> = -0.5 and a <Sz0 Sz1> that changed from
    run to run."""
    _require_julia(version)
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*3, T=0.0, itensor_version=version)
    tc.set_hamiltonian(_heis_pairs(tc, [(0, 1), (1, 2)]))
    tc.MBChain.maxm, tc.MBChain.nsweeps = 30, 20
    wf = _quiet(tc.get_gs)
    nrm = np.real(wf.dot(wf))
    e = np.real(wf.dot(tc.MBChain.hamiltonian*wf))/nrm
    zz = np.real(wf.dot((tc.Sz[0]*tc.Sz[1])*wf))/nrm
    assert e == pytest.approx(-1.0, abs=1e-7)
    assert zz == pytest.approx(-1.0/6.0, abs=1e-6)


@pytest.mark.parametrize("version", JULIA)
def test_generalized_and_non_hermitian_solves_leave_the_product_state(version):
    """The two other solves that started from the same product state:
    gs_energy_generalized(1+0.8*Sz0) on the even-site chain (-0.833333
    against the exact -1.496331) and the non-Hermitian gs_energy() on the
    J2-only chain plus 0.2j*Sz0 (-1.47+0.015j, not an eigenvalue, against
    -3.219401). NH-DMRG takes no noise on this backend, deliberately, so the
    start is all that rescues it."""
    _require_julia(version)
    pairs, n = EVEN6
    ed, H = _ed_matrix(pairs, n)
    A = np.asarray(ed.get_ED_obj().MO2matrix(1 + 0.8*ed.Sz[0]).todense())
    lam_ex = sla.eigh(H, A, eigvals_only=True)[0]
    sg = _chain(version, pairs, n)
    lam = np.real(_quiet(lambda: sg.gs_energy_generalized(1 + 0.8*sg.Sz[0])))
    assert lam == pytest.approx(lam_ex, abs=1e-6)
    pairs, n = J2_8
    nh = lambda c: 0.2j*c.Sz[0]
    _, H = _ed_matrix(pairs, n, extra=nh)
    w = np.linalg.eigvals(H)
    e_ex = w[np.argmin(w.real)]
    sn = _chain(version, pairs, n, extra=nh)
    e = complex(_quiet(sn.gs_energy))
    assert abs(e - e_ex) < 1e-4


# ---------------------------------------------------------- #11, #12

N4 = 4


def _ham4(sc):
    # unequal diagonal and a complex off-diagonal on every site, so that
    # neither rho/c^2 nor a transposed matrix can pass for rho
    return (_heis_pairs(sc, [(i, i + 1) for i in range(N4 - 1)])
            + 0.3*sc.Sz[0] + 0.2*sc.Sx[2] + 0.15*sc.Sy[3])


def _exact_rhos():
    ed = spinchain.Spin_Chain(["S=1/2"]*N4)
    ed.set_hamiltonian(_ham4(ed))
    out = []
    for i in range(N4):
        sx, sy, sz = [np.real(ed.vev(op[i], mode="ED")) for op in (ed.Sx, ed.Sy, ed.Sz)]
        # rho = 1/2 + 2 sum_a <S_a> S_a, in ITensor's (up, dn) order
        out.append(np.array([[0.5 + sz, sx - 1j*sy], [sx + 1j*sy, 0.5 - sz]]))
    return out


@pytest.mark.parametrize("version", JULIA + [_backend("python"), _backend(3)])
def test_get_rdm_is_the_density_matrix_of_the_ray_at_every_site(version):
    """Every site, the last one included (julia_live raised BoundsError
    there), matches the exact matrix, and after set_gs(c*s) and
    set_initial_wf(c*s) the matrix is still that of the ray: trace 1, not
    1/c^2 (0.25 at c=2 and 4.0 at c=0.5 on every backend)."""
    _require_julia(version)
    exact = _exact_rhos()
    sc = spinchain.Spin_Chain(["S=1/2"]*N4, itensor_version=version)
    sc.set_hamiltonian(_ham4(sc))
    sc.maxm, sc.nsweeps = 16, 12
    _quiet(sc.gs_energy)
    for i in range(N4):
        rho = np.asarray(sc.get_rdm(i=i))
        assert np.max(np.abs(rho - exact[i])) < 1e-6
    s = sc.get_gs().copy()
    s = s*(1.0/np.sqrt(np.real(s.dot(s))))
    for c in (2.0, 0.5):
        for setter in (sc.set_gs, sc.set_initial_wf):
            setter(s*c)
            for i in (1, N4 - 1):
                rho = np.asarray(sc.get_rdm(i=i))
                assert np.real(np.trace(rho)) == pytest.approx(1.0, abs=1e-10)
                assert np.max(np.abs(rho - exact[i])) < 1e-6


@pytest.mark.parametrize("version", JULIA)
def test_julia_reduced_dm_divides_by_the_norm(version):
    """The .jl formula itself, called on a scaled state directly, past the
    normalization get_rdm now does first: rho, not rho/c^2, at the first
    and the last site."""
    _require_julia(version)
    from dmrgpy.mpsjulialive import densitymatrix as dmjl
    exact = _exact_rhos()
    sc = spinchain.Spin_Chain(["S=1/2"]*N4, itensor_version=version)
    sc.set_hamiltonian(_ham4(sc))
    sc.maxm, sc.nsweeps = 16, 12
    _quiet(sc.gs_energy)
    s = sc.get_gs().copy()
    for c in (1.0, 3.0):
        for i in (0, N4 - 1):
            rho = dmjl.reduced_dm(sc, s*c, i)
            assert np.max(np.abs(rho - exact[i])) < 1e-6


# --------------------------------------------------------------- #23

@pytest.mark.parametrize("version", JULIA)
def test_julia_kpm_same_vector_test_is_relative(version):
    """kpm.jl's same_mps on vectors scaled to 1e-11: a vector is the same
    as itself and a cross pair is not, and two zero vectors take the full
    recursion (whose moments are zeros). Through the public correlator,
    C[eps*Sz0, eps*Sz3]/eps^2 at eps=1e-11 is C[Sz0,Sz3], which the
    absolute test answered with C[Sz3,Sz3], 2.07 of the peak off. The
    reference is the eps=1 call, the correlator being bilinear; run to run
    julia_live's band edge moves the curve by ~5e-4 of the peak."""
    _require_julia(version)
    from dmrgpy.mpsjulialive.juliasession import Main as Mainjl
    from dmrgpy.mpsjulialive.mpo import MPO
    L = 6
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=version)
    sc.set_hamiltonian(_heis_pairs(sc, [(i, i + 1) for i in range(L - 1)])
                       + 0.3*sc.Sz[0])
    sc.maxm, sc.nsweeps = 30, 10
    _quiet(sc.gs_energy)
    gs = sc.get_gs().jlmps
    eps = 1e-11
    a = Mainjl.apply_op(MPO(eps*sc.Sz[0], MBO=sc).jlmpo, gs, 30, 1e-12)
    b = Mainjl.apply_op(MPO(eps*sc.Sz[3], MBO=sc).jlmpo, gs, 30, 1e-12)
    assert bool(Mainjl.same_mps(a, a, 30, 1e-12))
    assert not bool(Mainjl.same_mps(a, b, 30, 1e-12))
    zero = Mainjl.mpstimesscalar(0.0, a)
    assert not bool(Mainjl.same_mps(zero, zero, 30, 1e-12))
    es = np.linspace(-0.5, 5.0, 200)
    kw = dict(delta=0.2, es=es)
    _, y_ab = _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[3]), **kw))
    _, y_bb = _quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[3], sc.Sz[3]), **kw))
    _, y = _quiet(lambda: sc.get_dynamical_correlator(
        name=(eps*sc.Sz[0], eps*sc.Sz[3]), **kw))
    y_ab, y_bb, y = np.asarray(y_ab), np.asarray(y_bb), np.asarray(y)/eps**2
    peak = np.max(np.abs(y_ab))
    assert np.max(np.abs(y_ab - y_bb))/peak > 1.0  # the pair tells them apart
    assert np.max(np.abs(y - y_ab))/peak < 2e-2


# ---------------------------------------------------------------- lead

@pytest.mark.parametrize("version", JULIA + [_backend("python"), _backend(3)])
def test_ed_requests_for_mps_only_quantities_raise_on_every_backend(version):
    """get_rdm and the bond entropy have no ED implementation. On julia_live
    get_rdm(mode="ED") returned the DMRG matrix and, with sc.mode="ED",
    both died with "'State' object has no attribute 'jlmps'"; they now
    raise NotImplementedError there as on v3 and "python". The site entropy
    has an ED route and keeps it, and the DMRG bond entropy still answers."""
    _require_julia(version)
    sc = spinchain.Spin_Chain(["S=1/2"]*N4, itensor_version=version)
    sc.set_hamiltonian(_ham4(sc))
    sc.maxm, sc.nsweeps = 16, 12
    _quiet(sc.gs_energy)
    s_bond = np.real(sc.get_bond_entropy(sc.get_gs(), 1, 2))
    with pytest.raises(NotImplementedError):
        sc.get_rdm(i=1, mode="ED")
    sc.mode = "ED"
    wf = sc.get_gs()
    with pytest.raises(NotImplementedError):
        sc.get_rdm(i=1)
    with pytest.raises(NotImplementedError):
        sc.get_bond_entropy(wf, 1, 2)
    s_site = np.real(sc.get_site_entropy(wf, 1))
    # the ED ground state's own bond-(1,2) entropy, from its Schmidt values:
    # sites 0..1 against 2..3 is a 4x4 reshape of the state vector in
    # either site order
    v = np.asarray(wf.v).ravel()
    p = np.linalg.svd(v.reshape(4, 4), compute_uv=False)**2
    p = p[p > 1e-14]/np.sum(p)
    assert s_bond == pytest.approx(-np.sum(p*np.log(p)), abs=1e-6)
    assert s_site > 0.0
