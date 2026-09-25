"""Regression tests for the `cpp` cluster of the 2026-09-25b hole hunt
(docs/audit_2026_09_25b_hole_hunt.md, findings 11, 23, 24, 25, 27, 28, 29).

Every one of them is a threshold or a normalization that was right only at
one scale, or only for one reading of a quantity, and every test below pins
the property that was wrong rather than a golden number:

- 11: get_rdm() divided the state by <psi|psi> instead of its norm, so a set
  state c*s gave rho/c^2. Pinned: rho is invariant under the ray, and is the
  ED density matrix.
- 23: the KPM same-vector shortcut compared the two (unnormalized) images
  on an absolute 1e-10, so any pair of small images was declared equal and
  the autocorrelator recursion answered a cross correlator. Pinned: the
  correlator is exactly bilinear in the two operators.
- 24: v3/"python" VUMPS stopped the H_AC/H_C Lanczos on a residual against
  max(1,|lambda|), neither free of the units nor of an energy offset, so the
  gauge mismatch floored above tol. Pinned: convergence and energy are the
  unit-scale ones under s*H and under an onsite constant.
- 25: the NH-DMRG SRTieBreak window was 1e-6*(1+|remin|), absolute in small
  units and growing with an offset, so the sweep followed the previous bond
  onto an excited eigenpair. Pinned on the C++ sessions directly (the
  Python entry's own unit scale would mask the small-units trigger): E0(s*H)
  = s*E0(H), and E0(H+c) = E0(H)+c at c=1e6.
- 27: v3's restarted Arnoldi stopped on absolute breakdown and residual
  tests, so iDMRG in small units returned an unconverged density at
  converged=True. Pinned: the density is scale-covariant down to s=1e-13.
- 28: the session's solver scale was read from the raw term list while the
  MPO's was read after AutoMPO merged it. Pinned: a list whose duplicate
  strings cancel gives the energy of the operator it sums to.
- 29: with verbose on, ITensor's log shows the energies of the unit-scaled
  operator. Pinned: a scaled solve announces its factor, an unscaled one
  does not.
"""

import numpy as np
import pytest

from dmrgpy import cppext, infinitechain, spinchain


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


MPS_BACKENDS = [_backend("python"), _backend(3), _backend(2)]
CPP = [_backend(3), _backend(2)]
needs_v3 = pytest.mark.skipif(not cppext.available(3), reason="needs v3")
needs_v2 = pytest.mark.skipif(not cppext.available(2), reason="needs v2")


def _heisenberg(sc, n):
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h


# -- 11: get_rdm of a set state ------------------------------------------------

@pytest.mark.parametrize("version", MPS_BACKENDS)
def test_rdm_of_a_set_state_is_invariant_under_the_ray(version):
    """rho(i) after set_gs(c*s) is the same matrix at every c and the ED one;
    the field terms make it non-trivial (unequal diagonal, off-diagonal)."""
    n = 4
    def ham(sc):
        return _heisenberg(sc, n) + 0.3*sc.Sz[0] + 0.2*sc.Sx[2]
    ref = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    ref.set_hamiltonian(ham(ref))
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.set_hamiltonian(ham(sc))
    sc.maxm, sc.nsweeps = 16, 12
    sc.gs_energy()
    s = sc.get_gs().copy()
    s = s*(1.0/np.sqrt(np.real(s.dot(s))))
    sx = np.array([[0, .5], [.5, 0]]); sy = np.array([[0, -.5j], [.5j, 0]])
    sz = np.diag([.5, -.5])
    for i in (1, n-1):
        ev = [np.real(ref.vev(op[i], mode="ED")) for op in (ref.Sx, ref.Sy, ref.Sz)]
        exact = 0.5*np.eye(2) + 2*(ev[0]*sx + ev[1]*sy + ev[2]*sz)
        sc.set_gs(s)
        rho1 = np.asarray(sc.get_rdm(i=i))
        assert min(np.max(np.abs(rho1-exact)), np.max(np.abs(rho1-exact.T))) < 1e-8
        for c in (2.0, 0.5):
            sc.set_gs(s*c)
            rho = np.asarray(sc.get_rdm(i=i))
            assert np.real(np.trace(rho)) == pytest.approx(1.0, abs=1e-10)
            assert np.max(np.abs(rho - rho1)) < 1e-10


@pytest.mark.parametrize("version", MPS_BACKENDS)
def test_rdm_after_an_unswept_wf0_is_normalized(version):
    """The gs_energy(wf0=x, reconverge=False) route reaches get_rdm too."""
    n = 4
    np.random.seed(2)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.set_hamiltonian(_heisenberg(sc, n) + 0.3*sc.Sz[0])
    sc.maxm, sc.nsweeps = 16, 12
    sc.gs_energy()
    s = sc.get_gs().copy()
    s = s*(1.0/np.sqrt(np.real(s.dot(s))))
    rho1 = np.asarray(sc.get_rdm(i=1))
    sc.gs_energy(wf0=s*3.0, reconverge=False)
    rho = np.asarray(sc.get_rdm(i=1))
    assert np.real(np.trace(rho)) == pytest.approx(1.0, abs=1e-10)
    assert np.max(np.abs(rho - rho1)) < 1e-10


# -- 23: the KPM same-vector shortcut ----------------------------------------------

_ES = np.linspace(-0.5, 5.0, 200)


def _kpm_chain(version, n=6):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = 30; sc.nsweeps = 10
    sc.set_hamiltonian(_heisenberg(sc, n) + 0.3*sc.Sz[0])
    return sc


@pytest.mark.parametrize("version", MPS_BACKENDS)
@pytest.mark.parametrize("eps", [1.2e-10, 1e-11])
def test_kpm_cross_correlator_is_bilinear_in_small_operators(version, eps):
    """C[eps*Sz0, eps*Sz3] = eps^2 C[Sz0,Sz3]. Below eps = 1e-10/||(Sz3-Sz0)|gs>||
    = 1.22e-10 the absolute test returned C[Sz3,Sz3] instead, 2.07 of the peak
    off. The pre-fix answer is kept as the thing it must NOT equal."""
    sc = _kpm_chain(version)
    _, yAB = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[3]), delta=0.2, es=_ES)
    _, yBB = sc.get_dynamical_correlator(name=(sc.Sz[3], sc.Sz[3]), delta=0.2, es=_ES)
    yAB, yBB = np.asarray(yAB), np.asarray(yBB)
    pk = np.max(np.abs(yAB))
    _, y = sc.get_dynamical_correlator(name=(eps*sc.Sz[0], eps*sc.Sz[3]),
                                       delta=0.2, es=_ES)
    y = np.asarray(y)/eps**2
    assert np.max(np.abs(y - yAB))/pk < 1e-8
    assert np.max(np.abs(y - yBB))/pk > 1.0


@pytest.mark.parametrize("version", MPS_BACKENDS)
def test_kpm_pair_with_physically_small_images(version):
    """The route that failed on the parent too: a raising-operator pair on a
    nearly saturated state, (S-_0, S+_0 + S+_3) on -sum Sz + 0.005 sum SxSx,
    whose images have norms ~1e-3 before the eps=3e-8 multiplies them.
    Anchored on the same call with the shortcut switched off."""
    n = 6
    def ham(sc):
        h = 0
        for i in range(n):
            h = h - sc.Sz[i]
        for i in range(n-1):
            h = h + 0.005*sc.Sx[i]*sc.Sx[i+1]
        return h
    sp = lambda c, i: c.Sx[i] + 1j*c.Sy[i]
    sm = lambda c, i: c.Sx[i] - 1j*c.Sy[i]
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = 30; sc.nsweeps = 10
    sc.set_hamiltonian(ham(sc))
    eps = 3e-8
    pair = (eps*sm(sc, 0), eps*(sp(sc, 0) + sp(sc, 3)))
    _, y = sc.get_dynamical_correlator(name=pair, delta=0.2, es=_ES)
    sc.kpm_accelerate = False
    _, y_full = sc.get_dynamical_correlator(name=pair, delta=0.2, es=_ES)
    y, y_full = np.asarray(y), np.asarray(y_full)
    assert np.max(np.abs(y - y_full)) < 1e-8*np.max(np.abs(y_full))


@pytest.mark.parametrize("version", MPS_BACKENDS)
def test_get_distribution_is_bilinear_in_small_operators(version):
    """general_kpm goes through the same shortcut (1.935 of the peak off)."""
    sc = _kpm_chain(version)
    H = _heisenberg(sc, 6) + 0.3*sc.Sz[0]
    _, d0 = sc.get_distribution(X=H, A=sc.Sz[0], B=sc.Sz[3], scale=5)
    eps = 1e-11
    _, d1 = sc.get_distribution(X=H, A=eps*sc.Sz[0], B=eps*sc.Sz[3], scale=5)
    d0 = np.asarray(d0); d1 = np.asarray(d1)/eps**2
    assert np.max(np.abs(d1 - d0)) < 1e-8*np.max(np.abs(d0))


# -- 24: VUMPS in small units and under an offset ---------------------------------

_TFIM_E0 = -0.440127030549   # D=8 grouped VUMPS on s=1, run to tol=1e-10


def _tfim(ic, s=1.0, c=0.0):
    h = s*(ic.SzC[0]*ic.SzR[0] + 0.8*ic.SxC[0])
    if c:
        h = h + 4*c*ic.SzC[0]*ic.SzC[0]   # 4c*Sz*Sz = c*Id on a spin-1/2 site
    return h


def _v3_vumps(h_of, D, maxiter=400):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
    ic.set_hamiltonian(h_of(ic))
    c = ic._make_cpp_chain()
    e0, conv, nit, gm = c.vumps_ground_state(
        ic._h_intra.to_terms(jordan_wigner_transform=False),
        ic._h_inter.to_terms(jordan_wigner_transform=False),
        D, 1e-10, maxiter, 4, 30)
    return np.real(e0), conv, nit, gm


@needs_v3
@pytest.mark.parametrize("s,c", [(1e-2, 0.0), (1e-8, 0.0), (1.0, 100.0), (1.0, -100.0)],
                         ids=["s=1e-2", "s=1e-8", "c=+100", "c=-100"])
def test_v3_grouped_vumps_converges_in_small_units_and_under_an_offset(s, c):
    """D=8 puts H_AC (n=128) on the Lanczos path. The gauge mismatch floored
    at 2.9e-9 (s=1e-2), 1.5e-3 (s=1e-8) and 2-3e-9 (c=+-100) on the old
    test, all converged=False."""
    e, conv, nit, gm = _v3_vumps(lambda ic: _tfim(ic, s, c), 8)
    assert conv and gm < 1e-10
    assert (e - c)/s == pytest.approx(_TFIM_E0, abs=1e-9)


@needs_v3
@pytest.mark.parametrize("s", [1e-2, 1e-4])
def test_v3_sequential_vumps_converges_in_small_units(s):
    """A reach-2 coupling routes to the sequential solver, whose local solves
    are Lanczos at every D. Anchored on its own s=1 energy."""
    def h_of(ic, s):
        return s*(ic.SzC[0]*ic.SzR[0] + 0.8*ic.SxC[0]
                  + 0.2*ic.SzC[0]*ic.get_operator("Sz", 0, group=2))
    e1, conv1, _, _ = _v3_vumps(lambda ic: h_of(ic, 1.0), 4)
    e, conv, nit, gm = _v3_vumps(lambda ic: h_of(ic, s), 4)
    assert conv1 and conv and gm < 1e-10
    assert e/s == pytest.approx(e1, abs=1e-9)


@pytest.mark.parametrize("s,c", [(1e-2, 0.0), (1.0, 100.0)], ids=["s=1e-2", "c=+100"])
def test_python_vumps_converges_in_small_units_and_under_an_offset(s, c):
    """The "python" twin: 7.2e-10 at s=1e-2 and 3.0e-9 at c=100 before."""
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.set_hamiltonian(_tfim(ic, s, c))
    ic.maxm = 8; ic.maxiter = 400; ic.etol = 1e-10
    np.random.seed(3)
    e = np.real(ic.gs_energy())
    assert ic.converged and ic._vumps_result.gauge_mismatch < 1e-10
    assert (e - c)/s == pytest.approx(_TFIM_E0, abs=1e-9)


def test_the_hamiltonian_unit_ignores_constants_and_follows_the_units():
    """The scale both VUMPS backends now measure residuals in: the largest
    coefficient of a non-constant term, whatever the constant is spelled as."""
    from dmrgpy.pyitensor import idmrg, vumps
    for s, c, want in ((1.0, 0.0, 1.0), (1e-6, 0.0, 1e-6), (1.0, 100.0, 1.0),
                       (1e-3, -7.0, 1e-3)):
        ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
        ic.set_hamiltonian(_tfim(ic, s, c))
        sites_uc, _ = idmrg._build_automaton(ic._h_intra.op, ic._h_inter.op,
                                             ic.site_types, ic.n_uc)
        u = vumps._hamiltonian_unit(ic._h_intra.op, ic._h_inter.op, sites_uc, ic.n_uc)
        assert u == pytest.approx(want, rel=1e-12)


# -- 25: the NH-DMRG SRTieBreak window (C++ sessions, called directly) ---------------

def _nh_ham(sc, n):
    return _heisenberg(sc, n) + 0.3j*sc.Sz[0]


def _nh_direct(version, n, maxm, s=1.0, c=0.0):
    """Chain::nhdmrg on s*H + c*Id with no Python-side unit scale, so the
    session's own Arnoldi and tie window are what is measured."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = maxm; sc.nsweeps = 10
    H = s*_nh_ham(sc, n)
    if c:
        H = H + c*sc.get_operator("Id", 0)
    sc.set_hamiltonian(H)
    sc._session.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
    sc._session.set_mpomaxm(max(sc.maxm, sc.mpomaxm))
    e, _, _ = sc._session.nhdmrg(H.to_terms(), H.get_dagger().to_terms(), 20, 2)
    return (complex(e) - c)/s


def _nh_ed(n):
    ref = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    Hm = ref.get_ED_obj().get_operator(_nh_ham(ref, n))
    Hm = Hm.toarray() if hasattr(Hm, "toarray") else np.asarray(Hm)
    ev = np.linalg.eigvals(Hm)
    return ev[np.argsort(ev.real)]


@needs_v2
@pytest.mark.parametrize("s,c", [(2e-6, 0.0), (1e-7, 0.0), (1.0, 1e6)],
                         ids=["s=2e-6", "s=1e-7", "c=1e6"])
def test_v2_nhdmrg_stays_on_the_lowest_level(s, c):
    """6 sites, gap 0.449: the window 1e-6*(1+|remin|) exceeded it below
    s = 2.2e-6 and above c ~ 5e5, and v2 returned levels #1 to #19 (0.449 to
    1.85 off) in every run measured."""
    ev = _nh_ed(6)
    e = _nh_direct(2, 6, 30, s, c)
    assert abs(e - ev[0]) < 1e-6


@needs_v3
@pytest.mark.parametrize("s", [2e-6, 1e-6])
def test_v3_nhdmrg_below_full_bond_dimension_is_scale_covariant(s):
    """v3 escapes on 6 sites at maxm=30; the record's v3 probe is 10 sites at
    maxm=8, where s=2e-6 landed on levels #1..#3. Anchored on the s=1 run of
    the same truncated calculation."""
    e1 = _nh_direct(3, 10, 8)
    e = _nh_direct(3, 10, 8, s)
    assert abs(e - e1) < 1e-6


# -- 27: v3 iDMRG in small units ------------------------------------------------------

@needs_v3
@pytest.mark.parametrize("s", [1e-11, 1e-13])
def test_v3_idmrg_density_is_scale_covariant(s):
    """1e-5..1e-3 relative off at s=1e-11..1e-12 and O(1) off (wrong sign at
    1e-14) from 1e-13 down, all at converged=True, before. etol is an energy
    and is scaled with the units."""
    def run(s):
        ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)
        ic.gs_method = "idmrg"
        ic.set_hamiltonian(s*(ic.SzC[0]*ic.SzR[0] + 0.8*ic.SxC[0]))
        ic.maxm = 16; ic.maxiter = 120; ic.etol = 1e-12*s
        return np.real(ic.gs_energy())/s, ic.converged
    e1, conv1 = run(1.0)
    e, conv = run(s)
    assert conv1 and conv
    assert e == pytest.approx(e1, rel=1e-9)


# -- 28: the solver scale is the MPO's ------------------------------------------------

@pytest.mark.parametrize("version", CPP)
@pytest.mark.parametrize("tag", ["Sz0", "Id"])
def test_cancelling_terms_do_not_change_the_solver_scale(version, tag):
    """s*H + X - X is s*H. With the scale read from the raw list (largest
    coefficient 1) and the MPO's from the merged AutoMPO (s), the solve ran
    unscaled on a unit-scaled MPO: 0.12 to 0.57 off at s=1e-10."""
    n = 6
    ref = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    ref.set_hamiltonian(_heisenberg(ref, n))
    e_ed = np.real(ref.gs_energy(mode="ED"))
    s = 1e-10
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = 30; sc.nsweeps = 10
    X = sc.Sz[0] if tag == "Sz0" else sc.get_operator("Id", 2)
    sc.set_hamiltonian(s*_heisenberg(sc, n) + X - X)
    assert np.real(sc.gs_energy())/s == pytest.approx(e_ed, abs=1e-9)


@pytest.mark.parametrize("version", CPP)
def test_cancelling_terms_excited_states(version):
    n = 6
    ref = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    Hm = ref.get_ED_obj().get_operator(_heisenberg(ref, n))
    Hm = Hm.toarray() if hasattr(Hm, "toarray") else np.asarray(Hm)
    ev = np.sort(np.linalg.eigvalsh(Hm))[:3]
    s = 1e-10
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = 30; sc.nsweeps = 10
    sc.set_hamiltonian(s*_heisenberg(sc, n) + sc.Sz[0] - sc.Sz[0])
    ex = np.sort(np.real(sc.get_excited(n=3)))/s
    assert np.max(np.abs(ex - ev)) < 1e-6


# -- 29: the verbose log of a scaled solve says so -------------------------------------

@pytest.mark.parametrize("version", CPP)
def test_verbose_log_announces_the_unit_scale(version, capfd):
    n = 6
    for s, announced in ((1.0, False), (0.5, True)):
        sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
        sc.maxm = 30; sc.nsweeps = 4; sc.verbose = 1
        sc.set_hamiltonian(s*_heisenberg(sc, n))
        capfd.readouterr()
        e = np.real(sc.gs_energy())
        out = capfd.readouterr().out
        assert ("2^1 = 2 times" in out) == announced
        assert "Energy after sweep" in out
        assert e/s == pytest.approx(-2.4935771339, abs=1e-8)
