"""Regression tests for the `construction` cluster of the 2026-09-25b hole
hunt (docs/audit_2026_09_25b_hole_hunt.md, findings 2 to 6), and for three
leads the 2026-09-25 record left open.

  #2  On every route that resolves to ED (sc.mode="ED", a mode="ED" call,
      v3's own fallback below 3 sites) gs_energy() passed none of its
      keywords on: gs_energy(wf0=x, reconverge=False) left the ED ground
      state on the chain, a typo was swallowed, and get_gs(wf0=x) raised
      TypeError from EDchain.get_gs().
  #3  mode.resolve_mode returned the chain's own mode ahead of the call's,
      so a chain whose mode was "DMRG" answered an explicit mode="ED" call,
      the cross-check, by DMRG.
  #4  sites.check_settings admitted nine attributes that are not settings,
      the Hamiltonian accumulators among them (hubbard=2.0 put 2*Id into the
      next set_hoppings() Hamiltonian).
  #5  On a chain whose ground state is current, gs_energy()/get_gs() read no
      keyword but wf0=, so a typo was swallowed there and raised on a fresh
      chain; with it the leads maxde= (unrefined on a current chain) and
      get_gs(best=True, **kwargs) (forwarded to a best_gs that took none).
  #6  On a chain with no Hamiltonian every reader failed with an exception
      that did not name the cause.
  Lead: an unknown keyword to Thermal_Spin_Chain named Spin_Chain().

Anchors are ED (numpy's eigvalsh/eigh of the ED matrix), the state object
itself (<x|O|x>/<x|x>), and the same call on a chain that is not current.
"""

import io
import contextlib
import warnings

import numpy as np
import pytest

from dmrgpy import (bosonchain, cppext, fermionchain, groundstate,
                    manybodychain, mixedchain, parafermionchain, sites,
                    spinchain, thermal)


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("python"), _backend(3)]


def _quiet(f, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f(*a, **k)


def _heis(sc, n=None):
    n = sc.ns if n is None else n
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h


def _exact(sc):
    """eigenvalues and eigenvectors of the chain's ED Hamiltonian matrix"""
    H = sc.get_ED_obj().get_hamiltonian()
    H = H.toarray() if hasattr(H, "toarray") else np.asarray(H)
    return np.linalg.eigh(H)


def _ray(sc, a, b):
    """|<a|b>|^2/(<a|a><b|b>), invariant under the scale and phase of each"""
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))


def _expect(x, op):
    return np.real(x.dot(op*x)/x.dot(x))


# ---------------------------------------------------------------- #2

def _field_chain(n, version, assign_ed):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 10
    sc.set_hamiltonian(_heis(sc) + 0.4*sc.Sz[0] + 0.3*sc.Sx[n-1])
    if assign_ed: sc.mode = "ED"
    return sc


# the routes that resolve to ED: the chain's own mode, v3's automatic
# fallback below three sites (nobody names ED), and a mode="ED" call
ED_ROUTES = [
    pytest.param(4, "python", True, "DMRG", id="python-chain-mode-ED"),
    pytest.param(2, 3, False, "DMRG", id="v3-2-sites-fallback",
                 marks=BACKENDS[1].marks),
    pytest.param(4, "python", False, "ED", id="python-call-mode-ED"),
]


@pytest.mark.parametrize("n,version,assign_ed,mode", ED_ROUTES)
def test_ed_route_takes_wf0_as_it_is(n, version, assign_ed, mode):
    """after gs_energy(wf0=x, reconverge=False) every reader measures x, as
    after set_gs(x); it used to measure the ED ground state (|<gs|x>|^2 =
    0.0732, vev(Sz0) -0.2287 against <x|Sz0|x> = 0.0390 on 4 sites). The
    energy returned is the one gs_energy() returns after set_gs(x) on ED,
    the lowest eigenvalue (the 2026-09-24c record's open choice, not
    decided here)"""
    np.random.seed(7)
    sc = _field_chain(n, version, assign_ed)
    _quiet(sc.gs_energy, mode=mode) # solved first
    x = _quiet(sc.random_state, mode=mode)
    assert x.mode == "ED"
    e = _quiet(sc.gs_energy, mode=mode, wf0=x, reconverge=False)
    assert _ray(sc, _quiet(sc.get_gs, mode=mode), x) == pytest.approx(1.0, abs=1e-12)
    assert np.real(_quiet(sc.vev, sc.Sz[0], mode=mode)) == pytest.approx(
        _expect(x, sc.Sz[0]), abs=1e-12)
    # the energy: the same as the set_gs() route's, on an identical chain
    twin = _field_chain(n, version, assign_ed)
    _quiet(twin.gs_energy, mode=mode)
    _quiet(twin.set_gs, x)
    assert e == pytest.approx(_quiet(twin.gs_energy, mode=mode), abs=1e-12)
    assert e == pytest.approx(_exact(sc)[0][0], abs=1e-10)


@pytest.mark.parametrize("n,version,assign_ed,mode", ED_ROUTES)
def test_ed_route_get_gs_takes_wf0(n, version, assign_ed, mode):
    """get_gs(wf0=x, reconverge=False) raised TypeError from
    EDchain.get_gs(); it now returns x, as on DMRG"""
    np.random.seed(7)
    sc = _field_chain(n, version, assign_ed)
    x = _quiet(sc.random_state, mode=mode)
    w = _quiet(sc.get_gs, mode=mode, wf0=x, reconverge=False)
    assert _ray(sc, w, x) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("n,version,assign_ed,mode", ED_ROUTES)
def test_ed_route_refuses_what_it_cannot_read(n, version, assign_ed, mode):
    """a misspelled keyword was swallowed on ED whether or not the chain
    was current, and an MPS passed as wf0= was dropped silently"""
    np.random.seed(7)
    sc = _field_chain(n, version, assign_ed)
    for current in (False, True):
        if current: _quiet(sc.gs_energy, mode=mode)
        for f in (sc.gs_energy, sc.get_gs):
            with pytest.raises(TypeError) as excinfo:
                _quiet(f, mode=mode, reconverg=False)
            assert "reconverg" in str(excinfo.value)
    dm = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    mps = _quiet(dm.random_state) # an MPS, from a DMRG chain
    with pytest.raises(TypeError) as excinfo:
        _quiet(sc.gs_energy, mode=mode, wf0=mps, reconverge=False)
    assert "ED state" in str(excinfo.value)


def test_ed_route_warm_start_ends_on_the_ground_state():
    """wf0=x as a start is where any sweep from x ends: the exact ground
    state, which replaces a state set by hand before, as the sweep does on
    DMRG; maxde= is met by the exact state and changes nothing"""
    np.random.seed(7)
    sc = _field_chain(4, "python", True)
    ev, vecs = _exact(sc)
    gs_sz0 = np.real(np.conj(vecs[:, 0]) @ (
        sc.get_ED_obj().get_operator(sc.Sz[0]) @ vecs[:, 0]))
    y = sc.random_state()
    sc.set_gs(y)
    assert np.real(sc.vev(sc.Sz[0])) == pytest.approx(_expect(y, sc.Sz[0]), abs=1e-12)
    e = sc.gs_energy(wf0=sc.random_state())
    assert e == pytest.approx(ev[0], abs=1e-10)
    assert np.real(sc.vev(sc.Sz[0])) == pytest.approx(gs_sz0, abs=1e-10)
    assert sc.gs_energy(maxde=1e-6, maxdepth=2) == pytest.approx(ev[0], abs=1e-10)


def test_fermionic_chain_ed_branch_reads_the_keywords():
    """Fermionic_Chain.gs_energy had its own mode="ED" branch, a third
    entry point that dropped the keywords the same way"""
    np.random.seed(3)
    fc = fermionchain.Fermionic_Chain(4, itensor_version="python")
    h = 0
    for i in range(3):
        h = h - fc.Cdag[i]*fc.C[i+1] - fc.Cdag[i+1]*fc.C[i]
    fc.set_hamiltonian(h + 0.3*fc.N[0])
    e0 = fc.gs_energy(mode="ED")
    assert e0 == pytest.approx(_exact(fc)[0][0], abs=1e-10)
    x = fc.random_state(mode="ED")
    fc.gs_energy(mode="ED", wf0=x, reconverge=False)
    assert np.real(fc.vev(fc.N[0], mode="ED")) == pytest.approx(
        _expect(x, fc.N[0]), abs=1e-12)
    with pytest.raises(TypeError):
        fc.gs_energy(mode="ED", reconverg=False)


# ---------------------------------------------------------------- #3

def _staggered(sc):
    return _heis(sc) + sum(0.2*(-1)**i*sc.Sz[i] for i in range(sc.ns))


@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("route", ["constructor", "assigned"])
def test_a_dmrg_mode_on_the_chain_does_not_override_a_mode_ed_call(version,
                                                                  route):
    """gs_energy(mode="ED") on a chain whose mode was "DMRG" returned the
    truncated DMRG energy, -3.6734578613 against the exact -3.7040879103 at
    maxm=2 on 8 sites, so the cross-check compared DMRG with DMRG"""
    np.random.seed(2)
    if route == "constructor":
        sc = spinchain.Spin_Chain(["S=1/2"]*8, itensor_version=version,
                                  mode="DMRG", maxm=2, nsweeps=4)
    else:
        sc = spinchain.Spin_Chain(["S=1/2"]*8, itensor_version=version)
        sc.mode, sc.maxm, sc.nsweeps = "DMRG", 2, 4
    sc.set_hamiltonian(_staggered(sc))
    ev, vecs = _exact(sc)
    ref = spinchain.Spin_Chain(["S=1/2"]*8, itensor_version="python")
    ref.set_hamiltonian(_staggered(ref))
    zz = ref.vev(ref.Sz[0]*ref.Sz[1], mode="ED")
    assert _quiet(sc.gs_energy, mode="ED") == pytest.approx(ev[0], abs=1e-10)
    assert np.real(_quiet(sc.vev, sc.Sz[0]*sc.Sz[1], mode="ED")) == pytest.approx(
        np.real(zz), abs=1e-10)
    # the DMRG answer really is truncated here, so the two differ
    assert _quiet(sc.gs_energy, mode="DMRG") - ev[0] > 1e-3


def test_an_ed_mode_on_the_chain_still_overrides_the_call():
    """documentation.md 4.3: DMRG unless self.mode forces ED"""
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version="python",
                              maxm=2, nsweeps=2)
    sc.set_hamiltonian(_staggered(sc))
    sc.mode = "ED"
    assert sc.get_mode(mode="DMRG") == "ED"
    assert sc.gs_energy(mode="DMRG") == pytest.approx(_exact(sc)[0][0], abs=1e-10)


def test_the_ed_only_finite_temperature_correlation_matrix_is_ed():
    """get_correlation_matrix(T>0), documented as ED only, asks for
    get_excited_states(mode="ED"); on a chain whose mode was "DMRG" it got
    truncated DMRG states, <n_i> 0.17 off the Fermi function"""
    def occupations(mode):
        fc = fermionchain.Fermionic_Chain(4, itensor_version="python",
                                          maxm=2, nsweeps=4)
        fc.mode = mode
        h = 0
        for i in range(3):
            h = h - fc.Cdag[i]*fc.C[i+1] - fc.Cdag[i+1]*fc.C[i]
        fc.set_hamiltonian(h + 0.2*fc.N[0] - 0.2*fc.N[3])
        return np.real(np.diag(_quiet(fc.get_correlation_matrix, T=1.0)))
    t = np.zeros((4, 4))
    for i in range(3): t[i, i+1] = t[i+1, i] = -1.0
    t[0, 0], t[3, 3] = 0.2, -0.2
    eps, U = np.linalg.eigh(t)
    fermi = (np.abs(U)**2) @ (1.0/(np.exp(eps/1.0)+1.0))
    assert np.max(np.abs(occupations(None) - fermi)) < 1e-10
    assert np.max(np.abs(occupations("DMRG") - fermi)) < 1e-10


def test_the_thermal_chain_leaves_mbchain_on_its_default_mode():
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*2, T=0.5, itensor_version="python")
    assert tc.mode is None and tc.MBChain.mode is None


# ---------------------------------------------------------------- #4

CONSTRUCTORS = {
    "Spin_Chain": lambda **k: spinchain.Spin_Chain(["S=1/2"]*4, **k),
    "Fermionic_Chain": lambda **k: fermionchain.Fermionic_Chain(4, **k),
    "Spinful_Fermionic_Chain":
        lambda **k: fermionchain.Spinful_Fermionic_Chain(2, **k),
    "Majorana_Chain": lambda **k: fermionchain.Majorana_Chain(4, **k),
    "Bosonic_Chain": lambda **k: bosonchain.Bosonic_Chain(3, maxnb=[3]*3, **k),
    "SpinBoson_Chain":
        lambda **k: bosonchain.SpinBoson_Chain(["B", "S=1/2"], **k),
    "Parafermionic_Chain":
        lambda **k: parafermionchain.Parafermionic_Chain(4, **k),
    "Mixed_Spin_Fermion_Chain": lambda **k: mixedchain.Mixed_Spin_Fermion_Chain(
        ["S=1/2", "fermion"], **k),
    "Thermal_Spin_Chain": lambda **k: thermal.Thermal_Spin_Chain(
        ["S=1/2"]*2, T=0.5, **k).MBChain,
}
ACCUMULATORS = ("hopping", "hubbard", "pairing", "exchange")
DEAD = ("fields", "resorder", "resordered_indexes", "hubbard_matrix", "fit_td")


@pytest.mark.parametrize("name", sorted(CONSTRUCTORS))
@pytest.mark.parametrize("key", ACCUMULATORS + DEAD)
def test_the_nine_names_that_are_not_settings_are_refused(name, key):
    """all nine were admitted, since the rule was 'a public attribute the
    chain has, less STATE'; the accumulators are refused as state, naming
    their setter, and the five nothing read are gone from the chain"""
    with pytest.raises(TypeError) as excinfo:
        CONSTRUCTORS[name](itensor_version="python", **{key: 1.0})
    assert key in str(excinfo.value)
    if key in ACCUMULATORS:
        assert "set_hamiltonian()" in str(excinfo.value)


def test_a_constructor_accumulator_no_longer_shifts_the_energy():
    """Fermionic_Chain(4, hubbard=2.0) + set_hoppings gave -0.2360679775,
    the free-fermion -2.2360679775 plus the constant 2*Id"""
    with pytest.raises(TypeError):
        fermionchain.Fermionic_Chain(4, itensor_version="python", hubbard=2.0)
    fc = fermionchain.Fermionic_Chain(4, itensor_version="python")
    fc.set_hoppings(lambda i, j: -1.0 if abs(i-j) == 1 else 0.0)
    t = np.diag([-1.0]*3, 1) + np.diag([-1.0]*3, -1)
    lev = np.linalg.eigvalsh(t)
    assert fc.gs_energy(mode="ED") == pytest.approx(lev[lev < 0].sum(), abs=1e-10)


def test_every_setting_is_settable_and_every_chain_attribute_is_classified():
    """SETTINGS is an allowlist: every name in it is an attribute every
    chain has, none of it is state, and every public attribute
    Many_Body_Chain.__init__ sets is either a setting or state, so a new
    attribute cannot become a constructor keyword by accident"""
    assert not sites.SETTINGS & set(sites.STATE)
    for name, build in CONSTRUCTORS.items():
        c = build(itensor_version="python")
        for k in sites.SETTINGS:
            assert hasattr(c, k) and not callable(getattr(c, k)), (name, k)
        for k in DEAD:
            assert not hasattr(c, k), (name, k)
    sc = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version="python")
    same = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version="python",
                                **{k: getattr(sc, k) for k in sites.SETTINGS})
    for k in sites.SETTINGS:
        assert getattr(same, k) == getattr(sc, k), k
    bare = manybodychain.Many_Body_Chain([2, 2, 2], itensor_version="python")
    public = {k for k in vars(bare) if not k.startswith("_")}
    public -= {"itensor_version", "path", "inipath"} # named, and initialize()'s
    assert public <= sites.SETTINGS | set(sites.STATE), sorted(
        public - sites.SETTINGS - set(sites.STATE))


def test_an_unknown_thermal_keyword_names_thermal_spin_chain():
    """the error named Spin_Chain(), the chain the wrapper builds"""
    with pytest.raises(TypeError) as excinfo:
        thermal.Thermal_Spin_Chain(["S=1/2"]*2, itensor_version="python",
                                   maxm=7, zz_bogus=1)
    assert str(excinfo.value).startswith("Thermal_Spin_Chain()")
    assert "zz_bogus" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        thermal.Thermal_Spin_Chain(["S=1/2"]*2, itensor_version="python",
                                   hamiltonian=None)
    assert str(excinfo.value).startswith("Thermal_Spin_Chain()")


# ---------------------------------------------------------------- #5

def _chain(version, n=10, maxm=6, nsweeps=1):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm, sc.nsweeps, sc.noise = maxm, nsweeps, 0.0
    sc.bond_ramp = False
    sc.set_hamiltonian(_heis(sc))
    return sc


TYPOS = ({"wf": None}, {"wf_0": None}, {"reconverg": False}, {"maxdepht": 2})


@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("kw", TYPOS, ids=lambda k: list(k)[0])
def test_a_typo_raises_on_a_current_chain_as_on_a_fresh_one(version, kw):
    """on a solved chain the typo came back with the stored energy and
    state; on a fresh one gs_energy_single() raises TypeError"""
    np.random.seed(11)
    sc = _chain(version, n=4)
    kw = dict(kw)
    if "wf" in kw or "wf_0" in kw: kw[list(kw)[0]] = sc.random_state()
    for current in (False, True):
        if current: _quiet(sc.gs_energy)
        for f in (sc.gs_energy, sc.get_gs):
            with pytest.raises(TypeError) as excinfo:
                _quiet(f, **kw)
            assert list(kw)[0] in str(excinfo.value)
    with pytest.raises(TypeError):
        _quiet(sc.get_excited, n=1, **kw) # forwarded into both


@pytest.mark.parametrize("version", [pytest.param("julia_live", id="julia_live")])
def test_julia_live_raises_on_a_current_chain_as_on_a_fresh_one(version):
    """_gs_energy_julia() takes no maxde= and raises on reconverge= without
    a state; on a solved julia_live chain both came back with the stored
    energy, so the stored answer there now holds for no keyword but wf0="""
    from _helpers import julia_available
    ok, reason = julia_available()
    if not ok:
        pytest.skip("requires a working juliacall/Julia toolchain: %s" % reason)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.maxm, sc.nsweeps = 10, 4
    sc.set_hamiltonian(_heis(sc))
    for current in (False, True):
        if current:
            e0 = _quiet(sc.gs_energy)
            assert sc.gs_energy() == e0 and sc.gs_energy(wf0=None) == e0
        for kw in ({"reconverge": True}, {"maxde": 1e-3}, {"reconverg": False}):
            for f in (sc.gs_energy, sc.get_gs):
                with pytest.raises(TypeError):
                    _quiet(f, **kw)
    # a switch onto julia_live, which builds no session, leaves none behind
    sp = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    assert sp._session is not None
    sp.setup_julia()
    assert sp._session is None


@pytest.mark.parametrize("version", BACKENDS)
def test_a_current_chain_still_answers_what_it_already_answers(version):
    """the short circuit is kept for the calls it answers: no keyword,
    wf0=None, reconverge=False, maxde=None, maxdepth= alone"""
    sc = _chain(version, n=4, maxm=10, nsweeps=4)
    e0 = _quiet(sc.gs_energy)
    gs = sc.get_gs()
    for kw in ({}, {"wf0": None}, {"reconverge": False}, {"maxde": None},
               {"maxdepth": 3}):
        assert sc.gs_energy(**kw) == e0, kw
        assert sc.get_gs(**kw) is gs, kw


@pytest.mark.parametrize("version", BACKENDS)
def test_maxde_refines_a_current_chain_as_a_fresh_one(version):
    """gs_energy(maxde=...) on a current chain returned the stored energy
    unrefined, -4.1431954920 against -4.2580352072 on a fresh chain (10
    sites, maxm=3, exact -4.2580352073); get_gs widens with it"""
    np.random.seed(3)
    fresh = _chain(version, maxm=3, nsweeps=4)
    ef = _quiet(fresh.gs_energy, maxde=1e-4)
    for entry in ("gs_energy", "get_gs"):
        np.random.seed(3)
        sc = _chain(version, maxm=3, nsweeps=4)
        e1 = _quiet(sc.gs_energy)
        stored = sc.get_gs()
        _quiet(getattr(sc, entry), maxde=1e-4)
        assert sc.e0 < e1 - 1e-3, entry # refined
        assert sc.e0 == pytest.approx(ef, abs=1e-6), entry
        assert sc.get_gs() is not stored


def test_reconverge_true_sweeps_from_the_stored_state():
    """gs_energy_single's docstring: reconverge=True overrides the cached
    energy. On a current chain the call returned the stored energy (moved
    0.000e+00); it now sweeps from the stored state, which is the same
    number as a sweep from a copy of it passed as wf0="""
    def solve(how):
        np.random.seed(11)
        sc = _chain("python")
        _quiet(sc.gs_energy)
        e1 = sc.e0
        if how == "reconverge": return e1, _quiet(sc.gs_energy, reconverge=True)
        return e1, _quiet(sc.gs_energy, wf0=sc.get_gs())
    (e1, er), (_, ew) = solve("reconverge"), solve("wf0")
    assert er < e1 - 1e-6 # it swept
    assert er == pytest.approx(ew, abs=1e-8)


@pytest.mark.parametrize("version", BACKENDS)
def test_the_non_hermitian_route_reaches_the_solver_too(version):
    """gs_energy(H=H2) on a current non-Hermitian chain returned the stored
    energy of the chain's own H; every keyword the stored answer does not
    answer now reaches gs_energy_nhdmrg(), which solves for H2 here (and,
    since the session cluster's fix, refuses H= with TypeError), while the
    calls the stored answer does answer still return it"""
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.maxm, sc.nsweeps = 10, 4
    h = _heis(sc) + 0.3j*sc.Sz[0] + 0.2*sc.Sx[1]
    sc.set_hamiltonian(h)
    e0 = _quiet(sc.gs_energy)
    gs = sc.get_gs()
    for kw in ({}, {"wf0": None}, {"reconverge": False}, {"maxdepth": 3}):
        assert sc.gs_energy(**kw) == e0, kw
        assert sc.get_gs(**kw) is gs, kw
    try:
        e2 = _quiet(sc.gs_energy, H=2*h)
    except TypeError:
        e2 = None # gs_energy_nhdmrg() refuses H= (the session cluster's fix)
    if e2 is not None:
        assert abs(e2 - 2*e0) < 1e-6 # H2's eigenvalue, not H's stored one


def test_best_takes_the_keywords_of_each_solve():
    """get_gs(best=True, **kwargs) forwarded to a best_gs that took none,
    so any keyword next to best=True raised 'best_gs() got an unexpected
    keyword argument'; the keywords now go to each solve, which reads them
    or raises on them, and wf0= is refused by name"""
    np.random.seed(5)
    sc = _chain("python", n=4, maxm=10, nsweeps=4)
    ev = _exact(sc)[0]
    w = _quiet(sc.get_gs, best=True, n=2, maxde=1e-4)
    assert _expect(w, sc.hamiltonian) == pytest.approx(ev[0], abs=1e-8)
    _quiet(sc.gs_energy) # takes the state best_gs() set: current again
    for current in (True, False):
        if not current: sc = _chain("python", n=4)
        assert groundstate.gs_is_current(sc) == current
        with pytest.raises(TypeError) as excinfo:
            _quiet(sc.get_gs, best=True, n=2, reconverg=False)
        assert "reconverg" in str(excinfo.value)
        with pytest.raises(TypeError) as excinfo:
            _quiet(sc.get_gs, best=True, n=2, wf0=sc.random_state())
        assert "best=True" in str(excinfo.value)


@pytest.mark.parametrize("version", BACKENDS)
def test_a_state_that_cannot_be_set_leaves_the_chain_as_it_was(version):
    """set_initial_wf() wrote computed_gs=False and gs_from_file=True before
    mark_injected(), which refuses a state with no norm (the session
    cluster's fix), so the ValueError left a half-set chain; the flags now
    come after it. Without that refusal the call simply succeeds"""
    sc = _chain(version, n=4, maxm=10, nsweeps=4)
    e0 = _quiet(sc.gs_energy)
    zero = 0.0*sc.get_gs()
    flags = (sc.computed_gs, sc.gs_from_file, sc.skip_dmrg_gs)
    for guess in (False, True):
        try:
            sc.set_initial_wf(zero, reconverge=guess)
        except ValueError:
            assert (sc.computed_gs, sc.gs_from_file, sc.skip_dmrg_gs) == flags
            assert groundstate.gs_is_current(sc) and sc.gs_energy() == e0


@pytest.mark.parametrize("version", BACKENDS)
def test_a_backend_switch_leaves_no_session_of_the_old_backend(version):
    """_switch_backend() left the previous backend's session on the chain
    wherever initialize() builds none: julia_live (see the julia_live test
    above) and a chain whose mode is "ED" """
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    old = sc._session
    assert old is not None
    sc.mode = "ED"
    if version == "python":
        if not cppext.available(3): pytest.skip("needs the v3 extension")
        sc.setup_cpp(3)
    else:
        sc.setup_python()
    assert sc._session is None


# ---------------------------------------------------------------- #6

NO_H = {
    "Spin_Chain": lambda: spinchain.Spin_Chain(["S=1/2"]*3, itensor_version="python"),
    "Fermionic_Chain": lambda: fermionchain.Fermionic_Chain(3, itensor_version="python"),
    "Spinful_Fermionic_Chain":
        lambda: fermionchain.Spinful_Fermionic_Chain(2, itensor_version="python"),
    "Bosonic_Chain":
        lambda: bosonchain.Bosonic_Chain(3, maxnb=[3]*3, itensor_version="python"),
    "Parafermionic_Chain":
        lambda: parafermionchain.Parafermionic_Chain(3, itensor_version="python"),
}


def _op(c):
    for name in ("Sz", "N", "Ntot", "Sig"):
        if hasattr(c, name): return getattr(c, name)[0]
    raise AttributeError("no site operator to measure on %s" % type(c).__name__)


READERS = {
    "gs_energy": lambda c, m: c.gs_energy(mode=m),
    "get_gs": lambda c, m: c.get_gs(mode=m),
    "vev": lambda c, m: c.vev(_op(c), mode=m),
    "get_excited": lambda c, m: c.get_excited(n=2, mode=m),
    "get_excited_states": lambda c, m: c.get_excited_states(n=2, mode=m),
    "get_gap": lambda c, m: c.get_gap(mode=m),
    "get_dynamical_correlator": lambda c, m: c.get_dynamical_correlator(
        mode=m, name=(_op(c), _op(c)), es=np.linspace(-1, 1, 5), delta=0.2),
}


@pytest.mark.parametrize("mode", ["DMRG", "ED"])
@pytest.mark.parametrize("reader", sorted(READERS))
@pytest.mark.parametrize("name", sorted(NO_H))
def test_a_chain_with_no_hamiltonian_says_so(name, reader, mode):
    """'NoneType' object has no attribute 'get_dagger' on DMRG, and 'No
    active exception to reraise', 'op', 'shape' or 'T' on ED, depending on
    the class; now the ValueError get_hamiltonian() gives"""
    c = NO_H[name]()
    with pytest.raises(ValueError, match="set_hamiltonian"):
        _quiet(READERS[reader], c, mode)


def test_ed_operators_and_states_need_no_hamiltonian():
    """the refusal is at the readers, not in get_ED_obj(), so an ED state
    or operator can still be built before the Hamiltonian is set"""
    sc = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version="python")
    x = sc.random_state(mode="ED")
    assert x.mode == "ED"
    sc.toMPO(sc.Sz[0], mode="ED")


def test_an_unknown_ed_spin_operator_is_named():
    """pychain/build.py ended an unknown-name lookup in a bare raise, 'No
    active exception to reraise' (the 2026-09-24 record's ISy lead)"""
    sc = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version="python")
    sc.set_hamiltonian(_heis(sc))
    with pytest.raises(ValueError, match="ISy"):
        sc.vev(sc.get_operator("ISy", 0), mode="ED")
