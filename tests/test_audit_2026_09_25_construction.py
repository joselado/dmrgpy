"""Regression tests for the `construction` cluster of the 2026-09-25 fixes
of the open items the 2026-09-24 records left behind.

  init-kwargs. Many_Body_Chain.__init__(**kwargs) handed its keywords to
      initialize(), which ignores them, so Spin_Chain(sites, maxm=4) ran at
      maxm=30 and a misspelled keyword was accepted without a word, on
      every model chain, since they all forward to it. A keyword naming a
      chain setting now takes effect as if assigned right after
      construction, mode= included, and anything else raises TypeError
      naming every offending key.
  get-gs-wf0. get_gs(wf0=x) on a chain whose ground state was current
      returned the stored state without reading x, where gs_energy(wf0=x)
      reads it; both short circuits now have the same condition.
  rootn-ij. On the lower-level get_dynamical_correlator_MB route ROOTN
      dropped i=/j= of a string name (C[Sz_0,Sz_0] for any sites) and
      accepted a misspelled keyword.

Anchors are ED and the explicit operator pair on chains of 4 to 10 sites.
"""

import io
import contextlib
import warnings

import numpy as np
import pytest

from dmrgpy import (bosonchain, cppext, fermionchain, mixedchain,
                    parafermionchain, spinchain, thermal)


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("python"), _backend(3)]


def _quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f()


def _heis(sc, B=0.0):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns):
        h = h + B*sc.Sz[i]
    return h


# ---------------------------------------------------------- init-kwargs

# every constructor that forwards its keywords to Many_Body_Chain.__init__,
# returning the Many_Body_Chain that holds the settings
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
SETTINGS = dict(maxm=7, nsweeps=3, kpmmaxm=11, noise=0.0, cutoff=1e-9,
                kpm_scale=0.6, tevol_method="TEBD", verbose=True)


@pytest.mark.parametrize("name", sorted(CONSTRUCTORS))
def test_a_constructor_setting_takes_effect(name):
    """every one of these read maxm=30, nsweeps=15, kpmmaxm=50"""
    c = CONSTRUCTORS[name](itensor_version="python", **SETTINGS)
    for k, v in SETTINGS.items():
        assert getattr(c, k) == v, k
    assert c.itensor_version == "python"


@pytest.mark.parametrize("name", sorted(CONSTRUCTORS))
def test_a_keyword_that_is_not_a_setting_raises_naming_every_one(name):
    """bogus_key=1 and kpm_nscale= (for kpm_n_scale) were accepted"""
    with pytest.raises(TypeError) as excinfo:
        CONSTRUCTORS[name](itensor_version="python", maxm=50, zz_bogus=1,
                           kpm_nscale=3)
    assert "kpm_nscale, zz_bogus" in str(excinfo.value)
    assert "maxm" not in str(excinfo.value).split(".")[0]


@pytest.mark.parametrize("key", ["_session", "_maxm", "gs_energy", "set_gs",
                                 "path"])
def test_private_names_methods_and_names_the_chain_lacks_are_not_settings(key):
    with pytest.raises(TypeError) as excinfo:
        spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python",
                             **{key: None})
    assert key in str(excinfo.value)


@pytest.mark.parametrize("key", ["hamiltonian", "conserved_sector", "wf0",
                                 "e0", "computed_gs", "ns",
                                 "use_ampo_hamiltonian"])
def test_the_chain_state_is_not_a_setting(key):
    """these exist on the chain, so the attribute test alone would pass
    them; each has its own entry point, which the message names"""
    with pytest.raises(TypeError) as excinfo:
        spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python",
                             **{key: None})
    assert key in str(excinfo.value)
    assert "not a setting" in str(excinfo.value)


def test_a_constructor_setting_goes_through_its_setter():
    """maxm=0 used to be dropped with every other keyword; the maxm
    property refuses it when assigned, and so it does at construction"""
    with pytest.raises(ValueError) as excinfo:
        spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python", maxm=0)
    assert "maxm" in str(excinfo.value)


def test_a_constructor_setting_is_the_same_as_assigning_it():
    """a 10-site chain asked for maxm=4 at construction ran at maxm=30 and
    returned the converged energy; it now returns the maxm=4 one, which is
    the number the same request made by assignment returns"""
    def solve(at_construction):
        np.random.seed(3)
        kw = dict(maxm=4, nsweeps=6) if at_construction else {}
        sc = spinchain.Spin_Chain(["S=1/2"]*10, itensor_version="python", **kw)
        if not at_construction: sc.maxm, sc.nsweeps = 4, 6
        sc.set_hamiltonian(_heis(sc))
        return sc.gs_energy(), sc.gs_energy(mode="ED")
    (e_ctor, exact), (e_set, _) = solve(True), solve(False)
    assert e_ctor == pytest.approx(e_set, abs=1e-12)
    assert e_ctor - exact > 1e-4 # the truncation is really in force


@pytest.mark.parametrize("name", sorted(CONSTRUCTORS))
def test_mode_at_construction_keeps_the_session(name):
    """mode= is applied after initialize(), as every other setting is, so
    mode="ED" at construction is sc.mode = "ED" afterwards: the session is
    built all the same. The first fix applied mode first, and the chain
    then had no session to go back to"""
    c = CONSTRUCTORS[name](itensor_version="python", mode="ED")
    assert c.mode == "ED" # was None
    assert c._session is not None


def test_mode_at_construction_can_be_left_again():
    """a chain built at mode="ED" answers by ED, and after sc.mode = None by
    DMRG, which is what assigning mode="ED" afterwards gives; the first fix
    raised AttributeError on the second read"""
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version="python",
                              mode="ED", maxm=20, nsweeps=10)
    assert sc.mode == "ED"
    sc.set_hamiltonian(_heis(sc, B=0.2))
    e_ed = sc.gs_energy()
    assert e_ed == pytest.approx(sc.gs_energy(mode="ED"), abs=1e-12)
    sc.mode = None
    assert _quiet(sc.gs_energy) == pytest.approx(e_ed, abs=1e-8)


@pytest.mark.parametrize("T", [0.0, 0.5])
def test_thermal_chain_with_mode_at_construction_still_solves(T):
    """Thermal_Spin_Chain writes its own mode onto MBChain in get_gs(); with
    mode="ED" applied before the session was built that raised
    AttributeError where the unfixed tree solved"""
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*3, T=T,
                                    itensor_version="python", mode="ED")
    h = 0
    for i in range(2):
        h = h + tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
    tc.set_hamiltonian(h)
    wf = _quiet(tc.get_gs)
    assert np.isfinite(np.real(tc.MBChain.vev(tc.Sz[0]*tc.Sz[1], wf=wf)))


def test_thermal_chain_mode_at_construction_is_the_mode_it_solves_with():
    """Thermal_Spin_Chain's get_gs() overwrote MBChain.mode with its own
    hardcoded "DMRG", so Thermal_Spin_Chain(..., mode="ED") ran DMRG; the
    wrapper now keeps the mode it was given and hands it to MBChain"""
    tc = thermal.Thermal_Spin_Chain(["S=1/2"]*3, T=0.0,
                                    itensor_version="python", mode="ED")
    h = 0
    for i in range(2):
        h = h + tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
    tc.set_hamiltonian(h)
    _quiet(tc.get_gs)
    assert tc.mode == "ED" and tc.MBChain.mode == "ED"
    assert thermal.Thermal_Spin_Chain(["S=1/2"]*2).mode == "DMRG"
    with pytest.raises(ValueError):
        thermal.Thermal_Spin_Chain(["S=1/2"]*2, mode="ed")


@pytest.mark.parametrize("name,key", [("Fermionic_Chain", "N"),
                                      ("Fermionic_Chain", "Cdag"),
                                      ("Parafermionic_Chain", "Sig")])
def test_what_the_model_class_builds_is_not_a_setting(name, key):
    """the operator lists a model class builds before calling
    Many_Body_Chain.__init__ exist when the keywords are checked, so the
    attribute test alone admitted them: Fermionic_Chain(4, N=5).N came
    back as 5 and Parafermionic_Chain(3, Sig=1) failed inside its own
    constructor"""
    with pytest.raises(TypeError) as excinfo:
        CONSTRUCTORS[name](itensor_version="python", **{key: 1})
    assert key in str(excinfo.value)
    assert "not a setting" in str(excinfo.value)


def test_the_v2_boson_refusal_does_not_offer_mode_ed():
    """mode="ED" at construction takes effect after the session, so it is
    no way out of a boson dimension itensor_version=2 lacks, and the
    message no longer says it is (the check runs before any extension is
    loaded, so this needs no compiled v2)"""
    with pytest.raises(ValueError) as excinfo:
        bosonchain.Bosonic_Chain(3, maxnb=[6]*3, itensor_version=2,
                                 mode="ED")
    assert "itensor_version=3" in str(excinfo.value)
    assert "or mode=" not in str(excinfo.value)


def test_an_unknown_mode_is_refused_at_construction():
    with pytest.raises(ValueError) as excinfo:
        spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python", mode="ed")
    assert "'ed'" in str(excinfo.value)


# ----------------------------------------------------------- get-gs-wf0

def _chain(version, n=6):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.set_hamiltonian(_heis(sc) + 0.3*sc.Sz[0])
    sc.maxm, sc.nsweeps = 20, 10
    return sc


def _energy(sc, w):
    return np.real(sc.vev(sc.hamiltonian, wf=w))/np.real(sc.overlap(w, w))


def _fidelity(sc, a, b):
    return abs(sc.overlap(a, b))**2/abs(sc.overlap(a, a)*sc.overlap(b, b))


@pytest.mark.parametrize("version", BACKENDS)
def test_get_gs_takes_wf0_as_it_is_on_a_current_chain(version):
    """get_gs(wf0=x, reconverge=False) on a solved chain returned the solved
    state, |<w|x>|^2 = 0.011 and e0 = E_0, where gs_energy() with the same
    keywords returns <x|H|x>"""
    np.random.seed(1)
    sc = _chain(version)
    sc.gs_energy()
    x = sc.random_state()
    w = sc.get_gs(wf0=x, reconverge=False)
    assert _fidelity(sc, w, x) == pytest.approx(1.0, abs=1e-10)
    assert sc.e0 == pytest.approx(_energy(sc, x), abs=1e-10)


@pytest.mark.parametrize("version", BACKENDS)
def test_get_gs_warm_starts_from_wf0_on_a_current_chain(version):
    """with reconverge unspecified x is a warm start: the call sweeps from
    it rather than returning the stored object, and lands on the ground
    state"""
    np.random.seed(1)
    sc = _chain(version)
    e0 = sc.gs_energy()
    gs = sc.get_gs()
    w = sc.get_gs(wf0=sc.random_state())
    assert w is not gs
    assert sc.e0 == pytest.approx(e0, abs=1e-8)
    assert _fidelity(sc, w, gs) == pytest.approx(1.0, abs=1e-8)


@pytest.mark.parametrize("version", BACKENDS)
def test_get_gs_without_wf0_still_returns_the_stored_state(version):
    """the short circuit itself is kept: a repeated read costs nothing"""
    sc = _chain(version)
    sc.gs_energy()
    gs = sc.get_gs()
    assert sc.get_gs() is gs
    assert sc.get_gs(wf0=None) is gs


# ------------------------------------------------------------- rootn-ij

@pytest.fixture(scope="module")
def heisenberg4():
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    sc.set_hamiltonian(_heis(sc))
    sc.maxm, sc.nsweeps = 20, 8
    return sc


ROOTN = dict(es=np.linspace(0.0, 3.0, 5), delta=0.4, N=4, nkry=12)


@pytest.mark.parametrize("i,j", [(1, 1), (0, 2)])
def test_rootn_lower_level_route_honours_the_sites(heisenberg4, i, j):
    """get_dynamical_correlator_MB(name="ZZ", i=, j=, submode="ROOTN") is
    the explicit [Sz_i, Sz_j] pair; it used to be [Sz_0, Sz_0] bit for
    bit, 4.1e-02 away on a 0.10 peak at (1,1)"""
    sc = heisenberg4
    y = np.asarray(sc.get_dynamical_correlator_MB(
        submode="ROOTN", name="ZZ", i=i, j=j, **ROOTN)[1])
    y_ij = np.asarray(sc.get_dynamical_correlator(
        submode="ROOTN", name=[sc.Sz[i], sc.Sz[j]], **ROOTN)[1])
    y_00 = np.asarray(sc.get_dynamical_correlator(
        submode="ROOTN", name=[sc.Sz[0], sc.Sz[0]], **ROOTN)[1])
    assert np.max(np.abs(y - y_ij)) < 1e-12
    assert np.max(np.abs(y - y_00)) > 1e-3


def test_rootn_rejects_an_unknown_keyword(heisenberg4):
    """nkyr= for nkry= returned the default-parameter spectrum"""
    sc = heisenberg4
    with pytest.raises(TypeError) as excinfo:
        sc.get_dynamical_correlator_MB(submode="ROOTN", name="ZZ", nkyr=3,
                                       **ROOTN)
    assert "nkyr" in str(excinfo.value)
