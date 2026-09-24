"""Regressions for the ground-state cluster of the 2026-09-24b hole hunt
(docs/audit_2026_09_24b_hole_hunt.md, findings 11 to 15).

One root: the pybind port kept the trigger of the file-based backend's
hand-off of a stored state to the C++ program (the `set_initial_wf(self.wf0)`
every correlator opened with) and dropped the hand-off itself, so a state
the caller set never reached the DMRG session.

- Finding 11. After `set_gs()` every ground-state-reading correlator
  submode measured the session's own solved state and then wrote it back
  over the one that was set. Pinned by the two pure members of the 3-site
  Heisenberg doublet, which mode="ED" separates by 0.4378 on a 0.445 peak:
  each submode must give each member's own density, and leave the member
  set, on the Python side and on the session.
- Finding 12. `set_initial_wf`/`set_initial_wf_guess` never reached the
  session either, and on "python" the explicit `gs_energy(wf0=x)` moved
  the caller's own `x`. Pinned by a spy on `session.set_wavefunction`, the
  overlap of the result with the target, <H> of the caller's state before
  and after, and the transverse-field Ising example's seeded branch.
- Finding 13. `get_kondo_spectrum(n_gs>1)` ran a hidden sweep from each
  member, and on v3 at a split below `delta` the upper member relaxed
  into the lower one in most runs. Pinned over several runs, since the
  defect showed in 16 of 24.
- Finding 14. After `n_gs>1`, `gs_energy()` was the last member's energy.
  Pinned at a split below `delta`, since eps=0 cannot see it.
- Finding 15. `n_gs>1` under `submode="EX"` returned the n_gs=1 value;
  EX now measures from the chain's own state reprojected onto its cached
  basis, and SECTOR, which cannot, raises.

Two more things the fix had to get right are pinned as well: the session's
band-edge cache, which the first KPM call after a push would otherwise fill
by sweeping the pushed state (seen only with a state that is not an
eigenstate), and the two cases that used to lean on the removed line, a
clone of a solved chain and a repeated correlator on a solved chain, which
must not re-sweep.
"""
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, groundstate, mps, spinchain
from dmrgpy.edtk.edchain import State


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("python"), _backend(3)]
ES = np.linspace(-1.0, 3.0, 21)
DELTA = 0.3


def heisenberg(version, n=3):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 10, 10
    return sc


def energy(sc, wf):
    return float(np.real(wf.aMb(sc.hamiltonian, wf)/wf.dot(wf)))


def fidelity(a, b):
    return abs(a.dot(b))**2/abs(a.dot(a)*b.dot(b))


def session_state(sc):
    return mps.MPS(MBO=sc, cpp_handle=sc._session.gs_wavefunction())


class _Doublet:
    """The 3-site Heisenberg ground doublet: exact members (dense eigh, no
    dmrgpy algebra) and each member's exact Lehmann density"""
    def __init__(self):
        sc = heisenberg("python")
        self.ed = sc.get_ED_obj()
        H = np.array(self.ed.get_hamiltonian().todense())
        Szt = np.array(self.ed.MO2matrix(sc.Sz[0] + sc.Sz[1] + sc.Sz[2]).todense())
        self.e, self.U = np.linalg.eigh(H)
        P = self.U[:, :2]
        w, V = np.linalg.eigh(P.conj().T @ Szt @ P)
        self.up, self.dn = P @ V[:, 1], P @ V[:, 0]
        self.sc = sc
        sc.get_gs(mode="ED")

    def lorentzian(self, member, A, B):
        Am = np.array(self.ed.MO2matrix(A).todense())
        Bm = np.array(self.ed.MO2matrix(B).todense())
        M = (member.conj() @ Am @ self.U)*(self.U.conj().T @ Bm @ member)
        x = ES[:, None] - (self.e - self.e[0])[None, :]
        return ((DELTA/np.pi)/(x**2 + DELTA**2)) @ M

    def ed_kpm(self, member, A, B):
        self.sc.set_gs(State(member, self.ed))
        return np.asarray(self.sc.get_dynamical_correlator(
                mode="ED", submode="KPM", name=[A, B], es=ES, delta=DELTA)[1])


@pytest.fixture(scope="module")
def doublet():
    return _Doublet()


def members(sc):
    """The two pure members of the doublet, built from the solved state as
    (Sz_tot -+ 1/2)|s>, which is exactly its 2Sz=+-1 component"""
    Szt = sc.Sz[0] + sc.Sz[1] + sc.Sz[2]
    s = sc.get_gs().copy()
    up = (Szt + 0.5)*s
    dn = (0.5 - 1*Szt)*s
    return Szt, up*(1/np.sqrt(up.dot(up).real)), dn*(1/np.sqrt(dn.dot(dn).real))


# ------------------------------------------------------------ finding 11

SUBMODES = [("KPM", {}, 5e-3), ("CVM", {}, 1e-6), ("ROOTN", {}, 1e-6),
            ("TD", dict(dt=0.05), 1e-4), ("TDZ", dict(dt=0.05), 5e-3),
            ("EX", dict(nex=8), 1e-6), ("CVM_explicit", {}, 1e-6)]


def operator_pair(sc, submode):
    """(S+_0, S-_2), whose densities in the two members differ at the level
    of the peak; CVM_explicit takes Hermitian pairs only, (S-_0^dag, S-_0)"""
    if submode == "CVM_explicit":
        Sm = sc.Sx[0] - 1j*sc.Sy[0]
        return [Sm.get_dagger(), Sm]
    return [sc.Sx[0] + 1j*sc.Sy[0], sc.Sx[2] - 1j*sc.Sy[2]]


@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("submode,extra,tol", SUBMODES,
                         ids=[s[0] for s in SUBMODES])
def test_every_submode_measures_the_member_set_gs_set(version, submode,
                                                      extra, tol, doublet):
    """KPM is held to mode="ED" KPM on the same member, the rest to the
    member's exact Lorentzian density; before the fix every one of them
    returned the solved state's curve for both members"""
    np.random.seed(1)
    sc = heisenberg(version)
    Szt, up, dn = members(sc)
    ref = doublet.ed_kpm if submode == "KPM" else doublet.lorentzian
    kw = dict(mode="DMRG", submode=submode, name=operator_pair(sc, submode),
              es=ES, delta=DELTA, **extra)
    got = {}
    for label, wf, exact in (("up", up, doublet.up), ("dn", dn, doublet.dn)):
        sc.set_gs(wf)
        y = np.asarray(sc.get_dynamical_correlator(**kw)[1])
        y_ref = ref(exact, *operator_pair(doublet.sc, submode))
        assert np.max(np.abs(y - y_ref)) < tol
        # the member it measured is the one set, and it is still set, on the
        # Python side and on the session, with its own energy
        assert fidelity(sc.get_gs(), wf) == pytest.approx(1.0, abs=1e-10)
        assert fidelity(session_state(sc), wf) == pytest.approx(1.0, abs=1e-10)
        assert sc.gs_energy() == pytest.approx(energy(sc, wf), abs=1e-12)
        got[label] = y
    # and the two members are told apart as mode="ED" tells them apart
    assert np.max(np.abs(got["up"] - got["dn"])) > 0.4


@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("setter", ["set_gs", "set_initial_wf"])
def test_injected_state_is_not_swept_by_the_first_kpm_call(version, setter):
    """KPM rescales with the session's band edges, whose lower one is filled
    lazily by a session gs_energy(skip_dmrg=True), which sweeps whenever the
    session holds no energy for its state -- and handing it a state drops
    that energy. With no correlator run before, the first KPM call after the
    hand-off used to sweep the pushed state in place (|<x|session>|^2 came
    out 0.989 on python and 0.987 on v3). The state is not an eigenstate, so
    any sweep moves it."""
    np.random.seed(3)
    sc = heisenberg(version)
    s = sc.get_gs()
    x = s + 0.3*(sc.Sx[0]*s)
    x = x*(1/np.sqrt(x.dot(x).real))
    e_x = energy(sc, x)
    assert e_x > sc.gs_energy() + 1e-3 # really not the ground state
    getattr(sc, setter)(x)
    sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=ES, delta=DELTA)
    assert fidelity(session_state(sc), x) == pytest.approx(1.0, abs=1e-12)
    assert fidelity(sc.get_gs(), x) == pytest.approx(1.0, abs=1e-12)
    assert sc.gs_energy() == pytest.approx(e_x, abs=1e-12)


# ------------------------------------------------------------ finding 12

class _Spy:
    """Stands in for the session and logs the calls that matter"""
    def __init__(self, sc):
        self.real, self.sc, self.log = sc._session, sc, []

    def __getattr__(self, name):
        attr = getattr(self.real, name)
        if name in ("set_wavefunction", "gs_energy", "set_hamiltonian",
                    "excited_states"):
            def logged(*a, **k):
                if name == "set_wavefunction":
                    self.log.append((name, mps.MPS(MBO=self.sc, cpp_handle=a[0])))
                else:
                    self.log.append((name, None))
                return attr(*a, **k)
            return logged
        return attr


def spy_on(sc):
    spy = _Spy(sc)
    sc._session = spy
    if sc._session_ham_cache is not None: # keep the send-cache hit
        sc._session_ham_cache = (spy, sc._session_ham_cache[1])
    return spy


@pytest.mark.parametrize("version", BACKENDS)
def test_warm_start_setters_reach_the_session(version):
    np.random.seed(1)
    sc = heisenberg(version)
    sc.noise, sc.cutoff = 1e-7, 1e-12
    Szt, up, dn = members(sc)
    sz = np.real(sc.vev(Szt))
    target = dn if sz > 0 else up # the member the solved state is far from
    spy = spy_on(sc)
    # set_initial_wf_guess: the next solve sweeps from the guess, which is an
    # exact eigenstate, so the result stays on it
    sc.set_initial_wf_guess(target)
    sc.gs_energy()
    pushed = [w for n, w in spy.log if n == "set_wavefunction"]
    assert len(pushed) == 1 and fidelity(pushed[0], target) == pytest.approx(1.0, abs=1e-12)
    assert "gs_energy" in [n for n, _ in spy.log] # it did sweep
    assert fidelity(sc.get_gs(), target) == pytest.approx(1.0, abs=1e-8)
    # set_initial_wf: taken as it is, no sweep at all
    spy.log.clear()
    sc.set_initial_wf(up if target is dn else dn)
    other = up if target is dn else dn
    sc.get_gs()
    names = [n for n, _ in spy.log]
    assert "gs_energy" not in names
    assert fidelity([w for n, w in spy.log if n == "set_wavefunction"][0], other) \
        == pytest.approx(1.0, abs=1e-12)
    assert fidelity(sc.get_gs(), other) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("version", BACKENDS)
def test_the_callers_state_is_not_mutated(version):
    """On "python" set_wavefunction stores the handle it is given, and the
    next sweep rewrites that MPS in place: gs_energy(wf0=x) used to move the
    caller's own x onto the ground state (<H> -0.987 -> -1.000)"""
    np.random.seed(2)
    sc = heisenberg(version)
    s = sc.get_gs()
    x = s + 0.3*(sc.Sx[0]*s)
    x = x*(1/np.sqrt(x.dot(x).real))
    keep = x.copy()
    e_x = energy(sc, x)
    sc.computed_gs = False # the only way the explicit route used to be read
    sc.gs_energy(wf0=x) # the explicit warm start: sweeps from x
    assert sc.gs_energy() < e_x - 1e-3
    sc.set_initial_wf_guess(x); sc.gs_energy()
    sc.set_initial_wf(x); sc.gs_energy()
    sc.set_gs(x)
    sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=ES, delta=DELTA)
    assert energy(sc, x) == pytest.approx(e_x, abs=1e-12)
    assert fidelity(x, keep) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("version", BACKENDS)
def test_explicit_wf0_is_read_on_a_current_chain(version):
    """Many_Body_Chain.gs_energy returned the stored energy before reading
    its keywords, so gs_energy(wf0=x) on a solved chain ignored x"""
    np.random.seed(1)
    sc = heisenberg(version)
    Szt, up, dn = members(sc)
    target = dn if np.real(sc.vev(Szt)) > 0 else up
    sc.gs_energy(wf0=target)
    assert fidelity(sc.get_gs(), target) == pytest.approx(1.0, abs=1e-8)


@pytest.mark.parametrize("version", BACKENDS)
def test_transverse_ising_warm_start_seeds_the_branch(version):
    """examples/topological/transverse_ising_model: the ferromagnetic guess
    seeds the symmetry-broken branch in the ordered phase. Without the
    hand-off the chain returned the symmetric state, Mz/n ~ 0.001 on the
    example's 40 sites where the seeded branch gives 0.4998. Sixteen sites
    at B=0.1 put the tunnelling split of the two branches near 1e-11, far
    below what a sweep resolves; at ten sites and B=0.2 it is about 1e-5
    and a sweep from the seed finds the symmetric state on its own."""
    np.random.seed(1)
    n = 16
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 8
    Mz = sum(sc.Sz)
    sc.set_hamiltonian(-1*Mz)
    wffe = sc.get_gs().copy()
    h = 0
    for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): h = h + 0.1*sc.Sx[i]
    sc.set_hamiltonian(h)
    sc.set_initial_wf_guess(wffe)
    assert np.real(sc.vev(Mz))/n > 0.45


# ------------------------------------------------------ findings 13, 14, 15

D, TP = 1e-3, 2*np.pi
EVS = np.array([-0.5e-3, 0.5e-3])
KONDO_ES = np.linspace(-30e-3, 30e-3, 3001)


def crossing_chain(version, eps):
    """The previous record's finding 12 chain: an S=1 impurity with
    D*Sz^2 + (D+eps)*Sz, |0> and |-1> split by eps, next to two
    field-polarized S=1/2 spectators"""
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=version)
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + (D + eps)*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc


def kondo_n_gs(sc, submode="KPM", **extra):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, d = sc.get_kondo_spectrum(EVS, site=0, T=0.0, order=2, mode="DMRG",
                                     submode=submode, delta=1e-4, es=KONDO_ES,
                                     n_gs=2, **extra)
    return d/TP


@pytest.mark.parametrize("version", [_backend(3), _backend(2)])
def test_n_gs_measures_each_member_unswept_and_restores_the_chain(version):
    """eps=1e-5 < delta=1e-4, so the split warning stays silent; the two-
    state average is 1.5 and the lower state alone gives 2.0. On v3 the
    upper member relaxed onto the lower one in 5 of 6 runs of this loop
    before the fix, and gs_energy() afterwards was the last member's"""
    for run in range(6):
        np.random.seed(100 + run)
        sc = crossing_chain(version, 1e-5)
        e_before = sc.gs_energy()
        wf_before = sc.get_gs()
        d = kondo_n_gs(sc)
        assert np.allclose(d, 1.5, atol=1e-3), (run, d)
        assert sc.gs_energy() == e_before # restored as a unit, bit for bit
        assert sc.get_gs() is wf_before


def test_n_gs_gs_energy_is_each_members_own_inside_the_loop(monkeypatch):
    """Finding 14's other half: while a member is measured, gs_energy() is
    that member's own energy, the convention the ED references use"""
    import dmrgpy.kondospectrumtk.secondorder_dc as secondorder_dc
    seen = []
    def record(chain, *a, **k):
        seen.append((chain.gs_energy(), energy(chain, chain.get_gs())))
        return np.zeros(len(EVS))
    monkeypatch.setattr(secondorder_dc, "second_order_dIdV_dc", record)
    np.random.seed(0)
    sc = crossing_chain("python", 1e-5)
    e_before = sc.gs_energy()
    kondo_n_gs(sc)
    assert len(seen) == 2
    for e_seen, e_member in seen:
        assert e_seen == pytest.approx(e_member, abs=1e-14)
    assert max(e for e, _ in seen) - min(e for e, _ in seen) == pytest.approx(1e-5, rel=1e-6)
    assert sc.gs_energy() == e_before


def test_n_gs_averages_under_ex():
    """EX measured from the lowest vector of its own cached basis, so n_gs=2
    returned the n_gs=1 value, anywhere in [1.1, 1.9] from run to run at the
    exact crossing. Measured from the chain's state reprojected onto that
    basis it is the average, the same on every run (1.512780: EX's
    Lorentzian lets the S+ pole at 2e-3 leak into the plateau)"""
    vals = []
    for seed in (100, 101):
        np.random.seed(seed)
        sc = crossing_chain("python", 0.0)
        sc.gs_energy()
        vals.append(kondo_n_gs(sc, submode="EX", nex=10))
    assert np.allclose(vals[0], 1.5, atol=0.02)
    assert np.allclose(vals[0], vals[1], rtol=0., atol=1e-6)


def test_ex_at_a_unique_ground_state_is_unchanged(monkeypatch):
    """At a unique ground state the chain's state is the basis's lowest
    vector, so the reprojection gives what the shipped prescription gave
    (the lowest vector itself, measured from its own energy), on the same
    cached basis. Three delta off the crossing, where EX gives 0.9287 and
    ED 1."""
    from dmrgpy import dcex
    np.random.seed(51)
    sc = crossing_chain("python", -3e-4)
    kw = dict(site=0, T=0.0, order=2, mode="DMRG", submode="EX", delta=1e-4,
              es=KONDO_ES, nex=10)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, new = sc.get_kondo_spectrum(EVS, **kw)
        monkeypatch.setattr(dcex, "_reference_coefficients",
                            lambda self, ws, C: np.eye(C.shape[1])[0])
        _, old = sc.get_kondo_spectrum(EVS, **kw) # same cached basis
    assert np.allclose(new/TP, 0.9287, atol=1e-3)
    assert np.max(np.abs(new - old)) < 1e-12*np.max(np.abs(old))


def test_ex_refuses_a_state_outside_its_basis():
    np.random.seed(4)
    sc = heisenberg("python", n=4)
    kw = dict(submode="EX", name=(sc.Sz[0], sc.Sz[0]), es=ES, delta=DELTA, nex=2)
    sc.get_dynamical_correlator(**kw) # builds the two-state basis
    s = sc.get_gs()
    x = s + 0.5*(sc.Sx[0]*s)
    sc.set_gs(x*(1/np.sqrt(x.dot(x).real)))
    with pytest.raises(ValueError, match="outside"):
        sc.get_dynamical_correlator(**kw)


@pytest.mark.parametrize("submode", ["SECTOR", "maxent"])
def test_n_gs_refuses_a_submode_that_does_not_read_the_state(submode):
    sc = crossing_chain("python", 0.0)
    with pytest.raises(NotImplementedError, match="n_gs"):
        kondo_n_gs(sc, submode=submode)


# ------------------------------------------- what the removed line did

@pytest.mark.parametrize("version", BACKENDS)
def test_a_clone_of_a_solved_chain_runs_a_correlator(version):
    """A clone used to keep computed_gs=True and a wf0 on the original's
    site set, with an empty session; without the removed line a correlator
    then raised "called before set_hamiltonian" on python and aborted the
    process on v3"""
    np.random.seed(5)
    sc = heisenberg(version, n=4) # a unique singlet ground state
    sc.get_gs()
    y0 = np.asarray(sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                                es=ES, delta=DELTA)[1])
    cl = sc.clone()
    assert not cl.computed_gs and cl.wf0 is None
    y1 = np.asarray(cl.get_dynamical_correlator(name=(cl.Sz[0], cl.Sz[0]),
                                                es=ES, delta=DELTA)[1])
    assert np.max(np.abs(y1 - y0)) < 5e-2*np.max(np.abs(y0))
    assert cl.gs_energy() == pytest.approx(sc.gs_energy(), abs=1e-8)


@pytest.mark.parametrize("version", BACKENDS)
@pytest.mark.parametrize("submode,extra", [("KPM", {}), ("TD", dict(dt=0.1))],
                         ids=["KPM", "TD"])
def test_a_repeated_correlator_on_a_solved_chain_does_not_resweep(version,
                                                                   submode, extra):
    np.random.seed(6)
    sc = heisenberg(version, n=4)
    kw = dict(submode=submode, name=(sc.Sz[0], sc.Sz[0]), es=ES, delta=DELTA, **extra)
    y0 = np.asarray(sc.get_dynamical_correlator(**kw)[1])
    spy = spy_on(sc)
    sweeps = []
    if version == "python": # count every DMRG sweep the session runs itself
        import dmrgpy.pyitensor.chain as pychain
        real = pychain.dmrg
        def counting(*a, **k):
            sweeps.append(1)
            return real(*a, **k)
        pychain.dmrg = counting
    try:
        y1 = np.asarray(sc.get_dynamical_correlator(**kw)[1])
    finally:
        if version == "python": pychain.dmrg = real
    assert spy.log == []
    assert sweeps == []
    assert np.allclose(y0, y1, rtol=0., atol=1e-12)
