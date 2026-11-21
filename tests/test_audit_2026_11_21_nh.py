"""Regressions for the "nh" cluster of docs/audit_2026_11_21_hole_hunt.md:
the level of smallest real part of a non-Hermitian Hamiltonian.

- finding 11: the default gs_energy() on DMRG (nhlevels.ground_pair) took
  the best converged Ritz pair within 1e-4*(c+|E-e_id|) of the smallest
  real part, so with a field of 5e-5 on the 6-site PT-symmetric XX chain it
  returned the conjugate pair 5e-5 above the lowest level, or the level
  1e-4 above it, with <Sz_total> 0 or +1 for -1;
- finding 12: ED ordered non-Hermitian eigenvalues by Re + 1e-3 Im, and
  returned a level of higher real part whenever its Im was lower by more
  than 1e3 times the difference;
- finding 13: algebra.biorthogonal_ground_state sorted h and dagger(h)
  with that one key, so on a smallest real part shared by a conjugate
  pair the left state belonged to the other member, and it raised, and
  the ED non-Hermitian KPM correlator with it.

Every anchor is np.linalg.eig or scipy.linalg.eig of the dense ED matrix.
The pre-fix orderings are kept inside the tests to show that each model
separates the two rules."""

import random

import numpy as np
import pytest
import scipy.linalg as sla
from scipy.sparse import csr_matrix

from dmrgpy import cppext, nhlevels, spinchain
from dmrgpy.algebra import algebra
from dmrgpy.algebra import kpm as kpmalg

VERSIONS = [
    pytest.param("python", id="python"),
    pytest.param(3, id="v3", marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
]

G = 0.3


def _old_key_index(ev):
    """The pre-fix ordering of algebra.sorteigen: Re(E - min Re) + 1e-3 Im."""
    w = ev - np.min(ev.real)
    return int(np.argmin(w.real + 1e-3 * w.imag))


def _pt_chain(n, eps, version="python"):
    """PT-symmetric XX chain, i G (-1)^i Sz_i, plus eps*Sz_total: at eps=0
    four levels share the smallest real part (E0 twice, E0 +- i gamma);
    the field splits them to E0-eps (Sz_total=-1), E0 +- i gamma (0) and
    E0+eps (+1)."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=version)
    sc.maxm = 16
    sc.nsweeps = 12
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
    for i in range(n):
        h = h + 1j * G * (-1)**i * sc.Sz[i]
    obs = 0
    for i in range(n):
        obs = obs + sc.Sz[i]
    if eps != 0.:
        h = h + eps * obs
    sc.set_hamiltonian(h)
    return sc, obs


def _impurity(n, eps, g):
    """-sum_{i>=1} Sz_i + (eps - i g) Sz_0: the smallest real part is
    -(n-1)/2 - eps/2 + i g/2, with Sz_0 = -1/2, and the level eps above it
    has Im lower by g."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
    h = 0
    for i in range(1, n):
        h = h - 1.0 * sc.Sz[i]
    h = h + (eps - 1j * g) * sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc, sc.Sz[0]


def _heis_pair(n=4, gamma=1.0):
    """Heisenberg chain with 0.3 Sz_total + i gamma Sz_0: the smallest real
    part is shared by the pair -1.242697 +- 0.291419i."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = (h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
             + sc.Sz[i] * sc.Sz[i + 1])
    for i in range(n):
        h = h + 0.3 * sc.Sz[i]
    h = h + 1j * gamma * sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc


def _dense(sc, op=None):
    ed = sc.get_ED_obj()
    m = ed.get_hamiltonian() if op is None else ed.get_operator(op)
    return m.toarray() if hasattr(m, "toarray") else np.array(m)


def _exact_lowest(hm, om):
    """(E, <O>) of the eigenvalue of smallest real part, unique here."""
    ev, vr = np.linalg.eig(hm)
    k = int(np.argmin(ev.real))
    v = vr[:, k]
    return ev, k, float(np.real(np.vdot(v, om @ v) / np.vdot(v, v)))


def _tied_lowest(w):
    """The index of smallest real part, the smallest Im among the real
    parts tied with it to 1e-9: the member of a conjugate pair the ordering
    takes."""
    tied = np.flatnonzero(w.real <= w.real.min() + 1e-9)
    return int(tied[np.argmin(w[tied].imag)])


@pytest.fixture
def no_fallback(monkeypatch):
    """Fail the test if the imaginary-time route hands over to NH-DMRG,
    which returns the conjugate pair on this chain in every run."""
    from dmrgpy import nhdmrg
    def refuse(*args, **kwargs):
        raise AssertionError("the imaginary-time route fell back to NH-DMRG")
    monkeypatch.setattr(nhdmrg, "gs_energy_nhdmrg", refuse)
    monkeypatch.setattr(nhdmrg, "nhdmrg", refuse)


# ------------------------------------------------------------ finding 11

def _pairs(e0, rows):
    """Accepted pairs (theta, absolute residual, c) at rscale=1, e_id=0,
    from (Re - e0, Im, relative residual) rows."""
    out = []
    for j, (dre, im, rel) in enumerate(rows):
        t = complex(e0 + dre, im)
        out.append((t, rel * (1.0 + abs(t)), j))
    return out


E0 = -1.39617641
# the accepted lists of two runs of the record's counter-repro (seeds 0 and
# 1 on "python"): the lowest level first, then the pair 5e-5 above it and
# the level 1e-4 above it, with their relative residuals
RUN_SEED0 = [(-1.371e-6, 0.0, 2.72e-6), (5.002e-5, -0.20121, 2.18e-7),
             (5.002e-5, 0.20121, 9.93e-8), (1.005e-4, 0.0, 1.75e-5)]
RUN_SEED1 = [(-2.210e-12, 0.0, 1.95e-7), (5.000e-5, -0.20121, 6.18e-9),
             (5.000e-5, 0.20121, 1.24e-8), (1.000e-4, 0.0, 0.0)]


@pytest.mark.parametrize("rows", [RUN_SEED0, RUN_SEED1], ids=["seed0", "seed1"])
@pytest.mark.parametrize("far", [False, True], ids=["accepted-only", "with-far-ritz"])
def test_selection_takes_the_smallest_real_part_among_converged_pairs(rows, far):
    """The rule of the fix takes the lowest level; the pre-fix window
    1e-4*(c+|E-e_id|) took the best converged pair of the four, a
    different level, on the same list."""
    acc = _pairs(E0, rows)
    ritz = [p[0] for p in acc] + ([E0 + 1.0, E0 + 2.5 - 0.4j] if far else [])
    pick, low = nhlevels._select_ground(acc, ritz, 1.0, 0.0, 1e-4)
    assert pick[2] == 0
    assert low is None
    # pre-fix: every pair in the window, then the smallest residual
    w = 1e-4 * (1.0 + abs(acc[0][0]))
    old = min([p for p in acc if p[0].real <= acc[0][0].real + w],
              key=lambda p: p[1])
    assert old[2] != 0


def test_selection_breaks_a_true_tie_by_residual_and_names_a_skipped_level():
    """Two Ritz values of one level (1e-12 apart) go to the better
    converged one; a pair below by 1e-3 whose residual is above the
    certificate tolerance is skipped and returned for the warning."""
    rows = [(-1e-3, 0.0, 5e-4), (0.0, 0.0, 3e-6), (1e-12, 0.0, 2e-9),
            (0.05, 0.1, 1e-9)]
    acc = _pairs(E0, rows)
    ritz = [p[0] for p in acc] + [E0 + 1.0]
    pick, low = nhlevels._select_ground(acc, ritz, 1.0, 0.0, 1e-4)
    assert pick[2] == 2
    assert low is not None and low[2] == 0


def test_selection_when_no_pair_reaches_the_certificate_tolerance():
    """cert_tol below every residual (gs_energy(tol=1e-30)): the best
    converged pair within 1e-4*(c+|E-e_id|) of the smallest real part,
    which the certificate and NH-DMRG then judge, and no warning."""
    acc = _pairs(E0, RUN_SEED0)
    pick, low = nhlevels._select_ground(acc, [p[0] for p in acc], 1.0, 0.0,
                                        1e-30)
    assert pick[2] == 2
    assert low is None


def test_tie_window_is_shift_invariant_and_scale_covariant():
    ev = np.array([E0, E0 + 5e-5 + 0.2j, E0 + 1.0])
    w = nhlevels._tie_window(ev)
    assert 1e-6 <= w < 2e-6
    assert nhlevels._tie_window(1e-6 * ev) == pytest.approx(1e-6 * w, rel=1e-6)
    assert nhlevels._tie_window(ev + 10.0) == pytest.approx(w, abs=1e-12)


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("version", VERSIONS)
def test_default_nh_gs_energy_is_the_level_of_smallest_real_part(version, seed,
                                                                 no_fallback):
    """The 6-site PT chain in a field 5e-5: E and <Sz_total> = -1 of the
    level of smallest real part, against np.linalg.eig. Before the fix 3
    of 4 seeded runs on "python" and 3 of 3 on v3 returned the pair 5e-5
    above it or the level 1e-4 above it."""
    np.random.seed(seed)
    random.seed(seed)
    sc, obs = _pt_chain(6, 5e-5, version)
    e = complex(sc.gs_energy())
    o = float(np.real(sc.vev(obs)))
    ev, k, o_true = _exact_lowest(_dense(sc), _dense(sc, obs))
    assert o_true == pytest.approx(-1.0, abs=1e-8)
    assert int(np.argmin(np.abs(ev - e))) == k
    assert abs(e - ev[k]) < 1e-5
    assert o == pytest.approx(o_true, abs=1e-3)
    assert sc._nh_last_residual[0] < 1e-4


# ------------------------------------------------------------ finding 12

ED_MODELS = {
    "impurity": lambda: _impurity(4, 1e-4, 0.5),
    "pt-field": lambda: _pt_chain(6, 5e-5),
}


@pytest.mark.parametrize("model", list(ED_MODELS))
def test_ed_ground_state_is_the_smallest_real_part(model):
    """gs_energy(mode="ED") and vev(mode="ED") on the level of smallest
    real part, where the pre-fix key took the one of lower Im above it."""
    sc, obs = ED_MODELS[model]()
    hm, om = _dense(sc), _dense(sc, obs)
    ev, k, o_true = _exact_lowest(hm, om)
    assert _old_key_index(ev) != k  # the model separates the two orderings
    e = complex(sc.gs_energy(mode="ED"))
    o = float(np.real(sc.vev(obs, mode="ED")))
    assert abs(e - ev[k]) < 1e-10
    assert o == pytest.approx(o_true, abs=1e-8)


def test_lowest_states_on_the_arpack_branch(monkeypatch):
    """The ARPACK branch of lowest_states and lowest_eigenvalues (above
    algebra.maxsize) orders the same way."""
    sc, obs = _pt_chain(6, 5e-5)
    hm = _dense(sc)
    ev, k, _ = _exact_lowest(hm, _dense(sc, obs))
    assert _old_key_index(ev) != k
    monkeypatch.setattr(algebra, "maxsize", 8)
    es, vs = algebra.lowest_states(csr_matrix(hm), n=1)
    assert abs(es[0] - ev[k]) < 1e-10
    v = np.asarray(vs[0]).reshape(-1)
    assert np.linalg.norm(hm @ v - ev[k] * v) < 1e-8 * np.linalg.norm(v)
    assert abs(algebra.lowest_eigenvalues(csr_matrix(hm), n=1)[0] - ev[k]) < 1e-10


@pytest.mark.parametrize("scale,offset", [(1.0, 0.0), (1e-8, 0.0), (1.0, 1e3)],
                         ids=["unit", "small-units", "offset"])
def test_sorteigen_orders_by_real_part_and_im_only_inside_a_tie(scale, offset):
    """Real part first, a tie at 1e-10 of the largest |E| ordered by Im,
    and any sequence of vectors carried along as a list (excited.py passes
    a list of MPS)."""
    eig = np.array([0.3 + 0j, -1.5 + 1e-4 - 0.25j, -1.5 + 0.25j,
                    -1.5 + 1e-14 + 0.1j]) * scale + offset
    labels = ["d", "a", "b", "c"]
    es, vs = algebra.sorteigen(eig, labels)
    assert isinstance(vs, list) and isinstance(es, np.ndarray)
    assert vs == ["c", "b", "a", "d"]
    assert np.array_equal(es, eig[[3, 2, 1, 0]])
    if scale == 1.0 and offset == 0.0:
        assert labels[_old_key_index(eig)] == "a"  # the pre-fix first


# ------------------------------------------------------------ finding 13

BIORTH_MODELS = {
    "pt-tied": lambda: _pt_chain(6, 0.0)[0],
    "heis4": _heis_pair,
}


def _check_pair(hm, e0, vr, vl):
    w = np.linalg.eigvals(hm)
    assert e0.real == pytest.approx(w.real.min(), abs=1e-10)
    nr, nl = np.linalg.norm(vr), np.linalg.norm(vl)
    assert np.linalg.norm(hm @ vr - e0 * vr) < 1e-8 * nr
    assert np.linalg.norm(hm.conj().T @ vl - np.conj(e0) * vl) < 1e-8 * nl
    assert np.conj(vl) @ vr == pytest.approx(1.0, abs=1e-10)
    return abs(np.conj(vl) @ vr) / (nr * nl)


@pytest.mark.parametrize("arpack", [False, True], ids=["dense", "arpack"])
@pytest.mark.parametrize("model", list(BIORTH_MODELS))
def test_biorthogonal_ground_state_on_a_conjugate_pair(model, arpack, monkeypatch):
    """The right and left states of one level, with the overlap scipy's
    pair of that level has (0.55 and 0.62), where the pre-fix pair raised:
    its left state belonged to the other member of the pair."""
    sc = BIORTH_MODELS[model]()
    hm = _dense(sc)
    ev, evd = np.linalg.eigvals(hm), np.linalg.eigvals(hm.conj().T)
    # pre-fix: each side by the old key, the left one on conj(E)
    assert abs(np.conj(evd[_old_key_index(evd)]) - ev[_old_key_index(ev)]) > 0.1
    if arpack:
        monkeypatch.setattr(algebra, "maxsize", 8)
    e0, vr, vl = algebra.biorthogonal_ground_state(csr_matrix(hm))
    cos = _check_pair(hm, e0, vr, vl)
    assert cos > 0.1
    w, VL, VR = sla.eig(hm, left=True, right=True)
    k = _tied_lowest(w)
    assert abs(e0 - w[k]) < 1e-8
    ref = abs(np.conj(VL[:, k]) @ VR[:, k]) / (
        np.linalg.norm(VL[:, k]) * np.linalg.norm(VR[:, k]))
    assert cos == pytest.approx(ref, abs=1e-6)


def test_biorthogonal_ground_state_on_a_degenerate_level():
    """A twofold level of h = Q (A x 1) Q^dagger: vl is the left vector of
    the whole eigenspace of largest overlap with vr, whose overlap is the
    one of A's own pair, whatever vector of the eigenspace vr is."""
    rng = np.random.default_rng(7)
    a = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
    q, _ = np.linalg.qr(rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12)))
    hm = q @ np.kron(a, np.eye(2)) @ q.conj().T
    e0, vr, vl = algebra.biorthogonal_ground_state(csr_matrix(hm))
    cos = _check_pair(hm, e0, vr, vl)
    w, VL, VR = sla.eig(a, left=True, right=True)
    k = int(np.argmin(w.real))
    assert abs(e0 - w[k]) < 1e-8
    ref = abs(np.conj(VL[:, k]) @ VR[:, k]) / (
        np.linalg.norm(VL[:, k]) * np.linalg.norm(VR[:, k]))
    assert cos == pytest.approx(ref, abs=1e-8)


@pytest.mark.parametrize("model", list(BIORTH_MODELS))
def test_ed_nh_kpm_correlator_on_a_conjugate_pair(model):
    """get_dynamical_correlator(mode="ED") on a non-Hermitian chain (the
    NH-KPM submode) answers, and equals the same recursion run on scipy's
    biorthogonal pair of the level of smallest real part."""
    sc = BIORTH_MODELS[model]()
    es = np.linspace(-0.5, 3.0, 20)
    delta, emax, n = 0.1, 20.0, 50
    x, y = sc.get_dynamical_correlator(mode="ED", name=(sc.Sz[1], sc.Sz[1]),
                                       delta=delta, es=es, E_max=emax, n=n)
    y = np.asarray(y)
    assert np.all(np.isfinite(y))
    ed = sc.get_ED_obj()
    h = ed.get_hamiltonian()
    a = ed.get_operator(sc.Sz[1])
    hm = h.toarray() if hasattr(h, "toarray") else np.array(h)
    w, VL, VR = sla.eig(hm, left=True, right=True)
    k = _tied_lowest(w)
    vr, vl = VR[:, k], VL[:, k]
    vl = vl / np.conj(np.conj(vl) @ vr)
    wfa = a @ vr
    wfb = algebra.dagger(a) @ vl
    ref = np.array([kpmalg.get_spec_kpm_nh(h, w[k] + z + 1j * delta, wfa, wfb,
                                           emax, n=n, kernel="jackson")
                    for z in es])
    assert np.max(np.abs(y - ref)) < 1e-8 * np.max(np.abs(ref))
