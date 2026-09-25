"""Regression tests for the `nh` cluster of the 2026-09-25b hole hunt,
findings 25, 26 and 27 (the NH-DMRG halves).

25. NH-DMRG's right solve picks its Ritz value through SRTieBreak, which
    treats every value with Re <= min Re + degtol as tied and follows the
    previous bond among them. degtol was 1e-6*(1+|min Re|): absolute at
    small |E| (1e-6/s wide in units of H for s*H) and growing with |E| (a
    constant offset c makes it 1e-6*c), so once it exceeded the real-part
    gap the sweep came back on a converged EXCITED eigenpair, which no
    residual certificate can reject. It is now measured against the Ritz
    spread, 1e-6*(max Re - min Re), with a roundoff floor, which is shift
    invariant and scale covariant (pyitensor/nhdmrg.py::_select_ritz; the
    C++ arnoldi_select_kbest follows in the cpp cluster's rebuild).
26. nhdmrg.py's eigen-residual certificate divided by 1 + |E|, which is
    neither scale covariant nor offset invariant, so below s ~ 2e-5, or
    next to an offset of ~1e4, every state passed it. It divides by
    c + |E - e_id| now (c = min(1, largest non-identity |coefficient|),
    e_id = summed identity coefficient), which is the old test for any H
    with no identity term and a largest coefficient of 1 or more.
27. The local Arnoldi's breakdown and restart tests are absolute, and the
    NH entry points handed the session the caller's own units. nhdmrg()
    and nhdmrg_generalized() now hand the session 2^k*H, the power of two
    mo_terms.h's unit_scale_up() uses, and divide the energy back.

The anchors are ED, the same calculation at scale 1, and linearity: no
golden number enters.
"""

import contextlib
import io
import math

import numpy as np
import pytest

from dmrgpy import cppext, nhdmrg, spinchain
from dmrgpy.pyitensor import nhdmrg as pnh


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


L = 6


def _ham(sc, n=L):
    """Heisenberg + 0.3j*Sz0: a real smallest-real-part level separated
    from the next one by a real-part gap of 0.4488 on 6 sites, 0.2911 on
    10, with a complex-conjugate pair right above it (the hunt's chain)."""
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h + 0.3j*sc.Sz[0]


_SPECTRA = {}


def _spectrum(n=L):
    """ED eigenvalues of _ham, sorted by real part."""
    if n not in _SPECTRA:
        ref = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
        m = ref.get_ED_obj().get_operator(_ham(ref, n))
        m = m.toarray() if hasattr(m, "toarray") else np.asarray(m)
        ev = np.linalg.eigvals(m)
        _SPECTRA[n] = ev[np.argsort(ev.real)]
    return _SPECTRA[n]


def _chain(version, n=L, maxm=30, nsweeps=10):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.maxm = maxm
    sc.nsweeps = nsweeps
    return sc


def _quiet(f, *args, **kwargs):
    """f(*args, **kwargs) and whatever it printed."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = f(*args, **kwargs)
    return out, buf.getvalue()


# ---------------------------------------------------------------- 25 --

def _old_select(evals, target):
    """The SRTieBreak window before the fix, kept as the reference the
    tests below must fail."""
    remin = evals.real.min()
    cand = np.flatnonzero(evals.real < remin + 1e-6*(1.0 + abs(remin)))
    return int(cand[np.argmin(np.abs(evals[cand] - target))])


@pytest.mark.parametrize("s,c", [(1e-6, 0.0), (1e-7, 0.0), (1e-12, 0.0),
                                 (1.0, 1e6), (1e-3, 1e5)])
def test_select_ritz_window_is_scale_covariant_and_shift_invariant(s, c):
    """The ED spectrum of the hunt's chain in units s and shifted by c, as
    a Ritz set, with the previous bond's value on level #3 (the case the
    sweep follows): the tie window must still exclude every level above
    the gap. The old window let levels #1..#3 in at s=1e-6 and c=1e6."""
    ev = _spectrum()[:12]
    ritz = s*ev + c
    target = s*ev[3] + c
    assert pnh._select_ritz(ritz, "SRTieBreak", target) == 0
    assert pnh._select_ritz(ev, "SRTieBreak", ev[3]) == 0 # unit scale
    if (s, c) in ((1e-6, 0.0), (1.0, 1e6)):
        assert _old_select(ritz, target) != 0


def test_select_ritz_lone_and_equal_values_stay_candidates():
    """A single Ritz value, an all-equal set and an all-zero set (the zero
    operator) each leave a candidate: the window is compared with <=."""
    assert pnh._select_ritz(np.array([2.0+0j]), "SRTieBreak", 5.0) == 0
    assert pnh._select_ritz(np.zeros(3, dtype=complex), "SRTieBreak", 0j) in (0, 1, 2)
    ev = np.array([1.0+0.5j, 1.0-0.5j, 1.0+0j])
    assert pnh._select_ritz(ev, "SRTieBreak", 1.0-0.4j) == 1


def test_select_ritz_roundoff_floor_keeps_a_split_degenerate_pair():
    """A Re-degenerate conjugate pair at a large offset comes out of the
    Hessenberg eigensolver split by roundoff of order eps*|E|, far above
    1e-6 of a spread of order 1: the floor keeps both members tied, so the
    sweep keeps following the previous bond's member instead of flipping
    between them (the reviewer's caveat 1)."""
    c = 1e12
    eps = np.finfo(float).eps
    ev = np.array([c - 1.0 + 0.1j, c - 1.0 + 10*eps*c - 0.1j, c + 0.5])
    assert ev[1].real > ev[0].real # the split survived
    assert pnh._select_ritz(ev, "SRTieBreak", ev[1]) == 1


def test_python_nh_session_small_units_lands_on_ground_level():
    """The window alone, at small units: the "python" session called with
    s*H itself (no unit scale at the entry), seeded. Before the fix this
    returned level #6, 1.06 off, at s=1e-7."""
    s = 1e-7
    ev = _spectrum()
    np.random.seed(7)
    sc = _chain("python")
    H = s*_ham(sc)
    sc.set_hamiltonian(H)
    sc._session.set_sweep_params(sc.maxm, sc.nsweeps, sc.cutoff, sc.noise)
    sc._session.set_mpomaxm(max(sc.maxm, sc.mpomaxm))
    e, _, _ = sc._session.nhdmrg(H.to_terms(), H.get_dagger().to_terms(), 20, 2)
    assert abs(complex(e)/s - ev[0]) < 1e-8


def test_python_nhdmrg_offset_lands_on_ground_level():
    """The window alone, at a constant offset, which no unit scale reaches
    (the largest coefficient is the offset): E0 - c must be the ground
    level. Before the fix, level #3, 0.482 off, at c=1e6. The tolerance
    is "python"'s recorded MPO truncation of an offset, 1.8/c, far below
    the 0.4488 gap."""
    c = 1e6
    ev = _spectrum()
    np.random.seed(7)
    sc = _chain("python")
    sc.set_hamiltonian(_ham(sc) + c)
    (e, _, _), _ = _quiet(nhdmrg.nhdmrg, sc)
    assert abs(complex(e) - c - ev[0]) < 1e-4


@pytest.mark.parametrize("version", [2, 3])
def test_cpp_nhdmrg_offset_lands_on_ground_level(version):
    """The same offset on v2 and v3, whose window is arnoldi_select_kbest's
    (the cpp cluster's half of finding 25): before, v2 at c=1e6 returned an
    excited level in 9 of 9 runs in the record. The unit scale at the NH
    entry is a no-op here (the largest coefficient is c), so this pins the
    C++ window alone, through the public driver."""
    if not cppext.available(version):
        pytest.skip("needs the compiled itensor_version=%d extension" % version)
    c = 1e6
    ev = _spectrum()
    sc = _chain(version)
    sc.set_hamiltonian(_ham(sc) + c)
    (e, _, _), _ = _quiet(nhdmrg.nhdmrg, sc)
    assert abs(complex(e) - c - ev[0]) < 1e-4


# ------------------------------------------------------------ 25, 27 --

@pytest.mark.parametrize("version,s", [
    (2, 1e-7), ("python", 1e-7), (3, 1e-13), (3, 1e-16), ("python", 1e-13)],
    ids=["v2-1e-7", "python-1e-7", "v3-1e-13", "v3-1e-16", "python-1e-13"])
def test_nh_gs_energy_is_scale_covariant(version, s):
    """E0(s*H)/s = E0(H) through the public gs_energy() of a non-Hermitian
    chain. Before: v2 1.47 to 1.85 off at 1e-7 (the SRTieBreak window of
    arnoldi_select_kbest, which the unit scale puts back at its
    calibrated scale), v3 0.24 to 2.6 off from 1e-13 (the Arnoldi's
    absolute breakdown and restart tests), "python" 1.06 off at 1e-7 and
    2.4 at 1e-13."""
    if version in (2, 3) and not cppext.available(version):
        pytest.skip("needs the compiled itensor_version=%d extension" % version)
    ev = _spectrum()
    if version == "python": np.random.seed(7)
    sc = _chain(version)
    sc.set_hamiltonian(s*_ham(sc))
    e0, _ = _quiet(sc.gs_energy)
    assert abs(complex(e0)/s - ev[0]) < 1e-8


@pytest.mark.skipif(not cppext.available(3), reason="needs the v3 extension")
@pytest.mark.parametrize("s", [2e-6, 1e-7])
def test_nhdmrg_v3_below_full_bond_dimension_small_units(s):
    """The default backend below full bond dimension (10 sites at maxm=8),
    where v3 did not escape the window: 0.29 to 0.32 off at s=2e-6 and 1.8
    to 2.9 off at 1e-7. Now the s=1 answer, the maxm=8 truncation value
    5.46e-5 above ED, against a real-part gap of 0.2911."""
    ev = _spectrum(10)
    sc = _chain(3, n=10, maxm=8)
    sc.set_hamiltonian(s*_ham(sc, 10))
    (e, _, _), _ = _quiet(nhdmrg.nhdmrg, sc, ntries=1)
    assert abs(complex(e)/s - ev[0]) < 1e-3


def _generalized_reference():
    """lambda0 of _ham against A = 1 + 0.2*Sz0, by dense linear algebra."""
    from scipy.linalg import eig
    ref = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
    hm = ref.get_ED_obj().get_operator(_ham(ref))
    am = ref.get_ED_obj().get_operator(1 + 0.2*ref.Sz[0])
    hm = hm.toarray() if hasattr(hm, "toarray") else np.asarray(hm)
    am = am.toarray() if hasattr(am, "toarray") else np.asarray(am)
    w = eig(hm, am, right=False)
    return w[np.argmin(w.real)]


@pytest.mark.parametrize("version", [_backend(3), _backend("python")])
def test_nh_generalized_is_scale_covariant(version):
    """lambda0(s*H, A)/s = lambda0(H, A) through gs_energy_generalized with
    A = 1 + 0.2*Sz0, at s=1e-13 where the local Arnoldi's absolute tests
    bit (v3 and "python" alike)."""
    lam_ref = _generalized_reference()
    s = 1e-13
    if version == "python": np.random.seed(7)
    sc = _chain(version)
    sc.set_hamiltonian(s*_ham(sc))
    lam, _ = _quiet(sc.gs_energy_generalized, 1 + 0.2*sc.Sz[0])
    assert abs(complex(lam)/s - lam_ref) < 1e-8


def test_nh_generalized_carries_lam0_into_the_solver_units(monkeypatch):
    """What the session receives: H and H^dagger at 2^k (A untouched) and
    the caller's lam0 in the same units, 2^k*lam0; and lambda is handed
    back in the caller's units. A lam0 left unscaled would start the
    self-consistency at a shift 2^k times too small."""
    s = 1e-7
    up = 2.0**24 # unit_scale_up(1e-7)
    np.random.seed(7)
    sc = _chain("python")
    sc.set_hamiltonian(s*_ham(sc))
    A = 1 + 0.2*sc.Sz[0]
    lam0 = s*(-2.5 - 0.1j)
    seen = {}
    orig = sc._session.nhdmrg_generalized
    def spy(th, thd, ta, krylovdim, restarts, lam0=None):
        seen.update(th=th, ta=ta, lam0=lam0)
        return orig(th, thd, ta, krylovdim, restarts, lam0=lam0)
    monkeypatch.setattr(sc._session, "nhdmrg_generalized", spy)
    (lam, _, _), _ = _quiet(nhdmrg.nhdmrg_generalized, sc, A, lam0=lam0,
                            ntries=1)
    assert seen["lam0"] == lam0*up
    assert [c for c, _ in seen["th"]] == [c*up for c, _ in (s*_ham(sc)).to_terms()]
    assert seen["ta"] == A.to_terms()
    assert abs(complex(lam)/s - _generalized_reference()) < 1e-8


@pytest.mark.parametrize("version", [_backend(3), _backend(2), _backend("python")])
def test_nh_solve_hands_back_the_callers_units(version):
    """The session solves 2^k*H; nothing the caller reads may carry the
    2^k. e0 must be the biorthogonal Rayleigh quotient of the stored pair
    under H itself, applied by the session afterwards from the chain's
    own, unscaled, terms."""
    s = 1e-7
    if version == "python": np.random.seed(7)
    sc = _chain(version)
    H = s*_ham(sc)
    sc.set_hamiltonian(H)
    e0, _ = _quiet(sc.gs_energy)
    psil, psir = sc.nh_left_wf, sc.wf0
    q = psil.dot(H*psir)/psil.dot(psir)
    assert abs(q - e0) < 1e-8*abs(e0)


# ---------------------------------------------------------------- 26 --

def test_residual_scale_is_the_old_denominator_at_unit_scale():
    """c + |E - e_id| is exactly 1 + |E| for any H with no identity term
    and a largest coefficient of 1 or more, and the unit scale is exactly
    1 there, so the certificate and the solve are the old ones."""
    sc = _chain("python")
    terms = _ham(sc).to_terms()
    c, e_id = nhdmrg._residual_scale(terms)
    assert c == 1.0 and e_id == 0.0
    assert nhdmrg._unit_scale_up(nhdmrg._max_abs_coef(terms)) == 1.0
    c, e_id = nhdmrg._residual_scale((1e-5*_ham(sc) + 7.0).to_terms())
    assert c == pytest.approx(1e-5) and e_id == pytest.approx(7.0)
    c, sh, sa = nhdmrg._generalized_residual_scale(terms,
            (1 + 0.2*sc.Sz[0]).to_terms())
    assert (c, sh, sa) == (1.0, 0.0, 1.0)


def test_unit_scale_up_matches_mo_terms():
    """mo_terms.h's unit_scale_up(): the power of two bringing a largest
    coefficient below 1 into [1,2), exactly 1 at 1 or more."""
    for cmax in (1.0, 1.5, 7.0, 1e6):
        assert nhdmrg._unit_scale_up(cmax) == 1.0
    for cmax in (0.0, -1.0, float("nan")):
        assert nhdmrg._unit_scale_up(cmax) == 1.0
    for cmax in (0.5, 0.3, 1e-7, 1e-13, 3e-300):
        up = nhdmrg._unit_scale_up(cmax)
        assert 1.0 <= cmax*up < 2.0
        assert math.frexp(up)[0] == 0.5 # a power of two


@pytest.mark.parametrize("s,c", [(1.0, 0.0), (1e-5, 0.0), (1.0, 1e4)],
                         ids=["s1", "s1e-5", "offset1e4"])
@pytest.mark.parametrize("version", [_backend("python"), _backend(3), _backend(2)])
def test_certificate_flags_an_unconverged_run_at_any_units(version, s, c):
    """A deliberately unconverged run (maxm=2, nsweeps=1; relative residual
    2.4e-2 to 4.0e-2) must warn at s=1, in small units and next to an
    offset. Before: it warned at s=1 only, and at s=1e-5 and c=1e4 it was
    accepted on its first attempt (certificate 1.2e-6 and 7.9e-6)."""
    if version == "python": np.random.seed(7)
    sc = _chain(version, maxm=2, nsweeps=1)
    sc.set_hamiltonian(s*_ham(sc) + c)
    _, out = _quiet(nhdmrg.nhdmrg, sc)
    assert "Warning: nhdmrg did not reach" in out


@pytest.mark.parametrize("version", [_backend("python"), _backend(3)])
def test_generalized_certificate_flags_an_unconverged_run_in_small_units(version):
    """nhdmrg_generalized's copy of the certificate, A = 1 + 0.2*Sz0: the
    same unconverged schedule must warn at s=1e-5 (before: 1 attempt, no
    warning, certificate 2e-6 at a relative residual of 6e-2 to 8e-2)."""
    if version == "python": np.random.seed(7)
    sc = _chain(version, maxm=2, nsweeps=1)
    sc.set_hamiltonian(1e-5*_ham(sc))
    _, out = _quiet(nhdmrg.nhdmrg_generalized, sc, 1 + 0.2*sc.Sz[0])
    assert "Warning: nhdmrg_generalized did not reach" in out


@pytest.mark.parametrize("version", [_backend(3), _backend(2), _backend("python")])
def test_certificate_passes_a_converged_run_in_small_units(version):
    """...and the relative test raises no false alarm on a converged run
    at the same units (or next to the same offset). At s=1e-14 and 1e-16
    "python"'s own H*psir, left in the caller's units, is off by 2 and 23
    per cent, which read a pair whose energy is exact to 5e-15 as a
    relative residual of 0.16 and 0.54 and warned: the certificate's
    algebra runs at the solve's unit scale for that reason."""
    for s, c in ((1e-5, 0.0), (1e-14, 0.0), (1e-16, 0.0), (1.0, 1e4)):
        if version == "python": np.random.seed(7)
        sc = _chain(version)
        sc.set_hamiltonian(s*_ham(sc) + c)
        _, out = _quiet(nhdmrg.nhdmrg, sc, ntries=1)
        assert "Warning" not in out, (s, c, out)
    if version == 2: return # no generalized NH solver on v2
    for s in (1e-14, 1e-16):
        if version == "python": np.random.seed(7)
        sc = _chain(version)
        sc.set_hamiltonian(s*_ham(sc))
        _, out = _quiet(nhdmrg.nhdmrg_generalized, sc, 1 + 0.2*sc.Sz[0],
                        ntries=1)
        assert "Warning" not in out, (s, out)
