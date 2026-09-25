"""Regression tests for the `ednum` cluster of the 2026-09-25b hole hunt
(docs/audit_2026_09_25b_hole_hunt.md): findings 15, 20, 21 and 22, and the
lead scale-arnolditk-absolute-stop.

All five are absolute thresholds in the numerical kernels under the ED
reference and the MPS Krylov helpers, reached once e7b1196 let a
Hamiltonian written in small units keep its terms:

* 15: mode="ED" submode="ED" at T=0 averages with equal weight over every
  level below an absolute dex=1e-5, and its sensitivity warning looked only
  at [dex/3, 3*dex], so a whole spectrum narrower than dex/3 came back as
  the infinite-temperature spectrum with no warning. The warning now also
  fires when the averaged manifold is wider than the broadening; no number
  changes.
* 20: the ED eigensolver above 2000 states (_deflated_lowest_hermitian)
  carried four absolute constants; in small units it smuggled a magnon in
  as a ground-multiplet copy or stalled for hours. It now runs on the
  matrix scaled to unit infinity norm when that norm is below 1.
* 21: submode="ROOTN" stopped its Lanczos basis on an absolute beta, so
  below a seed spread of 1e-12 (ED) or 1e-10 (MPS) the spectrum was one
  Lorentzian. beta is compared against ||H q_0|| now.
* 22: the ED State.normalize() took no tol and returned None silently; the
  arnolditk/power-method consumers of H|psi> hit the 1e-8 norm floor.
* lead: arnolditk's residual stop, warm-up stop and invariant-subspace test
  were absolute energies.

The anchors are the same calculation at unit scale and linearity in the
scale s of the Hamiltonian, so no golden number enters.
"""

import io
import contextlib
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, spinchain
from dmrgpy.edtk.dynamics import check_dex_sensitivity


def _heisenberg(sc, n, J=1.0):
    h = 0
    for i in range(n-1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                   + sc.Sz[i]*sc.Sz[i+1])
    return h


# ---------------------------------------------------------------------
# 15: the dex manifold of submode="ED"
# ---------------------------------------------------------------------

ES1 = np.linspace(-0.5, 4.0, 91)
D1 = 0.2


def _ed_correlator(s, J=1.0, field=0.3, **kw):
    """s*C_s(s*w) of Sz0 on the 6-site chain s*(J*Heisenberg + field*Sz0),
    es and delta scaled with it, and the RuntimeWarnings it raised."""
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version="python")
    sc.set_hamiltonian(s*(_heisenberg(sc, 6, J) + field*sc.Sz[0]))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _, y = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                es=s*ES1, delta=s*D1, mode="ED", submode="ED", **kw)
    msgs = [str(m.message) for m in w if issubclass(m.category, RuntimeWarning)]
    return s*np.asarray(y), msgs


@pytest.mark.parametrize("s", [3e-7, 1e-7, 1e-13])
def test_dex_above_the_whole_spectrum_is_no_longer_silent(s):
    """Every level of s*H below dex/3: the old guard saw no level near the
    cutoff and said nothing while the call returned the infinite-
    temperature spectrum, 0.738 of the peak off. It still returns that
    (no number changes), but says so now."""
    ref, msgs = _ed_correlator(1.0)
    assert msgs == []
    y, msgs = _ed_correlator(s)
    assert np.max(np.abs(y-ref))/np.max(np.abs(ref)) > 0.5  # still the dex average
    assert len(msgs) == 1 and "wide" in msgs[0] and "nex=64" in msgs[0]


@pytest.mark.parametrize("s", [1e-7, 1e-13])
def test_the_documented_remedy_is_exact_and_quiet(s):
    """dex chosen in the Hamiltonian's units: scale covariant, no warning."""
    ref, _ = _ed_correlator(1.0)
    y, msgs = _ed_correlator(s, dex=s*1e-5)
    assert np.max(np.abs(y-ref))/np.max(np.abs(ref)) < 1e-10
    assert msgs == []


def test_a_genuine_multiplet_is_averaged_quietly():
    """The 6-site ferromagnet: the dex average over its 7-fold S=3 ground
    multiplet is the intended use, and its width is roundoff."""
    _, msgs = _ed_correlator(1.0, J=-1.0, field=0.0)
    assert msgs == []


def test_width_guard_unit():
    quiet = np.array([0.0, 1e-15, 2e-15, 0.4, 1.0])
    loud = np.array([0.0, 1e-7, 2e-7, 3e-7])  # every level below dex/3
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        check_dex_sensitivity(quiet, 1e-5, delta=1e-3)
        check_dex_sensitivity(loud, 1e-5)  # no delta: the old test alone
        check_dex_sensitivity(loud, 1e-5, delta=1e-6)  # unresolved: fine
    with pytest.warns(RuntimeWarning, match="wide"):
        check_dex_sensitivity(loud, 1e-5, delta=1e-8)


# ---------------------------------------------------------------------
# 21: the Lanczos breakdown of submode="ROOTN"
# ---------------------------------------------------------------------

ESR = np.linspace(-0.5, 4.0, 31)


def _old_lanczos_basis(H, v0, k):
    """algebra/rootn.py's lanczos_basis before the fix, verbatim apart from
    its docstring: the absolute beta<1e-12 break."""
    n = v0.shape[0]
    q = v0/np.linalg.norm(v0)
    Q = [q]
    alphas = []
    betas = []
    beta = 0.0
    q_prev = np.zeros(n, dtype=complex)
    for j in range(k):
        w = H@Q[j] - beta*q_prev
        alpha = np.vdot(Q[j], w).real
        alphas.append(alpha)
        w = w - alpha*Q[j]
        for q_i in Q:
            w = w - np.vdot(q_i, w)*q_i
        beta = np.linalg.norm(w)
        if j == k-1: break
        if beta < 1e-12: break
        betas.append(beta)
        q_prev = Q[j]
        Q.append(w/beta)
    nb = len(Q)
    T = np.diag(alphas)
    for j in range(nb-1):
        T[j, j+1] = betas[j]
        T[j+1, j] = betas[j]
    Q = np.array(Q).T
    return Q, T


def _rootn(s, version="python", mode="ED", es=ESR, **kw):
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=version)
    sc.maxm = 40
    sc.nsweeps = 12
    sc.set_hamiltonian(s*(_heisenberg(sc, 6) + 0.3*sc.Sz[0]))
    _, y = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=s*es,
            delta=s*D1, mode=mode, submode="ROOTN", **kw)
    return s*np.asarray(y)


def _basis_sizes(monkeypatch):
    from dmrgpy.algebra import rootn
    sizes = []
    inner = rootn.lanczos_basis
    def wrapped(H, v0, k):
        Q, T = inner(H, v0, k)
        sizes.append(Q.shape[1])
        return Q, T
    monkeypatch.setattr(rootn, "lanczos_basis", wrapped)
    return sizes


@pytest.mark.parametrize("s", [1e-12, 1e-13, 1e-20])
def test_rootn_on_ed_is_scale_covariant(s, monkeypatch):
    """Below s=2e-12 the old break left the seed alone (basis 1 of 20) and
    the spectrum was a single pole, 1.85 of the peak off here."""
    ref = _rootn(1.0)
    sizes = _basis_sizes(monkeypatch)
    y = _rootn(s)
    assert min(sizes) == 20
    assert np.max(np.abs(y-ref))/np.max(np.abs(ref)) < 1e-5


def test_rootn_on_ed_is_bit_for_bit_the_old_one_at_unit_scale(monkeypatch):
    from dmrgpy.algebra import rootn
    new = _rootn(1.0)
    monkeypatch.setattr(rootn, "lanczos_basis", _old_lanczos_basis)
    old = _rootn(1.0)
    assert np.array_equal(new, old)


def test_breakdown_rule_is_relative_and_safe_at_zero():
    from dmrgpy.algebra.rootn import is_breakdown
    assert is_breakdown(0.0, 0.0, 1e-12)  # H q_0 = 0: an exact zero breaks
    # beta=1e-17 in units where ||H q_0||=1e-7: an ordinary direction (the
    # old absolute 1e-12 broke here), and a breakdown at unit scale
    assert not is_breakdown(1e-17, 1e-7, 1e-12)
    assert is_breakdown(1e-17, 1.0, 1e-12)


def _old_lanczos_basis_mps(self, Hmpo, v, nkry):
    """rootndmrg.py's _lanczos_basis_mps before the fix, verbatim apart from
    its docstring: the absolute beta<1e-10 break."""
    nrm = np.sqrt(v.dot(v).real)
    q = (1./nrm)*v
    Q = [q]
    alphas = []
    betas = []
    qprev = None
    beta = 0.0
    for it in range(nkry):
        w = Hmpo*Q[-1]
        if qprev is not None: w = w - beta*qprev
        alpha = Q[-1].dot(w).real
        alphas.append(alpha)
        w = w - alpha*Q[-1]
        for qi in Q: w = w - qi.dot(w)*qi
        if it == nkry-1: break
        beta = np.sqrt(abs(w.dot(w).real))
        if beta < 1e-10: break
        betas.append(beta)
        qprev = Q[-1]
        Q.append((1./beta)*w)
    return Q, np.array(alphas), np.array(betas)


V3 = pytest.mark.skipif(not cppext.available(3),
                        reason="needs the compiled v3 extension")


# the MPS route costs a truncated MPO application per Lanczos step, so it is
# checked on a short grid and a short recursion (the collapse is at the
# first step whatever N and nkry are)
ESV3 = np.linspace(-0.5, 4.0, 7)
V3KW = dict(version=3, mode="DMRG", es=ESV3, N=2, nkry=8)


@V3
def test_rootn_on_the_mps_route_is_scale_covariant():
    """The DMRG twin (rootndmrg.py) broke on an absolute 1e-10, a hundred
    times earlier than ED: 1.85 of the peak off at s=1e-10 on v3."""
    ref = _rootn(1.0, **V3KW)
    y = _rootn(1e-10, **V3KW)
    assert np.max(np.abs(y-ref))/np.max(np.abs(ref)) < 1e-8


@V3
def test_rootn_on_the_mps_route_is_the_old_one_at_unit_scale(monkeypatch):
    """Same chain, same ground state, both break rules: at unit scale no
    basis breaks early under either, so the curves are the same numbers."""
    from dmrgpy import rootndmrg
    sc = spinchain.Spin_Chain(["S=1/2"]*6, itensor_version=3)
    sc.maxm = 40
    sc.nsweeps = 12
    sc.set_hamiltonian(_heisenberg(sc, 6) + 0.3*sc.Sz[0])
    kw = dict(name=(sc.Sz[0], sc.Sz[0]), es=ESV3, delta=D1, mode="DMRG",
              submode="ROOTN", N=2, nkry=8)
    _, new = sc.get_dynamical_correlator(**kw)
    monkeypatch.setattr(rootndmrg, "_lanczos_basis_mps", _old_lanczos_basis_mps)
    _, old = sc.get_dynamical_correlator(**kw)
    assert np.array_equal(np.asarray(new), np.asarray(old))


# ---------------------------------------------------------------------
# 20: the ED eigensolver above 2000 states, in small units
# ---------------------------------------------------------------------

@pytest.fixture
def capped_arpack(monkeypatch):
    """Cap ARPACK at 1000 restarts, so a stalled deflated round raises
    ArpackNoConvergence in seconds instead of running for hours at the
    library's maxiter=1e6. Converging calls need far fewer."""
    from dmrgpy.algebra import algebra
    monkeypatch.setattr(algebra, "maxiter", 1000)
    return algebra


def _ferromagnet_levels(s, n=6, L=12):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
    sc.set_hamiltonian(s*_heisenberg(sc, L, J=-1.0))
    return np.real(sc.get_excited(n=n, mode="ED"))/s


@pytest.mark.parametrize("s", [2e-7, 1e-9])
def test_ground_multiplet_in_small_units_keeps_every_copy(s, capped_arpack):
    """12-site ferromagnet (dim 4096, the deflated ARPACK path), 13-fold
    ground multiplet at -2.75*s. The absolute partner window 1e-8 took the
    magnon, 0.034*s above, as a copy below s=2.9e-7: 8 of 8 calls wrong at
    s=2e-7 on the unfixed tree."""
    for _ in range(4):  # the failure went through ARPACK's start vector
        es = _ferromagnet_levels(s)
        assert len(es) == 6
        assert np.max(np.abs(es + 2.75)) < 1e-8


@pytest.mark.parametrize("s", [1e-5, 1e-6])
def test_deflated_rounds_converge_in_small_units(s, capped_arpack):
    """The parked levels sat near an absolute +10 over a spectrum of order
    s, and ARPACK could not converge next to them: 4 of 8 calls raised
    ArpackNoConvergence at each of these scales on the unfixed tree (with
    the same cap; uncapped, a hang of hours)."""
    for _ in range(6):
        es = _ferromagnet_levels(s)
        assert np.max(np.abs(es + 2.75)) < 1e-8


def _degenerate_matrix(dim=600, g=6, seed=3):
    rng = np.random.default_rng(seed)
    d = np.concatenate([-np.ones(g), np.linspace(-0.97, 1.0, dim-g)])
    Q, _ = np.linalg.qr(rng.normal(size=(dim, dim)))
    h = (Q*d)@Q.T
    return 0.5*(h+h.T)


@pytest.mark.parametrize("s", [1.0, 1e-7, 1e-9, 1e-13])
def test_deflated_solver_is_scale_covariant(s, capped_arpack):
    """An exactly 6-fold lowest level at -1 with the next at -0.97, straight
    into the solver: the magnon-for-copy swap fired in 3 of 4 calls at
    s=1e-7 on the unfixed tree."""
    h = _degenerate_matrix()
    for _ in range(3):
        es, vs = capped_arpack._deflated_lowest_hermitian(s*h, 6)
        assert len(es) == 6
        assert np.max(np.abs(es/s + 1.0)) < 1e-10


def test_unit_scale_runs_unscaled_and_small_units_at_unit_norm(monkeypatch):
    """At an infinity norm of 1 or more the routine is entered once and
    runs as it always did; below, it re-enters itself exactly once, on the
    matrix scaled by a power of two into [1,2)."""
    from dmrgpy.algebra import algebra
    seen = []
    inner = algebra._deflated_lowest_hermitian
    def wrapped(h, n):
        seen.append(algebra._infinity_norm(h))
        return inner(h, n)
    monkeypatch.setattr(algebra, "_deflated_lowest_hermitian", wrapped)
    h = _degenerate_matrix(dim=300)
    assert algebra._infinity_norm(h) >= 1.0
    algebra._deflated_lowest_hermitian(h, 2)
    assert len(seen) == 1
    seen.clear()
    s = 3e-9
    es, _ = algebra._deflated_lowest_hermitian(s*h, 2)
    assert len(seen) == 2 and 1.0 <= seen[1] < 2.0
    assert seen[1]/seen[0] == 2.0**round(np.log2(seen[1]/seen[0]))  # exact power of two
    assert np.max(np.abs(es/s + 1.0)) < 1e-10


# ---------------------------------------------------------------------
# 22: normalize() on a state written in small units
# ---------------------------------------------------------------------

def _small_state(version, mode):
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=version)
    sc.maxm = 16
    sc.nsweeps = 10
    if mode == "ED":
        sc.mode = "ED"
    sc.set_hamiltonian(_heisenberg(sc, 4))
    with contextlib.redirect_stdout(io.StringIO()):
        gs = sc.get_gs()
    return sc, gs


def test_ed_normalize_takes_the_documented_tol(capsys):
    """wf.normalize(tol=...) raised TypeError on mode="ED", and below the
    floor it returned None without the warning the MPS backends print."""
    sc, gs = _small_state("python", "ED")
    ref = (sc.Sx[0]*gs).normalize()
    x = (1e-9*sc.Sx[0])*gs
    capsys.readouterr()
    assert x.normalize() is None
    assert "not normalizable" in capsys.readouterr().out
    y = x.normalize(tol=0.)
    assert abs(y.dot(y) - 1.0) < 1e-12
    assert abs(abs(ref.dot(y)) - 1.0) < 1e-12
    assert abs(x.norm() - 5e-10) < 1e-15


def test_ed_normalize_does_not_turn_a_cancelled_state_into_noise():
    sc, gs = _small_state("python", "ED")
    with contextlib.redirect_stdout(io.StringIO()):
        assert (gs - gs).normalize(tol=0.) is None


@pytest.mark.parametrize("version,mode", [("python", "ED"), ("python", "DMRG")])
def test_normalize_image_divides_by_the_norm_itself(version, mode):
    from dmrgpy.algebra.krylov import normalize_image
    sc, gs = _small_state(version, mode)
    ref = (sc.Sx[0]*gs).normalize()
    y = normalize_image((1e-12*sc.Sx[0])*gs)
    assert abs(y.dot(y) - 1.0) < 1e-10
    assert abs(abs(ref.dot(y)) - 1.0) < 1e-10


# ---------------------------------------------------------------------
# lead: arnolditk's stops (and 22's consumer powermethod.estimate_radius)
# ---------------------------------------------------------------------

E0_ED_4 = None


def _ed_e0_4():
    global E0_ED_4
    if E0_ED_4 is None:
        sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
        sc.set_hamiltonian(_heisenberg(sc, 4))
        E0_ED_4 = sc.gs_energy(mode="ED")
    return E0_ED_4


def _arnoldi_e0(s, seed, **kw):
    from dmrgpy import mpsalgebra
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    sc.maxm = 16
    sc.nsweeps = 10
    h = s*_heisenberg(sc, 4)
    sc.set_hamiltonian(h)
    with contextlib.redirect_stdout(io.StringIO()):
        out = mpsalgebra.lowest_energy_arnoldi(sc, h, **kw)
    return np.real(np.array(out[0]).ravel()[0])/s


@pytest.mark.parametrize("seed", [11, 12])
@pytest.mark.parametrize("s", [1e-6, 1e-9, 1e-11])
def test_arnoldi_ground_state_is_scale_covariant(s, seed):
    """E0/s was 6.2e-4 and 4.0e-3 relative off at s=1e-6 (the residual stop
    was met in absolute units after one outer iteration), and at s=1e-9
    raised TypeError (estimate_radius lost its iterate to normalize()'s
    floor). The same seed at s=1 is within 1.5e-8 and 1.0e-7. Below
    s=1e-11 "python"'s own MPO product (s*H)|u> drifts (3e-5 relative at
    1e-12), which is not arnolditk's; see the ED-vector test below."""
    e1 = _arnoldi_e0(1.0, seed)
    es = _arnoldi_e0(s, seed)
    ed = _ed_e0_4()
    assert abs(e1-ed)/abs(ed) < 1e-6
    assert abs(es-ed)/abs(ed) < 1e-6


@pytest.fixture
def ed_vectors(monkeypatch):
    """Run arnolditk on exact ED vectors (State objects) rather than MPS,
    so that nothing but arnolditk's own tests can depend on the units."""
    from dmrgpy.algebra import arnolditk
    monkeypatch.setattr(arnolditk, "arnoldimode", "ED")


@pytest.mark.parametrize("s", [1e-9, 1e-13])
def test_arnoldi_on_exact_vectors_is_scale_covariant(s, ed_vectors):
    """With ED vectors the whole route is exactly linear in s: same seed,
    same numbers, down to 1e-13 (E0/s 4.3e-8 relative off ED at every s)."""
    e1 = _arnoldi_e0(1.0, 11)
    es = _arnoldi_e0(s, 11)
    assert es == pytest.approx(e1, rel=1e-10)
    assert abs(e1-_ed_e0_4())/abs(_ed_e0_4()) < 1e-6


def test_arnoldi_invariant_subspace_test_is_relative():
    """delta=s*1e-3 at s=1e-7, seed 11: two Krylov directions fell under
    the absolute beta<1e-8, were replaced by random vectors coupled with
    weight one, and the call returned E0/s=+0.351 for -1.616."""
    s = 1e-7
    e = _arnoldi_e0(s, 11, delta=s*1e-3)
    ed = _ed_e0_4()
    assert abs(e-ed)/abs(ed) < 1e-6


def test_energy_unit_is_one_at_ordinary_couplings():
    """At J=1, and at any largest coefficient above 1, every arnolditk stop
    is the old absolute one; below, it scales with the Hamiltonian."""
    from dmrgpy.algebra.arnolditk import energy_unit
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    h = _heisenberg(sc, 4) + 0.3*sc.Sz[0]
    assert energy_unit(h) == 1.0
    assert energy_unit(10.*h) == 1.0
    assert energy_unit(2.5e-9*h) == 2.5e-9


def test_estimate_radius_is_scale_covariant(ed_vectors):
    """It raised TypeError below s of about 1e-8 (normalize()'s floor)."""
    from dmrgpy.algebra import powermethod
    out = []
    for s in (1.0, 1e-12):
        np.random.seed(5)
        sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
        h = s*_heisenberg(sc, 4)
        out.append(powermethod.estimate_radius(sc, h)/s)
    assert out[1] == pytest.approx(out[0], rel=1e-12)


def test_krylov_hermiticity_test_is_relative_below_unit_scale():
    """The Krylov-matrix Hermiticity test was an absolute 1e-6, which every
    matrix in small units passed: a non-Hermitian one then went to eigh,
    which reads one triangle only (eigenvalues +-0.25 here for the true
    +-0.5)."""
    from dmrgpy.algebra import krylov
    herm = np.array([[1.0, 0.2j], [-0.2j, -0.5]])
    nonh = np.array([[0.0, 1.0], [0.25, 0.0]])
    for s in (1.0, 1e-9):
        assert krylov.is_hermitian_matrix(s*herm)
        assert not krylov.is_hermitian_matrix(s*nonh)
        es, _ = krylov.diagonalize(s*nonh)
        assert sorted(np.real(es)/s) == pytest.approx([-0.5, 0.5])
    assert krylov.is_hermitian_matrix(np.zeros((2, 2)))
    # at and above unit scale it is the old absolute test, bit for bit
    assert krylov.is_hermitian_matrix(np.array([[2.0, 0.5e-6], [0.0, 1.0]]))
