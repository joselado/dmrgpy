"""Regression tests for the `cvm` cluster of the 2026-09-25b hole hunt
(docs/audit_2026_09_25b_hole_hunt.md, findings 13, 14 and 16).

13. submode="CVM" stopped its linear solve on an absolute residual,
    cvm_tol=1e-5, while the right-hand side b = -eta*B|GS> carries the
    units of eta and of B. Wherever ||b|| was at or below 1e-5 -- a small
    operator, a Hamiltonian in small units, or eta <= 2e-5 at unit scale,
    get_kondo_spectrum's documented default delta=2e-6 included -- the
    start passed at iteration 0 and every frequency returned the flat
    eta*<AB>/pi. applyinverse_dmrg handed the same absolute cvm_tol to the
    session BiCGSTAB that submode="CVM_explicit" inverts with. Both are
    relative to the right-hand side now.
14. The early exits of that solve (cvm_patience, cvm_blowup) read the
    residual 2-norm, which exact conjugate gradient does not decrease
    monotonically, and returned the best iterate by that norm: at full
    bond dimension, on a line, the initial guess (1.5915e-04 against
    13.2636). They read the CG functional phi now, which exact CG lowers
    at every iteration, and they still stop a truncated solve early.
16. The A^dagger == B gate of CVM_explicit, of the non-Hermitian CVM/INV
    and of cvm_solver="variational" was an absolute test on a squared norm
    (Many_Body_Chain.is_zero_operator's 1e-4, EDchain's 1e-8 on
    Tr(D D^dagger)), so a small non-adjoint pair was admitted and answered
    with another pair's density. CVM_explicit and the non-Hermitian route
    now evaluate <GS|A (z-H)^-1 B|GS> for any pair, with no gate; the
    variational gate is canonical.is_dagger_pair first, and both
    is_zero_operator helpers are relative.

The anchors are mode="ED" (exact inverses, dense Lehmann sums) on the same
chain, and the same calculation at unit scale, so no golden number enters.
"""

import contextlib
import io
import warnings

import numpy as np
import pytest

from dmrgpy import cppext, cvm, spinchain
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc


def _backend(version):
    marks = ()
    if version in (2, 3):
        marks = pytest.mark.skipif(not cppext.available(version),
            reason="needs the compiled itensor_version=%d extension" % version)
    ident = "v%d" % version if version in (2, 3) else str(version)
    return pytest.param(version, id=ident, marks=marks)


BACKENDS = [_backend("python"), _backend(3), _backend(2)]


def _chain(version, L=6, s=1.0, maxm=40):
    """Open S=1/2 Heisenberg chain with 0.3*Sz0, times s; maxm above the
    exact bond dimension, so the solver is the only error left."""
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=version)
    sc.maxm, sc.nsweeps, sc.cvm_maxm = maxm, 12, maxm
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(s*(h + 0.3*sc.Sz[0]))
    return sc


def _corr(sc, A, B, es, delta, mode, submode):
    with contextlib.redirect_stdout(io.StringIO()):
        _x, y = sc.get_dynamical_correlator(name=(A, B), es=es, delta=delta,
                                            mode=mode, submode=submode)
    return np.asarray(y)


def _rel(y, ref):
    return np.max(np.abs(np.asarray(y) - ref))/np.max(np.abs(ref))


ES = np.linspace(0.2, 3.0, 5)
DELTA = 0.2


# ------------------------------------------------------------ finding 13

@pytest.mark.parametrize("version", BACKENDS)
def test_cvm_is_covariant_under_the_operator_scale(version):
    """C[eps*A, eps*B] = eps^2 C[A,B]: at eps=1e-4 every frequency was the
    flat eta*<Sz0 Sz0>/pi = 0.0159, 0.887 of the peak off."""
    sc = _chain(version)
    ref = _corr(sc, sc.Sz[0], sc.Sz[0], ES, DELTA, "ED", "INV")
    for eps in (1.0, 1e-4):
        sc = _chain(version)
        y = _corr(sc, eps*sc.Sz[0], eps*sc.Sz[0], ES, DELTA, "DMRG", "CVM")
        assert _rel(y/eps**2, ref) < 1e-4, eps


@pytest.mark.parametrize("version", BACKENDS)
def test_cvm_is_covariant_under_the_hamiltonian_units(version):
    """s*C_s(s*w; s*eta) = C_1(w) for H -> s*H: at s=1e-4 every frequency
    was 1.5915e-10, the whole peak off."""
    sc = _chain(version)
    ref = _corr(sc, sc.Sz[0], sc.Sz[0], ES, DELTA, "ED", "INV")
    s = 1e-4
    sc = _chain(version, s=s)
    y = s*_corr(sc, sc.Sz[0], sc.Sz[0], s*ES, s*DELTA, "DMRG", "CVM")
    assert _rel(y, ref) < 1e-4


@pytest.mark.parametrize("version", BACKENDS)
def test_cvm_explicit_is_covariant_under_the_operator_scale(version):
    """The same absolute tolerance reached submode="CVM_explicit" through
    applyinverse_dmrg: 0.12 to 0.37 of the peak off at eps=1e-4."""
    sc = _chain(version)
    ref = _corr(sc, sc.Sz[0], sc.Sz[0], ES, DELTA, "ED", "INV")
    eps = 1e-4
    y = _corr(sc, eps*sc.Sz[0], eps*sc.Sz[0], ES, DELTA, "DMRG",
              "CVM_explicit")
    assert _rel(y/eps**2, ref) < 1e-3


@pytest.mark.parametrize("version", BACKENDS)
def test_applyinverse_is_linear_in_the_right_hand_side(version):
    """A^-1 (s*wf) = s*A^-1 wf: the session BiCGSTAB compared an absolute
    cvm_tol with its residual, so a small wf passed after one step. (wf
    must not be an eigenvector of A, which one step inverts exactly.)"""
    sc = _chain(version)
    wf = sc.Sz[0]*sc.get_gs()
    M = -sc.hamiltonian + (sc.gs_energy() + 0.7 + 0.3j)
    x1 = sc.applyinverse(M, wf)
    s = 1e-6
    xs = sc.applyinverse(M, s*wf)
    d = x1 - (1./s)*xs
    assert np.sqrt(abs(d.dot(d))) < 1e-3*np.sqrt(abs(x1.dot(x1)))


def _kondo_chain(version, J=1e-3):
    """Three S=1/2 in eV units, J=1 meV and a 0.4 meV Zeeman term"""
    sc = spinchain.Spin_Chain(["1/2"]*3, itensor_version=version)
    sc.maxm, sc.nsweeps = 20, 10
    sc.set_hamiltonian(J*(sc.SS(0, 1) + sc.SS(1, 2))
                       + 0.4*J*(sc.Sz[0] + sc.Sz[1] + sc.Sz[2]))
    return sc


@pytest.mark.parametrize("version", [_backend("python"), _backend(3)])
def test_kondo_second_order_cvm_at_the_default_broadening(version):
    """get_kondo_spectrum(mode="DMRG", submode="CVM", order=2, T=0) at its
    documented delta=2e-6 returned 9.12e-09 against an exact 4.7124. The
    anchor is the exact Lehmann correlator on the SAME es grid, which
    shares the grid's discretization, so the grid can be coarse."""
    delta = 2e-6
    sc = _kondo_chain("python")
    with contextlib.redirect_stdout(io.StringIO()):
        E = np.sort(np.real(np.asarray(sc.get_excited(mode="ED", n=8))))
    lines = np.unique(np.round(E - E[0], 12))
    eVs = 1e-3*np.linspace(-3, 3, 7)
    es = [np.linspace(-10*delta, 3.2e-3, 40)]
    es += [l + delta*np.linspace(-4, 4, 9) for l in lines if l < 3.2e-3]
    es = np.sort(np.concatenate(es))
    with contextlib.redirect_stdout(io.StringIO()):
        ref = second_order_dIdV_dc(sc, 0, eVs, mode="ED", submode="ED",
                                   delta=delta, es=es)
        _, y = _kondo_chain(version).get_kondo_spectrum(
            eVs, site=0, T=0.0, order=2, mode="DMRG", submode="CVM",
            delta=delta, es=es)
    assert _rel(y, ref) < 1e-4


# ------------------------------------------------------------ finding 14

@pytest.mark.parametrize("version", BACKENDS)
def test_cvm_resolves_a_line_at_full_bond_dimension(version):
    """On a line at eta=2e-3 and 2e-4 the exact CG residual rose ~200x
    before it fell, cvm_blowup fired and the initial guess came back:
    1.5915e-04 against 13.2636, 1.5915e-05 against 132.633."""
    sc = _chain("python")
    with contextlib.redirect_stdout(io.StringIO()):
        E = np.sort(np.real(np.asarray(sc.get_excited(mode="ED", n=64))))
    w0 = E[2] - E[0]
    for delta in (2e-3, 2e-4):
        es = w0 + delta*np.array([-2., 0., 1.])
        ref = _corr(sc, sc.Sz[0], sc.Sz[0], es, delta, "ED", "INV")
        y = _corr(_chain(version), sc.Sz[0], sc.Sz[0], es, delta,
                  "DMRG", "CVM")
        assert _rel(y, ref) < 1e-5, delta


@pytest.mark.parametrize("version", [_backend(3)])
def test_cvm_on_the_test_suites_own_broadening(version):
    """10 sites at eta=0.05: cvm_patience fired at 7 of 121 points of
    linspace(0,3,121), up to 0.989 of the peak off. This is the worst of
    them, 3.97887e-03 against an exact 0.374193, the grid's peak."""
    es = np.array([0.75])
    sc = _chain("python", L=10)
    ref = _corr(sc, sc.Sz[0], sc.Sz[0], es, 0.05, "ED", "INV")
    y = _corr(_chain(version, L=10), sc.Sz[0], sc.Sz[0], es, 0.05,
              "DMRG", "CVM")
    assert _rel(y, ref) < 1e-4


def test_a_truncated_solve_still_stops_early_and_says_so():
    """What the exits were added for (a67228e): once cvm_maxm truncates
    the correction vector, the recurrence stops improving, and the solve
    has to stop well short of cvm_nit and warn rather than run the whole
    budget. A 12-site chain at cvm_maxm=4 is far below the bond dimension
    the correction vector needs. (The old loop stopped early here too,
    but its absolute warning test stayed silent.)"""
    sc = spinchain.Spin_Chain(["S=1/2"]*12, itensor_version=3) \
        if cppext.available(3) else None
    if sc is None: pytest.skip("needs the compiled itensor_version=3 extension")
    h = 0
    for i in range(11):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cvm_maxm = 40, 12, 4
    sc.get_gs()
    cvm._UNCONVERGED_WARNED.clear()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        C, _xc, nit, res = cvm.cvm_correction_vector(
            sc, sc.Sz[0], sc.Sz[0], 0.3, 0.15, tol=1e-5,
            max_it=int(sc.cvm_nit))
    assert nit < 300
    assert any("did not converge" in str(x.message) for x in w)


def test_unconverged_solve_warns_relative_to_the_right_hand_side():
    """The warning compared the absolute residual with 100*cvm_tol, so a
    residual never reduced below its start was silent when ||b|| was small.
    Here max_it=1 at eta=2e-3 on a line: one step cannot converge, and at
    ||b|| ~ 8e-4 the old test (best_res > 1e-3) could never fire."""
    sc = _chain("python")
    with contextlib.redirect_stdout(io.StringIO()):
        E = np.sort(np.real(np.asarray(sc.get_excited(mode="ED", n=64))))
    sc.get_gs()
    cvm._UNCONVERGED_WARNED.clear()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        C, _xc, nit, res = cvm.cvm_correction_vector(
            sc, sc.Sz[0], sc.Sz[0], E[2]-E[0], 2e-3, tol=1e-5, max_it=1)
    assert nit == 1
    assert any("did not converge" in str(x.message) for x in w)


# ------------------------------------------------------------ finding 16

@pytest.mark.parametrize("version", BACKENDS)
def test_cvm_explicit_returns_the_asked_pair_at_any_scale(version):
    """(eps*Sx0, eps*Sy1) has purely imaginary Lehmann weights, so the
    adjoint pair's density is minus this one's. CVM_explicit raised at
    eps=1 and returned the adjoint pair's curve at eps=1e-2 (2.004 of
    the peak off); it now evaluates the asked pair at both."""
    es = ES[:3]
    for eps in (1.0, 1e-2):
        sc = _chain(version)
        A, B = eps*sc.Sx[0], eps*sc.Sy[1]
        ref = _corr(sc, A, B, es, DELTA, "ED", "INV")
        y = _corr(sc, A, B, es, DELTA, "DMRG", "CVM_explicit")
        assert _rel(y, ref) < 1e-3, eps


def test_variational_gate_sends_a_small_non_adjoint_pair_to_cg():
    """cvm_solver="variational" reads only B, so (1e-2*Sz0, 1e-2*Sz3)
    came back as C[Sz3,Sz3], 2.124 of the peak off."""
    sc = _chain("python")
    sc.cvm_solver = "variational"
    eps = 1e-2
    A, B = eps*sc.Sz[0], eps*sc.Sz[3]
    assert not cvm._use_ddmrg(sc, A, B)
    assert cvm._use_ddmrg(sc, A, eps*sc.Sz[0]) # proven adjoint: variational
    ref = _corr(sc, A, B, ES, DELTA, "ED", "INV")
    y = _corr(sc, A, B, ES, DELTA, "DMRG", "CVM")
    assert _rel(y, ref) < 1e-4


def test_variational_gate_falls_back_to_a_relative_numerical_test():
    """4*Sx^3 = Sx on S=1/2 is an identity the canonical form cannot see,
    so the numerical test decides it, the same way at any scale."""
    sc = _chain("python")
    sc.cvm_solver = "variational"
    for eps in (1.0, 1e-6):
        A = eps*sc.Sx[0]
        B = 4*eps*sc.Sx[0]*sc.Sx[0]*sc.Sx[0]
        assert cvm._use_ddmrg(sc, A, B), eps
        assert not cvm._use_ddmrg(sc, A, eps*sc.Sy[0]), eps


@pytest.mark.parametrize("version", ["ED", "python", 3])
def test_is_zero_operator_is_scale_free(version):
    """Zero is a property of the operator, not of its units: the absolute
    thresholds called any operator below about 1e-2 (DMRG) or
    sqrt(2e-8/dim) (ED) zero."""
    sc = _chain("python" if version == "ED" else version, L=4)
    iz = sc.get_ED_obj().is_zero_operator if version == "ED" \
        else sc.is_zero_operator
    for eps in (1.0, 1e-2, 1e-6):
        assert not iz(eps*(sc.Sx[0] - sc.Sy[1])), eps
        assert not iz(eps*sc.Sz[2]), eps
        assert iz(eps*(sc.Sx[0]*sc.Sx[0] - 0.25)), eps  # Sx^2 = 1/4
        assert iz(eps*(4*sc.Sx[1]*sc.Sx[1]*sc.Sx[1] - sc.Sx[1])), eps
        assert iz(eps*sc.Sz[0] - eps*sc.Sz[0]), eps
    assert iz(0.*sc.Sz[0])


def _nh_chain(L):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h + 0.3*sc.Sz[0] + 0.1j*sc.Sz[1])
    return sc


@pytest.mark.parametrize("submode", ["CVM", "INV"])
def test_ed_non_hermitian_resolvent_returns_the_asked_pair(submode):
    """mode="ED" on a non-Hermitian chain raised for (eps*Sx0, eps*Sy1) at
    eps >= 3e-5 and returned the adjoint pair's density at 1e-5, 2.064 of
    the peak off. The anchor is the routine's own formula, <v|A (z-H)^-1
    B|v> with v the right ground state, evaluated densely."""
    sc = _nh_chain(6)
    ed = sc.get_ED_obj()
    dense = lambda m: np.asarray(m.todense()) if hasattr(m, "todense") \
        else np.asarray(m)
    H = dense(ed.MO2matrix(sc.hamiltonian))
    with contextlib.redirect_stdout(io.StringIO()):
        v = ed.get_gs().v
        e0 = ed.gs_energy()
    for eps in (1.0, 1e-5):
        A, B = eps*sc.Sx[0], eps*sc.Sy[1]
        av = dense(ed.MO2matrix(A.get_dagger())) @ v
        bv = dense(ed.MO2matrix(B)) @ v
        ref = []
        for e in ES:
            g = [np.conj(av) @ np.linalg.solve(
                    (e0+e+1j*d)*np.eye(H.shape[0]) - H, bv)
                 for d in (DELTA, -DELTA)]
            ref.append(0.5j*(g[0]-g[1])/np.pi)
        y = _corr(sc, A, B, ES, DELTA, "ED", submode)
        assert _rel(y, np.array(ref)) < 1e-10, eps
