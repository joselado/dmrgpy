"""Fermionic (Jordan-Wigner-strung) dynamical correlators on an infinite
chain: the IBC real-time window (`pyitensor/idmrg_window.py` and its
ITensor v3 port `Chain::td_dynamical_correlator_window`) and the
finite-window KPM reduction (`Infinite_Many_Body_Chain.kpm_finite`).

Why this module exists. When fermionic infinite chains were made to work
(`docs/idmrg_fermionic_infinite_chain_plan.md`), the Jordan-Wigner string
was threaded through the Hamiltonian (`idmrg._term_site_matrices` ->
`_classify_terms`) and, in a follow-up, through the *static* two-point
correlator. The dynamical correlators were never covered, and the window
path turned out to be the same failure class one layer further on: both
`apply_local_operator` (ket) and `snapshot_correlator` (bra) inserted a
bare `C`/`Cdag` matrix with no string at all, on both backends. Measured
on the model below at t=0, `<Cdag_x C_0>` came out +0.203 at x=2 against
an exact -0.001, and +0.101 with the wrong sign at x=3 -- a different
number entirely, not a less converged one. Separately, the bra applied
`A` where the conjugation in `_close_array_chain` needs `A^dag`, so the
documented `<psi|A_x ... B_0|psi>` was silently `<psi|A^dag_x ... B_0|psi>`
-- invisible for every Hermitian name (only "Sz" appears in the existing
tests and examples), exactly wrong for C/Cdag/Sp/Sm.

Oracles, in increasing order of what they can catch:

1. `S(x, t=0)` must equal the *static* `two_point_correlator` for every x,
   to machine precision. At t=0 the ket's and bra's semi-infinite strings
   cancel identically down to the finite between-the-endpoints string the
   static path threads, so this is an exact identity, not a convergence
   comparison -- and it pins every sign, ordering, sublattice and
   window-pad detail at once. Note the anticommutation sign: for two
   parity-odd operators at different sites, `<A_x B_0> = -<B_0 A_x>`, so
   the x>0 reference carries an explicit minus.
2. The same `S(x, 0)` against the exact free-fermion one-body density
   matrix -- an external, backend-independent number.
3. `S(x, t)` against the exact free-fermion Green function
   `<c^dag_x(t) c_0> = sum_l [e^{iht}]_{xl} P[l,0]`, which is what
   actually exercises the string *through* the time evolution.

The model is the gapped, dimerized (SSH-like) spinless chain already used
by `tests/test_infinite_chain.py`'s own fermionic tests -- `t3` couples
the two A sites of adjacent cells, i.e. a term whose endpoints have a site
strictly between them, so the string is load-bearing in the Hamiltonian as
well as in the observables. Gapped matters: the uniform half-filled chain
is critical and is this codebase's documented iDMRG convergence trap (see
`tests/test_idmrg_window_free_fermion.py`'s own module docstring).
"""
import numpy as np
import pytest
import scipy.linalg

from dmrgpy import cppext, infinitechain
from dmrgpy.pyitensor import idmrg_window

_FERMION_SITE_CODE = 0  # spinless fermion site (C/Cdag/F/N)

_T1, _T2, _T3 = 1.0, 0.4, 0.1


def _single_particle_matrix(t1, t2, t3, L=200, periodic=True):
    """The one-body hopping matrix h (H = sum_ij h_ij c^dag_i c_j) of the
    dimerized spinless chain, on a ring of L cells (or an open segment of
    the same, `periodic=False`), site index 2*cell+sublattice."""
    N = 2 * L
    h = np.zeros((N, N))
    for n in range(L):
        a, b, a2 = 2 * n, 2 * n + 1, (2 * n + 2) % N
        wraps = (2 * n + 2) >= N
        h[a, b] += t1; h[b, a] += t1
        if periodic or not wraps:
            h[b, a2] += t2; h[a2, b] += t2
            h[a, a2] += t3; h[a2, a] += t3
    return h


def _one_body_density_matrix(h):
    """P[i,j] = <c^dag_i c_j> with every negative single-particle level
    filled."""
    w, v = np.linalg.eigh(h)
    occ = v[:, w < 0]
    return (occ.conj() @ occ.T).real


_CHAINS = {}


def _fermionic_chain(itensor_version="python", maxm=20):
    """A converged, gapped, dimerized spinless infinite chain (n_uc=2).
    Cached per (backend, maxm): every test here reuses one iDMRG solve."""
    key = (itensor_version, maxm)
    if key in _CHAINS:
        return _CHAINS[key]
    ic = infinitechain.Infinite_Many_Body_Chain(
        [_FERMION_SITE_CODE] * 2, itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    ic.maxm = maxm
    ic.maxiter = 40
    C = [ic.get_operator("C", i, "C") for i in range(2)]
    Cd = [ic.get_operator("Cdag", i, "C") for i in range(2)]
    CR = [ic.get_operator("C", i, "R") for i in range(2)]
    CdR = [ic.get_operator("Cdag", i, "R") for i in range(2)]
    ic.set_hamiltonian(_T1 * (Cd[0] * C[1] + Cd[1] * C[0])
                       + _T2 * (Cd[1] * CR[0] + CdR[0] * C[1])
                       + _T3 * (Cd[0] * CR[0] + CdR[0] * C[0]))
    ic.gs_energy()
    _CHAINS[key] = ic
    return ic


def _static_reference(ic, opname_A, opname_B, x, n_uc=2):
    """`<A_x B_0>` from the *static* two_point_correlator, which only ever
    writes its own two operators left-to-right in increasing site order.

    x > 0 puts A at the larger site while leaving it first in the operator
    product, so for two parity-odd operators the reference picks up the
    anticommutation sign `A_x B_0 = -B_0 A_x`. x < 0 needs the sublattice
    position A actually sits on (`x % n_uc`), not 0."""
    odd = opname_A.startswith("C") and opname_B.startswith("C")
    if x > 0:
        val = complex(ic.correlator(opname_B, 0, opname_A, x))
        return -val if odd else val
    if x == 0:
        return complex(ic.correlator(opname_A, 0, opname_B, 0))
    return complex(ic.correlator(opname_A, x % n_uc, opname_B, -x))


@pytest.mark.parametrize("opname_A,opname_B", [("Cdag", "C"), ("C", "Cdag"),
                                                ("N", "N")])
def test_snapshot_at_t0_matches_the_static_correlator(opname_A, opname_B):
    """Oracle 1. Exact identity, both signs of x, across the window's own
    left/right pads (|x| reaches beyond the n_window=4 cells' own centre by
    design), for a fermionic pair in both orders and a parity-even pair
    that must be unaffected by any of this."""
    ic = _fermionic_chain()
    xs = list(range(-4, 5))
    _ts, xarr, S = idmrg_window.dynamical_correlator_td(
        ic._result, 4, opname_A, opname_B, dt=0.05, nt=1, cutoff=1e-10,
        maxdim=40, x_values=xs, connected=False)
    for ix, x in enumerate(xarr):
        assert S[0, ix] == pytest.approx(
            _static_reference(ic, opname_A, opname_B, int(x)), abs=1e-10), \
            "x={}".format(x)


def test_snapshot_at_t0_matches_the_static_correlator_for_a_spin_ladder_operator():
    """The bra-adjoint convention on its own, with no string anywhere near
    it: `Sp`/`Sm` are non-Hermitian but parity-EVEN, so this isolates the
    half of the fix the fermionic tests above cannot separate from the
    string. Applying `A` rather than `A^dag` to the bra silently computes
    <A^dag_x ... B_0>, which for this pair is <Sm_x Sm_0> -- a different,
    generally nonzero object, not a sign flip."""
    from dmrgpy import infinitechain as _ic
    ic = _ic.Infinite_Many_Body_Chain([2, 2], itensor_version="python")
    ic.gs_method = "idmrg"
    ic.maxm = 12
    ic.maxiter = 30
    Sx = [ic.get_operator("Sx", i, "C") for i in range(2)]
    Sy = [ic.get_operator("Sy", i, "C") for i in range(2)]
    Sz = [ic.get_operator("Sz", i, "C") for i in range(2)]
    SxR = [ic.get_operator("Sx", i, "R") for i in range(2)]
    SyR = [ic.get_operator("Sy", i, "R") for i in range(2)]
    SzR = [ic.get_operator("Sz", i, "R") for i in range(2)]
    ic.set_hamiltonian(Sx[0]*Sx[1] + Sy[0]*Sy[1] + Sz[0]*Sz[1]
                       + Sx[1]*SxR[0] + Sy[1]*SyR[0] + Sz[1]*SzR[0])
    ic.gs_energy()
    xs = list(range(-2, 3))
    _ts, xarr, S = idmrg_window.dynamical_correlator_td(
        ic._result, 4, "Sp", "Sm", dt=0.05, nt=1, cutoff=1e-10, maxdim=40,
        x_values=xs, connected=False)
    for ix, x in enumerate(xarr):
        assert S[0, ix] == pytest.approx(
            _static_reference(ic, "Sp", "Sm", int(x)), abs=1e-10), \
            "x={}".format(x)


def test_snapshot_at_t0_matches_the_exact_free_fermion_matrix():
    """Oracle 2: the same numbers against an external free-fermion
    reference, so the identity in the test above cannot be satisfied by two
    consistently-wrong implementations of the same convention."""
    ic = _fermionic_chain()
    P = _one_body_density_matrix(_single_particle_matrix(_T1, _T2, _T3))
    i0 = 2 * (P.shape[0] // 4)          # a bulk A site
    xs = list(range(-4, 5))
    _ts, xarr, S = idmrg_window.dynamical_correlator_td(
        ic._result, 4, "Cdag", "C", dt=0.05, nt=1, cutoff=1e-10, maxdim=40,
        x_values=xs, connected=False)
    for ix, x in enumerate(xarr):
        assert S[0, ix].real == pytest.approx(P[i0 + int(x), i0], abs=1e-6), \
            "x={}".format(x)


def test_time_evolution_matches_the_exact_free_fermion_green_function():
    """Oracle 3, the one that actually tests the string *through* TDVP:
    S(x,t) = <c^dag_x(t) c_0> = sum_l [e^{iht}]_{xl} P[l,0] for a quadratic
    h -- an exact, external number at every time, not just at t=0.

    The tolerance is set from the residual actually observed at this window
    size and bond dimension, per this repo's own practice, not from machine
    precision (t=0 above is the machine-precision anchor). Worst observed
    here is ~3.5e-2 at t=0.4, growing roughly linearly in t: that is the
    finite-window/finite-D error of the method itself, plus the string's
    own left-edge truncation (see `apply_local_operator`'s docstring).
    Before the vacuum normalization went in (see
    `dynamical_correlator_td`'s own docstring) the same numbers were off by
    up to 7e-1 -- a factor ~17 worse, and dominated by a spurious global
    phase rather than by anything physical."""
    ic = _fermionic_chain()
    h = _single_particle_matrix(_T1, _T2, _T3)
    P = _one_body_density_matrix(h)
    i0 = 2 * (h.shape[0] // 4)
    dt, nt = 0.1, 5
    xs = list(range(-2, 3))
    ts, xarr, S = idmrg_window.dynamical_correlator_td(
        ic._result, 6, "Cdag", "C", dt=dt, nt=nt, cutoff=1e-10, maxdim=60,
        x_values=xs, connected=False)
    for it, t in enumerate(ts):
        U = scipy.linalg.expm(1j * h * t)
        for ix, x in enumerate(xarr):
            exact = U[i0 + int(x), :] @ P[:, i0]
            assert S[it, ix] == pytest.approx(exact, abs=6e-2), \
                "x={}, t={}".format(x, t)


def test_dropping_the_string_would_be_a_different_number():
    """Guards the guard: with the parity test forced to False (i.e. the
    pre-fix behaviour restored) the same t=0 values move far beyond every
    tolerance above, so the tests here cannot be passing for some unrelated
    reason."""
    ic = _fermionic_chain()
    xs = list(range(1, 5))
    strung = idmrg_window.dynamical_correlator_td(
        ic._result, 4, "Cdag", "C", dt=0.05, nt=1, cutoff=1e-10, maxdim=40,
        x_values=xs, connected=False)[2][0]
    original = idmrg_window.is_fermionic
    idmrg_window.is_fermionic = lambda name: False
    try:
        bare = idmrg_window.dynamical_correlator_td(
            ic._result, 4, "Cdag", "C", dt=0.05, nt=1, cutoff=1e-10,
            maxdim=40, x_values=xs, connected=False)[2][0]
    finally:
        idmrg_window.is_fermionic = original
    assert max(abs(a - b) for a, b in zip(strung, bare)) > 0.1


def test_odd_total_parity_pair_is_rejected():
    """One fermionic operator against a parity-even one: the string can
    never close, so this is not a well-defined object on an infinite chain
    -- the same contract (and the same failure mode if it were allowed)
    `correlator` already enforces."""
    ic = _fermionic_chain()
    with pytest.raises(ValueError, match="odd total fermion parity"):
        idmrg_window.dynamical_correlator_td(
            ic._result, 4, "Cdag", "N", dt=0.05, nt=1, cutoff=1e-10,
            maxdim=40, x_values=[0, 1], connected=False)
    with pytest.raises(ValueError, match="odd total fermion parity"):
        ic.td_dynamical_correlator("N", 0, "C", n_window=4, dt=0.05, nt=1)


def test_local_expectation_rejects_a_fermionic_operator():
    """<C> is zero by symmetry; the stringless number this used to return
    for it is an artifact, so it raises instead."""
    ic = _fermionic_chain()
    window = idmrg_window.build_window(ic._result, 2)
    with pytest.raises(ValueError, match="parity-odd"):
        idmrg_window.local_expectation(window, ic._result, 2, "Cdag")
    # ...while a parity-even operator on the same fermionic sites is fine
    assert 0.0 < idmrg_window.local_expectation(
        window, ic._result, 2, "N").real < 1.0


@pytest.mark.skipif(not cppext.available(3),
                     reason="mpscpp3 extension not compiled")
def test_itensor_version3_window_matches_the_static_correlator_at_t0():
    """The same exact t=0 identity as the python backend's, on v3.

    This was an xfail until 2026-08-29, and it failed for a reason that
    has nothing to do with the Jordan-Wigner string this module is about:
    the v3 window tiled `idmrg_U_`, whose two ends live in bond bases
    minted by different iDMRG micro-steps, instead of the gauge-consistent
    unit cell that the same file's own `idmrg_onsite_expectation`/
    `idmrg_two_point_correlator` were fixed to use. Measured on a *spin*
    chain, with no fermions anywhere near it, the same identity missed by
    up to 7.4e-2 (x=+1: -0.0864 against an exact -0.1607) while x=0 stayed
    exact -- the signature of a bond-basis mismatch, not of a string. The
    string port itself was unconditionally correct on this backend the
    whole time (a line-for-line mirror of the python one, pinned by the
    tests above); `Chain::idmrg_build_window` now tiles
    `idmrg_cell_raw_`, and this passes."""
    ic = _fermionic_chain(3, maxm=12)
    xs = list(range(-2, 3))
    _ts, _xs, S = ic._session3.td_dynamical_correlator_window(
        4, "Cdag", "C", 0.1, 1, xs, 40, 1e-10, 50, False, 0)
    S = np.array(S).reshape(1, len(xs))
    for ix, x in enumerate(xs):
        assert S[0, ix] == pytest.approx(
            _static_reference(ic, "Cdag", "C", x), abs=1e-6), "x={}".format(x)


@pytest.mark.skipif(not cppext.available(3),
                     reason="mpscpp3 extension not compiled")
def test_itensor_version3_time_evolution_matches_the_green_function():
    """Oracle 3 on v3: the exact free-fermion Green function through TDVP.

    The t=0 test above pins the *gauge*; this one pins everything the
    evolution adds on top of it, against an external number rather than
    against the other backend. It is what caught the second half of the
    2026-08-29 window fix: moving the snapshot's closure onto the cell's
    transfer-matrix fixed points invalidated this backend's old
    `exp(+i*eshift*t)` phase correction (eshift is measured with the
    window's boundary legs *traced*, and the two closings see different
    energies), so S(x,t) was exact at t=0 and drifted linearly afterwards.
    v3 now co-evolves an unperturbed vacuum window and divides by its
    `<psi|psi(t)>`, exactly as `dynamical_correlator_td` does.

    Same tolerance as the python counterpart, and for the same reason:
    what is left is the finite-window/finite-D error of the method plus
    the string's own left-edge truncation, not the machinery."""
    ic = _fermionic_chain(3, maxm=20)
    h = _single_particle_matrix(_T1, _T2, _T3)
    P = _one_body_density_matrix(h)
    i0 = 2 * (h.shape[0] // 4)
    dt, nt = 0.1, 5
    xs = list(range(-2, 3))
    ts, xarr, S = ic._session3.td_dynamical_correlator_window(
        6, "Cdag", "C", dt, nt, xs, 60, 1e-10, 50, False, 0)
    S = np.array(S).reshape(nt, len(xs))
    for it, t in enumerate(ts):
        U = scipy.linalg.expm(1j * h * t)
        for ix, x in enumerate(xarr):
            exact = U[i0 + int(x), :] @ P[:, i0]
            assert S[it, ix] == pytest.approx(exact, abs=6e-2), \
                "x={}, t={}".format(x, t)


def test_kpm_finite_fermionic_matches_the_exact_free_fermion_sum_rule():
    """`kpm_finite` needs no string threading of its own -- it delegates to
    the ordinary finite-chain KPM path, whose operators go through
    `MultiOperator.to_terms()` and so carry the *global* (finite-chain)
    Jordan-Wigner form, the correct convention for the finite object it
    actually diagonalizes. This pins that, since nothing else does: the
    spectral function's own zeroth moment, integral C(w) dw = <0|A B|0>,
    must equal the exact <c^dag_i c_j> of the finite OPEN window
    `kpm_finite` builds -- a number that moves by O(1) if the string is
    dropped, and one this test computes from a plain N x N single-particle
    diagonalization rather than from dmrgpy.

    (Verified along the way, since it is what the sum rule rests on and it
    is nowhere stated: both the finite-chain KPM and the finite-chain ED
    dynamical correlator already carry parity-odd pairs correctly, and the
    operator order is the written one -- `name=(A,B)` integrates to
    `<A B>`, so `("Cdag", "C")` here and not `("C", "C")`, which is
    identically zero and would look like a broken path rather than the
    zero it correctly is.)"""
    ic = _fermionic_chain()
    n_window, r = 6, 3
    es, ys = ic.kpm_finite("Cdag", 0, "C", r, n_window=n_window, n=4000,
                            delta=2e-2, es=np.linspace(-6, 6, 4000))
    integral = np.trapezoid(np.array(ys).real, es)
    L = n_window                       # unit cells, 2 sites each
    h = _single_particle_matrix(_T1, _T2, _T3, L=L, periodic=False)
    P = _one_body_density_matrix(h)
    s_i = (n_window // 2) * 2          # kpm_finite's own central-cell site
    assert integral == pytest.approx(P[s_i, s_i + r], abs=3e-2)
