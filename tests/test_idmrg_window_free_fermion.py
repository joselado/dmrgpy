"""Cross-checks of the IBC-window real-time dynamical correlator
(pyitensor/idmrg_window.py) against an exact, non-interacting (free-
fermion) reference -- see `_free_fermion_reference.py`'s own module
docstring for the derivation and the sign/timing convention it follows.

This module found and fixed a real, previously-undiscovered bug:
`window.env_HL`/`env_HR` are not energy-baseline-subtracted (see
`window_total_energy`'s own docstring), so `window_tdvp_step`'s TDVP
evolution used to carry an arbitrary, macro-iteration-count-dependent
global phase `exp(-i*const*t)` -- confirmed directly, `S(x,t)` for the
*same* physical, equally-converged model came out essentially
uncorrelated between two iDMRG runs that only differed in how many
macro-iterations happened to run before convergence was checked. Fixed
by threading an `eshift` (the window's own ground-state total energy,
measured once before any perturbation) through `window_tdvp_step`,
restoring consistency with the rest of this codebase's own established
"TD" submode convention (`mpscpp3::quench_tdvp`'s `Hshift = H - EGS*Id`,
already used for ordinary finite chains) -- see `window_tdvp_step`'s own
docstring for the full derivation. The native ITensor v3 port
(`Chain::td_dynamical_correlator_window`) had the identical bug and was
fixed the same way (a post-hoc `exp(+i*eshift*t)` correction, `eshift`
measured via a throwaway t=0 `tdvp()` call -- see that method's own
comment in chain_session.h).

`test_window_tdvp_step_eshift_matches_exact_dense_evolution` is the
tightest, most direct validation of the fix itself: it builds the *exact*
dense Hamiltonian matrix for a 2-site window's own (env_HL, window,
env_HR) system (small enough -- local Hilbert space dimension
`chi_L*2*2*chi_R` at modest bond dimension -- for `np.linalg.eigh` and a
literal spectral-decomposition time evolution) and compares
`window_tdvp_step`'s own result against it directly, with no free-fermion
mapping or iDMRG-convergence-quality confound at all. This is the
strongest evidence the fix is not just "less wrong" but numerically
exact for whatever `env_HL`/`env_HR`/window state it is given.

The free-fermion comparison further downstream also reflects iDMRG's own
ground-state convergence quality (not just this module's own
correctness) -- confirmed directly: the *gapless*, exactly-half-filled
uniform XX chain (`J2=None`) has iDMRG state_overlap topping out around
0.87-0.96 even with many restart attempts (matching this codebase's own
documented gapless-iDMRG convergence caveats, e.g.
`test_idmrg_window.py`'s own docstrings), which shows up as a real,
growing-with-t discrepancy against the exact free-fermion answer despite
`window_tdvp_step` itself being verified exact above -- i.e. the
discrepancy is in *how well iDMRG's own U_list/env_HL/env_HR approximate
the true infinite chain*, not in this module's own TDVP code. Switching
to a *dimerized* (SSH-like) free-fermion chain (`J2 != J1`, still exactly
solvable the same way, see `_free_fermion_reference.py`) opens a gap and
gets state_overlap up to ~0.98-0.999 with the same modest bond dimension
-- used here for the quantitative comparison, with tolerances set from
the residual actually observed at that convergence quality, not
machine precision.

A LATER CORRECTION to that story, since it is the kind of note a reader
will otherwise trust: the residual observed at the time (~0.07) was not
mostly convergence quality. It was mostly a second, independent error --
a closure mismatch between the `eshift` measured by
`window_total_energy` (boundary legs traced) and the contraction
`snapshot_correlator` measures through (boundary legs closed by the
transfer-matrix fixed points), leaving a spurious factor of ~1.6 rad per
unit time on the evolved branch. It was found on 2026-08-29 by a sharper
oracle (a fermionic chain, where the same defect showed up as 3.0 rad/t
against an exact Green function) and is now cancelled exactly, by
dividing every `S(x,t)` by the vacuum amplitude measured through the
identical contraction -- see `idmrg_window.window_total_energy`'s own
docstring, point 3. The genuine convergence residual, now visible
underneath it, is ~1e-5 at these parameters. The convergence *ordering*
above (gapless chain worse than dimerized) is still real and still the
reason this module uses a dimerized model.
"""
import numpy as np
import pytest

from dmrgpy import spinchain, timedependent
from dmrgpy.pyitensor import idmrg, idmrg_window, kernels
from dmrgpy.pyitensor.mpsalgebra import inner

from _free_fermion_reference import FreeFermionXX


def test_free_fermion_reference_matches_ed():
    """Validates `_free_fermion_reference.py`'s own formula (not
    idmrg_window.py at all) against exact many-body ED on a small (8-site)
    open XX chain, at 8 sites' worth of `x` and several `t` -- the
    strongest possible check of the reference's own sign/timing
    conventions, since ED is exact regardless of system size for a chain
    this small. See `_free_fermion_reference.py`'s own docstring for why
    a complex conjugate is needed here (a sign convention internal to
    edtk/tdtk.py's `evolve()`, unrelated to idmrg_window.py).

    This used to additionally shift `ts_ed` by `+dt` to compensate for the
    ED backend's own "evolve-then-measure" off-by-one -- `cs[0]` was the
    value at `t=dt` while `ts[0]` said `t=0`. That was fixed at the source
    (every backend now measures before evolving, see timedependent.py's
    evolution_dmrg_DC docstring), so the shift is gone and `ts_ed` is used
    as-is; this test would fail by ~6e-3 if either the fix or this
    compensating shift were reverted on its own."""
    n, J = 8, 1.0
    sc = spinchain.Spin_Chain([2 for _ in range(n)])
    h = 0
    for i in range(n - 1):
        h = h + J * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1])
    sc.set_hamiltonian(h)
    ff = FreeFermionXX(n, J=J)

    dt, nt = 0.02, 15
    y0 = 3
    mean_y0 = sc.vev(sc.Sz[y0], mode="ED").real
    max_err = 0.0
    for x in range(n):
        ts_ed, cs_ed = timedependent.evolution_DC(
            sc, mode="ED", name=(sc.Sz[x], sc.Sz[y0]), nt=nt, dt=dt)
        mean_x = sc.vev(sc.Sz[x], mode="ED").real
        conn_ed = np.array(cs_ed) - mean_x * mean_y0

        conn_ff = np.conj(np.array([ff.sz_sz_connected(x, y0, t) for t in np.array(ts_ed)]))
        max_err = max(max_err, np.max(np.abs(conn_ed - conn_ff)))
    assert max_err < 1e-9


def _full_window_order(win):
    """(order_in, shape) for the *entire* window's own explicit tensors,
    as a single flat vector -- [left_link, s_1, ..., s_n, right_link].
    Rebuilt fresh on demand since site Index identities are freshly
    re-minted by every `window_tdvp_step` call (see this module's own
    per-step rebuilding elsewhere)."""
    n = win.mps.length()
    site_inds = [next(ind for ind in win.mps.A(i).inds if ind.hastags("Site"))
                 for i in range(1, n + 1)]
    order = [win.mps.A(1).inds[0]] + site_inds + [win.mps.A(n).inds[-1]]
    shape = tuple(ind.dim for ind in order)
    return order, shape


def _full_window_vector(win):
    order, _shape = _full_window_order(win)
    theta = win.mps.A(1)
    for i in range(2, win.mps.length() + 1):
        theta = theta * win.mps.A(i)
    return theta.transpose_to(order).reshape(-1)


def test_window_tdvp_step_eshift_matches_exact_dense_evolution():
    """The decisive check for the `eshift` fix itself: an n_window=3
    window's own (env_HL, window, env_HR) system is small enough to
    diagonalize exactly (dense `np.linalg.eigh` over the *whole* window,
    not just one two-site bond), giving a literal, closed-form
    e^{-i(H-EGS)t} to compare `window_tdvp_step`'s own full multi-step
    sweep against directly -- no free-fermion mapping, no iDMRG-
    convergence-quality confound. Deliberately n_window=3 (not 2): a
    2-site window has only one bond, whose own L/R *are* env_HL/env_HR
    directly with no `_extend_HL`/`_extend_HR` chaining at all -- it
    would silently miss a bug specific to the interior-environment
    extension `_all_left_environments_window`/`_all_right_environments_
    window` build and `window_tdvp_step`'s own multi-bond sweep actually
    exercises (confirmed as a real test-coverage gap, not hypothetical,
    by code review). If this passes tightly, the fix is numerically exact
    for whatever env_HL/env_HR/window state it is given (any *remaining*
    discrepancy in a free-fermion comparison downstream must come from
    iDMRG's own ground-state convergence quality, not this module's own
    TDVP code)."""
    field = 1.4
    best = None
    for seed in range(30):
        np.random.seed(seed)
        result = idmrg.idmrg_ground_state(
            [2], [[field, ["Sz", 0]]], [[1.0, ["Sx", 0], ["Sx", 1]]], 1,
            maxm=16, cutoff=1e-12, maxiter=150, etol=1e-14, niter=100, verbose=False)
        if not result.converged:
            continue  # unseeded local solve occasionally misses even the
                      # energy-density criterion in a given attempt -- skip
                      # it rather than let it win "best" on overlap alone
        if best is None or (result.state_overlap or 0) > (best.state_overlap or 0):
            best = result
        if (best.state_overlap or 0) > 0.999:
            break
    assert best is not None and best.converged

    n_window = 3  # n_uc=1: build_window rounds this up to a whole number of
    # tiling cells (2 unit cells each when n_uc=1, see
    # idmrg_window._window_cell), so the realized window is 4 sites / 3 bonds
    # -- still the interior-environment extension path this test is for (see
    # its own docstring), just one site wider than the 3 originally asked
    # for. Everything below reads the site count off the window itself.
    win = idmrg_window.build_window(best, n_window)
    n_sites = win.mps.length()

    # Build the exact dense n_window-site effective Hamiltonian: env_HL
    # contracted through every site's own MPO tensor, closed by env_HR --
    # kernels.make_matvec is a generic tensor-network contraction, not
    # hardcoded to a fixed piece count, so this generalizes
    # _window_two_site_heff's own 2-piece pattern to n_window pieces.
    order_in, shape = _full_window_order(win)
    site_inds_out = [ind.prime(1) for ind in order_in[1:-1]]
    order_out = [win.env_HL_bra] + site_inds_out + [win.env_HR_bra]
    pieces = [win.env_HL] + [win.mpo.A(i) for i in range(1, n_sites + 1)] + [win.env_HR]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    dim = int(np.prod(shape))
    Hdense = np.zeros((dim, dim), dtype=complex)
    basis = np.eye(dim, dtype=complex)
    for k in range(dim):
        Hdense[:, k] = matvec(basis[k])
    assert np.max(np.abs(Hdense - Hdense.conj().T)) < 1e-10

    eshift = idmrg_window.window_total_energy(win)
    Hshift = Hdense - eshift * np.eye(dim)
    evals, evecs = np.linalg.eigh(Hshift)

    # Perturb the center site and read off psi0.
    center = 2
    idmrg_window.apply_local_operator(win, best, center, "Sz")
    psi0 = _full_window_vector(win).astype(complex)

    dt, nt = 0.05, 6
    ts = [it * dt for it in range(nt)]
    psi_exact = []
    for t in ts:
        phase = np.exp(-1j * evals * t)
        psi_exact.append(evecs @ (phase * (evecs.conj().T @ psi0)))

    psi_tdvp = []
    for it in range(nt):
        psi_tdvp.append(_full_window_vector(win))
        if it < nt - 1:
            idmrg_window.window_tdvp_step(win, dt, cutoff=1e-12, maxdim=40, niter=100, eshift=eshift)

    for it in range(nt):
        assert np.linalg.norm(psi_exact[it] - psi_tdvp[it]) < 1e-4


def _dimerized_setup(n_window, J1=1.0, J2=0.2, maxm=24, attempts=30,
                      overlap_target=0.99, min_overlap=0.9, maxiter=80, niter=60):
    """Converged dimerized-XX iDMRG result with `state_overlap>min_overlap`
    that also actually passes `build_window(n_window)` -- a high
    `state_overlap` alone does not guarantee the periodic-wraparound bond
    dimensions match (see `build_window`'s own RuntimeError message,
    confirmed directly to trigger even at state_overlap>0.98 for this
    model), and a candidate that *does* pass `build_window` is not
    necessarily well-converged (confirmed directly: with a smaller
    `attempts` budget, the highest-overlap `build_window`-passing
    candidate was sometimes only ~0.7-0.8, since the search stops at the
    first passing candidate by overlap rank rather than requiring a
    floor) -- both conditions are checked here, over a large enough
    `attempts` budget that a qualifying candidate is reliably found
    (iDMRG's own local two-site solve is unseeded -- see idmrg.py's own
    `_local_two_site_solve` docstring -- so run-to-run convergence
    quality genuinely varies)."""
    hi = [[J1, ["Sx", 0], ["Sx", 1]], [J1, ["Sy", 0], ["Sy", 1]]]
    he = [[J2, ["Sx", 1], ["Sx", 2]], [J2, ["Sy", 1], ["Sy", 2]]]
    candidates = []
    for seed in range(attempts):
        np.random.seed(seed)
        result = idmrg.idmrg_ground_state([2, 2], hi, he, 2, maxm=maxm, cutoff=1e-12,
                                           maxiter=maxiter, etol=1e-14, niter=niter, verbose=False)
        candidates.append((result.state_overlap or 0, result))
        candidates.sort(key=lambda p: -p[0])
        if candidates[0][0] > overlap_target:
            break
    for ov, result in candidates:
        if ov < min_overlap:
            break  # sorted descending -- no further candidate can qualify
        try:
            idmrg_window.build_window(result, n_window)
            return result, J1, J2
        except RuntimeError:
            continue
    pytest.skip("no candidate iDMRG solve reached state_overlap>{} while "
                "also passing build_window's own self-consistency check "
                "in {} attempts".format(min_overlap, attempts))


def test_dynamical_correlator_td_reproducible_across_convergence():
    """The bug this module fixed: pre-fix, `S(x,t)` for the *same*
    physical, equally-converged model came out essentially uncorrelated
    between iDMRG runs that only differed in macro-iteration count (an
    arbitrary, non-reproducible run-to-run baseline in the unshifted
    env_HL/env_HR). Two independent, well-converged solves of the
    dimerized model here should now give closely-agreeing S(x,t), not
    just a coincidentally-agreeing e0."""
    n_window, dt, nt = 8, 0.05, 4
    result_a, J1, J2 = _dimerized_setup(n_window)
    result_b, _, _ = _dimerized_setup(n_window)
    assert result_a.e0 == pytest.approx(result_b.e0, abs=1e-3)

    ts_a, xs_a, S_a = idmrg_window.dynamical_correlator_td(
        result_a, n_window=n_window, opname_A="Sz", opname_B="Sz", dt=dt, nt=nt,
        cutoff=1e-10, maxdim=28, niter=30, x_values=[0, 1], connected=True)
    ts_b, xs_b, S_b = idmrg_window.dynamical_correlator_td(
        result_b, n_window=n_window, opname_A="Sz", opname_B="Sz", dt=dt, nt=nt,
        cutoff=1e-10, maxdim=28, niter=30, x_values=[0, 1], connected=True)
    # Loose tolerance (independent ground states, independent gauges) --
    # the point is that these are now close at all, not exact agreement.
    # Pre-fix, this max|diff| was order-1 (comparable to S itself); the
    # fixed version stays a small fraction of that.
    assert np.max(np.abs(S_a - S_b)) < 0.15


def test_dynamical_correlator_td_matches_dimerized_free_fermion():
    """Quantitative comparison against the exact free-fermion answer
    (`_free_fermion_reference.py`) for the dimerized XX chain, on a
    reference system (N=2000) far larger than many-body ED could reach.
    `t=0` is a tight check (the already-exact static correlator, no TDVP
    involved -- see `test_idmrg_window.py`'s own analogous check);
    `t>0` uses a tolerance set from the residual actually observed at
    this convergence quality. NOTE that until 2026-08-29 that residual was
    ~0.07 and was attributed, here and in this module's own docstring, to
    iDMRG's own ground-state approximation. That was wrong: it was
    dominated by a second, independent phase error (a closure mismatch
    between `eshift` and the correlator's own contraction, see
    `idmrg_window.window_total_energy`'s own docstring point 3), and it
    dropped to ~1e-5 the moment `dynamical_correlator_td` started
    dividing by the vacuum amplitude. `window_tdvp_step` itself was never
    at fault -- that part of the old note stands, and is verified exactly
    by `test_window_tdvp_step_eshift_matches_exact_dense_evolution`."""
    n_window, dt, nt = 8, 0.05, 5
    result, J1, J2 = _dimerized_setup(n_window, overlap_target=0.995)

    ts, xs, S = idmrg_window.dynamical_correlator_td(
        result, n_window=n_window, opname_A="Sz", opname_B="Sz", dt=dt, nt=nt,
        cutoff=1e-10, maxdim=40, niter=50, x_values=[-1, 0, 1], connected=True)

    ff = FreeFermionXX(2000, J=J1, J2=J2)
    x0 = 1000
    xs = list(xs)

    ix0 = xs.index(0)
    assert S[0, ix0].real == pytest.approx(0.25, abs=1e-6)
    assert S[0, ix0].imag == pytest.approx(0.0, abs=1e-6)

    max_err = 0.0
    for it, t in enumerate(ts):
        for ix, x in enumerate(xs):
            ff_val = ff.sz_sz_connected(x0 + int(x), x0, t)
            max_err = max(max_err, abs(S[it, ix] - ff_val))
    # 2e-3, three orders of magnitude below the bound this test used to
    # carry. The old 0.12 was set from an observed 0.069-0.070 residual
    # attributed to iDMRG's own ground-state convergence quality -- that
    # attribution was wrong. Almost all of it was a spurious global factor
    # from a closure mismatch: `eshift` (window_total_energy) is measured
    # with the window's boundary legs traced, while the correlator is
    # measured with them closed by the transfer-matrix fixed points, and
    # the two see different energies, so the evolved branch rotated at
    # ~1.6 rad per unit time relative to the un-evolved one (the arithmetic
    # checks out: 1.6 x t=0.2 x |S|~0.25 ~ 0.08). Since
    # dynamical_correlator_td started dividing by the vacuum amplitude
    # measured through the identical contraction (see its own docstring),
    # max_err here measures 1.1e-5 at t<=0.2 and 2.8e-5 out to t=0.9 --
    # i.e. the genuine convergence residual was always ~4 orders of
    # magnitude smaller than what this bound was absorbing. 2e-3 leaves
    # room for run-to-run variation in an unseeded iDMRG solve while still
    # catching any reintroduction of a phase error of that class.
    assert max_err < 2e-3
