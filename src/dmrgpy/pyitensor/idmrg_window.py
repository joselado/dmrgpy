"""Heterogeneous, infinite-boundary-condition (IBC) window construction for
iDMRG (pyitensor backend), following Milsted/Vanderstraeten et al.,
"Infinite boundary conditions for response functions and limit cycles in
iDMRG" (arXiv:1804.09163), Sec. V.1.

The paper's method computes real-time dynamical correlators of an infinite
chain by (1) freezing a uniform iMPS ground state everywhere except a
finite window of `n_window` unit cells, (2) capping the two window edges
with "the Hamiltonian terms outside the window, projected onto the reduced
D-dimensional Hilbert space of the left/right block", and (3) perturbing
and time-evolving only the window.

The key realization this module relies on: an ordinary iDMRG growth
environment (idmrg.py's own HL/HR, accumulated every macro-iteration) *is
already* exactly that projection -- no new environment solve is needed, only
plumbing to reuse it. `pyitensor/idmrg.py`'s `idmrg_ground_state` was
extended (see its own `env_window_boundary` comment and IDMRGResult's own
docstring) to snapshot HL/HR (and the automaton tensors W_bulk) as they
stood *entering* the last macro-iteration -- precisely the snapshot whose
own bond Indices match the *returned* U_list's own edge bonds (U_list[0]'s
left bond and U_list[n_uc-1]'s right bond respectively), so a window can be
built and capped with zero further approximation beyond what iDMRG's own
convergence already assumes for U_list itself.

This module builds and validates the capped window: an explicit, finite
`mpscontainer.MPS`/`MPO` pair (so `.A(i)`/`.length()`/etc. all work
normally), with two extra environment tensors (`.env_HL`/`.env_HR`) that
must cap its two open ends. Environment extension *through* the window
deliberately reuses idmrg.py's own explicit-Index `_extend_HL`/`_extend_HR`
(see `_extend_through_left`/`_extend_through_right` below) rather than
dmrg.py's generic, `_Chain`/`_link_at`-based `_extend_left`/`_extend_right`:
the latter infers a site's left/right Link purely by looking for a
same-Index neighbor in the *chain itself* (`_link_at`), which is correct
for an ordinary finite chain (whose two boundary sites truly have no
further leg at all, see mpscontainer.py's own docstring) but wrong here --
the window's own edge sites *do* carry one more leg each (the one
connecting to `env_HL`/`env_HR`, tensors that live outside the `_Chain`
object entirely), which `_link_at` cannot see and would silently treat as
absent. idmrg.py's own extend functions take the relevant Index objects
explicitly instead of inferring them, exactly sidestepping this -- fully
consistent, since the window's MPS/MPO tensors are themselves built (via
`_tile_periodic`) with the exact same per-site tensor shape idmrg.py's own
growth loop already uses.

Follow-up work (perturb the window, TDVP-evolve it, read off shifted
overlaps, Fourier transform -- Sec. V.1 steps 3-5 of the paper) will need
an analogous window-aware two-site sweep (mirroring tdvp.py's own
structure but built on these same explicit-Index extend primitives rather
than tdvp.py's own `_link_at`-based half-sweeps) -- out of scope here. This
module's own `window_energy_density` is the sanity check that validates
the construction so far: the window's own marginal energy density (a
finite difference between two window sizes -- see its own docstring for
why not the absolute total directly), computed via idmrg.py's independent
`_extend_HL`/`_extend_HR` code (not `idmrg_ground_state`'s own internal
accumulation this snapshot was taken from), must reproduce `result.e0`,
and must do so *without* the boundary contamination an ordinary
open-boundary window (`infinitechain.py`'s
`kpm_finite`/`_window_hamiltonian`) suffers from regardless of how large
`n_window` is."""

from . import idmrg as _idmrg_mod
from .mpscontainer import MPO, MPS
from .tensor import ITensor


def _tile_periodic(tensors_uc, boundary_left, boundary_right, n_window):
    """`n_window` periodic repeats of `tensors_uc` (a length-n_uc list of
    rank>=3 ITensors, each shaped (left Link, ...middle legs..., right
    Link), with consecutive tensors already sharing Link identity within
    one "natural" copy -- true of both idmrg.py's own U_list (built within
    a single macro-iteration's own sweep, see idmrg_ground_state's own
    comments) and W_bulk (built directly with shared boundary_idx Indices,
    see _build_periodic_mpo's own docstring)), returned as a flat length-
    n_window*len(tensors_uc) list.

    Every copy boundary (including the very first tensor's own left leg and
    the very last tensor's own right leg) gets an explicitly assigned Link
    Index rather than reusing `tensors_uc`'s own natural (period-n_uc,
    reused-by-identity-across-copies) boundary Indices -- reusing those
    directly across repeats would silently alias non-adjacent window sites
    onto the same Index object, exactly the class of bug idmrg.py's own
    `_relabel_pos`/`_project_channel` docstrings warn about at length.
    `boundary_left`/`boundary_right` are used as-is for the window's own
    two open ends (so they attach directly to whatever environment cap the
    caller supplies, e.g. IDMRGResult's env_HL_ket/env_HR_ket or
    env_HL_mpo/env_HR_mpo) -- every other copy-to-copy cut is a freshly
    minted Index (`.sim()` of the corresponding natural boundary Index, so
    the dimension is right by construction)."""
    n_uc = len(tensors_uc)
    out = []
    left_bound = boundary_left
    for c in range(n_window):
        last_copy = (c == n_window - 1)
        internal = [tensors_uc[p].inds[-1].sim() for p in range(n_uc - 1)]
        right_bound = boundary_right if last_copy else tensors_uc[-1].inds[-1].sim()
        legs_left = [left_bound] + internal
        legs_right = internal + [right_bound]
        for p in range(n_uc):
            T = tensors_uc[p]
            new_inds = (legs_left[p],) + T.inds[1:-1] + (legs_right[p],)
            out.append(ITensor(new_inds, T.array))
        left_bound = right_bound
    return out


class IBCWindow:
    """A finite, `n_window*n_uc`-site MPS (`.mps`) and matching MPO
    (`.mpo`) -- both ordinary `mpscontainer` chains (`.A(i)`/`.length()`/
    etc. all work normally) -- plus the converged environment caps
    (`.env_HL`/`.env_HL_bra` on the left, `.env_HR`/`.env_HR_bra` on the
    right) that must seed environment extension *through* this window for
    every subsequent computation (ground-state energy, perturbation, TDVP
    evolution, ...) to see the correct infinite-boundary physics. That
    extension needs idmrg.py's own explicit-Index `_extend_HL`/`_extend_HR`
    (see `_extend_through_left`/`_extend_through_right` below and this
    module's own docstring) -- *not* dmrg.py's/tdvp.py's generic
    `_Chain`/`_link_at`-based environment machinery, which cannot see the
    window's own edge sites' extra leg (the one connecting to `env_HL`/
    `env_HR`) and would silently mis-contract it. Not meant to be
    constructed directly -- see `build_window`."""

    def __init__(self, mps, mpo, env_HL, env_HL_bra, env_HR, env_HR_bra, n_uc, n_window):
        self.mps = mps
        self.mpo = mpo
        self.env_HL = env_HL
        self.env_HL_bra = env_HL_bra
        self.env_HR = env_HR
        self.env_HR_bra = env_HR_bra
        self.n_uc = n_uc
        self.n_window = n_window


def build_window(result, n_window):
    """Build an `n_window`-unit-cell IBCWindow from a converged
    `pyitensor.idmrg.IDMRGResult` (`itensor_version="python"` only --
    `result.U_list`/`result.W_bulk`/`result.env_*` are all pyitensor-only
    objects).

    `n_window` should be chosen large enough that the causal cone of
    whatever perturbation/time-evolution is later applied at the window's
    center never reaches its two open ends within the simulated time --
    same convergence caveat as `infinitechain.py`'s `kpm_finite` own
    `n_window` (check by increasing it and confirming results stop
    changing), except here the two ends are physically correct (IBC-capped,
    see this module's own docstring) rather than open-boundary artifacts,
    so a *much* smaller margin should suffice in practice than an
    open-boundary window of comparable accuracy would need."""
    if n_window < 1:
        raise ValueError("build_window: n_window must be >= 1")
    if result.env_HL_ket is None or result.env_HR_ket is None:
        raise RuntimeError(
            "build_window: result has no converged environment snapshot "
            "(env_HL/env_HR) to cap a window with -- idmrg_ground_state "
            "must complete at least one full macro-iteration first (always "
            "true given maxiter>=2 is already enforced there, so this "
            "should not happen in practice)")
    if n_window > 1 and result.U_list[-1].inds[-1].dim != result.U_list[0].inds[0].dim:
        # Tiling more than one copy requires U_list[n_uc-1]'s own right bond
        # to be reused (via a fresh Index of the same dimension, see
        # _tile_periodic) as U_list[0]'s own left bond in the next copy --
        # physically the *same* recurring cut of a truly periodic MPS, and
        # equal by construction for a genuinely self-consistent, converged
        # unit cell (the same equality _transfer_matrices/_compose already
        # rely on for the wraparound bond in onsite_expectation/
        # two_point_correlator). A mismatch here means U_list has not
        # actually settled into a self-consistent, translationally-
        # invariant state yet (bond dimension hasn't uniformly saturated
        # across the periodic boundary) -- confirmed directly: this can
        # happen even when idmrg_ground_state reports .converged=True
        # (a pure energy-density criterion) with a less-than-excellent
        # .state_overlap, exactly the gap IDMRGResult.state_overlap's own
        # docstring warns about. Check state_overlap and/or increase
        # maxiter/niter/maxm rather than treat this as a bug in this
        # module.
        raise RuntimeError(
            "build_window: U_list[-1]'s own right-bond dimension ({}) "
            "does not match U_list[0]'s own left-bond dimension ({}) -- "
            "a multi-copy window needs these to agree (see this "
            "function's own comment). This means the converged unit cell "
            "is not yet self-consistent across its own periodic boundary "
            "(check result.state_overlap, or increase maxiter/niter/maxm) "
            "-- not a bug in idmrg_window.py itself.".format(
                result.U_list[-1].inds[-1].dim, result.U_list[0].inds[0].dim))
    ket_tensors = _tile_periodic(result.U_list, result.env_HL_ket,
                                  result.env_HR_ket, n_window)
    mpo_tensors = _tile_periodic(result.W_bulk, result.env_HL_mpo,
                                  result.env_HR_mpo, n_window)
    mps = MPS(ket_tensors)
    mps.center = None  # not canonical relative to any single center yet
    mpo = MPO(mpo_tensors)
    mpo.center = None
    return IBCWindow(mps, mpo, result.env_HL, result.env_HL_bra,
                      result.env_HR, result.env_HR_bra, result.n_uc, n_window)


def _extend_through_left(window, up_to_site):
    """(HL, HL_bra) after extending window.env_HL through window sites
    1..up_to_site (inclusive), via idmrg.py's own explicit-Index
    `_extend_HL` (see this module's own docstring for why, not dmrg.py's
    `_extend_left`). Each site's own "old" (already-part-of-HL) ket leg and
    "new" (freshly dangling) ket leg are read directly off the window's own
    tensors -- known exactly, by construction (`_tile_periodic`), rather
    than inferred via any chain-boundary lookup."""
    ket, mpo = window.mps, window.mpo
    HL, HL_bra = window.env_HL, window.env_HL_bra
    left_ket_old = ket.A(1).inds[0]  # == window.env_HL's own ket leg
    for i in range(1, up_to_site + 1):
        U, W_p = ket.A(i), mpo.A(i)
        right_ket_new = U.inds[-1]
        HL, HL_bra = _idmrg_mod._extend_HL(HL, HL_bra, W_p, U, left_ket_old, right_ket_new)
        left_ket_old = right_ket_new
    return HL, HL_bra


def _extend_through_right(window, down_to_site):
    """Mirror of `_extend_through_left`: (HR, HR_bra) after extending
    window.env_HR through window sites n..down_to_site (inclusive, sweeping
    right to left), via idmrg.py's own `_extend_HR`."""
    ket, mpo = window.mps, window.mpo
    n = ket.length()
    HR, HR_bra = window.env_HR, window.env_HR_bra
    right_ket_old = ket.A(n).inds[-1]  # == window.env_HR's own ket leg
    for i in range(n, down_to_site - 1, -1):
        V, W_p = ket.A(i), mpo.A(i)
        left_ket_new = V.inds[0]
        HR, HR_bra = _idmrg_mod._extend_HR(HR, HR_bra, W_p, V, right_ket_old, left_ket_new)
        right_ket_old = left_ket_new
    return HR, HR_bra


def _relabel_index(T, old_ind, new_ind):
    """T with every occurrence of `old_ind` (matched by identity+prime
    level) replaced by `new_ind`. Safe here via plain identity matching
    (unlike idmrg.py's own positional `_relabel_pos`): every leg this
    module ever relabels is guaranteed distinct within its own tensor,
    since window tensors are always built with fresh per-copy Link legs
    (see `_tile_periodic`) rather than idmrg.py's own automaton tensors,
    which can have identical left/right legs for small n_uc."""
    new_inds = tuple(new_ind if ind == old_ind else ind for ind in T.inds)
    return ITensor(new_inds, T.array)


def window_total_energy(window):
    """Total (unnormalized-by-site-count) energy of the capped window,
    computed via idmrg.py's own independent `_extend_HL` code (called fresh
    here through the *whole* window, not reusing whatever accumulation
    `idmrg_ground_state` itself produced this snapshot from) -- a genuine
    cross-check of the whole construction (Index plumbing included), not a
    tautological re-derivation.

    This is *not* simply `e0 * (number of window sites)`, for two reasons,
    both confirmed directly (an earlier version of this function skipped
    both and reported energies off by 2-3 orders of magnitude for a
    Heisenberg test case -- e.g. -1014 instead of an expected ~-8.4 for a
    19-site configuration):

    1. `<window|window>` (using env_HL/env_HR and the window's own
       left-canonical tensors exactly as built, with *both* env_HL_ket and
       env_HR_ket left as free/dangling legs rather than contracted against
       any specific boundary vector) is *not* 1 -- it equals
       `dim(env_HR_ket)` exactly, by the isometry property of left-canonical
       tensors chained together (each site's own U being an isometry from
       (its left leg, physical leg) to its right leg collapses the sum over
       every index except the very last dangling one, leaving that
       dimension as the Frobenius norm-squared) -- this function divides by
       it, matching a Rayleigh-quotient convention (confirmed numerically:
       matches an independent Lanczos-Rayleigh-quotient computation using
       the exact same env_HL/env_HR/W_pL/W_pR to machine precision).
    2. Even after that fix, the *absolute* total still includes whatever
       (large, `n_window`-independent, and not independently exposed)
       energy env_HL/env_HR themselves represent from every macro-iteration
       already absorbed before the snapshot -- so this total is only
       meaningful as a *difference* between two window sizes (see
       `window_energy_density`), not compared directly against `e0` on its
       own."""
    n = window.mps.length()
    L, Lbra = _extend_through_left(window, n)
    # L's three dangling legs are (Lbra [a fresh bra-side mint from this
    # very last extension], env_HR_ket, env_HR_mpo) -- the ket/mpo legs
    # already match window.env_HR's own by construction (the window's last
    # site's own right legs *are* env_HR_ket/env_HR_mpo, see
    # _tile_periodic/build_window), only the bra leg needs relabeling
    # before the two can be multiplied down to a bare scalar.
    env_HR_relabeled = _relabel_index(window.env_HR, window.env_HR_bra, Lbra)
    closed = L * env_HR_relabeled
    norm = window.mps.A(n).inds[-1].dim  # == dim(env_HR_ket), see docstring
    return closed.scalar().real / norm


def window_energy_density(result, n_window):
    """The marginal energy *density* of growing the window by one more unit
    cell -- `(window_total_energy(n_window+n_uc) -
    window_total_energy(n_window)) / n_uc` -- which should reproduce
    `result.e0` (to within iDMRG's own convergence tolerance), for *any*
    `n_window` -- unlike an open-boundary window
    (`infinitechain.py`'s `kpm_finite`), which only approaches `e0` in the
    large-`n_window` limit and never removes edge contamination entirely
    (see that module's own docstring).

    This finite difference (mirroring `idmrg_ground_state`'s own
    `density = (energy - prev_energy) / (2*n_uc)` diagnostic) is the
    well-posed sanity check for this construction, not
    `window_total_energy` alone: the latter's own large, env_HL/env_HR-
    derived baseline (see its own docstring) is a fixed, `n_window`-
    independent constant that this difference exactly cancels. Confirmed
    numerically against several gapped test models (n_uc=1 and n_uc=2):
    this reproduces `e0` to within iDMRG's own residual convergence error
    (state_overlap not exactly 1 -- typically a few percent at maxm~20-24,
    matching how far state_overlap itself is from 1), converging cleanly
    (no further drift) by n_window~3-4, while `window_total_energy` divided
    by site count alone drifts with no clear limit. A genuinely gapless or
    poorly-converged `result` (low `state_overlap`, e.g. n_uc=1 gapless
    models -- a documented, pre-existing iDMRG limitation, see
    `idmrg_warmstart_correlator_fix`-era notes) will *not* reproduce `e0`
    well here either, for the same underlying reason static correlators
    already carry that caveat -- check `result.state_overlap` first.

    `n_window` counts *unit-cell copies* (see `build_window`'s own
    docstring -- an IBCWindow has `n_window*n_uc` physical sites), so
    growing by exactly one more copy (`n_window + 1`, not `n_window +
    n_uc`) adds exactly `n_uc` physical sites -- matching the `n_uc`
    divisor below. An earlier version of this function used `n_window +
    n_uc` here, silently adding `n_uc` *copies* (`n_uc*n_uc` physical
    sites) while still dividing by only `n_uc` -- invisible for n_uc=1
    (where n_uc*n_uc==n_uc), but confirmed directly to inflate every n_uc=2
    result by *exactly* a factor of 2 (density -2.278 instead of -1.139 for
    an easy, cleanly-converged alternating-field test model -- an exact
    ratio, not convergence noise, which is what exposed this as a real
    off-by-n_uc bug rather than an iDMRG convergence limitation)."""
    n_uc = result.n_uc
    win_a = build_window(result, n_window)
    win_b = build_window(result, n_window + 1)
    e_a = window_total_energy(win_a)
    e_b = window_total_energy(win_b)
    return (e_b - e_a) / n_uc
