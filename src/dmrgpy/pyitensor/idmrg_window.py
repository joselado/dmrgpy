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

`window_energy_density` is the sanity check that validates the static
construction: the window's own marginal energy density (a finite
difference between two window sizes -- see its own docstring for why not
the absolute total directly), computed via idmrg.py's independent
`_extend_HL`/`_extend_HR` code (not `idmrg_ground_state`'s own internal
accumulation this snapshot was taken from), must reproduce `result.e0`,
and must do so *without* the boundary contamination an ordinary
open-boundary window (`infinitechain.py`'s
`kpm_finite`/`_window_hamiltonian`) suffers from regardless of how large
`n_window` is.

== Real-time evolution of the window (Sec. V.1 steps 3-4 of the paper) ==

`window_tdvp_step` two-site-TDVP-evolves a window in place, exactly
mirroring `tdvp.py`'s own `tdvp_step` (a left-to-right half-sweep by dt/2,
then a right-to-left half-sweep by another dt/2) -- but `tdvp.py`'s own
sweep functions cannot be reused directly: they build each bond's local
effective Hamiltonian (`dmrg.py`'s `two_site_heff`/`one_site_heff`) and
environments (`_all_left_environments`/`_all_right_environments`,
`_extend_left`/`_extend_right`) via `mpsalgebra.py`'s `_link_at`, which
finds a site's left/right Link purely by looking for a same-Index
*neighbor in the chain itself* -- correct for an ordinary finite chain
(whose two boundary sites truly have no further leg, see
`mpscontainer.py`'s own docstring) but wrong for a window, whose edge
sites *do* carry one more leg each (connecting to `env_HL`/`env_HR`,
tensors outside the `_Chain` object entirely) that `_link_at` cannot see.
Confirmed directly: passing a window straight through `tdvp.py`'s own
sweep functions silently mis-contracts the bra/ket legs at the window's
own edge sites instead of raising.

The fix is not just a different environment seed (as `window_total_energy`
already needed): a window's own edge sites, unlike an ordinary chain's,
*always* have both a left and a right Link (`env_HL_ket`/`env_HR_ket` are
real Indices, never absent) -- so this module's own
`_window_two_site_heff`/`_window_one_site_heff` read a site's left/right
Link directly off its own tensor (`ket.A(i).inds[0]`/`.inds[-1]`, always
present, no `None`-guard needed) instead of via `_link_at`, and
`_all_left_environments_window`/`_all_right_environments_window` build
environments via idmrg.py's own explicit-Index `_extend_HL`/`_extend_HR`
(seeded from `window.env_HL`/`window.env_HR` at the two ends) rather than
`dmrg.py`'s `_link_at`-based ones. `tdvp.py`'s own `_lanczos_expm_multiply`
(a pure Krylov-propagator numerical primitive, chain-structure-agnostic)
and `svd.py`'s `svd` are reused unchanged; only the chain-structure-aware
pieces are window-specific.

`apply_local_operator`/`local_expectation` add the perturbation and
readout primitives (Sec. V.1 step 3). `local_expectation` needed its own
real fix, distinct from the TDVP one above: a plain observable's own
boundary weighting is *not* symmetric between the window's two ends,
because `U_list` is left-canonical -- see its own docstring.

== Shifted overlaps, S(x,t) (Sec. V.1 step 5, simplified to t1=0) ==

`dynamical_correlator_td`/`snapshot_correlator` reconstruct
`S(x,t)=<psi|A_x e^{-iHt}B_0|psi>` for *every* `x` from a *single* window
evolution (the paper's own headline efficiency result), by inserting
operator `A` directly into the *bra* side of the overlap at the shifted
absolute position `center+x`, rather than perturbing a second,
independently-evolved window and shifting one relative to the other
(mathematically equivalent for the `t1=0` case implemented here -- the
bra side is just the plain, un-evolved ground state, so there is nothing
to "evolve" on that side at all). `_padded_arrays` extends whichever
side's explicit tensor list falls short of `center+x` with more
(unevolved) `U_list` copies -- valid because, outside each window's own
causal cone, its tensors already equal `U_list` exactly. This is a
*simplification* of the paper's own Eq. 7 (which evolves a second window
backward by `t1` too, doubling the accessible total time for the same
TDVP cost) -- `dynamical_correlator_td`'s own docstring has the full
scope note. `_close_array_chain` is a plain-NumPy chain contraction
(mirroring idmrg.py's own `_transfer_matrices`/`_apply_transfer` style)
shared by `local_expectation` and this overlap machinery, since two
independently time-evolved chains have no shared ITensor Index
bookkeeping to exploit anyway."""

import numpy as np

from . import idmrg as _idmrg_mod
from . import kernels
from .mpscontainer import MPO, MPS
from .svd import svd
from .tdvp import _lanczos_expm_multiply
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


def _refresh_physical_legs(ket_tensors, mpo_tensors):
    """Mint a fresh, distinct physical Index for every window position,
    relabeling `ket_tensors[i]`'s single Site leg and `mpo_tensors[i]`'s
    matching unprimed/primed pair together (so ket and MPO still agree,
    site by site).

    `_tile_periodic` deliberately reuses the *same* physical Index across
    every copy of a given sublattice position (matching idmrg.py's own
    convention, e.g. `grow_by_mpo`'s own requirement that `W_bulk[p]`'s
    physical Index literally *be* `U_list[p]`'s own) -- safe there because
    nothing in idmrg.py's own machinery ever multiplies two same-position
    tensors together directly (each site's physical leg only ever meets
    its own W tensor at that same site, one site at a time via
    `_extend_HL`/`_extend_HR`). This module's own window TDVP machinery
    does not have that luxury: `_window_two_site_heff`'s `Ti * Tj` directly
    multiplies two *adjacent* window sites together -- for `n_uc=1` (every
    window site shares the identical physical Index, since there is only
    one sublattice position), this silently auto-contracts what should be
    two independent sites' physical legs into one, instead of keeping them
    as the two separate legs a genuine 2-site tensor needs. Confirmed
    directly: an n_uc=1 window's very first two-site TDVP step raised
    `transpose_to: ... is not a permutation of ...` before this fix --
    exactly the collision idmrg.py's own `_fresh_physical_copy` exists to
    avoid for its analogous p_L==p_R micro-step."""
    new_ket, new_mpo = [], []
    for U, W in zip(ket_tensors, mpo_tensors):
        old_phys = next(ind for ind in U.inds if ind.hastags("Site") and ind.plev == 0)
        new_phys = old_phys.sim()
        new_phys_p = new_phys.prime(1)
        U_inds = tuple(new_phys if ind == old_phys else ind for ind in U.inds)
        W_inds = tuple(
            new_phys if (ind.hastags("Site") and ind.plev == 0) else
            (new_phys_p if (ind.hastags("Site") and ind.plev == 1) else ind)
            for ind in W.inds)
        new_ket.append(ITensor(U_inds, U.array))
        new_mpo.append(ITensor(W_inds, W.array))
    return new_ket, new_mpo


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
    ket_tensors, mpo_tensors = _refresh_physical_legs(ket_tensors, mpo_tensors)
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


# -- real-time TDVP evolution of the window (Sec. V.1 steps 3-4) ----------

def _window_link_left(chain, i):
    """The Link Index at the left of `chain.A(i)` -- always present for a
    window's own tensors (unlike an ordinary open chain, whose site 1 has
    none, see mpscontainer.py's own docstring): a window's edge sites
    always carry env_HL_ket/env_HR_ket as a genuine leg. Read directly off
    the tensor rather than via mpsalgebra.py's `_link_at` (a same-Index
    neighbor lookup that would incorrectly return None at the window's own
    edges, see this module's own docstring)."""
    return chain.A(i).inds[0]


def _window_link_right(chain, i):
    """Mirror of `_window_link_left`."""
    return chain.A(i).inds[-1]


def _window_two_site_heff(L, Lbra, H, ket, i, R, Rbra):
    """Window analogue of dmrg.py's `two_site_heff`: the 2-site effective
    Hamiltonian at bond (i,i+1), built via the same shared
    `kernels.make_matvec` -- except the outer bond legs are always
    included (never `None`-guarded), since a window's own edge sites
    always carry a real Link (see this module's own docstring)."""
    Ti, Tj = ket.A(i), ket.A(i + 1)
    left_link = _window_link_left(ket, i)
    right_link = _window_link_right(ket, i + 1)
    s_i = next(ind for ind in Ti.inds if ind.hastags("Site"))
    s_j = next(ind for ind in Tj.inds if ind.hastags("Site"))

    order_in = [left_link, s_i, s_j, right_link]
    shape = tuple(ind.dim for ind in order_in)
    x0 = (Ti * Tj).transpose_to(order_in).reshape(-1)

    s_i_out, s_j_out = s_i.prime(1), s_j.prime(1)
    order_out = [Lbra, s_i_out, s_j_out, Rbra]
    H_i, H_j = H.A(i), H.A(i + 1)
    pieces = [p for p in (L, H_i, H_j, R) if p is not None]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    return matvec, order_in, shape, x0


def _window_one_site_heff(L, Lbra, H, ket, i, R, Rbra):
    """Window analogue of dmrg.py's `one_site_heff` -- the "backward
    evolution" piece of two-site TDVP (see tdvp.py's own module docstring
    for why it's needed), see `_window_two_site_heff`'s own docstring for
    why the bond legs are never `None`-guarded here."""
    T = ket.A(i)
    left_link = _window_link_left(ket, i)
    right_link = _window_link_right(ket, i)
    s_i = next(ind for ind in T.inds if ind.hastags("Site"))

    order_in = [left_link, s_i, right_link]
    shape = tuple(ind.dim for ind in order_in)
    x0 = T.transpose_to(order_in).reshape(-1)

    s_i_out = s_i.prime(1)
    order_out = [Lbra, s_i_out, Rbra]
    H_i = H.A(i)
    pieces = [p for p in (L, H_i, R) if p is not None]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    return matvec, order_in, shape, x0


def _all_left_environments_window(window):
    """{i: (L,Lbra)} for i=0..n-1 -- mirrors dmrg.py's own
    `_all_left_environments`'s indexing convention (env[i] = environment
    through site i, env[0] = the seed before any window site is absorbed)
    but built via idmrg.py's own explicit-Index `_extend_HL` (see this
    module's own docstring for why, not `_link_at`-based `_extend_left`),
    seeded from `window.env_HL`/`window.env_HL_bra` instead of the
    ordinary open-chain `(None, None)`. Used as the "before this
    half-sweep" static reference by `_half_sweep_rl_window` -- mirrors
    tdvp.py's own `_half_sweep_rl`, which likewise only ever needs the
    *other* side's environment precomputed, building its own side
    incrementally as the sweep progresses."""
    ket, H = window.mps, window.mpo
    n = ket.length()
    env = {0: (window.env_HL, window.env_HL_bra)}
    left_ket_old = ket.A(1).inds[0]
    for i in range(1, n):
        L_prev, Lbra_prev = env[i - 1]
        U, W_p = ket.A(i), H.A(i)
        right_ket_new = U.inds[-1]
        env[i] = _idmrg_mod._extend_HL(L_prev, Lbra_prev, W_p, U, left_ket_old, right_ket_new)
        left_ket_old = right_ket_new
    return env


def _all_right_environments_window(window):
    """Mirror of `_all_left_environments_window`: {i: (R,Rbra)} for
    i=n+1..2, seeded from `window.env_HR`/`window.env_HR_bra`. Used as the
    "before this half-sweep" static reference by `_half_sweep_lr_window`."""
    ket, H = window.mps, window.mpo
    n = ket.length()
    env = {n + 1: (window.env_HR, window.env_HR_bra)}
    right_ket_old = ket.A(n).inds[-1]
    for i in range(n, 1, -1):
        R_next, Rbra_next = env[i + 1]
        V, W_p = ket.A(i), H.A(i)
        left_ket_new = V.inds[0]
        env[i] = _idmrg_mod._extend_HR(R_next, Rbra_next, W_p, V, right_ket_old, left_ket_new)
        right_ket_old = left_ket_new
    return env


def _evolve_two_site_window(L, Lbra, H, ket, i, R, Rbra, tau, niter):
    """Forward (-i*tau) two-site evolution -- window analogue of tdvp.py's
    own `_evolve_two_site`, reusing its `_lanczos_expm_multiply` Krylov
    propagator unchanged (a pure numerical primitive, chain-structure-
    agnostic) against this module's own `_window_two_site_heff`."""
    matvec, order_in, shape, x0 = _window_two_site_heff(L, Lbra, H, ket, i, R, Rbra)
    evolved = _lanczos_expm_multiply(matvec, x0, -1j * tau, niter=niter)
    return ITensor(tuple(order_in), evolved.reshape(shape))


def _evolve_one_site_window(L, Lbra, H, ket, i, R, Rbra, tau, niter):
    """Backward (+i*tau) one-site evolution -- window analogue of tdvp.py's
    own `_evolve_one_site`, see `_evolve_two_site_window`'s own docstring."""
    matvec, order_in, shape, x0 = _window_one_site_heff(L, Lbra, H, ket, i, R, Rbra)
    evolved = _lanczos_expm_multiply(matvec, x0, 1j * tau, niter=niter)
    return ITensor(tuple(order_in), evolved.reshape(shape))


def _half_sweep_lr_window(window, tau, cutoff, maxdim, niter):
    """Window analogue of tdvp.py's own `_half_sweep_lr` -- identical
    algorithmic structure (forward-evolve each bond, SVD-split, backward-
    evolve the freshly-cut bond tensor before the next bond), built on
    this module's own window-aware environment/heff functions instead of
    tdvp.py's `_link_at`-based ones (see this module's own docstring for
    why those cannot be reused directly). One simplification unique to the
    explicit-Index approach: tdvp.py's own version must write a
    not-yet-evolved placeholder to site i+1 *before* calling `_extend_left`
    (which infers site i's own right leg via `_link_at`, needing site i+1
    to already carry a matching Link) -- `_extend_HL` here takes that leg
    as an explicit argument instead, so no placeholder write is needed.
    Mutates window.mps in place."""
    ket, H = window.mps, window.mpo
    n = ket.length()
    right_env = _all_right_environments_window(window)  # sites i+1..n+1, ket BEFORE this half-sweep
    left_env = {0: (window.env_HL, window.env_HL_bra)}
    left_ket_old = ket.A(1).inds[0]
    for i in range(1, n):
        L, Lbra = left_env[i - 1]
        R2, R2bra = right_env[i + 2]
        theta = _evolve_two_site_window(L, Lbra, H, ket, i, R2, R2bra, tau, niter)

        left_link = _window_link_left(ket, i)
        s_i = next(ind for ind in ket.A(i).inds if ind.hastags("Site"))
        U, S, V, spec = svd(theta, [left_link, s_i], cutoff=cutoff, maxdim=maxdim)
        ket.set_A(i, U)
        C = S * V
        ket.set_A(i + 1, C)

        right_ket_new = U.inds[-1]
        left_env[i] = _idmrg_mod._extend_HL(L, Lbra, H.A(i), U, left_ket_old, right_ket_new)
        left_ket_old = right_ket_new

        if i < n - 1:
            Lnew, Lnewbra = left_env[i]
            R2next, R2nextbra = right_env[i + 2]
            C_evolved = _evolve_one_site_window(Lnew, Lnewbra, H, ket, i + 1, R2next, R2nextbra, tau, niter)
            ket.set_A(i + 1, C_evolved)
    ket.center = n


def _half_sweep_rl_window(window, tau, cutoff, maxdim, niter):
    """Mirror of `_half_sweep_lr_window`, sweeping right to left -- window
    analogue of tdvp.py's own `_half_sweep_rl`."""
    ket, H = window.mps, window.mpo
    n = ket.length()
    left_env = _all_left_environments_window(window)  # sites 1..i-1, ket BEFORE this half-sweep
    right_env = {n + 1: (window.env_HR, window.env_HR_bra)}
    right_ket_old = ket.A(n).inds[-1]
    for i in range(n - 1, 0, -1):
        L2, L2bra = left_env[i - 1]
        R, Rbra = right_env[i + 2]
        theta = _evolve_two_site_window(L2, L2bra, H, ket, i, R, Rbra, tau, niter)

        right_link = _window_link_right(ket, i + 1)
        s_j = next(ind for ind in ket.A(i + 1).inds if ind.hastags("Site"))
        right_of_bond = [s_j, right_link]
        left_of_bond = [ind for ind in theta.inds if ind not in right_of_bond]
        U, S, V, spec = svd(theta, left_of_bond, cutoff=cutoff, maxdim=maxdim)
        ket.set_A(i + 1, V)
        C = U * S
        ket.set_A(i, C)

        left_ket_new = V.inds[0]
        right_env[i + 1] = _idmrg_mod._extend_HR(R, Rbra, H.A(i + 1), V, right_ket_old, left_ket_new)
        right_ket_old = left_ket_new

        if i > 1:
            L2prev, L2prevbra = left_env[i - 1]
            Rnew, Rnewbra = right_env[i + 1]
            C_evolved = _evolve_one_site_window(L2prev, L2prevbra, H, ket, i, Rnew, Rnewbra, tau, niter)
            ket.set_A(i, C_evolved)
    ket.center = 1


def window_tdvp_step(window, dt, cutoff, maxdim, niter=50):
    """One real-time step exp(-i*dt*H) on a capped window via two-site
    TDVP -- mirrors tdvp.py's own `tdvp_step` (a left-to-right half-sweep
    evolving by dt/2, then a right-to-left half-sweep by another dt/2),
    built on this module's own window-aware environment/heff functions
    (see this module's own docstring for why tdvp.py's own sweep functions
    cannot be reused directly on a window). Mutates window.mps in place."""
    tau = dt / 2.0
    _half_sweep_lr_window(window, tau, cutoff, maxdim, niter)
    _half_sweep_rl_window(window, tau, cutoff, maxdim, niter)
    return window


def apply_local_operator(window, result, site, opname):
    """Apply the named single-site operator (`result.sites_uc`'s own
    per-sublattice-type matrix convention, e.g. "Sz"/"Sx"/"Cdag" -- the
    same names `idmrg.py`'s own `_term_site_matrices` resolves via
    `site_type.matrix(name)`) to window site `site` (1-based), in place --
    Sec. V.1 step 3 of the paper (`A^dagger_0|psi>`/`B_0|psi>`).

    `mat` is applied to the ket's own physical axis with the same
    `M[in,out]` contraction convention idmrg.py's own `_op_transfer` uses
    (`Aop[l,o,r] = sum_i M[i,o] A[l,i,r]`) -- matching it matters for any
    non-symmetric operator (e.g. Cdag/C/Sp/Sm); it happens to be
    invisible for the symmetric ones (Sz, Sx, ...) this module's own
    tests exercise so far."""
    n_uc = result.n_uc
    p = (site - 1) % n_uc
    mat = result.sites_uc.site_type(p + 1).matrix(opname)
    T = window.mps.A(site)
    new_array = np.einsum('io,lir->lor', mat, T.array)
    window.mps.set_A(site, ITensor(T.inds, new_array))


def local_expectation(window, result, site, opname):
    """`<window|Op_site|window> / <window|window>`.

    *Not* a bare trace over both dangling boundary legs (env_HL_ket,
    env_HR_ket), unlike `window_total_energy`'s own energy-density
    convention: an earlier version of this function used exactly that
    (mirroring `inner(ket,ket)`'s own auto-contraction of the window's
    dangling legs, which -- since both operands share the *same* Index
    object there -- amounts to summing over *every* boundary basis state
    with equal weight) and it produced a visibly *non-uniform* profile
    even with no perturbation applied at all (confirmed directly: -0.484
    in the window's own interior, drifting to -0.10 near the right edge
    for a converged, translation-invariant ground state that should read
    exactly uniform everywhere). The bare trace is only correct on the
    *left*: idmrg.py's own `onsite_expectation`/`two_point_correlator`
    only ever close their own right side with the dominant right
    transfer-matrix fixed point (`_all_right_fixed_points`), *never* an
    analogous left one, precisely because `U_list` is left-canonical --
    the left-canonicality condition `sum_{l,s} conj(U[l,s,r']) U[l,s,r])
    = delta(r,r')` *is* the statement that the left fixed point is
    trivially the identity, needing no separate weighting. The window's
    own left boundary (env_HL_ket) inherits that same property (the
    window's tensors are literally `U_list`, away from any perturbation),
    but the *right* boundary does not -- it needs the same `rho_after`
    weighting `onsite_expectation` uses, evaluated at the window's own
    last sublattice position (always `n_uc-1`, since window site n has
    position `(n-1)%n_uc` for any n_window). Mirrors
    `_transfer_matrices`/`_apply_transfer`'s own plain-NumPy style
    directly (not idmrg.py's own functions verbatim, since those assume a
    *uniform* `U_list`, whereas the window's own tensors are generally
    non-uniform after a perturbation/time evolution)."""
    ket = window.mps
    n = ket.length()
    n_uc = result.n_uc
    p_last = (n - 1) % n_uc
    Es = _idmrg_mod._transfer_matrices(result.U_list, n_uc)
    rho_after, _eta = _idmrg_mod._all_right_fixed_points(Es, n_uc)
    rho_R = rho_after[p_last]

    p = (site - 1) % n_uc
    mat = result.sites_uc.site_type(p + 1).matrix(opname)

    arrays = [ket.A(i).array for i in range(1, n + 1)]
    op_arrays = list(arrays)
    op_arrays[site - 1] = np.einsum('io,lir->lor', mat, arrays[site - 1])
    norm = _close_array_chain(arrays, arrays, result, p_last).real
    val = _close_array_chain(arrays, op_arrays, result, p_last).real
    return val / norm


def _close_array_chain(bra_arrays, ket_arrays, result, p_right):
    """`Σ` over a chain of doubled (ket, conj(bra)) transfer steps, closed
    on the left by a bare trace (correct for a left-canonical `U_list` --
    see `local_expectation`'s own docstring) and on the right by the
    dominant right transfer-matrix fixed point at sublattice position
    `p_right` (idmrg.py's own `_all_right_fixed_points`, evaluated on the
    converged, unperturbed `result.U_list` -- the correct weighting for
    "everything beyond this chain", exactly as `onsite_expectation`/
    `two_point_correlator` already rely on).

    `bra_arrays`/`ket_arrays`: lists of `(chi_l,d,chi_r)` plain NumPy
    arrays of the same length, aligned site by site -- deliberately plain
    NumPy rather than ITensor objects, since a shifted overlap between two
    independently time-evolved windows has no shared Index bookkeeping to
    exploit anyway (every internal bond Index is freshly minted,
    independently, by each window's own TDVP sweep -- see
    `_half_sweep_lr_window`'s own docstring), so an explicit, positional
    contraction (mirroring idmrg.py's own `_transfer_matrices`/
    `_apply_transfer` style) is simpler and no less correct than trying to
    route this through ITensor's own identity-based auto-contraction."""
    E = None
    for Karr, Barr in zip(ket_arrays, bra_arrays):
        step = np.einsum('lir,LiR->lLrR', Karr, np.conj(Barr))
        E = step if E is None else np.einsum('lLrR,rRsS->lLsS', E, step)
    left_traced = np.einsum('llrR->rR', E)
    n_uc = result.n_uc
    Es = _idmrg_mod._transfer_matrices(result.U_list, n_uc)
    rho_after, _eta = _idmrg_mod._all_right_fixed_points(Es, n_uc)
    return np.einsum('rR,rR->', left_traced, rho_after[p_right % n_uc])


# -- Sec. V.1 steps 3-5: shifted overlaps, S(x,t) --------------------------

def _padded_arrays(window, result, extra_left, extra_right):
    """`window`'s own ket tensors as a plain list of `(chi_l,d,chi_r)`
    NumPy arrays, extended by `extra_left`/`extra_right` *unevolved*
    copies of `result.U_list` at each end.

    Valid (not a new approximation) because, away from a perturbation's
    own causal cone, the window's own explicit tensors already equal
    `U_list` exactly -- extending further with more `U_list` copies is
    exactly consistent with what the (conceptually infinite) window
    already represents. `extra_left` copies continue the periodic pattern
    *backward* from window site 1 (sublattice position `(-m) % n_uc` for
    the pad site `m` steps before site 1), `extra_right` continue it
    *forward* from window site n (position `(n+k) % n_uc` for the pad
    site `k` steps after site n) -- both via Python's own modulo, which
    already wraps negative operands into `[0, n_uc)` correctly. Consistent
    bond dimensions between the padding and the window's own edge sites
    are guaranteed by `build_window`'s own wraparound-dimension check
    (`U_list[-1]`'s right bond must equal `U_list[0]`'s left bond for a
    multi-copy window to exist at all)."""
    n_uc = result.n_uc
    n = window.mps.length()
    arrays = []
    for m in range(extra_left, 0, -1):
        p = (-m) % n_uc
        arrays.append(_idmrg_mod._to_array_lpr(result.U_list[p]))
    for i in range(1, n + 1):
        arrays.append(window.mps.A(i).array)
    for k in range(extra_right):
        p = (n + k) % n_uc
        arrays.append(_idmrg_mod._to_array_lpr(result.U_list[p]))
    return arrays


def snapshot_correlator(window_B, result, opname_A, x_values):
    """`{x: <ground_state| A_x |window_B>}` for every `x` in `x_values`,
    at whatever time `window_B` has already been evolved to -- this *is*
    `S(x,t) = <psi|A_x e^{-iHt}B_0|psi>` (Eq. 3-style Schrödinger-picture
    form) directly, with no extra `e^{iE0t}` phase correction needed:
    matching the rest of dmrgpy's own "TD" submode convention
    (`timedependent.py`'s `evolution_dmrg_DC`, which likewise reports this
    Schrödinger-picture matrix element as *the* correlator, with no
    Heisenberg-picture `e^{iE0t}` conversion) rather than the literal
    Heisenberg-picture `<psi|A(t)B(0)|psi>` the paper's own Eq. 3
    additionally converts to -- that conversion needs the *window's own*
    (finite, well-defined) ground-state energy as `E0`, not the
    infinite-system energy density the paper's own equation is written
    for, and only ever produces an overall, physically inert rigid shift
    of the resulting `S(k,omega)`'s own frequency axis, so it is omitted
    here for consistency with this codebase's own established convention.

    `x` is measured relative to `window_B`'s own center (the site
    `B_0` was applied to, `window_B.mps.length()//2 + 1` by the same
    convention `dynamical_correlator_td` uses to build it). The "bra" side
    is *always* the plain, unperturbed, un-evolved converged ground state
    (`result.U_list`, tiled to whatever absolute range `center+x` needs,
    via `_padded_arrays`-style padding built inline here rather than from
    an actual `IBCWindow` object, since the bra never needs a real window
    at all in this t1=0 simplification -- see this module's own docstring
    for the t1-nonzero extension this does not implement) with operator
    `opname_A` inserted at the single site `center+x`; the "ket" side is
    `window_B`'s own (generally non-uniform, evolved) explicit tensors,
    padded with unevolved `U_list` copies (`_padded_arrays`) whenever
    `center+x` falls outside `window_B`'s own explicit range."""
    n = window_B.mps.length()
    n_uc = result.n_uc
    center = n // 2 + 1
    out = {}
    for x in x_values:
        pos = center + x
        lo, hi = min(1, pos), max(n, pos)
        bra_arrays = []
        for i in range(lo, hi + 1):
            p = (i - 1) % n_uc
            arr = _idmrg_mod._to_array_lpr(result.U_list[p])
            if i == pos:
                mat = result.sites_uc.site_type(p + 1).matrix(opname_A)
                arr = np.einsum('io,lir->lor', mat, arr)
            bra_arrays.append(arr)
        ket_arrays = _padded_arrays(window_B, result, max(0, 1 - lo), max(0, hi - n))
        p_right = (hi - 1) % n_uc
        out[x] = _close_array_chain(bra_arrays, ket_arrays, result, p_right)
    return out


def dynamical_correlator_td(result, n_window, opname_A, opname_B, dt, nt,
                             cutoff, maxdim, niter=50, x_values=None):
    """Real-time dynamical correlator `S(x,t) = <psi|A_x e^{-iHt}B_0|psi>`
    for `x` in `x_values` and `t = 0, dt, 2dt, ..., (nt-1)*dt`, from a
    *single* window evolution -- Sec. V.1 of arXiv:1804.09163, simplified
    to `t1=0` (this is the naive Eq. 3 Schrödinger-picture form; see
    `snapshot_correlator`'s own docstring for why no `e^{iE0t}` conversion
    to the Heisenberg-picture `<A(t)B(0)>` form is applied, matching this
    codebase's own established "TD" submode convention) rather than the
    full two-branch trick of Eq. 7 (evolving a *second*, independent
    window backward by `t1` as well, which doubles the accessible total
    time `t1+t2` for the same number of TDVP steps on either branch) -- a
    documented, straightforward follow-up (`window_tdvp_step` already
    supports backward evolution via a negative `dt`; only the second
    branch and the corresponding overlap bookkeeping are missing), not
    attempted here.

    Even in this simplified form, the paper's own headline result already
    holds: every `x` in `x_values` comes from the *same* one window
    evolution (`window_B`), not one run per distance -- see
    `snapshot_correlator`'s own docstring for how a single run yields
    every `x`.

    `x_values` defaults to `range(-n_window*n_uc//4, n_window*n_uc//4 +
    1)` -- a conservative margin so `center+x` stays well inside
    `window_B`'s own causal-cone-limited interior for the `nt*dt` total
    time simulated here; same convergence caveat as `build_window`'s own
    `n_window` (check by increasing both and confirming `S(x,t)` stops
    changing).

    Returns `(ts, xs, S)`: `ts` (length `nt`), `xs` (sorted `x_values`),
    `S` (`nt` x `len(xs)` complex array, `S[it,ix] = S(xs[ix], ts[it])`)."""
    n_uc = result.n_uc
    if x_values is None:
        margin = max(1, n_window * n_uc // 4)
        x_values = range(-margin, margin + 1)
    xs = sorted(x_values)

    window_B = build_window(result, n_window)
    center = window_B.mps.length() // 2 + 1
    apply_local_operator(window_B, result, center, opname_B)

    ts = np.array([it * dt for it in range(nt)])
    S = np.zeros((nt, len(xs)), dtype=complex)
    for it in range(nt):
        snap = snapshot_correlator(window_B, result, opname_A, xs)
        for ix, x in enumerate(xs):
            S[it, ix] = snap[x]
        if it < nt - 1:
            window_tdvp_step(window_B, dt, cutoff=cutoff, maxdim=maxdim, niter=niter)
    return ts, np.array(xs), S
