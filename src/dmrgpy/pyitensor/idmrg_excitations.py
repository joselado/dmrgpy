"""Quasiparticle / tangent-space excitation ansatz (Haegeman et al.) on top
of a converged iDMRG ground state (idmrg.py's own IDMRGResult) -- gives a
momentum-resolved dispersion E(k) for the lowest "single-mode" excitation(s)
above the ground state, from which a scalar gap (min_k E(k)) can be
extracted. See idmrg.py's own module docstring for the ground-state growing
algorithm this builds on; this module only ever *reads* a converged
IDMRGResult (U_list, sites_uc) plus the original Hamiltonian term lists, it
never touches idmrg_ground_state itself.

== Why this, not the cheaper "diagonalize the converged superblock" idea ==

An infinite, translationally-invariant chain's excitations form a momentum-
resolved continuum/band, not a single discrete state the way a finite chain
has -- so "the first excited state" is the standard single-mode/quasiparticle
ansatz: a tangent-space vector

    |Phi_k(X)> = sum_n e^{ikn} |...A A B_n A A...>

with the SAME excitation tensor B (a function of a free matrix X) inserted at
every unit-cell position n, weighted by momentum k. This is the textbook
tangent-space excitation, not an approximation reusing the growing
algorithm's own truncated environments (which would only give a finite-size-
window-like estimate of the gap, not a genuine momentum-resolved one).

== Scope ==

- n_uc in {1, 2} (matching infinitechain.py's own overall restriction):
  multi-site unit cells are handled by first GROUPING the n_uc site tensors
  (and the n_uc automaton tensors) into one effective supersite of physical
  dimension d_g = prod(d_p) -- see `_group_unit_cell` -- then running a
  single, uniform (grouped-n_uc=1-shaped) construction. Momentum k is then
  per unit cell (the usual convention for multi-site-unit-cell tangent-space
  methods).
- Every 2-site bond term must have "reach" (supersite separation) exactly 1
  after grouping -- i.e. couple only adjacent supersites. This covers every
  ordinary nearest-neighbor-per-unit-cell Hamiltonian (all existing tests/
  examples), including the "skip-site R-suffix" case (reach-1 *after*
  grouping) -- but excludes a deliberately constructed longer-range term
  (e.g. `get_operator(name, i>=n_uc, group="R")` spanning 2 supersites),
  which is detected directly from the grouped automaton's own channel
  structure (see `_check_reach_one`) and rejected with NotImplementedError,
  mirroring idmrg.py's/infinitechain.py's existing "n_uc>2" style guards.
- itensor_version="python" (pyitensor) only.
- D==1 (the converged unit cell's bond dimension) ONLY -- i.e. an
  uncorrelated/product-like ground state, checked in
  `build_excitation_environment` and rejected with NotImplementedError
  otherwise. This is a deliberate, currently-necessary restriction, not a
  stylistic one -- see "KNOWN LIMITATION" immediately below.

== KNOWN LIMITATION: D>1 (genuinely entangled ground states) is unresolved ==

For D=1 (e.g. a field-polarized paramagnet, or any product-state ground
state), this module's dispersion E(k) has been checked to match the exact
free-fermion single-magnon dispersion of the XX+field model to ~14 digits
across the entire Brillouin zone, and every individual diagram in
`_h_eff_action` independently matches a from-scratch finite-ring tensor-
network contraction (bypassing all of this module's own machinery) to
~1e-13. Three real bugs were found and fixed this way (a missing
r-weighted trace in the Source_L/e_l energy-density computation -- see
`_regularized_environments`'s own docstring; a missing final
`_cap_right(..., r)` closure in 3 of the 7 diagrams -- see
`_h_eff_action`'s own docstring; and a numerically-singular metric
requiring the reduced-basis reformulation in `_reduce_metric`).

For a genuinely entangled ground state (D>1, e.g. the transverse-field
Ising or a dimerized Heisenberg chain), the resulting dispersion comes out
anomalously flat compared to the expected/reference answer (checked
against the exact TFIM free-fermion dispersion and against a large finite
chain's own DMRG gap), despite every individual diagram still
independently checking out against the same from-scratch ring
contraction, and H_eff(k) remaining exactly Hermitian. This was not
root-caused despite extensive investigation (ruled out: the metric's own
conditioning: `_reduce_metric` already handles this; Rh/Lh's own internal
consistency: verified both via the linear solve's own residual and via
independent direct iteration of their defining recursion; iDMRG's own
convergence quality: persisted even at state_overlap~0.9999). Given this,
`build_excitation_environment` explicitly rejects D>1 with
NotImplementedError rather than silently returning a dispersion that is
known, in at least some cases, to be wrong -- a deliberate, documented
scope limit (matching this module's own established style, e.g.
`_check_reach_one`) rather than a bug to route around quietly. Revisiting
this is a real, worthwhile follow-up, but needs a fresh, careful
re-derivation (or a from-scratch cross-check against an independent
tangent-space implementation) rather than further ad hoc debugging.

A second investigation pass found a much tighter, analytically-controlled
reproducer than the TFIM/dimerized-chain checks above, and used it to
narrow (not yet close) the search space:

- A fully-dimerized n_uc=2 Heisenberg chain (J_weak=0, i.e. decoupled
  singlet dimers) has an EXACT single-triplon excitation: dispersionless,
  E(k)=J_strong for every k (promoting one isolated dimer's singlet to a
  triplet costs exactly J_strong, independent of which dimer). This
  converges to D=1 after grouping (the inter-cell bond is exactly a
  product cut, zero Schmidt rank above 1), so it exercises this module's
  n_uc=2 grouping path while remaining inside the already-supported D=1
  scope -- and the code matches it to ~1e-16, confirming grouping itself
  is not the problem.
- Turning on a SMALL J_weak (perturbing away from that exactly solvable
  point) is exactly solvable too, to first order in degenerate
  perturbation theory: the triplon acquires nearest-dimer hopping
  amplitude -J_weak/4 (derived by hand from the explicit 4-spin matrix
  element of S_{b_n}.S_{a_{n+1}} between |t_n>|s_{n+1}> and
  |s_n>|t_{n+1}>), giving E(k) = J_strong - (J_weak/2)*cos(k). The
  *moment* J_weak becomes nonzero (D jumps straight from 1 to 4, not 1 to
  2 -- plausibly an SU(2)-multiplet effect, the leading correction mixing
  in a full triplet's worth of virtual admixture), the code's own k-
  dependence collapses to a few times 1e-4-1e-6 in absolute spread versus
  the ~J_weak (e.g. 0.02) spread this formula predicts -- a 2-4 order-of-
  magnitude suppression that reproduces even at J_weak=0.02, far inside
  where 1st-order perturbation theory should be essentially exact. Reran
  independently 3x (iDMRG's starting MPS is unseeded) and got the same
  eigenvalues/energies to 5-6 significant figures each time, ruling out
  "insufficiently converged/noisy U_list" as the explanation -- whatever
  is happening is a deterministic function of the (well-converged)
  environment, not iDMRG randomness.
- Splitting `_h_eff_action`'s diagram sum into the k-independent piece
  (diagrams 1/4/5) versus the momentum-carrying piece (6a/6b) on a random
  tangent vector shows the ratio between them collapsing by ~2 orders of
  magnitude between the D=1 XX case (comparable magnitude, ~0.18 vs 0.55)
  and any D>1 case tried (6a/6b suppressed to ~1% of the static piece or
  less) -- diagrams 2/3 (the *other* bond-touching-B terms, also k-
  independent) stayed comparably-sized to 6a/6b's *expected* scale in the
  dimer case, which if anything sharpens the puzzle: 2/3 and 6a/6b share
  the same mat_a/mat_b content and same order in J_weak, yet only the
  k-carrying pair collapses.
- Two new, so-far-unexplained numerical clues from the dimer case worth
  chasing first in any future session: (1) the dominant right fixed point
  r's eigenvalue spectrum comes out extremely skewed even at tiny J_weak
  (e.g. [6.3e-6, 6.3e-6, 6.3e-6, 0.9998] at J_weak=0.02) with the small
  eigenvalues sitting only ~6x above `_reduce_metric`'s default
  `rel_floor=1e-6` pruning threshold -- a numerically marginal regime the
  earlier TFIM check (eigenvalues either O(1) or exactly 0) never
  exercised; (2) Lh and Rh come out wildly asymmetric in norm (~1.7 vs
  ~2e-5 at J_weak=0.02) despite passing their own internal fixed-point
  residual check individually -- consistent with the *defining equation*
  itself being incomplete/asymmetric rather than a solve-time bug (the
  residual check can only catch the latter). The most promising concrete
  next step identified but not yet attempted: `IDMRGResult.env_HL`/
  `env_HR` (the growing algorithm's own converged MPO-environment
  tensors, already reused verbatim by `idmrg_window.py` -- see
  `IDMRGResult`'s own docstring) are an independent, already-correct (the
  energies/window dynamics built on them match ED) construction of
  exactly the same "background Hamiltonian environment" this module
  builds from scratch via `_regularized_environments` -- cross-checking
  Rh/Lh's F-channel content against theirs (after extracting the right
  automaton-channel slice and subtracting the same energy-density
  regularization) would show directly whether `_regularized_environments`
  is missing/misweighting a term, without having to re-derive the
  tangent-space formalism from scratch.

== Algorithm summary ==

1. Group the unit cell into one supersite (A, W) -- plain NumPy arrays, W
   still the S/F/pending automaton idmrg._build_periodic_mpo builds (see its
   own docstring), just re-expressed per supersite.
2. Null-space gauge fixing: V_L (isometry, V_L^dagger @ A_mat = 0,
   V_L^dagger @ V_L = I) via scipy.linalg.null_space -- the excitation
   tensor is B(X) = reshape(V_L @ X).
3. l = I (exact, A is left-canonical by construction -- grouping preserves
   this, see _group_unit_cell's own comment) and r = the ordinary dominant
   right fixed point of the transfer map (idmrg._dominant_right_fixed_point,
   reused directly).
4. Regularized environments Lh, Rh ("Hamiltonian to the left/right of a cut,
   energy density subtracted off so the sum is finite") -- solved as dense
   (D^2,D^2) linear systems, see `_regularized_environments`.
5. The momentum-dependent effective Hamiltonian H_eff(k)[X] -- built from a
   FINITE list of diagrams (no infinite geometric/momentum sum is actually
   needed here: the gauge condition V_L^dagger @ A = 0 makes every diagram
   where the ket differs from the bra by more than one adjacent unit cell
   vanish identically, since this module's Hamiltonians only have reach-1
   bonds -- see `_h_eff_action`'s own comment for the diagram list and the
   derivation of why only n_uc-separations 0 and +-1 ever contribute).
6. Generalized eigenproblem H_eff(k)[X] = E(k)*(X @ r) -- solved as a dense
   generalized eigenproblem at the (small) scale D*(d_g-1) x D.

This file's own numeric conventions ((ket,bra) index ordering, M[i,o]
operator-insertion convention) are deliberately built by reusing idmrg.py's
own already-validated primitives (_apply_transfer, _apply_transfer_from_left,
_op_transfer's M[i,o] convention, _dominant_right_fixed_point) rather than
re-deriving index conventions from scratch -- see the module-level comments
at each reuse site.
"""

import numpy as np
from scipy.linalg import null_space as _null_space

from . import idmrg

# Automaton channel positions, matching idmrg._build_periodic_mpo's own
# fixed convention (chans[p][0]="S", chans[p][1]="F") -- preserved exactly
# by grouping, since W_bulk[0]'s own left index and W_bulk[n_uc-1]'s own
# right index (both indexed by chans[0]) are the grouped W's own left/right
# axes verbatim (see _group_unit_cell).
_S_IDX = 0
_F_IDX = 1


def _group_unit_cell(U_list, W_bulk, n_uc):
    """Contract the n_uc per-sublattice left-canonical ket tensors U_list
    and automaton tensors W_bulk (both idmrg.py objects, from a converged
    IDMRGResult and idmrg._build_automaton respectively) into one effective
    supersite: A (D, d_g, D) and W (Dw, Dw, d_g, d_g), d_g = prod(d_p).

    A is exactly the matrix product of the n_uc site tensors -- still
    exactly left-canonical (sum_s A^s-dagger A^s = I_D): grouping two
    isometries U_0 (D,d0,chi), U_1 (chi,d1,D) preserves the isometry
    property (sum_{a,s0} conj(U_0[a,s0,m']) U_0[a,s0,m] = delta_{m,m'} lets
    the m-sum collapse first, leaving exactly U_1's own isometry condition
    summed over its own left index too -- verified by hand, generalizes to
    any n_uc by induction), so l=I_D exactly for the grouped chain too, no
    separate computation needed.

    W's channel structure (S=start, F=accumulator, pending bonds) is
    unchanged by grouping -- only the physical legs get Kronecker-combined,
    the channel (Link) legs are ordinarily contracted one bond at a time,
    with the *outer* (uncontracted) left/right channel axes of the grouped
    W being exactly W_bulk[0]'s own left axis and W_bulk[n_uc-1]'s own right
    axis (both indexed by the *same* channel list chans[0], since
    idmrg._build_periodic_mpo shares that Index across the whole periodic
    automaton) -- so _S_IDX/_F_IDX keep meaning what they mean for a single
    bulk tensor, regardless of n_uc."""
    A = idmrg._to_array_lpr(U_list[0])
    for p in range(1, n_uc):
        Ap = idmrg._to_array_lpr(U_list[p])
        A = np.einsum('lpr,rqs->lpqs', A, Ap)
        Dl, d0, d1, Dr = A.shape
        A = A.reshape(Dl, d0 * d1, Dr)

    W = W_bulk[0].array  # (Lchan, phys_in, phys_out, Mchan)
    for p in range(1, n_uc):
        Wp = W_bulk[p].array  # (Mchan, phys_in, phys_out, Rchan)
        W = np.einsum('Labm,mcdR->LacbdR', W, Wp)
        L, di0, di1, do0, do1, R = W.shape
        W = W.reshape(L, di0 * di1, do0 * do1, R)
    return A, W


def _check_reach_one(W):
    """Raise NotImplementedError if the grouped automaton W has any nonzero
    entry directly connecting two distinct "pending" channels (anything
    other than _S_IDX/_F_IDX) -- the signature of a bond term with reach>1
    (see idmrg._build_periodic_mpo's own "propagates one more site" branch,
    only reachable when a term's reach exceeds 1 unit cell). Checked
    directly on the grouped automaton's own array rather than by
    re-inspecting the original Hamiltonian term lists, so it also catches a
    deliberately constructed longer-range term (e.g.
    get_operator(name, i>=n_uc, group="R")) that idmrg.py's own automaton
    builder happily accepts but this module's diagrams (see
    _h_eff_action) do not support."""
    Dw = W.shape[0]
    pending = [c for c in range(Dw) if c not in (_S_IDX, _F_IDX)]
    for p in pending:
        for q in pending:
            if not np.allclose(W[p, :, :, q], 0):
                raise NotImplementedError(
                    "idmrg_excitations: a Hamiltonian term with reach>1 unit "
                    "cell (a bond spanning more than 2 adjacent supersites) "
                    "was detected -- the tangent-space excitation ansatz "
                    "implemented here only supports nearest-adjacent-unit-"
                    "cell (reach<=1) couplings.")


def _pending_channels(W):
    """[(mat_a, mat_b), ...] -- one pair per pending automaton channel,
    mat_a = the S->pending transition matrix (applied at the bond's own
    "earlier"/rel_a site), mat_b = the pending->F transition matrix (applied
    at the "later"/rel_b site) -- both read directly off the grouped
    automaton's array (W[chan_l, s_in, s_out, chan_r], idmrg._build_periodic_mpo's
    own axis convention), already in the M[i,o] convention idmrg._op_transfer
    uses (i=phys_in, o=phys_out)."""
    Dw = W.shape[0]
    out = []
    for p in range(Dw):
        if p in (_S_IDX, _F_IDX):
            continue
        out.append((W[_S_IDX, :, :, p], W[p, :, :, _F_IDX]))
    return out


def _onsite_matrix(W):
    """The direct S->F transition matrix -- this supersite's own onsite
    Hamiltonian content (M[i,o] convention), all-zero if there is none."""
    return W[_S_IDX, :, :, _F_IDX]


def _null_space_left(A):
    """V_L: (D*d_g, D*(d_g-1)) isometry spanning the null space of A_mat's
    own adjoint (A reshaped to (D*d_g, D)) -- i.e. V_L^dagger @ A_mat = 0,
    V_L^dagger @ V_L = I -- the standard left-gauge-fixing condition for the
    tangent-space excitation tensor B(X) = reshape(V_L @ X, A.shape)."""
    D, d_g, Dr = A.shape
    A_mat = A.reshape(D * d_g, Dr)
    V_L = _null_space(A_mat.conj().T)
    expected = D * d_g - Dr
    if V_L.shape[1] != expected:
        raise RuntimeError(
            "idmrg_excitations: the ground-state tensor's null space has "
            "dimension {}, expected {} (D*d_g - D) -- the left-canonical "
            "tensor is not full rank (a degenerate/non-injective ground "
            "state), which this module does not support.".format(
                V_L.shape[1], expected))
    return V_L


def _op_transfer_matrix(ket, bra, M=None):
    """E4[l,L,r,R] for an explicit ket/bra pair (D,d,D arrays, not
    necessarily the same tensor -- generalizes idmrg._op_transfer, which
    always used the same U_list[p] for both), with an explicit dense
    operator matrix M (M[i,o] convention, applied to the ket's physical leg
    only -- same convention as idmrg._op_transfer's own
    `einsum('io,lir->lor', M, A)`) -- M=None means plain identity/no
    operator. Reusing this exact convention (rather than a fresh derivation)
    is what lets every downstream contraction here feed straight into
    idmrg._apply_transfer/_apply_transfer_from_left unchanged."""
    if M is not None:
        ket = np.einsum('io,aic->aoc', M, ket)
    return np.einsum('lpr,LpR->lLrR', ket, np.conj(bra))


def _apply_op_ket(M, T):
    """T (D,d,D) with operator M (M[i,o] convention) applied to its
    physical leg -- M=None returns T unchanged. The single-tensor half of
    `_op_transfer_matrix` (no bra/transfer-tensor contraction), used to
    build B/A with an operator inserted while keeping bond legs open, for
    the diagrams in `_h_eff_action`."""
    if M is None:
        return T
    return np.einsum('io,aic->aoc', M, T)


def _cap_right(T, R):
    """T (D,d,D) with its right bond contracted against a (D,D) matrix R
    (ket-type contraction, i.e. R's *first* index matches T's right bond --
    the same (ket,bra) index convention idmrg._apply_transfer's own `rho`
    argument uses) -- used to close T's right bond against a fixed
    point/environment (r, R_h) while keeping T's own (left, phys) legs
    open."""
    return np.einsum('aoc,cb->aob', T, R)


def _cap_left(L, T):
    """T's left bond contracted against a (D,D) matrix L: L's *first*
    index matches T's left bond (ket-type), L's *second* index becomes the
    new left index -- confirmed by direct derivation (not just mirrored
    blindly from `_cap_right`): expanding
    tr[apply_transfer_from_left(E4(ket=A-with-mat_a, bra=A), l)] against
    idmrg.onsite_expectation's own established one-site formula the same
    way `_regularized_environments`' own e_l fix was found (see that
    function's docstring) shows the (ket,bra) output of
    idmrg._apply_transfer_from_left has its *ket*-type index first -- i.e.
    contracting L's first (ket) index against T's own (ket-type) left bond,
    NOT its second, is what correctly threads L (Lh, or a one-step
    apply_transfer_from_left result) into T's left bond. An earlier version
    of this function contracted L's second index instead -- a genuine
    transpose bug, caught by comparing the resulting dispersion against
    the exact TFIM/XX free-fermion dispersion (see the design plan's
    validation section) rather than assumed correct by symmetry with
    `_cap_right`, which is NOT actually L/R-symmetric here since L and T
    play structurally different roles."""
    return np.einsum('ba,boc->aoc', L, T)


def _dense_linear_map(D, action):
    """(D*D, D*D) dense matrix representing `action` (a (D,D) complex
    array -> (D,D) complex array linear map), built by applying `action` to
    each elementary basis matrix -- same "matvec, one basis vector at a
    time" style as this module's own H_eff construction (see
    `_build_H_eff_dense`), at the same (D^2,D^2) scale idmrg.py's own
    `_dominant_right_fixed_point` already uses for the ordinary transfer-
    matrix eigenproblem."""
    n = D * D
    mat = np.zeros((n, n), dtype=complex)
    basis = np.zeros((D, D), dtype=complex)
    for j in range(n):
        basis.flat[j] = 1.0
        mat[:, j] = action(basis).reshape(-1)
        basis.flat[j] = 0.0
    return mat


def _reduce_metric(r, rel_floor=1e-6):
    """r's eigendecomposition, restricted to eigenvalues above `rel_floor`
    times the largest one: (V, sqrt_eig, isqrt_eig), V (D,k) the kept
    eigenvectors, sqrt_eig/isqrt_eig (k,) their square roots/inverse square
    roots.

    The tangent-space norm N(X',X) = tr(X'^dagger @ X @ r) makes r the
    metric of the generalized eigenproblem H_eff(k)[X] = E*(X@r)
    (`excitation_energies`) -- solving that directly is numerically
    unusable whenever r is ill-conditioned, which is generic, not a corner
    case: r's eigenvalues are exactly the ground state's own entanglement
    spectrum (squared Schmidt values) across the cut, and any weakly-
    entangled (e.g. gapped, deep in a paramagnetic-like phase) ground state
    has a steeply decaying spectrum -- confirmed directly on a transverse-
    field Ising ground state (D=5), whose r came out with eigenvalues
    [0, 0, 0, 0.0073, 0.9927] (condition number ~1e10); solving the raw
    generalized eigenproblem there gave a dispersion that was flat and
    wrong by an O(1) amount at every momentum, not just imprecise.

    The fix (standard for a possibly-singular metric): substitute
    X = X_tilde @ diag(isqrt_eig) @ V^dagger (X_tilde: Dx x k) -- for X in
    the *dropped* subspace, |Phi_k(X)> has exactly zero norm (an
    unphysical, ill-defined direction, not merely a numerically awkward
    one), so restricting to the kept subspace discards nothing physical.
    In terms of X_tilde the norm becomes the plain Frobenius inner product
    (verified algebraically: tr(X'^dagger X r) = tr(X_tilde'^dagger
    X_tilde) after this substitution, using r = V @ diag(eig) @ V^dagger
    and V^dagger V = I on the kept subspace), turning the generalized
    eigenproblem into an ordinary Hermitian one, at the same time as
    shrinking its dimension from D to k (a bonus -- k is often much
    smaller than D, e.g. k=2 vs D=5 in the case above)."""
    herm = (r + r.conj().T) / 2
    evals, evecs = np.linalg.eigh(herm)
    floor = rel_floor * evals[-1] if evals.size and evals[-1] > 0 else 0.0
    keep = evals > max(floor, 0.0)
    V = evecs[:, keep]
    sqrt_eig = np.sqrt(evals[keep])
    isqrt_eig = 1.0 / sqrt_eig
    return V, sqrt_eig, isqrt_eig


class ExcitationEnvironment:
    """Everything the tangent-space excitation ansatz needs that does NOT
    depend on momentum k -- built once per converged IDMRGResult (expensive:
    a null-space computation plus two dense (D^2,D^2) linear solves) and
    reused across every `excitation_energies(k)` call, e.g. when scanning a
    dispersion or computing `excitation_gap`.

    `e_cell` is the unit cell's own energy density (see
    `_regularized_environments`'s own docstring for why this is computed
    self-consistently from tr(Source_R)/tr(Source_L) rather than taken from
    IDMRGResult.e0 directly) -- `excitation_energies` subtracts it from the
    raw H_eff(k) Rayleigh quotient. This subtraction is required, not
    optional: Lh/Rh (diagrams 4/5 in `_h_eff_action`) only regularize the
    Hamiltonian content strictly AWAY from the excitation tensor B itself,
    so B's own unit cell is never included in any energy-density
    subtraction anywhere else -- confirmed directly on an exactly solvable
    toy model (a pure onsite field, no bonds, D=1 product ground state):
    without this shift, a single-spin-flip momentum eigenstate came out at
    energy `h` instead of the exact `2h`, off by precisely one missed
    `e_cell` subtraction; subtracting it reproduces `2h` exactly.

    `r_V`/`r_isqrt` factor the tangent-space norm's own metric (r) onto its
    well-conditioned subspace -- see `_reduce_metric`'s own docstring for
    why solving the raw generalized eigenproblem H_eff(k)[X] = E*(X@r)
    directly is numerically unusable in general (r's eigenvalues track the
    ground state's own entanglement spectrum, which is often extremely
    skewed for a gapped/weakly-entangled state -- confirmed directly on a
    transverse-field Ising ground state, D=5, r's eigenvalues coming out as
    [0, 0, 0, 0.0073, 0.9927], a ~1e10 condition number)."""

    def __init__(self, A, W, D, d_g, V_L, l, r, Lh, Rh, e_cell):
        self.A = A
        self.W = W
        self.D = D
        self.d_g = d_g
        self.V_L = V_L
        self.l = l
        self.r = r
        self.Lh = Lh
        self.Rh = Rh
        self.e_cell = e_cell
        self.pending = _pending_channels(W)
        self.h1 = _onsite_matrix(W)
        self.r_V, self.r_sqrt, self.r_isqrt = _reduce_metric(r)


def _regularized_environments(A, W, r, l):
    """Lh, Rh (D,D): the regularized ("energy density subtracted off, so
    the semi-infinite sum is finite") left/right Hamiltonian environments --
    see this module's own docstring, point 4.

    Rh solves (I - T_A + r-outer-l)[Rh] = Source_R - e_cell*r, where
    T_A[X] = idmrg._apply_transfer(E_Id, X) (E_Id = the plain, no-operator
    transfer tensor of the grouped supersite) and
    Source_R = T_h1[r] + sum_bonds T_mat_a[T_mat_b[r]] (mat_a/mat_b: the
    S->pending/pending->F transition matrices, see `_pending_channels`) --
    i.e. exactly the one-site-onsite plus one-bond-away contributions a
    reach<=1 automaton can produce in a single step. The r-outer-l
    projector (P[X] = r*trace(l@X)) is the standard regularization that
    makes (I-T_A+P) invertible despite T_A's own dominant eigenvalue being
    exactly 1 along r -- see this module's docstring/the design plan for
    the derivation (tr(T_A[X])=tr(l@X) always, since l=I=sum_s A^s-dagger
    A^s, so (I-T_A) always maps into the traceless subspace).

    e_cell (the unit cell's own energy density, needed to regularize the
    divergent geometric sum) is deliberately computed HERE as tr(Source_R)
    itself (equivalently tr(Source_L), checked to agree), NOT taken from
    the IDMRGResult's own `e0` -- idmrg_ground_state's `e0` is a finite-
    difference of macro-iteration eigenvalues that is known to converge
    *before* U_list itself has fully settled into a self-consistent,
    translationally-invariant unit cell (see idmrg.py's own
    `_local_two_site_solve`/`IDMRGResult.state_overlap` docstrings, and
    tests/test_infinite_chain.py's `test_unit_cell_expectation_self_
    consistency`/`test_n_uc2_uniform_expectation_self_consistency`, which
    tolerate exactly this kind of gap -- confirmed directly here too: for a
    small test TFIM case, tr(Source_R) matched the *actual* <H_uc>/n_uc
    computed via idmrg.onsite_expectation/two_point_correlator to
    ~1e-14, while `result.e0` itself was off from that same actual value by
    ~4e-3). Using tr(Source_R) instead keeps Lh/Rh exactly, internally
    self-consistent with the U_list this whole module actually operates on,
    independent of that separate, pre-existing iDMRG convergence nuance.

    Lh mirrors this exactly via idmrg._apply_transfer_from_left instead
    (T_A^L[X] = apply_transfer_from_left(E_Id, X), P^L[X] = l*trace(r@X),
    Source_L = one-site onsite plus one-bond-away contributions built by
    *appending* to the right of l instead of prepending to the left of r)."""
    D = A.shape[0]
    pending = _pending_channels(W)
    h1 = _onsite_matrix(W)
    E_id = _op_transfer_matrix(A, A, None)

    def E_op(M):
        return _op_transfer_matrix(A, A, M)

    source_r = idmrg._apply_transfer(E_op(h1), r)
    for mat_a, mat_b in pending:
        inner = idmrg._apply_transfer(E_op(mat_b), r)
        source_r = source_r + idmrg._apply_transfer(E_op(mat_a), inner)

    source_l = idmrg._apply_transfer_from_left(E_op(h1), l)
    for mat_a, mat_b in pending:
        inner = idmrg._apply_transfer_from_left(E_op(mat_a), l)
        source_l = source_l + idmrg._apply_transfer_from_left(E_op(mat_b), inner)

    # e_r uses the plain trace (matches l=I: tr(T_A[X])=tr(X) always, see
    # this function's own docstring); e_l needs the trace closed with r
    # instead (tr(T_A^L[X]) is NOT tr(X) in general -- T_A^L's own
    # conserved functional is tied to r, not the identity, confirmed
    # directly against idmrg.onsite_expectation's own established formula:
    # closing a single apply_transfer_from_left(E_op(M), l) step with plain
    # trace reproduced a value systematically off from the true <M>, while
    # trace(r @ ...) matched to ~1e-14).
    e_r = np.trace(source_r)
    e_l = np.trace(r @ source_l)
    if abs(e_r - e_l) > 1e-6 * max(1.0, abs(e_r)):
        raise RuntimeError(
            "idmrg_excitations: tr(Source_R)={} and tr(Source_L)={} "
            "disagree -- both should equal the same unit-cell energy "
            "density, a mismatch signals a bug in the Source_R/Source_L "
            "construction.".format(e_r, e_l))
    e_cell = ((e_r + e_l) / 2).real

    def right_action(X):
        return X - idmrg._apply_transfer(E_id, X) + r * np.trace(l @ X)

    def left_action(X):
        return X - idmrg._apply_transfer_from_left(E_id, X) + l * np.trace(r @ X)

    rhs_r = (source_r - e_cell * r).reshape(-1)
    Mat_r = _dense_linear_map(D, right_action)
    Rh = np.linalg.solve(Mat_r, rhs_r).reshape(D, D)

    rhs_l = (source_l - e_cell * l).reshape(-1)
    Mat_l = _dense_linear_map(D, left_action)
    Lh = np.linalg.solve(Mat_l, rhs_l).reshape(D, D)

    # Internal consistency check (see this function's own docstring: the
    # regularized solution should satisfy the *un-regularized* fixed-point
    # equation exactly, since Source_R/Source_L are traceless after the
    # e_cell subtraction) -- catches an implementation bug in the dense
    # linear-map assembly, not a proof of the underlying physics, but a
    # real, load-bearing sanity check.
    resid_r = Rh - idmrg._apply_transfer(E_id, Rh) - (source_r - e_cell * r)
    if np.max(np.abs(resid_r)) > 1e-6 * max(1.0, np.max(np.abs(source_r))):
        raise RuntimeError(
            "idmrg_excitations: regularized right-environment fixed-point "
            "equation not satisfied (max residual {}) -- internal "
            "inconsistency in the environment construction.".format(
                np.max(np.abs(resid_r))))
    resid_l = Lh - idmrg._apply_transfer_from_left(E_id, Lh) - (source_l - e_cell * l)
    if np.max(np.abs(resid_l)) > 1e-6 * max(1.0, np.max(np.abs(source_l))):
        raise RuntimeError(
            "idmrg_excitations: regularized left-environment fixed-point "
            "equation not satisfied (max residual {}) -- internal "
            "inconsistency in the environment construction.".format(
                np.max(np.abs(resid_l))))
    return Lh, Rh, e_cell


def build_excitation_environment(result, h_intra_op, h_inter_op, n_uc, site_types):
    """Build an ExcitationEnvironment from a converged IDMRGResult
    (idmrg.idmrg_ground_state's own return value), the original Hamiltonian
    term lists (h_intra_op/h_inter_op, plain MultiOperator.op-format lists,
    exactly as idmrg_ground_state itself takes them -- infinitechain.py
    already stores these as self._h_intra.op/self._h_inter.op), and the raw
    per-site type codes (infinitechain.py's own self.site_types -- needed to
    rebuild the automaton via idmrg._build_automaton, which mints its own
    fresh SiteX from type codes rather than accepting result.sites_uc's
    already-built SiteType objects directly, exactly mirroring how
    idmrg_ground_state itself is always called)."""
    _, W_bulk = idmrg._build_automaton(h_intra_op, h_inter_op, site_types, n_uc)
    A, W = _group_unit_cell(result.U_list, W_bulk, n_uc)
    _check_reach_one(W)

    D, d_g, Dr = A.shape
    if Dr != D:
        raise RuntimeError(
            "idmrg_excitations: the converged unit cell's wraparound bond "
            "dimension is inconsistent (left dim {}, right dim {}) -- same "
            "failure mode idmrg._dominant_right_fixed_point already guards "
            "against, see its own comment -- try a different maxm/maxiter/"
            "etol combination for gs_energy().".format(D, Dr))
    if D != 1:
        raise NotImplementedError(
            "idmrg_excitations: the converged unit cell has bond dimension "
            "D={} (a genuinely entangled ground state) -- only D=1 "
            "(product-state-like ground states, e.g. a strongly-polarized "
            "paramagnet) is supported. This is a known, deliberate scope "
            "limit, not an oversight -- see this module's own docstring "
            "(\"KNOWN LIMITATION\" section) for what was tried and why D>1 "
            "is not enabled yet.".format(D))

    E_self = _op_transfer_matrix(A, A, None)
    r, eta = idmrg._dominant_right_fixed_point([E_self])
    r = (r + r.conj().T) / 2  # Hermitize defensively, mirrors svd._psd_sqrt_factor
    l = np.eye(D, dtype=complex)

    V_L = _null_space_left(A)
    Lh, Rh, e_cell = _regularized_environments(A, W, r, l)
    return ExcitationEnvironment(A, W, D, d_g, V_L, l, r, Lh, Rh, e_cell)


def _h_eff_action(k, X, env):
    """H_eff(k)[X] -- the momentum-dependent effective Hamiltonian acting on
    a tangent-space parameter X ((D*(d_g-1), D) matrix), returned in the
    same shape.

    Diagrams (Y accumulated in the (D,d_g,D) shape of B(X)=reshape(V_L@X),
    then projected once via V_L^dagger at the very end):

    1. Onsite term acting directly on B: cap_right(h1 applied to B, r).
    2/3. A reach-1 bond straddling B (at unit cell 0) and the *ordinary*
       ground-state tensor A at the adjacent cell (+1 or -1) -- the bond's
       "other half" is capped by one application of idmrg._apply_transfer(
       _from_left) against r/l, since beyond that single adjacent site
       everything is the plain converged background again.
    4/5. Background Hamiltonian content strictly to the right/left of B,
       with B itself untouched -- exactly Rh/Lh (built once, k-independent).
    6a/6b. The momentum-summed piece: a reach-1 bond connecting the bra's
       excitation at cell 0 to the SAME B (weighted e^{+-ik}) sitting on the
       KET side at cell +-1 instead, with an ordinary A at cell 0 on the ket
       side. This is the only place k enters.

    Every diagram's own "local, touched" part must be closed on BOTH sides
    -- the left with l, the right with r -- even where nothing acts, since
    that side still represents an infinite trivial background needing
    regularization, not simply "leave it open, it'll match B' directly".
    Diagrams 1/2/4/6a happen to end with a `_cap_right(...)` step that
    already performs this right-side closure explicitly; diagrams 3, 5 and
    6b instead end with `_cap_left(...)` (threading their own environment
    into position 0 from the left) and so need one further explicit
    `_cap_right(..., r)` -- easy to miss since omitting the *left*-side
    closure is harmless whenever it would have been closed with `l` (=
    Identity, a no-op, which is why leaving diagrams 1/2/4/6a's own left
    side implicit -- relying on the final V_L^dagger projection to match
    B' there directly -- is fine) but omitting the *right*-side closure
    with `r` (not the identity in general) is not a no-op at all. Found by
    direct numerical comparison against an exact, from-scratch finite-ring
    tensor-network contraction (bypassing all of this module's own
    machinery) after the raw (unfixed) diagrams 3/5/6b failed both that
    ring comparison and the required H_eff(k) Hermiticity identity -- see
    the design plan for the debugging trail; not something derivable from
    the diagrams' own shape alone.

    No genuine infinite/geometric momentum sum is needed for cell
    separations |n|>=2: this module's Hamiltonians have reach<=1 (checked by
    `_check_reach_one`), and the excitation tensor's own gauge condition
    (V_L^dagger @ A = 0) makes any diagram where the bra has B at cell 0 but
    the ket has a *plain* A there (no operator acting at cell 0 at all)
    vanish identically, regardless of what happens further away -- so a
    Hamiltonian term can only connect the bra's fixed excitation (cell 0) to
    a ket excitation at cell n!=0 by directly touching cell 0 itself, which
    (reach<=1) is only possible for n in {-1,0,+1}. This was confirmed by
    hand (not just assumed) before implementing -- see the design plan for
    the derivation."""
    D, d_g = env.D, env.d_g
    A, r, l = env.A, env.r, env.l
    Lh, Rh, h1 = env.Lh, env.Rh, env.h1

    B = (env.V_L @ X).reshape(D, d_g, D)

    Y = _cap_right(_apply_op_ket(h1, B), r)                    # diagram 1
    Y = Y + _cap_right(B, Rh)                                    # diagram 4
    Y = Y + _cap_right(_cap_left(Lh, B), r)                      # diagram 5

    phase_p, phase_m = np.exp(1j * k), np.exp(-1j * k)
    for mat_a, mat_b in env.pending:
        cap_b_r = idmrg._apply_transfer(_op_transfer_matrix(A, A, mat_b), r)
        Y = Y + _cap_right(_apply_op_ket(mat_a, B), cap_b_r)          # diagram 2

        cap_a_l = idmrg._apply_transfer_from_left(_op_transfer_matrix(A, A, mat_a), l)
        Y = Y + _cap_right(_cap_left(cap_a_l, _apply_op_ket(mat_b, B)), r)   # diagram 3

        cap_bB_r = idmrg._apply_transfer(_op_transfer_matrix(B, A, mat_b), r)
        Y = Y + phase_p * _cap_right(_apply_op_ket(mat_a, A), cap_bB_r)   # diagram 6a

        cap_aB_l = idmrg._apply_transfer_from_left(_op_transfer_matrix(B, A, mat_a), l)
        Y = Y + phase_m * _cap_right(
            _cap_left(cap_aB_l, _apply_op_ket(mat_b, A)), r)              # diagram 6b

    Y_mat = Y.reshape(D * d_g, D)
    return env.V_L.conj().T @ Y_mat


def _build_H_eff_dense(k, env):
    """Dense (Dx*k, Dx*k) matrix representing H_eff(k) in the *reduced*,
    well-conditioned tangent-space basis (see `_reduce_metric`'s own
    docstring for why the raw (Dx*D, Dx*D) basis is numerically unusable
    in general), Dx=D*(d_g-1), k = number of kept eigenvalues of r
    (env.r_V/env.r_sqrt/env.r_isqrt) -- built one basis vector at a time,
    same style as `_dense_linear_map`.

    For each reduced basis vector X_tilde (Dx,k): reconstruct the full
    X = X_tilde @ diag(r_isqrt) @ r_V^dagger (Dx,D), evaluate the existing,
    validated `_h_eff_action(k, X, env)` unchanged, then project the
    result back down the same way: Y_tilde = Y @ r_V @ diag(r_isqrt) --
    exactly the substitution `_reduce_metric` derives, applied to H_eff
    itself rather than just the norm."""
    D, d_g = env.D, env.d_g
    Dx = D * (d_g - 1)
    V, isqrt = env.r_V, env.r_isqrt
    kdim = V.shape[1]
    n = Dx * kdim
    H = np.zeros((n, n), dtype=complex)
    Xt = np.zeros((Dx, kdim), dtype=complex)
    for j in range(n):
        Xt.flat[j] = 1.0
        X = Xt @ (isqrt[:, None] * V.conj().T)
        Y = _h_eff_action(k, X, env)
        Yt = Y @ V @ np.diag(isqrt)
        H[:, j] = Yt.reshape(-1)
        Xt.flat[j] = 0.0
    return H


def excitation_energies(env, k, n=1):
    """The lowest `n` excitation energies (above the ground state) at
    momentum `k` (radians, per unit cell).

    The tangent-space norm <Phi_k(X')|Phi_k(X)> = tr(X'^dagger @ X @ r)
    (derived the same way as `_h_eff_action`'s own diagrams -- both bra/ket
    at the same, single insertion position, capped by l=I on the left and
    r on the right) makes the *generalized* eigenproblem
    H_eff(k)[X] = E_raw(k)*(X@r) the naive thing to solve -- but r is, in
    general, an ill-conditioned metric (see `_reduce_metric`'s own
    docstring), so this instead solves the equivalent *ordinary* Hermitian
    eigenproblem in the reduced basis `_build_H_eff_dense` builds. Finally
    subtracts `env.e_cell` from every raw eigenvalue -- see
    ExcitationEnvironment's own docstring for why this shift is required
    (B's own unit cell is never covered by the Lh/Rh background
    subtraction, so its own energy-density contribution has to be removed
    here instead)."""
    Hmat = _build_H_eff_dense(k, env)
    Hmat = (Hmat + Hmat.conj().T) / 2  # Hermitize (H is Hermitian; this is numerical noise cleanup)
    w = np.linalg.eigvalsh(Hmat)
    w = np.sort(w) - env.e_cell
    return w[:n]
