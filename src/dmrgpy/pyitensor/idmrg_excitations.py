"""Quasiparticle / tangent-space excitation ansatz (Haegeman et al.) on top
of a converged VUMPS ground state (vumps.py's own VUMPSResult) -- gives a
momentum-resolved dispersion E(k) for the lowest "single-mode" excitation(s)
above the ground state, from which a scalar gap (min_k E(k)) can be
extracted.

== Why this, not the cheaper "diagonalize the converged superblock" idea ==

An infinite, translationally-invariant chain's excitations form a momentum-
resolved continuum/band, not a single discrete state the way a finite chain
has -- so "the first excited state" is the standard single-mode/quasiparticle
ansatz: a tangent-space vector

    |Phi_k(X)> = sum_n e^{ikn} |...A_L A_L B_n(X) A_R A_R...>

with the SAME excitation tensor B (a function of a free matrix X) inserted at
every unit-cell position n, weighted by momentum k, and the ground state
written in the mixed gauge {A_L, A_R, C} (A_L left of the insertion, A_R
right of it, C the bond matrix between them in the k=0/no-insertion case).
This is the textbook mixed-gauge tangent-space excitation (Vanderstraeten,
Haegeman, Verstraete, "Tangent-space methods for uniform matrix product
states", arXiv:1810.07006, Sec. 6), not an approximation reusing a growing
algorithm's own truncated environments.

== Scope ==

- n_uc in {1, 2} (matching vumps.py's/infinitechain.py's own overall
  restriction): a multi-site unit cell is handled entirely inside
  vumps.py, which groups the n_uc sites into one effective supersite of
  physical dimension d_g=prod(d_p) before ever building AL/AR/C -- this
  module only ever sees the already-grouped, single-supersite result.
  Momentum k is per unit cell.
- Every bond term must have "reach" (supersite separation) exactly 1 after
  grouping -- checked once, inside `vumps.vumps_ground_state` itself
  (`idmrg_excitations._check_reach_one`, reused from this module), so a
  `VUMPSResult` reaching this module has already passed that check.
- itensor_version="python" (pyitensor), gs_method="vumps" only -- unlike
  the growing/infinite-size algorithm (pyitensor/idmrg.py), VUMPS solves
  directly at a fixed target bond dimension D in the thermodynamic limit
  and returns the mixed-gauge {AL,AR,C,GL,GR} representation this
  ansatz needs with no reconstruction step (see vumps.py's own module
  docstring).
- Any bond dimension D>=1 -- unlike an earlier, now-superseded version of
  this module (see "History" below), a genuinely entangled (D>1) ground
  state is fully supported.

== Algorithm (mirrors MPSKit.jl's own architecture) ==

Given AL, AR, C, GL, GR (vumps.VUMPSResult, already validated: GL/GR are
the ordinary, momentum-independent "Hamiltonian to the left/right of a cut"
environments, built the same way idmrg_excitations' own historical Lh/Rh
were, but from the mixed-gauge AL/AR rather than a single uniform-gauge
tensor) and the grouped automaton W (S/F/pending-channel finite-state
Hamiltonian, `idmrg._build_periodic_mpo`'s own construction):

1. Null-space gauge fixing: V_L (isometry, V_L^dagger @ AL_mat = 0,
   V_L^dagger @ V_L = I) via scipy.linalg.null_space -- the excitation
   tensor is B(X) = reshape(V_L @ X, AL.shape). No r-weighted metric
   reduction is needed here (unlike an earlier, uniform-gauge version of
   this module) -- the mixed-gauge tangent-space norm <Phi'|Phi> =
   tr(X'^dagger X) is trivially Euclidean, since gauge is split exactly at
   the excitation on both bra and ket (Vanderstraeten et al., Sec. 6.1).
2. GL_full/GR_full: the FULL channel-resolved (not just F/S-projected)
   background environments -- {S: I, F: GL, pending_p: Lvec_a_p} and
   {F: I, S: GR, pending_p: Rvec_b_p} respectively (Lvec_a/Rvec_b:
   vumps._precompute_bond_environments' own one-more-site bond content).
   `_full_channel_contraction` assembles "left_full * ket * H[site] *
   right_full" from these -- reusing it with (GL_full, X, GR_full) exactly
   reproduces vumps._h_ac_action(X,GL,GR,bond_envs,h1) (confirmed
   directly, to ~1e-16, before trusting it for the momentum-dependent
   pieces below).
3. GBL(k)/GBR(k): the momentum-summed, channel-resolved "the excitation
   already happened somewhere to the left/right" environments -- MPSKit's
   own `lBs`/`rBs` (src/environments/qp_envs.jl), built here as an
   explicit, self-consistent (period-1, since the unit cell is a single
   grouped supersite) linear fixed-point solve per momentum, mirroring how
   `_full_channel_contraction`'s own recursion for the *ordinary*
   background (GL_full/GR_full) already resums an infinite sum via a
   single regularized linear solve -- see `_build_GBL`/`_build_GBR`'s own
   docstrings for the exact recursion and the two subtle, previously
   error-prone points it depends on (the momentum-phase direction, and the
   E_RL/E_LR *mixed*-transfer fixed points used to regularize the k=0
   sector).
4. H_eff(k)[X] -- exactly 3 terms (MPSKit's own
   `_effective_excitation_local_apply`, `algorithms/excitation/
   quasiparticleexcitation.jl`): "B in center"
   (`_full_channel_contraction(GL_full, B, GR_full, ...)`), "B to the
   left" (`_full_channel_contraction(GBL, AR, GR_full, ...)`), "B to the
   right" (`_full_channel_contraction(GL_full, AL, GBR, ...)`) -- see
   `_h_eff_action`.
5. H_eff(k)[X] = lambda(k) * X, an ordinary (not generalized -- see point 1)
   Hermitian eigenproblem of size Dx*D (Dx=D*(d_g-1)), solved directly.
   `lambda_AC` (H_AC's own Rayleigh quotient on the converged AC=AL@C,
   *not* the ground state's physical energy density e_cell=lambda_AC-
   lambda_C -- see `ExcitationEnvironment`'s own docstring) is subtracted
   from every eigenvalue to give the excitation energy above the ground
   state.

This file's numeric conventions ((ket,bra) index ordering, M[i,o]
operator-insertion convention, the S/F/pending automaton channel layout)
are shared verbatim with vumps.py (which imports this module's primitives
as `idmrg_exc.*`) and idmrg.py (`_apply_transfer`,
`_apply_transfer_from_left`, `_dominant_right_fixed_point`,
`_dominant_left_fixed_point`).

== History ==

An earlier version of this module implemented a hand-derived, closed-form-
resummed diagram list on top of a single uniform-gauge tensor from
idmrg.py's growing algorithm, and was exact at D=1 but wrong (anomalously
suppressed dispersion) at D>1 despite eight independent investigation
passes (see `git log` on this file, and
`docs/idmrg_excitation_mpskit_port_plan.md`). That investigation's own
conclusion -- stop patching the diagram list, rewrite from scratch
mirroring MPSKit.jl's actual (not just published-equation) architecture,
on top of a genuine mixed-gauge {AL,AR,C} ground state -- is what this
version implements. Getting it working at D>1 additionally required two
further fixes, made while validating this module directly against an
independently-converged MPSKit.jl D=2 TFIM state (reproducible via
`docs/idmrg_excitation_mpskit_port_plan.md`'s own "Environment" section):
a real, pre-existing sign bug in vumps.py itself (`_solve_left_
environment`/`_solve_right_environment`/`_energy_density_and_source_from_
left`/`_right` were all missing a conjugate on the dominant fixed point
used to close an `apply_transfer_from_left`/`apply_transfer` output --
invisible at D=1, where that fixed point is always real, see vumps.py's
own docstrings at those functions for the derivation), and the
`lambda_AC`-not-`e_cell` renormalization point above. With both fixes, the
D=2 TFIM dispersion matches MPSKit's own D=2 result to 6 significant
figures at every momentum tested, and H_eff(k) is Hermitian to machine
precision, at D=1, 2 and 3.
"""

import numpy as np
from scipy.linalg import null_space as _null_space

from . import idmrg

# Automaton channel positions, matching idmrg._build_periodic_mpo's own
# fixed convention (chans[p][0]="S", chans[p][1]="F").
_S_IDX = 0
_F_IDX = 1


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
    builder happily accepts but this module's H_eff construction does not
    support. Called once, inside `vumps.vumps_ground_state`, so a
    `VUMPSResult` reaching this module's own `build_excitation_environment`
    has already passed this check."""
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
    necessarily the same tensor), with an explicit dense operator matrix M
    (M[i,o] convention, applied to the ket's physical leg only -- same
    convention as idmrg._op_transfer's own `einsum('io,lir->lor', M, A)`)
    -- M=None means plain identity/no operator. Reusing this exact
    convention (rather than a fresh derivation) is what lets every
    downstream contraction here feed straight into
    idmrg._apply_transfer/_apply_transfer_from_left unchanged."""
    if M is not None:
        ket = np.einsum('io,aic->aoc', M, ket)
    return np.einsum('lpr,LpR->lLrR', ket, np.conj(bra))


def _apply_op_ket(M, T):
    """T (D,d,D) with operator M (M[i,o] convention) applied to its
    physical leg -- M=None returns T unchanged."""
    if M is None:
        return T
    return np.einsum('io,aic->aoc', M, T)


def _cap_right(T, R):
    """T (D,d,D) with its right bond contracted against a (D,D) matrix R
    (ket-type contraction, i.e. R's *first* index matches T's right bond)
    -- used to close T's right bond against a fixed point/environment while
    keeping T's own (left, phys) legs open."""
    return np.einsum('aoc,cb->aob', T, R)


def _cap_left(L, T):
    """T's left bond contracted against a (D,D) matrix L: L's *first*
    index matches T's left bond (ket-type), L's *second* index becomes the
    new left index. NOT symmetric with `_cap_right` -- L and T play
    structurally different roles here (confirmed by direct derivation
    against idmrg.onsite_expectation's own one-site formula; see git
    history for the transpose bug this once caught)."""
    return np.einsum('ba,boc->aoc', L, T)


def _dense_linear_map(D, action):
    """(D*D, D*D) dense matrix representing `action` (a (D,D) complex
    array -> (D,D) complex array linear map), built by applying `action` to
    each elementary basis matrix -- same "matvec, one basis vector at a
    time" style used throughout this module and vumps.py."""
    n = D * D
    mat = np.zeros((n, n), dtype=complex)
    basis = np.zeros((D, D), dtype=complex)
    for j in range(n):
        basis.flat[j] = 1.0
        mat[:, j] = action(basis).reshape(-1)
        basis.flat[j] = 0.0
    return mat


def _mixed_fixed_points(E_mixed):
    """(r, l): the dominant right/left fixed points of an arbitrary mixed
    transfer tensor E_mixed (D,D,D,D), normalized so trace(l @ r) = 1 --
    reuses idmrg._dominant_right_fixed_point/_dominant_left_fixed_point
    directly (already-validated, generic eigenproblem code, not specific to
    a self-overlap transfer) rather than hand-deriving the fixed point in
    closed form (e.g. via C/C^dagger algebra) -- deliberately, since a
    previous investigation pass's hand-derived mixed-transfer fixed points
    were later found to disagree with MPSKit.jl's own ground truth in
    a subtle conjugate/transpose way (see this module's own "History"
    section) -- reusing the same numerical machinery idmrg.py already uses
    everywhere else avoids re-introducing that class of bug. Raises
    RuntimeError (via the reused functions' own degeneracy guard) if
    E_mixed's dominant eigenvalue is not well-separated, and if either
    fixed point's own eigenvalue is not close to 1 (the expected
    normalization for a converged VUMPS mixed-gauge state -- a large
    deviation signals `vumps_result` is not actually converged)."""
    D = E_mixed.shape[0]
    r, eta_r = idmrg._dominant_right_fixed_point([E_mixed])
    l, eta_l = idmrg._dominant_left_fixed_point([E_mixed])
    for name, eta in (("right", eta_r), ("left", eta_l)):
        if abs(abs(eta) - 1.0) > 1e-6:
            raise RuntimeError(
                "idmrg_excitations: the mixed transfer tensor's dominant {} "
                "eigenvalue is {} (expected magnitude 1) -- the VUMPS ground "
                "state this excitation environment was built from is not "
                "well converged.".format(name, eta))
    norm = np.trace(l @ r)
    l = l / norm
    return r, l


def _full_channel_contraction(left_full, ket, right_full, W, h1, pending):
    """sum over every nonzero automaton channel pair (c,c') of
    `left_full[c] -cap_left- (ket, with W[c,:,:,c'] applied if c!=c'-of-
    onsite-type) -cap_right- right_full[c']`, returned in (D,d_g,D) shape
    (matching `ket`'s own shape) -- the generic "left-channel-environment *
    ket * H[site] * right-channel-environment" contraction this module's
    H_eff(k) is built from three times (see this module's own docstring,
    algorithm step 4), for whichever `ket`/`left_full`/`right_full` triple
    a given term needs.

    `left_full`/`right_full`: dicts (or any `[]`-indexable-by-channel
    object) covering every channel this module's reach<=1 automaton can
    have -- {_S_IDX: ..., _F_IDX: ..., pending-channel-index: ...} for
    each pending channel (index = 2, 3, ... in `_pending_channels`' own
    enumeration order, matching W's own channel axis directly). Iterates
    the S->F (onsite h1), F->F (identity/background pass-through), S->S
    (identity/background pass-through) and S->pending_p/pending_p->F (each
    bond) transitions explicitly -- exactly the nonzero-transition set
    `_check_reach_one` already guarantees is exhaustive for any Hamiltonian
    this module accepts, so no other channel pair needs to be considered.

    With (left_full, right_full) = (GL_full, GR_full) -- the *ordinary*,
    momentum-independent background environments {S: I, F: GL,
    pending_p: Lvec_a_p} / {F: I, S: GR, pending_p: Rvec_b_p} -- and
    ket=X, this reproduces vumps._h_ac_action(X, GL, GR, bond_envs, h1)
    exactly (confirmed directly, to ~1e-16) -- the "B in center" H_eff
    term (see `_h_eff_action`) is exactly this call. The other two H_eff
    terms reuse this same helper with the momentum-dependent, channel-
    resolved GBL(k)/GBR(k) (`_build_GBL`/`_build_GBR`) standing in for
    GL_full/GR_full on whichever side the excitation has already passed."""
    S, F = _S_IDX, _F_IDX
    Y = _cap_left(left_full[S], _cap_right(_apply_op_ket(h1, ket), right_full[F]))
    Y = Y + _cap_left(left_full[F], _cap_right(ket, right_full[F]))
    Y = Y + _cap_left(left_full[S], _cap_right(ket, right_full[S]))
    for idx, (mat_a, mat_b) in enumerate(pending):
        p = idx + 2
        Y = Y + _cap_left(left_full[S],
                           _cap_right(_apply_op_ket(mat_a, ket), right_full[p]))
        Y = Y + _cap_left(left_full[p],
                           _cap_right(_apply_op_ket(mat_b, ket), right_full[F]))
    return Y


def _channel_resolvent(phase_factor, E_mixed, r_mixed, l_mixed, apply_fn, D):
    """A `solve(source) -> X` callable for (I - phase_factor*T + P)[X] =
    source, T[X] = apply_fn(E_mixed, X) (`apply_fn` is
    idmrg._apply_transfer_from_left for `_build_GBL`'s own use,
    idmrg._apply_transfer for `_build_GBR`'s) -- built once per momentum
    and mixed-transfer direction, reused across every channel this
    momentum's GBL(k)/GBR(k) solve needs (`_build_GBL`/`_build_GBR` each
    call this exactly once, for their own shared self-loop channel).

    T's own dominant eigenvalue sits exactly at 1 (E_mixed is always
    built with M=None -- the plain identity/background transfer -- for the
    two channels that ever hit this function's own self-loop, see
    `_build_GBL`/`_build_GBR`), so (I-phase_factor*T) is only singular
    exactly when phase_factor=1 (k=0 mod 2*pi for `_build_GBR`'s own
    direct `phase`, k=0 mod 2*pi for `_build_GBL`'s `1/phase` too, since
    |phase|=1 always) -- handled the same way idmrg_excitations' own
    (now-superseded) uniform-gauge resolvent handled its analogous k=0
    singularity: adding the r-outer-l-style projector P[X] = X0 *
    trace(Y0 @ X), X0 = T's own dominant *right* fixed point in the
    apply_fn sense (= `l_mixed`'s own defining fixed-point, NOT
    `r_mixed` -- `apply_fn`=apply_transfer_from_left's own forward-
    iterated fixed point is `_dominant_left_fixed_point`'s result, an easy
    mix-up since the two dominant-fixed-point functions' own "left"/
    "right" naming refers to which SIDE of the transfer matrix they solve,
    not which apply_fn direction they are the natural fixed point of), Y0
    = conj(r_mixed) (the *adjoint* map's own fixed point -- derived
    directly, and confirmed the conjugate is required by the same
    empirical D=2 cross-check that found vumps.py's own r_AL/l_AR
    conjugate bug, see this module's "History" section)."""
    at_zero = abs(phase_factor - 1.0) < 1e-10

    def action(X):
        out = X - phase_factor * apply_fn(E_mixed, X)
        if at_zero:
            out = out + l_mixed * np.trace(r_mixed.conj() @ X)
        return out

    Mat = _dense_linear_map(D, action)

    def solve(source):
        return np.linalg.solve(Mat, source.reshape(-1)).reshape(D, D)

    return solve


def _build_GBL(env, B, k):
    """GBL(k): the channel-resolved "the excitation has already happened
    somewhere to the left" environment (MPSKit's own `lBs`,
    `environments/qp_envs.jl`) -- a dict {S: (D,D), F: (D,D),
    pending_p: (D,D), ...}, solved self-consistently for this specific
    excitation tensor `B` and momentum `k` (radians per unit cell).

    Recursion (MPSKit's own, mirrored exactly -- see
    docs/idmrg_excitation_mpskit_port_plan.md's "Environment"/"Oracle run"
    sections for how this was cross-checked): for each target channel c',

        GBL[c'] = (1/phase) * [
            sum_c apply_transfer_from_left(E_RL(W[c,:,:,c']), GBL[c])
          + sum_c apply_transfer_from_left(E_B(W[c,:,:,c']), GL_full[c]) ]

    where E_RL(M) = _op_transfer_matrix(AR, AL, M) (the *mixed*, ket=AR/
    bra=AL, transfer -- "everything since the excitation is now in the
    AR-ket region, but the bra is always AL", the standard mixed-gauge
    tangent-space convention) and E_B(M) = _op_transfer_matrix(B, AL, M)
    (the excitation itself, inserted directly at this recursion step).
    Because this module's automaton only has S->S/F->F (identity pass-
    through), S->pending_p and pending_p->F nonzero transitions (any
    reach<=1 automaton, `_check_reach_one`), this decomposes into a
    solvable order: channel S (a genuine self-loop via S->S, solved via
    `_channel_resolvent`), then each pending_p (a *direct* formula from
    channel S, no self-loop, no solve needed -- only S transitions into a
    pending channel), then channel F (another self-loop, via F->F, its own
    inhomogeneous source built from the already-solved S and pending_p
    channels plus the S->F/pending_p->F transitions).

    A term easy to drop by pattern-matching too quickly from
    `_full_channel_contraction`'s own S->S/F->F self-terms: channel F's
    own source needs a `B` inserted directly into F's own identity pass-
    through slot, closed against `GL_full[F]=GL` (the *already*-
    accumulated ordinary background) -- i.e. "the excitation replaces
    exactly the site where the ordinary GL environment would otherwise
    grow by one more step". This is NOT optional/negligible: without it,
    D=1 still comes out exact (both terms vanish identically there -- V_L's
    own gauge condition, not this term's absence, is what keeps D=1
    correct either way) but D>1 fails Hermiticity by 10-30% (confirmed
    directly on the D=2 TFIM reproducer this module's own "History" section
    references); with it, Hermiticity is machine precision at every D
    tried. `_build_GBR` needs the mirror term (`GR_full[S]=GR` inserted
    into channel S's own identity self-loop) for the same reason."""
    D = env.D
    AL, AR = env.AL, env.AR
    h1, pending = env.h1, env.pending
    GL_full = env.GL_full
    phase = np.exp(1j * k)

    resolve = _channel_resolvent(1.0 / phase, env.E_RL, env.r_RL, env.l_RL,
                                  idmrg._apply_transfer_from_left, D)

    def E_RL(M):
        return _op_transfer_matrix(AR, AL, M)

    def E_B(M):
        return _op_transfer_matrix(B, AL, M)

    src_S = idmrg._apply_transfer_from_left(E_B(None), GL_full[_S_IDX])
    G_S = resolve(src_S / phase)

    G_pending = {}
    for idx, (mat_a, mat_b) in enumerate(pending):
        p = idx + 2
        term_bg = idmrg._apply_transfer_from_left(E_RL(mat_a), G_S)
        term_src = idmrg._apply_transfer_from_left(E_B(mat_a), GL_full[_S_IDX])
        G_pending[p] = (term_bg + term_src) / phase

    rhs_F = idmrg._apply_transfer_from_left(E_RL(h1), G_S)
    rhs_F = rhs_F + idmrg._apply_transfer_from_left(E_B(h1), GL_full[_S_IDX])
    rhs_F = rhs_F + idmrg._apply_transfer_from_left(E_B(None), GL_full[_F_IDX])
    for idx, (mat_a, mat_b) in enumerate(pending):
        p = idx + 2
        rhs_F = rhs_F + idmrg._apply_transfer_from_left(E_RL(mat_b), G_pending[p])
        rhs_F = rhs_F + idmrg._apply_transfer_from_left(E_B(mat_b), GL_full[p])
    G_F = resolve(rhs_F / phase)

    GBL = {_S_IDX: G_S, _F_IDX: G_F}
    GBL.update(G_pending)
    return GBL


def _build_GBR(env, B, k):
    """Mirror of `_build_GBL`: GBR(k), the channel-resolved "the excitation
    has already happened somewhere to the right" environment (MPSKit's own
    `rBs`). Recursion (note: MULTIPLIES by `phase`, `_build_GBL`'s own
    recursion DIVIDES -- confirmed against MPSKit's own source this is not
    a typo, both directions genuinely use the same `phase=exp(1j*k)`, just
    combined differently, see docs/idmrg_excitation_mpskit_port_plan.md's
    own "channel-resolved GBL/GBR recursion" section):

        GBR[c] = phase * [
            sum_c' apply_transfer(E_LR(W[c,:,:,c']), GBR[c'])
          + sum_c' apply_transfer(E_B(W[c,:,:,c']), GR_full[c']) ]

    where E_LR(M) = _op_transfer_matrix(AL, AR, M) (ket=AL/bra=AR, the
    mirror mixed transfer) and E_B(M) = _op_transfer_matrix(B, AR, M).
    Solve order mirrors `_build_GBL`'s (S<->F roles swapped, matching how
    GR_full's own channel roles are swapped relative to GL_full's -- see
    `build_excitation_environment`'s own docstring): channel F first (a
    self-loop via F->F, no incoming edges from any other channel), each
    pending_p next (a direct formula from channel F only), then channel S
    (a self-loop via S->S, with contributions from h1 (channel F),
    pending_p (via each mat_a) and its own identity self-term)."""
    D = env.D
    AL, AR = env.AL, env.AR
    h1, pending = env.h1, env.pending
    GR_full = env.GR_full
    phase = np.exp(1j * k)

    resolve = _channel_resolvent(phase, env.E_LR, env.r_LR, env.l_LR,
                                  idmrg._apply_transfer, D)

    def E_LR(M):
        return _op_transfer_matrix(AL, AR, M)

    def E_B(M):
        return _op_transfer_matrix(B, AR, M)

    src_F = idmrg._apply_transfer(E_B(None), GR_full[_F_IDX])
    G_F = resolve(src_F * phase)

    G_pending = {}
    for idx, (mat_a, mat_b) in enumerate(pending):
        p = idx + 2
        term_bg = idmrg._apply_transfer(E_LR(mat_b), G_F)
        term_src = idmrg._apply_transfer(E_B(mat_b), GR_full[_F_IDX])
        G_pending[p] = (term_bg + term_src) * phase

    rhs_S = idmrg._apply_transfer(E_LR(h1), G_F)
    rhs_S = rhs_S + idmrg._apply_transfer(E_B(h1), GR_full[_F_IDX])
    rhs_S = rhs_S + idmrg._apply_transfer(E_B(None), GR_full[_S_IDX])
    for idx, (mat_a, mat_b) in enumerate(pending):
        p = idx + 2
        rhs_S = rhs_S + idmrg._apply_transfer(E_LR(mat_a), G_pending[p])
        rhs_S = rhs_S + idmrg._apply_transfer(E_B(mat_a), GR_full[p])
    G_S = resolve(rhs_S * phase)

    GBR = {_F_IDX: G_F, _S_IDX: G_S}
    GBR.update(G_pending)
    return GBR


class ExcitationEnvironment:
    """Everything the tangent-space excitation ansatz needs that does NOT
    depend on momentum k -- built once per converged `vumps.VUMPSResult`
    and reused across every `excitation_energies(k)` call, e.g. when
    scanning a dispersion or computing `excitation_gap`.

    `lam_AC` is the constant subtracted from H_eff(k)'s raw eigenvalues to
    get the excitation energy above the ground state -- H_AC's own
    Rayleigh quotient <AC|H_AC[AC]>/<AC|AC> on the converged AC=AL@C, NOT
    the ground state's own physical energy density `e_cell` (also stored,
    for reference/diagnostics only, and equal to `lam_AC - lam_C`, the
    difference of H_AC's and H_C's own eigenvalues at the VUMPS fixed
    point -- the standard VUMPS energy-density identity). This is because
    H_eff(k)'s own "B in center" term (`_full_channel_contraction(GL_full,
    B, GR_full, ...)`) is structurally identical to H_AC[B] (same object,
    X=B instead of X=AC) -- the natural zero-point that makes the ground
    state itself (B=AC) have zero excitation energy under just that term
    is H_AC's own eigenvalue, not the ground state's physical energy
    density. Using `e_cell` here instead was this module's single largest
    remaining discrepancy against MPSKit's own D=2 TFIM oracle after every
    other piece (the GBL/GBR construction, the vumps.py conjugate fix) was
    already exact/Hermitian -- a `lambda_C` = `e_cell - lam_AC` sized
    (`lambda_C` here being H_C's own eigenvalue, generally comparable in
    scale to `e_cell` itself) constant shift, the same at every momentum,
    found by directly comparing `<AC|H_AC[AC]>/<AC|AC>` against the
    dispersion offset empirically before confirming the two matched to
    ~1e-5 relative."""

    def __init__(self, AL, AR, C, GL, GR, W, D, d_g, V_L, GL_full, GR_full,
                 E_RL, E_LR, r_RL, l_RL, r_LR, l_LR, e_cell, lam_AC):
        self.AL, self.AR, self.C = AL, AR, C
        self.GL, self.GR, self.W = GL, GR, W
        self.D, self.d_g = D, d_g
        self.V_L = V_L
        self.pending = _pending_channels(W)
        self.h1 = _onsite_matrix(W)
        self.GL_full, self.GR_full = GL_full, GR_full
        self.E_RL, self.E_LR = E_RL, E_LR
        self.r_RL, self.l_RL = r_RL, l_RL
        self.r_LR, self.l_LR = r_LR, l_LR
        self.e_cell = e_cell
        self.lam_AC = lam_AC


def build_excitation_environment(vumps_result):
    """Build an ExcitationEnvironment from a converged `vumps.VUMPSResult`
    (`vumps.vumps_ground_state`'s own return value) -- unlike an earlier,
    now-superseded version of this module (see this module's own "History"
    section), no separate Hamiltonian term lists / site types need to be
    passed in: AL, AR, C, GL, GR and the grouped automaton W are all
    already on `vumps_result`, and reach<=1 has already been checked
    inside `vumps_ground_state` itself."""
    from . import vumps as _vumps  # lazy: vumps.py imports this module at load time

    AL, AR, C = vumps_result.AL, vumps_result.AR, vumps_result.C
    GL, GR, W = vumps_result.GL, vumps_result.GR, vumps_result.W
    D, d_g = vumps_result.D, vumps_result.d_g
    pending = _pending_channels(W)
    h1 = _onsite_matrix(W)

    I = np.eye(D, dtype=complex)
    bond_envs = _vumps._precompute_bond_environments(AL, AR, pending)
    Lvec_a, Rvec_b = {}, {}
    for idx, (mat_a, mat_b, Lv, Rv) in enumerate(bond_envs):
        Lvec_a[idx + 2] = Lv
        Rvec_b[idx + 2] = Rv
    GL_full = {_S_IDX: I, _F_IDX: GL}
    GL_full.update(Lvec_a)
    GR_full = {_F_IDX: I, _S_IDX: GR}
    GR_full.update(Rvec_b)

    E_RL = _op_transfer_matrix(AR, AL, None)  # ket=AR, bra=AL
    E_LR = _op_transfer_matrix(AL, AR, None)  # ket=AL, bra=AR
    r_RL, l_RL = _mixed_fixed_points(E_RL)
    r_LR, l_LR = _mixed_fixed_points(E_LR)

    V_L = _null_space_left(AL)

    AC = np.einsum('lpm,mr->lpr', AL, C)
    bond_envs_full = [(mat_a, mat_b, GL_full[idx + 2], GR_full[idx + 2])
                       for idx, (mat_a, mat_b) in enumerate(pending)]
    H_AC = _vumps._h_ac_action(AC, GL, GR, bond_envs_full, h1)
    lam_AC = float((np.vdot(AC, H_AC) / np.vdot(AC, AC)).real)

    return ExcitationEnvironment(AL, AR, C, GL, GR, W, D, d_g, V_L,
                                  GL_full, GR_full, E_RL, E_LR,
                                  r_RL, l_RL, r_LR, l_LR,
                                  vumps_result.e_cell, lam_AC)


def _h_eff_action(k, X, env):
    """H_eff(k)[X] -- the momentum-dependent effective Hamiltonian acting
    on a tangent-space parameter X ((D*(d_g-1), D) matrix), returned in the
    same shape. Exactly the 3 terms this module's own docstring (algorithm
    step 4) and `ExcitationEnvironment`'s own docstring describe."""
    D, d_g = env.D, env.d_g
    AL, AR = env.AL, env.AR
    B = (env.V_L @ X).reshape(D, d_g, D)

    Y = _full_channel_contraction(env.GL_full, B, env.GR_full, env.W, env.h1, env.pending)

    GBL = _build_GBL(env, B, k)
    Y = Y + _full_channel_contraction(GBL, AR, env.GR_full, env.W, env.h1, env.pending)

    GBR = _build_GBR(env, B, k)
    Y = Y + _full_channel_contraction(env.GL_full, AL, GBR, env.W, env.h1, env.pending)

    Y_mat = Y.reshape(D * d_g, D)
    return env.V_L.conj().T @ Y_mat


def _build_H_eff_dense(k, env):
    """Dense (Dx*D, Dx*D) matrix representing H_eff(k), Dx=D*(d_g-1) --
    built one basis vector at a time, same style as `_dense_linear_map`.
    No metric reduction is needed (unlike the now-superseded uniform-gauge
    version of this module) -- see this module's own docstring, algorithm
    step 1, for why the mixed-gauge tangent-space norm is already the
    trivial Euclidean one in this (Dx, D)-shaped X basis."""
    D, d_g = env.D, env.d_g
    Dx = D * (d_g - 1)
    n = Dx * D
    H = np.zeros((n, n), dtype=complex)
    Xt = np.zeros((Dx, D), dtype=complex)
    for j in range(n):
        Xt.flat[j] = 1.0
        Y = _h_eff_action(k, Xt, env)
        H[:, j] = Y.reshape(-1)
        Xt.flat[j] = 0.0
    return H


def excitation_energies(env, k, n=1):
    """The lowest `n` excitation energies (above the ground state) at
    momentum `k` (radians, per unit cell) of the tangent-space/
    quasiparticle excitation ansatz.

    Solves the *ordinary* (not generalized -- see this module's own
    docstring, algorithm step 1) Hermitian eigenproblem H_eff(k)[X] =
    lambda*X directly, then subtracts `env.lam_AC` from every raw
    eigenvalue -- see `ExcitationEnvironment`'s own docstring for why this,
    not `env.e_cell`, is the correct renormalization constant."""
    Hmat = _build_H_eff_dense(k, env)
    Hmat = (Hmat + Hmat.conj().T) / 2  # Hermitize (H is Hermitian; this is numerical noise cleanup)
    w = np.linalg.eigvalsh(Hmat)
    w = np.sort(w) - env.lam_AC
    return w[:n]
