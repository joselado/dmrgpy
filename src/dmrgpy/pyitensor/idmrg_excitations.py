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
from scipy.linalg import lu_factor, lu_solve
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


# Every contraction in this section is written with np.tensordot/@ rather
# than np.einsum. These four helpers are the innermost loops of *every*
# infinite-chain code path (idmrg.py's growth loop, vumps.py's H_AC/H_C
# Lanczos matvecs, this module's own excitation ansatz), and np.einsum
# without optimize= runs numpy's own c_einsum kernel, which never dispatches
# to BLAS: profiling a single VUMPS macro-iteration on a 2-site native-
# spinful chain (D=8, d_g=16) spent 7.2s of its 12.0s total inside c_einsum
# across ~102k calls. Measured per-call, on those exact shapes: _apply_op_ket
# 142us -> 37us, _cap_right 75us -> 14us, _op_transfer_matrix 563us -> 98us,
# and the gap *grows* with D (at D=30 _apply_op_ket is 1647us -> 164us, 10x).
# The tensordot forms below are index-for-index the same contractions -- the
# einsum subscripts they replace are kept in each docstring.


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
        ket = _apply_op_ket(M, ket)
    # einsum('lpr,LpR->lLrR', ket, conj(bra)): contract only over the
    # physical leg p, then reorder the four open legs into (l,L,r,R).
    return np.tensordot(ket, np.conj(bra), axes=([1], [1])).transpose(0, 2, 1, 3)


def _apply_op_ket(M, T):
    """T (D,d,D) with operator M (M[i,o] convention) applied to its
    physical leg -- M=None returns T unchanged."""
    if M is None:
        return T
    # einsum('io,aic->aoc', M, T): T's physical leg against M's *in* leg,
    # M's *out* leg then moved back into the physical slot. Written as a
    # broadcasting matmul over T's leading (bond) axis rather than
    # np.tensordot: at these sizes np.tensordot's own Python-level shape
    # bookkeeping costs more than the arithmetic (measured 52us vs 26us per
    # call at D=8, d_g=16 -- and this is the single most-called helper here).
    return M.T @ T


def _cap_right(T, R):
    """T (D,d,D) with its right bond contracted against a (D,D) matrix R
    (ket-type contraction, i.e. R's *first* index matches T's right bond)
    -- used to close T's right bond against a fixed point/environment while
    keeping T's own (left, phys) legs open."""
    # einsum('aoc,cb->aob', T, R): a plain matrix product on T's last leg.
    return (T.reshape(-1, T.shape[2]) @ R).reshape(T.shape[0], T.shape[1], -1)


def _cap_left(L, T):
    """T's left bond contracted against a (D,D) matrix L: L's *first*
    index matches T's left bond (ket-type), L's *second* index becomes the
    new left index. NOT symmetric with `_cap_right` -- L and T play
    structurally different roles here (confirmed by direct derivation
    against idmrg.onsite_expectation's own one-site formula; see git
    history for the transpose bug this once caught)."""
    # einsum('ba,boc->aoc', L, T): contract L's *first* index against T's
    # left leg (see the docstring -- deliberately not the mirror of
    # _cap_right), leaving L's second index as the new left leg. One gemm
    # on T flattened to (left, phys*right), for the same reason
    # _apply_op_ket avoids np.tensordot (42us -> 17us at D=8, d_g=16).
    return (L.T @ T.reshape(T.shape[0], -1)).reshape(T.shape)


# (D*D)-sized linear systems below this size are solved by building the
# map densely and calling np.linalg.solve; above it, iteratively (GMRES on
# the same matvec). The dense build costs D*D applications of the action
# plus an O(D^6) solve, which is what made VUMPS impractical at the bond
# dimensions this module's own docstrings assume: profiling a D=12 run
# spent 5.6s of ~20s inside _dense_linear_map alone, and the cost grows
# like D^6. The threshold is deliberately low rather than "large enough to
# never trigger" -- a test monkeypatches it to 0 to check the two paths
# agree (see tests/test_infinite_chain.py).
_DENSE_SOLVE_MAX = 64

# Relative residual demanded of the iterative solve. This has to be far
# tighter than a typical GMRES default: the environments it computes feed
# straight into the reported energy density, which the infinite-chain tests
# compare at 1e-8..1e-10.
_ITERATIVE_SOLVE_RTOL = 1e-13

# The same dense-vs-iterative decision for `_channel_resolvent`, which needs
# its OWN, much higher threshold than `_solve_linear_map`'s -- the two solve
# genuinely different problems. `_solve_linear_map` is a one-shot solve, so
# there is nothing to amortize an O(D^6) factorization over and the dense
# build is pure overhead. `_channel_resolvent` returns a *reusable* solve
# callable that one `excitation_energies` call invokes ~2 per H_eff
# application; with the per-momentum cache (`_resolvents_for`) the LU
# factorization is paid once and every later solve is two triangular solves,
# which beats re-running GMRES to 1e-13 every time by a wide margin. So the
# threshold here is set by *memory*, not by flops: D*D=2048 is a
# (2048, 2048) complex factorization, ~67 MB.
_RESOLVENT_DENSE_MAX = 2048

# `excitation_energies` builds H_eff(k) densely at or below this dimension
# (dim = D*D*(d_g-1)) and uses an iterative eigensolver above it. Below it,
# assembling the matrix costs `dim` applications of `_h_eff_action` plus a
# dim^3 diagonalization, which is genuinely cheaper -- and exact -- for the
# small cases the test suite is full of; above it the assembly is what
# dominates, since the iterative solve's cost is set by ARPACK's own
# iteration count rather than by `dim`. That count is NOT small: measured on
# a D=12 TFIM chain (dim=144), n=2 took ~155-225 applications per momentum,
# not the "few tens" a first estimate suggests -- which is exactly why the
# threshold sits where it does rather than at 0.
#
# The true crossover depends on `n` (ARPACK needs more iterations for more
# eigenpairs) and so cannot be captured by one number: measured on TFIM,
# n=1 crosses over around dim~36-64 while n=2 still favours dense at
# dim=144. 256 is chosen to sit above the n=2 crossover, i.e. deliberately
# conservative -- the dense path is exact and cannot fail, so paying a
# little extra there is the safe direction to err.
# examples/idmrg/excitation_solver_scaling/main.py plots both curves.
_DENSE_EIG_MAX = 256

# Convergence tolerance handed to the iterative eigensolver, and the largest
# relative residual ||H x - lambda x|| / |lambda| accepted from it before
# falling back to the dense path.
#
# Unlike `_ITERATIVE_SOLVE_RTOL` (a linear solve, whose result feeds onward
# unchecked), this one is measured rather than argued: at D=12 asking ARPACK
# for machine precision (tol=0) cost 905 applications across a 4-momentum
# scan while tol=1e-10 cost 627 -- 30% cheaper for an answer that differed
# from the dense path by 1.15e-14 either way, i.e. no measurable accuracy at
# all. What actually guards accuracy here is the explicit residual check
# below, which is what a too-loose tolerance would trip.
_ITERATIVE_EIG_TOL = 1e-10
_ITERATIVE_EIG_RESIDUAL_MAX = 1e-7


def _solve_linear_map(D, action, rhs):
    """Solve `action(X) = rhs` for a (D,D) complex X, where `action` is a
    linear (D,D)->(D,D) map given as a Python callable.

    Small systems are solved exactly (dense build + LU); larger ones with
    GMRES on the same callable, falling back to the dense path if GMRES
    does not converge -- so this is never *less* robust than the dense
    solve it replaces, only faster when it works."""
    n = D * D
    if n > _DENSE_SOLVE_MAX:
        from scipy.sparse.linalg import LinearOperator, gmres
        op = LinearOperator((n, n), dtype=complex,
                             matvec=lambda x: action(x.reshape(D, D)).reshape(-1))
        try:
            x, info = gmres(op, rhs.reshape(-1), rtol=_ITERATIVE_SOLVE_RTOL,
                             atol=0.0, restart=min(n, 100), maxiter=200)
        except TypeError:  # scipy < 1.12 spells rtol "tol"
            x, info = gmres(op, rhs.reshape(-1), tol=_ITERATIVE_SOLVE_RTOL,
                             atol=0.0, restart=min(n, 100), maxiter=200)
        if info == 0:
            return x.reshape(D, D)
    Mat = _dense_linear_map(D, action)
    return np.linalg.solve(Mat, rhs.reshape(-1)).reshape(D, D)


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

    NOTHING here depends on the excitation tensor B -- only on the momentum
    and on the (k-independent) environment. That is what makes
    `_resolvents_for`'s per-momentum cache correct: one `excitation_energies`
    call applies H_eff(k) many times, each application with a different B,
    and every one of those applications used to rebuild BOTH of this
    momentum's resolvents from scratch. Since the dense build costs D*D
    applications of `action` plus an O(D^6) factorization, that was the
    dominant cost of the whole ansatz at any non-trivial bond dimension.

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

    # Dense path: build the map once, LU-factor it once, and reuse that
    # factorization for every solve this resolvent is asked for. The caller
    # (`_build_GBL`/`_build_GBR`) solves twice per call, and the resolvent
    # itself is now cached across the many `_h_eff_action` calls one
    # eigensolve makes (see `_resolvents_for`), so re-solving from scratch
    # each time was repeating an O(D^6) factorization for nothing.
    lu_cache = []

    def dense_solve(source):
        if not lu_cache:
            lu_cache.append(lu_factor(_dense_linear_map(D, action)))
        return lu_solve(lu_cache[0], source.reshape(-1)).reshape(D, D)

    if D * D <= _RESOLVENT_DENSE_MAX:
        return dense_solve

    # Above the threshold, never build the (D*D, D*D) matrix at all --
    # GMRES on the same callable, exactly as `_solve_linear_map` does (and
    # with the same tolerance, for the same reason: what this solves for
    # feeds straight into the reported dispersion). Falls back to the dense
    # path if GMRES does not converge, so this is never less robust than the
    # dense solve, only cheaper.
    from scipy.sparse.linalg import LinearOperator, gmres

    n = D * D
    op = LinearOperator((n, n), dtype=complex,
                         matvec=lambda x: action(x.reshape(D, D)).reshape(-1))

    def solve(source):
        rhs = source.reshape(-1)
        try:
            x, info = gmres(op, rhs, rtol=_ITERATIVE_SOLVE_RTOL, atol=0.0,
                             restart=min(n, 100), maxiter=200)
        except TypeError:  # scipy < 1.12 spells rtol "tol"
            x, info = gmres(op, rhs, tol=_ITERATIVE_SOLVE_RTOL, atol=0.0,
                             restart=min(n, 100), maxiter=200)
        if info == 0:
            return x.reshape(D, D)
        return dense_solve(source)

    return solve


def _resolvents_for(env, k):
    """(resolve_L, resolve_R): this momentum's two `_channel_resolvent`
    solve-callables, built on first use and cached on `env`.

    Both depend only on `k` and on the momentum-independent environment,
    never on the excitation tensor B (see `_channel_resolvent`'s own
    docstring), so one cache entry serves every `_h_eff_action` call at this
    momentum -- which is what turns the eigensolve's per-application cost
    from "rebuild and factorize two D^2-by-D^2 maps" into two triangular
    solves.

    Keyed on the momentum's exact float bit pattern rather than a rounded
    value: callers pass the same `k` object through a whole eigensolve, so
    exact equality always hits, and a near-miss would only cost a rebuild,
    never a wrong answer."""
    cache = env.resolvent_cache
    key = float(k)
    if key not in cache:
        phase = np.exp(1j * k)
        cache[key] = (
            _channel_resolvent(1.0 / phase, env.E_RL, env.r_RL, env.l_RL,
                                idmrg._apply_transfer_from_left, env.D),
            _channel_resolvent(phase, env.E_LR, env.r_LR, env.l_LR,
                                idmrg._apply_transfer, env.D),
        )
    return cache[key]


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

    resolve = _resolvents_for(env, k)[0]

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

    resolve = _resolvents_for(env, k)[1]

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
        # {k: (resolve_L, resolve_R)} -- see `_resolvents_for`. Lives on the
        # environment because the environment is itself cached per converged
        # ground state (infinitechain.py's `_excitation_env`), so a momentum
        # scan that revisits a k pays for its resolvents exactly once.
        self.resolvent_cache = {}
        # The same two caches for the spectral-weight machinery, filled
        # lazily by `_spectral_cache`/`_spectral_resolvents_for` -- lazily
        # rather than here, so a chain that only ever asks for a dispersion
        # never pays for either (see `_spectral_source_vector`).
        self._spectral_pieces = None
        self._spectral_resolvent_cache = {}


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


def _lowest_dense(k, env, n):
    """(w, X_list): the `n` lowest eigenpairs of H_eff(k), by building the
    matrix and diagonalizing it. Exact, and the cheaper of the two paths
    when `dim` is small -- see `_DENSE_EIG_MAX`."""
    D, d_g = env.D, env.d_g
    Dx = D * (d_g - 1)
    Hmat = _build_H_eff_dense(k, env)
    # Hermitize. H_eff is Hermitian analytically; what this removes is
    # numerical noise, measured at ~1e-16 relative for k away from 0 and
    # ~1e-11 at exactly k=0, where `_channel_resolvent`'s own regularizing
    # projector is armed and its residual inherits the ground state's own
    # VUMPS convergence tolerance.
    Hmat = (Hmat + Hmat.conj().T) / 2
    w, V = np.linalg.eigh(Hmat)
    order = np.argsort(w)[:n]
    return w[order], [V[:, j].reshape(Dx, D) for j in order]


def _lowest_iterative(k, env, n):
    """(w, X_list) via Lanczos on `_h_eff_action` directly, or None if the
    iterative solve did not produce eigenpairs good enough to trust.

    This is the whole point of not assembling H_eff(k): the dense build
    costs `dim` = D*D*(d_g-1) applications of `_h_eff_action`, while ARPACK
    needs a few hundred of them regardless of `dim` (measured, see
    `_ITERATIVE_EIG_TOL`) -- so this wins once `dim` is the larger of the
    two, and `_DENSE_EIG_MAX` is where that happens. Returning None rather
    than a bad answer keeps the same contract `_solve_linear_map` has -- the
    iterative path is never allowed to be *less* reliable than the dense one
    it replaces, only cheaper.

    Two details that are not optional:

    - `v0` is fixed and deterministic. ARPACK's default start vector is
      random, which would make a near-degenerate or gapless dispersion come
      out slightly differently run to run -- exactly the regime the n_uc=2
      Heisenberg cross-checks live in.
    - the eigenpairs are residual-checked before being returned, which also
      covers the Hermiticity question: Lanczos assumes a Hermitian operator,
      and if the H-vs-H^dagger asymmetry that `_lowest_dense` averages away
      were ever structural rather than noise, it would show up here as a
      residual this check rejects."""
    from scipy.sparse.linalg import ArpackError, ArpackNoConvergence
    from scipy.sparse.linalg import LinearOperator, eigsh

    D, d_g = env.D, env.d_g
    Dx = D * (d_g - 1)
    dim = Dx * D

    def matvec(x):
        return _h_eff_action(k, x.reshape(Dx, D), env).reshape(-1)

    op = LinearOperator((dim, dim), dtype=complex, matvec=matvec)
    v0 = np.ones(dim, dtype=complex) / np.sqrt(dim)
    try:
        w, V = eigsh(op, k=n, which="SA", v0=v0, tol=_ITERATIVE_EIG_TOL)
    except (ArpackNoConvergence, ArpackError):
        return None

    order = np.argsort(w)
    w, V = w[order], V[:, order]
    scale = max(np.max(np.abs(w)), 1.0)
    for j in range(len(w)):
        resid = np.linalg.norm(matvec(V[:, j]) - w[j] * V[:, j])
        if resid / scale > _ITERATIVE_EIG_RESIDUAL_MAX:
            return None
    return w, [V[:, j].reshape(Dx, D) for j in range(len(w))]


def excitation_energies(env, k, n=1, return_vectors=False):
    """The lowest `n` excitation energies (above the ground state) at
    momentum `k` (radians, per unit cell) of the tangent-space/
    quasiparticle excitation ansatz.

    Solves the *ordinary* (not generalized -- see this module's own
    docstring, algorithm step 1) Hermitian eigenproblem H_eff(k)[X] =
    lambda*X, then subtracts `env.lam_AC` from every raw eigenvalue -- see
    `ExcitationEnvironment`'s own docstring for why this, not `env.e_cell`,
    is the correct renormalization constant.

    Two solvers, picked by problem size (`_DENSE_EIG_MAX`): assemble
    H_eff(k) and call `eigh` for small ones, Lanczos on `_h_eff_action`
    directly for large ones, falling back to the dense path if the
    iterative solve does not converge to an acceptable residual. Both are
    exercised by the test suite (the threshold is monkeypatched, the same
    way `_DENSE_SOLVE_MAX` is). The dense path is not vestigial: `dim` can
    legitimately be 1 (a D=1 spin-1/2 chain), where `eigsh` cannot be used
    at all since ARPACK needs `nev < dim`.

    `return_vectors=True` additionally returns the tangent-space parameters
    X (a list of (D*(d_g-1), D) arrays, one per returned energy, in the same
    order). The excitation tensor itself is B = (V_L @ X).reshape(D, d_g, D)
    -- what a spectral weight <Psi|O|Phi_k(B)> is built from."""
    dim = env.D * env.D * (env.d_g - 1)
    n = min(n, dim)

    got = None
    if dim > _DENSE_EIG_MAX and n < dim:  # ARPACK needs nev < dim
        got = _lowest_iterative(k, env, n)
    if got is None:
        got = _lowest_dense(k, env, n)

    w, X_list = got
    w = np.asarray(w) - env.lam_AC
    if return_vectors:
        return w, X_list
    return w


# ---------------------------------------------------------------------------
# Spectral weights: <Phi_k(B)|O(k)|Psi>, the exact delta-peak weight of each
# excitation branch (docs/idmrg_improvement_plan.md item 2).
# ---------------------------------------------------------------------------

def _spectral_resolvent(phase_factor, E, X0, apply_fn, D):
    """A `solve(source) -> G` callable for the *transposed*, deflated
    geometric-series map

        (I - phase_factor*T + Pi)[G] = source,   T[G] = apply_fn(E, G)

    where `T`'s dominant eigenvalue is exactly 1 with right eigenvector the
    identity and left eigenvector `X0` (normalized to trace(X0)=1), and
    `Pi[G] = I * sum(X0 * G)` is that eigenvalue's own spectral projector
    written in the *bilinear* pairing `<X,Y> = sum_a X[a]*Y[a]` (no
    conjugation anywhere -- this whole helper lives on the transposed side
    of `_spectral_source_vector`'s own derivation, see its docstring).

    `Pi` is added at EVERY momentum, not only at the singular one. That is
    deliberate and is exactly what makes this well-conditioned near k=0:
    deflating replaces the identity direction's own eigenvalue 1-phase_factor
    (which vanishes at k=0, and is merely tiny just next to it) by
    2-phase_factor, leaving every other direction untouched. The resulting
    G differs from the un-deflated solution by a multiple of the identity,
    and `_spectral_source_vector` contracts that multiple into a multiple of
    AC, which its own final V_L^dagger projection annihilates exactly (V_L
    is built as the null space of AL_mat^dagger, and AC = AL_mat @ C) -- so
    the returned spectral weight is unchanged while the conditioning is
    bounded uniformly in k. Doing this only at k==0, the way
    `_channel_resolvent` has to, would leave a genuinely ill-conditioned
    solve at every small-but-nonzero momentum.

    Same dense-vs-iterative split (and same threshold) as
    `_channel_resolvent`: LU-factor the (D*D,D*D) map once and reuse it,
    or GMRES on the callable above `_RESOLVENT_DENSE_MAX`."""
    def action(G):
        return G - phase_factor * apply_fn(E, G) + np.eye(D) * np.sum(X0 * G)

    lu_cache = []

    def dense_solve(source):
        if not lu_cache:
            lu_cache.append(lu_factor(_dense_linear_map(D, action)))
        return lu_solve(lu_cache[0], source.reshape(-1)).reshape(D, D)

    if D * D <= _RESOLVENT_DENSE_MAX:
        return dense_solve

    from scipy.sparse.linalg import LinearOperator, gmres

    n = D * D
    op = LinearOperator((n, n), dtype=complex,
                         matvec=lambda x: action(x.reshape(D, D)).reshape(-1))

    def solve(source):
        rhs = source.reshape(-1)
        try:
            x, info = gmres(op, rhs, rtol=_ITERATIVE_SOLVE_RTOL, atol=0.0,
                             restart=min(n, 100), maxiter=200)
        except TypeError:  # scipy < 1.12 spells rtol "tol"
            x, info = gmres(op, rhs, tol=_ITERATIVE_SOLVE_RTOL, atol=0.0,
                             restart=min(n, 100), maxiter=200)
        if info == 0:
            return x.reshape(D, D)
        return dense_solve(source)

    return solve


def _spectral_cache(env):
    """The momentum-independent pieces `_spectral_source_vector` needs,
    built on first use and stashed on the environment: (AC, E_RR, E_LL,
    rho_R, rho_L).

    `rho_R = C^dagger C` and `rho_L = C C^dagger` are the two transfer
    matrices' own non-trivial fixed points, written in closed form rather
    than found by an eigensolve -- both are exact identities of the mixed
    canonical gauge (AC = AL@C = C@AR, AL left-orthonormal, AR right-
    orthonormal), confirmed directly: propagating C^dagger C forward
    through E_RR gives sum_p AC_p^dagger AC_p = C^dagger (sum_p AL_p^dagger
    AL_p) C = C^dagger C, and the mirror for C C^dagger through E_LL. The
    *other* fixed point of each is the identity on both sides, by the same
    orthonormality. Normalized to trace 1 so the deflation projector
    `_spectral_resolvent` builds from them is idempotent."""
    cached = env._spectral_pieces
    if cached is not None:
        return cached
    AL, AR, C = env.AL, env.AR, env.C
    AC = np.tensordot(AL, C, axes=([2], [0]))
    rho_R = C.conj().T @ C
    rho_L = C @ C.conj().T
    rho_R = rho_R / np.trace(rho_R)
    rho_L = rho_L / np.trace(rho_L)
    cached = (AC, _op_transfer_matrix(AR, AR), _op_transfer_matrix(AL, AL),
              rho_R, rho_L)
    env._spectral_pieces = cached
    return cached


def _spectral_resolvents_for(env, k):
    """(solve_R, solve_L): this momentum's two `_spectral_resolvent`
    callables, cached on the environment exactly the way `_resolvents_for`
    caches the H_eff(k) ones -- the matrices depend only on `k`, never on
    the operator, so one entry serves every operator asked for at this
    momentum (and every branch, since `_spectral_source_vector` needs the
    two solves once in total, not once per branch)."""
    cache = env._spectral_resolvent_cache
    key = float(k)
    if key not in cache:
        _, E_RR, E_LL, rho_R, rho_L = _spectral_cache(env)
        phase = np.exp(1j * k)
        cache[key] = (
            _spectral_resolvent(phase, E_RR, rho_R, idmrg._apply_transfer, env.D),
            _spectral_resolvent(1.0 / phase, E_LL, rho_L,
                                 idmrg._apply_transfer_from_left, env.D),
        )
    return cache[key]


def _spectral_source_vector(env, M, k):
    """v: the (Dx, D) matrix (Dx = D*(d_g-1)) representing the linear
    functional `X -> <Phi_k(V_L X) | O(k) | Psi>` as `trace(X^dagger @ v)`,
    for the supersite-local operator `M` (a d_g x d_g dense matrix in this
    module's own M[i,o] convention) and momentum `k` (radians per unit
    cell).

    == What is being computed ==

    With this module's own momentum convention |Phi_k(B)> = sum_n e^{ikn}
    |...AL AL B_n AR AR...> and O(k) = (1/sqrt(N)) sum_m e^{ikm} O_m, the
    overlap of the *normalized* ansatz state (norm sqrt(N * trace(X^dagger
    X)), the cross terms vanishing by the left gauge fixing) with O(k)|Psi>
    is exactly the per-site, N-independent quantity

        S(k) = sum_j e^{ikj} <B at 0 | O_j | Psi>

    (both factors of N cancel), which is what this returns in the X basis.
    Vanderstraeten, Haegeman & Verstraete, arXiv:1810.07006 Sec. 6.3
    ("dynamical correlations"), and Osborne & McCulloch, arXiv:2408.17117,
    for the same object written recursively.

    == The three regions ==

    Put the ket in mixed canonical form with its center at site 0 (AL to
    the left, AC at 0, AR to the right -- a valid gauge choice at any cut,
    see vumps.two_point_correlator's own docstring) and the bra with its
    excitation tensor B at site 0. Then, since AL is left- and AR is
    right-orthonormal, both semi-infinite background halves close to the
    identity and only three regions survive:

    - `j == 0`: the operator sits on the excitation itself, contributing
      `_apply_op_ket(M, AC)` directly.
    - `j > 0`: an (AR,AR) transfer chain to the right of the excitation,
      resummed as a geometric series `sum_{m>=0} (e^{ik} E_RR)^m` and
      capped by one operator-carrying transfer and the identity.
    - `j < 0`: the mirror, an (AL,AL) chain with `e^{-ik}`.

    Both series converge even at k=0 despite the transfer matrices' own
    unit dominant eigenvalue, because the source of each is orthogonal to
    the corresponding dominant eigenvector: both reduce to
    `sum conj(B) * AC = trace(X^dagger V_L^dagger AL_mat C) = 0`, which is
    the left gauge-fixing condition itself. The deflation
    `_spectral_resolvent` applies is therefore a conditioning device only,
    never a physical regularization -- see its own docstring.

    == Why the solves are transposed ==

    Written forwards, each geometric series takes a B-dependent source and
    the whole expression is one scalar, so producing `v` (i.e. keeping B's
    legs open) would need one linear solve per basis element of B --
    D^2*d_g solves per momentum. Transposing each solve instead (the map
    `_apply_transfer_from_left` and `_apply_transfer` are each other's
    transpose under the bilinear pairing, with the *same* E) moves the
    solve to the B-independent end of the contraction, so the whole source
    vector -- every branch's weight, and the sum rule below -- costs
    exactly two solves per momentum regardless of D, d_g or how many
    branches the caller wants.

    == Sum rule ==

    For a supersite-local O, O(k)|Psi> lies *exactly* inside the
    tangent space (B = O.A is itself a valid excitation tensor), so
    summing |<X_alpha, v>|^2 over a complete eigenbasis {X_alpha} of
    H_eff(k) -- i.e. `norm(v)**2` -- is the per-site *connected* static
    structure factor sum_r e^{ikr} (<O_0 O_r> - <O>^2). The connectedness
    is automatic and needs no explicit <O> subtraction: the disconnected
    piece lives entirely along the identity direction the transposed solves
    leave undetermined, and `V_L^dagger` annihilates it (see
    `_spectral_resolvent`). `spectral_weights(..., return_total=True)`
    exposes this number, and tests/test_infinite_chain_spectral.py checks
    it against an independent real-space correlator sum."""
    D, d_g = env.D, env.d_g
    AL, AR = env.AL, env.AR
    AC, E_RR, E_LL, _, _ = _spectral_cache(env)
    solve_R, solve_L = _spectral_resolvents_for(env, k)
    phase = np.exp(1j * k)
    I = np.eye(D, dtype=complex)

    # j == 0.
    v = _apply_op_ket(M, AC)

    # j > 0: cap = trace(apply_transfer_from_left(E_RR_O, Y)), so the
    # transposed chain starts from apply_transfer(E_RR_O, I).
    G = solve_R(idmrg._apply_transfer(_op_transfer_matrix(AR, AR, M), I))
    v = v + phase * _cap_right(AC, G)

    # j < 0: the mirror. cap = trace(apply_transfer(E_LL_O, Y)), so the
    # transposed chain starts from apply_transfer_from_left(E_LL_O, I).
    G = solve_L(idmrg._apply_transfer_from_left(_op_transfer_matrix(AL, AL, M), I))
    v = v + _cap_left(G, AC) / phase

    return env.V_L.conj().T @ v.reshape(D * d_g, D)


def spectral_weights(env, k, M, n=1, return_total=False):
    """(energies, weights): the lowest `n` excitation energies at momentum
    `k` (radians per unit cell) together with their exact delta-peak
    spectral weights for the supersite-local operator `M` (a d_g x d_g
    dense matrix, M[i,o] convention -- `vumps._embed_group_operator` builds
    one from a named operator on a chosen sub-site of the unit cell).

    `weights[a] = |<k,a| O(k) |Psi>|^2` for the normalized ansatz state
    `|k,a>` and `O(k) = N^{-1/2} sum_m e^{ikm} O_m`, i.e. exactly the
    residue of the delta peak that branch contributes to the zero-
    temperature dynamical structure factor

        S(k, w) = sum_a weights[a] * delta(w - energies[a])

    -- a *momentum-resolved* dynamical correlator directly in the
    thermodynamic limit, with no finite window and no broadening beyond
    whatever the caller chooses to apply afterwards. Both arrays are
    ordered by energy, matching `excitation_energies`' own ordering.

    `return_total=True` additionally returns the branch-complete total
    `sum_a |<k,a|O(k)|Psi>|^2` over EVERY branch, not just the `n`
    requested -- the per-site connected static structure factor at this
    momentum (see `_spectral_source_vector`'s own "Sum rule" section for
    why that identity is exact here). Comparing `weights.sum()` against it
    says how much of the k-resolved spectral weight the branches actually
    returned carry, which is the practical measure of whether a
    single-mode picture is adequate at that momentum. It costs nothing
    extra: it is the squared norm of the same source vector every weight is
    an inner product against.

    Two things to know before reading an individual weight:

    - **Within a degenerate multiplet the split between branches is
      basis-arbitrary**, since `numpy.linalg.eigh`/ARPACK pick an
      arbitrary basis of a degenerate eigenspace. Only the
      multiplet-summed weight is physical. This is not a corner case on a
      symmetric model: on the AKLT chain at D=2 the eight branches split
      into an SU(2) triplet and a quintuplet at every momentum, and the
      triplet's three individual weights come out ~0.08/0.84/0.08 while
      their sum is the whole `return_total`. At finite D on a model with
      no exact MPS ground state the multiplet is only *approximately*
      degenerate -- split by the variational error, not by machine
      precision: measured at 1.6e-5 on the Haldane chain at D=16, against
      a 1e-13 splitting at D=12 of the same chain, i.e. not even monotonic
      in D, since it tracks how well that particular VUMPS run converged.
      So group branches by the multiplet's own size (the physics), not by
      a fixed energy tolerance -- a 1e-6 tolerance silently reports one
      member's arbitrary share (0.43 instead of 0.97) on exactly that
      D=16 case.
    - **A near-zero fraction means the wrong branches were asked for, not
      that there is no response.** The operator can be forbidden by a
      selection rule from reaching the branches returned -- exactly what
      happens above the multiplet crossing on that same AKLT chain, where
      below k~0.9 the *lowest* branch is the quintuplet and `Sz` cannot
      reach it at all (weight ~1e-23), while the triplet carrying all of
      the weight sits above it. Raising `n` is the fix; `return_total` is
      what makes the situation visible rather than silent.

    The momentum sign convention is this module's own (|Phi_k> = sum_n
    e^{ikn}|...B_n...>); flipping it relabels k <-> -k and nothing else."""
    w, X_list = excitation_energies(env, k, n=n, return_vectors=True)
    v = _spectral_source_vector(env, M, k)
    weights = np.array([abs(np.vdot(X, v)) ** 2 for X in X_list])
    if return_total:
        return w, weights, float(np.vdot(v, v).real)
    return w, weights
