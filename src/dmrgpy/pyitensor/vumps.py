"""VUMPS (Variational Uniform Matrix Product States) ground-state solver
(Zauner-Stauber, Vanderstraeten, Fishman, Verstraete, Haegeman, "Variational
optimization algorithms for uniform matrix product states", arXiv:1701.07035;
see also Vanderstraeten, Haegeman, Verstraete, "Tangent-space methods for
uniform matrix product states", arXiv:1810.07006, Algorithm 4) -- an
alternative to idmrg.py's growing/infinite-size algorithm that solves
directly, at a FIXED target bond dimension D and in the thermodynamic limit,
for the mixed-gauge {AL, AR, C} fixed point of a single (super)site's
transfer/Hamiltonian environment, rather than growing a finite window and
taking its converged edge tensors as an approximation to the infinite chain.

== Why this, on top of idmrg.py's existing growing algorithm ==

idmrg_excitations.py's tangent-space excitation ansatz needs a
self-consistent mixed gauge {AL, AR, C} (AL @ C = C @ AR, C C^dagger = the
AL-transfer's own right fixed point, C^dagger C = the AR-transfer's own left
fixed point) to build its GBL/GBR environments correctly for a genuinely
entangled (D>1) ground state -- see that module's own module docstring
("History" section) for the multi-pass investigation this unblocked.
Reconstructing AR/C post-hoc from idmrg.py's single left-canonical tensor A
(via an eigendecomposition square root of A's own dominant right fixed
point r) was tried early on and found to produce a non-Hermitian H_eff(k)
-- a sign that a remaining, unisolated sign/gauge error persisted in that
reconstruction. VUMPS sidesteps the reconstruction entirely: AL, AR, C are
ALL primary variables of the algorithm itself, solved for jointly and
self-consistently, with GL/GR (this module's own name for
idmrg_excitations.py's own ordinary background environments, built the same
way but from AL/AR directly rather than from a single A/r) falling out as a
natural intermediate of every iteration -- there is no separate
reconstruction step left to get wrong. `build_excitation_environment` in
idmrg_excitations.py takes a converged `VUMPSResult` (this module's own
return value) directly, at any bond dimension D>=1.

A second, independent advantage over idmrg.py's growing algorithm: VUMPS
solves directly at a fixed target bond dimension D in the thermodynamic
limit, with no finite-window growth/truncation step -- idmrg_ground_state
instead grows a window indefinitely and *truncates* down to maxm at every
micro-step, so its own converged U_list is only ever an approximation to the
true D=maxm optimum (whatever residual truncation error survives from every
earlier, smaller-window growth step). VUMPS's {AL,AR,C} are the actual
D-dimensional variational optimum of the energy (up to Lanczos/outer-loop
convergence), not a truncated growth product.

== Scope ==
Same as idmrg_excitations.py: n_uc in {1,2} (grouped into one effective
supersite via `_group_automaton` -- there is no existing per-sublattice
ket tensor list to group here, VUMPS builds AL/AR from scratch), reach-1
bonds after grouping (idmrg_excitations._check_reach_one), pyitensor
engine only. Static
correlators (onsite_expectation/two_point_correlator-style) are not
implemented for a VUMPSResult -- only the converged energy density -- mirroring
the precedent Infinite_Many_Body_Chain already sets for itensor_version=3
(gs_energy works, vev/correlator raise NotImplementedError).

== Convergence robustness (honest scope note) ==

Plain single-attempt VUMPS from a random start turned out, during
development, to be genuinely unreliable for D>1 -- confirmed directly to
reach self-consistent (gauge_mismatch<tol) fixed points at an ENERGY WORSE
than the same model's own smaller-D result, which cannot happen for the
true variational optimum. `vumps_ground_state` mitigates this two ways:
(1) a D-ramp (warm-start each bond dimension from the previous one's own
best solution) plus multiple restarts at each step, and (2) a safety net
that explicitly enforces the variational principle across the ramp -- if a
given D's own restart search still lands worse than an already-known
smaller-D energy, a bounded extra budget of attempts -- warm-started from
that SAME already-good smaller-D solution each time (fresh noise per
attempt), not resampled fully at random -- specifically tries to beat it
before the ramp moves on (see `vumps_ground_state`'s own inline comment).

A large part of what used to make D>1 unreliable beyond what (1)/(2) alone
fix (an unguarded D=4 TFIM(g=1.5) run once landed at a wildly wrong, even
POSITIVE, energy despite a known-achievable D=1 reference around -1.56;
repeated calls at D=4 occasionally, roughly 1 in 10, landed ~10% away from
the exact answer even with both mitigations) turned out to be a genuine,
separate bug, not merely restart-search difficulty: `_solve_left_
environment`/`_solve_right_environment`/`_energy_density_and_source_from_
left`/`_right` were all closing an `apply_transfer_from_left`/
`apply_transfer` output against the WRONG (missing a conjugate) fixed
point -- see those functions' own docstrings for the derivation, and
idmrg_excitations.py's own "History" section for how this was found
(cross-checking the D>1 tangent-space excitation ansatz this module feeds
directly into against an independently-converged MPSKit.jl state).
Invisible at D=1 (where the relevant fixed point is always real), but a
real, D>1-only source of wrong energies on top of ordinary restart-search
difficulty. With the fix, the specific D=4 TFIM(g=1.5) case flagged above
now converges to the exact answer to ~1e-7 relative on 10/10 independent
`nrestarts=6` calls (previously confirmed to occasionally miss by ~10%) --
a measurable, direct improvement, though not a guarantee that restart
search itself can never land in a bad basin for a harder model/larger D
(the D-ramp/multi-restart/safety-net machinery above remains real,
load-bearing infrastructure, not vestigial). `.converged` is always
checked and reported honestly (never silently assumed, matching
idmrg_ground_state's own `.converged` discipline), and every fixed point
actually found does satisfy the VUMPS equations to the requested `tol`. A
caller needing a reliable D>2 result on a harder/less-tested model should
still call `vumps_ground_state` several times independently (not just via
one call's own `nrestarts`) and keep the lowest reported `e0`, the same
"rerun and take the best" discipline idmrg_ground_state's own
unseeded-random-MPS-initialization docstring already recommends.

== Algorithm (single-site VUMPS, Zauner-Stauber Algorithm 2) ==

Given a converged (AL, AR, C) that satisfies the mixed-gauge condition
exactly, GL (the ordinary, regularized background Hamiltonian environment
built from AL) and GR (mirror, from AR) close the two effective
Hamiltonians:

    H_AC[X] = onsite(X) + cap_right(X,GR) + cap_left(GL,X)
              + sum_bonds [ cap_right(mat_a(X), Rvec_b) + cap_left(Lvec_a, mat_b(X)) ]
    H_C[X]  = cap_right(X,GR) + cap_left(GL,X)
              + sum_bonds [ Lvec_a @ X @ Rvec_b ]

(Rvec_b/Lvec_a: the one-more-site content of a reach-1 bond, built by
applying mat_b to AR and closing with AR's own right fixed point I, or mat_a
to AL and closing with AL's own left fixed point I -- see
`_precompute_bond_environments`.) This is exactly the "B in center" term
idmrg_excitations.py's own `_full_channel_contraction` builds (with X in
place of B, and no momentum sum -- there is no momentum sum here, X sits
at a single, fixed position, connecting only to its two immediate
neighbors AL/AR), reusing the SAME cap_right/cap_left/apply_op_ket
primitives that module already validated (and, in the other direction,
`idmrg_excitations._full_channel_contraction(GL_full, X, GR_full, ...)`
reproduces `_h_ac_action` below exactly -- confirmed directly, to
~1e-16).

One VUMPS iteration, given the current (AL, AR, C):

1. r_AL := AL's own fresh dominant right fixed point (AL is always exactly
   left-canonical by construction -- see step 5 -- so its own left fixed
   point is exactly I, but its right fixed point is a genuine eigenproblem,
   generally != C C^dagger until convergence). l_AR mirrors this for AR.
2. e := the current energy-density estimate, averaged from AL's own and AR's
   own one-step Hamiltonian content (`_energy_density_and_source_from_left`/
   `_right`) -- the two individually agree only once AL/AR become mutually
   consistent, so early iterations use their average as a single, shared
   regularization target for both GL and GR (both must be regularized
   against the SAME e, or H_AC/H_C's own two background pieces would carry
   inconsistent zero-point shifts).
3. GL, GR: solve a regularized (I-T+projector) dense linear system,
   one-sided (AL alone for GL, AR alone for AR) -- see `_solve_left_
   environment`/`_solve_right_environment`.
4. Solve H_AC[AC] = lambda*AC and H_C[C] = mu*C for their LOWEST eigenpairs
   (Lanczos, dmrg._lanczos_ground_state, warm-started from the previous
   iteration's own AC/C).
5. Update AL, AR from the new AC, C via the standard least-squares
   (orthogonal-Procrustes) formula: AL = argmin_{AL isometric} ||AC-AL C||,
   solved by polar-decomposing AC = U_l P_l (U_l isometric) and C = U_c P_c
   (U_c unitary), giving AL = U_l U_c^dagger (and the mirror formula, right
   polar factors, for AR) -- see `_update_AL_AR`.

Convergence is measured as ||AC - AL@C|| + ||AC - C@AR|| (both -> 0 exactly
at a true VUMPS fixed point, by construction of step 5 -- this is the
standard "AC-C gauge mismatch" diagnostic most VUMPS implementations use).
"""

import numpy as np
from scipy.linalg import polar as _polar

from . import idmrg
from . import idmrg_excitations as idmrg_exc
from .dmrg import _lanczos_ground_state
from .index import Index
from .tensor import ITensor


def _group_automaton(W_bulk, n_uc):
    """The grouped, single-supersite automaton W (Dw,Dw,d_g,d_g) -- VUMPS
    has no existing per-sublattice ket tensor list to group (unlike
    idmrg.py's growing algorithm), so only the automaton itself needs
    grouping here; AL/AR are built directly at the grouped-supersite
    level from scratch (see `_random_initial_state`)."""
    W = W_bulk[0].array
    for p in range(1, n_uc):
        Wp = W_bulk[p].array
        W = np.einsum('Labm,mcdR->LacbdR', W, Wp)
        L, di0, di1, do0, do1, R = W.shape
        W = W.reshape(L, di0 * di1, do0 * do1, R)
    return W


def _random_initial_state(D, d_g):
    """(AL0, AR0, C0): a crude, mutually-INconsistent starting point for the
    VUMPS iteration -- AL0/AR0 are each individually exactly canonical
    (left/right respectively, verified by construction below), but are
    canonicalizations of the SAME raw random tensor from two different
    directions, so they generally do not yet satisfy AL0 @ C0 = C0 @ AR0 for
    any C0. This does not matter for an iterative fixed-point method: C0 is
    simply started at the identity, and the VUMPS iteration itself is what
    drives (AL,AR,C) toward mutual consistency (see this module's own
    docstring, step 5) -- there is no need for a consistent starting gauge.

    AL0 is obtained by reusing idmrg._canonicalize_periodic (the same
    two-sided fixed-point canonicalization apply_mpo/imps_sum already rely
    on) on a single random (D,d_g,D) tensor treated as a trivial n_uc=1
    periodic chain -- already-validated machinery, not a fresh derivation.
    AR0 is obtained by the standard reversal trick: canonicalizing the same
    raw tensor with its left/right bond legs swapped gives a LEFT-canonical
    tensor for the *reversed* chain, which is exactly a RIGHT-canonical
    tensor once swapped back (verified algebraically: if left-canonicality
    of the swapped tensor reads sum_{a,p} conj(T[a,p,m])T[a,p,m']=delta,
    swapping back exchanges which index pair is "left" vs "right" so this
    identity becomes exactly AR0's own right-canonicality condition
    sum_{p,r} AR0[l,p,r]conj(AR0[l',p,r])=delta).

    A fresh random tensor's own self-overlap transfer eigenvalue eta is
    generically far from 1 (confirmed directly: eta in the 4-15 range for
    small D*d_g here) -- _canonicalize_periodic's own G_left/G_right
    construction is a pure re-gauging (X_p/Y_p are individually
    trace-normalized to 1 regardless of eta, by _dominant_right/left_
    fixed_point's own convention), NOT a renormalization, so it silently
    reproduces the ORIGINAL eta in its output gram matrix instead of
    Identity whenever the input isn't already close to eta=1 -- confirmed
    directly (gram off from Identity by exactly eta, e.g. ~4.44 for a
    D=1 test case) and, checking every existing caller, never previously
    exercised: apply_mpo/imps_sum always feed it tensors built from an
    already-converged (eta=1, by left-canonical SVD construction)
    IDMRGResult. Dividing the raw random tensor by sqrt(its own eta) first
    (computed via idmrg._dominant_right_fixed_point) avoids this precondition
    gap entirely -- confirmed directly to bring the canonicalized output's own
    eta to 1 (to ~1e-15) for D in (1,3,4)."""
    A0 = _random_raw_tensor(D, d_g)
    AL0 = _canonicalize_raw(A0, D, d_g)
    A0_rev = np.transpose(A0, (2, 1, 0)).copy()
    AL0_rev = _canonicalize_raw(A0_rev, D, d_g)
    AR0 = np.transpose(AL0_rev, (2, 1, 0)).copy()

    C0 = np.eye(D, dtype=complex)
    return AL0, AR0, C0


def _random_raw_tensor(D, d_g, seed_AL=None, noise=0.05):
    """A raw (D,d_g,D) tensor to feed `_canonicalize_raw`: pure random
    noise if `seed_AL` is None, otherwise `seed_AL` (a smaller D_old<=D
    already-left-canonical tensor -- see `_grow_initial_state`) embedded
    into the top-left (D_old,d_g,D_old) block plus small random `noise` on
    every entry (including the embedded block, so the larger bond
    dimension's new directions aren't started exactly at zero, which would
    leave them exactly decoupled/unreachable by the D_old block's own
    already-converged content under a linear map like the transfer
    matrix)."""
    rng = np.random.default_rng()
    A0 = noise * (rng.standard_normal((D, d_g, D)) + 1j * rng.standard_normal((D, d_g, D)))
    if seed_AL is not None:
        D_old = seed_AL.shape[0]
        A0[:D_old, :, :D_old] += seed_AL
    return A0


def _canonicalize_raw(A0, D, d_g):
    """A0 (D,d_g,D), pre-normalized to its own eta=1 (see this function's
    own docstring in `_random_initial_state` -- the precondition
    idmrg._canonicalize_periodic silently assumes but does not itself
    enforce) and left-canonicalized via idmrg._canonicalize_periodic."""
    left, right, phys = (Index(D, tags="Link"), Index(D, tags="Link"),
                          Index(d_g, tags="Site"))
    _, eta0 = idmrg._dominant_right_fixed_point(
        idmrg._transfer_matrices([ITensor((left, phys, right), A0)], 1))
    A0 = A0 / np.sqrt(eta0.real)
    B0 = ITensor((left, phys, right), A0)
    U_list, _eta = idmrg._canonicalize_periodic([B0], 1, cutoff=0.0, maxdim=D)
    return idmrg._to_array_lpr(U_list[0])


def _grow_initial_state(D, d_g, AL_old, AR_old):
    """(AL0, AR0, C0) at bond dimension D, warm-started from a smaller,
    already-(locally-)converged D_old<D solution (AL_old, AR_old) --
    mirrors `_random_initial_state`'s own two-direction canonicalization
    trick, but seeding each raw tensor with the old solution embedded in
    its top-left block (plus small noise on every entry, see
    `_random_raw_tensor`) rather than starting from pure noise.

    This is the standard, much more reliable way to reach a large target D
    in practice -- confirmed directly this matters, not just nice-to-have:
    single-shot VUMPS from a PURE random D=4/D=6 start reliably landed on
    self-consistent (gauge_mismatch~1e-10) fixed points with an energy
    density *worse* than the same model's own D=2 result (impossible for
    the true variational optimum, whose energy can only improve, not
    worsen, as D grows) -- i.e. genuinely stuck in a bad basin, not merely
    slow. Ramping D up one step at a time from an already-good smaller
    solution (see `vumps_ground_state`'s own docstring) avoids searching
    the full, much larger D-dimensional random landscape from scratch."""
    A0 = _random_raw_tensor(D, d_g, seed_AL=AL_old)
    AL0 = _canonicalize_raw(A0, D, d_g)
    A0_rev = _random_raw_tensor(D, d_g, seed_AL=np.transpose(AR_old, (2, 1, 0)).copy())
    AL0_rev = _canonicalize_raw(A0_rev, D, d_g)
    AR0 = np.transpose(AL0_rev, (2, 1, 0)).copy()
    C0 = np.eye(D, dtype=complex)
    return AL0, AR0, C0


def _energy_density_and_source_from_left(AL, W, r_AL):
    """(e, source_l): source_l = Source_L(AL) (the one-more-site-to-the-
    left Hamiltonian content, closed against AL's own exact left fixed
    point I -- NOT regularized/summed yet, that is `_solve_left_
    environment`'s own job), e = trace(conj(r_AL) @ source_l) -- the energy
    density this content implies, weighted by AL's own (generally non-
    trivial) right fixed point r_AL.

    The conjugate on r_AL is required, not optional -- found during the
    iDMRG D>1 excitation work (see idmrg_excitations.py's own module
    docstring) by cross-checking against an independently-converged
    MPSKit.jl D=2 TFIM state fed directly into this formula: plain
    `trace(r_AL @ source_l)` reproduced the wrong energy density (~1.72
    instead of the exact ~2.601), while `trace(conj(r_AL) @ source_l)`
    matched to machine precision. Derivable directly too: `source_l` lives
    in the *output* space of `apply_transfer_from_left`, and the correct
    closing functional for that space is the fixed point of that map's own
    adjoint, which (for AL's self-transfer) is `conj(r_AL)`, not `r_AL`
    itself -- `r_AL` is instead the fixed point of the *other* direction,
    `idmrg._apply_transfer`, a distinct map for any non-Hermitian transfer
    tensor. This is invisible whenever r_AL happens to be real (e.g. every
    D=1 case, where a trace-1 1x1 fixed point is always real, and D>1
    cases lucky enough to converge to a real gauge) -- which is why this
    survived undetected through this module's own D=1/D>1 energy-bound
    test suite (test_vumps.py's D>1 checks only assert loose bounds, not
    exact agreement) until the excitation work above needed exact D>1
    energies to cross-check a Hermiticity/dispersion computation against."""
    D = AL.shape[0]
    pending = idmrg_exc._pending_channels(W)
    h1 = idmrg_exc._onsite_matrix(W)
    I = np.eye(D, dtype=complex)

    def E_op(M):
        return idmrg_exc._op_transfer_matrix(AL, AL, M)

    source_l = idmrg._apply_transfer_from_left(E_op(h1), I)
    for mat_a, mat_b in pending:
        inner = idmrg._apply_transfer_from_left(E_op(mat_a), I)
        source_l = source_l + idmrg._apply_transfer_from_left(E_op(mat_b), inner)
    e = float(np.trace(r_AL.conj() @ source_l).real)
    return e, source_l


def _energy_density_and_source_from_right(AR, W, l_AR):
    """Mirror of `_energy_density_and_source_from_left`: Source_R(AR)
    (closed against AR's own exact right fixed point I), energy
    trace(conj(l_AR) @ source_r) -- l-weighted since here it is `r` (=I)
    that is trivial and `l_AR` that is the generic fixed point, the reverse
    of the left-side function above. See that function's own docstring for
    why the conjugate on l_AR is required, not optional -- the same bug,
    mirrored."""
    D = AR.shape[0]
    pending = idmrg_exc._pending_channels(W)
    h1 = idmrg_exc._onsite_matrix(W)
    I = np.eye(D, dtype=complex)

    def E_op(M):
        return idmrg_exc._op_transfer_matrix(AR, AR, M)

    source_r = idmrg._apply_transfer(E_op(h1), I)
    for mat_a, mat_b in pending:
        inner = idmrg._apply_transfer(E_op(mat_b), I)
        source_r = source_r + idmrg._apply_transfer(E_op(mat_a), inner)
    e = float(np.trace(l_AR.conj() @ source_r).real)
    return e, source_r


def _solve_left_environment(AL, W, r_AL, source_l, e):
    """GL: solves (I - T_AL^L + P^L)[GL] = source_l - e*I, T_AL^L[X] =
    apply_transfer_from_left(E_id(AL,AL), X), P^L[X] = I*trace(conj(r_AL)@X)
    -- the same style of regularized (I-T+projector) fixed-point solve
    idmrg_excitations.py's own `_channel_resolvent` uses, specialized to
    l=I (always exactly true for AL). The conjugate on r_AL is required for the
    same reason `_energy_density_and_source_from_left`'s own `e` needs it
    -- see that function's docstring.

    `e` MUST be AL's own e_L = trace(conj(r_AL)@source_l) (`_energy_
    density_and_source_from_left`'s own return value) -- NOT some other
    estimate (e.g. an average with AR's own e_R). Confirmed directly this
    is load-bearing, not a style choice: P^L is built so that
    (I-T_AL^L+P^L) is invertible despite T_AL^L's own dominant eigenvalue
    sitting exactly at 1 along r_AL, but the RESULT is only the correct,
    bounded fixed-point solution when the right-hand side source_l-e*I has
    zero projection along that same r_AL-weighted null direction
    (trace(conj(r_AL)@(source_l-e*I))=0, i.e. e=e_L exactly) -- otherwise
    the residual (e_L-e) leaks into GL as a spurious, uncontrolled
    component along the offending eigendirection. An earlier version of
    this module regularized GL/GR with a shared e=(e_L+e_R)/2 (reasoning
    that both should eventually agree at convergence, so a shared
    regularization target seemed natural) -- confirmed directly to be the
    cause of a persistent, exact period-2 limit cycle in the outer VUMPS
    iteration (energy alternating between two fixed values indefinitely,
    gauge_mismatch never decreasing below ~2, for a Heisenberg n_uc=1 test
    case) instead of convergence: e_L and e_R genuinely differ before
    AL/AR become mutually consistent, so the shared average fed GL a
    nonzero-null-component RHS every iteration, injecting exactly the
    runaway/oscillatory behavior described above. Using each side's own e
    (this function's current form) resolved it."""
    D = AL.shape[0]
    I = np.eye(D, dtype=complex)
    E_id = idmrg_exc._op_transfer_matrix(AL, AL, None)

    def left_action(X):
        return X - idmrg._apply_transfer_from_left(E_id, X) + I * np.trace(r_AL.conj() @ X)

    rhs = (source_l - e * I).reshape(-1)
    Mat = idmrg_exc._dense_linear_map(D, left_action)
    return np.linalg.solve(Mat, rhs).reshape(D, D)


def _solve_right_environment(AR, W, l_AR, source_r, e):
    """Mirror of `_solve_left_environment`: GR from AR (r=I always),
    regularized against AR's own fresh dominant left fixed point l_AR (the
    conjugate on l_AR is required, mirroring `_solve_left_environment`'s
    own r_AL -- see that function's docstring)."""
    D = AR.shape[0]
    I = np.eye(D, dtype=complex)
    E_id = idmrg_exc._op_transfer_matrix(AR, AR, None)

    def right_action(X):
        return X - idmrg._apply_transfer(E_id, X) + I * np.trace(l_AR.conj() @ X)

    rhs = (source_r - e * I).reshape(-1)
    Mat = idmrg_exc._dense_linear_map(D, right_action)
    return np.linalg.solve(Mat, rhs).reshape(D, D)


def _precompute_bond_environments(AL, AR, pending):
    """[(mat_a, mat_b, Lvec_a, Rvec_b), ...] -- one entry per reach-1 bond
    channel, Lvec_a/Rvec_b built once per outer VUMPS iteration (they only
    depend on AL/AR/W, not on the H_AC/H_C trial vector) and reused across
    every Lanczos matvec, mirroring idmrg_excitations._build_H_eff_dense's
    own "build the resolvent once, reuse across basis vectors" pattern.
    Lvec_a = apply_transfer_from_left(op_transfer_matrix(AL,AL,mat_a), I) --
    "AL with mat_a applied, closed from the left with AL's own exact left
    fixed point I" -- lives on the bond immediately to the left of AC/C.
    Rvec_b mirrors this on the right, via AR/apply_transfer/I."""
    D = AL.shape[0]
    I = np.eye(D, dtype=complex)
    out = []
    for mat_a, mat_b in pending:
        Lvec_a = idmrg._apply_transfer_from_left(
            idmrg_exc._op_transfer_matrix(AL, AL, mat_a), I)
        Rvec_b = idmrg._apply_transfer(
            idmrg_exc._op_transfer_matrix(AR, AR, mat_b), I)
        out.append((mat_a, mat_b, Lvec_a, Rvec_b))
    return out


def _h_ac_action(X, GL, GR, bond_envs, h1):
    """H_AC[X] -- see this module's own docstring for the diagram list.
    X: (D,d_g,D). Every term either leaves both of X's own bond legs open
    (onsite, background) or closes exactly one of them against a
    precomputed one-site-away environment (the two bond terms) -- unlike
    idmrg_excitations._h_eff_action, there is no momentum sum and no
    "diagram must close both sides" bookkeeping: X sits at a single, fixed
    position, so its own legs are simply the natural output legs of H_AC's
    result, never something to trace out."""
    Y = idmrg_exc._apply_op_ket(h1, X) if h1 is not None else np.zeros_like(X)
    Y = Y + idmrg_exc._cap_right(X, GR)
    Y = Y + idmrg_exc._cap_left(GL, X)
    for mat_a, mat_b, Lvec_a, Rvec_b in bond_envs:
        Y = Y + idmrg_exc._cap_right(idmrg_exc._apply_op_ket(mat_a, X), Rvec_b)
        Y = Y + idmrg_exc._cap_left(Lvec_a, idmrg_exc._apply_op_ket(mat_b, X))
    return Y


def _h_c_action(C, GL, GR, bond_envs):
    """H_C[C] -- the "zero-site" effective Hamiltonian, C sitting exactly
    on the bond between the AL-region and the AR-region. Reuses
    idmrg_excitations._cap_right/_cap_left directly on C reshaped to
    (D,1,D) (a dummy unit physical leg) rather than re-deriving the
    equivalent plain-matrix contraction by hand -- avoids re-deriving
    _cap_left's own (validated, deliberately non-symmetric-with-_cap_right)
    index convention a second time."""
    D = C.shape[0]
    C3 = C.reshape(D, 1, D)
    Y = idmrg_exc._cap_right(C3, GR) + idmrg_exc._cap_left(GL, C3)
    for _mat_a, _mat_b, Lvec_a, Rvec_b in bond_envs:
        Y = Y + idmrg_exc._cap_right(idmrg_exc._cap_left(Lvec_a, C3), Rvec_b)
    return Y.reshape(D, D)


def _update_AL_AR(AC, C):
    """New (AL, AR) from the just-solved (AC, C), via the standard
    orthogonal-Procrustes least-squares update (see this module's own
    docstring, step 5): AL minimizes ||AC - AL@C|| over isometric AL,
    solved by polar-decomposing AC (as a (D*d_g,D) matrix, left polar:
    AC=U_l@P_l, U_l isometric) and C (square, polar: C=U_c@P_c, U_c
    unitary) and setting AL=U_l@U_c^dagger; AR mirrors this via AC's own
    RIGHT polar factor (AC as a (D,d_g*D) matrix, right polar: AC=P_r@U_r,
    U_r isometric rows) and AR=U_c^dagger@U_r."""
    D, d_g, _ = AC.shape
    U_l, _P_l = _polar(AC.reshape(D * d_g, D), side='right')
    U_c, _P_c = _polar(C, side='right')
    AL_new = (U_l @ U_c.conj().T).reshape(D, d_g, D)

    U_r, _P_r = _polar(AC.reshape(D, d_g * D), side='left')
    AR_new = (U_c.conj().T @ U_r).reshape(D, d_g, D)
    return AL_new, AR_new


def _gauge_mismatch(AC, C, AL, AR):
    """||AC - AL@C|| + ||AC - C@AR||, normalized by ||AC|| -- 0 exactly at
    a true VUMPS fixed point (see this module's own docstring); the
    standard convergence diagnostic for single-site VUMPS.

    Must be called with the SAME (AL, AR) that were used to build the H_AC/
    H_C this AC/C were solved from (i.e. the pre-update gauge, not the
    fresh `_update_AL_AR(AC, C)` refit) -- confirmed directly this is not
    merely a style choice: `_update_AL_AR` constructs its own AL/AR as the
    isometric factor *closest* to this exact (AC, C) by orthogonal-
    Procrustes least squares, so comparing AC against THAT freshly-fit
    AL/AR makes this diagnostic close to 0 near-trivially (exactly 0 for
    D=1, where any nonzero scalar's polar decomposition is exact) --
    regardless of whether AC/C have actually converged to the Hamiltonian's
    own ground state. Comparing against the gauge that was actually held
    FIXED while solving for this AC/C has no such trivial floor: it is 0
    only once the Lanczos-solved AC/C genuinely agree with the AL/AR that
    produced their own environment."""
    lhs1 = np.einsum('lpm,mr->lpr', AL, C)
    lhs2 = np.einsum('lm,mpr->lpr', C, AR)
    norm_ac = np.linalg.norm(AC)
    if norm_ac == 0:
        return 0.0
    return (np.linalg.norm(AC - lhs1) + np.linalg.norm(AC - lhs2)) / norm_ac


def _environments(AL, AR, W, pending):
    """(GL, GR, e_cell, bond_envs) from the current (AL, AR) -- the full
    per-iteration environment build (steps 1-3 of this module's own
    docstring), factored out so `vumps_ground_state` can call it once per
    iteration AND once more after the loop exits, to refresh GL/GR/e_cell
    against the FINAL (post-loop) AL/AR before returning (see
    `vumps_ground_state`'s own comment for why this extra call is
    needed: e_cell/GL/GR computed at the top of the last executed
    iteration reflect that iteration's *input* AL/AR, one update step
    behind the AL/AR actually returned)."""
    E_AL = idmrg_exc._op_transfer_matrix(AL, AL, None)
    r_AL, _ = idmrg._dominant_right_fixed_point([E_AL])
    r_AL = (r_AL + r_AL.conj().T) / 2

    E_AR = idmrg_exc._op_transfer_matrix(AR, AR, None)
    l_AR, _ = idmrg._dominant_left_fixed_point([E_AR])
    l_AR = (l_AR + l_AR.conj().T) / 2

    e_L, source_l = _energy_density_and_source_from_left(AL, W, r_AL)
    e_R, source_r = _energy_density_and_source_from_right(AR, W, l_AR)

    # GL/GR are each regularized against THEIR OWN side's energy estimate
    # (e_L for GL, e_R for GR), not a shared average -- see
    # _solve_left_environment's own docstring for why using a mismatched
    # `e` is a real, confirmed bug (a persistent period-2 limit cycle),
    # not merely a style choice. e_cell (the *reported* energy density,
    # used only for VUMPSResult.e0/the outer loop's own printout) is
    # still their average -- a reasonable single-number summary once AL/AR
    # are close to mutually consistent, meaningless as a regularization
    # target beforehand.
    GL = _solve_left_environment(AL, W, r_AL, source_l, e_L)
    GR = _solve_right_environment(AR, W, l_AR, source_r, e_R)
    bond_envs = _precompute_bond_environments(AL, AR, pending)
    e_cell = 0.5 * (e_L + e_R)
    return GL, GR, e_cell, bond_envs


class VUMPSResult:
    """Converged (or best-effort, if not `.converged`) VUMPS ground state:
    the mixed-gauge tensors AL/AR/C/AC (D,d_g,D)/(D,d_g,D)/(D,D)/(D,d_g,D)
    numpy arrays, at the GROUPED supersite level (d_g = prod of the unit
    cell's own physical dimensions -- see this module's own docstring,
    "Scope"), plus the converged environments GL/GR and the grouped
    automaton W itself (all of which idmrg_excitations.py's own
    ExcitationEnvironment machinery would need directly, with no
    reconstruction step, to build a corrected mixed-gauge excitation
    ansatz on top of this result -- see this module's own docstring).

    `e0` is the ground-state energy density per PHYSICAL site (dividing the
    per-unit-cell `e_cell` by n_uc), matching idmrg.IDMRGResult.e0's own
    convention exactly, so the two are directly comparable."""

    def __init__(self, sites_uc, n_uc, D, d_g, AL, AR, C, AC, GL, GR, W,
                 e_cell, converged, niter_done, gauge_mismatch):
        self.sites_uc = sites_uc
        self.n_uc = n_uc
        self.D = D
        self.d_g = d_g
        self.AL = AL
        self.AR = AR
        self.C = C
        self.AC = AC
        self.GL = GL
        self.GR = GR
        self.W = W
        self.e_cell = e_cell
        self.e0 = e_cell / n_uc
        self.converged = converged
        self.niter_done = niter_done
        self.gauge_mismatch = gauge_mismatch


def vumps_ground_state(site_types, h_intra_op, h_inter_op, n_uc, D,
                        tol=1e-10, maxiter=800, niter_lanczos=30,
                        nrestarts=4, verbose=False):
    """Run single-site VUMPS to convergence (or `maxiter` iterations) for a
    unit cell of `n_uc` sites (type codes `site_types`), Hamiltonian
    h_intra_op/h_inter_op (plain MultiOperator.op-format term lists, exactly
    as idmrg.idmrg_ground_state itself takes them), at FIXED target bond
    dimension `D` (see this module's own docstring for why this, unlike
    idmrg_ground_state's `maxm`, is not merely a truncation cap but the
    actual dimension the variational optimum is sought at). Returns a
    VUMPSResult; `.e0` is the converged ground-state energy density per
    physical site.

    Ramps the bond dimension from D'=1 up to the target D one step at a
    time, warm-starting each D' from the PREVIOUS step's own best solution
    (`_grow_initial_state`) rather than starting the target D at a fully
    random guess -- confirmed directly this matters, not merely
    nice-to-have: single-shot VUMPS from a pure random start at D>1
    reliably lands on self-consistent (gauge_mismatch~1e-10) fixed points
    with an ENERGY WORSE than the same model's own smaller-D result
    (impossible for the true variational optimum, whose energy can only
    improve as D grows) -- i.e. genuinely stuck in a bad basin, not merely
    slow to converge, and this got WORSE (not better) at larger D from a
    random start. `nrestarts` attempts (one warm-started from the previous
    step, the rest random) are tried at each intermediate D' (capped at 3,
    warm-starting already does most of the work there) and, in full, at the
    final target D -- keeping the best across all of them (preferring
    converged over non-converged, then lowest `e0`; a mid-iteration AL/AR
    can transiently hit a (near-)degenerate transfer-matrix spectrum,
    idmrg._dominant_right_fixed_point's own guard -- confirmed directly on
    a gapless Heisenberg chain, whose own SU(2) symmetry makes this a real,
    reachable state along some trajectories -- such an attempt is simply
    skipped, not fatal). Additionally, at every D' a variational-principle
    safety net spends a bounded extra attempt budget if the best result
    found so far is worse than an already-known smaller-D' energy (see the
    inline comment at that check, and this module's own "Convergence
    robustness" docstring section for why this is needed, not just
    defensive). None of this guarantees the true optimum (no restart/ramp
    strategy does), so a caller doing quantitative work should still
    cross-check by increasing `nrestarts` and/or `D` and confirming the
    reported `e0` stops changing, the same convergence discipline
    idmrg_ground_state's own `maxm`/`etol` already require."""
    if n_uc > 2:
        raise NotImplementedError(
            "vumps_ground_state: n_uc>2 (got {}) is not supported yet -- "
            "see this module's own docstring / idmrg.py's module "
            "docstring for the same n_uc<=2 restriction there".format(n_uc))
    if D < 1:
        raise ValueError("vumps_ground_state: D must be >= 1, got {}".format(D))
    if nrestarts < 1:
        raise ValueError("vumps_ground_state: nrestarts must be >= 1, got {}".format(nrestarts))

    sites_uc, W_bulk = idmrg._build_automaton(h_intra_op, h_inter_op, site_types, n_uc)
    W = _group_automaton(W_bulk, n_uc)
    idmrg_exc._check_reach_one(W)
    d_g = int(np.prod([sites_uc.dim(p + 1) for p in range(n_uc)]))
    pending = idmrg_exc._pending_channels(W)
    h1 = idmrg_exc._onsite_matrix(W)

    def one_attempt(D_cur, init):
        try:
            result = _vumps_single_run(sites_uc, n_uc, D_cur, d_g, W, pending, h1,
                                        tol, maxiter, niter_lanczos, verbose, init=init)
        except RuntimeError as exc:
            if verbose:
                print("vumps D={} attempt: failed ({})".format(D_cur, exc))
            return None
        if verbose:
            print("vumps D={} attempt: e0={} converged={}".format(
                D_cur, result.e0, result.converged))
        return result

    def pick_better(a, b):
        """a if a is a better VUMPSResult than b (converged preferred,
        then lower e0), else b. Either may be None."""
        if a is None:
            return b
        if b is None:
            return a
        if a.converged != b.converged:
            return a if a.converged else b
        return a if a.e0 < b.e0 else b

    best = None
    prev_AL = prev_AR = None
    best_e0_so_far = None
    for D_cur in range(1, D + 1):
        n_here = nrestarts if D_cur == D else min(nrestarts, 3)
        local_best = None
        for attempt_i in range(n_here):
            init = (_grow_initial_state(D_cur, d_g, prev_AL, prev_AR)
                    if attempt_i == 0 and prev_AL is not None else None)
            local_best = pick_better(one_attempt(D_cur, init), local_best)
        if local_best is None:
            raise RuntimeError(
                "vumps_ground_state: every attempt at D={} failed (degenerate "
                "transfer-matrix spectrum) -- try increasing nrestarts".format(D_cur))
        # Safety net: the variational principle forbids a LARGER bond
        # dimension from ever doing worse than an already-known smaller-D
        # result -- confirmed directly this is a real, reachable failure
        # mode of the restart search above, not just a hypothetical
        # (an unguarded run at D=4 was observed to land at a wildly wrong,
        # even POSITIVE energy for a Hamiltonian whose D=1 energy was
        # already known to be ~-1.56), not merely "close but imprecise".
        # When this happens, spend a bounded extra budget of attempts
        # specifically trying to beat the known-achievable reference
        # before accepting the ramp's own result -- warm-started from the
        # SAME already-good (prev_AL, prev_AR) each time (fresh noise per
        # attempt, via `_grow_initial_state`'s own internal RNG call), not
        # fully random: confirmed directly that fully-random extra
        # attempts at D_cur can still fail to recover from a bad local_best
        # even given a generous extra budget, whereas repeatedly perturbing
        # the SAME known-good smaller-D anchor searches the region actually
        # likely to contain a good D_cur solution instead of resampling the
        # full, much less favorable random landscape from scratch every
        # time -- see this module's "Convergence robustness" docstring
        # section.
        if best_e0_so_far is not None and local_best.e0 > best_e0_so_far + 1e-6:
            for _ in range(2 * nrestarts):
                if local_best.e0 <= best_e0_so_far + 1e-6:
                    break
                init = (_grow_initial_state(D_cur, d_g, prev_AL, prev_AR)
                        if prev_AL is not None else None)
                local_best = pick_better(one_attempt(D_cur, init), local_best)
        prev_AL, prev_AR = local_best.AL, local_best.AR
        best_e0_so_far = (local_best.e0 if best_e0_so_far is None
                           else min(best_e0_so_far, local_best.e0))
        best = local_best
    return best


def _vumps_single_run(sites_uc, n_uc, D, d_g, W, pending, h1,
                       tol, maxiter, niter_lanczos, verbose, init=None):
    """One VUMPS attempt -- see `vumps_ground_state`'s own docstring for
    why this is wrapped in a multi-restart/D-ramp driver rather than
    called directly. `init`: an optional (AL,AR,C) starting point (see
    `_grow_initial_state`); defaults to a fresh `_random_initial_state`."""
    AL, AR, C = init if init is not None else _random_initial_state(D, d_g)
    AC = np.einsum('lpm,mr->lpr', AL, C)

    converged = False
    mismatch = None
    dim_ac = D * d_g * D
    it = 0
    for it in range(maxiter):
        GL, GR, e_cell, bond_envs = _environments(AL, AR, W, pending)

        def matvec_ac(x, GL=GL, GR=GR, bond_envs=bond_envs):
            X = x.reshape(D, d_g, D)
            return _h_ac_action(X, GL, GR, bond_envs, h1).reshape(-1)

        _lam_ac, ac_vec = _lanczos_ground_state(
            matvec_ac, AC.reshape(-1), niter=min(niter_lanczos, dim_ac))
        AC_new = ac_vec.reshape(D, d_g, D)

        def matvec_c(x, GL=GL, GR=GR, bond_envs=bond_envs):
            X = x.reshape(D, D)
            return _h_c_action(X, GL, GR, bond_envs).reshape(-1)

        _lam_c, c_vec = _lanczos_ground_state(
            matvec_c, C.reshape(-1), niter=min(niter_lanczos, D * D))
        C_new = c_vec.reshape(D, D)

        # Compared against the CURRENT (AL, AR) -- the gauge that was held
        # fixed while building GL/GR/bond_envs above and solving for this
        # AC_new/C_new -- not the freshly-refit `_update_AL_AR(AC_new,
        # C_new)` output, which would make this diagnostic near-trivially
        # small regardless of true convergence (see _gauge_mismatch's own
        # docstring).
        mismatch = _gauge_mismatch(AC_new, C_new, AL, AR)

        if verbose:
            print("vumps iter {}: e0={} gauge_mismatch={}".format(
                it, e_cell / n_uc, mismatch))

        AC, C = AC_new, C_new
        AL, AR = _update_AL_AR(AC, C)
        if mismatch < tol:
            converged = True
            break

    if not converged and verbose:
        print("vumps_ground_state: reached maxiter={} without converging "
              "to tol={} (last gauge_mismatch={})".format(maxiter, tol, mismatch))

    # AL/AR were just updated from the (AC, C) actually returned, but
    # GL/GR/e_cell above are still one step stale (built from the loop
    # iteration's OWN starting AL/AR, i.e. the pre-update gauge -- see
    # `_environments`' own docstring) -- refresh them once more here so
    # the returned VUMPSResult's GL/GR/e0 are consistent with its own
    # AL/AR/C/AC, not the previous iteration's.
    GL, GR, e_cell, _bond_envs = _environments(AL, AR, W, pending)

    return VUMPSResult(sites_uc, n_uc, D, d_g, AL, AR, C, AC, GL, GR, W,
                        e_cell, converged, it + 1, mismatch)
