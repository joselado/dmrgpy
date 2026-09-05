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
engine only. Static correlators (`onsite_expectation`/`two_point_correlator`
below) ARE implemented directly on a VUMPSResult, computed from the mixed-
gauge {AC, AR} rather than idmrg.py's dominant-right-fixed-point
eigenproblem -- see those functions' own docstrings, and
Infinite_Many_Body_Chain.vev/correlator for the public dispatch (also
ported to itensor_version=3's own gs_method="vumps" path,
mpscpp3/chain_session.h's Chain::vumps_onsite_expectation/
vumps_two_point_correlator -- a line-for-line C++ port of this module's
own formula, see that header's own doc comment). `imps_sum` (see this
module's own "Summing two converged VUMPS iMPS" section docstring further
below) is the VUMPS-mixed-gauge analogue of idmrg.py's own `imps_sum` --
pyitensor only, same physical scope limit.

== Convergence robustness (honest scope note) ==

Before any of the below: the single most important thing about this
loop's convergence is that its inner eigensolves run under
`_lanczos_ground_state`'s RESIDUAL criterion (`residual_tol=tol/10`), not
its default eigenVALUE one. `gauge_mismatch` compares AC and C, two
independently-solved eigenVECTORS, and a value-based stopping rule leaves
an eigenvector accurate only to the square root of its tolerance -- so
under the default rule this loop floored at a mismatch of ~1e-6 and could
not reach `tol=1e-10` at ANY number of outer iterations, while reporting a
perfectly good energy the whole time. Measured on a D=8 TFIM chain,
switching to the residual criterion took `converged` from 0/3 runs to 3/3
and the whole driver from 32.6s to 6.7s (g=1.5) and 38.4s to 4.6s
(g=1.0, critical), with the energy unchanged to every printed digit. See
`_lanczos_ground_state`'s own docstring. A second, smaller piece of the
same problem: AC and C being solved independently means nothing ties
their SIGNS together, and a flipped pair makes `gauge_mismatch` read
~4 instead of ~1e-6 -- 11 of every 100 iterations, measured -- so both
the convergence test and the reported diagnostic were partly noise. Both
are fixed in `_vumps_single_run`; the paragraphs below describe what
remains after that.

Plain single-attempt VUMPS from a random start turned out, during
development, to be genuinely unreliable for D>1 -- confirmed directly to
reach self-consistent (gauge_mismatch<tol) fixed points at an ENERGY WORSE
than the same model's own smaller-D result, which cannot happen for the
true variational optimum. `vumps_ground_state` mitigates this two ways:
(1) a D-ramp (warm-start each bond dimension from the previous one's own
best solution -- by subspace expansion, see `_subspace_expand` and the
section comment above it; noise-padding only where expansion cannot supply
the directions) plus multiple restarts at each step, and (2) a safety net
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
from scipy.linalg import null_space

from . import idmrg
from . import idmrg_excitations as idmrg_exc
from .dmrg import _lanczos_ground_state
from .index import Index
from .tensor import ITensor
from .tensor import noPrime as _t_noPrime


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
    matrix).

    `np.random.default_rng()` here is OS-entropy-seeded, so genuinely
    non-deterministic run to run, independent of any `np.random.seed(n)` a
    caller/test sets. That was long suspected (flagged 2026-08-07) of being
    the cause of a ~33% flake in
    `tests/test_infinite_chain.py::test_excitation_energies_matches_exact_xx_dispersion`,
    on the trivial fully-polarized D=1 case. Measured 2026-08-16, it is not:
    seeding is only why the failures land in different places run to run.

    The cause is conditioning in the *excitation* solver when the state is
    represented at a larger bond dimension than it needs. Over 40
    independent runs of that test's configuration (maxm=4 for a state whose
    exact representation is D=1), the ground-state energy was always
    excellent -- |e0 error| <= 1.8e-9 -- while the excitation error ranged
    from 8e-13 to 1.6e-3, an amplification of ~1e6. A bond-dimension sweep
    isolates it (12 runs each): maxm=1 and 2 never failed (worst 1.8e-15 and
    1.8e-10), maxm=4 failed 2/12, maxm=8 failed 10/12, with |e0 error| <=
    8e-10 at every maxm. Padding an exactly-D=1 state out to larger D leaves
    near-null directions in the bond space, and `idmrg_excitations`'
    (1-E)-type inverses (`_channel_resolvent` and the `GL`/`GR` builds) are
    ill-conditioned there.

    So a seeding fix would make the failures reproducible, not make them go
    away. The rest of that conclusion -- "and the ground-state solver is not
    the thing to change" -- was WRONG, and is kept here rather than deleted
    because the way it was wrong is instructive. The conditioning
    amplification is real, but what it was amplifying was convergence error
    in the ground state itself: `_lanczos_ground_state` stopped on the
    lowest Ritz VALUE, which caps the eigenVECTOR at the square root of that
    tolerance (~1e-6), and the excitation ansatz's inverses turned that into
    ~1e-3. The ground-state solver was exactly the thing to change -- not
    its seeding, but its stopping criterion (see that function's own
    docstring, and `_vumps_single_run`'s use of `residual_tol`). With the
    state converged to machine precision the same amplification lands at
    ~1e-13 and the phenomenon is gone: the regression that used to pin the
    degradation now pins its absence, as
    `tests/test_infinite_chain.py::
    test_redundant_bond_dimension_no_longer_costs_excitation_accuracy`.

    Two things survive from the original analysis. The measurement itself
    stands (padding an exactly-D=1 state does leave near-null bond
    directions, and the inverses do amplify them -- it is the size of what
    they amplify that changed). And projecting those directions out of the
    excitation ansatz's linear solves remains the principled fix if this
    ever matters for a physical, genuinely-entangled case, where the input
    to the amplification is a real truncation error rather than solver
    noise; nothing here reseeds anything to paper over it."""
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


# == Subspace expansion (bond-dimension growth) ==============================
#
# How the D-ramp gets from an already-converged solution at D_old to a
# starting point at D_new > D_old. `_grow_initial_state` (below, kept as
# the fallback) answers this with noise: embed the old tensors in the
# top-left block of a larger random one and re-canonicalize. That works,
# but the new directions it hands VUMPS are *arbitrary* -- the optimizer
# then has to discover, by itself, which of them the Hamiltonian actually
# wants populated, and the "Convergence robustness" section of this
# module's own docstring is largely a record of what happens when it
# fails to (a D=4 solve landing above its own D=2 energy, the restart
# budget and the variational safety net that exist to catch that).
#
# `_subspace_expand` answers it the way ITensorInfiniteMPS.jl's own
# `subspace_expansion.jl` does, and for the same reason: pick the new
# directions to be the ones H itself couples the current state to most
# strongly. Concretely -- with NL/NR orthonormal bases of the null spaces
# of AL and AR (the directions the current isometries do NOT reach), and
# theta = AL.C.AR the current two-site state --
#
#     M = NL^dagger . (H2 theta) . NR^dagger,   M = U S V^dagger (SVD)
#
# and the k largest singular pairs give the new left/right tensor columns
# NL.U / NR.V. Because M's row/column spaces are orthogonal to AL/AR by
# construction, the enlarged tensors are still exactly isometric, and
# because C is embedded with zeros in the new block, the expanded state is
# *the same physical state* as the D_old solution -- expansion adds
# variational directions without moving the energy, so the ramp's next
# step provably starts no worse than where the previous one finished. That
# is exactly what the noise start cannot promise.
#
# H2 here is the genuine two-site effective Hamiltonian, background
# environments included (GL, GR, both onsite terms, the bond term between
# the two sites, and the two bond terms straddling the boundary to the
# neighbours) -- i.e. `_h_ac_action`'s own diagram list widened from one
# site to two, built from the same primitives. Upstream's own
# `InfiniteSum{ITensor}` fast path drops the environments and keeps only
# the local two-site term; the environments are cheap here (they are
# already built and sitting on the VUMPSResult) and this is a *direction*
# heuristic either way, so there is no reason to drop them.
#
# Measured, at nrestarts=1 (which is what isolates the warm start: with
# the default restart budget the driver's cost is dominated by the RANDOM
# attempts, which expansion does not touch). Total outer VUMPS iterations
# to walk the whole ramp, noise start -> expansion, median over 3-5
# independent runs (independent, not seeded: `_random_raw_tensor` draws
# from an OS-seeded `default_rng`, so `np.random.seed` does not reach it):
#
#   TFIM g=1.5        D=4     72 ->  71     (0.05s -> 0.04s)
#   TFIM g=1.5        D=8    246 -> 238     (1.16s -> 1.20s)
#   TFIM g=1.0        D=4    146 -> 148     (0.18s -> 0.17s)
#   TFIM g=1.0        D=8    212 -> 232     (0.83s -> 1.09s)
#   Heisenberg        D=4   1601 -> 1601    (1.67s -> 1.49s)
#   Heisenberg        D=8   4001 -> 2401    (28.96s -> 6.23s)
#   Heisenberg n_uc=2 D=4   2400 -> 2400    (2.25s -> 2.01s)
#   Heisenberg n_uc=2 D=8   3200 -> 3200    (8.32s -> 7.38s)
#
# Read that honestly: expansion is a wash on almost every row, and the one
# row where it is decisive is the hardest one -- the gapless Heisenberg
# chain at D=8, 4.6x faster and, more to the point, 5x more accurate AND
# reproducible (the noise start's runs scattered between 2.1e-3 and 4.0e-2
# from the Bethe-ansatz energy density; every expansion run landed at
# 2.1e-3). That is the shape a better warm start should have: it cannot
# help a problem that was never hard, and it stops the hard one from
# depending on luck.
#
# These numbers replace a much more flattering earlier set (563 -> 127,
# 2401 -> 907, ...) measured before `_lanczos_ground_state` grew its Ritz
# RESIDUAL criterion. Most of that apparent gain was against a broken
# baseline: with an eigenvector only accurate to ~1e-6, a poor starting
# point cost many extra outer iterations to recover from, so *any*
# improvement to the start looked large. With the eigenvectors solved
# properly the starting point matters much less, and the honest benefit is
# the table above. The lesson is worth keeping: measure an optimization
# against a baseline you have already checked is not itself broken.


def _h_two_site_action(theta, GL, GR, bond_envs, h1):
    """H2[theta] for a (D,d_g,d_g,D) two-site tensor -- the two-site
    widening of `_h_ac_action`'s own term list (see the section comment
    above). The two physical legs are contracted with operators
    individually, and the two bond legs are closed against GL/GR or against
    the one-more-site bond environments exactly as in the one-site case.

    Kept separate from the theta it is applied to (rather than folded into
    `_subspace_expand`) so it can be checked on its own: assembled densely
    over a basis, this map must come out Hermitian, which is what pins the
    operator-ordering conventions of the three bond terms below -- the
    intra-pair one in particular applies mat_a and mat_b to two DIFFERENT
    legs, so a transposed factor there survives every isometry/state
    check the expansion itself is subject to and shows up only here."""
    D, d_g = theta.shape[0], theta.shape[1]

    def as3(T):
        return T.reshape(D, d_g * d_g, T.shape[3])

    def un3(T):
        return T.reshape(T.shape[0], d_g, d_g, T.shape[2])

    # background environments: both physical legs ride along untouched, so
    # the (d_g,d_g) pair is just one fat physical leg here.
    Y = un3(idmrg_exc._cap_left(GL, as3(theta)))
    Y = Y + un3(idmrg_exc._cap_right(as3(theta), GR))
    # the two onsite terms
    if h1 is not None:
        Y = Y + np.einsum('io,liqr->loqr', h1, theta)
        Y = Y + np.einsum('io,lpir->lpor', h1, theta)
    for mat_a, mat_b, Lvec_a, Rvec_b in bond_envs:
        # the bond BETWEEN the two sites: no environment at all, both ends
        # of the term land inside theta itself.
        Y = Y + np.einsum('io,jk,lijr->lokr', mat_a, mat_b, theta)
        # the bond straddling theta's LEFT edge: mat_a sits on the site to
        # the left (already summed into Lvec_a), mat_b on theta's own first
        # physical leg.
        Y = Y + un3(idmrg_exc._cap_left(
            Lvec_a, as3(np.einsum('io,liqr->loqr', mat_b, theta))))
        # ... and the mirror, straddling theta's RIGHT edge.
        Y = Y + un3(idmrg_exc._cap_right(
            as3(np.einsum('io,lpir->lpor', mat_a, theta)), Rvec_b))
    return Y


def _null_space_left(AL):
    """(D, d_g, kL) orthonormal basis of the directions AL's own isometry
    does not reach: columns orthogonal to AL's columns under the (l,p)
    inner product, i.e. sum_{l,p} conj(AL[l,p,m]) NL[l,p,a] = 0."""
    D, d_g, _ = AL.shape
    M = AL.reshape(D * d_g, -1)
    N = null_space(M.conj().T)
    return N.reshape(D, d_g, N.shape[1])


def _null_space_right(AR):
    """(kR, d_g, D) mirror of `_null_space_left`: rows orthogonal to AR's
    own rows, sum_{p,r} NR[a,p,r] conj(AR[l,p,r]) = 0. The conjugate sits
    on AR (not on NR) because AR's right-canonicality is itself written
    with the conjugate on the second factor -- so the complement is the
    null space of conj(AR), not of AR."""
    D, d_g, _ = AR.shape
    M = AR.reshape(D, d_g * D)
    N = null_space(M.conj())
    return N.T.reshape(N.shape[1], d_g, D)


def _subspace_expand(result, pending, h1, D_new, cutoff=1e-12):
    """(AL, AR, C) at bond dimension D_new, grown from a converged
    `VUMPSResult` at a smaller D by subspace expansion (see the section
    comment above). Returns None if the expansion cannot supply any new
    direction at all (D_new <= D, an exhausted null space, or an H2 that
    has no weight outside the current subspace -- the last of which means
    the state is already exact at this D and there is nothing to grow
    into); the caller falls back to `_grow_initial_state` there.

    The returned dimension is D + k with k <= D_new - D: the null spaces
    have only D*(d_g-1) directions each and the SVD cutoff may discard
    more, so the caller must be prepared for a state SMALLER than it asked
    for and pad the remainder itself. It is never larger."""
    AL, AR, C = result.AL, result.AR, result.C
    D, d_g, _ = AL.shape
    k_want = D_new - D
    if k_want <= 0:
        return None

    NL = _null_space_left(AL)
    NR = _null_space_right(AR)
    if NL.shape[2] == 0 or NR.shape[0] == 0:
        return None

    bond_envs = _precompute_bond_environments(AL, AR, pending)
    theta = np.einsum('lpm,mn,nqr->lpqr', AL, C, AR)
    H2theta = _h_two_site_action(theta, result.GL, result.GR, bond_envs, h1)

    # <NL_a NR_b| H2 |theta>: both bras conjugate, hence conj on both
    # null-space factors.
    M = np.einsum('lpa,lpqr,bqr->ab', NL.conj(), H2theta, NR.conj())
    U, S, Vh = np.linalg.svd(M, full_matrices=False)
    keep = int(np.sum(S > cutoff * max(S[0], 1.0))) if S.size else 0
    keep = min(keep, k_want)
    if keep == 0:
        return None

    AL_add = np.einsum('lpa,ak->lpk', NL, U[:, :keep])
    AR_add = np.einsum('kb,bqr->kqr', Vh[:keep, :], NR)

    Dn = D + keep
    AL_new = np.zeros((Dn, d_g, Dn), dtype=complex)
    AL_new[:D, :, :D] = AL
    AL_new[:D, :, D:] = AL_add
    AR_new = np.zeros((Dn, d_g, Dn), dtype=complex)
    AR_new[:D, :, :D] = AR
    AR_new[D:, :, :D] = AR_add
    C_new = np.zeros((Dn, Dn), dtype=complex)
    C_new[:D, :D] = C
    return AL_new, AR_new, C_new


def _grow_initial_state(D, d_g, AL_old, AR_old):
    """(AL0, AR0, C0) at bond dimension D, warm-started from a smaller,
    already-(locally-)converged D_old<D solution (AL_old, AR_old) --
    mirrors `_random_initial_state`'s own two-direction canonicalization
    trick, but seeding each raw tensor with the old solution embedded in
    its top-left block (plus small noise on every entry, see
    `_random_raw_tensor`) rather than starting from pure noise.

    This is no longer the ramp's *first* choice: `_subspace_expand` (see
    the section comment above) picks the new directions from H itself
    rather than from noise; it is a wash on most models and decisive on
    the hardest one tested (see that comment's own table). What is left
    for this function is the two places expansion
    cannot serve -- padding out a ramp step the null space was too small to
    fill, and the variational safety net's extra attempts, which need a
    FRESH perturbation of the same anchor each time and so cannot use a
    deterministic expansion. Both remain load-bearing.

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

    return idmrg_exc._solve_linear_map(D, left_action, source_l - e * I)


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

    return idmrg_exc._solve_linear_map(D, right_action, source_r - e * I)


def _precompute_bond_environments(AL, AR, pending):
    """[(mat_a, mat_b, Lvec_a, Rvec_b), ...] -- one entry per reach-1 bond
    channel, Lvec_a/Rvec_b built once per outer VUMPS iteration (they only
    depend on AL/AR/W, not on the H_AC/H_C trial vector) and reused across
    every Lanczos matvec, mirroring idmrg_excitations._build_H_eff_dense's
    own "build the resolvent once, reuse across basis vectors" pattern.
    Lvec_a = apply_transfer_from_left(op_transfer_matrix(AL,AL,mat_a), I) --
    "AL with mat_a applied, closed from the left with AL's own exact left
    fixed point I" -- lives on the bond immediately to the left of AC/C.
    Rvec_b mirrors this on the right, via AR/apply_transfer/I.

    Both closures are done *without* ever materializing the (D,D,D,D)
    transfer matrix the formulas above name: closing E4 against the
    identity is just a trace over one of its index pairs, i.e. a direct
    contraction of the two rank-3 tensors over two legs
    (O(D^2 d D) work and O(D^2) memory, versus O(D^4) for both if E4 is
    built first). This runs once per bond channel per macro-iteration,
    and at D=30 the intermediate it avoids is 810k complex entries per
    bond."""
    out = []
    for mat_a, mat_b in pending:
        # apply_transfer_from_left(op_transfer_matrix(AL,AL,mat_a), I)
        #   = sum_l E4[l,l,r,R] = sum_{l,p} (mat_a AL)[l,p,r] conj(AL)[l,p,R]
        AL_op = idmrg_exc._apply_op_ket(mat_a, AL)
        Lvec_a = np.tensordot(AL_op, np.conj(AL), axes=([0, 1], [0, 1]))
        # apply_transfer(op_transfer_matrix(AR,AR,mat_b), I)
        #   = sum_r E4[l,L,r,r] = sum_{p,r} (mat_b AR)[l,p,r] conj(AR)[L,p,r]
        AR_op = idmrg_exc._apply_op_ket(mat_b, AR)
        Rvec_b = np.tensordot(AR_op, np.conj(AR), axes=([1, 2], [1, 2]))
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


def _multisite_ground_state(sites_uc, W_bulk, n_uc, D, tol, maxiter,
                             niter_lanczos, nrestarts, verbose):
    """`vumps_ground_state` for a cell too big to group: run
    `vumps_ms.ground_state` and wrap its result in the same `VUMPSResult`
    the grouped path returns.

    The grouped fields (`d_g`, the single `AL`/`AR`/`C`/`AC`, the grouped
    `W`) have no single-tensor meaning here -- the state IS a list of
    per-site tensors -- so they are filled with the per-site lists and
    `d_g` reports the product it *would* have had. Consumers that need the
    grouped form (`idmrg_excitations`' tangent-space ansatz, which builds
    on a single supersite) therefore have to be taught the multi-site form
    before they can be used at n_uc>2; `excitation_energies` raises for
    that case rather than silently reading a list as a tensor."""
    from . import vumps_ms
    W_list = [W_bulk[p].array for p in range(n_uc)]
    dims = [sites_uc.dim(p + 1) for p in range(n_uc)]
    out = vumps_ms.ground_state(W_list, dims, D, tol=tol, maxiter=maxiter,
                                 niter_lanczos=niter_lanczos,
                                 nrestarts=nrestarts, verbose=verbose)
    d_g = int(np.prod(dims))
    res = VUMPSResult(sites_uc, n_uc, D, d_g, out["AL"], out["AR"], out["C"],
                       out["AC"], out["GL"], out["GR"], W_list,
                       out["e_cell"], out["converged"], out["niter"],
                       out["mismatch"])
    res.multisite = True
    return res


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
    if D < 1:
        raise ValueError("vumps_ground_state: D must be >= 1, got {}".format(D))
    if nrestarts < 1:
        raise ValueError("vumps_ground_state: nrestarts must be >= 1, got {}".format(nrestarts))

    sites_uc, W_bulk = idmrg._build_automaton(h_intra_op, h_inter_op, site_types, n_uc)

    # n_uc > 2 goes to the sequential multi-site algorithm, which never
    # groups the cell and so costs LINEARLY rather than exponentially in
    # n_uc (see pyitensor/vumps_ms.py's own module docstring, and
    # arXiv:2003.01142 for why the grouped route is the wrong one to scale
    # up). n_uc <= 2 stays on the grouped path below purely so the values
    # this module has been validated against do not move; the two agree to
    # machine precision where both apply, which is what
    # tests/test_vumps_multisite.py checks.
    if n_uc > 2:
        return _multisite_ground_state(
            sites_uc, W_bulk, n_uc, D, tol, maxiter, niter_lanczos,
            nrestarts, verbose)

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

    def warm_start(prev, D_cur):
        """The ramp's first attempt at D_cur, grown from the previous
        step's own best result. Subspace expansion when it can supply the
        directions (see `_subspace_expand`), noise-padding for whatever it
        cannot -- an exhausted null space or an SVD that keeps fewer
        directions than the ramp step asks for both leave a state smaller
        than D_cur, which `_grow_initial_state` then pads the rest of the
        way. Never worse than the pure-noise start it replaces: the
        fallback IS that start."""
        expanded = _subspace_expand(prev, pending, h1, D_cur)
        if expanded is None:
            return _grow_initial_state(D_cur, d_g, prev.AL, prev.AR)
        AL_e, AR_e, C_e = expanded
        if AL_e.shape[0] < D_cur:
            return _grow_initial_state(D_cur, d_g, AL_e, AR_e)
        return AL_e, AR_e, C_e

    best = None
    prev_result = None
    prev_AL = prev_AR = None
    best_e0_so_far = None
    for D_cur in _d_ramp(D):
        n_here = nrestarts if D_cur == D else min(nrestarts, 3)
        local_best = None
        for attempt_i in range(n_here):
            init = (warm_start(prev_result, D_cur)
                    if attempt_i == 0 and prev_result is not None else None)
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
        prev_result = local_best
        prev_AL, prev_AR = local_best.AL, local_best.AR
        best_e0_so_far = (local_best.e0 if best_e0_so_far is None
                           else min(best_e0_so_far, local_best.e0))
        best = local_best
    return best


def _d_ramp(D):
    """The bond dimensions the driver above actually solves at on its way
    to D: 1, 2, 4, 8, ... , D (doubling, always ending exactly at D).

    Every step of this ramp is a full multi-restart VUMPS solve, and one
    solve costs roughly O(D^3), so ramping one integer at a time (as this
    function's predecessor did) makes the whole driver O(D^4) -- the ramp
    itself, not the solve at the target D, then dominates: 32 complete
    solves to reach D=8, and ~100 to reach D=30. Doubling reaches D=30 in
    6 steps whose combined cost is a small multiple of the last one's.

    What the ramp is *for* is unchanged (see this module's "Convergence
    robustness" docstring section): a pure random start at large D lands
    in a bad basin often enough to matter, so each step is warm-started
    from the previous one via `_grow_initial_state`, which embeds the old
    solution in the new tensor's top-left block plus noise and assumes
    nothing about how big the jump is."""
    ramp = []
    d = 1
    while d < D:
        ramp.append(d)
        d *= 2
    ramp.append(D)
    return ramp


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

        # residual_tol, not the default eigenVALUE test: the criterion
        # below is a norm difference between AC and C, two independently
        # solved eigenvectors, and the eigenvalue test caps eigenVECTOR
        # accuracy at ~sqrt(tol) -- which is why this loop used to floor
        # at a gauge mismatch of ~1e-6 and could not reach tol=1e-10 at
        # any number of iterations. See `_lanczos_ground_state`'s own
        # docstring for the measurement.
        _lam_ac, ac_vec = _lanczos_ground_state(
            matvec_ac, AC.reshape(-1), niter=min(niter_lanczos, dim_ac),
            residual_tol=tol / 10.0)
        AC_new = ac_vec.reshape(D, d_g, D)

        def matvec_c(x, GL=GL, GR=GR, bond_envs=bond_envs):
            X = x.reshape(D, D)
            return _h_c_action(X, GL, GR, bond_envs).reshape(-1)

        _lam_c, c_vec = _lanczos_ground_state(
            matvec_c, C.reshape(-1), niter=min(niter_lanczos, D * D),
            residual_tol=tol / 10.0)
        C_new = c_vec.reshape(D, D)

        # An eigenvector is only defined up to a phase, and AC and C were
        # just solved INDEPENDENTLY -- so nothing ties their phases
        # together, while `_gauge_mismatch` below compares them directly.
        # The Lanczos tridiagonal is real symmetric, so the freedom is a
        # sign rather than a general phase, and each solve is warm-started
        # from the previous iteration's own vector: aligning each new
        # vector with the one it started from is enough to keep the pair
        # consistent by induction (AC = AL@C exactly, at iteration 0).
        # Without this, a sign flip makes ||AC-AL@C|| read ~2||AC||
        # instead of ~0 -- measured on a D=8 TFIM chain, 11 of every 100
        # iterations reported a mismatch of 4.0 for a state whose true
        # mismatch was ~1e-6, so both the convergence test and the
        # `gauge_mismatch` on the returned VUMPSResult were a coin flip.
        if np.vdot(AC, AC_new).real < 0:
            AC_new = -AC_new
        if np.vdot(C, C_new).real < 0:
            C_new = -C_new

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


# == Static correlators ======================================================
#
# Computed directly from the converged mixed-gauge {AC, AR} rather than
# idmrg.py's dominant-right-fixed-point eigenproblem -- AC is EXACTLY the
# correctly normalized single-(super)site reduced state by construction of
# the mixed canonical gauge (Vanderstraeten, Haegeman, Verstraete,
# arXiv:1810.07006, Eq.(34): "locating the center site where the operator is
# acting ... everything to the left and right is contracted to the
# identity"), and AR is exactly right-orthonormal (sum_p AR_p AR_p^dagger =
# I), so both the left closure (via AC) and the right closure (via AR, for
# any operator strictly to the right of AC's own cell) are algebraically
# exact rather than needing a numerical eigen-solve the way idmrg.py's
# one-sided-isometric U_list does -- see `two_point_correlator`'s own
# docstring for the derivation (mirroring arXiv:1810.07006 Eq.(37)-(39)'s
# "uniform gauge" correlator formula, specialized to the mixed gauge).


def _embed_group_operator(sites_uc, n_uc, ops_by_pos):
    """The d_g x d_g grouped-supersite operator matrix (M[i,o] convention,
    same as idmrg._op_transfer/_onsite_matrix) obtained by placing each
    `ops_by_pos[p]` (a dense d_p x d_p M[i,o] matrix) at sub-site p of the
    n_uc-site unit cell and identity everywhere else, combined via
    `np.kron` in the SAME sequential (site 0 slowest-varying, site n_uc-1
    fastest-varying) order `_group_automaton` itself uses to build the
    grouped physical index -- confirmed directly against
    `_group_automaton`'s own `np.einsum('Labm,mcdR->LacbdR', ...)` +
    reshape construction, which is exactly what `np.kron`'s row-major
    flattening produces for a chain of per-site (in,out) matrices.
    `ops_by_pos` may name the SAME position at most once -- composing two
    operators at one physical site (e.g. `two_point_correlator`'s own r=0
    case) is the caller's job (an ordinary d_p x d_p matrix product
    before calling this), not something this function does."""
    mats = []
    for p in range(n_uc):
        if p in ops_by_pos:
            mats.append(ops_by_pos[p])
        else:
            d = sites_uc.dim(p + 1)
            mats.append(np.eye(d, dtype=complex))
    M = mats[0]
    for m in mats[1:]:
        M = np.kron(M, m)
    return M


def onsite_expectation(result, opname, p):
    """<opname> at sub-site p (0..n_uc-1) of the converged VUMPSResult's
    unit cell -- mirrors idmrg.onsite_expectation's own public signature
    and semantics, but computed directly from AC (arXiv:1810.07006 Eq.(34))
    rather than idmrg.py's dominant-right-fixed-point eigenproblem: AC
    already IS the correctly normalized single-(super)site reduced state
    by construction of the mixed canonical gauge, so no eigenproblem is
    needed here at all."""
    if not (0 <= p < result.n_uc):
        raise ValueError("onsite_expectation: p must be in 0..{} (n_uc-1), "
                          "got {!r}".format(result.n_uc - 1, p))
    if getattr(result, "multisite", False):
        from . import vumps_ms
        return vumps_ms.onsite_expectation(result.AC, result.sites_uc, opname, p)
    M = _embed_group_operator(
        result.sites_uc, result.n_uc,
        {p: result.sites_uc.site_type(p + 1).matrix(opname)})
    AC = result.AC
    AC_op = np.einsum('io,lir->lor', M, AC)
    val = np.einsum('lor,lor->', np.conj(AC), AC_op)
    norm = np.einsum('lir,lir->', np.conj(AC), AC).real
    return complex(val / norm)


def two_point_correlator(result, opname_i, p_i, opname_j, r):
    """<opname_i(site p_i) opname_j(site p_i + r)> of the converged
    VUMPSResult's infinite chain, r measured in physical sites (r>=0) --
    mirrors idmrg.two_point_correlator's own signature, r=0 same-site
    convention (M_j @ M_i, see that function's own docstring for why this
    order, not the reverse, and that this is only well-defined when
    opname_i/opname_j don't commute), and n_uc-periodicity, but built from
    the mixed-gauge {AC, AR} instead of idmrg.py's growing-algorithm
    dominant-fixed-point machinery.

    Places AC (the correctly normalized center) at the unit cell containing
    p_i and AR (right-orthonormal, `sum_p AR_p AR_p^dagger = I`) at every
    cell strictly to its right -- a valid mixed-canonical-gauge choice at
    ANY cut position, since AL@C=C@AR=AC (this module's own docstring).
    When both operators land in that same AC cell (r spans less than one
    unit cell), this reduces to onsite_expectation's own full-AC-
    contraction formula (Eq.(34)); otherwise AC's own right bond is left
    open, propagated through zero or more plain AR transfer tensors, has
    the second operator inserted via one more AR transfer tensor, and is
    closed by a direct trace -- exploiting AR's exact right-orthonormality
    to skip the dominant-right-fixed-point eigenproblem idmrg.
    two_point_correlator needs (confirmed algebraically: closing any number
    of further plain-AR transfer tensors after the last operator leaves the
    trace exactly invariant, since `sum_r E4_AR[l,L,r,r] = delta[l,L]` IS
    AR's own right-orthonormality condition -- the mixed-gauge analogue of
    arXiv:1810.07006 Eq.(37)-(39)'s "uniform gauge" correlator formula,
    which instead needs an explicit right fixed point `r` because a single
    uniform-gauge tensor A is not already right-canonical)."""
    if r < 0:
        raise ValueError("two_point_correlator: r must be >= 0")
    n_uc = result.n_uc
    if not (0 <= p_i < n_uc):
        raise ValueError("two_point_correlator: p_i must be in 0..{} "
                          "(n_uc-1), got {!r}".format(n_uc - 1, p_i))
    if getattr(result, "multisite", False):
        from . import vumps_ms
        return vumps_ms.two_point_correlator(
            result.AC, result.AR, result.sites_uc, n_uc,
            opname_i, p_i, opname_j, r)
    sites_uc = result.sites_uc
    AC = result.AC
    AR = result.AR
    norm = np.einsum('lir,lir->', np.conj(AC), AC).real

    cell_offset, p_j = divmod(p_i + r, n_uc)

    # Endpoint matrices, ordering sign, and whether a Jordan-Wigner string
    # is open between the two operators -- from the same helper that builds
    # the Hamiltonian's own 2-site terms, so the two cannot disagree about
    # the convention. Parity-even operators come back stringless and with
    # no endpoint factor, reproducing this function's pre-string behaviour
    # exactly. See idmrg.two_point_correlator's docstring.
    term = [1.0, [opname_i, p_i], [opname_j, p_i + r]]
    _rel, coef, mats, ferm = idmrg._term_site_matrices(term, sites_uc, n_uc)
    strung = ferm[0] if r > 0 else False

    def _F(p):
        return sites_uc.site_type(p + 1).matrix("F")

    def _string_ops(lo, hi):
        """{p: F} for every sub-site p with lo <= p < hi (empty if not
        strung) -- the part of the string living inside one unit cell."""
        return {p: _F(p) for p in range(lo, hi)} if strung else {}

    if cell_offset == 0:
        if p_j == p_i:
            M = _embed_group_operator(sites_uc, n_uc, {p_i: mats[0]})
        else:
            ops = {p_i: mats[0], p_j: mats[1]}
            ops.update(_string_ops(p_i + 1, p_j))   # string strictly between
            M = _embed_group_operator(sites_uc, n_uc, ops)
        AC_op = np.einsum('io,lir->lor', M, AC)
        val = np.einsum('lor,lor->', np.conj(AC), AC_op)
        return complex(coef * val / norm)

    # The string, when open, runs from just right of p_i to the end of the
    # AC cell, across every fully-crossed intermediate cell, and from the
    # start of the final cell up to just left of p_j.
    ops_i = {p_i: mats[0]}
    ops_i.update(_string_ops(p_i + 1, n_uc))
    Mi_embed = _embed_group_operator(sites_uc, n_uc, ops_i)
    AC_op = np.einsum('io,lir->lor', Mi_embed, AC)
    # Open right-bond object: bra/ket both AC, operator on the ket side
    # only -- this already IS the full left closure (AC's own left leg is
    # summed away here), leaving just the (ket-bond, bra-bond) legs open.
    X = np.einsum('lor,loR->rR', AC_op, np.conj(AC)) / norm

    if cell_offset > 1:
        # Every fully-crossed cell carries F on ALL of its sub-sites when
        # the string is open, and the plain transfer tensor otherwise.
        M_cross = (_embed_group_operator(sites_uc, n_uc, _string_ops(0, n_uc))
                   if strung else None)
        AR_cross = AR if M_cross is None else np.einsum('io,lir->lor', M_cross, AR)
        E_AR = np.einsum('lpr,LpR->lLrR', AR_cross, np.conj(AR))
        for _ in range(cell_offset - 1):
            X = idmrg._apply_transfer_from_left(E_AR, X)

    ops_j = {p_j: mats[1]}
    ops_j.update(_string_ops(0, p_j))               # string up to p_j
    Mj_embed = _embed_group_operator(sites_uc, n_uc, ops_j)
    AR_op = np.einsum('io,lir->lor', Mj_embed, AR)
    E_AR_op = np.einsum('lpr,LpR->lLrR', AR_op, np.conj(AR))
    X = idmrg._apply_transfer_from_left(E_AR_op, X)

    return complex(coef * np.trace(X))


# == Summing two converged VUMPS iMPS ========================================
#
# imps_sum(result_a, result_b) is the VUMPS/mixed-gauge analogue of
# idmrg.imps_sum -- see that function's own module-level docstring in
# idmrg.py ("Summing two converged iMPS" section) for the full physical
# derivation, which applies here unchanged. The construction:
#
# 1. Block-diagonal direct sum of result_a.AL and result_b.AL (grouped-
#    supersite tensors, bond dimension D_a+D_b) -- both are already
#    individually left-canonical (isometric) by construction of any
#    converged VUMPSResult, so this raw direct sum is ALREADY exactly
#    left-canonical too (a block-diagonal direct sum of two isometries is
#    itself an isometry: sum_p AL_sum_p^dagger AL_sum_p = block_diag(sum_p
#    AL_a_p^dagger AL_a_p, sum_p AL_b_p^dagger AL_b_p) = block_diag(I,I) =
#    I) -- unlike idmrg._periodic_direct_sum's own per-sublattice list
#    construction, no per-position loop is needed here: VUMPS already
#    works at the single grouped-supersite level (n_uc sites folded into
#    one d_g-dimensional site via `_group_automaton`), so there is only
#    ever one cut to sum across.
# 2. Re-canonicalize/truncate via `idmrg._canonicalize_periodic` (the same
#    two-sided fixed-point procedure idmrg.py's own apply_mpo/imps_sum
#    already use), called here on a trivial n_uc=1 "periodic chain" --
#    reused rather than skipped even though step 1's raw tensor is already
#    left-canonical, because this is also what applies the caller's
#    requested (cutoff, maxdim) truncation via a genuine two-sided
#    fixed-point SVD (not a naive per-block truncation), and, crucially,
#    is where the SAME degeneracy check idmrg.imps_sum relies on
#    (`idmrg._dominant_right_fixed_point`, called internally) fires for
#    the "two ordinary/tied-norm branches" case -- see below.
# 3. Complete the resulting left-canonical AL to the full mixed gauge
#    {AL, AR, C, AC} via `_complete_mixed_gauge` -- the standard "bringing
#    a uniform MPS to canonical form" construction (Vanderstraeten,
#    Haegeman, Verstraete, "Tangent-space methods for uniform matrix
#    product states", arXiv:1810.07006, Sec. 2.1, Eq. (9)-(17), specialized
#    to an already-left-canonical input): factor AL's own dominant right
#    transfer-matrix fixed point r = C C^dagger, then AR := C^-1 AL C and
#    AC := AL C (== C AR algebraically, by construction of AR -- not merely
#    approximately). This is NOT a further VUMPS energy-minimization step
#    (imps_sum's output is generally not an eigenstate of anything) --
#    purely bookkeeping to re-derive a valid, self-consistent mixed gauge
#    for whatever tensor step 2 produced, exactly as `_random_initial_
#    state`'s own two-direction canonicalization trick does for a fresh
#    random tensor (see that function's own docstring), except here C is
#    solved for properly (via the fixed-point equation) rather than left at
#    the identity, since there is no follow-up VUMPS iteration here to fix
#    up an inconsistent C the way there is for a fresh VUMPS run.
#
# PHYSICAL SCOPE -- identical to idmrg.imps_sum's own scope note (see
# idmrg.py's "Summing two converged iMPS" section docstring for the full
# derivation). Every converged VUMPSResult has AL exactly left-canonical
# AND AR exactly right-canonical by construction of the mixed gauge, so
# BOTH its left and right transfer eigenvalues are exactly 1 (not merely
# close) -- summing two ordinary VUMPSResults therefore produces a
# combined transfer matrix with a genuinely 2-fold degenerate dominant
# eigenvalue (block (a,a) and block (b,b) each contribute eigenvalue 1;
# see idmrg.py's `_dominant_eigenvalue_mixed`/`imps_overlap` machinery for
# why the cross (a,b) block generically has magnitude <1 instead, short of
# a genuine gauge-equivalence between the two states), and
# `idmrg._dominant_right_fixed_point`'s own degeneracy check (reused here
# unmodified, inside step 2's `_canonicalize_periodic` call) reliably
# raises RuntimeError there rather than silently collapsing to one
# arbitrary branch -- confirmed directly (see this module's test suite):
# the same "orthogonality catastrophe" arXiv:1810.07006 Eq.(22)-(24) and
# its own Sec. 2.1 remark on non-injective/"cat state" MPS describe in the
# uniform-MPS language directly. Only two states with a genuine per-site
# norm mismatch (e.g. one deliberately rescaled -- ordinary VUMPSResults
# never carry this on their own) have a well-posed sum, exactly mirroring
# idmrg.imps_sum's own worked example.


def _complete_mixed_gauge(AL):
    """(AR, C, AC): completes a left-canonical grouped-supersite tensor AL
    (D,d_g,D) to the full VUMPS mixed gauge -- the standard "bring a
    uniform MPS to canonical form" construction (Vanderstraeten, Haegeman,
    Verstraete, arXiv:1810.07006, Eq. (9)-(17)), specialized to an
    already-left-canonical input (that reference's own transform "L",
    which maps a general raw tensor to left-canonical form, is already the
    identity here): factor AL's own dominant right transfer-matrix fixed
    point r = C C^dagger (Hermitian PSD square root via `eigh`, mirroring
    `idmrg._psd_sqrt_factor`'s own reasoning for Hermitizing before the
    square root -- r is guaranteed Hermitian PSD in theory, but
    `idmrg._dominant_right_fixed_point`'s general (non-symmetric)
    eigensolver only returns it Hermitian up to numerical noise), then
    AR := C^-1 @ AL @ C and AC := AL @ C -- AC == C @ AR is then an exact
    algebraic identity (not merely approximate), since AR is defined in
    terms of C in the first place: C @ AR = C @ C^-1 @ AL @ C = AL @ C.

    AL's dominant right eigenvalue is guaranteed to equal exactly 1 (up to
    numerical noise), not merely close to it: AL is left-canonical by
    construction (sum_p AL_p^dagger AL_p = I exactly), which is precisely
    the statement that I is the transfer matrix's own LEFT fixed point
    with eigenvalue exactly 1 -- and a matrix's spectrum is identical
    whether obtained via its left or right eigenvectors (same operator, T
    and T^T share eigenvalues), so this is not an extra assumption, just
    that fact reused. `idmrg._dominant_right_fixed_point` raises
    RuntimeError if AL's own dominant eigenvalue is (near-)degenerate --
    see this module's own "Summing two converged VUMPS iMPS" section
    docstring above for why that is the expected, correct outcome when AL
    is itself already a block-diagonal direct sum of two equally-
    normalized branches, not a bug to route around here."""
    D, d_g, _ = AL.shape
    left, right, phys = (Index(D, tags="Link"), Index(D, tags="Link"),
                          Index(d_g, tags="Site"))
    Es = idmrg._transfer_matrices([ITensor((left, phys, right), AL)], 1)
    r, _eta = idmrg._dominant_right_fixed_point(Es)
    herm = (r + r.conj().T) / 2
    evals, evecs = np.linalg.eigh(herm)
    evals = np.clip(evals.real, 0.0, None)
    C = (evecs * np.sqrt(evals)[None, :]) @ evecs.conj().T
    Cinv = np.linalg.pinv(C)
    AR = np.einsum('ab,bpc,cd->apd', Cinv, AL, C)
    AC = np.einsum('apb,bc->apc', AL, C)
    return AR, C, AC


class UniformMPS:
    """A converged/summed VUMPS-mixed-gauge uniform iMPS with no
    ground-state-specific bookkeeping -- same shape as VUMPSResult
    (sites_uc, n_uc, D, d_g, AL, AR, C, AC) minus e0/GL/GR/W/converged/
    niter_done/gauge_mismatch, which have no meaning for imps_sum's output
    (not a Hamiltonian eigenstate, and not built from any particular
    automaton). onsite_expectation/two_point_correlator only ever read
    .sites_uc/.n_uc/.AC/.AR off their `result` argument, so they accept a
    UniformMPS directly, no changes needed there. `eta` is imps_sum's own
    norm diagnostic (mirrors idmrg.PeriodicMPS.eta)."""

    def __init__(self, sites_uc, n_uc, D, d_g, AL, AR, C, AC, eta):
        self.sites_uc = sites_uc
        self.n_uc = n_uc
        self.D = D
        self.d_g = d_g
        self.AL = AL
        self.AR = AR
        self.C = C
        self.AC = AC
        self.eta = eta


def imps_sum(result_a, result_b, cutoff=1e-12, maxdim=None):
    """Direct sum of two converged infinite MPS in VUMPS mixed gauge
    (`result_a`/`result_b`: VUMPSResult and/or UniformMPS, any combination
    -- same duck typing as idmrg.imps_sum), returning a new UniformMPS --
    the VUMPS analogue of idmrg.imps_sum. See this module's own "Summing
    two converged VUMPS iMPS" section docstring above for the construction
    and the physical scope limit (in particular, why summing two
    *ordinary* VUMPSResults, always individually normalized so eta=1 on
    both the left and right transfer eigenvalue, reliably raises
    RuntimeError rather than silently returning one arbitrary branch).

    Requires both states to share the same n_uc and grouped physical
    dimension d_g (mirrors idmrg.imps_sum's own per-sublattice physical-
    dimension check, specialized to VUMPS's single-grouped-supersite
    representation) -- raises ValueError otherwise. Bond dimensions need
    not match."""
    n_uc = result_a.n_uc
    if result_b.n_uc != n_uc:
        raise ValueError(
            "vumps.imps_sum: unit-cell size mismatch (result_a.n_uc={}, "
            "result_b.n_uc={})".format(n_uc, result_b.n_uc))
    if result_a.d_g != result_b.d_g:
        raise ValueError(
            "vumps.imps_sum: grouped physical dimension mismatch "
            "(result_a.d_g={}, result_b.d_g={})".format(
                result_a.d_g, result_b.d_g))
    d_g = result_a.d_g
    Da, Db = result_a.AL.shape[0], result_b.AL.shape[0]
    D = Da + Db
    raw = np.zeros((D, d_g, D), dtype=complex)
    raw[:Da, :, :Da] = result_a.AL
    raw[Da:, :, Da:] = result_b.AL

    left, right, phys = (Index(D, tags="Link"), Index(D, tags="Link"),
                          Index(d_g, tags="Site"))
    B = ITensor((left, phys, right), raw)
    AL_list, eta = idmrg._canonicalize_periodic([B], 1, cutoff, maxdim)
    AL_new = idmrg._to_array_lpr(AL_list[0])
    AR_new, C_new, AC_new = _complete_mixed_gauge(AL_new)
    D_new = AL_new.shape[0]
    return UniformMPS(result_a.sites_uc, n_uc, D_new, d_g,
                       AL_new, AR_new, C_new, AC_new, eta)


# == Applying a (bounded) MPO to a converged VUMPS iMPS ======================
#
# apply_mpo(result, W_bulk, cutoff=1e-12, maxdim=None) is the VUMPS-mixed-
# gauge analogue of idmrg.apply_mpo -- see that function's own "Applying a
# (bounded) MPO to the converged iMPS" section docstring in idmrg.py for the
# full derivation (Orus & Vidal, "Infinite time-evolving block decimation
# algorithm beyond unitary evolution", PRB 78, 155117 (2008): grow the
# converged unit cell by the MPO, then re-canonicalize/truncate the grown
# bond dimension back down via the standard two-sided fixed-point
# procedure) and the SAME scope restriction -- W_bulk must represent a
# genuinely *bounded* (non-extensive) periodic operator, the same tensor
# tiled once per unit-cell site with no unconditional "keep accumulating
# forever" self-loop; idmrg.py's own Hamiltonian automaton
# (`_build_periodic_mpo`) is explicitly out of scope here for the identical
# reason it is for idmrg.apply_mpo (see that module's own comment).
#
# `W_bulk` uses EXACTLY the same convention idmrg.apply_mpo itself takes: a
# list of n_uc ITensors, one per unit-cell sublattice position (0..n_uc-1),
# each rank-4 (Left, in, out, Right) with `in` a d_p x d_p physical Index
# matching sites_uc.si(p+1) -- so the identical W_bulk list built for one
# converged ground state can be fed to either idmrg.apply_mpo or this
# function to apply the same operator to the other backend's own converged
# state (e.g. to cross-check both against each other on the same physical
# operation, see examples/idmrg/vumps_apply_mpo/main.py).
#
# Construction: group W_bulk into a single Dw x Dw x d_g x d_g grouped-
# supersite MPO tensor via `_group_automaton` (already used, unmodified, to
# build VUMPS's own grouped Hamiltonian automaton -- grouping n_uc
# site-local MPO tensors into one is a purely structural contraction,
# independent of whether the result is a Hamiltonian automaton or an
# arbitrary bounded operator), then grow AL by it exactly the way
# idmrg.apply_mpo grows U_list (`idmrg.grow_by_mpo`, called on the trivial
# n_uc=1 "periodic chain" VUMPS already works at -- the same reduction
# `imps_sum` above uses for its own `_canonicalize_periodic` call),
# re-canonicalize/truncate (`idmrg._canonicalize_periodic`), and complete
# the resulting truncated left-canonical tensor to the full mixed gauge
# (`_complete_mixed_gauge`, this module's own "Summing two converged VUMPS
# iMPS" section machinery, reused verbatim) -- growing-then-truncating a
# single left-canonical tensor and re-deriving {AR,C,AC} from it afterwards
# is exactly the same bookkeeping step `imps_sum`'s own construction needs,
# for the same reason: neither operation is itself a VUMPS ground-state
# solve, so there is no outer iteration left to fix up an inconsistent
# gauge the way `vumps_ground_state`'s own iteration does.


def apply_mpo(result, W_bulk, cutoff=1e-12, maxdim=None):
    """Apply a periodic (bounded) MPO to the converged VUMPS iMPS `result`
    (a VUMPSResult or UniformMPS), returning a new UniformMPS representing
    W|psi> up to (cutoff, maxdim) truncation -- the VUMPS-mixed-gauge
    analogue of idmrg.apply_mpo. See this module's own "Applying a (bounded)
    MPO to a converged VUMPS iMPS" section docstring above for the
    construction, the W_bulk convention (identical to idmrg.apply_mpo's
    own), and the scope restriction (bounded/local periodic operators
    only)."""
    n_uc = result.n_uc
    d_g = result.d_g
    D = result.AL.shape[0]

    W_arr = _group_automaton(W_bulk, n_uc)
    Dw_l, d_in, d_out, Dw_r = W_arr.shape
    if d_in != d_g or d_out != d_g:
        raise ValueError(
            "vumps.apply_mpo: W_bulk's own grouped physical dimension "
            "({}, {}) does not match result's own d_g={}".format(
                d_in, d_out, d_g))

    left, right = Index(D, tags="Link"), Index(D, tags="Link")
    phys = Index(d_g, tags="Site")
    AL_tensor = ITensor((left, phys, right), result.AL)

    w_left, w_right = Index(Dw_l, tags="Link"), Index(Dw_r, tags="Link")
    W_tensor = ITensor((w_left, phys, phys.prime(1), w_right), W_arr)

    B = idmrg.grow_by_mpo([W_tensor], [AL_tensor], 1)
    B = [_t_noPrime(b, "Site") for b in B]
    AL_list_new, eta = idmrg._canonicalize_periodic(B, 1, cutoff, maxdim)
    AL_new = idmrg._to_array_lpr(AL_list_new[0])
    AR_new, C_new, AC_new = _complete_mixed_gauge(AL_new)
    D_new = AL_new.shape[0]
    return UniformMPS(result.sites_uc, n_uc, D_new, d_g,
                       AL_new, AR_new, C_new, AC_new, eta)
