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

A third investigation pass carried out exactly that cross-check, and it
splits cleanly in two: `Lh` is now independently confirmed correct; `Rh`
is not yet, and is the sharper suspect going forward.

- `IDMRGResult.env_HL`/`env_HR` are rank-3 (ket, bra, mpo-channel) tensors
  -- see `IDMRGResult`'s own docstring and `idmrg.py`'s `_extend_HL`/
  `_extend_HR`. `env_HL`'s own channel axis is a plain, unreversed
  automaton read (new sites are absorbed by right-multiplying the running
  channel-space product, i.e. appended in the *same* left-to-right order
  the S/F channel labels were defined for -- `_S_IDX`=identity/pass-
  through content, `_F_IDX`=the fully-summed accumulator, exactly as in
  any bulk `W` tensor). Concretely: `env_HL`'s `_F_IDX` channel, minus
  `e_cell * (result.niter_done - 1)` (the *unit-cell* count absorbed by
  the snapshot `env_window_boundary` captures -- note `e_cell` here is
  per unit cell, i.e. exactly `2*result.e0` for n_uc=2, *not* per physical
  site, an easy off-by-factor-of-n_uc mistake), reproduces `Lh` directly:
  exact to ~1e-16 at the exactly-solvable J_weak=0 dimer point (D=1), and
  to ~2e-4 relative at J_weak=0.02 (D=4) at niter_done=3, tightening
  further with more macro-iterations -- ordinary geometric convergence,
  nothing surprising. This is a genuine, independent confirmation that
  `Lh` (and so `_regularized_environments`'s left-side construction,
  `Source_L`/`left_action`) is correct -- not something the original
  internal residual check (also in `_regularized_environments`) could
  show, since that check only catches a solve-time bug, not a wrong
  defining equation.
- `env_HR` is built by *prepending* each newly-absorbed site (it grows
  outward from the cut in the opposite direction from `env_HL`), which
  swaps its own channel-axis reading: the newly-absorbed site's *right*
  leg attaches to the old `env_HR`, and its *left* leg becomes the new
  boundary -- the reverse of `env_HL`'s pattern. Worked out by hand (not
  just pattern-matched) what this does to the S/F labels: reading the
  resulting channel matrix product spatially (starting at the cut and
  reading outward) is still a perfectly standard, non-reversed automaton
  read -- but the *open* boundary axis is now a row (not column) index,
  so its `_S_IDX` entry (not `_F_IDX`) is the one meaning "enter this
  segment completely fresh", i.e. the fully-summed accumulator. Confirmed
  exactly at the D=1, J_weak=0 point, where the whole comparison reduces
  to scalars: raw `env_HR` channel `_S_IDX` equals `e_cell * (niter_done
  - 1)` to machine precision (matching `env_HL`'s own `_F_IDX` role,
  channel-swapped as derived), and channel `_F_IDX` equals 1 (the
  identity/pass-through content, matching `env_HL`'s own `_S_IDX`).
- This channel-swapped extraction does *not* reproduce `Rh` at D>1,
  though: at J_weak=0.02, `env_HR`'s own *pending* (bond-carrying, i.e.
  neither `_S_IDX` nor `_F_IDX`) channels come out with entries of size
  ~0.3-0.5 -- two to three orders of magnitude larger than `env_HL`'s own
  analogous pending-channel residuals (~0.01, shrinking to ~0.001 by
  niter_done=30, the expected dangling-half-a-J_weak-bond residual that
  ordinary geometric convergence in macro-iteration count should leave
  behind). Checked from niter_done=3 up through niter_done=30 (forced via
  etol=0 so the run doesn't stop early): `env_HL`'s pending-channel
  residuals shrink steadily with more macro-iterations as expected, but
  `env_HR`'s stay pinned at that same ~0.3-0.5 scale throughout -- so this
  is not simply "`env_HR` needs more macro-iterations to converge", it is
  a persistent, non-decaying discrepancy. It's concentrated exactly in
  the near-degenerate small-eigenvalue subspace of `r` (the 3 indices
  with eigenvalue ~6.3e-6 in the skewed-spectrum clue above, as opposed
  to the 1 dominant ~0.9998 index) -- plausibly (not confirmed) reflecting
  DMRG's own gauge freedom to rotate within that near-degenerate subspace
  rather than a genuinely missing term, but this was not disentangled
  before time ran out on this pass.
- Net effect on the overall investigation: `Lh` is now the *better-
  supported* of the two regularized environments (independently
  cross-checked against a completely different, already-validated
  construction), so suspicion should concentrate on `Rh`/`right_action`/
  `Source_R` (or on whatever makes diagrams 4/5/6a in `_h_eff_action`
  behave differently on the right side than diagrams 1/3/5/6b do on the
  left) rather than treating both sides as equally likely culprits, as
  the previous investigation pass's "wildly asymmetric in norm" clue
  alone was consistent with either side being the broken one.
- Two concrete next steps for whoever picks this up: (1) work out the
  proper "closing" formula for `env_HR`'s own dangling pending-channel
  content -- mirroring `Source_R`'s own `T_mat_a[T_mat_b[r]]` construction
  but through the row/column-reversed reading derived here -- so the
  cross-check above can actually complete on the `Rh` side instead of
  stalling on the pending channels; (2) cheaper first move: this
  investigation's own D=1 test case (the dimer, J_weak=0) has *no* bond
  content crossing the cut at all, so it never stress-tested the pending-
  channel-closing formula in the first place -- a D=1 case that *does*
  have a small nonzero dangling bond (e.g. a field-polarized XX chain
  detuned slightly below full saturation) would isolate exactly that
  formula in a setting simple enough to derive fully by hand, before
  returning to the genuinely D>1 dimer case where the pending-channel
  content is large and the near-degeneracy in `r` complicates things
  further.

A fourth investigation pass consulted the standard tangent-space excitation
literature directly (Haegeman, Pirvu, Weir, Cirac, Osborne, Verschelde,
Verstraete, "Variational matrix product ansatz for dispersion relations",
arXiv:1103.2286 -- the original single-tensor/"uniform gauge" derivation,
matching this module's own single-A construction; and Vanderstraeten,
Haegeman, Verstraete, "Tangent-space methods for uniform matrix product
states", arXiv:1810.07006, whose Sec. 6.3/Eq. (193)-(197) spell out the
mixed-gauge {A_L,A_R,C} version of the same construction with fully
explicit diagrams) to check `_h_eff_action`'s own diagram list against the
literature's. This found one genuine, provable gap -- now fixed -- but it
turned out to be insufficient to explain the D>1 anomaly on its own.

- Both references' full H_eff(k)/H_eff(p) formulas contain terms with a
  nested pseudo-inverse/geometric-series structure ((1-e^{ip}E)^{-1}
  sandwiched between two other environments -- Haegeman 2012's terms
  involving `(1-e^{+-ik}E)^{-1}`; the review's `L1`/`R1`, Eq. (186)/(187)
  and Eq. (195)/(196)) that this module's old diagrams 6a/6b did *not*
  have -- they were single terms (unit-cell separation exactly +-1), not a
  resummed sum over all separations n=+-1,+-2,....
- Re-deriving this module's own diagram 6a from scratch (in its own
  (ket,bra)/M[i,o] conventions, not by force-translating the mixed-gauge
  formula) confirms this gap is real: for ket-B at cell n>=2 (bond at
  cell (0,1), mat_a at 0 contributing to the output, mat_b at 1 applied to
  a *plain* A since B is further out), positions 2..n-1 carry the plain
  transfer map T_A with no operator at all -- a genuine geometric series
  in n, resummed via `_right_momentum_resolvent`'s (I - e^{ik}T_A)^{-1}
  (the same r-outer-l pseudo-inverse regularization `_regularized_
  environments` already uses at k=0, needed since T_A's own dominant
  eigenvalue 1 sits exactly at k=0). The fix (implemented, see
  `_h_eff_action`'s own "Diagram 6a's |n|>=2 tail" section) was verified
  two ways: (i) the closed-form resolvent solve reproduces an explicit
  truncated sum over n up to 400 to ~1e-17; (ii) it is the identity for
  D=1 (confirmed both numerically and by hand -- the tail's own "seed" is
  provably zero there, see the derivation in `_h_eff_action`), so every
  existing D=1 test (exact to ~1e-8-1e-14) is completely unaffected.
- Diagram 6b's own analogous tail (ket-B at cell n<=-2) was *also*
  re-derived from scratch and found to vanish identically for every D, not
  just D=1 -- its own "seed" contracts the exact (l,L)-pair-with-identity
  combination the left-gauge condition already forces to zero (see
  `_h_eff_action`'s docstring) -- so no analogous correction was needed or
  added there.
- However, this fix, despite being independently verified as internally
  consistent and mathematically real, does *not* resolve the D>1 anomaly:
  a controlled before/after comparison on the *same* converged n_uc=2
  dimerized-Heisenberg ground state (D=4, J_weak=0.02, the same reproducer
  as the third pass) shows the fix changes the lowest three (near-
  degenerate, S=1-triplet-like) eigenvalues by only ~1e-5 across the full
  Brillouin zone, versus the ~1e-2 (J_weak/2) amplitude the first-order
  degenerate-perturbation-theory formula predicts -- i.e. it adds real,
  nonzero, previously-entirely-missing k-dependence, but at 2-3 orders of
  magnitude too small a scale to be the (sole) explanation. Repeated at
  J_weak=0.1/0.3/0.5 (where D grows to 7/16/19 -- these larger-J_weak runs
  were not pushed to tight maxm/etol convergence, so are only a rough
  cross-check): the same qualitative pattern holds -- weak but nonzero,
  monotonically-increasing-from-k=0-to-k=pi k-dependence (the *correct*
  qualitative shape/sign of -cos(k)), but suppressed by roughly two orders
  of magnitude versus the naive J_weak-scale expectation.
- Net effect on the overall investigation: diagram 6a's resummation was a
  real, independent bug (now fixed, and left in -- it is strictly more
  correct than before and provably inert at D=1), but it is not *the*
  D>1 bug, or at least not the only one. Given that even the literature's
  own richer mixed-gauge H_eff(p) (Eq. 197) has roughly twice as many
  distinct diagram *types* as this module's 7 (onsite/1, diagrams 2/3,
  Lh/Rh/4/5, and 6a/6b) -- including combinations like `L1`/`R1` fed
  *into* a further `h`-touching term, which this module's single-A
  reformulation has no obvious analogue for -- the most promising
  remaining path is no longer "patch one more diagram in the existing
  single-tensor (uniform-gauge) derivation", which has now had two
  independently-confirmed real bugs fixed without closing the gap, but a
  from-scratch port to the mixed-gauge {A_L,A_R,C} formalism Sec. 6 of
  arXiv:1810.07006 actually uses (A_L to the left of B, A_R to the right,
  built from this module's own already-computed `A`/`r` via `C =
  cholesky(r)` or an eigendecomposition-based square root, A_R = C^-1 A_L
  C -- see the review's Algorithm 2/Eq. 15) and its Eq. (193)-(197)
  H_eff(p) verbatim -- a materially larger rewrite than any single-diagram
  patch tried so far, but one that reuses a fully worked-out, already-
  published diagram list end to end instead of re-deriving pieces of it by
  hand one investigation pass at a time.

A fifth investigation pass started exactly that mixed-gauge port -- and,
before writing any of this module's own code, cross-checked the review's
Eq. (193)-(197) directly against MPSKit.jl (github.com/QuantumKitHub/
MPSKit.jl, src/algorithms/excitation/{quasiparticleexcitation,exci_
transfer_system}.jl and src/environments/qp_envs.jl), the reference
open-source implementation by the same authors, since transcribing the
review's own tensor-network *diagrams* (image-only, no equation source
available) by eye had already once produced a plausible-looking but
ultimately unverifiable result (see the fourth pass). This paid off
immediately: MPSKit's own code confirmed the mixed-gauge H_eff, expressed
against a finite-state-automaton MPO (exactly this module's own W, S/F/
pending-channel structure -- not a raw dense two-site h), reduces (for a
single-site unit cell, matching this module's own already-grouped
supersite) to just 3 physical terms -- "B in center" (h touching B
directly, closed by the *full* automaton-channel-resolved left/right
environments GL/GR), "B to the left"/"B to the right" (h touching an
*adjacent* site, closed by a separately precomputed, resummed "B-
environment" GBL/GBR) -- plus an energy-density shift. This is a much more
tractable, unambiguous starting point than hand-parsing Eq. (197)'s own
12-term diagram list.

- Translating GL/GR/GBL/GBR into this module's own automaton channels
  confirmed the *structure* of the third pass's strongest suspicion
  exactly: `GL` (built from A_L, closed against A_L's own dominant right
  fixed point `r`) is genuinely different in kind from `GR` (must be built
  from A_R -- the *right*-canonical gauge tensor, C^-1 A_L C for
  C C^dagger = r -- closed against A_R's own dominant right fixed point,
  which is the *identity*, not `r`). This module's existing `Rh` (and the
  right-side closures inside diagrams 1/2/4) instead reuse A_L and `r`
  throughout (the *uniform*-gauge convention, internally consistent on its
  own terms -- see Haegeman 2012 and the "Algorithm summary" below, this
  is not a bug in that formulation) -- but is a *materially different*
  quantity from the mixed-gauge `Rh` a from-scratch port needs.
- This distinction is invisible at D=1 for a clean, checkable reason:
  `_dominant_right_fixed_point` normalizes `r` to trace 1 always (see its
  own docstring), and for D=1 a trace-1 1x1 matrix *is* the identity
  ([1]) -- so "close with r" and "close with the identity" are the exact
  same number whenever D=1, regardless of which formulation's closure
  convention is "correct" there. For D>1 they are provably different (e.g.
  `r`'s own diagonal in the J_weak=0.02 dimer reproducer is
  [0.9998,6e-6,6e-6,6e-6], nothing like the identity's [1,1,1,1] diagonal)
  -- confirmed directly: building the mixed-gauge A_R/C machinery (via an
  eigendecomposition square root of `r`) and verifying, to ~1e-15, that
  (i) A_R is genuinely right-canonical (sum_s A_R^s (A_R^s)^dagger = I),
  (ii) A_L C = C A_R holds, (iii) the *mixed* (A_L-ket, A_R-bra) transfer
  map's own dominant right fixed point is exactly C, and (iv) rebuilding
  `Rh` from A_R (closed against the identity instead of `r`, with the
  matching C^dagger-C-based projector regularization) changes its norm by
  five orders of magnitude (~2e-5 to ~1.7) relative to the existing A_L/r
  based one at J_weak=0.02, D=4 -- which is *exactly* the "Lh and Rh come
  out wildly asymmetric in norm (~1.7 vs ~2e-5)" clue the second
  investigation pass already flagged as suspicious but did not explain.
  This is a genuine, verified, quantitatively-confirmed finding, not
  speculation.
- However, a first full attempt at porting every diagram to this
  mixed-gauge convention (Lh via A_L/r unchanged; Rh via A_R/identity;
  diagrams 2/6a's right-side closures switched from A_L/r to A_R/identity;
  diagrams 3/6b symmetrically re-examined for which positions fall in the
  A_L-region vs the A_R-region relative to a *movable* excitation --
  concretely: for diagram 6b with the ket's B at cell n<=-2, the site
  between B and the output position 0 is *not* pure A_L-A_L as the fourth
  pass's tail derivation assumed, since positions right of B's own cell
  but left of 0 are a genuinely *mixed* (A_R-ket, A_L-bra) transfer segment
  -- so the fourth pass's "diagram 6b's tail is exactly zero" conclusion
  does not carry over unmodified to the mixed-gauge setting and was
  reworked using a second mixed fixed point, the (A_R-ket,A_L-bra)
  transfer's own dominant right fixed point, which the same C-based
  derivation shows equals C^dagger) did *not* reproduce the expected
  dispersion either -- the resulting H_eff(k) is not even Hermitian to
  better than ~1e-5 (versus ~1e-13 for this module's existing, shipped
  uniform-gauge diagrams), a clear, concrete sign that at least one more
  sign/conjugate/gauge-region error remains in the diagram 2/3/6a/6b
  mixed-gauge translation attempted here -- most likely in exactly the
  kind of place the fourth pass's own tail derivation shows is easy to get
  wrong (which of A_L/A_R/B sits at which position, for a diagram whose
  "far" position moves with n). This attempt was intentionally *not*
  merged into this module -- it is a strictly different (not strictly
  better) partial implementation than what is shipped, and merging an
  unverified, non-Hermitian H_eff would be worse than the current,
  honestly-scoped D>1 restriction.
- Net effect: the mixed-gauge Rh/A_R finding above is solid and worth
  reusing directly (it independently explains a specific, previously
  unexplained clue), but a *complete, correct* mixed-gauge H_eff still
  needs the same diagram-by-diagram individual validation (each diagram
  checked against a from-scratch finite-ring tensor-network contraction,
  the way the very first investigation pass validated this module's
  original uniform-gauge diagrams one at a time) rather than a single
  end-to-end dispersion comparison, which only tells you *that* something
  is still wrong, not *which* diagram. Concrete next step for whoever
  picks this up: build the same from-scratch finite-ring cross-check this
  module's own diagrams 1-6b were originally validated against (see this
  module's own git history / the design-plan references in the first
  investigation pass), but with A_L/A_R/C substituted in throughout, and
  check diagrams 2, 3, 6a, 6b one at a time (in that order -- 2 and 6a are
  the most-modified since they now use A_R/identity instead of A_L/r; 3
  and 6b are the most likely to have a subtle A_L-vs-A_R region error
  given diagram 6b's tail already needed correcting once above) rather
  than only checking the final assembled H_eff(k)'s Hermiticity/dispersion
  as this pass did.

A sixth investigation pass picked up exactly where the fifth left off --
`vumps.py`'s own module docstring flags "wiring its {AL,AR,C} output into a
corrected mixed-gauge excitation ansatz" as the natural next step -- and
made concrete progress, but did NOT close the gap either; D>1 is still
correctly rejected by `build_excitation_environment` below. Recorded here in
the same spirit as passes 1-5: real, verified sub-findings, not speculation,
even though the overall dispersion is still wrong.

- Two mixed-transfer fixed-point identities needed for a from-scratch
  mixed-gauge port were derived by hand AND confirmed numerically to
  machine precision on a converged `vumps.VUMPSResult` (D=2 TFIM):
  the (A_L-ket, A_R-bra) transfer map T[X] = sum_s A_L^s X (A_R^s)^dagger
  has C as its dominant (eigenvalue-1) RIGHT fixed point (immediate from
  A_L C = C A_R contracted against A_L's own left-isometry sum_s
  (A_L^s)^dagger A_L^s = I) and conj(C) as its dominant LEFT fixed point
  (same identity, conjugated) -- confirmed to ~1e-14 by direct
  diagonalization of the dense (D^2,D^2) transfer matrix, not merely
  algebra. The mirror pair, for the (A_R-ket, A_L-bra) map, is C^dagger
  (right) and C^T (left) -- derived the same way, via A_R's own
  right-isometry sum_s A_R^s(A_R^s)^dagger = I. These let a
  `_right_momentum_resolvent`-style (I - e^{ip}T)^{-1} regularized
  resolvent be built for the mixed transfer exactly the way the D=1-only
  diagrams above already do for the uniform transfer -- implemented as a
  `_mixed_momentum_resolvent(p, AL, AR, C)` prototype and independently
  verified (not just internally self-consistent) against an explicit
  truncated sum over n up to 150 unit cells, agreeing to ~1e-16 (mirroring
  how diagram 6a's own resummation was validated in pass 4).
- With VUMPS supplying AL, AR, C, GL, GR directly (GL/GR already the
  correctly regularized mixed-gauge background environments, reused
  verbatim from `vumps.py`'s own `_h_ac_action`/`_precompute_bond_
  environments` -- these are independently validated already, since VUMPS
  ground-state energies match known references), the mixed-gauge tangent-
  space norm collapses to the *trivial* Euclidean metric
  tr(X'^dagger X) -- unlike the uniform-gauge norm (`_reduce_metric`,
  weighted by the generally ill-conditioned r), because gauge is split
  exactly at the excitation for both bra and ket in the n=0 ("B in
  center") sector, and the norm's n!=0 cross terms vanish identically by
  the same null-space argument that kills diagram 6b's tail (confirmed
  by hand). This removes the whole `_reduce_metric` conditioning problem
  for a mixed-gauge port -- a genuine structural simplification, not just
  an implementation convenience.
- A prototype H_eff(p) assembled from exactly 3 pieces -- vumps.py's own
  `_h_ac_action(B, GL, GR, bond_envs, h1)` for the n=0 ("B in center")
  sector unchanged, plus two momentum-summed tails (mirroring diagrams
  6a/6b above, built from the same mixed resolvent) -- reproduces the
  D=1 dispersion exactly (both tails provably/numerically vanish at D=1,
  confirmed to ~1e-15 against the same XX-chain reference `excitation_
  energies` already matches), a necessary but not sufficient consistency
  check.
- At D>1 this prototype is measurably WRONG in the same qualitative way
  passes 2-5 already found (anomalously suppressed k-dependence, and
  H_eff(p) not Hermitian to better than ~1e-1 relative, worse than the
  ~1e-13 the D=1 diagrams achieve) -- but two SPECIFIC, previously
  unidentified missing diagrams were found and fixed along the way, each
  independently improving the numerical agreement with the exact D>1 TFIM
  free-fermion dispersion (eps(k) = 2*sqrt(J^2+h^2-2*J*h*cos k), checked
  directly against a D=2 VUMPS ground state at J=1, h=2.5):
  (i) an onsite h1 term at the bra's own excitation position also has a
  nonzero momentum-summed tail connecting it to a ket-B further to the
  RIGHT (n>=1) -- structurally identical to diagram 6a's tail (same seed,
  same resolvent) but with no mat_b/pending-channel step, since h1 has
  none; its own tail to the LEFT (n<=-1) is provably zero by the exact
  same argument that kills diagram 6b's tail (the seed is built with
  bra=A_L at the far position, forcing the same null-space contraction to
  vanish). Adding this alone (before finding (ii)) reduced the D=2 TFIM
  dispersion error at p=pi/4 from ~0.42 to ~0.07 (still wrong, but a
  measurable, real improvement, not noise).
  (ii) the OTHER reach-1 bond touching the bra's excitation position (the
  one straddling cells -1,0, whose mat_a sits at -1 -- background, and
  mat_b at 0) analogously has a nonzero tail connecting to a ket-B to the
  right, built from the SAME `Lvec_a` `vumps.py`'s own `_precompute_bond_
  environments` already computes (mat_a's own background content is
  n-independent, only mat_b's site-0 attachment needs to thread through
  the resolvent). Adding this on top of (i) changed the dispersion shape
  further (differences dropped to the ~0.2-0.5 range at J=1,h=2.5 rather
  than ~1 without either fix) but made the Hermiticity violation WORSE
  (~1e-1 relative versus ~1e-2 with only fix (i)) -- i.e. it is not simply
  "the missing piece", either an overall sign/normalization on it is
  wrong, or a further, still-unidentified counterpart diagram is needed
  before the two partially-cancel correctly. This was not resolved.
- An attempt at an independent, from-scratch brute-force cross-check
  (sweeping the actual automaton W tensor explicitly across a window,
  analogous to how `dynamics.py`'s KPM machinery or `idmrg.py`'s own
  HL/HR growth works) repeatedly tripped on a subtle, genuinely confusing
  point worth recording so a future attempt does not repeat it: a
  left-to-right MPO sweep (`apply_transfer_from_left`, one site at a
  time) through a UNIFORM A_L-A_L region does NOT preserve the trace of
  an arbitrary accumulated environment matrix from one step to the next
  (confirmed both by direct hand derivation and numerically) -- even
  though the trivial S-channel seed (identity) genuinely IS an exact
  fixed point of that same sweep (also confirmed, via A_L's own
  left-isometry). Trace-preservation under repeated identity-channel
  pass-through specifically needs the RIGHT-canonical property (which
  only A_R has), while the S-channel-stays-exactly-identity property
  needs the LEFT-canonical property (which only A_L has) -- these are two
  logically separate facts, easy to conflate, and a naive uniform
  automaton sweep silently gives a plausible-looking but wrong answer
  (confirmed directly: exactly 2x too large a per-site energy density on
  a D=2 TFIM background, NOT an obvious red flag at first glance) by
  implicitly assuming both hold simultaneously for whichever tensor sits
  in the window. Reusing `vumps.py`'s own already-validated GL/GR as
  the sweep's boundary conditions (rather than iterating a raw automaton
  sweep for many sites and hoping it converges) avoids this specific trap
  but was not itself successfully debugged to agreement with the diagram-
  based H_eff before time ran out on this pass -- a working, independent
  brute-force validator (matching the standard this module's D=1 diagrams
  were originally held to) is still the single most valuable next step,
  since Hermiticity/dispersion checks alone identify THAT something is
  still missing but not precisely WHICH diagram.
- Net effect: the mixed-gauge/VUMPS approach remains the right direction
  (D=1 exact, two genuine missing diagrams found and fixed, the metric
  problem structurally resolved, the needed fixed-point identities now
  proven rather than guessed) but is not yet complete -- D>1 stays
  unimplemented here. Concrete next steps for whoever picks this up:
  (1) get the brute-force W-sweep cross-check (using GL/GR as boundary
  conditions, not a raw many-site iteration) working and agreeing with
  the n=0 sector first (a pure consistency check against `_h_ac_action`,
  no new physics), THEN use it to test each of (i)/(ii) above and any
  further candidate diagrams one at a time, the same diagram-by-diagram
  discipline pass 1 established; (2) revisit whether (ii)'s sign/
  normalization is simply wrong (try its negation) before assuming a
  whole new diagram is missing; (3) systematically enumerate every
  remaining "attachment point to position 0" x "tail direction"
  combination against the null-space zero-tail argument (proven so far
  for: any tail seeded by bra=A_L at a negative far position) rather than
  continuing to find them one at a time by trial.

A seventh investigation pass picked up exactly the "brute-force W-sweep
cross-check" next step pass six left unfinished -- but replaced it with
something more useful: a from-scratch, dense-operator (no automaton, no
W tensor, no idmrg.py primitives at all) reproduction of the D=1 XX+field
matrix element <B'|H|Phi_p(B)> position-by-position (positions -3..+3
enumerated explicitly with raw Sx/Sy/Sz matrices, ket_at(n) swapped in by
hand), which pins down the *exact* target value pass six's own brute-force
attempt was trying to reach but never got working. This tool -- call it
the "position-explicit reproducer" -- is what let this pass finally
distinguish a genuinely correct diagram from a self-consistent-but-wrong
one (see below), which the previous six passes' own cross-checks (each
validated only against the *final* dispersion/Hermiticity, or against
another piece of this same module's own machinery) could not do.

- This pass first re-attempted pass five/six's own mixed-gauge diagram
  translation, but organized differently: instead of hand-transcribing
  Eq. (193)-(198) of Vanderstraeten/Haegeman/Verstraete (arXiv:1810.07006)
  from their own rendered *diagrams* (image-only, as pass four already
  found unreliable), this pass downloaded the actual PDF and read the
  diagrams directly at high zoom, AND cross-referenced MPSKit.jl's actual
  source (github.com/QuantumKitHub/MPSKit.jl,
  src/environments/qp_envs.jl's `left_excitation_transfer_system`/
  `right_excitation_transfer_system` and the `GBL`/`GBR` recursions
  feeding them) -- the same reference pass five used, but this time
  fetched to look at the *recursion*, not just the final 3-term formula.
  This resolved a real ambiguity pass four/five never pinned down: the
  paper's own `E_L^L`/`E_R^R`/`E_L^R`/`E_R^L` notation (Eq. 184-196) is
  "E_{bra-tensor}^{ket-tensor}" (subscript = which tensor type sits in
  the bra, superscript = which sits in the ket) -- confirmed by matching
  MPSKit's `TransferMatrix(right_gs.AR, O, left_gs.AL)`-style calls
  against the paper's own diagram connectivity at matching points, not
  guessed. Under this reading, `E_L^R` (used in the paper's own `L_1`/
  `L_B`, the "B already passed, propagate rightward" objects) is exactly
  this module's own (A_R-ket, A_L-bra) mixed transfer, and `E_R^L`
  (`R_1`/`R_B`) is the (A_L-ket, A_R-bra) one -- matching pass six's own
  fixed-point identities (C^dagger/C^T for the former, C/conj(C) for the
  latter) exactly, a genuine independent confirmation of that pass's own
  algebra.
- Attempting a *literal* implementation of the paper's own channel-
  resolved GBL/GBR recursion (mirroring MPSKit's `lBs`/`rBs`, i.e.
  building the FULL automaton-channel-resolved environment -- not just
  its already-regularized F-channel slice `vumps.GL`/`GR` -- and
  resumming it against a momentum phase via a dense pseudo-inverse
  (`numpy.linalg.lstsq`, chosen specifically to avoid hand-building a
  separate r-outer-l-style projector for each of the four mixed
  transfers, a real source of sign/index bugs in every previous pass)
  produced a construction that was internally self-consistent -- its own
  fixed-point equation checked out against an explicit truncated-geometric-
  series unrolling to ~1e-16, TWICE (once for the left-growing GBL, once
  for the right-growing GBR) -- but which the position-explicit reproducer
  showed was simply computing the WRONG quantity for the right-growing
  (GBR) half: internally consistent with its own (subtly wrong) recursion,
  not with the actual physics. Concretely: at D=1 (XX+field, same
  reproducer model as passes two-four), GBL's contribution matched the
  reproducer's own "B strictly left of the differentiated position"
  sector exactly (-0.4777+0.1478j vs the reproducer's -0.4777+0.1478j),
  but GBR's matched contribution came out identically zero against a
  target of -0.4777-0.1478j (the *complex conjugate* of GBL's own value --
  an evocative, unexplained near-symmetry of this specific reproducer
  model, not fully understood, but a strong hint the "correct" GBR
  should have been the same order of magnitude as GBL, not zero). No
  combination of swapping `apply_transfer`/`apply_transfer_from_left` or
  the summation-vs-target index in the recursion (4 combinations tried
  explicitly) fixed this without breaking GBL instead -- and removing the
  vumps-GL/GR content from the channel-resolved environments' own F-slice
  (a live hypothesis at the time, "double-counting the regularized
  background") changed nothing numerically (GL/GR were already ~0 for
  this trivial D=1 product-state test case, so this specific hypothesis
  was never actually exercised, and remains untested for a genuinely
  D=1-but-GL/GR-nonzero case). This construction (and the specific bug in
  it) was NOT resolved and was abandoned in favor of the approach below,
  rather than merged -- recorded here so a future pass does not re-spend
  the same effort re-discovering that the "obvious" GBR mirror of a
  working GBL is not simply obtained by swapping A_L<->A_R and
  apply_transfer<->apply_transfer_from_left throughout.
- Falling back to a direct, minimal port of THIS module's own already-
  validated (D=1-exact) diagram 6a/6b -- rather than the channel-resolved
  GBL/GBR abstraction -- worked immediately: diagram 6a's mixed-gauge
  translation is (per pending channel) `cap_right(apply_op_ket(mat_a,
  A_L), phase*cap_bB + phase**2*tail)` with `cap_bB = apply_transfer(
  op_transfer_matrix(B, A_R, mat_b), I)` (the n=1 term, closed with plain
  identity -- NOT `Rh`/`GR`, a real and easy mistake, see below), `tail =
  apply_transfer(op_transfer_matrix(A_L, A_R, mat_b), resolvent(seed))`,
  `seed = apply_transfer(op_transfer_matrix(B, A_R, None), I)`, and
  `resolvent` the (I - e^{ip} T)^{-1} pseudo-inverse of the (A_L-ket,
  A_R-bra) mixed transfer (regularized via plain `lstsq`, no explicit
  projector needed). Diagram 6b mirrors this with A_L/A_R swapped,
  `apply_transfer_from_left`, and phase e^{-ip}, closed via an explicit
  trailing `cap_right(..., I)` (matching this module's own existing D=1
  diagrams' documented "diagrams 3/5/6b need an explicit right-side
  closure" rule). This reproduces the position-explicit reproducer's own
  n=+-1-only values exactly (both components, not just one), and the
  resulting D=1 XX+field dispersion is exact to machine precision
  (~1e-15) AND exactly Hermitian (~1e-16) at every momentum tested --
  matching this module's own existing D=1 diagrams' accuracy, from a
  completely independent (VUMPS-based, mixed-gauge) code path. This is a
  genuine, useful confirmation that VUMPS's {A_L,A_R,C} plus this simpler
  6a/6b-style construction (not the channel-resolved one above) is
  internally sound at D=1 -- but D=1 alone cannot distinguish this from
  the FIVE previous passes' own D=1-exact-but-D>1-wrong attempts, since
  D=1 makes every tail/resolvent trivially vanish regardless of whether
  the underlying construction generalizes correctly.
- At D=2 (TFIM, J=1, h=2.5, the same reproducer passes two/six used),
  this minimal 6a/6b port alone is NOT Hermitian (~0.1-0.2 relative) and
  the dispersion is wrong by O(1) -- the same failure signature every
  prior pass hit. Adding pass six's own two extra diagrams (i) (an onsite
  h1 tail: `cap_right(apply_op_ket(h1, A_L), phase*resolvent(seed))`,
  i.e. exactly diagram 6a's own tail machinery with the mat_a/mat_b
  pending-channel step dropped since h1 has none) and (ii) (the OTHER
  reach-1 bond touching position 0: `cap_right(cap_left(Lvec_a,
  apply_op_ket(mat_b, A_L)), phase*resolvent(seed))`, reusing
  `vumps._precompute_bond_environments`'s own `Lvec_a` verbatim)
  reproduces pass six's own qualitative finding exactly, now from an
  independently-written implementation: (i) alone changes the dispersion
  (real, nonzero effect) but WORSENS Hermiticity (~0.1 -> ~1-2); adding
  (ii) on top narrows the dispersion error further (diff shrinks from
  O(1) to ~0.1-0.4 across the Brillouin zone) but worsens Hermiticity
  again (~1-2 -> ~2-3.5). Negating (i) alone, (ii) alone (holding the
  other fixed), and both together (four sign combinations in total, a
  strict superset of pass six's own "try (ii)'s negation" suggestion)
  were all tried explicitly here -- none improves on the untouched
  (i)+(ii), both positive, combination, and none comes close to
  Hermitian. This is an independent confirmation (different code, same
  physical conclusion) that (i) and (ii) are real, necessary corrections
  -- not artifacts of pass six's own particular implementation -- but
  that at least one more diagram (or a genuine error in how (i)/(ii)
  interact with each other or with 6a/6b's own tail) remains missing,
  ruling out "it's just a sign" as the explanation pass six's own next-
  steps list had flagged as the cheapest thing to try first.
- Net effect: no closer to a working D>1 ansatz than pass six, but two
  concrete, reusable assets for whoever picks this up next: (1) the
  position-explicit reproducer itself (enumerate ket_at(n) by hand with
  raw Sx/Sy/Sz matrices, no automaton) is a fast, unambiguous way to get
  the *exact* target value for any individual diagram's own n=+-1
  contribution at D=1, and should be built out to D=2 (TFIM) next --
  computationally harder (D=2 breaks the "everything is a scalar" shortcut
  D=1 enjoys, needing a genuine small-D finite-window MPO contraction) but
  the only way to pin down which of 6a/6b/(i)/(ii)/"something else" is
  wrong AT D>1 specifically, rather than only checking self-consistency
  (as this pass's abandoned GBL/GBR attempt shows can pass even when the
  underlying construction is wrong) or the final dispersion/Hermiticity
  (which identifies THAT something is wrong but, as six diagram-counting
  passes now confirm, not efficiently WHICH diagram); (2) the E_L^R/E_R^L
  notation key above (bra-subscript, ket-superscript), independently
  re-derived and cross-checked against MPSKit's own source rather than
  guessed from the paper's images alone, removes one previously-open
  translation ambiguity for any future attempt at the paper's own full
  Eq. (195)-(197) (this pass's own channel-resolved attempt shows that
  formula's structural shape is NOT simply "correct diagrams, wrong
  sign" -- something more fundamental about the GBR/rightward half's own
  recursion was wrong, and remains unidentified).

An eighth investigation pass began by consulting an independent model
(Fable) for a second opinion on strategy before touching any code again,
given seven passes' worth of "internally self-consistent but wrong"
results. Its recommendation, borne out below: stop patching diagrams and
first check whether the bug is in shared/background infrastructure that
is numerically DEGENERATE (and therefore invisible) at D=1 -- every
validation this module has ever used (D=1 exactness, D=1 finite-ring
contractions, pass seven's own D=1 position-explicit reproducer) is
structurally blind to exactly that class of bug, since T, the mixed-
transfer projector, and C all trivialize to scalars/identity at D=1.

- A cheap, decisive cross-check the Fable consultation prompted (not run
  before, despite being implied by data already in hand): pass six/seven's
  own mixed-gauge attempts already run on VUMPS's GL/GR -- which are
  independently validated via ground-state ENERGY agreement with known
  references (see vumps.py's own test suite), a completely different code
  path from idmrg.py's growing-algorithm env_HL/env_HR that the THIRD
  pass's own investigation flagged as having an unexplained, non-decaying
  env_HR residual. Since passes six/seven's mixed-gauge H_eff(p) is built
  entirely from VUMPS's own (AL,AR,C,GL,GR), NOT from idmrg.py's env_HL/
  env_HR at all, and still fails the exact same way (suppressed
  dispersion, non-Hermitian H_eff), this is strong evidence the third
  pass's env_HR anomaly -- while real, and still unexplained -- is NOT the
  root cause of the excitation bug: the bug survives switching to a
  completely independent, already-validated ground-state environment
  construction. This deprioritizes (does not close) that investigation
  thread and refocuses suspicion on the excitation-diagram/mixed-transfer
  machinery specifically, which both approaches DO share.
- This pass then rebuilt pass seven's own mixed-gauge H_eff(p)
  construction from scratch, directly from that pass's own prose
  description (6a/6b mixed-gauge port plus the sixth pass's (i)/(ii)
  extra diagrams), reusing vumps.py's already-validated `_h_ac_action`/
  GL/GR/bond_envs verbatim for the momentum-independent (n=0, "B in
  center") sector rather than re-deriving it. Re-validated against the
  SHIPPED, already-D=1-exact uniform-gauge `excitation_energies` (a
  different, independent D=1 ground state -- idmrg.py's growing
  algorithm, not VUMPS) on the same XX+field Hamiltonian: agrees to
  ~1e-16 and is exactly Hermitian (~1e-17-1e-25) at every momentum tested
  -- confirming this reconstruction faithfully reproduces what pass seven
  itself reported, before trusting anything it says about D>1.
- Two new candidate diagrams were then tried and DEFINITIVELY RULED OUT,
  closing off a real, previously-untested hypothesis rather than merely
  re-trying sign flips of existing pieces (which is all pass seven itself
  had time for): a leftward mirror of the sixth pass's onsite-h1 tail
  (already known to vanish, per pass six's own null-space argument -- kept
  here purely as an implementation sanity check) AND, new this pass, a
  leftward mirror of the sixth pass's own "(ii)" diagram (the OTHER
  reach-1-bond-touching-position-0 correction) -- structurally the
  natural leftward analogue nothing in passes six/seven records having
  actually tried. Checked in ISOLATION (not just "adding it changes
  nothing in the final dispersion", which could also mean a
  bug silently no-ops it) at D=2 TFIM (J=1, h=2.5, nonzero h1): both
  mirrors evaluate to a norm of ~1e-16-1e-17 on a random tangent vector,
  i.e. exactly zero, not merely small. This extends pass six's own
  null-space argument (proven so far only for "any tail seeded by
  bra=A_L at a negative far position") to cover this second diagram
  family too, and rules out "a missing mirror diagram" as the explanation
  for the remaining gap -- a genuine negative result, not merely another
  attempt that didn't help.
- A new, more granular data point (not reported by passes six/seven,
  which only ever checked the FULL i+ii-included construction's
  Hermiticity): at D=2 TFIM, diagrams 6a+6b ALONE (i.e. the mixed-gauge
  port of this module's own original, D=1-exact diagram list, with
  neither of the sixth pass's extra corrections) reproduce the classic,
  ORIGINAL failure signature first seen in the very first investigation
  pass almost exactly -- an anomalously FLAT dispersion (E(p) ranging only
  4.86-4.97 across the full Brillouin zone at J=1,h=2.5, versus the exact
  free-fermion band's own 3.00-6.99 range) -- with a comparatively SMALL
  Hermiticity defect (~2e-3 to 2e-2 relative, depending on momentum).
  Adding the sixth pass's (i)+(ii) corrections on top genuinely fixes the
  flatness -- E(p) now spans 2.5-6.5, the correct order of magnitude and
  qualitative shape (low near p=0, peaking near p~1.5-2, matching the
  exact band's own shape) -- but at the cost of a much LARGER Hermiticity
  defect (~0.15-0.31 relative), confirming passes six/seven's own
  qualitative finding (dispersion improves, Hermiticity worsens) with a
  cleaner, isolated before/after comparison than either pass's own
  end-to-end-only checks. This is a genuine, load-bearing new
  observation: it means (i)/(ii) supply real, necessary k-dependence (not
  spurious content an alternative sign would remove -- consistent with
  pass seven's own four-sign-combination sweep finding no improvement
  over the "both positive" combination), but are evidently missing a
  Hermitian-restoring partner of some kind that is NOT simply their own
  mirror image (now ruled out above) and NOT simply a sign flip (already
  ruled out by pass seven).
- Given this, the most likely remaining gap is structural rather than a
  further sign/mirror search: pass four's own literature comparison
  already flagged, and never followed up on, that the full published
  mixed-gauge H_eff(p) (Vanderstraeten/Haegeman/Verstraete arXiv:
  1810.07006 Eq. 197) has roughly twice as many distinct diagram TYPES as
  this module's simpler list, including terms where a resolvent-summed
  "L1"/"R1"-style object feeds INTO a further, separate Hamiltonian-
  touching term (a genuinely doubly-nested construction), which this
  module's diagrams 1-6b/(i)/(ii) have no analogue of at all -- every
  diagram implemented so far (in any of the eight passes) has exactly one
  H-term insertion. A term of that doubly-nested shape is a structurally
  different, larger addition than anything tried in passes five through
  eight (which only ever varied WHICH existing single-insertion diagrams
  were included, and with what sign), and is the most concrete lead this
  investigation currently has.
- Net effect: no closer to a working D>1 ansatz than pass seven (D>1
  remains correctly rejected by `build_excitation_environment` below),
  but two real, reusable results for whoever picks this up next: (1) the
  third pass's env_HR anomaly is now known NOT to be the root cause (via
  the VUMPS-bypass argument above), so it no longer needs chasing before
  further excitation-diagram work; (2) both natural "mirror" completions
  of the sixth pass's (i)/(ii) diagrams are now definitively, numerically
  ruled out (not merely untried), narrowing the search specifically
  toward a genuinely new, doubly-nested diagram type mirroring the
  published Eq. 197's own L1/R1-fed-into-h structure, rather than further
  recombination of the diagrams already tried across all eight passes so
  far.

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
   finite list of diagram *types* (see `_h_eff_action`'s own docstring),
   most of which genuinely only connect unit-cell separations 0 and +-1
   directly (the gauge condition V_L^dagger @ A = 0 kills anything
   further, given this module's Hamiltonians only have reach-1 bonds) --
   except diagram 6a, whose own tail (separations n>=2) is a genuine
   geometric series, resummed via `_right_momentum_resolvent`'s
   (I - e^{ik}T_A)^{-1} pseudo-inverse (diagram 6b's analogous tail
   vanishes identically instead, see `_h_eff_action`'s own docstring for
   why the two sides differ).
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


def _right_momentum_resolvent(k, env):
    """Solve operator for (I - e^{ik} T_A)[Y] = source, T_A[X] =
    idmrg._apply_transfer(E_id, X) -- the geometric-series resummation
    diagram 6a's "tail" (unit-cell separations |n|>=2, see `_h_eff_action`'s
    own docstring for the derivation of why this tail is generically
    nonzero for D>1) requires. Built once per k (not once per B/X, since
    the operator itself is k-dependent but B-independent) and reused across
    every basis vector `_build_H_eff_dense` feeds through `_h_eff_action`.

    At k=0 exactly, T_A's dominant eigenvalue 1 (along r) makes (I-T_A)
    singular, so this adds the same r-outer-l projector regularization
    `_regularized_environments`'s own `right_action` uses -- valid here
    too since the actual source fed in (see `_h_eff_action`'s diagram 6a)
    is guaranteed traceless against l by the very same left-gauge identity
    that makes diagram 6b's analogous tail vanish outright (see below) --
    confirmed directly: trace(l @ apply_transfer(op_transfer_matrix(B,A,
    None), r)) collapses to trace(r-weighted) of a quantity that is
    identically zero for *every* (row,col) entry by construction of
    `_null_space_left` (V_L^dagger @ A_mat = 0 kills the plain trace of
    _op_transfer_matrix(B,A,None) over its own (l,L) pair, for any closure
    on its (r,R) pair -- including the r-weighted one used here)."""
    D = env.D
    A, r, l = env.A, env.r, env.l
    E_id = _op_transfer_matrix(A, A, None)
    phase = np.exp(1j * k)
    at_zero = abs(phase - 1.0) < 1e-10

    def action(X):
        out = X - phase * idmrg._apply_transfer(E_id, X)
        if at_zero:
            out = out + r * np.trace(l @ X)
        return out

    Mat = _dense_linear_map(D, action)

    def solve(source):
        return np.linalg.solve(Mat, source.reshape(-1)).reshape(D, D)

    return solve


def _h_eff_action(k, X, env, resolvent):
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
       everything is the plain converged background again. B's own
       position never moves in this diagram, so there is no momentum sum
       to resum here -- k does not enter.
    4/5. Background Hamiltonian content strictly to the right/left of B,
       with B itself untouched -- exactly Rh/Lh (built once, k-independent).
    6a/6b. The momentum-summed piece: a reach-1 bond connecting the bra's
       excitation (fixed at cell 0 by the projection V_L^dagger applied at
       the very end) to the SAME B, weighted e^{+-ikn}, sitting on the KET
       side at cell n=+-1,+-2,... instead, with an ordinary A at every cell
       strictly between 0 and n. This is the only place k enters.

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

    == Diagram 6a's |n|>=2 tail (the fix for the D>1 dispersion bug) ==

    For |n|>=2, the site(s) strictly between B's own two positions (0 in
    the bra, n in the ket) carry an *unperturbed* background on both bra
    and ket -- i.e. plain repeated application of the ordinary transfer
    map T_A (no operator, no B) -- so summing over all n>=2 is a genuine
    geometric series, resummed here via `_right_momentum_resolvent`'s
    (I - e^{ik} T_A)^{-1} rather than truncated at the single n=1 term a
    prior version of this function stopped at (see git history). Crucially
    this tail is *not* generically zero for D>1: the "seed" it resums,
    apply_transfer(_op_transfer_matrix(B,A,None), r), leaves the (l,L)
    pair of _op_transfer_matrix(B,A,None) *open* (closed only against the
    metric r on its (r,R) pair), and the left-gauge condition
    (`_null_space_left`'s V_L^dagger @ A_mat = 0) only forces the *trace*
    (l=L, closed with plain identity) of that same tensor to vanish, not
    its off-diagonal (l!=L) entries -- so for D=1 (where (l,L) can only
    ever take their single, necessarily-diagonal value) the tail vanishes
    identically and this reduces exactly to the old, single-term diagram,
    but for D>1 it generically does not, which is exactly the observed
    failure signature (this module's own "KNOWN LIMITATION" section):
    exact to ~14 digits at D=1, anomalously flat/wrong at D>1.

    Diagram 6b's own analogous tail (ket-B at cell n<=-2) is, by contrast,
    provably zero for *every* D, not merely small -- no resummation is
    needed there. Its "seed" is apply_transfer_from_left(
    _op_transfer_matrix(B,A,None), l), which contracts l against exactly
    the (l,L) pair the left-gauge condition already forces to zero (l=I
    here, so this contraction *is*, verbatim, the plain (l=L) trace the
    gauge condition kills) -- so diagram 6b's single existing term is
    already exact, and no `resolvent`-based correction is added for it."""
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

        cap_bB_r = idmrg._apply_transfer(_op_transfer_matrix(B, A, mat_b), r)  # n=1
        seed = idmrg._apply_transfer(_op_transfer_matrix(B, A, None), r)      # n>=2 seed
        tail = idmrg._apply_transfer(_op_transfer_matrix(A, A, mat_b), resolvent(seed))
        Y = Y + _cap_right(_apply_op_ket(mat_a, A),
                            phase_p * cap_bB_r + phase_p**2 * tail)         # diagram 6a

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
    itself rather than just the norm.

    `_right_momentum_resolvent(k, env)` is built once here (it only
    depends on k, not on the trial vector X) and reused across all n
    basis-vector evaluations of `_h_eff_action`, rather than rebuilt (and
    re-factorized as a dense (D^2,D^2) linear solve) on every single call
    -- an O(D^6)-per-call cost that would otherwise dominate this
    function's own O(D^6) total cost by a factor of n."""
    D, d_g = env.D, env.d_g
    Dx = D * (d_g - 1)
    V, isqrt = env.r_V, env.r_isqrt
    kdim = V.shape[1]
    n = Dx * kdim
    resolvent = _right_momentum_resolvent(k, env)
    H = np.zeros((n, n), dtype=complex)
    Xt = np.zeros((Dx, kdim), dtype=complex)
    for j in range(n):
        Xt.flat[j] = 1.0
        X = Xt @ (isqrt[:, None] * V.conj().T)
        Y = _h_eff_action(k, X, env, resolvent)
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
