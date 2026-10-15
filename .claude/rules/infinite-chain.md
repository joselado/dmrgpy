---
paths:
  - "src/dmrgpy/infinitechain.py"
  - "src/dmrgpy/pyitensor/idmrg*.py"
  - "src/dmrgpy/pyitensor/vumps*.py"
  - "src/dmrgpy/mpscpp3/chain_session.h"
  - "tests/test_infinite*.py"
  - "tests/test_idmrg*.py"
  - "tests/test_vumps*.py"
  - "tests/test_lanczos_residual_criterion.py"
  - "examples/idmrg/**"
---

# Infinite chains: iDMRG and VUMPS

These are the notes CLAUDE.md used to carry under its "ITensor library"
header; they load when `infinitechain.py`, `pyitensor/idmrg*.py`,
`pyitensor/vumps*.py`, `mpscpp3/chain_session.h` or one of their tests is
read or edited. The short version of the hard rules is in CLAUDE.md's
"Infinite chains" section; this file is the reasoning and the measurements
behind them.

## Upstream reference

ITensor's own infinite-MPS package,
<https://github.com/ITensor/ITensorInfiniteMPS.jl>, implements iDMRG, VUMPS,
the infinite MPO/`InfiniteMPS` types and the transfer-matrix fixed-point
machinery that `infinitechain.py`, `pyitensor/idmrg.py`,
`pyitensor/vumps.py` and the `ic_*`/`vx_*` helpers in
`mpscpp3/chain_session.h` reimplement here. It is not vendored in this repo
and dmrgpy does not depend on it; consult it as the reference implementation
when an infinite-chain algorithm (unit-cell gauging, McCulloch wavefunction
prediction, environment energy subtraction, fixed-point conventions for
`<Sz>`/correlators) looks wrong or ambiguous, the same way `mpscppN/ITensor`
is consulted for the finite algorithms.

## Subspace expansion on the D-ramp

One thing has already been ported from upstream: `pyitensor/vumps.py`'s
`_subspace_expand` (after upstream's `subspace_expansion.jl`) is how the
`"python"` VUMPS D-ramp warm-starts each rung. It picks the new bond
directions from the SVD of `H` projected into `AL`/`AR`'s null spaces,
instead of the noise-padding `_grow_initial_state` still does on
`itensor_version=3`. The expansion cannot change the state or the energy
(both enlarged tensors stay exactly isometric, `C` is embedded with zeros in
the new block), only how fast the next rung converges: a wash on most
models, and 4.6x faster plus 5x more accurate (and run-to-run reproducible,
which the noise start is not) on the hardest one tested, a gapless `D=8`
Heisenberg chain. `pyitensor/vumps_ms.py` (the sequential solver) has the
per-bond version too, a loop over every bond of the cell, which is the shape
upstream's own `subspace_expansion(psi, H)` has. Two things generalize
rather than repeat: bond `n`'s two null spaces come from different tensors
(`AL[n]` and `AR[n+1]`), so `AR[n]`'s new directions come from bond `n-1`;
and the cell's last bond straddles the cell boundary, where the right
environment, automaton and tensor are site 0's of the next cell, by
periodicity the same `GR[0]`/`W_list[0]`/`AR[0]` this cell already holds.
Every bond grows by the same amount so the state stays uniform-`D`. Measured
the same way (`nrestarts=1`, median of 5): a wash on TFIM and on the easy
rows, 4.3x faster on gapless Heisenberg `n_uc=2 D=8` (1701 to 909
iterations, 21.4s to 5.0s), and on `n_uc=1 D=8` (where neither start
converges within `maxiter`) 1.4x faster with the accuracy spread about 5x
tighter at both ends. `itensor_version=3`'s own `vms_ground_state` has not
picked up subspace expansion; it noise-pads (`Chain::vms_grow_init`, the
per-site `vumps_grow_init`). It did for a while have neither: its ramp's
`reuse` test compared the previous rung's tensor size against the new `D`,
so it could never hold across a rung and every rung started from pure
noise. That is not a speed nicety on a model whose exact state is smaller
than the requested `D`, see "Redundant bond dimension" below. The residual
criterion, the AC/C phase alignment and the variational safety net are now
ported too.

## Couplings reaching further than one unit cell

The sequential solver is also the route for couplings reaching further than
one unit cell, which used to be rejected outright
(`infinitechain._canonicalize_hamiltonian` refused any term touching both
the previous and the next cell). A Hamiltonian with reach>1 has no grouped
reach-1 representation on that cell at all, so `vumps_ground_state` sends it
to `vumps_ms` regardless of `n_uc`; `idmrg_excitations.is_reach_one` read
off the grouped automaton is the test, `_check_reach_one` is now just the
raising wrapper around it, kept for the excitation ansatz. That is what
makes a long-range model affordable: range-R costs the sequential path R
extra automaton channels (linear), where the only previous route, rewriting
the chain on an `n_uc >= R` cell so every coupling is reach-1, costs the
grouped path a `d**R` supersite. `get_operator(name, i, group=c)` takes an
integer cell offset (`c = -1, 0, 1` being `L`, `C`, `R`) to write one, and
`_canonicalize_hamiltonian` canonicalizes a term's position by translating
it whole onto the cell its leftmost site lives in (well defined at any
reach, identical to the old three-cell rules everywhere those applied). Two
things worth knowing before touching this: `gs_method="idmrg"` needed no
change at all, since its growth loop consumes whatever `_build_periodic_mpo`
builds, and `_active_channels_at` has always carried one pending channel per
site of a term's reach, confirmed directly against exact answers, so do NOT
"fix" that by adding a guard; and `_window_hamiltonian` (`kpm_finite`'s
open-boundary tiling) tests the fit per term rather than dropping one fixed
copy, because a longer-ranged term runs off the end of the window sooner.
`itensor_version=3` has the same dispatch (`Chain::vumps_ground_state`'s
`use_multisite`, on `vumps_is_reach_one`) and agrees with `"python"` to
2e-13. `vev`/`correlator` follow along on it too, at any `n_uc` and any
reach. They briefly did not, and the reason is worth keeping: the C++ picks
the SOLVER on `n_uc>2 || reach>1` while `infinitechain.py` picked the READER
on `n_uc>2` alone, so a reach>1 chain on a short cell ran the sequential
solver and was then handed to the grouped reader, which had no snapshot,
even though `Chain::vms_onsite_expectation`/`vms_two_point_correlator`
could answer it and had existed since `b81b6a9`.
`Chain::vumps_onsite_expectation`/`vumps_two_point_correlator` now fall
through to the `vms_*` pair on `have_vms_snapshot_`, i.e. the snapshot that
exists decides the reader, and the Python side has no dispatch of its own
left to drift. `tests/test_infinite_long_range.py`,
`examples/idmrg/long_range_infinite_chain`.

## Redundant bond dimension

Redundant bond dimension used to be fatal on `itensor_version=3`. A gapped
model's exact state often needs fewer directions than the requested `maxm`
(a field-polarized chain needs exactly one), the extra ones then carry no
weight, and the transfer matrix picks up a decoupled unimodular block: a
genuinely DEGENERATE dominant eigenvalue, which looks identical to the one
thing `vx_check_perron_nondegenerate` exists to reject (a "cat state": two
branches with matched nonzero weight). So `gs_energy()` raised "every
attempt at D=... failed" for any such model at `maxm>1` on the sequential
solver, and the grouped one survived only by never landing exactly on the
degeneracy (its second eigenvalue was measured at 0.99996, just outside the
guard's 1e-9). Both environment builders now fall back to the fixed points
the state itself names, `C C^dag`/`C^dag C` (`Chain::vx_bond_fixed_points`),
exact under redundancy, and for a real cat state the same branch mixture an
unchecked eigensolver would give, where the eigensolver may instead return
an arbitrary single branch with its own wrong energy. Note what was tried
first and is the WRONG shape: a threshold on `C`'s weight spectrum to tell
redundancy from a cat state. The guard trips mid-convergence, where the
redundant direction is still on its way down, so the ratio to catch is a
moving number (2.7e-9, 1.2e-4 and 1.6e-2 on three cells of the same model)
with no defensible cutoff.

The pure-Python reference has since gone one step further, and the ORDERING
is the point: both its environment builders (`vumps_ms._cell_fixed_points`
on the sequential path, `vumps._transfer_fixed_points` on the grouped one)
prefer the bond candidate, `C C^dag` for the AL transfer's right fixed
point and `conj(C^dag C)` for the AR transfer's left one in this codebase's
`X[ket,bra]` ordering, whenever it reproduces itself under the transfer map
to a residual <= 1e-6, which in mixed canonical gauge is an exact algebraic
identity and therefore a yes/no test rather than a tuned threshold, so at
convergence no eigensolve runs at all; the guarded eigensolver runs only
when the residual says the gauge relation does not hold, and the two are
cross-checked against each other by residual. That is bond-candidate-FIRST,
and `itensor_version=3` got there on 2026-09-12: `Chain::vx_choose_fixed_point`,
one shared selection that both `vumps_build_environments` and
`vms_environments` go through, where each used to hold its own
`catch (ITError const&)` reaching `vx_bond_fixed_points` as a failure
fallback.
The difference between the two orderings was measurable and was not in the
C++'s favour: on a field-polarized reach-1 one-site cell at `maxm=4`,
`gs_energy()` on `itensor_version="python"` raised "every attempt at D=4
failed" in 7 of 20 runs before its own change and 0 of 20 after, while
`itensor_version=3` returned energies below the exact variational minimum
at a rate that moved between runs and builds (13 and 17 of 80 measured at
D=6 on the grouped cell, worst 2.3e-8; worst 2.8e-3 on the sequential one)
where `"python"` did so in 0 of 40. It is now 0 of 330 past 1e-12 across six
cells. NUMBERS CHANGE: any `itensor_version=3` VUMPS energy from before
2026-09-12 on a model whose exact state is smaller than the requested
`maxm` is not comparable, on either solver; those were wrong, not merely
less converged. Energies from runs that were already fine move only at
about 1e-15. Quote the threshold with any rate, as the known-issue file
itself instructs; see `docs/known_issue_v3_vumps_variational_floor.md`, now
FIXED, and `tests/test_vumps_redundant_bond_dimension.py`, whose narrow
`xfail` came out and which gained a repeated-run guard (a 6-in-30 defect is
invisible about 80% of the time to a one-run test).

Porting that ordering is what caught a latent defect in the primitive
itself, and the lesson generalizes: `vx_bond_fixed_points`' LEFT candidate
was `C^dag C` where the `X[ket,bra]` ordering needs `conj(C^dag C)`, the
transpose. Nothing had caught it because nothing had ever measured that
candidate; adding the residual read it immediately (4e-15 for the correct
orientation against 0.38 to 0.53 for the transpose, at convergence on the
polarized cell at D=6). It survived because the only two models this
function had been exercised on, a field-polarized chain, whose converged
`C` is real diagonal, and AKLT, whose `C` is real, both have a real
symmetric `C^dag C`, where the two orientations coincide. A fallback that is
only reached on failure is a fallback nothing validates.

## The Lanczos stopping criterion

Reading that led to a bigger, unrelated find, and it is the one to know
about: `pyitensor/dmrg.py`'s `_lanczos_ground_state` stops when the lowest
Ritz value stops moving. For a Hermitian operator the value error is
quadratic in the eigenvector error, so a value tolerance of 1e-12 caps the
eigenvector at about 1e-6: fine for finite DMRG, which reports energies, and
fatal for VUMPS, whose convergence criterion is a norm difference between
`AC` and `C`, two independently-solved eigenvectors. VUMPS's gauge mismatch
therefore floored at about 1e-6 and `tol=1e-10` was unreachable at any
number of outer iterations, while the energy was already correct, which is
exactly why nothing in `tests/` caught it. `_lanczos_ground_state` has an
opt-in `residual_tol` (Ritz residual `beta*|s_k|`, free from the
tridiagonal problem) that `vumps.py`/`vumps_ms.py` pass as `tol/10`; every
other caller is byte-identical. Alongside it, `AC`/`C` are sign-aligned
against the vectors they were warm-started from: solved independently,
their relative sign was arbitrary, and a flip made the mismatch read about
4 instead of about 1e-6 in 11 of every 100 iterations. Measured through the
public driver: `D=8` TFIM went from `converged` in 0/3 runs to 3/3, and
32.6s to 6.7s (g=1.5), 38.4s to 4.6s (critical g=1.0), energies unchanged
to every printed digit.

One consequence looked alarming and, measured, is not: finite DMRG's own
local solves carry the same ~1e-6 per-solve eigenvector cap, so anything
consuming a DMRG wavefunction rather than its energy was flagged as worth a
look. It has had one (2026-09-12) and the answer is that the cap does not
reach any consumer. On an 8-site Heisenberg chain at FULL bond dimension,
where the MPS can represent the exact state so the local solver is the only
error left, energies, `<Sz_i>`, two-site correlators, the bond entropy,
`<wf|wf>` and the first three excited states all agree with ED to
1e-13..1e-15, and forcing `residual_tol=1e-13` into finite DMRG's own two
`_lanczos_ground_state` call sites (`dmrg.py:421`/`:726`) moves every one of
them by <= 5e-15. Below full bond dimension the error is
truncation-dominated and the criterion makes no systematic difference either
way: swept over `maxm` in 4..32 against ED, the loose/tight difference in the
excited-state errors is non-monotonic noise (tight is slightly worse at
maxm=8 and 12), not an improvement. The mechanism is the sweep: each local
solve is warm-started from the previous sweep's own tensor and the sweep is
iterated to convergence, so a per-solve vector error contracts instead of
accumulating, which is exactly what VUMPS lacks, its criterion being a
difference between two independently solved eigenvectors with no outer
contraction. So do not "fix" finite DMRG's criterion on this reasoning; it
costs Krylov iterations and buys nothing measurable.

The C++ `itensor_version=3` port carried the eigenvalue criterion for a
while after this, so `Chain::vumps_ground_state` floored at about 1e-6 and
reported `converged=False` at `D>=8`, a documented divergence between the
two VUMPS backends. `Chain::vx_lanczos_ground_state` now has the same
`residual_tol` (both VUMPS drivers pass `tol/10`; every other caller is
byte-identical; since 2026-09-25b finding 24 its reference
`max(1,|lambda|)` is capped by the Hamiltonian's own unit, its largest
non-identity coefficient, on both backends, so that it holds in small units
and next to an offset) and `vx_align_phase` the same alignment, removing the
whole phase rather than just the sign, since a vector from the dense
`vx_dense_eig_max_` branch comes out of `zheev` with no fixed phase at all.
Re-measured through the public driver at `D=8`, `nrestarts=6`: `converged`
in 3/3 runs at both `g=1.5` and critical `g=1.0`, from 0/3.
`tests/test_lanczos_residual_criterion.py` covers both backends. See also
`docs/documentation.md`'s D-ramp section,
`tests/test_vumps_subspace_expansion.py`, and the measurement tables in
`vumps.py`'s own section comments.

That fix also overturned a previously-recorded diagnosis, which is worth
knowing before trusting a similar one: the tangent-space excitation ansatz's
accuracy loss at redundant bond dimension (excitation error 1.8e-15 at
`maxm=1` up to 3.5e-3 at `maxm=8` on an exactly-`D=1` state, and a ~33%
flake in the dispersion tests) had been investigated at length and
attributed to conditioning alone; `_random_raw_tensor`'s docstring in
`vumps.py` concluded "the ground-state solver is not the thing to change".
The amplification was real, but what it amplified was the ~1e-6
eigenvector error; the ground-state solver's stopping criterion was exactly
the thing to change. Re-measured, the error is now 1e-13 to 1e-12 across
`maxm` = 1..16 and no longer monotonic. The regression that pinned the
degradation now pins its absence
(`test_redundant_bond_dimension_no_longer_costs_excitation_accuracy`), as
its own predecessor's docstring had asked for; both docstrings keep the
superseded reasoning rather than deleting it.

## iDMRG on `itensor_version=3`

`Chain::idmrg_ground_state` is a full port of `pyitensor/idmrg.py`, static
observables included. It was for a long time only most of one, and the
three missing pieces are worth knowing about because they are exactly the
parts an energy-only test cannot see: (1) McCulloch's wavefunction
prediction (`theta = lambda_k . B_k . lambda_{k-1}^{-1} . A_k . lambda_k`)
as each micro-step's Krylov start, now `idmrg_wavefunction_prediction` (the
older heuristic it replaced, reusing a raw flat vector whenever its size
matches, survives as the fallback when the previous step's shapes do not
line up yet); (2) extraction of the converged unit cell from a single
micro-step's theta, now `idmrg_theta_cell` + `ic_canonicalize_cell`, kept on
the Chain as `idmrg_cell_`; tiling the per-micro-step `idmrg_U_` chain
instead (whose two ends live in bond bases minted by different micro-steps)
leaves the energy correct while putting `<Sz>` at -0.13 against an exact 0
on the XX chain; (3) the per-site energy baseline subtracted from both
growing environments each micro-step (`idmrg_subtract_energy_baseline`, the
reference `idmrg.h`'s `HL += -energy*IL`), with the density adding the shift
back. On top of (2), `bindings.cc` exposes
`idmrg_onsite_expectation`/`idmrg_two_point_correlator`/
`idmrg_local_excitation_gap`, so `gs_method="idmrg"` supports
`vev`/`correlator`/`local_excitation_gap` on `itensor_version=3`;
`td_dynamical_correlator_window`'s own connected-background subtraction
reads the same cell-based expectation value (it used to tile `idmrg_U_`,
i.e. it had the gauge bug (2) exists to fix).
`tests/test_infinite_chain.py::test_itensor_version3_matches_python_backend`
compares `vev` and a `correlator` sweep as well as `e0`, so the divergence
is guarded rather than merely documented, and
`tests/test_idmrg_correlator_v3.py` covers the surface directly. Still
`itensor_version="python"`-only: `local_excitation_gap(window>0)`
(`local_excitation_gap_windowed` is an explicit prototype, not stable API)
and the iMPS algebra over `IDMRGResult` (`imps_overlap`/`apply_mpo`/
`imps_sum`), exposed since 2026-10-15 as
`Infinite_Many_Body_Chain.overlap`/`apply_operator`/`summps` on
`itensor_version="python"` for both `gs_method`s (the VUMPS `Chain` methods
on v3 are not wired to them yet); `overlap` returns the fidelity per site
as the spectral radius of the mixed transfer matrix, not one eigenvalue,
since a state and its spin-flipped copy give a tied pair of opposite sign.

The iDMRG static observables on v3 are dense arrays, not ITensors: the
`ic_*` helpers in `chain_session.h`, the same design (and for the same
reasons) as the VUMPS `vx_*` ones next to them. Two deliberate departures
from a literal transcription of `idmrg.py`, both pure performance and
neither changing a returned number: the measurement-independent part
(transfer tensors + both families of fixed points) is built once per
ground-state run and cached (`ic_build_cache`) rather than recomputed inside
every `vev`/`correlator` call the way `idmrg.py` does; and a two-point
correlator applies the operator string to the right fixed point one site at
a time (O(chi^4) per site) instead of composing transfer matrices (O(chi^6)
per site). Measured on a critical Heisenberg chain at `maxm=24`: a single
`r=3` correlator went from ~4.1s to ~0.01s. Above `ic_dense_eig_max_` (=64,
matching `idmrg.py`'s own `_DENSE_EIG_MAX`) the transfer-matrix
eigenproblems switch from a dense `zgeev` on the whole `chi^2 x chi^2`
matrix to a matrix-free restarted Arnoldi (`ic_arnoldi_dominant`), the C++
counterpart of `idmrg.py`'s own ARPACK branch. Net effect against
`itensor_version="python"` on that same chain: the growth loop is no faster
(v3 runs ~1.5-2x slower, being ITensor-object-bound, see the benchmark note
in `docs/documentation.md`), while the observables are 5-500x faster.

## A dimerized iDMRG fixed point is not an extraction bug

On the S=1/2 Heisenberg chain, iDMRG at finite `maxm` gives alternating
bonds (n_uc=2, maxm=16: -0.4719 and -0.4143 for 3<SzSz>), and on n_uc=1
`correlator("Sz",0,"Sz",1)` reads -0.4137 against e=-0.4431. The 2026-10-15
audit checked the obvious suspect, a missing `lambda^-1` in `_theta_cell`/
`idmrg_theta_cell`, and it is not there: the growing chain's own central
bond, read straight off each micro-step's theta, alternates -0.4721/-0.4140
step to step, the cell reproduces both values, AKLT (maxm=2) and
Majumdar-Ghosh (two-site cell) tile exactly on both backends, and on n_uc=2
the mean bond equals e to 4e-7 (truncation; unchanged from etol 1e-10 to
1e-13, 1e-8 at maxm=48). The alternation shrinks with maxm. So it is the
iDMRG fixed point, and on n_uc=1 the density itself alternates between
macro-iterations, which is why `converged` stays False there. The response
is a warning, not code: `infinitechain._measure_cell_asymmetry` compares
the Schmidt spectra of the two bonds of the n_uc=1 cell
(`idmrg.cell_bond_weights`/`Chain::idmrg_cell_bond_weights`, gauge
invariant) and `gs_energy`/`vev`/`correlator` warn above an L1 distance of
1e-2 (dimerized: 0.4 to 0.8; uniform: 1e-9 to 3e-4). Do NOT average the
two cell positions for n_uc=1: that reports a bond no bond has.
`tests/test_audit_2026_10_15_infinite.py`.

## The `"python"` iDMRG transfer chain

`idmrg.py`'s transfer chain is `_TransferChain` (2026-09 audit, finding
#20): the per-position site tensors rather than materialized chi^4 arrays,
so `_CorrelatorEnv`, `_canonicalize_periodic`, `imps_overlap`, the window
environment and three of `vumps.py`'s own fixed-point sites allocate no
chi^4 array at all (a materialized chain is still built by the dense
eigensolve below `_DENSE_EIG_MAX`, by the ARPACK fallback and by
`_dominant_eigenvalue_mixed`, all through `_compose_chain`). Measured on a
2-site Heisenberg cell at maxm=64, seeded so the growth trajectory is
identical: the tracemalloc peak of the first `vev` went 517.8 MB to 5.8 MB
and of the growth loop 519.0 MB to 10.4 MB, peak RSS 810.8 MB to 304.1 MB,
growth 23.92 s to 5.23 s, with every returned number unchanged to every
printed digit (both routes on the same state agree to 4.4e-16 in the fixed
points and 5.8e-17 in the observables). The class is deliberately not
list-like, so a leftover `Es[p]` raises instead of quietly materializing
again.

## Open

`pyitensor/idmrg_excitations.py::_deflated_lanczos_run`'s `max(1,|lam|)`
reference is still absolute (2026-09-25b, left open). `"python"`'s Lanczos
at small units was fixed for VUMPS in that pass. `_grow_initial_state` on
`itensor_version=3` still noise-pads where `"python"` subspace-expands.
