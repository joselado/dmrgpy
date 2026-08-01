# Hybrid dynamical-correlator method — design plan

**Status:** proposed, not yet implemented. Written 2026-08-01 as a plan for
combining the existing dynamical-correlator submodes (KPM, TD, EX, CVM/INV,
ROOTN, TDZ, maxent) into a hybrid method that is more accurate or cheaper
than any single one alone. Nothing described here has been built; this
document is the plan to build from when the work is picked up.

## Terminology

`kpm`/`td`/`ex` map directly to `submode="KPM"/"TD"/"EX"` (both backends).
`inv` is a literal existing submode, but ED-only
(`src/dmrgpy/edtk/dynamics.py`, `submode="INV"`): a dense/sparse exact
inversion of `(z-H)` at each &omega;. DMRG has no `"INV"` — its nearest
analogue is `submode="CVM"` (`src/dmrgpy/cvm.py`), the iterative (CG)
version of the same resolvent, MPS-truncated instead of dense. ED also
exposes `submode="CVM"` as the CG version of its own `INV`
(`edtk/dynamics.py::dynamical_correlator_inv(mode="cv")` vs
`mode="full"`) — so "INV vs CVM" is really "exact dense solve vs
iterative solve of the identical linear system," already present in
miniature on the ED side. Also present: `ROOTN` (`rootndmrg.py`), `TDZ`
(`tdz.py`), `maxent` (`distribution.py`), and NH-KPM for non-Hermitian
`H` (`nonhermitian/kpm.py`).

## Cost/accuracy landscape

Every submode computes the same `G_AB(omega)`, but spends cost along a
different axis:

| Method | Cost driver | Shape | Per-point sharpness |
|---|---|---|---|
| **KPM** | # Chebyshev moments `N ~ scale/delta` | `O(N)`, shared across the *whole* window | Jackson-kernel broadened (trades resolution to kill Gibbs ringing) |
| **CVM** ("inv", iterative) | CG iterations x # frequency points | `O(n_pts*iters)`, independent of window width | limited only by `cvm_maxm`'s truncation floor |
| **INV** (ED, dense) | one `O(dim^3)` inverse per point | `O(n_pts*dim^3)` | numerically exact (only the artificial `eta`) |
| **EX** | # DMRG excited states `nex` | `O(nex)` once, reused for every omega | exact poles, but only within the truncated subspace |
| **TD/TDZ** | simulated time `T ~ damping_periods/delta` | `O(T/dt)`, amortized FFT over the window | Lorentzian(delta), no Gibbs ringing |

KPM's cost to resolve a feature of width `delta` inside bandwidth `W` is
`~W/delta` matvecs, paid once, covering the entire window. CVM's cost to
resolve `n_local` points around one feature is `~n_local*iters`,
independent of `W`, but paid per feature. A method that uses KPM to find
*where* structure is, then spends CVM's per-point budget only there,
costs `~N_coarse + n_features*n_local*iters` — below either extreme
whenever the spectrum is "mostly smooth with a few sharp features" (the
common case for gapped magnets/insulators, and even many gapless cases
away from the gap).

Existing prior art in this direction, worth knowing about before
extending it: `kpm_extrapolate` already boosts KPM moments post-hoc
(`examples/dynamical_correlator/dynamical_correlator_extrapolate_moments`);
`kpmdmrg.get_dynamical_correlator` even has a `deconvolve=` kwarg for
exactly this kind of post-hoc sharpening, but it is currently dead code
(accepted as a parameter, never read anywhere in the function body) — a
landing spot that already half-exists. `dynamical_correlator_sharpen`
(the example of that name) instead gets its sharpening entirely from
`kernel="lorentz"`, a single-method resolution trick, not a hybrid
across solvers.

## The core principle

Cheap-global methods (KPM, short TD) locate features but blur them
(kernel/damping broadening). Expensive-local methods (CVM, INV, ROOTN)
resolve a feature exactly (up to the user's own `delta`) but only pay off
per point. A hybrid should run the cheap method first to find *where* to
look, then spend the expensive method's budget only there.

## Phase 1 (recommended starting point): `submode="HYBRID"` — coarse locate, fine refine

Implement as a **generic combinator**, not a KPM/CVM-specific script: a
new function that calls the existing public
`get_dynamical_correlator(mode=..., submode=...)` twice — once
cheaply/globally, once expensively/locally — and splices the results.
Built entirely on top of the already mode-dispatching public API, so it
works on every backend either chosen submode already supports, with zero
new tensor-network code.

```python
def dynamical_correlator_hybrid(self, mode="DMRG", name=None, es=None,
        delta=1e-1, coarse="KPM", refine="CVM",
        coarse_kwargs=None, refine_kwargs=None,
        prominence=0.1, window_halfwidth=None, refine_npoints=15,
        refine_at=None, **kwargs):
    delta_coarse = (coarse_kwargs or {}).pop("delta", 4*delta)
    x_c, y_c = self.get_dynamical_correlator(mode=mode, submode=coarse,
            name=name, es=es, delta=delta_coarse, **(coarse_kwargs or {}))

    windows = find_peak_windows(x_c, y_c, prominence=prominence,
            halfwidth=window_halfwidth or 10*delta)
    if refine_at is not None:
        windows += [(w - (window_halfwidth or 10*delta),
                     w + (window_halfwidth or 10*delta)) for w in refine_at]

    y = y_c.copy()
    for lo, hi in merge_overlapping(windows):
        x_r = np.linspace(lo, hi, refine_npoints)
        _, y_r = self.get_dynamical_correlator(mode=mode, submode=refine,
                name=name, es=x_r, delta=delta, **(refine_kwargs or {}))
        mask = (x_c >= lo) & (x_c <= hi)
        y[mask] = blend(x_c[mask], y_c[mask], x_r, y_r)  # raised-cosine splice
    return x_c, y
```

(All numeric defaults above — `4*delta`, `10*delta`, `refine_npoints=15`
— are illustrative starting points, to be tuned once real timing/accuracy
data exists.)

- **File/integration**: new `src/dmrgpy/hybriddc.py` (sibling to
  `cvm.py`/`kpmdmrg.py`/`dcex.py`), intercepted in
  `manybodychain.py::get_dynamical_correlator` (currently around line
  466) *before* the `mode = self.get_mode(...)` DMRG/ED split, so it is
  mode-agnostic by construction — `coarse="KPM", refine="INV"` works on
  the ED backend for free, with the same code.
- **Peak detection**: `scipy.signal.find_peaks` on `|y_c|` with a
  `prominence` threshold; `refine_at=[...]` lets a caller force a window
  manually when the physics is already known (e.g. a known bound state),
  the same override-over-heuristic pattern `kpm_scale` already follows
  elsewhere in this codebase.
- **Validation**: mirror the existing pattern in
  `tests/test_dynamical_correlator.py`, which already has
  `test_cvm_dynamical_correlator_peak_matches_exact_gap`/
  `test_kpm_..._peak_matches_exact_gap` against the 4-site Heisenberg
  `HEISENBERG_4_GAP=0.658919` fixture and ED `"INV"` golden values; add
  `test_hybrid_dynamical_correlator_peak_matches_exact_gap` the same way.
  Add `examples/dynamical_correlator/dynamical_correlator_hybrid_kpm_cvm/main.py`
  following `dynamical_correlator_VS_ED`'s template, with a genuine
  **wall-clock-matched accuracy comparison** (hybrid vs plain KPM at the
  same total time, hybrid vs CVM-everywhere at the same accuracy) — that
  head-to-head is the actual point; a bare correctness check alone would
  not justify "cheaper."
- **Honesty**: document this as a heuristic accelerant (in the spirit of
  this codebase's existing "finite-window KPM is an approximation"
  disclosure for iDMRG), not a replacement for a manual CVM sweep on an
  unfamiliar spectrum — a real feature can hide below `prominence`, and a
  Gibbs-ringing artifact from the coarse pass can trigger a wasted refine
  window. Both are cheap to detect empirically (compare against a plain
  fine CVM sweep once per new physical regime).
- **Docs**: once implemented, add a `submode="HYBRID"` paragraph to
  `docs/user_guide.md`/`.tex` section 6, right after its existing
  "Choosing a method" paragraph, per CLAUDE.md's rule to document every
  new dynamical-correlator submode there.

## Phase 2: EX pole deflation feeding KPM/CVM's residual

For gapped systems with a few sharp poles plus a weaker continuum: run
`EX` first (cheap), then subtract each pole's own **closed-form**
Chebyshev moments (`mu_m^pole = w * T_m(eps_pole)`, `eps_pole` = the
pole's rescaled energy) from KPM's moments before reconstructing, so
KPM's finite moment budget is not spent reproducing peaks already known
exactly. Final answer = EX's exact poles + KPM's now-cleaner residual.

Needs a small refactor: `kpmdmrg.get_dynamical_correlator`'s C++ path
currently returns only the already-reconstructed spectrum, not the raw
`moments`, unlike `general_kpm_moments`, which already separates "get
moments" from "reconstruct from moments" for the generic-operator path —
that split is the template to follow. Main risk is the `scale`/`shift`
rescaling convention: this codebase has hit real bugs here before (the
KPM energy-truncation window-anchoring fix, see
`docs/user_guide.md` section 6's "Energy truncation" paragraph), so
budget real care for correctly mapping EX's exact `E_n` into the same
rescaled Chebyshev argument KPM uses, not just for the deflation math
itself (which is standard and low-risk). Build this after Phase 1, whose
window/splice helpers Phase 2 reuses for placing the EX poles back in.

## Phase 3 (optional layer on Phase 1): disagreement-triggered refinement

Replace or augment Phase 1's single-method peak-finder with a
cross-check: run coarse KPM and a short/cheap TD pass over the same
window; where they agree within tolerance, trust the cheap answer; where
they disagree, that disagreement *is* the refine trigger (feeds the same
`windows` list Phase 1 already builds). This turns "where do I refine"
from a heuristic into an actual two-estimator error signal — KPM's
kernel-truncation bias and TD's finite-time-window bias are different, so
agreement is real evidence. Best treated as a plug-in alternative
locator, not a separate pipeline.

## Phase 4 (stretch, higher risk): KPM-seeded CVM warm start

CVM's CG currently cold-starts from `b = -eta*B|GS>`. KPM's Chebyshev
vectors `T_m(H~)|B GS>` are a much better implicit approximation to the
resolvent — summed with the right `z`-dependent coefficients, they *are*
a partial reconstruction of the correction vector. Using that as CVM's
initial guess could cut CG iterations substantially. This is more
invasive than Phases 1-3: the C++ backends (`chain_session.h`, both
v2/v3) and pyitensor currently discard the per-moment vectors inside the
recursion, so surfacing them means touching the KPM recursion itself, not
just Python glue.

Flag explicitly against the prior finding recorded for this codebase
(`rootn_correction_vector_per_bond_sweep_fails` /
`dynamical_correlator_cvm_cg_limits`): warm-starting CVM's CG from a
*neighboring frequency's* converged correction vector was tried and
measurably hurt (truncated CG can stagnate at a much worse residual than
a cold start; see `cvm.py`'s own `dynamical_correlator` docstring). This
proposal is a *different* warm-start source (same frequency, KPM-derived,
not a neighboring frequency's CVM solution) — plausibly fine, but it
needs the same empirical rigor (iteration-count comparison on the same
test chains used before) before being trusted, not an assumption that it
inherits the earlier negative result just because both involve
"warm-starting CVM."

## Suggested build order

1. **Phase 1** — pure Python glue over already-working, already-tested
   public APIs, no C++/pyitensor changes required, and its
   window/splice/blend helpers are reused by every later phase.
2. **Phase 2** — next-highest additional value, given how many example
   models in this repo are gapped few-pole-plus-continuum systems
   (`examples/topological/`, `examples/kondo/`, `edge_mode_haldane/`).
3. **Phases 3-4** — once Phase 1 has real timing/accuracy data showing
   where cost is still being left on the table.
