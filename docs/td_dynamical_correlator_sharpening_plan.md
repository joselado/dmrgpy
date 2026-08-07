# TD dynamical-correlator lineshape sharpening — design plan

**Status:** Track A (pluggable damping/window family) implemented. Track B
(linear-prediction extrapolation) implemented. `dynamical_correlator`
(submode="TD")'s own defaults changed to `damping="exp"` (unchanged) +
`predict=True` (was `False`) based on the empirical comparison below.
Written 2026-08-07.

## Choosing the default combination

Once both tracks existed, the empirical question was which combination
to default `submode="TD"` to, so it "looks sharper like KPM" out of the
box rather than only on request. Measured on the 4-site Heisenberg chain
also used by `tests/test_dynamical_correlator.py`
(`HEISENBERG_4_GAP=0.658919`, `delta=0.05`, single isolated resonance):

| combination | FWHM | peak position |
|---|---|---|
| KPM (Jackson kernel), for reference | 0.082 | 0.6588 |
| `exp`, no predict (old default) | 0.205 | 0.6810 |
| `exp` + predict | **0.171** | **0.6564** |
| `gaussian` + predict | 0.207 | 0.6564 |
| `parzen` + predict | 0.173 | 0.6564 |

`predict=True` alone (any damping choice) fixes most of the peak-position
error (`0.6564` vs the exact `0.6588`, all three tied, vs `0.6810` for
the old undamped-window default). For FWHM, `exp`+predict and
`parzen`+predict are statistically tied for narrowest (`0.171` vs
`0.173`) and both clearly beat `gaussian`+predict (`0.207`, no better
than the old default) -- consistent with `gaussian`'s own larger
intrinsic FWHM at fixed `delta` (`2.35*delta` vs the Lorentzian's
`2*delta`, see `_damping_window`'s docstring) largely canceling out
prediction's narrowing. `exp` was kept as the paired default over
`parzen` for simplicity, since the two are statistically indistinguishable
here -- once `predict=True` extends the effective window by
`lp_extend_factor` (default 10x), the residual truncation-ringing effect
`parzen` specifically targets is already small relative to `exp`'s own
resolution at that much longer effective `Tmax`, so `parzen`'s extra
mechanism adds little on top of prediction alone. `gaussian`/`parzen`
remain available as opt-in choices (e.g. `gaussian` when far-tail
suppression matters more than peak height/width for a specific plot).

This comparison was only run for `submode="TD"` itself
(`dynamical_correlator`); `submode="TDZ"` and the internal `S(k,omega)`
reduction (`sxt_to_skomega`) share the same `_fourier_transform_correlator`
mechanics and accept the same kwargs, but keep their own prior defaults
(`damping="exp"`, `predict=False`) since the comparison wasn't repeated
for their different correlator shapes (TDZ's complex-time-reconstructed
`C(t)`, or `S(k,omega)`'s per-k reduction) -- revisit if/when that data
exists.

Neither knob touches the tail's *asymptotic* algebraic decay far from any
resonance (see `examples/dynamical_correlator/dynamical_correlator_td_sharpen`'s
right-hand panel): `"TD"` (any combination here) still sits several
orders of magnitude above `"KPM"`'s Jackson-kernel tail at `delta=0.05`
over a several-J-wide window. These knobs sharpen peaks and suppress
near-tail ringing; they don't close that specific gap, which is intrinsic
to windowing a single exponential/Gaussian/compact-support taper rather
than a polynomial-kernel reconstruction.

## Problem

`submode="TD"` (`timedependent.py::dynamical_correlator`) windows the raw
real-time correlator `C(t)` with `exp(-delta*t)` before the FFT
(`_fourier_transform_correlator`). Multiplying by a pure exponential is
mathematically equivalent to convolving the spectrum with a Lorentzian of
width `delta`, and Lorentzians have heavy algebraic
(`~1/(omega-omega0)^2`) tails. That is the "long tail" reported for "TD"
compared to KPM's default reconstruction, which uses the Jackson kernel
(`algebra/kpm.py::jackson_kernel`) -- specifically derived to decay much
faster away from a peak while still suppressing Gibbs ringing. `kpm.py`
already exposes `kernel="lorentz"` alongside `kernel="jackson"` for
comparison, confirming "Lorentzian tail" vs "Jackson-kernel tail" is
exactly the tradeoff at issue, just on the KPM side of the codebase
already.

This plan is narrower than `docs/hybrid_dynamical_correlator_plan.md`
(which combines whole submodes, e.g. KPM-coarse + CVM-refine): it only
tightens TD's own lineshape. The two are complementary -- the hybrid
plan's Phase 3 wants a cheap TD error signal to cross-check against KPM,
and a sharper TD directly improves that signal.

## Literature

- White & Affleck (2008); Barthel, White & Schollwöck, *"Spectral
  functions in one-dimensional quantum systems at T>0"*, PRB 79, 245101
  (2009), [arXiv:0901.2342](https://arxiv.org/abs/0901.2342) -- introduced
  **linear prediction**: fit `C(t)` with a linear autoregressive model on
  its computed tail, then extrapolate far beyond the directly simulated
  `Tmax` before windowing+FFT. The standard DMRG-dynamics fix for
  truncated-time-window artifacts, since it sharpens resolution without
  more real TDVP time (i.e. without more bond-dimension growth).
- Wolf, Justiniano, McCulloch & Schollwöck, *"Spectral functions and time
  evolution from the Chebyshev recursion"*, PRB 91, 115144 (2015),
  [arXiv:1501.07216](https://arxiv.org/abs/1501.07216) -- proves KPM/
  Chebyshev-recursion and real-time-evolution+FFT are the same
  calculation in a certain limit, and discusses linear prediction on both
  sides. The direct theoretical bridge between "why KPM looks sharper"
  and "how to fix TD."
- Tang, Jia, Moritz & Devereaux (2025), *"Improving Spectral Resolution
  from Real-time Evolution for Correlated Systems"*,
  [arXiv:2509.15539](https://arxiv.org/abs/2509.15539) -- a recent
  autoregressive/ML variant of the same idea; confirms extrapolation is
  still the state of the art for this problem.
- Weisse, Wellein, Alvermann & Fehske, *"The kernel polynomial method"*,
  Rev. Mod. Phys. 78, 275 (2006) -- the Jackson-kernel review already
  informally cited at `docs/user_guide.md`'s KPM section. Per followup
  literature, kernel-style damping **underestimates** peak height while
  linear prediction **overestimates** it -- the two techniques bound the
  true answer from opposite sides, a useful built-in sanity check now
  that both exist here.
- **TeNPy** (GPL-3.0, same license as dmrgpy) ships a `SpectralSimulation`
  class with built-in `"linear_prediction"` and `"gaussian windowing"`
  post-processing for exactly this time-domain-to-spectrum step
  ([docs](https://tenpy.readthedocs.io/en/v1.0.2/reference/tenpy.simulations.time_evolution.SpectralSimulation.html),
  [github](https://github.com/tenpy/tenpy)) -- confirms both tracks below
  are established practice, not experimental; consulted for API/design
  only, no code copied (dmrgpy's TD pipeline has its own shape). The
  windowed-FT literature for tensor-network correlators also documents a
  **Parzen window** for the same truncation-ringing problem, independent
  of the peak-broadening choice.

## Track A — pluggable window/damping family (implemented)

`timedependent.py::_damping_window(ts, delta, damping="exp")` replaces the
hardcoded `exp(-delta*t)` with a selectable family, the time-domain
analogue of `kpm.py`'s `kernel="jackson"/"lorentz"/"plain"` pattern:

- `"exp"` (default, unchanged behavior): Lorentzian tail.
- `"gaussian"`: `exp(-(delta*t)**2/2)`, Gaussian tail (`exp(-omega^2)`
  decay, much faster than Lorentzian's `1/omega^2`).
- `"parzen"`: Parzen window (compact support, vanishing derivative at
  `Tmax`) multiplied on top of the exponential, targeting the truncation-
  ringing artifact specifically rather than the peak-broadening shape.

Wired through `dynamical_correlator` (submode "TD"), `sxt_to_skomega`
(`S(k,omega)`), and `tdz.py::dynamical_correlator_tdz` (submode "TDZ") --
all three share `_fourier_transform_correlator`, so one change covers
every real-time-evolution-based submode. `damping` itself defaults to
`"exp"` everywhere (unchanged) -- only `dynamical_correlator`'s `predict`
default changed (see Track B and "Choosing the default combination"
above); `"gaussian"`/`"parzen"` stay opt-in.

## Track B — linear-prediction extrapolation (implemented)

`dynamicstk/linearprediction.py::linear_predict_extend` fits an
autoregressive (AR) model to the tail of the raw, undamped `C(t)` (an
SVD-regularized least-squares solve of the linear-prediction normal
equations, following Barthel et al.'s appendix -- plain Yule-Walker can
be ill-conditioned; uses only `numpy`, no new dependency) and extrapolates
`C(t)` to a longer effective `Tmax` before Track A's windowing+FFM runs.
Exposed as `predict=`/`lp_order=`/`lp_extend_factor=`/
`lp_fit_start_fraction=`/`lp_max_pole_radius=` on
`dynamical_correlator`/`dynamical_correlator_tdz`/`sxt_to_skomega`.
`dynamical_correlator` (submode "TD") defaults `predict=True` (see
"Choosing the default combination" above); `dynamical_correlator_tdz`
("TDZ") and `sxt_to_skomega` keep `predict=False`. Because extrapolation
happens on the raw signal before any deliberate damping, the extended
series can then use the same (or a smaller) `delta` for a sharper line at
no extra real TDVP cost -- the reason the literature above treats this as
the primary fix, with windowing (Track A) as a complementary, cheaper
partial remedy.

## Out of scope here

`docs/hybrid_dynamical_correlator_plan.md`'s Phase 3 (KPM/TD
disagreement-triggered refinement) stays in that separate plan; noting
only that a sharper TD (this plan) directly improves the quality of the
cross-check signal it depends on.
