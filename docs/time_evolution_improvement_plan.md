# Real-time evolution: what's worth doing next

Written 2026-09-26, after a review of `itensor_version="python"`'s real-time
evolution (`src/dmrgpy/pyitensor/tdvp.py`, `gse.py`, `tebd.py` and the
`quench_*`/`evolve_and_measure_*` loops in `pyitensor/chain.py`) against exact
propagation and against v3. `ROADMAP.md` tracks what each *backend* has; this
file tracks what the time-evolution *algorithms* could still gain. The same
shape as `docs/idmrg_improvement_plan.md`: ranked by value per unit of effort,
every anchor a real symbol in this checkout.

| # | Item | Effort | Status |
|---|------|--------|--------|
| 1 | `TDVP_GSE`: stop freezing the bond dimension after `tdvp_gse_sweeps` steps | ~1 day | open |
| 2 | v3 Krylov exponentiator: relative `ErrGoal` and `NormCutoff` | ~hours + rebuild | open |
| 3 | Reuse environments across TDVP half-sweeps (a few %) | ~hours | open |
| 4 | Controlled bond expansion TDVP (CBE-TDVP) | ~1–2 weeks | open |
| 5 | Fourth-order composition of the TDVP step | ~1 day | open, unmeasured |
| 6 | Expokit step-size formula for the Krylov sub-stepping | ~hours | open, low value |

What the review found right, and fixed, is recorded in `docs/user_guide.md`
§21 ("The 2026-09-26 review of `"python"` real-time evolution") and pinned by
`tests/test_pyitensor_time_evolution_review.py`: two-site and one-site TDVP are
exact to ~1e-11 at full bond dimension (spin-1/2 with long-range and DM terms,
spin-1, Jordan-Wigner fermions, bosons; real and complex dt) and track ED
below it exactly as v3 does. That is the baseline every item here should be
measured against.

---

## 1. `TDVP_GSE` stops growing the bond dimension after three steps

`tevol_method="TDVP_GSE"` runs a global subspace expansion
(`pyitensor/gse.py::global_subspace_expand`, Yang & White,
arXiv:2005.06104) only before the first `tdvp_gse_sweeps` steps (3 by
default, `manybodychain.py`), then one-site TDVP alone, which conserves the
bond dimension exactly. So on any evolution whose entanglement keeps
growing, the bond dimension freezes at whatever three expansions reached and
the error saturates, whatever `maxm` is. Measured on a 12-site Heisenberg
quench from the Néel state (dt=0.05, t up to 5, largest error of
<Sz_6 Sz_7 + Sz_0> against ED):

| maxm | `"TDVP"` | `"TDVP_GSE"` (python) | `"TDVP_GSE"` (v3) |
|---|---|---|---|
| 8  | 2.5e-2 | 2.3e-2 | 2.3e-2 |
| 16 | 1.0e-2 | 3.8e-3 | 2.5e-3 |
| 32 | 4.5e-4 | 1.1e-3 | 2.0e-3 |
| 64 | 8.7e-7 | 1.9e-3 | 2.2e-3 |

The schedule copies `mpscpp3/TDVP/sample/run.cc` (`if(n < 3) addBasis(...)`),
a demonstration. The package's own `mpscpp3/TDVP/README.md` recommends the
opposite once the bond dimension is large enough: turn the expansion off
**and switch to two-site TDVP**. Options, in order of simplicity:

- after `tdvp_gse_sweeps` steps, continue with `num_center=2` instead of 1
  (the README's recommendation; one line per loop);
- expand every `k` steps rather than only the first few;
- expand when needed, e.g. whenever a bond is below `maxm` and its smallest
  Schmidt value is above the cutoff.

Anchors, all of which must move together, since the loops are kept in
lockstep across backends: `pyitensor/chain.py::quench_tdvp_gse`/
`evolve_and_measure_tdvp_gse`, `mpscpp3/chain_session.h`'s
`quench_tdvp_gse`/`evolve_and_measure_tdvp_gse`, `mpsjulialive/tdvp.jl`, and
`tdz.py` (TDZ drives the same expansion step by step). This changes numbers
for every `TDVP_GSE` run, so it needs a user-guide note and a regression that
pins the non-saturation (for example the table above at maxm=64).

## 2. v3's Krylov exponentiator still has absolute thresholds

`pyitensor/tdvp.py::_lanczos_expm_multiply` is scale-free since 2026-09-26.
v3's counterpart, ITensor's `applyExp` (`mpscpp3/ITensor/itensor/
iterativesolvers.h`), is vendored and should stay unedited, but its two
thresholds are arguments: `ErrGoal` (1e-10, compared with an estimate that
carries `nrm`, the vector's norm) and `NormCutoff` (1e-7 on the Lanczos
beta, which carries the units of H). By reading, `TDVP/tdvp.h` forwards its
`args` to both `applyExp` calls (lines ~231 and ~263), so
`Chain::tdvp_step` could pass `"ErrGoal",1e-10*norm(psi)` and a
`"NormCutoff"` scaled by the Hamiltonian's own unit (its largest coefficient,
which `hscale_up_`'s machinery already computes). The 2026-09-25 record
measured the v3 side at 1.04e-4 of the peak for a TD correlator at an
operator scale of 1e-9, and 6.97e-6 at a Hamiltonian scale of 1e-6. Verify
the forwarding before relying on it; needs a rebuild of `mpscpp3`.

## 3. Reuse environments across TDVP half-sweeps

Each `tdvp_step` builds all right environments at the start of
`_half_sweep_lr` (`_all_right_environments`) and all left environments at the
start of `_half_sweep_rl` (`_all_left_environments`), although the preceding
half-sweep has just built exactly those, for the final tensors, one site at a
time. The left-to-right sweep's `left_env` dict is what the right-to-left
sweep needs; the right-to-left sweep's `right_env` dict is what the next
step's left-to-right sweep needs, as long as the caller changes nothing but
the center tensor in between (the `_renormalized` scaling touches site 1
only, which no right environment contains). Profiled on a 30-site
Heisenberg chain at maxm=64, building environments one site at a time
inside the sweeps is ~9% of a step and the two full rebuilds did not reach
the profile's top 25 entries (under ~4%), so this saves a few per cent: a
cheap, exact change, not a large one. Keep the bra/ket link relabelling (`dmrg.py::_relabel_bra_local`)
consistent: an environment is only valid while the link indices it was built
against are still the ones in the state.

## 4. Controlled bond expansion TDVP (CBE-TDVP)

J.-W. Li, A. Gleis and J. von Delft, arXiv:2208.10972 (Phys. Rev. Lett. 133,
026401 (2024)): one-site TDVP whose bond dimension grows on the fly, by
adding the directions of the two-site tangent space's orthogonal complement
that the Hamiltonian actually couples to, chosen by a cheap projection
rather than a full two-site SVD. It has one-site cost with two-site-like
accuracy, so it would replace both `"TDVP"` (two-site, whose local matvec is
a factor `d` more expensive than one-site's) and `"TDVP_GSE"` (whose MPO
applications for the Krylov vectors are its own cost). The gain is largest
at `d=4` (native Hubbard sites, spin-3/2) and for long-range Hamiltonians,
where two-site TDVP's projection error is largest. The same expansion step
serves ground-state DMRG (the authors' arXiv:2207.14712), so it could be
written once in `pyitensor/` and used by both `dmrg.py` and `tdvp.py`.
Largest item here, and a new `tevol_method`, so it needs the user-guide
section, a cross-check against `"TDVP"` and ED, and a decision on whether v3
gets a port.

## 5. Fourth-order composition of the TDVP step

`tdvp_step` is the symmetric second-order projector-splitting integrator (a
left-to-right half-sweep by dt/2 and its mirror). Composing three of them
with Yoshida's coefficients (H. Yoshida, Phys. Lett. A 150, 262 (1990)),
`dt1 = dt/(2-2^(1/3))`, `dt2 = -2^(1/3) dt/(2-2^(1/3))`, `dt1` again, raises
the order of the time-step error to four for three times the cost per step.
It does nothing for the projection and truncation errors, which usually
dominate below full bond dimension, so it pays only where the time-step
error is what limits accuracy (small systems at full bond dimension, long
times at large dt). The negative middle step is harmless for unitary
evolution but not for imaginary time (METTS, TDZ), so it would have to be
opt-in and refused there. Unmeasured; measure the error-versus-cost curve
before adopting it. See Paeckel et al., "Time-evolution methods for
matrix-product states", arXiv:1901.05824 (Ann. Phys. 411, 167998 (2019)),
for the error budget of every method in this file.

## 6. Expokit's step-size formula for the Krylov sub-stepping

When the `niter` budget runs out, `_lanczos_expm_multiply` now takes the
largest dyadic fraction of the step that converges on the basis already
built (`_largest_converged_substep`) and continues in sub-steps of that
size. Expokit (R. B. Sidje, ACM TOMS 24, 130 (1998)) instead predicts the
next step size from the error estimate's own scaling, which wastes fewer
matvecs: at |coeff| times the spectral width = 60 and `niter=50` the halving
takes 88 matvecs where one run of ~60 would do. It only matters for steps far
too long for the budget, which no ordinary run makes, hence the low priority.
