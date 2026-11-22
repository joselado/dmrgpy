---
paths:
  - "src/dmrgpy/dynamics.py"
  - "src/dmrgpy/dynamicstk/**"
  - "src/dmrgpy/kpmdmrg.py"
  - "src/dmrgpy/algebra/kpm.py"
  - "src/dmrgpy/algebra/kpmextrapolate.py"
  - "src/dmrgpy/timedependent.py"
  - "src/dmrgpy/timeevolution.py"
  - "src/dmrgpy/tdz.py"
  - "src/dmrgpy/taskdmrg.py"
  - "src/dmrgpy/distribution.py"
  - "src/dmrgpy/pyitensor/tdvp.py"
  - "src/dmrgpy/pyitensor/gse.py"
  - "src/dmrgpy/pyitensor/task.py"
  - "src/dmrgpy/mpsjulialive/dynamics.py"
  - "src/dmrgpy/mpsjulialive/timedependent.py"
  - "tests/test_audit_*kpm*.py"
  - "tests/test_audit_*realtime*.py"
  - "tests/test_audit_2026_10_task.py"
  - "tests/test_pyitensor_time_evolution_review.py"
  - "tests/test_nh_time_evolution.py"
  - "tests/test_time_evolution.py"
  - "tests/test_kpm*.py"
  - "tests/test_task*.py"
---

# Dynamical correlators and real-time evolution

These are durable rules the hole hunts left, with the mechanism and the
numbers; CLAUDE.md keeps one line per rule. Each rule names the test that
pins it. Do not undo one without reading its finding in `docs/audit_*.md`
first. The file loads when `dynamics.py`, `kpmdmrg.py`, `timedependent.py`, `taskdmrg.py`,
`pyitensor/tdvp.py` or a file under `dynamicstk/` is read or edited. The convention every
submode returns is in the module docstring of `src/dmrgpy/dynamics.py`, and
the rule that every correlator starts from `ground_state_on_session` is in
`.claude/rules/session.md`.

- **KPM broadening is calibrated.** Every KPM route goes through
  `algebra/kpm.py::polynomials_for_broadening`, which solves the Jackson
  kernel's resolution relation so that `delta` is the broadening it is
  everywhere else (FWHM `2*delta` at the band centre); the DMRG routes
  reconstruct from exactly those n moments, the ED route uses the DMRG
  rescaling so the two share the same `x` at the same physical energy (on
  the bandwidth-centred window only: under `kpm_energy_truncate`,
  `mode="ED"` raises), `kpm_n_scale` defaults to 1 and must be a positive
  integer, and
  every backend raises once a moment exceeds 1.5 times the exact bound
  `||vi|| ||vj||`. Two residuals are documented rather than removed: the
  width tightens as `sqrt(1-x^2)` away from the band centre (0.70 of the
  requested width for an intensive excitation of a long chain), and the
  Jackson line peaks about 1.5x higher than a Lorentzian of equal FWHM.
  `get_distribution()`'s own KPM path is deliberately not on the
  calibration (`test_audit_2026_09_24_kpm.py`, `test_audit_2026_09_24b_kpm.py`).
- **`submode="TD"`/`"TDZ"` return the two-sided transform** written in
  terms of forward runs only, `(F[(A,B)] + conj(F[(B^dagger,A^dagger)]))/2`
  (`timedependent.lehmann_density_from_one_sided`), collapsing to `Re F`
  exactly when `canonical.is_dagger_pair` proves `A = B^dagger`; it refuses
  rather than guesses for a name with no known adjoint. Taking `np.real`
  instead measured worse than leaving the transform whole. The Fourier sum
  is evaluated at each requested frequency, not interpolated off an FFT
  grid; `sxt_to_skomega`, the infinite-chain `S(k,omega)` reduction, and
  the two short-window direct callers stay on the FFT stage, and
  `sxt_to_skomega` on the raw transform, since it never sees the pair
  (`test_audit_2026_09_24_realtime.py`).
- **Real-time evolution on `"python"`** (`tdvp._lanczos_expm_multiply`)
  stops on `|coeff|*beta_k*|c_k| < 1e-10` relative to the vector, with
  sub-stepping when `niter` runs out; the host loop and the device path
  share one test, `_krylov_verdict`, and return their basis rather than a
  vector, so the sub-stepping above them cannot drift between the two. Keep
  it that way if either changes. Every `evolve_and_measure_*` loop on every
  backend restores the input norm `norm0` after each step, and
  `gse.global_subspace_expand` normalizes each Krylov vector
  (`test_pyitensor_time_evolution_review.py`). A non-Hermitian MPO takes
  `tdvp._arnoldi_expm_multiply` instead, chosen in `tdvp_step` by
  `mpsalgebra.mpo_is_hermitian(H)` and never by `dt` (METTS and TDZ pass an
  imaginary `dt` with a Hermitian H); the Lanczos route was 3.6e-3 off on a
  non-Hermitian H, since it drops `Im<q|A|q>` and the projections its
  reorthogonalization removes (`test_nh_time_evolution.py`). `TDVP_GSE` stops growing the
  bond dimension after `tdvp_gse_sweeps` steps, so a growing-entanglement
  quench saturates; recorded, not changed.
- **`mode="ED"` real-time routes are `expm_multiply`**, start from the
  cached ground state when called without `wf=`, and
  `evolve_and_measure`/`evolution_ABA` raise on an unknown keyword on DMRG,
  take `h=` on ED and default to `nt=1000` on both
  (`test_audit_2026_09_24c_realtime.py`).
- **TASK's error bar is two numbers, not one.** The two-site variance of
  a pole does not see the part of `B|GS>` outside the one-site tangent
  space, which is dropped before the recursion starts: zero for a sum of
  one-site operators and for a fermion operator on a state of definite
  parity, up to all of the weight for `C_i` on a parity mixture (a
  degenerate Majorana pair) at truncated `maxm`. `taskdmrg.get_poles`
  returns it as `missing_weight` and warns above 1e-6 of the weight
  (`test_audit_2026_10_task.py`).
