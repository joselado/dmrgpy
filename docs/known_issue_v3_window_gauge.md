# Known issue: the `itensor_version=3` IBC window measures in the wrong gauge

**Status: open, and the path is disabled** (found 2026-08-29, while
adding fermionic coverage to the infinite-chain dynamical correlators).
Affects `Infinite_Many_Body_Chain.td_dynamical_correlator` on
`itensor_version=3` only. The pure-Python backend
(`itensor_version="python"`) is unaffected and is exact on the same check.

Because a wrong number is worse than an error, both layers now refuse
rather than answer:

- `Infinite_Many_Body_Chain.td_dynamical_correlator` raises `RuntimeError`
  for `itensor_version=3`, naming `itensor_version="python"` and this
  document.
- `Chain::td_dynamical_correlator_window` itself throws
  `std::runtime_error` unless `Chain::set_allow_defective_window(true)`
  has been called — a deliberate, tests-only opt-in that dmrgpy's own API
  never uses. It exists so the tests that *pin* the defect (and whoever
  fixes it) can still run the code.

Undoing all of that — the raise, the C++ check, the opt-in, this file — is
the last step of the fix, not a separate cleanup.

## The symptom, and why it is not about fermions

`S(x, t=0)` produced by `Chain::td_dynamical_correlator_window` must equal
the static `correlator(...)` for the same operator pair and separation.
That is an exact identity, not a convergence statement: at `t=0` no
evolution has happened, so the window snapshot is just a two-point
function of the converged ground state. On `itensor_version="python"` it
holds to ~1e-15. On `itensor_version=3` it does not:

| model | pair | x | v3 `S(x,0)` | static `correlator` |
|---|---|---|---|---|
| spin-1/2 Heisenberg, n_uc=2, maxm=12 | Sz, Sz | 0 | +0.25000000 | +0.25000000 |
| " | " | +1 | -0.08642967 | -0.16074727 |
| " | " | +2 | -0.01154499 | +0.06045035 |
| " | " | -1 | -0.12969705 | -0.13457189 |
| dimerized spinless fermions, n_uc=2, maxm=12 | N, N | 0 | +0.51396559 | +0.51050298 |
| " | " | +1 | +0.28769256 | +0.02169647 |
| " | Cdag, C | +1 | -0.19434674 | -0.47769573 |
| Heisenberg, n_uc=2, maxm=10 (the example below) | Sz, Sz | -1 | +0.03552110 | -0.13209193 |
| " | " | +2 | -0.04213772 | +0.06037943 |

The last block is `examples/idmrg/td_dynamical_correlator_python_VS_v3/
main.py`, which now runs this identity on both backends and plots it: the
Python side lands at 3.3e-16 across every `x` while v3 misses by up to
1.7e-1 — on a model where the two backends' *energy densities* agree to
6.7e-11, which is exactly why an energy-only comparison never noticed.

The first block is a **spin** chain with a Hermitian, parity-even
operator: no Jordan-Wigner string is involved anywhere. So this is not the
fermionic-string bug that
`tests/test_idmrg_window_fermionic.py` was written for -- it is an
independent, pre-existing defect that the fermionic work merely had a
sharp enough oracle to expose. `x=0` is exact in every row, and the error
grows with `|x|`: the signature of two tensors being contracted in bond
bases that do not match, not of a wrong operator.

## Cause

`Chain::idmrg_build_window` (`mpscpp3/chain_session.h`) tiles `idmrg_U_`,
the raw per-micro-step iDMRG factors, and
`Chain::idmrg_window_snapshot_correlator` builds its bra from the same
`idmrg_U_`. That is exactly the construction the rest of this backend was
fixed away from: `idmrg_U_`'s two ends live in bond bases minted by
*different* micro-steps, so tiling it gives a chain whose energy is right
while its static observables are not (`CLAUDE.md`'s own iDMRG section
records the same failure at `<Sz> = -0.13` against an exact 0 on the XX
chain). `Chain::idmrg_onsite_expectation` /
`idmrg_two_point_correlator` were ported onto the gauge-consistent unit
cell (`idmrg_theta_cell` + `ic_canonicalize_cell`, kept as
`idmrg_cell_`) precisely for this reason, and
`pyitensor/idmrg_window.py`'s own `_window_cell` does the same on the
Python side -- the window path on this backend was simply never moved
over. The inconsistency is visible in the file today: the same method's
*background* subtraction already reads the cell-based
`idmrg_onsite_expectation`, while the snapshot next to it still tiles
`idmrg_U_`.

## Fix

Port `_window_cell`: tile `idmrg_cell_` (with its own per-position left
and right Index objects) in `idmrg_build_window`, build
`idmrg_window_snapshot_correlator`'s bra from the same cell, and take the
closing fixed points from the cell-based `ic_*` cache instead of the
`idmrg_U_`-based `idmrg_all_right_fixed_points()`. Nothing new has to be
invented -- every piece already exists on this backend for the static
observables; it is Index plumbing, not new physics.

`tests/test_idmrg_window_fermionic.py::
test_itensor_version3_window_matches_the_static_correlator_at_t0` is a
`strict=True` xfail pinning the defect: it flips to a failure-to-fail (and
so a visible test failure) the moment the gauge is fixed, which is when
this document should be deleted. Two other tests were rewritten to assert
the refusal rather than a success that had wrong numbers inside it
(`tests/test_idmrg_window_v3.py::test_td_dynamical_correlator_public_api_
refuses_v3` and `tests/test_infinite_chain.py::
test_td_dynamical_correlator_refuses_on_v3_backend`) — both say in their
docstrings what to restore, and that the restored version needs an oracle
this time. That is the real lesson here: the old tests asserted finiteness
and shape, and neither can see a gauge error.

## Scope

- Not affected: `itensor_version="python"` (`gs_method="idmrg"`) for the
  same feature; every static v3 observable (`vev`, `correlator`,
  `local_excitation_gap`), which already uses the cell; `kpm_finite`,
  which builds an ordinary finite chain and never touches this path.
- Affected: `td_dynamical_correlator(itensor_version=3)` for every
  operator, spin included. Its *shape* in `(k, omega)` is qualitatively
  right (which is why `tests/test_idmrg_window_v3.py`'s own loose
  consistency checks pass); the individual `S(x,t)` values are not
  quantitatively trustworthy, and the path now refuses to return them.
  The tests in that module opt in explicitly and assert only machinery
  properties — error handling, `eshift` measurement, shapes, run-to-run
  reproducibility — never that a value is physically right.
