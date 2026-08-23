# Fixed issue: `itensor_version=3` VUMPS did not return on a small-gap spinful model

Filed from the same external physics session (a triangulene heavy-fermion / Kondo-lattice
calculation) that reported the `local_excitation_gap` deflation bug — that one is
**fixed**, and so is this one — the file is kept as the record of what the failure looked
like, what caused it, and which regimes it silently affected.

## Summary

On an infinite spinful (c,f) Kondo chain (`Infinite_Many_Body_Chain([1,1],
itensor_version=3)`, two native spinful sites per cell), `gs_method="vumps"` does not
return an energy density even in the **free limit** (`U=J=0`) at `D=8`, `maxiter=60`,
`nrestarts=2`. Confirmed directly in this repository: the run was still at 100% CPU after
**10:22** (mm:ss) with no output and had to be killed. The `itensor_version="python"` VUMPS
returns on the identical model and parameters in **22.7 s** (`e0 = -0.19477349`,
`converged=False` at that iteration budget), and v3's own `gs_method="idmrg"` does `D=8`
on the same model in ~10 s.

So this is not "VUMPS is slow on this model" — the two VUMPS implementations differ by at
least a factor of ~25 on the same input, and the gap between them grows with wall time
since v3 had not returned at all.

## Why it matters beyond speed

`excitation_energies` / `excitation_gap` (the tangent-space quasiparticle ansatz, the
physically principled dispersion) **require** `gs_method="vumps"`. With v3's VUMPS
unusable on this Hamiltonian, the only remaining excitation estimator on that backend is
`local_excitation_gap` (which requires `gs_method="idmrg"`), and that is exactly the route
the reporting session fell back to before hitting the deflation bug. Together with
`docs/known_issue_idmrg_product_state_collapse.md`, no single (backend, solver) pair
covers the whole parameter range of this Hamiltonian.

## Cause (found, fixed)

`Chain::vumps_single_run` materialized the **full dense** `H_AC` and `H_C` every outer
iteration and diagonalized them with `zheev`, where `pyitensor/vumps.py` never forms
either matrix and runs `_lanczos_ground_state` on `_h_ac_action`/`_h_c_action` instead.
`vx_hermitian_ground_state`'s own comment states the assumption behind that choice
explicitly — "D and d_g are always small enough here ... negligible at the bond dimensions
this port targets".

That assumption holds for the models the port was validated on (TFIM/Heisenberg,
one-site unit cell, `D<=3`, so `d_g = 2` or `4` and `n = D*d_g*D <= 36`) and fails badly
for a two-site cell of *native spinful* sites, where `d_g = 4*4 = 16`:

| D | `n = D*d_g*D` | dense cost per outer iteration |
|---|---|---|
| 1 | 16 | trivial |
| 4 | 256 | 256 action calls + `zheev(256)` |
| 8 | **1024** | 1024 action calls + `zheev(1024)` |
| 16 | **4096** | 4096 action calls + `zheev(4096)` |

and that is paid per outer iteration (`maxiter`), per attempt (`nrestarts`), at every step
of the D-ramp. A verbose run confirmed the backend was not hung but grinding: it completed
`D=1` and `D=2` quickly, then `D=4`, and was still inside `D=8` when killed.

**Fixed** by mirroring the Python reference: `Chain::vx_lanczos_ground_state` is a direct
port of `pyitensor/dmrg.py`'s `_lanczos_ground_state` (Lanczos with full
reorthogonalization, early stop when the lowest Ritz value stabilizes) over flat complex
vectors, and `vumps_single_run` uses it for `H_AC`/`H_C` whenever `n` exceeds
`vx_dense_eig_max_` (64), keeping the exact dense path at or below it — so every case the
port was originally validated on stays on the code path it was validated against, and only
the regime that did not work at all changes. The Krylov dimension is `Infinite_Many_Body_
Chain.niter`, the same knob that feeds `pyitensor.vumps.vumps_ground_state`'s own
`niter_lanczos`.

### Measured effect of the fix

Same model, same parameters (`U=J=0`, `D=8`, `maxiter=60`, `nrestarts=2`), same machine:

| build | result |
|---|---|
| before | killed at 10:22, ramp still inside `D=8`, no energy |
| Lanczos above `n=256` | `e0 = -0.19515244` in 180.8 s |
| Lanczos above `n=64` (final) | `e0 = -0.19524372` in **21.7 s** |
| `itensor_version="python"` reference | `e0 = -0.19477349` in 22.7 s |

So the v3 backend now matches the Python reference's runtime on this model, at a slightly
*lower* (variationally better) energy. Both report `converged=False` at `maxiter=60`, which
is what the ~4e-4 spread between them reflects — neither has reached its fixed point at
that iteration budget, and raising `maxiter` is the caller's knob for that.

Note the D-ramp's missing "beat a known-good smaller-D energy" safety-net budget (recorded
in `vumps_ground_state`'s own doc comment as a deliberate simplification versus pyitensor's
driver) was **not** the cause — that net only ever *adds* attempts, so it cannot explain a
non-return. It remains a real difference between the two drivers, and a robustness gap
worth closing separately.

## Reproduction

`vumps_check.py` in the reporting session's own working notes; equivalently, build the
model with `build_infinite_kondo(J=0.0, U=0.0, maxm=8, backend=3, maxiter=60,
nrestarts=2)` (from `paper_figures/generate_idmrg_control.py` of that project), set
`gs_method="vumps"`, and call `gs_energy()`. Compare against `backend="python"`, which
returns in seconds.
