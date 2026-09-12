# Fixed issue: `itensor_version=3` VUMPS returned energies BELOW the exact variational minimum

**Status**: FIXED on 2026-09-12, on both solvers, by porting the Python side's
*ordering* as `Chain::vx_choose_fixed_point` -- see "The fix" at the bottom,
which also carries the before/after measurements and the second, latent defect
this turned up in `vx_bond_fixed_points` itself. The file is kept as the record
of what the failure looked like, because the reasoning that localized it is the
part worth keeping: everything above "The fix" is in the present tense as it was
written while the issue was open.

It affected `itensor_version=3` only, on `gs_method="vumps"` (the default for
`Infinite_Many_Body_Chain`), on both the grouped and the sequential solver.
`itensor_version="python"` was unaffected -- this was the C++ counterpart of the
2026-09 audit's finding #3, which was found and fixed only on the Python side.

## What happens

A variational method minimizes `<psi|H|psi>` over normalized states, so the
energy density it reports can be above the exact one (not converged) but never
below it. `Chain::vumps_ground_state` / `Chain::vms_ground_state` intermittently
return an energy below it, through the ordinary public
`Infinite_Many_Body_Chain.gs_energy()` and with nothing downstream flagging it.
(The Python failure this mirrors also reported `converged=True`; whether v3's own
flag does the same was not checked here, only the returned energy was.)

The model is the same field-polarized chain
`tests/test_vumps_redundant_bond_dimension.py` uses, where the exact answer is
known in closed form and the exact state has bond dimension 1:

```
H = -4.0 * sum_i Sz_i + 0.7 * sum_i Sz_i Sz_{i+reach}
exact energy density  e = -4.0/2 + 0.7/4 = -1.825   (<Sz> = 1/2 everywhere)
```

so any `e < -1.825` is impossible, and the size of the excursion is a direct
read of how wrong the environment was.

## What was measured

The script is reproduced in full at the bottom of this file, and every table row
and every Observed block below is one run of exactly that script -- no other
variant of it was used. Every run thread-pinned
(`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1`, `taskset -c 11`, one core) with this checkout's own
`src/` forced onto the path, on the `_dmrgcpp*.so` built 2026-09-12 13:28 -- i.e.
*after* the audit's own C++ fix (#7) and its rebuild, not on the stale
extension the finding was first seen on.

`ic.maxm = D`, `ic.maxiter = 300`, `ic.etol = 1e-12`, `ic.vumps_nrestarts = 2`.
"below (>1e-12)" counts any run returning `e < -1.825 - 1e-12`; "below (>1e-9)"
counts the subset that also exceeds the `abs=1e-9` tolerance the existing
regression test asserts at, i.e. the subset that would actually fail the suite.

GROUPED solver (`n_uc=1`, `reach=1`, the reach-1 one-site cell):

| backend | D | runs | below (>1e-12) | below (>1e-9) | raised | worst |
|---|---|---|---|---|---|---|
| 3 | 4 | 80 | 0 | 0 | 0 | -- |
| 3 | 6 | 80 | 13 | 2 | 0 | 2.33e-08 |
| 3 | 8 | 80 | 18 | 3 | 0 | 6.36e-08 |
| python | 6 | 40 | 0 | 0 | 1 | -- |
| python | 8 | 40 | 0 | 0 | 0 | -- |

SEQUENTIAL solver (`n_uc=1`, `reach=2` -- a one-site cell whose only coupling
reaches two sites, so `vumps_is_reach_one` is false and `use_multisite` routes
it to `vms_ground_state`):

| backend | D | runs | below (>1e-12) | below (>1e-9) | raised | worst |
|---|---|---|---|---|---|---|
| 3 | 4 | 30 | 6 | 2 | 0 | **6.13e-04** |
| 3 | 6 | 30 | 5 | 1 | 0 | **2.85e-03** |
| 3 | 8 | 30 | 3 | 2 | 0 | 6.73e-08 |
| python | 6 | 40 | 0 | 0 | 0 | -- |

The `raised` column is the script's own count of runs that ended in an
exception instead of a number. It is 0 everywhere except one `"python"`
grouped D=6 run in 40, which raised `numpy.linalg.LinAlgError`; a further 60
runs of that same cell raised none, so it is rare (~1 in 100) and it is a
different failure from the one this file is about -- not investigated here.

Independent measurements of the same cells, at different times, agree on the
behaviour and not on the rate. Besides the table above: one verifier measured
the sequential cell at D=6 reaching 1.38e-05 and D=8 7.28e-04 below over 6 runs
each (the numbers recorded in that test file's own `xfail` comment; the driver
and threshold behind them are not stated there); another, running exactly the
script below, got the
grouped cell at 10 of 80 (D=6, one past 1e-9, worst 3.13e-08) and 16 of 80 (D=8,
one past 1e-9, worst 2.60e-08), and the sequential cell at 5 of 30 (D=6, two
past 1e-9, worst 1.50e-03) and 1 of 30 (D=8, worst 3.61e-07). The counts move
run to run because nothing here is reproducible (see the note on the
eigensolver's start vector under "What a fix would involve") and
because the tail is heavy: the frequent excursions are 1e-12..1e-10 and the rare
ones are four to six orders of magnitude larger. Do not quote a single rate;
quote the threshold with it, and expect the >1e-12 and >1e-9 columns to differ
by roughly an order of magnitude in count.

## Why this is the same defect class as audit finding #3, not slow convergence

The obvious alternative explanation is that `etol=1e-12` simply does not pin the
energy read-off to better than ~1e-10, so the sub-1e-9 excursions are gauge
noise with a sign. Three things rule that out:

1. **`itensor_version="python"` never does it**, at the same `etol`, the same
   `maxiter`, the same model and the same D. On the grouped cell, 0 of 40 at
   D=6 and 0 of 40 at D=8 where v3 is 13 of 80 and 18 of 80; on the sequential
   cell, 0 of 40 at D=6 where v3 is 5 of 30. A tolerance floor would be a
   property of the criterion, not of the backend.
2. **The tail is not a tolerance.** 2.85e-03 below the exact minimum, on a model
   whose exact state is a product state, is nine orders of magnitude past any
   convergence explanation.
3. **The signature matches.** Finding #3's Python failure had exactly this
   shape: the state stayed right (`<Sz>` reads exactly 0.5 throughout) while the
   energy read-off went below the minimum, because the dominant fixed point of
   the transfer matrix is degenerate under redundant bond dimension and an
   arbitrary element of that degenerate subspace carries its own energy. The
   exact state here needs D=1, so every D>1 run is in exactly that regime --
   and every cell measured at D>1 shows excursions except the grouped one at
   D=4 (0 of 80). Note the size of the excursion does NOT grow monotonically
   with the surplus: on the sequential cell the worst tail is at D=4 and D=6
   (6.1e-04, 2.8e-03) and D=8 is three orders tighter, so the rate and the tail
   are properties of which element of the degenerate subspace a given run
   happens to land on, not of how much surplus there is.

## What the C++ does differently from the fixed Python

`mpscpp3/chain_session.h` already has the right primitive --
`Chain::vx_bond_fixed_points`, the fixed points the state's own `C` names. What
it does not have is the Python fix's *ordering*. Both C++ environment builders
are **eigensolver-first**:

```cpp
// vumps_build_environments, and vms_environments alongside it
try   { /* vx_dominant_right_fixed_point / ic_arnoldi_dominant + vx_check_perron_nondegenerate */ }
catch (ITError const&)
    {
    if (C.empty()) throw;
    auto [rb,lb] = vx_bond_fixed_points(C,D);   // reached ONLY on a guard trip
    r_AL = rb; l_AR = lb;
    }
```

so the bond candidate is a failure fallback: it is used when the degeneracy
guard *raises*, and not when the guard passes but the eigenvector it passed is
nonetheless an arbitrary element of a near-degenerate subspace. That is the case
this issue is about -- no exception is thrown, and the wrong fixed point is
accepted silently.

`pyitensor/vumps_ms.py::_cell_fixed_points` and
`pyitensor/vumps.py::_transfer_fixed_points` are **bond-candidate-first**: the
candidate `C C^dag` / `conj(C^dag C)` is *preferred* whenever it reproduces
itself under the transfer map to a residual <= 1e-6 -- in mixed canonical gauge
that is an exact algebraic identity, so it is a yes/no test rather than a tuned
threshold, and at convergence no eigensolve runs at all -- with the eigensolver
used only when the residual says the gauge relation does not hold, and the two
cross-checked against each other by residual when both are available.

## What a fix would involve

Port that ordering, not another guard. Concretely: give
`vumps_build_environments` and `vms_environments` a residual test for the
`vx_bond_fixed_points` candidate, take it when it passes, and fall through to
the eigensolver only when it does not -- then cross-check. Note what is
explicitly the WRONG shape, for the same reason it was wrong on the Python side
and is already recorded in `vx_bond_fixed_points`' own comment and in CLAUDE.md:
a threshold on `C`'s weight spectrum to tell redundancy from a genuine cat
state. That ratio is a moving number mid-convergence (2.7e-9, 1.2e-4 and 1.6e-2
on three cells of the same model) with no defensible cutoff. The residual test
is defensible because it is testing an identity, not a magnitude.

The Python side also pinned `eigs`' start vector (`v0`) while it was there,
which is what made its sequential solver reproducible run to run. The C++
eigensolver's own start is worth checking at the same time; the rate table above
is only meaningful because these runs are *not* reproducible.

## What is and is not guarded in the test suite

**Historical, as of while this was open.**
`tests/test_vumps_redundant_bond_dimension.py::test_sequential_solver_tolerates_redundant_bond_dimension`
applied a narrow, non-strict `xfail` for exactly `backend == 3 and (n_uc, reach)
== (1, 2) and D > 2`, which held the sequential half of this issue; a
deliberately narrow marker was chosen over a loosened tolerance so the other
eleven `(n_uc, reach, D)` combinations kept catching the defect they exist for.
That marker is **gone** -- all three of its cases xpass after the fix.

The GROUPED half was not guarded at all: `test_grouped_solver_tolerates_redundant_bond_dimension`
on `itensor_version=3` at D=6 and D=8 failed whenever an excursion happened to
land past its `abs=1e-9`, a small but real fraction of runs (2 of 80 and
3 of 80 in the measurements above, 1 of 80 and 1 of 80 in a verifier's). Both
are 0 of 80 now.

Two statements in that test file were stale in consequence, and both are
corrected there now: the grouped test's own "0 of 10 at each of D=2,4,6,8 on
both backends" (measured before the rebuild, and untrue of `itensor_version=3`
while this was open), and the module docstring's claim that "All four
environment builders ... now prefer the fixed points the state itself names",
which was precisely what this file measured the C++ half NOT doing. That
sentence is now true; it says so, and says since when.

What is newly guarded is the thing a single-run test could not catch. This
defect appeared in 6 of 30 runs on the sequential `(1,2)` cell at D=4 and 13 of
80 on the grouped `(1,1)` cell at D=6, so the two one-run tests above saw it
~20% and ~16% of the time.
`test_variational_bound_holds_over_repeated_runs` asserts the one-sided bound
over 8 runs of each, which would have caught it ~83% and ~74% of the time.

Two things about that test are deliberate and easy to get wrong if it is ever
rewritten. Its two cells run at *different* bond dimensions, because the two
halves failed at different ones -- the grouped cell is 0 of 80 at D=4 in the
table above and only starts excursing at D=6, so a grouped D=4 row would pin
nothing. And it is `itensor_version=3`-only: the bound is backend-independent,
but the `"python"` grouped D=6 cell carries the separate ~1-in-100
`LinAlgError` noted above, and eight runs of it per suite would be a 5-8% flake
in a test whose purpose is the opposite. The Python side's own guard is
`tests/test_audit_2026_09_pyitensor-infinite.py::test_sequential_vumps_never_returns_below_the_variational_minimum`.

## Where the code is

- `src/dmrgpy/mpscpp3/chain_session.h::vx_choose_fixed_point` -- the shared
  selection both builders now go through (the fix)
- `src/dmrgpy/mpscpp3/chain_session.h::vumps_build_environments` (grouped) and
  `::vms_environments` (sequential) -- which used to hold a
  `catch (ITError const&)` block each
- `src/dmrgpy/mpscpp3/chain_session.h::vx_bond_fixed_points` -- the primitive
  that was already there, and `::vx_fixed_point_residual`, which was not
- `src/dmrgpy/pyitensor/vumps.py::_transfer_fixed_points` and
  `src/dmrgpy/pyitensor/vumps_ms.py::_cell_fixed_points` -- the fixed reference
- `docs/audit_2026_09_hole_hunt.md` finding #3 -- the Python half, with its own
  reproduction

## Reproduction

```python
# floor.py -- any returned e < -1.825 is below the variational minimum
import sys
import numpy as np
from dmrgpy import infinitechain

FIELD, J = 4.0, 0.7
EXACT_E = -FIELD / 2.0 + J / 4.0          # -1.825


def polarized(n_uc, reach, D, backend):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc,
                                           itensor_version=backend)
    ic.maxm, ic.maxiter, ic.etol = D, 300, 1e-12
    ic.vumps_nrestarts = 2
    h = 0
    for i in range(n_uc):
        h = h - FIELD * ic.SzC[i]
        k = i + reach
        other = (ic.SzC[k] if k < n_uc
                 else ic.get_operator("Sz", k % n_uc, group=k // n_uc))
        h = h + J * ic.SzC[i] * other
    ic.set_hamiltonian(h)
    return ic


n_uc, reach, D, backend, nrun = (int(sys.argv[1]), int(sys.argv[2]),
                                 int(sys.argv[3]), sys.argv[4], int(sys.argv[5]))
backend = 3 if backend == "3" else backend
below = 0
below9 = 0
raised = 0
worst = 0.0
for r in range(nrun):
    try:
        e = polarized(n_uc, reach, D, backend).gs_energy()
    except Exception as exc:
        raised += 1
        print("run %3d RAISED %s" % (r, type(exc).__name__), flush=True)
        continue
    if e - EXACT_E < -1e-12:
        below += 1
        if e - EXACT_E < -1e-9:
            below9 += 1
        worst = min(worst, e - EXACT_E)
        print("run %3d e=%.16f BELOW by %.3e" % (r, e, EXACT_E - e), flush=True)
print("n_uc=%d reach=%d D=%d backend=%s runs=%d below=%d below9=%d raised=%d worst=%.3e"
      % (n_uc, reach, D, backend, nrun, below, below9, raised, worst))
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  PYTHONPATH=/path/to/dmrgpy/src taskset -c 11 python3 floor.py 1 1 6 3 80   # grouped
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  PYTHONPATH=/path/to/dmrgpy/src taskset -c 11 python3 floor.py 1 2 6 3 30   # sequential
```

Observed (grouped, D=6, 80 runs, 24.4 s) -- the run the grouped `3 | 6` table
row above is from:

```
run   3 e=-1.8250000000010664 BELOW by 1.066e-12
run   5 e=-1.8250000007764595 BELOW by 7.765e-10
run  20 e=-1.8250000000014475 BELOW by 1.448e-12
run  22 e=-1.8250000000227495 BELOW by 2.275e-11
run  23 e=-1.8250000000021145 BELOW by 2.115e-12
run  29 e=-1.8250000004137759 BELOW by 4.138e-10
run  46 e=-1.8250000000101481 BELOW by 1.015e-11
run  63 e=-1.8250000000446822 BELOW by 4.468e-11
run  68 e=-1.8250000019222643 BELOW by 1.922e-09
run  69 e=-1.8250000000717770 BELOW by 7.178e-11
run  71 e=-1.8250000232616630 BELOW by 2.326e-08
run  73 e=-1.8250000001444653 BELOW by 1.445e-10
run  77 e=-1.8250000000022939 BELOW by 2.294e-12
n_uc=1 reach=1 D=6 backend=3 runs=80 below=13 below9=2 raised=0 worst=-2.326e-08
```

Observed (sequential, D=6, 30 runs, 5.8 s -- note the 2.85e-03 outlier):

```
run   3 e=-1.8250000000015256 BELOW by 1.526e-12
run   8 e=-1.8250000000033524 BELOW by 3.352e-12
run   9 e=-1.8250000001923357 BELOW by 1.923e-10
run  10 e=-1.8278470533394651 BELOW by 2.847e-03
run  21 e=-1.8250000000215043 BELOW by 2.150e-11
n_uc=1 reach=2 D=6 backend=3 runs=30 below=5 below9=1 raised=0 worst=-2.847e-03
```

Observed (the same cells on `itensor_version="python"`, the fixed reference --
summary lines only; the one non-summary line printed across all three runs is
the `LinAlgError` noted under the tables):

```
run  17 RAISED LinAlgError
n_uc=1 reach=1 D=6 backend=python runs=40 below=0 below9=0 raised=1 worst=0.000e+00
n_uc=1 reach=1 D=8 backend=python runs=40 below=0 below9=0 raised=0 worst=0.000e+00
n_uc=1 reach=2 D=6 backend=python runs=40 below=0 below9=0 raised=0 worst=0.000e+00
```

(The last of those is the slow one -- 494.0 s for 40 runs, against 5.8 s for 30
runs of the same cell on `itensor_version=3`, which is the separate and expected
Python-vs-C++ gap, not part of this issue.)

---

## The fix (2026-09-12)

Ported the Python side's **ordering**, not another guard, exactly as "What a fix
would involve" above asked for.

`Chain::vx_choose_fixed_point` is now the single fixed-point selection both
environment builders go through, and it is bond-candidate-first: measure the
candidate's residual under the transfer map it is supposed to be a fixed point
of (`Chain::vx_fixed_point_residual`, the C++ analogue of
`vumps_ms._fixed_point_residual`), take it when that residual is at or below
`vx_bond_fp_residual_tol_ = 1e-6`, and only otherwise run the eigensolver --
then keep whichever of the two reproduces itself better. An empty `C` keeps the
pure eigensolver route byte-identical to before. The two `catch (ITError const&)`
blocks are gone; the eigensolver's own throw is now caught inside the shared
helper, where it still falls back to the candidate and still rethrows when there
is none.

That tolerance is not a tuning knob, which is the whole reason this shape was
chosen over the C-weight-spectrum threshold the earlier attempt tried: in mixed
canonical gauge the bond candidate being a fixed point is an exact algebraic
identity, so the residual is 0 to machine precision when the gauge relation
holds and O(0.1) when it does not. Measured directly on the polarized cell at
D=6 while converging: 4.4e-03 and 3.7e-03 mid-approach, then 2.9e-12, 3.4e-12,
3.8e-15, 6.5e-15 -- there is no band in between to calibrate against.

### A second, latent defect this turned up

`vx_bond_fixed_points`' LEFT candidate was `C^dag C` where this codebase's
`X[ket, bra]` index ordering needs `conj(C^dag C)` -- the transpose, `C^dag C`
being Hermitian. `pyitensor/vumps_ms.py::_bond_fixed_points` has always carried
the conjugate (`np.conj(C.conj().T @ C)`); the C++ port dropped it.

It was invisible for as long as nothing *measured* the candidate. Adding the
residual is what read it: at convergence on the polarized cell at D=6 the
correct orientation reproduces itself to 4e-15 (3e-15 on the grouped path) and
its transpose to 0.38-0.53, against the eigensolver's own 5e-16 -- so this was
not a near miss, and without fixing it the AR side's candidate would simply
never have been accepted and the fix would have been a half no-op. It survived
because the two models this function had ever been exercised on both have a real
symmetric `C^dag C` (a field-polarized chain's converged `C` is real diagonal,
AKLT's is real), where the two orientations coincide. Only a complex `C` tells
them apart -- and VUMPS's own random complex start produces one on every model,
so this mattered wherever the old `catch` fallback actually fired.

The same measurement also confirmed, as a by-product, that the sequential path's
dense `vms_cell_transfer` route and its matrix-free push-chain agree on what
"left fixed point" means: the dense eigenvector's residual under the push-chain
action is ~5e-16. Those are two independent pieces of code for the same map, and
the cross-check between them had never been made before.

### Measured effect

Exactly the script at the top of this section's "Reproduction", verbatim, on the
`_dmrgcpp*.so` rebuilt after the change, thread-pinned to one core the same way.
The "before" rows for the grouped D=6 and sequential D=6 cells were re-measured
on this machine immediately before the change (13 of 80 and 5 of 30, agreeing
with the table above); the other four "before" rows are that table's own.

| solver | cell | D | runs | below (>1e-12) before | after | worst before |
|---|---|---|---|---|---|---|
| grouped | n_uc=1 reach=1 | 4 | 80 | 0 | **0** | -- |
| grouped | n_uc=1 reach=1 | 6 | 80 | 13 | **0** | 2.33e-08 |
| grouped | n_uc=1 reach=1 | 8 | 80 | 18 | **0** | 6.36e-08 |
| sequential | n_uc=1 reach=2 | 4 | 30 | 6 | **0** | 6.13e-04 |
| sequential | n_uc=1 reach=2 | 6 | 30 | 5 | **0** | 2.85e-03 |
| sequential | n_uc=1 reach=2 | 8 | 30 | 3 | **0** | 6.73e-08 |

330 runs, 0 below the exact minimum by more than 1e-12, and 0 raised. Quote the
threshold with the rate, as this file has throughout: "0 of 330 past 1e-12".

Note what did **not** need changing. The eigensolver's start vector was the
other suspect named above, and it is not one: `ic_arnoldi_dominant` already
starts from the identity deterministically. The irreproducibility this file's
rate tables warn about comes from `vumps_random_init`, i.e. the random start
MPS -- which `pyitensor`'s driver has too, so the two backends are on the same
footing there and pinning `v0` would buy nothing. No parameter was added.

### Tests

- `tests/test_vumps_redundant_bond_dimension.py` -- 40 passed, and the 3
  previously-`xfail`ed cases xpass, so the marker was removed rather than
  loosened. Its new `test_variational_bound_holds_over_repeated_runs` is the
  repeated-run guard described under "What is and is not guarded".
- `tests/test_infinite_chain.py`, `test_infinite_long_range.py`,
  `test_lanczos_residual_criterion.py`, `test_vumps_subspace_expansion.py`,
  `test_idmrg_correlator_v3.py`, `test_audit_2026_09_pyitensor-infinite.py` --
  180 passed, 2 skipped. The AKLT rows matter most: a bond-dimension-2 exact
  state is genuinely entangled, so a wrongly-selected element of the degenerate
  subspace would carry its own energy rather than the ground state's, which is
  what says the candidate is the RIGHT fixed point and not merely a harmless
  one.
