# Known issue: `itensor_version=3` VUMPS returns energies BELOW the exact variational minimum

**Status**: NOT fixed. Guarded only on the sequential half, by a narrow
non-strict `xfail`; the grouped half is unguarded and intermittently fails the
suite (see "What is and is not guarded" below). Affects `itensor_version=3` only, on
`gs_method="vumps"` (the default for `Infinite_Many_Body_Chain`), on both the
grouped and the sequential solver. `itensor_version="python"` is unaffected --
this is the C++ counterpart of the 2026-09 audit's finding #3, which was found
and fixed only on the Python side.

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

`tests/test_vumps_redundant_bond_dimension.py::test_sequential_solver_tolerates_redundant_bond_dimension`
applies a narrow, non-strict `xfail` for exactly `backend == 3 and (n_uc, reach)
== (1, 2) and D > 2`, which holds the sequential half of this issue. That marker
should stay until this file says FIXED, and a deliberately narrow marker was
chosen over a loosened tolerance so the other eleven `(n_uc, reach, D)`
combinations keep catching the defect they exist for.

The GROUPED half is **not** guarded: `test_grouped_solver_tolerates_redundant_bond_dimension`
on `itensor_version=3` at D=6 and D=8 will fail whenever an excursion happens to
land past its `abs=1e-9`, which is a small but real fraction of runs (2 of 80 and
3 of 80 in the measurements above, 1 of 80 and 1 of 80 in a verifier's).

Two statements in that test file are stale in consequence, both outside the
scope of this document and both worth fixing there:

* its grouped test's own docstring (`test_grouped_solver_tolerates_redundant_
  bond_dimension`, one occurrence, at "0 of 10 at each of D=2,4,6,8 on both
  backends") -- measured before the rebuild, and no longer true of
  `itensor_version=3`, which is 2 of 80 at D=6 and 3 of 80 at D=8 past the
  `abs=1e-9` this test asserts at. The module docstring does not repeat that
  sentence; its own "It is 0 of 20 now" is about the `"python"` grouped D=4
  case and still holds.
* the module docstring's claim that "All four environment builders (`vumps.py`
  and `vumps_ms.py` on the Python side, `Chain::vx_*`'s grouped and sequential
  halves on the C++ one) now prefer the fixed points the state itself names" --
  which is the more serious of the two, because it is precisely what this file
  measures the C++ half NOT doing: on `itensor_version=3` the bond candidate is
  reached only from a `catch (ITError const&)`, i.e. it is a failure fallback,
  not a preference. Only the two Python builders prefer it.

## Where the code is

- `src/dmrgpy/mpscpp3/chain_session.h::vumps_build_environments` (grouped) and
  `::vms_environments` (sequential) -- the two `catch (ITError const&)` blocks
- `src/dmrgpy/mpscpp3/chain_session.h::vx_bond_fixed_points` -- the primitive
  that is already there
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
