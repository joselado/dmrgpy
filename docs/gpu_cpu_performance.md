# CPU vs GPU in `pyitensor`: what is measured, and when the device wins

Reference numbers for the pure-Python engine's array backend
(`pyitensor/backend.py`, `itensor_version="python"`). Everything here was
measured, not estimated; the design reasoning behind the port is in
`docs/pyitensor_gpu_port_plan.md`.

**Hardware**: NVIDIA H200 (141 GB) versus one Intel Xeon Gold 6248 core,
unless a row says "8 cores". Everything in this engine is complex128, so
these are FP64/ZGEMM numbers and they assume a data-centre GPU with real
double-precision throughput (V100/A100/H100/H200 class). On a consumer
card with 1/32-rate FP64 none of it transfers.

## The one-sentence version

The device wins on **bond dimension**, not on chain length or on which
calculation you run: below chi ~ 120 the CPU is faster, above it the GPU
pulls away without bound, and this is equally true for ground states and
for dynamical correlators.

Two exceptions, and the reason they are exceptions is worth reading before
you generalize from that. The four-point correlator wins on the device at
chi = 20, because `ctmode="batched"` gets its arithmetic-per-dispatch from
a *tuple batch* rather than from bond dimension. And the complex-time
correlator (`submode="TDZ"`) crosses over between chi = 30 and 60, because
it is TDVP-heavy and a two-site matvec at chi = 60 on 30 sites already
carries enough arithmetic per dispatch. The rule is really "the device
needs enough work per array operation", and chi is only the usual way to
get it -- a batch axis, or simply a bigger local tensor, buys the same
thing.

A caveat that applies to both, and to every number below: none of it holds
without `set_pad_bonds` + `set_jit`. Measured on TDZ, those two knobs are
worth 10.1x on their own -- more than every other optimization in this file
combined -- and eager on the device is *slower* than one CPU core at a
bond dimension where the padded, jitted configuration is 3.4x faster.

## End-to-end, warm (steady state, compile cost already paid)

KPM dynamical correlator, 30-site Heisenberg chain:

| kpmmaxm | CPU | GPU | speedup |
|---|---|---|---|
| 40 | 49.2 s | 365.0 s | 0.13x |
| 80 | 239.8 s | 504.3 s | 0.48x |
| 160 | 1775.5 s | 541.9 s | **3.28x** |
| 240 | 6325.3 s | 586.4 s | **10.79x** |

Ground state, 30-site 3-leg ladder (a model that actually needs large
chi -- see the warning below):

| maxm | CPU | GPU | speedup |
|---|---|---|---|
| 60 | 10.2 s | 25.1 s | 0.41x |
| 120 | 31.7 s | 23.0 s | 1.38x |
| 240 | 112.5 s | 23.9 s | **4.70x** |
| 480 | 465.1 s | 22.9 s | **20.27x** |

Read the *shape* of those columns, not just the ratios: GPU time is nearly
flat (365 -> 586 s; 25.1 -> 22.9 s) while CPU time explodes (49 -> 6325 s;
10 -> 465 s). The device is still dispatch-bound at the top of both
sweeps, so **both speedups are lower bounds** -- they keep growing with
chi.

Accuracy at every point above: ground-state energies agree with the host
run to <= 2.2e-11, KPM sum rules to <= 7.9e-07, and every spectrum
satisfies the exact zeroth-moment sum rule (int S_zz domega = 1/4 for a
spin-1/2 site) on both backends.

## Primitives, so you can predict your own case

Per-call cost, complex128, GPU speedup versus 1 core [versus 8 cores]:

| operation | chi=64 | chi=256 | chi=512 | chi=1024 |
|---|---|---|---|---|
| two-site matvec (contraction chain) | 5.5x [5.0x] | 246x [127x] | 480x [169x] | 688x [166x] |
| `eigh` (Gram route in `svd.py`) | 1.2x [1.0x] | 13x [4.2x] | 40x [10x] | 121x [27x] |
| `svd` (exact fallback) | **0.6x [0.6x]** | 2.7x [1.0x] | 6.9x [1.8x] | 24x [3.6x] |
| `qr` | 1.2x [1.0x] | 16x [7.7x] | 47x [18x] | 145x [40x] |

Two things follow. **Contractions are where the device is spectacular**,
which is why `svd.py`'s existing preference for the Gram+`eigh` route
matters here rather than being a mere CPU optimization: `svd` is the one
primitive the GPU is *worse* at below chi ~ 256. And **matrix size, not
call count, is what buys anything** -- see the dispatch floor next.

## The two costs that decide everything

**1. Per-call dispatch floor.** Every eager operation costs ~0.35 ms on
the device regardless of size (measured: torch ~0.09 ms, a `jax.jit`-ed
kernel ~0.07 ms). A calculation issuing thousands of small operations
cannot win, whatever the hardware. This is the whole reason for the
crossover at chi ~ 120-160, and the reason chain length *hurts* the GPU
(more sites = more operations, each paying the floor) while bond dimension
*helps* it (same call count, more arithmetic per call).

**2. Host<->device transfers.** One round trip for a theta-sized array:

| chi | matvec on device | one H2D+D2H round trip | |
|---|---|---|---|
| 64 | 0.34 ms | 0.31 ms | compute-bound |
| 512 | 1.52 ms | 4.31 ms | **transfer-bound (2.8x)** |
| 1024 | 8.31 ms | 21.96 ms | **transfer-bound (2.6x)** |

Above chi=64, shipping an array to the device and back costs *more* than
the arithmetic done on it. Hence the port's central rule: arrays stay
resident, and only scalars (Lanczos alpha/beta, energies, overlaps) and
the O(chi) singular-value vector ever come home. A design that converts
per call is unwinnable at any device speed -- confirmed the hard way,
since `kernels.py`'s older per-call JAX path made GPU runs **5-11x
slower** than plain NumPy.

## Cold versus warm: XLA compiles per shape

Eager JAX compiles a kernel per distinct (operation, shape), and DMRG
mints a new shape whenever a bond dimension changes. Measured on a small
CPU run: **672 compilations costing 18.4 s of a 29.1 s run**, and in a
bigger one 40676 compilations costing 87.6 s. It is a *one-time* tax per
shape, so:

| | cold | warm |
|---|---|---|
| ground state (n=8, maxm=30) | 44.9 s | 1.72 s |
| KPM (same chain) | 64.2 s | 4.93 s |

A script that does one calculation and exits pays the cold price; a sweep
inside one process pays it once. At kpmmaxm=160 on the device that is the
difference between 2.8x and 5.0x against the CPU.

**`backend.set_pad_bonds(K)`** attacks exactly this: freezing every bond
at K collapses the shape zoo (40676 -> 10407 compilations). On the device
that is a good trade, on the host a bad one:

| | cold speedup | warm cost |
|---|---|---|
| GPU, ground state | 1.44-1.84x | 8-30% |
| GPU, KPM | **1.87-3.07x** | 6-31% |
| CPU (either) | 1.64x | **2x slower** |

So: pad for one-shot device runs, don't pad for long warm runs or on the
host. It appends *zero* singular values after truncation has chosen what
to keep, so the state is unchanged (ground-state energies agree to
1.8e-15) -- but it does perturb later truncation *decisions* in a long
recursion (KPM spectra shift by ~3e-3 versus ~3e-4 run-to-run), so it is
exact in representation, not in trajectory. For one family of methods it
was not even the same algorithm: one-site TDVP can populate a padded zero
direction, which made a padded `TDVP_GSE` run one-site TDVP on the
manifold of bond dimension K instead of its Krylov expansion. That route is
exempt from padding now, see "`set_pad_bonds` used to change one-site TDVP
(fixed)" below.

Related: run KPM with `kpmmaxm == maxm`. Otherwise the ground-state solve
and the moment recursion have two separate shape families and every kernel
is compiled twice.

## Lowering the dispatch floor: `set_jit`

Padding stops the engine *minting* new shapes; jitting stops it
*dispatching* so many kernels. The floor is ~0.35 ms per eager operation
against ~0.07 ms for the same kernel under `jax.jit` (H200, Phase 0), so
what matters is how many operations one Lanczos iteration issues.
`backend.jit()` registers the four composites that account for nearly all
of them:

| composite | eager dispatches | where |
|---|---|---|
| planned matvec chain | 3-4 per step, 8-12 per call | `kernels.py::_matvec_chain` |
| contraction: transpose+reshape+matmul | 5 | `tensor.py::_contract_matmul` |
| Gram matrix + `eigh` + descending sort | 5 | `svd.py::_gram_spectrum` |
| Lanczos recurrence / reorthogonalization | 3 / 2 per basis vector | `dmrg.py` |

`backend.set_jit()` controls them: `"auto"` (the default) compiles exactly
when `set_pad_bonds` is set, `True`/`False` force it. The two knobs are
tied together on purpose -- `jax.jit` traces once per input *shape*, so
without padding the compile count rises instead of falling, and the
measurement below shows it doing exactly that.

n=6 Heisenberg, maxm=16, **JAX on CPU**, one pinned P-core
(`benchmarks/gpu/jit_speedup.py`); ground state (4 sweeps) and a TDVP
quench (8 steps), every row returning the identical value:

| case | GS cold | GS warm | kernels | TDVP cold | TDVP warm | kernels |
|---|---|---|---|---|---|---|
| eager, unpadded | 34.8 s | 0.66 s | - | 59.9 s | 2.10 s | - |
| padded, eager | 22.2 s | 1.04 s | - | 31.4 s | 3.54 s | - |
| **padded + jit** | **7.7 s** | 0.80 s | 86 | **11.1 s** | 2.41 s | 93 |
| jit, unpadded | 12.6 s | 0.48 s | 266 | 31.3 s | 1.32 s | 432 |

That is a host measurement; the floor the knob attacks is a device
property, so it is a lower bound. **On an H200** (n=20 chain, jobs
19951209/19951210, ratios against the eager/unpadded row):

| calculation | eager cold/warm | padded+jit | jit alone |
|---|---|---|---|
| ground state, maxm=32 | 127.5 s / 7.29 s | 18.6 / 6.20 (**6.8x** cold, 1.17x warm) | 48.2 / 5.65 (2.6x, 1.29x) |
| ground state, maxm=64 | 167.8 s / 7.15 s | 13.8 / 6.58 (**12.1x**, 1.09x) | 62.5 / 5.56 (2.7x, 1.28x) |
| TDVP quench, maxm=32 | 130.0 s / 14.82 s | 20.3 / 10.93 (**6.4x**, 1.36x) | 52.6 / 9.15 (2.5x, **1.62x**) |
| TDVP quench, maxm=64 | 129.0 s / 14.37 s | 17.7 / 11.72 (**7.3x**, 1.23x) | 54.1 / 9.09 (2.4x, **1.58x**) |

Three things follow, and the third matters most:

* **Cold: 6.4-12.1x on the device**, against 4.5-5.4x on the host -- the
  effect is about twice as large where it was predicted to be. Both knobs
  are needed: padding alone gives 2.4-3.8x, jit alone 2.4-2.7x, and jit
  alone traces 5-9x more kernels (897-1092 against 114-192) because
  nothing is holding the shapes still.
* **Warm: 1.1-1.6x, and the ranking flips.** Padding does real arithmetic
  on blocks known to be zero, so once every kernel is compiled it is jit
  *without* padding that wins (1.26-1.29x on a ground state, 1.58-1.62x
  on TDVP). So: `set_pad_bonds` + `set_jit` for a script that runs once
  and exits; `set_jit(True)` alone for a long sweep in one process.
* **It does not move the crossover.** The same sizes on one CPU core run
  in 0.79 s (ground state) and 1.24 s (TDVP) warm -- still ~7x faster
  than the best device configuration here. Lowering the floor makes the
  device far cheaper to *start*; it does not turn a small problem into a
  large one, and a 20-site chain at these bond dimensions is
  dispatch-bound whatever is done to the dispatch. The crossover is still
  a chi ~ 120-160 story.

`backend.compilations()` reports the traced-kernel count during a run.

One caveat on the ground-state rows, which is this file's own trap
(below) biting again: a uniform 20-site Heisenberg chain saturates well
before maxm=128, so its maxm=64 and maxm=128 runs are the same
calculation -- identical energies to 13 digits and flat times on both
devices (166-168 s, 2.22-2.23 s). Read them as dispatch-floor
measurements at fixed work, which is what they are, not as a chi sweep.

## Real-time evolution runs on the device too

TDVP (`timedependent.evolve_and_measure`, `submode="TD"` correlators) was
ported on 2026-08-26: the Krylov propagator's basis now stays resident,
with only alpha/beta per iteration coming home. Before that it kept the
basis in a NumPy buffer, i.e. one round trip per Krylov iteration per bond
per time step -- silently, since a transfer returns the same numbers.

This is the calculation the device suits best, for a physical reason: a
quench grows entanglement roughly linearly in time, so chi climbs into the
paying range by construction, while a 1D ground state's is capped by an
area law (see the trap below). Cross-backend agreement on a full
trajectory: 1.5e-14. The chi-swept GPU table for it has not been run yet;
`examples/backend_comparison/tdvp_cpu_VS_gpu` is the script that produces
it.

## Complex-time evolution (`submode="TDZ"`): removing synchronizations

`submode="TDZ"` (`tdz.py`, complex-time evolution + perturbative real-axis
reconstruction, arXiv:2311.10909) sits on the TDVP path above and adds
n_max+1 full-chain overlaps per time step. Profiled on the host it is ~90%
TDVP and ~8% those overlaps, but on a device the split that matters is a
different one: nothing in it is a large GEMM, so its cost is made of
per-call dispatch and, worse, of *host synchronizations*. JAX dispatch is
asynchronous -- the host runs ahead enqueueing kernels while the device
works, which is the only thing hiding the per-call floor -- and every
`float(bk.to_host(...))` drains that queue.

Three changes (2026-08-27) attack exactly that. None of them changes any
arithmetic, so all of them are verified by the spectrum being unchanged:

* **The Krylov exponentiator no longer synchronizes per iteration**
  (`tdvp.py`, `_lanczos_expm_device`). alpha and beta stay 0-d device
  arrays -- `w - alpha*q` and `w/beta` never needed a host value -- and a
  whole block of them comes home in one transfer. Since the stopping test
  cannot then be evaluated every iteration, the recursion *speculates*
  past its own stopping point and rolls back to the exact same Krylov
  dimension the per-iteration loop would have chosen, so the returned
  vector is identical rather than merely as accurate. At n=30 this takes a
  two-site TDVP step from ~1000 synchronizations to ~60. The next
  checkpoint is placed at the previous call's converged k, which is an
  excellent predictor here (Krylov convergence is set by the local
  effective Hamiltonian's spectral width, which barely moves between
  bonds), so the speculation usually wastes nothing.
* **The phi^(n) overlaps are batched** (`mpsalgebra.BatchedBras`). The
  bras are fixed for the whole run -- H^n(B|GS>), built once -- so they are
  stacked, conjugated and zero-padded to a common per-bond width once, and
  each step's n_max+1 overlaps become one sweep of two batched GEMMs per
  site instead of n_max+1 separate sweeps. Same arithmetic, a fifth of the
  dispatches. Padding is exact for the same reason `set_pad_bonds` is: a
  padded column contracts to zero.
* **`svd()` stopped round-tripping its S tensor.** The diagonal
  singular-value tensor was assembled with `np.diag()` from the host copy
  of the spectrum, so every factorization pushed a `keep x keep` matrix
  back across the bus -- at every bond, of every half-sweep, of every step,
  which makes it the most frequent O(chi^2) transfer in the engine. It is
  built from the device-resident spectrum now. The same call also fetched
  the O(chi) spectrum twice (once inside `_svd_truncated` for the
  truncation rule, once again for the `Spectrum`); it is fetched once and
  reused, halving the remaining per-SVD stalls.

The one synchronization deliberately left in place is the truncation rule
itself: `keep` is a cumulative sum with data-dependent branching over the
spectrum, worthless on a device and not comparable against a Python float
there, so O(chi) numbers come home per factorization by design.

### Does TDZ actually reach a bond dimension a device cares about?

Worth settling before reading any timing, because the answer was not what
this section originally predicted. The paper's headline is that the
complex-time contour keeps chi small -- it reports chi ~ 20-30 against
500-700 for real-time evolution -- and this port's crossover is chi ~
120-160, so the expectation written here first was that TDZ sits on the
wrong side of the line by construction.

Measured instead (n=30 Heisenberg, `alpha0=0.1`, `dt=0.1`, one core, the
largest bond dimension anywhere in the MPS):

| t | 0.1 | 1.1 | 2.1 | 3.1 | 4.1 | 5.1 | 6.1 | 6.6 |
|---|---|---|---|---|---|---|---|---|
| chi | 72 | 83 | 98 | 125 | 156 | 191 | 230 | **240 (cap)** |

It starts at 72 -- the state is `Sz|GS>`, so it inherits the ground
state's own bond dimension rather than starting at 1 -- crosses the
120-160 crossover at t ~ 3-4, and saturates a `maxm=240` cap at t = 6.6,
a third of the way through a `nt=100` run. So on this model, at this
alpha0, TDZ is *not* a small-chi calculation: the contour slows
entanglement growth, it does not prevent it.

Two consequences. Every bond dimension in the sweep below is genuinely
binding (at `maxm=30` and `60` the state is truncated from the very first
step, since it begins at 72), so no row is measuring padding tax on bond
dimension the state never reaches. And the paper's chi ~ 20-30 should be
read as a statement about its own impurity model and contour angle, not a
property of the method that transfers to every Hamiltonian -- a larger
`alpha0` damps harder and would push these numbers down.

### Measured, n=30

H200 against one Skylake Xeon core (`batch-skl`, not the Cascade Lake the
other tables in this file use -- see `tdz_cpu.sbatch`), `nt=100`,
`dt=0.1`, `alpha0=0.1`, `n_max=4`, one full TDZ correlator per cell.
**Warm** seconds (second run of the same size in the same process),
device configuration padded + jitted. The ground state is inside the timed
region and is ~5 s on the host, 13-47 s on the device, so the
correlator-only device advantage is slightly larger than these totals
show.

| maxm | CPU base | CPU both | GPU base | GPU +krylov | GPU +bras | GPU both | GPU both vs CPU |
|---|---|---|---|---|---|---|---|
| 30 | 98.6 | 99.9 | 135.5 | 128.6 | 130.5 | 124.7 | **0.80x** |
| 60 | 459.4 | 458.1 | 144.5 | 137.9 | 141.8 | 135.8 | **3.37x** |
| 120 | 2214.9 | 2105.4 | 200.6 | 157.5 | 197.6 | 151.5 | **13.9x** |
| 240 | not run | not run | 875.0 | 833.5 | 885.9 | 832.3 | -- |

Every cell agrees with its `base` to <= 4.2e-16, which is the claim the
whole exercise rests on: these are scheduling changes and they change no
number.

The crossover for TDZ sits between maxm 30 and 60 -- *below* the chi ~
120-160 this file quotes elsewhere, because TDZ is TDVP-heavy and a
two-site TDVP matvec at chi=60 on 30 sites already carries enough
arithmetic per dispatch. Read the CPU column's shape rather than only the
ratios: it grows ~4.8x per doubling of maxm while the device column grows
1.1x, 1.4x, 5.5x, so the ratio is still climbing at the top of the sweep
and 13.9x is a lower bound. The maxm=240 CPU row was not measured -- it
extrapolates to ~2.9 h per run and the job was stopped before reaching it.

### What the attribution actually says, including where it disappoints

Warm, device, `base` / `both`: **1.10x, 1.07x, 1.32x, 1.05x** at maxm
30/60/120/240. Two honest readings of that:

* **The Krylov synchronization removal is the change that pays**; the
  batched overlaps are ~neutral (1.00-1.02x on their own, occasionally
  slightly negative). That was predictable in advance and was not
  predicted: the overlaps are ~8% of the work, so by Amdahl a 5x
  improvement there is capped at 1.07x *however well it is done*. The
  batching was over-invested relative to its ceiling. It is kept because
  it costs nothing, it is exact, and it is the piece that would matter if
  the batch axis below is ever built -- but it is not what made this
  faster.
* **The spread 1.05-1.32x is not a trend.** `--reps 2` gives exactly one
  warm sample per cell, so few-percent differences between configurations
  are not resolvable and the 1.32x at maxm=120 may be an outlier rather
  than a peak. Treat the honest summary as "~1.1x, with one size showing
  1.3x", and re-run with more reps before quoting a curve.

### The knob that dominates all of this: pad + jit

Second pass, same job, eager (no padding, hence no jit), maxm=60:

| config | eager | padded + jitted |
|---|---|---|
| base | 1456.7 s | 144.5 s |
| both | 1694.4 s | 135.8 s |

**10.1x**, an order of magnitude more than anything else in this section --
and eager on the device is 3.2x *slower* than one CPU core at the same
size, against 3.4x faster with the knobs on. For this calculation
`set_pad_bonds` + `set_jit` is not a tuning option, it is the difference
between the device being worth using and not.

The second row is the finding that inverts a prediction. The
synchronization work was expected to help *most* in eager mode, where the
dispatch queue it unblocks is longest. It does the opposite: eager,
`both` is **0.86x**, i.e. 14% slower than `base`. The mechanism is the
speculation -- running past the Krylov stopping point costs a few extra
matvecs per call, and with no jit each of those is a full eager dispatch
at the ~0.35 ms floor, which outweighs the synchronizations saved. The
optimizations and the jit knob are therefore not independent: they
compose, and the speculative path assumes the fused kernels are there.
Anyone running this engine eager on a device should call
`tdvp.set_krylov_defer_sync(False)`.

The lever that would change that verdict is the one the four-point
correlator found and TDZ has not used yet: a genuine batch axis. Computing
C_ij(t) for many j from the *same* evolution -- the whole dynamical
structure factor rather than one operator pair -- would put ~N bras into
the batch that currently holds n_max+1 of them, which is the regime where
the next section's result says a device wins at chi=20. That is an API
change (dmrgpy has no multi-pair dynamical-correlator entry point today),
so it is named here as the follow-up rather than built.

## The four-point correlator: batching beats bond dimension

The one-sentence version above has one exception, and it is the most useful
result in this file: the four-point tensor
`<Cdag_i C_j Cdag_k C_l>` wins on the *device* at `chi = 20`, an order of
magnitude below the chi ~ 120-160 crossover everything else obeys. Nothing
about the hardware changed -- what changed is where the arithmetic-per-
dispatch comes from. Every other calculation here gets it from bond
dimension, so it needs a large chi to clear the floor. `ctmode="batched"`
(`pyitensor/fourpoint.py`) gets it from the *tuple batch* instead: the
`O(n^4)` tuples collapse onto a trie of a few dozen environment arrays,
each `(B, chi, chi)` with `B` up to `C(n,3)`, so a single GEMM carries tens
of thousands of tuples at any chi at all.

n=30 spinless fermionic chain, H200 against one Xeon core, full tensor:

| maxm | chi | ITensor v3 (C++) | `"sweep"` (host) | `"batched"` (host) | `"batched"` (device, cold) | `"batched"` (device, warm) |
|---|---|---|---|---|---|---|
| 20 | 20 | 24.7 s | 44.0 s | 1.89 s | 20.5 s | **0.97 s** |
| 40 | 40 | 40.0 s | -- | 9.44 s | 21.5 s | **0.95 s** |
| 80 | 69 | -- | -- | 22.8 s | 24.8 s | **0.98 s** |

Read the warm column vertically. It does not move: 0.97, 0.95, 0.98 s while
the host column grows 12x over the same range. The device is entirely
dispatch-bound here -- the GEMMs themselves are free on an H200 at these
sizes -- so warm device/host runs 1.95x, 9.90x, 23.34x purely because the
*host* gets slower, and every one of those figures is a lower bound. Against
the compiled C++ backend the same warm number is 25x at maxm=20 and 42x at
maxm=40. Agreement with the host result: 7.8e-16 to 1.3e-15.

At n=100 the same kernel enters a different regime, and both of the
sentences above stop being true in an interesting way. 10^8 tuples, 23.5M
distinct-index leaf values, a 1.6 GB output tensor:

| maxm | chi | `"batched"` (host) | `"batched"` (device, cold) | `"batched"` (device, warm) |
|---|---|---|---|---|
| 20 | 20 | 481.6 s | 103.7 s (4.65x) | **19.9 s (24.19x)** |
| 40 | 40 | 2171.8 s | 328.6 s (6.61x) | **79.1 s (27.44x)** |

Agreement with the host: 1.8e-15 and 2.8e-15. Two things changed:

* **The device wins cold too.** At n=30 a single one-shot run was 0.09-0.92x
  -- the compile cost swamped everything. Here the host needs 8 minutes at
  maxm=20, so the same compilation is amortized and even one tensor computed
  once is 4.65x ahead. The operating rule below is therefore not only
  "several tensors in one process" but also "one tensor, if n is large
  enough".
* **The device is no longer flat.** 19.9 s against 0.97 s at n=30 is 20.5x
  for 143x the leaf values, so the H200 is now doing real arithmetic between
  dispatches rather than idling. Part of that is the sweep being *blocked* at
  this size: one block would hold 6*C(100,3) = 941094 level-3 environments,
  6.0 GB at chi=20, so `fourpoint.py` splits the sweep on the first occupied
  site (exact -- environments with different first sites never merge) and
  picks a block width from a byte budget. Narrower blocks cost dispatches,
  which is why the width is a budget rather than one-block-per-site: at
  n=100, chi=20 it lands on 11 first-sites per block.

Extrapolating the compiled v3 backend by its measured n^4 scaling puts a
single n=100, maxm=20 tensor at ~51 min, i.e.\ ~6x slower than batched on
one core and ~150x slower than batched on the H200 -- but that is an
extrapolation from the n=30 row, not a measurement, and should be quoted as
one.

The cold column at n=30 is the cost of that flatness: ~20-25 s, near-constant,
and it is XLA compiling a kernel per array shape. This kernel mints a fresh
leading dimension `B` for every trie node at every site -- thousands of
distinct shapes -- and `set_pad_bonds` cannot help, because what varies is
the environment batch, not the bond. So the operating rule is the mirror of
the one in the cold-versus-warm section below: **the device pays off when
several tensors are computed in one process, or when chi is large enough
that ~22 s of compilation is small next to the host time** (at maxm=80 a
single cold run is already break-even at 0.92x). One tensor, one script,
small chi: stay on the host, where `"batched"` is still 13x faster than the
compiled v3 backend.

Two practical notes:

* the ground state does **not** have to be on the device. `_arrays_lpr`
  converts the MPS once, `O(n chi^2 d)`, so `backend.set_backend("jax")`
  immediately before the correlator call is enough -- which avoids paying
  DMRG's own below-crossover device tax entirely. Every row above did
  exactly that.
* this kernel is deliberately *not* jitted; see
  `docs/documentation.md`'s section on `pyitensor/fourpoint.py` for why the
  batch axis and `jax.jit` are mutually exclusive here.

## Consumer GPUs: measured on a GTX 1060 laptop

Every other table in this file was measured on an H200 and says so, with a
warning attached that on a consumer card with 1/32-rate FP64 none of it
transfers. That warning was an inference, not a measurement, and this
section is the measurement: **NVIDIA GeForce GTX 1060 Mobile (6 GB,
Pascal sm_61)** against the **Intel i7-7700HQ** (4 cores / 8 threads,
AVX2) in the same laptop, jax 0.11.1 on CUDA 13.0, threads pinned
(`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1`), `taskset` on the CPU runs.

The inference was right about the conclusion and wrong about the reason,
which matters for deciding what to do about it.

### The port is correct here; that was never the question

All 17 tests in `tests/test_pyitensor_gpu_backend.py` and
`tests/test_metts_gpu_backend.py` pass on the CUDA device, and end to end
the device reproduces the host *exactly*: ground-state energies agree to
all 10 printed digits at every bond dimension in the sweep below
(-13.8771534902, -13.8811961948, -13.8813580417, -13.8813610754,
-13.8813610980), and the KPM sum rule lands on 0.250000 (error 3.4e-07)
on both. Nothing about a consumer card makes the port less faithful. What
it changes is whether running it is worth doing.

### The GPU is not dispatch-bound here -- it is FP64-bound

This is the finding that separates this card from an H200, and it inverts
the advice in the rest of this file. Achieved GFLOP/s on the two-site
matvec (`benchmarks/gpu/gpu_microbench.py`, complex128, median of 5):

| chi | CPU 1 core | CPU 4 threads | GPU | GPU, complex64 |
|---|---|---|---|---|
| 128 | 33.1 | 60.8 | 92.8 | 896.8 |
| 256 | 41.0 | 90.4 | 122.6 | 1757.5 |
| 512 | 44.2 | 115.8 | 120.2 | 2513.9 |
| 1024 | 45.9 | 127.8 | 99.0 | 2990.4 |

The GTX 1060's FP64 peak is ~138 GFLOP/s (1/32 of its ~4.4 TFLOP/s FP32);
the i7-7700HQ's is ~218 GFLOP/s. The GPU column reaches **123 GFLOP/s,
89% of its own FP64 ceiling** -- it is flat out, with no headroom a
dispatch or scheduling fix could recover, while the CPU column reaches
59% of its ceiling. The `complex64` column is the proof: the same
kernel in single precision is **30.2x faster at chi=1024**, which is
this card's 1/32 FP64:FP32 ratio almost exactly. (That column is a
what-if -- the engine is complex128 throughout, see
`pyitensor/tensor.py` -- produced by `gpu_microbench.py --dtype
complex64`. Read it as a ceiling on what a mixed-precision port could
buy on such a card, not as a dmrgpy timing.)

So on this hardware the H200 story runs backwards. There, the device had
enormous arithmetic throughput and the whole problem was feeding it, so
every lever was about dispatch: pad the bonds, jit the composites, keep
arrays resident. Here the arithmetic throughput *is* the wall, the
per-call speedups stop growing with chi instead of climbing without
bound, and the ratios below are ceilings rather than lower bounds.

Per-call primitive ratios (GPU vs CPU, complex128) make the same point:

| operation | chi=128 | chi=256 | chi=512 | chi=1024 |
|---|---|---|---|---|
| matvec, vs 1 core | 2.80x | 2.99x | 2.72x | 2.16x |
| matvec, vs 4 threads | 1.53x | 1.36x | 1.04x | **0.78x** |
| `eigh`, vs 1 core | 1.01x | 1.73x | 3.73x | 3.85x |
| `svd`, vs 1 core | **0.10x** | **0.09x** | **0.09x** | **0.49x** |
| `qr`, vs 1 core | 0.79x | 1.21x | 1.85x | 1.92x |

The matvec row is the one to read: on an H200 it goes 5.5x -> 688x across
this range, here it *peaks at chi=256 and falls*. Against the whole CPU
it is already losing by chi=1024. `svd` is a rout at every size, an
order of magnitude worse -- which makes `svd.py`'s existing preference for
the Gram+`eigh` route load-bearing on this hardware rather than merely
nice.

### End to end: 3-leg ladder, n=24, warm seconds

`benchmarks/gpu/port_speedup.py --model ladder3 --n 24 --reps 2`. The GPU
column is the better of the two device configurations at each size
(padded+jitted up to maxm=30, eager above it -- see the padding note
below); the energies are identical across all four columns.

| maxm | CPU 1 core | CPU 4 threads | GPU | GPU vs 1 core | GPU vs 4 threads |
|---|---|---|---|---|---|
| 30 | 2.28 | 1.46 | 14.96 | 0.15x | 0.10x |
| 60 | 6.77 | 3.96 | 19.67 | 0.34x | 0.20x |
| 120 | 19.05 | 9.72 | 25.46 | 0.75x | 0.38x |
| 240 | 53.23 | 24.12 | 38.40 | **1.39x** | 0.63x |
| 360 | 93.08 | 41.92 | 76.79 | **1.21x** | 0.55x |

KPM dynamical correlator (n=16, kpmmaxm=40): 61.4 s on one core, 39.2 s on
four, **898 s cold on the device** -- 15x slower than one core, and the
calculation this file's H200 table has winning by 3.3x at kpmmaxm=160.
KPM issues thousands of small operations, so it is the worst case for a
card that is slow per FLOP *and* pays a dispatch floor.

Three things to take from the ground-state table:

* **The crossover against one core is maxm ~ 180**, not the chi ~ 120-160
  quoted elsewhere here -- close, but reached for a different reason, and
  it is a *peak* rather than a threshold: 1.39x at maxm=240 falls to 1.21x
  at 360, because the GPU is already at its FP64 ceiling while the CPU
  still has cache to lose. On an H200 the same column climbs to 20x.
* **Against the whole CPU there is no crossover in this range at all.**
  0.63x is the best cell in the last column. A laptop user does not run
  one core -- so on this machine the honest summary is that the GPU never
  wins a ground state, and the one-core column is there to show where the
  crossover *would* be if it did.
* **The 4-thread column is not free either.** Below maxm ~ 120 it is the
  BLAS-oversubscription case `CLAUDE.md` warns about; measured on a
  contaminated first pass (a leftover device job holding a core) the same
  4-thread sweep read 17.9 s at maxm=30 against 1.46 s clean, i.e. 12x.
  Kill everything else before timing on a 4-core laptop; the numbers above
  are from a re-run on an idle machine.

### `set_pad_bonds` used to pad the Hamiltonian MPO too (fixed)

Found here because 6 GB is small enough to turn a constant factor into a
hard failure. `mpobuilder.to_mpo` compresses the finite-state machine with
the same `position()`/`svd()` every MPS sweep uses, so `set_pad_bonds(K)`
padded the *operator's* bonds along with the state's: measured on a
12-site next-nearest-neighbour spin chain, `set_pad_bonds(60)` took the
Hamiltonian MPO from bond dimension 8 to 60.

That is not a small tax. The MPO bond `w` is a linear factor in the
two-site matvec's dominant O(chi^3 d^2 w) term and in every environment
tensor, and padding cannot buy anything back here, because an operator's
bonds never move -- it is built once and keeps its shape for the whole
run, so there was no shape churn to collapse. On this card it produced a
single 1.77 GiB allocation the allocator could not serve; XLA fell back to
a slower plan and a padded `maxm=60` ground state did not finish in 40
minutes, against 6.8 s on one CPU core. The same padding on an H200 is an
invisible constant factor, which is why it survived the original port.

`mpscontainer._Chain._pad_bonds` (True on `MPS`, False on `MPO`) now
exempts operators, via `backend.pad_bonds_suspended`. MPS bonds are still
frozen at K, which is what padding is for. Measured effect on this card,
same ladder3 ground state, padded configuration:

| maxm | before | after | peak device memory |
|---|---|---|---|
| 30 | 27.52 s | **14.96 s** | 4810 -> 1304 MiB |
| 60 | did not finish (>40 min) | **23.40 s** | OOM -> fits |

Energies are unchanged (identical to 10 digits), which is the claim the
change rests on: it removes arithmetic on blocks that were known zeros.
`tests/test_pad_bonds_mpo_exemption.py` pins it, on the host, where it
needs no GPU to check.

Note what the fixed padding does *not* do here: it wins at maxm=30 and
loses above it (65.9 s padded against 25.5 s eager at maxm=120, 281 s
against 38.4 s at maxm=240). Padding trades real arithmetic for shape
stability, and on a card this slow per FLOP the arithmetic is the
expensive half. **On a consumer GPU, pad only at small bond dimension**;
the opposite of the advice for an H200.

### `set_pad_bonds` used to change one-site TDVP (fixed)

Not specific to a consumer card, and found by the 2026-09-24 audit
(`docs/audit_2026_09_24_hole_hunt.md`, finding 16) rather than by a
benchmark. Padding appends zero singular values, which leaves the state
unchanged at every instant, and for a two-site method that is the whole
story, since the next SVD discards the zero directions again. A one-site
method does not truncate between its steps, and in the first
left-to-right half-sweep `qr_split`, a reduced QR, completes the padded
zero directions into live zero-weight basis vectors that one-site TDVP
then populates. So under `set_pad_bonds(K)` at the recommended K = maxm,
`tevol_method="TDVP_GSE"` added no direction at any bond of any call (27
of 27 bond steps at n=10, K=4, and 33 of 33 at n=12, K=8, against 1 to 4
per bond unpadded), because `gse._gse_bond_step` read the padded dimension
as the bond's rank and found no room left, and its bond growth came from
the QR completion rather than from the Krylov subspace. You can think of
the padded route as one-site TDVP on the manifold of bond dimension K: a
larger ansatz, not a wrong integrator, but not the method that was asked
for either.

The one-site route is now exempt from padding the way the MPO is:
`Chain.global_subspace_expand` and `Chain.tdvp_step(num_center=1)` run
under `backend.pad_bonds_suspended()`, `quench_tdvp_gse` and
`evolve_and_measure_tdvp_gse` strip the padding once at trajectory entry
(`_strip_bond_padding`, a lossless SVD sweep on the evolved state only, so
the caller's `wf` keeps its padding), and `_gse_bond_step` reads the true
rank from the spectrum, which keeps a direct padded caller of `gse.py`
right as well. Measured on an XXZ quench (Delta=0.7, hz=0.1, Neel start,
40 steps of dt=0.05), the largest distance between the padded and the
unpadded `<Sz_0>(t)` over the trajectory:

| chain | K | `tdvp_gse_sweeps` | before | after |
|---|---|---|---|---|
| n=10 | 4 (= maxm) | 0 | 0.4928 | 8.9e-16 |
| n=12 | 8 | 0 | 0.4929 | 1.3e-14 |
| n=10 | 4 | 3 | 4.086e-6 | 1.70e-7 |
| n=12 | 8 | 3 | 1.911e-7 | 4.7e-10 |

A `quench_tdvp_gse` correlator at n=8, K=4 went from 0.34 to 1.2e-15, and
on the Gram SVD route (K=24, above `_GRAM_MIN_DIM`) the padded and
unpadded runs now agree to 2.9e-15 and 1.3e-10 from a Neel start and to
1.2e-12 and 1.1e-12 from an entangled start. Unpadded runs are unchanged,
bit for bit at `tdvp_gse_sweeps` 0 and 3. Note the direction at
`tdvp_gse_sweeps=0`: the padded run used to sit 9.8e-5 from ED only
because it was one-site TDVP on the larger manifold, and it now sits
0.4929 from ED, like the unpadded frozen product state, which is the
correct one-site answer.

Two things are left as they are, deliberately. `Chain.tdvp_step` never
strips, because `submode="TDZ"` carries its wavefunction between calls and
a per-call strip would delete the expansion's zero-weight directions, so a
TDZ run with `TDVP_GSE` at `tdvp_gse_sweeps=0` under padding still starts
from the padded state (0.467 against unpadded in an emulated TDZ loop,
Neel start, n=8, K=4), while at `tdvp_gse_sweeps=3` its first expansion
strips it (to 8e-13). And under JAX with `set_jit("auto")` this route
again retraces once per bond dimension it grows through, the cost padding
was meant to remove, though it never kept frozen shapes anyway.
`tests/test_audit_2026_09_24_pyitensor.py` pins all of it on the host,
where it needs no GPU to check.

### The operating rule for a consumer card

Use the CPU. Concretely: on a 6 GB Pascal laptop the pure-Python engine is
faster on the host for every ground state up to at least maxm=360 and for
every KPM correlator, and the device's one-core crossover at maxm ~ 180
does not survive contact with the other three cores. A card with real
FP64 (V100/A100/H100/H200 class, 1/2-rate) is a different machine and the
rest of this file applies to it.

### Single precision does not rescue it (measured)

The `complex64` column above prices the *arithmetic* at ~30x, so the
obvious follow-up is to run the engine in single precision. It was tried
directly, with a temporary dtype switch in `backend.py`'s three array
constructors (`asarray`/`zeros`/`eye`), and it fails at every level that
matters:

| test | complex128 | complex64 |
|---|---|---|
| n=8 Heisenberg, full bond dimension (no truncation), error vs ED | 1.1e-14 | 3.1e-7 |
| n=16 Heisenberg, maxm=40, E0 | -6.911737 | -6.906817 (error 4.9e-3) |
| ladder3 n=24 maxm=120, host | E0 -13.8813580417 | **crash**: `SVD did not converge` |
| ladder3 n=24 maxm=120, GPU | E0 -13.8813580417, 24.7 s warm | **E0 = -20297.4**, 1628 s warm |

Read the rows in order, because each one fails in a different way:

* **Without truncation the floor is ~3e-7** -- float32 rounding, what
  anyone would expect, and on its own a usable if unimpressive accuracy.
* **With truncation it is four orders of magnitude worse.** `svd.py`'s
  Gram route diagonalizes `M M^dag`, i.e. the *squared* spectrum, so in
  float32 (eps ~ 1.2e-7) singular values below ~sqrt(eps) ~ 3.5e-4 are not
  resolved at all -- which is exactly where a truncation decision lives.
  The Lanczos tolerances (1e-12, and VUMPS's `residual_tol`) are likewise
  below what float32 can represent, so no solve can meet its own
  stopping criterion.
* **At a size worth putting on a device it diverges.** An energy of
  -20297 on a 24-site spin-1/2 ladder is not a poor answer but an
  impossible one -- the Lanczos basis loses orthogonality and the Rayleigh
  quotient runs away -- and **nothing raises**. On the host the same
  calculation at least dies loudly inside LAPACK.
* **It is 66x *slower* on the GPU, not 30x faster.** With the numerics
  broken every local solve runs to its iteration cap and the truncation
  keeps landing on the exact-`svd` fallback, which is the one primitive
  this card is worst at (0.09x against one CPU core, table above). The
  30x lives in the matvec; the engine stops being matvec-dominated the
  moment it stops converging.

So the switch was removed again rather than kept as an experimental knob:
a setting that returns -20297 without an error is a trap, not an option.
What single precision would actually require is genuine *mixed*
precision -- float32 only inside the two-site matvec, with the Lanczos
recurrence and reorthogonalization, the Gram/SVD truncation and every
reported quantity kept in float64 -- and that is a research port with an
unmeasured payoff, not a dtype flag. On this laptop it would also be
competing against a CPU that already wins.

## Device compatibility, calculation by calculation

Everything above is about *speed*. This section is about whether a
calculation runs on the JAX backend at all, whether it gives the NumPy
answer, and whether it really stays on the device -- three separate
questions, and a calculation can pass the first two and fail the third
with no visible symptom. Measured 2026-09-16 on the GTX 1060 laptop,
jax 0.11.1, tiny sizes (the point is coverage, not timing).
`tests/test_pyitensor_gpu_backend.py` already covered ground states,
static/KPM/TDZ correlators and TDVP; this sweep covers the rest.

| calculation | runs | vs NumPy (max abs diff) | on the device? |
|---|---|---|---|
| bond entanglement entropy | yes | 1.1e-13 | resident |
| excited states | yes | 4.9e-15 | resident |
| conserved-sector ground state | yes | 2.0e-14 | resident |
| `submode="SECTOR"` spectral function | yes | 2.2e-08 | resident |
| `submode="CVM"` correlator | yes | 8.7e-14 | resident |
| `submode="TD"` correlator | yes | 5.7e-15 | resident |
| TEBD evolution | yes | 1.2e-14 | resident |
| TDVP-GSE evolution | yes | 2.8e-11 | transfers only while bonds grow |
| four-point tensor, `ctmode="batched"` / `"full"` | yes | 1.6e-15 / 1.9e-15 | resident |
| MPS algebra (`exponential`, MPO application, `vev(npow=2)`) | yes | 7.3e-08 (see below) | resident |
| **iDMRG** (energy, `vev`, `correlator`, `local_excitation_gap`) | **no -> fixed** | <1e-12 energy, <4e-10 observables | **per growth step** |
| VUMPS (energy, `vev`, `correlator`) | yes | 2.5e-16 | **per iteration** |
| non-Hermitian DMRG (`nhdmrg`) | yes | 3.2e-14 | **per Arnoldi matvec** |

"vs NumPy" is meaningful because NumPy is bit-for-bit deterministic from a
fixed seed (re-running it gives a difference of exactly 0), so any nonzero
entry is the device's own roundoff carried through an iterative solve. The
two largest are not defects. The 7.3e-8 is `gs_energy_fluctuation` alone,
<H^2>-<H>^2 of a converged state: two ~6.2 numbers cancelling to ~1e-7 on
*both* backends (9.8e-8 host, 1.7e-7 device), i.e. noise at the precision
floor of the quantity itself; every other MPS-algebra entry agrees to
<=4e-12. The 2.2e-8 sits at the per-sector eigensolver's tolerance.

### iDMRG did not run on a device at all (fixed)

`pyitensor/idmrg.py` had never been ported to `backend.py`, and
`gs_energy()` with `gs_method="idmrg"` raised before its first growth step
on any device, twice over:

* `np.take(T.array, idx, axis=...)` in `_project_channel`. With an integer
  `idx` NumPy hands the call to the array's own `.take` method with its own
  `mode="raise"`, which `jnp.take` does not implement --
  `NotImplementedError`. Now `bk.xp().take`.
* `arr[idx] -= shift * arr[src]` in `_subtract_energy_baseline`. JAX arrays
  are immutable. Now `bk.setblock`, which is the same in-place write on
  NumPy.

Both are the traps `docs/documentation.md`'s GPU section lists; neither
changes a NumPy result (216 host iDMRG/infinite-chain tests pass
unchanged). `tests/test_pyitensor_gpu_compatibility.py` pins the fix --
it fails with the original `NotImplementedError` without it -- and pins
the other rows of the table above.

### Running is not the same as staying on the device

A NumPy free function that falls back to `__array__` copies the whole
tensor to the host and returns the right answer, so nothing in the table's
first two columns can see it. JAX can:
`jax.transfer_guard_device_to_host("disallow")` raises at the first
device-to-host transfer and `"log"` reports every one with its shape. The
column above comes from running each calculation under that guard with the
engine's *designed* synchronizations allowed (`backend.to_host`,
`backend.scalar`, `ITensor.scalar` -- see `backend.py`'s docstring), so
anything left is a transfer nobody intended. Whether one matters is
decided by whether it repeats, which is measured by doubling the iteration
count and counting again:

| calculation | undesigned transfers, 1x -> 2x iterations | what moves |
|---|---|---|
| TDVP-GSE (nt 10 -> 20) | 57 -> 57 | `gse.py`'s companion density matrix, only while bonds are still growing |
| iDMRG (20 -> 40 growth steps) | 159 -> 319 | site tensors (chi,d,chi), the chi x chi singular-value matrix, flat two-site Krylov vectors |
| VUMPS (15 -> 30 iterations) | 6,556 -> 13,573 | `C`/`AC` blocks, into host-side eigensolvers |
| NH-DMRG (4 -> 8 sweeps) | 13,734 -> 25,038 | flat two-site vectors, once per Arnoldi matvec |

TDVP-GSE is effectively resident: its count does not grow with the time
step count, because `_gse_bond_step` only expands while the bond dimension
has room to grow (it becomes a per-step cost again only in a run whose
bonds never saturate). The other three are correct on a device but are
**not device calculations**: every iteration ships O(chi^2 d^2) data across
the bus and back, which is the per-call pattern this file's transfer table
shows losing to the host above chi ~ 64 however fast the device is. The
origin is design, not a slip: NH-DMRG's `_arnoldi_smallest_real` is written
over flat host NumPy vectors on purpose (its docstring says so), exactly as
`tdvp.py`'s Krylov propagator was before its 2026-08-26 port; iDMRG and
VUMPS lean on SciPy/ARPACK eigensolvers, which only take host arrays. Each
would need the same kind of port the Krylov propagator got. None is
pressing on a consumer card, where the device loses anyway -- but on a
data-centre GPU, do not expect `gs_method="idmrg"`/`"vumps"` or `nhdmrg`
to speed up on the device until that port exists.

## A trap that produced a completely wrong conclusion

**Do not benchmark a ground state on a uniform 1D Heisenberg chain.** Its
entanglement is S ~ (1/3) log L, so the ground state converges at chi ~ 60
and `maxm` above that does nothing: measured at n=32, E0 is identical to
13 digits at maxm=120 and maxm=240, and both CPU and GPU times were flat
across the whole sweep. An early sweep on it reported "ground state: no
speedup, 0.11-0.19x", which was purely an artifact of benchmarking a knob
that had no effect. The same port on a 3-leg ladder gives 20.27x.

Convergence measured at n=24 (|E(240)-E(120)|), i.e. how much bond
dimension each model actually needs:

| model | \|dE\| | CPU time chi=60 -> 240 | verdict |
|---|---|---|---|
| Heisenberg chain | 2.3e-14 | flat (~1.9 s) | useless for this |
| 2-leg ladder | 3.9e-11 | 4.1 -> 12.2 s | saturates by chi~120 |
| J1-J2 at J2=0.35 | 6.9e-13 | flat | useless for this |
| J1-J2 at J2=0.5 | 0 (exact) | flat | **Majumdar-Ghosh: exact chi=2 dimer product state** |
| **3-leg ladder** | **3.0e-06** | **9.4 -> 102.3 s** | still unconverged at chi=480 |

`benchmarks/gpu/port_speedup.py --model` implements these.

## How to measure this yourself, without fooling yourself

* `benchmarks/gpu/gpu_microbench.py` -- the primitives above at your own
  shapes; imports no dmrgpy at all.
* `benchmarks/gpu/port_speedup.py` -- end-to-end GS + KPM, with
  `--model`, `--pad-bonds`, `--backends`, and cold/warm reported
  separately.
* `benchmarks/gpu/kpm_gpu_probe.py` -- where KPM's time actually goes.
* Pin BLAS threads (`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1`) for every
  timing, GPU runs included.
* On a hybrid-core host, pin cores too (`taskset`) and interleave the A/B
  in one process: an unpinned comparison here reported 1.3-1.7x for a
  change a pinned one measured at 1.02x.
* Beware cProfile on this workload -- it inflated a run 1.71x and did so
  *selectively*, toward call-heavy code, which produced a wrong conclusion
  about where KPM spends its time. `kpm_gpu_probe.py` now reports both the
  clean and the profiled time so the distortion is visible.
* Sanity-check any GPU number with the site's job-accounting tool: a
  "GPU" run that never touched the device shows 0% utilization.
* Utilization is not enough to show a calculation *stays* on the device.
  Run it under `jax.transfer_guard_device_to_host("log")` with
  `backend.to_host`/`backend.scalar`/`ITensor.scalar` wrapped in
  `"allow"`, and count the log lines at two iteration counts: a number
  that doubles is a per-iteration round trip (see the
  device-compatibility section).
