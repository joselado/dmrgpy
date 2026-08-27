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

The one exception, and the reason it is an exception, is worth reading
before you generalize from that: the four-point correlator wins on the
device at chi = 20, because `ctmode="batched"` gets its arithmetic-per-
dispatch from a *tuple batch* rather than from bond dimension. See "The
four-point correlator: batching beats bond dimension" below. The rule is
really "the device needs enough work per array operation", and chi is only
the usual way to get it.

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
exact in representation, not in trajectory.

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
