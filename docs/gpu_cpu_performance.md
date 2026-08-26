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
