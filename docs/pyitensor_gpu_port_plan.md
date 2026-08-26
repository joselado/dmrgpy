# Plan: porting `pyitensor` to GPU, and measuring what it buys

Written 2026-08-25 as a plan; **Phases 0-3 and 5 are now done and
measured** (2026-08-26), and the first two items of Phase 6's follow-up
list -- the dispatch floor and TDVP -- landed the same day, so this file
has become the design-and-history record. The results a *user* needs are in
`docs/gpu_cpu_performance.md`; what remains to do, and what deliberately
will not be done, is Sec. 9.

Status at a glance:

| | state |
|---|---|
| Phase 0, measure before porting | done -- gate passed at chi=64, Sec. 5b |
| Phase 1, `backend.py` + `tensor.py` + `svd.py` | done |
| Phase 2, device residency in DMRG | done, measured in Sec. 8 |
| Phase 3, the KPM path (incl. energy truncation) | done, measured in Sec. 8 |
| Phase 4, benchmark campaign | done (`benchmarks/gpu/`) |
| Phase 5, docs / tests / example | done |
| Phase 6, follow-ups | dispatch floor + TDVP **done**; METTS/tebd/gse open; iDMRG/VUMPS deliberately not planned -- Sec. 9 |

Site-specific operating detail (how these jobs were submitted on one
particular cluster) is deliberately kept in this checkout's untracked
local notes rather than here.

## 1. Why `pyitensor` is the right target

`pyitensor/` is the pure-Python/NumPy DMRG backend
(`itensor_version="python"`): a dense, always-complex128 reimplementation
of the ITensor v3 API subset that `mpscpp3/chain_session.h` uses. It has
**no block sparsity and no JIT**, which is exactly why it is slower than
compiled ITensor — and exactly why it is the backend a GPU helps most:
every contraction is one big dense `tensordot`, and every truncation is
one big dense factorization, with no block bookkeeping to get in the way.
The compiled C++ backends (`mpscpp2`/`mpscpp3`) are not candidates —
vendored upstream ITensor, treated as black boxes here.

The numeric surface is remarkably small, which is what makes this
tractable at all:

| file | lines | NumPy/SciPy call sites | role |
|---|---|---|---|
| `tensor.py` | 299 | **9** | the ITensor type: `zeros`, `asarray`, `transpose`, `tensordot`, `conj`, `eye`, one `np.number` isinstance |
| `svd.py` | 244 | 17 | the one factorization primitive (`svd`, Gram+`eigh` route, `qr`) every truncation in the package goes through |
| `kernels.py` | 445 | 18 | the planned transpose+reshape+matmul matvec chain used by DMRG's Lanczos and TDVP's Krylov step |
| `dmrg.py` | 681 | 27 | sweeps, environments, Lanczos |
| `tdvp.py` | 368 | 11 | Krylov exponentiation (shares the `kernels.py` matvec) |

Everything else (`idmrg.py` 98 sites, `vumps.py` 55, `vumps_ms.py` 53,
`metts.py` 33, `chain.py` 33, ...) sits *on top* of those primitives and
is explicitly out of scope for the first pass (see §7).

Nothing in `tensor.py`/`svd.py` mutates an array in place — every
operation returns a new array — so an immutable-array backend (JAX) is
not structurally excluded, and neither is any other.

## 2. What is already there (and one hazard)

`kernels.py` already carries an **opt-in JAX path** for the DMRG/TDVP
matvec: the whole contraction chain expressed as one fused
`jax.numpy.einsum`, JIT-compiled once at module level and keyed on a
deterministically-labelled subscript string so structurally identical
bonds share a compile-cache entry. `jax_enable_x64` is already set, so
complex128 is preserved.

**Both points below are now resolved** -- the auto-enable is gone and the
fused kernel is superseded by the resident-array port (Sec. 8) -- but the
reasoning is kept because it is why the port is shaped the way it is.

Two things follow:

* The fused-einsum kernel is a *ready-made* GPU matvec. We do not need to
  design one; we need to stop paying the per-call host↔device conversion
  around it.
* **Hazard to fix early.** `USE_JAX` defaults to
  `_detect_default_use_jax()`, which auto-enables JAX whenever
  `jax.devices()` reports any non-CPU device. Its docstring says this is
  "UNTESTED on real GPU/TPU hardware". On any cluster whose Python ships
  `jax_cuda12_plugin` (§3) it *will* fire, so a user running pyitensor on
  a GPU node silently gets a code path nobody has ever measured. Phase 0 either validates that heuristic or replaces it
  with an explicit opt-in.

Measured on CPU, the JAX path is 1.5–2.5x *slower* than plain NumPy,
because per-call `numpy`↔`jax` conversion and dispatch exceed the compute
saved at these tensor sizes. On a GPU that same conversion becomes a
literal H2D/D2H transfer per matvec — which is the whole design principle
of this port:

> **Arrays live on the device. Host syncs are limited to the Lanczos
> α/β scalars and the O(χ) singular-value vector `_truncate` needs.**

## 3. Environment assumptions

The measurements below were taken on an HPC GPU node; the site-specific
details (login, scratch layout, scheduler partitions, job scripts) live in
this checkout's untracked local notes rather than here, since they are
particular to one account and one cluster. What matters for reproducing
the numbers is only:

* Python 3.12 with `jax` 0.7.1 and its CUDA 12 plugin already installed --
  so there is a **zero-install GPU path**; `cupy` was *not* available, and
  `torch` 2.8 with CUDA 12.8 was.
* **FP64 hardware.** Everything in this engine is complex128 (ZGEMM), so
  the measurements assume a data-centre GPU with real double-precision
  throughput (V100/A100/H100/H200 class). On a consumer card with
  1/32-rate FP64 none of these numbers transfer.
* Because complex128 is ~4x the real-FLOP count of a same-shape real
  matmul, the crossover bond dimension sits higher than an FP32 ML
  workload would suggest.
* CPU baselines are run as separate jobs from GPU ones (many sites forbid
  CPU-only work on GPU partitions), BLAS-pinned per `CLAUDE.md`'s
  benchmarking section. The baseline CPU is therefore a different host
  from the GPU node's, which is acceptable for an "is the device worth it"
  comparison and is stated wherever a number is quoted.

## 4. Prior art (checked before designing)

* **cuDMRG** (ClarkResearchGroup) — a CuPy-based DMRG; the closest
  precedent for "swap the array module, keep the algorithm".
* **Renormalizer** — production TD-DMRG with a CuPy backend selected by
  an environment variable.
* **Cytnx** — tensor-network library whose selling point is moving a
  tensor to GPU without changing user code, i.e. device-as-attribute.
* **quimb + `autoray`** — the dispatch-design precedent: a thin layer
  that routes numpy-style calls to whatever array library holds the data.
  We deliberately do *not* take the dependency (`pyitensor` promises to
  work with only NumPy/SciPy) but we copy the idea in ~100 lines.
* **NVIDIA cuTensorNet / cuQuantum** — has MPS-DMRG primitives; relevant
  as a possible far-future backend, not for this port.
* **DMRG with Tensor Processing Units**, PRX Quantum 4, 010317 — the
  reference point for "accelerator DMRG only pays off above some χ".

Published GPU-vs-CPU speedups for MPS/SVD workloads cluster around
4x at χ≈256, ~4.4x at χ≈512, ~7x at χ≈2048 — i.e. **the win is real but
it is a large-χ win**, which is precisely the regime where dense
`pyitensor` is currently unusable. That is the honest framing for this
work: it does not make small chains faster, it makes large-χ chains
reachable.

## 5. Design

### 5.1 `pyitensor/backend.py` (new, ~100 lines)

A single module owning *all* array-namespace decisions:

```python
get_xp()                 # the active array module (numpy | cupy | jnp | torch-shim)
set_backend("numpy" | "cupy" | "jax" | "torch")   # explicit, process-wide
available()              # what is importable AND has a usable device
to_device(a) / to_host(a)  # the only two places a transfer may happen
```

Rules the module enforces, and which the rest of the package then gets
for free:

* NumPy stays the default. With the default backend, every code path
  must be **byte-identical** to today — same calls, same order, no shim
  overhead in the hot loop (resolve `xp` once per function, never per
  element).
* Only `to_host()` may pull data back. Anything that needs a Python
  float/complex goes through it explicitly, so every sync is greppable.
* dtype is complex128 everywhere, unconditionally, as today.

### 5.2 `tensor.py` — the ITensor type

Nine call sites become `xp` calls; `ITensor.array` simply holds a device
array. `__mul__`'s `tensordot`, `transpose`, `conj`, `zeros`, `eye`
all have direct equivalents in every candidate backend (torch needs
`permute`/`movedim` aliases — that is what the shim is for). `mul_plan()`
is pure index bookkeeping and never touches data: unchanged.

### 5.3 `svd.py` — the factorization primitive

`_svd_truncated()` already prefers the Gram route (`eigh` of the smaller
of `M M†`/`M† M`, then one matmul to recover the kept factor) over a full
`svd`, falling back to the exact SVD when the smallest kept singular
value is too small for the squared spectrum to resolve. That is
*fortunate*: dense complex128 `eigh` + `gemm` is much better GPU work
than a full SVD (cuSOLVER's complex128 SVD is the weakest primitive in
the set). The port keeps the same structure and adds one rule:

* `_truncate()` runs on the **host**, on the O(χ) singular-value vector
  only — it is a cumulative sum and a couple of comparisons, worthless on
  a GPU and full of data-dependent control flow. One D2H of χ floats per
  factorization is negligible against an O(χ³) factorization.
* The Gram/exact-SVD fallback threshold is a scalar comparison — also
  host-side, on values already transferred.

### 5.4 `dmrg.py` / `kernels.py` / `tdvp.py` — device residency

This is the phase that actually buys the speedup:

* Environment tensors (`_all_left_environments`/`_all_right_environments`)
  are built on device and stay there for the whole sweep.
* The Lanczos loop (`_lanczos_ground_state`) keeps its Krylov vectors on
  device. Per iteration it needs exactly two scalars on the host —
  `np.vdot(q, w).real` (α) and `np.linalg.norm(w)` (β) — plus the
  convergence test; the tiny tridiagonal eigenproblem
  (`_tridiag_ground_value`, an `eigvalsh` of an n×n matrix with n ≤ 30)
  stays on the host by design. Same for
  `_local_ground_state_penalized`'s overlap penalties.
* `kernels.py`: `_make_matvec_planned()`'s precomputed transpose+reshape+
  matmul plan is backend-agnostic already — it becomes the default GPU
  matvec. The JAX fused-einsum path becomes a *comparison* candidate
  rather than a default, and `_detect_default_use_jax()`'s silent
  auto-enable is replaced by an explicit opt-in unless Phase 0 vindicates
  it.
* TDVP comes nearly free: `_lanczos_expm_multiply` uses the same matvec
  builders and the same Krylov structure.

### 5.5 User-facing switch (as built)

Process-wide and explicit, in the established style of this package:

```python
from dmrgpy.pyitensor import backend
backend.set_backend("jax")        # default stays "numpy"
backend.set_pad_bonds(240)        # optional, see Sec. 8
```

No automatic device detection anywhere — a GPU must be asked for, never
inferred (see the §2 hazard, which was exactly that mistake made once
already). `set_backend` raises on an unknown name rather than falling back,
because a GPU run that quietly became a CPU run is a benchmark that lies.

Two departures from the plan as written. Only `"numpy"` and `"jax"` exist:
CuPy was never installed anywhere this ran, and torch turned out not to be
needed once the JAX dispatch floor was understood (Sec. 5b). And the
`Many_Body_Chain.setup_python(device="gpu")` convenience was not added --
the backend is process-wide state, so hanging it off a per-chain method
would imply a per-chain scope that does not exist.

## 5b. Phase 0 results (2026-08-25) -- GATE PASSED

Run on an NVIDIA H200 (141 GB, jax 0.7.1 + CUDA plugin, 1:50 wall, 74% GPU
utilization, 11.7 GB VRAM) against one Xeon Gold 6248 core. Raw JSON and
the figure come from `benchmarks/gpu/`'s own harness. The plan originally
targeted an A100; an H200 was what was free, and it is a better FP64 part,
so the measurement is on stronger hardware than planned rather than
weaker.

**Primitive microbenchmark, complex128, H200 vs one Xeon 6248 core (and vs
8 cores in parentheses):**

| op | chi=64 | chi=128 | chi=256 | chi=512 | chi=1024 |
|---|---|---|---|---|---|
| matvec | 5.5x (5.0x) | 36x (21x) | 246x (127x) | 480x (169x) | **688x (166x)** |
| eigh (Gram route) | 1.2x (1.0x) | 4.6x (2.4x) | 13x (4.2x) | 40x (10x) | 121x (27x) |
| svd (fallback) | 0.6x (0.6x) | 1.0x (0.7x) | 2.7x (1.0x) | 6.9x (1.8x) | 24x (3.6x) |
| qr | 1.2x (1.0x) | 6.2x (3.9x) | 16x (7.7x) | 47x (18x) | 145x (40x) |

The matvec numbers are not too good to be true: at chi=1024 the chain is
~3.4e11 real FLOPs, so 8.3 ms on the H200 is ~41 TFLOP/s (near its FP64
peak) and 5.7 s on one core is ~60 GFLOP/s (about right for one Cascade
Lake core). The GPU is simply delivering the FP64 throughput it advertises
against a serial baseline.

Also as predicted: **`svd` is the one primitive the GPU is bad at**
(slower than CPU below chi=256), which is exactly why `svd.py`'s
preference for the Gram+`eigh` route matters here rather than being a
mere CPU optimization. Below chi=256 the GPU matvec is flat at ~0.35 ms --
that is dispatch latency, not compute.

**The transfer result, which is the design-defining one:**

| chi | matvec on device | one H2D+D2H round trip | |
|---|---|---|---|
| 64 | 0.34 ms | 0.31 ms | compute-bound |
| 512 | 1.52 ms | 4.31 ms | **transfer-bound (2.8x)** |
| 1024 | 8.31 ms | 21.96 ms | **transfer-bound (2.6x)** |
| 2048 | 60.3 ms | 82.7 ms | **transfer-bound (1.4x)** |

Above chi=64, shipping theta to the device and back costs more than the
matvec itself. Any design that converts per call is therefore unwinnable
regardless of device speed -- the "arrays stay resident, only scalars come
back" principle of Sec. 5 is now measured rather than assumed.

**The USE_JAX hazard is confirmed and was fixed on the spot.** On the GPU
node `kernels.USE_JAX` did default to True, and unmodified pyitensor was
**5-11x slower** for it: 0.26s vs 2.91s (n=16, maxm=60), 0.48s vs 2.25s
(n=20, maxm=100), 0.53s vs 3.50s (n=24, maxm=200), energies agreeing to
8e-14. `_detect_default_use_jax()` now always returns False, following the
instruction its own docstring left for this outcome. This is independent
of the port: it was making every GPU-node user slower today.

**Verdict against the gate** (matvec crossover at chi <= 512 -> proceed):
crossover is at **chi=64**, the smallest size measured, on both the
1-thread and 8-thread baselines. Gate passed with a wide margin; Phase 1
is justified.

**One honest caveat carried into Phase 1.** The probe's own chains
(n=24, maxm=200) finish in ~0.5 s on CPU, i.e. real pyitensor runs at
these sizes never reach the chi where the GPU's advantage is dramatic.
The value of this port is not "make today's small runs faster" -- it is
"make large-chi runs, where dense pyitensor is currently unusable,
reachable". Phase 4's benchmark campaign must sweep chi, not just chain
length, or it will measure the wrong axis and understate the result.

## 5c. The KPM dynamical correlator is a different problem (measured)

Ground-state DMRG and the KPM dynamical correlator do *not* share a cost
profile, so the Sec. 5b verdict does not transfer to KPM and it gets its
own measurement (`benchmarks/gpu/kpm_gpu_probe.py`).

DMRG's sweep is dominated by the two-site matvec -- one big dense
contraction, the primitive the GPU wins most decisively (688x at
chi=1024). KPM's moment recursion (`pyitensor/chain.py`'s
`_kpm_moments_full`/`_kpm_moments_accelerated`) is instead a long serial
chain of

    a' = 2*(M a) - a_prev,   truncated back to kpmmaxm after every step

so each of its O(100-1000) moments costs one MPO application plus one MPS
sum, and *both* truncate through `svd.py` -- the primitive the GPU is
**worst** at (slower than CPU below chi=256 on an H200). Profiled on CPU
at n=12, kpmmaxm=40, ~150 moments: 78% of the run is the moment
recursion, 49% inside `svd()`, 25% in the Gram-route `eigh` alone, and
19% in `tensordot` spread over 17711 calls.

**A correction, twice over, because the first two versions of this
section were both wrong -- and the way they were wrong is worth keeping.**

Version 1 bucketed a cProfile run by file, found "55-63% host
bookkeeping", and derived an Amdahl ceiling of ~1.8x on any GPU port.
That was an attribution error: cProfile only sees *callables*, so
arithmetic invoked through operators (`am @ bm`, `mat.conj().T`,
`np.sqrt`, fancy slicing) has no Python frame and is charged to its
caller's self time. `svd.py`'s Gram and recovery matmuls -- the very work
a device accelerates -- were counted as bookkeeping.

Version 2 fixed the attribution but kept trusting the profiled *totals*,
reporting tensordot at 18-35% of the run and predicting that removing it
would pay accordingly. It did not: a pinned, interleaved A/B measured
1.02x, flat across kpmmaxm=40, 80 and 160. The reason is the second
distortion: **cProfile inflates this workload**, and inflates it
selectively. Measured directly, the same kpmmaxm=160 run takes 102.7 s
clean and 175.9 s under the profiler -- 1.71x -- because the run is
millions of small Python calls and cProfile charges each one. That
overhead lands precisely on the call-heavy buckets, so tensordot's
apparent share was inflated along with it.

Both effects push the apparent *Python-side* share up. Correcting for
them, and re-measuring after `ITensor.__mul__` stopped calling tensordot
(which by itself cut profiler inflation to 1.08x, since the new path
makes far fewer instrumented calls), the picture is:

| kpmmaxm | KPM wall (clean) | profiler inflation | eigh | pure bookkeeping | not separable |
|---|---|---|---|---|---|
| 40 | 5.8 s | 1.08x | 31% | **5%** | 47% |
| 80 | 21.0 s | 1.01x | 36% | **2%** | 55% |

The "not separable" column is `svd.py`/`tensor.py`/`chain.py` self time:
inline BLAS (the Gram matmul, the recovery matmul, the contraction matmul)
mixed with those files' own Python logic, which cProfile cannot split. But
the column that *is* clean -- Python-only overhead attributable to
`index.py` -- is **2-5%**, not 55-63%.

So the conclusion is the opposite of version 1's: **the KPM recursion is
overwhelmingly dense linear algebra** (~85-90% once the inline BLAS in the
inseparable column is accounted for), with `eigh` alone at 31-36% and
identifiable Python overhead in the low single digits. It is a legitimate
GPU target -- with the standing caveat that its factorization lean hits
`eigh`, which the H200 wins by only 1.0-2.4x at kpmmaxm-scale matrices
versus 4-27x at chi >= 256. No Amdahl ceiling is quoted: the honest
version of that number needs Phase 2's end-to-end device measurement.

Two things follow for the CPU side, both now done:

* `ITensor.__mul__` executes its contraction as an explicit
  transpose+reshape+matmul instead of calling `np.tensordot`. Isolated on
  the shapes the KPM recursion produces, that is 1.4-2.05x faster *with
  the plan built inside the call* (tensordot re-derives its own axes
  normalization and dispatches through `np.dot` every time), and it is
  bit-identical -- verified over a ten-shape zoo including scalars, vector
  dots and multi-axis permuted contractions, with outer products delegated
  back to `tensordot` because they alone rounded differently (5e-16).

  End to end it is worth **1.02x**, measured with a `taskset`-pinned,
  interleaved, multi-repetition A/B in one process against the same seeded
  random MPS start: 1.024x (n=10, kpmmaxm=40), 1.035x (n=12, kpmmaxm=40),
  1.021x (n=12, kpmmaxm=80), 1.020x (n=12, kpmmaxm=160), spectra
  bit-identical in every case. Keep it -- free, exact, and it benefits
  every contraction in the engine rather than only KPM -- but it is a 2%
  change, not a "cheaper win" that could substitute for the port.
* Three other candidates were measured and rejected: computing only the
  top-`keep` eigenvectors (`eigvalsh` + `eigh(subset_by_index=...)`) is
  0.38-0.88x, i.e. *slower*, because both passes pay the same O(m^3)
  tridiagonalization; scipy's alternative Hermitian drivers (`evd`,
  `evr`, `ev`) are all within 3% of numpy's `zheevd`; and the band-edge
  DMRG that costs ~18% of a first KPM call is already cached per chain
  (`_bandwidth_min`/`_bandwidth_max`).

Two methodological notes worth keeping, since both cost real time here:
`cProfile` percentages on this workload must be read against the clean
wall time it also distorts (the probe now measures both and prints the
inflation factor); and the first end-to-end numbers for the
`ITensor.__mul__` change looked like 1.3-1.7x purely because of CPU core
placement on a hybrid host (4 P-threads at 4.8 GHz, 10 E-cores at 2.4-3.8
GHz, DMRG ~2x slower on E silicon). Only a pinned, interleaved,
same-process A/B is trustworthy on this class of machine -- the same
caution `CLAUDE.md` gives for BLAS threads, applied to core placement.

**Tests for the KPM path exist now**, written so a port can be validated
against them (`tests/test_kpm_dynamical_correlator_python.py`):
pointwise agreement between pyitensor and the compiled v2/v3 backends
(measured 5.4e-12, asserted at 1e-8 -- far sharper than the pre-existing
peak-position check, since all three run the same recursion); the exact
zeroth-moment sum rule (integral of S_AA == mu_0 == 1/4 for Sz or Sx on a
spin-1/2 site, which is sensitive to every moment rather than just the
dominant pole); reproducibility across freshly built chains (~2e-14, so a
GPU-vs-CPU diff at 1e-10 is meaningful); and non-negativity of the
reconstructed spectral function.
`examples/dynamical_correlator/kpm_sum_rule_convergence` plots that sum
rule converging with kpmmaxm (2.7e-3 at kpmmaxm=10 down to 3.7e-6 at 80)
next to the lineshape against ED.

## 6. Phases (as planned; see the status table at the top and Sec. 9 for
where they actually landed)


### Phase 0 — measure before porting (one short GPU job, ≤30 min)

No dmrgpy code changes. Deliverable: a decision, not a port.

1. **Environment probe**: does `jax` see the GPU, does `torch.cuda`
   work, on both an A100 and a V100 (the `sm_70` question of §3).
2. **`benchmarks/gpu_microbench.py`** (new, standalone): complex128
   `tensordot`/ZGEMM, `eigh`, `qr`, `svd` at the shapes DMRG actually
   produces, sweeping χ ∈ {64, 128, 256, 512, 1024, 2048} with d=2 and
   MPO bond w=5, on each candidate backend, against a BLAS-pinned CPU
   baseline. Reports per-op time *and* measured H2D/D2H transfer cost.
   → gives the **crossover χ** and the **backend choice**.
3. **Unmodified pyitensor on a GPU node**: `gs_energy` at, say,
   n=20/maxdim=100 and n=24/maxdim=300, with `kernels.USE_JAX` forced
   True and forced False, profiled. Tests the §2 hazard directly and
   shows whether transfers dominate as predicted.

**Backend decision, settled (job 19936675, H200).** JAX, with the
requirement that hot kernels go through `jax.jit`. Eager `jnp` has a
per-call dispatch floor of ~0.35 ms that made it look 3-4x worse than
torch at small chi; jitting the same kernel removes it entirely:

| chi | jax (jit) | jax (eager) | torch |
|---|---|---|---|
| 32 | 0.090 ms | 0.345 ms | 0.100 ms |
| 64 | **0.072 ms** | 0.360 ms | 0.098 ms |
| 256 | 0.287 ms | 0.375 ms | 0.274 ms |
| 1024 | 8.31 ms | 8.37 ms | 7.99 ms |

So the eager-vs-torch gap was a property of how JAX was called, not of
JAX. At large chi all three converge (same cuBLAS ZGEMM underneath).
`kernels.py`'s existing fused-einsum-under-`jax.jit` design is therefore
the right shape for the port's matvec, and torch stays a measured
fallback rather than a needed one. The three namespace divergences torch
would have forced on the `xp` shim (`.size` is a method on a tensor,
`numpy` refuses a cuda tensor without `.cpu()`, `tensordot` spells its
axes `dims`) are recorded in `benchmarks/gpu/gpu_microbench.py`'s own
adapter; JAX needed none of them.

Backend decision rule, applied to Phase 0's numbers:

| candidate | install | API fit | note |
|---|---|---|---|
| **JAX** (`jnp`) | none — plugin already present | near-NumPy | immutability is fine (§1); fused-einsum kernel already written |
| **CuPy** | `pip install cupy-cuda12x` into a `$WRKDIR` venv — **needs approval** | true drop-in | cleanest `xp` story; the prior-art default |
| **torch** | none | needs a ~50-line adapter | strongest fallback if the other two disappoint |

Gate: **if the microbench shows no crossover below χ≈512 on an A100, stop
and report** rather than porting. pyitensor's realistic χ range is the
deciding fact, and it is cheap to establish first.

### Phase 1 — `backend.py` + `tensor.py` + `svd.py`

The `xp` abstraction lands with NumPy as the only backend initially, so
the whole existing test suite must pass unchanged and bit-identically.
Then the chosen GPU backend is registered and `tests/test_pyitensor_gpu.py`
(skipped without a device) checks agreement.

### Phase 2 — device residency in DMRG

Environments + Lanczos on device, syncs reduced to α/β. Success metric:
GPU utilization on a real `gs_energy` run is not transfer-bound
(`nvidia-smi` during the run; whatever job-accounting tool the site
provides afterwards).

### Phase 3 — TDVP, and truncation polish

`tdvp.py` (shares the matvec), plus whatever the Phase-2 profile says the
new bottleneck is — likely the factorization path in `svd.py`.

### Phase 4 — benchmark campaign

`benchmarks/run_benchmarks.py` already sweeps chain length across
backends and cross-checks against ED; extend it with a `--device`
axis rather than writing a second harness, and add χ as a sweep axis
(the GPU story is a χ story, not an L story). Report: speedup vs χ, vs GPU model (V100 / A100 /
H100), crossover point, and the accuracy cross-check at every point.

### Phase 5 — docs, tests, examples (repo conventions)

* `docs/documentation.{md,tex}` — architecture: a new device axis under
  backend dispatch. `docs/user_guide.{md,tex}` — the user-facing switch.
  Both formats kept in sync; `.tex` must still compile under `pdflatex`.
* `tests/test_pyitensor_gpu.py` — skips cleanly with no device.
* `examples/backend_comparison/pyitensor_gpu_VS_cpu/main.py` — **must
  plot** (CLAUDE.md rule): speedup vs χ, with the CPU and GPU curves
  overlaid, plus the ED/CPU agreement check.

### Phase 6 (explicitly deferred)

`idmrg.py`, `vumps.py`, `vumps_ms.py`, `metts.py`, `idmrg_window.py`,
`idmrg_excitations.py`, `nhdmrg.py`, `kpm_energy_truncation.py`. These
are ~5x the numeric surface of Phases 1–3 combined and they all sit on
top of `ITensor`/`svd`, so they inherit the device automatically for
anything expressed through those; only their *own* direct `np.` calls
(the counts in §1) need auditing. Do them profile-driven, one at a time,
or not at all — this plan is not a 13k-line rewrite.

## 7. Verification

* **Correctness**: GPU vs CPU agreement to ~1e-10, never bit-identical —
  different summation order in the reductions, plus DMRG is iterative.
  Ground-state energies additionally cross-checked against `mode="ED"` on
  small chains, which is the existing pattern in `tests/_helpers.py`.
* **No regression on CPU**: the full suite (`pytest tests -k "not
  julia_live"`, plus the Julia tests only if Julia code is touched —
  it isn't here) must pass with the default NumPy backend, and the
  `"python"` backend's timings must not move.
* **Benchmarks**: BLAS threads pinned (`MKL_NUM_THREADS=1
  OMP_NUM_THREADS=1`) in every job — mandatory on a shared node, see
  `CLAUDE.md`. GPU utilization checked with the site's job-accounting
  tool, so a "GPU" number that never touched the device is caught.

## 8. What the finished port measures (2026-08-26)

Full tables, primitives and pitfalls: `docs/gpu_cpu_performance.md`. The
short version, NVIDIA H200 against one Xeon Gold 6248 core, warm:

| calculation | small chi | large chi |
|---|---|---|
| KPM correlator, n=30 | 0.13x at kpmmaxm=40 | **10.79x at kpmmaxm=240** |
| ground state, 3-leg ladder n=30 | 0.41x at maxm=60 | **20.27x at maxm=480** |

The result is **not** "correlators yes, ground states no" -- it is "large
chi yes, small chi no", identically for both. GPU time is nearly flat
across each sweep (365 -> 586 s; 25.1 -> 22.9 s) while CPU time explodes
(49 -> 6325 s; 10 -> 465 s), so the device is still dispatch-bound at the
top of both ranges and both figures are lower bounds. Crossover is
chi ~ 120-160.

Three findings that were not in the plan and that changed how this is
used:

* **Chain length hurts, bond dimension helps.** The KPM crossover moved
  from kpmmaxm ~ 80 at n=16 to ~120 at n=30: more sites means more
  operations, each paying the ~0.35 ms eager-dispatch floor, while more
  bond dimension means more arithmetic per operation at constant call
  count. This is the opposite of the intuition that bigger systems favour
  a GPU.
* **Cold versus warm is a factor of 10, not a detail.** XLA compiles a
  kernel per (operation, shape) and DMRG mints a new shape whenever a bond
  dimension changes: 44.9 s cold versus 1.72 s warm on one small
  ground-state solve. `backend.set_pad_bonds(K)` freezes every bond at K
  and collapses that (40676 -> 10407 compilations), which on the device is
  worth 1.4-3.1x on a cold run for 6-31% on a warm one -- and on the host
  is a pure ~2x loss. So: pad for one-shot jobs, not for sweeps. This knob
  came from the user, not from this plan.
* **The ground-state benchmark in Sec. 5b's own gate was on the wrong
  model.** A uniform 1D Heisenberg chain converges at chi ~ 60 (E0
  identical to 13 digits from maxm=120 to 480), so an early sweep on it
  measured dispatch overhead only and reported "ground state: no speedup".
  The same port gives 20.27x on a 3-leg ladder.
  `benchmarks/gpu/port_speedup.py --model` now carries the working
  benchmark and the negative controls, J2/J1=0.5's exact chi=2 dimer
  product state included.

Correctness at every measured point: ground-state energies agree with the
host to <= 2.2e-11, KPM sum rules to <= 7.9e-07, every spectrum satisfies
the exact 1/4 zeroth-moment sum rule on both backends, and the
energy-truncation Krylov projection agrees to 1e-10 at the unit level.

## 9. What is left, in the order worth doing it

Ranked by measured payoff per unit of work, not by how much code each
touches. The counts below are direct `np.` call sites, and
"silent-transfer" counts `np.<free function>(array, ...)` patterns -- the
ones that copy a device array to the host per call with no error (Sec. 5's
first trap), i.e. the ones that make a port real work rather than a
namespace swap.

1. **[done 2026-08-26] Lower the dispatch floor: `jax.jit` on the hot
   kernels, now that `set_pad_bonds` can stabilize their shapes.** Not a
   port at all, and the highest-leverage item on this list. The floor is
   why nothing wins below chi ~ 120, and a jitted kernel measured
   0.072 ms against eager's 0.345 ms -- 5x. Lowering it does not speed up
   one calculation, it widens the set of problems for which the device is
   worth using at all. `jax.jit` was ruled out in Sec. 5.4 *because* it
   retraces per shape and DMRG's growing bonds defeat it; padding removes
   exactly that objection, so the two compose and neither is much use
   alone.

   *What landed.* `backend.jit(fn, static_argnums)` registers a hot
   composite and returns a wrapper that compiles it when jitting is on and
   calls the plain function otherwise -- so the NumPy path executes the
   same operations in the same order as before, bit-identically, which is
   the invariant Phase 1 set. `backend.set_jit()` decides: `"auto"` (the
   default) is on exactly when `set_pad_bonds` is, `True`/`False` force
   it. `backend.compilations()` reports the traced-kernel count, i.e. the
   shape zoo padding exists to collapse. Four composites cover everything
   the engine does per Lanczos iteration:

   | composite | eager dispatches | where |
   |---|---|---|
   | planned matvec chain | 3-4 per step, 8-12 per call | `kernels.py::_matvec_chain` |
   | contraction: transpose+reshape+matmul | 5 | `tensor.py::_contract_matmul` |
   | Gram matrix + `eigh` + descending sort | 5 | `svd.py::_gram_spectrum` |
   | Lanczos recurrence / reorthogonalization | 3 / 2 per basis vector | `dmrg.py` |

   The last one is a change of algorithm on the device only: the host
   keeps its modified-Gram-Schmidt loop (marginally more stable, free
   there, and bit-identical to the pre-port code), while a device backend
   projects against the whole basis at once, since the loop otherwise
   issues 2k tiny kernels per iteration -- at small chi, more than the
   matvec they orthogonalize. `make_matvec` also now refuses its two
   opt-in host paths (numba contractor, fused einsum) whenever
   `backend.is_device()`: both convert per call, so either would
   reinstate the very round trip this port removed.

   *Measured*, `benchmarks/gpu/jit_speedup.py`, n=6 Heisenberg, maxm=16,
   JAX **on CPU**, pinned to one P-core. Ground state (4 sweeps) and a
   TDVP quench (8 steps), every row returning the identical value:

   | case | GS cold | GS warm | kernels | TDVP cold | TDVP warm | kernels |
   |---|---|---|---|---|---|---|
   | eager, unpadded (the port as shipped) | 34.8 s | 0.66 s | - | 59.9 s | 2.10 s | - |
   | padded, eager | 22.2 s | 1.04 s | - | 31.4 s | 3.54 s | - |
   | **padded + jit** | **7.7 s** | 0.80 s | 86 | **11.1 s** | 2.41 s | 93 |
   | jit, unpadded | 12.6 s | 0.48 s | 266 | 31.3 s | 1.32 s | 432 |

   Read those honestly. The cold column is the point: 4.5x and 5.4x. The
   padded/unpadded jit rows show exactly the interaction the ranking
   predicted -- unpadded, jit traces 3-4.6x the kernels (266 vs 86, 432 vs
   93) and gives back most of the win, which is why `"auto"` ties the two
   knobs together.
   The warm column barely moves because on a *host* the eager dispatch
   floor is microseconds, not 0.35 ms -- the floor this item attacks is a
   device property, so the host numbers are a lower bound.

   *Measured on the device* (H200, jobs 19951209/19951210, n=20, same
   script; full table and caveats in `docs/gpu_cpu_performance.md`):
   padded+jit is **6.4-12.1x cold** against eager/unpadded, roughly twice
   the host effect, on both the ground state and a TDVP quench; padding
   alone gives 2.4-3.8x and jit alone 2.4-2.7x, so both knobs are needed
   (jit alone traces 897-1092 kernels against 114-192 padded). Warm it is
   1.1-1.6x and the ranking *flips* -- padding does real arithmetic on
   known-zero blocks, so once everything is compiled it is jit without
   padding that wins. Hence the operational rule: pad+jit for a script
   that runs once, jit alone for a long sweep in one process.

   What it does **not** do is move the crossover. One CPU core runs those
   same sizes in 0.79 s / 1.24 s warm, ~7x faster than the best device
   configuration there. Lowering the floor makes the device much cheaper
   to start; it does not make a small problem worth a device. The claim
   at the top of this item -- that this "widens the set of problems for
   which the device is worth using at all" -- is therefore *not*
   supported at n=20: what widened is the range of *session shapes*
   (one-shot scripts, cold runs) rather than the range of physics.
2. **[done 2026-08-26] TDVP (`tdvp.py`).** The best remaining physics
   target: real-time evolution grows entanglement with time, so chi climbs
   into the hundreds *by construction*, unlike a 1D ground state. It
   already shares `kernels.py`'s ported matvec, so this is the largest
   payoff per line changed of anything left.

   It turned out to be smaller than its 10-`np.`-site count suggested,
   and the count was measuring the wrong thing. Everything underneath the
   propagator was *already* resident -- `two_site_heff`/`one_site_heff`/
   `zero_site_heff` come from `dmrg.py` and go through
   `kernels.make_matvec`, the environments come from `dmrg.py`'s ported
   builders, the truncations from `svd.py`. The single unported piece was
   `_lanczos_expm_multiply`, the Krylov exponentiator, and it was
   host-bound in the way that costs most: `np.linalg.norm`, `np.vdot` and
   a preallocated `Q = np.empty((n, m))` buffer written column by column,
   i.e. a device->host->device round trip *per Krylov iteration, per bond,
   per time step*, silently -- it never raised, it just transferred.

   Now the basis stays where `v0` lives and exactly two numbers per
   iteration (alpha, beta) come home, the same split
   `_lanczos_ground_state` uses; the k x k projected tridiagonal matrix,
   its eigendecomposition and Saad's residual test stay on the host by
   design. The two storage layouts are behind `_KrylovBasis`: the host
   keeps its preallocated buffer and view-slices (immutable device arrays
   would make each column write an O(n m) copy), a device keeps a list and
   stacks the k built so far, O(n k), one dispatch. `tests/
   test_pyitensor_gpu_backend.py` asserts the residency directly through
   the returned array's *type*, because agreement alone cannot catch a
   regression to transferring -- a round trip returns the same numbers,
   only slower. Cross-backend trajectory agreement: 1.5e-14.
3. **METTS (`metts.py`, 31 / 0 / 3).** Probably the biggest absolute win
   in the library -- samples x time steps x Krylov iterations -- and it
   inherits TDVP's large chi, so it sits directly on (2), now done.
   Smaller than its count suggests, in the other direction from TDVP:
   most of those `np.` sites are RNG, pooled-variance statistics over
   samples, and a d x d single-site `eigh` for the collapse basis, all of
   which *should* stay on the host. Audit before porting.
4. **Cheap completeness: `tebd.py` (6 / 0 / 0), `gse.py` (7 / 0 / 0).**
   No traps in either. `kpm_energy_truncation.py` is already done.

Not worth doing now, with reasons rather than silence:

* **`idmrg.py` (97 np sites, 14 silent-transfer, 17 in-place) and
  `vumps.py` (54 / 18 / 4).** Those silent-transfer counts are the
  pathology that made the pre-existing per-call JAX path 5-11x slower, so
  this is a genuine porting job; and `mpscpp3`'s own iDMRG observables
  are already 5-500x faster than the Python ones, so the incentive is
  weak from both ends.
* **`nhdmrg.py` (22 / 3 / 6)** -- niche, and carries traps.
* **The ED path (`edtk/`)** -- small dense matrices by construction; it
  is the small-system fallback, and a device can only lose there.
* **More 1D ground-state work.** The physics caps it: an area law with
  log corrections keeps chi small, so ground states only benefit for
  ladder/2D-like geometries, which the current port already covers.

## 10. Running the benchmarks

`benchmarks/gpu/` holds the harness: `gpu_microbench.py` (primitives),
`pyitensor_gpu_probe.py` (unmodified pyitensor on a device),
`kpm_gpu_probe.py` (the KPM cost split), `port_speedup.py` (end-to-end
GS + KPM speedup, `--model`/`--pad-bonds`/`--backends`) and
`jit_speedup.py` (the dispatch-floor knobs: eager/padded/jitted x ground
state/TDVP, cold and warm, with the traced-kernel count). Each is a plain
script with `--help`; none of them needs a scheduler, and the job scripts
that submitted them on one particular cluster are intentionally not
tracked here.

Two rules that matter more than the scripts:

* pin BLAS threads (`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1`) for every
  timing, GPU runs included -- see `CLAUDE.md`;
* on a hybrid-core host, pin cores as well (`taskset`), and interleave the
  A/B in one process. An unpinned comparison here reported 1.3-1.7x for a
  change that a pinned, interleaved one measured at 1.02x.
