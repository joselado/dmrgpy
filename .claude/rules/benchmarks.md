---
paths:
  - "benchmarks/**"
  - "src/dmrgpy/blasthreads.py"
---

# Cross-backend performance benchmarking (`benchmarks/`)

These are the notes CLAUDE.md used to carry under "Cross-backend
performance benchmarking"; they load when a file under `benchmarks/` or
`src/dmrgpy/blasthreads.py` is read or edited. The rule on BLAS threads is
still in CLAUDE.md, since it holds for every timing.

**Before benchmarking anything, check BLAS thread counts.** DMRG is a great
many small dense linear-algebra calls (a two-site tensor at `maxm=30` on
spin-1/2 sites is 60x60), and at that size more BLAS threads is slower, not
faster: measured 1 thread against 2, `svd` 1.6x and `eigh` 2.3x slower. End
to end that alone is minor (~1.13x for `"python"`, nothing for v3), but
under oversubscription (a shared cluster node, or several dmrgpy runs in
parallel, where every process's BLAS thinks it owns the machine) it
dominates every other effect measured in this repository: on a 14-core host with another
job holding 10 cores, MKL's default thread count turned a 9.4s `"python"`
solve into 28.2s and a 1.6s v3 solve into 26.7s. Pin threads before
importing numpy (`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 ...`) when
timing anything, or timings across runs are not comparable.
`src/dmrgpy/blasthreads.py` has the full measurements and an opt-in
`limit()` context manager (needs `threadpoolctl`; the env-var route is more
effective, since some BLAS libraries fix threading on first use). dmrgpy
never changes thread counts on its own. The same cap applies to agents: at
most three Python processes at once across all subagents of a hunt or
workflow (`run3.sh` in the 2026-09-24b record's "Shared helpers").

`benchmarks/run_benchmarks.py` answers a different question than
`tests/`/`examples/`: not "is this backend correct" but "which backend is
fastest, and by how much, for each kind of calculation". It runs on every
backend available in the current environment (`v2`, `v3`, `"python"`,
`"julia_live"`, auto-detected via `cppext.available()` and a live Julia
import) and sweeps two families of calculation over a uniform S=1/2
Heisenberg chain: **finite**, swept over chain length (ground-state energy,
static correlator, KPM and TD dynamical correlators, first excited state,
bond entanglement entropy, real-time evolution after a quench, and the
ground state confined to the `Sz=0` sector), cross-checked against ED at
every size cheap enough for it (the ED entropy reference is an SVD of the
ED state vector, since `entropy.py` is `self._session`-only and has no ED
path); and **infinite**, swept over bond dimension (`infinitechain.py`):
energy density and `<Sz>` under both `gs_method="idmrg"` and
`gs_method="vumps"`, checked against the Bethe-ansatz density and against
`<Sz>=0`. The catalogue of what is covered is `benchmarks/calctk.py`'s
`CALCS` list; adding a capability is one entry there plus a label in
`reporttk.CALC_LABELS`. Calculations only some backends implement (the
sector, METTS, the non-Hermitian ground state, SECTOR and TASK, on v3,
`"python"` and `"julia_live"`; the infinite chain, on v3 and `"python"`)
are restricted to those backends and show up as blank-by-construction
rather than as failures in the coverage table. The sector section tabulates
the sector-against-dense ratio per size (block sparsity scales, its
bookkeeping does not), and the `<Sz>` sections catch a wrongly gauged iDMRG
unit cell, which an energy-only check cannot see.

It writes `benchmarks/output/report.tex` (tables, matplotlib plots and an
auto-generated "which backend is worth optimizing" paragraph) and compiles
it to `report.pdf` with `pdflatex` if available; `benchmarks/output/` is
gitignored. Run `python run_benchmarks.py --help` from `benchmarks/` for
the full option list; `--from-json results.json` regenerates the report
from a previous run's raw data. `--repeats>1` re-prepares each calculation
from a fresh chain per repeat, because the DMRG session caches its
ground-state energy (`groundstate.py`'s send-cache) and timing a second
`gs_energy()` on the same chain would record that backend as instantaneous.
The script always runs one untimed warm-up call on the Julia backend before
the timed sweep, since `juliacall` compiles each Julia function signature
the first time it is called in a process.
