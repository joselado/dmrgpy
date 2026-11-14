# Plan: faster Kondo spectra

This plan collects the speedups that the Opus re-read of 2026-11-13 found
in the Kondo commits of 2026-11-10 to 13 (workflow `wf_5422a87c-92a`, five
readers, read only). Every gain below is the reader's estimate from the
code and from timings already in the records. None of them is measured
yet. Each item says where the time goes, what to change, what the gain
should be, and whether a number moves. The three bugs and the test gaps
of the same re-read are already fixed (`27198eb`, `e8ac10a`).

## The baseline, measured

The light cases of `benchmarks/kondo_speedups/cases.py`, on HEAD of
2026-11-13 (`85ea31a`), one core of an Intel Xeon Gold 6248, BLAS
pinned, three or two repeats in fresh processes, and the cumulative time
of the functions the items aim at from one run under cProfile:

| case | time | where it goes |
|---|---|---|
| ED, 10 spins, T=0, 41 biases | 2.9 s | third-order Kondo term 1.66 s (F0 1.27 s), KondoSpectrum 1.20 s (eigh 0.63 s, eigenbasis 0.42 s) |
| ED, 10 spins, 1 K, 41 biases | 5.9 s | third-order Kondo term 4.28 s (the F spline 2.88 s), KondoSpectrum 1.21 s |
| ED, 9 spins, 20 K, 41 biases | 58.6 s | third-order Kondo term 55.7 s (coefficients 9.8 s, the F spline 2.7 s), potential term 1.6 s |
| TASK, order 2, nkry=100, 121 biases | 2.4 s | the projected spectrum 3.4 s of 3.5 s under the profiler |
| TASK, order 3, nkry=300, 0.3 K, 241 biases | 24.6 s | the build 14.4 s (block Krylov 11.5 s), third-order Kondo term 12.7 s |
| TASK, n_gs=2, a scan of 4 values of U | 20.3 s | get_excited_states 25.4 s of 32.4 s under the profiler |
| KRYLOV, 12 spins, 1 K, 61 biases | 53.7 s | the Krylov spaces 26.6 s, third-order Kondo term 26.1 s |
| LL on the CPU, 9 spins, T=0, 61 biases | 204.6 s | the flow, all of it |

The two heavy cases ran in a second job (21571380, same node type). ED
on 10 spins at 10 K with 241 biases took 2923 s and 2800 s in two runs,
with a peak of 6.0 GB, and the third run and the profile did not fit in
the 2 h of the job. The leading-log flow on 8 spins at 2 K with 61 biases
and `mem_gb=12` took 3209 s with a peak of 5.8 GB, 3153 s of it in the
flow and 36 s in the thermal kernel. The peak is half the budget, which
is what item 13 says: the memory model counts complex numbers where the
real route holds real ones.

So the ED sums (items 1, 2, 4), the manifold re-solve (item 10), the
Krylov spaces (item 11) and the leading-log flow (items 13 to 15) are
where the time is. The peak memory of every light case was below 0.9 GB.

We order the items by how much a user waits on the path they touch. The
ED route comes first, since it is the reference for every other route and
the only one at any T on small clusters.

## The ED Kondo route

These live in `kondospectrumtk/conductance.py` and `edkondo.py`, and also
serve `submode="KRYLOV"`, which feeds the same sums.

1. **Final states closed by the step.** The second- and third-order sums
   run over every final state f, also where the step Theta of eV - (E_f -
   E_i) is exactly 0 at T=0, or below Theta(-40) at T>0. Keeping only the
   states within max|eV| plus about 40 kT of the occupied ones makes the
   third-order Kondo term cost O(n_eV n_occ n_f dim) instead of
   O(n_eV n_occ dim^2), so a factor dim/n_f on that term, which is about
   60% of the 10-spin T=0 call of the hunt record. At T=0 no term changes;
   at T>0 the dropped terms are at most about 1e-13 of the result.
   `conductance.py:310`.
2. **The s->t exchange table is the transpose of the t->s one.** Both
   third-order terms build two dim^2 F tables per bias where one and its
   transpose suffice: about 1.7x on `third_order_kondo_dIdV` and 2x on
   the potential term, 30 s to 18 s at 12 spins on the Krylov route.
   Unchanged at T=0, where F0 is exactly even; spline rounding at T>0.
   `conductance.py:305`.
3. **Real arithmetic for a real H.** `eigh`, the three eigenbasis
   transforms and the Levi-Civita build run complex. Real arithmetic is
   about 3x on `eigh` and 4x on the rest, about 20 to 25% of a 10-spin
   call. Caching the transforms per ED object and site removes them from
   repeated calls at fixed H. Numbers move at roundoff, and inside a
   degenerate level the basis changes while every sum is invariant.
   `edkondo.py:72`.
4. **The blocked coefficient tensor is built once per group of biases.**
   With many biases at many occupied states it is built again for every
   group: about 12% at 512 states and 241 biases, up to half of the call
   at 917 of 1024 states. No number moves. `conductance.py:331`.
5. **The potential term is not chunked over biases**, and evaluates F on
   four arrays where two suffice by evenness. Chunking keeps the peak at
   about 0.1 GB, and the even form saves about 6% of a 9-spin 20 K call.
   `conductance.py:418`.

## The TASK route

These live in `pyitensor/task.py`, `taskdmrg.py` and `kondomethod.py`.

6. **order=2 builds three blocks it never reads**, the tangent-tangent
   P S P blocks: about 1.5x on the order=2 build at nkry=100. The result
   does not change. `task.py:1289`. This is the simplest item.
7. **"python" never takes the real basis.** pyitensor tensors are
   complex128 and the random start is complex, so a ground state of a
   real H comes in a complex gauge and `real_starts` returns early. A real
   gauge before the tangent space (as `nrgmps.real_gauge` does for the
   NRG hand-over), or a real start with real noise when H is real, should
   give 1.5 to 3x on the build and half the memory of the store. Numbers
   move by the imaginary part of the solved state, or by solver noise.
   `task.py:1102`.
8. **Sy and the Ritz transform run complex on real data.** About 15% of a
   T=0 build at nkry=100, and the transform grows as n^3. No number moves.
   `task.py:1296`.
9. **The Krylov loop computes the H overlaps twice** and gauges each
   vector once more than needed: 10 to 20% of `block_krylov` at nkry=300.
   Converged values do not move. `task.py:1221`.
10. **With n_gs>1 every call solves the manifold again.** A scan over
    biases, U or Jrho_s at n_gs=2 could reuse the members and their
    spectra: 10 to 20x on each repeat call. Repeated calls would then
    return the same manifold, where today they differ by solver noise (up
    to 3e-3 of the peak on a long chain at T=0). `kondomethod.py:510`.

## The Krylov thermal route

11. **`state_spectrum` does about twice the dense work it needs** per
    initial state: 1.5 to 2x on the Krylov part, 35 s of 65 s at 12 spins.
    `krylovkondo.py:198`.
12. **No convergence stop.** 3*nkry states are always built; a stop on
    the Ritz residual would cut up to (300/N)^2 per initial state, with dI/dV
    moving by the tolerance chosen. `krylovkondo.py:191`.

## The leading-log flow

13. **The CPU memory model counts complex numbers** while a real H runs
    real since `fce3382`, so `mem_gb` cuts about twice as many chunks as
    needed. `leadinglog.py:606`.
14. **The thermal kernel recomputes its D-independent parts** for every
    right-hand side and builds a 2 ms table: 1.3 to 3x at 7 to 8 spins at
    low T, nothing at 9 spins with many biases. `thermalkernel.py:84`.
15. **Sz block structure.** For a Hamiltonian that conserves Sz the four
    2d x 2d products and the hole sum are block sparse: 2 to 3x at 9 spins
    and T=0. `leadinglog.py:163`.
16. **On the GPU** the jitted right-hand side is compiled again for
    every flow (`leadinglog_device.py:346`), the tip sums re-read the whole
    operand for every bias (`:505`, 2 to 3x on the per-bias part at
    d>=4096), and the tip vertex keeps 20 stage copies (`:598`, about 1.6x
    at 13 spins, 0.5 K and 241 biases).

## Choices for the author

Three items change a default or a number, so they are not taken without a
pick: `su2=True` by default when H commutes with the total spin (a
(2S+1)-fold saving at T>0, 3.2x in states measured on 7 spins at 1.5 K);
the reuse of the n_gs manifold (item 10), which fixes the solver noise
between calls; and a default `mem_gb` for order="LL" on the CPU. On 9
spins with a transverse field one flow of 241 biases needs, by the
memory model, 2.8 GB at 0.3 K, 34 GB at 1 K and 182 GB at 2 K (complex
count, about half on the real route), so today a T>0 call above 8 spins
needs `mem_gb` by hand.

## How to measure

Every item is measured before and after on the call the record already
times: the 10-spin T=0 and the 9-spin 20 K `mode="ED"` calls of the
2026-11 hunt, the 20-site TASK chain at nkry=100 and 300, the 12-spin
Krylov chain at 1.5 K, and the 9-spin leading-log chain of the README.
Timings go to Triton as one Slurm array (one task per item and per
before/after), each with `-c 1` and `MKL_NUM_THREADS=1 OMP_NUM_THREADS=1`,
from `$WRKDIR/calculations/kondo_speedups`. The `#SBATCH` values are
checked against the Triton docs when the script is written, and the job
is submitted only after the author approves it. The local workstation runs
only the correctness checks, one at a time and below 3 GB.

## The ED route and the memory model, measured after the changes

Items 1, 3, 4, 5 and 6 are on master (`1e80085` to `31967be`, the record
of each in `.claude/rules/kondo.md`); item 2 was measured and not
applied, since after item 1 it saves only n_f/(2 dim) of the F entries
and gave 0 to 4%. The same cases as the baseline, before (`85ea31a`)
against after in one job on one core of Triton (21572809, three runs each,
the minimum):

| case | before | after | factor | largest change of the curve |
|---|---|---|---|---|
| ED, 10 spins, T=0, 41 biases | 3.20 s | 0.86 s | 3.7 | 5.8e-15 of the peak |
| ED, 10 spins, 1 K, 41 biases | 6.10 s | 2.37 s | 2.6 | 5.8e-15 |
| ED, 9 spins, 20 K, 41 biases | 58.4 s | 27.4 s | 2.1 | 2.0e-15 |
| ED, 10 spins, 10 K, 241 biases (one run each) | 2533 s | 802 s | 3.2 | 5.2e-15 |
| KRYLOV, 12 spins, 1 K, 61 biases | 51.3 s | 49.8 s | 1.03 | 2.3e-16 |

At T=0 the third-order term went from 1.84 s to 0.22 s and F0 from 1.38 s
to 0.17 s (item 1 keeps 99 of 1024 final states); at 20 K the coefficient
build went from 10.4 s to 1.4 s (item 3) and the rest of the third-order
term is the F tables of item 4. The peak memory went down in every case
(0.84 to 0.61 GB at 20 K, and 5.9 to 0.70 GB at 10 K with 241 biases, where
the hub blocks of item 4 take 796 of the 843 s of the call under the
profiler). The Krylov case does not gain, because its
spectra never take the real route and its time is in `state_spectrum`
(items 11 and 12).

Item 13 is applied too: `memory_model(d, m, real)` counts real numbers on
the real route. One flow of 241 biases (21572828), the peak above that of
the process before the call: 9 spins at T=0, 0.65 GB on the real route
against 0.65 GB of the new model (1.30 GB of the old count), and 1.32 GB
on the complex route against 1.30 GB; 8 spins at T=0, 0.24 and 0.47 GB
against 0.22 and 0.44 GB; 8 spins at 1 K with 55 occupied states and 61
biases, 1.78 GB on the real route against 1.99 GB. With `mem_gb=None` the
call now warns when the model puts the flow above what the process can
still take (the decider's pick of a warning, not a default budget).
