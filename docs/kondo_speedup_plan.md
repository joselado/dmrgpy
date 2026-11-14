# Plan: faster Kondo spectra

This plan collects the speedups that the Opus re-read of 2026-11-13 found
in the Kondo commits of 2026-11-10 to 13 (workflow `wf_5422a87c-92a`, five
readers, read only). Every gain in the items below is the reader's
estimate from the code and from timings already in the records. What was
measured after each change is in `.claude/rules/kondo.md` for the ED
route, and under "Measured on 2026-11-14" for the TASK route (items 6 to
10). Each item says where the time goes, what to change, what the gain
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

### Measured on 2026-11-14

Items 6 to 10 are done, in this order: 6, 8, 9, 7, 10, and two riders
(the residual of the block Krylov space on the device, the spectrum cache
in real storage). The times are local, one core (CPU 2 of a workstation
whose other cores ran other jobs), one BLAS thread, the minimum of two or
three repeats in one process with the TASK caches cleared before each,
on a ground state seeded so that every tree solves the same one. They
are not the Xeon numbers of the baseline table. The case numbers are
those of `run_case.py --list`.

| item | what it does | measured | curves |
|---|---|---|---|
| 6 | order=2 builds row and column 0 only | case 4: 2.57 s to 1.90 s | bit for bit |
| 8 | i S_y applied as a real operator, the Ritz transform real | v3, nkry=300, 0.3 K: build 7.56 s to 6.21 s (P S P 2.77 s to 1.50 s), call 15.9 s to 14.6 s; nkry=100: build 1.53 s to 1.44 s; "python" after item 7, nkry=300: projected spectrum 8.43 s to 7.35 s | 2e-16 of the peak |
| 9 | H overlaps off the projection pass, one gauge less | block_krylov 1.86 s to 1.69 s at nkry=100 and 9.48 s to 8.49 s at nkry=300 (complex basis), 5.07 s to 4.40 s at nkry=300 (real basis); case 4: 1.86 s to 1.66 s; case 5: 22.4 s to 21.4 s | converged poles 4e-15 of max abs Hq; curves 1.2e-15 (0.3 K), 7e-9 (T=0, loose poles) |
| 7 | the "python" state in a real gauge, H as a real MPO | build 12.97 s to 6.42 s at nkry=300 and 2.36 s to 1.39 s at nkry=100, store 32.1 to 16.0 MB; case 4: 1.75 s to 1.13 s; case 5: 21.3 s to 14.2 s; case 6: no change | moves: 3e-14 at nkry=300, 2e-7 to 1.5e-6 at nkry=100 and T=0 |
| 10 | the n_gs manifold kept on the chain (TASK only) | case 6: 17.8 s to 4.9 s, 4 solves and 8 builds to 1 and 2 | the first call bit for bit; later calls now take the first call's manifold (a new solve moved the curve by 4e-9 to 2.6e-8) |
| GPU rule | the residual of the block Krylov space formed on the backend | not timed; a test on JAX's CPU device reaches the branch and refuses a tangent vector on its way to the host | the same product on NumPy |
| cache | real and imaginary matrices kept in real storage | half the memory of an entry of a real basis (35 MB at nkry=400) | the same values handed out |

Not done: a bound of the spectrum cache in bytes (it stays a count of 8,
since a bound in bytes is more than a one-line change); the real gauge for
`get_poles` and `dynamical_correlator` with `submode="TASK"`, which go
through `taskdmrg.session_poles` and would move their numbers too; and
the timings on Triton.

Item 8 gives 1.06 to 1.22 on the build alone and less on a call, and
item 9 1.10 to 1.15 on the block Krylov space and 1.12 on case 4; both
are kept, since the rule to revert a change below 1.1 was read on the
function the item aims at, and together they are 1.3 on the real
"python" build at nkry=300. The plan's
estimates held for 6 and 7 (1.35 against 1.5 for 6, which the share of
the P S P loop bounds, and 2.0 to 2.2 on the build for 7) and for 9; for
8 the third-order sums did not get faster, although they now run in real
arithmetic, since the F tables dominate them at these sizes.

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
