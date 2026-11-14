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

### Measured on the Krylov route (2026-11-14)

On the local workstation, one core (CPU 2) and one BLAS thread, minimum of
a few runs; another session also ran Python on the machine, so the times
are good to about 5%. Case 7 is `krylov_n12_T1K`, and the new case
`krylov_n12_T1K_B0` is the same chain at zero field.

- **Item 11, done.** Three of the six dense products per initial state are
  gone (the explicit residual, H Q formed again, the projection of S-),
  and the matrices are built in Krylov coordinates. The Krylov basis is the
  same to the last bit, and dI/dV moved by 4.7e-16 of the peak on 10 spins
  and 5.9e-16 on 12 (eigsh seeded). `state_spectrum` went from 0.67 to
  0.60 s per state on 12 spins (1.12x). On the harness (`run_case.py
  --inline`, two runs each) case 7 went from 65.1 to 60.2 s (1.08x) and
  the zero-field case from 58.7 to 53.4 s (1.10x). The rest of the Krylov
  time is the loop that reorthogonalizes against each earlier block, a
  Python loop of small products. One block Gram-Schmidt
  over the whole space takes `state_spectrum` to 0.39 s (1.7x, case 7 1.2x
  to 1.4x), but it changes the basis at roundoff and moved dI/dV by
  2.9e-13 of the peak on 10 spins and 3.0e-13 on 12, above the 1e-13 this
  item was given. It was the branch `krylov-block-gram-schmidt`, one
  commit on item 11, and it is on master since 2c41a22 (the author's
  pick). The truncated-space test of item 11 then failed, since the
  matrices moved by 3.6e-13 of their largest element against the 1e-14 it
  asserted, and its bound is now 1e-12. For scale: the
  old code with its three start vectors in another order moves by 2.5e-14
  on 10 spins and 6.7e-12 on 12, and two runs of the same code on 12 spins
  differ by 9e-12, since SciPy's eigsh draws its start vector from fresh
  entropy.
- **The real route of the sums was already taken.** The Krylov spectra of
  a real H have had the form `conductance._spin_triplet` reads as real
  since `84090c0`; a test now pins it. That is why the ED speedups gave
  case 7 only 51.3 to 49.8 s: 90% of the third-order term of one Krylov
  spectrum (0.51 of 0.56 s on 12 spins) is the two exchange tables of F,
  one per direction, over n_f x dim with n_f = 248 of 301. Item 2 above,
  measured and not applied on the ED route where n_f is far below dim,
  would halve them here.
- **Item 12, not on the delivered branch.** The stop (tol = 1e-8 relative
  to onenorm(H), every Ritz state up to max|eV| + 40 kT above the initial
  state converged) is on the branch `krylov-residual-stop`. Where it fires
  it is exact enough: on 10 spins in the field of case 7 at 0.2 K and
  |eV| < 1 mV it stopped at 373 and 376 of 1024 states at nkry=400, 3.5x
  faster, and dI/dV moved by 1.8e-15 of the peak against the exhausted
  space (the acceptance asked for 6.8e-9). But at the default nkry=100 it
  fires on none of the chains measured: at 1 K and |eV| < 4 mV every one
  of the 29 states of the 10-spin chain was exhausted before its window
  converged, and the narrow window needs more than 300 states. There it
  costs 8% of `state_spectrum` on 12 spins (0.59 to 0.64 s per state) and
  saves nothing, so it is below the 1.1x of the brief. The default call
  differs from the exhausted space by 5.2e-13 of the peak on 10 spins at
  1 K, and by 2.7e-15 in the narrow window: at these sizes nkry=100 is
  already converged for dI/dV, long before every Ritz state of the window
  is. So item 12 stays on its branch, measured and not merged (the
  decider's pick). The measurement that would reopen it as an opt-in
  keyword is one at 12 spins or more with nkry from 200 to 400, where the
  Krylov space does not exhaust the Hilbert space. At 10 spins the 1024
  states are the whole space, which plain ED covers. The author picked on
  2026-11-14 to take it as an opt-in keyword instead, and it is on master
  as `tol=` (None by default), ported onto the block Gram-Schmidt loop,
  where a check needs no product in the Hilbert space. On 8 spins with
  `nkry=60` and `tol=1e-8` it stopped at 121 of 181 states at 0.1 K and
  |eV| < 0.5 meV and at 31 at T=0 and 0.3 mV, and dI/dV moved by 4e-12 and
  1e-11 of the peak.
- **`su2=None` by default, taken on the cut of the full sum** (the
  decider's pick, `ded5439` and `100eb84`). With the first cut the gates
  gave: (a) on against off at the default `p_cut` on 7 spins at 1.5 K
  agree to 2.1e-15 of the peak (44 states, 14 multiplets, the same
  levels); (b) on the zero-field 12-spin chain at 1 K on and off differ by
  2.4e-8 of the peak at nkry=100 and 200 alike, while nkry=100 and 200
  differ by 3.9e-12 (off) and 2.9e-12 (on). The cause was the weight cut,
  not the reduction: it kept a multiplet when (2S+1) exp(-x) was above
  `p_cut` times the largest weight, the full sum a state when exp(-x) is
  above `p_cut`, and two multiplets (S=2 and S=1 at 18.5 kT, exp(-x) of
  9.4e-9 and 8.4e-9) were kept by the first and dropped by the second;
  eigsh found every member of every level. Gate (a) passed only because
  no level of the 7-spin chain sits between the two cuts at the default
  `p_cut` (at `p_cut=1e-14` two multiplets do). The reduction
  now keeps a multiplet when exp(-x) is above `p_cut`, whatever its 2S+1,
  which moves `su2=True` results by about `p_cut` (2.4e-8 of the peak on
  the 12-spin chain, in the ledger). On that cut (a) gives 2.1e-15,
  (b) 2.5e-12, below the difference between nkry=100 and 200, and (c),
  auto against off at T>0, 2.5e-15 and 2.1e-15 on the 7-spin chain
  (`p_cut=1e-14` with nkry=40, and the defaults) and 2.3e-15 and 1.9e-15
  on the 3-site Hubbard chain with no Nf sector. Every `su2=True` test
  passes. The zero-field case `krylov_n12_T1K_B0` of the harness went
  from 53.4 s to 18.9 s (2.8x) with the default call. On a chain in a
  field the check fails at the first component of the total spin, 9 ms
  on 12 spins, and the call is the full sum as before.

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
pick. On 2026-11-14 a stand-in decider (a Fable agent, the author being
away) picked all three, and the author confirmed all three the same day:
`su2=True` by default when H commutes with the total spin (a
(2S+1)-fold saving at T>0, 3.2x in states measured on 7 spins at 1.5 K;
taken as `su2=None`, see the Krylov section);
the reuse of the n_gs manifold (item 10), which fixes the solver noise
between calls (taken, keyed on the content of the solved state); and a
default `mem_gb` for order="LL" on the CPU (not taken: a warning instead,
see the measured section). On 9
spins with a transverse field one flow of 241 biases needs, by the
memory model, 2.8 GB at 0.3 K, 34 GB at 1 K and 182 GB at 2 K (complex
count, about half on the real route), so today a T>0 call above 8 spins
needs `mem_gb` by hand.

The Krylov work of 2026-11-14 left three more picks, each measured in the
section above. The decider took `su2=None` by default together with the
weight cut of the full sum (2.8x on the zero-field 12-spin chain,
`su2=True` results move by about `p_cut`), and the author confirmed it.
The author picked the other two on 2026-11-14: the residual stop of
`krylov-residual-stop` is on master as an opt-in keyword, off by
default, so the default call keeps its cost and its numbers; and the
block Gram-Schmidt of `krylov-block-gram-schmidt` (1.7x on
`state_spectrum`, dI/dV moves by 3e-13 of the peak) is on master with
the truncated-space bound of the item 11 test loosened to 1e-12, still 9
times below the 9e-12 that two eigsh starts give on 12 spins.

Three more picks of the author, the same day. The rule to revert a change
below 1.1x is read on the function the item targets, not on the whole
call, so items 8, 9 and 11 stay (8 and 9 give 1.3 together on the real
"python" build at nkry=300, and 11 removes code). Item 2 is applied on
the Krylov route only, where n_f is close to dim, and kept only if it
passes 1.1x on the third-order term on Triton; on the ED route it stays
out. It is on master as `transpose_exchange=True`, which only
`krylov_kondo_dIdV` sets: on one Krylov spectrum of 10 spins (301
states) the third-order Kondo term went from 0.83 to 0.49 s at 1 K,
where every state is final, and from 0.19 to 0.17 s at T=0, where 150
are (one P-core of the workstation). On Triton it gave 1.49x on the
term of the 12-spin chain in a field (see the measurement of the author's
picks below), so it stays. dI/dV moved by 2e-16 of the peak at 1 K and
not at all at T=0. The saving is n_f/(2 dim) of the
F values, half at 1 K here and 41% at 248 of 301, so it is less than the
halving the plan wrote for the term. The 63 commits of the night were pushed on 2026-11-14.

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
still take (the decider's pick of a warning, not a default budget,
confirmed by the author).

## The TASK route, measured after the changes

Items 6 to 10 are on master (merged in `792279c`). The same cases as the
baseline, before (`85ea31a`) against after in one job on one core of
Triton (21576083, three or two runs each in fresh processes, the
minimum):

| case | before | after | factor | largest change of the curve |
|---|---|---|---|---|
| TASK, order 2, nkry=100, 121 biases | 2.33 s | 1.26 s | 1.85 | 4.9e-6 of the peak |
| TASK, order 3, nkry=300, 0.3 K, 241 biases | 20.8 s | 12.1 s | 1.72 | 2.7e-14 |
| TASK, n_gs=2, a scan of 4 values of U | 17.8 s | 4.8 s | 3.7 | 1.1e-5 |

Each run solves its own state, so at T=0 and for the n_gs manifold the
change between trees includes the solver noise of two separate solves
(two v3 solves of one chain differ by 3e-3 of the peak, see the rules
file); the agent's comparison on one seeded state gave 2.2e-7 for item 7
at T=0. In the scan the manifold is now solved once (get_excited_states
23.5 s to 6.2 s under the profiler) and the spectra are built twice
instead of eight times.

## The Krylov route, measured after the changes

Item 11 and the su2 default are on master (merged in `a3eb80a`), before
(`85ea31a`) against after in one job on one core of Triton (two runs
each, the minimum):

| case | before | after | factor | largest change of the curve |
|---|---|---|---|---|
| KRYLOV, 12 spins, 1 K, a field, 61 biases | 70.7 s | 69.5 s | 1.02 | 4.7e-16 of the peak |
| KRYLOV, 12 spins, 1 K, zero field, 61 biases | 63.4 s | 22.5 s | 2.8 | 3.5e-12 |

The node of this job was slower than that of the baseline (51.3 s for the
first case). With a field nothing but item 11 acts, and on this node it
gave 1.05 on `state_spectrum` and 1.02 on the call, below the 1.1 of the
rule; it stays because it removes code (the explicit residual and the
second product of H). At zero field the reduction to one member per
multiplet does the work, and the change against the full sum is the
2.5e-12 of gate (b).

## The author's Krylov picks, measured

The block Gram-Schmidt (`2c41a22`), the opt-in residual stop (`a9d8d06`,
off here) and item 2 on the Krylov route (`267a711`), against the pushed
master (`90f4ddb`), all four trees in one task on one core of Triton
(21578925, csl48, two runs each in fresh processes, the minimum; the
third-order Kondo term and `state_spectrum` from the profiled run):

| case | tree | call | `state_spectrum` | third-order Kondo term |
|---|---|---|---|---|
| KRYLOV, 12 spins, 1 K, a field | 90f4ddb | 52.5 s | 26.1 s | 25.3 s |
| | + block Gram-Schmidt | 44.7 s | 17.9 s | 25.3 s |
| | + residual stop, off | 44.3 s | 17.5 s | 25.4 s |
| | + item 2 | 36.4 s | 17.9 s | 17.0 s |
| KRYLOV, 12 spins, 1 K, zero field | 90f4ddb | 17.3 s | 8.6 s | 7.9 s |
| | + block Gram-Schmidt | 14.8 s | 5.9 s | 7.9 s |
| | + residual stop, off | 14.5 s | 5.9 s | 7.9 s |
| | + item 2 | 12.4 s | 5.9 s | 5.4 s |

The block Gram-Schmidt gives 1.46x on `state_spectrum` with a field and
1.46x at zero field (1.17x and 1.17x on the call), and item 2 gives 1.49x
and 1.45x on the third-order Kondo term (1.22x and 1.17x on the call), so
item 2 passes the 1.1 of the author's rule and stays. Together the call
goes from 52.5 to 36.4 s (1.44x) and from 17.3 to 12.4 s (1.40x). The stop
off costs nothing. The curve moves by 1.5e-12 of the peak with a field and
by 4.6e-13 at zero field against `90f4ddb`, the same to every digit in the
three new trees, so this is the change of basis of the block Gram-Schmidt
and item 2 adds nothing visible.
