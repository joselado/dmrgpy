# Audit, 2026-11: five-lens hole hunt and profile of the Kondo machinery

This record is the ninth hole hunt over the Python layer. It covers
`get_kondo_spectrum` and everything under `src/dmrgpy/kondospectrumtk/`,
and it looked for two kinds of problem in the same lenses: numbers that are
plausible and wrong, and time or memory that a route spends without need.
It ran on 2026-11-11 and 2026-11-12 on `553a824` plus the uncommitted work
that added `T>0` to `submode="TASK"`. That work was committed during the
hunt as `b0e923d`, with the same content, so the hunt audited `b0e923d`.
Both compiled extensions were rebuilt before the hunt, since they were
older than the `chain_session.h` of `733472b`. The Kondo tests passed on
that build (193 passed, `julia_live` not run, and 8 of 8 `julia_live` Kondo
tests passed separately).

Five hunters ran, one per lens, and every repro was executed. They returned
19 candidates, 10 of them about cost and 9 about wrong behaviour. Eight
reviewers, each briefed to refute, took one unit of one to four candidates
whose fixes touch the same files, and judged each candidate separately.
This is fewer reviewers than one per candidate, chosen to keep the workflow
small; every verdict below is still per candidate. The reviewers refuted 2,
confirmed 6 and narrowed 11. What they struck stays visible in the
reviewer paragraphs. A speedup counted only with a profile of the public
call at a realistic size, a prototype, old against new on the same inputs,
and a median of three timed runs, and the reviewers measured it again.

Heavy runs went to d22, a 12-core workstation, through `rund22.sh` (see
"Shared helpers"). That runner had a defect during the hunt: it copied the
run folder back after each run, and the copy replaced the `.out` file that
`tee` was writing with the empty copy on d22. So most outputs of d22 runs
are empty on disk. For those runs the record carries the decisive lines as
the agent quoted them in its report, and says so where an output is empty.
The local runs are spliced from disk, byte for byte. `<scratch>` is the
session's scratch folder, `<repo>` the checkout and `<d22>` the folder on
d22; outputs longer than 160 lines keep their first 120 and last 40 lines.
Reviewer paragraphs are their hand-back text.

This file is the evidence, not a task list. Each finding has a
`**Status**` line and keeps its repro.

The fixes landed on 2026-11-12 in six clusters that touch separate files,
each in its own git worktree, plus two commits by hand, one per cluster
(`49c23ca`, `f6ac1e3`, `fce3382`, `4d1abbc`, `2393bc0`, `76313ee`,
`317c39f`, `d098cd0`). The worktrees were cut from `553a824` and not from
`b0e923d`. Five of the six commits touch files that `b0e923d` did not
change; the sixth, on `kondomethod.py`, was merged by hand. All of them
were tested again on top of `b0e923d`. 16 of the 17 findings are fixed and
one is recorded and not applied. The regressions are
`tests/test_audit_2026_11_<cluster>.py`. Five fixes change returned
numbers, all by far less than the method's own error except the two-time
term below saturation; they are in `docs/numerical_changes.md`.

## The five lenses

| Lens | Brief |
|---|---|
| `task-route` | `taskkondo.py` and the Kondo use of `taskdmrg.py` and `pyitensor/task.py`: TASK against ED at saturation at T=0 and at the new T>0, the cold-sample guard, the unconverged-pole warning, absolute thresholds, the spectrum cache; profile of the TASK call |
| `ed-reference` | `edkondo.py`, `conductance.py`, `stepfunctions.py`, the formulas every route shares: the thermal F, grids, degeneracies, occupied states, scale covariance; profile of the ED route at 8 to 10 spins |
| `leading-log` | `leadinglog.py`: the flow against the bare third order and the two-pass reference, scale covariance, degenerate states, the coupling check, keywords; profile at 6 and 7 spins |
| `correlator-twotime` | `secondorder_dc.py`, `potentialdc.py`, `dmrgtwotime.py`, `twotime.py`, `edtwotimeref.py`: the correlator routes and the two-time TDVP Kondo term against ED, the default grids; profile of the two-time term |
| `dispatch` | `kondomethod.py`: every keyword on every branch, validation order, the new T>0 dispatch, the default submode per backend, the texts; repeated work across calls and `n_gs` members |

## Scope

Out of scope by construction: vendored ITensor and TDVP; the legacy quirks
of `.claude/rules/cpp-backend.md`; the open `docs/known_issue_*.md` items;
anything in an earlier audit record; what `ROADMAP.md` marks as absent;
the deliberate departure from the paper's eq. 25; `ysr.py` and `ysrtk/`,
`kondoscreening.py`, `kondocloud.py`, `wilsonchain.py` and `nrgthermo.py`,
which are other methods. `itensor_version="julia_live"` was read and not
run, for its JIT cost. The measured and rejected optimizations of
`.claude/rules/kondo.md` (the GPU in float64, complex64 in the leading-log
flow, numba threads in `_hole_sum`, smoothing TASK steps by their residual,
1 GB TASK batches) were out of bounds.

Every Python process ran through one of two wrappers that share one cap of
three concurrent processes over all agents and both machines, pin the BLAS
threads to one and put this checkout's `src` first:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 <script>
```

## Where the time goes

Each hunter profiled the public call of the routes its lens owns, before
any fix. The profiles are the baseline that the speedups below are
measured against.

- **task-route**, get_kondo_spectrum(mode=DMRG, submode=TASK), T=0, order=3, U=0.3, 121 biases, 'python', uniform Heisenberg chain J=1e-3 eV with a field, maxm=40, nkry=100, ground state pre-solved (1.0 s, excluded), 20 sites, d22, 1 core: 16 s in total. TaskKondoSpectrum.__init__ 15.5 s = 97% (all in pyitensor/task.py projected_spectrum 15.46: block_krylov 10.24 = 64%, the 3-operator P S P loop 5.2 = 32%). tottime: numpy tensordot 7.3 s (46%, 85k calls from apply_many), backend.project_rows 3.5 s (22%), reshape 1.6 s, _overlap_columns 1.6 s (10%), stepfunctions.F0 0.3 s. The conductance sums are under 3% at T=0, nkry=100. Nothing proposed here: 107 apply_many calls, duplicate starts already deflated, block Gram-Schmidt would save at most project_rows (22%) and is untested. Log: <scratch>/hunt/task-route-logs/05_profile_n20_T0.out
- **task-route**, get_kondo_spectrum(TASK) at T=0.3 K, order=3, U=0.3, 241 biases, nkry=300, 'python', 20-site disordered XXZ chain (J in 0.5 to 1.5e-3, staggered field 6e-4), maxm=40, 20 sites, dim 904, d22, 1 core: 39.2 s in total. conductance.second_order_dIdV 19.3 s = 49% of the call, of which stepfunctions.Theta 17.9 s (46%) and _theta_raw 0.6 s; third_order_kondo_dIdV 11.1 s = 28% (FBuilder.__call__ 10.3 s: scipy interpolate _evaluate 5.5 s, F0 3.7 s); TaskKondoSpectrum build 8.6 s = 22% (project_rows 1.7 s, tensordot 3.7 s cumulative share of the run, _overlap_columns 1.1 s). Log: .../task-route-logs/12_profile_thermal.out
- **ed-reference**, mode=ED, T=0, order 3, U=0.25, 41 biases, S=1/2 Heisenberg chain + fields, n=10 spins (1024 states), d22, 1 BLAS thread: 3.84 s in total. third_order_kondo_dIdV 1.80 s (47%), of which F0 on dim^2 points per bias 1.27 s (33%); KondoSpectrum.__init__ 1.46 s (38%): _eigenbasis_spin_operator 0.86 s (22%, two complex 1024^3 products per operator) and eigh 0.51 s (13%); second_order_dIdV 0.56 s (15%); peak RSS 1.80 GB. n=8 (256 states): 0.18 s total.
- **ed-reference**, mode=ED, T=1 K, order 3, U=0.25, 41 biases, n=10 spins (1024 states), n_occ=4, d22, 1 BLAS thread: 8.43 s in total. second_order_dIdV 3.19 s (38%), almost all in Theta on 41x1024x1024 arrays, 2.83 s (34%); FBuilder.__call__ 2.97 s (35%) of which scipy spline evaluation 1.46 s (17%) and the F0 it computes first for every point and then overwrites 1.23 s (15%); KondoSpectrum.__init__ 1.46 s (17%); FBuilder table build 0.24 s (3%); peak RSS 3.26 GB. n=8: 0.73 s total.
- **ed-reference**, mode=ED, T=20 K (all states occupied), order 3, U=0.25, 41 biases, n=9 spins (512 states), d22, 1 BLAS thread: 10.6 s in total. third_order_kondo_dIdV 8.10 s (76%): _triple_product_coefficients 3.61 s (34%), F evaluation 2.60 s (25%), occupied-state einsums; potential term 1.42 s (13%); peak RSS 3.48 GB. n=10 at 10 K (917 of 1024 occupied, 21 biases): 77.7 s, 22.9 GB.
- **leading-log**, get_kondo_spectrum(order="LL", mode="ED", T=0), current code, 61 biases, 1 core, 6 spins 1/2 (d=64, 1 occupied state), d22, one core: 2.85 s in total. 626 rhs evaluations, 52 solver steps. _tip_bracket 1.28 s tottime (45%), _bracket 0.78 s (27%, 1252 calls = four 2d x 2d products per rhs), rk_step own time 0.26 s (9%), einsum 0.12 s, _unblock 0.085 s, _block 0.058 s. Median of 3 runs: 3.05 s.
- **leading-log**, get_kondo_spectrum(order="LL", mode="ED", T=0), current code, 61 biases, 1 core, 7 spins 1/2 (d=128, 1 occupied state), d22, one core: 13.5 s in total. 614 rhs evaluations, 51 steps (41 accepted). _bracket 5.90 s tottime (44%), _tip_bracket 4.52 s (33%; inside it _hole_sum is 62% locally: 12.8 of 20.9 ms at d=128), rk_step own 1.06 s (8%), einsum 0.45 s, _unblock 0.35 s, _block 0.27 s, ndarray.copy (the _swap copies) 0.19 s. Median of 3: 14.3 to 14.6 s. _bracket runs at about 56 GFlop/s of complex BLAS on one core, so it is compute bound.
- **leading-log**, same public call, profile of the whole call (19_public_profile.py), 7 spins 1/2, d=128, d22, one core: 14.4 s in total. leading_log_dIdV 14.387 s of 14.391 s: the ED diagonalization is cached and costs nothing here, 99.97% of the public call is _joint_flow (rhs 12.92 s: _bracket 6.5 s, _tip_bracket 5.3 s, rk_step 1.15 s own).
- **leading-log**, prototype (real block-form state, real arithmetic, precomputed hole-kernel constants), same call, 7 spins 1/2, d=128, d22, one core: 5.52 s in total. 518 rhs evaluations. tip_bracket 3.14 s tottime (57%, mostly _hole_sum), _bracket 1.34 s (24%), rk_step 0.38 s, einsum 0.25 s. This profile (11_profile_new.out) was taken with the real/block change and BEFORE llhole was wired in; with the hole precompute the median is 4.02 s.
- **correlator-twotime**, two-time Kondo term (v3 TDVP), fine grid dt2=dtau=1000, n_t2_half=3, n_tau_half=15, 41 eV, 8 sites S=1/2 Heisenberg + field, maxm=40, local, 1 thread: 55.9 s in total. cProfile: tdvp_step (C++) 52.43 s = 94%; K_F cell averages (scipy quad) 2.07 s = 3.7%; 2670 MPS overlaps 0.69 s = 1.2%. 1908 tdvp_step calls = 270 per t2 row (6 ket trajectories for G + 3 for Gx, 30 steps each) + 18 for the three t2 trajectories. File 03_prof_v3_n8_fine.out
- **correlator-twotime**, two-time Kondo term (python TDVP), same grid, n_t2_half=1, 21 eV, 8 sites, maxm=40, d22, 1 thread: 33.8 s in total. cProfile (inflated; plain run is 18.0 s): pyitensor tdvp_step 32.34 s = 96%, of which _lanczos_expm_multiply 18.8 s (21216 calls), _half_sweep_lr/rl 14.7 s each; K_F negligible. File 07_prof_python_n8_d22.out
- **correlator-twotime**, two-time Kondo term (v3) on the grid the docstring and tests use: dt2=25/Gamma0/10, n_t2_half=10 (here 3), dtau=2*pi/2e-5/15, n_tau_half=15, 41 eV, 4 sites, n_t2_half=3, local, 1 thread: 227 s in total. cProfile: K_F cell-average quad 192.95 s = 85% (3.5e6 K_F calls, IntegrationWarning, eV*dt2 of hundreds of rad per cell); tdvp_step 32.66 s = 14%. File 03_prof_v3_n4.out
- **correlator-twotime**, KPM correlator route: second_order_dIdV_dc + third_order_potential_dIdV_dc (v3), delta=1e-5, 40001 es points, 201 eV, 8 sites, maxm=40, local, 1 thread: 17.4 s in total. kpm_dynamical_correlator (6 calls, Sx,Sy,Sz twice) 13.95 s = 80%; F0 kernel (6 calls) 2.34 s = 13%; KPM generate_profile 0.72 s = 4%; einsum 0.06 s. File 11_prof_kpm_v3_n8.out
- **correlator-twotime**, same KPM route on python, 8 sites, maxm=40, delta=1e-5, 40001 es, 201 eV, d22, 1 thread: 41.8 s in total. cProfile: pyitensor SVD/contract inside the six KPM correlators about 95% (svd 17.95 s cumulative, tensor __mul__ 9.2 s); F0 1.19 s = 2.8%. File 18_prof_kpm_python_n8_d22.out
- **dispatch**, get_kondo_spectrum mode=DMRG submode=TASK, T=0.2 K, order=3, n_gs=1, repeat call (TASK spectrum served from cache), 10 sites S=1/2, maxm=30, nkry=60, 61 biases, local, 1 thread: 1.04 s in total. terms 0.618 (59%); FBuilder.__init__ 0.422 (41%, all in _F_convolve); third_order_kondo_dIdV 0.326 (31%); second_order_dIdV 0.236 (23%) of which Theta 0.212 (20%); CubicSpline eval 0.185 (18%); F0 0.097. Whole call at T=0 repeat is 0.32 s, T>0 repeat 1.39 s. (01_profile_small.out)
- **dispatch**, same, n_gs=2 (Kramers doublet), T=0.3 K, order=3, repeat call, 13 sites S=1/2, maxm=24, nkry=60, 81 biases, local, 1 thread: 33 s in total. get_excited_states 15.76 (48%, rebuilt on every call: dmrg_excited 13.8); two TaskKondoSpectrum builds 14.73 (45%, session_spectrum -> projected_spectrum 14.3; block_krylov 9.8 of it, apply_many 9.06); the remainder is FBuilder and the conductance sums. Unprofiled repeats take 16.6 to 24.8 s. (04_ngs_profile.out)
- **dispatch**, same as the first row, conductance stage only (the TASK spectrum was cached by three earlier refused calls), 24 sites S=1/2, maxm=40, nkry=100, 121 biases, T=0.1 K, d22, 1 thread: 1.84 s in total. Theta 0.766 (42%); CubicSpline _evaluate 0.386 (21%); F0 0.239 (13%); FBuilder.__init__ 0.218 (12%). The first call of the same chain, which builds the TASK spectrum, takes 25.8 s (06b_fbuilder_cache.out), so FBuilder is 0.9% of a first call. The TASK build itself is not profiled here (another lens). (05_profile_d22.out)

## The findings at a glance

The two-time grid finding (`correlator-twotime-3`) and the guard's blind
spot (`dispatch-2`) were refuted and are under "Ruled out".

| # | Finding | Kind | Severity | Review | Status |
|---|---|---|---|---|---|
| 1 | The second-order sum runs over every initial state on TASK, where only the ground state is occupied: 6.5 GB and 49% of a T>0 call at 904 states | performance | MEDIUM | confirmed | FIXED `49c23ca` |
| 2 | At T>0 the unconverged-pole warning promises a change of 4.1e-3 where nkry=100 and 300 differ by 6.7e-8 | bug | LOW | narrowed | FIXED `d098cd0` |
| 3 | The same sum on mode="ED": 9.5 GB at 1024 states and 121 biases | performance | MEDIUM | narrowed | FIXED `49c23ca`, numbers change |
| 4 | FBuilder evaluates F0 on every point and then overwrites it inside the table: 17.7% of a T>0 ED call | performance | LOW | confirmed | FIXED `f6ac1e3` |
| 5 | The third-order coefficient tensor is built for every occupied state at once: 15 GB at 917 of 1024 states | performance | MEDIUM | narrowed | FIXED `49c23ca` |
| 6 | The ED eigenbasis spin operator is densified first: 14.6% of a T=0 ED call | performance | LOW | confirmed | FIXED `f6ac1e3` |
| 7 | Theta_prime loses its digits for 1e-4 < abs(x) < 1e-2 (0.0 for 1/6), reached by tabulate=False only | bug | LOW | narrowed | FIXED `f6ac1e3`, numbers change |
| 8 | The leading-log flow runs in complex arithmetic for a real Hamiltonian: 2x to 2.5x of the call | performance | MEDIUM | confirmed | FIXED `fce3382`, numbers change |
| 9 | _hole_sum recomputes D-independent abs and sqrt at every step: 1.15x, and 1.37x on top of finding 8 | performance | LOW | narrowed | FIXED `fce3382` |
| 10 | A tip site with no spin gives NaN at every bias with a false Kondo-scale warning, for an exact 1.5708 | bug | LOW | narrowed | FIXED `fce3382`, numbers change |
| 11 | The two-time G runs six ket trajectories per row that three bra trajectories replace: 1.8x to 3.1x | performance | MEDIUM | confirmed | FIXED `4d1abbc`, numbers change |
| 12 | es starting at 0 cuts the elastic line in half with no warning: 0.000 at zero bias for an exact 1.037 | bug | LOW | narrowed | FIXED `2393bc0` |
| 13 | An integer dt2 raises a UFuncTypeError after the first row of TDVP | bug | LOW | confirmed | FIXED `4d1abbc` |
| 14 | Sz is run twice and the F0 kernel is built on twice the rows: 4% of the call | performance | LOW | narrowed | NOT APPLIED |
| 15 | Every T>0 call rebuilds the 0.2 to 0.6 s FBuilder table | performance | LOW | narrowed | FIXED `f6ac1e3`, `317c39f` |
| 16 | Three texts say julia_live gets the correlator route by default; the code sends it to TASK | bug | LOW | narrowed | FIXED `76313ee` |
| 17 | mode="ed" with order="LL" is told LL is ED-only; n_gs is checked after setup | bug | LOW | narrowed | FIXED `76313ee` |


## Findings

### 1. The second-order sum runs over every initial state on TASK, where only the ground state is occupied: 6.5 GB and 49% of a T>0 call at 904 states

`performance` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `task-route` &middot; candidate `task-route-1`

**Status**: FIXED in `49c23ca` on 2026-11-12, one change with finding 3. `second_order_dIdV` sums over the occupied initial states, with `p_cut=P_CUT` as the third-order terms have (`p_cut=0` is the full sum), and builds its Theta arrays per block of `CHUNK_ENTRIES` = 4e6 entries. On TASK the result is identical, since p is one-hot. Measured on d22, one core, on the 20-site `"python"` chain at `nkry=300` (904 states), with the old module swapped in on the same cached spectrum: the public T=0 call 1.24 s to 0.87 s (3.6e-16 apart), the sums at 0.3 K 4.06 s to 1.95 s (3.7e-16), the peak memory of the second-order term at 41 biases 2.38 GB to 0.04 GB. Pinned by `tests/test_audit_2026_11_conductance.py` (`::test_second_order_p_cut_zero_is_the_full_sum`, `::test_second_order_chunks_do_not_matter`, `::test_second_order_default_drops_only_empty_states`, `::test_second_order_onehot_state_memory`).

The claim as reviewed: conductance.second_order_dIdV sums over every initial state although p is one-hot on the TASK route, so its (n_eV,dim,dim) temporaries cost O(n_eV dim^2); restricting to the occupied columns and chunking over biases gives identical numbers on TASK (p exactly zero off row 0) and cuts 2.4 s and 6.5 GB at dim 904 and 241 biases (hunter) to 0.04 s and 0.04 GB (my mock). Same change as ed-reference-1.

**Where**: src/dmrgpy/kondospectrumtk/conductance.py:139-171 (second_order_dIdV: eps, weight, x_ts, x_st over all i); conductance._occupied_states (line ~190) is already used by the third-order terms but not here. Reached from taskkondo via kondomethod.py:339 and from edkondo via mode='ED'.

The third-order terms sum over the thermally occupied initial states only (P_CUT = 1e-12, O(n_occ dim^2)). The second-order term does not: it forms x_ts = eV[:,None,None] - eps[None,:,:] with eps of shape (dim, dim) and the same for x_st, then evaluates Theta on both, so it costs O(n_eV dim^2) in time and memory although on TASK p is one-hot on row 0 (and on ED at low T only a few rows are occupied). Expected: cost O(n_eV dim n_occ). Observed on a 20-site chain, python backend, TASK, nkry=300 (dim 904): second_order_dIdV 2.43 s median at T=0 and 6.52 GB peak numpy memory (tracemalloc), against 0.013 s and 0.04 GB for the occupied-columns version (script 09). At T=0.3 K the thermal Theta is evaluated on the full (n_eV,dim,dim) tensor: 17.9 s in Theta alone, 19.3 s of a 39.2 s call (script 12). With 6001 biases the old code raises ArrayMemoryError for 36.5 GiB (script 08, first attempt), so a dense bias grid on a long chain is impossible without a patch. End to end, public call with the spectrum cached, median of 3 on d22: nkry=100 old 3.40 s, new 1.87 s (45% saved); nkry=300 old 30.65 s, new 11.34 s (63% saved); max relative difference of the returned dI/dV 3.7e-16 and 7.3e-16 (script 10). It survived because every test uses dim of a few tens, where the (n_eV,dim,dim) tensor is small, and the 2026-11-10 batching work profiled the Krylov build, not the sums; the build is only 22% of the call at nkry=300 and T>0. The restriction is exact at T=0 (the dropped columns have p=0, so they contribute exactly 0) and drops at T>0 only states with p < 1e-12 of the largest, which the third-order terms already drop.

The hunter's claim, before review: conductance.second_order_dIdV sums over every initial state f,i of the spectrum although p is one-hot on the TASK route, building (n_eV, dim, dim) arrays: at dim=904 and 241 biases it takes 6.5 GB of numpy memory and 2.4 s at T=0, and it is 49% of the public T>0 call (19.3 of 39.2 s); restricting to the occupied columns saves 63% of the call and gives the same numbers to 7e-16.

Decisive output, quoted by the hunter:

```
dim 904 TaskKondoSpectrum build 6.7s
old median 2.430s peak numpy mem 6.52 GB
new median 0.013s peak numpy mem 0.04 GB
max rel diff new/old 5.09399898660089e-16
nkry 100 ... cached-spectrum call, median of 3: old 3.40s new 1.87s max rel diff 3.665859961715946e-16
nkry 300 ... cached-spectrum call, median of 3: old 30.65s new 11.34s max rel diff 7.331719923431889e-16
profile: second_order_dIdV 19.256 s cumulative of 39.182 s; Theta 17.926 s tottime; numpy._core._exceptions._ArrayMemoryError: Unable to allocate 36.5 GiB for an array with shape (6001, 904, 904)
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`09_second_order_occ.py`:

```python
# conductance.second_order_dIdV builds (n_eV, dim, dim) arrays over ALL initial states although p is one-hot on TASK
import sys, time, tracemalloc, numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import conductance as C
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
n=20; nkry=int(sys.argv[1]) if len(sys.argv)>1 else 300; neV=int(sys.argv[2]) if len(sys.argv)>2 else 241
rng=np.random.default_rng(5)
sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
J=1e-3*rng.uniform(0.5,1.5,n-1); h=0
for i in range(n-1): h=h+J[i]*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1])
for i in range(n): h=h+6e-4*(-1)**i*sc.Sz[i]
sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=40,14,1e-10
sc.gs_energy()
t0=time.time(); ks=TaskKondoSpectrum(sc,3,nkry=nkry); tb=time.time()-t0
print("dim",ks.dim,"TaskKondoSpectrum build %.1fs"%tb,flush=True)
eV=np.linspace(-6e-3,6e-3,neV)
def new(ks,eVs,T0=1.0,U=0.0):
    eVs=np.asarray(eVs,dtype=float)
    occ=C._occupied_states(ks)
    M2=C._spin_matrix_elements_squared(ks)[:,occ]
    eps=ks.e[:,None]-ks.e[occ][None,:]
    weight=ks.p[occ][None,:]*M2
    spin=np.empty(len(eVs))
    for a in range(0,len(eVs),32):
        e=eVs[a:a+32]
        spin[a:a+32]=(np.einsum('fi,efi->e',weight,C._theta_raw(ks,e[:,None,None]-eps[None]))
                     +np.einsum('fi,efi->e',weight,C._theta_raw(ks,-e[:,None,None]-eps[None])))
    return 2*np.pi*T0**2*(spin+4*U**2*np.ones_like(eVs))
for name,f in (("old",C.second_order_dIdV),("new",new)):
    ts=[]
    for r in range(3):
        tracemalloc.start(); t0=time.time(); d=f(ks,eV,U=0.3); ts.append(time.time()-t0); pk=tracemalloc.get_traced_memory()[1]; tracemalloc.stop()
    print(name,"median %.3fs"%np.median(ts),"peak numpy mem %.2f GB"%(pk/1e9),flush=True)
    if name=="old": d0=d
print("max rel diff new/old",np.max(np.abs(d-d0))/np.max(np.abs(d0)))
t0=time.time(); d3=C.third_order_kondo_dIdV(ks,eV,-0.1,omega0=20e-3,Gamma0=5e-5); print("third_order_kondo %.2fs"%(time.time()-t0))
t0=time.time(); d3=C.third_order_potential_dIdV(ks,eV,-0.1,0.3,omega0=20e-3,Gamma0=5e-5); print("third_order_potential %.2fs"%(time.time()-t0))
```

`10_end_to_end.py`:

```python
# public get_kondo_spectrum(TASK), old vs second_order_dIdV restricted to the occupied initial states
import sys, time, numpy as np, warnings
warnings.simplefilter("ignore")
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import conductance as C
exec(open("09_second_order_occ.py").read().split("for name,f in")[0].split("eV=np.linspace")[1].join(["eV=np.linspace",""]) if False else "")
n=20; rng=np.random.default_rng(5)
sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
J=1e-3*rng.uniform(0.5,1.5,n-1); h=0
for i in range(n-1): h=h+J[i]*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1])
for i in range(n): h=h+6e-4*(-1)**i*sc.Sz[i]
sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=40,14,1e-10
sc.gs_energy()
def new(ks,eVs,T0=1.0,U=0.0):
    eVs=np.asarray(eVs,dtype=float); occ=C._occupied_states(ks)
    M2=C._spin_matrix_elements_squared(ks)[:,occ]; eps=ks.e[:,None]-ks.e[occ][None,:]
    weight=ks.p[occ][None,:]*M2; spin=np.empty(len(eVs))
    for a in range(0,len(eVs),32):
        e=eVs[a:a+32]
        spin[a:a+32]=(np.einsum('fi,efi->e',weight,C._theta_raw(ks,e[:,None,None]-eps[None]))+np.einsum('fi,efi->e',weight,C._theta_raw(ks,-e[:,None,None]-eps[None])))
    return 2*np.pi*T0**2*(spin+4*U**2*np.ones_like(eVs))
old=C.second_order_dIdV
eV=np.linspace(-6e-3,6e-3,241)
kw=dict(site=3,Jrho_s=-0.1,U=0.3,T=0.3,omega0=20e-3,Gamma0=5e-5,order=3)
for nkry in (100,300):
    t0=time.time(); sc.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=nkry,**kw); print("nkry",nkry,"first call (spectrum build included) %.1fs"%(time.time()-t0),flush=True)
    res={}
    for name,f in (("old",old),("new",new),("old",old),("new",new),("old",old),("new",new)):
        C.second_order_dIdV=f
        t0=time.time(); _,d=sc.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=nkry,**kw); res.setdefault(name,[]).append(time.time()-t0); res[name+"_d"]=d
    C.second_order_dIdV=old
    print("  cached-spectrum call, median of 3: old %.2fs new %.2fs"%(np.median(res["old"]),np.median(res["new"])),"max rel diff",np.max(np.abs(res["old_d"]-res["new_d"]))/np.max(np.abs(res["old_d"])),flush=True)
```

`12_profile_thermal.py`:

```python
# cProfile of the public call at T>0, nkry=300, 20-site disordered XXZ chain, python backend, 241 biases, order 3 with U
import time, cProfile, pstats, io, warnings, numpy as np
warnings.simplefilter("ignore")
from dmrgpy import spinchain
n=20; rng=np.random.default_rng(5)
sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
J=1e-3*rng.uniform(0.5,1.5,n-1); h=0
for i in range(n-1): h=h+J[i]*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1])
for i in range(n): h=h+6e-4*(-1)**i*sc.Sz[i]
sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=40,14,1e-10
sc.gs_energy()
eV=np.linspace(-6e-3,6e-3,241)
kw=dict(site=3,Jrho_s=-0.1,U=0.3,T=0.3,omega0=20e-3,Gamma0=5e-5,order=3)
pr=cProfile.Profile(); t0=time.time(); pr.enable()
sc.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=300,**kw)
pr.disable(); print("total %.1fs"%(time.time()-t0))
s=io.StringIO(); pstats.Stats(pr,stream=s).sort_stats("tottime").print_stats(10); print(s.getvalue()[:3000])
s=io.StringIO(); pstats.Stats(pr,stream=s).sort_stats("cumulative").print_stats(14); print(s.getvalue()[:3500])
```

Observed:

`09_second_order_occ.out`:

```
  from dmrgpy import spinchain
dim 904 TaskKondoSpectrum build 6.7s
old median 2.430s peak numpy mem 6.52 GB
new median 0.013s peak numpy mem 0.04 GB
max rel diff new/old 5.09399898660089e-16
third_order_kondo 4.86s
third_order_potential 0.02s
```

`10_end_to_end.out`:

```
nkry 100 first call (spectrum build included) 4.6s
  cached-spectrum call, median of 3: old 3.40s new 1.87s max rel diff 3.665859961715946e-16
nkry 300 first call (spectrum build included) 38.6s
  cached-spectrum call, median of 3: old 30.65s new 11.34s max rel diff 7.331719923431889e-16
```

`12_profile_thermal.out`:

```
total 39.2s
         8533467 function calls (8527331 primitive calls) in 39.182 seconds

   Ordered by: internal time
   List reduced from 581 to 10 due to restriction <10>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
      126   17.926    0.142   17.926    0.142 <d22>/src/dmrgpy/kondospectrumtk/stepfunctions.py:58(Theta)
      248    5.501    0.022    5.502    0.022 <d22 home>/venv/lib/python3.12/site-packages/scipy/interpolate/_interpolate.py:781(_evaluate)
      248    3.740    0.015    3.740    0.015 <d22>/src/dmrgpy/kondospectrumtk/stepfunctions.py:111(F0)
   243247    2.124    0.000    3.661    0.000 <d22 home>/venv/lib/python3.12/site-packages/numpy/_core/numeric.py:1001(tensordot)
     1818    1.682    0.001    1.693    0.001 <d22>/src/dmrgpy/pyitensor/backend.py:260(project_rows)
      304    1.128    0.004    1.147    0.004 <d22>/src/dmrgpy/pyitensor/task.py:1123(_overlap_columns)
      248    1.013    0.004   10.266    0.041 <d22>/src/dmrgpy/kondospectrumtk/stepfunctions.py:214(__call__)
      126    0.611    0.005   18.537    0.147 <d22>/src/dmrgpy/kondospectrumtk/conductance.py:114(_theta_raw)
        1    0.607    0.607    8.494    8.494 <d22>/src/dmrgpy/pyitensor/task.py:1253(projected_spectrum)
   800841    0.603    0.000    0.603    0.000 {method 'reshape' of 'numpy.ndarray' objects}



         8533467 function calls (8527331 primitive calls) in 39.182 seconds

   Ordered by: cumulative time
   List reduced from 581 to 14 due to restriction <14>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000   39.182   39.182 <d22>/src/dmrgpy/kondospectrumtk/kondomethod.py:23(get_kondo_spectrum)
        1    0.000    0.000   39.182   39.182 <d22>/src/dmrgpy/kondospectrumtk/kondomethod.py:278(_get_kondo_spectrum_dmrg)
        1    0.005    0.005   38.956   38.956 <d22>/src/dmrgpy/kondospectrumtk/kondomethod.py:323(terms)
        1    0.519    0.519   19.256   19.256 <d22>/src/dmrgpy/kondospectrumtk/conductance.py:139(second_order_dIdV)
      126    0.611    0.005   18.537    0.147 <d22>/src/dmrgpy/kondospectrumtk/conductance.py:114(_theta_raw)
      126   17.926    0.142   17.926    0.142 <d22>/src/dmrgpy/kondospectrumtk/stepfunctions.py:58(Theta)
        1    0.002    0.002   11.102   11.102 <d22>/src/dmrgpy/kondospectrumtk/conductance.py:206(third_order_kondo_dIdV)
        2    0.257    0.129   11.061    5.531 <d22>/src/dmrgpy/kondospectrumtk/conductance.py:265(one_direction)
      248    1.013    0.004   10.266    0.041 <d22>/src/dmrgpy/kondospectrumtk/stepfunctions.py:214(__call__)
        1    0.000    0.000    8.556    8.556 <d22>/src/dmrgpy/kondospectrumtk/taskkondo.py:68(__init__)
        1    0.000    0.000    8.548    8.548 <d22>/src/dmrgpy/taskdmrg.py:183(session_spectrum)
        1    0.000    0.000    8.498    8.498 <d22>/src/dmrgpy/pyitensor/task.py:1338(spectrum_from_arrays)
        1    0.607    0.607    8.494    8.494 <d22>/src/dmrgpy/pyitensor/task.py:1253(projected_spectrum)
        1    0.018    0.018    5.945    5.945 <d22>/src/dmrgpy/pyitensor/task.py:1141(block_krylov)
```

**Reviewer (CONFIRMED)**: Reviewed with ed-reference-1 as one change. I could not refute it. It is not documented as intended: conductance.py and kondo.md state the third-order sums run over occupied initial states only and say nothing for the second-order term. The TASK spectrum sets p to a one-hot vector (taskkondo.py:135-136), so the dropped columns have p exactly 0 and contribute exactly 0.0. That is why numbers_change is false here. I reproduced the new code at dim 904 and 241 biases at T=0 on a mock spectrum (0.04 s, 0.04 GB); the old code was skipped there because it would not fit under my memory cap, so I did not reproduce the 6.5 GB peak or the hunter's end-to-end 63% figure. The prototype gives identical results on real ED chains with T=0, a degenerate ground state (2 occupied states), a complex Hamiltonian, and T>0 with all states occupied: max relative difference 0 to 3.4e-15 (scripts 01 and 02). Struck sub-claim: none. Not re-measured: the hunter's end-to-end 63% / 45% saved figures from d22.

Reviewer's own evidence:

```
01_second.out: '904 241 0.0 new 0.04s peak 0.04 GB' then 'skip old (mem)'. 02_occ.out: '5 0.0 0.0 dim 32 nocc 2 maxrel 1.5568031143512647e-16'; '6 3.0 0.0 dim 64 nocc 64 maxrel 3.4e-15'
```

Reviewer's measurement: mock dim 904, 241 biases, T=0: old not runnable under my cap, new 0.04 s / 0.04 GB; max relative difference 0 to 3.4e-15 on ED chains, exactly 0 on TASK by construction; end-to-end percentages are the hunter's, not re-measured.

Reviewer's scripts and outputs:

`01_second.py`:

```python
import numpy as np, time, tracemalloc
from dmrgpy.kondospectrumtk import conductance as cd
from dmrgpy import spinchain
def new(ks, eVs, T0=1.0, U=0.0, p_cut=cd.P_CUT, nbytes=4e6):
    eVs=np.asarray(eVs,float); occ=cd._occupied_states(ks,p_cut)
    M2=cd._spin_matrix_elements_squared(ks)[:,occ]
    eps=ks.e[:,None]-ks.e[None,occ]; w=ks.p[None,occ]*M2
    spin=np.empty(len(eVs)); ch=max(1,int(nbytes)//w.size)
    for s in range(0,len(eVs),ch):
        v=eVs[s:s+ch]
        spin[s:s+len(v)]=(np.einsum('fi,efi->e',w,cd._theta_raw(ks,v[:,None,None]-eps[None]))
          +np.einsum('fi,efi->e',w,cd._theta_raw(ks,-v[:,None,None]-eps[None])))
    return 2*np.pi*T0**2*(spin+4*U**2*np.ones_like(eVs))
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
for case in ["heis n=8 T=0","heis n=8 T=1K","heis n=8 T=30K","degenerate GS T=0 (field-free S=1/2 dimer-less n=3 odd)","complex H n=6 T=2K"]:
    n=8
    sc=spinchain.Spin_Chain([2]*n)
    h=0
    if "n=3" in case:
        sc=spinchain.Spin_Chain([2]*3); n=3
    for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    if "complex" in case:
        n=6; sc=spinchain.Spin_Chain([2]*n); h=0
        for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]+0.3*(sc.Sx[i]*sc.Sy[i+1]-sc.Sy[i]*sc.Sx[i+1])
    h=h+0.05*sc.Sz[0]
    sc.set_hamiltonian(h)
    T=0. if "T=0" in case else float(case.split("T=")[-1].split("K")[0]) if "K" in case else 2.
    ks=KondoSpectrum(sc,0,T)
    eV=np.linspace(-0.01,0.01,41)
    a=cd.second_order_dIdV(ks,eV,U=0.2); b=new(ks,eV,U=0.2)
    print(case,"dim",ks.dim,"nocc",len(cd._occupied_states(ks)),"maxrel",np.max(np.abs(a-b)/np.abs(a)))
# size/time on mock
class K: pass
for dim,nb,T in [(904,241,0.),(1024,121,1.)]:
    ks=K(); ks.dim=dim; rng=np.random.default_rng(0)
    ks.e=np.sort(rng.uniform(0,0.2,dim)); ks.e-=ks.e[0]; ks.T=T; ks.kB=8.617e-5
    for nm in ['Sx','Sy','Sz']: setattr(ks,nm,(rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim)))*0.05)
    if T==0: ks.p=np.zeros(dim); ks.p[0]=1
    else:
        p=np.exp(-ks.e/(ks.kB*T)); ks.p=p/p.sum()
    eV=np.linspace(-0.01,0.01,nb)
    for nm,f in [('new',new),('old',cd.second_order_dIdV)]:
        if nm=='old' and dim*dim*nb*8*6>9e9: print('skip old (mem)');continue
        tracemalloc.start(); t=time.time(); r=f(ks,eV); dt=time.time()-t
        print(dim,nb,T,nm,"%.2fs peak %.2f GB"%(dt,tracemalloc.get_traced_memory()[1]/1e9)); tracemalloc.stop()
        if nm=='new': rn=r
        else: print(" maxrel",np.max(np.abs(r-rn)/np.abs(r)))
```

`02_occ.py`:

```python
import numpy as np, importlib.util, sys
sys.argv=['x']
src=open('01_second.py').read().split("from dmrgpy.kondospectrumtk.edkondo")[0]
exec(src)
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
for n,T,B in [(6,3.,0.0),(7,3.,0.0),(6,10.,0.2),(5,0.,0.0)]:
    sc=spinchain.Spin_Chain([2]*n); h=0
    for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h=h+B*1e-3*sc.Sz[0]; sc.set_hamiltonian(h)
    ks=KondoSpectrum(sc,0,T); eV=np.linspace(-0.01,0.01,41)
    a=cd.second_order_dIdV(ks,eV,U=0.2); b=new(ks,eV,U=0.2)
    print(n,T,B,"dim",ks.dim,"nocc",len(cd._occupied_states(ks)),"maxrel",np.max(np.abs(a-b)/np.abs(a)))
```

`01_second.out`:

```
heis n=8 T=0 dim 256 nocc 1 maxrel 0.0
heis n=8 T=1K dim 256 nocc 1 maxrel 0.0
heis n=8 T=30K dim 256 nocc 1 maxrel 0.0
degenerate GS T=0 (field-free S=1/2 dimer-less n=3 odd) dim 8 nocc 1 maxrel 0.0
complex H n=6 T=2K dim 64 nocc 1 maxrel 0.0
904 241 0.0 new 0.04s peak 0.04 GB
skip old (mem)
1024 121 1.0 new 0.30s peak 0.10 GB
1024 121 1.0 old 20.30s peak 9.53 GB
 maxrel 7.738594311995933e-12
```

`02_occ.out`:

```
6 3.0 0.0 dim 64 nocc 64 maxrel 3.417445875158219e-15
7 3.0 0.0 dim 128 nocc 128 maxrel 3.4174600487122648e-15
6 10.0 0.2 dim 64 nocc 64 maxrel 1.5852746424178214e-15
5 0.0 0.0 dim 32 nocc 2 maxrel 1.5568031143512647e-16
```

**Suggested fix (hunter)**: In second_order_dIdV take occ = _occupied_states(ks), build M2[:, occ], eps = ks.e[:, None] - ks.e[occ][None, :] and weight = ks.p[occ][None, :] * M2[:, occ], and evaluate the two Theta sums in chunks of about 32 biases so the temporary stays small (the prototype is the function `new` in 09_second_order_occ.py). Numbers: identical at T=0, and at T>0 within 1e-12 relative, the P_CUT the third-order terms use; measured 7e-16 here. Same patch helps mode='ED' on large Fock spaces. Also leave third_order_kondo_dIdV (11.1 s, 28% at T>0, FBuilder interpolation 5.5 s plus F0 3.7 s) for the conductance lens.

**Fix (reviewer)**: One patch in conductance.py shared with ed-reference-1: occ=_occupied_states(ks, p_cut); M2[:,occ]; eps=ks.e[:,None]-ks.e[None,occ]; weight=ks.p[None,occ]*M2; evaluate both Theta sums in bias chunks of max(1, 4e6//weight.size). Add a p_cut=P_CUT keyword like the third-order terms so p_cut=0 reproduces the full sum. Pin it in a test comparing against the old full-sum formula at p_cut=0 (T=0 degenerate, T>0) with the chunk forced to 1 and 3, plus a TASK-shaped test with p one-hot and a tracemalloc cap on the temporaries. The hunter fix is right.

### 2. At T>0 the unconverged-pole warning promises a change of 4.1e-3 where nkry=100 and 300 differ by 6.7e-8

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `task-route` &middot; candidate `task-route-2`

**Status**: FIXED in `d098cd0` on 2026-11-12, in the reviewer's form, which the author picked. At T>0 `convergence()` scales the weight of each loose pole in "step" by min(1, d/(5.4 kT)), with d = min(r, r^2/g) and g the distance to the nearest other Ritz energy; the third-nearest-pole filter stays, and T=0 is unchanged. Measured against `nkry=300` on one state at `nkry=100`, on d22 (`<scratch>/fix/taskguard/01_estimate.py`): gapped 20-site chain, 4.1e-3 to 3.0e-5 at 0.1 K (true change 6.5e-8) and to 1.1e-4 at 0.05 K (7.0e-7), so it no longer warns; uniform 24-site chain, 3.7e-2 to 1.9e-2 at 0.03 K (2.2e-3) and to 1.5e-2 at 0.05 K (8.7e-4), so it still warns. In every case measured the new size stays between the true change and the old size. Pinned by `tests/test_audit_2026_11_task.py` on a 12-site version of the gapped chain (`nkry=60` at 0.1 K: old 2.9e-3, new 3.1e-4, true 2.4e-6), with the old rule kept verbatim.

The claim as reviewed: At T>0 the TASK unconverged-pole warning is a loud false alarm on a gapped disordered chain: at T=0.1 K (20 sites, maxm=40, site 3, order=2, nkry=100) it promises a change of up to 4.1e-3 of the largest value while nkry=100 against nkry=300 differs by 6.7e-8 and two different solves differ by 3.7e-8, because a pole with residual r moves by about r^2/gap, which the criterion does not use. The T=0 half of the hunter's claim (overstatement by 2 to 4 orders, 5.8e-5 against 6.4e-3) is struck: it is an artifact of a uniform bias grid.

**Where**: src/dmrgpy/kondospectrumtk/taskkondo.py:208-262 (convergence: loose = residual > CONVERGED_RESIDUAL*width; step = max weight of loose poles / total; UNCONVERGED_STEP = 1e-3; THERMAL_WIDTH rule)

The warning's step is the largest weight of any pole with residual > 1e-8 of the Ritz bandwidth, divided by the inelastic weight, i.e. it assumes a flagged pole sits at a wrong energy by at least its own spacing. A Ritz energy error is of order residual^2 / gap, so a pole with residual 1e-4 of the width is far more accurate than a step at an arbitrary position; the criterion does not use the quadratic relation or the spacing. Chain: 20 sites, disordered XXZ (Jz=1.5J, J in 0.5 to 1.5e-3 eV, staggered field 6e-4), python, maxm=40, site 3, order=2 (the quantity the docstring names), 6001 biases so that no step falls between grid points. T=0: nkry=100 flags 113 of 160 poles, weight 0.027, step 6.4e-3 and warns; the curve differs from nkry=300 by 5.81e-5 of the peak. T=0.05 K: 91 of 160 flagged, step 4.1e-3, warns; difference 5.94e-7. T=0.1 K: 32 flagged, step 4.1e-3, warns; difference 5.09e-8. With the 241-point grid the differences were 3.8e-6, 2.4e-7, 3.4e-8. The T>0 spacing rule (third nearest pole > 5.4 kT) therefore keeps warning at T=0.1 K where the curve is converged to 5e-8, while at T=0.2 K and 0.3 K (no warning) the nkry=100 curve is within 5e-9 and 8e-10. It is a loud false alarm, not a silent wrong number, and the flagged-pole population is dominated by high-energy poles of small weight. It survives because the thresholds were fitted on one gapless uniform chain (24 sites, where the author measured 3e-3 at T=0), where the residual criterion happens to be of the right size. The docstring's 'relative to its largest value' is also normalized by the inelastic weight, not by the largest dI/dV (elastic part included).

The hunter's claim, before review: TaskKondoSpectrum.warn_unconverged/convergence overstate the error of unconverged poles by two to four orders of magnitude on a gapped disordered chain: the warning promises a change of up to 6.4e-3 (T=0) or 4.1e-3 (T=0.05 K) of the largest value, while nkry=100 against nkry=300 on a 6001-point bias grid differs by 5.8e-5 and 5.9e-7 of the peak (second order).

Decisive output, quoted by the hunter:

```
T 0.0 nkry 100 warned True {'poles': 160, 'unconverged': 113, 'weight': 0.0274, 'step': 0.006368}
  T 0.0 max|d100-d300|/max|d300| = 5.81e-05
T 0.05 nkry 100 warned True {... 'unconverged': 91, 'step': 0.004123}
  T 0.05 max|d100-d300|/max|d300| = 5.94e-07
T 0.1 nkry 100 warned True {... 'unconverged': 32, 'step': 0.004123}
  T 0.1 max|d100-d300|/max|d300| = 5.09e-08
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`11_dense_grid_order2.py`:

```python
# does the T>0 warning rule (third nearest pole > 5.4 kT) transfer to a chain it was not fitted on?
# 20-site disordered XXZ chain, gapped by a staggered field so the sample stays cold; nkry 100 vs 300
import sys, warnings, numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
from dmrgpy.kondospectrumtk import conductance as C
def new(ks,eVs,T0=1.0,U=0.0):
    eVs=np.asarray(eVs,dtype=float); occ=C._occupied_states(ks)
    M2=C._spin_matrix_elements_squared(ks)[:,occ]; eps=ks.e[:,None]-ks.e[occ][None,:]
    weight=ks.p[occ][None,:]*M2; spin=np.empty(len(eVs))
    for a in range(0,len(eVs),32):
        e=eVs[a:a+32]
        spin[a:a+32]=(np.einsum('fi,efi->e',weight,C._theta_raw(ks,e[:,None,None]-eps[None]))+np.einsum('fi,efi->e',weight,C._theta_raw(ks,-e[:,None,None]-eps[None])))
    return 2*np.pi*T0**2*(spin+4*U**2*np.ones_like(eVs))
C.second_order_dIdV=new   # only to fit the dense grid in memory; identical numbers (5e-16)
n=20; rng=np.random.default_rng(5)
sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
J=1e-3*rng.uniform(0.5,1.5,n-1)
h=0
for i in range(n-1): h=h+J[i]*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1])
for i in range(n): h=h+6e-4*(-1)**i*sc.Sz[i]+1e-4*rng.uniform(-1,1)*sc.Sz[i]
sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=40,14,1e-10
print("E0",sc.gs_energy(),flush=True)
eV=np.linspace(-6e-3,6e-3,6001)
kw=dict(site=3,Jrho_s=-0.1,U=0.,omega0=20e-3,Gamma0=5e-5,order=2)
curves={}
for T in (0.,0.05,0.1):
  for nkry in (100,300):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try: _,d=sc.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=nkry,T=T,**kw)
        except ValueError as ex: print("T",T,"nkry",nkry,"RAISED",str(ex)[:120],flush=True); d=None
    curves[T,nkry]=d
    ks=TaskKondoSpectrum(sc,3,nkry=nkry,T=T) if d is not None else None
    c=ks.convergence(6e-3) if ks else None
    if ks is not None and nkry==100:
        w_=np.abs(ks.Sx[:,0])**2+np.abs(ks.Sy[:,0])**2+np.abs(ks.Sz[:,0])**2; width=np.max(np.abs(ks.e)); m=(ks.residual>1e-8*width)&(w_>1e-14*w_.sum())&(ks.e<6e-3)
        print('   flagged poles',int(m.sum()),'residual/width: min %.1e max %.1e'%(ks.residual[m].min()/width,ks.residual[m].max()/width),'width',width,'max res^2/gap(1e-5)=%.1e eV'%(ks.residual[m].max()**2/1e-5))
    print("T",T,"nkry",nkry,"warned",any("not converged" in str(x.message) for x in w),c,flush=True)
  if curves[T,100] is not None and curves[T,300] is not None:
    d1,d3=curves[T,100],curves[T,300]
    print("  T",T,"max|d100-d300|/max|d300| = %.2e"%(np.max(np.abs(d1-d3))/np.max(np.abs(d3))),flush=True)
```

`06_thermal_width.py`:

```python
# does the T>0 warning rule (third nearest pole > 5.4 kT) transfer to a chain it was not fitted on?
# 20-site disordered XXZ chain, gapped by a staggered field so the sample stays cold; nkry 100 vs 300
import sys, warnings, numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
n=20; rng=np.random.default_rng(5)
sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
J=1e-3*rng.uniform(0.5,1.5,n-1)
h=0
for i in range(n-1): h=h+J[i]*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1])
for i in range(n): h=h+6e-4*(-1)**i*sc.Sz[i]+1e-4*rng.uniform(-1,1)*sc.Sz[i]
sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=40,14,1e-10
print("E0",sc.gs_energy(),flush=True)
eV=np.linspace(-6e-3,6e-3,241)
kw=dict(site=3,Jrho_s=-0.1,U=0.3,omega0=20e-3,Gamma0=5e-5,order=3)
curves={}
for T in (0.05,0.1,0.2,0.3):
  for nkry in (100,300):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try: _,d=sc.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=nkry,T=T,**kw)
        except ValueError as ex: print("T",T,"nkry",nkry,"RAISED",str(ex)[:120],flush=True); d=None
    curves[T,nkry]=d
    ks=TaskKondoSpectrum(sc,3,nkry=nkry,T=T) if d is not None else None
    c=ks.convergence(6e-3) if ks else None
    print("T",T,"nkry",nkry,"warned",any("not converged" in str(x.message) for x in w),c,flush=True)
  if curves[T,100] is not None and curves[T,300] is not None:
    d1,d3=curves[T,100],curves[T,300]
    print("  T",T,"max|d100-d300|/max|d300| = %.2e"%(np.max(np.abs(d1-d3))/np.max(np.abs(d3))),flush=True)
```

`07_residual_vs_curve.py`:

```python
# does the T>0 warning rule (third nearest pole > 5.4 kT) transfer to a chain it was not fitted on?
# 20-site disordered XXZ chain, gapped by a staggered field so the sample stays cold; nkry 100 vs 300
import sys, warnings, numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
n=20; rng=np.random.default_rng(5)
sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
J=1e-3*rng.uniform(0.5,1.5,n-1)
h=0
for i in range(n-1): h=h+J[i]*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1])
for i in range(n): h=h+6e-4*(-1)**i*sc.Sz[i]+1e-4*rng.uniform(-1,1)*sc.Sz[i]
sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=40,14,1e-10
print("E0",sc.gs_energy(),flush=True)
eV=np.linspace(-6e-3,6e-3,241)
kw=dict(site=3,Jrho_s=-0.1,U=0.3,omega0=20e-3,Gamma0=5e-5,order=3)
curves={}
for T in (0.,0.05):
  for nkry in (100,300):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try: _,d=sc.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=nkry,T=T,**kw)
        except ValueError as ex: print("T",T,"nkry",nkry,"RAISED",str(ex)[:120],flush=True); d=None
    curves[T,nkry]=d
    ks=TaskKondoSpectrum(sc,3,nkry=nkry,T=T) if d is not None else None
    c=ks.convergence(6e-3) if ks else None
    if ks is not None and nkry==100:
        w_=np.abs(ks.Sx[:,0])**2+np.abs(ks.Sy[:,0])**2+np.abs(ks.Sz[:,0])**2; width=np.max(np.abs(ks.e)); m=(ks.residual>1e-8*width)&(w_>1e-14*w_.sum())&(ks.e<6e-3)
        print('   flagged poles',int(m.sum()),'residual/width: min %.1e max %.1e'%(ks.residual[m].min()/width,ks.residual[m].max()/width),'width',width,'max res^2/gap(1e-5)=%.1e eV'%(ks.residual[m].max()**2/1e-5))
    print("T",T,"nkry",nkry,"warned",any("not converged" in str(x.message) for x in w),c,flush=True)
  if curves[T,100] is not None and curves[T,300] is not None:
    d1,d3=curves[T,100],curves[T,300]
    print("  T",T,"max|d100-d300|/max|d300| = %.2e"%(np.max(np.abs(d1-d3))/np.max(np.abs(d3))),flush=True)
```

Observed:

`11_dense_grid_order2.out`:

```
  from dmrgpy import spinchain
E0 -0.014975352636028269
   flagged poles 113 residual/width: min 1.2e-08 max 8.5e-02 width 0.009948906493900164 max res^2/gap(1e-5)=7.1e-02 eV
T 0.0 nkry 100 warned True {'poles': 160, 'unconverged': 113, 'weight': 0.027362352208516495, 'step': 0.006368439222184336}
T 0.0 nkry 300 warned False {'poles': 308, 'unconverged': 55, 'weight': 2.1957213286552133e-06, 'step': 3.249046990025471e-07}
  T 0.0 max|d100-d300|/max|d300| = 5.81e-05
   flagged poles 113 residual/width: min 1.2e-08 max 8.5e-02 width 0.009948906493900164 max res^2/gap(1e-5)=7.1e-02 eV
T 0.05 nkry 100 warned True {'poles': 160, 'unconverged': 91, 'weight': 0.017421532286787233, 'step': 0.004122583097322781}
T 0.05 nkry 300 warned False {'poles': 308, 'unconverged': 6, 'weight': 3.0066471371499707e-07, 'step': 1.665859861609411e-07}
  T 0.05 max|d100-d300|/max|d300| = 5.94e-07
   flagged poles 113 residual/width: min 1.2e-08 max 8.5e-02 width 0.009948906493900164 max res^2/gap(1e-5)=7.1e-02 eV
T 0.1 nkry 100 warned True {'poles': 160, 'unconverged': 32, 'weight': 0.01371268301475021, 'step': 0.004122583097322781}
T 0.1 nkry 300 warned False {'poles': 308, 'unconverged': 0, 'weight': 0.0, 'step': 0.0}
  T 0.1 max|d100-d300|/max|d300| = 5.09e-08
```

`06_thermal_width.out`:

```
  from dmrgpy import spinchain
E0 -0.014975352636028321
T 0.05 nkry 100 warned True {'poles': 159, 'unconverged': 91, 'weight': 0.01746228060957514, 'step': 0.004122711237482384}
T 0.05 nkry 300 warned False {'poles': 308, 'unconverged': 6, 'weight': 3.0066471371171635e-07, 'step': 1.6658598616116237e-07}
  T 0.05 max|d100-d300|/max|d300| = 2.51e-07
T 0.1 nkry 100 warned True {'poles': 159, 'unconverged': 36, 'weight': 0.013699619013533955, 'step': 0.004122711237482384}
T 0.1 nkry 300 warned False {'poles': 308, 'unconverged': 0, 'weight': 0.0, 'step': 0.0}
  T 0.1 max|d100-d300|/max|d300| = 3.42e-08
T 0.2 nkry 100 warned False {'poles': 159, 'unconverged': 2, 'weight': 2.2729480374743886e-07, 'step': 1.6643097529330925e-07}
T 0.2 nkry 300 warned False {'poles': 308, 'unconverged': 0, 'weight': 0.0, 'step': 0.0}
  T 0.2 max|d100-d300|/max|d300| = 4.99e-09
T 0.3 nkry 100 warned False {'poles': 159, 'unconverged': 0, 'weight': 0.0, 'step': 0.0}
T 0.3 nkry 300 warned False {'poles': 308, 'unconverged': 0, 'weight': 0.0, 'step': 0.0}
  T 0.3 max|d100-d300|/max|d300| = 8.47e-10
```

`07_residual_vs_curve.out`:

```
  from dmrgpy import spinchain
E0 -0.014975352636028283
   flagged poles 111 residual/width: min 1.4e-08 max 7.4e-02 width 0.009948906490797504 max res^2/gap(1e-5)=5.4e-02 eV
T 0.0 nkry 100 warned True {'poles': 159, 'unconverged': 111, 'weight': 0.027362331744125685, 'step': 0.006368493997328501}
T 0.0 nkry 300 warned False {'poles': 309, 'unconverged': 54, 'weight': 2.3010773159177406e-06, 'step': 3.249046990337223e-07}
  T 0.0 max|d100-d300|/max|d300| = 3.80e-06
   flagged poles 111 residual/width: min 1.4e-08 max 7.4e-02 width 0.009948906490797504 max res^2/gap(1e-5)=5.4e-02 eV
T 0.05 nkry 100 warned True {'poles': 159, 'unconverged': 89, 'weight': 0.01742062254154669, 'step': 0.004122523882453829}
T 0.05 nkry 300 warned False {'poles': 309, 'unconverged': 5, 'weight': 3.006645241429638e-07, 'step': 1.6658598616107402e-07}
  T 0.05 max|d100-d300|/max|d300| = 2.42e-07
```

**Reviewer (CONFIRMED (narrowed))**: I rebuilt the hunter's chain and measured on one state with TaskKondoSpectrum at nkry=100 and 300. The hunter's number at T=0 is the maximum over a uniform 6001-point grid with spacing 2 microeV. An unconverged pole at T=0 is a sharp Theta step, and a displacement much smaller than the grid spacing is invisible on the grid although its sup-norm effect is the full step height. Evaluating both curves at the biases between and at the pole energies of the two runs gives 2.18e-1 of the peak for nkry=100 against 300, and 3.41e-1 for two DMRG solves of the same chain at nkry=100 (nsweeps=14, maxm=40 against nsweeps=18, maxm=41, energies equal to 1e-17). The warned 6.4e-3 is therefore below the real sup-norm change, not above it, and the 'up to' wording is true at T=0 (and it says 'at the biases next to these poles'). On the uniform grid the same two comparisons give 1.15e-4 and 5.39e-5. At T>0 the thermal step is smooth with width 5.4 kT, and the hunter's claim survives: at T=0.1 K both the pole-resolved biases and the uniform grid give 6.7e-8 (nkry) and 3.5e-8 / 3.7e-8 (solve to solve), against a warned step of 4.1e-3. The criterion flags 34 of 159 poles by residual alone. So the warning is a bound that is loose by about 5 orders in this regime, from the second-order displacement, which the author's THERMAL_WIDTH fit (one 24-site gapless chain) does not capture. It is a warning, not a wrong number; no returned value changes. Not documented as a limitation: kondo.md states THERMAL_WIDTH is fitted on one chain, and the guide says 'once the thermal width is larger than the spacing... a smaller nkry is enough', which is the regime where the warning stays on. Sub-claim not checked: the docstring normalization (inelastic weight against largest dI/dV).

Reviewer's own evidence:

```
01_supnorm.out: 'T 0.0 conv100 {... unconverged 111, step 0.006368}'; 'T 0.0 points 4720 max|d100-d300|/max|d300| on pole-resolved biases = 2.18e-01'; 'uniform 6001 grid = 1.15e-04'; 'T 0.1 conv100 {... unconverged 34, step 0.004123}'; 'T 0.1 ... pole-resolved = 6.71e-08, uniform 6001 grid = 6.71e-08'. 02_solves.out: 'T 0.0 resolved solve-to-solve 3.41e-01 ... T 0.0 grid 5.39e-05', 'T 0.1 resolved 3.54e-08, grid 3.67e-08, warned step 0.00412'.
```

Reviewer's scripts and outputs:

`01_supnorm.py`:

```python
# Is the warned step a real sup-norm bound? Evaluate the second-order dI/dV of nkry=100 and 300 on one state
# at the biases BETWEEN the pole energies of the two runs (where a displaced step shows), not on a uniform grid.
import numpy as np, warnings
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
from dmrgpy.kondospectrumtk import conductance as C
def d2(ks,eVs):
    occ=C._occupied_states(ks); M2=C._spin_matrix_elements_squared(ks)[:,occ]
    eps=ks.e[:,None]-ks.e[occ][None,:]; wt=ks.p[occ][None,:]*M2; out=np.empty(len(eVs))
    for a in range(0,len(eVs),16):
        e=eVs[a:a+16]
        out[a:a+16]=(np.einsum('fi,efi->e',wt,C._theta_raw(ks,e[:,None,None]-eps[None]))+np.einsum('fi,efi->e',wt,C._theta_raw(ks,-e[:,None,None]-eps[None])))
    return 2*np.pi*out
n=20; rng=np.random.default_rng(5)
sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
J=1e-3*rng.uniform(0.5,1.5,n-1); h=0
for i in range(n-1): h=h+J[i]*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1])
for i in range(n): h=h+6e-4*(-1)**i*sc.Sz[i]+1e-4*rng.uniform(-1,1)*sc.Sz[i]
sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=40,14,1e-10
print("E0",sc.gs_energy(),flush=True)
for T in (0.,0.1):
    k1=TaskKondoSpectrum(sc,3,nkry=100,T=T); k3=TaskKondoSpectrum(sc,3,nkry=300,T=T)
    print("T",T,"conv100",k1.convergence(6e-3))
    E=np.sort(np.concatenate([k1.e,k3.e])); E=E[(E>0)&(E<6e-3)]
    kT=k1.kB*T
    pts=np.concatenate([E,(E[1:]+E[:-1])/2,E+1e-9,E-1e-9])
    pts=np.concatenate([pts,-pts]); pts=np.unique(pts)
    a=d2(k1,pts); b=d2(k3,pts)
    print("  points",len(pts),"max|d100-d300|/max|d300| on pole-resolved biases = %.2e"%(np.max(np.abs(a-b))/np.max(np.abs(b))))
    g=np.linspace(-6e-3,6e-3,6001); a=d2(k1,g); b=d2(k3,g)
    print("  uniform 6001 grid = %.2e"%(np.max(np.abs(a-b))/np.max(np.abs(b))))
```

`02_solves.py`:

```python
# solve-to-solve variation at fixed nkry=100 (the other half of the warning's claim), two different DMRG solves of one chain
import numpy as np
exec(open('01_supnorm.py').read().split("n=20;")[0])
def build(ns,maxm):
    n=20; rng=np.random.default_rng(5)
    sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
    J=1e-3*rng.uniform(0.5,1.5,n-1); h=0
    for i in range(n-1): h=h+J[i]*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1])
    for i in range(n): h=h+6e-4*(-1)**i*sc.Sz[i]+1e-4*rng.uniform(-1,1)*sc.Sz[i]
    sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=maxm,ns,1e-10
    sc.gs_energy(); return sc
A=build(14,40); B=build(18,41)
print("E0",A.gs_energy(),B.gs_energy())
g=np.linspace(-6e-3,6e-3,6001)
for T in (0.,0.1):
    ka=TaskKondoSpectrum(A,3,nkry=100,T=T); kb=TaskKondoSpectrum(B,3,nkry=100,T=T)
    E=np.sort(np.concatenate([ka.e,kb.e])); E=E[(E>0)&(E<6e-3)]
    pts=np.unique(np.concatenate([E,(E[1:]+E[:-1])/2])); pts=np.concatenate([pts,-pts])
    for name,x in (("resolved",pts),("grid",g)):
        a=d2(ka,x); b=d2(kb,x)
        print("T",T,name,"solve-to-solve %.2e"%(np.max(np.abs(a-b))/np.max(np.abs(b))),"warned step",ka.convergence(6e-3)['step'])
```

`01_supnorm.out`:

```
E0 -0.014975352636028245
T 0.0 conv100 {'poles': 159, 'unconverged': 111, 'weight': 0.027362344082563345, 'step': 0.006368541112196055}
  points 4720 max|d100-d300|/max|d300| on pole-resolved biases = 2.18e-01
  uniform 6001 grid = 1.15e-04
T 0.1 conv100 {'poles': 159, 'unconverged': 34, 'weight': 0.01367638606865277, 'step': 0.004123029395141609}
  points 4720 max|d100-d300|/max|d300| on pole-resolved biases = 6.71e-08
  uniform 6001 grid = 6.71e-08
```

`02_solves.out`:

```
E0 -0.0149753526360283 -0.014975352636028248
T 0.0 resolved solve-to-solve 3.41e-01 warned step 0.006368541471782321
T 0.0 grid solve-to-solve 5.39e-05 warned step 0.006368541471782321
T 0.1 resolved solve-to-solve 3.54e-08 warned step 0.004123140269227178
T 0.1 grid solve-to-solve 3.67e-08 warned step 0.004123140269227178
```

**Suggested fix (hunter)**: Estimate the pole displacement instead of flagging by residual alone: delta_e ~ residual^2 / (distance to the nearest other Ritz energy), capped at that distance, and let a pole count in 'step' only when delta_e exceeds the resolution that matters (the thermal width 5.4 kT at T>0, or Gamma0 at T=0). Normalize the step by the largest second-order dI/dV or by elastic plus inelastic weight, as the docstring says. Refit THERMAL_WIDTH and UNCONVERGED_STEP on at least two chains (gapless and gapped). Returned dI/dV numbers do not change; only which calls warn.

**Fix (reviewer)**: Change only TaskKondoSpectrum.convergence in taskkondo.py and only for T>0; leave T=0 untouched, since there the step is the true sup-norm of one pole and the warning understates the total change. For T>0 estimate each loose pole's displacement as d_k = min(res_k, res_k^2 / dist_k), with dist_k the distance to its nearest other Ritz energy, and count its weight in 'step' scaled by min(1, d_k/(THERMAL_WIDTH*kT)), which is the height of the thermal step difference. The hunter's fix is right in direction; I would not refit UNCONVERGED_STEP or normalize by the largest dI/dV in the same change. Pin it in tests/test_kondo_task_convergence.py with the gapped disordered chain: no warning at T=0.1 K and nkry=100, and the warning unchanged at T=0 and on the 24-site uniform chain. Returned dI/dV does not change. Alternative with no new estimate: say in the warning text that at T>0 the bound is loose when the poles are many and close.

### 3. The same sum on mode="ED": 9.5 GB at 1024 states and 121 biases

`performance` &middot; severity **MEDIUM** &middot; CONFIRMED (narrowed) &middot; lens `ed-reference` &middot; candidate `ed-reference-1`

**Status**: FIXED in `49c23ca`, the change of finding 1. `mode="ED"`, 10 spins 1/2, 41 biases, U=0.25, order 3, on d22: T=0 2.93 s to 2.45 s (5.0e-16 apart), T=1 K 7.01 s to 4.32 s; `second_order_dIdV` alone at 1 K 2.70 s and 2.90 GB to 0.03 s and 0.05 GB. NUMBERS CHANGE, on `mode="ED"` at T>0 only, through `P_CUT`: 9.8e-15 relative on that chain at 1 K, 5e-16 at T=0; `p_cut=0` gives the old sum back to 1e-13. `tests/test_ysr_bcs_second_order.py`, which pins `ysrtk/bcsperturbative.py` to this function at 1e-13 of the maximum, passes. Pinned by the tests of finding 1.

The claim as reviewed: conductance.second_order_dIdV builds (biases x dim x dim) arrays over all initial states, so ED at large dim and many biases needs far more memory than the occupied states warrant (mock dim 1024, 121 biases, T=1 K: old 20.3 s and 9.5 GB, restricted and chunked 0.30 s and 0.10 GB). Chunking alone is exact. Restricting to occupied states applies P_CUT (1e-12 relative) to the second-order term, which moves numbers by up to about 1e-11 relative (7.7e-12 measured on a spread-out p), so beyond roundoff but far below the Theta/F accuracy. numbers_change is true for ED at T>0, false at T=0 and on TASK.

**Where**: src/dmrgpy/kondospectrumtk/conductance.py:139-160 (second_order_dIdV: x_ts, x_st, _theta_raw temporaries, einsum over all f,i); the third-order terms next to it already chunk (conductance.py:240) and already restrict to occupied i (_occupied_states, P_CUT)

second_order_dIdV forms x_ts and x_st of shape (nbias, dim, dim) and evaluates Theta (sinh, cosh and three masked copies) on both, for every initial state i, although p_i is zero beyond n_occ states (4 of 1024 at 1 K, exactly the sets the third-order terms already select with _occupied_states). Memory therefore grows as nbias*dim^2 with a large constant: 3.23 GB at 41 biases, 8.71 GB at 121, 16.92 GB at 241 (n=10, T=1 K, fresh process each, script 10). On a 16 GB machine a 241-point scan of a 10-spin chain does not fit, while the third-order terms of the same call stay at a few hundred MB. Prototype (script proto.py, monkeypatch): restrict to the occupied columns and loop over bias chunks of 4e6 entries. Result: same sums, peak RSS 0.79 GB at every bias count; time 41 biases T=1 K median of 3 fresh processes (script 18, whole call including eigh) 7.56 to 4.93 s (34.7%), T=0 3.47 to 2.99 s (13.9%); 241 biases 55.2 to 21.7 s (single run). Numbers: not bit-identical, because P_CUT (1e-12 of the largest weight) is now also applied to the second-order term: max relative difference to the old public result 5.0e-16 at T=0 and 9.8e-15 at T=1 K (3.3e-13 of the peak in the 3 K run of script 07), so numbers_change is true only at that level. Why tests pass: they use 1 to 4 spin chains and 21 to 61 biases, where the arrays are kilobytes; the memory is a cost, not a wrong number.

The hunter's claim, before review: conductance.second_order_dIdV builds three unchunked (biases x dim x dim) arrays over ALL initial states, so a 241-bias ED spectrum of 10 spins 1/2 at T=1 K peaks at 16.92 GB (0.79 GB chunked and restricted to occupied states) and takes 55.2 s against 21.7 s; at 41 biases the public call saves 35% at T=1 K (7.56 to 4.93 s) and 14% at T=0 (3.47 to 2.99 s).

Decisive output, quoted by the hunter:

```
241 biases old      time 55.2 s, peak RSS 16.92 GB, sum 1243.262237
  241 biases second   time 21.7 s, peak RSS 0.79 GB, sum 1243.262237
  n=10 T=1 K 41 biases  old       median  7.56 s  (7.57 7.54 7.56)  saves  0.0%  max rel diff to old 0.0e+00
  n=10 T=1 K 41 biases  second    median  4.93 s  (4.94 4.93 4.92)  saves 34.7%  max rel diff to old 9.8e-15
  n=10 T=0 K 41 biases  second    median  2.99 s  (2.99 2.99 2.99)  saves 13.9%  max rel diff to old 5.0e-16
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`proto.py`:

```python
# prototypes (monkeypatch only), nothing in the repo is edited
import numpy as np
from scipy import sparse
from dmrgpy.kondospectrumtk import stepfunctions as sf, conductance as cd, edkondo
from dmrgpy.algebra import algebra

def Theta_new(x):
    """Theta(x) = (1-t-xt)/(1-t)^2, t=exp(-|x|), mirrored for x<0: one expm1, no sinh/cosh, no masks"""
    x = np.asarray(x, dtype=float)
    a = np.abs(x)
    em1 = np.expm1(-a)
    t = 1.0 + em1
    with np.errstate(divide="ignore", invalid="ignore"):
        pos = (-em1 - a*t)/(em1*em1)
    pos = np.where(a < 1e-4, 0.5 + a/6.0 - a**3/360.0, pos)
    return np.where(x >= 0, pos, 1.0 - pos)

def F_call_new(self, x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    if self._F_spline is None: return self._direct(x)
    inside = (x >= self._F_min) & (x <= self._F_max)
    if inside.all(): return self._F_spline(x)          # F0 only where the table does not reach
    out = np.empty_like(x)
    out[inside] = self._F_spline(x[inside])
    xo = x[~inside]
    out[~inside] = sf.F0(xo, omega0=self.omega0, Gamma0=self.Gamma0)
    return out

def second_order_new(ks, eVs, T0=1.0, U=0.0, p_cut=cd.P_CUT, nbytes=4e6):
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    M2 = cd._spin_matrix_elements_squared(ks)[:, occ]          # [f,i in occ]
    eps = ks.e[:, None] - ks.e[None, occ]
    weight = ks.p[None, occ]*M2
    spin = np.empty(len(eVs))
    chunk = max(1, int(nbytes)//weight.size)
    for s0 in range(0, len(eVs), chunk):
        v = eVs[s0:s0+chunk]
        spin[s0:s0+len(v)] = (np.einsum('fi,efi->e', weight, cd._theta_raw(ks, v[:, None, None] - eps[None]))
                              + np.einsum('fi,efi->e', weight, cd._theta_raw(ks, -v[:, None, None] - eps[None])))
    return 2*np.pi*T0**2*(spin + 4*(U**2)*np.ones_like(eVs))

def eigenbasis_new(self, name):
    edobj = self.chain.get_ED_obj()
    op = getattr(self.chain, name)[self.site]
    mat = edobj.get_operator(op)
    vs = self.vs
    if sparse.issparse(mat):   # sparse @ dense first (O(dim^2)), then ONE dense product instead of two
        mat = mat.tocsr().astype(complex)
        return np.conjugate(vs.T) @ (mat @ vs)
    mat = np.array(algebra.todense(mat), dtype=complex)
    return np.conjugate(vs.T)@mat@vs

def F0_new(x, omega0=20e-3, Gamma0=5e-6):
    x = np.asarray(x, dtype=float)
    ax = np.abs(x)
    return np.log((omega0 + ax)/np.sqrt(x*x + Gamma0*Gamma0))

_old = dict(F0=sf.F0, Theta=sf.Theta, call=sf.FBuilder.__call__, second=cd.second_order_dIdV, eig=edkondo.KondoSpectrum._eigenbasis_spin_operator)
def enable(which):
    if "theta" in which: sf.Theta = Theta_new; cd.Theta = Theta_new
    if "F" in which: sf.FBuilder.__call__ = F_call_new
    if "F0" in which: sf.F0 = F0_new; cd.F0 = F0_new
    if "second" in which: cd.second_order_dIdV = second_order_new
    if "eig" in which: edkondo.KondoSpectrum._eigenbasis_spin_operator = eigenbasis_new
def disable():
    sf.F0 = cd.F0 = _old["F0"]; sf.Theta = cd.Theta = _old["Theta"]; sf.FBuilder.__call__ = _old["call"]
    cd.second_order_dIdV = _old["second"]; edkondo.KondoSpectrum._eigenbasis_spin_operator = _old["eig"]

def third_order_kondo_blocked(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6, Fb=None, p_cut=cd.P_CUT,
                              block_bytes=2.0e8, w_bytes=3.0e8):
    """same sums as conductance.third_order_kondo_dIdV; the coefficient tensor is built for a block of initial
    states at a time (block_bytes), and the exchange diagram is accumulated into W[e,f,m] for a group of biases"""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    dim = ks.dim
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    e_fm = ks.e[:, None] - ks.e[None, :]
    nB = max(1, int(block_bytes//(16*dim*dim)))
    nG = max(1, int(w_bytes//(8*dim*dim)))          # biases per group (W accumulator)
    fchunk = max(1, int(4e6)//(dim*dim))
    def one_direction(v):
        out = np.zeros(len(v))
        for g0 in range(0, len(v), nG):
            vg = v[g0:g0+nG]
            W = np.zeros((len(vg), dim, dim))
            for b0 in range(0, len(occ), nB):
                ob = occ[b0:b0+nB]; p = ks.p[ob]
                coeff = np.imag(cd._triple_product_coefficients(ks, ob))/2.
                eps = ks.e[None, :] - ks.e[ob, None]
                for s0 in range(0, len(vg), max(1, int(4e6)//(dim*len(ob)))):
                    vc = vg[s0:s0+max(1, int(4e6)//(dim*len(ob)))]
                    Th = cd._theta_raw(ks, vc[:, None, None] - eps[None])
                    Fd = Fcall((vc[:, None, None] - eps[None]).ravel()).reshape(len(vc), len(ob), dim)
                    out[g0+s0:g0+s0+len(vc)] += np.einsum('i,ifm,eif,eim->e', p, coeff, Th, Fd, optimize=True)
                    W[s0:s0+len(vc)] += np.einsum('i,ifm,eif->efm', p, coeff, Th, optimize=True)
            for s0 in range(0, len(vg), fchunk):
                vc = vg[s0:s0+fchunk]
                Fx = Fcall((vc[:, None, None] - e_fm[None]).ravel()).reshape(len(vc), dim, dim)
                out[g0+s0:g0+s0+len(vc)] += np.einsum('efm,efm->e', W[s0:s0+len(vc)], Fx)
        return out
    total = one_direction(eVs) + one_direction(-eVs)
    return 4*np.pi*T0**2*Jrho_s*total

_old["third"] = cd.third_order_kondo_dIdV
_enable_old = enable
def enable(which):
    _enable_old(which)
    if "third" in which: cd.third_order_kondo_dIdV = third_order_kondo_blocked
_disable_old = disable
def disable():
    _disable_old(); cd.third_order_kondo_dIdV = _old["third"]
```

`10_memory_and_retime.py`:

```python
# (a) peak RSS and time of the public ED call at n=10, T=1 K, for 41 / 121 / 241 biases, old vs second-order chunked+occupied (fresh process each)
# (b) 5-repeat median at T=1 K, 41 biases, old vs F+second vs all4
import subprocess, sys, json
CHILD = r'''
import numpy as np, time, resource, sys, json
from dmrgpy import spinchain
import proto
which = tuple(w for w in sys.argv[1].split(",") if w); nb = int(sys.argv[2]); T = float(sys.argv[3]); reps = int(sys.argv[4])
n = 10
sc = spinchain.Spin_Chain(["1/2"]*n)
h = 0
for i in range(n-1):
    h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
sc.set_hamiltonian(h)
proto.enable(which)
eV = np.linspace(-0.01, 0.01, nb)
ts = []
for r in range(reps):
    t0 = time.perf_counter()
    _, d = sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
    ts.append(time.perf_counter()-t0)
print(json.dumps(dict(ts=ts, rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6, chk=float(np.sum(d)))))
'''
open("_child.py", "w").write(CHILD)
def run(which, nb, T, reps):
    out = subprocess.run([sys.executable, "-s", "_child.py", ",".join(which), str(nb), str(T), str(reps)], capture_output=True, text=True)
    return json.loads(out.stdout.strip().splitlines()[-1]) if out.returncode == 0 else dict(err=out.stderr[-300:])
import statistics
print("(a) memory, n=10, T=1 K")
for nb in (41, 121, 241):
    for name, which in (("old", ()), ("second", ("second",))):
        r = run(which, nb, 1.0, 1)
        print("  %3d biases %-7s " % (nb, name), ("time %.1f s, peak RSS %.2f GB, sum %.10g" % (r["ts"][0], r["rss"], r["chk"])) if "ts" in r else r, flush=True)
print("(b) 5-repeat median, n=10, T=1 K, 41 biases (one process per config)")
for name, which in (("old", ()), ("F+second", ("F", "second")), ("all4", ("F", "second", "eig", "F0"))):
    r = run(which, 41, 1.0, 5)
    print("  %-9s median %.2f s runs %s" % (name, statistics.median(r["ts"]), " ".join("%.2f" % t for t in r["ts"])), flush=True)
```

`18_clean_timing.py`:

```python
# clean design: a fresh process per (config, repeat), interleaved, timing the whole public call INCLUDING eigh. ED, n=10 (1024 states), 41 biases, U=0.25, order 3.
import subprocess, sys, json, statistics
CHILD = r'''
import numpy as np, time, sys, json
from dmrgpy import spinchain
import proto
n=int(sys.argv[1]); T=float(sys.argv[2]); nb=int(sys.argv[3]); which=tuple(w for w in sys.argv[4].split(",") if w)
sc = spinchain.Spin_Chain(["1/2"]*n)
h = 0
for i in range(n-1):
    h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
sc.set_hamiltonian(h)
proto.enable(which)
eV = np.linspace(-0.01, 0.01, nb)
t0=time.perf_counter()
_, d = sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
print(json.dumps(dict(t=time.perf_counter()-t0, d=d.tolist())))
'''
open("_child18.py","w").write(CHILD)
import numpy as np
def run(n,T,nb,which):
    o=subprocess.run([sys.executable,"-s","_child18.py",str(n),str(T),str(nb),",".join(which)],capture_output=True,text=True)
    return json.loads(o.stdout.strip().splitlines()[-1])
def study(n,T,nb,configs,reps=3):
    ts={k:[] for k,_ in configs}; ref=None; diffs={}
    for r in range(reps):
        for k,w in configs:
            res=run(n,T,nb,w); ts[k].append(res["t"])
            if ref is None: ref=np.array(res["d"])
            if k not in diffs: diffs[k]=np.max(np.abs(np.array(res["d"])-ref))/np.max(np.abs(ref))
    base=statistics.median(ts[configs[0][0]])
    for k,_ in configs:
        m=statistics.median(ts[k])
        print("  n=%d T=%.0f K %d biases  %-9s median %5.2f s  (%s)  saves %4.1f%%  max rel diff to old %.1e"%(n,T,nb,k,m," ".join("%.2f"%t for t in ts[k]),100*(1-m/base),diffs[k]),flush=True)
cfg=[("old",()),("second",("second",)),("F",("F",)),("eig",("eig",)),("F0",("F0",)),("all",("second","F","eig","F0"))]
for T in (0.0,1.0): study(10,T,41,cfg)
print("blocked third-order against the unpatched repo code, high occupation")
study(9,20.0,21,[("old",()),("third",("third",))])
study(9,5.0,21,[("old",()),("third",("third",))])
```

`07_proto_accuracy.py`:

```python
# old vs new on identical inputs (ED, 6 spins + S=1, both T), largest relative differences
import numpy as np
from dmrgpy import spinchain
import proto
x = np.concatenate([-np.logspace(-9,2.5,300), np.logspace(-9,2.5,300), [0.0]])
from dmrgpy.kondospectrumtk import stepfunctions as sf
from mpmath import mp, mpf, exp, expm1
mp.dps=40
def th_exact(t):
    t=mpf(t); 
    if t==0: return mpf(1)/2
    return mpf(1)/2 + (mp.sinh(t)-t)/(2*(mp.cosh(t)-1))
ref = np.array([float(th_exact(t)) for t in x])
old = sf.Theta(x); new = proto.Theta_new(x)
print("Theta vs 40-digit: max abs err old %.2e (at x=%.1e)  new %.2e (at x=%.1e)" % (np.max(np.abs(old-ref)), x[np.argmax(np.abs(old-ref))], np.max(np.abs(new-ref)), x[np.argmax(np.abs(new-ref))]))
def build():
    sc = spinchain.Spin_Chain(["1/2"]*5+["1"])
    h = 0
    for i in range(5):
        h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2] + 0.0002*sc.Sy[3]
    sc.set_hamiltonian(h); return sc
eV = np.linspace(-0.012, 0.012, 61); import warnings; warnings.simplefilter("ignore")
for T in (0.0, 0.5, 3.0, 40.0):
    sc = build()
    proto.disable(); _, a = sc.get_kondo_spectrum(eV, site=5, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
    for which in (("F",), ("F0",), ("second",), ("eig",), ("third",), ("F","F0","second","eig","third")):
        sc = build(); proto.disable(); proto.enable(which)
        _, b = sc.get_kondo_spectrum(eV, site=5, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
        print("T=%5.1f %-26s max|new-old| = %.2e   relative to peak %.2e" % (T, "+".join(which), np.max(np.abs(a-b)), np.max(np.abs(a-b))/np.max(np.abs(a))))
    proto.disable()
```

`06_profile_ed.py`:

```python
# end-to-end cProfile of mode="ED" get_kondo_spectrum, S=1/2 Heisenberg + field, 8 spins (256 states), 41 biases, U!=0, order 3
import numpy as np, cProfile, pstats, time, sys, resource
from dmrgpy import spinchain
n=int(sys.argv[1]) if len(sys.argv)>1 else 8
def build():
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
    sc.set_hamiltonian(h); return sc
eV = np.linspace(-0.01, 0.01, 41)
for T in (0.0, 1.0):
    sc = build()
    pr = cProfile.Profile(); t0=time.time(); pr.enable()
    sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, order=3, mode="ED")
    pr.disable(); tot=time.time()-t0
    print("\n=== n=%d spins (dim %d), T=%.1f K, 41 biases, U=0.25: total %.2f s, peak RSS %.2f GB" % (n, 2**n, T, tot, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6))
    st = pstats.Stats(pr); st.sort_stats("cumulative").print_stats(14)
```

Observed:

`10_memory_and_retime.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`18_clean_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`07_proto_accuracy.out`:

```
Theta vs 40-digit: max abs err old 8.17e-13 (at x=1.7e-04)  new 5.76e-13 (at x=-1.9e-04)
T=  0.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.0 F0                         max|new-old| = 1.78e-15   relative to peak 1.12e-16
T=  0.0 second                     max|new-old| = 3.55e-15   relative to peak 2.24e-16
T=  0.0 eig                        max|new-old| = 8.88e-16   relative to peak 5.59e-17
T=  0.0 third                      max|new-old| = 1.78e-15   relative to peak 1.12e-16
T=  0.0 F+F0+second+eig+third      max|new-old| = 5.33e-15   relative to peak 3.36e-16
T=  0.5 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.5 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.5 second                     max|new-old| = 5.33e-15   relative to peak 3.39e-16
T=  0.5 eig                        max|new-old| = 3.55e-15   relative to peak 2.26e-16
T=  0.5 third                      max|new-old| = 3.55e-15   relative to peak 2.26e-16
T=  0.5 F+F0+second+eig+third      max|new-old| = 7.11e-15   relative to peak 4.52e-16
T=  3.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  3.0 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  3.0 second                     max|new-old| = 5.10e-12   relative to peak 3.26e-13
T=  3.0 eig                        max|new-old| = 3.55e-15   relative to peak 2.27e-16
T=  3.0 third                      max|new-old| = 3.55e-15   relative to peak 2.27e-16
T=  3.0 F+F0+second+eig+third      max|new-old| = 5.10e-12   relative to peak 3.26e-13
T= 40.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T= 40.0 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T= 40.0 second                     max|new-old| = 5.33e-14   relative to peak 3.52e-15
T= 40.0 eig                        max|new-old| = 5.33e-15   relative to peak 3.52e-16
T= 40.0 third                      max|new-old| = 3.55e-15   relative to peak 2.35e-16
T= 40.0 F+F0+second+eig+third      max|new-old| = 5.51e-14   relative to peak 3.64e-15
```

`06_profile_ed_n10.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

**Reviewer (CONFIRMED (narrowed))**: Real, not documented, and the same defect as task-route-1. The difference between the two findings is exactly the P_CUT step. On TASK p is exactly zero off row 0, so nothing is dropped. On ED at T>0, states with p<1e-12*max(p) are dropped; each dropped column contributes at most 1e-12 of the largest, and there can be up to dim of them, so the bound is about dim*1e-12 relative. On my mock (dim 1024, T=1 K, spread energies) old and new differ by 7.7e-12 relative; the hunter measured 9.8e-15 on real chains with sharper occupations. At T=0 the difference is exactly 0 including a 2-fold degenerate ground state; on real ED chains with every state occupied it is 1.6e-15 to 3.4e-15 (einsum ordering). Narrowing: I re-measured only the function itself, not the public ED call, so the hunter's 35% and 14% public-call savings and the 16.92 GB at 241 biases are carried as theirs; my 9.5 GB at 121 biases scales consistently with it (three arrays of nbias*dim^2*8 bytes).

Reviewer's own evidence:

```
01_second.out: '1024 121 1.0 new 0.30s peak 0.10 GB'; '1024 121 1.0 old 20.30s peak 9.53 GB'; ' maxrel 7.738594311995933e-12'. 02_occ.out: '7 3.0 0.0 dim 128 nocc 128 maxrel 3.4174600487122648e-15'
```

Reviewer's measurement: function level, mock dim 1024, 121 biases, T=1 K, one timed call each: old 20.30 s / 9.53 GB, new 0.30 s / 0.10 GB; max relative difference 7.7e-12 (spread p), 1.6e-15 to 3.4e-15 (real chains, all occupied), 0 at T=0. Public-call percentages are the hunter's, not re-measured.

Reviewer's scripts and outputs:

`01_second.py`:

```python
import numpy as np, time, tracemalloc
from dmrgpy.kondospectrumtk import conductance as cd
from dmrgpy import spinchain
def new(ks, eVs, T0=1.0, U=0.0, p_cut=cd.P_CUT, nbytes=4e6):
    eVs=np.asarray(eVs,float); occ=cd._occupied_states(ks,p_cut)
    M2=cd._spin_matrix_elements_squared(ks)[:,occ]
    eps=ks.e[:,None]-ks.e[None,occ]; w=ks.p[None,occ]*M2
    spin=np.empty(len(eVs)); ch=max(1,int(nbytes)//w.size)
    for s in range(0,len(eVs),ch):
        v=eVs[s:s+ch]
        spin[s:s+len(v)]=(np.einsum('fi,efi->e',w,cd._theta_raw(ks,v[:,None,None]-eps[None]))
          +np.einsum('fi,efi->e',w,cd._theta_raw(ks,-v[:,None,None]-eps[None])))
    return 2*np.pi*T0**2*(spin+4*U**2*np.ones_like(eVs))
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
for case in ["heis n=8 T=0","heis n=8 T=1K","heis n=8 T=30K","degenerate GS T=0 (field-free S=1/2 dimer-less n=3 odd)","complex H n=6 T=2K"]:
    n=8
    sc=spinchain.Spin_Chain([2]*n)
    h=0
    if "n=3" in case:
        sc=spinchain.Spin_Chain([2]*3); n=3
    for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    if "complex" in case:
        n=6; sc=spinchain.Spin_Chain([2]*n); h=0
        for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]+0.3*(sc.Sx[i]*sc.Sy[i+1]-sc.Sy[i]*sc.Sx[i+1])
    h=h+0.05*sc.Sz[0]
    sc.set_hamiltonian(h)
    T=0. if "T=0" in case else float(case.split("T=")[-1].split("K")[0]) if "K" in case else 2.
    ks=KondoSpectrum(sc,0,T)
    eV=np.linspace(-0.01,0.01,41)
    a=cd.second_order_dIdV(ks,eV,U=0.2); b=new(ks,eV,U=0.2)
    print(case,"dim",ks.dim,"nocc",len(cd._occupied_states(ks)),"maxrel",np.max(np.abs(a-b)/np.abs(a)))
# size/time on mock
class K: pass
for dim,nb,T in [(904,241,0.),(1024,121,1.)]:
    ks=K(); ks.dim=dim; rng=np.random.default_rng(0)
    ks.e=np.sort(rng.uniform(0,0.2,dim)); ks.e-=ks.e[0]; ks.T=T; ks.kB=8.617e-5
    for nm in ['Sx','Sy','Sz']: setattr(ks,nm,(rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim)))*0.05)
    if T==0: ks.p=np.zeros(dim); ks.p[0]=1
    else:
        p=np.exp(-ks.e/(ks.kB*T)); ks.p=p/p.sum()
    eV=np.linspace(-0.01,0.01,nb)
    for nm,f in [('new',new),('old',cd.second_order_dIdV)]:
        if nm=='old' and dim*dim*nb*8*6>9e9: print('skip old (mem)');continue
        tracemalloc.start(); t=time.time(); r=f(ks,eV); dt=time.time()-t
        print(dim,nb,T,nm,"%.2fs peak %.2f GB"%(dt,tracemalloc.get_traced_memory()[1]/1e9)); tracemalloc.stop()
        if nm=='new': rn=r
        else: print(" maxrel",np.max(np.abs(r-rn)/np.abs(r)))
```

`02_occ.py`:

```python
import numpy as np, importlib.util, sys
sys.argv=['x']
src=open('01_second.py').read().split("from dmrgpy.kondospectrumtk.edkondo")[0]
exec(src)
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
for n,T,B in [(6,3.,0.0),(7,3.,0.0),(6,10.,0.2),(5,0.,0.0)]:
    sc=spinchain.Spin_Chain([2]*n); h=0
    for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h=h+B*1e-3*sc.Sz[0]; sc.set_hamiltonian(h)
    ks=KondoSpectrum(sc,0,T); eV=np.linspace(-0.01,0.01,41)
    a=cd.second_order_dIdV(ks,eV,U=0.2); b=new(ks,eV,U=0.2)
    print(n,T,B,"dim",ks.dim,"nocc",len(cd._occupied_states(ks)),"maxrel",np.max(np.abs(a-b)/np.abs(a)))
```

`01_second.out`:

```
heis n=8 T=0 dim 256 nocc 1 maxrel 0.0
heis n=8 T=1K dim 256 nocc 1 maxrel 0.0
heis n=8 T=30K dim 256 nocc 1 maxrel 0.0
degenerate GS T=0 (field-free S=1/2 dimer-less n=3 odd) dim 8 nocc 1 maxrel 0.0
complex H n=6 T=2K dim 64 nocc 1 maxrel 0.0
904 241 0.0 new 0.04s peak 0.04 GB
skip old (mem)
1024 121 1.0 new 0.30s peak 0.10 GB
1024 121 1.0 old 20.30s peak 9.53 GB
 maxrel 7.738594311995933e-12
```

`02_occ.out`:

```
6 3.0 0.0 dim 64 nocc 64 maxrel 3.417445875158219e-15
7 3.0 0.0 dim 128 nocc 128 maxrel 3.4174600487122648e-15
6 10.0 0.2 dim 64 nocc 64 maxrel 1.5852746424178214e-15
5 0.0 0.0 dim 32 nocc 2 maxrel 1.5568031143512647e-16
```

**Suggested fix (hunter)**: In second_order_dIdV take occ = _occupied_states(ks) and slice weight and eps to those columns, then evaluate both Theta arrays inside a loop over chunks of biases (chunk = max(1, 4e6//weight.size)), as third_order_kondo_dIdV already does. Document that P_CUT now also applies to the second-order term. Moves numbers by at most 1e-14 relative (3e-13 of the peak at 3 K).

**Fix (reviewer)**: Same patch as task-route-1. Record in docs/numerical_changes.md that the second-order term now uses P_CUT, with the 1e-11 bound; the P_CUT comment ('every term in this module') then becomes literally true. If zero change at T>0 is wanted, chunking without the occ restriction is exact but keeps O(dim^2) work per bias; not worth it. Hunter fix is right.

### 4. FBuilder evaluates F0 on every point and then overwrites it inside the table: 17.7% of a T>0 ED call

`performance` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `ed-reference` &middot; candidate `ed-reference-2`

**Status**: FIXED in `f6ac1e3` on 2026-11-12. `FBuilder.__call__` returns the spline when every point is inside the table and evaluates F0 only on the points outside it. Bit for bit the same. With findings 6 and 7 in the same commit, `mode="ED"` on 10 spins, 41 biases, on d22: 1 K 6.99 s to 5.37 s. Pinned by `tests/test_audit_2026_11_stepfunctions.py::test_fbuilder_call_equals_old_overwrite_form`.

The claim as reviewed: FBuilder.__call__ evaluates F0 on every point and then overwrites the points inside the table, so a T>0 ED call pays for a discarded F0 pass: 17.7% of the public call (14.00 to 11.52 s, n=10 spins, T=1 K, 41 biases), bit-identical.

**Where**: src/dmrgpy/kondospectrumtk/stepfunctions.py:214-220 (FBuilder.__call__: out = F0(x...) then out[inside] = spline)

At T=1 K with the default omega0 the table spans +-10*omega0 = 0.2 eV, so every argument of the exchange and direct logs (eV - e_f + e_m, a few tens of meV) is inside; the F0 pass (two logs per point, 60 calls on up to 1M points each) is entirely discarded: 1.23 s of 8.43 s in the cProfile at n=10. Prototype F_call_new: if inside.all() return the spline directly, else compute F0 only on the points outside. Fresh process per repeat, 3 repeats, whole call including eigh (script 18): T=1 K 7.56 to 6.31 s (16.5%); at T=0 the change does nothing (F0 is the real work there, 0.1%). Output is bit-identical (max difference 0.0 in script 07 at T=0.5, 3 and 40 K). A single-log form of F0, log((w0+|x|)/sqrt(x^2+G0^2)), was also tried and saves 2.4% at T=0 and 1.2% at T=1 K for a 1.3e-16 change: not worth it. A Theta rewritten with one expm1 gave no speedup (not proposed).

The hunter's claim, before review: FBuilder.__call__ evaluates the closed-form F0 on every requested point and then overwrites the points inside the table with the spline, so the T>0 route pays for a full F0 pass it discards: 16.5% of the whole public call (7.56 to 6.31 s, n=10 spins, T=1 K, 41 biases), bit-identical.

Decisive output, quoted by the hunter:

```
n=10 T=1 K 41 biases  old       median  7.56 s  (7.57 7.54 7.56)  saves  0.0%  max rel diff to old 0.0e+00
  n=10 T=1 K 41 biases  F         median  6.31 s  (6.32 6.21 6.31)  saves 16.5%  max rel diff to old 0.0e+00
  n=10 T=0 K 41 biases  F         median  3.47 s  (3.48 3.47 3.47)  saves  0.1%  max rel diff to old 0.0e+00
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`proto.py`:

```python
# prototypes (monkeypatch only), nothing in the repo is edited
import numpy as np
from scipy import sparse
from dmrgpy.kondospectrumtk import stepfunctions as sf, conductance as cd, edkondo
from dmrgpy.algebra import algebra

def Theta_new(x):
    """Theta(x) = (1-t-xt)/(1-t)^2, t=exp(-|x|), mirrored for x<0: one expm1, no sinh/cosh, no masks"""
    x = np.asarray(x, dtype=float)
    a = np.abs(x)
    em1 = np.expm1(-a)
    t = 1.0 + em1
    with np.errstate(divide="ignore", invalid="ignore"):
        pos = (-em1 - a*t)/(em1*em1)
    pos = np.where(a < 1e-4, 0.5 + a/6.0 - a**3/360.0, pos)
    return np.where(x >= 0, pos, 1.0 - pos)

def F_call_new(self, x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    if self._F_spline is None: return self._direct(x)
    inside = (x >= self._F_min) & (x <= self._F_max)
    if inside.all(): return self._F_spline(x)          # F0 only where the table does not reach
    out = np.empty_like(x)
    out[inside] = self._F_spline(x[inside])
    xo = x[~inside]
    out[~inside] = sf.F0(xo, omega0=self.omega0, Gamma0=self.Gamma0)
    return out

def second_order_new(ks, eVs, T0=1.0, U=0.0, p_cut=cd.P_CUT, nbytes=4e6):
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    M2 = cd._spin_matrix_elements_squared(ks)[:, occ]          # [f,i in occ]
    eps = ks.e[:, None] - ks.e[None, occ]
    weight = ks.p[None, occ]*M2
    spin = np.empty(len(eVs))
    chunk = max(1, int(nbytes)//weight.size)
    for s0 in range(0, len(eVs), chunk):
        v = eVs[s0:s0+chunk]
        spin[s0:s0+len(v)] = (np.einsum('fi,efi->e', weight, cd._theta_raw(ks, v[:, None, None] - eps[None]))
                              + np.einsum('fi,efi->e', weight, cd._theta_raw(ks, -v[:, None, None] - eps[None])))
    return 2*np.pi*T0**2*(spin + 4*(U**2)*np.ones_like(eVs))

def eigenbasis_new(self, name):
    edobj = self.chain.get_ED_obj()
    op = getattr(self.chain, name)[self.site]
    mat = edobj.get_operator(op)
    vs = self.vs
    if sparse.issparse(mat):   # sparse @ dense first (O(dim^2)), then ONE dense product instead of two
        mat = mat.tocsr().astype(complex)
        return np.conjugate(vs.T) @ (mat @ vs)
    mat = np.array(algebra.todense(mat), dtype=complex)
    return np.conjugate(vs.T)@mat@vs

def F0_new(x, omega0=20e-3, Gamma0=5e-6):
    x = np.asarray(x, dtype=float)
    ax = np.abs(x)
    return np.log((omega0 + ax)/np.sqrt(x*x + Gamma0*Gamma0))

_old = dict(F0=sf.F0, Theta=sf.Theta, call=sf.FBuilder.__call__, second=cd.second_order_dIdV, eig=edkondo.KondoSpectrum._eigenbasis_spin_operator)
def enable(which):
    if "theta" in which: sf.Theta = Theta_new; cd.Theta = Theta_new
    if "F" in which: sf.FBuilder.__call__ = F_call_new
    if "F0" in which: sf.F0 = F0_new; cd.F0 = F0_new
    if "second" in which: cd.second_order_dIdV = second_order_new
    if "eig" in which: edkondo.KondoSpectrum._eigenbasis_spin_operator = eigenbasis_new
def disable():
    sf.F0 = cd.F0 = _old["F0"]; sf.Theta = cd.Theta = _old["Theta"]; sf.FBuilder.__call__ = _old["call"]
    cd.second_order_dIdV = _old["second"]; edkondo.KondoSpectrum._eigenbasis_spin_operator = _old["eig"]

def third_order_kondo_blocked(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6, Fb=None, p_cut=cd.P_CUT,
                              block_bytes=2.0e8, w_bytes=3.0e8):
    """same sums as conductance.third_order_kondo_dIdV; the coefficient tensor is built for a block of initial
    states at a time (block_bytes), and the exchange diagram is accumulated into W[e,f,m] for a group of biases"""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    dim = ks.dim
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    e_fm = ks.e[:, None] - ks.e[None, :]
    nB = max(1, int(block_bytes//(16*dim*dim)))
    nG = max(1, int(w_bytes//(8*dim*dim)))          # biases per group (W accumulator)
    fchunk = max(1, int(4e6)//(dim*dim))
    def one_direction(v):
        out = np.zeros(len(v))
        for g0 in range(0, len(v), nG):
            vg = v[g0:g0+nG]
            W = np.zeros((len(vg), dim, dim))
            for b0 in range(0, len(occ), nB):
                ob = occ[b0:b0+nB]; p = ks.p[ob]
                coeff = np.imag(cd._triple_product_coefficients(ks, ob))/2.
                eps = ks.e[None, :] - ks.e[ob, None]
                for s0 in range(0, len(vg), max(1, int(4e6)//(dim*len(ob)))):
                    vc = vg[s0:s0+max(1, int(4e6)//(dim*len(ob)))]
                    Th = cd._theta_raw(ks, vc[:, None, None] - eps[None])
                    Fd = Fcall((vc[:, None, None] - eps[None]).ravel()).reshape(len(vc), len(ob), dim)
                    out[g0+s0:g0+s0+len(vc)] += np.einsum('i,ifm,eif,eim->e', p, coeff, Th, Fd, optimize=True)
                    W[s0:s0+len(vc)] += np.einsum('i,ifm,eif->efm', p, coeff, Th, optimize=True)
            for s0 in range(0, len(vg), fchunk):
                vc = vg[s0:s0+fchunk]
                Fx = Fcall((vc[:, None, None] - e_fm[None]).ravel()).reshape(len(vc), dim, dim)
                out[g0+s0:g0+s0+len(vc)] += np.einsum('efm,efm->e', W[s0:s0+len(vc)], Fx)
        return out
    total = one_direction(eVs) + one_direction(-eVs)
    return 4*np.pi*T0**2*Jrho_s*total

_old["third"] = cd.third_order_kondo_dIdV
_enable_old = enable
def enable(which):
    _enable_old(which)
    if "third" in which: cd.third_order_kondo_dIdV = third_order_kondo_blocked
_disable_old = disable
def disable():
    _disable_old(); cd.third_order_kondo_dIdV = _old["third"]
```

`18_clean_timing.py`:

```python
# clean design: a fresh process per (config, repeat), interleaved, timing the whole public call INCLUDING eigh. ED, n=10 (1024 states), 41 biases, U=0.25, order 3.
import subprocess, sys, json, statistics
CHILD = r'''
import numpy as np, time, sys, json
from dmrgpy import spinchain
import proto
n=int(sys.argv[1]); T=float(sys.argv[2]); nb=int(sys.argv[3]); which=tuple(w for w in sys.argv[4].split(",") if w)
sc = spinchain.Spin_Chain(["1/2"]*n)
h = 0
for i in range(n-1):
    h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
sc.set_hamiltonian(h)
proto.enable(which)
eV = np.linspace(-0.01, 0.01, nb)
t0=time.perf_counter()
_, d = sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
print(json.dumps(dict(t=time.perf_counter()-t0, d=d.tolist())))
'''
open("_child18.py","w").write(CHILD)
import numpy as np
def run(n,T,nb,which):
    o=subprocess.run([sys.executable,"-s","_child18.py",str(n),str(T),str(nb),",".join(which)],capture_output=True,text=True)
    return json.loads(o.stdout.strip().splitlines()[-1])
def study(n,T,nb,configs,reps=3):
    ts={k:[] for k,_ in configs}; ref=None; diffs={}
    for r in range(reps):
        for k,w in configs:
            res=run(n,T,nb,w); ts[k].append(res["t"])
            if ref is None: ref=np.array(res["d"])
            if k not in diffs: diffs[k]=np.max(np.abs(np.array(res["d"])-ref))/np.max(np.abs(ref))
    base=statistics.median(ts[configs[0][0]])
    for k,_ in configs:
        m=statistics.median(ts[k])
        print("  n=%d T=%.0f K %d biases  %-9s median %5.2f s  (%s)  saves %4.1f%%  max rel diff to old %.1e"%(n,T,nb,k,m," ".join("%.2f"%t for t in ts[k]),100*(1-m/base),diffs[k]),flush=True)
cfg=[("old",()),("second",("second",)),("F",("F",)),("eig",("eig",)),("F0",("F0",)),("all",("second","F","eig","F0"))]
for T in (0.0,1.0): study(10,T,41,cfg)
print("blocked third-order against the unpatched repo code, high occupation")
study(9,20.0,21,[("old",()),("third",("third",))])
study(9,5.0,21,[("old",()),("third",("third",))])
```

`07_proto_accuracy.py`:

```python
# old vs new on identical inputs (ED, 6 spins + S=1, both T), largest relative differences
import numpy as np
from dmrgpy import spinchain
import proto
x = np.concatenate([-np.logspace(-9,2.5,300), np.logspace(-9,2.5,300), [0.0]])
from dmrgpy.kondospectrumtk import stepfunctions as sf
from mpmath import mp, mpf, exp, expm1
mp.dps=40
def th_exact(t):
    t=mpf(t); 
    if t==0: return mpf(1)/2
    return mpf(1)/2 + (mp.sinh(t)-t)/(2*(mp.cosh(t)-1))
ref = np.array([float(th_exact(t)) for t in x])
old = sf.Theta(x); new = proto.Theta_new(x)
print("Theta vs 40-digit: max abs err old %.2e (at x=%.1e)  new %.2e (at x=%.1e)" % (np.max(np.abs(old-ref)), x[np.argmax(np.abs(old-ref))], np.max(np.abs(new-ref)), x[np.argmax(np.abs(new-ref))]))
def build():
    sc = spinchain.Spin_Chain(["1/2"]*5+["1"])
    h = 0
    for i in range(5):
        h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2] + 0.0002*sc.Sy[3]
    sc.set_hamiltonian(h); return sc
eV = np.linspace(-0.012, 0.012, 61); import warnings; warnings.simplefilter("ignore")
for T in (0.0, 0.5, 3.0, 40.0):
    sc = build()
    proto.disable(); _, a = sc.get_kondo_spectrum(eV, site=5, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
    for which in (("F",), ("F0",), ("second",), ("eig",), ("third",), ("F","F0","second","eig","third")):
        sc = build(); proto.disable(); proto.enable(which)
        _, b = sc.get_kondo_spectrum(eV, site=5, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
        print("T=%5.1f %-26s max|new-old| = %.2e   relative to peak %.2e" % (T, "+".join(which), np.max(np.abs(a-b)), np.max(np.abs(a-b))/np.max(np.abs(a))))
    proto.disable()
```

Observed:

`18_clean_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`07_proto_accuracy.out`:

```
Theta vs 40-digit: max abs err old 8.17e-13 (at x=1.7e-04)  new 5.76e-13 (at x=-1.9e-04)
T=  0.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.0 F0                         max|new-old| = 1.78e-15   relative to peak 1.12e-16
T=  0.0 second                     max|new-old| = 3.55e-15   relative to peak 2.24e-16
T=  0.0 eig                        max|new-old| = 8.88e-16   relative to peak 5.59e-17
T=  0.0 third                      max|new-old| = 1.78e-15   relative to peak 1.12e-16
T=  0.0 F+F0+second+eig+third      max|new-old| = 5.33e-15   relative to peak 3.36e-16
T=  0.5 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.5 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.5 second                     max|new-old| = 5.33e-15   relative to peak 3.39e-16
T=  0.5 eig                        max|new-old| = 3.55e-15   relative to peak 2.26e-16
T=  0.5 third                      max|new-old| = 3.55e-15   relative to peak 2.26e-16
T=  0.5 F+F0+second+eig+third      max|new-old| = 7.11e-15   relative to peak 4.52e-16
T=  3.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  3.0 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  3.0 second                     max|new-old| = 5.10e-12   relative to peak 3.26e-13
T=  3.0 eig                        max|new-old| = 3.55e-15   relative to peak 2.27e-16
T=  3.0 third                      max|new-old| = 3.55e-15   relative to peak 2.27e-16
T=  3.0 F+F0+second+eig+third      max|new-old| = 5.10e-12   relative to peak 3.26e-13
T= 40.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T= 40.0 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T= 40.0 second                     max|new-old| = 5.33e-14   relative to peak 3.52e-15
T= 40.0 eig                        max|new-old| = 5.33e-15   relative to peak 3.52e-16
T= 40.0 third                      max|new-old| = 3.55e-15   relative to peak 2.35e-16
T= 40.0 F+F0+second+eig+third      max|new-old| = 5.51e-14   relative to peak 3.64e-15
```

**Reviewer (CONFIRMED)**: I reproduced it with my own monkeypatch (old/new alternated in one process, 3 repeats, whole public get_kondo_spectrum(mode=ED) call including eigh). Output is bit-identical (reldiff 0.0) for real and complex Hamiltonians (Sy term) at T=0 and 1 K. At T=0 the change does nothing (5.94 vs 5.93 s), as the hunter said. At n=8 the saving is 13% (T=1, real) and 24% (complex, noisy). The behaviour is not documented as intended (the docstring says only that evaluation is O(1) per point). Nothing in kondo.md forbids it. The saving is a fixed cost per call that is spent on F0 of up to 1e6 points, so it grows with the number of biases and states and is irrelevant at T=0. The hunter's single-log F0 variant (2.4%) was not re-measured and rightly is not proposed.

Reviewer's own evidence:

```
cplx=False T=1 old median 14.00 [13.98, 14.0, 18.65] reldiff 0.0e+00
cplx=False T=1 F median 11.52 [11.51, 11.52, 11.71] reldiff 0.0e+00
cplx=False T=0 F median 5.93 (old 5.94)
n=8 cplx=True T=1: old 1.56, F 1.18, reldiff 0.0
```

Reviewer's measurement: own median old/new at n=10, T=1 K, 41 biases: 14.00 s / 11.52 s, saves 17.7%; T=0 saves 0.2%; max relative difference 0.0 (bitwise). Local host was loaded (slower than the hunter's d22 7.56 s) but the fraction agrees with 16.5%.

Reviewer's scripts and outputs:

`02_time.py`:

```python
import sys; sys.path.insert(0,"<repo>/src")
import numpy as np, time, statistics
from scipy import sparse
from dmrgpy import spinchain
from dmrgpy.algebra import algebra
from dmrgpy.kondospectrumtk import stepfunctions as sf, edkondo
def Fnew(self,x):
    x=np.atleast_1d(np.asarray(x,dtype=float))
    if self._F_spline is None: return self._direct(x)
    inside=(x>=self._F_min)&(x<=self._F_max)
    if inside.all(): return self._F_spline(x)
    out=np.empty_like(x); out[inside]=self._F_spline(x[inside]); out[~inside]=sf.F0(x[~inside],omega0=self.omega0,Gamma0=self.Gamma0); return out
def eignew(self,name):
    e=self.chain.get_ED_obj(); op=getattr(self.chain,name)[self.site]; mat=e.get_operator(op)
    if sparse.issparse(mat):
        return np.conjugate(self.vs.T)@(mat.tocsr().astype(complex)@self.vs)
    mat=np.array(algebra.todense(mat),dtype=complex); return np.conjugate(self.vs.T)@mat@self.vs
oF,oE=sf.FBuilder.__call__,edkondo.KondoSpectrum._eigenbasis_spin_operator
n=int(sys.argv[1])
def mk(cplx):
    sc=spinchain.Spin_Chain(["1/2"]*n); h=0
    for i in range(n-1): h=h+0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h=h+0.001*sc.Sz[0]+0.0003*sc.Sx[2]
    if cplx: h=h+0.0005*sc.Sy[1]
    sc.set_hamiltonian(h); return sc
eV=np.linspace(-0.01,0.01,41)
def run(sc,T):
    t=time.perf_counter(); r=sc.get_kondo_spectrum(eV,site=0,Jrho_s=-0.05,U=0.25,T=T,mode="ED")[1]; return time.perf_counter()-t,r
for cplx in (False,):
  sc=mk(cplx)
  for T in (0.,1.):
    res={}
    for rep in range(3):
        for k,(f,e) in {"old":(oF,oE),"F":(Fnew,oE),"eig":(oF,eignew)}.items():
            sf.FBuilder.__call__=f; edkondo.KondoSpectrum._eigenbasis_spin_operator=e
            t,r=run(sc,T); res.setdefault(k,[]).append((t,r))
    ref=res["old"][0][1]
    for k,v in res.items():
        m=statistics.median([a for a,_ in v]); d=max(np.max(abs(r-ref)) for _,r in v)/np.max(abs(ref))
        print("cplx=%s T=%g %-4s median %.2f %s reldiff %.1e"%(cplx,T,k,m,[round(a,2) for a,_ in v],d),flush=True)
```

`03_time10.out`:

```
cplx=False T=0 old  median 5.94 [7.3, 5.93, 5.94] reldiff 0.0e+00
cplx=False T=0 F    median 5.93 [6.05, 5.93, 5.93] reldiff 0.0e+00
cplx=False T=0 eig  median 5.07 [5.07, 5.07, 5.1] reldiff 0.0e+00
cplx=False T=1 old  median 14.00 [13.98, 14.0, 18.65] reldiff 0.0e+00
cplx=False T=1 F    median 11.52 [11.51, 11.52, 11.71] reldiff 0.0e+00
cplx=False T=1 eig  median 13.32 [13.03, 13.32, 17.36] reldiff 1.3e-16
```

`02_time.out`:

```
cplx=False T=0 old  median 0.46 [0.6, 0.46, 0.46] reldiff 0.0e+00
cplx=False T=0 F    median 0.46 [0.48, 0.46, 0.46] reldiff 0.0e+00
cplx=False T=0 eig  median 0.44 [0.44, 0.42, 0.45] reldiff 2.5e-16
cplx=False T=1 old  median 2.09 [2.13, 2.09, 1.74] reldiff 0.0e+00
cplx=False T=1 F    median 1.82 [1.85, 1.82, 1.8] reldiff 0.0e+00
cplx=False T=1 eig  median 2.00 [2.0, 1.85, 2.01] reldiff 2.5e-16
cplx=True T=0 old  median 0.39 [0.56, 0.39, 0.35] reldiff 0.0e+00
cplx=True T=0 F    median 0.36 [0.39, 0.35, 0.36] reldiff 0.0e+00
cplx=True T=0 eig  median 0.34 [0.34, 0.37, 0.26] reldiff 1.2e-16
cplx=True T=1 old  median 1.56 [1.29, 1.56, 1.98] reldiff 0.0e+00
cplx=True T=1 F    median 1.18 [1.15, 1.18, 1.59] reldiff 0.0e+00
cplx=True T=1 eig  median 1.65 [1.3, 1.65, 2.04] reldiff 0.0e+00
```

**Suggested fix (hunter)**: In FBuilder.__call__ build out only where needed: if inside.all(): return self._F_spline(x); otherwise allocate out, fill out[inside] from the spline and out[~inside] from F0. Returns identical numbers.

**Fix (reviewer)**: The hunter's fix is right. In FBuilder.__call__ (stepfunctions.py): if inside.all(): return self._F_spline(x); otherwise out = np.empty_like(x), out[inside] = spline, out[~inside] = F0(x[~inside]). Pin with a test in tests/test_kondo_spectrum.py that compares the new __call__ against the old overwrite form for an x array straddling the table edge (all inside, all outside, mixed) with exact equality.

### 5. The third-order coefficient tensor is built for every occupied state at once: 15 GB at 917 of 1024 states

`performance` &middot; severity **MEDIUM** &middot; CONFIRMED (narrowed) &middot; lens `ed-reference` &middot; candidate `ed-reference-3`

**Status**: FIXED in `49c23ca` on 2026-11-12, in the reviewer's form. The coefficients are built per block of occupied initial states under `COEFF_BYTES` = 200 MB, with both tunnelling directions inside each block, and the exchange diagram is accumulated per group of biases under `EXCHANGE_BYTES` = 600 MB. One block is the old code path. 9 spins, 20 K, all 512 states occupied, 41 biases, on d22: 3.22 GB to 0.59 GB and 8.05 s to 7.53 s; one block (10 spins, 1 K): 3.24 s against 3.23 s. Roundoff only (4e-15). The 917-state case at 10 spins and 10 K was not rerun. Pinned by `tests/test_audit_2026_11_conductance.py::test_third_order_blocked_equals_single_pass` and `::test_third_order_coefficient_memory_follows_the_budget`.

The claim as reviewed: third_order_kondo_dIdV builds the full (n_occ,dim,dim) complex coefficient tensor at once (16*n_occ*dim^2 bytes, 15 GB at n_occ 917 and dim 1024; the hunter measured a 22.9 GB peak), so high-T ED does not fit in 16 GB. Blocking over initial states fixes the memory with identical numbers (4e-15), but the hunter's prototype costs +36 to +46% time because it rebuilds the coefficient tensor once per tunneling direction. With both directions inside each block the cost is free: with all biases in one W group, -1% to -8% against old at lower memory; each extra bias group costs one more rebuild (+19% with several groups).

**Where**: src/dmrgpy/kondospectrumtk/conductance.py:189-203 (_triple_product_coefficients) and 214-231 (coeff = imag(...)/2 for all of occ, einsums over it)

The docstring says the third-order memory is bounded by chunking biases, but that only bounds the F arrays; coeff[i,f,m] has n_occ*dim^2 entries (complex during the einsum, then float) whatever the chunking. n_occ is the number of states within about 28 kT of the ground state, so it is 4 at 1 K but 176 of 1024 at 5 K and 917 of 1024 at 10 K (script 14). Measured peak RSS, second-order and F patched so they do not dominate: n=9: 0.72 GB (1 K, n_occ 4), 1.13 GB (5 K, 127), 3.52 GB (20 K, all 512); n=10: 0.85 GB, 4.82 GB (5 K, 176), 22.99 GB (10 K, 917). A 16 GB machine cannot run a 10-spin chain at 10 K (nor a 9-spin chain at about 100 K); the ED route is the only one that gives the thermal state, so this is where users go at T>0. Prototype third_order_kondo_blocked (proto.py): coefficients for blocks of initial states (200 MB), direct diagram contracted per block, exchange diagram accumulated into W[e,f,m] for a group of biases and contracted with Fx at the end. Measured against the unpatched repo code, fresh process, 3 repeats (script 18): n=9 T=20 K 21 biases 7.82 to 11.40 s (+46%), n=9 T=5 K 2.85 to 3.80 s (+33%); single runs on d22 (script 15, baseline had the other patches on): n=10 T=5 K 4.75 to 0.94 GB (14.8 to 21.1 s), n=10 T=10 K 22.92 to 0.96 GB (77.7 to 98.3 s), n=9 T=20 K 3.51 to 0.79 GB. Max relative difference to the old code 1.8e-15 (n=9), 7.8e-15 and 5.2e-15 (n=10). The extra time is the coefficient einsum re-run per bias group and the per-block einsum overhead; at 1 K (n_occ=4) nothing changes. Offered as a memory fix, not a speedup. Separate observation: at high occupation the coefficient einsum alone is 3.61 s of 10.6 s at n=9, T=20 K (profile 16); a direct expression of the six non-zero Levi-Civita terms is a lead, not tested.

The hunter's claim, before review: third_order_kondo_dIdV builds the full (n_occ, dim, dim) complex Levi-Civita coefficient tensor at once, so the ED route at 10 spins 1/2 and T=10 K (917 of 1024 states occupied) peaks at 22.9 GB; building it per block of initial states gives 0.96 GB with identical numbers (5e-15 relative) at a cost of about 26 to 46% more time.

Decisive output, quoted by the hunter:

```
n=10 dim=1024 T=  5.0 K: n_occ= 176 of 1024  time 14.3 s  peak RSS 4.82 GB (before call 0.42)
n=10 dim=1024 T= 10.0 K: n_occ= 917 of 1024  time 80.2 s  peak RSS 22.99 GB (before call 0.42)
n=10 dim=1024 T=10.0 K 21 biases | old   77.7 s 22.92 GB | blocked   98.3 s  0.96 GB | max rel diff 5.2e-15
  n=9 T=20 K 21 biases  old       median  7.82 s  (7.79 8.22 7.82)  saves  0.0%  max rel diff to old 0.0e+00
  n=9 T=20 K 21 biases  third     median 11.40 s  (11.40 11.51 11.35)  saves -45.8%  max rel diff to old 6.7e-16
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`proto.py`:

```python
# prototypes (monkeypatch only), nothing in the repo is edited
import numpy as np
from scipy import sparse
from dmrgpy.kondospectrumtk import stepfunctions as sf, conductance as cd, edkondo
from dmrgpy.algebra import algebra

def Theta_new(x):
    """Theta(x) = (1-t-xt)/(1-t)^2, t=exp(-|x|), mirrored for x<0: one expm1, no sinh/cosh, no masks"""
    x = np.asarray(x, dtype=float)
    a = np.abs(x)
    em1 = np.expm1(-a)
    t = 1.0 + em1
    with np.errstate(divide="ignore", invalid="ignore"):
        pos = (-em1 - a*t)/(em1*em1)
    pos = np.where(a < 1e-4, 0.5 + a/6.0 - a**3/360.0, pos)
    return np.where(x >= 0, pos, 1.0 - pos)

def F_call_new(self, x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    if self._F_spline is None: return self._direct(x)
    inside = (x >= self._F_min) & (x <= self._F_max)
    if inside.all(): return self._F_spline(x)          # F0 only where the table does not reach
    out = np.empty_like(x)
    out[inside] = self._F_spline(x[inside])
    xo = x[~inside]
    out[~inside] = sf.F0(xo, omega0=self.omega0, Gamma0=self.Gamma0)
    return out

def second_order_new(ks, eVs, T0=1.0, U=0.0, p_cut=cd.P_CUT, nbytes=4e6):
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    M2 = cd._spin_matrix_elements_squared(ks)[:, occ]          # [f,i in occ]
    eps = ks.e[:, None] - ks.e[None, occ]
    weight = ks.p[None, occ]*M2
    spin = np.empty(len(eVs))
    chunk = max(1, int(nbytes)//weight.size)
    for s0 in range(0, len(eVs), chunk):
        v = eVs[s0:s0+chunk]
        spin[s0:s0+len(v)] = (np.einsum('fi,efi->e', weight, cd._theta_raw(ks, v[:, None, None] - eps[None]))
                              + np.einsum('fi,efi->e', weight, cd._theta_raw(ks, -v[:, None, None] - eps[None])))
    return 2*np.pi*T0**2*(spin + 4*(U**2)*np.ones_like(eVs))

def eigenbasis_new(self, name):
    edobj = self.chain.get_ED_obj()
    op = getattr(self.chain, name)[self.site]
    mat = edobj.get_operator(op)
    vs = self.vs
    if sparse.issparse(mat):   # sparse @ dense first (O(dim^2)), then ONE dense product instead of two
        mat = mat.tocsr().astype(complex)
        return np.conjugate(vs.T) @ (mat @ vs)
    mat = np.array(algebra.todense(mat), dtype=complex)
    return np.conjugate(vs.T)@mat@vs

def F0_new(x, omega0=20e-3, Gamma0=5e-6):
    x = np.asarray(x, dtype=float)
    ax = np.abs(x)
    return np.log((omega0 + ax)/np.sqrt(x*x + Gamma0*Gamma0))

_old = dict(F0=sf.F0, Theta=sf.Theta, call=sf.FBuilder.__call__, second=cd.second_order_dIdV, eig=edkondo.KondoSpectrum._eigenbasis_spin_operator)
def enable(which):
    if "theta" in which: sf.Theta = Theta_new; cd.Theta = Theta_new
    if "F" in which: sf.FBuilder.__call__ = F_call_new
    if "F0" in which: sf.F0 = F0_new; cd.F0 = F0_new
    if "second" in which: cd.second_order_dIdV = second_order_new
    if "eig" in which: edkondo.KondoSpectrum._eigenbasis_spin_operator = eigenbasis_new
def disable():
    sf.F0 = cd.F0 = _old["F0"]; sf.Theta = cd.Theta = _old["Theta"]; sf.FBuilder.__call__ = _old["call"]
    cd.second_order_dIdV = _old["second"]; edkondo.KondoSpectrum._eigenbasis_spin_operator = _old["eig"]

def third_order_kondo_blocked(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6, Fb=None, p_cut=cd.P_CUT,
                              block_bytes=2.0e8, w_bytes=3.0e8):
    """same sums as conductance.third_order_kondo_dIdV; the coefficient tensor is built for a block of initial
    states at a time (block_bytes), and the exchange diagram is accumulated into W[e,f,m] for a group of biases"""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    dim = ks.dim
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    e_fm = ks.e[:, None] - ks.e[None, :]
    nB = max(1, int(block_bytes//(16*dim*dim)))
    nG = max(1, int(w_bytes//(8*dim*dim)))          # biases per group (W accumulator)
    fchunk = max(1, int(4e6)//(dim*dim))
    def one_direction(v):
        out = np.zeros(len(v))
        for g0 in range(0, len(v), nG):
            vg = v[g0:g0+nG]
            W = np.zeros((len(vg), dim, dim))
            for b0 in range(0, len(occ), nB):
                ob = occ[b0:b0+nB]; p = ks.p[ob]
                coeff = np.imag(cd._triple_product_coefficients(ks, ob))/2.
                eps = ks.e[None, :] - ks.e[ob, None]
                for s0 in range(0, len(vg), max(1, int(4e6)//(dim*len(ob)))):
                    vc = vg[s0:s0+max(1, int(4e6)//(dim*len(ob)))]
                    Th = cd._theta_raw(ks, vc[:, None, None] - eps[None])
                    Fd = Fcall((vc[:, None, None] - eps[None]).ravel()).reshape(len(vc), len(ob), dim)
                    out[g0+s0:g0+s0+len(vc)] += np.einsum('i,ifm,eif,eim->e', p, coeff, Th, Fd, optimize=True)
                    W[s0:s0+len(vc)] += np.einsum('i,ifm,eif->efm', p, coeff, Th, optimize=True)
            for s0 in range(0, len(vg), fchunk):
                vc = vg[s0:s0+fchunk]
                Fx = Fcall((vc[:, None, None] - e_fm[None]).ravel()).reshape(len(vc), dim, dim)
                out[g0+s0:g0+s0+len(vc)] += np.einsum('efm,efm->e', W[s0:s0+len(vc)], Fx)
        return out
    total = one_direction(eVs) + one_direction(-eVs)
    return 4*np.pi*T0**2*Jrho_s*total

_old["third"] = cd.third_order_kondo_dIdV
_enable_old = enable
def enable(which):
    _enable_old(which)
    if "third" in which: cd.third_order_kondo_dIdV = third_order_kondo_blocked
_disable_old = disable
def disable():
    _disable_old(); cd.third_order_kondo_dIdV = _old["third"]
```

`14_thermal_memory.py`:

```python
# third-order memory/time vs number of thermally occupied initial states (second-order patched so it does not dominate), fresh process each
import subprocess, sys, json
CHILD = r'''
import numpy as np, time, resource, sys, json, cProfile, pstats, io
from dmrgpy import spinchain
import proto
from dmrgpy.kondospectrumtk import conductance as cd
n=int(sys.argv[1]); T=float(sys.argv[2]); nb=int(sys.argv[3])
sc = spinchain.Spin_Chain(["1/2"]*n)
h = 0
for i in range(n-1):
    h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
sc.set_hamiltonian(h)
proto.enable(("second","F"))
eV = np.linspace(-0.01, 0.01, nb)
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
ks = KondoSpectrum(sc, 0, T); nocc = len(cd._occupied_states(ks))
rss0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6
t0=time.perf_counter()
_, d = sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
print(json.dumps(dict(nocc=nocc, dim=2**n, t=time.perf_counter()-t0, rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6, rss_before=rss0)))
'''
open("_child14.py","w").write(CHILD)
for n,T in ((9,1.0),(9,5.0),(9,20.0),(10,1.0),(10,5.0),(10,10.0)):
    out = subprocess.run([sys.executable,"-s","_child14.py",str(n),str(T),"21"],capture_output=True,text=True)
    if out.returncode: print(n,T,"FAILED",out.stderr[-200:]); continue
    r=json.loads(out.stdout.strip().splitlines()[-1])
    print("n=%d dim=%d T=%5.1f K: n_occ=%4d of %d  time %.1f s  peak RSS %.2f GB (before call %.2f)"%(n,r["dim"],T,r["nocc"],r["dim"],r["t"],r["rss"],r["rss_before"]),flush=True)
```

`15_blocked_third_order.py`:

```python
# third-order Kondo term with the coefficient tensor built per block of initial states: old vs new, accuracy local, time/memory on d22
import subprocess, sys, json
CHILD = r'''
import numpy as np, time, resource, sys, json
from dmrgpy import spinchain
import proto
n=int(sys.argv[1]); T=float(sys.argv[2]); nb=int(sys.argv[3]); which=tuple(w for w in sys.argv[4].split(",") if w)
sc = spinchain.Spin_Chain(["1/2"]*n)
h = 0
for i in range(n-1):
    h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
sc.set_hamiltonian(h)
proto.enable(("second","F")+which)
eV = np.linspace(-0.01, 0.01, nb)
t0=time.perf_counter()
_, d = sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
print(json.dumps(dict(t=time.perf_counter()-t0, rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6, d=d.tolist())))
'''
open("_child15.py","w").write(CHILD)
import numpy as np
cases = [(7,5.0,21),(9,5.0,21),(9,20.0,21),(10,5.0,21),(10,10.0,21)]
if len(sys.argv)>1 and sys.argv[1]=="small": cases=cases[:1]
for n,T,nb in cases:
    res={}
    for name,which in (("old",""),("blocked","third")):
        out = subprocess.run([sys.executable,"-s","_child15.py",str(n),str(T),str(nb),which],capture_output=True,text=True)
        res[name]=json.loads(out.stdout.strip().splitlines()[-1]) if out.returncode==0 else None
        if res[name] is None: print(n,T,name,"FAILED",out.stderr[-300:])
    if all(res.values()):
        a=np.array(res["old"]["d"]); b=np.array(res["blocked"]["d"])
        print("n=%2d dim=%4d T=%4.1f K %d biases | old %6.1f s %5.2f GB | blocked %6.1f s %5.2f GB | max rel diff %.1e"%(n,2**n,T,nb,res["old"]["t"],res["old"]["rss"],res["blocked"]["t"],res["blocked"]["rss"],np.max(np.abs(a-b))/np.max(np.abs(a))),flush=True)
```

`18_clean_timing.py`:

```python
# clean design: a fresh process per (config, repeat), interleaved, timing the whole public call INCLUDING eigh. ED, n=10 (1024 states), 41 biases, U=0.25, order 3.
import subprocess, sys, json, statistics
CHILD = r'''
import numpy as np, time, sys, json
from dmrgpy import spinchain
import proto
n=int(sys.argv[1]); T=float(sys.argv[2]); nb=int(sys.argv[3]); which=tuple(w for w in sys.argv[4].split(",") if w)
sc = spinchain.Spin_Chain(["1/2"]*n)
h = 0
for i in range(n-1):
    h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
sc.set_hamiltonian(h)
proto.enable(which)
eV = np.linspace(-0.01, 0.01, nb)
t0=time.perf_counter()
_, d = sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
print(json.dumps(dict(t=time.perf_counter()-t0, d=d.tolist())))
'''
open("_child18.py","w").write(CHILD)
import numpy as np
def run(n,T,nb,which):
    o=subprocess.run([sys.executable,"-s","_child18.py",str(n),str(T),str(nb),",".join(which)],capture_output=True,text=True)
    return json.loads(o.stdout.strip().splitlines()[-1])
def study(n,T,nb,configs,reps=3):
    ts={k:[] for k,_ in configs}; ref=None; diffs={}
    for r in range(reps):
        for k,w in configs:
            res=run(n,T,nb,w); ts[k].append(res["t"])
            if ref is None: ref=np.array(res["d"])
            if k not in diffs: diffs[k]=np.max(np.abs(np.array(res["d"])-ref))/np.max(np.abs(ref))
    base=statistics.median(ts[configs[0][0]])
    for k,_ in configs:
        m=statistics.median(ts[k])
        print("  n=%d T=%.0f K %d biases  %-9s median %5.2f s  (%s)  saves %4.1f%%  max rel diff to old %.1e"%(n,T,nb,k,m," ".join("%.2f"%t for t in ts[k]),100*(1-m/base),diffs[k]),flush=True)
cfg=[("old",()),("second",("second",)),("F",("F",)),("eig",("eig",)),("F0",("F0",)),("all",("second","F","eig","F0"))]
for T in (0.0,1.0): study(10,T,41,cfg)
print("blocked third-order against the unpatched repo code, high occupation")
study(9,20.0,21,[("old",()),("third",("third",))])
study(9,5.0,21,[("old",()),("third",("third",))])
```

`16_profile_highT.py`:

```python
# end-to-end cProfile of mode="ED" get_kondo_spectrum, S=1/2 Heisenberg + field, 8 spins (256 states), 41 biases, U!=0, order 3
import numpy as np, cProfile, pstats, time, sys, resource
from dmrgpy import spinchain
n=int(sys.argv[1]) if len(sys.argv)>1 else 8
def build():
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
    sc.set_hamiltonian(h); return sc
eV = np.linspace(-0.01, 0.01, 41)
for T in (20.0,):
    sc = build()
    pr = cProfile.Profile(); t0=time.time(); pr.enable()
    sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, order=3, mode="ED")
    pr.disable(); tot=time.time()-t0
    print("\n=== n=%d spins (dim %d), T=%.1f K, 41 biases, U=0.25: total %.2f s, peak RSS %.2f GB" % (n, 2**n, T, tot, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6))
    st = pstats.Stats(pr); st.sort_stats("cumulative").print_stats(12)
```

Observed:

`14_thermal_memory.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`15_blocked_third_order.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`18_clean_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`16_profile_highT.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

**Reviewer (CONFIRMED (narrowed))**: The memory defect is real, and the docstring's chunking only bounds the F arrays. The coefficient tensor is n_occ*dim^2 complex whatever the bias chunking, which exceeds 15 GB at n_occ 917, dim 1024 by arithmetic. I did not rerun that case; I measured dim 384 with 94 occupied states (old 0.50 GB, blocked 0.38 GB), which shows the direction only. The hunter's time cost survived on my mock (+36%, the same at 200 MB and 2 GB blocks), so block size alone does not make it free. Profiling (script 04) shows the coefficient einsum is 1.5 s of the 2.6 to 2.9 s old call, and the prototype's one_direction() rebuilds it for -v. Script 05 puts both directions inside the block loop with one W accumulator per direction: 2.83 s against old 2.86 s at 200 MB blocks, 2.64 s at 2 GB, max relative difference 4.6e-15 and 3.9e-15, peak 0.38 GB against 0.50 GB. The cost grows with the number of bias groups (two W accumulators of nbias*dim^2 floats against the budget): at dim 1024 and 241 biases W for two directions is about 4 GB, so a byte budget there means several groups and a rebuild per group; I did not measure that size. Struck: the claim that the 30 to 45% price is inherent. Not tested by me: the d22 single runs at n=10, and complex-H or sector cases (same einsum, new block axis only).

Reviewer's own evidence:

```
03_third.out: 'old 2.62s 0.50GB'; 'blocked 200000000.0 300000000.0 3.56s (+36%) 0.38GB maxrel 4.6e-15'; 'blocked 2000000000.0 3000000000.0 3.57s (+36%)'. 04_pieces.out: 'coeff 1.506 s'. 05_fixed.out: 'blocked2 block 200000000.0 2.83s (-1%) 0.38GB maxrel 4.64e-15'; 'blocked2 block 2000000000.0 2.64s (-8%)'; 'blocked2 small W group (several groups) 3.39s (+19%)'
```

Reviewer's measurement: mock dim 384, 94 occupied, 21 biases, T=20 K, median of 3: old 2.62 to 2.86 s / 0.50 GB; hunter blocked prototype 3.56 s (+36%) / 0.38 GB, independent of block size; both directions per block 2.64 to 2.83 s (-8% to -1%) / 0.38 GB; several W groups 3.39 s (+19%). Max relative difference 3.9e-15 to 4.6e-15, no change beyond roundoff. The 22.9 GB peak at n=10, T=10 K is the hunter's measurement, not rerun.

Reviewer's scripts and outputs:

`03_third.py`:

```python
import numpy as np, time, tracemalloc, sys
sys.path.insert(0,'<scratch>/hunt/ed-reference')
import proto
from dmrgpy.kondospectrumtk import conductance as cd
cd.second_order_dIdV=proto._old['second']
class K: pass
dim=384; rng=np.random.default_rng(1); ks=K(); ks.dim=dim
ks.e=np.sort(rng.uniform(0,0.2,dim)); ks.e-=ks.e[0]; ks.T=20.; ks.kB=8.617e-5
for nm in ['Sx','Sy','Sz']:
    a=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim)); setattr(ks,nm,(a+a.conj().T)*0.05)
p=np.exp(-ks.e/(ks.kB*ks.T)); ks.p=p/p.sum()
print("nocc",len(cd._occupied_states(ks)),"dim",dim)
eV=np.linspace(-0.01,0.01,21)
def run(f,**kw):
    tracemalloc.start(); ts=[]
    for _ in range(3):
        t=time.time(); r=f(ks,eV,-0.05,**kw); ts.append(time.time()-t)
    pk=tracemalloc.get_traced_memory()[1]/1e9; tracemalloc.stop(); return r,np.median(ts),pk
r0,t0,m0=run(cd.third_order_kondo_dIdV); print("old %.2fs %.2fGB"%(t0,m0))
for bb,wb in [(2e8,3e8),(2e9,3e9),(6e8,3e8)]:
    r,t,m=run(proto.third_order_kondo_blocked,block_bytes=bb,w_bytes=wb)
    print("blocked",bb,wb,"%.2fs (%+.0f%%) %.2fGB"%(t,100*(t/t0-1),m),"maxrel",np.max(np.abs(r-r0)/np.abs(r0)))
```

`04_pieces.py`:

```python
import numpy as np, time
exec(open('03_third.py').read().split("def run")[0])
occ=cd._occupied_states(ks); pp=ks.p[occ]; nv=len(eV)
def T(label,f,n=3):
    ts=[]
    for _ in range(n):
        t=time.time(); r=f(); ts.append(time.time()-t)
    print("%-28s %.3f s"%(label,np.median(ts))); return r
coeff=T("coeff",lambda:np.imag(cd._triple_product_coefficients(ks,occ))/2.)
eps=ks.e[None,:]-ks.e[occ,None]
Th=T("Theta",lambda:cd._theta_raw(ks,eV[:,None,None]-eps[None]))
Fd=np.random.rand(nv,len(occ),dim); Fx=np.random.rand(nv,dim,dim)
T("old direct einsum",lambda:np.einsum('i,ifm,eif,eim->e',pp,coeff,Th,Fd,optimize=True))
T("old exch einsum",lambda:np.einsum('i,ifm,eif,efm->e',pp,coeff,Th,Fx,optimize=True))
W=T("W einsum",lambda:np.einsum('i,ifm,eif->efm',pp,coeff,Th,optimize=True))
def Wmm():
    C=(pp[:,None,None]*coeff).transpose(1,0,2)  # f,i,m
    return np.matmul(Th.transpose(2,0,1),C).transpose(1,0,2)  # f,e,i @ f,i,m -> f,e,m
W2=T("W matmul",Wmm); print("maxdiff",np.abs(W-W2).max())
```

`05_fixed.py`:

```python
import numpy as np, time, tracemalloc
exec(open('03_third.py').read().split("def run")[0])
def blocked2(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6, Fb=None, p_cut=cd.P_CUT, block_bytes=2e8, w_bytes=6e8):
    eVs=np.asarray(eVs,float); occ=cd._occupied_states(ks,p_cut); dim=ks.dim
    Fcall=cd._get_F(ks,omega0,Gamma0,Fb); e_fm=ks.e[:,None]-ks.e[None,:]
    nB=max(1,int(block_bytes//(16*dim*dim)))
    nG=max(1,int(w_bytes//(2*8*dim*dim)))   # biases per group, both directions
    fch=max(1,int(4e6)//(dim*dim)); out=np.zeros(len(eVs))
    for g0 in range(0,len(eVs),nG):
        vg=eVs[g0:g0+nG]; vs=[vg,-vg]; Ws=[np.zeros((len(vg),dim,dim)) for _ in vs]; o=np.zeros(len(vg))
        for b0 in range(0,len(occ),nB):
            ob=occ[b0:b0+nB]; p=ks.p[ob]
            coeff=np.imag(cd._triple_product_coefficients(ks,ob))/2.
            eps=ks.e[None,:]-ks.e[ob,None]
            for v,W in zip(vs,Ws):
                ch=max(1,int(4e6)//(dim*len(ob)))
                for s0 in range(0,len(v),ch):
                    vc=v[s0:s0+ch]
                    Th=cd._theta_raw(ks,vc[:,None,None]-eps[None])
                    Fd=Fcall((vc[:,None,None]-eps[None]).ravel()).reshape(len(vc),len(ob),dim)
                    o[s0:s0+len(vc)]+=np.einsum('i,ifm,eif,eim->e',p,coeff,Th,Fd,optimize=True)
                    W[s0:s0+len(vc)]+=np.einsum('i,ifm,eif->efm',p,coeff,Th,optimize=True)
        for v,W in zip(vs,Ws):
            for s0 in range(0,len(v),fch):
                vc=v[s0:s0+fch]
                Fx=Fcall((vc[:,None,None]-e_fm[None]).ravel()).reshape(len(vc),dim,dim)
                o[s0:s0+len(vc)]+=np.einsum('efm,efm->e',W[s0:s0+len(vc)],Fx)
        out[g0:g0+len(vg)]=o
    return 4*np.pi*T0**2*Jrho_s*out
def run(f,**kw):
    tracemalloc.start(); ts=[]
    for _ in range(3):
        t=time.time(); r=f(ks,eV,-0.05,**kw); ts.append(time.time()-t)
    pk=tracemalloc.get_traced_memory()[1]/1e9; tracemalloc.stop(); return r,np.median(ts),pk
r0,t0,m0=run(cd.third_order_kondo_dIdV); print("old %.2fs %.2fGB"%(t0,m0))
for bb in [2e8,2e9]:
    r,t,m=run(blocked2,block_bytes=bb); print("blocked2 block",bb,"%.2fs (%+.0f%%) %.2fGB"%(t,100*(t/t0-1),m),"maxrel",np.max(np.abs(r-r0)/np.abs(r0)))
r,t,m=run(blocked2,block_bytes=2e8,w_bytes=6e7); print("blocked2 small W group (several groups) %.2fs (%+.0f%%) %.2fGB"%(t,100*(t/t0-1),m),np.max(np.abs(r-r0)/np.abs(r0)))
```

`03_third.out`:

```
nocc 94 dim 384
old 2.62s 0.50GB
blocked 200000000.0 300000000.0 3.56s (+36%) 0.38GB maxrel 4.641353040244497e-15
blocked 2000000000.0 3000000000.0 3.57s (+36%) 0.38GB maxrel 3.8677942002037475e-15
blocked 600000000.0 300000000.0 3.84s (+47%) 0.38GB maxrel 3.8677942002037475e-15
```

`04_pieces.out`:

```
nocc 94 dim 384
coeff                        1.506 s
Theta                        0.076 s
old direct einsum            0.102 s
old exch einsum              0.105 s
W einsum                     0.102 s
W matmul                     0.122 s
maxdiff 4.336808689942018e-19
```

`05_fixed.out`:

```
nocc 94 dim 384
old 2.86s 0.50GB
blocked2 block 200000000.0 2.83s (-1%) 0.38GB maxrel 4.641353040244497e-15
blocked2 block 2000000000.0 2.64s (-8%) 0.38GB maxrel 3.8677942002037475e-15
blocked2 small W group (several groups) 3.39s (+19%) 0.38GB 4.641353040244497e-15
```

**Suggested fix (hunter)**: Build the coefficient tensor for blocks of occupied initial states sized by a byte budget (about 200 MB), contract the direct diagram per block, and accumulate the exchange diagram into W[e,f,m] = sum_i p_i coeff Th for a group of biases sized by a second budget, contracting it with Fx once per group. Same sums, memory independent of n_occ, 1e-14 relative change; the price is 30 to 45% more time at high occupation, so keep the old single pass when n_occ*dim^2*16 bytes is below the budget.

**Fix (reviewer)**: In third_order_kondo_dIdV: pick a block of occupied initial states by a byte budget (about 200 MB of complex coefficients, nB = budget // (16*dim*dim)). For each bias group (sized so two W accumulators of nbias*dim^2 floats fit a second budget, about 600 MB), loop over blocks, build coeff once per block, and for both v and -v accumulate the direct diagram into the result and the exchange diagram into W[e,f,m]=sum_i p_i coeff Th; after the block loop contract each W with Fx in the existing bias chunks. This is one coefficient build per direction pair per group instead of the hunter's one per direction, which is the whole 36 to 46% cost. When nocc*dim^2*16 is under the budget it is one block and costs what the present code costs. Pin with a test comparing blocked against single-pass at p_cut=0 with tiny forced block and group budgets (T>0, several occupied states, complex H), plus a tracemalloc cap at dim of about 250 with all states occupied.

### 6. The ED eigenbasis spin operator is densified first: 14.6% of a T=0 ED call

`performance` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `ed-reference` &middot; candidate `ed-reference-4`

**Status**: FIXED in `f6ac1e3` on 2026-11-12. `conj(vs.T) @ (csr @ vs)` with the csr cast to complex, and the dense branch kept for a dense operator. T=0, 10 spins, 41 biases, on d22: 3.43 s to 3.05 s, 1.3e-16 apart. Pinned by `tests/test_audit_2026_11_stepfunctions.py::test_eigenbasis_spin_operator_spin_chain` and `::test_eigenbasis_spin_operator_hubbard_sector` (Nf=3 sector).

The claim as reviewed: KondoSpectrum._eigenbasis_spin_operator densifies the sparse spin operator and does two complex dim^3 products per operator; conj(vs.T) @ (csr @ vs) saves 14.6% of the whole ED call at T=0 (5.94 to 5.07 s, n=10) and 4.9% at T=1 K, for a roundoff (1.3e-16) change.

**Where**: src/dmrgpy/kondospectrumtk/edkondo.py:62-71

get_operator returns a complex csr_matrix; the code does algebra.todense, np.array(..., dtype=complex) and np.conjugate(vs.T)@mat@vs: 0.86 s for three operators at n=10 (22% of the T=0 call, profile 06). Prototype eigenbasis_new: np.conjugate(vs.T) @ (mat_csr @ vs), the sparse product being O(dim^2) and leaving one dense ZGEMM. Fresh process per repeat, 3 repeats, whole call including eigh (script 18): T=0 3.47 to 3.05 s (12.1%), T=1 K 7.56 to 7.11 s (5.9%); max relative difference 1.3e-16 (T=0) and 2.5e-16 (T=1 K). Combined with the three proposals above, all four changes together give 3.47 to 2.48 s at T=0 (28.5%) and 7.56 to 3.23 s at T=1 K (57.3%), max relative difference 3.8e-16 and 9.8e-15. Note: the first timing runs of scripts 08 and 09 were contaminated by other jobs on d22 and are not cited; script 18 is the clean one. eigh itself (0.51 s at n=10) is irreducible.

The hunter's claim, before review: KondoSpectrum._eigenbasis_spin_operator densifies the sparse spin operator and does two complex dim^3 products per operator; sparse-times-dense first and one product saves 12% of the whole call at T=0 (3.47 to 3.05 s, n=10) and 5.9% at T=1 K, for a 1.3e-16 change.

Decisive output, quoted by the hunter:

```
n=10 T=0 K 41 biases  eig       median  3.05 s  (3.05 3.05 3.05)  saves 12.1%  max rel diff to old 1.3e-16
  n=10 T=0 K 41 biases  all       median  2.48 s  (2.48 2.49 2.48)  saves 28.5%  max rel diff to old 3.8e-16
  n=10 T=1 K 41 biases  all       median  3.23 s  (3.23 3.23 3.23)  saves 57.3%  max rel diff to old 9.8e-15
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`proto.py`:

```python
# prototypes (monkeypatch only), nothing in the repo is edited
import numpy as np
from scipy import sparse
from dmrgpy.kondospectrumtk import stepfunctions as sf, conductance as cd, edkondo
from dmrgpy.algebra import algebra

def Theta_new(x):
    """Theta(x) = (1-t-xt)/(1-t)^2, t=exp(-|x|), mirrored for x<0: one expm1, no sinh/cosh, no masks"""
    x = np.asarray(x, dtype=float)
    a = np.abs(x)
    em1 = np.expm1(-a)
    t = 1.0 + em1
    with np.errstate(divide="ignore", invalid="ignore"):
        pos = (-em1 - a*t)/(em1*em1)
    pos = np.where(a < 1e-4, 0.5 + a/6.0 - a**3/360.0, pos)
    return np.where(x >= 0, pos, 1.0 - pos)

def F_call_new(self, x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    if self._F_spline is None: return self._direct(x)
    inside = (x >= self._F_min) & (x <= self._F_max)
    if inside.all(): return self._F_spline(x)          # F0 only where the table does not reach
    out = np.empty_like(x)
    out[inside] = self._F_spline(x[inside])
    xo = x[~inside]
    out[~inside] = sf.F0(xo, omega0=self.omega0, Gamma0=self.Gamma0)
    return out

def second_order_new(ks, eVs, T0=1.0, U=0.0, p_cut=cd.P_CUT, nbytes=4e6):
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    M2 = cd._spin_matrix_elements_squared(ks)[:, occ]          # [f,i in occ]
    eps = ks.e[:, None] - ks.e[None, occ]
    weight = ks.p[None, occ]*M2
    spin = np.empty(len(eVs))
    chunk = max(1, int(nbytes)//weight.size)
    for s0 in range(0, len(eVs), chunk):
        v = eVs[s0:s0+chunk]
        spin[s0:s0+len(v)] = (np.einsum('fi,efi->e', weight, cd._theta_raw(ks, v[:, None, None] - eps[None]))
                              + np.einsum('fi,efi->e', weight, cd._theta_raw(ks, -v[:, None, None] - eps[None])))
    return 2*np.pi*T0**2*(spin + 4*(U**2)*np.ones_like(eVs))

def eigenbasis_new(self, name):
    edobj = self.chain.get_ED_obj()
    op = getattr(self.chain, name)[self.site]
    mat = edobj.get_operator(op)
    vs = self.vs
    if sparse.issparse(mat):   # sparse @ dense first (O(dim^2)), then ONE dense product instead of two
        mat = mat.tocsr().astype(complex)
        return np.conjugate(vs.T) @ (mat @ vs)
    mat = np.array(algebra.todense(mat), dtype=complex)
    return np.conjugate(vs.T)@mat@vs

def F0_new(x, omega0=20e-3, Gamma0=5e-6):
    x = np.asarray(x, dtype=float)
    ax = np.abs(x)
    return np.log((omega0 + ax)/np.sqrt(x*x + Gamma0*Gamma0))

_old = dict(F0=sf.F0, Theta=sf.Theta, call=sf.FBuilder.__call__, second=cd.second_order_dIdV, eig=edkondo.KondoSpectrum._eigenbasis_spin_operator)
def enable(which):
    if "theta" in which: sf.Theta = Theta_new; cd.Theta = Theta_new
    if "F" in which: sf.FBuilder.__call__ = F_call_new
    if "F0" in which: sf.F0 = F0_new; cd.F0 = F0_new
    if "second" in which: cd.second_order_dIdV = second_order_new
    if "eig" in which: edkondo.KondoSpectrum._eigenbasis_spin_operator = eigenbasis_new
def disable():
    sf.F0 = cd.F0 = _old["F0"]; sf.Theta = cd.Theta = _old["Theta"]; sf.FBuilder.__call__ = _old["call"]
    cd.second_order_dIdV = _old["second"]; edkondo.KondoSpectrum._eigenbasis_spin_operator = _old["eig"]

def third_order_kondo_blocked(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6, Fb=None, p_cut=cd.P_CUT,
                              block_bytes=2.0e8, w_bytes=3.0e8):
    """same sums as conductance.third_order_kondo_dIdV; the coefficient tensor is built for a block of initial
    states at a time (block_bytes), and the exchange diagram is accumulated into W[e,f,m] for a group of biases"""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    dim = ks.dim
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    e_fm = ks.e[:, None] - ks.e[None, :]
    nB = max(1, int(block_bytes//(16*dim*dim)))
    nG = max(1, int(w_bytes//(8*dim*dim)))          # biases per group (W accumulator)
    fchunk = max(1, int(4e6)//(dim*dim))
    def one_direction(v):
        out = np.zeros(len(v))
        for g0 in range(0, len(v), nG):
            vg = v[g0:g0+nG]
            W = np.zeros((len(vg), dim, dim))
            for b0 in range(0, len(occ), nB):
                ob = occ[b0:b0+nB]; p = ks.p[ob]
                coeff = np.imag(cd._triple_product_coefficients(ks, ob))/2.
                eps = ks.e[None, :] - ks.e[ob, None]
                for s0 in range(0, len(vg), max(1, int(4e6)//(dim*len(ob)))):
                    vc = vg[s0:s0+max(1, int(4e6)//(dim*len(ob)))]
                    Th = cd._theta_raw(ks, vc[:, None, None] - eps[None])
                    Fd = Fcall((vc[:, None, None] - eps[None]).ravel()).reshape(len(vc), len(ob), dim)
                    out[g0+s0:g0+s0+len(vc)] += np.einsum('i,ifm,eif,eim->e', p, coeff, Th, Fd, optimize=True)
                    W[s0:s0+len(vc)] += np.einsum('i,ifm,eif->efm', p, coeff, Th, optimize=True)
            for s0 in range(0, len(vg), fchunk):
                vc = vg[s0:s0+fchunk]
                Fx = Fcall((vc[:, None, None] - e_fm[None]).ravel()).reshape(len(vc), dim, dim)
                out[g0+s0:g0+s0+len(vc)] += np.einsum('efm,efm->e', W[s0:s0+len(vc)], Fx)
        return out
    total = one_direction(eVs) + one_direction(-eVs)
    return 4*np.pi*T0**2*Jrho_s*total

_old["third"] = cd.third_order_kondo_dIdV
_enable_old = enable
def enable(which):
    _enable_old(which)
    if "third" in which: cd.third_order_kondo_dIdV = third_order_kondo_blocked
_disable_old = disable
def disable():
    _disable_old(); cd.third_order_kondo_dIdV = _old["third"]
```

`18_clean_timing.py`:

```python
# clean design: a fresh process per (config, repeat), interleaved, timing the whole public call INCLUDING eigh. ED, n=10 (1024 states), 41 biases, U=0.25, order 3.
import subprocess, sys, json, statistics
CHILD = r'''
import numpy as np, time, sys, json
from dmrgpy import spinchain
import proto
n=int(sys.argv[1]); T=float(sys.argv[2]); nb=int(sys.argv[3]); which=tuple(w for w in sys.argv[4].split(",") if w)
sc = spinchain.Spin_Chain(["1/2"]*n)
h = 0
for i in range(n-1):
    h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2]
sc.set_hamiltonian(h)
proto.enable(which)
eV = np.linspace(-0.01, 0.01, nb)
t0=time.perf_counter()
_, d = sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T, mode="ED")
print(json.dumps(dict(t=time.perf_counter()-t0, d=d.tolist())))
'''
open("_child18.py","w").write(CHILD)
import numpy as np
def run(n,T,nb,which):
    o=subprocess.run([sys.executable,"-s","_child18.py",str(n),str(T),str(nb),",".join(which)],capture_output=True,text=True)
    return json.loads(o.stdout.strip().splitlines()[-1])
def study(n,T,nb,configs,reps=3):
    ts={k:[] for k,_ in configs}; ref=None; diffs={}
    for r in range(reps):
        for k,w in configs:
            res=run(n,T,nb,w); ts[k].append(res["t"])
            if ref is None: ref=np.array(res["d"])
            if k not in diffs: diffs[k]=np.max(np.abs(np.array(res["d"])-ref))/np.max(np.abs(ref))
    base=statistics.median(ts[configs[0][0]])
    for k,_ in configs:
        m=statistics.median(ts[k])
        print("  n=%d T=%.0f K %d biases  %-9s median %5.2f s  (%s)  saves %4.1f%%  max rel diff to old %.1e"%(n,T,nb,k,m," ".join("%.2f"%t for t in ts[k]),100*(1-m/base),diffs[k]),flush=True)
cfg=[("old",()),("second",("second",)),("F",("F",)),("eig",("eig",)),("F0",("F0",)),("all",("second","F","eig","F0"))]
for T in (0.0,1.0): study(10,T,41,cfg)
print("blocked third-order against the unpatched repo code, high occupation")
study(9,20.0,21,[("old",()),("third",("third",))])
study(9,5.0,21,[("old",()),("third",("third",))])
```

Observed:

`18_clean_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`07_proto_accuracy.out`:

```
Theta vs 40-digit: max abs err old 8.17e-13 (at x=1.7e-04)  new 5.76e-13 (at x=-1.9e-04)
T=  0.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.0 F0                         max|new-old| = 1.78e-15   relative to peak 1.12e-16
T=  0.0 second                     max|new-old| = 3.55e-15   relative to peak 2.24e-16
T=  0.0 eig                        max|new-old| = 8.88e-16   relative to peak 5.59e-17
T=  0.0 third                      max|new-old| = 1.78e-15   relative to peak 1.12e-16
T=  0.0 F+F0+second+eig+third      max|new-old| = 5.33e-15   relative to peak 3.36e-16
T=  0.5 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.5 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  0.5 second                     max|new-old| = 5.33e-15   relative to peak 3.39e-16
T=  0.5 eig                        max|new-old| = 3.55e-15   relative to peak 2.26e-16
T=  0.5 third                      max|new-old| = 3.55e-15   relative to peak 2.26e-16
T=  0.5 F+F0+second+eig+third      max|new-old| = 7.11e-15   relative to peak 4.52e-16
T=  3.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  3.0 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T=  3.0 second                     max|new-old| = 5.10e-12   relative to peak 3.26e-13
T=  3.0 eig                        max|new-old| = 3.55e-15   relative to peak 2.27e-16
T=  3.0 third                      max|new-old| = 3.55e-15   relative to peak 2.27e-16
T=  3.0 F+F0+second+eig+third      max|new-old| = 5.10e-12   relative to peak 3.26e-13
T= 40.0 F                          max|new-old| = 0.00e+00   relative to peak 0.00e+00
T= 40.0 F0                         max|new-old| = 0.00e+00   relative to peak 0.00e+00
T= 40.0 second                     max|new-old| = 5.33e-14   relative to peak 3.52e-15
T= 40.0 eig                        max|new-old| = 5.33e-15   relative to peak 3.52e-16
T= 40.0 third                      max|new-old| = 3.55e-15   relative to peak 2.35e-16
T= 40.0 F+F0+second+eig+third      max|new-old| = 5.51e-14   relative to peak 3.64e-15
```

**Reviewer (CONFIRMED)**: Reproduced with my own prototype that keeps the dense fallback when get_operator does not return a sparse matrix. T=0 n=10: 5.94 to 5.07 s (14.6%); T=1 K: 14.00 to 13.32 s (4.9%); relative difference 1.3e-16 (real) and 1.2e-16 (complex H with Sy). Tested on a path the hunter did not use: a spinful Hubbard chain (4 sites, Nf=4 sector restriction, order 3, U=0.1) at T=0 and T=1 K gives reldiff 0.0 and 1.6e-16, so the sector-restricted operators go through the same code. Not documented as intended; kondo.md only says edkondo takes its spin matrices from get_operator, which is unchanged. Gain is modest and grows as dim^3 against the O(dim^2) sparse product, but eigh (dim^3) stays, so it never exceeds a fraction of the call.

Reviewer's own evidence:

```
cplx=False T=0 eig median 5.07 [5.07, 5.07, 5.1] reldiff 0.0e+00 (old 5.94)
cplx=False T=1 eig median 13.32 reldiff 1.3e-16 (old 14.00)
Hubbard Nf=4 sector: T=0 reldiff 0.0, T=1 reldiff 1.647e-16
```

Reviewer's measurement: own median old/new, n=10, 41 biases: T=0 5.94/5.07 s (14.6%), T=1 K 14.00/13.32 s (4.9%); max relative difference 1.3e-16. Does not move numbers beyond roundoff.

Reviewer's scripts and outputs:

`02_time.py`:

```python
import sys; sys.path.insert(0,"<repo>/src")
import numpy as np, time, statistics
from scipy import sparse
from dmrgpy import spinchain
from dmrgpy.algebra import algebra
from dmrgpy.kondospectrumtk import stepfunctions as sf, edkondo
def Fnew(self,x):
    x=np.atleast_1d(np.asarray(x,dtype=float))
    if self._F_spline is None: return self._direct(x)
    inside=(x>=self._F_min)&(x<=self._F_max)
    if inside.all(): return self._F_spline(x)
    out=np.empty_like(x); out[inside]=self._F_spline(x[inside]); out[~inside]=sf.F0(x[~inside],omega0=self.omega0,Gamma0=self.Gamma0); return out
def eignew(self,name):
    e=self.chain.get_ED_obj(); op=getattr(self.chain,name)[self.site]; mat=e.get_operator(op)
    if sparse.issparse(mat):
        return np.conjugate(self.vs.T)@(mat.tocsr().astype(complex)@self.vs)
    mat=np.array(algebra.todense(mat),dtype=complex); return np.conjugate(self.vs.T)@mat@self.vs
oF,oE=sf.FBuilder.__call__,edkondo.KondoSpectrum._eigenbasis_spin_operator
n=int(sys.argv[1])
def mk(cplx):
    sc=spinchain.Spin_Chain(["1/2"]*n); h=0
    for i in range(n-1): h=h+0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h=h+0.001*sc.Sz[0]+0.0003*sc.Sx[2]
    if cplx: h=h+0.0005*sc.Sy[1]
    sc.set_hamiltonian(h); return sc
eV=np.linspace(-0.01,0.01,41)
def run(sc,T):
    t=time.perf_counter(); r=sc.get_kondo_spectrum(eV,site=0,Jrho_s=-0.05,U=0.25,T=T,mode="ED")[1]; return time.perf_counter()-t,r
for cplx in (False,):
  sc=mk(cplx)
  for T in (0.,1.):
    res={}
    for rep in range(3):
        for k,(f,e) in {"old":(oF,oE),"F":(Fnew,oE),"eig":(oF,eignew)}.items():
            sf.FBuilder.__call__=f; edkondo.KondoSpectrum._eigenbasis_spin_operator=e
            t,r=run(sc,T); res.setdefault(k,[]).append((t,r))
    ref=res["old"][0][1]
    for k,v in res.items():
        m=statistics.median([a for a,_ in v]); d=max(np.max(abs(r-ref)) for _,r in v)/np.max(abs(ref))
        print("cplx=%s T=%g %-4s median %.2f %s reldiff %.1e"%(cplx,T,k,m,[round(a,2) for a,_ in v],d),flush=True)
```

`04_hub.py`:

```python
import sys; sys.path.insert(0,"<repo>/src")
import numpy as np
exec(open("02_time.py").read().split("oF,oE=")[0])
from dmrgpy import fermionchain
oE=edkondo.KondoSpectrum._eigenbasis_spin_operator
for cls in ("Spinful_Fermionic_Chain",):
    fc=getattr(fermionchain,cls)(4) if hasattr(fermionchain,cls) else None
    h=0
    for i in range(3):
        for s in ("up","dn"): pass
    h=sum(-0.01*(fc.Cdagup[i]*fc.Cup[i+1]+fc.Cdagdn[i]*fc.Cdn[i+1]+fc.Cdagup[i+1]*fc.Cup[i]+fc.Cdagdn[i+1]*fc.Cdn[i]) for i in range(3))+sum(0.02*fc.Nup[i]*fc.Ndn[i] for i in range(4))
    fc.set_hamiltonian(h); fc.set_conserved_sector(Nf=4)
    eV=np.linspace(-.01,.01,11)
    for T in (0.,1.):
        edkondo.KondoSpectrum._eigenbasis_spin_operator=oE
        a=fc.get_kondo_spectrum(eV,site=0,Jrho_s=-0.05,U=0.1,T=T,mode="ED")[1]
        edkondo.KondoSpectrum._eigenbasis_spin_operator=eignew
        b=fc.get_kondo_spectrum(eV,site=0,Jrho_s=-0.05,U=0.1,T=T,mode="ED")[1]
        print(T,np.max(abs(a-b))/np.max(abs(a)))
```

`03_time10.out`:

```
cplx=False T=0 old  median 5.94 [7.3, 5.93, 5.94] reldiff 0.0e+00
cplx=False T=0 F    median 5.93 [6.05, 5.93, 5.93] reldiff 0.0e+00
cplx=False T=0 eig  median 5.07 [5.07, 5.07, 5.1] reldiff 0.0e+00
cplx=False T=1 old  median 14.00 [13.98, 14.0, 18.65] reldiff 0.0e+00
cplx=False T=1 F    median 11.52 [11.51, 11.52, 11.71] reldiff 0.0e+00
cplx=False T=1 eig  median 13.32 [13.03, 13.32, 17.36] reldiff 1.3e-16
```

`04_hub.out`:

```
0.0 0.0
1.0 1.647298392733542e-16
```

**Suggested fix (hunter)**: Replace the dense triple product by np.conjugate(vs.T) @ (sparse_op.tocsr() @ vs), keeping algebra.todense only for an operator that is not sparse. Roundoff-level change.

**Fix (reviewer)**: The hunter's fix is right with one care: keep the dense branch for a non-sparse operator (np.array(algebra.todense(mat), dtype=complex)) and cast the csr to complex before the product. Replace edkondo.py:62-71 body accordingly. Existing tests tests/test_kondo_spectrum.py and tests/test_kondo_spectrum_spinful.py (sector case) cover it; add an assertion that Sx/Sy/Sz in the eigenbasis equal the old dense formula to 1e-12 on a spin chain and on a Nf-sector Hubbard chain.

### 7. Theta_prime loses its digits for 1e-4 < abs(x) < 1e-2 (0.0 for 1/6), reached by tabulate=False only

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `ed-reference` &middot; candidate `ed-reference-5`

**Status**: FIXED in `f6ac1e3` on 2026-11-12. `Theta_prime` is the closed form for 0.05 <= |x| <= 40 and a Taylor series below, within 1e-12 of a 50-digit value at every tested x. The FBuilder docstring now states its accuracy for omega0 >~ kT only, and says the table and the direct convolution are 2e-4 of F off for omega0 << kT; the wider table and finer u grid of the hunter's fix were not taken, since they move numbers in a regime outside the paper's. NUMBERS CHANGE, at T>0 only: 7.6e-13 relative on `mode="ED"`, 10 spins, 1 K, from table nodes between 1e-2 and 0.05 kT where the old form lost up to 4e-8. Pinned by `tests/test_audit_2026_11_stepfunctions.py::test_theta_prime_keeps_digits` (skipped without mpmath) and `::test_theta_prime_unit_normalized`.

The claim as reviewed: stepfunctions.Theta_prime loses its digits for 1e-4<|x|<1e-2 (returns 0.0 at x=1e-4 for an exact 1/6, 0.1645 at 3e-4, 0.16698 at 1e-3), which makes FBuilder(tabulate=False) 8.1e-3 off (7% of F=0.1156) at x/kT=1e-4, T=300 K, omega0=2 meV; no production path reaches the band (the table nodes and the u grid at the default npts=4001 have none in it: 0 nodes), so this is a latent defect of the validation path and of finer npts. Separately, the docstring's accuracy statements ('~1e-8 against the direct convolution') fail for omega0 much smaller than kT: table and direct are both 2.1e-4 off (0.2% of F) at x/kT=0.05, T=300 K, omega0=2 meV, because the u step does not resolve the ln(w0+|y|) kink; this is outside the paper's regimes (1 to 20 K, omega0 20 to 200 meV).

**Where**: src/dmrgpy/kondospectrumtk/stepfunctions.py:82-95 (Theta_prime mid branch), 125-160 (_F_convolve K_at and u grid), 163-210 (FBuilder docstring, tabulate=False, npts)

Theta_prime's mid branch is (x sinh(x)/2 - cosh(x) + 1)/(cosh(x)-1)^2, a difference of two numbers of size x^2/2 divided by (x^2/2)^2, and the series branch only takes over below 1e-4: the absolute error is about 4e-16/x^4, so Theta_prime(1e-4) = 0.0 (exact 0.16667), 1.1e-2 off at 3e-4, 3.1e-4 at 1e-3, 3e-8 at 1e-2 (script 05, 30-digit reference). Theta itself is fine (8e-13 absolute at worst, same as a one-expm1 form). Consequences, all anchored on a 25-digit quadrature: (a) _F_convolve evaluates K_at = Theta_prime(-x/kT)/kT at the requested point and the exact-subtraction trick (K-K_at)*L - K_at*int(L) only cancels K_at if Simpson integrates L exactly, which it does not (that is why the subtraction exists), so a garbage K_at leaks: FBuilder(...,tabulate=False) at T=300 K, omega0=2 meV: x/kT=1e-4 error +8.1e-03, 3e-4 +9.2e-05 (script 19), while the table is at -5.5e-06 there. tabulate=False is the documented way to validate the spline, and a random-x comparison in script 11 reported 8.0e-03 'table error' that is the direct path's. (b) Production escapes because no u node (spacing 0.03 kT at npts=4001) or table node (0.05 kT) falls in the band; finer npts do: npts=64001 gives +6.4e-07 and 256001 gives -1.1e-04 at F(0)=5.16 (script 13), against -3.5e-09 at the default. (c) Separately, the docstring's 'table reproduces the direct convolution to ~1e-8' and 'converged in npts to ~1e-6' hold at 1 K with omega0=20 meV (3.5e-9 to 4e-8 against mpmath, script 05) but not when omega0 << kT: at T=300 K, omega0=2 meV the table and direct are both 2.1e-4 off (0.2% of F=0.1156) at x/kT=0.05 because the 0.03 kT step does not resolve the ln(w0+|y|) kink (w0/kT=0.08), and at T=10 K, omega0=2 meV, x=-40.09 kT (just beyond the table, which ends at 40 kT when 40 kT > 10 omega0) the returned F0 is 2.2e-4 off (0.4% of F=0.056) since the thermal broadening there is not negligible (script 12). At Gamma0=5e-8 the table is 1e-5 off at x/kT=-2.19 (script 12). None of this reaches an existing test: tests/test_kondo_spectrum.py checks Theta_prime on a 21-point grid from -5 to 5 with atol 1e-4, and the paper figures use 1 K to 20 K with omega0 of 20 to 200 meV. Size in dI/dV: at or below 1e-4 relative for the default regime, so LOW.

The hunter's claim, before review: stepfunctions.Theta_prime loses all its digits for 1e-4 < |x| < 1e-2 (0.0 returned at x=1e-4 for an exact 1/6), which makes FBuilder(tabulate=False) off by 8.1e-3 (7% of F=0.1156) at x/kT=1e-4 and makes finer npts worse (F(0) error -3.5e-9 at npts=4001, +6.4e-7 at 64001, -1.1e-4 at 256001); the production table never lands in the band, and the FBuilder accuracy statements in its docstring measure 2e-4 to 2e-5 in some regimes.

Decisive output, quoted by the hunter:

```
Theta_prime(1e-4) = 0.0  (exact 1/6 = 0.16666666666666666 )
   x/kT     F(mpmath)   table err   tabulate=False err
    1e-04    0.115620   -5.5e-06    +8.1e-03
    3e-04    0.115620   -5.5e-06    +9.2e-05
  npts=  4001  error -3.5e-09
  npts= 64001  error +6.4e-07
  npts=256001  error -1.1e-04
  T= 10.000 w0=0.002 G0=5e-04 x/kT=-40.09  mpmath 0.05639355  table -2.20e-04  direct -2.03e-12 (errors vs mpmath)
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`05_F_mpmath.py`:

```python
# high-precision reference for F(x,T) and for Theta' accuracy
import numpy as np, mpmath as mp
from dmrgpy.kondospectrumtk.stepfunctions import Theta_prime, FBuilder, F0, _F_convolve
mp.mp.dps=30
kB=mp.mpf('8.617333262e-5')
def thp(t):  # exact Theta'(t) in mp
    t=mp.mpf(t)
    if abs(t)<mp.mpf('1e-6'): return mp.mpf(1)/6-t**2/60
    return (t*mp.sinh(t)/2-mp.cosh(t)+1)/(mp.cosh(t)-1)**2
xs=np.concatenate([np.logspace(-5,-1,9),[0.5,2,10,30,39.9]])
ex=np.array([float(thp(t)) for t in xs])
print("Theta_prime(x) vs 30-digit: x, current, exact, abs err")
for t,e in zip(xs,ex): print("  %9.2e %+.6e %+.6e %9.2e"%(t,Theta_prime(np.array([t]))[0],e,abs(Theta_prime(np.array([t]))[0]-e)))
def Fmp(x,T,w0,G0):
    kT=kB*mp.mpf(T); x=mp.mpf(x); w0=mp.mpf(w0); G0=mp.mpf(G0)
    f=lambda u: thp(u/kT)/kT*(mp.log(w0+abs(x+u))-mp.log(x*x+2*x*u+u*u+G0*G0)/2)
    W=60*kT
    pts=[-W,-x-100*G0,-x-10*G0,-x-G0,-x,-x+G0,-x+10*G0,-x+100*G0,W]
    pts=sorted(set(p for p in pts if -W<=p<=W))
    return mp.quad(f,pts)
print()
print("F(0,T): mpmath vs FBuilder(table) vs _F_convolve(direct 4001/8001 pts)")
for T,w0,G0 in [(1.0,0.02,5e-6),(1.0,0.02,5e-8),(0.1,0.02,5e-6),(5.0,0.02,5e-6),(30.,0.02,5e-6)]:
    ref=float(Fmp(0,T,w0,G0))
    fb=FBuilder(T,omega0=w0,Gamma0=G0)
    d4=_F_convolve(np.array([0.]),8.617333262e-5*T,w0,G0,npts=4001)[0]
    d8=_F_convolve(np.array([0.]),8.617333262e-5*T,w0,G0,npts=8001)[0]
    d16=_F_convolve(np.array([0.]),8.617333262e-5*T,w0,G0,npts=64001)[0]
    print("  T=%5.1f w0=%.3f G0=%.0e ref %.8f table %.8f (%.1e) direct4001 %.8f (%.1e) direct8001 (%.1e) direct64001 (%.1e)"%(T,w0,G0,ref,fb(np.array([0.]))[0],fb(np.array([0.]))[0]-ref,d4,d4-ref,d8-ref,d16-ref))
```

`12_table_vs_mpmath.py`:

```python
import numpy as np, mpmath as mp
from dmrgpy.kondospectrumtk.stepfunctions import FBuilder, _F_convolve
mp.mp.dps=25
kB=mp.mpf('8.617333262e-5')
def thp(t):
    t=mp.mpf(t)
    if abs(t)<mp.mpf('1e-6'): return mp.mpf(1)/6-t**2/60
    return (t*mp.sinh(t)/2-mp.cosh(t)+1)/(mp.cosh(t)-1)**2
def Fmp(x,T,w0,G0):
    kT=kB*mp.mpf(T); x=mp.mpf(x); w0=mp.mpf(w0); G0=mp.mpf(G0)
    f=lambda u: thp(u/kT)/kT*(mp.log(w0+abs(x+u))-mp.log((x+u)**2+G0*G0)/2)
    W=60*kT
    pts=[-W,-x-100*G0,-x-10*G0,-x-G0,-x,-x+G0,-x+10*G0,-x+100*G0,-x-w0,-x+w0,-x-5*w0,-x+5*w0,W]
    pts=sorted(set(p for p in pts if -W<=p<=W))
    return mp.quad(f,pts,maxdegree=8)
cases=[(300.,2e-3,5e-6,[0.0,0.5,2.0,-1.0]),(10.,2e-3,5e-4,[-40.09,0.0,5.0]),(1.0,0.2,5e-8,[-2.19,0.0,3.0])]
for T,w0,G0,xk in cases:
    kT=8.617333262e-5*T; fb=FBuilder(T,omega0=w0,Gamma0=G0)
    for k in xk:
        x=k*kT
        ref=float(Fmp(x,T,w0,G0)); tab=fb(np.array([x]))[0]; dr=fb._direct(np.array([x]))[0]
        print("T=%6.1f w0=%.3f G0=%.0e x/kT=%+7.2f  mpmath %.8f  table %+.2e  direct %+.2e (errors vs mpmath)"%(T,w0,G0,k,ref,tab-ref,dr-ref),flush=True)
```

`13_direct_theta_prime_cancellation.py`:

```python
# FBuilder(tabulate=False) and _F_convolve at |x| ~ 1e-4..1e-2 kT: Theta_prime cancellation enters the kernel value K_at
import numpy as np, mpmath as mp
from dmrgpy.kondospectrumtk.stepfunctions import FBuilder, Theta_prime
mp.mp.dps=25
kB=mp.mpf('8.617333262e-5')
def thp(t):
    t=mp.mpf(t)
    if abs(t)<mp.mpf('1e-6'): return mp.mpf(1)/6-t**2/60
    return (t*mp.sinh(t)/2-mp.cosh(t)+1)/(mp.cosh(t)-1)**2
def Fmp(x,T,w0,G0):
    kT=kB*mp.mpf(T); x=mp.mpf(x); w0=mp.mpf(w0); G0=mp.mpf(G0)
    f=lambda u: thp(u/kT)/kT*(mp.log(w0+abs(x+u))-mp.log((x+u)**2+G0*G0)/2)
    W=60*kT
    pts=sorted(set(p for p in [-W,-x-100*G0,-x-10*G0,-x-G0,-x,-x+G0,-x+10*G0,-x+100*G0,-x-w0,-x+w0,W] if -W<=p<=W))
    return mp.quad(f,pts,maxdegree=8)
T,w0,G0=1.0,20e-3,5e-6
kT=8.617333262e-5*T
fb=FBuilder(T,omega0=w0,Gamma0=G0); fd=FBuilder(T,omega0=w0,Gamma0=G0,tabulate=False)
print("Theta_prime(1e-4) =",Theta_prime(np.array([1e-4]))[0]," (exact 1/6 =",1/6,")")
print("T=1 K, omega0=20 meV, Gamma0=5 ueV; error against 25-digit quadrature")
print("   x/kT      F(mpmath)    table err    tabulate=False err")
for k in (1e-4,3e-4,1e-3,3e-3,1e-2,1e-1,1.0):
    x=k*kT; ref=float(Fmp(x,T,w0,G0))
    print("  %7.0e  %10.6f  %+9.1e   %+9.1e"%(k,ref,fb(np.array([x]))[0]-ref,fd(np.array([x]))[0]-ref),flush=True)
print("effect of npts on F(0) (default 4001):")
from dmrgpy.kondospectrumtk.stepfunctions import _F_convolve
ref0=float(Fmp(0,T,w0,G0))
for npts in (1001,4001,8001,16001,64001,256001):
    print("  npts=%6d  error %+.1e"%(npts,_F_convolve(np.array([0.]),kT,w0,G0,npts=npts)[0]-ref0))
```

`19_theta_prime_band_anchor.py`:

```python
# which side is wrong in the Theta_prime cancellation band: 25-digit quadrature at 300 K, omega0=2 meV, Gamma0=5 ueV
import numpy as np, mpmath as mp
from dmrgpy.kondospectrumtk.stepfunctions import FBuilder
mp.mp.dps=25
kB=mp.mpf('8.617333262e-5')
def thp(t):
    t=mp.mpf(t)
    if abs(t)<mp.mpf('1e-6'): return mp.mpf(1)/6-t**2/60
    return (t*mp.sinh(t)/2-mp.cosh(t)+1)/(mp.cosh(t)-1)**2
def Fmp(x,T,w0,G0):
    kT=kB*mp.mpf(T); x=mp.mpf(x); w0=mp.mpf(w0); G0=mp.mpf(G0)
    f=lambda u: thp(u/kT)/kT*(mp.log(w0+abs(x+u))-mp.log((x+u)**2+G0*G0)/2)
    W=60*kT
    pts=sorted(set(p for p in [-W,-x-100*G0,-x-10*G0,-x-G0,-x,-x+G0,-x+10*G0,-x+100*G0,-x-w0,-x+w0,-x-5*w0,-x+5*w0,W] if -W<=p<=W))
    return mp.quad(f,pts,maxdegree=8)
T,w0,G0=300.,2e-3,5e-6
kT=8.617333262e-5*T
fb=FBuilder(T,omega0=w0,Gamma0=G0); fd=FBuilder(T,omega0=w0,Gamma0=G0,tabulate=False)
print("T=300 K, omega0=2 meV, Gamma0=5 ueV; errors against 25-digit quadrature")
print("   x/kT     F(mpmath)   table err   tabulate=False err")
for k in (1e-4,3e-4,1e-3,3e-3,1e-2,0.05):
    x=k*kT; ref=float(Fmp(x,T,w0,G0))
    print("  %7.0e  %9.6f  %+9.1e   %+9.1e"%(k,ref,fb(np.array([x]))[0]-ref,fd(np.array([x]))[0]-ref),flush=True)
```

`11_table_and_grids.py`:

```python
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.stepfunctions import FBuilder, F0
kB=8.617333262e-5
rng=np.random.default_rng(1)
print("table (spline) vs direct convolution at 3000 random x, max abs / max rel-to-F")
for T,w0,G0 in [(1e-3,0.02,5e-6),(0.1,0.02,5e-6),(1.0,0.02,5e-6),(1.0,0.2,5e-8),(10.,0.002,5e-4),(300.,0.002,5e-6),(300.,0.2,5e-6)]:
    kT=kB*T; fb=FBuilder(T,omega0=w0,Gamma0=G0)
    top=max(40*kT,10*w0)
    x=np.concatenate([rng.uniform(-45*kT,45*kT,1000), rng.choice([-1,1],1000)*np.exp(rng.uniform(np.log(kT),np.log(top*1.5),1000)), rng.uniform(-3*w0,3*w0,1000)])
    a=fb(x); d=fb._direct(x)
    j=np.argmax(np.abs(a-d))
    print("  T=%7.3f w0=%.3f G0=%.0e  max|tab-dir| %.2e at x/kT=%+.2f ; F range [%.2f, %.2f]"%(T,w0,G0,np.abs(a-d).max(),x[j]/kT,d.min(),d.max()))
# grids at the public API
n=9
sc=spinchain.Spin_Chain(["1/2"]*n)
h=0
for i in range(n-1): h=h+0.004*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h+0.001*sc.Sz[0]+0.0003*sc.Sx[2])
for T in (0.0,1.0):
    eV=np.linspace(-0.01,0.01,31)
    _,ref=sc.get_kondo_spectrum(eV,site=0,Jrho_s=-0.05,U=0.25,T=T,mode="ED")
    perm=rng.permutation(len(eV)); eVp=np.concatenate([eV[perm],eV[:5]])
    _,got=sc.get_kondo_spectrum(eVp,site=0,Jrho_s=-0.05,U=0.25,T=T,mode="ED")
    print("T=%.0f shuffled + duplicated grid: max|got-ref[perm]| = %.2e"%(T,max(np.abs(got[:31]-ref[perm]).max(),np.abs(got[31:]-ref[:5]).max())))
    _,one=sc.get_kondo_spectrum(eV[7:8],site=0,Jrho_s=-0.05,U=0.25,T=T,mode="ED")
    print("   one bias point: %.12g vs %.12g"%(one[0],ref[7]))
    try:
        r=sc.get_kondo_spectrum(0.003,site=0,Jrho_s=-0.05,U=0.25,T=T,mode="ED"); print("   scalar eV ->",r)
    except Exception as e: print("   scalar eV raises",type(e).__name__,e)
```

Observed:

`05_F_mpmath.out`:

```
Theta_prime(x) vs 30-digit: x, current, exact, abs err
   1.00e-05 +1.666667e-01 +1.666667e-01  3.11e-11
   3.16e-05 +1.666667e-01 +1.666667e-01  9.73e-14
   1.00e-04 +0.000000e+00 +1.666667e-01  1.67e-01
   3.16e-04 +1.776357e-01 +1.666667e-01  1.10e-02
   1.00e-03 +1.669775e-01 +1.666667e-01  3.11e-04
   3.16e-03 +1.666664e-01 +1.666665e-01  9.73e-08
   1.00e-02 +1.666650e-01 +1.666650e-01  3.15e-08
   3.16e-02 +1.666500e-01 +1.666500e-01  2.73e-10
   1.00e-01 +1.665001e-01 +1.665001e-01  2.39e-12
   5.00e-01 +1.625613e-01 +1.625613e-01  7.63e-15
   2.00e+00 +1.133284e-01 +1.133284e-01  5.55e-17
   1.00e+01 +3.632736e-04 +3.632736e-04  1.63e-19
   3.00e+01 +2.620134e-12 +2.620134e-12  4.04e-28
   3.99e+01 +1.779465e-16 +1.779465e-16  2.47e-32

F(0,T): mpmath vs FBuilder(table) vs _F_convolve(direct 4001/8001 pts)
  T=  1.0 w0=0.020 G0=5e-06 ref 5.16494600 table 5.16494599 (-3.5e-09) direct4001 5.16494599 (-3.5e-09) direct8001 (-1.2e-09) direct64001 (6.4e-07)
  T=  1.0 w0=0.020 G0=5e-08 ref 5.19471863 table 5.19471859 (-3.7e-08) direct4001 5.19471859 (-3.7e-08) direct8001 (-6.0e-09) direct64001 (9.1e-07)
  T=  0.1 w0=0.020 G0=5e-06 ref 7.21378924 table 7.21378924 (-1.3e-10) direct4001 7.21378924 (-1.3e-10) direct8001 (1.9e-10) direct64001 (6.4e-07)
  T=  5.0 w0=0.020 G0=5e-06 ref 3.61261200 table 3.61261196 (-4.1e-08) direct4001 3.61261196 (-4.1e-08) direct8001 (-3.8e-09) direct64001 (6.8e-07)
  T= 30.0 w0=0.020 G0=5e-06 ref 2.00182983 table 2.00182979 (-3.7e-08) direct4001 2.00182979 (-3.7e-08) direct8001 (-5.7e-09) direct64001 (6.1e-07)
```

`12_table_vs_mpmath.out`:

```
T= 300.0 w0=0.002 G0=5e-06 x/kT=  +0.00  mpmath 0.11561988  table -5.50e-06  direct -5.50e-06 (errors vs mpmath)
T= 300.0 w0=0.002 G0=5e-06 x/kT=  +0.50  mpmath 0.11344553  table -1.08e-06  direct -1.08e-06 (errors vs mpmath)
T= 300.0 w0=0.002 G0=5e-06 x/kT=  +2.00  mpmath 0.08668680  table -7.73e-06  direct -7.73e-06 (errors vs mpmath)
T= 300.0 w0=0.002 G0=5e-06 x/kT=  -1.00  mpmath 0.10725208  table +9.06e-06  direct +9.06e-06 (errors vs mpmath)
T=  10.0 w0=0.002 G0=5e-04 x/kT= -40.09  mpmath 0.05639355  table -2.20e-04  direct -2.03e-12 (errors vs mpmath)
T=  10.0 w0=0.002 G0=5e-04 x/kT=  +0.00  mpmath 0.86572325  table +9.75e-11  direct +9.75e-11 (errors vs mpmath)
T=  10.0 w0=0.002 G0=5e-04 x/kT=  +5.00  mpmath 0.45034522  table +1.70e-08  direct +1.70e-08 (errors vs mpmath)
T=   1.0 w0=0.200 G0=5e-08 x/kT=  -2.19  mpmath 7.11751136  table +1.03e-05  direct -1.57e-08 (errors vs mpmath)
T=   1.0 w0=0.200 G0=5e-08 x/kT=  +0.00  mpmath 7.48960783  table -3.65e-08  direct -3.65e-08 (errors vs mpmath)
T=   1.0 w0=0.200 G0=5e-08 x/kT=  +3.00  mpmath 6.87308651  table +1.12e-08  direct +1.12e-08 (errors vs mpmath)
```

`13_direct_theta_prime_cancellation.out`:

```
Theta_prime(1e-4) = 0.0  (exact 1/6 = 0.16666666666666666 )
T=1 K, omega0=20 meV, Gamma0=5 ueV; error against 25-digit quadrature
   x/kT      F(mpmath)    table err    tabulate=False err
    1e-04    5.164946   -3.5e-09    -7.6e-06
    3e-04    5.164946   -3.5e-09    -9.9e-08
    1e-03    5.164946   -3.4e-09    +2.4e-08
    3e-03    5.164945   -3.0e-09    +3.3e-08
    1e-02    5.164937   +1.9e-09    +6.9e-08
    1e-01    5.164070   +1.1e-08    +1.1e-08
    1e+00    5.080318   +8.1e-08    +8.1e-08
effect of npts on F(0) (default 4001):
  npts=  1001  error -2.5e-06
  npts=  4001  error -3.5e-09
  npts=  8001  error -1.2e-09
  npts= 16001  error -7.0e-10
  npts= 64001  error +6.4e-07
  npts=256001  error -1.1e-04
```

`19_theta_prime_band_anchor.out`:

```
T=300 K, omega0=2 meV, Gamma0=5 ueV; errors against 25-digit quadrature
   x/kT     F(mpmath)   table err   tabulate=False err
    1e-04   0.115620   -5.5e-06    +8.1e-03
    3e-04   0.115620   -5.5e-06    +9.2e-05
    1e-03   0.115620   -5.3e-06    +3.0e-05
    3e-03   0.115620   -3.9e-06    +1.0e-04
    1e-02   0.115619   +1.1e-05    +2.1e-04
    5e-02   0.115598   +2.1e-04    +2.1e-04
```

`11_table_and_grids.out`:

```
table (spline) vs direct convolution at 3000 random x, max abs / max rel-to-F
  T=  0.001 w0=0.020 G0=5e-06  max|tab-dir| 2.68e-08 at x/kT=-2308315.82 ; F range [0.07, 8.29]
  T=  0.100 w0=0.020 G0=5e-06  max|tab-dir| 2.28e-08 at x/kT=+0.57 ; F range [0.06, 7.21]
  T=  1.000 w0=0.020 G0=5e-06  max|tab-dir| 3.41e-07 at x/kT=+1.07 ; F range [0.06, 5.16]
  T=  1.000 w0=0.200 G0=5e-08  max|tab-dir| 2.05e-05 at x/kT=-2.19 ; F range [0.06, 7.49]
  T= 10.000 w0=0.002 G0=5e-04  max|tab-dir| 2.20e-04 at x/kT=-40.09 ; F range [0.04, 0.87]
  T=300.000 w0=0.002 G0=5e-06  max|tab-dir| 8.00e-03 at x/kT=+0.00 ; F range [0.00, 0.12]
  T=300.000 w0=0.200 G0=5e-06  max|tab-dir| 9.51e-05 at x/kT=+77.55 ; F range [0.06, 2.00]
T=0 shuffled + duplicated grid: max|got-ref[perm]| = 0.00e+00
   one bias point: 6.44646347568 vs 6.44646347568
   scalar eV raises IndexError too many indices for array: array is 0-dimensional, but 1 were indexed
T=1 shuffled + duplicated grid: max|got-ref[perm]| = 4.44e-16
   one bias point: 6.43699549707 vs 6.43699549707
   scalar eV raises IndexError too many indices for array: array is 0-dimensional, but 1 were indexed
```

**Reviewer (CONFIRMED (narrowed))**: Theta_prime values reproduced against a 40-digit mpmath reference: 0.0 vs 0.1666667 at 1e-4, 0.16448 at 3e-4, 0.166978 at 1e-3, 3e-8 off at 1e-2. Direct path at T=300 K, omega0=2 meV, x=1e-4 kT gives 0.12370 against table 0.11561, matching the hunter's 8.1e-3. Production: at T=1 K the spline grid has 0 nodes with |x|/kT in (1e-4,1e-2) and the u grid (spacing 0.03 kT, node at 0) has 0, so the cancellation never enters the shipped numbers. I re-measured sub-claim (c) against mpmath quadrature: error 2.1e-4 for both table and direct at x/kT=0.05, shrinking to 1e-6 at x/kT=0.5 and 8e-6 at 2, so the kink-resolution error is confined to |x| below omega0 when omega0<kT, and is real. Sub-claims not re-measured: npts=64001/256001 error growth (b) and the beyond-table 2.2e-4 case at T=10 K, omega0=2 meV, Gamma0=5e-8 (I accept them as plausible from the mechanism but did not reproduce them). The size of the defect in dI/dV is at most 2e-4 relative to F in an unphysical regime for this method, and zero in the tested regime. It is not documented as a known boundary. The docstring statement should be narrowed (valid for omega0 >~ kT), which is the cheapest fix.

Reviewer's own evidence:

```
Theta_prime(1e-4)=0.0 exact 0.1666666665; (3e-4) 0.16447748 vs 0.16666667; (1e-3) 0.16697752 vs 0.16666665
nodes in band: 0 ; u in band 0
direct [0.12370187] table [0.11561438] (T=300, omega0=2e-3, x=1e-4 kT)
0.05 0.11559785760071856 table 0.00021131 direct 0.00021131
0.5 table -1.08e-06 direct -1.08e-06
```

Reviewer's scripts and outputs:

`01_tp.py`:

```python
import sys; sys.path.insert(0,"<repo>/src")
import numpy as np, time
from dmrgpy.kondospectrumtk import stepfunctions as sf
for x in [1e-4,3e-4,1e-3,1e-2]:
    ex=(x*np.sinh(x)/2-np.cosh(x)+1)  # naive
    # exact via series/ high precision
    import mpmath as mp; mp.mp.dps=40; X=mp.mpf(x)
    e=(X*mp.sinh(X)/2-mp.cosh(X)+1)/(mp.cosh(X)-1)**2
    print(x, sf.Theta_prime(np.array([x]))[0], float(e))
# production table nodes vs band
T=1.0;kT=8.617333262e-5*T
fb=sf.FBuilder(T); print("nodes in band:", np.sum((np.abs(fb._F_spline.x/kT)>1e-4)&(np.abs(fb._F_spline.x/kT)<1e-2)))
u=np.linspace(-60*kT,60*kT,4001)/kT; print("u in band",np.sum((abs(u)>1e-4)&(abs(u)<1e-2)))
# ed-ref-2 timing
x=np.random.default_rng(0).uniform(-.05,.05,1_000_000)
t=time.time();a=fb(x);print("call",time.time()-t)
t=time.time();F0=sf.F0(x);print("F0 only",time.time()-t)
# direct F accuracy at x/kT=1e-4 
fd=sf.FBuilder(300,omega0=2e-3,tabulate=False); k=8.617333262e-5*300
print("direct",fd(np.array([1e-4*k])), sf.FBuilder(300,omega0=2e-3)(np.array([1e-4*k])))
```

`06_acc.py`:

```python
import sys; sys.path.insert(0,"<repo>/src")
import numpy as np, mpmath as mp
from dmrgpy.kondospectrumtk import stepfunctions as sf
mp.mp.dps=30
kB=8.617333262e-5; T=300.; kT=kB*T; w0=2e-3; G=5e-6
def Kp(v): # Theta'(v)
    v=mp.mpf(v); return (v*mp.sinh(v)/2-mp.cosh(v)+1)/(mp.cosh(v)-1)**2 if v!=0 else mp.mpf(1)/6
def Fex(x):
    f=lambda u: Kp(u/kT)/kT*(mp.log(w0+abs(x+u))-mp.log(mp.mpf(x+u)**2+G**2)/2)
    pts=[-60*kT,-x-1e-3*kT if False else -x,60*kT]
    pts=sorted(set([-60*kT,-w0-x,w0-x,-x,60*kT, -x-G,-x+G]))
    pts=[p for p in pts if abs(p)<=60*kT]
    return mp.quad(f,pts,maxdegree=10)
fb=sf.FBuilder(T,omega0=w0,Gamma0=G); fd=sf.FBuilder(T,omega0=w0,Gamma0=G,tabulate=False)
for r in (0.05,0.5,2.0):
    x=r*kT; e=float(Fex(x)); print(r,e,"table",fb(np.array([x]))[0]-e,"direct",fd(np.array([x]))[0]-e)
```

`01_tp.out`:

```
0.0001 0.0 0.1666666665
0.0003 0.16447748225935538 0.16666666516666667
0.001 0.16697751512375347 0.166666650000001
0.01 0.16666496849996815 0.1666650000099206
nodes in band: 0
u in band 0
call 0.23351502418518066
F0 only 0.03144669532775879
direct [0.12370187] [0.11561438]
```

`06_acc.out`:

```
0.05 0.11559785760071856 table 0.0002113105188499692 direct 0.0002113105188499831
0.5 0.11344552645096022 table -1.0795250156431235e-06 direct -1.0795250174194804e-06
2.0 0.08668679724971438 table -7.725061330884908e-06 direct -7.72506133005224e-06
```

**Suggested fix (hunter)**: Write Theta_prime as [x cosh(x/2) - 2 sinh(x/2)] / (4 sinh^3(x/2)) with a Taylor series for |x| < 0.05 (or use the series to |x| = 0.1, error 1e-12 there), so the kernel is accurate everywhere; evaluate K_at from it. Make the table extend to where the thermal correction is below 1e-8 (e.g. at least 100 kT) or blend into F0 only beyond that, and state the accuracy of FBuilder as a function of omega0/kT (resolve the kink by refining u near 0 when omega0 < kT). Changes F by under 1e-5 at default parameters and up to 2e-4 for omega0 << kT; numbers change only there.

**Fix (reviewer)**: The hunter's Theta_prime rewrite is right: use [x cosh(x/2) - 2 sinh(x/2)] / (4 sinh^3(x/2)) (or the 1/6 - x^2/60 + x^4/... series up to |x|~0.05-0.1) in stepfunctions.py Theta_prime; this changes no production number (no node in the band) and fixes tabulate=False, which is the validation path. Pin with a test comparing Theta_prime to a mpmath-computed value at x in {1e-4, 3e-4, 1e-3, 1e-2} to 1e-10 (tests/test_kondo_spectrum.py currently checks only a 21-point grid with atol 1e-4). Do NOT adopt the hunter's wider changes (table out to 100 kT, refining u for omega0<kT) in the same pass: they move numbers and are for an out-of-regime case. For (c), narrow the FBuilder docstring accuracy sentence to omega0 >~ kT instead of changing the grid.

### 8. The leading-log flow runs in complex arithmetic for a real Hamiltonian: 2x to 2.5x of the call

`performance` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `leading-log` &middot; candidate `leading-log-1`

**Status**: FIXED in `fce3382` on 2026-11-12. `_joint_flow` keeps the state in real block form when the block form of the bare vertex is real to `REAL_TOL` = 1e-12, and runs the old complex code otherwise; the `Gamma0 > 0` check is kept, and the step control is unchanged. With finding 9, on d22, one core, 61 biases: 6 spins 2.80 s to 0.99 s, 7 spins 13.52 s to 4.04 s; a field along y takes the complex route, 2.83 s to 2.41 s (the table only). NUMBERS CHANGE: 1.2e-12 (6 spins) and 2.7e-11 (7 spins) relative, from the integrator's step sequence (614 to 518 right-hand-side evaluations at 7 spins), inside `rtol=1e-10`; roundoff on the complex route; `_TWO_PASS` at 1e-9 holds. Pinned by `tests/test_audit_2026_11_leadinglog.py` (`::test_the_block_form_is_real_for_a_real_hamiltonian_only`, `::test_real_and_complex_routes_agree`, `::test_a_field_along_y_takes_the_complex_route_and_gives_the_rotated_values`, `::test_the_real_route_gives_the_flowed_vertex_of_the_complex_one`, `::test_gamma0_must_be_positive_on_the_real_route`, `::test_a_stopped_flow_on_the_real_route_is_the_complex_one`).

The claim as reviewed: For a real Hamiltonian the leading-log flow can run in real block arithmetic, which makes get_kondo_spectrum(order="LL", T=0) 1.98x faster at 6 spins (2.79 s to 1.41 s) and 2.46x at 7 spins (13.57 s to 5.52 s) on d22, with values moving by 1e-12 to 6e-11 as integrator step-sequence movement inside rtol=1e-10, not roundoff.

**Where**: src/dmrgpy/kondospectrumtk/leadinglog.py:312-387 (_joint_flow: state kept as (4,d,d) A-matrices, _block/_unblock at every rhs call, complex y), :133-138 (_bracket), :190-217 (_tip_bracket splits P and Q into real and imaginary halves, so a real problem pays the 8m-column hole sum and two matmuls)

WHAT THE CODE DOES. _joint_flow stores G, A_on and the tip vertex as complex (4,d,d) matrices, converts to block form _block() and back _unblock() on every right-hand side, and multiplies in complex. The block matrix [[A0+A3, A1-iA2],[A1+iA2, A0-A3]] is real whenever Sx and Sz are real and Sy is imaginary in the eigenbasis, which is the case for every real Hamiltonian (08_realness.py: KondoSpectrum on a tilted-field chain and on the degenerate B=0 Heisenberg chain gives max|Im Sx|, max|Im Sz|, max|Re Sy| all exactly 0.0; with a By field they are 0.3). The kernel is real, the products and _swap keep a real block real, so the whole flow stays real: the four 2d x 2d complex products become real products (a zgemm is 4x a dgemm), P.imag and Q.imag vanish so the hole sum has 4m columns instead of 8m, and the state vector is half the size. The prototype (llreal.py) also keeps the state in block form, so _block/_unblock disappear from the right-hand side; it falls back to the current complex code when the bare block has an imaginary part above 1e-12 relative.

WHAT IS EXPECTED. Same equations. EVIDENCE THE ALGEBRA IS IDENTICAL: one right-hand side evaluated by the current code and by the prototype on the same mid-flow state (7 spins, state after 25 solver steps) differs by 2.14e-15 relative (16_equiv_timing.out); the imaginary part of the old rhs in block form is exactly 0.0.

WHAT MOVES. The solver's error norm sees a different state layout, so the step sequence changes slightly (41 to 38 accepted steps, 614 to 518 rhs evaluations at 7 spins): about 16% of the gain is fewer evaluations and 2.6/1.16 = 2.2x is per-evaluation. The values move by 2.7e-11 relative at 7 spins and 1.2e-12 at 6 spins, which is step-sequence movement inside the integration tolerance. Against an rtol=1e-12 reference at 7 spins the old code is 3.06e-11 off and the prototype 9.77e-12 off (6 spins: both 1.51e-11), so the new numbers are not worse. By the rule of this audit that makes it numbers_change=true, at the 1e-11 level; if a byte-for-byte roundoff-only change is wanted, keep the current layout of y and make only the inside of rhs real (not measured separately).

SIZES. Stopped flows (S=1/2 in a field at Jrho_s=-0.2, dimer in a field at -0.2): same NaN mask, D_stop equal to 1e-10 relative, values 5e-11 to 1.2e-10. S=1 anisotropic at -0.15 and the 4 and 5 spin chains: 3e-11 and 1.4e-11. Occupied manifold of two states (B=0 trimer, real H, doublet ground state): 7.9e-11 on dI/dV, 3.7e-12 on info["coupling"], and tip_vertex(cols=None) agrees to 1.9e-11 (17_degenerate_real_and_probes.out). The complex fallback (By and Bz field dimer, complex H) agrees to 6e-12 (09_real_check.out).

TIMINGS (median of 3, d22, one core, 61 biases, pinned threads, 16_equiv_timing.out for 7 spins, 10_timing.out for 6 spins; 10_timing.out was produced before llreal got the hole precompute, so it is the real+block number): 6 spins old 3.05 s, new 1.49 s (2.05x); 7 spins old 14.30 s, new 5.50 s (2.60x). Through the public get_kondo_spectrum call at 7 spins with both changes (19_public_profile.out): 14.70 s to 4.25 s, 3.46x. The ED part of the call is cached and is 0.03% of the time.

WHY NO TEST SEES IT. Nothing is wrong: the tests pin values to rel=1e-9 (_TWO_PASS) and every model there is far below the cost where this matters; the real structure is invisible to them.

The hunter's claim, before review: The leading-log flow runs in complex arithmetic although, for a real Hamiltonian, every vertex is a real 2d x 2d matrix in block form: keeping the state in block form and in real arithmetic makes get_kondo_spectrum(order="LL") 2.05x faster at 6 spins (3.05 s to 1.49 s) and 2.60x at 7 spins (14.30 s to 5.50 s), a saving of 51 to 62% of the whole public call, with values moving by 1e-12 to 3e-11.

Decisive output, quoted by the hunter:

```
rhs, USE_PRE=False: max |new-old| / max|old| = 2.14e-15
n=7 d=128 old                  median 14.30 s (runs 13.95, 14.48, 14.30)  steps 41  rhs evals 614
n=7 d=128 real+block           median 5.50 s (runs 5.56, 5.47, 5.50)  steps 38  rhs evals 518
  real+block           speedup 2.60x  |x-old|max rel 2.70e-11   |x-ref1e-12|max rel 9.77e-12   (old vs ref1e-12: 3.06e-11)
n=6 d=64  old median 3.05 s ; real+block median 1.49 s (2.05x), max rel diff vs old 1.17e-12
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`llreal.py`:

```python
"""Prototype (not in the repo): leading-log flow kept in block form, in real arithmetic whenever the
vertices are real in block form (real H and Sx, Sz real, Sy imaginary in the eigenbasis, the usual case).
Same equations as leadinglog._joint_flow; monkeypatched in with install()."""
import numpy as np
from scipy.integrate import DOP853
from scipy.optimize import brentq
from dmrgpy.kondospectrumtk import leadinglog as L
import llhole
USE_PRE = True

def block_is_real(S, tol=1e-12):
    """block form of the bare vertex: [[A0+A3, A1-iA2],[A1+iA2, A0-A3]] real?"""
    d = S.shape[-1]
    B = L._block(np.concatenate([np.zeros((1, d, d), complex), S]))
    return np.max(np.abs(B.imag)) <= tol*max(np.max(np.abs(B)), 1e-300)

def tip_bracket_b(Gb, Ab, cols, w, e, D, Dlow, omega0, Gamma0, real, pre=None):
    d, m, n = len(e), len(cols), len(w)
    G4 = Gb.reshape(2, d, 2, d); A4 = Ab.reshape(2, d, 2, d)
    P = np.einsum('iakc,kcnj->jianc', G4, A4[:, :, :, cols]).reshape(m, 4*d, d)
    x = w[None, None, :] - (e[None, :, None] - e[cols][:, None, None])
    Kp = L._kernel(x, D, omega0, Gamma0) if Dlow < 0. else L._kernel_integral(x, Dlow, D, omega0, Gamma0)
    Q = np.einsum('kanc,ickj->ainjc', A4, G4[:, :, :, cols]).reshape(d, 4*m, d)
    if real:
        part = np.matmul(P, Kp)
        hole = np.empty((n, d, 4*m))
        if pre is not None and Dlow < 0.: llhole.hole_sum_pre(np.ascontiguousarray(Q), pre[0], pre[1], D, hole)
        else: L._hole_sum(np.ascontiguousarray(Q), w, e, D, Dlow, omega0, Gamma0, hole)
        hole = hole.reshape(n, d, 2, 2, m)
    else:
        part = (np.matmul(np.ascontiguousarray(P.real), Kp) + 1j*np.matmul(np.ascontiguousarray(P.imag), Kp))
        Qr = np.ascontiguousarray(np.concatenate([Q.real, Q.imag], axis=1))
        hole = np.empty((n, d, 8*m))
        L._hole_sum(Qr, w, e, D, Dlow, omega0, Gamma0, hole)
        hole = (hole[:, :, :4*m] + 1j*hole[:, :, 4*m:]).reshape(n, d, 2, 2, m)
    part = part.reshape(m, 2, d, 2, n).transpose(4, 1, 2, 3, 0)
    return part - hole.transpose(0, 2, 1, 3, 4)

def _coupling_block(B):
    d = B.shape[-1]//2
    a1 = (B[:d, d:] + B[d:, :d])/2.; a2 = (B[:d, d:] - B[d:, :d])/2.; a3 = (B[:d, :d] - B[d:, d:])/2.
    return max(np.max(np.abs(a1)), np.max(np.abs(a2)), np.max(np.abs(a3)))

def joint_flow(ks, omegas, cols, Jrho_s, omega0, Gamma0, rtol, atol, pole_threshold):
    S = L._spin_ops(ks); d = S.shape[-1]
    e = np.ascontiguousarray(ks.e, dtype=float); dE = e[:, None] - e[None, :]
    w = np.ascontiguousarray(omegas, dtype=float); cols = np.asarray(cols, dtype=int)
    n, m = len(w), len(cols)
    Dmax, Dmin = L._energy_range(ks, omega0, w, Gamma0)
    real = block_is_real(S)
    dt = float if real else complex
    Gb0 = L._block(L._bare_vertex(S, L.SIGN*Jrho_s, 0.)); Ab0 = L._block(L._bare_vertex(S, 1., 0.))
    T0 = Ab0.reshape(2, d, 2, d)[:, :, :, cols]
    nb = 4*d*d
    pre = llhole.precompute(w, e, omega0, Gamma0) if (real and USE_PRE and 16*n*d*d < 256e6) else None
    y0 = np.concatenate([Gb0.ravel(), Ab0.ravel(), np.broadcast_to(T0, (n, 2, d, 2, m)).ravel()])
    y0 = (y0.real if real else y0).astype(dt)
    smax = np.max(np.abs(S))/2.; ref = pole_threshold*smax
    def pack(y):
        return y[:2*nb].reshape(2, 2*d, 2*d)
    def rhs(s, y):
        D = np.exp(s)
        if not real: y = y.view(complex)
        X = pack(y)
        Kf = L._kernel(dE, D, omega0, Gamma0); Kb = np.tile(Kf, (2, 2))
        Gb, Ab = X[0], X[1]; GKb = Gb*Kb
        out = np.empty(len(y), dtype=dt)
        out[:nb] = L._bracket(Gb, GKb).ravel()
        out[nb:2*nb] = L._bracket(GKb, Ab).ravel()
        out[2*nb:] = tip_bracket_b(Gb, Ab, cols, w, e, D, -1., omega0, Gamma0, real, pre).ravel()
        out *= -D
        return out if real else out.view(float)
    yy0 = y0 if real else y0.view(float)
    sample_coupling = lambda y: _coupling_block(pack(y if real else y.view(complex))[0])
    solver = DOP853(rhs, np.log(Dmax), yy0, np.log(Dmin), rtol=rtol, atol=atol)
    g_old = sample_coupling(yy0); max_coupling = g_old; y, D_stop = solver.y, None
    while solver.status == "running":
        message = solver.step()
        if solver.status == "failed": raise RuntimeError(str(message))
        y = solver.y; g_new = sample_coupling(y)
        if g_old <= ref and g_new >= ref:
            dense = solver.dense_output()
            s_stop = brentq(lambda s: ref - sample_coupling(dense(s)), solver.t_old, solver.t, xtol=4*L._EPS, rtol=4*L._EPS)
            y = dense(s_stop); D_stop = float(np.exp(s_stop)); max_coupling = max(max_coupling, sample_coupling(y)); break
        max_coupling = max(max_coupling, g_new); g_old = g_new
    yc = y if real else y.view(complex)
    Tb = yc[2*nb:].reshape(n, 2, d, 2, m).astype(complex)
    if D_stop is not None and D_stop > Dmin:
        X = pack(yc)
        Tb = Tb + tip_bracket_b(X[0], X[1], cols, w, e, D_stop, Dmin, omega0, Gamma0, real)
    return L._cut_columns(Tb), {"D_stop": D_stop, "max_coupling": float(max_coupling/smax)}

def install():
    L._joint_flow = joint_flow
```

`llhole.py`:

```python
import numpy as np
from numba import njit
@njit(cache=True, fastmath=True)
def hole_sum_pre(Q, Rr, Wt, D, out):
    """out[k,a,l] = sum_c Q[a,l,c]*(1/(D+Rr[a,k,c]) - 1/(D+Wt[a,k,c])), Rr = sqrt(x^2+G^2), Wt = omega0+|x|, x=w[k]-e[a]+e[c]"""
    d, nk, _ = Rr.shape; nl = Q.shape[1]
    row = np.empty(d)
    for a in range(d):
        for k in range(nk):
            for c in range(d):
                r = Rr[a, k, c]; wt = Wt[a, k, c]
                row[c] = (wt - r)/((D + r)*(D + wt))
            for l in range(nl):
                acc = 0.
                for c in range(d): acc += Q[a, l, c]*row[c]
                out[k, a, l] = acc
def precompute(w, e, omega0, Gamma0):
    x = (w[None, :, None] - e[:, None, None]) + e[None, None, :]
    ax = np.abs(x)
    return np.ascontiguousarray(np.sqrt(ax*ax + Gamma0*Gamma0)), np.ascontiguousarray(omega0 + ax)
```

`16_equiv_timing.py`:

```python
"""(1) one right-hand side, old vs prototype, on the same state; (2) steps/evals split; (3) timing old / real+block / real+block+holepre,
median of 3, plus the rtol=1e-12 reference at the same size. Run on d22, one core."""
import numpy as np, warnings, time, sys
from common import chain
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
rec = {}
Orig = L.DOP853
class Spy(Orig):
    def __init__(self, fun, t0, y0, *a, **k):
        super().__init__(fun, t0, y0, *a, **k); rec["fun"] = fun; rec["y0"] = y0.copy(); rec["t0"] = t0; rec["n"] = 0; rec["self"] = self
    def step(self):
        r = super().step(); rec["n"] += 1; rec["nfev"] = self.nfev; return r
L.DOP853 = Spy
# the solver classes are referenced as L.DOP853 inside both flows (llreal imported DOP853 itself): patch it there too
llreal.DOP853 = Spy
n = int(sys.argv[1]); ks = chain(n); d = ks.dim
eVs = np.linspace(-4e-3, 4e-3, 61) + 1.3e-5
L.leading_log_dIdV(chain(3), np.array([1e-3, 2e-3]), -0.05, U=0.25); llreal.install(); llreal.USE_PRE = False; L.leading_log_dIdV(chain(3), np.array([1e-3, 2e-3]), -0.05, U=0.25)
# (1) rhs equivalence at a mid-flow state of the OLD solver
L._joint_flow = old
L.leading_log_dIdV(ks, eVs, -0.05, U=0.25)
sv = rec["self"]
# rebuild the old solver state after 25 steps
L.DOP853 = Orig; llreal.DOP853 = Orig
# reuse rec["fun"] (old rhs closure) with a state from a short re-run
s_ = rec["t0"]; fun_old = rec["fun"]; y = rec["y0"]
solver = Orig(fun_old, s_, y, np.log(1e-30), rtol=1e-10, atol=1e-14)
for _ in range(25): solver.step()
ymid, smid = solver.y, solver.t
nx = 8*d*d; yc = ymid.view(complex)
w = np.unique(np.concatenate([eVs, -eVs])); nb = 4*d*d; nw = len(w)
G = L._block(yc[:nx].reshape(2, 4, d, d)[0]); A = L._block(yc[:nx].reshape(2, 4, d, d)[1]); T = yc[nx:]
ynew = np.concatenate([G.ravel().real, A.ravel().real, T.real])
print("imag parts of the state in block form (should be ~0): G %.1e A %.1e T %.1e" % (np.max(np.abs(G.imag)), np.max(np.abs(A.imag)), np.max(np.abs(T.imag))))
out_old = fun_old(smid, ymid).view(complex)
Gr = L._block(out_old[:nx//2].reshape(4, d, d)); Ar = L._block(out_old[nx//2:nx].reshape(4, d, d))
ref_vec = np.concatenate([Gr.ravel(), Ar.ravel(), out_old[nx:]])
# new rhs: call through a fresh llreal flow's closure
rec2 = {}
class Spy2(Orig):
    def __init__(self, fun, t0, y0, *a, **k):
        super().__init__(fun, t0, y0, *a, **k); rec2["fun"] = fun
llreal.DOP853 = Spy2
L._joint_flow = llreal.joint_flow
for pre in (False, True):
    llreal.USE_PRE = pre
    L.leading_log_dIdV(ks, eVs, -0.05, U=0.25)  # builds closure (runs the whole flow, fine)
    out_new = rec2["fun"](smid, ynew)
    print("rhs, USE_PRE=%s: max |new-old| / max|old| = %.2e (imag part of old rhs in block form: %.1e)" % (pre, np.max(np.abs(out_new - ref_vec.real))/np.max(np.abs(ref_vec)), np.max(np.abs(ref_vec.imag))))
llreal.DOP853 = Spy; L.DOP853 = Spy
# (2)+(3)
def run(f, **kw):
    L._joint_flow = f
    t = time.perf_counter(); o = L.leading_log_dIdV(ks, eVs, -0.05, U=0.25, **kw); return o, time.perf_counter()-t, rec["n"], rec["nfev"]
res = {}
for label, f, pre in (("old", old, False), ("real+block", llreal.joint_flow, False), ("real+block+holepre", llreal.joint_flow, True)):
    llreal.USE_PRE = pre; ts = []
    for r in range(3):
        o, t, steps, nfev = run(f); ts.append(t)
    res[label] = (np.median(ts), o)
    print("n=%d d=%d %-20s median %.2f s (runs %s)  steps %d  rhs evals %d" % (n, d, label, np.median(ts), ", ".join("%.2f" % x for x in ts), steps, nfev), flush=True)
ref12, t12, s12, f12 = run(old, rtol=1e-12, atol=1e-16)
print("reference old rtol=1e-12: %.2f s, steps %d, evals %d" % (t12, s12, f12))
for k, (t, o) in res.items():
    print("  %-20s speedup %.2fx  |x-old|max rel %.2e   |x-ref1e-12|max rel %.2e" % (k, res["old"][0]/t, np.max(np.abs(o-res["old"][1])/np.abs(res["old"][1])), np.max(np.abs(o-ref12)/np.abs(ref12))))
```

`10_timing.py`:

```python
"""old vs real/block-form flow (llreal), and old vs rtol=1e-8, same inputs, median of 3, one core (d22)."""
import numpy as np, warnings, time, sys
from common import chain
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
def run(f, ks, eVs, **kw):
    L._joint_flow = f
    t = time.perf_counter(); out = L.leading_log_dIdV(ks, eVs, -0.05, U=0.25, **kw); return out, time.perf_counter()-t
L.leading_log_dIdV(chain(3), np.array([1e-3, 2e-3]), -0.05, U=0.25); llreal.install(); L.leading_log_dIdV(chain(3), np.array([1e-3, 2e-3]), -0.05, U=0.25)
for n in [int(a) for a in sys.argv[1:]] or (6, 7):
    ks = chain(n); eVs = np.linspace(-4e-3, 4e-3, 61) + 1.3e-5
    res = {}
    for label, f, kw in (("old", old, {}), ("real+block", llreal.joint_flow, {}), ("old rtol1e-8", old, dict(rtol=1e-8, atol=1e-12)), ("real+block rtol1e-8", llreal.joint_flow, dict(rtol=1e-8, atol=1e-12))):
        ts = []
        for r in range(3):
            out, t = run(f, ks, eVs, **kw); ts.append(t)
        res[label] = (np.median(ts), out)
        print("n=%d d=%d  %-22s median %.2f s (runs %s)" % (n, ks.dim, label, np.median(ts), ", ".join("%.2f" % t for t in ts)), flush=True)
    ref = res["old"][1]
    ref12, _ = run(old, ks, eVs, rtol=1e-12, atol=1e-16) if n <= 6 else (None, 0)
    for k, (t, o) in res.items():
        print("   %-22s speedup %.2fx  max rel diff vs old %.2e" % (k, res["old"][0]/t, np.max(np.abs(o-ref)/np.abs(ref))) + ("   vs rtol1e-12: %.2e" % np.max(np.abs(o-ref12)/np.abs(ref12)) if ref12 is not None else ""))
```

`08_realness.py`:

```python
"""Is the eigenbasis of a real-H chain real in KondoSpectrum (Sx,Sz real, Sy imaginary)? Heisenberg B=0 (degenerate) and tilted field."""
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
for label, field in (("tilted Bx,Bz", 1), ("B=0 Heisenberg (degenerate)", 0), ("By field (complex H)", 2)):
    n = 5
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = sum(1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 0.8*sc.Sz[i]*sc.Sz[i+1]) for i in range(n-1))
    if field == 1: h = h + sum(1e-4*(1+0.1*i)*sc.Sz[i] + 3e-5*sc.Sx[i] for i in range(n))
    if field == 2: h = h + sum(1e-4*sc.Sy[i] + 1e-4*sc.Sz[i] for i in range(n))
    sc.set_hamiltonian(h)
    ks = KondoSpectrum(sc, 0, 0.)
    sc_ = max(np.max(np.abs(x)) for x in (ks.Sx, ks.Sy, ks.Sz))
    print("%-30s d=%d  max|Im Sx|=%.1e  max|Im Sz|=%.1e  max|Re Sy|=%.1e  (scale %.2f) vs dtype %s" % (label, ks.dim, np.max(np.abs(ks.Sx.imag)), np.max(np.abs(ks.Sz.imag)), np.max(np.abs(ks.Sy.real)), sc_, ks.vs.dtype))
```

`15_real_stop.py`:

```python
"""prototype (real block form + precomputed hole sum) vs current code on STOPPED flows in the real path, and S=1 anisotropic."""
import numpy as np, warnings
from common import spin, ks_from
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
S = spin(.5); I2 = np.eye(2)
cases = {"S=1/2 Bz=0.3meV g=-0.2 (stops)": (ks_from(0.3e-3*S[2], S), -0.2, np.array([-1e-2, -1e-3, 1e-5, 3e-4, 1e-2])),
         "dimer Bz g=-0.2 (stops)": (ks_from(1e-3*sum(np.kron(x, x) for x in S) + 1e-4*np.kron(S[2], I2)+3e-5*np.kron(I2, S[0]), [np.kron(x, I2) for x in S]), -0.2, np.array([-1e-2, -3e-3, 1e-5, 2e-3, 1.2e-2]))}
S1 = spin(1.)
cases["S=1 D=-1meV B=(1,0,2) g=-0.15"] = (ks_from(-1e-3*S1[2]@S1[2] + 2*5.788e-5*(S1[0]+2*S1[2]), S1), -0.15, np.linspace(-4e-3, 4e-3, 9) + 1.3e-5)
for name, (ks, g, eVs) in cases.items():
    L._joint_flow = old; a, ia = L.leading_log_dIdV(ks, eVs, g, U=.25, return_info=True)
    llreal.install(); b, ib = L.leading_log_dIdV(ks, eVs, g, U=.25, return_info=True)
    ok = ~np.isnan(a)
    print("%-34s real=%s  nan old/new %d/%d same mask %s  max rel diff %.1e  D_stop %s vs %s" % (name, llreal.block_is_real(L._spin_ops(ks)), np.isnan(a).sum(), np.isnan(b).sum(), np.array_equal(np.isnan(a), np.isnan(b)), np.max(np.abs(a[ok]-b[ok])/np.abs(a[ok])) if ok.any() else 0, ia["D_stop"], ib["D_stop"]))
```

`17_degenerate_real_and_probes.py`:

```python
"""(a) prototype vs current with nocc>1 (B=0 trimer, real H, doublet GS; random-rotated doublet) and tip_vertex(cols=None);
(b) probes: info['coupling'] under a rotation inside the degenerate ground doublet; coupling when the occupied states carry no spin;
(c) keywords at the public call."""
import numpy as np, warnings
from types import SimpleNamespace
from common import spin, ks_from
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
S = spin(.5); I2 = np.eye(2)
def at(a, i, n=3):
    o = np.eye(1)
    for j in range(n): o = np.kron(o, S[a] if j == i else I2)
    return o
Htri = 1e-3*sum(at(a, i)@at(a, i+1) for i in range(2) for a in range(3))
ks = ks_from(Htri, [at(a, 1) for a in range(3)])
eVs = np.linspace(-4e-3, 4e-3, 17) + 1.3e-5
print("trimer nocc =", int(np.sum(ks.p > 0)), " real block:", llreal.block_is_real(L._spin_ops(ks)))
L._joint_flow = old; a, ia = L.leading_log_dIdV(ks, eVs, -0.05, U=0.25, return_info=True)
llreal.USE_PRE = True; L._joint_flow = llreal.joint_flow; b, ib = L.leading_log_dIdV(ks, eVs, -0.05, U=0.25, return_info=True)
print("(a) nocc=2: max rel diff dIdV %.1e, coupling %.1e" % (np.max(np.abs(a-b)/np.abs(a)), np.max(np.abs(ia["coupling"]-ib["coupling"]))))
# tip_vertex with cols=None
w = np.array([-3e-3, 1e-4, 2e-3]); Us = np.array([.1, .2, -.3])
L._joint_flow = old; A0, _ = L.tip_vertex(ks, w, Us, -0.05)
L._joint_flow = llreal.joint_flow; A1, _ = L.tip_vertex(ks, w, Us, -0.05)
print("(a) tip_vertex(cols=None) shapes", A0.shape, A1.shape, " max abs diff / max abs: %.1e" % (np.max(np.abs(A0-A1))/np.max(np.abs(A0))))
# (b1) rotation inside the doublet: ks has p on the 2 lowest states; rotate those two eigenstates by a random unitary
L._joint_flow = old
rng = np.random.default_rng(5)
occ = np.nonzero(ks.p > 0)[0]
def rotate(ks, th, ph):
    U = np.eye(ks.dim, dtype=complex)
    c, s = np.cos(th), np.sin(th)
    U[np.ix_(occ, occ)] = np.array([[c, -s*np.exp(-1j*ph)], [s*np.exp(1j*ph), c]])
    return SimpleNamespace(e=ks.e, p=ks.p, dim=ks.dim, T=0., kB=ks.kB, Sx=U.conj().T@ks.Sx@U, Sy=U.conj().T@ks.Sy@U, Sz=U.conj().T@ks.Sz@U)
for g in (-0.05, -0.1):
    base = None
    for th, ph in ((0, 0), (0.4, 0.3), (np.pi/4, 1.1), (1.2, 2.5)):
        k2 = rotate(ks, th, ph)
        o, info = L.leading_log_dIdV(k2, eVs, g, U=0.25, return_info=True)
        if base is None: base = (o, info["coupling"])
        print("(b1) g=%+.2f rotation(th=%.2f,ph=%.2f): max rel diff dIdV %.1e   coupling max %.4f  max |coupling - coupling(unrotated)| %.1e   NaNs %d" % (g, th, ph, np.nanmax(np.abs(o-base[0])/np.abs(base[0])), info["coupling"].max(), np.max(np.abs(info["coupling"]-base[1])), np.isnan(o).sum()))
# (b2) occupied states with no spin on the tip site: tip on an empty orbital-like site: occupied column of S is zero
d = 4
Sz = np.zeros((d, d), complex); Sx = np.zeros((d, d), complex); Sy = np.zeros((d, d), complex)
Sx[1:3, 1:3] = S[0]; Sy[1:3, 1:3] = S[1]; Sz[1:3, 1:3] = S[2]
ks0 = SimpleNamespace(e=np.array([0., 1e-3, 1.5e-3, 3e-3]), p=np.array([1., 0, 0, 0]), dim=4, T=0., kB=8.6e-5, Sx=Sx, Sy=Sy, Sz=Sz)
with warnings.catch_warnings(record=True) as w_:
    warnings.simplefilter("always")
    o, info = L.leading_log_dIdV(ks0, np.array([-2e-3, 5e-4, 2e-3]), -0.05, U=0.25, return_info=True)
print("(b2) zero-spin occupied state: dIdV", o, " coupling", info["coupling"], " D_stop", info["D_stop"], " warnings:", [str(x.message)[:60] for x in w_])
ks00 = SimpleNamespace(e=np.array([0., 1e-3]), p=np.array([.5, .5]), dim=2, T=0., kB=8.6e-5, Sx=0*S[0], Sy=0*S[1], Sz=0*S[2])
with warnings.catch_warnings(record=True) as w_:
    warnings.simplefilter("always")
    o, info = L.leading_log_dIdV(ks00, np.array([-2e-3, 5e-4, 2e-3]), -0.05, U=0.25, return_info=True)
print("(b2') S=0 entirely: dIdV", o, " coupling", info["coupling"], " D_stop", info["D_stop"], " warnings:", [str(x.message)[:90] for x in w_])
```

`19_public_profile.py`:

```python
"""Public call get_kondo_spectrum(order='LL', T=0) on 7 spins 1/2 (128 states), 61 biases: cProfile of the current code, then old vs prototype medians of 3 through the public call."""
import numpy as np, warnings, time, cProfile, pstats, io
from dmrgpy import spinchain
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
n = 7
sc = spinchain.Spin_Chain(["1/2"]*n)
h = sum(1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 0.8*sc.Sz[i]*sc.Sz[i+1]) for i in range(n-1))
h = h + sum(1e-4*(1+0.1*i)*sc.Sz[i] + 3e-5*sc.Sx[i] for i in range(n))
sc.set_hamiltonian(h)
eVs = np.linspace(-4e-3, 4e-3, 61) + 1.3e-5
call = lambda: sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, U=0.25, T=0., order="LL")
call()  # compile + cache the diagonalization
pr = cProfile.Profile(); t0 = time.time(); pr.enable(); a = call()[1]; pr.disable(); tot = time.time()-t0
print("public call, current code: %.2f s" % tot)
s = io.StringIO(); pstats.Stats(pr, stream=s).sort_stats("cumtime").print_stats(14); print(s.getvalue()[:3600])
old = L._joint_flow
def timed(f):
    L._joint_flow = f; ts = []
    for _ in range(3):
        t = time.perf_counter(); o = call()[1]; ts.append(time.perf_counter()-t)
    return np.median(ts), o, ts
llreal.USE_PRE = True
llreal.install(); L.leading_log_dIdV(chain3 := None, [], 0) if False else None
t_old, o_old, r_old = timed(old); t_new, o_new, r_new = timed(llreal.joint_flow)
print("public call median of 3: old %.2f s %s, prototype %.2f s %s, speedup %.2fx, max rel diff %.1e" % (t_old, np.round(r_old, 2), t_new, np.round(r_new, 2), t_old/t_new, np.max(np.abs(o_old-o_new)/np.abs(o_old))))
```

Observed:

`16_equiv_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`10_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`08_realness.out`:

```
tilted Bx,Bz                   d=32  max|Im Sx|=0.0e+00  max|Im Sz|=0.0e+00  max|Re Sy|=0.0e+00  (scale 0.49) vs dtype complex128
B=0 Heisenberg (degenerate)    d=32  max|Im Sx|=0.0e+00  max|Im Sz|=0.0e+00  max|Re Sy|=0.0e+00  (scale 0.50) vs dtype complex128
By field (complex H)           d=32  max|Im Sx|=3.3e-01  max|Im Sz|=2.1e-01  max|Re Sy|=4.5e-01  (scale 0.45) vs dtype complex128
```

`15_real_stop.out`:

```
S=1/2 Bz=0.3meV g=-0.2 (stops)     real=True  nan old/new 3/3 same mask True  max rel diff 1.2e-10  D_stop 0.00021579239647816011 vs 0.00021579239625258066
dimer Bz g=-0.2 (stops)            real=True  nan old/new 2/2 same mask True  max rel diff 5.5e-11  D_stop 0.00011244635155873904 vs 0.00011244635150224161
S=1 D=-1meV B=(1,0,2) g=-0.15      real=True  nan old/new 0/0 same mask True  max rel diff 3.3e-11  D_stop None vs None
```

`17_degenerate_real_and_probes.out`:

```
trimer nocc = 2  real block: True
(a) nocc=2: max rel diff dIdV 7.9e-11, coupling 3.7e-12
(a) tip_vertex(cols=None) shapes (3, 4, 8, 8) (3, 4, 8, 8)  max abs diff / max abs: 1.9e-11
(b1) g=-0.05 rotation(th=0.00,ph=0.00): max rel diff dIdV 0.0e+00   coupling max 0.0682  max |coupling - coupling(unrotated)| 0.0e+00   NaNs 0
(b1) g=-0.05 rotation(th=0.40,ph=0.30): max rel diff dIdV 2.5e-11   coupling max 0.0682  max |coupling - coupling(unrotated)| 1.5e-12   NaNs 0
(b1) g=-0.05 rotation(th=0.79,ph=1.10): max rel diff dIdV 5.1e-11   coupling max 0.0682  max |coupling - coupling(unrotated)| 2.2e-12   NaNs 0
(b1) g=-0.05 rotation(th=1.20,ph=2.50): max rel diff dIdV 2.1e-11   coupling max 0.0682  max |coupling - coupling(unrotated)| 1.1e-12   NaNs 0
(b1) g=-0.10 rotation(th=0.00,ph=0.00): max rel diff dIdV 0.0e+00   coupling max 0.2126  max |coupling - coupling(unrotated)| 0.0e+00   NaNs 0
(b1) g=-0.10 rotation(th=0.40,ph=0.30): max rel diff dIdV 6.7e-11   coupling max 0.2126  max |coupling - coupling(unrotated)| 8.1e-12   NaNs 0
(b1) g=-0.10 rotation(th=0.79,ph=1.10): max rel diff dIdV 1.2e-10   coupling max 0.2126  max |coupling - coupling(unrotated)| 1.4e-11   NaNs 0
(b1) g=-0.10 rotation(th=1.20,ph=2.50): max rel diff dIdV 7.4e-11   coupling max 0.2126  max |coupling - coupling(unrotated)| 8.8e-12   NaNs 0
(b2) zero-spin occupied state: dIdV [1.57079633 1.57079633 1.57079633]  coupling [nan nan nan]  D_stop None  warnings: ['invalid value encountered in divide']
(b2') S=0 entirely: dIdV [nan nan nan]  coupling [nan nan nan]  D_stop 230000000.00000015  warnings: ['invalid value encountered in scalar divide', 'invalid value encountered in divide', 'leading_log_dIdV: the sample coupling reached order one at D = 2.3e+08 (the Kondo scale is']
```

`19_public_profile.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`09_real_check.out`:

```
4 max rel diff 1.4479535289506981e-11 coupling diff 5.81756864903582e-13 0.0646683621787231 0.06466836217812141
5 max rel diff 3.3206355763397695e-11 coupling diff 1.3085088568232095e-12 0.06428689750398003 0.06428689750271242
dimer ByBz (complex path) g=-0.05 real?False 6.183811507308057e-12 None None 0 0
dimer ByBz (complex path) g=-0.2 real?False 1.5959981521439282e-12 0.00010198288116540182 0.0001019828811680509 4 4
```

**Reviewer (CONFIRMED)**: I tried to refute it along five lines. (1) Documented or rejected: kondo.md lists the GPU, complex64, numba threads in _hole_sum and the one-solve flow, and does not mention real arithmetic. The "Do NOT" items (no symmetrizing G(0), no onshell_flow plus per-bias flow) are untouched, because the prototype keeps the same equations and does not symmetrize. (2) Probe: I used my own public-call harness (02_timing.py, d22, one core, pinned threads, median of 3), and it reproduces the hunter's sizes within 4%. (3) Anchor: I compared both codes against an rtol=1e-12, atol=1e-16 reference of the old code on 8 models, each with 41 biases (03_inputs.py), and the new code is no further from the reference than the old one. 5x1/2: old-ref 1.8e-11, new-ref 2.4e-11. 4x1/2 B=0 degenerate: 3.1e-11 and 2.7e-11. Hubbard: 3.0e-11 and 3.1e-11. The 5x1/2 case with a nearly real By=1e-9 is the one exception, at new-ref 4.6e-11, still well under rtol=1e-10. The NaN masks are identical in all cases. (4) Untested paths: By=3e-5 (complex eigenbasis) falls back to the complex route and gives new-old 3.3e-11. By=1e-9 is complex above the 1e-12 tolerance and also falls back. The degenerate B=0 chain, S=1 with Dz, and the interleaved and native spinful Hubbard chains all agree to between 3e-12 and 6e-11. The strong-coupling case (g=-0.2) is NaN everywhere in both codes with the same mask. (5) The pinned test file tests/test_kondo_leading_log.py with the prototype installed: 21 of 22 pass. The one failure, test_gamma0_must_be_positive, is a prototype omission, because llreal.joint_flow dropped the `if not Gamma0 > 0.: raise ValueError` of the original. The port must keep it. Struck: nothing in the claim. I only correct the "roundoff or approximation" question: the 1e-11 movement is neither. It is the integrator's step sequence changing because the error norm sees a different state layout, and it stays inside the integration tolerance. Against the tight reference the new code is no worse. The hunter's side claim that block-reality holds for every real H also holds, because a real eigenbasis gives real Sx and Sz and imaginary Sy.

Reviewer's own evidence:

```
02_timing.out: n=6 old 2.79 [2.79 2.82 2.79] | real+block 1.41 [1.41 1.41 1.4] (1.98x, rel diff 1.2e-12); n=7 old 13.57 | real+block 5.52 (2.46x, rel diff 2.7e-11). 03_inputs.out: '5x1/2 By=3e-5 (complex) NaNs old/new 0/0 same mask True | new-old 3.3e-11 | old-ref 1.3e-11 new-ref 2.2e-11'; '4x1/2 B=0 degenerate ... new-old 9.5e-12 | old-ref 3.1e-11 new-ref 2.7e-11'; 'hubbard native 3 sites Nf=3 ... new-old 3.0e-12'. 04_pinned_tests.out: '1 failed, 21 passed' (the failure is the missing Gamma0 check in the prototype).
```

Reviewer's measurement: Own medians of 3 on d22 through the public call, T=0: 6 spins 2.79 s to 1.41 s (1.98x, 49% saved), 7 spins 13.57 s to 5.52 s (2.46x, 59% saved). Max relative difference of dI/dV old to new: 1.2e-12 at 6 spins, 2.7e-11 at 7 spins, up to 6.3e-11 on the other inputs. The new code is as close to the rtol=1e-12 reference as the old one.

Reviewer's scripts and outputs:

`02_timing.py`:

```python
"""Own timing: public get_kondo_spectrum(order='LL') old vs real+block (+hole precompute), median of 3, 7 spins 1/2 and 6 spins, real H."""
import numpy as np, warnings, time, sys
from dmrgpy import spinchain
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
def build(n, by=0.):
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = sum(1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 0.8*sc.Sz[i]*sc.Sz[i+1]) for i in range(n-1))
    h = h + sum(1e-4*(1+0.1*i)*sc.Sz[i] + 3e-5*sc.Sx[i] + by*sc.Sy[i] for i in range(n))
    sc.set_hamiltonian(h); return sc
eVs = np.linspace(-4e-3, 4e-3, 61) + 1.3e-5
def med(sc, f):
    L._joint_flow = f; ts = []
    for _ in range(3):
        t = time.perf_counter(); o = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, U=0.25, T=0., order="LL")[1]; ts.append(time.perf_counter()-t)
    return np.median(ts), o, ts
for n in (6, 7):
    sc = build(n); sc.get_kondo_spectrum(eVs[:3], site=0, Jrho_s=-0.05, U=0.25, T=0., order="LL")
    llreal.install(); L._joint_flow = old
    to, oo, ro = med(sc, old)
    llreal.USE_PRE = False; tr, orr, rr = med(sc, llreal.joint_flow)
    llreal.USE_PRE = True; tp, op, rp = med(sc, llreal.joint_flow)
    print("n=%d old %.2f %s | real+block %.2f %s (%.2fx, rel diff %.1e) | +pre %.2f %s (%.2fx, rel diff %.1e)" % (n, to, np.round(ro,2), tr, np.round(rr,2), to/tr, np.max(np.abs(orr-oo)/np.abs(oo)), tp, np.round(rp,2), to/tp, np.max(np.abs(op-oo)/np.abs(oo))), flush=True)
```

`03_inputs.py`:

```python
"""Old vs real-block prototype on inputs the hunter did not use, through the public call."""
import numpy as np, warnings, sys
sys.path.insert(0, ".")
import llreal
from dmrgpy import spinchain, fermionchain
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow; llreal.install(); new = L._joint_flow
def spinc(sp, n, by=0., bz=1e-4, bx=3e-5, aniso=0.8, Dz=0.):
    sc = spinchain.Spin_Chain([sp]*n)
    h = sum(1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + aniso*sc.Sz[i]*sc.Sz[i+1]) for i in range(n-1))
    h = h + sum(bz*sc.Sz[i] + bx*sc.Sx[i] + by*sc.Sy[i] + Dz*sc.Sz[i]*sc.Sz[i] for i in range(n))
    sc.set_hamiltonian(h); return sc
def hub(n, nf):
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = sum(-2e-3*(fc.Cdag[2*i+s]*fc.C[2*(i+1)+s] + fc.Cdag[2*(i+1)+s]*fc.C[2*i+s]) for i in range(n-1) for s in (0,1))
    h = h + sum(1e-3*fc.Nup[i]*fc.Ndn[i] for i in range(n)) + 1e-4*fc.Sz[0]
    fc.set_hamiltonian(h); fc.set_conserved_sector(Nf=nf); return fc
def hubn(n, nf):
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = sum(-2e-3*(fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagup[i+1]*fc.Cup[i] + fc.Cdagdn[i]*fc.Cdn[i+1] + fc.Cdagdn[i+1]*fc.Cdn[i]) for i in range(n-1))
    h = h + sum(1e-3*fc.Nup[i]*fc.Ndn[i] for i in range(n)) + 1e-4*fc.Sz[0]
    fc.set_hamiltonian(h); fc.set_conserved_sector(Nf=nf); return fc
eVs = np.linspace(-4e-3, 4e-3, 41) + 1.3e-5
cases = [("5x1/2 By=0", lambda: spinc("1/2",5), -0.05),
         ("5x1/2 By=3e-5 (complex)", lambda: spinc("1/2",5,by=3e-5), -0.05),
         ("5x1/2 By=1e-9 (nearly real)", lambda: spinc("1/2",5,by=1e-9), -0.05),
         ("4x1/2 B=0 degenerate", lambda: spinc("1/2",4,bz=0,bx=0), -0.05),
         ("3x1 Dz", lambda: spinc("1",3,Dz=2e-4), -0.05),
         ("5x1/2 strong g=-0.2 (stops?)", lambda: spinc("1/2",5), -0.2),
         ("hubbard 3 sites Nf=3", lambda: hub(3,3), -0.05),
         ("hubbard native 3 sites Nf=3", lambda: hubn(3,3), -0.05)]
for name, mk, g in cases:
    try:
        sc = mk(); out = {}
        tight = lambda *a: old(*a[:6], 1e-12, 1e-16, *a[8:])
        for lab, f, kw in (("old", old, {}), ("new", new, {}), ("ref", tight, {})):
            L._joint_flow = f
            out[lab] = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=g, U=0.25, T=0., order="LL", **kw)[1]
        ks = sc._kondo_ks if hasattr(sc, "_kondo_ks") else None
        mask = ~np.isnan(out["old"])
        same = (np.isnan(out["old"]) == np.isnan(out["new"])).all()
        rd = lambda a, b: np.nanmax(np.abs(a-b)/np.abs(b)) if mask.any() else np.nan
        print("%-32s NaNs old/new %d/%d same mask %s | new-old %.1e | old-ref %.1e new-ref %.1e" % (name, np.isnan(out["old"]).sum(), np.isnan(out["new"]).sum(), same, rd(out["new"], out["old"]), rd(out["old"], out["ref"]), rd(out["new"], out["ref"])), flush=True)
    except Exception as ex:
        print(name, "ERR", type(ex).__name__, str(ex)[:150], flush=True)
    L._joint_flow = old
```

`04_pinned_tests.py`:

```python
"""Run tests/test_kondo_leading_log.py with the real-block prototype installed (and the hole precompute)."""
import sys, pytest
sys.path.insert(0, ".")
sys.path.insert(0, "<repo>/src")
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
llreal.install(); llreal.USE_PRE = True
assert L._joint_flow is llreal.joint_flow
sys.exit(pytest.main(["-q", "-p", "no:cacheprovider", "<repo>/tests/test_kondo_leading_log.py"]))
```

`02_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`03_inputs.out`:

```
5x1/2 By=0                       NaNs old/new 0/0 same mask True | new-old 4.2e-11 | old-ref 1.8e-11 new-ref 2.4e-11
5x1/2 By=3e-5 (complex)          NaNs old/new 0/0 same mask True | new-old 3.3e-11 | old-ref 1.3e-11 new-ref 2.2e-11
5x1/2 By=1e-9 (nearly real)      NaNs old/new 0/0 same mask True | new-old 6.3e-11 | old-ref 1.8e-11 new-ref 4.6e-11
4x1/2 B=0 degenerate             NaNs old/new 0/0 same mask True | new-old 9.5e-12 | old-ref 3.1e-11 new-ref 2.7e-11
3x1 Dz                           NaNs old/new 0/0 same mask True | new-old 4.2e-12 | old-ref 1.6e-11 new-ref 1.2e-11
5x1/2 strong g=-0.2 (stops?)     NaNs old/new 41/41 same mask True | new-old nan | old-ref nan new-ref nan
hubbard 3 sites Nf=3             NaNs old/new 0/0 same mask True | new-old 3.0e-12 | old-ref 3.0e-11 new-ref 3.1e-11
hubbard native 3 sites Nf=3      NaNs old/new 0/0 same mask True | new-old 3.0e-12 | old-ref 3.0e-11 new-ref 3.1e-11
```

`04_pinned_tests.out`:

```
.....F................                                                   [100%]
=================================== FAILURES ===================================
_________________________ test_gamma0_must_be_positive _________________________

    def test_gamma0_must_be_positive():
        S = _spin(0.5)
        ks = SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=0.,
                             kB=8.617333262e-5, Sx=S[0], Sy=S[1], Sz=S[2])
>       with pytest.raises(ValueError):
             ^^^^^^^^^^^^^^^^^^^^^^^^^
E       Failed: DID NOT RAISE <class 'ValueError'>

<repo>/tests/test_kondo_leading_log.py:176: Failed
=============================== warnings summary ===============================
tests/test_kondo_leading_log.py::test_gamma0_must_be_positive
  <scratch>/review/leading-log/llreal.py:75: RuntimeWarning: divide by zero encountered in log
    solver = DOP853(rhs, np.log(Dmax), yy0, np.log(Dmin), rtol=rtol, atol=atol)

tests/test_kondo_leading_log.py::test_gamma0_must_be_positive
  <repo>/src/dmrgpy/kondospectrumtk/leadinglog.py:87: RuntimeWarning: divide by zero encountered in divide
    return 1./(D + np.sqrt(ax*ax + Gamma0*Gamma0)) - 1./(D + omega0 + ax)

tests/test_kondo_leading_log.py::test_gamma0_must_be_positive
  <scratch>/review/leading-log/llreal.py:66: RuntimeWarning: invalid value encountered in multiply
    Gb, Ab = X[0], X[1]; GKb = Gb*Kb

tests/test_kondo_leading_log.py::test_gamma0_must_be_positive
  <repo>/src/dmrgpy/kondospectrumtk/leadinglog.py:138: RuntimeWarning: invalid value encountered in matmul
    return np.matmul(Xb, Yb) - _swap(np.matmul(_swap(Yb), _swap(Xb)))

tests/test_kondo_leading_log.py::test_get_kondo_spectrum_order_LL
  <repo>/src/dmrgpy/manybodychain.py:24: PendingDeprecationWarning: the matrix subclass is not the recommended way to represent matrices or deal with linear algebra (see https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html). Please adjust your code to use regular ndarray.
    one = np.matrix(np.identity(3))

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
=========================== short test summary info ============================
FAILED ../../../../../../..<repo>/tests/test_kondo_leading_log.py::test_gamma0_must_be_positive
1 failed, 21 passed, 5 warnings in 5.71s
```

**Suggested fix (hunter)**: In _joint_flow, test once whether the block form of the bare vertex is real (_block of the stacked [0,Sx,Sy,Sz] has no imaginary part above about 1e-13 of its largest entry). If it is, keep y as the real concatenation of the two 2d x 2d block matrices of G and A_on and the (n,2,d,2,m) block tip array, write the right-hand side directly in block form (the rate of G is _bracket(Gb, Gb*tile(Kf,(2,2))), no _unblock/_block), call np.matmul on real arrays in _tip_bracket (no P.imag, Q.imag, and _hole_sum with the 4m real columns), and convert to the (4,d,d) form only in _cut_columns at the end. Keep the present complex code as the fallback for a complex H (a By field, a complex eigenbasis). It changes returned numbers at the 1e-11 level (integrator step sequence); the equivalence of the algebra is 2e-15 per evaluation.

**Fix (reviewer)**: The hunter's fix is right: in _joint_flow test once that _block of the bare vertex is real (relative 1e-12), keep y real in block form, and keep the present complex code as the fallback. Three points for the port. Keep the Gamma0>0 ValueError, which the prototype dropped and which test_gamma0_must_be_positive pins. Take sample_coupling in block form, as the prototype's _coupling_block does. Do not touch the step control. Pin it in tests/test_kondo_leading_log.py with a real-path model against _TWO_PASS (rel=1e-9 stays valid), a By-field model that goes through the complex fallback and agrees with the real route's neighbour value, and a block-reality assertion. The ledger entry in docs/numerical_changes.md should say 1e-11, from the step sequence.

### 9. _hole_sum recomputes D-independent abs and sqrt at every step: 1.15x, and 1.37x on top of finding 8

`performance` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `leading-log` &middot; candidate `leading-log-2`

**Status**: FIXED in `fce3382` on 2026-11-12. `_hole_tables` builds r and omega0+|x| once per flow, and `_hole_sum_table` does one division per element on both routes; the on-the-fly `_hole_sum` stays for the closed-form tail and above `HOLE_TABLE_BYTES` = 256 MB, so the 9-spin call with 241 biases is unchanged, as the method's README now says. One core, no threads. Alone, on the complex route at 6 spins: 2.83 s to 2.41 s, 9e-16 apart. Pinned by `tests/test_audit_2026_11_leadinglog.py::test_the_table_sum_is_the_sum_on_the_fly` and `::test_the_table_is_off_above_the_cap_and_the_numbers_do_not_move`.

The claim as reviewed: _hole_sum recomputes abs and sqrt for every (bias, final state, intermediate state) at every right-hand-side evaluation although none of it depends on D. Precomputing r and omega0+|x| once per flow makes it 1.17x faster at 6 spins (2.82 s to 2.41 s) and 1.15x at 7 spins (13.52 s to 11.71 s) on the current complex flow, with a difference of 7e-16 to 9e-16, and 1.37x on top of leading-log-1 (5.52 s to 4.04 s at 7 spins). The 256 MB gate (16*n*d*d bytes) switches it off at 9 spins (and at 8 spins with 241 biases), exactly where a call costs minutes, so it is a gain for the cheap sizes only.

**Where**: src/dmrgpy/kondospectrumtk/leadinglog.py:159-187 (_hole_sum, called from :215 for every rhs evaluation), :346-356 (rhs)

WHAT THE CODE DOES. _hole_sum(Q,w,e,D,...) computes, for each (a,k), row[c] = (omega0+ax-r)/((D+r)(D+omega0+ax)) with ax=|w[k]-e[a]+e[c]| and r=sqrt(ax^2+Gamma0^2). Only D changes along the flow, ax and r are constants, but are recomputed 518 to 626 times (abs and sqrt per element, d*n*d = 2.0e6 elements at 7 spins with 122 distinct biases). After leading-log-1 this routine is about 57% of the call (profile 11: tip_bracket 3.14 s of 5.52 s, most of it the hole sum).

PROTOTYPE (llhole.py): precompute R[a,k,c]=r and Wt[a,k,c]=omega0+ax once per flow, a numba kernel then does one division per element. Micro (12_hole_micro.out, d22, median of 21, pinned): d=64, n=122: 1.36 ms to 0.55 ms (2.49x); d=128, n=122: 5.45 ms to 2.46 ms (2.21x); max relative difference 2.3e-16. Extra memory 2 n d^2 doubles: 8 MB at d=64 and 32 MB at d=128. SCOPE: the saving was measured on top of leading-log-1 (real path), end to end 7 spins 5.50 s to 4.02 s (1.37x), 6 spins 1.49 s to about 1.13 s; alone on the current complex code it is a smaller share (the hole sum is about 28% of the old rhs). It must be gated by size: at 8 spins with 241 biases the two arrays are 0.5 GB and at 9 spins with 241 biases 2.0 GB (the README's 1.7 GB peak run), so above roughly 256 MB the existing on-the-fly kernel has to stay (the prototype gates at 16*n*d*d < 256e6). This does not thread anything and does not use the GPU or complex64.

WHY NO TEST SEES IT. Pure speed, numbers equal to roundoff.

The hunter's claim, before review: _hole_sum recomputes abs, sqrt and two additions for every (bias, final state, intermediate state) at every right-hand-side evaluation although they do not depend on D; precomputing r=sqrt(x^2+Gamma0^2) and omega0+|x| once per flow makes the sum 2.2 to 2.5x faster and the whole 7-spin call a further 1.37x (5.50 s to 4.02 s) on top of leading-log-1, with a difference of 2e-16.

Decisive output, quoted by the hunter:

```
d=64 n=122: _hole_sum 1.36 ms, precomputed 0.55 ms (2.49x), max rel diff 1.9e-16, extra memory 8 MB
d=128 n=122: _hole_sum 5.45 ms, precomputed 2.46 ms (2.21x), max rel diff 2.3e-16, extra memory 32 MB
n=7 d=128 real+block+holepre   median 4.02 s (runs 4.05, 4.01, 4.02)  steps 38  rhs evals 518
  real+block+holepre   speedup 3.55x  |x-old|max rel 2.70e-11   |x-ref1e-12|max rel 9.77e-12
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`llhole.py`:

```python
import numpy as np
from numba import njit
@njit(cache=True, fastmath=True)
def hole_sum_pre(Q, Rr, Wt, D, out):
    """out[k,a,l] = sum_c Q[a,l,c]*(1/(D+Rr[a,k,c]) - 1/(D+Wt[a,k,c])), Rr = sqrt(x^2+G^2), Wt = omega0+|x|, x=w[k]-e[a]+e[c]"""
    d, nk, _ = Rr.shape; nl = Q.shape[1]
    row = np.empty(d)
    for a in range(d):
        for k in range(nk):
            for c in range(d):
                r = Rr[a, k, c]; wt = Wt[a, k, c]
                row[c] = (wt - r)/((D + r)*(D + wt))
            for l in range(nl):
                acc = 0.
                for c in range(d): acc += Q[a, l, c]*row[c]
                out[k, a, l] = acc
def precompute(w, e, omega0, Gamma0):
    x = (w[None, :, None] - e[:, None, None]) + e[None, None, :]
    ax = np.abs(x)
    return np.ascontiguousarray(np.sqrt(ax*ax + Gamma0*Gamma0)), np.ascontiguousarray(omega0 + ax)
```

`12_hole_micro.py`:

```python
import numpy as np, time
from common import chain
import llhole
from dmrgpy.kondospectrumtk import leadinglog as L
def med(f, r=21):
    ts = []
    for _ in range(r):
        t = time.perf_counter(); f(); ts.append(time.perf_counter()-t)
    return np.median(ts)
for nspin in (6, 7):
    ks = chain(nspin); d = ks.dim; e = np.ascontiguousarray(ks.e)
    w = np.unique(np.concatenate([np.linspace(-4e-3, 4e-3, 61) + 1.3e-5]*1 + [-np.linspace(-4e-3, 4e-3, 61) - 1.3e-5])); n = len(w)
    rng = np.random.default_rng(0); Q = rng.normal(size=(d, 4, d)); D = 3e-4; w0, G0 = 20e-3, 5e-6
    o1 = np.empty((n, d, 4)); o2 = np.empty((n, d, 4))
    L._hole_sum(Q, w, e, D, -1., w0, G0, o1)
    R, Wt = llhole.precompute(w, e, w0, G0); llhole.hole_sum_pre(Q, R, Wt, D, o2)
    t1 = med(lambda: L._hole_sum(Q, w, e, D, -1., w0, G0, o1)); t2 = med(lambda: llhole.hole_sum_pre(Q, R, Wt, D, o2))
    print("d=%d n=%d: _hole_sum %.2f ms, precomputed %.2f ms (%.2fx), max rel diff %.1e, extra memory %.0f MB" % (d, n, t1*1e3, t2*1e3, t1/t2, np.max(np.abs(o1-o2))/np.max(np.abs(o1)), (R.nbytes+Wt.nbytes)/1e6))
```

`16_equiv_timing.py`:

```python
"""(1) one right-hand side, old vs prototype, on the same state; (2) steps/evals split; (3) timing old / real+block / real+block+holepre,
median of 3, plus the rtol=1e-12 reference at the same size. Run on d22, one core."""
import numpy as np, warnings, time, sys
from common import chain
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
rec = {}
Orig = L.DOP853
class Spy(Orig):
    def __init__(self, fun, t0, y0, *a, **k):
        super().__init__(fun, t0, y0, *a, **k); rec["fun"] = fun; rec["y0"] = y0.copy(); rec["t0"] = t0; rec["n"] = 0; rec["self"] = self
    def step(self):
        r = super().step(); rec["n"] += 1; rec["nfev"] = self.nfev; return r
L.DOP853 = Spy
# the solver classes are referenced as L.DOP853 inside both flows (llreal imported DOP853 itself): patch it there too
llreal.DOP853 = Spy
n = int(sys.argv[1]); ks = chain(n); d = ks.dim
eVs = np.linspace(-4e-3, 4e-3, 61) + 1.3e-5
L.leading_log_dIdV(chain(3), np.array([1e-3, 2e-3]), -0.05, U=0.25); llreal.install(); llreal.USE_PRE = False; L.leading_log_dIdV(chain(3), np.array([1e-3, 2e-3]), -0.05, U=0.25)
# (1) rhs equivalence at a mid-flow state of the OLD solver
L._joint_flow = old
L.leading_log_dIdV(ks, eVs, -0.05, U=0.25)
sv = rec["self"]
# rebuild the old solver state after 25 steps
L.DOP853 = Orig; llreal.DOP853 = Orig
# reuse rec["fun"] (old rhs closure) with a state from a short re-run
s_ = rec["t0"]; fun_old = rec["fun"]; y = rec["y0"]
solver = Orig(fun_old, s_, y, np.log(1e-30), rtol=1e-10, atol=1e-14)
for _ in range(25): solver.step()
ymid, smid = solver.y, solver.t
nx = 8*d*d; yc = ymid.view(complex)
w = np.unique(np.concatenate([eVs, -eVs])); nb = 4*d*d; nw = len(w)
G = L._block(yc[:nx].reshape(2, 4, d, d)[0]); A = L._block(yc[:nx].reshape(2, 4, d, d)[1]); T = yc[nx:]
ynew = np.concatenate([G.ravel().real, A.ravel().real, T.real])
print("imag parts of the state in block form (should be ~0): G %.1e A %.1e T %.1e" % (np.max(np.abs(G.imag)), np.max(np.abs(A.imag)), np.max(np.abs(T.imag))))
out_old = fun_old(smid, ymid).view(complex)
Gr = L._block(out_old[:nx//2].reshape(4, d, d)); Ar = L._block(out_old[nx//2:nx].reshape(4, d, d))
ref_vec = np.concatenate([Gr.ravel(), Ar.ravel(), out_old[nx:]])
# new rhs: call through a fresh llreal flow's closure
rec2 = {}
class Spy2(Orig):
    def __init__(self, fun, t0, y0, *a, **k):
        super().__init__(fun, t0, y0, *a, **k); rec2["fun"] = fun
llreal.DOP853 = Spy2
L._joint_flow = llreal.joint_flow
for pre in (False, True):
    llreal.USE_PRE = pre
    L.leading_log_dIdV(ks, eVs, -0.05, U=0.25)  # builds closure (runs the whole flow, fine)
    out_new = rec2["fun"](smid, ynew)
    print("rhs, USE_PRE=%s: max |new-old| / max|old| = %.2e (imag part of old rhs in block form: %.1e)" % (pre, np.max(np.abs(out_new - ref_vec.real))/np.max(np.abs(ref_vec)), np.max(np.abs(ref_vec.imag))))
llreal.DOP853 = Spy; L.DOP853 = Spy
# (2)+(3)
def run(f, **kw):
    L._joint_flow = f
    t = time.perf_counter(); o = L.leading_log_dIdV(ks, eVs, -0.05, U=0.25, **kw); return o, time.perf_counter()-t, rec["n"], rec["nfev"]
res = {}
for label, f, pre in (("old", old, False), ("real+block", llreal.joint_flow, False), ("real+block+holepre", llreal.joint_flow, True)):
    llreal.USE_PRE = pre; ts = []
    for r in range(3):
        o, t, steps, nfev = run(f); ts.append(t)
    res[label] = (np.median(ts), o)
    print("n=%d d=%d %-20s median %.2f s (runs %s)  steps %d  rhs evals %d" % (n, d, label, np.median(ts), ", ".join("%.2f" % x for x in ts), steps, nfev), flush=True)
ref12, t12, s12, f12 = run(old, rtol=1e-12, atol=1e-16)
print("reference old rtol=1e-12: %.2f s, steps %d, evals %d" % (t12, s12, f12))
for k, (t, o) in res.items():
    print("  %-20s speedup %.2fx  |x-old|max rel %.2e   |x-ref1e-12|max rel %.2e" % (k, res["old"][0]/t, np.max(np.abs(o-res["old"][1])/np.abs(res["old"][1])), np.max(np.abs(o-ref12)/np.abs(ref12))))
```

`13_timing_pre.py`:

```python
"""old vs real/block-form flow (llreal), and old vs rtol=1e-8, same inputs, median of 3, one core (d22)."""
import numpy as np, warnings, time, sys
from common import chain
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
def run(f, ks, eVs, **kw):
    L._joint_flow = f
    t = time.perf_counter(); out = L.leading_log_dIdV(ks, eVs, -0.05, U=0.25, **kw); return out, time.perf_counter()-t
L.leading_log_dIdV(chain(3), np.array([1e-3, 2e-3]), -0.05, U=0.25); llreal.install(); L.leading_log_dIdV(chain(3), np.array([1e-3, 2e-3]), -0.05, U=0.25)
for n in [int(a) for a in sys.argv[1:]] or (6, 7):
    ks = chain(n); eVs = np.linspace(-4e-3, 4e-3, 61) + 1.3e-5
    res = {}
    for label, f, kw in (("old", old, {}), ("real+block+holepre", llreal.joint_flow, {}), ("old rtol1e-8", old, dict(rtol=1e-8, atol=1e-12)), ("real+block+holepre rtol1e-8", llreal.joint_flow, dict(rtol=1e-8, atol=1e-12))):
        ts = []
        for r in range(3):
            out, t = run(f, ks, eVs, **kw); ts.append(t)
        res[label] = (np.median(ts), out)
        print("n=%d d=%d  %-22s median %.2f s (runs %s)" % (n, ks.dim, label, np.median(ts), ", ".join("%.2f" % t for t in ts)), flush=True)
    ref = res["old"][1]
    ref12, _ = run(old, ks, eVs, rtol=1e-12, atol=1e-16) if n <= 6 else (None, 0)
    for k, (t, o) in res.items():
        print("   %-22s speedup %.2fx  max rel diff vs old %.2e" % (k, res["old"][0]/t, np.max(np.abs(o-ref)/np.abs(ref))) + ("   vs rtol1e-12: %.2e" % np.max(np.abs(o-ref12)/np.abs(ref12)) if ref12 is not None else ""))
```

Observed:

`12_hole_micro.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`16_equiv_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`13_timing_pre.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

**Reviewer (CONFIRMED (narrowed))**: The hole sum is real work and the saving is real, but the hunter's headline 1.37x is conditional on leading-log-1 being applied, while the claim as written reads as standalone. Standalone on the unmodified code I measure 15 to 17% (1.15x to 1.17x), and the hole sum takes 27 to 34% of the call, not the 57% the hunter profiled after leading-log-1. The gain is also gated by memory: the two tables take 16*n*d^2 bytes (n = number of distinct incoming energies, 122 for 61 biases), which is 32 MB at 7 spins but about 0.5 GB at 9 spins with 61 biases and 2 GB at 241 biases, so above the 256 MB cap the on-the-fly kernel stays. The size the audit flagged as expensive (9 spins, 392 s on d22) therefore gets nothing, and the saving is a few seconds at 6 to 8 spins. Not threaded, not GPU, not complex64, so it does not hit any 'Do NOT' item in kondo.md. Numbers: standalone the end-to-end change is 7e-16 to 9e-16, which is roundoff, so numbers_change=false holds when it is applied alone. The hunter's 2.7e-11 belongs to leading-log-1. Struck: the implication that 1.37x is available at the large sizes, and the 57% share.

Reviewer's own evidence:

```
05_holepre_alone.out: 'n=6 current code: orig 2.82 s, +hole precompute 2.41 s (1.17x), max rel diff 7.0e-16 ; fraction of the call in _hole_sum 34%'; 'n=7 ... orig 13.52 s, +hole precompute 11.71 s (1.15x), max rel diff 9.4e-16 ; fraction ... 27%'. Combined with leading-log-1 (02_timing.out): 'n=7 ... +pre 4.04 (3.36x, rel diff 2.7e-11)' against real+block 5.52 s.
```

Reviewer's measurement: Own medians of 3 on d22 through the public call. Standalone on the current code: 6 spins 2.82 s to 2.41 s (1.17x, 15% saved), 7 spins 13.52 s to 11.71 s (1.15x, 13% saved), max relative difference 9.4e-16. On top of leading-log-1 at 7 spins: 5.52 s to 4.04 s (1.37x, 27% saved). Inside the gate only: no change at 9 spins.

Reviewer's scripts and outputs:

`05_holepre_alone.py`:

```python
"""Hole precompute alone on the CURRENT complex flow (no real-block change), public call, median of 3, 6 and 7 spins; plus the fraction of the call spent in _hole_sum."""
import numpy as np, warnings, time, sys
sys.path.insert(0, ".")
from dmrgpy import spinchain
import llhole
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
orig = L._hole_sum
cache = {}
spent = [0.]
def timed_orig(*a):
    t = time.perf_counter(); orig(*a); spent[0] += time.perf_counter()-t
def pre_hole(Q, w, e, D, Dlow, omega0, Gamma0, out):
    if Dlow >= 0.: return orig(Q, w, e, D, Dlow, omega0, Gamma0, out)
    key = (len(w), len(e), float(w[0]), float(w[-1]), float(omega0), float(Gamma0))
    if key not in cache: cache.clear(); cache[key] = llhole.precompute(w, e, omega0, Gamma0)
    R, W = cache[key]; llhole.hole_sum_pre(Q, R, W, D, out)
def build(n):
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = sum(1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 0.8*sc.Sz[i]*sc.Sz[i+1]) for i in range(n-1))
    h = h + sum(1e-4*(1+0.1*i)*sc.Sz[i] + 3e-5*sc.Sx[i] for i in range(n))
    sc.set_hamiltonian(h); return sc
eVs = np.linspace(-4e-3, 4e-3, 61) + 1.3e-5
call = lambda sc: sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, U=0.25, T=0., order="LL")[1]
for n in (6, 7):
    sc = build(n); call(sc)
    res = {}
    for lab, f in (("orig", orig), ("pre", pre_hole)):
        L._hole_sum = f; ts = []
        for _ in range(3):
            t = time.perf_counter(); o = call(sc); ts.append(time.perf_counter()-t)
        res[lab] = (np.median(ts), o, ts)
    L._hole_sum = timed_orig; spent[0] = 0.; t = time.perf_counter(); call(sc); tot = time.perf_counter()-t
    print("n=%d current code: orig %.2f s, +hole precompute %.2f s (%.2fx), max rel diff %.1e ; fraction of the call in _hole_sum %.0f%%" % (
        n, res["orig"][0], res["pre"][0], res["orig"][0]/res["pre"][0], np.max(np.abs(res["pre"][1]-res["orig"][1])/np.abs(res["orig"][1])), 100*spent[0]/tot), flush=True)
```

`02_timing.py`:

```python
"""Own timing: public get_kondo_spectrum(order='LL') old vs real+block (+hole precompute), median of 3, 7 spins 1/2 and 6 spins, real H."""
import numpy as np, warnings, time, sys
from dmrgpy import spinchain
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
def build(n, by=0.):
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = sum(1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 0.8*sc.Sz[i]*sc.Sz[i+1]) for i in range(n-1))
    h = h + sum(1e-4*(1+0.1*i)*sc.Sz[i] + 3e-5*sc.Sx[i] + by*sc.Sy[i] for i in range(n))
    sc.set_hamiltonian(h); return sc
eVs = np.linspace(-4e-3, 4e-3, 61) + 1.3e-5
def med(sc, f):
    L._joint_flow = f; ts = []
    for _ in range(3):
        t = time.perf_counter(); o = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, U=0.25, T=0., order="LL")[1]; ts.append(time.perf_counter()-t)
    return np.median(ts), o, ts
for n in (6, 7):
    sc = build(n); sc.get_kondo_spectrum(eVs[:3], site=0, Jrho_s=-0.05, U=0.25, T=0., order="LL")
    llreal.install(); L._joint_flow = old
    to, oo, ro = med(sc, old)
    llreal.USE_PRE = False; tr, orr, rr = med(sc, llreal.joint_flow)
    llreal.USE_PRE = True; tp, op, rp = med(sc, llreal.joint_flow)
    print("n=%d old %.2f %s | real+block %.2f %s (%.2fx, rel diff %.1e) | +pre %.2f %s (%.2fx, rel diff %.1e)" % (n, to, np.round(ro,2), tr, np.round(rr,2), to/tr, np.max(np.abs(orr-oo)/np.abs(oo)), tp, np.round(rp,2), to/tp, np.max(np.abs(op-oo)/np.abs(oo))), flush=True)
```

`05_holepre_alone.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`02_timing.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

**Suggested fix (hunter)**: Build R and Wt once at the top of _joint_flow (a helper like llhole.precompute(w,e,omega0,Gamma0), shape (d,n,d), C-contiguous) when 2*n*d*d*8 bytes is below a cap of a few hundred MB, and pass them to a variant of _hole_sum that reads them with one division per element; keep the present on-the-fly _hole_sum for the closed-form tail call (Dlow>=0) and for sizes above the cap. The hole sum itself changes only by roundoff (2e-16); the end-to-end 2.7e-11 is the same step-sequence movement as leading-log-1, since the two changes were timed together.

**Fix (reviewer)**: Worth doing after leading-log-1 and as a second step: build R and Wt in _joint_flow once (shape (d,n,d), C-contiguous), pass them to a hole-sum variant that does one division per element, keep the on-the-fly _hole_sum for the Dlow>=0 closed-form tail and above the cap. Its cap is a memory trade, so make the cap a module constant (HOLE_TABLE_BYTES) and say in the README that at 9 spins nothing changes. Pin it with a test that evaluates both kernels on one (Q,w,e,D) at roundoff (1e-14 relative). Do not describe it as a speedup of the 9-spin call.

### 10. A tip site with no spin gives NaN at every bias with a false Kondo-scale warning, for an exact 1.5708

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `leading-log` &middot; candidate `leading-log-3`

**Status**: FIXED in `fce3382` on 2026-11-12, by the reviewer's three edits rather than an early return: the stop needs a strict crossing, `max_coupling` is divided by the sample coupling only when it is positive, and a zero bare vertex gives a ratio of 0. NUMBERS CHANGE, for an impurity with no spin on its occupied states only: NaN with a false warning to 2*pi*4*U^2, 1.5708 at U=0.25, with no warning. Pinned by `tests/test_audit_2026_11_leadinglog.py::test_a_state_with_no_spin_gives_the_potential_scattering` and `::test_an_empty_orbital_in_a_sector_gives_the_second_order`.

The claim as reviewed: When the tip site carries no spin in any state (a vacuum or a filled orbital in a conserved sector, S=0 everywhere), order="LL" returns NaN at all biases with the false warning "the sample coupling reached order one at D = 2.2e+08 (the Kondo scale is above Gamma0)", plus two divide RuntimeWarnings, while order=2 and order=3 return the exact potential-scattering value 2*pi*4*U^2 = 1.5708. When S vanishes only on the occupied columns, dI/dV is right but info["coupling"] is NaN at every bias with a divide warning. It is a corner of the input (a site with no spin has no Kondo effect to compute) and nothing in the docstring, the README or kondo.md names it as outside the method.

**Where**: src/dmrgpy/kondospectrumtk/leadinglog.py:340-341 (ref = pole_threshold*smax with smax = 0) and :371 (g_old <= ref and g_new >= ref is true for 0 <= 0 at the first step); :463-466 (bare = max(bare, BARE_FLOOR*max(bare)) = 0, so flow/bare is 0/0)

WHAT THE CODE DOES. The stop criterion compares max|G_ab| with pole_threshold*max|S_ab|/2. With S = 0 on every state both are 0, so `g_old <= ref and g_new >= ref` holds at the very first step, brentq returns the start, D_stop = Dmax, and every bias within 10 D_stop of anything is NaN. Separately, when only the occupied columns of S vanish (S acts on excited states only) bare is zero over [f,i] where it matters, the ratio flow/bare is 0/0 and info["coupling"] is NaN at every bias with a RuntimeWarning, while the dI/dV itself is right (1.5708 = pure potential scattering, 17_degenerate_real_and_probes.out (b2)). EXPECTED: the exchange vertex is identically zero, so the flow does nothing and dI/dV = 2*pi*4*U^2, which is what order=2 and order=3 give. REPRO through the public call: Spinful_Fermionic_Chain(2) in the vacuum sector Nf=0 (the tip orbital is empty, Sx=Sy=Sz=0 in the only state). The same stand-in with Sx=Sy=Sz=0 on a 2-level ks gives NaN too. SIZE: 100% of the biases (3 of 3) NaN, against an exact 1.5708. WHY TESTS PASS: every test model has a nonzero spin on the tip site; an empty or doubly occupied orbital in a sector is the only way to reach it, and the other two orders hide it because they only multiply by S. It is a corner of the input space, not a physical regime (a spinless site has no Kondo effect to compute), which is why it is LOW: it is silent only in that the NaN is attributed to the Kondo scale.

The hunter's claim, before review: A tip site that carries no spin in any state makes order="LL" return NaN at every bias with a false "Kondo scale is above Gamma0" warning (D_stop = Dmax = 2.2e8), where order=2 and order=3 return the exact potential-scattering value 2*pi*4*U^2 = 1.5708.

Decisive output, quoted by the hunter:

```
order=2 dIdV=[1.57079633 1.57079633 1.57079633]  warnings=[]
order=3 dIdV=[1.57079633 1.57079633 1.57079633]  warnings=[]
order=LL dIdV=[nan nan nan]  warnings=['invalid value encountered in divide', 'invalid value encountered in scalar divide', 'leading_log_dIdV: the sample coupling reached order one at D = 2.2e+08']
expected pure potential scattering: 2*pi*4*U^2 = 1.5707963267948966
(b2 zero spin on occupied columns only) dIdV [1.57 1.57 1.57]  coupling [nan nan nan]
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`18_zero_spin_public.py`:

```python
"""get_kondo_spectrum on a chain whose tip site carries no spin in any state (vacuum sector, Nf=0): the exact answer is pure potential scattering, 2*pi*4*U^2 per direction sum; keywords probe."""
import numpy as np, warnings
from dmrgpy import fermionchain
fc = fermionchain.Spinful_Fermionic_Chain(2)
h = sum(-5e-3*(fc.Cdag[s]*fc.C[2+s] + fc.Cdag[2+s]*fc.C[s]) for s in (0, 1))
fc.set_hamiltonian(h); fc.set_conserved_sector(Nf=0)
eV = np.array([-2e-3, 5e-4, 2e-3])
for order in (2, 3, "LL"):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _, d = fc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=0., order=order)
    print("order=%s dIdV=%s  warnings=%s" % (order, d, sorted(set(str(x.message)[:70] for x in w))))
print("expected pure potential scattering: 2*pi*4*U^2 =", 2*np.pi*4*0.25**2)
# keywords at the public call
from dmrgpy import spinchain
sc = spinchain.Spin_Chain(["1/2", "1/2"]); sc.set_hamiltonian(1e-3*(sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1]) + 1e-4*sc.Sz[0])
for kw in (dict(), dict(T=0., Jrho=-0.05), dict(T=0., nkry=50), dict(T=0., rtol=1e-6), dict(T=0., kB=1.0), dict(T=0., Gamma0=-1e-6)):
    try:
        _, d = sc.get_kondo_spectrum(np.array([5e-4]), site=0, Jrho_s=-0.05, U=0.25, order="LL", **kw)
        print("kw", kw, "->", d)
    except Exception as ex: print("kw", kw, "->", type(ex).__name__, str(ex)[:90])
```

`17_degenerate_real_and_probes.py`:

```python
"""(a) prototype vs current with nocc>1 (B=0 trimer, real H, doublet GS; random-rotated doublet) and tip_vertex(cols=None);
(b) probes: info['coupling'] under a rotation inside the degenerate ground doublet; coupling when the occupied states carry no spin;
(c) keywords at the public call."""
import numpy as np, warnings
from types import SimpleNamespace
from common import spin, ks_from
import llreal
from dmrgpy.kondospectrumtk import leadinglog as L
warnings.simplefilter("ignore")
old = L._joint_flow
S = spin(.5); I2 = np.eye(2)
def at(a, i, n=3):
    o = np.eye(1)
    for j in range(n): o = np.kron(o, S[a] if j == i else I2)
    return o
Htri = 1e-3*sum(at(a, i)@at(a, i+1) for i in range(2) for a in range(3))
ks = ks_from(Htri, [at(a, 1) for a in range(3)])
eVs = np.linspace(-4e-3, 4e-3, 17) + 1.3e-5
print("trimer nocc =", int(np.sum(ks.p > 0)), " real block:", llreal.block_is_real(L._spin_ops(ks)))
L._joint_flow = old; a, ia = L.leading_log_dIdV(ks, eVs, -0.05, U=0.25, return_info=True)
llreal.USE_PRE = True; L._joint_flow = llreal.joint_flow; b, ib = L.leading_log_dIdV(ks, eVs, -0.05, U=0.25, return_info=True)
print("(a) nocc=2: max rel diff dIdV %.1e, coupling %.1e" % (np.max(np.abs(a-b)/np.abs(a)), np.max(np.abs(ia["coupling"]-ib["coupling"]))))
# tip_vertex with cols=None
w = np.array([-3e-3, 1e-4, 2e-3]); Us = np.array([.1, .2, -.3])
L._joint_flow = old; A0, _ = L.tip_vertex(ks, w, Us, -0.05)
L._joint_flow = llreal.joint_flow; A1, _ = L.tip_vertex(ks, w, Us, -0.05)
print("(a) tip_vertex(cols=None) shapes", A0.shape, A1.shape, " max abs diff / max abs: %.1e" % (np.max(np.abs(A0-A1))/np.max(np.abs(A0))))
# (b1) rotation inside the doublet: ks has p on the 2 lowest states; rotate those two eigenstates by a random unitary
L._joint_flow = old
rng = np.random.default_rng(5)
occ = np.nonzero(ks.p > 0)[0]
def rotate(ks, th, ph):
    U = np.eye(ks.dim, dtype=complex)
    c, s = np.cos(th), np.sin(th)
    U[np.ix_(occ, occ)] = np.array([[c, -s*np.exp(-1j*ph)], [s*np.exp(1j*ph), c]])
    return SimpleNamespace(e=ks.e, p=ks.p, dim=ks.dim, T=0., kB=ks.kB, Sx=U.conj().T@ks.Sx@U, Sy=U.conj().T@ks.Sy@U, Sz=U.conj().T@ks.Sz@U)
for g in (-0.05, -0.1):
    base = None
    for th, ph in ((0, 0), (0.4, 0.3), (np.pi/4, 1.1), (1.2, 2.5)):
        k2 = rotate(ks, th, ph)
        o, info = L.leading_log_dIdV(k2, eVs, g, U=0.25, return_info=True)
        if base is None: base = (o, info["coupling"])
        print("(b1) g=%+.2f rotation(th=%.2f,ph=%.2f): max rel diff dIdV %.1e   coupling max %.4f  max |coupling - coupling(unrotated)| %.1e   NaNs %d" % (g, th, ph, np.nanmax(np.abs(o-base[0])/np.abs(base[0])), info["coupling"].max(), np.max(np.abs(info["coupling"]-base[1])), np.isnan(o).sum()))
# (b2) occupied states with no spin on the tip site: tip on an empty orbital-like site: occupied column of S is zero
d = 4
Sz = np.zeros((d, d), complex); Sx = np.zeros((d, d), complex); Sy = np.zeros((d, d), complex)
Sx[1:3, 1:3] = S[0]; Sy[1:3, 1:3] = S[1]; Sz[1:3, 1:3] = S[2]
ks0 = SimpleNamespace(e=np.array([0., 1e-3, 1.5e-3, 3e-3]), p=np.array([1., 0, 0, 0]), dim=4, T=0., kB=8.6e-5, Sx=Sx, Sy=Sy, Sz=Sz)
with warnings.catch_warnings(record=True) as w_:
    warnings.simplefilter("always")
    o, info = L.leading_log_dIdV(ks0, np.array([-2e-3, 5e-4, 2e-3]), -0.05, U=0.25, return_info=True)
print("(b2) zero-spin occupied state: dIdV", o, " coupling", info["coupling"], " D_stop", info["D_stop"], " warnings:", [str(x.message)[:60] for x in w_])
ks00 = SimpleNamespace(e=np.array([0., 1e-3]), p=np.array([.5, .5]), dim=2, T=0., kB=8.6e-5, Sx=0*S[0], Sy=0*S[1], Sz=0*S[2])
with warnings.catch_warnings(record=True) as w_:
    warnings.simplefilter("always")
    o, info = L.leading_log_dIdV(ks00, np.array([-2e-3, 5e-4, 2e-3]), -0.05, U=0.25, return_info=True)
print("(b2') S=0 entirely: dIdV", o, " coupling", info["coupling"], " D_stop", info["D_stop"], " warnings:", [str(x.message)[:90] for x in w_])
```

Observed:

`18_zero_spin_public.out`:

```
order=2 dIdV=[1.57079633 1.57079633 1.57079633]  warnings=[]
order=3 dIdV=[1.57079633 1.57079633 1.57079633]  warnings=[]
order=LL dIdV=[nan nan nan]  warnings=['invalid value encountered in divide', 'invalid value encountered in scalar divide', 'leading_log_dIdV: the sample coupling reached order one at D = 2.2e+08']
expected pure potential scattering: 2*pi*4*U^2 = 1.5707963267948966
kw {} -> NotImplementedError order="LL" is T=0 only for now (pass T=0); see docs/new_methods/kondo_leading_log/README.m
kw {'T': 0.0, 'Jrho': -0.05} -> TypeError get_kondo_spectrum(mode="ED") got unexpected keyword argument(s): Jrho. They are either mo
kw {'T': 0.0, 'nkry': 50} -> TypeError get_kondo_spectrum(mode="ED") got unexpected keyword argument(s): nkry. They are either mo
kw {'T': 0.0, 'rtol': 1e-06} -> TypeError get_kondo_spectrum(mode="ED") got unexpected keyword argument(s): rtol. They are either mo
kw {'T': 0.0, 'kB': 1.0} -> [1.42855277]
kw {'T': 0.0, 'Gamma0': -1e-06} -> ValueError Gamma0 must be positive, it cuts the logs
```

`17_degenerate_real_and_probes.out`:

```
trimer nocc = 2  real block: True
(a) nocc=2: max rel diff dIdV 7.9e-11, coupling 3.7e-12
(a) tip_vertex(cols=None) shapes (3, 4, 8, 8) (3, 4, 8, 8)  max abs diff / max abs: 1.9e-11
(b1) g=-0.05 rotation(th=0.00,ph=0.00): max rel diff dIdV 0.0e+00   coupling max 0.0682  max |coupling - coupling(unrotated)| 0.0e+00   NaNs 0
(b1) g=-0.05 rotation(th=0.40,ph=0.30): max rel diff dIdV 2.5e-11   coupling max 0.0682  max |coupling - coupling(unrotated)| 1.5e-12   NaNs 0
(b1) g=-0.05 rotation(th=0.79,ph=1.10): max rel diff dIdV 5.1e-11   coupling max 0.0682  max |coupling - coupling(unrotated)| 2.2e-12   NaNs 0
(b1) g=-0.05 rotation(th=1.20,ph=2.50): max rel diff dIdV 2.1e-11   coupling max 0.0682  max |coupling - coupling(unrotated)| 1.1e-12   NaNs 0
(b1) g=-0.10 rotation(th=0.00,ph=0.00): max rel diff dIdV 0.0e+00   coupling max 0.2126  max |coupling - coupling(unrotated)| 0.0e+00   NaNs 0
(b1) g=-0.10 rotation(th=0.40,ph=0.30): max rel diff dIdV 6.7e-11   coupling max 0.2126  max |coupling - coupling(unrotated)| 8.1e-12   NaNs 0
(b1) g=-0.10 rotation(th=0.79,ph=1.10): max rel diff dIdV 1.2e-10   coupling max 0.2126  max |coupling - coupling(unrotated)| 1.4e-11   NaNs 0
(b1) g=-0.10 rotation(th=1.20,ph=2.50): max rel diff dIdV 7.4e-11   coupling max 0.2126  max |coupling - coupling(unrotated)| 8.8e-12   NaNs 0
(b2) zero-spin occupied state: dIdV [1.57079633 1.57079633 1.57079633]  coupling [nan nan nan]  D_stop None  warnings: ['invalid value encountered in divide']
(b2') S=0 entirely: dIdV [nan nan nan]  coupling [nan nan nan]  D_stop 230000000.00000015  warnings: ['invalid value encountered in scalar divide', 'invalid value encountered in divide', 'leading_log_dIdV: the sample coupling reached order one at D = 2.3e+08 (the Kondo scale is']
```

**Reviewer (CONFIRMED (narrowed))**: Documented? No: kondo.md and the method README (section 7) discuss a dimer ground state with no spin (an elastic vertex of zero, still a spin on other states) and the D_stop NaN mask, but not S=0 on the whole spectrum. It is not a NotImplementedError boundary either. Probe: I ran my own public call (01_zero_spin.py) on Spinful_Fermionic_Chain(2) in the Nf=0 sector and in Nf=4 (both orbitals full), and both give NaN for LL and 1.5708 for orders 2 and 3, so the reproduction holds for empty and for filled orbitals. Mechanism: smax=0 makes ref=0, so `g_old <= ref and g_new >= ref` is true at the first step and brentq returns the start (D_stop=Dmax), and the bare floor gives bare=0 so flow/bare is 0/0. The right answer: the exchange vertex is identically zero and the potential scattering of the tip does not flow (direct and exchange cancel, documented in tip_vertex), so the flow is trivial and dI/dV equals the order-2 value, as the other orders give. Severity stays LOW: no spin means no physical regime for this method, and a NaN is not a wrong number. It is still a false diagnosis. Struck: nothing, but I replace the hunter's fix, which skips the flow with an early return, by a smaller one.

Reviewer's own evidence:

```
01_zero_spin.out: '2 [1.57079633 1.57079633 1.57079633] []'; 'LL [nan nan nan] [invalid value encountered in divide, invalid value encountered in scalar divide, leading_log_dIdV: the sample coupling reached order one at D]'; 'Nf=4 LL [nan nan nan]'. With the minimal patch (06_zero_spin_fix.out): 'S=0 everywhere, patched: dIdV [1.57079633 1.57079633 1.57079633] coupling [0. 0. 0.] D_stop None warnings []'; 'occupied column zero, patched: ... coupling [0. 0. 0.]'; '22 passed'.
```

Reviewer's scripts and outputs:

`01_zero_spin.py`:

```python
import numpy as np, warnings
from dmrgpy import fermionchain
fc = fermionchain.Spinful_Fermionic_Chain(2)
h = sum(-5e-3*(fc.Cdag[s]*fc.C[2+s] + fc.Cdag[2+s]*fc.C[s]) for s in (0, 1))
fc.set_hamiltonian(h); fc.set_conserved_sector(Nf=0)
eV = np.array([-2e-3, 5e-4, 2e-3])
for order in (2, 3, "LL"):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _, d = fc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=0., order=order)
    print(order, d, sorted(set(str(x.message)[:60] for x in w)))
# a real-life variant: half filled single site 'spinless-like' sector Nf=2 on 2 sites => site 0 spin nonzero. Nf=4 doubly occ both
fc.set_conserved_sector(Nf=4)
for order in (2,"LL"):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _, d = fc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=0., order=order)
    print("Nf=4",order, d, sorted(set(str(x.message)[:60] for x in w)))
```

`06_zero_spin_fix.py`:

```python
"""Prototype of the minimal fix (strict crossing, guarded ratios) by exec of a patched copy of leadinglog.py into its namespace; check the zero-spin case and that the regression tests still pass."""
import numpy as np, warnings, sys, re
from dmrgpy.kondospectrumtk import leadinglog as L
from types import SimpleNamespace
src = open(L.__file__).read()
def rep(a, b):
    global src
    assert src.count(a) == 1, a
    src = src.replace(a, b)
rep("if g_old <= ref and g_new >= ref:", "if g_old <= ref and g_new > ref:")
rep("max_coupling/smax)}", "max_coupling/(smax if smax > 0. else 1.))}")
rep("flow/bare[None]", "flow/np.where(bare > 0., bare, np.inf)[None]")
src = src.replace("cache=True", "cache=False"); src = src.replace("from .conductance", "from dmrgpy.kondospectrumtk.conductance")
src = src.replace("np.maximum(bare, BARE_FLOOR*np.max(bare))", "np.maximum(bare, BARE_FLOOR*np.max(bare))")
exec(compile(src, "leadinglog_patched", "exec"), L.__dict__)
d = 4
Z = np.zeros((d, d), complex)
S = [np.zeros((2, 2), complex) for _ in range(3)]
Sx = Z.copy(); Sy = Z.copy(); Sz = Z.copy()
ks = SimpleNamespace(e=np.array([0., 1e-3, 1.5e-3, 3e-3]), p=np.array([1., 0, 0, 0]), dim=4, T=0., kB=8.6e-5, Sx=Sx, Sy=Sy, Sz=Sz)
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    o, info = L.leading_log_dIdV(ks, np.array([-2e-3, 5e-4, 2e-3]), -0.05, U=0.25, return_info=True)
print("S=0 everywhere, patched: dIdV", o, "coupling", info["coupling"], "D_stop", info["D_stop"], "warnings", [str(x.message)[:60] for x in w])
print("expected", 2*np.pi*4*0.25**2)
# occupied column zero only
spin = np.array([[0,.5],[.5,0]]); 
Sx = Z.copy(); Sx[1:3,1:3] = spin; Sy = Z.copy(); Sy[1,2] = -.5j; Sy[2,1] = .5j; Sz = Z.copy(); Sz[1,1] = .5; Sz[2,2] = -.5
ks.Sx, ks.Sy, ks.Sz = Sx, Sy, Sz
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    o, info = L.leading_log_dIdV(ks, np.array([-2e-3, 5e-4, 2e-3]), -0.05, U=0.25, return_info=True)
print("occupied column zero, patched: dIdV", o, "coupling", info["coupling"], "warnings", [str(x.message)[:60] for x in w])
import pytest
sys.exit(pytest.main(["-q", "-p", "no:cacheprovider", "<repo>/tests/test_kondo_leading_log.py"]))
```

`01_zero_spin.out`:

```
2 [1.57079633 1.57079633 1.57079633] []
3 [1.57079633 1.57079633 1.57079633] []
LL [nan nan nan] ['invalid value encountered in divide', 'invalid value encountered in scalar divide', 'leading_log_dIdV: the sample coupling reached order one at D']
Nf=4 2 [1.57079633 1.57079633 1.57079633] []
Nf=4 LL [nan nan nan] ['invalid value encountered in divide', 'invalid value encountered in scalar divide', 'leading_log_dIdV: the sample coupling reached order one at D']
```

`06_zero_spin_fix.out`:

```
S=0 everywhere, patched: dIdV [1.57079633 1.57079633 1.57079633] coupling [0. 0. 0.] D_stop None warnings []
expected 1.5707963267948966
occupied column zero, patched: dIdV [1.57079633 1.57079633 1.57079633] coupling [0. 0. 0.] warnings []
......................                                                   [100%]
=============================== warnings summary ===============================
tests/test_kondo_leading_log.py::test_get_kondo_spectrum_order_LL
  <repo>/src/dmrgpy/manybodychain.py:24: PendingDeprecationWarning: the matrix subclass is not the recommended way to represent matrices or deal with linear algebra (see https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html). Please adjust your code to use regular ndarray.
    one = np.matrix(np.identity(3))

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
22 passed, 1 warning in 10.81s
```

**Suggested fix (hunter)**: When smax == 0 (the stacked spin operators are identically zero) skip the flow and return the bare vertex (dI/dV = 2*pi*T0^2*4*U^2 summed over both directions, i.e. the order=2 value), with info D_stop=None and coupling = 0. Make the crossing test strict (g_new > ref) so a vertex that is not above the threshold never triggers a stop, and guard the coupling ratio with where(bare > 0) so a zero bare vertex gives ratio 0 instead of 0/0. Changes returned numbers only in this degenerate case (NaN to 1.5708).

**Fix (reviewer)**: The hunter's early return is not needed. Three edits in src/dmrgpy/kondospectrumtk/leadinglog.py give the right answer without a special branch. In _joint_flow, make the crossing strict, `g_old <= ref and g_new > ref`, so a vertex that is not above the threshold never triggers a stop. In _joint_flow's returned info, divide max_coupling by `smax if smax > 0. else 1.`. In leading_log_dIdV, take the ratio against `np.where(bare > 0., bare, np.inf)`, so a zero bare vertex gives ratio 0, not 0/0. I prototyped exactly these (06_zero_spin_fix.py): S=0 gives dI/dV=1.5708, coupling 0, D_stop None, no warnings, and tests/test_kondo_leading_log.py still passes 22 of 22. Pin it with a test on a stand-in ks with S=0 (and one with S zero on the occupied column only) that expects 2*pi*4*U^2 and no warnings, plus the public Nf=0 repro. Also say in the docstring that a site without spin gives the potential-scattering term.

### 11. The two-time G runs six ket trajectories per row that three bra trajectories replace: 1.8x to 3.1x

`performance` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `correlator-twotime` &middot; candidate `correlator-twotime-1`

**Status**: FIXED in `4d1abbc` on 2026-11-12. G is built from the three bra trajectories of S_l|GS>, run once before the t2 loop; Gx keeps its three trajectories per row, so a row costs three TDVP runs instead of nine. Public `get_kondo_spectrum(order=3, submode="KPM")`, 8 sites, `maxm=60`, 21 biases: v3 72.3 s to 26.5 s (local), `"python"` 56.5 s to 21.8 s (d22). NUMBERS CHANGE below saturation: at saturation the term moves by 3.5e-11 (v3) and 3.4e-13 (`"python"`); below it the bra evolved backward and the ket evolved forward are two TDVP approximations, and the term moves by the truncation error, 8.7e-4 and 7.3e-5 of the ED peak on v3 (8 sites, `maxm` 4 and 8) and 1.8e-3 and 1.5e-4 on `"python"` (10 sites). The error against the exact term is the same size before and after: 1.42e-2 to 1.48e-2 and 2.8e-4 to 2.1e-4 on v3, 1.22e-2 to 1.24e-2 and 5.6e-4 to 6.0e-4 on `"python"`, four points, neither closer. Pinned by `tests/test_audit_2026_11_twotime.py::test_G_from_bra_trajectories_equals_ket_trajectories` (the old construction kept in the test, v3 and `"python"`, real and complex field).

The claim as reviewed: The two-time Kondo term's G block runs six ket TDVP trajectories per t2 row that three bra trajectories of S_l|GS> (computed once) replace, since <S_l GS|U(tau)|branch> = <U(-tau) S_l GS|branch>. The term is unchanged to 1e-11 relative at saturation, and the call is 1.8x (v3, 5 sites) to 3.1x (python, 5 sites) faster in my measurement, against the hunter's 2.2x to 2.5x at 8 sites.

**Where**: src/dmrgpy/kondospectrumtk/dmrgtwotime.py:179-198 (the j,k loop in _levi_civita_coeff_G_batches_dmrg that calls _tdvp_trajectory(chain, Hop, branch, dtau, n_tau_half) six times per t2 checkpoint)

The code evolves branch = S_k psi_j(t2) along tau for each of the six (j,k) pairs with j != k, at every t2 row, and overlaps each checkpoint with the fixed reference S_l|GS>. But <S_l GS| exp(-iH'tau) |branch> = <exp(+iH'tau) S_l GS | branch>, and exp(+iH'tau) S_l|GS> does not depend on t2. So three trajectories of S_l|GS> (one per l, both time directions, computed once for the whole calculation) replace 6 x n_t2 trajectories, and each G row costs only 6 x (2 n_tau_half+1) MPS overlaps. The Gx half still needs its three trajectories per row (its bra psi_l(-s) depends on s), so a row goes from 9 trajectories (270 steps) to 3 (90 steps). Profiled share of tdvp_step in the call is 94% (v3, 8 sites) and 96% (python, 8 sites), so the saving is about 60% of the whole public call, and it grows with n_t2_half because the three bra trajectories are a fixed cost. Numbers: at saturation (8 sites, maxm 60 above the exact bond dimension) G differs from the old G by 7.6e-11 (max |G| 0.75), Gx is bit identical, and the eV-swept term differs by 4.1e-12 relative to its peak. At truncation (python, 10 sites, maxm=8) G differs by 7.4e-3 and the term by 5.3e-5 relative, and against ED the new G is slightly closer (max |G-G_ED| 7.3e-3 against 9.7e-3 for the old one), so the truncation error is of the same order, neither route is better in kind. Medians of 3 repeats on identical inputs, threads pinned: v3 4 sites 15.64 s to 6.21 s (2.52x); v3 8 sites 38.29 s to 15.39 s (2.49x); python on d22 8 sites 18.00 s to 8.29 s (2.17x); python 10 sites maxm=8 23.8 s to 11.0 s (2.17x, one run).

The hunter's claim, before review: The two-time Kondo term runs six TDVP ket trajectories per t2 row for G that are not needed: G can be assembled from three bra trajectories computed once, cutting the TDVP steps from 270 to 90 per row and the call by 2.17x on python (d22, 8 sites, 18.0 s to 8.29 s) and 2.49x on v3 (8 sites, 38.3 s to 15.4 s), with the term unchanged to 4e-12 at saturation.

Decisive output, quoted by the hunter:

```
v3 n=8: max |G_old-G_new| 7.577812020918761e-11 max|G| 0.7500000000000036 | max |Gx_old-Gx_new| 0.0 | term max rel diff 4.055259112424148e-12 | old median s 38.29 | new median s 15.39 | speedup 2.4886301443122547
python d22 n=8: term max rel diff 3.1735900892862974e-13 | old median s 18.0 | new median s 8.29 | speedup 2.1731590778854812
truncated python n=10 maxm=8: max |G_old-G_new| 0.007356 term max rel diff 5.33e-05; vs ED: max|G_old-G_ED| 0.009738868292891696 max|G_new-G_ED| 0.00729242562805752
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`proto_g_bra.py`:

```python
"""Prototype: G(t2,tau) from three bra trajectories computed once, not six ket
trajectories per t2. <Sl gs|U(tau)|Sk psi_j(t2)> = <U(-tau) Sl gs | Sk psi_j(t2)>."""
import numpy as np
from dmrgpy.kondospectrumtk import dmrgtwotime as dt
from dmrgpy.kondospectrumtk.dmrgtwotime import _tdvp_trajectory, _EPS3, _AXES

def batches_new(chain, site, dt2, n_t2_half, dtau, n_tau_half, t2_batch=1):
    ops = {name: getattr(chain, name) for name in _AXES}
    E0 = chain.gs_energy()
    Hop = chain.toMPO(chain.hamiltonian - E0)
    gs = chain.get_gs()
    ref_wf = {l: ops[l][site]*gs for l in _AXES}
    # once: Sl|GS> at times -tau, tau on the tau grid (index n <-> tau_n = -tau_n of the walk)
    bra = {}
    for l in _AXES:
        _, wfs = _tdvp_trajectory(chain, Hop, ref_wf[l], dtau, n_tau_half)
        bra[l] = wfs[::-1]          # entry n is the state at time -tau_n
    t2_times, t2_wfs_by_j = {}, {}
    for j in _AXES:
        times, wfs = _tdvp_trajectory(chain, Hop, ops[j][site]*gs, dt2, n_t2_half)
        t2_times[j] = times; t2_wfs_by_j[j] = wfs
    n_t2 = len(t2_times[_AXES[0]])
    ntau = 2*n_tau_half + 1
    tau_grid = dtau*np.arange(-n_tau_half, n_tau_half+1)
    for it2 in range(n_t2):
        t2 = t2_times[_AXES[0]][it2]
        G_row = np.zeros((1, ntau), dtype=complex)
        for jj, j in enumerate(_AXES):
            psi = t2_wfs_by_j[j][it2]
            for kk, k in enumerate(_AXES):
                if jj == kk: continue
                branch = ops[k][site]*psi
                for ll, l in enumerate(_AXES):
                    c = _EPS3[jj, kk, ll]
                    if c == 0.: continue
                    G_row[0, :] += c*np.array([b.dot(branch) for b in bra[l]])
        mirror = n_t2 - 1 - it2
        kets = {(j, k): ops[k][site]*t2_wfs_by_j[j][mirror] for j in _AXES for k in _AXES if j != k}
        Gx_row = np.zeros((1, ntau), dtype=complex)
        for ll, l in enumerate(_AXES):
            _, bras = _tdvp_trajectory(chain, Hop, t2_wfs_by_j[l][mirror], dtau, n_tau_half)
            bras = bras[::-1]
            for jj, j in enumerate(_AXES):
                for kk, k in enumerate(_AXES):
                    c = _EPS3[jj, kk, ll]
                    if c == 0.: continue
                    Gx_row[0, :] += c*np.array([b.dot(kets[(j, k)]) for b in bras])
        yield np.array([t2]), tau_grid, G_row, Gx_row
```

`06_g_bra_check.py`:

```python
import sys, time, warnings, numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0,".")
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import dmrgtwotime as dt
import proto_g_bra as P
G=2.0; MUB=5.7883818066e-5
v=sys.argv[1] if sys.argv[1]=="python" else int(sys.argv[1]); n=int(sys.argv[2]); nt2=int(sys.argv[3]); rep=int(sys.argv[4])
def build():
    sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version=v)
    h=G*MUB*10.*sc.Sz[0]
    for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h); sc.maxm=int(sys.argv[5]) if len(sys.argv)>5 else 60; sc.nsweeps=12
    return sc
old=dt._levi_civita_coeff_G_batches_dmrg
eV=np.linspace(-4e-3,4e-3,21)
kw=dict(omega0=2e-3,Gamma0=2e-4,dt2=1000.,n_t2_half=nt2,dtau=1000.,n_tau_half=15)
sc=build(); sc.get_gs()
# rows
ro=list(old(sc,0,1000.,nt2,1000.,15)); rn=list(P.batches_new(sc,0,1000.,nt2,1000.,15))
print("max |G_old-G_new|",max(np.abs(a[2]-b[2]).max() for a,b in zip(ro,rn)),"max|G|",max(np.abs(a[2]).max() for a in ro))
print("max |Gx_old-Gx_new|",max(np.abs(a[3]-b[3]).max() for a,b in zip(ro,rn)))
res={}
for name,f in (("old",old),("new",P.batches_new)):
    ts=[]
    for r in range(rep):
        dt._levi_civita_coeff_G_batches_dmrg=f
        t=time.time(); out=dt.two_time_kondo_term_dmrg(sc,0,eV,**kw); ts.append(time.time()-t)
    res[name]=(out,sorted(ts)[len(ts)//2],ts)
dt._levi_civita_coeff_G_batches_dmrg=old
print("term max rel diff",np.abs(res["old"][0]-res["new"][0]).max()/np.abs(res["old"][0]).max(),"peak",np.abs(res["old"][0]).max())
for k in res: print(k,"median s",round(res[k][1],2),"all",[round(x,2) for x in res[k][2]])
print("speedup",res["old"][1]/res["new"][1])
```

`10_trunc_vs_ed.py`:

```python
import sys, warnings, numpy as np
warnings.filterwarnings("ignore"); sys.path.insert(0,".")
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import dmrgtwotime as dt
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.edtwotimeref import _levi_civita_coeff_G_chunk
import proto_g_bra as P
G=2.0; MUB=5.7883818066e-5
n=10; maxm=int(sys.argv[1])
sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version="python")
h=G*MUB*10.*sc.Sz[0]
for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=maxm; sc.nsweeps=12; sc.get_gs()
ks=KondoSpectrum(sc,0,0.)
tau=1000.*np.arange(-15,16)
ro=list(dt._levi_civita_coeff_G_batches_dmrg(sc,0,1000.,1,1000.,15)); rn=list(P.batches_new(sc,0,1000.,1,1000.,15))
eo=en=0.
for a,b in zip(ro,rn):
    ref=_levi_civita_coeff_G_chunk(ks,a[0],tau)[0]
    eo=max(eo,np.abs(a[2][0]-ref).max()); en=max(en,np.abs(b[2][0]-ref).max())
print("maxm",maxm,"GS energy err vs ED",sc.gs_energy()-ks.e0 if hasattr(ks,'e0') else None)
print("max|G_old-G_ED|",eo,"max|G_new-G_ED|",en)
```

`03_prof_twotime.py`:

```python
import sys, time, cProfile, pstats, numpy as np, warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain
G=2.0; MUB=5.7883818066e-5
v=sys.argv[1] if sys.argv[1]=="python" else int(sys.argv[1]); n=int(sys.argv[2]); nt2=int(sys.argv[3]); nev=int(sys.argv[4])
sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version=v)
h=G*MUB*10.*sc.Sz[0]
for i in range(n-1):
    h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=10
sc.get_gs()
eV=np.linspace(-4e-3,4e-3,nev)
from dmrgpy.kondospectrumtk.dmrgtwotime import two_time_kondo_term_dmrg
pr=cProfile.Profile(); t=time.time(); pr.enable()
r=two_time_kondo_term_dmrg(sc,0,eV,omega0=2e-3,Gamma0=5e-6,dt2=float(sys.argv[5]),n_t2_half=nt2,dtau=float(sys.argv[6]),n_tau_half=15)
pr.disable(); tot=time.time()-t
print("version",v,"n",n,"nt2half",nt2,"nev",nev,"total",round(tot,2))
pstats.Stats(pr).sort_stats("cumulative").print_stats(14)
pstats.Stats(pr).sort_stats("tottime").print_stats(8)
```

Observed:

`06_g_bra_v3_n4.out`:

```
max |G_old-G_new| 2.6360281894337346e-11 max|G| 0.7499999999999997
max |Gx_old-Gx_new| 0.0
term max rel diff 1.0704378963806254e-12 peak 1.1260394125836568
old median s 15.64 all [15.64, 15.63, 15.75]
new median s 6.21 all [6.83, 6.19, 6.21]
speedup 2.5192859761077098
```

`06_g_bra_v3_n8.out`:

```
max |G_old-G_new| 7.577812020918761e-11 max|G| 0.7500000000000036
max |Gx_old-Gx_new| 0.0
term max rel diff 4.055259112424148e-12 peak 0.9955504432202749
old median s 38.29 all [38.73, 38.29, 36.74]
new median s 15.39 all [15.39, 16.17, 15.27]
speedup 2.4886301443122547
```

`08_g_bra_python_n8_d22.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`09_g_bra_python_n10_trunc_d22.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`10_trunc_vs_ed.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`07_prof_python_n8_d22.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`03_prof_v3_n8_fine.out`:

```
version 3 n 8 nt2half 3 nev 41 total 55.94
         1371740 function calls (1333728 primitive calls) in 55.944 seconds

   Ordered by: cumulative time
   List reduced from 162 to 18 due to restriction <18>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000   55.944   55.944 <repo>/src/dmrgpy/kondospectrumtk/dmrgtwotime.py:220(two_time_kondo_term_dmrg)
        1    0.005    0.005   55.944   55.944 <repo>/src/dmrgpy/kondospectrumtk/twotime.py:163(kondo_term_from_two_time)
        8    0.001    0.000   53.822    6.728 <repo>/src/dmrgpy/kondospectrumtk/dmrgtwotime.py:260(batches)
        8    0.020    0.003   53.822    6.728 <repo>/src/dmrgpy/kondospectrumtk/dmrgtwotime.py:151(_levi_civita_coeff_G_batches_dmrg)
       66    0.010    0.000   52.945    0.802 <repo>/src/dmrgpy/kondospectrumtk/dmrgtwotime.py:89(_tdvp_trajectory)
      132    0.037    0.000   52.712    0.399 <repo>/src/dmrgpy/kondospectrumtk/dmrgtwotime.py:115(_walk)
     1908   52.428    0.027   52.428    0.027 {built-in method dmrgpy.mpscpp3._dmrgcpp.tdvp_step}
38290/287    1.694    0.000    2.070    0.007 <repo>/src/dmrgpy/kondospectrumtk/twotime.py:110(K_F)
      574    0.002    0.000    2.044    0.004 <home>
      574    0.002    0.000    2.041    0.004 <home>
       82    0.016    0.000    1.512    0.018 {built-in method scipy.integrate._quadpack._qagpe}
    24108    0.053    0.000    1.297    0.000 <repo>/src/dmrgpy/kondospectrumtk/twotime.py:134(<lambda>)
    13608    0.029    0.000    0.716    0.000 <repo>/src/dmrgpy/kondospectrumtk/twotime.py:136(<lambda>)
     2670    0.002    0.000    0.690    0.000 <repo>/src/dmrgpy/mps.py:26(dot)
     2670    0.002    0.000    0.688    0.000 <repo>/src/dmrgpy/manybodychain.py:1035(overlap)
     2670    0.003    0.000    0.686    0.000 <repo>/src/dmrgpy/mpsalgebra.py:169(overlap)
     2670    0.003    0.000    0.675    0.000 <repo>/src/dmrgpy/mpsalgebra.py:183(overlap_dmrg)
     2670    0.669    0.000    0.669    0.000 {built-in method dmrgpy.mpscpp3._dmrgcpp.overlap}
```

**Reviewer (CONFIRMED)**: I tried to break it on inputs the hunter did not use. The identity needs only a Hermitian Hop, and Hop = H - E0 is Hermitian. I tested a complex Hamiltonian (field along 0.6 Sx + 0.8 Sy, which makes the MPO complex) on v3 and python, and a spin-1 chain on python. G agrees to 4e-12 to 5e-11 against |G| of 0.75 (2.0 for spin 1), and Gx is bit identical because that block is untouched. The mirror and index logic is right: tdvp walk entry n is the state at time -tau_n, so the reversed list gives U(-tau_n) ref, matching the old overlap at tau_n. Nothing in kondo.md or CLAUDE.md forbids it. The rejected optimizations listed there (GPU float64, complex64, numba threads, TASK batching) are unrelated. The only mention is the '1.5x the TDVP cost' note for Gx, which becomes stale after the fix: a row goes from 9 trajectories to 3 plus the 3 for Gx, i.e. Gx then costs as much as G. The hunter's truncated-regime number (term moves 5e-5 relative at maxm=8) I did not re-run. It is a different but equally valid TDVP error, not a regression. I did not test T>0, n_gs>1 or the fermion chains. The change sits below the n_gs and class loops and reads only Sx/Sy/Sz[site], the chain's Hamiltonian and gs, so I do not expect a difference. My timings are lower in speedup than the hunter's at v3 because my chain is 5 sites with n_t2_half=4: the three bra trajectories are a fixed cost, so the gain grows with n_t2_half. The v3 median was taken on a loaded shared machine.

Reviewer's own evidence:

```
01_check.out: '3 4 xy 1/2 G 4.56e-11 |G| 0.75 Gx 0.0'; 'python 4 xy 1/2 G 1.34e-11 |G| 0.75 Gx 0.0'; 'python 3 z 1 G 2.84e-12 |G| 2.0 Gx 0.0'. 02_time.out: '3 rel diff 8.6e-12 old 27.17 new 14.78 x 1.84'; 'python rel diff 7.6e-13 old 57.7 new 18.6 x 3.11'.
```

Reviewer's measurement: Own medians of 3 repeats of two_time_kondo_term_dmrg, 5 sites, n_t2_half=4, n_tau_half=15, maxm=60, threads pinned: v3 27.17 s to 14.78 s (1.84x, 46% saved), python 57.7 s to 18.6 s (3.1x, 68% saved, local, loaded machine, noisy). Max relative difference of the term 8.6e-12 (v3) and 7.6e-13 (python). I did not re-time on d22 or at 8 sites.

Reviewer's scripts and outputs:

`01_check.py`:

```python
import sys, time, warnings, numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0,".")
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import dmrgtwotime as dt
import proto_g_bra as P
MUB=5.7883818066e-5
def build(v,n,field,spin="1/2"):
    sc=spinchain.Spin_Chain([spin]*n, itensor_version=v)
    h=0
    for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    if field=="z": h=h+2*MUB*10*sc.Sz[0]
    if field=="xy": h=h+2*MUB*10*(0.6*sc.Sx[0]+0.8*sc.Sy[0])+1e-4*sc.Sz[1]
    if field=="none": h=h+1e-6*sc.Sz[0]
    sc.set_hamiltonian(h); sc.maxm=60; sc.nsweeps=14
    return sc
for v,n,field,spin in [(3,4,"z","1/2"),(3,4,"xy","1/2"),("python",4,"xy","1/2"),("python",3,"z","1")]:
    sc=build(v,n,field,spin); sc.get_gs()
    ro=list(dt._levi_civita_coeff_G_batches_dmrg(sc,0,1000.,2,1000.,6))
    rn=list(P.batches_new(sc,0,1000.,2,1000.,6))
    print(v,n,field,spin,"G",max(np.abs(a[2]-b[2]).max() for a,b in zip(ro,rn)),"|G|",max(np.abs(a[2]).max() for a in ro),"Gx",max(np.abs(a[3]-b[3]).max() for a,b in zip(ro,rn)),flush=True)
```

`02_time.py`:

```python
import sys, time, warnings, numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0,".")
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import dmrgtwotime as dt
import proto_g_bra as P
MUB=5.7883818066e-5
old=dt._levi_civita_coeff_G_batches_dmrg
eV=np.linspace(-4e-3,4e-3,21)
kw=dict(omega0=2e-3,Gamma0=2e-4,dt2=1000.,n_t2_half=4,dtau=1000.,n_tau_half=15)
for v in (3,"python"):
    sc=spinchain.Spin_Chain(["1/2"]*5, itensor_version=v)
    h=2*MUB*10*sc.Sz[0]
    for i in range(4): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h); sc.maxm=60; sc.nsweeps=14; sc.get_gs()
    res={}
    for name,f in (("old",old),("new",P.batches_new)):
        ts=[]
        for r in range(3):
            dt._levi_civita_coeff_G_batches_dmrg=f
            t=time.time(); out=dt.two_time_kondo_term_dmrg(sc,0,eV,**kw); ts.append(time.time()-t)
        res[name]=(out,sorted(ts)[1])
    dt._levi_civita_coeff_G_batches_dmrg=old
    print(v,"rel diff",np.abs(res["old"][0]-res["new"][0]).max()/np.abs(res["old"][0]).max(),"old",res["old"][1],"new",res["new"][1],"x",res["old"][1]/res["new"][1],flush=True)
```

`proto_g_bra.py`:

```python
"""Prototype: G(t2,tau) from three bra trajectories computed once, not six ket
trajectories per t2. <Sl gs|U(tau)|Sk psi_j(t2)> = <U(-tau) Sl gs | Sk psi_j(t2)>."""
import numpy as np
from dmrgpy.kondospectrumtk import dmrgtwotime as dt
from dmrgpy.kondospectrumtk.dmrgtwotime import _tdvp_trajectory, _EPS3, _AXES

def batches_new(chain, site, dt2, n_t2_half, dtau, n_tau_half, t2_batch=1):
    ops = {name: getattr(chain, name) for name in _AXES}
    E0 = chain.gs_energy()
    Hop = chain.toMPO(chain.hamiltonian - E0)
    gs = chain.get_gs()
    ref_wf = {l: ops[l][site]*gs for l in _AXES}
    # once: Sl|GS> at times -tau, tau on the tau grid (index n <-> tau_n = -tau_n of the walk)
    bra = {}
    for l in _AXES:
        _, wfs = _tdvp_trajectory(chain, Hop, ref_wf[l], dtau, n_tau_half)
        bra[l] = wfs[::-1]          # entry n is the state at time -tau_n
    t2_times, t2_wfs_by_j = {}, {}
    for j in _AXES:
        times, wfs = _tdvp_trajectory(chain, Hop, ops[j][site]*gs, dt2, n_t2_half)
        t2_times[j] = times; t2_wfs_by_j[j] = wfs
    n_t2 = len(t2_times[_AXES[0]])
    ntau = 2*n_tau_half + 1
    tau_grid = dtau*np.arange(-n_tau_half, n_tau_half+1)
    for it2 in range(n_t2):
        t2 = t2_times[_AXES[0]][it2]
        G_row = np.zeros((1, ntau), dtype=complex)
        for jj, j in enumerate(_AXES):
            psi = t2_wfs_by_j[j][it2]
            for kk, k in enumerate(_AXES):
                if jj == kk: continue
                branch = ops[k][site]*psi
                for ll, l in enumerate(_AXES):
                    c = _EPS3[jj, kk, ll]
                    if c == 0.: continue
                    G_row[0, :] += c*np.array([b.dot(branch) for b in bra[l]])
        mirror = n_t2 - 1 - it2
        kets = {(j, k): ops[k][site]*t2_wfs_by_j[j][mirror] for j in _AXES for k in _AXES if j != k}
        Gx_row = np.zeros((1, ntau), dtype=complex)
        for ll, l in enumerate(_AXES):
            _, bras = _tdvp_trajectory(chain, Hop, t2_wfs_by_j[l][mirror], dtau, n_tau_half)
            bras = bras[::-1]
            for jj, j in enumerate(_AXES):
                for kk, k in enumerate(_AXES):
                    c = _EPS3[jj, kk, ll]
                    if c == 0.: continue
                    Gx_row[0, :] += c*np.array([b.dot(kets[(j, k)]) for b in bras])
        yield np.array([t2]), tau_grid, G_row, Gx_row
```

`01_check.out`:

```
3 4 z 1/2 G 3.928955292672254e-12 |G| 0.7500000000000011 Gx 0.0
3 4 xy 1/2 G 4.562734861217345e-11 |G| 0.7500000000000013 Gx 0.0
python 4 xy 1/2 G 1.3376054584725167e-11 |G| 0.7499999999999979 Gx 0.0
python 3 z 1 G 2.8379948610919545e-12 |G| 2.000000000000001 Gx 0.0
```

`02_time.out`:

```
3 rel diff 8.575113021032527e-12 old 27.169398546218872 new 14.783569097518921 x 1.8378105021187763
python rel diff 7.569166944769097e-13 old 57.6772358417511 new 18.560473203659058 x 3.1075304605046634
```

**Suggested fix (hunter)**: In _levi_civita_coeff_G_batches_dmrg, before the t2 loop compute bra[l] = reversed(_tdvp_trajectory(chain, Hop, ref_wf[l], dtau, n_tau_half)[1]) for the three l (entry n is the state at time -tau_n), and in the j,k loop replace the per-pair _tdvp_trajectory call by G_row += c * [b.dot(branch) for b in bra[l]], with l the third axis of (j,k). Keep the Gx block as it is. At saturation the numbers move by roundoff (1e-11 in G, 4e-12 in the term); at truncated maxm the term moves by 5e-5 relative, which is the truncation error of the TDVP trajectories either way. A further saving for Gx exists only when dtau is a multiple of dt2 (the bras psi_l(-s-tau) then lie on one long trajectory), which is a restriction on the grids and is left as a lead.

**Fix (reviewer)**: The hunter's fix is right. In _levi_civita_coeff_G_batches_dmrg (src/dmrgpy/kondospectrumtk/dmrgtwotime.py), before the t2 loop compute bra[l] = _tdvp_trajectory(chain, Hop, ref_wf[l], dtau, n_tau_half)[1][::-1] for the three l, and in the j,k loop replace the per-pair trajectory with G_row += c*[b.dot(branch) for b in bra[l]]. tau_grid_row must then be built from dtau*np.arange(-n_tau_half, n_tau_half+1) because the trajectory call that set it is gone (proto_g_bra.py does this). Keep the Gx block. Pin it with a test that compares the new G against a reference built from the old six-trajectory route at saturation on a 4-site chain (atol 1e-9), on v3 and python, with a complex field. Update the kondo.md '1.5x the TDVP cost' sentence and the function docstring ('no cheap way to batch').

### 12. es starting at 0 cuts the elastic line in half with no warning: 0.000 at zero bias for an exact 1.037

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `correlator-twotime` &middot; candidate `correlator-twotime-2`

**Status**: FIXED in `2393bc0` on 2026-11-12, in the reviewer's form. `secondorder_dc._cumulative_theta0_weight` warns when the correlator at the lowest `es` is above 1e-2 of its maximum, which means the lower edge cuts a line, and when the largest `es` is below max|eV|. The guide's `es` item says both. No number changes. Pinned by `tests/test_audit_2026_11_dcroutes.py` (`::test_es_from_zero_warns_and_cuts_the_elastic_line`, `::test_covered_grid_does_not_warn`, `::test_singlet_with_es_from_zero_does_not_warn`, `::test_es_below_max_eV_warns`).

The claim as reviewed: get_kondo_spectrum(mode=DMRG, correlator submode) at order=2, or order=3 with U=0, never checks the documented es precondition of second_order_dIdV_dc: with es[0]=0 (linspace(0,W,N)) the elastic line at w=0 is half cut and the zero-bias value is 0.000 against an exact 1.037 (max error 1.037 on a 4.71 peak; 4-site Heisenberg chain, 10 T on site 0, v3 KPM, delta=2e-5), with no warning; the user guide states only the max|eV| side of the requirement. A missing check plus a guide gap, not a wrong value under the documented contract.

**Where**: src/dmrgpy/kondospectrumtk/secondorder_dc.py:_cumulative_theta0_weight (cumulative_trapezoid starting at the first es point, np.interp clamping at both ends) and second_order_dIdV_dc; the only coverage check in the folder is the sum-rule warning in potentialdc.py:third_order_potential_dIdV_dc, which runs only for order=3 with U!=0

The cumulative weight starts at 0 at the lowest es point and np.interp clamps outside the grid, so the elastic line at w=0 (weight 0.165 of 0.75 on this chain, and it is what sets the zero-bias plateau) is dropped when es starts at 0 or above, and every transition above max(es) is dropped when max(es) < max|eV|. The docstring says es must start several delta below 0 and pass max|eV|, but nothing checks it for the second-order term, and get_kondo_spectrum(order=2) never reaches the potential-term sum-rule warning. Chain: 4 site Heisenberg impurity chain, B=10 T on site 0, J=1 meV, T=0, v3, KPM delta=2e-5. Anchor mode="ED" (T=0) at the same eV. Existing tests use es that start at -1e-3 and run past the last line, so they never leave the grid short. The Casimir sum rule S(S+1) is not what to test for the second-order term (es may legitimately stop at max|eV|), the exact precondition is es.min() < -several delta and es.max() >= max|eV| + several delta.

The hunter's claim, before review: second_order_dIdV_dc returns a wrong second-order spectrum with no warning when es does not start below w=0 or does not reach max|eV|: es starting at 0 (the natural linspace(0, W, N)) gives 0.0 at zero bias against an exact 1.037 (the whole elastic plateau lost, max error 1.037 on a 6.97 peak), and an es top below max|eV| drops 0.23 of a 4.71 peak.

Decisive output, quoted by the hunter:

```
es lower edge -0.001 maxdiff 0.423 value at eV=0: got 1.037 ED 1.037 nwarn 0
es lower edge 0.0 maxdiff 1.037 value at eV=0: got 0.0 ED 1.037 nwarn 0
es lower edge 0.0005 maxdiff 1.228 value at eV=0: got 0.0 ED 1.037 nwarn 0
es top 0.002 maxdiff 0.23058119843072422 peak 4.712388980384691 nwarn 0
es top 0.003 maxdiff 0.018460923364479687 peak 4.712388980384691 nwarn 0   (the 0.42 and 0.018 floors are the delta=2e-5 broadening)
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`14_es_lower.py`:

```python
import warnings, numpy as np
from dmrgpy import spinchain
G=2.0; MUB=5.7883818066e-5
sc=spinchain.Spin_Chain(["1/2"]*4, itensor_version=3)
h=G*MUB*10.*sc.Sz[0]
for i in range(3): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12
eV=np.linspace(-4e-3,4e-3,41)
ed=sc.get_kondo_spectrum(eV,order=2,U=0.,mode="ED",T=0)[1]
for lo in (-1e-3,-1e-4,0.0,5e-4):
    es=np.linspace(lo,1e-2,int((1e-2-lo)/5e-6))
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter("always")
        r=sc.get_kondo_spectrum(eV,order=2,U=0.,mode="DMRG",T=0,submode="KPM",delta=2e-5,es=es)[1]
    print("es lower edge",lo,"maxdiff",round(np.abs(r-ed).max(),3),"value at eV=0: got",round(r[20],3),"ED",round(ed[20],3),"nwarn",len(W))
```

`04_es_clamp.py`:

```python
import warnings, numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
G=2.0; MUB=5.7883818066e-5
sc=spinchain.Spin_Chain(["1/2"]*4, itensor_version=3)
h=G*MUB*10.*sc.Sz[0]
for i in range(3): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12
ks=KondoSpectrum(sc,0,0.)
w=ks.e; wt=sum(np.abs(S[:,0])**2 for S in (ks.Sx,ks.Sy,ks.Sz))
print("transitions with weight>1e-3:", [(round(e*1e3,3),round(x,3)) for e,x in zip(w,wt) if x>1e-3])
eV=np.linspace(-5e-3,5e-3,21)
ed=sc.get_kondo_spectrum(eV,order=2,U=0.,mode="ED",T=0)[1]
for top in (2e-3,3e-3,6e-3):
    es=np.linspace(-1e-3,top,int((top+1e-3)/1e-5))
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter("always")
        r=sc.get_kondo_spectrum(eV,order=2,U=0.,mode="DMRG",T=0,submode="KPM",delta=2e-5,es=es)[1]
    print("es top",top,"maxdiff",np.abs(r-ed).max(),"peak",ed.max(),"nwarn",len(W))
```

Observed:

`14_es_lower.out`:

```
es lower edge -0.001 maxdiff 0.423 value at eV=0: got 1.037 ED 1.037 nwarn 0
es lower edge -0.0001 maxdiff 0.423 value at eV=0: got 1.037 ED 1.037 nwarn 0
es lower edge 0.0 maxdiff 1.037 value at eV=0: got 0.0 ED 1.037 nwarn 0
es lower edge 0.0005 maxdiff 1.228 value at eV=0: got 0.0 ED 1.037 nwarn 0
```

`04_es_clamp.out`:

```
transitions with weight>1e-3: [(np.float64(0.0), np.float64(0.165)), (np.float64(0.423), np.float64(0.028)), (np.float64(0.984), np.float64(0.017)), (np.float64(1.225), np.float64(0.017)), (np.float64(1.236), np.float64(0.096)), (np.float64(1.61), np.float64(0.058)), (np.float64(1.814), np.float64(0.33)), (np.float64(1.963), np.float64(0.002)), (np.float64(2.561), np.float64(0.01)), (np.float64(2.606), np.float64(0.025)), (np.float64(3.166), np.float64(0.001))]
es top 0.002 maxdiff 0.23058119843072422 peak 4.712388980384691 nwarn 0
es top 0.003 maxdiff 0.018460923364479687 peak 4.712388980384691 nwarn 0
es top 0.006 maxdiff 0.018199191205371612 peak 4.712388980384691 nwarn 0
```

**Reviewer (CONFIRMED (narrowed))**: Reproduced with my own script. es lower edge -1e-3 and -1e-4 give 1.037 at eV=0 (equal to ED); -2e-5 (one delta) gives 0.958; 0.0 gives 0.000, max error 1.037; nwarn 0 everywhere. The requirement is documented in the second_order_dIdV_dc docstring and the 2026-09-24b Ruled out line says it holds, so no numerical defect is new. What survives: (a) no check exists for order=2, because the 2026-10 finding-2 sum-rule warning lives only in third_order_potential_dIdV_dc, which does not run for order=2 or U=0, and its own text already says es must include w=0; (b) docs/guide/kondo.md says only that es must cover transitions below max|eV|, not that it must start below w=0, so a user following the guide writes linspace(0,W,N) and loses the zero-bias plateau; (c) the loss is the full elastic weight. Struck: the top-side claim (0.23 of 4.71) is a plain violation of an explicit documented requirement, and in my run (es top 0.002 for max|eV|=0.004) the error hid under the 0.42 delta floor, so I did not re-establish 0.23. Struck: the hunter es-only lower check es.min() > -3*delta, since at B=0 (singlet ground state, no elastic weight) es from 0 is harmless (max error 0.025, same as es from -1e-3), so an es-only rule would warn falsely there.

Reviewer's own evidence:

```
01_es_edge.out: 'es -0.001 0.01 maxdiff 0.423 eV=0: 1.037 ED 1.037 nwarn 0' / 'es -2e-05 0.01 maxdiff 0.384 eV=0: 0.958' / 'es 0.0 0.01 maxdiff 1.037 eV=0: 0.000 ED 1.037 nwarn 0'. 02_edge_probe.out: 'B 10.0 es 0.0 0.01 S[0]/Smax 1 maxdiff vs ED 1.037' against 'B 0.0 es 0.0 0.01 S[0]/Smax 2.59e-07 maxdiff vs ED 0.025'.
```

Reviewer's scripts and outputs:

`01_es_edge.py`:

```python
import sys, warnings, numpy as np
sys.path.insert(0,"<repo>/src")
import dmrgpy; print(dmrgpy.__file__)
from dmrgpy import spinchain
G=2.0; MUB=5.7883818066e-5
sc=spinchain.Spin_Chain(["1/2"]*4, itensor_version=3)
h=G*MUB*10.*sc.Sz[0]
for i in range(3): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12
eV=np.linspace(-4e-3,4e-3,41)
ed=sc.get_kondo_spectrum(eV,order=2,U=0.,mode="ED",T=0)[1]
for lo,hi in ((-1e-3,1e-2),(-1e-4,1e-2),(-2e-5,1e-2),(0.0,1e-2),(-1e-3,2e-3),(-1e-3,3e-3)):
    es=np.linspace(lo,hi,int((hi-lo)/5e-6))
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter("always")
        r=sc.get_kondo_spectrum(eV,order=2,U=0.,mode="DMRG",T=0,submode="KPM",delta=2e-5,es=es)[1]
    print("es",lo,hi,"maxdiff",round(np.abs(r-ed).max(),3),"eV=0:",round(r[20],3),"ED",round(ed[20],3),"nwarn",len(W),"peak",round(ed.max(),3))
```

`02_edge_probe.py`:

```python
import sys, numpy as np
sys.path.insert(0,"<repo>/src")
from dmrgpy import spinchain
G=2.0; MUB=5.7883818066e-5
def build(B):
    sc=spinchain.Spin_Chain(["1/2"]*4, itensor_version=3)
    h=G*MUB*B*sc.Sz[0]
    for i in range(3): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12
    return sc
eV=np.linspace(-4e-3,4e-3,41)
for B in (10.,0.):
    sc=build(B)
    ed=sc.get_kondo_spectrum(eV,order=2,U=0.,mode="ED",T=0)[1]
    for lo,hi in ((-1e-3,1e-2),(0.0,1e-2),(-1e-3,2e-3)):
        es=np.linspace(lo,hi,int((hi-lo)/5e-6))
        x,S=sc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=(sc.Sz[0].get_dagger(),sc.Sz[0]),delta=2e-5,es=es)
        S=np.real(S)
        r=sc.get_kondo_spectrum(eV,order=2,U=0.,mode="DMRG",T=0,submode="KPM",delta=2e-5,es=es)[1]
        print("B",B,"es",lo,hi,"S[0]/Smax %.3g"%(S[0]/S.max()),"maxdiff vs ED %.3f"%np.abs(r-ed).max(),"eV=0 %.3f ED %.3f"%(r[20],ed[20]))
```

`01_es_edge.out`:

```
<repo>/src/dmrgpy/__init__.py
es -0.001 0.01 maxdiff 0.423 eV=0: 1.037 ED 1.037 nwarn 0 peak 4.712
es -0.0001 0.01 maxdiff 0.423 eV=0: 1.037 ED 1.037 nwarn 0 peak 4.712
es -2e-05 0.01 maxdiff 0.384 eV=0: 0.958 ED 1.037 nwarn 0 peak 4.712
es 0.0 0.01 maxdiff 1.037 eV=0: 0.0 ED 1.037 nwarn 0 peak 4.712
es -0.001 0.002 maxdiff 0.419 eV=0: 1.037 ED 1.037 nwarn 0 peak 4.712
es -0.001 0.003 maxdiff 0.422 eV=0: 1.037 ED 1.037 nwarn 0 peak 4.712
```

`02_edge_probe.out`:

```
B 10.0 es -0.001 0.01 S[0]/Smax 0 maxdiff vs ED 0.423 eV=0 1.037 ED 1.037
B 10.0 es 0.0 0.01 S[0]/Smax 1 maxdiff vs ED 1.037 eV=0 0.000 ED 1.037
B 10.0 es -0.001 0.002 S[0]/Smax 0 maxdiff vs ED 0.419 eV=0 1.037 ED 1.037
B 0.0 es -0.001 0.01 S[0]/Smax 0 maxdiff vs ED 0.026 eV=0 0.000 ED 0.000
B 0.0 es 0.0 0.01 S[0]/Smax 2.59e-07 maxdiff vs ED 0.025 eV=0 0.000 ED 0.000
B 0.0 es -0.001 0.002 S[0]/Smax 0 maxdiff vs ED 0.054 eV=0 0.000 ED 0.000
```

**Suggested fix (hunter)**: In second_order_dIdV_dc raise (or at least warn) when es.min() > -3*delta or es.max() < max|eV| + 3*delta, and in the potential term the same lower check; the checks use only es, delta and eVs, so they cost nothing. Returned numbers do not change for any grid that satisfies them.

**Fix (reviewer)**: In secondorder_dc._cumulative_theta0_weight, after sorting, use the data: warn (RuntimeWarning) when S[0] exceeds 1e-2 of S.max() (the lower edge cuts a line; the discriminator is clean, 1 against 2.6e-7), and warn when x.max() < max|eV| (only necessary condition for the top side; it can warn with no loss when no line lies above, which is acceptable). This needs eVs, so pass it in. Do not use the hunter es.min() > -3*delta rule alone. Add a sentence to docs/guide/kondo.md for order=2 saying es must start several delta below 0 and pass max|eV|. Pin with a test in tests/test_audit_2026_10_kondo.py style: the 4-site chain above at es from 0 warns, and es from -1e-3 and the B=0 singlet case do not. Values do not change.

### 13. An integer dt2 raises a UFuncTypeError after the first row of TDVP

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `correlator-twotime` &middot; candidate `correlator-twotime-4`

**Status**: FIXED in `4d1abbc` on 2026-11-12. `two_time_kondo_term_dmrg` casts `dt2` and `dtau` to float and validates the four grid arguments before any TDVP step; `kondo_term_from_two_time` casts too. No number changes. Pinned by `tests/test_audit_2026_11_twotime.py::test_integer_step_equals_float_step` and `::test_bad_grid_raises_before_any_tdvp_step`.

The claim as reviewed: An integer dt2 (dt2=1000 instead of 1000.0) in get_kondo_spectrum(order=3, mode="DMRG", submode="KPM") raises UFuncTypeError, but only after the first row of TDVP trajectories has been computed. The same call with 1000.0 returns a curve.

**Where**: src/dmrgpy/kondospectrumtk/dmrgtwotime.py:269 (t2_grid = dt2*np.arange(...) is an integer array for an int dt2) and twotime.py:221-223 (weights = np.full(len(t2_chunk), dt2) is then int64 and weights[...] *= 0.5 fails)

dt2 = t2_grid[1]-t2_grid[0] inherits the integer dtype, np.full(..., dt2) makes an int array and the in-place half weight for the end points cannot be cast. It is a loud failure, not a wrong number, but it fires after the whole first row (nine TDVP trajectories plus the t2 trajectories) has been paid for, and 1000 is the natural way to type the argument. The ED reference builds its grid with linspace and never hits it. Anchor: the same call with 1000.0 returns a curve.

The hunter's claim, before review: An integer dt2 (dt2=1000 instead of 1000.0) in get_kondo_spectrum(order=3, mode="DMRG", submode="KPM") raises UFuncTypeError, but only after the t2 trajectories and the first row of TDVP have been computed.

Decisive output, quoted by the hunter:

```
dt2 1000.0 ok [5.093 2.27  1.492 2.27  5.093]
dt2 1000 -> UFuncTypeError Cannot cast ufunc 'multiply' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`17_int_dt2.py`:

```python
import numpy as np, warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain
sc=spinchain.Spin_Chain(["1/2"]*3, itensor_version="python")
h=2*5.7883818066e-5*10*sc.Sz[0]+1e-3*(sc.Sx[0]*sc.Sx[1]+sc.Sy[0]*sc.Sy[1]+sc.Sz[0]*sc.Sz[1])+1e-3*(sc.Sx[1]*sc.Sx[2]+sc.Sy[1]*sc.Sy[2]+sc.Sz[1]*sc.Sz[2])
sc.set_hamiltonian(h)
eV=np.linspace(-2e-3,2e-3,5)
for dt2 in (1000.0, 1000):
    try:
        r=sc.get_kondo_spectrum(eV,Jrho_s=-0.1,order=3,mode="DMRG",T=0,submode="TD" if False else "KPM",delta=1e-4,es=np.linspace(-1e-3,1e-2,500),
              dt2=dt2,n_t2_half=2,dtau=1000.,n_tau_half=5,omega0=2e-3,Gamma0=2e-4)
        print("dt2",repr(dt2),"ok",np.round(r[1],3))
    except Exception as e: print("dt2",repr(dt2),"->",type(e).__name__,e)
```

Observed:

`17_int_dt2.out`:

```
dt2 1000.0 ok [5.093 2.27  1.492 2.27  5.093]
dt2 1000 -> UFuncTypeError Cannot cast ufunc 'multiply' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
```

**Reviewer (CONFIRMED)**: I reproduced it. It is a loud failure, not a wrong number, and nothing documents an integer dt2 as unsupported. The cause is as stated: weights = np.full(len(t2_chunk), dt2) is int64 when dt2 is int, and weights[...] *= 0.5 cannot be cast. The row is already paid for because the batch generator has run by then. The claim is real but minor, and I struck no sub-claim.

Reviewer's own evidence:

```
03_int.out: 'dt2 1000.0 ok [5.093 2.27 1.492 2.27 5.093]' and 'dt2 1000 -> UFuncTypeError Cannot cast ufunc multiply output from dtype(float64) to dtype(int64) with casting rule same_kind'.
```

Reviewer's scripts and outputs:

`03_int.py`:

```python
import numpy as np, warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain
sc=spinchain.Spin_Chain(["1/2"]*3, itensor_version="python")
h=2*5.7883818066e-5*10*sc.Sz[0]+1e-3*(sc.Sx[0]*sc.Sx[1]+sc.Sy[0]*sc.Sy[1]+sc.Sz[0]*sc.Sz[1])+1e-3*(sc.Sx[1]*sc.Sx[2]+sc.Sy[1]*sc.Sy[2]+sc.Sz[1]*sc.Sz[2])
sc.set_hamiltonian(h)
eV=np.linspace(-2e-3,2e-3,5)
for dt2 in (1000.0, 1000):
    try:
        r=sc.get_kondo_spectrum(eV,Jrho_s=-0.1,order=3,mode="DMRG",T=0,submode="TD" if False else "KPM",delta=1e-4,es=np.linspace(-1e-3,1e-2,500),
              dt2=dt2,n_t2_half=2,dtau=1000.,n_tau_half=5,omega0=2e-3,Gamma0=2e-4)
        print("dt2",repr(dt2),"ok",np.round(r[1],3))
    except Exception as e: print("dt2",repr(dt2),"->",type(e).__name__,e)
```

`03_int.out`:

```
dt2 1000.0 ok [5.093 2.27  1.492 2.27  5.093]
dt2 1000 -> UFuncTypeError Cannot cast ufunc 'multiply' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
```

**Suggested fix (hunter)**: Convert dt2 and dtau to float at the top of two_time_kondo_term_dmrg (and weights = np.full(len(t2_chunk), float(dt2)) in kondo_term_from_two_time), and validate the four grid arguments (positive, n_t2_half and n_tau_half positive integers) before any TDVP starts. No number changes.

**Fix (reviewer)**: The hunter's fix is right. Cast dt2 and dtau to float at the top of two_time_kondo_term_dmrg, and also use np.full(..., float(dt2)) in kondo_term_from_two_time (twotime.py), so the ED and any other caller are covered. Validate the four grid arguments (positive steps, positive integer half-widths) before the first TDVP step. Pin it with a test that calls with dt2=int on a 3-site python chain and compares against the float call, plus a ValueError test for a zero n_t2_half. No numbers change.

### 14. Sz is run twice and the F0 kernel is built on twice the rows: 4% of the call

`performance` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `correlator-twotime` &middot; candidate `correlator-twotime-5`

**Status**: NOT APPLIED. The redundancy is real, but it is 4% of the public call, and the hunter's form makes the second-order term depend on how the pair is spelled at 1e-8. The reviewer's cheaper form, a per-call memo of the Sz run and the F0 kernel built once by evenness, moves no number and is left here for the next pass. Re-measured on 2026-11-13 on HEAD, after finding 11 made the two-time term faster, with the reviewer.s `06_fraction.py` below (v3, one thread, this workstation): the public call takes 100.1 s, the pair 14.8 s, a fraction of 0.15, so the hunter.s 2.35x on the pair would now save 8% of the call, not 4%.

The claim as reviewed: With a correlator submode (not the default TASK) at order=3 and U!=0, Sz is computed twice and the second-order term and the potential term use three distinct spectral functions through six KPM runs, with the F0 kernel built on twice the rows needed; the redundancy is real but the dc pair is only 6% of the public call (15.7 s of 255.8 s, 8 spins, v3, KPM, delta=1e-5, 40001 es points, 21 biases, guide time grids) so sharing saves about 4% of it, not 2.35x of the call. The 5.2e-9 shift of the second-order term is KPM recursion noise of the MPS route, not a different quantity.

**Where**: src/dmrgpy/kondospectrumtk/secondorder_dc.py, the `weight()` closure of `second_order_dIdV_dc` (weights of S-, S+, Sz; lines 79-99 when recorded) and potentialdc.py:101-168 (_convolved_F0_weight per op, kernel built on [eVs,-eVs])

0.5*w(S-)+0.5*w(S+)+w(Sz) equals w(Sx)+w(Sy)+w(Sz) moment by moment, so the second-order term is the Theta0 cumulative of S_tot(w)=sum_k S_kk(w), and the potential term is the F0 convolution of the same S_tot. Currently each of the six calls to get_dynamical_correlator is a separate KPM run (Sz is run twice, with identical moments). F0 is even, so the kernel for -eV is minus the kernel for +eV, and h(eV)-h(-eV) collapses to weighted(eV)*(Theta0(eV)+Theta0(-eV)): half the rows, and the kernel is op-independent, so three builds become one. Measured on 8 sites, delta=1e-5, 40001 es points, 201 eV: second 5.2e-9 of the peak, potential 3.9e-15 (the 5.2e-9 is identical on v3 and python to six digits, so it is not truncation noise; it comes from the complex S+- names and is far below delta effects). Honest fraction: it only applies when the potential term is requested, and there the two-time term costs far more (56 s for a 7 row grid at 8 sites on v3, minutes to hours for a converged grid), so the saving is about 10 s of a call that is dominated by candidate 1's term; for order=2 nothing changes.

The hunter's claim, before review: With order=3 and U!=0 on a correlator submode, the six KPM correlators (Sx,Sy,Sz for the potential term; S-,S+,Sz for the second-order term) carry only three distinct spectral functions, and the F0 kernel is rebuilt three times on twice the rows needed; sharing them takes the second-order plus potential pair from 17.67 s to 7.52 s on v3 (2.35x) and from 21.25 s to 10.18 s on python (2.09x, d22), with the potential term identical to 4e-15 and the second-order term moved by 5.2e-9 of its peak.

Decisive output, quoted by the hunter:

```
v3: second max abs diff 3.64e-08 rel to peak 5.2e-09 | potential max abs diff 3.0e-15 rel to peak 3.9e-15 | old median s 17.67 | new median s 7.52 | speedup 2.349
python d22: second rel to peak 5.2e-09 | potential rel to peak 5.2e-15 | old median s 21.25 | new median s 10.18 | speedup 2.087
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`proto_kpm_share.py`:

```python
"""Prototype: one S_tot(w)=sum_k S_kk(w) from Sx,Sy,Sz, shared by the second-order
cumulative weight and the potential convolution; the F0 kernel on nev rows only
(F0 is even, so the -eV rows are minus the +eV rows)."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from dmrgpy.kondospectrumtk.stepfunctions import F0, Theta0

def s_tot(chain, site, mode, submode, delta, es, **kw):
    S = 0.; x0 = None
    for op in (chain.Sx[site], chain.Sy[site], chain.Sz[site]):
        x, s = chain.get_dynamical_correlator(mode=mode, submode=submode,
                 name=(op.get_dagger(), op), delta=delta, es=es, **kw)
        x = np.asarray(x, float); o = np.argsort(x, kind="stable")
        S = S + np.asarray(s)[o].real; x0 = x[o]
    return x0, S

def both_terms(chain, site, eVs, Jrho_s, U, T0, omega0, Gamma0, mode, submode, delta, es, **kw):
    eVs = np.asarray(eVs, float)
    x, S = s_tot(chain, site, mode, submode, delta, es, **kw)
    cum = cumulative_trapezoid(S, x, initial=0.)
    second = 2*np.pi*T0**2*(np.interp(eVs, x, cum) + np.interp(-eVs, x, cum) + 4*U**2)
    dx = np.diff(x); w = np.zeros_like(x); w[1:] += dx/2; w[:-1] += dx/2
    ws = w*S; weighted = np.empty(len(eVs))
    for a in range(0, len(eVs), 16):               # blocks, bounded memory
        e = eVs[a:a+16, None]
        K = F0((e - x[None, :]).ravel(), omega0=omega0, Gamma0=Gamma0).reshape(len(e), -1) \
          - F0((e + x[None, :]).ravel(), omega0=omega0, Gamma0=Gamma0).reshape(len(e), -1)
        weighted[a:a+16] = K @ ws
    # h(eV)-h(-eV) with h = weighted*Theta0 and weighted odd in eV
    pot = weighted*(Theta0(eVs) + Theta0(-eVs))
    return second, 4*np.pi*T0**2*Jrho_s*U*pot
```

`12_kpm_share_check.py`:

```python
import sys, time, warnings, numpy as np
warnings.filterwarnings("ignore"); sys.path.insert(0,".")
from dmrgpy import spinchain
import proto_kpm_share as P
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
G=2.0; MUB=5.7883818066e-5
v=sys.argv[1] if sys.argv[1]=="python" else int(sys.argv[1]); n=int(sys.argv[2]); delta=float(sys.argv[3]); nes=int(sys.argv[4]); nev=int(sys.argv[5])
sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version=v)
h=G*MUB*10.*sc.Sz[0]
for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12; sc.get_gs()
eV=np.linspace(-5e-3,5e-3,nev); es=np.linspace(-2e-2,2e-2,nes)
def old():
    a=second_order_dIdV_dc(sc,0,eV,mode="DMRG",submode="KPM",delta=delta,es=es,U=0.3)
    b=third_order_potential_dIdV_dc(sc,0,eV,-0.1,0.3,mode="DMRG",submode="KPM",delta=delta,es=es)
    return a,b
def new(): return P.both_terms(sc,0,eV,-0.1,0.3,1.0,20e-3,5e-6,"DMRG","KPM",delta,es)
res={}
for name,f in (("old",old),("new",new)):
    ts=[]
    for r in range(3):
        t=time.time(); out=f(); ts.append(time.time()-t)
    res[name]=(out,sorted(ts)[1],ts)
for i,nm in enumerate(("second","potential")):
    d=np.abs(res["old"][0][i]-res["new"][0][i]).max(); print(nm,"max abs diff",d,"rel to peak",d/np.abs(res["old"][0][i]).max())
for k in res: print(k,"median s",round(res[k][1],2),[round(x,2) for x in res[k][2]])
print("speedup",res["old"][1]/res["new"][1])
```

`11_prof_kpm.py`:

```python
import sys, time, cProfile, pstats, warnings, numpy as np
warnings.filterwarnings("ignore")
from dmrgpy import spinchain
G=2.0; MUB=5.7883818066e-5
v=sys.argv[1] if sys.argv[1]=="python" else int(sys.argv[1]); n=int(sys.argv[2]); delta=float(sys.argv[3]); nes=int(sys.argv[4]); nev=int(sys.argv[5])
sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version=v)
h=G*MUB*10.*sc.Sz[0]
for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12; sc.get_gs()
eV=np.linspace(-5e-3,5e-3,nev); es=np.linspace(-2e-2,2e-2,nes)
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
pr=cProfile.Profile(); t=time.time(); pr.enable()
a=second_order_dIdV_dc(sc,0,eV,mode="DMRG",submode="KPM",delta=delta,es=es,U=0.3)
b=third_order_potential_dIdV_dc(sc,0,eV,-0.1,0.3,mode="DMRG",submode="KPM",delta=delta,es=es)
pr.disable(); print("version",v,"n",n,"delta",delta,"nes",nes,"nev",nev,"total",round(time.time()-t,2))
pstats.Stats(pr).sort_stats("tottime").print_stats(10)
np.save("11_b_%s.npy"%sys.argv[1],b)
```

Observed:

`12_kpm_share_v3_n8.out`:

```
second max abs diff 3.640152623773929e-08 rel to peak 5.219354596231663e-09
potential max abs diff 2.9976021664879227e-15 rel to peak 3.906970207738117e-15
old median s 17.67 [18.21, 17.67, 17.34]
new median s 7.52 [7.52, 7.52, 7.54]
speedup 2.349005127642054
```

`19_kpm_share_python_n8_d22.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`11_prof_kpm_v3_n8.out`:

```
version 3 n 8 delta 1e-05 nes 40001 nev 201 total 17.44
         5686 function calls (5634 primitive calls) in 17.443 seconds

   Ordered by: internal time
   List reduced from 276 to 10 due to restriction <10>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        6   13.954    2.326   13.954    2.326 {built-in method dmrgpy.mpscpp3._dmrgcpp.kpm_dynamical_correlator}
        6    2.342    0.390    2.342    0.390 <repo>/src/dmrgpy/kondospectrumtk/stepfunctions.py:111(F0)
        6    0.720    0.120    0.734    0.122 <repo>/src/dmrgpy/algebra/kpm.py:840(generate_profile)
        3    0.320    0.107   10.151    3.384 <repo>/src/dmrgpy/kondospectrumtk/potentialdc.py:43(_convolved_F0_weight)
        3    0.064    0.021    0.064    0.021 {built-in method numpy._core._multiarray_umath.c_einsum}
        6    0.014    0.002    0.014    0.002 <repo>/src/dmrgpy/algebra/kpm.py:948(jackson_kernel)
        1    0.007    0.007   10.158   10.158 <repo>/src/dmrgpy/kondospectrumtk/potentialdc.py:101(third_order_potential_dIdV_dc)
       18    0.003    0.000    0.003    0.000 {built-in method numpy._core._multiarray_umath.interp}
        6    0.003    0.001    0.747    0.124 <repo>/src/dmrgpy/kpmdmrg.py:240(dynamical_correlator_from_moments)
       36    0.001    0.000    0.001    0.000 {built-in method numpy.array}
```

`18_prof_kpm_python_n8_d22.out`:

```
version python n 8 delta 1e-05 nes 40001 nev 201 total 41.78
         80068635 function calls (80068576 primitive calls) in 41.589 seconds

   Ordered by: internal time
   List reduced from 522 to 10 due to restriction <10>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
   131663    3.685    0.000    4.771    0.000 <d22 home>/venv/lib/python3.12/site-packages/numpy/linalg/_linalg.py:1667(svd)
   131663    2.168    0.000   17.951    0.000 <d22>/src/dmrgpy/pyitensor/svd.py:209(svd)
   493713    1.598    0.000    9.232    0.000 <d22>/src/dmrgpy/pyitensor/tensor.py:194(__mul__)
   135289    1.552    0.000    3.600    0.000 <d22>/src/dmrgpy/pyitensor/svd.py:33(_truncate)
   474958    1.356    0.000    2.139    0.000 <d22>/src/dmrgpy/pyitensor/tensor.py:123(_contract_matmul_impl)
  1138108    1.353    0.000    1.353    0.000 {method 'reduce' of 'numpy.ufunc' objects}
  1180387    1.348    0.000    2.666    0.000 <d22>/src/dmrgpy/pyitensor/tensor.py:147(__init__)
  2659133    1.283    0.000    2.401    0.000 <d22>/src/dmrgpy/pyitensor/tensor.py:25(_find)
        6    1.188    0.198    1.188    0.198 <d22>/src/dmrgpy/kondospectrumtk/stepfunctions.py:111(F0)
   474961    1.127    0.000    3.492    0.000 <d22>/src/dmrgpy/pyitensor/tensor.py:66(contract_arrays)
```

**Reviewer (CONFIRMED (narrowed))**: Why 5.2e-9: 0.5 w(S-)+0.5 w(S+) equals w(Sx)+w(Sy) exactly as bilinear moments. I checked on mode=ED (same Lorentzian route): the two spellings agree to 3.8e-16 of the peak. On v3 KPM at n=8, delta=1e-5, 40001 points they differ by 2.0e-8 of the S peak pointwise (cumulative 1.3e-8 of its total, 5.2e-9 of the dI/dV peak in the hunter run); at n=6, delta=2e-5 the difference is 1.2e-11 and the integrals agree to 12 digits. A repeated S- run reproduces itself exactly (0.0), so it is deterministic, but it depends on how the operator is spelled, which is the truncation/rounding of the Chebyshev recursion on a different vector, growing with the number of moments. It is far below the 0.4 delta floor against ED, so it is not an accuracy change; numbers move beyond roundoff, at 1e-8. The F0 evenness claim holds (F0(x)-F0(-x) max 0.0). Scope: the default submode on python and v3 is TASK, which never calls these functions, so this applies to explicit submode=KPM/CVM only; for order=2 nothing is saved (3 runs either way); the two-time Kondo term dominates the call (about 240 s of 256). Not re-verified: complex Hamiltonian, n_gs>1, spinful chains, python backend (I did not repeat the d22 timing).

Reviewer's own evidence:

```
06_fraction.out: 'second 7.6 s  potential 8.1 s  sum 15.7 / public order=3 U=0.3 call 255.8 s ; dc pair fraction 0.06 ; saving at 2.35x = 0.04 of call'. 05_ed_anchor.out: 'DMRG pm-xy rel 1.955e-08 / ED pm-xy rel 3.788e-16'. 04_why_5e9_n8.out: 'repeat Sm vs Sm: 0.000e+00 / cum diff max 6.575e-09 rel to cum max 1.315e-08'.
```

Reviewer's measurement: Own single run on v3, 8 spins: dc pair 15.7 s of a 255.8 s public call (6%); at the hunter 2.35x on the pair the saving is 4% of the call. Not repeated 3 times: the fraction, not the pair speedup, is what bounds it. Max relative difference of the second-order term 5.2e-9 of peak (hunter), 2e-8 of S pointwise (mine); potential term roundoff.

Reviewer's scripts and outputs:

`03_why_5e9.py`:

```python
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0,"<repo>/src")
from dmrgpy import spinchain
G=2.0; MUB=5.7883818066e-5
v=3 if sys.argv[1]=="3" else "python"
n=6
sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version=v)
h=G*MUB*10.*sc.Sz[0]
for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12; sc.get_gs()
es=np.linspace(-2e-2,2e-2,4001); d=2e-5
def S(op):
    x,s=sc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=(op.get_dagger(),op),delta=d,es=es)
    return np.asarray(x,float),np.asarray(s)
Sm=sc.Sx[0]-1j*sc.Sy[0]; Sp=sc.Sx[0]+1j*sc.Sy[0]
x,a=S(Sm); _,b=S(Sp); _,c=S(sc.Sx[0]); _,e=S(sc.Sy[0]); _,f=S(sc.Sz[0])
pm=0.5*a+0.5*b; xy=c+e
print("imag parts max: Sm %.2e Sp %.2e Sx %.2e Sy %.2e"%tuple(np.abs(q.imag).max() for q in (a,b,c,e)))
print("max |pm-xy| real: %.3e  rel to max S %.3e"%(np.abs(pm.real-xy.real).max(), np.abs(pm.real-xy.real).max()/np.abs(xy.real).max()))
print("integral pm %.12f xy %.12f"%(np.trapezoid(pm.real,x),np.trapezoid(xy.real,x)))
print("integral Sx %.12f Sy %.12f Sm %.12f Sp %.12f"%tuple(np.trapezoid(q.real,x) for q in (c,e,a,b)))
k=np.argmax(np.abs(pm.real-xy.real)); print("worst at w=%.5g, pm %.6g xy %.6g"%(x[k],pm.real[k],xy.real[k]))
```

`04_why_5e9_n8.py`:

```python
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0,"<repo>/src")
from dmrgpy import spinchain
from scipy.integrate import cumulative_trapezoid as ct
G=2.0; MUB=5.7883818066e-5
n=8
sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version=3)
h=G*MUB*10.*sc.Sz[0]
for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12; sc.get_gs()
es=np.linspace(-2e-2,2e-2,40001); d=1e-5
def S(op):
    x,s=sc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=(op.get_dagger(),op),delta=d,es=es)
    return np.asarray(x,float),np.asarray(s).real
Sm=sc.Sx[0]-1j*sc.Sy[0]; Sp=sc.Sx[0]+1j*sc.Sy[0]
x,a=S(Sm); _,b=S(Sp); _,c=S(sc.Sx[0]); _,e=S(sc.Sy[0]); _,a2=S(Sm)
print("repeat Sm vs Sm: %.3e"%np.abs(a-a2).max())
pm=0.5*a+0.5*b; xy=c+e
print("S diff max %.3e  rel %.3e (Smax %.3g)"%(np.abs(pm-xy).max(),np.abs(pm-xy).max()/xy.max(),xy.max()))
cp=ct(pm,x,initial=0.); cx=ct(xy,x,initial=0.)
print("cum diff max %.3e  rel to cum max %.3e"%(np.abs(cp-cx).max(),np.abs(cp-cx).max()/cx.max()))
print("totals pm %.10f xy %.10f"%(cp[-1],cx[-1]))
print("min of S: pm %.3e xy %.3e"%(pm.min(),xy.min()))
```

`05_ed_anchor.py`:

```python
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0,"<repo>/src")
from dmrgpy import spinchain
G=2.0; MUB=5.7883818066e-5
n=8
sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version=3)
h=G*MUB*10.*sc.Sz[0]
for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12; sc.get_gs()
es=np.linspace(-2e-2,2e-2,40001); d=1e-5
def S(op,mode):
    x,s=sc.get_dynamical_correlator(mode=mode,submode="KPM" if mode=="DMRG" else "ED",name=(op.get_dagger(),op),delta=d,es=es)
    return np.asarray(s).real
Sm=sc.Sx[0]-1j*sc.Sy[0]; Sp=sc.Sx[0]+1j*sc.Sy[0]
for mode in ("DMRG","ED"):
    pm=0.5*S(Sm,mode)+0.5*S(Sp,mode); xy=S(sc.Sx[0],mode)+S(sc.Sy[0],mode)
    print(mode,"pm-xy rel %.3e"%(np.abs(pm-xy).max()/xy.max()))
    if mode=="DMRG": P,X=pm,xy
    else: print("DMRG(pm) vs ED-lorentz not comparable shape")
```

`06_fraction.py`:

```python
import sys, time, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0,"<repo>/src")
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
G=2.0; MUB=5.7883818066e-5
n=8
sc=spinchain.Spin_Chain(["1/2"]*n, itensor_version=3)
h=G*MUB*10.*sc.Sz[0]
for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12; sc.get_gs()
eV=np.linspace(-4e-3,4e-3,21); es=np.linspace(-2e-2,2e-2,40001); d=1e-5
t=time.time(); second_order_dIdV_dc(sc,0,eV,mode="DMRG",submode="KPM",delta=d,es=es,U=0.3); t2=time.time()-t
t=time.time(); third_order_potential_dIdV_dc(sc,0,eV,-0.1,0.3,omega0=2e-3,Gamma0=2e-3,mode="DMRG",submode="KPM",delta=d,es=es); tp=time.time()-t
print("second %.1f s  potential %.1f s  sum %.1f"%(t2,tp,t2+tp),flush=True)
t=time.time()
sc.get_kondo_spectrum(eV,Jrho_s=-0.1,U=0.3,T=0.,order=3,omega0=2e-3,Gamma0=2e-3,mode="DMRG",submode="KPM",delta=d,es=es,dt2=150.,n_t2_half=20,dtau=150.,n_tau_half=20)
tt=time.time()-t
print("public order=3 U=0.3 call %.1f s ; dc pair fraction %.2f ; saving at 2.35x = %.2f of call"%(tt,(t2+tp)/tt,(t2+tp)/tt*(1-1/2.35)))
```

`03_why_5e9.out`:

```
imag parts max: Sm 0.00e+00 Sp 0.00e+00 Sx 0.00e+00 Sy 0.00e+00
max |pm-xy| real: 7.456e-08  rel to max S 1.182e-11
integral pm 0.499999989223 xy 0.499999989223
integral Sx 0.249999994612 Sy 0.249999994612 Sm 0.084780990638 Sp 0.915218987808
worst at w=0.00194, pm 25.6647 xy 25.6647
imag parts max: Sm 2.47e-12 Sp 9.29e-12 Sx 4.48e-12 Sy 1.26e-12
max |pm-xy| real: 7.457e-08  rel to max S 1.182e-11
integral pm 0.499999989223 xy 0.499999989223
integral Sx 0.249999994612 Sy 0.249999994612 Sm 0.084780990638 Sp 0.915218987808
worst at w=0.00194, pm 25.6647 xy 25.6647
```

`04_why_5e9_n8.out`:

```
repeat Sm vs Sm: 0.000e+00
S diff max 2.211e-04  rel 1.955e-08 (Smax 1.13e+04)
cum diff max 6.575e-09  rel to cum max 1.315e-08
totals pm 0.5000000683 xy 0.5000000683
min of S: pm -5.365e-07 xy -2.152e-06
```

`05_ed_anchor.out`:

```
DMRG pm-xy rel 1.955e-08
ED pm-xy rel 3.788e-16
DMRG(pm) vs ED-lorentz not comparable shape
```

`06_fraction.out`:

```
second 7.6 s  potential 8.1 s  sum 15.7
public order=3 U=0.3 call 255.8 s ; dc pair fraction 0.06 ; saving at 2.35x = 0.04 of call
```

**Suggested fix (hunter)**: Compute S_tot once from the Sx, Sy, Sz correlators in get_kondo_spectrum's DMRG branch (or a small cache keyed on the chain state and the correlator arguments) and hand it to both term functions; in the potential term build the F0 kernel once, on the nev positive biases only, and use that h(eV)-h(-eV) reduces to weighted(eV). Second-order values move by 5e-9 of the peak, potential by roundoff.

**Fix (reviewer)**: The hunter fix (share S_tot from Sx,Sy,Sz and halve the F0 rows) is correct but not worth the 4%, and it makes the second-order term spelling-dependent at 1e-8 and breaks the symmetry with order=2, which keeps S-,S+,Sz. A cheaper fix that moves no number: a per-call memo in terms() keyed on (operator name, delta, es id) so Sz is run once (5 runs of 6), and build the F0 kernel once by evenness on the nev positive biases in potentialdc (exactly equal, F0 is even to 0.0); pin with an old-against-new equality test at 1e-14 on a 4-site chain. Apply it only if the author wants the 4%; recording it as low-value is reasonable.

### 15. Every T>0 call rebuilds the 0.2 to 0.6 s FBuilder table

`performance` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `dispatch` &middot; candidate `dispatch-1`

**Status**: FIXED on 2026-11-12 in two steps: `stepfunctions.cached_fbuilder`, a `functools.lru_cache` of four FBuilder objects keyed on float arguments, in `f6ac1e3`, and its three call sites (the ED route, the T>0 TASK route, the fallback in `conductance._get_F`) in `317c39f`. The object is shared, and its docstring says it must not be mutated. No number changes. Pinned by `tests/test_audit_2026_11_stepfunctions.py::test_cached_fbuilder_is_shared_and_keyed` and `::test_get_kondo_spectrum_builds_the_thermal_table_once` (two equal ED calls at 1 K build one table, a different omega0 a second one).

The claim as reviewed: Every T>0 get_kondo_spectrum call with order=3 rebuilds the FBuilder (0.63 s locally on a loaded host, 0.21 s on d22; 2399 spline nodes, about 100 KB), a pure function of (T, omega0, Gamma0, kB); this matters only for repeated calls at fixed T on a cheap or cached spectrum (12.3% of a repeat TASK call on 24 sites by the hunter's d22 measurement), and is 0.9% of a first call. On the ED route it is negligible against the per-call eigh.

**Where**: src/dmrgpy/kondospectrumtk/kondomethod.py:270 (_get_kondo_spectrum_ed), :319-322 (_get_kondo_spectrum_dmrg, TASK branch); src/dmrgpy/kondospectrumtk/stepfunctions.py:191 (FBuilder.__init__)

The FBuilder is built once per get_kondo_spectrum call, shared across the two third-order terms and across n_gs members, but not across calls. A field, Jrho_s, U or bias sweep on one chain at fixed T, omega0, Gamma0 and kB repeats the same 4001-point convolution on a 1601+800 point grid each time. FBuilder holds only scalars and a CubicSpline that are read-only after __init__, so a module-level memo keyed on the float values of its arguments is content-keyed and does not touch the cache rule of kondo.md (that rule is about the TASK spectrum, whose key stays the digest of the exported arrays). Prototype in 06b_fbuilder_cache.py: a wrapper that stores the instances in a dict, patched over stepfunctions.FBuilder, which both call sites import at call time. Old against new on identical inputs (24 sites, python, T=0.1 K, TASK spectrum cached): the returned dI/dV is bitwise identical (max relative difference 0.0). Nine alternating repeats on d22, one thread: old median 1.823 s (range 1.78 to 1.87), new median 1.598 s (range 1.57 to 1.79, the first of the new repeats includes the one build); the bare FBuilder construction is 0.214 s in 7 repeats, matching the saving of 0.225 s. That is 12.3% of a cached repeat call and 0.87% of the 25.8 s first call; on a 10-site chain locally it is 0.42 of 1.04 s. Worth doing only for sweeps over cheap chains or for repeated calls. The first draft of this measurement (06_fbuilder_cache.py, noisy host) is superseded by 06b and should not be cited.

The hunter's claim, before review: Every T>0 call with order=3 rebuilds the thermal Kondo-log table FBuilder (0.214 s on d22, 0.42 s locally), a pure function of (T, omega0, Gamma0, kB), which is 12.3% of a repeat call on 24 sites but only 0.9% of a first call.

Decisive output, quoted by the hunter:

```
bare FBuilder construction, 7 repeats: [0.213 0.214 0.213 0.213 0.214 0.214 0.214] median 0.214
old repeats [1.786 1.832 1.823 1.824 1.873 1.795 1.818 1.84  1.782] median 1.823
new repeats [1.792 1.604 1.603 1.577 1.622 1.598 1.587 1.591 1.567] median 1.598
max |new-old|/max|old| = 0.000e+00
median saving 0.225 s = 12.3% of a cached repeat call, 0.87% of the first call
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`06b_fbuilder_cache.py`:

```python
import time, warnings, numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import stepfunctions as sf
n=24
sc=spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
h=0
for i in range(n-1):
    J=0.6e-3 if i==0 else 1e-3
    h=h+J*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h=h+3e-4*sum(sc.Sz[i] for i in range(n))
sc.set_hamiltonian(h); sc.maxm,sc.nsweeps=40,16
eV=np.linspace(-6e-3,6e-3,121)
kw=dict(site=0,Jrho_s=-0.1,U=0.3,omega0=20e-3,Gamma0=5e-5,order=3,mode="DMRG",submode="TASK",nkry=100,T=0.1)
Orig=sf.FBuilder; _cache={}
def Cached(T, omega0=20e-3, Gamma0=5e-6, kB=8.617333262e-5, npts=4001, width=60., grid_span=10., tabulate=True):
    key=(float(T),float(omega0),float(Gamma0),float(kB),npts,width,grid_span,tabulate)
    if key not in _cache: _cache[key]=Orig(T,omega0=omega0,Gamma0=Gamma0,kB=kB,npts=npts,width=width,grid_span=grid_span,tabulate=tabulate)
    return _cache[key]
def call():
    t=time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); r=sc.get_kondo_spectrum(eV,**kw)[1]
    return time.time()-t, r
t=time.time(); sc.gs_energy(); print("DMRG solve %.1f s"%(time.time()-t))
t0,_=call(); print("first call (builds the TASK spectrum, T>0) %.2f s"%t0,flush=True)
fb=[]
for k in range(7):
    t=time.time(); Orig(0.1,omega0=20e-3,Gamma0=5e-5); fb.append(time.time()-t)
print("bare FBuilder construction, 7 repeats:",np.round(fb,3),"median %.3f"%np.median(fb))
old=[];new=[]
for k in range(9):
    sf.FBuilder=Orig; dt,ro=call(); old.append(dt)
    sf.FBuilder=Cached; dt,rn=call(); new.append(dt)
print("old repeats",np.round(old,3),"median %.3f"%np.median(old))
print("new repeats",np.round(new,3),"median %.3f"%np.median(new))
print("max |new-old|/max|old| = %.3e"%(np.max(np.abs(rn-ro))/np.max(np.abs(ro))))
print("median saving %.3f s = %.1f%% of a cached repeat call, %.2f%% of the first call"%(np.median(old)-np.median(new),100*(1-np.median(new)/np.median(old)),100*(np.median(old)-np.median(new))/t0))
```

`01_profile_small.py`:

```python
import cProfile, pstats, time, numpy as np, warnings
from dmrgpy import spinchain
KB=8.617333262e-5
def chain(n=10, maxm=30, version="python", seed=1):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    h = 0
    for i in range(n-1):
        J = (0.6e-3 if i==0 else 1e-3)
        h = h + J*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h = h + 4e-4*sum(sc.Sz[i] for i in range(n))
    sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = maxm, 20
    return sc
sc = chain()
eV = np.linspace(-6e-3,6e-3,61)
kw = dict(site=0, Jrho_s=-0.1, U=0.3, omega0=20e-3, Gamma0=5e-5, order=3, mode="DMRG", submode="TASK", nkry=60)
sc.gs_energy()
for T in (0.,0.2):
    for rep in range(2):
        t=time.time()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.get_kondo_spectrum(eV, T=T, **kw)
        print("T",T,"rep",rep,"%.2f s"%(time.time()-t))
pr = cProfile.Profile(); pr.enable()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sc.get_kondo_spectrum(eV, T=0.2, **kw)
pr.disable()
pstats.Stats(pr).sort_stats("cumulative").print_stats(18)
```

Observed:

`06b_fbuilder_cache.out`: (empty file: the d22 runner of this hunt copied the run folder back after each run and replaced the `.out` that `tee` was writing with the empty remote copy; the decisive lines are the agent's verbatim quote above)

`01_profile_small.out`:

```
T 0.0 rep 0 2.27 s
T 0.0 rep 1 0.32 s
T 0.2 rep 0 1.62 s
T 0.2 rep 1 1.39 s
         52890 function calls in 1.040 seconds

   Ordered by: cumulative time
   List reduced from 391 to 18 due to restriction <18>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    1.040    1.040 <repo>/src/dmrgpy/kondospectrumtk/kondomethod.py:23(get_kondo_spectrum)
        1    0.000    0.000    1.040    1.040 <repo>/src/dmrgpy/kondospectrumtk/kondomethod.py:278(_get_kondo_spectrum_dmrg)
        1    0.001    0.001    0.618    0.618 <repo>/src/dmrgpy/kondospectrumtk/kondomethod.py:323(terms)
        1    0.000    0.000    0.422    0.422 <repo>/src/dmrgpy/kondospectrumtk/stepfunctions.py:191(__init__)
        1    0.001    0.001    0.420    0.420 <repo>/src/dmrgpy/kondospectrumtk/stepfunctions.py:211(_direct)
        1    0.211    0.211    0.419    0.419 <repo>/src/dmrgpy/kondospectrumtk/stepfunctions.py:125(_F_convolve)
        1    0.003    0.003    0.326    0.326 <repo>/src/dmrgpy/kondospectrumtk/conductance.py:206(third_order_kondo_dIdV)
        2    0.008    0.004    0.319    0.159 <repo>/src/dmrgpy/kondospectrumtk/conductance.py:265(one_direction)
        8    0.024    0.003    0.307    0.038 <repo>/src/dmrgpy/kondospectrumtk/stepfunctions.py:214(__call__)
        1    0.010    0.010    0.236    0.236 <repo>/src/dmrgpy/kondospectrumtk/conductance.py:139(second_order_dIdV)
        6    0.012    0.002    0.224    0.037 <repo>/src/dmrgpy/kondospectrumtk/conductance.py:114(_theta_raw)
        6    0.212    0.035    0.212    0.035 <repo>/src/dmrgpy/kondospectrumtk/stepfunctions.py:58(Theta)
        8    0.000    0.000    0.185    0.023 <home>
        8    0.000    0.000    0.185    0.023 <home>
        8    0.185    0.023    0.185    0.023 <home>
        1    0.114    0.114    0.114    0.114 <repo>/src/dmrgpy/kondospectrumtk/stepfunctions.py:138(L)
        8    0.097    0.012    0.098    0.012 <repo>/src/dmrgpy/kondospectrumtk/stepfunctions.py:111(F0)
        1    0.000    0.000    0.093    0.093 <home>
```

**Reviewer (CONFIRMED (narrowed))**: I checked what I could independently: the build time is independent of T (0.632, 0.633, 0.634 s at 0.1, 1, 5 K), the object holds two scalars, a kT and a CubicSpline that nothing mutates after __init__, so a content-keyed memo is safe and does not touch kondo.md's cache rule (that rule concerns the TASK spectrum under a digest of the exported arrays; an FBuilder has no chain state). Call sites: kondomethod.py:270 and :322 build one per call, conductance.py:127 builds when Fb is None, and tests construct their own (test_kondo_spectrum_paper_fig7.py, test_kondo_spectrum.py) with tabulate default, so a factory leaves them unchanged. I did not re-run the 24-site TASK timing (d22, TASK solve cost), so the 12.3% figure rests on the hunter's output; it is consistent with the construction time I measured, and the saving per call equals the build time by construction, with identical output because the object is the same. Narrowing: the saving is a fixed 0.2 to 0.6 s per call, relevant only to sweeps over cached TASK spectra; for ED and for first calls it is noise.

Reviewer's own evidence:

```
T=0.1 0.6316 ; T=1.0 0.6330 ; T=5.0 0.6339 (s, FBuilder construction median of 5)
tracemalloc peak 384 MB during construction, retained 99480 bytes, 2399 nodes
```

Reviewer's measurement: construction median 0.63 s locally (loaded host), independent of T; hunter's d22 whole-call 1.823 s to 1.598 s (12.3%) not re-run by me; output bitwise identical by construction (same object). Numbers do not change.

Reviewer's scripts and outputs:

`05_fb.py`:

```python
import sys; sys.path.insert(0,"<repo>/src")
import time, numpy as np
from dmrgpy.kondospectrumtk import stepfunctions as sf
for T in (0.1,1.,5.):
    t=[]
    for k in range(5):
        s=time.perf_counter(); sf.FBuilder(T,omega0=20e-3,Gamma0=5e-5); t.append(time.perf_counter()-s)
    print(T,np.median(t))
import tracemalloc; tracemalloc.start(); fb=sf.FBuilder(1.); print(tracemalloc.get_traced_memory(), fb._F_spline.x.size)
```

`05_fb.out`:

```
0.1 0.6316174009989481
1.0 0.6329885119921528
5.0 0.6339259830128867
(99480, 384081839) 2399
```

**Suggested fix (hunter)**: Memoize FBuilder behind a small lru_cache (maxsize about 4) keyed on (float(T), float(omega0), float(Gamma0), float(kB), npts, width, grid_span, tabulate), exposed as a function in stepfunctions.py that both call sites use; keep the class for tabulate=False tests. Does not change returned numbers (bitwise equal in the measurement).

**Fix (reviewer)**: The hunter's fix is right in shape (a small lru_cache factory keyed on float(T), omega0, Gamma0, kB, npts, width, grid_span, tabulate, maxsize about 4, used at kondomethod.py:270, :322 and conductance.py:127). Note the construction peaks at 384 MB transiently, so keep maxsize small only for the retained object, which is 100 KB, and that is fine. Better: do ed-reference-2 first, since it removes a cost from every call, and add the memo only if sweeps over cached TASK spectra matter. Pin with a test that two calls with the same arguments return the same object and that a changed omega0 or kB returns a different one.

### 16. Three texts say julia_live gets the correlator route by default; the code sends it to TASK

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `dispatch` &middot; candidate `dispatch-3`

**Status**: FIXED in `76313ee` on 2026-11-12. The comment, the two docstring places and the implicit-default TypeError say TASK on `"python"`, 3 and `"julia_live"` and the correlator route (order=2 only) on v2; the TypeError builds its list from `taskdmrg.BACKENDS`. The sentence in `.claude/rules/kondo.md` was corrected with the integration. Pinned by `tests/test_audit_2026_11_dispatch.py::test_implicit_default_error_names_every_task_backend`.

The claim as reviewed: The comment (kondomethod.py ~288-290), the get_kondo_spectrum docstring (lines ~106-108 and ~185) and the implicit-default TypeError text (~314-315) say TASK is the default on python and 3 only and that julia_live gets the correlator route, while the code (BACKENDS = ("python",3,"julia_live"), taskdmrg.py:19) sends julia_live to TASK; docs/guide/kondo.md and ROADMAP.md already say all three, so the stale texts are the three in kondomethod.py plus the phrase 'the default on both' in .claude/rules/kondo.md. Text only, no number wrong.

**Where**: src/dmrgpy/kondospectrumtk/kondomethod.py:106-108, 185, 287-292, 314-315; src/dmrgpy/taskdmrg.py:19

Evidence by reading, julia_live not run. kondomethod.py:291-292 sets submode = "TASK" if self.itensor_version in BACKENDS else "KPM", and taskdmrg.py:19 is BACKENDS = ("python",3,"julia_live") (printed in 02_dispatch_checks.out). The comment at 288-290 says "v2 and julia_live have no TASK, so they get the correlator route", the docstring says "default TASK on itensor_version=python and 3 ... KPM on the others" and the TypeError at 314 says "TASK is the default on itensor_version=\"python\" and 3". The new T>0 message at 355-358 does list julia_live among the TASK backends, which is what the code does, so the file contradicts itself. Consequence for a julia_live user: a call without submode= takes the TASK route (no es/delta/grid, raises on them, T>0 allowed with the guard) rather than the KPM route the docstring describes, and a user who passes delta= or es= with the default gets a TypeError that blames a default the text says does not apply to them. No number is wrong.

The hunter's claim, before review: The default submode comment, the docstring and the TypeError text of get_kondo_spectrum say TASK is the default on python and 3 and that julia_live gets the correlator route, but the code sends julia_live to TASK because taskdmrg.BACKENDS contains it.

Decisive output, quoted by the hunter:

```
BACKENDS ('python', 3, 'julia_live')
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`02_dispatch_checks.py`:

```python
import numpy as np, warnings
from dmrgpy import spinchain, taskdmrg
print("BACKENDS", taskdmrg.BACKENDS)
sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
h=0
for i in range(3):
    h = h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h+4e-4*sum(sc.Sz[i] for i in range(4))
sc.set_hamiltonian(h); sc.maxm, sc.nsweeps, sc.cutoff = 16, 20, 0.
eV=np.linspace(-4e-3,4e-3,41)
kw=dict(site=0,Jrho_s=-0.1,U=0.3,omega0=20e-3,Gamma0=5e-5,order=3)
def run(label, f):
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r=f()
        print(label,"->", np.round(r[1][[0,20,40]],5) if hasattr(r,"__len__") else r, [str(x.message)[:60] for x in w])
    except Exception as e:
        print(label,"-> RAISES",type(e).__name__,str(e)[:110])
# kB handling at T>0 on ED
run("ED T=2 kB default", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=2.,**kw))
run("ED T=2 kB=-8.6e-5", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=2.,kB=-8.617333262e-5,**kw))
run("ED T=2 kB=0", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=2.,kB=0.,**kw))
run("ED T=nan", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=float("nan"),**kw))
run("ED T=-1", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=-1.,**kw))
run("TASK T=0.2 kB=-8.6e-5", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.2,kB=-8.617333262e-5,nkry=30,**kw))
run("TASK T=0.2 kB=0", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.2,kB=0.,nkry=30,**kw))
run("TASK T=0 kB=-1 (inert)", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,kB=-1.,nkry=30,**kw))
# mode typo with LL
run("mode=ed order=LL", lambda: sc.get_kondo_spectrum(eV,mode="ed",T=0.,**dict(kw,order="LL")))
run("mode=ed order=3", lambda: sc.get_kondo_spectrum(eV,mode="ed",T=0.,**kw))
# KPM with nkry
run("KPM nkry", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,submode="KPM",nkry=30,delta=1e-4,es=np.linspace(-.03,.03,100),**dict(kw,order=2)))
# order=2 inert args
a=sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,nkry=30,**dict(kw,order=2))[1]
b=sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,nkry=30,**dict(kw,order=2,Jrho_s=5.,omega0=1.,Gamma0=1e-2))[1]
print("order=2 insensitive to Jrho_s/omega0/Gamma0:", np.max(np.abs(a-b)))
# bool order / order=True
run("order=True", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,order=True)))
run("order=3.0", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,order=3.0)))
run("order=np.int64(3)", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,order=np.int64(3))))
run("site=-1 ED", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,site=-1)))
run("site=-1 TASK", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,nkry=30,**dict(kw,site=-1)))
run("site=4", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,site=4)))
run("T=1 string", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T="1",nkry=30,**kw))
run("empty eV TASK T>0", lambda: sc.get_kondo_spectrum(np.array([]),mode="DMRG",T=0.1,nkry=30,**kw))
```

Observed:

`02_dispatch_checks.out`:

```
BACKENDS ('python', 3, 'julia_live')
ED T=2 kB default -> [8.87539 4.42909 8.67066] []
ED T=2 kB=-8.6e-5 -> RAISES ValueError `x` must contain only finite values.
ED T=2 kB=0 -> RAISES ValueError Geometric sequence cannot include zero
ED T=nan -> RAISES ValueError `x` must contain only finite values.
ED T=-1 -> RAISES ValueError KondoSpectrum requires T>=0
TASK T=0.2 kB=-8.6e-5 -> RAISES ValueError `x` must contain only finite values.
TASK T=0.2 kB=0 -> RAISES ValueError Geometric sequence cannot include zero
TASK T=0 kB=-1 (inert) -> [8.88649 2.26195 8.66271] []
mode=ed order=LL -> RAISES NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kondo_leading_log/README.md
mode=ed order=3 -> RAISES ValueError mode must be "ED" or "DMRG"
KPM nkry -> RAISES TypeError Unexpected keyword argument(s) for the KPM dynamical correlator: nkry. Solver parameters are attributes of the
order=2 insensitive to Jrho_s/omega0/Gamma0: 0.0
order=True -> RAISES ValueError order must be 2, 3 or "LL"
order=3.0 -> [8.88649 2.26195 8.66271] []
order=np.int64(3) -> [8.88649 2.26195 8.66271] []
site=-1 ED -> [8.88649 2.26195 8.66271] []
site=-1 TASK -> [8.88649 2.26195 8.66271] []
site=4 -> RAISES IndexError list index out of range
T=1 string -> RAISES TypeError '>=' not supported between instances of 'str' and 'float'
empty eV TASK T>0 -> RAISES IndexError index 0 is out of bounds for axis 0 with size 0
```

**Reviewer (CONFIRMED (narrowed))**: Checked by reading (julia_live not run, per brief). The code line `submode = "TASK" if self.itensor_version in BACKENDS else "KPM"` and BACKENDS contain julia_live. Not intended or documented as otherwise: ROADMAP.md row for the DMRG third-order route says the default submode is TASK on all three, and docs/guide/kondo.md says TASK is the default on 3, python and julia_live and KPM stays on v2. So the code and the user-facing docs agree and only the in-code texts and the rule file are stale. The hunter's consequence (a julia_live user passing delta=/es= gets a TypeError that blames a default the text says does not apply) is real but only cosmetic: the error still names submode="KPM" as the fix. Narrowing: the guide is already correct, so the fix is not in the guide, and the T>0 message already lists julia_live.

Reviewer's own evidence:

```
kondomethod.py: `submode = "TASK" if self.itensor_version in BACKENDS else "KPM"`; taskdmrg.py:19 `BACKENDS = ("python",3,"julia_live")`; comment 'v2 and julia_live have no TASK, so they get the correlator route'; docs/guide/kondo.md: 'the default on itensor_version=3, "python" and "julia_live"'.
```

Reviewer's scripts and outputs:


`02_dispatch_checks.out`:

```
BACKENDS ('python', 3, 'julia_live')
ED T=2 kB default -> [8.87539 4.42909 8.67066] []
ED T=2 kB=-8.6e-5 -> RAISES ValueError `x` must contain only finite values.
ED T=2 kB=0 -> RAISES ValueError Geometric sequence cannot include zero
ED T=nan -> RAISES ValueError `x` must contain only finite values.
ED T=-1 -> RAISES ValueError KondoSpectrum requires T>=0
TASK T=0.2 kB=-8.6e-5 -> RAISES ValueError `x` must contain only finite values.
TASK T=0.2 kB=0 -> RAISES ValueError Geometric sequence cannot include zero
TASK T=0 kB=-1 (inert) -> [8.88649 2.26195 8.66271] []
mode=ed order=LL -> RAISES NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kondo_leading_log/README.md
mode=ed order=3 -> RAISES ValueError mode must be "ED" or "DMRG"
KPM nkry -> RAISES TypeError Unexpected keyword argument(s) for the KPM dynamical correlator: nkry. Solver parameters are attributes of the
order=2 insensitive to Jrho_s/omega0/Gamma0: 0.0
order=True -> RAISES ValueError order must be 2, 3 or "LL"
order=3.0 -> [8.88649 2.26195 8.66271] []
order=np.int64(3) -> [8.88649 2.26195 8.66271] []
site=-1 ED -> [8.88649 2.26195 8.66271] []
site=-1 TASK -> [8.88649 2.26195 8.66271] []
site=4 -> RAISES IndexError list index out of range
T=1 string -> RAISES TypeError '>=' not supported between instances of 'str' and 'float'
empty eV TASK T>0 -> RAISES IndexError index 0 is out of bounds for axis 0 with size 0
```

**Suggested fix (hunter)**: Rewrite the three texts (and the kondo.md sentence "the default on both") to say TASK is the default on python, 3 and julia_live, and KPM on v2 only. Changes no numbers.

**Fix (reviewer)**: Hunter fix is right. Edit the comment, the docstring lines ~106-108 and ~185 and the TypeError string in kondomethod.py to say TASK is the default on python, 3 and julia_live and KPM on v2 only (and that on v2 the correlator route covers order=2 only), and the 'the default on both' sentence in .claude/rules/kondo.md. Optionally build the TypeError text from taskdmrg.BACKENDS so it cannot drift. A test asserting the message contains julia_live (or iterates BACKENDS) would pin it.

### 17. mode="ed" with order="LL" is told LL is ED-only; n_gs is checked after setup

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `dispatch` &middot; candidate `dispatch-4`

**Status**: FIXED in `76313ee` on 2026-11-12. `mode` is validated first, then the LL checks; `n_gs` is validated at the top of `_get_kondo_spectrum_dmrg`, before the submode, the T>0 FBuilder and any solve. The commit was written on a tree without the T>0 route; the integration merged it onto `b0e923d` and extended its test to T=0.5. Pinned by `tests/test_audit_2026_11_dispatch.py::test_mistyped_mode_with_LL_is_a_value_error_naming_the_modes`, `::test_LL_on_DMRG_is_still_not_implemented` and `::test_n_gs_is_checked_before_anything_is_built` (FBuilder, TaskKondoSpectrum, KondoSpectrum and `get_gs` made to fail, at T=0 and T=0.5).

The claim as reviewed: In get_kondo_spectrum, the order="LL" and mode != "ED" refusal runs before mode is validated, so mode="ed" (a typo) with order="LL" raises NotImplementedError saying LL is ED-only instead of the ValueError 'mode must be "ED" or "DMRG"' that order=3 gets; and n_gs is validated at the end of _get_kondo_spectrum_dmrg, after the 0.4 s FBuilder is built, but only on the T>0 order=3 TASK route (and after the other DMRG setup), so a bad n_gs at T=0 wastes nothing measurable. Both raise loudly and no number is wrong.

**Where**: src/dmrgpy/kondospectrumtk/kondomethod.py:214-221 (LL mode check before the mode check at 249-255), :318-322 against :402-406 (FBuilder before n_gs validation)

get_kondo_spectrum(order="LL", mode="ed") raises NotImplementedError 'order="LL" is implemented on mode="ED" only', because `order == "LL" and mode != "ED"` runs before the `else: raise ValueError('mode must be "ED" or "DMRG"')`; with order=3 the same typo gets the correct ValueError (02_dispatch_checks.out). This is the 4.10 pattern (one raise serving unsupported and mistyped). Separately, n_gs=0 raises ValueError only after FBuilder(T, ...) has been tabulated (0.40 s locally, 12_cheap_checks.out: 'n_gs=0 -> ValueError after 0.40 s'); no number is wrong. Both are loud.

The hunter's claim, before review: Two validations run after the branch they qualify: n_gs is checked after the 0.4 s FBuilder is built, and order="LL" is refused before mode is validated, so mode="ed" with order="LL" is told LL is implemented on mode="ED" only.

Decisive output, quoted by the hunter:

```
mode=ed order=LL -> RAISES NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kondo_leading_log/README.md
mode=ed order=3 -> RAISES ValueError mode must be "ED" or "DMRG"
n_gs=0 -> ValueError after 0.40 s
```

Repro, hunter (every script ran through `run3.sh` locally or `rund22.sh` on d22, see "Shared helpers"):

`02_dispatch_checks.py`:

```python
import numpy as np, warnings
from dmrgpy import spinchain, taskdmrg
print("BACKENDS", taskdmrg.BACKENDS)
sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
h=0
for i in range(3):
    h = h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
h = h+4e-4*sum(sc.Sz[i] for i in range(4))
sc.set_hamiltonian(h); sc.maxm, sc.nsweeps, sc.cutoff = 16, 20, 0.
eV=np.linspace(-4e-3,4e-3,41)
kw=dict(site=0,Jrho_s=-0.1,U=0.3,omega0=20e-3,Gamma0=5e-5,order=3)
def run(label, f):
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r=f()
        print(label,"->", np.round(r[1][[0,20,40]],5) if hasattr(r,"__len__") else r, [str(x.message)[:60] for x in w])
    except Exception as e:
        print(label,"-> RAISES",type(e).__name__,str(e)[:110])
# kB handling at T>0 on ED
run("ED T=2 kB default", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=2.,**kw))
run("ED T=2 kB=-8.6e-5", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=2.,kB=-8.617333262e-5,**kw))
run("ED T=2 kB=0", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=2.,kB=0.,**kw))
run("ED T=nan", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=float("nan"),**kw))
run("ED T=-1", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=-1.,**kw))
run("TASK T=0.2 kB=-8.6e-5", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.2,kB=-8.617333262e-5,nkry=30,**kw))
run("TASK T=0.2 kB=0", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.2,kB=0.,nkry=30,**kw))
run("TASK T=0 kB=-1 (inert)", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,kB=-1.,nkry=30,**kw))
# mode typo with LL
run("mode=ed order=LL", lambda: sc.get_kondo_spectrum(eV,mode="ed",T=0.,**dict(kw,order="LL")))
run("mode=ed order=3", lambda: sc.get_kondo_spectrum(eV,mode="ed",T=0.,**kw))
# KPM with nkry
run("KPM nkry", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,submode="KPM",nkry=30,delta=1e-4,es=np.linspace(-.03,.03,100),**dict(kw,order=2)))
# order=2 inert args
a=sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,nkry=30,**dict(kw,order=2))[1]
b=sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,nkry=30,**dict(kw,order=2,Jrho_s=5.,omega0=1.,Gamma0=1e-2))[1]
print("order=2 insensitive to Jrho_s/omega0/Gamma0:", np.max(np.abs(a-b)))
# bool order / order=True
run("order=True", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,order=True)))
run("order=3.0", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,order=3.0)))
run("order=np.int64(3)", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,order=np.int64(3))))
run("site=-1 ED", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,site=-1)))
run("site=-1 TASK", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T=0.,nkry=30,**dict(kw,site=-1)))
run("site=4", lambda: sc.get_kondo_spectrum(eV,mode="ED",T=0.,**dict(kw,site=4)))
run("T=1 string", lambda: sc.get_kondo_spectrum(eV,mode="DMRG",T="1",nkry=30,**kw))
run("empty eV TASK T>0", lambda: sc.get_kondo_spectrum(np.array([]),mode="DMRG",T=0.1,nkry=30,**kw))
```

`12_cheap_checks.py`:

```python
import numpy as np, warnings
from dmrgpy import spinchain
KB=8.617333262e-5
n=4
def mk(version="python", nn=n):
    sc=spinchain.Spin_Chain(["S=1/2"]*nn,itensor_version=version)
    h=0
    for i in range(nn-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    h=h+4e-4*sum(sc.Sz[i] for i in range(nn))
    sc.set_hamiltonian(h); sc.maxm,sc.nsweeps,sc.cutoff=16,20,0.
    return sc
eV=np.linspace(-4e-3,4e-3,41)
kw=dict(site=0,Jrho_s=-0.1,U=0.3,omega0=20e-3,Gamma0=5e-5,order=3)
sc=mk(); T=0.3
def rel(a,b): return np.max(np.abs(a-b))/np.max(np.abs(b))
_,ed=sc.get_kondo_spectrum(eV,mode="ED",T=T,**dict(kw,order=2))
_,tk=sc.get_kondo_spectrum(eV,mode="DMRG",T=T,nkry=60,**dict(kw,order=2))
print("order=2 T=0.3 TASK vs ED:",rel(tk,ed))
# kB=1 with T = kT
kT=KB*T
_,a=sc.get_kondo_spectrum(eV,mode="DMRG",T=T,nkry=60,**kw)
_,b=sc.get_kondo_spectrum(eV,mode="DMRG",T=kT,kB=1.,nkry=60,**kw)
_,c=sc.get_kondo_spectrum(eV,mode="ED",T=kT,kB=1.,**kw)
_,d=sc.get_kondo_spectrum(eV,mode="ED",T=T,**kw)
print("TASK kB=1,T=kT vs default:",rel(b,a)," ED:",rel(c,d)," TASK vs ED:",rel(a,d))
# short-chain fallback with kB=1
s1=mk("python",1)
_,a=s1.get_kondo_spectrum(eV,mode="DMRG",T=T,**kw); _,b=s1.get_kondo_spectrum(eV,mode="DMRG",T=kT,kB=1.,**kw)
print("short-chain fallback kB=1 vs default:",rel(b,a))
# LL refusals as reached
for k in (dict(mode="DMRG",T=0.),dict(mode="ED",T=0.3),dict(mode="DMRG",T=0.3)):
    try: sc.get_kondo_spectrum(eV,**dict(kw,order="LL"),**k)
    except Exception as ex: print("LL",k,"->",type(ex).__name__,str(ex)[:70])
# n_gs invalid is checked after FBuilder is built (precondition after the work)
import time
t=time.time()
try: sc.get_kondo_spectrum(eV,mode="DMRG",T=T,nkry=60,n_gs=0,**kw)
except Exception as ex: print("n_gs=0 ->",type(ex).__name__,"after %.2f s"%(time.time()-t))
t=time.time()
try: sc.get_kondo_spectrum(eV,mode="DMRG",T=T,nkry=60,n_gs=0,submode="TASK",**kw)
except Exception as ex: pass
```

Observed:

`02_dispatch_checks.out`:

```
BACKENDS ('python', 3, 'julia_live')
ED T=2 kB default -> [8.87539 4.42909 8.67066] []
ED T=2 kB=-8.6e-5 -> RAISES ValueError `x` must contain only finite values.
ED T=2 kB=0 -> RAISES ValueError Geometric sequence cannot include zero
ED T=nan -> RAISES ValueError `x` must contain only finite values.
ED T=-1 -> RAISES ValueError KondoSpectrum requires T>=0
TASK T=0.2 kB=-8.6e-5 -> RAISES ValueError `x` must contain only finite values.
TASK T=0.2 kB=0 -> RAISES ValueError Geometric sequence cannot include zero
TASK T=0 kB=-1 (inert) -> [8.88649 2.26195 8.66271] []
mode=ed order=LL -> RAISES NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kondo_leading_log/README.md
mode=ed order=3 -> RAISES ValueError mode must be "ED" or "DMRG"
KPM nkry -> RAISES TypeError Unexpected keyword argument(s) for the KPM dynamical correlator: nkry. Solver parameters are attributes of the
order=2 insensitive to Jrho_s/omega0/Gamma0: 0.0
order=True -> RAISES ValueError order must be 2, 3 or "LL"
order=3.0 -> [8.88649 2.26195 8.66271] []
order=np.int64(3) -> [8.88649 2.26195 8.66271] []
site=-1 ED -> [8.88649 2.26195 8.66271] []
site=-1 TASK -> [8.88649 2.26195 8.66271] []
site=4 -> RAISES IndexError list index out of range
T=1 string -> RAISES TypeError '>=' not supported between instances of 'str' and 'float'
empty eV TASK T>0 -> RAISES IndexError index 0 is out of bounds for axis 0 with size 0
```

`12_cheap_checks.out`:

```
order=2 T=0.3 TASK vs ED: 2.0520391785564677e-05
TASK kB=1,T=kT vs default: 0.0  ED: 0.0  TASK vs ED: 2.6448305466596026e-05
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
short-chain fallback kB=1 vs default: 0.0
LL {'mode': 'DMRG', 'T': 0.0} -> NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kond
LL {'mode': 'ED', 'T': 0.3} -> NotImplementedError order="LL" is T=0 only for now (pass T=0); see docs/new_methods/kondo_
LL {'mode': 'DMRG', 'T': 0.3} -> NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kond
n_gs=0 -> ValueError after 0.40 s
```

**Reviewer (CONFIRMED (narrowed))**: Read the code: kondomethod.py checks `order == "LL" and mode != "ED"` (line ~214) before the if mode == "ED" / elif "DMRG" / else ValueError chain (~249-255), so any mode string other than "ED" with LL gets the NotImplementedError. This is the development.md 4.10 pattern (one raise serving unsupported and mistyped), though the wrong message is still loud and names mode="ED", so a user can recover. For n_gs, the FBuilder in the TASK branch is built only `if T > 0. and order == 3` (line ~319-322) and the n_gs check is at ~402, so the 0.40 s the hunter measured is reproduced only at T>0, order=3; at T=0 or order=2 the check is already cheap. The 0.4 s is a once-per-call cost on a mistyped argument, so the cost sub-claim is real but trivial. The hunter's phrasing implied it applies to every route; narrowed to T>0, order=3 (the ED-route FBuilder at line 270 is not reached from the DMRG n_gs check at all).

Reviewer's own evidence:

```
kondomethod.py ~214: `if order == "LL" and mode != "ED": raise NotImplementedError(...)` precedes `elif mode == "DMRG": ... else: raise ValueError("mode must be \"ED\" or \"DMRG\"")`; ~321-322 `Fb = FBuilder(...)` guarded by `T > 0. and order == 3`; ~402-406 n_gs check after it.
```

Reviewer's scripts and outputs:


`02_dispatch_checks.out`:

```
BACKENDS ('python', 3, 'julia_live')
ED T=2 kB default -> [8.87539 4.42909 8.67066] []
ED T=2 kB=-8.6e-5 -> RAISES ValueError `x` must contain only finite values.
ED T=2 kB=0 -> RAISES ValueError Geometric sequence cannot include zero
ED T=nan -> RAISES ValueError `x` must contain only finite values.
ED T=-1 -> RAISES ValueError KondoSpectrum requires T>=0
TASK T=0.2 kB=-8.6e-5 -> RAISES ValueError `x` must contain only finite values.
TASK T=0.2 kB=0 -> RAISES ValueError Geometric sequence cannot include zero
TASK T=0 kB=-1 (inert) -> [8.88649 2.26195 8.66271] []
mode=ed order=LL -> RAISES NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kondo_leading_log/README.md
mode=ed order=3 -> RAISES ValueError mode must be "ED" or "DMRG"
KPM nkry -> RAISES TypeError Unexpected keyword argument(s) for the KPM dynamical correlator: nkry. Solver parameters are attributes of the
order=2 insensitive to Jrho_s/omega0/Gamma0: 0.0
order=True -> RAISES ValueError order must be 2, 3 or "LL"
order=3.0 -> [8.88649 2.26195 8.66271] []
order=np.int64(3) -> [8.88649 2.26195 8.66271] []
site=-1 ED -> [8.88649 2.26195 8.66271] []
site=-1 TASK -> [8.88649 2.26195 8.66271] []
site=4 -> RAISES IndexError list index out of range
T=1 string -> RAISES TypeError '>=' not supported between instances of 'str' and 'float'
empty eV TASK T>0 -> RAISES IndexError index 0 is out of bounds for axis 0 with size 0
```

`12_cheap_checks.out`:

```
order=2 T=0.3 TASK vs ED: 2.0520391785564677e-05
TASK kB=1,T=kT vs default: 0.0  ED: 0.0  TASK vs ED: 2.6448305466596026e-05
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
short-chain fallback kB=1 vs default: 0.0
LL {'mode': 'DMRG', 'T': 0.0} -> NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kond
LL {'mode': 'ED', 'T': 0.3} -> NotImplementedError order="LL" is T=0 only for now (pass T=0); see docs/new_methods/kondo_
LL {'mode': 'DMRG', 'T': 0.3} -> NotImplementedError order="LL" is implemented on mode="ED" only; see docs/new_methods/kond
n_gs=0 -> ValueError after 0.40 s
```

**Suggested fix (hunter)**: Validate mode (membership in ("ED","DMRG")) first, then the LL checks; move the n_gs validation to the top of _get_kondo_spectrum_dmrg before Fb is built. Changes no numbers.

**Fix (reviewer)**: Hunter fix is right. In get_kondo_spectrum, test `mode not in ("ED","DMRG")` first (ValueError), then the LL checks. Move the n_gs validation to the top of _get_kondo_spectrum_dmrg (before submode resolution and the FBuilder). Pin with a test that mode="ed", order="LL" raises ValueError, and that n_gs=0 with T>0 raises before FBuilder is called (monkeypatch FBuilder to fail).

## Ruled out

Two candidates were refuted. They are kept here, with the reviewer's evidence, so the next hunt does not raise them again.

**Refuted, `correlator-twotime-3`** (lens `correlator-twotime`): The two-time Kondo term has no check that dt2 and dtau satisfy the Nyquist condition of the spectrum: the grid used in the docstring and the tests (dt2=5e5, dtau=2e4 on a 3.2 meV spectral span, Nyquist dt about 968) returns a plausible curve 253% off the exact ED Lehmann term, 207% at dt2=5000, and only converges below dt2 about 1000 (6.9%, 2.3%, 1.2% at dt2=1000, 500, 250).

Reviewer: I reproduced the numbers: 2.53 and 0.0691 of the exact peak, with the span 3.245 meV and Nyquist dt of 968. The size is real, but the finding fails on three other counts.

(1) Documented. The docstring of two_time_kondo_term_dmrg says the grid has no default because none is both correct and practical. It requires dt2 ~ 1/omega0 and a t2 range of about 25/Gamma0. It tells the user to converge the grid explicitly against edtwotimeref and to treat coarse results as 'indicative only'. It calls the test grid (n_t2_half=10, n_tau_half=15) 'a small, fast, deliberately-coarse example -- not a recommended production resolution'. It also says a badly under-resolved grid 'returns a finite but silently, badly wrong result instead of erroring (confirmed directly)'. docs/guide/kondo.md says the grids 'must resolve 1/omega0 and 1/Gamma0' and that 'Too coarse a grid gives a wrong result without any error, so compare with mode="ED" on a small chain first'.

(2) Wrong reading of the claim. The route computes the Riemann-sum construction of the third-order term on the grid it is given. The test is named 'matches grid_consistent_ed_reference' on purpose: it pins that the TDVP machinery reproduces the ED two-time G on the same grid, not that the grid is converged. The hunter's own run shows the DMRG machinery is not at fault. The test grid dt2=25/Gamma0/10=5e5 and the hunter's dt2=1.25e4 both violate the documented dt2 ~ 1/omega0 = 500 by 25 to 1000 times. Not satisfying the documented rule is a misuse, not a defect. The statement that the 'documented grid' is aliased is false: the docstring documents it as a non-production example.

(3) The only surviving part is a gap in the documentation. The stated rule (1/omega0, 1/Gamma0) does not name the spectral span of S_k|GS> (dt below pi/span) as a separate requirement. On this chain dt2 = 1/omega0 = 500 passes Nyquist (968) anyway, and converging against ED is already instructed. An optional warning would be a feature, not a bug fix. I filed no narrowed finding. I also did not re-run the 5000/2000, 500/500 and 250/250 rows, so those three figures are the hunter's and I did not verify them.

```
01_repro.out: 'spectral span of Sk|GS> (meV): 3.2451177239731597   Nyquist dt = pi/span = 968.098208080841'; 'exact term (ED Lehmann) max 0.8657556323069852'; row '12500 10 20000 15 2.53'; row '1000 125 500 100 0.0691'. Docstring text in dmrgtwotime.py lines 241-253 and docs/guide/kondo.md lines 167-171 as quoted above.
```

`01_repro.py`:

```python
import warnings, numpy as np
warnings.filterwarnings("ignore")
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk import conductance
from dmrgpy.kondospectrumtk.twotime import kondo_term_from_two_time
from dmrgpy.kondospectrumtk.edtwotimeref import _levi_civita_coeff_G_chunk
G=2.0; MUB=5.7883818066e-5
sc=spinchain.Spin_Chain(["1/2"]*4, itensor_version=3)
h=G*MUB*10.*sc.Sz[0]
for i in range(3): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
ks=KondoSpectrum(sc,0,0.)
print("spectral span of Sk|GS> (meV):",1e3*(ks.e.max()-ks.e.min()),"  Nyquist dt = pi/span =",np.pi/(ks.e.max()-ks.e.min()))
omega0,Gamma0=2e-3,2e-4
eV=np.linspace(-4e-3,4e-3,9)
exact=conductance.third_order_kondo_dIdV(ks,eV,1.0,T0=1.,omega0=omega0,Gamma0=Gamma0)/(4*np.pi)
print("exact term (ED Lehmann) max",np.abs(exact).max())
def run(dt2,n2,dtau,n_tau):
    t2=float(dt2)*np.arange(-n2,n2+1); tau=float(dtau)*np.arange(-n_tau,n_tau+1)
    def b():
        for i in range(len(t2)):
            c=t2[i:i+1]
            yield c,_levi_civita_coeff_G_chunk(ks,c,tau),_levi_civita_coeff_G_chunk(ks,c,tau,exchange=True)
    return kondo_term_from_two_time(t2,tau,b(),eV,omega0,Gamma0)
print("exact third_order_kondo_dIdV via edkondo; ED two-time on same grid")
print("   dt2  n_t2_half  dtau n_tau_half   max|term-exact|/max|exact|")
for dt2,n2,dtau,nt in ((1.25e4,10,2e4,15),(1000,125,500,100)):
    r=run(dt2,n2,dtau,nt)
    print("%8.0f %5d %8.0f %5d      %.3g"%(dt2,n2,dtau,nt,np.abs(r-exact).max()/np.abs(exact).max()),flush=True)
```

`01_repro.out`:

```
spectral span of Sk|GS> (meV): 3.2451177239731597   Nyquist dt = pi/span = 968.098208080841
exact term (ED Lehmann) max 0.8657556323069852
exact third_order_kondo_dIdV via edkondo; ED two-time on same grid
   dt2  n_t2_half  dtau n_tau_half   max|term-exact|/max|exact|
   12500    10    20000    15      2.53
    1000   125      500   100      0.0691
```

**Refuted, `dispatch-2`** (lens `dispatch`): The T>0 TASK guard cannot see a thermally occupied state that H conserves away from the Krylov space of S_site|GS>, so it lets through a sample with 10% of its weight outside the ground state and returns dI/dV 4.5% of the peak away from mode="ED", with no warning.

Reviewer: Reproduced on v3: ED levels [0, 2.8e-4, 9.8e-4, 1.3e-3], TASK poles [9.8e-4, 1.3e-3, 1.1e-2], T=1.482 K, ED weight outside GS 0.101, max|TASK-ED|/max|ED| = 0.045, warnings 0. This agrees with the hunter. It is the documented limitation, stated in two places: the get_kondo_spectrum docstring ('The check sees the states the spin of the site reaches from the ground state, not the others') and docs/guide/kondo.md ('An excitation that this spin does not reach is not checked'). kondo.md also says to keep _check_cold built on the Krylov poles and connected_to_ground_state, and says the tangent space holds only the excitations of the ground state. The control shows nothing lies beyond the boundary: the same spectator with an XY coupling to the site is reached, the guard sees the 3.0e-4 pole and raises or lowers T to 0.457 K, and TASK matches ED to 0.000 at the guard edge. On python the Krylov start noise leaks the unreached state in on my run (TASK poles [3.0e-4...], guard catches it), the hunter's python run was blind, v3 is blind in both; that is a property of the same boundary (no guarantee either way), not a second defect. The unmeasured consequence (a Hubbard chain with no Nf sector, ED grand-canonical over N) is speculative and was not measured by the hunter; struck. What remains is only an optional doc improvement: a measured size next to the sentence. Nothing in the code is wrong.

```
03_spectator_v3.out: 'Ising (spectator conserved) | ED levels [0. 0.00028 0.00098 0.0013 ] | TASK poles [0.00098 0.0013 0.01098]', 'T=1.482 K, ED weight outside GS 0.101', 'max|TASK-ED|/max|ED| = 0.045, warnings 0'; 'with XY (spectator reached) | ... TASK poles [0.000302 ...]', 'ED weight outside GS 0.000', 'max|TASK-ED|/max|ED| = 0.000'. 03_spectator.out (python): the Ising case already sees the 2.8e-4 pole, guard correct.
```

`03_spectator.py`:

```python
# own repro of dispatch-2, plus the control: the same spectator coupled by an XX term so the spin reaches its flip
import numpy as np, warnings
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import taskkondo
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
import sys; VER=sys.argv[1] if sys.argv[1]=="python" else int(sys.argv[1])
KB=8.617333262e-5
eV=np.linspace(-2.5e-3,2.5e-3,251)
kw=dict(site=0,Jrho_s=-0.1,U=0.0,omega0=20e-3,Gamma0=2e-5,order=2)
for label,Jxy in (("Ising (spectator conserved)",0.),("with XY (spectator reached)",3e-4)):
    sc=spinchain.Spin_Chain(["S=1/2"]*3,itensor_version=VER)
    H=1e-3*sc.Sz[0]+6e-4*sc.Sz[0]*sc.Sz[1]+2e-5*sc.Sz[1]+1e-2*sc.Sz[2]
    if Jxy: H=H+Jxy*(sc.Sx[0]*sc.Sx[1]+sc.Sy[0]*sc.Sy[1])
    sc.set_hamiltonian(H); sc.maxm,sc.nsweeps,sc.cutoff=8,20,0.
    e=np.sort(np.linalg.eigvalsh(sc.get_ED_obj().get_operator(sc.hamiltonian).todense())); e-=e[0]
    ks=TaskKondoSpectrum(sc,0,nkry=100)
    poles=np.sort(ks.e[1:][ks.e[1:]>1e-9*ks.e.max()])
    T=0.9*poles[0]/(KB*np.log(1/taskkondo.THERMAL_OCCUPATION)); kT=KB*T
    p=np.exp(-e/kT); p/=p.sum()
    print(label,"| ED levels",np.round(e[:4],6),"| TASK poles",np.round(poles[:3],6))
    print("  T=%.3f K, ED weight outside GS %.3f"%(T,1-p[0]))
    _,ed=sc.get_kondo_spectrum(eV,mode="ED",T=T,**kw)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            _,tk=sc.get_kondo_spectrum(eV,mode="DMRG",submode="TASK",nkry=100,T=T,**kw)
            print("  max|TASK-ED|/max|ED| = %.3f, warnings %d"%(np.max(np.abs(tk-ed))/np.max(np.abs(ed)),len(w)))
        except ValueError as ex: print("  RAISED:",str(ex)[:110])
```

`03_spectator_v3.out`:

```
Ising (spectator conserved) | ED levels [0.      0.00028 0.00098 0.0013 ] | TASK poles [0.00098 0.0013  0.01098]
  T=1.482 K, ED weight outside GS 0.101
  max|TASK-ED|/max|ED| = 0.045, warnings 0
with XY (spectator reached) | ED levels [0.       0.000302 0.001025 0.001322] | TASK poles [0.000302 0.001025 0.001322]
  T=0.457 K, ED weight outside GS 0.000
  max|TASK-ED|/max|ED| = 0.000, warnings 0
```

`03_spectator.out`:

```
Ising (spectator conserved) | ED levels [0.      0.00028 0.00098 0.0013 ] | TASK poles [0.00028 0.00098 0.0013 ]
  T=0.423 K, ED weight outside GS 0.000
  max|TASK-ED|/max|ED| = 0.000, warnings 0
with XY (spectator reached) | ED levels [0.       0.000302 0.001025 0.001322] | TASK poles [0.000302 0.001025 0.001322]
  T=0.457 K, ED weight outside GS 0.000
  max|TASK-ED|/max|ED| = 0.000, warnings 0
```

Checked by the hunters and found clean, with their numbers:

- **task-route**: TASK(T>0) against mode='ED' at saturation, order 3 with U=0.2, 3-site Hubbard b=0.1, both spinful classes, 'python' and v3, kB=1: kT = gap/45 gives 8.3e-13 (interleaved python), 3.9e-11 (native python), 4.9e-12 (interleaved v3), 2.6e-14 (native v3); kT = gap/12 gives 3.95e-6 on all four, matching the Boltzmann factor exp(-12) = 6e-6 of the dropped occupation. (04_fermion_T.out)
- **task-route**: The three raises of _check_cold fire as intended in everything I ran: lowest pole occupied (0.82, 0.94, 0.036, 0.016, 0.003 Boltzmann factors all raised against the 1e-3 bound); at charge gap 0.01 (occupation 2.4e-4) a 4-site Hubbard with a near-degenerate N=3 state passed with 5.2e-4 of ED, consistent with the stated bound. (02, 03)
- **task-route**: Cache key of session_spectrum: the digest covers the site arrays of the state, of H and of every operator (so site and operator change it) and the settings (nkry, pad, expand); T, kB and the bias are not inputs of the spectrum, so a call with another T or bias correctly reuses it. Read-only review, plus 10_end_to_end.py repeats with identical outputs on the cached path.
- **task-route**: Absolute thresholds in energy: degeneracy_tol (relative to max|e|), the 1e-14 weight cut (relative to total), CONVERGED_RESIDUAL (relative to max|e|), connected_to_ground_state tol (relative to the largest element), _check_cold near and THERMAL_OCCUPATION (relative to kT), 1e-300 floors only. No absolute scale found in taskkondo.py; everything ran with J=1e-3 eV and with kB=1 units.
- **task-route**: Duplicate Sx/Sy/Sz starts and their daggers: block_krylov already deflates them and sizes the store by the surviving count (code and profile: 107 apply_many calls for nkry=100), so the 6 starts cost no more than 3.
- **task-route**: convergence() cost at T>0 (the apart loop with a sort per pole) is negligible: not visible in the top 30 of the profile at dim 904.
- **ed-reference**: Rounding at zero bias at T=0: a ground doublet (3 and 5 spins 1/2) with eigenvalue noise injected at 1e-17, 1e-15 and 1e-13 gives dI/dV(-1e-9, 0, +1e-9) identical to 7 digits (3.337933 and 1.842207, script 01), because Theta0(eps)+Theta0(-eps)=1 and the coefficient sums are symmetric in the degenerate pair. A split of 5e-11 eV (above the 1e-9 relative tolerance, here 2e-11 eV) drops dI/dV(0) from 2.0088 to 1.7168 while +-1e-9 keeps 2.0088: that is the exact T=0 step of a split doublet, not a rounding artefact (script 02).
- **ed-reference**: Scale covariance of the whole ED result (H, eV, omega0, Gamma0 and T times s with kB fixed, s = 1e-3, 0.37, 1e3, U=0.25, S=1/2 plus S=1 chain with fields, T = 0, 0.3, 4 K): max relative change 6.7e-16 at T=0, 1.9e-11 at 0.3 K and 1.3e-11 at 4 K (script 04). P_CUT, degeneracy_tol, the FBuilder grids and the Theta series cut are all relative or dimensionless; no absolute threshold in eV found in edkondo.py, conductance.py or stepfunctions.py (Gamma0 is a parameter, not a threshold).
- **ed-reference**: T0 meaning: T0=2 against T0=1 gives exactly 4 (max |b/a-4| = 0.0) for order 2 at T=1 K, order 2 with U at T=0, order 3 with U at T=1 K and T=0, and order="LL" with U at T=0 (script 20); conductance.py, secondorder_dc.py, potentialdc.py and leadinglog.py all use 2*pi*T0^2 (second order, LL) and 4*pi*T0^2*Jrho_s (third order).
- **ed-reference**: T to 0 continuity of the public ED spectrum (5-spin chain, ground doublet, six biases away from the steps): max |d(T)-d(0)| = 7.7e-8 at 1e-9 K, 1.6e-8 at 1e-6 K, 1.1e-7 at 1e-3 K, 1.1e-5 at 1e-2 K (script 20); the 1e-8 floor is the F table accuracy.
- **ed-reference**: FBuilder table against a 30-digit quadrature at the default regime: F(0) error 3.5e-9 (1 K), 1.3e-10 (0.1 K), 4.1e-8 (5 K), 3.7e-8 (30 K) at omega0=20 meV, Gamma0=5 ueV (script 05); the T->0 limit approaches F0 smoothly (0.01 K: |F0-F| 7.7e-2 only within kT of x=0). Theta against d/dx[x/(1-e^-x)] to 2.8e-6 relative at x=-30 (cancellation of 0.5 - 0.5, 1e-16 absolute); Theta' integrates to 1.0000048 over [-80,80] (the Theta_prime band, candidate 5). Script 03's scipy-quad reference was itself wrong by 5e-4 and is not cited.
- **ed-reference**: Grids: a shuffled and duplicated bias list gives the same values per point (max difference 0.0 at T=0, 4.4e-16 at T=1 K, n=9 where the F chunk is 15 biases), a single bias point equals the same point inside the grid, and a non-uniform grid is per-point independent in the sums read (script 11). A scalar eV (0-d) raises IndexError 'too many indices for array' on both temperatures: visible, cheap, not reported.
- **ed-reference**: Occupied-state cut: p > 1e-12*max(p) is relative; n_occ = 4, 127 or 176, 512 or 917 states at n=9/10 for 1 K, 5 K, 10 to 20 K (script 14); it is the memory driver of candidate 3, not a wrong number.
- **ed-reference**: Second-order U term: 4*U^2 is added as a constant with Theta(eV)+Theta(-eV)=1, which also holds at T=0 (Theta0(0)+Theta0(-0)=1); the potential term is odd in eV with F even so it is exactly 0 at eV=0 (read and checked through the parity of the code, paper tests pin the values).
- **leading-log**: First order against conductance.py order 3 (second+third+potential, U=0.25) on models the tests do not use, 01_first_order_generic.out: random Hermitian triples d=4 nondegenerate (g^2 ratio 3.98, odd ratio 8.22), d=5 with a 2-fold ground state (4.17, 8.22), S=1 D=-1 meV B=0 doublet (4.11, 8.01), S=3/2 with D and Bx (4.16, 8.03), spin-1/2 trimer B=0 doublet with tip in the middle (4.12, 8.02) and at the end (4.17, 8.03). The expected values are 4 and 8. Only the 3-fold random row is off in the odd part (see leads).
- **leading-log**: Spinful fermion chains, both classes, 3 sites, Nf=3 sector and unrestricted Fock space (14_hubbard.out): |LL-order3| halves by 4.05 to 4.08 when Jrho_s is halved, same numbers on the interleaved and native class, no NaN.
- **leading-log**: Scale covariance (H, eV, omega0, Gamma0 times s from 1e-6 to 1e6) on the three test models, 02_scale.out: largest relative change 5.3e-11, no NaN, no warning.
- **leading-log**: Closed form of the free S=1/2 at B=0 for both signs of Jrho_s and for U=0, 0.25, -0.4 at Jrho_s=+0.3, +0.05, -0.02, -0.10 (04_closed_form_signs.out): at most 1.3e-7 relative from an independent quadrature, no stop.
- **leading-log**: U entering as A(U)=A(0)+U: confirmed by the closed form above at three values of U, and by the repo test of the cancellation.
- **leading-log**: Rotation inside a degenerate ground doublet (17_degenerate_real_and_probes.out): dI/dV changes by at most 1.2e-10 and info['coupling'] by at most 1.4e-11 (maximum coupling 0.0682 at g=-0.05 and 0.2126 at g=-0.10, no NaN), so the open-transition coupling check is basis independent in that case.
- **leading-log**: Degenerate ground state convention: leading_log_dIdV takes the occupied set from the same KondoSpectrum.p (equal weights below 1e-9 of the bandwidth) as order=3, so the two share one convention by construction; the nocc=2 and nocc=3 models above agree to first order.
- **leading-log**: Keywords of get_kondo_spectrum(order="LL"): a typo (Jrho=), nkry=, rtol= raise TypeError; Gamma0<=0 raises ValueError; the default T=1.0 raises NotImplementedError with 'pass T=0'; kB is accepted and unread at T=0 (harmless); order='ll' raises ValueError; mode='DMRG' raises NotImplementedError (18_zero_spin_public.out).
- **leading-log**: The flow's error control is not dominated by the tip quadrature: the tip block is about half of the state vector and the step counts of the real prototype (38) and the current code (41) are close; the high-D region costs 7 of 41 steps (lead above).
- **leading-log**: GPU, complex64 and numba threading of _hole_sum were not touched, per the brief.
- **correlator-twotime**: Scale covariance of the KPM second-order route: H, eV, es, delta all multiplied by 1000 give the same max difference from ED to 13 digits (0.7663215808164 against 0.7663215808164 at delta=5e-5, 0.000273163622519 at delta=2e-6), so no absolute threshold acts on the correlator route when delta and es follow H (script 02).
- **correlator-twotime**: v2, v3 and python KPM at delta=2e-6 on the 4-site chain agree with ED order=2, U=0.3 to 2.330e-4 against a peak of 6.974 (all three equal to 8 digits), so the earlier lead of v2 being 0.39 off on a 3-site Hubbard chain is broadening, not a v2 defect (script 15; at delta=5e-5 all three give 0.764, at 1e-5 0.088).
- **correlator-twotime**: KPM potential term against conductance.third_order_potential_dIdV: 0.071 of a 0.827 peak at delta=5e-5 (broadening), and the potential routine reproduces the same values with the shared-S_tot formula to 4e-15 (script 05, 12).
- **correlator-twotime**: The two-time construction itself: G and Gx against the ED two-time reference agree to 7.6e-11 (v3) and 3.2e-11 (python) with Gx bit identical between the old and the prototype, field of 10 T on the site, exchange diagram included; the grid errors in candidate 3 are the grid, since the ED reference run on the same grid shows them too.
- **correlator-twotime**: Dispatch of the frame diff in kondomethod.py: T default is 1.0 K, so mode='DMRG' with a correlator submode and no T= raises ValueError naming TASK (loud), v2 with the implicit submode takes KPM and raises the same way at T>0, NaN T is rejected by 'not T >= 0.'; no silent path found.
- **correlator-twotime**: Spectral weights below a threshold at the potential term: with es past the top of the spectrum the sum-rule warning in potentialdc does fire only through the Casimir comparison; no case of a wrong potential term without warning was found when es is wide enough.
- **dispatch**: T>0 order=2 TASK against mode=ED, 4-site chain, T=0.3 K: max relative difference 2.05e-05 (12_cheap_checks.out), no Fb is needed or built.
- **dispatch**: kB=1 with T given as kT against the default kB: bitwise equal (0.0) on TASK, ED and the short-chain fallback; kB is consumed by FBuilder, TaskKondoSpectrum and the guard consistently.
- **dispatch**: Short-chain ED fallback boundary (09_short_chain_boundary.out): v3 ns=2 and python ns=1 return the exact thermal ED spectrum at any T with a printed fallback note (max difference 0.0 against mode=ED, warm T=2 K included), while v3 ns=3 and python ns=2 raise at the same warm T; at a cold T the ns=3 TASK result agrees with ED to 1.0e-07 (v3, gap/kT=15.5) and 8e-11 (python ns=2, gap/kT=23). The switch is visible by a print, not by a warning, and the two sides are consistent where the guard passes.
- **dispatch**: Guard strength on dense low spectra (07, 08): the guard bounds only the lowest pole, yet at 0.99 of its edge the neglected excited occupation Z-1 is 0.001 and the error against ED is 6.1e-04 (8 sites) and 5.5e-04 (10 sites) of the peak, 6.1e-05 on a 9-site chain with a dense low band; at 0.5 of the edge 6e-07. The number of states below 5 kT was 0 in every case, so the sum over unreached poles does not accumulate in these models.
- **dispatch**: T>0 with n_gs=2 on a Kramers doublet agrees with ED to 1e-9 in the existing tests (test_kondo_task_convergence.py) and each member builds its own TaskKondoSpectrum with n_gs; the FBuilder is built once per call and shared by the members (confirmed in kondomethod.py:318-322 and in the 13-site profile, 1 FBuilder call).
- **dispatch**: order=2 ignores Jrho_s, omega0 and Gamma0 on TASK (difference 0.0 when Jrho_s=5, omega0=1, Gamma0=1e-2 are passed); this is the documented meaning of order=2 and the same on ED.
- **dispatch**: Unknown or grid keywords: submode="KPM" with nkry= raises TypeError naming nkry; ED raises TypeError on unread keywords; TASK raises on es/delta/dt2/n_t2_half/dtau/n_tau_half; nkry validation runs before the short-chain fallback; the non-TASK DMRG routes raise at T>0 with a message that truthfully lists python, 3 and julia_live (BACKENDS) as the TASK backends. kB, T0, omega0, Gamma0 are consumed on every TASK branch; on the non-TASK DMRG routes kB is accepted and inert because T must be 0.
- **dispatch**: T validation: T=-1 raises on ED and DMRG, T=nan raises (via FBuilder, 'x must contain only finite values'), T as a string raises TypeError; order=True raises, order=3.0 and np.int64(3) are accepted and give the order=3 result; negative site indexes work like the list index (site=-1 equals the mirror site on a symmetric chain) and site=4 on 4 sites raises IndexError.

## Left open, and new leads

None of these was reviewed. Leads from the hunters:

- **task-route**: Block Gram-Schmidt in task.block_krylov.add(): project_rows is 22% of the T=0 20-site call (3.5 of 16 s) and runs two gemv passes per vector; projecting the whole block against the existing store with two gemm passes and orthogonalizing sequentially only inside the block could cut it, at roundoff changes to Ritz energies. Not prototyped.
- **task-route**: third_order_kondo_dIdV at T>0 is 11.1 s = 28% of the nkry=300 call (FBuilder.__call__ 10.3 s: scipy interpolate _evaluate 5.5 s, F0 3.7 s) after the second-order fix; belongs to the conductance lens, untouched.
- **task-route**: When a DMRG start lands in a mixture of two electron numbers (Hubbard 4 sites, N=4 and N=3 states 0.003 to 0.007 apart, no Nf sector), _check_cold raises, but the message says 'the lowest excitation the spin of site 0 reaches' when the cause is a charge mixture that S_site conserves; the message does not point at set_conserved_sector(Nf=). Not a wrong number.
- **task-route**: _check_cold window near = 1e-3 kT is narrower than a truncated-maxm partner splitting (1.5e-9 eV per its own comment) once T < about 17 mK at J of order meV, so a Kramers doublet with n_gs=2 raises 'lowest excitation ... Boltzmann factor 1' at very low T. By design loud; not run.
- **task-route**: Blind spot of the T>0 guard by construction: only states S_site reaches from the ground state are seen (documented in the docstring). I could not build a case where DMRG converged cleanly to the ground state and an unreached low state was occupied: a mirror-symmetric 5-site chain with an odd-parity state 1.7e-6 eV above the ground state made the guard raise (the converged state mixed them), and the charge-mixture runs raised as well. A clean exact degeneracy with an unreached partner of inequivalent spectrum (n_gs=1) was not constructed.
- **task-route**: v3 was not profiled (d22 has no compiled v3 here); only 'python' timings above.
- **ed-reference**: Parameter validation is inconsistent between T=0 and T>0 (script 17, read-only judgement: not a hole by the standard): Gamma0=0 gives NaN at eV=0 at T=0 (a true log divergence) but a cryptic scipy ValueError ('y must contain only finite values') at T>0; omega0<0 gives all-NaN at T=0 and the same ValueError at T>0; Gamma0<0 is silently the same as Gamma0>0 (only Gamma0^2 enters); omega0=0 is accepted and finite except at the zero-bias point. A one-line check in get_kondo_spectrum would make both temperatures say the same thing.
- **ed-reference**: Levi-Civita coefficient build (conductance._triple_product_coefficients) is 3.61 s of 10.6 s at n=9, T=20 K (34%), via einsum with the 27-entry epsilon tensor over complex arrays. Writing the six non-zero terms explicitly, and using real arithmetic for Sx and Sz when the chain is real, could roughly halve it; not prototyped.
- **ed-reference**: The TASK route on a chain too short for the session (kondomethod.py:331-336) falls back to the full thermal KondoSpectrum at T>0, while on a longer chain the same call keeps the sample in its ground state and raises through the guard when that is not the thermal state; the fallback is the more exact answer, but the same call means two different physical states depending on the length. Read only, not run (the task-route lens owns the diff).
- **ed-reference**: degeneracy_tol at T=0 is relative to the bandwidth e[-1]-e[0] after the ground state is shifted to 0, so it is shift invariant by construction; invariance of the ED result under a constant shift of H was not tested (can only break above about 1e7 bandwidths).
- **ed-reference**: Stale docstring text in conductance.third_order_kondo_dIdV: 'building it tabulates _band_integral over hundreds of adaptive-quadrature points' (there is no _band_integral; the table is built by _F_convolve with Simpson, 0.24 s at default). Also FBuilder.__call__ docstring: '1e6 points in 0.15 s' not re-measured here.
- **ed-reference**: p_cut (P_CUT=1e-12) acts on p_i relative to max(p) so it is relative, but the coefficient tensor memory scales with n_occ (candidate 3); the cut is a memory knob as well as an accuracy knob, and a user-facing p_cut keyword does not exist on get_kondo_spectrum.
- **leading-log**: Real-path generalization: a complex H whose eigenbasis is real up to a diagonal gauge, or a By field, falls back to the complex code. A gauge fixing of the eigenvectors (make the largest component of each real) could widen the real fast path; not tried.
- **leading-log**: Start the flow at Dmax of about 1e2 to 1e3 times the scale and add the first-order closed-form tail (kernel integrated from Dmax' to infinity times the bare bracket): 7 of 41 accepted steps at 5 spins lie above D/scale=100 (06_steps.out), so about 15% of the evaluations; not prototyped, the error is O(g * omega0/Dmax') and needs a measurement.
- **leading-log**: The default rtol=1e-10, atol=1e-14 is 9 orders beyond the accuracy of the method; rtol=1e-8 halves the steps (42 to 25 at 5 spins, 1.6 to 1.7x end to end at 6 and 7 spins) but moves numbers 5e-10 to 2.7e-9 and the pinned _TWO_PASS test is rel=1e-9, so it is not offered as a candidate; if the author accepts a looser default, the tests' tolerance has to follow.
- **leading-log**: get_kondo_spectrum(order="LL") cannot reach rtol, atol, pole_threshold, stop_margin, p_cut or return_info (the ED route refuses unknown keywords with TypeError, which is loud), so D_stop, max_coupling and info['coupling'] are unreachable from the public call: the user sees only warnings. Probed in 18_zero_spin_public.out; recorded as a limitation, not a silent wrong number.
- **leading-log**: The 3-fold degenerate random model in 01 has |odd part|/|third| = 0.20 above the repo test's 0.05 bound; the odd-part ratio is still 8.87 and the g^2 ratio 4.32. I read it as higher-order terms from operators that are not spins (random Hermitian triples), not as a defect, but did not separate the two.
- **leading-log**: The dense-spectrum NaN mask (every bias NaN when a stopped flow meets a dense spectrum) was left as recorded and not re-measured.
- **leading-log**: The sample-coupling criterion max|G_ab| >= smax is basis dependent in principle; on the SU(2)-symmetric free S=1 and S=3/2 triplet/quartet it gave D_stop equal to 1e-9 relative under random unitary rotations in the degenerate manifold (03_basis_stop.out), because G stays proportional to S there. An anisotropic degenerate manifold with a stopped flow was not tried.
- **correlator-twotime**: twotime.K_F cell averages: at the docstring/test grid (dt2 about 1e6, eV*dt2 of hundreds of rad per cell) 85% of the call (193 of 227 s on 4 sites, 3.5e6 K_F calls, IntegrationWarning 'probably divergent') goes to two scipy quad calls (real and imaginary separately) per near cell per eV per chunk; the integrand exp(i eV u) R(u) has an eV independent real R, so one complex Gauss quadrature vectorized over eV, or the QAWO weight, would remove it. Not prototyped, and that grid is aliased anyway (candidate 3), so the saving matters only for users who follow the recipe. Also K_F is rebuilt per chunk, a single point at a time, with numpy overhead per call.
- **correlator-twotime**: Gx on the two-time route: the bras psi_l(-s-tau) lie on one long trajectory when dtau is a multiple of dt2, which would remove the remaining three trajectories per row (another 3x on the TDVP part); needs a restriction on the grids, not prototyped.
- **correlator-twotime**: v3 prints 'warning: applyExp not converged in 50 steps' dozens of times on the docstring grid at 8 sites (dtau=2e4 with H of order 1e-3 eV, a step phase of tens of rad) and the result is used silently; known open applyExp lead, here with a concrete trigger, not quantified for its effect on G.
- **correlator-twotime**: submode='CVM' through the Kondo correlator route took 394 s for 4001 es points on 4 sites (v3), printing one line per energy; accuracy matched KPM at equal delta (second order 0.81 against 0.76 off ED at delta=5e-5, potential 0.17 against 0.07), so only the cost is odd. submode EX raised 'asked for 20 states, but Hilbert space has dimension 16' without n=; TD, TDZ and ROOTN were not run to the end (timeouts), their agreement with ED on the second-order and potential terms is open.
- **correlator-twotime**: twotime.kondo_term_from_two_time marks the end points with np.isclose(t2_chunk, t2_first) at the default rtol=1e-5, so for more than 1e5 half-points a neighbouring row would be halved too; unreachable with TDVP rows, only matters for the ED reference grids.
- **correlator-twotime**: n_gs>1 with the two-time term was not compared to ED: edtwotimeref._levi_civita_coeff_G_chunk reads column 0 of the eigenbasis only, so it has no manifold average to compare against on a degenerate chain.
- **correlator-twotime**: The lead 'secondorder_dc.py docstring claims 0.2% at every bias point, measures 0.46%' and the Nf-sector, nkry, JAX-memory, MPO-compression and leading-log NaN leads are outside the correlator and two-time files or were not reached.
- **dispatch**: n_gs>1 rebuilds get_excited_states on every call: 15.76 of 33.0 s (48%) on a 13-site doublet at n_gs=2, python, and the fresh members also defeat the TASK spectrum cache, so both member builds repeat (14.7 s). Caching the members under a content key (digest of the exported ground state and H) would save up to about half of a repeat call, but it changes which basis of the manifold is averaged between calls (basis-invariant only for an exactly degenerate manifold), so it needs numbers_change=true and a prototype; not built. Check SPECTRUM_CACHE_SIZE >= n_gs first.
- **dispatch**: Theta is 42% (0.77 s) of a cached 24-site T>0 repeat call, with the F spline evaluation 21% and F0 13%: conductance.py/stepfunctions.py, outside this lens (the Theta array is evaluated 10 times in _theta_raw).
- **dispatch**: Routes not profiled: mode=ED order=3 at T>0, and the KPM/correlator and two-time DMRG routes (T=0 only). The ED timing in the first draft of 06 was noise and is not cited.
- **dispatch**: Empty eV on the TASK route at T>0 raises IndexError (index 0 is out of bounds) from inside the conductance sums, not a message about eV; get_kondo_spectrum guards `if len(eV)` for the warning only. Not checked at T=0 or on ED.
- **dispatch**: kB<=0 with T>0 is not validated: kB=-8.6e-5 and kB=0 both raise ValueError from deep inside FBuilder ('x must contain only finite values', 'Geometric sequence cannot include zero') on ED and TASK. Loud but the message does not name kB.
- **dispatch**: The n_gs over-declaration at T>0 (n_gs=3 on a doublet) averages the third member, a gapped state, with weight 1/3 instead of its Boltzmann weight; the only signal is the RuntimeWarning when the split exceeds Gamma0. Not run.
- **dispatch**: A degenerate partner that differs from the DMRG state on two or more sites and is not in the one-site tangent space would not be counted by _check_cold; my attempt (accidental degeneracy of up-up-up-up with up-down-down-up in an Ising chain) always produced a superposition whose tangent space did contain the partner, so the guard raised 6 of 6 runs (python and v3); a pure product start (set_gs) was not tried.

Found by the reviewers while reviewing something else:

- **second-order-sums**: The third-order coefficient einsum _triple_product_coefficients ('jkl,ifl,fmk,mij->ifm') takes 1.5 s of a 2.6 to 2.9 s third_order_kondo_dIdV call at dim 384 with 94 occupied states (script 04_pieces), over half of the call at T>0. A direct expression of the six non-zero Levi-Civita terms as batched matmuls could remove most of it. Untested lead, the hunter listed it too. Evidence: <scratch>/review/second-order-sums/04_pieces.out
- **stepfunctions-edkondo**: Peak memory of FBuilder construction is 384 MB (tracemalloc, T=1 K, default npts=4001) for a table that retains 100 KB: _F_convolve builds a chunk of int(2e7)//npts rows by npts floats several times over (y, integrand, log temporaries). A smaller chunk would cut the transient peak with no change in numbers; not reviewed further.
- **leading-log**: The leading-log prototype in the hunter's llreal.py dropped the `Gamma0 > 0` ValueError of _joint_flow, so tests/test_kondo_leading_log.py::test_gamma0_must_be_positive fails with it (it is a prototype gap and the repo is not affected). A port of leading-log-1 must keep that check at the top of the new flow. Evidence: <scratch>/review/leading-log/04_pinned_tests.out ('1 failed, 21 passed').
- **task-guard**: TaskKondoSpectrum.convergence normalizes 'step' by the inelastic weight sum_k|<a|S_k|GS>|^2, while the docstring and warning text say 'of the largest value' of dI/dV; not checked here whether the two differ materially (elastic part and U terms add to the largest dI/dV).
- **task-guard**: Python backend TASK start noise leaks a conserved-spectator state into the Krylov space on some runs (my python run of the hunter's Ising model saw the 2.8e-4 pole, the hunter's python run and v3 did not), so the T>0 guard outcome on a conserved label depends on the backend; same documented boundary, but the guard result is not reproducible across runs there.

## Shared helpers

`run3.sh`, the local runner. A worktree sets `DMRGPY_SRC=<worktree>/src`.

```bash
#!/bin/bash
# Runs one python command under a workflow-wide cap of three concurrent runs.
# Usage, from the folder holding the script:
#   <scratch>/hunt/run3.sh [timeout_seconds] script.py [args] 2>&1 | tee script.out
#   (pytest too:  <scratch>/hunt/run3.sh 900 -m pytest <args>)
# The first argument is a timeout in seconds when it is an integer (default 900).
# It waits (checking every 5 s) until one of three lock slots is free, then
# runs python3 with BLAS threads pinned to 1 and this checkout's src first on
# PYTHONPATH (or DMRGPY_SRC=<worktree>/src when set). The slot is released when the run exits, also when it is killed.
LOCKDIR=<scratch>/hunt/locks
TO=900
if [[ "$1" =~ ^[0-9]+$ ]]; then TO=$1; shift; fi
while true; do
  for i in 1 2 3; do
    exec 9>"$LOCKDIR/slot$i.lock"
    if flock -n 9; then
      MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
        MPLBACKEND=Agg PYTHONPATH=${DMRGPY_SRC:-<repo>/src} \
        timeout $TO python3 "$@"
      rc=$?
      exec 9>&-
      exit $rc
    fi
    exec 9>&-
  done
  sleep 5
done
```

`rund22.sh`, the d22 runner, as it is after the hunt. During the hunt the two `rsync` lines of the run folder had no `--exclude '*.out'`, which is how the outputs of d22 runs were lost (see the introduction).

```bash
#!/bin/bash
# Runs one python command on d22 (<d22 host>) for heavy work,
# under the SAME workflow-wide cap of three concurrent runs as run3.sh (one global
# cap over both machines). Usage, from the folder holding the script:
#   <scratch>/hunt/rund22.sh [timeout_seconds] script.py [args] 2>&1 | tee script.out
# First it syncs this checkout's src/ (without the compiled mpscpp2/mpscpp3, so
# itensor_version=2/3 fall back to ED there: use "python", ED or leading-log on
# d22, and run v3 locally with run3.sh), then copies the current folder to
# d22:<d22>/runs/<folder name>/, runs python3 there pinned to
# one BLAS thread on the P-cores (taskset 0-11), and copies the folder back, so
# files the script writes appear here. A worktree sets DMRGPY_SRC=<worktree>/src
# and DMRGPY_SRCNAME=src-<cluster>, so it syncs to its own remote copy. Default timeout 540 s; keep one call
# under the 10-minute Bash limit and split longer work.
if [ -n "$DMRGPY_SRC" ] && [ -z "$DMRGPY_SRCNAME" ]; then
  echo "rund22.sh: DMRGPY_SRC is set, so set DMRGPY_SRCNAME=src-<cluster> too" >&2; exit 2
fi
HOST=<d22 host>
LOCKDIR=<scratch>/hunt/locks
TO=540
if [[ "$1" =~ ^[0-9]+$ ]]; then TO=$1; shift; fi
NAME=$(basename "$PWD")
REMOTE=<d22>/runs/$NAME
while true; do
  for i in 1 2 3; do
    exec 9>"$LOCKDIR/slot$i.lock"
    if flock -n 9; then
      ssh -o BatchMode=yes $HOST "mkdir -p <d22>/${DMRGPY_SRCNAME:-src} $REMOTE" &&
      rsync -a --delete --exclude 'mpscpp2/' --exclude 'mpscpp3/' --exclude '__pycache__/' \
        ${DMRGPY_SRC:-<repo>/src}/ $HOST:<d22>/${DMRGPY_SRCNAME:-src}/ &&
      rsync -a --exclude '*.out' ./ $HOST:$REMOTE/ &&
      ssh -o BatchMode=yes $HOST "cd $REMOTE && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<d22>/${DMRGPY_SRCNAME:-src} timeout $TO taskset -c 0-11 <d22 python> -s $(printf '%q ' "$@")"
      rc=$?
      # never copy *.out back: tee is writing the local one, and a copy
      # would replace it with the remote file
      rsync -a --exclude '*.out' $HOST:$REMOTE/ ./
      exec 9>&-
      exit $rc
    fi
    exec 9>&-
  done
  sleep 5
done
```
