# Audit, 2026-09: eight-lens hole hunt

Findings from an eight-lens automated audit of the Python layer, run 2026-09-12
on `43d1a35` (clean tree, both compiled extensions current). Each lens hunted one
class of problem; every finding below carries a repro that was actually executed
and its verbatim output, and every one was then handed to an independent reviewer
whose brief was to *refute* it. `REFUTED` findings are not reproduced here.

This file is the evidence, not a task list -- same convention as
`audit_2026_08_hole_hunt.md`: it records what was observed and how to reproduce
it, so a fix (or a decision that the behaviour is intended after all) does not
have to re-derive any of it. Nothing here is fixed yet; mark entries as they are.

## The eight lenses

| Lens | Brief |
|---|---|
| `python-backend-parity` | pyitensor vs ED vs v3 parity on the public method surface |
| `wavefunction-consumers` | Consumers of the DMRG wavefunction under the ~1e-6 Lanczos eigenvector cap |
| `dispatch-matrix` | Backend x method dispatch matrix, and dispatch sites added since the last audit |
| `pyitensor-performance` | pyitensor hot spots and superlinear scaling |
| `cpp-v3-completeness` | mpscpp3 completeness against mpscpp2 and against the pyitensor reference |
| `ed-and-operators` | ED backend, MultiOperator algebra and Jordan-Wigner correctness |
| `recent-commits` | Regressions and stale callers from the last month of commits |
| `docs-examples-drift` | Documentation, docstring and example drift against the code |

## Scope

Out of scope by construction, and excluded from every lens's brief: vendored
ITensor (`mpscpp2/ITensor/`, `mpscpp3/ITensor/`); the deliberately reproduced
legacy bugs CLAUDE.md lists (`evoloperator`'s z^3/6 term on `H2`, the `"moise"`
key, the unreachable `"tevol_fit_td"` branch); the `docs/known_issue_*.md`
items; anything already recorded in `audit_2026_08_hole_hunt.md`; and gaps
`ROADMAP.md` marks as absent. `itensor_version="julia_live"` was not executed
(juliacall JIT cost), so no finding below is about the Julia backend.

Every repro in this file was run with threads pinned and the checkout's own
`src/` forced onto the path, which is how they should be re-run:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 <script>
```

## Findings


### 1. pyitensor's TDVP integrator is first-order in dt, not second-order: every real-time evolution on itensor_version="python" carries an O(dt) error where itensor_version=3 is exact

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `python-backend-parity`

**Where**: `src/dmrgpy/pyitensor/tdvp.py:_half_sweep_lr / _half_sweep_rl (tdvp_step, ~line 300); reached from pyitensor/chain.py:846 quench_tdvp, :870 evolve_and_measure_tdvp, :723 tdvp_step; dispatched at src/dmrgpy/timedependent.py:169,248 and src/dmrgpy/tdz.py:168`

pyitensor/tdvp.py's module docstring and tdvp_step()'s own docstring both state the integrator is "a left-to-right half-sweep evolving by dt/2, then a right-to-left half-sweep evolving by another dt/2 -- mirrors mpscpp3/chain_session.h's tdvp_step()", i.e. a Strang-symmetric projector splitting. Two-site TDVP composed that way is second-order, and at full bond dimension (no truncation) it is EXACT by the Lubich-Oseledets property -- which is precisely what itensor_version=3 delivers (error 1e-8 -> 4e-11 as dt shrinks, i.e. dt-independent to 10 digits). itensor_version="python" instead shows an error that halves exactly when dt halves: first order. Richardson-extrapolating python's first-order sequence reproduces v3's value to 5 digits, so it is converging to the right answer, just one order too slowly.

Isolation actually performed: not truncation (identical with cutoff=1e-14, 1e-24 and 0.0; ground-state bond dims are the full 2,4,8,4,2); not maxdim (64 at n=6); not the Krylov stopping rule (errgoal forced from 1e-10 to 1e-16 changes nothing); not one-vs-two-site (num_center=1 and 2 give the same error to 2 digits); not the MPO or the environments (norm and <H> are conserved to 1e-10 along the trajectory, and python's TEBD and MPO-Taylor integrators on the same chain are bit-identical to v3's -- so the shared MPO builder, apply_mpo, environments and Hamiltonian are all fine); not the state (identical value at t=0 and identical dt->0 limit). n<=4 is exact on both backends, n>=5 breaks. That combination -- symplectic (norm and energy exactly conserved) but first-order -- is the signature of the RL half-sweep not being the exact adjoint of the LR one.

Blast radius: this is the default integrator (tevol_method="TDVP") on the backend a `pip install dmrgpy` defaults to (cppext.default_backend() returns "python" with no compiled extension). It reaches evolve_and_measure(), evolution_ABA(), evolution_DC(), get_dynamical_correlator(submode="TD") and (submode="TDZ"), tdz.py's complex-time stepper, and metts_dynamical_correlator.

Why no test caught it: tests/test_time_evolution.py's TDVP-vs-golden test uses n=2 (where two-site TDVP is trivially exact); its other python-backend evolution tests use TEBD (which is correct); test_td_dynamical_correlator_is_dt_independent asserts rel=0.05, ~3x looser than the effect; examples/time_evolution/tdvp_VS_ED_time_evolution/main.py asserts `sc.itensor_version==3` explicitly and so never exercises "python".

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/python-backend-parity && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 repro_tdvp_order.py

--- repro_tdvp_order.py (n=6 S=1/2 Heisenberg + staggered Sz + uniform Sx, maxm=64, cutoff=1e-14; quench correlator <Sz0(t)Sz0> evaluated at fixed total time T=2.0, nt=T/dt) ---
for dt in (0.2,0.1,0.05,0.025):
    for b in ("ED",3,"python"):
        sc=mk(b); sc.tevol_method="TDVP"; sc.get_gs()
        ts,cs=timedependent.evolution_DC(sc,mode=sc.get_mode(),name=(sc.Sz[0],sc.Sz[0]),nt=int(round(T/dt))+1,dt=dt)
```

Observed:

```
C(t=2.0) for the quench correlator <Sz0(t) Sz0>, full bond dim (maxm=64, exact)
dt     ED (exact)                     itensor_version=3              itensor_version='python'
0.2    (0.046375655+0.025851436j)     (0.046375644+0.025851433j) err=1.14e-08 (0.046955621+0.03808968j) err=1.23e-02
0.1    (0.046375646+0.025851434j)     (0.046375644+0.025851433j) err=2.62e-09 (0.046624126+0.03201687j) err=6.17e-03
0.05   (0.046375644+0.025851433j)     (0.046375644+0.025851433j) err=2.16e-10 (0.046489704+0.028940268j) err=3.09e-03
0.025  (0.046375644+0.025851433j)     (0.046375644+0.025851433j) err=4.11e-11 (0.046430128+0.027396689j) err=1.55e-03

--- and, via evolve_and_measure from a non-eigenstate (repro_tdvp_order2.py) ---
evolve_and_measure(<Sz0>) starting from wf = Sz[0]|gs> (NOT an eigenstate), T=1.0, maxm=64 (exact)
dt=0.2    v3=(-0.0224190323+0j)     python=(-0.0308425871-0j)      |diff|=8.424e-03
dt=0.1    v3=(-0.0224190323+0j)     python=(-0.0266667366+0j)      |diff|=4.248e-03
dt=0.05   v3=(-0.0224190323+0j)     python=(-0.0245485207+0j)      |diff|=2.129e-03
dt=0.025  v3=(-0.0224190323-0j)     python=(-0.0234847198+0j)      |diff|=1.066e-03

--- ruling out every other cause (diag_f1.py) ---
== (a) which integrator on python? N=6, T=1, ref(v3/exact)=-0.0224190323
python TDVP   dt .2/.1/.05 -> [-0.03084259 -0.02666674 -0.02454852]  errs ['8.42e-03', '4.25e-03', '2.13e-03']
python TEBD   dt .2/.1/.05 -> [-0.02270601 -0.02249067 -0.02243693]  errs ['2.87e-04', '7.16e-05', '1.79e-05']
python MPO    dt .2/.1/.05 -> [-0.02131903 -0.02035494 -0.02170778]  errs ['1.10e-03', '2.06e-03', '7.11e-04']
v3     TEBD   dt .2/.1/.05 -> [-0.02270601 -0.02249067 -0.02243693]  errs ['2.87e-04', '7.16e-05', '1.79e-05']
v3     MPO    dt .2/.1/.05 -> [-0.02131903 -0.02035494 -0.02170778]  errs ['1.10e-03', '2.06e-03', '7.11e-04']
== (b) N sweep on python TDVP (dt=0.2 vs 0.05)
N=4 python -0.020633481 / -0.020633481  (spread 1.57e-14) | v3 -0.020633481 / -0.020633481 (spread 6.78e-12)
N=5 python -0.031390709 / -0.024936152  (spread 6.45e-03) | v3 -0.022775121 / -0.022775121 (spread 3.40e-13)
N=6 python -0.030842587 / -0.024548521  (spread 6.29e-03) | v3 -0.022419032 / -0.022419032 (spread 7.73e-12)
N=7 python -0.030881737 / -0.024476923  (spread 6.40e-03) | v3 -0.022348545 / -0.022348545 (spread 6.06e-12)
== (c) norm & energy conservation along python TDVP, N=6 dt=0.1
  step 0  <psi|psi>=0.250000000000  <H>=-2.277702034585
  step 5  <psi|psi>=0.249999999825  <H>=-2.277702036179

--- cutoff / Krylov / num_center are all irrelevant ---
cutoff=1e-14 dt=0.2 err=8.424e-03 | cutoff=1e-24 err=8.424e-03 | cutoff=0 err=8.424e-03
errgoal=1e-10 dt=0.2 err=8.424e-03 | errgoal=1e-16 dt=0.2 err=8.424e-03
num_center=1 dt=0.2 err=8.450e-03 | num_center=2 dt=0.2 err=8.424e-03
```

**Expected**: Two-site TDVP composed as LR(dt/2) then RL(dt/2) is a symmetric (Strang) projector splitting and is second-order in dt; at bond dimensions large enough that no SVD truncation occurs it is exact (Lubich & Oseledets projector-splitting exactness), which is what itensor_version=3 demonstrates here (dt-independent to 4e-11). python should therefore reproduce ED/v3 to the same ~1e-10 at any of these dt, and at worst show an error that falls as dt^2. Observing an error that falls exactly as dt^1, on the same chain, same state, same maxm, same cutoff, with norm and energy exactly conserved, is a defect in the sweep composition, not a tolerance.

**Suggested fix**: Compare pyitensor/tdvp.py's `_half_sweep_lr`/`_half_sweep_rl` against mpscpp3/TDVP/tdvp.h step by step: the RL half-sweep must be the exact adjoint of the LR one (same sequence of forward two-site / backward one-site local flows, reversed), and the empirical first-order behaviour says it is not. Add a regression test that pins the order rather than a value: evolve `Sz[0]|gs>` on an n>=6 chain at maxm >= 2^(n/2) for a fixed total time at dt and dt/2 and assert the two agree to ~1e-8 (v3 satisfies this trivially; python currently differs by 4e-3). Until it is fixed, `tevol_method="TEBD"` on "python" is the correct workaround for nearest-neighbour Hamiltonians -- it is bit-identical to v3's TEBD.

**Reviewer (CONFIRMED)**: Re-ran repro_tdvp_order2.py thread-pinned: |v3-python| = 8.424e-03 / 4.248e-03 / 2.129e-03 / 1.066e-03 at dt = .2/.1/.05/.025 -- digit-for-digit the pasted output. Tried to explain it away as DMRG/tolerance error and could not: I then built a SELF-CONTAINED order test (order_selfcontained.py in my scratch dir) that uses no v3 and no ED reference at all -- python TDVP against its own dt=T/512 limit on the same n=6 chain gives err 8.34e-3, 4.16e-3, 2.05e-3, 9.82e-4, 4.50e-4 with successive ratios 2.00, 2.04, 2.08, 2.18. That is O(dt^1) on the backend's own terms, and the linear extrapolation of that sequence (err ~ 0.0417*dt, so the T/512 reference itself carries ~8e-5) lands the dt->0 limit on v3/ED's -0.0224190, i.e. python converges to the right answer one order too slowly. Bond dimensions are full at n=6 (2,4,8,4,2 vs maxm=64) so this is not truncation, and n=4 being dt-independent to 1.6e-14 rules out a time-grid/measurement-offset alternative (that would be n-independent). Not vendored code (pyitensor/tdvp.py is dmrgpy's own), not a deliberately-reproduced legacy bug, absent from docs/audit_2026_08_hole_hunt.md and from ROADMAP.md (line 81 lists two-site TDVP as fully implemented on all three backends).

WHERE I DISAGREE WITH THE HUNTER: the suggested location ("the RL half-sweep is not the exact adjoint of the LR one") is unverified and probably wrong. Reading _half_sweep_lr/_half_sweep_rl, the composition is already palindromic -- F(1,2) B(2) F(2,3) ... B(n-1) F(n-1,n) | F(n-1,n) B(n-1) ... B(2) F(1,2) -- and the environment choices are exact mirrors (LR's backward step at site i+1 uses the freshly-extended left_env[i] plus the untouched right_env[i+2]; RL's backward step at site i uses the untouched left_env[i-1] plus the freshly-extended right_env[i+1]). The hunter's own data also contradicts the attribution: num_center=1 runs _half_sweep_lr_onesite/_half_sweep_rl_onesite, entirely different functions, and shows the same 8.4e-3. The suspect set is therefore the machinery both paths share -- _lanczos_expm_multiply, _extend_left/_extend_right and the environment builders, or the driver loop between steps in pyitensor/chain.py -- none of which TEBD or MPO-Taylor use, which is consistent with those two being bit-identical to v3. Localizing further is the fixer's job; the defect itself is real and reproducible.


---


### 2. pyitensor's session keeps the previous Hamiltonian's ground state as its DMRG start, so a parameter sweep on one chain is silently trapped at a stale eigenstate

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `dispatch-matrix`

**Where**: `src/dmrgpy/pyitensor/chain.py:475-511 (set_hamiltonian / gs_energy); src/dmrgpy/manybodychain.py:702-718 (restart); src/dmrgpy/mpscpp3/chain_session.h:518-553 (same retention, latent)`

`pyitensor/chain.py::set_hamiltonian` invalidates `_solve_H_cache`, `_wf0_energy` and both bandwidth caches but leaves `self.wf0` — the previous Hamiltonian's converged MPS — in place, and `gs_energy()` then warm-starts DMRG from it (`floor_dim = _max_link_dim(self.wf0)`, `if self.wf0 is None: ...`). The Python layer believes otherwise: `Many_Body_Chain.set_hamiltonian` calls `restart()`, which sets `self.wf0 = None` and whose own comment says it "promises a genuinely cold recalculation". Nothing propagates that to the session. When the retained state happens to be an exact eigenstate of the NEW Hamiltonian — a sign flip, or any sweep that passes through a field-polarized phase, since a fully polarized product state is an exact eigenstate of the Heisenberg chain at every field — the variational solve is stationary and never moves. pyitensor's DMRG has no noise term at all (documented in `pyitensor/dmrg.py`'s module docstring and `chain.py:1613-1615`), so unlike the compiled backends it has no escape mechanism; `mpscpp3/chain_session.h::set_hamiltonian_mpo` retains `wf0_`/`have_wf0_` exactly the same way, so v3 escapes empirically (via noise) rather than structurally. This is the DEFAULT backend for a `pip install dmrgpy` (`cppext.default_backend()`), and the loop `for B in [...]: sc.set_hamiltonian(h0+B*Sz); sc.gs_energy()` is the textbook workflow. Supporting evidence that this was previously seen and misattributed: `mpsalgebra.py:29-38` works around `bandwidth()` returning <=0 with the comment that it "can occasionally underestimate a highly degenerate operator's spectral width depending on the random initial wavefunction" — it is not randomness, it is this.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 3 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/dispatch-matrix/final_A.py    # and .../sweep.py (v3 vs python field sweep), .../bw.py, .../fixtest.py
```

Observed:

```
# final_A.py (n=6 Heisenberg, maxm=40, nsweeps=40, noise=1e-2, itensor_version="python")
default backend chosen: 3
B=4 e0 = -10.75
B=0 e0 = 1.25  (maxm=40, nsweeps=40, noise=1e-2)
B=0 exact (ED)  = -2.493577
B=0 bandwidth() = 0.0

# sweep.py: same field sweep on ONE chain object, both backends
--- itensor_version = 3 --- field sweep on ONE chain object
  B= 4.0  sweep-on-one-chain e0=  -10.750000   fresh-chain e0=  -10.750000   ED=  -10.750000
  B= 1.0  sweep-on-one-chain e0=   -3.001995   fresh-chain e0=   -3.001995   ED=   -3.001995
  B= 0.0  sweep-on-one-chain e0=   -2.493577   fresh-chain e0=   -2.493577   ED=   -2.493577
--- itensor_version = python --- field sweep on ONE chain object
  B= 4.0  sweep-on-one-chain e0=  -10.750000   fresh-chain e0=  -10.750000   ED=  -10.750000
  B= 1.0  sweep-on-one-chain e0=   -1.750000   fresh-chain e0=   -3.001995   ED=   -3.001995
  B= 0.5  sweep-on-one-chain e0=   -0.250000   fresh-chain e0=   -2.501995   ED=   -2.501995
  B= 0.0  sweep-on-one-chain e0=    1.250000   fresh-chain e0=   -2.493577   ED=   -2.493577

# bw.py (n=4): Many_Body_Chain.bandwidth(h) per backend
2 fresh-clone e0,e1 = -1.616025 -0.75  same-clone a0,a1 = -1.616025 -0.75  bandwidth()= 2.366025
3 fresh-clone e0,e1 = -1.616025 -0.75  same-clone a0,a1 = -1.616025 -0.75  bandwidth()= 2.366025
python fresh-clone e0,e1 = -1.616025 -0.75  same-clone a0,a1 = -1.616025 1.616025  bandwidth()= 0.0
exact: min -1.616025 max 0.75 bandwidth 2.366025

# fixtest.py: the ONLY change is `sc._session.wf0 = None` before each gs_energy
B=4.0  with session wf0 cleared: e0=-10.750000
B=0.0  with session wf0 cleared: e0=-2.493577
bandwidth with a cleared session start:
  e0=-2.493577 e1=-1.250000 bandwidth=3.743577
```

**Expected**: Every energy in the sweep should match the fresh-chain / ED reference, as it does on itensor_version=3: B=0 must give -2.493577, not +1.25 (which is the ferromagnetic MAXIMUM of the same Hamiltonian — wrong sign, not merely under-converged). `bandwidth(h)` must give 2.366025 at n=4 (exact max-min from a dense eigvalsh of the same operator), not 0.0. Proof it is the retained start state and not convergence: the error survives maxm=40/nsweeps=40/noise=1e-2, and clearing `sc._session.wf0` alone recovers -2.493577 and a bandwidth of 3.743577 (= 1.25 - (-2.493577), exact for n=6).

**Suggested fix**: Clear the session's start state when the Hamiltonian it was converged against is replaced: in `pyitensor/chain.py::set_hamiltonian` (and `set_hamiltonian_mpo`) set `self.wf0 = None` alongside `self._wf0_energy = None`, and do the same for `wf0_`/`have_wf0_` in both `mpscppN/chain_session.h::set_hamiltonian_mpo` so the two backends stop diverging by luck. This does not cost legitimate warm starts: `groundstate.py::gs_energy_single` sends a deliberate `set_initial_wf()` state via `self._session.set_wavefunction(wf0.cpp_handle)` at line 150, AFTER the `set_hamiltonian`/`set_hamiltonian_mpo` call at lines 139/148, so an explicitly requested warm start is re-applied on top. If cross-Hamiltonian reuse is wanted as an optimisation, it has to be opt-in and must at minimum re-randomise/perturb, since pyitensor has no noise term to escape a stationary start.

**Reviewer (CONFIRMED)**: Tried hard to refute; could not. I re-ran final_A.py, sweep.py, bw.py and fixtest.py myself (thread-pinned, PYTHONPATH forced to this checkout, dmrgpy.__file__ printed to rule out the site-packages symlink pitfall) and got the hunter's numbers verbatim: itensor_version="python", N=6, sweep on ONE chain gives B=1.0 -> -1.750000, B=0.5 -> -0.250000, B=0.0 -> +1.250000, against fresh-chain/ED -3.001995 / -2.501995 / -2.493577; itensor_version=3 on the identical script tracks ED at every point. NOT convergence error: +1.25 is the FERROMAGNETIC MAXIMUM of the 6-site Heisenberg chain (5 bonds x 0.25), i.e. the all-down state retained from the B=4 solve, which is an exact eigenstate of the B=0 Hamiltonian and therefore a stationary point no variational sweep can leave; the error survived maxm=40/nsweeps=40/noise=1e-2 (I ran that exact configuration). Mechanism verified in source, not just by outcome: pyitensor/chain.py:475-483 set_hamiltonian clears _solve_H_cache, _wf0_energy, _bandwidth_min and _bandwidth_max but NOT self.wf0, and gs_energy (line 495-500) then warm-starts from it (floor_dim = _max_link_dim(self.wf0); `if self.wf0 is None` never fires). The same file sets wf0=None at lines 179 and 375, so clearing it is an idiom the file already has and set_hamiltonian simply omits. manybodychain.restart() clears self.wf0 and _session_ham_cache with the comment "restart() promises a genuinely cold recalculation", and groundstate.py does re-send the terms (the send-cache key includes them), so the Python layer's cold-start contract is real and the session silently ignores it. Causal proof: clearing sc._session.wf0 alone (fixtest.py, no other change) restores -2.493577 and a bandwidth of 3.743577. Not audit #7: that finding was the computed_gs gate short-circuiting repeated gs_energy() calls with an UNCHANGED Hamiltonian, was fixed, and its mechanism is absent here (set_hamiltonian is called, computed_gs is cleared, the Hamiltonian is re-sent). Blast radius is not hypothetical: Many_Body_Chain.bandwidth (manybodychain.py:691-696) IS literally the two-set_hamiltonian-on-one-clone pattern, and on my run it returned -0.0 on "python" against an exact 2.366025 (v2 and v3 both returned 2.366025) — the same-clone probe showed gs_energy(-h) coming back as +1.616025, the MAXIMUM, i.e. stuck at the previous solve's state. Caveat on the write-up's scope, not its verdict: I verified the C++ side only by outcome (v3 tracks ED on the same sweep), not by reading mpscpp3/chain_session.h, so the claim that v3 retains wf0_ identically and escapes "by luck via noise" is untested by me; the finding stands without it.

**Independently reproduced by the orchestrator** (6-site Heisenberg + uniform field `B*Sz`, `itensor_version="python"`, one chain object reused across the sweep vs. a fresh chain per point vs. ED):

```
  B     sweep(reused chain)     fresh chain          ED
  0.0       -2.4935771339      -2.4935771339      -2.4935771339
  1.0       -2.4935771339      -3.0019953569      -3.0019953569
  2.0       -2.4935771339      -4.7500000000      -4.7500000000
  3.0       -7.7500000000      -7.7500000000      -7.7500000000
```

The reused chain returns the B=0 energy for B=1 and B=2 -- variationally *above* the true
ground state, so it is a trapped start rather than a stale cached number; at B=3 the
polarized state is reachable from the trapped one and the coincidence hides the bug.


---


### 3. pyitensor's sequential VUMPS returns silently wrong energies (converged=True, 29% below the exact variational minimum) or raises, at any maxm above the state's own bond dimension

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `recent-commits`

**Where**: `src/dmrgpy/pyitensor/vumps_ms.py:262-298 (_cell_fixed_points), :806-833 (ground_state's restart/safety-net loop); reached from src/dmrgpy/pyitensor/vumps.py:892 (_multisite_ground_state) and :979-1001 (the n_uc>2 / reach>1 dispatch)`

`_cell_fixed_points` gets the cell transfer matrix's dominant fixed point with a bare `scipy.sparse.linalg.eigs(op, k=1, which='LM')`. When the requested bond dimension exceeds what the state actually needs (a field-polarized chain needs D=1), the extra directions carry no weight and the transfer matrix acquires a decoupled unimodular block, i.e. a genuinely DEGENERATE dominant eigenvalue. Two things then go wrong: (1) ARPACK fails to converge and the exception is swallowed by `ground_state`'s `attempt()`, so every restart at that rung dies and the driver raises "every attempt at D=... failed"; (2) when it does return, the eigenvector is an arbitrary element of the degenerate subspace whose trace can be ~0, and the normalization `r_AL = r_AL/tr_r` (guarded only by `abs(tr_r) > 1e-300`, which is no guard at all) blows the environment up, so `e_L`/`e_R` come back as garbage. `ground_state`'s "variational safety net" only retries when `local["e_cell"] > best_e + 1e-6`, so an energy *below* the previous rung sails through and wins `better()` (which prefers converged, then LOWEST e_cell). The result is reported with `converged=True`. Crucially the STATE is fine -- `<Sz>` still reads exactly 0.5 -- only the energy read-off is broken, so nothing downstream flags it. This is the exact failure class a2eb46e (in-window) fixed on `itensor_version=3` with `Chain::vx_bond_fixed_points` (fall back to the fixed points the state itself names, C C^dag / C^dag C); the pure-Python side never got it. 71ba8eb (in-window) made this the MANDATORY route for any Hamiltonian with reach>1 at any n_uc, and it was already the only route for n_uc>2. `"python"` is the constructor default for `Infinite_Many_Body_Chain` (infinitechain.py:298) and the pip-install default backend, so this is the out-of-the-box path. `tests/test_vumps_redundant_bond_dimension.py::test_sequential_solver_tolerates_redundant_bond_dimension` exercises the sequential path at D=2 only, which passes.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/recent-commits && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 vumps_conv.py   # and vumps_freq.py, vumps_nuc3.py, vumps_diag.py in the same directory
```

Observed:

```
# vumps_conv.py -- default backend ("python"), 3-site cell, reach-1, maxm=6.
# Exact: fully polarized product state, e/site = -FIELD/2 + J/4 = -1.825.
run 0: e0=-2.3473307291666665  converged=True  exact=-1.825  (delta=-0.522)
   <Sz> = (0.5+0j)

# vumps_freq.py -- 1-site cell, reach-2 coupling (forced onto the sequential
# solver by 71ba8eb's dispatch), maxm=8, 10 independent runs per backend:
backend=python  exact=-1.825  9/10 wrong-or-raised in 62s
    ['RuntimeError', 'RuntimeError', -110365641.40041503, 'RuntimeError', 'RuntimeError', -1.825, 'RuntimeError', 'RuntimeError', 'RuntimeError', 'RuntimeError']
backend=3       exact=-1.825  0/10 wrong-or-raised in 2s
    [-1.825, -1.825, -1.825, -1.825, -1.825, -1.825, -1.825, -1.82499992, -1.825, -1.825]

# vumps_nuc3.py -- plain reach-1 Hamiltonian on a 3-site cell (no long-range
# term needed; n_uc>2 alone routes to the sequential solver), 8 runs:
backend=python  n_uc=3 reach=1 D=6 (10s): [-1.82926432, -110671909.8125, -1.825, -1.825, -1.825, -1.825, -1.825, 'RuntimeError']
backend=python  n_uc=3 reach=1 D=8 (42s): [-1.825, -1.825, 'RuntimeError', -1.825, -1.825, 'RuntimeError', 'RuntimeError', -1.825]
backend=3       n_uc=3 reach=1 D=6 (18s): [-1.825, -1.825, -1.825, -1.825, -1.825, -1.825, -1.825, -1.825]
backend=3       n_uc=3 reach=1 D=8 (65s): [-1.825, -1.825, -1.825, -1.825, -1.825, -1.825, -1.825, -1.825]

# vumps_diag.py -- the swallowed exception behind every 'RuntimeError' above:
  single_run raised: ArpackNoConvergence ARPACK error -1: No convergence (5001 iterations, 0/1 eigenvectors converged)
  single_run raised: ArpackNoConvergence ARPACK error -1: No convergence (5001 iterations, 0/1 eigenvectors converged)
  single_run raised: ArpackNoConvergence ARPACK error -1: No convergence (5001 iterations, 0/1 eigenvectors converged)
  single_run raised: ArpackNoConvergence ARPACK error -1: No convergence (5001 iterations, 0/1 eigenvectors converged)
FINAL RuntimeError vumps_ms.ground_state: every attempt at D=8 failed -- try increasing nrestarts
```

**Expected**: The model is exactly solvable: H = -4 sum_i Sz_i + 0.7 sum_i Sz_i Sz_{i+r} has the fully polarized product state as its ground state, so e/site = -4/2 + 0.7/4 = -1.825 exactly, at every maxm >= 1. Values of -2.347, -1.82926 and -1.1e8 are all BELOW that, which the variational principle forbids for a normalized state -- so they are not 'less converged', they are wrong. `itensor_version=3` returns -1.825 in 26/26 runs across the same (n_uc, reach, D) grid, which is what the pure-Python backend should do too; `tests/test_vumps_redundant_bond_dimension.py`'s own module docstring states this is fixed on both backends.

**Suggested fix**: Port `Chain::vx_bond_fixed_points`'s fallback into `_cell_fixed_points`: when `eigs` raises (ArpackNoConvergence) or the dominant eigenvalue is degenerate / the returned eigenvector's trace is not safely nonzero, use the fixed points the state itself names -- `C C^dag` for the AL transfer's right fixed point and `C^dag C` for the AR transfer's left one -- which are exact under redundancy and give the same branch mixture in a genuine cat state. Replace the `abs(tr_r) > 1e-300` normalization guard with a real tolerance relative to `norm(r_AL)` and raise rather than divide when it fails. Separately, extend `ground_state`'s safety net to reject an `e_cell` implausibly far BELOW the previous rung (not just above it), so a broken environment cannot win `better()` with `converged=True`. Extend `test_sequential_solver_tolerates_redundant_bond_dimension` from D=2 to D in {2,4,6,8} and assert the value, not just the absence of an exception.

**Reviewer (CONFIRMED)**: I tried three ways to refute this and all failed.

(a) Re-ran the hunter's vumps_conv.py verbatim (thread-pinned, PYTHONPATH into the repo src): 'run 4: e0=-5929246569643.829 converged=True exact=-1.825 (delta=-5.93e+12), <Sz> = (0.5+0j)'. Worse than the pasted -2.347, same shape.

(b) Convergence-error explanation: impossible. H = -4 sum Sz + 0.7 sum Sz_i Sz_{i+1} is diagonal in the Sz product basis (classical Ising + field); the minimum over all configurations is the polarized one at -2 + 0.175 = -1.825 exactly, so every returned value below that is forbidden for a normalized state, not 'less converged'. <Sz> comes back exactly 0.5, so the state is right and only the energy read-off is broken, as claimed.

(c) Tuned-down-parameters explanation: refuted by my own independent script (scratchpad/verify-recent-commits/v1.py) at the library DEFAULT vumps_nrestarts=4, 3-site cell, reach-1, D=6, 12 fresh runs: 5/12 bad -- e0 = -1.8275009, -1.9825033, -1.9007161, -5.7458333, -1.8250327, every one converged=True, <Sz>=0.5, every one below -1.825. At D=8 (8 runs) I also reproduced the other half: 'run 0 RAISED RuntimeError: vumps_ms.ground_state: every attempt at D=8 failed -- try increasing nrestarts' plus e0=-6.45e9 converged=True.

Mechanism check (diag.py, monkeypatching vumps_ms._cell_fixed_points and logging on a bad run): the returned r_AL has unit trace by construction but max|r_AL| = 6.4e14 with cond = 9.6e16 on the run that produced e=-7.04e10 -- i.e. the eigensolver returned a near-traceless element of the degenerate subspace and the `r_AL = r_AL/tr_r` normalization blew it up, exactly the mechanism described. (One bad-run sample; the script breaks at the first one.)

Scope: not a known issue -- the opposite is asserted in-window. docs/documentation.md:1527 says itensor_version="python", 'whose vumps_ms._cell_fixed_points has no such guard at all -- returned the exact energy', and tests/test_vumps_redundant_bond_dimension.py only exercises the sequential python path at D=2 (and AKLT at D=3,4), whose own docstring concedes vumps_ms._cell_fixed_points 'makes the arbitrary pick'. Reached out of the box: "python" is the constructor default and the pip-install default backend, gs_method="vumps" is the default, and in-window 71ba8eb made this the mandatory route for any reach>1 Hamiltonian. Code reading confirms ground_state's safety net only fires on e_cell ABOVE best_e + 1e-6 and better() prefers converged-then-lowest, so a blown-up energy wins.


---


### 4. ED boson occupation-number projectors are |k><0|+|k><1|+... , not |k><k|: bc.D[i][k] under mode="ED" returns negative "probabilities" that do not sum to 1

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `ed-and-operators`

**Where**: `src/dmrgpy/pyboson/boson.py:40 (also 36-42); surfaced through src/dmrgpy/bosonchain.py:26 (Bosonic_Chain.D) and :31-34 (D0..D3)`

`BosonChain.create_operators` builds the occupation-number projectors with

    ops = ids[i]*0.0          # a (d,d) dense numpy array
    op  = ops*0.0
    op[n] = 1.0               # <-- sets the whole ROW n, not the element (n,n)

On a 2-D numpy array `op[n] = 1.0` assigns the entire row, so the single-site operator stored under the name "N"+str(n) is sum_m |n><m| rather than the projector |n><n|. It is not Hermitian, it is not idempotent, and it changes the boson number. Every DMRG backend builds the same name correctly as a genuine projector (`pyitensor/sites/boson.py::_boson_ops` does `build_matrix(dim,[(k+1,k+1,1.0)])`; mpscpp3 uses ITensor's own BosonFourSite), so ED and DMRG disagree.

Why it has stayed invisible: every off-diagonal piece |n><m| (m!=n) changes the total boson number, and a boson-number-conserving Hamiltonian has a ground state of definite total N, so all of that garbage has exactly zero expectation value. The first run below (a number-conserving hopping+interaction chain) shows ED and DMRG agreeing to 1e-9 on all three projectors — which is precisely why no existing test catches it. Adding a single number-breaking drive term 0.5*(A_i + Adag_i) exposes it immediately.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/ed-and-operators/boson_public2.py
# and, for the bare single-site matrix:
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/ed-and-operators/boson_proj.py
```

Observed:

```
boson_public2.py  (3 sites, maxnb=[3,3,3], hopping + N^2 + 0.5*(A+Adag) drive):
E0 ED  = -5.248282665570111
E0 DMRG= -5.248282665570105
P(n=0) site0  ED=-0.0081902253  DMRG=0.1176555998  diff=1.26e-01
P(n=1) site0  ED=-0.0273662843  DMRG=0.5442488879  diff=5.72e-01
P(n=2) site0  ED=0.0938331174  DMRG=0.3380955124  diff=2.44e-01
sum_k P   ED=0.0582766077  DMRG=1.0000000000
<N>_0     ED=1.2204399126  DMRG=1.2204399126
sum k*P_k ED=0.1602999504  DMRG=1.2204399126

boson_proj.py  (single-site N0 embedded on a [3,3] chain, printed as the 9x9 many-body matrix):
N0 on site 0 ... hermitian? False  trace= (3+0j)
N1 on site 0 ... hermitian? False  trace= (3+0j)
N2 on site 0 ... hermitian? False  trace= (3+0j)
N0 rank/nnz: 9
[[1 0 0 1 0 0 1 0 0]
 [0 1 0 0 1 0 0 1 0]
 [0 0 1 0 0 1 0 0 1]
 [0 0 0 0 0 0 0 0 0]
 ... ]

(For comparison, the same script with the drive term removed gives ED==DMRG to 1e-9 on every P(n) and sum_k P = 1.0000000000022 -- the number-conserving case that hides the bug.)
```

**Expected**: bc.D[i][k] is documented in docs/user_guide.md:48 as the occupation projector $\hat n_i^{(k)}=\lvert k\rangle\langle k\rvert$. Its expectation value must be a probability: in [0,1], summing to 1 over k, with sum_k k*P_k == <N_i>. The DMRG column satisfies all three (sum = 1.0000000000, sum k*P_k = 1.2204399126 = <N>_0 exactly); the ED column satisfies none (two negative entries, sum 0.058, sum k*P_k = 0.160 against <N>_0 = 1.220). The single-site matrix must be diag(0,..,1,..,0), and the printed N0 shows it is instead the all-ones row 0.

**Suggested fix**: src/dmrgpy/pyboson/boson.py:40 -- `op[n] = 1.0` should be `op[n,n] = 1.0`. Add a regression test asserting sum_k <D[i][k]> == 1 and sum_k k*<D[i][k]> == <N_i> on a chain whose Hamiltonian is NOT boson-number conserving (a number-conserving one passes even with the bug in place).

**Reviewer (CONFIRMED)**: Re-ran both repros thread-pinned; output matches the hunter's to the digit.

boson_proj.py: N0 on a [3,3] chain prints as the all-ones row-0 pattern ([[1,0,0,1,0,0,1,0,0],[0,1,0,0,1,0,0,1,0],[0,0,1,0,0,1,0,0,1],0...]), hermitian?=False, trace=3 for every one of N0/N1/N2. boson_public2.py: E0 ED=-5.248282665570111 vs DMRG=-5.248282665570109 (so the Hamiltonian itself is fine), but P(n=0)=-0.0081902253 / P(n=1)=-0.0273662843 / P(n=2)=0.0938331174 under ED against 0.1177/0.5442/0.3381 under DMRG; sum_k P = 0.0583 (ED) vs 1.0000000000 (DMRG); sum_k k*P_k = 0.1603 (ED) vs 1.2204399126 = <N>_0 (DMRG).

Attempts to explain it away, all failed: (a) not a convergence artifact -- the defect is visible in the bare single-site matrix before any solver runs, and E0 agrees to 2e-15; (b) not a one2many transposition -- the printed 9x9 is unambiguously sum_m |0><m| on site 0, not diag; (c) not intended design -- docs/user_guide.md:48 states bc.D[i][k] is the projector |k><k|, pyitensor/sites/boson.py::_boson_ops builds build_matrix(dim,[(k+1,k+1,1.0)]) i.e. a genuine diagonal projector, and mpscpp3/get_sites.h routes to ITensor's BosonFourSite with MaxOcc; (d) not out of scope -- src/dmrgpy/pyboson/boson.py is dmrgpy's own ED code, and the item is not in docs/audit_2026_08_hole_hunt.md, any docs/known_issue_*.md, or CLAUDE.md's list of deliberately-reproduced legacy bugs.

Root cause confirmed by reading src/dmrgpy/pyboson/boson.py:36-42: `ops = ids[i]*0.0` is a (d,d) array, so `op[n] = 1.0` assigns the whole row n instead of the element (n,n). grep over tests/ and examples/ finds nothing that exercises bc.D at all, which is consistent with the bug surviving.


---


### 5. The same dynamical-correlator submode name computes different quantities under mode="ED" and mode="DMRG": submode="ED" silently drops Im of the Lehmann weight, and DMRG/CVM + ROOTN + the docstrings sit off the convention the default KPM path and ED/INV/CVM actually compute

> Retitled by the orchestrator: the reviewer confirmed the observation but disputed the
> hunter's diagnosis in the original title, which was: *mode="ED" computes a different dynamical correlator from the one documented, and from mode="DMRG": submode="ED" drops Im of the Lehmann weight, submode="CVM"/"INV" keeps the complex weight instead of -Im G/pi; three mutually inconsistent conventions, all coinciding only when A=B^dag*

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `ed-and-operators`

**Where**: `src/dmrgpy/edtk/dynamics.py:273 (submode="ED"), :356 (the T>0 sum, same shape), :409 (submode="INV"/"CVM"), :159-177 (submode="ROOTN"); src/dmrgpy/cvm.py:239 (mode="DMRG", submode="CVM"); definition in docs/user_guide.md:1148 and :1155-1156`

docs/user_guide.md:1148 defines one quantity and says every submode computes it:

    G_AB(w) = <GS| A (w - H + E0 + i*delta)^{-1} B |GS>,   S_AB(w) = -(1/pi) Im G_AB(w)
    "All of the submodes below compute G_AB(omega) (or equivalently S_AB(omega))"

With M_n = <0|A|n><n|B|0> and D = (w-Dn)^2+delta^2, the documented answer is
    S_AB = sum_n [ Re(M_n)*delta - Im(M_n)*(w-Dn) ] / (pi*D).
Three different things are actually returned:

* `edtk/dynamics.py:273` (`submode="ED"`, the "exact Lehmann sum" everything else is validated against) returns `-out.imag/(2*pi)` of `dynamical_sum`, which is sum_n Re(M_n)*delta/(pi*D) -- the *real part only*. The Im(M_n) dispersive term is dropped.
* `edtk/dynamics.py:409` (`submode="INV"` and `"CVM"` under mode="ED") returns `1j*(g1-g2)/2/pi` = sum_n M_n*delta/(pi*D) -- the full *complex* delta-weight, a third quantity again.
* `edtk/dynamics.py:159-177` (`submode="ROOTN"`) and `cvm.py:239` (`mode="DMRG", submode="CVM"`, `-G.imag/np.pi`) return the documented -Im G/pi.

All three coincide exactly when M_n is real -- i.e. when B = A^dag, or when the Hamiltonian is real so the eigenvectors can be chosen real. That is every case in tests/ (test_tompo_correlator_operators.py:112/118 use `Amo, Amo.get_dagger()`; :142 uses `Sz[0],Sz[1]` on a real H), which is why nothing catches it. They diverge for the ordinary off-diagonal Green's function <c^dag_i ... c_j>, i != j, on a chain with complex hoppings -- exactly what `fermionchain.get_gr` and any momentum-resolved A(k,w) is built out of.

Note the rootn docstring (edtk/dynamics.py:169-172) explicitly claims it "Follows the same ... convention ... as dynamical_correlator_ED/dynamical_correlator_inv above, so results are directly comparable to submode='ED'/'CVM'/'INV'". Measured, that claim is false in both directions.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/ed-and-operators/probe8.py 2>&1 | grep -v '^CVM in E'
# 4-site Fermionic_Chain, random complex Hermitian hopping + NN interaction,
# A = Cdag[0], B = C[2], es in [-1,6], delta=0.15
```

Observed:

```
exact weights: max|Re|=0.233561 max|Im|=0.573439
  cand Sigma M_n L (complex)      peak=0.619179
  cand Re[Sigma M_n L]            peak=0.233561
  cand -(1/pi)Im G^R (retarded)   peak=0.425685
  ED  /ED     max|Im|=0.0000 -> Re[Sigma M_n L]:1.67e-16  -(1/pi)Im G^R (retarded):2.88e-01  Sigma M_n L (complex):5.73e-01
  ED  /INV    max|Im|=0.5734 -> Sigma M_n L (complex):6.87e-16  Re[Sigma M_n L]:5.73e-01  -(1/pi)Im G^R (retarded):5.75e-01
  ED  /CVM    max|Im|=0.5734 -> Sigma M_n L (complex):6.94e-16  Re[Sigma M_n L]:5.73e-01  -(1/pi)Im G^R (retarded):5.75e-01
  ED  /ROOTN  max|Im|=0.0000 -> -(1/pi)Im G^R (retarded):6.11e-16  Re[Sigma M_n L]:2.88e-01  Sigma M_n L (complex):5.75e-01
  DMRG/CVM    max|Im|=0.0000 -> -(1/pi)Im G^R (retarded):6.05e-07  Re[Sigma M_n L]:2.88e-01  Sigma M_n L (complex):5.75e-01
max|CVM(ED) - CVM(DMRG)| = 0.5746   (peak of the correlator = 0.4257)
max|submodeED(ED) - CVM(DMRG)| = 0.2879

(Independently reproduced on a 4-site spin chain with a DM term, A=S+_0, B=S+_2:
  ED  /ED     -> Re(AB):2.91e-15   retarded(AB):9.33e-02
  ED  /INV    -> AB complex:1.51e-15
  ED  /ROOTN  -> retarded(AB):1.58e-15
  DMRG/CVM    -> retarded(AB):3.21e-06 )
```

**Expected**: Per docs/user_guide.md:1148 every submode must return S_AB(w) = -(1/pi) Im G_AB(w), which is the `-(1/pi)Im G^R (retarded)` column. ROOTN and mode="DMRG" submode="CVM" match it to 6e-16 and 6e-7. The ED backend does not: `submode="ED"` is off by up to 0.288 on a peak of 0.426 (68%) and `mode="ED", submode="CVM"` is off by 0.575 (135%) -- and the same submode name "CVM" returns two different quantities depending on whether mode is "ED" or "DMRG", differing by 0.5746, more than the correlator's own peak. Since `mode="ED", submode="ED"` is what tests/ and the docs both call the exact reference, any future cross-backend check of a complex-weight correlator compares two different observables.

**Suggested fix**: Pick the documented convention (-(1/pi)Im G_AB) and make the ED paths produce it. In `dynamical_correlator_ED` (edtk/dynamics.py:207-273) and `dynamical_sum`, accumulate the full retarded resolvent instead of the +i*delta minus -i*delta difference: replace the kernel with `M_n/(w+1j*delta-Dn)` and return `-out.imag/np.pi`. `dynamical_correlator_finite_T`/`dynamical_sum_thermal` (:301-385) have the identical structure and need the same change. In `dynamical_correlator_inv` (:388-409), return `-o.imag/np.pi` rather than the complex `o/np.pi` (the `mode="cv"` branch already does two solves at +-delta, so both parts are available). Fix the false comparability claim in the rootn docstring (:169-172). Also add a regression test that uses an asymmetric pair (A=Cdag_0, B=C_2) on a Hamiltonian with complex hoppings -- with B=A^dag, as every current test does, all three conventions agree and the bug is unobservable.

**Reviewer (CONFIRMED)**: The INCONSISTENCY reproduces exactly; the hunter's diagnosis of which side is wrong does not survive, so read the fix direction with care.

Re-ran probe8.py thread-pinned, identical numbers: ED/ED matches Re[Sigma M_n L] to 1.67e-16; ED/INV and ED/CVM match the complex Sigma M_n L to 6.9e-16; ED/ROOTN and DMRG/CVM match -(1/pi)Im G^R to 6.1e-16 and 6.6e-11. max|CVM(ED)-CVM(DMRG)| = 0.5746 against a correlator peak of 0.4257. Algebra checked by hand: dynamical_sum's 1/(w+i.delta-D) - 1/(w-i.delta-D) = -2i.delta/D, so -out.imag/(2pi) = Re(M_n).delta/(pi.D), while cvm.py:239's -G.imag/pi = [Re(M).delta - Im(M)(w-D)]/(pi.D). Same-submode-name cross-backend divergence is therefore real, and algebra/rootn.py:86's docstring claim that ROOTN is "directly comparable to submode='ED'/'CVM'/'INV'" is measurably false in both directions.

Two corrections to the framing, both from checks I ran myself:

1. It is TWO conventions, not three. Re[ED/INV] == ED/ED to 1e-16 (visible in the probe8 table itself): submode="ED" is just INV/CVM with .imag taken of a complex accumulator.

2. The hunter's "expected" column is the wrong target. I ran a kernel-independent sum rule (es in [-12,18], 3000 pts, delta=0.05): <A B> = 0.11052176 - 0.27135291j from ED, and integral of y gives DMRG/KPM = 0.11052123 - 0.27135161j (exact), ED/INV = 0.11028 - 0.27077j (window truncation only), ED/ED = 0.11028 (real part only), ED/ROOTN = 0.13145 (principal-value tails outside the window). So the DEFAULT submode KPM computes the complex Lehmann density Sigma M_n delta(w-Dn), i.e. i(G^R-G^A)/2pi, which is what ED/INV/CVM already return; -(1/pi)Im G^R (DMRG/CVM, ROOTN) equals that only for real M_n and is the Hermitian-pair specialization the docs wrote as if general. Adopting the suggested fix would move the default KPM path and INV/CVM away from what they correctly compute. The defect is that submode="ED" silently drops Im(M_n) and that ROOTN/DMRG-CVM/the docs disagree with the rest -- the resolution needs a convention decision plus a docs correction, not the mechanical -out.imag -> -o.imag/pi change proposed.

Ruled out as out of scope: not in docs/audit_2026_08_hole_hunt.md (#1 is the non-Hermitian submode dispatch, #10 is _fourier_transform_correlator's missing 1/pi on submode TD/TDZ -- both distinct), not a known_issue_*.md item, not one of CLAUDE.md's deliberately-reproduced legacy bugs.

**Reviewer on severity**: Keep high, but for a narrower reason than stated: the defensible high-severity core is (a) the same submode name "CVM" returning two quantities differing by more than the correlator peak across mode=, (b) submode="ED" -- the documented exact reference -- dropping Im(M_n) that the default KPM path keeps, and (c) a docstring asserting comparability that is false. The claim that mode="DMRG" submode="CVM" is the correct target is not supported; by the sum rule it is DMRG/CVM and ROOTN that sit off the house convention.


---


### 6. vev(op, npow=n) silently ignores npow on every ED route, returning <op> instead of <op^n>

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `docs-examples-drift` &middot; independently found by `python-backend-parity`, `dispatch-matrix`, `docs-examples-drift`

**Where**: `src/dmrgpy/edtk/edchain.py:226 (EDchain.vev(self,op,T=0.,**kwargs)); dispatched from src/dmrgpy/manybodychain.py:762-764; documented at docs/user_guide.md:558 and docs/documentation.md:3562`

`Many_Body_Chain.vev(MO, npow=n)` is documented as computing <X^n> by repeated MPO application. On the DMRG path (`vev.py::multi_vev`) `npow` is forwarded into `self._session.vev(..., npow=int(npow))` and honoured. On the ED path the call is `self.get_ED_obj().vev(MOf, **kwargs)`, and `EDchain.vev` has signature `(self, op, T=0., **kwargs)` — at T=0 it computes `algebra.braket_wAw(wf0, MO2matrix(op))` and never looks at `kwargs`, so `npow` is swallowed and the first moment is returned for every n. Confirmed that this is the single site for both statistics: `Spin_Chain.get_ED_obj()` returns a `pychain.build.Spin_chain` and `Fermionic_Chain.get_ED_obj()` an `pyfermion.mbfermion.MBFermion`, and `inspect.getsourcefile(type(o).vev)` resolves to edchain.py:226 for both. documentation.md §4.6 explicitly records that the same kwarg raises a plain TypeError on julia_live "rather than silently misbehaving" — the ED backend is the one that silently misbehaves, and that is recorded nowhere.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 5 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/npow.py   # 4-site Heisenberg; prints vev(h,npow=1..3) under mode='DMRG' and mode='ED'
```

Observed:

```
E0 (ED)            = -1.6160254037844393
exact <H^2> = E0^2 = 2.61153810567666
DMRG vev(h,npow=1) = (-1.6160254037844388+0j)
DMRG vev(h,npow=2) = (2.6115381056766576+0j)
ED   vev(h,npow=1) = (-1.6160254037844384+0j)
ED   vev(h,npow=2) = (-1.6160254037844384+0j)   <-- should be E0^2
ED   vev(h,npow=3) = (-1.6160254037844384+0j)   <-- should be E0^3 = -4.220311921724574
```

**Expected**: `sc.vev(h, npow=2, mode="ED")` must return <H^2>. For an exact ED eigenstate that is E0^2 = 2.6115381056766600, which is exactly what the DMRG path returns (2.6115381056766576). Either value is verifiable without dmrgpy: the ED ground state is an eigenstate, so <H^n> = E0^n identically. Downstream this is not a cosmetic difference — `gs_energy_fluctuation()` builds sqrt(|<H^2>-<H>^2|) out of it, and on any ED route that becomes sqrt(|<H>-<H>^2|), a number with no physical meaning (see the companion finding).

**Suggested fix**: Give `EDchain.vev` a real `npow` parameter rather than letting it fall into `**kwargs`: at T=0, build the matrix once and apply it n times to wf0 (or use `np.linalg.matrix_power`-free repeated `op @ v`), i.e. `v = wf0; for _ in range(npow): v = op@v; return np.vdot(wf0, v)` with `npow=0` returning 1.0 to match `multi_vev`'s own early return. If a fix is not wanted immediately, the minimum is to raise on `npow != 1` here so ED matches julia_live's honest TypeError instead of returning a wrong number.

**Reviewer (CONFIRMED)**: Re-ran my own script (4-site Heisenberg, threads pinned, taskset -c 7). Output matched the hunter's to the digit: DMRG npow=1/2/3 = -1.6160254037844388 / 2.6115381056766607 / -4.220311921724573 (exactly E0, E0^2, E0^3), while ED npow=1/2/3 all returned -1.6160254037844384. I extended the check to a 4-site spinless Fermionic_Chain: E0=-2.23606797749979, E0^2=5.0, ED vev(h,npow=2)=-2.236067977499789 -- same failure, so it is not spin-specific. inspect.getsourcefile/getsourcelines on type(get_ED_obj()).vev resolves to edtk/edchain.py:226 for BOTH pychain.build.Spin_chain and pyfermion.mbfermion.MBFermion, confirming the single shared site; that signature is (self,op,T=0.,**kwargs) and the T==0 branch never reads kwargs. Tried to explain it away three ways and could not: (a) not a convergence artefact -- the ED state is an exact eigenstate and the DMRG path returns the exact power on the same chain; (b) not documented as DMRG-only -- user_guide.md:558 and user_guide.tex:642 present `sc.vev(h, npow=2)` as general, and the model-profile table two lines below says these 'take the same mode=/**kwargs as vev'; (c) not already recorded -- grep over docs/audit_2026_08_hole_hunt.md and docs/known_issue_*.md finds no npow hit at all. documentation.md:3562 does discuss npow, and it cuts the hunter's way: it records that julia_live raises 'a plain TypeError ... rather than silently misbehaving', with no mention of the ED path, which does silently misbehave.

**Independently reproduced by the orchestrator**, and one addition. On a 4-site
Heisenberg chain (default backend v3):

```
ED   npow=1 5.55e-17   vev(A*A) 0.2500   vev(A,npow=2) 5.55e-17
DMRG npow=1 -2.52e-12  vev(A*A) 0.2500   vev(A,npow=2) 0.2500
```

so the ED route returns <A> for npow=2 while the explicit product <A^2> is correct --
`npow=` is dropped, not mis-evaluated.

The addition concerns how `gs_energy_fluctuation` is reached, because the corroborating
titles below are imprecise about it. `Many_Body_Chain.gs_energy_fluctuation(self,**kwargs)`
(`manybodychain.py:1041`) accepts `**kwargs` and **never consumes them** -- it calls
`self.vev(h)` and `self.vev(h,npow=2)` with no `mode=`. So `mode="ED"` passed to it is
silently swallowed (a textbook instance of documentation.md 4.10's "a `**kwargs` with no
consumer"), and the ED route is reached only via `self.mode` or an automatic fallback:

```
mode="ED" kwarg on a v3 chain : 4.712160915387242e-08   (identical to mode="DMRG" -- kwarg dropped)
self.mode="ED", 4 sites       : 2.056103963680119       (exact answer: 0)
ns=2, automatic v3->ED fallback: 1.1456439237389602      (exact answer: 0, and nobody opted in)
```

The last line is the one that matters: the `ns<3` fallback is automatic, so a 2-site chain
on the default backend returns a garbage fluctuation with no user action at all.

<details><summary>Corroborating report from lens <code>python-backend-parity</code>: mode="ED" silently ignores vev(..., npow=n) and returns <A> for every n, so gs_energy_fluctuation(mode="ED") returns 2.30 for an exact eigenstate instead of 0</summary>

**Where**: `src/dmrgpy/edtk/edchain.py:226 (EDchain.vev(self, op, T=0., **kwargs) -- npow is swallowed by **kwargs); consumers src/dmrgpy/manybodychain.py:1041 gs_energy_fluctuation and src/dmrgpy/groundstate.py:163 (the maxde refinement loop)`

Many_Body_Chain.vev(MO, mode=..., npow=n) is documented and honoured on itensor_version 2/3/"python" (the DMRG path forwards npow into Chain::vev, which computes <A^n>). The ED path does not: manybodychain.py:751 routes to `self.get_ED_obj().vev(MOf, **kwargs)`, and EDchain.vev's signature is `vev(self, op, T=0., **kwargs)` -- `npow` lands in **kwargs and is never read. grep confirms the string "npow" appears nowhere under edtk/. So the exact reference a user cross-checks DMRG against silently returns the wrong moment, with no warning.

The public symptom is gs_energy_fluctuation(), whose whole implementation is sqrt(|<H^2> - <H>^2|). On ED that becomes sqrt(|<H> - <H>^2|), a number with no physical meaning: on a 4-site chain it returns 2.3023368536 for a state that is an exact eigenvector of H, where the correct answer is 0 (v3 gives 2.1e-8, python 3.7e-8). Since gs_energy_fluctuation is the documented way to check convergence, "my DMRG fluctuation is 1e-8 but ED says 2.3" is a maximally confusing result. groundstate.py:163's get_gs(maxde=...) refinement loop uses the same <H^2>-<H>^2 idiom.

Note this is a distinct defect from any of the 21 findings in docs/audit_2026_08_hole_hunt.md -- "npow" and "gs_energy_fluctuation" appear zero times in that file.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/python-backend-parity && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 repro_npow.py

--- repro_npow.py (n=4 S=1/2 Heisenberg + staggered 0.3*Sz, maxm=20, nsweeps=20) ---
A = sc.Sz[0]*sc.Sz[2]
print(sc.vev(A), sc.vev(A,npow=2), sc.gs_energy_fluctuation())
# A^2 = (Sz0 Sz2)^2 = (1/4)(1/4) = 1/16 exactly for spin-1/2, so <A^2> must be 0.0625
```

Observed:

```
ED       <A>=0.1232719263  <A^2>(npow=2)=0.1232719263   gs_energy_fluctuation=2.3023368536
3        <A>=0.1232719263  <A^2>(npow=2)=0.0625000000   gs_energy_fluctuation=0.0000000211
python   <A>=0.1232719263  <A^2>(npow=2)=0.0625000000   gs_energy_fluctuation=0.0000000371
exact <A^2> for spin-1/2 = 1/16 = 0.0625
```

**Expected**: vev(A, npow=2, mode="ED") must return <A^2> = 0.0625 (an exact identity for spin-1/2: (Sz_i Sz_j)^2 = 1/16), matching both DMRG backends, and gs_energy_fluctuation() on an exact eigenstate must be ~0 on every solver. An unsupported kwarg should at minimum raise rather than be silently discarded -- this is the same "dropped keyword argument" class the 2026-08 audit was chartered against.

**Suggested fix**: Give EDchain.vev an explicit `npow=1` parameter and apply it: build the operator matrix once and raise it to the power (`op = MO2matrix(op); M = numpy.linalg.matrix_power(M, npow)`) before `algebra.braket_wAw`, mirroring Chain::vev's own semantics. Add a cross-backend regression asserting vev(Sz[i]*Sz[j], npow=2) == 0.0625 on all of ED / 3 / "python", and one asserting gs_energy_fluctuation() < 1e-6 on all three.

**Reviewer (CONFIRMED)**: The code is dispositive before any run: edtk/edchain.py:226 is `def vev(self,op,T=0.,**kwargs)` and its T==0 branch is `return algebra.braket_wAw(wf0, self.MO2matrix(op))` -- npow lands in **kwargs and is never read; `grep -rn npow src/dmrgpy/edtk/` returns nothing, while vev.py:23-34 forwards `npow=int(npow)` into the session for the DMRG backends. manybodychain.py:751 routes mode=="ED" to `self.get_ED_obj().vev(MOf,**kwargs)`, so gs_energy_fluctuation's sqrt(|<H^2>-<H>^2|) degenerates to sqrt(|<H>-<H>^2|) there.

Re-ran repro_npow.py: ED <A>=0.1232719263 and <A^2>(npow=2)=0.1232719263 against the exact spin-1/2 identity (Sz_i Sz_j)^2 = 1/16 = 0.0625, which v3 and python both return exactly; gs_energy_fluctuation came out ED 2.3023368536 vs v3 5.16e-08 / python 5.99e-08 (the DMRG values differ from the hunter's pasted 2.11e-8/3.71e-8 only in iterative noise, as expected). Not convergence error -- ED is an exact diagonalization and the wrong value is exactly <A>, not an approximation of <A^2>. docs/user_guide.md:558 documents `e2 = sc.vev(h, npow=2)` and the section states these calls take the same `mode=` as vev, so the ED path is documented, not undefined behaviour. "npow" and "gs_energy_fluctuation" appear nowhere in docs/audit_2026_08_hole_hunt.md, nowhere in ROADMAP.md, and in no docs/known_issue_*.md.

</details>

<details><summary>Corroborating report from lens <code>dispatch-matrix</code>: EDchain.vev silently discards npow=, so vev(A,npow=n) returns <A> instead of <A^n> and gs_energy_fluctuation() returns a garbage number on every ED path</summary>

**Where**: `src/dmrgpy/edtk/edchain.py:226-234 (vev); src/dmrgpy/manybodychain.py:761-764 (the ED dispatch that forwards it); src/dmrgpy/manybodychain.py:1041-1046 (gs_energy_fluctuation)`

`EDchain.vev(self,op,T=0.,**kwargs)` forwards `**kwargs` only into the `T>0` branch (`thermal_vev_ex`). At T=0 — the default, and what `Many_Body_Chain.vev(...)`'s ED dispatch always reaches — `npow=` lands in `**kwargs` and is dropped on the floor; the function returns `braket_wAw(wf0, MO2matrix(op))`, i.e. <A>, whatever power was asked for. The DMRG side honours it (`vev.multi_vev` passes `npow=int(npow)` to `self._session.vev`). This is exactly §4.10's "`**kwargs` with no consumer" shape, in the one module the 2026-08 audit's fix for `kpmdmrg` did not reach. The user-visible consequence is `Many_Body_Chain.gs_energy_fluctuation()`, which computes `sqrt(|<H^2> - <H>^2|)`: on ED it evaluates `sqrt(|<H> - <H>^2|)`, a dimensionally meaningless quantity that is large and plausible-looking rather than ~0. It fires on an explicit `mode="ED"` chain AND — with no user action at all — on `itensor_version=3` below 3 sites, where `mode.py`'s documented ns<3 fallback silently routes to ED: `itensor_version=2` and `"python"` return 0.0 for the same 2-site chain while v3 returns 1.145644.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 3 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/dispatch-matrix/npow.py    # and matrix.py 2 / matrix.py 4 for the fallback row
```

Observed:

```
# npow.py, n=4 Heisenberg, itensor_version=3, same chain before/after sc.mode="ED"
DMRG  vev(Sz0,npow=1) = (-0+0j)
DMRG  vev(Sz0,npow=2) = (0.25+0j)
DMRG  gs_energy_fluctuation = 0.0
ED    vev(Sz0,npow=1) = 0j
ED    vev(Sz0,npow=2) = 0j
ED    vev(Sz0,npow=7) = 0j
ED    gs_energy_fluctuation = 2.056104

# matrix.py 2 -- 2-site chain, NO mode= set by the user; v3 falls back to ED automatically
=== N=2 backend=v2 ===
  gs_energy_fluctuation              val 0j
=== N=2 backend=v3 ===
  gs_energy_fluctuation              val (1.145644+0j)
=== N=2 backend=python ===
  gs_energy_fluctuation              val 0j
```

**Expected**: `vev(Sz0,npow=2)` must be 0.25 on every backend (Sz^2 = 1/4 identically for an S=1/2 site), not 0.0 — 0.0 is precisely <Sz0>, i.e. the npow=1 answer. `gs_energy_fluctuation()` on a converged ground state must be ~0 on every backend, as v2/python/DMRG all report; 2.056104 is exactly sqrt(|-1.616025 - (-1.616025)^2|), confirming <H^2> was replaced by <H>. Either the ED path implements npow (it trivially can: matrix-power the assembled operator, or apply it n times to wf0), or it must raise rather than answer a different question.

**Suggested fix**: Give `EDchain.vev` an explicit `npow=1` parameter and honour it at T=0 — `A = self.MO2matrix(op)`, then `braket_wAw(wf0, np.linalg.matrix_power(A,npow))` or n successive applications to the state vector; keep `**kwargs` for the thermal branch only. Failing that, raise `NotImplementedError` on `npow!=1` in the ED branch of `Many_Body_Chain.vev`, and note that this silently corrupts `gs_energy_fluctuation()` on the automatic v3/ns<3 fallback where the caller never chose ED.

**Reviewer (CONFIRMED)**: Reproduced and independently cross-checked. Source is unambiguous: edtk/edchain.py:226-234, `def vev(self,op,T=0.,**kwargs)` — the T==0 branch computes braket_wAw(wf0, MO2matrix(op)) and never touches **kwargs; only the T>0 branch forwards them to thermal_vev_ex. My run of npow.py (n=4 Heisenberg, v3): DMRG vev(Sz0,npow=2)=0.25 while the same chain with sc.mode="ED" gives 0j for npow=1, npow=2 AND npow=7 — identical for every power, which is the signature of the argument being dropped rather than of any numerical issue. The hunter's own exact-reference block in npow.py crashed (ValueError in np.vdot, a flaw in their reference script, not in the finding), so I checked analytically instead, which is stronger: for an S=1/2 site Sz^2 = 1/4 identically as an operator, so <Sz0^2> = 0.25 in ANY state — 0.0 cannot be right at any convergence. Likewise ED gs_energy_fluctuation = 2.056104 is exactly sqrt(|-1.616025 - (-1.616025)^2|) = sqrt(4.227561), i.e. <H^2> was replaced by <H>, confirming the substitution numerically. The automatic-fallback leg reproduces on my own script too (no hunter code): a 2-site Heisenberg chain with NO mode= set by the user gives gs_energy_fluctuation = 1.145644 on itensor_version=3 (mode.py prints its ns<3 -> ED fallback) versus 1.49e-08 on v2 and 1.05e-08 on "python"; 1.145644 = sqrt(|-0.75 - 0.75^2|), the same substitution. Not out of scope: not vendored ITensor, not one of CLAUDE.md's three deliberately-reproduced legacy bugs, not in docs/known_issue_*.md, and not in docs/audit_2026_08_hole_hunt.md (I grepped it for npow — no hit). The DMRG side does honour npow (vev.multi_vev passes npow=int(npow)), so this is a genuine cross-solver divergence on a documented kwarg.

</details>

<details><summary>Corroborating report from lens <code>docs-examples-drift</code>: gs_energy_fluctuation() returns a meaningless number on every ED route, including the automatic ns<3 fallback nobody opted into</summary>

**Where**: `src/dmrgpy/manybodychain.py:1041 (gs_energy_fluctuation(self,**kwargs)); documented at docs/user_guide.md:576-581`

`gs_energy_fluctuation(self,**kwargs)` accepts `**kwargs` and forwards none of it: its body is `h=self.get_hamiltonian(); e=self.vev(h); e2=self.vev(h,npow=2); return sqrt(|e2-e^2|)`. Two consequences. (i) `mode="ED"` is silently swallowed — the value returned is byte-identical to the no-kwarg DMRG call, while setting the attribute `sc.mode="ED"` (which `vev` reads through `get_mode`) gives a different answer, proving the kwarg never reached the dispatcher. (ii) Whenever the call does reach ED — explicit `sc.mode="ED"`, or one of mode.py's automatic fallbacks: `itensor_version=3` with `ns<3`, or no compiled extension for an explicitly named C++ version — the ignored `npow` (companion finding) turns the formula into sqrt(|<H>-<H>^2|). On a 2-site chain, where v3 falls back to ED with no user opt-in at all, `gs_energy_fluctuation()` returns 1.1456 for a state that is an exact eigenstate and must give ~0. The user guide documents this method as 'a measure of how sharply the DMRG/ED state is an eigenstate' with the formula sqrt(<H^2>-<H>^2), and `examples/groundstate/GS_enforce_maximum_fluctuation` uses it as a convergence criterion.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 5 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/fluct.py ; MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 5 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/fluct2.py
```

Observed:

```
# fluct.py (6-site Heisenberg, maxm=1 to make DMRG and ED differ visibly)
maxm=1 DMRG  gs_energy()          = -1.5577939788496395
maxm=1 gs_energy(mode='ED')       = -2.4935771338879236
maxm=1 gs_energy_fluctuation()          = 0.022827801908630284
maxm=1 gs_energy_fluctuation(mode='ED') = 0.022827801908630284
mode attr 'ED' gs_energy_fluctuation()  = 2.9515257167330686

# fluct2.py (2-site chain -> mode.py routes v3 to ED automatically, no kwarg passed)
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
get_mode -> ED
E0 = -0.75  exact singlet = -0.75
gs_energy_fluctuation() = 1.1456439237389602  (exact eigenstate -> should be ~0)
sqrt(|E0-E0^2|)        = 1.14564392373896
```

**Expected**: `gs_energy_fluctuation(mode="ED")` should equal the `sc.mode="ED"` result (it does not: 0.0228 vs 2.95), and on the 2-site chain the answer must be ~0 to machine precision, since the ED ground state of the 2-site Heisenberg model is an exact eigenstate at E0=-0.75. The printed 1.1456439237389602 is exactly sqrt(|E0 - E0^2|) = 1.14564392373896, which identifies the mechanism: <H^2> came back as <H>.

**Suggested fix**: Forward the kwargs: `h=self.get_hamiltonian(); e=self.vev(h,**kwargs); e2=self.vev(h,npow=2,**kwargs)`. That alone fixes (i). (ii) additionally needs the `npow` fix in `EDchain.vev`; until that lands, the ED value stays wrong even with the kwargs forwarded, so both changes belong in the same commit. A regression test asserting `gs_energy_fluctuation()<1e-8` on a 2-site chain would have caught this and costs milliseconds.

**Reviewer (CONFIRMED)**: Re-ran both halves with my own script. Part (i): on a 6-site S=1/2 chain at maxm=1, gs_energy_fluctuation() and gs_energy_fluctuation(mode='ED') returned the identical 6.524258637127171e-07, while setting sc.mode='ED' on a freshly built identical chain gave 2.9515257167330686 -- the kwarg demonstrably never reaches the dispatcher, which matches the source: the body is h=self.get_hamiltonian(); e=self.vev(h); e2=self.vev(h,npow=2), forwarding none of **kwargs. Part (ii): on a 2-site chain, get_mode() prints the ns<3 fallback banner and returns 'ED'; E0=-0.75 (exact singlet) and gs_energy_fluctuation()=1.1456439237389602, against sqrt(|E0-E0^2|)=1.14564392373896 -- agreement to 15 digits identifies the mechanism as <H^2> coming back as <H>, i.e. the F1 root cause. One discrepancy with the pasted output, which I judge immaterial: my maxm=1 DMRG fluct() was 6.5e-7 where the hunter pasted 0.0228 (maxm=1 DMRG lands on different product-state minima run to run, and that number is not the claim). Every load-bearing observation matched exactly. Not a convergence artefact (the 2-site ED state is an exact eigenstate; the correct answer is 0), not in the audit or known_issue files, and the shipped example examples/groundstate/GS_enforce_maximum_fluctuation/main.py does use this method as the y-axis of a convergence plot. F1 and F2 share one root cause and should be judged as one fix.

</details>


---


### 7. itensor_version=3's Chain::tdvp_step hardcodes DoNormalize=true, so complex-time evolution cannot decay and submode="TDZ" returns 36% too much spectral weight (peak 2.8x too high); itensor_version="python" is correct

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `python-backend-parity`

**Where**: `src/dmrgpy/mpscpp3/chain_session.h:1601-1614 (Chain::tdvp_step, the {"DoNormalize",true} arg); consumed by src/dmrgpy/tdz.py:168-175`

Chain::tdvp_step() passes {"DoNormalize",true} to ITensor's tdvp() unconditionally. For real-time evolution that is harmless (TDVP is norm-preserving anyway, and quench_tdvp/evolve_and_measure_tdvp restore the norm explicitly afterwards). For COMPLEX time -- the one thing this method exists to support, per its own 20-line comment ("dt may be any complex number ... complex time evolution (TDZ) share this same code path unchanged") -- the norm genuinely decays, and forcing it to 1 every step destroys exactly the damping the TDZ contour is built on. pyitensor's tdvp_step does not normalize, so "python" gets it right.

Measured directly: one complex-dz step on a state of norm^2 = 0.25 returns norm^2 = 1.0000000000 on v3 and 0.2109991178 on "python". Downstream, the exact sum rule int S(w) dw = <Sz0 Sz0> = 0.25 is satisfied by KPM/EX/TD/SECTOR and by python's TDZ (0.2495), while v3's TDZ returns 0.3407 with a peak of 0.9287 against the correct 0.3295. Both v3 TDZ integrators are affected (tevol_method="TDVP" and "TDVP_GSE" give the identical wrong number, since tdz.py:170 routes through the same tdvp_step with num_center=1).

This relocates two things already on record. docs/audit_2026_08_hole_hunt.md finding #10 states "TDZ's own residual factor of 2 lives in its reconstruction and is NOT addressed here" -- but tdz.py's reconstruction is shared Python code and "python" gets the right answer with it, so the reconstruction is not the culprit. And finding #9's "Where a fix goes" explicitly advises "not by flipping tdvp_step's DoNormalize, since tdz.py's complex-time path and metts_vev depend on its current semantics" -- the tdz.py half of that rationale is exactly backwards: tdz.py is broken BY the current semantics. (The metts_vev half stands: METTS normalizes after each beta/2 imaginary-time step by design, so that caller does want it.)

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/python-backend-parity && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 repro_tdz_norm.py && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 repro_tdz_sumrule.py

--- repro_tdz_norm.py (n=6 Heisenberg+0.3*Sz, maxm=30) ---
dz = 0.2*np.exp(-1j*0.5)                      # a TDZ-style complex contour increment
wf = sc.toMPO(sc.Sz[0])*sc.wf0                # norm^2 = 0.25
h  = sc._session.tdvp_step(Hop.cpp_handle, wf.cpp_handle, dz)

--- repro_tdz_sumrule.py (n=6, es=linspace(-20,20,4001), delta=0.1, maxm=20, nsweeps=12) ---
# exact sum rule: int S(w) dw = <Sz0 Sz0> = 0.25
```

Observed:

```
itensor_version=3         |psi|^2 before=0.2500000000  after one complex-dz tdvp_step=1.0000000000
itensor_version=python    |psi|^2 before=0.2500000000  after one complex-dz tdvp_step=0.2109991178

exact sum rule int S(w)dw = <Sz0 Sz0> = 0.25000000000000006
itensor_version=3        KPM   integral=0.250000  peak=0.7473 at w=0.490
itensor_version=3        EX    integral=0.246602  peak=0.3616 at w=0.490
itensor_version=3        TD    integral=0.249474  peak=0.3277 at w=0.520
itensor_version=3        TDZ   integral=0.340676  peak=0.9287 at w=0.520
itensor_version=python   KPM   integral=0.250000  peak=0.7473 at w=0.490
itensor_version=python   EX    integral=0.246602  peak=0.3616 at w=0.490
itensor_version=python   TD    integral=0.249495  peak=0.3300 at w=0.520
itensor_version=python   TDZ   integral=0.249495  peak=0.3295 at w=0.520

(widening: both v3 TDZ integrators are equally wrong)
=== F3 widening: TDZ with tevol_method=TDVP_GSE (exact sum rule 0.25)
  v=3        TDVP      integral=0.340676 peak=0.9287
  v=3        TDVP_GSE  integral=0.340676 peak=0.9287
  v=python   TDVP      integral=0.249495 peak=0.3295
  v=python   TDVP_GSE  integral=0.249473 peak=0.3277
```

**Expected**: ROADMAP.md:96 lists TDZ as fully implemented on all three backends, and docs/user_guide.md fixes the codebase-wide convention that every submode returns the same S_AB(omega). A complex-time TDVP step must preserve the norm decay that the contour introduces -- that decay IS the damping the TDZ reconstruction inverts. Expected: v3's tdvp_step on a norm^2=0.25 state with a complex dz returns something below 0.25 (python returns 0.2109991178), and v3's TDZ satisfies the sum rule at 0.2495 with a peak of ~0.33, agreeing with its own TD/KPM/EX and with python's TDZ.

**Suggested fix**: Do not simply flip DoNormalize to false -- metts_vev()'s imaginary-time steps rely on the current semantics. Either (a) add a `bool normalize=true` parameter to Chain::tdvp_step and have bindings.cc / tdz.py's two call sites pass false, or (b) set DoNormalize=false in tdvp_step and add the explicit psi.normalize() to the callers that want it (quench_tdvp and evolve_and_measure_tdvp already restore norm0 explicitly; metts_vev would need one added). Then extend tests/test_dynamical_correlator.py with a sum-rule assertion for TDZ -- int S(w) dw == <A B> to a few percent -- on both itensor_version=3 and "python"; the existing tests only check peak position, which the factor does not move.

**Reviewer (CONFIRMED)**: Reproduced both halves: one complex-dz tdvp_step on a norm^2=0.25 state returns 1.0000000000 on v3 and 0.2109991178 on python; the n=6 sum rule gives v3 TDZ integral 0.340676 / peak 0.9287 against python TDZ 0.249495 / peak 0.3295 with the exact rule at 0.25, while KPM/EX/TD agree with python on both backends. {"DoNormalize",true} is in mpscpp3/chain_session.h (dmrgpy's own code), not in vendored mpscpp3/TDVP/, so it is not out of scope on that axis.

I then ran the decisive control the hunter did not: v3 with tevol_method="MPO", which routes tdz.py's else-branch to the NON-normalizing evolve_taylor_step instead of tdvp_step. v3 TDZ integral = 0.249508, peak 0.4807 -- identical to python's 0.249508/0.4800 and to the exact 0.25 on the same chain. Same backend, same reconstruction, same Hamiltonian; the only thing changed is whether the stepper normalizes. That pins the cause to DoNormalize and leaves nothing to explain away. (It also means tevol_method="MPO" is a working v3 TDZ workaround, narrowing the hunter's "both v3 integrators are affected" to TDVP and TDVP_GSE only.)

SCOPE, since this brushes the audit: docs/audit_2026_08_hole_hunt.md #10 does record an undiagnosed TDZ residual on the default backend and conjectures it "lives in tdz.py's reconstruction", and #9 advises against flipping DoNormalize "since tdz.py's complex-time path ... depend[s] on its current semantics". I re-ran at the audit's OWN configuration (n=4, Heisenberg + 0.3*Sz, dt/defaults, es=linspace(-20,20,4001), delta=0.1): v3 TDZ 0.340842, python TDZ 0.249508. So the audit's diagnosis is disproved by measurement, the residual is v3-only, and the tdz.py half of #9's rationale is exactly backwards (the metts_vev half stands). The audit's "factor of exactly 2" was arithmetic on its own pi-inflated number, not a measurement, so today's 1.36x is not a discrepancy needing explanation. I am not calling this out of scope: the symptom was on record, the cause and the parity split are new, and the recorded guidance is wrong in a way that would misdirect a fixer.

**Reviewer on severity**: medium-high rather than high. It is a wrong number from a documented public API on the default backend, which argues for high; against that, the symptom was already on record in audit #10 as a known-unfixed TDZ inaccuracy, so this sharpens an existing entry rather than opening a new hole, and a working workaround (tevol_method="MPO") exists on the affected backend.


---


### 8. idmrg_window._close_array_chain composes the full (chi,chi,chi,chi) transfer chain: O(n*chi^6) where propagating the boundary matrix is O(n*chi^3*d) — 97% of td_dynamical_correlator

`optimization` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `pyitensor-performance`

**Where**: `src/dmrgpy/pyitensor/idmrg_window.py:964-1022 (_close_array_chain), hot via :1077 snapshot_correlator -> :1169 dynamical_correlator_td -> :1331 dynamical_correlator_komega; called from src/dmrgpy/infinitechain.py:1342 td_dynamical_correlator`

`_close_array_chain` builds a running rank-4 transfer tensor E and grows it site by site with `E = np.einsum('lLrR,rRsS->lLsS', E, step)`. That composition is O(chi^6) per site (and, with numpy's default `optimize=False`, runs as an unoptimised C loop rather than BLAS), and it is done twice per call — once for the measured chain and once again for the `E_id` ground-state calibration denominator. But the composed E is only ever used as `close(E4) = einsum('rR,rR->', einsum('lL,lLrR->rR', l, E4), rho_after[p_right])`, i.e. it is immediately closed on both ends. Propagating the left boundary matrix `l` through the chain one site at a time instead (two `tensordot`s per site, O(chi^3 d)) computes the identical scalar without ever materialising a rank-4 object. Measured exponent of the shipped routine over maxm=6,8,12,16 is ~chi^5.3 (all-points log-log slope; top-two-point ratio chi^5.8); a standalone head-to-head of the two contraction orders gives chi^6.1 for compose against chi^3-ish for propagate. Note `mpscpp3/chain_session.h::idmrg_close_array_chain` (line 10522) is a literal port with the same six-deep loop nest, so the v3 backend shares the defect — out of this lens, but relevant to whoever fixes it. Practical impact: `td_dynamical_correlator`'s own defaults are `nt=200, maxdim=60` and `Infinite_Many_Body_Chain.maxm` defaults to 30, and at maxm=16 with nt=4 the call already takes ~20 s of which 97% is this closure.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s9_window.py 16 4 6   # cProfile of td_dynamical_correlator
# chi sweep: same script with first arg 6, 8, 12, 16
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s10_close.py   # compose vs propagate, isolated
# end-to-end, _close_array_chain replaced at runtime (repo untouched):
for p in shipped patched; do MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s16_window_patch.py 16 4 6 $p; done
```

Observed:

```
## td_dynamical_correlator maxm=16 nt=4 n_window=6 wall=19.94s
   _close_array_chain              ncalls=24     tot=0.007 cum=19.270
   snapshot_correlator             ncalls=8      tot=0.006 cum=19.277
   <built-in method numpy._core._multiarray_umath.c_einsum> ncalls=1323 tot=18.838 cum=18.838

_close_array_chain cum by maxm:  6 -> 0.101s | 8 -> 0.435s | 12 -> 3.617s | 16 -> 19.270s
  (all-points slope log(19.27/0.101)/log(16/6) = 5.35)

 chi  nsites   compose(s)   propagate(s)   speedup   |diff|
   8     12      0.00522       0.000093     56.2x   1.89e-15
  12     12      0.05768       0.000101    569.7x   1.92e-15
  16     12      0.30226       0.000130   2325.2x   1.70e-14
  24     12      3.36436       0.000224  15033.7x   4.65e-16
  32     12     25.42211       0.000466  54571.6x   2.72e-15

returned <class 'tuple'> [(200,), (800,), (200, 800)]
maxm=16 nt=4 nwin=6  shipped   15.847 s
maxm=16 nt=4 nwin=6  PATCHED   0.458 s
S(k,w) shape (200, 800) max |a| 0.1174437922930987 max rel diff 5.060539318354355e-13
(an earlier, more loaded pair of the same runs: shipped 21.312 s vs PATCHED 0.673 s, 31.7x)
```

**Expected**: Closing a chain of doubled transfer steps against a left and a right boundary matrix is O(n*chi^3*d); the rank-4 intermediate is never needed. The module already knows this shape — `idmrg.py::_apply_site_transfer`'s docstring spells out the identical re-association ('building the transfer tensor costs O(chi^4 d) and applying it another O(chi^4), whereas contracting rho in first makes both halves O(chi^3 d). Exact, not approximate -- only the contraction order changes'), and `vumps.py::_precompute_bond_environments`'s docstring makes the same argument for its own closures. The runtime-patched version returns S(k,omega) agreeing with the shipped one to 5.1e-13 relative, which is ordinary floating-point reassociation noise, so the win is free.

**Suggested fix**: Rewrite `_close_array_chain` to propagate the boundary: X = l_before[p_left % n_cell]; for each (K, B) in zip(ket_arrays, bra_arrays): tmp = np.tensordot(X, K, axes=([0],[0])); X = np.tensordot(tmp, np.conj(B), axes=([0,1],[0,1])); then close with np.einsum('rR,rR->', X, rho_after[p_right % n_cell]). Apply the same loop to the E_id calibration chain (and see finding 4 — that denominator depends only on (result, p_left, len(ket_arrays)) and can be memoized outright). The same transformation applies to `mpscpp3/chain_session.h::idmrg_close_array_chain`.

**Reviewer (CONFIRMED)**: Re-ran everything myself, core-pinned (taskset -c 8, MKL/OMP/OPENBLAS/NUMEXPR=1, PYTHONPATH=<repo>/src). Code reading first: the rank-4 E built in the loop is consumed ONLY by the local close() at the end of _close_array_chain — no other consumer — so the re-association is unconditionally available, and the ket/bra pair can have different bond dims without breaking the propagate form (X carries (l_ket,l_bra) and stays rectangular). Ran the hunter's s16_window_patch.py at maxm=16 nt=4 nwin=6: shipped 22.49 s vs PATCHED 0.71 s, and I loaded the two saved S(k,omega) arrays myself: max abs diff 2.18e-13 on max|a|=0.117, i.e. 1.9e-12 relative — floating-point reassociation noise, not a changed answer. I then wrote my own split-variant harness (verify-pyitensor-performance/v_window_split.py) that patches _close_array_chain with four independent combinations (compose/propagate x cached/uncached) so finding 1 and finding 4 are separated: compose+cache 20.83 s, propagate-only (no cache) 1.254 s, propagate+cache 0.713 s, shipped 22.56 s — so the contraction order alone is the ~20x, independently of the caching. Output of propagate+cache vs shipped: 6.1e-13 max relative. Scaling: closure-only cost (shipped minus propagate+cache) is 0.525 s at maxm=8, 3.70 s at 12, 21.8 s at 16 → log-log exponent 5.4, matching the hunter's 5.35. Isolated head-to-head (s10_close.py) reproduced: 37.8x at chi=8 up to 37704x at chi=32, |diff| ~1e-15. Scope checks came back clean: nothing in docs/audit_2026_08_hole_hunt.md, the three known_issue_*.md files, or docs/pip_install_and_pyitensor_performance_plan.md / idmrg_improvement_plan.md covers _close_array_chain's contraction order; td_dynamical_correlator is a documented public Infinite_Many_Body_Chain method (not a prototype), gated to gs_method='idmrg' but with its own defaults nt=200, maxdim=60, i.e. far past the nt=4/maxm=16 setting that already costs 22 s. One sub-claim I did NOT verify and am not inheriting: that mpscpp3/chain_session.h:10522 idmrg_close_array_chain is a literal port of the same loop nest — out of this lens, unread.


---


### 9. idmrg_excitations._op_transfer_matrix returns a transposed (non-contiguous) view, so every idmrg._apply_transfer/_apply_transfer_from_left .reshape() silently copies the whole chi^4 tensor — 41% of a VUMPS ground-state solve

`optimization` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `pyitensor-performance`

**Where**: `src/dmrgpy/pyitensor/idmrg_excitations.py:246-259 (_op_transfer_matrix, the `.transpose(0, 2, 1, 3)`); consumed by src/dmrgpy/pyitensor/idmrg.py:2018 (_apply_transfer) and :2815 (_apply_transfer_from_left); hot via src/dmrgpy/pyitensor/vumps.py:652/693 (_solve_left_environment/_solve_right_environment) and :585/628 (_energy_density_and_source_from_left/right)`

`_op_transfer_matrix` ends with `np.tensordot(ket, np.conj(bra), axes=([1],[1])).transpose(0, 2, 1, 3)`. The transpose produces a stride-permuted view that is not C-contiguous. Both consumers then do `E4.reshape(chi*chi, -1)`, and numpy cannot reshape a non-contiguous array without materialising a full copy — so every single application of that transfer tensor pays an extra chi^4 complex memcpy. In `_solve_left_environment`/`_solve_right_environment` the tensor is built ONCE and then applied inside the iterative `_solve_linear_map` many thousands of times, so the copy is paid per iteration for a tensor that never changes. cProfile attributes 3.76 s of a 9.07 s D=24 VUMPS ground-state solve to `ndarray.reshape` called from exactly those two functions (1.895 s from `_apply_transfer`, 1.869 s from `_apply_transfer_from_left`), i.e. 41% of the whole solve is memcpy. `gs_method="vumps"` is the DEFAULT for `Infinite_Many_Body_Chain`, so this is on the default path. The repo has already found and fixed this exact bug class once, in the finite-chain matvec: `kernels.py:398` says outright 'since a transposed array is generally not contiguous, that `.reshape` forces a full copy of the operator tensor each time, not a view' — `_op_transfer_matrix` is the same mistake, unfixed. `idmrg_excitations`' own public surfaces (`excitation_energies`, `dynamical_structure_factor`, `spectral_weights`) route through the same helper and very likely pay the same cost; I did not time them, so treat that extension as unverified.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s13_reshape.py 24   # cProfile, reshape attributed to its callers
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s14_stride.py   # per-application micro-benchmark
for d in 16 24 32; do for p in shipped patched; do MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s15_patch.py $d $p; done; done   # runtime np.ascontiguousarray patch, min of 3
```

Observed:

```
TOTAL 9.07 s  e0=-0.4431330329722514
FUNC ~ 0 <method 'reshape' of 'numpy.ndarray' objects> ncalls 265459 tot 3.825
    caller idmrg.py:2018(_apply_transfer)              (41115, 41115, 1.895152605, ...)
    caller idmrg.py:2815(_apply_transfer_from_left)    (40590, 40590, 1.869433572, ...)
    caller numeric.py:997(tensordot)                   (13482, 13482, 0.011392150, ...)

  D   E4 contiguous?   _apply_transfer(view)  _apply_transfer(ascontig)  copy-overhead  site-route
    8   False             0.000006 s            0.000002 s          2.9x        0.000008 s (1x)
   16   False             0.000064 s            0.000016 s          4.0x        0.000014 s (4x)
   24   False             0.000404 s            0.000082 s          4.9x        0.000027 s (15x)
   32   False             0.002990 s            0.000629 s          4.8x        0.000050 s (60x)
   48   False             0.028741 s            0.004985 s          5.8x        0.000150 s (191x)

D=16  shipped   2.096 s   e0=-0.44311587653365
D=16  PATCHED   1.564 s   e0=-0.44311587654291
D=24  shipped   9.656 s   e0=-0.44313303298424
D=24  PATCHED   5.316 s   e0=-0.44313303298469
D=32  shipped  39.751 s   e0=-0.44313904633464
D=32  PATCHED  15.276 s   e0=-0.44313904633538
```

**Expected**: A transfer tensor built once and applied N times should be laid out so the N applications are pure BLAS. `np.ascontiguousarray` alone gives 1.34x / 1.82x / 2.60x end-to-end at D=16/24/32 (the win grows with D, as it must). The e0 values agree to ~1e-11, which is within the shipped code's own run-to-run variation (two independent shipped D=24 runs already differed by 1.3e-11) — `ascontiguousarray` cannot change arithmetic, since the reshape was copying to a contiguous buffer anyway; it only stops doing so per application.

**Suggested fix**: Minimal: return `np.ascontiguousarray(np.tensordot(...).transpose(0, 2, 1, 3))` from `_op_transfer_matrix`, or equivalently build it as `np.einsum('lpr,LpR->lLrR', ket, np.conj(bra))` (which idmrg.py's own `_transfer_matrices` already does and which yields a contiguous result). Better: in `_solve_left_environment`/`_solve_right_environment`, where the operator is the plain identity transfer of AL/AR, drop `E_id` entirely and call `idmrg._apply_site_transfer(AL, None, X)` — O(D^3 d) instead of O(D^4), measured 15x at D=24 and 191x at D=48 per application, with no O(D^4) allocation at all.

**Reviewer (CONFIRMED)**: Checked the mechanism independently of the hunter's harness before touching their scripts: for D=8 and D=24, ie._op_transfer_matrix(A,A,None) gives C_CONTIGUOUS=False with strides (8192,128,1024,16) / (221184,384,9216,16), and np.shares_memory(E4.reshape(D*D,-1), E4) is False — the reshape really is a full chi^4 copy on every application. The einsum form idmrg._transfer_matrices uses ('lpr,LpR->lLrR') is C_CONTIGUOUS=True and agrees to 4e-15, so the contiguous alternative is free. Profile reproduced (s13_reshape.py 24): of a 15.21 s D=24 VUMPS solve, 2.754 s of reshape is attributed to idmrg.py:2018(_apply_transfer) and 2.715 s to idmrg.py:2815(_apply_transfer_from_left) = 36% (hunter: 41% of 9.07 s; my box is ~1.5x slower overall, the ratio is what matters). End-to-end runtime patch (s15_patch.py, min of 3): D=16 3.314 s -> 2.571 s (1.29x), D=24 13.720 s -> 7.961 s (1.72x), e0 agreeing to ~1e-11 — matching the hunter's 1.34x/1.82x. Per-application micro (s14_stride.py): view vs ascontiguous 4.0x at D=16/24, 5.3x at D=48. Scope: this is NOT already recorded. docs/idmrg_fermionic_infinite_chain_plan.md Step 5 item 1 is the einsum->tensordot rewrite of exactly this helper (the code comment confirms it shipped: 563us->98us at D=8), and my contiguity check shows the einsum form was contiguous while the tensordot+transpose form is not — so the non-contiguity is a regression introduced BY that recorded fix, and postdates the doc. The same plan's note about vumps._environments materialising E4 only to close it against I was separately addressed (_precompute_bond_environments now closes without materialising); _solve_left_environment/_solve_right_environment were not. One constraint on the suggested 'better' fix, not a refutation: only a right-action _apply_site_transfer exists (idmrg.py:2243); _solve_left_environment would need a left-action mirror written. The minimal np.ascontiguousarray fix stands on its own and is what I measured.


---


### 10. promote_to_dense() on itensor_version="python" silently re-solves the ground state unconstrained, returning the GLOBAL ground-state energy instead of the sector's -- the exact opposite of what its own docstring guarantees

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `python-backend-parity`

**Where**: `src/dmrgpy/manybodychain.py:530-601 (promote_to_dense, the docstring guarantee and the `self._session_ham_cache = None` at the end); mechanism in src/dmrgpy/pyitensor/chain.py:475-484 (set_hamiltonian clears _wf0_energy) + :485-489 (gs_energy's skip_dmrg short circuit)`

Many_Body_Chain.promote_to_dense()'s docstring states verbatim: "The Hamiltonian and the band-edge/iDMRG/VUMPS caches are dropped (they were built on the QN indices) and re-sent on the next call that needs them; the ground-state energy and wavefunction are kept, so a bare `gs_energy()` afterwards returns the sector's energy rather than re-solving unconstrained. Call `restart()` if an unconstrained re-solve is what you want."

On a Hubbard chain confined to Nf=3 (where the sector energy -1.9408140222 differs from the global one -2.3399130755, unlike the spin-Sz=0 and spinless-Nf=3 cases where they coincide and the bug is invisible), mode="ED" and itensor_version=3 honour that contract; itensor_version="python" returns the global ground state, and a subsequent vev of the particle number reads 2 instead of 3. No warning is printed.

Mechanism, confirmed by timing: promote_to_dense() sets self._session_ham_cache = None, so the next gs_energy() re-sends the Hamiltonian, which invalidates the session's energy cache on BOTH backends -- the post-promote call takes 0.063s (v3) / 0.103s (python) against 0.0000s for a genuinely cached repeat, i.e. neither backend actually honours the "no re-solve" promise. v3 is rescued by physics: its DMRG restart begins from the promoted Nf=3 state and, because H commutes with N, cannot leave that sector, so it lands back on -1.9408140222. pyitensor's solver does leave it. So the documented guarantee is unimplemented on both, and only bites as a wrong number on "python".

The blast radius is anything that goes through get_gs()/gs_energy() after promotion -- which is the documented workflow itself ("promote after the expensive sector-confined part of the calculation", then apply a charge-changing operator). Passing the saved wf0 explicitly still works (self.wf0 and self.e0 are preserved correctly on both backends, verified), so the failure is confined to implicit re-derivation.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/python-backend-parity && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 repro_promote.py

--- repro_promote.py (Spinful_Fermionic_Chain n=3, U=1.7 Hubbard, maxm=64, nsweeps=30) ---
sc.set_conserved_sector(Nf=3); e_sec = sc.gs_energy()
sc.promote_to_dense()
e_after = sc.gs_energy(); q_after = sc.vev(sum(sc.N)).real
```

Observed:

```
ED        sector Nf=3: e0=-1.9408140222 <N>=3.000000  | after promote_to_dense: e0=-1.9408140222 <N>=3.000000
3         sector Nf=3: e0=-1.9408140222 <N>=3.000000  | after promote_to_dense: e0=-1.9408140222 <N>=3.000000
python    sector Nf=3: e0=-1.9408140222 <N>=3.000000  | after promote_to_dense: e0=-2.3399130755 <N>=2.000000
(global, unconstrained, ED) e0 = -2.3399130755495743

--- mechanism: neither backend actually skips the re-solve (final2.py) ---
=== F4 mechanism: is gs_energy() after promote_to_dense a re-solve?
  3        sector e0=-1.9408140222  repeat call 0.0000s  after promote 0.0633s -> e0=-1.9408140222
  python   sector e0=-1.9408140222  repeat call 0.0000s  after promote 0.1026s -> e0=-2.3399130755
```

**Expected**: Per promote_to_dense()'s own docstring, gs_energy() immediately after promotion must return the sector's energy (-1.9408140222) on every sector-capable backend, and must not re-solve (the ED and v3 answers, and the timing of a genuinely cached call, both confirm that is the intended semantics). It should also be instant, not a 0.06-0.10s re-sweep.

**Suggested fix**: Make the guarantee real rather than accidental. In pyitensor/chain.py::gs_energy, when skip_dmrg is set and self.wf0 is not None but self._wf0_energy is None, return inner(wf0,H,wf0)/inner(wf0,wf0) instead of falling through to a full dmrg() (mpscpp3's Chain::gs_energy would want the same treatment so it stops re-sweeping too). Alternatively have Many_Body_Chain.promote_to_dense re-prime the session after the forced Hamiltonian re-send (set_wavefunction(self.wf0) plus restoring the energy cache). Add a regression on the Hubbard Nf=3 case above -- it is the only one of the three sector shapes where the sector and global energies differ, so it is the only one that can see this.

**Reviewer (CONFIRMED)**: Re-ran repro_promote.py: ED and itensor_version=3 both return e0=-1.9408140222 with <N>=3.000000 after promote_to_dense, while python returns -2.3399130755 with <N>=2.000000 -- which is exactly the unconstrained ED global ground state (-2.3399130755495743). Digit-for-digit the pasted output.

I checked whether the hunter misread the intended design, and they did not: the guarantee is verbatim in promote_to_dense's own docstring (manybodychain.py ~570) AND repeated in docs/user_guide.md:514-519 -- "the ground-state energy and wavefunction are kept -- a bare gs_energy() afterwards therefore still returns the sector's energy rather than re-solving; call restart() if an unconstrained re-solve is what you want" -- and user_guide.md:540 states promote_to_dense is available on itensor_version=3 AND "python" alike, so python is squarely inside the promise. Mechanism verified in source: promote_to_dense sets self._session_ham_cache=None, the next gs_energy re-sends the Hamiltonian, pyitensor/chain.py:475-483 set_hamiltonian clears self._wf0_energy, and chain.py:485-489 gates the skip on `self._wf0_energy is not None`, so skip_dmrg falls through to a full dmrg() from the promoted state -- which, python's sector being a charge PENALTY on the variational solve rather than a structural confinement (pyitensor/sector.py), is now free to leave the sector. v3 is rescued only because its DMRG cannot leave a sector H commutes with. Nothing in ROADMAP.md, the audit, or docs/known_issue_*.md covers it. The hunter's supporting timing claim (neither backend actually skips the re-solve; 0.063s/0.103s vs 0.0000s) I did not re-time -- reported, not independently verified, and the verdict does not rest on it.


---


### 11. submode="CVM_explicit" returns exactly 2x the spectral function it is documented to share with submode="CVM", on every backend, plus an np.abs() that destroys the sign of negative-weight correlators

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `python-backend-parity`

**Where**: `src/dmrgpy/nonhermitian/dynamics.py:7-30 (dynamical_correlator_cvm_explicit, the `return es, np.abs(1j*outz/np.pi)` at line 28), re-exported at src/dmrgpy/cvm.py:315 and dispatched at src/dmrgpy/dynamics.py:60`

docs/user_guide.md:1425 documents submode="CVM_explicit" as "Computes the same G_AB(omega) as submode='CVM'" and recommends it precisely as the cross-check for a suspicious CVM curve. It returns exactly twice that, at every frequency, on itensor_version=3 and on "python" alike (measured ratio 2.00011 / 2.00044 / 2.00047 / 2.00214 against the exact ED/EX Lehmann reference; CVM, ROOTN and ED's own CVM all sit at 1.0001).

The arithmetic is visible in one line. The routine forms outz = f(w,+delta) - f(w,-delta) = G(w+i*delta) - G(w-i*delta) = 2i*Im G(w), then returns 1j*outz/pi = -2*Im G(w)/pi = 2*S_AB(w). The factor 1/2 that turns the two-sided difference into Im G is missing.

A second defect in the same return: np.abs() is applied unconditionally. That is invisible for a diagonal correlator (S_AB >= 0), but it silently flips the sign of any off-diagonal pair whose spectral weight is negative -- and S_AB for A != B routinely is. The same function is also what a Hermitian chain reaches (it prints "Non Hermitian mode in dynamical correlator" even for a Hermitian H, because dynamics.py routes CVM_explicit into nonhermitian/dynamics.py unconditionally), so both defects are on the default Hermitian path, not only the non-Hermitian one.

This is backend-independent shared Python -- I found it via the parity sweep rather than because it is a python-vs-v3 divergence. "factor of 2" in a CVM context appears nowhere in docs/audit_2026_08_hole_hunt.md (its finding #10's factor-2 discussion is about TDZ), and ROADMAP.md:95 lists CVM_explicit as implemented on both C++/python.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/python-backend-parity && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 repro_cvmexp.py

--- repro_cvmexp.py (n=6 S=1/2 Heisenberg + staggered 0.3*Sz, maxm=30, nsweeps=20, es=linspace(0.2,4,8), delta=0.3) ---
for v,sub in [("ED","EX"),("ED","CVM"),(3,"CVM"),(3,"ROOTN"),(3,"CVM_explicit"),
              ("python","CVM"),("python","CVM_explicit")]:
    x,y = sc.get_dynamical_correlator(submode=sub, name=(sc.Sz[0],sc.Sz[0]), es=es, delta=0.3)
```

Observed:

```
ED       EX             [0.1084689 0.0459946 0.086541  0.0525315]   ratio to ED/EX = [1. 1. 1. 1.]
ED       CVM            [0.1084749 0.0460046 0.086561  0.0525879]   ratio to ED/EX = [1.00006 1.00022 1.00023 1.00107]
3        CVM            [0.1084749 0.0460046 0.086561  0.0525879]   ratio to ED/EX = [1.00006 1.00022 1.00023 1.00107]
3        ROOTN          [0.1084749 0.0460046 0.086561  0.0525879]   ratio to ED/EX = [1.00006 1.00022 1.00023 1.00107]
Non Hermitian mode in dynamical correlator
3        CVM_explicit   [0.2169497 0.0920093 0.1731231 0.1051757]   ratio to ED/EX = [2.00011 2.00044 2.00047 2.00214]
python   CVM            [0.1084749 0.0460046 0.086561  0.0525879]   ratio to ED/EX = [1.00006 1.00022 1.00023 1.00107]
Non Hermitian mode in dynamical correlator
python   CVM_explicit   [0.2169497 0.0920097 0.1731231 0.105176 ]   ratio to ED/EX = [2.00011 2.00045 2.00047 2.00215]
```

**Expected**: CVM_explicit must return the same numbers as CVM/ROOTN/EX -- 0.1084749, 0.0460046, 0.086561, 0.0525879 here -- since user_guide.md §6 states it computes the same G_AB(omega) by a more literal transcription and offers it specifically as a cross-check. -Im G(w)/pi from a two-sided resolvent difference is (i/(2*pi))*(G(w+i*eta) - G(w-i*eta)), i.e. half what the code computes.

**Suggested fix**: In nonhermitian/dynamics.py::dynamical_correlator_cvm_explicit change the return to `return es, 0.5j*outz/np.pi` and drop the unconditional np.abs() (keep the complex/real convention the other submodes use, so a negative-weight off-diagonal correlator survives). Add a cross-submode assertion to tests/test_dynamical_correlator.py that CVM_explicit matches CVM to a few percent on a Hermitian chain -- the existing tests never compare their amplitudes.

**Reviewer (CONFIRMED)**: Re-ran repro_cvmexp.py: CVM_explicit / (ED EX reference) = 2.00011, 2.00044, 2.00047, 2.00215 on itensor_version=3 and 2.00011, 2.00045, 2.00047, 2.00215 on "python", while CVM, ROOTN and ED's own CVM all sit at 1.0001-1.0011. The arithmetic is unambiguous and matches: nonhermitian/dynamics.py:27-28 forms outz = f(w,+delta) - f(w,-delta) = G(w+i*eta) - G(w-i*eta) = 2i*Im G, then returns 1j*outz/pi = -2*Im G/pi = 2*S_AB. docs/user_guide.md:1425 promises verbatim that it "Computes the same G_AB(omega) as submode='CVM'" and recommends it as the cross-check for a suspicious CVM curve, so the claimed expectation is the documented one.

I checked three refutation routes and all failed. (a) Not a regression from the audit's fix: commit 1b87543 touched only src/dmrgpy/dynamics.py (dispatch), cvm.py:315 is a bare re-export of the nonhermitian function and its git history shows no removed 0.5. (b) The audit's "max|EX-CVM_explicit| = 2.79e-08" (line 119) does NOT contradict this -- reading lines 100-135, that control ran on a NON-Hermitian H (h + 0.3j*Sz[0]) in exactly the regime where every submode was literally the same function object, so the agreement is trivial and says nothing about the amplitude on a Hermitian chain. "factor of 2" in a CVM context appears nowhere in the audit, and ROADMAP.md:95 marks CVM_explicit implemented on v3/python. (c) Not convergence: the ratio is 2.000 to four digits at four separate frequencies.

The np.abs() SUB-CLAIM IS REFUTED, however, and should be dropped from the finding. The routine raises unless A^dagger == B (`is_zero_operator(A.get_dagger()-B)` at line 19); with A^dag = B the weight is -Im <psi|(w+E0+i*eta-H)^-1|psi>/pi with |psi> = B|GS>, i.e. (eta/pi)*<psi|((w+E0-H)^2+eta^2)^-1|psi> >= 0 for Hermitian H, so np.abs is a no-op and cannot flip any sign. The "off-diagonal pair whose spectral weight is negative" case the hunter describes is rejected by the guard before it can reach the abs. On a genuinely non-Hermitian H the two-sided difference is complex and abs is lossy, but that is arguably deliberate for the path that prints "Non Hermitian mode". The factor of 2 carries the finding on its own.


---


### 12. Parafermionic_Chain.get_dynamical_correlator ignores self.mode="ED" and dispatches to DMRG anyway, aborting the whole process with SIGABRT on itensor_version=3

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `python-backend-parity`

**Where**: `src/dmrgpy/parafermionchain.py:40-45 (the get_dynamical_correlator override -- `mode="DMRG"` default with no self.get_mode() consultation), vs src/dmrgpy/manybodychain.py:909 which does call get_mode()`

Every dispatching method on Many_Body_Chain resolves the solver through self.get_mode(mode=mode), so an enforced self.mode="ED" wins over the per-call default. Parafermionic_Chain overrides get_dynamical_correlator with `def get_dynamical_correlator(self, mode="DMRG", **kwargs)` and branches on the raw `mode` argument only -- self.mode is never read. So on a chain that has been put into ED mode (or that mode.py would route to ED for any of its other reasons), the call still goes to the DMRG backend.

Because the ED path never pushes the Hamiltonian into the C++ session, that lands in Chain::kpm_dynamical_correlator with have_H_ false, and ITensor's Error() calls abort(): the user's entire Python process dies with SIGABRT and a core dump, uncatchable from Python. On itensor_version="python" the same dispatch raises a catchable RuntimeError instead, which is how the divergence surfaced in the cross-backend sweep. Passing mode="ED" explicitly works fine on both, so the ED implementation exists and is correct -- only the implicit route is broken.

This is the same shape as audit finding #14 (get_distribution/get_distribution_moments never calling get_mode()) but a different method, a different class, and a worse symptom (process abort rather than AttributeError). "Parafermionic" appears zero times in docs/audit_2026_08_hole_hunt.md. I checked the other model classes: this is the only get_dynamical_correlator/get_gs/gs_energy/vev override in the tree that skips get_mode().

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/python-backend-parity && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 repro_pf.py ; echo "shell saw exit: $?"

--- repro_pf.py (Parafermionic_Chain(4, Z=3), H = sum Sig_i Sigd_{i+1} + 0.4 sum Tau_i + h.c.) ---
sc.set_hamiltonian(h); sc.maxm=64; sc.nsweeps=25
sc.mode = "ED"
print(sc.get_mode())                       # -> ED
sc.get_dynamical_correlator(mode="ED", name=(sc.Tau[1],sc.Tau[1]), es=..., delta=0.3)   # works
sc.get_dynamical_correlator(          name=(sc.Tau[1],sc.Tau[1]), es=..., delta=0.3)   # aborts
```

Observed:

```
itensor_version 3 mode ED get_mode ED
gs_energy(ED) -4.367430075242429
now calling get_dynamical_correlator with mode='ED' explicitly...
ED ok [5.58537966e-05-9.79464992e-18j 4.81727617e-04-2.39614256e-17j]
now calling get_dynamical_correlator with NO mode (self.mode='ED' set)...
From line 2558, file chain_session.h

Chain::kpm_dynamical_correlator called before set_hamiltonian

Chain::kpm_dynamical_correlator called before set_hamiltonian
timeout: the monitored command dumped core

(same chain on itensor_version="python" -- catchable, but still the wrong dispatch:)
  File ".../dmrgpy/pyitensor/chain.py", line 1413, in kpm_dynamical_correlator
    raise RuntimeError("Chain.kpm_dynamical_correlator called before set_hamiltonian")
RuntimeError: Chain.kpm_dynamical_correlator called before set_hamiltonian
```

**Expected**: With self.mode="ED" set on the chain, sc.get_dynamical_correlator(...) must run the ED implementation and return the same array the explicit mode="ED" call returns ([5.585e-05, 4.817e-04, ...]), exactly as Many_Body_Chain.get_dynamical_correlator does for every other model class. Under no circumstances should a public API call abort the user's process.

**Suggested fix**: In parafermionchain.py's override, resolve the solver the way the base class does -- `mode = self.get_mode(mode=mode)` as the first line -- so an enforced self.mode (and mode.py's own DMRG->ED fallbacks: no compiled extension, v3 with ns<3) are honoured. Add a test asserting that `Parafermionic_Chain(4).mode="ED"; get_dynamical_correlator(...)` equals the explicit-mode result. Independently worth doing: the `if have_H_` guards in chain_session.h that reach ITensor's Error() could throw std::runtime_error instead, so a mis-dispatch anywhere else surfaces as a Python exception rather than a core dump.

**Reviewer (CONFIRMED)**: Re-ran repro_pf.py thread-pinned and reproduced the abort exactly: get_mode() prints "ED", gs_energy(ED) = -4.367430075242429, the explicit mode="ED" correlator returns [5.58537966e-05, 4.81727617e-04], and the same call with no mode= (self.mode="ED" set) prints "From line 2558, file chain_session.h / Chain::kpm_dynamical_correlator called before set_hamiltonian" and dumps core. Digit-for-digit the pasted output.

The code is dispositive: parafermionchain.py:40 is `def get_dynamical_correlator(self,mode="DMRG",**kwargs)` branching on the raw argument, while manybodychain.py:911 -- the method it shadows -- opens with `mode = self.get_mode(mode=mode)`. I verified the hunter's "only one in the tree" claim with `grep -rn "def get_dynamical_correlator" src/dmrgpy/*.py`: the only model-class override is parafermionchain.py's, and it is the only one skipping get_mode(). That also means mode.py's own DMRG->ED fallbacks (extension not compiled, v3 with ns<3) are bypassed on this class. Scope: this is a different method, class and symptom from audit finding #14 (get_distribution/get_distribution_moments), "Parafermion" appears nowhere in the audit, and it is not ROADMAP-marked or in any docs/known_issue_*.md. Not attributable to vendored ITensor either -- the abort is ITensor's Error() being reached by a dmrgpy-side mis-dispatch, not an ITensor defect.


---


### 13. pyitensor's applyMPO is an unorthogonalized zip-up: truncating at a provably lossless bond dimension loses 9.4e-5, so applyoperator/vev(npow>1)/gs_energy_fluctuation are wrong at default settings on itensor_version="python"

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `wavefunction-consumers`

**Where**: `src/dmrgpy/pyitensor/mpsalgebra.py:453-472 (_apply_chain's SVD loop), :477-489 (applyMPO docstring), src/dmrgpy/pyitensor/chain.py:1567-1575 (_apply_mpo/_apply_mpo_with), consumers at chain.py:583,589,668,721,771,818-819,827,843,856-857,908-909,964-965,1001,1036,1423-1536,1670,1761-1769,1855-1883`

`applyMPO(K, x)` is implemented by `_apply_chain`, a left-to-right zip-up: both K and X are put at `center=1`, then at each bond the running product `piece = leftover*K.A(i)*X.A(i)` is SVD-truncated to `cutoff`/`maxdim`. The right environment of that SVD is *not* orthogonal — K and X being individually right-canonical does not make the product K.X right-canonical — so the truncation is not the optimal (2-norm-minimizing) truncation, and discarded singular weight does not measure the actual error. The docstring nevertheless claims the routine is "exact up to cutoff/maxdim regardless of x0's value, so ignoring it is correctness-preserving", and it silently discards the variational seed `x0` that would make it so. Measured on an 8-site chain where the exact H|psi> has Schmidt rank exactly [2,4,8,16,8,4,2] — i.e. truncating to maxdim=16 removes literally nothing — `applyMPO(maxdim=16)` still loses 9.4e-05 in 2-norm, while maxdim=32 gives 1.6e-14. Because `_apply_mpo` passes the chain's own `self.maxm`, every public consumer inherits this: `sc.applyoperator` (`||H|gs>||^2` off by 8.8e-09 at maxm=16 vs 2.5e-14 on v3), `sc.vev(MO, npow>1)`, and `sc.gs_energy_fluctuation()`, which at stock defaults (maxm=30, nsweeps=15) on a 10-site Heisenberg chain reports 5.06e-05 for a state that is an exact eigenstate — 420x the v3 value of 1.19e-07, which is just the double-precision cancellation floor of subtracting 18.13 from 18.13. The user guide advertises that number as "a measure of how sharply the DMRG/ED state is an eigenstate", so a user tuning maxm by watching it would conclude the state is converged to 1e-4 when it is converged to 1e-15. Building the operator explicitly instead (`sc.vev(h*h)`, which needs no MPO application) gives 18.130863826435 against v3's 18.130863826376, confirming the loss is in the application, not the MPO or the state. The same primitive backs the KPM Chebyshev recursion (`_apply_mpo_with` at chain.py:1423-1536, 1855-1883, hundreds of applications per correlator), the MPO-Taylor `exponential_apply`/`evolve_taylor_step`/`quench` paths, and CVM/`apply_inverse` — not separately measured, but they apply the same non-optimal truncation repeatedly, where the error compounds. Default TDVP time evolution is NOT affected (it Krylov-propagates `two_site_heff` directly).

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/wavefunction-consumers/repro_f2.py   # and repro_f2b.py for the mechanism (both standalone)
```

Observed:

```
$ python3 repro_f2.py
python  <H>=-4.258035207272  vev(h,npow=2)=18.130863823805  vev(h*h)=18.130863826435  gs_energy_fluctuation=5.0600e-05
3       <H>=-4.258035207251  vev(h,npow=2)=18.130863826190  vev(h*h)=18.130863826376  gs_energy_fluctuation=1.1921e-07

$ python3 repro_f2b.py
exact Schmidt rank of H|psi> per bond: [2, 4, 8, 16, 8, 4, 2]
applyMPO(maxdim=16  )  || H|psi> - exact ||_2 = 9.395e-05
applyMPO(maxdim=32  )  || H|psi> - exact ||_2 = 1.608e-14
applyMPO(maxdim=None)  || H|psi> - exact ||_2 = 1.608e-14
```

**Expected**: applyMPO at a bond dimension that is provably sufficient must be exact to machine precision, as ITensor's own applyMPO is (5.9e-14 at the identical maxm=16 on itensor_version=3). gs_energy_fluctuation() on a state that is an exact eigenstate to 1e-15 must return the cancellation floor (~1e-7 for <H^2>~18), not 5e-05. Known because (a) the exact Schmidt spectrum of H|psi> was computed densely and caps at the physical bound at every bond, and (b) raising maxdim above the cap, with cutoff held at 0, drops the error by 10 orders.

**Suggested fix**: In `_apply_chain`, do the zip-up contraction at full bond dimension (`maxdim=None`, `cutoff=0`) and then apply the requested `cutoff`/`maxdim` with a proper `position()` sweep, which truncates in an orthogonal gauge — or honor `x0` and run the variational fit the parameter exists for. Either way correct the `applyMPO` docstring: it is currently a zip-up, not "exact up to cutoff/maxdim". Separately, `gs_energy_fluctuation` should not go through `vev(npow=2)` at all — `vev(h*h)` builds H^2 as an MPO and needs no MPO application (2.1e-14 error at the same settings).

**Reviewer (CONFIRMED)**: I re-ran both repros thread-pinned and reproduced them to the printed digit (repro_f2: python fluct 5.0601e-05 vs v3 2.5288e-07; repro_f2b: 9.395e-05 at maxdim=16, 1.03e-14 at 32/None). I then tried three ways to explain it away and failed on all three.

(a) Not DMRG convergence error. repro_f2b holds ONE DMRG state fixed and only varies applyMPO's maxdim, so no state-quality argument applies. I strengthened the losslessness claim the hunter only asserted via a rank count: I printed the 17th singular value of the exact dense H|psi> at every bond and it is identically 0.0e+00 at all seven bonds, so maxdim=16 discards literally zero weight, yet the zip-up loses 9.4e-05.

(b) The proposed fix is real and available. I monkey-free-tested it: applyMPO(...,maxdim=None) followed by position(n) then position(1,maxdim=md) (truncation in an orthogonal gauge) gives 1.007e-14 at maxdim=16, 20 and 32 alike. So the mechanism is exactly as diagnosed: the zip-up keeps the top singular vectors of `piece`, whose right environment (a product of two individually right-canonical chains, which is not right-canonical) is non-orthogonal, so at the first bond where 2*D_prev exceeds maxdim (bond 5 here) it can discard directions the final state actually needs.

(c) Not a deliberate documented simplification. mpsalgebra.py's module docstring asserts the opposite invariant -- 'the standard, correct way to compress a tensor-train sum/product... dmrgpy only ever observes final numerical results bounded by Cutoff/MaxDim' -- and _apply_chain's own docstring claims the interleaved and every-site-first variants 'both reach the same fixed point'. That last claim is true (both truncate in a non-orthogonal gauge) but the invariant it is used to justify is not. applyMPO's docstring's 'exact up to cutoff/maxdim regardless of x0's value, so ignoring it is correctness-preserving' is the claim this refutes. Nothing about it is recorded in docs/known_issue_*.md, docs/audit_2026_08_hole_hunt.md or ROADMAP.md. Parity check: mpscpp3's apply_operator uses the same maxm_/cutoff_ as pyitensor's _apply_mpo, so the backend comparison is apples-to-apples.

Sweeping maxm (my own script, both backends, maxm in 16/30/40/64/128) does shrink the python number (4.8e-3, 5.06e-5, 3.3e-6, then a ~5.5e-6 floor), but that does not rescue it: the failure at a mathematically sufficient bond dimension is an algorithmic defect, not a tolerance.

Unverified parts of the finding, flagged as such: the compounding claim for KPM/MPO-Taylor/CVM is stated as not measured and I did not measure it either; and the reported 'applyoperator ||H|gs>||^2 off by 8.8e-09 at maxm=16' did not reproduce on my chain (I measured 2.3e-05 for python vs 1.1e-08 for v3 in <H^2> at maxm=16 on a 10-site chain), so that particular number should not be quoted. Also, my v3 gs_energy_fluctuation came out 1.0e-07 to 2.5e-07 across runs rather than the reported 1.19e-07 -- that is randomMPS-start noise on the v3 backend, not a refutation.

**Reviewer on severity**: Keep medium, but the stated 'expected' is wrong on one point and should not be quoted. The finding says gs_energy_fluctuation 'must return the cancellation floor (~1e-7)' on python. It cannot, for a reason independent of applyMPO: I densified the pyitensor DMRG state and computed its exact fluctuation ||H|psi> - E|psi>|| with no MPO application at all -- 8.364e-06 at maxm=30 and 5.918e-06 at maxm=128, i.e. the ~1e-6 Lanczos eigenvector cap CLAUDE.md already documents. So the inflation attributable to applyMPO at stock defaults is ~6x (reported 5.06e-05 vs the state's true 8.4e-06), not the 420x the finding gets by comparing against v3 -- that comparison conflates a better-converged v3 state with a better MPO application. Medium still stands on the algorithmic defect plus the two docstrings that assert the violated invariant.


---


### 14. Many_Body_Chain.exponential() never reaches its DMRG exponential path for any multi-site Hamiltonian, because MultiOperator.is_hermitian() false-rejects every two-site term; it silently returns an unnormalized 2-term Taylor series instead, 4% wrong at z=1 and unbounded

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `wavefunction-consumers`

**Where**: `src/dmrgpy/mpsalgebra.py:9-17 (the dispatch and its fallback), src/dmrgpy/multioperator.py:65-72 (is_hermitian/is_antihermitian), src/dmrgpy/manybodychain.py:837-839 (the public method)`

`mpsalgebra.exponential(self,h,wf)` dispatches on the *symbolic* `MultiOperator.is_hermitian()` / `is_antihermitian()`, which test `(self -/+ self.get_dagger()).simplify() == 0`. `simplify()` does not recognize that `get_dagger()`'s operator-order reversal is a no-op for factors on different sites, so `Sx[i]*Sx[j]+Sy[i]*Sy[j]+Sz[i]*Sz[j]` — the single most common Hamiltonian shape in this library — reports `is_hermitian() == False` while the chain's own numerical probe `sc.is_hermitian(h)` correctly reports True. Both branches therefore fail and control falls into the `else`, which prints "Warning, using 3rd order taylor expansion mode" and returns `wf + h*wf + h*h*wf/2` — a 2-term (mislabelled 3rd-order) Taylor truncation with no step subdivision, no normalization, and no dependence on ||h||, instead of `exponential_dmrg`'s converged expansion with its `nt0 = bandwidth*nt` sub-steps. `sc.exponential(h, wf)` is a documented public method (`docs/user_guide.md:210`, "e^{h}|psi>") and `examples/time_evolution/exponential_EV/main.py` is built entirely around comparing it between DMRG and ED — that example happens to escape the bug only because its operator `sum(sc.Sx)` is a sum of *single-site* terms, for which `is_hermitian()` does work. Replace it with any two-site H and the DMRG/ED comparison the example makes diverges: 9.6e-04 at z=0.25, 6.7e-03 at z=0.5, 4.1e-02 at z=1.0, growing without bound in z (and, used in a time-evolution loop, compounding per step). ED mode is unaffected (it never consults the symbolic test), as is `timeevolution.evolve_WF`'s DMRG branch, which calls `exponential_dmrg` directly. NOTE on prior art: `infinitechain.py:540-555` documents the `is_hermitian()` false-negative itself as a known, deliberately-not-fixed limitation *of that module's own Hamiltonian check*; the consequence in `mpsalgebra.exponential` — a public method that silently returns a badly wrong number on all three backends — is not recorded anywhere (not in docs/audit_2026_08_hole_hunt.md, ROADMAP.md, or any docs/known_issue_*.md).

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/wavefunction-consumers/repro_f1.py   # standalone
```

Observed:

```
heisenberg.is_hermitian() = False    sc.is_hermitian(h) = True
Warning, using 3rd order taylor expansion mode
z=0.25  <gs|e^{zH}|gs>  ED=1.20623025  DMRG=1.20507812  relerr=9.551e-04
Warning, using 3rd order taylor expansion mode
z=0.50  <gs|e^{zH}|gs>  ED=1.45499141  DMRG=1.44531250  relerr=6.652e-03
Warning, using 3rd order taylor expansion mode
z=1.00  <gs|e^{zH}|gs>  ED=2.11700002  DMRG=2.03125000  relerr=4.051e-02
```

**Expected**: `sc.exponential(z*h, wf)` on a Hermitian Hamiltonian should take the `exponential_dmrg` path and agree with `mode="ED"` to DMRG's own accuracy (the sibling single-site case in examples/time_evolution/exponential_EV does exactly that). It should never silently substitute a 2-term Taylor series whose error is O((z||h||)^3) with no convergence control; if a path genuinely is unavailable it should raise, not print a warning and return a number.

**Suggested fix**: Gate `mpsalgebra.exponential` on the chain's numerical `self.is_hermitian(h)` (already used everywhere else in the codebase, e.g. excited.py:100, and correct here) rather than the symbolic `h.is_hermitian()`, with `self.is_hermitian(1j*h)` for the anti-Hermitian case; and make the remaining `else` branch `raise NotImplementedError` instead of printing a warning and returning the unconverged Taylor truncation. Fixing `MultiOperator.is_hermitian()`'s `simplify()` to canonically order factors on distinct sites would fix this and `infinitechain.py`'s documented gap at the same time.

**Reviewer (CONFIRMED)**: Reproduced verbatim, including the 'Warning, using 3rd order taylor expansion mode' print and relerr 9.551e-04 / 6.652e-03 / 4.051e-02. I replaced the hunter's ED reference with an independent scipy.linalg.expm on the dense 16x16 Hamiltonian: it agrees with mode='ED' to 8 digits, so the reference is sound and the DMRG dispatch really is 4% off at z=1.

Scope: the ROOT CAUSE is documented in-source twice -- infinitechain.py:540-555 ('a pre-existing, general limitation, not specific to this module... Left as a known, documented gap') and docs/documentation.md:1743-1749, which routes around the same false-negative for kpm_finite. Neither is one of the three out-of-scope markers in my brief (no known_issue_*.md entry, nothing in docs/audit_2026_08_hole_hunt.md -- I grepped 'exponential' there and got zero hits -- and nothing in ROADMAP.md), and neither records this consequence. Also confirmed the example that exercises this method, examples/time_evolution/exponential_EV/main.py, escapes only because sum(sc.Sx) is single-site: I checked Mx.is_hermitian() symbolically and it is True.

One thing the hunter got wrong, and it matters for anyone acting on this. The finding's 'expected' is that routing to exponential_dmrg would agree with ED. It would not. Calling mpsalgebra.exponential_dmrg(sc, z*h, wf) directly (bypassing the symbolic gate) on the default v3 backend gives 0.829 / 0.687 / 0.472 against the dense reference 1.206 / 1.455 / 2.117 -- relative errors 3.1e-01, 5.3e-01, 7.8e-01. The cause is a sign convention, not truncation: exponential_dmrg sets tau = complex(-dt.real, dt.imag) = -1 for its dt=1.0 Hermitian call, chain_session.h's custom_exp(H,z) builds I + zH + z^2 H^2/2 = e^{zH}, so the DMRG path computes e^{-h} while edchain.exponential does algebra.expm(h) = e^{+h} and user_guide.md:210 documents e^{h}. I confirmed the sign independently on an Mz eigenstate (<Mz>=-2): at z=0.5 DMRG returns 2.71828155 = e^{+1} while ED returns 0.36787944 = e^{-1}. The exponential_EV example cannot see this because e^{+zMx} and e^{-zMx} have identical expectation in a Z-polarized state (I checked: both 1.13169806 at z=0.5).

So the finding is real -- a public, documented method silently returns an uncontrolled 2-term Taylor truncation for the most common Hamiltonian shape in the library -- but its suggested fix (gate on the numerical self.is_hermitian(h)) would, applied as written today, replace a 4% error with a 31-78% sign-flipped one. Fix the tau sign first, or the 'fix' is a regression.


---


### 15. mode.py's automatic v3 ns<3 -> ED fallback yields a chain whose get_gs() returns an ED State that overlap/aMb/get_rdm/get_distribution cannot consume, because those modules branch on self.mode instead of get_mode()

`hole` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `dispatch-matrix` &middot; independently found by `cpp-v3-completeness`

**Where**: `src/dmrgpy/mpsalgebra.py:52-63 (overlap/overlap_aMb read self.mode); src/dmrgpy/densitymatrix.py:4-13 (reduced_dm reads neither); src/dmrgpy/manybodychain.py:1128-1141 (random_state reads self.mode); src/dmrgpy/mode.py:88-92 (the fallback)`

`mode.py::resolve_mode` falls back to ED for `itensor_version==3 and self.ns<3` without touching `self.mode`, which stays `None`. Several modules never call `get_mode()` and instead test `self.mode` directly (`if self.mode is not None: mode = self.mode`) or not at all, so they take the DMRG branch on a chain whose ground state is now an ED `State`. The result is not a clear refusal but an `AttributeError: 'State' object has no attribute 'cpp_handle'` several frames deep — and, worse, two public methods on the SAME chain now disagree about what kind of object it holds: `get_gs()` returns a `State` while `random_state()` returns an `MPS`, so `sc.overlap(sc.get_gs(), sc.random_state())` cannot work under any branch. This is audit finding #14's class ("the only public entry points that never call get_mode()"), fixed there for `get_distribution*` and still open for the whole `mpsalgebra`/`densitymatrix` family. `itensor_version=2` answers all of these correctly at the same size, so the divergence is purely the fallback's. `get_distribution(n=...)` shows the second half of the problem: the ED route audit #14 installed leads to `distribution_kpm()`, whose signature does not accept the same keywords the DMRG route does.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 3 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/dispatch-matrix/matrix.py 2
```

Observed:

```
=== N=2 backend=v3 ===          # user set only itensor_version=3; mode is None
  get_gs                             obj State
  get_rdm(i=0)                       EXC AttributeError: 'State' object has no attribute 'cpp_handle' | depth=4 at densitymatrix.py:12
  get_distribution                   EXC TypeError: distribution_kpm() got an unexpected keyword argument 'n' | depth=5 at distribution.py:14
  overlap(gs,gs)                     EXC AttributeError: 'State' object has no attribute 'cpp_handle' | depth=5 at mpsalgebra.py:68
  aMb                                EXC AttributeError: 'State' object has no attribute 'cpp_handle' | depth=6 at mpsalgebra.py:87
  random_state                       obj MPS

=== N=2 backend=v2 ===          # identical chain, identical calls
  get_gs                             obj MPS
  get_rdm(i=0)                       arr(2, 2) [ 0.5+0.j -0. +0.j -0. +0.j]
  get_distribution                   seq[2] [-0.61875 -0.61256 -0.60638]
  overlap(gs,gs)                     val (1+0j)
  aMb                                val (-0.75+0j)
  random_state                       obj MPS
```

**Expected**: Once mode.py has decided this chain is answered by ED, every method must follow that decision: `overlap`, `aMb`, `get_rdm` and `random_state` should all take their ED branch (each of which exists and works — the explicit `sc.mode="ED"` column answers overlap, aMb, get_site_entropy and get_pair_entropy correctly at n=4), so the results should match the v2 column. At minimum the failure should be a named refusal naming the fallback, not an `AttributeError` on an attribute of an object the user never asked for.

**Suggested fix**: Replace the `if self.mode is not None: mode = self.mode` idiom in `mpsalgebra.py` (`overlap`, `overlap_aMb`, `exponential`) and the unguarded `self._session` dereference in `densitymatrix.reduced_dm` and `manybodychain.random_state` with `mode = self.get_mode(mode=mode)`, the way `gs_energy`/`get_gs`/`vev`/`get_dynamical_correlator` already do — `get_mode()` is the only thing that sees the ns<3 and extension-not-compiled fallbacks. Separately, align `distribution_kpm`'s signature with the DMRG entry point so audit #14's ED route accepts the documented `n=`/`X=`.

**Reviewer (CONFIRMED)**: Core claim reproduced on my own script (not the hunter's matrix.py). A 2-site Heisenberg chain with only itensor_version=3 set and self.mode left None: get_mode() returns ED, gs_energy() correctly returns -0.75, get_gs() returns an edtk State — and then overlap(gs,gs), aMb(gs,h,gs) and get_rdm(i=0) all die with `AttributeError: 'State' object has no attribute 'cpp_handle'`, while random_state() on the SAME chain returns an MPS. v2 at the identical size answers all of them (overlap (1+0j), aMb (-0.75+0j), rdm [0.5,0]), so the divergence is the fallback's alone. Source confirms the idiom: mpsalgebra.py:52-63 uses `if self.mode is not None: mode = self.mode`, which honours an explicit sc.mode="ED" but is blind to mode.py's automatic fallbacks — get_mode()/resolve_mode() (mode.py:88-92) returns "ED" without ever writing self.mode. I checked that the ED branches these calls should reach genuinely work: on an n=4 chain with explicit sc.mode="ED", overlap -> (1+0j), aMb -> (-1.616025+0j), random_state -> State. Not already recorded: audit #14 covers only get_distribution/get_distribution_moments (both now route through get_mode) and its reviewer explicitly carved the entropy leg out as a different root cause; overlap/aMb/random_state appear nowhere in that document (grepped).

**Reviewer on severity**: Agree with medium for the core, but two sub-claims in the write-up should be struck and do not survive on their own. (1) get_rdm: my run shows densitymatrix.reduced_dm fails with the same AttributeError even under an EXPLICIT sc.mode="ED" on an n=4 chain — it reads neither self.mode nor get_mode() and there is no ED branch to route to. That puts it in the same already-carved-out class as entropy.py ("_session-only, no ED path of its own", CLAUDE.md), not in the self.mode-idiom class, so the finding's "each of which exists and works" is false for this one leg and adding get_mode() there would fix nothing. (2) get_distribution(n=...): I read the signatures — `n=` is not a parameter of the DMRG route either; kpmdmrg.general_kpm_moments merely swallows it in **kwargs and ignores it, while edtk/distribution.py's distribution_kpm(wf0,X,scale,delta,xs) rejects it loudly. So "the ED route does not accept the keywords the DMRG route does" is wrong as stated — DMRG silently ignores an undocumented kwarg. Severity of what remains (overlap/aMb/random_state, plus the get_gs->State / random_state->MPS split on one chain) is medium as reported.

<details><summary>Corroborating report from lens <code>cpp-v3-completeness</code>: get_rdm() has no ED path and no get_mode() guard: crashes with the exact AttributeError the 2026-08 audit fixed for its neighbours, and rejects mode= entirely</summary>

**Where**: `src/dmrgpy/densitymatrix.py:12 (reduced_dm), src/dmrgpy/manybodychain.py:1118-1122 (get_rdm)`

`densitymatrix.reduced_dm` branches on `itensor_version=="julia_live"` and otherwise goes straight to `self._session.reduced_dm(wf.cpp_handle, i+1)`. There is no `get_mode()` call and no ED branch, so every path on which mode.py routes the chain to ED hands an ED `State` to session-only code. The v3-specific automatic fallback is the live one: CLAUDE.md documents that `mode.py::get_mode()` silently routes `itensor_version==3` with `self.ns<3` to ED (ITensor v3's two-site dmrg aborts below 3 sites), and `get_rdm` is the entry point that sweep missed. The failure string is verbatim the one `get_distribution`'s own fix comment quotes as the symptom it exists to prevent (manybodychain.py:957-965: "get_mode(), like every neighbouring entry point: without it the DMRG branch ran even when mode.py had routed everything else to ED ... failing with the opaque 'State' object has no attribute 'cpp_handle'"). Separately, `get_rdm` does not accept `mode=` at all, so the documented cross-validation idiom `get_rdm(mode="ED")` is a TypeError, unlike essentially every neighbouring method. Note this is NOT ROADMAP-marked absent: ROADMAP.md:74 lists "Bond entanglement entropy / reduced density matrix" as available on all three DMRG backends with "backend-agnostic dispatch", and its columns are v3/pyitensor/Julia — ED is not a column, and the ED `EDchain`/`MBFermion` objects have no `_session` at all. The same file's `reduced_dm_projective` DOES work under ED routing (it is what `get_site_entropy`/`get_mutual_information` use, and those pass), but it is not a drop-in: its basis is the projector list `[N, Cdag]`/`[Sz+1/2, S+]`, so its diagonal comes out in the opposite order from ITensor's (measured on a 6-site fermion chain in sector mode: v3 `get_rdm(i=1)` diagonal `[0.425147, 0.574853]` = (empty, occupied), while `<N_1>` = 0.574853 sits first in the projective basis). Eigenvalues agree — site entropy matched ED to 2.7e-13 — but the matrix does not.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/cpp-v3-completeness && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 -u final_repro.py    # section A; the script builds spinchain.Spin_Chain(["S=1/2"]*2, itensor_version=3) with a Heisenberg bond and calls gs_energy/vev/get_site_entropy/get_rdm/get_bond_entropy
```

Observed:

```
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
A v3 ns=2 (auto ED)  gs_energy         OK [-0.75+0.j]
A v3 ns=2 (auto ED)  vev               OK [0.+0.j]
A v3 ns=2 (auto ED)  get_site_entropy  OK [0.693147+0.j]
A v3 ns=2 (auto ED)  get_rdm           AttributeError: 'State' object has no attribute 'cpp_handle'
A v3 ns=2            get_rdm(mode=ED)  TypeError: reduced_dm() got an unexpected keyword argument 'mode'
```

**Expected**: Either the one-site reduced density matrix (the ED object can produce it — `reduced_dm_projective` on the same chain does, and the entropies derived from it match ED to 2.7e-13), or a `NotImplementedError` naming the reason, which is the established pattern its two neighbours already follow: `get_distribution_moments(mode='ED')` raises "get_distribution_moments has no ED implementation ... Note mode.py routes to ED on its own when the requested C++ extension is unavailable, or for itensor_version=3 on a chain with fewer than 3 sites", and `metts_vev` on the same 2-site v3 chain raises "metts_vev with itensor_version=3 can't handle a chain this short (n=2 < 3 sites)". Both of those were produced by the same run, so the pattern is live in this tree; `get_rdm` simply was not converted.

**Suggested fix**: Route `get_rdm` through `get_mode()` like its neighbours: `mode = self.get_mode(mode=mode)`; on `"DMRG"` keep the current `self._session.reduced_dm(...)`, on `"ED"` raise `NotImplementedError` naming the reason (the v3 ns<3 / uncompiled-extension routing) in the wording `get_distribution_moments` already uses. Also give `get_rdm`/`reduced_dm` a `mode=` kwarg so `get_rdm(mode='ED')` is a routed request rather than a TypeError. Do NOT silently substitute `reduced_dm_projective` — its basis ordering differs (see description); if it is used, transpose/reorder to ITensor's (empty, occupied) convention and say so in the docstring.

**Reviewer (CONFIRMED)**: Re-ran the hunter's final_repro.py section A thread-pinned (MKL/OMP/OPENBLAS=1, taskset -c 3, PYTHONPATH=.../src). Every line of the pasted output matched verbatim: gs_energy/vev/get_site_entropy OK on the auto-ED-routed 2-site v3 chain, get_rdm -> AttributeError: 'State' object has no attribute 'cpp_handle', get_rdm(mode='ED') -> TypeError: reduced_dm() got an unexpected keyword argument 'mode'.

I tried three refutation angles and none held.
(a) Out of scope? Not vendored ITensor, not a deliberately-reproduced legacy bug, not in a docs/known_issue_*.md, and not recorded as an audit item -- the only mention in docs/audit_2026_08_hole_hunt.md is line 2182's passing 'Unaffected' remark inside finding #16 ('mode="ED" chains have no get_rdm path at all'), which records the state without treating it as a finding. ROADMAP.md's columns really are v3/pyitensor/Julia (checked lines 68-76), so ED is not a marked gap either way.
(b) Legitimate absence with a legitimate error? No: grep over src/dmrgpy/edtk/ finds no reduced_dm at all, so there is no ED implementation, and the established in-tree pattern for exactly that situation is a NotImplementedError naming the routing -- manybodychain.py:975-985's get_distribution_moments does precisely this, and get_distribution above it carries the comment quoting this very AttributeError as the symptom get_mode() exists to prevent. So the neighbour-precedent claim checks out in the current tree.
(c) Corner case only? The opposite -- I found a more mainstream trigger than the hunter reported. On an ordinary 4-site Heisenberg Spin_Chain(itensor_version=3) with the compiled extension present, setting the documented sc.mode="ED" and calling get_rdm(i=0) gives the same AttributeError, while the same chain without sc.mode returns [[0.5,0],[0,0.5]] correctly. densitymatrix.reduced_dm never consults self.mode, so the plain documented way to force ED is enough; the v3 ns<3 route is just one more door.

One discrepancy, in the hunter's disfavour as a reporter but not as a finding: their observed_output omits the get_bond_entropy line. In my run 'A v3 ns=2 (auto ED) get_bond_entropy' fails with the identical AttributeError, so the get_mode() sweep missed at least one entry point beyond the one the finding names. Incomplete paste, not a mismatch.

</details>


---


### 16. mpsalgebra.applyoperator/summps test type(wf)==np.ndarray for their ED branch, but ED wavefunctions are edtk.edchain.State, so both are dead on the ED backend and end in a bare raise

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `dispatch-matrix`

**Where**: `src/dmrgpy/mpsalgebra.py:90-95 (applyoperator), 111-116 (summps); contrast :99-107 (applyinverse, which checks State and works)`

`applyoperator(self,A,wf)` and `summps(self,wf1,wf2)` dispatch on the runtime type of the wavefunction: `if type(wf)==mps.MPS: mode="DMRG" ; elif type(wf)==np.ndarray: mode="ED" ; else: raise`. The ED backend has not handed out bare `np.ndarray` wavefunctions — `EDchain.get_gs()` returns an `edtk.edchain.State`. `applyinverse`, ten lines below, imports `State` and tests for it correctly, which is the same-file proof the other two are stale. So both ED branches are unreachable and every call on an ED chain falls into the bare `else: raise`, producing `RuntimeError: No active exception to reraise` with no indication of what went wrong. This fires on an explicit `sc.mode="ED"` (so it is not merely a consequence of finding C) and again on every automatic ED fallback. It is the same `else`-doing-two-jobs / stale-precondition shape §4.10 catalogues.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 3 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/dispatch-matrix/final_D.py
```

Observed:

```
type(get_gs()) = State  is ndarray? False
  applyoperator  -> RuntimeError: No active exception to reraise
  summps         -> RuntimeError: No active exception to reraise
  applyinverse   -> State
```

**Expected**: `applyoperator` should return the ED state `A|wf>` via `self.get_ED_obj().applyoperator(A,wf)` (the branch is written and merely unreachable), and `summps` should return `wf1 + wf2` — `applyinverse`, which makes the identical dispatch with the correct type test, returns a `State` on the same chain and the same call site. A backend-unsupported case should in any event raise a named `NotImplementedError`, never a bare `raise` with no active exception.

**Suggested fix**: In `mpsalgebra.applyoperator` and `mpsalgebra.summps`, import `State` from `.edtk.edchain` and test `type(wf)==State` (or better `isinstance`), exactly as `applyinverse` at line 99 already does; and replace the trailing bare `raise` in all three with an explicit `TypeError` naming the wavefunction type received and the two it accepts.

**Reviewer (CONFIRMED)**: Unarguable in source and reproduced. mpsalgebra.py:90-95 and :111-116 both dispatch `elif type(wf)==np.ndarray: mode="ED"`, while applyinverse ten lines below (:99-107) imports State from .edtk.edchain and tests `type(wf)==State` — same file, same dispatch shape, one of them updated and two not. EDchain.get_gs (edchain.py:206-212) returns `State(self.get_gs_array(),self)`, so the ndarray branch can never be taken by anything the ED backend hands out. I re-ran final_D.py myself on an n=4 chain with explicit sc.mode="ED": type(get_gs()) = State, isinstance ndarray False; applyoperator -> RuntimeError: No active exception to reraise; summps -> same; applyinverse -> State. No numerics involved, so convergence/tolerance objections do not apply. Not vendored, not a CLAUDE.md legacy-reproduced bug, not in docs/known_issue_*.md, and absent from docs/audit_2026_08_hole_hunt.md (grepped for applyoperator/summps).

**Reviewer on severity**: Suggest LOW rather than medium. I grepped every internal caller (`grep -rn "applyoperator|summps" src/dmrgpy --include=*.py`): the only ones are kpmdmrg.py:208/210 and mps.py:38/169, all of which are on the DMRG path and pass an MPS; mpsjulialive has its own separate implementations. So no internal code path can reach the dead branch, and the blast radius is a user calling sc.applyoperator/sc.summps directly on an ED chain (explicit mode="ED", or an automatic fallback) — a loud crash, not a wrong number, with the working ED implementation sitting one line away. The defect itself is certain; only its reach is smaller than 'medium' implies.


---


### 17. Conserved sector + mode="ED": get_correlation_matrix/_entropy/_eigenvalues crash — the 2026-08 audit's own fix for finding #11 picked a session-only dmmode under a premise ("ED is unreachable") that has since become false

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `cpp-v3-completeness`

**Where**: `src/dmrgpy/entropytk/correlationentropy.py:60-68 (dmmode default) and :93/:235 (cpp_correlation_matrix -> self._session)`

`get_correlation_matrix_zeroT` resolves `dmmode=None` as: sector set -> `"full"`, otherwise -> `"fast"`. `"full"` is `cpp_correlation_matrix(wf)`, i.e. `wf.MBO._session.correlation_matrix(wf.cpp_handle)` — session-only. So on a sector-mode chain answered by ED the call dies with `'MBFermion' object has no attribute '_session'`, and everything layered on it in entanglement.py goes with it (`get_correlation_eigenvalues`, `get_correlation_entropy`, `get_correlation_entropy_density`, `get_correlated_orbitals`, `get_correlated_density`, `MPS.get_correlation_entropy_from_wf`). The sharp part is that enabling a sector BREAKS a call that works without one: the same chain, same `mode='ED'`, no sector, returns the matrix fine via the `"fast"` default. This is a regression from the audit's own remedy for finding #11 (docs/audit_2026_08_hole_hunt.md:1440-1519, "Status: FIXED -- dmmode now defaults against the chain's state instead of being hardcoded to 'fast'"). That finding's Affects line states the premise explicitly: "ED is unreachable because a sector-mode chain deliberately refuses to fall back." That premise no longer holds — CLAUDE.md now records that `mode="ED"` implements the whole sector API by a third mechanism (edtk/edchain.py restricting every operator to the sector submatrix, `tests/test_sector_conservation_ed.py`) and that "falling back to ED ... is now correct rather than forbidden". The fix's own suggested-remedy text offered `"full" (or "explicit")` as interchangeable; `"explicit"` is backend-agnostic and `"full"` is not, and the one that was chosen is the one that cannot answer on ED. docs/audit_2026_08_hole_hunt.md already carries a precedent section for exactly this shape ("A regression this audit's own fix introduced"), and the regression test pinned for #11 evidently exercises the DMRG backends only.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/cpp-v3-completeness && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 -u final_repro.py    # section B; 6-site Fermionic_Chain(itensor_version=3), t + 0.8*N N, f.mode='ED', with and without f.set_conserved_sector(Nf=3)
```

Observed:

```
B ED no sector       get_correlation_matrix  OK [0.58884 -0.j 0.426785-0.j 0.484375-0.j 0.484375-0.j]
B ED sector Nf=3     get_correlation_matrix  AttributeError: 'MBFermion' object has no attribute '_session'
B ED sector dmmode=explicit                  OK [0.58884 -0.j 0.426785-0.j 0.484375-0.j 0.484375-0.j]
B ED sector Nf=3     get_correlation_entropy AttributeError: 'MBFermion' object has no attribute '_session'
B ED sector Nf=3     get_correlation_eigenvalues AttributeError: 'MBFermion' object has no attribute '_session'
B DMRG(v3) sector    get_correlation_matrix  OK [0.58884 -0.j 0.426785-0.j 0.484375-0.j 0.484375-0.j]
```

**Expected**: The sector-restricted correlation matrix — the same numbers the v3 DMRG backend and the explicit dmmode both return on the identical chain, `[0.58884, 0.426785, 0.484375, 0.484375, ...]` on the diagonal (this ground state lies in Nf=3, so the sector and non-sector answers coincide, which is what makes the three OK lines above a valid oracle). The audit itself established this is computable in the sector ("dmmode='explicit' and dmmode='full' reproduce the ED reference to 1e-12 (v3) / 2e-9 (pyitensor)"), and CLAUDE.md documents `mode="ED"` as a first-class sector backend, so nothing here is a legitimate absence.

**Suggested fix**: Resolve the dmmode default against the resolved MODE as well as the sector, at correlationentropy.py:60-68: keep `"full"` only when `self.get_mode(**kwargs) == "DMRG"`, and use `"explicit"` when the chain will be answered by ED (confirmed above to return the right matrix on exactly this chain). A one-line alternative that is also correct is to make the sector default `"explicit"` unconditionally — it is backend-agnostic, it is what mpsjulialive/mps.py already hardcodes for the same reason, and the audit's own fix note listed it as an accepted option. Extend the #11 regression test to the `mode="ED"` sector case, which is what would have caught this.

**Reviewer (CONFIRMED)**: Re-ran section B thread-pinned; all six lines matched the pasted output exactly, including the numbers. The oracle structure is what makes this decisive rather than a convergence artefact: the same 6-site chain returns the identical diagonal [0.58884, 0.426785, 0.484375, 0.484375] under (i) mode='ED' with no sector, (ii) mode='ED' + sector + dmmode='explicit', and (iii) v3 DMRG + sector -- and only mode='ED' + sector with the default dmmode raises AttributeError: 'MBFermion' object has no attribute '_session'. Three independent routes agreeing to the printed digits rules out 'it is merely DMRG convergence error', and it rules out 'the quantity is not computable in this sector'.

I verified the regression narrative rather than taking it on trust. correlationentropy.py:60-68 does resolve dmmode=None to 'full' whenever conserved_sector is truthy, and 'full' is cpp_correlation_matrix at :232-235, i.e. wf.MBO._session.correlation_matrix(...) -- session-only, unconditionally. The audit fix that introduced that default is 1b87543 (2026-08-30, 'Fix the 21 holes a five-lens cross-backend audit found'), and mode="ED" sector support arrived afterwards in d62a306 (2026-09-02, 'Let mode="ED" target a conserved sector too'); git merge-base --is-ancestor confirms the ordering. So finding #11's stated premise ('ED is unreachable because a sector-mode chain deliberately refuses to fall back') was true when written and has since lapsed, exactly as claimed. The audit's own remedy text at :1519 does offer 'full (or "explicit")' as interchangeable, and the backend-agnostic one is the one that was not chosen.

Also verified the hunter's assertion that the pinned regression test covers DMRG only: tests/test_audit_2026_08_regressions.py::test_sector_correlation_matrix_default_mode_works builds hopping_chain(n=4) and calls get_correlation_matrix() with no mode= anywhere in the test, so nothing in it would have caught the ED route. Not vendored, not legacy-reproduced, not in a known_issue doc, and not an open audit item -- #11 is marked FIXED, and this is a distinct failure its fix created.


---


### 18. Bosonic_Chain with any maxnb != 4 on itensor_version=2 aborts the whole process (SIGABRT, exit 134) instead of the documented "silently uses the fixed 4-level site"

`hole` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `ed-and-operators` &middot; independently found by `recent-commits`

**Where**: `src/dmrgpy/mpscpp2/get_sites.h:110-112 (dmrgpy's own header, not vendored ITensor); src/dmrgpy/bosonchain.py:15-20; docs/user_guide.md:117-129`

`bosonchain.Bosonic_Chain` encodes the per-site dimension as the site type code 100+maxnb[i]. `mpscpp2/get_sites.h`'s `SpinX` constructor accepts the single literal code 104 (`else if (nm==104) sites.set(i,BosonFourSite(i))`) and falls through to `Error(format("SpinX cannot read index of size "))` for anything else. ITensor's `Error()` calls `abort()`, so this is an uncatchable SIGABRT that takes the interpreter down -- the same shape as audit finding #21 (maxm<=0 on v2/v3).

docs/user_guide.md:124-127 describes the v2 behaviour as "`itensor_version=2` and the Julia backend still only understand the single fixed 4-level boson site regardless of what `maxnb` requests, so a non-default `maxnb` should be run under `itensor_version=3` ... for DMRG/ED results to actually agree", and bosonchain.py:18-19 repeats it. That wording describes a silently-wrong-but-finite answer; the reality is a core dump. Nothing in mode.py or Bosonic_Chain.__init__ guards it, even though the exact same class of precondition (v3 with ns<3, missing extension) is guarded there. maxnb=4 works correctly on v2 (verified), so the boundary is sharp.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/ed-and-operators/p3c.py 2 ; echo "exit=$?"
# the script builds Bosonic_Chain(3,maxnb=[3,3,3]), calls setup_cpp(version=2),
# then gs_energy(mode="DMRG"); run it with argument 3 / python / ED for the others
```

Observed:

```
--- maxnb=[3,3,3], itensor_version=2
/bin/bash: line 1: 2273354 Aborted                 (core dumped) ... python3 p3c.py 2
exit=134
From line 112, file get_sites.h

SpinX cannot read index of size 

SpinX cannot read index of size 

--- same chain, other backends (all fine)
ED   E0=-5.248282665570
3       E0=-5.248282665570
python  E0=-5.248282665570

--- maxnb=[4,4,4], itensor_version=2 (the one dimension v2 accepts)
ED   E0=-7.469145256243
2       E0=-7.469145256243
```

**Expected**: Either a Python-level `ValueError` naming the unsupported combination (what mode.py does for every other backend precondition), or the documented silent fall-back to the 4-level site. An uncatchable abort that kills the interpreter is neither, and it contradicts docs/user_guide.md:124-127 and bosonchain.py:18-19, both of which promise v2 "understands the single fixed 4-level boson site regardless of what maxnb requests".

**Suggested fix**: Guard it on the Python side, where it can raise: in `Bosonic_Chain.__init__` (and in `Many_Body_Chain.setup_cpp`, which can switch an existing chain to v2 after construction), raise ValueError when `itensor_version==2` and any `maxnb[i] != 4`, naming itensor_version=3/"python"/mode="ED" as the working alternatives -- the same pattern `Mixed_Spin_Fermion_Chain.__init__` (mixedchain.py:86-92) already uses to reject v2. Then correct docs/user_guide.md:124-127 from "only understand the single fixed 4-level boson site" to "reject a non-default maxnb". The identical claim is made there about the Julia backend; I could not test it (julia_live was out of bounds for this audit) but it is worth checking for the same discrepancy.

**Reviewer (CONFIRMED)**: Re-ran p3c.py 2 thread-pinned: "From line 112, file get_sites.h / SpinX cannot read index of size" then "Aborted (core dumped)", exit=134. Boundary verified sharp in the other direction with my own script: Bosonic_Chain(3,maxnb=[4,4,4]) + setup_cpp(version=2) runs to completion (ED -6.5015138, v2 DMRG -6.4754524, exit=0).

Code confirms the mechanism: src/dmrgpy/mpscpp2/get_sites.h:110 accepts only the literal `nm==104`, everything else falls to `Error(format("SpinX cannot read index of size "))`, and ITensor's Error() calls abort() so nothing is catchable from Python. This is dmrgpy's own header, not vendored ITensor. No guard exists in mode.py or Bosonic_Chain.__init__ even though mixedchain.py:86-92 already uses exactly the Python-side ValueError pattern for rejecting v2.

The documentation discrepancy is real: docs/user_guide.md:124-127 says v2 "still only understand the single fixed 4-level boson site regardless of what maxnb requests", which promises a finite (if wrong) answer, not a core dump. Not out of scope: audit finding #21 establishes uncatchable SIGABRT on unvalidated input as an in-scope class, and this is a different trigger from that one (maxm<=0), so it is not already recorded.

One stale-comment note the hunter did not flag: src/dmrgpy/bosonchain.py:15-20's comment lists "pyitensor" alongside mpscpp2 as only understanding code 104, but pyitensor/sites/boson.py::get_boson_site() builds an arbitrary-dimension site and p3c.py python returns the correct -5.248282665570 -- that half of the comment is wrong too.

**Reviewer on severity**: Medium is the ceiling, not a floor: reaching this needs an explicit itensor_version=2 / setup_cpp(version=2) on a combination the user guide already tells you to avoid (the default backend is 3). The finding is real because the failure mode is an uncatchable process abort and the docs affirmatively describe a different one.

<details><summary>Corroborating report from lens <code>recent-commits</code>: Bosonic_Chain with a non-default maxnb aborts the whole process (SIGABRT, uncatchable) on itensor_version=2, where the user guide says it merely degrades to a 4-level site</summary>

**Where**: `src/dmrgpy/mpscpp2/get_sites.h:112; docs/user_guide.md:119-129 (text rewritten by in-window commit 3593a2c)`

`Bosonic_Chain(n, maxnb=[...])` encodes the local dimension as site-type code 100+dim. `mpscpp2/get_sites.h`'s `SpinX(std::vector<int>)` only knows the literal code 104, and its `else` branch calls ITensor's `Error(...)`, which calls `abort()` -- so the process dies with a core dump that no Python `try/except` can intercept, the moment `setup_cpp(version=2)` builds the session. The user guide (rewritten by 3593a2c, in window) says v2 and Julia "still only understand the single fixed 4-level boson site regardless of what maxnb requests", which reads as silently-degraded-but-running and tells users to prefer v3 "for DMRG/ED results to actually agree" -- i.e. it describes a numerical caveat where the reality is process termination. The error message itself is also broken: `format("SpinX cannot read index of size ")` has no `%d`, so the abort does not even print the offending code (visible in the output below as a trailing blank).

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/recent-commits && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 -u -c "
from dmrgpy import bosonchain
bc = bosonchain.Bosonic_Chain(3,maxnb=[6,6,6])
print('constructed, backend', bc.itensor_version)
try:
    bc.setup_cpp(version=2); print('setup_cpp(2) returned')
except Exception as e:
    print('caught',type(e).__name__,e)
print('still alive')"
```

Observed:

```
constructed, backend 3
From line 112, file get_sites.h

SpinX cannot read index of size 

SpinX cannot read index of size 
timeout: the monitored command dumped core

# neither 'caught ...' nor 'still alive' is ever reached
```

**Expected**: Either a catchable Python exception naming the unsupported site-type code and the backend (the pattern CLAUDE.md says the v3 sector code adopted precisely because ITensor's `Error()` is uncatchable), or -- matching what the user guide currently promises -- construction of the fixed 4-level boson site with a warning. Silently taking down the caller's interpreter is neither.

**Suggested fix**: Replace `Error(format("SpinX cannot read index of size "))` in `mpscpp2/get_sites.h:112` with a `throw std::invalid_argument` that includes the offending code and the site index (and fix the missing `%d`), so pybind11 surfaces it as a Python exception. Cheaper and backend-agnostic: validate the site-type codes in `sites.py::initialize` / `Many_Body_Chain.setup_cpp()` against what the requested backend supports, before constructing the session. Either way, correct the user_guide.md paragraph to say v2 rejects a non-default maxnb rather than degrading it.

**Reviewer (CONFIRMED)**: Re-ran the repro: Bosonic_Chain(3, maxnb=[6,6,6]) builds site codes [106,106,106]; setup_cpp(version=2) prints 'From line 112, file get_sites.h / SpinX cannot read index of size ' twice and the process dies with exit code 134 (SIGABRT, core dumped). Neither the `except Exception` handler nor the trailing print is ever reached, so it is genuinely uncatchable from Python, and the missing %d in format("SpinX cannot read index of size ") is visible as the trailing blank in my output too.

Not out of scope as vendored ITensor: mpscpp2/get_sites.h is dmrgpy's own site dispatch (only 2,0,1,3,4,5,6,104,-2,-3 are handled; the else calls ITensor's Error()). CLAUDE.md itself records that ITensor's Error() calls abort() and that the v3 sector code therefore throws std::invalid_argument instead -- the repo's own established shape for exactly this, so the finding's fix is not an invention. The doc-mismatch half is also real and in-window: the user_guide.md paragraph ('still only understand the single fixed 4-level boson site regardless of what maxnb requests ... for DMRG/ED results to actually agree') was written by 3593a2c (2026-09-06) and reads as numerical degradation, not process death.

**Reviewer on severity**: low rather than medium. Reaching it requires two deliberate non-default steps at once (a non-default maxnb plus an explicit setup_cpp(version=2) / itensor_version attribute switch away from the v3 default) onto a backend the user guide already steers away from -- and per finding #3 the constructor kwarg route is not even available on this class. I would not go below low: an uncatchable SIGABRT in a library call is worse than a wrong number.

</details>


---


### 19. Bosonic_Chain and Parafermionic_Chain never cache their ED object, so every mode="ED" call rebuilds the full sparse Hamiltonian and re-diagonalizes: 37x / 585x slower for a per-site vev sweep

`optimization` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `ed-and-operators`

**Where**: `src/dmrgpy/bosonchain.py:47-56 (Bosonic_Chain.get_ED_obj) and :110-118 (SpinBoson_Chain.get_ED_obj); src/dmrgpy/parafermionchain.py:35-37 (Parafermionic_Chain.get_ED_obj); contrast src/dmrgpy/fermionchain.py:138-153`

`Many_Body_Chain` has a caching protocol for the ED backend: `has_ED_obj`/`ED_obj`, invalidated by `restart()` (manybodychain.py:706). `Fermionic_Chain.get_ED_obj` implements it correctly.

`Bosonic_Chain.get_ED_obj` writes `self.ed_obj = out` (lower case, a different attribute from the `self.ED_obj` the protocol uses) and **never sets `self.has_ED_obj = True`**, so its `if not self.has_ED_obj` guard is always true and the object is rebuilt on every single call. `Parafermionic_Chain.get_ED_obj` has no cache at all -- it unconditionally constructs a fresh `parafermion.Parafermion_Chain(self)`.

Rebuilding the object throws away not just the sparse operator dictionary (one2many over every site x every operator name) but also the ED object's own `computed_gs`/`Diagonalized_Hamiltonian` caches, so every `vev(..., mode="ED")`, `gs_energy(mode="ED")`, `get_excited(mode="ED")` pays a full ground-state solve. A per-site observable sweep -- `get_density()`, `get_density_fluctuation()`, any correlator sweep -- therefore costs ns full ED solves instead of one. This is pure waste: the answers are identical, and the fermionic chain already demonstrates the fix.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/ed-and-operators/cache.py
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/ed-and-operators/cache2.py
# cache2.py: min of 5 repeats, one full per-site <N_i> (resp. <Tau_i>) sweep
# via the public API vs the same sweep on a single ED object fetched once
```

Observed:

```
cache.py:
Bosonic        get_ED_obj() identical across calls: False
Fermionic      get_ED_obj() identical across calls: True
Parafermionic  get_ED_obj() identical across calls: False

cache2.py (min of 5, threads pinned):
Bosonic(4 sites,dim 256)        4 vevs via public API = 0.6896 s ;  same via one cached ED obj = 0.0188 s  (36.7x)
Parafermion(5 sites,dim 243)    5 vevs via public API = 0.8666 s ;  same via one cached ED obj = 0.0015 s  (585.3x)
```

**Expected**: `sc.get_ED_obj() is sc.get_ED_obj()` should be True for every chain class between `restart()` calls, as it already is for `Fermionic_Chain` -- that is what the `has_ED_obj`/`ED_obj` protocol in Many_Body_Chain exists for. A per-site vev sweep should cost one ED solve plus ns cheap operator assemblies (the 0.0188 s / 0.0015 s column), not ns full solves.

**Suggested fix**: bosonchain.py:47-56 -- store into `self.ED_obj`, set `self.has_ED_obj = True`, and read `self.ED_obj` on the cached branch (the current code writes `self.ed_obj` and never flips the flag, so both halves are dead). parafermionchain.py:35-37 -- add the same guard: `if self.has_ED_obj: return self.ED_obj`, then build, `self._apply_sector_to_ed(obj)` for consistency with the other classes, store, set the flag. `restart()` already clears `has_ED_obj`, so correctness after `set_hamiltonian` is preserved (verified: reusing a chain across two different Hamiltonians gives bit-identical energies to a fresh chain on all three classes). SpinBoson_Chain.get_ED_obj has the same `self.ed_obj` typo but is an unfinished stub.

**Reviewer (CONFIRMED)**: Both repros re-ran thread-pinned. cache.py: get_ED_obj() identical across calls = False (Bosonic), True (Fermionic), False (Parafermionic); 1 vev vs 5 vevs scales 4.62x (Bosonic) and 5.63x (Parafermionic), i.e. every single vev pays a full rebuild+solve. cache2.py on my box: Bosonic(4 sites, dim 256) 0.3255 s vs 0.0086 s (37.9x), Parafermion(5 sites, dim 243) 0.4246 s vs 0.0007 s (577.1x) -- same ratios as reported, absolute times faster.

Code reading confirms both mechanisms exactly as described: bosonchain.py:47-56 stores `self.ed_obj` (lower case) and never sets `self.has_ED_obj`, so the `else: return self.ed_obj` branch is dead code and the `if not self.has_ED_obj` guard is permanently true; parafermionchain.py:36-38 has no cache at all. fermionchain.py:138-153 implements the protocol correctly (`self.ED_obj` + `self.has_ED_obj = True`), and manybodychain.py:238-239 declares both fields.

I checked the one thing that could have made the current shape the correct one -- stale-cache risk -- and it does not hold: set_hamiltonian() calls restart() by default (manybodychain.py:660,684) and set_conserved_sector() calls restart() at its end (manybodychain.py:424), and restart() clears has_ED_obj (:706). So caching is safe by exactly the invalidation path Fermionic_Chain already relies on.

Minor overstatement worth noting: the "cached" column in cache2.py pre-warms get_gs_array() outside the timer, so 37x/577x is the full sweep-vs-sweep gap rather than a per-call saving; the substantive claim (ns full ED solves instead of one) is carried by cache.py's linear 1-to-5 scaling and is correct.


---


### 20. idmrg._dominant_fixed_point's ARPACK matvec applies materialised chi^4 transfer tensors (O(chi^4)/iteration, chi^4 memory) although _apply_site_transfer in the same file does it in O(chi^3 d)

`optimization` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `pyitensor-performance`

**Where**: `src/dmrgpy/pyitensor/idmrg.py:2344-2420 (_dominant_fixed_point's matvec), :1985-2010 (_transfer_matrices), :2278-2283 (_CorrelatorEnv.__init__); reached from infinitechain.vev/correlator with gs_method="idmrg"`

`_CorrelatorEnv.__init__` unconditionally builds `Es = _transfer_matrices(cell, n_cell)` — one (chi,chi,chi,chi) array per cell site — and hands them to `_dominant_fixed_point`, whose ARPACK matvec applies each with `_apply_transfer` (a chi^2 x chi^2 gemv, O(chi^4) work and chi^4 bytes streamed per site per iteration). The module already contains `_apply_site_transfer(A, M, rho)`, whose own docstring states it computes 'the same quantity `_apply_transfer(_op_transfer_mat(...), rho)` computes, re-associated ... both halves O(chi^3 d)' and which is exact rather than approximate — but it is used only by the operator-string walk in `two_point_correlator`, not by the fixed-point solves that dominate the env build. Consequence: the first `vev`/`correlator` on a `gs_method="idmrg"` chain scales as ~chi^4 (0.69 s at maxm=32, 4.50 s at 48, 11.55 s at 64, 540 MB peak allocation at 64 — one E4 tensor is 268 MB at chi=64 and 1.36 GB at chi=96), with 95% of that inside `_apply_transfer`/`_apply_transfer_from_left`. Solving the identical eigenproblem matrix-free over the cell tensors, with the same ARPACK settings (k=2, tol=0, ncv=40, same v0), is 9.2x faster at chi=48 and 20.1x at chi=64 per eigenproblem — and the env solves two of them (right and left families) plus the `_transfer_matrices` build, none of which would be needed.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s3_idmrg.py     # maxm sweep of growth / first vev / later vev / correlator
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s4_env.py 48     # cProfile of the first vev
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s6_env_win.py 64   # same eigenproblem, matrix-free
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s5_micro.py   # per-application micro-benchmark
```

Observed:

```
maxm=  8  growth=0.172s  vev1(env build)=0.011s  vev2=0.0001s  corr(r=3)=0.0003s
maxm= 16  growth=0.367s  vev1(env build)=0.038s  vev2=0.0001s  corr(r=3)=0.0004s
maxm= 24  growth=0.725s  vev1(env build)=0.130s  vev2=0.0002s  corr(r=3)=0.0007s
maxm= 32  growth=1.722s  vev1(env build)=0.716s  vev2=0.0003s  corr(r=3)=0.0011s

### vev first call, maxm= 48   (total 2.983 s)
   _dominant_fixed_point   ncalls=2     cum=2.900
   _apply_transfer         ncalls=307   tot=1.512
   _apply_transfer_from_left ncalls=305 tot=1.329

maxm=48  growth=6.34s  FIRST vev (env build)=4.496s  peak-alloc=171.8 MB
   tiled cell: 2 sites, chi=48 ; one E4 tensor = 84.9 MB
   matrix-free right fixed point alone: 0.151s (eta=1.000591941722)
   shipped _dominant_right_fixed_point:  1.400s (eta=1.000591941722)  -> 9.2x
maxm=64  FIRST vev (env build)=11.549s  peak-alloc=540.2 MB
   tiled cell: 2 sites, chi=64 ; one E4 tensor = 268.4 MB
   matrix-free right fixed point alone: 0.277s (eta=1.000728272866)
   shipped _dominant_right_fixed_point:  5.569s (eta=1.000728272866)  -> 20.1x

 chi   d   E4-route(s)   site-route(s)   speedup   max|diff|
  32   2     0.000836       0.000107      7.8x   1.01e-16  (E4 mem 16.8 MB)
  64   2     0.014850       0.000599     24.8x   6.07e-17  (E4 mem 268.4 MB)
  96   2     0.074505       0.002107     35.4x   5.21e-17  (E4 mem 1359.0 MB)
```

**Expected**: The ARPACK branch exists precisely to avoid the O(chi^6) dense route, and its own comment says 'the iterative matvec is O(n_uc*chi^4)'. Given `_apply_site_transfer` in the same file, that matvec can be O(n_uc*chi^3*d), which is what an ITensorInfiniteMPS-style implementation does: the fixed-point solve should never materialise a rank-4 transfer tensor. The dominant eigenvalue returned by the two routes matches to 12 printed digits and the per-application outputs to ~1e-16, so this is a pure contraction-order change.

**Suggested fix**: Give `_dominant_fixed_point` the cell tensors (and, where present, the per-site operator matrices) rather than pre-built Es, and let its `matvec` loop `_apply_site_transfer(A, M, X)` over them. Keep `_transfer_matrices` only for the `n <= _DENSE_EIG_MAX` dense fallback and for `_compose`; above the threshold, `_CorrelatorEnv` then never allocates a chi^4 array at all (O(chi^2 d) instead of 540 MB at chi=64). `_all_right_fixed_points`/`_all_left_fixed_points`' own cheap propagation steps convert the same way.

**Reviewer (CONFIRMED)**: Reproduced all three legs. s4_env.py 48: the first vev costs 3.536 s, of which _apply_transfer 1.796 s + _apply_transfer_from_left 1.531 s = 94% (hunter claimed 95%), inside 2 _dominant_fixed_point calls / 344 ARPACK iterations. s6_env_win.py: at maxm=48 the shipped _dominant_right_fixed_point takes 1.416 s vs 0.084 s matrix-free (16.8x) with eta identical to all 12 printed digits (1.000591941719); at maxm=64, 4.042 s vs 0.210 s (19.3x), eta 1.000728333741 both ways, peak alloc 540.2 MB with one E4 tensor at 268.4 MB. s5_micro.py per-application: 10.9x at chi=32, 25.2x at 48, 33.1x at 64, 48.0x at 96, max|diff| ~1e-16 — a pure contraction-order change, exact as _apply_site_transfer's own docstring states. s3_idmrg.py sweep reproduced: vev1 0.007/0.024/0.119/0.599 s at maxm 8/16/24/32 against growth 0.105/0.221/0.480/1.222 s. The ARPACK comment's own 'O(n_uc*chi^4)' is an acknowledged cost, not a design decision that rules the O(chi^3 d) route out — the same file ships the cheaper helper and calls it exact. Nothing in the audit/known-issue/plan docs covers it. Two fix constraints I found (neither refutes the finding, both bound the suggested_fix): _apply_site_transfer(A,M,rho) uses A for both ket and bra, while _transfer_matrices supports a mixed bra_list that imps_overlap (idmrg.py:2692) relies on, so a matrix-free path needs an explicit bra argument; and _all_right_fixed_points/_all_left_fixed_points' own propagation steps consume Es too, so the cell has to be threaded through them as well.


---


### 21. user_guide's MPS/operator-algebra table claims all nine calls take mode=/**kwargs; only one of the nine actually honours it

`hole` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `docs/user_guide.md:202-214 and docs/user_guide.tex:235-253 ('Each takes the same mode=/**kwargs as the rest of the API'); code at src/dmrgpy/manybodychain.py:824-856 and src/dmrgpy/mpsalgebra.py:50,90,110,119,179`

The §2 'MPS and operator algebra' table lists nine primitives and prefaces them with 'Each takes the same `mode=`/`**kwargs` as the rest of the API' (identical sentence in the .tex). Executed against a 4-site Heisenberg chain, three of the nine raise TypeError on `mode=`: `scale_mps` (manybodychain.py:849 takes only `(self,x,wf)` — no kwargs at all), `operator_norm` and `is_zero_operator` (both land in mpsalgebra.py:179 `operator_norm(self,op,ntries=5,simplify=True)`, which has neither `mode` nor `**kwargs`). Three more accept `mode=` and discard it: `applyoperator`/`summps`/`applyinverse` (mpsalgebra.py:90/110/99) take `**kwargs`, never read them, and decide DMRG-vs-ED from `type(wf)` instead — so `mode="ED"` on an MPS argument silently runs DMRG. Two (`overlap`, `aMb`) route on `mode=` but then hand the MPS objects to the ED object, which just calls `.dot`, i.e. an MPS inner product again. Only `trace` behaves as the sentence promises. This is documentation advertising a uniform kwarg surface that does not exist, which is the same class as the audit's 'kwargs with no consumer'.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 5 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/algebra_mode.py
```

Observed:

```
chain mode: None itensor_version: 3
OK    overlap(a,b,mode='ED') -> complex (1.0000000000000009+0j)
OK    aMb(a,A,b,mode='ED') -> complex (6.854660068722485e-13+0j)
OK    applyoperator(A,wf,mode='ED') -> MPS 
OK    summps(wf,wf,mode='ED') -> MPS 
FAIL  scale_mps(2.0,wf,mode='ED') -> TypeError: Many_Body_Chain.scale_mps() got an unexpected keyword argument 'mode'
FAIL  operator_norm(A,mode='ED') -> TypeError: operator_norm() got an unexpected keyword argument 'mode'
FAIL  is_zero_operator(A,mode='ED') -> TypeError: operator_norm() got an unexpected keyword argument 'mode'
OK    trace(A,mode='ED') -> complex128 0j
OK    exponential(A,wf,mode='ED') -> MPS
```

**Expected**: Either every call in the table accepts and honours `mode=` (what both documents state), or the sentence names the subset that does. A reader following the guide and writing `sc.operator_norm(h, mode="ED")` gets a TypeError four frames deep with no hint that the guide is wrong; a reader writing `sc.applyoperator(A, wf, mode="ED")` gets no error at all and the wrong backend.

**Suggested fix**: Cheapest correct fix is on the code side, since three of the four dispatchers already have the branch: give `scale_mps` and `operator_norm`/`is_zero_operator` a `**kwargs` passthrough, and make `applyoperator`/`summps`/`applyinverse` prefer an explicit `mode=` (falling back to the current `type(wf)` inference when none is given) rather than dropping it. If that is more churn than wanted, replace the blanket sentence in both .md and .tex with 'overlap, aMb, exponential and trace take mode=; the remainder infer the backend from the wavefunction type'.

**Reviewer (CONFIRMED)**: I executed all eleven table rows myself rather than rerunning the hunter's nine (4-site Heisenberg, real MPS ground state). Four raise TypeError on mode=: applyinverse ('applyinverse_dmrg() got an unexpected keyword argument mode'), scale_mps, operator_norm, is_zero_operator. Note this CORRECTS the hunter, who listed applyinverse among the silent-ignore group -- it actually raises, so the count of hard failures is four, not three. Of the rest: applyoperator and summps take **kwargs and decide on type(wf); exponential's signature has mode='DMRG' but line 6 immediately overwrites it with mode = wf.mode, so it is accept-and-discard too; overlap and aMb do dispatch on mode=; inverse_trace honours mode= genuinely; trace forwards **kwargs into toMPO. So the guide's blanket sentence 'Each takes the same mode=/**kwargs as the rest of the API' (user_guide.md:202 and user_guide.tex:235) is false for at least six of eleven rows and hard-fails for four. Could not explain it away: the type(wf) inference may well be the intended design, but that does not make the documented sentence true, and a reader typing sc.operator_norm(h, mode='ED') gets a TypeError with no hint. Nothing about this is in the audit or known_issue docs.


---


### 22. A mistyped ctmode/dmmode/fpmode string produces 'RuntimeError: No active exception to reraise' instead of naming the valid options

`hole` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `src/dmrgpy/entropytk/correlationentropy.py:261 (ctmode), :96 (dmmode), :84 (non-fermion chain); src/dmrgpy/fermionicparity.py:73 (fpmode); documented enumerations at docs/user_guide.md:1040-1073 (five ctmode values), :1138 (ctmode=None resolver), :2340-2341 (fpmode)`

Three dispatchers over a documented string enumeration end their if/elif chain in a bare `raise` with no active exception, so Python turns a user typo into `RuntimeError: No active exception to reraise` — a message that names neither the argument nor the accepted values. `get_four_correlation_tensor(ctmode=...)` is the worst of the three because the user guide makes an explicit promise about it ('Passing a `ctmode` explicitly is still a hard request — it raises rather than silently falling back if that method isn't available'): that promise is kept for a valid-but-unavailable mode (the per-mode helpers raise a clear ValueError), and broken for a misspelling. `get_correlation_matrix(dmmode=...)` additionally `print`s 'fasst not recognized' to stdout before raising, so the diagnostic exists but is not in the exception. This is exactly documentation.md §4.10's 'an else serving both "unsupported here" and "you typo'd it"' class, which audit finding #12 fixed at `get_dynamical_correlator(name=)` and left standing here. I enumerated all 40 bare-`raise` sites under src/dmrgpy (excluding ITensor/) and executed these three, chosen because each gates a string option the user guide documents by name.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 5 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/bareraise.py   # 4-site Fermionic_Chain; calls get_four_correlation_tensor(ctmode='sweeep'), get_correlation_matrix(dmmode='fasst'), get_fermionic_parity(fpmode='fulll')
```

Observed:

```
--- get_four_correlation_tensor(ctmode='typo') ---
   RuntimeError: No active exception to reraise
--- get_correlation_matrix(dmmode='typo') ---
fasst not recognized
   RuntimeError: No active exception to reraise
--- get_fermionic_parity(fpmode='typo') ---
   RuntimeError: No active exception to reraise
```

**Expected**: A ValueError naming the argument and the accepted values, e.g. ValueError("get_four_correlation_tensor: ctmode='sweeep' is not recognized; expected one of 'explicit', 'full', 'sweep', 'fold', 'batched', or None to auto-select") — the shape already used elsewhere in the same file (get_four_correlation_tensor_fold/_batched both raise informative ValueErrors) and the shape audit #12 applied when it fixed the identical pattern in get_dynamical_correlator.

**Suggested fix**: Replace the three bare `raise` statements with ValueErrors listing the valid strings. While in correlationentropy.py, the neighbouring bare `raise` at :84 (reached when `operators is None` on a non-fermionic chain, after a `print("Unrecognized type",type(self))`) deserves the same treatment — it is the same defect one branch away.

**Reviewer (CONFIRMED)**: Reproduced all three with my own script on a 4-site Fermionic_Chain ground state: get_four_correlation_tensor(ctmode='sweeep') -> RuntimeError: No active exception to reraise; get_correlation_matrix(dmmode='fasst') -> prints 'fasst not recognized' to stdout then the same RuntimeError; get_fermionic_parity(fpmode='fulll') -> same RuntimeError. Read the three sites and all three are bare `raise` outside any except block (correlationentropy.py:261 and :96, fermionicparity.py:73), plus a fourth one branch away at correlationentropy.py:84. Scope check: the audit does record this exact pattern class (documentation.md 4.10, and audit finding #12 fixed it at get_dynamical_correlator), and audit #11 touches dmmode, but for a different defect (dmmode='fast' defaulting into a ValueError under a conserved sector, marked FIXED). These three typo-path sites are not recorded anywhere, so the finding is not out of scope -- only its taxonomy is already named.

**Reviewer on severity**: low rather than medium. No wrong number and no silent misbehaviour: the call does raise, so the user guide's stated promise ('a ctmode explicitly is a hard request -- it raises rather than silently falling back') is literally kept. The defect is message quality on a typo path, which is a papercut, not a correctness or a doc-drift failure.


---


### 23. Spin_Chain.get_hamiltonian() raises TypeError: 'int' object is not iterable for any chain without set_hamiltonian() — the dead self.exchange path 43d1a35 fixed in meanfield.py but left here

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `src/dmrgpy/spinchain.py:136-152 (Spin_Chain.get_hamiltonian, the `else: # conventional way` branch, iterating `self.exchange` at :144 and `self.fields` at :147); dead attributes initialized at src/dmrgpy/manybodychain.py:76-77`

`Spin_Chain.get_hamiltonian()` returns `self.hamiltonian` when one was set, and otherwise falls back to a 'conventional way' branch that builds the operator from `self.exchange` and `self.fields`. Those two attributes were populated by `Spin_Chain.set_exchange()`, a builder that has been removed; `Many_Body_Chain.__init__` now leaves both as the integer 0, so the fallback cannot run — it dies on `for c in self.exchange`. This is the same corpse commit 43d1a35 diagnosed and repaired in `meanfield.py` ('every call raised TypeError: int object is not iterable, on any chain, while the user guide documented it as working'); the fix was not applied to this sibling site, nor to `pychainwrapper.py::old2ampo` (which reads the same two attributes and additionally references a bare undefined name `fields` at :18-19, so it would NameError even if `self.exchange` were a list). `gs_energy_fluctuation()` is a public method that calls `get_hamiltonian()` unconditionally and so inherits the failure.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 5 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/getham.py   # Spin_Chain(['S=1/2']*4) with no set_hamiltonian(), then get_hamiltonian() and gs_energy_fluctuation()
```

Observed:

```
hamiltonian attr: None
exchange: 0  fields: 0
get_hamiltonian FAIL: TypeError: 'int' object is not iterable
gs_energy_fluctuation FAIL: TypeError: 'int' object is not iterable
```

**Expected**: Either the documented behaviour (return the Hamiltonian as a MultiOperator) or, since there is no builder left that can populate `self.exchange`, a clear ValueError of the kind meanfield.py now raises: 'this chain has no Hamiltonian yet; call set_hamiltonian() first'. A bare TypeError about integers is not actionable for a caller who never touched `exchange`.

**Suggested fix**: Delete the unreachable `else` branch of `Spin_Chain.get_hamiltonian` and replace it with `raise ValueError("this chain has no Hamiltonian yet; call set_hamiltonian() first")`, mirroring `meanfield.decompose_spin_hamiltonian`. Then drop `self.exchange`/`self.fields` from `Many_Body_Chain.__init__` (manybodychain.py:76-77) and remove `pychainwrapper.old2ampo` plus the `h = h + self.exchange` line in `update_hamiltonian` (manybodychain.py:819), which are the remaining readers.

**Reviewer (CONFIRMED)**: Reproduced: Spin_Chain(['S=1/2']*4) with no set_hamiltonian() -> hamiltonian attr None, exchange 0, fields 0, and both get_hamiltonian() and gs_energy_fluctuation() die with TypeError: 'int' object is not iterable. Verified the branch is genuinely dead, not merely unused: grep -rn 'set_exchange' over src/dmrgpy returns only comments (spinchain.py:66-68, :89-93 and meanfield.py:29 all record that the builder was removed and that self.exchange is now permanently the integer 0) -- there is no surviving code path that can populate it. Scope check that could have flipped this: git show 43d1a35 -s --format=%B is entirely about meanfield.py and says nothing about Spin_Chain.get_hamiltonian or pychainwrapper.old2ampo, so this sibling site is not a recorded deferral. One correction to the hunter's write-up: manybodychain.py:819's `h = h + self.exchange` is harmless (0 is the additive identity for a MultiOperator sum), so only spinchain.py:136-152 and the caller-less pychainwrapper.old2ampo are actually broken.

**Reviewer on severity**: low rather than medium. It is only reachable by calling get_hamiltonian() before set_hamiltonian(), which is a caller error either way; the outcome is an unhelpful exception rather than a wrong number, and nothing in the user guide documents get_hamiltonian() as a way to obtain a Hamiltonian that was never set (unlike meanfield, which the guide did document as working).


---


### 24. examples/boson_models/v2_VS_v3_boson/main.py -- an assert-carrying regression example is non-deterministic on a clean tree (2 of 5 runs fail for the reviewer, 5 of 7 for the hunter), from DMRG under-convergence that affects v2 and v3 symmetrically at its default nsweeps=15/maxm=30

> Retitled by the orchestrator: the reviewer confirmed the observation but disputed the
> hunter's diagnosis in the original title, which was: *examples/boson_models/v2_VS_v3_boson/main.py — an assert-carrying regression example fails 5 of 7 runs, with ITensor v2 landing up to 0.025 above the ED/pure-Python energy*

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `examples/boson_models/v2_VS_v3_boson/main.py:70 (the U-sweep assert, tol=1e-2); the divergent backend is itensor_version=2 on Bosonic_Chain`

This is one of the v2_VS_v3_* cross-backend regression scripts CLAUDE.md points at as 'the fastest way to check a mpscpp3 change didn't diverge from mpscpp2's numerics'. Run seven times on a clean tree with threads pinned, its U-sweep assert failed five times. The printed table identifies which side is wrong: ED and the pure-Python backend agree with each other to 5e-11 at every U, itensor_version=3 tracks them to 1e-3..4e-3 (inside the 1e-2 tolerance), and itensor_version=2 is the one that drifts — e.g. U=0.5 v2=-0.5788227 against ED=-0.6033631 (0.0245 above the true ground state), U=0.2 v2=-2.2601163 against ED=-2.2759299 (0.0158 above). Because the energies are always *above* ED's, this is v2's boson DMRG failing to converge rather than a tolerance being merely tight, and because both v2 and v3 start from unseeded random-ish states the failing U point moves between runs. Net effect: a regression script whose assert is noise, so a real future v2/v3 divergence in the boson backend would not be distinguishable from its background failure rate.

Repro:

```bash
cd /home/joselado/Documents/programs/dmrgpy/examples/boson_models/v2_VS_v3_boson && for i in 1 2 3 4; do echo "##### RUN $i"; MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src MPLBACKEND=Agg taskset -c 5 timeout 200 python3 main.py 2>&1 | grep -E '^U =|AssertionError|TEST PASSED'; done
```

Observed:

```
# batch run (exit status recorded by the driver script): EXIT=1 boson_models/v2_VS_v3_boson
U = 0.5   v2 = -0.5788227193947524   v3 = -0.6033628792737404   ED = -0.6033630654449267   python = -0.603363065444918
AssertionError: U=0.5: v2 vs v3 disagree by 0.0245402 (tol=0.01)

# four further runs:
##### RUN 1
U = 0.2   v2 = -2.2601162902491643   v3 = -2.2759081611049825   ED = -2.2759298805426926   python = -2.275929880538307
AssertionError: U=0.2: v2 vs v3 disagree by 0.0157919 (tol=0.01)
##### RUN 2   (all five U points passed)
##### RUN 3
AssertionError: U=0.1: v2 vs v3 disagree by 0.0102322 (tol=0.01)
##### RUN 4
AssertionError: U=0.2: v2 vs ED disagree by 0.0151103 (tol=0.01)
# tally over all seven executions: 5 failures, 2 passes
```

**Expected**: An assert-carrying example is a pass/fail regression guard, so it must pass deterministically on a clean tree. Since ED and 'python' agree to 5e-11 here, the exact answer is not in doubt: v2's Bosonic_Chain DMRG should reach it too. The script's own header calls itself a 'Regression test'.

**Suggested fix**: Diagnose the v2 side first rather than widening the tolerance — the errors are one-sided (always above ED), which points at v2's boson DMRG stopping short; raising nsweeps/maxm/noise inside get_energy() and re-measuring would tell you whether it is convergence or something structural in mpscpp2's 4-level boson site. If it is convergence, pin nsweeps/maxm explicitly in the script (cheap at n=6) and keep tol=1e-2; if v2 cannot reach 1e-2 on this model at any reasonable settings, that is itself worth recording in the script's comment block rather than leaving an assert that fires most of the time.

**Reviewer (CONFIRMED)**: CONFIRMED on the observable claim, with the diagnosis in the title REFUTED. First ruled out a stale extension (my MEMORY.md's usual explanation): mpscpp2/_dmrgcpp*.so is dated 2026-09-05 16:54 against newest mpscpp2 source 2026-08-30, and mpscpp3's .so is 2026-09-05 21:34 against chain_session.h 2026-09-05 21:33 -- both fresh, so the flakiness is not a build artefact. Ran the script as committed five times (pinned, MPLBACKEND=Agg): the U-sweep assert fired in runs 1 and 3 and passed in 2, 4, 5 -- 2/5, lower than the hunter's 5/7 but decisively non-deterministic on a clean tree. What I could not reproduce is the attribution. In BOTH of my failures the outlier was itensor_version=3, not 2: run 1 U=0.2 gave v3=-2.2616090 against ED=-2.2759299 and v2=-2.2756110; run 3 U=0.4 gave v3=-0.9427030 against ED=-0.9576023 and v2=-0.9575679. That directly contradicts 'v3 tracks them to 1e-3..4e-3 ... v2 is the one that drifts'. I then measured the mechanism: 4 repeats per point at U=0.2/0.4/0.5, worst error against ED at the script's defaults (nsweeps=15, maxm=30) was 6.6e-4/8.9e-4/2.3e-4 for v2 and 1.1e-3/1.4e-3/4.8e-4 for v3; at nsweeps=40/maxm=100 both collapse to 1.5e-5/1.1e-12/1.8e-7 (v2) and 5.2e-5/1.1e-12/3.2e-7 (v3). So this is DMRG under-convergence from a random start, affecting BOTH compiled backends symmetrically, with a tail that occasionally exceeds the script's tol=1e-2. The brief's 'tighten and it shrinks -> REFUTED' rule does not dispose of this one, because the claim is about the committed assert, not about a number being wrong: a regression guard that fires ~40% of the time on a clean tree cannot distinguish a future v2/v3 divergence from its own background rate, and the fix is two lines (pin nsweeps/maxm inside get_energy).

**Reviewer on severity**: Severity medium is defensible, but the title and description must be rewritten: the failing backend is not specifically v2, the '0.025 above' and '5 of 7' figures did not reproduce (I saw 0.014/0.015 and 2 of 5), and the cause is shared under-convergence at the script's default sweep parameters rather than anything structural in mpscpp2's boson site.


---


### 25. The user guide's SpinBoson_Chain constructor example uses a site label the code rejects, and the class silently ignores its own maxnb= argument

`hole` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `docs/user_guide.md:50 (SpinBoson_Chain(["boson","S=1/2",...])); code at src/dmrgpy/bosonchain.py:60-67 (get_site accepts only "B"/"B4") and :80-83 (__init__(self,sitesin,n=None,maxnb=None) — both n and maxnb are dead, the two lines that used them are commented out)`

Two drifts in one class, both introduced or left by the recent documentation pass (3593a2c lists SpinBoson_Chain under 'Undocumented API, all verified working before being written up'). (1) The §1 model table gives the constructor call literally as `SpinBoson_Chain(["boson","S=1/2",...])`. `bosonchain.get_site` accepts only `"B"` or `"B4"` for a boson location, so the documented call dies in the constructor with `RuntimeError: No active exception to reraise` — the bare-raise problem of finding 4 met at the first line a reader following the guide would type. The shipped example `examples/fermion_models/spinboson_chain/main.py` uses `"B"`, confirming which spelling is real. This is also an md/tex divergence: docs/user_guide.tex:66 covers the same class without giving any example call, so only the .md is wrong. (2) `SpinBoson_Chain.__init__` takes `n=None, maxnb=None` and consumes neither — the two lines that would have applied maxnb are commented out immediately below the signature. A caller asking for larger bosons gets the 4-level site silently: sites come back as [104,104] with maxnb=[8,8], where the sibling `Bosonic_Chain(2, maxnb=[8,8])` correctly gives [108,108]. The guide's own boson section documents maxnb threading through to the site type code 100+maxnb for `Bosonic_Chain`, so the two classes look interchangeable in this respect and are not.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 5 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/sbc2.py
```

Observed:

```
--- exactly as docs/user_guide.md line 50 writes it ---
FAIL: RuntimeError: No active exception to reraise
--- the label the code actually accepts ---
OK sites = [104, 2]
--- maxnb= is accepted and ignored (D0..D3 only, site code 104) ---
maxnb=[8,8] -> sites = [104, 104]  (104 == 4 levels)
Bosonic_Chain(2,maxnb=[8,8]) -> sites = [108, 108]
```

**Expected**: The documented snippet must construct: `SpinBoson_Chain(["B","S=1/2"])` works and gives sites [104,2]. And a constructor argument named `maxnb` should either change the local dimension (as it does on Bosonic_Chain, giving site code 100+maxnb) or not be in the signature at all — accepting it and returning 4-level sites is the 'kwarg with no consumer' shape, here with a silently wrong Hilbert space rather than a silently ignored option.

**Suggested fix**: In docs/user_guide.md:50 change `["boson","S=1/2",...]` to `["B","S=1/2",...]` (the .tex needs no change, but adding the same concrete example there would remove the divergence in the other direction). In bosonchain.py either wire `maxnb` through `get_site` the way `Bosonic_Chain` does (the site code is just 100+maxnb, and get_site already returns 104 for "B4") or delete `n=`/`maxnb=` from the signature so the call fails loudly. Fixing the bare `raise` in `get_site` at :66 to a ValueError naming the accepted labels would have made this a one-line diagnosis instead of a RuntimeError.

**Reviewer (CONFIRMED)**: Reproduced both halves. (1) SpinBoson_Chain(['boson','S=1/2']) exactly as user_guide.md:50 writes it dies with RuntimeError: No active exception to reraise, from the bare raise in bosonchain.get_site at :66, which accepts only 'B' and 'B4'; SpinBoson_Chain(['B','S=1/2']) works and gives sites [104, 2]. The shipped example examples/fermion_models/spinboson_chain/main.py uses the 'B' spelling, and user_guide.tex:66 gives no example call at all, so the md/tex divergence claim also holds. (2) SpinBoson_Chain(['B','B'], maxnb=[8,8]) returns sites [104, 104] with no maxnb attribute set, while Bosonic_Chain(2, maxnb=[8,8]) correctly gives [108, 108]; read bosonchain.py:80-83 and the two lines that would consume maxnb are commented out directly under the signature, so n= and maxnb= are dead parameters. This is a wrong Hilbert space handed back silently, which is worse than the kwarg-with-no-consumer shape documentation.md 4.10 names. Nothing in the audit, known_issue or ROADMAP files covers it.


---


### 26. get_excited_states(n=1, mode="DMRG") returns a plain list of energies while every other (n, mode) returns an ndarray, so get_gs_manifold(n=1) raises TypeError on DMRG but works on ED

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `wavefunction-consumers`

**Where**: `src/dmrgpy/excited.py:110-112 (the n==1 Hermitian branch) vs :101-107 (the sibling non-Hermitian n==1 branch, which returns np.array), consumer src/dmrgpy/groundstate.py:389`

`excited.get_excited_states` has two `n==1` short circuits. The non-Hermitian one at line 107 returns `(np.array([e0]), [self.wf0.copy()])`; the Hermitian one at line 112 returns `([e0], [w])` — a bare Python list. Every other path (`get_excited_states_dmrg` at :55, the purify branch, and the whole ED route) returns an ndarray, so `n=1` under `mode="DMRG"` is the single inconsistent case in the API. `groundstate.get_gs_manifold` (line 389) does `es[np.abs(es-e0)<tol]`, which on a list raises `TypeError: unsupported operand type(s) for -: 'list' and 'float'` — and it raises from *inside* the recursion, after the "Recalling with 2 states" print, so the failure appears at n=1 even when the caller asked for more. `Many_Body_Chain.get_gs_manifold(n=1)` is public and is reached indirectly by `fidelity.get_fidelity`'s default `fmode="derivative"` branch (fidelity.py:31) whenever its `n` resolves to 1. `mode="ED"` returns 1 state correctly on the same chain, so this is a DMRG-only break of a capability that exists and works.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/wavefunction-consumers/repro_f3.py   # standalone
```

Observed:

```
get_excited_states(n=1,mode=ED   ) -> type(es)=ndarray
get_excited_states(n=2,mode=ED   ) -> type(es)=ndarray
get_excited_states(n=1,mode=DMRG ) -> type(es)=list
get_excited_states(n=2,mode=DMRG ) -> type(es)=ndarray
Recalling with  2 states
get_gs_manifold(n=1,mode=ED   ) -> 1 states
get_gs_manifold(n=1,mode=DMRG ) RAISED TypeError("unsupported operand type(s) for -: 'list' and 'float'")
```

**Expected**: `get_excited_states` should return the same type for every (n, mode) — an ndarray of energies plus a list of wavefunctions — as its own non-Hermitian n==1 sibling three lines above already does, and `get_gs_manifold(n=1, mode="DMRG")` should return the one ground state, exactly as `mode="ED"` does on the same chain.

**Suggested fix**: In `src/dmrgpy/excited.py:112`, change `return ([e0],[w])` to `return (np.array([e0]),[w])`, matching line 107. Optionally also make `groundstate.get_gs_manifold` defensive with `es = np.array(es)` at its top.

**Reviewer (CONFIRMED)**: Reproduced exactly, line for line including the 'Recalling with 2 states' print ordering. Read the code: excited.py:107 returns (np.array([e0]),[...]) for the non-Hermitian n==1 short circuit and excited.py:112 returns ([e0],[w]) for the Hermitian one three lines below -- two sibling branches of the same function disagreeing on return type, which is a defect independent of any consumer. groundstate.get_gs_manifold:389 does es[np.abs(es-e0)<tol], which raises on a list; es[0] on the preceding line happens to work, which is why the failure surfaces one line late.

The only refutation angle I could construct is that get_gs_manifold(n=1) is outside the method's contract -- one state cannot establish degeneracy. That angle dies on the fact that mode='ED' answers the identical call correctly (it recurses to n=2 and returns 1 state), so the capability exists and only the DMRG route breaks. Not recorded in the audit doc, ROADMAP, or any known_issue file. The suggested one-line fix is correct and matches the sibling branch.


---


### 27. get_bond_entropy's guard is off by one (b>self.ns rather than b>=self.ns), so an out-of-range site index kills the whole process with an uncatchable SIGABRT on both C++ backends

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `dispatch-matrix`

**Where**: `src/dmrgpy/entropy.py:29-36 (compute_entropy_single); src/dmrgpy/entropy.py:38-42 (bond_entropy); src/dmrgpy/manybodychain.py:877-879`

`compute_entropy_single` validates `if b<1 or b>self.ns: raise`, but sites are indexed 0..ns-1 and the valid bond range is 1..ns-1 — `compute_entropy`'s own loop three lines above uses `range(1,self.ns)`, i.e. the module already knows the right bound. `bond_entropy(self,wf,i,j)` forms `b = max(i,j)`, so `get_bond_entropy(wf, ns-1, ns)` — an out-of-range site index, the kind of GIGO the 2026-08 audit's findings #8 and #21 cover at other sites — slips past the guard and reaches ITensor, which calls `abort()`. The process dies with a core dump; `try/except Exception` cannot catch it, so a parameter loop or a test session is terminated outright rather than reporting a bad index. `itensor_version="python"` raises a plain `IndexError` for the same call, which is the behaviour the C++ backends should have too.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 3 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/dispatch-matrix/bond.py 3    # n=4 Heisenberg chain, loops get_bond_entropy(wf,b-1,b) for b=1..4
```

Observed:

```
### itensor_version=2
  b=2 -> 0.319368
  b=3 -> 0.693147
From line 101, file itensor_operators.cc
Default constructed ITensor in product
Default constructed ITensor in product
timeout: the monitored command dumped core

### itensor_version=3
  b=2 -> 0.319368
  b=3 -> 0.693147
From line 940, file itensor.cc
Default constructed ITensor in product
Default constructed ITensor in product
timeout: the monitored command dumped core

### itensor_version=python
  b=1 -> 0.693147
  b=2 -> 0.319368
  b=3 -> 0.693147
  b=4 -> EXC IndexError: list index out of range
```

**Expected**: `get_bond_entropy(wf, 3, 4)` on a 4-site chain names site 4, which does not exist — it must raise a catchable Python error naming the out-of-range index on every backend (as "python" does), not abort the interpreter. The guard that is already present is clearly meant to do this and simply has the wrong comparison.

**Suggested fix**: Change `if b<1 or b>self.ns: raise` in `entropy.compute_entropy_single` to `if b<1 or b>=self.ns: raise IndexError(...)` naming `b` and the valid range 1..ns-1, matching `compute_entropy`'s own `range(1,self.ns)`; and have `bond_entropy` validate both `i` and `j` against `0..ns-1` before forming `max(i,j)`, so the message points at the site the caller actually passed.

**Reviewer (CONFIRMED)**: Reproduced on all three backends by re-running bond.py under timeout, one backend per process. itensor_version=2: b=1..3 print 0.693147/0.319368/0.693147, then b=4 prints ITensor's "From line 101, file itensor_operators.cc / Default constructed ITensor in product" and the shell reports the command dumped core. itensor_version=3: same, from itensor.cc:940. itensor_version="python": b=4 raises a catchable IndexError. So the cross-backend divergence and the uncatchable abort are both real. The guard is dmrgpy's own code, not vendored: entropy.py:31 reads `if b<1 or b>self.ns: raise`, three lines below compute_entropy's own `for i in range(1,self.ns)` — the module already encodes the correct upper bound and the guard contradicts it; bond_entropy (entropy.py:38-42) forms b=max(i,j) with no validation of i or j. Scope checked: audit #16's reviewer wrote that this abort path is "unreachable from the public API (mps.get_entropy loops b in range(0,ns-1) and get_bond_entropy(i,j=i+1) caps compute_entropy_single at b=ns-1)" — that statement silently assumes a valid i, and this finding corrects the record by passing i=ns-1,j=ns, which is exactly the kind of out-of-range site index audit #8 and #21 treat as in-scope. #21 (unvalidated maxm<=0 -> uncatchable SIGABRT on v2/v3) is the same class and was FIXED by validating on assignment, so this is a new instance of an already-accepted defect class rather than something the project has decided to tolerate.


---


### 28. self.mode is returned unvalidated by resolve_mode, so a mistyped mode string makes get_gs() return None and gs_energy() raise 'No active exception to reraise'

`hole` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `dispatch-matrix`

**Where**: `src/dmrgpy/mode.py:105-110 (resolve_mode returns self.mode unchecked); src/dmrgpy/manybodychain.py:1064-1074 (get_gs has no else); src/dmrgpy/manybodychain.py:1078-1087 (gs_energy's bare raise)`

`resolve_mode` does `if self.mode is not None: return self.mode` before it ever checks the string against the two legal values — the validation `if mode in ["ED","DMRG"] ... else: print(...); raise` applies only to the *call argument*, never to the chain attribute. `Many_Body_Chain.mode` is the documented way to pin a solver (`sc.mode = "ED"`), so a lowercase or misspelled assignment is a realistic user error. `get_gs()` then falls off the end of its `if/elif` and returns `None` — the caller gets a null wavefunction with no error at all — while `gs_energy()` hits a bare `raise` with no active exception and reports `RuntimeError: No active exception to reraise`, which names neither the attribute nor the typo. This is the same §4.10 shape (a dispatch decision taken before the information that qualifies it) as the `else`-serving-two-jobs cases the audit fixed elsewhere.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 3 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/dispatch-matrix/lower.py    # sets sc.mode = "ed" (lowercase) on a 4-site chain
```

Observed:

```
get_gs() -> None
gs_energy() -> RuntimeError No active exception to reraise
```

**Expected**: An unrecognized `self.mode` should raise a `ValueError` naming the offending value and the accepted set ("ED", "DMRG", or None), at the point the attribute is read — the same treatment `timedependent.check_tevol_method()` now gives `tevol_method` after audit finding #2. `get_gs()` must never return `None` silently: a returned `None` propagates into `overlap`/`vev`/entropy calls and fails much further away from the cause.

**Suggested fix**: In `mode.py::resolve_mode`, validate the attribute where it is read: `if self.mode is not None: if self.mode not in ("ED","DMRG"): raise ValueError(...); return self.mode`. Independently, give `Many_Body_Chain.get_gs` a trailing `else: raise ValueError(...)` so it can never fall through to an implicit `None`, and replace `gs_energy`'s bare `raise` with the same explicit error.

**Reviewer (CONFIRMED)**: Reproduced verbatim: running lower.py (4-site chain, sc.mode="ed" lowercase) gives `get_gs() -> None` and `gs_energy() -> RuntimeError No active exception to reraise`. Source matches exactly: mode.py:105-107 does `if self.mode is not None: return self.mode` BEFORE the `if mode in ["ED","DMRG"]` check, which therefore only ever validates the call argument, never the attribute; manybodychain.get_gs (:1064-1074) has an `if mode=="DMRG" / elif mode=="ED"` with no else, so it falls off the end and implicitly returns None; gs_energy (:1078-1087) ends in a bare `raise` with no active exception. The silent None is the load-bearing half — it propagates into overlap/vev/entropy and fails far from the cause. sc.mode is a documented user-facing knob (CLAUDE.md: 'Most public methods accept a mode="DMRG"|"ED" kwarg'; the audit's own repro scripts set sc.mode="ED" directly), so a typo is a realistic input rather than a contrived one. In scope: the codebase already treats this exact shape as a defect — audit #2 fixed an unvalidated tevol_method string (timedependent.check_tevol_method) and audit #14's fix replaced get_distribution_moments' bare raise with a named NotImplementedError. No entry for resolve_mode/self.mode validation exists in docs/audit_2026_08_hole_hunt.md (grepped).


---


### 29. kpmdmrg.general_kpm_moments still has a bare `raise` for a missing X, so get_distribution() reports "RuntimeError: No active exception to reraise"

`hole` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `cpp-v3-completeness`

**Where**: `src/dmrgpy/kpmdmrg.py:197`

`general_kpm_moments` opens with `if X is None: raise` — a bare re-raise with no active exception, so the caller sees `RuntimeError: No active exception to reraise` instead of a message naming the required argument. `Many_Body_Chain.get_distribution()` forwards straight into it (manybodychain.py:967 -> distribution.py:11 -> kpmdmrg.py:252 -> :197), and `get_distribution` itself documents no required argument in its docstring ("Return the distribution of an operator's spectrum"), so a caller who omits `X=` gets no usable signal. This is the same defect class the 2026-08 audit fixed at kpmdmrg.py:104 (finding #12, the `name` type check) and at two adjacent sites; this sibling in the same file was not swept. A second instance of the same class sits at entropytk/correlationentropy.py:84 (`print("Unrecognized type",type(self)); raise`), reached by calling `get_correlation_matrix()` on a non-fermionic chain — it produced the identical RuntimeError on a Spin_Chain in my sweep.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/cpp-v3-completeness && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 -u final_repro.py    # section C: 4-site Heisenberg Spin_Chain(itensor_version=3), sc.get_distribution(es=np.linspace(-1,1,5))
```

Observed:

```
C get_distribution() with no X   RuntimeError: No active exception to reraise

(full traceback from an earlier run of the same call, scratchpad/s4_dist.py:)
  File "/home/joselado/Documents/programs/dmrgpy/src/dmrgpy/kpmdmrg.py", line 197, in general_kpm_moments
    if X is None: raise
                  ^^^^^
RuntimeError: No active exception to reraise
```

**Expected**: A `TypeError`/`ValueError` naming the missing argument, e.g. "general_kpm_moments: X= (the operator whose spectral distribution is wanted) is required" — the wording the audit's own fix to kpmdmrg.py:104 adopted for the neighbouring check in this file.

**Suggested fix**: Replace `if X is None: raise` with an explicit `raise TypeError("get_distribution/general_kpm_moments: X= is required (the operator whose spectral distribution to compute)")`, and do the same for correlationentropy.py:84 (`raise TypeError("get_correlation_matrix: no default operators for chain type "+type(self).__name__+"; pass operators=")`).

**Reviewer (CONFIRMED)**: Reproduced verbatim in section C: sc.get_distribution(es=...) with no X= gives RuntimeError: No active exception to reraise, from kpmdmrg.py:197's `if X is None: raise`. I also reproduced the second instance the finding mentions in passing -- get_correlation_matrix() on a 4-site Spin_Chain prints 'Unrecognized type <class dmrgpy.spinchain.Spin_Chain>' and then raises the same RuntimeError from correlationentropy.py:84.

Checked the two refutation angles that could apply. (a) Already recorded? Audit finding #12 (docs/audit_2026_08_hole_hunt.md:1554-1557) names kpmdmrg.py:104 specifically -- the `name` type check -- and is marked fixed; line 197 is a different check in the same function family and was not swept. The user_guide's 'Was/Now' table likewise records only the name="ZZ" case as fixed. So this instance is live and unrecorded, not a re-report. (b) Is X actually required, i.e. is the bare raise merely guarding an impossible call? No -- I confirmed get_distribution(X=sc.Sz[0], es=...) returns a spectrum fine, and there is no default operator anywhere on the path (distribution.get_distribution -> general_kpm -> general_kpm_moments), while get_distribution's docstring documents no required argument. So a caller who omits X gets a message that names neither the function nor the argument.

Worth flagging for triage: this is off the stated lens. A bare `raise` in backend-agnostic Python has nothing to do with mpscpp3 completeness against mpscpp2/pyitensor. Lens mismatch is not one of the brief's refutation criteria, so it survives as a real (if minor) defect, but it should not be counted as evidence about the v3 backend.


---


### 30. Bosonic_Chain and Parafermionic_Chain constructors still reject itensor_version=, the kwarg deb6bf8 added to the fermionic subclasses and the user guide tells boson users to pass

`hole` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `recent-commits`

**Where**: `src/dmrgpy/bosonchain.py:10 (`def __init__(self,n,maxnb=None)`), :80 (`SpinBoson_Chain.__init__(self,sitesin,n=None,maxnb=None)`), src/dmrgpy/parafermionchain.py:7 (`def __init__(self,n,Z=3)`), src/dmrgpy/spinfermionchain.py:6; contrast src/dmrgpy/fermionchain.py:13/207/231/466/647 (all `**kwargs`)`

In-window commit deb6bf8 ("Forward **kwargs from the fermionic subclass constructors") gave every `fermionchain` class a `**kwargs` passthrough so `itensor_version=` reaches `Many_Body_Chain.__init__`. `Bosonic_Chain`, `SpinBoson_Chain`, `Parafermionic_Chain` and `Spin_Fermion_Hamiltonian` were left with fixed signatures, so the kwarg is a TypeError there. This collides with the boson paragraph the in-window doc-sync commit wrote (user_guide.md:119-129), which tells users a non-default `maxnb` "should be run under itensor_version=3 ... or \"python\"" -- the natural spelling of which is the constructor kwarg. The workaround (construct, then `setup_python()`/`setup_cpp(version=3)`) exists but is undocumented for these classes, and on the v2 route it is also the thing that trips finding #2's abort. It is not a wrong-answer bug: with no kwarg the classes go through `Many_Body_Chain.__init__(sites)` and correctly pick up `cppext.default_backend()`.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/recent-commits && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 -c "
from dmrgpy import bosonchain, parafermionchain
for f,lab in ((lambda: bosonchain.Bosonic_Chain(3,maxnb=[6,6,6],itensor_version='python'),'Bosonic_Chain'),
              (lambda: parafermionchain.Parafermionic_Chain(4,itensor_version='python'),'Parafermionic_Chain')):
    try: f(); print(lab,'accepted')
    except TypeError as e: print(lab,'->',type(e).__name__,e)"
```

Observed:

```
Bosonic_Chain -> TypeError Bosonic_Chain.__init__() got an unexpected keyword argument 'itensor_version'
Parafermionic_Chain -> TypeError Parafermionic_Chain.__init__() got an unexpected keyword argument 'itensor_version'

# for contrast, the same sweep over every chain class with no kwarg, under a
# simulated missing C++ extension (cppext._backends[2]=_backends[3]=None):
Spin_Chain                         itensor_version='python'
Fermionic_Chain                    itensor_version='python'
Spinful_Fermionic_Chain            itensor_version='python'
Spinful_Fermionic_Chain_Native     itensor_version='python'
Bosonic_Chain                      itensor_version='python'
Parafermionic_Chain                itensor_version='python'
Mixed_Spin_Fermion_Chain           itensor_version='python'
```

**Expected**: `Bosonic_Chain(3, maxnb=[6,6,6], itensor_version="python")` should construct a chain on the pure-Python backend, the way `Fermionic_Chain(4, itensor_version="python")` does since deb6bf8 and the way user_guide.md:127 implies.

**Suggested fix**: Add `**kwargs` to `Bosonic_Chain.__init__`, `SpinBoson_Chain.__init__`, `Parafermionic_Chain.__init__` and `Spin_Fermion_Hamiltonian.__init__` and forward it to the `Many_Body_Chain.__init__` call already in each body -- the same one-line change deb6bf8 made in fermionchain.py. Note `Parafermionic_Chain` calls `self.get_operator(...)` before `Many_Body_Chain.__init__`, so forward the kwargs at the existing init call site rather than reordering.

**Reviewer (CONFIRMED)**: Reproduced exactly: Bosonic_Chain(3, maxnb=[6,6,6], itensor_version='python') and Parafermionic_Chain(4, itensor_version='python') both raise TypeError ('unexpected keyword argument'), while Fermionic_Chain(4, itensor_version='python') and Spin_Chain(['S=1/2']*4, itensor_version='python') construct and report itensor_version='python'. Signatures verified in source: bosonchain.py `def __init__(self,n,maxnb=None)`, SpinBoson_Chain `(self,sitesin,n=None,maxnb=None)`, parafermionchain.py `def __init__(self,n,Z=3)`, spinfermionchain.py `Spin_Fermion_Hamiltonian.__init__(self,sites)` -> super().__init__(len(sites)) with nothing forwarded; fermionchain.py uses **kwargs throughout. In-window: deb6bf8 is 2026-08-26 and did exactly this for the three fermionic subclasses. The hunter's own caveat is accurate -- with no kwarg these classes still pick up cppext.default_backend() correctly, so there is no wrong answer here, only an inconsistent surface. The note about Parafermionic_Chain calling get_operator before Many_Body_Chain.__init__ also checks out in the source.


---


### 31. Infinite_Many_Body_Chain.get_operator's docstring says gs_method="idmrg" raises for couplings reaching past one unit cell; it does not, and the answer is exact

`hole` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `recent-commits`

**Where**: `src/dmrgpy/infinitechain.py:485-487 (get_operator docstring), contradicted by :504-534 (_require_reach_one's own docstring) and by its only call site, :919`

`get_operator(name, i, group=<int>)` -- the API 71ba8eb added so a user can write a coupling longer than one unit cell -- documents that integer form as "supported by set_hamiltonian and by gs_method=\"vumps\" ...; gs_method=\"idmrg\" and excitation_energies/excitation_gap are reach-1 only and raise for it". That is wrong for iDMRG. `_require_reach_one` is called from exactly one place (line 919, `excitation_energies`), and its own docstring says the opposite in as many words ("gs_method=\"idmrg\"'s growth loop consumes whatever automaton _build_periodic_mpo hands it, which has always carried one pending channel per site of a term's reach"), as does CLAUDE.md ("gs_method=\"idmrg\" needed no change at all ... so do NOT 'fix' that by adding a guard"). Executed on a reach-2 chain, iDMRG returns the exact answer on both backends. Consequence: a user reading the docstring at the point of use rewrites their model on an n_uc>=R cell -- the exact cost 71ba8eb existed to remove -- or avoids iDMRG entirely.

Repro:

```bash
cd /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/recent-commits && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 reach2.py   # H = -Sz_i + 0.5 Sz_i Sz_{i+1} + 0.25 Sz_i Sz_{i+2} on a 1-site cell
```

Observed:

```
backend=python  idmrg  e=-0.31249999999999983  <Sz>=(0.5+0j)
backend=3       idmrg  e=-0.31249999999999867  <Sz>=(0.5+0j)
backend=3       vumps  e=-0.3125000000000001  <Sz>=(0.5+0j)
# exact for the polarized state: -1*(0.5) + 0.5*0.25 + 0.25*0.25 = -0.3125
# (backend=python vumps is finding #1 and raised here)
```

**Expected**: The docstring should say what the code does: both ground-state methods handle any finite reach, and only the tangent-space excitation ansatz (excitation_energies/excitation_gap) is reach-1 -- which is precisely what `_require_reach_one`'s docstring, thirty lines below, already says correctly.

**Suggested fix**: In `get_operator`'s docstring replace "`gs_method=\"idmrg\"` and `excitation_energies`/`excitation_gap` are reach-1 only and raise for it" with "both `gs_method` values handle any finite reach; only `excitation_energies`/`excitation_gap` are reach-1 and raise for it (see `_require_reach_one`)". Check docs/user_guide.md and docs/documentation.md for the same claim while there.

**Reviewer (CONFIRMED)**: Executed my own version (v4.py) on the reach-2 polarized chain, 1-site cell, maxm=8: gs_method='idmrg' returns e=-0.31249999999999994 on backend 'python' and -0.3124999999999999 on backend 3, against the exact -0.3125, with <Sz>=0.5 on both -- no raise on either. Only excitation_energies raises (NotImplementedError from _require_reach_one), which is what the docstring should have said.

Source checks: _require_reach_one has exactly one caller, infinitechain.py:919 (excitation_energies/excitation_gap); its own docstring states 'Both GROUND-STATE algorithms handle any finite reach and are not gated by this'; the module-level comment at infinitechain.py:191-195 says the same correctly; docs/user_guide.md:2803 says the same correctly; CLAUDE.md says 'gs_method="idmrg" needed no change at all ... do NOT fix that by adding a guard'. So get_operator's docstring is the lone dissenter. `git log -L 484,487:src/dmrgpy/infinitechain.py` pins the wrong sentence to in-window commit 71ba8eb, i.e. it is stale text introduced by the very commit that removed the restriction. Docstring-only, no behavioural impact -- low is right.


---


### 32. idmrg_window re-solves the call-invariant transfer-matrix fixed points on every _close_array_chain / local_expectation call — 25 full rebuilds for a 4-step td_dynamical_correlator

`optimization` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `pyitensor-performance`

**Where**: `src/dmrgpy/pyitensor/idmrg_window.py:986-1000 (inside _close_array_chain) and :949-951 (inside local_expectation)`

`_close_array_chain` calls `_idmrg_mod._transfer_matrices(cell, n_cell)`, then `_all_right_fixed_points(Es, n_cell)`, then `_all_left_fixed_points(Es, n_cell)` — on every invocation. None of those three depends on the call's arguments: they are functions of `result` alone (the converged unit cell). `local_expectation` builds another copy on its own before calling `_close_array_chain` twice. Since `snapshot_correlator` calls `_close_array_chain` once per x value and `dynamical_correlator_td` calls `snapshot_correlator` once per time step, a `td_dynamical_correlator(nt=200, x_values=[...5 values...])` run re-solves ~2000 pairs of transfer-matrix eigenproblems that could be solved once. The `E_id` calibration chain in the same function is likewise a pure function of (result, p_left, chain length) and is recomputed per call. Both idmrg.py (`_CorrelatorEnv` / `_correlator_env`, whose docstring says rebuilding per call 'made a seven-point correlator sweep cost seven-plus full eigensolves instead of one') and mpscpp3 (`iw_build_cache`) already have the cache this module lacks. At maxm=16 with nt=4 this is 0.47 s of a 19.9 s call — small only because finding 1's chi^6 closure dwarfs it; once that is fixed it becomes the dominant remaining cost, and it carries finding 3's chi^4 scaling.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/pyitensor-performance/s9_window.py 16 4 6
```

Observed:

```
## td_dynamical_correlator maxm=16 nt=4 n_window=6 wall=19.94s
   _close_array_chain         ncalls=24     tot=0.007 cum=19.270
   _transfer_matrices         ncalls=25     tot=0.000 cum=0.026
   _all_right_fixed_points    ncalls=25     tot=0.000 cum=0.224
   _all_left_fixed_points     ncalls=25     tot=0.000 cum=0.217
   _dominant_fixed_point      ncalls=74     tot=0.005 cum=0.830
   eigs                       ncalls=74     tot=0.010 cum=0.813
   snapshot_correlator        ncalls=12     tot=0.008 cum=10.192   (nt=6 variant, same shape)

(with the same maxm=12/nt=6 run: _transfer_matrices ncalls=37, _all_right_fixed_points ncalls=37,
 _all_left_fixed_points ncalls=37 — one per _close_array_chain call, 36 of them)
```

**Expected**: These quantities are properties of the converged `IDMRGResult`, not of the snapshot being measured, so they should be built once and memoized on the result — exactly as `idmrg._correlator_env` does for the identical objects on the finite-observable path, and as `mpscpp3`'s `iw_build_cache` does for the v3 window. The `E_id` denominator should be memoized on (result, p_left, nsites) alongside them.

**Suggested fix**: Add an `_IWEnv`-style cache on the IDMRGResult holding `(cell, n_cell, rho_after, l_before)`, invalidated by `cell is not cached_cell` the way `_correlator_env` already does it, and have both `_close_array_chain` and `local_expectation` read it. Memoize the `close(E_id)` calibration value per (p_left, len(ket_arrays)) in the same object. (The runtime patch in s16_window_patch.py does exactly this together with finding 1's fix and reproduces S(k,omega) to 5.1e-13.)

**Reviewer (CONFIRMED)**: Counted the calls directly rather than reading them off a profile: I wrapped idmrg._transfer_matrices with a counter in my own harness (v_window_split.py) and the shipped run of td_dynamical_correlator(maxm=16, nt=4, n_window=6, 5 x-values) makes exactly 25 calls — 24 _close_array_chain invocations plus the one local_expectation does on its own — while the memoized variant makes 2. Confirmed the quantities really are call-invariant: caching them changes S(k,omega) by 6.1e-13 relative across the whole nt=4 run, and _close_array_chain's three rebuilt objects (_transfer_matrices / _all_right_fixed_points / _all_left_fixed_points on _window_cell(result)) take no argument derived from the call. Isolated the cost, which is what makes the 'low' severity honest but also shows it is not zero: with finding 1's contraction fix applied, propagate-without-cache is 1.254 s and propagate-with-cache 0.713 s — the caching is 43% of the remaining time, exactly the hunter's 'small only because finding 1 dwarfs it; once that is fixed it becomes the dominant remaining cost'. Prior art for the fix is in the repo: idmrg._correlator_env/_CorrelatorEnv memoizes the identical objects with an `env.cell is not cell` invalidation, and its docstring records the same defect being fixed on the finite-observable path. Not covered by any docs/audit_* or known_issue_* item.


---


### 33. examples/readme_examples/energy_VS_length/main.py is an empty stub: it computes no energy, sweeps no length, prints nothing and plots nothing

`hole` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `examples/readme_examples/energy_VS_length/main.py (whole file, 21 lines)`

The whole script builds one fixed 30-site S=1/2 Heisenberg chain, calls `set_hamiltonian`, then does `import matplotlib.pyplot as plt` and ends. There is no `gs_energy()` call, no loop over length, no `print`, and no plotting call — it exits 0 having produced no output at all. CLAUDE.md describes `examples/readme_examples/` as mirroring the snippets shown in README.md, but README.md has no 'energy versus length' section (its section list runs from 'Ground state energy of an S=1/2 spin chain' through 'Ground state of an infinite chain with iDMRG' with no such entry), so this directory mirrors nothing. It also violates the 'examples should plot' rule twice over: it imports pyplot without using it, and the quantity its own name promises is a sweep.

Repro:

```bash
cd /home/joselado/Documents/programs/dmrgpy/examples/readme_examples/energy_VS_length && cat main.py && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src MPLBACKEND=Agg taskset -c 3 python3 main.py ; echo "EXIT=$?" ; grep -n 'energy_VS_length\|energy vs' /home/joselado/Documents/programs/dmrgpy/README.md
```

Observed:

```
# (full file, after the sys.path preamble)
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(30)] # spins in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # create the Hamiltonian


import matplotlib.pyplot as plt

# (end of file -- 21 lines total)

# run: no stdout at all
EXIT=0
# grep over README.md for a matching section: no matches
```

**Expected**: Either a working example that sweeps chain length and plots E0 (and/or E0/n) versus n — which is what the directory name advertises and what the '105 scripts that had quietly settled for print-only' audit was about — or deletion of the directory, since README.md has no snippet for it to mirror.

**Suggested fix**: Add the missing loop: `for n in [4,8,12,16,20,24,30]` building the chain, collecting `sc.gs_energy()`, then `plt.plot(ns, es/np.array(ns))` for the energy density with the Bethe-ansatz value 0.25-ln2 as a horizontal reference. That is a few seconds of compute at these sizes and turns a dead directory into the regression the name promises.

**Reviewer (CONFIRMED)**: Read the whole file and ran it pinned with MPLBACKEND=Agg: zero stdout, EXIT=0. The file is exactly the sys.path preamble, a fixed 30-site Heisenberg chain, set_hamiltonian(h), a bare `import matplotlib.pyplot as plt`, and end of file -- no gs_energy() call, no length loop, no print, no plotting call. grep -ci over README.md for 'energy_VS_length|energy vs length|versus length' returns 0, so the directory mirrors no README snippet, which is what examples/readme_examples/ exists for per CLAUDE.md. Nothing to explain away: this is not a convergence or tolerance question and is recorded in no audit or ROADMAP entry.


---


### 34. Two examples compute a full sequence and only print it, against CLAUDE.md's 'examples should plot' rule — both already import pyplot and never call it

`hole` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `examples/groundstate/energy_fluctuation/main.py:21-26 (7-point maxm sweep); examples/utilities/multioperator_density/main.py:32-39 (6-site density profile, DMRG and ED)`

A grep over all 286 example scripts for 'imports matplotlib but never calls a drawing method' returned 13 hits. Ten of them are the legitimate single-scalar exception CLAUDE.md allows (a pair of ED-vs-DMRG ground-state energies and nothing else). The other three are real violations: `readme_examples/energy_VS_length` (reported separately, it is broken outright), plus these two. `groundstate/energy_fluctuation` loops maxm over [1,2,5,10,20,30,40] and prints an energy and a fluctuation for each — a textbook convergence-versus-bond-dimension curve, and the obvious axis is already the loop variable. `utilities/multioperator_density` builds a 6-site density profile twice, once under mode='DMRG' and once under mode='ED', and prints both lists — a site axis plus a backend overlay, both already computed. Both scripts carry `import matplotlib.pyplot as plt` at the top and never touch it, which is the signature CLAUDE.md describes ('usually because the script already computed an array or a comparison and simply never plotted it').

Repro:

```bash
cd /home/joselado/Documents/programs/dmrgpy/examples && for f in $(find . -name '*.py'); do if grep -qE 'import matplotlib|import pylab' $f && ! grep -qE '\.(plot|scatter|imshow|errorbar|semilogy|semilogx|loglog|bar|barh|fill_between|pcolormesh|hist|step|contourf|axhline|axvline|stem|plot_surface)\(' $f; then echo "NOPLOTCALL $f"; fi; done
```

Observed:

```
NOPLOTCALL ./fermion_models/spinless_fermions/main.py
NOPLOTCALL ./fermion_models/fermionic_energy/main.py
NOPLOTCALL ./fermion_models/majorana_chain/main.py
NOPLOTCALL ./fermion_models/charge_gap/main.py
NOPLOTCALL ./finite_temperature/finite_temperature_ground_state/main.py
NOPLOTCALL ./magnetization/total_spin/main.py
NOPLOTCALL ./readme_examples/energy_VS_length/main.py
NOPLOTCALL ./magnetization/constrain_density/main.py
NOPLOTCALL ./groundstate/energy_fluctuation/main.py
NOPLOTCALL ./utilities/multioperator_density/main.py
NOPLOTCALL ./utilities/multioperator_vev/main.py
NOPLOTCALL ./spin_models/parafermion_energy/main.py
NOPLOTCALL ./spin_models/bilinear_biquadratic/main.py

# the two with a real sequence:
# groundstate/energy_fluctuation/main.py
for maxm in [1,2,5,10,20,30,40]: # loop over bond dimension
    ...
    print("Energy",e,"fluctuation",de,"for bond dimension",maxm)
# utilities/multioperator_density/main.py
print("Density DMRG",[fc.vev(di,mode="DMRG").real for di in den])
print("Density ED",[fc.vev(di,mode="ED").real for di in den])
```

**Expected**: Per CLAUDE.md, every examples/*/*/main.py that produces a sequence of values should end with a matplotlib plot of it; the single-scalar exception explicitly does not apply when 'an obvious cheap axis was sitting right there'. Here the axes are the loop variable (maxm) and the site index, both already in hand.

**Suggested fix**: energy_fluctuation: collect the sweep into lists and `plt.semilogy(maxms, des)` with E(maxm) on a twin axis — the fluctuation falling with bond dimension is the whole point of the script. multioperator_density: `plt.plot(range(n), den_dmrg, marker='o', label='DMRG')` and the same for ED, which also makes the backend agreement visible instead of leaving it to the reader to diff two printed lists. Note that energy_fluctuation's printed numbers are themselves affected by the gs_energy_fluctuation findings above if the chain ever routes to ED, so re-check it after those are fixed.

**Reviewer (CONFIRMED)**: Read both files in full rather than trusting the grep. groundstate/energy_fluctuation/main.py imports pyplot at line 5, loops maxm over [1,2,5,10,20,30,40] computing e and de, prints them, and ends -- the loop variable is a ready-made x-axis and nothing is drawn. utilities/multioperator_density/main.py imports pyplot, builds a 6-site density profile and prints it twice, once mode='DMRG' and once mode='ED' -- a site axis plus a backend overlay, already computed, never plotted. Both match CLAUDE.md's own description of the 105-script class ('already computed an array or a comparison and simply never plotted it'), and neither qualifies for the single-scalar exception, which CLAUDE.md explicitly says does not apply when an obvious cheap axis is sitting there. I did not re-audit all 286 scripts, so I am vouching only for these two, which are the substance of the finding.


---


### 35. get_four_correlation_tensor's public docstring names the wrong default-resolver order: it claims 'sweep' first and never mentions 'batched' or 'fold'

`hole` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `src/dmrgpy/entropytk/correlationentropy.py:237-246 (the docstring); the actual resolver is _four_correlation_tensor_default_ctmode at :436-495`

`get_four_correlation_tensor(wf, ctmode=None)` is the public dispatcher and its docstring describes the auto-selection as '"sweep" whenever it applies (itensor_version in (3,"python"), non-native-spinful fermionic sites), else "full" ... else "explicit"' — a three-way choice. Two more modes were added since: `"batched"` (pyitensor-only, and now the *first* thing the resolver tries, measured 15-28x faster than sweep) and `"fold"` (for native spinful sites under v3). The resolver's own docstring at :436 is correct and lists all five, and docs/user_guide.md:1135-1141 is correct too ('"batched" whenever it applies ... else "sweep" ... else "full" ... else the always-correct "explicit" fallback'), so the public docstring is the one artefact left describing the pre-batched behaviour. A caller reading `help(get_four_correlation_tensor)` to decide whether to pass an explicit ctmode is told the default is 'sweep' when it is in fact 'batched' on the pure-Python backend.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONPATH=/home/joselado/Documents/programs/dmrgpy/src taskset -c 5 python3 /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/ctdoc.py   # prints the docstring, then asks the resolver what it picks on a 4-site Fermionic_Chain under 'python' and under 3
```

Observed:

```
PUBLIC DOCSTRING of get_four_correlation_tensor:
Return the correlation tensor as <Cdag_i C_j Cdag_k C_l>.

    ctmode=None (the default) auto-selects the fastest method actually
    available for this wavefunction's backend/chain type -- "sweep"
    whenever it applies (itensor_version in (3,"python"), non-native-
    spinful fermionic sites), else "full" whenever it applies, else the
    always-correct but slowest "explicit" fallback (see
    _four_correlation_tensor_default_ctmode()). Passing a ctmode
    explicitly is still a hard request: it raises rather than silently
    falling back if that method isn't available for this wavefunction.
itensor_version='python' -> resolver picks 'batched'
itensor_version=3 -> resolver picks 'sweep'
```

**Expected**: The docstring should describe the order the resolver implements — batched, then sweep, then fold (native spinful), then full, then explicit — matching both `_four_correlation_tensor_default_ctmode`'s own docstring and user_guide.md §5, which already say exactly that.

**Suggested fix**: Replace the three-mode sentence with the five-mode order, or shorten it to 'auto-selects the fastest method available for this wavefunction — see _four_correlation_tensor_default_ctmode() for the order' so there is one place to keep in sync rather than two.

**Reviewer (CONFIRMED)**: Read both docstrings and the resolver body rather than taking the hunter's word. The public docstring at correlationentropy.py:238-246 describes a three-way choice ('sweep' whenever it applies, else 'full', else 'explicit') and never names 'batched' or 'fold'. The resolver _four_correlation_tensor_default_ctmode at :436 documents five modes with batched first, and its code confirms it: the first branch returns 'batched' for itensor_version=='python' with four_correlation_tensor_batched on the session, then 'sweep' for version 3/'python', then 'fold' for native spinful, then 'full', then 'explicit'. The dispatcher at :250-263 does have live elif branches for both 'fold' and 'batched'. user_guide.md:1130-1142 is already correct and names batched/fold, so the public docstring is the sole stale artefact. Not a convergence question and not recorded anywhere.


---


### 36. documentation.md's VUMPS section reports its measurements 'at the default nrestarts=6'; the default is 4 and has never been 6

`hole` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `docs-examples-drift`

**Where**: `docs/documentation.md:505 (and the matching passage in docs/documentation.tex); code at src/dmrgpy/infinitechain.py:369 (self.vumps_nrestarts = 4) and src/dmrgpy/pyitensor/vumps.py:924 (nrestarts=4)`

The 'Let VUMPS converge' section reports its before/after timing table with the phrase 'Measured through the public driver at the default `nrestarts=6`'. `Infinite_Chain.__init__` sets `self.vumps_nrestarts = 4`, and `pyitensor.vumps.vumps_ground_state`'s own signature default is also 4. `git log -L369,369:src/dmrgpy/infinitechain.py` shows the attribute was introduced at 4 in b5593ec and never changed, so this is not a default that moved after the measurement was written — the sentence was wrong when written. The number matters because the table is the evidence for the residual-criterion fix (0/3 -> 3/3 converged), and a reader reproducing it at the real default 4 is running a different experiment from the one tabulated. Note the user guide's own VUMPS snippet is fine: it writes `ic.vumps_nrestarts = 6` as an explicit setting, not as the default.

Repro:

```bash
cd /home/joselado/Documents/programs/dmrgpy && sed -n '505,506p' docs/documentation.md && grep -n 'vumps_nrestarts' src/dmrgpy/infinitechain.py && grep -n 'nrestarts=4' src/dmrgpy/pyitensor/vumps.py && git log -L369,369:src/dmrgpy/infinitechain.py --oneline | head -12
```

Observed:

```
# docs/documentation.md:505
Measured through the public driver at the default `nrestarts=6`, before

# code
src/dmrgpy/infinitechain.py:369:        self.vumps_nrestarts = 4  # gs_method="vumps" only: independent
src/dmrgpy/pyitensor/vumps.py:924:                        nrestarts=4, verbose=False):

# git log -L on that line
b5593ec Add VUMPS ground-state solver; fix iDMRG excitation diagram 6a's |n|>=2 tail
+        self.vumps_nrestarts = 4  # gs_method="vumps" only: independent
(no later change to this line)
```

**Expected**: Either 'at nrestarts=6' (dropping the word 'default', since 6 is what the measurement used), or the default raised to 6 in infinitechain.py and vumps.py if 6 is genuinely the recommended setting the text implies.

**Suggested fix**: Change 'at the default `nrestarts=6`' to 'at `nrestarts=6` (the default is 4)' in both documentation.md and documentation.tex. Line 526's '`D=8`, `nrestarts=6`' in the same section is already phrased correctly and needs no change.

**Reviewer (CONFIRMED)**: Verified the three facts independently. docs/documentation.md:505 reads exactly 'Measured through the public driver at the default `nrestarts=6`, before', and docs/documentation.tex:773 carries the same sentence, so both formats are wrong. The code default is 4 in all three places I could find: infinitechain.py:369 (self.vumps_nrestarts = 4), pyitensor/vumps.py:924 and pyitensor/vumps_ms.py:772 (both nrestarts=4 in the signature). I grepped for any backend-specific override of 6 and found none -- the only other '6' is a prose mention in vumps.py:126. git log -L369,369:src/dmrgpy/infinitechain.py shows the line introduced as 4 in b5593ec with no later change, so this is not a default that moved after the text was written. Low severity but factually indisputable.


---


## Refuted


- **mpscpp3 bindings.cc exposes Chain::idmrg_local_excitation_gap with no caller anywhere in the repo (Python only ever uses the _detail variant)** (lens `cpp-v3-completeness`) -- Every fact the hunter reports is accurate -- I checked all of them -- but none of them adds up to a defect, and certainly not to a completeness hole.

Facts verified: `grep -rn '\.idmrg_local_excitation_gap\b' --include=*.py` over src/tests/examples/benchmarks returns 0; infinitechain.py:1153 calls only the _detail form; and I enumerated all 78 .def names in mpscpp3/bindings.cc and cross-grepped each -- the only three with no Python caller are idmrg_local_excitation_gap, vms_onsite_expectation and vms_two_point_correlator, and the latter two are indeed reached from C++ (chain_session.h:3957-3958 and 4042, vumps_onsite_expectation/vumps_two_point_correlator falling through on !have_vumps_snapshot_ && have_vms_snapshot_), so under the hunter's own C++-inclusive reading their 'only one' claim holds. I am not refuting on a factual quibble.

What does not follow is the conclusion. (1) The lens is mpscpp3 *incompleteness* relative to mpscpp2/pyitensor; an extra exposed method is the opposite of a missing one. (2) Nothing crashes, no number is wrong, no supported call path is affected -- the hunter concedes 'Harmless'. (3) The C++ method is a one-line convenience overload, `idmrg_local_excitation_gap(int niter) const { return idmrg_local_excitation_gap_detail(niter).gap; }` (chain_session.h:3357-3358), carrying its own pybind docstring that already states its limitations ('a cheap, momentum-less cross-check, not a variationally optimal excited state'). Exposing both a plain and a detail form of the same primitive is ordinary API design, and the name is the one the whole surrounding machinery is referred to by throughout chain_session.h's comments, tests/test_infinite_chain.py:1396, tests/test_idmrg_correlator_v3.py:4 and examples/idmrg/idmrg_correlator_python_VS_v3/main.py:9 -- removing it would break those references' referent for zero gain. (4) The one substantive harm alleged -- that it 'bypasses the _warn_if_growth_missed_local_ground_state check' -- does not distinguish it from anything else: that check lives in infinitechain.py, so reaching into the private `_session3` attribute for *any* method bypasses it, along with the rest of the Python layer. There is no route to this binding that is not already a deliberate descent below the public API.

Acting on the suggested fix would be churn with a real downside (removing a documented public binding to satisfy a grep), which is the tell that this is not a finding.



## What each lens covered, and what it ruled out


### `python-backend-parity`

I walked the public `Many_Body_Chain` surface three ways (`itensor_version="python"` vs `itensor_version=3` vs `mode="ED"`) on 4-8 site chains with maxm set large enough to be exact (`maxm >= 2^(n/2)`, cutoff 1e-12..0, nsweeps 20-30), so any disagreement above ~1e-10 is a defect rather than convergence. Phase A (S=1/2, n=6): gs_energy, vev (incl. npow), site/bond/pair correlators, get_correlator, get_correlation_matrix (fast/explicit/full), bond/site/pair entropy, mutual information, get_rdm, get_excited/get_excited_states/get_gap, gs_energy_fluctuation, random_state, and get_dynamical_correlator across KPM/TD/TDZ/CVM/CVM_explicit/ROOTN/EX/SECTOR. Phase B (reduced probe set): Fermionic_Chain n=6, Spinful_Fermionic_Chain (Hubbard) n=3, S=1 n=5, S=3/2 n=4, Bosonic_Chain n=4, Parafermionic_Chain Z=3 n=4. Phase C: set_conserved_sector + promote_to_dense/promote_mps on Sz=0 spin, Nf=3 spinless fermions, Nf=3 Hubbard. Phase D: TDZ, SECTOR, METTS, excited-state orthogonality. Scratch scripts and pickled result tables are in /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/python-backend-parity/ (NOTES.md summarizes).

The headline is that the premise of the lens mostly does not bite: on the *static* surface the pure-Python backend is excellent. python vs v3 vs ED agree to <= 1e-11 for ground-state energy, every vev and correlator, the correlation matrix in all working dmmodes, correlation entropy, site/pair/bond entropy, reduced density matrices, excited-state energies and gaps, the sector-confined ground state and its charge, and submodes KPM/EX/CVM/ROOTN/SECTOR -- across all six site types tested. Where python diverges it does so in exactly two places, both time-evolution-related, and the parity sweep also turned up three defects that are *not* python's (one on v3, two shared/ED). I did not run julia_live (per instructions) and did not run the full pytest suite. I could not root-cause F1 to a specific line: I ruled out cutoff, maxdim, Krylov stopping rule, num_center, the initial state and the MPO (python's TEBD and MPO-Taylor integrators are bit-identical to v3's on the same chain, and norm and energy are conserved to 1e-10), which localizes it to the LR/RL half-sweep composition in pyitensor/tdvp.py but no further. An imaginary-time discriminator was inconclusive (dominated by finite-beta truncation, err plateaued at 5.4e-5 for dt=0.2/0.1/0.05). I deliberately did not report: ED-vs-DMRG KPM amplitude differences (python==v3 to 1e-12, ED's KPM has its own defaults), get_rdm/get_bond_entropy under mode="ED" (documented `_session`-only), reduced_dm_projective NotImplementedError on S=1/S=3/2/boson/parafermion (raises cleanly and identically on all three backends), and a 1-sigma metts_vev deviation at nsamples=40.


### `wavefunction-consumers`

**The lens premise is refuted, and the real ~1e-6 mechanism is a different one.** CLAUDE.md records that `pyitensor/dmrg.py::_lanczos_ground_state`'s eigenVALUE stopping rule caps finite DMRG's own MPS tensors at ~1e-6. Measured directly (n=8 Heisenberg, maxm=20, MPS densified and compared against the ED eigenvector, phase-aligned): the converged wavefunction error is 5.0e-08 after 1 sweep, 2.7e-11 after 3, 8.7e-13 after 6, 6.7e-13 after 10 — six orders below the claimed cap. Mechanism: every local solve warm-starts from the current 2-site tensor and always takes >=2 Krylov steps, returning a Ritz vector strictly better than its input, so the outer sweep is a contraction and the per-solve stopping cap never becomes the fixed point's error. VUMPS differs only because its criterion compares two *independently* solved eigenvectors. Monkeypatching `residual_tol` through (`scratchpad/.../ex3.py`) moves the Heisenberg number from 2.8e-11 to 1.2e-13 at 4 sweeps — a modest win, not a lifted floor — and changes nothing at all on TFIM. What *does* produce a ~1e-6 wavefunction error is the default SVD `cutoff=1e-12` (`manybodychain.py:165`): on TFIM g=0.5, n=8, vector error runs 5.52e-07 / 3.46e-08 / 6.30e-10 / 3.01e-12 at cutoff 1e-12 / 1e-14 / 1e-16 / 0, i.e. ~sqrt(cutoff), while dE stays ~7e-13 throughout. This affects all three backends identically (v2, v3 and "python" all sit at dE=6.9e-13 on that model), so the `vx_lanczos` residual fix is irrelevant to finite DMRG and the mechanism sentence in CLAUDE.md should say "default cutoff", not "Lanczos". Practically, though, no local-observable consumer visibly pays for it: on the same cutoff-limited TFIM state, `<Sx_i>` and `<Sz_i>` still agree with ED to 3e-12 on every backend, because the discarded directions barely couple to local operators.

**Clean negatives (measured, n=8 staggered-field Heisenberg whose DMRG state is exact to ~2e-12, so anything above ~1e-11 would be the consumer's own).** `vev(Sz_i)` 1.2e-14 (python) / 5.3e-12 (v3) / 4.9e-12 (v2); `<Sz_i Sz_j>` 3.9e-14 / 2.5e-12 / 2.0e-12; `get_bond_entropy` 2.5e-14 / 1.0e-11; `get_site_entropy` 2.5e-14 / 1.0e-11; `get_pair_entropy` 1.8e-14 / 1.0e-11; `get_mutual_information` 1.4e-14 / 9.1e-12; `get_rdm` 1.9e-13 / 5.3e-12; `densitymatrix.reduced_dm_projective` 6.8e-15 / 2.8e-12; `overlap`, `summps`, `scale_mps`, `applyinverse` all 1e-14 or better; `get_excited_states(n=3)` energies 1e-15 and `|<e0|ek>|` 1e-15..5e-15 (the overlap-penalty method does *not* compound anything measurable); `effectivehamiltonian.get_manifold` correctly refuses a degenerate cut; `fidelity.get_fidelity` (default `fmode="derivative"`, n>=2) agrees ED 4.4e-14 vs DMRG 1.5e-13; `timedependent.evolve_and_measure` started from the DMRG ground state under its own H holds `<Sz_0>` constant to 2.7e-15 (python) / 8.9e-16 (v3) — v2 drifts 3.7e-3 over t=1, but that is its MPO-Taylor integrator, already the subject of audit finding #2, not the wavefunction. `degeneracy.get_gs_degeneracy` returns 1.0 correctly. Two further low-value executed observations not raised as findings: `get_correlation_matrix()` on a `Spin_Chain` dies with a bare `raise` -> `RuntimeError("No active exception to reraise")` after printing "Unrecognized type ...Spin_Chain" (same shape as audit #12, different site), and `wavefunction.py`'s `Wavefunction` class is dead code (`np.random()` would TypeError on construction; nothing in src/tests/examples ever instantiates it). I did not chase `gap.py::sector_gap`'s `e1 = wf1.overlap(h1*wf1)` (penalized H rather than `h0`) because I had no executed repro for it.

**Coverage/limits.** Everything above is `mode="DMRG"` on `itensor_version` "python", 3 and 2, n<=10, ED cross-check by explicit kron. `julia_live` was excluded per instructions. I did not reach `topology.py` (it operates on Bloch Hamiltonians, not chain MPS) or `dmtk`/`mpsalgebratk/disentangle.py`. All scratch scripts and the self-contained repros are in /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/wavefunction-consumers/ (repro_f1.py, repro_f2.py, repro_f2b.py, repro_f3.py are standalone; ex*.py import lib.py from that directory).


### `dispatch-matrix`

**Part (a), the matrix.** I drove ~34 public `Many_Body_Chain` methods against four columns (`itensor_version=2`, `=3`, `"python"`, and `itensor_version=3` + `sc.mode="ED"`) on a uniform S=1/2 Heisenberg chain at n=4 and again at n=2 (the size that trips mode.py's v3→ED fallback), recording return type/value/exception class/traceback depth per cell (`scratchpad/dispatch-matrix/matrix.py`, `n2.log`). At n=4 the three DMRG columns agree to ~1e-6 on everything they answer, with two exceptions I chased down (findings A and B). The n=2 v3 row — the silent-fallback row — is where the dispatch damage is concentrated (finding C): the same chain hands back an ED `State` from `get_gs()` but an `MPS` from `random_state()`, and `overlap`/`aMb`/`get_rdm` then fail with `AttributeError: 'State' object has no attribute 'cpp_handle'` because `mpsalgebra.py`/`densitymatrix.py` branch on `self.mode` (None, since the fallback was automatic) rather than `self.get_mode()` — audit finding #14's class, fixed for `get_distribution*` and not for these. `itensor_version=2` answers every one of those at n=2. I ruled out several suspects empirically: the sector + ns<3 combination is correct (the ED fallback really does target the sector: Sz=N polarized sector gives 0.25/0.5 against a global −0.75/−1.0), `sectordc._check_backend` raises clearly on v2, on ns<3, on a non-Hermitian H and on `mode="ED"`, and `gs_energy_generalized` raises a clear `NotImplementedError` on v2 and on ED. `get_correlation_matrix`/`get_correlation_entropy` fail identically on all four columns for a spin chain ("Unrecognized type … RuntimeError: No active exception to reraise", `correlationentropy.py:84`) — a bare `raise` for a genuinely unsupported model, so backend-independent and not a dispatch finding.

**Part (b), sites since 2026-08-30.** Read against §4.10: `sectordc.py`'s gating (executed — every refusal path is an explicit, informative raise, and its comment correctly explains why `mode.py`'s sector guard cannot fire for it); `meanfield.py`'s rewrite in 43d1a35 (read — `get_new()` does `sc0.copy()` per iteration, so it is immune to finding A and forwards `mode=`/`**kwargs` to both `gs_energy` and `get_magnetization`; `mixmode` validates its own name); the `7815df1` default-backend change (grepped every subclass constructor — none hardcodes an `itensor_version` default, `default_backend()` is honoured, and `sites.py` creates a session for `"python"` unconditionally); `infinitechain.py`'s `gs_method`/`itensor_version` cross-dispatch (read only — every unsupported combination ends in a named `NotImplementedError`/`ValueError` that prints both coordinates, and `td_dynamical_correlator` rejects stray `kwargs` explicitly; `kpm_finite`'s hardcoded `"python"` window chain is documented in its own docstring). I did **not** run the mechanical `git diff 1b87543..HEAD` grep for the four §4.10 shapes — the matrix run found more per minute. Finding B is a §4.10 "`**kwargs` with no consumer" instance that survived in `edtk/edchain.py`, and finding D is a stale type test rather than a dispatch-order problem.

**Could not reach.** `itensor_version="julia_live"` was out of scope to execute. By reading: `mpsalgebra.overlap_dmrg` and `overlap_aMb_dmrg_MO` have no `julia_live` branch and dereference `self._session`, which `sites.py::initialize` only creates for `itensor_version in (2,3,"python")` — so `Many_Body_Chain.overlap()`/`aMb()` on a Julia chain look like they should fail on a `None`/absent session even though `mpsjulialive/mps.py::MPS` has a working `overlap`/`dot` of its own. Unverified; I did not run it. Also unmeasured: I found no optimization-class finding worth reporting under this lens — the dispatch layers themselves are cheap, and the one caching structure I read (`groundstate.solver_key`/`_session_ham_cache`) is correct.


### `pyitensor-performance`

I profiled the pure-Python backend with cProfile (threads pinned, min-of-N wall timings) across: finite `gs_energy` (L=6..14), excited states (n=2..6), TDVP real-time evolution and `submode="TD"` (nt axis), KPM moments (n axis), METTS `metts_vev` (nsamples axis), bond entropy, `get_correlation_matrix` (all four dmmodes), the four-point tensor (L=4..10), and the infinite-chain paths — iDMRG growth + static observables (maxm=8..64), VUMPS ground state (D=16..32) and `td_dynamical_correlator` (maxm=6..24). Per the plan, exponents were fitted on wall time for the unbounded axes (chi/D, nt, nsamples) and cross-checked structurally on cProfile `ncalls` where the site cap made wall time uninformative. Every finding below was confirmed by a runtime monkeypatch or a standalone re-implementation that reproduces the shipped number and is timed against it; the repository was never modified.

Clean (ruled out, so the next auditor need not redo them): finite `gs_energy` is linear in L with properly incremental left/right environments in `_dmrg_one_sweep`; `tdvp._lanczos_expm_multiply` stops on Saad's residual rather than running the full `niter`, and TDVP/METTS/KPM are all linear in nt/nsamples/n_moments with no per-step MPO rebuild; excited states are linear in n; the four-point tensor's default `ctmode="sweep"` scales ~L^1.7, not L^4; `mpobuilder.to_mpo` is called exactly twice per `gs_energy` (the `is_hermitian` double-build is real but only ~20% at these sizes and is already recorded as deliberately deferred). One thing measured but not elevated: `Chain.correlation_matrix` (reached as `dmmode="full"`, which is the *default* in conserved-sector mode) builds one full MPO per (i,j) pair — 36 `to_mpo` calls at L=8, O(L^2) builds — but `mpscpp3/chain_session.h::correlation_matrix` has the identical shape, so it is a shared design choice rather than a pyitensor-specific miss, and the default non-sector `dmmode="fast"` path is fine. I did not reach `idmrg_excitations`' own public entry points (`excitation_energies`, `dynamical_structure_factor`, `spectral_weights`) with timings; they route through the same `_op_transfer_matrix` as finding 2, so the same copy cost very likely applies there, but I have marked that part of finding 2 explicitly as unverified.


### `cpp-v3-completeness`

What I covered. (a) mpscpp3 vs mpscpp2: v3's pybind surface is a strict superset of v2's (42 defs vs 82; no v2-only method), and for every shared method the `py::arg` names and defaults are byte-identical except `itensor_smoke_test`. Numerically I cross-checked, on 5-6 site chains against ED and pyitensor: gs_energy, vev profiles, `correlation_matrix`, site/pair entropy, mutual information, correlation entropy, `get_excited`/`get_excited_states` (incl. a degenerate pair), gs energy fluctuation, gs degeneracy, `general_kpm` moments, operator norm/`aMb`, `exponential_apply`, and the dynamical correlator in submodes KPM/CVM/EX/TD. Everything agrees (v2-v3 ≤1.6e-10 for CVM, ≤5.6e-7 for EX; all ED cross-checks ≤6e-12). The two v2-v3 gaps I saw are both documented and out of scope: TD differs by 7.6e-4 (TDVP on v3 vs MPO-Taylor on v2) and `exponential_apply` is the deliberately-reproduced 3rd-order-Taylor/`H2` bug on both.

(b) `ic_*`/`vx_*` vs `pyitensor/idmrg.py` / `vumps.py`: I could not find a static observable that diverges. Cross-checked v3 vs "python" for both `gs_method`s on a gapped dimerized n_uc=2 XXZ+transverse-field cell (unique state, so the comparison is well-posed): `e0`, `vev` of Sx/Sy/Sz at both sublattices, and `correlator` for SzSz/SxSx/SySy/SxSz/SzSx at r=0..4 — idmrg agrees to ≤1.2e-7, vumps to ≤3e-15. Repeated on a **mixed S=1/2 + S=1 cell** with p_i=0 and p_i=1 (idmrg ≤2.6e-8, vumps ≤2e-14), and on `local_excitation_gap(window=0)` (≤1.2e-6 at maxm=8). Perron selection, the degenerate-dominant-eigenvalue guard and the `C C†` redundancy fallback are all present in the C++ (`vx_perron_reorder`, `vx_check_perron_nondegenerate`, `vx_bond_fixed_points`). I also ran the one test-parametrize gap I found: `test_infinite_chain.py::test_fermionic_correlator_matches_free_fermion_exact` covers `("python","idmrg")`, `("python","vumps")`, `(3,"vumps")` but **not `(3,"idmrg")`**, while the energy test one function above it (`test_free_fermion_energy_density_matches_band_integral`) does cover `(3,"idmrg")` — the exact "energy pinned, observable not" shape the brief warns about. I ran the missing combination against the exact free-fermion one-body density matrix: it passes to 5.6e-9. So that gap is in the pin, not in the code; one parametrize entry closes it.

Timing: v3 `gs_energy` is *faster* than v2, not slower (n=30, maxm=60, min of 3, one core, threads pinned: v2 14.56s vs v3 3.35s), so I have no optimization finding to report. What I could not reach: `julia_live` (forbidden), `itensor_version=3` VUMPS excitation energies on a strongly anisotropic XYZ+field model (both backends refuse it identically with an unconverged mixed-transfer-tensor error, so no divergence to see), and `td_dynamical_correlator` at t>0 (ROADMAP records an existing python-VS-v3 example that already compares every x). The four findings below are all in the Python dispatch layer that decides *which* backend answers — which is where the lens's remaining holes turned out to live.


### `ed-and-operators`

**Covered and ruled out (all executed).** Fermionic sign conventions are clean everywhere I could reach. Spinless `Fermionic_Chain`: 21 probes (`<Cdag_i C_j>` both orderings, `<C_i C_j>`, `<Cdag_i Cdag_j>`, `<C_i Cdag_j>`, mixed 3-factor `N_k Cdag_i C_j`, and six 4-point products) on a 4-site chain with a random complex Hermitian hopping matrix, an antisymmetric pairing matrix and a nearest-neighbour interaction, against an independent dense Jordan-Wigner reference I wrote from scratch — ED agrees to 1e-16..1e-18 and the compiled v3 DMRG to 1e-13..1e-9. Spinful: the same style of probe (`<Cdagup_0 Cdn_2>`, `<Cup_0 Cdn_2>`, `<Cdagup_0 Cdagdn_2>`, `<Sx_0 Sx_2>`, `<Delta_0 Delta_2^dag>`, a 4-point on-site-pair correlator) on both `Spinful_Fermionic_Chain` (interleaved) and `Spinful_Fermionic_Chain_Native` (Electron sites) against the same dense reference built on flat modes 2i=up/2i+1=dn — ED exact to 1e-16 on both classes; the native DMRG's ~1e-6 deviations track its own E0 error of 6e-6 (bond truncation at maxm=30 on a dim-64 chain), i.e. convergence, not a sign. Z3 parafermions with a chiral (complex-phase) hopping and transverse field: `<Tau_0>`, `<Sig_0 SigDag_2>`, `<Chi_0 Chid_2>`, `<Psi_0 Psid_3>` ED vs v3 DMRG agree to 5e-7..5e-6 (DMRG convergence). Also verified clean: a second `set_hamiltonian()` correctly drops the cached ED object on fermion/spin/boson chains (no accumulation into `MBFermion.h`); dimension-1 and dimension-2 ED sectors (`set_conserved_sector(Nf=0)`, `(Nf=1)`, `(Nf=ns)`) return correct energies and vevs without crashing; `multioperatortk/charge.py::charge_components` correctly splits `Cdag` and `C+Cdag` into definite-Nf pieces.

**What I could not reach / did not pursue.** `mixedchain.Mixed_Spin_Fermion_Chain` has no ED backend at all (its own docstring says so), so there was nothing to audit on the ED side. `itensor_version="julia_live"` was forbidden, so the user guide's identical claim about the Julia backend and non-default `maxnb` is untested. I noticed but did not pursue two things: `edtk/dynamics.py::dynamical_correlator_kpm`'s effective broadening is not the `delta=` the caller passes (it sets the polynomial count instead), so `mode="ED", submode="KPM"` — the *default* submode — differs from every other submode by a resolution-dependent factor that I could not cleanly separate from KPM's inherent approximation at 4 sites, and I make no claim about it; and `EDchain.sector_nonconserving`'s symmetry-violation threshold is *relative* to the largest matrix element (`tol*max|data|`), so a small charge-violating term riding along with a large conserving one would be silently projected rather than refused — plausible but I did not construct a case. `pyzn/zn.py`, `pyfermion/states.py::four2many`, `sympymultioperator.get_dagger` and `pyboson.SpinBosonChain` are all unreferenced dead code / self-declared stubs and are excluded deliberately.


### `recent-commits`

Covered: every commit since 2026-08-25 (`git log --since=2026-08-25`, 45 commits), reading the diffs of 43d1a35 (meanfield rewrite), 7815df1 (default backend), e448699 (pyitensor MPO finite-state machine), 54580b0/3593a2c (release + doc sync), deb6bf8, 71ba8eb and a2eb46e (VUMPS long-range dispatch + v3 redundant-bond-dimension fix). All work ran from /tmp/claude-1000/.../scratchpad/recent-commits with threads pinned and PYTHONPATH forced to this checkout; nothing in the repo was touched.

RULED OUT (executed, clean): (a) The e448699 MPO automaton — 600 randomized cases (n=1..5, spin-1/2/spin-1/spinless-fermion/electron sites, 1-4 factors per term, deliberately *unsorted* factor order, exact duplicates, terms differing only by coefficient, zero and complex coefficients, JW-strung odd-parity terms); `_automaton_mpo`, `to_mpo` and `AutoMPO.dense_matrix()` agree to <1e-10 in every case. The only mismatches were `_sum_of_term_mpos` on 1-site chains, which is the known pre-existing bug the rewrite fixed. Both "load-bearing rules" (coefficient on the transition into F; accumulate into F vs assign structural) hold: structural writes never target column `col_F` and no term's row is ever the F row. (b) 43d1a35 meanfield — `tests/test_meanfield.py` 8/8 in 1.4s; the documented default (DMRG, itensor_version=3) path runs end to end; `mixmode="broyden"` runs; grep over examples/, docs/*.md, docs/*.tex, README.md, tests/, src/ found no surviving caller of `set_fields`/`set_swave_pairing`/`set_exchange` other than deliberate comments. (c) 54580b0 packaging — built the wheel in a scratch copy: 0 `mpscpp*` entries, `dmrgpy/juliapkg.json` + `mpsjulialive/*.jl` + `algebra/fortran/*.f90` present, `compileall` clean (2 pre-existing SyntaxWarnings), and importing the unzipped wheel with the checkout off PYTHONPATH gives the correct "defaults to python" warning, `default_backend()=="python"`, and DMRG==ED to 1e-15. (d) Four assert-carrying examples from CLAUDE.md (tdvp_VS_ED_time_evolution, backend_switch_consistency, static_correlator_VS_ED, entanglement_entropy_VS_ED) all print TEST PASSED. (e) 7815df1 default-backend sweep — with `cppext._backends[2]=_backends[3]=None`, all 7 chain classes resolve to `"python"`, and Spinful_Fermionic_Chain_Native / Bosonic_Chain(maxnb=5) / Parafermionic_Chain(Z=3) all agree with ED to ~1e-13 on that backend.

Demoted to notes rather than reported: `meanfield.spinchain_meanfield(..., maxite=0)` raises `UnboundLocalError: cannot access local variable 'error'` (for-else references a name the loop never bound); several stale comments left by the window's commits — `Spinful_Fermionic_Chain_Native`'s docstring says pyitensor "does not implement" a native spinful site (it does, and agrees with ED to 1e-13), `bosonchain.py:17` says pyitensor "still only understand[s] the plain 104 code" (`pyitensor/sites/siteset.py:35` routes 102..199), and `pyproject.toml`'s trailing NOTE still says "dmrgpy transparently falls back to ED", which 7815df1 changed. CLAUDE.md also contradicts itself on sector+ED ("mode.py raise rather than fall back to ED" vs "falling back to ED is now correct"); `mode.py:32` shows only non-(3,"python") DMRG raises, so the second statement is the code.

Could not reach: `itensor_version="julia_live"` (excluded by brief), and I did not run the full suite.


### `docs-examples-drift`

Covered, mechanically first then by execution. (a) docs/user_guide.{md,tex} + docs/documentation.{md,tex} against the code: an AST/regex sweep matched every `.method(` mentioned in either document against every name defined under src/ (zero missing names — the 3593a2c/43d1a35 user-guide audits closed that class); a second sweep matched every documented `kwarg=` against the real `inspect.signature` of `Many_Body_Chain`; a third compared every "default `x=v`" claim in prose against every default value parsed out of every `def` in the tree. That third sweep is what produced findings 9 and 10. Both `.tex` files compile: `pdflatex -interaction=nonstopmode` gives 0 errors, 0 undefined references, 2 overfull hboxes in user_guide.tex and 0 in documentation.tex (the same counts 3593a2c recorded), so I report no LaTeX breakage. A full backtick-identifier md-vs-tex parity diff was too noisy to be useful (verbatim/lstlisting blocks defeat it); the one md/tex divergence I do report (finding 11) came out of executing the doc's own snippet instead. (b) src/ docstrings and kwargs: an AST pass over every FunctionDef outside ITensor/ listing parameters and `**kwargs` never referenced in the body (69 hits), plus `grep -rnE '^\s*raise\s*$'` (40 bare-raise sites). I executed the top candidates from both lists; findings 1, 2, 4, 9 and 11 come from there. (c) examples/: 30 scripts run with the mandatory thread pinning, `PYTHONPATH` to this checkout and `MPLBACKEND=Agg` — the assert-carrying VS_ED regression templates CLAUDE.md names, the v2_VS_v3_* comparisons, the sector/idmrg/long-range/backend-switch scripts touching recently-changed code, and all 13 of `readme_examples/`. Two non-zero exits: `boson_models/v2_VS_v3_boson` (finding 6) and `readme_examples/edge_correlator`, the latter only my own 120 s timeout on a legitimately slow n=20 S=1 KPM run, not a defect. Plot-rule check was a grep for scripts importing matplotlib but never calling a drawing method (13 hits; 10 are the legitimate single-scalar ED-vs-DMRG exception, 3 are findings 7 and 8).

Not reached / ruled out: nothing was run on `itensor_version="julia_live"` (per instructions), so the Julia mirror modules' own drift is unaudited — note `mpsjulialive/dynamics.py:129` accepts and drops a documented `deconvolve=`, which I did not execute and therefore do not report. I checked each finding against docs/audit_2026_08_hole_hunt.md (including its Refuted section), ROADMAP.md and the docs/known_issue_*.md files; none of the eleven is already recorded — finding 1 is adjacent to audit #4's "kwarg with no consumer" class and to documentation.md §4.6's note that `npow=` raises on julia_live "rather than silently misbehaving", which is precisely what it does on ED, and finding 4 is the same `else`-serving-two-purposes shape audit #12 fixed at one site and left at three others. The repository is read-only and was left byte-clean: `python clean.py` from the repo root, then reverting one PNG the examples overwrote and deleting two they created; `git status --porcelain` is empty. All scratch under /tmp/claude-1000/-home-joselado-Documents-programs-dmrgpy/26f5986a-45f8-441f-88b8-2d3f4e95216b/scratchpad/docs-examples-drift/ (FINDINGS.md, examples.log, boson.log, and one .py per repro).

