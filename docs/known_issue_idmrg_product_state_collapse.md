# Fixed issue: iDMRG collapsed into a product state and reported `converged=True`

Filed from the same external physics session (a triangulene heavy-fermion / Kondo-lattice
calculation) that reported the `local_excitation_gap` deflation bug — that one is **fixed**
(see `Chain::idmrg_local_excitation_gap`'s own comment and
`tests/test_idmrg_local_gap_deflation.py`); so is this one, and the file is kept as the
record of what the failure looked like and what it took to break it.

## Summary

On an infinite spinful (c,f) Kondo chain (`Infinite_Many_Body_Chain([1,1],
itensor_version=3)`, `gs_method="idmrg"`, `maxm=16`, `maxiter=200`), sweeping a uniform
chemical potential `mu` past roughly `mu = -0.20` makes `Chain::idmrg_ground_state`
return the **empty product state** — `n_tot = 0`, every correlator identically 0,
`e/site = +0.064171 = U/8` exactly (the `U*(Nup-1/2)*(Ndn-1/2)` term on an empty f site)
— and report `converged=True` after ~0 s, against ~3500 s for a healthy point at the same
bond dimension. Symmetrically, large positive `mu` returns the fully-filled product state
(`n_tot = 4`, `e/site = -0.540133`, again the exact product-state value to 6 digits).

Three independent arguments that this is wrong rather than physics, all from the
reporting session:

1. **Exact free-fermion reference.** At `U=J=0` the band integral gives `n_tot = 0.656`
   at `mu=-0.30` (the lower band spans `[-0.517, -0.019]` after the mid-gap shift, so
   occupied states certainly exist). A repulsive `U` and an AF `J` cannot empty a band the
   free model fills.
2. **Variational.** The half-filled state — whose energy at any `mu` is known exactly,
   since a uniform shift only moves it by `-mu*n/2` — gives `e = +0.034655` at `mu=-0.30`,
   below the `+0.064171` the solver returned. The solver is beaten by a state already in
   hand.
3. **The other solver finds a much better state at identical parameters.** `pyitensor`
   VUMPS (`itensor_version="python"`, `gs_method="vumps"`, same `maxm=16`) returns
   `n_tot = 1.1696` (`n_c = 0.7169`, `n_f = 0.4527`) at `e/site = +0.005052` — 59 meV/site
   below the growing algorithm's answer. So the product state is an attractor of this
   *algorithm*, not of the physics.

## Why it is not a convergence-tolerance problem

The obvious first guess is the convergence test: a product state is an exact fixed point,
so the energy-density finite difference between consecutive macro-iterations is 0 and
`etol` trips immediately — hence `converged=True` in ~0 s. But setting `etol=0`, so the
test can never trip, changes nothing except the flag: the answer is identical and the run
still finishes in 2-3 s rather than ~3500 s.

That runtime *is* the diagnosis. Once the state collapses to a product state its bond
dimension collapses to 1; the SVD at every subsequent micro-step then has a single nonzero
singular value, and nothing in the growth loop can grow `chi` back. It is an absorbing
state, not slow convergence.

This is the same family as the documented random-start metastability of this backend
(`mpscpp3/get_sites.h` builds every site with `ConserveQNs=false` and
`chain_session.h`'s `default_mps()` starts from `randomMPS`), but with a sharper failure
mode: an unlucky start is one thing, an absorbing state the algorithm cannot leave is
another.

## It appears stochastic in the random start, not deterministic at a given mu

Worth knowing before anyone validates a fix against this: an independent reproduction here
at `mu=-0.30`, `maxm=16`, `maxiter=200` — the same point the reporting session got the
empty state at, in ~0 s — did **not** collapse. It ran normally for over half an hour, the
same order as a healthy point at that bond dimension.

That is consistent with the mechanism rather than against it: `default_mps()` starts from
`randomMPS` with no seeding, so whether a given run falls into the absorbing product state
is a property of that run's start, not of `mu` alone. It does mean "reproduce at
`mu=-0.30`" is not a reliable test, and that any fix has to be validated over repeated
runs (or with a seeded start) rather than a single one.

## The fix: White noise, applied on demand

Implemented in BOTH backends -- the `"python"` reference
(`pyitensor/idmrg.py`'s `_noise_perturbed_split`) and its ITensor v3 C++ port
(`Chain::idmrg_noisy_isometry`) -- because both collapse identically on the model
below, so this was a shared gap in the algorithm rather than a porting bug.

**It takes two halves, and the second is the one that does the work here.**

1. *Enlarge the basis.* White's density-matrix perturbation, modelled on ITensor's own
   `LocalOp::deltaRho`: decompose `rho + noise*drho` instead of taking a plain SVD, with
   `drho` built from the individual MPO channel components `W^a|theta>`. Using the
   channels rather than `H|theta>` is essential -- `H|theta>` is *parallel* to `|theta>`
   at an eigenstate and would open no new direction at all, while a channel that creates
   a particle takes the vacuum to a one-particle state. This loop extends `HL` and `HR`
   in the same micro-step, so `U` and `V` come from two independently perturbed density
   matrices (a finite sweep needs only one direction); they are then no longer an SVD
   pair of `theta`, which nothing here requires.
2. *Unpin the solver.* Necessary but **not sufficient**: measured directly, the
   density-matrix term alone grew `chi` away from 1 and the state still never moved
   (energy identically 0 for every macro-iteration). A product state like the vacuum is
   an exact eigenstate of the FULL Hamiltonian, so it stays an exact eigenvector of the
   local effective Hamiltonian in any basis, however enriched -- and a Krylov solve
   started exactly on an eigenvector breaks down on its first step and hands the same
   vector back. A random admixture of weight `O(noise)` on the local solve's start vector
   breaks that pinning, and the enriched basis is then what gives the freed solve
   somewhere lower to go.

**The schedule is demand-driven, not a fixed number of iterations.** Two measurements
forced this. First, escape depends on *duration* rather than magnitude: `noise=1e-4` over
12 macro-iterations escapes while `noise=0.1` over 10 does not, magnitude-independent
across three orders of magnitude. (The physical reading: the growing chain has to reach a
length at which a particle-carrying state beats the vacuum, and the noise only has to
hold the door open until it does.) Second, a fixed "noise for the first N iterations"
default perturbed *healthy* runs enough to regress a free-fermion correlator test by
1.3e-6 against its 1e-6 tolerance.

So noise arms only while the state genuinely is a product state, keyed on the purity of
the **noise-free** reduced state (`tr rho^2`, exactly 1 for a product state and strictly
less otherwise), and runs for `_NOISE_TAIL` (6) further macro-iterations after that
clears so it cannot switch off the instant the basis opens and let the state fall
straight back. An already-entangled model never arms it and is bit-for-bit unaffected.
`noise_iters` (40) caps the total, so a model whose ground state genuinely IS a product
state -- a field-polarized chain -- stops re-arming and still finishes on a clean,
unperturbed tail. The `etol` test is suppressed whenever noise is active, so a perturbed
state can never be what a run reports as converged; that also removes the
"converged=True in 0 s" pathology by construction.

Both knobs are exposed as `Infinite_Many_Body_Chain.noise` / `.noise_iters`;
`noise = 0` disables the mechanism entirely and restores the pre-fix numerics exactly.

## Verification

`tests/test_idmrg_product_state_trap.py`, on both backends. The model there is the
cheapest one that reproduces the trap deterministically *and* has a closed-form answer:
a dimerized spinless-fermion chain (`t1=1.0`, `t2=0.4`) at a uniform `mu=+1.25`, where
the band integral gives `e/site = -0.01648450` at `n/site = 0.166092`.

| | before | after |
|---|---|---|
| `"python"` | `e=0`, `n=0`, 0.0 s, 4/4 runs | `e=-0.016251`, `n=0.146` |
| `itensor_version=3` | `e=0`, `n=0`, 0.0 s, 6/6 runs | `e=-0.016481`, `n=0.158` |
| exact | | `-0.01648450`, `0.166092` |

The third test in that file is the case the fix could plausibly damage: a field-polarized
XX chain whose true ground state IS a product state, so it arms the noise and keeps
re-arming it. It still returns `<Sz> = 0.5000000000` and `<Sz(0)Sz(3)> = 0.2500000000`
exactly, because the cap ends the noisy phase and the clean tail snaps the state back.

## What was NOT done

The reporting session also floated a `chi=1` detector and letting `gs_energy()` accept an
initial state so a caller could walk a parameter adiabatically. Neither is implemented:
with the trap itself fixed, a detector has much less to catch, and initial-state seeding
is a separate API question rather than a bug fix.

## Practical scope today

`itensor_version=3` + `gs_method="idmrg"` is reliable on the half-filled plateau for this
model (it reproduces `n_c = 1.0009`, `n_f = 0.9991` and matches a finite-chain correlator
benchmark) but cannot reach hole doping beyond about `mu = -0.20`. `pyitensor` VUMPS does
reach it. Combined with `docs/known_issue_v3_vumps_no_return.md`, no single
(backend, solver) pair covers the whole parameter range of this Hamiltonian.
