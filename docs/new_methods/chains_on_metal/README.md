# Spin chains on metals: several baths, and an Anderson level

This folder documents two extensions of `ysr.spin_chain_on_metal`, the
Heisenberg chain on a metal whose STM spectrum the TASK and NRG routes of
`ysr.get_stm_spectrum` compute. The first puts several spins of the chain
on the metal, each with a Wilson chain of its own
(`spin_chain_on_metal(ns, n, site=[...])`). The second makes the spin on
the metal an Anderson level, an orbital with a charge, instead of a Kondo
spin (`ysr.anderson_spin_chain`). The tests are
`tests/test_ysr_chains_on_metal.py`, the examples
`examples/kondo/spin_dimer_two_metals` and
`examples/kondo/anderson_level_on_spin_chain`, and the scripts of the
runs on Triton are in `runs/`.

Both came out of the question whether the one-crossing approximation
(OCA) of Jacob's STM calculations (arXiv:1501.02449, arXiv:2104.02503)
could be adapted to tensor networks. The answer was no: OCA over the
states of a whole chain is not size consistent when several sites touch
the metal, and truncating those states removes the virtual charge states
that renormalize the spin excitations. What OCA is used for, charge
fluctuations of a level on a metal and several sites on a metal, is done
here with the metal inside the MPS, which is size consistent by
construction and controlled by the bond dimension and the length of the
Wilson chains.

## 1. The models

Several spins on the metal. With `site` a list, every spin of the list
is Kondo coupled, `JK S_i.s_0^(i)`, to the first site of its own Wilson
chain, the same chain for all, and `JK` may differ from spin to spin.
The baths are independent. This is the model of spins far apart on the
scale of the Fermi wavelength; on a common metal the number of channels
is the rank of the coupling matrix of the spins to the metal (Eickhoff
and Anders, arXiv:2008.09034), which this model does not contain.

An Anderson level. The spin `site` of the chain becomes an orbital d,

    H_d = ed N_d + U n_{d up} n_{d dn} + tD sum_s (d^dag_s f_{0 s} + h.c.),

that hops onto its own Wilson chain, and the exchange `Js S_d.S_i` with
its neighbours acts on S_d = d^dag (sigma/2) d, so it is there when the
level holds one electron and vanishes when it is empty or doubly
occupied. The other spins stay charge frozen. The tip puts its electron
on the level itself, `X^dag_s = d^dag_s + q f^dag_{0 s}`, and at q=0
dI/dV is the spectral function of the level summed over spin, the ideal
STM limit of Jacob (arXiv:1805.10624). The hybridization is
Gamma = pi rho tD^2 with rho = 1/(2D), and on the mean-energy mesh
`tD*sqrt(lambda_factor(Lambda))` on the chain stands for tD in the
continuum, the factor of `island_chain`'s t_c.

## 2. Where the Wilson chains go

A cut of the MPS carries the entropy of what it separates, so the order
of the sites sets the bond dimension a number needs, not the number. We
compared four orders (`ysr.METAL_ORDERS` and one by hand):

- linear: all spins, then the chains one after the other;
- comb: every chain right after its spin;
- sides: the spins in the middle, the chains of the left half reversed
  on the left end and those of the right half on the right end, stacked
  with the chain of the spin nearest the end innermost;
- for two baths, sides is the order with one chain at each end.

The measurement: 5 spins, Js = 0.5, JK = 1, Lambda = 2, 13 Wilson sites
per chain, v3 with 20 sweeps, energy above the lowest one found for the
same model (the script is `runs/order_probe.py`).

| baths (spins) | order | maxm 32 | maxm 64 | maxm 128 |
|---|---|---|---|---|
| 1 (2) | linear | 3.0e-4 | 6.1e-6 | 0 |
| 1 (2) | comb | 9.9e-3 | 1.3e-3 | 4.3e-5 |
| 2 (1, 3) | sides | 4.6e-4 | 1.9e-5 | 0 |
| 2 (1, 3) | comb | 1.2e-2 | 1.1e-3 | 1.2e-4 |
| 2 (1, 3) | linear | 2.5e-2 | 2.5e-3 | 3.5e-4 |
| 3 (0, 2, 4) | comb | 2.3e-2 | 5.1e-3 | 0 |
| 3 (0, 2, 4) | sides | 2.9e-2 | 6.9e-3 | 8.1e-4 |
| 3 (0, 2, 4) | linear | 1.5e-1 | 4.9e-2 | 1.0e-2 |

With four baths (spins 0, 1, 3, 4) at maxm 64, sides is 4.1e-3 above
comb and linear 9.7e-2. So linear, whose cuts inside a chain also cut
every spin from the chains still to come, grows fast with the number of
baths; sides halves that growth; comb does not grow with the number of
baths, since a cut inside a chain carries that chain and one bond of the
spin chain, but every spin bond it jumps costs (4.8i of
`docs/development.md` measured 2 to 4 times the bond dimension for one
chain).

JK = 1 is strong coupling (rho J = 0.5 on the chain). At JK = 0.3 sides
leads with three and four baths too:

| baths | order | maxm 64 | maxm 128 | largest entropy at 128 |
|---|---|---|---|---|
| 3 (0, 2, 4) | sides | 1.4e-3 | 0 | 1.93 |
| 3 (0, 2, 4) | comb | 1.9e-3 | 2.1e-4 | 2.15 |
| 4 (0, 1, 3, 4) | sides | 1.4e-3 | 0 | 1.79 |
| 4 (0, 1, 3, 4) | comb | 2.6e-3 | 6.2e-4 | 2.35 |

Comb has a second problem. At JK = 0 the exact energy is the Heisenberg
chain plus free Wilson chains. With three baths at maxm 128 comb stopped
0.177 above it and with four baths 0.36 above it, while sides was within
5e-9 (three baths) and 6.4e-9 (four), and linear within 5e-9 (three).
With nothing coupling the spins to their chains, the spin bond that
jumps a chain is carried by no bond of the start state, and two-site
DMRG does not build it.

So the default for several baths is sides, at any number of baths. Comb
leads only at JK = 1 with three and four baths, which is strong
coupling, and it stalls where a spin is weakly coupled, a stall that a
sweep in maxm does not show; an explicit `order="comb"` with a spin of
zero coupling warns. The deficits of sides are deficits of convergence,
which fall with maxm.

## 3. TASK with several baths

The one-site tangent space of TASK misses part of the vertex X^dag|GS>
when there are several baths, also where the state is exact. Each chain
keeps its own charge, since only exchange couples it to its spin, so a
block that holds a whole chain has a definite charge in the ground state
and its Schmidt basis holds only that charge; the vertex adds an
electron to the chain and moves the spins at the same time, which a
change of one site with the other bases fixed cannot do. Measured
against ED at saturation (v3, maxm 256 to 512, TASK error over the peak):

| chain | order | tip on the first spin | tip on the second |
|---|---|---|---|
| 2 spins, 3 sites per chain | linear | 1.7e-4 | 9.6e-2 |
| 2 spins, 3 sites per chain | sides | 9.6e-2 | 9.6e-2 |
| 2 spins, 3 sites per chain | comb | 7.8e-3 | 9.6e-2 |
| 3 spins (0, 2), 3 sites | linear | 6.2e-6 | 3.0e-2 |
| 3 spins (0, 2), 3 sites | sides | 3.0e-2 | 3.0e-2 |

with a larger Krylov space changing nothing (nkry 100 and 400 gave the
same 8.3e-3 on a comb chain), and KPM on the same state right (2.6e-4
of ED's KPM, the two tips equal to 5e-12 by mirror symmetry). The bond
expansion of TASK, `expand_bonds=8, expand_rounds=3`, adds the missing
directions: 1e-11 to 3e-12 for both tips in sides and comb on the
first chain of the table. `expand_bonds=4, expand_rounds=2` was not
enough (4e-3 to 5e-3). It is not a complete cure in every order: in
linear the tip on the first spin keeps 1.7e-4 of the peak, the same at
nkry 400 and with `expand_bonds=16`, while `pad_bonds="full"` is exact
(3e-12). So the missing part is not only in X^dag|GS> and its first
Lanczos step, which is what the expansion adds, but also further along
the recursion.

On a long chain the start vector itself leaves the tangent space. The
weight of X^dag|GS> that TASK misses (`taskdmrg.get_poles`, its
`missing_weight` over the sum rule), on 2 spins with Js = 0.3, JK =
0.3 A_Lambda, Lambda = 2, maxm 128, against the length n of each Wilson
chain:

| n | one bath (under spin 0 or 1) | two baths, sides | two baths, linear, tip spin first |
|---|---|---|---|
| 3, 5 | 0 | 0 | 0 |
| 9 | 0 | 1e-3 | 0 |
| 15 | 0 | 0.622 | 0 |
| 21 | | 0.95 (0.104 of 0.110) | 1e-14 (energy 1.2e-5 above sides) |

With 8 spins, two of them on a chain of 9 Wilson sites, maxm 128, so that
the block of spins (2^8 states) is cut: tip-first linear misses 4.5e-7
to 4.8e-7 with the tip on the end spin and 1.8e-6 to 3.5e-6 with the tip
on spin 3 of 8, of each of the four vertices.

and with `expand_bonds=8, expand_rounds=3` at n = 21 it is 0.018 at
maxm 128 and 0.020 at maxm 64. The vertex of the tip on spin j acts on
spin j and on the first site of its chain. In the linear order with the
tip spin listed first, that site comes right after a block that holds
all the spins and nothing else, whose left basis holds every spin state
the vertex makes; in sides the chain of the tip sits at the end of the
MPS, and the bond between that whole chain, of fixed charge, and the
spins cannot carry an electron added to the chain together with the
flip of the spin.

So the order that is best for the energy is not the one that is right
for the spectrum. What the code does, decided on 2026-11-18 by Fable,
the reviewer the author named for the night:

- The default order stays sides.
- Plain TASK on several Kondo baths raises before any solve, unless the
  order is linear and the tip spin is the first of `site=`; the message
  names that rebuild (and the reordering of a per-spin `JK` with it) and
  the bond expansion.
- Linear with the tip spin first runs with a warning, until three things
  are measured: no missing weight at n = 21, the same on 5 or more spins
  with 2^ns above maxm, and dI/dV within 1e-2 of the peak of the
  variational correction vector. With three baths or more the warning
  stays and names the cost of the linear order in the energy.
- With bond expansion the call runs in any order, and raises when the
  vertices miss more than `ysr.MISSING_RAISE` = 0.1 of their sum rule (a
  choice between the 0.02 of 8/3 and the 0.62 to 0.95 of no expansion,
  not a measurement).
- Anderson levels keep a warning: the hopping moves charge onto the
  chain, so the argument does not carry over, and nothing was measured.
- The expansion stays opt-in: it leaves 0.018 at n = 21, costs 2612 s at
  maxm 64 on the dimer, and at maxm 256 on three and four baths its
  tangent space passes the cap of 4e6 (`task.MAX_PADDED_DIMENSION`) and
  TASK raises (three of the four maxm 256 chain runs stopped there).

The correction vector is not a safe reference here either. With the
default conjugate gradient (`cvm_solver="cg"`) the correction vector of
these vertices stopped with residuals of 1e-2 to 5e-2 at positive bias
even on 8 sites at saturation (maxm 128, no truncation), and dI/dV was
0.115 of the peak off ED; `cvm_solver="variational"` (dynamical DMRG,
valid since the pair is (X, X^dag)) was 5.8e-9 off and 2.8 times
faster. The references of section 5 are variational
(`runs/reference.py`).

## 4. The Anderson level: checks

- One level and no other spin is `wilsonchain.anderson_chain`: the same
  energy and the same spectral function to roundoff on ED.
- At large U the level is a Kondo spin (Schrieffer and Wolff): at
  JK = 8 tD^2/U fixed, the low-energy levels approach those of
  `spin_chain_on_metal` and the dI/dV at |eV| < 0.6 approaches
  64 (tD/U)^2 times the Kondo one, the electron-like vertex d^dag
  becoming the composite one. On 2 spins and 3 Wilson sites the ratio
  over 64 (tD/U)^2 is 0.61, 0.75, 0.85, 0.92 and 0.96 at U = 4, 8, 16,
  32 and 64, the shape 15% and 7.9% off at U = 32 and 64 and the levels
  8.5e-3 and 4.5e-3 off, every error falling as 1/U.
- TASK against ED at saturation, the level in the middle of three spins,
  off particle-hole symmetry and with q = 0.2 (test).
- The NRG route (head and tail of `nrgmps`) with a complete hand-over
  basis and no truncation in the tail is the ED spectrum of the chain,
  for one level and for a level next to a spin, off particle-hole
  symmetry and with q (test).
- The Friedel sum rule. One level at ed = -U/2, U = 0.5, Gamma = 0.05
  (tD on the mean mesh at Lambda = 4), a head of 12 Wilson sites, the
  default hand-over (m = 4) and tail (keep 800): pi Gamma rho_d(0) per
  spin is 1.064 on one mesh and 1.000 averaged over two (maxm 160),
  0.998 over four (maxm 200). On the two meshes rho_d has fallen to 0.84
  of its height at w = 1e-3 and to 0.29 at 1e-2, against a Haldane T_K
  of about 2.6e-3. One mesh takes 10 to 20 minutes on one core: of 434 s
  on a head of 10 sites at maxm 128, 323 s are the TASK poles of the
  head, 98 s the six DMRG solves of `head_state` (three charges, each
  with its seed) and 11 s the NRG tail.

## 5. Runs

The scripts are in `runs/`; the night of 2026-11-18 ran them on one core
each, locally and as Slurm arrays on Triton.

An Anderson level at the end of a chain of spins, NRG route
(`runs/anderson.py`; U = 0.5, Gamma = 0.05, Js = 0.01, Lambda = 4, a head
of 12 Wilson sites at maxm 200, four meshes). pi Gamma rho_d(0) per spin:

| spins | ed | pi Gamma rho_d(0) | rho_d(+0.01) / rho_d(-0.01) |
|---|---|---|---|
| 1 | -U/2 | 0.998 | 1 |
| 3 | -U/2 | 1.0006 | 1 |
| 3 | -U/3 | 0.989 | 1.12 |
| 3 | -U/6 | 0.877 | 1.61 |

Away from ed = -U/2 the height falls below the unitary value, as the
Friedel sum rule pi Gamma rho_d(0) = sin^2(pi n_d / 2) says for an
occupation below one, and the resonance leans to positive bias. The
Kondo spin of `spin_chain_on_metal` at the Schrieffer-Wolff coupling of
the symmetric level, JK = 8 tD^2/U, gives A_X(0) = 1.354 (one spin) and
1.350 (three spins); times 64 (tD/U)^2 that predicts rho_d(0) = 12.74 for
one spin, against the 12.71 the level gives, 0.3% apart at U/Gamma = 10.

The spin dimer on two metals (`runs/dimer.py`, `runs/dimer_linear.py`,
`runs/promote.py`; Js = 0.3, rho J = 0.15, Lambda = 2, 21 Wilson sites
per chain, 44 sites, the tip on spin 0, Lorentzian half-width 0.02). The
TASK routes on one state each, at maxm 128 and nkry 100:

| layout | expansion | time | weight missed (of 0.44) | peak |
|---|---|---|---|---|
| sides | none | 1091 s | 0.417 | 0.018 |
| sides | 8/3 | 2110 s | 7.7e-3 | 0.431 |
| linear, tip first | none | 1052 s | 4e-15 | 0.395 |
| linear, tip first | 8/3 | 1977 s | 4e-15 | 0.436 |

So the linear order holds the vertex with no expansion, as the table of
section 3 said, but holding the start vector is not enough: the
expansion still moves the peak of the linear curve by 10%, through the
Krylov recursion, which leaves the one-site tangent space too. Across
settings the peak of sides with 8/3 is 0.449, 0.432 and 0.429 at maxm
64, 128 (nkry 100) and 128 (nkry 200), and at half-width 0.01 it is
0.790, 0.723 and 0.664, so the sharp line is not converged in nkry. The
energy of the linear order is 1.2e-5 above sides at maxm 128 and 3.5e-3
at maxm 64. The CG correction vector gave peaks of 0.128 and 0.353 at
maxm 128 and 256 and is not a reference (section 3).
The Js sweep, sides with 8/3 against linear tip-first, maxm 128, nkry
200: the two routes agree to 1.4% of the peak up to Js = 0.01 (the
Kondo peak of each spin, its height 0.83 at Js = 0 falling to 0.53 at
Js = 0.01), and differ by 6%, 24%, 23% and 46% of the peak at Js =
0.03, 0.1, 0.2 and 0.5, where the zero-bias value agrees (0.24, 0.05,
0.018, 0.004) and the lines of the singlet-triplet step move. Both are
even in the bias to a few percent (the model is particle-hole symmetric:
up to 0.047 of the peak for sides with 8/3, 0.014 for linear), so the
symmetry does not tell which is right. On a Wilson chain at Lambda = 2
with a Lorentzian of 0.02 the curves are combs of discrete lines, which
hides the steps at eV = +-Js. At Js = 0 linear stalls 1.4e-2 above sides
in energy (nothing couples the second spin to anything but its own
chain, which sits behind the first one), while the spectrum under the
first spin is the same to 0.4%; at Js = 0.003 the two energies are
1.9e-5 apart.
Against an exact reference (sparse ED of 8 sites, 2 spins on chains of
3 Wilson sites, the ground state non-degenerate, each vertex by a Lanczos
continued fraction, `runs/exact_small.py`), with nkry 200 and a
Lorentzian of 0.02, the error over the peak:

| Js | route | maxm 16 | maxm 32 | maxm 64 |
|---|---|---|---|---|
| 0.03 | linear, tip first | 7.3e-4 | 9.0e-6 | 1.4e-7 |
| 0.03 | sides, 8/3 | 3.5e-10 | 3.5e-10 | 3.5e-10 |
| 0.1 | linear, tip first | 9.3e-3 | 1.3e-4 | 1.2e-7 |
| 0.1 | sides, 8/3 | 2.8e-9 | 2.0e-10 | 2.0e-10 |
| 0.3 | linear, tip first | 9.8e-3 | 1.2e-3 | 3.6e-6 |
| 0.3 | sides, 8/3 | 2.0e-8 | 1.5e-10 | 1.5e-10 |

Both converge with the bond dimension, and on a chain this short eight
directions per bond fill most of the space, so the expansion is exact
almost at once. With 4 Wilson sites per chain (10 sites) the ground
state is degenerate between charges, as `anderson_chain` warns for an
even chain, and no single reference state exists. Which route is right
on the long chain, where they differ, is what the variational
references on 44 sites were meant to answer; they failed (section 5,
after the convergence table).

Chains with every spin on its own metal (`runs/chain.py`,
`runs/chain_linear.py`; Js = 0.3, rho J = 0.15, 15 Wilson sites per chain,
maxm 128). With three spins (48 sites) sides with 8/3 and linear
tip-first give the same zero-bias dI/dV under the end spin (0.263 and
0.262) and under the middle one (0.019 and 0.020), and differ by 38% and
47% of the peak in the lines away from zero bias; with four spins (64
sites) the zero-bias values are 0.011 and 0.013. The energy of the linear
order is 7.2e-4 (three chains) and 1.5e-3 (four) above sides at this
bond dimension, so here the linear curves are the less converged ones.
The middle spin, whose exchange with its metal acts as ferromagnetic in
the ground state of the chain, has almost no zero-bias weight, as for one
metal in 4.8i of `docs/development.md`. With sides and comb at maxm 256
the expanded tangent space passed the cap of TASK and the runs stopped.
The spectral moments are a reference that needs no correction vector:
m_k = <GS|X^dag (H-E0)^k X|GS> from the MPS (`runs/moments.py`, X|GS>
and H X|GS> applied at maxm 128), against the same sums over the TASK
poles, on the 44-site dimer at Js = 0.3 (electron part of the first
vertex):

| route | m0 exact / TASK | m1 | m2 |
|---|---|---|---|
| linear, tip first | 0.109856 / 0.109856 | 0.080364 / 0.080364 | 0.069497 / 0.069497 |
| sides, 8/3 | 0.109818 / 0.107849 | 0.080369 / 0.079323 | 0.069511 / 0.068887 |

So the linear order with the tip first holds the first three moments to
all digits, and sides with the expansion misses 1.8%, 1.3% and 0.9% of
them, the weight it leaves out. This is what Fable's decision rests on.
It does not fix the lines: the expansion still moves the peak of the
linear curve by 10% with nothing missing, through higher moments, and
the variational references on Triton (`runs/reference.py`) had not
finished by the end of the night (they did not give a reference, see
below).

Convergence in the bond dimension, the dimer at Js = 0.3, Lorentzian of
0.02 (both bias signs counted in each window):

| route | maxm | dI/dV(0) | weight, \|w\| < 0.15 | 0.15 to 0.35 | 0.35 to 0.6 |
|---|---|---|---|---|---|
| sides, 8/3 | 64 | 0.0094 | 0.0034 | 0.050 | 0.418 |
| sides, 8/3 | 128 | 0.0097 | 0.0040 | 0.058 | 0.354 |
| sides, 8/3 | 256 | 0.0097 | 0.0036 | 0.055 | 0.279 |
| linear, tip first | 64 | 0.0098 | 0.0036 | 0.060 | 0.330 |
| linear, tip first | 128 | 0.0100 | 0.0037 | 0.058 | 0.294 |

So the low-energy part, the zero-bias value and the weight below
|eV| = 0.35, where the singlet-triplet step of Js = 0.3 lies, agrees
between the routes and the bond dimensions to about 10%, while the lines
above 0.35, near the band edge of the Wilson chain, still move by a
third from maxm 64 to 256 (as maxm grows, linear at 128 is closer to
sides at 256 than sides at 128 is). The peak heights and line positions
up there are not converged by either route at these bond dimensions,
and the differences of 6% to 46% between the routes in the Js sweep
are in those lines.

The variational references did not give a reference. Of the 21 tasks of
Triton array 21688905, the 7 chunks of sides at maxm 128 and 2 chunks of
linear at maxm 128 finished, in 17 to 34 h each, and the 12 others hit
the 36 h limit (all of sides at maxm 256). The finished curves are not
a spectral function: dI/dV is negative at 8 of the 21 biases of sides,
-7.4e-5 at zero bias against 0.0097 from both TASK routes, and its
largest value is 0.0042. The two DDMRG calls at each bias (the
electron and the hole part) are each negative, of order 1e-6 to 1e-4,
so the variational functional was not minimized at 44 sites and maxm
128, while on 8 sites it is 5.8e-9 off ED (section 3). As a result
Fable's condition (iii) cannot be tested with these runs, and the
warning on the linear tip-first route stays. What is left as a
reference on 44 sites is the moments above, which both routes satisfy
to the weight they keep. `compare_ref.py` on Triton reads the outputs
in `$WRKDIR/calculations/chains_on_metal/out`.

## 6. A chain between two leads, normal or superconducting

With `site=[0, ns-1]` in the sides order the chain sits between two
semi-infinite leads, and since 2026-11-20 every lead can be its own
metal: `D`, `Lambda` and `Delta` are a number or one value per spin of
the list, `Delta` being the uniform onsite pairing of `kondo_chain` on
every site of that lead. Exchange moves no charge in or out of a lead,
so the phase of a pairing is a gauge and the model has no Josephson
coupling; a complex `Delta` raises.

Checks, all in `tests/test_ysr_chains_on_metal.py`: one paired lead
under one spin is `kondo_chain` term for term, `subgap_energies` gives
the same E_ysr on both to 1e-10, and on v3 and "python" its DMRG E_ysr,
`parity_ground_state` and the TASK STM spectrum agree with ED (1e-7,
1e-7 and 1e-4 of the peak), so this one layout is not refused; at JK = 0 the ED energy is the
Heisenberg chain plus two BdG chains with different Delta, D and Lambda,
to 1e-10, and `fc.lead_qp_gap` is the BdG gap of each; a reordered
`site` list with its values reordered is the same model in all three
orders (spins 0 and 1 of 3, which are not a mirror pair), and the same
list with the values not reordered is not; Delta_R -> -Delta_R leaves
the six lowest levels unchanged to 1e-10; DMRG against ED to 1e-7 on v3
and "python".

The STM spectrum with pairing is not verified on DMRG. On 2 spins with
3 Wilson sites per lead (Delta = 0.15 and 0.3, JK = 1, maxm 256), TASK
with `expand_bonds=8, expand_rounds=3` in the sides order was 3e-4 to
6e-4 of the peak off ED on both backends, against 1e-8 at Delta = 0, and
plain TASK with the tip first in the linear order was 1.6e-3 off on
"python" and 0.89 on v3. The ground state was exact to 8e-13 on both,
with an elastic weight of 1e-25 and no degeneracy, so the cause is in
the TASK space; it was not diagnosed. All three STM entry points
(`ysr.get_stm_spectrum`, `zaverage.stm_poles`,
`surrogate.get_stm_spectrum`) raise on DMRG for these chains, and
`mode="ED", submode="CVM"` runs. On ED a mirror-symmetric junction can
have a degenerate ground level whose members give different spectra (3
spins between two equal leads with Delta = 0.3 at JK = 2: the tips on
the two end spins were 0.13 of the peak apart with one member), and
`get_stm_spectrum` averages over the level, after which they agree to
below 1e-8 of the peak.

The example `examples/kondo/spin_chain_between_metal_and_superconductor`
puts a normal metal on the left and a superconductor on the right of 4
spins (8 Wilson sites per lead, Lambda = 4, JK = 0.5, Js = 0.3). The
binding <S_end . s_0> of the right spin to its lead goes from -0.1227 at
Delta = 0 to -0.0338 at Delta = 0.5, and that of the left spin from
-0.1227 to -0.1217; maxm 64 and 128 agree to five digits.

## 7. What is not in it

- A common metal: the baths of different spins are independent.
- Energies resolved by the parity of each lead, and a charge path
  through the chain (the Josephson effect).
- The STM spectrum with superconducting leads on DMRG.
- T > 0 with several baths: TASK is a T = 0 method, and the NRG tail
  takes one Wilson chain.
- The NRG tail with several baths, which is the interleaved chains of
  milestone N6 of `docs/nrg_mps_kondo_plan.md`.
- Charge fluctuations on the spins other than the Anderson levels.
