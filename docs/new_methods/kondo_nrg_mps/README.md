# The STM spectrum of several spins by NRG, with a DMRG head

This folder documents the numerical renormalization group (NRG) route to
the STM spectrum of a spin cluster on a metal: `nrgspectral.py` (the
spectral NRG with the complete basis, and the exact start shell of the
spins), `nrgmps.py` (the shells read off a DMRG state, and the hand-over
from a DMRG head to an NRG tail) and `ysr.get_stm_spectrum(fc, es,
submode="NRG")`, the public entry. The plan it follows is
`docs/nrg_mps_kondo_plan.md`, milestones N0 to N4 and N7. The tests are
`tests/test_nrgspectral.py` and `tests/test_nrgmps.py`, the examples
`examples/kondo/kondo_resonance_nrg` and
`examples/kondo/spin_chain_nrg_head_tail`.

## 1. The question

The STM spectrum of one Kondo spin is a solved problem for NRG. Can we
compute it for a chain of spins on a metal, where exact NRG has to keep
every state of the spins while the energy of a shell is above their
bandwidth, and where TASK (the tangent-space Krylov route of
`get_stm_spectrum`) does not resolve the Kondo resonance below about
T_K/3?

The answer is yes up to at least 6 spins, where the exact route is enough:
the spins and the first Wilson site are diagonalized exactly and the NRG
goes on from there. Beyond that, DMRG solves the spins and the first
Wilson sites (the head), and the NRG goes on from the bond basis of that
state (the tail). The tail gives the spectrum below the energy of the
hand-over, and TASK on the head gives it above. This second route agrees
with the exact one on 3 spins; at 7 and 8 spins it runs and gives the
expected physics, but it was not converged in maxm (section 6).

## 2. The spectral NRG

The model is the one of `ysr.spin_chain_on_metal`, ns spins 1/2 with a
Heisenberg exchange Js, one of them Kondo coupled by JK to the first site
of a Wilson chain, and for ns = 1 the one of `ysr.kondo_chain` at
Delta = 0. The quantity is

    A(w) = sum_sigma sum_ab (p_a + p_b) |<a|X^dag_sigma|b>|^2 delta(w - E_a + E_b),

with X^dag_sigma the composite tunneling vertex of
`ysr.tunneling_operators`. At T = 0 this is dI/dV with a normal tip, in
units of 2 pi T0^2 per unit of tip density of states, the convention of
`get_stm_spectrum`.

The iteration is the one of `nrgthermo`: shells labelled by the bath
charge Q and 2S^z, one Wilson site added at a time with the
Jordan-Wigner sign (-1)^Q, the lowest `keep` states kept and the cut moved
past degenerate partners. Two things are new. Every shell keeps its
discarded states, so that the discarded states of all shells and all
states of the last one form a complete basis (Anders and Schiller,
cond-mat/0604517). And the density matrix is the full one of Weichselbaum
and von Delft (cond-mat/0607497): every discarded state of shell n
carries its Boltzmann weight with the degeneracy 4^(N-n) of the sites
after it, and a backward pass gives the reduced density matrix R^n of
every shell. The Lehmann sum then runs over the pairs of states of each
shell that are not both kept, with the weight

    X^dag_ab [ (X^dag R^n)_ab + (R^n X^dag)_ab ],

the electron and the hole part of the anticommutator (Peters, Pruschke and
Anders, cond-mat/0607494). The weights add up to
sum_sigma Tr[rho {X_sigma, X^dag_sigma}] exactly, at any `keep`. At
T = 0 rho is the equal mixture of the ground manifold of the last shell,
the states within `DEG_TOL` = 1e-4 of omega_N. The tolerance is not
1e-8 because on an even chain the fixed point has a fourfold ground
manifold (a zero mode), split by an irrelevant correction of 1e-7 omega_N
at 30 sites and 2.4e-8 at 32.

The poles are binned on a logarithmic mesh (200 bins per decade) and
broadened with the log-Gaussian kernel of `ysrtk.zaverage.broaden_poles`,
z-averaged over twisted meshes (`nrgspectral.stm_spectrum`). Two rules
are needed on top of the kernel.

- At T = 0 the curve below the last shell is that of the kernel, not of
  the spectrum: below omega_N it falls to 0.86 of the plateau. Since the
  stop rule ends the chain past the fixed point, where the spectrum no
  longer changes, `nrgspectral.curve` reads every |w| below
  `FLOOR` = 30 times omega_N at that energy. A twist average has to read
  all meshes at one energy, 30 times the largest omega_N: read at their
  own omega_N, the meshes are all at the same place of their comb and the
  average gives 1.93 instead of 1.00.
- At T > 0 the NRG keeps the weight inside |w| < T but not its shape. On
  the free chain at T = 3.5e-4 the poles carry 250 times too much weight
  below 1e-7, while the weight inside |w| < T is right to 0.15%. The
  kernel with omega0 = T then gives a dip to 0.68 of a flat density inside
  |w| < T and a bump just above it, for the exact density of states as
  well. What STM measures at T > 0 is dI/dV = int dw A(w) [-f'(w - eV)],
  with the Fermi function of the tip, and this average over the thermal
  window is what `nrgspectral.thermal_didv` returns (omega0 = T/4). On the
  free chain it is the exact density of states of the same chains to
  0.3% at every bias.

The tail stops at the fixed point (`stop_tol`): once the lowest levels of
two consecutive shells, in units of their scale, agree with those of the
shells two sites before, `stop_extra` more shells are added. At T > 0 it
stops at the first shell below T/100.

## 3. Several spins in the start shell

`nrgspectral.spin_cluster_shell` diagonalizes the spins and the first
Wilson site exactly, in blocks of (Q, 2S^z), up to 10 spins. The keep it
needs was measured at Js = 1e-3 D, rho JK = 0.15, Lambda = 4, one mesh,
as the largest change of A(w) on 1e-9 to 0.3 against keep = 3200, in
units of the peak of A. The peak is the unit because an even chain is a
singlet whose A(w) vanishes at low w, where a relative change is a ratio
of two near-zero numbers (a first version of this table, relative point
by point, gave 100% at 4 spins and 5e4% at 6 that way):

| ns | keep 400 | keep 800 | keep 1600 | time at 1600 |
|---|---|---|---|---|
| 1 | 0.20% | 0.68% | 0.02% | 45 s |
| 2 | 0.75% | 0.16% | 0.11% | 37 s |
| 3 | 0.88% | 0.86% | 0.09% | 49 s |
| 4 | 1.0% | 0.14% | 0.16% | 38 s |
| 5 | 17% | 0.98% | 0.15% | 56 s |
| 6 | 106% | 17% | 1.5% | 46 s |

keep = 3200 itself takes 230 to 430 s per mesh. At Js = 1e-2 D three spins
behave the same (0.75% at keep = 400). So up to 4 spins keep = 800 is
enough, the keep needed grows from 5 spins on, and at 6 spins keep = 1600
is still 1.5% of the peak away (5.6% relative where A is above 1% of its
peak). Where the exact route stops is between 6 spins and the next one we
did not measure; the head of section 4 is the route beyond. Below the
keep that converges, the result can be qualitatively wrong, not just less
accurate (106% of the peak at 6 spins and keep = 400). keep = 1600 is
the default of the examples with several spins.

Beyond 6 spins the reference is keep = 6400, on one core of Triton (Xeon
Gold 6248), 12 to 14 min and about 30 GB per run. At Js = 1e-3 D, as in
the table, keep = 3200 is 0.16% of the peak off at 6 spins, 1.0% at 7 and
16% at 8, and keep = 1600 is 1.4%, 17% and 111% off. At Js = 1e-2 D on 7
spins keep = 3200 is 0.19% off with the metal under spin 3 and 0.31%
under the end spin. So the exact route holds at 7 spins with keep = 3200,
and at 8 spins it does not hold even at keep = 6400.

The shortcut of the plan, the lowest states of the bare spin chain as the
impurity, is uncontrolled, as the plan expected. With 3 spins at
Js = 1e-3 D and the metal on the end spin at rho JK = 0.1, all 8 states
give a resonance whose half width is 1.3e-5, and the projected ground
doublet gives 3.2e-7, 42 times smaller. The states that are cut take part
in the Kondo flow between D and Js. `spin_cluster_shell(project=k)` keeps
the shortcut for this test only.

## 4. The shells of a DMRG state

A DMRG ground state of the chain is an NRG in disguise (Verstraete et
al., cond-mat/0504305): the left isometries up to the Wilson site k map
the block onto chi states, which play the part of the kept states of
shell k (Saberi, Weichselbaum and von Delft, arXiv:0804.0193). `nrgmps`
writes any operator of the block in that basis by contracting the MPO the
backend exports (the arrays TASK uses) with the left isometries, so the
site conventions and the Jordan-Wigner strings are those of the backend.
Three things are needed for this to be a shell.

- A graded basis. N and S^z of the block are diagonalized together in the
  left basis, and their eigenvalues have to be integers to 1e-8. A plain
  solve of a doublet is an arbitrary member of it, so `head_state` solves
  in the sector (Nf, 2S^z), with the warm start that keeps the spin sites
  singly occupied on v3, and promotes the state to dense indices. In a
  field it scans the sectors.
- A real basis. The "python" backend returns a state that is real up to a
  phase in a complex gauge; `real_gauge` takes a real orthonormal basis of
  the span of each left basis and checks that the span is closed under
  conjugation (residual at roundoff).
- The block Hamiltonian, the sum of the terms of H with every factor in
  the block.

The gate of the plan was whether this basis holds the low excited states
of the block. On v3, 20 Wilson sites, rho JK = 0.15, the levels of the
block Hamiltonian in the left basis against the NRG shells at keep = 1600:

- One spin, maxm = 64: the lowest 10 levels to 1e-4 to 1e-3 of omega_k,
  about 20 of the 64 states within 1e-2 omega_k. maxm = 128 and 256: the
  bond stops at 108 to 142 states (the cutoff 1e-12, not maxm), 63 to 88
  levels within 1e-2 omega_k up to the site 10, the lowest 10 to 3e-5.
- Three spins (Js = 1e-2 D), maxm = 64: the lowest 10 levels off by 0.02
  to 0.5 omega_k. maxm = 256: 72, 27, 88 and 35 levels within 1e-2 omega_k
  at the sites 2, 6, 10 and 14, the lowest 10 to 2e-4.
- Where the left basis is complete (the first Wilson sites of a short
  chain) the levels are those of NRG to 4e-15, on both backends.

So the basis holds a few tens of low states well, and the number grows
with maxm until the cutoff of the Schmidt weights stops it.

## 5. The hand-over and the tail

`handover_shell` turns the left basis at the Wilson site k0 into the start
shell of the NRG: the block Hamiltonian, the labels, f of the site k0 for
the next hopping, the tip vertices, and their anticommutator projected as
one operator, for the sum rule. The tail is then the iteration of section
2 on the Wilson sites after k0 of the same mesh.

With a complete basis at the hand-over (`pad_bonds="full"`) the head and
the tail are the exact NRG of the same chain to 5e-14, at T = 0 and
T > 0, on "python" and on 3, with the bath on the end spin and on the
first spin of two. With the basis DMRG gives, against the exact NRG of
the same chain on one mesh below omega_k0/2 (head of 12 Wilson sites,
the hand-over m sites before its end, keep = 800 in the tail and 1600 in
the reference):

| spins | m | maxm | chi | largest error | mean error |
|---|---|---|---|---|---|
| 1 | 4 | 128 | 121 | 1.5% | 0.1% |
| 3 | 3 | 128 | 64 | 4.6% | -3.1% |
| 3 | 4 | 256 | 146 | 1.8% | -0.1% |
| 3 | 5 (head 14) | 512 | 224 | 1.9% | -0.1% |

The largest error of 1.5 to 1.9% does not move with chi above about 120.
It is mostly the different keep of the two runs: at 3 spins keep = 800 is
1.0% off keep = 3200 (section 3), and two chain lengths differ by up to
0.7% at a fixed energy on one mesh. With keep = 1600 on both sides the
head and the tail agree with the exact route to 8e-4 below 1e-3 (the
example below). A hand-over basis of 64 states is not enough. So
At T > 0 the same head (maxm = 160, 146 states at the hand-over, keep =
1600 on both sides) gives dI/dV with the Fermi function of the tip at
T = 1e-4, about omega_k0/100, to 0.93% of the exact route below
omega_k0/2, 0.02% on average.

A head with an odd number of electrons is one member of a doublet, and
its left basis is not symmetric under the flip of all spins: on 3 spins
12 of the partner blocks (Q, 2Sz) and (Q, -2Sz) of the hand-over have
different sizes, while the blocks of equal size agree exactly. The
average of partner spectra at B = 0 (`flip`) restores most of the
symmetry, and it is needed: with the metal under the middle spin, an
effectively ferromagnetic and so marginal flow, the tail is 3.8% off the
exact route with it (mean -2.3%) and 59% without it; under the end spin
both are below 0.3% apart. A head with an even number of electrons
(2Sz = 0) is symmetric to 4e-8.

The hand-over basis of a head with a moment is now built from both
members of the doublet (`symmetric="auto"`, the default): the left basis
of the state joined with that of R|GS>, R = exp(-i pi S^y_tot) applied
site by site, which is the other member exactly (same energy, opposite
S^z). The join is closed under R exactly: the basis of (Q, -2Sz) is the
image of that of (Q, 2Sz), and the one of 2Sz = 0 is cut by the
eigenvalues of an overlap matrix that commutes with R. This matters. A
first version that dropped the nearly dependent directions (within 1e-3
of the span) label by label did not keep the closure, and the tail
without the average was still 60% off; we read this as the split of the
partner levels acting on the free moment as a field, which the flow
makes relevant (the size of the split was not measured). Measured on
the same 3 spins (head of 12 Wilson sites, maxm = 160, keep = 1600 in the
tail and in the exact route, one mesh, below w_join/3):

| metal under | basis | flip average | max error / peak | max relative | mean relative |
|---|---|---|---|---|---|
| middle spin | plain (145 states) | yes | 3.4% | 12.6% | -2.9% |
| middle spin | plain | no | 3.1% | 60% | -23% |
| middle spin | joined (284 states) | yes or no | 1.5% | 5.6% | -0.4% |
| end spin | plain or joined | yes or no | 1.3 to 1.4% | 1.7 to 1.9% | -0.1% |

The exact route gives the same curve with and without the average (to
1e-4 of the peak), so the average was only ever standing in for the
basis. The joined basis doubles the hand-over basis and takes the tail
from about 48 s to 70 to 90 s; the DMRG head and TASK are unchanged.

m = 4 is the default, and the head needs a maxm of at least about 160:
the bond at the hand-over stops near 150 states at the cutoff 1e-12 even
with maxm = 256 or 512. The example `spin_chain_nrg_head_tail` (3 spins,
maxm = 160, nkry = 30, one mesh) agrees with the exact route to 8e-4
with one member of the doublet in the hand-over basis and to 2.9e-4 with
both (the default since 2026-11-14), both
below 1e-3, in 6 minutes on one core.

## 6. The joined spectrum

`head_tail_poles` keeps the tail poles below w_join (default omega_k0, the
scale of the hand-over shell) and adds the TASK poles of the head state
above it (`head_task_poles`, one block Krylov space from P X^dag_sigma|GS>
and P X_sigma|GS> of both spins). The tail cannot see the part of
X^dag|GS> outside the hand-over basis, at high energy (`missing_weight`):
1.4% of the weight on a head of 7 Wilson sites, 80% on a head of 10 with
59 states at the hand-over. TASK supplies it above w_join. Nothing
conserves the weight across the cut: the joined poles carry 0.6120 of a
sum rule of 0.6138 on the head state (one spin, rho J = 0.2, head of 10
sites, m = 3, maxm = 128), and 1.1% less than it on the head of 7 sites
at nkry = 80.

The whole curve was compared with the exact route on 3 spins (Js =
1e-2 D, rho JK = 0.15, the metal under the end spin), averaged over two
meshes with b = 0.6: a head of 12 Wilson sites at maxm = 160 with the
hand-over 4 sites before its end (146 and 159 states, w_join = 1e-2 and
2e-2), keep = 800 in the tail and 1600 in the exact route. The largest
relative error and its mean, on three ranges of |eV|:

| nkry | below w_join/2 | w_join/2 to 2 w_join | above 2 w_join | time, 2 meshes |
|---|---|---|---|---|
| 30 | 1.5%, 0.1% | 8%, -2% | 12%, -2.5% | 476 s |
| 90 | 1.1%, 0.2% | 5%, -3% | 3.8%, -0.5% | 1268 s |

Below the cut the error is that of the tail, here mostly the keep of 800
against 1600 (section 5). Above it the error is the Krylov resolution of
TASK, which falls with nkry; 90 is the default. Around the cut the two
routes put their weight on different lines, and the joined curve is 3%
low on average there. The joined poles carry the sum rule of the head
state to +0.3% and -0.2% at nkry = 90 (+0.5% and -1.2% at 30).

Step 9 of the plan's ladder takes the route past 3 spins. It was run on 7
and 8 spins (Js = 1e-2 D, rho JK = 0.15, one mesh, b = 0.8, the tail
alone, keep = 1600), with the hand-over at the Wilson site 7 (scale
1e-2) and m = 4, 5 and 6 head sites after it, maxm 256 to 768, so the
hand-over basis grows from about 160 to 420 states. The size of that
basis is set by m more than by maxm. Between the two largest bases A(w)
changes by 0.00% of its peak on 7 spins with the metal under the end
spin, 0.37% under spin 3, and by 0.07% and 0.12% on 8 spins under spin 7
and spin 4. The bases of m = 4 (160 to 210 states) are 2 to 20% off the
largest one, so the about 150 states of section 5 are not enough on
these chains.

On 7 spins the exact route holds (section 3), and head and tail agree
with it at keep = 6400. With the metal under the end spin (a moment that
the metal screens, a Kondo peak of 4.2) the difference is 1.6% of the peak
below 1e-5, 0.4% from 1e-5 to 1e-4 and at most 0.3% up to 3e-3. Under
spin 3 (a triplet whose exchange with the metal acts as ferromagnetic, no
Kondo peak) it is 0.07% below 1e-3 and 1.3% from 1e-3 to 3e-3. The 1.6%
does not move with the hand-over basis (1.67% at 168 states, 1.56% at
422), so it is probably the keep of the tail, 1600 against 6400, which
was not measured. From 3e-3 up, near and above w_join, the tail alone is
11% and 62% off, and there TASK on the head takes over. On 8 spins the
exact route does not hold, so there only the convergence of the head
route is measured; the chain is a singlet and A vanishes below 1e-4, as
it should.

The same ladder on 10 spins, with the metal under the end spin and under
its neighbour, converges more slowly with the neighbour: between the two
largest hand-over bases (402 and 414 states, 327 and 333) A changes by
0.47% and 1.35% of its peak, and the next ones (340, 292) are 0.78% and
2.3% off. Both chains are singlets, and A vanishes at low energy (below
1e-10 of its peak at 1e-9 and about 1e-5 at 1e-5), so these are runs
without a Kondo resonance, and nothing beyond 7 spins checks the
resonance itself. Each run of the head took 1.5 to 4.5 min on one core.

## 7. Checks

- Untruncated, the poles are those of a dense diagonalization of the same
  model, one spin and several, with a field, U and a twist, at T = 0 and
  T > 0 (1e-13 on the resolvent).
- They are the ED STM spectrum of `kondo_chain` and `spin_chain_on_metal`
  (`get_stm_spectrum(mode="ED", submode="ED")`), whose spins are
  penalized fermion sites (2e-14).
- With three coupled Wilson sites and three decoupled ones, the truncated
  result is the untruncated one of the short chain to 1e-15, with up to
  46% of rho on the discarded states of the earlier shells.
- The sum rule and Tr rho = 1 hold to 1e-12 at small keep.
- At J = 0 the spectrum is (U^2/2 + 3/8) times the density of states per
  spin of the first Wilson site.
- One spin at rho J = 0.2, four meshes, b = 0.4: the T-matrix
  t(w) = pi^2 rho J^2 A(w)/2 is 1 to 1.5% from w = 0 to 1e-4 (Wilson's
  T_K = 3.49e-3), and dI/dV at zero bias gives 0.989 at T = 0.1 T_K and
  0.720 at T = T_K.
- The same spin in a field B on the impurity: at B = T_K the top stays
  flat at 0.90, and at B = 10 T_K the resonance splits into two peaks at
  +-9.06 T_K with 0.21 at zero bias, as in Costi's NRG.
- An independent review recomputed the untruncated and decoupled-tail
  cases with its own dense diagonalization in a different Jordan-Wigner
  order (27 cases, 3e-11 at most) and found six defects, all fixed and
  pinned by tests: the real gauge of the "python" state, the Zeeman sector
  of the head, the broken end of a 120-site Wilson chain, the absolute
  energies in the Boltzmann weights at T below 1e-12, a field below 1e-10
  averaged away, and the definition of the missing weight.

## 8. Cost

One core, BLAS on one thread. The exact route costs 40 to 56 s per mesh
at keep = 1600 on a Wilson chain of 44 sites, for one to six spins, and
230 to 430 s at keep = 3200; its cost is set by keep, not by the number of
spins, up to the dense first shell of 4 2^ns states. On the head route the
NRG tail is cheap (9 s at keep = 800 from a hand-over of 126 states), and
the time goes to the head: on 3 spins and 12 Wilson sites at maxm = 128
the sector solve took 56 s and TASK 115 s for 30 blocks of four, twice
that for 60, and at maxm = 160 the full joined spectrum took about 4 min
per mesh at nkry = 30 and 10 min at nkry = 90. TASK grows as chi^3,
which is why `head_task_poles` builds one block Krylov space for the four
starts (three times faster than four separate recursions on a short
chain). A coarse Krylov space is not enough above w_join: nkry = 30 was
12% off there, nkry = 90 3.8% (section 6).

## 9. What the route does not contain

- Sharp features at finite energy: the logarithmic mesh broadens a line
  at w by a fraction of w, so the spin excitation steps near eV = Js are
  those of TASK on the head.
- More than one channel: the metal is under one spin.
- A temperature above the scale of the hand-over: TASK on the head is a
  T = 0 method, so T has to stay far below w_join (a warning above
  w_join/10). Above it, the exact route of section 3 is the only one.
- A superconducting substrate: the (Q, 2S^z) labels need a conserved
  charge.
- A reference beyond 7 spins: at 8 spins the exact route is not converged
  at keep = 6400, and the head route is checked by its convergence only
  (section 6).
- A hand-over below the floor of the MPO. ITensor's svdMPO (v3) cuts a
  squared singular value below an absolute 1e-13 at unit scale, so a
  Wilson hopping below about 3e-7 of the unit of H leaves the MPO of the
  head, and the sites after it are nearly decoupled from the chain (on 30
  Wilson sites at Lambda = 4 the MPO link dimension is 7 down to
  |t| = 6e-7, 6 at 3.0e-7 and 3 from 1.5e-7 on). So the depth floor of
  the plan's N2 is a scale, not a number of sites. Measured on one spin
  (maxm 128 to 512, cutoff 1e-12 and 1e-14, 30 or 90 sweeps) and three
  spins (maxm 256 and 512), the lowest levels of the block agree with
  NRG to below 8e-4 omega_k down to the deepest shell whose hopping is
  at least 3.5e-7 (the Wilson site 21 at Lambda = 4, 32 at 2.5 and 17 at
  6), and they are off by 0.36 to 1.8 omega_k from the first shell whose
  hopping is 3.0e-7 or less, whatever maxm, cutoff and sweeps. On
  "python" the builder truncates its return sweep at a relative 1e-14 and
  drops channels one bond earlier at Lambda = 4, and the solve of a head
  of 30 Wilson sites failed the real gauge (residual 1.07e-4 against
  1e-4). `nrgmps.lost_wilson_bonds` reads the dropped bonds off the MPO,
  and `head_tail_poles` warns for a hand-over whose block holds one (at
  Lambda = 6 this warns one shell before the measured floor). The
  hand-overs of this record are at the Wilson site 7, scale 1e-2, far
  above it. A per-bond cutoff of the MPO would remove the floor; it is
  the bond-local truncation that the 2026-09-25 audit left open.
- `head_tail_poles` also warns when the hand-over basis holds fewer than
  `MIN_HANDOVER` = 100 states.

## 10. Papers

- Verstraete, Weichselbaum, Schollwoeck, Cirac, von Delft,
  cond-mat/0504305; Saberi, Weichselbaum, von Delft, arXiv:0804.0193:
  NRG states are MPS, and the DMRG state of a Wilson chain holds them.
- Anders and Schiller, cond-mat/0604517: the complete basis.
- Weichselbaum and von Delft, cond-mat/0607497; Peters, Pruschke, Anders,
  cond-mat/0607494: the full density matrix and the spectral sum.
- Costi, cond-mat/0004302: the T-matrix of the Kondo model.
- Bulla, Costi, Pruschke, cond-mat/0701105: the NRG review.
