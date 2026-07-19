# DMRGPY Physics User Guide

This guide explains what DMRGPY computes from a **physics** point of
view — the quantities each method returns, the formula behind them, and
when to reach for which method. It intentionally says nothing about the
solver backends (DMRG vs ED, C++ vs Python vs Julia) or code layout —
see `docs/documentation.md` for that. Every `Many_Body_Chain` subclass
(`Spin_Chain`, `Fermionic_Chain`, `Bosonic_Chain`, ...) exposes the same
physics-facing API described here, `mode="DMRG"|"ED"` selectable on
almost every call so a result can be cross-checked against an exact
reference on small systems.

## Contents

1. [Physical models and Hilbert spaces](#1-physical-models-and-hilbert-spaces)
2. [Building a Hamiltonian and observables](#2-building-a-hamiltonian-and-observables)
3. [Ground-state properties](#3-ground-state-properties)
4. [Excited states and energy gaps](#4-excited-states-and-energy-gaps)
5. [Entanglement and quantum information](#5-entanglement-and-quantum-information)
6. [Dynamical (frequency-dependent) correlators](#6-dynamical-frequency-dependent-correlators)
7. [Real-time dynamics: quenches](#7-real-time-dynamics-quenches)
8. [Density of states](#8-density-of-states)
9. [Finite temperature](#9-finite-temperature)
10. [Topological invariants](#10-topological-invariants)
11. [Mean-field decoupling](#11-mean-field-decoupling)
12. [Fidelity susceptibility and quantum phase transitions](#12-fidelity-susceptibility-and-quantum-phase-transitions)
13. [Ground-state degeneracy](#13-ground-state-degeneracy)
14. [Reduced density matrices and operator distributions](#14-reduced-density-matrices-and-operator-distributions)
15. [Post-processing tools](#15-post-processing-tools)
16. [Worked-example cookbook](#16-worked-example-cookbook)

## 1. Physical models and Hilbert spaces

Every model is a chain of $n$ local Hilbert spaces $\mathcal H=\bigotimes_{i=1}^n \mathcal H_i$, on which a Hamiltonian and observables are built out of local operators.

| Chain class | Local Hilbert space | Key operators |
|---|---|---|
| `spinchain.Spin_Chain` | spin-$S$, $S\in\{\tfrac12,1,\tfrac32,2,\tfrac52,3\}$ per site | $S^x_i,S^y_i,S^z_i$ |
| `fermionchain.Fermionic_Chain` | spinless fermion (occupied/empty) | $c_i,c_i^\dagger,n_i=c_i^\dagger c_i$, Jordan-Wigner string $F_i$ |
| `fermionchain.Majorana_Chain` | Majorana fermion | Majorana operators built from `Fermionic_Chain` |
| `fermionchain.Spinful_Fermionic_Chain` | spin-$\tfrac12$ fermion (4 states: $0,\uparrow,\downarrow,\uparrow\downarrow$), built from two interleaved spinless sites per physical site | $c_{i\sigma},c^\dagger_{i\sigma},n_{i\sigma}$, plus derived $S^x_i,S^y_i,S^z_i=\tfrac12(n_{i\uparrow}-n_{i\downarrow})$, onsite pairing $\Delta_i=\tfrac12 c_{i\uparrow}c_{i\downarrow}$ |
| `bosonchain.Bosonic_Chain` | truncated boson Fock space, $n_i\in\{0,\ldots,n_{\max}\}$ (default $n_{\max}=4$) | $a_i,a_i^\dagger,n_i$, occupation projectors $\hat n_i^{(k)}=\lvert k\rangle\langle k\rvert$ |
| `parafermionchain.Parafermionic_Chain` | $\mathbb Z_N$ parafermion (clock model), $N\in\{2,3,4\}$ | clock/shift operators $\sigma_i,\tau_i$ and composite parafermion operators $\chi_i,\psi_i$ built as $\tau$-string $\times\sigma_i$ |

Spinful fermionic chains are built by *interleaving* two spinless
fermionic sites per physical site (site $2i$ = spin up, site $2i+1$ =
spin down) rather than by a genuinely 4-dimensional local space, so that
the same Jordan-Wigner machinery used for spinless fermions applies
unchanged; `Spinful_Fermionic_Chain` wraps this bookkeeping for you.

## 2. Building a Hamiltonian and observables

Any Hamiltonian or observable is written exactly as its second-quantized
or spin-operator expression, using operator lists indexed by site, e.g.
for the spin-$\tfrac12$ Heisenberg chain

$$H=\sum_{i=1}^{n-1} J\,\mathbf S_i\cdot\mathbf S_{i+1}=\sum_{i=1}^{n-1} J\left(S_i^xS_{i+1}^x+S_i^yS_{i+1}^y+S_i^zS_{i+1}^z\right)$$

```python
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(30)]
sc = spinchain.Spin_Chain(spins)
h = 0
for i in range(29):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
```

Any sum of products of these operators is a valid `MultiOperator`, so
essentially any local (or long-range) Hamiltonian on the given lattice
can be built term by term: single-ion anisotropy $D(S_i^z)^2$,
Dzyaloshinskii-Moriya terms $\mathbf D\cdot(\mathbf S_i\times\mathbf
S_{i+1})$, biquadratic exchange $(\mathbf S_i\cdot\mathbf S_{i+1})^2$,
hopping $t\sum_i(c_i^\dagger c_{i+1}+\text{h.c.})$, Hubbard interaction
$U\sum_i n_{i\uparrow}n_{i\downarrow}$, and so on — see §16 for concrete
Hamiltonians.

## 3. Ground-state properties

**Ground-state energy**

$$E_0=\langle\mathrm{GS}|H|\mathrm{GS}\rangle=\min_{|\psi\rangle}\frac{\langle\psi|H|\psi\rangle}{\langle\psi|\psi\rangle}$$

```python
e0 = sc.gs_energy()          # E0
wf = sc.get_gs()             # the |GS> wavefunction itself
```

**Expectation values and moments.** For any operator $O$ built the same
way as $H$,

$$\langle O\rangle=\langle\mathrm{GS}|O|\mathrm{GS}\rangle,\qquad \langle O^n\rangle=\langle\mathrm{GS}|O^n|\mathrm{GS}\rangle$$

```python
mz = [sc.vev(sc.Sz[i]).real for i in range(n)]        # local magnetization profile
e2 = sc.vev(h, npow=2)                                  # <H^2>, e.g. for fluctuations
```

**Energy fluctuation** (a measure of how sharply the DMRG/ED state is an eigenstate, and physically the variance of $H$ in the prepared state):

$$\delta E=\sqrt{\langle H^2\rangle-\langle H\rangle^2}$$

```python
de = sc.gs_energy_fluctuation()
```

**Static two-point correlators.** `sc.vev(sc.Sz[0]*sc.Sz[i])` gives
$\langle S^z_0 S^z_i\rangle$ directly; `correlator.get_correlator`
provides shorthand names for common correlators over a list of site
pairs, e.g. `"SS"` for the full dot product

$$\langle\mathbf S_i\cdot\mathbf S_j\rangle=\langle S_i^xS_j^x\rangle+\langle S_i^yS_j^y\rangle+\langle S_i^zS_j^z\rangle$$

or fermionic correlators like `"cdc"` ($\langle c_i^\dagger c_j\rangle$), `"density"`/`"densitydensity"` ($\langle n_in_j\rangle$), and pairing correlators `"delta"`/`"deltadeltad"`. These static correlators are the equal-time limit of the dynamical correlators in §6, and are what you Fourier-transform to get an equal-time structure factor $S(q)=\sum_{i,j}e^{iq(i-j)}\langle O_iO_j\rangle$.

## 4. Excited states and energy gaps

**Low-lying spectrum.**

```python
es = sc.get_excited(n=4)          # [E0, E1, E2, E3]
gap = sc.get_gap()                # E1 - E0
```

`get_excited_states(n, purify=True)` additionally returns the
wavefunctions; `purify=True` re-diagonalizes $H$ in the Gram-Schmidt
orthogonalized subspace spanned by the raw excited-state MPS, correcting
for near-degenerate states that DMRG's iterative solvers can otherwise
mix.

**Sector (charge) gaps.** For a conserved quantity $A$ with $[H,A]=0$
(e.g. total particle number $\hat N=\sum_i n_i$, or total $S^z$), the gap
to the lowest state with $\langle A\rangle$ shifted by $d$ from the
ground-state sector is obtained by adding a Lagrange-multiplier penalty
and increasing $\lambda$ until the constraint is satisfied:

$$H_\lambda=H+\lambda\big(A-\langle A\rangle_0-d\big)^2,\qquad \Delta_A(d)=E_\lambda-E_0$$

```python
gap_charge = fc.get_charge_gap(d=2)   # gap to add/remove a pair of particles
```

This is exactly how a superconducting/charging gap is extracted in a
finite fermionic chain: $d=2$ probes the energy cost of adding a pair of
particles (Cooper-pair-like excitation) rather than a single particle,
which avoids odd/even-parity finite-size artifacts.

**Example: the Haldane gap.** For the spin-1 Heisenberg chain, the
famous Haldane gap between the (nearly four-fold degenerate, due to
fractionalized $S=\tfrac12$ edge states) ground-state manifold and the
first bulk excitation is

```python
es = sc.get_excited(n=6)
haldane_gap = es[4] - es[0]     # skip the 4 near-degenerate edge states
```

## 5. Entanglement and quantum information

**Entanglement entropy of a real-space bipartition.** Cutting the chain
at bond $(i,i+1)$, the Schmidt decomposition of the ground state is
$|\mathrm{GS}\rangle=\sum_\alpha\lambda_\alpha|\alpha\rangle_L\otimes|\alpha\rangle_R$
and the von Neumann entanglement entropy of either half is

$$S_{i}=-\sum_\alpha\lambda_\alpha^2\log\lambda_\alpha^2=-\mathrm{Tr}\big[\rho_L\log\rho_L\big],\qquad \rho_L=\mathrm{Tr}_R|\mathrm{GS}\rangle\langle\mathrm{GS}|$$

```python
s = wf.get_bond_entropy(i, i+1)
```

`get_site_entropy(i)`/`get_pair_entropy(i,j)` instead build a reduced
density matrix from local *projectors* (e.g. $S^z=\pm\tfrac12$ for
spin-$\tfrac12$, occupation number for fermions) rather than a full MPS
bond cut, which lets you ask about the entanglement of a single site or
a pair of (possibly non-adjacent) sites with the rest of the system.

**Mutual information** between two sites (or subsystems) $i,j$:

$$I(i,j)=S_i+S_j-S_{ij}$$

```python
mi = wf.get_mutual_information(1, 2)
```

quantifies total (classical + quantum) correlation between $i$ and $j$,
and is often used to map out effective couplings or emergent
degrees of freedom (e.g. between the two edge spins of a Haldane chain).

**CFT central charge.** For a critical (gapless) chain, the
entanglement entropy of a cut at position $\ell$ in an open chain of
length $L$ follows the Calabrese-Cardy formula:

$$S(\ell)=\frac{c}{6}\log\!\left[\frac{2L}{\pi a}\sin\frac{\pi\ell}{L}\right]+\text{const.}$$

DMRGPY computes $S(\ell)$ at every bond and least-squares fits this
formula to extract the central charge $c$ — the universal number that
identifies the underlying conformal field theory (e.g. $c=\tfrac12$ for
the critical Ising chain, $c=1$ for a free boson / XX chain):

```python
c = wf.get_CFT_central_charge()
```

**Single-particle correlation matrix and orbital entanglement.** For
fermionic systems, the correlation matrix

$$C_{ij}=\langle c_i^\dagger c_j\rangle$$

is directly accessible (`sc.get_correlation_matrix()`), and diagonalizing
it, $C=U\,n\,U^\dagger$, gives natural-orbital occupations
$n_\alpha\in[0,1]$ (`get_correlated_orbitals`) and an orbital-resolved
entanglement entropy

$$S=-\sum_\alpha\Big[n_\alpha\log n_\alpha+(1-n_\alpha)\log(1-n_\alpha)\Big]$$

(`get_correlation_entropy`). `get_highorder_correlation_matrix` and
`get_four_correlation_tensor` give the corresponding two-particle
correlators $\langle c_i^\dagger c_j^\dagger c_l c_k\rangle$ and
$\langle c_i^\dagger c_j c_k^\dagger c_l\rangle$.

## 6. Dynamical (frequency-dependent) correlators

The central quantity is a retarded correlator (Green's function) between
two operators $A,B$,

$$G_{AB}(\omega)=\langle\mathrm{GS}|A\,\frac{1}{\omega-H+E_0+i\delta}\,B|\mathrm{GS}\rangle,\qquad S_{AB}(\omega)=-\frac1\pi\,\mathrm{Im}\,G_{AB}(\omega)$$

where $\delta$ is a small broadening. Choosing $A=B=S^z_i$ at the same
site gives the local dynamical spin structure factor $S^{zz}_{ii}(\omega)$
(what a local probe like NMR/ESR couples to); choosing $A=S^z_i$,
$B=S^z_j$ at different sites and Fourier-transforming over $i-j$ gives
the momentum-resolved dynamical structure factor $S(q,\omega)$ measured
in inelastic neutron scattering. All five submodes below compute
$G_{AB}(\omega)$ (or equivalently $S_{AB}(\omega)$); they differ only in
*how*, and therefore in what energy range/resolution/cost trade-off they
offer:

```python
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                      name=(sc.Sz[0], sc.Sz[0]))
```

**`submode="KPM"` — Kernel Polynomial Method.** $H$ is linearly rescaled
into $\tilde H=(H-b)/a$ so its full spectrum lies in $[-1,1]$
(`kpm_scale`, default 0.7, sets how much margin is left inside $[-1,1]$),
then $G_{AB}$'s spectral function is expanded in Chebyshev polynomials
$T_m(\tilde H)$:

$$\mu_m=\langle\mathrm{GS}|A\,T_m(\tilde H)\,B|\mathrm{GS}\rangle,\qquad S_{AB}(\tilde\omega)=\frac{1}{\pi\sqrt{1-\tilde\omega^2}}\left[\mu_0+2\sum_{m=1}^{N-1}g_m\,\mu_m\,T_m(\tilde\omega)\right]$$

with Jackson-kernel damping coefficients $g_m$ suppressing Gibbs
ringing from the truncation at $N$ moments (`kpm_n_scale` scales $N$,
i.e.\ the energy resolution, relative to the rescaled bandwidth). This is
the default, general-purpose method: robust, works across the whole
spectrum at once, cost grows only linearly with the number of moments.

**`submode="CVM"` — correction-vector method.** Instead of a global
polynomial expansion, this solves directly for the correction vector at
one frequency $\omega$ at a time, via the positive-definite linear system

$$\big[(H-\omega-E_0)^2+\eta^2\big]\,x_c=-\eta\,B|\mathrm{GS}\rangle$$

(solved by conjugate gradient in MPS form), from which

$$x=(\omega+E_0+i\eta-H)^{-1}B|\mathrm{GS}\rangle=i\,x_c+\frac{H-\omega-E_0}{\eta}\,x_c,\qquad G_{AB}(\omega)=\langle\mathrm{GS}|A|x\rangle$$

Here $\eta$ (`delta`) is the artificial broadening that regularizes the
resolvent at a real frequency. CVM is more accurate at a single targeted
frequency/energy window (e.g. zooming in on a sharp resonance) than a
global KPM expansion, at the cost of re-solving the linear system for
every $\omega$ on the requested grid.

**`submode="TD"` — time-dependent DMRG.** Real-time evolution gives the
correlator directly in the time domain,

$$C(t)=\langle\mathrm{GS}|A(t)B(0)|\mathrm{GS}\rangle=e^{iE_0t}\langle\mathrm{GS}|A\,e^{-iHt}\,B|\mathrm{GS}\rangle$$

which is then windowed with an exponential damping factor
$e^{-\delta t}$ (equivalent to a Lorentzian broadening of width $\delta$
in frequency) and Fourier transformed,

$$S_{AB}(\omega)=\frac1\pi\,\mathrm{Re}\!\int_0^{T}\!dt\;e^{i\omega t}\,C(t)\,e^{-\delta t}$$

The total simulated time $T$ (`damping_periods`/$\delta$) must be long
enough that the damping has suppressed truncation ringing by $t=T$.
Frequency resolution is set by $T$ (via the usual $\Delta\omega\sim
1/T$ time-frequency uncertainty), so this method is best when you want
fine resolution over a *narrow* frequency window, at the cost of a real
dynamical simulation (bond dimension grows with entanglement generated
during the evolution).

**`submode="TDZ"` — complex-time evolution (Cao, Lu, Stoudenmire &
Parcollet, arXiv:2311.10909).** Real-time evolution (`"TD"` above) grows
entanglement, so the MPS bond dimension needed for a given accuracy
grows with the simulated time $T$. TDZ instead evolves along a complex
time contour

$$z(t,\alpha_0)=\int_0^t e^{-i\alpha_0 f(t')}\,dt',\qquad f(t)=e^{-t\omega_0},\qquad \omega_0=2\pi/t_{\max}$$

Since $\mathrm{Im}\,z(t,\alpha_0)<0$ for $\alpha_0>0$, this progressively
damps high-energy content as it evolves, so the bond dimension needed
for a given accuracy grows far more slowly than under real-time
evolution alone (the original paper reports $\chi\sim20$–$30$ vs
$\chi\sim500$–$700$ for comparable accuracy on the Anderson impurity
model). The true real-time ($\alpha_0=0$) correlator is then recovered
order by order via a perturbative Taylor expansion in $\alpha_0$ around
the simulated contour,

$$C(t,0)\approx\phi^{(0)}(t,\alpha_0)+\sum_{n=1}^{n_{\max}}g^{(n)}(t,\alpha_0)$$

where $\phi^{(n)}(t,\alpha_0)=\langle H^n B|\mathrm{GS}\rangle\cdot
|\psi(t,\alpha_0)\rangle$ (precomputed once per $n$, reused as a fixed
overlap target at every time step) and $g^{(n)}$ are explicit
combinatorial expressions in $\phi^{(1..n)}$ and the pure contour
integrals $J^{(n)}(t,\alpha_0)=-i\,\partial^n_{\alpha_0}z(t,\alpha_0)$
(see the paper's Appendix B; this implementation hardcodes $n\le4$,
which the paper finds always suffices for $\alpha_0\lesssim0.3$). The
reconstructed $C(t,0)$ is then windowed/Fourier-transformed exactly as
in `"TD"`.

```python
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="TDZ",
                                      name=(sc.Sz[0], sc.Sz[0]),
                                      alpha0=0.1, n_max=4, dt=0.05)
```

`alpha0` is the contour angle parameter (larger reduces the bond
dimension needed further, but requires a larger `n_max` to reconstruct
the real axis accurately); `n_max` (≤4) is the reconstruction order;
`dt`/`tmax`/`nt` set the underlying time step/duration exactly as in
`"TD"`. Uses two-site TDVP when available (`itensor_version` 3 or
`"python"`, `tevol_method="TDVP"`, the paper's own setup) or falls back
to the MPO-Taylor propagator otherwise (`itensor_version=2`, which has
no TDVP) — the same TDVP-vs-Taylor choice `"TD"` already makes. Current
scope: only the "greater" branch of the correlator is computed (the same
simplification `"TD"` itself already makes), so this is best used the
same way as `"TD"`: high-resolution work in a narrow frequency window,
now reachable at a lower bond-dimension cost for a given simulated time.

**`submode="EX"` — exact diagonalization in a truncated DMRG subspace.**
Builds $A$, $B$, $H$ explicitly in the subspace spanned by the lowest
`nex` (orthogonalized) DMRG excited states, then evaluates the exact
Lehmann sum in that subspace,

$$G_{AB}(\omega)=\sum_{n=1}^{n_{ex}}\frac{\langle\mathrm{GS}|A|n\rangle\langle n|B|\mathrm{GS}\rangle}{\omega-E_n+i\delta}$$

i.e.\ a small, explicit sum over poles at the computed excited-state
energies $E_n$, each with residue given by the transition matrix
elements. Cheap and exact *within* the truncated subspace; only as good
as how many/which excited states were computed.

**`submode="maxent"` — maximum-entropy reconstruction.** Reconstructs a
positive-definite spectral function from a finite set of moments
$\langle(H-E_0)^k\rangle$ using a maximum-entropy method
(`distribution.get_distribution_maxent`), rather than a Chebyshev
expansion — useful when positivity of the reconstructed $S(\omega)$
matters more than matching KPM's polynomial-expansion artifacts.

**Choosing a method:** KPM (default) for a first look at the full
spectrum; CVM or TD when you need high resolution in a specific,
narrow frequency window; TDZ instead of TD when that window also needs
long simulated times/low frequencies, where TD's real-time bond-dimension
growth becomes limiting; EX when a handful of excited states already
capture the physics (e.g. a small gapped system); maxent when you want a
guaranteed-positive reconstruction from limited moment data (e.g.\
combined with finite-temperature ED, see §9).

## 7. Real-time dynamics: quenches

Beyond frequency-domain correlators, DMRGPY directly simulates real-time
unitary evolution $|\psi(t)\rangle=e^{-iHt}|\psi(0)\rangle$ and measures
an observable along the way:

```python
from dmrgpy import timedependent
wf0 = sc.get_gs()                          # prepare |GS> of H0
sc.set_hamiltonian(h1)                     # quench to a different H
(ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0],
                                             nt=200, dt=1e-2, wf=wf0)
```

giving $\langle\psi(t)|S^z_0|\psi(t)\rangle$ as a function of time after
a Hamiltonian quench $H_0\to H_1$ — the standard non-equilibrium quantum
quench setup: prepare the ground state of one Hamiltonian (e.g.\ a
symmetry-broken/Néel-ordered $H_0$), then let it evolve under a
different $H_1$ (e.g.\ the isotropic Heisenberg point) and watch
observables relax, oscillate, or thermalize. `evolution_ABA` similarly
lets you apply an operator $A$ as an instantaneous local quench (e.g.\
flip a spin, or add a particle) before evolving and measuring $B(t)$.
`timeevolution.imaginary_exponential` computes autocorrelation functions
directly, $\langle\psi_0|e^{iHt}|\psi_0\rangle$, without a separate
measurement operator.

## 8. Density of states

**Many-body density of states.** The full many-body spectral density

$$\rho(E)=\sum_n\delta(E-E_n)$$

(sum over *all* eigenstates of $H$, not resolved by any operator) is
obtained via the same KPM machinery as §6 applied directly to $H$
itself:

```python
from dmrgpy import dos
es, rho = dos.get_dos(sc)
```

**Single-particle local density of states.** For fermionic chains, the
physically distinct quantity usually meant by "DOS" (e.g.\ as measured
by scanning tunneling spectroscopy) is the local single-particle
spectral function at site $i$,

$$A_i(\omega)=-\frac1\pi\,\mathrm{Im}\,G^R_{ii}(\omega),\qquad G^R_{ii}(\omega)\ \text{built from}\ \langle c_i(t)c_i^\dagger(0)\rangle\ (\omega>0)\ \text{and}\ \langle c_i^\dagger(t)c_i(0)\rangle\ (\omega<0)$$

i.e.\ the particle-addition and particle-removal (electron/hole)
branches of the local Green's function concatenated across $\omega=0$.
There is no ready-made helper for this (the module this section used to
point to, `fermiondos.py`, relied on a `get_dynamical_correlator` calling
convention that no longer exists and was removed) -- build the two
branches directly with `get_dynamical_correlator` instead, using
`(fc.C[i], fc.Cdag[i])` for the particle-removal branch and
`(fc.Cdag[i], fc.C[i])` for the particle-addition branch, then
concatenate.

## 9. Finite temperature

Thermal (mixed-state, finite-$T$) expectation values are obtained via
**purification**: each physical site is paired with an ancilla site, and
the maximally-entangled state of every physical-ancilla pair (e.g.\ the
singlet-forming Heisenberg coupling $\mathbf S_i^{\rm phys}\cdot\mathbf
S_i^{\rm anc}$ as the "Hamiltonian" preparing that state) is exactly the
purification of the infinite-temperature ($T=\infty$) physical density
matrix $\rho\propto\mathbb 1$. Imaginary-time evolving this purified
state under the *physical* Hamiltonian $H$,

$$|\Psi(\beta)\rangle=e^{-\beta H/2}|\Psi(0)\rangle,\qquad \beta=1/T$$

(stepped via repeated small first-order updates
$|\Psi(\beta+\Delta\beta)\rangle\propto(1-\Delta\beta\,H)|\Psi(\beta)\rangle$,
renormalized at each step) and tracing out the ancillas gives exactly the
thermal (Gibbs) density matrix of the physical chain,
$\rho(T)=\mathrm{Tr}_{\rm anc}|\Psi(\beta)\rangle\langle\Psi(\beta)|\propto e^{-\beta H}$:

```python
from dmrgpy import thermal
tc = thermal.Thermal_Spin_Chain(spins, T=0.1)
```

`T=0` recovers ordinary ground-state DMRG. For small systems ($n\lesssim
14$), `get_correlation_matrix(T=...)` instead computes the exact thermal
average directly by brute-force ED, $\rho=\sum_nZ^{-1}e^{-E_n/T}|n\rangle\langle n|$ — a useful cross-check of the purification approach.

## 10. Topological invariants

**Many-body Berry phase.** For a ground state that depends on an
adiabatic parameter $k$ threaded around a closed loop (e.g.\ inserted
flux), the discretized Berry phase

$$\gamma=\arg\prod_{k}\langle\psi(k)|\psi(k+\delta k)\rangle$$

quantized (typically in units of $\pi$) signals a topologically
nontrivial ground state / obstruction to adiabatic continuity — the
many-body generalization of a Zak/Berry phase, computed by running DMRG
at a discrete set of parameter points and chaining ground-state
overlaps.

**Single-particle Berry phase / Wilson loop.** For a translationally
invariant single-particle Hamiltonian $H(k)$ with a set of occupied
bands, the non-Abelian Wilson loop over one Brillouin-zone circuit is

$$W=\det\Big[\textstyle\prod_k U(k,k+\delta k)\Big],\qquad U_{mn}(k,k+\delta k)=\langle u_m(k)|u_n(k+\delta k)\rangle$$

with $m,n$ running over occupied bands; $\gamma=\arg W$ is the
polarization/Zak phase of that band manifold (and, combined with a scan
over a second momentum direction, the ingredient for a Chern number).

## 11. Mean-field decoupling

For a spin Hamiltonian with exchange couplings $J_{ij}$, a
self-consistent mean-field (Weiss-field) decoupling replaces the
two-body exchange term with a one-body field,

$$H_{\rm MF}=\sum_i\mathbf h_i\cdot\mathbf S_i,\qquad \mathbf h_i=(1-p)\sum_jJ_{ij}\langle\mathbf S_j\rangle$$

iterated to self-consistency: solve for $\langle\mathbf S_i\rangle$ under
$H_{\rm MF}$ (by DMRG), rebuild $\mathbf h_i$ from the new expectation
values, mix old/new fields, and repeat until $\max_i|\Delta\langle\mathbf
S_i\rangle|$ falls below a tolerance. The parameter $p\in[0,1]$
interpolates between pure mean-field theory ($p=0$) and keeping a
fraction $p$ of the original many-body exchange treated exactly
alongside the self-consistent field (useful for hybrid
mean-field-plus-fluctuations treatments of, e.g., magnetically ordered
phases where pure MF overestimates order).

```python
from dmrgpy import meanfield
meanfield.spinchain_meanfield(sc, p=0.0)
```

## 12. Fidelity susceptibility and quantum phase transitions

For a Hamiltonian $H(\lambda)=H_0+\lambda H_1$ depending on a tuning
parameter $\lambda$, the fidelity susceptibility measures how sharply
the ground state changes as $\lambda$ varies — and diverges at a quantum
phase transition, where the gap to the first excited state closes:

$$\chi(\lambda)=\sum_{n\neq0}\frac{|\langle 0|H_1|n\rangle|^2}{(E_n-E_0)^2+\delta^2}\qquad(\text{perturbative form, }\texttt{fmode="PT"})$$

`get_fidelity`'s **default** is not this perturbative form but a
non-perturbative estimator (`fmode="derivative"`) that computes $\chi$
directly from finite differences of the ground-state overlap matrix
between nearby $\lambda$ values, which also handles a (near-)degenerate
ground-state manifold via a smooth gauge choice:

```python
from dmrgpy import fidelity
chi = fidelity.get_fidelity(sc, h0, h1, lam, n=3) # fmode="derivative" (default)
chi_pt = fidelity.get_fidelity(sc, h0, h1, lam, n=3, fmode="PT") # perturbative form above
```

Scanning $\lambda$ and plotting $\chi(\lambda)$, a peak (sharpening and
diverging with system size) locates a quantum critical point without
needing to know its universality class in advance.

## 13. Ground-state degeneracy

Exact and near (e.g.\ symmetry-protected, or finite-size-split
topological) ground-state degeneracies are estimated by a narrow,
super-Gaussian-broadened level count around the ground-state energy
(`degeneracy.py`'s `gs_degeneracy_simple`/`eigenvalue_degeneracy`):

$$g(E_0)\approx\sum_i\exp\!\left[-\left(\frac{(E_i-E_0)^2}{\delta}\right)^2\right]$$

Note the quartic falloff in $(E_i-E_0)$ (not a plain Gaussian) — this
makes the window considerably narrower than $\delta$ would suggest for a
standard Gaussian, so `delta` needs to be picked accordingly (typically
larger than the target energy resolution) to count near-degenerate
levels rather than only exactly-degenerate ones.

summed over a growing number of low-lying computed eigenstates $E_i$
until the count converges — a value near an integer $g$ signals a
$g$-fold degenerate (or near-degenerate, at the working precision
$\delta$) ground-state manifold, as expected e.g.\ for the four
edge-state-split ground states of an open Haldane chain, or a
symmetry-broken ordered phase.

```python
from dmrgpy import degeneracy
g = sc.get_gs_degeneracy()
```

## 14. Reduced density matrices and operator distributions

**Reduced density matrix.** The single-site reduced density matrix
$\rho_i=\mathrm{Tr}_{j\neq i}|\mathrm{GS}\rangle\langle\mathrm{GS}|$ is
available directly (`sc.get_rdm(i=0)`), the basic object entanglement
entropies, local observables, and further post-processing (e.g.\ local
susceptibilities) are built from.

**Operator distributions.** More generally, the full probability
distribution of an arbitrary Hermitian operator $X$ (not just its
expectation value) in the ground state,

$$P(x)=\langle\mathrm{GS}|\,\delta(X-x)\,|\mathrm{GS}\rangle$$

is computed via the same KPM machinery as the dynamical correlators of
§6 (`sc.get_distribution`), or reconstructed from a finite set of raw
moments $\langle X^k\rangle$ via maximum entropy
(`get_distribution_maxent`). Useful for e.g.\ full counting statistics of
a conserved charge, or distinguishing a sharply peaked (well-defined
quantum number) ground state from a broadly spread one.

## 15. Post-processing tools

- **Analytic continuation** (`analyticcontinuation.py`): Padé
  continuation of a correlator known on the imaginary/complex-frequency
  axis (e.g.\ from a Matsubara-like or complex-shifted CVM calculation,
  §6's `submode="CVMimag"`) to the real frequency axis, where the
  physical spectral function lives.
- **Function fitting** (`functionfit.py`): a generic multi-start Powell
  minimizer used e.g.\ to fit the Calabrese-Cardy entropy formula in §5.
- **Finite-size extrapolation** (`extrapolate.py`): polynomial
  extrapolation in $1/L$ of a size-dependent quantity $y(L)$ toward the
  thermodynamic limit $L\to\infty$ — standard practice for extracting
  bulk quantities (energy density, order parameters, gaps) from finite
  DMRG chains.
- **Maximum-entropy reconstruction** (`reconstruct.py`): reconstructs a
  positive spectral function from a truncated moment expansion,
  underlying both the `"maxent"` dynamical-correlator submode and
  `get_distribution_maxent`.

## 16. Worked-example cookbook

**Central charge of a critical transverse-field Ising chain**
($H=\sum_iS_i^zS_{i+1}^z+\tfrac12\sum_iS_i^x$ at the critical field,
expected $c=\tfrac12$):

```python
sc.maxm = 200          # larger bond dimension: needed at criticality
wf = sc.get_gs()
print(wf.get_CFT_central_charge())
```

**Haldane gap of the spin-1 chain** (see §4):

```python
spins = ["S=1" for i in range(n)]
sc = spinchain.Spin_Chain(spins)
# ... build Heisenberg h ...
es = sc.get_excited(n=6)
print("Haldane gap:", es[4]-es[0])
```

**Charge gap of a Hubbard chain** (see §4):

```python
print("Single-particle gap:", fc.get_gap())
print("Charge (pair) gap:", fc.get_charge_gap(d=2))
```

**Fidelity susceptibility across the Ising transition** (see §12):

```python
h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))   # Ising coupling
h1 = 2*sum(sc.Sx[i] for i in range(n))                # transverse field
for lam in lambdas:
    chis.append(fidelity.get_fidelity(sc, h0, h1, lam, n=3))
```

**Momentum- and frequency-resolved dynamical structure factor**
$S(q,\omega)$ (see §6), by combining the site-resolved KPM correlator
with a lattice Fourier transform:

```python
Sqw = {}
for i in range(n):
    for j in range(n):
        x, y = sc.get_dynamical_correlator(submode="KPM", name=(sc.Sz[i], sc.Sz[j]))
        Sqw[(i, j)] = (x, y)          # combine with sum_ij e^{iq(i-j)} S_ij(w) offline
```
