# dI/dV of a spin chain on a superconducting substrate

Plan, 2026-11-06.

**Decision (2026-11-06): start with the many-body route on one spin,
milestones S1 to S3** (pairing Wilson chain, composite operator, parity
helper), checked against the NRG numbers of arXiv:1008.1734. S0 and the
surrogate chain come after.

**Status (2026-11-07): S0 to S4, the surrogate half of S5, S6 and the
guide sections of S7 are implemented.** The code is `src/dmrgpy/ysr.py`
(one spin of any S with anisotropy), `src/dmrgpy/ysrtk/zaverage.py` (the
continuum above the gap), `src/dmrgpy/ysrtk/surrogate.py` (several spins),
`src/dmrgpy/ysrtk/bcsperturbative.py` (S0) and
`src/dmrgpy/ysrtk/ysrcurrent.py` (the YSR width of S6),
`src/dmrgpy/ysrtk/paritystate.py` (a ground state of definite parity for
the STM spectrum), `src/dmrgpy/ysrtk/subgapscan.py` (warm-started scans)
and the representative-energy mesh `mesh="zp"` of `wilsonchain.py`, with
`tests/test_ysr*.py`, ten examples under `examples/kondo/` and the guide
sections from "Spin excitations on a superconductor" on. The details and the measurements
are in `docs/development.md` 4.8f and 4.8g. In short:

- One spin 1/2: the pairing Wilson chain is exact against the BCS star,
  the classical YSR energy agrees to 5e-4 of Delta with a field v on the
  chain read as v/A_Lambda, the STM spectrum agrees with ED to 1e-11
  (TASK). Averaged over four twists z, Delta_c/T_K = 4.41 at rho J = 0.1
  against the NRG value 3.7, and 4.33 to 5.02 at rho J = 0.08 to 0.2. At
  rho J = 0.1 the absolute scales depend on the mesh: Delta_c and T_K are
  29% and 8% above the paper on the mean-energy mesh, 0.7% and 16.5% below
  it on the Zitko-Pruschke mesh, while the ratio is 4.40 on both. So the
  offset of the ratio is not the mesh; the first suspect is T_K from the
  local susceptibility instead of Wilson's chi_imp. The rho J -> 0 limit
  is 3.8 to 4.1 depending on the fit, so these data do not fix it.
- Parity: near the singlet-doublet transition a plain solve can stop at a
  mixture of the two parities. `parity=True` solves inside the parity
  sector, and the spectrum is then exact against ED (8 sites: one pole at
  0.0190312 Delta, elastic weight 0).
- Scans: a scan over Delta or Dz, or a bisection of Delta_c, starts each
  point from the states of the one before, about 2x faster per point on
  26 sites at the same energies to 1e-12.
- Higher spin: 2S native sites bound by a ferromagnetic Hund coupling,
  so that the Sz sectors work and S_imp^2 = S(S+1) is exact. Averaged over
  z = 0 and 1/2, Delta_c(S)/Delta_c(1/2) = 1.57, 2.72, 5.25 for S = 1,
  3/2, 2, against 1.59, 2.78, 5.35 in Table I. The renormalized anisotropy
  of eq. 2 agrees to about 0.01.
- Continuum: z-averaging with a log-Gaussian in xi = sqrt(w^2 - Delta^2)
  gives the BCS density of states to 4e-6 at nz=4, b=0.4. The Kondo
  resonance below about T_K/3 is not resolved at nkry=300.
- Several spins: the surrogate bath of arXiv:2307.11646 reproduces the
  published fit tables for L = 1 to 8, agrees with the Wilson chain to
  0.06 Delta at L=4 for Delta >= 0.1 D, and has no Kondo scale for
  Delta << D. Dimers by ED and chains of 4 and 8 spins on d22 give the
  YSR band with the weights of the open-chain modes.
- S0: Ternes' second order with Dynes densities of states on sample and
  tip gives back `second_order_dIdV` at Delta=0 to roundoff, the SIS peak
  weight to 0.2%, and the coherence peaks of Fe-OEP-Cl on Pb at 2Delta and
  2Delta + 2D. The substrate relaxation rate carries the BCS coherence
  factor and vanishes below 2Delta at T=0. No current-driven occupations.
- S6: the rate model of arXiv:1502.05048, with single-particle and
  Andreev currents, reduces to `ysr.tip_convolution` at weak coupling and
  takes eps, u^2, v^2 from the many-body poles. For a spin-1/2 doublet the
  crossover conductances are right only up to a factor of order 2.

Two choices differ from the plan below. The impurity is native spinful
sites with one electron each, not a `Mixed_Spin_Fermion_Chain` spin,
because only the native class has Sz quantum numbers; the parity sectors
are then values of 2Sz, with no parity operator. And Lambda=4 replaces
Lambda=2, since at Lambda=2 the energy differences of order Delta were not
converged at maxm=128. Open: the Wilson-ladder half of S5 (two channels
for a dimer), a long-chain high-spin dI/dV, the rho J -> 0 limit of
Delta_c/T_K (needs Wilson's chi_imp and more twists), a Lambda=2 run on
the zp mesh, current-driven spin occupations in S0, the degeneracy of a
spin-1/2 doublet in the rate model of S6, the slower "python" sector solve
of coupled surrogate spins, `parity=None` (a sector solve only when the
elastic weight trips), and an adaptive stop of the warm sweeps.

## 1. The question

We want the STM spectrum dI/dV(V) of a magnetic atom, or a chain of them,
adsorbed on a superconductor. The spectrum should contain both the Kondo
physics that `get_kondo_spectrum` gives on a normal metal and the physics
of the superconducting gap: the coherence peaks at the gap edge, the
Yu-Shiba-Rusinov (YSR) bound states inside the gap, the spin excitation
steps moved out to Delta plus the excitation energy, and the
screened-to-free quantum phase transition when T_K and Delta are
comparable.

The short answer from the literature is that no paper writes this as one
formalism for a spin chain. Every piece exists separately (see "What the
literature gives"), and the composite-operator route of
`docs/kondo_higher_order_plan.md` (its Route A) is the formalism that joins
them. That route is exact in the substrate coupling and lowest order in
the tip coupling, so the superconductor enters only through the bath
Hamiltonian, and the tip enters only through the operator whose spectral
function we compute. Route A is itself not implemented yet, so its
milestones A1 and A2 are prerequisites here and are not repeated.

## 2. The formalism

Let us start with one spin S on a substrate. The tip tunnels into the
substrate electron psi0 at the site of the spin, either directly (a
potential ratio U) or with a spin flip, through the vertex

    X_sigma^dag = sum_sigma' psi0^dag_sigma' [U delta + sigma.S/2]_{sigma' sigma}.

At lowest order in the tip coupling T0, with a normal tip of flat density
of states at T=0, the current is a golden rule between the exact
eigenstates of spin plus substrate, so

    dI/dV(V) ~ T0^2 sum_sigma A_sigma(eV),
    A_sigma(w) = sum_n |<n|X_sigma^dag|GS>|^2 delta(w - (E_n - E_GS)),  w > 0,

and the removal part with X_sigma for w < 0. This is a dynamical
correlator of a composite operator, which dmrgpy already computes for any
`MultiOperator`. Jacob, arXiv:1805.10624, states the same ideal-STM limit
for a molecule on a normal substrate; the check that the composite
operator gives back Ternes' second and third order is in
`docs/kondo_higher_order_plan.md` ("The check that licenses the composite
operator", agreement to 3e-17 and below 4e-8).

The superconductor enters only through the bath. For a BCS band

    H_band = sum_k eps_k c^dag_k c_k + Delta sum_k (c^dag_{k up} c^dag_{-k dn} + h.c.),

which is the model of Zitko, Bodensiek and Pruschke, arXiv:1008.1734.
The pairing term is diagonal in energy and the same for both spins, so the
Wilson map, a real orthogonal transformation of the energy shells applied
identically to both spins, leaves sum_q a^dag_{q up} a^dag_{q dn}
invariant. As a result the Wilson chain of a BCS band is the normal Wilson
chain plus a uniform onsite pairing Delta (f^dag_{n up} f^dag_{n dn} + h.c.)
on every chain site, and `wilsonchain.wilson_chain` needs one new term,
not a new discretization. This is a one-line argument, not a quoted
result; milestone S1 checks it numerically against the BCS bath Green
function.

What this gives, in one clause each: the YSR poles appear as discrete
peaks inside the gap with their electron-hole weight asymmetry (u^2
against v^2); the coherence peaks come from the Bogoliubov continuum of
the chain; the spin steps appear at Delta + eps_fi through the sigma.S
part of X; the Kondo resonance appears when Delta goes to zero; and the
quantum phase transition shows up as a change in ground state parity.
What it does not give: the width of a YSR peak at T=0, which in
experiment comes from quasiparticle relaxation (Ruby et al.,
arXiv:1502.05048) and here is a phenomenological broadening; Andreev
processes, which are order T0^4 and matter only at high tip conductance
(same paper); and finite temperature on the DMRG side.

A superconducting tip, which is what most experiments use for resolution,
is a convolution of the same A(w) with the tip BCS density of states,

    I(V) ~ integral dw rho_tip(w - eV) A(w) [f(w - eV) - f(w)],

so the YSR peaks move to +-(Delta_tip + eps_YSR). This is post-processing
and changes nothing in the many-body part.

## 3. A cheaper perturbative step

Before the many-body route there is a cheap one. Ternes' second order is a
golden rule into the substrate with a flat density of states. Replacing it
by the BCS density of states turns every inelastic step at eps_fi into a
coherence peak at eV = Delta + eps_fi, which is the "protection of excited
spin states" experiment of Heinrich, Braun, Pascual and Franke,
arXiv:1405.6984. This reuses `conductance.py::second_order_dIdV` with the
step function replaced by a convolution with rho_BCS, works at any T, and
costs nothing. It has no YSR states by construction, because a YSR state
is a bound state of the exchange coupling to all orders. The third order
with a BCS bath, the log cut off by Delta, is a possible second step, but
it fails exactly where the YSR state leaves the gap edge, so it is not
where the effort should go.

## 4. Geometry: from one spin to a chain

A full Wilson chain per spin is the reference for one spin. For a chain of
L spins it is too expensive: L x (n+1) sites with the exchange acting
across a whole Wilson chain in the MPS ordering. Two cheaper baths exist.

- The surrogate model of Baran, Frost and Paaske, arXiv:2307.11646: a few
  discrete BCS levels per impurity, with energies and couplings fitted to
  the Matsubara hybridization function. Already 1 level (the zero
  bandwidth limit) puts the singlet-doublet transition roughly right, and
  2 to 4 levels agree with NRG on the subgap energies. Odd numbers of
  levels contain a level at zero energy, which gets the singlet energy
  right at weak coupling, even numbers get the strong coupling side right.
  This is the bath that makes a chain feasible, since each spin then
  carries only 2 to 4 fermion sites. The authors point to public codes
  (FermPy by Frost, "Surrogate models" by Baran), which the
  physics-literature rule says to check before writing our own fit.
- A ladder: the spin chain on one leg, a one-dimensional BCS tight-binding
  chain on the other. It gives the YSR hybridization between neighbors
  with a real distance dependence (in one dimension, not the 2D or 3D of a
  real substrate), and it is directly a DMRG geometry.

For the surrogate chain the coupling between the baths of neighboring
spins is a hopping between their levels, which in a real substrate decays
with distance as sin(k_F r)/(k_F r) e^(-r/xi). In the surrogate it is a
parameter, and the plan says so wherever a chain result is shown.

## 5. What the code has to handle

- Delta breaks the particle number. So `set_conserved_sector(Nf=)` and
  `submode="SECTOR"` are out, and the routes for A(w) are TASK and KPM.
  Sz is still conserved at zero transverse field, so `Sz=` stays usable on
  `itensor_version=3`.
- The ground state changes fermion parity at the YSR transition, and
  dmrgpy has no parity sector. The YSR energy is E0(odd) - E0(even),
  taken from the ground state and the first excited state (overlap
  penalty), with each one's parity read off <(-1)^N>, since either can be
  the lower one. The parity operator is a product over all sites, so it
  should be built as a product MPO and not summed.
- `Mixed_Spin_Fermion_Chain` has no ED backend. The ED reference is the
  explicit Fock-times-spin construction of
  `docs/kondo_higher_order/composite_check.py` with Delta added to its
  discrete bath, which is exactly a surrogate model. It must be checked
  first that a pairing term `Cdagup[i]*Cdagdn[i]` builds on that class
  with a spin site inside the Jordan-Wigner string.
- `anderson_chain`'s rule "take n odd" is a normal-state rule. With Delta
  the degeneracy is set by the parity of the ground state, not by n, and
  the warning has to follow that.
- The gap cuts the Wilson chain. Sites with Lambda^(-n/2) D well below
  Delta only renormalize the gapped continuum, so the chain can stop
  there, which makes the superconducting case cheaper than the normal one
  for the same Lambda. Near the transition the ground state can depend on
  the discretization (z), as 1008.1734 warns, so results there need two z
  values.

## 6. Validation ladder

Each step is a test before the next one starts.

1. The bath alone: the local Green function of the first site of the
   pairing Wilson chain against the BCS one, gap and coherence peaks,
   after log-Gaussian broadening.
2. Classical spin limit: eps_YSR/Delta = (1 - alpha^2)/(1 + alpha^2),
   alpha = pi rho J S/2, for a static exchange field (no spin flips), on
   the Wilson chain.
3. S=1/2 Kondo spin: the transition at Delta_c/T_K = 3.7 with Wilson's
   definition of T_K (Table I of 1008.1734). The value 0.3 quoted
   elsewhere is T_K/Delta with a different T_K, so the comparison fixes
   the definition first.
4. Higher spin with anisotropy: Delta_c/T_K for S = 1, 3/2, 2 (5.9, 10.3,
   19.8 in the same table) and the subgap multiplets against D.
5. The surrogate bath against the full Wilson chain on one spin, at 1 to
   4 levels, both on subgap energies and on the dI/dV weights.
6. Two spins: YSR dimers against Schmid et al., arXiv:2203.15011, and the
   singlet limit of 1008.1734 section VII (subgap states pushed to the
   gap edge).
7. The chain: against the effective t-J model of Steiner, Mora, Franke
   and von Oppen, arXiv:2107.00031, as a check of trends only, since that
   model is itself an approximation.

## 7. Milestones

- S0. The perturbative step of "A cheaper perturbative step": BCS density
  of states in `second_order_dIdV`, normal or superconducting tip. Test
  against a direct golden-rule sum on a discrete BCS bath.
- S1. `wilson_chain(..., Delta=)` and `anderson_chain(..., Delta=)`, the
  pairing term, and ladder step 1.
- S2. Route A1 and A2 of `docs/kondo_higher_order_plan.md` (impurity spin
  on `Mixed_Spin_Fermion_Chain`, composite operator, TASK), run with
  Delta. Ladder steps 2 to 4. ED reference from `composite_check.py`
  with pairing.
- S3. Parity of a state and E0(odd) - E0(even) as a small helper, used by
  steps 3 and 4.
- S4. The surrogate bath (`wilsonchain.surrogate_bcs_bath` or similar),
  the fit to the hybridization function, ladder step 5.
- S5. Dimers and chains on the surrogate bath and on the ladder, ladder
  steps 6 and 7, with the tip at each spin.
- S6. Superconducting tip convolution and a phenomenological YSR width.
- S7. Public entry, docs page under `docs/guide/`, examples with plots.

## 8. What the literature gives

- Zitko, Bodensiek, Pruschke, arXiv:1008.1734: NRG of a high-spin Kondo
  impurity with anisotropy in a BCS host, Delta_c/T_K by spin and channel
  number. The model and the reference numbers for S2.
- Baran, Frost, Paaske, arXiv:2307.11646: the surrogate BCS bath, ED with
  a few levels, benchmarked against NRG and QMC. The bath for S4.
- von Oppen and Franke, arXiv:2105.06651: YSR states of real transition
  metal impurities, higher spin, several channels, crystal field.
- Ruby et al., arXiv:1502.05048: single-electron against Andreev tunneling
  into YSR states with a superconducting tip, and the relaxation rates
  that set the peak widths.
- Heinrich, Braun, Pascual, Franke, arXiv:1405.6984: spin excitations of
  molecules on Pb, steps moved by the gap, long excited-state lifetimes.
- Steiner, Mora, Franke, von Oppen, arXiv:2107.00031: an effective model
  for dilute YSR chains, with tunneling spectra.
- Schmid et al., arXiv:2203.15011: quantum YSR dimers, the change of
  ground state multiplicity with the binding of a quasiparticle.
- Trivini et al., arXiv:2207.00617: pair excitations of a quantum spin on
  a proximitized superconductor, peaks outside the gap at strong exchange.
- Experiments that combine YSR and spin-flip excitations in one spectrum:
  arXiv:1811.11591 (FePc), arXiv:2006.09459 (Fe on Nb).
- Codes: NRG Ljubljana (github.com/rokzitko/nrgljubljana) has
  superconducting bands and is the reference for S2; the surrogate
  authors' FermPy and "Surrogate models" repositories are the reference
  for S4.
