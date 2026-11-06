# Beyond third order in the STM Kondo spectrum

Plan, 2026-11-06. Started in the session of 2026-11-05, which ran out of
memory before writing anything down (a dense `eigh` of a 38610-state
sector); this file collects what that session found and the checks run
since.

**Status (2026-11-06): Route B is chosen and implemented on the ED
spectrum at T=0, see `docs/new_methods/kondo_leading_log/README.md`.**
Route A is not implemented. The key identity behind the
non-perturbative route, the composite-operator form of dI/dV, is checked
on an exact discrete bath at zeroth order in J (3e-17) and at first order
in J (below 4e-8), see "The check that licenses the composite operator".

## 1. The question

`get_kondo_spectrum` is the third-order theory of Ternes, NJP 17 063016
(2015), arXiv:1505.04430. Second order is the golden rule for a tip
electron that scatters on the impurity spin, and third order adds one
exchange scattering with a sample electron, which gives the logs
F(eV - eps) at every spin excitation. The question is how to go beyond
that in a way that is physically meaningful, with the aim of computing
dI/dV spectra where third order fails.

Let us first say why the bare fourth order is not the answer. Each new
order adds one more power of the same log, rho J ln(omega0/|x|). So the
expansion parameter is rho J ln(omega0/|x|), not rho J, and it passes 1
near T_K. A fourth-order term adds nothing where the log is large, which
is exactly where third order already fails, and it improves things only
where third order was already adequate. What one wants is either a
resummation of the logs (the leading logs sum to a pole at T_K), or a
calculation that does not expand in J at all.

## 2. Two routes

There are two physically meaningful routes. They are complementary, not
alternatives.

**Route A, non-perturbative: dI/dV as the spectral function of a
composite operator.** Ternes' model has two small couplings: the tip
coupling T0 and the sample exchange J. The experiment is in the tunneling
limit, so T0 is small in the experiment, while J is what gives rise to the Kondo
physics. The theory we want is therefore exact in J and lowest order in
T0. At order T0^2 the current is a golden rule, but between the exact
eigenstates of the sample with the impurity, not the bare ones. The tip
electron enters through the tunneling vertex

    X_l^dag = sum_l' psi0^dag_l' tau(l', l),   tau = sigma.S/2 + U,

where psi0 is the sample electron at the impurity site, S the impurity
spin and U the potential scattering ratio. So at T=0 and a flat tip
density of states

    dI/dV(V) = 2 pi T0^2 * 2 * sum_l A_l(eV),
    A_l(w)   = sum_n |<n|X_l^dag|GS>|^2 delta(w - (E_n - E_GS))   (eV > 0),

with the removal part, X_l instead of X_l^dag, for eV < 0, and the factor
2 the spin-average ("SA") factor of `conductance.py`'s docstring. This is
a dynamical correlator of a composite operator (a fermion at the first
Wilson-chain site times the impurity spin) in the ground state of
impurity plus Wilson chain, which dmrgpy already knows how to compute.
Expanding A_l in J gives back the second order at J^0 and the third order
at J^1, and the next orders after that. Beyond that it contains what no
finite order has: the Kondo resonance below T_K, the renormalized
anisotropy and g-factor, the Korringa broadening of the steps. Jacob,
arXiv:1805.10624, states the same ideal-STM-limit result for the
molecule-substrate spectral function; the composite operator for exchange
tunneling we derive ourselves, since no source was found that writes it.

What it does not contain, in one clause each: the nonequilibrium
populations driven by the tip (spin pumping, Loth et al. Nat. Phys. 6,
340), which are order T0^4 and need the tip current to compete with
relaxation; finite temperature on the DMRG side (ED handles any T, DMRG
needs purification or METTS); and the continuum, which a Wilson chain
replaces by a logarithmic set of levels, so the chain must be long enough
for its lowest energy scale to be below T_K and below the smallest
splitting of interest, and the delta peaks must be broadened (or the
integrated current compared, which needs no broadening).

**Route B, perturbative: leading-log resummation.** Keep the structure of
`conductance.py` and replace the bare Jrho F(x) by a frequency-dependent
renormalized coupling g_fm(w), one per transition of the impurity
spectrum, flowing with logs that are cut off by the splittings and by a
decoherence rate Gamma (the Korringa rate, of order (rho J)^2 times the
splitting), instead of the fixed Gamma0. For S=1/2 at B=0 this sums
J + J^2 F + ... into J / (1 - J F). It reuses every object
`conductance.py` already builds, works at any T as cheaply as the third
order, and fails at T_K by construction, but it is controlled where it
applies, which third order is not near the thresholds. The structure is
Rosch, Paaske, Kroha, Woelfle (cond-mat/0202404, cond-mat/0408506), the
decoherence cutoff Paaske et al. (cond-mat/0401180), a multi-level
anisotropic spin Hoerig, Schuricht, Andergassen (1110.6103), and the
controlled two-loop check Schoeller and Reininghaus (0902.1446).

Route A is the reference for Route B: where both apply, the resummed
logs must agree with the composite-operator spectrum at weak coupling.

## 3. The check that licenses the composite operator

The script is `docs/kondo_higher_order/composite_check.py`. It builds the
impurity (S=1, D=-1 meV, B=(1,0,2) T, g=2, U=0.25) coupled by
J psi0^dag (sigma/2) psi0 . S to a discrete particle-hole symmetric bath
with three positive levels and their mirror images (12 fermion modes,
sectors of 2772 and 2376 states), diagonalizes the N and N+1 sectors
exactly and compares the cumulative weight of X_l^dag below V, which is
the current I(V) of the discrete model, with `conductance.py`'s second-
and third-order coefficients where the continuum F is replaced by its
discrete counterpart sum_q w_q / (eps_q - x). The bias points sit midway
between thresholds, so the comparison needs no broadening.

The runs were made on d22 (`/l/ladovj1/kondo_ho`), with the derivative in
J taken as a central difference.

- Zeroth order: the weight agrees with `second_order_dIdV`'s structure to
  3e-17, against a scale of 0.06, so the operator and the spin-average
  factor 2 are right.
- First order at U=0: agreement to 3.3e-8 against a scale of 0.065, in
  every window. So the composite operator reproduces the third-order
  Kondo term, including the exchange-diagram argument F(eV - (e_f - e_m)).
- First order at U=0.25: agreement to 4e-9 and 3.5e-8 on two bath
  discretizations, once one more term is added (below). The residuals
  are the dJ^2 error of the finite difference (dJ=2e-5); dJ=2e-3 and
  2e-4 gave the same mismatch without the term, so it is not a
  nonlinearity.

The term is the interference of the second-order spin amplitude with the
third-order amplitude that has U at the tunneling vertex and the spin at
the exchange vertex,

    U * sum_a |<f|S_a|i>|^2 * [F_electron(eV) - F_hole(eV)],   any f,

where F_electron and F_hole are the direct (electron-like) and exchange
(hole-like) logs, both at the incoming energy. `conductance.py` keeps
only the f=i potential loop. With Ternes' symmetric F0 the two logs are
equal and this term cancels exactly, so `conductance.py` is right within
its own convention. On a real band, and on a discrete bath, the two
differ at order eV/D with no log, the same order at which F0 differs
from the paper's defining integral eq. 20. This does not contradict
`tests/test_kondo_spectrum_tmatrix.py`, which agrees with `conductance.py`
to roundoff with U in the vertices: that test replaces the band integral
over the empty and over the full half by the same symmetric F0, so it
cannot see this term by construction. On a discrete bath the term
can be large: with a near-degenerate pair of final states 9 microeV apart
it gave a 0.26 error in the cumulative weight before it was added. For
the plan this means two things. The composite operator contains the
third order exactly, but it is not the F0 convention, so the comparison
with the continuum third order (step 3 of the ladder) holds at leading
log only. And a non-log piece of order eV/D in the U term is a real
property of the model, which the third-order code drops on purpose.

## 4. Validation ladder

Each step is a test before the next one starts.

1. J^0 of the composite spectral function against `second_order_dIdV`,
   on the discrete bath (done, 3e-17).
2. J^1 against `third_order_kondo_dIdV` plus `third_order_potential_dIdV`
   on the discrete bath, with the U term that cancels under F0 added
   (done, below 4e-8). This also checks the exchange-diagram argument
   F(eV - (e_f - e_m)), which `conductance.py` takes from a brute-force
   T-matrix and not from the paper.
3. The Wilson chain at weak J against the continuum third order: the
   integrated current I(V) of A against the integral of
   `get_kondo_spectrum(order=3)`, with omega0 set by the chain's band
   half-width D. The constant under the log differs between a flat band
   of half-width D and Ternes' F0, so this agrees at leading log only,
   and the plan is to fix the map of omega0 to D here, once.
4. S=1/2, B=0, at increasing J: the zero-bias peak against the universal
   Kondo spectral function, and with U the Fano shape (Zitko,
   arXiv:1110.4561).
5. S=2 with anisotropy, Fe on CuN, in fields along three directions,
   against the NRG of Zitko and Pruschke, arXiv:1002.4082.
6. Route B against Route A in the weak-coupling window, and both against
   third order where the logs are small.

## 5. Implementation milestones, Route A

- A1. The model builder: an impurity spin with anisotropy and field,
  coupled by J to the first site of a Wilson chain from
  `wilsonchain.wilson_chain`, on `Mixed_Spin_Fermion_Chain` (a genuine
  spin-S site next to spinful fermion sites, `itensor_version=3` and
  `"python"`). That class has no ED backend, so the ED reference of steps
  1 and 2 stays the explicit Fock-times-spin construction of the check
  script, moved into a test.
- A2. The composite operator X_l^dag as a `MultiOperator`
  (Cdagup[1]*Sz[0] and so on, with the U term), and its spectral function
  through the existing correlator submodes. The ones to try first are
  TASK (a Krylov spectrum, exact at saturation, already wired into the
  Kondo route) and SECTOR (the N+1 sector is exactly where X^dag|GS>
  lands). KPM is the fallback. Test: step 3 of the ladder.
- A3. Broadening of a logarithmic discretization: log-Gaussian broadening
  of the delta peaks, or z-averaging over several discretizations, as NRG
  does. Until then, compare integrated currents.
- A4. The public entry: a new order or submode of `get_kondo_spectrum`,
  or a separate method on the Wilson-chain side (see "Decisions").
- A5. Finite temperature: ED on small baths at any T first; DMRG later
  through purification, if it is needed.

## 6. Implementation milestones, Route B

- B1. Derive the one-loop flow of the couplings for a general impurity
  spectrum in the eigenbasis that `KondoSpectrum` already holds, with
  frequency-dependent couplings and the decoherence cutoff, and check
  that its expansion to J^2 reproduces the existing third-order
  coefficients (`_triple_product_coefficients`, the potential loop).
- B2. Implement it as `order="ll"` (or similar) of `get_kondo_spectrum`
  on the ED route; the DMRG routes through TASK inherit it, since TASK
  feeds `conductance.py` a `KondoSpectrum`-shaped object.
- B3. Validate against Route A, step 6.

## 7. Papers to read

From a survey of papers on 2026-11-05. It read abstracts and search
snippets only, no equations; "useful for" means where to look, not
verified content.

- cond-mat/0202404, Rosch, Paaske, Kroha, Woelfle, PRL 90, 076804
  (2003): perturbative RG with frequency-dependent couplings for a Kondo
  dot at large V and B; the template for Route B.
- cond-mat/0307365, Paaske, Rosch, Woelfle, PRB 69, 155330 (2004):
  nonequilibrium perturbation theory for the Kondo dot in a field.
- cond-mat/0401180, Paaske, Rosch, Kroha, Woelfle, PRB 70, 155301 (2004):
  decoherence rates that cut off the logs; the source of Gamma in Route B.
- cond-mat/0408506, Rosch, Paaske, Kroha, Woelfle, JPSJ 74, 118 (2005):
  the longer account of the perturbative RG.
- 0902.1446, Schoeller and Reininghaus, PRB 80, 045117 (2009): two-loop
  real-time RG for the anisotropic Kondo model in a field; the controlled
  reference for Route B.
- 1110.6103, Hoerig, Schuricht, Andergassen: two-loop RG for a spin-1
  Kondo dot with anisotropy; closest to a multi-level spin-S impurity.
  Also 1402.0479 (Hoerig, Mora, Schuricht, PRB 89, 165411) for arbitrary S.
- 1201.6295, Pletyukhov and Schoeller, PRL 108, 260601 (2012), and
  1405.3150: RG from weak to strong coupling, within a few percent of NRG.
- 1104.2764, Hurley, Baadji, Sanvito, PRB 84, 115435 (2011): a third-order
  self-energy with NEGF for Co and Fe on CuN; the one independent
  third-order theory to cross-check against.
- 1805.10624, Jacob, J. Phys.: Condens. Matter 30, 354003 (2018): dI/dV
  in the ideal STM limit as the equilibrium spectral function; the
  statement Route A rests on.
- 1507.08474, Jacob and Fernandez-Rossier, EPJB 89, 210 (2016), and
  1712.01005, Jacob, PRB 97, 075428 (2018): OCA for anisotropic spins,
  split Kondo peaks, renormalized excitation energies; benchmarks.
- 1002.4082, Zitko and Pruschke, NJP 12, 063040 (2010): NRG for a spin-2
  anisotropic Kondo model, Fe on CuN, fields along three directions;
  benchmark for step 5.
- 1110.4561, Zitko, PRB 84, 195116 (2011): NRG universal Kondo spectral
  function with a Fano parameter, S=1/2; benchmark for step 4.
- 1112.6288, Korytar, Lorente, Gauyacq, PRB 85, 125434 (2012): sharp
  peaks at the inelastic thresholds with their own energy scale, which
  survive above T_K.
- 1401.7272, Delgado, Hirjibehedin, Fernandez-Rossier, Surf. Sci. 630,
  337 (2014): renormalized anisotropy and g-factor at second order in
  rho J; experiment 1407.1542, Oberg et al.
- 1103.3676, Delgado and Fernandez-Rossier, PRB 84, 045439 (2011):
  cotunneling theory, the second-order baseline.
- 1308.4524, Zhang et al., Nat. Commun. 4, 2110 (2013): spin-1/2 radical
  on Au(111) fitted with the third-order theory in T and B.
- 1701.03587, Ternes, Prog. Surf. Sci. 92, 83 (2017): review of single
  spins, Kondo correlations, dimers and chains.

No paper was found that takes Ternes' theory to fourth order, or that
computes STM dI/dV of a spin impurity with DMRG on a Wilson chain.

## 8. Decisions for the author

- Which route first: A (all orders, T=0, needs the chain and broadening)
  or B (cheap, any T, perturbative).
- The first target model: S=1/2 Kondo (the cleanest benchmark) or the
  anisotropic S=1 and S=2 cases the third order is used for today.
- Where Route A lives in the API: inside `get_kondo_spectrum` or as its
  own method on the Wilson-chain side.
