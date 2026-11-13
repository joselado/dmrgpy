# An NRG route to the Kondo spectrum of several spins, with matrix product states

Plan, 2026-11-13. The numbers quoted below are from the existing
docstrings and from the papers of the last section.

**Status (2026-11-13): N0 to N4 and N7 are implemented**, in
`nrgspectral.py`, `nrgmps.py` and `ysr.get_stm_spectrum(submode="NRG")`.
The gate N2 passed: a DMRG bond basis of about 150 states holds the low
states of the block well enough that the head and the tail give the exact
NRG of 3 spins to the 2% two NRG runs differ by, while 64 states were 3%
low. The exact route stops near 6 spins, as estimated in section 4 (keep
1600 is 6% off there). Ladder step 9 is open: at 7 and 8 spins the head
route runs and gives the expected physics, but maxm = 256 was not enough
to converge it (6% and 15% changes from 128). A flip-symmetric hand-over
for a head with a moment is open too. The measurements are in
`docs/new_methods/kondo_nrg_mps/README.md`; what follows is the plan as
written.

**Decision (2026-11-13): the bath on one spin.** The target is the
one-channel geometry of `spin_chain_on_metal`, milestones N0 to N4 and N7.
N6 (two and three channels) and the fork network stay out for now.

## 1. The question

Can we compute the Kondo spectrum the way the numerical renormalization
group (NRG) does it, with several spins instead of one, and with matrix
product states (MPS) instead of exact diagonalization?

The short answer is yes, subject to one gate that has to be measured first
(milestone N2). The way we do this is by splitting the Wilson chain in two. DMRG solves the spins and the first Wilson sites, the head.
An NRG iteration then continues from the bond basis of that state, the
tail: add one site, diagonalize, keep the low states, and repeat until the
spectrum stops changing. The head has no separation of energy scales, so it
needs the entanglement truncation of DMRG. The tail has it, so there the
energy truncation of NRG is correct and costs seconds per site.

This is not a new formalism. Every piece is published (see "Papers and
codes"). What we did not find in our searches is the combination applied to
a spin chain on a metal with an STM spectrum as the output.

## 2. What NRG is in MPS language, and what the code already has

NRG takes the impurity and the first site of the Wilson chain, diagonalizes,
keeps the lowest D states, adds the next site, and repeats. The hopping
falls as Lambda^(-n/2), so each new site is a small perturbation on the kept
states. This is what makes the truncation by energy correct.

The kept states of NRG are an MPS (Verstraete et al., cond-mat/0504305).
The matrix that rotates "kept states times new site" into the new kept
states is the MPS tensor of that site. So the DMRG ground state of
`ysr.kondo_chain` or `ysr.spin_chain_on_metal` is already the variational
form of the NRG ground state. Saberi, Weichselbaum and von Delft
(arXiv:0804.0193) showed that such a state also contains the NRG kept
states of every shell, when its bond dimension is large enough.

What the code does not have is the tower of excited states of each shell.
We have three related pieces:

- `nrgthermo.py` is an NRG with exact diagonalization for one spin 1/2. It
  gives thermodynamics only (T chi_imp, S_imp, Wilson's T_K). It drops the
  eigenvectors and the discarded states at every step, and it carries only
  f_N and S^z_imp.
- `ysr.get_stm_spectrum` and `ysrtk/zaverage.py` compute the spectrum from
  one DMRG state of a chain of fixed length, with one global Krylov space
  (TASK). The `zaverage` docstring says what limits it: at 21 sites and
  `nkry=300` the lowest pole sits near T_K/10 and carries all the weight
  below it. The Kondo resonance below about T_K/3 is not resolved, and the
  unitary value is not tested.
- `kondoscreening.py` reads ground-state expectation values on a Wilson
  chain. There is no NRG iteration in it.

So the fixed length is not the real limit of the present route. The limit
is that one Krylov space has to converge at every energy scale at once. NRG
solves one scale per site instead, and this is the piece we add.

The one-site tangent space of TASK and the complete basis of NRG are the
same construction (Gleis, Li, von Delft, arXiv:2207.13161): the discarded
states of each site times a basis of the rest of the chain. The difference
is how H is diagonalized in that space, globally by Krylov or shell by
shell. This is why the two routes can share code and can be compared on one
state.

## 3. The observable

With the substrate in the Hamiltonian, the tip is the only perturbation.
dI/dV is then the spectral function A(w) of the composite operator
X^dag_sigma of `ysr.tunneling_operators`, the T-matrix of the Kondo model
(Costi, cond-mat/0004302). This is the convention of `get_stm_spectrum`,
and the new route returns the same quantity in the same units.

The new route does not go through `conductance.py`. That code is Ternes'
perturbation theory, with the substrate put in by hand as `Jrho_s`. An
eigenbasis that already contains the Wilson chain, used there with a
nonzero `Jrho_s`, would count the substrate twice. No third-order term is
needed: the bath produces the Kondo logarithm to all orders.

## 4. Where NRG with exact diagonalization stops for several spins

Let us take ns spins with exchange Js, and the bath on one spin. The spin
energies are far below the band edge, Js/D of order 1e-3 for real atoms.
While the energy of the current shell is above the bandwidth of the spins,
every spin state is a low-energy state, and none can be cut. NRG then has
to keep all 2^ns spin states, each with its own tower of bath states.

For one spin, `nrgthermo` needs `keep=400` at Lambda=4, meaning roughly 200
bath states per impurity state. The same count gives about 3000 states for
4 spins and 13000 for 6. So we expect the exact route to stop near 6 spins.
This is an estimate; the first milestones measure `keep` against ns and
Js/D.

One shortcut looks natural and is not controlled: take the lowest states of
the bare spin chain from DMRG, and use them as the NRG impurity. The states
that are cut take part in the Kondo flow between D and the cut. A one-loop
estimate gives a correction of order rho J ln(D/E_cut) to the coupling,
which is not small for D/E_cut = 1000, and T_K depends exponentially on the
coupling. The plan has a direct test of this (milestone N1, 3 spins, all 8
states against the projected doublet). We do not build on this shortcut.

## 5. The method: a DMRG head and an NRG tail

The chain keeps the order of `spin_chain_on_metal`: the spins, then the
Wilson sites. We choose a bond n0 inside the Wilson chain, at a shell whose
energy is below the scale of the spins.

1. The head. DMRG solves the spins and the Wilson sites up to n0 plus m
   more. This is the present calculation, on a short chain. The left bond
   basis at n0 contains the spins and the Kondo flow from D down to the
   scale of n0. Its size is the smallest of three numbers: `maxm`, the
   4^m states of the m sites to its right (16 for m=2, 256 for m=4), and
   what the `cutoff` leaves. A state with a small weight in the ground
   state is cut, also when its energy is low. So m and `cutoff` are
   parameters of the method next to `maxm` and n0.
2. The hand-over. We write five things in that left basis: the block
   Hamiltonian H_n0 with its ground energy subtracted, the charge N and
   the S^z of the block, the two matrices of f_n0 that the next hopping
   needs, and the tip operators X^dag_sigma. N gives the fermion parity
   for the Jordan-Wigner sign of the next hopping, and N and S^z together
   give the (Q, 2Sz) labels of `nrgthermo`.
3. The tail. From here on it is NRG: add a site, diagonalize a matrix of
   size 4 maxm, keep the lowest states, carry f and X. The kept and
   discarded states of every shell are stored, so we have the complete
   basis of Anders and Schiller (cond-mat/0604517).
4. The spectrum. A backward pass gives the reduced density matrix of every
   shell, and the Lehmann sum over kept and discarded states gives the
   poles of A(w) (Weichselbaum and von Delft, cond-mat/0607497; Peters,
   Pruschke, Anders, cond-mat/0607494). `zaverage.broaden_poles` and the
   twist average are used as they are.

You can think of the head as a renormalized impurity: the spins, dressed by
the high-energy part of the bath, with `maxm` states. Schwarz et al.
(arXiv:1708.06315) use the same idea in the other direction, NRG first and
an MPS afterwards.

Above the scale of n0 the spectrum comes from TASK on the head. We expect
the Krylov space to converge there, since the head is short and its lowest
scale is high, but this is not measured: on `spin_chain_on_metal` two
orders of the sites differed by 9e-2 of the peak above |eV|=1 at one
`nkry`, and this is still open. Below the scale of n0 the spectrum comes
from the tail. In a window around that scale both are valid, and their
agreement there is a check that needs no outside reference. Nothing
conserves the weight across the cut, so the joined spectrum does not keep
the sum rule by construction, and the ladder measures it.

"Add sites until it converges" is then exact in the sense of NRG. The
energies of each shell, in units of the hopping of that shell, go to a
fixed point with period two in the number of sites. The tail stops when two
shells of the same parity agree to a relative tolerance. At T=0 the
spectrum below the last shell does not change any more. At T>0 the tail
stops at the shell of the temperature, so a low temperature comes with no
extra machinery, as long as T is below the scale of n0. This holds only if
the basis at n0 passes the gate below, since at T>0 the excited kept
states carry weight too.

A deep tail needs the symmetry labels. At B=0 `nrgthermo` sets the spectra
of the blocks (Q, 2Sz) and (Q, -2Sz) to their mean, and without this a
chain past about 46 sites at Lambda=4 gave a wrong T_K. The tail inherits
this fix through the labels of the hand-over.

Three errors have to be measured apart, each against its own reference.

- The basis error: do the `maxm` states at n0 contain the low eigenstates
  of H_n0? The reference is the exact NRG of milestone N1 on a few spins,
  where we can compute the weight of every NRG kept state inside the DMRG
  basis, against `maxm`. Saberi et al. find that the needed bond dimension
  is larger for the kept states than for the ground state alone. Their
  rule D' = 2^m sqrt(D) is for a chain with spin up and down unfolded, and
  does not transfer to our d=4 sites.
- The missing weight: X^dag changes the charge by one, and X^dag times the
  ground state can leave the basis. TASK already measures this as
  `missing_weight`, and `expand_bonds` adds the missing directions.
- The shell-diagonal error: NRG neglects the coupling between shells. The
  reference is TASK on the same exported state, never on a second solve,
  since two solves of one chain differ by 3e-3 of the peak.

There is also a depth floor. DMRG converges the total energy to an absolute
accuracy, and a Wilson hopping at shell 40 with Lambda=4 is 1e-12. So a bond
basis deep in the chain is not reliable at the scale of its own shell.
`CLAUDE.md` lists the same effect as open for bonds far below the largest
coefficient. This is the reason the tail is an NRG iteration and not more
DMRG sites: each NRG step works relative to its own ground energy, as the
`offset` of `nrgthermo` does. The plan measures the deepest n0 at which the
DMRG block spectrum still agrees with `nrgthermo` in relative terms.

If the basis error cannot be closed with `maxm` and `expand_bonds`, the
next step is an MPS that holds several states with one extra index, and
sweeps that lower the sum of their energies (Pizorn and Verstraete,
arXiv:1102.1401). Nothing in the code holds several states in one MPS
today, so this is new machinery. It enters only if the measurement asks
for it.

## 6. Geometry: the bath on one spin, or on every spin

With the bath on one spin there is one channel, and everything above
applies as written.

With a bath on every spin, the number of channels is the rank of the
coupling matrix of the spins to the metal (Eickhoff and Anders,
arXiv:2008.09034). For atoms a few lattice constants apart this is close to
ns. In exact NRG each channel multiplies the size of a shell by 4. In a
chain MPS the channels have to be interleaved, and the bond dimension grows
exponentially with their number. Interleaved Wilson chains exist for this
(Mitchell et al., arXiv:1308.1903; Stadler et al., arXiv:1602.02182). Our
estimate is that they carry us to about three channels; we did not read
this number off the two papers. Beyond that, the geometry that scales is one Wilson
chain hanging from each spin, a fork tensor network (Holzner et al.,
arXiv:0804.0550; Bauernfeind et al., arXiv:1612.05587). `pyitensor` has
only chains, so this is a separate project with its own plan.

## 7. Milestones

- N0. Spectral NRG with exact diagonalization for one spin 1/2
  (`nrgspectral.py`, on the internals of `nrgthermo`). It stores the
  rotation of every shell, keeps the kept and discarded states, carries
  X^dag_sigma and the spin, does the backward pass, and returns poles at
  T=0 and T>0. The start shell is an argument: energies, labels, f
  matrices and operator matrices. N1 and N3 reuse it through that argument.
- N1. Several spins in the start shell, exact: ns spins with Js and the
  bath on `site`. Measures `keep` against ns and Js/D, and does the 3-spin
  test of the shortcut of "Where NRG with exact diagonalization stops".
- N2. Shells read off a DMRG state. For each bond, the block Hamiltonian in
  the left basis, relative to its ground energy, and its spectrum. This is
  the flow diagram of an MPS. Measures the basis error and the depth floor.
  This is the gate for everything after it.
- N3. The hand-over and the tail: export H_n0, N, S^z, f_n0 and X^dag at
  bond n0 from "python" and v3 (the export of `taskdmrg.session_spectrum`
  is the model), then N0's iteration. The labels come from diagonalizing N
  and S^z of the block in the left basis. A first test without labels is
  fine on a shallow tail only, since a matrix of size 4 maxm is small, and
  the parity is needed in any case.
- N4. The spectrum: tail below the scale of n0, TASK on the head above it,
  the check in the overlap window, the twist average for
  `spin_chain_on_metal` (`zaverage` is tied to `kondo_chain` today), and
  the stop rule on the fixed point.
- N5. Only if N2 fails: several states in one MPS and variational sweeps.
- N6. Two and three channels on interleaved chains. The fork network is not
  part of this plan.
- N7. Public entry (`get_stm_spectrum(..., submode="NRG")` or its own
  function), a page in `docs/guide/`, a folder in `docs/new_methods/`,
  examples with plots.

New `pyitensor` code takes its arrays from `backend.xp()`, as the rule for
new methods asks. The tail is dense linear algebra on small matrices and
stays on the host.

## 8. Validation ladder

1. N0 without truncation on 3 to 5 Wilson sites: the poles are those of a
   full ED of `kondo_chain`, to roundoff, at T=0 and T>0.
2. J=0: the spectrum is U^2/2 times the density of states, as
   `zaverage_stm_spectrum` gives with `continuum_norm`.
3. The sum rule of A(w), which the complete basis keeps exactly in N0 and
   N1.
4. One spin, Kondo regime: the unitary value at w=0 that the `zaverage`
   docstring records as untested, and a width set by Wilson's T_K of
   `nrgthermo`. The tolerance is whatever the measurement supports at
   Lambda=4 with the twist average; we do not promise it in advance.
5. One spin in a field: the splitting of the resonance, against Costi.
6. N1 without truncation against a full ED of `spin_chain_on_metal`.
7. N3 against N1 on 3 to 5 spins, at T=0 and T>0, as a function of `maxm`,
   `cutoff`, n0 and the number m of head sites after n0.
8. N4: tail against TASK in the overlap window, on one exported state, with
   the convergence of TASK in `nkry` on the head. The integrated weight on
   each side of the cut from both routes, and the total weight of the
   joined spectrum against the sum rule, as measured numbers.
9. A chain beyond the exact route (10 to 15 spins): convergence in `maxm`,
   n0 and `keep`, the singlet against triplet ground states of
   `spin_chain_on_metal` (bath on an even or an odd spin), and the
   resonance present only where the ground state has a moment.

## 9. What the method does not give

- Sharp features at finite energy. The logarithmic mesh broadens a line at
  energy w by a fraction of w, so the spin excitation steps at eV near Js
  stay with TASK on the head. A mesh that is linear in a window and
  logarithmic below it (Schwarz et al.) is a later option.
- More than one channel before N6.
- High temperature for many spins. T has to stay below the scale of n0.
  Above it, the exact route of N1 is the only one, up to its limit in ns.
- A superconducting substrate, which is not in this plan. NRG in a BCS
  host exists (arXiv:1008.1734, the reference of `ysr.py`). Here the
  pairing breaks Q, so the (Q, 2Sz) blocks of `nrgthermo` become parity
  and Sz, and the present sector route of `ysr.py` stays for now.

## 10. Papers and codes

- Verstraete, Weichselbaum, Schollwoeck, Cirac, von Delft,
  cond-mat/0504305: NRG states are MPS, and a variational version.
- Saberi, Weichselbaum, von Delft, arXiv:0804.0193: the DMRG ground state
  of a Wilson chain contains the NRG states; block Hamiltonian per site,
  flow diagrams from an MPS. The model for N2.
- Pizorn and Verstraete, arXiv:1102.1401: NRG as a variational method for
  an MPS with a state index; cost D^3, as NRG. The model for N5.
- Anders and Schiller, cond-mat/0604517: the complete basis of discarded
  states. Weichselbaum and von Delft, cond-mat/0607497, and Peters,
  Pruschke, Anders, cond-mat/0607494: spectral functions that keep the sum
  rule. The model for N0.
- Costi, cond-mat/0004302: the T-matrix of the Kondo model by NRG, with a
  field. Reference curves for the ladder.
- Bulla, Costi, Pruschke, cond-mat/0701105: the NRG review.
- Weichselbaum, arXiv:1209.2062: NRG in tensor-network form.
- Gleis, Li, von Delft, arXiv:2207.13161: kept and discarded spaces of an
  MPS; arXiv:2207.14712: controlled bond expansion, already in
  `pyitensor/cbe.py`.
- Schwarz, Weymann, von Delft, Weichselbaum, arXiv:1708.06315: NRG states
  as one site of an MPS (the renormalized impurity), and the mesh that is
  logarithmic outside a window and linear inside.
- Eickhoff and Anders, arXiv:2008.09034: several impurities map to a number
  of channels equal to the rank of the coupling matrix; NRG up to 5 sites.
- Mitchell, Galpin, Wilson-Fletcher, Logan, Bulla, arXiv:1308.1903, and
  Stadler, Mitchell, von Delft, Weichselbaum, arXiv:1602.02182:
  interleaved Wilson chains.
- Holzner, Weichselbaum, von Delft, arXiv:0804.0550, and Bauernfeind et
  al., arXiv:1612.05587: one chain per channel in a star or fork network.
- Peters, arXiv:1103.5837: Lehmann spectra from the DMRG basis for three
  coupled Anderson impurities.
- Danu, Vojta, Grover, Assaad, arXiv:2204.00029, and arXiv:1903.08622:
  spin chains on a metal by quantum Monte Carlo, with a bath on every
  spin. Choi et al., arXiv:1904.09941: the review of the experiments.
- NRG Ljubljana (github.com/rokzitko/nrgljubljana): the reference NRG code,
  with complete-basis spectral functions. Its repository states GPL
  version 2 or later and a Zenodo record states AGPL version 3 or later,
  so the license file has to be read before any code is mirrored. dmrgpy
  is GPL-3.0-or-later. We use it as an outside cross-check for N0 at most.
- QSpace and MuNRG (Weichselbaum; Lee): we did not check their license or
  availability.

A spin cluster or a spin chain as the impurity of an NRG or of a DMRG on a
Wilson chain, with an STM spectrum as output, was not found in these
searches. This is a statement about the searches, not about the field.
