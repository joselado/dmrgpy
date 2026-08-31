# Sector-resolved dynamical correlators and spectral functions

Plan, 2026-08-31. Revised after an independent review (Sections 3, 4, 7,
8, 10 and 13 carry its corrections; the review's own findings are marked
"[review]").

**Status: Phases 1-3 implemented** (`src/dmrgpy/sectordc.py`,
`src/dmrgpy/multioperatortk/charge.py`,
`tests/test_sector_spectral_function.py`, three examples under
`examples/dynamical_correlator/sector_*`, documented in
`docs/user_guide.{md,tex}` and `docs/documentation.{md,tex}`). Phase 4
(Krylov enrichment of the target-sector basis) is not implemented and is
not needed for any use case so far: where the sum is truncated, the
`captured` diagnostic says so and `nex` fixes it. Two Phase-1 items were
deliberately dropped as the code took shape: the Loewdin
rediagonalization of Section 5's step 7 and the `rediagonalize=` kwarg
that would have exposed it. Neither was implemented, so nothing here
claims to have measured its effect -- it was deferred because the failure
it guards against (states within the target sector that are not
mutually orthogonal) is already voiced twice: by the per-state
`<H^2>-<H>^2` fluctuation warning, and by the captured weight exceeding
1, which is now warned about explicitly (Section 8). If either fires
routinely on a real problem, step 7 is the fix to implement.

## 1. The idea

`set_conserved_sector` (mpscpp3's QN block-sparse mode, and
`itensor_version="python"`'s charge-graded equivalent) confines a whole
calculation to one quantum-number sector. That makes a second, quite
different way of computing dynamical correlators available, one dmrgpy
does not have today:

    C_AB(w) = sum_n <0|A|n> <n|B|0> * delta_broadened(w - (E_n - E_0))

where the intermediate states |n> are *not* excited states of the full
Hilbert space (what `submode="EX"`/`dcex.py` computes) but the lowest
eigenstates **of the single sector that B|0> actually lives in**.

For an electron spectral function on a chain with N particles, B =
c_i^dag sends the N-particle ground state into the N+1 sector, and the
only states that can contribute are the N+1 eigenstates; the hole part
lives entirely in N-1. For a spin correlator on an Sz=0 ground state,
S^+ reaches Sz=+1 and S^- reaches Sz=-1 (2Sz units: +/-2), S^z stays in
Sz=0, and a two-magnon operator would reach Sz=+/-2.

Why this is worth having, next to KPM/CVM/TD/EX:

* Each intermediate state is a genuine, individually converged DMRG
  eigenstate of a *smaller* Hilbert space, and orthogonality to the
  other states in the sector is exact by construction rather than
  enforced by the overlap penalty. The excited-state search that
  `submode="EX"` relies on is much better conditioned inside a sector.
* No broadening/resolution tradeoff, no Chebyshev expansion order, no
  time window and no Fourier windowing artifact: the poles come out at
  their exact (converged) energies with exact weights. The only
  approximation is truncation of the sum at `nex` states, and that
  truncation is *measurable* (Section 8's sum rule).
* It resolves which quantum-number channel each feature belongs to,
  which none of the existing submodes do. A spin-1 chain's S(q,w) split
  into its Sz channels, or a Hubbard chain's spectral function split
  into upper/lower Hubbard band, come out separately.

The cost: it only ever gives the *low-energy* part of the response, up
to whatever the `nex` lowest states in the sector cover. It is the right
tool for a few sharp low-lying excitations and the wrong one for a broad
continuum, which is exactly the opposite tradeoff to KPM.

## 2. Literature check

The method is textbook rather than novel: targeting excited states in a
fixed quantum-number sector and summing the Lehmann representation over
them is how DMRG spectral functions were computed before correction
vector / KPM / time-evolution methods took over, and it is still used
where individual sharp excitations matter. Representative:

* Jeckelmann, *DMRG methods for momentum- and frequency-resolved
  dynamical correlation functions*, arXiv:0808.2620 -- review; Section
  on the Lehmann/spectral-weight route and its limits versus the
  correction vector.
* Peters, *Spectral functions for single- and multi-impurity models
  using DMRG*, arXiv:1103.5837 -- explicit "peak spectrum" from Lehmann
  weights of DMRG-targeted states.
* The selection-rule structure ("the ground state is in the total
  spin-0 sector, so S^z connects only to spin 0/1 excitations") is
  standard in the dynamical-structure-factor literature, e.g.
  arXiv:2501.13059 (J1-J2 spin-1 chain bound states).

No permissively licensed reference implementation to mirror was found:
ITensor(s).jl ships no sector-resolved Lehmann helper, and MPSKit's
excitation machinery is the infinite-system quasiparticle ansatz, a
different object. So this is written from scratch against dmrgpy's own
sector API, with the references above only fixing conventions.

## 3. Load-bearing facts this design rests on (all verified in-tree)

1. `Chain::promote_mps` (mpscpp3/chain_session.h:487, bindings.cc:144)
   and `pyitensor/chain.py::promote_mps` rebase a wavefunction onto
   **`dense_sites_` / `self.dense_sites`** -- the chain's *own original*
   site set, fixed at construction (`chain_session.h:290`,
   `pyitensor/chain.py:167`) and never rebuilt by
   `set_conserved_sector`, which only replaces `sites_`. So states
   promoted out of *different* sectors of the same session land on the
   *same* dense indices and can be contracted against each other and
   against operators built later. This is the single fact the whole
   method depends on.
2. **[review] That fact is already regression-tested**, on both
   backends: `tests/test_sector_promotion.py:195
   ::test_states_promoted_from_different_sectors_are_comparable` (and
   its `test_sector_promotion_python.py` twin) solves Nf and Nf-1 in
   sequence on one chain, promotes both, and asserts the photoemission
   weight `|<gs_{N-1}|c_i|gs_N>|^2` against ED site by site. That is
   precisely this method's core matrix element. Bond indices are a
   non-issue for the same reason: independently produced MPS carry
   distinct link indices, and `innerC(wf1,A,wf2)`
   (`chain_session.h:1085-1094`) contracts each against its own dagger.
3. **[review] Jordan-Wigner is not a backend convention question here.**
   dmrgpy applies its own JW strings in Python before any backend sees
   the terms (`multioperator.py:331`, anchored at site 1), identically
   for A and B and for every backend. A different string anchoring would
   contribute a total-parity factor to each of `<0|A|n>` and `<n|B|0>`,
   which cancels in their product between states of definite N.
4. `Many_Body_Chain.promote_mps(wf)` (manybodychain.py:491) exposes the
   conversion per wavefunction, on both backends, with a clear error on
   a pre-sector compiled extension.
5. Clearing the sector (`set_conserved_sector()`, no args) sets
   `sites_ = dense_sites_` -- so after clearing, previously promoted
   handles are consistent with everything built from then on. The
   charge-changing operator B can only be built *after* that point
   (`sector_terms` rejects it inside a sector), which fixes the step
   order in Section 5.
6. Excited states inside a sector work on both backends:
   `Chain::excited_states` starts from `default_mps()`, which routes to
   `sector_mps()` when a sector is set (chain_session.h:4994), and
   pyitensor's `excited_states` runs `_check_in_sector` per state.
7. Index identity is per session. A fresh `Chain` mints fresh Indices,
   so MPS handles are **not** portable between chain objects --
   everything below must happen inside one session.
8. `Chain::op_charge` (chain_session.h:4585) and
   `pyitensor/sector.py::op_charge` both infer a single-site operator's
   charge from its dense matrix elements, returning "no definite charge"
   for e.g. `Sx`. The Python one is pure Python and usable from any
   backend, which is what the dispatch layer needs.

## 4. Core API

A new submode, not a new top-level method, so the whole existing
cross-check surface (`mode="ED"` reference, KPM/TD/CVM comparison,
`examples/dynamical_correlator/*_VS_ED`) applies unchanged:

```python
(x, y) = fc.get_dynamical_correlator(submode="SECTOR",
                                     name=(fc.C[0], fc.Cdag[0]),
                                     nex=20, delta=0.05,
                                     es=np.linspace(-2., 8., 800))
```

**Return convention, pinned. [review]** The existing submodes do not
agree with each other on dtype: `edtk/dynamics.py`'s ED reference
returns a real array (`-Im[...]/(2 pi)`), while `dcex.py:112-114` and
the KPM path return a complex one. `SECTOR` returns **complex, exactly
`dcex`'s `1j*(adv-ret)/(2 pi)`** -- it is the DMRG-side sibling of
`submode="EX"` and must be interchangeable with it. Its real part is the
ED convention. This matters beyond dtype when `A != B^dag` (the weight
`<0|A|n><n|B|0>` is then genuinely complex), so the docs state it and
the ED cross-check compares against `.real`.

New kwargs, all optional:

* `nex=20` -- states per target sector, **capped at the exact dimension
  of the target sector [review]**, computed with the same reachable-
  charge DP that `sector_state_plan` already runs. Without the cap the
  default overshoots small sectors (the N+/-1 sectors of this plan's own
  n=8 validation chains, or an Sz=+2 sector of a short spin chain): the
  overlap-penalty search then hunts states that do not exist, burns most
  of the runtime, floods the user with fluctuation warnings, and gets
  quietly rescued by the Loewdin cut -- a plausible-looking spectrum
  produced the slow, alarming way. Phase 1, not later.
* `sector=None` -- the reference sector Q of the ground state. Default:
  the chain's own `conserved_sector` if set; otherwise inferred by
  measuring the conserved charges on an unconstrained ground state and
  rounding, with a warning if the measurement is not close to an
  integer (a genuinely charge-superposed ground state has no sector
  interpretation, and this method cannot be applied to it).
* `conserve=None` -- which quantum numbers to conserve. Default
  **[review]**: every quantity the site types offer *that the
  Hamiltonian actually conserves*, tested with Section 6's own helper
  (`multioperator_charge(H) == 0` per quantity). Taking everything on
  offer instead fails correctly but late, deep in the Hamiltonian send,
  on e.g. a Hubbard chain with spin-orbit coupling.
* `return_poles=False` -- also return the raw `(energies, weights)`
  Lehmann poles instead of only the broadened curve. This is what makes
  the method more than another broadened spectrum. **[review]** Note in
  the docs that inside a degenerate multiplet the *individual* pole
  weights are basis-dependent and therefore run-dependent (only the
  multiplet sum is invariant), and that truncating `nex` in the middle
  of a multiplet drops a run-dependent share of its weight.

Error cases, all raised with the offending operator named:

* `A` or `B` has no definite charge (`Sx`, `C+Cdag`): `ValueError`
  pointing at Phase 2 (decomposition) below. This message has to be
  excellent -- `name=(Sx,Sx)` is the single most common spin-correlator
  call in the examples tree, and it errors until Phase 3.
* `charge(A) + charge(B) != 0`: the correlator vanishes identically by
  symmetry. Raise rather than return zeros -- silently returning a zero
  spectrum for a typo'd operator pair is exactly the kind of dispatch
  hole `docs/audit_2026_08_hole_hunt.md` catalogues.

**Gating must be explicit, inside the new branch. [review]** The earlier
draft said this could reuse `mode.py`'s sector messages; it cannot. The
sector lives on an internal clone (Section 5), so the *user's* chain has
no sector set and `mode.py:27`'s `if sector:` never fires. Therefore:
add `"SECTOR"` to `dynamics.py`'s `SUBMODES` tuple, and in the SECTOR
branch itself reject `mode="ED"` (it would otherwise reach
`edtk/dynamics.py` and die with a generic "no ED implementation"),
`itensor_version=2` (`dynamics.py`'s DMRG branch admits 2, 3 and
"python"), `julia_live` (which forwards the submode to `dynamicsjl`),
`ns<3` on v3, and a non-Hermitian Hamiltonian -- verifying, not
assuming, that the existing non-Hermitian raise-list already covers the
last of these. No ED fallback anywhere: a fallback would answer with the
*global* excited states, a different calculation.

## 5. Pipeline

All of it inside **one session on a deepcopy of the user's chain** (the
`gap.py` pattern): the pipeline ends with the sector cleared, which must
not be a side effect on the caller's object, and `__deepcopy__` already
rebuilds a fresh session and re-applies the sector.

1. Determine `names` (conserved quantities) and the reference sector Q
   (Section 4).
2. Determine `dq = charge(B)` over `names` -- a *tuple*, one integer per
   conserved quantity (Section 6). Target sector Q' = Q + dq. Check
   `charge(A) == -dq`.
3. `clone.set_conserved_sector(**Q)`; `gs_energy()` -> `E0`, `wf0`.
   `wf0_dense = clone.promote_mps(wf0)`.
4. `clone.set_conserved_sector(**Q')`; `get_excited_states(n=nex,
   purify=False)` -> `E_n`, `wfs`. `wfs_dense = [promote_mps(w) for w in
   wfs]`. (The existing `_warn_unconverged_excited_states` fluctuation
   warning fires here for free.)
5. `clone.set_conserved_sector()` -- back to dense. Every handle from
   steps 3-4 is still valid (fact 5), and only now can the
   charge-changing A/B be built at all.
6. Matrix elements with `clone.aMb(wf_n, B, wf0_dense)` and
   `clone.aMb(wf0_dense, A, wf_n)` -- `aMb` sandwiches directly, with no
   intermediate compressed MPS, which is why `dcex.py` switched to it.
7. Optional refinement, reusing `dcex.py`'s machinery verbatim: build
   `S`, `H` in the promoted basis, Loewdin-orthogonalize (drop
   directions below `1e-8*max`), rediagonalize. The basis is
   near-linearly-dependent for the same reasons as in `dcex`, and H is
   block diagonal across sectors so this stays inside the target sector.
   **[review]** Send H as a prebuilt `StaticOperator`/MPO here:
   `Chain::overlap_aMb` (chain_session.h:1090) rebuilds the operator MPO
   from terms on *every* call, so the nex^2 Hamiltonian sandwiches would
   otherwise rebuild H's MPO nex^2 times. `mpsalgebra.py:71` already
   special-cases `StaticOperator`, so this is a call-site choice, not
   new machinery.
8. Lehmann sum in the convention pinned in Section 4.

Two target sectors (`Q+dq` and `Q-dq`, i.e. particle and hole) are only
needed by the Layer-2 spectral function of Section 7 -- one
`get_dynamical_correlator(name=(A,B))` call has one target sector, fixed
by `B`, exactly as the existing correlator convention implies.

**Global phases need no fixing. [review]** Each converged `|n>` carries
an arbitrary, genuinely run-dependent global phase (v3 starts from
`sector_mps`, a random in-sector state). It cancels exactly in the
product `<0|A|n><n|B|0>` -- the bra contributes `e^{-i theta_n}`, the ket
`e^{+i theta_n}` -- *provided both factors are computed from the same
promoted handle*, which steps 3-6 guarantee. Same for `|0>`'s own phase.
Worth stating in the docstring, since it is the first thing a reader
worries about.

## 6. Charge inference

A backend-agnostic Python helper (new `multioperatortk/charge.py`, or
extend `pyitensor/sector.py` and import from it -- pyitensor is pure
Python and always importable regardless of `itensor_version`):

```python
def multioperator_charge(sites, MO, names):  # -> tuple[int] or None
```

For each term of the `MultiOperator`, sum `op_charge` over its factors;
the operator has a definite charge only if every term with a nonzero
coefficient agrees. Returns `None` otherwise. Built on
`pyitensor.sites`' `state_charge`, so it needs only the chain's site-type
codes, not a live session. Also used for the `conserve=` default
(Section 4) by applying it to the Hamiltonian.

Multi-quantity from day one: on a spinful chain `Cdagup` carries
`(Nf=+1, Sz=+1)`, and the Sz component is what makes the spin-resolved
spectral function of Section 7 fall out for free.

## 7. Layer 2: the physics-facing spectral functions

Two convenience methods that call Layer 1 twice and assemble the
standard object, since "give me A(w)" should not require the user to
know that the hole part is a separate sector solve:

**`get_spectral_function(i, j=None, spin=None, ...)`** on fermionic
chains (spinless and spinful):

    A_ij(w) = sum_n <0|c_i|n^{N+1}><n^{N+1}|c_j^dag|0> L(w - (E_n^{N+1}-E0))
            + sum_m <0|c_j^dag|m^{N-1}><m^{N-1}|c_i|0> L(w + (E_m^{N-1}-E0))

the greater part at `w = E_n^{N+1}-E0` and the lesser part mirrored to
`w = -(E_m^{N-1}-E0)`, on one frequency axis.

**[review] The chemical potential is not cosmetic here.** In these raw
energies the particle poles sit at `w >= E_0^{N+1}-E_0^N` and the hole
poles at `w <= -(E_0^{N-1}-E_0^N)`; both can land on the *same* side of
zero for any Hamiltonian with an explicit chemical-potential term or
attractive interactions, so "particle at positive w, hole at negative w"
is only true when mu is near zero. The method therefore shifts by
`mu = (E_0^{N+1} - E_0^{N-1})/2` **by default** (`shift="mu"`, with
`shift=None` for the raw axis), which is the convention the formula
above is usually written in, and always returns the unshifted pole
energies alongside.

**`get_spin_spectral_function(i, j=None, ...)`** on spin chains and
spinful fermion chains: transverse channel from the two sectors
`Sz +/- 2` (2Sz units) via `(S^-, S^+)` and `(S^+, S^-)`, longitudinal
channel from `(S^z, S^z)` inside `Sz = 0`; returns the channels
separately and their sum, because the whole point of the method is that
the channel decomposition is available. **[review]** In the longitudinal
channel the target sector *is* the reference sector, so `n=0` is the
ground state itself and the sum contains an elastic `w=0` pole of weight
`<0|Sz_i|0><0|Sz_j|0>`. That is correct (the ED reference contains it
too, its state sum running over everything), but users expecting the
*connected* structure factor will be surprised, so it is documented and
a `connected=True` switch subtracts it (default `False`, i.e. the raw
Lehmann sum that matches the ED reference term for term).

Both accept a list of sites (or all sites) so `S(q,w)` / `A(k,w)` follow
from one Fourier transform of the site-resolved poles, which is much
cheaper here than with KPM: the sector solves are shared across all
`(i,j)` pairs -- one set of `nex` states per sector serves the entire
matrix of site pairs. That sharing is a real, method-specific win and
should be in the design from the start (cache keyed on
`(sector, nex, hamiltonian)`, like `dcex.get_cached_excited_states`).

## 8. Diagnostics (not optional)

Truncating the Lehmann sum at `nex` states is the method's only
approximation, and it is directly measurable:

    captured = sum_n |<n|B|0>|^2 / <0|B^dag B|0>

**[review]** The denominator must be evaluated as
`clone.vev(B.get_dagger()*B, wf=wf0_dense)` -- explicitly on the
promoted ground state. A bare `vev()` at that point in the pipeline
would find `wf0=None` and `computed_gs=False` (step 5's
`set_conserved_sector()` calls `restart()`) and silently *re-solve
unconstrained*, i.e. normalize the sum rule against the wrong state.

The same number above 1 is the one silent failure this method has: more
weight than B|gs> holds means the target sector's states are not as
mutually orthogonal as the sum assumes, and some of it is being counted
twice. Within a sector they are still separated by the overlap-penalty
search, which can return two overlapping states without saying so, so
this is warned about explicitly as well.

Return `captured`, and warn when it drops below a threshold (say 0.9) --
"this spectrum accounts for 62% of the spectral weight of B|0>, increase
nex" is the difference between a usable method and a silently wrong one.
Same for the `<H^2>-<H>^2` fluctuation warning already emitted per state.

## 9. Backend gating

`itensor_version in (3, "python")` and DMRG only, `ns >= 3` on v3,
Hermitian only -- enforced explicitly in the new branch, for the reason
given at the end of Section 4. `julia_live` and v2 raise
`NotImplementedError` naming the two supported backends.

## 10. Validation

* **Exact reference**: `mode="ED"`, `submode="ED"` on n=8-12 chains --
  the full Lehmann sum. The `SECTOR` result must match it peak for peak
  in the region covered by the captured weight, comparing `.real`
  against ED's real convention (Section 4). **[review]** The validation
  chains must have a **unique** ground state (even-length
  antiferromagnetic chains): `dynamical_correlator_ED` averages over a
  near-degenerate ground manifold with equal `1/nex` weights below its
  `dex` cutoff, while SECTOR uses one DMRG ground state, so on a
  degenerate chain the two disagree by construction rather than by bug.
* **[review] The assembled A(w) has no single-call ED reference.** ED
  `name=(C_i,Cdag_j)` gives only the particle poles; the hole part is a
  separate `(Cdag_j,C_i)` call whose axis must be mirrored before
  summing. The validation script performs the same two-call assembly on
  the ED side -- otherwise it "validates" the particle half against the
  whole.
* **Sector-restricted ED**: reuse `tests/test_sector_conservation.py`'s
  pattern (*restrict* the ED Hamiltonian to the sector, never filter
  full-spectrum eigenvectors by `<N>` -- degenerate eigenvalues come
  back as arbitrary superpositions) to check the per-sector energies and
  matrix elements individually, not just the assembled curve.
* **Free fermions**: `tests/_free_fermion_reference.py` gives the exact
  single-particle spectrum of a tight-binding chain; its spectral
  function is a set of delta functions at the single-particle levels
  with known weights. An analytic test independent of ED.
* **Cross-submode**: agreement with `submode="KPM"` and `submode="TD"`
  on a moderate chain, within KPM's own broadening -- the check that
  matters for a user deciding between methods.
* **Both backends**: v3 vs `"python"` on the same small chain, the
  standard pattern for anything sector-related.
* **Hubbard**: Mott gap visible in `get_spectral_function`, and the
  gap value matching `E_0^{N+1} + E_0^{N-1} - 2 E_0^N` computed
  directly. A physics-level test the numerics cannot fake.

**[review] Write this one first, before any C++ is touched**: an n=6-8
tight-binding chain on `itensor_version="python"` with `nex` = the exact
sector dimension, comparing pole energies *and weights* against
sector-restricted ED. It exercises promotion, cross-sector `aMb`, the JW
signs, the Lehmann assembly and the `nex` cap in seconds. Its second
assertion must use a **non-self-adjoint pair** -- `(S^+, S^-)` on a
short spin chain -- because `dcex.py:61`'s `A = name[0].get_dagger()`
placement is invisible for pairs like `(C,Cdag)`, and a misplaced
conjugation there is exactly the silent convention slip this family of
code has been bitten by before.

Tests: `tests/test_sector_spectral_function.py`. Examples (each ending
in a plot, per CLAUDE.md): a Hubbard chain A(w) with the Mott gap; a
Heisenberg chain S(q,w) with the des Cloizeaux-Pearson lower edge; a
`SECTOR` vs `KPM` vs `ED` overlay; a `nex` convergence sweep showing the
captured-weight diagnostic rising toward 1.

## 11. Documentation

New subsection in `docs/user_guide.{md,tex}` under the dynamical
correlator section (the physics, the sum rule, the return convention,
when to prefer this over KPM), and in `docs/documentation.{md,tex}` for
the dispatch and the promote/clear pipeline (this adds a dispatch
branch, and `documentation.md` 4.10 is explicitly about getting those
right). `.tex` recompiled with `pdflatex`.

## 12. Phases

1. **Core**: charge inference, `submode="SECTOR"`, the `nex` sector-
   dimension cap, diagnostics, explicit gating, the free-fermion
   pure-Python test first, then ED and v3 tests, one example, docs.
2. **Layer 2**: `get_spectral_function` / `get_spin_spectral_function`,
   shared sector-solve cache, site-resolved sweeps and the k-space
   Fourier transform, Hubbard/Heisenberg examples.
3. **Charge-indefinite operators**: decompose `Sx` -> `(S^+ + S^-)/2`
   etc. into charge-homogeneous components, pair them up (only
   `dq_A = -dq_B` survives), sum contributions over the two sectors.
   Makes `name=(Sx,Sx)` -- the most common spin correlator in the
   examples tree -- work without the user rewriting it in ladder
   operators.
4. **Krylov enrichment (optional, do not let it creep earlier)**: extend
   the target-sector basis with `B|0>, H B|0>, H^2 B|0>, ...` and
   rediagonalize, so the sum rule is exhausted by construction rather
   than by luck. Must happen entirely in dense mode: `B|0>` can only be
   built after promotion, and there is no demote-to-QN path. Cheap on
   top of Section 5's step 7, which already solves a generalized
   eigenproblem in a non-orthogonal basis.

## 13. Known risks

* **Sum-rule truncation is fundamental**, not a bug: for a gapless
  continuum the low-lying states carry a small fraction of the weight
  and `nex` would have to be huge. The diagnostic makes this visible;
  the docs must say plainly that this method is for sharp low-lying
  excitations.
* **[review] `nex` past the target sector's dimension** is the single
  most likely practical failure, because it fires on this plan's own
  validation chains and does not crash: it produces a slow,
  warning-spewing run whose spectrum still looks plausible, and reads as
  "the method doesn't converge". Hence the Phase-1 cap in Section 4.
* **Excited-state quality inside a sector** still uses the
  overlap-penalty method, whose failure mode (a spurious stationary
  point) is documented in `pyitensor/dmrg.py`. Better conditioned in a
  sector, not immune. The fluctuation warning is the guard.
* **Cross-sector energy differences** are differences of separately
  converged variational energies. Both sides must be converged to the
  same standard, and the physical quantity is a small difference of two
  large numbers -- the tests must check the *gap*, not just each energy.
* **Degenerate multiplets** make individual pole weights run-dependent
  (Section 4, `return_poles`); only multiplet sums are meaningful.
* **v3 QN mode is off by default** for a reason (`get_sites.h`); the
  sector path is the opt-in one, and its per-block bookkeeping only pays
  off with size. For the small chains most examples use, this method is
  not necessarily faster than KPM -- it is *sharper*, which is the point.
