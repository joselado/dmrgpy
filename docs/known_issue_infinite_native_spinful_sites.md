# RESOLVED: `Infinite_Many_Body_Chain` with native spinful (site_type=1) sites

**Status: fixed on 2026-08-22.** Both halves of this report were real, and both have
been addressed; the diagnosis in the "What was tried" section below turned out to be
partly wrong about the *causes*, so read this header before the original text. The
report is kept in full because its reproduction case is now a regression test
(`tests/test_infinite_chain.py::test_interacting_spinful_cell_backends_agree`) and its
measurements are the before-numbers for the performance work.

## What was actually wrong, and what changed

**1. `itensor_version=3` -- the term-count error.** The report guessed that the C++
term-inspector counts operator *factors* rather than distinct site *indices*. It does
not: `idmrg_classify_terms` groups factors by site correctly, and
`infinitechain.py::_canonicalize_hamiltonian` already enforced the documented
"1- and 2-site terms" contract on the raw terms. The extra sites were created
*afterwards*, by `infinitechain.py` itself: it passed the C++ backend
`MultiOperator.to_terms()`, i.e. the **finite-chain** Jordan-Wigner form, whose strings
run from site 1 of the chain and put an explicit `F` factor on every intervening site.
An infinite chain has no site 1. Dumping the failing Hamiltonian's own terms showed all
16 inter-cell terms in forms like
`(-0.258, [('Adagup',1), ('F',2), ('F',1), ('Aup',3)])` -- three or four sites each.

This was also worse than an over-strict check: `idmrg_classify_terms` hardcoded
`carry_ferm=false` and `idmrg_build_row` always propagated `Id` instead of `F`, both
documented in their own comments as "fermionic terms not supported yet ... silently
wrong for a fermionic one". So merely relaxing the site count would have returned wrong
numbers. What the report saw as "works for the finite chain but not the infinite one"
was really "works whenever every fermionic term connects *adjacent* sites" -- the
global string then collapses to exactly the right endpoint-`F` composition with nothing
in between, which is the only case the spin-only test suite could ever exercise.

Fixed by giving `MultiOperator.to_terms()` a `jordan_wigner_transform=False` opt-out,
using it at both `infinitechain.py` v3 call sites, and porting
`pyitensor/idmrg.py::_term_site_matrices`' local, translation-invariant string threading
into `mpscpp3/chain_session.h::idmrg_classify_terms` (fermionic reordering sign,
endpoint `F`, `carry_ferm` wired through `idmrg_build_row`, odd-parity rejection). One
implementation covers both `gs_method="idmrg"` and `gs_method="vumps"`, which share that
classifier.

**2. `itensor_version="python"` -- the wall-clock.** Profiling (which the report asked
for) found nothing site-type-specific. Native spinful sites only make two general costs
bite harder, via `d_g = 4*4 = 16` instead of 4:

- `np.einsum` without `optimize=` in the innermost helpers (`_apply_op_ket`,
  `_op_transfer_matrix`, `_cap_right`, `_cap_left`): numpy's `c_einsum` never dispatches
  to BLAS. One D=8 macro-iteration spent 7.2s of 12.0s there across ~102k calls. Now
  written as `matmul`/`tensordot`.
- Two `O(chi^6)` dense steps per VUMPS iteration -- composing the unit-cell transfer
  matrix and `np.linalg.eig`-ing it for the dominant fixed point, and building the
  environment solve densely by applying its action `chi^2` times. At `chi=12` the `eig`
  alone was 14.7s of a 20s run. Now ARPACK/GMRES on matvecs, with the dense route kept
  as an exact fallback.
- A `D`-ramp that stepped one integer at a time, making the driver `O(D^4)`: 32 complete
  VUMPS solves to reach `D=8`, ~100 to reach `D=30`. Now doubling (`1,2,4,...,D`),
  mirrored in the C++ port.

Measured on the report's own kind of case: the gapped dimerized spinless chain at D=12
went from 19.9s to 1.0s (~19x), and VUMPS on the report's D=8 spinful case from 12.8s to
1.6s. The C++ VUMPS got the geometric ramp but not the Krylov fixed point (no ARPACK
there), so it now scales worse in `D` than the Python one -- see
`docs/documentation.md`.

**3. Coverage.** `tests/test_infinite_chain.py` grew end-to-end fermionic tests on both
backends and both solvers, anchored on the exact free-fermion band integral rather than
golden values, plus the `site_type=1` case and this report's own interacting cell.

## Where it stands now

- The report's exact Hamiltonian runs on both backends. At `maxm=8, maxiter=320` the two
  agree to 8e-7 (`python` -0.15885471, `v3` -0.15885551); `v3`+`gs_method="idmrg"` at the
  report's own `maxm=30, maxiter=60` takes ~250s and gives -0.16615.
- On *free* fermions, where an exact answer exists, every backend/solver combination
  reproduces the band integral to ~4e-10 (`examples/idmrg/fermionic_infinite_chain`).
- The interacting model itself converges slowly: its energy still wanders at the 1e-3
  level between `maxiter` settings, and the two backends' independent random starts land
  in slightly different places at moderate `maxm`. That is a property of this Kondo-lattice-like
  model, not of the fixed machinery -- but it means the workaround adopted below (finite
  chains) remains a reasonable choice for quantitative work on it.
- **Cross-site fermionic correlators now work too** (fixed immediately after the above, same
  day). `<C^dag_i C_j>` with `j != i` needs the Jordan-Wigner string on every site strictly
  between the two operators; no correlator path had one, which was the same failure class as
  the Hamiltonian bug, one layer up and in a separate code path. All three paths
  (`"python"`+idmrg, `"python"`+vumps, v3+vumps) now thread it, taking their endpoint
  matrices and their "is a string open" flag from the same helper that builds the
  Hamiltonian's own 2-site terms. Verified against the exact free-fermion one-body density
  matrix over r=0..8 including the alternating sign structure: 4.8e-08 (python/idmrg),
  1.5e-06 (python/vumps), 1.9e-06 (v3/vumps). A pair of *odd total fermion parity* (e.g.
  `Cdag` with `N`) still raises -- its string can never close, and the quantity vanishes
  identically in any parity-conserving state.

---

# Original report (2026-08-2x)

## Motivation / what was being attempted

A Kondo-lattice-like chain: one itinerant ("c") and one correlated ("f") spinful fermion
orbital per unit cell, with Hubbard `U1` on the f-orbital, a Heisenberg-like Kondo exchange
`J`, and a hybridization-like `tc` term (`(N_f - 1/2) * (spin-current between c and f)`).

`Spinful_Fermionic_Chain` (finite chains) represents each spinful orbital as **two**
interleaved spinless-fermion sites (up/down), so the physical 2-orbital (c,f) unit cell
needs **4** tensor-network sites. `Infinite_Many_Body_Chain` caps unit cells at `n_uc<=2`
(documented, `infinitechain.py`'s own constructor comment: "unit cells with more than 2
sites ... not supported yet ... needs a genuine redesign").

The natural workaround: use **native spinful sites** (`site_type=1`, `ElectronSite`/
`HubbardSite`, dim=4 per site — see `Spinful_Fermionic_Chain_Native` for the finite-chain
analogue) instead. That needs only **1** tensor-network site per orbital, so the (c,f) unit
cell fits exactly into `n_uc=2`. `pyitensor/sites/siteset.py`'s `TYPE_CODE_TO_SITE` maps
`site_type=1` to `ElectronSite` generically, and it's used directly by both `pyitensor/chain.py`
(finite) and `pyitensor/idmrg.py`/`vumps.py` (infinite) via the same `SiteX` class — so nothing
in the infinite framework's own site-construction code path rejects it. But nobody has actually
exercised this combination: `tests/test_infinite_chain.py`'s fermionic coverage (the
`_FERMION_SITE_CODE = 0` block, ~line 1466) only tests **spinless** fermion sites in the
infinite framework, never `site_type=1`.

## What was tried, and what happened

### 1. Free-fermion sanity check, `itensor_version="python"` (default), `gs_method="vumps"`

A 2-site (n_uc=2) native-spinful infinite chain, free-fermion hopping only (no interaction),
no chemical potential — i.e. the simplest possible case. **Hung**: killed after 24 CPU-minutes
with zero output (not even the pre-Hamiltonian numpy reference calculation had printed, though
that's instant — stdout block-buffering hid this, a separate minor confound). In hindsight this
specific test case likely landed at a gapless/metallic filling (no `mu` term to pin a definite
insulating filling), and VUMPS is documented elsewhere in this codebase to converge very slowly
for gapless cases (`examples/idmrg/local_excitation_gap/main.py`'s own docstring references
this). So this particular run may not indict the site-type=1 path specifically — see run 2.

### 2. Full interacting model, at a point known to be gapped, `itensor_version="python"`, `gs_method="vumps"`

Same 2-site chain, now with the actual `U1`/`tc`/`J`/`mu` interaction terms, at parameters
independently confirmed gapped via finite-chain DMRG (`J=0.2, tc=0.2, mu=0.1`, gap~0.04-0.05 in
the finite chain). `maxm=30`. Even with `maxiter` capped as low as **2**, did not finish within
120s (`time timeout 120 python3 ...` — exit 124). A finite chain of comparable physical content
(`Spinful_Fermionic_Chain`, N=10 unit cells = 40 sites, `maxm=40, nsweeps=30`) computes a full
ground state *and* first excited state in ~163-300s. So a single micro-iteration on a 2-site
*infinite* native-spinful chain is far more expensive than an entire finite calculation ~10x
its size — points at a real per-iteration cost problem specific to this site-type/algorithm
combination on the `"python"` backend, not merely "needs more iterations to converge."

### 3. Same interacting model, `itensor_version=3` (compiled C++ backend)

Same Hamiltonian, `ic = infinitechain.Infinite_Many_Body_Chain([1,1], itensor_version=3)`.
Failed **immediately** (~0.1s):
```
From line 3932, file chain_session.h
Chain::idmrg_ground_state: a term touches more than 2 distinct sites -- only 1- and 2-site terms are supported
```
Note this fired even though `gs_method` was left at `"vumps"` (the default) — the error text
names `idmrg_ground_state`, suggesting v3's vumps path still routes through (or shares a
term-count check with) the idmrg ground-state code. The likely cause: interaction terms built
as products like `(N_f - 1/2) * (Cdag_f * C_c)` are 3-operator products that only touch **2**
distinct site *indices* (both `N_f` and one `C`/`Cdag` factor live on the f site) — but the C++
term-inspector apparently counts something closer to *operator factors* (or fails to collapse
same-site factors into one local operator matrix before counting), rather than *distinct site
indices*, which is what the class's own Python-level docs promise support for (1- and 2-site
terms, meaning terms touching only C, or C and one adjacent cell). The finite-chain path
(`Spinful_Fermionic_Chain`, same style of 3-operator-product term) has no trouble with this
exact term shape, which points at the *infinite*-chain C++ term-canonicalizer specifically, not
`chain_session.h`'s general AutoMPO handling.

## Suggested reproduction (self-contained, needs `DMRGROOT` set)

```python
import os, sys
sys.path.append(os.environ["DMRGROOT"])
import numpy as np
from dmrgpy import infinitechain

ta, tb, alpha = -0.258, 0.011, -0.0545j   # free-fermion part, any reasonable values
U1, tc, J, mu = 0.513, 0.2, 0.2, 0.1       # interacting part

ic = infinitechain.Infinite_Many_Body_Chain([1, 1], itensor_version=3)  # or "python"
ic.gs_method = "vumps"
ic.maxm = 30
ic.maxiter = 60

CupC=[ic.get_operator("Cup",i,"C") for i in range(2)]; CdagupC=[ic.get_operator("Cdagup",i,"C") for i in range(2)]
CdnC=[ic.get_operator("Cdn",i,"C") for i in range(2)]; CdagdnC=[ic.get_operator("Cdagdn",i,"C") for i in range(2)]
CupR=[ic.get_operator("Cup",i,"R") for i in range(2)]; CdagupR=[ic.get_operator("Cdagup",i,"R") for i in range(2)]
CdnR=[ic.get_operator("Cdn",i,"R") for i in range(2)]; CdagdnR=[ic.get_operator("Cdagdn",i,"R") for i in range(2)]
NupC=[ic.get_operator("Nup",i,"C") for i in range(2)]; NdnC=[ic.get_operator("Ndn",i,"C") for i in range(2)]
SxC=[0.5*CdagupC[i]*CdnC[i]+0.5*CdagdnC[i]*CupC[i] for i in range(2)]
SyC=[-0.5j*CdagupC[i]*CdnC[i]+0.5j*CdagdnC[i]*CupC[i] for i in range(2)]
SzC=[0.5*NupC[i]-0.5*NdnC[i] for i in range(2)]

H = 0
for flavor in ["up","dn"]:
    Cdc = CdagupC if flavor=="up" else CdagdnC; Cc = CupC if flavor=="up" else CdnC
    Cdr = CdagupR if flavor=="up" else CdagdnR; Cr = CupR if flavor=="up" else CdnR
    H = H + ta*Cdc[0]*Cr[0] + ta*Cdr[0]*Cc[0]
    H = H + tb*Cdc[1]*Cr[1] + tb*Cdr[1]*Cc[1]
    H = H + alpha*Cdc[0]*Cr[1] + np.conj(alpha)*Cdr[1]*Cc[0]
    H = H + alpha*Cdc[1]*Cr[0] + np.conj(alpha)*Cdr[0]*Cc[1]
c, f = 0, 1
Tup = -1j*(CdagupC[f]*CupC[c] - CdagupC[c]*CupC[f])
Tdn = -1j*(CdagdnC[f]*CdnC[c] - CdagdnC[c]*CdnC[f])
H = H + U1/2*(NupC[f]-0.5)*(NdnC[f]-0.5)
H = H - tc/2*(NupC[f]-0.5)*Tdn + tc/2*(NdnC[f]-0.5)*Tup
H = H + J*(SxC[f]*SxC[c] + SyC[f]*SyC[c] + SzC[f]*SzC[c])
H = H + mu*NupC[c] + mu*NdnC[c]

ic.set_hamiltonian(H)
e0 = ic.gs_energy()   # itensor_version=3: fails instantly with the term-count error above
                        # itensor_version="python": does not return within 120s even at maxiter=2
```

## What would make this usable

1. `itensor_version=3`: fix (or at least correctly document) the term-touches-2-sites check
   for the infinite growth loop so multi-operator-product terms confined to 2 site *indices*
   (regardless of factor count) pass, matching the finite-chain behavior and the class's own
   documented 1-/2-site-term contract.
2. `itensor_version="python"`: profile a single VUMPS macro-iteration on a native-spinful
   (`site_type=1`) `n_uc=2` chain to find the actual per-iteration bottleneck (dense 4-dim local
   tensors? Jordan-Wigner string overhead specific to `ElectronSite`'s `F` operator? something
   not vectorized for this site type?) — since a finite chain ~10x larger does full DMRG in
   comparable wall-clock, this looks like a fixable inefficiency rather than an inherent cost.
3. Either way, `tests/test_infinite_chain.py`'s fermionic block should grow a `site_type=1`
   case alongside its existing `site_type=0` one, so this combination has real coverage instead
   of silently falling through untested.

## Workaround adopted in the calling project

Dropped infinite DMRG for this model; using finite chains instead (both open and periodic
boundary conditions — PBC specifically to get a bulk gap uncontaminated by edge states, since
open chains here host genuine near-zero-energy edge modes).
