---
paths:
  - "src/dmrgpy/mpscpp2/**"
  - "src/dmrgpy/mpscpp3/**"
  - "src/dmrgpy/cppext.py"
  - "src/dmrgpy/mps.py"
  - "src/dmrgpy/multioperatortk/staticoperator.py"
  - "tests/test_audit_2026_09_25b_cpp.py"
---

# The in-process pybind11 extension (the C++ DMRG backend, v2 and v3)

This is the section CLAUDE.md used to carry as "In-process pybind11
extension"; it loads when a file under `mpscpp2/` or `mpscpp3/` is read or
edited. The vendored ITensor copies under `mpscppN/ITensor/` are upstream
code and are not read or modified here; the dmrgpy-owned C++ is
`chain_session.h`, `bindings.cc`, `check_task.h`, `get_sites.h` and
`mo_terms.h` in each `mpscppN/`, plus `mpscpp3/TDVP/`.

`mpscpp2/bindings.cc` and `mpscpp3/bindings.cc` each compile a pybind11
extension module (`mpscpp2/_dmrgcpp*.so` and `mpscpp3/_dmrgcpp*.so`, built
by the `pybind` target in each backend's own `Makefile`) that runs DMRG
entirely in-process against its own vendored ITensor copy: no task files, no
operator files, no subprocess, no temp directory. `mpscpp3` is a
line-by-line port of `mpscpp2` to the ITensor v3 API (IQTensor/IQIndex
merged into a single ITensor/Index type with QN blocks living on the Index;
see the porting-notes comment at the top of `mpscpp3/chain_session.h` for
the full v2 to v3 API mapping) exposing the identical `Chain` class and
Python surface, not an independent reimplementation. There is no
file-based/subprocess design for either version; that was removed entirely
from `mpscpp2`, and `mpscpp3` never had one.

Both `mpscppN/Makefile`s list `$(wildcard *.h)` among the extension's
prerequisites. Without that, editing `chain_session.h`, where the actual
DMRG implementation lives, leaves `make pybind` reporting "Nothing to be
done" while the `.so` keeps the old code, which is a confusing way to test a
change that was never compiled. A C++ change is not tested until `make
pybind` has rebuilt the `.so` and the test ran against the rebuilt one.

## The session and its pieces

- **`cppext.py`** lazily imports whichever compiled extension a chain needs
  and caches whether it loaded (`get_backend(version)`/`available(version)`,
  `version` is `2` or `3`).
- **`mpscppN/chain_session.h`**'s `Chain` class is the stateful, in-process
  session: one `Chain` instance per `Many_Body_Chain` (held as
  `self._session`, created in `sites.py::initialize` whenever
  `itensor_version in (2,3)` and the corresponding extension is available),
  with a method per DMRG task (`gs_energy`, `vev`, `apply_operator`, KPM/CVM
  dynamical correlators, `correlation_matrix`, `reduced_dm`,
  `excited_states`, time evolution, ...), built directly against ITensor's
  `AutoMPO`/`MPS`/`MPO` APIs. `check_task.h`, `get_sites.h` and `mo_terms.h`
  are the only other `mpscppN/*.h` headers it depends on; everything else it
  needs is self-contained in `chain_session.h`.
- **`mpscppN/mo_terms.h`** and `MultiOperator.to_terms()`
  (`multioperator.py`) convert a `MultiOperator` into a plain list of
  `(coefficient, [(opname, site), ...])` tuples that `Chain` turns directly
  into an ITensor `HTerm`/`AutoMPO`.
- **The `cpp_handle` pattern**: `mps.py::MPS` and
  `multioperatortk/staticoperator.py::StaticOperator` both carry an opaque
  `cpp_handle` attribute holding the pybind11-wrapped ITensor object living
  in C++ memory. These classes have no other representation now, so
  `cpp_handle` is always set for `itensor_version in (2,3)` chains (it is
  `None` only for objects that were never meant to touch the C++ backend,
  e.g. `mpsjulialive`'s own `MPS`/`MPO` classes are separate types
  entirely).
- Because `cpp_handle` objects have no pickle/deepcopy support,
  `MPS.copy()`/`StaticOperator.copy()` use a shallow `copy.copy()` (never
  `deepcopy()`). `Many_Body_Chain.__deepcopy__` never copies `_session`
  directly (a live C++ session has no well-defined "copy" semantics);
  instead it builds a fresh `Chain` for the clone when the original had one,
  which is what `clone()`/`bandwidth()`/`lowest_eigenvalue()` rely on.
- **Both backends can be loaded in the same process** (needed for anything
  that directly compares `itensor_version=2` against `=3`, e.g. the
  `v2_VS_v3_*`-named examples in each theme folder, or
  `examples/dynamical_correlator/dynamical_correlator_v2_VS_v3`) only
  because `mpscpp3/bindings.cc` registers its pybind11 types with
  `py::module_local()`. Without it, `mpscpp2` and `mpscpp3` define their
  own, ABI-incompatible `Chain`/`MPS`/`MPO` C++ types (same unqualified,
  unnamespaced names, different vendored ITensor copies), and libstdc++'s
  cross-DSO RTTI equality falls back to comparing the mangled type-name
  string, not the address, so pybind11 (keyed on `std::type_index`) sees
  the exact same "type" already registered the moment both extensions are
  imported into one process, and aborts with `generic_type: type "Chain" is
  already registered!`. `module_local()` keeps mpscpp3's registration in a
  module-scoped table instead of the shared global one, avoiding the
  collision regardless of import order (confirmed both ways). This is
  applied only on the mpscpp3 side; mpscpp2 is untouched and still uses the
  ordinary global registration.

## Scale, thresholds and reproduced quirks

- **Both backends hand ITensor their operators at unit scale.** Three
  thresholds inside ITensor are absolute (svdMPO's `Cutoff` of 1e-13 on
  squared singular values, its 1e-14 coefficient skip, Davidson's 1e-10
  residual below which a Krylov direction is randomized), so `build_mpo()`
  goes through `mo_terms.h`'s `to_mpo_unit()` (`unit_scale_up` is the
  factor) and every `dmrg()` on the Hamiltonian through
  `solver_hamiltonian()`, both at the power of two that
  brings a largest coefficient below 1 into [1,2), divided back exactly; at
  a largest coefficient of 1 or more they run the unscaled code byte for
  byte, and `set_hamiltonian_mpo` resets the scale to 1. The session's
  `hscale_up_` is the factor `to_mpo_unit()` actually applied to the merged
  AutoMPO (its `up_out`), not a maximum over the raw terms, which cancelling
  duplicates made wrong (2026-09-25b finding 28). A new MPO built from an
  AutoMPO, or a new `dmrg()` call, should go through the same two, or it
  reintroduces the small-units failure (2026-09-25 record, item 2). This
  fixes the Hamiltonian's units but not a hierarchy inside it: a bond whose
  strongest term is far below the largest coefficient still loses channels
  on v2, v3 and `"python"`, held by strict xfails. The cause at small units
  is ITensor's own `toMPO`, whose svdMPO truncates on an absolute cutoff,
  plus Davidson's absolute 1e-10 randomization threshold below it, not the
  `"python"` MPO builder.
- **NH-DMRG is not a `dmrg()` call.** It is unit-scaled once, at the Python
  entry in `nhdmrg.py`, for all three session backends, so do NOT add a
  second scale inside `Chain::nhdmrg`.
- **NH-DMRG's Ritz tie window** is `1e-6*(max Re - min Re) +
  100*eps*max|ev|`, compared with `<=`, identically in `pyitensor/nhdmrg.py`
  and both `arnoldi_select_kbest`; keep the two spellings identical. Its
  certificate divides by `c + |E - e_id|`.
- **`Chain::vx_lanczos_ground_state`'s residual reference** `max(1,|lambda|)`
  is capped by the Hamiltonian's own unit, its largest non-identity
  coefficient, so that it holds in small units and next to an offset
  (2026-09-25b finding 24); the `"python"` VUMPS has the same cap.
- **`verbose` v2/v3 runs print the solver scale**, one line before the solve.
- **Pre-existing bugs in the original design are deliberately reproduced,
  not fixed**, in both `chain_session.h`s (see comments at the call sites
  for details): `evoloperator()`'s z^3/6 term multiplies `H2` instead of the
  computed `H3`; the old file-based backend's `"noise"` tasks.in key was
  actually read back as `"moise"`, and `"tevol_fit"` was written but
  `"tevol_fit_td"` was read, making that fitting branch unreachable. Both
  are reproduced as hardcoded values rather than silently "fixed". If
  ITensor results look numerically "close but not exact" to a hand
  derivation, check here before assuming a porting bug. The hole-hunt
  process treats these three as out of scope by construction.
- **Still unedited after the 2026-09-25 rebuild**: the
  `vx_deterministic_start` comment in `mpscpp3/chain_session.h`, recorded as
  open by the 2026-09-24 audit, is left to a pass that rebuilds.

## `mpscpp3`-specific behaviour

- **Every site is built with `ConserveQNs=false` by default**
  (`mpscpp3/get_sites.h`), and `Chain`'s DMRG starting state is an actual
  `randomMPS(sites_,maxm_)`, not a plain product state (`chain_session.h`'s
  `default_mps()`). This is not cosmetic: v2's own plain `MPS(sites_)` (no
  `InitState`) ends up, in practice, performing an unconstrained search of
  the full Hilbert space rather than one pinned to a single total-Sz or
  particle-number sector, confirmed by cross-checking against the compiled
  v2 extension with a symmetry-breaking field term, where v2 correctly finds
  a fully-polarized ground state. ITensor v3 is stricter about
  QN-conserving tensors needing a well-defined flux from the start, and a
  naive port (QN on, product-state start) gets DMRG permanently stuck at
  that starting product state's energy (DMRG trapped at an exact
  eigenstate, reproducible even with a large noise term) instead of raising
  an error. Turning QN off and starting from an actual random MPS reproduces
  v2's real behavior at the cost of the QN block-sparsity speedup for
  `itensor_version=3`: a genuine performance and memory tradeoff of this
  backend, not a bug to "fix" later. Note the "permanently stuck"
  observation is specific to the all-first-basis-state product start it was
  measured with: that state (all spins up) is the unique state of its own
  Sz=N/2 sector, so no noise term can move it and its energy is the correct
  answer for that sector. It says nothing about QN mode from a non-trivial
  in-sector start, which is what the opt-in sector mode below actually does.
- **Opt-in conserved-sector (QN) mode.**
  `Many_Body_Chain.set_conserved_sector(Nf=6)`/`(Sz=0)`/`(Nf=8,Sz=0)` (off
  by default, no-args clears it) confines the whole calculation to one
  quantum-number sector: `Chain::set_conserved_sector` rebuilds `sites_`
  via `SpinX(site_types,conserved)` with QN-carrying indices and drops
  everything built on the old site set. The same public API (including
  `promote_to_dense`/`promote_mps`) is also implemented by
  `itensor_version="python"`, by a different mechanism (dense storage with
  a charge grading on the site indices plus a charge penalty on the
  variational solves, see `pyitensor/sector.py`'s module docstring and the
  pure-Python-backend section of `docs/development.md`), and by
  `mode="ED"` by a third (`edtk/edchain.py`: a conserved charge is diagonal
  in the ED product basis, so a sector is a set of basis states and every
  assembled operator is restricted to that submatrix; `promote_to_dense`/
  `promote_mps` are deliberately not part of it, since an ED re-solve is
  cheap). Everything below is v3-specific. It pays off with size, and only
  with size: measured through the Python API on Heisenberg chains, sector
  against dense is 0.6x (n=20, maxm=60), 2.0x (n=40, maxm=100), 4.2x (n=60,
  maxm=200), since block sparsity scales and its per-block bookkeeping does
  not. Do not quote a single speedup number. Three ITensor behaviors shape
  the implementation, each confirmed directly and none of them a clean
  failure:
  - `randomMPS(SiteSet)` refuses to guess a sector under QNs and
    `randomMPS(InitState,m>1)` does not exist here (`mps.cc`), so the start
    state is a normalized `sum()` of random in-sector product states
    (`sector_mps`), with the per-site assignment coming from an exact DP
    over reachable partial charges (`sector_state_plan`; greedy dead-ends
    on a Hubbard chain at fixed `Nf` and `Sz`).
  - Building a non-QN-definite operator (`Sx` on Sz-conserving sites)
    aborts the process inside `ITensor::set`, and AutoMPO aborts on a
    flux-violating term even when its coefficient cancels exactly ("Index
    does not contain given QN block"), which is why the
    `Sx*Sx+Sy*Sy+Sz*Sz` Heisenberg Hamiltonian needs `mo_terms.h`'s
    `expand_xy_terms`+`combine_terms` normalization here rather than being
    left to AutoMPO, and why `op_charge()` infers an operator's charge from
    its dense matrix elements instead of building it on QN indices. Every
    terms-to-MPO path routes through `mpo_from_terms`/`ampo_from_terms` so
    the check covers `vev`/correlators/KPM vertices too; both are
    pass-through when no sector is set, keeping the default path
    byte-identical.
  - ITensor's own `Error()` calls `abort()` (`util/error.h`), so nothing it
    raises is catchable from Python; the sector code throws
    `std::invalid_argument`/`std::runtime_error` instead.
  Python side: the sector lives on the chain, is re-applied wherever a
  session is rebuilt (`sites.py::initialize`, `__deepcopy__`), is part of
  `groundstate.py`'s send-cache key, and `mode.py` raises on DMRG with
  `itensor_version=2` or `"julia_live"`, the two backends with no quantum
  numbers, rather than letting one of them answer with the global ground
  state; falling back to ED (missing extension, short chain) is correct,
  since ED restricts to the sector itself. `tests/test_sector_conservation.py`
  cross-checks every sector against ED restricted to that sector. Note the
  reference must restrict (the charge operator is diagonal in the ED basis,
  so a sector is a submatrix), not filter full-spectrum eigenvectors by
  <N>: sector-degenerate eigenvalues come back as arbitrary superpositions.
- **`applyMPO()` does not restore the standard unprimed physical index** the
  way v2's `exactApplyMPO`/`fitApplyMPO` implicitly did. It just contracts
  whichever leg of the MPO happens to match the input MPS's own (possibly
  already-primed) physical index, so the next `sum()` of two
  independently-produced results can hit a hard
  `ITensor::operator+=: different index structure` abort. Every MPO application in
  `mpscpp3/chain_session.h` goes through a private `apply_mpo()` wrapper
  that always calls `.noPrime("Site")` on the result before returning,
  restoring that invariant globally instead of patching individual call
  sites. Likewise, `nmultMPO(A,B)` (used to build `H*H` for
  `custom_exp()`/`evoloperator()`) requires `nmultMPO(A,prime(B))` when `A`
  and `B` share the same physical indices (any two operators built on the
  same chain always do), where v2 tolerated the unprimed form directly, and
  the result's output leg then needs `.mapPrime(2,1)` to get back to the
  standard single-application convention before it can be summed against
  any other operator (see `mult_mpo()`'s comment, and ITensor v3's own
  `unittest/mpo_test.cc` "nmultMPO" section, which documents exactly this
  usage).
- **Plain `inner()` throws on complex MPS/MPO** ("Cannot use inner(...) with
  complex MPS/MPO, use innerC(...) instead"), unlike v2's `overlap()` (which
  just silently returns the real part regardless). Since `MOTerm::coef` is
  always `Cplx`, any MPO/MPS built from Python-supplied terms is
  complex-typed even when every coefficient's imaginary part is exactly 0;
  confirmed directly, the KPM dynamical correlator path aborted here in
  `same_mps()` until every "Re[<x|y>]" use in `chain_session.h` was changed
  to `innerC(...).real()`.
- **Real-time MPS evolution defaults to TDVP.** `mpscpp3/TDVP/` (vendored
  from ITensor's own TDVP repo, see its `README.md`) provides a proper
  two-site TDVP integrator, wired into `chain_session.h` as
  `Chain::quench_tdvp()`/`Chain::evolve_and_measure_tdvp()` (via the private
  `tdvp_step()` helper) alongside the pre-existing `Chain::quench()`/
  `Chain::evolve_and_measure()`, which apply a hand-rolled 2nd-order Taylor
  expansion of `exp(-i dt H)` as an MPO (`evoloperator()`) instead.
  `Many_Body_Chain.tevol_method` (`manybodychain.py`, default `"TDVP"`)
  picks between them in `timedependent.py`'s `evolution_dmrg_DC()`/
  `evolve_and_measure_dmrg()`; every one of those dispatches tests
  `itensor_version in (3,"python")`, so `"TDVP"` applies on the C++ side
  only at `itensor_version==3` (mpscpp2 has no `TDVP/` and no `_tdvp`
  methods at all), and any other combination silently falls back to the
  MPO-Taylor path, which remains the only option for `itensor_version=2`.
  `Chain::tdvp_step` takes `DoNormalize` from the step (`dt.imag()==0.0`),
  not hardcoded: forcing unit norm after a COMPLEX-time step deletes exactly
  the decay `tdz.py`'s `submode="TDZ"` contour is built on, which is what
  made v3's TDZ return too much spectral weight (the 2026-09 audit measured
  36% too much, 0.340676 against an exact 0.25, finding #7; the post-fix
  0.249473 re-measures, the "before" figure does not now that the extension
  has been rebuilt). Only the two-site TDVP algorithm is used
  (`NumCenter=2`); the global subspace expansion machinery in
  `TDVP/basisextension.h` (mainly useful for one-site TDVP or long-range
  Hamiltonians) is not wired in, since two-site TDVP already grows bond
  dimension via SVD like two-site DMRG.
- **`evolve_and_measure_*` loops restore the input norm `norm0` after each
  step**, as every `quench_*` loop always did; v3's copy scales by
  `norm0/nrm` rather than calling `MPS::normalize()`, which aborts below an
  absolute 1e-20 (2026-09-26 review).
- **In `Chain::global_subspace_expand` a `Cutoff=0` sweep is not "no
  cutoff"**: ITensor v3's `truncate()` discards exactly-zero weights at
  `Cutoff=0`, which is what the expansion adds next to a site pinned to one
  basis vector, so the sweeps run only when a bond exceeds `maxm`, with
  `Cutoff=-1` (2026-09-24 third pass, finding 14; this was also the
  product-state failure the TDVP_GSE example had recorded as isolated).
- **The two-site `dmrg()` aborts the process below 3 sites** (`SIGABRT`,
  "LocalOp is default constructed", an internal check in
  `itensor/mps/localop.h`), for any chain with fewer than 3 physical sites
  and independent of physics type. `mode.py` routes `itensor_version==3`
  with `self.ns<3` to ED for that reason; see CLAUDE.md's "Running /
  verifying changes".
- **The iDMRG and VUMPS helpers in `chain_session.h`** (`idmrg_*`, `ic_*`,
  `vx_*`, `vms_*`) are covered in `.claude/rules/infinite-chain.md`, which
  loads together with this file when `mpscpp3/chain_session.h` is edited.
