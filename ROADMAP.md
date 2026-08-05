# DMRGPY Roadmap: Backend Capability Matrix

DMRGPY exposes the same Python API (`Many_Body_Chain` and friends) over
several interchangeable computational backends, selected via
`itensor_version=`. This document tracks **what each backend can
actually do today** and where the gaps are. See `CLAUDE.md` for the
architectural background (`mode.py` dispatch, the `cpp_handle` pattern,
etc.) — this file is a feature-by-feature status board, not a design doc.

## The three backends this roadmap covers

- **ITensor v3 (`itensor_version=3`, the default)** — in-process pybind11
  extension (`mpscpp3/`) over a vendored copy of real ITensor v3 (C++).
  Fastest backend for most iterative methods once compiled; needs a C++
  toolchain + LAPACK/BLAS (`python install.py`).
- **pyitensor (`itensor_version="python"`)** — a from-scratch, pure
  Python/NumPy/SciPy reimplementation of the same ITensor v3 API subset
  (`pyitensor/`). Zero compiler dependency, always available. Slower than
  compiled ITensor per-operation, but is where most *new* algorithms in
  this codebase land first (easier to prototype/debug than C++), and has
  actually overtaken `itensor_version=3` for several iterative workloads
  after the kernel/contraction-order optimization pass (see
  `pyitensor/kernels.py`).
- **Julia (`itensor_version="julia_live"`)** — a live, in-process Julia
  session (`mpsjulialive/`, via `pyjulia`) driving real ITensors.jl. Its
  own parallel, hand-written set of modules mirroring the top-level ones
  — a feature only exists here if someone has explicitly ported it; there
  is no automatic fallback into this backend for anything.

Two other backends exist but are out of scope for this roadmap's main
matrix: **ITensor v2** (`itensor_version=2`, the original/legacy compiled
backend — same session API as v3 but missing everything built on v3's
`TDVP/` module: no TDVP, no TEBD, no generalized eigenproblem, no
NH-DMRG, no METTS) and **ED** (`mode="ED"` or the automatic fallback in
`mode.py` when a requested C++ extension isn't compiled) — exact
diagonalization, always available, the correctness reference every other
backend is checked against, but exponential in system size.

Legend: ✅ implemented and wired into the public API · 🟡 implemented but
restricted (see note) · ❌ not implemented (raises `NotImplementedError`
or is simply absent) · — not meaningful for this backend/method combo.

## 1. Ground state / core DMRG

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| Ground-state energy/wavefunction (`gs_energy`/`get_gs`) | ✅ | ✅ | ✅ | All three run real two-site DMRG. v3 aborts (`SIGABRT`) for chains with <3 sites — worked around by falling back to ED automatically (`mode.py`), not a v3-native fix. |
| Excited states, overlap-penalty method | ✅ | ✅ | ✅ | Julia has its own implementation (`mpsjulialive/excited.py`), not a shared code path. |
| Excited states, non-Hermitian (Arnoldi/IRAM, `arpacktk`) | ✅ | ✅ | ✅ | Backend-agnostic — built only on generic `MultiOperator * MPS` application, so it rides on whichever backend supplies that. |
| Generalized eigenvalue DMRG (Lagrange-multiplier trick, Hermitian) | ✅ | ✅ | ✅ | Julia via `mpsjulialive/generalized.jl`'s `get_gs_generalized` (same algorithm against ITensorMPS.jl's own `dmrg()`/`Sweeps`/`add()`). `groundstate.gs_energy_generalized` still raises for `itensor_version=2`. No ED fallback either (there is no ED implementation of this method at all). |
| Non-Hermitian DMRG (biorthogonal, complex energies) | ✅ | ✅ | ✅ | Julia is the one backend that isn't a port of ITensorNHDMRG.jl — it calls the real package (`mpsjulialive/nhdmrg.jl`), with two corrections on top: its left vector solves the *transpose* equation (needs a conjugation), and its left solve isn't anchored to the right solve's eigenvalue, so a conjugate pair tied for the smallest real part needs an `exp(i*theta)*H` tie-break. Both are invisible on a complex-*symmetric* H — see `tests/test_nh_dmrg.py::nh_asymmetric_hopping_chain`. |
| Non-Hermitian generalized eigenproblem | ✅ | ✅ | ✅ | Same story as above, one level up (`gs_energy_generalized_nhdmrg`; `mpsjulialive/generalized.jl`'s `get_gs_generalized_nhdmrg` wraps the same outer loop around real ITensorNHDMRG.jl sweeps). v2 still excluded. |
| iDMRG (infinite chain, ground state) | ✅ | ✅ | ❌ | `infinitechain.py` hard-restricts to `itensor_version in ("python", 3)` at construction time; v2 and Julia have no port at all. |
| iDMRG vev / two-point correlator | ❌ | ✅ | ❌ | pyitensor only — v3's iDMRG session has no equivalent to `two_point_correlator`/`onsite_expectation` yet. |
| iDMRG excitation ansatz (quasiparticle, tangent-space) | ❌ | 🟡 | ❌ | pyitensor only, and only `D=1` bond dimension (`D>1` gives an unresolved anomalous dispersion — deliberately descoped, see `NotImplementedError` for `D>1`). |
| iDMRG local excitation gap (deflation-based, `D`-capable) | ❌ | ✅ | ❌ | pyitensor only (`idmrg.py::local_excitation_gap`/`local_excitation_gap_windowed`). |
| iDMRG dynamical correlator (finite-window KPM reduction) | ❌ | ✅ | ❌ | `infinitechain.py::kpm_finite` is gated to `itensor_version="python"` only. |
| iDMRG real-time window dynamics (IBC-style TDVP, `td_dynamical_correlator`) | ✅ | ✅ | ❌ | Both wired (`idmrg_window.py` for pyitensor, `mpscpp3/chain_session.h`'s `td_dynamical_correlator_window` for v3); Phase 2+ generic sweep machinery beyond the current window construction is still pyitensor-only groundwork (`idmrg_window.py`'s own Phase 0-1 notes). |
| iDMRG cat-state superposition (`imps_sum` across distinct branches) | ❌ | ❌ | ❌ | Two *identical* branches sum fine (`imps_sum`); a true superposition of two physically distinct symmetry-broken ground states needs new correlator machinery not yet written for any backend. |

## 2. Expectation values, static correlators, entanglement

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| `vev` (operator expectation value) | ✅ | ✅ | ✅ | |
| Two-point static correlators | ✅ | ✅ | ✅ | |
| All-pairs correlation matrix (native accelerated) | ✅ | ✅ | 🟡 | Julia's `MPS` class doesn't expose a dedicated accelerated method the way `pyitensor.chain.Chain.correlation_matrix`/v3's `Chain::correlation_matrix` do — falls back to the generic per-pair loop. |
| Four-point correlator, `ctmode="explicit"` (backend-agnostic loop) | ✅ | ✅ | ✅ | Always correct, always available, slowest — the universal fallback. |
| Four-point correlator, `ctmode="full"` (native AutoMPO per-tuple, incl. spinful-native sites) | ✅ | ✅ | ❌ | v2 also has this; Julia does not. |
| Four-point correlator, `ctmode="sweep"` (single-sweep, environment-reuse, ITensorCorrelators.jl-style) | ✅ | ✅ | ❌ | The fastest mode when it applies; no native-spinful-site counterpart yet even on v3/pyitensor. |
| Bond entanglement entropy / reduced density matrix | ✅ | ✅ | ✅ | Backend-agnostic dispatch (`entanglement.py` has no `itensor_version` branching at all — every backend's session exposes the same primitive). |
| Topological invariants (Berry phase, etc.) | ✅ | ✅ | ✅ | Built purely on `get_gs()` + MPS overlap — generic across every backend. |

## 3. Time evolution

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| Real-time evolution, TDVP (two-site) | ✅ | ✅ | ✅ | Julia has its own implementation (`mpsjulialive/timedependent.py` + `tdvp.jl`), and now honors `self.tevol_method` rather than ignoring it — see the TDVP_GSE row; anything outside `("TDVP","TDVP_GSE")` raises there instead of silently running plain TDVP. |
| Real-time evolution, TDVP + global subspace expansion (`tevol_method="TDVP_GSE"`) | ✅ | ✅ | ✅ | Julia needs no port: ITensorMPS.jl ships both pieces (`expand(psi,H; alg="global_krylov")`, citing the same Yang-White paper, and `tdvp(...; nsite=1)`), so `mpsjulialive/tdvp.jl` only wires them into the same expand-then-step structure. Unlike v3 it handles a bond-dimension-1 (product-state) start fine — see `examples/time_evolution/tdvp_gse_julia_VS_ED_time_evolution`. |
| Real-time evolution, TEBD (2nd-order Trotter gates) | ✅ | ✅ | ❌ | Restricted to `itensor_version in (3,"python")`; v2 has no TDVP module at all so was never a candidate either. Julia would need a Trotter-gate primitive (bond Hamiltonians built from the term list, with explicit fermionic-sign handling that AutoMPO does for free elsewhere) — the one piece of this group still missing there. |
| Real-time evolution, MPO-Taylor (`tevol_method="MPO"`, the pre-TDVP legacy path) | ✅ | — | — | The only option for v2; still selectable on v3 by explicitly setting `tevol_method="MPO"`. |
| Complex-time evolution / TDZ damping (arXiv:2311.10909) | ✅ | ✅ | ✅ | Julia needed exactly one new primitive (`advance_complex_time_step`); the rest of `tdz.py` is backend-agnostic MPS/MPO algebra. |

## 4. Dynamical correlators (`get_dynamical_correlator(submode=...)`)

| submode | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| `KPM` (Kernel Polynomial Method, Chebyshev) | ✅ | ✅ | ✅ | Julia's is a separate implementation (`mpsjulialive/dynamics.py`), Hermitian-only. |
| `KPM` energy truncation (Holzner Sec. III-B) | ✅ | ✅ | ❌ | Net *slower* on both backends where it's measured, kept for regimes/paper-fidelity, not a default. |
| Non-Hermitian KPM (NHKPM.jl-ported coupled recursion) | ✅ | ✅ | ❌ | Auto-selected instead of plain KPM whenever the Hamiltonian is non-Hermitian, for backends that support it; ED also has this path. |
| `CVM` (correction-vector / conjugate-gradient resolvent) | ✅ | ✅ | ✅ | Julia's CVM needed no new primitive beyond `set_sweep_params`-equivalent bookkeeping — already backend-agnostic MPS/MPO algebra. |
| `CVM_explicit`, `CVMimag` (analytic continuation variant) | ✅ | ✅ | ❌ | Only plain `CVM` is wired for Julia so far. |
| `TDZ` (complex-time evolution correlator) | ✅ | ✅ | ✅ | |
| `TD` (real-time evolution correlator) | ✅ | ✅ | ✅ | Dispatches through `timedependent.py`, which has its own Julia branch. |
| `ROOTN` (root-N Krylov correction-vector) | ✅ | ✅ | ❌ | Per-bond-sweep shortcut confirmed impossible on both C++ and pyitensor (needs real multi-target state-averaging) — not attempted for Julia. |
| `EX` (exact spectral decomposition over excited states) | ✅ | ✅ | ✅ | Non-Hermitian-capable; entirely backend-agnostic (`self.get_excited_states()` + generic MPS algebra). |
| `maxent` (max-entropy positive-definite reconstruction) | ✅ | ✅ | ✅ | Same — backend-agnostic post-processing over `power_vev` moments. |

## 5. Finite-temperature methods

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| Ancilla purification + imaginary-time annealing (`thermal.py`) | ✅ | ✅ | ✅ | Backend-agnostic — just wraps a doubled-site `Spin_Chain` and calls its generic `get_gs`/evolve methods. |
| Exact Boltzmann sum over ED excited states (`thermalvev.py`) | — | — | — | ED-only by construction (needs the full spectrum); not a DMRG-backend feature. |
| METTS thermal average (`metts_vev`) | ✅ | ✅ | ❌ | Needs imaginary-time TDVP, so v2 (no `TDVP/`) and Julia (no port) are both excluded. `njobs>1` multiprocess pooling is pyitensor-only (no equivalent lever for a single live v3/Julia session). |
| METTS dynamical correlator (finite-T real-time correlator) | ✅ | ✅ | ❌ | Same restriction as `metts_vev`. |

## 6. Physical models (Hilbert-space coverage)

Model coverage is largely backend-independent — `spinchain.py`,
`fermionchain.py`, `spinfermionchain.py`, `bosonchain.py`,
`parafermionchain.py` all build their operators as backend-agnostic
`MultiOperator`s, so any new model works with any backend that supports
the *method* being called on it (per the tables above). The one
model-specific exception:

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| Native spinful-fermion sites (`Spinful_Fermionic_Chain_Native`, real `ElectronSite`) | ✅ | ❌ | ❌ | Only v3 wires up a genuine per-site `Cup/Cdn/Cdagup/Cdagdn` site type; the other backends only offer the two-sites-per-orbital ("interleaved") representation. Slower for iterative GS/TDVP/KPM, faster for the 4-point tensor — a real perf tradeoff, not a strict upgrade, so both representations stay available where implemented. |

## What's missing, roughly ranked by how load-bearing it'd be

1. **Julia: iDMRG.** No port at all (`infinitechain.py` rejects
   `itensor_version="julia_live"` outright). Given how much iDMRG
   machinery (excitation ansatz, window dynamics, warm-start correlators)
   has landed on pyitensor+v3 recently, a from-scratch Julia port is a
   large undertaking, not a quick win — likely not worth it unless a
   concrete use case needs Julia specifically for infinite systems.
2. **Julia: TEBD and METTS.** TEBD needs a Trotter-gate evolution
   primitive (and, for fermionic models, the explicit Jordan-Wigner sign
   handling inside a two-site gate that ITensor's AutoMPO does for free
   on the MPO path); METTS needs an imaginary-time sampler on top of the
   imaginary-time TDVP step `advance_complex_time_step` already provides.
   Both are incremental rather than architectural gaps — TDVP_GSE, the
   third member of this group, is done.
3. **v3/pyitensor: four-point `ctmode="sweep"` for native-spinful sites.**
   The fast environment-reuse algorithm has no flavor-resolved
   counterpart yet; `ctmode="full"` remains the practical choice for
   `Spinful_Fermionic_Chain_Native`.
4. **iDMRG excitation ansatz beyond `D=1`.** Known to be broken (anomalous
   flat dispersion), not just unimplemented — needs real debugging, not
   just porting effort, before `D>1` can be un-blocked.
5. **iDMRG cat-state superpositions.** No backend supports summing two
   *physically distinct* symmetry-broken iDMRG ground states; needs new
   correlator machinery from scratch.
6. **v2 (legacy) parity.** Deliberately not being chased — v2 predates
   `TDVP/` entirely and is kept mainly as the historical QN-conserving
   reference implementation, not a target for new features.

## Pointers for extending a backend

- **v3 (C++)**: add/extend a method on `Chain` in
  `src/dmrgpy/mpscpp3/chain_session.h`, expose it through
  `src/dmrgpy/mpscpp3/bindings.cc`, then wire the Python-side dispatch
  (usually a `self.itensor_version==3` branch next to the existing
  `"python"`/`"julia_live"` ones in the relevant top-level module).
- **pyitensor**: add the method to the relevant `pyitensor/*.py` module
  (`dmrg.py`, `tdvp.py`, `idmrg.py`, `metts.py`, `nhdmrg.py`, ...) and to
  `pyitensor/chain.py`'s `Chain` class, which is the thing every
  top-level dispatch module actually calls into for `itensor_version=
  "python"`.
- **Julia**: add a mirrored module under `src/dmrgpy/mpsjulialive/`
  (see `dynamics.py`/`timedependent.py` for the pattern of "reuse
  backend-agnostic MPS/MPO algebra wherever possible, write a Julia
  primitive only for the one piece that genuinely needs it"), then add
  an `itensor_version=="julia_live"` branch at the relevant top-level
  dispatch point. A `.jl` file holding the actual algorithm goes next to
  it and gets appended to `juliasession.py`'s `initialize()` file list
  (order only matters up to first *use*, since Julia resolves globals at
  call time). `generalized.jl`+`generalized.py` / `nhdmrg.jl`+`nhdmrg.py`
  are the most recent worked example of the whole shape, including the
  case where the Julia ecosystem already has the algorithm and the work
  is reconciling *its* conventions with dmrgpy's rather than porting
  anything.

After adding a capability to a backend, update this file's matrix and,
if it's a new physics-facing method, `docs/user_guide.md`/`.tex` per
`CLAUDE.md`'s documentation policy.
