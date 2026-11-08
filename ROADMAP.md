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
  session (`mpsjulialive/`, via `juliacall`) driving real ITensors.jl. Its
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
| Controlled bond expansion ground state (`dmrg_method="CBE"`) | ✅ | ✅ | ✅ | One-site DMRG with the shrewd-selection expansion: `pyitensor/cbe.py`, `mpscpp3/cbe.h`, `mpsjulialive/cbe.jl`. Not in a conserved sector on v3 and Julia. On Julia it is about 7x slower than two-site `dmrg()` at 24 sites. |
| Excited states, overlap-penalty method | ✅ | ✅ | ✅ | Julia has its own implementation (`mpsjulialive/excited.py`), not a shared code path. |
| Excited states, non-Hermitian (Arnoldi/IRAM, `arpacktk`) | ✅ | ✅ | ✅ | Backend-agnostic — built only on generic `MultiOperator * MPS` application, so it rides on whichever backend supplies that. |
| Generalized eigenvalue DMRG (Lagrange-multiplier trick, Hermitian) | ✅ | ✅ | ✅ | Julia via `mpsjulialive/generalized.jl`'s `get_gs_generalized` (same algorithm against ITensorMPS.jl's own `dmrg()`/`Sweeps`/`add()`). `groundstate.gs_energy_generalized` still raises for `itensor_version=2`. No ED fallback either (there is no ED implementation of this method at all). |
| Non-Hermitian DMRG (biorthogonal, complex energies) | ✅ | ✅ | ✅ | Julia is the one backend that isn't a port of ITensorNHDMRG.jl — it calls the real package (`mpsjulialive/nhdmrg.jl`), with two corrections on top: its left vector solves the *transpose* equation (needs a conjugation), and its left solve isn't anchored to the right solve's eigenvalue, so a conjugate pair tied for the smallest real part needs an `exp(i*theta)*H` tie-break (anchored on the untied run's own eigenvalue, since the rotation can otherwise re-target a different one). Both are invisible on a complex-*symmetric* H — see `tests/test_nh_dmrg.py::nh_asymmetric_hopping_chain`. |
| Non-Hermitian generalized eigenproblem | ✅ | ✅ | ✅ | Same story as above, one level up (`gs_energy_generalized_nhdmrg`; `mpsjulialive/generalized.jl`'s `get_gs_generalized_nhdmrg` wraps the same outer loop around real ITensorNHDMRG.jl sweeps). v2 still excluded. |
| Conserved sectors (`set_conserved_sector(Nf=..)`/`(Sz=..)`, `promote_to_dense`, `promote_mps`) | ✅ | ✅ | ✅ | Julia uses ITensors.jl QN site indices (`mpsjulialive/sector.py`, `sites_qn.jl`), so the confinement is structural as on v3 and needs no penalty. QN sites exist for spinless fermions, `Electron`, S=1/2, S=1 and bosons; another site type in a sector raises. The charge bookkeeping and the in-sector start state are pyitensor's (`pyitensor/sector.py`), so an operator that breaks the sector raises the same `ValueError`. Checked against ED restricted to the sector in `tests/test_sector_conservation.py` and `tests/test_sector_promotion.py`. |
| `gs_energy` keywords (`maxde=`, `reconverge=True`, `wf0=`, the `restart=False` warm start) | ✅ | ✅ | ✅ | Julia runs the same `maxde=` doubling loop as `gs_energy_single`, keeps the stored state under `maxde=`, and `stored_answer_holds` uses the session backends' rule. `set_hamiltonian(H, restart=False)` makes the stored state the next solve's start; `restart()` is still a cold start. Only a cold solve marks `e0` as the KPM lower band edge. See `tests/test_audit_2026_09_25b_session.py`. |
| Bond-dimension ramp (`bond_ramp`, `bond_ramp_start`) | ✅ | ✅ | ✅ | Julia's `get_gs.jl::make_sweeps_ramped` is the sweep-by-sweep counterpart of pyitensor's `_make_sweeps_ramped`, warm-start floor included (`tests/test_bond_dimension_ramp.py`). |
| `applyinverse` (`A^-1|wf>`) | ✅ | ✅ | ✅ | Julia: `mpsjulialive/bicgstab.jl`, a line-by-line port of pyitensor's `_bicstab`, with sums compressed by SVD rather than by `summps`. |
| `exponential` (`exp(h)|wf>`) | ✅ | ✅ | ✅ | Julia: `nt0` second-order Taylor steps in `mpsalgebra.jl`, the same operator as pyitensor's `_custom_exp` and v3's `custom_exp`, built exactly and compressed at a relative cutoff of 1e-28. |
| `clone()`/`copy()` of a chain | ✅ | ✅ | ✅ | A Julia clone shares the chain's Julia site indices, so states of the clone and of the original contract with each other (`submode="SECTOR"`, `bandwidth()`). |
| iDMRG (infinite chain, ground state) | ✅ | ✅ | ✅ | `infinitechain.py` hard-restricts to `itensor_version in ("python", 3)` at construction time; v2 and Julia have no port at all. Both `gs_method` values are now on v3 too: `"idmrg"` (`Chain::idmrg_ground_state`) and `"vumps"` (`Chain::vumps_ground_state`, C++ port of `pyitensor/vumps.py` — dense LAPACK-based linear algebra, not ITensor tensor-network objects, see that method's own doc comment for why; the default on BOTH backends since 2026-08-08). Cross-checked directly against `itensor_version="python"`'s own VUMPS at `D=1,2,3` on TFIM and Heisenberg (n_uc=1 and 2), matching to ~1e-10 or tighter — see `tests/test_vumps_v3.py`. Julia since 2026-11-08: `gs_method="vumps"` only, `n_uc<=2` with nearest-cell couplings, a dense-array port of `pyitensor/vumps.py` in `mpsjulialive/vumps.jl` (KrylovKit); `gs_method="idmrg"`, longer reach and `n_uc>2` raise. |
| iDMRG vev / two-point correlator | ✅ | ✅ | 🟡 | Both backends cover both `gs_method`s. pyitensor: `"idmrg"` (`idmrg.py`, dominant-right-fixed-point eigenproblem on the gauge-consistent cell of `IDMRGResult`) and `"vumps"` (`vumps.py`, exact contraction on the mixed-gauge `{AC, AR}`, following arXiv:1810.07006 Eq.(34)/(37)-(39)). v3: `"vumps"` through `Chain::vumps_onsite_expectation`/`vumps_two_point_correlator` (a C++ port of the same AC/AR formula, falling through to the `vms_*` pair for the sequential solver) and `"idmrg"` through `Chain::idmrg_onsite_expectation`/`idmrg_two_point_correlator` (dense arrays on the extracted cell `idmrg_cell_`, `ic_*` helpers). `Infinite_Many_Body_Chain.vev`/`correlator` dispatch to all four combinations. Cross-checked against `itensor_version="python"` in `tests/test_vumps_correlator_v3.py`, `tests/test_idmrg_correlator_v3.py` and `tests/test_infinite_chain.py::test_itensor_version3_matches_python_backend`. Julia: `vev`/`correlator` on its VUMPS state (`n_uc<=2`), not on `gs_method="idmrg"`. |
| iDMRG excitation ansatz (quasiparticle, tangent-space) | ✅ | ✅ | ❌ | `Chain::vumps_excitation_energies` (C++ port of `pyitensor/idmrg_excitations.py`, same dense-array approach as `vumps_ground_state`) — requires `gs_method="vumps"` on both backends. Cross-checked directly against `itensor_version="python"` across a full momentum scan at `D=1,2,3` on TFIM (n_uc=1) and Heisenberg (n_uc=2), matching to ~1e-10 (TFIM, gapped) or ~1e-4..1e-7 (Heisenberg, gapless/critical — both backends land on slightly different local optima of the same non-convex restart search, not a discrepancy in the algorithm itself) — see `tests/test_vumps_excitations_v3.py`. Any `D>=1` is supported, matching pyitensor. Both of the solver improvements the pure-Python implementation carries are on v3 too since 2026-09-22 (item 10 below): the channel resolvents are built once per momentum and kept with their LU factorization rather than rebuilt inside every application of `H_eff(k)`, and `H_eff(k)` is solved by Lanczos on its action above `vumps_h_eff_dense_max_` instead of being assembled. A 3-momentum scan at `D=16` on TFIM went from 164.1 s to 1.8 s and at `D=10` on the `n_uc=2` Heisenberg cell from 15.6 s to 1.2 s, the latter on a state that came out the same either way, its dispersion agreeing to 1.2e-11. |
| iDMRG excitation spectral weights / `S(k,w)` | ✅ | ✅ | ❌ | `Infinite_Many_Body_Chain.spectral_weights`/`dynamical_structure_factor` (`pyitensor/idmrg_excitations.py`'s `spectral_weights`/`_spectral_source_vector`) give each quasiparticle branch's exact delta-peak residue `|<k,a|O(k)|Psi>|^2`, i.e. a momentum-resolved dynamical structure factor directly in the thermodynamic limit — no finite window, unlike `kpm_finite`/`td_dynamical_correlator`. `return_total=True` adds the branch-complete total, which is exactly the connected static structure factor (a one-site operator applied to a uniform MPS lies exactly in the tangent space), so `weights.sum()/total` measures how much of the response the single-mode picture carries. Ported to v3 as `Chain::vumps_spectral_weights` (`gs_method="vumps"`, `n_uc<=2`): `vx_lanczos_lowest` now hands back its converged eigenvectors (and the dense branch the `zheev` ones), and the mixed-transfer source vector and the two deflated resolvents are dense-array ports of `_spectral_source_vector`/`_spectral_resolvent`; `vumps_excitation_energies` is untouched. The two backends agree to 9e-13 in energies, multiplet sums and totals on the dimerized Heisenberg cell at D=4, and to ~1e-11 on the TFIM at D=1 to 3 (set by the two ground states, not by the port). Validated four ways: the exact product-state limit, the static sum rule against an independent real-space `correlator` sum, the f-sum rule, and — sharpest — the S=1 **AKLT** point, whose ground state is exactly a D=2 MPS, where the closed-form structure factor and Arovas-Auerbach-Haldane's single-mode dispersion (PRL 60, 531 (1988), which the first moment must reproduce exactly) both come out to ~5e-15, and the SU(2) triplet/quintuplet split with its ΔS selection rule emerges rather than being imposed. On the pure Haldane chain the k=π gap converges monotonically to the literature 0.4104789 while the sum rule holds identically at every D. See `tests/test_infinite_chain_spectral.py`, `tests/test_infinite_chain_spectral_v3.py` (the same checks on v3 plus the cross-check against `"python"`), `examples/idmrg/dynamical_structure_factor_tfim/main.py` and `examples/idmrg/haldane_structure_factor/main.py`. |
| iDMRG local excitation gap (deflation-based, `D`-capable) | 🟡 | ✅ | ❌ | `local_excitation_gap(window=0)` runs on both backends (`idmrg.py::local_excitation_gap` and its C++ port `Chain::idmrg_local_excitation_gap`, `gs_method="idmrg"`). `window>0` (`local_excitation_gap_windowed`, a prototype) is pyitensor only. |
| iDMRG dynamical correlator (finite-window KPM reduction) | ✅ | ✅ | ❌ | `infinitechain.py::kpm_finite` builds its open-boundary window chain on the backend of the infinite chain (`itensor_version=3` on v3, `"python"` otherwise) and runs the ordinary finite KPM on it. Window spectra of the two backends agree in `tests/test_infinite_chain_algebra_v3.py`. |
| iDMRG real-time window dynamics (IBC-style TDVP, `td_dynamical_correlator`) | ✅ | ✅ | ❌ | Both wired (`idmrg_window.py` for pyitensor, `mpscpp3/chain_session.h`'s `td_dynamical_correlator_window` for v3); Phase 2+ generic sweep machinery beyond the current window construction is still pyitensor-only groundwork (`idmrg_window.py`'s own Phase 0-1 notes). **Fermionic (parity-odd) operator pairs are supported since 2026-08-29** on both backends — the Jordan-Wigner string is threaded across the window on the ket before the evolution and across the bra at measurement, pinned by the exact `S(x,t=0) == correlator(...)` identity and by the free-fermion Green function `<c†_x(t) c_0>` at `t>0` (`tests/test_idmrg_window_fermionic.py`, `examples/idmrg/fermionic_dynamical_correlator/main.py`); an odd-total-parity pair raises. **Four pre-existing bugs were found and fixed on the way**, all operator-independent. (1) The bra applied `A` where the contraction's own conjugation needs `A†` (silently `<A†_x ... B_0>`, invisible for Hermitian names). (2) Every `S(x,t)` carried a spurious global factor, because the window Hamiltonian's energy baseline is measured with the window's boundary legs *traced* while correlators close them with the transfer-matrix fixed points — now cancelled on **both** backends by co-evolving an unperturbed vacuum window and dividing by its `<ψ|ψ(t)>` through the identical contraction (the dimerized-XX free-fermion residual went from ~0.07, previously blamed on iDMRG convergence, to ~1e-5; v3's own `exp(+i*eshift*t)` estimate of the same factor, and with it the whole `eshift` measurement, is gone). (3) v3's window tiled the raw per-micro-step `idmrg_U_` factors instead of the gauge-consistent cell every other v3 static observable uses, missing the t=0 identity by up to 1.7e-1 even on a spin chain, on a model whose energy density the two backends agree on to 6.7e-11 — shape and finiteness looked fine throughout, which is how it survived; it now tiles the *raw* (un-re-gauged) cell `idmrg_cell_raw_`, whose outer legs are literally the environment caps, with the closures moved onto that cell's own left/right transfer fixed points plus a ground-state calibration denominator. (4) A latent index-stride bug in `idmrg_close_array_chain` indexed the bra's left axis with the bra's *right* extent — silent for as long as every tiled tensor was square, and wrong the moment it is not (a dimerized spinless chain at `maxm=20` converges to a (18,20)/(20,18) cell, where it put `S(x,0)` at 0.83 against an exact 0.51). Both backends are now exact on the t=0 identity to ~1e-15 and agree on `S(x,t)`; see `examples/idmrg/td_dynamical_correlator_python_VS_v3/main.py`, which evaluates it at *every* `x` — comparing at `x=0` alone, the one point a gauge error cannot touch, is what let (3) survive in the first place. |
| iDMRG cat-state superposition (`imps_sum` across distinct branches) | ❌ | ❌ | ❌ | Two branches with a genuine per-site norm mismatch sum fine, on both of pyitensor's ground-state representations (`idmrg.imps_sum` on `IDMRGResult`'s periodic `U_list`, and, since 2026-08-08, `vumps.imps_sum` on `VUMPSResult`'s mixed gauge `{AL,AR,C,AC}` — same construction and physical scope, ported per arXiv:1810.07006 Eq.(9)-(17)'s canonical-form algorithm, see `tests/test_vumps_imps_sum.py`/`examples/idmrg/vumps_imps_sum/main.py`); a true cat-state superposition of two physically distinct, EQUALLY-normalized symmetry-broken ground states (the degenerate-eta case both `imps_sum`s reliably raise `RuntimeError` on) needs new correlator machinery not yet written for any backend. |
| VUMPS `imps_sum` (norm-mismatch direct sum), ported to `itensor_version=3` | ✅ | ✅ | ❌ | `Chain::vumps_imps_sum` (`mpscpp3/chain_session.h`) is a dense-array port of `pyitensor/vumps.py`'s own `imps_sum` (block-diagonal direct sum of `this` Chain's own converged snapshot `AL` with a second, externally-supplied `AL_b`, re-canonicalized/truncated via the new `vx_canonicalize_n1` — the n_uc=1 specialization of `idmrg._canonicalize_periodic`, VUMPS always groups down to one supersite — then completed to the full mixed gauge via `vumps_complete_mixed_gauge`), reliably throwing for two ordinary equally-normalized states exactly like pyitensor's own degeneracy check. It is now reached through the public `Infinite_Many_Body_Chain.summps` on `itensor_version=3` with `gs_method="vumps"` (grouped cells only, the same scope and errors as `"python"`), which loads the result into a fresh Chain with `Chain::vumps_load_uniform_state` so `vev`/`correlator` read it; `gs_method="idmrg"` on v3 has no iMPS algebra, since its cell cannot be exported or loaded back. Cross-checked directly against `itensor_version="python"`'s own `vumps.imps_sum` at the exact D=1 field-polarized point and at genuinely entangled D=2,3 TFIM (dominant-branch `vev`/`correlator` matching to ~1e-6) — see `tests/test_vumps_imps_sum_v3.py`. |
| VUMPS `apply_mpo` (bounded periodic-MPO application), ported to `itensor_version=3` | ✅ | ✅ | ❌ | `Chain::vumps_apply_mpo` (`mpscpp3/chain_session.h`) is a dense-array port of `pyitensor/vumps.py`'s own `apply_mpo` (grow `this` Chain's own converged `AL` by a caller-supplied bounded MPO via the new `vx_grow_by_mpo_n1`, reusing `vumps_group_automaton` unmodified to group a multi-site `W_bulk` — the same grouping routine VUMPS's own Hamiltonian automaton uses, since the grouping contraction itself is agnostic to what it's grouping — then `vx_canonicalize_n1`/`vumps_complete_mixed_gauge`, shared with `vumps_imps_sum` immediately above). It is now reached through the public `Infinite_Many_Body_Chain.apply_operator` (and `overlap`) on `itensor_version=3` with `gs_method="vumps"`, grouped cells only, with the same bounded-operator restriction and `.eta` as pyitensor; `gs_method="idmrg"` on v3 raises NotImplementedError. Cross-checked directly against `itensor_version="python"`'s own `vumps.apply_mpo`: the exact D=1 field-polarized Pauli-X flip, unitary invariants (`<Sz>` flip, `<Sz(0)Sz(r)>` unchanged, `eta~1`) at genuinely entangled D=2,3 TFIM, and a genuinely bond-growing `chi_W>1` two-site gate tiled once per an n_uc=2 unit cell — see `tests/test_vumps_apply_mpo_v3.py` and `examples/idmrg/vumps_apply_mpo_v3_VS_python/main.py`. |

## 2. Expectation values, static correlators, entanglement

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| `vev` (operator expectation value) | ✅ | ✅ | ✅ | Julia takes `wf=`, `npow=` and `T=` with the semantics of `vev.py`'s `multi_vev` (the last factor an exact `<a|M|b>`, `T>0` raises naming `mode="ED"` and `metts_vev()`). |
| `gs_energy_fluctuation` / energy variance | ✅ | ✅ | ✅ | All three go through `vev.energy_variance`, `||(H-<H>)|psi>||^2` with the `maxm` cap lifted. |
| Chain-level MPS algebra (`overlap`, `aMb`, `applyoperator`, `summps`, `random_state(orthogonal=..)`) | ✅ | ✅ | ✅ | Julia states are recognized by `wavefunction_mode` through their `jlmps` attribute; `aMb` is the exact `inner(a',M,b)`. |
| MPO algebra (`toMPO`: `M*N`, `M+N`, scalars, `get_dagger()`, `trace()`) | ✅ | ✅ | ✅ | Julia: `mpsjulialive/mpo_ops.jl`, the product cut at `mpomaxm`. `M*N` is the operator product with `N` acting first on every backend, as for `MultiOperator` and `mode="ED"` (before 2026-11-08 v2, v3 and pyitensor gave `N*M`); `tests/test_julia_live_mpo_algebra.py` pins it on all four. |
| Two-point static correlators | ✅ | ✅ | ✅ | |
| All-pairs correlation matrix (native accelerated) | ✅ | ✅ | 🟡 | Julia's `dmmode="full"` is ITensorMPS's `correlation_matrix(psi,"Cdag","C")` (`mpsjulialive/correlations.jl`), on chains of spinless fermions only; other site types use the generic per-pair loop. |
| Four-point correlator, `ctmode="explicit"` (backend-agnostic loop) | ✅ | ✅ | ✅ | Always correct, always available, slowest — the universal fallback. |
| Four-point correlator, `ctmode="full"` (native AutoMPO per-tuple, incl. spinful-native sites) | ✅ | ✅ | 🟡 | v2 also has this. Julia (`gs_extra.jl::four_correlation_tensor`, one exact MPO per tuple) is the default there on chains of spinless fermions only; other site types raise and use `ctmode="explicit"`. |
| Four-point correlator, `ctmode="sweep"` (single-sweep, environment-reuse, ITensorCorrelators.jl-style) | ✅ | ✅ | ✅ | The fastest mode when it applies; no native-spinful-site counterpart yet even on v3/pyitensor. Julia: `mpsjulialive/fourpoint.jl`, spinless and interleaved spinful fermions, the default `ctmode` there. |
| Bond entanglement entropy / reduced density matrix | ✅ | ✅ | ✅ | Backend-agnostic dispatch (`entanglement.py` has no `itensor_version` branching at all — every backend's session exposes the same primitive). |
| Topological invariants (Berry phase, etc.) | ✅ | ✅ | ✅ | Built purely on `get_gs()` + MPS overlap — generic across every backend. |

## 3. Time evolution

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| Real-time evolution, TDVP (two-site) | ✅ | ✅ | ✅ | Julia has its own implementation (`mpsjulialive/timedependent.py` + `tdvp.jl`), and now honors `self.tevol_method` rather than ignoring it — see the TDVP_GSE row; anything outside `("TDVP","TDVP_GSE")` raises there instead of silently running plain TDVP. The local exponentiator is chosen by `chain.is_hermitian(h)`: Lanczos for a Hermitian H, Arnoldi otherwise. |
| Real-time evolution, TDVP + global subspace expansion (`tevol_method="TDVP_GSE"`) | ✅ | ✅ | ✅ | Julia needs no port: ITensorMPS.jl ships both pieces (`expand(psi,H; alg="global_krylov")`, citing the same Yang-White paper, and `tdvp(...; nsite=1)`), so `mpsjulialive/tdvp.jl` only wires them into the same expand-then-step structure. It handles a bond-dimension-1 (product-state) start fine — see `examples/time_evolution/tdvp_gse_julia_VS_ED_time_evolution` — and so does v3 since the 2026-09-24 third pass fixed its edge truncation (`docs/audit_2026_09_24c_hole_hunt.md`, finding 14). |
| Real-time evolution, one-site TDVP with controlled bond expansion (`tevol_method="TDVP_CBE"`) | ✅ | ✅ | ✅ | Julia in `mpsjulialive/cbe.jl` (KrylovKit Arnoldi for both half sweeps, so complex time steps work), wired into `evolve_and_measure`, `submode="TD"` and TDZ. |
| Real-time evolution, TEBD (2nd-order Trotter gates) | ✅ | ✅ | ✅ | v2 has no TDVP module at all so was never a candidate. Julia's `mpsjulialive/tebd.jl` builds each bond's gate from ITensor outer products (`op(name,s)`, no dense-matrix kron/reshape) and exponentiates via ITensors.jl's own tensor `exp()` (the same ITensor-native mechanism as `mpscpp3/tebd.h`'s `BondGate`, unlike pyitensor's manual `scipy.linalg.expm`). Unlike the other two backends, it resolves Jordan-Wigner from scratch off the *raw* `"C"/"Cdag"` term list (`mpo.py`'s `text_mpo()`), not `MultiOperator.to_terms()`'s pre-dressed `"A"/"Adag"/"F"` form — real ITensors.jl's builtin `"Fermion"` site type only defines `op("C",..)`/`op("Cdag",..)`/`op("F",..)`, not dmrgpy's own `"A"/"Adag"` names. |
| Real-time evolution, MPO-Taylor (`tevol_method="MPO"`, the pre-TDVP legacy path) | ✅ | — | — | The only option for v2; still selectable on v3 by explicitly setting `tevol_method="MPO"`. |
| Complex-time evolution / TDZ damping (arXiv:2311.10909) | ✅ | ✅ | ✅ | Julia needed exactly one new primitive (`advance_complex_time_step`); the rest of `tdz.py` is backend-agnostic MPS/MPO algebra. |

## 4. Dynamical correlators (`get_dynamical_correlator(submode=...)`)

| submode | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| `KPM` (Kernel Polynomial Method, Chebyshev) | ✅ | ✅ | ✅ | Julia takes the shared dispatch of `dynamics.py`/`kpmdmrg.py`; only its moments are its own (`mpsjulialive/dynamics.py::kpm_moments` + `kpm.jl`), so `get_dynamical_correlator_moments`, `get_gr` and `get_distribution` work there too. |
| `KPM` energy truncation (Holzner Sec. III-B) | ✅ | ✅ | ❌ | Net *slower* on both backends where it's measured, kept for regimes/paper-fidelity, not a default. Julia raises `NotImplementedError` naming the backends that have it. |
| Non-Hermitian KPM (NHKPM.jl-ported coupled recursion) | ✅ | ✅ | ✅ | Auto-selected instead of plain KPM whenever the Hamiltonian is non-Hermitian, for backends that support it; ED also has this path. Julia's recursion is `mpsjulialive/nhkpm.jl`, on the left state the NH solve keeps (`nh_left_wf`). |
| `CVM` (correction-vector / conjugate-gradient resolvent) | ✅ | ✅ | ✅ | Julia's CVM needed no new primitive beyond `set_sweep_params`-equivalent bookkeeping — already backend-agnostic MPS/MPO algebra. |
| `CVM` variational solver (`cvm_solver="variational"`, Jeckelmann dynamical DMRG, A=B^dagger) | ✅ | ✅ | ✅ | Two-site sweep minimizing `<x|M|x>-2Re<b|x>` with `M=(H-E0-omega)^2+eta^2` as a squared MPO: `pyitensor/ddmrg.py` and, on v3, `Chain::ddmrg_correction_vector` (`mpscpp3/chain_session.h`, local CG on the two-site block, solver scale `hscale_up_` applied to H, omega, eta and b alike). v2 and Julia use the global CG. Cross-checked against `"python"` and ED in `tests/test_ddmrg_v3.py`. Julia: `mpsjulialive/ddmrg.jl`, the same two-site sweep, on the chain's own ground state. |
| `CVM_explicit`, `CVMimag` (analytic continuation variant) | ✅ | ✅ | ✅ | Julia through the `applyinverse` of `bicgstab.jl`, which also gives the non-Hermitian `CVM`. `CVMimag`'s continuation to the real axis is ill-conditioned on every MPS backend (about 0.1 of the peak); its resolvent grid matches ED. |
| `TDZ` (complex-time evolution correlator) | ✅ | ✅ | ✅ | |
| `TD` (real-time evolution correlator) | ✅ | ✅ | ✅ | Dispatches through `timedependent.py`, which has its own Julia branch. |
| `ROOTN` (root-N Krylov correction-vector) | ✅ | ✅ | ✅ | Backend-agnostic MPS algebra (`rootndmrg.py`), reached on Julia through the shared dispatch. The per-bond-sweep shortcut is confirmed impossible on both C++ and pyitensor (needs real multi-target state-averaging). |
| `SECTOR` (Lehmann sum over per-sector eigenstates, `get_spectral_function`, `get_spin_spectral_function`) | ✅ | ✅ | ✅ | Julia solves the sectors on a clone that shares the caller's dense site set, and assembles the poles with the Julia MPS's own `dot`/`aMb`. Checked against the restricted ED Lehmann sum and the spin sum rules in `tests/test_sector_spectral_function.py`. |
| `TASK` (tangent-space Krylov poles, arXiv:2510.07279) | ✅ | ✅ | ✅ | Julia exports the ground-state MPS and the MPO site tensors as arrays (`mpsjulialive/arrays.jl`) and runs `pyitensor/task.py` on them, as v3 does; `pad_bonds`, `expand_bonds` and `missing_weight` are the same code. |
| Third-order Kondo spectrum, DMRG route (`get_kondo_spectrum(order=3, mode="DMRG")`) | ✅ | ✅ | ✅ | Julia's two-time construction calls `tdvp.jl`'s `tdvp_step`; the default `submode` is `"TASK"` on all three. Checked against ED in `tests/test_kondo_spectrum_dmrgtwotime.py` and `tests/test_kondo_spectrum_task.py`. |
| `EX` (exact spectral decomposition over excited states) | ✅ | ✅ | ✅ | Non-Hermitian-capable; entirely backend-agnostic (`self.get_excited_states()` + generic MPS algebra). |
| `maxent` (max-entropy positive-definite reconstruction) | ✅ | ✅ | ✅ | Same — backend-agnostic post-processing over `power_vev` moments. |

## 5. Finite-temperature methods

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| Ancilla purification + imaginary-time annealing (`thermal.py`) | ✅ | ✅ | ✅ | Backend-agnostic — just wraps a doubled-site `Spin_Chain` and calls its generic `get_gs`/evolve methods. |
| Exact Boltzmann sum over ED excited states (`thermalvev.py`) | — | — | — | ED-only by construction (needs the full spectrum); not a DMRG-backend feature. |
| METTS thermal average (`metts_vev`) | ✅ | ✅ | ✅ | Needs imaginary-time TDVP, so v2 (no `TDVP/`) is excluded. Julia's `mpsjulialive/metts.jl` is a value-level port of pyitensor/metts.py's algorithm (per-site eigendecomposition + sequential-sampling collapse, both via plain Julia arrays rather than ITensor-level primitives) that reuses `tdvp.jl`'s `tdvp_step` unchanged for the imaginary-time evolution -- no new evolution primitive needed, only the sampling/collapse driver. `njobs>1` multiprocess pooling is pyitensor-only (no equivalent lever for a single live v3/Julia session -- `julia_live` raises for `njobs>1`, same as v3). |
| METTS dynamical correlator (finite-T real-time correlator) | ✅ | ✅ | ✅ | Julia's `mpsjulialive/metts.jl`'s `metts_dynamical_correlator` reuses the same METTS Markov chain and `tdvp.jl`'s `tdvp_step` (now with a real `dt`, generalizing the same `-im*dt` convention `metts_vev` already exploits with an imaginary one) -- no new evolution primitive needed, and (unlike the C++ port) no manual norm-restoration workaround either, since this backend's `tdvp_step` already passes `normalize=false` for every caller. `tdvp_cutoff`/`tdvp_maxdim` are wired through here too (Julia's `tdvp_step` takes them per call, unlike v3's hardcoded restriction to this chain's own cutoff/maxm); `njobs>1` remains python-only. |

## 6. Physical models (Hilbert-space coverage)

Model coverage is largely backend-independent — `spinchain.py`,
`fermionchain.py`, `spinfermionchain.py`, `bosonchain.py`,
`parafermionchain.py` all build their operators as backend-agnostic
`MultiOperator`s, so any new model works with any backend that supports
the *method* being called on it (per the tables above). The
model-specific rows:

| Capability | v3 | pyitensor | Julia | Notes |
|---|---|---|---|---|
| Native spinful-fermion sites (`Spinful_Fermionic_Chain_Native`, real `ElectronSite`) | ✅ | ✅ | ✅ | A genuine per-site `Cup/Cdn/Cdagup/Cdagdn` site type: v3's `HubbardSite`, pyitensor's `ElectronSite` (`pyitensor/sites/fermion.py`) and ITensors.jl's `"Electron"` site, basis `Emp, Up, Dn, UpDn`. Julia gets the raw names and ITensorMPS threads the Jordan-Wigner strings, up before down on a site. Slower for iterative GS/TDVP/KPM than the interleaved `Spinful_Fermionic_Chain`, faster for the 4-point tensor — a real perf tradeoff, not a strict upgrade, so both representations stay available. |
| Higher spins (S=3/2, 2, 5/2), Z2/Z3/Z4 parafermions, bosons of any `maxnb` | ✅ | ✅ | ✅ | Julia builds every site code of `pyitensor/sites/siteset.py`; `mpsjulialive/sites_custom.jl` defines what ITensors.jl does not ship (the higher spins, the clock `Sig`/`Tau`, the boson projectors, `Sp/Sm/ISy`). An unknown code raises naming it. See `tests/test_julia_live_site_types.py`, `tests/test_boson_chain.py`, `tests/test_parafermion_chain.py`. |
| Mixed spin-fermion chain (`Mixed_Spin_Fermion_Chain`) | ✅ | ✅ | ✅ | v2 is refused (it aborts on a string crossing a spin site). On Julia ITensorMPS puts an `F` on every site a fermionic term crosses, and `sites_custom.jl` makes `F` the identity on spin sites. A conserved sector on a chain with both kinds is ED only on every backend. |
| Parity-odd fermion operators as MPOs (a bare `C_i`, Green's functions) | ✅ | ✅ | ✅ | Julia's OpSum refuses parity-odd fermionic terms, so `mpsjulialive/mpo.py::get_MPO` sends them in the Jordan-Wigner form of `MultiOperator.to_terms()`, with `A/Adag/F` defined as plain local matrices in `fermion_ops.jl`. Parity-even operators keep the raw names. |

## What's missing, roughly ranked by how load-bearing it'd be

1. **Julia: what is still missing.** The finite-chain surface of
   `julia_live` is now close to v3 and pyitensor (the rows above). What
   it still lacks:
   - **Infinite chains beyond VUMPS:** `gs_method="idmrg"`, the
     excitation ansatz, spectral weights, the finite windows and the iMPS
     algebra still raise on julia_live, and its VUMPS covers `n_uc<=2`
     with nearest-cell couplings only.
   - **Deferred:** `get_krylov_spectrum` (CLIK), `Lindblad_Chain`, the
     Holzner KPM energy truncation (`kpm_energy_truncate`), YSR
     (`ysr.py`), and the rest of the
     infinite-chain stack: iDMRG, the excitation ansatz, spectral
     weights, window dynamics and the iMPS algebra
     (`infinitechain.py` rejects `julia_live` at construction).
     iDMRG in particular is a large port, worth it only for a concrete
     use case that needs Julia on an infinite system.
   - **Site-type limits:** `dmmode="full"` and `ctmode="full"` cover
     spinless fermions only, and conserved sectors cover spinless
     fermions, `Electron`, S=1/2, S=1 and bosons, not the higher spins or
     the parafermions.
2. **Julia: real time and finite temperature.** Done. TDVP, TDVP_GSE,
   TEBD (`mpsjulialive/tebd.jl`, Jordan-Wigner threaded off the raw
   `"C"/"Cdag"` term list, gates exponentiated by ITensors.jl's own
   `exp()`), TDZ and METTS (`metts_vev` and `metts_dynamical_correlator`,
   `mpsjulialive/metts.jl`) all run there. TEBD on v2 was never a
   candidate, since v2 has no TDVP module at all.
3. **v3/pyitensor: four-point `ctmode="sweep"` for native-spinful sites.**
   The fast environment-reuse algorithm has no flavor-resolved
   counterpart yet; `ctmode="full"` remains the practical choice for
   `Spinful_Fermionic_Chain_Native`.
4. **iDMRG excitation ansatz beyond `D=1`.** Now implemented — see
   `pyitensor/idmrg_excitations.py`'s own module docstring ("History"
   section) for the full account of what it took: after eight
   investigation passes patching the old, single-uniform-gauge-tensor
   diagram list (documented in that module's own git history and
   `docs/idmrg_excitation_mpskit_port_plan.md`) failed to fix `D>1`, the
   ansatz was rewritten from scratch on top of a VUMPS ground state's own
   mixed-gauge `{AL,AR,C,GL,GR}`, mirroring MPSKit.jl's actual
   architecture (channel-resolved `GBL(k)`/`GBR(k)`, a 3-term `H_eff(k)`)
   rather than its published equations alone. This also surfaced and
   fixed two further, previously-undiscovered bugs: a real, pre-existing
   missing-conjugate sign bug in `vumps.py` itself (`_solve_left_
   environment`/`_solve_right_environment`/`_energy_density_and_source_
   from_left`/`_right`, invisible at `D=1` where the relevant fixed point
   is always real — this also measurably improves `gs_method="vumps"`'s
   own `D>1` ground-state convergence reliability, previously documented
   as a known, unresolved robustness gap in `vumps.py`'s own "Convergence
   robustness" section), and the correct excitation-energy renormalization
   constant (`H_AC`'s own eigenvalue `lambda_AC`, not the ground state's
   physical energy density `e_cell`). The `D=2` transverse-field Ising
   dispersion now matches an independently-converged MPSKit.jl result to
   6 significant figures, and `H_eff(k)` is Hermitian to machine precision
   at `D=1`, `2` and `3` — see `examples/idmrg/excitation_gap_tfim/
   main.py`. **Update:** both `vumps_ground_state` and this excitation
   ansatz are now also ported to `itensor_version=3` (`Chain::
   vumps_ground_state`/`Chain::vumps_excitation_energies`,
   `mpscpp3/chain_session.h`) — see item 7 below.
5. **iDMRG cat-state superpositions.** No backend supports summing two
   *physically distinct, equally-normalized* symmetry-broken iDMRG ground
   states; needs new correlator machinery from scratch. (Summing two
   states with a genuine per-site norm mismatch is supported, on both
   pyitensor ground-state representations — `idmrg.imps_sum` and, since
   2026-08-08, `vumps.imps_sum` — see item 7's own row in the matrix
   above.)
6. **v2 (legacy) parity.** Deliberately not being chased — v2 predates
   `TDVP/` entirely and is kept mainly as the historical QN-conserving
   reference implementation, not a target for new features.
7. **iDMRG VUMPS ground state + excitation ansatz, ported to
   `itensor_version=3`.** `Chain::vumps_ground_state`/
   `Chain::vumps_excitation_energies` (`mpscpp3/chain_session.h`) close
   the gap items 1 and 4's own table rows used to flag (v3 previously had
   neither `gs_method="vumps"` nor the tangent-space excitation ansatz at
   all). Unlike `Chain::idmrg_ground_state` (built from genuine ITensor
   tensor-network objects, mirroring pyitensor's own `IDMRGResult`
   line-for-line), this port is deliberately built from plain dense
   row-major `std::vector<Cplx>` arrays closed over LAPACK
   (`itensor::zgeev_wrapper`/`zheev_wrapper`/`zgesv_wrapper`/
   `zgesvd_wrapper`) rather than ITensor `Index`/`ITensor` objects — see
   `chain_session.h`'s own comment on `VumpsResult` for why: `D` and
   `d_g` (the grouped supersite's own physical dimension) are always
   small in this feature's scope (`n_uc<=2`, reach<=1 bonds), so exact
   dense linear algebra is both simpler and far less risky than
   re-deriving VUMPS's/the excitation ansatz's already extremely subtle
   fixed-point and gauge bookkeeping (see item 4 above — eight
   independent investigation passes were needed to get this right in
   pure Python the first time) against ITensor v3's `Index`/`IndexSet`
   machinery a second time. One deliberate simplification relative to
   pyitensor's own `vumps_ground_state` driver: the D-ramp/multi-restart
   strategy is ported, but not the "spend extra attempts to beat an
   already-known smaller-D energy" safety-net budget (see
   `Chain::vumps_ground_state`'s own comment) — a caller doing
   quantitative work at a harder/larger `D` should still call it several
   times independently and keep the lowest reported `e0`, exactly as
   pyitensor's own docstring already recommends.

   Validated directly against `itensor_version="python"`'s own
   already-MPSKit.jl-validated implementation (not just internal
   self-consistency, per this project's own literature-check standing
   feedback): ground-state energy density matches to ~1e-10 or tighter at
   `D=1,2,3` on TFIM (`n_uc=1`) and the uniform Heisenberg chain
   (`n_uc=2`); the excitation dispersion matches across a full momentum
   scan to the same tolerance on the gapped TFIM case, and to ~1e-4..1e-7
   on the gapless/critical Heisenberg case (both backends' non-convex
   restart searches can land on slightly different local optima there —
   not a discrepancy in the ported algorithm itself). See
   `tests/test_vumps_v3.py`/`tests/test_vumps_excitations_v3.py` and
   `examples/idmrg/vumps_excitation_v3_VS_python/main.py`.

   **A genuine, previously-undetected bug this audit also found and
   fixed in the ALREADY-SHIPPED `Chain::idmrg_ground_state`** (not new
   code this item added): `idmrg_build_row` was putting a site's own
   onsite/field-operator term on the automaton's `F,F` self-loop transition
   instead of the direct `S,F` transition — exactly the bug
   `pyitensor/idmrg.py`'s own `_build_periodic_mpo` docstring documents
   finding and fixing on the Python side (see that docstring's own two
   confirmed reproducers: a bare-field Hamiltonian with no bond terms
   silently drops the field entirely, since `W` stays block-diagonal
   between `{S}` and `{F,pending...}`; a Hamiltonian with both bond and
   field terms instead diverges exponentially in the number of macro-
   iterations, since every further-absorbed site while `F` is already hot
   re-adds the onsite content on top of an already-summed total). This
   was never caught by the existing test suite because no
   `itensor_version=3` `idmrg_ground_state` test exercises an onsite/field
   term (every `tests/test_infinite_chain.py` v3 case is a pure
   Heisenberg/XXX-style bond-only Hamiltonian) — confirmed directly while
   building this item's own VUMPS test cases (a TFIM Hamiltonian, which
   does have a field term, is the very first thing this port was
   validated against). Fixed at the shared root cause (`idmrg_build_row`,
   used by both `idmrg_ground_state` and the new `vumps_ground_state`),
   not patched per call site.

8. **iDMRG VUMPS static correlators, ported to `itensor_version=3`.**
   `Chain::vumps_onsite_expectation`/`Chain::vumps_two_point_correlator`
   (`mpscpp3/chain_session.h`) close the correlator half of the gap item
   7 above left open (that item ported the ground-state solver and the
   excitation ansatz, but explicitly not `vev`/`correlator`). At the time,
   `gs_method="idmrg"` still had no correlator support at all on this
   backend; item 9 below closes that half too. Same deliberate architecture as item 7's own
   port: a line-for-line dense-array translation of `pyitensor/vumps.py`'s
   `onsite_expectation`/`two_point_correlator` (`AC`'s exact normalization
   and `AR`'s exact right-orthonormality mean neither the left nor the
   right closure needs an eigensolve), plus one new private helper,
   `Chain::vumps_embed_group_operator` (the C++ analogue of
   `_embed_group_operator`'s `np.kron`-based embedding, reusing the same
   per-sublattice dense operator matrices `idmrg_op_dense` already
   provides for `idmrg_ground_state`'s own automaton construction).
   `Infinite_Many_Body_Chain.vev`/`correlator` now route to it whenever
   `itensor_version=3` and `self.gs_method=="vumps"`. Validated directly
   against `itensor_version="python"`'s own VUMPS correlators: the same
   exact `D=1` closed-form cases both backends already had (a
   field-polarized product state, a decoupled Heisenberg singlet dimer),
   plus a direct cross-check against `itensor_version="python"` on TFIM
   at `D=2,3`, matching to ~1e-14 or tighter across the ground-state
   energy, `vev`, and a multi-`r` correlator scan — see
   `tests/test_vumps_correlator_v3.py` and
   `examples/idmrg/vumps_correlator_v3_VS_python/main.py`.

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

9. **iDMRG growing-algorithm static correlators and local gap, ported to
   `itensor_version=3`.** `Chain::idmrg_onsite_expectation`/
   `Chain::idmrg_two_point_correlator`/`Chain::idmrg_local_excitation_gap`
   (`mpscpp3/chain_session.h`) close the last piece item 8 above left
   open, and with it the whole `gs_method="idmrg"` half of the v3 iDMRG
   surface. They rest on three prerequisites from `pyitensor/idmrg.py`
   that the original C++ port never got and that this item ported first:
   McCulloch's wavefunction prediction as each micro-step's Krylov start
   (`idmrg_wavefunction_prediction`), extraction of the gauge-consistent
   unit cell from a single micro-step's theta (`idmrg_theta_cell` +
   `ic_canonicalize_cell`), and the per-site energy baseline subtracted
   from both growing environments (`idmrg_subtract_energy_baseline`) —
   the unit cell in particular is a hard prerequisite, since tiling the
   raw per-micro-step factors instead leaves the energy correct while
   getting correlators wrong by orders of magnitude (and, on the XX
   chain, by a sign). The observables themselves are a dense-array
   translation of `idmrg.py`'s transfer-matrix machinery (`ic_*`, same
   architecture as items 7/8's `vx_*`), with two performance departures
   that do not change any returned number: the measurement-independent
   fixed points are cached per ground-state run instead of resolved per
   call, and an operator string is applied to the fixed point site by
   site (O(chi^4)) instead of composed as transfer matrices (O(chi^6)).
   `td_dynamical_correlator_window`'s own connected-background
   subtraction was switched onto the same cell-based expectation value,
   fixing a latent gauge bug there. Validated against exact closed-form
   values (a field-polarized product state; the `<H_uc> = n_uc*e0`
   translational-invariance identity, which is what actually guards the
   gauge), against `itensor_version="python"` on gapped spin and
   fermionic chains including Jordan-Wigner-strung correlators, and on
   the local gap — see `tests/test_idmrg_correlator_v3.py`,
   `tests/test_infinite_chain.py::test_itensor_version3_matches_python_backend`
   (which now compares the converged *state*, not just `e0`), and
   `examples/idmrg/idmrg_correlator_python_VS_v3/main.py`. Still
   `"python"`-only: `local_excitation_gap(window>0)` (a prototype, not
   stable API) and the iMPS algebra over `IDMRGResult`
   (`imps_overlap`/`apply_mpo`/`imps_sum`), exposed as
   `overlap`/`apply_operator`/`summps` for `gs_method="idmrg"` on
   `"python"` only. On v3 the same three public methods run for
   `gs_method="vumps"` (`tests/test_infinite_chain_algebra_v3.py`).

10. **The excitation ansatz's solver improvements, mirrored on
   `itensor_version=3`.** Done, both halves, on 2026-09-22.
   `pyitensor/idmrg_excitations.py` caches the momentum-dependent channel
   resolvents per momentum instead of rebuilding them inside every
   `_h_eff_action` call (they depend only on `k` and the
   momentum-independent environment, never on the excitation tensor `B`),
   and solves `H_eff(k)` by Lanczos on its action rather than by assembling
   it once the eigenproblem exceeds `_DENSE_EIG_MAX`;
   `mpscpp3/chain_session.h` now has both.
   `Chain::vumps_exc_resolvents` holds the two resolvents of the momentum
   most recently asked for, and `Chain::vx_resolvent_build`/
   `vx_resolvent_solve` keep that map's LU factorization, so the hundreds
   of solves one eigensolve asks for share one factorization instead of
   rebuilding and refactorizing a `D^2`-by-`D^2` map per application (the
   factorization is hand-rolled, `vx_lu_factor`/`vx_lu_solve`: ITensor
   wraps `zgesv` only, which factors and solves in one call and keeps its
   pivots to itself). `Chain::vx_lanczos_lowest` then solves `H_eff(k)` on
   its action above `vumps_h_eff_dense_max_`, falling back to the dense
   path whenever it cannot vouch for what it found. Measured over a
   3-momentum scan, threads pinned to one core, before against after: a
   `D=16` TFIM chain went from 164.1 s to 1.8 s and an `n_uc=2`
   Heisenberg chain at `D=10` from 15.6 s to 1.2 s, the resolvent cache
   carrying 164.1 s to 10.3 s and 15.6 s to 5.4 s of that and the Lanczos
   solver the rest. Forcing the two solvers against each other on one
   converged state, which is the only way to price the second half by
   itself, it is 5.3x at `dim=256` and 6.0x at `dim=432`. Numbers are
   unchanged: dense against forced Lanczos on the same state agrees to
   4.6e-11 over `D=2,3,4` at `n=1,2,3` on both models, the new build
   differs from the pre-change one by at most 1.3e-10 on TFIM at
   `D=1,2,3` against that build's own run-to-run scatter of 1.2e-10, and
   the agreement with `itensor_version="python"` over the same scan is
   2.7e-10, the figure the pre-change build gave too. At `D=16` the
   converged state itself moves by ~3e-5 from run to run, so that point
   prices the change and does not check it.
   Two things there are worth knowing before touching it. The threshold is
   64, not the 256 `_DENSE_EIG_MAX` uses, because the crossover is not
   transferable between the two backends and was measured again here: one
   application of `H_eff(k)` solves four channel resolvents, so Lanczos
   already wins at `dim=36`, the smallest size tried (see the measurement
   table at `vumps_h_eff_dense_max_`), and 64 is above that while keeping
   every model the port was validated on, `dim<=48`, on the exact dense
   path. And a single-vector Lanczos cannot return a degenerate eigenvalue
   twice, which the `n_uc=2` Heisenberg cell does away from `k=0`: asked
   for the lowest three at `k=0.37`, `D=2`, it returned three distinct
   eigenvalues where two of them should have been the same one, an error
   of 1.0, with every returned value a genuine eigenpair that no residual
   test rejects. `vx_lanczos_lowest` therefore runs one deflated Lanczos
   per eigenvalue, each from its own generic start vector, and
   `tests/test_vumps_excitations_v3.py` pins the pair of solvers against
   each other on exactly that model. The `"python"` side, which asked
   ARPACK for all `n` pairs at once from one constant start vector, was
   recorded here as not exposed, on a measurement at `D=2` (dimension
   12) that could not discriminate, since there the Krylov basis is the
   whole space. It was exposed above its dense threshold: on the same
   cell at `maxm=10` (dimension 300) the single call was wrong in 7 of 16
   calls over four momenta and `n=1` to 4, dropping one copy of a
   degenerate level. Since the 2026-09-24 audit
   (`docs/audit_2026_09_24_hole_hunt.md`, finding 15) it runs one deflated
   Lanczos per value for `n>=2`, a port of `vx_lanczos_lowest`, and is 0
   of 16 wrong there, so the deflation was a fix owed to both backends.
   See `docs/idmrg_improvement_plan.md`, which also ranks what else is
   worth doing on infinite chains.

After adding a capability to a backend, update this file's matrix and,
if it's a new physics-facing method, its page under `docs/guide/` per
`CLAUDE.md`'s documentation policy.
