---
paths:
  - "src/dmrgpy/multioperator.py"
  - "src/dmrgpy/multioperatortk/canonical.py"
  - "src/dmrgpy/nhdmrg.py"
  - "src/dmrgpy/nhlevels.py"
  - "src/dmrgpy/cvm.py"
  - "src/dmrgpy/rootndmrg.py"
  - "src/dmrgpy/kpmdmrg.py"
  - "src/dmrgpy/nrgmps.py"
  - "src/dmrgpy/algebra/**"
  - "src/dmrgpy/pyitensor/mpobuilder.py"
  - "src/dmrgpy/pyitensor/nhdmrg.py"
  - "src/dmrgpy/mpscpp2/mo_terms.h"
  - "src/dmrgpy/mpscpp3/mo_terms.h"
  - "src/dmrgpy/mpscpp2/chain_session.h"
  - "src/dmrgpy/mpscpp3/chain_session.h"
  - "tests/test_audit_2026_09_25_scale.py"
  - "tests/test_audit_2026_09_25b_cvm.py"
  - "tests/test_audit_2026_09_25b_nh.py"
  - "tests/test_audit_2026_09_25b_cpp.py"
  - "tests/test_audit_2026_09_25b_ednum.py"
  - "tests/test_audit_2026_09_25b_mpobuilder.py"
  - "tests/test_audit_2026_11_21_mpo_cutoff.py"
  - "tests/test_audit_2026_11_21_nrg.py"
---

# Thresholds are relative to a scale of the problem

These are durable rules the hole hunts left, with the mechanism and the
numbers; CLAUDE.md keeps one line per rule. Each rule names the test that
pins it. Do not undo one without reading its finding in `docs/audit_*.md`
first. The file loads when a file that holds one of these thresholds is read or edited.

- **Every threshold is relative to a scale of the problem it acts on.**
  Where a rule exists it is: multiply by the power of two that brings a
  largest coefficient below 1 into [1,2), and leave everything
  byte-identical at 1 or more. It covers v2/v3's `to_mpo_unit`/
  `solver_hamiltonian` (`hscale_up_`, read from the merged AutoMPO),
  `pyitensor`'s `to_mpo`, the NH-DMRG entry in `nhdmrg.py` (once, for all
  three session backends, so do NOT also scale inside `Chain::nhdmrg`) and
  ED's deflated eigensolver. The others: `multioperator.clean_threshold` is
  1e-12 relative to the largest coefficient (in `canonical_dict`, to the
  larger of that and the summed signature magnitudes); CVM stops relative
  to `||b||` and its early exits read the CG functional; NH-DMRG's Ritz tie
  window is `1e-6*(max Re - min Re) + 100*eps*max|ev|` with `<=`,
  identically in `pyitensor/nhdmrg.py`, in both `arnoldi_select_kbest`
  and `nh_ritz_order`, and in `nhlevels._tie_window` over the Ritz values
  of the imaginary-time projection; the VUMPS Lanczos reference is capped
  by the Hamiltonian's unit; the KPM same-vector shortcut and ROOTN's
  Lanczos breakdown compare against the vectors' own norms
  (`test_audit_2026_09_25_scale.py`,
  `test_audit_2026_09_25b_{cvm,nh,cpp,ednum,mpobuilder}.py`). The NRG
  tail of `nrgmps` has the same number of sites at every band D
  (`max_sites`, `OMEGA_MIN` relative to D), and its hand-over check is
  1e-12 of the largest hopping of the head (`test_audit_2026_11_21_nrg.py`).
  A weak term is covered by a cutoff of the MPO build that is one per bond
  (`mo_terms.h`'s `bond_cutoff_args` on v2/v3, `mpobuilder._bond_cutoff_factor`
  on `"python"`). The cutoff is lowered by m^2, floored at (1e-14)^2. m is
  the smaller of two ratios to the largest coefficient: the strongest
  coefficient crossing the weakest bond, and the weakest multi-site
  coefficient below 1e-5. A term's span starts after the Jordan-Wigner
  strings that cancel in pairs, and an operator with m = 1 is built
  exactly as before. Do NOT edit the vendored svdMPO for this
  (`test_audit_2026_09_25_scale.py`, `test_audit_2026_11_21_mpo_cutoff.py`).
