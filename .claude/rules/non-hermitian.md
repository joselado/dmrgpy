---
paths:
  - "src/dmrgpy/nhdmrg.py"
  - "src/dmrgpy/nhlevels.py"
  - "src/dmrgpy/nonhermitian/**"
  - "src/dmrgpy/pyitensor/nhdmrg.py"
  - "src/dmrgpy/mpsjulialive/nhdmrg.py"
  - "src/dmrgpy/mpsjulialive/nhdmrg.jl"
  - "src/dmrgpy/algebra/algebra.py"
  - "tests/test_audit_*nh*.py"
  - "tests/test_nh_*.py"
---

# Non-Hermitian Hamiltonians

These are durable rules the hole hunts left, with the mechanism and the
numbers; CLAUDE.md keeps one line per rule. Each rule names the test that
pins it. Do not undo one without reading its finding in `docs/audit_*.md`
first. The file loads when `nhdmrg.py`, `nhlevels.py`, `pyitensor/nhdmrg.py` or a file under
`nonhermitian/` is read or edited. The Ritz tie window is in `.claude/rules/thresholds.md`.

- **NH-DMRG sweeps the transposed left state.** Inside every v2, v3 and
  `"python"` sweep the left state is ITensorNHDMRG.jl's own,
  H^T|wl> = E|wl>, solved against the conjugate of the adjoint MPO and
  conjugated into dmrgpy's adjoint psil on the way out, and each bond is
  solved by a Krylov-Schur restart that keeps 3/5 of the Ritz vectors;
  `restarts` counts thick restarts. With the adjoint left state the
  shared isometry has to hold r and conj(r) of a complex-symmetric H, and
  on the 10-site PT-symmetric XX chain every attempt stopped at a
  non-eigenpair, Re E up to 0.03 below the spectrum
  (`test_audit_2026_11_15_nh.py`).
- **A non-Hermitian `gs_energy()` on DMRG is imaginary time first.**
  `nh_gs_method` (a setting, in `method_key`) defaults to
  `"imaginary_time"`: `nhlevels.py` on v3 and `"python"` without a
  sector, its pair checked by `nhdmrg()`'s own certificate and stored by
  `nhdmrg.store_nh_pair`; a failed check runs `nhdmrg()` and keeps the
  pair of smaller residual, a raised route and every other backend go to
  `gs_energy_nhdmrg`; `mode="ED"` is untouched. The pair taken is the
  smallest real part among the accepted Ritz pairs below the certificate's
  tolerance, and only a tie inside the Ritz tie window goes by residual
  (`test_audit_2026_11_21_nh.py`). Its
  grid is in units of the largest coefficient off the identity, since an
  offset as the unit leaves nothing projected. `get_nh_levels()` returns
  the accepted Ritz pairs (`test_nh_levels.py`).
- **ED orders a non-Hermitian spectrum by real part first**:
  `algebra.eigen_order`, behind `sorteigen`, groups real parts at 1e-10 of
  the largest |E| and orders by Im only inside a group;
  `biorthogonal_ground_state` takes its left state by its eigenvalue, never
  by a sort of `dagger(h)` (`test_audit_2026_11_21_nh.py`).
