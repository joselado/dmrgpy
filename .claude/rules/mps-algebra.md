---
paths:
  - "src/dmrgpy/mpsalgebra.py"
  - "src/dmrgpy/mpsalgebratk/**"
  - "src/dmrgpy/pyitensor/mpsalgebra.py"
  - "src/dmrgpy/algebra/inverse.py"
  - "src/dmrgpy/mpsjulialive/bicgstab.jl"
  - "src/dmrgpy/mpsjulialive/mpsalgebra.jl"
  - "tests/test_bicgstab*.py"
  - "tests/test_audit_2026_11_21_exact-arith.py"
---

# MPS algebra: the linear solves

These are durable rules the hole hunts left, with the mechanism and the
numbers; CLAUDE.md keeps one line per rule. Each rule names the test that
pins it. Do not undo one without reading its finding in `docs/audit_*.md`
first. The file loads when `mpsalgebra.py` or a file under `mpsalgebratk/` is read or edited.

- **A BiCGSTAB solve checks its answer.** After BiCGSTAB, `applyinverse`
  on DMRG measures ||Ax-b||/||b|| once with the maxm cap lifted
  (`mpsalgebra.check_bicgstab`), warns above 10 delta and raises above 1;
  before, v3 returned a vector 3e64 times the right-hand side on 10 sites
  at maxm=12, silently. Do NOT loosen it to quiet a test; the divergence
  of the recurrence itself is still open
  (`test_bicgstab_residual_check.py`). On `"python"`, v3 without a
  conserved sector and `julia_live` the check forms A x and A x - b with
  nothing truncated (`residual_norm`, `mps_residual_norm`); the default
  density-matrix routes overstate a converged residual next to a level up
  to 2e6 times. A sum whose coefficients cancel passes
  `sum_mps_weighted(method="exact")`, as `nhlevels._combine` does; do NOT
  make it the default, the ESR Krylov route sums hundreds of states
  (`test_audit_2026_11_21_exact-arith.py`). The recurrence itself applies
  and sums exactly, truncating once per step, on `"python"` and on v3
  without a sector (the author's pick for v3): next to a level it converges
  where the density-matrix routes stopped at 0.62, and where `maxm` limits
  the solve it costs 2 to 4 times more and diverges in the same cases, so
  the divergence is not from the route of the apply.
