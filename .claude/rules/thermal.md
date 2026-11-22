---
paths:
  - "src/dmrgpy/thermal.py"
  - "src/dmrgpy/vevtk/metts*.py"
  - "src/dmrgpy/pyitensor/metts.py"
  - "src/dmrgpy/mpsjulialive/metts.py"
  - "src/dmrgpy/mpsjulialive/metts.jl"
  - "tests/test_audit_*thermal*.py"
  - "tests/test_audit_*metts*.py"
  - "tests/test_metts*.py"
  - "tests/test_thermal*.py"
---

# Finite temperature: purification and METTS

These are durable rules the hole hunts left, with the mechanism and the
numbers; CLAUDE.md keeps one line per rule. Each rule names the test that
pins it. Do not undo one without reading its finding in `docs/audit_*.md`
first. The file loads when `thermal.py` or a METTS file is read or edited.

- **`Thermal_Spin_Chain.anneal()` is order-8 Taylor steps** of
  `exp(-dtau(H-E_mid))` at `dtau*W <= 2`, about `2*beta*W` factor
  applications (`test_audit_2026_09_25b_thermal.py`).
- **`Thermal_Spin_Chain.get_gs()` returns its stored state only under the
  key it was purified under** (`thermal.purification_key`: T, the
  `anneal_*` settings, `mode`, the backend and `groundstate.solver_key`),
  read before any solve. Do NOT go back to a short circuit on `computed_gs`
  alone, or to a setter on T alone (`test_audit_2026_11_21_thermal.py`).
- **A collapse basis METTS cannot sample with is refused** by both calls on
  every backend (`mettsvev._check_basis`): a non-Hermitian name, or names
  all commuting on each site with a charge or parity H conserves, which
  refuses every fermion chain (the author's pick). Do NOT go back to a list
  of diagonal names (`test_audit_2026_10_23*.py`, `test_audit_2026_11_21_metts.py`).
