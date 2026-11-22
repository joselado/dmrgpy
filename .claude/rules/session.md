---
paths:
  - "src/dmrgpy/manybodychain.py"
  - "src/dmrgpy/groundstate.py"
  - "src/dmrgpy/mode.py"
  - "src/dmrgpy/sites.py"
  - "src/dmrgpy/vev.py"
  - "src/dmrgpy/excited.py"
  - "src/dmrgpy/entropy.py"
  - "src/dmrgpy/entropytk/**"
  - "src/dmrgpy/mpsjulialive/groundstate.py"
  - "src/dmrgpy/mpsjulialive/get_gs.jl"
  - "tests/test_audit_*session*.py"
  - "tests/test_audit_*groundstate*.py"
  - "tests/test_audit_*construction*.py"
  - "tests/test_audit_*fluctuation*.py"
  - "tests/test_audit_2026_09_25b_julia.py"
  - "tests/test_audit_2026_10_23*.py"
  - "tests/test_audit_2026_11_21_corrmat.py"
  - "tests/test_restart_drops_session_state.py"
---

# Ground state, session, dispatch and static observables

These are durable rules the hole hunts left, with the mechanism and the
numbers; CLAUDE.md keeps one line per rule. Each rule names the test that
pins it. Do not undo one without reading its finding in `docs/audit_*.md`
first. The file loads when `manybodychain.py`, `groundstate.py`, `mode.py`, `sites.py`, `vev.py` or a file
under `entropytk/` is read or edited.

- **`restart()` is a cold start on every backend.** It calls the session's
  `clear_wavefunction()` (v2, v3, `"python"`), so `set_hamiltonian(H)`
  solves from a fresh random MPS, and only `restart=False` or an explicit
  start keeps a warm one. A v3 warm start cannot leave the total-Sz sector
  of its state when H conserves Sz, and a field sweep on one chain stayed
  in an old sector (`test_restart_drops_session_state.py`).
- **A state becomes the chain's normalized.** `groundstate.unit_copy` is
  applied in `mark_injected`, in the `wf0=` copies and in `set_gs`'s ED
  branch, and a zero state raises; an explicit `wf=` to
  `evolve_and_measure`/`evolution_ABA` is still evolved as given
  (`test_audit_2026_09_25b_session.py`).
- **Every correlator starts from `ground_state_on_session`**, with
  `groundstate.mark_injected` marking a state the caller supplied, so that
  a DMRG correlator after `set_gs()`, `set_initial_wf` or
  `set_initial_wf_guess` measures that state and not the session's own
  solved one (on `julia_live` too, since 2026-09-25); the mark's third field,
  `supplied`, read through `groundstate.state_supplied()`, is what SECTOR
  uses to refuse a set state that is not its reference sector's ground
  state (`test_audit_2026_09_24b_groundstate.py`,
  `test_audit_2026_09_24c_session.py`).
- **`groundstate.stored_answer_holds` is the one condition** under which
  `gs_energy`/`get_gs` return a stored answer; a keyword the stored answer
  would ignore (`maxde=`, `reconverge=True`, `wf0=`) raises or re-solves.
  After `gs_energy_generalized`, `e0` is the state's own energy `<wg|H|wg>`
  (biorthogonal on the NH route) and lambda is kept as `lam_generalized`
  (`test_audit_2026_09_25b_session.py`).
- **The session's lower band edge, `minimum_energy()`, never sweeps a state
  the session holds without its energy**: it solves from a fresh start and
  puts the state back (`test_audit_2026_09_24c_session.py`).
- **`vev.energy_variance` is how an energy variance is computed**,
  `||(H-<H>)|psi>||^2` with the `maxm` cap lifted. `<H^2>-<H>^2` through
  `vev(npow=2)` is a small difference of two numbers of order E0^2 whose
  error is set by the truncation of H|psi>, and H^2 as an MPO is O(L^3) on
  `"python"` (`test_audit_2026_09_24c_fluctuation.py`).
- **A chain constructor keyword is a setting assigned after `initialize()`**,
  exactly as an assignment afterwards, `mode` included (`sites.check_settings`,
  the explicit allowlist `sites.SETTINGS`); anything else, a Hamiltonian
  accumulator such as `hubbard=` included, raises
  (`test_audit_2026_09_25_construction.py`,
  `test_audit_2026_09_25b_construction.py`).
- **`mode.resolve_mode` lets only a chain pinned to `"ED"` override a call's
  `mode=`**; a `"DMRG"` pin is the same as None, and both names are
  validated before any fallback can return (`test_audit_2026_09_25b_session.py`).
- **`julia_live` starts every DMRG solve from a random MPS** of link
  dimension min(maxm,10), with the chain's noise on the first half of the
  sweeps, so its solves vary run to run by solver noise and a warm 16-site
  solve is about 40% slower than before; the KPM window there trusts `e0`
  as the lower band edge only when `groundstate.mark_lower_edge` says a
  plain solve produced it (`test_audit_2026_09_25b_julia.py`).
- **`gs_energy_generalized` refuses a metric it cannot use**: `ValueError`
  for a non-Hermitian A, `RuntimeError` once `<psi|A|psi>` is not positive,
  on v3 through a C++ throw, never ITensor's aborting `Error()`; an NH-DMRG
  state that did not converge warns again on every reuse
  (`groundstate.warn_unconverged`). Entropies of a site or a pair are those
  of the normalized state, with a relative eigenvalue floor
  (`test_audit_2026_10_23_leads.py`).
- **The DMRG `vev` has no `**kwargs`**: an unknown keyword raises, and
  `T>0` raises naming `mode="ED"` and `metts_vev()`
  (`test_audit_2026_10_23.py`).
- **One predicate, `correlationentropy._full_applies`, picks the default
  `dmmode`** in a sector and outside: `"full"` only for the chain's own C
  on spinless sites with a one-sweep session, else `"explicit"` in a sector
  and `"fast"` outside. Do NOT add a condition to one branch of
  `_default_dmmode` only (`test_audit_2026_11_21_corrmat.py`).
