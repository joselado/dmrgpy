# Plan: make fermionic `Infinite_Many_Body_Chain` work (v3 correctness, `"python"` speed)

**Status: implemented on 2026-08-22** (steps 1-6; Step 5.2's geometric ramp landed on both
backends, and an extra Krylov pass beyond the plan removed the `O(chi^6)` dense fixed
point/environment solves on the `"python"` side, and cross-site fermionic correlators were
threaded afterwards too). The outcome and the before/after numbers are in the commit
messages on `idmrg-fermionic-infinite`; the user-facing result is in
`docs/user_guide.md`'s iDMRG section and `docs/documentation.md`'s backend-dispatch
section. Kept only as a record of the investigation -- the bug report it was written
against has been removed now that it is fixed.

Addresses an external bug report (since removed, now that it is fixed) about
`Infinite_Many_Body_Chain` with native spinful (`site_type=1`) sites being unusable: it
failed instantly on `itensor_version=3` and did not finish on `itensor_version="python"`.
That report had two
independent halves — an immediate hard error on `itensor_version=3`, and impractical
wall-clock on `itensor_version="python"` — plus a test-coverage gap. All three were
re-confirmed in this checkout before writing this plan; the diagnosis below replaces the
report's own guesses about causes.

## What is actually broken (confirmed, not conjectured)

### A. `itensor_version=3`: the infinite path feeds it the *finite-chain* Jordan-Wigner form

`infinitechain.py:477/478` and `:490/491` pass `self._h_intra.to_terms()` /
`self._h_inter.to_terms()` to `Chain::idmrg_ground_state` / `Chain::vumps_ground_state`.
`MultiOperator.to_terms()` (`multioperator.py:158`) applies the **global, finite-chain**
`jordan_wigner()` — every fermionic operator at site *s* gets an explicit `F` factor on
every site `1..s-1`. Dumping the actual terms for the report's model:

```
inter: 16 terms, 16 of them touch >2 distinct sites, e.g.
   (-0.258, [('Adagup',1), ('F',2), ('F',1), ('Aup',3)])
   ( 0.011, [('F',1), ('Adagup',2), ('F',3), ('F',2), ('F',1), ('Aup',4)])
```

`idmrg_classify_terms` (`mpscpp3/chain_session.h:3917`) groups factors by site and errors
when more than 2 *sites* are touched — the explicit string sites are counted, hence the
report's error. Two consequences the report did not identify:

1. The report's guess ("the C++ inspector counts operator *factors* rather than distinct
   site *indices*") is **wrong** — it does group by site correctly, and `_canonicalize_
   hamiltonian` (`infinitechain.py:95`) already enforces the documented ≤2-distinct-sites
   contract on the *raw* terms. The extra sites are created afterwards, by `to_terms()`.
2. The bug is broader than an over-strict check. `idmrg_classify_terms` hardcodes
   `carry_ferm=false` (`chain_session.h:3952`) and `idmrg_build_row` always propagates
   `eye()` instead of `F` on a pending channel (`:4097`) — both files already say so in
   their own comments (`:246`, `:2448`, `:4016`: "fermionic terms not supported yet …
   silently wrong for a fermionic one — no detection/guard"). So today's v3 infinite path:
   - **works by accident** when every fermionic term connects *adjacent* sites (the global
     string then reduces to the correct endpoint-`F` composition with nothing in between),
     which is why the existing spin-only tests never noticed;
   - **hard-errors** as soon as a fermionic term has reach ≥ 2 (the report's case);
   - would be **silently wrong** if the site check were merely relaxed, because the
     pending channel would still carry `Id` where the string needs `F`.

Meanwhile `itensor_version="python"` passes the **raw** `.op` terms and does its own
*local, translation-invariant* Jordan-Wigner inside `pyitensor/idmrg.py::_term_site_
matrices` (endpoint `F`, `carry_ferm` on intervening sites, fermionic reorder sign,
odd-total-parity rejection) — validated against `autompo.HTerm.resolve()` by
`tests/test_infinite_chain.py::test_term_site_matrices_matches_autompo_resolve`. The two
backends are simply running two different JW conventions, and only one of them is right
for an infinite chain.

### B. `itensor_version="python"`: cost is structural, not mysterious

Profiled the report's model (native spinful, `n_uc=2`, `gs_method="vumps"`) at `maxm=8,
maxiter=1`: **12.0 s**, of which **7.2 s is `numpy.core._multiarray_umath.c_einsum`**
across **101 874 calls**. Two multiplicative factors explain the report's "doesn't finish
at `maxiter=2`":

1. **`np.einsum` without `optimize`** in the shared hot helpers — `idmrg_excitations.py::
   _apply_op_ket` (29 799 calls), `_op_transfer_matrix` (6 528), `_cap_right`/`_cap_left`,
   via `vumps.py::_h_ac_action`. `c_einsum` never dispatches to BLAS. Measured against
   `tensordot`/`matmul` equivalents on the exact shapes involved:

   | helper | einsum | tensordot/matmul | speedup |
   |---|---|---|---|
   | `_apply_op_ket` (D=8, d_g=16) | 142 µs | 37 µs | 3.8x |
   | `_cap_right` (D=8) | 75 µs | 14 µs | 5.2x |
   | `_op_transfer_matrix` (D=8) | 563 µs | 98 µs | 5.7x |
   | `_apply_op_ket` (D=30, d_g=16) | 1647 µs | 164 µs | **10.0x** |

   The gap *grows* with D, so it is worst exactly at the report's `maxm=30`.
2. **The D-ramp is linear.** `vumps.py::vumps_ground_state` (`:717`) runs a full
   multi-restart VUMPS solve at *every* `D_cur = 1..D` (3 restarts below D, 4 at D, plus a
   bounded variational-safety-net budget). Confirmed: 32 `_vumps_single_run` calls for
   `maxm=8, maxiter=1`. Total cost is then ~`0.75·D⁴` rather than ~`D³`; at `D=30` that is
   a ~10x multiplier over a geometric ramp.

Native spinful sites make both factors bite harder simply because `d_g = 4·4 = 16`
(vs 4 for spinless `n_uc=2`), not because of anything site-type-specific. **Nothing has
established that the `"python"` backend is *correct* for `site_type=1`** — the report never
got a converged number out of it. Treat it as *slow and unvalidated*, not "slow but working".

### C. No end-to-end fermionic infinite-chain test exists on any backend

`tests/test_infinite_chain.py`'s every `itensor_version=3` test is an `Infinite_Spin_Chain`;
its `_FERMION_SITE_CODE = 0` block is unit-level coverage of `_term_site_matrices` only. No
test drives a fermionic Hamiltonian through `gs_energy()` on either backend. This is why A
survived.

## Fix

### Step 1 — v3: stop JW-ing the terms before handing them to C++

Give `MultiOperator.to_terms()` an opt-out (`to_terms(jordan_wigner=True)`), keeping
`_filter_small` and the 0→1-based shift, and use `to_terms(jordan_wigner=False)` at both
`infinitechain.py` call sites (idmrg and vumps branches). Rationale for this over relaxing
the site check and collapsing pure-`F` string sites: it makes the C++ path a genuine port of
the Python reference (this repo's stated design for `mpscpp3`), keeps exactly one JW
convention alive in the infinite framework, and does not depend on the fragile fact that the
global strings happen to cancel between a term's two endpoints.

Guard the plumbing change for free: for spin models JW is a no-op, so assert in a test that
`to_terms()` and `to_terms(jordan_wigner=False)` are identical for the existing spin
Hamiltonians — the whole existing v3 spin suite then covers the refactor.

### Step 2 — v3: port the local JW into `idmrg_classify_terms`

Mirror `pyitensor/idmrg.py::_term_site_matrices` line-for-line in
`mpscpp3/chain_session.h::idmrg_classify_terms` (both `idmrg_ground_state` at `:2508` and
`vumps_ground_state` at `:2768` call it, so one change fixes both `gs_method`s):

- fermionic reorder sign when sorting factors into site order (`is_fermionic` = name starts
  with `'C'`, matching `ITensor autompo.cc` and `pyitensor/sites/base.py:59`);
- per-site parity `is_site_fermionic`, and the endpoint-`F` factor when the carried parity
  differs (`combined = F @ combined`, exact order — this is the `F@C == -C` sign class the
  Python side already had to fix once);
- set `IdmrgBond::carry_ferm` from the carry between the two touched sites instead of
  `false`;
- reject odd total fermion parity with the same error text as the Python side;
- keep the existing >2-distinct-sites error — with raw terms it now fires only for genuinely
  3-site terms, which is the documented contract.

Same-site factor composition (`M_last @ … @ M_0`) already matches the Python reference and
needs no change. `idmrg_op_dense` resolves raw names via `sites_.op(name, site)`; verified
`HubbardSite`/`ElectronSite` (`get_sites.h:122`, `ITensor/itensor/mps/sites/electron.h:154`)
provides `Cup`/`Cdagup`/`Cdn`/`Cdagdn`/`Nup`/`Ndn`/`F`.

### Step 3 — v3: propagate the string in the MPO builder

`idmrg_build_row` (`chain_session.h:4097`): use the site's own `F` matrix instead of `eye()`
on a pending→pending transition when `bonds[…].carry_ferm` is true. Update the three
"fermionic terms not supported yet" comments (`:246`, `:2448`, `:4016`) to describe what is
now supported.

Rebuild with `make pybind` in `mpscpp3/` (the `$(wildcard *.h)` prerequisite covers
`chain_session.h`; see CLAUDE.md).

### Step 4 — tests (the actual coverage gap)

Add to `tests/test_infinite_chain.py`:

1. **Spinless, `n_uc=2`, reach-2 hopping** — a dimerized (SSH-like) chain plus a small
   next-nearest hopping `t₂·C†₀C₂`. This is the minimal shape that hard-errors today and the
   only one that exercises `carry_ferm=true`. It is free-fermion and gapped, so the exact
   energy density is an analytic Bloch-band integral (`∫dk/2π` over the negative band) — a
   backend-independent anchor, not a golden value. Assert v3 (both `gs_method`s) and
   `"python"` against it.
2. **Native spinful (`site_type=1`), free** — the same Bloch Hamiltonian for both spin
   species; the exact energy density is precisely **2×** the spinless one. This is what
   closes the open question in §B about whether `"python"` is correct for `site_type=1`, and
   it is the *only* check on `ElectronSite`'s on-site spin convention (`electron.h` defines
   `Cdn = Fup·Adn`, i.e. an intra-site string that no spinless site type has) — note that in
   the test so nobody later "simplifies" it away as redundant with test 1.
3. **The report's interacting cell** — v3 vs `"python"` agreement at small `maxm`, plus a
   loose cross-check against a finite `Spinful_Fermionic_Chain_Native` chain's energy per
   site. Reuse the repro verbatim from the known-issue doc.

Keep the parameters small (D ≤ 8–12); tests 1–2 are gapped by construction, avoiding the
gapless slow-convergence trap the report hit.

### Step 5 — `"python"` performance

Two separate commits, in this order (do **not** bundle them):

1. **einsum → `tensordot`/`matmul`** in `idmrg_excitations.py::_apply_op_ket`,
   `_op_transfer_matrix`, `_cap_right`, `_cap_left` (shared by iDMRG, VUMPS *and* the
   excitation ansatz, so every infinite-chain path benefits at once). Purely mechanical,
   numerically identical up to floating-point associativity; the measured 4–10x is the
   single highest-value change. Re-profile afterwards and record the new hot spots.
   While here: `vumps.py::_environments` (`:470`) builds a full `D²×D²` `E4` per bond per
   iteration via `_op_transfer_matrix` only to immediately contract it against `I` —
   contract directly and skip materializing it.
2. **Geometric D-ramp** (`1,2,4,…,D`) in `vumps_ground_state` instead of `1..D`, keeping the
   restart count and the variational safety net unchanged. ~10x at `D=30`, but it changes
   convergence behaviour that `vumps.py`'s own docstring documents real failures around, so
   it needs its own commit and a full run of the existing spin tests (plus the new fermionic
   ones) before it is trusted. Mechanically the warm start already allows arbitrary jumps —
   `_grow_initial_state` (`vumps.py:304`) embeds the old solution in the new tensor's
   top-left block plus noise and assumes nothing about `D - D_old` — but its docstring's
   robustness argument is specifically about ramp *granularity* (a pure-random large-D start
   was observed landing in a basin worse than the same model's D=2 result), and a `D→2D` jump
   leaves proportionally more noise-filled space. So treat the geometric ramp as the primary
   attempt and a hybrid (geometric up to `D/2`, unit steps after) as the named fallback if the
   spin suite's converged energies degrade.

Target: the report's model at `maxm=30` should complete a macro-iteration in seconds, not
tens of minutes.

### Step 6 — examples and docs

- New `examples/idmrg/fermionic_infinite_chain/main.py`: sweep `t₂` (or the dimerization),
  **plot** energy density from v3, `"python"`, and the analytic band integral on one axis
  (CLAUDE.md: examples plot, they do not merely print/assert).
- Retire the originating bug report once it is fixed (its two cause-guesses were both
  wrong; the confirmed diagnosis is in §A/§B above and its reproduction case survives as
  `tests/test_infinite_chain.py::test_interacting_spinful_cell_backends_agree`).
- Update `docs/user_guide.{md,tex}` (fermionic infinite chains are now supported; state the
  ≤2-distinct-sites and even-parity contracts) and `docs/documentation.{md,tex}` (the v3
  infinite path consumes raw, non-JW terms and threads JW itself, unlike every finite v3
  call site). Keep `.md`/`.tex` in sync; re-check `pdflatex`.

## Non-goals

- `n_uc ≤ 2` stays — lifting it is the redesign `infinitechain.py`'s own constructor comment
  describes, unrelated to this issue.
- Fermionic **tangent-space excitations** (`idmrg_excitations.py`, `vumps_excitation_
  energies`) stay untested and out of scope.
- **Fermionic cross-site two-point correlators are out of scope and are currently silently
  wrong on both backends.** Checked directly: `pyitensor/vumps.py::two_point_correlator`,
  `pyitensor/idmrg.py::two_point_correlator` and `Chain::vumps_two_point_correlator`
  (`chain_session.h:2941`) contain no `is_fermionic`/`F`/carry handling at all — they insert
  the operator into the transfer matrix with a plain `'io,lir->lor'` contraction, so a
  `<C†_i C_j>` with `j>i` is missing the Jordan-Wigner string on every intervening site.
  Same failure class as §A, one layer up, but a *separate* code path from
  `idmrg_classify_terms`: Steps 1–3 do not fix it. In scope here: ground-state energy
  density, single-site `vev`, and correlators of parity-even (spin/boson/`N`-type)
  operators. Add a guard so a fermionic opname pair raises instead of returning a wrong
  number, and track the string threading separately.
- `mpscpp2` (ITensor v2) is uninvolved — it has no infinite-chain path at all.

## Ordering / risk

Steps 1–3 are one coherent change (v3 is broken-or-accidental until all three land); Step 4
gates it. Step 5.1 is independent and low-risk — it can land first if a quick win is wanted.
Step 5.2 is the only change with real behavioural risk.
