# Numerical changes ledger

This file lists what raises an error and which results changed in each audit.
It was section 21 of `docs/user_guide.md`. Read it before comparing a number
with a result computed before 2026-09-26. Section numbers mentioned below
refer to the user guide as it was then. The full records are in
`docs/audit_*.md`.

## The 2026-08 audit

A cross-backend audit in August 2026 (`docs/audit_2026_08_hole_hunt.md`,
which records every reproduction) went looking for calls that silently did
something other than what was asked. Most of what it found is now either
fixed or loud. Two of the fixes **change numbers** and are worth knowing
about if you have results from before them:

- **`submode="TD"` and `submode="TDZ"` were a factor of $\pi$ too large.**
  The Fourier transform behind them omitted the $1/\pi$ of the convention
  $S_{AB}(\omega)=\frac1\pi\mathrm{Re}\int_0^T\!dt\,e^{i\omega t}C(t)w_\delta(t)$
  that this guide documents and that every other submode already followed,
  and gave the $t=0$ sample full rather than half weight. On a 4-site chain
  the exact sum rule $\int S(\omega)d\omega=\langle AB\rangle=0.25$ came
  out as 1.28. Peak positions and widths were never affected -- the error
  is uniform in $\omega$ -- so only absolute spectral weight changes.
  `sxt_to_skomega` and the infinite-chain $S(k,\omega)$ reduction share the
  same transform and inherit the fix.
- **`itensor_version=3` real-time evolution from a non-unit-norm state.**
  `evolution_ABA` (whose start state is $A|\mathrm{gs}\rangle$) and
  `evolve_and_measure(wf=...)` were rescaled by $1/\|\psi\|^2$ for every
  $t>0$ under `tevol_method="TDVP"`/`"TDVP_GSE"`, while $C(0)$ stayed
  correct. Against ED on a 5-site chain the error went from 0.176 to
  1.7e-9. `"TEBD"`, `"MPO"`, `itensor_version=2` and `"python"` were never
  affected.

Beyond those, the following now raise where they used to be silent:

| Call | Was | Now |
|---|---|---|
| `sc.maxm = m` with `m<1` | SIGABRT (v2/v3), wrong energy (`"python"`) | `ValueError` |
| changing `maxm`/`nsweeps` between `gs_energy()` calls | returned the first, less converged energy | re-solves |
| `sc.tevol_method = "tdvp"` (any unrecognized name) | ran the legacy MPO-Taylor integrator | `ValueError` |
| `submode=` on a non-Hermitian Hamiltonian | every submode returned the CVM resolvent | dispatches; `NotImplementedError` for the Hermitian-only ones |
| unknown keyword to `submode="KPM"` | silently discarded | `TypeError` |
| `name="ZZ"` with `submode="KPM"`/`"EX"`, or any `mode="ED"` | `RuntimeError: No active exception to reraise` | works |
| `submode="CVMimag"`/`"maxent"` | `exit()` -- terminated your process | `NotImplementedError` |
| `kpm_energy_truncate=True` on `itensor_version=2` | silently ignored | `NotImplementedError` |
| `kpm_energy_truncate=True` with far too small a `kpm_scale` | plausible spectrum at the wrong energy | `RuntimeError` in the worst regime (see `docs/known_issue_kpm_energy_truncation_window.md`) |
| an operator on an out-of-range site (`"python"`) | treated as the identity | `IndexError` |
| `get_excited_states(n=...)` beyond the Hilbert space | ground on for minutes | `ValueError` |
| `Thermal_Spin_Chain(..., T<0)` | treated as `T=0` | `ValueError` |
| `Thermal_Spin_Chain(..., itensor_version=...)` | dropped | honored |
| `get_distribution*(mode="ED")` | `AttributeError` deep inside | dispatches, or `NotImplementedError` |
| `vev(C[i]*C[i])` | `AttributeError` on every DMRG backend | `0` (as ED always answered) |
| `Many_Body_Chain.evolution()` | `AttributeError` (it called a function that does not exist) | removed |
| `set_hoppings()`, `set_hubbard()`, their `_spinful`/`_spinless` variants, `set_hoppings_MB()`, `set_hubbard_MB()`, `set_pairings_MB()`, `update_hamiltonian()` | built the Hamiltonian as $\sum_{ij}f(i,j)A_iB_j$ over both orderings, keeping every zero term | removed; write the sum and call `set_hamiltonian()` |
| `Fermionic_Chain.hamiltonian_free()`, `gs_energy_free()`, `get_correlator_free()`, `get_gr_free()` | `AttributeError` (they read `self.hoppings`, which nothing sets) | removed |
| `gs_energy_generalized()` on a chain below 2 sites (`"python"`) | the Rayleigh quotient of an untouched random state (-0.3049 for an exact -0.5) | `RuntimeError` |

And these combinations now work where they used to fail:

- **conserved-sector mode + `get_rdm()`** on `itensor_version=3`, which
  returned uninitialized heap memory (non-Hermitian, different every run);
- **conserved-sector mode + a repeated-site operator product** such as the
  four-point correlator's `Cdag_0 C_1 Cdag_0 C_1`, which aborted the whole
  process inside ITensor and now evaluates to 0 as it should;
- **conserved-sector mode + site/pair entropy, mutual information, and the
  default `get_correlation_matrix()`**, all of which raised;
- **`get_rdm(i=ns-1)`** (the last site) on `itensor_version="python"`, and
  on `"julia_live"` since the 2026-09-25b pass, where it raised
  `BoundsError`.

One thing the audit's own `name=` fix broke, since fixed: centralizing
`name=` resolution made every submode go through the same normalizer,
whose explicit-pair branch admitted only `MultiOperator`s -- so
`name=(A,B)` with the already-built operators `toMPO()` returns started
raising `ValueError` on the paths that had always accepted them
(`mode="ED"`, and `submode="EX"`). They are accepted again, per element
and including products of them; see *What `name=` accepts* in §6 for
which submodes can and cannot take one.

**The ED Lehmann sum is now exact (2026-09-02).** `mode="ED",
submode="ED"` -- the reference every other submode is validated against --
truncated its final-state sum to `emu < max(es)`, comparing *absolute*
eigenvalues against a frequency measured from the ground state. Three
consequences, all reproduced: adding a constant to `H` changed the answer
(1.1e-3 at `H+3*Id` on a 6-site Heisenberg chain) although a constant can
move no pole; a spectrum sitting entirely above `max(es)` truncated itself
away and raised numba's `"zero-size array to reduction"` from deep inside
the sum; and the result depended on how wide an `es` window the caller
happened to ask for. `dynamical_correlator_finite_T` (the `T>0` route,
and the reference for the METTS correlator) carried the same line.

Both now sum over the full spectrum. **For almost every existing caller
this changes nothing measurable**: whenever the ground-state energy is
comfortably negative -- the ordinary case, and why this survived since
2025-11 -- the old threshold retained everything that mattered, and the
difference is float noise (1.4e-17 on the 5-orbital IETS atom of
`tests/test_atom_iets.py`, 8.3e-13 on a 6-site Heisenberg chain). It
differs only where the old code was wrong. Nor is it slower: the `T=0`
sum reads only `nex` rows of one matrix-element array and `nex` columns
of the other, so only those are now built, which is cheaper than the full
`n x n` products the truncated version was computing. The `T>0` route
*does* cost more, and unavoidably so: it sums over every initial state
with its own Boltzmann weight, so every row and column is genuinely read
and no such shortcut exists. That is O(dim^2) in the sum, the same order
as the dense diagonalization it already performs -- 0.47 s on a
1024-dimensional Hilbert space, against 0.004 s for a 6-site chain.

### The 2026-09 audit

A second cross-backend audit (`docs/audit_2026_09_hole_hunt.md`) recorded
36 confirmed findings, all of them now fixed. The last to land was the
iDMRG matrix-free item, which stood as PARTIAL for a while because its
compute half arrived without its memory half; the second pass closed on
2026-09-22, and what it buys an infinite-chain run is in §18. It also
left two *open items*, found while
re-measuring the record rather than during the hunt itself, and both of
those were closed on 2026-09-22: the real-time submodes sat off the
correlator convention (O1) and `delta` meant a different width under
`submode="KPM"` on each solver (O2). The list below is the
physics-facing half — what changes for a user who never reads the
architecture docs.

**Results that are not comparable across this change.** Read this the way
you read the two 2026-08 items above.

- **`submode="CVM_explicit"` returned exactly twice the correct
  spectrum**, at every frequency. The submode is DMRG-only (`mode="ED"`
  raises `NotImplementedError`), and the factor of 2 was the same on
  every DMRG backend, the fix being one backend-agnostic line.
  Verified directly on the 6-site staggered Heisenberg chain:
  the first four points came back as
  `[0.2169497, 0.0920095, 0.1731231, 0.1051759]` against `submode="CVM"`'s
  `[0.1084749, 0.0460046, 0.086561, 0.0525879]`, i.e. a ratio of 2 to six
  digits; they now agree to 5.7e-07. Any peak height, integrated weight
  or figure produced with `"CVM_explicit"` before this is a factor of 2
  too large. This is the one correlator change that hits Hermitian pairs
  as well.
- **Three dynamical-correlator routes moved onto the house convention**
  of §6 (the complex Lehmann density) from $-\frac1\pi\mathrm{Im}\,G^R$:
  `mode="DMRG" submode="CVM"`, `submode="ROOTN"` under both `mode=`
  values, and `mode="ED" submode="ED"` (at $T=0$ and $T>0$), which is
  the exact reference every other submode is validated against. **Where
  every $M_n$ is real nothing moves at all** — the two conventions are
  then identical (measured to 1.1e-16), and that covers every
  `name=(A,B)` example in this guide, Hermitian pairs and the
  $S(q,\omega)$ sweep's $(S^z_i,S^z_j)$ alike; §6 states the criterion
  and why $A=B^\dagger$ is not it. For a pair with complex $M_n$ the curves move by more than
  their own peak: on a 4-site complex-hopping chain with
  $A=c^\dagger_0$, $B=c_2$ at $\delta=0.15$ the two conventions differ by
  up to 5.720e-01 against a peak of 0.61. `submode="ED"` also returns a complex array
  now where it returned `float64` before (same real part, to the digit;
  the imaginary part it used to drop peaks at 0.5682 on that pair), so a
  caller that assumed a float array should take `.real`. `submode="KPM"`
  (the default), `"INV"`, `"CVM"` under `mode="ED"`, `"EX"`, `"TD"`,
  `"TDZ"` and `"SECTOR"` were untouched *by that change*; `"TD"`,
  `"TDZ"` and `"KPM"` moved later, on 2026-09-22, for the two reasons in
  the following two entries.
  `mode="ED" submode="ROOTN"` and `mode="DMRG" submode="ROOTN"` are also
  roughly **2x slower** now, for the reason given in §6.
- **`submode="TD"` and `submode="TDZ"` now return the density too**,
  where they used to return the complex one-sided transform it is
  assembled from (the audit's open item O1, closed 2026-09-22). §6's
  `"TD"` entry carries the identity and the cost rule; what matters for
  a result you already have is that the array moved, and that the
  imaginary part moved most. On the audit's 4-site complex-hopping chain
  ($A=c^\dagger_0$, $B=c_2$, $\delta=0.4$, `dt=0.1`, exact peak 0.2313)
  the distance to the exact Lehmann density went from 1.16e-01 to
  2.57e-04 under `mode="DMRG"`, from 2.76e-01 to 2.57e-04 under
  `mode="ED"`, and from 1.13e-01 to 3.31e-02 for `"TDZ"`, a residual
  this entry used to put down to the contour and which was in fact the
  frequency stage's FFT grid; since the 2026-09-24 audit (finding 7, in
  the subsection below) `"TDZ"` lands 5.4e-04 from exact and `"TD"`
  3.2e-05. On a Hermitian pair the
  imaginary part used to run at 60% of the peak (1.11e-01 against
  0.1869 on a 4-site Heisenberg chain at $\delta=0.3$) and is now
  exactly zero, with the curve itself landing 2.82e-04 from exact. Taking
  `.real` of the result, which this guide used to instruct, is now a
  no-op on a real-weight pair and *wrong* on a complex-weight one.
  `mode="ED"` moved further than `mode="DMRG"` because that route
  additionally read the operator pair in the opposite order, returning
  $S_{BA}$ where every other route returns $S_{AB}$; it was invisible
  whenever $A$ and $B$ are the same operator, which is every example
  here. A pair that is not provably its own adjoint now costs two
  evolutions rather than one.
- **`delta` under `submode="KPM"` is now the same broadening it is in
  the resolvent submodes**, FWHM $=2\delta$ at the band centre, on both
  solvers (the audit's open item O2, closed 2026-09-22), and
  `kpm_n_scale`'s default changes from 3 to 1 because the base moment
  count it multiplies is now the calibrated one. Neither solver realized
  the intended width before, and they missed it differently: `mode="ED"`
  rescaled onto a window three times the bandwidth and counted moments
  by its own rule, giving about $0.93\delta$, while `mode="DMRG"` gave
  about $1.6\delta$. Every KPM curve therefore changes width, and the
  two solvers stop disagreeing: measured with $A=B=S^z_0$ on three
  chains, $\max|\mathrm{ED}-\mathrm{DMRG}|$ against the resolvent peak
  went from 95 to 124 per cent down to 4 to 7 per cent, and to 0.07 to
  0.2 per cent on `itensor_version="python"` once the 2026-09-24 audit
  cut the DMRG routes to the calibrated moment count (finding 4, in the
  subsection below). Peak heights
  change with the width and the sum rule does not, which is exactly why
  nothing caught this: both routes integrated to 0.250000 against an
  exact 0.250000 before and after. A script that was tuning `delta` by
  eye against a KPM curve, or carrying `kpm_n_scale=3` explicitly, gets
  a different spectrum now. `get_distribution()`'s own KPM path is
  deliberately left off the calibration (§14).
- **`vev(op, npow=n)` on an ED route** returned $\langle O\rangle$ for
  every $n$; it now returns $\langle O^n\rangle$. Measured on a 4-site
  Heisenberg chain, `vev(Sz[0], npow=2, mode="ED")` went from 5.6e-17
  (i.e. $\langle S^z_0\rangle$) to 0.25.
- **`gs_energy_fluctuation` on an ED route** inherited that and also
  ignored `mode=`. On the same chain `sc.mode="ED"` reported 2.0561 and
  now reports 0.0; a 2-site chain on the *default* backend (which
  `mode.py` routes to ED by itself) reported 1.1456 and now reports
  1.05e-08. Anything using it as a convergence criterion behaves
  differently on those routes. DMRG-route values are bit-unchanged.
- **`exponential(h, wf)` on the DMRG backends** — see §2. Measured as
  $\langle\mathrm{GS}|e^{zH}|\mathrm{GS}\rangle$ on a 4-site Heisenberg
  chain under `itensor_version="python"`, against the ED value: 0.6776
  vs 0.6676 at $z=0.25$, 0.5184 vs 0.4457 at $z=0.5$, 0.6897 vs 0.1987
  at $z=1$ (1.5%, 16%, 247%). It now agrees with ED to ~1e-8 at all
  three.
- **Boson occupation projectors under `mode="ED"`** were wrong on any
  Hamiltonian that does not conserve the total boson number. On a
  3-site, 4-level chain with an $a_i+a^\dagger_i$ drive, site 0's
  $P(n)$ read `[0.0042, 0.0054, -0.048, 0.0706]`, summing to 0.032 —
  a negative probability included. It now reads
  `[0.0322, 0.2393, 0.4732, 0.2553]`, summing to 1.0000 and matching
  DMRG. A number-*conserving* Hamiltonian is unaffected to machine
  precision, which is why nothing caught this; $\langle n_i\rangle$,
  the energies and every DMRG backend were always right.
- **`SpinBoson_Chain(..., maxnb=[...])`** silently built the 4-level
  site regardless: `SpinBoson_Chain(["B","S=1/2"], maxnb=[6,None])` gave
  sites `[104, 2]` and now gives `[106, 2]`. Earlier results at a
  non-default `maxnb` are results at `maxnb=4`. With no `maxnb=` nothing
  changes.
- **`itensor_version="python"`**, several ways. Its `applyMPO` used to
  truncate in a gauge where discarding the smallest singular values is
  not the optimal truncation, so every repeated-application consumer
  moved at the 1e-5..1e-9 level: `applyoperator`, `vev(..., npow>1)`,
  `gs_energy_fluctuation` (5.06e-05 → 6.27e-06 on a 10-site Heisenberg
  chain at `maxm=30` — the old number was mostly the MPO-application
  error rather than the state's), the KPM Chebyshev recursion and
  therefore every KPM correlator and `get_distribution`, CVM,
  `applyinverse`, and the MPO-Taylor evolution paths. The new values are
  the ones that match `itensor_version=3` and ED. The product is now
  formed exactly and truncated once afterwards, which is more work per
  application; the cost is workload-dependent, so budget for it rather
  than expecting a fixed factor. Separately, its TDVP real-time
  evolution did not canonicalize the state before the first half-sweep,
  so a trajectory whose start state was not already in that gauge took
  its first step under an unprojected generator; the audit record
  (`docs/audit_2026_09_hole_hunt.md`, finding #1) carries the
  reproduction and the mechanism. Which start states escaped it is not
  recorded there and is not claimed here. And a chain object reused across
  several `set_hamiltonian()` calls could keep a start state belonging to
  the previous Hamiltonian. Both are fixed; no
  `itensor_version` other than `"python"` was affected by any of this.
- **`promote_to_dense()` now keeps the sector's energy on
  `itensor_version="python"`** — the guarantee §3 states ("a bare
  `gs_energy()` afterwards still returns the sector's energy rather than
  re-solving") was true on `itensor_version=3` and `mode="ED"` but not
  there: promotion invalidated the session's cached energy, and
  pyitensor's unconstrained re-solve then left the sector. On a 3-site
  Hubbard chain confined to `Nf=3` the audit recorded a post-promotion
  energy of $-2.3399130755$ with $\langle N\rangle=2$ against the
  sector's own $-1.9408140222$ with $\langle N\rangle=3$; it now returns
  the sector value, re-measured here. Any `"python"` result that called
  `gs_energy()`/`get_gs()` after `promote_to_dense()` — the documented
  workflow — is a different number than it was. The other backends are
  unchanged, and a chain where the sector and global ground states
  coincide never saw it.

**And these now raise, or now work, where they used to be silent:**

| Call | Was | Now |
|---|---|---|
| `Bosonic_Chain(n, maxnb=[...])` with `maxnb != 4` on `itensor_version=2` | SIGABRT inside ITensor | `ValueError` naming the working backends |
| `Bosonic_Chain`/`SpinBoson_Chain`/`Parafermionic_Chain`/`Spin_Fermion_Hamiltonian(..., itensor_version=...)` | `TypeError: unexpected keyword argument` | honored |
| `SpinBoson_Chain(["boson", ...])`, or any unknown site label | `RuntimeError: No active exception to reraise` | `ValueError` listing the accepted labels |
| `SpinBoson_Chain(..., maxnb=[...])` of the wrong length, or `n=` contradicting `len(sitesin)` | silently dropped | `ValueError` |
| `SpinBoson_Chain` under `mode="ED"` | failure inside an unimplemented stub | `NotImplementedError` naming the DMRG alternatives |
| `vev(..., npow=n)` with `n<0`, or with `n!=1` at `T>0`, **on an ED route** | accepted silently, answering with $\langle O\rangle$ | `ValueError` / `NotImplementedError` (the DMRG routes still answer) |
| `gs_energy_fluctuation(npow=...)` | accepted and misused | `TypeError` |
| `exponential(h, wf)` with `h` neither Hermitian nor anti-Hermitian, on DMRG | a warning and a number | `NotImplementedError` |
| `applyoperator`/`summps`/`applyinverse`/`scale_mps`/`exponential` with a `mode=` naming the other backend | ran the wavefunction's own backend anyway (or a bare `raise`) | `TypeError` saying which backend the state belongs to |
| `get_rdm()` on any ED route | `AttributeError` several frames deep | `NotImplementedError` (on `"julia_live"` since the 2026-09-25b pass) |
| the bond entanglement entropy on any ED route | `AttributeError` several frames deep | `NotImplementedError` pointing at `get_site_entropy`/`get_pair_entropy` (on `"julia_live"` since the 2026-09-25b pass) |
| `get_bond_entropy(wf, i, j)` with a site out of range | the audit recorded a SIGABRT on the C++ backends | `IndexError` |
| `get_hamiltonian()` before `set_hamiltonian()` | `'int' object is not iterable`, or `takes 0 positional arguments` | `ValueError` naming the fix |
| a misspelled `ctmode`/`dmmode`/`fpmode` | `RuntimeError: No active exception to reraise` | `ValueError` listing the accepted values |
| a misspelled `basis=` in `get_correlation_matrix` | the electron-basis matrix, silently | `ValueError` |
| `submode="CVM_explicit"` with $A^\dagger\neq B$ | a `print` then `RuntimeError: No active exception to reraise` | `NotImplementedError` naming the submode; since the 2026-09-25b pass it returns $C[A,B]$ |
| `get_distribution()` with no `X=` | `RuntimeError: No active exception to reraise` | `TypeError` naming `X=` |
| `MultiOperator / "a"`, `*` by an ndarray, or `multioperator.obj2MO` on an object it does not recognize | `RuntimeError: No active exception to reraise` | `TypeError` naming what was passed |
| `sc.mode = "dmrg"` (any unrecognized spelling) | accepted, and silently decided by the automatic fallbacks | `ValueError` |
| a conserved sector under `mode="ED"` + the `get_correlation_*` family | `AttributeError` | works, via the backend-agnostic `dmmode="explicit"` |
| `get_dynamical_correlator()` on `Parafermionic_Chain` without an explicit `mode=` | the audit recorded a process abort | works, and obeys `self.mode` and the automatic ED fallbacks like every other chain class |

`Parafermionic_Chain` lost its own `get_dynamical_correlator` override in
the process — it was the base-class method minus the mode resolution and
minus the `name="..."` string resolution — so both the documented string
form of `name=` and `mode=`/`self.mode` work there now.
`Parafermionic_Chain.test(ntries=...)` is honored too, as
`Spin_Chain.test(ntries=...)` already was.

### The 2026-09-24 audit

A third hunt (`docs/audit_2026_09_24_hole_hunt.md`, five lenses over the
commits `1c2606d..765b537` that closed the 2026-09 audit's open items)
recorded 16 confirmed findings. Every code defect among them is fixed,
with one residual kept by design, `submode="TDZ"` under `set_pad_bonds`
at `tdvp_gse_sweeps=0` still being started from the padded state (§19),
and finding 6, which was documentation only, is closed by §6's
`"KPM"` entry.
Several of the findings corrected what the entries above say, the
`"TDZ"` residual and the KPM agreement between the two solvers among
them, and those entries now point here. Each item below names its
finding, and the record carries the reproduction that was actually run.

**Results that are not comparable across this change.**

- **`submode="TDZ"`, and `submode="TD"` at `predict=False`, were
  evaluated on an FFT grid** of spacing $2\pi/(n_t\,dt)$, about
  $1.05\delta$ at the default `damping_periods=6`, and interpolated
  linearly onto `es`, meaning that a Lorentzian of half-width $\delta$
  got about one sample per $\delta$. The damped sum is now evaluated at
  each requested frequency. On the 2026-09 audit's complex-hopping chain
  at $\delta=0.4$ (exact peak 0.2313) the distance to the exact density
  goes from 3.314e-02 to 5.408e-04 for `"TDZ"` and for `"TD"` at
  `predict=False`, and from 2.571e-04 to 3.169e-05 for `"TD"` at its
  default, and the returned curves move by up to 14 per cent of the peak
  for the first two (19 per cent at $\delta=0.25$) and by 2.89e-04 for
  the third. `sxt_to_skomega`, the infinite-chain $S(k,\omega)$
  reduction, stayed on the FFT stage and returned the same bits as
  before (§6, finding 7), until the second pass below.
- **Every DMRG KPM spectrum moves by roughly $2/N$ onto the ED one**,
  $N$ being the calibrated moment count, since the DMRG routes
  reconstructed from $N+2$ moments. On two decoupled Heisenberg dimers
  the band-centre peak goes from 0.666761 to 0.620591 at $\delta=0.2$ and
  from 1.266336 to 1.220236 at $\delta=0.1$, on `itensor_version=3` and
  `"python"` alike, both now equal to the ED peak, and the three
  residuals of the 2026-09 record's O2 table go from about 1.3e-02 to
  2.55e-04, 2.09e-04 and 3.50e-04 on `"python"`. `mode="ED"` does not
  move (§6, finding 4).
- **`get_distribution(mode="ED")` moves by exactly the factor `scale`**:
  on a 4-site open Heisenberg chain with $X=S^z_0$ at the default
  `scale=10` and $\delta=0.05$ the total weight goes from 0.1000 to
  1.0000 and the peak from 0.517316 to 5.173159, and with `xs=` the
  imaginary part, which was a copy of the real one, goes from
  $\max|\mathrm{Im}\,y|=4.96$ to exactly 0. Nothing moves on any DMRG
  backend (§14, finding 3).
- **`ISy` in the first operator of a correlator** returned exactly minus
  the correlator under `"KPM"`, `"EX"`, `"TD"` and `"TDZ"`, since
  `get_dagger()` treated this anti-Hermitian name as Hermitian. On a
  4-site Heisenberg chain in a field, the KPM curve of
  `(get_operator("ISy",0), Sz[1])` goes from exactly minus the curve of
  `(1j*Sy[0], Sz[1])` to identical to it, and its first moment from
  +0.079470 to -0.079470 (exact -0.079685); `sc.is_hermitian(ISy)` is
  now `False` and `sc.is_hermitian(1j*ISy)` `True`, where both were the
  other way round. Only operators that contain `ISy` move (§6,
  finding 2).
- **An operator that mixes pre-Jordan-Wigner `C`-type names with bare
  `A`-type ones** was canonicalized with the wrong sign, so `simplify()`
  changed its value and the Hermiticity, zero and adjoint-pair proofs
  could be false. Such terms are now left as written (§2). A `"TD"`
  pair of that kind on a 4-site complex-hopping chain goes from
  3.616e-02 off exact, with one evolution, to 3.063e-05 with two, on
  `"python"` and v3 alike, and no dmrgpy chain class builds such a term,
  so no ordinary Hamiltonian moves (finding 1).
- **`mode="ED"` `submode="TD"` on a degenerate ground state** built the
  two halves of the density on two different randomly chosen ground
  states: three identical calls on a degenerate 3-site Heisenberg chain
  spread by 2.03e-01 against a 1.768e-01 peak. Both halves now use the
  cached ground state, so the calls are identical and agree with
  `"KPM"`, `"INV"`, `"CVM"` and `"ROOTN"`, though still not with
  `submode="ED"`'s average over the manifold, a split between the two ED
  conventions that predates this fix. `evolution_ABA(mode="ED")` and
  `evolve_and_measure(mode="ED")` without `wf=` change the same way, and
  a non-degenerate ground state moves by at most 3e-11 (finding 8).
- **The Kondo potential-interference term on a non-uniform `es`**
  weighted every point with the grid's first spacing. On a single
  $S=1/2$ at 10 T the peak on a sinh-mapped grid goes from 68.2332 to
  0.6704 against an exact 0.6709, and a uniform grid moves by 2.4e-08
  (§17, finding 13).
- **`excitation_energies(k, n)` with `n>=2` on `itensor_version="python"`**,
  and `spectral_weights`/`dynamical_structure_factor` with it, could
  drop one copy of a degenerate level of $H_{\rm eff}(k)$ above the dense
  threshold (dimension 256) and return the next distinct level in its
  place. At $k=0.37$ on the `n_uc=2` critical Heisenberg cell at
  `maxm=10`, `[0.289896072 0.291597913 0.313428946]` becomes
  `[0.289896072 0.289896072 0.291597913]`, and the $S^x$ and $S^z$ weight
  fractions go from 0.058 to 0.184014 and from 0.339 to 0.041282, both
  now equal to the dense path's. `n=1`, `excitation_gap` and dimensions
  up to 256 are byte-identical (§18, finding 15).
- **`tevol_method="TDVP_GSE"` under `set_pad_bonds(K)`** on
  `itensor_version="python"` ran a different algorithm, one-site TDVP on
  the manifold of bond dimension $K$, its bond growth coming from a QR
  completion of the padded zero directions rather than from the Krylov
  expansion. Padded runs through `quench_tdvp_gse` and
  `evolve_and_measure_tdvp_gse` now follow the unpadded ones: on a
  10-site XXZ quench at $K=$`maxm`$=4$ the two went from 0.4928 apart to
  8.9e-16 at `tdvp_gse_sweeps=0` and from 4.1e-06 to 1.7e-07 at
  `tdvp_gse_sweeps=3`. Note the direction at `tdvp_gse_sweeps=0`, where
  the padded run used to sit 9.8e-05 from ED only because the larger
  manifold held the answer, and now sits 0.4929 from it, like the
  unpadded frozen product state, which is the correct one-site answer
  (§19, finding 16).

**And these now raise, or now work, where they used to be silent:**

| Call | Was | Now |
|---|---|---|
| `kpm_n_scale` not a positive integer, on `"python"`, `mode="ED"` or `"julia_live"` | rounded down silently (1.5 gave 1x, anything below 1 the 16-moment floor) | `TypeError`/`ValueError` naming it, as v2/v3 already did; `kpm_n_scale=2.0` raises too |
| a misspelled keyword to `submode="TDZ"` (`alpha=`, `nmax=`) | ignored, the default-parameter spectrum bit for bit | `TypeError` |
| a `toMPO()` operator to `submode="TDZ"` | a crash several frames deep | `TypeError` naming `toMPO` |
| `get_dynamical_correlator_MB(name="ZZ", i=1, j=1)` under `"TD"`/`"TDZ"` | $C[S^z_0,S^z_0]$, the sites dropped | the sites asked for (the public route was never affected; `"ROOTN"` on this route followed on 2026-09-25) |
| `get_kondo_spectrum(mode="ED")` with a keyword it does not read | ignored, so `Jrho=` for `Jrho_s=` removed the Kondo peak | `TypeError` naming the unknown keywords |
| `get_kondo_spectrum(mode="DMRG")` at a degenerate ground state | one arbitrary member's spectrum | the same by default; `n_gs=g` gives the equal-weight average over the manifold that `mode="ED"` takes |
| the potential term with an `es` that misses spectral weight | a silently low spectrum | `RuntimeWarning` from the sum rule |

### The 2026-09-24 second pass

A fourth hunt (`docs/audit_2026_09_24b_hole_hunt.md`, five lenses over
the single commit `30200a4` that closed the third) recorded 18 confirmed
findings, twelve of them older than that commit and reached by probing
next to it, and every one is fixed. Each item below names its finding;
the record carries the reproduction that was actually run.

**Results that are not comparable across this change.**

- **Every infinite-chain $S(k,\omega)$**, on v3 and `"python"`, was
  mirrored in frequency and interpolated off an FFT grid. It is now
  conjugated after the spatial sum and evaluated directly at each
  frequency: on the transverse-field paramagnet the $k=\pi$ magnon moves
  from $\omega=-1.045$, with 0.928 of its weight below zero, to $+0.905$
  against $\varepsilon=0.900$, and on a single magnon at $\delta T=6$ the
  curve goes from 21 per cent of the peak off the closed form to 2e-14
  (§18, findings 7 and 8).
- **`evolve_and_measure` and `evolution_ABA` on `mode="ED"` ran backwards
  in time**, and on every DMRG backend returned the conjugate of
  $\langle O\rangle$. Larmor precession of a +x chain in a field now gives
  $\langle S^y_0\rangle(t=1)=+0.4207$ on ED, where it gave $-0.4207$, and
  $O=S^z_0+iS^x_0$ on an eigenstate now gives $+0.5i$ on DMRG, where it
  gave $-0.5i$. A real Hamiltonian, a real start and a Hermitian
  observable with real matrix elements, which is what every earlier test
  measured, do not move (§7, findings 9 and 10).
- **A state set by hand now reaches the solver.** After `set_gs()` every
  DMRG correlator measured the session's own solved state and then
  overwrote the chain's with it (TD on a pure doublet member goes from
  2.04e-01 to 3.3e-06 off exact), and `set_initial_wf`/
  `set_initial_wf_guess` never reached the session, so the transverse-field
  Ising example plotted $M_z/n\approx0.001$ in the ordered phase where the
  seeded branch gives 0.4998 (findings 11 and 12).
- **`get_kondo_spectrum(mode="DMRG", n_gs>1)`** re-swept each member before
  measuring it, so on v2 and v3 an excited member split below `delta`
  relaxed into the lower one (1 of 6 runs at the average 1.5, now 6 of 6),
  left `gs_energy()` at the last member's energy afterwards, and was inert
  under `submode="EX"`, which now measures from the chain's own state
  (1.1392 to 1.8711 across seeds, now 1.512780 on every one; §17,
  findings 13 to 15).
- **Smaller moves**: the Kondo terms on a non-increasing `es` (3.035 to
  0.03183 off exact on a refined block appended after a coarse one,
  finding 17); a padded ground state evolved with `set_pad_bonds` cleared
  under `TDVP_GSE` at `tdvp_gse_sweeps=0` (0.4929 to 1.2e-15 off the
  unpadded run, finding 18); `mpsalgebra.disentangle_manifold` on a
  Hermitian operator the canonical proof cannot see (a non-orthonormal
  basis, Gram error 0.32, now 1e-15, and proven operators move at the
  gauge level, finding 1).

**And these now raise where they used to be silent:**

| Call | Was | Now |
|---|---|---|
| KPM with a pole outside the rescaled window (`kpm_scale` below 1/2) | DMRG spectra up to 109 times the true peak, ED integrals up to 1e6 | `RuntimeError` on every backend, the guard being 1.5 times the exact moment bound (findings 2 and 3) |
| `kpm_energy_truncate=True` on `mode="ED"` | the untruncated curve | `NotImplementedError`, also through `mode.py`'s fallback (finding 4) |
| `kpm_n_scale=1.5` on a non-Hermitian `mode="ED"` KPM | `TypeError`, where `mode="DMRG"` accepted it | accepted on both, neither reading it (finding 5) |
| a misspelled key, a method or `itensor_version` in `kpm_finite`'s `window_chain_kwargs` | ignored, the default spectrum bit for bit | `TypeError` (finding 6) |
| `name=(A,B)` with `i=`/`j=` in `get_dynamical_correlator` and every wrapper of it | the sites silently dropped | `TypeError` (finding 16) |
| `n_gs>1` under `submode="SECTOR"` | the single-state value | `NotImplementedError` (finding 15) |

### The 2026-09-24 third pass

A fifth hunt (`docs/audit_2026_09_24c_hole_hunt.md`, four lenses over the
single commit `867e2b4` that closed the second pass) recorded 18 findings,
six of them from that commit and twelve older, and every one is fixed; two
needed a rebuild of the compiled extensions. Each item below names its
finding.

**Results that are not comparable across this change.**

- **`gs_energy_fluctuation()` below full bond dimension, on every DMRG
  backend**, was set by truncating $H|\psi\rangle$ to `maxm`. It is now
  $\lVert(H-\langle H\rangle)|\psi\rangle\rVert$ with an uncapped
  application: on an 8-site Heisenberg chain solved and measured at
  `maxm=3` it goes from 6.7e-03 to 0.344 on `"python"` and from 1.5e-02 to
  0.347 on v3, and an exact state measured at `maxm=3` from 1.64 (v3) to
  1e-11. `mode="ED"` subtracts $\langle H\rangle$ first too, so its value at
  an exact eigenstate goes from the 1e-07 roundoff floor to 1e-15. With it,
  `gs_energy(maxde=...)` stops later, since it used to stop on an
  under-reported number, and returns the refined energy it leaves on the
  chain, -3.374933 where it returned the unrefined -3.279373; a
  correlator afterwards measures that refined state instead of re-solving
  at the original `maxm` (§3, findings 6 and 7).
- **After `set_gs()` of a state off the ground manifold** every submode on
  both modes measures it from its own energy. DMRG KPM and ED KPM, CVM,
  INV, ROOTN and TD move by $E_x-E_0$ (0.3 on the 3-site chain in a field,
  where the elastic line sat at +0.3), and ED `submode="ED"`, which read no
  state at all, now gives the set state's density, 2.834 away from what it
  returned on a 2.835 peak (§6, findings 1 and 2).
- **`promote_to_dense()` keeps the state the chain holds.** A state set in
  a sector and promoted with no read in between reads -0.957107, its own
  energy, where it read the sector ground energy -1.616025, and a carried
  state that is not the sector ground state is kept (-1.131526) where
  `"python"` re-swept it out of the sector to the global ground state
  -1.857107 (§3, findings 3 and 4).
- **`set_hamiltonian(H2, restart=False)`** now drops the stored ground state,
  so `gs_energy()`, `vev()`, `get_excited()` and the direct KPM moments answer
  for H2: -1.780099 where they gave H1's -1.616025 on a 4-site chain (finding
  5).
- **Weakly non-Hermitian operators.** The Hermiticity probe is scale-free,
  so an anti-Hermitian part below about 1e-2 in absolute size is no longer
  called Hermitian: a weak loss term keeps its decay rate (Im $E_0$ =
  -0.002134, where v2 and v3 returned 0 and `"python"` a real part 3 to 18
  per cent off), a non-Hermitian chain written in small units is solved as
  one, and `disentangle_manifold` diagonalizes such an operator rather than
  its Hermitian part (§2, finding 12).
- **`"python"` operators whose coefficients are all below about 2e-7**
  were built as a different operator, a lone `1e-7*Sz0` as the zero MPO;
  `vev(1e-7*Sz0)` goes from 0 to -1.8936e-08 (finding 13).
- **v3 `TDVP_GSE` from a start whose site 0 is one local basis vector** (any
  ladder operator or projector on site 0, or a product state) left site 0
  frozen: `evolution_ABA(A=C_0)` goes from 4.7e-01 to 1e-07 off exact, the
  public (Cdag_0, C_0) TD spectrum's peak from 0.081 to 0.47, and a Néel
  start is exact in 8 of 8 runs where up to 5 of 8 failed (§7, finding 14).
- **Every `mode="ED"` real-time route** (`evolve_and_measure`,
  `evolution_ABA`, `submode="TD"`) uses an exact propagator, and moves by
  1e-08 to 1e-06 on the small chains the tests use and up to 1e-03 at
  `dt=0.2` on an 8-site chain with an extensive energy (§7, finding 17).

**Faster, with no number changed.** The first read after `set_gs()` no
longer solves the upper band edge and an energy fluctuation it discarded
(7.1 s on a 24-site `"python"` chain for a 0.04 s result), and `"python"`'s
upper band edge caps its local Krylov dimension at 20, costing about one
ground-state solve where it cost 5 to 9 on a Heisenberg chain (findings 9
and 10). A malformed KPM call raises before any ground-state work again,
and a SECTOR call makes no solve on the caller's chain (finding 11).

**And these now raise, or now work, where they used to do something else:**

| Call | Was | Now |
|---|---|---|
| a misspelled keyword to `evolve_and_measure`/`evolution_ABA` on DMRG | ignored, the run at the defaults | `TypeError`, as on ED (finding 15) |
| `h=` to `evolve_and_measure`/`evolution_ABA` on ED | `TypeError`, "multiple values for argument 'h'" | evolves under `h`, as on DMRG (finding 16) |
| `submode="SECTOR"` after `set_gs()` of a state that is not its sector's ground state | the ground state's spectrum | `NotImplementedError` (finding 8) |
| `disentangle_manifold` on an ED manifold | `AttributeError` for every operator | the eigenbasis (finding 18) |
| an unknown keyword to `gs_energy_fluctuation()` on DMRG | forwarded and dropped | `TypeError` (finding 7) |

### The 2026-09-25 open items

A sixth pass (`docs/audit_2026_09_25_open_items.md`) took ten of the items
the five hunts had left open or recorded as unreviewed leads, reproduced
each on `8dd2198`, fixed it and handed every fix to a reviewer briefed to
refute it; one needed a rebuild of both compiled extensions. Each item
below names its entry in that record.

**Results that are not comparable across this change.**

- **A setting passed to a chain constructor** now takes effect: a 12-site
  Heisenberg chain built with `maxm=4, nsweeps=10` on v3 goes from
  -5.1420906326, the `maxm=30` answer, to -5.1323602278, the number the
  same request made by assignment returns, and a constructor `mode="ED"`
  answers every read that does not pass `mode=` itself by ED (the KPM ZZ
  peak of a 4-site Heisenberg chain in a field of 0.2, 0.205144 on
  `"python"` and 0.205151 on v3, both DMRG, becomes 0.205190). The same
  holds for `Thermal_Spin_Chain(..., mode="ED")`, which ran DMRG (§1,
  init-kwargs).
- **v2 and v3 Hamiltonians in small units** (largest coefficient below
  about $10^{-6}$): the ground state, the excited states, both band edges,
  the KPM spectra and `gs_energy_generalized` now match the unit-scale
  calculation, where a 6-site Heisenberg chain written as $sH$ returned the
  Néel $-1.25$ for $E_0/s$ against $-2.4936$ from $s\approx4\times10^{-7}$
  (v3) and $10^{-7}$ (v2) down. Every v2 or v3 calculation whose largest
  coefficient is below 1 now runs on $2^kH$, which moved nothing beyond the
  run-to-run noise at ordinary scales (§2, small-units).
- **Operators with coefficients below $10^{-8}$** keep their terms on every
  backend: `vev(eps*Sz0)/eps` at `eps=1e-9` goes from exactly 0 to
  -0.189361 on a 6-site chain in a field, and $E_0/s$ of a Hamiltonian at
  $s=10^{-8}$ from 0 to the right answer on ED and `"python"`. A term at or
  below $10^{-12}$ of its operator's largest coefficient is now dropped
  where the absolute rule kept it when it was above $10^{-8}$, and an
  anti-Hermitian part between about $10^{-10}$ and $10^{-8}$ of the largest
  coefficient now makes `gs_energy()` take the non-Hermitian route
  (§2, clean-threshold).
- **`"julia_live"` warm starts**: `set_initial_wf_guess()` of an exact
  doublet member goes from $|\langle t|\mathrm{gs}\rangle|^2=0.0008$ (0.0008
  to 0.68 over runs) to 1.0000, and `set_initial_wf(x)` followed by
  `gs_energy()` from the solved $-1.0$ to $\langle x|H|x\rangle$ (§3,
  julia-warm-start).
- **`get_gs(wf0=x, reconverge=False)` on a solved chain** returns x with
  `e0` = $\langle x|H|x\rangle$ (0.0388284989 where it returned the stored
  $-2.5231886435$ on a 6-site `"python"` chain), and every later reader
  measures x; `get_gs(wf0=x)` sweeps from x (§3, get-gs-wf0). Since the
  2026-09-25b pass the state returned and measured is $x/\lVert x\rVert$,
  and the call works on the routes that resolve to ED, where it raised
  `TypeError`.
- **Readers after a correlator that followed `gs_energy_generalized()`**
  measure the generalized state: on a Hermitian chain never solved before,
  `e0` goes from $-2.493577$ to $\lambda=-3.597994$ and $\langle
  S^z_0\rangle$ from 0 to $-0.4801$ (6 sites, $A=1+0.8S^z_0$); on a
  non-Hermitian chain with no correlator before the generalized solve,
  solved first or not, `e0` goes from $-1.596396$ to
  $\lambda=-2.17581-0.213083i$ (§3, generalized-cache). Since the
  2026-09-25b pass `e0` is the generalized state's own energy instead of
  $\lambda$, which is kept as `lam_generalized` (next subsection).
- **`tc.MBChain` after `Thermal_Spin_Chain.get_gs()` at $T>10^{-5}$**:
  on a 3-site chain at $T=1$, `gs_energy()` goes from $-2.25$ to $-0.416388$
  and the KPM sum rule from $-0.166370$ to $-0.069291$ on the DMRG
  backends, and on `mode="ED"` `vev(Sz0*Sz1)` from 0 to $-0.069398$ and
  `gs_energy()` from $-2.25$ to $-1.0$, the lowest eigenvalue (§9,
  thermal-bypass). The 2026-09-25b stepper moved these again, to
  $-0.428230$ for `gs_energy()`, $-0.071372$ for `vev(Sz0*Sz1)` (both
  exact) and $-0.071331$ (DMRG) and $-0.071188$ (ED) for the sum rule, and
  the switch is now at $T>10^{-5}W$, $W$ the spectral width.
- **`gs_energy(wf0=x, reconverge=False)` on a non-Hermitian chain** is
  $\langle x|H|x\rangle$ on `"python"`, v2 and v3, where it returned the
  NH-DMRG $-1.596396$ (§6, nh-injected).
- **`submode="ROOTN"` on the lower-level route** (`get_dynamical_correlator_MB`
  with a string name) honours `i=`/`j=`: at (1,1) on a 4-site `"python"`
  chain it goes from the (0,0) curve, $4.07\times10^{-2}$ off ED's (1,1) on
  a 0.103 peak, to $4\times10^{-16}$ (rootn-ij).

**And these now raise, where they used to do something else:**

| Call | Was | Now |
|---|---|---|
| a chain constructor keyword that is not a setting (a misspelling, the chain's state, a model-built operator list) | ignored | `TypeError` naming every offending key (init-kwargs) |
| `maxm=0` or an unknown `mode=` at construction | ignored | `ValueError`, as on assignment (init-kwargs) |
| an unknown keyword to `submode="ROOTN"` | ignored | `TypeError` (rootn-ij) |
| `get_gs(best=True, wf0=x)` on a solved chain | the stored state | `TypeError`, as on a fresh chain (get-gs-wf0) |
| non-Hermitian `submode="KPM"` after `set_gs()` | a spectrum built on the left state of an earlier solve | `RuntimeError` naming the missing left state (nh-injected) |
| `gs_energy(wf0=x)` on a non-Hermitian chain | the NH-DMRG energy, x dropped | `TypeError`; `reconverge=False` takes x (nh-injected) |

**And these now work:** on `"julia_live"`, `set_gs()` and
`gs_energy(wf0=x)` on a solved chain, which raised `AttributeError` and
`TypeError`; `gs_energy_generalized()` on a v3 MPO Hamiltonian, which raised
`AttributeError` (-3.59799449, the `MultiOperator` route's value); and the
first correlator after a plain NH solve no longer repeats the solve.

### The 2026-09-25b hole hunt

A sixth hunt (`docs/audit_2026_09_25b_hole_hunt.md`, four lenses over the
single commit `e7b1196` that fixed the open items above) recorded 29
findings and two unreviewed leads, and its fix pass found one more
(finding 30). Three of the 29 came from that commit and the rest are
older, reached by probing next to it, and seventeen are one family: a
threshold in absolute energy units, which the relative term cutoff of
`e7b1196` let a Hamiltonian or an operator written in small units reach
for the first time. Every one is fixed, and both compiled extensions were
rebuilt. Each item below names its finding; the record carries the
reproduction that was actually run, and every number is from it.

**Results that are not comparable across this change.**

- **A state whose norm is not one**, set with `set_gs()`,
  `set_initial_wf()` or `set_initial_wf_guess()`, or passed as `wf0=` with
  `reconverge=False`, is now normalized once, where it becomes the chain's,
  on every backend. Before, `e0`, the session's `vev()` and the energy
  fluctuation divided by $\langle x|x\rangle$ while KPM, CVM, TD, TDZ,
  ROOTN, `evolve_and_measure()` and, on `mode="ED"` and `"julia_live"`,
  `vev()` did not, so each of those now drops by exactly
  $1/\langle x|x\rangle$. On a 4-site Heisenberg chain at $x=2s$, $s$ its
  ground state: the KPM sum rule of $(S^z_0,S^z_0)$ goes from 0.999995 to
  0.249999, `evolve_and_measure(Sz0*Sz1)` at $t=0$ and ED's `vev(Sz0*Sz1)`
  from $-0.910684$ to $-0.227671$, `gs_energy_fluctuation()` on ED and
  `"julia_live"` from 9.696152 to 0 ($8\times10^{-9}$ on `"julia_live"`),
  and `"julia_live"`'s
  `get_excited_states(n=2)` from $[-6.464102, -0.957107]$ to
  $[-1.616025, -0.957107]$. `get_gs()` hands back $x/\lVert x\rVert$, and
  `get_rdm()` is the density matrix of the ray on all four backends, where
  it was $\rho/c^2$ for a state of norm $c$ (trace 0.25 at $c=2$). An
  explicit `wf=` to `evolve_and_measure`/`evolution_ABA` is left as given
  (findings 1 and 11).
- **Readers after `gs_energy(wf0=x, reconverge=False)` on a route that
  resolves to ED** (the chain's `mode="ED"`, a `mode="ED"` call, or v3's own
  fallback below three sites) measure x, where the ED ground state stayed
  on the chain: on a 4-site Heisenberg chain with $0.4S^z_0+0.3S^x_3$,
  `vev(Sz0)` goes from $-0.2286763037$ to $\langle x|S^z_0|x\rangle$,
  $-0.1089088536$ and $0.0390101221$ for the random states drawn on
  `"python"` and v3. The returned energy stays the lowest eigenvalue
  (finding 2).
- **A chain whose `mode` is `"DMRG"` no longer overrides a `mode="ED"`
  call**, so the ED cross-check is ED again: `gs_energy(mode="ED")` on an
  8-site Heisenberg chain in a staggered field at `maxm=2` goes from
  $-3.6734578613$ to the exact $-3.7040879103$, and
  `get_correlation_matrix(T=1)`, which is ED only, from occupations up to
  0.149 off the Fermi function to $5.6\times10^{-16}$, and from 4.2 s to
  0.1 s (finding 3).
- **`gs_energy(reconverge=True)` and `gs_energy(maxde=...)` on a solved
  chain** reach the solver, with `get_gs()` alike: the first sweeps from the
  stored state ($-4.2548001917$ to $-4.2548333556$ on a 10-site Heisenberg
  chain at `maxm=6` and one sweep, `"python"`), and the second refines
  ($-4.1431954920$ to $-4.2580352072$ at `maxm=3`, exact $-4.2580352073$)
  (finding 5, and the 2026-09-25 record's `maxde=` lead).
- **`"julia_live"` solves on a Hamiltonian whose couplings skip a site**
  leave the product state they were trapped in (§3): the ground energy goes
  from $-0.5$ to $-1.0$ on the Heisenberg chain on the even sites of 6 and
  from $-1.5$ to $-3.232051$ on an 8-site chain with next-nearest-neighbour
  exchange only, `Thermal_Spin_Chain` at $T=0$ from
  $\langle H\rangle=-0.5$ to $-1.0$, `gs_energy_generalized` from
  $-0.833333$ to $-1.496331$, the non-Hermitian solve from
  $-1.4997+0.0994i$ to $-3.219401$, and the lowest excited state from
  $-3.190943$ to $-3.232051$, all now exact. Every other `"julia_live"`
  solve moves within its run-to-run noise, and a warm 16-site solve costs
  about 40 per cent more (finding 7).
- **`gs_energy()` and `e0` after `gs_energy_generalized()`** are the
  generalized state's own energy (§4), where they were $\lambda$: on a
  4-site Heisenberg chain in a field of 0.3, from $-2.162930$ to $-1.385986$
  for $A=1+0.8S^z_0$, from $-0.808013$ to $-1.616025$ for $A=2\,\mathrm{Id}$
  and from $-1.115761$ to $-1.565921$ for $A=1.5+0.4S^z_0$, on `"python"`,
  v3 and `"julia_live"`, and on a non-Hermitian chain the biorthogonal
  energy of the pair. $\lambda$ is still the return value and is kept as
  `sc.lam_generalized`. Every KPM, CVM, ROOTN and TDZ spectrum afterwards
  moves rigidly by $\lambda-E_{w_g}$ onto the lines TD and EX already gave
  ($-0.777$ for $A=1+0.8S^z_0$, $+0.808$ for $A=2\,\mathrm{Id}$), and on
  `"julia_live"` the KPM window takes its lower edge from a solve of $H$,
  so the calls that raised "KPM moments diverging" return, and the
  $A=1+0.8S^z_0$ line's height goes from 1.2248 to 1.3037
  (findings 8 and 10).
- **`submode="CVM"` whenever $\eta\lVert B|\mathrm{GS}\rangle\rVert$ is
  small, and on exact solves at small $\eta$** (§6): on a 6-site Heisenberg
  chain with $0.3S^z_0$, every frequency of $(S^z_0,S^z_0)$ at an operator
  scale of $10^{-4}$ was the flat $\eta\langle AB\rangle/\pi$ (0.0159 after
  dividing by $\epsilon^2$) and is now the ED curve, and the same for
  $sH$ at $s=10^{-4}$; `get_kondo_spectrum(mode="DMRG", submode="CVM",
  order=2, T=0)` at its default `delta=2e-6` on three spins at
  $J=10^{-3}$ goes from about $10^{-8}$, flat, to the exact curve on the
  same grid (4.6981 at the ends, 0.6763 in the middle); the on-line point
  $\omega_0=0.52710639$ of that 6-site chain goes from
  $1.5915\times10^{-4}$ to 13.2636 at $\eta=2\times10^{-3}$ and from
  $1.5915\times10^{-5}$ to 132.6333 at $\eta=2\times10^{-4}$, as ED; and on a
  10-site chain on a 121-point grid, v3, 10 wrong points (up to 0.989 of
  the peak) become none at $\eta=0.05$, and 45 become none at $\eta=0.02$.
  At unit scale the tolerance tightens by $1/\lVert b\rVert$ (tenfold at
  $\eta=0.2$), so ordinary curves move within the old tolerance, the
  6-site one from $6.8\times10^{-7}$ to $6.7\times10^{-9}$ of the peak on
  `"python"`. A truncated solve still stops early and warns, and is not
  converged either way: a 20-site chain at `maxm=cvm_maxm=30` gives
  0.30333897 at $\omega=0.3$ where it gave the flat 0.01193662. The
  121-point grid takes 1214 s where it took 677 s, and a truncated point
  about twice as long. `applyinverse` stops relative to the norm of the
  vector it inverts against, so `CVM_explicit` on
  $(\epsilon S^z_0,\epsilon S^z_0)$ is $1.4\times10^{-6}$ of the peak off at
  every $\epsilon$ on v3, where it was 0.37 off at $\epsilon=10^{-4}$
  (findings 13 and 14).
- **Non-adjoint pairs under `CVM_explicit`, the non-Hermitian CVM and INV,
  and `cvm_solver="variational"`** return their own $C[A,B]$ (§6), where an
  absolute adjoint gate answered small pairs with another pair's density:
  `CVM_explicit` on $(10^{-2}S^x_0,10^{-2}S^y_1)$ from 2.004 of the peak off
  to $3.5\times10^{-5}$, the variational solver on
  $(10^{-2}S^z_0,10^{-2}S^z_3)$ from 2.124 to $5.9\times10^{-8}$, and ED's
  non-Hermitian CVM/INV on $(10^{-5}S^x_0,10^{-5}S^y_1)$ from 2.064 to
  $1.4\times10^{-15}$ (finding 16).
- **Every finite-temperature `Thermal_Spin_Chain` result** (§9). At
  $T=0.5$, on 10 sites ($\langle H\rangle$, `"python"` DMRG at `maxm=64`)
  from $-2.777892$ to $-3.166391$ (exact $-3.166396$), on 12 sites from
  $-3.287097$ to $-3.849203$ (exact $-3.849213$), on 6 sites on `mode="ED"`
  from $-1.674494$ to $-1.800863$, and on 3 sites from $-0.754223$ to
  $-0.769459$ (exact $-0.769460$), now the same at every scale and under a
  constant offset, where an offset of 20 gave 0.393848. On that 3-site
  chain in a field $B$ at $T=B$, $\langle S^z_{\rm tot}\rangle$ goes from $-0.089480$,
  $-0.013179$ and $-0.002545$ at $B=10^{-2}$, $5\times10^{-3}$ and
  $10^{-3}$ to $-0.231057$ at all three (exact $-0.231059$). A float $T$
  above 5 returns ($-0.071119$ at $T=5.5$) where it raised
  `ZeroDivisionError`, and a numpy scalar above 5 gets its temperature
  ($-0.055408$ at `np.float64(7.0)`) where it got $T=\infty$. The
  2026-09-25 thermal-bypass numbers move with it, `tc.MBChain.gs_energy()`
  at $T=1$ on 3 sites from $-0.416388$ to $-0.428230$ and `vev(Sz0*Sz1)`
  from $-0.069398$ to $-0.071372$. On `"julia_live"` a 4-site chain at
  $T=0.5$ comes out $2.3\times10^{-5}$ off where the other backends are
  $2\times10^{-8}$ off, since its band-edge solves returned a narrower
  width (not examined further). The cost grows with $\beta W$: 8.6 s where
  it was 4.9 s on 12 sites at $T=0.5$, and 19.8 s where it was 0.7 s for
  the field chain at $B=10^{-3}$ (findings 17 to 19).
- **`get_excited(n, mode="ED")` above 2000 states in small units** keeps
  every copy of a degenerate level: on a 12-site ferromagnet at
  $s=2\times10^{-7}$ the six-fold ground level came back with a magnon at
  $-2.715926s$ in place of a copy of $-2.75s$ in 8 of 8 calls, and the
  calls that stalled at $s=10^{-5}$ and $10^{-6}$ return. Below an infinity
  norm of 1 the levels also move at roundoff (finding 20).
- **`submode="ROOTN"` in small units** keeps its Lanczos basis, where it
  collapsed to the seed alone and returned one Lorentzian: the peak of
  $sC_s$ on a 6-site chain goes from 0.393371 to the unit-scale 0.157205 at
  $s=10^{-12}$ on ED and $s=10^{-10}$ on v3 (finding 21).
- **A KPM cross-correlator of two small operators** is no longer read as
  an autocorrelator: $C[\epsilon S^z_0,\epsilon S^z_3]/\epsilon^2$ at
  $\epsilon=1.2\times10^{-10}$ and $10^{-11}$ was $C[S^z_3,S^z_3]$, 2.070 of
  the peak off, and is now $C[S^z_0,S^z_3]$ to $10^{-14}$ on v2, v3 and
  `"python"` and to $7\times10^{-5}$, its run-to-run noise, on
  `"julia_live"`; `get_distribution` likewise, from 1.935 to $10^{-14}$
  (finding 23).
- **VUMPS in small units or next to a large constant**, on v3 and
  `"python"`, reports `converged=True` where it could not: TFIM at $D=8$
  from $s=0.1$ down to $10^{-8}$ and at constant offsets of $+10$, $+100$
  and $-100$ on v3; below about $s=10^{-6}$ the energy moves too, $e_0/s$
  from $-0.440126707917$ to $-0.440127030549$ at $s=10^{-8}$
  (finding 24).
- **v3 iDMRG in small units**: the energy density of TFIM is
  $4.8\times10^{-12}$ off the closed form at every $s$ down to $10^{-14}$,
  where it was 0.12 off at $s=10^{-13}$ and of the wrong sign at $10^{-14}$
  (finding 27).
- **NH-DMRG in small units or next to a large constant**, on v2, v3 and
  `"python"`, lands on the lowest level where it returned an excited one:
  $E_0/s$ on a 6-site Heisenberg chain with $0.3iS^z_0$ at $s=10^{-7}$ goes
  from $-0.9904-0.0568i$ to $-2.4610410218$ on v2, $E_0-c$ at
  $c=5\times10^5$ from $-1.9919-0.1109i$ to $-2.4610374186$ on `"python"`,
  and v3 at $s=10^{-13}$ from 0.570 off to $7\times10^{-15}$. At unit scale
  the tie window between Ritz values is about three times wider than it
  was, which moved nothing measured (findings 25 and 27).
- **v2 and v3 Hamiltonians whose duplicate terms cancel** are solved at the
  scale of the MPO they build (§2): $sH+S^z_0-S^z_0$ at $s=10^{-10}$ from
  0.09 to 0.53 off in $E_0/s$ on a 6-site Heisenberg chain to
  $5\times10^{-15}$. On v3 the XX chain at $J=1$, whose realified
  coefficients are 0.5, now also solves at scale 2, moving within its
  run-to-run noise (finding 28).
- **`"python"` operators whose largest coefficient is below 1** are built
  at unit scale (§2): a Heisenberg MPO at $s=10^{-7}$ goes from
  $1.1\times10^{-9}$ off relative to itself ($2.5\times10^{-4}$ at
  $s=2^{-40}$) to $10^{-15}$ at every scale, on numpy with OpenBLAS, and
  every `"python"` `MultiOperator*MPS` in small units moves with it,
  $\lVert(sH)|u\rangle\rVert/s$ reading 0.897242774928 at every $s$ down to
  $10^{-16}$, where it read 1.432 at $10^{-16}$ (finding 30).
- **`mpsalgebra.lowest_energy_arnoldi` in small units**, a non-default
  route: $E_0/s$ on a 4-site Heisenberg chain at $s=10^{-6}$ goes from
  $-1.6150192976$ to $-1.6160253791$, and at $s=10^{-7}$ with
  `delta=s*1e-3` from $+0.3510023868$ to $-1.6160254038$ (lead
  `scale-arnolditk-absolute-stop`).

**New attributes.** `sc.lam_generalized` holds $\lambda$ after
`gs_energy_generalized()`; `tc.anneal_step` and `tc.anneal_order` set the
thermal stepper (§9); `sites.SETTINGS` is the explicit list of the 36
settings a chain constructor accepts (§1); and on `mode="ED"` the state
gains `normalize(tol=)` and `norm()`, as on the MPS backends (§2).

**Quieter changes.** With `verbose` on, v2 and v3 print one line before
each solve that runs at $2^kH$, saying that the energies logged below it
are $2^k$ times those of the operator (finding 29); the returned energies
never were. NH-DMRG's convergence certificate is relative to the
Hamiltonian's units and constant, so a run it used to accept silently in
small units or next to an offset now retries up to `ntries` times and
warns with its "best relative residual" (finding 26). `State.normalize()`
on `mode="ED"` prints the MPS backends' warning where it was silent below
its floor (finding 22).

**And these now raise, or now work, where they used to do something
else:**

| Call | Was | Now |
|---|---|---|
| a constructor keyword naming one of nine attributes that are not settings (`hopping`, `hubbard`, `pairing`, `exchange`, `fields`, `resorder`, `resordered_indexes`, `hubbard_matrix`, `fit_td`) | accepted; `hubbard=2.0` added $2\,\mathrm{Id}$ to the Hamiltonian | `TypeError` naming the key, `Thermal_Spin_Chain` included; the last five are gone from the chain (finding 4) |
| a misspelled keyword (`wf=`, `reconverg=`) to `gs_energy()`, `get_gs()` or `get_excited(n=1)` on a solved chain | swallowed, the stored answer returned | `TypeError`, as on a fresh chain (finding 5) |
| a reader on a chain with no Hamiltonian | an error not naming the cause | `ValueError` saying to call `set_hamiltonian()` (finding 6) |
| `gs_energy(H=H2)` on a non-Hermitian chain | H2's pair stored as the chain's state | `TypeError` naming `nhdmrg(H=...)` (finding 9) |
| `get_gs(wf0=x, reconverge=False)` on a route that resolves to ED | `TypeError` | $x/\lVert x\rVert$ (finding 2) |
| an MPS as `wf0=`, or a misspelled keyword, to `gs_energy()` on a route that resolves to ED | ignored | `TypeError` (finding 2) |
| a set state whose norm is zero or not finite | `None` stored as the state | `ValueError` (finding 1) |
| an unknown single-site name on an ED spin chain (`ISy`, `Adag`) | `RuntimeError: No active exception to reraise` | `ValueError` naming it (finding 6) |
| `mode="ED"`, `submode="ED"` with the averaged manifold wider than `delta` | silent | `RuntimeWarning` (finding 15) |
| `Thermal_Spin_Chain` with a negative or NaN `T` set after construction | the ground state | `ValueError` (finding 19) |
| `get_rdm(i=ns-1)` on `"julia_live"` | `BoundsError` | the matrix (finding 12) |
| `get_rdm()` or the bond entropy under an ED mode on `"julia_live"` | the DMRG matrix, or `AttributeError` | `NotImplementedError`, as on the other backends (lead `session-julia-ed-guard-carveout-rdm-bond-entropy`) |
| `CVM_explicit`, or the non-Hermitian CVM/INV, on a non-adjoint pair | `NotImplementedError` | $C[A,B]$ (finding 16) |

Results from before this pass are not comparable wherever a number above
moved; everywhere else the calculation is the one it was.

### The 2026-09-26 Kondo exchange diagram

Not an audit: a check of the Kondo spectrum against the paper found that
its eq. 25 puts the third-order Kondo term's exchange-diagram log at
$eV+\epsilon_{im}$ where the second-order T-matrix of its own
Hamiltonian puts it at $eV-(\epsilon_f-\epsilon_m)$ (§17, "Where the
exchange log sits"). `get_kondo_spectrum` now uses the latter on both
`mode="ED"` and `mode="DMRG"`. Every third-order spectrum with an
inelastic transition out of an occupied state moves, by up to 5–20% of
the third-order term at $T=0$ and a few per cent of the total at 1 K
(the Fig. 7d 10 T overshoots from 1.226/1.183 to 1.248/1.203); zero-bias
values, a free $S=1/2$ at $B=0$, `order=2` and the potential term do not.
Results from before are not comparable. The lower-level
`kondospectrumtk.twotime.kondo_term_from_two_time` now takes
`(t, G, Gx)` triples and raises on the old `(t, G)` pairs, and its
kernel `K_W` is replaced by the one-sided `K_F`.

### The 2026-09-26 review of `"python"` real-time evolution

Not an audit: a check of `itensor_version="python"`'s real-time evolution
against exact propagation. The integrators themselves were right. Two-site
and one-site TDVP reproduce $e^{-iH\,dt}$ to $10^{-11}$ at full bond
dimension on spin-1/2 (long-range and Dzyaloshinskii-Moriya), spin-1,
Jordan-Wigner fermion and boson chains, at real and complex $dt$, and below
full bond dimension they track ED as closely as v3 does. Three things
around them changed numbers:

- **The Krylov exponentiator** behind every `"python"` TDVP route (TD,
  TDZ, `evolve_and_measure`, `evolution_ABA`, METTS, the infinite-chain
  window, the Kondo two-time construction) now stops when
  $|z|\,\beta_k\,|[e^{zT_k}]_{k1}|<10^{-10}$ for $e^{zH}$, an error
  relative to the evolved vector and free of units. It used to stop on
  $\lVert v\rVert\,\beta_k\,|[e^{zT_k}]_{k1}|$, which carries the units of
  $H$ and the norm of the state, and it returned a step its 50-vector
  budget could not take unconverged (0.31 off at $|z|$ times the spectral
  width equal to 100); such a step is now split. This closes the "absolute
  Krylov error goal" the 2026-09-25 record left open, on `"python"`:
  $C[\epsilon S^z_0,\epsilon S^z_3]/\epsilon^2$ under `submode="TD"` was
  $8.8\times10^{-2}$ off $C[S^z_0,S^z_3]$ at $\epsilon=10^{-8}$ and $0.24$
  at $10^{-9}$, and is $10^{-11}$ at both, and a Néel quench under $sH$ at
  $dt/s$ was $0.48$ off ED at $s=10^{-12}$ and is at the $2.9\times10^{-6}$
  of $s=1$. At ordinary units results move below $10^{-10}$. v3's
  `applyExp` is vendored ITensor and keeps its absolute goal.
- **`evolve_and_measure` and `evolution_ABA`** on `"python"` (TDVP,
  TDVP_GSE, TEBD) and on v3 TEBD now put the state back to its initial
  norm after every step, as v3's TDVP and every `"julia_live"` loop
  already did. Real-time evolution conserves the norm, so what the
  truncation removes from it is error, and it used to scale every later
  $\langle\psi(t)|O|\psi(t)\rangle$: on a 12-site Néel quench at
  `maxm=8`, $\langle\psi|\psi\rangle$ fell to 0.982 by $t=5$ under TDVP and
  to 0.970 under TEBD, and $\langle H\rangle$ drifted by
  $4.9\times10^{-2}$ where v3 drifted by $1.8\times10^{-4}$ (now the
  same $1.8\times10^{-4}$ on both). Nothing moves where nothing truncates,
  and the `quench` loops behind `submode="TD"` always renormalized.
- **`tevol_method="TDVP_GSE"` on `"python"`** now normalizes each Krylov
  vector $H^k|\psi\rangle$ before choosing the new directions, as v3's
  expansion does. The choice used to depend on $\lVert H\rVert$: that
  8-site Néel quench was $1.0\times10^{-4}$ off ED at $s=1$ and
  $1.3\times10^{-3}$ at every $s\le10^{-4}$, and is $2.7\times10^{-5}$ at
  all of them.

Results from before are not comparable where these apply. The
regressions are in `tests/test_pyitensor_time_evolution_review.py`. One
thing was measured and left as it is, on every backend: under
`TDVP_GSE` the bond dimension stops growing after `tdvp_gse_sweeps` steps
(3 by default), so a quench whose entanglement keeps growing saturates.
On the 12-site Néel quench the largest error up to $t=5$ stays near
$2\times10^{-3}$ from `maxm=16` to 64, where `"TDVP"` reaches $9\times10^{-7}$ at
`maxm=64`. Raise `tdvp_gse_sweeps` to the number of steps, or use
`"TDVP"`, for such a run.

### The 2026-10-14 documentation pass

These bugs were found while the user guide was rewritten. The regressions
are in `tests/test_docs_pass_2026_10_14.py`. Four fixes change numbers:

- **DMRG on `itensor_version=2` and `3` now does 4 Davidson steps per
  two-site solve**, where it did ITensor's default of 2. With 2 steps a
  chain whose Hamiltonian conserves a charge the sites do not (QNs are off
  by default) could stop short of its ground state at the default
  `nsweeps=15`. A 4-site spinless chain with hopping and a
  nearest-neighbor interaction gave between -1.6215 and -1.7015600 over
  repeated runs, against ED -1.7015621187. It now gives the ED value to
  1e-15 in every run, and a 30-site chain of the same kind moved by about
  1e-6. Every v2 and v3 ground-state energy can move a little. A sweep
  costs about 1.5 times more. `"python"` already solved each two-site
  problem to convergence and is unchanged.
- **The chain returned by `meanfield.spinchain_meanfield` reports the
  mean-field energy.** Its Hamiltonian now includes the constant
  $-(1-p)\sum_{ij}J_{ij}\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle$.
  Before, its `gs_energy()` counted the decoupled exchange twice: -2.5
  against -1.25 for a 6-site Heisenberg chain in the Néel state. The
  state and the magnetization do not change.
- **`get_fidelity(fmode="PT", n=...)` uses `n`.** The sum used to run over
  the default number of states of `get_excited_states()` whatever `n` was:
  2 states on DMRG and 10 on `mode="ED"`. Now it is `n` states on both,
  with the default `n=2`. An ED result computed with the default changes.
- **`wf.get_dm(inds=...)`** used the first `len(inds)` sites whatever
  `inds` said. Only `inds` of the form `range(k)` were right.

These used to fail and now work:

| Call | Was | Now |
|---|---|---|
| `Parafermionic_Chain(n, Z=2)` | process abort on v2/v3 (also with `mode="ED"`), `ValueError` on `"python"` | works on ED, v2, v3 and `"python"` (the transverse-field Ising chain) |
| `get_charge_gap(mode="ED")` | `TypeError` | works, agrees with DMRG |
| `get_dynamical_correlator(mode="ED", submode="ED", dex=...)` after `set_gs()` | `TypeError` | works; `dex` has no effect on a set state |
| `get_full_SS_correlator(mode="ED")` | left the last state of the ground manifold set on the chain | the chain is left as it was |
| infinite chain `gs_energy()` that does not converge | only `converged=False` | also a `RuntimeWarning` with the reason |

The last one matters for a J1-J2 chain at $J_2=0.4$ on a one-site cell.
Its ground state is dimerized, so no uniform state of one site is a fixed
point. On `"python"` VUMPS at `maxm=4` it gave +0.006, -0.001, +0.154 and
-0.050 in four runs, and -0.380 at `maxm=30`, all unconverged. On a
two-site cell it converges to -0.3801177 at `maxm=4`.

### The 2026-10-15 effective Hamiltonian

- **`get_heff` and the two fits in `effectivehamiltonian.py` take a
  `projector=`, the model space of the effective Hamiltonian.** With it,
  the $n$ lowest eigenstates are projected on the model space and Löwdin
  orthonormalized, and the fit is $H_\mathrm{eff}=\sum_nE_n|\tilde n\rangle\langle\tilde n|$
  (des Cloizeaux), which has exactly the eigenvalues $E_n$.
  `effectivehamiltonian.single_occupancy_projector(fc)` is the singly
  occupied space, where the spin operators of a Hubbard chain act. Without
  `projector=` the call is unchanged and returns the same numbers as
  before: the operators are fitted inside the span of the eigenstates.
- **The exchange of the Hubbard dimer changes when the projector is
  passed**, which the user guide and its figure now do. Without it the
  singlet carries doubly occupied weight, where $\mathbf S_0\cdot\mathbf S_1=0$,
  and the fit gives $J=(E_T-E_S)/(1/4+3w/4)$, with $w$ the singly occupied
  weight of the singlet. With it, $J=E_T-E_S$, on ED and on DMRG:

| $U/t$ | without projector | with projector | $E_T-E_S$ |
|---|---|---|---|
| 1 | 2.18109 | 1.56155 | 1.56155 |
| 4 | 0.93064 | 0.82843 | 0.82843 |
| 20 | 0.19949 | 0.19804 | 0.19804 |

- The check for linearly dependent states in `orthonormalize` is now
  relative to the largest eigenvalue of the Gram matrix (1e-10), not a
  test for a non-positive one. A manifold larger than the model space
  raises a `ValueError` that names the model space.

### The 2026-10-15 infinite chains

No returned number moves. Two kinds of state now give a `RuntimeWarning`
when an observable is read off them
(`tests/test_audit_2026_10_15_infinite.py`).

- **iDMRG on a one-site cell can end in a dimerized state.** On the
  spin-1/2 Heisenberg chain the central bond of the growing chain
  alternates between a strong and a weak bond from one step to the next
  ($3\langle S^z_0S^z_1\rangle=-0.4721$ and $-0.4140$ at `maxm=16`), so the
  fixed point of iDMRG at finite `maxm` has a period of two sites. The
  extraction of the cell is exact: the cell reproduces both bonds of the
  growing chain, AKLT at `maxm=2` and Majumdar-Ghosh on a two-site cell
  give exact correlators on both backends, and on `n_uc=2` the mean of the
  two bonds is the energy per site to $4\times10^{-7}$ at `maxm=16`, the
  truncation error (unchanged from `etol=1e-10` to `1e-13`). On `n_uc=1`,
  `correlator("Sz",0,"Sz",1)` reads only the first bond of the cell
  ($3\langle S^z_0S^z_1\rangle=-0.4137$ against $e=-0.4431$ at `maxm=16`),
  and the energy per site alternates between macro-iterations, so
  `converged` stays `False`. `gs_energy`, `vev` and `correlator` now warn
  when the Schmidt spectra of the two bonds of the cell differ, with
  their entanglement entropies (0.909 and 0.968 at `maxm=16`). The
  alternation of the two-site cell decreases with `maxm`: bonds $-0.4970$
  and $-0.3887$ at `maxm=8`, $-0.4719$ and $-0.4143$ at 16, $-0.4607$ and
  $-0.4256$ at 32. VUMPS on the same two-site cell at `maxm=16` finds a
  lower energy ($-0.44311147$ against $-0.44308908$), with nearly equal
  bonds and Neel order instead.
- **An unconverged VUMPS state gives wrong correlators while its energy
  looks right.** On the Heisenberg chain on one site at `maxm=16` both
  backends give $e=-0.44284628$, while the sum of the bond correlators is
  $-0.391$ on `itensor_version=3` and $-0.467$ and $-0.434$ in two runs
  on `"python"`, with $\langle S^z\rangle\neq0$ and a gauge mismatch of
  1.5. VUMPS does not converge on this one-site cell; the reason was not
  investigated. `vev`, `correlator` and `excitation_energies`
  now warn when `converged` is `False`, with the gauge mismatch.

### The 2026-10-15 fidelity susceptibility

- **`get_fidelity(fmode="derivative")` on a ground state that is
  degenerate within `delta` changed.** It returned $|\det U_d|$, with $U_d$
  the second difference of the overlap matrix of the manifold, which scales
  as $\chi^g$ for $g$ states. It now returns the second difference of
  $|\det U|$, the same form as for a single state, which is the sum of
  $\chi$ over the $g$ states:
  $\sum_{a\in\mathrm{GS}}\sum_{m\notin\mathrm{GS}}|\langle m|H_1|a\rangle|^2/(E_m-E_0)^2$.
  With a single ground state nothing changes. Open Ising chain,
  $h_0=4\sum S^zS^z$, $h_1=2\sum S^x$, `dl=0.02`, $\chi/n$:

| case | before | now |
|---|---|---|
| ED, $n=12$, $\lambda=0.6$, `delta=1e-4` (one state) | 0.19304 | 0.19304 |
| ED, $n=12$, $\lambda=0.6$, `delta=1e-2` (doublet, split 2.8e-3) | 0.3454 | 0.3421 (= 0.1930 + 0.1490, even plus odd) |
| ED, $n=8$, $\lambda=0$, exact doublet, total $\chi$ | 0.7657 | 1.74996 (exact 1.75) |
| DMRG v3, $n=24$, $\lambda=0.6$, default `delta` | 0.4357 | 0.2694 |

- The two states of the doublet do not have the same $\chi$ on an open
  chain, so the new value is not twice the value of one state. A call with
  `delta` below the splitting, on DMRG, still gives a meaningless number
  (55 at $n=24$, $\lambda=0.6$, `delta=1e-12`): DMRG lands on an arbitrary
  superposition of the two states, and that path did not change.
- An unknown `fmode` raises `ValueError`; it returned `None`.

### The 2026-10-22 tangent-space Krylov route

- **`get_kondo_spectrum(mode="DMRG")` without a `submode` changed on
  `itensor_version=3` and `"python"`.** The default was `"KPM"`, with the
  second-order and potential terms from the broadened dynamical correlator
  and the Kondo term from two-time TDVP. It is now `"TASK"`, the exact sums
  of `mode="ED"` over the tangent-space Krylov spectrum of the ground state
  (arXiv:2510.07279), with no frequency or time grid. The result is exact
  when `maxm` saturates the Hilbert space (`cutoff=0`) and `nkry` exhausts
  the Krylov space. On `itensor_version=2` and `"julia_live"` the default
  is still `"KPM"`.

| case | before (`"KPM"`) | now (`"TASK"`) |
|---|---|---|
| 3-site chain of the Kondo guide page, against ED | 0.5% off, 119 s | 1e-15, 0.02 s |
| S=1 impurity and 4 spins, saturated, against ED | broadened by `delta` | below 1e-9 |
| 10-site S=1 impurity chain, saturated, `nkry=400` | not run (two-time TDVP too slow) | 2e-11 |

- With the implicit default, the grid keywords of the old route (`es`,
  `delta`, `dt2`, `n_t2_half`, `dtau`, `n_tau_half`) raise `TypeError`,
  whose message says to pass `submode="KPM"`. An explicit `submode="KPM"`
  gives exactly what the default gave before.
- On a chain too short for the session (v3 below 3 sites, `"python"` on 1
  site), the TASK route uses the exact T=0 ED spectrum, as the other DMRG
  routes do.
- `get_dynamical_correlator(submode="TASK")` is new on both backends;
  nothing it replaces changed.

## The 2026-10 audit

`docs/audit_2026_10_hole_hunt.md`. No fix changed a number that was
returned before. `taskdmrg.get_poles` gained the `missing_weight` entry and
a warning; the potential term's coverage warning now fires on the spinful
fermion chains; `TaskKondoSpectrum` warns when handed a state above the
ground state. `get_kondo_spectrum(order=3)` with a correlator submode on a
chain too short for the session returns the exact ED Kondo term where it
raised `AttributeError`, and raises `NotImplementedError` on v2 where it
raised `AttributeError`. A non-integer `nkry` raises where it was cut to an
integer on the TASK correlator.

## The 2026-10-23 audit

`docs/audit_2026_10_23_hole_hunt.md`. Two calls return different numbers
now, and in both the old number was wrong:

| Call | Chain | Before | After |
|---|---|---|---|
| `excited_vev(Sz[0]*Sz[1], n=3)`, every backend and `mode="ED"` | 4-site XXZ, Jz=0.7 | the ground-state value, once, whatever `n` | one value per state, equal to ED to 1e-5 |
| `get_bond_entropy` of `2*wf` on v2 and v3 | 8-site XXZ, Jz=0.7 | 2.08 to 3.52 above the exact entropy | exact to 1e-10 |

The bond entropy of a unit-norm state, every solver output, is unchanged.
Four calls raise where they returned a wrong number: `vev(O, T>0)` on DMRG
(`NotImplementedError`; it returned the T=0 value), `vev` with an unknown
keyword on DMRG (`TypeError`), `Spin_Chain.get_correlator(name=...)`
(`TypeError`; it returned <Sz Sz> for every name), and `metts_vev` with a
non-Hermitian basis operator, or with an all-diagonal basis whose charge H
conserves (`ValueError`; `("N",)` on a hopping chain returned the value of
one particle-number sector with an error bar of 0). `Fermionic_Chain.excited_vev(mode="ED")`
works where it raised a matmul `ValueError`.

Afterwards, while fixing the broken examples. `submode="maxent"`,
`submode="CVMimag"`, `get_distribution_maxent` and
`reconstruct.reconstruct_distribution` return a result where they raised
`NotImplementedError`, since `dmrgpy.maxenttk` and `dmrgpy.padetk` now
exist. `"CVMimag"` returns `-Im G/pi`; its code had `-2 Im G/pi`, twice the
weight of every other route, which no earlier result could have used since
the route never ran. `Spinon_Chain.gs_energy()` returns the projected
energy (-2.4935771 on 6 sites, equal to the Heisenberg ED) where every call
raised `AttributeError`.

Afterwards, the reviewed leads (the record's last section, findings 6 to
11). Two calls return different numbers now, and in both the old number was
wrong:

| Call | Chain | Before | After |
|---|---|---|---|
| `get_site_entropy`, `get_pair_entropy`, `get_mutual_information`, every backend and `mode="ED"` | any state with a site or pair density-matrix eigenvalue at or below 1e-6, e.g. a 4-site chain in a staggered field with XY coupling 0.003 | the -p ln p of each such eigenvalue dropped: 5.625e-7 | exact: 8.657e-6, at most about 1.4e-5 higher per dropped eigenvalue |
| `gs_energy` and everything after it on `"python"` | two equivalent spins on a `SpinBoson_Chain` (exact swap symmetry inside the first two-site block) | an excited level, -4.470573 in 11 of 12 runs | the ground state, -5.917370 |

Every other fresh `"python"` solve agrees with the old one to its
convergence: the random part that now enters each Lanczos start is there
on the first sweep only, from a generator of its own, so numpy's global
stream and every seeded start are where they were. A warm start takes no
noise and is unchanged. Calls that raise now
where they returned a number: `get_pair_entropy`/`get_mutual_information`
are scale invariant where they raised (or, at norm^2 within 1e-3 of 1, were
up to 3e-4 off per 1e-3); `gs_energy_generalized` with a non-Hermitian
metric (`ValueError`) or one whose `<psi|A|psi>` goes to 0 or below
(`RuntimeError`; A=-1 returned the largest eigenvalue, 2.4936 on a 6-site
Heisenberg chain, and on v3 a singular A aborted the interpreter); and
`metts_vev` with an all-diagonal basis on a chain whose conserved charge H
keeps (`("Nupdn",)`, `("F",)`, `("N2",)`, ...: one charge sector, error bar
0). An NH-DMRG state that did not converge prints its warning again on
every reuse.

## The 2026-10-25 Davidson solver on `"python"`

The two-site solve of `"python"` DMRG is now a Davidson solver with a
diagonal preconditioner, at least 4 steps per bond, as v3 does. Before, it
was Lanczos with up to 200 steps per bond. A ground state in a conserved
sector is much faster: on a 4-cell c-f chain on native spinful sites,
sector (8, 0), `maxm=300`, 20 sweeps on one core, 323 s before and about
20 s now (15 s on v3). Excited states keep a larger budget, 20 steps with
a 1e-8 residual stop. The spectral-bound solve is unchanged.

| quantity | before | now |
|---|---|---|
| E0 of the sector chain above, against v3 (-4.081749584520) | 7e-11 | 5e-12 |
| energy variance of an exact gapped product state (Néel state of a staggered field) | 7e-15 | 2e-12 |
| energy variance of a 16-site Heisenberg chain at `maxm=64` | 9.8e-6 | 9.7e-6 |

So the energy is unchanged at its convergence. The eigenvector of an exact
state is converged a little less, to about 1e-6 in norm, since each bond
stops once its energy is stable. A truncation-limited state is unchanged.

## The 2026-10-29 non-Hermitian real-time evolution on `"python"`

Real-time evolution on `"python"` with a non-Hermitian Hamiltonian was
wrong. Every TDVP route exponentiated the local generator with a Lanczos
recursion, and that recursion is correct only for a Hermitian operator.
This covers two-site TDVP, `tevol_method="TDVP_GSE"`, and `Chain.tdvp_step`
as TDZ and the two-time Kondo construction use it. Now a non-Hermitian MPO
takes an Arnoldi exponentiator. A Hermitian Hamiltonian still takes the
Lanczos one, so its results do not change. Measured on a 4-site Heisenberg
chain with $0.5 S^x_0 + 0.7i S^z_0$ added, $dt=0.05$, 40 steps, as the
largest error of $\langle S^z_2(t)\rangle$ against the exact normalized
evolution:

| backend | before | now |
|---|---|---|
| `"python"`, two-site TDVP | 3.6e-3 | 9.0e-7 |
| `itensor_version=3`, two-site TDVP | 9.0e-7 | 9.0e-7 |

The remaining 9e-7 is the time-step error that the two backends share.
`mode="ED"` returns the unnormalized $\langle\psi|O|\psi\rangle$ for a
non-Hermitian $H$, while both DMRG backends normalize the state after each
step. So ED and DMRG still answer different quantities there; this is not
changed yet.

## The 2026-10-30 Lindbladian spectrum on `"python"`

`Lindblad_Chain.get_spectrum` with random starts now begins its second
branch from a traceless state, as in arXiv:2509.07709. The exact evolution
conserves the trace, so this branch no longer contains the steady state and
resolves the slowest decaying modes instead. A seeded start
(`start=(psi, phi)`) is already traceless, so its numbers do not change.
Measured on the driven, dephased 4-site XXZ chain at the defaults, as the
distance of each ED eigenvalue to the nearest Krylov one:

| eigenvalue | before | now |
|---|---|---|
| steady state | 1.0e-6 | 1.3e-6 |
| gap mode, $-0.1439$ | 2.2e-4 | 8.9e-5 |
| $-0.2281 - 1.6373i$ | 1.9e-4 | 6.2e-5 |
| $-0.2281 + 1.6373i$ | 2.3e-4 | 3.0e-4 |

The magnetization profile of the steady state is unchanged, 2.9e-6 off ED.

`get_krylov_spectrum` on a Hermitian generator no longer runs the contour
probe and takes the angle $\pi/8$ directly. The probe could only return
$\pi/2$ there, so the angle, and every number, is the same unless the
probe misfired on an unconverged slowest eigenvalue. On a non-Hermitian
generator the probe now bisects its angle grid instead of scanning it. The
two give the same angle, since instability grows monotonically with the
angle; they could differ only where the heuristic test itself is not
monotonic.

The non-Hermitian Arnoldi exponentiator now stops on its error estimate
relative to the result rather than to the input vector. A unitary step is
unchanged. A dissipative step that shrinks the vector now takes more
iterations. On the 5-site Bose-Hubbard Lindbladian at dt=0.5 the old test
stopped some local steps after one iteration with a result 5e-19 of the
input. `Lindblad_Chain.get_spectrum` also projects the trace out of its
second branch after every step and sub-steps its TDVP when a step does not
conserve the trace. On chains where neither triggers, as on the driven
4-site XXZ chain, every number above is unchanged. On the seeded 3-site
Bose-Hubbard chain the occupations go from 4.8e-8 to 1.6e-8 off ED, and on
the 4-site one each of the slowest eight eigenvalues comes closer to ED,
by factors 1.1 to 10. The 5-site check gave no steady state at all before (top
eigenvalue +0.485). It now gives $\lambda_1$ 7e-4 off ED and the
occupations 1e-3 off, with spurious modes still in the list; see
`docs/clik_mps_lindblad_plan.md`.
