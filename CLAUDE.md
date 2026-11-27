# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository. It holds what every session needs: how to build, how
to test, the architecture, and the rules the audits taught. The detail that
only matters for one part of the codebase is in a rule file under
`.claude/rules/`, loaded only when a file it covers is read, written or
edited:

- `cpp-backend.md`: the pybind11 extension, `mpscpp2/` and `mpscpp3/`.
- `pyitensor.md`: the pure-Python backend, its MPO builder and its GPU rule.
- `infinite-chain.md`: iDMRG and VUMPS on both backends.
- `kondo.md`: `kondospectrumtk/`, YSR, the NRG routes and ESR-STM.
- `multioperator.md`: the canonical form and the Hermiticity proof.
- `session.md`, `thresholds.md`, `dynamics.md`, `non-hermitian.md`,
  `thermal.md`, `mps-algebra.md`: the bodies of the durable rules of the
  audits, by the part of the code each one guards.
- `install.md`, `packaging.md`, `benchmarks.md`: `install.py` and
  `installtk/`, the PyPI wheel, and `benchmarks/`.

A rule file loads when the Read, Write or Edit tool touches a matching
file, not when a shell command prints it, so a session that reads with
`cat` or `sed` opens the matching one by hand before touching those files;
a pointer at the place each used to live says where it went. The history of
every fix, with its repro and its reviewer, is in `docs/audit_*.md`, and
the ledger of every result that changed is `docs/numerical_changes.md`.

## What this is

DMRGPY is a Python library for quasi-one-dimensional spin chains, fermionic
systems, parafermion and bosonic models using matrix product states (MPS).
It exposes a single Python API (`src/dmrgpy/`) backed by interchangeable
solvers per calculation:

- **DMRG**, executed in-process via a pybind11 extension over ITensor
  (C++, either v2 or v3, see below), via a pure-Python reimplementation
  (`itensor_version="python"`), or via a live Julia/ITensors.jl backend
- **ED** (exact diagonalization), a pure-Python fallback used for small
  systems, for cross-checking DMRG results, and automatically whenever the
  C++ extension is not compiled

The Python layer is what this project actually develops. The C++/ITensor
backends under `src/dmrgpy/mpscpp2/ITensor/` (ITensor v2) and
`src/dmrgpy/mpscpp3/ITensor/` (ITensor v3) are vendored copies of the
ITensor library maintained upstream: treat them as black boxes, not code to
read or modify here. Read them only if some feature does not work with the
`mpscpp2`/`mpscpp3` code and you need to figure out why.

## Install / build

There is a `pyproject.toml`, but it is only for building the PyPI
distribution (see "Packaging / PyPI" below); it is not how this repo is
developed. For development the project is used from `src/` via a path added
to `PYTHONPATH` (or a symlink into `site-packages`), which is what
`install.py` sets up; there is no `setup.py` and no editable install.

```bash
python install.py                    # compiles ITensor v3 (mpscpp3) + its pybind11 extension
python install.py --itensor-version=2    # compiles ITensor v2 (mpscpp2) instead
python install.py --itensor-version=both # compiles both v2 and v3 backends
python install.py --gpp=g++-6        # use a specific compiler (needs g++ >= 6, LAPACK/BLAS)
python install.py --doctor           # only the requirement checks, no multi-minute ITensor build
python install.py --extension-only   # after a C++ edit: rebuild only the extension, not ITensor
python install.py --build-only       # build and self-test, change nothing outside this checkout
python -m dmrgpy.doctor              # which backends load in this session, and how each was built
python install_julia.py              # alternative: Julia backend instead of/alongside C++
```

`install.py` first checks every build requirement by a real trial compile
and link (`installtk/requirements.py`), and only then compiles the vendored
ITensor and the `pybind` target of the backend that `--itensor-version`
selects (`installtk/install2.py`). The compiler and the LAPACK/BLAS are
each found by a cascade of candidates, a conda compiler first under a conda
Python, and `pybind11` is installed with `pip` if it is missing. The
build ends with a self-test in a fresh process, and ITensor is rebuilt
only when its sources or the configuration changed (`--force`). The
details are in `.claude/rules/install.md`, which loads with `install.py`,
`installtk/`, `cppext.py`, `buildinfo.py` and `doctor.py`. There is no
lint config, and the only CI is the build of the user guide.

A C++ change is not tested until the `.so` has been rebuilt
(`python install.py --extension-only`); dmrgpy warns at import when the
`.so` is older than its sources, and when one is there and does not load.
To test the installer itself, use a worktree and `--build-only`: a plain
`python install.py` rebuilds the working `.so` and links into
site-packages.

## Packaging / PyPI

`pyproject.toml` builds a **pure-Python** distribution (`python -m build`):
the ED backend, the `itensor_version="python"` backend and the Julia
bridge, with no C++ at all. Users who want the compiled backends clone the
repo and run `install.py`. The release steps and what is load-bearing in
`pyproject.toml` are in `.claude/rules/packaging.md`, which loads with
`pyproject.toml`, `release_pypi.py`, `.gitignore` and `juliapkg.json`.
Three of its rules matter outside those files:

- `numba` is a core dependency, not an extra: `algebra/kpmextrapolate.py`
  imports it at module level, on the import path of every chain.
- Every module in the wheel must byte-compile, and the wheel says
  `OS Independent`: the ED and `"python"` backends run on Windows, so a
  shell command in their code path breaks that.
- `.gitignore` has an unanchored `julia*` pattern, so a new file or folder
  whose name starts with `julia` is ignored unless it is un-ignored, as
  `src/dmrgpy/juliapkg.json` is.

## Running / verifying changes

`tests/` is a real pytest suite. Run it with `python run_tests.py` from the repo root (a thin wrapper: `pytest tests`
from the repo root, then `clean.py`) or plain `pytest tests`.
`tests/conftest.py` inserts this checkout's own `src/` at the front of
`sys.path`, so pytest always exercises the current working tree's code,
never a stale `site-packages`/symlinked checkout (see the worktree pitfall
below). `tests/_helpers.py`'s `energy_ed_v2_v3`/`vev_ed_v2_v3` are the
standard way most tests cross-check ED against DMRG on the C++ backends (a
`versions=` kwarg lets a test opt a specific `itensor_version` out; it
covers `(2, 3)` only, and the pure-Python backend is tested by the dozens
of files that parametrize their own `itensor_version` instead).
`tests/reference_data.py` holds shared golden values pinned to a specific
historical commit (its own module docstring explains how to regenerate
them); several test files inline their own golden values as local
constants instead, which is equally fine for values used in one place. When
a test's assertion is against exact DMRG/ED agreement (not a golden
regression value), prefer `pytest.approx(..., abs=1e-6)`-style tolerances
over exact equality, since DMRG is iterative. The audits' regressions live
in `tests/test_audit_<date>_<cluster>.py`, one file per fix cluster; each
`**Status**` line in the matching `docs/audit_*.md` record names the test
that pins it.

**`itensor_version="julia_live"` tests only need to run when a change
touches Julia code** (`mpsjulialive/`, any `*.jl` file, or a shared entry
point like `manybodychain.py`/`infinitechain.py` that dispatches into that
backend). These tests carry a large, mostly fixed cost unrelated to what
they are testing: `juliacall`'s JIT compiles each distinct Julia function
signature the first time it is called in a process, and this compile cost
dominates wall time regardless of how small the test's own parameters are.
Confirmed directly by isolating a `julia_live` `metts_dynamical_correlator`
call down to `nsamples=1`, `nwarmup=1`, `nt=1`: it still took ~73s, almost
entirely JIT, against ~27s of genuine algorithmic work for the same call at
its normal (nsamples=40) parameters once the function was already warm.
Reducing test parameters further cannot buy back that fixed cost, and it is
paid independently by every distinct Julia function signature exercised
across the whole suite, which is why the full suite's slowest handful of
tests are almost all `[julia_live]`-parametrized. For a change that does not
touch Julia code at all, skip them locally with
`pytest tests -k "not julia_live"` (the parametrize id `julia_live` appears in every such test's
node id); run the full suite including `julia_live` before finishing work
that does touch Julia code, or before a release.

The `examples/` directory (100+ self-contained scripts, one per physical
model or feature, organized into thematic subfolders) is a second,
complementary regression surface: broader and more exploratory than
`tests/` (most scripts just print and plot rather than assert), but
covering more ground (every physical model, every dynamical-correlator
submode, every backend combination including `"python"`). Some
`examples/*/*/main.py` scripts do use real `assert`s and are meant as
regression tests in the same spirit as `tests/` (e.g.
`time_evolution/tdvp_VS_ED_time_evolution`,
`backend_comparison/backend_switch_consistency`,
`finite_temperature/thermal_purification_VS_exact`,
`staticcorrelators/static_correlator_VS_ED`,
`entanglement/entanglement_entropy_VS_ED`,
`dynamical_correlator/dynamical_correlator_VS_ED`). To run one:

```bash
cd examples/readme_examples/static_correlator_S12 && python <script>.py
```

`examples/readme_examples/` mirrors the snippets shown in `README.md`.
`examples/backend_comparison/` and the `v2_VS_v3_*`-named examples inside
each theme folder directly compare the ITensor v2 and v3 C++ backends
against each other (same script, both `itensor_version`s, small systems),
the fastest way to check a `mpscpp3` change did not diverge from
`mpscpp2`'s numerics. One rule if you use one as a template: pin the sweep
schedule inside the script instead of inheriting the library defaults
(`boson_models/v2_VS_v3_boson` inherited them and its `assert` fired
non-deterministically; it now pins `nsweeps=80`/`maxm=100`, 2026-09 audit
#24). After running examples, `python clean.py` recursively removes
generated working directories (`.mpsfolder`, `.pychainfolder`,
`.dmrgfolder`) and stray `ERROR`/`*.OUT` files from the tree.

**Examples should plot, not just print/assert.** What sets `examples/`
apart from `tests/` is that a human is expected to look at the result, so
every `examples/*/*/main.py` ends with a `matplotlib` plot of whatever it
computed whenever the script produces a sequence of values; the scripts
that also carry a regression `assert` keep it and plot afterwards. A script
whose entire output is a single scalar is the one exception, and a last
resort: sweep a cheap axis (chain length, a coupling, time, bond dimension,
backend against backend) instead. An audit found 105 scripts that had
settled for print-only output with such an axis already computed.
**Examples are short and simple** (the author's rule since 2026-10-23): 40
lines, 50 at most, written in plain English the way a condensed matter
physicist who is not a programmer would write them, about the physics and the
one capability shown, nothing else. The `new-example` skill
(`.claude/skills/new-example/SKILL.md`) has the full rules, the two
approved models, the path line, and how to run a script against this
worktree.

**Worktree/symlink pitfall**: `site-packages/dmrgpy` is typically a symlink
into one specific checkout's `src/dmrgpy/` (usually the primary working
directory, not whichever worktree you are in;
`python3 -c 'import dmrgpy; print(dmrgpy.__file__)'` shows which). Since every
`examples/*/*/main.py` does `sys.path.append(...)` (append, not insert),
that symlinked package resolves before the appended path when running a
script directly with plain `python3 main.py` from a different worktree,
silently testing the wrong checkout's code with no error. Force the intended
checkout explicitly with
`PYTHONPATH=/path/to/this/worktree/src python3 main.py`, or
`sys.path.insert(0, ...)` in a one-off driver script.
`tests/conftest.py` avoids this problem entirely for anything run through
pytest.

**A real, confirmed `mpscpp3` bug**: ITensor v3's two-site `dmrg()` (see
`chain_session.h`'s `dmrg_args()`) aborts the whole process (`SIGABRT`,
`"LocalOp is default constructed"`, an internal ITensor v3 check in
`itensor/mps/localop.h`) for any chain with fewer than 3 physical sites,
independent of physics type. `mode.py`'s `get_mode()` falls back to ED
automatically for `itensor_version==3` with `self.ns<3` (same mechanism as
the "extension not compiled" fallback), so this no longer crashes; a few
older tests under `tests/` predate that fix and still manually exclude v3
for 2-site chains via `_helpers.py`'s `versions=` kwarg (the dimer tests in
`test_spin_chain.py`/`test_fermion_chain.py`, `test_entanglement.py`'s
Bell-pair test). Those exclusions are redundant defense-in-depth now, not
something to "fix" by removing.

### Cross-backend performance benchmarking

**Pin BLAS threads before timing anything**
(`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 ...`, set before numpy is
imported), or timings across runs are not comparable. DMRG is a great many
small dense linear-algebra calls, where more threads are slower, and under
oversubscription this dominates every other effect: on a 14-core host with
another job holding 10 cores, a 9.4s `"python"` solve became 28.2s and a
1.6s v3 solve 26.7s. dmrgpy never changes thread counts on its own. The
same cap applies to agents: at most three Python processes at once across
all subagents of a hunt or workflow (`run3.sh` in the 2026-09-24b record's
"Shared helpers").

`benchmarks/run_benchmarks.py` answers "which backend is fastest, and by
how much, for each kind of calculation" on every backend available, and
writes `benchmarks/output/report.pdf`. Its catalogue is the `CALCS` list of
`benchmarks/calctk.py`. The full notes are in `.claude/rules/benchmarks.md`,
which loads with `benchmarks/` and `src/dmrgpy/blasthreads.py`.

## The audits, and what to know before touching the code

Ten hole hunts, a fix pass and one review ran between 2026-08 and
2026-11-21, each hunt a multi-lens search for numbers that are plausible and
wrong, with every finding reproduced, handed to a reviewer briefed to refute
it, and fixed:
`docs/audit_2026_08_hole_hunt.md` (21 findings), `docs/audit_2026_09_hole_hunt.md`
(36), `docs/audit_2026_09_24_hole_hunt.md` (16), `docs/audit_2026_09_24b_hole_hunt.md`
(18), `docs/audit_2026_09_24c_hole_hunt.md` (18), `docs/audit_2026_09_25_open_items.md`
(ten open items, a fix pass rather than a hunt) and
`docs/audit_2026_09_25b_hole_hunt.md` (29, all fixed on 2026-09-26 together
with its leads), followed by the 2026-09-26 review of `"python"` real-time
evolution (regressions in `tests/test_pyitensor_time_evolution_review.py`,
measurements in `docs/numerical_changes.md`, leads ranked in
`docs/time_evolution_improvement_plan.md`), and by
`docs/audit_2026_10_hole_hunt.md` (7, over the TASK correlator and the
Kondo spectrum on fermion chains, all fixed on 2026-10-23), and by
`docs/audit_2026_10_23_hole_hunt.md` (four bugs and one performance item,
over the least-covered areas rather than a commit; the bugs fixed the same
day, and its nine leads reviewed the same day into six more findings, all
fixed), and by `docs/audit_2026_11_hole_hunt.md` (17 findings over the
Kondo machinery, bugs and measured hot spots together, 16 fixed on
2026-11-12; its rules are in `.claude/rules/kondo.md`), and by
`docs/audit_2026_11_21_hole_hunt.md` (22 findings over the 158 commits
since that one, hunted by read-only agents with every repro run by the
lead, all fixed on 2026-11-21; three only in part: finding 1 on DMRG is a
guard rather than a closing stage, finding 3's DMRG half warns, and
finding 4 warns on DMRG and keeps the free-spin band on the Krylov
routes). One more record is not a hunt:
`docs/audit_2026_11_15_hole_hunt.md` holds four defects found by the
2026-11-15 backend benchmark, one fixed at the size it was found at and
three open (see "Left open" below). Four things to
read before changing the kind of code they cover:

- **`docs/numerical_changes.md`** before comparing any number against
  a result from before 2026-09-26. It is the ledger of every fix that moved
  numbers rather than behaviour, by call and by size, and most of the fixes
  did: the TD/TDZ pi factor and the KPM calibration, every injected state
  (normalized now), every constructor keyword (applied now), every
  finite-temperature result, every `gs_energy_fluctuation`, every `CVM`
  spectrum, and more. The records' `**Status**` lines carry the same
  entries finding by finding.
- **`docs/development.md` section 4.10** before adding a dispatch. Almost
  every finding of the first hunt was a dispatch decision taken before the
  information that should inform it: a short circuit ahead of a cache key,
  a precondition ahead of the branch it qualifies, an `else` serving both
  "unsupported here" and "you typo'd it", a `**kwargs` with no consumer.
- **`src/dmrgpy/dynamics.py`'s module docstring** before adding a
  correlator submode. It is the single normative statement of the
  convention: every submode returns the complex Lehmann density
  `i(G^R-G^A)/(2pi)`, unchanged wherever the Lehmann weights
  `M_n = <GS|A|n><n|B|GS>` are real (`Im M_n == 0` is the discriminant,
  which a Hermitian pair implies but does not exhaust).
- **The `hole-hunt` skill** (`.claude/skills/hole-hunt/SKILL.md`) before
  running or extending an audit; the record shape and the evidence standard
  are fixed there.

The durable rules the hunts left, one line each with the test that pins it.
The body of a rule, with the mechanism and the numbers, is in the rule file
named above its group. Open that file and the finding before undoing one,
by hand if the session reads with `cat` or `sed`.

Ground state, session, dispatch and static observables
(`.claude/rules/session.md`):

- `restart()` is a cold start on every backend; only `restart=False` or an
  explicit start keeps a warm state (`test_restart_drops_session_state.py`).
- A state becomes the chain's normalized (`groundstate.unit_copy`), and a
  zero state raises (`test_audit_2026_09_25b_session.py`).
- Every correlator starts from `ground_state_on_session`, so it measures a
  state the caller set (`groundstate.mark_injected`) and not the session's
  own (`test_audit_2026_09_24b_groundstate.py`,
  `test_audit_2026_09_24c_session.py`).
- `groundstate.stored_answer_holds` is the one condition under which
  `gs_energy`/`get_gs` return a stored answer
  (`test_audit_2026_09_25b_session.py`).
- `minimum_energy()` never sweeps a state the session holds without its
  energy (`test_audit_2026_09_24c_session.py`).
- An energy variance is `vev.energy_variance`, never `<H^2>-<H>^2` through
  `vev(npow=2)` (`test_audit_2026_09_24c_fluctuation.py`).
- A chain constructor keyword is a setting of the allowlist
  `sites.SETTINGS`, assigned after `initialize()`; anything else raises
  (`test_audit_2026_09_25_construction.py`,
  `test_audit_2026_09_25b_construction.py`).
- `mode.resolve_mode` lets only a chain pinned to `"ED"` override a call's
  `mode=` (`test_audit_2026_09_25b_session.py`).
- The DMRG `vev` has no `**kwargs`: an unknown keyword raises, and so does
  `T>0` (`test_audit_2026_10_23.py`).
- `gs_energy_generalized` refuses a metric it cannot use, by a Python
  exception on every backend, and the entropy of a site or a pair is that
  of the normalized state (`test_audit_2026_10_23_leads.py`).
- `julia_live` starts every DMRG solve from a random MPS, so its solves
  vary run to run by solver noise (`test_audit_2026_09_25b_julia.py`).
- One predicate, `correlationentropy._full_applies`, picks the default
  `dmmode` in a sector and outside; do NOT add a condition to one branch
  only (`test_audit_2026_11_21_corrmat.py`).

Thresholds (`.claude/rules/thresholds.md`):

- Every threshold is relative to a scale of the problem it acts on. An
  operator whose largest coefficient is below 1 is brought into [1,2) by a
  power of two and is left byte-identical at 1 or more; do NOT also scale
  inside `Chain::nhdmrg` (`test_audit_2026_09_25_scale.py`,
  `test_audit_2026_09_25b_{cvm,nh,cpp,ednum,mpobuilder}.py`).
- A weak term is covered by a cutoff of the MPO build that is one per
  bond; do NOT edit the vendored svdMPO for this
  (`test_audit_2026_11_21_mpo_cutoff.py`).

Dynamical correlators and real-time evolution
(`.claude/rules/dynamics.md`):

- KPM broadening is calibrated: every KPM route goes through
  `algebra/kpm.py::polynomials_for_broadening`, so `delta` is the
  broadening it is everywhere else; `get_distribution()` is deliberately
  not on it (`test_audit_2026_09_24_kpm.py`, `test_audit_2026_09_24b_kpm.py`).
- `submode="TD"`/`"TDZ"` return the two-sided transform from forward runs
  only (`timedependent.lehmann_density_from_one_sided`); do not take
  `np.real` instead (`test_audit_2026_09_24_realtime.py`).
- Real-time evolution on `"python"` has one stopping test,
  `_krylov_verdict`, shared by the host loop and the device path; every
  `evolve_and_measure_*` loop restores the input norm; a non-Hermitian MPO
  takes the Arnoldi route, chosen by `mpsalgebra.mpo_is_hermitian(H)` and
  never by `dt` (`test_pyitensor_time_evolution_review.py`,
  `test_nh_time_evolution.py`).
- `mode="ED"` real-time routes are `expm_multiply` and start from the
  cached ground state when called without `wf=`
  (`test_audit_2026_09_24c_realtime.py`).
- TASK's error bar is two numbers: the two-site variance of a pole and
  `missing_weight`, the part of `B|GS>` outside the one-site tangent space
  (`test_audit_2026_10_task.py`).

Non-Hermitian Hamiltonians (`.claude/rules/non-hermitian.md`):

- NH-DMRG sweeps the transposed left state, H^T|wl> = E|wl>, on v2, v3 and
  `"python"`; do not go back to the adjoint left state
  (`test_audit_2026_11_15_nh.py`).
- A non-Hermitian `gs_energy()` on DMRG is imaginary time first
  (`nh_gs_method`, `nhlevels.py`), with its pair checked by the certificate
  of `nhdmrg()` (`test_audit_2026_11_21_nh.py`, `test_nh_levels.py`).
- ED orders a non-Hermitian spectrum by real part first
  (`algebra.eigen_order`), and takes the left state by its eigenvalue
  (`test_audit_2026_11_21_nh.py`).

Finite temperature (`.claude/rules/thermal.md`):

- `Thermal_Spin_Chain.anneal()` is order-8 Taylor steps at `dtau*W <= 2`
  (`test_audit_2026_09_25b_thermal.py`).
- `Thermal_Spin_Chain.get_gs()` returns its stored state only under
  `thermal.purification_key`; do NOT go back to a short circuit on
  `computed_gs` alone (`test_audit_2026_11_21_thermal.py`).
- A collapse basis METTS cannot sample with is refused
  (`mettsvev._check_basis`), which refuses every fermion chain (the
  author's pick); do NOT go back to a list of diagonal names
  (`test_audit_2026_10_23*.py`, `test_audit_2026_11_21_metts.py`).

MPS algebra (`.claude/rules/mps-algebra.md`):

- A BiCGSTAB solve checks its answer (`mpsalgebra.check_bicgstab`): it
  warns above 10 delta and raises above 1. Do NOT loosen it to quiet a
  test, and do NOT make `sum_mps_weighted(method="exact")` the default
  (`test_bicgstab_residual_check.py`, `test_audit_2026_11_21_exact-arith.py`).

Pure-Python backend (`.claude/rules/pyitensor.md`):

- `"python"` DMRG mixes noise into its Lanczos starts (`_noisy_start`), on
  the first sweep of a fresh start only; do NOT remove it and do NOT widen
  it (`test_audit_2026_10_23_leads.py`).

Kondo, NRG and ESR-STM (`.claude/rules/kondo.md`):

- The NRG tail of `nrgmps` has the same number of sites at every band D,
  and its hand-over check is relative to the largest hopping of the head
  (`test_audit_2026_11_21_nrg.py`).
- The ESR Krylov basis is closed under the Korringa operators on ED, and
  both Krylov routes end with a guard that raises naming `submode=None`
  (`test_audit_2026_11_21_esr_krylov.py`).
- The ESR band follows the drive on ED (`esrstm._ed_band`); DMRG warns
  instead, and the Krylov routes keep the free-spin band
  (`test_audit_2026_11_21_esr_window.py`).
- The Krylov Kondo stop reads residual intervals: a Ritz state is in the
  window when theta - r <= top (`test_audit_2026_11_21_kondo_krylov.py`).

**Out of scope by construction** for every hunt: the vendored ITensor and
TDVP folders, the three legacy quirks `.claude/rules/cpp-backend.md` lists
as deliberately reproduced (`evoloperator()`'s z^3/6 term on `H2`, the
`"moise"` key, the unreachable `"tevol_fit_td"` branch), the open
`docs/known_issue_*.md` items, and what `ROADMAP.md` marks as absent by
design.

Left open and recorded as such, in the records' "Left open" and "New leads"
sections: three of the four findings of
`docs/audit_2026_11_15_hole_hunt.md`, from the 2026-11-15 backend
benchmark and not yet fixed (the `"python"` and `"julia_live"` Sz=0
sector solves going wrong between 56 and 80 sites at `maxm=128`, and
`"python"`'s `randomMPS` overflowing from about 300 sites at
`maxm=128`), and the lead of its finding 3, fixed at 10 sites: on the
PT-symmetric XX chain, several eigenvalues with one real part, one
NH-DMRG attempt in three still stalls at 12 sites and the attempts stall
at 14 (ITensorNHDMRG.jl does worse), while the imaginary-time route,
now `get_nh_levels()` and the default of `gs_energy()`, finds the
distinct levels of smallest real part up to 16 sites (9 of 9 there on
`"python"`, 8 of 9 on v3 at maxm=30; its stored pair passes the
certificate at maxm=64, at maxm=30 neither route does and the better
pair is kept with a warning);
the `kpm_energy_truncate` window problem outside its worst regime
(`docs/known_issue_kpm_energy_truncation_window.md`); the sliver just below
`kpm_scale=1/2` where the elastic line is distorted while the moments are
still under the bound; the truncated CVM regime (warned, about 2x the cost,
still unconverged); `Thermal_Spin_Chain`'s cost and its variational band
width; the MPO-set Hamiltonian and
small-unit real-time evolution; v3's Krylov exponentiator's absolute error
goal (`applyExp`) and v3's run-to-run `emax` noise, which is the whole of
the ED-against-v3 KPM residual; `kpm_finite`'s `window_chain_kwargs` check;
the reason for
`julia_live`'s slowdown and its KPM guard raising `juliacall.JuliaError`
rather than `RuntimeError`; `ShiftInv`'s `applyinverse(delta=)`; and the
short `nt` defaults of `td_dynamical_correlator`. The variational
correction vector (`cvm_solver="variational"`) raises from about 16 sites
on a Kondo chain, where its compressed MPO of (H-E0-omega)^2+eta^2 is not
positive; the fix is open
(`docs/known_issue_ddmrg_compressed_operator.md`).
The fixes of `docs/audit_2026_11_21_hole_hunt.md` left open: the ESR
DMRG guard runs only when H conserves a charge that S^+_k changes; the
ESR Krylov routes keep the free-spin band and neither grow nor warn; the
`julia_live` side of the exact BiCGSTAB check was not measured; the cost
of the exact check and of the ESR closing stage at real size was not
measured; and v3's exact `bicstab` recurrence costs 2 to 4 times more
where `maxm` limits the solve, with the divergence under truncation still
open.

## Hard problems go to Fable

When a task feels too hard for Opus, say so and suggest giving it to a
Fable agent (the Agent tool with `model: "fable"`), rather than pushing on
or settling for a weaker answer (the author, 2026-11-26: "if something
feel too hard for opus, suggest giving it to fable"). Typical cases are a
derivation that decides a threshold or a bound, a numerical failure with
no clear mechanism, and the design of a new algorithm. The suggestion
names the task and its cost; the author decides. Such an agent keeps the
usual cap of at most three Python processes at once across all agents.

## Architecture

### Entry point classes

Each physical model is a thin subclass of `Many_Body_Chain`
(`manybodychain.py`), which holds all shared state (bond dimensions, sweep
count, cutoffs, mode, etc.) and implements the generic operations
(`set_hamiltonian`, `gs_energy`, `vev`, `get_gs`, correlators, time
evolution, entanglement, ...). Model modules just define the local Hilbert
space and convenience operators:

- `spinchain.py` -> `Spin_Chain` (aliased directly to `Many_Body_Chain`)
- `fermionchain.py` -> `Fermionic_Chain` (spinless fermions, `C`/`Cdag`/`N`/Fermi-string `F`)
- `spinfermionchain.py`, `bosonchain.py`, `parafermionchain.py`: same pattern for other statistics/models

### Operator representation: `MultiOperator`

Hamiltonians and observables are built as `MultiOperator` objects
(`multioperator.py`), sums of products of named single-site operators (e.g.
`Sx[i]*Sx[j]`). This is a backend-agnostic intermediate representation: the
same `MultiOperator` is later either converted to ITensor terms for the
C++/Julia/`"python"` DMRG backends (`MultiOperator.to_terms()`), or
converted to a sparse matrix via `multioperator.MO2matrix` for the ED
backend (`edtk/edchain.py`). `multioperatortk/` holds supporting machinery
(Jordan-Wigner strings for fermionic operators, static/long-range operator
construction, the canonical form below).

**Canonical form and the Hermiticity proof.**
`multioperatortk/canonical.py` rewrites every term into one spelling, so
that terms which are the same operator collect, and
`MultiOperator.simplify()`/`is_zero()`/`is_hermitian()`/`is_antihermitian()`
all read off it. The proof is **one-sided**: an empty canonical form of
`H-H^dagger` means Hermitian, and a nonempty one means not proven, so
`Many_Body_Chain.is_hermitian()` falls back to its random-witness probe; do
not turn that fallback into a rejection. Do NOT widen the same-site
reordering beyond two diagonal factors, and do NOT give the `A`-type names
no parity. The reasons are in `.claude/rules/multioperator.md`, which loads
with `multioperator.py` and `multioperatortk/`.

### Backend dispatch (DMRG vs ED, C++ vs Julia)

`mode.py` decides, per call, whether a calculation runs via DMRG or ED
(`get_mode`/`resolve_mode`): DMRG unless `self.mode` pins the chain to ED,
or unless `self.itensor_version` is `2` or `3` and the corresponding
pybind11 extension is not compiled (see `cppext.available(version)`), in
which case it falls back to ED with a printed note. Two more automatic
fallbacks exist, both structural: ITensor v3's two-site `dmrg()` aborts
below 3 sites (see "A real, confirmed `mpscpp3` bug" above), and
`"python"`'s two-site DMRG has no update to make on 1 site, where `dmrg()`
returned `None` before the fallback existed, so `itensor_version==3` with
`ns<3` and `"python"` with `ns<2` route to ED. `itensor_version="python"`
never falls back for a missing extension (it has none).
`groundstate.py::gs_energy_generalized` needs its own copy of the short-chain
guard (it has no ED fallback to be routed to), and there it raises, because
the symptom was worse than `None`: the outer self-consistent iteration
returned the Rayleigh quotient of a state no sweep ever touched, a silently
wrong number (-0.3049 for an exact -0.5); that guard sits before the
non-Hermitian dispatch, since NH-DMRG escapes v3's abort by never calling
`dmrg()` but its own sweep is two-site too. A conserved sector
(`set_conserved_sector`) makes `get_mode` raise on DMRG with
`itensor_version=2`, the one backend with no quantum numbers, rather
than answer with the global ground state; falling back to
ED is correct there, since ED restricts to the sector itself.

**Which backend a chain gets when the caller names none** is
`cppext.default_backend()`, not `DEFAULT_ITENSOR_VERSION` directly: the
default C++ version when that extension is compiled, and `"python"`
otherwise. That fallback is what makes `pip install dmrgpy` usable; before
it every default-backend chain in a pip install silently ran ED, which
cannot reach the sizes an MPS backend exists for. Only the implicit choice
moves: an explicit `itensor_version=3` with no extension still falls back
to ED, and `mode="ED"` is untouched.
`Many_Body_Chain.__init__`/`Mixed_Spin_Fermion_Chain.__init__` take
`itensor_version=None` to mean "pick one for me". See
`tests/test_default_backend_without_cpp.py`. Most public methods on
`Many_Body_Chain` accept a `mode="DMRG"|"ED"` kwarg so results can be
cross-validated between solvers.

- **DMRG, C++ (`itensor_version=2` or `3`, the default)**: entirely
  in-process through a pybind11 extension over a vendored ITensor, one
  `Chain` session per `Many_Body_Chain` held as `self._session`. Both
  versions expose the identical `Chain` session API and public
  `Many_Body_Chain` surface; `itensor_version` only picks which compiled
  extension (`mpscpp2._dmrgcpp` or `mpscpp3._dmrgcpp`, via
  `cppext.get_backend(version)`) backs the session.
  `cppext.DEFAULT_ITENSOR_VERSION` (currently `3`) is the single source of
  truth for the default, threaded through `Many_Body_Chain.__init__`,
  `setup_cpp()`, `get_backend()`/`available()`, and `install.py`'s own
  `--itensor-version` default (`installtk/__init__.py`'s matching
  constant); change it there, not per call site.
  `Many_Body_Chain.setup_cpp(version=2)` switches an existing chain between
  them. Everything else about this backend, the `cpp_handle` pattern, the
  unit-scale rule, the three deliberately reproduced legacy quirks, the
  `mpscpp3`-specific behaviour (QNs off by default and a random-MPS start,
  the opt-in conserved-sector mode, `apply_mpo()`'s `.noPrime("Site")`,
  `innerC`, TDVP as the default real-time method, `py::module_local()`), is
  in `.claude/rules/cpp-backend.md`, formerly this file's "In-process
  pybind11 extension" section.
- **DMRG, pure Python (`itensor_version="python"`)**: `pyitensor/`
  (`Many_Body_Chain.setup_python()` switches an existing chain to it) is a
  from-scratch, pure-Python/NumPy/SciPy reimplementation of exactly the
  ITensor v3 API subset `mpscpp3/chain_session.h` uses: Index/ITensor tensor
  core, the same eleven site types, an AutoMPO term compiler with its own
  from-scratch Jordan-Wigner threading, MPS/MPO algebra, two-site
  Lanczos-based DMRG (ground and excited states via the overlap-penalty
  method), two-site Krylov-based TDVP, and a `pyitensor/chain.py::Chain`
  exposing the identical method surface as the compiled backends'
  `self._session`, so every `itensor_version in (2,3)` call site elsewhere
  treats `"python"` as a third, always-available option wherever such a
  site was updated to include it. Exists so DMRG/TDVP work with zero
  compiler/pybind11 dependency, at the cost of being substantially slower
  than compiled ITensor (no block-sparsity, no JIT); see
  `pyitensor/__init__.py`'s docstring and the module docstrings throughout
  `pyitensor/` for the specific simplifications taken and why each one does
  not affect dmrgpy's own results. The MPO builder
  (`pyitensor/mpobuilder.py`, a finite-state machine over the terms'
  partial products) has two rules that are easy to get wrong, and a new
  `pyitensor` method takes its array namespace from `pyitensor/backend.py`
  so that it also runs on the GPU (the author's rule of 2026-10-24); both
  are in `.claude/rules/pyitensor.md`, which loads with `pyitensor/`.
- **DMRG, Julia (`itensor_version="julia_live"`)**: a live, in-process
  Julia session (`mpsjulialive/`, via `juliacall`/`juliasession.py`) with
  its own parallel set of modules (`mpsjulialive/groundstate.py`, `vev.py`,
  `mps.py`, `mpo.py`, ...) mirroring the top-level ones but talking to
  Julia's ITensors.jl instead of the C++ pybind11 extension.
  `Many_Body_Chain.setup_julia()` switches a chain to this mode. Since
  2026-11-08 its finite-chain surface is close to v3 and `"python"`
  (sectors, SECTOR, TASK, NH-KPM, the explicit resolvents, MPO and MPS
  algebra, every site type); `ROADMAP.md`'s Julia column is the list, and
  `docs/development.md` 4.6 the architecture. A parity-odd operator (a bare
  `C`) goes to OpSum as the Jordan-Wigner-dressed `A`/`Adag`/`F` names of
  `MultiOperator.to_terms()` (`fermion_ops.jl`), since OpSum refuses a
  parity-odd fermionic term; everything else keeps the raw names, whose
  Jordan-Wigner string ITensors.jl threads itself.
  `itensor_version="julia"` (a separate, older subprocess-based Julia path
  via `juliarun.py`) is not reachable through the normal public API today
  and should be treated as legacy/inert.
- **ED path**: `edtk/edchain.py` (`EDchain`) builds dense/sparse operators
  directly in Python/NumPy/SciPy (`pyfermion/`, `pyspin/`, `pyboson/`,
  `pyzn/` provide the per-statistics many-body operator construction,
  `edtk/one2many.py` promotes single-site operators to the full Hilbert
  space) and diagonalizes with `scipy.sparse.linalg`.

**Conserved sectors.** `Many_Body_Chain.set_conserved_sector(Nf=6)`/
`(Sz=0)`/`(Nf=8,Sz=0)` (off by default, no-args clears it) confines the
whole calculation to one quantum-number sector, by four mechanisms behind
one API: QN-carrying ITensor indices on `itensor_version=3` and, through a
second, QN-carrying ITensors.jl site set beside the dense one
(`mpsjulialive/sites_qn.jl`, `sector.py`), on `"julia_live"`; dense storage
with a charge grading plus a charge penalty on `"python"`
(`pyitensor/sector.py`), and a restriction of every assembled operator to
the basis states carrying the charge on `mode="ED"` (`edtk/edchain.py`,
where `promote_to_dense`/`promote_mps` are deliberately not implemented
since an ED re-solve is cheap). The sector lives on the chain, is re-applied
wherever a session is rebuilt (`sites.py::initialize`, `__deepcopy__`) and
is part of `groundstate.py`'s send-cache key. It pays off with size, and
only with size (0.6x at n=20, 4.2x at n=60 on Heisenberg chains); do not
quote a single speedup number. `tests/test_sector_conservation.py` and
`tests/test_sector_conservation_ed.py` cross-check every sector against ED
restricted to that sector; the reference must restrict, not filter
full-spectrum eigenvectors by `<N>`, since sector-degenerate eigenvalues
come back as arbitrary superpositions. See `docs/development.md` 4.3a and
the C++ rule file for the three ITensor behaviours that shape the v3
implementation.

### Supporting `*tk` packages

Functionality is generally split into a top-level module (the public method)
and a `<topic>tk/` package with the implementation, e.g.
`correlator.py`/`fermionchaintk/staticcorrelator.py`,
`dynamics.py`/`dynamicstk/`, `entropy.py`/`entropytk/`,
`mpsalgebra.py`/`mpsalgebratk/`, `entanglement.py`. When changing behavior
for a specific feature, look for both the thin dispatch module and its `tk`
counterpart.

### KPM / dynamical correlators

Dynamical correlators and generic operator distributions are computed with
the Kernel Polynomial Method (`kpmdmrg.py`, `Chain::kpm_dynamical_correlator`/
`Chain::general_kpm` in `mpscpp2/chain_session.h`/`mpscpp3/chain_session.h`
on the C++ DMRG side, `algebra/kpm.py` on the ED side) rather than direct
spectral decomposition, since exact diagonalization of the full spectrum is
infeasible for large chains. The convention is in the module docstring of
`dynamics.py`, and the calibration in `.claude/rules/dynamics.md`.

`submode="SECTOR"` (`sectordc.py`) is the exception to that shape, and the
only submode that does not compute a spectrum on the chain it is called on:
it solves the reference sector and the one sector `B|gs>` lands in (N+1 for
`Cdag`, Sz+2 for `S+`) on an internal clone, promotes both to dense indices
and contracts them, summing an explicit Lehmann representation over
converged per-sector eigenstates. It rests entirely on `promote_mps`
rebasing onto the chain's original `dense_sites_`, which
`set_conserved_sector` never rebuilds
(`tests/test_sector_promotion.py::test_states_promoted_from_different_sectors_are_comparable`).
`get_spectral_function`/`get_spin_spectral_function` assemble the
particle/hole and Sz channels on top, and a charge-indefinite operator
(`Sx`) is split into `(S+ + S-)/2` and summed channel by channel
(`multioperatortk/charge.py`) rather than rejected. Two things to know
before touching it: the gating is explicit inside the submode (`mode.py`'s
sector guard keys on the caller's `conserved_sector`, which is empty here,
so it never fires), and the per-sector solves are cached with the clone
that produced them (`_sector_states_cache`) because ITensor Index identity
is per session; that cache is also what makes a site sweep, and so
`S(q,w)`/`A(k,w)`, cost about one site's correlator. See
`docs/sector_resolved_spectral_function_plan.md` for the design and
`tests/test_sector_spectral_function.py` for what is pinned.

### STM/Kondo tunneling spectra (`kondospectrumtk/`)

`get_kondo_spectrum` (on `Spin_Chain` and both spinful fermion chains,
`kondospectrumtk/kondomethod.py`) is the third-order perturbation theory of
Ternes (NJP 17 063016, arXiv:1505.04430); `mode="ED"` is the reference
implementation at any T and `mode="DMRG"` reproduces it at T=0 via the
dynamical correlator and a two-time TDVP construction. Every prefactor was
pinned against the paper's Fig. 7, except that the Kondo term's exchange
diagram deliberately departs from the paper's eq. 25, where a brute-force
T-matrix puts the log elsewhere; do NOT "restore" eq. 25 to get back on the
figure. The notes, formerly this file's section of the same name, are in
`.claude/rules/kondo.md`.

`get_esr_spectrum` (`esrstm.py`, `esrtk/`, the same three chains) is the
ESR-STM line shape: the rate equation between eigenstates with Korringa
relaxation and a DC plus homodyne readout of the tip-site spin, on the ED
eigenbasis or on the lowest DMRG eigenstates, never on TASK or a
correlator (excited states are populated at 1 K and 10 to 30 GHz). Plan
and derivation in `docs/esr_stm_plan.md`, architecture in
`docs/development.md` 4.8l, rules in `.claude/rules/kondo.md`.

### Julia vs C++ backend

The two DMRG backends are independent implementations, not a shared
protocol: the C++ path runs in-process through the pybind11 extension,
while `itensor_version="julia_live"` (`mpsjulialive/`) drives a live
Julia/ITensors.jl session via `juliacall` with its own mirrored set of
modules. A feature missing on one side simply is not implemented there yet
(check the relevant `mpsjulialive/*.py` file, or lack thereof) rather than
falling back automatically; the only automatic fallbacks are the DMRG to ED
ones `mode.py` makes (missing extension, too few sites), none of which
substitutes one DMRG backend for another.

### Infinite chains (iDMRG and VUMPS)

`infinitechain.py` is the public API, `pyitensor/idmrg.py`, `vumps.py` and
`vumps_ms.py` the `"python"` implementation, and the `idmrg_*`/`ic_*`/
`vx_*`/`vms_*` helpers in `mpscpp3/chain_session.h` the `itensor_version=3`
one; the upstream reference is ITensor's own
<https://github.com/ITensor/ITensorInfiniteMPS.jl>, not vendored here. The
rules, with the reasoning and the measurements in
`.claude/rules/infinite-chain.md` (formerly this file's notes under the
"ITensor library" header):

- Both environment builders on both backends try the bond candidate first
  (`C C^dag` and `conj(C^dag C)` in this codebase's `X[ket,bra]` ordering,
  `vumps._transfer_fixed_points`, `vumps_ms._cell_fixed_points`,
  `Chain::vx_choose_fixed_point`), a yes/no residual test, and run the
  guarded eigensolver only when it fails; a threshold on `C`'s weight
  spectrum is the WRONG shape, since the ratio to catch is a moving number
  mid-convergence. The LEFT candidate is the transpose, `conj(C^dag C)`,
  and only a measured residual ever caught that.
- VUMPS's Lanczos stops on the Ritz residual (`residual_tol`, passed as
  `tol/10` by `vumps.py`/`vumps_ms.py` and by `Chain::vumps_ground_state`),
  since an eigenvalue criterion caps the eigenvector at ~1e-6 and the
  VUMPS criterion is a difference of two independently solved eigenvectors;
  `AC`/`C` are phase-aligned against their warm starts. Do not "fix" finite
  DMRG's criterion on this reasoning: measured, its sweep contracts the
  per-solve error and the tight criterion moves nothing above 5e-15.
- Any Hamiltonian with reach>1 goes to the sequential solver (`vumps_ms`,
  `Chain::vumps_ground_state`'s `use_multisite`), and the snapshot that
  exists decides the reader; `gs_method="idmrg"`'s `_active_channels_at`
  carries one pending channel per site of a term's reach by design, so do
  NOT add a guard there.
- `idmrg.py`'s transfer chain is `_TransferChain`, deliberately not
  list-like, so a leftover `Es[p]` raises instead of materializing a chi^4
  array; v3's static observables come from the extracted unit cell
  (`idmrg_cell_`), never from tiling the per-micro-step `idmrg_U_`.
- Any `itensor_version=3` VUMPS energy from before 2026-09-12 on a model
  whose exact state is smaller than the requested `maxm` is not comparable
  (`docs/known_issue_v3_vumps_variational_floor.md`, FIXED;
  `tests/test_vumps_redundant_bond_dimension.py`, with a repeated-run guard).

## Documentation

There are three documents. `docs/development.md` is the development
document: the architecture (backends, dispatch, directory layout). The
user guide covers the physics-facing API (what calculation each method
performs, with the formula behind it), kept brief and physics-first, with
no history in it. `docs/numerical_changes.md` is the ledger of what raises
and what changed numbers.

**The user guide is a website.** Its only source is `docs/guide/*.md`, one
page per topic, built with MkDocs Material from `mkdocs.yml` (repo root,
`docs_dir: docs/guide`, so nothing else under `docs/` is published) and
deployed to <https://joselado.github.io/dmrgpy/> by
`.github/workflows/guide.yml` on every push to `master` that touches it.
The PDF is generated from the same pages (`docs/build_guide_pdf.py`, pandoc
and xelatex, also run by the workflow); there is no handwritten `.tex` any
more, and `docs/user_guide.md` is only a pointer. Math is `$...$`/`$$...$$`
(MathJax through `pymdownx.arithmatex`); a `|` inside math in a table row
is written `\vert`/`\lvert`/`\rvert`, or the table splits. Pages link to
each other by page name (`[Ground state](ground-state.md)`), never by a
section number. A new page needs an entry in `nav`, or it is neither in
the menu nor in the PDF. `mkdocs build` is strict and fails on a broken
link or anchor; run it (`pip install mkdocs==1.6.1 mkdocs-material==9.7.7`)
before committing a guide change. Every figure in the guide is computed:
`docs/figures/<id>/compute.py` writes `data.npz` (committed), `plot.py`
reads only that and writes the light and dark PNGs into `docs/guide/img/`;
`docs/figures/README.md` says how to recompute, replot and run on a remote
host, and `docs/figures/plan.md` is the physics each figure checks.
Under each figure a collapsed block shows its `compute.py` and `plot.py`,
copied in at build time by `pymdownx.snippets` with `base_path:
docs/figures`, a path relative to the current directory, so `mkdocs build`
runs from the repo root.

Whenever a new feature is implemented (a new physical model, a new
`Many_Body_Chain` method, a new dynamical-correlator submode, a new
post-processing tool) update the matching page under `docs/guide/`, and
`docs/development.md` too if the change affects architecture or
backend dispatch rather than just physics; the
`docs-sync` skill (`.claude/skills/docs-sync/SKILL.md`) is the procedure.

**The user guide's style is fixed.** Its readers are condensed matter
physicists who want to compute things, not to learn how the code works. Each
capability gets a short physics explanation (the quantity, its formula, what
it is used for, what measures it), then the code, then the few keywords that
control the physics or the accuracy as a one-line-per-item list. Short plain
sentences, one idea each. No audit history, no "used to", no old numbers, no
timings, no implementation internals (caches, penalties, block sparsity,
stopping rules), and a backend caveat only when a number would otherwise be
wrong, in one sentence. Every snippet in it has been run. The
"Ground state" page (`docs/guide/ground-state.md`, formerly section 3,
"Ground-state properties") is the reference: a new or edited page matches
it in length per topic and in tone. History and changed numbers go to
`docs/numerical_changes.md`, never into the guide.

## Git workflow

Working directly on `master` is fine at all times, unless the user says
otherwise for a given task or the change is a really major one with real
potential to conflict with other concurrent development (in which case use
a feature branch/PR instead).

When a feature branch is warranted, new work (a new branch, a new PR) should
start from `master`, not from another feature/staging branch, unless
explicitly told otherwise. This repo previously reused a long-lived branch
(`mpo-algebra-tdvp-gse`) as a staging area for several stacked PRs; once
`master` and the staging branch diverged (each gained commits the other
lacked), a PR built on the staging branch picked up unrelated content and
produced conflicting/misleading diffs against `master`. Branching from
`master` (`git fetch origin master && git checkout -b <branch>
origin/master`) avoids this.
