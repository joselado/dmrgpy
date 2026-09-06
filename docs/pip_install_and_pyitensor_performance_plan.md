# Making `pip install dmrgpy` usable: the ED fallback, and `to_mpo`'s O(L^4)

**Status: both fixed.** See "What was done" at the end of this file for
the implementation, the measurements and the tests. Everything in the two
problem statements below is kept in the present tense as it was written,
because the diagnoses are what the fixes were built on.

Two independent problems, both reported from outside this repo (a session
writing the README of a course whose notebooks build `Spin_Chain`/
`Fermionic_Chain` with the default backend and call `gs_energy(mode="DMRG")`
at 40-100+ sites, and which wants `pip install dmrgpy` to be the recommended
path). Everything below was measured in this checkout at
commit `a2eb46e`, not taken on report. The reporter's original absolute
timings turned out to be wrong and have been retracted by them; the
scaling problem they pointed at is real, and is the subject of problem 2.
See "Settled: no regression, and where the reported numbers came from".

The two problems compound: a pip install silently gives you ED (problem 1),
and the backend it *should* give you instead is too slow at the sizes those
notebooks use (problem 2). Fixing only the first would route users onto a
backend that cannot do the job.

---

## Problem 1: a pip install silently falls back to ED, not to the pure-Python MPS backend

Confirmed. `mode.py`'s `get_mode` returns `"ED"` whenever the requested
C++ extension is missing:

```python
# src/dmrgpy/mode.py:62-64
if not cppext.available(self.itensor_version):
    print("C++ extension not compiled, using default ED routines")
    return "ED" # use exact diagonalization
```

and `cppext.py:32-33` sets `DEFAULT_ITENSOR_VERSION = 3` with the comment
that the `"python"` backend is "never selected by default". The wheel ships
no C++ (see CLAUDE.md's "Packaging / PyPI"), so **every default-backend chain
in a pip install runs exact diagonalization**, which cannot reach the sizes
the pyitensor backend exists to serve. The package does warn at import, but
a warning is not a working default: scripts written against the compiled
backend keep running and silently change algorithm.

### What a fix has to respect

The obvious one-line change (make `get_mode` return `"DMRG"` and let it
route to pyitensor) is not enough, because `get_mode` returns a *mode*
(`"ED"`/`"DMRG"`) while the backend choice lives in a different attribute,
`self.itensor_version`. The flip therefore belongs where the version is
chosen -- `Many_Body_Chain.__init__`, or `setup_cpp()` -- switching
`itensor_version` from `3` to `"python"` when `cppext.available(3)` is
False, and downgrading the import-time ED warning to a one-line note.

Three things in `mode.py` must NOT be swept into the same change:

* **The `ns < 3` fallback.** `itensor_version==3` with fewer than 3 sites
  aborts the whole process inside ITensor v3 (`"LocalOp is default
  constructed"`), so `get_mode` routes it to ED. That is a real mpscpp3 bug
  and is v3-specific: the pure-Python backend has no such limit, so a chain
  that is *already* on `"python"` should stay on DMRG here rather than
  inheriting v3's workaround.
* **The conserved-sector guard**, which raises rather than falling back --
  answering with the global ground state instead of the requested sector
  would be silently wrong.
* **ED as a deliberate choice.** `mode="ED"` and small-system
  cross-checking are the point of that backend; only the *implicit*
  fallback should change.

### Scope note

This is a behaviour change to the default backend of every chain in the
package, so it wants its own commit and a test that a chain built with no
C++ available reports `"python"` rather than ED. `tests/` has no coverage of
the `itensor_version="python"` path at all today (CLAUDE.md records this:
`_helpers.py`'s `versions=` only covers `(2, 3)`), which is worth fixing in
the same change since the wheel's *primary* DMRG backend would then be the
only untested one.

---

## Problem 2: `to_mpo` is O(L^4) and is 95% of `gs_energy` at L=100

Measured here, Heisenberg S=1/2 open chain, `Sx Sx + Sy Sy + Sz Sz`
nearest-neighbour, defaults (`maxm=30`, `nsweeps=15`),
`itensor_version="python"`, one `gs_energy(mode="DMRG")` call,
`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1`:

| L | wall time | ratio vs previous |
|---|---|---|
| 20 | 1.08 s | |
| 40 | 8.29 s | 7.7x for 2x L |
| 60 | 33.0 s | 4.0x for 1.5x L |
| 100 | 239 s | 7.2x for 1.67x L |

From L=40 to L=100 (2.5x) the time grows 28.8x, i.e. **~L^3.7**, where fixed
`maxm` and fixed `nsweeps` should give ~L.

`cProfile` on the same call says exactly where it goes:

| | L=40 | L=100 |
|---|---|---|
| `gs_energy` total | 9.73 s | 245.5 s |
| `mpobuilder.to_mpo` (2 calls) | 6.19 s (64%) | **232.7 s (95%)** |
| — of which `svd._svd_truncated` | 4.94 s | 186.8 s |
| `manybodychain.is_hermitian` | 5.77 s (59%) | 217.7 s (89%) |
| the actual DMRG (`pyitensor.dmrg.dmrg`) | 3.35 s | ~13 s |

**The DMRG is fine.** 3.35 s -> ~13 s for 2.5x the sites is essentially
linear, which is what the algorithm should do. The entire scaling problem
is MPO construction, and it is worth stating plainly because it is
counter-intuitive: at L=100 this backend spends 95% of a "ground state
calculation" building the Hamiltonian and 5% solving it.

### Root cause

`pyitensor/mpobuilder.py`'s `to_mpo` does not port ITensor's automaton
MPO-compression algorithm (a deliberate simplification, documented in its
own module docstring and in `pyitensor/__init__.py`). Instead every `HTerm`
becomes its own exact bond-dimension-1 MPO, all `T` of them are
block-diagonally concatenated by `mpsalgebra.sum_many()`, and a single
bidirectional truncating sweep compresses the result.

The *final* bond dimension is correct -- the docstring records 39 -> 5 on a
nearest-neighbour Heisenberg chain at N=14, and that is the constant 5 it
should be. **The cost is the intermediate.** The concatenated MPO carries
bond dimension ~`T` = 3(L-1) before compression, so the compression sweep
runs `O(L)` truncating SVDs on matrices of size `O(L)`, i.e. `O(L * L^3) =
O(L^4)`. That matches the measured L^3.7 and matches `_svd_truncated` being
76% of total wall time at L=100.

Note this refines, and partly corrects, the first guess made about it: the
problem is not that the MPO bond dimension grows with L (it does not), but
that the *route to* the correct constant-5 MPO goes through an object whose
bond dimension does.

### A second, independent multiplier: the MPO is built twice, and cached never

`pyitensor/chain.py`'s `_mpo()` calls `to_mpo` fresh on every invocation --
there is no cache:

```python
# src/dmrgpy/pyitensor/chain.py:393
def _mpo(self, terms, cutoff=_BUILD_CUTOFF, maxdim=None):
    return to_mpo(self._ampo(terms), cutoff=cutoff, ...)
```

so a single `gs_energy` builds the Hamiltonian MPO **twice** (`ncalls = 2`
in both profiles): once for `groundstate.gs_energy`'s unconditional
`self.is_hermitian(self.hamiltonian)` check, once for the DMRG itself. That
is why `is_hermitian` costs 89% of `gs_energy` at L=100 -- not because the
check is expensive in itself, but because it rebuilds the MPO.

`mpsalgebra.is_hermitian` has already been optimized once for exactly this
symptom (its docstring records cutting it from 31-38% of `gs_energy` by
building its witness state at a small fixed bond dimension instead of
`self.maxm`). That fix capped the *witness*; it did not touch the *MPO
build*, which is what now dominates.

### Fixes, cheapest first

1. **Cache the MPO on the `Chain`**, keyed on the terms (and the sector,
   which `_ampo` already folds in). Roughly a free 2x on every `gs_energy`,
   and more for any workflow that calls several observables on one chain.
   Small, local, testable; does not change any number.
2. **Do not rebuild an MPO for the Hermiticity check.** `is_hermitian`
   needs to know whether `op - op.get_dagger()` is exactly zero. With (1)
   in place this mostly falls out; without it, the check could compare
   terms symbolically at the `MultiOperator` level and never build an MPO
   at all, which would be better still.
3. **Port ITensor's automaton MPO compiler**, which constructs the
   bond-dimension-5 MPO directly instead of reaching it by compressing a
   bond-dimension-3(L-1) one. This is the actual fix for the scaling and
   the only one that gets L=100 to seconds. It is real work -- it is the
   piece `pyitensor` explicitly chose not to port -- but it is also
   well-specified, self-contained, and has a reference implementation in
   `mpscpp3/ITensor`. Intermediate option if that is too much: compress
   incrementally (fold in terms in batches, compressing as you go) so the
   intermediate bond dimension stays O(1) rather than O(L).

(1) and (2) are a constant-factor win and leave the wall in place; only (3)
changes the exponent. For the course's L=100 target, (3) is required.

### What was ruled out

* Not the MPO's final bond dimension -- it compresses to 5 correctly.
* Not `set_hamiltonian`, which is ~0.00 s at every L measured (it stores a
  `MultiOperator`; the MPO is built later, inside `gs_energy`).
* Not BLAS thread oversubscription -- but do not read that as "threads
  don't matter here", because on the reporting machine they were worth
  more than everything else combined. On this (otherwise idle) host,
  unpinned vs pinned was 2.04 s vs 1.08 s at L=20 and 12.7 s vs 8.3 s at
  L=40, i.e. ~1.5-2x. On theirs, unpinned gave 25.6 s at L=20 against
  1.15 s pinned -- a factor of ~22 -- and L=40 did not finish in ten
  minutes against ~9 s pinned. That is the regime CLAUDE.md's benchmarking
  note describes, where oversubscription dominates every other effect in
  the file. **Pin threads before timing anything on this backend**;
  every number above is pinned. See `src/dmrgpy/blasthreads.py`.
* Not the DMRG sweeps, per the table above.

---

## Settled: no regression, and where the reported numbers came from

The first version of this note flagged a 20x discrepancy between the
reporter's timings (L=12: 4.0 s, L=20: 25.6 s, L=40: did not finish in ten
minutes) and this checkout's (0.28 s / 1.08 s / 8.3 s), and asked whoever
picked it up to re-measure against the installed wheel before doing
anything -- because a since-fixed regression would have made the right
action a release rather than a port.

**That has been done, and there is no regression.** The PyPI 0.1.1 wheel and
`src@c78909e` are identical, on the reporting machine, with threads pinned:

| | L=12 | L=20 | L=40 |
|---|---|---|---|
| wheel 0.1.1 | 0.29 s | 1.15 s | 8.77 s |
| src `c78909e` | 0.30 s | 1.16 s | 8.60 s |

with `E/L` agreeing to every printed digit, and both matching the numbers
measured here. The original timings were BLAS thread oversubscription on an
unpinned run (and a 10-minute command cap killing L=40), and have been
retracted.

So: **the wheel is not behind, a release fixes nothing here, and the
priority of fix (3) stands as written.** The scaling exponent (~L^3.7) was
always a property of the algorithm and is unaffected either way.

The practical consequence for the course that prompted this: L=40 in ~9 s
is fine for exercises, and L=100 in minutes is workable for self-study but
borderline for a live session -- which is exactly the range fix (3) would
move.

## Reproduction

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 - <<'EOF'
import time
from dmrgpy import spinchain
for L in (20, 40, 60, 100):
    sc = spinchain.Spin_Chain(["S=1/2"] * L, itensor_version="python")
    H = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(L - 1))
    sc.set_hamiltonian(H)
    t = time.time(); e = sc.gs_energy(mode="DMRG")
    print(L, e / L, time.time() - t)
EOF
```

`e/L` converges toward -0.4431 (Bethe ansatz, with the open-boundary
finite-size correction); the values here run -0.4341 (L=20) to -0.4413
(L=100), which is the expected approach and confirms the calculation is
correct at every size -- this is a performance problem only. Swap
`cProfile` around the `gs_energy` call, sorted by cumulative time, to see
the `to_mpo` split.

For problem 1, run anything with the default backend in an environment
where `cppext.available(3)` is False and observe "C++ extension not
compiled, using default ED routines".


---

## What was done

Both problems are fixed, in that order reversed: the exponent first, since
flipping the default (problem 1) before fixing the scaling (problem 2)
would have routed every pip user onto a backend that takes minutes at
L=100.

### Problem 2 — `to_mpo` is now a finite-state machine, not a compression

Fix (3) from the list above, taken directly rather than via the
intermediate options. `pyitensor/mpobuilder.py::to_mpo` assembles the MPO
as a finite-state machine over the terms' partial products -- the shape
ITensor's own `toMPO(...,{"Exact",true})` produces -- instead of
concatenating T bond-dimension-1 MPOs and compressing. The states at a
bond are `I` (term not started), `F` (finished) and one per *distinct*
left-partial product among the terms straddling that bond; sharing partial
states between terms is what compresses, and it creates no spurious paths
because sharing a state means the prefixes are identical. Two rules are
load-bearing: the coefficient goes on the transition *into* `F` (so terms
differing only by a coefficient still share partial states), and
transitions into `F` accumulate while structural transitions are assigned
(identical terms trace the same path -- coefficients must sum, structure
must not be doubled).

The build is now O(L^2) rather than O(L^4). The machine itself is O(L);
what is left of the square is `HTerm.resolve()` spelling each of the O(L)
terms out over all L sites -- cheap (small dense per-site matrices, ~30k
of them at L=100), nowhere near dominant at the sizes measured, but it is
the next wall if one ever appears.

The two bidirectional truncating sweeps are **kept**, now purely to honour
the caller's `cutoff`/`maxdim` and to squeeze out redundancy prefix-sharing
cannot see. They are cheap because the bond dimension going into them is
now O(1) instead of O(T).

`to_mpo` in isolation, nearest-neighbour Heisenberg, min of 5, threads
pinned (note the machine was under load from an unrelated job throughout,
so treat the ratios as the result and the absolutes as an upper bound):

| L | old | new | speedup |
|---|---|---|---|
| 20 | 0.125 s | 0.027 s | 4.6x |
| 40 | 0.663 s | 0.029 s | 22.7x |
| 60 | 2.96 s | 0.119 s | 24.8x |
| 100 | 24.6 s | 0.235 s | **105x** |

End to end, `gs_energy(mode="DMRG")` at L=100 went from 239 s to ~42 s on
this host; what remains is the DMRG itself, which was always linear and
was never the problem. Below L≈40 the end-to-end win is small, because
there the MPO build was never the dominant cost.

Fixes (1) and (2) from the list above -- caching the MPO on the `Chain`,
and not rebuilding one for the Hermiticity check -- were **not** done.
They were a constant-factor 2x on top of an O(L^4) build; against an O(L)
build they buy a fraction of a second at L=100 and are not worth the cache
invalidation surface. Worth revisiting only if a profile says otherwise.

The previous construction is kept as `_sum_of_term_mpos` and is now the
independent reference the machine is checked against, over a zoo of term
shapes (long-range, gaps in the support, several factors per site, complex
coefficients, mixed local dimensions, bare-coefficient terms, spinless and
spinful fermions with their Jordan-Wigner strings, odd fermion parity,
n=1 and n=2 chains) in `tests/test_mpo_automaton_builder.py`. The real
reference there is `AutoMPO.dense_matrix()`, which Kronecker-multiplies
the same per-site matrices without ever forming an MPO and so shares no
code with either builder.

**A bug found on the way, worth recording.** On a **one-site** chain the
old builder silently kept only the *last* term: `to_mpo` had no bonds to
sweep there, so `sum_many`'s concatenation was returned as-is, and
`0.8*Sz + 0.6*Sx` came back as `0.6*Sx` alone. Through the public API that
made a 1-site Hamiltonian a different operator -- and, being non-Hermitian
by accident, sent it to NH-DMRG. Regression-tested.

### Problem 1 — the default backend, not the mode

Implemented where the note said it belonged: at the point the *version* is
chosen, not in `get_mode`. `Many_Body_Chain.__init__` (and
`Mixed_Spin_Fermion_Chain.__init__`) now take `itensor_version=None`
meaning "pick one for me", resolved through the new
`cppext.default_backend()` -- `DEFAULT_ITENSOR_VERSION` when that
extension is compiled, `"python"` when it is not. The import-time warning
in `dmrgpy/__init__.py` no longer says chains fall back to ED, because
they no longer do; it now says they default to the pure-Python backend and
that compiling the C++ one is a speed choice.

All three things the note said must not be swept into this were kept:

* **The `ns<3` fallback** stays scoped to `itensor_version==3`. A chain
  resolved to `"python"` skips it, as intended.
* **The conserved-sector guard** still raises rather than falling back.
* **`mode="ED"`** is untouched, and so is an *explicit* `itensor_version=3`
  on a machine with no extension -- that still falls back to ED, because a
  caller who named a version asked for that backend specifically.

`tests/test_default_backend_without_cpp.py` pins all of it, simulating the
missing extension by emptying `cppext._backends` so the real decision is
exercised on a machine where the extension is built. That file is also the
first coverage `tests/` has had of the `itensor_version="python"` dispatch
path at all, which the note flagged as worth fixing in the same change
given the wheel's primary DMRG backend would otherwise be the only
untested one.

### A third, smaller fix that fell out of the first

Making the 1-site MPO correct exposed that `itensor_version="python"`
cannot do DMRG on a 1-site chain *at all*: pyitensor's DMRG is two-site
(`dmrg.py::_dmrg_one_sweep` sweeps `for i in range(1, n)`), so the sweep
body never runs and `dmrg()` returns the `energy = None` it started with.
Before the MPO fix this was masked -- the wrong 1-site operator happened to
be non-Hermitian, so the call went to NH-DMRG and returned a wrong number
instead of `None`.

`mode.py::resolve_mode` now routes `itensor_version=="python"` with
`ns < 2` to ED, the same mechanism as the existing `itensor_version==3`
with `ns < 3` fallback next to it. The threshold is 2, not 3: a two-site
chain has exactly one two-site update and solves correctly, checked
against ED.

`gs_energy_generalized` needed its own copy of the guard, in
`groundstate.py`, because it has no ED fallback to be routed to -- and
there the symptom was worse than `None`: the outer self-consistent
iteration still returns a lambda, the Rayleigh quotient of a state no
sweep ever touched, so a 1-site chain answered **-0.3049 for an exact
-0.5**. Unlike the `itensor_version==3` guard a few lines below it, this
one sits *before* the non-Hermitian dispatch: NH-DMRG escapes ITensor
v3's short-chain abort because it never calls `dmrg()`, but its own sweep
is two-site as well, so it does not escape this one (checked -- same
wrong-number symptom).

None of these three 1-site bugs was reachable from the default backend
before, on a machine with a compiled extension: `itensor_version=3` with
`ns<3` was already routed to ED. They mattered because problem 1's fix
makes `"python"` the default wherever there is no extension, which is
exactly where a 1-site chain would now meet them.
