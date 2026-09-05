# Making `pip install dmrgpy` usable: the ED fallback, and `to_mpo`'s O(L^4)

Two independent problems, both reported from outside this repo (a session
writing the README of a course whose notebooks build `Spin_Chain`/
`Fermionic_Chain` with the default backend and call `gs_energy(mode="DMRG")`
at 40-100+ sites, and which wants `pip install dmrgpy` to be the recommended
path). Neither is fixed. Everything below was measured in this checkout at
commit `a2eb46e`, not taken on report -- and one of the two headline claims
did not reproduce, so read the "What did not reproduce" section before
trusting the reporter's numbers.

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
* Not BLAS thread oversubscription, though that is a real secondary effect:
  unpinned vs pinned was 2.04 s vs 1.08 s at L=20 and 12.7 s vs 8.3 s at
  L=40, i.e. ~1.5-2x, no more. Timings above are all pinned. See
  `src/dmrgpy/blasthreads.py`.
* Not the DMRG sweeps, per the table above.

---

## What did not reproduce

The reporter measured L=12: 4.0 s, L=20: 25.6 s, L=40: did not finish in 10
minutes. Here the same script gives L=12: 0.28 s, L=20: 1.08 s, L=40: 8.3 s
-- roughly 20x faster at L=20 and finishing at L=40. Thread pinning explains
at most 2x of that.

The untested differences are that they ran **PyPI dmrgpy 0.1.1 in a clean
venv** while this was run from `src/` at `a2eb46e`, and that it was a
different machine. Since 0.1.1 was tagged at `7cbe80e`, a fix landing
between `7cbe80e` and now would show up exactly like this. **Anyone picking
this up should first re-measure against the installed wheel**, because if
the gap is a since-fixed regression then the priority of (3) above drops
sharply, and the real action is a release rather than a port. The scaling
exponent measured here (~L^3.7) is a property of the algorithm and would not
change either way.

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
