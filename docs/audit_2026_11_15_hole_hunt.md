# Audit, 2026-11-15: four findings from the backend benchmark

These four defects were not found by a hunt but by the cross-backend
benchmark of 2026-11-15 (`benchmarks/run_benchmarks.py --adaptive`, run as a cluster array job
on one Intel Xeon Gold 6248 core per task, `maxm=128`,
`nsweeps=10`, `cutoff=1e-10`). That run compares v3, `"python"` and
`"julia_live"` on the same chains, and checks every answer against exact
diagonalization where it is cheap and against the other backends above
that. Each defect was then reproduced locally on `e55a400` plus
uncommitted work that does not reach any of the code paths below (the
benchmark harness, the right-environment reuse in `pyitensor/dmrg.py`'s
`dmrg_excited`, and a new `method="state_averaged"` route behind a keyword
in `excited.py`). Every repro below was executed and its output is spliced
from disk; RuntimeWarning lines are removed from the outputs and nothing
else is. Unlike a hunt, no finding has been handed to an independent
reviewer yet: they were recorded on the author's request, to be addressed
in a later session, so each one is `OPEN` and its cause is marked as
measured or as a reading where it is one.

This file is the evidence, not a task list. Each finding keeps its repro
and gains a `**Status**` line when it is acted on.

## Scope

The four code paths below, on the backends named. Out of scope by
construction, as in every record: the vendored ITensor and TDVP folders,
the legacy quirks of `.claude/rules/cpp-backend.md`, the open
`docs/known_issue_*.md` items and what earlier records hold. Every repro
runs as

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 <scratch>/<script>
```

with `<scratch>/common.py` putting `<repo>/src` and `<repo>/benchmarks`
first on `sys.path` and refusing to run on any other checkout. The chains
are the benchmark's own builders in `benchmarks/calctk.py`: `make_chain`
is `H = 2 sum_i S_i.S_{i+1}` on S=1/2, `make_nh_chain` the PT-symmetric XX
chain `sum_i (Sx_i Sx_{i+1} + Sy_i Sy_{i+1}) + i 0.3 sum_i (-1)^i Sz_i`,
`make_spinless_chain` spinless fermions with hopping 1 and
`(N_i-1/2)(N_{i+1}-1/2)`.

```python
import os, sys
REPO = "<repo>"
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "benchmarks"))
import numpy as np  # noqa: F401  (before dmrgpy)
import dmrgpy
assert os.path.realpath(dmrgpy.__file__).startswith(os.path.join(REPO, "src")), dmrgpy.__file__
HEAVY = dict(maxm=128, nsweeps=10, cutoff=1e-10)
```

## Findings

### 1. `"python"` returns a positive ground-state energy, +36.42 for about -70.5, in the Sz=0 sector of the 80-site Heisenberg chain

`bug` &middot; severity **HIGH** &middot; OPEN &middot; found by the benchmark

**Where**: `src/dmrgpy/pyitensor/sector.py` (start state `sector_mps`,
charge penalty), `src/dmrgpy/pyitensor/chain.py::_default_mps`, the
sector branch of `gs_energy`.

With `set_conserved_sector(Sz=0)`, `maxm=128`, `nsweeps=10` and
`cutoff=1e-10`, the 56-site chain gives -49.26125979, equal to v3 and to
the dense solve. The 80-site chain gives +36.42182885 locally, with no
warning and no exception; the v3 sector solve gives -70.530474 on the same
chain. On the cluster run the same call raised instead, "the converged state is not
in the requested sector (<(Q-q)^2> = 1.23e-06)". So whether the state
leaves the sector and is caught, or comes back with a positive energy, is
decided by the draw; a positive energy for a Heisenberg antiferromagnet is
not a ground state in either sector. Not yet known: whether the penalty
energy is part of the reported number (its sign would fit), at what size
between 56 and 80 it starts, and whether a smaller `maxm` moves it.

**Expected**: the v3 value, -70.530474, as at every size up to 56.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 <scratch>/r3_python_sector.py
```
```python
# python Sz=0 sector ground state on the benchmark chain, heavy settings
from common import np, HEAVY
import time, calctk
for n in (56, 80):
    np.random.seed(1234)
    sc = calctk.make_chain(n, "python", HEAVY)
    sc.set_conserved_sector(Sz=0)
    t0 = time.perf_counter()
    try:
        e = sc.gs_energy(mode="DMRG")
        print(f"N={n:3d} gs_energy={e:.8f}  ({time.perf_counter()-t0:.1f} s)")
    except Exception as ex:
        print(f"N={n:3d} raised after {time.perf_counter()-t0:.1f} s: {ex!r}")
```

Observed:

```
N= 56 gs_energy=-49.26125979  (41.3 s)
N= 80 gs_energy=36.42182885  (48.8 s)
```

The cluster benchmark (its `sector` task) gave, for the same
call, `-49.26126` on all three backends at N=56 and, at N=80, v3
`-70.530474` and `"python"`: `RuntimeError('gs_energy: the converged state
is not in the requested sector (<(Q-q)^2> = 1.23e-06). Raise the charge
penalty with set_sector_penalty(), or the sweep count/bond dimension.')`.

**Suggested fix**: first find whether the +36.42 is H plus the penalty, by
evaluating `<H>` and `<(Q-q)^2>` on the returned state; then whether the
start state (a sum of random in-sector product states at `maxm`) is
already the problem at this length, by solving from a single in-sector
product state.

### 2. `"julia_live"` converges to a wrong state in the Sz=0 sector from N=56 on, -48.09 for -49.26, and to positive energies from N=224

`bug` &middot; severity **HIGH** &middot; OPEN &middot; found by the benchmark

**Where**: `src/dmrgpy/mpsjulialive/sector.py`, `sites_qn.jl` and the
QN-carrying DMRG call the sector route makes.

The returned state is in the sector (`<Sz_tot>` 2e-15) and its `<H>`
agrees with the reported energy, so the number is not misreported: the
solve converges to a state that is not the ground state. At N=40 all three
backends agree with v3 to 1e-8. In the benchmark the julia_live energy then
drifts away from v3 with N, and is positive from N=224 (v3 `-198.153978`,
julia_live `16.177386`; at N=640 `-566.84188` against `222.190641`), at a
fraction of v3's time, so the fast julia_live sector times above N=40 in
the benchmark are not the cost of a ground state.

**Expected**: v3's -49.26125980 at N=56.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 <scratch>/r1_julia_sector.py
```
```python
# julia_live Sz=0 sector ground state against v3 on the benchmark chain
from common import np, HEAVY
import calctk
for n in (40, 56):
    for iv in (3, "julia_live"):
        np.random.seed(1234)
        sc = calctk.make_chain(n, iv, HEAVY)
        sc.set_conserved_sector(Sz=0)
        e = sc.gs_energy(mode="DMRG")
        wf = sc.get_gs()
        h = sc.hamiltonian
        eh = wf.dot(h * wf).real if hasattr(wf, "dot") else float("nan")
        sz = sum(sc.Sz)
        print(f"N={n:3d} {str(iv):10s} gs_energy={e:.8f}  <H>={eh:.8f}  <Sz_tot>={sc.vev(sz).real:+.2e}")
```

Observed:

```
N= 40 3          gs_energy=-35.08294658  <H>=-35.08294652  <Sz_tot>=+5.55e-17
N= 40 julia_live gs_energy=-35.08294658  <H>=-35.08294650  <Sz_tot>=+0.00e+00
N= 56 3          gs_energy=-49.26125980  <H>=-49.26125964  <Sz_tot>=+0.00e+00
N= 56 julia_live gs_energy=-48.09088891  <H>=-48.09088881  <Sz_tot>=+2.41e-15
```

**Suggested fix**: compare the julia_live start state of the sector route
with v3's (a random QN MPS at link dimension min(maxm,10) against the
product-state sum of the other backends), and run the same QN-carrying
ITensors.jl solve by hand on 56 sites with more sweeps and with noise, to
tell a start-state problem from a solver setting.

### 3. v3 and `"python"` `nhdmrg()` settle on a fixed point that is not an eigenpair when four eigenvalues share the smallest real part, already at N=10

`bug` &middot; severity **MEDIUM** &middot; OPEN &middot; found by the benchmark

**Where**: `src/dmrgpy/pyitensor/nhdmrg.py` (`_arnoldi_smallest_real`,
`_select_ritz`, the fidelity truncation in `_nhdmrg_one_sweep`),
`src/dmrgpy/mpscpp3/chain_session.h` (`Chain::nhdmrg`,
`arnoldi_smallest_real`), the retry loop of `src/dmrgpy/nhdmrg.py`.

On the 10-site PT-symmetric XX chain the four lowest eigenvalues share
`Re E = -2.566782` (a real eigenvalue twice and the pair
`-2.566782 +- 0.264096i`). Both backends warn ("nhdmrg did not reach the
residual tolerance after 5 tries"), so this is not silent, but the value
they return is not an eigenvalue: relative residuals 0.29 to 0.44, and
`Re E` up to 0.032 below the lowest real part, which a user reading only
the energy takes for a ground state. `"julia_live"`, which runs
ITensorNHDMRG.jl itself, returns -2.566782271488228 on the same chain
(cluster benchmark, its `nh_gs` task).

What was measured about the cause (`x2`, `x3`, `x4` below): a larger
Krylov space helps and is not enough (`krylovdim=40`: residual 1.6e-2);
more restarts, more sweeps and the `exp(i theta) H` rotation that the
julia_live driver uses for real-part ties do not help; replacing the local
Arnoldi of the python backend by a fully converged ARPACK solve
(`which="SR"`, six roots, `tol=1e-12`) brings the residual from 0.33 to
1.3e-2, and then 10, 30 and 60 sweeps return the same pair to eight digits.
So there are two parts. The local Arnoldi restarts on a single Ritz vector
and loses the cluster; and with a converged local solve the sweep still
has a fixed point that is not an eigenpair, which points at the step where
the ports differ from ITensorNHDMRG.jl on purpose, the shared "fidelity"
isometry of the left and right states (a reading, not yet measured).

**Expected**: an eigenvalue of smallest real part, as julia_live returns.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 <scratch>/r2_nh.py
```
```python
# nhdmrg() on the PT-symmetric XX chain, i 0.3 (-1)^i Sz_i, against full ED
from common import np, HEAVY
import calctk
n = 10
ref = calctk.make_nh_chain(n, "python", HEAVY)
Hd = np.array(ref.get_ED_obj().get_hamiltonian().todense())
ev = np.linalg.eigvals(Hd)
print("full ED, lowest real parts:", np.round(np.sort_complex(ev[np.argsort(ev.real)][:4]), 6))
for params in (HEAVY, dict(maxm=30, nsweeps=10, cutoff=1e-10)):
    for iv in ("python", 3):
        np.random.seed(1234)
        sc = calctk.make_nh_chain(n, iv, params)
        e, psil, psir = sc.nhdmrg()
        r = sc.hamiltonian * psir - e * psir
        res = abs(r.dot(r)) ** 0.5 / abs(psir.dot(psir)) ** 0.5
        print(f"maxm={params['maxm']:3d} {str(iv):6s} E={complex(e):.8f}  ||H r - E r||/||r||={res:.2e}"
              f"  dist to spectrum={np.min(abs(ev - e)):.2e}")
```

Observed:

```
full ED, lowest real parts: [-2.566782+0.j       -2.566782+0.264096j -2.566782-0.j
 -2.566782-0.264096j]
maxm=128 python E=-2.58848584+0.04292640j  ||H r - E r||/||r||=3.08e-01  dist to spectrum=4.81e-02
maxm=128 3      E=-2.53892971-0.11404201j  ||H r - E r||/||r||=4.36e-01  dist to spectrum=1.17e-01
maxm= 30 python E=-2.59077859-0.04949930j  ||H r - E r||/||r||=3.26e-01  dist to spectrum=5.50e-02
maxm= 30 3      E=-2.59583469+0.04370225j  ||H r - E r||/||r||=2.88e-01  dist to spectrum=5.25e-02
```

What cures it, `x2_nh_knobs.py` (output as printed):

```python
# What cures nhdmrg() on the 10-site PT chain (four eigenvalues share Re E)?
from common import np
import warnings, time, calctk
warnings.simplefilter("ignore")
n = 10
P = dict(maxm=30, nsweeps=10, cutoff=1e-10)
ref = calctk.make_nh_chain(n, "python", P)
ev = np.linalg.eigvals(np.array(ref.get_ED_obj().get_hamiltonian().todense()))
emin = ev.real.min()
def run(iv, theta=0.0, nsweeps=10, **kw):
    np.random.seed(1234)
    sc = calctk.make_nh_chain(n, iv, dict(P, nsweeps=nsweeps))
    h = sc.hamiltonian
    if theta:
        sc.set_hamiltonian(np.exp(1j * theta) * h)
    t0 = time.perf_counter()
    try:
        e, psil, psir = sc.nhdmrg(**kw)
    except Exception as ex:
        return f"raised {ex!r}"[:90]
    e = e / np.exp(1j * theta)
    r = h * psir - e * psir
    res = abs(r.dot(r)) ** 0.5 / abs(psir.dot(psir)) ** 0.5
    return (f"E={complex(e):.6f} Re-Emin={e.real-emin:+.1e} res={res:.1e} "
            f"dist={np.min(abs(ev-e)):.1e} ({time.perf_counter()-t0:.0f} s)")
for iv in ("python", 3):
    print(iv, "default            ", run(iv))
    print(iv, "krylovdim=40       ", run(iv, krylovdim=40))
    print(iv, "restarts=6         ", run(iv, restarts=6))
    print(iv, "nsweeps=30         ", run(iv, nsweeps=30))
    for th in (1e-2, 1e-3):
        print(iv, f"exp(i {th:g}) H     ", run(iv, theta=th))
```
```
python default             E=-2.590779-0.049499j Re-Emin=-2.4e-02 res=3.3e-01 dist=5.5e-02 (24 s)
python krylovdim=40        E=-2.566919-0.264041j Re-Emin=-1.4e-04 res=1.6e-02 dist=1.5e-04 (42 s)
python restarts=6          E=-2.583256-0.056254j Re-Emin=-1.6e-02 res=2.8e-01 dist=5.9e-02 (50 s)
python nsweeps=30          E=-2.565193+0.260702j Re-Emin=+1.6e-03 res=1.4e-01 dist=3.7e-03 (70 s)
python exp(i 0.01) H      E=-2.616340+0.270307j Re-Emin=-5.0e-02 res=5.3e-01 dist=5.0e-02 (24 s)
python exp(i 0.001) H      E=-2.590685-0.046608j Re-Emin=-2.4e-02 res=3.0e-01 dist=5.2e-02 (23 s)
3 default             E=-2.599250+0.020800j Re-Emin=-3.2e-02 res=3.0e-01 dist=3.9e-02 (49 s)
3 krylovdim=40        E=-2.566919+0.264041j Re-Emin=-1.4e-04 res=1.6e-02 dist=1.5e-04 (73 s)
3 restarts=6          E=-2.592548-0.037162j Re-Emin=-2.6e-02 res=2.2e-01 dist=4.5e-02 (63 s)
3 nsweeps=30          E=-2.616332+0.270314j Re-Emin=-5.0e-02 res=5.3e-01 dist=5.0e-02 (94 s)
3 exp(i 0.01) H      E=-2.616338+0.270304j Re-Emin=-5.0e-02 res=5.3e-01 dist=5.0e-02 (32 s)
3 exp(i 0.001) H      E=-2.593846+0.045073j Re-Emin=-2.7e-02 res=3.1e-01 dist=5.3e-02 (31 s)
```

A converged local solver, `x3_nh_arpack.py` (output as printed):

```python
# Does a converged local eigensolver cure python nhdmrg() on the PT chain?
from common import np
import warnings, time, calctk
import scipy.sparse.linalg as sla
import dmrgpy.pyitensor.nhdmrg as pnh
warnings.simplefilter("ignore")
orig = pnh._arnoldi_smallest_real
calls = dict(n=0)
def arpack(matvec, x0, krylovdim, restarts, sel="SR", target=0.0):
    calls["n"] += 1
    x0 = np.asarray(x0, dtype=complex)
    d = x0.size
    if d <= 8:
        M = np.column_stack([matvec(c) for c in np.eye(d, dtype=complex)])
        ev, V = np.linalg.eig(M)
    else:
        A = sla.LinearOperator((d, d), matvec=matvec, dtype=complex)
        k = min(6, d - 2)
        ev, V = sla.eigs(A, k=k, which="SR", v0=x0, ncv=min(d, max(2 * k + 1, krylovdim)),
                         tol=1e-12, maxiter=2000)
    i = pnh._select_ritz(ev, sel, target)
    return complex(ev[i]), V[:, i]
n = 10
for maxm in (30, 128):
    P = dict(maxm=maxm, nsweeps=10, cutoff=1e-10)
    ref = calctk.make_nh_chain(n, "python", P)
    ev = np.linalg.eigvals(np.array(ref.get_ED_obj().get_hamiltonian().todense()))
    emin = ev.real.min()
    for name, solver in (("arnoldi", orig), ("arpack ", arpack)):
        pnh._arnoldi_smallest_real = solver
        np.random.seed(1234)
        sc = calctk.make_nh_chain(n, "python", P)
        t0 = time.perf_counter()
        e, psil, psir = sc.nhdmrg()
        r = sc.hamiltonian * psir - e * psir
        res = abs(r.dot(r)) ** 0.5 / abs(psir.dot(psir)) ** 0.5
        print(f"maxm={maxm:3d} {name} E={complex(e):.8f} Re-Emin={e.real-emin:+.1e} "
              f"res={res:.1e} dist={np.min(abs(ev-e)):.1e} ({time.perf_counter()-t0:.0f} s)")
pnh._arnoldi_smallest_real = orig
```
```
maxm= 30 arnoldi E=-2.59077859-0.04949930j Re-Emin=-2.4e-02 res=3.3e-01 dist=5.5e-02 (17 s)
maxm= 30 arpack  E=-2.56694731-0.00003899j Re-Emin=-1.7e-04 res=1.3e-02 dist=1.7e-04 (75 s)
maxm=128 arnoldi E=-2.58848584+0.04292640j Re-Emin=-2.2e-02 res=3.1e-01 dist=4.8e-02 (19 s)
maxm=128 arpack  E=-2.56694796-0.00004170j Re-Emin=-1.7e-04 res=1.2e-02 dist=1.7e-04 (65 s)
```

The same with 30 and 60 sweeps, `x4_nh_arpack_sweeps.py`:

```
maxm= 30 nsweeps=30 arpack  E=-2.56694731-0.00003899j Re-Emin=-1.7e-04 res=1.3e-02 dist=1.7e-04 (194 s)
maxm= 30 nsweeps=60 arpack  E=-2.56694731-0.00003899j Re-Emin=-1.7e-04 res=1.3e-02 dist=1.7e-04 (419 s)
```

**Suggested fix**: port the left/right truncation and biorthogonalization
of ITensorNHDMRG.jl's "onesided" route (with `biorthoalg="biorthoblock"`)
into both ports, keep the tie-break anchor, and replace the single-vector
restart of the local Arnoldi by a thick (Krylov-Schur) restart; check on
this chain at N=10 against full ED and against julia_live, and on the
existing `tests/test_nh_dmrg.py` models.

### 4. `"python"`'s random start MPS overflows at maxm=128 from about 300 sites, and the solve raises `LinAlgError('SVD did not converge')` before its first sweep

`bug` &middot; severity **LOW** &middot; OPEN &middot; found by the benchmark

**Where**: `src/dmrgpy/pyitensor/mpsalgebra.py::randomMPS` (lines
1032-1050), which every `"python"` solve without a sector starts from
(`chain.py::_default_mps`).

`randomMPS` draws Gaussian tensors at full bond dimension, canonicalizes
them in one left-to-right `position(n)` sweep and normalizes only at the
end. The norm piles up on the moving center, a factor of order
`sqrt(d*m)` per site, so at `m=128` it passes the double range after a few
hundred sites; the overflow warnings come first and the SVD of the next
center then fails. In the benchmark the correlation-matrix case raised at
N=448 after 10 s (at N=320 the same case took 541 s), and locally
`randomMPS` already raises at N=320 with this seed. The failure is loud,
and only long chains at large `maxm` reach it.

**Expected**: a finite, normalized random MPS at any length.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 <scratch>/r4_random_mps.py
```
```python
# pyitensor's random start MPS at maxm=128: is it finite, and does the
# 448-site correlation-matrix case fail before any sweep?
from common import np, HEAVY
import time, calctk
from dmrgpy.pyitensor.mpsalgebra import randomMPS
for n in (320, 448):
    np.random.seed(1234)
    fc = calctk.make_spinless_chain(n, "python", HEAVY)
    sites = fc._session.sites
    psi = randomMPS(sites, 128)
    finite = all(np.all(np.isfinite(np.asarray(psi.A(i).array))) for i in range(1, n + 1))
    print(f"N={n:3d} randomMPS(m=128) all entries finite: {finite}")
np.random.seed(1234)
fc = calctk.make_spinless_chain(448, "python", HEAVY)
t0 = time.perf_counter()
try:
    fc.gs_energy(mode="DMRG")
    print("N=448 gs_energy finished")
except Exception as ex:
    print(f"N=448 gs_energy raised after {time.perf_counter()-t0:.1f} s: {ex!r}")
```

Observed (the traceback, repo paths shortened to `<repo>`):

```
Traceback (most recent call last):
  File "<scratch>/r4_random_mps.py", line 10, in <module>
    psi = randomMPS(sites, 128)
  File "<repo>/src/dmrgpy/pyitensor/mpsalgebra.py", line 1049, in randomMPS
    mps.position(n)  # lossless left-to-right sweep: exact canonicalization
    ~~~~~~~~~~~~^^^
  File "<repo>/src/dmrgpy/pyitensor/mpscontainer.py", line 105, in position
    self._shift_right(self.center, cutoff, maxdim)
    ~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "<repo>/src/dmrgpy/pyitensor/mpscontainer.py", line 124, in _shift_right
    U, S, V, spec = svd(T, left_inds, cutoff=cutoff, maxdim=maxdim)
                    ~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "<repo>/src/dmrgpy/pyitensor/svd.py", line 235, in svd
    U, S, S_host, Vh, keep, discarded = _svd_truncated(mat, cutoff, maxdim,
                                        ~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^
                                                       mindim)
                                                       ^^^^^^^
  File "<repo>/src/dmrgpy/pyitensor/svd.py", line 202, in _svd_truncated
    U, S, Vh = bk.xp().linalg.svd(mat, full_matrices=False)
               ~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "<site-packages>/numpy/linalg/_linalg.py", line 1841, in svd
    u, s, vh = gufunc(a, signature=signature)
               ~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^
  File "<site-packages>/numpy/linalg/_linalg.py", line 153, in _raise_linalgerror_svd_nonconvergence
    raise LinAlgError("SVD did not converge")
numpy.linalg.LinAlgError: SVD did not converge
```

**Suggested fix**: divide the center by its norm at every step of the
canonicalizing sweep (or build the left-canonical tensors directly), which
leaves the state itself unchanged, so no number moves.
