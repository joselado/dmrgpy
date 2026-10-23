# Audit, 2026-10-23: seven-lens hole hunt over the whole Python layer

This record is the eighth hole hunt over the Python layer. Unlike the last
five, it is not scoped to a commit: each lens took one area that the earlier
hunts covered least. It ran on 2026-10-23 on `c54bbfa` (clean tree apart
from three untracked example PNGs; both compiled extensions current, no
source in `mpscpp2/` or `mpscpp3/` newer than its `.so`) as seven hunters,
one per lens, and one reviewer per candidate briefed to refute it, twenty
agents at most, in about 25 minutes. Every repro was executed through a
runner that pins the BLAS threads, puts this checkout's `src` on the path
and allows three Python processes at once across all agents. The hunters
returned seven candidates; the reviewers refuted none, narrowed four, and
lowered the severity of three. Two candidates from different lenses were the
same defect and are merged into finding 2. Three lenses found nothing, and
what they checked is under "Ruled out".

This file is the evidence, not a task list. Each entry has a `**Status**`
line saying what was done, and keeps its repro.

The four bugs were fixed on 2026-10-23 in the same session, by hand. The
regressions are `tests/test_audit_2026_10_23.py`: 15 of its 18 tests fail
on `c54bbfa`, and the other three are controls. No fix changes a number that
a correct call returned before: every fix either raises where a wrong number
came out, or replaces a wrong number by the right one, as recorded in
`docs/numerical_changes.md`. The performance finding is recorded and not
applied.

## The seven lenses

| Lens | Brief |
|---|---|
| `exotic-statistics` | Bosonic, parafermionic and mixed spin-fermion chains against ED or a dense reference: energies, vev, static and dynamical correlators |
| `mps-algebra-entanglement` | Overlaps, operator application, sums of states, normalization, excited states and entanglement against ED, on v2, v3 and `"python"` |
| `finite-temperature` | Purification, METTS and every `T=` keyword against exact thermal averages from ED |
| `infinite-chain` | iDMRG and VUMPS beyond the S=1/2 Heisenberg chain: spin 1, XX in a field, free fermions, dimerized and longer-range couplings, observables other than the energy |
| `static-observables` | vev and static correlators with complex, non-Hermitian, same-site and Jordan-Wigner operators, spinful chains, after `set_gs` |
| `nonhermitian-generalized` | NH-DMRG and `gs_energy_generalized`: which eigenvalue is selected, PT-symmetric chains, the exceptional point, left and right states |
| `performance` | The largest avoidable cost of `"python"` or ED on an ordinary 12 to 30 site calculation, with a measured prototype |

## Scope

Out of scope by construction: the vendored ITensor and TDVP folders; the
three legacy quirks `.claude/rules/cpp-backend.md` lists as deliberately
reproduced; the open `docs/known_issue_*.md` items; anything already in an
earlier `docs/audit_*.md` record, including their "Left open" and "New
leads" sections; and what `ROADMAP.md` marks as absent by design.
`itensor_version="julia_live"` was out of scope, for its JIT cost. Every
repro below ran as

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 <script>
```

## The findings at a glance

| # | Finding | Severity | Status |
|---|---|---|---|
| 1 | METTS with `basis_ops=("N",)` on a number-conserving fermion chain never leaves the particle-number sector of its first sample: <Ntot> 3, 1, 1, 4 by seed, stderr 0, for an exact 1.7811 | HIGH | FIXED (raises) |
| 2 | DMRG `vev()` ignores every keyword: `vev(O, T=2.0)` returns the T=0 value -0.2277 for an exact -0.0345; `excited_vev(O, n=3)` returns the ground-state value on both modes | MEDIUM | FIXED |
| 3 | `Spin_Chain.get_correlator(name="SS")` ignores `name` and returns <Sz Sz>: -0.236 for an exact <S.S> of -0.664 | MEDIUM | FIXED (raises) |
| 4 | v2 and v3 bond entropies of a state that is not unit norm are wrong, up to 3.5 too high at norm^2 4, while `"python"` is exact | LOW | FIXED |
| 5 | `"python"` QR through `numpy.linalg.qr` instead of `scipy.linalg.qr`: about 1.3x per call, about 1.1x on a KPM correlator | LOW, performance | Open |
| 6 | `get_pair_entropy` and `get_mutual_information` do not normalize: they raise on `2*wf`, and at norm^2 1.0004 return 1.1e-4 too much | LOW | FIXED |
| 7 | `entropy_dm` drops eigenvalues below an absolute 1e-6: site entropy 5.6e-7 for an exact 8.66e-6 | LOW | FIXED |
| 8 | An unconverged NH-DMRG state is stored and reused silently, residual 0.1 to 0.4 | LOW to MEDIUM | FIXED (warns on reuse) |
| 9 | `gs_energy_generalized` never checks its metric: A=-1 returns the largest eigenvalue, a singular A aborts v3 | LOW to MEDIUM | FIXED (raises) |
| 10 | METTS with a diagonal basis whose own sum is not conserved stays in one charge sector: <Ntot> 6, 4, 2 by seed, stderr 0, for 2.74 | MEDIUM | FIXED (raises) |
| 11 | `"python"` DMRG keeps the symmetry sector of its first local solve: -4.470573 for -5.917370 on two equivalent spins | HIGH | FIXED |

## Findings

### 1. METTS with `basis_ops=("N",)` on a spinless fermion chain that conserves N returns the value of one particle-number sector, chosen by the seed, with an error bar of 0

`bug` &middot; severity **HIGH** &middot; CONFIRMED (narrowed) &middot; lens `finite-temperature`

**Status**: FIXED on 2026-10-23. `vevtk/mettsvev.py::_check_basis` runs
before the backend on `"python"` and v3. It raises when a basis operator is
not Hermitian (the collapse calls `eigh` on it), and when every basis
operator is diagonal and H conserves `sum_i O_i` for one of them, tested on
a random witness state as `mpsalgebra.is_hermitian` does (the canonical form
cannot prove `[H,N]=0` for hopping, since that rests on same-site identities
such as `N C = 0`). The same guard catches `("Sz",)` alone on a Heisenberg
chain. `("Sz","Sx")` and a diagonal basis on a chain that does not conserve
it (transverse-field Ising with `("Sz",)`) still run. For a spinless fermion
chain that conserves N there is no valid basis, so METTS now refuses it and
the message points to `mode="ED"`. Pinned by
`tests/test_audit_2026_10_23.py::test_metts_refuses_*` and the control
`test_metts_accepts_basis_that_breaks_the_charge`. A call that raises now
returned a wrong number before.

**Where**: `src/dmrgpy/vevtk/mettsvev.py::metts_vev`; the collapse in
`src/dmrgpy/pyitensor/metts.py` and `mpscpp3`'s `Chain::metts_vev`.

The collapse projects each site onto an eigenstate of the basis operator.
With only `N`, every sample is a product state of definite particle number,
and imaginary-time evolution under a number-conserving H keeps that number,
so the Markov chain never leaves the sector of its random first sample. The
default basis `("Sz","Sx")` raises on a fermion site, so a user is pushed to
`("N",)`, the one choice that gives a silently wrong answer. The variance of
`Ntot` is exactly zero, so the reported error bar is zero too.

**Expected**: the grand-canonical thermal average, <Ntot> = 1.7811 and a
variance of 0.726 at T=1, or an error.

Repro:

`hunt_finite-temperature/06_metts_fermion.py`:

```python
import numpy as np
from dmrgpy import fermionchain
n=4
for v in ("python",3):
  for seed in (1,2,3,4):
    sc=fermionchain.Fermionic_Chain(n,itensor_version=v)
    sc.set_hamiltonian(sum(sc.Cdag[i]*sc.C[i+1]+sc.Cdag[i+1]*sc.C[i] for i in range(n-1))+0.3*sum(sc.N[i] for i in range(n)))
    Ntot=sum(sc.N[i] for i in range(n))
    ex=sc.vev(Ntot,mode="ED",T=1.0).real
    m=sc.metts_vev([Ntot,Ntot*Ntot],T=1.0,nsamples=300,nwarmup=20,seed=seed,basis_ops=("N",))
    print(v,"seed",seed,"exact <Ntot>=%.4f"%ex,"METTS <Ntot>=%.4f +- %.4f"%(m[0][0].real,m[0][1]),"METTS <Ntot^2>-<Ntot>^2=%.4f"%(m[1][0].real-m[0][0].real**2))
# exact variance
sc=fermionchain.Fermionic_Chain(n,mode="ED")
sc.set_hamiltonian(sum(sc.Cdag[i]*sc.C[i+1]+sc.Cdag[i+1]*sc.C[i] for i in range(n-1))+0.3*sum(sc.N[i] for i in range(n)))
N=sum(sc.N[i] for i in range(n))
print("exact <N^2>-<N>^2 =",(sc.vev(N*N,T=1.0)-sc.vev(N,T=1.0)**2).real)
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 06_metts_fermion.py
```

Observed:

```
python seed 1 exact <Ntot>=1.7811 METTS <Ntot>=3.0000 +- 0.0000 METTS <Ntot^2>-<Ntot>^2=-0.0000
python seed 2 exact <Ntot>=1.7811 METTS <Ntot>=1.0000 +- 0.0000 METTS <Ntot^2>-<Ntot>^2=0.0000
python seed 3 exact <Ntot>=1.7811 METTS <Ntot>=1.0000 +- 0.0000 METTS <Ntot^2>-<Ntot>^2=0.0000
python seed 4 exact <Ntot>=1.7811 METTS <Ntot>=4.0000 +- 0.0000 METTS <Ntot^2>-<Ntot>^2=-0.0000
3 seed 1 exact <Ntot>=1.7811 METTS <Ntot>=4.0000 +- 0.0000 METTS <Ntot^2>-<Ntot>^2=0.0000
3 seed 2 exact <Ntot>=1.7811 METTS <Ntot>=0.0000 +- 0.0000 METTS <Ntot^2>-<Ntot>^2=0.0000
3 seed 3 exact <Ntot>=1.7811 METTS <Ntot>=2.0000 +- 0.0000 METTS <Ntot^2>-<Ntot>^2=0.0000
3 seed 4 exact <Ntot>=1.7811 METTS <Ntot>=2.0000 +- 0.0000 METTS <Ntot^2>-<Ntot>^2=0.0000
exact <N^2>-<N>^2 = 0.7259518997798393
```


The reviewer's probe of the other basis choices:

`review_finite-temperature_1/02_basis.py`:

```python
from dmrgpy import fermionchain
n=4
sc=fermionchain.Fermionic_Chain(n,itensor_version="python")
sc.set_hamiltonian(sum(sc.Cdag[i]*sc.C[i+1]+sc.Cdag[i+1]*sc.C[i] for i in range(n-1)))
for b in [("Sz","Sx"),("N",),("N","Cdag"),("N","Sx"),("Sx",),("N","C")]:
  try:
    m=sc.metts_vev(sc.N[0],T=1.0,nsamples=20,nwarmup=5,seed=1,basis_ops=b); print(b,m)
  except Exception as e: print(b,"ERR",type(e).__name__,str(e)[:120])
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 02_basis.py
```

Observed:

```
('Sz', 'Sx') ERR ValueError Operator 'Sz' not recognized for site type FermionSite
('N',) (np.complex128(0.9177943746241498+0j), np.float64(0.01410160769875215))
('N', 'Cdag') (np.complex128(0.6098642502046879+0j), np.float64(0.06376662067268915))
('N', 'Sx') ERR ValueError Operator 'Sx' not recognized for site type FermionSite
('Sx',) ERR ValueError Operator 'Sx' not recognized for site type FermionSite
('N', 'C') (np.complex128(0.9177943746241498+0j), np.float64(0.01410160769875215))
```


**Reviewer (CONFIRMED, narrowed)**: re-ran the repro and got the same
output byte for byte. No audit record, known issue, guide page or docstring
says METTS needs a basis that breaks the conserved charge. The result
changes with the seed and the variance is exactly 0, so this is not
convergence. Narrowed: `("N",)` is not the only basis accepted, since
`("N","Cdag")` and `("N","C")` also run, because nothing checks that a basis
operator is Hermitian; `("N","Cdag")` gave an unvalidated 0.61 +- 0.06.

**Suggested fix**: raise when the basis cannot change a charge that H
conserves, and check that each basis operator is Hermitian.

### 2. On every DMRG backend `vev()` ignores every keyword, so `vev(O, T=2.0)` returns the T=0 value, and `excited_vev(O, n=3)` returns the ground-state value on both modes

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED (two candidates, merged) &middot; lenses `finite-temperature`, `static-observables`

**Status**: FIXED on 2026-10-23. `vev.multi_vev` takes `wf`, `npow` and `T`
and no `**kwargs`, so an unknown keyword raises `TypeError`, and a nonzero
`T` raises `NotImplementedError` naming `mode="ED"` and `metts_vev()`.
`Many_Body_Chain.excited_vev_MB(MO, n=2, mode="DMRG")` now takes the states
of `get_excited_states(mode=mode, n=n)` and returns `<n|MO|n>` for each, on
both modes; `Fermionic_Chain` uses the same method, so its ED route no longer
fails with a matmul error; `vev.excited_vev` is removed. NUMBERS CHANGE for
`excited_vev`: on the 4-site XXZ chain (Jz=0.7) it returned the ground-state
value three times and now returns one value per state, equal to ED to 1e-5.
Pinned by `test_dmrg_vev_refuses_temperature_and_unknown_keywords` and
`test_excited_vev_is_per_state`.

**Where**: `src/dmrgpy/vev.py::multi_vev` (`**kwargs` never read),
`vev.excited_vev`, `manybodychain.py` `vev()` and `excited_vev()`,
`fermionchain.py::excited_vev_spinless`.

`Many_Body_Chain.vev` forwards its keywords to `multi_vev`, which takes
`**kwargs` and never reads them. Only `EDchain.vev` consumes `T`. This is
the shape `docs/development.md` section 4.10 warns about: a `**kwargs` with
no consumer. `excited_vev` passed `excited=True` into the same sink, so it
returned the ground-state vev whatever `n` was, and its `mode=` was dropped
too.

**Expected**: the thermal value -0.0345 or an error for `T=2`; an error for
an unknown keyword; one value per excited state for `excited_vev`.

Repro (`finite-temperature`):

`hunt_finite-temperature/03_vev_T.py`:

```python
import numpy as np
from dmrgpy import spinchain
n=4
for v in (2,3,"python"):
    sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1] for i in range(n-1)))
    O=sc.Sz[0]*sc.Sz[1]
    print(v,"vev T=0:",sc.vev(O).real,"vev T=2 (default mode):",sc.vev(O,T=2.0).real,
          "vev bogus=7:",sc.vev(O,bogus=7).real if v==3 else "")
# independent exact value, plain numpy
sx=np.array([[0,1],[1,0]])/2; sy=np.array([[0,-1j],[1j,0]])/2; sz=np.diag([.5,-.5]); I=np.eye(2)
def op(m,i):
    out=np.eye(1)
    for k in range(n): out=np.kron(out,m if k==i else I)
    return out
H=sum(op(m,i)@op(m,i+1) for i in range(n-1) for m in (sx,sy,sz))
e,V=np.linalg.eigh(H); O=op(sz,0)@op(sz,1)
for T in (2.0,):
    w=np.exp(-(e-e[0])/T); w/=w.sum()
    print("exact thermal T=%g:"%T,np.sum(w*np.einsum('ji,jk,ki->i',V.conj(),O,V)).real,"exact T=0:",(V[:,0].conj()@O@V[:,0]).real)
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 03_vev_T.py
```

Observed:

```
2 vev T=0: -0.2276709006307321 vev T=2 (default mode): -0.2276709006307321 vev bogus=7: 
3 vev T=0: -0.22767090063074086 vev T=2 (default mode): -0.22767090063074086 vev bogus=7: -0.22767090063074086
python vev T=0: -0.2276709006307395 vev T=2 (default mode): -0.2276709006307395 vev bogus=7: 
exact thermal T=2: -0.0344746681073044 exact T=0: -0.22767090063073975
```


Repro (`static-observables`):

`hunt_static-observables/04_kw.py`:

```python
import numpy as np
from dmrgpy import spinchain
n=4
sc=spinchain.Spin_Chain([2]*n,itensor_version=3)
h=0
for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
O=sc.Sz[0]*sc.Sz[1]
print("ED T=0     ",sc.vev(O,mode="ED").real)
print("ED T=2     ",sc.vev(O,mode="ED",T=2.0).real)
print("DMRG       ",sc.vev(O).real)
print("DMRG T=2   ",sc.vev(O,T=2.0).real)
print("DMRG bogus ",sc.vev(O,bogus=3).real)
print("DMRG excited=False/ npow=2 ok;  mode typo:")
try: print(sc.vev(O,mode="dmrg"))
except Exception as e: print("raised",e)
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 04_kw.py
```

Observed:

```
ED T=0      -0.2276709006307398
ED T=2      -0.0344746681073044
DMRG        -0.22767090063074466
DMRG T=2    -0.22767090063074466
DMRG bogus  -0.22767090063074466
DMRG excited=False/ npow=2 ok;  mode typo:
raised unrecognized mode 'dmrg': expected one of 'ED' or 'DMRG'
```


`hunt_static-observables/05_excited.py`:

```python
import numpy as np
from dmrgpy import spinchain, fermionchain
n=4
sc=spinchain.Spin_Chain([2]*n,itensor_version=3)
h=0
for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
O=sc.Sz[0]*sc.Sz[1]
print("gs vev DMRG    ",sc.vev(O).real)
for k in (dict(),dict(n=3)):
    try: print("excited_vev DMRG",k,np.ravel(sc.excited_vev(O,**k)))
    except Exception as e: print("excited_vev DMRG raised",k,type(e).__name__,e)
    try: print("excited_vev ED  ",k,np.ravel(sc.excited_vev(O,mode="ED",**k)))
    except Exception as e: print("excited_vev ED raised",k,type(e).__name__,e)
fc=fermionchain.Fermionic_Chain(n,itensor_version=3)
fc.set_hamiltonian(sum(fc.Cdag[i]*fc.C[i+1]+fc.Cdag[i+1]*fc.C[i] for i in range(n-1)))
for md in ("DMRG","ED"):
    try: print("ferm excited_vev",md,np.ravel(fc.excited_vev(fc.N[0],mode=md,n=3)))
    except Exception as e: print("ferm raised",md,type(e).__name__,e)
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 05_excited.py
```

Observed:

```
gs vev DMRG     -0.22767090063073872
excited_vev DMRG {} [-0.2276709+0.j]
excited_vev ED   {} [-0.2276709+0.j]
excited_vev DMRG {'n': 3} [-0.2276709+0.j]
excited_vev ED   {'n': 3} [-0.2276709+0.j]
ferm excited_vev DMRG [0.5+0.j]
ferm raised ED ValueError matmul: dimension mismatch with signature (n,k=16),(k=3,1?)->(n,1?)
```


**Reviewer (CONFIRMED, severity lowered from HIGH to MEDIUM)**: both
reviewers reproduced the numbers. The guide shows `T=` only together with
`mode="ED"`, so a DMRG `T=` is not a documented promise, but nothing rejects
it either, and the answer is wrong by a factor 6.6. `excited_vev` is
undocumented and untested, and was wrong on every route.

**Suggested fix**: no `**kwargs` in `multi_vev`, raise on a nonzero `T`,
and compute `excited_vev` from the excited states.

### 3. `Spin_Chain.get_correlator(pairs, name="SS")` ignores `name` and returns <Sz_i Sz_j>, -0.236 for an exact <S.S> of -0.664

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED (narrowed) &middot; lens `static-observables`

**Status**: FIXED on 2026-10-23. `Many_Body_Chain.get_correlator` raises
`TypeError` on `name=` and keeps its old default (<Sz Sz> on a spin chain,
<Cdag C> on a spinless fermion chain), so no call without `name` changes.
Delegating to `correlator.get_correlator` was not done, since it would
change the no-name default from <Sz Sz> to <S.S>. Pinned by
`test_get_correlator_refuses_name`.

**Where**: `src/dmrgpy/manybodychain.py::get_correlator`.

The method is a deprecated stub. It builds `Sz[i]*Sz[j]` for a spin chain
and passes `**kwargs` on to `vev()`, which swallowed `name` (finding 2). The
ED route uses the same stub, so a DMRG against ED comparison cannot see it.
`correlator.get_correlator`, which does read `name`, is only reached from
the fermion chains.

**Expected**: <S.S> for `name="SS"`, or an error.

Repro:

`hunt_static-observables/02_name_ignored.py`:

```python
import numpy as np
from dmrgpy import spinchain
n=4
sc=spinchain.Spin_Chain([2]*n,itensor_version=3)
h=0
for i in range(n-1):
    h=h+sc.Sx[i]*sc.Sx[i+1]+0.5*sc.Sy[i]*sc.Sy[i+1]+1.5*sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
pairs=[(0,1),(1,3)]
print(" name   get_correlator(name=...)   explicit vev (ED)")
for nm,mk in [("SS",lambda i,j: sc.Sx[i]*sc.Sx[j]+sc.Sy[i]*sc.Sy[j]+sc.Sz[i]*sc.Sz[j]),
              ("XX",lambda i,j: sc.Sx[i]*sc.Sx[j]),("YY",lambda i,j: sc.Sy[i]*sc.Sy[j]),("ZZ",lambda i,j: sc.Sz[i]*sc.Sz[j])]:
    got=np.real(sc.get_correlator(pairs=pairs,name=nm))
    ref=np.real([sc.vev(mk(i,j),mode="ED") for i,j in pairs])
    print(nm,got,ref)
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 02_name_ignored.py
```

Observed:

```
 name   get_correlator(name=...)   explicit vev (ED)
SS [-0.23617323  0.10955944] [-0.66366708  0.19132062]
XX [-0.23617323  0.10955944] [-0.21677815  0.05672955]
YY [-0.23617323  0.10955944] [-0.2107157   0.02503162]
ZZ [-0.23617323  0.10955944] [-0.23617323  0.10955944]
```


**Reviewer (CONFIRMED, severity lowered to MEDIUM)**: reproduced the
numbers exactly; `XX` and `YY` also return the `ZZ` values. No test, guide
page or README calls it on a spin chain, and it prints a deprecation notice,
hence MEDIUM. Struck: the hunter's guess that the tests pass because they
use `ZZ`.

**Suggested fix**: raise on `name=`, keeping the old default.

### 4. On v2 and v3 the bond entropy of a state that is not unit norm is wrong, up to 3.5 too high for a state scaled by 2, while `"python"` is exact

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `mps-algebra-entanglement`

**Status**: FIXED on 2026-10-23. `entropy.compute_entropy_single`
normalizes the state before the compiled session cuts it, raises on a zero
state, and leaves a unit-norm state untouched. After the fix the repro
below gives zero error on v2, v3 and `"python"` for the scaled state.
NUMBERS CHANGE only for a state that was not unit norm, where the old number
was wrong: on the 8-site XXZ chain (Jz=0.7) at norm^2 4, the bond entropies
were 2.08 to 3.52 above the exact ones and now agree to 1e-10. Pinned by
`test_bond_entropy_is_scale_invariant`.

**Where**: `src/dmrgpy/entropy.py::compute_entropy_single`.

The compiled sessions take `-sum w ln w` over the raw Schmidt weights, so
the result depends on the norm. `get_site_entropy` normalizes and
`"python"` normalizes, so only v2 and v3 bond entropies were affected.
Ground states are unit norm, so only a state the caller builds (`2*wf`,
`Sx[0]*wf`, a sum of states) reaches it.

**Expected**: the entropy of the normalized state, equal on all backends.

Repro (before the fix):

`hunt_mps-algebra-entanglement/03_ent.py`:

```python
import numpy as np
from dmrgpy import spinchain
def ent(v,b):
    s=np.linalg.svd(v.reshape(2**b,-1),compute_uv=False)**2
    s=s/s.sum(); s=s[s>1e-14]; return -np.sum(s*np.log(s))
n=8
for ver in (2,3,"python"):
  for cut in (None,1e-14):
    sc=spinchain.Spin_Chain([2]*n,itensor_version=ver)
    h=0
    for i in range(n-1):
        h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+0.7*sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm=60; sc.nsweeps=20
    if cut: sc.cutoff=cut
    sc.gs_energy(); wf=sc.get_gs(); wfe=sc.get_gs(mode="ED")
    ve=np.array(wfe.v).ravel(); ve=ve/np.linalg.norm(ve)
    eb=np.array([ent(ve,b) for b in range(1,n)])
    mb=np.array([wf.get_bond_entropy(b-1,b) for b in range(1,n)])
    print(ver,"cutoff",sc.cutoff,"gs bond err",np.round(mb-eb,5))
    w2=2.0*wf
    print("   scaled x2 bond entropies", np.round(np.array([w2.get_bond_entropy(b-1,b) for b in range(1,n)])-eb,4),"norm2",w2.dot(w2).real)
    print("   get_entropy()",wf.get_entropy(),"mean ed",eb.mean())
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 03_ent.py
```

Observed:

```
2 cutoff 1e-12 gs bond err [ 0. -0. -0. -0. -0. -0.  0.]
   scaled x2 bond entropies [2.0794 3.5196 1.9355 3.2148 1.9355 3.5196 2.0794] norm2 3.999999999999995
   get_entropy() 0.5866415841989673 mean ed 0.5866415842008028
2 cutoff 1e-14 gs bond err [ 0. -0. -0. -0. -0. -0. -0.]
   scaled x2 bond entropies [2.0794 3.5196 1.9355 3.2148 1.9355 3.5196 2.0794] norm2 4.0000000000000036
   get_entropy() 0.5866415842000802 mean ed 0.5866415842008028
3 cutoff 1e-12 gs bond err [ 0. -0. -0. -0. -0. -0.  0.]
   scaled x2 bond entropies [2.0794 3.5196 1.9355 3.2148 1.9355 3.5196 2.0794] norm2 3.9999999999999947
   get_entropy() 0.5866415841996426 mean ed 0.5866415842008028
3 cutoff 1e-14 gs bond err [-0. -0. -0. -0. -0. -0.  0.]
   scaled x2 bond entropies [2.0794 3.5196 1.9355 3.2148 1.9355 3.5196 2.0794] norm2 4.0
   get_entropy() 0.5866415841999311 mean ed 0.5866415842008028
python cutoff 1e-12 gs bond err [ 0. -0.  0. -0.  0. -0.  0.]
   scaled x2 bond entropies [ 0. -0.  0. -0.  0. -0.  0.] norm2 3.9999999999999973
   get_entropy() 0.5866415842007761 mean ed 0.5866415842008028
python cutoff 1e-14 gs bond err [-0. -0. -0. -0. -0. -0.  0.]
   scaled x2 bond entropies [-0. -0. -0. -0. -0. -0.  0.] norm2 4.0
   get_entropy() 0.5866415842007271 mean ed 0.5866415842008028
```


**Reviewer (CONFIRMED, severity lowered from MEDIUM to LOW)**: reproduced
it at two cutoffs. Lowered for the same reason the 2026-09-25b record kept
`get_rdm`'s norm finding LOW: only a caller-built state reaches it. Struck:
a 0.0041 `"python"` error in the hunter's second probe, which was the probe
reading the ED vector with the wrong site order.

**Suggested fix**: normalize in `compute_entropy_single`.

### 5. On `"python"`, `qr_split` uses `numpy.linalg.qr`, about 1.3x slower per call than `scipy.linalg.qr`, which is worth about 1.1x on a KPM correlator

`performance` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `performance`

**Status**: open. The gain is small and inside the run-to-run spread of a
loaded host, so it was not applied.

**Where**: `src/dmrgpy/pyitensor/svd.py::qr_split`, called by the exact
zip-up sweep of `applyMPO` in the KPM recursion.

The hunter's profile put about 42% of a 16-site KPM correlator on
`"python"` in QR. Replacing `numpy.linalg.qr` with
`scipy.linalg.qr(mode="economic", check_finite=False)` gave identical
moments and a median 1.19x end to end on four runs each, with a spread of
10 to 18 s per run, wider than the effect.

Repro (microbenchmark, hunter):

`hunt_performance/09_qr.py`:

```python
import numpy as np,scipy.linalg as sl,timeit
rng=np.random.default_rng(0)
A=rng.standard_normal((500,250))+1j*rng.standard_normal((500,250))
def t(f): return min(timeit.repeat(f,number=10,repeat=5))/10*1e3
print("np.qr reduced   %.1f ms"%t(lambda:np.linalg.qr(A)))
print("sp.qr economic  %.1f ms"%t(lambda:sl.qr(A,mode='economic',check_finite=False)))
print("sp.qr econ ow   %.1f ms"%t(lambda:sl.qr(A.copy(),mode='economic',overwrite_a=True,check_finite=False)))
print("np.qr r only    %.1f ms"%t(lambda:np.linalg.qr(A,mode='r')))
print("svd gesdd       %.1f ms"%t(lambda:np.linalg.svd(A,full_matrices=False)))
B=np.asfortranarray(A)
print("sp.qr fortran ow %.1f ms"%t(lambda:sl.qr(B.copy(order='F'),mode='economic',overwrite_a=True,check_finite=False)))
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 09_qr.py
```

Observed:

```
np.qr reduced   26.0 ms
sp.qr economic  22.0 ms
sp.qr econ ow   20.4 ms
np.qr r only    9.0 ms
svd gesdd       41.4 ms
sp.qr fortran ow 14.0 ms
```


The prototype, end to end:

`hunt_performance/11_proto.py`:

```python
import time,numpy as np,scipy.linalg as sl,statistics
from dmrgpy import spinchain
from dmrgpy.pyitensor import svd as S, backend as bk
n=16
sc=None
def run():
    global sc
    if sc is not None: return go()
    sc=spinchain.Spin_Chain([2]*n,itensor_version="python")
    h=0
    for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.gs_energy()
    return go()
def go():
    t=time.time()
    r=sc.get_dynamical_correlator(name=(sc.Sz[n//2],sc.Sz[n//2]),submode="KPM",delta=0.2,es=np.linspace(0,3,20))
    return time.time()-t,np.array(r[1])
def fastqr(a,mode="reduced"):
    return sl.qr(np.array(a,order="F"),mode="economic",overwrite_a=True,check_finite=False)
class XP:  # proxy so only qr is swapped
    def __init__(s,x): s.x=x
    def __getattr__(s,k): return getattr(s.x,k)
orig=np.linalg.qr
res={}
for rep in range(4):
    for name in ("orig","fast"):
        np.linalg.qr = orig if name=="orig" else fastqr
        t,y=run(); res.setdefault(name,[]).append(t)
        if rep==0: res[name+"y"]=y
for k in("orig","fast"): print(k,["%.2f"%x for x in res[k]],"median %.2f"%statistics.median(res[k]))
print("speedup %.2fx maxdiff %.2e"%(statistics.median(res["orig"])/statistics.median(res["fast"]),abs(res["origy"]-res["fasty"]).max()))
```

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 11_proto.py
```

Observed:

```
orig ['17.64', '14.49', '17.42', '10.01'] median 15.96
fast ['14.51', '15.04', '10.66', '12.37'] median 13.44
speedup 1.19x maxdiff 0.00e+00
```


**Reviewer (CONFIRMED, narrowed)**: re-ran the microbenchmark twice on a
host with load 4.7 on 14 cores: numpy 23.9 and 21.4 ms, scipy economic 18.3
and 16.0 ms, so 1.3x, not the hunter's 1.85x. The Fortran-ordered,
overwritten copy the hunter proposed adds nothing over plain
`scipy.linalg.qr`. End to end at most about 1.1x.

**Suggested fix**: `scipy.linalg.qr(mode="economic", check_finite=False)`
in `qr_split` on the numpy backend, if a 1.1x is wanted.

## Ruled out

- `exotic-statistics`: bosonic chains (uniform and non-uniform `maxnb`),
  parafermions Z=2, 3, 4, the mixed spin-fermion chain against a dense
  Jordan-Wigner reference, and the spin-boson chain against a dense
  reference: energies, vev, KPM dynamical correlators, `get_excited`,
  conserved sectors and `evolution_ABA` agree on v3 and `"python"`, within
  the documented KPM residual where KPM is involved.
- `mps-algebra-entanglement`: overlaps of complex states, site and pair
  entropies, mutual information and `get_excited_states` (degenerate levels
  included) agree with ED on v2, v3 and `"python"`.
- `infinite-chain`: spin-1 Heisenberg, XX in a field (gapless and
  polarized), free fermions and dimerized chains, energies and correlators,
  iDMRG and VUMPS on v3 and `"python"`, agree with their exact anchors
  within truncation.
- `nonhermitian-generalized`: NH-DMRG ground and excited states on
  interacting Hatano-Nelson chains, and `gs_energy_generalized` for six
  metrics and scales from 1e-3 to 1e3, agree with ED.
- `finite-temperature`: `Thermal_Spin_Chain` agrees with ED to about 1e-7
  on S=1/2, S=1 and mixed chains; `get_dynamical_correlator(T=...)` on DMRG
  raises.

## New leads, not reviewed

These were reviewed later the same day; what came of each is under "The
leads, reviewed" at the end of this record.

- NH-DMRG past the exceptional point (PT-symmetric chain, g=1.5, n=6) prints
  that the residual tolerance was not reached (relative residual 0.2 to
  0.37) but stores the value, and later calls use it.
- `gs_energy_generalized` on v3 with a non-Hermitian or indefinite metric A
  aborts the whole process (SIGABRT) instead of raising.
- Inside the 4-fold degenerate lowest level of a PT-symmetric chain at
  g=0.3, DMRG returns a real member and ED a complex one, and `vev(Sz_i)`
  differs by backend.
- METTS on the spinful fermion chain and on bosons, where a charge is
  conserved too, was not probed; the guard of finding 1 covers them when
  the basis is all diagonal.
- `MPS.get_entropy(b=None)` averages bonds arithmetically and
  `entropy.compute_entropy(b=None)` geometrically.
- `get_pair_entropy` and `get_mutual_information` raise on a state that is
  not unit norm while `get_site_entropy` normalizes it.
- `entropy_dm`'s 1e-6 eigenvalue floor is absolute, not relative.
- The iDMRG Hubbard probe did not finish in time; `"python"` at U=0 gave
  -1.27093 against -4/pi = -1.27324 at `maxm=30`, probably truncation.
- The spin-boson chain gives near-degenerate ground-state energies off by
  up to 5e-5 on v3, not rerun with tighter sweeps.

## Found afterwards, by the check of the user guide

`evolve_and_measure` on DMRG without `wf=` took `self.wf0`, which is None
on a chain that was never solved and right after `set_hamiltonian`, so the
session failed with `AttributeError: 'NoneType' object has no attribute
'cpp_handle'`, while `evolution_ABA` solves the ground state in the same
case. **Status**: FIXED on 2026-10-23, `timedependent.evolve_and_measure_dmrg`
calls `self.get_gs()` when `wf` is None. Pinned by
`test_evolve_and_measure_solves_the_ground_state_when_none_is_given`, which
also records that v2's MPO-Taylor integrator drifts by about 4e-6 per step
at `dt=0.05` on the stationary ground state. No number changes: the call
raised before.

The same check found two more things.

- `metts_vev` and `metts_dynamical_correlator` on v3 in a conserved sector
  did not return within 3 minutes, where the guide promises an error and
  `"python"` raises one (`Chain._reject_sector`). The probable cause was
  finding 1's basis check, which builds a random witness state before the
  backend can refuse. **Status**: FIXED on 2026-10-23, both routes refuse a
  sector first, with the same `ValueError` on every backend. Pinned by
  `test_metts_refuses_a_conserved_sector`.
- `get_CFT_central_charge()` fitted a linear model with a minimizer from a
  random start. It found the same minimum as the exact solution (identical
  to 6 digits on XX chains of 8, 20 and 40 sites), so the -0.63 seen on the
  8-site Heisenberg chain is the true fit, driven by the alternating
  even/odd term, not noise. **Status**: the fit is now a linear least
  squares (`entropy.central_charge`), deterministic, same numbers; the
  guide says the result can be negative on a short chain.

The guide's remark that some real-time routes still use absolute
thresholds was checked against the records and kept: small-unit real-time
evolution and v3's `applyExp` absolute error goal are listed as open in
`CLAUDE.md`.

## Found afterwards, while fixing the six broken examples

The example rewrite left six scripts that failed as their originals did.

- `Spinon_Chain.get_gs` took the `(energies, states)` pair `mpsarnoldi`
  returns as the state, so every call failed on `wf.dot`. **Status**:
  FIXED, it unpacks the pair; the projected energy equals the Heisenberg ED
  energy to 1e-6 on 6 sites on v3 and `"python"`
  (`test_spinon_chain_reproduces_the_heisenberg_energy`).
- `idmrg.imps_sum` refused the example's states with a "non-left-canonical
  tensor" error at 2e-4 against a 1e-4 check. The example stopped iDMRG on
  the gapless Heisenberg chain after 30 iterations, and that state is not a
  fixed point; with 100 the sum returns the Heisenberg state, <Sz_0 Sz_1>
  equal to 3e-16. **Status**: the example converges its states, and the
  error message now names an unconverged input state as a cause, not only a
  bug in the code.
- `dmrgpy.maxenttk` and `dmrgpy.padetk`, imported by `"maxent"`,
  `"CVMimag"`, `get_distribution_maxent` and `reconstruct_distribution`,
  were never in the repository. **Status**: written, and tested on exact
  cases (`test_maxent_recovers_a_density_of_its_own_family`,
  `test_pade_reproduces_a_rational_function_near_the_axis`). Running
  `"CVMimag"` for the first time showed that it returned twice the weight
  of every other route (0.4929 against 0.2449 for ED on a 6-site chain),
  from a `-2 Im G/pi` where the convention is `-Im G/pi`. **Status**: FIXED,
  0.2465 now (`test_cvmimag_has_the_weight_of_the_exact_spectrum`). The
  continuation from 20 points two units above the axis smooths the peaks,
  8% low on that chain, as the guide says.

## The leads, reviewed

The nine leads above were taken up later on 2026-10-23, on `743c254`
(clean tree apart from the same three untracked PNGs, both extensions
current), as a second pass of the same process: three hunters, one per
group of leads (`nonhermitian-leads` for leads 1 to 3, `entropy-leads` for
5 to 7, `misc-leads` for 4, 8 and 9), and one reviewer per group of
candidates briefed to refute them, through the same three-slot runner.
Lead 3 and lead 8 did not reproduce, lead 9 did not reproduce as stated
but turned up finding 11, lead 5 is real but in code nothing calls, and
the other five are findings 6 to 10 below. `julia_live` was out of scope
again. The fixes landed in the same session, and the regressions are
`tests/test_audit_2026_10_23_leads.py`.

Not findings:

- Lead 3, the degenerate level of the PT-symmetric chain at g=0.3. The
  level holds four states with the same real part, -1.396126 (twice real,
  and +-0.201207i). Every DMRG state on v2, v3 and `"python"` is an
  eigenstate to a residual of 1e-14, and every `vev(Sz_i)` lies inside the
  range that `Sz_i` takes over the span of the level (+-0.2716 on site 0),
  as does ED's; the projector onto the level gives `<Sz_i>=0` for every
  site. Which member a solve returns is a free choice inside a tied level,
  and so it varies from run to run.
- Lead 8, the iDMRG Hubbard energy at U=0. It is truncation: `"python"`
  iDMRG gives -1.270737, -1.272652 and -1.272931 at `maxm` 30, 60 and 100
  for -4/pi = -1.273240, and VUMPS and v3 iDMRG at `maxm=30` agree with it
  to 2e-4. Every run carries the dimerized-cell warning.
- Lead 9 as stated, the spin-boson energies on v3. v3 matches a dense
  diagonalization, -5.91736957, to 7e-14 at `maxm` 20 and above; the 5e-5
  was `maxm=8` truncation (5e-7 there, 4e-9 at 12). Running it on
  `"python"` turned up finding 11.
- Lead 5, the two means over bonds. `MPS.get_entropy(b=None)` takes the
  arithmetic mean and agrees with the exact one to 1e-10 on all three
  backends (0.403891 on a dimerized 8-site XXZ chain).
  `entropy.compute_entropy(b=None)` took the geometric mean, 0.145644 on
  the same chain and exactly 0 as soon as one bond is unentangled, but
  nothing called it. **Status**: removed, together with its `gmean`
  (`test_compute_entropy_with_a_geometric_mean_is_gone`).

### 6. `get_pair_entropy` and `get_mutual_information` do not normalize the state: they raise on `2*wf`, and at norm^2 1.0004 they return 1.1e-4 too much

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `entropy-leads`

**Status**: FIXED, `entropy._unit_state` divides the state by its norm in
`pair_entropy` and `site_entropy` and raises `ValueError` on a zero state
(where `site_entropy` failed on the `None` that `MPS.normalize` returns).
Pinned by `test_pair_entropy_and_mutual_information_are_scale_invariant`
and `test_entropy_of_a_zero_state_raises`. No number of a unit-norm state
changes.

**Where**: `src/dmrgpy/entropy.py`, `pair_entropy`, `site_entropy`,
`mutual_information`, and the trace check of `entropy_dm`.

`site_entropy` normalized the state and `pair_entropy` did not, so the
mutual information, which calls both, raised on any state that is not unit
norm, on every backend and in both modes. Inside the 1e-3 trace tolerance
of `entropy_dm` the unnormalized pair density matrix went through, and the
entropy came out off by about 3e-4 per 1e-3 of norm error. Every solver
output is unit norm, so only a state built by the caller reaches it, the
same reason finding 4 is LOW.

**Expected**: the entropy of the ray, independent of the scale of the state,
as the bond entropy after finding 4 and the site entropy already were.

Repro:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 03_norm.py
```

`03_norm.py`:

```python
import numpy as np
from dmrgpy import spinchain
n=6
for ver in (2,3,"python"):
  for mode in ("DMRG","ED"):
    sc=spinchain.Spin_Chain([2]*n,itensor_version=ver)
    h=sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+0.7*sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
    sc.set_hamiltonian(h); sc.maxm=60; sc.nsweeps=20
    sc.gs_energy(mode=mode); wf=sc.get_gs(mode=mode)
    w2=2.0*wf
    print(ver,mode,"norm2 of 2*wf",np.real(w2.dot(w2)))
    for name,f1,f2 in [("site",lambda w:sc.get_site_entropy(w,1),None),
                       ("pair",lambda w:sc.get_pair_entropy(w,1,3),None),
                       ("MI",lambda w:sc.get_mutual_information(w,1,3),None)]:
        a=f1(wf)
        try: b=f1(w2)
        except Exception as e: b="RAISES %s: %s"%(type(e).__name__,str(e)[:70])
        print("  ",name,"unit:",a," x2:",b)
```

Observed:

```
2 DMRG norm2 of 2*wf 4.000000000000002
   site unit: 0.6931471805599452  x2: 0.6931471805599452
   pair unit: 1.2822012840095116  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000002+2.08
   MI unit: 0.10409307711037874  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000002+2.08
2 ED norm2 of 2*wf 4.0000000000000036
   site unit: 0.6931471805599454  x2: 0.6931471805599454
   pair unit: 1.2822012840095312  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000003+0j),
   MI unit: 0.10409307711035964  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000003+0j),
3 DMRG norm2 of 2*wf 4.0000000000000036
   site unit: 0.6931471805599452  x2: 0.6931471805599452
   pair unit: 1.2822012840095154  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000004-7.61
   MI unit: 0.10409307711037497  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000004-7.61
3 ED norm2 of 2*wf 4.0000000000000036
   site unit: 0.6931471805599454  x2: 0.6931471805599454
   pair unit: 1.2822012840095312  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000003+0j),
   MI unit: 0.10409307711035964  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000003+0j),
python DMRG norm2 of 2*wf 3.999999999999994
   site unit: 0.6931471805599452  x2: 0.6931471805599452
   pair unit: 1.282201284009529  x2: RAISES ValueError: this density matrix is not normalized (trace = (3.9999999999999907-1.6
   MI unit: 0.10409307711036142  x2: RAISES ValueError: this density matrix is not normalized (trace = (3.9999999999999907-1.6
python ED norm2 of 2*wf 4.0000000000000036
   site unit: 0.6931471805599454  x2: 0.6931471805599454
   pair unit: 1.2822012840095312  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000003+0j),
   MI unit: 0.10409307711035964  x2: RAISES ValueError: this density matrix is not normalized (trace = (4.000000000000003+0j),
```

The silent window is the tail of `04_floor.py` (finding 7): the pair entropy
of sites 1 and 3 of the same chain, exact 1.2822012840, is 1.2823140958 at
norm^2 1.0004 and 1.2824267701 at 1.0008, and raises from 1.0014 on.

**Reviewer (CONFIRMED, narrowed)**: the raise is real and the same on every
backend and in both modes, but it is loud; what is silent is the 1e-3
window. Narrowed from a wrong number by default to a robustness
inconsistency.

**Suggested fix**: normalize once in `pair_entropy`, as `site_entropy` does.

### 7. `entropy_dm` drops every eigenvalue below an absolute 1e-6: a site with a Schmidt weight of 5.6e-7 has entropy 5.6e-7 for an exact 8.66e-6

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `entropy-leads`

**Status**: FIXED, the floor is `1e-14` times the largest eigenvalue, which
only removes roundoff and negative eigenvalues. NUMBERS CHANGE for a site or
pair entropy with a density-matrix eigenvalue at or below 1e-6, upward by
at most about 1.4e-5 per such eigenvalue (5.625e-7 to 8.657e-6 on the chain
below). Pinned by `test_entropy_dm_keeps_small_weights`,
`test_entropy_dm_ignores_roundoff_eigenvalues` and
`test_site_entropy_of_a_nearly_product_state_against_ed`.

**Where**: `src/dmrgpy/entropy.py`, `entropy_dm`, `ds = ds[ds>1e-6]`.

A dropped weight p loses its -p ln p, at most 1.4e-5 at p=1e-6, and the
survivor 1-p still gives about p, so a nearly product site returns about p
instead of p(1+|ln p|). The density matrices are at most 4x4 (site and pair
entropies, S=1/2 and spinless fermions), so the absolute error is bounded
by about 4e-5; the relative error is about 94% for a site whose weight is
below the floor. Any state with real entanglement has eigenvalues far above
1e-6, which is why nothing caught it. The bond entropy does not go through
`entropy_dm`.

**Expected**: -sum p ln p over every eigenvalue that is not roundoff.

Repro: a 4-site chain in a staggered field of 1 with a weak XY coupling eps,
against the SVD of the ED vector,

`05_site.py`:

```python
import numpy as np
from dmrgpy import spinchain
n=4
for eps in (0.004,0.003,0.002,0.001):
    sc=spinchain.Spin_Chain([2]*n,itensor_version=3)
    h=0
    for i in range(n-1): h=h+eps*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1])
    for i in range(n): h=h+(1.0 if i%2==0 else -1.0)*sc.Sz[i]
    sc.set_hamiltonian(h)
    wfe=sc.get_gs(mode="ED"); v=np.array(wfe.v).ravel(); v/=np.linalg.norm(v)
    s=np.linalg.svd(v.reshape(2,-1),compute_uv=False)**2   # site 0
    Sx=-np.sum(s*np.log(s)); 
    print("eps",eps,"small eig %.2e"%s.min(),"exact S %.4e"%Sx,"dmrgpy get_site_entropy(site0) %.4e"%sc.get_site_entropy(wfe,0),"rel err %.2f"%((sc.get_site_entropy(wfe,0)-Sx)/Sx))
```

Observed:

```
eps 0.004 small eig 1.00e-06 exact S 1.4815e-05 dmrgpy get_site_entropy(site0) 9.9999e-07 rel err -0.93
eps 0.003 small eig 5.62e-07 exact S 8.6573e-06 dmrgpy get_site_entropy(site0) 5.6250e-07 rel err -0.94
eps 0.002 small eig 2.50e-07 exact S 4.0504e-06 dmrgpy get_site_entropy(site0) 2.5000e-07 rel err -0.94
eps 0.001 small eig 6.25e-08 exact S 1.0993e-06 dmrgpy get_site_entropy(site0) 6.2500e-08 rel err -0.94
```

**Reviewer (CONFIRMED, narrowed)**: reproduced on `diag(1-p,p)` directly,

```python
import numpy as np
from dmrgpy.entropy import entropy_dm
# exact: site with weights (1-p,p)
for p in (5.6e-7,1e-6,2e-6,1e-5):
    dm=np.diag([1-p,p]).astype(complex)
    ex=-(1-p)*np.log(1-p)-p*np.log(p)
    print(p,"exact",ex,"entropy_dm",entropy_dm(dm),"rel err",entropy_dm(dm)/ex-1)
```

```
5.6e-07 exact 8.621384113039921e-06 entropy_dm 5.599998432382783e-07 rel err -0.9350452507514109
1e-06 exact 1.4815510057992861e-05 entropy_dm 9.99999500028589e-07 rel err -0.9325032012995667
2e-06 exact 2.824472475475381e-05 entropy_dm 2.824472475475381e-05 rel err 0.0
1e-05 exact 0.00012512920464949012 entropy_dm 0.00012512920464949012 rel err 0.0
```

and narrowed to a small absolute error that is large only relative to
entropies of order 1e-5. The floor also guarded `log(0)` and negative
eigenvalues, which a relative floor keeps doing.

**Suggested fix**: a relative floor of order machine precision.

### 8. An NH-DMRG solve that does not converge is stored and reused silently: past the exceptional point, residual 0.1 to 0.4, energy up to 0.18 from every eigenvalue

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED (narrowed) &middot; lens `nonhermitian-leads`

**Status**: FIXED, the stored state carries its residual
(`self._gs_unconverged`, keyed to the state object) and
`groundstate.warn_unconverged` prints the warning again whenever
`gs_energy()`, `get_gs()` or a session reader hands that state out;
`set_gs()` or a converged solve clears it. Kept a warning rather than an
error, since the guide says `nhdmrg` warns. Pinned by
`test_unconverged_nhdmrg_state_warns_on_every_reuse` and
`test_converged_nhdmrg_state_is_reused_quietly`. No number changes.

**Where**: `src/dmrgpy/nhdmrg.py`, the end of `nhdmrg()` and
`gs_energy_nhdmrg()`; `groundstate.stored_answer_holds` and
`ground_state_on_session`, which hand out the stored state.

At g=1.5 the PT-symmetric chain `XX + 1j*g*(-1)^i Sz_i` on 6 sites is past
its exceptional point: every eigenvalue is on the imaginary axis (largest
real part 7e-15), so the smallest real part is a 64-fold tie. NH-DMRG at the
default effort does not converge on any backend, prints one warning, and
stores its best attempt as the ground state with `computed_gs=True`. Every
later `gs_energy()`, `vev()` and `get_gs()` takes it through
`stored_answer_holds` without a word. At `maxm=200`, `nsweeps=40`,
`krylovdim=60` v3 and `"python"` converge with no warning.

**Expected**: a state that did not converge is not passed off as a solved
one; at least the warning is repeated wherever it is used.

Repro (the hunter's run; the solves are unseeded, so the reviewer's rerun
below gives other numbers):

`01_past_ep.py`:

```python
import numpy as np, io, contextlib
from dmrgpy import spinchain
n=6; g=1.5
def build(version):
    sc=spinchain.Spin_Chain(["S=1/2"]*n)
    if version=="python": sc.setup_python()
    else: sc.setup_cpp(version=version)
    h=0
    for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]
    for i in range(n): h=h+1j*g*(-1)**i*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc,h
sc,h=build(3)
M=sc.get_ED_obj().MO2matrix(h); M=M.toarray() if hasattr(M,'toarray') else np.asarray(M)
ev=np.linalg.eigvals(M); ev=ev[np.argsort(ev.real)]
print("ED lowest eigenvalues (by real part):",np.round(ev[:6],5))
print("ED exceptional?: min real",ev[0].real)
def cap(f,*a,**k):
    b=io.StringIO()
    with contextlib.redirect_stdout(b): r=f(*a,**k)
    return r,b.getvalue().strip()
for v in [2,3,"python"]:
    sc,h=build(v)
    try:
        e,out=cap(sc.gs_energy)
    except Exception as ex:
        print(v,"RAISED",repr(ex)); continue
    print("== backend",v,"first gs_energy:",e)
    print("   stdout:",out or "<nothing>")
    e2,out2=cap(sc.gs_energy)
    print("   second gs_energy:",e2,"stdout:",out2 or "<nothing>")
    sz,out3=cap(sc.vev,sc.Sz[0]) if False else cap(sc.vev,sc.Sz[0])
    print("   vev(Sz0):",sz,"stdout:",out3 or "<nothing>")
    # get_gs
    w,out4=cap(sc.get_gs); print("   get_gs stdout:",out4 or "<nothing>")
    # residual of the stored right state, via ED
    wf=np.asarray(sc.wf0.get_dense()) if hasattr(sc.wf0,'get_dense') else None
    if wf is not None:
        r=M@wf-e2*wf
        print("   ||H psi - e psi||/||psi|| with stored e0:",np.linalg.norm(r)/np.linalg.norm(wf))
    print("   e0 vs nearest ED eigenvalue dist:",np.min(abs(ev-e2)))
```

Observed:

```
ED lowest eigenvalues (by real part): [-0.+2.68268j -0.+2.68268j -0.+2.56355j -0.-2.56355j -0.-1.64841j
 -0.-2.56355j]
ED exceptional?: min real -3.9948309959493555e-15
== backend 2 first gs_energy: (-0.16960700291160882-2.7960662376200385j)
   stdout: Warning: nhdmrg did not reach the residual tolerance after 5 tries (best relative residual 0.3002605082577219); consider raising nsweeps, maxm or krylovdim
   second gs_energy: (-0.16960700291160882-2.7960662376200385j) stdout: <nothing>
   vev(Sz0): (-0.45874671597937533+0j) stdout: <nothing>
   get_gs stdout: <nothing>
   e0 vs nearest ED eigenvalue dist: 0.17728746551800345
== backend 3 first gs_energy: (-0.03726574074563059-1.6069223417159817j)
   stdout: Warning: nhdmrg did not reach the residual tolerance after 5 tries (best relative residual 0.21950766115370612); consider raising nsweeps, maxm or krylovdim
   second gs_energy: (-0.03726574074563059-1.6069223417159817j) stdout: <nothing>
   vev(Sz0): (-0.3892993011941425+2.7755575615628914e-17j) stdout: <nothing>
   get_gs stdout: <nothing>
   e0 vs nearest ED eigenvalue dist: 0.05576736448104086
== backend python first gs_energy: (-0.04544335802891021+4.050045723467864j)
   stdout: Warning: nhdmrg did not reach the residual tolerance after 5 tries (best relative residual 0.15994803990221412); consider raising nsweeps, maxm or krylovdim
   second gs_energy: (-0.04544335802891021+4.050045723467864j) stdout: <nothing>
   vev(Sz0): (0.47045292570303615+1.3877787807814457e-17j) stdout: <nothing>
   get_gs stdout: <nothing>
   e0 vs nearest ED eigenvalue dist: 0.045548277843159186
```

**Reviewer (CONFIRMED, narrowed)**: the rerun of the same script gave

```
ED lowest eigenvalues (by real part): [-0.+2.68268j -0.+2.68268j -0.+2.56355j -0.-2.56355j -0.-1.64841j
 -0.-2.56355j]
ED exceptional?: min real -3.9948309959493555e-15
== backend 2 first gs_energy: (-4.527255606809075e-05-2.8481255875494482j)
   stdout: Warning: nhdmrg did not reach the residual tolerance after 5 tries (best relative residual 0.10185255431874717); consider raising nsweeps, maxm or krylovdim
   second gs_energy: (-4.527255606809075e-05-2.8481255875494482j) stdout: <nothing>
   vev(Sz0): (-0.4585421664262995+1.3877787807814457e-17j) stdout: <nothing>
   get_gs stdout: <nothing>
   e0 vs nearest ED eigenvalue dist: 0.0004447058682909525
== backend 3 first gs_energy: (-0.1033295699283077+2.467490741769158j)
   stdout: Warning: nhdmrg did not reach the residual tolerance after 5 tries (best relative residual 0.374559481095739); consider raising nsweeps, maxm or krylovdim
   second gs_energy: (-0.1033295699283077+2.467490741769158j) stdout: <nothing>
   vev(Sz0): (0.4304682867518357-2.7755575615628914e-17j) stdout: <nothing>
   get_gs stdout: <nothing>
   e0 vs nearest ED eigenvalue dist: 0.14108488822271964
== backend python first gs_energy: (-0.05000520501898697+2.873561548312339j)
   stdout: Warning: nhdmrg did not reach the residual tolerance after 5 tries (best relative residual 0.23510660910558964); consider raising nsweeps, maxm or krylovdim
   second gs_energy: (-0.05000520501898697+2.873561548312339j) stdout: <nothing>
   vev(Sz0): (0.3582581429780045-6.938893903907228e-18j) stdout: <nothing>
   get_gs stdout: <nothing>
   e0 vs nearest ED eigenvalue dist: 0.056304616574155156
```

so on v2 the energy came 4e-4 from an eigenvalue at a state residual of
0.10; the sound claim is a stored state of relative residual 0.1 to 0.4 and
an energy off by up to 0.18. The one print is the documented contract
(`docs/guide/excited-states.md`: "nhdmrg warns when it does not
converge"); what is not documented is that the attempt is then reused in
silence. The reviewer advised repeating the warning rather than raising.

**Suggested fix**: mark the stored state as unconverged and warn on reuse.

### 9. `gs_energy_generalized` never checks its metric: A=-1 returns the largest eigenvalue, an indefinite or non-Hermitian A a number that is no eigenvalue, and a singular A aborts the interpreter on v3

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED (narrowed) &middot; lens `nonhermitian-leads`

**Status**: FIXED. `groundstate.gs_energy_generalized` raises `ValueError`
for a metric that `is_hermitian` does not prove or witness Hermitian, before
any backend; the guard in the sweep loop of `pyitensor/dmrg.py` and of
`Chain::gs_energy_generalized` requires `<psi|A|psi>` to be positive, not
only nonzero, and raises `RuntimeError`; and the four ITensor `Error()`
calls of v3's `gs_energy_generalized` and `nhdmrg_generalized`, which call
`abort()`, are `std::runtime_error` throws. A positive definite metric is
untouched. Pinned by `test_generalized_refuses_a_non_hermitian_metric`,
`test_generalized_refuses_a_metric_that_is_not_positive` (the same v3
singular case aborted with rc=134 in the reviewer's subprocess repro below;
the test itself was first run after the rebuild) and
`test_generalized_with_a_positive_metric_matches_ed`.

**Where**: `src/dmrgpy/groundstate.py`, `gs_energy_generalized`;
`src/dmrgpy/pyitensor/dmrg.py`, `dmrg_generalized`;
`src/dmrgpy/mpscpp3/chain_session.h`, `Chain::gs_energy_generalized` and
`Chain::nhdmrg_generalized`.

The docstring and the guide ask for a Hermitian positive definite A, and
nothing checked it. The only guard was `|<psi|A|psi>|>1e-14`, so a negative
`<psi|A|psi>` went through, and on v3 the guard itself was an ITensor
`Error()`, which calls `abort()` and kills the interpreter, while
`"python"` raises `RuntimeError` for the same input. Measured on a 6-site
Heisenberg chain against `scipy.linalg.eig(H, A)`:

| A | exact lowest | v3 | `"python"` |
|---|---|---|---|
| `1+0.5*Sz0` | -2.9431 | -2.9431 | -2.9431 |
| `0.2+Sz0+Sz1` (indefinite) | -11.9885 | -1.3593 | -1.3656 |
| `1+1j*Sz0` (not Hermitian) | -1.7092-0.7299i | -1.84 to -1.86 | -1.37 to -1.93 |
| `-1` | -1.25 | 2.4936 | 2.4936 |
| `Sz0+0.5` (singular) | | SIGABRT | `RuntimeError` |

Repro, one case per subprocess (`python3 03_gen_case.py <version> <kind>`):

`03_gen_case.py`:

```python
import sys, numpy as np
from dmrgpy import spinchain
v=sys.argv[1]; kind=sys.argv[2]
v = v if v=="python" else int(v)
n=6
sc=spinchain.Spin_Chain(["S=1/2"]*n)
if v=="python": sc.setup_python()
else: sc.setup_cpp(version=v)
h=0
for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
if kind=="nonherm":   A=1+0.3*(sc.Sx[0]+1j*sc.Sy[0])
elif kind=="indef":   A=0.2+sc.Sz[0]*1.0+sc.Sz[1]  # eigenvalues of A: 0.2-1.5 ... indefinite
elif kind=="negdef":  A=-1+0.0*sc.Sz[0]
elif kind=="imag":    A=1+1j*sc.Sz[0]
elif kind=="sing":    A=sc.Sz[0]+0.5
elif kind=="zero":    A=0*sc.Sz[0]+1e-30*sc.Sz[1]
elif kind=="nonhermH": A=1+0.5*sc.Sz[0]; sc.set_hamiltonian(h+0.3j*sc.Sz[1])
elif kind=="pos":     A=1+0.5*sc.Sz[0]
try:
    lam=sc.gs_energy_generalized(A)
    print("RESULT",v,kind,lam)
except BaseException as e:
    print("RAISED",v,kind,type(e).__name__,str(e)[:150])
# ED reference for A spectrum
M=np.asarray(sc.get_ED_obj().MO2matrix(sc.multioperator if False else A).todense()) if hasattr(sc.get_ED_obj().MO2matrix(A),'todense') else None
if M is not None: print("A herm:",np.allclose(M,M.conj().T),"eig A:",np.round(np.sort(np.linalg.eigvals(M).real)[[0,-1]],3))
```

Observed (the reviewer's run; `rc=134` is the abort):

```
--- 3 sing
Chain::gs_energy_generalized: <psi|A|psi> collapsed to ~0 (or NaN) mid-iteration (A may not be positive definite, or the sweep drove psi toward A's near-null-space) -- cannot form a meaningful generalized Rayleigh quotient

Chain::gs_energy_generalized: <psi|A|psi> collapsed to ~0 (or NaN) mid-iteration (A may not be positive definite, or the sweep drove psi toward A's near-null-space) -- cannot form a meaningful generalized Rayleigh quotient
/tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/822cbd6d-5d6c-44ee-a178-6671f2bbbf05/scratchpad/run3.sh: line 6: 2465326 Aborted                 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 "$@"
rc=134
--- 3 zero
Chain::gs_energy_generalized: <psi|A|psi> collapsed to ~0 (or NaN) mid-iteration (A may not be positive definite, or the sweep drove psi toward A's near-null-space) -- cannot form a meaningful generalized Rayleigh quotient

Chain::gs_energy_generalized: <psi|A|psi> collapsed to ~0 (or NaN) mid-iteration (A may not be positive definite, or the sweep drove psi toward A's near-null-space) -- cannot form a meaningful generalized Rayleigh quotient
/tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/822cbd6d-5d6c-44ee-a178-6671f2bbbf05/scratchpad/run3.sh: line 6: 2465334 Aborted                 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 "$@"
rc=134
--- 3 indef
RESULT 3 indef -1.35930959495233
A herm: True eig A: [-0.8  1.2]
rc=0
--- 3 negdef
RESULT 3 negdef 2.4935771338879267
A herm: True eig A: [-1. -1.]
rc=0
--- 3 imag
RESULT 3 imag -1.8398396978810243
A herm: False eig A: [1. 1.]
rc=0
--- 3 pos
RESULT 3 pos -2.9430868812000304
A herm: True eig A: [0.75 1.25]
rc=0
--- python sing
RAISED python sing RuntimeError dmrg_generalized: <psi|A|psi> collapsed to ~0 (or NaN) mid-iteration (A may not be positive definite, or the sweep drove psi toward A's near-null-spac
A herm: True eig A: [0. 1.]
rc=0
--- python zero
RAISED python zero RuntimeError dmrg_generalized: <psi|A|psi> collapsed to ~0 (or NaN) mid-iteration (A may not be positive definite, or the sweep drove psi toward A's near-null-spac
A herm: True eig A: [-0.  0.]
rc=0
--- python indef
RESULT python indef -1.3655665680440299
A herm: True eig A: [-0.8  1.2]
rc=0
--- python negdef
RESULT python negdef 2.493577133887925
A herm: True eig A: [-1. -1.]
rc=0
--- python imag
RESULT python imag -1.9302872259580341
A herm: False eig A: [1. 1.]
rc=0
--- python pos
RESULT python pos -2.9430868812000246
A herm: True eig A: [0.75 1.25]
rc=0
```

and the exact generalized eigenvalues, `05_gen_ed.py`:

```python
import numpy as np, scipy.linalg as sla
from dmrgpy import spinchain
n=6
sc=spinchain.Spin_Chain(["S=1/2"]*n)
h=0
for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
Hm=np.asarray(sc.get_ED_obj().MO2matrix(h).todense())
for name,A,res in [("pos 1+0.5Sz0",1+0.5*sc.Sz[0],{"v3":-2.943086881200029,"py":-2.9430868812000264}),
                   ("indef 0.2+Sz0+Sz1",0.2+sc.Sz[0]+sc.Sz[1],{"v3":-1.3593178726249613,"py":-1.3655665809159678}),
                   ("nonherm 1+0.3 S+",1+0.3*(sc.Sx[0]+1j*sc.Sy[0]),{"v3":-2.4805830169451135,"py":-2.4935771338878863}),
                   ("imag 1+1j Sz0",1+1j*sc.Sz[0],{"v3":-1.860025954148969,"py":-1.3703071598479315}),
                   ("negdef -1",-1+0*sc.Sz[0],{"v3":2.4935771338879253,"py":2.4935771338879253})]:
    Am=np.asarray(sc.get_ED_obj().MO2matrix(A).todense())
    w=sla.eig(Hm,Am,right=False); w=w[np.isfinite(w)]; w=w[np.argsort(w.real)]
    print(name,"ED generalized eigenvalues lowest/highest by real part:",np.round(w[0],4),np.round(w[-1],4),"| dmrgpy",res)
    for k,v in res.items(): print("    dist of",k,"result to nearest ED eigenvalue:",round(np.min(abs(w-v)),5))
```

```
pos 1+0.5Sz0 ED generalized eigenvalues lowest/highest by real part: (-2.9431+0j) (1.6667+0j) | dmrgpy {'v3': -2.943086881200029, 'py': -2.9430868812000264}
    dist of v3 result to nearest ED eigenvalue: 0.0
    dist of py result to nearest ED eigenvalue: 0.0
indef 0.2+Sz0+Sz1 ED generalized eigenvalues lowest/highest by real part: (-11.9885+0j) (5.3529+0j) | dmrgpy {'v3': -1.3593178726249613, 'py': -1.3655665809159678}
    dist of v3 result to nearest ED eigenvalue: 0.03244
    dist of py result to nearest ED eigenvalue: 0.02619
nonherm 1+0.3 S+ ED generalized eigenvalues lowest/highest by real part: (-2.4936+0j) (1.2501+0.0001j) | dmrgpy {'v3': -2.4805830169451135, 'py': -2.4935771338878863}
    dist of v3 result to nearest ED eigenvalue: 0.01299
    dist of py result to nearest ED eigenvalue: 0.0
imag 1+1j Sz0 ED generalized eigenvalues lowest/highest by real part: (-1.7092-0.7299j) (1-0.5j) | dmrgpy {'v3': -1.860025954148969, 'py': -1.3703071598479315}
    dist of v3 result to nearest ED eigenvalue: 0.74529
    dist of py result to nearest ED eigenvalue: 0.50251
negdef -1 ED generalized eigenvalues lowest/highest by real part: (-1.25+0j) (2.4936+0j) | dmrgpy {'v3': 2.4935771338879253, 'py': 2.4935771338879253}
    dist of v3 result to nearest ED eigenvalue: 0.0
    dist of py result to nearest ED eigenvalue: 0.0
```

**Reviewer (CONFIRMED, narrowed)**: the abort is real and deterministic,
and only at that one site; the two `Error()` calls ahead of the loop sit
behind Python guards. Non-Hermitian and indefinite metrics do not abort.
The wrong numbers are garbage in, since nothing promises a result for such
an A, so the claim is narrowed to missing validation. A proof of positive
definiteness is not cheap and was not asked for; the positive guard catches
the negative and the indefinite cases that drive `<psi|A|psi>` below 0.

**Suggested fix**: check Hermiticity at the Python entry, require a positive
`<psi|A|psi>`, and throw a C++ exception instead of `Error()`.

After the fix, the same script on every case:

```
RAISED 3 indef RuntimeError Chain::gs_energy_generalized: <psi|A|psi> = -0.797818 is not positive mid-iteration (A is not positive definite, or the sweep drove psi toward A's nea
RAISED python indef RuntimeError dmrg_generalized: <psi|A|psi> = -0.7983099575762932 is not positive mid-iteration (A is not positive definite, or the sweep drove psi toward A's near-
RAISED 3 negdef RuntimeError Chain::gs_energy_generalized: <psi|A|psi> = -1 is not positive mid-iteration (A is not positive definite, or the sweep drove psi toward A's near-null-
RAISED python negdef RuntimeError dmrg_generalized: <psi|A|psi> = -1.0000000000000004 is not positive mid-iteration (A is not positive definite, or the sweep drove psi toward A's near-
RAISED 3 imag ValueError gs_energy_generalized needs a Hermitian positive definite metric A, and this A is not Hermitian
RAISED python imag ValueError gs_energy_generalized needs a Hermitian positive definite metric A, and this A is not Hermitian
RAISED 3 sing RuntimeError Chain::gs_energy_generalized: <psi|A|psi> = 0 is not positive mid-iteration (A is not positive definite, or the sweep drove psi toward A's near-null-s
RAISED python sing RuntimeError dmrg_generalized: <psi|A|psi> = 3.925231146708865e-16 is not positive mid-iteration (A is not positive definite, or the sweep drove psi toward A's nea
RESULT 3 pos -2.9430868812000326
RESULT python pos -2.943086881200026
```

The indefinite metric is caught here because the sweep drives
`<psi|A|psi>` below 0; an indefinite A whose sweep stays on the positive
side would still go through, which is the part a cheap check cannot
promise.

### 10. METTS with a diagonal collapse basis whose own sum H does not conserve still never leaves its first charge sector: `("Nupdn",)` on a Hubbard chain gives <Ntot> 6, 4, 2 by seed with stderr 0, for an exact 2.74

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED (narrowed) &middot; lens `misc-leads`

**Status**: FIXED (raises). `_check_basis` tests, besides each basis name's
own sum, every conserved quantity the chain offers
(`get_sector_charge_operators`: `Nf` and `Sz` on the spinful chains, `Nb`
on bosons), and `Nupdn` joins `canonical._DIAGONAL`, where it was missing,
so it is no longer taken for a charge-changing operator. Pinned by
`test_metts_refuses_a_diagonal_basis_on_a_charge_conserving_chain`, five
bases on two chains and two backends. No number changes; these calls raise
now.

**Where**: `src/dmrgpy/vevtk/mettsvev.py`, `_check_basis`;
`src/dmrgpy/multioperatortk/canonical.py`, `_DIAGONAL`; the collapse is
`pyitensor/metts.py`'s `_eigenbasis` and v3's `Chain::metts_vev`.

Finding 1's guard refused an all-diagonal basis only when H conserves
`sum_i O_i` for one of the basis names. A diagonal name whose own sum is not
conserved (`F` on a Hubbard chain, the occupation projectors `N0`, `N1`,
`N2` on bosons) passed it, and so did `Nupdn`, which the name table did not
know to be diagonal. The collapse diagonalizes each name with `eigh`, which
returns the computational basis vectors inside a degenerate eigenspace of a
diagonal matrix, so every collapsed product state still has a definite
particle number, and the imaginary-time evolution keeps it. The guide
promises that the guard covers this case ("The default `basis_ops` avoids
this"), and on `Spinful_Fermionic_Chain_Native` and `Bosonic_Chain` the
default `("Sz","Sx")` does not exist, so an explicit basis, and on these
sites a diagonal one, is the only way to call it.

**Expected**: a refusal, as for `("N",)`.

Repro, 4-site native Hubbard chain (t=1, U=4) and 4-site Bose-Hubbard chain
(`maxnb=4`, U=2), T=1, against `vev(..., T=1, mode="ED")`:

`01_metts_charge.py`:

```python
import numpy as np
from dmrgpy import fermionchain, bosonchain
T=1.0
def run(name, mk, ops, bases, v):
    for b in bases:
        out=[]
        for seed in (1,2,3,4):
            sc=mk(v)
            obs=ops(sc)
            try:
                r=sc.metts_vev(obs,T=T,nsamples=60,nwarmup=10,seed=seed,basis_ops=b)
            except Exception as e:
                out=str(e)[:90]; break
            out.append(tuple(round(float(np.real(x[0])),3)+0 for x in r)+(round(float(np.real(r[0][1])),3),))
        print(name,"v=",v,"basis",b,out,flush=True)
def hub(v):
    n=4
    sc=fermionchain.Spinful_Fermionic_Chain_Native(n,itensor_version=v)
    h=0
    for i in range(n-1):
        h=h-(sc.Cdagup[i]*sc.Cup[i+1]+sc.Cdagdn[i]*sc.Cdn[i+1]); 
        h=h-(sc.Cdagup[i+1]*sc.Cup[i]+sc.Cdagdn[i+1]*sc.Cdn[i])
    for i in range(n): h=h+4.0*sc.Nup[i]*sc.Ndn[i]
    sc.set_hamiltonian(h); return sc
def hubops(sc): return [sum(sc.Ntot), sum(sc.Sz)*sum(sc.Sz)]
def bose(v):
    n=4
    sc=bosonchain.Bosonic_Chain(n,maxnb=[4]*n,itensor_version=v)
    h=0
    for i in range(n-1): h=h+sc.Adag[i]*sc.A[i+1]+sc.Adag[i+1]*sc.A[i]
    for i in range(n): h=h+2.0*sc.N[i]*sc.N[i]
    sc.set_hamiltonian(h); return sc
def bops(sc): return [sum(sc.N), sum(sc.N)*sum(sc.N)]
# ED reference
for nm,mk,ops in (("hub",hub,hubops),("bose",bose,bops)):
    sc=mk("python"); sc.mode="ED"
    print(nm,"ED T=1:",[np.real(sc.vev(o,T=T,mode="ED")) for o in ops(sc)],flush=True)
for v in ("python",3):
    run("hub",hub,hubops,[("Sz","Sx"),("Sz",),("Nupdn",),("F",),("Nupdn","Sz")],v)
    run("bose",bose,bops,[("N",),("N2",),("N0",),("N0","N1")],v)
```

Observed (each tuple is <Ntot>, <Sz_tot^2> or <Ntot^2>, and the error bar
of the first, per seed):

```
hub ED T=1: [np.float64(2.7375381737677134), np.float64(0.5302420479869359)]
bose ED T=1: [np.float64(0.8235453428927778), np.float64(1.3371388373689657)]
hub v= python basis ('Sz', 'Sx') Operator 'Sx' not recognized for site type ElectronSite
hub v= python basis ('Sz',) metts_vev: H conserves sum_i Sz_i and every basis operator is diagonal, so the samples nev
hub v= python basis ('Nupdn',) [(6.0, 0.0, 0.0), (4.0, 1.0, 0.0), (2.0, 0.0, 0.0), (6.0, 1.0, 0.0)]
hub v= python basis ('F',) [(5.0, 0.25, 0.0), (5.0, 0.25, 0.0), (5.0, 2.25, 0.0), (4.0, 0.0, 0.0)]
hub v= python basis ('Nupdn', 'Sz') [(6.0, 0.0, 0.0), (4.0, 1.0, 0.0), (2.0, 0.0, 0.0), (6.0, 1.0, 0.0)]
bose v= python basis ('N',) metts_vev: H conserves sum_i N_i and every basis operator is diagonal, so the samples neve
bose v= python basis ('N2',) [(8.0, 64.0, 0.0), (4.0, 16.0, 0.0), (2.0, 4.0, 0.0), (10.0, 100.0, 0.0)]
bose v= python basis ('N0',) [(5.0, 25.0, 0.0), (5.0, 25.0, 0.0), (3.0, 9.0, 0.0), (6.0, 36.0, 0.0)]
bose v= python basis ('N0', 'N1') [(5.0, 25.0, 0.0), (5.0, 25.0, 0.0), (3.0, 9.0, 0.0), (6.0, 36.0, 0.0)]
hub v= 3 basis ('Sz', 'Sx') Operator "Sx" name not recognized
hub v= 3 basis ('Sz',) metts_vev: H conserves sum_i Sz_i and every basis operator is diagonal, so the samples nev
hub v= 3 basis ('Nupdn',) [(7.0, 0.25, 0.0), (0.0, 0.0, 0.0), (5.0, 0.25, 0.0), (4.0, 0.0, 0.0)]
hub v= 3 basis ('F',) [(2.0, 0.0, 0.0), (4.0, 4.0, 0.0), (4.0, 1.0, 0.0), (4.0, 0.0, 0.0)]
hub v= 3 basis ('Nupdn', 'Sz') [(7.0, 0.25, 0.0), (0.0, 0.0, 0.0), (5.0, 0.25, 0.0), (4.0, 0.0, 0.0)]
bose v= 3 basis ('N',) metts_vev: H conserves sum_i N_i and every basis operator is diagonal, so the samples neve
bose v= 3 basis ('N2',) [(7.0, 49.0, 0.0), (12.0, 144.0, 0.0), (3.0, 9.0, 0.0), (6.0, 36.0, 0.0)]
bose v= 3 basis ('N0',) [(1.0, 1.0, 0.0), (12.0, 144.0, 0.0), (5.0, 25.0, 0.0), (6.0, 36.0, 0.0)]
bose v= 3 basis ('N0', 'N1') [(1.0, 1.0, 0.0), (12.0, 144.0, 0.0), (5.0, 25.0, 0.0), (6.0, 36.0, 0.0)]
```

**Reviewer (CONFIRMED, narrowed)**: reproduced on `"python"`,

```python
import numpy as np
from dmrgpy import fermionchain, bosonchain
n=4
def hub():
    sc=fermionchain.Spinful_Fermionic_Chain_Native(n,itensor_version="python")
    h=0
    for i in range(n-1):
        h=h-(sc.Cdagup[i]*sc.Cup[i+1]+sc.Cdagdn[i]*sc.Cdn[i+1]); 
        h=h-(sc.Cdagup[i+1]*sc.Cup[i]+sc.Cdagdn[i+1]*sc.Cdn[i])
    for i in range(n): h=h+4.0*sc.Nup[i]*sc.Ndn[i]
    sc.set_hamiltonian(h); return sc
sc=hub(); sc.mode="ED"
print("ED Ntot",np.real(sc.vev(sum(sc.Ntot),T=1.0,mode="ED")))
for b in [("Nupdn",),("F",),("Nup",),("Ntot",),("Nup","Ndn")]:
  for seed in (1,2,3):
    sc=hub()
    try:
        r=sc.metts_vev(sum(sc.Ntot),T=1.0,nsamples=40,nwarmup=5,seed=seed,basis_ops=b)
        print(b,seed,r)
    except Exception as e: print(b,seed,"RAISE",str(e)[:70]); break
def bose():
    sc=bosonchain.Bosonic_Chain(n,maxnb=[4]*n,itensor_version="python")
    h=0
    for i in range(n-1): h=h+sc.Adag[i]*sc.A[i+1]+sc.Adag[i+1]*sc.A[i]
    for i in range(n): h=h+2.0*sc.N[i]*sc.N[i]
    sc.set_hamiltonian(h); return sc
sc=bose(); print("ED N",np.real(sc.vev(sum(sc.N),T=1.0,mode="ED")))
for b in [("N2",),("N",)]:
  for seed in (1,2,3):
    sc=bose()
    try: print(b,seed,sc.metts_vev(sum(sc.N),T=1.0,nsamples=40,nwarmup=5,seed=seed,basis_ops=b))
    except Exception as e: print(b,seed,"RAISE",str(e)[:70]); break
```

```
ED Ntot 2.7375381737677134
('Nupdn',) 1 (np.complex128(5.999999999999993+0j), np.float64(2.473606439217309e-16))
('Nupdn',) 2 (np.complex128(3.9999999999999973+0j), np.float64(2.3152124556859323e-16))
('Nupdn',) 3 (np.complex128(1.9999999999999973+0j), np.float64(1.0958969660475593e-16))
('F',) 1 (np.complex128(4.999999999999997+0j), np.float64(2.5441511085923908e-16))
('F',) 2 (np.complex128(4.999999999999995+0j), np.float64(1.6677051211084999e-16))
('F',) 3 (np.complex128(4.9999999999999964+0j), np.float64(3.1482263771856656e-16))
('Nup',) 1 RAISE metts_vev: H conserves sum_i Nup_i and every basis operator is diagona
('Ntot',) 1 RAISE metts_vev: H conserves sum_i Ntot_i and every basis operator is diagon
('Nup', 'Ndn') 1 RAISE metts_vev: H conserves sum_i Nup_i and every basis operator is diagona
ED N 0.8235453428927778
('N2',) 1 (np.complex128(8+0j), np.float64(5.526579057060682e-16))
('N2',) 2 (np.complex128(4+0j), np.float64(7.625825054161494e-17))
('N2',) 3 (np.complex128(2+0j), np.float64(1.5498323523479857e-16))
('N',) 1 RAISE metts_vev: H conserves sum_i N_i and every basis operator is diagonal,
```

No random rotation inside a degenerate eigenspace is implemented, and no
known issue covers this. Narrowed to explicit, unusual bases: the default
basis is not affected, and where it does not exist it fails loudly. The
reviewer advised testing the chain's conserved charges rather than the
names' sums, and reading diagonality off the site matrix rather than a
name table; the fix takes the first and adds the one missing name to the
table.

**Suggested fix**: refuse an all-diagonal basis whenever H conserves a
diagonal charge of the site type.

### 11. `"python"` DMRG keeps the symmetry sector of its first local solve: two equivalent spins on a spin-boson chain give -4.470573 for an exact -5.917370, in 11 of 12 runs, at every `maxm` and noise

`bug` &middot; severity **HIGH** &middot; CONFIRMED &middot; lens `misc-leads` (found while checking lead 9)

**Status**: FIXED, `pyitensor/dmrg.py`'s `_noisy_start` mixes a random
vector of relative weight `noise` into every Lanczos start of the first
sweep of a fresh start (the chain's default `noise=1e-7` puts a component
of norm 3e-4 in), which the module docstring used to say was not needed.
18 of 18 runs give the ground state now. Three narrower versions were tried
against the full suite first. Noise on every sweep of the schedule's first
half left a short run's eigenvectors less converged (a TDVP-GSE
trajectory test moved from 1.7e-7 to 1.03e-6, past its 1e-6 bound). Noise
on a warm start moved the state inside the degenerate doublet the caller
had chosen a member of (three warm-start tests, overlap 0.99999988). And
drawing from numpy's global stream shifted every seeded random start after
it, which flipped a seeded KPM test. A component of norm `noise` instead of
`sqrt(noise)` did not escape the trap (5 of 12 runs). The generator is
reseeded at the start of every `dmrg()` call, so the same solve twice is
the same bit for bit (two `test_unpadded_runs_are_*` tests compare two runs
exactly). One older test,
`test_padded_tdvp_gse_follows_the_unpadded_trajectory`, now sets
`noise=0` on its preparation chain: the noisy first sweep reaches the same
Neel state, to 1e-15 in every `<Sz_i>`, but fills the padded bond
directions differently, and that test's gauge-sensitive comparison of a
padded and an unpadded run moved from 5e-8 to 1.7e-6, past its 1e-6
bound. NUMBERS CHANGE on
`"python"` for a Hamiltonian with an exact symmetry between two sites that
share the first two-site block, from an excited level to the ground state;
other fresh solves agree with the old ones to their convergence, and warm
starts are untouched. Pinned by
`test_python_dmrg_leaves_the_sector_of_its_first_solve`.

**Where**: `src/dmrgpy/pyitensor/dmrg.py`, `_local_ground_state` and
`_dmrg_one_sweep`, and the module docstring's paragraph on the noise term.

The chain is two S=1/2 spins with a Heisenberg coupling, both coupled to
the first of four bosonic sites through `(Sx_0+Sx_1)(A_2+A_2^dagger)`, so H
commutes exactly with the swap of the two spins. The random start has the
generic weight 0.75 on the swap-symmetric sector, but both spins sit inside
the first two-site block, where the left environment is trivial, so that
block's effective Hamiltonian commutes with the swap and its lowest
eigenvector lies entirely in one sector, chosen by a local energy computed
against a still random right basis. Every later Lanczos solve starts from
the current state, and a Krylov space grown from a vector in a conserved
sector never leaves it. `"python"` had no noise term (`Sweeps.noise` was
accepted and unused), so nothing mixed the sectors back. The answer is an
exact eigenstate, with an energy variance of 2.7e-10, so no check fires.
v2 and v3 apply noise and find the ground state.

**Expected**: -5.917370, the lowest eigenvalue of a dense diagonalization
(the chain has no ED backend).

Repro:

`03_spinboson.py`:

```python
import numpy as np
from dmrgpy import bosonchain
def build(v,mode=None):
    ns,nb=2,4
    bc=bosonchain.SpinBoson_Chain(["S=1/2"]*ns+["B"]*nb,itensor_version=v)
    h=0
    for i in range(ns,ns+nb-1): h=h+1.0*bc.Adag[i]*bc.A[i+1]
    for i in range(ns,ns+nb): h=h+1.0*(bc.N[i]-1.0)*(bc.N[i]-1.0)
    h=h+bc.Sx[0]*bc.Sx[1]+bc.Sy[0]*bc.Sy[1]+bc.Sz[0]*bc.Sz[1]
    h=h+(bc.Sx[0]+bc.Sx[1])*(bc.A[2]+bc.Adag[2])
    h=h+h.get_dagger()
    bc.set_hamiltonian(h); return bc
bc=build(3)
try:
    print("ED E0",bc.gs_energy(mode="ED"))
    ed=bc.get_ED_obj(); 
    import scipy.sparse.linalg as sl
except Exception as e: print("ED fail",repr(e))
# dense reference built independently
import scipy.sparse as sp
def dense_ref():
    sx=np.array([[0,1],[1,0]])/2; sy=np.array([[0,-1j],[1j,0]])/2; sz=np.diag([.5,-.5]); I2=np.eye(2)
    d=4; a=np.diag(np.sqrt(np.arange(1,d)),1); ad=a.T; N=np.diag(np.arange(d)); Ib=np.eye(d)
    ops=[2,2,d,d,d,d]
    def emb(o,k):
        m=np.eye(1)
        for j,dj in enumerate(ops): m=np.kron(m,o if j==k else np.eye(dj))
        return m
    H=0
    for i in range(2,5): H=H+emb(ad,i)@emb(a,i+1)
    for i in range(2,6): H=H+(emb(N,i)-np.eye(1024))@(emb(N,i)-np.eye(1024))
    for o in (sx,sy,sz): H=H+emb(o,0)@emb(o,1)
    H=H+(emb(sx,0)+emb(sx,1))@(emb(a,2)+emb(ad,2))
    H=H+H.conj().T
    return np.linalg.eigvalsh(H)
w=dense_ref(); print("dense ref lowest 5",w[:5])
for v in (3,"python"):
  for (maxm,ns_,noise) in ((20,10,None),(40,20,None),(100,40,None),(100,40,1e-5),(8,10,None),(12,10,None)):
    b=build(v); b.maxm=maxm; b.nsweeps=ns_
    if noise is not None: b.noise=noise
    e=b.gs_energy(); print(v,"maxm",maxm,"nsweeps",ns_,"noise",noise,"E",e,"err vs dense",e-w[0],flush=True)
```

Observed:

```
ED fail NotImplementedError('the ED backend for SpinBoson_Chain is not implemented yet; use mode="DMRG" (itensor_version=3 or "python") on this chain, or Bosonic_Chain for a purely bosonic one')
dense ref lowest 5 [-5.91736957 -5.91736957 -4.77760767 -4.77760767 -4.47057317]
3 maxm 20 nsweeps 10 noise None E -5.917369571447391 err vs dense 7.283063041541027e-14
3 maxm 40 nsweeps 20 noise None E -5.917369571447391 err vs dense 7.283063041541027e-14
3 maxm 100 nsweeps 40 noise None E -5.917369571447382 err vs dense 8.171241461241152e-14
3 maxm 100 nsweeps 40 noise 1e-05 E -5.917369571447395 err vs dense 6.838973831690964e-14
3 maxm 8 nsweeps 10 noise None E -5.917369033987776 err vs dense 5.374596874219151e-07
3 maxm 12 nsweeps 10 noise None E -5.91736956715765 err vs dense 4.2898138374880546e-09
python maxm 20 nsweeps 10 noise None E -4.470573168577817 err vs dense 1.446796402869646
python maxm 40 nsweeps 20 noise None E -4.470573168577821 err vs dense 1.4467964028696425
python maxm 100 nsweeps 40 noise None E -4.470573168577824 err vs dense 1.4467964028696398
python maxm 100 nsweeps 40 noise 1e-05 E -4.470573168577829 err vs dense 1.4467964028696345
python maxm 8 nsweeps 10 noise None E -4.470272915670787 err vs dense 1.4470966557766767
python maxm 12 nsweeps 10 noise None E -5.917369571345399 err vs dense 1.0206413492142019e-10
```

and the symmetry discriminator, a coupling of 0.999 on the second spin:

`09_sb_asym.py`:

```python
import numpy as np
from dmrgpy import bosonchain
def run(v,g1,verbose=False,maxm=20,ns=10,order="SSBBBB"):
    b=bosonchain.SpinBoson_Chain(["S=1/2"]*2+["B"]*4,itensor_version=v)
    h=b.Sx[0]*b.Sx[1]+b.Sy[0]*b.Sy[1]+b.Sz[0]*b.Sz[1]+(b.Sx[0]+g1*b.Sx[1])*(b.A[2]+b.Adag[2])
    for i in range(2,6): h=h+(b.N[i]-1.0)*(b.N[i]-1.0)
    for i in range(2,5): h=h+b.Adag[i]*b.A[i+1]
    h=h+h.get_dagger(); b.set_hamiltonian(h); b.maxm=maxm; b.nsweeps=ns
    return b.gs_energy()
for g1 in (1.0,0.999,0.9):
    print("python g1=",g1,[round(float(run("python",g1)),5) for _ in range(4)],flush=True)
print("v3     g1=0.9",[round(float(run(3,0.9)),5) for _ in range(2)])
```

```
python g1= 1.0 [-5.91737, -4.47057, -4.47057, -4.47057]
python g1= 0.999 [-5.9153, -5.9153, -5.9153, -5.9153]
python g1= 0.9 [-5.71058, -5.71058, -5.71058, -5.71058]
v3     g1=0.9 [-5.71058, -5.71058]
```

**Reviewer (CONFIRMED)**: the dense reference uses the same conventions,
and v3 agrees with it to 7e-14. The reviewer located the lock at the first
local solve, measuring the swap-symmetric weight of the start and of the
solution there on a smaller chain (2 spins and 2 bosons, exact -3.659464):

`01_probe.py`:

```python
import numpy as np
from dmrgpy import bosonchain
from dmrgpy.pyitensor import dmrg as D
def build(v,nb=2,g1=1.0):
    b=bosonchain.SpinBoson_Chain(["S=1/2"]*2+["B"]*nb,itensor_version=v)
    h=0
    for i in range(2,2+nb-1): h=h+b.Adag[i]*b.A[i+1]
    for i in range(2,2+nb): h=h+(b.N[i]-1.0)*(b.N[i]-1.0)
    h=h+b.Sx[0]*b.Sx[1]+b.Sy[0]*b.Sy[1]+b.Sz[0]*b.Sz[1]
    h=h+(b.Sx[0]+g1*b.Sx[1])*(b.A[2]+b.Adag[2])
    h=h+h.get_dagger()
    b.set_hamiltonian(h); b.maxm=20; b.nsweeps=10; return b
orig=D._lanczos_ground_state
log=[]
def wrap(matvec,v0,*a,**k):
    r=orig(matvec,v0,*a,**k); log.append((v0.size,float(np.real(r[0])))); return r
D._lanczos_ground_state=wrap
# swap parity of start
from dmrgpy.pyitensor import mpsalgebra as M
b=build("python")
s=b._session
print([n for n in dir(s) if 'mps' in n.lower() or 'wf' in n.lower()][:30])
for seed in range(6):
    np.random.seed(seed); log.clear()
    b=build("python"); e=b.gs_energy()
    print(seed,e,[round(x[1],5) for x in log][:12],"...",log[-1])
```

`02_sector.py`:

```python
import numpy as np
exec(open('01_probe.py').read().split("orig=")[0])
orig=D._lanczos_ground_state
rec=[]
def wrap(matvec,v0,*a,**k):
    r=orig(matvec,v0,*a,**k)
    rec.append((v0.size,float(np.real(r[0])),r[1].copy(),v0.copy())); return r
D._lanczos_ground_state=wrap
def sym_weight(v,chiR):
    t=v.reshape(2,2,-1)   # first solve: no left link, order s0,s1,R
    sw=t.transpose(1,0,2); p=0.5*(t+sw)
    return float(np.sum(abs(p)**2)/np.sum(abs(t)**2))
for g1 in (1.0,0.999):
  for seed in range(6):
    np.random.seed(seed); rec.clear()
    b=build("python",g1=g1); e=b.gs_energy()
    n0,e0,ev,x0=rec[0]
    print("g1",g1,"seed",seed,"E",round(e,6),"first-solve dim",n0,"E_loc",round(e0,5),
          "start swap-sym weight",round(sym_weight(x0,0),3),"solution swap-sym weight",round(sym_weight(ev,0),6))
```

```
g1 1.0 seed 0 E -3.659464 first-solve dim 16 E_loc -0.03681 start swap-sym weight 0.75 solution swap-sym weight 1.0
g1 1.0 seed 1 E -2.328427 first-solve dim 16 E_loc 0.01479 start swap-sym weight 0.75 solution swap-sym weight 0.0
g1 1.0 seed 2 E -2.328427 first-solve dim 16 E_loc 0.00969 start swap-sym weight 0.75 solution swap-sym weight 0.0
g1 1.0 seed 3 E -3.659464 first-solve dim 16 E_loc 0.00455 start swap-sym weight 0.75 solution swap-sym weight 1.0
g1 1.0 seed 4 E -2.328427 first-solve dim 16 E_loc 0.01797 start swap-sym weight 0.75 solution swap-sym weight 0.0
g1 1.0 seed 5 E -2.328427 first-solve dim 16 E_loc 0.06049 start swap-sym weight 0.75 solution swap-sym weight 0.0
g1 0.999 seed 0 E -3.657403 first-solve dim 16 E_loc -0.03674 start swap-sym weight 0.75 solution swap-sym weight 1.0
g1 0.999 seed 1 E -3.657403 first-solve dim 16 E_loc 0.01479 start swap-sym weight 0.75 solution swap-sym weight 0.0
g1 0.999 seed 2 E -3.657403 first-solve dim 16 E_loc 0.00969 start swap-sym weight 0.75 solution swap-sym weight 0.0
g1 0.999 seed 3 E -3.657403 first-solve dim 16 E_loc 0.00455 start swap-sym weight 0.75 solution swap-sym weight 1.0
g1 0.999 seed 4 E -3.657403 first-solve dim 16 E_loc 0.01797 start swap-sym weight 0.75 solution swap-sym weight 0.0
g1 0.999 seed 5 E -3.657403 first-solve dim 16 E_loc 0.06049 start swap-sym weight 0.75 solution swap-sym weight 0.0
```

The solution's weight is exactly 1 or 0, never in between, and decides the
final energy. Narrowed from the hunter's "the random start cannot leave the
sector": the start has the other sector's weight, and the solver discards
it. The failure does not need the boson hopping (hop 0 gives -1.5 for an
exact -2.41246).

**Suggested fix**: a noise term on `"python"` that follows the sweep
schedule.

After the fix:

`sb_after.py`:

```python
import numpy as np
from dmrgpy import bosonchain
def run(v, nb=4, hop=1.0, maxm=20, nsweeps=10):
    b = bosonchain.SpinBoson_Chain(["S=1/2"]*2+["B"]*nb, itensor_version=v)
    h = b.Sx[0]*b.Sx[1]+b.Sy[0]*b.Sy[1]+b.Sz[0]*b.Sz[1]+(b.Sx[0]+b.Sx[1])*(b.A[2]+b.Adag[2])
    for i in range(2, 2+nb): h = h+(b.N[i]-1.0)*(b.N[i]-1.0)
    for i in range(2, 1+nb): h = h+hop*b.Adag[i]*b.A[i+1]
    h = h+h.get_dagger(); b.set_hamiltonian(h); b.maxm = maxm; b.nsweeps = nsweeps
    return b.gs_energy()
print("exact -5.917370 (nb 4, hop 1), -2.412460 (nb 4, hop 0)")
print("python nb4 hop1", [round(float(run("python")), 6) for _ in range(12)])
print("python nb4 hop0", [round(float(run("python", hop=0.0)), 6) for _ in range(6)])
print("v3     nb4 hop1", [round(float(run(3)), 6) for _ in range(2)])
```

```
exact -5.917370 (nb 4, hop 1), -2.412460 (nb 4, hop 0)
python nb4 hop1 [-5.91737, -5.91737, -5.91737, -5.91737, -5.91737, -5.91737, -5.91737, -5.91737, -5.91737, -5.91737, -5.91737, -5.91737]
python nb4 hop0 [-2.412455, -2.412455, -2.412455, -2.412455, -2.412455, -2.412455]
v3     nb4 hop1 [-5.91737, -5.91737]
```
