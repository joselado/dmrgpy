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
