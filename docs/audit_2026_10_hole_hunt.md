# Audit, 2026-10: five-lens hole hunt over 8d12cb6 and 5a861ba

This record is the seventh hole hunt over the Python layer. It is scoped to two
commits: `8d12cb6`, which added `submode="TASK"`, the tangent-space Krylov
dynamical correlator of arXiv:2510.07279, and made it the default Kondo route
on v3 and `"python"`, and `5a861ba`, which moved `get_kondo_spectrum` into
`kondospectrumtk/kondomethod.py` and put it on the two spinful fermion chains.
It ran on 2026-10-23 on `5a861ba` (clean tree apart from three untracked
example PNGs, both compiled extensions current: `make -n pybind` has nothing
to do in `mpscpp3`) as five hunters, one per lens, and one reviewer per
candidate briefed to refute it. Every repro was executed. The hunters returned
four candidates with a measured wrong number or a lost safeguard and a handful
of loud or cosmetic items; the reviewers refuted none, narrowed every one, and
struck several sub-claims, which stay visible below. Two lenses, the v3 and
`"python"` state export and the TASK dispatch, found no defect, and what they
ruled out is under "Ruled out" so that the next hunt does not repeat it.

No compiled snapshot of a parent commit was built. Both features are new in
these two commits, so a defect in them is new by construction; the two
findings that are older than the commits (finding 4's crash, and the
underlying v2 limitation) were checked against `Spin_Chain`, which had the
same route before, and say so.

This file is the evidence, not a task list. Each entry gains a `**Status**`
line when it is acted on, and keeps its repro.

All seven findings were fixed on 2026-10-23 in the same session, by hand
rather than by parallel agents, since the three clusters (`task`:
`taskdmrg.py`, `dynamics.py`; `kondo`: `potentialdc.py`, `taskkondo.py`,
`kondomethod.py`; `docs`: the example, the guide, `fermionchain.py`) are
small. The regressions are `tests/test_audit_2026_10_task.py` and
`tests/test_audit_2026_10_kondo.py`; 17 of their 19 tests fail on
`5a861ba`, and the other two are controls. No fix changes a number that
was returned before, so `docs/numerical_changes.md` only records that.

## The five lenses

| Lens | Brief |
|---|---|
| `task-core` | The numerics of `pyitensor/task.py` on `"python"`: sum rules, gauge, rank-deficient and degenerate states, complex input, Lanczos thresholds, off-diagonal pairs, the two-site variance, against ED |
| `task-v3-export` | How v3 and `"python"` hand their state to the TASK code: unit scale of `H_`, index order, QN blocks, conserved sectors, injected states |
| `task-dispatch` | `dynamics.py` and `taskdmrg.py` around `submode="TASK"`: keyword consumers, non-Hermitian H, `mode=` pins, short chains, the Lehmann-density convention |
| `kondo-task` | `get_kondo_spectrum(mode="DMRG")` through TASK: `taskkondo.py`, the snapping of near-zero energies, `n_gs`, the ED fallback, the two-site impurity spin of the interleaved class |
| `kondo-fermion` | `5a861ba` outside TASK: the impurity spin on both fermion classes, the KPM/EX/two-time routes on fermion chains, sectors on the ED route, the Hubbard-to-Heisenberg example |

## Scope

Out of scope by construction: vendored ITensor and TDVP; the three legacy
quirks of `.claude/rules/cpp-backend.md`; the open `docs/known_issue_*.md`
items; anything in an earlier audit record; what `ROADMAP.md` marks as absent.
`itensor_version="julia_live"` is out of scope, since TASK does not exist
there and its JIT cost would dominate.

Every script ran through a wrapper that caps the hunt at three concurrent
Python processes, pins threads and forces this checkout's `src`:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 <script>
```

Scripts and outputs are spliced in from disk, not retyped. `<scratch>` stands
for the session's scratch folder and `<repo>` for the checkout, the only edit
made to either. Reviewer paragraphs are condensed from their hand-back text.

## The findings at a glance

| # | Finding | Severity | Origin |
|---|---|---|---|
| 1 | TASK loses the weight of a fermion operator on a parity-mixed state at truncated `maxm`, up to 100%, with a small variance | MEDIUM | `8d12cb6` |
| 2 | The potential term's `es` coverage warning never fires on the fermion chains | LOW | `5a861ba` |
| 3 | The TASK Kondo route has a free detector of a wrong ground state, a negative Ritz energy, and does not read it | LOW | `8d12cb6` |
| 4 | `order=3` with an explicit correlator submode crashes on a short chain and on v2, and a comment promises the v2 route | LOW | older |
| 5 | The Hubbard-to-Heisenberg example draws a $(t/U)^2$ line through points that are not yet asymptotic | LOW | `5a861ba` |
| 6 | `nkry` is cut to an integer on the TASK correlator and not checked on the Kondo ED fallback; a scalar `es` gives an opaque IndexError | LOW | `8d12cb6` |
| 7 | Two docstrings are stale: `dynamics.py` calls TASK `"python"`-only, the native class calls v2 unavailable | LOW | both |

## Findings

### 1. TASK loses up to 100% of the spectral weight of `C_i`/`Cdag_i` when the ground state mixes fermion parities and `maxm` is truncated, and the two-site variance stays small

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED (narrowed) &middot; lens `task-core`

**Status**: FIXED (diagnostic) on 2026-10-23. `taskdmrg.get_poles` returns
`missing_weight`, the sum rule `<A B>` of the same state minus the elastic
term minus the pole weights, and raises a RuntimeWarning when it is above
1e-6 of the weight; the docstring and `docs/guide/correlator-methods.md`
say what the variance does not see. The weights themselves are unchanged,
so no number changes. Pinned by `tests/test_audit_2026_10_task.py`
(the odd-term chain of `13_parity.py` on `"python"` and 3, 0.107 missing
of a sum rule of 0.426 at site 9, warns; a parity-definite Kitaev state at
`maxm=2` and a one-site spin operator at `maxm=4` stay below 1e-8; a bond
operator warns). The regression does not use the topological Kitaev chain:
which mixture of its degenerate pair DMRG reaches depends on the random
start, and in the full suite v3 once reached one that kept the weight, so
the test uses a ground state that mixes the parities by itself.

**Where**: `src/dmrgpy/pyitensor/task.py`, `TangentSpace.project_operator`
(the start vector `P B|GS>`) and `poles`; `two_site_variance` is the error
bar the docs promise; `taskdmrg.get_poles` returns no other check. Reached on
`itensor_version="python"` and 3 through `get_dynamical_correlator(submode="TASK")`.

TASK starts its Lanczos from `P B|GS>`, the projection of `B|GS>` onto the
one-site tangent space. The docstring says this is exact for a sum of
one-site operators. A fermion `C_i = F_0...F_{i-1} a_i` is a product over
`i+1` sites. On a state of definite parity the left reduced density matrix commutes
with the parity of the left block, so the string maps the kept left space
into itself and the projection is still exact; that is why every
number-conserving test passes. On a parity mixture this commutation breaks,
and the part of `B|GS>` outside the tangent space is dropped. The variance measures the part
of `H|alpha>` outside the two-site space, not this, so it does not see it.
The guide (`docs/guide/correlator-methods.md`) calls the variance "an error
bar for each pole", and nothing mentions Jordan-Wigner strings or parity.

A parity mixture is a legitimate ground state of the topological Kitaev
chain, whose even and odd ground states are split by far less than the DMRG
resolution; DMRG returns an arbitrary superposition (`<P>` about -2.6e-5
below), with E0 within 2.4e-7 of ED. An explicit odd term in H gives the
same loss.

**Expected**: the weights plus the elastic term sum to `<GS|A B|GS>` of the
same MPS, or the call reports how much is missing.

Repro, hunter (`14_kitaev.py`, with `common.py`):

```python
import numpy as np
from dmrgpy import spinchain, fermionchain
from dmrgpy import taskdmrg
def heis(n=6, hz=0.3, maxm=64, J=1.0, extra=None, cutoff=0., dm=0., hy=0.0, fieldsite=0):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
        if dm: h = h + dm*(sc.Sx[i]*sc.Sy[i+1]-sc.Sy[i]*sc.Sx[i+1])
    if hz: h = h + hz*sc.Sz[fieldsite]
    if hy: h = h + hy*sc.Sy[fieldsite]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = maxm, 20, cutoff
    return sc
def ed_lehmann(chain, A, B):
    ed = chain.get_ED_obj()
    h = np.array(ed.get_hamiltonian().todense())
    e, U = np.linalg.eigh(h)
    Ae = U.conj().T @ np.array(ed.MO2matrix(A).todense()) @ U
    Be = U.conj().T @ np.array(ed.MO2matrix(B).todense()) @ U
    return e, U, Ae, Be
def poles(chain, A, B, nkry=100, variance=True):
    from dmrgpy.groundstate import ground_state_on_session
    ground_state_on_session(chain)
    from dmrgpy.multioperatortk import canonical
    Adag = None if canonical.is_dagger_pair(A,B) else A.get_dagger().to_terms()
    ws,wt,el,var,dim = taskdmrg.session_poles(chain,A.to_terms(),B.to_terms(),Adag,nkry,variance=variance)
    return ws,wt,el,var,dim
def report(chain, A, B, nkry=100, tag=""):
    ws,wt,el,var,dim = poles(chain,A,B,nkry)
    e,U,Ae,Be = ed_lehmann(chain,A,B)
    M = Ae[0,:]*Be[:,0]
    tot_ed = M.sum()  # <A B>
    print(tag,"dim",dim,"npoles",len(ws))
    print("  sum task + elastic:",wt.sum()+el," ED <AB>:",tot_ed," elastic",el, "ED elastic",M[0])
    # spectral first moment <A (H-E0) B>
    m1_task = (wt*ws).sum(); m1_ed = (M*(e-e[0])).sum()
    print("  first moment task",m1_task," ED",m1_ed)
    print("  max variance",None if var is None else np.max(var))
    return ws,wt,el,var
```

```python
from common import *
n=10
fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
h=0
t,D,mu=1.0,0.8,0.5
for i in range(n-1):
    h = h - t*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]) + D*(fc.C[i]*fc.C[i+1]+fc.Cdag[i+1]*fc.Cdag[i])
for i in range(n): h = h - mu*fc.N[i]
fc.set_hamiltonian(h)
ed=fc.get_ED_obj(); H=np.array(ed.get_hamiltonian().todense()); w=np.linalg.eigvalsh(H); print("ED lowest",w[:3])
for maxm in (2,4,8,16):
  fc.maxm, fc.nsweeps, fc.cutoff = maxm, 40, 0.
  fc.clean() if hasattr(fc,"clean") else None
  e0=fc.gs_energy()
  par=fc.vev(sum((1-2*fc.N[i]) for i in range(1))*1) if False else None
  print("maxm",maxm,"E0",e0)
  for i in (0,5,9):
    ws,wt,el,var,dim=poles(fc,fc.Cdag[i],fc.C[i],100)
    nvev=fc.vev(fc.Cdag[i]*fc.C[i]).real
    t_=(wt.sum()+el).real
    print("  site",i,"task %.6f vev %.6f missing %.1f%%  elastic %.4g (<Cdag><C>) weighted var %.3g"%(t_,nvev,100*(nvev-t_)/nvev,el.real,(abs(wt)*var).sum()))
```

Observed:

```
ED lowest [-10.82630556 -10.82630509  -9.35658415]
maxm 2 E0 -10.826294769767687
  site 0 task 0.630483 vev 0.630483 missing 0.0%  elastic 0.2315 (<Cdag><C>) weighted var 0.0234
  site 5 task 0.007391 vev 0.566340 missing 98.7%  elastic 2.659e-08 (<Cdag><C>) weighted var 0.00227
  site 9 task 0.000000 vev 0.630477 missing 100.0%  elastic 1.849e-11 (<Cdag><C>) weighted var 2.01e-08
maxm 4 E0 -10.826305320961078
  site 0 task 0.630476 vev 0.630476 missing 0.0%  elastic 0.2315 (<Cdag><C>) weighted var 0.0635
  site 5 task 0.475535 vev 0.566338 missing 16.0%  elastic 1.353e-06 (<Cdag><C>) weighted var 0.123
  site 9 task 0.016006 vev 0.630476 missing 97.5%  elastic 5.967e-10 (<Cdag><C>) weighted var 0.00486
maxm 8 E0 -10.826305325700543
  site 0 task 0.630476 vev 0.630476 missing 0.0%  elastic 0.2315 (<Cdag><C>) weighted var 0.000443
  site 5 task 0.566338 vev 0.566338 missing 0.0%  elastic 1.351e-06 (<Cdag><C>) weighted var 0.000523
  site 9 task 0.630247 vev 0.630476 missing 0.0%  elastic 6.045e-10 (<Cdag><C>) weighted var 0.000464
maxm 16 E0 -10.826305325709342
  site 0 task 0.630476 vev 0.630476 missing 0.0%  elastic 0.2315 (<Cdag><C>) weighted var 3.67e-10
  site 5 task 0.566338 vev 0.566338 missing 0.0%  elastic 1.351e-06 (<Cdag><C>) weighted var 4.78e-10
  site 9 task 0.630476 vev 0.630476 missing 0.0%  elastic 6.042e-10 (<Cdag><C>) weighted var 4.17e-10
```

The same through the public API, by the reviewer, on both backends
(`01.py`):

```python
import numpy as np, sys
from dmrgpy import fermionchain
def build(it,mu,D=0.8,t=1.0,n=10,maxm=4,seed=None):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=it)
    h=0
    for i in range(n-1):
        h = h - t*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]) + D*(fc.C[i]*fc.C[i+1]+fc.Cdag[i+1]*fc.Cdag[i])
    for i in range(n): h = h - mu*fc.N[i]
    fc.set_hamiltonian(h); fc.maxm,fc.nsweeps,fc.cutoff=maxm,40,0.
    return fc
def run(it,mu,maxm):
    fc=build(it,mu,maxm=maxm)
    e=fc.gs_energy(); ed=fc.gs_energy(mode="ED")
    P=1
    for i in range(10): P=P*(1-2*fc.N[i])
    par=fc.vev(P).real
    print("== it",it,"mu",mu,"maxm",maxm,"E0",e,"ED",ed,"parity<P>",par)
    es=np.linspace(-1,14,3000)
    for i in (0,5,9):
        A,B=fc.Cdag[i],fc.C[i]
        try:
            _,y,p=fc.get_dynamical_correlator(submode="TASK",name=(A,B),es=es,delta=0.05,return_poles=True)
        except Exception as ex:
            print(" err",repr(ex)[:200]); return
        _,yed=fc.get_dynamical_correlator(mode="ED",submode="ED",name=(A,B),es=es,delta=0.05)
        print(" site",i,"TASK %.5f ED %.5f vev %.5f var*w %.3g"%(np.trapezoid(y.real,es),np.trapezoid(np.asarray(yed).real,es),fc.vev(A*B).real,(abs(p["weight"])*p["variance"]).sum()))
for it,mu,mm in [("python",0.5,4),(3,0.5,4),(2,0.5,4),("python",3.0,2),("python",3.0,4),(3,3.0,2)]:
    run(it,mu,mm)
```

```
== it python mu 0.5 maxm 4 E0 -10.826305320961072 ED -10.826305560725821 parity<P> -2.626815496188846e-05
 site 0 TASK 0.62378 ED 0.62378 vev 0.63048 var*w 0.0635
 site 5 TASK 0.47258 ED 0.56260 vev 0.56634 var*w 0.123
 site 9 TASK 0.01590 ED 0.62378 vev 0.63048 var*w 0.00486
== it 3 mu 0.5 maxm 4 E0 -10.826305320958912 ED -10.826305560725821 parity<P> -2.5446809169811335e-05
 site 0 TASK 0.62378 ED 0.62378 vev 0.63048 var*w 0.0632
 site 5 TASK 0.47255 ED 0.56260 vev 0.56634 var*w 0.123
 site 9 TASK 0.01578 ED 0.62378 vev 0.63048 var*w 0.00483
== it 2 mu 0.5 maxm 4 E0 -10.82630532095894 ED -10.826305560725821 parity<P> -2.544680908129448e-05
 err NotImplementedError('get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got itensor_version=2)')
== it python mu 3.0 maxm 2 E0 -31.012791366621812 ED -31.013153354298773 parity<P> 0.9999999996842441
 site 0 TASK 0.97431 ED 0.97426 vev 0.97981 var*w 0.000362
 site 5 TASK 0.95117 ED 0.95110 vev 0.95669 var*w 0.000415
 site 9 TASK 0.97428 ED 0.97426 vev 0.97978 var*w 0.000362
== it python mu 3.0 maxm 4 E0 -31.01315327302097 ED -31.013153354298773 parity<P> 0.9999999996842441
 site 0 TASK 0.97426 ED 0.97426 vev 0.97976 var*w 1.78e-07
 site 5 TASK 0.95110 ED 0.95110 vev 0.95662 var*w 1.96e-07
 site 9 TASK 0.97426 ED 0.97426 vev 0.97976 var*w 1.78e-07
== it 3 mu 3.0 maxm 2 E0 -31.012791363224835 ED -31.013153354298773 parity<P> 1.00000000000033
 site 0 TASK 0.97431 ED 0.97426 vev 0.97981 var*w 0.000362
 site 5 TASK 0.95117 ED 0.95110 vev 0.95669 var*w 0.000415
 site 9 TASK 0.97428 ED 0.97426 vev 0.97978 var*w 0.000362
```

And on 14 sites (`02.py`, the parity column is a placeholder 0, removed
after an n=24 attempt ran out of memory):

```python
import numpy as np
from dmrgpy import fermionchain
n=14
for mu,maxm in [(0.5,4),(0.5,8)]:
    fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
    h=0
    for i in range(n-1):
        h = h - (fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]) + 0.8*(fc.C[i]*fc.C[i+1]+fc.Cdag[i+1]*fc.Cdag[i])
    for i in range(n): h = h - mu*fc.N[i]
    fc.set_hamiltonian(h); fc.maxm,fc.nsweeps,fc.cutoff=maxm,40,0.
    e=fc.gs_energy()
    P=1
    break_=0
    pass
    print("== n",n,"mu",mu,"maxm",maxm,"E0",e,"<P>",0)
    es=np.linspace(-1,16,3400)
    for i in (0,7,13):
        A,B=fc.Cdag[i],fc.C[i]
        _,y,p=fc.get_dynamical_correlator(submode="TASK",name=(A,B),es=es,delta=0.05,return_poles=True)
        print(" site",i,"TASK %.5f vev %.5f var*w %.3g"%(np.trapezoid(y.real,es),fc.vev(A*B).real,(abs(p["weight"])*p["variance"]).sum()))
```

```
== n 14 mu 0.5 maxm 4 E0 -15.50348459324436 <P> 0
 site 0 TASK 0.62389 vev 0.63048 var*w 0.0543
 site 7 TASK 0.42371 vev 0.56634 var*w 0.0753
 site 13 TASK 0.00000 vev 0.63048 var*w 8.8e-08
== n 14 mu 0.5 maxm 8 E0 -15.50348460119361 <P> 0
 site 0 TASK 0.62389 vev 0.63048 var*w 0.00199
 site 7 TASK 0.56268 vev 0.56634 var*w 0.00241
 site 13 TASK 0.62263 vev 0.63048 var*w 0.00789
```

The number-conserving control loses nothing (`12_jw_sumrule.py`), and an odd
term on an otherwise number-conserving chain loses up to 25% (`13_parity.py`,
not rerun by the reviewer):

```python
from common import *
n=10
fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
h=0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] + 1.0*fc.N[i]*fc.N[i+1]
h = h + 0.3*fc.N[0]
fc.set_hamiltonian(h)
fc.maxm, fc.nsweeps, fc.cutoff = 4, 40, 0.
print("E0",fc.gs_energy(), "ED",fc.gs_energy(mode="ED"))
print("site  sum_rule_task  <Cdag C>=<N> (vev)   missing   frac   maxvar   | hole channel: sum, <C Cdag>=1-<N>, missing")
for i in range(n):
    ws,wt,el,var,dim=poles(fc,fc.Cdag[i],fc.C[i],100)
    nvev=fc.vev(fc.N[i]).real
    ws2,wt2,el2,var2,dim2=poles(fc,fc.C[i],fc.Cdag[i],100)
    t=(wt.sum()+el).real; t2=(wt2.sum()+el2).real
    print(i,"%.6f %.6f  %.6f %.1f%%  %.3g | %.6f %.6f %.6f %.1f%%"%(t,nvev,nvev-t,100*(nvev-t)/nvev,var.max(),t2,1-nvev,1-nvev-t2,100*(1-nvev-t2)/(1-nvev)))
# spin comparison: one-site op
```

```
E0 -5.188583075816922 ED -5.198205910561557
site  sum_rule_task  <Cdag C>=<N> (vev)   missing   frac   maxvar   | hole channel: sum, <C Cdag>=1-<N>, missing
0 0.301048 0.301048  0.000000 0.0%  0.987 | 0.698952 0.698952 0.000000 0.0%
1 0.526016 0.526016  0.000000 0.0%  0.986 | 0.473984 0.473984 0.000000 0.0%
2 0.285077 0.285077  0.000000 0.0%  0.986 | 0.714923 0.714923 0.000000 0.0%
3 0.463689 0.463689  0.000000 0.0%  0.986 | 0.536311 0.536311 0.000000 0.0%
4 0.390192 0.390192  0.000000 0.0%  0.987 | 0.609808 0.609808 0.000000 0.0%
5 0.347208 0.347208  0.000000 0.0%  0.981 | 0.652792 0.652792 0.000000 0.0%
6 0.493675 0.493675  0.000000 0.0%  0.985 | 0.506325 0.506325 0.000000 0.0%
7 0.295713 0.295713  0.000000 0.0%  0.979 | 0.704287 0.704287 0.000000 0.0%
8 0.474593 0.474593  0.000000 0.0%  0.983 | 0.525407 0.525407 -0.000000 -0.0%
9 0.422789 0.422789  0.000000 0.0%  0.985 | 0.577211 0.577211 -0.000000 -0.0%
```

```python
from common import *
n=10
fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
h=0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] + 1.0*fc.N[i]*fc.N[i+1]
h = h + 0.3*fc.N[0] + 0.5*(fc.Cdag[4]+fc.C[4])
fc.set_hamiltonian(h)
for maxm in (4,):
  fc.maxm, fc.nsweeps, fc.cutoff = maxm, 40, 0.
  print("maxm",maxm,"E0",fc.gs_energy(), "ED",fc.gs_energy(mode="ED"))
  print("site  task  vev<Cdag C>   missing  frac  maxvar(Ritz, <1 relative) weighted-var | hole")
  for i in (0,2,5,7,9):
    ws,wt,el,var,dim=poles(fc,fc.Cdag[i],fc.C[i],100)
    nvev=fc.vev(fc.Cdag[i]*fc.C[i]).real
    ws2,wt2,el2,var2,dim2=poles(fc,fc.C[i],fc.Cdag[i],100)
    nv2=fc.vev(fc.C[i]*fc.Cdag[i]).real
    t=(wt.sum()+el).real; t2=(wt2.sum()+el2).real
    print(i,"%.6f %.6f %.6f %.1f%% var wtd %.3g | %.6f %.6f %.1f%%"%(t,nvev,nvev-t,100*(nvev-t)/nvev,(abs(wt)*var).sum(),t2,nv2,100*(nv2-t2)/nv2))
```

```
maxm 4 E0 -5.374245579143566 ED -5.380600875062999
site  task  vev<Cdag C>   missing  frac  maxvar(Ritz, <1 relative) weighted-var | hole
0 0.323252 0.323252 0.000000 0.0% var wtd 0.0204 | 0.676748 0.676748 0.0%
2 0.323481 0.323481 0.000000 0.0% var wtd 0.0349 | 0.676519 0.676519 0.0%
5 0.330095 0.335800 0.005706 1.7% var wtd 0.0302 | 0.650850 0.664200 2.0%
7 0.238861 0.294880 0.056019 19.0% var wtd 0.0179 | 0.636659 0.705120 9.7%
9 0.318625 0.426006 0.107381 25.2% var wtd 0.0183 | 0.382113 0.573994 33.4%
```

**Reviewer (CONFIRMED, narrowed)**: reproduced through the public API on
`"python"` and v3 (site 9: 0.0159 and 0.0158 against the sum rule 0.6305).
The anchor is the `vev` sum rule of the same MPS; the ED integral is a
parity eigenstate and agrees with the sum rule only at site 0, so ED is not
a clean anchor for the lost fraction. Narrowed: "0% at `maxm>=8`" holds for
n=10 only, since at n=14 `maxm=8` still loses 1.2% at the last site; the
mechanism fits the parity-definite control (mu=3, `maxm=2`, no loss) and the
hunter's dense contraction, but the reviewer did not inspect `P B|GS>`
directly; the variance is not miscomputed, it is the documented error bar
that fails to cover the loss. Severity medium, leaning high: the output looks
like a normal spectrum, and Majorana and superconducting chains are what a
`C`/`Cdag` user runs.

**Suggested fix**: return the missing weight, `vev(A*B) - elastic -
sum(weights)`, in the `get_poles` dict and warn when it is not small; it
also catches the documented loss of a bond operator (3.4% for `Sx4*Sx5` at
`maxm=4`). For a fermion operator on a state whose parity expectation is
not +-1, warn or refuse and ask for a ground state projected onto one parity
sector, which for the degenerate Kitaev pair is the physical answer anyway.
The component of `B|GS>` outside the tangent space cannot be added as a start
vector, since the recursion has no representation for it.

### 2. The potential term's coverage warning for a too-narrow `es` never fires on the spinful fermion chains, so the term comes back up to 99% wrong in silence

`bug` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `kondo-fermion`

**Status**: FIXED on 2026-10-23. `_spin_casimir` keeps S(S+1) on a spin
site and measures `sum_k <S_k S_k>` with `chain.vev` on any other chain, so
the coverage warning fires on both fermion classes; the warning still says
S(S+1) on a spin site. No returned number changes. Pinned by
`tests/test_audit_2026_10_kondo.py` (the measured sum rule against the ED
spectrum, below 3/4 on a Hubbard site; a narrow `es` warns, a wide one does
not).

**Where**: `src/dmrgpy/kondospectrumtk/potentialdc.py:85-93`
(`_spin_casimir` returns None unless the chain is a `Spin_Chain`) and the
skip of the warning when it is None (about line 142); reached from
`kondomethod.py:307-311` with `mode="DMRG"`, an explicit correlator submode,
`U!=0` and `order=3`.

The potential-interference term needs `es` to reach past the top of the
`S_k` spectrum, and `potentialdc` checks this with a sum rule against
S(S+1). On a Hubbard site there is no fixed S(S+1), so the check returns
None and is skipped. Before `5a861ba` only `Spin_Chain` reached this route.
The default TASK route does not read `es`, which is why only an explicit
`submode="KPM"` (or another correlator submode) is exposed.

**Expected**: the same RuntimeWarning as on a spin chain.

Repro, hunter (`03_potential_sumrule.py`, the internal function):

```python
import sys; sys.stdout.reconfigure(line_buffering=True)
import numpy as np, warnings
from dmrgpy import fermionchain, spinchain
from dmrgpy.kondospectrumtk import conductance, potentialdc
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
def hubbard(n,U,b,iv=3):
    fc=fermionchain.Spinful_Fermionic_Chain_Native(n,itensor_version=iv)
    h=0
    for i in range(n-1):
        h=h-(fc.Cdagup[i]*fc.Cup[i+1]+fc.Cdagdn[i]*fc.Cdn[i+1])
    h=h+h.get_dagger()
    for i in range(n): h=h+U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)+b*fc.Sz[i]
    fc.set_hamiltonian(h); fc.maxm=64; fc.nsweeps=20; fc.cutoff=0.
    return fc
fc=hubbard(3,2.,0.3)
eV=np.linspace(-1.2,1.2,13)
ks=KondoSpectrum(fc,1,0.)
print("top of spectrum from GS:",ks.e[-1])
ref=conductance.third_order_potential_dIdV(ks,eV,-0.05,0.3,T0=1.,omega0=4.,Gamma0=1e-2)
for top in (0.8,12.):
    es=np.linspace(-0.3,top,800)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        d=potentialdc.third_order_potential_dIdV_dc(fc,1,eV,-0.05,0.3,T0=1.,omega0=4.,Gamma0=1e-2,submode="KPM",delta=0.05,es=es)
    print("es top",top,"max|ref|",np.abs(ref).max(),"max|diff|",np.abs(d-ref).max(),"nwarnings",len(w))
# same on a Spin_Chain for contrast: does the sum rule warn there
sc=spinchain.Spin_Chain(["S=1/2"]*3)
h=0
for i in range(2): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
h=h+0.3*sum(sc.Sz); sc.set_hamiltonian(h)
es=np.linspace(-0.3,0.6,400)
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    potentialdc.third_order_potential_dIdV_dc(sc,1,eV,-0.05,0.3,T0=1.,omega0=4.,Gamma0=1e-2,submode="KPM",delta=0.05,es=es)
print("spin chain narrow es: nwarnings",len(w))
# weight above es top from ED
w_tot=0; w_hi=0
occ=np.where(ks.p>0)[0]
Xi=np.stack([ks.Sx,ks.Sy,ks.Sz],-1)
for i in occ:
    wts=ks.p[i]*np.sum(np.abs(Xi[:,i,:])**2,axis=1)
    en=ks.e-ks.e[i]
    w_tot+=wts.sum(); w_hi+=wts[en>0.8+0.1].sum()
print("fraction of sum_k S_kk weight above es top 0.8: ",w_hi/w_tot)
```

Observed:

```
top of spectrum from GS: 6.790178748749579
es top 0.8 max|ref| 0.11267238289417131 max|diff| 0.11173705507734331 nwarnings 0
es top 12.0 max|ref| 0.11267238289417131 max|diff| 0.00033555620674763365 nwarnings 0
spin chain narrow es: nwarnings 1
fraction of sum_k S_kk weight above es top 0.8:  0.9620831522254859
```

Reviewer, through the public `get_kondo_spectrum` on both classes
(`01_public.py`):

```python
import sys; sys.stdout.reconfigure(line_buffering=True)
import numpy as np, warnings
from dmrgpy import fermionchain
def hubbard(cls,n,U,b,iv=3):
    fc=cls(n,itensor_version=iv)
    h=0
    for i in range(n-1):
        h=h-(fc.Cdagup[i]*fc.Cup[i+1]+fc.Cdagdn[i]*fc.Cdn[i+1])
    h=h+h.get_dagger()
    for i in range(n): h=h+U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)+b*fc.Sz[i]
    fc.set_hamiltonian(h); fc.maxm=64; fc.nsweeps=20; fc.cutoff=0.
    return fc
eV=np.linspace(-1.2,1.2,13)
kw=dict(site=1,Jrho_s=-0.05,U=0.3,T=0.,order=3,omega0=4.,Gamma0=1e-2,T0=1.)
for name,cls,iv in (("native",fermionchain.Spinful_Fermionic_Chain_Native,3),("interleaved",fermionchain.Spinful_Fermionic_Chain,3)):
    fc=hubbard(cls,3,2.,0.3,iv)
    try:
        x,ref=fc.get_kondo_spectrum(eV,mode="ED",**kw)
    except Exception as e:
        print(name,"ED failed",repr(e)); continue
    print(name,"ED max|dIdV|",np.abs(ref).max())
    for top in (0.8,12.):
        es=np.linspace(-0.3,top,800)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                x,d=fc.get_kondo_spectrum(eV,mode="DMRG",submode="KPM",delta=0.05,es=es,dt2=1.,n_t2_half=20,dtau=1.,n_tau_half=20,**kw)
            except Exception as e:
                print(name,top,"fail",repr(e)); continue
        print(name,"es top",top,"max|diff| full",np.abs(d-ref).max(),"warnings",[str(i.message)[:60] for i in w])
```

```
native ED max|dIdV| 2.516470213431383
native es top 0.8 max|diff| full 0.19167316950141355 warnings []
native es top 12.0 max|diff| full 0.13192128127507408 warnings []
interleaved ED max|dIdV| 2.516470213431383
interleaved es top 0.8 max|diff| full 0.19167316950059998 warnings []
interleaved es top 12.0 max|diff| full 0.13192128127635527 warnings []
```

**Reviewer (CONFIRMED, narrowed)**: no warning on either class through the
public route. The 99% is the potential term alone. In the full dI/dV, where
the Kondo term dominates (peak 2.52), the narrow grid costs roughly 0.06 to
0.11, about 3 to 4.5%; the wide-grid residual of 0.132 comes from the
reviewer's coarse time grids, and the two grids were not run against each
other directly. The general need for a wide `es` is documented (guide
`kondo.md`, the `kondomethod.py` docstring, `.claude/rules/kondo.md`), but
not that fermion chains get no warning. Severity low to medium.

**Suggested fix**: compare against the measured sum rule
`sum_a <GS|S_a S_a|GS>` (a `vev` of `Sx^2+Sy^2+Sz^2` at the site), not a
constant: on a Hubbard site it is below 3/4 by the empty and doubly occupied
weight. Reword the warning, which says S(S+1). No returned number changes.

### 3. The TASK Kondo route hands a wrong ground state to the sums without a word, although a negative Ritz energy proves it is wrong

`diagnostic` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `kondo-task`

**Status**: FIXED on 2026-10-23. `TaskKondoSpectrum` raises a
RuntimeWarning when a Ritz energy survives below the snapping window. No
number changes. Pinned by `tests/test_audit_2026_10_kondo.py` (a state set
with `set_gs` above the ground state warns, the solved ground state does
not).

**Where**: `src/dmrgpy/kondospectrumtk/taskkondo.py`, `TaskKondoSpectrum.__init__`
(the energies `e` are used as returned).

On a 5-site Heisenberg chain with a field `B=1e-8` splitting the Kramers
doublet, DMRG converges to a mixture of the two members; `gs_energy()` is
itself off by about the splitting. The TASK spectrum then contains a Ritz
energy below the state it started from, `e[1]=-4.3e-9`, and the Kondo curve
at `n_gs=1` misses ED by an order-one amount. A negative Ritz energy means
the tangent space reaches a state below the "ground state", which is an
exact proof that it is not one.

Repro, hunter (`05_snap.py`, and `06_snap2.py`, which `exec`s the builders
of `01_sumrule.py`, not reproduced):

```python
import numpy as np
from dmrgpy import spinchain
# odd chain, Kramers doublet split by a field B on all sites; ED vs TASK n_gs=1, 2
def chain(n,B,big):
    sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
    h=0
    for i in range(n-1):
        h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    for i in range(n): h=h+B*sc.Sz[i]
    if big: h=h+ 1e-3*sc.Sz[0]*sc.Sz[1]*0 # no-op
    sc.set_hamiltonian(h); sc.maxm=64; sc.nsweeps=30; sc.cutoff=0.
    return sc
eV=np.linspace(-1e-4,1e-4,41)
kw=dict(site=0,Jrho_s=-0.1,U=0.,T=0.,omega0=20e-3,Gamma0=5e-8,order=3)
for B in [0,1e-12,1e-11,1e-10,1e-8,1e-6]:
    c=chain(5,B,0)
    _,ed=c.get_kondo_spectrum(eV,mode="ED",**kw)
    for ng in [1,2]:
        _,d=c.get_kondo_spectrum(eV,mode="DMRG",nkry=50,n_gs=ng,**kw)
        print("B %.0e n_gs %d dev %.3e"%(B,ng,np.max(np.abs(d-ed))/np.max(np.abs(ed))))
```

```
B 0e+00 n_gs 1 dev 1.687e-09
B 0e+00 n_gs 2 dev 2.806e-10
B 1e-12 n_gs 1 dev 1.243e-09
B 1e-12 n_gs 2 dev 2.527e-11
B 1e-11 n_gs 1 dev 1.523e-08
B 1e-11 n_gs 2 dev 8.102e-01
B 1e-10 n_gs 1 dev 2.663e-07
B 1e-10 n_gs 2 dev 8.102e-01
B 1e-08 n_gs 1 dev 1.621e+00
B 1e-08 n_gs 2 dev 8.099e-01
B 1e-06 n_gs 1 dev 1.187e-13
B 1e-06 n_gs 2 dev 7.607e-01
```

```python
import numpy as np
exec(open('05_snap.py').read().split("eV=np.linspace")[0])
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
for B in [1e-9,1e-8,1e-7,1e-6]:
    c=chain(5,B,0)
    ks=TaskKondoSpectrum(c,0,nkry=50)
    ek=KondoSpectrum(c,0,0.)
    Sz=sum(c.Sz[i] for i in range(5))
    print("B",B,"TASK e[:4]",ks.e[:4],"ED e[:4]",ek.e[:4],"ED p",ek.p[:3])
    print("  <Stot z> dmrg %.6f"%np.real(c.vev(Sz)),"gs E err %.2e"%(c.gs_energy()-c.gs_energy(mode="ED")))
    print("  TASK <Sz_0> %.6f"%np.real(ks.Sz[0,0]),"ED %.6f"%np.real(ek.Sz[0,0]))
```

```
B 1e-09 TASK e[:4] [0.00000000e+00 5.41275513e-10 7.20779243e-04 7.20780243e-04] ED e[:4] [0.00000000e+00 9.99999991e-10 7.20779472e-04 7.20780472e-04] ED p [1. 0. 0.]
  <Stot z> dmrg -0.270638 gs E err 2.29e-10
  TASK <Sz_0> -0.138476 ED -0.255833
B 1e-08 TASK e[:4] [ 0.00000000e+00 -4.34309990e-09  7.20772300e-04  7.20782301e-04] ED e[:4] [0.00000000e+00 9.99999999e-09 7.20779472e-04 7.20789472e-04] ED p [1. 0. 0.]
  <Stot z> dmrg 0.217150 gs E err 7.17e-09
  TASK <Sz_0> 0.111105 ED -0.255833
B 1e-07 TASK e[:4] [0.00000000e+00 1.00000000e-07 7.20779472e-04 7.20879472e-04] ED e[:4] [0.00000000e+00 1.00000000e-07 7.20779472e-04 7.20879472e-04] ED p [1. 0. 0.]
  <Stot z> dmrg -0.500000 gs E err -3.04e-18
  TASK <Sz_0> -0.255833 ED -0.255833
B 1e-06 TASK e[:4] [0.00000000e+00 1.00000000e-06 7.20779472e-04 7.21779472e-04] ED e[:4] [0.00000000e+00 1.00000000e-06 7.20779472e-04 7.21779472e-04] ED p [1. 0. 0.]
  <Stot z> dmrg -0.500000 gs E err 6.51e-19
  TASK <Sz_0> -0.255833 ED -0.255833
```

Reviewer (`r2_routes.py`):

```python
import numpy as np
from dmrgpy import spinchain
def chain(n,B,sweeps=30):
    sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
    h=0
    for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    for i in range(n): h=h+B*sc.Sz[i]
    sc.set_hamiltonian(h); sc.maxm=64; sc.nsweeps=sweeps; sc.cutoff=0.
    return sc
eV=np.linspace(-1e-4,1e-4,41)
kw=dict(site=0,Jrho_s=-0.1,U=0.,T=0.,omega0=20e-3,Gamma0=5e-8,order=2)
def dev(a,b): return np.max(np.abs(a-b))/np.max(np.abs(b))
for B in [1e-8,1e-6]:
  for sw in [30]:
    c=chain(5,B,sw)
    _,ed=c.get_kondo_spectrum(eV,mode="ED",**kw)
    Sz=sum(c.Sz[i] for i in range(5))
    print("B %.0e sweeps %d gsE err %.2e <Sz> %.4f"%(B,sw,c.gs_energy()-c.gs_energy(mode="ED"),np.real(c.vev(Sz))))
    for sm in ["TASK","EX","KPM"]:
        try:
            ex=dict(nkry=50) if sm=="TASK" else dict(delta=2e-6,es=np.linspace(-2e-2,2e-2,2001))
            _,d=c.get_kondo_spectrum(eV,mode="DMRG",submode=sm,**ex,**kw)
            print("   ",sm,"dev vs ED %.3e"%dev(d,ed))
        except Exception as e: print("   ",sm,type(e).__name__,str(e)[:100])
print("gamma0 in meV-units: J=1e-3, Gamma0=5e-8 ; B=1e-8 vs Gamma0 ratio",1e-8/5e-8)
```

```
B 1e-08 sweeps 30 gsE err 4.69e-09 <Sz> -0.0315
    TASK dev vs ED 8.315e-11
    EX dev vs ED 2.957e+00
    KPM dev vs ED 6.486e+00
B 1e-06 sweeps 30 gsE err 4.34e-19 <Sz> -0.5000
    TASK dev vs ED 6.754e-10
    EX dev vs ED 2.528e+00
    KPM dev vs ED 5.180e+00
gamma0 in meV-units: J=1e-3, Gamma0=5e-8 ; B=1e-8 vs Gamma0 ratio 0.2
```

**Reviewer (CONFIRMED, narrowed)**: the wrong state is DMRG
non-convergence at an accidental near-degeneracy, which every route inherits
equally and which the `n_gs` caveat of the docstring documents. The 162% is
one random start: the same script gives 1.5e-8 at `B=1e-11` and 2.7e-7 at
`1e-10`, and an `order=2` run gave 133%. `B/Gamma0=0.2`, so the regime is
below the kernel width. The framing "a TASK bug" is struck. What is specific
to TASK is the free detector, `min(e) < -tol`, which catches the mixed
doublet; it does not catch the second case the hunter raised, a `"python"`
Nf-sector Hubbard chain off by 2.5e-4 from a ground-state energy error of
2.4e-10 (there `min(e)=0`), whose amplification was not compared against
the other routes.

**Suggested fix**: warn in `TaskKondoSpectrum` when a Ritz energy lies below
`-degeneracy_tol*width`. No number changes.

### 4. `get_kondo_spectrum(order=3)` with an explicit correlator submode raises `AttributeError` on a chain below three sites and on v2, and the comment in `kondomethod.py` promises that v2 runs that route

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `kondo-fermion`

**Status**: FIXED on 2026-10-23. On a chain too short for the session the
Kondo term of the correlator route is the exact one of `conductance.py` on
the T=0 ED spectrum, as the correlators there already were; `order=3` with
a correlator submode raises NotImplementedError on any backend other than
`"python"` and 3; the comment is corrected. A call that crashed now returns
a number, and no returned number changes. Pinned by
`tests/test_audit_2026_10_kondo.py` (order 3 minus order 2 on a 2-site v3
chain is the ED Kondo term to 1e-10; v2 raises). The wiring test
`test_kondo_spectrum.py::test_get_kondo_spectrum_dmrg_combines_terms_correctly`
runs on a one-site chain and now tells it the session is there.

**Where**: `kondospectrumtk/dmrgtwotime.py` (calls `chain._session.tdvp_step`
after the correlator fell back to ED for ns<3); `kondomethod.py:244-248`
(the comment "v2 and julia_live have no TASK, so they keep the correlator and
two-time route"; v2 has no `tdvp_step`).

Loud, with no wrong number. A 2-site `Spin_Chain` on v3 fails the same way,
so the crash is older than `5a861ba`; that commit makes the Hubbard dimer, the
example's own case, a natural way to reach it. The default TASK route on v3
and `"python"` falls back to ED for a short chain correctly.

Repro, hunter (`05_twotime_short.py`, `10_v2.py`):

```python
import sys; sys.stdout.reconfigure(line_buffering=True)
import numpy as np, warnings
from dmrgpy import fermionchain
def hubbard(cls,n,U,b,iv):
    fc=cls(n,itensor_version=iv)
    h=0
    for i in range(n-1):
        h=h-(fc.Cdagup[i]*fc.Cup[i+1]+fc.Cdagdn[i]*fc.Cdn[i+1])
    h=h+h.get_dagger()
    for i in range(n): h=h+U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)+b*fc.Sz[i]
    fc.set_hamiltonian(h); fc.maxm=64; fc.nsweeps=20; fc.cutoff=0.
    return fc
eV=np.linspace(-1.5,1.5,7)
es=np.linspace(-0.3,6,400)
tw=dict(dt2=0.2,n_t2_half=20,dtau=0.2,n_tau_half=20)
KW=dict(site=0,Jrho_s=-0.05,kB=1.,omega0=4.,Gamma0=0.2,order=3,T=0.)
for cls,n,iv in ((fermionchain.Spinful_Fermionic_Chain_Native,2,3),(fermionchain.Spinful_Fermionic_Chain_Native,2,"python"),
                 (fermionchain.Spinful_Fermionic_Chain_Native,3,3),(fermionchain.Spinful_Fermionic_Chain,2,3),
                 (fermionchain.Spinful_Fermionic_Chain_Native,3,2)):
    fc=hubbard(cls,n,4.,0.3,iv)
    ed=fc.get_kondo_spectrum(eV,**KW)[1]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            d=fc.get_kondo_spectrum(eV,mode="DMRG",submode="KPM",delta=0.05,es=es,**tw,**KW)[1]
        print(cls.__name__,n,iv,"maxdiff",np.abs(d-ed).max(),"scale",np.abs(ed).max())
    except BaseException as e:
        print(cls.__name__,n,iv,"ERR",type(e).__name__,str(e)[:200])
```

```
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
Spinful_Fermionic_Chain_Native 2 3 ERR AttributeError 'Fermionic_State' object has no attribute 'cpp_handle'
Spinful_Fermionic_Chain_Native 2 python maxdiff 0.4996508568666149 scale 4.676159083217802
Spinful_Fermionic_Chain_Native 3 3 maxdiff 0.11897284576786471 scale 4.592844482187626
Spinful_Fermionic_Chain 2 3 maxdiff 0.49965085686623245 scale 4.676159083217802
Spinful_Fermionic_Chain_Native 3 2 ERR AttributeError 'dmrgpy.mpscpp2._dmrgcpp.Chain' object has no attribute 'tdvp_step'
```

```python
import sys; sys.stdout.reconfigure(line_buffering=True)
import numpy as np, warnings
from dmrgpy import fermionchain, spinchain
exec(open("07_ngs.py").read().split("eV=")[0].split("import numpy as np, warnings")[1])
eV=np.linspace(-1.5,1.5,13); es=np.linspace(-0.3,8,500)
tw=dict(dt2=0.3,n_t2_half=5,dtau=0.3,n_tau_half=5)
for cls in (fermionchain.Spinful_Fermionic_Chain,fermionchain.Spinful_Fermionic_Chain_Native):
  for order in (2,3):
    try:
        fc=hubbard(cls,3,4.,0.3,2)
        ed=fc.get_kondo_spectrum(eV,site=1,T=0.,kB=1.,order=order,Jrho_s=-.05,omega0=4.,Gamma0=.2)[1]
        kw=dict(tw) if order==3 else {}
        d=fc.get_kondo_spectrum(eV,site=1,T=0.,kB=1.,order=order,Jrho_s=-.05,omega0=4.,Gamma0=.2,mode="DMRG",delta=0.05,es=es,**kw)[1]
        print(cls.__name__,"v2 order",order,"diff",np.abs(d-ed).max(),ed.max())
    except BaseException as e:
        print(cls.__name__,"v2 order",order,"ERR",type(e).__name__,str(e)[:160])
# spin chain v2 order 3 for contrast
sc=spinchain.Spin_Chain(["S=1/2"]*3,itensor_version=2)
h=0
for i in range(2): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h+0.3*sum(sc.Sz))
try:
    sc.get_kondo_spectrum(eV,site=1,T=0.,kB=1.,order=3,Jrho_s=-.05,omega0=4.,Gamma0=.2,mode="DMRG",delta=0.05,es=es,**tw); print("spin v2 order3 ok")
except BaseException as e: print("Spin_Chain v2 order3 ERR",type(e).__name__,str(e)[:160])
```

```
Spinful_Fermionic_Chain v2 order 2 diff 0.38804149656505116 3.0295943695357903
Spinful_Fermionic_Chain v2 order 3 ERR AttributeError 'dmrgpy.mpscpp2._dmrgcpp.Chain' object has no attribute 'tdvp_step'
Spinful_Fermionic_Chain_Native v2 order 2 diff 0.3880414965725505 3.0295943695357903
Spinful_Fermionic_Chain_Native v2 order 3 ERR AttributeError 'dmrgpy.mpscpp2._dmrgcpp.Chain' object has no attribute 'tdvp_step'
Spin_Chain v2 order3 ERR AttributeError 'dmrgpy.mpscpp2._dmrgcpp.Chain' object has no attribute 'tdvp_step'
```

Reviewer (`02_spin2_v2native.py`):

```python
import sys; sys.stdout.reconfigure(line_buffering=True)
import numpy as np, warnings
from dmrgpy import spinchain, fermionchain
eV=np.linspace(-1.5,1.5,7); es=np.linspace(-0.3,8,200)
tw=dict(dt2=0.3,n_t2_half=3,dtau=0.3,n_tau_half=3)
sc=spinchain.Spin_Chain(["S=1/2"]*2,itensor_version=3)
h=sc.Sx[0]*sc.Sx[1]+sc.Sy[0]*sc.Sy[1]+sc.Sz[0]*sc.Sz[1]
sc.set_hamiltonian(h)
try:
    sc.get_kondo_spectrum(eV,site=0,T=0.,kB=1.,order=3,Jrho_s=-.05,omega0=4.,Gamma0=.2,mode="DMRG",submode="KPM",delta=0.05,es=es,**tw); print("Spin_Chain n=2 v3 ok")
except BaseException as e: print("Spin_Chain n=2 v3 ERR",type(e).__name__,str(e)[:160])
# v2 native: what does DMRG do
fc=fermionchain.Spinful_Fermionic_Chain_Native(3,itensor_version=2)
h=0
for i in range(2): h=h-(fc.Cdagup[i]*fc.Cup[i+1]+fc.Cdagdn[i]*fc.Cdn[i+1])
h=h+h.get_dagger()
for i in range(3): h=h+4*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)+0.3*fc.Sz[i]
fc.set_hamiltonian(h)
print("sites",[type(s).__name__ for s in fc.sites][:3] if hasattr(fc,'sites') else None, getattr(fc,'site_names',None), getattr(fc,'sites',None) and fc.sites[:3])
fc.maxm=4; fc.nsweeps=10; fc.cutoff=1e-12
print("DMRG gs maxm=4",fc.gs_energy(mode="DMRG"),"ED",fc.gs_energy(mode="ED"))
print("Sz0 vev DMRG",fc.vev(fc.Sz[0]*fc.Sz[1],mode="DMRG"),"ED",fc.vev(fc.Sz[0]*fc.Sz[1],mode="ED"))
```

```
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
Spin_Chain n=2 v3 ERR AttributeError 'State' object has no attribute 'cpp_handle'
sites ['int', 'int', 'int'] None [1, 1, 1]
DMRG gs maxm=4 -4.386067977499788 ED -4.386067977499783
Sz0 vev DMRG (-0.13090169944006+0j) ED (-0.1309016994374949+0j)
```

**Reviewer (CONFIRMED)**: the 2-site `Spin_Chain` raises the same
`AttributeError`; v2 `order=3` raises on both fermion classes and on
`Spin_Chain`; only `order=2` runs on v2.

**Suggested fix**: route the two-time Kondo term to ED when the correlator
did, and raise a clear `NotImplementedError` for `order=3` on v2; correct
the comment.

### 5. The Hubbard-to-Heisenberg example draws a $(t/U)^2$ guide line through points at U/t=4 to 32 that fall by 2.3 to 2.9 per doubling, not 4

`docs` &middot; severity **LOW** &middot; CONFIRMED (narrowed) &middot; lens `kondo-fermion`

**Status**: FIXED on 2026-10-23. The example sweeps the deviation up to
U/t=256 and anchors the guide line there, and its header and
`docs/guide/kondo.md` say the deviation falls as $(t/U)^2$ once U/t is a
few tens. The example was run; its printed deviations reach the ratio 4
per doubling at the largest U/t. No test pins the sweep.

**Where**: `examples/kondo/kondo_third_order_hubbard_VS_heisenberg/main.py:16-18`
and its guide line (about line 118, anchored at the last point);
`docs/guide/kondo.md:192`; the `5a861ba` commit message.

The docs say the difference is "of order $(t/U)^2$", which is an asymptotic
statement and is true: the reviewer extended the sweep to U/t=256, where the
ratio per doubling is 4.00 at both fields. At the plotted U/t the points are
not yet asymptotic, and the line drawn through the last one misses the first
three. The ratio 6.66 at U/t=32 is grid quantisation of step positions.

Reviewer (`01_asym.py`, the example's own chain extended):

```python
import sys; sys.stdout.reconfigure(line_buffering=True)
import numpy as np
src=open("<repo>/examples/kondo/kondo_third_order_hubbard_VS_heisenberg/main.py").read()
pre=src.split("x = np.linspace(-3.")[0].replace("import matplotlib.pyplot as plt","").replace("import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')","")
exec(pre)
n=6
for N in (240,2400):
  x=np.linspace(-3.,3.,N)
  for b in (0.,1.):
    hs=spectrum(heisenberg(n,1.,b=b),x,1.)
    prev=None
    for Ut in (4.,8.,16.,32.,64.,128.,256.):
        J=4*t**2/Ut
        d=np.mean(np.abs(spectrum(hubbard(n,Ut*t,b=b*J),x,J)-hs))
        print("grid",N,"b",b,"U/t",Ut,"dev %.6g"%d,"ratio",None if prev is None else round(prev/d,2)); prev=d
```

```
grid 240 b 0.0 U/t 4.0 dev 0.116849 ratio None
grid 240 b 0.0 U/t 8.0 dev 0.0517791 ratio 2.26
grid 240 b 0.0 U/t 16.0 dev 0.0176509 ratio 2.93
grid 240 b 0.0 U/t 32.0 dev 0.00264979 ratio 6.66
grid 240 b 0.0 U/t 64.0 dev 0.000775398 ratio 3.42
grid 240 b 0.0 U/t 128.0 dev 0.00016204 ratio 4.79
grid 240 b 0.0 U/t 256.0 dev 4.05303e-05 ratio 4.0
grid 240 b 1.0 U/t 4.0 dev 0.0590835 ratio None
grid 240 b 1.0 U/t 8.0 dev 0.0244679 ratio 2.41
grid 240 b 1.0 U/t 16.0 dev 0.00929073 ratio 2.63
grid 240 b 1.0 U/t 32.0 dev 0.00210128 ratio 4.42
grid 240 b 1.0 U/t 64.0 dev 0.000485023 ratio 4.33
grid 240 b 1.0 U/t 128.0 dev 0.000121467 ratio 3.99
grid 240 b 1.0 U/t 256.0 dev 3.03798e-05 ratio 4.0
../run3.sh: line 7: 2134240 Killed                  MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 "$@"
```

**Reviewer (CONFIRMED, narrowed)**: the physics claim stands, the plot over
the plotted range does not. The dimer test
(`tests/test_kondo_spectrum_spinful.py`, U=20 and 40) is a real pin; no test
covers the chain sweep.

**Suggested fix**: extend the sweep to U/t of 64 to 256, or anchor the line
at the largest U/t and say "for U/t of a few tens and above".

### 6. `nkry` is cut to an integer on the TASK correlator and not checked at all on the Kondo ED fallback, and a scalar `es` fails with an opaque IndexError

`bug` &middot; severity **LOW** &middot; lenses `task-dispatch`, `kondo-task`

**Status**: FIXED on 2026-10-23. `nkry` must be a positive integer in
`taskdmrg.get_poles` and, ahead of the short-chain branch, in
`kondomethod.py`; `es=None` and a 2-D `es` raise ValueError and a scalar
`es` is one frequency. Numbers change only for a caller who passed a
non-integer `nkry`, which now raises. Pinned by
`tests/test_audit_2026_10_task.py` and `tests/test_audit_2026_10_kondo.py`.

**Where**: `taskdmrg.py`, `get_poles` (`int(nkry)` after a `nkry<1` check)
and `dynamical_correlator` (`es[:,None]`); `kondomethod.py:256` pops `nkry`,
which only `taskkondo.py:37` validates, and the short-chain branch never
builds a `TaskKondoSpectrum`.

`nkry=2.7` on the correlator silently runs a 2-step Krylov space; `nkry=-5`,
`"x"` or `2.5` on a short-chain Kondo call return the ED curve, and raise on
a longer chain. On the fallback the value is inert and the ED curve exact, so
no number is wrong there. `es=None` or a scalar `es` raise
`IndexError: too many indices`.

Repro (`01_kw.py`, task-dispatch hunter):

```python
import numpy as np, traceback
from dmrgpy import spinchain, fermionchain
def mk(v, n=6):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h=0
    for i in range(n-1): h = h + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); return sc
es=np.linspace(-1,5,7)
def tryit(label,f):
    try:
        r=f(); y=r[1] if isinstance(r,tuple) else r
        print(label,"OK",np.round(np.asarray(y)[:4],4))
    except Exception as e:
        print(label,"RAISE",type(e).__name__,str(e)[:100])
for v in ["python",3,2]:
    sc=mk(v)
    print("== version",v)
    pair=(sc.Sz[1],sc.Sz[1])
    tryit("base",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=0.2))
    tryit("kernel=lorentz",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=0.2,kernel="lorentz"))
    tryit("bogus",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=0.2,bogus=3))
    tryit("n=50 (kpm kw)",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=0.2,n=50))
    tryit("nex",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=0.2,nex=3))
    tryit("string ZZ",lambda: sc.get_dynamical_correlator(submode="TASK",name="ZZ",i=1,j=1,es=es,delta=0.2))
    tryit("mode=ED",lambda: sc.get_dynamical_correlator(submode="TASK",mode="ED",name=pair,es=es,delta=0.2))
    tryit("delta=0",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=0))
    tryit("nkry=0",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=.2,nkry=0))
    tryit("nkry=2.7",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=.2,nkry=2.7))
    tryit("noname",lambda: sc.get_dynamical_correlator(submode="TASK",es=es,delta=.2))
    tryit("wf",lambda: sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=.2,wf=None))
sc=mk("python"); sc.mode="ED"
tryit("pinned ED",lambda: sc.get_dynamical_correlator(submode="TASK",name=(sc.Sz[1],sc.Sz[1]),es=es,delta=0.2))
sc=mk(3,2)
tryit("v3 ns=2",lambda: sc.get_dynamical_correlator(submode="TASK",name=(sc.Sz[1],sc.Sz[1]),es=es,delta=0.2))
sc=mk("python",1)
tryit("py ns=1",lambda: sc.get_dynamical_correlator(submode="TASK",name=(sc.Sz[0],sc.Sz[0]),es=es,delta=0.2))
```

```
== version python
base OK [0.0037+0.j 0.0199+0.j 0.132 +0.j 0.0906+0.j]
kernel=lorentz RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'kernel'
bogus RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'bogus'
n=50 (kpm kw) RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'n'
nex RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'nex'
string ZZ OK [0.0037+0.j 0.0199+0.j 0.132 +0.j 0.0906+0.j]
mode=ED RAISE NotImplementedError get_dynamical_correlator: submode='TASK' has no ED implementation
delta=0 RAISE ValueError delta must be positive, got 0
nkry=0 RAISE ValueError nkry must be a positive integer, got 0
nkry=2.7 OK [0.0037+0.j 0.0193+0.j 0.0622+0.j 0.0503+0.j]
noname RAISE ValueError name must be a documented correlator string (with i=/j=) or a pair of operators, got None
wf RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'wf'
== version 3
base OK [0.0037+0.j 0.0199+0.j 0.132 +0.j 0.0906+0.j]
kernel=lorentz RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'kernel'
bogus RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'bogus'
n=50 (kpm kw) RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'n'
nex RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'nex'
string ZZ OK [0.0037+0.j 0.0199+0.j 0.132 +0.j 0.0906+0.j]
mode=ED RAISE NotImplementedError get_dynamical_correlator: submode='TASK' has no ED implementation
delta=0 RAISE ValueError delta must be positive, got 0
nkry=0 RAISE ValueError nkry must be a positive integer, got 0
nkry=2.7 OK [0.0037+0.j 0.0193+0.j 0.0622+0.j 0.0503+0.j]
noname RAISE ValueError name must be a documented correlator string (with i=/j=) or a pair of operators, got None
wf RAISE TypeError dynamical_correlator() got an unexpected keyword argument 'wf'
== version 2
base RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
kernel=lorentz RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
bogus RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
n=50 (kpm kw) RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
nex RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
string ZZ RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
mode=ED RAISE NotImplementedError get_dynamical_correlator: submode='TASK' has no ED implementation
delta=0 RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
nkry=0 RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
nkry=2.7 RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
noname RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
wf RAISE NotImplementedError get_dynamical_correlator: submode="TASK" is implemented on itensor_version="python" and 3 only (got 
pinned ED RAISE NotImplementedError get_dynamical_correlator: submode='TASK' has no ED implementation
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
v3 ns=2 RAISE NotImplementedError get_dynamical_correlator: submode='TASK' has no ED implementation
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
py ns=1 RAISE TypeError the ED spin chain cannot build an operator from 0 (int): expected a MultiOperator, an EDOperator or 
```

The `es` crash and the non-Hermitian refusal (`05_misc.py`, task-dispatch
hunter):

```python
import numpy as np
from dmrgpy import spinchain, fermionchain
es=np.linspace(-1,6,29); d=0.15
def tryit(label,f):
    try:
        r=f(); print(label,"OK",type(r).__name__,[np.shape(a) for a in r] if isinstance(r,tuple) else "")
        return r
    except Exception as e:
        print(label,"RAISE",type(e).__name__,str(e)[:140])
for v in ["python",3]:
    print("==",v)
    # complex hopping
    n=4
    fc = fermionchain.Fermionic_Chain(n, itensor_version=v)
    rng=np.random.RandomState(3)
    h=0
    for i in range(n-1):
        t=np.exp(1j*rng.rand()*3)*(0.5+rng.rand())
        h = h + t*fc.Cdag[i]*fc.C[i+1]+np.conj(t)*fc.Cdag[i+1]*fc.C[i]
    h = h + 0.4*fc.N[0]*fc.N[1]+0.3*fc.N[0]-0.2*fc.N[3]
    fc.set_hamiltonian(h); fc.maxm=64; fc.nsweeps=20
    for lab,pair in [("Cdag0 C2",(fc.Cdag[0],fc.C[2])),("C1 Cdag3",(fc.C[1],fc.Cdag[3])),("Cdag1 C1",(fc.Cdag[1],fc.C[1]))]:
        ed = fc.get_dynamical_correlator(mode="ED",submode="ED",name=pair,es=es,delta=d)[1]
        t = fc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=d,nkry=100)[1]
        print("complex hop",lab,"max|ED|=%.3f dev=%.2e maxIm ED %.3f"%(abs(ed).max(),abs(t-ed).max(),abs(ed.imag).max()))
    # non-Hermitian
    n=4
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h=0
    for i in range(n-1): h = h + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.3j*sc.Sz[0]
    sc.set_hamiltonian(h)
    tryit("nonherm TASK",lambda: sc.get_dynamical_correlator(submode="TASK",name=(sc.Sz[1],sc.Sz[1]),es=es,delta=d))
    # es=None
    sc = spinchain.Spin_Chain([2]*4, itensor_version=v)
    h=0
    for i in range(3): h = h + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    tryit("es=None",lambda: sc.get_dynamical_correlator(submode="TASK",name=(sc.Sz[1],sc.Sz[1]),es=None,delta=d))
    tryit("es scalar",lambda: sc.get_dynamical_correlator(submode="TASK",name=(sc.Sz[1],sc.Sz[1]),es=1.0,delta=d))
    tryit("es 2d",lambda: sc.get_dynamical_correlator(submode="TASK",name=(sc.Sz[1],sc.Sz[1]),es=np.ones((2,2)),delta=d))
    tryit("return_poles via ED",lambda: sc.get_dynamical_correlator(mode="ED",submode="ED",name=(sc.Sz[1],sc.Sz[1]),es=es,delta=d,return_poles=True))
    tryit("DMRG return_poles w/ variance=False",lambda: sc.get_dynamical_correlator(submode="TASK",name=(sc.Sz[1],sc.Sz[1]),es=es,delta=d,return_poles=True,variance=False))
    # degenerate gs: n=5 Heisenberg
    sc = spinchain.Spin_Chain([2]*5, itensor_version=v)
    h=0
    for i in range(4): h = h + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm=64; sc.nsweeps=20
    for lab,pair in [("Sz0Sz0",(sc.Sz[0],sc.Sz[0])),("Sx0Sx0",(sc.Sx[0],sc.Sx[0]))]:
        ed = sc.get_dynamical_correlator(mode="ED",submode="ED",name=pair,es=es,delta=d)[1]
        t = sc.get_dynamical_correlator(submode="TASK",name=pair,es=es,delta=d,nkry=100)[1]
        print("degenerate",lab,"max|ED|=%.3f dev=%.2e"%(abs(ed).max(),abs(t-ed).max()), "<Sz_0>=%.3f"%sc.vev(sc.Sz[0]).real)
```

```
== python
complex hop Cdag0 C2 max|ED|=0.296 dev=1.53e-13 maxIm ED 0.171
complex hop C1 Cdag3 max|ED|=0.244 dev=1.46e-13 maxIm ED 0.097
complex hop Cdag1 C1 max|ED|=0.532 dev=1.49e-13 maxIm ED 0.000
nonherm TASK RAISE NotImplementedError get_dynamical_correlator: submode='TASK' assumes a Hermitian Hamiltonian and has no non-Hermitian implementation (KPM/CVM/CVM_explicit do, a
es=None RAISE IndexError too many indices for array: array is 0-dimensional, but 1 were indexed
es scalar RAISE IndexError too many indices for array: array is 0-dimensional, but 1 were indexed
es 2d RAISE ValueError operands could not be broadcast together with shapes (2,1,2) (1,4) 
return_poles via ED RAISE TypeError dynamical_correlator_ED() got an unexpected keyword argument 'return_poles'
DMRG return_poles w/ variance=False OK tuple [(29,), (29,), ()]
degenerate Sz0Sz0 max|ED|=0.185 dev=2.66e-15 <Sz_0>=0.246
degenerate Sx0Sx0 max|ED|=0.185 dev=8.60e-16 <Sz_0>=0.246
== 3
complex hop Cdag0 C2 max|ED|=0.296 dev=2.46e-14 maxIm ED 0.171
complex hop C1 Cdag3 max|ED|=0.244 dev=2.87e-14 maxIm ED 0.097
complex hop Cdag1 C1 max|ED|=0.532 dev=4.88e-14 maxIm ED 0.000
nonherm TASK RAISE NotImplementedError get_dynamical_correlator: submode='TASK' assumes a Hermitian Hamiltonian and has no non-Hermitian implementation (KPM/CVM/CVM_explicit do, a
es=None RAISE IndexError too many indices for array: array is 0-dimensional, but 1 were indexed
es scalar RAISE IndexError too many indices for array: array is 0-dimensional, but 1 were indexed
es 2d RAISE ValueError operands could not be broadcast together with shapes (2,1,2) (1,4) 
return_poles via ED RAISE TypeError dynamical_correlator_ED() got an unexpected keyword argument 'return_poles'
DMRG return_poles w/ variance=False OK tuple [(29,), (29,), ()]
degenerate Sz0Sz0 max|ED|=0.185 dev=1.46e-13 <Sz_0>=0.210
degenerate Sx0Sx0 max|ED|=0.185 dev=7.41e-14 <Sz_0>=0.210
```

Repro (`07_short.py`, kondo-task hunter):

```python
import numpy as np
from dmrgpy import spinchain
def mk(n,ver):
    sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=ver)
    h=0
    for i in range(n-1): h=h+1e-3*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h); return sc
eV=np.linspace(-3e-3,3e-3,21)
kw=dict(site=0,Jrho_s=-0.1,U=0.,T=0.,omega0=20e-3,Gamma0=5e-5,order=3)
for ver,n in [("python",1),(3,2),("python",3)]:
    c=mk(n,ver)
    for extra in [dict(nkry=-5),dict(nkry="x"),dict(nkry=2.5),dict(n_gs=2)]:
        try:
            _,d=c.get_kondo_spectrum(eV,mode="DMRG",**kw,**extra)
            print(ver,n,extra,"OK max %.4g"%np.max(np.abs(d)))
        except Exception as e: print(ver,n,extra,type(e).__name__,str(e)[:80])
```

```
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
python 1 {'nkry': -5} TypeError the ED spin chain cannot build an operator from 0 (int): expected a MultiOperato
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
python 1 {'nkry': 'x'} TypeError the ED spin chain cannot build an operator from 0 (int): expected a MultiOperato
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
python 1 {'nkry': 2.5} TypeError the ED spin chain cannot build an operator from 0 (int): expected a MultiOperato
pyitensor's two-site DMRG can't handle a chain this short (n=1 < 2 sites), using default ED routines
python 1 {'n_gs': 2} NotImplementedError n_gs>1 sets each member of the ground-state manifold as the ground state of this
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
3 2 {'nkry': -5} OK max 8.226
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
3 2 {'nkry': 'x'} OK max 8.226
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
3 2 {'nkry': 2.5} OK max 8.226
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
3 2 {'n_gs': 2} NotImplementedError n_gs>1 sets each member of the ground-state manifold as the ground state of this
python 3 {'nkry': -5} ValueError nkry must be a positive integer, got -5
python 3 {'nkry': 'x'} ValueError nkry must be a positive integer, got 'x'
python 3 {'nkry': 2.5} ValueError nkry must be a positive integer, got 2.5
python 3 {'n_gs': 2} OK max 7.363
```

**Reviewer (CONFIRMED, narrowed)** for the Kondo half: reproduced exactly
(`r1_short.py`), and narrowed to "a typo is not reported on a short chain",
with no numerical consequence. The correlator's `int(nkry)` and the `es`
crash were checked by reading `taskdmrg.py`, not by a separate reviewer.

**Suggested fix**: one positive-integer check of `nkry` in `get_poles` and
ahead of the `get_mode` branch in `kondomethod.py`, next to the `n_gs`
check; `np.atleast_1d(es)` with a 1-D check.

### 7. Two docstrings are stale: `dynamics.py` says TASK is `"python"`-only, and `Spinful_Fermionic_Chain_Native` says v2 cannot run DMRG on it

`docs` &middot; severity **LOW** &middot; CONFIRMED &middot; lenses `task-dispatch`, `kondo-fermion`

**Status**: FIXED on 2026-10-23, both docstrings corrected. No test.

**Where**: `src/dmrgpy/dynamics.py:97-99` (TASK is on `"python"` and 3,
`taskdmrg.BACKENDS`); `src/dmrgpy/fermionchain.py:258-262` ("mpscpp2/get_sites.h
never registered a Hubbard/Electron site type"; it does, at lines 43, 80 and
106).

**Reviewer (CONFIRMED)** for the native class: v2 runs genuine DMRG on the
native site, `gs_energy` -4.386067977499788 against ED -4.386067977499783 on
a 3-site chain (in finding 4's `02_spin2_v2native.out`); the sub-claim of
possible wrong numbers on v2 is struck. The `dynamics.py` line was checked
by reading.

**Suggested fix**: correct both docstrings.

## Ruled out

Each with the lens and the anchor, so the next hunt does not repeat it.

- TASK at saturation, every pair tested (one-site, `Sp`/`Sm`, bond, `S(q)`, `Sy`, off-diagonal, complex H, fermion pairs with complex hopping, S=1, mixed, boson, both spinful classes): ED to 1e-12 or better, elastic term included identically (`task-core`, `task-dispatch`).
- Fermion `C`/`Cdag` on a number-conserving chain at truncated `maxm`: sum rule exact at every `maxm` (`task-core`).
- `mixed_gauge` on 300 random, zeroed and oversized bonds, real and complex: isometry and reconstruction hold; product-state and dimerized ground states exact (`task-core`, `task-v3-export`).
- `two_site_variance`, `apply` and `block_krylov` against dense constructions: 3e-15, 2e-15, 6e-14 (`task-core`).
- Unit scale: `hamiltonian_mpo` returns `H_` in physical units; v3 agrees with ED from J=1 down to J=1e-12; `"python"` at J<=1e-9 inherits the known Lanczos floor of `pyitensor/dmrg.py` (`docs/audit_2026_09_25_open_items.md`) (`task-v3-export`).
- Conserved sectors: conserving operators match restricted ED; charge-changing ones raise a clear ValueError; after `promote_to_dense()` they match full ED. On `"python"` the charge penalty is not in `self.H`; its tangent space also lists zero-weight out-of-sector poles at negative omega (weight 9e-25), which v3 does not (`task-v3-export`).
- Injected states (`set_gs` of an excited state), a Hamiltonian changed after a solve, a ground state solved by ED first: TASK measures what the chain holds (`task-v3-export`, `task-dispatch`).
- Non-Hermitian H, `mode="ED"`, an ED pin, v2, unknown keywords: all raise (`task-dispatch`).
- The interleaved class's two-site impurity spin: documented in the `kondomethod.py` docstring (0.33 against 0.16 relative error at `maxm=4`, equal at `maxm=8`) (`kondo-task`).
- Snapping of near-zero Ritz energies: only splittings below 1e-9 of the bandwidth, and ED snaps the same (`kondo-task`).
- The third-order Kondo term is about 1.3 times less accurate than the second-order one at truncated `maxm`, consistent with the docstring (`kondo-task`).
- `Sx/Sy/Sz[site]` on both fermion classes is exactly c^dagger sigma c/2 of orbital `site`; the non-TASK routes address it through the MultiOperator, never a session index (two-time grid DMRG against ED to 5e-15 and 1e-12); ED Nf sectors exact at T=0 and canonical at T>0 as documented; an Sz sector raises; the move of `Spin_Chain`'s method into the mixin is behaviour-identical (`kondo-fermion`).

## New leads, not reviewed

- `"python"` with an Nf sector amplifies a ground-state energy error of 2.4e-10 into a 2.5e-4 error of the saturated TASK Kondo spectrum (finding 3's reviewer); not compared against the other routes.
- The default `nkry=100` is not exhausted on larger tangent spaces and nothing says so (4.1e-5 on a 10-site chain at `maxm=40`, 5e-13 at `nkry=300`): a risk, not a measured failure.
- v3's fully polarized sector gives an empty tangent space; the public `get_dynamical_correlator` on it was not run.
- v2 `order=2` KPM on a 3-site Hubbard chain differs from ED by 0.39 on a peak of 3.03, the same on both classes; probably broadening at `delta=0.05`, not checked.
