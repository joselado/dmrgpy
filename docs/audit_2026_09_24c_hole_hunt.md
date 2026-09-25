# Audit, 2026-09-24 (third pass): four-lens hole hunt of the fix pass `867e2b4`

Findings from a four-lens automated audit of the Python layer, run 2026-09-24 on
`867e2b4` (clean tree, both compiled extensions current, `make -n pybind`
reporting nothing to do in `mpscpp2` and `mpscpp3`), scoped to one commit,
`867e2b4`, the fix pass that closed the eighteen findings of
`audit_2026_09_24b_hole_hunt.md`. The reason for scoping a hunt to a single fix
commit is the measurement of the two records before this one: seven or eight of
the sixteen holes of `audit_2026_09_24_hole_hunt.md` came in with `765b537`, the
commit that closed the hunt before it, and the second pass then found eighteen
more reached from `30200a4`. Each lens took one of `867e2b4`'s four fix clusters.
Every finding below carries a repro that was actually executed and its verbatim
output, and every one was then handed to an independent reviewer whose brief was
to *refute* it; seven were turned up by a reviewer while reviewing a neighbouring
candidate and got a reviewer of their own before entering the record. Of the 22
candidates reviewed, one was refuted, on scope: `factor>1` in
`_fourier_transform_correlator` builds its oversampled grid with the wrong
spacing, which is real and reproduces to every digit, but is the 2026-09-24
record's own "New lead" and was left untouched by `867e2b4`, so it is not
reproduced here. Three pairs of the surviving 21 turned out to be one defect
reached from two sides, and each pair is one entry below carrying both reviews,
so the record has 18 findings. Fourteen of the 21 came back narrowed, and each
entry keeps its struck sub-claims visible next to what survived.

This file is the evidence, not a task list, the same convention as the four
earlier records: it records what was observed and how to reproduce it, so that a
fix, or a decision that the behaviour is intended after all, does not have to
re-derive any of it. Fixed entries gain a `**Status**` line under their
classification line and are kept rather than deleted, since the repro doubles as
the regression check. None is fixed yet.

What the hunt says about `867e2b4` itself, in one paragraph: the fixes hold on the
ground each one was written for (the nulls are collected in "Ruled out"), and the
new holes sit at the edges of the injected-state machinery that the commit built.
The injection mark counts only while `self.wf0` is the marked object, so the copy
`promote_to_dense` makes retires it (finding 3); the injected state carries its
own energy into CVM, ROOTN, TD, TDZ and EX while the KPM axis stays on the band
edge (finding 1); `ground_state_on_session` went in ahead of every KPM argument
check (finding 11); the band-edge pre-fill costs an upper-edge solve that on
`"python"` is five to nine full solves (findings 9 and 10); and
`disentangle_manifold`, moved onto the chain's probe, meets an ED state that has
no probe (finding 18) and a probe whose threshold is absolute (finding 12). Six of
the eighteen came in with `867e2b4`, whole or in part; twelve are older and were
reached by probing next to it, and each entry says which. The oldest date from
2020. The widest in reach is finding 7: `gs_energy_fluctuation()`, the documented
convergence diagnostic, under-reports the variance by 10 to 51 times in the
ordinary solve-then-measure workflow on every DMRG backend, and the `maxde=`
keyword that enforces it returns states over their own tolerance.

## The four lenses

| Lens | Brief |
|---|---|
| `groundstate` | Does the injected-state machinery of `867e2b4` (`mark_injected`, `pending_injection`, `gs_is_current`, `ground_state_on_session`, the send-cache key, `gs_energy_single`'s unswept take with e0 = <wf\|H\|wf>, the band-edge pre-fill, `__deepcopy__`'s reset, `gs_energy(wf0=x)`) answer the question asked in every state a chain can be in, and do the `n_gs` snapshot and `dcex`'s projection leave the chain and the session consistent? |
| `kpm` | Is the 1.5\*\|\|vi\|\|\*\|\|vj\|\| moment bound a bound on every route that now carries it, computed from the vectors the recursion uses, silent on every correct run, and did moving the `kpm_n_scale` and `kpm_energy_truncate` checks into ED's Hermitian branch and the `i=None`/`j=None` defaults leave a route answering a different question? |
| `realtime` | After ED `evolution_ABC` moved to e^{-iHt}, DMRG `evolve_and_measure` stopped conjugating and `sxt_to_skomega` moved to a conjugate after the spatial sum on the direct sum, is every real-time consumer on one convention, on every integrator and backend, on non-Hermitian and fermionic pairs and complex Hamiltonians? |
| `misc` | Do the four misc fixes hold beyond the case each was written for: `disentangle_manifold` on the chain's probe with `eigh` of the symmetrized matrix, `kpm_finite`'s key rejection, the Kondo grids sorted with `S`, and `_strip_bond_padding` keyed on the state? |

## Scope

Out of scope by construction: the eighteen findings of
`audit_2026_09_24b_hole_hunt.md` as originally stated (a new defect in how one of
them was fixed is exactly in scope), together with that record's "Ruled out",
"New leads" and left-open items; every entry of the three records before it and
their New leads (which is what refuted the one candidate above); the open
`docs/known_issue_*.md` items, in particular the `kpm_energy_truncate` window
problem; the legacy bugs `CLAUDE.md` says are deliberately reproduced
(`evoloperator`'s z^3/6 term on `H2`, the `"moise"` key, the unreachable
`"tevol_fit_td"` branch); vendored ITensor; gaps `ROADMAP.md` marks as absent; and
the two documented KPM residuals (the `sqrt(1-x^2)` width profile and the Jackson
line shape). A finding that lies outside `867e2b4` but was reached by probing next
to it is reported with its origin named, since a hole the suite passes through is
worth the record whatever commit brought it in; the at-a-glance table carries that
origin as a column. One entry, finding 14, widens a failure that an example
comment, `docs/documentation.md` and `ROADMAP.md` record as an isolated,
intermittent edge case; it is in scope because none of the three says how wide it
is, and it is marked as such.

`itensor_version="julia_live"` was out of scope except by reading: its juliacall
JIT costs over a minute per function signature, which a three-process cap cannot
absorb, and `867e2b4`'s two Julia changes (`mpsjulialive/kpm.jl`,
`mpsjulialive/timedependent.py`) are pinned by `julia_live`-parametrized tests of
the second pass. Where an entry names `julia_live`, it says "by reading".

As a baseline, the four regression files of the second pass were run on the frame
before any lens reported:

```bash
pytest tests/test_audit_2026_09_24b_groundstate.py tests/test_audit_2026_09_24b_kpm.py tests/test_audit_2026_09_24b_realtime.py tests/test_audit_2026_09_24b_misc.py -q -p no:cacheprovider -k "not julia_live"
```

which gave 144 passed, 3 deselected in 111 seconds, so every finding below is a
hole the suite passes through rather than code already known to be red. The broad
39-minute selection the second pass ran as its baseline was not repeated.

From the start the whole hunt, all agents together, was capped at three
concurrent Python processes, and every script went through a small wrapper,
`run3.sh` (in "Shared helpers"), which waits for one of three `flock` slots and
then runs

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 <script> [args]
```

from the folder holding the script, so every repro below is shown as the
`run3.sh` line that was actually typed. Agents did not build, edit the tree or
run pytest during the hunt. Where an entry dates a defect by running an older
commit, the script ran against a read-only `git archive` extraction of that
commit in the scratch folder, on the pure-Python backend or ED only, since an
archive carries no compiled extension. The scripts and their outputs are spliced
into this file from disk rather than retyped, because the scratch folders they
ran in do not survive the session. The only edits are path substitutions:
`<repo>` for the checkout, `<scratch>` for the hunt's scratch folder and
`<python>` for the interpreter's prefix (in one profiler listing). The outputs drop three kinds of line and nothing else: the
wrapper's own `[run3] slot N acquired` notice (105 lines), the
missing-extension `UserWarning` (and the source line Python prints under it) that
an archived older tree prints on import (8 warnings), and ITensor's per-step
`TDVP_GSE` progress lines, 714 of them ("maxLinkDim after global subspace
expansion", "Global subspace expansion: cputime", "In applyExp, ..."), which are
kept in the one probe that prints them on purpose, finding 14's
`04_v3_gse_verbose`. A script already spliced under an earlier finding is
pointed to rather than repeated, and a reviewer's rerun of a hunter's script
unchanged is shown as a pointer plus its own output. Reviewers returned their
reports as text, so the reviewer paragraphs below are condensed from those
reports while the scripts and outputs are theirs verbatim.

## The findings at a glance

Origin says where the defect came from: `867e2b4` for what the fix pass itself
introduced, a commit or a period for what is older and was reached by probing
next to it. The last column names the fix cluster.

| # | Finding | Severity | Origin | Numbers change on fix | Fix cluster |
|---|---|---|---|---|---|
| 1 | after `set_gs` of a non-ground state, KPM measures from the solved E0 and every other DMRG submode from the state's own energy | LOW to MEDIUM | `867e2b4` (DMRG half), older (ED half) | only after `set_gs` of a state off the ground manifold | session |
| 2 | ED `submode="ED"` reads no state, neither `set_gs` nor `wf0=` | LOW | `0fd42d4` (2020), exposed by `867e2b4` | only after `set_gs` or with `wf0=` | session |
| 3 | a pending injection is dropped by `promote_to_dense` | LOW to MEDIUM | pre-`867e2b4` | only for a setter followed directly by a promotion | session |
| 4 | after `promote_to_dense` the next read re-sweeps any carried state that is not the sector ground state | MEDIUM | `867e2b4` (`"python"`), the 2026-09 finding 10 fix (v3) | only when the carried state is not the sector ground state | session |
| 5 | after `set_hamiltonian(restart=False)`, every read but the public correlator answers for the old Hamiltonian | LOW to MEDIUM | `7c0a71b`, `1b87543` | yes, after `restart=False` | session |
| 6 | `gs_energy(maxde=...)` returns the unrefined energy, and the next correlator discards the refinement | MEDIUM | `4731b5a` (2020) | yes, for every `maxde=` caller | fluctuation |
| 7 | `gs_energy_fluctuation()` measures the truncation of H\|psi> at `maxm`, 10 to 51 times low in the ordinary workflow | MEDIUM | `4731b5a` (2020), carried by `a8233ab`, `8fc78c6`, `68c96eb` | yes, every fluctuation below full bond dimension | fluctuation |
| 8 | `submode="SECTOR"` ignores `set_gs` | LOW to MEDIUM | `ed00449` | no (a wrong spectrum becomes a raise), or yes under a state check | session |
| 9 | the injected-state band-edge pre-fill costs an upper-edge solve on the first read | LOW (performance) | `867e2b4` | no | session (clean fix needs a rebuild) |
| 10 | `"python"`'s upper band edge costs 5 to 9 full solves on a Heisenberg chain | LOW to MEDIUM (performance) | `68c96eb`, `695c452`, `ea34f1a` | no (Emax moves at 1e-10) | pyitensor |
| 11 | `ground_state_on_session` runs ahead of every KPM argument check, and ahead of SECTOR | LOW | `867e2b4` | no | session |
| 12 | the DMRG Hermiticity probe's threshold is absolute, so a weak anti-Hermitian part is called Hermitian | MEDIUM | `40e526e` (2022); the `disentangle_manifold` half `867e2b4` | yes, on weakly non-Hermitian operators | hermiticity |
| 13 | `"python"`'s MPO builder turns operators with all coefficients below about 2e-7 into different operators | LOW to MEDIUM | `e448699` | yes, `"python"` small-coefficient operators | pyitensor |
| 14 | v3 `TDVP_GSE` loses its expansion at the left edge for a start pinned at site 0 | MEDIUM | `19bbee0` | yes, v3 `TDVP_GSE` from site-0 ladder or projector starts | realtime (needs a rebuild) |
| 15 | DMRG `evolve_and_measure`/`evolution_ABA` swallow every unknown keyword | LOW | `f54b3c2` (2020) | no (validation) | realtime |
| 16 | ED `evolve_and_measure`/`evolution_ABA` raise on `h=`, and default to a different `nt` | LOW | `f54b3c2`, `6e16aab` (2020) | no (a raise becomes a result) | realtime |
| 17 | the ED propagator is RK45 at default tolerances, its error set by the absolute energy | LOW to MEDIUM | `f54b3c2` (2020) | yes, ED real-time at 1e-8 to 1e-3 | realtime |
| 18 | `disentangle_manifold` raises on every ED manifold | LOW | `867e2b4` | no (a raise becomes a result) | hermiticity |

The clusters and their files: `session` (1 to 5, 8, 9, 11: `groundstate.py`'s
`set_gs` ED branch, `_take_injected_state` and `solver_key` docstring,
`manybodychain.py`'s `promote_to_dense` and `set_hamiltonian`, `dynamics.py`,
`kpmdmrg.py`, `edtk/dynamics.py`, `sectordc.py`, `dcex.py`'s docstring,
`spinchain.py`'s `n_gs` docstring and snapshot); `fluctuation` (6 and 7:
`groundstate.py`'s `maxde` block, `manybodychain.gs_energy_fluctuation`,
`vev.py`, `examples/groundstate/GS_enforce_maximum_fluctuation`, the user guide's
fluctuation paragraph); `hermiticity` (12 and 18: `mpsalgebra.is_hermitian`,
`algebra/algebra.py`'s two ED checks, `mpsalgebratk/disentangle.py`,
`edtk/edchain.py`); `pyitensor` (10 and 13: `pyitensor/chain.py::_maximum_energy`,
`pyitensor/dmrg.py`, `pyitensor/mpobuilder.py`); `realtime` (14 to 17:
`mpscpp3/chain_session.h::global_subspace_expand`, `timedependent.py`'s two
dispatchers, `mpsjulialive/timedependent.py`, `edtk/tdtk.py`). They are not fully
file-disjoint: `session` and `fluctuation` both edit `groundstate.py`, in
different functions, so those two go one after the other. Four pairs must land
together, because each pair shares one fix or one convention: 1 with 2 (the ED
frequency origin after `set_gs` must be one choice), 3 with 4 (one re-mark in
`promote_to_dense`), 5 with 11 (where the ground-state gate sits), and 6 with 7
(the `maxde` block). Finding 14 needs a rebuild of the v3 extension, and the clean
fix of finding 9 a rebuild of both.

## Findings

### 1. After `set_gs` of a state whose energy differs from the solved ground energy, `submode="KPM"` on every DMRG backend puts every line at E_n - E_0 of the solved ground state while CVM, CVM_explicit, ROOTN, TD, TDZ and EX put them at E_n - E_x, so the default submode's spectrum is rigidly shifted by E_x - E_0 (0.300 on a line at 1.000 on a 3-site chain, with the elastic line at +0.3), and `mode="ED"` uses three origins at once

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `groundstate`

**Where**: `src/dmrgpy/kpmdmrg.py:204` (the reconstructed axis,
`xs + (emin+emax)/2 - emin`, with `emin` the session's lower band edge, i.e. the
solved E_0); `pyitensor/chain.py`'s `_scaled_hamiltonian` and
`_scaled_hamiltonian_gs_anchored` and the same band edge in both
`chain_session.h`'s `kpm_dynamical_correlator`; `groundstate.py:200-202`
(`_take_injected_state` sets `e0 = <x|H|x>`, which CVM, ROOTN, TD, TDZ and EX
read); `groundstate.py:634-637` (the ED branch of `set_gs` leaves `ED_obj.e0` at
the solved energy, read by `edtk/dynamics.py:108-120`); `spinchain.py:250-254`,
the `n_gs` paragraph of `get_kondo_spectrum`, the only place the split is
written down.

`867e2b4` made a state set by hand reach the session unswept, with `e0 =
<x|H|x>`, and `set_gs`'s docstring now promises "its own energy <wf|H|wf>" to
every consumer. The submodes that take their frequency origin from the chain's
`e0` (CVM, CVM_explicit, ROOTN, TD, TDZ, EX) therefore measure x from E_x. KPM
does not read `e0`: it rescales H on the session's band edges and reconstructs
the axis relative to `emin`, which is the solved ground energy whatever state is
being measured. So one chain state gets two origins depending on the submode,
and on the default one the elastic |1> -> |1> line sits at +0.3 instead of 0,
while `gs_energy()` reports -0.85 and KPM measures from -1.15. `mode="ED"` has a
third arrangement: KPM, CVM, INV, ROOTN and TD measure the set state from the
solved E_0 (consistently with `gs_energy(mode="ED")`, which also reports -1.15),
EX measures it from its own energy, and `submode="ED"` does not read the state at
all (finding 2). It survived because every regression of `867e2b4` sets a member
of a degenerate doublet or of a pair split below `delta`, where E_x = E_0 to
within the split and the origins coincide, which is also the case the `n_gs`
docstring bounds; and because on the hunter's own 3-site chain the Sz0-Sz0 pair
cannot tell "|1> from E_1" apart from "|0> from E_0" (the two exact curves have
the same lines), so a pass on that pair says nothing about which state was read.

**Expected**: one frequency origin for one state across every submode and both
modes. On the 3-site Heisenberg chain in a uniform Bz=0.3 (E = -1.15, -0.85,
-0.15, 0.05), measured from |1> with its own energy, (S+_0,S-_0) has lines at
-0.3, 0.7 and 1.2; measured from E_0 they sit at 0.0, 1.0 and 1.5. The state's own
energy is what `gs_energy()` reports on the DMRG backends, what EX uses on both
modes, and what `set_gs`'s docstring promises.

Repro, from the hunter (`02` compares every submode after `set_gs(|1>)` against
the exact Lehmann sums from both origins; `07` runs the same sequence on a
read-only archive of the parent `e5049b0`, where every DMRG submode measured the
session's own solved state, the second pass's finding 11):

```bash
cd <scratch>/groundstate && <scratch>/run3.sh 02_setgs_excited_origin.py 2>&1 | tee 02_setgs_excited_origin.out
```

`groundstate/02_setgs_excited_origin.py`:

```python
# set_gs(|1>), |1> a NON-degenerate excited eigenstate: where does each
# submode put the lines?  Exact: Lehmann sum from |1>, with D_n = E_n - E_1
# (the state's own energy) or D_n = E_n - E_0 (the solved ground energy).
import numpy as np, io, contextlib
from dmrgpy import spinchain

def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h

n, B, delta = 3, 0.3, 0.05
es = np.linspace(-0.6,2.4,601)
def peak(y,lo=0.6):   # position of the strongest inelastic line (omega>lo)
    m = es>lo
    return es[m][np.argmax(np.real(y)[m])]
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()

ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
def dense(m):
    m = ed.get_ED_obj().MO2matrix(m)
    return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
E,V = np.linalg.eigh(dense(ham(ed,B)))
Sz0 = dense(ed.Sz[0])
M = np.abs(V.conj().T@Sz0@V[:,1])**2
print("E0..3 exact:",np.round(E[:4],6))
exact = {}
for ref,lab in ((E[1],"own"),(E[0],"solved")):
    D = E-ref
    exact[lab] = sum(M[k]*delta/np.pi/((es-D[k])**2+delta**2) for k in range(len(E)))
    print("exact, D_n=E_n-%s: strongest inelastic line at %.3f, peak %.4f"
          %("E_1" if lab=="own" else "E_0",peak(exact[lab]),np.max(exact[lab][es>0.6])))

ed.get_gs(mode="ED")  # sets the ED e0 (gs_energy alone does not)
eE,wE = ed.get_excited_states(n=2,mode="ED")
ed.set_gs(wE[1])
for sub in ("KPM","CVM","TD"):
    x,y = quiet(lambda: ed.get_dynamical_correlator(mode="ED",submode=sub,
                name=(ed.Sz[0],ed.Sz[0]),delta=delta,es=es))
    print("mode=ED     submode=%-5s line at %.3f"%(sub,peak(y)))

for v in ("python",3):
    np.random.seed(2)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,B))
    e0 = sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=2))
    x1 = ww[1]
    print("v=%s solved E0=%.6f, E(|1>)=%.6f"%(v,e0,ee[1]))
    for sub in ("KPM","CVM","ROOTN","TD","EX"):
        sc.set_gs(x1)
        kw = dict(delta=delta,es=es)
        if sub=="EX": kw["nex"]=8
        x,y = quiet(lambda: sc.get_dynamical_correlator(submode=sub,
                    name=(sc.Sz[0],sc.Sz[0]),**kw))
        dv = {k:np.max(np.abs(np.real(y)-exact[k])) for k in exact}
        print("v=%-6s submode=%-5s line at %.3f  max|C-exact_own|=%.3f max|C-exact_solved|=%.3f  gs_energy()=%.6f"
              %(v,sub,peak(y),dv["own"],dv["solved"],sc.gs_energy()))
```

`groundstate/02_setgs_excited_origin.out`:


```
E0..3 exact: [-1.15 -0.85 -0.15  0.05]
exact, D_n=E_n-E_1: strongest inelastic line at 1.000, peak 0.5358
exact, D_n=E_n-E_0: strongest inelastic line at 1.300, peak 0.5358
mode=ED     submode=KPM   line at 1.300
mode=ED     submode=CVM   line at 1.300
mode=ED     submode=TD    line at 1.300
v=python solved E0=-1.150000, E(|1>)=-0.850000
v=python submode=KPM   line at 1.300  max|C-exact_own|=1.224 max|C-exact_solved|=0.538  gs_energy()=-0.850000
v=python submode=CVM   line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=python submode=ROOTN line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=python submode=TD    line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=python submode=EX    line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=3 solved E0=-1.150000, E(|1>)=-0.850000
v=3      submode=KPM   line at 1.300  max|C-exact_own|=1.224 max|C-exact_solved|=0.538  gs_energy()=-0.850000
v=3      submode=CVM   line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=3      submode=ROOTN line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=3      submode=TD    line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=3      submode=EX    line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
```

```bash
cd <scratch>/groundstate && <scratch>/run3.sh 07_parent_compare.py 2>&1 | tee 07_parent_compare.out
```

`groundstate/07_parent_compare.py`:

```python
# The same two scenarios as 02 and 03, on the PARENT of 867e2b4 (e5049b0,
# extracted read-only with git archive into ./parent_src), "python" backend
# only (the archive has no compiled extensions).  What did they return
# before the injected-state machinery?
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),"parent_src","src"))
import numpy as np, io, contextlib
import dmrgpy; print(dmrgpy.__file__)
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def heis(sc,B=0.0):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
# scenario 03: set_gs(x) in a sector, promote_to_dense(), no read between
np.random.seed(3)
sc = spinchain.Spin_Chain(["S=1/2"]*4,itensor_version="python")
sc.maxm=20; sc.nsweeps=12
sc.set_hamiltonian(heis(sc)); sc.set_conserved_sector(Sz=0)
e0 = sc.gs_energy()
ee,ww = quiet(lambda: sc.get_excited_states(n=2)); x = ww[1]
sc.set_gs(x); sc.promote_to_dense()
b = sc.gs_energy(); wf = sc.get_gs(); xd = sc.promote_mps(x)
ov = abs(sc.overlap(xd,wf))**2/(abs(sc.overlap(xd,xd))*abs(sc.overlap(wf,wf)))
print("parent, scenario 03: gs_energy()=%.6f |<x|wf0>|^2=%.4f <Sz0Sz1>=%.6f"%(b,ov,sc.vev(sc.Sz[0]*sc.Sz[1]).real))
# scenario 02: set_gs(|1>), non-degenerate excited state, line positions
n,B,delta = 3,0.3,0.05
es = np.linspace(-0.6,2.4,601)
def peak(y,lo=0.6):
    m = es>lo; return es[m][np.argmax(np.real(y)[m])]
np.random.seed(2)
sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
sc.maxm=20; sc.nsweeps=12; sc.set_hamiltonian(heis(sc,B)); sc.gs_energy()
ee,ww = quiet(lambda: sc.get_excited_states(n=2)); x1 = ww[1]
for sub in ("KPM","CVM","TD"):
    sc.set_gs(x1)
    x,y = quiet(lambda: sc.get_dynamical_correlator(submode=sub,name=(sc.Sz[0],sc.Sz[0]),delta=delta,es=es))
    print("parent, scenario 02: submode=%-4s line at %.3f  gs_energy()=%.6f"%(sub,peak(y),sc.gs_energy()))
```

`groundstate/07_parent_compare.out`:


```
<scratch>/groundstate/parent_src/src/dmrgpy/__init__.py
parent, scenario 03: gs_energy()=-1.616025 |<x|wf0>|^2=0.0000 <Sz0Sz1>=-0.227671
parent, scenario 02: submode=KPM  line at 1.000  gs_energy()=-1.150000
parent, scenario 02: submode=CVM  line at 1.000  gs_energy()=-1.150000
parent, scenario 02: submode=TD   line at 1.000  gs_energy()=-1.150000
```


**Reviewer (CONFIRMED, NARROWED)**: the hunter's repro reproduces line for line.
Because Sz0-Sz0 cannot separate the two origins on this chain (the exact curves
coincide, and the DMRG CVM rows are 0.000 from both), the reviewer added
(S+_0,S-_0), whose exact lines are [-0.3, 0.7, 1.2] from |1>'s own energy, [0.0,
1.0, 1.5] from E_0 and [1.2] for the ground state. On `"python"`, v3 and v2 (v2
not in the hunter's run), CVM, CVM_explicit, ROOTN, TD and EX give [-0.3, 0.7,
1.2] with |C - own| = 0.000 and |C - gs| = 2.834, so they measure the injected
state from its own energy; TDZ agrees to 0.010, its method error; KPM gives [0.0,
1.0, 1.5] for S+S- and [0.3, 1.3, 1.8] for SzSz, the elastic line at +0.3 and the
emission line at 0. On ED after `set_gs(|1>)`, `gs_energy(mode="ED")` and
`ED_obj.e0` are both -1.150000; KPM, CVM, INV, ROOTN and TD give the curve from
E_0 (CVM |C - solved| = 0.000), EX the curve from E_1 (|C - own| = 0.000) and
`submode="ED"` the ground state's (lines [1.2], |C - gs| = 0.000). Nothing
documents the split for a non-degenerate state: the `n_gs` paragraph bounds it by
the degeneracy split and `_take_injected_state`'s docstring speaks of "a member of
a degenerate ground manifold". Dating: `07` shows the split came with `867e2b4`,
since the parent measured the solved state in every submode. Struck, each
explicitly:

- "`mode="ED"` puts them at E_n - E_0 in every submode": on ED, EX measures |1>
  from its own energy and `submode="ED"` ignores `set_gs`; only KPM, CVM, INV,
  ROOTN and TD are on E_n - E_0.
- An ED counterpart of the `gs_energy()` disagreement: ED reports -1.15 and
  measures from -1.15, so it is self-consistent and differs from DMRG only in
  convention. Only DMRG reports one energy and measures KPM from another.
- The hunter's Sz0-Sz0 "max|C - exact_own| = 0.000" as evidence that the
  injected state is measured; the claim stands on the S+S- measurement instead.
- Scope widened rather than struck: the DMRG split holds on v2, and for
  CVM_explicit and TDZ, not only for the four submodes and two backends named.

On severity the reviewer puts it at LOW to MEDIUM: it needs a state outside the
ground manifold, a niche use, but within that use it is silent, on the default
submode, and of unbounded size E_x - E_0; the ED half, a convention mismatch, is
LOW.

```bash
cd <scratch>/reviews/groundstate_1 && <scratch>/run3.sh 01_repro_setgs_excited_origin.py 2>&1 | tee 01_repro_setgs_excited_origin.out
```

`reviews/groundstate_1/01_repro_setgs_excited_origin.py` is `groundstate/02_setgs_excited_origin.py` (finding 1) unchanged, rerun; its output, `reviews/groundstate_1/01_repro_setgs_excited_origin.out`:


```
E0..3 exact: [-1.15 -0.85 -0.15  0.05]
exact, D_n=E_n-E_1: strongest inelastic line at 1.000, peak 0.5358
exact, D_n=E_n-E_0: strongest inelastic line at 1.300, peak 0.5358
mode=ED     submode=KPM   line at 1.300
mode=ED     submode=CVM   line at 1.300
mode=ED     submode=TD    line at 1.300
v=python solved E0=-1.150000, E(|1>)=-0.850000
v=python submode=KPM   line at 1.300  max|C-exact_own|=1.224 max|C-exact_solved|=0.538  gs_energy()=-0.850000
v=python submode=CVM   line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=python submode=ROOTN line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=python submode=TD    line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=python submode=EX    line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=3 solved E0=-1.150000, E(|1>)=-0.850000
v=3      submode=KPM   line at 1.300  max|C-exact_own|=1.224 max|C-exact_solved|=0.538  gs_energy()=-0.850000
v=3      submode=CVM   line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=3      submode=ROOTN line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=3      submode=TD    line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
v=3      submode=EX    line at 1.000  max|C-exact_own|=0.000 max|C-exact_solved|=0.689  gs_energy()=-0.850000
```

```bash
cd <scratch>/reviews/groundstate_1 && <scratch>/run3.sh 02_elastic_and_ed.py 2>&1 | tee 02_elastic_and_ed.out
```

`reviews/groundstate_1/02_elastic_and_ed.py`:

```python
# Reviewer probe for groundstate_1.
# (a) The elastic line and the de-excitation line after set_gs(|1>), |1> a
#     non-degenerate excited eigenstate: under the state's own origin the
#     elastic weight sits at w=0 and the |1>->|0> emission at w=E0-E1=-0.3;
#     under the solved-E0 origin they sit at +0.3 and 0.
# (b) Every submode, on ED, "python", v3 and v2, against three exact curves:
#     |1> from E_1 (own), |1> from E_0 (solved), and |0> from E_0 (the
#     solved ground state itself, i.e. set_gs ignored).
# (c) gs_energy() on both solvers after set_gs.
import numpy as np, io, contextlib, warnings
warnings.simplefilter("ignore")
from dmrgpy import spinchain

def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h

n, B, delta = 3, 0.3, 0.05
es = np.linspace(-0.8,2.4,641)
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()

ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
def dense(m):
    m = ed.get_ED_obj().MO2matrix(m)
    return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
E,V = np.linalg.eigh(dense(ham(ed,B)))
print("E exact:",np.round(E,6))
Szt = dense(ed.Sz[0]+ed.Sz[1]+ed.Sz[2])
print("Sz_tot of |0>,|1>: %.3f %.3f"%(np.real(V[:,0].conj()@Szt@V[:,0]),np.real(V[:,1].conj()@Szt@V[:,1])))

def exact_curves(A,B_):
    Am, Bm = dense(A), dense(B_)
    out = {}
    for lab,st,ref in (("own",V[:,1],E[1]),("solved",V[:,1],E[0]),("gs",V[:,0],E[0])):
        M = (st.conj()@Am@V)*(V.conj().T@Bm@st)
        D = E-ref
        out[lab] = sum(M[k]*delta/np.pi/((es-D[k])**2+delta**2) for k in range(len(E)))
    return out

def lines(y,thr=0.05):
    """positions of local maxima of Re y in es above thr*max"""
    y = np.real(y); m = np.max(y)
    i = [k for k in range(1,len(y)-1) if y[k]>=y[k-1] and y[k]>y[k+1] and y[k]>thr*m]
    return [round(es[k],3) for k in i]

pairs = {"SzSz":lambda c:(c.Sz[0],c.Sz[0]), "SpSm":lambda c:(c.Sx[0]+1j*c.Sy[0],c.Sx[0]-1j*c.Sy[0])}
for pl,pf in pairs.items():
    ex = exact_curves(*pf(ed))
    for k in ex: print("exact %s %-6s lines at %s"%(pl,k,lines(ex[k])))

def report(tag,sub,y,ex):
    dv = {k:np.max(np.abs(np.real(y)-np.real(ex[k]))) for k in ex}
    print("%-9s %-12s lines %-34s |C-own|=%.3f |C-solved|=%.3f |C-gs|=%.3f"
          %(tag,sub,lines(y),dv["own"],dv["solved"],dv["gs"]))

# ---- ED
ed.get_gs(mode="ED")
eE,wE = ed.get_excited_states(n=2,mode="ED")
ed.set_gs(wE[1])
print("ED: gs_energy(mode=ED) after set_gs(|1>) = %.6f ; ED_obj.e0 = %.6f"
      %(ed.gs_energy(mode="ED"),ed.get_ED_obj().e0))
for pl,pf in pairs.items():
    ex = exact_curves(*pf(ed))
    for sub in ("KPM","CVM","INV","ROOTN","TD","ED","EX"):
        kw = dict(delta=delta,es=es)
        if sub=="EX": kw["nex"]=8
        try:
            x,y = quiet(lambda: ed.get_dynamical_correlator(mode="ED",submode=sub,
                        name=pf(ed),**kw))
            report("ED "+pl,sub,y,ex)
        except Exception as err:
            print("ED %s %s raised %s: %s"%(pl,sub,type(err).__name__,str(err)[:120]))

# ---- DMRG
for v in ("python",3,2):
    np.random.seed(2)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,B))
    e0 = sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=2))
    x1 = ww[1]
    print("v=%s solved E0=%.6f, E(|1>)=%.6f"%(v,e0,ee[1]))
    for pl,pf in pairs.items():
        ex = exact_curves(*pf(ed))
        for sub in ("KPM","CVM","CVM_explicit","ROOTN","TD","TDZ","EX"):
            sc.set_gs(x1)
            kw = dict(delta=delta,es=es)
            if sub=="EX": kw["nex"]=8
            if sub in ("TD","TDZ"): kw["dt"]=0.05
            try:
                x,y = quiet(lambda: sc.get_dynamical_correlator(submode=sub,
                            name=pf(sc),**kw))
                report("v=%s %s"%(v,pl),sub,y,ex)
            except Exception as err:
                print("v=%s %s %s raised %s: %s"%(v,pl,sub,type(err).__name__,str(err)[:120]))
        print("   gs_energy() after = %.6f"%sc.gs_energy())
```

`reviews/groundstate_1/02_elastic_and_ed.out`:


```
E exact: [-1.15 -0.85 -0.15  0.05  0.15  0.35  0.65  0.95]
Sz_tot of |0>,|1>: -0.500 0.500
exact SzSz own    lines at [np.float64(0.0), np.float64(1.0), np.float64(1.5)]
exact SzSz solved lines at [np.float64(0.3), np.float64(1.3), np.float64(1.8)]
exact SzSz gs     lines at [np.float64(0.0), np.float64(1.0), np.float64(1.5)]
exact SpSm own    lines at [np.float64(-0.3), np.float64(0.7), np.float64(1.2)]
exact SpSm solved lines at [np.float64(0.0), np.float64(1.0), np.float64(1.5)]
exact SpSm gs     lines at [np.float64(1.2)]
ED: gs_energy(mode=ED) after set_gs(|1>) = -1.150000 ; ED_obj.e0 = -1.150000
ED SzSz   KPM          lines [np.float64(0.3), np.float64(1.3), np.float64(1.8)] |C-own|=1.225 |C-solved|=0.538 |C-gs|=1.225
ED SzSz   CVM          lines [np.float64(0.3), np.float64(1.3), np.float64(1.8)] |C-own|=0.689 |C-solved|=0.000 |C-gs|=0.689
ED SzSz   INV          lines [np.float64(0.3), np.float64(1.3), np.float64(1.8)] |C-own|=0.689 |C-solved|=0.000 |C-gs|=0.689
ED SzSz   ROOTN        lines [np.float64(0.3), np.float64(1.3), np.float64(1.8)] |C-own|=0.689 |C-solved|=0.000 |C-gs|=0.689
ED SzSz   TD           lines [np.float64(0.3), np.float64(1.3), np.float64(1.8)] |C-own|=0.689 |C-solved|=0.000 |C-gs|=0.689
ED SzSz   ED           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
ED SzSz   EX           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
ED SpSm   KPM          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=6.042 |C-solved|=3.295 |C-gs|=6.128
ED SpSm   CVM          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=2.755 |C-solved|=0.000 |C-gs|=2.833
ED SpSm   INV          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=2.755 |C-solved|=0.000 |C-gs|=2.833
ED SpSm   ROOTN        lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=2.755 |C-solved|=0.000 |C-gs|=2.833
ED SpSm   TD           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=2.755 |C-solved|=0.000 |C-gs|=2.833
ED SpSm   ED           lines [np.float64(1.2)]                  |C-own|=2.834 |C-solved|=2.833 |C-gs|=0.000
ED SpSm   EX           lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=python solved E0=-1.150000, E(|1>)=-0.850000
v=python SzSz KPM          lines [np.float64(0.3), np.float64(1.3), np.float64(1.8)] |C-own|=1.224 |C-solved|=0.538 |C-gs|=1.224
v=python SzSz CVM          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=python SzSz CVM_explicit lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=python SzSz ROOTN        lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=python SzSz TD           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=python SzSz TDZ          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.015 |C-solved|=0.687 |C-gs|=0.015
v=python SzSz EX           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
   gs_energy() after = -0.850000
v=python SpSm KPM          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=6.038 |C-solved|=3.291 |C-gs|=6.124
v=python SpSm CVM          lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=python SpSm CVM_explicit lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=python SpSm ROOTN        lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=python SpSm TD           lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=python SpSm TDZ          lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.010 |C-solved|=2.748 |C-gs|=2.827
v=python SpSm EX           lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
   gs_energy() after = -0.850000
v=3 solved E0=-1.150000, E(|1>)=-0.850000
v=3 SzSz  KPM          lines [np.float64(0.3), np.float64(1.3), np.float64(1.8)] |C-own|=1.224 |C-solved|=0.538 |C-gs|=1.224
v=3 SzSz  CVM          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=3 SzSz  CVM_explicit lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=3 SzSz  ROOTN        lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=3 SzSz  TD           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=3 SzSz  TDZ          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.015 |C-solved|=0.687 |C-gs|=0.015
v=3 SzSz  EX           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
   gs_energy() after = -0.850000
v=3 SpSm  KPM          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=6.038 |C-solved|=3.291 |C-gs|=6.124
v=3 SpSm  CVM          lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=3 SpSm  CVM_explicit lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=3 SpSm  ROOTN        lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=3 SpSm  TD           lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=3 SpSm  TDZ          lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.010 |C-solved|=2.748 |C-gs|=2.827
v=3 SpSm  EX           lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
   gs_energy() after = -0.850000
v=2 solved E0=-1.150000, E(|1>)=-0.850000
v=2 SzSz  KPM          lines [np.float64(0.3), np.float64(1.3), np.float64(1.8)] |C-own|=1.224 |C-solved|=0.538 |C-gs|=1.224
v=2 SzSz  CVM          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=2 SzSz  CVM_explicit lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=2 SzSz  ROOTN        lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
v=2 SzSz  TD           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.002 |C-solved|=0.689 |C-gs|=0.002
v=2 SzSz  TDZ          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.016 |C-solved|=0.687 |C-gs|=0.016
v=2 SzSz  EX           lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=0.000 |C-solved|=0.689 |C-gs|=0.000
   gs_energy() after = -0.850000
v=2 SpSm  KPM          lines [np.float64(0.0), np.float64(1.0), np.float64(1.5)] |C-own|=6.038 |C-solved|=3.291 |C-gs|=6.124
v=2 SpSm  CVM          lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=2 SpSm  CVM_explicit lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=2 SpSm  ROOTN        lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
v=2 SpSm  TD           lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.002 |C-solved|=2.755 |C-gs|=2.834
v=2 SpSm  TDZ          lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.010 |C-solved|=2.749 |C-gs|=2.828
v=2 SpSm  EX           lines [np.float64(-0.3), np.float64(0.7), np.float64(1.2)] |C-own|=0.000 |C-solved|=2.755 |C-gs|=2.834
   gs_energy() after = -0.850000
```

```bash
cd <scratch>/reviews/groundstate_1 && <scratch>/run3.sh 03_fix_safety.py 2>&1 | tee 03_fix_safety.out
```

`reviews/groundstate_1/03_fix_safety.py`:

```python
# Reviewer probe for groundstate_1: does the suggested fix hold?
# (a) DMRG half: at a unique ground state, is the band edge emin the session
#     returns equal to self.e0 (so shifting by e0 instead of emin moves no
#     default number)?  "python", v3, v2, and "python"/v3 under
#     kpm_energy_truncate=True.  After set_gs(|1>), reconstruct the same
#     moments with the axis shifted by self.e0 and read the lines.
# (b) ED half: simulate the suggested ED fix (ED_obj.e0 = <v|H|v>) by setting
#     the attribute, and see whether ED KPM then follows, given that e0 is
#     also ED KPM's window edge (m = H - e0, width = top of m).
import numpy as np, io, contextlib, warnings
warnings.simplefilter("ignore")
from dmrgpy import spinchain, kpmdmrg

def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, B, delta = 3, 0.3, 0.05
es = np.linspace(-0.8,2.4,641)
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def lines(y,thr=0.05):
    y = np.real(y); m = np.max(y)
    i = [k for k in range(1,len(y)-1) if y[k]>=y[k-1] and y[k]>y[k+1] and y[k]>thr*m]
    return [float(round(es[k],3)) for k in i]
def weight(y): return float(np.trapezoid(np.real(y),es))

ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
def dense(m):
    m = ed.get_ED_obj().MO2matrix(m)
    return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
E,V = np.linalg.eigh(dense(ham(ed,B)))
pair = lambda c:(c.Sx[0]+1j*c.Sy[0],c.Sx[0]-1j*c.Sy[0])
Am,Bm = [dense(o) for o in pair(ed)]
M = (V[:,1].conj()@Am@V)*(V.conj().T@Bm@V[:,1])
own = sum(M[k]*delta/np.pi/((es-(E[k]-E[1]))**2+delta**2) for k in range(len(E)))
print("exact own-origin S+S- lines",lines(own),"weight in es %.4f, sum M_n %.4f"
      %(weight(own),np.real(M.sum())))

# (a) DMRG
for v,trunc in (("python",False),(3,False),(2,False),("python",True),(3,True)):
    np.random.seed(2)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,B))
    sc.kpm_energy_truncate = trunc
    e0 = sc.gs_energy()
    mus,emin,emax,scale,nn,dl = quiet(lambda: kpmdmrg.dynamical_correlator_moments(
            sc,name=pair(sc),delta=delta))
    print("v=%-6s truncate=%-5s unique GS: emin-e0 = %.3e  (emin %.9f, e0 %.9f)"
          %(v,trunc,emin-e0,emin,e0))
    if trunc: continue
    ee,ww = quiet(lambda: sc.get_excited_states(n=2)); x1 = ww[1]
    sc.set_gs(x1)
    mus,emin,emax,scale,nn,dl = quiet(lambda: kpmdmrg.dynamical_correlator_moments(
            sc,name=pair(sc),delta=delta))
    print("   after set_gs(|1>): gs_energy()=%.6f emin=%.6f"%(sc.gs_energy(),emin))
    x,y = kpmdmrg.dynamical_correlator_from_moments(mus,emin,emax,scale,nn,es,delta=dl)
    shift = sc.gs_energy() - emin   # simulated fix: origin at e0, window unchanged
    x,yf = kpmdmrg.dynamical_correlator_from_moments(mus,emin,emax,scale,nn,es+shift,delta=dl)
    print("   current  KPM lines %s weight %.4f"%(lines(y),weight(y)))
    print("   e0-shift KPM lines %s weight %.4f"%(lines(yf),weight(yf)))

# (b) ED, suggested fix simulated
ed.get_gs(mode="ED")
eE,wE = ed.get_excited_states(n=2,mode="ED")
ed.set_gs(wE[1])
for lab,e0set in (("as shipped (e0 = E_0)",None),("simulated fix (e0 = E_1)",E[1])):
    if e0set is not None: ed.get_ED_obj().e0 = e0set
    for sub in ("KPM","CVM"):
        x,y = quiet(lambda: ed.get_dynamical_correlator(mode="ED",submode=sub,
                    name=pair(ed),delta=delta,es=es))
        print("ED %-25s %-4s lines %s weight %.4f max|C-own| %.3f"
              %(lab,sub,lines(y),weight(y),np.max(np.abs(np.real(y)-np.real(own)))))
```

`reviews/groundstate_1/03_fix_safety.out`:


```
exact own-origin S+S- lines [-0.3, 0.7, 1.2] weight in es 0.8088, sum M_n 0.8333
v=python truncate=False unique GS: emin-e0 = 0.000e+00  (emin -1.150000000, e0 -1.150000000)
   after set_gs(|1>): gs_energy()=-0.850000 emin=-1.150000
   current  KPM lines [0.0, 1.0, 1.5] weight 0.8333
   e0-shift KPM lines [-0.3, 0.7, 1.2] weight 0.8333
v=3      truncate=False unique GS: emin-e0 = 0.000e+00  (emin -1.150000000, e0 -1.150000000)
   after set_gs(|1>): gs_energy()=-0.850000 emin=-1.150000
   current  KPM lines [0.0, 1.0, 1.5] weight 0.8333
   e0-shift KPM lines [-0.3, 0.7, 1.2] weight 0.8333
v=2      truncate=False unique GS: emin-e0 = 0.000e+00  (emin -1.150000000, e0 -1.150000000)
   after set_gs(|1>): gs_energy()=-0.850000 emin=-1.150000
   current  KPM lines [0.0, 1.0, 1.5] weight 0.8333
   e0-shift KPM lines [-0.3, 0.7, 1.2] weight 0.8333
v=python truncate=True  unique GS: emin-e0 = 0.000e+00  (emin -1.150000000, e0 -1.150000000)
v=3      truncate=True  unique GS: emin-e0 = 0.000e+00  (emin -1.150000000, e0 -1.150000000)
ED as shipped (e0 = E_0)     KPM  lines [0.0, 1.0, 1.5] weight 0.8333 max|C-own| 6.042
ED as shipped (e0 = E_0)     CVM  lines [0.0, 1.0, 1.5] weight 0.8135 max|C-own| 2.755
ED simulated fix (e0 = E_1)  KPM  lines [-0.295, 0.7] weight 0.5763 max|C-own| 9.790
ED simulated fix (e0 = E_1)  CVM  lines [-0.3, 0.7, 1.2] weight 0.8088 max|C-own| 0.000
```


**Suggested fix**: on the DMRG side, shift the reconstructed KPM axis by the
state's energy rather than by `emin`, keeping the rescaling window on the band
edges. `dynamical_correlator_from_moments` has no `self`, so the origin has to be
passed in, and it has two callers, `kpmdmrg.py:89` and `mpsjulialive/
dynamics.py:190`. The reviewer measured that this moves no default number: at a
unique ground state `emin - e0` is exactly 0.000e+00 on `"python"`, v3 and v2,
and also under `kpm_energy_truncate=True` on `"python"` and v3; after
`set_gs(|1>)` the shifted axis gives the exact own-origin lines with the weight
kept (0.8333, the exact sum of M_n). On the ED side the hunter's suggestion,
setting `ED_obj.e0 = <v|H|v>`, is wrong: in `edtk/dynamics.py::
dynamical_correlator_kpm` `e0` is both the frequency origin and the bottom of the
rescaling window, and the key of `_kpm_emax_cache`, so moving it pushes the E_0
pole to x = -0.952, just outside the evaluation domain, and ED KPM loses a third
of its weight (0.5763 of 0.8333, maximum deviation 9.79, measured). The two
numbers have to be split: the window stays on the lowest eigenvalue of H, and the
state's energy goes in separately, as the origin for CVM, INV, ROOTN and TD and as
an `es` shift for KPM. Whether `gs_energy(mode="ED")` should then report <v|H|v>
after `set_gs`, as DMRG does, is a choice to make with it. The origin chosen here
is the one finding 2 has to use too. The `n_gs` docstring's KPM caveat goes away
once both halves land. NUMBERS CHANGE only after `set_gs`, `set_initial_wf` or
`set_initial_wf_guess` of a state whose energy differs from the solved ground
energy.

### 2. `mode="ED"` `submode="ED"` reads no state at all, neither the one `set_gs` set nor an explicit `wf0=`, so after `set_gs` of a non-degenerate excited eigenstate it returns the ground state's spectrum, 2.834 off the set state's exact density on a 2.835 peak, while CVM, INV, ROOTN, TD and KPM on the same ED chain read the set state and EX reads it from its own energy

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `groundstate` (turned up by the reviewer of finding 1)

**Where**: `src/dmrgpy/edtk/dynamics.py:109-111` (the `submode=="ED"` branch
passes `emu, vs` to `dynamical_correlator_ED`, whose signature
`(h, a0, b0, delta, emu, vs, dex, es)` has no state argument; the local `wf0`
computed two lines above is never forwarded).

`submode="ED"` is the exact full-spectrum Lehmann route, and its docstring and
the user guide's "dex cutoff" section define it as the equal-weight average over
the eigenstates below `dex`, measured from E_0. It builds that manifold from the
eigenvectors alone, so it ignores whatever state the chain holds and whatever
`wf0=` the caller passes. For a degenerate ground manifold this is the "split
between the two ED conventions" that the 2026-09-24 record's finding 8 left open,
and there an average is a defensible convention. For a non-degenerate state it is
not a convention at all: the state is simply not read. It survived because no
shipped caller sets a state and then calls `submode="ED"` (`atomtk/iets.py`,
`manybodychain.py:850`, the Kondo second-order checks and the METTS reference all
call it on a freshly solved chain), and because
`test_every_submode_measures_the_member_set_gs_set` covers only the DMRG
backends.

**Expected**: the density of the state the chain holds, or of the explicit
`wf0=`, as every other ED submode and every DMRG backend now return; or a
refusal.

Repro, from the reviewer of finding 1, who found it (`02_elastic_and_ed` compares
every ED submode after `set_gs(|1>)` against three exact Lehmann curves):

`reviews/groundstate_1/02_elastic_and_ed.py` and its output are spliced under finding 1 above.


**Reviewer (CONFIRMED, NARROWED)**: reproduced with an ED-only script at HEAD. On
the 3-site chain at Bz=0.3 after `set_gs(wE[1])`, `submode="ED"` is 0.000 from
the |0> curve and 2.834 from |1>'s own, on a 2.835 peak; EX is 0.000 from |1>'s
own and CVM 0.000 from |1> measured from E_0. The size is an exact Lehmann sum,
so nothing needs converging. The same mechanism is behind the recorded
degenerate case: at B=0, after `set_gs` of the pure Sz=+1/2 member of the ground
doublet, `submode="ED"` is 0.000 from the manifold average and 1.417 from the
member, whose peak is 2.835, while EX and CVM give the member to 0.000. On the
parent `e5049b0` (a read-only archive, ED only) `submode="ED"` gives identical
numbers, and there EX agreed with it bit for bit: `867e2b4` moved EX onto the
chain's state, as the second pass's finding 15 intended, and so exposed the older
gap without introducing a wrong behaviour; the branch has never read a state
(`0fd42d4`, 2020-03-21). With no `set_gs` at all and an explicit `wf0=wE[1]`
(|<wE[1]|V1>|^2 = 1.000000), `submode="ED"` still returns the |0> curve while
CVM, INV and ROOTN take the explicit state. This makes the second pass's
finding 11 Status sentence, "`mode="ED"` honours the same `set_gs` to 3.3e-6",
false for this submode. Struck, each explicitly:

- The parenthetical on which origin the other ED submodes use after
  `set_gs(excited)`: that is finding 1, not a property of this submode.
- KPM's large deviations in the reviewer-of-1's table as evidence: at
  `delta=0.05` they are dominated by the documented Jackson line shape and width
  profile.
- "The split against EX dates at least from `867e2b4`" as a regression: the
  submode is unchanged across `867e2b4`; the gap dates from `0fd42d4`.
- "Every other ED submode measures |1>" is kept, but for a degenerate member it
  is the split already recorded; what is new is the non-degenerate case and the
  ignored explicit `wf0=`.

```bash
cd <scratch>/reviews/groundstate_1_new1 && <scratch>/run3.sh 01_ed_submode_state.py 2>&1 | tee 01_ed_submode_state.out
```

`reviews/groundstate_1_new1/01_ed_submode_state.py`:

```python
# Reviewer probe for groundstate_1_new1, ED only, HEAD (867e2b4).
# (a) the hunter's case: set_gs(|1>), |1> a non-degenerate excited eigenstate,
#     on a 3-site Heisenberg chain with Bz=0.3: submode ED against EX and CVM.
# (b) the recorded case: B=0, degenerate ground doublet, set_gs(one member):
#     submode ED against CVM/EX and against the member's density and the
#     manifold average -- is this the same mechanism as the "split between
#     the two ED conventions" recorded in the 2026-09-24 record's finding 8?
# (c) does submode ED read ANY state it is given: it has no wf0 route at all.
import numpy as np, io, contextlib, warnings, inspect
warnings.simplefilter("ignore")
import dmrgpy; print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain
from dmrgpy.edtk import dynamics as edd

def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, delta = 3, 0.05
es = np.linspace(-0.8,2.4,641)
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()

def setup(B):
    ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
    def dense(m):
        m = ed.get_ED_obj().MO2matrix(m)
        return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
    E,V = np.linalg.eigh(dense(ham(ed,B)))
    return ed,dense,E,V

def curve(dense,E,st,ref,A,B_):
    Am,Bm = dense(A),dense(B_)
    M = (st.conj()@Am@V_[0])*(V_[0].conj().T@Bm@st)
    D = E-ref
    return sum(M[k]*delta/np.pi/((es-D[k])**2+delta**2) for k in range(len(E)))

def run(ed,sub,pair):
    kw = dict(delta=delta,es=es)
    if sub=="EX": kw["nex"]=8
    x,y = quiet(lambda: ed.get_dynamical_correlator(mode="ED",submode=sub,name=pair,**kw))
    return np.array(y)

# (a) hunter's case
ed,dense,E,V = setup(0.3); V_=[V]
print("(a) Bz=0.3 E:",np.round(E,4))
ed.get_gs(mode="ED")
eE,wE = ed.get_excited_states(n=2,mode="ED")
ed.set_gs(wE[1])
pair = (ed.Sx[0]+1j*ed.Sy[0], ed.Sx[0]-1j*ed.Sy[0])
ref = {"|1> own E1": curve(dense,E,V[:,1],E[1],*pair),
       "|1> from E0": curve(dense,E,V[:,1],E[0],*pair),
       "|0> from E0": curve(dense,E,V[:,0],E[0],*pair)}
print("   peak of |1> own curve %.3f, of |0> curve %.3f"%(np.max(np.real(ref["|1> own E1"])),np.max(np.real(ref["|0> from E0"]))))
for sub in ("ED","EX","CVM","KPM"):
    y = run(ed,sub,pair)
    print("   SpSm %-4s "%sub+"  ".join("max|y-%s|=%.3f"%(k,np.max(np.abs(y-r))) for k,r in ref.items()))

# (b) degenerate doublet, B=0
ed,dense,E,V = setup(0.0); V_=[V]
print("(b) B=0 E:",np.round(E,4))
ed.get_gs(mode="ED")
Szt = dense(ed.Sz[0]+ed.Sz[1]+ed.Sz[2])
# pure Sz=+1/2 member of the ground doublet
P = V[:,:2]; w,U = np.linalg.eigh(P.conj().T@Szt@P); up = P@U[:,1]
from dmrgpy.edtk.edchain import State
upst = State(up,ed.get_ED_obj())
ed.set_gs(upst)
pair = (ed.Sx[0]+1j*ed.Sy[0], ed.Sx[0]-1j*ed.Sy[0])
rm = curve(dense,E,up,E[0],*pair)
ra = 0.5*(curve(dense,E,V[:,0],E[0],*pair)+curve(dense,E,V[:,1],E[0],*pair))
print("   peak member %.3f, peak manifold average %.3f, max|member-average| %.3f"
      %(np.max(np.real(rm)),np.max(np.real(ra)),np.max(np.abs(rm-ra))))
for sub in ("ED","EX","CVM","KPM"):
    y = run(ed,sub,pair)
    print("   SpSm %-4s max|y-member|=%.3f max|y-manifold avg|=%.3f"%(sub,np.max(np.abs(y-rm)),np.max(np.abs(y-ra))))

# (c) the ED submode's own signature: which arguments can carry a state
print("(c) dynamical_correlator_ED signature:",inspect.signature(edd.dynamical_correlator_ED))
print("    dynamical_correlator_finite_T signature:",inspect.signature(edd.dynamical_correlator_finite_T))
```

`reviews/groundstate_1_new1/01_ed_submode_state.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
(a) Bz=0.3 E: [-1.15 -0.85 -0.15  0.05  0.15  0.35  0.65  0.95]
   peak of |1> own curve 2.835, of |0> curve 1.061
   SpSm ED   max|y-|1> own E1|=2.834  max|y-|1> from E0|=2.833  max|y-|0> from E0|=0.000
   SpSm EX   max|y-|1> own E1|=0.000  max|y-|1> from E0|=2.755  max|y-|0> from E0|=2.834
   SpSm CVM  max|y-|1> own E1|=2.755  max|y-|1> from E0|=0.000  max|y-|0> from E0|=2.833
   SpSm KPM  max|y-|1> own E1|=6.042  max|y-|1> from E0|=3.295  max|y-|0> from E0|=6.128
(b) B=0 E: [-1.  -1.  -0.  -0.   0.5  0.5  0.5  0.5]
   peak member 2.835, peak manifold average 1.418, max|member-average| 1.417
   SpSm ED   max|y-member|=1.417 max|y-manifold avg|=0.000
   SpSm EX   max|y-member|=0.000 max|y-manifold avg|=1.417
   SpSm CVM  max|y-member|=0.000 max|y-manifold avg|=1.417
   SpSm KPM  max|y-member|=3.327 max|y-manifold avg|=4.744
(c) dynamical_correlator_ED signature: (h, a0, b0, delta=0.02, emu=None, vs=None, dex=1e-05, es=array([-1.00000000e+00, -9.81636060e-01, -9.63272120e-01, -9.44908180e-01,
       -9.26544240e-01, -9.08180301e-01, -8.89816361e-01, -8.71452421e-01,
       -8.53088481e-01, -8.34724541e-01, -8.16360601e-01, -7.97996661e-01,
       -7.79632721e-01, -7.61268781e-01, -7.42904841e-01, -7.24540902e-01,
       -7.06176962e-01, -6.87813022e-01, -6.69449082e-01, -6.51085142e-01,
       -6.32721202e-01, -6.14357262e-01, -5.95993322e-01, -5.77629382e-01,
       -5.59265442e-01, -5.40901503e-01, -5.22537563e-01, -5.04173623e-01,
       -4.85809683e-01, -4.67445743e-01, -4.49081803e-01, -4.30717863e-01,
       -4.12353923e-01, -3.93989983e-01, -3.75626043e-01, -3.57262104e-01,
       -3.38898164e-01, -3.20534224e-01, -3.02170284e-01, -2.83806344e-01,
       -2.65442404e-01, -2.47078464e-01, -2.28714524e-01, -2.10350584e-01,
       -1.91986644e-01, -1.73622705e-01, -1.55258765e-01, -1.36894825e-01,
       -1.18530885e-01, -1.00166945e-01, -8.18030050e-02, -6.34390651e-02,
       -4.50751252e-02, -2.67111853e-02, -8.34724541e-03,  1.00166945e-02,
        2.83806344e-02,  4.67445743e-02,  6.51085142e-02,  8.34724541e-02,
        1.01836394e-01,  1.20200334e-01,  1.38564274e-01,  1.56928214e-01,
        1.75292154e-01,  1.93656093e-01,  2.12020033e-01,  2.30383973e-01,
        2.48747913e-01,  2.67111853e-01,  2.85475793e-01,  3.03839733e-01,
        3.22203673e-01,  3.40567613e-01,  3.58931553e-01,  3.77295492e-01,
        3.95659432e-01,  4.14023372e-01,  4.32387312e-01,  4.50751252e-01,
        4.69115192e-01,  4.87479132e-01,  5.05843072e-01,  5.24207012e-01,
        5.42570952e-01,  5.60934891e-01,  5.79298831e-01,  5.97662771e-01,
        6.16026711e-01,  6.34390651e-01,  6.52754591e-01,  6.71118531e-01,
        6.89482471e-01,  7.07846411e-01,  7.26210351e-01,  7.44574290e-01,
        7.62938230e-01,  7.81302170e-01,  7.99666110e-01,  8.18030050e-01,
        8.36393990e-01,  8.54757930e-01,  8.73121870e-01,  8.91485810e-01,
        9.09849750e-01,  9.28213689e-01,  9.46577629e-01,  9.64941569e-01,
        9.83305509e-01,  1.00166945e+00,  1.02003339e+00,  1.03839733e+00,
        1.05676127e+00,  1.07512521e+00,  1.09348915e+00,  1.11185309e+00,
        1.13021703e+00,  1.14858097e+00,  1.16694491e+00,  1.18530885e+00,
        1.20367279e+00,  1.22203673e+00,  1.24040067e+00,  1.25876461e+00,
        1.27712855e+00,  1.29549249e+00,  1.31385643e+00,  1.33222037e+00,
        1.35058431e+00,  1.36894825e+00,  1.38731219e+00,  1.40567613e+00,
        1.42404007e+00,  1.44240401e+00,  1.46076795e+00,  1.47913189e+00,
        1.49749583e+00,  1.51585977e+00,  1.53422371e+00,  1.55258765e+00,
        1.57095159e+00,  1.58931553e+00,  1.60767947e+00,  1.62604341e+00,
        1.64440735e+00,  1.66277129e+00,  1.68113523e+00,  1.69949917e+00,
        1.71786311e+00,  1.73622705e+00,  1.75459098e+00,  1.77295492e+00,
        1.79131886e+00,  1.80968280e+00,  1.82804674e+00,  1.84641068e+00,
        1.86477462e+00,  1.88313856e+00,  1.90150250e+00,  1.91986644e+00,
        1.93823038e+00,  1.95659432e+00,  1.97495826e+00,  1.99332220e+00,
        2.01168614e+00,  2.03005008e+00,  2.04841402e+00,  2.06677796e+00,
        2.08514190e+00,  2.10350584e+00,  2.12186978e+00,  2.14023372e+00,
        2.15859766e+00,  2.17696160e+00,  2.19532554e+00,  2.21368948e+00,
        2.23205342e+00,  2.25041736e+00,  2.26878130e+00,  2.28714524e+00,
        2.30550918e+00,  2.32387312e+00,  2.34223706e+00,  2.36060100e+00,
        2.37896494e+00,  2.39732888e+00,  2.41569282e+00,  2.43405676e+00,
        2.45242070e+00,  2.47078464e+00,  2.48914858e+00,  2.50751252e+00,
        2.52587646e+00,  2.54424040e+00,  2.56260434e+00,  2.58096828e+00,
        2.59933222e+00,  2.61769616e+00,  2.63606010e+00,  2.65442404e+00,
        2.67278798e+00,  2.69115192e+00,  2.70951586e+00,  2.72787980e+00,
        2.74624374e+00,  2.76460768e+00,  2.78297162e+00,  2.80133556e+00,
        2.81969950e+00,  2.83806344e+00,  2.85642738e+00,  2.87479132e+00,
        2.89315526e+00,  2.91151920e+00,  2.92988314e+00,  2.94824708e+00,
        2.96661102e+00,  2.98497496e+00,  3.00333890e+00,  3.02170284e+00,
        3.04006678e+00,  3.05843072e+00,  3.07679466e+00,  3.09515860e+00,
        3.11352254e+00,  3.13188648e+00,  3.15025042e+00,  3.16861436e+00,
        3.18697830e+00,  3.20534224e+00,  3.22370618e+00,  3.24207012e+00,
        3.26043406e+00,  3.27879800e+00,  3.29716194e+00,  3.31552588e+00,
        3.33388982e+00,  3.35225376e+00,  3.37061770e+00,  3.38898164e+00,
        3.40734558e+00,  3.42570952e+00,  3.44407346e+00,  3.46243740e+00,
        3.48080134e+00,  3.49916528e+00,  3.51752922e+00,  3.53589316e+00,
        3.55425710e+00,  3.57262104e+00,  3.59098497e+00,  3.60934891e+00,
        3.62771285e+00,  3.64607679e+00,  3.66444073e+00,  3.68280467e+00,
        3.70116861e+00,  3.71953255e+00,  3.73789649e+00,  3.75626043e+00,
        3.77462437e+00,  3.79298831e+00,  3.81135225e+00,  3.82971619e+00,
        3.84808013e+00,  3.86644407e+00,  3.88480801e+00,  3.90317195e+00,
        3.92153589e+00,  3.93989983e+00,  3.95826377e+00,  3.97662771e+00,
        3.99499165e+00,  4.01335559e+00,  4.03171953e+00,  4.05008347e+00,
        4.06844741e+00,  4.08681135e+00,  4.10517529e+00,  4.12353923e+00,
        4.14190317e+00,  4.16026711e+00,  4.17863105e+00,  4.19699499e+00,
        4.21535893e+00,  4.23372287e+00,  4.25208681e+00,  4.27045075e+00,
        4.28881469e+00,  4.30717863e+00,  4.32554257e+00,  4.34390651e+00,
        4.36227045e+00,  4.38063439e+00,  4.39899833e+00,  4.41736227e+00,
        4.43572621e+00,  4.45409015e+00,  4.47245409e+00,  4.49081803e+00,
        4.50918197e+00,  4.52754591e+00,  4.54590985e+00,  4.56427379e+00,
        4.58263773e+00,  4.60100167e+00,  4.61936561e+00,  4.63772955e+00,
        4.65609349e+00,  4.67445743e+00,  4.69282137e+00,  4.71118531e+00,
        4.72954925e+00,  4.74791319e+00,  4.76627713e+00,  4.78464107e+00,
        4.80300501e+00,  4.82136895e+00,  4.83973289e+00,  4.85809683e+00,
        4.87646077e+00,  4.89482471e+00,  4.91318865e+00,  4.93155259e+00,
        4.94991653e+00,  4.96828047e+00,  4.98664441e+00,  5.00500835e+00,
        5.02337229e+00,  5.04173623e+00,  5.06010017e+00,  5.07846411e+00,
        5.09682805e+00,  5.11519199e+00,  5.13355593e+00,  5.15191987e+00,
        5.17028381e+00,  5.18864775e+00,  5.20701169e+00,  5.22537563e+00,
        5.24373957e+00,  5.26210351e+00,  5.28046745e+00,  5.29883139e+00,
        5.31719533e+00,  5.33555927e+00,  5.35392321e+00,  5.37228715e+00,
        5.39065109e+00,  5.40901503e+00,  5.42737896e+00,  5.44574290e+00,
        5.46410684e+00,  5.48247078e+00,  5.50083472e+00,  5.51919866e+00,
        5.53756260e+00,  5.55592654e+00,  5.57429048e+00,  5.59265442e+00,
        5.61101836e+00,  5.62938230e+00,  5.64774624e+00,  5.66611018e+00,
        5.68447412e+00,  5.70283806e+00,  5.72120200e+00,  5.73956594e+00,
        5.75792988e+00,  5.77629382e+00,  5.79465776e+00,  5.81302170e+00,
        5.83138564e+00,  5.84974958e+00,  5.86811352e+00,  5.88647746e+00,
        5.90484140e+00,  5.92320534e+00,  5.94156928e+00,  5.95993322e+00,
        5.97829716e+00,  5.99666110e+00,  6.01502504e+00,  6.03338898e+00,
        6.05175292e+00,  6.07011686e+00,  6.08848080e+00,  6.10684474e+00,
        6.12520868e+00,  6.14357262e+00,  6.16193656e+00,  6.18030050e+00,
        6.19866444e+00,  6.21702838e+00,  6.23539232e+00,  6.25375626e+00,
        6.27212020e+00,  6.29048414e+00,  6.30884808e+00,  6.32721202e+00,
        6.34557596e+00,  6.36393990e+00,  6.38230384e+00,  6.40066778e+00,
        6.41903172e+00,  6.43739566e+00,  6.45575960e+00,  6.47412354e+00,
        6.49248748e+00,  6.51085142e+00,  6.52921536e+00,  6.54757930e+00,
        6.56594324e+00,  6.58430718e+00,  6.60267112e+00,  6.62103506e+00,
        6.63939900e+00,  6.65776294e+00,  6.67612688e+00,  6.69449082e+00,
        6.71285476e+00,  6.73121870e+00,  6.74958264e+00,  6.76794658e+00,
        6.78631052e+00,  6.80467446e+00,  6.82303840e+00,  6.84140234e+00,
        6.85976628e+00,  6.87813022e+00,  6.89649416e+00,  6.91485810e+00,
        6.93322204e+00,  6.95158598e+00,  6.96994992e+00,  6.98831386e+00,
        7.00667780e+00,  7.02504174e+00,  7.04340568e+00,  7.06176962e+00,
        7.08013356e+00,  7.09849750e+00,  7.11686144e+00,  7.13522538e+00,
        7.15358932e+00,  7.17195326e+00,  7.19031720e+00,  7.20868114e+00,
        7.22704508e+00,  7.24540902e+00,  7.26377295e+00,  7.28213689e+00,
        7.30050083e+00,  7.31886477e+00,  7.33722871e+00,  7.35559265e+00,
        7.37395659e+00,  7.39232053e+00,  7.41068447e+00,  7.42904841e+00,
        7.44741235e+00,  7.46577629e+00,  7.48414023e+00,  7.50250417e+00,
        7.52086811e+00,  7.53923205e+00,  7.55759599e+00,  7.57595993e+00,
        7.59432387e+00,  7.61268781e+00,  7.63105175e+00,  7.64941569e+00,
        7.66777963e+00,  7.68614357e+00,  7.70450751e+00,  7.72287145e+00,
        7.74123539e+00,  7.75959933e+00,  7.77796327e+00,  7.79632721e+00,
        7.81469115e+00,  7.83305509e+00,  7.85141903e+00,  7.86978297e+00,
        7.88814691e+00,  7.90651085e+00,  7.92487479e+00,  7.94323873e+00,
        7.96160267e+00,  7.97996661e+00,  7.99833055e+00,  8.01669449e+00,
        8.03505843e+00,  8.05342237e+00,  8.07178631e+00,  8.09015025e+00,
        8.10851419e+00,  8.12687813e+00,  8.14524207e+00,  8.16360601e+00,
        8.18196995e+00,  8.20033389e+00,  8.21869783e+00,  8.23706177e+00,
        8.25542571e+00,  8.27378965e+00,  8.29215359e+00,  8.31051753e+00,
        8.32888147e+00,  8.34724541e+00,  8.36560935e+00,  8.38397329e+00,
        8.40233723e+00,  8.42070117e+00,  8.43906511e+00,  8.45742905e+00,
        8.47579299e+00,  8.49415693e+00,  8.51252087e+00,  8.53088481e+00,
        8.54924875e+00,  8.56761269e+00,  8.58597663e+00,  8.60434057e+00,
        8.62270451e+00,  8.64106845e+00,  8.65943239e+00,  8.67779633e+00,
        8.69616027e+00,  8.71452421e+00,  8.73288815e+00,  8.75125209e+00,
        8.76961603e+00,  8.78797997e+00,  8.80634391e+00,  8.82470785e+00,
        8.84307179e+00,  8.86143573e+00,  8.87979967e+00,  8.89816361e+00,
        8.91652755e+00,  8.93489149e+00,  8.95325543e+00,  8.97161937e+00,
        8.98998331e+00,  9.00834725e+00,  9.02671119e+00,  9.04507513e+00,
        9.06343907e+00,  9.08180301e+00,  9.10016694e+00,  9.11853088e+00,
        9.13689482e+00,  9.15525876e+00,  9.17362270e+00,  9.19198664e+00,
        9.21035058e+00,  9.22871452e+00,  9.24707846e+00,  9.26544240e+00,
        9.28380634e+00,  9.30217028e+00,  9.32053422e+00,  9.33889816e+00,
        9.35726210e+00,  9.37562604e+00,  9.39398998e+00,  9.41235392e+00,
        9.43071786e+00,  9.44908180e+00,  9.46744574e+00,  9.48580968e+00,
        9.50417362e+00,  9.52253756e+00,  9.54090150e+00,  9.55926544e+00,
        9.57762938e+00,  9.59599332e+00,  9.61435726e+00,  9.63272120e+00,
        9.65108514e+00,  9.66944908e+00,  9.68781302e+00,  9.70617696e+00,
        9.72454090e+00,  9.74290484e+00,  9.76126878e+00,  9.77963272e+00,
        9.79799666e+00,  9.81636060e+00,  9.83472454e+00,  9.85308848e+00,
        9.87145242e+00,  9.88981636e+00,  9.90818030e+00,  9.92654424e+00,
        9.94490818e+00,  9.96327212e+00,  9.98163606e+00,  1.00000000e+01]))
    dynamical_correlator_finite_T signature: (h, a0, b0, T, delta=0.02, emu=None, vs=None, es=array([-1.00000000e+00, -9.81636060e-01, -9.63272120e-01, -9.44908180e-01,
       -9.26544240e-01, -9.08180301e-01, -8.89816361e-01, -8.71452421e-01,
       -8.53088481e-01, -8.34724541e-01, -8.16360601e-01, -7.97996661e-01,
       -7.79632721e-01, -7.61268781e-01, -7.42904841e-01, -7.24540902e-01,
       -7.06176962e-01, -6.87813022e-01, -6.69449082e-01, -6.51085142e-01,
       -6.32721202e-01, -6.14357262e-01, -5.95993322e-01, -5.77629382e-01,
       -5.59265442e-01, -5.40901503e-01, -5.22537563e-01, -5.04173623e-01,
       -4.85809683e-01, -4.67445743e-01, -4.49081803e-01, -4.30717863e-01,
       -4.12353923e-01, -3.93989983e-01, -3.75626043e-01, -3.57262104e-01,
       -3.38898164e-01, -3.20534224e-01, -3.02170284e-01, -2.83806344e-01,
       -2.65442404e-01, -2.47078464e-01, -2.28714524e-01, -2.10350584e-01,
       -1.91986644e-01, -1.73622705e-01, -1.55258765e-01, -1.36894825e-01,
       -1.18530885e-01, -1.00166945e-01, -8.18030050e-02, -6.34390651e-02,
       -4.50751252e-02, -2.67111853e-02, -8.34724541e-03,  1.00166945e-02,
        2.83806344e-02,  4.67445743e-02,  6.51085142e-02,  8.34724541e-02,
        1.01836394e-01,  1.20200334e-01,  1.38564274e-01,  1.56928214e-01,
        1.75292154e-01,  1.93656093e-01,  2.12020033e-01,  2.30383973e-01,
        2.48747913e-01,  2.67111853e-01,  2.85475793e-01,  3.03839733e-01,
        3.22203673e-01,  3.40567613e-01,  3.58931553e-01,  3.77295492e-01,
        3.95659432e-01,  4.14023372e-01,  4.32387312e-01,  4.50751252e-01,
        4.69115192e-01,  4.87479132e-01,  5.05843072e-01,  5.24207012e-01,
        5.42570952e-01,  5.60934891e-01,  5.79298831e-01,  5.97662771e-01,
        6.16026711e-01,  6.34390651e-01,  6.52754591e-01,  6.71118531e-01,
        6.89482471e-01,  7.07846411e-01,  7.26210351e-01,  7.44574290e-01,
        7.62938230e-01,  7.81302170e-01,  7.99666110e-01,  8.18030050e-01,
        8.36393990e-01,  8.54757930e-01,  8.73121870e-01,  8.91485810e-01,
        9.09849750e-01,  9.28213689e-01,  9.46577629e-01,  9.64941569e-01,
        9.83305509e-01,  1.00166945e+00,  1.02003339e+00,  1.03839733e+00,
        1.05676127e+00,  1.07512521e+00,  1.09348915e+00,  1.11185309e+00,
        1.13021703e+00,  1.14858097e+00,  1.16694491e+00,  1.18530885e+00,
        1.20367279e+00,  1.22203673e+00,  1.24040067e+00,  1.25876461e+00,
        1.27712855e+00,  1.29549249e+00,  1.31385643e+00,  1.33222037e+00,
        1.35058431e+00,  1.36894825e+00,  1.38731219e+00,  1.40567613e+00,
        1.42404007e+00,  1.44240401e+00,  1.46076795e+00,  1.47913189e+00,
        1.49749583e+00,  1.51585977e+00,  1.53422371e+00,  1.55258765e+00,
        1.57095159e+00,  1.58931553e+00,  1.60767947e+00,  1.62604341e+00,
        1.64440735e+00,  1.66277129e+00,  1.68113523e+00,  1.69949917e+00,
        1.71786311e+00,  1.73622705e+00,  1.75459098e+00,  1.77295492e+00,
        1.79131886e+00,  1.80968280e+00,  1.82804674e+00,  1.84641068e+00,
        1.86477462e+00,  1.88313856e+00,  1.90150250e+00,  1.91986644e+00,
        1.93823038e+00,  1.95659432e+00,  1.97495826e+00,  1.99332220e+00,
        2.01168614e+00,  2.03005008e+00,  2.04841402e+00,  2.06677796e+00,
        2.08514190e+00,  2.10350584e+00,  2.12186978e+00,  2.14023372e+00,
        2.15859766e+00,  2.17696160e+00,  2.19532554e+00,  2.21368948e+00,
        2.23205342e+00,  2.25041736e+00,  2.26878130e+00,  2.28714524e+00,
        2.30550918e+00,  2.32387312e+00,  2.34223706e+00,  2.36060100e+00,
        2.37896494e+00,  2.39732888e+00,  2.41569282e+00,  2.43405676e+00,
        2.45242070e+00,  2.47078464e+00,  2.48914858e+00,  2.50751252e+00,
        2.52587646e+00,  2.54424040e+00,  2.56260434e+00,  2.58096828e+00,
        2.59933222e+00,  2.61769616e+00,  2.63606010e+00,  2.65442404e+00,
        2.67278798e+00,  2.69115192e+00,  2.70951586e+00,  2.72787980e+00,
        2.74624374e+00,  2.76460768e+00,  2.78297162e+00,  2.80133556e+00,
        2.81969950e+00,  2.83806344e+00,  2.85642738e+00,  2.87479132e+00,
        2.89315526e+00,  2.91151920e+00,  2.92988314e+00,  2.94824708e+00,
        2.96661102e+00,  2.98497496e+00,  3.00333890e+00,  3.02170284e+00,
        3.04006678e+00,  3.05843072e+00,  3.07679466e+00,  3.09515860e+00,
        3.11352254e+00,  3.13188648e+00,  3.15025042e+00,  3.16861436e+00,
        3.18697830e+00,  3.20534224e+00,  3.22370618e+00,  3.24207012e+00,
        3.26043406e+00,  3.27879800e+00,  3.29716194e+00,  3.31552588e+00,
        3.33388982e+00,  3.35225376e+00,  3.37061770e+00,  3.38898164e+00,
        3.40734558e+00,  3.42570952e+00,  3.44407346e+00,  3.46243740e+00,
        3.48080134e+00,  3.49916528e+00,  3.51752922e+00,  3.53589316e+00,
        3.55425710e+00,  3.57262104e+00,  3.59098497e+00,  3.60934891e+00,
        3.62771285e+00,  3.64607679e+00,  3.66444073e+00,  3.68280467e+00,
        3.70116861e+00,  3.71953255e+00,  3.73789649e+00,  3.75626043e+00,
        3.77462437e+00,  3.79298831e+00,  3.81135225e+00,  3.82971619e+00,
        3.84808013e+00,  3.86644407e+00,  3.88480801e+00,  3.90317195e+00,
        3.92153589e+00,  3.93989983e+00,  3.95826377e+00,  3.97662771e+00,
        3.99499165e+00,  4.01335559e+00,  4.03171953e+00,  4.05008347e+00,
        4.06844741e+00,  4.08681135e+00,  4.10517529e+00,  4.12353923e+00,
        4.14190317e+00,  4.16026711e+00,  4.17863105e+00,  4.19699499e+00,
        4.21535893e+00,  4.23372287e+00,  4.25208681e+00,  4.27045075e+00,
        4.28881469e+00,  4.30717863e+00,  4.32554257e+00,  4.34390651e+00,
        4.36227045e+00,  4.38063439e+00,  4.39899833e+00,  4.41736227e+00,
        4.43572621e+00,  4.45409015e+00,  4.47245409e+00,  4.49081803e+00,
        4.50918197e+00,  4.52754591e+00,  4.54590985e+00,  4.56427379e+00,
        4.58263773e+00,  4.60100167e+00,  4.61936561e+00,  4.63772955e+00,
        4.65609349e+00,  4.67445743e+00,  4.69282137e+00,  4.71118531e+00,
        4.72954925e+00,  4.74791319e+00,  4.76627713e+00,  4.78464107e+00,
        4.80300501e+00,  4.82136895e+00,  4.83973289e+00,  4.85809683e+00,
        4.87646077e+00,  4.89482471e+00,  4.91318865e+00,  4.93155259e+00,
        4.94991653e+00,  4.96828047e+00,  4.98664441e+00,  5.00500835e+00,
        5.02337229e+00,  5.04173623e+00,  5.06010017e+00,  5.07846411e+00,
        5.09682805e+00,  5.11519199e+00,  5.13355593e+00,  5.15191987e+00,
        5.17028381e+00,  5.18864775e+00,  5.20701169e+00,  5.22537563e+00,
        5.24373957e+00,  5.26210351e+00,  5.28046745e+00,  5.29883139e+00,
        5.31719533e+00,  5.33555927e+00,  5.35392321e+00,  5.37228715e+00,
        5.39065109e+00,  5.40901503e+00,  5.42737896e+00,  5.44574290e+00,
        5.46410684e+00,  5.48247078e+00,  5.50083472e+00,  5.51919866e+00,
        5.53756260e+00,  5.55592654e+00,  5.57429048e+00,  5.59265442e+00,
        5.61101836e+00,  5.62938230e+00,  5.64774624e+00,  5.66611018e+00,
        5.68447412e+00,  5.70283806e+00,  5.72120200e+00,  5.73956594e+00,
        5.75792988e+00,  5.77629382e+00,  5.79465776e+00,  5.81302170e+00,
        5.83138564e+00,  5.84974958e+00,  5.86811352e+00,  5.88647746e+00,
        5.90484140e+00,  5.92320534e+00,  5.94156928e+00,  5.95993322e+00,
        5.97829716e+00,  5.99666110e+00,  6.01502504e+00,  6.03338898e+00,
        6.05175292e+00,  6.07011686e+00,  6.08848080e+00,  6.10684474e+00,
        6.12520868e+00,  6.14357262e+00,  6.16193656e+00,  6.18030050e+00,
        6.19866444e+00,  6.21702838e+00,  6.23539232e+00,  6.25375626e+00,
        6.27212020e+00,  6.29048414e+00,  6.30884808e+00,  6.32721202e+00,
        6.34557596e+00,  6.36393990e+00,  6.38230384e+00,  6.40066778e+00,
        6.41903172e+00,  6.43739566e+00,  6.45575960e+00,  6.47412354e+00,
        6.49248748e+00,  6.51085142e+00,  6.52921536e+00,  6.54757930e+00,
        6.56594324e+00,  6.58430718e+00,  6.60267112e+00,  6.62103506e+00,
        6.63939900e+00,  6.65776294e+00,  6.67612688e+00,  6.69449082e+00,
        6.71285476e+00,  6.73121870e+00,  6.74958264e+00,  6.76794658e+00,
        6.78631052e+00,  6.80467446e+00,  6.82303840e+00,  6.84140234e+00,
        6.85976628e+00,  6.87813022e+00,  6.89649416e+00,  6.91485810e+00,
        6.93322204e+00,  6.95158598e+00,  6.96994992e+00,  6.98831386e+00,
        7.00667780e+00,  7.02504174e+00,  7.04340568e+00,  7.06176962e+00,
        7.08013356e+00,  7.09849750e+00,  7.11686144e+00,  7.13522538e+00,
        7.15358932e+00,  7.17195326e+00,  7.19031720e+00,  7.20868114e+00,
        7.22704508e+00,  7.24540902e+00,  7.26377295e+00,  7.28213689e+00,
        7.30050083e+00,  7.31886477e+00,  7.33722871e+00,  7.35559265e+00,
        7.37395659e+00,  7.39232053e+00,  7.41068447e+00,  7.42904841e+00,
        7.44741235e+00,  7.46577629e+00,  7.48414023e+00,  7.50250417e+00,
        7.52086811e+00,  7.53923205e+00,  7.55759599e+00,  7.57595993e+00,
        7.59432387e+00,  7.61268781e+00,  7.63105175e+00,  7.64941569e+00,
        7.66777963e+00,  7.68614357e+00,  7.70450751e+00,  7.72287145e+00,
        7.74123539e+00,  7.75959933e+00,  7.77796327e+00,  7.79632721e+00,
        7.81469115e+00,  7.83305509e+00,  7.85141903e+00,  7.86978297e+00,
        7.88814691e+00,  7.90651085e+00,  7.92487479e+00,  7.94323873e+00,
        7.96160267e+00,  7.97996661e+00,  7.99833055e+00,  8.01669449e+00,
        8.03505843e+00,  8.05342237e+00,  8.07178631e+00,  8.09015025e+00,
        8.10851419e+00,  8.12687813e+00,  8.14524207e+00,  8.16360601e+00,
        8.18196995e+00,  8.20033389e+00,  8.21869783e+00,  8.23706177e+00,
        8.25542571e+00,  8.27378965e+00,  8.29215359e+00,  8.31051753e+00,
        8.32888147e+00,  8.34724541e+00,  8.36560935e+00,  8.38397329e+00,
        8.40233723e+00,  8.42070117e+00,  8.43906511e+00,  8.45742905e+00,
        8.47579299e+00,  8.49415693e+00,  8.51252087e+00,  8.53088481e+00,
        8.54924875e+00,  8.56761269e+00,  8.58597663e+00,  8.60434057e+00,
        8.62270451e+00,  8.64106845e+00,  8.65943239e+00,  8.67779633e+00,
        8.69616027e+00,  8.71452421e+00,  8.73288815e+00,  8.75125209e+00,
        8.76961603e+00,  8.78797997e+00,  8.80634391e+00,  8.82470785e+00,
        8.84307179e+00,  8.86143573e+00,  8.87979967e+00,  8.89816361e+00,
        8.91652755e+00,  8.93489149e+00,  8.95325543e+00,  8.97161937e+00,
        8.98998331e+00,  9.00834725e+00,  9.02671119e+00,  9.04507513e+00,
        9.06343907e+00,  9.08180301e+00,  9.10016694e+00,  9.11853088e+00,
        9.13689482e+00,  9.15525876e+00,  9.17362270e+00,  9.19198664e+00,
        9.21035058e+00,  9.22871452e+00,  9.24707846e+00,  9.26544240e+00,
        9.28380634e+00,  9.30217028e+00,  9.32053422e+00,  9.33889816e+00,
        9.35726210e+00,  9.37562604e+00,  9.39398998e+00,  9.41235392e+00,
        9.43071786e+00,  9.44908180e+00,  9.46744574e+00,  9.48580968e+00,
        9.50417362e+00,  9.52253756e+00,  9.54090150e+00,  9.55926544e+00,
        9.57762938e+00,  9.59599332e+00,  9.61435726e+00,  9.63272120e+00,
        9.65108514e+00,  9.66944908e+00,  9.68781302e+00,  9.70617696e+00,
        9.72454090e+00,  9.74290484e+00,  9.76126878e+00,  9.77963272e+00,
        9.79799666e+00,  9.81636060e+00,  9.83472454e+00,  9.85308848e+00,
        9.87145242e+00,  9.88981636e+00,  9.90818030e+00,  9.92654424e+00,
        9.94490818e+00,  9.96327212e+00,  9.98163606e+00,  1.00000000e+01]))
```

```bash
cd <scratch>/reviews/groundstate_1_new1 && <scratch>/run3.sh 02_parent_ed_split.py 2>&1 | tee 02_parent_ed_split.out
```

`reviews/groundstate_1_new1/02_parent_ed_split.py`:

```python
# Reviewer probe for groundstate_1_new1: the same ED-only measurements as 01,
# on the PARENT of 867e2b4 (e5049b0, the read-only git-archive extraction the
# groundstate hunter left in ../../groundstate/parent_src). Did submode ED
# read the set state before 867e2b4, and which side of the split was EX on?
import sys, os
sys.path.insert(0, "<scratch>/groundstate/parent_src/src")
import numpy as np, io, contextlib, warnings
warnings.simplefilter("ignore")
import dmrgpy; print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain
from dmrgpy.edtk.edchain import State

def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, delta = 3, 0.05
es = np.linspace(-0.8,2.4,641)
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def setup(B):
    ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
    def dense(m):
        m = ed.get_ED_obj().MO2matrix(m)
        return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
    E,V = np.linalg.eigh(dense(ham(ed,B)))
    return ed,dense,E,V
def curve(dense,E,V,st,ref,A,B_):
    Am,Bm = dense(A),dense(B_)
    M = (st.conj()@Am@V)*(V.conj().T@Bm@st)
    return sum(M[k]*delta/np.pi/((es-(E[k]-ref))**2+delta**2) for k in range(len(E)))
def run(ed,sub,pair):
    kw = dict(delta=delta,es=es)
    if sub=="EX": kw["nex"]=8
    x,y = quiet(lambda: ed.get_dynamical_correlator(mode="ED",submode=sub,name=pair,**kw))
    return np.array(y)

ed,dense,E,V = setup(0.3)
ed.get_gs(mode="ED")
eE,wE = ed.get_excited_states(n=2,mode="ED")
ed.set_gs(wE[1])
pair = (ed.Sx[0]+1j*ed.Sy[0], ed.Sx[0]-1j*ed.Sy[0])
ref = {"|1> own E1": curve(dense,E,V,V[:,1],E[1],*pair),
       "|1> from E0": curve(dense,E,V,V[:,1],E[0],*pair),
       "|0> from E0": curve(dense,E,V,V[:,0],E[0],*pair)}
print("(a) parent, Bz=0.3, set_gs(|1>)")
for sub in ("ED","EX","CVM"):
    y = run(ed,sub,pair)
    print("   SpSm %-4s "%sub+"  ".join("max|y-%s|=%.3f"%(k,np.max(np.abs(y-r))) for k,r in ref.items()))

ed,dense,E,V = setup(0.0)
ed.get_gs(mode="ED")
Szt = dense(ed.Sz[0]+ed.Sz[1]+ed.Sz[2])
P = V[:,:2]; w,U = np.linalg.eigh(P.conj().T@Szt@P); up = P@U[:,1]
ed.set_gs(State(up,ed.get_ED_obj()))
pair = (ed.Sx[0]+1j*ed.Sy[0], ed.Sx[0]-1j*ed.Sy[0])
rm = curve(dense,E,V,up,E[0],*pair)
ra = 0.5*(curve(dense,E,V,V[:,0],E[0],*pair)+curve(dense,E,V,V[:,1],E[0],*pair))
print("(b) parent, B=0, set_gs(Sz=+1/2 member)")
for sub in ("ED","EX","CVM"):
    y = run(ed,sub,pair)
    print("   SpSm %-4s max|y-member|=%.3f max|y-manifold avg|=%.3f"%(sub,np.max(np.abs(y-rm)),np.max(np.abs(y-ra))))
```

`reviews/groundstate_1_new1/02_parent_ed_split.out`:


```
dmrgpy from <scratch>/groundstate/parent_src/src/dmrgpy/__init__.py
(a) parent, Bz=0.3, set_gs(|1>)
   SpSm ED   max|y-|1> own E1|=2.834  max|y-|1> from E0|=2.833  max|y-|0> from E0|=0.000
   SpSm EX   max|y-|1> own E1|=2.834  max|y-|1> from E0|=2.833  max|y-|0> from E0|=0.000
   SpSm CVM  max|y-|1> own E1|=2.755  max|y-|1> from E0|=0.000  max|y-|0> from E0|=2.833
(b) parent, B=0, set_gs(Sz=+1/2 member)
   SpSm ED   max|y-member|=1.417 max|y-manifold avg|=0.000
   SpSm EX   max|y-member|=0.000 max|y-manifold avg|=1.417
   SpSm CVM  max|y-member|=0.000 max|y-manifold avg|=1.417
```

```bash
cd <scratch>/reviews/groundstate_1_new1 && <scratch>/run3.sh 03_explicit_wf0.py 2>&1 | tee 03_explicit_wf0.out
```

`reviews/groundstate_1_new1/03_explicit_wf0.py`:

```python
# Reviewer probe for groundstate_1_new1, HEAD: does submode ED read an
# EXPLICIT wf0= (the named parameter of edtk/dynamics.get_dynamical_correlator,
# forwarded by Many_Body_Chain.get_dynamical_correlator's mode="ED" branch),
# with no set_gs at all?  Same 3-site Bz=0.3 chain, wf0 = |1>.
import numpy as np, io, contextlib, warnings
warnings.simplefilter("ignore")
import dmrgpy; print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain

def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, delta, B = 3, 0.05, 0.3
es = np.linspace(-0.8,2.4,641)
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
def dense(m):
    m = ed.get_ED_obj().MO2matrix(m)
    return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
E,V = np.linalg.eigh(dense(ham(ed,B)))
def curve(st,ref,A,B_):
    Am,Bm = dense(A),dense(B_)
    M = (st.conj()@Am@V)*(V.conj().T@Bm@st)
    return sum(M[k]*delta/np.pi/((es-(E[k]-ref))**2+delta**2) for k in range(len(E)))
ed.get_gs(mode="ED") # without it CVM raises TypeError from e0=None (the recorded lead)
eE,wE = ed.get_excited_states(n=2,mode="ED")
print("|<wE[1]|V1>|^2 = %.6f"%abs(np.vdot(V[:,1],wE[1].v))**2)
pair = (ed.Sx[0]+1j*ed.Sy[0], ed.Sx[0]-1j*ed.Sy[0])
ref = {"|1> own E1": curve(V[:,1],E[1],*pair),
       "|1> from E0": curve(V[:,1],E[0],*pair),
       "|0> from E0": curve(V[:,0],E[0],*pair)}
for sub in ("ED","CVM","INV","ROOTN"):
    x,y = quiet(lambda: ed.get_dynamical_correlator(mode="ED",submode=sub,
                name=pair,delta=delta,es=es,wf0=wE[1]))
    y = np.array(y)
    print("   wf0=|1> SpSm %-5s "%sub+"  ".join("max|y-%s|=%.3f"%(k,np.max(np.abs(y-r))) for k,r in ref.items()))
```

`reviews/groundstate_1_new1/03_explicit_wf0.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
|<wE[1]|V1>|^2 = 1.000000
   wf0=|1> SpSm ED    max|y-|1> own E1|=2.834  max|y-|1> from E0|=2.833  max|y-|0> from E0|=0.000
   wf0=|1> SpSm CVM   max|y-|1> own E1|=2.755  max|y-|1> from E0|=0.000  max|y-|0> from E0|=2.833
   wf0=|1> SpSm INV   max|y-|1> own E1|=2.755  max|y-|1> from E0|=0.000  max|y-|0> from E0|=2.833
   wf0=|1> SpSm ROOTN max|y-|1> own E1|=2.755  max|y-|1> from E0|=0.000  max|y-|0> from E0|=2.833
```


**Suggested fix**: of "read `wf0`" and "raise when a state has been set", the
first, since every other ED submode and all three DMRG backends now read the
state. The shape already exists in `867e2b4`'s own `dcex` fix: project the
chain's state, or the explicit `wf0`, onto the eigenbasis, c = vs^H wf0, use it as
the single initial vector with M_n built the way `dynamical_sum` does for one
row, and put the origin at the state's own energy, sum |c|^2 E, the same choice
finding 1 makes; keep the `dex` manifold average only when no state was injected
or passed, so every default number stays bit-identical. The ED side cannot tell
an injected state from its own solve today (`groundstate.set_gs`'s ED branch sets
only `ED_obj.computed_gs=True` and `ED_obj.wf0`), so it needs a mark of its own,
cleared by any solve and mirroring `groundstate.mark_injected`, plus forwarding
an explicit `wf0=` into the branch. Short of that, a `NotImplementedError` when an
explicit `wf0=` reaches `submode="ED"` would at least stop the keyword from being
dropped. NUMBERS CHANGE only after `set_gs` or with an explicit `wf0=`.

### 3. `set_gs`, `set_initial_wf` or `set_initial_wf_guess` inside a conserved sector, followed by `promote_to_dense()` with no ground-state read in between, silently discards the state on `"python"` and v3: the first reader afterwards returns the session's own sector ground state (-1.616025 against <x|H|x> = -0.957107, overlap 0.0000)

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED &middot; lens `groundstate`

**Where**: `src/dmrgpy/manybodychain.py:636` (`self.wf0 =
self.promote_mps(self.wf0)` replaces the marked object); `groundstate.py:129`
(`pending_injection`'s identity test, `mark[1] is not self.wf0`);
`pyitensor/chain.py:393-432` and `mpscpp3` `Chain::promote_to_dense` (the
session promotes and carries its own solved state, which x never reached).

`867e2b4`'s injection mark holds the injected object and counts only while
`self.wf0` is that object, and its docstring lists what retires it: `restart()`,
a solve, a backend switch. `promote_to_dense` is not on that list, but it replaces
`self.wf0` by a promoted copy, so a mark still pending when the promotion runs is
retired silently, and since x never reached the session, the session promotes its
own sector ground state and the chain adopts it. `promote_to_dense`'s docstring
reads "Leave conserved-sector mode while keeping the state computed in it". It
survived because `tests/test_sector_promotion*.py` promote only after a solve,
and the `867e2b4` tests never combine an injected state with a sector. The older
half: on the parent `e5049b0` x is lost the same way (overlap 0.0000), there
because `set_gs` never reached the session at all (the second pass's finding 11).

**Expected**: after `set_gs(x)` and `promote_to_dense()`, `gs_energy()` =
<x|H|x> = -0.957107 (the Sz=0 triplet member of the 4-site Heisenberg chain) and
`get_gs()` = x.

Repro, from the hunter (`03` with and without a read between the setter and the
promotion; `07`, spliced under finding 1, dates it):

```bash
cd <scratch>/groundstate && <scratch>/run3.sh 03_setgs_then_promote.py 2>&1 | tee 03_setgs_then_promote.out
```

`groundstate/03_setgs_then_promote.py`:

```python
# set_gs(x) inside a conserved sector, then promote_to_dense() with no read in
# between: which state does the chain hold afterwards?  x = the sector's
# first excited state (the Sz=0 triplet member of a 4-site Heisenberg chain).
import numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
n = 4
for v in ("python",3):
  for read_first in (True,False):
    np.random.seed(3)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc))
    sc.set_conserved_sector(Sz=0)
    e0 = sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=2))
    x = ww[1]
    ex = (sc.aMb(x,sc.hamiltonian,x)/sc.overlap(x,x)).real
    sc.set_gs(x)
    if read_first: sc.gs_energy()   # hand x to the session before promoting
    sc.promote_to_dense()
    b = sc.gs_energy()
    wf = sc.get_gs()
    xd = sc.promote_mps(x)
    ov = abs(sc.overlap(xd,wf))**2/(abs(sc.overlap(xd,xd))*abs(sc.overlap(wf,wf)))
    s01 = sc.vev(sc.Sz[0]*sc.Sz[1]).real
    print("v=%-6s sector E0=%.6f <x|H|x>=%.6f | set_gs(x);%s promote_to_dense(): gs_energy()=%.6f |<x|wf0>|^2=%.4f <Sz0Sz1>=%.6f"
          %(v,e0,ex," gs_energy();" if read_first else "             ",b,ov,s01))
```

`groundstate/03_setgs_then_promote.out`:


```
v=python sector E0=-1.616025 <x|H|x>=-0.957107 | set_gs(x); gs_energy(); promote_to_dense(): gs_energy()=-0.957107 |<x|wf0>|^2=1.0000 <Sz0Sz1>=-0.176777
v=python sector E0=-1.616025 <x|H|x>=-0.957107 | set_gs(x);              promote_to_dense(): gs_energy()=-1.616025 |<x|wf0>|^2=0.0000 <Sz0Sz1>=-0.227671
v=3      sector E0=-1.616025 <x|H|x>=-0.957107 | set_gs(x); gs_energy(); promote_to_dense(): gs_energy()=-0.957107 |<x|wf0>|^2=1.0000 <Sz0Sz1>=-0.176777
v=3      sector E0=-1.616025 <x|H|x>=-0.957107 | set_gs(x);              promote_to_dense(): gs_energy()=-1.616025 |<x|wf0>|^2=0.0000 <Sz0Sz1>=-0.227671
```

`groundstate/07_parent_compare.py` and its output are spliced under finding 1 above.


**Reviewer (CONFIRMED)**: reproduced exactly on both backends. With a field
-0.9*sum Sz added so that the sector ground state (-1.616025), the global one (ED
-1.857107, Sz_tot=1) and x (-0.957107) are three different energies, the loss
shows in all 18 combinations of setter (`set_gs`, `set_initial_wf`,
`set_initial_wf_guess`), first reader (`gs_energy`, `get_gs`, `vev`) and backend:
e0=-1.616025, overlap 0.0000, Sz_tot=0, so the chain returns the session's
promoted sector ground state, not an unconstrained re-solve. With a
non-eigenstate x = normalize(|1> + 0.6|0>) (<x|H|x> = -1.131526), the claim holds
too (-1.616025, overlap 0.2647 = |<x|GS>|^2), no promotion keeps x, and
re-marking `promote_mps(x)` as injected after the promotion restores it exactly
on both backends and both seeds. Nothing in `docs/known_issue_*` or `ROADMAP.md`
mentions promotion, so this is an unintended retire, not a documented boundary.
Struck, each explicitly:

- The implicit counterpart, stated in the hunter's `why_tests_pass` and in its
  ruled-out list, that a read between the setter and `promote_to_dense()` hides
  the defect: that holds only for an exact eigenstate. For a non-eigenstate x the
  read-first path loses it too, on `"python"` to the global ground state
  (-1.857107, Sz_tot=1) and on v3 to the sector one; that is finding 4.
- `groundstate.py:86-90` (`gs_is_current`, then the sector-key mismatch) as part
  of the mechanism: `mark_injected` already sets `computed_gs=False`, so
  `gs_is_current` fails at its first test and the sector key plays no part. The
  loss is entirely `pending_injection`'s identity test together with x never
  having reached the session.

There is a public workaround, justified by reading: promote first, then
`set_gs(sc.promote_mps(x))`, which is the reviewer's re-mark case.

```bash
cd <scratch>/reviews/groundstate_2 && <scratch>/run3.sh 01_hunter_repro.py 2>&1 | tee 01_hunter_repro.out
```

`reviews/groundstate_2/01_hunter_repro.py` is `groundstate/03_setgs_then_promote.py` (finding 3) unchanged, rerun; its output, `reviews/groundstate_2/01_hunter_repro.out`:


```
v=python sector E0=-1.616025 <x|H|x>=-0.957107 | set_gs(x); gs_energy(); promote_to_dense(): gs_energy()=-0.957107 |<x|wf0>|^2=1.0000 <Sz0Sz1>=-0.176777
v=python sector E0=-1.616025 <x|H|x>=-0.957107 | set_gs(x);              promote_to_dense(): gs_energy()=-1.616025 |<x|wf0>|^2=0.0000 <Sz0Sz1>=-0.227671
v=3      sector E0=-1.616025 <x|H|x>=-0.957107 | set_gs(x); gs_energy(); promote_to_dense(): gs_energy()=-0.957107 |<x|wf0>|^2=1.0000 <Sz0Sz1>=-0.176777
v=3      sector E0=-1.616025 <x|H|x>=-0.957107 | set_gs(x);              promote_to_dense(): gs_energy()=-1.616025 |<x|wf0>|^2=0.0000 <Sz0Sz1>=-0.227671
```

```bash
cd <scratch>/reviews/groundstate_2 && <scratch>/run3.sh 02_field_and_readers.py 2>&1 | tee 02_field_and_readers.out
```

`reviews/groundstate_2/02_field_and_readers.py`:

```python
# Attack groundstate_2.  A field B*sum Sz puts the global ground state in
# Sz!=0, so the sector ground state, the global ground state and x (the
# sector's first excited state) are three different energies.  For each
# setter (set_gs, set_initial_wf, set_initial_wf_guess) and each first reader
# after promote_to_dense() (gs_energy, get_gs, vev), which state does the
# chain hold?  Control: the same with a read before promoting.
import numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
n,B = 4,0.9
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h - B*sc.Sz[i]
    return h
# exact references
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed))
print("ED global E0 = %.6f"%ed.gs_energy(mode="ED"))
def fresh(v,seed):
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc)); sc.set_conserved_sector(Sz=0)
    return sc
for v in ("python",3):
  for setter in ("set_gs","set_initial_wf","set_initial_wf_guess"):
    for reader in ("gs_energy","get_gs","vev"):
      for read_first in (False,True):
        sc = fresh(v,3)
        e0 = sc.gs_energy()
        ee,ww = quiet(lambda: sc.get_excited_states(n=2)); x = ww[1]
        ex = (sc.aMb(x,sc.hamiltonian,x)/sc.overlap(x,x)).real
        getattr(sc,setter)(x)
        if read_first: sc.gs_energy()
        sc.promote_to_dense()
        if reader=="gs_energy": sc.gs_energy()
        elif reader=="get_gs": sc.get_gs()
        else: sc.vev(sc.Sz[0])
        wf = sc.wf0
        xd = sc.promote_mps(x)
        ov = abs(sc.overlap(xd,wf))**2/(abs(sc.overlap(xd,xd))*abs(sc.overlap(wf,wf)))
        eh = (sc.aMb(wf,sc.hamiltonian,wf)/sc.overlap(wf,wf)).real
        szt = sum(sc.vev(sc.Sz[i]).real for i in range(n))
        print("v=%-6s %-20s read_first=%-5s reader=%-9s sectorE0=%.6f <x|H|x>=%.6f | e0=%.6f <wf|H|wf>=%.6f |<x|wf>|^2=%.4f Sz_tot=%.4f"
              %(v,setter,read_first,reader,e0,ex,sc.e0.real if sc.e0 is not None else np.nan,eh,ov,szt))
```

`reviews/groundstate_2/02_field_and_readers.out`:


```
ED global E0 = -1.857107
v=python set_gs               read_first=False reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_gs               read_first=True  reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=python set_gs               read_first=False reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_gs               read_first=True  reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=python set_gs               read_first=False reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_gs               read_first=True  reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=python set_initial_wf       read_first=False reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_initial_wf       read_first=True  reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=python set_initial_wf       read_first=False reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_initial_wf       read_first=True  reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=python set_initial_wf       read_first=False reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_initial_wf       read_first=True  reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=python set_initial_wf_guess read_first=False reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_initial_wf_guess read_first=True  reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=python set_initial_wf_guess read_first=False reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_initial_wf_guess read_first=True  reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=python set_initial_wf_guess read_first=False reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=python set_initial_wf_guess read_first=True  reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=3      set_gs               read_first=False reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=0.0000
v=3      set_gs               read_first=True  reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      set_gs               read_first=False reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=3      set_gs               read_first=True  reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      set_gs               read_first=False reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=3      set_gs               read_first=True  reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      set_initial_wf       read_first=False reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=-0.0000
v=3      set_initial_wf       read_first=True  reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=3      set_initial_wf       read_first=False reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=0.0000
v=3      set_initial_wf       read_first=True  reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      set_initial_wf       read_first=False reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=0.0000
v=3      set_initial_wf       read_first=True  reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=3      set_initial_wf_guess read_first=False reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=0.0000
v=3      set_initial_wf_guess read_first=True  reader=gs_energy sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=3      set_initial_wf_guess read_first=False reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=0.0000
v=3      set_initial_wf_guess read_first=True  reader=get_gs    sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=3      set_initial_wf_guess read_first=False reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.0000 Sz_tot=0.0000
v=3      set_initial_wf_guess read_first=True  reader=vev       sectorE0=-1.616025 <x|H|x>=-0.957107 | e0=-0.957107 <wf|H|wf>=-0.957107 |<x|wf>|^2=1.0000 Sz_tot=0.0000
```

```bash
cd <scratch>/reviews/groundstate_2 && <scratch>/run3.sh 03_non_eigenstate_and_fixes.py 2>&1 | tee 03_non_eigenstate_and_fixes.out
```

`reviews/groundstate_2/03_non_eigenstate_and_fixes.py`:

```python
# Attack the control and the suggested fixes of groundstate_2 with an x that
# is NOT an eigenstate: x = normalize(|1> + 0.6|0>) inside Sz=0, where |0>,|1>
# are the sector's two lowest states (the triplet |1> of the hunter's repro
# is an exact eigenstate, which a DMRG sweep from it cannot leave except by
# roundoff, so the hunter's "read first" control could not tell "x kept"
# from "x swept and stayed").  Field B*sum Sz: global GS at Sz_tot=1.
#   noprom    : set_gs(x); gs_energy()                   (in-sector reference)
#   hunter    : set_gs(x); promote; gs_energy()          (the claim)
#   readfirst : set_gs(x); gs_energy(); promote; gs_energy()  (= fix "take the injection first")
#   remark    : set_gs(x); promote; re-mark promote_mps(x) as injected; gs_energy()  (= fix "re-mark")
import time, numpy as np, io, contextlib
from dmrgpy import spinchain, groundstate
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
n,B = 4,0.9
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h - B*sc.Sz[i]
    return h
t0 = time.time()
for v in ("python",3):
  for case in ("noprom","hunter","readfirst","remark"):
    for seed in (3,4):
      np.random.seed(seed)
      sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
      sc.maxm=20; sc.nsweeps=12
      sc.set_hamiltonian(ham(sc)); sc.set_conserved_sector(Sz=0)
      e0 = sc.gs_energy()
      ee,ww = quiet(lambda: sc.get_excited_states(n=2))
      x = (ww[1] + 0.6*ww[0]).normalize()
      ex = (sc.aMb(x,sc.hamiltonian,x)/sc.overlap(x,x)).real
      sc.set_gs(x)
      if case=="noprom":
          sc.gs_energy()
      elif case=="hunter":
          sc.promote_to_dense(); sc.gs_energy()
      elif case=="readfirst":
          sc.gs_energy(); sc.promote_to_dense(); sc.gs_energy()
      elif case=="remark":
          sc.promote_to_dense()
          groundstate.mark_injected(sc, sc.promote_mps(x), reconverge=False)
          sc.gs_energy()
      wf = sc.wf0
      xd = x if case=="noprom" else sc.promote_mps(x)
      ov = abs(sc.overlap(xd,wf))**2/(abs(sc.overlap(xd,xd))*abs(sc.overlap(wf,wf)))
      eh = (sc.aMb(wf,sc.hamiltonian,wf)/sc.overlap(wf,wf)).real
      szt = sum(sc.vev(sc.Sz[i]).real for i in range(n))
      print("v=%-6s %-9s seed=%d sectorE0=%.6f <x|H|x>=%.6f | e0=%.6f <wf|H|wf>=%.6f |<x|wf>|^2=%.4f Sz_tot=%.4f"
            %(v,case,seed,e0,ex,float(np.real(sc.e0)),eh,ov,szt), flush=True)
print("elapsed %.1fs"%(time.time()-t0))
```

`reviews/groundstate_2/03_non_eigenstate_and_fixes.out`:


```
v=python noprom    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=python noprom    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=python hunter    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
v=python hunter    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
v=python readfirst seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=python readfirst seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=python remark    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=python remark    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      noprom    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      noprom    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      hunter    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=0.0000
v=3      hunter    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=0.0000
v=3      readfirst seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=0.0000
v=3      readfirst seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
v=3      remark    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=3      remark    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
elapsed 7.7s
```

```bash
cd <scratch>/reviews/groundstate_2 && <scratch>/run3.sh 04_v3_resweep_after_promote.py 2>&1 | tee 04_v3_resweep_after_promote.out
```

`reviews/groundstate_2/04_v3_resweep_after_promote.py`:

```python
# Side probe: the 2026-09 record's finding 10 says v3's post-promotion
# re-sweep "cannot leave that sector" because H commutes with the charge.
# Probe 02 saw v3 land on the global Sz_tot=1 ground state once after
# promotion.  How often does v3 leave, (a) on the plain documented workflow
# (sector solve; promote_to_dense(); gs_energy()), (b) after a sector sweep
# started from the triplet x (set_initial_wf_guess(x); gs_energy(); promote;
# get_gs())?  Also "python" for (a), which has the _promoted_gs bridge.
import time, numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
n,B = 4,0.9
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h - B*sc.Sz[i]
    return h
t0 = time.time()
for v,case,seeds in ((3,"plain",range(12)),("python","plain",range(6)),(3,"guess_x",range(12))):
  res = []
  for seed in seeds:
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc)); sc.set_conserved_sector(Sz=0)
    e0 = sc.gs_energy()
    if case=="guess_x":
      ee,ww = quiet(lambda: sc.get_excited_states(n=2)); x = ww[1]
      sc.set_initial_wf_guess(x); sc.gs_energy()
    eb = float(np.real(sc.e0))
    sc.promote_to_dense()
    ea = sc.gs_energy()
    szt = sum(sc.vev(sc.Sz[i]).real for i in range(n))
    res.append((seed,eb,ea,szt))
    print("v=%-6s %-8s seed=%2d before promote e0=%.6f  after promote gs_energy()=%.6f Sz_tot=%.4f"%(v,case,seed,eb,ea,szt),flush=True)
  left = sum(1 for r in res if abs(r[3])>1e-3)
  print("  -> v=%s %s: left Sz=0 in %d of %d runs"%(v,case,left,len(res)))
print("elapsed %.1fs"%(time.time()-t0))
```

`reviews/groundstate_2/04_v3_resweep_after_promote.out`:


```
v=3      plain    seed= 0 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      plain    seed= 1 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      plain    seed= 2 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 3 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 4 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      plain    seed= 5 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 6 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 7 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 8 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 9 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed=10 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed=11 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
  -> v=3 plain: left Sz=0 in 0 of 12 runs
v=python plain    seed= 0 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 1 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 2 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 3 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 4 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 5 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
  -> v=python plain: left Sz=0 in 0 of 6 runs
v=3      guess_x  seed= 0 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=0.0000
v=3      guess_x  seed= 1 before promote e0=-0.957107  after promote gs_energy()=-1.857097 Sz_tot=1.0000
v=3      guess_x  seed= 2 before promote e0=-0.957107  after promote gs_energy()=-1.032239 Sz_tot=1.0040
v=3      guess_x  seed= 3 before promote e0=-0.957107  after promote gs_energy()=-1.857107 Sz_tot=1.0000
v=3      guess_x  seed= 4 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      guess_x  seed= 5 before promote e0=-0.957107  after promote gs_energy()=-1.857107 Sz_tot=1.0000
v=3      guess_x  seed= 6 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      guess_x  seed= 7 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=-0.0000
v=3      guess_x  seed= 8 before promote e0=-0.957107  after promote gs_energy()=-1.857107 Sz_tot=1.0000
v=3      guess_x  seed= 9 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=-0.0000
v=3      guess_x  seed=10 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=0.0000
v=3      guess_x  seed=11 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
  -> v=3 guess_x: left Sz=0 in 5 of 12 runs
elapsed 8.5s
```


**Suggested fix**: in `Many_Body_Chain.promote_to_dense`, read `pending =
groundstate.pending_injection(self)` before promoting, and after `self.wf0 =
self.promote_mps(self.wf0)` re-mark the promoted copy with the same mode (skip or
reconverge), so the next read hands the promoted x to the session. The hunter's
"simpler alternative", calling `get_gs()` at the top of `promote_to_dense`, is
wrong and should not be taken: it is exactly the read-first path of finding 4,
where the session holds x with no energy and re-sweeps it after the promotion.
The reviewer's wider version covers this finding and finding 4 together: whenever
the chain holds a state it considers its own (`gs_is_current`, or a pending
mark), re-mark the promoted `wf0` as injected, with the pending mode if one was
pending and skip otherwise, so the next read goes through `_take_injected_state`
and takes it unswept with e0 = <wf|H|wf> on both backends. The regression should
use a non-eigenstate x and a field that separates the sector ground state from
the global one; with the triplet and no field neither a re-sweep nor a sector
escape can show. NUMBERS CHANGE only for a setter followed directly by
`promote_to_dense`: `gs_energy()` from the sector ground energy to <x|H|x>
(-1.616025 to -0.957107 here) and every observable afterwards from the sector
ground state to x (<Sz0Sz1> -0.227671 to -0.176777).

### 4. After `promote_to_dense()`, the next ground-state read re-sweeps whatever state the chain carries instead of keeping it, on `"python"` whenever the session holds no cached energy for it (every injected state) and on v3 always, so a carried state that is not the sector ground state is replaced, and on `"python"` one far from it leaves the sector for the global ground state (-1.857107 at Sz_tot=1 against the kept -1.131526)

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `groundstate` (turned up by the reviewer of finding 3)

**Where**: `src/dmrgpy/groundstate.py::_take_injected_state` (the session
receives the state through `set_wavefunction`, which drops its energy:
`pyitensor/chain.py:643-646` sets `_wf0_energy=None`), so
`pyitensor/chain.py::promote_to_dense`'s `_promoted_gs` bridge carries
`energy=None` and `gs_energy(skip_dmrg=True)` at `pyitensor/chain.py:613` falls
through to a dense sweep; on v3, `mpscpp3/chain_session.h:505-520`
(`promote_to_dense`) carries `have_wf0_energy_`, but `set_hamiltonian_mpo` at
`:567` clears it on the re-send that `manybodychain.py`'s `promote_to_dense`
forces (`self._session_ham_cache = None`), so that carry is dead code on v3.

`promote_to_dense` promises to keep the state computed in the sector, and it does
for the one state that cannot tell: the sector ground state, which a re-sweep on
the dense sites lands back on. Any other carried state is replaced by whatever the
unconstrained sweep reaches. On `"python"` the trigger is a state the session holds
without a cached energy, which since `867e2b4` is every state taken unswept
through `set_gs`, `set_initial_wf` or `gs_energy(wf0=x, reconverge=False)`; on v3
it is every state, because the 2026-09 audit's finding 10 fixed this re-sweep on
`"python"` only. The plain solve-then-promote workflow is unaffected, which is why
the existing `test_promote_to_dense_keeps_the_sector_energy`, which carries the
sector ground state, cannot see it; a KPM `Sx` correlator after the promotion,
the call promotion exists for, then measures the replaced state.

**Expected**: the carried state, kept: after `set_gs(x); gs_energy();
promote_to_dense(); gs_energy()` with x not an eigenstate (<x|H|x> = -1.131526) on
a 4-site Heisenberg chain in a field -0.9*sum Sz at Sz=0, the last read returns
-1.131526 with overlap 1.0000.

Repro, from the reviewer of finding 3, who found it (both scripts are spliced
under finding 3: the `readfirst` rows of `03_non_eigenstate_and_fixes` and the
seed counts of `04_v3_resweep_after_promote`):

`reviews/groundstate_2/03_non_eigenstate_and_fixes.py` and its output are spliced under finding 3 above.

`reviews/groundstate_2/04_v3_resweep_after_promote.py` and its output are spliced under finding 3 above.


**Reviewer (CONFIRMED, NARROWED)**: both scripts rerun digit for digit (the v3
seed counts differ, since numpy's seed does not reach v3's C++ RNG): in the
`readfirst` rows `"python"` gives -1.857107 with overlap 0.0000 and Sz_tot=1.0000,
v3 gives -1.616025 with overlap 0.2647, and the `noprom` and `remark` rows keep x
exactly on both. The noise term is not the cause (v3 escapes at `noise=0.0` too,
every escape starting from a sweep trapped at the excited eigenstate -0.957107),
and a realistic restore does not trigger it: saving the converged sector ground
state and restoring it with `set_gs` or `set_initial_wf`, then read, promote,
read, stays at -1.616025 in Sz=0 in 6 of 6 runs on both backends. How far the
state has to be for `"python"` to leave the sector: with x = |0> + eps|1>,
`"python"` returns the sector ground state at eps = 0.01 and 0.1 and the global
one, Sz_tot=1, at eps = 0.6 and 3.0, while v3 returns the sector ground state at
every eps with overlap 1/(1+eps^2); so the state is lost for every eps > 0 on both
backends, and leaving the sector is the `"python"`-only consequence at larger eps,
which a KPM `Sx` correlator as the reader then measures. Dating, on a read-only
archive of `867e2b4^`: the parent's read after `set_gs(x)` already gives the
sector ground state rather than x and stays there after the promotion, while the
audited tree gives x and then the global ground state, so the out-of-sector
result on `"python"` came with `867e2b4`. Struck, each explicitly:

- That the promoted state "leaves the sector" whenever it is not the sector ground
  state: being replaced is universal, leaving is conditional (`"python"` far from
  the sector ground state; v3 only from a sweep trapped at an excited
  eigenstate).
- The rate "5 of 12" for v3 `guess_x`: it is run-to-run, 1 to 5 of 12 over four
  runs.
- That the v3 re-sweep after promotion is new: it is the 2026-09 finding 10
  mechanism, fixed on `"python"` only; what is new is its consequence, reachable
  only since `867e2b4` made injected states and warm starts reach the session.

This contradicts the 2026-09 record's finding 10 statement that v3 "cannot leave
that sector", which is true only for the sector ground state.

```bash
cd <scratch>/reviews/groundstate_2_new1 && <scratch>/run3.sh 01_hunter_03_rerun.py 2>&1 | tee 01_hunter_03_rerun.out
```

`reviews/groundstate_2_new1/01_hunter_03_rerun.py` is `reviews/groundstate_2/03_non_eigenstate_and_fixes.py` (finding 3) unchanged, rerun; its output, `reviews/groundstate_2_new1/01_hunter_03_rerun.out`:


```
v=python noprom    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=python noprom    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=python hunter    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
v=python hunter    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
v=python readfirst seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=python readfirst seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=python remark    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=python remark    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      noprom    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=3      noprom    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      hunter    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=0.0000
v=3      hunter    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=0.0000
v=3      readfirst seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
v=3      readfirst seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=0.0000
v=3      remark    seed=3 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=3      remark    seed=4 sectorE0=-1.616025 <x|H|x>=-1.131526 | e0=-1.131526 <wf|H|wf>=-1.131526 |<x|wf>|^2=1.0000 Sz_tot=0.0000
elapsed 6.5s
```

```bash
cd <scratch>/reviews/groundstate_2_new1 && <scratch>/run3.sh 02_hunter_04_rerun.py 2>&1 | tee 02_hunter_04_rerun.out
```

`reviews/groundstate_2_new1/02_hunter_04_rerun.py` is `reviews/groundstate_2/04_v3_resweep_after_promote.py` (finding 3) unchanged, rerun; its output, `reviews/groundstate_2_new1/02_hunter_04_rerun.out`:


```
v=3      plain    seed= 0 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      plain    seed= 1 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 2 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      plain    seed= 3 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 4 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 5 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 6 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed= 7 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      plain    seed= 8 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      plain    seed= 9 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      plain    seed=10 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      plain    seed=11 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
  -> v=3 plain: left Sz=0 in 0 of 12 runs
v=python plain    seed= 0 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 1 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 2 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 3 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 4 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=python plain    seed= 5 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
  -> v=python plain: left Sz=0 in 0 of 6 runs
v=3      guess_x  seed= 0 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=-0.0000
v=3      guess_x  seed= 1 before promote e0=-0.957107  after promote gs_energy()=-1.857107 Sz_tot=1.0000
v=3      guess_x  seed= 2 before promote e0=-1.593878  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      guess_x  seed= 3 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=-0.0000
v=3      guess_x  seed= 4 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=0.0000
v=3      guess_x  seed= 5 before promote e0=-1.616025  after promote gs_energy()=-1.616025 Sz_tot=0.0000
v=3      guess_x  seed= 6 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=0.0000
v=3      guess_x  seed= 7 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=0.0000
v=3      guess_x  seed= 8 before promote e0=-0.957107  after promote gs_energy()=-1.857107 Sz_tot=1.0000
v=3      guess_x  seed= 9 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=0.0000
v=3      guess_x  seed=10 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=0.0000
v=3      guess_x  seed=11 before promote e0=-0.957107  after promote gs_energy()=-0.957107 Sz_tot=0.0000
  -> v=3 guess_x: left Sz=0 in 2 of 12 runs
elapsed 6.0s
```

```bash
cd <scratch>/reviews/groundstate_2_new1 && <scratch>/run3.sh 03_realistic_and_noise.py 2>&1 | tee 03_realistic_and_noise.out
```

`reviews/groundstate_2_new1/03_realistic_and_noise.py`:

```python
# Reviewer probe for groundstate_2_new1.
# (a) A realistic carried state: the chain's own converged sector ground
#     state, saved and restored with set_gs()/set_initial_wf(), then read,
#     promoted and read again.  Does "python" leave the sector here too, or
#     only for a hand-built non-eigenstate?
# (b) v3's escape after set_initial_wf_guess(x_triplet): is it the DMRG
#     noise term?  Default noise (1e-7) against noise=0, 12 runs each.
import time, numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
n,B = 4,0.9
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h - B*sc.Sz[i]
    return h
def build(v,seed,noise=None):
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    if noise is not None: sc.noise = noise
    sc.set_hamiltonian(ham(sc)); sc.set_conserved_sector(Sz=0)
    return sc
def szt(sc): return sum(sc.vev(sc.Sz[i]).real for i in range(n))
t0 = time.time()
print("(a) restore the sector ground state, read, promote, read")
for v in ("python",3):
  for setter in ("set_gs","set_initial_wf"):
    left = 0
    for seed in range(6):
      sc = build(v,seed)
      e0 = sc.gs_energy()
      x = sc.get_gs().copy()
      getattr(sc,setter)(x)
      e1 = sc.gs_energy()          # the read that takes the injection
      sc.promote_to_dense()
      e2 = sc.gs_energy()
      s = szt(sc); left += abs(s)>1e-3
      print("v=%-6s %-15s seed=%d sectorE0=%.6f read=%.6f after promote=%.6f Sz_tot=%.4f"
            %(v,setter,seed,e0,e1,e2,s),flush=True)
    print("  -> v=%s %s: left Sz=0 in %d of 6"%(v,setter,left))
print("(b) v3 set_initial_wf_guess(triplet); gs_energy(); promote; gs_energy()")
for noise in (None,0.0):
  left = 0; trapped = 0
  for seed in range(12):
    sc = build(3,seed,noise)
    sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=2)); x = ww[1]
    sc.set_initial_wf_guess(x); eb = sc.gs_energy()
    trapped += abs(eb+0.957107)<1e-5
    sc.promote_to_dense(); ea = sc.gs_energy(); s = szt(sc)
    left += abs(s)>1e-3
    print("noise=%s run=%2d before=%.6f after=%.6f Sz_tot=%.4f"%(sc.noise,seed,eb,ea,s),flush=True)
  print("  -> noise=%s: trapped at -0.957107 in %d of 12, left Sz=0 in %d of 12"
        %("default" if noise is None else noise,trapped,left))
print("elapsed %.1fs"%(time.time()-t0))
```

`reviews/groundstate_2_new1/03_realistic_and_noise.out`:


```
(a) restore the sector ground state, read, promote, read
v=python set_gs          seed=0 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_gs          seed=1 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_gs          seed=2 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_gs          seed=3 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_gs          seed=4 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_gs          seed=5 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
  -> v=python set_gs: left Sz=0 in 0 of 6
v=python set_initial_wf  seed=0 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_initial_wf  seed=1 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_initial_wf  seed=2 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_initial_wf  seed=3 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_initial_wf  seed=4 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=python set_initial_wf  seed=5 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
  -> v=python set_initial_wf: left Sz=0 in 0 of 6
v=3      set_gs          seed=0 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
v=3      set_gs          seed=1 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
v=3      set_gs          seed=2 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
v=3      set_gs          seed=3 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=3      set_gs          seed=4 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
v=3      set_gs          seed=5 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
  -> v=3 set_gs: left Sz=0 in 0 of 6
v=3      set_initial_wf  seed=0 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
v=3      set_initial_wf  seed=1 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=-0.0000
v=3      set_initial_wf  seed=2 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
v=3      set_initial_wf  seed=3 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
v=3      set_initial_wf  seed=4 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
v=3      set_initial_wf  seed=5 sectorE0=-1.616025 read=-1.616025 after promote=-1.616025 Sz_tot=0.0000
  -> v=3 set_initial_wf: left Sz=0 in 0 of 6
(b) v3 set_initial_wf_guess(triplet); gs_energy(); promote; gs_energy()
noise=1e-07 run= 0 before=-0.957107 after=-1.857107 Sz_tot=1.0000
noise=1e-07 run= 1 before=-1.616025 after=-1.616025 Sz_tot=0.0000
noise=1e-07 run= 2 before=-1.616025 after=-1.616025 Sz_tot=-0.0000
noise=1e-07 run= 3 before=-1.616025 after=-1.616025 Sz_tot=0.0000
noise=1e-07 run= 4 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
noise=1e-07 run= 5 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
noise=1e-07 run= 6 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
noise=1e-07 run= 7 before=-0.957107 after=-0.957107 Sz_tot=0.0000
noise=1e-07 run= 8 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
noise=1e-07 run= 9 before=-1.616025 after=-1.616025 Sz_tot=0.0000
noise=1e-07 run=10 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
noise=1e-07 run=11 before=-0.957107 after=-1.857097 Sz_tot=1.0000
  -> noise=default: trapped at -0.957107 in 8 of 12, left Sz=0 in 2 of 12
noise=0.0 run= 0 before=-0.957107 after=-1.857107 Sz_tot=1.0000
noise=0.0 run= 1 before=-0.957107 after=-0.957107 Sz_tot=0.0000
noise=0.0 run= 2 before=-0.957107 after=-0.957107 Sz_tot=0.0000
noise=0.0 run= 3 before=-0.957107 after=-0.957107 Sz_tot=0.0000
noise=0.0 run= 4 before=-1.616025 after=-1.616025 Sz_tot=0.0000
noise=0.0 run= 5 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
noise=0.0 run= 6 before=-0.957107 after=-0.957107 Sz_tot=0.0000
noise=0.0 run= 7 before=-0.957107 after=-0.957107 Sz_tot=0.0000
noise=0.0 run= 8 before=-0.957107 after=-0.957107 Sz_tot=0.0000
noise=0.0 run= 9 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
noise=0.0 run=10 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
noise=0.0 run=11 before=-0.957107 after=-0.957107 Sz_tot=-0.0000
  -> noise=0.0: trapped at -0.957107 in 11 of 12, left Sz=0 in 1 of 12
elapsed 11.3s
```

```bash
cd <scratch>/reviews/groundstate_2_new1 && <scratch>/run3.sh 04_origin_parent.py parent 2>&1 | tee 04_origin_parent.out
```

`reviews/groundstate_2_new1/04_origin_parent.py`:

```python
# Origin check for the "python" half: the readfirst sequence
#   set_gs(x); gs_energy(); promote_to_dense(); gs_energy()
# on the parent of 867e2b4 (git archive 867e2b4^, extracted next to this
# script, "python" backend only since the archive has no compiled .so),
# then on the audited tree (whatever PYTHONPATH the wrapper sets), in two
# subprocess-free passes selected by argv.
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
which = sys.argv[1] if len(sys.argv)>1 else "parent"
if which=="parent": sys.path.insert(0, os.path.join(here,"parent","src"))
import numpy as np, io, contextlib
import dmrgpy
print("tree:", which, "->", os.path.dirname(dmrgpy.__file__))
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
n,B = 4,0.9
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h - B*sc.Sz[i]
    return h
for case in ("readfirst","noread"):
  for seed in (3,4):
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc)); sc.set_conserved_sector(Sz=0)
    e0 = sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=2))
    x = (ww[1] + 0.6*ww[0]).normalize()
    ex = (sc.aMb(x,sc.hamiltonian,x)/sc.overlap(x,x)).real
    sc.set_gs(x)
    e1 = sc.gs_energy() if case=="readfirst" else float("nan")
    sc.promote_to_dense()
    e2 = sc.gs_energy()
    wf = sc.get_gs(); xd = sc.promote_mps(x)
    ov = abs(sc.overlap(xd,wf))**2/(abs(sc.overlap(xd,xd))*abs(sc.overlap(wf,wf)))
    s = sum(sc.vev(sc.Sz[i]).real for i in range(n))
    print("%-6s %-9s seed=%d <x|H|x>=%.6f read=%.6f after promote=%.6f |<x|wf>|^2=%.4f Sz_tot=%.4f"
          %(which,case,seed,ex,e1,e2,ov,s),flush=True)
```

`reviews/groundstate_2_new1/04_origin_parent.out`:


```
tree: parent -> <scratch>/reviews/groundstate_2_new1/parent/src/dmrgpy
parent readfirst seed=3 <x|H|x>=-1.131526 read=-1.616025 after promote=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
parent readfirst seed=4 <x|H|x>=-1.131526 read=-1.616025 after promote=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
parent noread    seed=3 <x|H|x>=-1.131526 read=nan after promote=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
parent noread    seed=4 <x|H|x>=-1.131526 read=nan after promote=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
```

```bash
cd <scratch>/reviews/groundstate_2_new1 && <scratch>/run3.sh 05_origin_current.py current 2>&1 | tee 05_origin_current.out
```

`reviews/groundstate_2_new1/05_origin_current.py` is `reviews/groundstate_2_new1/04_origin_parent.py` (finding 4) unchanged, run with the argument `current`; its output, `reviews/groundstate_2_new1/05_origin_current.out`:


```
tree: current -> <repo>/src/dmrgpy
current readfirst seed=3 <x|H|x>=-1.131526 read=-1.131526 after promote=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
current readfirst seed=4 <x|H|x>=-1.131526 read=-1.131526 after promote=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
current noread    seed=3 <x|H|x>=-1.131526 read=nan after promote=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
current noread    seed=4 <x|H|x>=-1.131526 read=nan after promote=-1.616025 |<x|wf>|^2=0.2647 Sz_tot=-0.0000
```

```bash
cd <scratch>/reviews/groundstate_2_new1 && <scratch>/run3.sh 06_python_eps_and_reader.py 2>&1 | tee 06_python_eps_and_reader.out
```

`reviews/groundstate_2_new1/06_python_eps_and_reader.py`:

```python
# How far from an eigenstate must the carried state be for "python" to
# leave the sector?  x = normalize(|0> + eps|1>), |0>,|1> the sector's two
# lowest states; set_gs(x); gs_energy(); promote_to_dense(); then a reader.
# Readers: gs_energy(), and a KPM correlator of the charge-changing Sx
# (the call promotion exists for), after which the chain's state is read.
import numpy as np, io, contextlib, time
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
n,B = 4,0.9
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h - B*sc.Sz[i]
    return h
t0 = time.time()
for v in ("python",3):
 for reader in ("gs_energy","kpm_Sx"):
  for eps in (0.0,0.01,0.1,0.6,3.0):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc)); sc.set_conserved_sector(Sz=0)
    sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=2))
    x = (ww[0] + eps*ww[1]).normalize()
    ex = (sc.aMb(x,sc.hamiltonian,x)/sc.overlap(x,x)).real
    sc.set_gs(x); sc.gs_energy()
    sc.promote_to_dense()
    if reader=="gs_energy": sc.gs_energy()
    else: quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sx[0],sc.Sx[0]),
                                                   es=np.linspace(-1,4,50),delta=0.2))
    wf = sc.wf0; xd = sc.promote_mps(x)
    ov = abs(sc.overlap(xd,wf))**2/(abs(sc.overlap(xd,xd))*abs(sc.overlap(wf,wf)))
    eh = (sc.aMb(wf,sc.hamiltonian,wf)/sc.overlap(wf,wf)).real
    s = sum(sc.vev(sc.Sz[i]).real for i in range(n))
    print("v=%-6s reader=%-9s eps=%-4s <x|H|x>=%.6f | chain e0=%.6f <wf|H|wf>=%.6f |<x|wf>|^2=%.4f Sz_tot=%.4f"
          %(v,reader,eps,ex,float(np.real(sc.e0)),eh,ov,s),flush=True)
print("elapsed %.1fs"%(time.time()-t0))
```

`reviews/groundstate_2_new1/06_python_eps_and_reader.out`:


```
v=python reader=gs_energy eps=0.0  <x|H|x>=-1.616025 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=python reader=gs_energy eps=0.01 <x|H|x>=-1.615960 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.9999 Sz_tot=0.0000
v=python reader=gs_energy eps=0.1  <x|H|x>=-1.609501 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.9901 Sz_tot=0.0000
v=python reader=gs_energy eps=0.6  <x|H|x>=-1.441606 | chain e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=python reader=gs_energy eps=3.0  <x|H|x>=-1.022999 | chain e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=python reader=kpm_Sx    eps=0.0  <x|H|x>=-1.616025 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=1.0000 Sz_tot=0.0000
v=python reader=kpm_Sx    eps=0.01 <x|H|x>=-1.615960 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.9999 Sz_tot=0.0000
v=python reader=kpm_Sx    eps=0.1  <x|H|x>=-1.609501 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.9901 Sz_tot=0.0000
v=python reader=kpm_Sx    eps=0.6  <x|H|x>=-1.441606 | chain e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=python reader=kpm_Sx    eps=3.0  <x|H|x>=-1.022999 | chain e0=-1.857107 <wf|H|wf>=-1.857107 |<x|wf>|^2=0.0000 Sz_tot=1.0000
v=3      reader=gs_energy eps=0.0  <x|H|x>=-1.616025 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=3      reader=gs_energy eps=0.01 <x|H|x>=-1.615960 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.9999 Sz_tot=0.0000
v=3      reader=gs_energy eps=0.1  <x|H|x>=-1.609501 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.9901 Sz_tot=-0.0000
v=3      reader=gs_energy eps=0.6  <x|H|x>=-1.441606 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.7353 Sz_tot=0.0000
v=3      reader=gs_energy eps=3.0  <x|H|x>=-1.022999 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.1000 Sz_tot=0.0000
v=3      reader=kpm_Sx    eps=0.0  <x|H|x>=-1.616025 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=1.0000 Sz_tot=-0.0000
v=3      reader=kpm_Sx    eps=0.01 <x|H|x>=-1.615960 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.9999 Sz_tot=-0.0000
v=3      reader=kpm_Sx    eps=0.1  <x|H|x>=-1.609501 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.9901 Sz_tot=-0.0000
v=3      reader=kpm_Sx    eps=0.6  <x|H|x>=-1.441606 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.7353 Sz_tot=0.0000
v=3      reader=kpm_Sx    eps=3.0  <x|H|x>=-1.022999 | chain e0=-1.616025 <wf|H|wf>=-1.616025 |<x|wf>|^2=0.1000 Sz_tot=0.0000
elapsed 9.8s
```


**Suggested fix**: finding 3's re-mark, in its wider form. At the end of
`Many_Body_Chain.promote_to_dense`, after `self.wf0 = self.promote_mps(self.wf0)`,
call `groundstate.mark_injected(self, self.wf0, reconverge=False)`, so the next
read takes the promoted state unswept with e0 = <wf|H|wf> on every backend
(-1.131526, overlap 1.0000 on both, measured). It needs no rebuild, fixes v3,
which the 2026-09 finding 10 fix never reached, and makes `pyitensor`'s
`_promoted_gs` bridge redundant. Three conditions: mark only when the state was
current before the session promote (`gs_is_current` held), so a stale state is
not frozen as current; mark when an injection was pending too, keeping its mode
(a pending `set_initial_wf_guess` stays "reconverge"), which is finding 3; and
extend `test_promote_to_dense_keeps_the_sector_energy` with a `set_gs` of a
non-lowest state on both backends. A `"python"`-only patch (computing <H> in
`pyitensor`'s `promote_to_dense` when `_wf0_energy` is None) would leave v3 losing
the state and should be rejected. The hunter's `03` labels `readfirst` as the fix
"take the injection first"; it is the failing path. NUMBERS CHANGE only when the
carried state is not the sector ground state: on the chain above, from -1.857107
(`"python"`) or -1.616025 (v3) to -1.131526.

### 5. After `set_hamiltonian(H2, restart=False)` on a solved chain, `gs_is_current` still holds, so `gs_energy()`, `get_gs()`/`vev()`, `get_excited()` and `get_dynamical_correlator_moments()` keep answering for H1 on v2, v3 and `"python"` (-1.616025 against an exact -1.780099, 9.2 per cent) while the public correlator alone re-solves, and `gs_energy()` on one chain reads either number depending on whether a correlator ran in between

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED &middot; lenses `groundstate` and `kpm` (found by both; the `kpm` half turned up by the reviewer of finding 11)

**Where**: `src/dmrgpy/groundstate.py:43-52` (`solver_key` omits the
Hamiltonian, and its docstring rests on the premise that changing it "goes
through `set_hamiltonian()`, which resets `computed_gs` outright", false under
`restart=False`) and `:55-90` (`gs_is_current`); `manybodychain.py:1235`
(`gs_energy` returns `self.e0`) and `:1212` (`get_gs` returns `self.wf0`);
`vev.py:30`; `excited.py:47` (`get_gs`, then `session.excited_states` on the
session still holding H1, then rediagonalized with H2 in H1's states);
`kpmdmrg.py:130` (`dynamical_correlator_moments` reads `get_gs()` and the
session's H), reached from `get_dynamical_correlator_moments` and, by reading,
`fermionchain.get_gr`, which also adds the stale `self.e0`; `groundstate.py:
301-305`, the one reader that tests `hamiltonian_on_session`; `dcex.py:14-15` and
`:182` (the cache "cleared by `set_hamiltonian()`").

`restart=False` is undocumented, and the only meaning it can sensibly have is "keep
the previous state as a warm start". The readers do not agree on that: the public
`get_dynamical_correlator`, through `ground_state_on_session`, sees that the
session's Hamiltonian is not the chain's and re-solves; every other reader goes
through `gs_is_current`, whose key has no Hamiltonian, and answers with H1's
energy, state and band edges. `get_excited(n=2)` returns [-1.748511, -1.298528],
which are exactly the Ritz values of H2 in the span of H1's lowest four
eigenstates. `867e2b4`'s Status line for the second pass's finding 12, "a stored
state whose Hamiltonian is no longer the session's is solved again", holds for the
correlator path only. It survived because nothing in `tests/`, `examples/` or
`benchmarks/` passes `restart=False`; the only in-tree caller,
`mpsjulialive/dynamics.py`, calls `restart()` itself first.

**Expected**: ED of H2 (4-site Heisenberg plus 0.8*Sx on site 0, a term that does
not commute with H1, so no warm start can be trapped): E0, E1 = -1.780099,
-1.321921 and <Sx_0> = -0.349221, and every read on the chain agreeing with it.

Repro, from the `groundstate` hunter (`04` with a non-commuting H2, `01` with a
uniform field, which commutes):

```bash
cd <scratch>/groundstate && <scratch>/run3.sh 04_restart_false_order.py 2>&1 | tee 04_restart_false_order.out
```

`groundstate/04_restart_false_order.py`:

```python
# set_hamiltonian(H2, restart=False) on a solved chain, H2 NOT commuting with
# H1 (a transverse field on site 0 only), so no warm start can be trapped in an H1
# eigenstate.  What do gs_energy()/get_gs()/get_excited() answer, before and
# after one correlator?  Anchor: ED of H2.
import numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc,Bx):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + Bx*sc.Sx[0]
    return h
n, Bx = 4, 0.8
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,Bx))
print("ED H2: E0,E1 =",np.round(ed.get_excited(n=2,mode="ED"),6),
      " <Sx_0> =",round(ed.vev(ed.Sx[0],mode="ED").real,6))
ed1 = spinchain.Spin_Chain(["S=1/2"]*n); ed1.set_hamiltonian(ham(ed1,0.0))
print("ED H1: E0 =",round(ed1.gs_energy(mode="ED"),6))
es = np.linspace(-0.5,3,60)
for v in ("python",3,2):
    np.random.seed(4)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,0.0)); sc.gs_energy()
    sc.set_hamiltonian(ham(sc,Bx),restart=False)
    a = sc.gs_energy(); sa = sc.vev(sc.Sx[0]).real
    b = quiet(lambda: sc.get_excited(n=2))
    quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]),es=es,delta=0.2))
    c = sc.gs_energy(); sc_ = sc.vev(sc.Sx[0]).real
    print("v=%-6s restart=False: gs_energy()=%.6f <Sx_0>=%.6f get_excited(n=2)=%s | after one KPM call: gs_energy()=%.6f <Sx_0>=%.6f"
          %(v,a,sa,np.round(b,6),c,sc_))
# dcex's cached basis across restart=False (python): loud or silent?
np.random.seed(5)
sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
sc.maxm=20; sc.nsweeps=12
sc.set_hamiltonian(ham(sc,0.0))
quiet(lambda: sc.get_dynamical_correlator(submode="EX",nex=4,name=(sc.Sz[0],sc.Sz[0]),es=es,delta=0.2))
sc.set_hamiltonian(ham(sc,Bx),restart=False)
try:
    quiet(lambda: sc.get_dynamical_correlator(submode="EX",nex=4,name=(sc.Sz[0],sc.Sz[0]),es=es,delta=0.2))
    print("EX after restart=False: returned (cached basis reused silently)")
except ValueError as e:
    print("EX after restart=False: ValueError:",str(e)[:120])
```

`groundstate/04_restart_false_order.out`:


```
ED H2: E0,E1 = [-1.780099 -1.321921]  <Sx_0> = -0.349221
ED H1: E0 = -1.616025
v=python restart=False: gs_energy()=-1.616025 <Sx_0>=0.000000 get_excited(n=2)=[-1.748511 -1.298528] | after one KPM call: gs_energy()=-1.780099 <Sx_0>=-0.349221
v=3      restart=False: gs_energy()=-1.616025 <Sx_0>=-0.000000 get_excited(n=2)=[-1.748511 -1.298528] | after one KPM call: gs_energy()=-1.780099 <Sx_0>=-0.349221
v=2      restart=False: gs_energy()=-1.616025 <Sx_0>=-0.000000 get_excited(n=2)=[-1.748511 -1.298528] | after one KPM call: gs_energy()=-1.780099 <Sx_0>=-0.349221
EX after restart=False: ValueError: submode='EX' measures from the chain's ground state expressed in its cached excited-state basis, and that state lies out
```

```bash
cd <scratch>/groundstate && <scratch>/run3.sh 01_restart_false.py 2>&1 | tee 01_restart_false.out
```

`groundstate/01_restart_false.py`:

```python
# set_hamiltonian(H2, restart=False) after a solve: which reads answer for H2?
import numpy as np
from dmrgpy import spinchain
import dmrgpy; print(dmrgpy.__file__)

def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h

n = 4
ed = spinchain.Spin_Chain(["S=1/2"]*n)
ed.set_hamiltonian(ham(ed,1.0))
es_ed = ed.get_excited(n=3,mode="ED")
szt = sum(ed.Sz)
print("ED H2: E0..2 =",np.round(es_ed,6),"  <Sz_tot> =",np.round(ed.vev(szt,mode="ED").real,6))
ed1 = spinchain.Spin_Chain(["S=1/2"]*n); ed1.set_hamiltonian(ham(ed1,0.0))
print("ED H1: E0..2 =",np.round(ed1.get_excited(n=3,mode="ED"),6))

for v in ("python",3,2):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm = 20; sc.nsweeps = 10
    sc.set_hamiltonian(ham(sc,0.0))
    e1 = sc.gs_energy()
    sc.set_hamiltonian(ham(sc,1.0),restart=False)
    a = sc.gs_energy()
    b = sc.get_excited(n=3)
    c = sc.vev(sum(sc.Sz)).real
    es = np.linspace(-0.5,3,50)
    sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]),es=es,delta=0.2)
    d = sc.gs_energy()
    f = sc.get_excited(n=3)
    g = sc.vev(sum(sc.Sz)).real
    print("v=%-6s H1 E0=%.6f | after restart=False: gs_energy=%.6f get_excited=%s <Sz_tot>=%.6f"
          %(v,e1,a,np.round(b,6),c))
    print("          after one KPM correlator:      gs_energy=%.6f get_excited=%s <Sz_tot>=%.6f"
          %(d,np.round(f,6),g))
```

`groundstate/01_restart_false.out`:


```
<repo>/src/dmrgpy/__init__.py
ED H2: E0..2 = [-1.957107 -1.616025 -1.25    ]   <Sz_tot> = -1.0
ED H1: E0..2 = [-1.616025 -0.957107 -0.957107]
v=python H1 E0=-1.616025 | after restart=False: gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -1.069016] <Sz_tot>=-0.000000
          after one KPM correlator:      gs_energy=-1.957107 get_excited=[-1.957107 -1.616025 -1.25    ] <Sz_tot>=-1.000000
v=3      H1 E0=-1.616025 | after restart=False: gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -1.197311] <Sz_tot>=0.000000
          after one KPM correlator:      gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -1.25    ] <Sz_tot>=-0.000000
v=2      H1 E0=-1.616025 | after restart=False: gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -0.957107] <Sz_tot>=0.000000
          after one KPM correlator:      gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -1.25    ] <Sz_tot>=-0.000000
```


and from the reviewer of finding 11, who reached the same defect through the direct
KPM entry points (`07` on an 8-site chain with a staggered field, `08` the same on
a read-only archive of `30200a4`):

```bash
cd <scratch>/reviews/kpm_2 && <scratch>/run3.sh 07_direct_kpm_after_restart_false.py 2>&1 | tee 07_direct_kpm_after_restart_false.out
```

`reviews/kpm_2/07_direct_kpm_after_restart_false.py`:

```python
# Reviewer probe, a lead next to kpm_2: the direct KPM callers
# (get_dynamical_correlator_moments, which fermionchain.get_gr's
# kpmdmrg.get_dynamical_correlator shares) call only self.get_gs(), not
# ground_state_on_session, so after set_hamiltonian(H2, restart=False) on a
# solved chain do they answer with H1?  H2 = H1 + a staggered field, whose
# ground state differs from H1's. Reference: a fresh chain on H2.
import sys
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain

L = 8
def chain(v):
    sc = spinchain.Spin_Chain([2]*L, itensor_version=v)
    sc.maxm = 30; sc.nsweeps = 8; sc.kpmmaxm = 30
    return sc
def ham(sc, hs):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L):
        h = h + hs*(-1)**i*sc.Sz[i]
    return h

for v in ("python", 3):
    ref = chain(v); ref.set_hamiltonian(ham(ref, 1.0))
    mus_ref, emin_r, emax_r, *_ = ref.get_dynamical_correlator_moments(name=(ref.Sz[0], ref.Sz[0]), delta=0.2)
    e_ref = ref.e0
    # direct route after restart=False
    sc = chain(v); sc.set_hamiltonian(ham(sc, 0.0)); e1 = sc.gs_energy()
    sc.set_hamiltonian(ham(sc, 1.0), restart=False)
    mus_d, emin_d, emax_d, *_ = sc.get_dynamical_correlator_moments(name=(sc.Sz[0], sc.Sz[0]), delta=0.2)
    e_d = sc.e0
    # public route after restart=False
    sc2 = chain(v); sc2.set_hamiltonian(ham(sc2, 0.0)); sc2.gs_energy()
    sc2.set_hamiltonian(ham(sc2, 1.0), restart=False)
    x2, y2 = sc2.get_dynamical_correlator(name=(sc2.Sz[0], sc2.Sz[0]), delta=0.2, es=np.linspace(0, 3, 31))
    x3, y3 = ref.get_dynamical_correlator(name=(ref.Sz[0], ref.Sz[0]), delta=0.2, es=np.linspace(0, 3, 31))
    n = min(len(mus_ref), len(mus_d))
    print("v=%-6s E(H1)=%.6f E(H2) fresh=%.6f | direct after restart=False: e0=%.6f emax=%.4f (fresh %.4f) max|mu-mu_ref|=%.3e | public after restart=False: e0=%.6f max|y-y_ref|=%.3e (peak %.4f)" % (
        v, e1, e_ref, e_d, emax_d, emax_r, np.max(np.abs(np.array(mus_d[:n])-np.array(mus_ref[:n]))),
        sc2.e0, np.max(np.abs(y2-y3)), np.max(np.abs(y3))))
    sys.stdout.flush()
```

`reviews/kpm_2/07_direct_kpm_after_restart_false.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
v=python E(H1)=-3.374933 E(H2) fresh=-6.327539 | direct after restart=False: e0=-3.374933 emax=1.7500 (fresh 3.3486) max|mu-mu_ref|=3.598e-01 | public after restart=False: e0=-6.327539 max|y-y_ref|=1.679e-13 (peak 0.7446)
v=3      E(H1)=-3.374933 E(H2) fresh=-6.327539 | direct after restart=False: e0=-3.374933 emax=1.7500 (fresh 3.3486) max|mu-mu_ref|=3.598e-01 | public after restart=False: e0=-6.327539 max|y-y_ref|=1.004e-13 (peak 0.7446)
```

```bash
cd <scratch>/reviews/kpm_2 && <scratch>/run3.sh 08_direct_kpm_after_restart_false_30200a4.py 2>&1 | tee 08_direct_kpm_after_restart_false_30200a4.out
```

`reviews/kpm_2/08_direct_kpm_after_restart_false_30200a4.py`:

```python
# Reviewer probe: 07_direct_kpm_after_restart_false.py against ./old30200a4 (read-only git archive of 30200a4), itensor_version="python" only.
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "old30200a4", "src"))
import sys
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain

L = 8
def chain(v):
    sc = spinchain.Spin_Chain([2]*L, itensor_version=v)
    sc.maxm = 30; sc.nsweeps = 8; sc.kpmmaxm = 30
    return sc
def ham(sc, hs):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L):
        h = h + hs*(-1)**i*sc.Sz[i]
    return h

for v in ("python",):
    ref = chain(v); ref.set_hamiltonian(ham(ref, 1.0))
    mus_ref, emin_r, emax_r, *_ = ref.get_dynamical_correlator_moments(name=(ref.Sz[0], ref.Sz[0]), delta=0.2)
    e_ref = ref.e0
    # direct route after restart=False
    sc = chain(v); sc.set_hamiltonian(ham(sc, 0.0)); e1 = sc.gs_energy()
    sc.set_hamiltonian(ham(sc, 1.0), restart=False)
    mus_d, emin_d, emax_d, *_ = sc.get_dynamical_correlator_moments(name=(sc.Sz[0], sc.Sz[0]), delta=0.2)
    e_d = sc.e0
    # public route after restart=False
    sc2 = chain(v); sc2.set_hamiltonian(ham(sc2, 0.0)); sc2.gs_energy()
    sc2.set_hamiltonian(ham(sc2, 1.0), restart=False)
    x2, y2 = sc2.get_dynamical_correlator(name=(sc2.Sz[0], sc2.Sz[0]), delta=0.2, es=np.linspace(0, 3, 31))
    x3, y3 = ref.get_dynamical_correlator(name=(ref.Sz[0], ref.Sz[0]), delta=0.2, es=np.linspace(0, 3, 31))
    n = min(len(mus_ref), len(mus_d))
    print("v=%-6s E(H1)=%.6f E(H2) fresh=%.6f | direct after restart=False: e0=%.6f emax=%.4f (fresh %.4f) max|mu-mu_ref|=%.3e | public after restart=False: e0=%.6f max|y-y_ref|=%.3e (peak %.4f)" % (
        v, e1, e_ref, e_d, emax_d, emax_r, np.max(np.abs(np.array(mus_d[:n])-np.array(mus_ref[:n]))),
        sc2.e0, np.max(np.abs(y2-y3)), np.max(np.abs(y3))))
    sys.stdout.flush()
```

`reviews/kpm_2/08_direct_kpm_after_restart_false_30200a4.out`:


```
dmrgpy from <scratch>/reviews/kpm_2/old30200a4/src/dmrgpy/__init__.py
v=python E(H1)=-3.374933 E(H2) fresh=-6.327539 | direct after restart=False: e0=-3.374933 emax=1.7500 (fresh 3.3486) max|mu-mu_ref|=3.598e-01 | public after restart=False: e0=-6.327539 max|y-y_ref|=1.420e-13 (peak 0.7446)
```


**Reviewer of the `groundstate` candidate (CONFIRMED)**: the hunter's `04`
reproduces to every printed digit on all three backends. The anchor is a numpy
Kronecker build of H1 and H2 (lowest levels [-1.780099, -1.321921, -0.90355],
<Sx_0> = -0.349221), and the Ritz values of H2 in the span of H1's lowest four
eigenstates are exactly [-1.748511, -1.298528], confirming the `excited.py`
mechanism. Four sites at `maxm=20` is full bond dimension, and after the
correlator's re-solve every backend lands on exact H2 values, so nothing here is
an unconverged schedule. Measured rather than read: on the `restart=False` chain
`get_dynamical_correlator_moments` agrees with a fresh H1 chain to 3.0e-15
(`"python"`) and 5.1e-13 (v3), with emin -1.616025, emax 0.75 and n=31, against
emax 1.15 and n=38 for a fresh H2 chain, and it does not even trigger the
re-solve. On a read-only archive of `867e2b4^` the same sequence prints the
identical lines, so both the stale reads and the order dependence predate
`867e2b4`. `dcex` is loud with a non-commuting H2 (the span check raises) and
silent with a commuting one (H2's ground state lies inside H1's cached basis), where
it returns a result from H1's basis rediagonalized with H2; that differs from a
fresh H2 chain by 5.8e-1 on a 0.52 peak, but is two different `nex=4`
truncations, the stale one closer to ED, so it is "silently a different basis",
not "wrong"; what survives there is the false "cleared by `set_hamiltonian()`".
With a commuting field, after the re-solve `"python"` moves to -1.957107 while v3
and v2 stay at -1.616025, H1's singlet being an exact H2 eigenstate that traps the
warm start. Struck, each explicitly:

- The v3 value of the third `get_excited` level in the commuting run (-1.197311):
  run-dependent; a rerun gave -0.957107. The level is wrong on every backend
  against -1.25, but the specific value does not reproduce.
- Any tie of the stale reads to `867e2b4`: the parent gives identical output; what
  `867e2b4` contributed is a Status line and a docstring that overclaim.
- The suggested fix's implied promise that a warm re-solve makes every reader
  right: on v2 and v3 a commuting H2 traps the warm start at H1's singlet, so the
  fix buys consistency between readers, not correctness, whenever H1 and H2 share
  an eigenstate.
- "`ground_state_on_session`'s own extra test then becomes redundant": it also
  catches a rebuilt session with an unchanged Hamiltonian.

```bash
cd <scratch>/reviews/groundstate_3 && <scratch>/run3.sh 01_hunter_repro.py 2>&1 | tee 01_hunter_repro.out
```

`reviews/groundstate_3/01_hunter_repro.py` is `groundstate/04_restart_false_order.py` (finding 5) unchanged, rerun; its output, `reviews/groundstate_3/01_hunter_repro.out`:


```
ED H2: E0,E1 = [-1.780099 -1.321921]  <Sx_0> = -0.349221
ED H1: E0 = -1.616025
v=python restart=False: gs_energy()=-1.616025 <Sx_0>=0.000000 get_excited(n=2)=[-1.748511 -1.298528] | after one KPM call: gs_energy()=-1.780099 <Sx_0>=-0.349221
v=3      restart=False: gs_energy()=-1.616025 <Sx_0>=-0.000000 get_excited(n=2)=[-1.748511 -1.298528] | after one KPM call: gs_energy()=-1.780099 <Sx_0>=-0.349221
v=2      restart=False: gs_energy()=-1.616025 <Sx_0>=0.000000 get_excited(n=2)=[-1.748511 -1.298528] | after one KPM call: gs_energy()=-1.780099 <Sx_0>=-0.349221
EX after restart=False: ValueError: submode='EX' measures from the chain's ground state expressed in its cached excited-state basis, and that state lies out
```

```bash
cd <scratch>/reviews/groundstate_3 && <scratch>/run3.sh 02_anchor_moments_excited.py 2>&1 | tee 02_anchor_moments_excited.out
```

`reviews/groundstate_3/02_anchor_moments_excited.py`:

```python
# Reviewer probe for groundstate_3.
# (a) independent numpy anchor: exact H2 spectrum, and the Ritz values of H2 in
#     the span of H1's lowest 4 eigenstates (what get_excited(n=2) would return
#     if the session computed n+2=4 states of H1 and rediagonalized H2 in them).
# (b) get_dynamical_correlator_moments on the restart=False chain: H1 or H2 moments?
# (c) get_excited(n=2) after a correlator has re-solved: does it become right?
import numpy as np, io, contextlib
import dmrgpy
from dmrgpy import spinchain
print(dmrgpy.__file__)
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc,Bx):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h + Bx*sc.Sx[0]
n, Bx = 4, 0.8
# ---- (a) numpy anchor
sx = np.array([[0,1],[1,0]])/2; sy = np.array([[0,-1j],[1j,0]])/2; sz = np.diag([0.5,-0.5])
def op(o,i):
    m = np.array([[1.0]])
    for k in range(n): m = np.kron(m, o if k==i else np.eye(2))
    return m
H1 = sum(op(s,i)@op(s,i+1) for i in range(n-1) for s in (sx,sy,sz))
H2 = H1 + Bx*op(sx,0)
e1,v1 = np.linalg.eigh(H1); e2,v2 = np.linalg.eigh(H2)
print("numpy H1 lowest 5:",np.round(e1[:5],6))
print("numpy H2 lowest 3:",np.round(e2[:3],6), " <Sx_0> in H2 GS:",round(np.real(v2[:,0].conj()@op(sx,0)@v2[:,0]),6))
P = v1[:,:4]
print("Ritz values of H2 in span(H1 lowest 4):",np.round(np.linalg.eigvalsh(P.conj().T@H2@P),6))
# ---- (b),(c)
for v in ("python",3):
    np.random.seed(4)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v); sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,0.0)); sc.gs_energy()
    sc.set_hamiltonian(ham(sc,Bx),restart=False)
    m_stale = quiet(lambda: sc.get_dynamical_correlator_moments(name=(sc.Sz[0],sc.Sz[0]),delta=0.2))
    e_after_mom = sc.gs_energy()
    quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]),es=np.linspace(-0.5,3,40),delta=0.2))
    e_after_kpm = sc.gs_energy()
    ex_after = quiet(lambda: sc.get_excited(n=2))
    refs = {}
    for lab,B in (("H1",0.0),("H2",Bx)):
        np.random.seed(4)
        f = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v); f.maxm=20; f.nsweeps=12
        f.set_hamiltonian(ham(f,B)); f.gs_energy()
        refs[lab] = quiet(lambda: f.get_dynamical_correlator_moments(name=(f.Sz[0],f.Sz[0]),delta=0.2))
    mus = np.array(m_stale[0])
    for lab in ("H1","H2"):
        r = refs[lab]; k = min(len(mus),len(r[0]))
        print("v=%-6s moments(restart=False chain) vs fresh %s chain: max|dmu|=%.3e  emin %.6f vs %.6f  emax %.6f vs %.6f  n %d vs %d"
              %(v,lab,np.max(np.abs(mus[:k]-np.array(r[0])[:k])),m_stale[1],r[1],m_stale[2],r[2],m_stale[4],r[4]))
    print("v=%-6s gs_energy() after moments=%.6f, after one KPM correlator=%.6f, get_excited(n=2) after it=%s"
          %(v,e_after_mom,e_after_kpm,np.round(ex_after,6)))
```

`reviews/groundstate_3/02_anchor_moments_excited.out`:


```
<repo>/src/dmrgpy/__init__.py
numpy H1 lowest 5: [-1.616025 -0.957107 -0.957107 -0.957107 -0.25    ]
numpy H2 lowest 3: [-1.780099 -1.321921 -0.90355 ]  <Sx_0> in H2 GS: -0.349221
Ritz values of H2 in span(H1 lowest 4): [-1.748511 -1.298528 -0.824621 -0.615685]
v=python moments(restart=False chain) vs fresh H1 chain: max|dmu|=2.998e-15  emin -1.616025 vs -1.616025  emax 0.750000 vs 0.750000  n 31 vs 31
v=python moments(restart=False chain) vs fresh H2 chain: max|dmu|=3.022e-01  emin -1.616025 vs -1.780099  emax 0.750000 vs 1.150000  n 31 vs 38
v=python gs_energy() after moments=-1.616025, after one KPM correlator=-1.780099, get_excited(n=2) after it=[-1.780099 -1.321921]
v=3      moments(restart=False chain) vs fresh H1 chain: max|dmu|=5.064e-13  emin -1.616025 vs -1.616025  emax 0.750000 vs 0.750000  n 31 vs 31
v=3      moments(restart=False chain) vs fresh H2 chain: max|dmu|=3.022e-01  emin -1.616025 vs -1.780099  emax 0.750000 vs 1.150000  n 31 vs 38
v=3      gs_energy() after moments=-1.616025, after one KPM correlator=-1.780099, get_excited(n=2) after it=[-1.780099 -1.321921]
```

```bash
cd <scratch>/reviews/groundstate_3 && <scratch>/run3.sh 03_ex_commuting.py 2>&1 | tee 03_ex_commuting.out
```

`reviews/groundstate_3/03_ex_commuting.py`:

```python
# Reviewer probe for groundstate_3: submode="EX" after set_hamiltonian(H2,restart=False)
# with H2 = H1 + B*Sz_total, which COMMUTES with H1, so H2's ground state (the
# Sz=-1 member of H1's triplet) can lie inside H1's cached nex=4 basis and pass
# dcex's span check.  Loud (ValueError) or silent, and if silent, right or wrong?
# Anchor: EX on a fresh H2 chain at the same nex, and mode="ED" submode="ED".
import numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, B, nex = 4, 1.0, 4
es = np.linspace(-1.0,2.0,301); d = 0.05
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
x,yed = ed.get_dynamical_correlator(mode="ED",submode="ED",name=(ed.Sx[0],ed.Sx[0]),es=es,delta=d)
def peaks(y):
    y = np.real(y); idx = [k for k in range(1,len(y)-1) if y[k]>y[k-1] and y[k]>y[k+1] and y[k]>0.05*np.max(y)]
    return [(round(es[k],3),round(y[k],3)) for k in idx]
print("ED  H2 Sx0Sx0 peaks (omega,height):",peaks(yed))
for v in ("python",):
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v); sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,0.0)); sc.gs_energy()
    quiet(lambda: sc.get_dynamical_correlator(submode="EX",nex=nex,name=(sc.Sx[0],sc.Sx[0]),es=es,delta=d))
    sc.set_hamiltonian(ham(sc,B),restart=False)
    try:
        x,ys = quiet(lambda: sc.get_dynamical_correlator(submode="EX",nex=nex,name=(sc.Sx[0],sc.Sx[0]),es=es,delta=d))
        print("v=%s EX after restart=False: returned; gs_energy()=%.6f <Sz_tot>=%.6f"%(v,sc.gs_energy(),sc.vev(sum(sc.Sz)).real))
        print("   peaks:",peaks(ys))
    except ValueError as e:
        ys = None; print("v=%s EX after restart=False: ValueError: %s"%(v,str(e)[:100]))
    np.random.seed(7)
    f = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v); f.maxm=20; f.nsweeps=12
    f.set_hamiltonian(ham(f,B)); f.gs_energy()
    x,yf = quiet(lambda: f.get_dynamical_correlator(submode="EX",nex=nex,name=(f.Sx[0],f.Sx[0]),es=es,delta=d))
    print("v=%s EX fresh H2 chain: gs_energy()=%.6f peaks: %s"%(v,f.gs_energy(),peaks(yf)))
    if ys is not None:
        print("   max|EX(restart=False) - EX(fresh H2)| = %.3e on a peak of %.3e"%(np.max(np.abs(np.real(ys)-np.real(yf))),np.max(np.real(yf))))
```

`reviews/groundstate_3/03_ex_commuting.out`:


```
ED  H2 Sx0Sx0 peaks (omega,height): [(np.float64(0.34), np.float64(0.527)), (np.float64(0.71), np.float64(0.143)), (np.float64(1.0), np.float64(0.587)), (np.float64(1.71), np.float64(0.066))]
v=python EX after restart=False: returned; gs_energy()=-1.957107 <Sz_tot>=-1.000000
   peaks: [(np.float64(0.34), np.float64(0.525)), (np.float64(1.0), np.float64(0.583))]
v=python EX fresh H2 chain: gs_energy()=-1.957107 peaks: [(np.float64(0.34), np.float64(0.523)), (np.float64(0.71), np.float64(0.126))]
   max|EX(restart=False) - EX(fresh H2)| = 5.765e-01 on a peak of 5.234e-01
```

```bash
cd <scratch>/reviews/groundstate_3 && <scratch>/run3.sh 04_origin_pre867e2b4.py 2>&1 | tee 04_origin_pre867e2b4.out
```

`reviews/groundstate_3/04_origin_pre867e2b4.py`:

```python
# Reviewer probe for groundstate_3: the same sequence on the tree BEFORE 867e2b4
# (git archive of e5049b0 = 867e2b4^, python backend only, no C++ in the archive),
# to date the stale reads and the order dependence.
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),"tree_30200a4","src"))
import numpy as np, io, contextlib
import dmrgpy
from dmrgpy import spinchain
print(dmrgpy.__file__)
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc,Bx):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h + Bx*sc.Sx[0]
n, Bx = 4, 0.8
np.random.seed(4)
sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python"); sc.maxm=20; sc.nsweeps=12
sc.set_hamiltonian(ham(sc,0.0)); sc.gs_energy()
sc.set_hamiltonian(ham(sc,Bx),restart=False)
a = sc.gs_energy(); sa = sc.vev(sc.Sx[0]).real
b = quiet(lambda: sc.get_excited(n=2))
quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]),es=np.linspace(-0.5,3,40),delta=0.2))
c = sc.gs_energy(); s2 = sc.vev(sc.Sx[0]).real
print("pre-867e2b4 python: restart=False gs_energy()=%.6f <Sx_0>=%.6f get_excited(n=2)=%s | after one KPM call gs_energy()=%.6f <Sx_0>=%.6f"
      %(a,sa,np.round(b,6),c,s2))
```

`reviews/groundstate_3/04_origin_pre867e2b4.out`:


```
<scratch>/reviews/groundstate_3/tree_30200a4/src/dmrgpy/__init__.py
pre-867e2b4 python: restart=False gs_energy()=-1.616025 <Sx_0>=0.000000 get_excited(n=2)=[-1.748511 -1.298528] | after one KPM call gs_energy()=-1.780099 <Sx_0>=-0.349221
```

```bash
cd <scratch>/reviews/groundstate_3 && <scratch>/run3.sh 05_hunter_commuting.py 2>&1 | tee 05_hunter_commuting.out
```

`reviews/groundstate_3/05_hunter_commuting.py` is `groundstate/01_restart_false.py` (finding 5) unchanged, rerun; its output, `reviews/groundstate_3/05_hunter_commuting.out`:


```
<repo>/src/dmrgpy/__init__.py
ED H2: E0..2 = [-1.957107 -1.616025 -1.25    ]   <Sz_tot> = -1.0
ED H1: E0..2 = [-1.616025 -0.957107 -0.957107]
v=python H1 E0=-1.616025 | after restart=False: gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -1.069016] <Sz_tot>=-0.000000
          after one KPM correlator:      gs_energy=-1.957107 get_excited=[-1.957107 -1.616025 -1.25    ] <Sz_tot>=-1.000000
v=3      H1 E0=-1.616025 | after restart=False: gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -0.957107] <Sz_tot>=-0.000000
          after one KPM correlator:      gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -1.25    ] <Sz_tot>=-0.000000
v=2      H1 E0=-1.616025 | after restart=False: gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -0.957107] <Sz_tot>=-0.000000
          after one KPM correlator:      gs_energy=-1.616025 get_excited=[-1.957107 -1.616025 -1.25    ] <Sz_tot>=-0.000000
```


**Reviewer of the `kpm` candidate (CONFIRMED, NARROWED)**: `07` reproduces
exactly on both backends: after `restart=False` the direct route gives e0 =
-3.374933 (H1's) against a fresh H2 value of -6.327539, emax 1.7500 against
3.3486, and moments 3.598e-01 off, while the public route is within 3.5e-13
(`"python"`) and 9.2e-14 (v3) of the fresh curve on a 0.7446 peak. On one chain
after `restart=False`, `gs_is_current` is True while `hamiltonian_on_session` is
False, bare `gs_energy()` returns -3.374933 and `vev(Sz[0])` 0.000000 against a
fresh -0.466067; after one public correlator every consumer is right. The runs are
deterministic, and the fresh H2 chain is the anchor. The origin is older than
`30200a4`: `restart=` dates from `7c0a71b` (2020) and `gs_is_current` from
`1b87543`. Struck, each explicitly:

- The location "the direct KPM entry points reach only `get_gs()`, not
  `ground_state_on_session`" as the defect: that is a symptom; the cause is
  `gs_is_current` and `solver_key`'s false premise, which bare `gs_energy()` and
  `vev()` read just the same.
- "`867e2b4`'s helper exists to fix it and was not applied to these callers" as
  the remedy: adding `ground_state_on_session` to the two direct callers would
  leave `gs_energy()`, `vev()` and every other `get_gs()` consumer answering with
  H1.

```bash
cd <scratch>/reviews/kpm_2_new2 && <scratch>/run3.sh 01_direct_kpm_after_restart_false.py 2>&1 | tee 01_direct_kpm_after_restart_false.out
```

`reviews/kpm_2_new2/01_direct_kpm_after_restart_false.py` is `reviews/kpm_2/07_direct_kpm_after_restart_false.py` (finding 5) unchanged, rerun; its output, `reviews/kpm_2_new2/01_direct_kpm_after_restart_false.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
v=python E(H1)=-3.374933 E(H2) fresh=-6.327539 | direct after restart=False: e0=-3.374933 emax=1.7500 (fresh 3.3486) max|mu-mu_ref|=3.598e-01 | public after restart=False: e0=-6.327539 max|y-y_ref|=3.460e-13 (peak 0.7446)
v=3      E(H1)=-3.374933 E(H2) fresh=-6.327539 | direct after restart=False: e0=-3.374933 emax=1.7500 (fresh 3.3486) max|mu-mu_ref|=3.598e-01 | public after restart=False: e0=-6.327539 max|y-y_ref|=9.215e-14 (peak 0.7446)
```

```bash
cd <scratch>/reviews/kpm_2_new2 && <scratch>/run3.sh 02_every_consumer_after_restart_false.py 2>&1 | tee 02_every_consumer_after_restart_false.out
```

`reviews/kpm_2_new2/02_every_consumer_after_restart_false.py`:

```python
# Reviewer probe for kpm_2_new2: is the stale answer after
# set_hamiltonian(H2, restart=False) specific to the direct KPM route, or does
# every get_gs()/gs_energy() consumer see it?  Same chain/params as the
# hunter's 07. On ONE chain after the switch, in order: bare gs_energy(),
# vev(Sz[0]), vev(H2) (energy of the stored state under H2), the direct
# moments call, then the public correlator, then gs_energy() again.
# Reference: a fresh chain on H2.
import sys
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain
from dmrgpy import groundstate

L = 8
def chain(v):
    sc = spinchain.Spin_Chain([2]*L, itensor_version=v)
    sc.maxm = 30; sc.nsweeps = 8; sc.kpmmaxm = 30
    return sc
def ham(sc, hs):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L):
        h = h + hs*(-1)**i*sc.Sz[i]
    return h

for v in ("python", 3):
    ref = chain(v); ref.set_hamiltonian(ham(ref, 1.0))
    e_ref = ref.gs_energy(); sz_ref = ref.vev(ref.Sz[0]).real
    mus_ref, _, emax_r, *_ = ref.get_dynamical_correlator_moments(name=(ref.Sz[0], ref.Sz[0]), delta=0.2)
    sc = chain(v); sc.set_hamiltonian(ham(sc, 0.0)); e1 = sc.gs_energy()
    sz1 = sc.vev(sc.Sz[0]).real
    sc.set_hamiltonian(ham(sc, 1.0), restart=False)
    print("v=%s fresh H2: E=%.6f <Sz0>=%.6f emax=%.4f | H1 solved: E=%.6f <Sz0>=%.6f" % (v, e_ref, sz_ref, emax_r, e1, sz1))
    print("  gs_is_current after restart=False:", groundstate.gs_is_current(sc),
          " hamiltonian_on_session:", groundstate.hamiltonian_on_session(sc))
    e_a = sc.gs_energy()
    print("  bare gs_energy()            = %.6f" % e_a)
    print("  vev(Sz[0])                  = %.6f" % sc.vev(sc.Sz[0]).real)
    print("  vev(H2) on stored state     = %.6f" % sc.vev(ham(sc, 1.0)).real)
    mus_d, _, emax_d, *_ = sc.get_dynamical_correlator_moments(name=(sc.Sz[0], sc.Sz[0]), delta=0.2)
    n = min(len(mus_ref), len(mus_d))
    print("  direct moments: e0=%.6f emax=%.4f max|mu-mu_ref|=%.3e" % (
        sc.e0, emax_d, np.max(np.abs(np.array(mus_d[:n])-np.array(mus_ref[:n])))))
    x2, y2 = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), delta=0.2, es=np.linspace(0, 3, 31))
    x3, y3 = ref.get_dynamical_correlator(name=(ref.Sz[0], ref.Sz[0]), delta=0.2, es=np.linspace(0, 3, 31))
    print("  public correlator: e0=%.6f max|y-y_ref|=%.3e (peak %.4f)" % (sc.e0, np.max(np.abs(y2-y3)), np.max(np.abs(y3))))
    print("  gs_energy() after public    = %.6f   vev(Sz[0]) = %.6f" % (sc.gs_energy(), sc.vev(sc.Sz[0]).real))
    mus_d2, _, emax_d2, *_ = sc.get_dynamical_correlator_moments(name=(sc.Sz[0], sc.Sz[0]), delta=0.2)
    print("  direct moments after public: emax=%.4f max|mu-mu_ref|=%.3e" % (
        emax_d2, np.max(np.abs(np.array(mus_d2[:n])-np.array(mus_ref[:n])))))
    sys.stdout.flush()
```

`reviews/kpm_2_new2/02_every_consumer_after_restart_false.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
v=python fresh H2: E=-6.327539 <Sz0>=-0.466067 emax=3.3486 | H1 solved: E=-3.374933 <Sz0>=-0.000000
  gs_is_current after restart=False: True  hamiltonian_on_session: False
  bare gs_energy()            = -3.374933
  vev(Sz[0])                  = -0.000000
  vev(H2) on stored state     = -3.374933
  direct moments: e0=-3.374933 emax=1.7500 max|mu-mu_ref|=3.598e-01
  public correlator: e0=-6.327539 max|y-y_ref|=9.337e-14 (peak 0.7446)
  gs_energy() after public    = -6.327539   vev(Sz[0]) = -0.466067
  direct moments after public: emax=3.3486 max|mu-mu_ref|=1.307e-13
v=3 fresh H2: E=-6.327539 <Sz0>=-0.466067 emax=3.3486 | H1 solved: E=-3.374933 <Sz0>=-0.000000
  gs_is_current after restart=False: True  hamiltonian_on_session: False
  bare gs_energy()            = -3.374933
  vev(Sz[0])                  = -0.000000
  vev(H2) on stored state     = -3.374933
  direct moments: e0=-3.374933 emax=1.7500 max|mu-mu_ref|=3.598e-01
  public correlator: e0=-6.327539 max|y-y_ref|=6.439e-14 (peak 0.7446)
  gs_energy() after public    = -6.327539   vev(Sz[0]) = -0.466067
  direct moments after public: emax=3.3486 max|mu-mu_ref|=4.069e-14
```


**Suggested fix**: put it at the event rather than in `solver_key`, which both
reviewers prefer, since fingerprinting `hamiltonian.to_terms()` inside
`_gs_solver_key` would add an O(terms) cost to `gs_is_current`, the hottest path
in the API. In `set_hamiltonian`, whatever `restart` is, set `computed_gs=False`
and drop `_dcex_excited_cache` and `has_ED_obj`, all keyed on the Hamiltonian,
while keeping `wf0` and the session's state as the warm start; that makes
`solver_key`'s docstring premise true as written, corrects the `dcex` docstring and
its error text, and also covers the ED reads of the second pass's recorded
`restart=False` TD lead. The non-skip, non-injected path of `gs_energy_single`
never reads `self.wf0`: it re-sends the Hamiltonian and sweeps from the session's
own state, so what `restart=False` keeps is the v2/v3 session warm start, the path
measured right here to 1e-13. That path runs with `skip=self.skip_dmrg_gs=True` and
is correct only because the session's `set_hamiltonian` invalidates its energy
cache, which a test should pin. Do not remove `ground_state_on_session`'s
`hamiltonian_on_session` test. The alternative, making `gs_is_current` also
require `hamiltonian_on_session`, needs checking first: after `set_gs()` the skip
path may not record a Hamiltonian send, which would make every later
`gs_energy()` re-sweep and could break
`test_a_repeated_correlator_on_a_solved_chain_does_not_resweep`. The
`set_hamiltonian` docstring should then say what `restart=False` promises, a warm
start, which on v2/v3 stays trapped in an H1 eigenstate that is also an H2
eigenstate. The regression should pin `gs_energy()`, `vev()`, `get_excited()`
and `get_dynamical_correlator_moments()` against ED of a non-commuting H2 right
after `restart=False`, on `"python"` and v3. This finding lands with finding 11,
which moves `ground_state_on_session`. NUMBERS CHANGE only after
`restart=False`: `gs_energy()`, `vev()`, `get_excited()` and the direct KPM
moments move from the H1 answer to the H2 one (-1.616025 to -1.780099 here);
correlators do not move.

### 6. `gs_energy(maxde=...)` returns the energy of the first, unrefined solve (-3.279373 against a stored -3.374932 and an exact -3.374933, 2.8 per cent, on an 8-site chain at `maxm=3`; 0.17 per cent in the example's own setting) and the next dynamical correlator throws the refinement away, re-solving at the original `maxm`

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `groundstate`

**Where**: `src/dmrgpy/groundstate.py:377-404` (the `maxde` block: the recursive
`gs_energy_single` at `:396`, the key restored at `:403`, `return out` at `:404`);
`groundstate.py:301-305` (`ground_state_on_session` re-solves because the send
cache is keyed on the doubled `maxm`); `examples/groundstate/
GS_enforce_maximum_fluctuation/main.py` (prints and plots the returned, unrefined
energy).

The `maxde` block solves, measures the fluctuation, and if it is above the
tolerance doubles `maxm` and calls itself, then restores the parameters. Two
things go wrong on the way out. The value returned is `out`, the first solve's,
while `self.e0` and `self.wf0` hold the refined state. And the restore puts
`_gs_solver_key` back on the original `maxm` while the send cache stays keyed on
the doubled one, so the next correlator's `ground_state_on_session` sees a
Hamiltonian that is not on the session and re-solves at the original `maxm`,
discarding the refinement; static `vev()` between the two calls reads the refined
state, so the two routes disagree about which state the chain holds. It survived
because no file under `tests/` uses `maxde`, and the one example asserts nothing.

**Expected**: ED E0 = -3.37493260. `gs_energy(maxde=1e-4)` should return the
energy of the state it leaves on the chain, and a correlator afterwards should
measure that state.

Repro, from the hunter:

```bash
cd <scratch>/groundstate && <scratch>/run3.sh 05_maxde_return.py 2>&1 | tee 05_maxde_return.out
```

`groundstate/05_maxde_return.py`:

```python
# gs_energy(maxde=...): the fluctuation-driven refinement at doubled maxm.
# What does the call return, what does the chain hold afterwards, and what
# does the next correlator measure from?  Anchor: ED.
import numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
n = 8
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed))
print("ED E0 = %.8f"%ed.gs_energy(mode="ED"))
for v in ("python",3,2):
    np.random.seed(6)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=3; sc.nsweeps=8
    sc.set_hamiltonian(ham(sc))
    r = quiet(lambda: sc.gs_energy(maxde=1e-4))
    e_after = sc.gs_energy()
    h = sc.hamiltonian; wf = sc.get_gs()
    ewf = (sc.aMb(wf,h,wf)/sc.overlap(wf,wf)).real
    quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]),es=np.linspace(0,3,30),delta=0.3))
    e_corr = sc.gs_energy()
    print("v=%-6s maxm=3: gs_energy(maxde=1e-4) returned %.8f | then gs_energy()=%.8f, <wf0|H|wf0>=%.8f | after one KPM call gs_energy()=%.8f"
          %(v,r,e_after,ewf,e_corr))
```

`groundstate/05_maxde_return.out`:


```
ED E0 = -3.37493260
v=python maxm=3: gs_energy(maxde=1e-4) returned -3.27937301 | then gs_energy()=-3.37493206, <wf0|H|wf0>=-3.37493206 | after one KPM call gs_energy()=-3.27937301
v=3      maxm=3: gs_energy(maxde=1e-4) returned -3.27925299 | then gs_energy()=-3.37493206, <wf0|H|wf0>=-3.37493206 | after one KPM call gs_energy()=-3.00000000
v=2      maxm=3: gs_energy(maxde=1e-4) returned -3.27925299 | then gs_energy()=-3.37493206, <wf0|H|wf0>=-3.37493206 | after one KPM call gs_energy()=-3.00000000
```


**Reviewer (CONFIRMED)**: reproduced with a control, a fresh chain with the same
seed and parameters and no `maxde`: on `"python"` the plain `gs_energy()` at
`maxm=3` gives -3.27937301 and `gs_energy(maxde=1e-4)` returns -3.27937301 with
two retries and `self.e0` = -3.37493206; on v3 and v2 the pair is -3.27925299 and
-3.27925299. So the returned number is exactly the first solve's, and the
refinement is real (ED <Sz0Sz1> = -0.22041272, the refined state -0.22041326).
After one KPM call the chain holds a different state: on `"python"` overlap
0.856164 with the refined one, E = -3.27937301, <Sz0Sz1> = -0.22303364; on v3 and
v2 overlap 0.764908, E = -3.00000000, <Sz0Sz1> = -0.25000000, a warm start from
the truncated `maxm=12` state trapped at a dimer product state, recorded as an
observation rather than a separate defect. At the example's own setting (10-site
S=1, `maxm=10`, v3), `maxde=1e-1` needs no retry, and `maxde=1e-3` does one and
returns -12.87252493 against a stored -12.89448122 and ED -12.89456013, so the
defect survives at a realistic `maxm` with a size of 0.17 per cent. On a read-only
archive of `e5049b0`, `"python"` gives identical numbers to every digit: `867e2b4`
moved the trigger (`ground_state_on_session` instead of `set_initial_wf(self.wf0)`)
and not the behaviour; `return out` is in `4731b5a` (2020). Nothing in the
repository's documentation covers the return value; `user_guide.md:699` presents
`get_gs(maxde=...)` as a convergence criterion. Struck, each explicitly:

- No substantive sub-claim. The origin is narrowed to "pre-existing, trigger moved
  in `867e2b4`".
- "The example prints the returned energy next to a fluctuation computed from the
  refined state": that fluctuation is evaluated after `maxm` is restored, so it is
  not a clean property of the refined state either (finding 7); both printed
  numbers are off.

```bash
cd <scratch>/reviews/groundstate_4 && <scratch>/run3.sh 01_maxde_return_review.py 2>&1 | tee 01_maxde_return_review.out
```

`reviews/groundstate_4/01_maxde_return_review.py`:

```python
# Review of groundstate_4: what gs_energy(maxde=...) returns, what the chain
# holds, and what a dynamical correlator afterwards measures from.
# Control: the same seed and parameters WITHOUT maxde, to identify the
# returned number as the first, unrefined solve.  Anchor: ED.
import numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
n = 8
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed))
e_ed = ed.gs_energy(mode="ED")
zz_ed = ed.vev(ed.Sz[0]*ed.Sz[1],mode="ED").real
print("ED E0 = %.8f   ED <Sz0Sz1> = %.8f"%(e_ed,zz_ed))
def make(v):
    np.random.seed(6)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=3; sc.nsweeps=8
    sc.set_hamiltonian(ham(sc))
    return sc
for v in ("python",3,2):
    ctrl = make(v); e_ctrl = ctrl.gs_energy()
    sc = make(v)
    log = io.StringIO()
    with contextlib.redirect_stdout(log): r = sc.gs_energy(maxde=1e-4)
    nretry = log.getvalue().count("Energy fluctuation")
    print("v=%s: plain gs_energy() at maxm=3 = %.8f ; gs_energy(maxde=1e-4) returned %.8f ; retries=%d ; self.e0=%.8f ; maxm now %d"
          %(v,e_ctrl,r,nretry,sc.e0,sc.maxm))
    wf_ref = sc.get_gs().copy()
    ewf = (sc.aMb(wf_ref,sc.hamiltonian,wf_ref)/sc.overlap(wf_ref,wf_ref)).real
    zz = sc.vev(sc.Sz[0]*sc.Sz[1]).real
    de = sc.gs_energy_fluctuation()
    print("   after: gs_energy()=%.8f <wf|H|wf>=%.8f vev(Sz0Sz1)=%.8f fluct=%.3e"%(sc.gs_energy(),ewf,zz,de))
    es = np.linspace(0,3,30)
    quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]),es=es,delta=0.3))
    wf_new = sc.get_gs()
    ov = abs(sc.overlap(wf_ref,wf_new))**2/abs(sc.overlap(wf_ref,wf_ref)*sc.overlap(wf_new,wf_new))
    enew = (sc.aMb(wf_new,sc.hamiltonian,wf_new)/sc.overlap(wf_new,wf_new)).real
    print("   after one KPM call: gs_energy()=%.8f <wf|H|wf>=%.8f |<refined|now>|^2=%.6f vev(Sz0Sz1)=%.8f"
          %(sc.gs_energy(),enew,ov,sc.vev(sc.Sz[0]*sc.Sz[1]).real))
```

`reviews/groundstate_4/01_maxde_return_review.out`:


```
ED E0 = -3.37493260   ED <Sz0Sz1> = -0.22041272
v=python: plain gs_energy() at maxm=3 = -3.27937301 ; gs_energy(maxde=1e-4) returned -3.27937301 ; retries=2 ; self.e0=-3.37493206 ; maxm now 3
   after: gs_energy()=-3.37493206 <wf|H|wf>=-3.37493206 vev(Sz0Sz1)=-0.22041326 fluct=9.790e-01
   after one KPM call: gs_energy()=-3.27937301 <wf|H|wf>=-3.27937301 |<refined|now>|^2=0.856164 vev(Sz0Sz1)=-0.22303364
v=3: plain gs_energy() at maxm=3 = -3.27925299 ; gs_energy(maxde=1e-4) returned -3.27925299 ; retries=2 ; self.e0=-3.37493206 ; maxm now 3
   after: gs_energy()=-3.37493206 <wf|H|wf>=-3.37493206 vev(Sz0Sz1)=-0.22041323 fluct=1.636e+00
   after one KPM call: gs_energy()=-3.00000000 <wf|H|wf>=-3.00000000 |<refined|now>|^2=0.764908 vev(Sz0Sz1)=-0.25000000
v=2: plain gs_energy() at maxm=3 = -3.27925299 ; gs_energy(maxde=1e-4) returned -3.27925299 ; retries=2 ; self.e0=-3.37493206 ; maxm now 3
   after: gs_energy()=-3.37493206 <wf|H|wf>=-3.37493206 vev(Sz0Sz1)=-0.22041327 fluct=9.790e-01
   after one KPM call: gs_energy()=-3.00000000 <wf|H|wf>=-3.00000000 |<refined|now>|^2=0.764908 vev(Sz0Sz1)=-0.25000000
```

```bash
cd <scratch>/reviews/groundstate_4 && <scratch>/run3.sh 02_fluct_truncation_and_example.py 2>&1 | tee 02_fluct_truncation_and_example.out
```

`reviews/groundstate_4/02_fluct_truncation_and_example.py`:

```python
# (a) gs_energy_fluctuation / vev(npow=2) on a FIXED, essentially exact state
#     as the chain's maxm changes: is it a property of the state, or of the
#     truncation of H|psi> at maxm?  (8-site S=1/2 Heisenberg, state solved at
#     maxm=16 = full bond dimension; variance should be ~0.)
# (b) the example's own setting, 10-site S=1 at maxm=10: returned energy vs
#     the stored refined one.
import numpy as np, io, contextlib
from dmrgpy import spinchain
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
n = 8
for v in ("python",3):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=16; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc))
    e = sc.gs_energy(); wf = sc.get_gs().copy()
    h = sc.hamiltonian
    row = []
    for m in (3,6,12,16,40,100):
        sc.maxm = m
        e1 = sc.vev(h,wf=wf).real; e2 = sc.vev(h,wf=wf,npow=2).real
        row.append("maxm=%d: %.3e"%(m,np.sqrt(abs(e2-e1**2))))
    print("(a) v=%s E=%.10f  sqrt(|<H^2>-<H>^2|) of the SAME state at "%(v,e)+" | ".join(row))
spins = ["S=1"]*10
ed = spinchain.Spin_Chain(spins); ed.set_hamiltonian(ham(ed))
print("(b) ED E0 (10-site S=1) = %.8f"%ed.gs_energy(mode="ED"))
for v in (3,):
    for maxde in (1e-1,1e-3):
        np.random.seed(2)
        sc = spinchain.Spin_Chain(spins,itensor_version=v); sc.maxm=10
        sc.set_hamiltonian(ham(sc))
        log = io.StringIO()
        with contextlib.redirect_stdout(log): r = sc.gs_energy(maxde=maxde)
        print("(b) v=%s maxde=%g: returned %.8f, stored e0 %.8f, retries %d, then gs_energy_fluctuation()=%.3e"
              %(v,maxde,r,sc.e0,log.getvalue().count("Energy fluctuation"),sc.gs_energy_fluctuation()))
```

`reviews/groundstate_4/02_fluct_truncation_and_example.out`:


```
(a) v=python E=-3.3749325987  sqrt(|<H^2>-<H>^2|) of the SAME state at maxm=3: 9.790e-01 | maxm=6: 8.325e-02 | maxm=12: 1.171e-03 | maxm=16: 1.460e-07 | maxm=40: 1.460e-07 | maxm=100: 1.460e-07
(a) v=3 E=-3.3749325987  sqrt(|<H^2>-<H>^2|) of the SAME state at maxm=3: 1.636e+00 | maxm=6: 1.180e-01 | maxm=12: 1.310e-03 | maxm=16: 1.520e-07 | maxm=40: 1.520e-07 | maxm=100: 1.520e-07
(b) ED E0 (10-site S=1) = -12.89456013
(b) v=3 maxde=0.1: returned -12.87252493, stored e0 -12.87252493, retries 0, then gs_energy_fluctuation()=1.896e-02
(b) v=3 maxde=0.001: returned -12.87252493, stored e0 -12.89448122, retries 1, then gs_energy_fluctuation()=1.113e+00
```

```bash
cd <scratch>/reviews/groundstate_4 && <scratch>/run3.sh 03_maxde_parent.py 2>&1 | tee 03_maxde_parent.out
```

`reviews/groundstate_4/03_maxde_parent.py`:

```python
# 01 on the PARENT of 867e2b4 (e5049b0, hunter's read-only git archive), python only
import sys; sys.path.insert(0,'<scratch>/groundstate/parent_src/src')
import dmrgpy; print(dmrgpy.__file__)
# Review of groundstate_4: what gs_energy(maxde=...) returns, what the chain
# holds, and what a dynamical correlator afterwards measures from.
# Control: the same seed and parameters WITHOUT maxde, to identify the
# returned number as the first, unrefined solve.  Anchor: ED.
import numpy as np, io, contextlib
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
n = 8
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed))
e_ed = ed.gs_energy(mode="ED")
zz_ed = ed.vev(ed.Sz[0]*ed.Sz[1],mode="ED").real
print("ED E0 = %.8f   ED <Sz0Sz1> = %.8f"%(e_ed,zz_ed))
def make(v):
    np.random.seed(6)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=3; sc.nsweeps=8
    sc.set_hamiltonian(ham(sc))
    return sc
for v in ("python",):
    ctrl = make(v); e_ctrl = ctrl.gs_energy()
    sc = make(v)
    log = io.StringIO()
    with contextlib.redirect_stdout(log): r = sc.gs_energy(maxde=1e-4)
    nretry = log.getvalue().count("Energy fluctuation")
    print("v=%s: plain gs_energy() at maxm=3 = %.8f ; gs_energy(maxde=1e-4) returned %.8f ; retries=%d ; self.e0=%.8f ; maxm now %d"
          %(v,e_ctrl,r,nretry,sc.e0,sc.maxm))
    wf_ref = sc.get_gs().copy()
    ewf = (sc.aMb(wf_ref,sc.hamiltonian,wf_ref)/sc.overlap(wf_ref,wf_ref)).real
    zz = sc.vev(sc.Sz[0]*sc.Sz[1]).real
    de = sc.gs_energy_fluctuation()
    print("   after: gs_energy()=%.8f <wf|H|wf>=%.8f vev(Sz0Sz1)=%.8f fluct=%.3e"%(sc.gs_energy(),ewf,zz,de))
    es = np.linspace(0,3,30)
    quiet(lambda: sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]),es=es,delta=0.3))
    wf_new = sc.get_gs()
    ov = abs(sc.overlap(wf_ref,wf_new))**2/abs(sc.overlap(wf_ref,wf_ref)*sc.overlap(wf_new,wf_new))
    enew = (sc.aMb(wf_new,sc.hamiltonian,wf_new)/sc.overlap(wf_new,wf_new)).real
    print("   after one KPM call: gs_energy()=%.8f <wf|H|wf>=%.8f |<refined|now>|^2=%.6f vev(Sz0Sz1)=%.8f"
          %(sc.gs_energy(),enew,ov,sc.vev(sc.Sz[0]*sc.Sz[1]).real))
```

`reviews/groundstate_4/03_maxde_parent.out`:


```
<scratch>/groundstate/parent_src/src/dmrgpy/__init__.py
ED E0 = -3.37493260   ED <Sz0Sz1> = -0.22041272
v=python: plain gs_energy() at maxm=3 = -3.27937301 ; gs_energy(maxde=1e-4) returned -3.27937301 ; retries=2 ; self.e0=-3.37493206 ; maxm now 3
   after: gs_energy()=-3.37493206 <wf|H|wf>=-3.37493206 vev(Sz0Sz1)=-0.22041326 fluct=9.790e-01
   after one KPM call: gs_energy()=-3.27937301 <wf|H|wf>=-3.27937301 |<refined|now>|^2=0.856164 vev(Sz0Sz1)=-0.22303364
```


**Suggested fix**: `return self.e0` at the end of the `maxde` block. The hunter's
alternative, returning the recursive call's result, fixes nothing at depth 2,
where the inner call returns its own first solve. Both of the hunter's options for
making the refinement stick are wrong in shape: recording the refined parameters
as the solver key makes the very next `gs_energy()` re-solve at the restored
`maxm`, and re-sending the Hamiltonian under the original key drops the session's
energy and band-edge caches, so the next KPM sweeps the refined state at
`maxm=3` anyway. What fits is to re-key the send cache without re-sending once the
parameters are restored, `self._session_ham_cache = (self._session,
_send_key(self)[0])`: the terms on the session are identical, the sweep parameters
are pushed on every call anyway, and the only MPO field in the key,
`max(maxm, mpomaxm)`, is the same before and after the retry at the default
`mpomaxm=5000`. The regression should pin the returned value against `self.e0`,
the overlap with the refined state after one KPM call, and the correlator against
ED. Lands with finding 7, which rewrites the same block. NUMBERS CHANGE for every
`maxde=` caller: the returned value moves to the refined energy (-3.2794 to
-3.3749 on the 8-site chain at `maxm=3`), correlators after such a call measure the
refined state, and the example's printed energies change.

### 7. `gs_energy_fluctuation()` and `vev(H, npow=2)` truncate H|psi> at the chain's current `maxm` on every DMRG backend, so the ordinary solve-then-measure workflow under-reports the variance 10 to 51 times (6.7e-3 against a true 0.344), a state wider than `maxm` over-reports it at order one (1.64 against about 1e-7), and `gs_energy(maxde=...)` stops early and returns states 2.3 times over the tolerance it was asked for

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `groundstate` (turned up by the reviewer of finding 6)

**Where**: `src/dmrgpy/vev.py:31-34` (`multi_vev` pushes `self.maxm` to the
session before `vev(npow)`); `pyitensor/chain.py:707-709` (`_apply_mpo` with
`maxdim=self.maxm`, `68c96eb`); `mpscpp3/chain_session.h:1110-1112`
(`apply_mpo` with `MaxDim maxm_`, `8fc78c6`); `mpscpp2/chain_session.h:416`
(`exactApplyMPO` with `Maxm maxm_`, `a8233ab`); consumers
`manybodychain.py:1155` (`gs_energy_fluctuation`) and `groundstate.py:378-379`
(the `maxde` loop); `docs/user_guide.md:705-714`, the fluctuation paragraph.

The route computes <psi|H|trunc(H|psi>)> - <H>^2 with H|psi> truncated to the
chain's `maxm` at the moment of measurement, and that number is not the variance
in either direction. A state converged at bond dimension chi has its residual
(H-E)|psi> mostly outside the bond-chi manifold, so truncating H|psi> back to chi
throws away the very quantity being measured and the diagnostic under-reports;
a state whose own bond dimension exceeds `maxm` gets H|psi> cut below its own
rank and the diagnostic over-reports at order one. The case the library produces
by itself is `gs_energy(maxde=...)`, which refines at 2x, 4x, ... `maxm` and then
restores `maxm`. The design is pre-pybind (the `mpscpp2` comment says it mirrors
the file-based `vev.h`) and was carried into every current backend. The
2026-09 audit's finding 13 touched this route on `"python"` only (a zip-up
accuracy loss) and suggested that `gs_energy_fluctuation` not go through
`vev(npow=2)` at all, which was not adopted. It survived because the
user guide's measurement of the fluctuation floor was made near full bond
dimension (10-site S=1/2 at `maxm=30`, full being 32), the one regime where the
route is right.

**Expected**: the variance <H^2> - <H>^2 of the state, independent of `maxm`.
`vev(h*h)`, which builds H^2 as one MPO and applies no MPS truncation, gives it,
and agrees with the npow route evaluated at `maxm=300` (no truncation at all) to
four digits on both backends.

Repro, from the reviewer of finding 6, who found it (`02`, part a, evaluates the
fluctuation of one fixed, essentially exact 8-site state at six `maxm`; part b is
the example's own setting; both scripts are spliced under finding 6):

`reviews/groundstate_4/02_fluct_truncation_and_example.py` and its output are spliced under finding 6 above.

`reviews/groundstate_4/01_maxde_return_review.py` and its output are spliced under finding 6 above.


**Reviewer (CONFIRMED, NARROWED)**: the hunter's script reruns to every printed
digit. The truth is `vev(h*h)`, which agrees with the untruncated npow route to
four digits (0.2625 and 0.2625 on v3, 0.2653 and 0.2653 on `"python"`, 0.02283 and
0.02283 on the v3 refined state) and gives 1.26e-07 (`"python"`), 4.2e-08 (v3)
and 6.0e-08 (v2) on the exact 8-site state, where the npow route at `maxm=3`
reports 9.790e-01, 1.636e+00 and 9.805e-01; the v2 leg is now measured. The
regime the hunter did not test, a state solved and measured at the same `maxm`,
under-reports in every case:

| backend, `maxm` | reported | true |
|---|---|---|
| `"python"`, 3 | 6.725e-03 | 3.436e-01 (51x) |
| `"python"`, 6 | 3.996e-03 | 6.726e-02 (17x) |
| v3, 3 | 1.492e-02 | 3.474e-01 (23x) |
| v3, 6 | 6.795e-03 | 6.766e-02 (10x) |

and on the 10-site S=1 chain at `maxm=10`, 1.896e-02 against 2.625e-01 (v3) and
2.727e-02 against 2.653e-01 (`"python"`). At `maxde=1e-3` the loop stopped at
`maxm=20` reading 2.604e-04 per site (v3) and 1.672e-04 (`"python"`), where the
true values are 2.283e-03 and 2.290e-03, so the returned state is 2.3 times over
its tolerance. At `maxde=1e-6` the loop stops at `maxm=80` on a state exact to
every ED digit, true fluctuation 2.52e-05, while `gs_energy_fluctuation()` at the
restored `maxm=10` reports 1.694e+00. So all three points plotted by
`examples/groundstate/GS_enforce_maximum_fluctuation` are wrong: 1.9e-2, 1.113 and
1.694 against true values of 0.263, 0.0228 and 2.5e-5. Nothing in `CLAUDE.md`,
the known-issue files, `ROADMAP.md` or the four records covers truncation below
the rank of H|psi>. Struck, each explicitly:

- "The two backends also disagree with each other at small `maxm`" as an
  independent defect: truncating below a state's rank has no canonical answer,
  so two correct `applyMPO` implementations differ legitimately; kept as a
  symptom.
- The hunter's sqrt(bandwidth*dE) bound as the anchor: superseded by the measured
  `vev(h*h)` truths.
- The framing that the defect needs a state wider than `maxm`: at the solve's own
  `maxm` the route under-reports 10 to 51 times, and that is where the `maxde`
  criterion is evaluated.
- "`user_guide.md:705-714` attributes the fluctuation floor to Lanczos accuracy
  only", narrowed: the text is right for the regime it was measured in; what is
  wrong is its silence below full bond dimension.
- "The `maxde` criterion is measured the same way at each retry's `maxm`" is
  kept, now measured rather than read, with the consequence reversed: the loop
  under-enforces.

```bash
cd <scratch>/reviews/groundstate_4_new1 && <scratch>/run3.sh 01_hunter_repro.py 2>&1 | tee 01_hunter_repro.out
```

`reviews/groundstate_4_new1/01_hunter_repro.py` is `reviews/groundstate_4/02_fluct_truncation_and_example.py` (finding 6) unchanged, rerun; its output, `reviews/groundstate_4_new1/01_hunter_repro.out`:


```
(a) v=python E=-3.3749325987  sqrt(|<H^2>-<H>^2|) of the SAME state at maxm=3: 9.790e-01 | maxm=6: 8.325e-02 | maxm=12: 1.171e-03 | maxm=16: 1.460e-07 | maxm=40: 1.460e-07 | maxm=100: 1.460e-07
(a) v=3 E=-3.3749325987  sqrt(|<H^2>-<H>^2|) of the SAME state at maxm=3: 1.636e+00 | maxm=6: 1.180e-01 | maxm=12: 1.310e-03 | maxm=16: 1.460e-07 | maxm=40: 1.460e-07 | maxm=100: 1.460e-07
(b) ED E0 (10-site S=1) = -12.89456013
(b) v=3 maxde=0.1: returned -12.87252493, stored e0 -12.87252493, retries 0, then gs_energy_fluctuation()=1.896e-02
(b) v=3 maxde=0.001: returned -12.87252493, stored e0 -12.89448122, retries 1, then gs_energy_fluctuation()=1.113e+00
```

```bash
cd <scratch>/reviews/groundstate_4_new1 && <scratch>/run3.sh 02_anchor_and_regimes.py 2>&1 | tee 02_anchor_and_regimes.out
```

`reviews/groundstate_4_new1/02_anchor_and_regimes.py`:

```python
# Independent anchor for the fluctuation: vev(h*h) builds H^2 as ONE MPO
# (bond <= 25, mpomaxm=5000) and takes <psi|H^2|psi> with npow=1, so no
# MPS truncation happens anywhere.  Compare it with the npow=2 route
# (H|psi> truncated at self.maxm) in three regimes:
#  1. state bond 16 (8-site, full) measured at maxm=3/6 (chi_psi > maxm)
#  2. state solved at maxm=10 and measured at maxm=10 (chi_psi <= maxm)
#  3. maxde-refined state (chi_psi = 20) measured at restored maxm=10
import numpy as np, io, contextlib
from dmrgpy import spinchain
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
def fl_npow(sc,wf,m):
    old = sc.maxm; sc.maxm = m
    e = sc.vev(sc.hamiltonian,wf=wf).real
    e2 = sc.vev(sc.hamiltonian,wf=wf,npow=2).real
    sc.maxm = old
    return np.sqrt(abs(e2-e**2))
def fl_hh(sc,wf,hh):
    e = sc.vev(sc.hamiltonian,wf=wf).real
    e2 = sc.vev(hh,wf=wf).real
    return np.sqrt(abs(e2-e**2)), e
def maxbond(sc,wf):
    try:
        return max(wf.cpp_handle.bond_dimensions())
    except Exception:
        return "?"
print("== 1. 8-site S=1/2, state solved at maxm=16 (exact), measured at small maxm")
for v in ("python",3,2):
    np.random.seed(1)
    sc = spinchain.Spin_Chain(["S=1/2"]*8,itensor_version=v)
    sc.maxm=16; sc.nsweeps=12
    h = ham(sc); sc.set_hamiltonian(h); hh = h*h
    sc.gs_energy(); wf = sc.get_gs().copy()
    t,e = fl_hh(sc,wf,hh)
    print("v=%-6s E=%.10f  truth sqrt(|<H*H>-E^2|)=%.3e | npow route: maxm=3 %.3e  maxm=6 %.3e  maxm=16 %.3e"
          %(v,e,t,fl_npow(sc,wf,3),fl_npow(sc,wf,6),fl_npow(sc,wf,16)))
spins = ["S=1"]*10
print("== 2./3. 10-site S=1 Heisenberg, maxm=10")
for v in (3,"python"):
    for maxde in (None,1e-1,1e-3):
        np.random.seed(2)
        sc = spinchain.Spin_Chain(spins,itensor_version=v); sc.maxm=10
        h = ham(sc); sc.set_hamiltonian(h); hh = h*h
        log = io.StringIO()
        with contextlib.redirect_stdout(log): sc.gs_energy(maxde=maxde)
        inner = [l for l in log.getvalue().splitlines() if "Energy fluctuation" in l]
        wf = sc.get_gs()
        t,e = fl_hh(sc,wf,hh)
        print("v=%-6s maxde=%-5s E=%.8f  gs_energy_fluctuation()=%.3e  npow@maxm=10 %.3e  npow@maxm=300 %.3e  truth(H*H) %.3e  loop-log %s"
              %(v,maxde,e,sc.gs_energy_fluctuation(),fl_npow(sc,wf,10),fl_npow(sc,wf,300),t,inner))
```

`reviews/groundstate_4_new1/02_anchor_and_regimes.out`:


```
== 1. 8-site S=1/2, state solved at maxm=16 (exact), measured at small maxm
v=python E=-3.3749325987  truth sqrt(|<H*H>-E^2|)=1.264e-07 | npow route: maxm=3 9.790e-01  maxm=6 8.325e-02  maxm=16 1.460e-07
v=3      E=-3.3749325987  truth sqrt(|<H*H>-E^2|)=4.215e-08 | npow route: maxm=3 1.636e+00  maxm=6 1.180e-01  maxm=16 1.520e-07
v=2      E=-3.3749325987  truth sqrt(|<H*H>-E^2|)=5.960e-08 | npow route: maxm=3 9.805e-01  maxm=6 8.326e-02  maxm=16 8.429e-08
== 2./3. 10-site S=1 Heisenberg, maxm=10
v=3      maxde=None  E=-12.87252493  gs_energy_fluctuation()=1.896e-02  npow@maxm=10 1.896e-02  npow@maxm=300 2.625e-01  truth(H*H) 2.625e-01  loop-log []
v=3      maxde=0.1   E=-12.87252493  gs_energy_fluctuation()=1.896e-02  npow@maxm=10 1.896e-02  npow@maxm=300 2.625e-01  truth(H*H) 2.625e-01  loop-log []
v=3      maxde=0.001 E=-12.89448122  gs_energy_fluctuation()=1.113e+00  npow@maxm=10 1.113e+00  npow@maxm=300 2.283e-02  truth(H*H) 2.283e-02  loop-log ['Energy fluctuation =  0.0018958122517855418 10']
v=python maxde=None  E=-12.87233087  gs_energy_fluctuation()=2.727e-02  npow@maxm=10 2.727e-02  npow@maxm=300 2.653e-01  truth(H*H) 2.653e-01  loop-log []
v=python maxde=0.1   E=-12.87233087  gs_energy_fluctuation()=2.727e-02  npow@maxm=10 2.727e-02  npow@maxm=300 2.653e-01  truth(H*H) 2.653e-01  loop-log []
v=python maxde=0.001 E=-12.89448263  gs_energy_fluctuation()=1.116e+00  npow@maxm=10 1.116e+00  npow@maxm=300 2.290e-02  truth(H*H) 2.290e-02  loop-log ['Energy fluctuation =  0.002727295809732486 10']
```

```bash
cd <scratch>/reviews/groundstate_4_new1 && <scratch>/run3.sh 03_consistent_regime_and_maxde_check.py 2>&1 | tee 03_consistent_regime_and_maxde_check.out
```

`reviews/groundstate_4_new1/03_consistent_regime_and_maxde_check.py`:

```python
# (A) consistent regime on the 8-site chain: state SOLVED at maxm=m and
#     measured at the same maxm=m, npow route vs vev(h*h) truth.
# (B) what the maxde loop itself read: the refined 10-site S=1 state
#     measured by the npow route at the retry's own maxm=20, per site,
#     against the truth per site and the requested maxde=1e-3.
# (C) the example's third point, maxde=1e-6, on v3 (per-site loop values).
import numpy as np, io, contextlib, time
from dmrgpy import spinchain
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
def fl_npow(sc,wf,m):
    old = sc.maxm; sc.maxm = m
    e = sc.vev(sc.hamiltonian,wf=wf).real
    e2 = sc.vev(sc.hamiltonian,wf=wf,npow=2).real
    sc.maxm = old
    return np.sqrt(abs(e2-e**2))
def truth(sc,wf,hh):
    e = sc.vev(sc.hamiltonian,wf=wf).real
    return np.sqrt(abs(sc.vev(hh,wf=wf).real-e**2)), e
print("== (A) 8-site S=1/2, ED E0=-3.3749325987; solved and measured at the same maxm")
for v in ("python",3):
    for m in (3,6):
        np.random.seed(1)
        sc = spinchain.Spin_Chain(["S=1/2"]*8,itensor_version=v)
        sc.maxm=m; sc.nsweeps=12
        h = ham(sc); sc.set_hamiltonian(h); hh = h*h
        sc.gs_energy(); wf = sc.get_gs()
        t,e = truth(sc,wf,hh)
        print("v=%-6s maxm=%d E=%.8f  gs_energy_fluctuation()=%.3e  truth=%.3e  ratio truth/reported=%.1f"
              %(v,m,e,sc.gs_energy_fluctuation(),t,t/sc.gs_energy_fluctuation()))
spins = ["S=1"]*10
print("== (B) 10-site S=1 (ED E0=-12.89456013), maxm=10, maxde=1e-3: what the loop read")
for v in (3,"python"):
    np.random.seed(2)
    sc = spinchain.Spin_Chain(spins,itensor_version=v); sc.maxm=10
    h = ham(sc); sc.set_hamiltonian(h); hh = h*h
    with contextlib.redirect_stdout(io.StringIO()): sc.gs_energy(maxde=1e-3)
    wf = sc.get_gs(); t,e = truth(sc,wf,hh)
    print("v=%-6s refined E=%.8f  npow@maxm=20 per site=%.3e (<= maxde, so the loop stopped)  truth per site=%.3e  requested maxde=1e-3"
          %(v,e,fl_npow(sc,wf,20)/10,t/10))
print("== (C) v3, maxde=1e-6 (the example's third point)")
np.random.seed(2)
sc = spinchain.Spin_Chain(spins,itensor_version=3); sc.maxm=10
h = ham(sc); sc.set_hamiltonian(h); hh = h*h
log = io.StringIO(); t0=time.time()
with contextlib.redirect_stdout(log): sc.gs_energy(maxde=1e-6)
inner = [l.split("=")[1].strip() for l in log.getvalue().splitlines() if "Energy fluctuation" in l]
wf = sc.get_gs(); t,e = truth(sc,wf,hh)
print("E=%.8f  loop reads per site (value, maxm): %s  final gs_energy_fluctuation()@maxm=10=%.3e  truth=%.3e  (%.0fs)"
      %(e,inner,sc.gs_energy_fluctuation(),t,time.time()-t0))
```

`reviews/groundstate_4_new1/03_consistent_regime_and_maxde_check.out`:


```
== (A) 8-site S=1/2, ED E0=-3.3749325987; solved and measured at the same maxm
v=python maxm=3 E=-3.27937301  gs_energy_fluctuation()=6.725e-03  truth=3.436e-01  ratio truth/reported=51.1
v=python maxm=6 E=-3.37329763  gs_energy_fluctuation()=3.996e-03  truth=6.726e-02  ratio truth/reported=16.8
v=3      maxm=3 E=-3.27925299  gs_energy_fluctuation()=1.492e-02  truth=3.474e-01  ratio truth/reported=23.3
v=3      maxm=6 E=-3.37328923  gs_energy_fluctuation()=6.795e-03  truth=6.766e-02  ratio truth/reported=10.0
== (B) 10-site S=1 (ED E0=-12.89456013), maxm=10, maxde=1e-3: what the loop read
v=3      refined E=-12.89448122  npow@maxm=20 per site=2.604e-04 (<= maxde, so the loop stopped)  truth per site=2.283e-03  requested maxde=1e-3
v=python refined E=-12.89448263  npow@maxm=20 per site=1.672e-04 (<= maxde, so the loop stopped)  truth per site=2.290e-03  requested maxde=1e-3
== (C) v3, maxde=1e-6 (the example's third point)
E=-12.89456013  loop reads per site (value, maxm): ['0.0018958122840179927 10', '0.00026043811918661966 20', '1.7603494959469614e-05 40']  final gs_energy_fluctuation()@maxm=10=1.694e+00  truth=2.522e-05  (1s)
```


**Suggested fix**: compute the second moment as `self.vev(h*h)` in
`manybodychain.gs_energy_fluctuation` and in the `maxde` loop at
`groundstate.py:378-379`. That builds one MPO of bond dimension at most D_H^2,
needs no MPS truncation, works on every backend including ED, and agrees with the
untruncated route to four digits. Raising `maxm` before measuring would fix only
the over-reporting side, and the under-reporting side is the one the ordinary
workflow hits. For a general `vev(op, npow=2)` the untruncated quantity exists on
both compiled sides (on v3, `innerC(A,x,B,y)` = <Ax|By>, `mpo.h:473`, valid as
<A^2> for Hermitian A, or an application with no `MaxDim` and only the cutoff; on
`pyitensor`, `applyMPO` with `maxdim=None`); for `npow>2`, apply without `MaxDim`
or document that the result truncates at `self.maxm`, and at minimum warn when the
state's bond dimension exceeds it. The user guide's fluctuation paragraph should
say that the number below full bond dimension is set by `maxm`. The example
then needs its plot rechecked. Lands with finding 6. A reviewer-found companion,
that the `maxde` loop compares a per-site fluctuation against `maxde` while
`gs_energy_fluctuation()` returns the total, got no reviewer of its own and is in
"New leads". NUMBERS CHANGE for every `gs_energy_fluctuation()` and
`vev(H, npow=2)` below full bond dimension, on every DMRG backend, by the factors
in the table above, and for every `maxde=` run, which will now stop later.

### 8. `submode="SECTOR"` never reads the state `set_gs` put on the chain, so after `set_gs` of anything but the global ground state it returns the global ground state's spectrum with no warning (1.0485 off the set state's exact curve on a 1.0616 peak, and 2.8339 off for a member that is the lowest of its own sector) while CVM on the same chain returns the set state's to 0.0000

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED &middot; lens `groundstate`

**Where**: `src/dmrgpy/sectordc.py:285-288` (`_sector_states`:
`clone.set_conserved_sector(**reference); clone.gs_energy(); wf0 =
clone.get_gs()`); `dynamics.py:250` (the SECTOR dispatch, with no test for an
injected state); `spinchain.py`'s `_N_GS_SUBMODES` and its raise, the only guard.

SECTOR solves the reference sector and the target sector on an internal clone, and
the clone's ground state is reset by `__deepcopy__` since `867e2b4`. Its reference
sector comes from `sector=`, else the chain's `conserved_sector`, else a charge
measured on the clone's own unconstrained re-solve, so neither the reference charge
nor the |GS> ever comes from the chain's state. `867e2b4` made `set_gs` reach every
other submode and recognised SECTOR's behaviour in the comment above
`_N_GS_SUBMODES`, but guarded only `get_kondo_spectrum(n_gs>1)`, so a plain SECTOR
call after `set_gs` is silent. The user guide now says, at `user_guide.md:4245`,
that a state set by hand reaches the solver, and its SECTOR entry writes |GS> with
no caveat where the EX entry above it names "the chain's own state, the one
`set_gs()` set if any". It survived because
`test_every_submode_measures_the_member_set_gs_set` lists KPM, CVM,
CVM_explicit, ROOTN, TD, TDZ and EX and leaves SECTOR out.

**Expected**: either the correlator of the state the chain holds (|2>, the upper
doublet's Sz=-1/2 member of the 3-site Heisenberg chain in Bz=0.3, E=-0.15), as
every other submode gives, or a refusal.

Repro, from the hunter:

```bash
cd <scratch>/groundstate && <scratch>/run3.sh 09_setgs_sector_submode.py 2>&1 | tee 09_setgs_sector_submode.out
```

`groundstate/09_setgs_sector_submode.py`:

```python
# set_gs(|2>), |2> the second excited eigenstate of the 3-site Heisenberg
# chain in a field (the upper doublet's Sz=-1/2 member, E=-0.15, the same
# Sz sector as the ground state |0>).  Which state does each submode measure?
# Exact Lehmann sums from |2> (own energy) and from |0>, as references.
import numpy as np, io, contextlib, warnings
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore"); return f()
def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, B, delta = 3, 0.3, 0.05
es = np.linspace(-1.5,2.4,781)
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
def dense(m):
    m = ed.get_ED_obj().MO2matrix(m)
    return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
E,V = np.linalg.eigh(dense(ham(ed,B))); Sz0 = dense(ed.Sz[0])
def lehmann(k):
    M = np.abs(V.conj().T@Sz0@V[:,k])**2; D = E-E[k]
    return sum(M[m]*delta/np.pi/((es-D[m])**2+delta**2) for m in range(len(E)))
ref2, ref0 = lehmann(2), lehmann(0)
print("E0..3 exact:",np.round(E[:4],6))
for v in ("python",3):
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,B)); sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=3))
    x = ww[2]
    print("v=%s E(|2>) from DMRG = %.6f, <Sz_tot> = %.4f"%(v,ee[2],(sc.aMb(x,sum(sc.Sz),x)/sc.overlap(x,x)).real))
    for sub in ("CVM","SECTOR"):
        sc.set_gs(x)
        kw = dict(delta=delta,es=es)
        if sub=="SECTOR": kw["nex"]=4
        y = quiet(lambda: sc.get_dynamical_correlator(submode=sub,name=(sc.Sz[0],sc.Sz[0]),**kw))[1]
        print("v=%-6s submode=%-6s max|C - C_exact(|2>)| = %.4f   max|C - C_exact(|0>)| = %.4f   (peak %.4f)"
              %(v,sub,np.max(np.abs(np.real(y)-ref2)),np.max(np.abs(np.real(y)-ref0)),np.max(ref2)))
```

`groundstate/09_setgs_sector_submode.out`:


```
E0..3 exact: [-1.15 -0.85 -0.15  0.05]
v=python E(|2>) from DMRG = -0.150000, <Sz_tot> = -0.5000
v=python submode=CVM    max|C - C_exact(|2>)| = 0.0000   max|C - C_exact(|0>)| = 1.0485   (peak 1.0616)
v=python submode=SECTOR max|C - C_exact(|2>)| = 1.0485   max|C - C_exact(|0>)| = 0.0000   (peak 1.0616)
v=3 E(|2>) from DMRG = -0.150000, <Sz_tot> = -0.5000
v=3      submode=CVM    max|C - C_exact(|2>)| = 0.0000   max|C - C_exact(|0>)| = 1.0485   (peak 1.0616)
v=3      submode=SECTOR max|C - C_exact(|2>)| = 1.0485   max|C - C_exact(|0>)| = 0.0000   (peak 1.0616)
```


**Reviewer (CONFIRMED)**: reproduces to every printed digit on both backends; the
anchor is an exact Lehmann sum over `np.linalg.eigh` of the ED matrix, the
injected level is non-degenerate, and a 3-site chain at `maxm=20` is exact. The
chain is left undisturbed (after SECTOR, |<2|gs>|^2 = 1.000000 and a following CVM
still measures |2>), so the wrong number is confined to the SECTOR curve. With
(S-_0,S+_0), whose poles differ between the two doublet members, `set_gs(|1>)`,
the Sz=+1/2 member and the lowest state of its own sector, then SECTOR gives 0.0000
from |0> and 2.8339 from |1> on both backends, identical to SECTOR with no
`set_gs` at all, while CVM gives |1> to 0.0000; passing `sector={"Sz":1}` gives
|1>, but only because |1> is the lowest of that sector, and after `set_gs(|2>)` no
`sector=` recovers |2>. The origin is `ed00449`; `867e2b4` is where it became a
mismatch. The mark cannot be the guard: `pending_injection(sc)` is not None right
after `set_gs`, but is None after one CVM call, before SECTOR is reached, while the
chain still holds the injected state. Struck, each explicitly:

- The hunter's `why_tests_pass` claim that when the injected state is the lowest
  of its own sector SECTOR's re-solve lands on it: it does not, since the
  reference sector is measured on the clone's unconstrained re-solve; a doublet
  test would see the defect as long as it does not pass `sector=`.
- The first form of the suggested fix, refusing when `pending_injection(self)` is
  not None at the SECTOR branch: it cannot fire, since `ground_state_on_session`
  runs before the dispatch and `_take_injected_state` clears the mark, and a check
  moved above it misses `set_gs` -> CVM -> SECTOR.

The reviewer puts it closer to LOW: the wrong answer is always the global ground
state's correct spectrum, and `set_gs` followed by SECTOR is unusual; it is silent,
and the user guide tells the reader every DMRG correlator honours `set_gs`.

```bash
cd <scratch>/reviews/groundstate_5 && <scratch>/run3.sh 01_hunter_repro.py 2>&1 | tee 01_hunter_repro.out
```

`reviews/groundstate_5/01_hunter_repro.py` is `groundstate/09_setgs_sector_submode.py` (finding 8) unchanged, rerun; its output, `reviews/groundstate_5/01_hunter_repro.out`:


```
E0..3 exact: [-1.15 -0.85 -0.15  0.05]
v=python E(|2>) from DMRG = -0.150000, <Sz_tot> = -0.5000
v=python submode=CVM    max|C - C_exact(|2>)| = 0.0000   max|C - C_exact(|0>)| = 1.0485   (peak 1.0616)
v=python submode=SECTOR max|C - C_exact(|2>)| = 1.0485   max|C - C_exact(|0>)| = 0.0000   (peak 1.0616)
v=3 E(|2>) from DMRG = -0.150000, <Sz_tot> = -0.5000
v=3      submode=CVM    max|C - C_exact(|2>)| = 0.0000   max|C - C_exact(|0>)| = 1.0485   (peak 1.0616)
v=3      submode=SECTOR max|C - C_exact(|2>)| = 1.0485   max|C - C_exact(|0>)| = 0.0000   (peak 1.0616)
```

```bash
cd <scratch>/reviews/groundstate_5 && <scratch>/run3.sh 02_sector_attacks.py 2>&1 | tee 02_sector_attacks.out
```

`reviews/groundstate_5/02_sector_attacks.py`:

```python
# Attacks on groundstate_5.  Same chain as the hunter (3-site Heisenberg, Bz=0.3).
# (a) set_gs(|1>), the Sz=+1/2 ground doublet member, lowest of its OWN sector:
#     the hunter says SECTOR's re-solve would land on it; does it?
# (b) with sector={"Sz":+1} passed explicitly after set_gs(|1>).
# (c) after set_gs(|2>) + SECTOR, is the chain's state still |2> (no side effect),
#     i.e. does a following CVM without a fresh set_gs still give |2>'s spectrum?
import numpy as np, io, contextlib, warnings
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore"); return f()
def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, B, delta = 3, 0.3, 0.05
es = np.linspace(-1.5,2.4,781)
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
def dense(m):
    m = ed.get_ED_obj().MO2matrix(m)
    return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
E,V = np.linalg.eigh(dense(ham(ed,B))); Sz0 = dense(ed.Sz[0]); Szt = dense(sum(ed.Sz))
print("E exact:",np.round(E,6))
print("Sz_tot exact:",np.round([ (V[:,k].conj()@Szt@V[:,k]).real for k in range(len(E))],3))
def lehmann(k):
    M = np.abs(V.conj().T@Sz0@V[:,k])**2; D = E-E[k]
    return sum(M[m]*delta/np.pi/((es-D[m])**2+delta**2) for m in range(len(E)))
refs = {k:lehmann(k) for k in range(3)}
def report(tag,y):
    d = {k:np.max(np.abs(np.real(y)-refs[k])) for k in refs}
    print("  %-40s dev from |0>: %.4f  |1>: %.4f  |2>: %.4f"%(tag,d[0],d[1],d[2]))
for v in ("python",3):
    print("v=%s"%v)
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,B)); sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=3))
    for k in range(3):
        x=ww[k]; print("  DMRG member %d: E=%.6f Sz_tot=%.4f"%(k,ee[k],(sc.aMb(x,sum(sc.Sz),x)/sc.overlap(x,x)).real))
    kw = dict(delta=delta,es=es)
    dc = lambda sub,**e: quiet(lambda: sc.get_dynamical_correlator(submode=sub,name=(sc.Sz[0],sc.Sz[0]),**kw,**e))[1]
    # (a)
    sc.set_gs(ww[1])
    report("(a) set_gs(|1>) CVM", dc("CVM"))
    sc.set_gs(ww[1])
    report("(a) set_gs(|1>) SECTOR", dc("SECTOR",nex=4))
    # (b)
    sc.set_gs(ww[1])
    report("(b) set_gs(|1>) SECTOR sector={Sz:+1}", dc("SECTOR",nex=4,sector={"Sz":1}))
    # (c)
    sc.set_gs(ww[2])
    report("(c) set_gs(|2>) SECTOR", dc("SECTOR",nex=4))
    g = sc.get_gs()
    print("  (c) after SECTOR: |<|2>|chain gs>|^2 = %.6f, gs_energy() = %.6f"
          %(abs(sc.overlap(ww[2],g))**2/abs(sc.overlap(g,g)), sc.gs_energy()))
    report("(c) then CVM with no new set_gs", dc("CVM"))
```

`reviews/groundstate_5/02_sector_attacks.out`:


```
E exact: [-1.15 -0.85 -0.15  0.05  0.15  0.35  0.65  0.95]
Sz_tot exact: [-0.5  0.5 -0.5 -1.5  0.5 -0.5  0.5  1.5]
v=python
  DMRG member 0: E=-1.150000 Sz_tot=-0.5000
  DMRG member 1: E=-0.850000 Sz_tot=0.5000
  DMRG member 2: E=-0.150000 Sz_tot=-0.5000
  (a) set_gs(|1>) CVM                      dev from |0>: 0.0000  |1>: 0.0000  |2>: 1.0485
  (a) set_gs(|1>) SECTOR                   dev from |0>: 0.0000  |1>: 0.0000  |2>: 1.0485
  (b) set_gs(|1>) SECTOR sector={Sz:+1}    dev from |0>: 0.0000  |1>: 0.0000  |2>: 1.0485
  (c) set_gs(|2>) SECTOR                   dev from |0>: 0.0000  |1>: 0.0000  |2>: 1.0485
  (c) after SECTOR: |<|2>|chain gs>|^2 = 1.000000, gs_energy() = -0.150000
  (c) then CVM with no new set_gs          dev from |0>: 1.0485  |1>: 1.0485  |2>: 0.0000
v=3
  DMRG member 0: E=-1.150000 Sz_tot=-0.5000
  DMRG member 1: E=-0.850000 Sz_tot=0.5000
  DMRG member 2: E=-0.150000 Sz_tot=-0.5000
  (a) set_gs(|1>) CVM                      dev from |0>: 0.0000  |1>: 0.0000  |2>: 1.0485
  (a) set_gs(|1>) SECTOR                   dev from |0>: 0.0000  |1>: 0.0000  |2>: 1.0485
  (b) set_gs(|1>) SECTOR sector={Sz:+1}    dev from |0>: 0.0000  |1>: 0.0000  |2>: 1.0485
  (c) set_gs(|2>) SECTOR                   dev from |0>: 0.0000  |1>: 0.0000  |2>: 1.0485
  (c) after SECTOR: |<|2>|chain gs>|^2 = 1.000000, gs_energy() = -0.150000
  (c) then CVM with no new set_gs          dev from |0>: 1.0485  |1>: 1.0485  |2>: 0.0000
```

```bash
cd <scratch>/reviews/groundstate_5 && <scratch>/run3.sh 03_sector_other_sector.py 2>&1 | tee 03_sector_other_sector.out
```

`reviews/groundstate_5/03_sector_other_sector.py`:

```python
# Probe 02's (a)/(b) could not tell |0> from |1>: Sz0-Sz0 is spin-flip invariant
# and the Zeeman term does not move excitation energies inside a sector, so both
# doublet members give the same curve.  Redo them with (Sm0,Sp0), which from |0>
# (Sz=-1/2) lands in Sz=+1/2 and from |1> (Sz=+1/2) lands in Sz=+3/2: different
# poles.  Question: after set_gs(|1>), the lowest state of its OWN sector, does
# SECTOR (default reference sector) measure |1> (as the hunter's why_tests_pass
# predicts) or |0>?  And does sector={"Sz":1} recover |1>?
import numpy as np, io, contextlib, warnings
from dmrgpy import spinchain
from dmrgpy.operatornames import name2MO
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore"); return f()
def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, B, delta = 3, 0.3, 0.05
es = np.linspace(-1.5,2.4,781)
ed = spinchain.Spin_Chain(["S=1/2"]*n); ed.set_hamiltonian(ham(ed,B))
def dense(m):
    m = ed.get_ED_obj().MO2matrix(m)
    return np.array(m.todense()) if hasattr(m,"todense") else np.array(m)
E,V = np.linalg.eigh(dense(ham(ed,B)))
Spd = dense(name2MO("Sp",ed)[0]); Smd = dense(name2MO("Sm",ed)[0])
def lehmann(k):
    M = (V[:,k].conj()@Smd@V) * (V.conj().T@Spd@V[:,k]); D = E-E[k]
    return np.real(sum(M[m]*delta/np.pi/((es-D[m])**2+delta**2) for m in range(len(E))))
refs = {k:lehmann(k) for k in range(3)}
print("peaks:", {k:round(float(np.max(np.abs(refs[k]))),4) for k in refs},
      " max|ref0-ref1| = %.4f"%np.max(np.abs(refs[0]-refs[1])))
def report(tag,y):
    d = {k:np.max(np.abs(np.real(y)-refs[k])) for k in refs}
    print("  %-44s dev from |0>: %.4f  |1>: %.4f  |2>: %.4f"%(tag,d[0],d[1],d[2]))
for v in ("python",3):
    print("v=%s"%v)
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,B)); sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=3))
    print("  DMRG members E:", np.round(np.real(ee),6))
    Sp, Sm = name2MO("Sp",sc), name2MO("Sm",sc)
    kw = dict(delta=delta,es=es)
    dc = lambda sub,**e: quiet(lambda: sc.get_dynamical_correlator(submode=sub,name=(Sm[0],Sp[0]),**kw,**e))[1]
    report("no set_gs, SECTOR", dc("SECTOR",nex=4))
    sc.set_gs(ww[1]); report("set_gs(|1>) CVM", dc("CVM"))
    sc.set_gs(ww[1]); report("set_gs(|1>) SECTOR", dc("SECTOR",nex=4))
    sc.set_gs(ww[1]); report("set_gs(|1>) SECTOR sector={Sz:+1}", dc("SECTOR",nex=4,sector={"Sz":1}))
    sc.set_gs(ww[2]); report("set_gs(|2>) CVM", dc("CVM"))
    sc.set_gs(ww[2]); report("set_gs(|2>) SECTOR", dc("SECTOR",nex=4))
```

`reviews/groundstate_5/03_sector_other_sector.out`:


```
peaks: {0: 2.8351, 1: 1.061, 2: 2.1232}  max|ref0-ref1| = 2.8339
v=python
  DMRG members E: [-1.15 -0.85 -0.15]
  no set_gs, SECTOR                            dev from |0>: 0.0000  |1>: 2.8339  |2>: 2.8193
  set_gs(|1>) CVM                              dev from |0>: 2.8339  |1>: 0.0000  |2>: 2.1228
  set_gs(|1>) SECTOR                           dev from |0>: 0.0000  |1>: 2.8339  |2>: 2.8193
  set_gs(|1>) SECTOR sector={Sz:+1}            dev from |0>: 2.8339  |1>: 0.0000  |2>: 2.1228
  set_gs(|2>) CVM                              dev from |0>: 2.8193  |1>: 2.1228  |2>: 0.0000
  set_gs(|2>) SECTOR                           dev from |0>: 0.0000  |1>: 2.8339  |2>: 2.8193
v=3
  DMRG members E: [-1.15 -0.85 -0.15]
  no set_gs, SECTOR                            dev from |0>: 0.0000  |1>: 2.8339  |2>: 2.8193
  set_gs(|1>) CVM                              dev from |0>: 2.8339  |1>: 0.0000  |2>: 2.1228
  set_gs(|1>) SECTOR                           dev from |0>: 0.0000  |1>: 2.8339  |2>: 2.8193
  set_gs(|1>) SECTOR sector={Sz:+1}            dev from |0>: 2.8339  |1>: 0.0000  |2>: 2.1228
  set_gs(|2>) CVM                              dev from |0>: 2.8193  |1>: 2.1228  |2>: 0.0000
  set_gs(|2>) SECTOR                           dev from |0>: 0.0000  |1>: 2.8339  |2>: 2.8193
```

```bash
cd <scratch>/reviews/groundstate_5 && <scratch>/run3.sh 04_injection_mark_taken.py 2>&1 | tee 04_injection_mark_taken.out
```

`reviews/groundstate_5/04_injection_mark_taken.py`:

```python
# Can the hunter's first fix option (refuse when groundstate.pending_injection(self)
# is not None) see the case?  set_gs(|2>), then CVM, then SECTOR with no new set_gs:
# print the mark just before SECTOR.  Also: the data a session-local check would
# need, the chain's own <Sz_tot> and e0 on its own session, against SECTOR's
# reference-sector energy (info from return_poles is not used; e0 of |0> is exact -1.15).
import numpy as np, io, contextlib, warnings
from dmrgpy import spinchain, groundstate
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore"); return f()
def ham(sc,B):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + B*sc.Sz[i]
    return h
n, B, delta = 3, 0.3, 0.05
es = np.linspace(-1.5,2.4,781)
for v in ("python",3):
    np.random.seed(7)
    sc = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    sc.maxm=20; sc.nsweeps=12
    sc.set_hamiltonian(ham(sc,B)); sc.gs_energy()
    ee,ww = quiet(lambda: sc.get_excited_states(n=3))
    sc.set_gs(ww[2])
    print("v=%-6s right after set_gs: pending_injection is None? %s"%(v, groundstate.pending_injection(sc) is None))
    quiet(lambda: sc.get_dynamical_correlator(submode="CVM",name=(sc.Sz[0],sc.Sz[0]),delta=delta,es=es))
    print("v=%-6s after CVM, before SECTOR: pending_injection is None? %s"%(v, groundstate.pending_injection(sc) is None))
    print("v=%-6s chain e0 = %.6f, chain <Sz_tot> = %.4f (session-local, no cross-session overlap)"
          %(v, sc.gs_energy(), complex(sc.vev(sum(sc.Sz))).real))
```

`reviews/groundstate_5/04_injection_mark_taken.out`:


```
v=python right after set_gs: pending_injection is None? False
v=python after CVM, before SECTOR: pending_injection is None? True
v=python chain e0 = -0.150000, chain <Sz_tot> = -0.5000 (session-local, no cross-session overlap)
v=3      right after set_gs: pending_injection is None? False
v=3      after CVM, before SECTOR: pending_injection is None? True
v=3      chain e0 = -0.150000, chain <Sz_tot> = -0.5000 (session-local, no cross-session overlap)
```


**Suggested fix**: the chain needs a persistent record that its current state was
supplied by the caller rather than solved: set when `_take_injected_state` hands
the state over, retired by the next real solve, by `restart()` and by
`set_hamiltonian()`, and added to the `n_gs` snapshot tuple in `spinchain.py` so
the Kondo `finally` restores it. The check belongs in `sectordc._prepare`, next to
`_check_backend`, so it covers direct `sectordc` callers too. Rather than a blanket
refusal, the reviewer proposes a session-local test that returns the right curve
where one exists: take the reference charges from the chain's own state (`vev` of
the total charge) rather than the clone's re-solve, compare the chain's `e0`
(<wf|H|wf>) with the clone's reference-sector ground energy, and proceed when they
agree to a tolerance and that ground state is non-degenerate, raising otherwise.
That gives the right curve for `set_gs(|1>)`, refuses `set_gs(|2>)`, and saves the
clone's unconstrained solve on the ordinary path. The user guide's sentence at
`:4245` and the SECTOR entry should say which state SECTOR measures. No number
moves under the refusal; under the charge-from-the-chain test, a set lowest-of-its-
sector state gets its own curve instead of the global ground state's.

### 9. The band-edge pre-fill `867e2b4` put in `_take_injected_state` turns the first read after `set_gs`/`set_initial_wf` into an upper-edge -H solve, plus a full solve on an unsolved chain, plus an energy fluctuation it throws away: 7.1 s on `"python"` and 0.42 s on v3 on a 24-site chain, against 0.037 s and 0.010 s for the <x|H|x> the read returns

`performance` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `groundstate`

**Where**: `src/dmrgpy/groundstate.py:195-198` (`_take_injected_state`'s
Hermitian `session.excited_states(1,1.0,False)`); `pyitensor/chain.py:672-697`
(`excited_states`: `gs_energy()` when `wf0` is None, `_bandwidth`, and
`_energy_fluctuation` for the one state).

The pre-fill is deliberate: `_take_injected_state`'s docstring records that
without it the first KPM call after a push sweeps the pushed state in place (to
overlap 0.989 on `"python"` and 0.987 on v3), because the session fills its lower
band edge lazily through `gs_energy(skip_dmrg=True)`. What the docstring and the
second pass's Status line misstate is the cost. The docstring says it "costs
nothing when the edges are cached already", while the discarded energy fluctuation
remains and is paid on every injection; the Status line says "one reduced solve
when the edges are missing", while on `"python"` that reduced solve is four to five
full solves (finding 10) and an unsolved chain adds a full solve on top. A solve
caches only the lower edge, so a solved chain that never ran KPM or excited states
pays the upper edge too. It survived because the `867e2b4` tests inject on 3-site
chains, where every piece is too cheap to see, and the only number they pin is the
energy, which the pre-fill does not touch.

**Expected**: the read returns <x|H|x> (it does, exactly) at roughly the cost of
computing it.

Repro, from the hunter:

```bash
cd <scratch>/groundstate && <scratch>/run3.sh 06_fill_cost.py 2>&1 | tee 06_fill_cost.out
```

`groundstate/06_fill_cost.py`:

```python
# Cost of taking an injected state "as it is": set_initial_wf(x) (or
# set_gs(x)) then gs_energy() on a chain whose session holds no energy.
# The documented work is one <x|H|x>; what runs is the band-edge fill.
import numpy as np, io, contextlib, time
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
n = 24
for v in ("python",3):
    ts_take, ts_amb, ts_solve, ts_second = [],[],[],[]
    for rep in range(3):
        np.random.seed(10+rep)
        src = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
        src.maxm=30; src.nsweeps=6
        src.set_hamiltonian(ham(src))
        x = src.random_mps()
        t=time.perf_counter(); ex = (src.aMb(x,src.hamiltonian,x)/src.overlap(x,x)).real; ts_amb.append(time.perf_counter()-t)
        src.set_initial_wf(x)       # fresh chain: the session has never solved
        t=time.perf_counter(); e = quiet(lambda: src.gs_energy()); ts_take.append(time.perf_counter()-t)
        x2 = src.random_mps()
        src.set_gs(x2)              # second injection, same chain, edges cached now
        t=time.perf_counter(); quiet(lambda: src.gs_energy()); ts_second.append(time.perf_counter()-t)
        ref = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
        ref.maxm=30; ref.nsweeps=6; ref.set_hamiltonian(ham(ref))
        t=time.perf_counter(); quiet(lambda: ref.gs_energy()); ts_solve.append(time.perf_counter()-t)
    print("v=%-6s n=%d: gs_energy() after set_initial_wf(x) returned %.6f (= <x|H|x> %.6f); median time %.3f s, "
          "against %.4f s for <x|H|x> itself and %.3f s for a full ground-state solve; a second set_gs(x2)+gs_energy() %.3f s"
          %(v,n,e,ex,np.median(ts_take),np.median(ts_amb),np.median(ts_solve),np.median(ts_second)))
```

`groundstate/06_fill_cost.out`:


```
v=python n=24: gs_energy() after set_initial_wf(x) returned 0.100309 (= <x|H|x> 0.100309); median time 8.418 s, against 0.0322 s for <x|H|x> itself and 1.029 s for a full ground-state solve; a second set_gs(x2)+gs_energy() 0.181 s
v=3      n=24: gs_energy() after set_initial_wf(x) returned 0.008513 (= <x|H|x> 0.008513); median time 0.469 s, against 0.0118 s for <x|H|x> itself and 0.249 s for a full ground-state solve; a second set_gs(x2)+gs_energy() 0.092 s
```


**Reviewer (CONFIRMED, NARROWED)**: re-measured with a breakdown on the 24-site
Heisenberg chain (`maxm=30`, `nsweeps=6`, two seeds): on `"python"` the injected
read is a median 7.066 s against 0.0370 s for <x|H|x>, a second injection 0.179 s,
and the pieces on a fresh session are send H 0.022 s, full solve 0.899 s,
upper-edge -H solve 4.309 s, energy fluctuation 0.139 s; on v3 the read is 0.419 s
against 0.0103 s, the second injection 0.082 s, and the pieces 0.002, 0.222, 0.161
and 0.024 s. The returned energy equals <x|H|x> in every run. A `"python"` chain
solved first (1.442 s) that never ran KPM still pays 7.787 s on its first
injection, v3 0.202 s against a 0.261 s solve; after a KPM call has filled both
edges an injection costs 0.228 s on `"python"` and 0.036 s on v3, against 0.031 s
and 0.0106 s for <x|H|x>, which is the discarded fluctuation. A profile of
`_maximum_energy` led to finding 10. Struck, each explicitly:

- "260x the documented <x|H|x>" as a broken promise: `set_initial_wf`'s "Nothing
  is computed here" describes the setter, which is true; what is false is the
  private docstring's "costs nothing when the edges are cached already" and the
  Status line's "one reduced solve".
- "On a chain whose session holds no energy" as the trigger: a solved chain that
  never ran KPM pays essentially the same on `"python"` (7.79 s).
- The hunter's `why_tests_pass` claim that on the tests' pre-solved chains "the
  band edges are already cached": a solve caches only the lower edge.
- The exact figures 8.4 s and 0.47 s: replaced by the re-measured 7.1 s and 0.42 s.

```bash
cd <scratch>/reviews/groundstate_6 && <scratch>/run3.sh 01_breakdown.py 2>&1 | tee 01_breakdown.out
```

`reviews/groundstate_6/01_breakdown.py`:

```python
# Reviewer probe for groundstate_6: reproduce the hunter's timing and split
# the fill (session.excited_states(1,1.0,False)) into its parts on a fresh
# 24-site Heisenberg chain: send H, full solve, upper edge (-H solve) and
# energy fluctuation, against <x|H|x>.
import numpy as np, io, contextlib, time
from dmrgpy import spinchain, groundstate
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
def T(f):
    t=time.perf_counter(); r=quiet(f); return time.perf_counter()-t, r
n = 24
for v in ("python",3):
    rows = []
    for rep in range(2):
        np.random.seed(10+rep)
        # (a) the hunter's route, end to end
        a = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
        a.maxm=30; a.nsweeps=6; a.set_hamiltonian(ham(a))
        x = a.random_mps()
        t_amb,ex = T(lambda: (a.aMb(x,a.hamiltonian,x)/a.overlap(x,x)).real)
        a.set_initial_wf(x)
        t_take,e = T(lambda: a.gs_energy())
        x2 = a.random_mps(); a.set_gs(x2)
        t_second,_ = T(lambda: a.gs_energy())
        # (b) the same pieces one by one on a fresh chain's session
        b = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
        b.maxm=30; b.nsweeps=6; b.set_hamiltonian(ham(b))
        groundstate._session_parameters(b)
        t_send,_ = T(lambda: groundstate.send_hamiltonian(b))
        t_gs,_   = T(lambda: b._session.gs_energy())
        t_edge,_ = T(lambda: b._session.excited_states(1,1.0,False)) # upper edge + fluctuation
        t_fluc,_ = T(lambda: b._session.excited_states(1,1.0,False)) # fluctuation only (edges cached)
        rows.append((t_take,t_amb,t_second,t_send,t_gs,t_edge-t_fluc,t_fluc))
        print("v=%-6s rep=%d e=%.6f <x|H|x>=%.6f equal=%s"%(v,rep,e,ex,abs(e-ex)<1e-12))
    m = np.median(np.array(rows),axis=0)
    print("v=%-6s n=%d medians: take-injected %.3f s | <x|H|x> %.4f s | second injection %.3f s"%(v,n,m[0],m[1],m[2]))
    print("          parts: send H %.3f s, full solve %.3f s, upper-edge -H solve %.3f s, energy fluctuation %.3f s (sum %.3f s)"
          %(m[3],m[4],m[5],m[6],m[3]+m[4]+m[5]+m[6]))
```

`reviews/groundstate_6/01_breakdown.out`:


```
v=python rep=0 e=-0.167349 <x|H|x>=-0.167349 equal=True
v=python rep=1 e=0.060280 <x|H|x>=0.060280 equal=True
v=python n=24 medians: take-injected 7.066 s | <x|H|x> 0.0370 s | second injection 0.179 s
          parts: send H 0.022 s, full solve 0.899 s, upper-edge -H solve 4.309 s, energy fluctuation 0.139 s (sum 5.369 s)
v=3      rep=0 e=0.053183 <x|H|x>=0.053183 equal=True
v=3      rep=1 e=-0.081149 <x|H|x>=-0.081149 equal=True
v=3      n=24 medians: take-injected 0.419 s | <x|H|x> 0.0103 s | second injection 0.082 s
          parts: send H 0.002 s, full solve 0.222 s, upper-edge -H solve 0.161 s, energy fluctuation 0.024 s (sum 0.410 s)
```

```bash
cd <scratch>/reviews/groundstate_6 && <scratch>/run3.sh 02_upper_edge_profile.py 2>&1 | tee 02_upper_edge_profile.out
```

`reviews/groundstate_6/02_upper_edge_profile.py`:

```python
# Reviewer probe for groundstate_6: why the "reduced-effort" upper-edge
# solve (pyitensor Chain._maximum_energy, 5 sweeps at maxdim<=20) costs
# several times a full 6-sweep maxm=30 ground-state solve on "python".
import numpy as np, io, contextlib, time, cProfile, pstats
from dmrgpy import spinchain, groundstate
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
n = 24
np.random.seed(10)
b = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version="python")
b.maxm=30; b.nsweeps=6; b.set_hamiltonian(ham(b))
groundstate._session_parameters(b); groundstate.send_hamiltonian(b)
s = b._session
t=time.perf_counter(); e0=s.gs_energy(); t_gs=time.perf_counter()-t
pr = cProfile.Profile()
t=time.perf_counter(); pr.enable(); emax=s._maximum_energy(); pr.disable(); t_max=time.perf_counter()-t
print("full solve %.3f s (E0=%.6f); upper edge %.3f s (Emax=%.6f, exact ferromagnetic top %.6f)"%(t_gs,e0,t_max,emax,0.25*(n-1)))
st = pstats.Stats(pr); st.sort_stats("cumulative"); st.print_stats(18)
```

`reviews/groundstate_6/02_upper_edge_profile.out`:


```
full solve 1.995 s (E0=-10.453786); upper edge 7.543 s (Emax=5.750000, exact ferromagnetic top 5.750000)
         1469260 function calls in 7.540 seconds

   Ordered by: cumulative time
   List reduced from 186 to 18 due to restriction <18>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    7.543    7.543 <repo>/src/dmrgpy/pyitensor/chain.py:1818(_maximum_energy)
        1    0.000    0.000    7.502    7.502 <repo>/src/dmrgpy/pyitensor/dmrg.py:476(dmrg)
        5    0.003    0.001    7.502    1.500 <repo>/src/dmrgpy/pyitensor/dmrg.py:446(_dmrg_one_sweep)
      230    0.004    0.000    7.197    0.031 <repo>/src/dmrgpy/pyitensor/dmrg.py:389(_local_ground_state)
      230    3.199    0.014    7.081    0.031 <repo>/src/dmrgpy/pyitensor/dmrg.py:118(_lanczos_ground_state)
     6999    0.027    0.000    1.914    0.000 <repo>/src/dmrgpy/pyitensor/dmrg.py:73(_tridiag_ground_value)
    15433    0.018    0.000    1.688    0.000 <repo>/src/dmrgpy/pyitensor/backend.py:331(wrapper)
     6999    1.450    0.000    1.552    0.000 <python>/lib/python3.14/site-packages/numpy/linalg/_linalg.py:1264(eigvalsh)
     6999    0.007    0.000    1.525    0.000 <repo>/src/dmrgpy/pyitensor/kernels.py:426(matvec)
     6999    1.174    0.000    1.508    0.000 <repo>/src/dmrgpy/pyitensor/kernels.py:433(_matvec_chain_impl)
     7229    0.123    0.000    0.341    0.000 <repo>/src/dmrgpy/pyitensor/dmrg.py:63(_tridiag_matrix)
    76933    0.327    0.000    0.327    0.000 {method 'reshape' of 'numpy.ndarray' objects}
    21247    0.086    0.000    0.159    0.000 <python>/lib/python3.14/site-packages/numpy/lib/_twodim_base_impl.py:258(diag)
      369    0.010    0.000    0.142    0.000 <repo>/src/dmrgpy/pyitensor/tensor.py:337(contract_many)
      230    0.003    0.000    0.132    0.001 <repo>/src/dmrgpy/pyitensor/dmrg.py:427(_apply_local_update)
      253    0.013    0.000    0.128    0.001 <repo>/src/dmrgpy/pyitensor/svd.py:189(svd)
      322    0.112    0.000    0.118    0.000 <python>/lib/python3.14/site-packages/numpy/linalg/_linalg.py:1514(eigh)
      230    0.001    0.000    0.112    0.000 <repo>/src/dmrgpy/pyitensor/dmrg.py:281(_extend_right)
```

```bash
cd <scratch>/reviews/groundstate_6 && <scratch>/run3.sh 03_solved_chain.py 2>&1 | tee 03_solved_chain.out
```

`reviews/groundstate_6/03_solved_chain.py`:

```python
# Reviewer probe for groundstate_6: the fill on a chain that WAS solved
# first (session holds its energy) but never ran KPM, so only the upper
# band edge is missing; and on a chain whose edges a KPM call filled.
import numpy as np, io, contextlib, time
from dmrgpy import spinchain
def quiet(f):
    with contextlib.redirect_stdout(io.StringIO()): return f()
def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
def T(f):
    t=time.perf_counter(); r=quiet(f); return time.perf_counter()-t, r
n = 24
for v in ("python",3):
    np.random.seed(20)
    a = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    a.maxm=30; a.nsweeps=6; a.set_hamiltonian(ham(a))
    t_solve,e0 = T(lambda: a.gs_energy())
    x = a.random_mps()
    t_amb,ex = T(lambda: (a.aMb(x,a.hamiltonian,x)/a.overlap(x,x)).real)
    a.set_gs(x); t1,e1 = T(lambda: a.gs_energy())
    x2 = a.random_mps(); a.set_gs(x2); t2,_ = T(lambda: a.gs_energy())
    print("v=%-6s solved chain (solve %.3f s, E0=%.5f), never ran KPM: first set_gs+gs_energy %.3f s (e=%.6f, <x|H|x>=%.6f), second %.3f s, <x|H|x> alone %.4f s"
          %(v,t_solve,e0,t1,e1,ex,t2,t_amb))
    # a chain whose edges a KPM call has already filled
    b = spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=v)
    b.maxm=30; b.nsweeps=6; b.set_hamiltonian(ham(b))
    quiet(lambda: b.gs_energy())
    t_kpm,_ = T(lambda: b.get_dynamical_correlator(name=(b.Sz[0],b.Sz[0]),es=np.linspace(-0.5,3,20),delta=0.5,submode="KPM"))
    y = b.random_mps(); b.set_gs(y); t3,_ = T(lambda: b.gs_energy())
    print("v=%-6s after a KPM call (%.3f s): set_gs+gs_energy %.3f s"%(v,t_kpm,t3))
```

`reviews/groundstate_6/03_solved_chain.out`:


```
v=python solved chain (solve 1.442 s, E0=-10.45379), never ran KPM: first set_gs+gs_energy 7.787 s (e=0.089277, <x|H|x>=0.089277), second 0.453 s, <x|H|x> alone 0.0309 s
v=python after a KPM call (37.988 s): set_gs+gs_energy 0.228 s
v=3      solved chain (solve 0.261 s, E0=-10.45379), never ran KPM: first set_gs+gs_energy 0.202 s (e=-0.032464, <x|H|x>=-0.032464), second 0.162 s, <x|H|x> alone 0.0106 s
v=3      after a KPM call (3.821 s): set_gs+gs_energy 0.036 s
```


**Suggested fix**: the hunter's second option, a `set_wavefunction(wf, energy)`
that keeps the energy so no pre-fill is needed, is wrong: on `"python"`
`_minimum_energy` is `gs_energy(skip_dmrg=True)`, which returns `_wf0_energy` when
set, and v3's `minimum_energy` reads the same cache, so the lower band edge would
become <x|H|x>, about 0 for a random x against E0 = -10.45 here, the KPM window
would exclude the spectrum below x's energy where x has weight, and the moment
guard would raise. The docstring's "the edges come out as the Hamiltonian's, not as
the injected state's energy" is a requirement. Filling the edges only on the routes
that read them has the right shape but is fragile, since every session method that
reaches `minimum_energy()`/`bandwidth()` internally would need its own pre-fill.
The clean fix is to decouple the lower-edge fill from `wf0`, the way the upper edge
already is: `_minimum_energy` and `Chain::minimum_energy` use the cached energy
when the session has one and otherwise run a reduced solve on a fresh
`default_mps()`, never sweeping `wf0`; then no pre-fill is needed, the cost lands
on the first KPM call, and the discarded fluctuation goes with `excited_states(1)`.
That needs a C++ edit in both `chain_session.h` and a rebuild. The docstring's cost
accounting and the Status line should be corrected either way. Finding 10 is most
of the `"python"` cost. No number changes.

### 10. On `"python"`, the reduced-effort upper band edge `_maximum_energy()` costs 5 to 9 full ground-state solves on a 24-site Heisenberg chain (5.3 to 8.3 s against 0.90 to 1.0 s; v3's costs 0.75 of its own), because the value-criterion Lanczos stalls at the SU(2)-symmetric ferromagnetic top, up to 193 matvecs per local solve where 20 give Emax to 1e-10

`performance` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `groundstate` (turned up by the reviewer of finding 9)

**Where**: `src/dmrgpy/pyitensor/chain.py:1818-1832` (`_maximum_energy`, DMRG on
-H at 5 sweeps and `maxdim` 20, cached per Hamiltonian), through
`pyitensor/dmrg.py:118` (`_lanczos_ground_state`) and the `max(niter, 200)` floor
at `dmrg.py:421`.

The upper band edge is the KPM window's top and is meant to be cheap, which is why
it runs a short schedule at a small bond dimension. On a Heisenberg chain the top
of the spectrum is the degenerate ferromagnetic multiplet, and there the local
Lanczos, which stops when the lowest Ritz value stops moving and is floored at 200
iterations, runs far past what the value needs. It is paid once per Hamiltonian:
the first KPM or excited-state call on each point of a `"python"` parameter scan,
and since `867e2b4` the first injection on a chain whose edges are not filled
(finding 9).

**Expected**: a bound solve at most as expensive as a ground-state solve, as on v3
(0.75 of its own), for an Emax exact to 1e-10.

Repro, from the reviewer of finding 9, who found it (both scripts are spliced
under finding 9: `02_upper_edge_profile` profiles `_maximum_energy`, 230 local
solves and 6999 Lanczos iterations; `01_breakdown` times it without the profiler):

`reviews/groundstate_6/02_upper_edge_profile.py` and its output are spliced under finding 9 above.

`reviews/groundstate_6/01_breakdown.py` and its output are spliced under finding 9 above.


**Reviewer (CONFIRMED, NARROWED)**: every probe on a fresh `"python"` session, 24
sites, `maxm=30`, `nsweeps=6`, with a matvec counter wrapped around
`_lanczos_ground_state` inside the probe process. The full solve takes 0.998,
0.984 and 0.904 s (E0 = -10.453786), matvecs per local solve by sweep 20.4, 11.7,
7.0, 4.2, 3.5, 3.5; the upper edge takes 5.690, 8.257 and 7.939 s (Emax = 5.750000,
the exact 0.25*(n-1)), matvecs per local solve by sweep 53.9, 92.8 to 109.4, 7.6 to
21.6, 2.0, 2.0, a ratio of 5.7 to 8.8 per repetition and 6.3 at the minimum of
three. Two controls: +H on the identical 5-sweep, `maxdim=20` schedule from a
`randomMPS(30)` costs 0.842 to 1.108 s, so the schedule is not the cause; and a
Zeeman field that makes the top unique brings the -H solve to 0.268 s (hz=0.5) and
0.242 s (hz=1.0), with Emax exact, so the stall is the SU(2)-symmetric top. With
the default rule, 28 of 230 local solves pass 100 matvecs, the worst is 193, none
reaches the 200 floor. `tol=1e-6` gives 0.78 to 0.84 s with Emax low by 3.4e-5 to
8.2e-5; a cap at 20 matvecs gives 0.95 to 0.98 s with Emax within 6e-11 to 1.3e-10.
On v3 the hunter's 0.7 is confirmed (minimum ratio 0.75). On a solved `"python"`
chain that never ran KPM, the first `set_gs(x)` + `gs_energy()` takes 5.635 s. No
known-issue file, `ROADMAP.md` entry or earlier record covers the band-edge cost.
Struck, each explicitly:

- "Costs 4 to 5 times a full solve": re-measured without a profiler it is 5.7 to
  8.8 per repetition, so the record gives 5 to 9.
- "Every `"python"` KPM call and every excited-state call pays this":
  `_bandwidth_max` is cached, cleared only by `set_hamiltonian` and
  `_forget_everything_built_on_sites`; what survives is the first such call per
  Hamiltonian, once per point in a scan.
- "About 30 Lanczos iterations per local solve": an average that hides the shape,
  53, then 93 to 109, then 8 to 22, then 2 and 2.
- The hunter's untested guess about the degenerate multiplet: now measured as the
  observable (a field removes the stall, +H on the same schedule does not have it);
  whether it is the degeneracy itself or the nearby one-magnon levels was not told
  apart, and the finding does not depend on it.

```bash
cd <scratch>/reviews/groundstate_6_new1 && <scratch>/run3.sh 01_cost_and_iterations.py 2>&1 | tee 01_cost_and_iterations.out
```

`reviews/groundstate_6_new1/01_cost_and_iterations.py`:

```python
# Reviewer of groundstate_6_new1: re-measure the cost of the "python" upper
# band edge (Chain._maximum_energy) against the full ground-state solve,
# min of 3 on fresh sessions, and count Lanczos matvecs per local solve,
# per sweep, for both. Then two controls on the same 24-site chain:
#  (C) the SAME 5-sweep maxdim-20 schedule on +H instead of -H;
#  (B) the upper edge of Heisenberg plus a uniform field hz*sum(Sz), which
#      Zeeman-splits the 25-fold ferromagnetic top into a unique state.
import numpy as np, io, contextlib, time
from dmrgpy import spinchain, groundstate
from dmrgpy.pyitensor import dmrg as D
from dmrgpy.pyitensor.mpsalgebra import randomMPS

rec = {"sweep": -1, "log": []}
_orig_lz = D._lanczos_ground_state
def counting_lz(matvec, v0, **kw):
    c = [0]
    def mv(x):
        c[0] += 1
        return matvec(x)
    out = _orig_lz(mv, v0, **kw)
    rec["log"].append((rec["sweep"], c[0]))
    return out
D._lanczos_ground_state = counting_lz
_orig_sw = D._dmrg_one_sweep
def counting_sw(*a, **k):
    rec["sweep"] += 1
    return _orig_sw(*a, **k)
D._dmrg_one_sweep = counting_sw

def reset():
    rec["sweep"] = -1; rec["log"] = []
def per_sweep():
    lg = rec["log"]
    ns = max(s for s, _ in lg) + 1
    return [np.mean([c for s, c in lg if s == k]) for k in range(ns)], \
           sum(c for _, c in lg), len(lg)

def ham(sc, hz=0.0):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    if hz:
        for i in range(sc.ns):
            h = h + hz*sc.Sz[i]
    return h

def session(n, hz=0.0, seed=0):
    np.random.seed(seed)
    b = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    b.maxm = 30; b.nsweeps = 6; b.set_hamiltonian(ham(b, hz))
    groundstate._session_parameters(b); groundstate.send_hamiltonian(b)
    return b._session

def timed(f):
    reset()
    t = time.perf_counter(); p = time.process_time()
    with contextlib.redirect_stdout(io.StringIO()):
        r = f()
    return time.perf_counter()-t, time.process_time()-p, r

n = 24
print("n=%d S=1/2 Heisenberg, maxm=30 nsweeps=6 (upper edge: 5 sweeps, maxdim 20)" % n)
rows = {"gs": [], "edge": [], "plusH": []}
for rep in range(3):
    s = session(n, seed=rep)
    w, cpu, e0 = timed(lambda: s.gs_energy())
    ps, tot, nc = per_sweep()
    rows["gs"].append((w, cpu))
    print("rep %d full solve   wall %.3f cpu %.3f E0=%.6f  matvecs/local solve per sweep %s (total %d over %d solves)"
          % (rep, w, cpu, e0, np.round(ps, 1).tolist(), tot, nc))
    w, cpu, em = timed(lambda: s._maximum_energy())
    ps, tot, nc = per_sweep()
    rows["edge"].append((w, cpu))
    print("rep %d upper edge   wall %.3f cpu %.3f Emax=%.6f matvecs/local solve per sweep %s (total %d over %d solves)"
          % (rep, w, cpu, em, np.round(ps, 1).tolist(), tot, nc))
    np.random.seed(100+rep)
    psi = randomMPS(s.sites, 30)
    w, cpu, ep = timed(lambda: D.dmrg(psi, s.H, s._make_sweeps(ns=5, maxdim=20)))
    ps, tot, nc = per_sweep()
    rows["plusH"].append((w, cpu))
    print("rep %d +H, edge sch wall %.3f cpu %.3f E=%.6f   matvecs/local solve per sweep %s (total %d over %d solves)"
          % (rep, w, cpu, ep, np.round(ps, 1).tolist(), tot, nc))
m = {k: np.min(np.array(v), axis=0) for k, v in rows.items()}
print("min of 3 wall: full %.3f, edge %.3f (ratio %.2f), +H on edge schedule %.3f"
      % (m["gs"][0], m["edge"][0], m["edge"][0]/m["gs"][0], m["plusH"][0]))
print("min of 3 cpu : full %.3f, edge %.3f (ratio %.2f), +H on edge schedule %.3f"
      % (m["gs"][1], m["edge"][1], m["edge"][1]/m["gs"][1], m["plusH"][1]))
r = np.array(rows["edge"])[:, 0]/np.array(rows["gs"])[:, 0]
print("per-rep wall ratios edge/full: %s" % np.round(r, 2).tolist())

for hz in (0.5, 1.0):
    s = session(n, hz=hz, seed=7)
    w, cpu, em = timed(lambda: s._maximum_energy())
    ps, tot, nc = per_sweep()
    print("hz=%.1f upper edge wall %.3f cpu %.3f Emax=%.6f (exact %.6f) matvecs/local solve per sweep %s (total %d)"
          % (hz, w, cpu, em, 0.25*(n-1)+hz*n/2, np.round(ps, 1).tolist(), tot))
```

`reviews/groundstate_6_new1/01_cost_and_iterations.out`:


```
n=24 S=1/2 Heisenberg, maxm=30 nsweeps=6 (upper edge: 5 sweeps, maxdim 20)
rep 0 full solve   wall 0.998 cpu 0.997 E0=-10.453786  matvecs/local solve per sweep [20.4, 11.7, 7.0, 4.2, 3.5, 3.5] (total 2315 over 276 solves)
rep 0 upper edge   wall 5.690 cpu 5.679 Emax=5.750000 matvecs/local solve per sweep [53.9, 92.8, 7.6, 2.0, 2.0] (total 7284 over 230 solves)
rep 0 +H, edge sch wall 0.842 cpu 0.841 E=-10.453784   matvecs/local solve per sweep [21.5, 6.0, 4.3, 4.3, 4.3] (total 1863 over 230 solves)
rep 1 full solve   wall 0.984 cpu 0.983 E0=-10.453786  matvecs/local solve per sweep [20.5, 11.8, 7.0, 4.2, 3.5, 3.5] (total 2324 over 276 solves)
rep 1 upper edge   wall 8.257 cpu 8.252 Emax=5.750000 matvecs/local solve per sweep [54.0, 99.2, 9.9, 2.0, 2.0] (total 7684 over 230 solves)
rep 1 +H, edge sch wall 1.108 cpu 1.108 E=-10.453784   matvecs/local solve per sweep [21.5, 6.8, 4.3, 4.3, 4.3] (total 1897 over 230 solves)
rep 2 full solve   wall 0.904 cpu 0.904 E0=-10.453786  matvecs/local solve per sweep [20.5, 11.2, 7.0, 4.2, 3.5, 3.5] (total 2296 over 276 solves)
rep 2 upper edge   wall 7.939 cpu 7.934 Emax=5.750000 matvecs/local solve per sweep [52.9, 109.4, 21.6, 2.1, 2.0] (total 8649 over 230 solves)
rep 2 +H, edge sch wall 1.028 cpu 1.028 E=-10.453784   matvecs/local solve per sweep [21.7, 6.0, 4.3, 4.3, 4.3] (total 1870 over 230 solves)
min of 3 wall: full 0.904, edge 5.690 (ratio 6.29), +H on edge schedule 0.842
min of 3 cpu : full 0.904, edge 5.679 (ratio 6.28), +H on edge schedule 0.841
per-rep wall ratios edge/full: [5.7, 8.39, 8.78]
hz=0.5 upper edge wall 0.268 cpu 0.268 Emax=11.750000 (exact 11.750000) matvecs/local solve per sweep [14.3, 2.1, 2.0, 2.0, 2.0] (total 1028)
hz=1.0 upper edge wall 0.242 cpu 0.242 Emax=17.750000 (exact 17.750000) matvecs/local solve per sweep [12.0, 2.0, 2.0, 2.0, 2.0] (total 923)
```

```bash
cd <scratch>/reviews/groundstate_6_new1 && <scratch>/run3.sh 02_fix_and_v3.py 2>&1 | tee 02_fix_and_v3.out
```

`reviews/groundstate_6_new1/02_fix_and_v3.py`:

```python
# Reviewer of groundstate_6_new1: (1) how the upper-edge matvecs are
# distributed (does any local solve hit the 200-iteration floor?);
# (2) does a looser Lanczos stopping rule for this BOUND solve alone keep
# Emax and remove the cost (tol 1e-6, and a 20-matvec cap); (3) v3's own
# edge/full ratio, min of 3; (4) what the first set_gs + gs_energy costs
# on a solved "python" chain after 867e2b4.
import numpy as np, io, contextlib, time
from dmrgpy import spinchain, groundstate
from dmrgpy.pyitensor import dmrg as D

MODE = {"tol": None, "niter": None}
LOG = []
_orig_lz = D._lanczos_ground_state
def lz(matvec, v0, niter=30, tol=1e-12, residual_tol=None):
    c = [0]
    def mv(x):
        c[0] += 1
        return matvec(x)
    if MODE["tol"] is not None: tol = MODE["tol"]
    if MODE["niter"] is not None: niter = MODE["niter"]
    out = _orig_lz(mv, v0, niter=niter, tol=tol, residual_tol=residual_tol)
    LOG.append(c[0])
    return out
D._lanczos_ground_state = lz

def ham(sc):
    h = 0
    for i in range(sc.ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return h
def chain(v, n, seed):
    np.random.seed(seed)
    b = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    b.maxm = 30; b.nsweeps = 6; b.set_hamiltonian(ham(b))
    return b
def session(v, n, seed):
    b = chain(v, n, seed)
    groundstate._session_parameters(b); groundstate.send_hamiltonian(b)
    return b._session
def T(f):
    t = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        r = f()
    return time.perf_counter()-t, r

n = 24
exact = 0.25*(n-1)
for label, tol, niter in (("default (tol 1e-12, floor 200)", None, None),
                          ("tol 1e-6", 1e-6, None),
                          ("cap 20 matvecs", None, 20)):
    ws = []; errs = []; hist = None
    for rep in range(3):
        s = session("python", n, rep)
        MODE["tol"] = None; MODE["niter"] = None
        T(lambda: s.gs_energy())
        MODE["tol"] = tol; MODE["niter"] = niter; LOG.clear()
        w, em = T(lambda: s._maximum_energy())
        MODE["tol"] = None; MODE["niter"] = None
        ws.append(w); errs.append(exact-em)
        if rep == 0:
            a = np.array(LOG)
            hist = "max %d, solves >=100: %d, >=199: %d of %d" % (a.max(), (a >= 100).sum(), (a >= 199).sum(), len(a))
    print("%-32s wall min %.3f s (per rep %s), exact-Emax per rep %s | rep0 matvecs %s"
          % (label, min(ws), np.round(ws, 3).tolist(), ["%.1e" % e for e in errs], hist))

# v3, the hunter's split: excited_states(1,1.0,False) twice, the second
# with the edges cached
r = []
for rep in range(3):
    s = session(3, n, 10+rep)
    t_gs, _ = T(lambda: s.gs_energy())
    t1, _ = T(lambda: s.excited_states(1, 1.0, False))
    t2, _ = T(lambda: s.excited_states(1, 1.0, False))
    r.append((t_gs, t1-t2))
r = np.array(r)
print("v3: full solve per rep %s, upper edge per rep %s, min ratio edge/full %.2f"
      % (np.round(r[:, 0], 3).tolist(), np.round(r[:, 1], 3).tolist(), r[:, 1].min()/r[:, 0].min()))

# first injection on a solved "python" chain that never ran KPM
b = chain("python", n, 20)
T(lambda: b.gs_energy())
x = b.random_mps()
t_amb, ex = T(lambda: (b.aMb(x, b.hamiltonian, x)/b.overlap(x, x)).real)
b.set_gs(x)
t1, e1 = T(lambda: b.gs_energy())
x2 = b.random_mps(); b.set_gs(x2)
t2, _ = T(lambda: b.gs_energy())
print("python solved chain: first set_gs+gs_energy %.3f s (e=%.6f, <x|H|x>=%.6f), second %.3f s, <x|H|x> alone %.4f s"
      % (t1, e1, ex, t2, t_amb))
```

`reviews/groundstate_6_new1/02_fix_and_v3.out`:


```
default (tol 1e-12, floor 200)   wall min 5.316 s (per rep [5.316, 7.037, 7.634]), exact-Emax per rep ['2.8e-11', '3.5e-11', '3.6e-11'] | rep0 matvecs max 193, solves >=100: 28, >=199: 0 of 230
tol 1e-6                         wall min 0.782 s (per rep [0.844, 0.825, 0.782]), exact-Emax per rep ['3.4e-05', '4.4e-05', '8.2e-05'] | rep0 matvecs max 47, solves >=100: 0, >=199: 0 of 230
cap 20 matvecs                   wall min 0.951 s (per rep [0.968, 0.951, 0.981]), exact-Emax per rep ['9.7e-11', '6.0e-11', '1.3e-10'] | rep0 matvecs max 20, solves >=100: 0, >=199: 0 of 230
v3: full solve per rep [0.216, 0.217, 0.351], upper edge per rep [0.166, 0.162, 0.334], min ratio edge/full 0.75
python solved chain: first set_gs+gs_energy 5.635 s (e=-0.011702, <x|H|x>=-0.011702), second 0.175 s, <x|H|x> alone 0.0299 s
```


**Suggested fix**: a small matvec cap, about 20, on the bound solve only, passed
to `_local_ground_state` so that it bypasses the `max(niter, 200)` floor at
`dmrg.py:421`; the floor stays for every other caller, since its own comment
records that the excited-state penalty method and spin-3/2 chains need it.
Measured, the cap gives 0.95 s at 1e-10 against 5.3 s. Loosening `tol` to 1e-6 is
as cheap (0.78 s) but pays a variational underestimate of 3e-5 to 8e-5, which the
comment at `chain.py:1820` permits but the cap does not need. Even capped, the
bound solve costs about one full solve, because `_make_sweeps` has no bond-dimension
ramp while `gs_energy` ramps from 10 to 30, and it starts at `randomMPS(self.maxm)`
before truncating to 20; v3 reaches 0.75 because ITensor's Davidson runs `niter=2`
per local solve. No number changes beyond Emax at the 1e-10 level.

### 11. `867e2b4` put `ground_state_on_session` at the top of `get_dynamical_correlator`, ahead of every KPM argument check, so on v2, v3 and `"python"` a malformed KPM call runs a full ground-state solve before it raises (where `30200a4` and `mode="ED"` raise with none), and the first `submode="SECTOR"` call on an unsolved chain pays one global solve that SECTOR never reads

`bug` &middot; severity **LOW** &middot; CONFIRMED (the SECTOR half CONFIRMED, NARROWED) &middot; lens `kpm` (the SECTOR half turned up by this finding's reviewer)

**Where**: `src/dmrgpy/dynamics.py:199-202` (`ground_state_on_session` ahead of the
dispatch), against the checks in `kpmdmrg.py:109-151` (unknown-keyword
`TypeError`, `delta<0`, `validate_kpm_n_scale`, the v2 `kpm_energy_truncate`
refusal) and SECTOR's own solves on a clone (`sectordc.py:286-292`, the clone
made at `:314`); the texts that are now false: `kpmdmrg.py:124` ("validated here,
once, ahead of both session calls and of the ground state"),
`algebra/kpm.py::validate_kpm_n_scale`'s docstring ("before any ground-state
work") and the second pass's ruled-out line "rejected before the ground state on
every Hermitian route".

`30200a4`'s first line in the dispatcher was `set_initial_wf(self.wf0)`, which only
reset flags, and the KPM route validated its arguments before its own
`get_gs()`. `867e2b4` replaced that line by `ground_state_on_session`, which
solves, so the solve now happens before any submode has looked at its arguments.
Two things follow. A malformed KPM call solves and then raises. And SECTOR, which
solves both of its sectors, and measures its reference charge, on an internal
clone and never reads the caller's state, now also solves the caller's chain on
its first call. The same shape reaches SECTOR on v2 (a solve, then SECTOR's own
`NotImplementedError`) and the non-Hermitian KPM's missing-`E_max` `ValueError`.
It survived because the only before-the-ground-state spy in the suite
(`GroundStateSpy` in `tests/test_audit_2026_09_24b_kpm.py`) patches ED's
`get_gs_array`, and the DMRG tests check only that the error is raised.

**Expected**: a malformed call rejected before any ground-state work, as on
`mode="ED"` and as at `30200a4`; and a SECTOR call that makes no solve on the
caller's chain, as at `30200a4`.

Repro, from the hunter (24-site Heisenberg chain, `maxm=40`, `nsweeps=6`):

```bash
cd <scratch>/kpm && <scratch>/run3.sh 06_dmrg_check_order.py 2>&1 | tee 06_dmrg_check_order.out
```

`kpm/06_dmrg_check_order.py`:

```python
# On the DMRG routes, is kpm_n_scale (and the other KPM argument checks:
# unknown keyword, delta<0, kpm_energy_truncate on v2) still rejected before
# any ground-state work, as the previous record's ruled-out list says?  867e2b4
# put groundstate.ground_state_on_session at the top of
# dynamics.get_dynamical_correlator, ahead of kpmdmrg's checks.
import time
from dmrgpy import spinchain
from dmrgpy import groundstate

calls = []
_orig = groundstate.gs_energy
def spy(*a, **k):
    calls.append(1)
    return _orig(*a, **k)
groundstate.gs_energy = spy

L = 24
def fresh(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40; sc.nsweeps = 6
    return sc

cases = [("kpm_n_scale=1.5", dict(attr=("kpm_n_scale", 1.5), kw={})),
         ("unknown keyword deltaa=", dict(attr=None, kw={"deltaa": 0.1})),
         ("delta=-0.1", dict(attr=None, kw={"delta": -0.1})),
         ("kpm_energy_truncate on v2", dict(attr=("kpm_energy_truncate", True), kw={}))]
for v in ("python", 3, 2):
    for label, c in cases:
        if label.startswith("kpm_energy_truncate") and v != 2: continue
        sc = fresh(v)
        if c["attr"]: setattr(sc, *c["attr"])
        del calls[:]
        t0 = time.time()
        try:
            sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), **c["kw"])
            res = "returned"
        except Exception as e:
            res = "%s: %s" % (type(e).__name__, str(e)[:60])
        print("v=%-6s %-26s -> %s | gs_energy calls=%d computed_gs=%s e0=%s %.2fs" % (
            v, label, res, len(calls), sc.computed_gs, sc.e0, time.time()-t0))
# the ED route, for contrast (4 sites)
sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
sc.set_hamiltonian(sum(sc.Sz[i]*sc.Sz[i+1]+sc.Sx[i]*sc.Sx[i+1] for i in range(3)))
sc.kpm_n_scale = 1.5
ed = sc.get_ED_obj()
try:
    sc.get_dynamical_correlator(mode="ED", name=(sc.Sz[0], sc.Sz[0]))
except Exception as e:
    print("ED kpm_n_scale=1.5 ->", type(e).__name__, "| ED computed_gs =", ed.computed_gs)
```

`kpm/06_dmrg_check_order.out`:


```
v=python kpm_n_scale=1.5            -> TypeError: kpm_n_scale must be a positive integer (it multiplies the ca | gs_energy calls=1 computed_gs=True e0=-10.453785758037014 1.86s
v=python unknown keyword deltaa=    -> TypeError: Unexpected keyword argument(s) for the KPM dynamical correla | gs_energy calls=1 computed_gs=True e0=-10.453785758036968 1.31s
v=python delta=-0.1                 -> ValueError: delta must be >= 0, got -0.1 | gs_energy calls=1 computed_gs=True e0=-10.453785758037114 1.21s
v=3      kpm_n_scale=1.5            -> TypeError: kpm_n_scale must be a positive integer (it multiplies the ca | gs_energy calls=1 computed_gs=True e0=-10.453785758033138 0.80s
v=3      unknown keyword deltaa=    -> TypeError: Unexpected keyword argument(s) for the KPM dynamical correla | gs_energy calls=1 computed_gs=True e0=-10.453785758031206 0.40s
v=3      delta=-0.1                 -> ValueError: delta must be >= 0, got -0.1 | gs_energy calls=1 computed_gs=True e0=-10.453785758031024 0.37s
v=2      kpm_n_scale=1.5            -> TypeError: kpm_n_scale must be a positive integer (it multiplies the ca | gs_energy calls=1 computed_gs=True e0=-10.453785758032371 0.88s
v=2      unknown keyword deltaa=    -> TypeError: Unexpected keyword argument(s) for the KPM dynamical correla | gs_energy calls=1 computed_gs=True e0=-10.45378575803365 0.91s
v=2      delta=-0.1                 -> ValueError: delta must be >= 0, got -0.1 | gs_energy calls=1 computed_gs=True e0=-10.453785758034087 0.92s
v=2      kpm_energy_truncate on v2  -> NotImplementedError: kpm_energy_truncate is implemented for itensor_version=3 and | gs_energy calls=1 computed_gs=True e0=-10.453785758034112 0.93s
ED kpm_n_scale=1.5 -> TypeError | ED computed_gs = False
```


**Reviewer (CONFIRMED)**: rerun at 16 sites, `maxm=30`, `nsweeps=6`: on every
backend each malformed call raises after exactly one `gs_energy` call with
`computed_gs` True afterwards (`"python"` 0.42 to 0.56 s, v3 0.14 s, v2 0.25 to
0.29 s, v2 including `kpm_energy_truncate`). On a read-only archive of `30200a4`
all three raise with zero `gs_energy` calls in 0.00 to 0.02 s, so the ordering was
true there and `867e2b4` broke it; the second pass's ruled-out item re-enters for
exactly that reason. Softening, not a strike: a corrected retry on the same chain
makes zero further calls, so interactively the solve is brought forward rather
than wasted, and is lost only by a script that dies on the raise or moves on. TD
was never ahead of its solve (`timedependent.dynamical_correlator` calls
`get_gs()` before reading its keywords at both commits), so the shape is
KPM-specific among the Hermitian submodes. A correct non-Hermitian KPM call costs
one caller solve at both commits and returns identical values: the solve is moved,
not duplicated. Nothing struck.

```bash
cd <scratch>/reviews/kpm_2 && <scratch>/run3.sh 01_check_order_current.py 2>&1 | tee 01_check_order_current.out
```

`reviews/kpm_2/01_check_order_current.py`:

```python
# Reviewer probe for kpm_2, on the tree under audit (867e2b4).
# (a) the hunter's four malformed KPM calls, counting gs_energy calls and
#     session sweeps before the raise, on "python", v3 and v2;
# (b) the retry pattern: after the raise, fix the argument and call again,
#     counting gs_energy calls and whether the session sweeps again;
# (c) a correct call from a fresh chain, for the reference cost.
import sys, time
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain, groundstate

calls = []
_orig = groundstate.gs_energy
def spy(*a, **k):
    calls.append(1)
    return _orig(*a, **k)
groundstate.gs_energy = spy

L = 16
def fresh(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 30; sc.nsweeps = 6; sc.kpmmaxm = 20
    return sc

es = [0.0, 0.5, 1.0]
cases = [("kpm_n_scale=1.5", ("kpm_n_scale", 1.5), {}, ("kpm_n_scale", 1)),
         ("unknown kw deltaa=", None, {"deltaa": 0.1}, None),
         ("delta=-0.1", None, {"delta": -0.1}, None),
         ("kpm_energy_truncate v2", ("kpm_energy_truncate", True), {}, ("kpm_energy_truncate", False))]
for v in ("python", 3, 2):
    for label, attr, kw, fix in cases:
        if label.startswith("kpm_energy_truncate") and v != 2: continue
        sc = fresh(v)
        if attr: setattr(sc, *attr)
        del calls[:]
        t0 = time.time()
        try:
            sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=es, **kw)
            res = "returned"
        except Exception as e:
            res = type(e).__name__
        t1 = time.time() - t0
        n1 = len(calls); cg = sc.computed_gs
        # retry with the argument corrected
        if fix: setattr(sc, *fix)
        del calls[:]
        t0 = time.time()
        sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=es, delta=0.2)
        t2 = time.time() - t0
        print("v=%-6s %-24s -> %-19s gs_energy=%d computed_gs=%s %.2fs | retry: gs_energy=%d %.2fs" % (
            v, label, res, n1, cg, t1, len(calls), t2))
    # reference: a correct call on a fresh chain, and a correct call after gs_energy()
    sc = fresh(v); del calls[:]; t0 = time.time()
    sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=es, delta=0.2)
    print("v=%-6s fresh correct call: gs_energy=%d %.2fs" % (v, len(calls), time.time()-t0))
    sys.stdout.flush()
```

`reviews/kpm_2/01_check_order_current.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
v=python kpm_n_scale=1.5          -> TypeError           gs_energy=1 computed_gs=True 0.56s | retry: gs_energy=0 7.50s
v=python unknown kw deltaa=       -> TypeError           gs_energy=1 computed_gs=True 0.50s | retry: gs_energy=0 7.57s
v=python delta=-0.1               -> ValueError          gs_energy=1 computed_gs=True 0.42s | retry: gs_energy=0 5.65s
v=python fresh correct call: gs_energy=1 7.38s
v=3      kpm_n_scale=1.5          -> TypeError           gs_energy=1 computed_gs=True 0.14s | retry: gs_energy=0 0.92s
v=3      unknown kw deltaa=       -> TypeError           gs_energy=1 computed_gs=True 0.14s | retry: gs_energy=0 0.90s
v=3      delta=-0.1               -> ValueError          gs_energy=1 computed_gs=True 0.14s | retry: gs_energy=0 0.89s
v=3      fresh correct call: gs_energy=1 1.08s
v=2      kpm_n_scale=1.5          -> TypeError           gs_energy=1 computed_gs=True 0.29s | retry: gs_energy=0 1.53s
v=2      unknown kw deltaa=       -> TypeError           gs_energy=1 computed_gs=True 0.26s | retry: gs_energy=0 1.51s
v=2      delta=-0.1               -> ValueError          gs_energy=1 computed_gs=True 0.25s | retry: gs_energy=0 1.55s
v=2      kpm_energy_truncate v2   -> NotImplementedError gs_energy=1 computed_gs=True 0.25s | retry: gs_energy=0 1.53s
v=2      fresh correct call: gs_energy=1 1.80s
```

```bash
cd <scratch>/reviews/kpm_2 && <scratch>/run3.sh 02_check_order_30200a4.py 2>&1 | tee 02_check_order_30200a4.out
```

`reviews/kpm_2/02_check_order_30200a4.py`:

```python
# Reviewer probe for kpm_2: the same malformed KPM calls on the tree as it
# was at 30200a4 (the parent of 867e2b4's fixes, extracted read-only with
# git archive into ./old30200a4), itensor_version="python" only, since the
# extracted tree has no compiled extension. The anchor for "before 867e2b4
# it raised with no solve".
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "old30200a4", "src"))
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain, groundstate

calls = []
_orig = groundstate.gs_energy
def spy(*a, **k):
    calls.append(1)
    return _orig(*a, **k)
groundstate.gs_energy = spy

L = 16
def fresh(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 30; sc.nsweeps = 6; sc.kpmmaxm = 20
    return sc

es = [0.0, 0.5, 1.0]
cases = [("kpm_n_scale=1.5", ("kpm_n_scale", 1.5), {}),
         ("unknown kw deltaa=", None, {"deltaa": 0.1}),
         ("delta=-0.1", None, {"delta": -0.1})]
v = "python"
for label, attr, kw in cases:
    sc = fresh(v)
    if attr: setattr(sc, *attr)
    del calls[:]
    t0 = time.time()
    try:
        sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), es=es, **kw)
        res = "returned"
    except Exception as e:
        res = type(e).__name__
    print("30200a4 v=%-6s %-24s -> %-19s gs_energy=%d computed_gs=%s %.2fs" % (
        v, label, res, len(calls), sc.computed_gs, time.time()-t0))
```

`reviews/kpm_2/02_check_order_30200a4.out`:


```
dmrgpy from <scratch>/reviews/kpm_2/old30200a4/src/dmrgpy/__init__.py
30200a4 v=python kpm_n_scale=1.5          -> TypeError           gs_energy=0 computed_gs=False 0.02s
30200a4 v=python unknown kw deltaa=       -> TypeError           gs_energy=0 computed_gs=False 0.00s
30200a4 v=python delta=-0.1               -> ValueError          gs_energy=0 computed_gs=False 0.00s
```

```bash
cd <scratch>/reviews/kpm_2 && <scratch>/run3.sh 03_sector_nh_extra_solve.py 2>&1 | tee 03_sector_nh_extra_solve.out
```

`reviews/kpm_2/03_sector_nh_extra_solve.py`:

```python
# Reviewer probe for kpm_2, tree under audit (867e2b4). Does the eager
# ground_state_on_session at the top of dynamics.get_dynamical_correlator
# also cost a solve on CORRECT calls of submodes that never read the
# caller's ground state?  submode="SECTOR" solves on an internal clone;
# the non-Hermitian KPM calls self.gs_energy() itself.
# gs_energy calls are counted per chain object: the caller's chain versus
# any other chain (the SECTOR clones).
import sys, os, time
if len(sys.argv) > 1 and sys.argv[1] == "old":
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "old30200a4", "src"))
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain, groundstate
from dmrgpy.operatornames import name2MO

calls = []
_orig = groundstate.gs_energy
def spy(self, *a, **k):
    calls.append(id(self))
    return _orig(self, *a, **k)
groundstate.gs_energy = spy

def heis(n, v, nh=0.0):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    if nh: h = h + 1j*nh*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm = 30; sc.nsweeps = 6; sc.kpmmaxm = 20
    return sc

versions = ("python",) if (len(sys.argv) > 1 and sys.argv[1] == "old") else ("python", 3, 2)
es = np.linspace(-0.5, 3.0, 50)
for v in versions:
    sc = heis(8, v)
    Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)
    del calls[:]
    t0 = time.time()
    try:
        x, y = sc.get_dynamical_correlator(submode="SECTOR", name=(Sm[0], Sp[0]),
                                           nex=8, es=es, delta=0.2, quiet=True)
        res = "returned, max y=%.6f" % np.max(np.real(y))
    except Exception as e:
        res = "%s: %s" % (type(e).__name__, str(e)[:50])
    mine = sum(1 for c in calls if c == id(sc))
    print("v=%-6s SECTOR 8 sites -> %s | gs_energy on caller=%d on others=%d caller computed_gs=%s %.2fs" % (
        v, res, mine, len(calls)-mine, sc.computed_gs, time.time()-t0))
    sys.stdout.flush()
for v in versions:
    if v == 2: continue
    sc = heis(6, v, nh=0.2)
    del calls[:]
    t0 = time.time()
    try:
        out = sc.get_dynamical_correlator(submode="KPM", name=(sc.Sz[1], sc.Sz[1]),
                                          es=np.linspace(-1, 3, 20), delta=0.3)
        res = "returned"
    except Exception as e:
        res = "%s: %s" % (type(e).__name__, str(e)[:50])
    mine = sum(1 for c in calls if c == id(sc))
    print("v=%-6s NH-KPM 6 sites -> %s | gs_energy on caller=%d on others=%d %.2fs" % (
        v, res, mine, len(calls)-mine, time.time()-t0))
    sys.stdout.flush()
```

`reviews/kpm_2/03_sector_nh_extra_solve.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
v=python SECTOR 8 sites -> returned, max y=0.360100 | gs_energy on caller=1 on others=3 caller computed_gs=True 3.62s
v=3      SECTOR 8 sites -> returned, max y=0.360301 | gs_energy on caller=1 on others=3 caller computed_gs=True 1.91s
v=2      SECTOR 8 sites -> NotImplementedError: get_dynamical_correlator(submode="SECTOR") needs c | gs_energy on caller=1 on others=0 caller computed_gs=True 0.03s
v=python NH-KPM 6 sites -> ValueError: E_max (an upper bound for the spectral radius of t | gs_energy on caller=1 on others=0 0.18s
v=3      NH-KPM 6 sites -> ValueError: E_max (an upper bound for the spectral radius of t | gs_energy on caller=1 on others=0 0.14s
```

```bash
cd <scratch>/reviews/kpm_2 && <scratch>/run3.sh 04_sector_nh_extra_solve_30200a4.py 2>&1 | tee 04_sector_nh_extra_solve_30200a4.out
```

`reviews/kpm_2/04_sector_nh_extra_solve_30200a4.py`:

```python
# Reviewer probe for kpm_2: 03_sector_nh_extra_solve.py run against the
# tree as it was at 30200a4 (./old30200a4, extracted read-only with git
# archive), itensor_version="python" only.
import sys, os, runpy
here = os.path.dirname(os.path.abspath(__file__))
sys.argv = [os.path.join(here, "03_sector_nh_extra_solve.py"), "old"]
runpy.run_path(sys.argv[0], run_name="__main__")
```

`reviews/kpm_2/04_sector_nh_extra_solve_30200a4.out`:


```
dmrgpy from <scratch>/reviews/kpm_2/old30200a4/src/dmrgpy/__init__.py
v=python SECTOR 8 sites -> returned, max y=0.360100 | gs_energy on caller=0 on others=3 caller computed_gs=False 5.26s
v=python NH-KPM 6 sites -> ValueError: E_max (an upper bound for the spectral radius of t | gs_energy on caller=0 on others=0 0.03s
```

```bash
cd <scratch>/reviews/kpm_2 && <scratch>/run3.sh 05_nhkpm_correct_call.py 2>&1 | tee 05_nhkpm_correct_call.out
```

`reviews/kpm_2/05_nhkpm_correct_call.py`:

```python
# Reviewer probe for kpm_2: a CORRECT non-Hermitian KPM call (E_max given)
# on "python" and v3, counting gs_energy calls on the caller's chain and the
# returned curve, to tell a moved solve (1 call) from an extra one (2).
# With argument "old", runs against ./old30200a4 ("python" only).
import sys, os, time
OLD = len(sys.argv) > 1 and sys.argv[1] == "old"
if OLD:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "old30200a4", "src"))
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain, groundstate

calls = []
_orig = groundstate.gs_energy
def spy(self, *a, **k):
    calls.append(id(self))
    return _orig(self, *a, **k)
groundstate.gs_energy = spy

def heis(n, v, nh):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 1j*nh*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm = 30; sc.nsweeps = 6; sc.kpmmaxm = 20
    return sc

for v in (("python",) if OLD else ("python", 3)):
    sc = heis(6, v, 0.2)
    del calls[:]
    t0 = time.time()
    x, y = sc.get_dynamical_correlator(submode="KPM", name=(sc.Sz[1], sc.Sz[1]),
                                       es=np.linspace(-1, 3, 6), delta=0.3,
                                       E_max=6.0, n=60)
    mine = sum(1 for c in calls if c == id(sc))
    print("v=%-6s NH-KPM 6 sites correct call | gs_energy on caller=%d on others=%d %.2fs | y[:3]=%s" % (
        v, mine, len(calls)-mine, time.time()-t0, np.round(y[:3], 6)))
    sys.stdout.flush()
```

`reviews/kpm_2/05_nhkpm_correct_call.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
v=python NH-KPM 6 sites correct call | gs_energy on caller=1 on others=0 11.96s | y[:3]=[ 3.85000e-04+0.001112j -2.25530e-02+0.15554j  -5.76246e-01+0.126186j]
v=3      NH-KPM 6 sites correct call | gs_energy on caller=1 on others=0 4.50s | y[:3]=[ 3.85000e-04+0.001112j -2.25530e-02+0.15554j  -5.76246e-01+0.126186j]
```

```bash
cd <scratch>/reviews/kpm_2 && <scratch>/run3.sh 06_nhkpm_correct_call_30200a4.py 2>&1 | tee 06_nhkpm_correct_call_30200a4.out
```

`reviews/kpm_2/06_nhkpm_correct_call_30200a4.py`:

```python
# Reviewer probe for kpm_2: 05_nhkpm_correct_call.py against ./old30200a4.
import sys, os, runpy
here = os.path.dirname(os.path.abspath(__file__))
sys.argv = [os.path.join(here, "05_nhkpm_correct_call.py"), "old"]
runpy.run_path(sys.argv[0], run_name="__main__")
```

`reviews/kpm_2/06_nhkpm_correct_call_30200a4.out`:


```
dmrgpy from <scratch>/reviews/kpm_2/old30200a4/src/dmrgpy/__init__.py
v=python NH-KPM 6 sites correct call | gs_energy on caller=1 on others=0 13.37s | y[:3]=[ 3.85000e-04+0.001112j -2.25530e-02+0.15554j  -5.76246e-01+0.126186j]
```


**Reviewer of the SECTOR half (CONFIRMED, NARROWED)**: reproduced on an 8-site
Heisenberg chain, `name=(Sm[0],Sp[0])`, `nex=8`, `delta=0.2`, counting
`gs_energy` per chain object: on a fresh chain the caller makes 1 solve and the
clones 3, max y = 0.360100 (`"python"`, 3.08 s in total, 0.11 s of it in
`ground_state_on_session`) and 0.360301 (v3, 1.84 s, 0.03 s); at `30200a4` the
caller makes 0 and the clones 3, with the same curve. A caller that solved first
pays nothing, a repeated call pays nothing (the `_sector_states_cache`), a
different site pays nothing. By reading, `sectordc.py` reads nothing of the
caller's state (`_measure_charge` measures `clone.vev`), and `get_spectral_function`
/`get_spin_spectral_function` call `sectordc.sector_poles` directly and never reach
the dispatcher. At size (`maxm=40`, `nsweeps=8`, `nex=6`) the extra solve is 0.46 s
of 80.29 s on `"python"` at n=14 (0.6 per cent) and 0.44 s of 12.15 s on v3 at n=24
(3.6 per cent). Struck, each explicitly:

- "Every correct SECTOR call": only the first call on a chain whose ground state is
  not current pays.
- "`867e2b4` also makes clones reset their ground state, so the caller's solve is
  not consumed" as a causal link: the solve is waste because SECTOR's pipeline never
  reads the caller's state at all.
- The v2 row (a `NotImplementedError` after a caller solve): a solve ahead of a
  refusal, which belongs to the first half of this entry.
- Any implication that `get_spectral_function`/`get_spin_spectral_function` pay
  it: by reading they never reach the dispatcher.

```bash
cd <scratch>/reviews/kpm_2_new1 && <scratch>/run3.sh 01_repro_and_attacks.py 2>&1 | tee 01_repro_and_attacks.out
```

`reviews/kpm_2_new1/01_repro_and_attacks.py`:

```python
# Reviewer of kpm_2_new1. (a) The hunter's repro (8-site Heisenberg,
# SECTOR, name=(Sm0,Sp0), nex=8), python and v3. (b) The same call on a
# chain already solved by the caller. (c) Two consecutive SECTOR calls on
# an unsolved chain. (d) A second, different SECTOR call (another site)
# after the first. gs_energy counted per chain object, and the time spent
# inside ground_state_on_session measured separately.
import sys, time
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain, groundstate
from dmrgpy.operatornames import name2MO

calls = []
_orig = groundstate.gs_energy
def spy(self, *a, **k):
    calls.append(id(self))
    return _orig(self, *a, **k)
groundstate.gs_energy = spy
tgos = [0.0]
_orig_gos = groundstate.ground_state_on_session
def spy_gos(self):
    t = time.time(); r = _orig_gos(self); tgos[0] += time.time() - t; return r
groundstate.ground_state_on_session = spy_gos

def heis(n, v):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 30; sc.nsweeps = 6; sc.kpmmaxm = 20
    return sc

es = np.linspace(-0.5, 3.0, 50)
def sector_call(sc, site=0):
    Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)
    del calls[:]; tgos[0] = 0.0
    t0 = time.time()
    x, y = sc.get_dynamical_correlator(submode="SECTOR", name=(Sm[site], Sp[site]),
                                       nex=8, es=es, delta=0.2, quiet=True)
    mine = sum(1 for c in calls if c == id(sc))
    return (np.max(np.real(y)), mine, len(calls)-mine, time.time()-t0, tgos[0])

fmt = "  max y=%.6f | gs_energy caller=%d others=%d | total %.2fs, in ground_state_on_session %.2fs"
for v in ("python", 3):
    print("=== itensor_version=%r" % (v,))
    sc = heis(8, v)
    print("(a) unsolved caller, first SECTOR call:"); print(fmt % sector_call(sc))
    print("(c) same chain, same call again:"); print(fmt % sector_call(sc))
    print("(d) same chain, site 3:"); print(fmt % sector_call(sc, 3))
    sc = heis(8, v)
    e = sc.gs_energy()
    print("(b) caller solved first (e0=%.8f), then SECTOR:" % e); print(fmt % sector_call(sc))
    sys.stdout.flush()
```

`reviews/kpm_2_new1/01_repro_and_attacks.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
=== itensor_version='python'
(a) unsolved caller, first SECTOR call:
  max y=0.360100 | gs_energy caller=1 others=3 | total 3.08s, in ground_state_on_session 0.11s
(c) same chain, same call again:
  max y=0.360100 | gs_energy caller=0 others=1 | total 0.15s, in ground_state_on_session 0.00s
(d) same chain, site 3:
  max y=0.265237 | gs_energy caller=0 others=1 | total 0.15s, in ground_state_on_session 0.00s
(b) caller solved first (e0=-3.37493260), then SECTOR:
  max y=0.360100 | gs_energy caller=0 others=3 | total 2.96s, in ground_state_on_session 0.00s
=== itensor_version=3
(a) unsolved caller, first SECTOR call:
  max y=0.360301 | gs_energy caller=1 others=3 | total 1.84s, in ground_state_on_session 0.03s
(c) same chain, same call again:
  max y=0.360301 | gs_energy caller=0 others=1 | total 0.05s, in ground_state_on_session 0.00s
(d) same chain, site 3:
  max y=0.268602 | gs_energy caller=0 others=1 | total 0.05s, in ground_state_on_session 0.00s
(b) caller solved first (e0=-3.37493260), then SECTOR:
  max y=0.360301 | gs_energy caller=0 others=3 | total 1.81s, in ground_state_on_session 0.00s
```

```bash
cd <scratch>/reviews/kpm_2_new1 && <scratch>/run3.sh 02_share_at_size.py 2>&1 | tee 02_share_at_size.out
```

`reviews/kpm_2_new1/02_share_at_size.py`:

```python
# Reviewer of kpm_2_new1: the wall-time share of the caller's extra
# global solve in a first SECTOR call on an unsolved chain, at a size
# larger than the hunter's 8 sites. Heisenberg, name=(Sm0,Sp0), nex=6.
import sys, time
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain, groundstate
from dmrgpy.operatornames import name2MO

calls = []
_orig = groundstate.gs_energy
def spy(self, *a, **k):
    calls.append(id(self))
    return _orig(self, *a, **k)
groundstate.gs_energy = spy
tgos = [0.0]
_orig_gos = groundstate.ground_state_on_session
def spy_gos(self):
    t = time.time(); r = _orig_gos(self); tgos[0] += time.time() - t; return r
groundstate.ground_state_on_session = spy_gos

def heis(n, v):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40; sc.nsweeps = 8
    return sc

es = np.linspace(-0.5, 3.0, 50)
for v, n in (("python", 14), (3, 24)):
    sc = heis(n, v)
    Sp, Sm = name2MO("Sp", sc), name2MO("Sm", sc)
    del calls[:]; tgos[0] = 0.0
    t0 = time.time()
    x, y = sc.get_dynamical_correlator(submode="SECTOR", name=(Sm[0], Sp[0]),
                                       nex=6, es=es, delta=0.2, quiet=True)
    tt = time.time() - t0
    mine = sum(1 for c in calls if c == id(sc))
    print("v=%-6s n=%d: max y=%.6f | gs_energy caller=%d others=%d | total %.2fs, "
          "caller solve %.2fs = %.1f%%" % (v, n, np.max(np.real(y)), mine,
          len(calls)-mine, tt, tgos[0], 100*tgos[0]/tt))
    sys.stdout.flush()
```

`reviews/kpm_2_new1/02_share_at_size.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
v=python n=14: max y=0.286387 | gs_energy caller=1 others=3 | total 80.29s, caller solve 0.46s = 0.6%
v=3      n=24: max y=0.198101 | gs_energy caller=1 others=3 | total 12.15s, caller solve 0.44s = 3.6%
```


**Suggested fix**: take `ground_state_on_session` off the top of the dispatcher
and call it at each session-driving submode's own `self.get_gs()` site, after that
submode's argument checks: for KPM, replace `self.get_gs()` in
`kpmdmrg.dynamical_correlator_moments` by `ground_state_on_session(self)`, and let
SECTOR and the non-Hermitian branch not call it at all. It must replace those
`get_gs()` calls rather than simply be deleted, since what it adds over `get_gs()`
(`hamiltonian_on_session`/`send_hamiltonian` for `restart=False`, and the
injected-state hand-off of the second pass's findings 11 and 12) has to stay on
every submode that drives the session. Done inside `kpmdmrg`, it also reaches the
direct callers `get_dynamical_correlator_moments` and `fermionchain.get_gr`; it
does not by itself fix finding 5, whose cause is `gs_is_current`, so the two land
together. A KPM-only pre-dispatch validator, the hunter's first option, would
duplicate four checks already in `kpmdmrg` and leave SECTOR and the non-Hermitian
route as they are. The narrowest SECTOR-only edit is `if submode != "SECTOR":
ground_state_on_session(self)`. A regression spy on `groundstate.gs_energy` for
the DMRG routes is missing and should come with it. No number changes; one side
effect goes away, that the caller's chain comes out of a SECTOR call with
`computed_gs=True`.

### 12. The DMRG Hermiticity probe `Many_Body_Chain.is_hermitian` thresholds the unnormalized ||(A-A^dag)w||^2 at an absolute 1e-4, so any anti-Hermitian part below about 1e-2 in absolute size is called Hermitian: `gs_energy()` then drops a weak decay rate whole on v2 and v3 (Im E0 = -0.002134 returned as 0) and returns a real part 3 to 18 per cent off on `"python"`, and since `867e2b4` `disentangle_manifold` diagonalizes the operator's Hermitian part instead of the operator (eigen-residual 0.702 of max|ma| against 1e-17)

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lenses `kpm` and `misc` (found by both independently)

**Where**: `src/dmrgpy/mpsalgebra.py:362` (`return not norm>1e-4`, `40e526e`,
2022) and the comment at `:355-357` saying `applyoperator()` normalizes its
result, which is false; every gate that reads the verdict:
`mpsalgebratk/disentangle.py:15`/`:26` (new in `867e2b4`), `groundstate.py:195`,
`:408` and `:526`, `mpsalgebra.py:52`/`:75`, `dynamics.py:203` and `:286`,
`sectordc.py:118`, `manybodychain.py:786` (the verdict's cache); the ED side
decides on two other absolute tests, `algebra/algebra.py:428-436` (1e-8 on
||h-h^dag||_F^2, read at `edtk/dynamics.py:29`) and `:356-360` (1e-6 max-abs,
read in `lowest_states`).

The chain's probe runs the canonical proof first and, when it does not land, applies
A - A^dag to a unit random witness and compares ||(A-A^dag)w||^2 against 1e-4.
The witness is unit-normalized and the result is not, so the threshold is 1e-4 in
squared energy units of whatever units the operator is written in, and the
Hermitian/non-Hermitian dispatch becomes a property of the units rather than of
the operator. Two regimes fall through: a whole operator written small (a
non-Hermitian Hamiltonian rescaled by s <= 1e-2), and an O(1) Hamiltonian with a
weak anti-Hermitian part, which is exactly weak dissipation, where the decay rate
is the answer. The verdict routes `gs_energy()` to Hermitian DMRG, the correlator
dispatch to the Hermitian KPM, and, since `867e2b4`, `disentangle_manifold` to
`eigh` of the symmetrized representation, where before the commit a non-Hermitian
A went to the bare proof, got False and was diagonalized correctly by `eig`. Near
the threshold the verdict is witness-dependent. It survived because every
non-Hermitian operator in `tests/` and in the earlier records is O(1) in its own
units, where the absolute and relative criteria agree, and the one pinned
non-Hermitian `disentangle` case, `Sx0 + 1j*Sy0 + 0.3*Sz0 + 0.1*Sz1`, has a probe
norm of order 1.

**Expected**: a verdict that does not depend on units. Anchors: E0(s*H) = s*E0(H),
with the exact E0 = -1.242697 +- 0.291419i for H = Heis4 + 0.3*Sz_tot + 1j*Sz0
(a conjugate pair from `np.linalg.eigvals`), which ED returns at every s and every
DMRG backend at s=1; the open Hatano-Nelson chain (n=6, g=0.6), similar to the XX
chain with hopping sqrt(1-g^2), so its E0/s is exactly -1.3975836830; and
`disentangle_manifold` on the full 8-dimensional space of 3 spins returning right
eigenvectors of A = eps*(Sz0 + S+0) + 0.3*eps*Sz1, an eigen-residual at roundoff.

Repro, from the `kpm` hunter (`03` the probe's verdicts against ED's, `04`/`05`
rescaled models against ED, `17` three fresh chains per point on every backend):

```bash
cd <scratch>/kpm && <scratch>/run3.sh 03_herm_verdicts.py 2>&1 | tee 03_herm_verdicts.out
```

`kpm/03_herm_verdicts.py`:

```python
# Which Hermiticity verdict does each mode take for H = Heisenberg + c*1j*Sz0
# (4 sites)?  mode="ED" decides on edtk/dynamics.py's is_hermitian(h), an
# absolute 1e-8 on ||h-h^dag||_F^2; mode="DMRG" on Many_Body_Chain.is_hermitian
# (symbolic proof, then a random-witness probe with a 1e-4 threshold).  Also:
# does the probe's applyoperator normalize (the docstring says so)?
import numpy as np
from dmrgpy import spinchain
from dmrgpy.algebra import algebra

L = 4
for v in ("python", 3):
    print("itensor_version", v)
    for c in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 1e-4, 3e-5, 1e-5, 1e-6):
        sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
        h = 0
        for i in range(L-1):
            h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        h = h + c*1j*sc.Sz[0]
        sc.set_hamiltonian(h)
        dmrg_verdict = sc.is_hermitian(sc.hamiltonian)
        hm = sc.get_ED_obj().get_hamiltonian()
        ed_verdict = algebra.is_hermitian(hm)
        fro2 = np.sum(np.abs((hm - hm.conj().T).toarray())**2) if hasattr(hm, "toarray") else np.sum(np.abs(hm-hm.conj().T)**2)
        print("  c=%.0e  ED is_hermitian=%s (||h-h^dag||_F^2=%.2e)  DMRG is_hermitian=%s" % (
            c, ed_verdict, fro2, dmrg_verdict))
    # probe normalization
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    sc.set_hamiltonian(sum(sc.Sz[i]*sc.Sz[i+1] for i in range(L-1)))
    wf = sc.random_mps()
    for c in (1.0, 1e-2, 1e-4):
        w2 = (c*sc.Sz[0])*wf
        print("  ||(%.0e*Sz0)|rand>||^2 = %.3e  (<rand|rand>=%.3f)" % (c, (w2.dot(w2)).real, (wf.dot(wf)).real))
```

`kpm/03_herm_verdicts.out`:


```
itensor_version python
  c=1e-01  ED is_hermitian=False (||h-h^dag||_F^2=1.60e-01)  DMRG is_hermitian=False
  c=3e-02  ED is_hermitian=False (||h-h^dag||_F^2=1.44e-02)  DMRG is_hermitian=False
  c=1e-02  ED is_hermitian=False (||h-h^dag||_F^2=1.60e-03)  DMRG is_hermitian=True
  c=3e-03  ED is_hermitian=False (||h-h^dag||_F^2=1.44e-04)  DMRG is_hermitian=True
  c=1e-03  ED is_hermitian=False (||h-h^dag||_F^2=1.60e-05)  DMRG is_hermitian=True
  c=1e-04  ED is_hermitian=False (||h-h^dag||_F^2=1.60e-07)  DMRG is_hermitian=True
  c=3e-05  ED is_hermitian=False (||h-h^dag||_F^2=1.44e-08)  DMRG is_hermitian=True
  c=1e-05  ED is_hermitian=True (||h-h^dag||_F^2=1.60e-09)  DMRG is_hermitian=True
  c=1e-06  ED is_hermitian=True (||h-h^dag||_F^2=1.60e-11)  DMRG is_hermitian=True
  ||(1e+00*Sz0)|rand>||^2 = 2.500e-01  (<rand|rand>=1.000)
  ||(1e-02*Sz0)|rand>||^2 = 2.500e-05  (<rand|rand>=1.000)
  ||(1e-04*Sz0)|rand>||^2 = 2.500e-09  (<rand|rand>=1.000)
itensor_version 3
  c=1e-01  ED is_hermitian=False (||h-h^dag||_F^2=1.60e-01)  DMRG is_hermitian=False
  c=3e-02  ED is_hermitian=False (||h-h^dag||_F^2=1.44e-02)  DMRG is_hermitian=False
  c=1e-02  ED is_hermitian=False (||h-h^dag||_F^2=1.60e-03)  DMRG is_hermitian=False
  c=3e-03  ED is_hermitian=False (||h-h^dag||_F^2=1.44e-04)  DMRG is_hermitian=True
  c=1e-03  ED is_hermitian=False (||h-h^dag||_F^2=1.60e-05)  DMRG is_hermitian=True
  c=1e-04  ED is_hermitian=False (||h-h^dag||_F^2=1.60e-07)  DMRG is_hermitian=True
  c=3e-05  ED is_hermitian=False (||h-h^dag||_F^2=1.44e-08)  DMRG is_hermitian=True
  c=1e-05  ED is_hermitian=True (||h-h^dag||_F^2=1.60e-09)  DMRG is_hermitian=True
  c=1e-06  ED is_hermitian=True (||h-h^dag||_F^2=1.60e-11)  DMRG is_hermitian=True
  ||(1e+00*Sz0)|rand>||^2 = 2.500e-01  (<rand|rand>=1.000)
  ||(1e-02*Sz0)|rand>||^2 = 2.500e-05  (<rand|rand>=1.000)
  ||(1e-04*Sz0)|rand>||^2 = 2.500e-09  (<rand|rand>=1.000)
```

```bash
cd <scratch>/kpm && <scratch>/run3.sh 04_herm_scale.py 2>&1 | tee 04_herm_scale.out
```

`kpm/04_herm_scale.py`:

```python
# The DMRG Hermiticity probe is an absolute threshold, ||(H-H^dag)|w>||^2 >
# 1e-4 on a unit witness, so the same non-Hermitian model written in other
# energy units (H -> s*H) takes the Hermitian branch on DMRG while ED, whose
# own absolute threshold is 1e-8 on ||h-h^dag||_F^2, still sees it as
# non-Hermitian.  What does each route then return for gs_energy and the KPM
# correlator?  Model: 4-site Heisenberg + 0.3*Sz_tot + 1j*gamma*Sz0, gamma=1,
# in units of J, scaled by s.  Scale covariance is the anchor: E0(s*H)=s*E0(H)
# and C_{sH}(s*w) = C_H(w)/s.
import numpy as np
from dmrgpy import spinchain

L = 4
gamma = 1.0
def build(s, v):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L):
        h = h + 0.3*sc.Sz[i]
    h = h + 1j*gamma*sc.Sz[0]
    sc.set_hamiltonian(s*h)
    sc.maxm = 20; sc.nsweeps = 10
    return sc

es = np.linspace(-0.5, 4.0, 200)
ref = {}
for v in ("python", 3):
    for s in (1.0, 1e-2, 1e-3):
        sc = build(s, v)
        hd = sc.is_hermitian(sc.hamiltonian)
        e_ed = sc.gs_energy(mode="ED")
        try:
            e_dm = sc.gs_energy()
        except Exception as e:
            e_dm = "raised %s" % type(e).__name__
        print("v=%s s=%.0e DMRG is_hermitian=%s  E0/s: ED %s  DMRG %s" % (
            v, s, hd, np.round(np.complex128(e_ed)/s, 6),
            np.round(np.complex128(e_dm)/s, 6) if not isinstance(e_dm, str) else e_dm))
        name = (sc.Sz[1], sc.Sz[1])
        # KPM on each mode, default kwargs (as a user would call it)
        for mode in ("DMRG", "ED"):
            try:
                x, y = sc.get_dynamical_correlator(mode=mode, name=name,
                                                   delta=0.1*s, es=es*s)
                y = np.asarray(y)*s
                print("   KPM mode=%s returned: max|C|=%.4f  int C dw=%.4f" % (
                    mode, np.max(np.abs(y)), np.trapezoid(y.real, es)))
                if s == 1e-3: ref[(v, mode)] = y
            except Exception as e:
                print("   KPM mode=%s raised %s: %s" % (mode, type(e).__name__, str(e)[:70]))
        if s == 1e-3:
            # ED non-Hermitian KPM with E_max given (the route ED takes), same scale
            x, y = sc.get_dynamical_correlator(mode="ED", name=name, delta=0.1*s,
                                               es=es*s, E_max=20.*s, n=200)
            y = np.asarray(y)*s
            ref[(v, "EDnh")] = y
            print("   ED NH-KPM (E_max given): max|C|=%.4f  int Re C dw=%.4f" % (
                np.max(np.abs(y)), np.trapezoid(y.real, es)))
    # the Hermitian part alone, gamma=0, for the curve DMRG would give if it
    # simply dropped the non-Hermitian term
    sc0 = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h0 = 0
    for i in range(L-1):
        h0 = h0 + sc0.Sx[i]*sc0.Sx[i+1] + sc0.Sy[i]*sc0.Sy[i+1] + sc0.Sz[i]*sc0.Sz[i+1]
    for i in range(L):
        h0 = h0 + 0.3*sc0.Sz[i]
    sc0.set_hamiltonian(1e-3*h0)
    x, y0 = sc0.get_dynamical_correlator(mode="ED", submode="KPM", name=(sc0.Sz[1], sc0.Sz[1]),
                                         delta=1e-4, es=es*1e-3)
    y0 = np.asarray(y0)*1e-3
    if (v, "DMRG") in ref:
        yd = ref[(v, "DMRG")]
        print("  v=%s s=1e-3: max|DMRG KPM - ED NH-KPM| = %.4f, max|DMRG KPM - gamma=0 ED KPM| = %.4f, peak gamma=0 %.4f" % (
            v, np.max(np.abs(yd-ref[(v, "EDnh")])), np.max(np.abs(yd-y0)), np.max(np.abs(y0))))
```

`kpm/04_herm_scale.out`:


```
v=python s=1e+00 DMRG is_hermitian=False  E0/s: ED (-1.242697-0.291419j)  DMRG (-1.242697+0.291419j)
   KPM mode=DMRG raised ValueError: E_max (an upper bound for the spectral radius of the Hamiltonian) must
   KPM mode=ED raised ValueError: E_max (an upper bound for the spectral radius of the Hamiltonian) must
v=python s=1e-02 DMRG is_hermitian=True  E0/s: ED (-1.242697-0.291419j)  DMRG (-1.777037+0j)
   KPM mode=DMRG raised RuntimeError: KPM moments diverging: scaled spectrum outside [-1,1] (a Chebyshev mom
   KPM mode=ED raised ValueError: E_max (an upper bound for the spectral radius of the Hamiltonian) must
v=python s=1e-03 DMRG is_hermitian=True  E0/s: ED (-1.242697-0.291419j)  DMRG (-2.003153+0j)
   KPM mode=DMRG raised RuntimeError: KPM moments diverging: scaled spectrum outside [-1,1] (a Chebyshev mom
   KPM mode=ED raised ValueError: E_max (an upper bound for the spectral radius of the Hamiltonian) must
Traceback (most recent call last):
  File "<scratch>/kpm/04_herm_scale.py", line 54, in <module>
    x, y = sc.get_dynamical_correlator(mode="ED", name=name, delta=0.1*s,
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                       es=es*s, E_max=20.*s, n=200)
                                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "<repo>/src/dmrgpy/manybodychain.py", line 1032, in get_dynamical_correlator
    return edobj.get_dynamical_correlator(**kwargs)
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^
  File "<repo>/src/dmrgpy/edtk/edchain.py", line 293, in get_dynamical_correlator
    return dynamics.get_dynamical_correlator(self,**kwargs)
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^
  File "<repo>/src/dmrgpy/edtk/dynamics.py", line 87, in get_dynamical_correlator
    return dynamical_correlator_nhkpm_ed(self,name=name,**kwargs)
  File "<repo>/src/dmrgpy/nonhermitian/kpm.py", line 79, in dynamical_correlator_nhkpm_ed
    e0,vr,vl = algebra.biorthogonal_ground_state(h)
               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^
  File "<repo>/src/dmrgpy/algebra/algebra.py", line 339, in biorthogonal_ground_state
    raise ValueError("Right and left ground states are (near) "
            "orthogonal, <vl|vr> = "+str(norm)+" - biorthogonal "
            "normalization failed")
ValueError: Right and left ground states are (near) orthogonal, <vl|vr> = (-1.3877787807814457e-15+4.996003610813204e-16j) - biorthogonal normalization failed
```

```bash
cd <scratch>/kpm && <scratch>/run3.sh 05_herm_scale_anchor.py 2>&1 | tee 05_herm_scale_anchor.out
```

`kpm/05_herm_scale_anchor.py`:

```python
# Anchor for 04 with plain numpy: the exact spectrum of H = Heis4 + 0.3*Sz_tot
# + 1j*gamma*Sz0 (dense eig of the ED matrix), the ground state of its
# Hermitian part, and what the DMRG route returns once the model is written in
# units where the probe calls it Hermitian (s*H).  Then the weak-gamma case,
# where the moment guard does not fire: does DMRG KPM return a curve, and how
# far is it from the gamma=0 curve (i.e. does the non-Hermitian term enter)?
import numpy as np
from dmrgpy import spinchain

L = 4
def build(s, gamma, v="python"):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L):
        h = h + 0.3*sc.Sz[i]
    hr = h
    h = h + 1j*gamma*sc.Sz[0]
    sc.set_hamiltonian(s*h)
    sc.maxm = 20; sc.nsweeps = 10
    return sc, hr

sc, hr = build(1.0, 1.0)
H = sc.get_ED_obj().get_hamiltonian().toarray()
ev = np.linalg.eigvals(H)
ev = ev[np.argsort(ev.real)]
print("exact eigenvalues of H (gamma=1), lowest real parts:", np.round(ev[:4], 6))
Hr = (H + H.conj().T)/2
print("lowest eigenvalue of the Hermitian part (H+H^dag)/2:", np.round(np.linalg.eigvalsh(Hr)[0], 6))
for v in ("python", 3):
    for s in (1.0, 1e-2, 1e-3):
        sc, _ = build(s, 1.0, v)
        print("v=%s s=%.0e is_hermitian=%s DMRG E0/s=%s" % (
            v, s, sc.is_hermitian(sc.hamiltonian), np.round(np.complex128(sc.gs_energy())/s, 6)))

es = np.linspace(-0.5, 4.0, 200)
s = 1e-3
for gamma in (0.0, 0.03, 0.1):
    sc, _ = build(s, gamma)
    herm = sc.is_hermitian(sc.hamiltonian)
    try:
        x, y = sc.get_dynamical_correlator(name=(sc.Sz[1], sc.Sz[1]), delta=0.1*s, es=es*s)
        y = np.asarray(y)*s
        out = "returned, max|C|=%.4f int=%.4f" % (np.max(np.abs(y)), np.trapezoid(y.real, es))
    except Exception as e:
        y = None
        out = "raised %s" % type(e).__name__
    if gamma == 0.0: y0 = y
    diff = "" if (y is None or y0 is None or gamma == 0.0) else " max|C-C(gamma=0)|=%.4f" % np.max(np.abs(y-y0))
    Hm = sc.get_ED_obj().get_hamiltonian().toarray()/s
    ev = np.linalg.eigvals(Hm); ev = ev[np.argsort(ev.real)]
    print("gamma=%.2f s=1e-3 DMRG is_hermitian=%s E0_DMRG/s=%.6f exact E0=%s  KPM %s%s" % (
        gamma, herm, np.real(sc.gs_energy())/s, np.round(ev[0], 6), out, diff))
```

`kpm/05_herm_scale_anchor.out`:


```
exact eigenvalues of H (gamma=1), lowest real parts: [-1.242697+0.291419j -1.242697-0.291419j -1.219455-0.456421j
 -0.619455+0.456421j]
lowest eigenvalue of the Hermitian part (H+H^dag)/2: -1.616025
v=python s=1e+00 is_hermitian=False DMRG E0/s=(-1.242697+0.291419j)
v=python s=1e-02 is_hermitian=True DMRG E0/s=(-1.540751+0j)
v=python s=1e-03 is_hermitian=True DMRG E0/s=(-1.741346+0j)
v=3 s=1e+00 is_hermitian=False DMRG E0/s=(-1.242697-0.291419j)
v=3 s=1e-02 is_hermitian=False DMRG E0/s=(-1.242697+0.291419j)
v=3 s=1e-03 is_hermitian=True DMRG E0/s=(-1.184871+0j)
gamma=0.00 s=1e-3 DMRG is_hermitian=True E0_DMRG/s=-1.616025 exact E0=(-1.616025+0j)  KPM returned, max|C|=0.5528 int=0.2500
gamma=0.03 s=1e-3 DMRG is_hermitian=True E0_DMRG/s=-1.615745 exact E0=(-1.615745-0j)  KPM returned, max|C|=0.5564 int=0.2499 max|C-C(gamma=0)|=0.0044
gamma=0.10 s=1e-3 DMRG is_hermitian=True E0_DMRG/s=-1.612903 exact E0=(-1.612903+0j)  KPM returned, max|C|=0.5626 int=0.2482 max|C-C(gamma=0)|=0.0475
```

```bash
cd <scratch>/kpm && <scratch>/run3.sh 17_probe_consolidation.py 2>&1 | tee 17_probe_consolidation.out
```

`kpm/17_probe_consolidation.py`:

```python
# Consolidation of 03-05: H = Heis4 + 0.3*Sz_tot + 1j*Sz0 written in units s.
# Exact E0 from np.linalg.eigvals of the ED matrix (lowest real part, a
# conjugate pair); gs_energy()/s on each DMRG backend, three fresh chains per
# point, with the probe's verdict; ED for reference; the KPM correlator's
# outcome on DMRG.
import numpy as np
from dmrgpy import spinchain

L = 4
es = np.linspace(-0.5, 4.0, 200)
def build(s, v):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L):
        h = h + 0.3*sc.Sz[i]
    h = h + 1j*sc.Sz[0]
    sc.set_hamiltonian(s*h)
    sc.maxm = 20; sc.nsweeps = 10
    return sc

sc = build(1.0, "python")
ev = np.linalg.eigvals(sc.get_ED_obj().get_hamiltonian().toarray())
ev = ev[np.argsort(ev.real)]
print("exact lowest eigenvalues:", np.round(ev[:3], 6))
for v in ("python", 3, 2):
    for s in (1.0, 1e-2, 1e-3):
        verdicts, es0 = [], []
        for rep in range(3):
            sc = build(s, v)
            verdicts.append(sc.is_hermitian(sc.hamiltonian))
            try:
                es0.append("%s" % np.round(np.complex128(sc.gs_energy())/s, 4))
            except Exception as e:
                es0.append("raised %s" % type(e).__name__)
        eed = np.complex128(sc.gs_energy(mode="ED"))/s
        try:
            sc.get_dynamical_correlator(name=(sc.Sz[1], sc.Sz[1]), delta=0.1*s, es=es*s)
            kp = "returned"
        except Exception as e:
            kp = "%s: %s" % (type(e).__name__, str(e)[:45])
        print("v=%-6s s=%.0e probe=%s DMRG E0/s=%s  ED E0/s=%s  KPM: %s" % (
            v, s, verdicts, es0, np.round(eed, 4), kp))
```

`kpm/17_probe_consolidation.out`:


```
exact lowest eigenvalues: [-1.242697+0.291419j -1.242697-0.291419j -1.219455-0.456421j]
v=python s=1e+00 probe=[False, False, False] DMRG E0/s=['(-1.2427+0.2914j)', '(-1.2427+0.2914j)', '(-1.2427-0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=python s=1e-02 probe=[True, True, True] DMRG E0/s=['(-1.6263+0j)', '(-1.9858+0j)', '(-2.0343+0j)']  ED E0/s=(-1.2427-0.2914j)  KPM: RuntimeError: KPM moments diverging: scaled spectrum outsid
v=python s=1e-03 probe=[True, True, True] DMRG E0/s=['(-2.0363+0j)', '(-2.0417+0j)', '(-1.7919+0j)']  ED E0/s=(-1.2427-0.2914j)  KPM: RuntimeError: KPM moments diverging: scaled spectrum outsid
v=3      s=1e+00 probe=[False, False, False] DMRG E0/s=['(-1.2427+0.2914j)', '(-1.2427+0.2914j)', '(-1.2427-0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=3      s=1e-02 probe=[False, True, False] DMRG E0/s=['(-1.2427+0.2914j)', '(-1.2204+0j)', '(-1.2427+0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=3      s=1e-03 probe=[True, True, True] DMRG E0/s=['(-1.2087+0j)', '(-1.2014+0j)', '(-1.2187+0j)']  ED E0/s=(-1.2427-0.2914j)  KPM: RuntimeError: KPM moments diverging: scaled spectrum outsid
v=2      s=1e+00 probe=[False, False, False] DMRG E0/s=['(-1.2427-0.2914j)', '(-1.2427-0.2914j)', '(-1.2427-0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=2      s=1e-02 probe=[False, True, False] DMRG E0/s=['(-1.2427-0.2914j)', '(-1.2548+0j)', '(-1.2427-0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=2      s=1e-03 probe=[True, True, True] DMRG E0/s=['(-1.2155+0j)', '(-1.2115+0j)', '(-1.2254+0j)']  ED E0/s=(-1.2427-0.2914j)  KPM: RuntimeError: KPM moments diverging: scaled spectrum outsid
```


and from the `misc` hunter (`01` the probe norm and `disentangle_manifold` against
`eig`, `02`/`08` `gs_energy()` on the rescaled Hatano-Nelson chain on `"python"`
and v3):

```bash
cd <scratch>/misc && <scratch>/run3.sh 01_probe_threshold.py 2>&1 | tee 01_probe_threshold.out
```

`misc/01_probe_threshold.py`:

```python
"""The chain's Hermiticity probe ends in `return not norm>1e-4`, with norm
= ||(A - A^dag)|w>||^2 for a random witness w. Nothing rescales by the size
of A, so a small enough non-Hermitian A is called Hermitian. 867e2b4 routed
disentangle_manifold onto that probe, and its Hermitian branch then
diagonalizes (ma+ma^dag)/2, i.e. the Hermitian part of A, instead of A.

Anchor: A = eps*(Sz0 + S+0) + 0.3*eps*Sz1 is non-Hermitian by construction
(A - A^dag = eps*(S+0 - S-0) = 2i*eps*Sy0 != 0) and its exact right
eigenvectors on the full 4-dim Hilbert space of 2 spins are those of the
4x4 matrix ma = <w_i|A|w_j> on an orthonormal basis; the pre-867e2b4
decision (bare proof -> False -> eig) returned exactly those."""
import warnings
import numpy as np
import scipy.linalg as dlg
from dmrgpy import spinchain, mpsalgebra
from dmrgpy.mpsalgebratk.disentangle import get_representation

warnings.simplefilter("ignore")
sc = spinchain.Spin_Chain([2]*3, itensor_version="python")
h = sc.Sx[0]*sc.Sx[1] + 0.6*sc.Sy[1]*sc.Sy[2] + 0.3*sc.Sz[0]*sc.Sz[1] \
    + 0.7*sc.Sx[0] + 0.45*sc.Sz[1] + 0.2*sc.Sy[2] + 0.33*sc.Sz[2]
sc.set_hamiltonian(h)
sc.maxm, sc.nsweeps = 8, 20
np.random.seed(0)
es, wfs = sc.get_excited_states(n=8)   # full 8-dim Hilbert space of 3 spins
G = np.array([[a.dot(b) for b in wfs] for a in wfs])
print("manifold Gram error %.2e" % np.max(np.abs(G - np.eye(8))))

def probe_norm(A, seed):
    np.random.seed(seed)
    old = sc.maxm; sc.maxm = min(old, 8)
    try:
        w = sc.random_mps()
        d = (A - A.get_dagger())*w
        return (d.dot(d)).real
    finally:
        sc.maxm = old

def eig_residual(wfs, out, A):
    """max over output states of ||(A - lambda) v|| / ||v|| in the
    manifold's own coordinates, lambda the Rayleigh quotient"""
    ma = get_representation(wfs, A)
    C = np.array([[a.dot(b) for b in out] for a in wfs])  # C[i,j] = <w_i|out_j>
    worst = 0.
    for j in range(C.shape[1]):
        v = C[:, j]
        lam = (v.conj() @ ma @ v)/(v.conj() @ v)
        worst = max(worst, np.linalg.norm(ma @ v - lam*v)/np.linalg.norm(v))
    return worst, np.max(np.abs(ma))

print("%8s %6s %12s %10s %12s %12s" % ("eps", "proof", "probe_norm", "probe",
                                        "resid_new", "resid_eig"))
for eps in [1.0, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 1e-4]:
    A = eps*(sc.Sz[0] + sc.Sx[0] + 1j*sc.Sy[0]) + 0.3*eps*sc.Sz[1]
    pn = probe_norm(A, 7)
    np.random.seed(7)
    verdict = sc.is_hermitian(A)
    np.random.seed(7)
    out = mpsalgebra.disentangle_manifold(wfs, A)
    r_new, amax = eig_residual(wfs, out, A)
    # the pre-867e2b4 decision, rebuilt: bare proof -> False -> eig
    ma = get_representation(wfs, A)
    ev, V = dlg.eig(ma)
    r_old = max(np.linalg.norm(ma @ V[:, j] - ev[j]*V[:, j]) for j in range(8))
    print("%8.0e %6s %12.3e %10s %12.3e %12.3e   (rel. to max|ma|=%.2e: %.3f)"
          % (eps, A.is_hermitian(), pn, verdict, r_new, r_old, amax, r_new/amax))

# seed dependence of the verdict at the flip region
for eps in [3e-2, 1e-2]:
    votes = []
    for seed in range(20):
        A = eps*(sc.Sz[0] + sc.Sx[0] + 1j*sc.Sy[0]) + 0.3*eps*sc.Sz[1]
        np.random.seed(seed)
        votes.append(sc.is_hermitian(A))
    print("eps=%.0e: probe says Hermitian in %d of 20 seeds" % (eps, sum(votes)))
```

`misc/01_probe_threshold.out`:


```
WARNING, state is not normalizable. Returning None
WARNING, state is not normalizable. Returning None
manifold Gram error 6.59e-15
     eps  proof   probe_norm      probe    resid_new    resid_eig
   1e+00  False    1.000e+00      False    4.066e-15    1.420e-15   (rel. to max|ma|=7.13e-01: 0.000)
   1e-01  False    1.000e-02      False    4.354e-16    9.958e-17   (rel. to max|ma|=7.13e-02: 0.000)
   3e-02  False    9.000e-04      False    1.145e-16    3.065e-17   (rel. to max|ma|=2.14e-02: 0.000)
   1e-02  False    1.000e-04       True    5.000e-03    1.052e-17   (rel. to max|ma|=7.13e-03: 0.702)
   3e-03  False    9.000e-06       True    1.500e-03    3.959e-18   (rel. to max|ma|=2.14e-03: 0.702)
   1e-03  False    1.000e-06       True    5.000e-04    1.026e-18   (rel. to max|ma|=7.13e-04: 0.702)
   1e-04  False    1.000e-08       True    5.000e-05    1.232e-19   (rel. to max|ma|=7.13e-05: 0.702)
eps=3e-02: probe says Hermitian in 0 of 20 seeds
eps=1e-02: probe says Hermitian in 18 of 20 seeds
```

```bash
cd <scratch>/misc && <scratch>/run3.sh 02_probe_gs_energy.py 2>&1 | tee 02_probe_gs_energy.out
```

`misc/02_probe_gs_energy.py`:

```python
"""The same absolute probe threshold, one level up: gs_energy() dispatches
Hermitian DMRG vs NH-DMRG on self.is_hermitian(self.hamiltonian), and the
proof cannot see a non-Hermitian H, so the probe decides. A Hamiltonian
written in small units (eV for a meV-scale model, as the Kondo code does)
has ||(H-H^dag)|w>||^2 below 1e-4 whatever its relative non-Hermiticity.

Anchor: the Hatano-Nelson (non-reciprocal XX) chain under open boundaries,
H = s * sum_i [(1+g) S+_i S-_{i+1} + (1-g) S-_i S+_{i+1}] / 2,
is similar to the Hermitian XX chain with hopping sqrt(1-g^2), so its
spectrum is real and E0(H) = sqrt(1-g^2) * E0(Hermitian part) exactly;
checked against a dense eig of the 2^n matrix below."""
import warnings
import numpy as np
from dmrgpy import spinchain

warnings.simplefilter("ignore")
n, g = 6, 0.6
sp = np.array([[0, 1], [0, 0]], dtype=complex); sm = sp.T.copy()
def site_op(o, i):
    out = np.array([[1.0+0j]])
    for k in range(n):
        out = np.kron(out, o if k == i else np.eye(2))
    return out
Hd = sum((1+g)*site_op(sp, i) @ site_op(sm, i+1) + (1-g)*site_op(sm, i) @ site_op(sp, i+1)
         for i in range(n-1))/2
ev = np.linalg.eigvals(Hd)
e_exact = ev[np.argmin(ev.real)]
eh = np.linalg.eigvalsh((Hd + Hd.conj().T)/2)[0]
print("exact E0(H)/s = %.10f%+.1ej, Hermitian part E0/s = %.10f, ratio %.6f (sqrt(1-g^2) = %.6f)"
      % (e_exact.real, e_exact.imag, eh, e_exact.real/eh, np.sqrt(1-g*g)))

print("%8s %6s %8s %16s %12s" % ("s", "proof", "probe", "gs_energy()/s", "rel.err"))
for s in [1.0, 1e-1, 1e-2, 1e-3]:
    sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
    Sp = [sc.Sx[i] + 1j*sc.Sy[i] for i in range(n)]
    Sm = [sc.Sx[i] - 1j*sc.Sy[i] for i in range(n)]
    h = 0
    for i in range(n-1):
        h = h + s*((1+g)*Sp[i]*Sm[i+1] + (1-g)*Sm[i]*Sp[i+1])/2
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 10
    np.random.seed(1)
    herm = sc.is_hermitian(sc.hamiltonian)
    np.random.seed(1)
    e = sc.gs_energy()
    print("%8.0e %6s %8s %16.10f %12.3e" % (s, h.is_hermitian(), herm, np.real(e)/s,
                                           abs(e/s - e_exact)/abs(e_exact)))
```

`misc/02_probe_gs_energy.out`:


```
exact E0(H)/s = -1.3975836830+0.0e+00j, Hermitian part E0/s = -1.7469796037, ratio 0.800000 (sqrt(1-g^2) = 0.800000)
       s  proof    probe    gs_energy()/s      rel.err
   1e+00  False    False    -1.3975836830    4.290e-15
   1e-01  False    False    -1.3975836830    3.655e-15
   1e-02  False     True    -1.5008764226    7.391e-02
   1e-03  False     True    -1.4136759840    1.151e-02
```

```bash
cd <scratch>/misc && <scratch>/run3.sh 08_probe_gs_energy_v3.py 2>&1 | tee 08_probe_gs_energy_v3.out
```

`misc/08_probe_gs_energy_v3.py`:

```python
"""v3 copy of 02. The same absolute probe threshold, one level up: gs_energy() dispatches
Hermitian DMRG vs NH-DMRG on self.is_hermitian(self.hamiltonian), and the
proof cannot see a non-Hermitian H, so the probe decides. A Hamiltonian
written in small units (eV for a meV-scale model, as the Kondo code does)
has ||(H-H^dag)|w>||^2 below 1e-4 whatever its relative non-Hermiticity.

Anchor: the Hatano-Nelson (non-reciprocal XX) chain under open boundaries,
H = s * sum_i [(1+g) S+_i S-_{i+1} + (1-g) S-_i S+_{i+1}] / 2,
is similar to the Hermitian XX chain with hopping sqrt(1-g^2), so its
spectrum is real and E0(H) = sqrt(1-g^2) * E0(Hermitian part) exactly;
checked against a dense eig of the 2^n matrix below."""
import warnings
import numpy as np
from dmrgpy import spinchain

warnings.simplefilter("ignore")
n, g = 6, 0.6
sp = np.array([[0, 1], [0, 0]], dtype=complex); sm = sp.T.copy()
def site_op(o, i):
    out = np.array([[1.0+0j]])
    for k in range(n):
        out = np.kron(out, o if k == i else np.eye(2))
    return out
Hd = sum((1+g)*site_op(sp, i) @ site_op(sm, i+1) + (1-g)*site_op(sm, i) @ site_op(sp, i+1)
         for i in range(n-1))/2
ev = np.linalg.eigvals(Hd)
e_exact = ev[np.argmin(ev.real)]
eh = np.linalg.eigvalsh((Hd + Hd.conj().T)/2)[0]
print("exact E0(H)/s = %.10f%+.1ej, Hermitian part E0/s = %.10f, ratio %.6f (sqrt(1-g^2) = %.6f)"
      % (e_exact.real, e_exact.imag, eh, e_exact.real/eh, np.sqrt(1-g*g)))

print("%8s %6s %8s %16s %12s" % ("s", "proof", "probe", "gs_energy()/s", "rel.err"))
for s in [1.0, 1e-1, 1e-2, 1e-3]:
    sc = spinchain.Spin_Chain([2]*n, itensor_version=3)
    Sp = [sc.Sx[i] + 1j*sc.Sy[i] for i in range(n)]
    Sm = [sc.Sx[i] - 1j*sc.Sy[i] for i in range(n)]
    h = 0
    for i in range(n-1):
        h = h + s*((1+g)*Sp[i]*Sm[i+1] + (1-g)*Sm[i]*Sp[i+1])/2
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 20, 10
    np.random.seed(1)
    herm = sc.is_hermitian(sc.hamiltonian)
    np.random.seed(1)
    e = sc.gs_energy()
    print("%8.0e %6s %8s %16.10f %12.3e" % (s, h.is_hermitian(), herm, np.real(e)/s,
                                           abs(e/s - e_exact)/abs(e_exact)))
```

`misc/08_probe_gs_energy_v3.out`:


```
exact E0(H)/s = -1.3975836830+0.0e+00j, Hermitian part E0/s = -1.7469796037, ratio 0.800000 (sqrt(1-g^2) = 0.800000)
       s  proof    probe    gs_energy()/s      rel.err
   1e+00  False    False    -1.3975836830    3.654e-15
   1e-01  False    False    -1.3975836830    3.336e-15
   1e-02  False     True    -1.4437449925    3.303e-02
   1e-03  False     True    -1.2631516711    9.619e-02
```


**Reviewer of the `kpm` candidate (CONFIRMED, NARROWED)**: `17` reproduces within
run-to-run noise: at s=1 the probe says False everywhere and DMRG gives E0/s =
-1.2427 +- 0.2914i on all three backends; at s=1e-3 it says True everywhere, and
`"python"` gives -1.78, -1.5788, -2.0909, v3 -1.2151, -1.2364, -1.2189, v2
-1.2758, -1.2152, -1.2326, all real; the `"python"` values reach below the lowest
eigenvalue of the Hermitian part (about -1.616), so they are not even Rayleigh
quotients of H. On the probe's own `maxm=8` witness, ||dh w||^2 for s*(Heis4 +
1j*Sz0) is exactly s^2 (1e6, 1, 1e-6) on all three backends with ||w||^2 =
1.000000. Away from the knife edge the verdict is witness-dependent (1j*c*(Sz0+Sz1),
20 probes per point: `"python"` True in 18, 9 and 4 of 20 at c = 6e-3, 7e-3, 8e-3,
v3 in 20, 9 and 6). In natural units with a first-order decay rate (Heis4 +
1.5*Sz_tot + 1j*c*Sz0, exact E0 = -2.457105 - 0.002134i at c=5e-3), the probe says
True on every backend, v3 and v2 return -2.457105 + 0i (the decay rate lost),
`"python"` -2.887505 and -2.537763 (18 and 3 per cent off, below the exact lowest
real part); at c=2e-2 every backend is exact. The Hermitian KPM on that chain
returns silently, on `"python"` off by the full peak height (0.612 on a 0.61 peak),
while `mode="ED"` takes the non-Hermitian branch. The counterweight: when the exact
E0 is real (no first-order expectation value of the anti-Hermitian term), the
Hermitian solver is right to 1e-6 or better. The two ED checks disagree with each
other on the same matrix at c=3e-6 (`is_hermitian` True, `ishermitian` False).
Struck, each explicitly:

- "The same Hamiltonian written in smaller energy units" as the defect: rescaling
  is one instance; the mechanism is an absolute threshold in whatever units H is
  written in, so a J=1 chain with a weak dissipative term is misdispatched too.
- "45 to 64 per cent in the real part on `"python"`": over two sessions 21 to 68
  per cent; v3 -1.20 to -1.29 and v2 -1.21 to -1.28.
- "KPM on the misdispatched chains is not silent, it raises the moment guard":
  only on the rescaled chain; in natural units the Hermitian KPM returns silently.
- "The verdict flips between repeated runs at s=1e-2": there ||dh w||^2 is 1e-4
  identically, the threshold itself, so that flip is roundoff; genuine witness
  dependence exists regardless (9 of 20 at c=7e-3).
- An implied general harm: when the exact E0 is real, the answer is right.

```bash
cd <scratch>/reviews/kpm_1 && <scratch>/run3.sh 01_repro_consolidation.py 2>&1 | tee 01_repro_consolidation.out
```

`reviews/kpm_1/01_repro_consolidation.py` is `kpm/17_probe_consolidation.py` (finding 12) unchanged, rerun; its output, `reviews/kpm_1/01_repro_consolidation.out`:


```
exact lowest eigenvalues: [-1.242697+0.291419j -1.242697-0.291419j -1.219455-0.456421j]
v=python s=1e+00 probe=[False, False, False] DMRG E0/s=['(-1.2427-0.2914j)', '(-1.2427+0.2914j)', '(-1.2427-0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=python s=1e-02 probe=[True, True, True] DMRG E0/s=['(-2.0874+0j)', '(-1.5035+0j)', '(-2.0628+0j)']  ED E0/s=(-1.2427-0.2914j)  KPM: RuntimeError: KPM moments diverging: scaled spectrum outsid
v=python s=1e-03 probe=[True, True, True] DMRG E0/s=['(-1.78+0j)', '(-1.5788+0j)', '(-2.0909+0j)']  ED E0/s=(-1.2427-0.2914j)  KPM: RuntimeError: KPM moments diverging: scaled spectrum outsid
v=3      s=1e+00 probe=[False, False, False] DMRG E0/s=['(-1.2427+0.2914j)', '(-1.2427-0.2914j)', '(-1.2427-0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=3      s=1e-02 probe=[False, True, False] DMRG E0/s=['(-1.2427+0.2914j)', '(-1.2934+0j)', '(-1.2427+0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=3      s=1e-03 probe=[True, True, True] DMRG E0/s=['(-1.2151+0j)', '(-1.2364+0j)', '(-1.2189+0j)']  ED E0/s=(-1.2427-0.2914j)  KPM: RuntimeError: KPM moments diverging: scaled spectrum outsid
v=2      s=1e+00 probe=[False, False, False] DMRG E0/s=['(-1.2427-0.2914j)', '(-1.2427-0.2914j)', '(-1.2427-0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=2      s=1e-02 probe=[False, True, False] DMRG E0/s=['(-1.2427-0.2914j)', '(-1.2129+0j)', '(-1.2427-0.2914j)']  ED E0/s=(-1.2427-0.2914j)  KPM: ValueError: E_max (an upper bound for the spectral radius
v=2      s=1e-03 probe=[True, True, True] DMRG E0/s=['(-1.2758+0j)', '(-1.2152+0j)', '(-1.2326+0j)']  ED E0/s=(-1.2427-0.2914j)  KPM: RuntimeError: KPM moments diverging: scaled spectrum outsid
```

```bash
cd <scratch>/reviews/kpm_1 && <scratch>/run3.sh 02_attacks.py 2>&1 | tee 02_attacks.out
```

`reviews/kpm_1/02_attacks.py`:

```python
# Reviewer attacks on kpm_1.
# (a) natural units, s=1: a weak anti-Hermitian term i*c*Sz0 on the J=1 chain.
#     What does the Hermitian solver return when the probe misdispatches, and
#     how far is that from the exact lowest eigenvalue?
# (b) witness dependence: dh = 2i*c*(Sz0+Sz1), whose square is not a constant,
#     so the probe's verdict near threshold depends on the random witness.
# (c) the suggested relative test: ||dh w||^2 / ||H w||^2 on an unproven
#     Hermitian H (1j*Sx0*Sy0 = -Sz0/2 on spin-1/2) and on the non-Hermitian
#     one, at three scales, on every backend.
import numpy as np
from dmrgpy import spinchain
from dmrgpy import mpsalgebra

L = 4
def heis(sc, field=0.3):
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L):
        h = h + field*sc.Sz[i]
    return h

def chain(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    sc.maxm = 20; sc.nsweeps = 10
    return sc

def exact(sc):
    ev = np.linalg.eigvals(sc.get_ED_obj().get_hamiltonian().toarray())
    return ev[np.argsort(ev.real)][0]

print("(a) s=1, H = Heis4 + 0.3*Sz_tot + 1j*c*Sz0")
for v in ("python", 3, 2):
    for c in (5e-3, 2e-3):
        out = []
        for rep in range(3):
            sc = chain(v)
            sc.set_hamiltonian(heis(sc) + 1j*c*sc.Sz[0])
            herm = sc.is_hermitian(sc.hamiltonian)
            e = np.complex128(sc.gs_energy())
            out.append((herm, e))
        e0 = exact(sc)
        print("  v=%-6s c=%.0e exact=%s | %s" % (v, c, np.round(e0, 6),
              ["herm=%s E=%s err=%.1e" % (h, np.round(e, 6), abs(e-e0)) for h, e in out]))

print("(b) witness dependence, H = Heis4 + 1j*c*(Sz0+Sz1), 20 probes per point")
for v in ("python", 3):
    for c in (6e-3, 7e-3, 8e-3):
        sc = chain(v)
        H = heis(sc) + 1j*c*(sc.Sz[0]+sc.Sz[1])
        sc.set_hamiltonian(H)
        verdicts = [mpsalgebra.is_hermitian(sc, H) for _ in range(20)]
        print("  v=%-6s c=%.0e  True in %d of 20" % (v, c, sum(verdicts)))

print("(c) ||dh w||^2 and ||H w||^2 on the probe's own witness (maxm=8)")
for v in ("python", 3, 2):
    for label, extra in (("unproven Hermitian 1j*Sx0*Sy0", lambda sc: 1j*sc.Sx[0]*sc.Sy[0]),
                         ("non-Hermitian 1j*Sz0", lambda sc: 1j*sc.Sz[0])):
        for s in (1e3, 1.0, 1e-3):
            sc = chain(v)
            H = s*(heis(sc) + extra(sc))
            proven = H.is_hermitian()
            dh = H - H.get_dagger()
            sc.maxm = 8
            w = sc.random_mps()
            a = dh*w; b = H*w
            na = a.dot(a).real if a is not None else 0.0
            nb = b.dot(b).real
            print("  v=%-6s %-31s s=%.0e proven=%s ||w||^2=%.6f ||dh w||^2=%.3e ||H w||^2=%.3e ratio=%.3e" % (
                v, label, s, proven, w.dot(w).real, na, nb, na/nb))
```

`reviews/kpm_1/02_attacks.out`:


```
(a) s=1, H = Heis4 + 0.3*Sz_tot + 1j*c*Sz0
  v=python c=5e-03 exact=(-1.616018-0j) | ['herm=True E=(-1.616018+0j) err=2.2e-16', 'herm=True E=(-1.616018+0j) err=6.7e-16', 'herm=True E=(-1.616018+0j) err=8.9e-16']
  v=python c=2e-03 exact=(-1.616024-0j) | ['herm=True E=(-1.616024+0j) err=4.7e-15', 'herm=True E=(-1.616024+0j) err=4.7e-15', 'herm=True E=(-1.616024+0j) err=4.7e-15']
  v=3      c=5e-03 exact=(-1.616018-0j) | ['herm=True E=(-1.616018+0j) err=1.1e-08', 'herm=True E=(-1.616018+0j) err=1.1e-08', 'herm=True E=(-1.616018+0j) err=2.1e-09']
  v=3      c=2e-03 exact=(-1.616024-0j) | ['herm=True E=(-1.616024+0j) err=7.5e-08', 'herm=True E=(-1.616024+0j) err=6.1e-07', 'herm=True E=(-1.616024+0j) err=1.4e-07']
  v=2      c=5e-03 exact=(-1.616018-0j) | ['herm=True E=(-1.616017+0j) err=4.0e-07', 'herm=True E=(-1.616017+0j) err=1.0e-06', 'herm=True E=(-1.616016+0j) err=1.4e-06']
  v=2      c=2e-03 exact=(-1.616024-0j) | ['herm=True E=(-1.616024+0j) err=1.2e-07', 'herm=True E=(-1.616024+0j) err=4.7e-08', 'herm=True E=(-1.616024+0j) err=8.5e-08']
(b) witness dependence, H = Heis4 + 1j*c*(Sz0+Sz1), 20 probes per point
  v=python c=6e-03  True in 18 of 20
  v=python c=7e-03  True in 9 of 20
  v=python c=8e-03  True in 4 of 20
  v=3      c=6e-03  True in 20 of 20
  v=3      c=7e-03  True in 9 of 20
  v=3      c=8e-03  True in 6 of 20
(c) ||dh w||^2 and ||H w||^2 on the probe's own witness (maxm=8)
  v=python unproven Hermitian 1j*Sx0*Sy0   s=1e+03 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=4.426e+05 ratio=0.000e+00
  v=python unproven Hermitian 1j*Sx0*Sy0   s=1e+00 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=5.953e-01 ratio=0.000e+00
  v=python unproven Hermitian 1j*Sx0*Sy0   s=1e-03 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=4.895e-07 ratio=0.000e+00
  v=python non-Hermitian 1j*Sz0            s=1e+03 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e+06 ||H w||^2=1.398e+06 ratio=7.155e-01
  v=python non-Hermitian 1j*Sz0            s=1e+00 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e+00 ||H w||^2=1.181e+00 ratio=8.465e-01
  v=python non-Hermitian 1j*Sz0            s=1e-03 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e-06 ||H w||^2=6.803e-07 ratio=1.470e+00
  v=3      unproven Hermitian 1j*Sx0*Sy0   s=1e+03 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=5.432e+05 ratio=0.000e+00
  v=3      unproven Hermitian 1j*Sx0*Sy0   s=1e+00 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=6.465e-01 ratio=0.000e+00
  v=3      unproven Hermitian 1j*Sx0*Sy0   s=1e-03 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=6.908e-07 ratio=0.000e+00
  v=3      non-Hermitian 1j*Sz0            s=1e+03 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e+06 ||H w||^2=1.086e+06 ratio=9.205e-01
  v=3      non-Hermitian 1j*Sz0            s=1e+00 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e+00 ||H w||^2=6.619e-01 ratio=1.511e+00
  v=3      non-Hermitian 1j*Sz0            s=1e-03 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e-06 ||H w||^2=8.248e-07 ratio=1.212e+00
  v=2      unproven Hermitian 1j*Sx0*Sy0   s=1e+03 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=6.969e+05 ratio=0.000e+00
  v=2      unproven Hermitian 1j*Sx0*Sy0   s=1e+00 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=4.517e-01 ratio=0.000e+00
  v=2      unproven Hermitian 1j*Sx0*Sy0   s=1e-03 proven=False ||w||^2=1.000000 ||dh w||^2=0.000e+00 ||H w||^2=4.919e-07 ratio=0.000e+00
  v=2      non-Hermitian 1j*Sz0            s=1e+03 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e+06 ||H w||^2=7.167e+05 ratio=1.395e+00
  v=2      non-Hermitian 1j*Sz0            s=1e+00 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e+00 ||H w||^2=7.472e-01 ratio=1.338e+00
  v=2      non-Hermitian 1j*Sz0            s=1e-03 proven=False ||w||^2=1.000000 ||dh w||^2=1.000e-06 ||H w||^2=7.222e-07 ratio=1.385e+00
```

```bash
cd <scratch>/reviews/kpm_1 && <scratch>/run3.sh 03_first_order_decay.py 2>&1 | tee 03_first_order_decay.out
```

`reviews/kpm_1/03_first_order_decay.py`:

```python
# Natural units, s=1, a weak anti-Hermitian term whose ground-state expectation
# is first order: H = Heis4 + 1.5*Sz_tot + 1j*c*Sz0, ground state fully
# polarized down, so exact E0 ~ -2.25 - 1j*c/2. Does the misdispatched
# Hermitian solver lose the decay rate Im E0?
import numpy as np
from dmrgpy import spinchain

L = 4
for v in ("python", 3, 2):
    for c in (5e-3, 2e-2):
        out = []
        for rep in range(2):
            sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
            h = 0
            for i in range(L-1):
                h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(L):
                h = h + 1.5*sc.Sz[i]
            h = h + 1j*c*sc.Sz[0]
            sc.set_hamiltonian(h)
            sc.maxm = 20; sc.nsweeps = 10
            herm = sc.is_hermitian(sc.hamiltonian)
            out.append((herm, np.complex128(sc.gs_energy())))
        ev = np.linalg.eigvals(sc.get_ED_obj().get_hamiltonian().toarray())
        e0 = ev[np.argsort(ev.real)][0]
        print("v=%-6s c=%.0e exact=%s | %s" % (v, c, np.round(e0, 6),
              ["herm=%s E=%s" % (hh, np.round(e, 6)) for hh, e in out]))
```

`reviews/kpm_1/03_first_order_decay.out`:


```
v=python c=5e-03 exact=(-2.457105-0.002134j) | ['herm=True E=(-2.887505+0j)', 'herm=True E=(-2.537763+0j)']
v=python c=2e-02 exact=(-2.457083-0.008536j) | ['herm=False E=(-2.457083-0.008536j)', 'herm=False E=(-2.457083-0.008536j)']
v=3      c=5e-03 exact=(-2.457105-0.002134j) | ['herm=True E=(-2.457105+0j)', 'herm=True E=(-2.457105+0j)']
v=3      c=2e-02 exact=(-2.457083-0.008536j) | ['herm=False E=(-2.457083-0.008536j)', 'herm=False E=(-2.457083-0.008536j)']
v=2      c=5e-03 exact=(-2.457105-0.002134j) | ['herm=True E=(-2.457105+0j)', 'herm=True E=(-2.457105+0j)']
v=2      c=2e-02 exact=(-2.457083-0.008536j) | ['herm=False E=(-2.457083-0.008536j)', 'herm=False E=(-2.457083-0.008536j)']
```

```bash
cd <scratch>/reviews/kpm_1 && <scratch>/run3.sh 04_kpm_natural_units.py 2>&1 | tee 04_kpm_natural_units.out
```

`reviews/kpm_1/04_kpm_natural_units.py`:

```python
# Does the misdispatched Hermitian KPM raise in natural units, or return?
# 03's chain (s=1, field 1.5, 1j*c*Sz0 at c=5e-3, probe True), delta=0.1.
# Reference curves: the same chain at c=0 (Hermitian, DMRG KPM) and the
# c=5e-3 chain on mode="ED" (whose own check calls it non-Hermitian).
# Plus the two ED checks at c=3e-6, where they may disagree.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.algebra import algebra

L = 4
es = np.linspace(-0.5, 4.0, 200)
def build(v, c):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(L):
        h = h + 1.5*sc.Sz[i]
    if c != 0: h = h + 1j*c*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm = 20; sc.nsweeps = 10
    return sc

for v in ("python", 3):
    ref = build(v, 0.0)
    _, y0 = ref.get_dynamical_correlator(name=(ref.Sz[1], ref.Sz[1]), delta=0.1, es=es)
    sc = build(v, 5e-3)
    herm = sc.is_hermitian(sc.hamiltonian)
    try:
        _, y = sc.get_dynamical_correlator(name=(sc.Sz[1], sc.Sz[1]), delta=0.1, es=es)
        print("v=%-6s c=5e-3 probe=%s KPM returned: max|y|=%.4f  max|y - y(c=0)|=%.3e  max|Im y|=%.3e" % (
            v, herm, np.max(np.abs(y)), np.max(np.abs(y - y0)), np.max(np.abs(np.imag(y)))))
    except Exception as e:
        print("v=%-6s c=5e-3 probe=%s KPM raised %s: %s" % (v, herm, type(e).__name__, str(e)[:80]))
    try:
        _, ye = sc.get_dynamical_correlator(name=(sc.Sz[1], sc.Sz[1]), delta=0.1, es=es, mode="ED")
        print("        mode=ED returned: max|y_ED - y_DMRG(c=0)|=%.3e" % np.max(np.abs(ye - y0)))
    except Exception as e:
        print("        mode=ED raised %s: %s" % (type(e).__name__, str(e)[:80]))

sc = build("python", 3e-6)
m = sc.get_ED_obj().get_hamiltonian()
print("ED checks at c=3e-6: algebra.is_hermitian (edtk/dynamics)=%s  algebra.ishermitian (lowest_states)=%s" % (
    algebra.is_hermitian(m), algebra.ishermitian(m)))
```

`reviews/kpm_1/04_kpm_natural_units.out`:


```
v=python c=5e-3 probe=True KPM returned: max|y|=0.4038  max|y - y(c=0)|=6.117e-01  max|Im y|=7.876e-05
        mode=ED raised ValueError: E_max (an upper bound for the spectral radius of the Hamiltonian) must be provid
v=3      c=5e-3 probe=True KPM returned: max|y|=0.6118  max|y - y(c=0)|=1.694e-04  max|Im y|=1.675e-04
        mode=ED raised ValueError: E_max (an upper bound for the spectral radius of the Hamiltonian) must be provid
ED checks at c=3e-6: algebra.is_hermitian (edtk/dynamics)=True  algebra.ishermitian (lowest_states)=False
```


**Reviewer of the `misc` candidate (CONFIRMED, NARROWED)**: `01` and `02` reproduce
to every digit; the probe norm is exactly eps^2, because Sy^2 = 1/4 on spin-1/2,
so the flip at eps=1e-2 is deterministic, and from eps=1e-2 down the new
`disentangle` residual is 0.702 of max|ma| while `eig` gives 1.052e-17 to
1.232e-19. `git show 867e2b4^:` confirms the old branch took `eig` there, the bare
proof being False at every eps. `08` on v3 did not reproduce the hunter's numbers
(4.004e-02 and 3.151e-02 against 3.303e-02 and 9.619e-02), since v3 starts from an
unseeded `randomMPS`; the defect itself reproduces on both. On the Hatano-Nelson
chain at s=1e-2, seeds 1, 2 and 3 give relative errors of 7.391e-02, 7.972e-04 and
5.610e-02 at `nsweeps=10` and 3.597e-02, 1.601e-02 and 2.021e-03 at `nsweeps=30`:
Hermitian DMRG on a non-Hermitian MPO has no fixed answer, so the error is silent
but arbitrary in size. An O(1) XX chain with a loss term i*gamma*Sz2 at gamma=3e-3
is called Hermitian, and `gs_energy()` returns -1.7765282438 + 0i against an exact
-1.7765279985 - 6.42e-05i; at gamma >= 1e-2 NH-DMRG is exact to 1e-15. Hatano-Nelson
at s=1, g=3e-3 is misdispatched too, but there Hermitian DMRG happens to be right
(1.3e-15), part of why nothing caught it. Dividing by the largest raw coefficient
separates cleanly: exactly 1.0 for the hunter's ladder at every eps from 1 to 1e-6,
gamma^2 for the loss ladder, exactly 0.0 for Hermitian operators the proof cannot
see (`1j*Sx0*Sy0`, and Heisenberg plus `i*Sx1*Sy1` plus `i*Sy3*Sz3`) at scale 1 and
1e-6; at eps=1e-7 the raw witness norm is exactly 0, the MPO application's own
floor. Struck, each explicitly:

- The fixed sizes 7.4 per cent (`"python"`) and 9.6 per cent (v3): the size is
  arbitrary, 8.0e-4 to 7.4e-2 across seeds and sweep counts.
- "Any non-Hermitian operator of overall scale below about 1e-2" as the boundary:
  it is the absolute size of the anti-Hermitian part, so weak loss on an O(1)
  Hamiltonian falls through too.
- That both `wf is None -> True` exits "lean the same way" as a consequence:
  nothing executed reaches them.
- The list of other consumers as observed consequences: they read the same
  verdict, but only `disentangle_manifold` and `gs_energy()` were executed.

```bash
cd <scratch>/reviews/misc_1 && <scratch>/run3.sh 01_probe_threshold.py 2>&1 | tee 01_probe_threshold.out
```

`reviews/misc_1/01_probe_threshold.py` is `misc/01_probe_threshold.py` (finding 12) unchanged, rerun; its output, `reviews/misc_1/01_probe_threshold.out`:


```
WARNING, state is not normalizable. Returning None
WARNING, state is not normalizable. Returning None
manifold Gram error 6.59e-15
     eps  proof   probe_norm      probe    resid_new    resid_eig
   1e+00  False    1.000e+00      False    4.066e-15    1.420e-15   (rel. to max|ma|=7.13e-01: 0.000)
   1e-01  False    1.000e-02      False    4.354e-16    9.958e-17   (rel. to max|ma|=7.13e-02: 0.000)
   3e-02  False    9.000e-04      False    1.145e-16    3.065e-17   (rel. to max|ma|=2.14e-02: 0.000)
   1e-02  False    1.000e-04       True    5.000e-03    1.052e-17   (rel. to max|ma|=7.13e-03: 0.702)
   3e-03  False    9.000e-06       True    1.500e-03    3.959e-18   (rel. to max|ma|=2.14e-03: 0.702)
   1e-03  False    1.000e-06       True    5.000e-04    1.026e-18   (rel. to max|ma|=7.13e-04: 0.702)
   1e-04  False    1.000e-08       True    5.000e-05    1.232e-19   (rel. to max|ma|=7.13e-05: 0.702)
eps=3e-02: probe says Hermitian in 0 of 20 seeds
eps=1e-02: probe says Hermitian in 18 of 20 seeds
```

```bash
cd <scratch>/reviews/misc_1 && <scratch>/run3.sh 02_probe_gs_energy.py 2>&1 | tee 02_probe_gs_energy.out
```

`reviews/misc_1/02_probe_gs_energy.py` is `misc/02_probe_gs_energy.py` (finding 12) unchanged, rerun; its output, `reviews/misc_1/02_probe_gs_energy.out`:


```
exact E0(H)/s = -1.3975836830+0.0e+00j, Hermitian part E0/s = -1.7469796037, ratio 0.800000 (sqrt(1-g^2) = 0.800000)
       s  proof    probe    gs_energy()/s      rel.err
   1e+00  False    False    -1.3975836830    4.290e-15
   1e-01  False    False    -1.3975836830    3.655e-15
   1e-02  False     True    -1.5008764226    7.391e-02
   1e-03  False     True    -1.4136759840    1.151e-02
```

```bash
cd <scratch>/reviews/misc_1 && <scratch>/run3.sh 03_probe_gs_energy_v3.py 2>&1 | tee 03_probe_gs_energy_v3.out
```

`reviews/misc_1/03_probe_gs_energy_v3.py` is `misc/08_probe_gs_energy_v3.py` (finding 12) unchanged, rerun; its output, `reviews/misc_1/03_probe_gs_energy_v3.out`:


```
exact E0(H)/s = -1.3975836830+0.0e+00j, Hermitian part E0/s = -1.7469796037, ratio 0.800000 (sqrt(1-g^2) = 0.800000)
       s  proof    probe    gs_energy()/s      rel.err
   1e+00  False    False    -1.3975836830    3.654e-15
   1e-01  False    False    -1.3975836830    3.654e-15
   1e-02  False     True    -1.3416248009    4.004e-02
   1e-03  False     True    -1.3535425511    3.151e-02
```

```bash
cd <scratch>/reviews/misc_1 && <scratch>/run3.sh 04_attacks.py 2>&1 | tee 04_attacks.out
```

`reviews/misc_1/04_attacks.py`:

```python
"""Reviewer attacks on misc_1.

A. What the proof and a rescaled probe say at small eps: the proof on
   A = eps*(Sz0+S+0)+0.3*eps*Sz1, the chain probe on A, and the chain probe on
   A/max|coefficient| (the canonical coefficients), which is the same operator
   up to a positive real factor, so has the same Hermiticity.
B. Is the Hatano-Nelson gs_energy error at s=1e-2 a stable size or seed/sweep
   noise? ("python", seeds 1..3, nsweeps 10 and 30.)
C. Does the absolute threshold also misclassify an O(1) Hamiltonian with a weak
   anti-Hermitian part? Hatano-Nelson at s=1, small g, and an XX chain with a
   small imaginary onsite potential (loss); compare gs_energy against a dense eig.
"""
import warnings
import numpy as np
from dmrgpy import spinchain
from dmrgpy.multioperatortk import canonical

warnings.simplefilter("ignore")

# ---------------- A ----------------
sc = spinchain.Spin_Chain([2]*3, itensor_version="python")
sc.set_hamiltonian(sc.Sz[0]*sc.Sz[1])
print("A. proof / probe / probe on A rescaled to max|coef|=1")
print("%8s %6s %6s %10s %14s" % ("eps", "proof", "probe", "probe_resc", "||dA w||^2"))
for eps in [1.0, 1e-1, 1e-2, 1e-3, 1e-6, 1e-9]:
    A = eps*(sc.Sz[0] + sc.Sx[0] + 1j*sc.Sy[0]) + 0.3*eps*sc.Sz[1]
    cmax = max(abs(c) for c in canonical.canonical_dict(A).values()) \
        if canonical.canonical_dict(A) else 1.0
    np.random.seed(3)
    p = sc.is_hermitian(A)
    np.random.seed(3)
    pr = sc.is_hermitian((1.0/cmax)*A)
    np.random.seed(3)
    w = sc.random_mps()
    d = (A - A.get_dagger())*w
    nd = None if d is None else d.dot(d).real
    print("%8.0e %6s %6s %10s %14s" % (eps, A.is_hermitian(), p, pr,
                                       "None" if nd is None else "%.3e" % nd))

# ---------------- dense helpers ----------------
n = 6
sp = np.array([[0, 1], [0, 0]], dtype=complex); sm = sp.T.copy()
szm = np.diag([0.5, -0.5]).astype(complex)
def site_op(o, i):
    out = np.array([[1.0+0j]])
    for k in range(n):
        out = np.kron(out, o if k == i else np.eye(2))
    return out

def hn_dense(g, s):
    return s*sum((1+g)*site_op(sp, i) @ site_op(sm, i+1) + (1-g)*site_op(sm, i) @ site_op(sp, i+1)
                 for i in range(n-1))/2

def hn_chain(g, s, nsweeps, maxm=20):
    c = spinchain.Spin_Chain([2]*n, itensor_version="python")
    Sp = [c.Sx[i] + 1j*c.Sy[i] for i in range(n)]
    Sm = [c.Sx[i] - 1j*c.Sy[i] for i in range(n)]
    h = 0
    for i in range(n-1):
        h = h + s*((1+g)*Sp[i]*Sm[i+1] + (1-g)*Sm[i]*Sp[i+1])/2
    c.set_hamiltonian(h)
    c.maxm, c.nsweeps = maxm, nsweeps
    return c

def e0_exact(Hd):
    ev = np.linalg.eigvals(Hd)
    return ev[np.argmin(ev.real)]

# ---------------- B ----------------
g, s = 0.6, 1e-2
ex = e0_exact(hn_dense(g, 1.0))
print("\nB. Hatano-Nelson n=6 g=0.6 s=1e-2 ('python'), exact E0/s = %.10f" % ex.real)
print("%5s %8s %6s %16s %10s" % ("seed", "nsweeps", "probe", "gs_energy()/s", "rel.err"))
for nsw in [10, 30]:
    for seed in [1, 2, 3]:
        c = hn_chain(g, s, nsw)
        np.random.seed(seed)
        herm = c.is_hermitian(c.hamiltonian)
        np.random.seed(seed)
        e = c.gs_energy()
        print("%5d %8d %6s %16.10f %10.3e" % (seed, nsw, herm, np.real(e)/s,
                                             abs(np.real(e)/s - ex.real)/abs(ex.real)))

# ---------------- C ----------------
print("\nC. O(1) Hamiltonians with a weak anti-Hermitian part ('python', seed 1, nsweeps 20)")
print("%-28s %6s %6s %22s %22s %10s" % ("model", "proof", "probe", "gs_energy()", "exact E0", "|diff|"))
for gg in [1e-1, 3e-2, 1e-2, 3e-3]:
    c = hn_chain(gg, 1.0, 20)
    np.random.seed(1); herm = c.is_hermitian(c.hamiltonian)
    np.random.seed(1); e = c.gs_energy()
    ex = e0_exact(hn_dense(gg, 1.0))
    print("%-28s %6s %6s %22s %22s %10.3e" % ("Hatano-Nelson g=%.0e" % gg, c.hamiltonian.is_hermitian(),
          herm, "%.10f%+.2ej" % (np.real(e), np.imag(e)), "%.10f%+.2ej" % (ex.real, ex.imag), abs(e-ex)))
for gam in [1e-1, 3e-2, 1e-2, 3e-3]:
    c = spinchain.Spin_Chain([2]*n, itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + c.Sx[i]*c.Sx[i+1] + c.Sy[i]*c.Sy[i+1]
    h = h + 0.3*c.Sz[0] + 1j*gam*c.Sz[2]
    c.set_hamiltonian(h); c.maxm, c.nsweeps = 20, 20
    Hd = sum((site_op(sp, i) @ site_op(sm, i+1) + site_op(sm, i) @ site_op(sp, i+1))/2
             for i in range(n-1)) + 0.3*site_op(szm, 0) + 1j*gam*site_op(szm, 2)
    ex = e0_exact(Hd)
    np.random.seed(1); herm = c.is_hermitian(c.hamiltonian)
    np.random.seed(1); e = c.gs_energy()
    print("%-28s %6s %6s %22s %22s %10.3e" % ("XX + i*%.0e*Sz2 (loss)" % gam, c.hamiltonian.is_hermitian(),
          herm, "%.10f%+.2ej" % (np.real(e), np.imag(e)), "%.10f%+.2ej" % (ex.real, ex.imag), abs(e-ex)))
```

`reviews/misc_1/04_attacks.out`:


```
A. proof / probe / probe on A rescaled to max|coef|=1
     eps  proof  probe probe_resc     ||dA w||^2
   1e+00  False  False      False      1.000e+00
   1e-01  False  False      False      1.000e-02
   1e-02  False   True      False      1.000e-04
   1e-03  False   True      False      1.000e-06
   1e-06  False   True      False      1.000e-12
   1e-09   True   True       True      0.000e+00

B. Hatano-Nelson n=6 g=0.6 s=1e-2 ('python'), exact E0/s = -1.3975836830
 seed  nsweeps  probe    gs_energy()/s    rel.err
    1       10   True    -1.5008764226  7.391e-02
    2       10   True    -1.3964695212  7.972e-04
    3       10   True    -1.4759814037  5.610e-02
    1       30   True    -1.4478517034  3.597e-02
    2       30   True    -1.4199584493  1.601e-02
    3       30   True    -1.3947593602  2.021e-03

C. O(1) Hamiltonians with a weak anti-Hermitian part ('python', seed 1, nsweeps 20)
model                         proof  probe            gs_energy()               exact E0     |diff|
Hatano-Nelson g=1e-01         False  False -1.7382227586-1.08e-16j -1.7382227586-1.27e-19j  9.104e-15
Hatano-Nelson g=3e-02         False  False -1.7461932859-6.26e-17j -1.7461932859+0.00e+00j  4.485e-16
Hatano-Nelson g=1e-02         False  False -1.7468922526-4.99e-17j -1.7468922526-1.91e-19j  6.217e-15
Hatano-Nelson g=3e-03         False   True -1.7469717423+0.00e+00j -1.7469717423-4.92e-30j  1.332e-15
XX + i*1e-01*Sz2 (loss)       False  False -1.7742328250-2.16e-03j -1.7742328250-2.16e-03j  3.114e-15
XX + i*3e-02*Sz2 (loss)       False  False -1.7763236974-6.42e-04j -1.7763236974-6.42e-04j  2.223e-16
XX + i*1e-02*Sz2 (loss)       False  False -1.7765071362-2.14e-04j -1.7765071362-2.14e-04j  1.561e-15
XX + i*3e-03*Sz2 (loss)       False   True -1.7765282438+0.00e+00j -1.7765279985-6.42e-05j  6.419e-05
```

```bash
cd <scratch>/reviews/misc_1 && <scratch>/run3.sh 05_normalized_probe.py 2>&1 | tee 05_normalized_probe.out
```

`reviews/misc_1/05_normalized_probe.py`:

```python
"""Separation measured for the fix: the probe's witness norm ||(A-A^dag)w||^2
divided by cmax^2, cmax = max |coefficient| over the raw terms of A (before
differencing, so no clean_threshold floor), for non-Hermitian operators at any
scale or weakness and for Hermitian operators the canonical proof cannot see."""
import warnings
import numpy as np
from dmrgpy import spinchain

warnings.simplefilter("ignore")
sc = spinchain.Spin_Chain([2]*6, itensor_version="python")
sc.set_hamiltonian(sc.Sz[0]*sc.Sz[1])

def ratio(A, seed=3):
    cmax = max(abs(t[0]) for t in A.op)
    np.random.seed(seed)
    old = sc.maxm; sc.maxm = 8
    try:
        w = sc.random_mps()
        d = (A - A.get_dagger())*w
        nd = d.dot(d).real
    finally:
        sc.maxm = old
    return A.is_hermitian(), nd, cmax, nd/cmax**2

cases = []
for eps in [1.0, 1e-2, 1e-4, 1e-6, 1e-7]:
    cases.append(("NH  eps*(Sz0+S+0)+0.3eps*Sz1, eps=%.0e" % eps,
                  eps*(sc.Sz[0] + sc.Sx[0] + 1j*sc.Sy[0]) + 0.3*eps*sc.Sz[1]))
xx = 0
for i in range(5):
    xx = xx + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
for gam in [1e-2, 3e-3, 1e-4, 1e-6]:
    cases.append(("NH  XX+0.3Sz0+i*%.0e*Sz2" % gam, xx + 0.3*sc.Sz[0] + 1j*gam*sc.Sz[2]))
for s in [1.0, 1e-6]:
    cases.append(("H   %.0e*1j*Sx0*Sy0 (= -Sz0/2)" % s, s*1j*sc.Sx[0]*sc.Sy[0]))
    heis = 0
    for i in range(5):
        heis = heis + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    heis = heis + 1j*sc.Sx[1]*sc.Sy[1] + 1j*sc.Sy[3]*sc.Sz[3]
    cases.append(("H   %.0e*(Heis + i Sx1Sy1 + i Sy3Sz3)" % s, s*heis))
print("%-42s %6s %12s %10s %12s" % ("operator", "proof", "||dA w||^2", "cmax", "ratio"))
for name, A in cases:
    pr, nd, cmax, r = ratio(A)
    print("%-42s %6s %12.3e %10.3e %12.3e" % (name, pr, nd, cmax, r))
```

`reviews/misc_1/05_normalized_probe.out`:


```
operator                                    proof   ||dA w||^2       cmax        ratio
NH  eps*(Sz0+S+0)+0.3eps*Sz1, eps=1e+00     False    1.000e+00  1.000e+00    1.000e+00
NH  eps*(Sz0+S+0)+0.3eps*Sz1, eps=1e-02     False    1.000e-04  1.000e-02    1.000e+00
NH  eps*(Sz0+S+0)+0.3eps*Sz1, eps=1e-04     False    1.000e-08  1.000e-04    1.000e+00
NH  eps*(Sz0+S+0)+0.3eps*Sz1, eps=1e-06     False    1.000e-12  1.000e-06    1.000e+00
NH  eps*(Sz0+S+0)+0.3eps*Sz1, eps=1e-07     False    0.000e+00  1.000e-07    0.000e+00
NH  XX+0.3Sz0+i*1e-02*Sz2                   False    1.000e-04  1.000e+00    1.000e-04
NH  XX+0.3Sz0+i*3e-03*Sz2                   False    9.000e-06  1.000e+00    9.000e-06
NH  XX+0.3Sz0+i*1e-04*Sz2                   False    1.000e-08  1.000e+00    1.000e-08
NH  XX+0.3Sz0+i*1e-06*Sz2                   False    1.000e-12  1.000e+00    1.000e-12
H   1e+00*1j*Sx0*Sy0 (= -Sz0/2)             False    0.000e+00  1.000e+00    0.000e+00
H   1e+00*(Heis + i Sx1Sy1 + i Sy3Sz3)      False    0.000e+00  1.000e+00    0.000e+00
H   1e-06*1j*Sx0*Sy0 (= -Sz0/2)             False    0.000e+00  1.000e-06    0.000e+00
H   1e-06*(Heis + i Sx1Sy1 + i Sy3Sz3)      False    0.000e+00  1.000e-06    0.000e+00
```


**Suggested fix**: make the probe scale-free by rescaling the operator itself before
the probe, which both reviewers measured and which costs no extra application:
take cmax = max|c| over the raw terms of the operator, before differencing and
before `canonical_dict` (so the `clean_threshold` floor is removed too), apply
dA = (op - op^dag)/cmax to the witness, and report Hermitian only when ||dA w||^2
falls below a roundoff-level tolerance of about 1e-20; running
`canonical.is_hermitian` on op/cmax makes the proof scale-free as well. Dividing
the norm afterwards is weaker, since it cannot see what the application has already
floored to zero (eps=1e-7). A denominator from a second application, or from the
coefficients, must keep the numerator a witness norm: a coefficient ratio for dA
would call every unproven Hermitian operator non-Hermitian, since its dA has
nonzero coefficients by construction. Guard 0/0 for an operator that annihilates
the witness and treat it as Hermitian, keep the rule that the fallback is never
turned into a rejection, and correct the comment at `mpsalgebra.py:355-357`.
Independently, `disentangle_manifold` should not trust the verdict when `ma`
contradicts it: if ||ma - ma^dag|| exceeds about 1e-8*||ma||, take `eig`; `ma` is in
hand, so the check is free and would have caught the `867e2b4` half whatever the
probe says. On the ED side one shared relative test should replace both
`algebra.is_zero_matrix`'s 1e-8 on ||h-h^dag||_F^2 and `ishermitian`'s 1e-6 max-abs,
since the two already disagree with each other. Finding 18 lands with this, in the
same file. NUMBERS CHANGE only for non-Hermitian operators whose anti-Hermitian part
is below about 1e-2 in absolute size, and every such number was wrong: the
Hatano-Nelson chain at s=1e-2 moves from -1.5009 (`"python"`) and -1.4437 (v3) to
the exact -1.3976, the weak-loss chain from a real -1.7765282 to -1.7765280 -
6.42e-5i, and `disentangle_manifold` returns `eig`'s eigenvectors again. Hermitian
operators do not move.

### 13. On `"python"`, `to_mpo`'s first truncating sweep runs on the uncanonicalized automaton, so any operator whose coefficients are all at or below about 2e-7 in absolute size becomes a different operator: a lone `1e-7*Sz0` becomes the zero MPO (`vev` exactly 0, the KPM correlator identically 0 with no raise), and a Heisenberg Hamiltonian at scale 3e-7 is built 89 to 95 per cent wrong

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `kpm`

**Where**: `src/dmrgpy/pyitensor/mpobuilder.py:280-292` (`to_mpo`; the two
`result.position(..., cutoff=cutoff)` sweeps at `:290-291`), called with
`pyitensor/chain.py:41` `_BUILD_CUTOFF = 1e-14` through `Chain._mpo`
(`chain.py:486`), i.e. every operator the `"python"` session builds: `vev`,
`apply_operator`, the KPM vertices at `chain.py:1591-1592`, the Hamiltonian
itself.

The automaton is exact, and the two sweeps after it exist only to honour `cutoff`
and `maxdim`. The first sweep runs left to right on the machine as built, before
any canonicalization, so it weighs each channel from the left, where the identity
channel is a string of identities whose weight dominates the first SVD even though
that channel is dead once the right boundary is reached. A relative discarded
weight of 1e-14 is then measured against the identity string, and a small term's
channel falls below it: the discarded fraction for a lone eps*Sz is
eps^2*||Sz||_F^2/||I||_F^2 = eps^2/4, kept at 3e-7, lost at 2e-7, the same at every
L and every site. The moment guard reads a zero operator as zero moments, so KPM
returns an identically zero spectrum with no raise. It came in with `e448699`
("Build the pyitensor MPO as a finite-state machine, not a compression"): the
reference `_sum_of_term_mpos` at the same cutoff keeps both cases exactly. It
survived because `tests/test_mpo_automaton_builder.py`'s term zoo, like every test
operator, has no coefficient below about 1e-3.

**Expected**: vev(eps*A) = eps*vev(A) and C[eps*A, eps*B] = eps^2*C[A,B] exactly, as
v2, v3 and ED give, and as the automaton itself gives before the sweeps.

Repro, from the hunter (`10` to `13` locate the threshold and compare backends, `14`
compares the builders and rescales a whole Hamiltonian):

```bash
cd <scratch>/kpm && <scratch>/run3.sh 10_small_coefficients.py 2>&1 | tee 10_small_coefficients.out
```

`kpm/10_small_coefficients.py`:

```python
# Follow-up of 09 (iii): C[eps*Sz0, eps*Sz3] is not eps^2*C[Sz0,Sz3] at
# eps=1e-9.  Is the returned correlator zero, on which modes, and is it the
# MultiOperator's absolute clean_threshold=1e-8 on coefficients?  Also the
# same threshold on a Hamiltonian term and on vev.
import numpy as np
from dmrgpy import spinchain, multioperator

L = 6
es = np.linspace(-0.5, 5.0, 300)
sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
h = 0
for i in range(L-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h + 0.3*sc.Sz[0])
sc.maxm = 30; sc.nsweeps = 10
print("clean_threshold =", multioperator.clean_threshold)
for eps in (1e-6, 1e-8, 2e-8, 1e-9):
    A = eps*sc.Sz[0]
    print("eps=%.0e  terms of eps*Sz0: %s" % (eps, A.to_terms()))
    for mode in ("DMRG", "ED"):
        x, y = sc.get_dynamical_correlator(mode=mode, name=(eps*sc.Sz[0], eps*sc.Sz[3]), delta=0.1, es=es)
        x, y1 = sc.get_dynamical_correlator(mode=mode, name=(sc.Sz[0], sc.Sz[3]), delta=0.1, es=es)
        print("   %-4s max|C[eps Sz0,eps Sz3]|/eps^2 = %.4f   max|C[Sz0,Sz3]| = %.4f   vev(eps*Sz0)/eps = %.4f  vev(Sz0) = %.4f" % (
            mode, np.max(np.abs(y))/eps**2, np.max(np.abs(y1)),
            np.real(sc.vev(eps*sc.Sz[0], mode=mode))/eps, np.real(sc.vev(sc.Sz[0], mode=mode))))
```

`kpm/10_small_coefficients.out`:


```
clean_threshold = 1e-08
eps=1e-06  terms of eps*Sz0: [((1e-06+0j), [('Sz', 1)])]
   DMRG max|C[eps Sz0,eps Sz3]|/eps^2 = 0.5049   max|C[Sz0,Sz3]| = 0.5049   vev(eps*Sz0)/eps = -0.1894  vev(Sz0) = -0.1894
   ED   max|C[eps Sz0,eps Sz3]|/eps^2 = 0.5052   max|C[Sz0,Sz3]| = 0.5052   vev(eps*Sz0)/eps = -0.1894  vev(Sz0) = -0.1894
eps=1e-08  terms of eps*Sz0: []
   DMRG max|C[eps Sz0,eps Sz3]|/eps^2 = 0.0000   max|C[Sz0,Sz3]| = 0.5049   vev(eps*Sz0)/eps = 0.0000  vev(Sz0) = -0.1894
   ED   max|C[eps Sz0,eps Sz3]|/eps^2 = 0.0000   max|C[Sz0,Sz3]| = 0.5052   vev(eps*Sz0)/eps = 0.0000  vev(Sz0) = -0.1894
eps=2e-08  terms of eps*Sz0: [((2e-08+0j), [('Sz', 1)])]
   DMRG max|C[eps Sz0,eps Sz3]|/eps^2 = 0.0000   max|C[Sz0,Sz3]| = 0.5049   vev(eps*Sz0)/eps = 0.0000  vev(Sz0) = -0.1894
   ED   max|C[eps Sz0,eps Sz3]|/eps^2 = 0.5052   max|C[Sz0,Sz3]| = 0.5052   vev(eps*Sz0)/eps = -0.1894  vev(Sz0) = -0.1894
eps=1e-09  terms of eps*Sz0: []
   DMRG max|C[eps Sz0,eps Sz3]|/eps^2 = 0.0000   max|C[Sz0,Sz3]| = 0.5049   vev(eps*Sz0)/eps = 0.0000  vev(Sz0) = -0.1894
   ED   max|C[eps Sz0,eps Sz3]|/eps^2 = 0.0000   max|C[Sz0,Sz3]| = 0.5052   vev(eps*Sz0)/eps = 0.0000  vev(Sz0) = -0.1894
```

```bash
cd <scratch>/kpm && <scratch>/run3.sh 11_small_coefficients_where.py 2>&1 | tee 11_small_coefficients_where.out
```

`kpm/11_small_coefficients_where.py`:

```python
# Where does the DMRG route lose a 2e-8 coefficient that the MultiOperator
# itself keeps (10: ED right, DMRG zero)?  Walk the operator through each
# step both routes take, then call the session directly.
import numpy as np
from dmrgpy import spinchain, multioperator

L = 6
sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
h = 0
for i in range(L-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h + 0.3*sc.Sz[0])
sc.maxm = 30; sc.nsweeps = 10
wf = sc.get_gs()
for eps in (2e-8, 1e-7, 1e-6):
    A = eps*sc.Sz[0]
    B = multioperator.obj2MO(A, name="vev_multioperator")
    print("eps=%.0e to_terms %s | obj2MO %s | dagger %s | simplify %s" % (
        eps, A.to_terms(), B.to_terms(), A.get_dagger().to_terms(), A.simplify().to_terms()))
    c = sc._session.vev(A.to_terms(), wf.cpp_handle, npow=1)
    c1 = sc._session.vev(sc.Sz[0].to_terms(), wf.cpp_handle, npow=1)
    w2 = A*wf
    print("   session.vev(eps*Sz0)/eps = %.6f  session.vev(Sz0) = %.6f  ||(eps*Sz0)|gs>||/eps = %.6f" % (
        np.real(c)/eps, np.real(c1), np.sqrt(abs(w2.dot(w2)))/eps))
```

`kpm/11_small_coefficients_where.out`:


```
eps=2e-08 to_terms [((2e-08+0j), [('Sz', 1)])] | obj2MO [((2e-08+0j), [('Sz', 1)])] | dagger [((2e-08+0j), [('Sz', 1)])] | simplify [((2e-08+0j), [('Sz', 1)])]
   session.vev(eps*Sz0)/eps = 0.000000  session.vev(Sz0) = -0.189361  ||(eps*Sz0)|gs>||/eps = 0.000000
eps=1e-07 to_terms [((1e-07+0j), [('Sz', 1)])] | obj2MO [((1e-07+0j), [('Sz', 1)])] | dagger [((1e-07+0j), [('Sz', 1)])] | simplify [((1e-07+0j), [('Sz', 1)])]
   session.vev(eps*Sz0)/eps = 0.000000  session.vev(Sz0) = -0.189361  ||(eps*Sz0)|gs>||/eps = 0.000000
eps=1e-06 to_terms [((1e-06+0j), [('Sz', 1)])] | obj2MO [((1e-06+0j), [('Sz', 1)])] | dagger [((1e-06+0j), [('Sz', 1)])] | simplify [((1e-06+0j), [('Sz', 1)])]
   session.vev(eps*Sz0)/eps = -0.189361  session.vev(Sz0) = -0.189361  ||(eps*Sz0)|gs>||/eps = 0.500000
```

```bash
cd <scratch>/kpm && <scratch>/run3.sh 12_small_coefficients_backends.py 2>&1 | tee 12_small_coefficients_backends.out
```

`kpm/12_small_coefficients_backends.py`:

```python
# 11 showed the "python" session returns 0 for vev(eps*Sz0) at eps<=1e-7 while
# the MultiOperator keeps the term.  Locate it (the MPO's own norm, before and
# after to_mpo's truncating sweeps) and check v3/v2 at the same eps, on vev and
# on the KPM correlator C[eps*Sz0, eps*Sz3]/eps^2 against C[Sz0,Sz3].
import numpy as np
from dmrgpy import spinchain
from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor import mpobuilder

L = 6
es = np.linspace(-0.5, 5.0, 300)
for v in ("python", 3, 2):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h + 0.3*sc.Sz[0])
    sc.maxm = 30; sc.nsweeps = 10
    x, yref = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[3]), delta=0.1, es=es)
    ref_vev = np.real(sc.vev(sc.Sz[0]))
    for eps in (1e-5, 1e-6, 3e-7, 1e-7, 3e-8):
        x, y = sc.get_dynamical_correlator(name=(eps*sc.Sz[0], eps*sc.Sz[3]), delta=0.1, es=es)
        y = np.asarray(y)/eps**2
        print("v=%-6s eps=%.0e vev(eps*Sz0)/eps=%.5f (ref %.5f)  max|C/eps^2|=%.4f  max|C/eps^2-C[Sz0,Sz3]|=%.2e" % (
            v, eps, np.real(sc.vev(eps*sc.Sz[0]))/eps, ref_vev, np.max(np.abs(y)), np.max(np.abs(y-yref))))
    if v == "python":
        from dmrgpy.pyitensor.mpsalgebra import inner
        sess = sc._session
        psi = sc.get_gs().cpp_handle
        for eps in (3e-7, 1e-7):
            terms = (eps*sc.Sz[0]).to_terms()
            raw = mpobuilder._automaton_mpo(sess._ampo(terms))
            full = sess._mpo(terms)
            print("   python eps=%.0e <gs|M|gs>/eps: automaton (before to_mpo's sweeps) %.5f, after to_mpo %.5f" % (
                eps, np.real(inner(psi, raw, psi))/eps/np.real(inner(psi, psi)),
                np.real(inner(psi, full, psi))/eps/np.real(inner(psi, psi))))
```

`kpm/12_small_coefficients_backends.out`:


```
v=python eps=1e-05 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=3.44e-15
v=python eps=1e-06 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=2.41e-15
v=python eps=3e-07 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=2.46e-15
v=python eps=1e-07 vev(eps*Sz0)/eps=0.00000 (ref -0.18936)  max|C/eps^2|=0.0000  max|C/eps^2-C[Sz0,Sz3]|=5.05e-01
v=python eps=3e-08 vev(eps*Sz0)/eps=0.00000 (ref -0.18936)  max|C/eps^2|=0.0000  max|C/eps^2-C[Sz0,Sz3]|=5.05e-01
   python eps=3e-07 <gs|M|gs>/eps: automaton (before to_mpo's sweeps) -0.18936, after to_mpo -0.18936
   python eps=1e-07 <gs|M|gs>/eps: automaton (before to_mpo's sweeps) -0.18936, after to_mpo 0.00000
v=3      eps=1e-05 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=4.22e-15
v=3      eps=1e-06 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=3.44e-15
v=3      eps=3e-07 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=2.72e-15
v=3      eps=1e-07 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=2.50e-15
v=3      eps=3e-08 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=3.66e-15
v=2      eps=1e-05 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=4.00e-15
v=2      eps=1e-06 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=2.33e-15
v=2      eps=3e-07 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=2.91e-15
v=2      eps=1e-07 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=1.75e-15
v=2      eps=3e-08 vev(eps*Sz0)/eps=-0.18936 (ref -0.18936)  max|C/eps^2|=0.5049  max|C/eps^2-C[Sz0,Sz3]|=4.72e-15
```

```bash
cd <scratch>/kpm && <scratch>/run3.sh 13_small_coefficients_scaling.py 2>&1 | tee 13_small_coefficients_scaling.out
```

`kpm/13_small_coefficients_scaling.py`:

```python
# How far does the "python" MPO builder's loss of small operators reach?
# (a) vev(eps*Sz_i)/eps against vev(Sz_i), for chain length L and site i;
# (b) the same term inside a sum, vev(Sz_j + eps*Sz_i) - vev(Sz_j) against
#     eps*vev(Sz_i).  v3 as the reference (it kept 3e-8 in 12).
import numpy as np
from dmrgpy import spinchain

for L in (6, 12, 20):
    for v in ("python", 3):
        sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
        h = 0
        for i in range(L-1):
            h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        sc.set_hamiltonian(h + 0.3*sc.Sz[0])
        sc.maxm = 30; sc.nsweeps = 8
        sc.gs_energy()
        out = []
        for site in (0, L//2):
            ref = np.real(sc.vev(sc.Sz[site]))
            row = []
            for eps in (1e-3, 1e-4, 1e-5, 1e-6, 3e-7, 1e-7):
                row.append("%.0e:%.3f" % (eps, np.real(sc.vev(eps*sc.Sz[site]))/eps/ref))
            out.append("site %d ratio vev(eps Sz)/(eps vev(Sz)) %s" % (site, " ".join(row)))
        ref0 = np.real(sc.vev(sc.Sz[0])); base = np.real(sc.vev(sc.Sz[1]))
        row = []
        for eps in (1e-5, 1e-6, 1e-7, 1e-8):
            row.append("%.0e:%.3f" % (eps, (np.real(sc.vev(sc.Sz[1] + eps*sc.Sz[0])) - base)/(eps*ref0)))
        print("L=%2d v=%-6s %s | %s" % (L, v, out[0], out[1]))
        print("              in a sum, (vev(Sz1+eps Sz0)-vev(Sz1))/(eps vev(Sz0)) %s" % " ".join(row))
```

`kpm/13_small_coefficients_scaling.out`:


```
L= 6 v=python site 0 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:-0.000 | site 3 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:0.000
              in a sum, (vev(Sz1+eps Sz0)-vev(Sz1))/(eps vev(Sz0)) 1e-05:1.000 1e-06:1.000 1e-07:-0.000 1e-08:-0.000
L= 6 v=3      site 0 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:1.000 | site 3 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:1.000
              in a sum, (vev(Sz1+eps Sz0)-vev(Sz1))/(eps vev(Sz0)) 1e-05:1.000 1e-06:1.000 1e-07:1.000 1e-08:-0.000
L=12 v=python site 0 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:-0.000 | site 6 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:-0.000
              in a sum, (vev(Sz1+eps Sz0)-vev(Sz1))/(eps vev(Sz0)) 1e-05:1.000 1e-06:1.000 1e-07:-0.000 1e-08:-0.000
L=12 v=3      site 0 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:1.000 | site 6 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:1.000
              in a sum, (vev(Sz1+eps Sz0)-vev(Sz1))/(eps vev(Sz0)) 1e-05:1.000 1e-06:1.000 1e-07:1.000 1e-08:-0.000
L=20 v=python site 0 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:-0.000 | site 10 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:-0.000
              in a sum, (vev(Sz1+eps Sz0)-vev(Sz1))/(eps vev(Sz0)) 1e-05:1.000 1e-06:1.000 1e-07:-0.000 1e-08:-0.000
L=20 v=3      site 0 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:1.000 | site 10 ratio vev(eps Sz)/(eps vev(Sz)) 1e-03:1.000 1e-04:1.000 1e-05:1.000 1e-06:1.000 3e-07:1.000 1e-07:1.000
              in a sum, (vev(Sz1+eps Sz0)-vev(Sz1))/(eps vev(Sz0)) 1e-05:1.000 1e-06:1.000 1e-07:1.000 1e-08:-0.000
```

```bash
cd <scratch>/kpm && <scratch>/run3.sh 14_small_units.py 2>&1 | tee 14_small_units.out
```

`kpm/14_small_units.py`:

```python
# Does the "python" loss of small operators reach a whole Hamiltonian written
# in small energy units (s*H), and the KPM correlator on it?  Anchor: scale
# covariance, E0(s*H) = s*E0(H), and ED.  Also: is it the finite-state-machine
# builder (e448699) or older, i.e. does the reference _sum_of_term_mpos lose
# the lone term too?
import numpy as np
from dmrgpy import spinchain
from dmrgpy.pyitensor import mpobuilder
from dmrgpy.pyitensor.mpsalgebra import inner

L = 6
es = np.linspace(-0.5, 5.0, 300)
for v in ("python", 3):
    for s in (1.0, 1e-6, 3e-7, 1e-7):
        sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
        h = 0
        for i in range(L-1):
            h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        sc.set_hamiltonian(s*h)
        sc.maxm = 30; sc.nsweeps = 8
        try:
            e = sc.gs_energy()/s
        except Exception as ex:
            e = "raised %s" % type(ex).__name__
        eed = sc.gs_energy(mode="ED")/s
        try:
            x, y = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]), delta=0.1*s, es=es*s)
            kp = "max|C|*s=%.4f" % (np.max(np.abs(np.asarray(y)))*s)
        except Exception as ex:
            kp = "KPM raised %s: %s" % (type(ex).__name__, str(ex)[:50])
        print("v=%-6s s=%.0e E0/s DMRG=%s ED=%.6f  %s" % (v, s, e, eed, kp))

sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version="python")
h = sum(sc.Sz[i]*sc.Sz[i+1] for i in range(L-1))
sc.set_hamiltonian(h)
psi = sc.get_gs().cpp_handle
sess = sc._session
for eps in (3e-7, 1e-7):
    ampo = sess._ampo((eps*sc.Sz[0]).to_terms())
    old = mpobuilder._sum_of_term_mpos(ampo, cutoff=1e-14, maxdim=None) if "cutoff" in mpobuilder._sum_of_term_mpos.__code__.co_varnames else mpobuilder._sum_of_term_mpos(ampo)
    new = mpobuilder.to_mpo(ampo, cutoff=1e-14)
    new0 = mpobuilder.to_mpo(ampo, cutoff=0.0)
    n = np.real(inner(psi, psi))
    print("eps=%.0e <Sz0> via old builder %.5f, to_mpo(cutoff=1e-14) %.5f, to_mpo(cutoff=0) %.5f" % (
        eps, np.real(inner(psi, old, psi))/n/eps, np.real(inner(psi, new, psi))/n/eps, np.real(inner(psi, new0, psi))/n/eps))
```

`kpm/14_small_units.out`:


```
v=python s=1e+00 E0/s DMRG=-2.4935771338879213 ED=-2.493577  max|C|*s=0.6259
v=python s=1e-06 E0/s DMRG=-2.493577133887088 ED=-2.493577  max|C|*s=0.6259
v=python s=3e-07 E0/s DMRG=-0.7499999999999989 ED=-2.493577  max|C|*s=1.7814
v=python s=1e-07 E0/s DMRG=-0.7499999999999993 ED=-2.493577  max|C|*s=1.7814
v=3      s=1e+00 E0/s DMRG=-2.4935771338879245 ED=-2.493577  max|C|*s=0.6259
v=3      s=1e-06 E0/s DMRG=-2.4935771326684724 ED=-2.493577  max|C|*s=0.6259
v=3      s=3e-07 E0/s DMRG=-1.2499999956335732 ED=-2.493577  max|C|*s=1.7388
v=3      s=1e-07 E0/s DMRG=-1.249999886820754 ED=-2.493577  max|C|*s=1.7388
eps=3e-07 <Sz0> via old builder 0.50000, to_mpo(cutoff=1e-14) 0.50000, to_mpo(cutoff=0) 0.50000
eps=1e-07 <Sz0> via old builder 0.50000, to_mpo(cutoff=1e-14) 0.00000, to_mpo(cutoff=0) 0.50000
```


**Reviewer (CONFIRMED, NARROWED)**: the headline reproduces through the public API
on a 6-site Heisenberg chain plus 0.3*Sz0 at pinned `maxm=40`, `nsweeps=12`: on
`"python"`, vev(eps*Sz0)/eps is -0.189361 at eps=3e-7 and 0.000000 at 1e-7 and
3e-8, max|C/eps^2| 0.5049 and then 0.0000; on v3 and v2 at the same eps, -0.189361
and 0.5049/0.5051 (within 1.6e-15 to 3.9e-15 of C[Sz0,Sz3]); `mode="ED"` -0.189361.
With the solver taken out, four builds against `AutoMPO.dense_matrix()`, the
Kronecker reference the suite itself uses: (A) a lone eps*Sz on site 1 or 4 of L=6
has relative error about 1e-15 at eps = 1e-6 and 3e-7 and 1.0 at 2e-7, 1e-7, 3e-8
and 1e-10 under `to_mpo(1e-14)`, and 1e-15 or below on every row under `to_mpo(0)`,
the old builder at 1e-14 and a reordered sweep (exact left to right, then the
truncating right to left); (B) Heisenberg scaled by s has relative error 8.9e-01
with bond dimension 3 at s = 3e-7 and 1e-7 on L=6, and 9.5e-01 on L=12, against
1e-15 and bond dimension 5 at s = 1 and 1e-6 and under the other three builds; (C)
Sz2 + eps*Sz1 loses the small term at eps=1e-7 under every build alike, the
relative cutoff doing its job and older than `e448699`; (D) the reordered sweep
keeps the Heisenberg MPO at [4,5,...,5,4] on L=14, with and without a field,
identical to production. Nothing in `CLAUDE.md`, the known-issue files, `ROADMAP.md`
or the earlier records covers it. Struck, each explicitly:

- The in-sum sub-claim, that a 1e-7 term next to a coefficient-1 term is dropped on
  `"python"` and kept on v3: true, but the old builder drops it too, so does the fix,
  and v3 drops it at 1e-8; it is a one-decade threshold difference by design, a note
  at most.
- "Relative to the largest term inside a sum" in the numbers-change statement, for
  the same reason.
- The hunter's "confounded" label on the whole-Hamiltonian row: row (B) shows the
  `"python"` MPO itself is wrong at s <= 3e-7, so its -0.75 is this defect; only the
  v2/v3 small-unit failure is unexplained, and is in "New leads".
- "About 1e-7 or less", narrowed: a lone term is lost at 2e-7 and kept at 3e-7; a
  whole Heisenberg Hamiltonian is lost at 3e-7 and kept at 1e-6.

```bash
cd <scratch>/reviews/kpm_3 && <scratch>/run3.sh 01_public_repro.py 2>&1 | tee 01_public_repro.out
```

`reviews/kpm_3/01_public_repro.py`:

```python
# Reviewer repro of kpm_3 through the public API: vev(eps*Sz0)/eps and the
# KPM correlator C[eps*Sz0, eps*Sz3]/eps^2 on python, v3, v2 and mode="ED",
# on a pinned 6-site Heisenberg chain with a small field (unique GS).
import numpy as np
from dmrgpy import spinchain

L = 6
es = np.linspace(-0.5, 5.0, 300)
for v in ("python", 3, 2):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=v)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h + 0.3*sc.Sz[0])
    sc.maxm = 40; sc.nsweeps = 12
    x, yref = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[3]), delta=0.1, es=es)
    ref = np.real(sc.vev(sc.Sz[0]))
    ref_ed = np.real(sc.vev(sc.Sz[0], mode="ED"))
    for eps in (3e-7, 1e-7, 3e-8):
        x, y = sc.get_dynamical_correlator(name=(eps*sc.Sz[0], eps*sc.Sz[3]), delta=0.1, es=es)
        y = np.asarray(y)/eps**2
        out = "v=%-6s eps=%.0e vev/eps=%.6f (ref %.6f, ED %.6f) max|C/eps^2|=%.4f dev=%.2e" % (
            v, eps, np.real(sc.vev(eps*sc.Sz[0]))/eps, ref, ref_ed,
            np.max(np.abs(y)), np.max(np.abs(y-np.asarray(yref))))
        if v == "python":
            out += "  | ED vev/eps=%.6f" % (np.real(sc.vev(eps*sc.Sz[0], mode="ED"))/eps)
        print(out, flush=True)
```

`reviews/kpm_3/01_public_repro.out`:


```
v=python eps=3e-07 vev/eps=-0.189361 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.5049 dev=4.18e-15  | ED vev/eps=-0.189361
v=python eps=1e-07 vev/eps=0.000000 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.0000 dev=5.05e-01  | ED vev/eps=-0.189361
v=python eps=3e-08 vev/eps=0.000000 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.0000 dev=5.05e-01  | ED vev/eps=-0.189361
v=3      eps=3e-07 vev/eps=-0.189361 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.5049 dev=3.03e-15
v=3      eps=1e-07 vev/eps=-0.189361 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.5049 dev=1.61e-15
v=3      eps=3e-08 vev/eps=-0.189361 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.5049 dev=3.94e-15
v=2      eps=3e-07 vev/eps=-0.189361 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.5051 dev=3.72e-15
v=2      eps=1e-07 vev/eps=-0.189361 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.5051 dev=3.72e-15
v=2      eps=3e-08 vev/eps=-0.189361 (ref -0.189361, ED -0.189361) max|C/eps^2|=0.5051 dev=2.72e-15
```

```bash
cd <scratch>/reviews/kpm_3 && <scratch>/run3.sh 02_dense_builder.py 2>&1 | tee 02_dense_builder.out
```

`reviews/kpm_3/02_dense_builder.py`:

```python
# Solver-free check of pyitensor's to_mpo against AutoMPO.dense_matrix(), the
# independent Kronecker reference, for (A) a lone eps*Sz term, (B) a whole
# Heisenberg Hamiltonian in units s, (C) Sz2 + eps*Sz1 inside a sum.  Four
# builds: production to_mpo(cutoff=1e-14), to_mpo(cutoff=0), the pre-e448699
# reference _sum_of_term_mpos(cutoff=1e-14), and a candidate fix done here
# without editing the repo: an exact left-to-right sweep first, then the
# truncating right-to-left sweep (operator-Schmidt truncation on a canonical
# MPO, so a channel dead from the right has zero weight).
import numpy as np
from dmrgpy.pyitensor import mpobuilder as mb
from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyitensor.tensor import contract_many

SPIN_HALF = 2


def mpo_dense(mpo, sites):
    n = mpo.length()
    T = contract_many([mpo.A(i) for i in range(1, n + 1)])
    si = [sites.si(i) for i in range(1, n + 1)]
    arr = np.asarray(T.transpose_to([i.prime(1) for i in si] + si))
    dim = int(np.prod([i.dim for i in si]))
    return arr.reshape(dim, dim)


def bond_dims(mpo):
    return [mpo.A(i).inds[-1].dim for i in range(1, mpo.length())]


def fixed(ampo, cutoff=1e-14):
    r = mb._automaton_mpo(ampo)
    r.position(r.length(), cutoff=0.0)
    r.position(1, cutoff=cutoff)
    return r


def heis(n, s):
    return [(s, [(op, i), (op, i + 1)]) for i in range(1, n) for op in ("Sx", "Sy", "Sz")]


builds = [("to_mpo(1e-14)", lambda a: mb.to_mpo(a, cutoff=1e-14)),
          ("to_mpo(0)", lambda a: mb.to_mpo(a, cutoff=0.0)),
          ("old(1e-14)", lambda a: mb._sum_of_term_mpos(a, cutoff=1e-14)),
          ("fixed", fixed)]


def report(label, n, terms):
    sites = SiteX([SPIN_HALF]*n)
    a = AutoMPO.from_terms(sites, terms)
    ref = a.dense_matrix()
    nr = np.linalg.norm(ref)
    parts = []
    for name, b in builds:
        m = b(a)
        err = np.linalg.norm(mpo_dense(m, sites) - ref)/nr
        parts.append("%s relerr=%.1e bd=%s" % (name, err, max(bond_dims(m))))
    print("%-34s " % label + " | ".join(parts), flush=True)


L = 6
for site in (1, 4):
    for eps in (1e-6, 3e-7, 2e-7, 1e-7, 3e-8, 1e-10):
        report("A lone eps*Sz%d L=%d eps=%.0e" % (site, L, eps), L, [(eps, [("Sz", site)])])
for n in (6, 12):
    for s in (1.0, 1e-6, 3e-7, 1e-7):
        report("B s*Heis L=%d s=%.0e" % (n, s), n, heis(n, s))
for eps in (1e-6, 1e-7, 1e-8):
    report("C Sz2+eps*Sz1 L=%d eps=%.0e" % (L, eps), L, [(1.0, [("Sz", 2)]), (eps, [("Sz", 1)])])
for nm, tt in (("D Heis L=14", heis(14, 1.0)), ("D Heis+0.3Sz1 L=14", heis(14, 1.0) + [(0.3, [("Sz", 1)])])):
    a = AutoMPO.from_terms(SiteX([SPIN_HALF]*14), tt)
    print(nm, "bond dims (no dense at this size):",
          " | ".join("%s %s" % (b0, bond_dims(b(a))) for b0, b in builds), flush=True)
```

`reviews/kpm_3/02_dense_builder.out`:


```
A lone eps*Sz1 L=6 eps=1e-06       to_mpo(1e-14) relerr=1.2e-15 bd=1 | to_mpo(0) relerr=1.2e-15 bd=1 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=1.2e-15 bd=1
A lone eps*Sz1 L=6 eps=3e-07       to_mpo(1e-14) relerr=8.0e-16 bd=1 | to_mpo(0) relerr=8.2e-16 bd=2 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=8.0e-16 bd=1
A lone eps*Sz1 L=6 eps=2e-07       to_mpo(1e-14) relerr=1.0e+00 bd=1 | to_mpo(0) relerr=1.1e-15 bd=2 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=1.3e-15 bd=1
A lone eps*Sz1 L=6 eps=1e-07       to_mpo(1e-14) relerr=1.0e+00 bd=1 | to_mpo(0) relerr=1.1e-15 bd=2 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=1.3e-15 bd=1
A lone eps*Sz1 L=6 eps=3e-08       to_mpo(1e-14) relerr=1.0e+00 bd=1 | to_mpo(0) relerr=1.2e-15 bd=1 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=1.2e-15 bd=1
A lone eps*Sz1 L=6 eps=1e-10       to_mpo(1e-14) relerr=1.0e+00 bd=1 | to_mpo(0) relerr=1.1e-15 bd=1 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=1.1e-15 bd=1
A lone eps*Sz4 L=6 eps=1e-06       to_mpo(1e-14) relerr=3.0e-16 bd=1 | to_mpo(0) relerr=3.0e-16 bd=1 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=3.0e-16 bd=1
A lone eps*Sz4 L=6 eps=3e-07       to_mpo(1e-14) relerr=3.5e-16 bd=1 | to_mpo(0) relerr=4.3e-16 bd=2 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=3.5e-16 bd=1
A lone eps*Sz4 L=6 eps=2e-07       to_mpo(1e-14) relerr=1.0e+00 bd=1 | to_mpo(0) relerr=3.9e-16 bd=1 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=3.9e-16 bd=1
A lone eps*Sz4 L=6 eps=1e-07       to_mpo(1e-14) relerr=1.0e+00 bd=1 | to_mpo(0) relerr=3.9e-16 bd=1 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=3.9e-16 bd=1
A lone eps*Sz4 L=6 eps=3e-08       to_mpo(1e-14) relerr=1.0e+00 bd=1 | to_mpo(0) relerr=6.2e-16 bd=1 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=6.2e-16 bd=1
A lone eps*Sz4 L=6 eps=1e-10       to_mpo(1e-14) relerr=1.0e+00 bd=1 | to_mpo(0) relerr=7.1e-16 bd=1 | old(1e-14) relerr=0.0e+00 bd=1 | fixed relerr=7.1e-16 bd=1
B s*Heis L=6 s=1e+00               to_mpo(1e-14) relerr=8.4e-16 bd=5 | to_mpo(0) relerr=8.4e-16 bd=5 | old(1e-14) relerr=2.6e-15 bd=5 | fixed relerr=8.4e-16 bd=5
B s*Heis L=6 s=1e-06               to_mpo(1e-14) relerr=9.5e-16 bd=5 | to_mpo(0) relerr=9.5e-16 bd=5 | old(1e-14) relerr=2.9e-15 bd=5 | fixed relerr=9.5e-16 bd=5
B s*Heis L=6 s=3e-07               to_mpo(1e-14) relerr=8.9e-01 bd=3 | to_mpo(0) relerr=1.4e-15 bd=5 | old(1e-14) relerr=5.2e-15 bd=5 | fixed relerr=1.4e-15 bd=5
B s*Heis L=6 s=1e-07               to_mpo(1e-14) relerr=8.9e-01 bd=3 | to_mpo(0) relerr=9.7e-16 bd=5 | old(1e-14) relerr=2.6e-15 bd=5 | fixed relerr=9.7e-16 bd=5
B s*Heis L=12 s=1e+00              to_mpo(1e-14) relerr=1.3e-15 bd=5 | to_mpo(0) relerr=1.3e-15 bd=5 | old(1e-14) relerr=3.7e-15 bd=5 | fixed relerr=1.3e-15 bd=5
B s*Heis L=12 s=1e-06              to_mpo(1e-14) relerr=1.9e-15 bd=5 | to_mpo(0) relerr=1.9e-15 bd=5 | old(1e-14) relerr=3.8e-15 bd=5 | fixed relerr=1.9e-15 bd=5
B s*Heis L=12 s=3e-07              to_mpo(1e-14) relerr=9.5e-01 bd=3 | to_mpo(0) relerr=1.4e-15 bd=5 | old(1e-14) relerr=3.8e-15 bd=5 | fixed relerr=1.4e-15 bd=5
B s*Heis L=12 s=1e-07              to_mpo(1e-14) relerr=9.5e-01 bd=3 | to_mpo(0) relerr=1.3e-15 bd=5 | old(1e-14) relerr=4.2e-15 bd=5 | fixed relerr=1.3e-15 bd=5
C Sz2+eps*Sz1 L=6 eps=1e-06        to_mpo(1e-14) relerr=1.4e-15 bd=2 | to_mpo(0) relerr=1.4e-15 bd=2 | old(1e-14) relerr=3.5e-16 bd=2 | fixed relerr=1.4e-15 bd=2
C Sz2+eps*Sz1 L=6 eps=1e-07        to_mpo(1e-14) relerr=1.0e-07 bd=1 | to_mpo(0) relerr=9.4e-16 bd=2 | old(1e-14) relerr=1.0e-07 bd=1 | fixed relerr=1.0e-07 bd=1
C Sz2+eps*Sz1 L=6 eps=1e-08        to_mpo(1e-14) relerr=1.0e-08 bd=1 | to_mpo(0) relerr=1.3e-15 bd=2 | old(1e-14) relerr=1.0e-08 bd=1 | fixed relerr=1.0e-08 bd=1
D Heis L=14 bond dims (no dense at this size): to_mpo(1e-14) [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4] | to_mpo(0) [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4] | old(1e-14) [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4] | fixed [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4]
D Heis+0.3Sz1 L=14 bond dims (no dense at this size): to_mpo(1e-14) [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4] | to_mpo(0) [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4] | old(1e-14) [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4] | fixed [4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4]
```


**Suggested fix**: run the first sweep exact, then truncate on the way back:
`_automaton_mpo`, then `position(n, cutoff=0)`, then `position(1, cutoff=cutoff,
maxdim=maxdim)`. The exact first sweep left-canonicalizes the machine, so the
truncating return sweep sees true operator-Schmidt values, in which a channel dead
from the right carries zero weight. Measured as the "fixed" build: relative error
1e-15 on every row of (A) and (B), bond dimension unchanged on L=14, and the
in-sum truncation of (C) unchanged, as it should be. The hunter's first option,
contracting the boundary vectors into the edge tensors, misdiagnoses the cause:
they are already contracted, and the truncation sees the identity string's weight
because it measures from the left. Rescaling by 1/max|coef| would cure the lone-term
and uniform-scale cases but treats the symptom. The regression should add a lone
1e-7 term and a 1e-7-scaled Hamiltonian to `TERM_ZOO`, checked against
`dense_matrix()` at a relative tolerance rather than `abs=1e-12`, and the existing
zoo should be rerun against the reordered sweep. NUMBERS CHANGE only on
`itensor_version="python"`, and only for operators whose coefficients are all at
or below about 2e-7 to 3e-7: `vev(1e-7*Sz0)` from 0 to -1.8936e-8, C[1e-7*Sz0,
1e-7*Sz3] from 0 to 1e-14 times the 0.5049-peak curve on the 6-site chain, and the
MPO of a Hamiltonian at scale 3e-7 from 89 to 95 per cent wrong to exact. Whether
`gs_energy()` at that scale then comes out right was not measured: v2 and v3 fail
there too, by a mechanism that is not the MPO builder (in "New leads").

### 14. On v3, `tevol_method="TDVP_GSE"` loses the whole Krylov expansion at the left edge before the first one-site update when site 0 is pinned to one local basis state, which every ladder operator or projector on site 0 produces, so site 0 stays frozen for the whole trajectory and the public TD correlator at the default i=j=0 of a fermion chain is off by 1.7 times its own peak (0.081 against an exact 0.476)

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `realtime`

**Where**: `src/dmrgpy/mpscpp3/chain_session.h:1942` (`Chain::
global_subspace_expand`), `:1497` (`quench_tdvp_gse`, behind `evolution_DC` and
`submode="TD"`), `:1526` (`evolve_and_measure_tdvp_gse`, behind
`evolve_and_measure` and `evolution_ABA`); dispatched from `timedependent.py:190`
and `:290`.

When site 0 of the state being evolved is exactly one local basis vector, v3's
global subspace expansion adds directions (`addBasis` reports a larger maximum link
dimension) and then loses them at the edge bond before the one-site TDVP update
runs, so one-site TDVP can never move site 0. If site 0 is pinned to the first
local basis state (Emp for spinless fermions, Up for S=1/2: `C_0`, `1-N_0`,
`S+_0`, `1/2+Sz_0`), site 0 stays frozen for the whole trajectory and the result
equals `tdvp_gse_sweeps=0` exactly; if it is pinned to the second (Occ, Dn: `N_0`,
`S-_0`, `1/2-Sz_0`), the trajectory is degraded rather than frozen. Bulk-site
starts and the right edge are unaffected. The same failure for a product-state
start is recorded as "a real, isolated bug", intermittent (5 of 8 at 0.02 to
0.04), in `examples/time_evolution/tdvp_gse_VS_ED_time_evolution/main.py:28-41`,
`docs/documentation.md` around line 4449 and `ROADMAP.md` line 82; none says how
wide it is, and here it is deterministic. The origin is `19bbee0` (2026-07-20). It
survived because no test runs v3's `quench_tdvp_gse` against ED at all; the v3
`evolve_and_measure_tdvp_gse` tests start from `Sx_0|GS>`, and `Sx` is invertible,
so it never pins a site; the product-state example mixes an XX+YY coupling into
`h0` to dodge the failure; the fermionic TD tests run the default
`tevol_method="TDVP"`; C(0), and with it every sum rule, is exact, since it is
measured before the first step; and the failure is left/right asymmetric, so a
test on the last site passes.

**Expected**: two references that agree: the exact series from ED matrices by
dense `eigh`, and `"python"`'s `TDVP_GSE` at the same settings, exact to 1e-11 to
1e-14 on every row here, as is v3's own two-site TDVP.

Repro, from the hunter (`08` every integrator on a complex-hopping chain, `09`
`tdvp_gse_sweeps` and `maxm` scans, `10` to `12` which starts fail, with the Schmidt
ranks from the ED vector, `13` the public TD spectrum):

```bash
cd <scratch>/realtime && <scratch>/run3.sh 08_quench_integrators_complex.py 2>&1 | tee 08_quench_integrators_complex.out
```

`realtime/08_quench_integrators_complex.py`:

```python
"""Lens realtime, lead 8: the TD correlator series evolution_DC(mode="DMRG")
(the conjugation 867e2b4 kept) on every integrator, with a COMPLEX
Hamiltonian and a non-Hermitian, fermionic pair, against a hand-built
exact series.

Same chain as 02 (4 spinless fermions, hopping -t e^{i phi}, 0.3 N_0).
Pair name=(A,B) = (Cdag_0, C_2): the house series is
  C(t) = <GS|A e^{+i(H-E0)t} B|GS> = sum_n <GS|A|n><n|B|GS> e^{+i D_n t},
what evolution_DC documents for both ED and DMRG.
"""
import numpy as np
from dmrgpy import fermionchain, cppext, timedependent

n, hop, phi, v0 = 4, 1.0, 0.4, 0.3
nt, dt = 40, 0.05
ts = dt*np.arange(nt)

def build(version, method):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    fc.maxm, fc.nsweeps = 20, 20
    fc.tevol_method = method
    h = 0
    for j in range(n - 1):
        h = h + (-hop*np.exp(1j*phi))*fc.Cdag[j+1]*fc.C[j]
        h = h + (-hop*np.exp(-1j*phi))*fc.Cdag[j]*fc.C[j+1]
    h = h + v0*fc.N[0]
    fc.set_hamiltonian(h)
    return fc

fc = build("python", "TDVP")
ed = fc.get_ED_obj()
H = ed.get_operator(fc.hamiltonian).toarray()
A = ed.get_operator(fc.Cdag[0]).toarray()
B = ed.get_operator(fc.C[2]).toarray()
w, v = np.linalg.eigh(H)
g = v[:, 0]
a_n = np.conj(v).T@(A.conj().T@g)      # <n|A^dag|GS> = conj(<GS|A|n>)
b_n = np.conj(v).T@(B@g)                # <n|B|GS>
M = np.conj(a_n)*b_n
print("max|Im M_n| = %.3f (complex weights)" % np.max(np.abs(M.imag)))
ref = np.array([np.sum(M*np.exp(1j*(w - w[0])*t)) for t in ts])
bwd = np.array([np.sum(M*np.exp(-1j*(w - w[0])*t)) for t in ts])
print("anchor: max|fwd-bwd| = %.3f, max|fwd-conj(fwd)| = %.3f"
      % (np.max(np.abs(ref-bwd)), np.max(np.abs(ref-np.conj(ref)))))

rows = [("ED", "python", "TDVP")]
rows += [("DMRG", "python", m) for m in ("TDVP", "TDVP_GSE", "TEBD", "AUTO", "MPO")]
if cppext.available(3):
    rows += [("DMRG", 3, m) for m in ("TDVP", "TDVP_GSE", "TEBD", "AUTO", "MPO")]
if cppext.available(2):
    rows += [("DMRG", 2, m) for m in ("MPO",)]
for mode, version, method in rows:
    fc = build(version, method)
    try:
        _t, cs = timedependent.evolution_DC(fc, mode=mode,
                     name=(fc.Cdag[0], fc.C[2]), nt=nt, dt=dt)
    except Exception as exc:
        print("%-4s %-7s %-9s raised %s: %s" % (mode, version, method,
              type(exc).__name__, str(exc)[:100])); continue
    cs = np.asarray(cs)
    print("%-4s %-7s %-9s  max|cs-ref|=%.2e  max|cs-bwd|=%.2e  "
          "max|cs-conj(ref)|=%.2e" % (mode, version, method,
          np.max(np.abs(cs-ref)), np.max(np.abs(cs-bwd)),
          np.max(np.abs(cs-np.conj(ref)))))
```

`realtime/08_quench_integrators_complex.out`:


```



max|Im M_n| = 0.142 (complex weights)
anchor: max|fwd-bwd| = 0.279, max|fwd-conj(fwd)| = 0.613
ED   python  TDVP       max|cs-ref|=5.18e-10  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG python  TDVP       max|cs-ref|=2.76e-14  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG python  TDVP_GSE   max|cs-ref|=1.00e-14  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG python  TEBD       max|cs-ref|=5.88e-05  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG python  AUTO       max|cs-ref|=5.88e-05  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG python  MPO        max|cs-ref|=2.25e-04  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG 3       TDVP       max|cs-ref|=2.45e-13  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG 3       TDVP_GSE   max|cs-ref|=2.30e-01  max|cs-bwd|=3.55e-01  max|cs-conj(ref)|=4.69e-01
DMRG 3       TEBD       max|cs-ref|=5.88e-05  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG 3       AUTO       max|cs-ref|=5.88e-05  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG 3       MPO        max|cs-ref|=2.25e-04  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
DMRG 2       MPO        max|cs-ref|=2.25e-04  max|cs-bwd|=2.79e-01  max|cs-conj(ref)|=6.13e-01
```

```bash
cd <scratch>/realtime && <scratch>/run3.sh 09_v3_quench_gse_discriminate.py 2>&1 | tee 09_v3_quench_gse_discriminate.out
```

`realtime/09_v3_quench_gse_discriminate.py`:

```python
"""Lens realtime, 09: what makes itensor_version=3 evolution_DC under
tevol_method="TDVP_GSE" (Chain::quench_tdvp_gse) miss the exact TD series
by 0.23 in 08, where "python"'s TDVP_GSE is at 1e-14 and v3's own
evolve_and_measure_tdvp_gse (02) at 6e-13?

Every row: v3 and "python" evolution_DC(mode="DMRG") under TDVP_GSE
against the exact Lehmann series sum_n <GS|A|n><n|B|GS> e^{+i D_n t}
built from ED's matrices by dense eigh (as in 08). Varies the chain
(complex / real hopping, spin chain), the pair (non-Hermitian fermionic,
Hermitian N_0,N_0), tdvp_gse_sweeps, and maxm.
"""
import numpy as np, functools
print = functools.partial(print, flush=True)
from dmrgpy import fermionchain, spinchain, timedependent

nt, dt = 40, 0.05
ts = dt*np.arange(nt)

def ferm(version, phi, n=4):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    h = 0
    for j in range(n - 1):
        h = h + (-np.exp(1j*phi))*fc.Cdag[j+1]*fc.C[j]
        h = h + (-np.exp(-1j*phi))*fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0]
    fc.set_hamiltonian(h)
    return fc

def spin(version, n=6):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc

def exact(ch, A, B):
    ed = ch.get_ED_obj()
    H = ed.get_operator(ch.hamiltonian).toarray()
    Am = ed.get_operator(A).toarray(); Bm = ed.get_operator(B).toarray()
    w, v = np.linalg.eigh(H); g = v[:, 0]
    M = np.conj(np.conj(v).T@(Am.conj().T@g))*(np.conj(v).T@(Bm@g))
    return np.array([np.sum(M*np.exp(1j*(w - w[0])*t)) for t in ts])

cases = [
    ("fermion phi=0.4, (Cdag0,C2)", lambda v: ferm(v, 0.4), lambda c: (c.Cdag[0], c.C[2])),
    ("fermion phi=0,   (Cdag0,C2)", lambda v: ferm(v, 0.0), lambda c: (c.Cdag[0], c.C[2])),
    ("fermion phi=0.4, (N0,N0)",    lambda v: ferm(v, 0.4), lambda c: (c.N[0], c.N[0])),
    ("fermion phi=0.4, (Cdag0,C0)", lambda v: ferm(v, 0.4), lambda c: (c.Cdag[0], c.C[0])),
    ("spin Heisenberg, (Sz0,Sz0)",  spin,                  lambda c: (c.Sz[0], c.Sz[0])),
    ("spin Heisenberg, (Sx0,Sx2)",  spin,                  lambda c: (c.Sx[0], c.Sx[2])),
]
for label, mk, pair in cases:
    ref = None
    for version in (3, "python"):
        for sweeps, maxm in ((3, 20), (0, 20), (nt, 20), (3, 64)):
            ch = mk(version)
            ch.maxm, ch.nsweeps = maxm, 20
            ch.tevol_method = "TDVP_GSE"
            ch.tdvp_gse_sweeps = sweeps
            A, B = pair(ch)
            if ref is None:
                ref = exact(ch, A, B)
            _t, cs = timedependent.evolution_DC(ch, mode="DMRG", name=(A, B),
                                                nt=nt, dt=dt)
            cs = np.asarray(cs)
            err = np.abs(cs - ref)
            print("%-30s v=%-6s gse_sweeps=%2d maxm=%2d  max|cs-exact|=%.2e"
                  " (scale %.3f)  |cs0-ref0|=%.1e  err at t=%.2f: %.2e"
                  % (label, version, sweeps, maxm, err.max(),
                     np.abs(ref).max(), err[0], ts[5], err[5]))
            if version == "python" and (sweeps, maxm) != (3, 20):
                break
```

`realtime/09_v3_quench_gse_discriminate.out`:


```



fermion phi=0.4, (Cdag0,C2)    v=3      gse_sweeps= 3 maxm=20  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=6.7e-13  err at t=0.25: 2.13e-02
fermion phi=0.4, (Cdag0,C2)    v=3      gse_sweeps= 0 maxm=20  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=5.2e-14  err at t=0.25: 2.13e-02








































fermion phi=0.4, (Cdag0,C2)    v=3      gse_sweeps=40 maxm=20  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=1.0e-12  err at t=0.25: 2.13e-02



fermion phi=0.4, (Cdag0,C2)    v=3      gse_sweeps= 3 maxm=64  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=3.9e-13  err at t=0.25: 2.13e-02
fermion phi=0.4, (Cdag0,C2)    v=python gse_sweeps= 3 maxm=20  max|cs-exact|=3.36e-14 (scale 0.324)  |cs0-ref0|=3.4e-14  err at t=0.25: 3.32e-14
fermion phi=0.4, (Cdag0,C2)    v=python gse_sweeps= 0 maxm=20  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=1.3e-14  err at t=0.25: 2.13e-02



fermion phi=0,   (Cdag0,C2)    v=3      gse_sweeps= 3 maxm=20  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=5.1e-13  err at t=0.25: 2.13e-02
fermion phi=0,   (Cdag0,C2)    v=3      gse_sweeps= 0 maxm=20  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=3.1e-14  err at t=0.25: 2.13e-02








































fermion phi=0,   (Cdag0,C2)    v=3      gse_sweeps=40 maxm=20  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=3.6e-13  err at t=0.25: 2.13e-02



fermion phi=0,   (Cdag0,C2)    v=3      gse_sweeps= 3 maxm=64  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=1.6e-15  err at t=0.25: 2.13e-02
fermion phi=0,   (Cdag0,C2)    v=python gse_sweeps= 3 maxm=20  max|cs-exact|=8.00e-14 (scale 0.324)  |cs0-ref0|=8.0e-14  err at t=0.25: 7.96e-14
fermion phi=0,   (Cdag0,C2)    v=python gse_sweeps= 0 maxm=20  max|cs-exact|=2.30e-01 (scale 0.324)  |cs0-ref0|=1.3e-14  err at t=0.25: 2.13e-02



fermion phi=0.4, (N0,N0)       v=3      gse_sweeps= 3 maxm=20  max|cs-exact|=1.19e-11 (scale 0.408)  |cs0-ref0|=4.7e-13  err at t=0.25: 4.49e-13
fermion phi=0.4, (N0,N0)       v=3      gse_sweeps= 0 maxm=20  max|cs-exact|=3.95e-01 (scale 0.408)  |cs0-ref0|=8.4e-13  err at t=0.25: 1.15e-02








































fermion phi=0.4, (N0,N0)       v=3      gse_sweeps=40 maxm=20  max|cs-exact|=1.27e-11 (scale 0.408)  |cs0-ref0|=4.0e-13  err at t=0.25: 4.07e-13



fermion phi=0.4, (N0,N0)       v=3      gse_sweeps= 3 maxm=64  max|cs-exact|=1.18e-11 (scale 0.408)  |cs0-ref0|=6.3e-13  err at t=0.25: 5.52e-13
fermion phi=0.4, (N0,N0)       v=python gse_sweeps= 3 maxm=20  max|cs-exact|=2.95e-13 (scale 0.408)  |cs0-ref0|=3.2e-14  err at t=0.25: 3.26e-14
fermion phi=0.4, (N0,N0)       v=python gse_sweeps= 0 maxm=20  max|cs-exact|=3.95e-01 (scale 0.408)  |cs0-ref0|=4.7e-14  err at t=0.25: 1.15e-02



fermion phi=0.4, (Cdag0,C0)    v=3      gse_sweeps= 3 maxm=20  max|cs-exact|=8.44e-02 (scale 0.408)  |cs0-ref0|=2.7e-13  err at t=0.25: 1.10e-03
fermion phi=0.4, (Cdag0,C0)    v=3      gse_sweeps= 0 maxm=20  max|cs-exact|=8.44e-02 (scale 0.408)  |cs0-ref0|=6.6e-14  err at t=0.25: 1.10e-03








































fermion phi=0.4, (Cdag0,C0)    v=3      gse_sweeps=40 maxm=20  max|cs-exact|=8.44e-02 (scale 0.408)  |cs0-ref0|=2.4e-12  err at t=0.25: 1.10e-03



fermion phi=0.4, (Cdag0,C0)    v=3      gse_sweeps= 3 maxm=64  max|cs-exact|=8.44e-02 (scale 0.408)  |cs0-ref0|=2.6e-13  err at t=0.25: 1.10e-03
fermion phi=0.4, (Cdag0,C0)    v=python gse_sweeps= 3 maxm=20  max|cs-exact|=2.59e-14 (scale 0.408)  |cs0-ref0|=2.6e-14  err at t=0.25: 2.57e-14
fermion phi=0.4, (Cdag0,C0)    v=python gse_sweeps= 0 maxm=20  max|cs-exact|=8.44e-02 (scale 0.408)  |cs0-ref0|=2.3e-14  err at t=0.25: 1.10e-03



spin Heisenberg, (Sz0,Sz0)     v=3      gse_sweeps= 3 maxm=20  max|cs-exact|=2.22e-11 (scale 0.250)  |cs0-ref0|=1.1e-16  err at t=0.25: 5.00e-14
spin Heisenberg, (Sz0,Sz0)     v=3      gse_sweeps= 0 maxm=20  max|cs-exact|=2.32e-11 (scale 0.250)  |cs0-ref0|=8.9e-16  err at t=0.25: 1.09e-13








































spin Heisenberg, (Sz0,Sz0)     v=3      gse_sweeps=40 maxm=20  max|cs-exact|=2.15e-11 (scale 0.250)  |cs0-ref0|=1.1e-15  err at t=0.25: 1.30e-13



spin Heisenberg, (Sz0,Sz0)     v=3      gse_sweeps= 3 maxm=64  max|cs-exact|=2.16e-11 (scale 0.250)  |cs0-ref0|=5.6e-16  err at t=0.25: 1.19e-13
spin Heisenberg, (Sz0,Sz0)     v=python gse_sweeps= 3 maxm=20  max|cs-exact|=1.54e-12 (scale 0.250)  |cs0-ref0|=3.9e-16  err at t=0.25: 2.43e-15
spin Heisenberg, (Sz0,Sz0)     v=python gse_sweeps= 0 maxm=20  max|cs-exact|=6.42e-15 (scale 0.250)  |cs0-ref0|=7.5e-16  err at t=0.25: 8.21e-16



spin Heisenberg, (Sx0,Sx2)     v=3      gse_sweeps= 3 maxm=20  max|cs-exact|=3.04e-11 (scale 0.116)  |cs0-ref0|=4.6e-13  err at t=0.25: 2.61e-13
spin Heisenberg, (Sx0,Sx2)     v=3      gse_sweeps= 0 maxm=20  max|cs-exact|=3.28e-11 (scale 0.116)  |cs0-ref0|=1.3e-13  err at t=0.25: 5.05e-13








































spin Heisenberg, (Sx0,Sx2)     v=3      gse_sweeps=40 maxm=20  max|cs-exact|=3.25e-11 (scale 0.116)  |cs0-ref0|=3.2e-13  err at t=0.25: 9.85e-13



spin Heisenberg, (Sx0,Sx2)     v=3      gse_sweeps= 3 maxm=64  max|cs-exact|=3.21e-11 (scale 0.116)  |cs0-ref0|=2.4e-13  err at t=0.25: 4.86e-13
spin Heisenberg, (Sx0,Sx2)     v=python gse_sweeps= 3 maxm=20  max|cs-exact|=2.87e-12 (scale 0.116)  |cs0-ref0|=4.6e-16  err at t=0.25: 5.15e-13
spin Heisenberg, (Sx0,Sx2)     v=python gse_sweeps= 0 maxm=20  max|cs-exact|=1.32e-13 (scale 0.116)  |cs0-ref0|=1.1e-15  err at t=0.25: 1.03e-14
```

```bash
cd <scratch>/realtime && <scratch>/run3.sh 10_v3_gse_parity_odd.py 2>&1 | tee 10_v3_gse_parity_odd.out
```

`realtime/10_v3_gse_parity_odd.py`:

```python
"""Lens realtime, 10: is itensor_version=3's global subspace expansion a
no-op on a fermion-parity-odd state, on both GSE routes?

09 showed Chain::quench_tdvp_gse on a fermion-odd pair returning, bit for
bit in its error, what "python" returns at tdvp_gse_sweeps=0 (no expansion
at all), while a parity-even pair (N0,N0) and every spin pair expand fine.
Here, on a 6-site chain (bond dimension not full, so one-site TDVP without
expansion is visibly wrong), both v3 GSE routes:
  evolution_DC(name=(Cdag_0, C_2))           -> Chain::quench_tdvp_gse
  evolution_ABA(A=C_0 or Cdag_0, B=N_1)       -> Chain::evolve_and_measure_tdvp_gse
at tdvp_gse_sweeps = 0 and 3, against exact (ED matrices, dense eigh,
e^{-iHt} by hand), with "python" alongside.
"""
import numpy as np, functools
import scipy.linalg as sla
print = functools.partial(print, flush=True)
from dmrgpy import fermionchain, timedependent

n, nt, dt = 6, 40, 0.05
ts = dt*np.arange(nt)

def ferm(version, sweeps):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    fc.maxm, fc.nsweeps = 40, 20
    fc.tevol_method = "TDVP_GSE"; fc.tdvp_gse_sweeps = sweeps
    h = 0
    for j in range(n - 1):
        h = h + (-np.exp(0.4j))*fc.Cdag[j+1]*fc.C[j]
        h = h + (-np.exp(-0.4j))*fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0] + 0.5*fc.N[2]*fc.N[3]
    fc.set_hamiltonian(h)
    return fc

fc = ferm("python", 3)
ed = fc.get_ED_obj()
H = ed.get_operator(fc.hamiltonian).toarray()
w, v = np.linalg.eigh(H); g = v[:, 0]
print("E0=%.6f gap=%.4f" % (w[0], w[1]-w[0]))
mat = lambda o: ed.get_operator(o).toarray()
def dc_ref(A, B):
    M = np.conj(np.conj(v).T@(mat(A).conj().T@g))*(np.conj(v).T@(mat(B)@g))
    return np.array([np.sum(M*np.exp(1j*(w - w[0])*t)) for t in ts])
def aba_ref(A, B):
    psi0 = mat(A)@g; Bm = mat(B); c = np.conj(v).T@psi0
    return np.array([np.vdot(v@(np.exp(-1j*w*t)*c), Bm@(v@(np.exp(-1j*w*t)*c)))
                     for t in ts])

for version in (3, "python"):
    for sweeps in (0, 3):
        ch = ferm(version, sweeps)
        _t, cs = timedependent.evolution_DC(ch, mode="DMRG",
                        name=(ch.Cdag[0], ch.C[2]), nt=nt, dt=dt)
        ref = dc_ref(ch.Cdag[0], ch.C[2])
        print("v=%-6s sweeps=%d evolution_DC (Cdag0,C2): max|err|=%.2e (scale %.3f)"
              % (version, sweeps, np.max(np.abs(np.asarray(cs)-ref)), np.abs(ref).max()))
        ch = ferm(version, sweeps)
        _t, cs = timedependent.evolution_DC(ch, mode="DMRG",
                        name=(ch.N[0], ch.N[0]), nt=nt, dt=dt)
        ref = dc_ref(ch.N[0], ch.N[0])
        print("v=%-6s sweeps=%d evolution_DC (N0,N0):    max|err|=%.2e (scale %.3f)"
              % (version, sweeps, np.max(np.abs(np.asarray(cs)-ref)), np.abs(ref).max()))
        for Aname in ("C", "Cdag"):
            ch = ferm(version, sweeps)
            A = ch.C[0] if Aname == "C" else ch.Cdag[0]
            _t, y = timedependent.evolution_ABA(ch, A=A, B=ch.N[1],
                                                 mode="DMRG", nt=nt, dt=dt)
            ref = aba_ref(A, ch.N[1])
            print("v=%-6s sweeps=%d evolution_ABA A=%-4s_0 B=N_1: max|err|=%.2e (scale %.3f)"
                  % (version, sweeps, Aname, np.max(np.abs(np.asarray(y)-ref)),
                     np.abs(ref).max()))
```

`realtime/10_v3_gse_parity_odd.out`:


```
E0=-3.313785 gap=0.3511
v=3      sweeps=0 evolution_DC (Cdag0,C2): max|err|=2.01e-01 (scale 0.307)
v=3      sweeps=0 evolution_DC (N0,N0):    max|err|=3.97e-01 (scale 0.418)
v=3      sweeps=0 evolution_ABA A=C   _0 B=N_1: max|err|=2.87e-02 (scale 0.171)
v=3      sweeps=0 evolution_ABA A=Cdag_0 B=N_1: max|err|=4.00e-02 (scale 0.519)



v=3      sweeps=3 evolution_DC (Cdag0,C2): max|err|=2.01e-01 (scale 0.307)



v=3      sweeps=3 evolution_DC (N0,N0):    max|err|=2.30e-07 (scale 0.418)



v=3      sweeps=3 evolution_ABA A=C   _0 B=N_1: max|err|=2.87e-02 (scale 0.171)



v=3      sweeps=3 evolution_ABA A=Cdag_0 B=N_1: max|err|=5.54e-11 (scale 0.519)
v=python sweeps=0 evolution_DC (Cdag0,C2): max|err|=2.01e-01 (scale 0.307)
v=python sweeps=0 evolution_DC (N0,N0):    max|err|=3.97e-01 (scale 0.418)
v=python sweeps=0 evolution_ABA A=C   _0 B=N_1: max|err|=2.87e-02 (scale 0.171)
v=python sweeps=0 evolution_ABA A=Cdag_0 B=N_1: max|err|=4.00e-02 (scale 0.519)
v=python sweeps=3 evolution_DC (Cdag0,C2): max|err|=1.52e-14 (scale 0.307)
v=python sweeps=3 evolution_DC (N0,N0):    max|err|=2.30e-07 (scale 0.418)
v=python sweeps=3 evolution_ABA A=C   _0 B=N_1: max|err|=6.84e-14 (scale 0.171)
v=python sweeps=3 evolution_ABA A=Cdag_0 B=N_1: max|err|=1.07e-13 (scale 0.519)
```

```bash
cd <scratch>/realtime && <scratch>/run3.sh 11_v3_gse_which_states.py 2>&1 | tee 11_v3_gse_which_states.out
```

`realtime/11_v3_gse_which_states.py`:

```python
"""Lens realtime, 11: which states does itensor_version=3's global subspace
expansion fail to expand? 10 showed C_0|GS> not expanded while
Cdag_0|GS> is, on both v3 GSE routes. Discriminate: filling (mu), the
spin analog (S-_0 vs S+_0), gse_cutoff and krylov order.

Observable: evolution_ABA(A, B) under tevol_method="TDVP_GSE" against
exact (ED matrices, dense eigh, e^{-iHt} by hand), and the bond dimensions
of the state after the run (return_wf) to see whether it grew.
"""
import numpy as np, functools
print = functools.partial(print, flush=True)
from dmrgpy import fermionchain, spinchain, timedependent

n, nt, dt = 6, 40, 0.05
ts = dt*np.arange(nt)

def ferm(version, mu, sweeps=3, cutoff=1e-8, kord=3):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    fc.maxm, fc.nsweeps = 40, 20
    fc.tevol_method = "TDVP_GSE"; fc.tdvp_gse_sweeps = sweeps
    fc.tdvp_gse_cutoff = cutoff; fc.tdvp_gse_krylov_order = kord
    h = 0
    for j in range(n - 1):
        h = h - fc.Cdag[j+1]*fc.C[j] - fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0] + 0.5*fc.N[2]*fc.N[3] - mu*sum(fc.N[j] for j in range(n))
    fc.set_hamiltonian(h)
    return fc

def spin(version, sweeps=3):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    sc.maxm, sc.nsweeps = 40, 20
    sc.tevol_method = "TDVP_GSE"; sc.tdvp_gse_sweeps = sweeps
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc

def aba_err(ch, A, B):
    ed = ch.get_ED_obj()
    H = ed.get_operator(ch.hamiltonian).toarray()
    w, v = np.linalg.eigh(H); g = v[:, 0]
    psi0 = ed.get_operator(A).toarray()@g; Bm = ed.get_operator(B).toarray()
    c = np.conj(v).T@psi0
    ref = np.array([np.vdot(v@(np.exp(-1j*w*t)*c), Bm@(v@(np.exp(-1j*w*t)*c)))
                    for t in ts])
    _t, y = timedependent.evolution_ABA(ch, A=A, B=B, mode="DMRG", nt=nt, dt=dt)
    return np.max(np.abs(np.asarray(y) - ref)), np.abs(ref).max()

for mu in (0.0, 0.8, -0.8):
    for Aname in ("C", "Cdag"):
        for version in (3, "python"):
            ch = ferm(version, mu)
            ng = ch.vev(sum(ch.N[j] for j in range(n))).real
            A = ch.C[0] if Aname == "C" else ch.Cdag[0]
            e, s = aba_err(ch, A, ch.N[1])
            print("fermion mu=%+.1f <N>=%.3f A=%-4s_0 v=%-6s max|err|=%.2e (scale %.3f)"
                  % (mu, ng, Aname, version, e, s))
for Aname in ("Sm", "Sp"):
    for version in (3, "python"):
        ch = spin(version)
        A = (ch.Sx[0] - 1j*ch.Sy[0]) if Aname == "Sm" else (ch.Sx[0] + 1j*ch.Sy[0])
        e, s = aba_err(ch, A, ch.Sz[1])
        print("spin A=%s_0 v=%-6s max|err|=%.2e (scale %.3f)" % (Aname, version, e, s))
for cutoff, kord in ((1e-8, 3), (1e-14, 3), (0.0, 3), (1e-8, 5)):
    ch = ferm(3, 0.0, cutoff=cutoff, kord=kord)
    e, s = aba_err(ch, ch.C[0], ch.N[1])
    print("fermion mu=0 A=C_0 v=3 gse_cutoff=%g krylov=%d max|err|=%.2e (scale %.3f)"
          % (cutoff, kord, e, s))
```

`realtime/11_v3_gse_which_states.out`:


```



fermion mu=+0.0 <N>=3.000 A=C   _0 v=3      max|err|=2.87e-02 (scale 0.171)
fermion mu=+0.0 <N>=3.000 A=C   _0 v=python max|err|=5.41e-14 (scale 0.171)



fermion mu=+0.0 <N>=3.000 A=Cdag_0 v=3      max|err|=4.51e-03 (scale 0.519)
fermion mu=+0.0 <N>=3.000 A=Cdag_0 v=python max|err|=1.47e-13 (scale 0.519)



fermion mu=+0.8 <N>=4.000 A=C   _0 v=3      max|err|=3.00e-01 (scale 0.314)
fermion mu=+0.8 <N>=4.000 A=C   _0 v=python max|err|=4.06e-07 (scale 0.314)



fermion mu=+0.8 <N>=4.000 A=Cdag_0 v=3      max|err|=8.28e-04 (scale 0.284)
fermion mu=+0.8 <N>=4.000 A=Cdag_0 v=python max|err|=3.63e-10 (scale 0.284)



fermion mu=-0.8 <N>=2.000 A=C   _0 v=3      max|err|=4.09e-03 (scale 0.020)
fermion mu=-0.8 <N>=2.000 A=C   _0 v=python max|err|=9.00e-13 (scale 0.020)



fermion mu=-0.8 <N>=2.000 A=Cdag_0 v=3      max|err|=1.48e-02 (scale 0.618)
fermion mu=-0.8 <N>=2.000 A=Cdag_0 v=python max|err|=3.63e-07 (scale 0.618)



spin A=Sm_0 v=3      max|err|=1.06e-03 (scale 0.168)
spin A=Sm_0 v=python max|err|=2.50e-12 (scale 0.168)



spin A=Sp_0 v=3      max|err|=6.86e-02 (scale 0.275)
spin A=Sp_0 v=python max|err|=4.30e-14 (scale 0.275)



fermion mu=0 A=C_0 v=3 gse_cutoff=1e-08 krylov=3 max|err|=2.87e-02 (scale 0.171)



fermion mu=0 A=C_0 v=3 gse_cutoff=1e-14 krylov=3 max|err|=2.87e-02 (scale 0.171)



fermion mu=0 A=C_0 v=3 gse_cutoff=0 krylov=3 max|err|=4.69e-02 (scale 0.171)



fermion mu=0 A=C_0 v=3 gse_cutoff=1e-08 krylov=5 max|err|=2.87e-02 (scale 0.171)
```

```bash
cd <scratch>/realtime && <scratch>/run3.sh 12_v3_gse_edge_vs_bulk.py 2>&1 | tee 12_v3_gse_edge_vs_bulk.out
```

`realtime/12_v3_gse_edge_vs_bulk.py`:

```python
"""Lens realtime, 12: is the itensor_version=3 TDVP_GSE failure of 10/11
the documented "bond-dimension-1 start" edge case reached through a local
operator on the EDGE site (which leaves A|GS> with Schmidt rank 1 across
the first bond), or does a bulk operator trigger it too?

Same 6-site spinless chain as 11 (mu=0.8, <N>=4) and a 6-site Heisenberg
chain; evolution_ABA(A, B) and evolution_DC under TDVP_GSE, v3 and
"python", against exact (ED matrices, dense eigh, e^{-iHt} by hand).
The Schmidt rank of A|GS> across every bond is printed from the ED vector.
"""
import numpy as np, functools
print = functools.partial(print, flush=True)
from dmrgpy import fermionchain, spinchain, timedependent

n, nt, dt = 6, 40, 0.05
ts = dt*np.arange(nt)

def ferm(version):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    fc.maxm, fc.nsweeps = 40, 20
    fc.tevol_method = "TDVP_GSE"
    h = 0
    for j in range(n - 1):
        h = h - fc.Cdag[j+1]*fc.C[j] - fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0] + 0.5*fc.N[2]*fc.N[3] - 0.8*sum(fc.N[j] for j in range(n))
    fc.set_hamiltonian(h)
    return fc

def spin(version):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    sc.maxm, sc.nsweeps = 40, 20
    sc.tevol_method = "TDVP_GSE"
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc

def ranks(vec):
    out = []
    for b in range(1, n):
        s = np.linalg.svd(vec.reshape(2**b, -1), compute_uv=False)
        out.append(int(np.sum(s > 1e-10*s[0])))
    return out

def run(ch, A, B, kind):
    ed = ch.get_ED_obj()
    H = ed.get_operator(ch.hamiltonian).toarray()
    w, v = np.linalg.eigh(H); g = v[:, 0]
    Am = ed.get_operator(A).toarray(); Bm = ed.get_operator(B).toarray()
    if kind == "ABA":
        psi0 = Am@g; c = np.conj(v).T@psi0
        ref = np.array([np.vdot(v@(np.exp(-1j*w*t)*c), Bm@(v@(np.exp(-1j*w*t)*c)))
                        for t in ts])
        _t, y = timedependent.evolution_ABA(ch, A=A, B=B, mode="DMRG", nt=nt, dt=dt)
        start = psi0
    else:   # evolution_DC name=(A,B): the evolved state is A^dagger|GS>
        M = np.conj(np.conj(v).T@(Am.conj().T@g))*(np.conj(v).T@(Bm@g))
        ref = np.array([np.sum(M*np.exp(1j*(w - w[0])*t)) for t in ts])
        _t, y = timedependent.evolution_DC(ch, mode="DMRG", name=(A, B), nt=nt, dt=dt)
        start = Am.conj().T@g
    return np.max(np.abs(np.asarray(y) - ref)), np.abs(ref).max(), ranks(start)

cases = [
    ("fermion ABA A=C_0  B=N_1", ferm, lambda c: (c.C[0], c.N[1]), "ABA"),
    ("fermion ABA A=C_2  B=N_1", ferm, lambda c: (c.C[2], c.N[1]), "ABA"),
    ("fermion ABA A=C_5  B=N_1", ferm, lambda c: (c.C[5], c.N[1]), "ABA"),
    ("fermion DC (Cdag_2,C_2) ", ferm, lambda c: (c.Cdag[2], c.C[2]), "DC"),
    ("fermion DC (Cdag_2,C_3) ", ferm, lambda c: (c.Cdag[2], c.C[3]), "DC"),
    ("spin ABA A=Sp_0 B=Sz_1  ", spin, lambda c: (c.Sx[0] + 1j*c.Sy[0], c.Sz[1]), "ABA"),
    ("spin ABA A=Sp_2 B=Sz_1  ", spin, lambda c: (c.Sx[2] + 1j*c.Sy[2], c.Sz[1]), "ABA"),
    ("spin DC (Sm_2,Sp_2)     ", spin, lambda c: (c.Sx[2] - 1j*c.Sy[2], c.Sx[2] + 1j*c.Sy[2]), "DC"),
]
for label, mk, pair, kind in cases:
    for version in (3, "python"):
        ch = mk(version)
        A, B = pair(ch)
        e, s, r = run(ch, A, B, kind)
        print("%s v=%-6s max|err|=%.2e (scale %.3f)  Schmidt ranks of the "
              "evolved start: %s" % (label, version, e, s, r))
```

`realtime/12_v3_gse_edge_vs_bulk.out`:


```



fermion ABA A=C_0  B=N_1 v=3      max|err|=3.00e-01 (scale 0.314)  Schmidt ranks of the evolved start: [2, 4, 4, 2, 1]
fermion ABA A=C_0  B=N_1 v=python max|err|=4.06e-07 (scale 0.314)  Schmidt ranks of the evolved start: [2, 4, 4, 2, 1]



fermion ABA A=C_2  B=N_1 v=3      max|err|=2.25e-11 (scale 0.300)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
fermion ABA A=C_2  B=N_1 v=python max|err|=1.02e-11 (scale 0.300)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]



fermion ABA A=C_5  B=N_1 v=3      max|err|=1.50e-07 (scale 0.482)  Schmidt ranks of the evolved start: [1, 2, 4, 4, 2]
fermion ABA A=C_5  B=N_1 v=python max|err|=1.50e-07 (scale 0.482)  Schmidt ranks of the evolved start: [1, 2, 4, 4, 2]



fermion DC (Cdag_2,C_2)  v=3      max|err|=4.64e-11 (scale 0.637)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
fermion DC (Cdag_2,C_2)  v=python max|err|=2.90e-12 (scale 0.637)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]



fermion DC (Cdag_2,C_3)  v=3      max|err|=5.15e-11 (scale 0.339)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
fermion DC (Cdag_2,C_3)  v=python max|err|=9.47e-12 (scale 0.339)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]



spin ABA A=Sp_0 B=Sz_1   v=3      max|err|=6.86e-02 (scale 0.275)  Schmidt ranks of the evolved start: [2, 4, 4, 2, 1]
spin ABA A=Sp_0 B=Sz_1   v=python max|err|=3.77e-14 (scale 0.275)  Schmidt ranks of the evolved start: [2, 4, 4, 2, 1]



spin ABA A=Sp_2 B=Sz_1   v=3      max|err|=1.42e-10 (scale 0.150)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
spin ABA A=Sp_2 B=Sz_1   v=python max|err|=1.69e-13 (scale 0.150)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]



spin DC (Sm_2,Sp_2)      v=3      max|err|=4.70e-11 (scale 0.571)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
spin DC (Sm_2,Sp_2)      v=python max|err|=2.43e-14 (scale 0.571)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
```

```bash
cd <scratch>/realtime && <scratch>/run3.sh 13_v3_gse_public_td_spectrum.py 2>&1 | tee 13_v3_gse_public_td_spectrum.out
```

`realtime/13_v3_gse_public_td_spectrum.py`:

```python
"""Lens realtime, 13: what the itensor_version=3 TDVP_GSE failure of 12
does to the public spectrum, get_dynamical_correlator(submode="TD"), at an
edge-site pair, the default i=j=0 of every string name.

Reference: the house Lehmann density sum_n M_n delta/pi/((w-D_n)^2+delta^2),
M_n = <GS|A|n><n|B|GS>, from ED matrices by dense eigh. Rows: v3 and
"python" under TDVP_GSE and TDVP, and mode="ED", all at predict=False so
the frequency stage is the plain damped sum.
"""
import numpy as np, functools
print = functools.partial(print, flush=True)
from dmrgpy import fermionchain, spinchain

n, dt, delta = 6, 0.1, 0.2
es = np.linspace(-1.0, 6.0, 701)

def ferm(version, method):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    fc.maxm, fc.nsweeps = 40, 20
    fc.tevol_method = method
    h = 0
    for j in range(n - 1):
        h = h - fc.Cdag[j+1]*fc.C[j] - fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0] + 0.5*fc.N[2]*fc.N[3] - 0.8*sum(fc.N[j] for j in range(n))
    fc.set_hamiltonian(h)
    return fc

def spin(version, method):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    sc.maxm, sc.nsweeps = 40, 20
    sc.tevol_method = method
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc

def lehmann(ch, A, B):
    ed = ch.get_ED_obj()
    H = ed.get_operator(ch.hamiltonian).toarray()
    w, v = np.linalg.eigh(H); g = v[:, 0]
    Am = ed.get_operator(A).toarray(); Bm = ed.get_operator(B).toarray()
    M = np.conj(np.conj(v).T@(Am.conj().T@g))*(np.conj(v).T@(Bm@g))
    D = w - w[0]
    return (M[None, :]*delta/np.pi/((es[:, None] - D[None, :])**2 + delta**2)).sum(1)

cases = [("fermion (Cdag_0, C_0)", ferm, lambda c: (c.Cdag[0], c.C[0])),
         ("spin (S+_0, S-_0)    ", spin, lambda c: (c.Sx[0] + 1j*c.Sy[0], c.Sx[0] - 1j*c.Sy[0]))]
for label, mk, pair in cases:
    ref = None
    for mode, version, method in [("ED", "python", "TDVP"), ("DMRG", 3, "TDVP_GSE"),
                                  ("DMRG", 3, "TDVP"), ("DMRG", "python", "TDVP_GSE")]:
        ch = mk(version, method)
        A, B = pair(ch)
        if ref is None:
            ref = lehmann(ch, A, B)
        e, y = ch.get_dynamical_correlator(mode=mode, submode="TD", name=(A, B),
                                           es=es, delta=delta, dt=dt, predict=False)
        y = np.asarray(y)
        i = np.argmax(np.abs(ref))
        print("%s %-4s v=%-6s %-8s  max|y-exact|/peak = %.3e   at the exact peak"
              " w=%.3f: %.4f vs %.4f   integral %.4f vs %.4f"
              % (label, mode, version, method, np.max(np.abs(y - ref))/np.max(np.abs(ref)),
                 es[i], y[i].real, ref[i].real, np.trapezoid(y.real, es), np.trapezoid(ref.real, es)))
```

`realtime/13_v3_gse_public_td_spectrum.out`:


```
fermion (Cdag_0, C_0) ED   v=python TDVP      max|y-exact|/peak = 2.760e-03   at the exact peak w=0.140: 0.4744 vs 0.4756   integral 0.6742 vs 0.6739



fermion (Cdag_0, C_0) DMRG v=3      TDVP_GSE  max|y-exact|/peak = 1.735e+00   at the exact peak w=0.140: 0.0808 vs 0.4756   integral 0.6777 vs 0.6739
fermion (Cdag_0, C_0) DMRG v=3      TDVP      max|y-exact|/peak = 2.760e-03   at the exact peak w=0.140: 0.4744 vs 0.4756   integral 0.6742 vs 0.6739
fermion (Cdag_0, C_0) DMRG v=python TDVP_GSE  max|y-exact|/peak = 2.757e-03   at the exact peak w=0.140: 0.4744 vs 0.4756   integral 0.6742 vs 0.6739
spin (S+_0, S-_0)     ED   v=python TDVP      max|y-exact|/peak = 2.963e-03   at the exact peak w=0.430: 0.2852 vs 0.2860   integral 0.3502 vs 0.3500



spin (S+_0, S-_0)     DMRG v=3      TDVP_GSE  max|y-exact|/peak = 6.066e-03   at the exact peak w=0.430: 0.2856 vs 0.2860   integral 0.3503 vs 0.3500
spin (S+_0, S-_0)     DMRG v=3      TDVP      max|y-exact|/peak = 2.963e-03   at the exact peak w=0.430: 0.2852 vs 0.2860   integral 0.3502 vs 0.3500
spin (S+_0, S-_0)     DMRG v=python TDVP_GSE  max|y-exact|/peak = 2.963e-03   at the exact peak w=0.430: 0.2852 vs 0.2860   integral 0.3502 vs 0.3500
```


**Reviewer (CONFIRMED, NARROWED)**: `12` and `13` reproduce to the printed digit
(fermion `A=C_0`: v3 3.00e-01 against `"python"` 4.06e-07 on a 0.314 scale; `A=C_2`
2.24e-11; `A=C_5` 1.50e-07, equal to `"python"`; spin `A=Sp_0` 6.86e-02 against
2.23e-14; the public (Cdag_0,C_0) TD spectrum 1.735 of its peak, a height of 0.0808
against 0.4756 at w=0.140, an integral of 0.6777 against 0.6739, where ED, v3 TDVP
and `"python"` TDVP_GSE are all at 2.76e-03), but spin (S+_0,S-_0) is only 6.066e-03
against 2.963e-03 for the others, although its evolved state also has a rank-1 first
bond. A site-0 operator zoo, with the site-0 density as the observable and the
final site-0 entropy against the exact one at T=2, sorts it: frozen for `C_0`
(peak-to-peak <N_0> 0.000 against 4.741e-01, entropy 0 against 0.6330), `1-N_0`,
`Sp_0` and `1/2+Sz_0`; degraded but moving for `N_0` (1.82e-02), `Sm_0`
(4.46e-04) and `1/2-Sz_0` (1.18e-02), where `"python"` is at 5.8e-08, 5.9e-13 and
1.5e-09; and rescued by any admixture (`C_0+0.3C_2`, 1.26e-10). A verbose single
step shows `addBasis` raising the maximum link dimension to 7 with no warning at the
edge bond, after which the TDVP sweep keeps 1,2,4,4,2, the state's own ranks, with a
site-0 entropy of 2.2e-16; for `C_5` it keeps 2,4,7,4,2. So the expansion is
computed and then lost between `addBasis` and the first TDVP update; which operation
drops it could not be isolated, since no link-dimension getter is exposed on the
MPS handle. On 10 sites (`maxm=60`, B=N_4), bulk `A=C_4` gives 1.15e-09 at
`tdvp_gse_sweeps=3` against 1.04e-05 at 0, so the v3 expansion works and matters in
the bulk, while edge `A=C_0` gives 3.01e-03 at both, with a site-0 entropy of 0
against 0.6922. `git log -S global_subspace_expand` confirms `19bbee0`. Struck, each
explicitly:

- "Does no subspace expansion at all": `addBasis` adds directions, and bulk
  expansion on v3 works; what holds is that the expansion is lost before the first
  update on an edge-pinned start.
- "Rank 1 on the bond next to site 0" as the discriminant: pinning to the first
  local basis state freezes, pinning to the second degrades, both with a rank-1 edge
  bond.
- "Every ladder or projector operator on site 0" producing the large failure: all
  are affected, but only the first-basis-state ones freeze.
- The quoted sizes as general: they belong to the hunter's 6-site chains; on 10
  sites with B=N_4 the frozen branch is 3.01e-3 on a 0.314 scale.

```bash
cd <scratch>/reviews/realtime_1 && <scratch>/run3.sh 01_repro_edge_vs_bulk.py 2>&1 | tee 01_repro_edge_vs_bulk.out
```

`reviews/realtime_1/01_repro_edge_vs_bulk.py` is `realtime/12_v3_gse_edge_vs_bulk.py` (finding 14) unchanged, rerun; its output, `reviews/realtime_1/01_repro_edge_vs_bulk.out`:


```
fermion ABA A=C_0  B=N_1 v=3      max|err|=3.00e-01 (scale 0.314)  Schmidt ranks of the evolved start: [2, 4, 4, 2, 1]
fermion ABA A=C_0  B=N_1 v=python max|err|=4.06e-07 (scale 0.314)  Schmidt ranks of the evolved start: [2, 4, 4, 2, 1]
fermion ABA A=C_2  B=N_1 v=3      max|err|=2.24e-11 (scale 0.300)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
fermion ABA A=C_2  B=N_1 v=python max|err|=1.06e-11 (scale 0.300)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
fermion ABA A=C_5  B=N_1 v=3      max|err|=1.50e-07 (scale 0.482)  Schmidt ranks of the evolved start: [1, 2, 4, 4, 2]
fermion ABA A=C_5  B=N_1 v=python max|err|=1.50e-07 (scale 0.482)  Schmidt ranks of the evolved start: [1, 2, 4, 4, 2]
fermion DC (Cdag_2,C_2)  v=3      max|err|=4.10e-11 (scale 0.637)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
fermion DC (Cdag_2,C_2)  v=python max|err|=2.89e-12 (scale 0.637)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
fermion DC (Cdag_2,C_3)  v=3      max|err|=5.21e-11 (scale 0.339)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
fermion DC (Cdag_2,C_3)  v=python max|err|=9.46e-12 (scale 0.339)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
spin ABA A=Sp_0 B=Sz_1   v=3      max|err|=6.86e-02 (scale 0.275)  Schmidt ranks of the evolved start: [2, 4, 4, 2, 1]
spin ABA A=Sp_0 B=Sz_1   v=python max|err|=2.23e-14 (scale 0.275)  Schmidt ranks of the evolved start: [2, 4, 4, 2, 1]
spin ABA A=Sp_2 B=Sz_1   v=3      max|err|=1.42e-10 (scale 0.150)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
spin ABA A=Sp_2 B=Sz_1   v=python max|err|=1.68e-13 (scale 0.150)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
spin DC (Sm_2,Sp_2)      v=3      max|err|=4.48e-11 (scale 0.571)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
spin DC (Sm_2,Sp_2)      v=python max|err|=3.92e-13 (scale 0.571)  Schmidt ranks of the evolved start: [2, 4, 4, 4, 2]
```

```bash
cd <scratch>/reviews/realtime_1 && <scratch>/run3.sh 02_repro_public_td_spectrum.py 2>&1 | tee 02_repro_public_td_spectrum.out
```

`reviews/realtime_1/02_repro_public_td_spectrum.py` is `realtime/13_v3_gse_public_td_spectrum.py` (finding 14) unchanged, rerun; its output, `reviews/realtime_1/02_repro_public_td_spectrum.out`:


```
fermion (Cdag_0, C_0) ED   v=python TDVP      max|y-exact|/peak = 2.760e-03   at the exact peak w=0.140: 0.4744 vs 0.4756   integral 0.6742 vs 0.6739
fermion (Cdag_0, C_0) DMRG v=3      TDVP_GSE  max|y-exact|/peak = 1.735e+00   at the exact peak w=0.140: 0.0808 vs 0.4756   integral 0.6777 vs 0.6739
fermion (Cdag_0, C_0) DMRG v=3      TDVP      max|y-exact|/peak = 2.760e-03   at the exact peak w=0.140: 0.4744 vs 0.4756   integral 0.6742 vs 0.6739
fermion (Cdag_0, C_0) DMRG v=python TDVP_GSE  max|y-exact|/peak = 2.757e-03   at the exact peak w=0.140: 0.4744 vs 0.4756   integral 0.6742 vs 0.6739
spin (S+_0, S-_0)     ED   v=python TDVP      max|y-exact|/peak = 2.963e-03   at the exact peak w=0.430: 0.2852 vs 0.2860   integral 0.3502 vs 0.3500
spin (S+_0, S-_0)     DMRG v=3      TDVP_GSE  max|y-exact|/peak = 6.066e-03   at the exact peak w=0.430: 0.2856 vs 0.2860   integral 0.3503 vs 0.3500
spin (S+_0, S-_0)     DMRG v=3      TDVP      max|y-exact|/peak = 2.963e-03   at the exact peak w=0.430: 0.2852 vs 0.2860   integral 0.3502 vs 0.3500
spin (S+_0, S-_0)     DMRG v=python TDVP_GSE  max|y-exact|/peak = 2.963e-03   at the exact peak w=0.430: 0.2852 vs 0.2860   integral 0.3502 vs 0.3500
```

```bash
cd <scratch>/reviews/realtime_1 && <scratch>/run3.sh 03_site0_operator_zoo.py 2>&1 | tee 03_site0_operator_zoo.out
```

`reviews/realtime_1/03_site0_operator_zoo.py`:

```python
"""Review of realtime_1, 03: which site-0 operators break itensor_version=3
TDVP_GSE, and is site 0 frozen when it breaks?

Same 6-site chains as the hunter's 12. For each A acting on site 0 (and a
few controls), evolution_ABA(A, B) on v3 and "python" TDVP_GSE, against
exact (ED matrices, dense eigh). B is the site-0 density (N_0 / Sz_0), so a
frozen site 0 shows up as a flat trajectory. With return_wf=True the final
state's site-0 entropy (= the entanglement across the first bond) is
compared with the exact one at t = nt*dt.
"""
import numpy as np, functools
print = functools.partial(print, flush=True)
from dmrgpy import fermionchain, spinchain, timedependent

n, nt, dt = 6, 40, 0.05
ts = dt*np.arange(nt)

def ferm(version):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    fc.maxm, fc.nsweeps = 40, 20
    fc.tevol_method = "TDVP_GSE"
    h = 0
    for j in range(n - 1):
        h = h - fc.Cdag[j+1]*fc.C[j] - fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0] + 0.5*fc.N[2]*fc.N[3] - 0.8*sum(fc.N[j] for j in range(n))
    fc.set_hamiltonian(h)
    return fc

def spin(version):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    sc.maxm, sc.nsweeps = 40, 20
    sc.tevol_method = "TDVP_GSE"
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    return sc

def ent(p):
    p = p[p > 1e-14]
    return float(-np.sum(p*np.log(p)))

def site_ent_exact(ed, ch, vec):
    # site-0 reduced density matrix from ED operators, independent of the
    # basis ordering: rho_ab = <vec| |b><a|_0 |vec>, built from N_0/Sz_0 and
    # the site-0 ladder operator
    vec = vec/np.linalg.norm(vec)
    return vec

def run(ch, A, B, P, L):
    """P: site-0 projector-like diagonal op (N_0 or Sz_0+1/2), L: site-0
    lowering-type op, used to build the exact site-0 density matrix."""
    ed = ch.get_ED_obj()
    H = ed.get_operator(ch.hamiltonian).toarray()
    w, v = np.linalg.eigh(H); g = v[:, 0]
    Am = ed.get_operator(A).toarray(); Bm = ed.get_operator(B).toarray()
    Pm = ed.get_operator(P).toarray(); Lm = ed.get_operator(L).toarray()
    psi0 = Am@g
    c = np.conj(v).T@psi0
    ref = np.array([np.vdot(v@(np.exp(-1j*w*t)*c), Bm@(v@(np.exp(-1j*w*t)*c))) for t in ts])
    psiT = v@(np.exp(-1j*w*nt*dt)*c); psiT = psiT/np.linalg.norm(psiT)
    p1 = np.vdot(psiT, Pm@psiT).real; off = np.vdot(psiT, Lm@psiT)
    rho = np.array([[p1, off], [np.conj(off), 1 - p1]])
    Sex = ent(np.linalg.eigvalsh(rho))
    _t, y, wf = timedependent.evolution_ABA(ch, A=A, B=B, mode="DMRG", nt=nt, dt=dt,
                                            return_wf=True)
    y = np.asarray(y)
    Sdm = ch.get_site_entropy(wf, 0)
    return np.max(np.abs(y - ref)), np.abs(ref).max(), np.ptp(y.real), np.ptp(ref.real), Sdm, Sex

fcases = [("C_0", lambda c: c.C[0]), ("Cdag_0", lambda c: c.Cdag[0]),
          ("N_0", lambda c: c.N[0]), ("1-N_0", lambda c: 1 - c.N[0]),
          ("C_1", lambda c: c.C[1]), ("C_0+0.3C_2", lambda c: c.C[0] + 0.3*c.C[2])]
scases = [("Sp_0", lambda c: c.Sx[0] + 1j*c.Sy[0]), ("Sm_0", lambda c: c.Sx[0] - 1j*c.Sy[0]),
          ("1/2+Sz_0", lambda c: 0.5 + c.Sz[0]), ("1/2-Sz_0", lambda c: 0.5 - c.Sz[0]),
          ("Sp_1", lambda c: c.Sx[1] + 1j*c.Sy[1])]
for label, mkA in fcases:
    for version in (3, "python"):
        ch = ferm(version)
        e, s, pty, ptr, Sdm, Sex = run(ch, mkA(ch), ch.N[0], ch.N[0], ch.C[0])
        print("fermion A=%-11s B=N_0  v=%-6s max|err|=%.2e (scale %.3f)  ptp<N_0>(t) dmrg %.3e exact %.3e"
              "  S_site0(T) dmrg %.4f exact %.4f" % (label, version, e, s, pty, ptr, Sdm, Sex))
for label, mkA in scases:
    for version in (3, "python"):
        ch = spin(version)
        e, s, pty, ptr, Sdm, Sex = run(ch, mkA(ch), ch.Sz[0], 0.5 + ch.Sz[0], ch.Sx[0] - 1j*ch.Sy[0])
        print("spin    A=%-11s B=Sz_0 v=%-6s max|err|=%.2e (scale %.3f)  ptp<Sz_0>(t) dmrg %.3e exact %.3e"
              "  S_site0(T) dmrg %.4f exact %.4f" % (label, version, e, s, pty, ptr, Sdm, Sex))
```

`reviews/realtime_1/03_site0_operator_zoo.out`:


```
fermion A=C_0         B=N_0  v=3      max|err|=4.74e-01 (scale 0.474)  ptp<N_0>(t) dmrg 0.000e+00 exact 4.741e-01  S_site0(T) dmrg -0.0000 exact 0.6330
fermion A=C_0         B=N_0  v=python max|err|=1.08e-07 (scale 0.474)  ptp<N_0>(t) dmrg 4.741e-01 exact 4.741e-01  S_site0(T) dmrg 0.6330 exact 0.6330
fermion A=Cdag_0      B=N_0  v=3      max|err|=2.05e-04 (scale 0.290)  ptp<N_0>(t) dmrg 1.770e-02 exact 1.775e-02  S_site0(T) dmrg 0.2365 exact 0.2371
fermion A=Cdag_0      B=N_0  v=python max|err|=3.11e-11 (scale 0.290)  ptp<N_0>(t) dmrg 1.775e-02 exact 1.775e-02  S_site0(T) dmrg 0.2371 exact 0.2371
fermion A=N_0         B=N_0  v=3      max|err|=1.82e-02 (scale 0.710)  ptp<N_0>(t) dmrg 2.524e-01 exact 2.622e-01  S_site0(T) dmrg 0.6013 exact 0.6042
fermion A=N_0         B=N_0  v=python max|err|=5.84e-08 (scale 0.710)  ptp<N_0>(t) dmrg 2.622e-01 exact 2.622e-01  S_site0(T) dmrg 0.6042 exact 0.6042
fermion A=1-N_0       B=N_0  v=3      max|err|=2.66e-01 (scale 0.266)  ptp<N_0>(t) dmrg 0.000e+00 exact 2.663e-01  S_site0(T) dmrg 0.0000 exact 0.3249
fermion A=1-N_0       B=N_0  v=python max|err|=9.28e-08 (scale 0.266)  ptp<N_0>(t) dmrg 2.663e-01 exact 2.663e-01  S_site0(T) dmrg 0.3249 exact 0.3249
fermion A=C_1         B=N_0  v=3      max|err|=1.99e-11 (scale 0.314)  ptp<N_0>(t) dmrg 6.265e-02 exact 6.265e-02  S_site0(T) dmrg 0.6930 exact 0.6930
fermion A=C_1         B=N_0  v=python max|err|=2.50e-12 (scale 0.314)  ptp<N_0>(t) dmrg 6.265e-02 exact 6.265e-02  S_site0(T) dmrg 0.6930 exact 0.6930
fermion A=C_0+0.3C_2  B=N_0  v=3      max|err|=1.26e-10 (scale 0.449)  ptp<N_0>(t) dmrg 4.122e-01 exact 4.122e-01  S_site0(T) dmrg 0.6106 exact 0.6106
fermion A=C_0+0.3C_2  B=N_0  v=python max|err|=6.90e-11 (scale 0.449)  ptp<N_0>(t) dmrg 4.122e-01 exact 4.122e-01  S_site0(T) dmrg 0.6106 exact 0.6106
spin    A=Sp_0        B=Sz_0 v=3      max|err|=6.97e-02 (scale 0.316)  ptp<Sz_0>(t) dmrg 9.437e-16 exact 6.966e-02  S_site0(T) dmrg -0.0000 exact 0.3621
spin    A=Sp_0        B=Sz_0 v=python max|err|=4.96e-14 (scale 0.316)  ptp<Sz_0>(t) dmrg 6.966e-02 exact 6.966e-02  S_site0(T) dmrg 0.3621 exact 0.3621
spin    A=Sm_0        B=Sz_0 v=3      max|err|=4.46e-04 (scale 0.184)  ptp<Sz_0>(t) dmrg 2.620e-02 exact 2.575e-02  S_site0(T) dmrg 0.2685 exact 0.2652
spin    A=Sm_0        B=Sz_0 v=python max|err|=5.90e-13 (scale 0.184)  ptp<Sz_0>(t) dmrg 2.575e-02 exact 2.575e-02  S_site0(T) dmrg 0.2652 exact 0.2652
spin    A=1/2+Sz_0    B=Sz_0 v=3      max|err|=2.08e-01 (scale 0.184)  ptp<Sz_0>(t) dmrg 5.551e-16 exact 2.079e-01  S_site0(T) dmrg 0.0000 exact 0.6801
spin    A=1/2+Sz_0    B=Sz_0 v=python max|err|=2.61e-09 (scale 0.184)  ptp<Sz_0>(t) dmrg 2.079e-01 exact 2.079e-01  S_site0(T) dmrg 0.6801 exact 0.6801
spin    A=1/2-Sz_0    B=Sz_0 v=3      max|err|=1.18e-02 (scale 0.316)  ptp<Sz_0>(t) dmrg 3.163e-01 exact 3.252e-01  S_site0(T) dmrg 0.6928 exact 0.6918
spin    A=1/2-Sz_0    B=Sz_0 v=python max|err|=1.48e-09 (scale 0.316)  ptp<Sz_0>(t) dmrg 3.252e-01 exact 3.252e-01  S_site0(T) dmrg 0.6918 exact 0.6918
spin    A=Sp_1        B=Sz_0 v=3      max|err|=1.49e-11 (scale 0.188)  ptp<Sz_0>(t) dmrg 3.330e-02 exact 3.330e-02  S_site0(T) dmrg 0.1633 exact 0.1633
spin    A=Sp_1        B=Sz_0 v=python max|err|=2.32e-13 (scale 0.188)  ptp<Sz_0>(t) dmrg 3.330e-02 exact 3.330e-02  S_site0(T) dmrg 0.1633 exact 0.1633
```

```bash
cd <scratch>/reviews/realtime_1 && <scratch>/run3.sh 04_v3_gse_verbose.py > 04_v3_gse_verbose.out 2>&1
```

`reviews/realtime_1/04_v3_gse_verbose.py`:

```python
"""Review of realtime_1, 04: what addBasis reports per bond on
itensor_version=3, with verbose on, for one GSE step from C_0|GS> (site 0
pinned to the empty state, which freezes in 03), from Cdag_0|GS> (pinned
occupied, which does not freeze) and from C_5|GS> (right edge, works).
ITensor bond b is between ITensor sites b and b+1, i.e. dmrgpy sites b-1, b.
"""
import sys, numpy as np, functools
print = functools.partial(print, flush=True)
from dmrgpy import fermionchain, timedependent

n = 6
def ferm():
    fc = fermionchain.Fermionic_Chain(n, itensor_version=3)
    fc.maxm, fc.nsweeps = 40, 20
    fc.tevol_method = "TDVP_GSE"
    h = 0
    for j in range(n - 1):
        h = h - fc.Cdag[j+1]*fc.C[j] - fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0] + 0.5*fc.N[2]*fc.N[3] - 0.8*sum(fc.N[j] for j in range(n))
    fc.set_hamiltonian(h)
    return fc

for label, mkA in [("C_0", lambda c: c.C[0]), ("Cdag_0", lambda c: c.Cdag[0]),
                   ("C_5", lambda c: c.C[5])]:
    fc = ferm()
    wf = fc.get_gs()
    fc.tdvp_gse_sweeps = 1
    fc.verbose = True
    print("===== A=%s: one GSE step, verbose" % label)
    sys.stdout.flush()
    _t, y, wfT = timedependent.evolution_ABA(fc, A=mkA(fc), B=fc.N[0], mode="DMRG",
                                             nt=1, dt=0.05, wf=wf, return_wf=True)
    sys.stdout.flush()
    fc.verbose = False
    print("===== A=%s done; site-0 entropy after the step %.3e" % (label, fc.get_site_entropy(wfT, 0)))
```

`reviews/realtime_1/04_v3_gse_verbose.out`:


```
===== A=C_0: one GSE step, verbose
norm(psi1)=1.00000000000000000000
maxLinkDim(psi1) = 4
norm(psi2)=1.00000000000000000000
maxLinkDim(psi2) = 4
warning: at bond 6, already reach maximum bond dimension.
warning: at bond 5, already reach maximum bond dimension.

maxLinkDim after global subspace expansion = 7
Global subspace expansion: cputime = 0.00163s, walltime = 0.00262s
Sweep=1, HS=1, Bond=1/5
In applyExp, number of matrix-vector multiplies: 1
In applyExp, number of matrix-vector multiplies: 1
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
Sweep=1, HS=1, Bond=2/5
In applyExp, number of matrix-vector multiplies: 2
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=1, Bond=3/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4
Sweep=1, HS=1, Bond=4/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4
Sweep=1, HS=1, Bond=5/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=1, Bond=6/5
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
Sweep=1, HS=2, Bond=6/5
In applyExp, number of matrix-vector multiplies: 2
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
Sweep=1, HS=2, Bond=5/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=2, Bond=4/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4
Sweep=1, HS=2, Bond=3/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4

    vN Entropy at center bond b=3 = 0.686670577611
    Eigs at center bond b=3: 0.5568 0.4432 
Sweep=1, HS=2, Bond=2/5
In applyExp, number of matrix-vector multiplies: 2
In applyExp, number of matrix-vector multiplies: 1
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=2, Bond=1/5
In applyExp, number of matrix-vector multiplies: 1
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
    Largest link dim during sweep 1/1 was 4
    Largest truncation error: 0
    Energy after sweep 1/1 is -4.882654859370
    Sweep 1/1 CPU time = 0.00462s (Wall time = 0.00463s)
===== A=C_0 done; site-0 entropy after the step 2.220e-16
===== A=Cdag_0: one GSE step, verbose
norm(psi1)=1.00000000000000000000
maxLinkDim(psi1) = 2
norm(psi2)=1.00000000000000000000
maxLinkDim(psi2) = 2
warning: at bond 6, already reach maximum bond dimension.

maxLinkDim after global subspace expansion = 4
Global subspace expansion: cputime = 0.00241s, walltime = 0.00241s
Sweep=1, HS=1, Bond=1/5
In applyExp, number of matrix-vector multiplies: 1
In applyExp, number of matrix-vector multiplies: 1
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
Sweep=1, HS=1, Bond=2/5
In applyExp, number of matrix-vector multiplies: 2
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=1, Bond=3/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4
Sweep=1, HS=1, Bond=4/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=3
Sweep=1, HS=1, Bond=5/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=1, Bond=6/5
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
Sweep=1, HS=2, Bond=6/5
In applyExp, number of matrix-vector multiplies: 2
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
Sweep=1, HS=2, Bond=5/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=2, Bond=4/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=3
Sweep=1, HS=2, Bond=3/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4

    vN Entropy at center bond b=3 = 0.101159594699
    Eigs at center bond b=3: 0.9792 0.0208 
Sweep=1, HS=2, Bond=2/5
In applyExp, number of matrix-vector multiplies: 2
In applyExp, number of matrix-vector multiplies: 1
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=2, Bond=1/5
In applyExp, number of matrix-vector multiplies: 1
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
    Largest link dim during sweep 1/1 was 4
    Largest truncation error: 0
    Energy after sweep 1/1 is -5.204754934857
    Sweep 1/1 CPU time = 0.00437s (Wall time = 0.00437s)
===== A=Cdag_0 done; site-0 entropy after the step 1.110e-16
===== A=C_5: one GSE step, verbose
norm(psi1)=0.99999999999999988898
maxLinkDim(psi1) = 4
norm(psi2)=1.00000000000000000000
maxLinkDim(psi2) = 4

maxLinkDim after global subspace expansion = 11
Global subspace expansion: cputime = 0.00174s, walltime = 0.00273s
Sweep=1, HS=1, Bond=1/5
In applyExp, number of matrix-vector multiplies: 2
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=1, Bond=2/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4
Sweep=1, HS=1, Bond=3/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=7
Sweep=1, HS=1, Bond=4/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4
Sweep=1, HS=1, Bond=5/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=1, Bond=6/5
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
Sweep=1, HS=2, Bond=6/5
In applyExp, number of matrix-vector multiplies: 2
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=1
Sweep=1, HS=2, Bond=5/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
Sweep=1, HS=2, Bond=4/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4
Sweep=1, HS=2, Bond=3/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=7

    vN Entropy at center bond b=3 = 0.641268010870
    Eigs at center bond b=3: 0.7218 0.2680 0.0072 0.0030 
Sweep=1, HS=2, Bond=2/5
In applyExp, number of iterations: 3
In applyExp, number of matrix-vector multiplies: 4
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=4
Sweep=1, HS=2, Bond=1/5
In applyExp, number of matrix-vector multiplies: 2
    Truncated to Cutoff=1.0E-12, Min_dim=1, Max_dim=40
    Trunc. err=0.0E+00, States kept: dim=2
    Largest link dim during sweep 1/1 was 7
    Largest truncation error: 0
    Energy after sweep 1/1 is -4.757609409656
    Sweep 1/1 CPU time = 0.00526s (Wall time = 0.00526s)
===== A=C_5 done; site-0 entropy after the step 6.249e-01
```

```bash
cd <scratch>/reviews/realtime_1 && <scratch>/run3.sh 05_bulk_start_longer_chain.py 2>&1 | tee 05_bulk_start_longer_chain.out
```

`reviews/realtime_1/05_bulk_start_longer_chain.py`:

```python
"""Review of realtime_1, 05: does itensor_version=3 TDVP_GSE actually expand
on a bulk start when the chain is long enough that the bond dimensions do
not saturate? (On the 6-site chain of 03 a bulk start works, but there the
ranks are close to full, so any directions at all would do.)

10-site spinless chain, same Hamiltonian shape as the hunter's; A=C_4 (bulk)
and A=C_0 (edge), B=N_4, evolution_ABA against exact. Rows: v3 and "python"
TDVP_GSE at tdvp_gse_sweeps 3 and 0, and v3 two-site TDVP.
"""
import numpy as np, functools
print = functools.partial(print, flush=True)
from dmrgpy import fermionchain, timedependent

n, nt, dt = 10, 40, 0.05
ts = dt*np.arange(nt)

def ferm(version, method, sweeps):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    fc.maxm, fc.nsweeps = 60, 20
    fc.tevol_method = method
    fc.tdvp_gse_sweeps = sweeps
    h = 0
    for j in range(n - 1):
        h = h - fc.Cdag[j+1]*fc.C[j] - fc.Cdag[j]*fc.C[j+1]
    h = h + 0.3*fc.N[0] + 0.5*fc.N[4]*fc.N[5] - 0.8*sum(fc.N[j] for j in range(n))
    fc.set_hamiltonian(h)
    return fc

ref_cache = {}
for site in (4, 0):
    for version, method, sweeps in [(3, "TDVP_GSE", 3), (3, "TDVP_GSE", 0), (3, "TDVP", 0),
                                    ("python", "TDVP_GSE", 3), ("python", "TDVP_GSE", 0)]:
        ch = ferm(version, method, sweeps)
        A, B = ch.C[site], ch.N[4]
        if site not in ref_cache:
            ed = ch.get_ED_obj()
            H = ed.get_operator(ch.hamiltonian).toarray()
            w, v = np.linalg.eigh(H); g = v[:, 0]
            c = np.conj(v).T@(ed.get_operator(A).toarray()@g)
            Bm = ed.get_operator(B).toarray()
            ref_cache[site] = np.array([np.vdot(v@(np.exp(-1j*w*t)*c), Bm@(v@(np.exp(-1j*w*t)*c)))
                                        for t in ts])
        ref = ref_cache[site]
        _t, y, wf = timedependent.evolution_ABA(ch, A=A, B=B, mode="DMRG", nt=nt, dt=dt,
                                                return_wf=True)
        y = np.asarray(y)
        print("A=C_%d v=%-6s %-8s gse_sweeps=%d  max|err|=%.2e (scale %.3f)  S_site0(T)=%.4f"
              % (site, version, method, sweeps, np.max(np.abs(y - ref)), np.abs(ref).max(),
                 ch.get_site_entropy(wf, 0)))
```

`reviews/realtime_1/05_bulk_start_longer_chain.out`:


```
A=C_4 v=3      TDVP_GSE gse_sweeps=3  max|err|=1.15e-09 (scale 0.344)  S_site0(T)=0.6714
A=C_4 v=3      TDVP_GSE gse_sweeps=0  max|err|=1.04e-05 (scale 0.344)  S_site0(T)=0.6714
A=C_4 v=3      TDVP     gse_sweeps=0  max|err|=6.51e-10 (scale 0.344)  S_site0(T)=0.6714
A=C_4 v=python TDVP_GSE gse_sweeps=3  max|err|=1.03e-09 (scale 0.344)  S_site0(T)=0.6714
A=C_4 v=python TDVP_GSE gse_sweeps=0  max|err|=1.27e-05 (scale 0.344)  S_site0(T)=0.6714
A=C_0 v=3      TDVP_GSE gse_sweeps=3  max|err|=3.01e-03 (scale 0.314)  S_site0(T)=0.0000
A=C_0 v=3      TDVP_GSE gse_sweeps=0  max|err|=3.01e-03 (scale 0.314)  S_site0(T)=0.0000
A=C_0 v=3      TDVP     gse_sweeps=0  max|err|=4.41e-08 (scale 0.314)  S_site0(T)=0.6922
A=C_0 v=python TDVP_GSE gse_sweeps=3  max|err|=1.58e-08 (scale 0.314)  S_site0(T)=0.6922
A=C_0 v=python TDVP_GSE gse_sweeps=0  max|err|=3.01e-03 (scale 0.314)  S_site0(T)=-0.0000
```


**Suggested fix**: the mechanism is not settled, and the fixer's first step is to
print `linkDims(phi)` after each `position()` sweep on `C_0|GS>` against `C_5|GS>`.
The loss happens between `addBasis` returning and the first one-site update, where
the only operations are `global_subspace_expand`'s two `phi.position(...,
{"Cutoff",0.0,"MaxDim",maxm_})` sweeps and TDVP's own `psi.position(1)`, a no-op
there. That points at the two sweeps: they move the centre by SVD, and ITensor v3's
`truncate()` runs `while(truncerr+P(n) <= cutoff*scale && n >= mindim)`, which
discards exactly-zero weights even at `Cutoff=0`; the comment in `chain_session.h`
claiming `Cutoff=0` "only trims once the count exceeds maxm_" is false for exact
zeros, and a site that is exactly one basis vector plausibly yields exact zeros for
one basis ordering and roundoff for the other. The other candidate is that the
edge bond's added direction ends up decoupled, so TDVP's own `Cutoff=1e-12` SVD
drops it. The likely cure is in dmrgpy's own code rather than in the vendored
header: skip both `position()` sweeps unless `maxLinkDim(phi) > maxm_`, since
`addBasis` already returns the centre at site 1 with right-orthonormal expanded
tensors, or, when truncation is needed, use a negative `Cutoff` so exact zeros
survive, or better, cap per bond inside the expansion the way `pyitensor/gse.py`
does with `bond_maxdim`; `pyitensor/gse.py` is exact on every row here and is the
reference. The hunter's fallback of one two-site step whenever a link has dimension
1 is a workaround, not a fix: it swaps integrators mid-trajectory and keys on the
wrong discriminant. The regression should run v3 `TDVP_GSE` against exact from edge
starts pinned to both basis states (`C_0` and `Cdag_0` on a half-filled chain,
`S+_0` and `S-_0`), through `evolution_ABA` and through `quench_tdvp_gse` via
`submode="TD"`, including the site-0 entropy check, which discriminates sharply.
Needs a rebuild of the v3 extension. The example comment, the `documentation.md`
paragraph and the `ROADMAP.md` line should point here. NUMBERS CHANGE only for
v3 with `tevol_method="TDVP_GSE"` from a site-0-pinned start, towards exact: the
(Cdag_0,C_0) TD peak from 0.081 to about 0.474 on the 6-site chain, and
`evolution_ABA(A=C_0)` from 3.0e-01 to about 4e-7 off exact.

### 15. `evolve_and_measure(mode="DMRG")` and `evolution_ABA(mode="DMRG")` take a `**kwargs` that nothing reads on v2, v3 and `"python"`, so a misspelled `DT=0.2` runs at `dt=1e-2` without a word (the last point at t=0.19 instead of the intended 3.8, <Sx_0> = 0.4910 against -0.3955 there) while `mode="ED"` raises on the same call

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `realtime`

**Where**: `src/dmrgpy/timedependent.py:215-216` (the `evolve_and_measure_dmrg`
signature; `kwargs` reaches only the `julia_live` branch at `:262`), `:310-317`
(`evolution_ABA` forwards `**kwargs` into it); `mpsjulialive/timedependent.py:
104-105` (the same signature, by reading).

This is `documentation.md` section 4.10's "`**kwargs` with no consumer". The DMRG
body never reads `kwargs` outside the `julia_live` branch, and the `julia_live`
function names it only on its signature line, so a typo of `dt`, `nt` or
`tevol_method` is dropped and the run goes ahead at the defaults. `867e2b4` edited
this function's return and docstring and left the signature alone. It survived
because every test and example spells the keywords right, and because the returned
`ts` shows the `dt` actually used, so a caller who plots against the returned `ts`
sees the right curve over the wrong window.

**Expected**: Larmor precession, H = B sum Sz from the +x state, so <Sx_0>(t) =
cos(Bt)/2, and a `TypeError` on an unknown keyword, as `mode="ED"` gives and as
`evolution_DC` gives on DMRG.

Repro, from the hunter:

```bash
cd <scratch>/realtime && <scratch>/run3.sh 06_evolve_kwargs_swallowed.py 2>&1 | tee 06_evolve_kwargs_swallowed.out
```

`realtime/06_evolve_kwargs_swallowed.py`:

```python
"""Lens realtime, lead 6: evolve_and_measure_dmrg(..., **kwargs) has no
consumer of **kwargs on the compiled and "python" backends, so a
misspelled keyword is dropped and the run proceeds at the defaults
(dt=1e-2, nt=1000), where ED (evolution_ABC has no **kwargs) raises.
evolution_ABA(mode="DMRG") forwards its **kwargs into the same function.

Anchor: Larmor precession, H = B sum Sz, start |+x...+x>, so
<Sx_0>(t) = cos(Bt)/2 exactly; the caller asks for dt=0.2 (spelled DT=).
"""
import numpy as np
from dmrgpy import spinchain, timedependent, cppext

B, n = 1.0, 3
def chain(version):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    sc.maxm, sc.nsweeps = 10, 10
    sc.set_hamiltonian(sum(-sc.Sx[i] for i in range(n)))
    return sc

for version in ["python"] + ([3] if cppext.available(3) else []) + \
               ([2] if cppext.available(2) else []):
    for mode in ("DMRG", "ED"):
        sc = chain(version)
        wf = sc.get_gs(mode=mode)          # +x polarized
        sc.set_hamiltonian(B*sum(sc.Sz[i] for i in range(n)))
        for label, kw in [("DT=0.2 (typo of dt)", dict(nt=20, DT=0.2)),
                          ("n_t=20 (typo of nt)", dict(dt=0.2, n_t=20)),
                          ("tevol_method='TEBD' as a kwarg",
                           dict(nt=20, dt=0.2, tevol_method="TEBD"))]:
            try:
                ts, y = timedependent.evolve_and_measure(
                    sc, operator=sc.Sx[0], wf=wf, mode=mode, **kw)
            except TypeError as exc:
                print("[%s %s] %-32s raises TypeError: %s"
                      % (version, mode, label, str(exc)[:70]))
                continue
            y = np.asarray(y)
            intended = 0.2*np.arange(20)
            m = min(len(y), 20)
            print("[%s %s] %-32s returned len=%d, ts[1]=%.3f, ts[-1]=%.2f;"
                  " <Sx0>(t=%.1f) returned %.4f, cos(Bt)/2 at the intended"
                  " t = %.4f" % (version, mode, label, len(y), ts[1], ts[-1],
                                  intended[m-1], y[m-1].real,
                                  np.cos(B*intended[m-1])/2))
```

`realtime/06_evolve_kwargs_swallowed.out`:


```
[python DMRG] DT=0.2 (typo of dt)              returned len=20, ts[1]=0.010, ts[-1]=0.19; <Sx0>(t=3.8) returned 0.4910, cos(Bt)/2 at the intended t = -0.3955
[python DMRG] n_t=20 (typo of nt)              returned len=1000, ts[1]=0.200, ts[-1]=199.80; <Sx0>(t=3.8) returned -0.3955, cos(Bt)/2 at the intended t = -0.3955
[python DMRG] tevol_method='TEBD' as a kwarg   returned len=20, ts[1]=0.200, ts[-1]=3.80; <Sx0>(t=3.8) returned -0.3955, cos(Bt)/2 at the intended t = -0.3955
[python ED] DT=0.2 (typo of dt)              raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 
[python ED] n_t=20 (typo of nt)              raises TypeError: evolution_ABC() got an unexpected keyword argument 'n_t'. Did you mean
[python ED] tevol_method='TEBD' as a kwarg   raises TypeError: evolution_ABC() got an unexpected keyword argument 'tevol_method'
[3 DMRG] DT=0.2 (typo of dt)              returned len=20, ts[1]=0.010, ts[-1]=0.19; <Sx0>(t=3.8) returned 0.4910, cos(Bt)/2 at the intended t = -0.3955
[3 DMRG] n_t=20 (typo of nt)              returned len=1000, ts[1]=0.200, ts[-1]=199.80; <Sx0>(t=3.8) returned -0.3955, cos(Bt)/2 at the intended t = -0.3955
[3 DMRG] tevol_method='TEBD' as a kwarg   returned len=20, ts[1]=0.200, ts[-1]=3.80; <Sx0>(t=3.8) returned -0.3955, cos(Bt)/2 at the intended t = -0.3955
[3 ED] DT=0.2 (typo of dt)              raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 
[3 ED] n_t=20 (typo of nt)              raises TypeError: evolution_ABC() got an unexpected keyword argument 'n_t'. Did you mean
[3 ED] tevol_method='TEBD' as a kwarg   raises TypeError: evolution_ABC() got an unexpected keyword argument 'tevol_method'
[2 DMRG] DT=0.2 (typo of dt)              returned len=20, ts[1]=0.010, ts[-1]=0.19; <Sx0>(t=3.8) returned 0.4910, cos(Bt)/2 at the intended t = -0.3955
[2 DMRG] n_t=20 (typo of nt)              returned len=1000, ts[1]=0.200, ts[-1]=199.80; <Sx0>(t=3.8) returned -0.3851, cos(Bt)/2 at the intended t = -0.3955
[2 DMRG] tevol_method='TEBD' as a kwarg   returned len=20, ts[1]=0.200, ts[-1]=3.80; <Sx0>(t=3.8) returned -0.3851, cos(Bt)/2 at the intended t = -0.3955
[2 ED] DT=0.2 (typo of dt)              raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 
[2 ED] n_t=20 (typo of nt)              raises TypeError: evolution_ABC() got an unexpected keyword argument 'n_t'. Did you mean
[2 ED] tevol_method='TEBD' as a kwarg   raises TypeError: evolution_ABC() got an unexpected keyword argument 'tevol_method'
```


**Reviewer (CONFIRMED)**: reproduced on `"python"`, v3 and v2 alike:
`evolve_and_measure(mode="DMRG", nt=20, DT=0.2)` returns 20 points with `dt` 0.010
and `ts[-1]` = 0.19, whose last value 0.4910 equals cos(0.19)/2 to every printed
digit, against -0.3955 at the intended t=3.8; `evolution_ABA` gives the identical
result; the correctly spelled `dt=0.2` returns -0.3955 on `"python"` and v3 (-0.3851
on v2, its MPO-Taylor step error, not part of this finding); `mode="ED"` raises
`TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean
'dt'?` for both entry points. No hidden consumer, and no caller relies on the
swallow: every call in `src/`, `examples/` and `tests/` passes only named
parameters. Nothing struck; one correction to the origin, which is `f54b3c2`
(2020-02-21), not `6e16aab`, which only added `h=`; and one clarification the
hunter concedes, that 0.4910 is the exact value at the reported t=0.19, so the
defect is the dropped keyword and the wrong window, not a wrong number on the
returned grid.

```bash
cd <scratch>/reviews/realtime_2 && <scratch>/run3.sh 01_kwargs_review.py 2>&1 | tee 01_kwargs_review.out
```

`reviews/realtime_2/01_kwargs_review.py`:

```python
"""Reviewer probe for realtime_2: do evolve_and_measure/evolution_ABA on
mode="DMRG" drop a misspelled keyword and run at the defaults, while
mode="ED" raises?  Anchor: Larmor precession, H = B sum Sz from the +x
product state, <Sx_0>(t) = cos(Bt)/2 exactly.

Also probes: (a) the defaults each mode falls back to (DMRG nt=1000 vs ED
nt=100), (b) keywords the DMRG route documents (h=, return_wf=) on ED.
"""
import numpy as np
import dmrgpy
from dmrgpy import spinchain, timedependent, cppext
print("dmrgpy from", dmrgpy.__file__)

B, n = 1.0, 3
def chain(version, mode):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=version)
    sc.maxm, sc.nsweeps = 10, 10
    sc.set_hamiltonian(sum(-sc.Sx[i] for i in range(n)))
    wf = sc.get_gs(mode=mode)      # +x polarized
    sc.set_hamiltonian(B*sum(sc.Sz[i] for i in range(n)))
    return sc, wf

def show(tag, f):
    try:
        out = f()
    except Exception as exc:
        print("%-44s raises %s: %s" % (tag, type(exc).__name__, str(exc)[:80]))
        return
    ts, y = out[0], np.asarray(out[1])
    print("%-44s len=%d dt_used=%.3f ts[-1]=%.2f y[-1]=%.4f exact(ts[-1])=%.4f"
          " exact(0.2*(len-1))=%.4f" % (tag, len(y), ts[1]-ts[0], ts[-1],
          y[-1].real, np.cos(B*ts[-1])/2, np.cos(B*0.2*(len(y)-1))/2))

versions = ["python"] + ([3] if cppext.available(3) else []) + \
           ([2] if cppext.available(2) else [])
for v in versions:
    for mode in ("DMRG", "ED"):
        sc, wf = chain(v, mode)
        p = "[%s %s]" % (v, mode)
        show(p+" em nt=20 DT=0.2", lambda: timedependent.evolve_and_measure(
            sc, operator=sc.Sx[0], wf=wf, mode=mode, nt=20, DT=0.2))
        show(p+" em nt=20 dt=0.2 (correct)", lambda: timedependent.evolve_and_measure(
            sc, operator=sc.Sx[0], wf=wf, mode=mode, nt=20, dt=0.2))
        show(p+" ABA nt=20 DT=0.2", lambda: timedependent.evolution_ABA(
            sc, A=None, B=sc.Sx[0], wf=wf, mode=mode, nt=20, DT=0.2))
        show(p+" em defaults (no nt, no dt)", lambda: timedependent.evolve_and_measure(
            sc, operator=sc.Sx[0], wf=wf, mode=mode))
        show(p+" em h=H (documented on DMRG)", lambda: timedependent.evolve_and_measure(
            sc, operator=sc.Sx[0], wf=wf, mode=mode, nt=5, dt=0.2, h=sc.hamiltonian))
        show(p+" em return_wf=False", lambda: timedependent.evolve_and_measure(
            sc, operator=sc.Sx[0], wf=wf, mode=mode, nt=5, dt=0.2, return_wf=False))
```

`reviews/realtime_2/01_kwargs_review.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
[python DMRG] em nt=20 DT=0.2                len=20 dt_used=0.010 ts[-1]=0.19 y[-1]=0.4910 exact(ts[-1])=0.4910 exact(0.2*(len-1))=-0.3955
[python DMRG] em nt=20 dt=0.2 (correct)      len=20 dt_used=0.200 ts[-1]=3.80 y[-1]=-0.3955 exact(ts[-1])=-0.3955 exact(0.2*(len-1))=-0.3955
[python DMRG] ABA nt=20 DT=0.2               len=20 dt_used=0.010 ts[-1]=0.19 y[-1]=0.4910 exact(ts[-1])=0.4910 exact(0.2*(len-1))=-0.3955
[python DMRG] em defaults (no nt, no dt)     len=1000 dt_used=0.010 ts[-1]=9.99 y[-1]=-0.4222 exact(ts[-1])=-0.4222 exact(0.2*(len-1))=0.1520
[python DMRG] em h=H (documented on DMRG)    len=5 dt_used=0.200 ts[-1]=0.80 y[-1]=0.3484 exact(ts[-1])=0.3484 exact(0.2*(len-1))=0.3484
[python DMRG] em return_wf=False             len=5 dt_used=0.200 ts[-1]=0.80 y[-1]=0.3484 exact(ts[-1])=0.3484 exact(0.2*(len-1))=0.3484
[python ED] em nt=20 DT=0.2                  raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 'dt'?
[python ED] em nt=20 dt=0.2 (correct)        len=20 dt_used=0.200 ts[-1]=3.80 y[-1]=-0.3955 exact(ts[-1])=-0.3955 exact(0.2*(len-1))=-0.3955
[python ED] ABA nt=20 DT=0.2                 raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 'dt'?
[python ED] em defaults (no nt, no dt)       len=100 dt_used=0.010 ts[-1]=0.99 y[-1]=0.2743 exact(ts[-1])=0.2743 exact(0.2*(len-1))=0.2907
[python ED] em h=H (documented on DMRG)      raises TypeError: evolve_and_measure() got multiple values for argument 'h'
[python ED] em return_wf=False               raises TypeError: evolution_ABC() got an unexpected keyword argument 'return_wf'
[3 DMRG] em nt=20 DT=0.2                     len=20 dt_used=0.010 ts[-1]=0.19 y[-1]=0.4910 exact(ts[-1])=0.4910 exact(0.2*(len-1))=-0.3955
[3 DMRG] em nt=20 dt=0.2 (correct)           len=20 dt_used=0.200 ts[-1]=3.80 y[-1]=-0.3955 exact(ts[-1])=-0.3955 exact(0.2*(len-1))=-0.3955
[3 DMRG] ABA nt=20 DT=0.2                    len=20 dt_used=0.010 ts[-1]=0.19 y[-1]=0.4910 exact(ts[-1])=0.4910 exact(0.2*(len-1))=-0.3955
[3 DMRG] em defaults (no nt, no dt)          len=1000 dt_used=0.010 ts[-1]=9.99 y[-1]=-0.4222 exact(ts[-1])=-0.4222 exact(0.2*(len-1))=0.1520
[3 DMRG] em h=H (documented on DMRG)         len=5 dt_used=0.200 ts[-1]=0.80 y[-1]=0.3484 exact(ts[-1])=0.3484 exact(0.2*(len-1))=0.3484
[3 DMRG] em return_wf=False                  len=5 dt_used=0.200 ts[-1]=0.80 y[-1]=0.3484 exact(ts[-1])=0.3484 exact(0.2*(len-1))=0.3484
[3 ED] em nt=20 DT=0.2                       raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 'dt'?
[3 ED] em nt=20 dt=0.2 (correct)             len=20 dt_used=0.200 ts[-1]=3.80 y[-1]=-0.3955 exact(ts[-1])=-0.3955 exact(0.2*(len-1))=-0.3955
[3 ED] ABA nt=20 DT=0.2                      raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 'dt'?
[3 ED] em defaults (no nt, no dt)            len=100 dt_used=0.010 ts[-1]=0.99 y[-1]=0.2743 exact(ts[-1])=0.2743 exact(0.2*(len-1))=0.2907
[3 ED] em h=H (documented on DMRG)           raises TypeError: evolve_and_measure() got multiple values for argument 'h'
[3 ED] em return_wf=False                    raises TypeError: evolution_ABC() got an unexpected keyword argument 'return_wf'
[2 DMRG] em nt=20 DT=0.2                     len=20 dt_used=0.010 ts[-1]=0.19 y[-1]=0.4910 exact(ts[-1])=0.4910 exact(0.2*(len-1))=-0.3955
[2 DMRG] em nt=20 dt=0.2 (correct)           len=20 dt_used=0.200 ts[-1]=3.80 y[-1]=-0.3851 exact(ts[-1])=-0.3955 exact(0.2*(len-1))=-0.3955
[2 DMRG] ABA nt=20 DT=0.2                    len=20 dt_used=0.010 ts[-1]=0.19 y[-1]=0.4910 exact(ts[-1])=0.4910 exact(0.2*(len-1))=-0.3955
[2 DMRG] em defaults (no nt, no dt)          len=1000 dt_used=0.010 ts[-1]=9.99 y[-1]=-0.4222 exact(ts[-1])=-0.4222 exact(0.2*(len-1))=0.1520
[2 DMRG] em h=H (documented on DMRG)         len=5 dt_used=0.200 ts[-1]=0.80 y[-1]=0.3458 exact(ts[-1])=0.3484 exact(0.2*(len-1))=0.3484
[2 DMRG] em return_wf=False                  len=5 dt_used=0.200 ts[-1]=0.80 y[-1]=0.3458 exact(ts[-1])=0.3484 exact(0.2*(len-1))=0.3484
[2 ED] em nt=20 DT=0.2                       raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 'dt'?
[2 ED] em nt=20 dt=0.2 (correct)             len=20 dt_used=0.200 ts[-1]=3.80 y[-1]=-0.3955 exact(ts[-1])=-0.3955 exact(0.2*(len-1))=-0.3955
[2 ED] ABA nt=20 DT=0.2                      raises TypeError: evolution_ABC() got an unexpected keyword argument 'DT'. Did you mean 'dt'?
[2 ED] em defaults (no nt, no dt)            len=100 dt_used=0.010 ts[-1]=0.99 y[-1]=0.2743 exact(ts[-1])=0.2743 exact(0.2*(len-1))=0.2907
[2 ED] em h=H (documented on DMRG)           raises TypeError: evolve_and_measure() got multiple values for argument 'h'
[2 ED] em return_wf=False                    raises TypeError: evolution_ABC() got an unexpected keyword argument 'return_wf'
```

```bash
cd <scratch>/reviews/realtime_2 && <scratch>/run3.sh 02_h_keyword_on_ed.py 2>&1 | tee 02_h_keyword_on_ed.out
```

`reviews/realtime_2/02_h_keyword_on_ed.py`:

```python
"""Reviewer side probe: h= is a live keyword of evolve_and_measure_dmrg
(the evolution Hamiltonian, independent of the chain's own), but on
mode="ED" both dispatchers pass self.hamiltonian and then **kwargs, so h=
collides.  Anchor: chain Hamiltonian H1 = -sum Sx (its ground state is the
+x state), evolution Hamiltonian H2 = B sum Sz, <Sx_0>(t) = cos(Bt)/2.
"""
import numpy as np
from dmrgpy import spinchain, timedependent, cppext

B, n, nt, dt = 1.0, 3, 5, 0.2
exact = np.cos(B*dt*np.arange(nt))/2
for v in ["python"] + ([3] if cppext.available(3) else []):
    for mode in ("DMRG", "ED"):
        sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
        sc.maxm, sc.nsweeps = 10, 10
        sc.set_hamiltonian(sum(-sc.Sx[i] for i in range(n)))
        wf = sc.get_gs(mode=mode)
        H2 = B*sum(sc.Sz[i] for i in range(n))
        for name, f in [
            ("evolve_and_measure(h=H2)", lambda: timedependent.evolve_and_measure(
                sc, operator=sc.Sx[0], wf=wf, mode=mode, nt=nt, dt=dt, h=H2)),
            ("evolution_ABA(h=H2)", lambda: timedependent.evolution_ABA(
                sc, A=None, B=sc.Sx[0], wf=wf, mode=mode, nt=nt, dt=dt, h=H2))]:
            try:
                ts, y = f()
                print("[%s %s] %-26s max|y - cos(Bt)/2| = %.2e" % (
                    v, mode, name, np.max(np.abs(np.asarray(y) - exact))))
            except Exception as exc:
                print("[%s %s] %-26s raises %s: %s" % (
                    v, mode, name, type(exc).__name__, str(exc)[:80]))
```

`reviews/realtime_2/02_h_keyword_on_ed.out`:


```
[python DMRG] evolve_and_measure(h=H2)   max|y - cos(Bt)/2| = 1.83e-15
[python DMRG] evolution_ABA(h=H2)        max|y - cos(Bt)/2| = 3.55e-15
[python ED] evolve_and_measure(h=H2)   raises TypeError: evolve_and_measure() got multiple values for argument 'h'
[python ED] evolution_ABA(h=H2)        raises TypeError: dmrgpy.edtk.timedependent.evolution_ABA() got multiple values for keyword argume
[3 DMRG] evolve_and_measure(h=H2)   max|y - cos(Bt)/2| = 2.22e-16
[3 DMRG] evolution_ABA(h=H2)        max|y - cos(Bt)/2| = 2.22e-16
[3 ED] evolve_and_measure(h=H2)   raises TypeError: evolve_and_measure() got multiple values for argument 'h'
[3 ED] evolution_ABA(h=H2)        raises TypeError: dmrgpy.edtk.timedependent.evolution_ABA() got multiple values for keyword argume
```


**Suggested fix**: drop `**kwargs` from `evolve_and_measure_dmrg` in
`timedependent.py` and in `mpsjulialive/timedependent.py`; the `julia_live`
dispatch at `timedependent.py:262` already names every argument, so it only loses
its trailing `**kwargs`. `evolution_ABA` can keep its own, since both of its
consumers then reject unknown names. ED's keyword set is not the reference for which
keywords are valid, since ED rejects `h=`, which DMRG honours (finding 16), so the
two are best fixed together. No returned number changes for a correctly spelled
call; a misspelled one raises.

### 16. On `mode="ED"`, `evolve_and_measure` and `evolution_ABA` raise `TypeError` ("got multiple values for argument 'h'") on `h=`, a keyword the DMRG route accepts and honours, because both ED dispatchers pass `self.hamiltonian` next to a forwarded `**kwargs`; and the default `nt` is 1000 on DMRG against 100 on ED

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `realtime` (turned up by the reviewer of finding 15)

**Where**: `src/dmrgpy/timedependent.py:208-211` (`evolve_and_measure`'s ED branch
passes `h` positionally and then `**kwargs`) and `:318-321` (`evolution_ABA`'s ED
branch passes `h=self.hamiltonian` and then `**kwargs`); `edtk/timedependent.py:8`
(`evolution_ABC`'s `nt=100`) against `timedependent.py:215` (`nt=1000`).

On DMRG, `h=H2` evolves under H2 exactly, and it is forwarded to `julia_live` too;
on ED the same call dies on the collision. The ED layer itself honours an arbitrary
`h` (called directly with `h=H2` while the chain still holds H1, it is exact to
5.54e-09), so the defect sits only in the two dispatcher lines. The ED dispatch
with a fixed `h` dates from `f54b3c2` (2020-02-21), and the collision became
reachable when `6e16aab` (2020-08-25) added `h=None` to `evolve_and_measure_dmrg`.
It is the section 4.10 shape again: the dispatcher fixes a value and then forwards
a `**kwargs` that can carry the same name. It survived because no test, example or
document passes `h=` to either function; the documented quench recipe swaps the
chain's Hamiltonian with `set_hamiltonian` instead, and that works on both modes.

**Expected**: `evolve_and_measure(mode="ED", h=H2)` evolving under H2 as DMRG does,
and one default `nt` on both modes.

Repro, from the reviewer of finding 15, who found it (both scripts spliced under
finding 15):

`reviews/realtime_2/02_h_keyword_on_ed.py` and its output are spliced under finding 15 above.

`reviews/realtime_2/01_kwargs_review.py` and its output are spliced under finding 15 above.


**Reviewer (CONFIRMED, NARROWED)**: rerun on `"python"`, v3 and v2 with H1 = -sum Sx,
so the start is the +x eigenstate, and H2 = sum Sz, closed form cos(Bt)/2 at
`nt=5`, `dt=0.2`. On DMRG with `h=H2` the error is 6.68e-16 and 8.33e-16 on
`"python"`, 1.67e-16 on v3 and 2.57e-03 on v2 (its MPO-Taylor stepper), and without
`h` <Sx_0> stays at 0.5 to 3e-16, so `h=` is live on DMRG and not ignored. On ED
both functions raise on all three backends. The `set_hamiltonian(H2)` workaround is
exact to 5.54e-09 on ED. At default arguments DMRG returns 1000 points to `ts[-1]` =
9.99 and ED 100 points to 0.99. Struck, each explicitly:

- "The evolution Hamiltonian keyword the DMRG route documents": `h=` is in the
  signature and honoured, but no docstring, user guide, architecture document or
  README mentions it; it is accepted and honoured, not documented.
- The observation that ED also rejects `return_wf=False`: by its own docstring
  `return_wf` wraps a `cpp_handle` and is DMRG-only, so a named `TypeError` there
  is a documented boundary.
- "The ED reference cannot be run for a quench from one Hamiltonian's ground state
  under another" is kept only with its own qualifier: the swap is the documented
  route and it works, so this is an inconvenience, not a missing capability.

```bash
cd <scratch>/reviews/realtime_2_new1 && <scratch>/run3.sh 01_h_keyword_ed.py 2>&1 | tee 01_h_keyword_ed.out
```

`reviews/realtime_2_new1/01_h_keyword_ed.py`:

```python
"""Reviewer repro of realtime_2_new1.  Chain H1 = -sum Sx (GS = +x state),
evolution Hamiltonian H2 = B sum Sz, <Sx_0>(t) = cos(Bt)/2; under H1 the
+x state is an eigenstate, so <Sx_0> would stay 0.5 (discriminating).
Also: (a) the ED workaround of swapping the chain Hamiltonian, (b) default
nt on both routes, (c) whether ED's own evolution_ABC signature ever had h
as a free keyword (i.e. whether any ED caller could pass h=)."""
import inspect
import numpy as np
import dmrgpy
from dmrgpy import spinchain, timedependent, cppext
from dmrgpy.edtk import timedependent as tded
print("dmrgpy from", dmrgpy.__file__)
print("ED evolution_ABC signature:", inspect.signature(tded.evolution_ABC))
print("DMRG evolve_and_measure_dmrg signature:",
      inspect.signature(timedependent.evolve_and_measure_dmrg))

B, n, nt, dt = 1.0, 3, 5, 0.2
ts0 = dt*np.arange(nt)
exact = np.cos(B*ts0)/2
def build(v):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    sc.maxm, sc.nsweeps = 10, 10
    sc.set_hamiltonian(sum(-sc.Sx[i] for i in range(n)))
    return sc
for v in ["python"] + ([3] if cppext.available(3) else []) + ([2] if cppext.available(2) else []):
    for mode in ("DMRG", "ED"):
        sc = build(v)
        wf = sc.get_gs(mode=mode)
        H2 = B*sum(sc.Sz[i] for i in range(n))
        calls = [
            ("evolve_and_measure(h=H2)", lambda: timedependent.evolve_and_measure(
                sc, operator=sc.Sx[0], wf=wf, mode=mode, nt=nt, dt=dt, h=H2)),
            ("evolution_ABA(h=H2)", lambda: timedependent.evolution_ABA(
                sc, A=None, B=sc.Sx[0], wf=wf, mode=mode, nt=nt, dt=dt, h=H2)),
            ("evolve_and_measure(no h)", lambda: timedependent.evolve_and_measure(
                sc, operator=sc.Sx[0], wf=wf, mode=mode, nt=nt, dt=dt))]
        for name, f in calls:
            try:
                ts, y = f()
                print("[%s %s] %-26s max|y-cos(Bt)/2|=%.2e max|y-0.5|=%.2e" % (
                    v, mode, name, np.max(np.abs(np.asarray(y)-exact)),
                    np.max(np.abs(np.asarray(y)-0.5))))
            except Exception as exc:
                print("[%s %s] %-26s raises %s: %s" % (
                    v, mode, name, type(exc).__name__, str(exc)[:90]))
        # (a) workaround: keep the H1 ground state, swap the chain Hamiltonian
        sc.set_hamiltonian(H2)
        try:
            ts, y = timedependent.evolve_and_measure(
                sc, operator=sc.Sx[0], wf=wf, mode=mode, nt=nt, dt=dt)
            print("[%s %s] %-26s max|y-cos(Bt)/2|=%.2e" % (
                v, mode, "swap H then evolve", np.max(np.abs(np.asarray(y)-exact))))
        except Exception as exc:
            print("[%s %s] swap H then evolve raises %s: %s" % (
                v, mode, type(exc).__name__, str(exc)[:90]))
    # (b) default nt on both routes
    for mode in ("DMRG", "ED"):
        sc = build(v)
        wf = sc.get_gs(mode=mode)
        ts, y = timedependent.evolve_and_measure(sc, operator=sc.Sx[0], wf=wf, mode=mode)
        ts2, y2 = timedependent.evolution_ABA(sc, B=sc.Sx[0], wf=wf, mode=mode)
        print("[%s %s] defaults: em len=%d ts[-1]=%.2f  ABA len=%d ts[-1]=%.2f" % (
            v, mode, len(ts), ts[-1], len(ts2), ts2[-1]))
```

`reviews/realtime_2_new1/01_h_keyword_ed.out`:


```
dmrgpy from <repo>/src/dmrgpy/__init__.py
ED evolution_ABC signature: (self, h, A=None, B=None, C=None, wf=None, nt=100, dt=0.01)
DMRG evolve_and_measure_dmrg signature: (self, operator=None, nt=1000, h=None, dt=0.01, wf=None, return_wf=False, **kwargs)
[python DMRG] evolve_and_measure(h=H2)   max|y-cos(Bt)/2|=6.68e-16 max|y-0.5|=1.52e-01
[python DMRG] evolution_ABA(h=H2)        max|y-cos(Bt)/2|=8.33e-16 max|y-0.5|=1.52e-01
[python DMRG] evolve_and_measure(no h)   max|y-cos(Bt)/2|=1.52e-01 max|y-0.5|=6.67e-16
[python DMRG] swap H then evolve         max|y-cos(Bt)/2|=6.68e-16
[python ED] evolve_and_measure(h=H2)   raises TypeError: evolve_and_measure() got multiple values for argument 'h'
[python ED] evolution_ABA(h=H2)        raises TypeError: dmrgpy.edtk.timedependent.evolution_ABA() got multiple values for keyword argument 'h'
[python ED] evolve_and_measure(no h)   max|y-cos(Bt)/2|=1.52e-01 max|y-0.5|=3.29e-08
[python ED] swap H then evolve         max|y-cos(Bt)/2|=5.54e-09
[python DMRG] defaults: em len=1000 ts[-1]=9.99  ABA len=1000 ts[-1]=9.99
[python ED] defaults: em len=100 ts[-1]=0.99  ABA len=100 ts[-1]=0.99
[3 DMRG] evolve_and_measure(h=H2)   max|y-cos(Bt)/2|=1.67e-16 max|y-0.5|=1.52e-01
[3 DMRG] evolution_ABA(h=H2)        max|y-cos(Bt)/2|=1.67e-16 max|y-0.5|=1.52e-01
[3 DMRG] evolve_and_measure(no h)   max|y-cos(Bt)/2|=1.52e-01 max|y-0.5|=3.33e-16
[3 DMRG] swap H then evolve         max|y-cos(Bt)/2|=1.67e-16
[3 ED] evolve_and_measure(h=H2)   raises TypeError: evolve_and_measure() got multiple values for argument 'h'
[3 ED] evolution_ABA(h=H2)        raises TypeError: dmrgpy.edtk.timedependent.evolution_ABA() got multiple values for keyword argument 'h'
[3 ED] evolve_and_measure(no h)   max|y-cos(Bt)/2|=1.52e-01 max|y-0.5|=3.29e-08
[3 ED] swap H then evolve         max|y-cos(Bt)/2|=5.54e-09
[3 DMRG] defaults: em len=1000 ts[-1]=9.99  ABA len=1000 ts[-1]=9.99
[3 ED] defaults: em len=100 ts[-1]=0.99  ABA len=100 ts[-1]=0.99
[2 DMRG] evolve_and_measure(h=H2)   max|y-cos(Bt)/2|=2.57e-03 max|y-0.5|=1.54e-01
[2 DMRG] evolution_ABA(h=H2)        max|y-cos(Bt)/2|=2.57e-03 max|y-0.5|=1.54e-01
[2 DMRG] evolve_and_measure(no h)   max|y-cos(Bt)/2|=1.59e-01 max|y-0.5|=7.71e-03
[2 DMRG] swap H then evolve         max|y-cos(Bt)/2|=2.57e-03
[2 ED] evolve_and_measure(h=H2)   raises TypeError: evolve_and_measure() got multiple values for argument 'h'
[2 ED] evolution_ABA(h=H2)        raises TypeError: dmrgpy.edtk.timedependent.evolution_ABA() got multiple values for keyword argument 'h'
[2 ED] evolve_and_measure(no h)   max|y-cos(Bt)/2|=1.52e-01 max|y-0.5|=3.29e-08
[2 ED] swap H then evolve         max|y-cos(Bt)/2|=5.54e-09
[2 DMRG] defaults: em len=1000 ts[-1]=9.99  ABA len=1000 ts[-1]=9.99
[2 ED] defaults: em len=100 ts[-1]=0.99  ABA len=100 ts[-1]=0.99
```

```bash
cd <scratch>/reviews/realtime_2_new1 && <scratch>/run3.sh 02_ed_layer_takes_h.py 2>&1 | tee 02_ed_layer_takes_h.out
```

`reviews/realtime_2_new1/02_ed_layer_takes_h.py`:

```python
# A=None replaced by identity: the top-level evolution_ABA converts None, tded does not (first run crashed on that, a probe error)
"""Does the ED layer itself honour an arbitrary evolution Hamiltonian, so
that the fix is only in the dispatcher?  Call edtk.timedependent directly
with h=H2 while the chain (and its EDchain) still hold H1."""
import numpy as np
from dmrgpy import spinchain, multioperator
from dmrgpy.edtk import timedependent as tded
B, n, nt, dt = 1.0, 3, 5, 0.2
exact = np.cos(B*dt*np.arange(nt))/2
sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
sc.set_hamiltonian(sum(-sc.Sx[i] for i in range(n)))
wf = sc.get_gs(mode="ED")
H2 = B*sum(sc.Sz[i] for i in range(n))
ed = sc.get_ED_obj()
ts, y = tded.evolve_and_measure(ed, H2, operator=sc.Sx[0], wf=wf, nt=nt, dt=dt)
print("tded.evolve_and_measure(ed, H2)  max|y-cos(Bt)/2| = %.2e" % np.max(np.abs(y-exact)))
ts, y = tded.evolution_ABA(ed, h=H2, A=multioperator.identity(), B=sc.Sx[0], wf=wf, nt=nt, dt=dt)
print("tded.evolution_ABA(ed, h=H2)     max|y-cos(Bt)/2| = %.2e" % np.max(np.abs(y-exact)))
ts, y = tded.evolve_and_measure(ed, sc.hamiltonian, operator=sc.Sx[0], wf=wf, nt=nt, dt=dt)
print("tded.evolve_and_measure(ed, H1)  max|y-0.5|       = %.2e" % np.max(np.abs(y-0.5)))
```

`reviews/realtime_2_new1/02_ed_layer_takes_h.out`:


```
tded.evolve_and_measure(ed, H2)  max|y-cos(Bt)/2| = 5.54e-09
tded.evolution_ABA(ed, h=H2)     max|y-cos(Bt)/2| = 5.54e-09
tded.evolve_and_measure(ed, H1)  max|y-0.5|       = 3.29e-08
```


**Suggested fix**: in the two ED branches of `timedependent.py`, `h =
kwargs.pop("h", None)`, then `if h is None: h = self.hamiltonian`, then forward.
The explicit `None` test matters, since on DMRG `h=None` means the chain's
Hamiltonian and ED should read it the same way, and a `MultiOperator`'s truthiness
is not that test. For `nt`, put `nt=1000, dt=1e-2` on the two public dispatchers
and forward them explicitly to both modes rather than changing `evolution_ABC`'s
default, since `evolution_DC` keeps its own `nt=100` and does not go through
`evolution_ABC`. Pins: `h=H2` on ED against cos(Bt)/2 and against DMRG, and equal
`len(ts)` across modes at defaults. Lands with finding 15. No number changes for a
call that worked; a default-argument ED call would return the longer trajectory.

### 17. The ED real-time propagator is `solve_ivp` RK45 at scipy's default rtol=1e-3/atol=1e-6, one call per dt, neither exact nor unitary, with an error set by dt times the absolute energy of the evolved state: a constant +20 added to H moves `evolve_and_measure(mode="ED")` from 2.1e-07 to 5.6e-04 off exact at the tests' own dt=0.1, and at dt=0.2 on an 8-site chain ED is 1.2 per cent off where `"python"` TDVP is at 7e-9

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `realtime`

**Where**: `src/dmrgpy/edtk/tdtk.py:26-33` (`solve_ivp(f, tspan, v0,
method="RK45", t_eval=...)` with default tolerances); consumers
`edtk/timedependent.py:54-55` (`evolution_ABC`, behind `evolve_and_measure` and
`evolution_ABA` on ED) and `:128` (`evolution_DC`, `submode="TD"` on ED).

ED is the correctness reference of the ED-versus-DMRG real-time tests, and its
propagator is an adaptive Runge-Kutta at scipy's default tolerances. What RK45
resolves is the absolute frequency of the state, so its error scales with dt times
|E|, not with the bandwidth: `evolution_ABC` integrates the unshifted H, so a
constant offset, which is physically inert, degrades it, while `evolution_DC`
integrates H - E0 and is correspondingly less affected. A chemical potential or a
Hubbard U supplies such an offset without anyone adding one on purpose. `867e2b4`
made `evolution_ABC` call it with `-Hop` and left the tolerance alone; the
propagator dates from `f54b3c2` (2020-02-21). It survived because the tests use
small, low-energy chains at dt of 0.05 to 0.1, where RK45 lands at 1e-8 to 1e-6,
below the tolerances they assert, and no test moves the energy origin; the
2026-09-24 record's reviewer called ED "exact at this level rather than exact by
construction" in passing, and neither it nor the second pass made a finding of it.

**Expected**: the exact e^{-iHt} from a dense eigendecomposition of ED's own H,
which is independent of both backends and invariant under a constant offset of H.

Repro, from the hunter:

```bash
cd <scratch>/realtime && <scratch>/run3.sh 03_ed_rk45_accuracy.py 2>&1 | tee 03_ed_rk45_accuracy.out
```

`realtime/03_ed_rk45_accuracy.py`:

```python
"""Lens realtime, lead 3: ED's propagator is solve_ivp RK45 at scipy's
default rtol=1e-3/atol=1e-6, one call per dt. How far is ED, the reference
every DMRG evolution test is held to, from the exact e^{-iHt}, as dt and
the spectral width grow?

Anchor: e^{-iHt} by dense eigendecomposition of ED's own H matrix.
Two ED routes: evolve_and_measure (evolution_ABC) and evolution_DC (the
series behind submode="TD"), each at its caller's dt.
"""
import numpy as np
from dmrgpy import spinchain, timedependent

def chain(n, J, hz):
    sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        h = h + hz*(0.3 + 0.1*i)*sc.Sz[i] + 0.4*sc.Sx[i]
    sc.set_hamiltonian(h)
    return sc

for (n, J, hz) in [(8, 1.0, 1.0), (8, 3.0, 3.0)]:
    sc = chain(n, J, hz)
    ed = sc.get_ED_obj()
    H = ed.get_operator(sc.hamiltonian).toarray()
    w, v = np.linalg.eigh(H)
    print("n=%d J=%.1f hz=%.1f  bandwidth=%.2f" % (n, J, hz, w[-1] - w[0]))
    g = v[:, 0]
    # evolve_and_measure from a product-like start: Sx_0 |GS> normalized
    wf = sc.get_gs(mode="ED")
    A = ed.get_operator(sc.Sz[0]).toarray()
    Bop = ed.get_operator(sc.Sz[3]).toarray()
    for dt, nt in [(0.05, 200), (0.1, 100), (0.2, 50), (0.5, 20)]:
        ts = dt*np.arange(nt)
        # exact, same unnormalized convention as evolution_ABA
        psi0 = A@ed.get_gs_array()
        c = v.conj().T@psi0
        ref = np.array([np.vdot(v@(np.exp(-1j*w*t)*c), Bop@(v@(np.exp(-1j*w*t)*c)))
                        for t in ts])
        _t, y = timedependent.evolution_ABA(sc, A=sc.Sz[0], B=sc.Sz[3],
                                             mode="ED", nt=nt, dt=dt)
        y = np.asarray(y)
        print("  evolution_ABA ED dt=%.2f T=%.1f  max|y-exact|=%.2e  "
              "(scale %.3f)" % (dt, dt*nt, np.max(np.abs(y - ref)),
                                np.max(np.abs(ref))))
    # evolution_DC, the TD series, at TD's default dt=0.1 and nt=600
    for dt, nt in [(0.1, 600), (0.2, 300)]:
        ts = dt*np.arange(nt)
        _t, cs = timedependent.evolution_DC(sc, mode="ED", name=(sc.Sz[3], sc.Sz[0]),
                                            nt=nt, dt=dt)
        cs = np.asarray(cs)
        e0 = w[0]
        # <GS|A e^{+i(H-E0)t} B|GS> for name=(A,B)=(Sz3,Sz0)
        gA = v.conj().T@(Bop@g)   # <n|Sz3|GS>^*... build explicitly below
        a_n = (v.conj().T@(Bop.conj().T@g)).conj()  # <GS|Sz3|n>
        b_n = v.conj().T@(A@g)                       # <n|Sz0|GS>
        ref = np.array([np.sum(a_n*b_n*np.exp(1j*(w - e0)*t)) for t in ts])
        print("  evolution_DC ED  dt=%.2f T=%.0f  max|cs-exact|=%.2e  "
              "(scale %.3f)  |cs(T)-exact(T)|=%.2e"
              % (dt, dt*nt, np.max(np.abs(cs - ref)), np.max(np.abs(ref)),
                 abs(cs[-1] - ref[-1])))
```

`realtime/03_ed_rk45_accuracy.out`:


```
n=8 J=1.0 hz=1.0  bandwidth=8.63
  evolution_ABA ED dt=0.05 T=10.0  max|y-exact|=2.55e-08  (scale 0.031)
  evolution_ABA ED dt=0.10 T=10.0  max|y-exact|=7.56e-08  (scale 0.031)
  evolution_ABA ED dt=0.20 T=10.0  max|y-exact|=2.69e-06  (scale 0.030)
  evolution_ABA ED dt=0.50 T=10.0  max|y-exact|=6.68e-05  (scale 0.030)
  evolution_DC ED  dt=0.10 T=60  max|cs-exact|=2.86e-08  (scale 0.105)  |cs(T)-exact(T)|=2.03e-08
  evolution_DC ED  dt=0.20 T=60  max|cs-exact|=1.07e-07  (scale 0.105)  |cs(T)-exact(T)|=7.05e-08
n=8 J=3.0 hz=3.0  bandwidth=24.30
  evolution_ABA ED dt=0.05 T=10.0  max|y-exact|=1.86e-05  (scale 0.069)
  evolution_ABA ED dt=0.10 T=10.0  max|y-exact|=4.64e-05  (scale 0.068)
  evolution_ABA ED dt=0.20 T=10.0  max|y-exact|=8.46e-04  (scale 0.068)
  evolution_ABA ED dt=0.50 T=10.0  max|y-exact|=1.06e-03  (scale 0.061)
  evolution_DC ED  dt=0.10 T=60  max|cs-exact|=3.09e-06  (scale 0.106)  |cs(T)-exact(T)|=2.21e-06
  evolution_DC ED  dt=0.20 T=60  max|cs-exact|=5.50e-05  (scale 0.106)  |cs(T)-exact(T)|=3.22e-05
```


**Reviewer (CONFIRMED, NARROWED)**: the hunter's script reruns to identical numbers
(the route is deterministic): on the bandwidth-24.30 chain `evolution_ABA` on ED is
1.86e-05, 4.64e-05, 8.46e-04 and 1.06e-03 off at dt = 0.05, 0.1, 0.2 and 0.5, and
`evolution_DC` 3.09e-06 and 5.50e-05 at dt = 0.1 and 0.2. Calling `tdtk.evolve`
directly on H + c*1 at a fixed bandwidth of 8.63, from a seeded random state over
T=10, the state error at dt=0.1 is 1.88e-06 at c=-0.52 (spectrum centred), 2.04e-06
at c=0, 2.49e-02 at c=+20 and 1.55e-01 at c=+60, with a norm error of -3.0e-03 at
c=+20. Through the public route, `evolve_and_measure(mode="ED")` on <Sz_2>(t) at
dt=0.1 with h + c*identity (the assembled matrices checked to differ by exactly
c*1) is 2.12e-07 off at c=0, 1.71e-04 at c=10 and 5.64e-04 at c=20, the last two
above the abs=1e-4 that `test_time_evolution.py` asserts, at a dt the tests use. At
`evolution_ABA`'s default dt=0.01 the error is 1.28e-08 with a norm drift of
3.2e-07. At dt=0.2 the relative norm drift is 8.97e-04, much smaller than the 1.2
per cent observable error, so phase error dominates. `"python"` TDVP at `maxm=32`
(exact for n=8) with dt=0.2 is 7.32e-09 off where ED is 8.46e-04; at n=10 ED is
1.97e-04 off at dt=0.1 and 1.11e-03 at dt=0.2. Struck, each explicitly:

- The error "grows with dt times bandwidth": at fixed bandwidth a constant shift
  takes it from 2e-06 to 2.5e-02, so the controlling scale is dt times the absolute
  energy, |E0|-like for `evolve_and_measure`/`evolution_ABA` and bandwidth-like only
  for `evolution_DC`.
- The hunter's `why_tests_pass`, "no test varies dt times bandwidth": what no test
  varies is the energy origin.
- "Every real-time test is anchored to ED": several compare v3 with `"python"` or
  with closed forms instead.
- "ED reproduces e^{-iHt} to machine precision" as a documented expectation: nothing
  documents it; the finding stands on the size of the error and on its dependence on
  the energy origin.

```bash
cd <scratch>/reviews/realtime_4 && <scratch>/run3.sh 01_hunter_repro.py 2>&1 | tee 01_hunter_repro.out
```

`reviews/realtime_4/01_hunter_repro.py` is `realtime/03_ed_rk45_accuracy.py` (finding 17) unchanged, rerun; its output, `reviews/realtime_4/01_hunter_repro.out`:


```
n=8 J=1.0 hz=1.0  bandwidth=8.63
  evolution_ABA ED dt=0.05 T=10.0  max|y-exact|=2.55e-08  (scale 0.031)
  evolution_ABA ED dt=0.10 T=10.0  max|y-exact|=7.56e-08  (scale 0.031)
  evolution_ABA ED dt=0.20 T=10.0  max|y-exact|=2.69e-06  (scale 0.030)
  evolution_ABA ED dt=0.50 T=10.0  max|y-exact|=6.68e-05  (scale 0.030)
  evolution_DC ED  dt=0.10 T=60  max|cs-exact|=2.86e-08  (scale 0.105)  |cs(T)-exact(T)|=2.03e-08
  evolution_DC ED  dt=0.20 T=60  max|cs-exact|=1.07e-07  (scale 0.105)  |cs(T)-exact(T)|=7.05e-08
n=8 J=3.0 hz=3.0  bandwidth=24.30
  evolution_ABA ED dt=0.05 T=10.0  max|y-exact|=1.86e-05  (scale 0.069)
  evolution_ABA ED dt=0.10 T=10.0  max|y-exact|=4.64e-05  (scale 0.068)
  evolution_ABA ED dt=0.20 T=10.0  max|y-exact|=8.46e-04  (scale 0.068)
  evolution_ABA ED dt=0.50 T=10.0  max|y-exact|=1.06e-03  (scale 0.061)
  evolution_DC ED  dt=0.10 T=60  max|cs-exact|=3.09e-06  (scale 0.106)  |cs(T)-exact(T)|=2.21e-06
  evolution_DC ED  dt=0.20 T=60  max|cs-exact|=5.50e-05  (scale 0.106)  |cs(T)-exact(T)|=3.22e-05
```

```bash
cd <scratch>/reviews/realtime_4 && <scratch>/run3.sh 02_mechanism_and_defaults.py 2>&1 | tee 02_mechanism_and_defaults.out
```

`reviews/realtime_4/02_mechanism_and_defaults.py`:

```python
"""Reviewer probe for realtime_4.

(a) Mechanism: is the RK45 error set by the bandwidth or by max|E|*dt?
    Call tdtk.evolve directly on H + c*I at fixed bandwidth; a constant c
    changes only a global phase of the exact propagator, so any dependence
    of the error on c is the solver resolving an absolute frequency.
(b) Size at the public defaults: evolution_ABA(mode="ED") at its default
    dt=0.01 on the bandwidth-24 chain, and the norm drift <psi(t)|psi(t)>.
(c) Dimension: the same chain at n=10.
(d) Is the DMRG backend it anchors more or less accurate? "python" TDVP at
    full bond dimension (maxm=32 >= 2^(8/2)) on the same case, dt=0.2.
"""
import numpy as np
from dmrgpy import spinchain, timedependent, multioperator
from dmrgpy.edtk.tdtk import evolve

def chain(n, J, hz):
    sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        h = h + hz*(0.3 + 0.1*i)*sc.Sz[i] + 0.4*sc.Sx[i]
    sc.set_hamiltonian(h)
    return sc

def exact_aba(ed, w, v, A, Bop, ts):
    psi0 = A@ed.get_gs_array()
    c = v.conj().T@psi0
    return np.array([np.vdot(v@(np.exp(-1j*w*t)*c), Bop@(v@(np.exp(-1j*w*t)*c)))
                     for t in ts])

# (a) shift test on the low-bandwidth chain (bandwidth 8.63)
sc = chain(8, 1.0, 1.0)
ed = sc.get_ED_obj()
Hs = ed.get_operator(sc.hamiltonian).tocsc()
H = Hs.toarray()
w, v = np.linalg.eigh(H)
print("(a) n=8 J=hz=1  E0=%.3f Emax=%.3f bandwidth=%.2f" % (w[0], w[-1], w[-1]-w[0]))
rng = np.random.default_rng(1)
psi = rng.normal(size=H.shape[0]) + 1j*rng.normal(size=H.shape[0])
psi /= np.linalg.norm(psi)
from scipy.sparse import identity
for c in [-0.5*(w[0]+w[-1]), 0.0, 20.0, 60.0]:
    Hc = Hs + c*identity(H.shape[0], format="csc")
    for dt, nt in [(0.1, 100), (0.2, 50)]:
        x = psi.copy()
        for _ in range(nt):
            x = evolve(x, -Hc, t=dt, dt=dt)
        ex = v@(np.exp(-1j*(w + c)*dt*nt)*(v.conj().T@psi))
        print("   shift c=%+7.2f  max|E+c|=%6.2f  dt=%.1f T=10  |psi_RK-psi_exact|=%.2e  norm-1=%+.2e"
              % (c, max(abs(w[0]+c), abs(w[-1]+c)), dt, np.linalg.norm(x-ex),
                 np.linalg.norm(x)-1))

# (b) defaults on the bandwidth-24 chain
sc = chain(8, 3.0, 3.0)
ed = sc.get_ED_obj()
H = ed.get_operator(sc.hamiltonian).toarray()
w, v = np.linalg.eigh(H)
print("(b) n=8 J=hz=3  E0=%.3f Emax=%.3f bandwidth=%.2f" % (w[0], w[-1], w[-1]-w[0]))
A = ed.get_operator(sc.Sz[0]).toarray()
Bop = ed.get_operator(sc.Sz[3]).toarray()
for dt, nt in [(0.01, 1000), (0.2, 50)]:
    ts = dt*np.arange(nt)
    ref = exact_aba(ed, w, v, A, Bop, ts)
    _t, y = timedependent.evolution_ABA(sc, A=sc.Sz[0], B=sc.Sz[3], mode="ED", nt=nt, dt=dt)
    y = np.asarray(y)
    _t, nrm = timedependent.evolution_ABA(sc, A=sc.Sz[0], B=multioperator.identity(),
                                          mode="ED", nt=nt, dt=dt)
    nrm = np.asarray(nrm).real
    print("   evolution_ABA ED dt=%.2f T=10  max|y-exact|=%.2e scale %.3f  "
          "norm drift max|<psi|psi>(t)/<psi|psi>(0)-1|=%.2e"
          % (dt, np.max(np.abs(y-ref)), np.max(np.abs(ref)), np.max(np.abs(nrm/nrm[0]-1))))

# (d) python TDVP at full bond dimension, same case, dt=0.2
sc.maxm = 32
sc.tevol_method = "TDVP"
wf = sc.get_gs()
dt, nt = 0.2, 50
ts = dt*np.arange(nt)
ref = exact_aba(ed, w, v, A, Bop, ts)
_t, yd = timedependent.evolution_ABA(sc, A=sc.Sz[0], B=sc.Sz[3], mode="DMRG", nt=nt, dt=dt)
yd = np.asarray(yd)
print("(d) python TDVP maxm=32 dt=0.20  max|y-exact|=%.2e  (DMRG GS energy err %.1e)"
      % (np.max(np.abs(yd-ref)), abs(sc.gs_energy()-w[0])))

# (c) n=10
sc = chain(10, 3.0, 3.0)
ed = sc.get_ED_obj()
H = ed.get_operator(sc.hamiltonian).toarray()
w, v = np.linalg.eigh(H)
print("(c) n=10 J=hz=3  E0=%.3f Emax=%.3f bandwidth=%.2f" % (w[0], w[-1], w[-1]-w[0]))
A = ed.get_operator(sc.Sz[0]).toarray()
Bop = ed.get_operator(sc.Sz[3]).toarray()
for dt, nt in [(0.1, 100), (0.2, 50)]:
    ts = dt*np.arange(nt)
    ref = exact_aba(ed, w, v, A, Bop, ts)
    _t, y = timedependent.evolution_ABA(sc, A=sc.Sz[0], B=sc.Sz[3], mode="ED", nt=nt, dt=dt)
    y = np.asarray(y)
    print("   evolution_ABA ED dt=%.2f T=10  max|y-exact|=%.2e scale %.3f"
          % (dt, np.max(np.abs(y-ref)), np.max(np.abs(ref))))
```

`reviews/realtime_4/02_mechanism_and_defaults.out`:


```
(a) n=8 J=hz=1  E0=-3.791 Emax=4.840 bandwidth=8.63
   shift c=  -0.52  max|E+c|=  4.32  dt=0.1 T=10  |psi_RK-psi_exact|=1.88e-06  norm-1=-5.22e-07
   shift c=  -0.52  max|E+c|=  4.32  dt=0.2 T=10  |psi_RK-psi_exact|=7.05e-05  norm-1=-1.47e-05
   shift c=  +0.00  max|E+c|=  4.84  dt=0.1 T=10  |psi_RK-psi_exact|=2.04e-06  norm-1=-4.43e-07
   shift c=  +0.00  max|E+c|=  4.84  dt=0.2 T=10  |psi_RK-psi_exact|=6.33e-05  norm-1=-1.05e-05
   shift c= +20.00  max|E+c|= 24.84  dt=0.1 T=10  |psi_RK-psi_exact|=2.49e-02  norm-1=-2.97e-03
   shift c= +20.00  max|E+c|= 24.84  dt=0.2 T=10  |psi_RK-psi_exact|=3.75e-02  norm-1=-3.93e-03
   shift c= +60.00  max|E+c|= 64.84  dt=0.1 T=10  |psi_RK-psi_exact|=1.55e-01  norm-1=-9.28e-03
   shift c= +60.00  max|E+c|= 64.84  dt=0.2 T=10  |psi_RK-psi_exact|=1.52e-01  norm-1=-1.13e-02
(b) n=8 J=hz=3  E0=-11.070 Emax=13.233 bandwidth=24.30
   evolution_ABA ED dt=0.01 T=10  max|y-exact|=1.28e-08 scale 0.069  norm drift max|<psi|psi>(t)/<psi|psi>(0)-1|=3.18e-07
   evolution_ABA ED dt=0.20 T=10  max|y-exact|=8.46e-04 scale 0.068  norm drift max|<psi|psi>(t)/<psi|psi>(0)-1|=8.97e-04
(d) python TDVP maxm=32 dt=0.20  max|y-exact|=7.32e-09  (DMRG GS energy err 7.1e-14)
(c) n=10 J=hz=3  E0=-14.311 Emax=18.206 bandwidth=32.52
   evolution_ABA ED dt=0.10 T=10  max|y-exact|=1.97e-04 scale 0.058
   evolution_ABA ED dt=0.20 T=10  max|y-exact|=1.11e-03 scale 0.057
```

```bash
cd <scratch>/reviews/realtime_4 && <scratch>/run3.sh 03_offset_public_and_fix.py 2>&1 | tee 03_offset_public_and_fix.out
```

`reviews/realtime_4/03_offset_public_and_fix.py`:

```python
"""Reviewer probe 3 for realtime_4.

(a) A constant added to H changes nothing physical in <psi(t)|O|psi(t)>,
    so evolve_and_measure(mode="ED") should not move. Measure it through
    the public route on the low-bandwidth chain at dt=0.1, the regime the
    tests use, with H and H + c*1 for c = 0, 10, 20.
(b) The suggested fix: expm_multiply(-1j*dt*H) step by step against the
    exact propagator on the same shifted operators.
"""
import numpy as np
from scipy.sparse import identity
from scipy.sparse.linalg import expm_multiply
from dmrgpy import spinchain, timedependent, multioperator

def chain(n, J, hz, c):
    sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        h = h + hz*(0.3 + 0.1*i)*sc.Sz[i] + 0.4*sc.Sx[i]
    if c != 0.0:
        h = h + c*multioperator.identity()
    sc.set_hamiltonian(h)
    return sc

n, dt, nt = 8, 0.1, 100
sc0 = chain(n, 1.0, 1.0, 0.0)
ed0 = sc0.get_ED_obj()
H0 = ed0.get_operator(sc0.hamiltonian).tocsc()
w, v = np.linalg.eigh(H0.toarray())
# a real product start: all up, then rotate site 0 to +x (a quench)
d = H0.shape[0]
rng = np.random.default_rng(3)
wf = rng.normal(size=d) + 0j
wf /= np.linalg.norm(wf)
O = ed0.get_operator(sc0.Sz[2]).toarray()
ts = dt*np.arange(nt)
c0 = v.conj().T@wf
ref = np.array([np.vdot(v@(np.exp(-1j*w*t)*c0), O@(v@(np.exp(-1j*w*t)*c0))) for t in ts]).real
print("exact <Sz_2>(t) range [%.4f, %.4f]" % (ref.min(), ref.max()))
for c in [0.0, 10.0, 20.0]:
    sc = chain(n, 1.0, 1.0, c)
    Hc = sc.get_ED_obj().get_operator(sc.hamiltonian).tocsc()
    print(" c=%5.1f  ||H_ED(c) - H_ED(0) - c*1|| = %.1e" %
          (c, abs(Hc - H0 - c*identity(d, format="csc")).max()))
    _t, y = timedependent.evolve_and_measure(sc, operator=sc.Sz[2], nt=nt, dt=dt,
                                             wf=wf.copy(), mode="ED")
    y = np.asarray(y).real
    print("   evolve_and_measure ED dt=%.1f T=10  max|y-exact|=%.2e" % (dt, np.max(np.abs(y-ref))))
    # suggested fix on the same operator
    x = wf.copy()
    ys = []
    for _ in range(nt):
        ys.append(np.vdot(x, O@x).real)
        x = expm_multiply(-1j*dt*Hc, x)
    print("   expm_multiply        dt=%.1f T=10  max|y-exact|=%.2e  norm-1=%+.1e"
          % (dt, np.max(np.abs(np.array(ys)-ref)), np.linalg.norm(x)-1))
```

`reviews/realtime_4/03_offset_public_and_fix.out`:


```
exact <Sz_2>(t) range [-0.0863, 0.0567]
 c=  0.0  ||H_ED(c) - H_ED(0) - c*1|| = 0.0e+00
   evolve_and_measure ED dt=0.1 T=10  max|y-exact|=2.12e-07
   expm_multiply        dt=0.1 T=10  max|y-exact|=7.22e-16  norm-1=+2.2e-16
 c= 10.0  ||H_ED(c) - H_ED(0) - c*1|| = 0.0e+00
   evolve_and_measure ED dt=0.1 T=10  max|y-exact|=1.71e-04
   expm_multiply        dt=0.1 T=10  max|y-exact|=7.91e-16  norm-1=+2.4e-15
 c= 20.0  ||H_ED(c) - H_ED(0) - c*1|| = 0.0e+00
   evolve_and_measure ED dt=0.1 T=10  max|y-exact|=5.64e-04
   expm_multiply        dt=0.1 T=10  max|y-exact|=1.06e-15  norm-1=+2.0e-15
```


**Suggested fix**: replace RK45 inside `scipy_evolution` by
`scipy.sparse.linalg.expm_multiply(1j*t*h, psi)`, keeping `evolve()`'s e^{+iht}
sign, so `evolution_DC` and the `-Hop` in `evolution_ABC` stay as they are. Measured
on the shifted operators: 7.2e-16, 7.9e-16 and 1.1e-15 off exact at c = 0, 10 and 20,
with norm error around 2e-15. Tightening `solve_ivp` to rtol=1e-10/atol=1e-12 only
moves the threshold, since the error still scales with dt times |E|, and costs more
steps than `expm_multiply` at large |E|dt; subtracting a reference energy in
`evolution_ABC` would remove most of the offset sensitivity but leave a non-unitary
integrator. When the dimension is a few thousand or less and `nt` is large, one
`eigh` per call and phase multiplication is exact and cheaper, but `expm_multiply`
is the general choice. The regression should pin offset invariance directly:
`evolve_and_measure(mode="ED")` with h and with h + 20*identity agreeing to about
1e-12. NUMBERS CHANGE on every `mode="ED"` real-time route (`evolve_and_measure`,
`evolution_ABA`, `submode="TD"`), towards exact: by 1e-8 to 1e-6 on the small chains
the tests use, 8e-4 at dt=0.2 and 1e-3 at dt=0.5 on the bandwidth-24 chain, and
5.6e-4 at dt=0.1 under a +20 offset.

### 18. `disentangle._is_hermitian` falls back to the bare proof only when a state's `MBO` is None, but an ED `State` carries its EDchain, which has no `is_hermitian`, so since `867e2b4` `disentangle_manifold` raises `AttributeError` on every ED manifold for every operator, including a proven `Sz0` it handled to 3.3e-15 before

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `misc`

**Where**: `src/dmrgpy/mpsalgebratk/disentangle.py:14-16` (`mbo =
getattr(wfs[0],"MBO",None); if mbo is not None: return mbo.is_hermitian(A)`);
`edtk/edchain.py:335-343` (a `State` stores its EDchain as `MBO`) and `:272`
(`get_excited_states` builds `State`s on the EDchain).

`867e2b4` moved `disentangle_manifold` onto the chain's Hermiticity probe, with a
fallback to the bare proof "when a state carries no chain". An ED `State` does carry
one, just not a `Many_Body_Chain`: its `MBO` is `pychain.build.Spin_chain` for
spins and `MBFermion` for fermions, neither of which has `is_hermitian`. So the call
raises for every operator, proven Hermitian, unprovable Hermitian and non-Hermitian
alike, whether the ED states come from `mode="ED"` as a keyword or from
`sc.mode="ED"` on the chain. It survived because all three pinned manifolds in
`tests/test_audit_2026_09_24b_misc.py` are `"python"` MPS, and the fallback test
`test_states_without_a_chain_fall_back_to_the_bare_proof` uses a synthetic class with
`MBO = None`, a shape no real state has; `disentangle_manifold` has no in-tree
caller in `src/` or `examples/`.

**Expected**: on an ED manifold (the 8 eigenstates of a generic 3-spin Hamiltonian
from `get_excited_states(n=8, mode="ED")`), an orthonormal basis diagonalizing A,
as on the `"python"` and v3 manifolds (Gram error 3e-15 to 1e-14, off-diagonal of A
4e-16 to 2e-15).

Repro, from the hunter (`03` the three backends, `04` the same ED manifold with the
pre-`867e2b4` decision patched back in):

```bash
cd <scratch>/misc && <scratch>/run3.sh 03_disentangle_backends.py 2>&1 | tee 03_disentangle_backends.out
```

`misc/03_disentangle_backends.py`:

```python
"""disentangle_manifold's _is_hermitian takes wfs[0].MBO.is_hermitian(A).
Which object is MBO on each backend's states? An ED State carries the
EDchain, not the Many_Body_Chain. Run the finding-1 operator
(1j*Sx0*Sy0, exactly -Sz0/2, unprovable) and a proven one (Sz0) on an ED
manifold, a v3 manifold and a "python" manifold, and report Gram error and
the off-diagonal of A on the output, or the exception."""
import warnings, traceback
import numpy as np
from dmrgpy import spinchain, mpsalgebra
from dmrgpy.mpsalgebratk.disentangle import get_representation

warnings.simplefilter("ignore")

def manifold(backend):
    iv = 3 if backend == "v3" else "python"
    sc = spinchain.Spin_Chain([2]*3, itensor_version=iv)
    h = sc.Sx[0]*sc.Sx[1] + 0.6*sc.Sy[1]*sc.Sy[2] + 0.3*sc.Sz[0]*sc.Sz[1] \
        + 0.7*sc.Sx[0] + 0.45*sc.Sz[1] + 0.2*sc.Sy[2] + 0.33*sc.Sz[2]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 8, 20
    np.random.seed(0)
    if backend == "ED":
        es, wfs = sc.get_excited_states(n=8, mode="ED")
    else:
        es, wfs = sc.get_excited_states(n=8)
    return sc, wfs

for backend in ["ED", "v3", "python"]:
    sc, wfs = manifold(backend)
    print("backend %-6s state type %s, MBO type %s" %
          (backend, type(wfs[0]).__name__, type(getattr(wfs[0], "MBO", None)).__name__))
    for label, A in [("Sz0 (proven)", sc.Sz[0]), ("1j*Sx0*Sy0 (unproven)", 1j*sc.Sx[0]*sc.Sy[0])]:
        try:
            np.random.seed(7)
            out = mpsalgebra.disentangle_manifold(wfs, A)
            G = np.array([[a.dot(b) for b in out] for a in out])
            ma = get_representation(out, A)
            print("   %-24s Gram err %.3e, off-diag of A %.3e" %
                  (label, np.max(np.abs(G - np.eye(len(out)))),
                   np.max(np.abs(ma - np.diag(np.diag(ma))))))
        except Exception as ex:
            print("   %-24s raised %s: %s" % (label, type(ex).__name__, ex))
```

`misc/03_disentangle_backends.out`:


```
backend ED     state type State, MBO type Spin_chain
   Sz0 (proven)             raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
   1j*Sx0*Sy0 (unproven)    raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
WARNING, state is not normalizable. Returning None
WARNING, state is not normalizable. Returning None
backend v3     state type MPS, MBO type Spin_Chain
   Sz0 (proven)             Gram err 3.998e-15, off-diag of A 1.264e-15
   1j*Sx0*Sy0 (unproven)    Gram err 3.331e-15, off-diag of A 4.480e-16
WARNING, state is not normalizable. Returning None
WARNING, state is not normalizable. Returning None
backend python state type MPS, MBO type Spin_Chain
   Sz0 (proven)             Gram err 6.439e-15, off-diag of A 1.961e-15
   1j*Sx0*Sy0 (unproven)    Gram err 1.155e-14, off-diag of A 1.182e-15
```

```bash
cd <scratch>/misc && <scratch>/run3.sh 04_disentangle_ed_prefix.py 2>&1 | tee 04_disentangle_ed_prefix.out
```

`misc/04_disentangle_ed_prefix.py`:

```python
"""Companion of 03: on the same ED manifold, the pre-867e2b4 decision
(the bare proof, A.is_hermitian()) rebuilt by patching _is_hermitian, to
show the ED route ran before 867e2b4 and is broken by it."""
import warnings
import numpy as np
from dmrgpy import spinchain, mpsalgebra
from dmrgpy.mpsalgebratk import disentangle
from dmrgpy.mpsalgebratk.disentangle import get_representation

warnings.simplefilter("ignore")
sc = spinchain.Spin_Chain([2]*3, itensor_version="python")
h = sc.Sx[0]*sc.Sx[1] + 0.6*sc.Sy[1]*sc.Sy[2] + 0.3*sc.Sz[0]*sc.Sz[1] \
    + 0.7*sc.Sx[0] + 0.45*sc.Sz[1] + 0.2*sc.Sy[2] + 0.33*sc.Sz[2]
sc.set_hamiltonian(h)
es, wfs = sc.get_excited_states(n=8, mode="ED")
print("ED state MBO is", type(wfs[0].MBO), "has is_hermitian:", hasattr(wfs[0].MBO, "is_hermitian"))

def report(tag, A):
    out = mpsalgebra.disentangle_manifold(wfs, A)
    G = np.array([[a.dot(b) for b in out] for a in out])
    ma = get_representation(out, A)
    print("  %-34s Gram err %.3e, off-diag of A %.3e" %
          (tag, np.max(np.abs(G - np.eye(len(out)))), np.max(np.abs(ma - np.diag(np.diag(ma))))))

orig = disentangle._is_hermitian
disentangle._is_hermitian = lambda w, A: A.is_hermitian()   # pre-867e2b4 decision
print("pre-867e2b4 decision (bare proof):")
report("Sz0 (proven)", sc.Sz[0])
report("1j*Sx0*Sy0 (unproven)", 1j*sc.Sx[0]*sc.Sy[0])
disentangle._is_hermitian = orig
print("867e2b4 decision:")
for tag, A in [("Sz0 (proven)", sc.Sz[0]), ("1j*Sx0*Sy0 (unproven)", 1j*sc.Sx[0]*sc.Sy[0])]:
    try:
        report(tag, A)
    except Exception as ex:
        print("  %-34s raised %s: %s" % (tag, type(ex).__name__, ex))
```

`misc/04_disentangle_ed_prefix.out`:


```
ED state MBO is <class 'dmrgpy.pychain.build.Spin_chain'> has is_hermitian: False
pre-867e2b4 decision (bare proof):
  Sz0 (proven)                       Gram err 3.331e-15, off-diag of A 5.689e-16
  1j*Sx0*Sy0 (unproven)              Gram err 3.402e-01, off-diag of A 8.505e-02
867e2b4 decision:
  Sz0 (proven)                       raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
  1j*Sx0*Sy0 (unproven)              raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
```


**Reviewer (CONFIRMED)**: both scripts reproduce: the ED manifold's states are
`State` with `MBO` of type `Spin_chain`, and `Sz0` and `1j*Sx0*Sy0` both raise
`AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'`, while v3
(Gram 4.885e-15 and 3.997e-15) and `"python"` (6.439e-15 and 1.155e-14) pass; with
the old decision patched in, `Sz0` gives Gram 3.331e-15 and `1j*Sx0*Sy0` 3.402e-01,
the second pass's finding 1 still live there. The same raise appears with
`sc.mode="ED"` set on the chain, on a fermionic ED manifold (`MBFermion`, `N0` and
`C0+Cdag0`), and for the genuinely non-Hermitian `Sx0+1j*Sy0`. A dense exact test,
M = `MO2matrix(A)` and ||M - M^H|| against ||M||, separates the cases cleanly
(exactly 0 for `Sz0`, `1j*Sx0*Sy0`, `N0` and `C0+Cdag0`; 2.828 against ||M|| = 2.0
for `Sx0+1j*Sy0`), and with it `1j*Sx0*Sy0` gives Gram 3.775e-15 and the fermion
operators 6.7e-16 and 8.9e-16. The behaviour is not documented, and an
`AttributeError` is not a `NotImplementedError` naming a restriction. Nothing
struck. The second pass's finding 1 Status, "falling back to the bare proof when a
state carries no chain", needs an addendum, since ED states never reach that
fallback.

```bash
cd <scratch>/reviews/misc_2 && <scratch>/run3.sh 01_disentangle_backends.py 2>&1 | tee 01_disentangle_backends.out
```

`reviews/misc_2/01_disentangle_backends.py` is `misc/03_disentangle_backends.py` (finding 18) unchanged, rerun; its output, `reviews/misc_2/01_disentangle_backends.out`:


```
backend ED     state type State, MBO type Spin_chain
   Sz0 (proven)             raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
   1j*Sx0*Sy0 (unproven)    raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
WARNING, state is not normalizable. Returning None
WARNING, state is not normalizable. Returning None
backend v3     state type MPS, MBO type Spin_Chain
   Sz0 (proven)             Gram err 4.885e-15, off-diag of A 9.254e-16
   1j*Sx0*Sy0 (unproven)    Gram err 3.997e-15, off-diag of A 3.244e-16
WARNING, state is not normalizable. Returning None
WARNING, state is not normalizable. Returning None
backend python state type MPS, MBO type Spin_Chain
   Sz0 (proven)             Gram err 6.439e-15, off-diag of A 1.961e-15
   1j*Sx0*Sy0 (unproven)    Gram err 1.155e-14, off-diag of A 1.182e-15
```

```bash
cd <scratch>/reviews/misc_2 && <scratch>/run3.sh 02_disentangle_ed_prefix.py 2>&1 | tee 02_disentangle_ed_prefix.out
```

`reviews/misc_2/02_disentangle_ed_prefix.py` is `misc/04_disentangle_ed_prefix.py` (finding 18) unchanged, rerun; its output, `reviews/misc_2/02_disentangle_ed_prefix.out`:


```
ED state MBO is <class 'dmrgpy.pychain.build.Spin_chain'> has is_hermitian: False
pre-867e2b4 decision (bare proof):
  Sz0 (proven)                       Gram err 3.331e-15, off-diag of A 5.689e-16
  1j*Sx0*Sy0 (unproven)              Gram err 3.402e-01, off-diag of A 8.505e-02
867e2b4 decision:
  Sz0 (proven)                       raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
  1j*Sx0*Sy0 (unproven)              raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
```

```bash
cd <scratch>/reviews/misc_2 && <scratch>/run3.sh 03_ed_routes_and_fix.py 2>&1 | tee 03_ed_routes_and_fix.out
```

`reviews/misc_2/03_ed_routes_and_fix.py`:

```python
"""Reviewer probe for misc_2. (a) Other ways of landing on an ED manifold:
a chain with sc.mode="ED" (no mode= kwarg), a fermionic chain's ED states,
and get_gs(mode="ED") single-state list. (b) Whether the 867e2b4 code
really raises there and what the pre-867e2b4 bare proof gave. (c) The
hunter's proposed exact dense test, monkeypatched in, on proven, unproven
Hermitian and genuinely non-Hermitian operators."""
import warnings
import numpy as np
from dmrgpy import spinchain, fermionchain, mpsalgebra
from dmrgpy.mpsalgebratk import disentangle
from dmrgpy.mpsalgebratk.disentangle import get_representation

warnings.simplefilter("ignore")

def report(tag, wfs, A):
    try:
        out = mpsalgebra.disentangle_manifold(wfs, A)
        G = np.array([[a.dot(b) for b in out] for a in out])
        ma = get_representation(out, A)
        print("  %-30s Gram err %.3e, off-diag %.3e" % (tag,
              np.max(np.abs(G - np.eye(len(out)))),
              np.max(np.abs(ma - np.diag(np.diag(ma))))))
    except Exception as ex:
        print("  %-30s raised %s: %s" % (tag, type(ex).__name__, ex))

# spin chain with mode="ED" on the chain
sc = spinchain.Spin_Chain([2]*3, itensor_version="python")
h = sc.Sx[0]*sc.Sx[1] + 0.6*sc.Sy[1]*sc.Sy[2] + 0.3*sc.Sz[0]*sc.Sz[1] \
    + 0.7*sc.Sx[0] + 0.45*sc.Sz[1] + 0.2*sc.Sy[2] + 0.33*sc.Sz[2]
sc.set_hamiltonian(h)
sc.mode = "ED"
es, wfs = sc.get_excited_states(n=8)
print("spin chain, sc.mode='ED': state", type(wfs[0]).__name__, "MBO", type(wfs[0].MBO).__name__,
      "has is_hermitian", hasattr(wfs[0].MBO, "is_hermitian"))
spinops = [("Sz0 (proven)", sc.Sz[0]), ("1j*Sx0*Sy0 (unproven)", 1j*sc.Sx[0]*sc.Sy[0]),
           ("Sx0+1j*Sy0 (non-Herm)", sc.Sx[0] + 1j*sc.Sy[0])]
for t, A in spinops: report(t, wfs, A)

# fermionic chain, ED
fc = fermionchain.Fermionic_Chain(3, itensor_version="python")
hf = sum(fc.Cdag[i]*fc.C[i+1] for i in range(2))
hf = hf + hf.get_dagger() + 0.3*fc.N[0] + 0.7*fc.N[0]*fc.N[1]
fc.set_hamiltonian(hf)
es2, wfs2 = fc.get_excited_states(n=8, mode="ED")
print("fermion chain ED: state", type(wfs2[0]).__name__, "MBO", type(wfs2[0].MBO).__name__,
      "has is_hermitian", hasattr(wfs2[0].MBO, "is_hermitian"))
fops = [("N0 (proven)", fc.N[0]), ("C0+Cdag0 (proven?)", fc.C[0]+fc.Cdag[0])]
for t, A in fops: report(t, wfs2, A)

print("pre-867e2b4 bare proof on the spin ED manifold:")
orig = disentangle._is_hermitian
disentangle._is_hermitian = lambda w, A: A.is_hermitian()
for t, A in spinops: report(t, wfs, A)

def dense_test(w, A):
    mbo = getattr(w[0], "MBO", None)
    if mbo is not None and hasattr(mbo, "is_hermitian"): return mbo.is_hermitian(A)
    if mbo is not None and hasattr(mbo, "MO2matrix"):
        M = mbo.MO2matrix(A)
        M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
        d = np.linalg.norm(M - M.conj().T); n = np.linalg.norm(M)
        print("    dense test: ||M-M^H|| = %.3e, ||M|| = %.3e" % (d, n))
        return d <= 1e-12*max(n, 1.0)
    return A.is_hermitian()
print("hunter's dense-matrix fix, monkeypatched, spin ED manifold:")
disentangle._is_hermitian = dense_test
for t, A in spinops: report(t, wfs, A)
print("hunter's dense-matrix fix, fermion ED manifold:")
for t, A in fops: report(t, wfs2, A)
disentangle._is_hermitian = orig
```

`reviews/misc_2/03_ed_routes_and_fix.out`:


```
spin chain, sc.mode='ED': state State MBO Spin_chain has is_hermitian False
  Sz0 (proven)                   raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
  1j*Sx0*Sy0 (unproven)          raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
  Sx0+1j*Sy0 (non-Herm)          raised AttributeError: 'Spin_chain' object has no attribute 'is_hermitian'
fermion chain ED: state State MBO MBFermion has is_hermitian False
  N0 (proven)                    raised AttributeError: 'MBFermion' object has no attribute 'is_hermitian'
  C0+Cdag0 (proven?)             raised AttributeError: 'MBFermion' object has no attribute 'is_hermitian'
pre-867e2b4 bare proof on the spin ED manifold:
  Sz0 (proven)                   Gram err 3.331e-15, off-diag 5.689e-16
  1j*Sx0*Sy0 (unproven)          Gram err 3.402e-01, off-diag 8.505e-02
  Sx0+1j*Sy0 (non-Herm)          Gram err 1.000e+00, off-diag 1.056e-08
hunter's dense-matrix fix, monkeypatched, spin ED manifold:
    dense test: ||M-M^H|| = 0.000e+00, ||M|| = 1.414e+00
  Sz0 (proven)                   Gram err 3.331e-15, off-diag 5.689e-16
    dense test: ||M-M^H|| = 0.000e+00, ||M|| = 7.071e-01
  1j*Sx0*Sy0 (unproven)          Gram err 3.775e-15, off-diag 2.762e-16
    dense test: ||M-M^H|| = 2.828e+00, ||M|| = 2.000e+00
  Sx0+1j*Sy0 (non-Herm)          Gram err 1.000e+00, off-diag 1.056e-08
hunter's dense-matrix fix, fermion ED manifold:
    dense test: ||M-M^H|| = 0.000e+00, ||M|| = 2.000e+00
  N0 (proven)                    Gram err 6.661e-16, off-diag 3.520e-16
    dense test: ||M-M^H|| = 0.000e+00, ||M|| = 2.828e+00
  C0+Cdag0 (proven?)             Gram err 8.882e-16, off-diag 4.996e-16
```


**Suggested fix**: not `hasattr(mbo, "is_hermitian")` with a fall back to the bare
proof, which would bring back the second pass's finding 1 on ED (Gram 0.3402 for
`1j*Sx0*Sy0`). Add an `is_hermitian(self, A)` method to `edtk/edchain.py::EDchain`,
next to `is_zero_operator`, `applyoperator` and `overlap`, since the EDchain already
mirrors the `Many_Body_Chain` method surface; that keeps `_is_hermitian` and its
docstring true and serves any later caller that asks the `MBO`. It should build the
dense matrix `M = self.MO2matrix(A)` (`.toarray()` when sparse) and test ||M - M^H||
at roundoff relative to ||M||, which needs nothing from `get_dagger()` and so stays
exact for names off the parity table, such as the parafermion `Sig`. It should not
reuse `algebra.is_hermitian` as it stands, whose `is_zero_matrix` uses the absolute
1e-8 of finding 12. The regression should build its manifold from
`get_excited_states(mode="ED")` on a spin and a fermion chain, with a proven, an
unprovable Hermitian and a non-Hermitian operator. Lands with finding 12. No
returned number changes: a raise becomes a result, the same as before `867e2b4` for
a proven operator and orthonormal instead of 0.340 off for an unprovable Hermitian
one.

## Ruled out

Executed and found clean, so the next pass need not redo them.

- **groundstate**: `get_kondo_spectrum(mode="DMRG", n_gs=2)` leaves the chain and
  the session consistent: on the crossing chain (eps=1e-5) |<wf0|session>|^2 =
  1.000000000000 afterwards, `gs_energy()` is unchanged to every digit, and a KPM
  and a CVM correlator after the call equal the ones before it bit for bit, on
  `"python"` and v3. The injected state reaches the session unswept: after
  `set_gs` of a non-degenerate excited eigenstate, CVM, CVM_explicit, ROOTN, TD and
  EX on `"python"`, v3 and v2 measure that state from its own energy (established
  by finding 1's reviewer on (S+_0,S-_0), since the hunter's Sz0-Sz0 cannot tell the
  two states apart on that chain), and `gs_energy()` stays <x|H|x> = -0.850000
  through five correlator calls. `gs_energy()` after `set_initial_wf(x)` on a
  never-solved 24-site chain returns exactly <x|H|x> on `"python"` and v3. `dcex`'s
  cached basis across `set_hamiltonian(restart=False)` raises from its 1e-6 span
  check when the new ground state lies outside the old basis (with a commuting H2 it
  is silent instead; see finding 5).
- **kpm**: the 1.5\*\|\|vi\|\|\*\|\|vj\|\| guard never fires on a correct run: on
  `"python"`, 8 sites, `kpmmaxm` 3, 4, 8 and 16, `delta` 0.1 and 0.03 (133 and 442
  moments), auto (Sz0,Sz0) with acceleration on and off, cross (Sz0,Sz3) and complex
  (Sp0,Sm0), the largest |mu|/bound recorded inside the loop is 0.9859; on v3 and
  v2, 14 sites, `kpmmaxm` 8, 16 and 32, `delta` 0.1 and 0.02 (240 and 1201
  moments), 0.9313, equal to the exact-moment ratio to four digits at
  `kpmmaxm=32`. The bound is computed from the vectors the recursion uses on every
  route (after the `kpmmaxm` truncation on DMRG; B\|0> and A^dag\|0> on ED; `wf` in
  `distribution_kpm`). The energy-truncated loop on correct runs stays at 0.8446 to
  0.9451 on `"python"` and returns on v3 (a cross pair at `kpm_scale=0.6` raises the
  truncation guard, not the moment guard). A charge-conserving pair (Sz0,Sz2) at
  Sz=0 in a conserved sector gives 0.7288 on `"python"` and v3, spectra within
  1.9e-7 and 7.8e-13 of exact, and a charge-changing (C0,Cdag0) at Nf=2 is refused
  with `ValueError` before any moment. `get_distribution` on both modes and all
  three DMRG backends gives line weights within 1.3e-4 of the exact projector
  weights (the A=/B= form within 7e-6); `scale=0.9` raises on every route and
  `scale=1.2` returns 0.9993 (DMRG) and 0.9996 (ED). Zero and near-zero start
  vectors on an exactly polarized chain neither raise nor give NaN. `_same_mps`'s
  absolute 1e-10 is unreachable through operator coefficients (at eps=1e-6 the
  cross pair is exact to 6e-15; below about 1e-7 the operator is gone before any
  recursion, finding 13 and `clean_threshold`). The `i=None`/`j=None` defaults are
  bit for bit right on DMRG and ED, and no in-tree wrapper overrides
  `get_dynamical_correlator`. NH-DMRG's ground energy coming back as the conjugate
  of ED's is a degenerate conjugate pair at the bottom of the spectrum, not a
  defect.
- **realtime**: `867e2b4`'s direction and conjugation fixes hold on every
  integrator for `evolution_ABA` with a non-unit-norm start, a complex Hamiltonian
  and a non-Hermitian fermionic observable (4 sites, hopping -e^{0.4i},
  \|\|Cdag_0\|GS>\|\|^2 = 0.592, B = Cdag_1 C_2, against a hand-written e^{-iHt}):
  ED 2.4e-10, `"python"` TDVP and TDVP_GSE 1e-13, TEBD and AUTO 7.7e-05, MPO
  1.3e-04, v3 the same, v2 MPO 1.3e-04, every row 0.226 from the backward
  trajectory and 0.301 from its conjugate. The TD series `evolution_DC`, whose
  conjugation `867e2b4` kept, holds on every integrator but v3 `TDVP_GSE`
  (finding 14) for a complex-weight pair (Cdag_0,C_2), max\|Im M_n\| = 0.142. The
  whole infinite-chain TD pipeline on a chiral ferromagnet (a DM term in a field,
  single-magnon closed form eps(q) = h + J(1-cos q) + Dm sin q) puts the (Sx,Sx)
  peaks at eps(k), not eps(-k), on `"python"` and v3 (k=+pi/2: 2.4962 and 2.4949
  against 2.5000; k=-pi/2: 1.5040 and 1.5052 against 1.5000), and the tangent-space
  excitation ansatz agrees on the momentum label. The connected-background
  subtraction with a complex static value and A != B^dagger ((S+,Sz), <S+> =
  0.322i) gives 0 to 5e-6 (`"python"`) and 1e-16 (v3) at x != 0, t=0, where a
  conjugated background would read 0.246i. The fermionic IBC window on a complex
  Hamiltonian runs forward with the unconjugated S(x,t) that `sxt_to_skomega`
  assumes. By grep and reading, nothing in `kondospectrumtk/` calls
  `evolve_and_measure`, `evolution_ABA`, `evolution_ABC` or ED `evolve()`
  (`dmrgtwotime.py` drives `_session.tdvp_step` directly), so the two sign flips
  cannot reach the Kondo third-order term.
- **misc**: the strip's other caller, `quench_tdvp_gse`, gives the same
  <Sz0(t)Sz0> from a padded-then-unpadded ground state as from a never-padded one,
  to 1.5e-15 and 1.9e-15 at cutoff 1e-10 (sweeps 0 and 3) and 1.5e-15 and 1.2e-15
  at cutoff 0. By reading, the lossless sweep's cutoff=0 drops a singular value only
  when the tail weight is exactly 0, so the copy-and-compare cannot discard a small
  nonzero one. Every chain setting the KPM route reads exists on a fresh chain, all
  of them passed together at their defaults to `kpm_finite`'s `window_chain_kwargs`
  give the default spectrum bit for bit, and the lazily defaulted `getattr`
  settings in `src/` are all set in `__init__`. The Kondo grid sort holds on the
  DMRG KPM route too (a coarse grid plus a duplicated fine block, identical to the
  sorted grid to 0.000e+00 on peaks 4.714 and 14.271, permuted, as a list and
  descending). `disentangle_manifold` on v3 and `"python"` manifolds gives Gram
  errors 3.3e-15 to 1.2e-14 for `1j*Sx0*Sy0` and `Sz0`, and no consumer relied on
  the old `eig` basis.

## New leads, not reviewed

Observed along the way and not handed to a reviewer of their own, so they are leads
rather than findings.

- The `maxde` loop compares a per-site fluctuation (`de/self.ns`,
  `groundstate.py:380-381`, `4731b5a`) against `maxde`, while
  `gs_energy_fluctuation()` returns the total, and neither the user guide nor the
  docstring says `maxde` is per site; `examples/groundstate/
  GS_enforce_maximum_fluctuation` plots the total against the per-site request on
  one axis (at `maxm=10` on v3 the loop reads 1.896e-03 per site where
  `gs_energy_fluctuation()` prints 1.896e-02 for the same state). Turned up by
  finding 7's reviewer, with executed evidence, but given no reviewer of its own.
- After `set_gs(x)` on a non-Hermitian chain the NH-KPM route pairs x with the
  stale `nh_left_wf` of the last NH-DMRG solve (or raises `AttributeError` if there
  was none), and `gs_energy(wf0=x)` on a non-Hermitian chain reaches
  `gs_energy_nhdmrg`, which drops `wf0` as an unknown keyword (by reading).
- `get_gs(wf0=x)` on a current chain returns the stored state without reading x
  (`manybodychain.py:1212`), the shape `867e2b4` fixed for `gs_energy(wf0=x)` only
  (by reading).
- Any solver-parameter change after `set_gs(x)` (`cutoff`, `noise`, `maxm`) makes
  `gs_is_current` False and the next read re-solves from x, discarding the
  caller's state; arguably the contract, but an injected state does not survive a
  harmless-looking tweak (by reading).
- `Thermal_Spin_Chain.get_gs` assigns `MBChain.wf0` and `MBChain.hamiltonian`
  directly (`thermal.py:59-62`), so a correlator on `MBChain` finds
  `hamiltonian_on_session` False and re-solves over the annealed state; the same
  before `867e2b4` (by reading).
- `gs_energy_generalized` re-sends the Hamiltonian with `session.set_hamiltonian`
  without updating `_session_ham_cache` (`groundstate.py:578`), so on a chain whose
  cache is empty the next correlator may re-solve a plain ground state over the
  generalized one, against its own CAVEAT (by reading).
- The v3 -3.000000 after `maxde` in finding 6, a re-solve at `maxm=3` warm-started
  from a truncated `maxm=12` state landing exactly on a dimer product energy, was
  not investigated.
- Every DMRG backend returns a wrong ground state for a Hamiltonian written in
  small energy units, on v2 and v3 by a mechanism other than finding 13: 6-site
  Heisenberg s*H gives `gs_energy()`/s of -1.351 and -1.320 on v3 at s=5e-7 (exact
  -2.493577), -1.25 on v3 and -1.747 on v2 at 3e-7, -1.25 on both at 1e-7, with
  the returned state's <wf|H|wf> not matching the returned energy on v3 at 5e-7;
  -1.25 is the Neel ZZ energy, which suggests the XX+YY terms vanish in the vendored
  AutoMPO's compression (`isZero` 1E-13, `toMPO` Cutoff 1E-13, by reading) or an
  absolute Davidson `ErrGoal`; not located.
- `multioperator.clean_threshold = 1e-8` drops every term with |coef| <= 1e-8 on
  every backend including ED: `(1e-8*Sz0).to_terms() == []`, `vev` and correlator
  exactly 0, and `gs_energy()`/s = 0 at s=1e-8 on all three DMRG backends; absolute,
  origin `593b394` or earlier.
- The second pass's finding 3 Status sentence that a harsh `kpmmaxm` truncation
  raises is a 4-site accident: on 8 sites `kpmmaxm=3` returns spectra 84 to 181 per
  cent wrong at moment ratios 0.85 to 0.99 with no raise. Not a defect of the
  guard, which claims only to catch growth, but the Status line and the guard's
  comment overstate what it catches.
- The moment guard's message blames the band edge or `kpm_scale` when the cause is
  a non-Hermitian Hamiltonian misdispatched to the Hermitian KPM (finding 12).
- ED's `algebra.biorthogonal_ground_state` raised `ValueError` (<vl|vr> about
  1e-15) on the conjugate-pair degenerate ground state at s=1e-3 while it passed at
  s=1: left and right eigenvectors picked from different members of the pair.
- In the energy-truncated loops of `"python"` and v3 the start vector and its first
  product are not truncated, so mu0 and mu1 carry the out-of-window weight every
  later moment has lost (by reading; `kpm_energy_truncate` territory).
- `fermionchaintk/dynamicalcorrelator.py` has no importer and carries the typo
  `name=(mi,mj)**kwargs` in `getd`; `pychainwrapper.py:76` and
  `infinitechain.py:1404` call `dynamics.`/`kpmdmrg.get_dynamical_correlator`
  directly and bypass the public method's `i=`/`j=` check (by reading).
- `julia_live` `TDVP_GSE` from an edge-site ladder start was not run (out of scope);
  ITensorMPS's own `expand(alg="global_krylov")` may or may not share finding 14.
- `_damped_sum_at` zeroes every `es` outside the FFT band [-pi/dt, pi/dt), which at
  dt=0.5 empties the default window above 6.28; the docstring states it.
- `sxt_to_skomega` on a pair with A != B^dagger returns the one-sided transform of
  conj(S(k,t)), whose real part carries a dispersive term; documented, and not
  measured against a two-sided density.
- `disentangle_manifold` now runs the chain's random-witness probe, which draws
  from the global `np.random` state and overwrites the chain's one-entry
  `_is_hermitian_cache`, so a seeded calculation sees different random numbers
  afterwards and the next `gs_energy()` re-probes the Hamiltonian (by reading).
- The probe's witness bond dimension is not the temporary `self.maxm = min(maxm, 8)`
  it sets: `randommps.random_mps_dummy` calls `self._session.random_mps()`, which on
  `"python"` uses the session's stored `maxm` (by reading).

## Statements this hunt overturns

Each of these is written as established and is contradicted by an entry above;
when the corresponding fix lands, the older text should gain a pointer here.

- The second pass's finding 12 Status, "a stored state whose Hamiltonian is no
  longer the session's is solved again" (true for the correlator path only);
  `groundstate.solver_key`'s docstring premise that changing the Hamiltonian "goes
  through `set_hamiltonian()`, which resets `computed_gs` outright"; `dcex.py:14-15`
  and its error text at `:182`, that the cache is "cleared by `set_hamiltonian()`"
  (finding 5).
- The second pass's finding 11 Status, "`mode="ED"` honours the same `set_gs` to
  3.3e-6" (false for `submode="ED"`, finding 2); `set_gs`'s docstring promise of
  "its own energy <wf|H|wf>" to every consumer (not KPM's axis, finding 1);
  `user_guide.md:4245`, "A state set by hand now reaches the solver" (not SECTOR,
  finding 8, nor ED `submode="ED"`, finding 2).
- `_take_injected_state`'s docstring, that the pre-fill "costs nothing when the
  edges are cached already", and the second pass's finding 11 Status, "one reduced
  solve when the edges are missing" (findings 9 and 10).
- `kpmdmrg.py:124`, "validated here, once, ahead of both session calls and of the
  ground state"; `validate_kpm_n_scale`'s docstring, "before any ground-state
  work"; the second pass's ruled-out line "rejected before the ground state on every
  Hermitian route" (finding 11).
- The comment at `mpsalgebra.py:355-357`, that `applyoperator()` normalizes its
  result (finding 12).
- The second pass's finding 1 Status, "falling back to the bare proof when a state
  carries no chain" (an ED state never reaches it, finding 18; and the probe it now
  trusts is absolute, finding 12).
- `promote_to_dense`'s docstring, "keeping the state computed in it" (findings 3
  and 4); the 2026-09 record's finding 10 statement that v3 "cannot leave that
  sector" (finding 4).
- The "isolated" `TDVP_GSE` failure of
  `examples/time_evolution/tdvp_gse_VS_ED_time_evolution/main.py:28-41`,
  `docs/documentation.md` around line 4449 and `ROADMAP.md` line 82, and the comment
  in `mpscpp3/chain_session.h` that `Cutoff=0` "only trims once the count exceeds
  maxm_", which ITensor v3's `truncate()` contradicts for exact zeros (finding 14).
- `user_guide.md:705-714`'s fluctuation paragraph, silent that below full bond
  dimension the number is set by `maxm`, and the plot of
  `examples/groundstate/GS_enforce_maximum_fluctuation` (findings 6 and 7).

## Shared helpers

The three-slot wrapper every script ran through, with the same substitutions as the
splices:

`run3.sh`:

```bash
#!/bin/bash
# Runs one python script under a workflow-wide cap of three concurrent runs.
# Usage, from the folder holding the script:
#   <this file> NN_slug.py [args] 2>&1 | tee NN_slug.out
# It waits (checking every 5 s) until one of three lock slots is free, then runs
# the script with threads pinned and this checkout's src on PYTHONPATH. The
# slot is released when the script exits, including when it is killed.
LOCKDIR=<scratch>/locks
while true; do
  for i in 1 2 3; do
    exec 9>"$LOCKDIR/slot$i.lock"
    if flock -n 9; then
      echo "[run3] slot $i acquired" >&2
      MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
        MPLBACKEND=Agg PYTHONPATH=<repo>/src \
        python3 "$@"
      rc=$?
      exec 9>&-
      exit $rc
    fi
    exec 9>&-
  done
  sleep 5
done
```

The `echo` line is the one change from the second pass's copy: it tells a Bash
timeout spent waiting for a slot from one spent running.
