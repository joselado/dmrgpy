# Audit, 2026-09-24 (second pass): five-lens hole hunt of the fix pass `30200a4`

Findings from a five-lens automated audit of the Python layer, run 2026-09-24 on
`6398577` (clean tree, both compiled extensions current, `make -n pybind`
reporting nothing to do in `mpscpp2` and `mpscpp3`), scoped to one commit,
`30200a4`, the fix pass that closed the sixteen findings of
`audit_2026_09_24_hole_hunt.md` (`6398577` on top of it is documentation only,
and `30200a4` itself touched no C++ and no Julia file). The reason for scoping a
hunt to a single fix commit is the previous record's own measurement: seven or
eight of its sixteen holes came in with `765b537`, the one commit that closed
the hunt before it. Each lens took one of `30200a4`'s five fix clusters. Every
finding below carries a repro that was actually executed and its verbatim
output, and every one was then handed to an independent reviewer whose brief
was to *refute* it; four of them were turned up by a reviewer while reviewing a
neighbouring candidate and got a reviewer of their own before entering the
record. None was refuted outright, and almost every one came back narrowed, so
each entry keeps its struck sub-claims visible next to what survived.

This file is the evidence, not a task list, the same convention as the three
earlier records: it records what was observed and how to reproduce it, so that
a fix, or a decision that the behaviour is intended after all, does not have to
re-derive any of it. Fixed entries gain a `**Status**` line under their
classification line and are kept rather than deleted, since the repro doubles
as the regression check.

All 18 have now been addressed: every entry carries a `**Status**` line
reading FIXED, the regressions live in four files,
`tests/test_audit_2026_09_24b_<cluster>.py` for `groundstate`, `kpm`,
`realtime` and `misc`, and finding 3 needed a rebuild of both compiled
extensions. The four files pass together (144 passed). Two things are left
open and say so in their Status lines: the sliver just below
`kpm_scale=1/2` that no moment bound sees (findings 2 and 3) and
`julia_live`'s warm start (finding 12). The fixes that changed numbers
rather than behaviour are listed together in `CLAUDE.md`'s paragraph for
this audit.

What the hunt says about `30200a4` itself, in one paragraph: the fixes hold on
the ground each one was written for (the nulls are collected in "Ruled out"),
and the new holes sit where a fix met a neighbour it did not look at. The
previous record's finding 12 fix, the `n_gs=` average, installs states through a
round trip in `dynamics.py` that re-solves them; the `kpm_n_scale` check was
placed one frame ahead of the branch that decides whether the value is read; the
direct Fourier sum was wired in everywhere except the one reduction it was most
needed in; the padding strip reads the flag rather than the state. Twelve of the
eighteen findings are older than `30200a4` and were reached by probing next to
it, and each entry says which. The oldest two, ED real-time evolution running
backwards (finding 9) and DMRG `evolve_and_measure` returning the conjugate
(finding 10), date from 2020: neither can show on the real Hermitian observables
every existing test measures, and on a non-Hermitian observable with real matrix
elements, S+ or a hopping, the two cancel, so ED and DMRG agree while both are
wrong.

## The five lenses

| Lens | Brief |
|---|---|
| `operators` | Can the `_LOCAL_LADDER` rule in `canonical.py` and `get_dagger()`'s `_dagger_phase` (`ISy` to `-ISy`) prove something false, refuse something they used to prove, or change an operator's value through any consumer of the dagger or of the canonical form? |
| `kpm` | Does every KPM route reconstruct from exactly the calibrated n moments, normalize the same way and validate `kpm_n_scale` before dispatch, and did anything that consumed the old n+2 count, the old ED distribution normalization or the old rounding quietly break? |
| `realtime` | Is the direct per-frequency Fourier sum right on every branch that reaches it, and are ED TD on the cached state, `evolution_ABA`/`evolve_and_measure` on ED without `wf=`, TDZ's named parameters and `lehmann_density_from_one_sided`'s `i=`/`j=` right on every pair shape and in every state the chain can be in? |
| `kondo` | Are `mode="ED"`'s keyword check, `mode="DMRG"`'s `n_gs=` average, the potential term's trapezoid weights and the `es` coverage requirement right, and do they leave the chain and the two modes consistent with each other? |
| `pyitensor` | Is the deflated Lanczos per value of the excitation ansatz right on every spectrum shape, and does the `set_pad_bonds` exemption of the one-site TDVP route reach exactly the routes it should? |

## Scope

Out of scope by construction: the sixteen findings of
`audit_2026_09_24_hole_hunt.md` as originally stated (a new defect in how one of
them was fixed is exactly in scope), together with that record's "New leads" and
left-open items; every entry of `audit_2026_08_hole_hunt.md` and
`audit_2026_09_hole_hunt.md`; the open `docs/known_issue_*.md` items, in
particular the `kpm_energy_truncate` window problem; the legacy bugs `CLAUDE.md`
says are deliberately reproduced (`evoloperator`'s z^3/6 term on `H2`, the
`"moise"` key, the unreachable `"tevol_fit_td"` branch); vendored ITensor; gaps
`ROADMAP.md` marks as absent; and the two documented KPM residuals (the
`sqrt(1-x^2)` width profile and the Jackson line shape). A finding that lies
outside `30200a4` but was reached by probing next to it is reported with its
origin named, since a hole the suite passes through is worth the record whatever
commit brought it in; the at-a-glance table carries that origin as a column.

`itensor_version="julia_live"` was in scope for exactly one planned probe in the
`kpm` lens, whether `mpsjulialive/dynamics.py::_kpm_dynamical_correlator`
reconstructs from the same n as v3 (it does, and its curve is within 8.4e-12 of
v3's), and for one reviewer run on finding 10 (it carries the same conjugation);
everything else excludes it, and its KPM guard (finding 3) and time-evolution
return (finding 10) are covered by reading.

As a baseline, the tests nearest the commit were run on the frame before any
lens reported:

```bash
pytest tests -q -p no:cacheprovider -k "not julia_live and (audit_2026_09 or kpm or kondo or canonical or multioperator or excitation or pad_bonds or tdz or time or dynamical or distribution or correlator or dagger or hermit or idmrg or vumps or infinite)"
```

gave 957 passed, 2 skipped, 478 deselected in 39 minutes 30 seconds, so every
finding below is a hole the suite passes through rather than code already known
to be red.

Every repro was run as

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 <script> [args]
```

from the folder holding the script. Midway through the review phase the run was
capped at three concurrent Python processes across all agents together, and from
then on every script went through a small wrapper, `run3.sh` (in "Shared
helpers"), which waits for one of three `flock` slots and sets exactly the
environment above; scripts started before the cap ran directly, at most two per
agent. No finding makes a timing claim, so the change affects no number. The
scripts and their outputs are spliced into this file from disk rather than
retyped, because the scratch folders they ran in do not survive the session;
`<repo>` stands for the checkout and `<scratch>` for the scratch folder, the only
edit made to either. The outputs drop four kinds of line and nothing else: the
repeated notice "ITensor v3's two-site DMRG can't handle a chain this short ...,
using default ED routines" that `mode.py` prints on every call below three sites,
the per-frequency "CVM in E = ... CG iterations" progress lines, the
byte-compile `SyntaxWarning`s of files an AST scan parsed, and the
missing-extension `UserWarning` an archived older tree prints. Where one script
was run on several backends, one output is spliced and the others are named.
Reviewers returned their reports as text, since the harness refuses report files
from subagents, so the reviewer paragraphs below are condensed from those
reports while the scripts and outputs are theirs verbatim.

## The findings at a glance

Origin says where the defect came from: `30200a4` for what the fix pass itself
introduced, a commit or a period for what is older and was reached by probing next
to it. The fix clusters are named here; unlike the previous hunt's they are not
fully file-disjoint, since `manybodychain.py` is touched by `ground-state` (the
state setters) and by `kpm` (`get_dynamical_correlator`), and `pyitensor/chain.py`
by `kpm` (the moment guard) and by `misc` (the padding strip), so those pairs are
taken one after the other. Finding 3 is the only one that needs a rebuild of the
compiled extensions.

| # | Finding | Severity | Origin | Numbers change on fix | Fix cluster |
|---|---|---|---|---|---|
| 1 | `disentangle_manifold` takes `eig` on the one-sided proof, a non-orthonormal basis for an unproven Hermitian operator | LOW | `fe3580e` or earlier | only for an unproven Hermitian operator | misc |
| 2 | ED KPM has no divergence check since it adopted the DMRG window | MEDIUM | `765b537` | no (validation; a wrong spectrum becomes a raise) | kpm |
| 3 | the DMRG moment guard is 1e3 to 5e3 times looser than the exact bound and not scale-invariant | MEDIUM | `a67228e`, `695c452` | no (validation); needs a rebuild | kpm |
| 4 | ED KPM ignores `kpm_energy_truncate` | LOW | `51357d8`, `765b537` | no (validation) | kpm |
| 5 | the ED `kpm_n_scale` check sits ahead of the non-Hermitian branch | LOW | `30200a4` | no | kpm |
| 6 | `kpm_finite`'s `window_chain_kwargs` accepts and ignores a misspelled key | LOW | `d22a576` | no (validation) | misc |
| 7 | `sxt_to_skomega` is still on the FFT-plus-interpolation stage, 20 per cent at a converged window | MEDIUM | `30200a4` | yes, towards exact | realtime |
| 8 | every infinite-chain S(k,omega) is mirrored in frequency | MEDIUM | the `idmrg_window` port | yes, every infinite-chain S(k,omega) | realtime |
| 9 | ED `evolve_and_measure`/`evolution_ABA` run backwards in time | MEDIUM | `f54b3c2` (2020) | yes, on complex H, start or observable | realtime |
| 10 | DMRG `evolve_and_measure` returns the conjugate of <O> | LOW to MEDIUM | `9e33f32` (2020) | yes, by twice Im<O> | realtime |
| 11 | `set_gs()` is ignored by every DMRG correlator and reverted by it | LOW to MEDIUM | `a8233ab` | only after `set_gs` | ground-state |
| 12 | `set_initial_wf`/`set_initial_wf_guess` never reach the session; the TFIM example plots the wrong branch | MEDIUM | `a8233ab` | yes, for every caller of the two setters | ground-state |
| 13 | `n_gs>1` re-sweeps each member, the excited one relaxes | MEDIUM | `30200a4` | only on v2/v3 at a split below delta | ground-state |
| 14 | after `n_gs>1`, `gs_energy()` is the last member's energy | LOW to MEDIUM | `30200a4` | only between an `n_gs` call and the next correlator | ground-state |
| 15 | `n_gs` is inert under `submode="EX"` | LOW to MEDIUM | `30200a4` | no (validation) now; yes under the reprojection | ground-state |
| 16 | the public correlator drops `i=`/`j=` next to an operator pair | LOW | `1b87543`, `30200a4` | no (validation) | kpm |
| 17 | both Kondo terms assume an increasing `es` | LOW | pre-`30200a4` | only on non-increasing grids | misc |
| 18 | the padding strip reads the flag, not the state | LOW | `30200a4` | only padded-then-unpadded | misc |

The clusters and their files: `ground-state` (11 to 15: `dynamics.py`,
`groundstate.py`, the setters in `manybodychain.py`, `spinchain.py`, `dcex.py`);
`kpm` (2 to 5 and 16: `algebra/kpm.py`, `edtk/dynamics.py`, `manybodychain.py`'s
`get_dynamical_correlator`, `pyitensor/chain.py`'s moment guard, both
`chain_session.h`, `mpsjulialive/kpm.jl`); `realtime` (7 to 10:
`timedependent.py`, `edtk/timedependent.py`, `mpsjulialive/timedependent.py`,
`tests/test_infinite_chain.py`, `examples/idmrg/td_dynamical_correlator`); `misc`
(1, 6, 17, 18: `mpsalgebratk/disentangle.py`, `infinitechain.py`'s `kpm_finite`,
`kondospectrumtk/`, `pyitensor/chain.py`'s padding strip). Two pairs must land
together: 9 with 10 (either alone makes ED and DMRG disagree on S+ where today they
agree by cancellation), and 7 with 8 (same function, same test); 11 to 14 share one
fix at `dynamics.py:190` and `gs_energy_single`.

## Findings

### 1. `disentangle_manifold` chooses between `eigh` and `eig` on the one-sided Hermiticity proof, so a Hermitian operator the proof cannot see is diagonalized as a general matrix and the "disentangled" basis comes back non-orthonormal by an O(1), roundoff-dependent amount (2.7e-02 to 7.9e-01 measured) whenever the operator is degenerate on the manifold

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `operators`

**Status**: FIXED. `disentangle_manifold` decides on `wfs[0].MBO.is_hermitian(A)`, the chain's probe (proof first, then the random witness), falling back to the bare proof when a state carries no chain (`disentangle._is_hermitian`), and on the Hermitian branch takes `eigh` of `(ma+ma^dagger)/2`; the `eig` branch stays for a non-Hermitian operator, which is bit-identical. Pinned by `tests/test_audit_2026_09_24b_misc.py::test_unproven_hermitian_operator_gives_an_orthonormal_eigenbasis` (spin, fermion, parafermion: Gram error below 1e-12, exact levels, A diagonal on the output), `::test_bare_proof_would_have_taken_eig`, `::test_proven_operator_output_is_unchanged` (columns up to a phase, or eigenprojectors, to 1e-12), `::test_non_hermitian_operator_still_takes_eig` and `::test_states_without_a_chain_fall_back_to_the_bare_proof`. NUMBERS CHANGE: for `C0*Adag1+Cdag0*A1` on the hunter's 2-site fermion manifold the Gram error goes from 3.237e-01 to 1.177e-15 and the level error from 5.640e-02 to 1.665e-15; for `1j*Sx0*Sy0` from 3.546e-01 to 2.220e-15 and from 2.042e-01 to 1.554e-15; for `Sig0+SigDag0` on a 2-site Z3 chain from 1.65e-01 to 1.87e-14. Because the whole Hermitian branch is now symmetrized, a proven operator moves too, at the gauge level: a phase flip on one state of `0.3*N0+0.1*N1`, a rotation inside the degenerate eigenspace of `-0.5*Sz0` (eigenprojectors unchanged), and on truncated Heisenberg manifolds a correction at the truncation level for `-0.5*Sz_tot` (eigenvalues by 9.8e-5 at L=8, `maxm=6`), since plain `eigh` read one triangle of an `ma` that is Hermitian only to that level; the output of `0.3*Sz0+0.1*Sz1` is unchanged to 1e-15. See findings 12 and 18 of `audit_2026_09_24c_hole_hunt.md`: an ED state carries its `EDchain` as `MBO`, so it never reached the bare-proof fallback and raised `AttributeError` instead, and the chain's probe this fix relies on compared against an absolute 1e-4, so a small non-Hermitian operator was diagonalized as its Hermitian part; `EDchain.is_hermitian` exists now and the probe is scale-free.

**Where**: `src/dmrgpy/mpsalgebratk/disentangle.py:12` (the branch on the bare
`A.is_hermitian()`), reached only as `dmrgpy.mpsalgebra.disentangle_manifold`
(the import at `mpsalgebra.py:422`); the ED twin
`algebra/algebra.py::disentangle_manifold` takes `eigh` unconditionally.

The function represents `A` on a manifold of states, `ma = <w_i|A|w_j>`, and
diagonalizes `ma`; `scipy.linalg.eig` does not orthogonalize inside a degenerate
eigenspace, `eigh` does. The branch is decided by `MultiOperator.is_hermitian()`,
the canonical-form proof, which is one-sided by design: "not proven" is not
"false", which is why every other Hermiticity gate in `src/` goes through
`Many_Body_Chain.is_hermitian`, proof first and a random-witness probe after it.
This is the only branch in `src/` decided by the bare proof. Three families of
Hermitian operator reach the `eig` branch: same-site identities (`1j*Sx0*Sy0`,
exactly `-Sz0/2`), terms mixing a C-type and an A-type name, whose proof
`30200a4` withdrew (the previous record's finding 1), and every operator on a
`Parafermionic_Chain`, whose names are off `_PARITY`. The basis still spans the
manifold and still consists of eigenvectors of `ma`; it is the orthonormality that
is lost. It survived because the natural disentangling operators (`Sz`, `N`) are
diagonal and proven, because `eig` of a Hermitian matrix with distinct eigenvalues
is orthogonal to roundoff, and because nothing in the tree calls the function
(last touched in `fe3580e`, before `canonical.py` existed).

**Expected**: an orthonormal eigenbasis of `ma`, since `ma` is Hermitian (to
4e-16 in every case below).

Repro, from the hunter (the manifold is orthonormal to 1.6e-15; "state is not
normalizable" comes from `get_excited_states`' purify step dropping two surplus
states on a 4-dimensional space):

```bash
cd <scratch>/operators && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 08_disentangle_unproven.py
```

`operators/08_disentangle_unproven.py`:

```python
"""mpsalgebratk/disentangle.py picks eigh or eig on MultiOperator.is_hermitian(),
the one-sided symbolic proof, so a Hermitian operator the proof cannot see
is diagonalized as a general matrix, and eig does not orthogonalize inside a
degenerate eigenspace. Manifold: the four eigenstates of a generic 2-site
Hamiltonian (full bond dimension, orthonormal to 1e-15). The disentangling
operator has a two-fold degenerate spectrum on it.

Anchors: (i) the representation matrix ma = <w_i|A|w_j> is Hermitian to
roundoff, so its eigenbasis is orthonormal and the output of the step must
be an orthonormal basis; (ii) the output spans the same space as the input,
so the Hamiltonian's representation on it, diagonalized, must give back the
input energies es exactly; (iii) the same call with A spelled so that the
proof lands.

Part 1, the 30200a4 half: T = C0*Adag1 + Cdag0*A1 (= A0 Adag1 + Adag0 A1, a
hopping, spectrum -1,0,0,1), proven Hermitian before the fix and refused
after it; "old grading" stubs _mixes_representations for that one call.
Part 2, pre-existing: A = 1j*Sx0*Sy0, exactly -Sz0/2 on spin-1/2, whose
Hermiticity rests on a same-site identity the canonical form does not model."""
import numpy as np, warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain, fermionchain, mpsalgebra
from dmrgpy.mpsalgebratk.disentangle import get_representation
from dmrgpy.multioperatortk import canonical

def run(label, chain, wfs, es, A):
    ma = get_representation(wfs, A)
    out = mpsalgebra.disentangle_manifold(wfs, A)
    G = np.array([[a.dot(b) for b in out] for a in out])
    Hrep = get_representation(out, chain.hamiltonian)
    eH = np.sort(np.linalg.eigvals(Hrep).real)
    print("%-26s A.is_hermitian()=%-5s ||ma-ma^H||=%.1e (||ma||=%.3f) | "
          "max|<out_i|out_j>-d_ij|=%.3e | eig(H on out)=%s vs es=%s, max err %.3e"
          % (label, A.is_hermitian(), np.linalg.norm(ma-ma.conj().T), np.linalg.norm(ma),
             np.max(np.abs(G-np.eye(len(out)))), np.round(eH, 6), np.round(np.real(es), 6),
             np.max(np.abs(eH-np.sort(np.real(es))))))

# Part 1
fc = fermionchain.Fermionic_Chain(2, itensor_version="python")
hf = 0.8*fc.Cdag[0]*fc.C[1] + 0.8*fc.Cdag[1]*fc.C[0] + 0.5*(fc.C[0]*fc.C[1] + fc.Cdag[1]*fc.Cdag[0]) \
     + 0.3*fc.N[0] - 0.4*fc.N[1]
fc.set_hamiltonian(hf); fc.maxm = 8; fc.nsweeps = 20
np.random.seed(0)
es, wfs = fc.get_excited_states(n=4)
T = fc.C[0]*fc.Adag[1] + fc.Cdag[0]*fc.A[1]
run("T, after 30200a4", fc, wfs, es, T)
orig = canonical._mixes_representations
canonical._mixes_representations = lambda f, p: False
try: run("T, old grading", fc, wfs, es, T)
finally: canonical._mixes_representations = orig
run("A0*Adag1+Adag0*A1 (same T)", fc, wfs, es, fc.A[0]*fc.Adag[1] + fc.Adag[0]*fc.A[1])

# Part 2
sc = spinchain.Spin_Chain([2]*2, itensor_version="python")
h = sc.Sx[0]*sc.Sx[1] + 0.6*sc.Sy[0]*sc.Sy[1] + 0.3*sc.Sz[0]*sc.Sz[1] \
    + 0.7*sc.Sx[0] + 0.45*sc.Sz[1] + 0.2*sc.Sy[0]
sc.set_hamiltonian(h); sc.maxm = 8; sc.nsweeps = 20
np.random.seed(0)
es, wfs = sc.get_excited_states(n=4)
S = np.array([[a.dot(b) for b in wfs] for a in wfs])
print("spin manifold: es=%s, max|<w_i|w_j>-d_ij| = %.1e" % (np.round(np.real(es), 8), np.max(np.abs(S-np.eye(4)))))
run("1j*Sx0*Sy0", sc, wfs, es, 1j*sc.Sx[0]*sc.Sy[0])
run("-0.5*Sz0 (same operator)", sc, wfs, es, -0.5*sc.Sz[0])
```

`operators/08_disentangle_unproven.out`:

```
WARNING, state is not normalizable. Returning None
WARNING, state is not normalizable. Returning None
T, after 30200a4           A.is_hermitian()=False ||ma-ma^H||=4.4e-16 (||ma||=1.414) | max|<out_i|out_j>-d_ij|=3.237e-01 | eig(H on out)=[-0.923212 -0.565039  0.396097  0.823212] vs es=[-0.923212 -0.552494  0.452494  0.823212], max err 5.640e-02
T, old grading             A.is_hermitian()=True  ||ma-ma^H||=4.4e-16 (||ma||=1.414) | max|<out_i|out_j>-d_ij|=1.530e-15 | eig(H on out)=[-0.923212 -0.552494  0.452494  0.823212] vs es=[-0.923212 -0.552494  0.452494  0.823212], max err 1.887e-15
A0*Adag1+Adag0*A1 (same T) A.is_hermitian()=True  ||ma-ma^H||=4.4e-16 (||ma||=1.414) | max|<out_i|out_j>-d_ij|=1.530e-15 | eig(H on out)=[-0.923212 -0.552494  0.452494  0.823212] vs es=[-0.923212 -0.552494  0.452494  0.823212], max err 1.887e-15
WARNING, state is not normalizable. Returning None
WARNING, state is not normalizable. Returning None
spin manifold: es=[-0.75062831 -0.10217921  0.14580064  0.70700688], max|<w_i|w_j>-d_ij| = 1.6e-15
1j*Sx0*Sy0                 A.is_hermitian()=False ||ma-ma^H||=3.1e-16 (||ma||=0.500) | max|<out_i|out_j>-d_ij|=3.546e-01 | eig(H on out)=[-0.546406 -0.118028  0.140409  0.743598] vs es=[-0.750628 -0.102179  0.145801  0.707007], max err 2.042e-01
-0.5*Sz0 (same operator)   A.is_hermitian()=True  ||ma-ma^H||=3.1e-16 (||ma||=0.500) | max|<out_i|out_j>-d_ij|=2.665e-15 | eig(H on out)=[-0.750628 -0.102179  0.145801  0.707007] vs es=[-0.750628 -0.102179  0.145801  0.707007], max err 1.471e-15
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced to the digit, and showed the
mechanism needs no dmrgpy (a random exactly-Hermitian 4x4 with spectrum
[-1,0,0,1] gives an `eig` Gram error of 0.32 to 0.53 over six draws, `eigh`
5e-16). Found the third family: on a 2-site Z3 chain `Sig0+SigDag0` gives 0.165
and `Sig0*SigDag1+Sig1*SigDag0` 0.476. On realistic input, the singlet plus
triplet of an open Heisenberg chain with `A = sum_i 1j*Sx_i*Sy_i = -Sz_tot/2`,
the provable spelling gives the identical `ma` and comes out orthonormal, so the
proof decides, not the matrix:

`reviews/operators_O1/04_submanifold.py`:

```python
"""Does the O(1) Gram error survive on the intended input, a low-lying
SUB-manifold of a longer chain, rather than the whole Hilbert space of a
2-site chain? On a sub-manifold the degeneracy of ma = <w_i|A|w_j> is exact
only to the accuracy of the manifold's span, and eig's non-orthogonality
inside a split pair scales like roundoff/split.

(a) open Heisenberg chain, L sites: singlet + triplet (n=4),
    A = sum_i 1j*Sx_i*Sy_i = -Sz_tot/2, spectrum on the manifold {-1/2,0,0,1/2}.
(b) Z3 clock chain in its ordered phase: the three quasi-degenerate ground
    states (n=3), A = Sig0+SigDag0, spectrum about m*{2,-1,-1}.
Seed 0 for the manifold; sweeps pinned.
"""
import numpy as np, warnings, scipy.linalg as dlg
warnings.filterwarnings("ignore")
from dmrgpy import spinchain, parafermionchain, mpsalgebra
from dmrgpy.mpsalgebratk.disentangle import get_representation

def gerr(G): return np.max(np.abs(G-np.eye(len(G))))
def gram(ws): return np.array([[a.dot(b) for b in ws] for a in ws])

def report(label, chain, wfs, es, A, n):
    eED = np.sort(np.real(chain.get_excited(n=n, mode="ED")))
    ma = get_representation(wfs, A)
    ev = np.sort(np.linalg.eigvalsh((ma+ma.conj().T)/2))
    split = np.min(np.diff(ev))
    _, Ve = dlg.eig(ma)
    out = mpsalgebra.disentangle_manifold(wfs, A)
    print("%-34s manifold Gram %.1e, max|es-ED| %.1e | proof=%s | spec(ma)=%s, closest pair split %.1e"
          % (label, gerr(gram(wfs)), np.max(np.abs(np.sort(np.real(es))-eED)), A.is_hermitian(),
             np.round(ev, 5), split))
    print("%-34s Gram(out) err %.2e, Gram(eig coeffs) err %.2e, roundoff/split %.1e"
          % ("", gerr(gram(out)), gerr(Ve.conj().T@Ve), 1e-16*np.linalg.norm(ma)/split))

for L, maxm in ((8, 16), (8, 6), (12, 12)):
    sc = spinchain.Spin_Chain([2]*L, itensor_version="python")
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1] for i in range(L-1))
    sc.set_hamiltonian(h); sc.maxm = maxm; sc.nsweeps = 20
    np.random.seed(0)
    es, wfs = sc.get_excited_states(n=4)
    A = sum(1j*sc.Sx[i]*sc.Sy[i] for i in range(L))
    report("Heisenberg L=%d maxm=%d" % (L, maxm), sc, wfs, es, A, 4)

for L, maxm in ((5, 9), (5, 4)):
    pc = parafermionchain.Parafermionic_Chain(L, itensor_version="python")
    h = -sum(pc.Sig[i]*pc.Sigd[i+1] + pc.Sig[i+1]*pc.Sigd[i] for i in range(L-1)) \
        - 0.2*sum(pc.Tau[i] + pc.Taud[i] for i in range(L))
    pc.set_hamiltonian(h); pc.maxm = maxm; pc.nsweeps = 20
    np.random.seed(0)
    es, wfs = pc.get_excited_states(n=3)
    report("Z3 clock L=%d maxm=%d" % (L, maxm), pc, wfs, es, pc.Sig[0] + pc.Sigd[0], 3)
```

`reviews/operators_O1/04_submanifold.out`:

```
Heisenberg L=8 maxm=16             manifold Gram 6.1e-15, max|es-ED| 1.2e-14 | proof=False | spec(ma)=[-0.5 -0.   0.   0.5], closest pair split 2.5e-16
                                   Gram(out) err 7.89e-01, Gram(eig coeffs) err 7.89e-01, roundoff/split 2.8e-01
Heisenberg L=8 maxm=6              manifold Gram 4.4e-04, max|es-ED| 7.5e-03 | proof=False | spec(ma)=[-4.9929e-01 -1.2000e-04  3.4000e-04  4.9919e-01], closest pair split 4.6e-04
                                   Gram(out) err 1.21e-01, Gram(eig coeffs) err 1.21e-01, roundoff/split 1.5e-13
Heisenberg L=12 maxm=12            manifold Gram 2.3e-06, max|es-ED| 5.7e-05 | proof=False | spec(ma)=[-0.5     -0.      -0.       0.49999], closest pair split 1.4e-06
                                   Gram(out) err 3.66e-01, Gram(eig coeffs) err 3.66e-01, roundoff/split 5.0e-11
Z3 clock L=5 maxm=9                manifold Gram 6.3e-15, max|es-ED| 2.1e-14 | proof=False | spec(ma)=[-0.98478 -0.98476  1.96954], closest pair split 2.5e-05
                                   Gram(out) err 2.54e-10, Gram(eig coeffs) err 2.25e-11, roundoff/split 9.8e-12
Z3 clock L=5 maxm=4                manifold Gram 8.0e-05, max|es-ED| 6.4e-04 | proof=False | spec(ma)=[-0.98482 -0.9847   1.96952], closest pair split 1.3e-04
                                   Gram(out) err 7.95e-05, Gram(eig coeffs) err 4.57e-12, roundoff/split 1.9e-12
```

Struck: the hunter's "the Hamiltonian represented on that basis misses the exact
levels by up to 5.6e-02 (and -0.5464 against -0.7506)" is not an independent
defect, since the generalized problem `eig(Hrep, G)` gives the exact levels back
to 1.4e-15 (the span is preserved); `30200a4` as the cause, since the branch has
misrouted unproven Hermitian operators since it was written, and under the false
proof `30200a4` removed it was worse (on the anti-Hermitian `X =
C1*Adag0+Cdag1*A0` the old grading took `eigh` and returned states with an
eigen-residual of 1.00, i.e. not eigenvectors of X at all); and "seed 0 gives
3.2e-01" as a size, replaced by "O(1), draw-dependent".

**Suggested fix**: decide on `wfs[0].MBO.is_hermitian(A)`, the chain's probe,
falling back to the bare proof when `MBO` is `None`, and take `eigh` of
`(ma+ma^dagger)/2`, since on a truncated manifold `ma` is non-Hermitian at the
truncation level. The hunter's alternative, `np.allclose(ma, ma^dagger)`, is
wrong on exactly that input: on L=8 at `maxm=6` and L=12 at `maxm=12` it sends
the Hermitian `A` back to `eig` (0.121 and 0.366), where the probe gives `eigh`
at 2e-16 to 7e-16 (`reviews/operators_O1/05_truncated.py`). The `eig` branch must
stay, for a non-normal `ma`. Returned states change only for a Hermitian `A` the
proof cannot see.

### 2. Since `765b537`, `mode="ED"` KPM rescales like the DMRG routes but carries no moment-divergence check, so whenever the ground state or the top of the band falls outside [-1,1] and carries weight it returns an unflagged wrong spectrum: a spurious 0.80 peak where the exact density is zero at `kpm_scale=0.49`, integrals of 102.6, 2.84e4 and -4.69e6 against a sum rule of 0.25 further down

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `kpm` &middot; origin `765b537`

**Status**: FIXED. `algebra/kpm.py` gains `KPM_MOMENT_BOUND_FACTOR = 1.5` and a shared `check_kpm_moment(mu, bound)`, and `get_moments_vivj_python` checks every moment against 1.5*||vi||*||vj||, raising the DMRG loops' `RuntimeError` message; both of its callers, `edtk/dynamics.py`'s Hermitian KPM branch and `edtk/distribution.py`, go through `dm_vivj_energy`, and the non-Hermitian KPM uses `get_mu_n_nh`, which the bound does not reach. The "sliver" just below 1/2 that no moment bound sees is left open, as the reviewer measured it. Pinned by `tests/test_audit_2026_09_24b_kpm.py::test_ed_kpm_below_half_raises_where_it_returned_garbage` (0.49 and 0.45), `::test_ed_kpm_below_half_still_exact_without_elastic_weight` and `::test_ed_kpm_at_the_default_scale_is_unchanged`. No returned number of a correct run changes: on the field chain at `kpm_scale` 0.49 (integral 0.2467, a spurious 0.80 peak) and 0.45 (integral 102.6) ED now raises, the fieldless chain at 0.45 still returns 0.24997, 0.7 and 0.55 are bit-identical, and ED `get_distribution` with too small a `scale` (integrals of -8.07e+07 at 0.9 and 5.24e+15 at 0.5 on X=Sx0+Sx1) now raises.

**Where**: `src/dmrgpy/algebra/kpm.py::get_moments_vivj_python` (the ED moment
recursion, no check), reached from `edtk/dynamics.py::dynamical_correlator_kpm`
and `edtk/distribution.py::distribution_kpm`; the DMRG loops it should match are
`pyitensor/chain.py::_check_kpm_moment` and `check_kpm_moment` in both
`chain_session.h` (finding 3 is about their own threshold). Reached through
`mode="ED"`, `sc.mode="ED"`, and `mode.py`'s `ns<3` fallback on
`itensor_version=3` even when the call says `mode="DMRG"`.

Since O2 the ED route places E0 and Emax at x = -+1/(2*kpm_scale), the DMRG
rescaling, so below `kpm_scale=1/2` both sit outside [-1,1], where T_k grows as
cosh(k*arccosh|x|). Before `765b537` the ED window was 3*max(|E0|,|Emax-E0|),
safe whatever `kpm_scale` says, and ED did not read `kpm_scale` at all. The user
guide's KPM entry promises that "the moment recursion detects the resulting
exponential divergence and aborts with an explicit error rather than returning a
silently wrong spectrum", in the general `submode="KPM"` paragraph and with no
backend named; the "safe ~0.5 floor" (`manybodychain.py:261-264`,
`user_guide.md:1634`) is stated and not enforced on ED; and
`known_issue_kpm_energy_truncation_window.md`'s "the ED path is unaffected",
written about energy truncation, was true before `765b537` and is not now. It
survived because every ED KPM test runs at the default `kpm_scale=0.7`, and the
only test below 0.5, `test_kpm_divergence_guard_catchable.py`, runs `mode="DMRG"`.

**Expected**: a raise, as the DMRG loops do below about 0.48 on this chain, or a
correct spectrum.

Repro, from the hunter, on a 4-site S=1/2 Heisenberg chain plus 0.2*Sz_0, pair
(Sz_0,Sz_0), delta=0.1:

```bash
cd <scratch>/kpm && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 05_ed_kpm_scale_below_half.py
```

`kpm/05_ed_kpm_scale_below_half.py`:

```python
"""Probe 05: mode="ED" submode="KPM" at kpm_scale below 1/2.

Since O2 the ED route rescales exactly like the DMRG ones, so the ground
state sits at x0 = -1/(2*kpm_scale), OUTSIDE [-1,1] once kpm_scale < 1/2.
The DMRG routes carry a moment-divergence check (check_kpm_moment /
_check_kpm_moment); the ED route (algebra/kpm.py::get_moments_vivj_python)
has none. Hypothesis: ED returns a finite, unflagged, wrong spectrum there,
while DMRG raises. Anchor: the exact Lehmann sum <0|A delta(w-H+E0) B|0>
from np.linalg.eigh, whose weight is <0|AB|0> exactly, plus the same
Lehmann sum broadened with a Gaussian of the calibrated FWHM (only for
orientation; the sum rule and the peak magnitude are the discriminants)."""
import numpy as np
np.random.seed(5)
from dmrgpy import spinchain


def build(v="python", ns=4):
    np.random.seed(5)
    sc = spinchain.Spin_Chain(["S=1/2"]*ns, itensor_version=v)
    h = 0
    for i in range(ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 30
    return sc


def lehmann(sc, A, B):
    ed = sc.get_ED_obj()
    H = np.array(ed.get_hamiltonian().todense())
    e, V = np.linalg.eigh(H); gs = V[:, 0]
    Am = np.array(ed.MO2matrix(A).todense()); Bm = np.array(ed.MO2matrix(B).todense())
    w = (np.conjugate(V.T)@(np.conjugate(Am.T)@gs)).conj()*(np.conjugate(V.T)@(Bm@gs))
    return e-e[0], w, e


es = np.linspace(-0.5, 4.0, 901)
delta = 0.1
sc0 = build()
om, w, e = lehmann(sc0, sc0.Sz[0], sc0.Sz[0])
W = e[-1]-e[0]
print("exact: W=%.6f  sum rule <0|Sz0 Sz0|0> = %.6f  weight at omega=0 = %.6f"
      % (W, np.sum(w).real, np.sum(w[np.abs(om) < 1e-9]).real))
for ks in (0.7, 0.55, 0.5, 0.45, 0.4, 0.3):
    for mode, v in (("ED", "python"), ("DMRG", "python"), ("DMRG", 3)):
        sc = build(v)
        sc.kpm_scale = ks
        try:
            x, y = sc.get_dynamical_correlator(mode=mode, submode="KPM",
                                               name=(sc.Sz[0], sc.Sz[0]), delta=delta, es=es)
            y = np.real(np.asarray(y))
            res = ("returns: integral=%.4e (exact 0.2500)  max|y|=%.3e  finite=%s"
                   % (np.trapezoid(y, x), np.max(np.abs(y)), np.all(np.isfinite(y))))
        except Exception as ex:
            res = "RAISES %s: %s" % (type(ex).__name__, str(ex)[:70])
        print("kpm_scale=%.2f  x0=%+.3f  %-4s v=%-7r %s" % (ks, -1/(2*ks), mode, v, res), flush=True)

print("2-site itensor_version=3 chain, no mode= (mode.py falls back to ED),"
      " kpm_scale=0.3 with kpm_energy_truncate=True as in the user guide's example")
sc = build(3, ns=2)
om2, w2, e2 = lehmann(sc, sc.Sz[0], sc.Sz[0])
sc.kpm_scale = 0.3
sc.kpm_energy_truncate = True
print("  get_mode() ->", sc.get_mode(), "  exact sum rule %.6f" % np.sum(w2).real)
try:
    x, y = sc.get_dynamical_correlator(submode="KPM", name=(sc.Sz[0], sc.Sz[0]), delta=delta, es=es)
    y = np.real(np.asarray(y))
    print("  returns: integral=%.4e  max|y|=%.3e  finite=%s" % (np.trapezoid(y, x), np.max(np.abs(y)), np.all(np.isfinite(y))))
except Exception as ex:
    print("  RAISES %s: %s" % (type(ex).__name__, str(ex)[:90]))
```

`kpm/05_ed_kpm_scale_below_half.out`:

```
exact: W=2.478269  sum rule <0|Sz0 Sz0|0> = 0.250000  weight at omega=0 = 0.014530
kpm_scale=0.70  x0=-0.714  ED   v='python' returns: integral=2.5000e-01 (exact 0.2500)  max|y|=7.682e-01  finite=True
kpm_scale=0.70  x0=-0.714  DMRG v='python' returns: integral=2.5000e-01 (exact 0.2500)  max|y|=7.680e-01  finite=True
kpm_scale=0.70  x0=-0.714  DMRG v=3       returns: integral=2.5000e-01 (exact 0.2500)  max|y|=7.683e-01  finite=True
kpm_scale=0.55  x0=-0.909  ED   v='python' returns: integral=2.4969e-01 (exact 0.2500)  max|y|=7.966e-01  finite=True
kpm_scale=0.55  x0=-0.909  DMRG v='python' returns: integral=2.4999e-01 (exact 0.2500)  max|y|=7.961e-01  finite=True
kpm_scale=0.55  x0=-0.909  DMRG v=3       returns: integral=2.4999e-01 (exact 0.2500)  max|y|=7.962e-01  finite=True
kpm_scale=0.50  x0=-1.000  ED   v='python' returns: integral=2.3547e-01 (exact 0.2500)  max|y|=8.255e-01  finite=True
kpm_scale=0.50  x0=-1.000  DMRG v='python' returns: integral=2.3582e-01 (exact 0.2500)  max|y|=8.251e-01  finite=True
kpm_scale=0.50  x0=-1.000  DMRG v=3       returns: integral=2.3582e-01 (exact 0.2500)  max|y|=8.251e-01  finite=True
kpm_scale=0.45  x0=-1.111  ED   v='python' returns: integral=1.0260e+02 (exact 0.2500)  max|y|=8.377e+03  finite=True
kpm_scale=0.45  x0=-1.111  DMRG v='python' RAISES RuntimeError: KPM moments diverging: scaled spectrum outside [-1,1] (band-edge estim
kpm_scale=0.45  x0=-1.111  DMRG v=3       RAISES RuntimeError: KPM moments diverging: scaled spectrum outside [-1,1] (band-edge estim
kpm_scale=0.40  x0=-1.250  ED   v='python' returns: integral=2.8415e+04 (exact 0.2500)  max|y|=2.910e+06  finite=True
kpm_scale=0.40  x0=-1.250  DMRG v='python' RAISES RuntimeError: KPM moments diverging: scaled spectrum outside [-1,1] (band-edge estim
kpm_scale=0.40  x0=-1.250  DMRG v=3       RAISES RuntimeError: KPM moments diverging: scaled spectrum outside [-1,1] (band-edge estim
kpm_scale=0.30  x0=-1.667  ED   v='python' returns: integral=-4.6948e+06 (exact 0.2500)  max|y|=4.320e+08  finite=True
kpm_scale=0.30  x0=-1.667  DMRG v='python' RAISES RuntimeError: KPM moments diverging: scaled spectrum outside [-1,1] (band-edge estim
kpm_scale=0.30  x0=-1.667  DMRG v=3       RAISES RuntimeError: KPM moments diverging: scaled spectrum outside [-1,1] (band-edge estim
2-site itensor_version=3 chain, no mode= (mode.py falls back to ED), kpm_scale=0.3 with kpm_energy_truncate=True as in the user guide's example
  get_mode() -> ED   exact sum rule 0.250000
  returns: integral=7.2673e+01  max|y|=4.739e+03  finite=True
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced digit for digit; ran the Python
tree of `765b537^` from a `git archive` extract, where ED returns 0.249999 at
every `kpm_scale` from 0.7 down to 0.3, and `git log -S "half =
width*kpm_scale"` names `765b537` alone, so this is a regression of the O2 fix,
inside the previous hunt's window and absent from its record (whose "Ruled out"
checked that `kpm_scale` is pushed to ED, not that the push is safe):

`reviews/kpm_A/03_pre765_ed.out`:

```
  import numpy as np, dmrgpy
dmrgpy from <scratch>/reviews/kpm_A/pre765/src/dmrgpy/__init__.py
pre-765b537 ED kpm_scale=0.70 integral=0.249999 max|y|=1.5523
pre-765b537 ED kpm_scale=0.50 integral=0.249999 max|y|=1.5523
pre-765b537 ED kpm_scale=0.45 integral=0.249999 max|y|=1.5523
pre-765b537 ED kpm_scale=0.40 integral=0.249999 max|y|=1.5523
pre-765b537 ED kpm_scale=0.30 integral=0.249999 max|y|=1.5523
```

Narrowed on the mechanism: below 1/2 the answer goes wrong only when a pole
outside the window carries weight, either real weight (here the elastic weight
<Sz_0>^2 = 0.0145 at E0) or roundoff grown past about n*arccosh(1/(2 kpm_scale))
= 75; on the fieldless chain, where Sz_0|gs> has no weight at either end, ED is
exact to 6e-14 at 0.45 and 0.40. And the dangerous rows are the quiet ones just
below 1/2, not the conspicuous 1e2 to 1e6 ones: at 0.49 the sum rule reads
0.2467, one per cent off, while the curve carries a 0.80 peak at omega=0.10
where the exact density is zero (97 per cent of the true peak), and the DMRG
routes are silent there too (finding 3):

`reviews/kpm_A/04_silent_band_field.py`:

```python
"""Review probe 04: the hunter's field chain (4-site Heisenberg + 0.2*Sz_0),
scanned finely just below kpm_scale=1/2. For each kpm_scale: the ED KPM
integral and max|y|, max|mu_k|/(||vi|| ||vj||) (<= 1 exactly when every
weighted pole is inside [-1,1]), whether the DMRG loops' bound
|mu_k| > 1e3*(||vi|| ||vj||+1) would fire, max|y_ED - y_ref| against the
same pipeline fed the exact in-window moments (so the out-of-window pole
simply dropped), and what v3 and "python" DMRG do at the same settings."""
import numpy as np
from reviewlib import build, spectrum, exact_moments, Spy

es = np.linspace(-0.5, 4.0, 901); delta = 0.1
sc0 = build()
e, w, nrm = spectrum(sc0, sc0.Sz[0], sc0.Sz[0])
print("W=%.6f  sum rule %.6f  weight at omega=0 %.6f  ||vi||*||vj||=%.4f"
      % (e[-1]-e[0], np.sum(w).real, np.sum(w[np.abs(e-e[0]) < 1e-9]).real, nrm))
for ks in (0.52, 0.50, 0.495, 0.49, 0.485, 0.48, 0.47, 0.46, 0.45):
    sc = build(); sc.kpm_scale = ks
    with Spy() as sp:
        x, y = sc.get_dynamical_correlator(mode="ED", submode="KPM",
                    name=(sc.Sz[0], sc.Sz[0]), delta=delta, es=es)
        mus = sp.rec[-1]; n = len(mus)
        ref, xn = exact_moments(e, w, ks, n)
        sp.sub = lambda nn: ref
        sc2 = build(); sc2.kpm_scale = ks
        _, yref = sc2.get_dynamical_correlator(mode="ED", submode="KPM",
                    name=(sc2.Sz[0], sc2.Sz[0]), delta=delta, es=es)
    y = np.real(np.asarray(y)); yref = np.real(np.asarray(yref))
    ratio = np.max(np.abs(mus))/nrm
    fires = np.max(np.abs(mus)) > 1e3*(nrm+1.0)
    line = ("kpm_scale=%.3f x0=%+.4f n=%3d ED: integral=%.4e max|y|=%.3e "
            "max|mu|/(|vi||vj|)=%.3e 1e3-bound fires=%s | ref integral=%.4f ref peak=%.4f "
            "max|y-yref|/refpeak=%.3e"
            % (ks, -1/(2*ks), n, np.trapezoid(y, x), np.max(np.abs(y)), ratio, fires,
               np.trapezoid(yref, x), np.max(np.abs(yref)),
               np.max(np.abs(y-yref))/np.max(np.abs(yref))))
    print(line, flush=True)
    for v in ("python", 3):
        sc = build(v); sc.kpm_scale = ks
        try:
            x, yd = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                        name=(sc.Sz[0], sc.Sz[0]), delta=delta, es=es)
            yd = np.real(np.asarray(yd))
            print("    DMRG v=%-7r returns integral=%.4e max|y|=%.3e max|y-yref|/refpeak=%.3e"
                  % (v, np.trapezoid(yd, x), np.max(np.abs(yd)),
                     np.max(np.abs(yd-yref))/np.max(np.abs(yref))), flush=True)
        except Exception as ex:
            print("    DMRG v=%-7r RAISES %s: %s" % (v, type(ex).__name__, str(ex)[:40]), flush=True)
```

`reviews/kpm_A/04_silent_band_field.out`:

```
W=2.478269  sum rule 0.250000  weight at omega=0 0.014530  ||vi||*||vj||=0.2500
kpm_scale=0.520 x0=-0.9615 n= 48 ED: integral=2.3978e-01 max|y|=8.195e-01 max|mu|/(|vi||vj|)=1.000e+00 1e3-bound fires=False | ref integral=0.2398 ref peak=0.8195 max|y-yref|/refpeak=2.181e-14
    DMRG v='python' returns integral=2.4989e-01 max|y|=8.193e-01 max|y-yref|/refpeak=3.252e-01
    DMRG v=3       returns integral=2.4989e-01 max|y|=8.193e-01 max|y-yref|/refpeak=3.252e-01
kpm_scale=0.500 x0=-1.0000 n= 46 ED: integral=2.3547e-01 max|y|=8.255e-01 max|mu|/(|vi||vj|)=1.000e+00 1e3-bound fires=False | ref integral=0.2355 ref peak=0.8255 max|y-yref|/refpeak=1.937e-14
    DMRG v='python' returns integral=2.3582e-01 max|y|=8.251e-01 max|y-yref|/refpeak=5.762e-02
    DMRG v=3       returns integral=2.3582e-01 max|y|=8.251e-01 max|y-yref|/refpeak=5.762e-02
kpm_scale=0.495 x0=-1.0101 n= 45 ED: integral=2.3573e-01 max|y|=8.187e-01 max|mu|/(|vi||vj|)=1.473e+01 1e3-bound fires=False | ref integral=0.2355 ref peak=0.8182 max|y-yref|/refpeak=1.038e-01
    DMRG v='python' returns integral=2.3035e-01 max|y|=1.931e+00 max|y-yref|/refpeak=2.360e+00
    DMRG v=3       returns integral=2.3035e-01 max|y|=1.931e+00 max|y-yref|/refpeak=2.360e+00
kpm_scale=0.490 x0=-1.0204 n= 45 ED: integral=2.4665e-01 max|y|=8.358e-01 max|mu|/(|vi||vj|)=2.073e+02 1e3-bound fires=False | ref integral=0.2354 ref peak=0.8289 max|y-yref|/refpeak=9.686e-01
    DMRG v='python' returns integral=2.8310e-01 max|y|=5.672e+00 max|y-yref|/refpeak=6.843e+00
    DMRG v=3       returns integral=2.8310e-01 max|y|=5.672e+00 max|y-yref|/refpeak=6.843e+00
kpm_scale=0.485 x0=-1.0309 n= 44 ED: integral=2.9683e-01 max|y|=4.562e+00 max|mu|/(|vi||vj|)=1.246e+03 1e3-bound fires=False | ref integral=0.2354 ref peak=0.8215 max|y-yref|/refpeak=5.553e+00
    DMRG v='python' returns integral=4.3928e-01 max|y|=2.947e+01 max|y-yref|/refpeak=3.588e+01
    DMRG v=3       returns integral=4.3931e-01 max|y|=2.948e+01 max|y-yref|/refpeak=3.588e+01
kpm_scale=0.480 x0=-1.0417 n= 44 ED: integral=5.4676e-01 max|y|=2.300e+01 max|mu|/(|vi||vj|)=6.849e+03 1e3-bound fires=True | ref integral=0.2354 ref peak=0.8325 max|y-yref|/refpeak=2.762e+01
    DMRG v='python' RAISES RuntimeError: KPM moments diverging: scaled spectrum o
    DMRG v=3       RAISES RuntimeError: KPM moments diverging: scaled spectrum o
kpm_scale=0.470 x0=-1.0638 n= 43 ED: integral=3.5999e+00 max|y|=2.574e+02 max|mu|/(|vi||vj|)=8.837e+04 1e3-bound fires=True | ref integral=0.2354 ref peak=0.8365 max|y-yref|/refpeak=3.077e+02
    DMRG v='python' RAISES RuntimeError: KPM moments diverging: scaled spectrum o
    DMRG v=3       RAISES RuntimeError: KPM moments diverging: scaled spectrum o
kpm_scale=0.460 x0=-1.0870 n= 42 ED: integral=2.2188e+01 max|y|=1.725e+03 max|mu|/(|vi||vj|)=6.858e+05 1e3-bound fires=True | ref integral=0.2354 ref peak=0.8409 max|y-yref|/refpeak=2.052e+03
    DMRG v='python' RAISES RuntimeError: KPM moments diverging: scaled spectrum o
    DMRG v=3       RAISES RuntimeError: KPM moments diverging: scaled spectrum o
kpm_scale=0.450 x0=-1.1111 n= 41 ED: integral=1.0260e+02 max|y|=8.377e+03 max|mu|/(|vi||vj|)=3.788e+06 1e3-bound fires=True | ref integral=0.2354 ref peak=0.8458 max|y-yref|/refpeak=9.904e+03
    DMRG v='python' RAISES RuntimeError: KPM moments diverging: scaled spectrum o
    DMRG v=3       RAISES RuntimeError: KPM moments diverging: scaled spectrum o
```

(The DMRG "max|y-yref|" at 0.52 and 0.50 is the reference carrying ED's own
|x| <= 0.95 cutoff, a lead below, not a DMRG error.) Struck: "every DMRG route
raises at the same settings", which holds only at the sampled points; and the
realistic trigger as stated, since the user guide's truncation example as written
gives 1.85e3 on a 2-site v3 dimer, the 72.67 being a construction on it (the ED
half of that is finding 4).

**Suggested fix**: the exact Cauchy-Schwarz bound, not the DMRG loops' 1e3*(b+1):
raise in `get_moments_vivj_python` when max_k |mu_k| > c*||vi||*||vj|| with c
between 1.5 and 2, which also covers `get_distribution(mode="ED")` with too small a
`scale`. ED carries no MPS noise to allow for. Measured: every correct run sits at
a ratio of at most 1.000, every run with an error of 1.1 per cent or more at 2.2
or above, and the one casualty is a fieldless run at 0.45, delta=0.025, ratio 17,
error 2.6e-4 (one delta step further it is 1.1e4 times the peak); no moment bound
sees the sub-percent contamination at 0.4995 (ratio 1.000, error 1.75e-3). The
hunter's two alternatives are both wrong: the 1e3*(b+1) bound misses 0.485 to
0.495, and an up-front raise at 1/(2 kpm_scale) > 1 rejects the exact fieldless
answers at 0.40 to 0.49.

`reviews/kpm_A/09_damped_statistic.py`:

```python
"""Review probe 09: which moment statistic separates the right ED answers from
the wrong ones? For each case: max|mu_k|/(||vi|| ||vj||) (raw), the same with
the Jackson factors g_k the reconstruction applies (damped), and the pointwise
error against the exact in-window reference of reviewlib. A guard wants a
statistic that is <= 1 on every correct run and > 1 as soon as the error is
visible."""
import numpy as np
from reviewlib import build, spectrum, exact_moments, Spy
from dmrgpy.algebra import kpm

def run(field, ks, delta, es=np.linspace(-0.5, 4.0, 901)):
    sc = build(field=field); sc.kpm_scale = ks
    e, w, nrm = spectrum(sc, sc.Sz[0], sc.Sz[0])
    with Spy() as sp:
        x, y = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=(sc.Sz[0], sc.Sz[0]), delta=delta, es=es)
        mus = sp.rec[-1]; n = len(mus); ref, _ = exact_moments(e, w, ks, n); sp.sub = lambda nn: ref
        sc2 = build(field=field); sc2.kpm_scale = ks
        _, yr = sc2.get_dynamical_correlator(mode="ED", submode="KPM", name=(sc2.Sz[0], sc2.Sz[0]), delta=delta, es=es)
    y = np.real(y); yr = np.real(yr)
    g = kpm.jackson_kernel(np.ones(n))
    print("field=%.1f ks=%.4f delta=%.3f n=%3d raw=%.3e damped=%.3e  max|y-yref|/refpeak=%.3e"
          % (field, ks, delta, n, np.max(np.abs(mus))/nrm, np.max(np.abs(g*mus))/nrm,
             np.max(np.abs(y-yr))/np.max(np.abs(yr))), flush=True)

for c in ((0.2, 0.7, 0.1), (0.0, 0.45, 0.05), (0.0, 0.40, 0.05), (0.0, 0.45, 0.03),
          (0.2, 0.499, 0.1), (0.2, 0.498, 0.1), (0.2, 0.495, 0.1), (0.2, 0.49, 0.1), (0.2, 0.485, 0.1), (0.2, 0.48, 0.1),
          (0.0, 0.45, 0.025), (0.0, 0.40, 0.03), (0.0, 0.45, 0.02), (0.0, 0.35, 0.2), (0.0, 0.35, 0.1), (0.0, 0.30, 0.2)):
    run(*c)
```

`reviews/kpm_A/09_damped_statistic.out`:

```
field=0.2 ks=0.7000 delta=0.100 n= 64 raw=1.000e+00 damped=1.000e+00  max|y-yref|/refpeak=1.944e-14
field=0.0 ks=0.4500 delta=0.050 n= 79 raw=1.000e+00 damped=1.000e+00  max|y-yref|/refpeak=5.752e-14
field=0.0 ks=0.4000 delta=0.050 n= 70 raw=1.000e+00 damped=1.000e+00  max|y-yref|/refpeak=6.132e-14
field=0.0 ks=0.4500 delta=0.030 n=131 raw=1.000e+00 damped=1.000e+00  max|y-yref|/refpeak=2.292e-09
field=0.2 ks=0.4990 delta=0.100 n= 46 raw=1.243e+00 damped=1.000e+00  max|y-yref|/refpeak=4.847e-03
field=0.2 ks=0.4980 delta=0.100 n= 46 raw=2.359e+00 damped=1.000e+00  max|y-yref|/refpeak=1.106e-02
field=0.2 ks=0.4950 delta=0.100 n= 45 raw=1.473e+01 damped=1.000e+00  max|y-yref|/refpeak=1.038e-01
field=0.2 ks=0.4900 delta=0.100 n= 45 raw=2.073e+02 damped=1.566e+00  max|y-yref|/refpeak=9.686e-01
field=0.2 ks=0.4850 delta=0.100 n= 44 raw=1.246e+03 damped=6.039e+00  max|y-yref|/refpeak=5.553e+00
field=0.2 ks=0.4800 delta=0.100 n= 44 raw=6.849e+03 damped=2.336e+01  max|y-yref|/refpeak=2.762e+01
field=0.0 ks=0.4500 delta=0.025 n=158 raw=1.700e+01 damped=1.000e+00  max|y-yref|/refpeak=2.643e-04
field=0.0 ks=0.4000 delta=0.030 n=117 raw=3.277e+04 damped=1.045e+00  max|y-yref|/refpeak=6.077e-01
field=0.0 ks=0.4500 delta=0.020 n=197 raw=2.416e+09 damped=3.879e+04  max|y-yref|/refpeak=1.129e+04
field=0.0 ks=0.3500 delta=0.200 n= 16 raw=2.217e+00 damped=1.000e+00  max|y-yref|/refpeak=2.254e-01
field=0.0 ks=0.3500 delta=0.100 n= 31 raw=5.824e+02 damped=2.755e+00  max|y-yref|/refpeak=4.048e+00
field=0.0 ks=0.3000 delta=0.200 n= 16 raw=2.018e+02 damped=1.922e+00  max|y-yref|/refpeak=5.637e+00
```

### 3. The moment-divergence guard of every DMRG KPM route, `|mu_k| > 1e3*(||vi|| ||vj|| + 1)`, is 1e3 to 5e3 times looser than the exact bound and not scale-invariant, so v2, v3, `"python"` and `julia_live` return spectra up to 109 times the true peak without raising below `kpm_scale=1/2`, and scaling the operators by 1e-2 disables it over the whole band

`bug` &middot; severity **MEDIUM** for the +1 term, **LOW to MEDIUM** for the band below 1/2 &middot; CONFIRMED, NARROWED &middot; lens `kpm` (found by the reviewer of finding 2) &middot; origin `a67228e`, `695c452`

**Status**: FIXED, with a rebuild of both extensions. The threshold is 1.5*||vi||*||vj||, with no +1, in `pyitensor/chain.py::_check_kpm_moment` (now the shared check), in `check_kpm_moment` of both `chain_session.h` (the v3 copy covering the plain and the energy-truncated loops) and in `mpsjulialive/kpm.jl`; every accelerated loop now checks both moments it appends per step. `general_kpm` keeps the bound, since every in-tree X there is Hermitian and `scale_operator` places it in [-0.8,0.8]. The sliver below 1/2 stays open. One claim of the reviewer did not hold: truncation can push a moment past the bound when it is harsh, 10.05 at `kpmmaxm=2` and 2.60 at `kpmmaxm=3` on 4 sites at `kpm_scale=0.55`, on runs whose spectra are 683 and 55 per cent of the peak wrong, so the guard now raises there too (every run with an error below 10 per cent stays at a ratio of 1.0000 or less); its message then names the band edge and `kpm_scale` rather than `kpmmaxm`, a lead. `julia_live` raises `juliacall.JuliaError` rather than `RuntimeError`, as it did before. Pinned by `tests/test_audit_2026_09_24b_kpm.py::test_dmrg_kpm_raises_in_the_band_the_old_guard_let_through` (three `kpm_scale` values on `"python"`, v3 and v2), `::test_dmrg_kpm_full_recursion_raises_too` (with `julia_live`), `::test_dmrg_guard_is_scale_invariant`, `::test_dmrg_kpm_correct_runs_do_not_raise` and `::test_dmrg_kpm_truncated_correct_runs_do_not_raise` (8 sites, `kpmmaxm=8`, auto and cross pairs, acceleration on and off); before the rebuild the v2 and v3 cases failed, 10 of 15. No returned number of a correct run changes: at 0.495, 0.49 and 0.485, where `"python"` returned curves 1.931, 5.672 and 29.47 high (v3 1.879, 5.672, 29.49; v2 1.930, 5.666, 29.49), and at 0.45 with operators scaled by 1e-2, every backend now raises; 0.7 and 0.55 on `"python"` are identical, and v3 moves only by its documented run-to-run emax noise.

**Where**: `src/dmrgpy/pyitensor/chain.py::_check_kpm_moment`; `check_kpm_moment`
in `src/dmrgpy/mpscpp2/chain_session.h` and `src/dmrgpy/mpscpp3/chain_session.h`
(called from both the plain and the energy-truncated loops);
`src/dmrgpy/mpsjulialive/kpm.jl:22` (read, not run). In the accelerated loops only
`out[-1]` is checked, so the second moment appended per step never is.

The exact Cauchy-Schwarz bound is |mu_k| <= ||vi|| ||vj|| whenever the rescaled
spectrum lies inside [-1,1], since |T_k(x)| <= 1 there. The shipped threshold is
1e3 to 5e3 times looser, and the +1 turns it into an absolute threshold of about
1e3 whenever ||vi|| ||vj|| < 1. Below `kpm_scale=1/2` on the untruncated midpoint
window, E0 sits at x0 = -1/(2 kpm_scale) by construction; the DMRG backends then
compute exactly the right moments (they match the exact all-pole moments at the
run's own rescaling to 1e-12 on v2, v3 and `"python"`) and reconstruct a spectrum
with a spurious feature at the elastic line or the band top, often negative, and do
not raise until the ratio crosses the threshold. The band's edge is where the
out-of-window weight fraction times cosh(n*arccosh(1/(2 kpm_scale))) reaches
1e3(1+1/bound), so it moves with n and with the operator. The constants were
introduced in `a67228e` as a detector for an underestimated emax at 0.7 and kept
in `695c452` when `|mu0|` became `||vi|| ||vj||`; neither commit calibrates them,
and the comment's reason for the +1 ("meaningful when both norms are tiny") does not
hold, moments being bilinear in vi and vj. The regime is below the documented
"safe ~0.5 floor" and is reached only by a caller who lowers `kpm_scale` with
`kpm_energy_truncate=False`; the user guide's divergence promise quoted under
finding 2 carries no such qualifier.

**Expected**: a raise wherever the reconstructed spectrum is not the calibrated
one, as the user guide promises.

Evidence, from the reviewer who owned this candidate (v3 rows shown; `"python"`
and v2 print the same shape):

`reviews/kpm_A2/01_band_scan.py`:

```python
"""kpm_A2 probe 01: the silent band below kpm_scale=1/2, re-measured on the
DMRG backends themselves (v2, v3, "python"), pinned schedule (maxm=30,
nsweeps=30, numpy seeded), 4-site S=1/2 Heisenberg chain + 0.2*Sz_0, pair
(Sz_0, Sz_0), at delta = 0.2, 0.1, 0.05.

Per run: raw ratio max_k|mu_k|/(||vi|| ||vj||) over the n moments returned
(exact norms), fidelity of the DMRG moments against the exact all-pole
moments at the run's own rescaling, and the reconstructed spectrum against
(a) the same pipeline fed exact in-window moments (yref) and (b) a Lorentzian
of the exact spectrum. Part 2: the same pair scaled by 1e-2 on both operators,
to see the +1 in the threshold at work."""
import numpy as np
import a2lib as L

es = np.linspace(-0.5, 4.0, 901)
sc0 = L.heis("python")
e, w, bnd = L.exact_weights(sc0, sc0.Sz[0], sc0.Sz[0])
W = e[-1]-e[0]
print("W=%.6f sum rule=%.6f elastic weight=%.6f (fraction %.4f) ||vi|| ||vj||=%.6f"
      % (W, np.sum(w).real, np.sum(w[np.abs(e-e[0]) < 1e-9]).real,
         np.sum(w[np.abs(e-e[0]) < 1e-9]).real/np.sum(w).real, bnd))
for delta in (0.2, 0.1, 0.05):
    ylz = L.lorentz(e, w, es, delta)
    print("--- delta=%.2f  Lorentzian exact peak %.4f" % (delta, ylz.max()))
    for ks in (0.7, 0.52, 0.50, 0.499, 0.497, 0.495, 0.49, 0.485, 0.48, 0.47, 0.46):
        for v in ("python", 3, 2):
            sc = L.heis(v)
            sc.kpm_scale = ks
            tag = "delta=%.2f ks=%.3f v=%-6r" % (delta, ks, v)
            try:
                mus, emin, emax, scale, n, d = sc.get_dynamical_correlator_moments(
                    name=(sc.Sz[0], sc.Sz[0]), delta=delta)
            except Exception as ex:
                print(tag, "RAISES %s: %s" % (type(ex).__name__, str(ex)[:34]), flush=True)
                continue
            mus = np.asarray(mus)
            x0 = (e[0]-(emin+emax)/2.)*scale
            ref_all = L.exact_moments_allpoles(e, w, emin, emax, scale, n)
            ref_in = L.exact_moments_inwindow(e, w, emin, emax, scale, n)
            y = np.real(L.reconstruct(mus, emin, emax, scale, n, es, delta))
            yref = np.real(L.reconstruct(ref_in, emin, emax, scale, n, es, delta))
            ratio = np.max(np.abs(mus))/bnd
            fid = np.max(np.abs(mus-ref_all))/np.max(np.abs(ref_all))
            err = np.abs(y-yref)
            iw = np.argmax(err)
            print("%s n=%3d x0=%+.4f emax-err=%.1e ratio=%.3e fid=%.1e | max|y|=%.3f "
                  "refpeak=%.3f err/refpeak=%.3e at omega=%.3f | c=2 fires=%s"
                  % (tag, n, x0, (e[-1]-emax), ratio, fid, np.max(np.abs(y)), yref.max(),
                     err[iw]/yref.max(), es[iw], np.max(np.abs(mus)) > 2*bnd), flush=True)

print("=== part 2: pair (1e-2*Sz_0, 1e-2*Sz_0), delta=0.1")
eps = 1e-2
for ks in (0.48, 0.47, 0.46, 0.45):
    for v in ("python", 3):
        sc = L.heis(v)
        sc.kpm_scale = ks
        A = eps*sc.Sz[0]
        e2, w2, b2 = L.exact_weights(sc, A, A)
        tag = "ks=%.3f v=%-6r" % (ks, v)
        try:
            mus, emin, emax, scale, n, d = sc.get_dynamical_correlator_moments(
                name=(A, A), delta=0.1)
        except Exception as ex:
            print(tag, "RAISES %s: %s" % (type(ex).__name__, str(ex)[:34]), flush=True)
            continue
        mus = np.asarray(mus)
        ref_in = L.exact_moments_inwindow(e2, w2, emin, emax, scale, n)
        y = np.real(L.reconstruct(mus, emin, emax, scale, n, es, 0.1))
        yref = np.real(L.reconstruct(ref_in, emin, emax, scale, n, es, 0.1))
        print("%s bound=%.3e threshold=%.1f max|mu|=%.3e ratio=%.3e | integral y/eps^2=%.4e "
              "(exact 0.25) err/refpeak=%.3e" % (tag, b2, 1e3*(b2+1), np.max(np.abs(mus)),
              np.max(np.abs(mus))/b2, np.trapezoid(y, es)/eps**2,
              np.max(np.abs(y-yref))/yref.max()), flush=True)
```

`reviews/kpm_A2/01_band_scan.out`:

```
W=2.478269 sum rule=0.250000 elastic weight=0.014530 (fraction 0.0581) ||vi|| ||vj||=0.250000
--- delta=0.20  Lorentzian exact peak 0.2513
delta=0.20 ks=0.700 v='python' n= 32 x0=-0.7143 emax-err=-5.6e-16 ratio=1.000e+00 fid=1.9e-13 | max|y|=0.390 refpeak=0.390 err/refpeak=2.258e-13 at omega=-0.015 | c=2 fires=False
delta=0.20 ks=0.700 v=3      n= 32 x0=-0.7143 emax-err=9.0e-04 ratio=1.000e+00 fid=1.1e-13 | max|y|=0.390 refpeak=0.390 err/refpeak=1.254e-13 at omega=-0.015 | c=2 fires=False
delta=0.20 ks=0.700 v=2      n= 32 x0=-0.7143 emax-err=5.2e-04 ratio=1.000e+00 fid=3.3e-13 | max|y|=0.390 refpeak=0.390 err/refpeak=3.721e-13 at omega=-0.015 | c=2 fires=False
delta=0.20 ks=0.520 v='python' n= 24 x0=-0.9615 emax-err=-5.6e-16 ratio=1.000e+00 fid=2.2e-13 | max|y|=0.418 refpeak=0.418 err/refpeak=6.233e-13 at omega=-0.025 | c=2 fires=False
delta=0.20 ks=0.520 v=3      n= 24 x0=-0.9615 emax-err=1.5e-05 ratio=1.000e+00 fid=6.0e-14 | max|y|=0.418 refpeak=0.418 err/refpeak=8.422e-14 at omega=-0.025 | c=2 fires=False
delta=0.20 ks=0.520 v=2      n= 24 x0=-0.9615 emax-err=8.1e-04 ratio=1.000e+00 fid=3.4e-13 | max|y|=0.418 refpeak=0.418 err/refpeak=9.214e-13 at omega=-0.025 | c=2 fires=False
delta=0.20 ks=0.500 v='python' n= 23 x0=-1.0000 emax-err=-5.6e-16 ratio=1.000e+00 fid=2.1e-13 | max|y|=0.422 refpeak=0.422 err/refpeak=1.033e-12 at omega=0.015 | c=2 fires=False
delta=0.20 ks=0.500 v=3      n= 23 x0=-1.0000 emax-err=2.1e-05 ratio=1.000e+00 fid=2.1e-13 | max|y|=0.422 refpeak=0.422 err/refpeak=4.357e-13 at omega=0.015 | c=2 fires=False
delta=0.20 ks=0.500 v=2      n= 23 x0=-1.0000 emax-err=6.5e-05 ratio=1.000e+00 fid=1.7e-13 | max|y|=0.422 refpeak=0.422 err/refpeak=6.300e-13 at omega=0.015 | c=2 fires=False
delta=0.20 ks=0.499 v='python' n= 23 x0=-1.0020 emax-err=-5.6e-16 ratio=1.000e+00 fid=2.5e-13 | max|y|=0.423 refpeak=0.423 err/refpeak=6.830e-01 at omega=0.015 | c=2 fires=False
delta=0.20 ks=0.499 v=3      n= 23 x0=-1.0020 emax-err=8.5e-04 ratio=1.000e+00 fid=4.0e-13 | max|y|=0.423 refpeak=0.423 err/refpeak=6.829e-01 at omega=0.015 | c=2 fires=False
delta=0.20 ks=0.499 v=2      n= 23 x0=-1.0020 emax-err=3.4e-04 ratio=1.000e+00 fid=4.2e-13 | max|y|=0.423 refpeak=0.423 err/refpeak=6.830e-01 at omega=0.015 | c=2 fires=False
delta=0.20 ks=0.497 v='python' n= 23 x0=-1.0060 emax-err=-5.6e-16 ratio=1.000e+00 fid=5.7e-13 | max|y|=0.425 refpeak=0.425 err/refpeak=6.571e-01 at omega=0.020 | c=2 fires=False
delta=0.20 ks=0.497 v=3      n= 23 x0=-1.0060 emax-err=7.8e-05 ratio=1.000e+00 fid=2.8e-13 | max|y|=0.425 refpeak=0.425 err/refpeak=6.571e-01 at omega=0.020 | c=2 fires=False
delta=0.20 ks=0.497 v=2      n= 23 x0=-1.0060 emax-err=8.1e-04 ratio=1.000e+00 fid=1.3e-12 | max|y|=0.425 refpeak=0.425 err/refpeak=6.568e-01 at omega=0.020 | c=2 fires=False
delta=0.20 ks=0.495 v='python' n= 23 x0=-1.0101 emax-err=-5.6e-16 ratio=1.207e+00 fid=8.5e-13 | max|y|=0.427 refpeak=0.427 err/refpeak=6.045e-01 at omega=0.025 | c=2 fires=False
delta=0.20 ks=0.495 v=3      n= 23 x0=-1.0101 emax-err=8.1e-04 ratio=1.211e+00 fid=1.6e-12 | max|y|=0.427 refpeak=0.427 err/refpeak=6.042e-01 at omega=0.025 | c=2 fires=False
delta=0.20 ks=0.495 v=2      n= 23 x0=-1.0101 emax-err=1.5e-03 ratio=1.214e+00 fid=3.7e-13 | max|y|=0.427 refpeak=0.427 err/refpeak=6.038e-01 at omega=0.025 | c=2 fires=False
delta=0.20 ks=0.490 v='python' n= 22 x0=-1.0204 emax-err=-5.6e-16 ratio=1.632e+00 fid=1.9e-12 | max|y|=0.415 refpeak=0.415 err/refpeak=5.149e-01 at omega=0.060 | c=2 fires=False
delta=0.20 ks=0.490 v=3      n= 22 x0=-1.0204 emax-err=1.7e-05 ratio=1.632e+00 fid=3.5e-14 | max|y|=0.415 refpeak=0.415 err/refpeak=5.149e-01 at omega=0.060 | c=2 fires=False
delta=0.20 ks=0.490 v=2      n= 22 x0=-1.0204 emax-err=8.1e-05 ratio=1.632e+00 fid=1.4e-12 | max|y|=0.415 refpeak=0.415 err/refpeak=5.149e-01 at omega=0.060 | c=2 fires=False
delta=0.20 ks=0.485 v='python' n= 22 x0=-1.0309 emax-err=-5.6e-16 ratio=4.656e+00 fid=1.9e-12 | max|y|=0.498 refpeak=0.420 err/refpeak=1.188e+00 at omega=0.070 | c=2 fires=True
delta=0.20 ks=0.485 v=3      n= 22 x0=-1.0309 emax-err=3.8e-07 ratio=4.656e+00 fid=2.7e-12 | max|y|=0.498 refpeak=0.420 err/refpeak=1.188e+00 at omega=0.070 | c=2 fires=True
delta=0.20 ks=0.485 v=2      n= 22 x0=-1.0309 emax-err=7.8e-04 ratio=4.660e+00 fid=2.0e-13 | max|y|=0.498 refpeak=0.420 err/refpeak=1.188e+00 at omega=0.070 | c=2 fires=True
delta=0.20 ks=0.480 v='python' n= 22 x0=-1.0417 emax-err=-5.6e-16 ratio=1.152e+01 fid=1.7e-12 | max|y|=1.003 refpeak=0.426 err/refpeak=2.358e+00 at omega=0.075 | c=2 fires=True
delta=0.20 ks=0.480 v=3      n= 22 x0=-1.0417 emax-err=8.2e-06 ratio=1.152e+01 fid=1.4e-13 | max|y|=1.003 refpeak=0.426 err/refpeak=2.358e+00 at omega=0.075 | c=2 fires=True
delta=0.20 ks=0.480 v=2      n= 22 x0=-1.0417 emax-err=1.5e-03 ratio=1.152e+01 fid=1.3e-12 | max|y|=1.003 refpeak=0.426 err/refpeak=2.358e+00 at omega=0.075 | c=2 fires=True
delta=0.20 ks=0.470 v='python' n= 22 x0=-1.0638 emax-err=-5.6e-16 ratio=4.992e+01 fid=1.7e-12 | max|y|=3.106 refpeak=0.437 err/refpeak=7.105e+00 at omega=0.095 | c=2 fires=True
delta=0.20 ks=0.470 v=3      n= 22 x0=-1.0638 emax-err=6.2e-04 ratio=4.992e+01 fid=7.9e-13 | max|y|=3.110 refpeak=0.437 err/refpeak=7.113e+00 at omega=0.095 | c=2 fires=True
delta=0.20 ks=0.470 v=2      n= 22 x0=-1.0638 emax-err=6.7e-06 ratio=4.992e+01 fid=1.6e-12 | max|y|=3.106 refpeak=0.437 err/refpeak=7.105e+00 at omega=0.095 | c=2 fires=True
delta=0.20 ks=0.460 v='python' n= 21 x0=-1.0870 emax-err=-5.6e-16 ratio=1.142e+02 fid=1.7e-12 | max|y|=5.726 refpeak=0.430 err/refpeak=1.331e+01 at omega=0.120 | c=2 fires=True
delta=0.20 ks=0.460 v=3      n= 21 x0=-1.0870 emax-err=6.7e-05 ratio=1.142e+02 fid=3.9e-12 | max|y|=5.727 refpeak=0.430 err/refpeak=1.331e+01 at omega=0.120 | c=2 fires=True
delta=0.20 ks=0.460 v=2      n= 21 x0=-1.0870 emax-err=1.7e-03 ratio=1.142e+02 fid=1.3e-12 | max|y|=5.750 refpeak=0.431 err/refpeak=1.336e+01 at omega=0.120 | c=2 fires=True
--- delta=0.10  Lorentzian exact peak 0.4850
delta=0.10 ks=0.700 v='python' n= 64 x0=-0.7143 emax-err=-5.6e-16 ratio=1.000e+00 fid=2.2e-13 | max|y|=0.768 refpeak=0.768 err/refpeak=2.247e-13 at omega=-0.005 | c=2 fires=False
delta=0.10 ks=0.700 v=3      n= 64 x0=-0.7143 emax-err=2.7e-05 ratio=1.000e+00 fid=4.1e-13 | max|y|=0.768 refpeak=0.768 err/refpeak=4.203e-13 at omega=-0.005 | c=2 fires=False
delta=0.10 ks=0.700 v=2      n= 64 x0=-0.7143 emax-err=7.7e-05 ratio=1.000e+00 fid=2.2e-13 | max|y|=0.768 refpeak=0.768 err/refpeak=1.782e-13 at omega=-0.005 | c=2 fires=False
delta=0.10 ks=0.520 v='python' n= 48 x0=-0.9615 emax-err=-5.6e-16 ratio=1.000e+00 fid=2.5e-13 | max|y|=0.819 refpeak=0.819 err/refpeak=5.491e-13 at omega=-0.005 | c=2 fires=False
delta=0.10 ks=0.520 v=3      n= 48 x0=-0.9615 emax-err=1.5e-06 ratio=1.000e+00 fid=2.2e-13 | max|y|=0.819 refpeak=0.819 err/refpeak=5.383e-13 at omega=-0.005 | c=2 fires=False
delta=0.10 ks=0.520 v=2      n= 48 x0=-0.9615 emax-err=6.6e-05 ratio=1.000e+00 fid=4.3e-13 | max|y|=0.819 refpeak=0.819 err/refpeak=1.146e-12 at omega=-0.005 | c=2 fires=False
delta=0.10 ks=0.500 v='python' n= 46 x0=-1.0000 emax-err=-5.6e-16 ratio=1.000e+00 fid=2.1e-13 | max|y|=0.825 refpeak=0.825 err/refpeak=1.903e-13 at omega=0.015 | c=2 fires=False
delta=0.10 ks=0.500 v=3      n= 46 x0=-1.0000 emax-err=1.6e-04 ratio=1.000e+00 fid=3.2e-13 | max|y|=0.825 refpeak=0.825 err/refpeak=2.328e-13 at omega=0.675 | c=2 fires=False
delta=0.10 ks=0.500 v=2      n= 46 x0=-1.0000 emax-err=2.7e-05 ratio=1.000e+00 fid=9.5e-14 | max|y|=0.825 refpeak=0.825 err/refpeak=4.615e-14 at omega=0.660 | c=2 fires=False
delta=0.10 ks=0.499 v='python' n= 46 x0=-1.0020 emax-err=-5.6e-16 ratio=1.243e+00 fid=4.6e-13 | max|y|=0.828 refpeak=0.828 err/refpeak=1.956e-01 at omega=0.015 | c=2 fires=False
delta=0.10 ks=0.499 v=3      n= 46 x0=-1.0020 emax-err=9.5e-06 ratio=1.243e+00 fid=1.4e-12 | max|y|=0.828 refpeak=0.828 err/refpeak=1.956e-01 at omega=0.015 | c=2 fires=False
delta=0.10 ks=0.499 v=2      n= 46 x0=-1.0020 emax-err=7.7e-04 ratio=1.247e+00 fid=8.9e-13 | max|y|=0.828 refpeak=0.828 err/refpeak=1.956e-01 at omega=0.015 | c=2 fires=False
delta=0.10 ks=0.497 v='python' n= 46 x0=-1.0060 emax-err=-5.6e-16 ratio=4.770e+00 fid=1.3e-12 | max|y|=0.873 refpeak=0.831 err/refpeak=1.050e+00 at omega=0.020 | c=2 fires=True
delta=0.10 ks=0.497 v=3      n= 46 x0=-1.0060 emax-err=1.7e-05 ratio=4.770e+00 fid=1.5e-12 | max|y|=0.873 refpeak=0.831 err/refpeak=1.050e+00 at omega=0.020 | c=2 fires=True
delta=0.10 ks=0.497 v=2      n= 46 x0=-1.0060 emax-err=6.3e-05 ratio=4.771e+00 fid=3.7e-12 | max|y|=0.873 refpeak=0.831 err/refpeak=1.050e+00 at omega=0.020 | c=2 fires=True
delta=0.10 ks=0.495 v='python' n= 45 x0=-1.0101 emax-err=-5.6e-16 ratio=1.473e+01 fid=1.5e-12 | max|y|=1.931 refpeak=0.818 err/refpeak=2.362e+00 at omega=0.025 | c=2 fires=True
delta=0.10 ks=0.495 v=3      n= 45 x0=-1.0101 emax-err=4.4e-04 ratio=1.474e+01 fid=9.7e-13 | max|y|=1.930 refpeak=0.818 err/refpeak=2.360e+00 at omega=0.025 | c=2 fires=True
delta=0.10 ks=0.495 v=2      n= 45 x0=-1.0101 emax-err=4.1e-04 ratio=1.474e+01 fid=1.9e-13 | max|y|=1.930 refpeak=0.818 err/refpeak=2.360e+00 at omega=0.025 | c=2 fires=True
delta=0.10 ks=0.490 v='python' n= 45 x0=-1.0204 emax-err=-5.6e-16 ratio=2.073e+02 fid=1.2e-12 | max|y|=5.672 refpeak=0.828 err/refpeak=6.847e+00 at omega=0.050 | c=2 fires=True
delta=0.10 ks=0.490 v=3      n= 45 x0=-1.0204 emax-err=2.5e-05 ratio=2.073e+02 fid=1.2e-12 | max|y|=5.672 refpeak=0.828 err/refpeak=6.847e+00 at omega=0.050 | c=2 fires=True
delta=0.10 ks=0.490 v=2      n= 45 x0=-1.0204 emax-err=6.7e-04 ratio=2.073e+02 fid=1.2e-12 | max|y|=5.669 refpeak=0.828 err/refpeak=6.842e+00 at omega=0.050 | c=2 fires=True
delta=0.10 ks=0.485 v='python' n= 44 x0=-1.0309 emax-err=-5.6e-16 ratio=1.246e+03 fid=1.6e-12 | max|y|=29.474 refpeak=0.821 err/refpeak=3.590e+01 at omega=0.060 | c=2 fires=True
delta=0.10 ks=0.485 v=3      n= 44 x0=-1.0309 emax-err=5.9e-06 ratio=1.246e+03 fid=2.5e-12 | max|y|=29.474 refpeak=0.821 err/refpeak=3.590e+01 at omega=0.060 | c=2 fires=True
delta=0.10 ks=0.485 v=2      n= 44 x0=-1.0309 emax-err=1.5e-03 ratio=1.246e+03 fid=1.2e-12 | max|y|=29.447 refpeak=0.822 err/refpeak=3.584e+01 at omega=0.060 | c=2 fires=True
delta=0.10 ks=0.480 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.10 ks=0.480 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.10 ks=0.480 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.10 ks=0.470 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.10 ks=0.470 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.10 ks=0.470 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.10 ks=0.460 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.10 ks=0.460 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.10 ks=0.460 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
--- delta=0.05  Lorentzian exact peak 0.9592
delta=0.05 ks=0.700 v='python' n=128 x0=-0.7143 emax-err=-5.6e-16 ratio=1.000e+00 fid=2.5e-13 | max|y|=1.522 refpeak=1.522 err/refpeak=2.252e-13 at omega=0.000 | c=2 fires=False
delta=0.05 ks=0.700 v=3      n=128 x0=-0.7143 emax-err=3.2e-03 ratio=1.000e+00 fid=2.7e-13 | max|y|=1.524 refpeak=1.524 err/refpeak=2.355e-13 at omega=0.000 | c=2 fires=False
delta=0.05 ks=0.700 v=2      n=128 x0=-0.7143 emax-err=6.5e-03 ratio=1.000e+00 fid=1.4e-13 | max|y|=1.525 refpeak=1.525 err/refpeak=1.154e-13 at omega=1.365 | c=2 fires=False
delta=0.05 ks=0.520 v='python' n= 95 x0=-0.9615 emax-err=-5.6e-16 ratio=1.000e+00 fid=2.5e-13 | max|y|=1.603 refpeak=1.603 err/refpeak=5.255e-13 at omega=0.000 | c=2 fires=False
delta=0.05 ks=0.520 v=3      n= 95 x0=-0.9615 emax-err=1.4e-06 ratio=1.000e+00 fid=4.8e-13 | max|y|=1.603 refpeak=1.603 err/refpeak=1.076e-12 at omega=0.000 | c=2 fires=False
delta=0.05 ks=0.520 v=2      n= 95 x0=-0.9615 emax-err=8.3e-05 ratio=1.000e+00 fid=7.2e-14 | max|y|=1.603 refpeak=1.603 err/refpeak=7.383e-14 at omega=0.000 | c=2 fires=False
delta=0.05 ks=0.500 v='python' n= 92 x0=-1.0000 emax-err=-5.6e-16 ratio=1.000e+00 fid=4.6e-13 | max|y|=1.630 refpeak=1.630 err/refpeak=1.003e-13 at omega=1.370 | c=2 fires=False
delta=0.05 ks=0.500 v=3      n= 92 x0=-1.0000 emax-err=2.3e-05 ratio=1.000e+00 fid=7.4e-13 | max|y|=1.630 refpeak=1.630 err/refpeak=2.612e-13 at omega=1.370 | c=2 fires=False
delta=0.05 ks=0.500 v=2      n= 92 x0=-1.0000 emax-err=1.1e-03 ratio=1.000e+00 fid=1.1e-13 | max|y|=1.630 refpeak=1.630 err/refpeak=8.193e-14 at omega=1.365 | c=2 fires=False
delta=0.05 ks=0.499 v='python' n= 91 x0=-1.0020 emax-err=-5.6e-16 ratio=8.992e+00 fid=2.9e-13 | max|y|=1.616 refpeak=1.616 err/refpeak=1.557e-01 at omega=0.015 | c=2 fires=True
delta=0.05 ks=0.499 v=3      n= 91 x0=-1.0020 emax-err=7.0e-05 ratio=8.993e+00 fid=2.4e-12 | max|y|=1.617 refpeak=1.616 err/refpeak=1.557e-01 at omega=0.015 | c=2 fires=True
delta=0.05 ks=0.499 v=2      n= 91 x0=-1.0020 emax-err=2.3e-04 ratio=8.995e+00 fid=2.7e-12 | max|y|=1.617 refpeak=1.617 err/refpeak=1.557e-01 at omega=0.015 | c=2 fires=True
delta=0.05 ks=0.497 v='python' n= 91 x0=-1.0060 emax-err=-5.6e-16 ratio=5.701e+02 fid=1.1e-12 | max|y|=15.198 refpeak=1.625 err/refpeak=9.355e+00 at omega=0.020 | c=2 fires=True
delta=0.05 ks=0.497 v=3      n= 90 x0=-1.0060 emax-err=2.0e-02 ratio=5.104e+02 fid=4.0e-12 | max|y|=13.658 refpeak=1.616 err/refpeak=8.451e+00 at omega=0.020 | c=2 fires=True
delta=0.05 ks=0.497 v=2      n= 91 x0=-1.0060 emax-err=2.2e-05 ratio=5.701e+02 fid=3.8e-12 | max|y|=15.198 refpeak=1.625 err/refpeak=9.355e+00 at omega=0.020 | c=2 fires=True
delta=0.05 ks=0.495 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.495 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.495 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.490 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.490 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.490 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.485 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.485 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.485 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.480 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.480 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.480 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.470 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.470 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.470 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.460 v='python' RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.460 v=3      RAISES RuntimeError: KPM moments diverging: scaled spec
delta=0.05 ks=0.460 v=2      RAISES RuntimeError: KPM moments diverging: scaled spec
=== part 2: pair (1e-2*Sz_0, 1e-2*Sz_0), delta=0.1
ks=0.480 v='python' bound=2.500e-05 threshold=1000.0 max|mu|=1.712e-01 ratio=6.849e+03 | integral y/eps^2=1.1994e+00 (exact 0.25) err/refpeak=1.435e+02
ks=0.480 v=3      bound=2.500e-05 threshold=1000.0 max|mu|=1.712e-01 ratio=6.849e+03 | integral y/eps^2=1.1994e+00 (exact 0.25) err/refpeak=1.435e+02
ks=0.470 v='python' bound=2.500e-05 threshold=1000.0 max|mu|=2.209e+00 ratio=8.837e+04 | integral y/eps^2=7.7973e+00 (exact 0.25) err/refpeak=1.211e+03
ks=0.470 v=3      bound=2.500e-05 threshold=1000.0 max|mu|=2.209e+00 ratio=8.837e+04 | integral y/eps^2=7.7237e+00 (exact 0.25) err/refpeak=1.191e+03
ks=0.460 v='python' bound=2.500e-05 threshold=1000.0 max|mu|=1.714e+01 ratio=6.858e+05 | integral y/eps^2=4.0882e+01 (exact 0.25) err/refpeak=6.745e+03
ks=0.460 v=3      bound=2.500e-05 threshold=1000.0 max|mu|=1.714e+01 ratio=6.858e+05 | integral y/eps^2=4.0879e+01 (exact 0.25) err/refpeak=6.744e+03
ks=0.450 v='python' bound=2.500e-05 threshold=1000.0 max|mu|=9.471e+01 ratio=3.788e+06 | integral y/eps^2=1.7208e+02 (exact 0.25) err/refpeak=3.193e+04
ks=0.450 v=3      bound=2.500e-05 threshold=1000.0 max|mu|=9.471e+01 ratio=3.788e+06 | integral y/eps^2=1.7208e+02 (exact 0.25) err/refpeak=3.193e+04
```

**Reviewer (CONFIRMED, NARROWED)**: the band's size in the table above, and its
edge per chain: at delta=0.2 on the same chain silent down to 0.42, where the curve
is 109 times the peak (raising at 0.40); at delta=0.05 silent at 0.498 and 0.497;
8 sites silent at 0.498 and 0.497 (1.8 to 11 times); a 6-site spinless-fermion
(N_0,N_0) pair, elastic fraction 0.56, silent from 0.499 to 0.497 at 2.7 to 61
times; the fieldless chain, whose out-of-window weight sits at the top, 10.2 times
at 0.35 and at 0.30 (the user guide's truncation snippet's value, without
truncation); and operators scaled by 1e-2 give 172*eps^2 against 0.25 at 0.45,
about 700 times off, with no raise (`reviews/kpm_A2/02_band_true_anchor.py`,
`03_what_the_band_looks_like.py`). Struck: "fires only from about 0.48 down", a
single chain's number; and the natural worry that MPS truncation needs a looser
tolerance on DMRG, which measurement removes: 80 correct 12-site runs at `kpmmaxm`
8, 16 and 50, `kpm_scale` 0.7 and 0.55, delta 0.1 and 0.05, `kpm_accelerate` on
and off, auto and cross pairs, on v3 and `"python"`, errors up to 30 per cent from
truncation, and not one moment above ratio 1.000000 (truncation only loses norm);
healthy energy-truncated runs peak at 0.967 (`"python"`) and 0.977 (v3):

`reviews/kpm_A2/04_correct_run_ratio.py`:

```python
"""kpm_A2 probe 04: how large does max_k |mu_k|/(||vi|| ||vj||) get on CORRECT
DMRG runs, where MPS truncation actually bites? This is what sets how tight a
moment bound can be on the DMRG backends.

12-site S=1/2 chains, ground state at maxm=40/nsweeps=20 (pinned, seeded):
  H12f  Heisenberg + 0.2 Sz_0, pairs (Sz_0,Sz_0) and cross (Sz_0,Sz_5)
  H12s  Heisenberg + staggered 0.5*(-1)^i Sz_i, pair (Sz_5,Sz_5)
kpm_scale 0.7 and 0.55, delta 0.1 (and 0.05 on a subset), kpmmaxm 8/16/50,
kpm_accelerate on/off (the cross pair always runs the full recursion). The
bound is 0.25 exactly for every pair here (Sz^2 = 1/4 on S=1/2), so the ratio
is max|mu|/0.25. 'tail' is the max over the second half of the moments, where
truncation drift would show. Correctness is measured, not assumed: err is
max|y_DMRG - y_ED|/max(y_ED) against mode="ED" at the same kpm_scale/delta."""
import numpy as np
import a2lib as L

BOUND = 0.25


def build(v, kind):
    if kind == "H12f":
        return L.heis(v, ns=12, field=0.2, maxm=40, nsweeps=20)
    return L.heis(v, ns=12, field=0.0, stagger=0.5, maxm=40, nsweeps=20)


cases = [("H12f", 0, 0), ("H12f", 0, 5), ("H12s", 5, 5)]
ed_cache = {}
worst = {}
for kind, i, j in cases:
    for delta, kss, mms in ((0.1, (0.7, 0.55), (8, 16, 50)), (0.05, (0.55,), (8, 50))):
        for ks in kss:
            key = (kind, i, j, delta, ks)
            sc = build("python", kind)
            sc.kpm_scale = ks
            es = np.linspace(-0.5, 6.0, 1301)
            xe, ye = sc.get_dynamical_correlator(mode="ED", submode="KPM",
                        name=(sc.Sz[i], sc.Sz[j]), delta=delta, es=es)
            ye = np.real(np.asarray(ye))
            for v in ("python", 3):
                for mm in mms:
                    for acc in ((True, False) if i == j else (False,)):
                        sc = build(v, kind)
                        sc.kpm_scale = ks
                        sc.kpmmaxm = mm
                        sc.kpm_accelerate = acc
                        tag = "%s (%d,%d) delta=%.2f ks=%.2f v=%-6r kpmmaxm=%2d acc=%d" % (
                            kind, i, j, delta, ks, v, mm, acc)
                        try:
                            mus, emin, emax, scale, n, d = sc.get_dynamical_correlator_moments(
                                name=(sc.Sz[i], sc.Sz[j]), delta=delta)
                        except Exception as ex:
                            print(tag, "RAISES %s: %s" % (type(ex).__name__, str(ex)[:40]), flush=True)
                            continue
                        mus = np.asarray(mus)
                        y = np.real(L.reconstruct(mus, emin, emax, scale, n, es, delta))
                        r = np.abs(mus)/BOUND
                        err = np.max(np.abs(y-ye))/np.max(np.abs(ye))
                        print("%s n=%3d ratio max=%.6f tail=%.6f last=%.6f | err vs ED=%.3e"
                              % (tag, n, r.max(), r[n//2:].max(), r[-1], err), flush=True)
                        worst[(v, acc)] = max(worst.get((v, acc), 0.0), r.max())
print("worst ratio per (backend, accelerate):", {str(k): round(x, 6) for k, x in worst.items()})
```

`reviews/kpm_A2/04_correct_run_ratio.out`:

```
H12f (0,0) delta=0.10 ks=0.70 v='python' kpmmaxm= 8 acc=1 n=207 ratio max=0.999985 tail=0.552108 last=0.438282 | err vs ED=1.358e-01
H12f (0,0) delta=0.10 ks=0.70 v='python' kpmmaxm= 8 acc=0 n=207 ratio max=0.999985 tail=0.505689 last=0.382250 | err vs ED=8.803e-02
H12f (0,0) delta=0.10 ks=0.70 v='python' kpmmaxm=16 acc=1 n=207 ratio max=1.000000 tail=0.498587 last=0.427380 | err vs ED=1.244e-03
H12f (0,0) delta=0.10 ks=0.70 v='python' kpmmaxm=16 acc=0 n=207 ratio max=1.000000 tail=0.499615 last=0.425567 | err vs ED=1.448e-03
H12f (0,0) delta=0.10 ks=0.70 v='python' kpmmaxm=50 acc=1 n=207 ratio max=1.000000 tail=0.500002 last=0.427158 | err vs ED=6.853e-04
H12f (0,0) delta=0.10 ks=0.70 v='python' kpmmaxm=50 acc=0 n=207 ratio max=1.000000 tail=0.500002 last=0.427158 | err vs ED=6.853e-04
H12f (0,0) delta=0.10 ks=0.70 v=3      kpmmaxm= 8 acc=1 n=207 ratio max=0.999985 tail=0.559369 last=0.419295 | err vs ED=1.315e-01
H12f (0,0) delta=0.10 ks=0.70 v=3      kpmmaxm= 8 acc=0 n=207 ratio max=0.999985 tail=0.511205 last=0.367375 | err vs ED=8.426e-02
H12f (0,0) delta=0.10 ks=0.70 v=3      kpmmaxm=16 acc=1 n=207 ratio max=1.000000 tail=0.503889 last=0.416585 | err vs ED=1.180e-03
H12f (0,0) delta=0.10 ks=0.70 v=3      kpmmaxm=16 acc=0 n=207 ratio max=1.000000 tail=0.500166 last=0.425209 | err vs ED=1.307e-03
H12f (0,0) delta=0.10 ks=0.70 v=3      kpmmaxm=50 acc=1 n=207 ratio max=1.000000 tail=0.500479 last=0.426236 | err vs ED=6.156e-04
H12f (0,0) delta=0.10 ks=0.70 v=3      kpmmaxm=50 acc=0 n=207 ratio max=1.000000 tail=0.500513 last=0.426170 | err vs ED=6.105e-04
H12f (0,0) delta=0.10 ks=0.55 v='python' kpmmaxm= 8 acc=1 n=163 ratio max=0.999985 tail=0.645719 last=0.256549 | err vs ED=7.480e-02
H12f (0,0) delta=0.10 ks=0.55 v='python' kpmmaxm= 8 acc=0 n=163 ratio max=0.999985 tail=0.655056 last=0.324549 | err vs ED=5.938e-02
H12f (0,0) delta=0.10 ks=0.55 v='python' kpmmaxm=16 acc=1 n=163 ratio max=1.000000 tail=0.632362 last=0.263875 | err vs ED=1.610e-03
H12f (0,0) delta=0.10 ks=0.55 v='python' kpmmaxm=16 acc=0 n=163 ratio max=1.000000 tail=0.631540 last=0.264365 | err vs ED=1.589e-03
H12f (0,0) delta=0.10 ks=0.55 v='python' kpmmaxm=50 acc=1 n=163 ratio max=1.000000 tail=0.632458 last=0.262381 | err vs ED=1.618e-03
H12f (0,0) delta=0.10 ks=0.55 v='python' kpmmaxm=50 acc=0 n=163 ratio max=1.000000 tail=0.632458 last=0.262381 | err vs ED=1.618e-03
H12f (0,0) delta=0.10 ks=0.55 v=3      kpmmaxm= 8 acc=1 n=163 ratio max=0.999985 tail=0.657148 last=0.242104 | err vs ED=7.545e-02
H12f (0,0) delta=0.10 ks=0.55 v=3      kpmmaxm= 8 acc=0 n=163 ratio max=0.999985 tail=0.664473 last=0.317347 | err vs ED=6.361e-02
H12f (0,0) delta=0.10 ks=0.55 v=3      kpmmaxm=16 acc=1 n=163 ratio max=1.000000 tail=0.643998 last=0.280934 | err vs ED=1.285e-03
H12f (0,0) delta=0.10 ks=0.55 v=3      kpmmaxm=16 acc=0 n=163 ratio max=1.000000 tail=0.634733 last=0.271033 | err vs ED=1.243e-03
H12f (0,0) delta=0.10 ks=0.55 v=3      kpmmaxm=50 acc=1 n=163 ratio max=1.000000 tail=0.630400 last=0.265413 | err vs ED=1.427e-03
H12f (0,0) delta=0.10 ks=0.55 v=3      kpmmaxm=50 acc=0 n=163 ratio max=1.000000 tail=0.632139 last=0.262861 | err vs ED=1.588e-03
H12f (0,0) delta=0.05 ks=0.55 v='python' kpmmaxm= 8 acc=1 n=326 ratio max=0.999985 tail=0.675140 last=0.085558 | err vs ED=1.369e-01
H12f (0,0) delta=0.05 ks=0.55 v='python' kpmmaxm= 8 acc=0 n=326 ratio max=0.999985 tail=0.698599 last=0.017662 | err vs ED=1.050e-01
H12f (0,0) delta=0.05 ks=0.55 v='python' kpmmaxm=50 acc=1 n=326 ratio max=1.000000 tail=0.695174 last=0.027421 | err vs ED=1.451e-03
H12f (0,0) delta=0.05 ks=0.55 v='python' kpmmaxm=50 acc=0 n=326 ratio max=1.000000 tail=0.695174 last=0.027421 | err vs ED=1.451e-03
H12f (0,0) delta=0.05 ks=0.55 v=3      kpmmaxm= 8 acc=1 n=326 ratio max=0.999985 tail=0.686667 last=0.076698 | err vs ED=1.404e-01
H12f (0,0) delta=0.05 ks=0.55 v=3      kpmmaxm= 8 acc=0 n=326 ratio max=0.999985 tail=0.717927 last=0.006099 | err vs ED=1.172e-01
H12f (0,0) delta=0.05 ks=0.55 v=3      kpmmaxm=50 acc=1 n=326 ratio max=1.000000 tail=0.697479 last=0.034291 | err vs ED=1.196e-03
H12f (0,0) delta=0.05 ks=0.55 v=3      kpmmaxm=50 acc=0 n=326 ratio max=1.000000 tail=0.698967 last=0.039391 | err vs ED=9.967e-04
H12f (0,5) delta=0.10 ks=0.70 v='python' kpmmaxm= 8 acc=0 n=207 ratio max=0.442294 tail=0.407339 last=0.008252 | err vs ED=3.074e-01
H12f (0,5) delta=0.10 ks=0.70 v='python' kpmmaxm=16 acc=0 n=207 ratio max=0.438497 tail=0.319759 last=0.139258 | err vs ED=1.539e-02
H12f (0,5) delta=0.10 ks=0.70 v='python' kpmmaxm=50 acc=0 n=207 ratio max=0.428338 tail=0.314513 last=0.136843 | err vs ED=8.037e-04
H12f (0,5) delta=0.10 ks=0.70 v=3      kpmmaxm= 8 acc=0 n=207 ratio max=0.458846 tail=0.407240 last=0.014149 | err vs ED=3.252e-01
H12f (0,5) delta=0.10 ks=0.70 v=3      kpmmaxm=16 acc=0 n=207 ratio max=0.437837 tail=0.320070 last=0.139788 | err vs ED=1.437e-02
H12f (0,5) delta=0.10 ks=0.70 v=3      kpmmaxm=50 acc=0 n=207 ratio max=0.429471 tail=0.315008 last=0.138417 | err vs ED=3.872e-04
H12f (0,5) delta=0.10 ks=0.55 v='python' kpmmaxm= 8 acc=0 n=163 ratio max=0.495813 tail=0.495813 last=0.118785 | err vs ED=2.342e-01
H12f (0,5) delta=0.10 ks=0.55 v='python' kpmmaxm=16 acc=0 n=163 ratio max=0.487054 tail=0.487054 last=0.121626 | err vs ED=1.188e-02
H12f (0,5) delta=0.10 ks=0.55 v='python' kpmmaxm=50 acc=0 n=163 ratio max=0.485898 tail=0.485898 last=0.129430 | err vs ED=1.324e-03
H12f (0,5) delta=0.10 ks=0.55 v=3      kpmmaxm= 8 acc=0 n=163 ratio max=0.449369 tail=0.444745 last=0.049109 | err vs ED=2.319e-01
H12f (0,5) delta=0.10 ks=0.55 v=3      kpmmaxm=16 acc=0 n=163 ratio max=0.483338 tail=0.483338 last=0.124258 | err vs ED=1.221e-02
H12f (0,5) delta=0.10 ks=0.55 v=3      kpmmaxm=50 acc=0 n=163 ratio max=0.483508 tail=0.483508 last=0.134302 | err vs ED=8.024e-04
H12f (0,5) delta=0.05 ks=0.55 v='python' kpmmaxm= 8 acc=0 n=326 ratio max=0.495813 tail=0.427105 last=0.305565 | err vs ED=2.839e-01
H12f (0,5) delta=0.05 ks=0.55 v='python' kpmmaxm=50 acc=0 n=326 ratio max=0.485898 tail=0.441917 last=0.033234 | err vs ED=1.573e-03
H12f (0,5) delta=0.05 ks=0.55 v=3      kpmmaxm= 8 acc=0 n=326 ratio max=0.449939 tail=0.325121 last=0.286280 | err vs ED=3.017e-01
H12f (0,5) delta=0.05 ks=0.55 v=3      kpmmaxm=50 acc=0 n=326 ratio max=0.485213 tail=0.442854 last=0.037694 | err vs ED=1.413e-03
H12s (5,5) delta=0.10 ks=0.70 v='python' kpmmaxm= 8 acc=1 n=268 ratio max=1.000000 tail=0.807962 last=0.599134 | err vs ED=8.913e-03
H12s (5,5) delta=0.10 ks=0.70 v='python' kpmmaxm= 8 acc=0 n=268 ratio max=1.000000 tail=0.823798 last=0.582372 | err vs ED=9.983e-03
H12s (5,5) delta=0.10 ks=0.70 v='python' kpmmaxm=16 acc=1 n=268 ratio max=1.000000 tail=0.830062 last=0.542524 | err vs ED=8.442e-04
H12s (5,5) delta=0.10 ks=0.70 v='python' kpmmaxm=16 acc=0 n=268 ratio max=1.000000 tail=0.830128 last=0.542381 | err vs ED=8.459e-04
H12s (5,5) delta=0.10 ks=0.70 v='python' kpmmaxm=50 acc=1 n=268 ratio max=1.000000 tail=0.829779 last=0.542011 | err vs ED=8.524e-04
H12s (5,5) delta=0.10 ks=0.70 v='python' kpmmaxm=50 acc=0 n=268 ratio max=1.000000 tail=0.829779 last=0.542011 | err vs ED=8.524e-04
H12s (5,5) delta=0.10 ks=0.70 v=3      kpmmaxm= 8 acc=1 n=267 ratio max=1.000000 tail=0.809086 last=0.109407 | err vs ED=8.732e-03
H12s (5,5) delta=0.10 ks=0.70 v=3      kpmmaxm= 8 acc=0 n=268 ratio max=1.000000 tail=0.812319 last=0.578298 | err vs ED=9.310e-03
H12s (5,5) delta=0.10 ks=0.70 v=3      kpmmaxm=16 acc=1 n=268 ratio max=1.000000 tail=0.830515 last=0.547780 | err vs ED=4.655e-04
H12s (5,5) delta=0.10 ks=0.70 v=3      kpmmaxm=16 acc=0 n=268 ratio max=1.000000 tail=0.830366 last=0.545226 | err vs ED=6.267e-04
H12s (5,5) delta=0.10 ks=0.70 v=3      kpmmaxm=50 acc=1 n=268 ratio max=1.000000 tail=0.829788 last=0.542112 | err vs ED=8.449e-04
H12s (5,5) delta=0.10 ks=0.70 v=3      kpmmaxm=50 acc=0 n=268 ratio max=1.000000 tail=0.829785 last=0.542086 | err vs ED=8.469e-04
H12s (5,5) delta=0.10 ks=0.55 v='python' kpmmaxm= 8 acc=1 n=210 ratio max=1.000000 tail=0.800704 last=0.035193 | err vs ED=6.507e-03
H12s (5,5) delta=0.10 ks=0.55 v='python' kpmmaxm= 8 acc=0 n=210 ratio max=1.000000 tail=0.814442 last=0.025968 | err vs ED=5.867e-03
H12s (5,5) delta=0.10 ks=0.55 v='python' kpmmaxm=16 acc=1 n=210 ratio max=1.000000 tail=0.810118 last=0.071963 | err vs ED=2.043e-03
H12s (5,5) delta=0.10 ks=0.55 v='python' kpmmaxm=16 acc=0 n=210 ratio max=1.000000 tail=0.810064 last=0.071991 | err vs ED=2.042e-03
H12s (5,5) delta=0.10 ks=0.55 v='python' kpmmaxm=50 acc=1 n=210 ratio max=1.000000 tail=0.809939 last=0.072489 | err vs ED=2.042e-03
H12s (5,5) delta=0.10 ks=0.55 v='python' kpmmaxm=50 acc=0 n=210 ratio max=1.000000 tail=0.809939 last=0.072489 | err vs ED=2.042e-03
H12s (5,5) delta=0.10 ks=0.55 v=3      kpmmaxm= 8 acc=1 n=210 ratio max=1.000000 tail=0.812218 last=0.037117 | err vs ED=6.399e-03
H12s (5,5) delta=0.10 ks=0.55 v=3      kpmmaxm= 8 acc=0 n=210 ratio max=1.000000 tail=0.819215 last=0.023770 | err vs ED=5.565e-03
H12s (5,5) delta=0.10 ks=0.55 v=3      kpmmaxm=16 acc=1 n=210 ratio max=1.000000 tail=0.809920 last=0.071937 | err vs ED=2.038e-03
H12s (5,5) delta=0.10 ks=0.55 v=3      kpmmaxm=16 acc=0 n=210 ratio max=1.000000 tail=0.810071 last=0.071870 | err vs ED=2.025e-03
H12s (5,5) delta=0.10 ks=0.55 v=3      kpmmaxm=50 acc=1 n=210 ratio max=1.000000 tail=0.809930 last=0.072478 | err vs ED=2.041e-03
H12s (5,5) delta=0.10 ks=0.55 v=3      kpmmaxm=50 acc=0 n=210 ratio max=1.000000 tail=0.809751 last=0.072238 | err vs ED=2.019e-03
H12s (5,5) delta=0.05 ks=0.55 v='python' kpmmaxm= 8 acc=1 n=421 ratio max=1.000000 tail=0.865135 last=0.219986 | err vs ED=1.394e-02
H12s (5,5) delta=0.05 ks=0.55 v='python' kpmmaxm= 8 acc=0 n=421 ratio max=1.000000 tail=0.834907 last=0.184954 | err vs ED=8.961e-03
H12s (5,5) delta=0.05 ks=0.55 v='python' kpmmaxm=50 acc=1 n=421 ratio max=1.000000 tail=0.854481 last=0.175032 | err vs ED=1.582e-03
H12s (5,5) delta=0.05 ks=0.55 v='python' kpmmaxm=50 acc=0 n=421 ratio max=1.000000 tail=0.854481 last=0.175032 | err vs ED=1.582e-03
H12s (5,5) delta=0.05 ks=0.55 v=3      kpmmaxm= 8 acc=1 n=421 ratio max=1.000000 tail=0.878623 last=0.267887 | err vs ED=1.233e-02
H12s (5,5) delta=0.05 ks=0.55 v=3      kpmmaxm= 8 acc=0 n=421 ratio max=1.000000 tail=0.837908 last=0.192388 | err vs ED=9.716e-03
H12s (5,5) delta=0.05 ks=0.55 v=3      kpmmaxm=50 acc=1 n=420 ratio max=1.000000 tail=0.846686 last=0.372539 | err vs ED=4.212e-03
H12s (5,5) delta=0.05 ks=0.55 v=3      kpmmaxm=50 acc=0 n=421 ratio max=1.000000 tail=0.854463 last=0.175046 | err vs ED=1.581e-03
worst ratio per (backend, accelerate): {"('python', True)": np.float64(1.0), "('python', False)": np.float64(1.0), '(3, True)': np.float64(1.0), '(3, False)': np.float64(1.0)}
```

Also narrowed: no bound closes the sliver just below 1/2, where the elastic line
is already mangled at a ratio of 1.2 to 1.6 (at 0.499, delta=0.1, a -0.162 dip and
a low-energy integral of -0.0008 against the elastic weight +0.0145, the
inelastic part right to 9e-5).

**Suggested fix**: replace `1e3*(bound+1.0)` with `c*bound`, c=1.5, in all four
guards, including the truncated v3 loops (a rebuild of both extensions), and
check both moments appended per accelerated step. A trial of c=2 monkeypatched
into `"python"` (`reviews/kpm_A2/05_tight_guard_trial.py`) had zero false
positives on 17 correct runs, including the exact fieldless answers below 1/2 and
the test suite's own truncated setting (ratio 0.96), and raised on every band row
with an error of 2 per cent or more; the one behaviour change is that the
fieldless chain at 0.45, delta=0.025 now raises where the shipped guard returned a
2 per cent wrong spectrum. `general_kpm` shares these loops, where c*bound holds
for a Hermitian X only; `nhkpm_moments` does not call the guard. Beyond the
bound: on the midpoint window x0 does not depend on the noisy emax, so raise when
the weight that will grow, |<vj|gs><gs|vi>|, exceeds a small fraction of the bound
at `kpm_scale < 0.5`, pointing at `kpm_energy_truncate`. Qualify the user guide's
promise in both formats.

### 4. `mode="ED"` KPM never reads `kpm_energy_truncate`, so under the flag ED silently computes the centred-window spectrum while v3 and `"python"` switch to the ground-state-anchored window, which makes two O2 sentences false and, through `mode.py`'s fallback, turns the user guide's truncation example on a 2-site v3 chain into finding 2's divergence

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `kpm` &middot; origin `51357d8` (the flag), made a stated contradiction by `765b537`

**Status**: FIXED. `edtk/dynamics.py::get_dynamical_correlator` computes `herm = is_hermitian(h)` once and, on `submode=="KPM" and herm`, after the T>0 guard and before `get_gs_array()`, raises `NotImplementedError` under `kpm_energy_truncate`, naming `mode.py`'s fallback in the message; `Many_Body_Chain.get_dynamical_correlator` pushes the flag across. The non-Hermitian KPM reads neither the flag nor `kpm_scale` on either mode and still returns. The exact projected anchored window on ED is left as a feature for later. Pinned by `tests/test_audit_2026_09_24b_kpm.py::test_ed_kpm_refuses_energy_truncation_before_the_ground_state` (a spy on `get_gs_array`), `::test_two_site_v3_fallback_refuses_energy_truncation_by_name` and `::test_non_hermitian_ed_kpm_ignores_both_flags`. No number changes: the ED curve under the flag was identical to the flag-off one and now raises. The known-issue sentence and the O2 sentences in `dynamics.py` and `CLAUDE.md` are corrected.

**Where**: `src/dmrgpy/manybodychain.py:985-986` (the ED push of
`get_dynamical_correlator` forwards `kpm_scale` and `kpm_n_scale` only) and
`src/dmrgpy/edtk/dynamics.py::dynamical_correlator_kpm` (always the
bandwidth-centred window).

With the flag set the ED curve is bit for bit the flag-off curve. Nothing
documents the boundary: the user guide, README and ROADMAP scope the flag to v3
and `"python"` and document the v2 raise, and `kpmdmrg.py:142-151` already turned
exactly this accepted-and-ignored shape into a `NotImplementedError` on v2. Under
the flag two O2 statements are false: `dynamics.py`'s "they now share both" and
`CLAUDE.md`'s "the ED route adopted the DMRG rescaling so the two share the same x
at the same physical energy".

**Expected**: either the anchored window on ED, or a raise naming the
unsupported combination.

Repro, from the hunter, part (a) of the script (parts (b) and (c) belong to
finding 6 and to "Ruled out"):

```bash
cd <scratch>/kpm && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 06_truncate_flag_on_ed_and_remaining_routes.py
```

`kpm/06_truncate_flag_on_ed_and_remaining_routes.py`:

```python
"""Probe 06.
(a) kpm_energy_truncate=True at a SAFE kpm_scale: does mode="ED" answer the
    same calculation as mode="DMRG"? Anchor: exact Chebyshev moments from
    dense numpy on each of the two windows (bandwidth-centred and
    ground-state-anchored, both with the calibrated n), reconstructed with
    a from-scratch Jackson kernel.
(b) kpm_finite: window_chain_kwargs with a misspelled key.
(c) the remaining kpm_n_scale readers: Kondo second-order DMRG, fermionchain
    get_gr, kpm_finite via window_chain_kwargs."""
import numpy as np
np.random.seed(3)
from dmrgpy import spinchain, fermionchain, infinitechain
from dmrgpy.algebra.kpm import polynomials_for_broadening


def build(v="python"):
    np.random.seed(3)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=v)
    h = 0
    for i in range(3):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 30
    return sc


def jackson_curve(sc, A, B, es, delta, window, kpm_scale):
    ed = sc.get_ED_obj()
    H = np.array(ed.get_hamiltonian().todense())
    e, V = np.linalg.eigh(H); gs = V[:, 0]; E0 = e[0]; W = e[-1]-e[0]
    if window == "centred":
        half = W*kpm_scale; centre = E0 + W/2.
    else:  # ground-state anchored, pyitensor/chain.py::_scaled_hamiltonian_gs_anchored
        wp = 1.0-0.025/2.; half = W*kpm_scale/(2*wp); centre = E0 + wp*half
    n = polynomials_for_broadening(half, delta)
    Hs = (H - centre*np.eye(len(H)))/half
    Am = np.array(ed.MO2matrix(A).todense()); Bm = np.array(ed.MO2matrix(B).todense())
    vb = Bm@gs; va = np.conjugate(Am.T)@gs
    t0, t1 = vb.copy(), Hs@vb; mus = [np.vdot(va, t0), np.vdot(va, t1)]
    for k in range(2, n):
        t0, t1 = t1, 2*Hs@t1 - t0; mus.append(np.vdot(va, t1))
    mus = np.array(mus); N = n; k = np.arange(N)
    g = ((N-k+1)*np.cos(np.pi*k/(N+1)) + np.sin(np.pi*k/(N+1))/np.tan(np.pi/(N+1)))/(N+1)
    xs = (es + E0 - centre)/half
    ok = np.abs(xs) < 0.99
    y = np.zeros(len(es))
    T = np.cos(np.outer(np.arccos(xs[ok]), k)); w = np.where(k == 0, 1.0, 2.0)
    y[ok] = np.real(T@(g*w*mus))/(np.pi*np.sqrt(1-xs[ok]**2))/half
    return y, n


print("(a) kpm_energy_truncate=True, safe kpm_scale")
es = np.linspace(0.0, 3.0, 601)
delta = 0.1
for ks in (0.7, 1.2):
    ref_c, n_c = jackson_curve(build(), build().Sz[0], build().Sz[0], es, delta, "centred", ks)
    ref_a, n_a = jackson_curve(build(), build().Sz[0], build().Sz[0], es, delta, "anchored", ks)
    out = {}
    for mode, v in (("ED", "python"), ("DMRG", "python"), ("DMRG", 3)):
        sc = build(v); sc.kpm_scale = ks; sc.kpm_energy_truncate = True
        x, y = sc.get_dynamical_correlator(mode=mode, submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
                                           delta=delta, es=es)
        y = np.real(np.asarray(y)); out[(mode, v)] = y
        print("  kpm_scale=%.1f %-4s v=%-7r  max|y-anchored ref (n=%d)|=%.2e  max|y-centred ref (n=%d)|=%.2e  peak=%.4f"
              % (ks, mode, v, n_a, np.max(np.abs(y-ref_a)), n_c, np.max(np.abs(y-ref_c)), np.max(y)))
    d = np.max(np.abs(out[("ED", "python")]-out[("DMRG", "python")]))
    print("  kpm_scale=%.1f  max|ED - DMRG(python)| = %.3e = %.1f%% of the DMRG peak %.4f"
          % (ks, d, 100*d/np.max(out[("DMRG", "python")]), np.max(out[("DMRG", "python")])))

print("(b) kpm_finite, window_chain_kwargs")
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.set_hamiltonian(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0])
esw = np.linspace(-0.5, 3.0, 351)
curves = {}
for label, wk in (("default", dict()),
                  ("kpm_n_scale=3", dict(kpm_n_scale=3)),
                  ("kpm_nscale=3 (typo)", dict(kpm_nscale=3)),
                  ("kpmscale=0.9 (typo)", dict(kpmscale=0.9))):
    np.random.seed(3)
    wk = dict(maxm=20, nsweeps=10, **wk)
    try:
        x, y = ic.kpm_finite("Sz", 0, "Sz", 0, n_window=6, window_chain_kwargs=wk, delta=0.2, es=esw)
        curves[label] = np.real(np.asarray(y))
        print("  %-22s returns, peak=%.5f" % (label, np.max(curves[label])))
    except Exception as ex:
        print("  %-22s RAISES %s: %s" % (label, type(ex).__name__, str(ex)[:80]))
for label in ("kpm_n_scale=3", "kpm_nscale=3 (typo)", "kpmscale=0.9 (typo)"):
    if label in curves:
        print("  max|%s - default| = %.2e" % (label, np.max(np.abs(curves[label]-curves["default"]))))
for val in (1.5, 0):
    try:
        ic.kpm_finite("Sz", 0, "Sz", 0, n_window=6,
                      window_chain_kwargs=dict(maxm=20, nsweeps=10, kpm_n_scale=val), delta=0.2, es=esw)
        print("  window kpm_n_scale=%r returns" % (val,))
    except Exception as ex:
        print("  window kpm_n_scale=%r RAISES %s (names it: %s)" % (val, type(ex).__name__, "kpm_n_scale" in str(ex)))

print("(c) remaining kpm_n_scale readers at kpm_n_scale=1.5")
sc = spinchain.Spin_Chain(["S=1"]*3, itensor_version="python")
h = sc.Sx[0]*sc.Sx[1]+sc.Sy[0]*sc.Sy[1]+sc.Sz[0]*sc.Sz[1] + 0.3*sc.Sz[0]*sc.Sz[0]
h = h + sc.Sx[1]*sc.Sx[2]+sc.Sy[1]*sc.Sy[2]+sc.Sz[1]*sc.Sz[2]
sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = 20, 10
sc.kpm_n_scale = 1.5
try:
    sc.get_kondo_spectrum(np.linspace(-1, 1, 5), site=0, T=0.0, order=2, mode="DMRG",
                          delta=0.05, es=np.linspace(-0.3, 3.0, 200))
    print("  Kondo order=2 mode=DMRG returns")
except Exception as ex:
    print("  Kondo order=2 mode=DMRG RAISES %s (names it: %s), gs_done=%s"
          % (type(ex).__name__, "kpm_n_scale" in str(ex), sc.computed_gs))
fc = fermionchain.Fermionic_Chain(3, itensor_version="python")
fc.set_hamiltonian(fc.Cdag[0]*fc.C[1]+fc.Cdag[1]*fc.C[0]+fc.Cdag[1]*fc.C[2]+fc.Cdag[2]*fc.C[1])
fc.maxm, fc.nsweeps = 20, 10
fc.kpm_n_scale = 1.5
try:
    fc.get_gr(delta=0.1, es=np.linspace(-2, 2, 41))
    print("  fermionchain get_gr returns")
except Exception as ex:
    print("  fermionchain get_gr RAISES %s (names it: %s), gs_done=%s"
          % (type(ex).__name__, "kpm_n_scale" in str(ex), fc.computed_gs))
```

`kpm/06_truncate_flag_on_ed_and_remaining_routes.out`:

```
(a) kpm_energy_truncate=True, safe kpm_scale
  kpm_scale=0.7 ED   v='python'  max|y-anchored ref (n=32)|=2.59e+08  max|y-centred ref (n=64)|=1.57e-05  peak=0.7682
  kpm_scale=0.7 DMRG v='python'  max|y-anchored ref (n=32)|=2.59e+08  max|y-centred ref (n=64)|=1.70e+00  peak=0.7247
  kpm_scale=0.7 DMRG v=3        max|y-anchored ref (n=32)|=2.59e+08  max|y-centred ref (n=64)|=1.70e+00  peak=0.7247
  kpm_scale=0.7  max|ED - DMRG(python)| = 1.700e+00 = 234.5% of the DMRG peak 0.7247
  kpm_scale=1.2 ED   v='python'  max|y-anchored ref (n=56)|=3.73e-01  max|y-centred ref (n=110)|=1.58e-14  peak=0.7372
  kpm_scale=1.2 DMRG v='python'  max|y-anchored ref (n=56)|=5.07e-03  max|y-centred ref (n=110)|=3.68e-01  peak=0.8705
  kpm_scale=1.2 DMRG v=3        max|y-anchored ref (n=56)|=4.38e-03  max|y-centred ref (n=110)|=3.69e-01  peak=0.8715
  kpm_scale=1.2  max|ED - DMRG(python)| = 3.683e-01 = 42.3% of the DMRG peak 0.8705
(b) kpm_finite, window_chain_kwargs
  default                returns, peak=0.30669
  kpm_n_scale=3          returns, peak=0.90702
  kpm_nscale=3 (typo)    returns, peak=0.30669
  kpmscale=0.9 (typo)    returns, peak=0.30669
  max|kpm_n_scale=3 - default| = 6.00e-01
  max|kpm_nscale=3 (typo) - default| = 0.00e+00
  max|kpmscale=0.9 (typo) - default| = 0.00e+00
  window kpm_n_scale=1.5 RAISES TypeError (names it: True)
  window kpm_n_scale=0 RAISES ValueError (names it: True)
(c) remaining kpm_n_scale readers at kpm_n_scale=1.5
  Kondo order=2 mode=DMRG RAISES TypeError (names it: True), gs_done=False
  fermionchain get_gr RAISES TypeError (names it: True), gs_done=False
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced exactly (0.368, 42.3 per cent of
0.8705, on `"python"` and v3), confirmed the hunter's anchored reference is exactly
`pyitensor/chain.py::_scaled_hamiltonian_gs_anchored` and its C++ twin, and that
truncation is inactive at `kpm_scale=1.2` (3.2e-14 against
`kpm_truncate_nsweeps=0`). Struck: the 42 per cent as a defect size. At 1.2 both
curves are exact reconstructions on their own windows and agree on every pole
position and weight; the gap is the documented sqrt(1-x^2) resolution profile
read on two windows (the elastic pole sits at x=-0.9875 in the anchored window,
height ratio 6.16 against a predicted 5.77; the inelastic peak ratio 1.180 against
1.164). Struck also: the hunter's `kpm_scale=0.7` row as evidence, since there the
truncated DMRG curve dips to -1.457 just below the window top and ED is the sane
one (a known-issue lead, below). Severity MEDIUM to LOW: the wrong number appears
only where this meets finding 2, which the record must not count twice.

`reviews/kpm_B/03_what_the_42_percent_is.py`:

```python
"""Review of kpm finding B: what the 42 per cent ED-vs-DMRG gap under the flag
is made of. Same 4-site chain, kpm_scale=1.2, delta=0.1.
- split the gap into the elastic pole at omega=0 (|<GS|Sz0|GS>|^2, nonzero
  because 0.2*Sz_0 breaks SU(2)) and the inelastic part;
- compare the height ratio at each pole with the ratio of sqrt(1-x^2) of the
  two windows, the documented Jackson resolution profile;
- locate the DMRG-vs-anchored-reference residual (5e-3);
- account for the DMRG sum-rule deficit under the flag: how much of the
  elastic line lies below x=-0.99, the reconstruction grid's lower edge."""
import numpy as np
from math import erf
from dmrgpy import spinchain


def build(v="python"):
    np.random.seed(3)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=v)
    h = 0
    for i in range(3):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 30
    return sc


ks, delta = 1.2, 0.1
sc = build()
ed = sc.get_ED_obj()
H = np.array(ed.get_hamiltonian().todense()); e, V = np.linalg.eigh(H)
E0, W = e[0], e[-1]-e[0]
Am = np.array(ed.MO2matrix(sc.Sz[0]).todense())
M = np.abs(V.conj().T@(Am@V[:, 0]))**2
print("elastic weight |<GS|Sz0|GS>|^2 = %.5f, <GS|Sz0|GS> = %.5f" % (M[0], np.real(V[:, 0].conj()@Am@V[:, 0])))

es = np.linspace(-0.2, 3.0, 3201)
curves = {}
for mode, v in (("ED", "python"), ("DMRG", "python")):
    s = build(v); s.kpm_scale = ks; s.kpm_energy_truncate = True
    x, y = s.get_dynamical_correlator(mode=mode, submode="KPM", name=(s.Sz[0], s.Sz[0]), delta=delta, es=es)
    curves[mode] = np.real(np.asarray(y))
d = np.abs(curves["ED"]-curves["DMRG"])
pk = curves["DMRG"].max()
for lo, hi in ((-0.2, 0.2), (0.2, 3.0)):
    m = (es >= lo) & (es < hi)
    i = np.argmax(d*m)
    print("  omega in [%.1f,%.1f): max|ED-DMRG| = %.4f (%.1f%% of the DMRG peak %.4f) at omega=%.3f"
          % (lo, hi, d[i], 100*d[i]/pk, pk, es[i]))

wp = 1-0.025/2
wins = {"centred": (W*ks, E0+W/2), "anchored": (W*ks/(2*wp), E0+wp*W*ks/(2*wp))}
for i in (0, int(np.argmax(M[1:]))+1):
    om = e[i]-E0
    f = {}
    for k, (half, c) in wins.items():
        xx = (e[i]-c)/half; f[k] = np.sqrt(1-xx*xx)
    j = np.argmin(np.abs(es-om))
    # local maximum near the pole
    sl = slice(max(j-40, 0), j+40)
    hE, hD = curves["ED"][sl].max(), curves["DMRG"][sl].max()
    print("  pole omega=%.4f weight=%.4f: height DMRG/ED = %.3f, profile ratio sqrt(1-x_c^2)/sqrt(1-x_a^2) = %.3f/%.3f = %.3f"
          % (om, M[i], hD/hE, f["centred"], f["anchored"], f["centred"]/f["anchored"]))

# DMRG-vs-anchored-reference residual location: reuse script 01's reference
import importlib.util, sys, os
spec = importlib.util.spec_from_file_location("r01", os.path.join(os.path.dirname(os.path.abspath(__file__)), "01_flag_on_ed_vs_dmrg.py"))
src = open(spec.origin).read().split("delta = 0.1\nks = 1.2")[0]
ns = {}; exec(src, ns)
ref_a, n_a = ns["jackson_curve"](build(), es, delta, "anchored", ks)
r = np.abs(curves["DMRG"]-ref_a)
i = np.argmax(r)
print("  max|DMRG - anchored ref| = %.2e at omega=%.4f; restricted to omega>0.2: %.2e"
      % (r[i], es[i], np.max(r[es > 0.2])))

# elastic line below the reconstruction grid's lower edge x=-0.99
half = wins["anchored"][0]; xa = (E0-wins["anchored"][1])/half
N = n_a
# Jackson line of a delta at x: sigma ~ pi*sqrt(1-x^2)/N in x (Weisse et al.)
sig = np.pi*np.sqrt(1-xa*xa)/N
frac_below = 0.5*(1+erf((-0.99-xa)/(sig*np.sqrt(2))))
print("  anchored: E0 at x=%.4f, grid edge -0.99 is %.4f below it, elastic line sigma_x=%.4f -> fraction below edge %.3f"
      % (xa, xa+0.99, sig, frac_below))
print("  predicted weight lost = %.5f; measured integral deficit on es in [-0.2,3]: ED %.5f, DMRG %.5f"
      % (frac_below*M[0], M.sum()-np.trapezoid(curves["ED"], es), M.sum()-np.trapezoid(curves["DMRG"], es)))
```

`reviews/kpm_B/03_what_the_42_percent_is.out`:

```
elastic weight |<GS|Sz0|GS>|^2 = 0.01453, <GS|Sz0|GS> = -0.12054
  omega in [-0.2,0.2): max|ED-DMRG| = 0.3974 (45.6% of the DMRG peak 0.8706) at omega=-0.003
  omega in [0.2,3.0): max|ED-DMRG| = 0.1340 (15.4% of the DMRG peak 0.8706) at omega=0.672
  pole omega=0.0000 weight=0.0145: height DMRG/ED = 6.158, profile ratio sqrt(1-x_c^2)/sqrt(1-x_a^2) = 0.909/0.158 = 5.767
  pole omega=0.6779 weight=0.1504: height DMRG/ED = 1.180, profile ratio sqrt(1-x_c^2)/sqrt(1-x_a^2) = 0.982/0.843 = 1.164
  max|DMRG - anchored ref| = 6.44e-03 at omega=-0.0010; restricted to omega>0.2: 5.87e-04
  anchored: E0 at x=-0.9875, grid edge -0.99 is 0.0025 below it, elastic line sigma_x=0.0088 -> fraction below edge 0.389
  predicted weight lost = 0.00565; measured integral deficit on es in [-0.2,3]: ED 0.00005, DMRG 0.00562
```

**Suggested fix**: raise `NotImplementedError` in the ED push when
`self.kpm_energy_truncate` is set and the submode is KPM, naming the fallback in
the message the way `get_distribution_moments`' ED branch does, since the caller
may have asked for DMRG. No test or example breaks (every truncation test uses a
separate unflagged ED chain on `submode="INV"`). Cost: a 2-site v3 fallback that
today returns a correct centred curve at 0.7 or 1.2 would refuse instead, the v2
precedent. Gate it on the Hermitian branch, since the non-Hermitian KPM reads
neither the flag nor `kpm_scale` on either mode. The exact projected anchored
window on ED is a feature for later, and the reviewer's prototype of it
(`reviews/kpm_B/04_kpm_scale_07_row.py`) is the right test anchor. Correct the
known-issue sentence and the two O2 sentences.

### 5. `30200a4` placed the `mode="ED"` `kpm_n_scale` check one frame ahead of the Hermiticity branch that decides whether the value is read, so on a non-Hermitian Hamiltonian ED rejects `kpm_n_scale=1.5` while DMRG accepts it, and neither reads it

`bug` &middot; severity **LOW** &middot; CONFIRMED &middot; lens `kpm` &middot; origin `30200a4`

**Status**: FIXED, with finding 4 in the same branch. The `kpm_n_scale` check sits in `edtk/dynamics.py`'s Hermitian KPM branch, after the T>0 guard and before the ground state, so a T>0 call gets its own error back; the pre-check in `manybodychain.py` is gone and its comment points there, and `validate_kpm_n_scale`'s docstring lists the new call site. Pinned by `tests/test_audit_2026_09_24b_kpm.py::test_non_hermitian_ed_kpm_ignores_both_flags` and `::test_hermitian_ed_kpm_n_scale_raises_before_the_ground_state` (1.5 and 0). No number changes: a non-Hermitian ED KPM at `kpm_n_scale=1.5` used to raise `TypeError` and now returns the spectrum at 1 exactly (difference 0.0).

**Where**: `src/dmrgpy/manybodychain.py` (the pre-check on
`kwargs.get("submode","KPM")=="KPM"` in the ED push, added by `30200a4`), one
frame before `edtk/dynamics.py`'s non-Hermitian branch into
`nonhermitian/kpm.py::dynamical_correlator_nhkpm_ed`; the DMRG check sits inside
`kpmdmrg.dynamical_correlator_moments`, after `dynamics.py`'s non-Hermitian branch
has returned.

The fix's own rule, stated four times (`validate_kpm_n_scale`'s docstring,
`test_kpm_n_scale_is_only_checked_where_it_is_read`, `documentation.md:4664`, the
previous record's finding 5 Status, "the same reach as on `mode="DMRG"`"), is that
the value is checked only where it is read. Non-Hermitian KPM takes `n=` and never
reads `kpm_n_scale` on either mode (`nonhermitian/kpm.py`, `nhkpm_moments` in
`pyitensor/chain.py` and `mpscpp3`), so ED's reach is wider than DMRG's. It is
the section 4.10 shape, a precondition ahead of the branch it qualifies. The
rejection also fires on the fallback routes (`sc.mode="ED"`, v3 on `ns<3`), so a
caller who asked for DMRG on a 2-site non-Hermitian chain gets the ED rejection.
It survived because every `kpm_n_scale` test uses a Hermitian chain.

**Expected**: the same answer on both modes; since nothing reads the value there,
both accept it.

Repro, from the hunter, part (b) of the script (part (a) is a null in "Ruled
out"):

```bash
cd <scratch>/kpm && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 03_distribution_nondiagonal_and_nh.py
```

`kpm/03_distribution_nondiagonal_and_nh.py`:

```python
"""Probe 03.
(a) get_distribution on a NON-diagonal X (the regression tests use Sz only),
    mode="ED" both branches (own grid, xs=) and mode="DMRG" on "python":
    total weight, max|Im|, and per-line weights against |<v|GS>|^2 from
    np.linalg.eigh.
(b) Non-Hermitian KPM: is kpm_n_scale read at all? Spectra at
    kpm_n_scale = 1 and 3, mode="ED" and mode="DMRG" ("python", v3)."""
import numpy as np
np.random.seed(11)
from dmrgpy import spinchain


def heis(v="python", nonherm=False):
    np.random.seed(11)
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=v)
    h = 0
    for i in range(3):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0] + 0.1*sc.Sx[3]
    if nonherm: h = h + 0.3j*sc.Sz[1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 30
    return sc


def exact_lines(sc, X):
    ed = sc.get_ED_obj()
    H = np.array(ed.get_hamiltonian().todense())
    _e, V = np.linalg.eigh(H); gs = V[:, 0]
    xm = np.array(ed.MO2matrix(X).todense())
    lam, u = np.linalg.eigh(xm)
    w = np.abs(np.conjugate(u.T)@gs)**2
    vals = np.unique(np.round(lam, 8))
    return vals, np.array([np.sum(w[np.abs(lam-v) < 1e-6]) for v in vals])


def line_weights(x, y, vals):
    x = np.asarray(x, dtype=float); y = np.real(np.asarray(y))
    mids = np.concatenate([[x[0]-1], (vals[1:]+vals[:-1])/2, [x[-1]+1]])
    return np.array([np.trapezoid(y[(x >= mids[k]) & (x < mids[k+1])],
                                  x[(x >= mids[k]) & (x < mids[k+1])]) for k in range(len(vals))])


print("(a) get_distribution, non-diagonal X")
sc = heis()
for label, X in (("Sx0+Sx1", sc.Sx[0]+sc.Sx[1]), ("Sz0*Sz1+0.3*Sx2", sc.Sz[0]*sc.Sz[1]+0.3*sc.Sx[2])):
    vals, wex = exact_lines(sc, X)
    for scale in (2.0, None):
        x, y = sc.get_distribution(mode="ED", X=X, delta=0.02, scale=scale)
        wk = line_weights(x, y, vals)
        print("  ED  grid %-16s scale=%-4s total=%.6f max|Im|=%.1e  max|w_line-exact|=%.1e  lines=%s exact=%s"
              % (label, scale, np.trapezoid(np.real(y), x), np.max(np.abs(np.imag(y))),
                 np.max(np.abs(wk-wex)), np.round(wk, 5), np.round(wex, 5)))
        xs = np.linspace(vals[0]-0.5, vals[-1]+0.5, 4001)
        x2, y2 = sc.get_distribution(mode="ED", X=X, delta=0.02, scale=scale, xs=xs)
        wk2 = line_weights(x2, y2, vals)
        print("  ED  xs=  %-16s scale=%-4s total=%.6f max|Im|=%.1e  max|w_line-exact|=%.1e"
              % (label, scale, np.trapezoid(np.real(y2), x2), np.max(np.abs(np.imag(y2))),
                 np.max(np.abs(wk2-wex))))
    x, y = sc.get_distribution(mode="DMRG", X=X, delta=0.02)
    wk = line_weights(x, y, vals)
    print("  DMRG python %-16s        total=%.6f max|Im|=%.1e  max|w_line-exact|=%.1e  lines=%s"
          % (label, np.trapezoid(np.real(y), x), np.max(np.abs(np.imag(y))),
             np.max(np.abs(wk-wex)), np.round(wk, 5)))

print("(b) non-Hermitian KPM: kpm_n_scale 1 vs 3")
es = np.linspace(0.0, 3.0, 13)
for mode, v in (("ED", "python"), ("DMRG", "python"), ("DMRG", 3)):
    ys = []
    for ns in (1, 3):
        sc = heis(v, nonherm=True)
        assert not sc.is_hermitian(sc.hamiltonian)
        sc.kpm_n_scale = ns
        _x, y = sc.get_dynamical_correlator(mode=mode, submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
                                            delta=0.2, es=es, E_max=10, n=60)
        ys.append(np.asarray(y))
    print("  NH-KPM mode=%-4s v=%-6r  max|y(kpm_n_scale=3)-y(kpm_n_scale=1)| = %.1e   (peak %.4f)"
          % (mode, v, np.max(np.abs(ys[1]-ys[0])), np.max(np.abs(ys[0]))))
```

`kpm/03_distribution_nondiagonal_and_nh.out`:

```
(a) get_distribution, non-diagonal X
  ED  grid Sx0+Sx1          scale=2.0  total=1.000000 max|Im|=0.0e+00  max|w_line-exact|=3.8e-06  lines=[0.0255  0.94214 0.03236] exact=[0.02549 0.94215 0.03236]
  ED  xs=  Sx0+Sx1          scale=2.0  total=1.000000 max|Im|=0.0e+00  max|w_line-exact|=3.8e-06
  ED  grid Sx0+Sx1          scale=None total=1.000000 max|Im|=0.0e+00  max|w_line-exact|=3.8e-06  lines=[0.0255  0.94214 0.03236] exact=[0.02549 0.94215 0.03236]
  ED  xs=  Sx0+Sx1          scale=None total=1.000000 max|Im|=0.0e+00  max|w_line-exact|=3.8e-06
  DMRG python Sx0+Sx1                 total=1.000000 max|Im|=1.4e-15  max|w_line-exact|=8.7e-06  lines=[0.0255  0.94214 0.03236]
  ED  grid Sz0*Sz1+0.3*Sx2  scale=2.0  total=1.000000 max|Im|=0.0e+00  max|w_line-exact|=1.3e-04  lines=[0.42545 0.52735 0.02187 0.02532] exact=[0.42544 0.52748 0.02176 0.02532]
  ED  xs=  Sz0*Sz1+0.3*Sx2  scale=2.0  total=0.999999 max|Im|=0.0e+00  max|w_line-exact|=1.3e-04
  ED  grid Sz0*Sz1+0.3*Sx2  scale=None total=0.999999 max|Im|=0.0e+00  max|w_line-exact|=1.4e-04  lines=[0.42546 0.52735 0.02187 0.02532] exact=[0.42544 0.52748 0.02176 0.02532]
  ED  xs=  Sz0*Sz1+0.3*Sx2  scale=None total=0.999999 max|Im|=0.0e+00  max|w_line-exact|=1.4e-04
  DMRG python Sz0*Sz1+0.3*Sx2         total=0.999986 max|Im|=8.8e-15  max|w_line-exact|=3.4e-04  lines=[0.42548 0.52715 0.02203 0.02533]
(b) non-Hermitian KPM: kpm_n_scale 1 vs 3
  NH-KPM mode=ED   v='python'  max|y(kpm_n_scale=3)-y(kpm_n_scale=1)| = 0.0e+00   (peak 57.3752)
  NH-KPM mode=DMRG v='python'  max|y(kpm_n_scale=3)-y(kpm_n_scale=1)| = 0.0e+00   (peak 57.3752)
  NH-KPM mode=DMRG v=3       max|y(kpm_n_scale=3)-y(kpm_n_scale=1)| = 9.3e-13   (peak 57.3752)
```

**Reviewer (CONFIRMED)**: nothing struck; extended to the fallback routes and
confirmed the v3 difference (9.3e-13) is run-to-run noise (2.6e-13 on a rerun at
the same value); on the parent `c65e39e`, extracted to scratch, `mode="ED"`
accepts every value and returns the value bit-identical to HEAD's at 1, so the
behaviour is new in `30200a4`:

`reviews/kpm_C/01_nh_kpm_validator_reach.py`:

```python
"""Review of kpm finding C.
(a) Which routes raise on an invalid kpm_n_scale for a NON-Hermitian H
    under submode="KPM" (explicit and default), mode="ED" against
    mode="DMRG" on "python" and v3, and the fallback-to-ED routes
    (sc.mode="ED"; v3 on 2 sites).
(b) Is kpm_n_scale read by NH-KPM at all: spectra at kpm_n_scale=1 and 3,
    and on v3 the run-to-run floor at the SAME kpm_n_scale, so a nonzero
    difference can be attributed.
(c) Hermitian control: ED and DMRG both raise on the Hermitian KPM route,
    and ED's other submodes accept (the pinned "where read" rule)."""
import numpy as np, warnings
warnings.filterwarnings("ignore")
from dmrgpy import spinchain, cppext
import dmrgpy
print("dmrgpy from", dmrgpy.__file__, " v3 compiled:", cppext.available(3))

VALUES = [1.5, 2.0, 0, -1, True, np.float64(2.0), np.True_, np.int64(2), 1]
ES = np.linspace(0.0, 3.0, 7)

def chain(n=4, v="python", nonherm=True, seed=11):
    np.random.seed(seed)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    h = h + 0.2*sc.Sz[0] + 0.1*sc.Sx[n-1]
    if nonherm: h = h + 0.3j*sc.Sz[1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 30
    return sc

def outcome(f):
    try:
        f(); return "ok"
    except Exception as e:
        return "%s%s" % (type(e).__name__, "" if "kpm_n_scale" in str(e) else "[" + str(e)[:60] + "]")

print("\n(a) non-Hermitian H, submode KPM, kpm_n_scale value -> outcome")
routes = [
  ("ED  explicit submode", lambda: chain(), dict(mode="ED", submode="KPM")),
  ("ED  default submode",  lambda: chain(), dict(mode="ED")),
  ("DMRG python",          lambda: chain(), dict(mode="DMRG", submode="KPM")),
  ("DMRG v3",              lambda: chain(v=3), dict(mode="DMRG", submode="KPM")),
  ("DMRG v3 2-site (->ED)",lambda: chain(n=2, v=3), dict(mode="DMRG", submode="KPM")),
]
def ed_chain():
    sc = chain(); sc.mode = "ED"; return sc
routes.append(("sc.mode=ED, mode=DMRG", ed_chain, dict(mode="DMRG", submode="KPM")))
for label, make, kw in routes:
    row = []
    for val in VALUES:
        sc = make(); assert not sc.is_hermitian(sc.hamiltonian)
        sc.kpm_n_scale = val
        row.append(outcome(lambda: sc.get_dynamical_correlator(
            name=(sc.Sz[0], sc.Sz[0]), delta=0.2, es=ES, E_max=10, n=60, **kw)))
    print("  %-24s " % label + "  ".join("%r:%s" % (v, r) for v, r in zip(VALUES, row)), flush=True)

print("\n(b) NH-KPM spectra, kpm_n_scale 1 vs 3, and same-value floor")
for mode, v in (("ED", "python"), ("DMRG", "python"), ("DMRG", 3)):
    ys = {}
    for tag, ns in (("1a", 1), ("3", 3), ("1b", 1)):
        sc = chain(v=v); sc.kpm_n_scale = ns
        _x, y = sc.get_dynamical_correlator(mode=mode, submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
                                            delta=0.2, es=ES, E_max=10, n=60)
        ys[tag] = np.asarray(y)
    print("  mode=%-4s v=%-8r |y(3)-y(1a)|=%.1e  |y(1b)-y(1a)|=%.1e  peak=%.4f"
          % (mode, v, np.max(np.abs(ys["3"]-ys["1a"])), np.max(np.abs(ys["1b"]-ys["1a"])),
             np.max(np.abs(ys["1a"]))), flush=True)
# the parameter the route does read: n
sc = chain(); _x, y60 = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
                                                    delta=0.2, es=ES, E_max=10, n=60)
_x, y120 = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
                                       delta=0.2, es=ES, E_max=10, n=120)
print("  control: ED NH-KPM n=60 vs n=120  max|diff|=%.3e (the route does read n=)" % np.max(np.abs(np.asarray(y120)-np.asarray(y60))))

print("\n(c) Hermitian H, kpm_n_scale=1.5")
for label, v, kw in (("ED KPM", "python", dict(mode="ED", submode="KPM")),
                     ("ED default", "python", dict(mode="ED")),
                     ("ED submode=ED", "python", dict(mode="ED", submode="ED")),
                     ("ED submode=INV", "python", dict(mode="ED", submode="INV")),
                     ("DMRG python KPM", "python", dict(mode="DMRG", submode="KPM")),
                     ("DMRG v3 KPM", 3, dict(mode="DMRG", submode="KPM"))):
    sc = chain(v=v, nonherm=False); sc.kpm_n_scale = 1.5
    print("  %-18s %s" % (label, outcome(lambda: sc.get_dynamical_correlator(
        name=(sc.Sz[0], sc.Sz[0]), delta=0.2, es=ES, **kw))), flush=True)
```

`reviews/kpm_C/01_nh_kpm_validator_reach.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py  v3 compiled: True

(a) non-Hermitian H, submode KPM, kpm_n_scale value -> outcome
  ED  explicit submode     1.5:TypeError  2.0:TypeError  0:ValueError  -1:ValueError  True:TypeError  np.float64(2.0):TypeError  np.True_:TypeError  np.int64(2):ok  1:ok
  ED  default submode      1.5:TypeError  2.0:TypeError  0:ValueError  -1:ValueError  True:TypeError  np.float64(2.0):TypeError  np.True_:TypeError  np.int64(2):ok  1:ok
  DMRG python              1.5:ok  2.0:ok  0:ok  -1:ok  True:ok  np.float64(2.0):ok  np.True_:ok  np.int64(2):ok  1:ok
  DMRG v3                  1.5:ok  2.0:ok  0:ok  -1:ok  True:ok  np.float64(2.0):ok  np.True_:ok  np.int64(2):ok  1:ok
  DMRG v3 2-site (->ED)    1.5:TypeError  2.0:TypeError  0:ValueError  -1:ValueError  True:TypeError  np.float64(2.0):TypeError  np.True_:TypeError  np.int64(2):ok  1:ok
  sc.mode=ED, mode=DMRG    1.5:TypeError  2.0:TypeError  0:ValueError  -1:ValueError  True:TypeError  np.float64(2.0):TypeError  np.True_:TypeError  np.int64(2):ok  1:ok

(b) NH-KPM spectra, kpm_n_scale 1 vs 3, and same-value floor
  mode=ED   v='python' |y(3)-y(1a)|=0.0e+00  |y(1b)-y(1a)|=0.0e+00  peak=55.8234
  mode=DMRG v='python' |y(3)-y(1a)|=0.0e+00  |y(1b)-y(1a)|=0.0e+00  peak=55.8234
  mode=DMRG v=3        |y(3)-y(1a)|=3.4e-13  |y(1b)-y(1a)|=2.6e-13  peak=55.8234
  control: ED NH-KPM n=60 vs n=120  max|diff|=6.295e+01 (the route does read n=)

(c) Hermitian H, kpm_n_scale=1.5
  ED KPM             TypeError
  ED default         TypeError
  ED submode=ED      ok
  ED submode=INV     ok
  DMRG python KPM    TypeError
  DMRG v3 KPM        TypeError
```

A cosmetic side effect: on a Hermitian chain at T=0.1 with the default submode
the pre-check's `TypeError` pre-empts the `NotImplementedError` that says T>0
needs `submode='ED'`, the error the caller can act on.

**Suggested fix**: move the check into `edtk/dynamics.py`, after the T>0 guard
and before `get_gs_array()`, gated on the branch's own `is_hermitian(h)` computed
once into a local; delete the pre-check in `manybodychain.py`. Dropping the
pre-check outright loses "before any ground-state work" (Hermitian ED KPM then
raises only after the ground state and the band-edge solve), and gating it in
`manybodychain.py` on `self.is_hermitian` uses a different predicate from the
branch. Tried in a scratch copy of the tree (`reviews/kpm_C/fix_move/`), it
passes `tests/test_audit_2026_09_24_kpm.py` (50 passed, as does the drop variant,
so the existing tests pin what the check does, not where). Add a non-Hermitian
mirror of `test_kpm_n_scale_is_only_checked_where_it_is_read` and a spy on
`get_gs_array` pinning the Hermitian raise before the ground state; say in the
non-Hermitian KPM entry that its count is `n=`. No number changes.

### 6. `kpm_finite`'s `window_chain_kwargs` is a bare `setattr` loop on a temporary chain, so a misspelled key is stored where nothing reads it and the call returns the default-settings spectrum bit for bit: `kpm_nscale=3` gives the default where `kpm_n_scale=3` moves the peak by 196 per cent

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `kpm` &middot; origin `d22a576`

**Status**: FIXED. `kpm_finite` rejects `itensor_version` and `mode` by name, with a message pointing at the hardcoded `"python"` window, and then raises one `TypeError` naming, sorted, every key the temporary chain does not hold as a setting: a missing attribute, a private name or a method (a method passes `hasattr` too, and `gs_energy=3` would have overwritten it); both checks run before any DMRG. `deconvolve` is gone from the docstring. Pinned by `tests/test_audit_2026_09_24b_misc.py::test_misspelled_window_key_raises` (`kpm_nscale`, `kpmscale`, `max_m`, `nsweep`), `::test_every_unknown_window_key_is_named_sorted`, `::test_backend_keys_are_rejected_by_name`, `::test_method_or_private_window_key_raises` and `::test_correct_window_key_still_moves_the_spectrum`. No number changes: the correct `kpm_n_scale=3` gives the peak 0.90702 before and after on the reviewer's Heisenberg chain, and every in-tree caller passes only `maxm`, `nsweeps` or `kpmmaxm`.

**Where**: `src/dmrgpy/infinitechain.py:1368-1369`, `for k, v in
(window_chain_kwargs or {}).items(): setattr(wc, k, v)`, on a `Many_Body_Chain`
the method builds and discards.

The docstring defines the legal set as "the same attributes an ordinary finite
chain exposes (manybodychain.py's Many_Body_Chain.__init__)", so an unknown name
is outside the contract, not an alias; and since the chain carrying the stray
attribute is gone when the call returns, the caller has nothing on which to
notice it. The same method raises `TypeError` on `kpm_n_scale=3` as a call
keyword, while it drops `kpm_nscale=3` in the dict silently. `30200a4` touched
only the value path, which is why a correctly spelled `kpm_n_scale=1.5` in the
dict now raises. Evidence: part (b) of finding 4's script above, and the
reviewer's sweep over four intended/typo pairs and three broadenings:

`reviews/kpm_D/01_window_kwargs_typos.py`:

```python
"""Review of kpm lens candidate D: kpm_finite's window_chain_kwargs is a bare
setattr loop, so a misspelled key is accepted and ignored.
(1) the hunter's part (b), then the same at delta=0.1 and 0.35, with more
    typo/intended pairs (kpm_n_scale, kpm_scale, maxm);
(2) seed control: is bit-identity of typo vs default meaningful?"""
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.set_hamiltonian(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0])
esw = np.linspace(-0.5, 3.0, 351)

def run(wk, delta, seed=3):
    np.random.seed(seed)
    x, y = ic.kpm_finite("Sz", 0, "Sz", 0, n_window=6, window_chain_kwargs=wk,
                         delta=delta, es=esw)
    return np.real(np.asarray(y))

base = dict(maxm=20, nsweeps=10)
pairs = [  # (intended, typo)
    ("kpm_n_scale=3", dict(kpm_n_scale=3), "kpm_nscale=3", dict(kpm_nscale=3)),
    ("kpm_scale=0.9", dict(kpm_scale=0.9), "kpmscale=0.9", dict(kpmscale=0.9)),
    ("maxm=2",        dict(maxm=2),        "max_m=2",      dict(max_m=2)),
    ("nsweeps=1",     dict(nsweeps=1),     "nsweep=1",     dict(nsweep=1)),
]
for delta in (0.2, 0.1, 0.35):
    print("delta=%.2f" % delta)
    y0 = run(dict(base), delta)
    y0b = run(dict(base), delta)
    y0s = run(dict(base), delta, seed=11)
    print("  default peak=%.5f   same-seed rerun max|d|=%.2e   seed=11 rerun max|d|=%.2e"
          % (np.max(y0), np.max(np.abs(y0b-y0)), np.max(np.abs(y0s-y0))))
    for li, wi, lt, wt in pairs:
        yi = run(dict(base, **wi), delta)
        yt = run(dict(base, **wt), delta)
        print("  intended %-14s peak=%.5f max|y-default|=%.3e (%.1f%% of default peak) | "
              "typo %-13s peak=%.5f max|y-default|=%.2e bit-identical=%s"
              % (li, np.max(yi), np.max(np.abs(yi-y0)), 100*np.max(np.abs(yi-y0))/np.max(y0),
                 lt, np.max(yt), np.max(np.abs(yt-y0)), np.array_equal(yt, y0)))
```

`reviews/kpm_D/01_window_kwargs_typos.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
delta=0.20
  default peak=0.30669   same-seed rerun max|d|=0.00e+00   seed=11 rerun max|d|=1.38e-14
  intended kpm_n_scale=3  peak=0.90702 max|y-default|=6.003e-01 (195.7% of default peak) | typo kpm_nscale=3  peak=0.30669 max|y-default|=0.00e+00 bit-identical=True
  intended kpm_scale=0.9  peak=0.28579 max|y-default|=2.111e-02 (6.9% of default peak) | typo kpmscale=0.9  peak=0.30669 max|y-default|=0.00e+00 bit-identical=True
  intended maxm=2         peak=0.29783 max|y-default|=1.910e-01 (62.3% of default peak) | typo max_m=2       peak=0.30669 max|y-default|=0.00e+00 bit-identical=True
  intended nsweeps=1      peak=0.30669 max|y-default|=3.963e-09 (0.0% of default peak) | typo nsweep=1      peak=0.30669 max|y-default|=0.00e+00 bit-identical=True
delta=0.10
  default peak=0.61319   same-seed rerun max|d|=0.00e+00   seed=11 rerun max|d|=3.16e-14
  intended kpm_n_scale=3  peak=1.82490 max|y-default|=1.212e+00 (197.6% of default peak) | typo kpm_nscale=3  peak=0.61319 max|y-default|=0.00e+00 bit-identical=True
  intended kpm_scale=0.9  peak=0.57156 max|y-default|=4.163e-02 (6.8% of default peak) | typo kpmscale=0.9  peak=0.61319 max|y-default|=0.00e+00 bit-identical=True
  intended maxm=2         peak=0.58352 max|y-default|=5.784e-01 (94.3% of default peak) | typo max_m=2       peak=0.61319 max|y-default|=0.00e+00 bit-identical=True
  intended nsweeps=1      peak=0.61319 max|y-default|=7.917e-09 (0.0% of default peak) | typo nsweep=1      peak=0.61319 max|y-default|=0.00e+00 bit-identical=True
delta=0.35
  default peak=0.18250   same-seed rerun max|d|=0.00e+00   seed=11 rerun max|d|=7.80e-15
  intended kpm_n_scale=3  peak=0.53186 max|y-default|=3.497e-01 (191.6% of default peak) | typo kpm_nscale=3  peak=0.18250 max|y-default|=0.00e+00 bit-identical=True
  intended kpm_scale=0.9  peak=0.16870 max|y-default|=1.404e-02 (7.7% of default peak) | typo kpmscale=0.9  peak=0.18250 max|y-default|=0.00e+00 bit-identical=True
  intended maxm=2         peak=0.17574 max|y-default|=7.310e-02 (40.1% of default peak) | typo max_m=2       peak=0.18250 max|y-default|=0.00e+00 bit-identical=True
  intended nsweeps=1      peak=0.18250 max|y-default|=2.364e-09 (0.0% of default peak) | typo nsweep=1      peak=0.18250 max|y-default|=0.00e+00 bit-identical=True
```

**Reviewer (CONFIRMED, NARROWED)**: not recorded before (the 2026-09 reviewer read
`kpm_finite` and cleared only its hardcoded `"python"` window); the only
user-facing `setattr` loop in `src/`. Struck: the 0.60 as the defect's size (it is
the knob's effect, 192 to 198 per cent of the default peak for `kpm_n_scale=3`,
40 to 94 per cent for `maxm=2`, about 7 per cent for `kpm_scale=0.9`); and
"elsewhere the codebase raises on unknown call keywords" as the comparison, since
`Many_Body_Chain.__init__(**kwargs)` itself drops every keyword
(`Spin_Chain(['1/2']*4, maxm=50).maxm` is 30), which the 2026-08 record noted at
its line 1723 as a reviewer sharpening that never became a finding (a lead below).
A correctly spelled `itensor_version=3` in the dict is accepted and ignored too: it
passes `hasattr`, is set after the session is built, and the window runs on
`"python"` bit-identical to the default while reporting 3.

**Suggested fix**: raise `TypeError` naming any key for which `hasattr(wc, k)` is
false, before the `setattr`, and reject `itensor_version` and `mode` by name with a
message pointing at the hardcoded `"python"` window. All 27 documented or consumed
knobs pass the `hasattr` test (every one is created by `__init__`). No number
changes for correct keys. Also remove `deconvolve` from the docstring's list of
forwarded keywords, since `deconvolve=True` now raises.

### 7. `sxt_to_skomega` is hard-wired to the FFT-plus-interpolation stage that finding 7 of the previous record replaced, so every infinite-chain `td_dynamical_correlator` returns S(k,omega) as a straight line between FFT nodes: at a converged window, delta*T = 6, 20 per cent of the peak off with heights down to 0.79 of exact, where the direct sum on the same series is at 0.25 per cent

`bug` &middot; severity **MEDIUM** (latent for every checked-in caller, real for any caller who converges the window) &middot; CONFIRMED, NARROWED &middot; lens `realtime` &middot; origin `30200a4` (the hard-wiring and its docstring)

**Status**: FIXED, together with finding 8. `sxt_to_skomega` no longer passes `_evaluation="fft"`, so every infinite-chain S(k,omega) goes through the direct per-frequency sum; the switch stays, and `_fourier_transform_correlator`'s docstring now says no production route calls the FFT stage (checked by a grep of `src/`). The two short-window callers call the public `td_dynamical_correlator(..., x_values=[0], ks=[0.0])` instead of the transform directly. `test_sxt_to_skomega_stays_on_the_fft_stage` is replaced by `tests/test_audit_2026_09_24_realtime.py::test_sxt_to_skomega_matches_the_damped_trapezoid_sum` (delta*T = 6, relative error below 1e-11, and the FFT stage must miss by more than 0.1), with `tests/test_audit_2026_09_24b_realtime.py::test_sxt_to_skomega_is_the_direct_sum_at_every_window` on four windows. The short `nt` defaults are unchanged; raising them is the previous record's lead, not this finding. NUMBERS CHANGE for every infinite-chain S(k,omega): on a single magnon (L=16, dt=0.1, nt=1200, delta=0.05) the FFT stage was 2.103e-01 of the peak off the closed form with heights down to 0.797, and the direct sum is 2.0e-14 off; on the reviewer's saved TFIM window data (Sx,Sx) the k=pi line moves from +1.045 to +0.905 against eps = 0.900.

**Where**: `src/dmrgpy/timedependent.py:768` (`sxt_to_skomega` passes
`_evaluation="fft"` whatever `nt` it is given), reached from
`Infinite_Many_Body_Chain.td_dynamical_correlator` on v3 (`infinitechain.py:1542`)
and `pyitensor/idmrg_window.py:1461` on `"python"`; the docstring of
`_fourier_transform_correlator` (`timedependent.py:658-671`); the two
short-window callers `tests/test_infinite_chain.py:1091` and
`examples/idmrg/td_dynamical_correlator/main.py:121`, which also pass
`_evaluation="fft"`.

Finding 7 of the previous record measured that interpolating the damped transform
off an FFT grid of spacing 2*pi/(nt*dt) is 14 to 19 per cent of the peak wrong,
and `30200a4` replaced that stage by a direct per-frequency sum everywhere except
here, deferring `sxt_to_skomega` and the two short-window callers "until that
regime is measured". The deferral is not conditioned on delta*T: a caller who
converges the window (delta*T = 6, the finite-chain default) still gets a node
spacing of 1.047*delta and the full interpolation error, which falls only as
(pi/(delta*T))^2, so it reaches 1 per cent at delta*T of about 31. There is no
public way out: `_evaluation`, `predict` and `damping` are accepted on neither
route, and `factor` refines dt, not nt*dt. The docstring says two false things,
written in the same commit: that `sxt_to_skomega` is the switch's "one caller",
and that the two short-window callers "get the direct evaluation". It survived
because `test_sxt_to_skomega_stays_on_the_fft_stage` pins the FFT bits and no test
compares S(k,omega) with anything exact.

**Expected**: the direct sum, which on every series below equals the exact damped
trapezoid sum to roundoff.

Repro, from the hunter, on a single-magnon S(x,t) with a closed-form transform:

```bash
cd <scratch>/realtime && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 02_sxt_to_skomega_grid.py
```

`realtime/02_sxt_to_skomega_grid.py`:

```python
"""Lens realtime, lead 1b: sxt_to_skomega is hard-wired to
_evaluation="fft" whatever the caller's nt, so a caller who converges the
window (delta*T = 6, the finite-chain default) still gets the finding-7
FFT-grid interpolation error in S(k,omega).

Anchor, no dmrgpy in it: a single-magnon S(x,t) = (1/L) sum_q e^{iqx}
e^{+i eps(q) t} on L sites, eps(q) = 2 - 2 cos q, so that on the q grid
S(k,t) = sum_x e^{-ikx} S(x,t) = e^{+i eps(k) t} exactly and the one-sided
damped transform is (1/pi)(1 - e^{-(delta+i(w-eps))T})/(delta+i(w-eps)),
whose infinite-T limit has the density peak 1/(pi*delta) at w = eps(k)."""
import numpy as np
from dmrgpy import timedependent as td

L = 16
xs = np.arange(L) - L//2
qs = 2*np.pi*np.arange(L)/L
qs = np.where(qs > np.pi, qs - 2*np.pi, qs)
eps = 2 - 2*np.cos(qs)
ks = np.sort(qs)
epsk = 2 - 2*np.cos(ks)

def run(dt, nt, delta, es):
    ts = dt*np.arange(nt)
    S = np.exp(1j*np.outer(ts, eps)) @ (np.exp(1j*np.outer(qs, xs))/L)
    _k, _e, Skw = td.sxt_to_skomega(ts, xs, S, dt, ks=ks, es=es, delta=delta)
    T = ts[-1]
    z = delta + 1j*(es[None, :] - epsk[:, None])
    ex_T = (1 - np.exp(-z*T))/z/np.pi
    ex_inf = 1/z/np.pi
    peak = 1/(np.pi*delta)
    # the same series through the direct stage, for comparison
    Sd = np.array([td._fourier_transform_correlator(ts, S @ np.exp(-1j*k*xs),
                   dt, es=es, delta=delta)[1] for k in ks])
    pos_err = []
    pos_err_d = []
    for ik in range(len(ks)):
        pos_err.append(abs(es[np.argmax(Skw[ik].real)] - epsk[ik]))
        pos_err_d.append(abs(es[np.argmax(Sd[ik].real)] - epsk[ik]))
    print("dt=%.2f nt=%d delta=%.3f  delta*T=%.2f  FFT grid=%.3f=%.2f*delta"
          % (dt, nt, delta, delta*T, 2*np.pi/(nt*dt), 2*np.pi/(nt*dt)/delta))
    print("  sxt_to_skomega (fft):  max|S-exact_T|/peak = %.3e   "
          "max|ReS-exact_inf|/peak = %.3e   max|peak pos - eps(k)| = %.4f"
          % (np.max(np.abs(Skw - ex_T))/peak,
             np.max(np.abs(Skw.real - ex_inf.real))/peak, max(pos_err)))
    print("  direct stage:          max|S-exact_T|/peak = %.3e   "
          "max|ReS-exact_inf|/peak = %.3e   max|peak pos - eps(k)| = %.4f"
          % (np.max(np.abs(Sd - ex_T))/peak,
             np.max(np.abs(Sd.real - ex_inf.real))/peak, max(pos_err_d)))
    print("  min over k of the peak height, fft / exact_inf: %.3f"
          % (np.min(np.max(Skw.real, axis=1))/np.max(ex_inf.real)))

es = np.linspace(-1, 6, 1401)            # 0.005 spacing
# td_dynamical_correlator's own defaults: dt=0.1, nt=200 (and sxt_to_skomega's delta=5e-2)
run(0.1, 200, 5e-2, es)
# the same delta, window converged the way the finite-chain TD default does it
run(0.1, int(6/5e-2/0.1), 5e-2, es)
# the example's delta=0.15 at nt=40 dt=0.05, and at a converged window
run(0.05, 40, 0.15, es)
run(0.05, int(6/0.15/0.05), 0.15, es)
```

`realtime/02_sxt_to_skomega_grid.out`:

```
dt=0.10 nt=200 delta=0.050  delta*T=1.00  FFT grid=0.314=6.28*delta
  sxt_to_skomega (fft):  max|S-exact_T|/peak = 5.024e-01   max|ReS-exact_inf|/peak = 8.691e-01   max|peak pos - eps(k)| = 0.1522
  direct stage:          max|S-exact_T|/peak = 3.353e-04   max|ReS-exact_inf|/peak = 3.697e-01   max|peak pos - eps(k)| = 0.0022
  min over k of the peak height, fft / exact_inf: 0.145
dt=0.10 nt=1200 delta=0.050  delta*T=6.00  FFT grid=0.052=1.05*delta
  sxt_to_skomega (fft):  max|S-exact_T|/peak = 2.098e-01   max|ReS-exact_inf|/peak = 2.118e-01   max|peak pos - eps(k)| = 0.0228
  direct stage:          max|S-exact_T|/peak = 2.521e-04   max|ReS-exact_inf|/peak = 2.489e-03   max|peak pos - eps(k)| = 0.0022
  min over k of the peak height, fft / exact_inf: 0.795
dt=0.05 nt=40 delta=0.150  delta*T=0.29  FFT grid=3.142=20.94*delta
  sxt_to_skomega (fft):  max|S-exact_T|/peak = 2.302e-01   max|ReS-exact_inf|/peak = 9.585e-01   max|peak pos - eps(k)| = 1.2346
  direct stage:          max|S-exact_T|/peak = 2.700e-04   max|ReS-exact_inf|/peak = 7.464e-01   max|peak pos - eps(k)| = 0.0022
  min over k of the peak height, fft / exact_inf: 0.082
dt=0.05 nt=800 delta=0.150  delta*T=5.99  FFT grid=0.157=1.05*delta
  sxt_to_skomega (fft):  max|S-exact_T|/peak = 2.106e-01   max|ReS-exact_inf|/peak = 2.130e-01   max|peak pos - eps(k)| = 0.0778
  direct stage:          max|S-exact_T|/peak = 1.876e-04   max|ReS-exact_inf|/peak = 2.493e-03   max|peak pos - eps(k)| = 0.0022
  min over k of the peak height, fft / exact_inf: 0.790
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced the hunter's scripts bit for bit,
then decomposed the error: the direct stage equals an independent closed form of
the damped trapezoid sum to 1e-13, and `sxt_to_skomega` equals `np.interp` of that
same closed form at the FFT nodes to 7.5e-14, so kernel sign, 1/pi, damping and
endpoint weights all match and the 20 per cent is the interpolation alone, with
peak heights at the worst-case linear-interpolation bound:

`reviews/realtime_RA/11_interp_decomposition.py`:

```python
"""Reviewer R-A, probe 11: is sxt_to_skomega's error at a converged window
exactly the linear interpolation of the correct damped trapezoid sum on the
FFT grid, or does the FFT stage also differ in window/damping/weights?

Anchor built independently of the hunter's: the single-magnon S(x,t) on L
sites, S(k,t) = e^{+i eps(k) t} exactly on the q grid, and its damped
trapezoid sum in closed form (a geometric series),

  T(w) = (dt/pi) [ sum_{j<n} r^j - (1 + r^{n-1})/2 ],  r = e^{(i(eps-w)-delta) dt},

which is what _damped_sum_at should return to rounding. Then (b) the FFT
stage should be my own np.interp of T(w) at the FFT nodes 2 pi m/(n dt),
and (c) the worst-case peak of a Lorentzian of HWHM delta sampled at spacing
s and interpolated linearly is 1/(1+(s/2delta)^2) of its maximum.

Run through the public route's own defaults: es=None (window [-1,10], 800
points), factor=1, damping="exp", predict=False."""
import numpy as np
from dmrgpy import timedependent as td

L = 16
xs = np.arange(L) - L//2
qs = 2*np.pi*np.arange(L)/L
qs = np.where(qs > np.pi, qs - 2*np.pi, qs)
eps_q = 2 - 2*np.cos(qs)
ks = np.sort(qs)
epsk = 2 - 2*np.cos(ks)


def trapezoid_closed_form(w, eps, delta, dt, n):
    r = np.exp((1j*(eps - w) - delta)*dt)
    geo = (1 - r**n)/(1 - r)
    return dt/np.pi*(geo - 0.5*(1 + r**(n - 1)))


def run(label, dt, nt, delta):
    ts = dt*np.arange(nt)
    S = np.exp(1j*np.outer(ts, eps_q)) @ (np.exp(1j*np.outer(qs, xs))/L)
    _k, es, Skw = td.sxt_to_skomega(ts, xs, S, dt, ks=ks, delta=delta)
    grid = 2*np.pi/(nt*dt)
    m = np.arange(int(np.floor(-1.0/grid)) - 1, int(np.ceil(10.0/grid)) + 2)
    nodes = grid*m
    err_cf, err_interp, ratios, ratios_d, ratios_exact = [], [], [], [], []
    for ik, k in enumerate(ks):
        Skt = S @ np.exp(-1j*k*xs)
        _e, gd = td._fourier_transform_correlator(ts, Skt, dt, es=es,
                                                  delta=delta)
        cf = trapezoid_closed_form(es, epsk[ik], delta, dt, nt)
        err_cf.append(np.max(np.abs(gd - cf)))
        cfn = trapezoid_closed_form(nodes, epsk[ik], delta, dt, nt)
        mine = np.interp(es, nodes, cfn.real) + 1j*np.interp(es, nodes,
                                                             cfn.imag)
        err_interp.append(np.max(np.abs(Skw[ik] - mine)))
        # peak heights on the same es, against the infinite-T density
        # 1/(pi*delta) at w = eps(k)
        ratios.append(np.max(Skw[ik].real)*np.pi*delta)
        ratios_d.append(np.max(gd.real)*np.pi*delta)
        ratios_exact.append(np.max(cf.real)*np.pi*delta)
    s = grid/delta
    print("%s: dt=%.3f nt=%d delta=%.3f delta*T=%.3f grid=%.4f=%.3f*delta"
          % (label, dt, nt, delta, delta*ts[-1], grid, s))
    print("  max|direct stage - closed-form trapezoid sum|       = %.2e"
          % max(err_cf))
    print("  max|sxt_to_skomega - np.interp(closed form @ nodes)| = %.2e"
          % max(err_interp))
    print("  peak height / (1/(pi delta)), min over k: fft %.3f  direct "
          "%.3f  closed-form trapezoid on es %.3f"
          % (min(ratios), min(ratios_d), min(ratios_exact)))
    print("  worst-case linear-interp bound 1/(1+(s/2)^2) times the "
          "window's (1-e^{-delta T}) = %.3f"
          % ((1/(1 + (s/2)**2))*(1 - np.exp(-delta*ts[-1]))))


run("defaults (delta*T=1)", 0.1, 200, 0.05)
run("default dt, nt; delta=0.3", 0.1, 200, 0.3)
run("default dt, delta; nt=1200", 0.1, 1200, 0.05)
run("example heatmap call", 0.05, 40, 0.15)
```

`reviews/realtime_RA/11_interp_decomposition.out`:

```
defaults (delta*T=1): dt=0.100 nt=200 delta=0.050 delta*T=0.995 grid=0.3142=6.283*delta
  max|direct stage - closed-form trapezoid sum|       = 2.99e-14
  max|sxt_to_skomega - np.interp(closed form @ nodes)| = 1.78e-14
  peak height / (1/(pi delta)), min over k: fft 0.145  direct 0.629  closed-form trapezoid on es 0.629
  worst-case linear-interp bound 1/(1+(s/2)^2) times the window's (1-e^{-delta T}) = 0.058
default dt, nt; delta=0.3: dt=0.100 nt=200 delta=0.300 delta*T=5.970 grid=0.3142=1.047*delta
  max|direct stage - closed-form trapezoid sum|       = 4.46e-15
  max|sxt_to_skomega - np.interp(closed form @ nodes)| = 2.81e-15
  peak height / (1/(pi delta)), min over k: fft 0.797  direct 0.997  closed-form trapezoid on es 0.997
  worst-case linear-interp bound 1/(1+(s/2)^2) times the window's (1-e^{-delta T}) = 0.783
default dt, delta; nt=1200: dt=0.100 nt=1200 delta=0.050 delta*T=5.995 grid=0.0524=1.047*delta
  max|direct stage - closed-form trapezoid sum|       = 1.19e-13
  max|sxt_to_skomega - np.interp(closed form @ nodes)| = 7.51e-14
  peak height / (1/(pi delta)), min over k: fft 0.792  direct 0.980  closed-form trapezoid on es 0.980
  worst-case linear-interp bound 1/(1+(s/2)^2) times the window's (1-e^{-delta T}) = 0.783
example heatmap call: dt=0.050 nt=40 delta=0.150 delta*T=0.293 grid=3.1416=20.944*delta
  max|direct stage - closed-form trapezoid sum|       = 4.74e-15
  max|sxt_to_skomega - np.interp(closed form @ nodes)| = 2.14e-15
  peak height / (1/(pi delta)), min over k: fft 0.082  direct 0.254  closed-form trapezoid on es 0.254
  worst-case linear-interp bound 1/(1+(s/2)^2) times the window's (1-e^{-delta T}) = 0.002
```

On real v3 infinite-window data (TFIM paramagnet, SxSx, delta*T = 5.97) the FFT
stage is 0.197 of the peak off the direct stage, heights down to 0.813, lines
moved by up to 0.155 against half a node spacing of 0.157, while the direct stage
puts the line within 0.009 of the exact dispersion
(`reviews/realtime_RA/15_real_window_on_spectrum.py`; on which side of zero that
line sits is finding 8). The previous record's reason for deferring, that "the
direct stage moved the test's TD peak to the window edge", is wrong: on the test's
own seeded data the direct stage puts the line at omega = -2.75, below the test's
window, and the FFT's straight line had been hiding it
(`reviews/realtime_RA/17_test_call_wide_window.py`). Narrowed: at the short
windows every checked-in caller uses (delta*T at most 1) this is the previous
record's recorded lead, now measured rather than new (the direct stage is closer
to the delta-broadened transform in at least 39 of 40 seeded series at every
window from 0.29 to 6, but neither stage is converged there, and what that regime
needs is a longer `nt`); the example's "peaks up to 1.23 off at 8 per cent of the
exact height" is the synthetic anchor at the example's parameters, not a
measurement of the example; the test's vacuity (a two-node straight line passes
its >0.5 ratio for any KPM peak in [-1,3.81]) predates `30200a4`, and the test
also passes on the direct stage.

**Suggested fix**: drop `_evaluation="fft"` in `sxt_to_skomega` and correct the
docstring's two sentences; keep the switch, which
`test_direct_evaluation_is_the_fft_on_its_own_grid` uses as its reference; replace
`test_sxt_to_skomega_stays_on_the_fft_stage` by a pin of `sxt_to_skomega` against
the closed form at delta*T = 6, to about 1e-12. Raising the short `nt` defaults is
the recorded lead's fix and a separate decision (it triples the default cost). The
KPM cross-check test can only be made meaningful once finding 8 settles which side
of zero the line belongs on. NUMBERS CHANGE for every infinite-chain S(k,omega).

### 8. Every infinite-chain S(k,omega), on both backends, is mirrored in frequency: the lines sit at omega = -D_n where every finite route puts them at +D_n, so at the default `window=[-1,10]` the returned array is mostly the tail of lines lying below the window (9.5 to 53 per cent of a line's weight inside it), and the example's KPM cross-check passes only because of the mirror

`bug` &middot; severity **MEDIUM** (every value of the feature, one feature) &middot; CONFIRMED &middot; lens `realtime` (found by the reviewer of finding 7) &middot; origin the `idmrg_window` port

**Status**: FIXED, together with finding 7. `sxt_to_skomega` conjugates the momentum series after the spatial sum, `Skt = np.conj(S@phase)`, the finite route's conjugate-at-return at the one step both routes share, and its docstring and `dynamical_correlator_komega`'s now say so; nothing below the reduction changed. The two direct callers now go through the public route, so they carry the same conjugation, and both, the example's assertion included, moved to the local Sx,Sx anchor on the TFIM paramagnet. The fermionic example gained one sentence, and its code is unchanged. Conjugating before the sum passes the magnon's closed form, which is parity-symmetric, so the chiral ring is what pins the choice. Pinned by `tests/test_audit_2026_09_24b_realtime.py::test_sxt_to_skomega_momentum_label_and_sign_on_a_chiral_ring` and the rewritten `tests/test_infinite_chain.py::test_td_dynamical_correlator_agrees_qualitatively_with_kpm_finite` (both peaks in [0.9, 1.9], TD weight on omega<0 below 0.15, correlation with KPM above 0.8; passed three times on unseeded iDMRG). NUMBERS CHANGE for every infinite-chain S(k,omega) on both backends: on the TFIM paramagnet through the public route the `"python"` lines go from -2.095, -1.570 and -1.045 to +1.910, +1.495 and +0.905 at k = 0, pi/2 and pi (eps = 1.900, 1.487, 0.900), and the weight below zero from 0.971, 0.960 and 0.928 to 0.034, 0.047 and 0.083; v3 gives the same new lines. The example's cross-check now lands at TD +1.120 against KPM +1.060, and its tracked heatmap `td_dynamical_correlator_Skw.png` is regenerated.

**Where**: `src/dmrgpy/timedependent.py::sxt_to_skomega`, the reduction shared by
`Infinite_Many_Body_Chain.td_dynamical_correlator` on v3 and
`pyitensor/idmrg_window.py::dynamical_correlator_komega` on `"python"`; the two
direct callers of `_fourier_transform_correlator` that bypass it
(`tests/test_infinite_chain.py::test_td_dynamical_correlator_agrees_qualitatively_with_kpm_finite`
and the cross-check block of `examples/idmrg/td_dynamical_correlator/main.py`).

S(x,t) itself is right on both backends, <psi|A_x e^{-i(H-E0)t} B_0|psi>, as
`chain_session.h:246` and the `idmrg_window` docstrings say, and it is pinned
against an `expm` free-fermion oracle (`tests/test_idmrg_window_fermionic.py`).
The mirror enters in the reduction: the finite TD route evolves forward and
conjugates at return (`timedependent.py:201`), which is what puts its lines at
+D_n under the e^{-i omega t} kernel (ED gets the same by evolving with e^{+iHt},
finding 9's `evolve()`), and `idmrg_window` ported the energy shift H - E0 from
that route, "matching this codebase's established TD convention", but not the
conjugation. Nothing documents negative frequencies, and several places assert
the finite convention: the user guide's `td_dynamical_correlator` paragraph (the
reduction uses "the same damping/FFT tail submode="TD" uses", "take the real part
if you want the density", the density of `dynamics.py` with D_n = E_n - E_0 > 0),
`dynamical_correlator_komega`'s "delta means the same thing here as in every
other dynamical-correlator submode", and every default and example window. It
survived because every checked-in assertion compares shapes on symmetric grids or
compares the wrong quantities (below), and because the FFT stage of finding 7 had
been flattening the curve at the window edge, where the previous record's finding 7
Status took the line's absence for an artifact of the direct stage.

**Expected**: the lines at +D_n, as on every finite route.

Evidence, from the reviewer of finding 7 (the side lead that became this
candidate), on a TFIM paramagnet, SxSx, v3 infinite-window data:

`reviews/realtime_RA/15_real_window_on_spectrum.out`:

```
SxSx on es in [-2.5, 0.5]: max|fft-direct|/max|direct| = 0.197 (Re part 0.157)
   per-k peak height fft/direct in [0.813, 1.000], median 0.975
   |peak position fft - direct| up to 0.155 (grid/2 = 0.157)
   |direct peak - (-eps(k))| up to 0.009; |fft peak - (-eps(k))| up to 0.151   (es spacing 0.005)
   Re weight inside the default window [-1,10], as a fraction of the whole line's, per k, min/max: 0.076 / 0.614
SzSz on es in [-4.5, 0.5]: max|fft-direct|/max|direct| = 0.215 (Re part 0.175)
   per-k peak height fft/direct in [0.719, 1.000], median 0.982
   |peak position fft - direct| up to 1.155 (grid/2 = 0.157)
```

**Reviewer (CONFIRMED)**: built a method-independent anchor. For A = B = Sx the
slope of Im S(0,t) at t=0 is -<Sx(H-E0)Sx> = 0.7<Sz> exactly, by the double
commutator, which fixes the sign of the time series with no Fourier transform in
it; the infinite route has the e^{-i(H-E0)t} sign and the finite route's returned
series the e^{+i(H-E0)t} one, and the infinite route's local line sits at -1.115
where `kpm_finite` puts the same line at +1.055 (the shapes correlate at +0.978
with KPM's mirror image and at -0.400 with KPM itself). The ground-state phase is
applied correctly: the lines sit at -eps(k) to within 0.015 on both backends, with
no energy-density offset.

`reviews/realtime_RA2/01_python_infinite_anchor.py`:

```python
"""Reviewer R-A2, probe 01: where does the infinite-chain TD route put the
single-magnon line of the TFIM paramagnet, on itensor_version="python"?

Model H = 1.4 sum Sz + sum Sx Sx (tests/test_infinite_chain.py's model),
Pauli form 0.7 sz + 0.25 sx sx, deep paramagnet. Sx creates exactly one
quasiparticle, eps(k) = 2 sqrt(0.5525 + 0.35 cos k): 1.9 at k=0, 0.9 at
k=pi. Excitation energies D_n = E_n - E_0 > 0, so the house density
(dynamics.py) has support at omega = +eps(k) only.

Three independent reads:
 (a) the raw time series itself, no Fourier transform: for A = B = Sx,
     S(0,t) = sum_n |M_n| e^{-/+ i D_n t}, so d Im S(0,t)/dt at t=0 is
     -/+ <Sx (H-E0) Sx> = -/+ (-0.7 <Sz>) exactly (double commutator,
     only the field term fails to commute with Sx). The sign says which
     time convention the series is on, with no kernel involved.
 (b) S(k,omega) with the direct-stage transform (the default
     _evaluation="direct", so line positions are not the FFT
     interpolation R-A recorded), and the public td_dynamical_correlator
     (FFT stage) on a two-signed es.
 (c) kpm_finite on the same model, the same es, local r=0 and the TD
     route's local x=0 row: positive frequency by construction."""
import numpy as np
np.random.seed(0)
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_window
from dmrgpy.timedependent import _fourier_transform_correlator

ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
ic.gs_method = "idmrg"
ic.maxm, ic.maxiter, ic.etol, ic.niter = 20, 300, 1e-12, 150
ic.set_hamiltonian(1.4*ic.SzC[0] + ic.SxC[0]*ic.SxR[0])
e0 = ic.gs_energy()
sz = complex(ic.vev("Sz", 0)).real
print("e0 = %.12f   <Sz> = %.8f" % (e0, sz))

dt, nt, delta, n_window = 0.1, 120, 0.3, 10
xs_in = list(range(-4, 5))
ts, xs, S = idmrg_window.dynamical_correlator_td(
    ic._result, n_window, "Sx", "Sx", dt, nt, 1e-10, 30, niter=30,
    x_values=xs_in)
ts, xs, S = np.array(ts), np.array(xs), np.array(S)
i0 = list(xs).index(0)

# (a) small-t slope of Im S(0,t)
slope_exact = -0.7*sz          # <Sx (H-E0) Sx> > 0
slope_num = (S[1, i0].imag - S[0, i0].imag)/dt
print("(a) S(0,0) = %.6f%+.6fi   (static <Sx Sx> = %.6f)"
      % (S[0, i0].real, S[0, i0].imag, complex(ic.correlator("Sx", 0, "Sx", 0)).real))
print("    dIm S(0,t)/dt at t=0: numerical %.5f ; exact for e^{-i(H-E0)t} "
      "is %.5f, for e^{+i(H-E0)t} is %+.5f" % (slope_num, -slope_exact,
                                               slope_exact))

# (b) direct-stage S(k,omega), the transform done here by hand exactly as
# sxt_to_skomega does it, minus the _evaluation="fft"
es = np.linspace(-3, 3, 1201)
eps = lambda k: 2*np.sqrt(0.5525 + 0.35*np.cos(k))
ks = np.array([0.0, np.pi/2, np.pi])
print("(b) direct stage, delta=%.2f, delta*T=%.1f" % (delta, delta*ts[-1]))
for k in ks:
    Skt = S @ np.exp(-1j*k*xs)
    _e, g = _fourier_transform_correlator(ts, Skt, dt, es=es, delta=delta)
    y = g.real
    wneg = np.trapezoid(np.abs(y[es < 0]), es[es < 0])/np.trapezoid(np.abs(y), es)
    print("    k=%.3f  argmax Re S at w=%+.3f   (+eps=%+.3f, -eps=%+.3f)   "
          "|Re| weight on w<0 / total = %.3f" % (k, es[np.argmax(y)], eps(k),
                                                 -eps(k), wneg))

# the public route, FFT stage, two-signed es
kk, ee, Sk = ic.td_dynamical_correlator(
    "Sx", 0, "Sx", n_window=n_window, dt=dt, nt=nt, x_values=xs_in,
    maxdim=30, cutoff=1e-10, niter=30, ks=ks, es=es, delta=delta)
print("    public td_dynamical_correlator (FFT stage), same run parameters:")
for ik, k in enumerate(ks):
    y = Sk[ik].real
    wneg = np.trapezoid(np.abs(y[ee < 0]), ee[ee < 0])/np.trapezoid(np.abs(y), ee)
    print("    k=%.3f  argmax Re S at w=%+.3f   |Re| weight on w<0 / total = %.3f"
          % (k, ee[np.argmax(y)], wneg))
# how much of the line is inside the default window [-1,10], per k
kd, ed, Skd = ic.td_dynamical_correlator(
    "Sx", 0, "Sx", n_window=n_window, dt=dt, nt=nt, x_values=xs_in,
    maxdim=30, cutoff=1e-10, niter=30, ks=ks, delta=delta)
print("    default es (window=[-1,10]) spans [%.2f, %.2f]" % (ed[0], ed[-1]))
for ik, k in enumerate(ks):
    wwin = np.trapezoid(np.abs(Skd[ik].real), ed)
    wall = np.trapezoid(np.abs(Sk[ik].real), ee)
    print("    k=%.3f  |Re| weight inside default window / on [-3,3] = %.3f"
          % (k, wwin/wall))

# (c) local spectra: TD route at x=0 against kpm_finite r=0
_e, gloc = _fourier_transform_correlator(ts, S[:, i0], dt, es=es, delta=delta)
es_k, ykpm = ic.kpm_finite("Sx", 0, "Sx", 0, n_window=12,
                            window_chain_kwargs=dict(maxm=16, nsweeps=5),
                            delta=delta, es=es)
ykpm = np.asarray(ykpm).real
ytd = gloc.real
print("(c) local x=0:  kpm_finite argmax at w=%+.3f ;  TD route argmax at "
      "w=%+.3f  (local band is [0.9,1.9])" % (es_k[np.argmax(ykpm)],
                                             es[np.argmax(ytd)]))
def corr(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float(a @ b/np.sqrt((a @ a)*(b @ b)))
print("    shape correlation Re TD(w) with KPM(w) = %+.3f ; with KPM(-w) = %+.3f"
      % (corr(ytd, ykpm), corr(ytd, ykpm[::-1])))
print("    KPM weight on w<0 / total = %.3f ; TD |Re| weight on w<0 / total = %.3f"
      % (np.trapezoid(np.abs(ykpm[es < 0]), es[es < 0])/np.trapezoid(np.abs(ykpm), es),
         np.trapezoid(np.abs(ytd[es < 0]), es[es < 0])/np.trapezoid(np.abs(ytd), es)))
np.savez("01_S.npz", ts=ts, xs=xs, S=S, es=es, ykpm=ykpm, ytd=ytd)
```

`reviews/realtime_RA2/01_python_infinite_anchor.out`:

```
e0 = -0.722505349665   <Sz> = -0.48365300
(a) S(0,0) = 0.250000-0.000000i   (static <Sx Sx> = 0.250000)
    dIm S(0,t)/dt at t=0: numerical -0.33730 ; exact for e^{-i(H-E0)t} is -0.33856, for e^{+i(H-E0)t} is +0.33856
(b) direct stage, delta=0.30, delta*T=3.6
    k=0.000  argmax Re S at w=-1.910   (+eps=+1.900, -eps=-1.900)   |Re| weight on w<0 / total = 0.966
    k=1.571  argmax Re S at w=-1.495   (+eps=+1.487, -eps=-1.487)   |Re| weight on w<0 / total = 0.953
    k=3.142  argmax Re S at w=-0.905   (+eps=+0.900, -eps=-0.900)   |Re| weight on w<0 / total = 0.916
    public td_dynamical_correlator (FFT stage), same run parameters:
    k=0.000  argmax Re S at w=-2.095   |Re| weight on w<0 / total = 0.971
    k=1.571  argmax Re S at w=-1.570   |Re| weight on w<0 / total = 0.960
    k=3.142  argmax Re S at w=-1.045   |Re| weight on w<0 / total = 0.928
    default es (window=[-1,10]) spans [-1.00, 10.00]
    k=0.000  |Re| weight inside default window / on [-3,3] = 0.095
    k=1.571  |Re| weight inside default window / on [-3,3] = 0.182
    k=3.142  |Re| weight inside default window / on [-3,3] = 0.527
(c) local x=0:  kpm_finite argmax at w=+1.055 ;  TD route argmax at w=-1.115  (local band is [0.9,1.9])
    shape correlation Re TD(w) with KPM(w) = -0.400 ; with KPM(-w) = +0.978
    KPM weight on w<0 / total = 0.000 ; TD |Re| weight on w<0 / total = 0.943
```

v3 prints the same slope and lines at -1.915, -1.475 and -0.910
(`02_v3_infinite_anchor.out`); the finite routes on the same model land positive
(`03_finite_chain_anchor.out`: DMRG TD at +1.100, DMRG KPM at +1.060, ED TD at
+1.100). The test's KPM peak near omega=0 is the elastic line, since `kpm_finite`
does not subtract <A><B>, so the test compares a disconnected elastic line with a
connected mirrored continuum and passes under both signs and both stages; the
example's cross-check passes today only because of the mirror and fires once the
sign is fixed, its TD peak moving to the two-magnon continuum at +2.1 to +2.7
(`06_test_assertion_under_fix.out`, `07_example_assertion_under_fix.out`). Struck:
the neighbour's exemption for the fermionic example. For (Cdag_c, C_c) the house
convention has D_n = E_n(N-1) - E_0 = -eps_l > 0, so its band belongs at positive
omega; it only looked consistent because negative omega reads as a removal energy
on a symmetric grid, and since both of its panels go through `sxt_to_skomega` they
move together and its comparison survives the fix unchanged
(`04_fermionic_example_sign.out`).

**Suggested fix**: in `sxt_to_skomega`, conjugate after the spatial sum,
`Skt = np.conj(S@phase)`, the finite route's own conjugate-at-return at the one
step the two routes share, reaching both backends at once. Conjugating each x
series before the sum is the transform of the conjugate at -k, indistinguishable on
a parity-symmetric model and wrong on a chiral one:

`reviews/realtime_RA2/08_fix_momentum_label.py`:

```python
"""Reviewer R-A2, probe 08: where must the conjugation of a fix go, after
the spatial sum (np.conj(S @ phase)) or before it (np.conj(S) @ phase)?
The two coincide on a parity-symmetric model and differ by k -> -k on a
chiral one, so the discriminating anchor is a chiral free-fermion ring
(complex first-neighbour hopping t e^{i phi}, band minimum off k=0, Fermi
sea off centre), with an exact reference built from the Lehmann sum.

S(x,t) is built exactly as the IBC window documents it,
<0|A_x e^{-i(H-E0)t} B_0|0> = sum_n M_n(x) e^{-i D_n t}, for the pair
(A,B) = (Cdag_x, C_0): |n> = c_q|0>, D_q = -eps_q > 0,
M_q(x) = conj(v[x,q]) v[0,q] n_q. The reference is the house density of
the momentum series, sum_x e^{-ikx} sum_q M_q(x) delta/(pi((w-D_q)^2+delta^2)),
i.e. S(k,t) = sum_x e^{-ikx} S(x,t) with its own Lehmann lines at +D_q.
The whole ring is summed, so each k picks exactly its own states.
Errors are normalized by the largest reference value over all k shown
(a k with no occupied state has a reference of zero)."""
import numpy as np
from dmrgpy.timedependent import _fourier_transform_correlator

N, t, phi, mu = 40, 1.0, 0.7, -0.4
h = np.zeros((N, N), dtype=complex)
for x in range(N):
    h[(x+1) % N, x] += -t*np.exp(1j*phi)
    h[x, (x+1) % N] += -t*np.exp(-1j*phi)
h -= mu*np.eye(N)
eps, v = np.linalg.eigh(h)
occ = eps < 0
D = -eps[occ]
xs = np.arange(N)
M = (np.conj(v[:, occ])*v[0, occ][None, :])       # M[x, q]
print("occupied levels: %d of %d, D_q in [%.3f, %.3f]" % (occ.sum(), N,
                                                         D.min(), D.max()))

dt, nt, delta = 0.05, 1600, 0.2                    # delta*T = 16
ts = dt*np.arange(nt)
S = (M[None, :, :]*np.exp(-1j*D[None, None, :]*ts[:, None, None])).sum(-1)
es = np.linspace(-3.5, 3.5, 1401)
lor = delta/np.pi/((es[:, None] - D[None, :])**2 + delta**2)    # (w, q)

ms = (-7, -5, -3, 3, 5, 7)
refs = {}
for m in ms:
    k = 2*np.pi*m/N
    refs[m] = (lor*(np.exp(-1j*k*xs) @ M)[None, :]).sum(1)
scale = max(np.max(np.abs(r)) for r in refs.values())
worst = {"as is": 0.0, "conj after sum": 0.0, "conj before sum": 0.0}
for m in ms:
    k = 2*np.pi*m/N
    phase = np.exp(-1j*k*xs)
    ref = refs[m]
    variants = {"as is": S @ phase,
                "conj after sum": np.conj(S @ phase),
                "conj before sum": np.conj(S) @ phase}
    line = "k=%+.3f ref max %.3f at w=%+.3f |" % (k, np.max(np.abs(ref))/scale,
                                                   es[np.argmax(np.abs(ref))])
    for name, Skt in variants.items():
        _e, g = _fourier_transform_correlator(ts, Skt, dt, es=es, delta=delta)
        err = np.max(np.abs(g.real - ref.real))/scale
        worst[name] = max(worst[name], err)
        line += " %s: argmax %+.3f, err %.1e |" % (
            name, es[np.argmax(np.abs(g.real))], err)
    print(line)
print("worst max|Re - ref| / scale over the six k:",
      ", ".join("%s %.2e" % kv for kv in worst.items()))
```

`reviews/realtime_RA2/08_fix_momentum_label.out`:

```
occupied levels: 18 of 40, D_q in [0.054, 1.595]
k=-1.100 ref max 1.000 at w=+1.440 | as is: argmax -1.440, err 1.0e+00 | conj after sum: argmax +1.440, err 8.4e-06 | conj before sum: argmax +0.705, err 1.0e+00 |
k=-0.785 ref max 1.000 at w=+1.595 | as is: argmax -1.595, err 1.0e+00 | conj after sum: argmax +1.595, err 8.4e-06 | conj before sum: argmax +0.345, err 1.0e+00 |
k=-0.471 ref max 1.000 at w=+1.550 | as is: argmax -1.550, err 1.0e+00 | conj after sum: argmax +1.550, err 8.4e-06 | conj before sum: argmax +0.380, err 9.7e-01 |
k=+0.471 ref max 1.000 at w=+0.380 | as is: argmax -0.380, err 9.3e-01 | conj after sum: argmax +0.380, err 8.4e-06 | conj before sum: argmax +1.550, err 9.7e-01 |
k=+0.785 ref max 0.000 at w=+0.355 | as is: argmax -0.345, err 5.0e-16 | conj after sum: argmax +0.345, err 3.5e-16 | conj before sum: argmax +1.595, err 1.0e+00 |
k=+1.100 ref max 0.000 at w=+0.695 | as is: argmax -0.705, err 7.2e-16 | conj after sum: argmax +0.705, err 1.3e-16 | conj before sum: argmax +1.440, err 1.0e+00 |
worst max|Re - ref| / scale over the six k: as is 9.96e-01, conj after sum 8.45e-06, conj before sum 1.00e+00
```

Do not touch anything below the reduction (`dynamical_correlator_td`,
`snapshot_correlator`, `Chain::td_dynamical_correlator_window`), since S(x,t) is
right and pinned. Ship it with finding 7 (same function, same test). The two
direct callers need the conjugation themselves, or a shared one-series helper, and
both should move to the Sx,Sx local anchor; the example's assertion changes in the
same commit, and the fermionic example gains one sentence. Pin it with a
DMRG-free test on the chiral ring (the only model that sees the k label) and a sign
test on the local Sx,Sx anchor (weight on omega<0 below 0.1, today 0.943). Correct
the user guide's S(k,omega) paragraph in both formats, `documentation.md`'s "both
agree on where the dominant spectral weight sits" (true only by accident), the
closing paragraph of `dynamics.py`, and the docstrings of `sxt_to_skomega` and
`dynamical_correlator_komega`. NUMBERS CHANGE for every infinite-chain S(k,omega)
on both backends, the fermionic example's panels included.

### 9. `mode="ED"` `evolve_and_measure` and `evolution_ABA` evolve with e^{+iHt}, so ED returns the time-reversed trajectory, which a real Hamiltonian, a real start and a real observable cannot tell apart: Larmor precession comes out as <Sy>(t) = -sin(Bt)/2 against +sin(Bt)/2, a flux ring's density runs the wrong way round, and ED is the reference every DMRG evolution test is held to

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `realtime` &middot; origin `f54b3c2` (2020)

**Status**: FIXED, together with finding 10. `evolution_ABC` advances both of its states with `evolve(..., -Hop, ...)`, i.e. e^{-iHt}, and its docstring says so and why the minus sits there; `evolve()` and `scipy_evolution` are unchanged, since `evolution_DC` relies on e^{+iHt}. The comment at `examples/time_evolution/time_evolution_ABA/main.py:13` now reads <X|A^dagger e^{iHt} B e^{-iHt} A|X>. Pinned by `tests/test_audit_2026_09_24b_realtime.py::test_larmor_precession_runs_forward` (Sx, Sy and S+ against cos, sin and e^{+iBt}, on ED, `"python"` and v3) and `::test_flux_ring_circulates_forward`; against the committed code 24 of that file's 28 tests fail. NUMBERS CHANGE for `evolve_and_measure(mode="ED")` and `evolution_ABA(mode="ED")` wherever the Hamiltonian breaks time reversal, the start state is complex or the observable has imaginary matrix elements: Larmor precession (3 sites, B=1, +x start) <Sy_0>(t=1.0) goes from -0.4207 to +0.4207 (closed form +0.4207), and on the 3-site flux ring (phi = pi/6, mu = 3) <N_1>(t=1.5) from 0.1024 to 0.8413. Every existing test, which measures Sz or N under a real Hamiltonian, is unchanged.

**Where**: `src/dmrgpy/edtk/timedependent.py::evolution_ABC` (lines 35-36), which
advances both of its states with `edtk/tdtk.py::scipy_evolution`, integrating
dpsi/dt = +1j*h@psi (`tdtk.py:28`); reached by `evolve_and_measure(mode="ED")`
and `evolution_ABA(mode="ED")`. Not `evolution_DC`, which uses the same
`evolve()` on purpose (below).

The user guide's section 7 says |psi(t)> = e^{-iHt}|psi(0)> and that the call
returns <psi(t)|O|psi(t)>, and every DMRG backend follows it. ED returns
<psi(-t)|O|psi(-t)>, which equals the documented value only when the Hamiltonian,
the start state and the observable are all real and the observable Hermitian,
because then psi(-t) = psi(t)^*. Every ED-versus-DMRG evolution test and example
measures `Sz[i]` or `N[i]` from a real state under a real Hamiltonian, or from a
complex state whose pieces a conserved charge keeps apart, so none can see the
direction. The one place that writes the backward form, the comment at
`examples/time_evolution/time_evolution_ABA/main.py:13`, is stale against the
normative text. `30200a4` edited `evolution_ABC` (for finding 8 of the previous
record) without touching the sign.

**Expected**: the forward trajectory.

Repro, from the hunter, on the seed-3 complex-hopping fermion chain of
`tests/test_audit_2026_09_correlator-conventions.py`, against an independent kron
anchor:

```bash
cd <scratch>/realtime && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 04_ed_evolution_direction.py
```

`realtime/04_ed_evolution_direction.py`:

```python
"""Lens realtime, lead 2: edtk/tdtk.scipy_evolution integrates
dpsi/dt = +i H psi, i.e. psi(t) = e^{+iHt} psi, and edtk/timedependent.
evolution_ABC (which backs evolution_ABA(mode="ED") and
evolve_and_measure(mode="ED")) advances both of its states with it, so it
should return <O>(-t) where the DMRG side returns <O>(t). For a real
Hamiltonian with a real start the two coincide, which is every ED-versus-
DMRG evolution test; on a complex-hopping chain they do not.

Anchor, no dmrgpy in it: the audit's seed-3 complex-hopping chain built
again by hand as a Jordan-Wigner kron, psi(t) = expm(-iHt) Cdag_0|GS>, and
<psi(t)|N_2|psi(t)> at +t and at -t."""
import numpy as np
import scipy.linalg as sla
from dmrgpy import fermionchain, timedependent

n = 4
a = np.array([[0, 1], [0, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)

def cmat(i):
    out = np.array([[1.0 + 0j]])
    for k in range(n):
        out = np.kron(out, Z if k < i else (a if k == i else I2))
    return out

def build(itensor_version):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    rng = np.random.RandomState(3)
    t = rng.random((n, n)) + 1j * rng.random((n, n))
    t = t + t.conj().T
    h = 0
    for i in range(n):
        for j in range(n):
            h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1):
        h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 30, 12
    fc.tevol_method = "TDVP"
    return fc, t

fc, t = build("python")
C = [cmat(i) for i in range(n)]
Cd = [x.conj().T for x in C]
Nn = [Cd[i] @ C[i] for i in range(n)]
H = sum(t[i, j] * Cd[i] @ C[j] for i in range(n) for j in range(n))
H = H + sum(0.8 * Nn[i] @ Nn[i + 1] for i in range(n - 1))
e, U = np.linalg.eigh(H)
print("ground state gap e1-e0 = %.4f (unique)" % (e[1] - e[0]))
g = U[:, 0]
nt, dt = 40, 0.05
ts = dt * np.arange(nt)
psi0 = Cd[0] @ g

def exact(sign):
    out = []
    for tt in ts:
        p = sla.expm(-1j * sign * H * tt) @ psi0
        out.append(np.vdot(p, Nn[2] @ p))
    return np.array(out)

ex_p, ex_m = exact(+1), exact(-1)
print("exact: max|<N2>(t) - <N2>(-t)| = %.4e  (the size of the effect)"
      % np.max(np.abs(ex_p - ex_m)))

_t, y_ed = timedependent.evolution_ABA(fc, A=fc.Cdag[0], B=fc.N[2],
                                        mode="ED", nt=nt, dt=dt)
y_ed = np.asarray(y_ed)
print("mode=ED               : max|y - exact(+t)| = %.3e   max|y - exact(-t)| = %.3e"
      % (np.max(np.abs(y_ed - ex_p)), np.max(np.abs(y_ed - ex_m))))
for v in ("python", 3):
    fcv, _ = build(v)
    _t, y = timedependent.evolution_ABA(fcv, A=fcv.Cdag[0], B=fcv.N[2],
                                         mode="DMRG", nt=nt, dt=dt)
    y = np.asarray(y)
    print("mode=DMRG, itensor=%-6s: max|y - exact(+t)| = %.3e   max|y - exact(-t)| = %.3e"
          % (v, np.max(np.abs(y - ex_p)), np.max(np.abs(y - ex_m))))

# the same through evolve_and_measure with an explicit start: the ED
# ground state itself, measured on a current-like operator, J = i(Cd0 C1 - Cd1 C0)
J = 1j * (fc.Cdag[0] * fc.C[1] - fc.Cdag[1] * fc.C[0])
Jm = 1j * (Cd[0] @ C[1] - Cd[1] @ C[0])
phi0 = (Cd[1] + Cd[2]) @ g
phi0 = phi0 / np.linalg.norm(phi0)
ex = {s: np.array([np.vdot(sla.expm(-1j*s*H*tt) @ phi0, Jm @ (sla.expm(-1j*s*H*tt) @ phi0))
                   for tt in ts]) for s in (+1, -1)}
ed = fc.get_ED_obj()
from dmrgpy.edtk.edchain import State
# express phi0 in the ED object's own basis via its own operators
gs_ed = ed.get_gs_array()
Ad = np.array(ed.MO2matrix(fc.Cdag[1] + fc.Cdag[2]).todense())
phi_ed = Ad @ gs_ed
phi_ed = phi_ed / np.linalg.norm(phi_ed)
_t, yj = timedependent.evolve_and_measure(fc, operator=J, mode="ED", nt=nt,
                                           dt=dt, wf=State(phi_ed, ed))
yj = np.asarray(yj)
print("evolve_and_measure(mode=ED), <J01>(t): max|y - exact(+t)| = %.3e   "
      "max|y - exact(-t)| = %.3e   (max|exact| = %.3f)"
      % (np.max(np.abs(yj - ex[+1])), np.max(np.abs(yj - ex[-1])),
         np.max(np.abs(ex[+1]))))
```

`realtime/04_ed_evolution_direction.out`:

```
ground state gap e1-e0 = 0.1219 (unique)
exact: max|<N2>(t) - <N2>(-t)| = 4.3520e-01  (the size of the effect)
mode=ED               : max|y - exact(+t)| = 4.352e-01   max|y - exact(-t)| = 3.618e-08
mode=DMRG, itensor=python: max|y - exact(+t)| = 1.031e-13   max|y - exact(-t)| = 4.352e-01
mode=DMRG, itensor=3     : max|y - exact(+t)| = 1.230e-06   max|y - exact(-t)| = 4.352e-01
evolve_and_measure(mode=ED), <J01>(t): max|y - exact(+t)| = 1.144e+00   max|y - exact(-t)| = 2.599e-07   (max|exact| = 0.717)
```

and on a real Hamiltonian with a real start, Larmor precession of a
+x-polarized chain:

`realtime/07_precession_direction.py`:

```python
"""Lens realtime, lead 2b: the same time-direction question on a REAL
Hamiltonian with a REAL start, where it still shows, through an operator
with imaginary matrix elements. Larmor precession has a closed form: for
H = B sum_i Sz_i and a start polarized along +x, the Heisenberg equations
give <Sx_i>(t) = cos(Bt)/2 and <Sy_i>(t) = +sin(Bt)/2 (dSy/dt = i[H,Sy] =
B Sx). The quench is written the way examples/time_evolution/* write it:
ground state of h0 = -sum Sx, then set_hamiltonian(h1) and
evolve_and_measure from that state."""
import numpy as np
from dmrgpy import spinchain, timedependent

n, B, nt, dt = 3, 1.0, 40, 0.1
ts = dt*np.arange(nt)
exact_sy = 0.5*np.sin(B*ts)

def run(v, mode):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    sc.maxm, sc.nsweeps = 10, 8
    sc.tevol_method = "TDVP"
    h0 = 0
    for i in range(n): h0 = h0 - sc.Sx[i]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs(mode=mode)
    h1 = 0
    for i in range(n): h1 = h1 + B*sc.Sz[i]
    sc.set_hamiltonian(h1)
    _t, y = timedependent.evolve_and_measure(sc, operator=sc.Sy[0], nt=nt,
                                              dt=dt, wf=wf, mode=mode)
    _t, x = timedependent.evolve_and_measure(sc, operator=sc.Sx[0], nt=nt,
                                              dt=dt, wf=wf, mode=mode)
    return np.asarray(y).real, np.asarray(x).real

for v, mode in (("python", "ED"), ("python", "DMRG"), (3, "DMRG")):
    y, x = run(v, mode)
    print("mode=%-4s itensor=%-6s  <Sx0>: max|y-cos/2| = %.2e   "
          "<Sy0>: max|y-(+sin/2)| = %.2e   max|y-(-sin/2)| = %.2e   <Sy0>(t=1.0) = %+.4f"
          % (mode, v, np.max(np.abs(x-0.5*np.cos(B*ts))),
             np.max(np.abs(y-exact_sy)), np.max(np.abs(y+exact_sy)), y[10]))
print("closed form <Sy0>(t=1.0) = %+.4f" % exact_sy[10])
```

`realtime/07_precession_direction.out`:

```
mode=ED   itensor=python  <Sx0>: max|y-cos/2| = 2.40e-08   <Sy0>: max|y-(+sin/2)| = 1.00e+00   max|y-(-sin/2)| = 1.39e-08   <Sy0>(t=1.0) = -0.4207
mode=DMRG itensor=python  <Sx0>: max|y-cos/2| = 1.81e-14   <Sy0>: max|y-(+sin/2)| = 1.43e-14   max|y-(-sin/2)| = 1.00e+00   <Sy0>(t=1.0) = +0.4207
mode=DMRG itensor=3       <Sx0>: max|y-cos/2| = 5.55e-16   <Sy0>: max|y-(+sin/2)| = 6.66e-16   max|y-(-sin/2)| = 1.00e+00   <Sy0>(t=1.0) = +0.4207
closed form <Sy0>(t=1.0) = +0.4207
```

**Reviewer (CONFIRMED)**: nothing struck. Built two anchors of its own with the
Schrodinger sign fixed explicitly, after ruling out that ED's `Sy` or its `H`
were sign-flipped (`[Sx,Sy] = +i Sz` on ED's own matrices, and `H_ED` is exactly
B sum Sz): free spins, where `tdtk.evolve` on ED's own H matches e^{+iHt} to 1.8e-7
and misses e^{-iHt} by 0.69 (`reviews/realtime_RB/03_larmor_own_anchor.py`), and
one fermion on a 3-site ring with flux pi/6, where the Hamiltonian breaks time
reversal, the observable N_k is real, and time reversal swaps N_1 and N_2 (v2's
0.157 is its MPO-Taylor error on the forward side; the 1e-6 ED residuals are
`solve_ivp`'s tolerances):

`reviews/realtime_RB/04_flux_ring_own_anchor.py`:

```python
"""Reviewer anchor B: one fermion on a 3-site ring with a flux.

H = -t sum_j (e^{i phi} Cdag_{j+1} C_j + h.c.) + mu sum_j N_j, mu = 3 > 2t,
so the vacuum is the unique ground state (E0 = 0, gap mu - 2t = 1) and
Cdag_0|vac> is one particle on site 0. In the one-particle sector the
Jordan-Wigner strings only cross empty sites, so the 3x3 matrix
h_sp[j+1,j] = -t e^{i phi}, h_sp[j,j+1] = -t e^{-i phi}, h_sp[j,j] = mu
is exact, and <N_k>(t) = |(expm(-1j*h_sp*t) e_0)_k|^2 with the sign of the
Schrodinger equation fixed here. The flux makes the particle circulate one
way, so N_1(t) != N_2(t) and time reversal swaps them: a TR-breaking H
seen through a real observable. No dmrgpy in the anchor."""
import numpy as np
import scipy.linalg as sla
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import fermionchain, timedependent

n, tt, phi, mu = 3, 1.0, np.pi/6, 3.0
nt, dt = 40, 0.1
ts = dt*np.arange(nt)
hsp = np.diag([mu]*n).astype(complex)
for j in range(n):
    hsp[(j+1) % n, j] += -tt*np.exp(1j*phi)
    hsp[j, (j+1) % n] += -tt*np.exp(-1j*phi)
e0 = np.zeros(n, complex); e0[0] = 1
def exact(sign):
    return np.array([np.abs(sla.expm(-1j*sign*hsp*t)@e0)**2 for t in ts])  # (nt, n)
fw, bw = exact(+1), exact(-1)
print("anchor: max|N1(t)-N2(t)| forward = %.3f ; max|N_k(t)-N_k(-t)| = %.3f"
      % (np.max(np.abs(fw[:, 1]-fw[:, 2])), np.max(np.abs(fw-bw))))

def build(v):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=v)
    h = 0
    for j in range(n):
        h = h - tt*np.exp(1j*phi)*fc.Cdag[(j+1) % n]*fc.C[j]
        h = h - tt*np.exp(-1j*phi)*fc.Cdag[j]*fc.C[(j+1) % n]
    for j in range(n): h = h + mu*fc.N[j]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 20, 10
    fc.tevol_method = "TDVP"
    return fc

rows = [("ED", "python")] + [("DMRG", v) for v in ("python", 3, 2)]
for mode, v in rows:
    fc = build(v)
    print("%-4s v=%-6s E0 = %+.2e" % (mode, v, fc.gs_energy(mode=mode)))
    for k in (1, 2):
        _t, y = timedependent.evolution_ABA(fc, A=fc.Cdag[0], B=fc.N[k], mode=mode, nt=nt, dt=dt)
        y = np.asarray(y)
        print("%-4s v=%-6s <N_%d>: max|y-exact(+t)| = %.2e   max|y-exact(-t)| = %.2e   y(t=1.5)=%.4f  (+t: %.4f, -t: %.4f)"
              % (mode, v, k, np.max(np.abs(y-fw[:, k])), np.max(np.abs(y-bw[:, k])),
                 y[15].real, fw[15, k], bw[15, k]))
```

`reviews/realtime_RB/04_flux_ring_own_anchor.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
anchor: max|N1(t)-N2(t)| forward = 1.000 ; max|N_k(t)-N_k(-t)| = 1.000
ED   v=python E0 = +0.00e+00
ED   v=python <N_1>: max|y-exact(+t)| = 1.00e+00   max|y-exact(-t)| = 4.45e-06   y(t=1.5)=0.1024  (+t: 0.8413, -t: 0.1024)
ED   v=python <N_2>: max|y-exact(+t)| = 1.00e+00   max|y-exact(-t)| = 2.68e-06   y(t=1.5)=0.8413  (+t: 0.1024, -t: 0.8413)
DMRG v=python E0 = +9.91e-16
DMRG v=python <N_1>: max|y-exact(+t)| = 2.04e-14   max|y-exact(-t)| = 1.00e+00   y(t=1.5)=0.8413  (+t: 0.8413, -t: 0.1024)
DMRG v=python <N_2>: max|y-exact(+t)| = 2.34e-14   max|y-exact(-t)| = 1.00e+00   y(t=1.5)=0.1024  (+t: 0.1024, -t: 0.8413)
DMRG v=3      E0 = +1.79e-25
DMRG v=3      <N_1>: max|y-exact(+t)| = 4.88e-15   max|y-exact(-t)| = 1.00e+00   y(t=1.5)=0.8413  (+t: 0.8413, -t: 0.1024)
DMRG v=3      <N_2>: max|y-exact(+t)| = 5.44e-15   max|y-exact(-t)| = 1.00e+00   y(t=1.5)=0.1024  (+t: 0.1024, -t: 0.8413)
DMRG v=2      E0 = +1.54e-25
DMRG v=2      <N_1>: max|y-exact(+t)| = 1.57e-01   max|y-exact(-t)| = 1.04e+00   y(t=1.5)=0.8181  (+t: 0.8413, -t: 0.1024)
DMRG v=2      <N_2>: max|y-exact(+t)| = 1.45e-01   max|y-exact(-t)| = 1.07e+00   y(t=1.5)=0.1533  (+t: 0.1024, -t: 0.8413)
```

Measured, not asserted, that `evolve()` itself must not change sign:
`evolution_DC`'s e^{-i omega t} kernel is built on e^{+iHt}, and flipping it puts
ED TD at 0.226 from an exact reference on complex weights, against 3.2e-5 as it is
(`reviews/realtime_RB/07_evolution_DC_needs_plus_iHt.py`). Widened the NUMBERS
CHANGE clause: ED also moves on a complex start state and on a non-Hermitian
observable with a complex expectation value, even under a real Hamiltonian.

**Suggested fix**: in `evolution_ABC` alone, advance both states with
`evolve(..., -Hop, ...)`; conjugating the returned array instead cannot work, since
on the flux ring the observable is real. Under an in-process patch every anchor
lands forward and the 20 consumer tests pass
(`reviews/realtime_RB/05_fix_under_patch.py`, `06_consumer_tests_under_patch.out`).
Correct the stale comment and `evolution_ABC`'s docstring, which leaves the sign of
U open. Land it in the same commit as finding 10, because either fix alone makes
ED and DMRG disagree on S+ where today they agree by cancellation. Regressions:
the Larmor <Sy_0> and the flux-ring <N_1>, both closed forms. NUMBERS CHANGE for
`evolve_and_measure(mode="ED")` and `evolution_ABA(mode="ED")` wherever the
Hamiltonian breaks time reversal, the start state is complex, or the observable
has imaginary matrix elements.

### 10. DMRG `evolve_and_measure` conjugates its result on return, so on every backend and every integrator it returns <psi(t)|O^dagger|psi(t)> instead of <psi(t)|O|psi(t)>: a bond current reads with the wrong sign, and for O = Sz_0 + i*Sx_0 on an eigenstate at t=0 it returns -0.5i where `vev()` on the same state returns +0.5i

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `realtime` (found by the reviewer of finding 9) &middot; origin `9e33f32` (2020)

**Status**: FIXED, together with finding 9. `evolve_and_measure_dmrg` returns `cs` unconjugated in both of its branches, in `timedependent.py` and in `mpsjulialive/timedependent.py`; `evolution_dmrg_DC` keeps its conjugation, and the docstring now states the return as <psi(t)|O|psi(t)> and says why the correlator route conjugates and this one does not. Pinned by `tests/test_audit_2026_09_24b_realtime.py::test_evolve_and_measure_on_an_eigenstate_is_vev` (ED, and `"python"` and v3 on TDVP, TDVP_GSE, TEBD, AUTO and MPO, v2 on MPO, `julia_live` on TDVP and TEBD), `::test_evolve_and_measure_return_wf_branch_is_not_conjugated` and `::test_evolution_aba_on_an_eigenstate_is_vev`. NUMBERS CHANGE for `evolve_and_measure(mode="DMRG")` and `evolution_ABA(mode="DMRG")` on any observable with a complex expectation value, by exactly twice its imaginary part, on v2, v3, `"python"` and `julia_live`: on the +x state of -sum Sx (3 sites), O = Sz_0 + i*Sx_0 goes from -0.5i to +0.5i, and ED stays at +0.5i. Hermitian observables do not move.

**Where**: `src/dmrgpy/timedependent.py:288-289` (`evolve_and_measure_dmrg`
returns `cs.real-1j*cs.imag`, in both its plain and its `return_wf` branch) and the
same line at `src/dmrgpy/mpsjulialive/timedependent.py:149-150`; reached by
`evolve_and_measure(mode="DMRG")` and `evolution_ABA(mode="DMRG")`.

Every session method measures <psi|O|psi> (`inner(psi,A,psi)` on `"python"`,
`innerC` on v3, v2's MPO loop, Julia's `inner(psi,Aop,psi)`), and the wrapper
conjugates. The line was copied in `9e33f32` from `evolution_dmrg`'s return, where
the conjugation carries the correlator convention (it still does, at line 201),
and the pybind port `d8f7cc5` carried it over. On a Hermitian observable the
imaginary part is roundoff, so the flip is invisible, which is why no test sees it.
It contradicts the user guide's section 7, the Julia loop's own comment ("the
returned correlator <psi(t)|Aop|psi(t)>"), and `vev()` on the same state, so the
library disagrees with itself at t=0 with no propagator involved.

**Expected**: <psi(t)|O|psi(t)>.

Evidence, from the reviewer who owned this candidate: part A is an eigenstate
anchor with no propagator (the +x state of H = -sum Sx on 3 sites, O = Sz_0 +
1j*Sx_0, exactly +0.5j at every t); part B an anchor on bare Pauli matrices; part C
ED as it is and under finding 9's fix; part D free fermions, <Cdag_0 C_1> against
a single-particle anchor:

`reviews/realtime_RB2/01_conj_anchor.py`:

```python
"""Reviewer of R-B2: does evolve_and_measure(mode="DMRG") return
conj(<psi(t)|O|psi(t)>) for a non-Hermitian O?

Part A, no propagator in the anchor at all. H = -sum_i Sx_i on n=3 spins,
whose unique ground state is the +x product state, an eigenstate of H, so
<O>(t) is constant for every t whatever the sign of the propagator. For
O = Sz_0 + 1j*Sx_0 the exact value is 0 + 1j*(1/2) = +0.5j, by hand.
Printed for every tevol_method on "python" and v3, MPO on v2, together with
the raw session correlator (the list before timedependent.py's return
line), sc.vev(O) on the same state, and mode="ED".

Part B, Larmor: H1 = B sum Sz_i, start +x, own numpy anchor on Pauli
matrices with psi(t) = expm(-1j*H*t) psi0 written explicitly (no dmrgpy),
<S+_0>(t) and <Sy_0>(t).

Part C, the consequence sub-claim: mode="ED" as it is and under an
in-process copy of the R-B fix (evolution_ABC fed -Hop, identical to the
neighbouring reviewer's rb_patch.py), next to DMRG unpatched.

Part D, free fermions: Cdag_0 C_1 after a quench on n=4, own single-particle
anchor G(t) = conj(U) G0 U^T with U = expm(-1j*h1*t)."""
import sys
import numpy as np
import scipy.linalg as sla
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain, fermionchain, timedependent
from dmrgpy.edtk import timedependent as tded
from dmrgpy.edtk.tdtk import evolve
from dmrgpy.edtk.edchain import State

np.set_printoptions(precision=6, suppress=True)

def fmt(z):
    return "%+.6f%+.6fj" % (z.real, z.imag)

# ---------------------------------------------------------------- Part A
print("=== Part A: eigenstate anchor, O = Sz_0 + 1j*Sx_0, exact +0.5j at all t")
n, nt, dt = 3, 6, 0.1
methods = {"python": ("TDVP", "TDVP_GSE", "TEBD", "AUTO", "MPO"),
           3: ("TDVP", "TDVP_GSE", "TEBD", "AUTO", "MPO"),
           2: ("MPO", "TDVP")}
for v in ("python", 3, 2):
    for meth in methods[v]:
        sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
        sc.maxm, sc.nsweeps = 10, 10
        sc.tevol_method = meth
        h0 = 0
        for i in range(n): h0 = h0 - sc.Sx[i]
        sc.set_hamiltonian(h0)
        wf = sc.get_gs(mode="DMRG")
        O = sc.Sz[0] + 1j*sc.Sx[0]
        vv = complex(sc.vev(O))
        try:
            _t, y = timedependent.evolve_and_measure(sc, operator=O, nt=nt, dt=dt, wf=wf)
        except Exception as e:
            print("v=%-6s %-8s evolve_and_measure raised %r" % (v, meth, e))
            continue
        y = np.asarray(y)
        # the raw session list, before the return line conjugates it
        s = sc._session
        if v in (3, "python") and meth in ("TDVP", "AUTO"):
            raw, _w = s.evolve_and_measure_tdvp(sc.hamiltonian.to_terms(), O.to_terms(),
                                                wf.cpp_handle, nt, dt)
        elif v in (3, "python") and meth == "TEBD":
            raw, _w = s.evolve_and_measure_tebd(sc.hamiltonian.to_terms(), O.to_terms(),
                                                wf.cpp_handle, nt, dt)
        elif v in (3, "python") and meth == "TDVP_GSE":
            raw, _w = s.evolve_and_measure_tdvp_gse(sc.hamiltonian.to_terms(), O.to_terms(),
                                                    wf.cpp_handle, nt, dt, sc.tdvp_gse_sweeps,
                                                    sc.tdvp_gse_krylov_order, sc.tdvp_gse_cutoff)
        else:
            raw, _w = s.evolve_and_measure(sc.hamiltonian.to_terms(), O.to_terms(),
                                           wf.cpp_handle, nt, dt, False)
        raw = np.asarray(raw)
        print("v=%-6s %-8s vev=%s  raw session[0]=%s  evolve_and_measure[0]=%s  "
              "max|y-(+0.5j)|=%.1e  max|y-(-0.5j)|=%.1e"
              % (v, meth, fmt(vv), fmt(raw[0]), fmt(y[0]),
                 np.max(np.abs(y-0.5j)), np.max(np.abs(y+0.5j))))
# ED on the same chain
sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
h0 = 0
for i in range(n): h0 = h0 - sc.Sx[i]
sc.set_hamiltonian(h0)
O = sc.Sz[0] + 1j*sc.Sx[0]
_t, y = timedependent.evolve_and_measure(sc, operator=O, nt=nt, dt=dt, mode="ED")
y = np.asarray(y)
print("mode=ED           vev(ED)=%s  evolve_and_measure[0]=%s  max|y-(+0.5j)|=%.1e  max|y-(-0.5j)|=%.1e"
      % (fmt(complex(sc.vev(O, mode="ED"))), fmt(y[0]), np.max(np.abs(y-0.5j)), np.max(np.abs(y+0.5j))))

# ---------------------------------------------------------------- Part B/C
print("=== Part B: Larmor, own numpy anchor on Pauli matrices")
B, nt, dt = 1.0, 40, 0.1
ts = dt*np.arange(nt)
sx = 0.5*np.array([[0, 1], [1, 0]], dtype=complex)
sy = 0.5*np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = 0.5*np.array([[1, 0], [0, -1]], dtype=complex)
sp = sx + 1j*sy
print("anchor algebra ||[sx,sy]-i sz|| = %.1e" % np.linalg.norm(sx@sy-sy@sx-1j*sz))
p0 = np.array([1, 1], dtype=complex)/np.sqrt(2)
H1 = B*sz
ps = [sla.expm(-1j*H1*t)@p0 for t in ts]
a_sp = np.array([np.vdot(p, sp@p) for p in ps])
a_sy = np.array([np.vdot(p, sy@p) for p in ps]).real
a_szix = np.array([np.vdot(p, (sz+1j*sx)@p) for p in ps])
print("anchor: max|<S+>-e^{+iBt}/2| = %.1e   max|<Sy>-sin(Bt)/2| = %.1e   <S+>(t=1) = %s"
      % (np.max(np.abs(a_sp-0.5*np.exp(1j*B*ts))), np.max(np.abs(a_sy-0.5*np.sin(B*ts))), fmt(a_sp[10])))

def larmor_chain(v, n=3):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    sc.maxm, sc.nsweeps = 10, 10
    sc.tevol_method = "TDVP" if v != 2 else "MPO"
    h0 = 0
    for i in range(n): h0 = h0 - sc.Sx[i]
    sc.set_hamiltonian(h0)
    return sc

def h_field(sc, n=3):
    h1 = 0
    for i in range(n): h1 = h1 + B*sc.Sz[i]
    return h1

obs = lambda sc: (("S+", sc.Sx[0]+1j*sc.Sy[0], a_sp), ("Sy", sc.Sy[0], a_sy),
                  ("Sz+iSx", sc.Sz[0]+1j*sc.Sx[0], a_szix))
res = {}
for v in ("python", 3, 2):
    sc = larmor_chain(v)
    wf = sc.get_gs(mode="DMRG")
    sc.set_hamiltonian(h_field(sc))
    for lab, op, a in obs(sc):
        _t, y = timedependent.evolve_and_measure(sc, operator=op, nt=nt, dt=dt, wf=wf)
        y = np.asarray(y)
        res[(v, lab)] = y
        print("DMRG v=%-6s <%-6s>: max|y-anchor| = %.2e   max|y-conj(anchor)| = %.2e   y(t=1) = %s"
              % (v, lab, np.max(np.abs(y-a)), np.max(np.abs(y-np.conj(a))), fmt(y[10])))

print("=== Part C: mode=ED as it is and under the R-B fix, same Larmor setup")

def evolution_ABC_fixed(self, h, A=None, B=None, C=None, wf=None, nt=100, dt=0.01):
    # copy of the neighbouring reviewer's rb_patch.evolution_ABC_fixed
    nt = int(nt)
    Aop = self.get_operator(A); Bop = self.get_operator(B); Cop = self.get_operator(C)
    ts = np.array([dt*ii for ii in range(nt)])
    Hop = self.get_operator(h)
    if wf is None: wf = self.get_gs_array()
    elif type(wf) == State: wf = wf.v
    wfA = Aop@wf; wfC = Cop@wf
    cs = []
    for it in range(nt):
        cs.append(np.conjugate(wfC)@Bop@wfA)
        wfA = evolve(wfA, -Hop, t=dt, dt=dt)
        wfC = evolve(wfC, -Hop, t=dt, dt=dt)
    return ts, np.array(cs)

_orig = tded.evolution_ABC
for patch in ("as is", "R-B fix"):
    tded.evolution_ABC = _orig if patch == "as is" else evolution_ABC_fixed
    sc = larmor_chain("python")
    ed = sc.get_ED_obj()
    psi0 = ed.get_gs_array()
    sc.set_hamiltonian(h_field(sc))
    ed = sc.get_ED_obj()
    for lab, op, a in obs(sc):
        _t, y = timedependent.evolve_and_measure(sc, operator=op, nt=nt, dt=dt,
                                                wf=State(psi0, ed), mode="ED")
        y = np.asarray(y)
        yd = res[("python", lab)]
        print("ED %-7s <%-6s>: max|y-anchor| = %.2e  max|y-conj(anchor)| = %.2e  max|ED-DMRG(python)| = %.2e"
              % (patch, lab, np.max(np.abs(y-a)), np.max(np.abs(y-np.conj(a))), np.max(np.abs(y-yd))))
tded.evolution_ABC = _orig

# ---------------------------------------------------------------- Part D
print("=== Part D: free fermions, Cdag_0 C_1 after a quench, own single-particle anchor")
nf, nt, dt = 4, 30, 0.1
ts = dt*np.arange(nt)
mu0 = np.array([0.6, -0.6, 0.6, -0.6])
h0m = np.zeros((nf, nf)); h1m = np.zeros((nf, nf))
for i in range(nf-1):
    h0m[i, i+1] = h0m[i+1, i] = -1.0
    h1m[i, i+1] = h1m[i+1, i] = -1.0
h0m += np.diag(mu0)
h1m += np.diag([0.0, 0.3, 0.0, 0.0])
e, P = np.linalg.eigh(h0m)
print("h0 single-particle levels", e, "(no zero level, unique Fock ground state)")
occ = P[:, e < 0]
G0 = np.conj(occ)@occ.T      # G0[k,l] = <c_k^dag c_l>
Gt = []
for t in ts:
    U = sla.expm(-1j*h1m*t)
    Gt.append(np.conj(U)@G0@U.T)
a01 = np.array([g[0, 1] for g in Gt])
print("anchor <Cdag0 C1>(0) = %s, max|Im| over t = %.3f, <Cdag0 C1>(t=1) = %s"
      % (fmt(a01[0]), np.max(np.abs(a01.imag)), fmt(a01[10])))

def fchain(v):
    fc = fermionchain.Fermionic_Chain(nf, itensor_version=v)
    fc.maxm, fc.nsweeps = 20, 12
    fc.tevol_method = "TDVP" if v != 2 else "MPO"
    return fc

def build(fc, hm):
    h = 0
    for i in range(nf):
        for j in range(nf):
            if abs(hm[i, j]) > 0: h = h + hm[i, j]*fc.Cdag[i]*fc.C[j]
    return h

resf = {}
for v in ("python", 3):
    fc = fchain(v)
    fc.set_hamiltonian(build(fc, h0m))
    wf = fc.get_gs(mode="DMRG")
    op = fc.Cdag[0]*fc.C[1]
    v0 = complex(fc.vev(op))
    fc.set_hamiltonian(build(fc, h1m))
    _t, y = timedependent.evolve_and_measure(fc, operator=op, nt=nt, dt=dt, wf=wf)
    y = np.asarray(y)
    resf[v] = y
    print("DMRG v=%-6s <Cdag0 C1>: vev(t=0)=%s  max|y-anchor| = %.2e  max|y-conj(anchor)| = %.2e  y(t=1)=%s"
          % (v, fmt(v0), np.max(np.abs(y-a01)), np.max(np.abs(y-np.conj(a01))), fmt(y[10])))
for patch in ("as is", "R-B fix"):
    tded.evolution_ABC = _orig if patch == "as is" else evolution_ABC_fixed
    fc = fchain("python")
    fc.set_hamiltonian(build(fc, h0m))
    psi0 = fc.get_ED_obj().get_gs_array()
    fc.set_hamiltonian(build(fc, h1m))
    ed = fc.get_ED_obj()
    op = fc.Cdag[0]*fc.C[1]
    _t, y = timedependent.evolve_and_measure(fc, operator=op, nt=nt, dt=dt,
                                            wf=State(psi0, ed), mode="ED")
    y = np.asarray(y)
    print("ED %-7s <Cdag0 C1>: max|y-anchor| = %.2e  max|y-conj(anchor)| = %.2e  max|ED-DMRG(python)| = %.2e"
          % (patch, np.max(np.abs(y-a01)), np.max(np.abs(y-np.conj(a01))), np.max(np.abs(y-resf["python"]))))
tded.evolution_ABC = _orig
```

`reviews/realtime_RB2/01_conj_anchor.out`:

```

maxLinkDim after global subspace expansion = 1
Global subspace expansion: cputime = 0.00100s, walltime = 0.00100s

maxLinkDim after global subspace expansion = 1
Global subspace expansion: cputime = 0.000954s, walltime = 0.000955s

maxLinkDim after global subspace expansion = 1
Global subspace expansion: cputime = 0.000869s, walltime = 0.000869s

maxLinkDim after global subspace expansion = 1
Global subspace expansion: cputime = 0.000867s, walltime = 0.000866s

maxLinkDim after global subspace expansion = 1
Global subspace expansion: cputime = 0.000899s, walltime = 0.000898s

maxLinkDim after global subspace expansion = 1
Global subspace expansion: cputime = 0.000962s, walltime = 0.000962s
dmrgpy from <repo>/src/dmrgpy/__init__.py
=== Part A: eigenstate anchor, O = Sz_0 + 1j*Sx_0, exact +0.5j at all t
v=python TDVP     vev=+0.000000+0.500000j  raw session[0]=+0.000000+0.500000j  evolve_and_measure[0]=+0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=6.0e-16
v=python TDVP_GSE vev=+0.000000+0.500000j  raw session[0]=+0.000000+0.500000j  evolve_and_measure[0]=+0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=2.6e-15
v=python TEBD     vev=+0.000000+0.500000j  raw session[0]=+0.000000+0.500000j  evolve_and_measure[0]=+0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=7.8e-16
v=python AUTO     vev=-0.000000+0.500000j  raw session[0]=-0.000000+0.500000j  evolve_and_measure[0]=-0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=6.0e-16
v=python MPO      vev=-0.000000+0.500000j  raw session[0]=-0.000000+0.500000j  evolve_and_measure[0]=-0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=6.0e-04
v=3      TDVP     vev=-0.000000+0.500000j  raw session[0]=-0.000000+0.500000j  evolve_and_measure[0]=-0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=3.1e-16
v=3      TDVP_GSE vev=+0.000000+0.500000j  raw session[0]=+0.000000+0.500000j  evolve_and_measure[0]=+0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=4.3e-16
v=3      TEBD     vev=+0.000000+0.500000j  raw session[0]=+0.000000+0.500000j  evolve_and_measure[0]=+0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=1.2e-15
v=3      AUTO     vev=-0.000000+0.500000j  raw session[0]=-0.000000+0.500000j  evolve_and_measure[0]=-0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=9.7e-16
v=3      MPO      vev=-0.000000+0.500000j  raw session[0]=-0.000000+0.500000j  evolve_and_measure[0]=-0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=6.0e-04
v=2      MPO      vev=+0.000000+0.500000j  raw session[0]=+0.000000+0.500000j  evolve_and_measure[0]=+0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=6.0e-04
v=2      TDVP     vev=+0.000000+0.500000j  raw session[0]=+0.000000+0.500000j  evolve_and_measure[0]=+0.000000-0.500000j  max|y-(+0.5j)|=1.0e+00  max|y-(-0.5j)|=6.0e-04
mode=ED           vev(ED)=-0.000000+0.500000j  evolve_and_measure[0]=-0.000000+0.500000j  max|y-(+0.5j)|=5.9e-09  max|y-(-0.5j)|=1.0e+00
=== Part B: Larmor, own numpy anchor on Pauli matrices
anchor algebra ||[sx,sy]-i sz|| = 0.0e+00
anchor: max|<S+>-e^{+iBt}/2| = 2.3e-16   max|<Sy>-sin(Bt)/2| = 1.7e-16   <S+>(t=1) = +0.270151+0.420735j
DMRG v=python <S+    >: max|y-anchor| = 1.00e+00   max|y-conj(anchor)| = 1.68e-14   y(t=1) = +0.270151-0.420735j
DMRG v=python <Sy    >: max|y-anchor| = 1.07e-14   max|y-conj(anchor)| = 1.07e-14   y(t=1) = +0.420735-0.000000j
DMRG v=python <Sz+iSx>: max|y-anchor| = 1.00e+00   max|y-conj(anchor)| = 1.70e-14   y(t=1) = +0.000000-0.270151j
DMRG v=3      <S+    >: max|y-anchor| = 1.00e+00   max|y-conj(anchor)| = 5.74e-16   y(t=1) = +0.270151-0.420735j
DMRG v=3      <Sy    >: max|y-anchor| = 5.00e-16   max|y-conj(anchor)| = 5.00e-16   y(t=1) = +0.420735+0.000000j
DMRG v=3      <Sz+iSx>: max|y-anchor| = 1.00e+00   max|y-conj(anchor)| = 1.96e-15   y(t=1) = -0.000000-0.270151j
DMRG v=2      <S+    >: max|y-anchor| = 1.00e+00   max|y-conj(anchor)| = 5.70e-03   y(t=1) = +0.269017-0.421653j
DMRG v=2      <Sy    >: max|y-anchor| = 4.94e-03   max|y-conj(anchor)| = 4.94e-03   y(t=1) = +0.421653+0.000000j
DMRG v=2      <Sz+iSx>: max|y-anchor| = 1.00e+00   max|y-conj(anchor)| = 3.54e-03   y(t=1) = -0.000146-0.269017j
=== Part C: mode=ED as it is and under the R-B fix, same Larmor setup
ED as is   <S+    >: max|y-anchor| = 1.00e+00  max|y-conj(anchor)| = 2.66e-08  max|ED-DMRG(python)| = 2.66e-08
ED as is   <Sy    >: max|y-anchor| = 1.00e+00  max|y-conj(anchor)| = 1.00e+00  max|ED-DMRG(python)| = 1.00e+00
ED as is   <Sz+iSx>: max|y-anchor| = 2.40e-08  max|y-conj(anchor)| = 1.00e+00  max|ED-DMRG(python)| = 1.00e+00
ED R-B fix <S+    >: max|y-anchor| = 2.66e-08  max|y-conj(anchor)| = 1.00e+00  max|ED-DMRG(python)| = 1.00e+00
ED R-B fix <Sy    >: max|y-anchor| = 1.39e-08  max|y-conj(anchor)| = 1.39e-08  max|ED-DMRG(python)| = 1.39e-08
ED R-B fix <Sz+iSx>: max|y-anchor| = 2.40e-08  max|y-conj(anchor)| = 1.00e+00  max|ED-DMRG(python)| = 1.00e+00
=== Part D: free fermions, Cdag_0 C_1 after a quench, own single-particle anchor
h0 single-particle levels [-1.725698 -0.861374  0.861374  1.725698] (no zero level, unique Fock ground state)
anchor <Cdag0 C1>(0) = +0.370094+0.000000j, max|Im| over t = 0.219, <Cdag0 C1>(t=1) = +0.326148-0.161149j
DMRG v=python <Cdag0 C1>: vev(t=0)=+0.370094+0.000000j  max|y-anchor| = 4.38e-01  max|y-conj(anchor)| = 5.85e-14  y(t=1)=+0.326148+0.161149j
DMRG v=3      <Cdag0 C1>: vev(t=0)=+0.370094+0.000000j  max|y-anchor| = 4.38e-01  max|y-conj(anchor)| = 7.20e-13  y(t=1)=+0.326148+0.161149j
ED as is   <Cdag0 C1>: max|y-anchor| = 4.38e-01  max|y-conj(anchor)| = 1.24e-07  max|ED-DMRG(python)| = 1.24e-07
ED R-B fix <Cdag0 C1>: max|y-anchor| = 1.24e-07  max|y-conj(anchor)| = 4.38e-01  max|ED-DMRG(python)| = 4.38e-01
```

**Reviewer (CONFIRMED, NARROWED)**: the raw session list and `vev` are right on
every row and the wrapper flips the sign, on `"python"`, v2, v3 and `julia_live`
(the latter in `reviews/realtime_RB2/03_aba_and_julia.py`, which also covers
`evolution_ABA`, the only consumer in `src/`). Narrowed the neighbour's "ED and
DMRG agree on S+ only because both are wrong": that holds only for observables
with real matrix elements (S+, Cdag_0 C_1), where ED's e^{+iHt} happens to produce
the same conjugate; for O = Sz_0 + 1j*Sx_0 the two disagree today by 1.00 on a
value of 0.5, and ED is the right one, so this defect is observable on its own.
No test, example or benchmark measures a non-Hermitian observable through either
call.

**Suggested fix**: return `cs` unconjugated in both return statements, in both
files; leave `evolution_dmrg_DC`'s conjugation, the session methods and the ED
side alone. With the fix applied in-process the seven test files that call the two
functions pass (112 passed, the largest imaginary part returned 1.94e-16;
`reviews/realtime_RB2/02_consumer_tests_under_fix.out`). Land it with finding 9
(with both in, every row of part C and D lands on its anchor). Regression: part
A's eigenstate anchor, `evolve_and_measure(...)[0] == vev(O)` on every backend
and on ED. Say in the docstring why `evolution_dmrg_DC` conjugates and this one
does not. NUMBERS CHANGE for any observable with a complex expectation value, by
exactly twice its imaginary part.

### 11. On every DMRG backend, `set_gs()` reaches only the Python-side `wf0`, and every ground-state-reading correlator submode (KPM, CVM, ROOTN, TD, TDZ) returns the density of the session's own solved state and then overwrites `wf0` with it: 0.204 off the requested member's exact density on a 0.445 peak, while `mode="ED"` honours the same `set_gs` to 3.3e-6

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `realtime` &middot; origin `a8233ab` (the pybind port)

**Status**: FIXED, with findings 12 to 14 in one change. The public `set_gs`, `set_initial_wf` and `set_initial_wf_guess` go through `groundstate.mark_injected()`, a mark that holds the injected object and counts only while `self.wf0` is that object, so a later solve or `restart()` retires it; `gs_is_current` is False while a mark is pending; `gs_energy_single` hands the session a detached copy of an injected state and, under skip, takes it unswept with e0 = <wf|H|wf>. `dynamics.py:190` is now `groundstate.ground_state_on_session()`, which makes no session call when the stored state is current and the Hamiltonian is on the session, and otherwise solves and sends. Two things the reviewers did not see had to go in with it: the session fills its lower band edge lazily through `gs_energy(skip_dmrg=True)`, which after a push would have swept the pushed state in place (|<x|session>|^2 0.989 on `"python"`, 0.987 on v3, measured), so the edges are filled before the push (`session.excited_states(1,1.0,False)`, no C++ change; one reduced solve when the edges are missing); and `__deepcopy__` now resets the clone's ground state, since the copied `wf0` carries the original session's site indices. Pinned by `tests/test_audit_2026_09_24b_groundstate.py::test_every_submode_measures_the_member_set_gs_set` (KPM, CVM, CVM_explicit, ROOTN, TD, TDZ and EX, on `"python"` and v3), `::test_injected_state_is_not_swept_by_the_first_kpm_call`, `::test_a_clone_of_a_solved_chain_runs_a_correlator` and `::test_a_repeated_correlator_on_a_solved_chain_does_not_resweep` (zero session calls, zero sweeps); against the committed code 37 of that file's 40 tests fail. NUMBERS CHANGE only after `set_gs` or `set_initial_wf`: on the 3-site Heisenberg chain after `set_gs` of the pure 2Sz=+1 member, TD goes from 2.036e-01 to 3.3e-06 off the member's exact density on `"python"` (from 1.505e-01 on v3), CVM and ROOTN from 2.036e-01 to 1e-15, and KPM from 2.787e-01 to 9.3e-04 off ED's KPM (from 4.538e-01 on v3), the up/down difference going from 0 to ED's 1.317. See findings 1, 2, 8, 9 and 11 of `audit_2026_09_24c_hole_hunt.md`: KPM kept measuring a set state from the solved E0 while the other submodes took its own energy, `mode="ED"` `submode="ED"` and SECTOR never read the set state, the band-edge pre-fill cost an upper-edge solve and a discarded fluctuation on every first read, and `ground_state_on_session` at the top of the dispatcher ran ahead of the KPM argument checks; all five are fixed.

**Where**: `src/dmrgpy/dynamics.py:190` (`self.set_initial_wf(self.wf0)` at the
top of every DMRG correlator, which sets `computed_gs=False`),
`src/dmrgpy/groundstate.py:175` (`gs_energy_single` calls
`session.set_wavefunction` only when handed a `wf0=` argument), `:178`
(`session.gs_energy(skip_dmrg=True)` returns the session's own state) and `:185`
(overwrites `self.wf0` with it). KPM and TD read the session's state directly and
would ignore `set_gs` even without line 190; CVM, ROOTN and TDZ read the
Python-side `wf0` and are defeated by line 190 alone.

Before the pybind port, line 190 was the trigger of a hand-off:
`taskdmrg.py::write_tasks` wrote `gs_from_file`/`starting_file_gs`, and
`mpscpp2/get_gs.h` returned that state unswept under `skip_dmrg_gs`, to the KPM
correlator (`dyncorr.h`) and to time evolution. The port kept the trigger and
dropped the hand-off. `gs_is_current`'s docstring (`groundstate.py:55`) says an
injected state "is returned unconditionally ... re-solving over it would discard
exactly what they asked to use", which is contradicted, and its premise, that such
a state carries no `_gs_solver_key`, is false (the key after `set_gs` is the
solve's own). The effect is that `vev`, `get_gs` and `evolution_ABA(mode="DMRG")`
honour `set_gs` only until the first correlator call, after which an unrelated
earlier call has changed their result. It matters wherever one member of a
degenerate manifold is chosen by hand; its size is the whole distance between the
requested member's density and the solved state's, set by whichever superposition
the random start reached (on v3 `np.random.seed` does not reach `randomMPS`, so it
varies between runs).

**Expected**: the density in the state the caller set, as on `mode="ED"`.

Repro, from the hunter, on the 3-site Heisenberg doublet with `set_gs` of the pure
2Sz=+1 member:

```bash
cd <scratch>/realtime && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 08d_td_after_set_gs_pure_member.py
```

`realtime/08d_td_after_set_gs_pure_member.py`:

```python
"""Lens realtime, lead 4b sized: set_gs() a PURE member of the 3-site
Heisenberg doublet on mode="DMRG" and on mode="ED", then TD of (Sp_0,
Sm_2), against the exact density in that member (dense eigh, no dmrgpy
algebra in the reference). The DMRG member is built from the solved state
as (Sz_tot + 1/2)|s>, normalized, which is exactly its 2Sz=+1 component."""
import numpy as np
np.random.seed(1)
from dmrgpy import spinchain
from dmrgpy.edtk.edchain import State

n = 3
ES = np.linspace(-1, 3, 41); delta = 0.3
kw = dict(submode="TD", es=ES, delta=delta, dt=0.05)

def chain(v):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 10, 10
    return sc

sc = chain("python")
ed = sc.get_ED_obj()
H = np.array(ed.get_hamiltonian().todense())
Sp0 = np.array(ed.MO2matrix(sc.Sx[0] + 1j*sc.Sy[0]).todense())
Sm2 = np.array(ed.MO2matrix(sc.Sx[2] - 1j*sc.Sy[2]).todense())
Szt = np.array(ed.MO2matrix(sc.Sz[0] + sc.Sz[1] + sc.Sz[2]).todense())
e, U = np.linalg.eigh(H)
P = U[:, :2]
w, V = np.linalg.eigh(P.conj().T @ Szt @ P)
up = P @ V[:, 1]
M = (up.conj() @ Sp0 @ U) * (U.conj().T @ Sm2 @ up)
ref_up = ((delta/np.pi)/((ES[:, None]-(e-e[0])[None, :])**2+delta**2)) @ M
print("exact density in the 2Sz=+1 member: peak %.4f" % np.max(np.abs(ref_up)))

sc.get_gs(mode="ED")  # solve first: set_gs on a fresh ED object leaves e0=None
sc.set_gs(State(up, ed))
y = np.asarray(sc.get_dynamical_correlator(mode="ED", name=[sc.Sx[0] + 1j*sc.Sy[0],
        sc.Sx[2] - 1j*sc.Sy[2]], **kw)[1])
print("mode=ED   after set_gs(up member): max|y - exact| = %.3e" % np.max(np.abs(y - ref_up)))

sd = chain("python")
sd.get_gs()
Szd = sd.Sz[0] + sd.Sz[1] + sd.Sz[2]
wu = (Szd + 0.5)*sd.get_gs()
wu = wu*(1/np.sqrt(wu.dot(wu).real))
sd.set_gs(wu)
print("DMRG state after set_gs: <Sz_tot> = %+.4f" % sd.vev(Szd).real)
nm = [sd.Sx[0] + 1j*sd.Sy[0], sd.Sx[2] - 1j*sd.Sy[2]]
y = np.asarray(sd.get_dynamical_correlator(mode="DMRG", name=nm, **kw)[1])
print("mode=DMRG after set_gs(up member): max|y - exact| = %.3e   (DMRG state now <Sz_tot> = %+.4f)"
      % (np.max(np.abs(y - ref_up)), sd.vev(Szd).real))
sd.set_gs(wu)
y = np.asarray(sd.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=nm,
        es=ES, delta=delta, dt=0.05)[1])
print("mode=DMRG TDZ after set_gs(up member): max|y - exact| = %.3e" % np.max(np.abs(y - ref_up)))
```

`realtime/08d_td_after_set_gs_pure_member.out`:

```
exact density in the 2Sz=+1 member: peak 0.4446
mode=ED   after set_gs(up member): max|y - exact| = 3.323e-06
DMRG state after set_gs: <Sz_tot> = +0.5000
mode=DMRG after set_gs(up member): max|y - exact| = 2.036e-01   (DMRG state now <Sz_tot> = +0.0349)
mode=DMRG TDZ after set_gs(up member): max|y - exact| = 2.044e-01
```

**Reviewer (CONFIRMED, NARROWED)**: widened to every ground-state-reading submode
with a discriminator that needs no reference: `set_gs` the up member, run the
submode, `set_gs` the down member and run again; ED separates them, DMRG returns
the solved state's curve both times (v3 and v2 print the same columns in
`04_every_submode_after_set_gs_v3.out` and `_v2.out`):

`reviews/realtime_RC/04_every_submode_after_set_gs.py`:

```python
"""Reviewer probe: does ANY mode="DMRG" dynamical-correlator submode honour
set_gs()?  Discriminator that needs no reference and no method-accuracy
argument: set_gs the pure 2Sz=+1 member of the 3-site Heisenberg doublet,
run submode S, then set_gs the 2Sz=-1 member and run S again.  The two
members' (S+_0, S-_2) densities differ at the level of the peak, so a
submode that honours set_gs returns two different curves, and one that
ignores it returns the same curve twice.  The ED columns give the size of
that difference on mode="ED" (exact, set_gs honoured), and the last
column checks that the DMRG curve is the one of the *solved* state,
y(solved) computed on the same chain before any set_gs.
Backend from argv[1]: python | 3 | 2."""
import sys
import numpy as np
np.random.seed(1)
from dmrgpy import spinchain
from dmrgpy.edtk.edchain import State

v = sys.argv[1]
v = v if v == "python" else int(v)
n = 3
ES = np.linspace(-1, 3, 21); delta = 0.3


def chain():
    sc = spinchain.Spin_Chain([2]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 10, 10
    return sc


KW = {"KPM": dict(delta=delta), "CVM": dict(delta=delta),
      "CVM_explicit": dict(delta=delta), "ROOTN": dict(delta=delta),
      "TD": dict(delta=delta, dt=0.05), "TDZ": dict(delta=delta, dt=0.05),
      "EX": dict(delta=delta, nex=6)}

# exact members, for the ED columns
se = chain()
ed = se.get_ED_obj()
H = np.array(ed.get_hamiltonian().todense())
Szt = np.array(ed.MO2matrix(se.Sz[0] + se.Sz[1] + se.Sz[2]).todense())
e, U = np.linalg.eigh(H)
P = U[:, :2]
w, V = np.linalg.eigh(P.conj().T @ Szt @ P)
up_e, dn_e = P @ V[:, 1], P @ V[:, 0]
se.get_gs(mode="ED")
nm_e = [se.Sx[0] + 1j*se.Sy[0], se.Sx[2] - 1j*se.Sy[2]]

sd = chain()
print("itensor_version=%r, mode actually used: %s" % (v, sd.get_mode(mode="DMRG")))
sd.get_gs()
Szd = sd.Sz[0] + sd.Sz[1] + sd.Sz[2]
s0 = sd.get_gs().copy()
print("solved state <Sz_tot> = %+.4f" % sd.vev(Szd).real)
up = (Szd + 0.5)*s0; up = up*(1/np.sqrt(up.dot(up).real))
dn = (0.5 - 1*Szd)*s0; dn = dn*(1/np.sqrt(dn.dot(dn).real))
nm = [sd.Sx[0] + 1j*sd.Sy[0], sd.Sx[2] - 1j*sd.Sy[2]]

print("%-13s %12s %12s %12s %14s" % ("submode", "ED |up-dn|", "DMRG |up-dn|",
      "|DMRG_up-solv|", "<Sz> after call"))
for sub in ["KPM", "CVM", "CVM_explicit", "ROOTN", "TD", "TDZ", "EX"]:
    kw = dict(submode=sub, es=ES, **KW[sub])
    try:
        if sub in ("KPM", "CVM", "TD", "EX", "ROOTN"):
            se.set_gs(State(up_e, ed))
            yeu = np.asarray(se.get_dynamical_correlator(mode="ED", name=nm_e, **kw)[1])
            se.set_gs(State(dn_e, ed))
            yed = np.asarray(se.get_dynamical_correlator(mode="ED", name=nm_e, **kw)[1])
            ded = "%.3e" % np.max(np.abs(yeu - yed))
        else:
            ded = "n/a"
        sd.set_gs(s0)
        ys = np.asarray(sd.get_dynamical_correlator(mode="DMRG", name=nm, **kw)[1])
        sd.set_gs(up)
        yu = np.asarray(sd.get_dynamical_correlator(mode="DMRG", name=nm, **kw)[1])
        szu = sd.vev(Szd).real
        sd.set_gs(dn)
        yd = np.asarray(sd.get_dynamical_correlator(mode="DMRG", name=nm, **kw)[1])
        print("%-13s %12s %12.3e %12.3e %+14.4f" % (sub, ded, np.max(np.abs(yu - yd)),
              np.max(np.abs(yu - ys)), szu))
    except Exception as ex:
        print("%-13s raised %s: %s" % (sub, type(ex).__name__, str(ex)[:150]))
```

`reviews/realtime_RC/04_every_submode_after_set_gs_python.out`:

```
itensor_version='python', mode actually used: DMRG
solved state <Sz_tot> = +0.0349
submode         ED |up-dn| DMRG |up-dn| |DMRG_up-solv| <Sz> after call
KPM              1.318e+00    0.000e+00    0.000e+00        +0.0349
CVM              4.378e-01    0.000e+00    0.000e+00        +0.0349
Non Hermitian mode in dynamical correlator
CVM_explicit  raised NotImplementedError: get_dynamical_correlator: submode='CVM_explicit' is only implemented for a Hermitian operator pair, A^dagger == B (it inverts (z-H) against the single
ROOTN            4.378e-01    0.000e+00    0.000e+00        +0.0349
TD               4.378e-01    0.000e+00    0.000e+00        +0.0349
TDZ                    n/a    0.000e+00    0.000e+00        +0.0349
EX               0.000e+00    0.000e+00    0.000e+00        +0.0349
```

0.4651 times 0.4378 is the hunter's 0.2036 exactly; once the session itself is
handed the member, DMRG TD separates the members by ED's 0.4378 on both backends,
so the ignoring is not a limitation of the method. Struck: "TD and TDZ" (it is
every submode that reads the ground state, so the claim is about `set_gs`); "the
ED/DMRG split that finding 8's fix opened" (ED KPM, CVM and ROOTN honoured
`set_gs` before `30200a4`; finding 8 only brought ED TD in line); EX, which reads
its own reference on both modes (finding 15). KPM's half was recorded in the
previous record's finding 12 Status ("set_gs alone leaves KPM on the solved
state"); new are the revert of the Python-side state, the uniformity across
submodes, and the mechanism.

**Suggested fix**: both of the hunter's variants are wrong. "Have `gs_energy_single`
send `self.wf0` when it has no solver key" rests on the false premise above; "have
`set_gs` push to the session" works for exact eigenstates only, since
`set_wavefunction` drops the session's cached energy and `skip_dmrg=True` then
means "sweep", so the next correlator call relaxes an injected state (finding 13
is exactly this, done by hand), and on `"python"` it also mutates the caller's own
MPS (`reviews/realtime_RC/05_fix_variant_push_to_session.py`). Restore the
pre-port contract: the public `set_gs`/`set_initial_wf`/`set_initial_wf_guess`
mark the state as injected (the solver's own trailing calls at `groundstate.py:185`,
`:376`, `:390` and `nhdmrg.py:297` must not); `gs_energy_single` pushes a real copy
of an injected state to the session and, under `skip_dmrg_gs`, takes e0 =
<wf|H|wf> rather than the session's cached energy, as the pre-port
`get_gs_energy` did; `reconverge=True` sweeps from it; and line 190 becomes "solve
only when the state is not current, and make sure the session has the Hamiltonian"
(`send_hamiltonian`), since a clone of a solved chain keeps `computed_gs=True` with
an empty session and a bare deletion crashes it (finding 13's reviewer). One fix
closes findings 11, 12, 13 and 14. NUMBERS CHANGE only after `set_gs` or
`set_initial_wf`.

### 12. `set_initial_wf` and `set_initial_wf_guess`, the documented way to re-enter `gs_energy()` with a state in hand, never reach the DMRG session on any backend, so the warm start does nothing and the transverse-field Ising example, which seeds the ferromagnetic branch at every field, plots Mz/n of 0.001 across the ordered phase on v3 and v2 where the seeded branch gives up to 0.4998

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `realtime` (found by the reviewer of finding 11) &middot; origin `a8233ab` (the pybind port)

**Status**: FIXED, with finding 11's change. `set_initial_wf_guess` sweeps from a detached copy of the guess, `set_initial_wf` returns the state unswept, and the explicit `gs_energy(wf0=x)` pushes a copy too, so the caller's `x` is no longer mutated on `"python"`; `Many_Body_Chain.gs_energy` no longer returns its stored energy ahead of a `wf0=` argument, and a stored state whose Hamiltonian is no longer the session's is solved again. `set_initial_wf_guess`'s docstring and `tests/test_bond_dimension_ramp.py`'s module docstring are corrected. `julia_live` has the same shape (`get_gs_dmrg` starts from `random_state()` unless given `wf0=`) and is left as it is, recorded here as open. (2026-09-25: fixed, see `audit_2026_09_25_open_items.md`, item 4.) Pinned by `tests/test_audit_2026_09_24b_groundstate.py::test_warm_start_setters_reach_the_session`, `::test_the_callers_state_is_not_mutated`, `::test_explicit_wf0_is_read_on_a_current_chain` and `::test_transverse_ising_warm_start_seeds_the_branch`. NUMBERS CHANGE for every caller of the two setters: on v3, `set_initial_wf_guess(target)` then `gs_energy()` goes from |<target|result>|^2 = 0.4520 to 1.0000; the transverse-field Ising example at n=40 with its shipped schedule goes from Mz/n = +0.00085 to +0.4998 at B = +-0.026 and from 0.0000 to 0.3764 at B = -0.436, and its `_map` sibling from 0.2773 to 0.3701 at B = 0.436 and from 0.0000 to 0.0387 at B = 0.462. See finding 5 of `audit_2026_09_24c_hole_hunt.md`: "a stored state whose Hamiltonian is no longer the session's is solved again" held for the public correlator only, and every other reader after `set_hamiltonian(restart=False)` answered for the old Hamiltonian; `set_hamiltonian` now resets `computed_gs` whatever `restart` is.

**Where**: `Many_Body_Chain.set_initial_wf`/`set_initial_wf_guess`
(`manybodychain.py`, around line 1133), which set `self.wf0` and the pre-port
flags `gs_from_file`/`skip_dmrg_gs`, whose only reader now is `gs_energy_single`'s
`skip_dmrg`; `groundstate.py:174-175` hands a state to the session only on an
explicit `wf0=`; `:185` overwrites the Python-side guess. Dependents:
`examples/topological/transverse_ising_model/main.py` and
`transverse_ising_model_map/main.py`.

The two setters fail differently. After `set_initial_wf_guess` a sweep does run,
from the wrong start: on v2 and v3 from the session's retained previous solution,
on `"python"` from that too, or from a fresh random MPS after the Hamiltonian's
terms change (the 2026-09 record's finding 2 fix clears the session state there).
After `set_initial_wf` nothing is swept and the session returns its cached energy
and state, which is finding 11's surface. The explicit `gs_energy(wf0=wf)` is the
one working route, and only once `computed_gs` has been reset, since
`Many_Body_Chain.gs_energy` returns `self.e0` on `gs_is_current()` before it reads
`**kwargs`. `docs/user_guide.md:361` (and `user_guide.tex:433`) name
`set_initial_wf` as a way to re-enter `gs_energy()` with a wavefunction in hand;
`set_initial_wf_guess`'s docstring ("Set the initial guess, and perform the DMRG
GS calculation") is false twice; `docs/audit_2026_09_hole_hunt.md:212` states that
`gs_energy_single` re-applies a `set_initial_wf()` state through
`set_wavefunction`, which the code did not do when that was written; and
`tests/test_bond_dimension_ramp.py`'s module docstring names "an MPS handed in by
set_initial_wf", while its test relies only on the session keeping its previous
state. Nothing in `tests/` goes through the documented route. `julia_live` has
the same shape by reading (`get_gs_dmrg` starts from `random_state()` unless given
`wf0=`).

**Expected**: the next solve starts from the guess (`set_initial_wf_guess`) or
returns it unswept (`set_initial_wf`).

Evidence, from the reviewer who owned this candidate, a spy on the session logging
every `set_wavefunction` call and the overlap of the start state with the target
just before each solve (v3 shown; `"python"` and v2 in `_python.out` and `_v2.out`
show the same pattern):

`reviews/realtime_RC2/01_warm_start_reaches_session.py`:

```python
"""R-C2 refutation probe 1: does set_initial_wf_guess(wf) / set_initial_wf(wf)
reach the DMRG session as the start state of the next solve?

3-site open S=1/2 Heisenberg chain, ground doublet at E=-1 (Sz_tot=+-1/2).
target = the pure doublet member the solved state has LESS weight in, built
as (Sz_tot -+ 1/2)|s>, normalized: an exact eigenstate, so a sweep that
really starts from it stays on it (<Sz_tot> = +-1/2).

Spy: every session.set_wavefunction call is logged, and immediately before
each session.gs_energy call the overlap |<target|session start>|^2 is read
off session.gs_wavefunction() (safe: only after the first solve, so the
session always has a wf0 and v2/v3 cannot abort).
Backend from argv[1]: python | 2 | 3."""
import sys
import numpy as np
np.random.seed(1)
import dmrgpy
from dmrgpy import spinchain, mps
print("dmrgpy from", dmrgpy.__file__)

v = sys.argv[1]
v = v if v == "python" else int(v)
n = 3
sd = spinchain.Spin_Chain([2]*n, itensor_version=v)
h = 0
for i in range(n-1):
    h = h + sd.Sx[i]*sd.Sx[i+1] + sd.Sy[i]*sd.Sy[i+1] + sd.Sz[i]*sd.Sz[i+1]
sd.set_hamiltonian(h)
sd.maxm, sd.nsweeps, sd.noise, sd.cutoff = 10, 10, 1e-7, 1e-12
Szd = sd.Sz[0] + sd.Sz[1] + sd.Sz[2]

s0 = sd.get_gs().copy()
sz0 = sd.vev(Szd).real
print("itensor_version=%r: solved E = %.12f, <Sz_tot> = %+.4f" % (v, sd.e0, sz0))
shift = -0.5 if sz0 > 0 else +0.5   # pick the member the solved state is far from
target = (Szd + shift)*s0
target = target*(1/np.sqrt(target.dot(target).real))
sz_t = +shift  # (Sz_tot+shift) annihilates the member at -shift
print("target member: <Sz_tot> = %+.4f, E = %.12f, |<target|solved>|^2 = %.4f"
      % (sd.vev(Szd, wf=target).real, sd.vev(h, wf=target).real,
         abs(target.dot(s0))**2))

log = []
real = sd._session


def w(handle):
    return abs(target.dot(mps.MPS(MBO=sd, cpp_handle=handle)))**2


class Spy:
    def __getattr__(self, name):
        attr = getattr(real, name)
        if name == "set_wavefunction":
            def f(hd):
                log.append("  session.set_wavefunction(|<target|wf>|^2=%.4f)" % w(hd))
                return attr(hd)
            return f
        if name == "gs_energy":
            def g(*a, **k):
                log.append("  session.gs_energy(%s) starts from |<target|start>|^2=%.4f"
                           % (k, w(real.gs_wavefunction())))
                return attr(*a, **k)
            return g
        return attr


spy = Spy()
sd._session = spy
sd._session_ham_cache = (spy, sd._session_ham_cache[1])  # keep the send-cache hit


def report(label):
    wf = sd.get_gs()
    print("%s: E = %.12f, <Sz_tot> = %+.4f, |<target|result>|^2 = %.4f"
          % (label, sd.e0, sd.vev(Szd).real, abs(target.dot(wf))**2))
    for line in log: print(line)
    log.clear()


sd.set_initial_wf_guess(target)
sd.gs_energy()
report("(a) set_initial_wf_guess(target); gs_energy()")

sd.set_initial_wf(target)
sd.get_gs()
report("(b) set_initial_wf(target); get_gs()")

sd.set_initial_wf(target)
sd.gs_energy()
report("(c) set_initial_wf(target); gs_energy()")

sd.computed_gs = False
sd.gs_energy(wf0=target.copy())
report("(d) gs_energy(wf0=target), the explicit argument")
print("expected if the guess were the start: <Sz_tot> = %+.4f" % sz_t)
```

`reviews/realtime_RC2/01_warm_start_reaches_session_v3.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
itensor_version=3: solved E = -1.000000000000, <Sz_tot> = -0.4074
target member: <Sz_tot> = +0.5000, E = -1.000000000000, |<target|solved>|^2 = 0.0926
(a) set_initial_wf_guess(target); gs_energy(): E = -1.000000000000, <Sz_tot> = -0.4074, |<target|result>|^2 = 0.0926
  session.gs_energy({'skip_dmrg': False}) starts from |<target|start>|^2=0.0926
(b) set_initial_wf(target); get_gs(): E = -1.000000000000, <Sz_tot> = -0.4074, |<target|result>|^2 = 0.0926
  session.gs_energy({'skip_dmrg': True}) starts from |<target|start>|^2=0.0926
(c) set_initial_wf(target); gs_energy(): E = -1.000000000000, <Sz_tot> = -0.4074, |<target|result>|^2 = 0.0926
  session.gs_energy({'skip_dmrg': True}) starts from |<target|start>|^2=0.0926
(d) gs_energy(wf0=target), the explicit argument: E = -1.000000000000, <Sz_tot> = +0.5000, |<target|result>|^2 = 1.0000
  session.set_wavefunction(|<target|wf>|^2=1.0000)
  session.gs_energy({'skip_dmrg': True}) starts from |<target|start>|^2=1.0000
expected if the guess were the start: <Sz_tot> = +0.5000
```

and the example itself, reduced to a script that runs the same Hamiltonian and
field loop three ways (as written, with the explicit `wf0=`, and with no guess),
at the example's own n=40, maxm=30, nsweeps=15:

`reviews/realtime_RC2/02_tfim_reduced.py`:

```python
"""R-C2 probe 2: does the missing hand-off change what
examples/topological/transverse_ising_model{,_map}/main.py compute?

Same Hamiltonian and loop shape as the examples (one chain, set_hamiltonian
per field, warm-start guess = the -Mz ferromagnet solved first), three
variants each on its own chain:
  as_written : set_initial_wf_guess(wffe) then vev(Mz)   (the example)
  explicit   : gs_energy(wf0=wffe.copy()) then vev(Mz)    (what it means)
  no_guess   : nothing, then vev(Mz)
Schedule pinned explicitly (not inherited).
argv: backend bmin bmax nb n maxm nsweeps"""
import sys
import numpy as np
np.random.seed(1)
import dmrgpy
from dmrgpy import spinchain
print("dmrgpy from", dmrgpy.__file__)

v = sys.argv[1]; v = v if v == "python" else int(v)
bmin, bmax, nb, n, maxm, nsweeps = (float(sys.argv[2]), float(sys.argv[3]),
        int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]))
bs = np.linspace(bmin, bmax, nb)
print("backend=%r n=%d maxm=%d nsweeps=%d noise=1e-7 cutoff=1e-12 ramp=on(10,0.5,0.1)"
      % (v, n, maxm, nsweeps))


def run(variant):
    np.random.seed(1)  # same RNG stream per variant: a no-op guess must give bit-identical output
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    sc.maxm, sc.nsweeps, sc.noise, sc.cutoff = maxm, nsweeps, 1e-7, 1e-12
    sc.bond_ramp, sc.bond_ramp_start = True, 10
    sc.bond_ramp_fraction, sc.bond_ramp_noise_decay = 0.5, 0.1

    def geth(b):
        h = 0
        for i in range(n-1): h = h - sc.Sz[i]*sc.Sz[i+1]
        for i in range(n): h = h + b*sc.Sx[i]
        return h
    Mz = 0
    for i in range(n): Mz = Mz + sc.Sz[i]
    sc.set_hamiltonian(-Mz); wffe = sc.get_gs().copy()
    out = []
    for b in bs:
        sc.set_hamiltonian(geth(b))
        if variant == "as_written": sc.set_initial_wf_guess(wffe)
        elif variant == "explicit": sc.gs_energy(wf0=wffe.copy())
        mz = sc.vev(Mz).real/n
        out.append((mz, sc.e0))
    return out


res = {k: run(k) for k in ("as_written", "explicit", "no_guess")}
print("%6s | %10s %10s %10s | %16s %16s %16s" % ("B", "Mz/n writ", "Mz/n expl",
      "Mz/n none", "E writ", "E expl", "E none"))
for j, b in enumerate(bs):
    a, e, z = res["as_written"][j], res["explicit"][j], res["no_guess"][j]
    print("%+6.3f | %+10.5f %+10.5f %+10.5f | %16.10f %16.10f %16.10f"
          % (b, a[0], e[0], z[0], a[1], e[1], z[1]))
A = np.array([r[0] for r in res["as_written"]])
E = np.array([r[0] for r in res["explicit"]])
Z = np.array([r[0] for r in res["no_guess"]])
eA = np.array([r[1] for r in res["as_written"]])
eE = np.array([r[1] for r in res["explicit"]])
print("max|Mz writ - Mz expl| = %.3e   max|Mz writ - Mz none| = %.3e (bit-identical: %s)"
      % (np.max(np.abs(A-E)), np.max(np.abs(A-Z)), bool(np.all(A == Z))))
print("max|E writ - E expl| = %.3e" % np.max(np.abs(eA-eE)))
```

`reviews/realtime_RC2/02_tfim_v3_full_m1to1.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
backend=3 n=40 maxm=30 nsweeps=15 noise=1e-7 cutoff=1e-12 ramp=on(10,0.5,0.1)
     B |  Mz/n writ  Mz/n expl  Mz/n none |           E writ           E expl           E none
-1.000 |   -0.00000   -0.00000   -0.00000 |   -21.2379959283   -21.2379959283   -21.2379959283
-0.949 |   -0.00000   +0.00000   -0.00000 |   -20.2816979119   -20.2816979119   -20.2816979119
-0.897 |   +0.00000   +0.00000   -0.00000 |   -19.3338297099   -19.3338297099   -19.3338297099
-0.846 |   -0.00000   +0.00000   -0.00000 |   -18.3960807356   -18.3960807356   -18.3960807356
-0.795 |   -0.00000   -0.00000   -0.00000 |   -17.4706455944   -17.4706455944   -17.4706455944
-0.744 |   +0.00000   +0.00000   -0.00000 |   -16.5604410529   -16.5604410529   -16.5604410529
-0.692 |   +0.00000   -0.00000   +0.00000 |   -15.6694565529   -15.6694565529   -15.6694565529
-0.641 |   -0.00000   -0.00000   -0.00000 |   -14.8033594197   -14.8033594197   -14.8033594197
-0.590 |   -0.00000   +0.00000   +0.00000 |   -13.9706477835   -13.9706477835   -13.9706477835
-0.538 |   +0.00000   +0.00000   +0.00000 |   -13.1852348968   -13.1852348968   -13.1852348968
-0.487 |   +0.00000   +0.00000   +0.00000 |   -12.4738059085   -12.4738059085   -12.4738059085
-0.436 |   -0.00000   +0.36524   -0.00000 |   -11.8749888311   -11.8748308715   -11.8749888311
-0.385 |   -0.00005   +0.43452   -0.00179 |   -11.3760526411   -11.3760498537   -11.3760526410
-0.333 |   -0.00006   +0.45748   -0.00190 |   -10.9554971392   -10.9554971267   -10.9554971392
-0.282 |   +0.00104   +0.47240   -0.00087 |   -10.6043824162   -10.6043824162   -10.6043824162
-0.231 |   +0.00106   +0.48277   -0.00089 |   -10.3174473036   -10.3174473036   -10.3174473036
-0.179 |   +0.00100   +0.49008   -0.00084 |   -10.0912199155   -10.0912199155   -10.0912199155
-0.128 |   +0.00101   +0.49511   -0.00085 |    -9.9233398879    -9.9233398879    -9.9233398879
-0.077 |   +0.00102   +0.49828   -0.00085 |    -9.8122270641    -9.8122270641    -9.8122270641
-0.026 |   +0.00102   +0.49981   -0.00085 |    -9.7569045426    -9.7569045426    -9.7569045426
+0.026 |   +0.00102   +0.49981   -0.00085 |    -9.7569045426    -9.7569045426    -9.7569045426
+0.077 |   +0.00102   +0.49828   -0.00085 |    -9.8122270641    -9.8122270641    -9.8122270641
+0.128 |   +0.00101   +0.49511   -0.00085 |    -9.9233398879    -9.9233398879    -9.9233398879
+0.179 |   +0.00100   +0.49008   -0.00084 |   -10.0912199155   -10.0912199155   -10.0912199155
+0.231 |   +0.00098   +0.48277   -0.00082 |   -10.3174473036   -10.3174473036   -10.3174473036
+0.282 |   +0.00091   +0.47240   -0.00076 |   -10.6043824162   -10.6043824162   -10.6043824162
+0.333 |   +0.00088   +0.45748   -0.00074 |   -10.9554971392   -10.9554971267   -10.9554971392
+0.385 |   -0.00050   +0.43452   +0.00063 |   -11.3760526411   -11.3760498538   -11.3760526411
+0.436 |   -0.00000   +0.37046   +0.00000 |   -11.8749888311   -11.8748218286   -11.8749888311
+0.487 |   -0.00000   +0.00000   -0.00000 |   -12.4738059085   -12.4738059085   -12.4738059085
+0.538 |   -0.00000   -0.00000   -0.00000 |   -13.1852348968   -13.1852348968   -13.1852348968
+0.590 |   +0.00000   +0.00000   -0.00000 |   -13.9706477835   -13.9706477835   -13.9706477835
+0.641 |   -0.00000   -0.00000   -0.00000 |   -14.8033594197   -14.8033594197   -14.8033594197
+0.692 |   -0.00000   -0.00000   -0.00000 |   -15.6694565529   -15.6694565529   -15.6694565529
+0.744 |   -0.00000   +0.00000   -0.00000 |   -16.5604410529   -16.5604410529   -16.5604410529
+0.795 |   +0.00000   -0.00000   +0.00000 |   -17.4706455944   -17.4706455944   -17.4706455944
+0.846 |   -0.00000   +0.00000   +0.00000 |   -18.3960807356   -18.3960807356   -18.3960807356
+0.897 |   +0.00000   +0.00000   +0.00000 |   -19.3338297099   -19.3338297099   -19.3338297099
+0.949 |   +0.00000   +0.00000   -0.00000 |   -20.2816979119   -20.2816979119   -20.2816979119
+1.000 |   -0.00000   +0.00000   +0.00000 |   -21.2379959283   -21.2379959283   -21.2379959283
max|Mz writ - Mz expl| = 4.988e-01   max|Mz writ - Mz none| = 1.947e-03
max|E writ - E expl| = 1.670e-04
```

**Reviewer (CONFIRMED, NARROWED)**: v2 prints the same (`02_tfim_v2_full_m1to1.out`,
0.4988); on `"python"` at n=16 the sign in the ordered phase is a per-field coin
flip, the two branches being near-degenerate for |b| below about 0.46
(`02_tfim_python_n16_m1to1_reseeded.out`, and `03_rng_trace.py`, which shows the
as-written and no-guess runs take identical code paths); `transverse_ising_model_map`
is right by accident on v3, starting at b=0 from the all-up state the previous
solve retained, and off by at most 0.083 near the transition
(`02_tfim_v3_full_0to1.out`). The explicit route's energy is up to 1.7e-4 higher
near |b| = 0.44, as it must be: the symmetry-broken branch sits above the finite
chain's symmetric state near the transition, so the as-written numbers are the
symmetric state's zero, not a wrong energy, and what changes is the branch the
example explicitly seeds. Struck: "the user guide documents
`set_initial_wf_guess`" (it names `set_initial_wf`; both are affected); "returns
whatever the session's own random start converges to" (it sweeps from the retained
state); and `best_gs` leaving the best of n solves on the Python side only, which
holds at 6e-9 and becomes a separate lead. Kept as its own finding rather than
folded into finding 11: the root line is shared, the broken contract (a solve that
does not start where it was told) and its consequence are not.

**Suggested fix**: finding 11's. Pushing `self.wf0` on every `gs_energy_single`
would break the cache hit `dynamics.py:190` relies on and alias the caller's MPS on
`"python"` (`gs_energy(wf0=x)` already mutates `x` there,
`reviews/realtime_RC2/05_explicit_route_aliasing.py`), so the push must be of a copy
and only of a state the caller injected. With it both examples become right without
editing them. Correct `set_initial_wf_guess`'s docstring and the 2026-09 record's
line 212. NUMBERS CHANGE for every caller of the two setters, the two examples
included.

### 13. `get_kondo_spectrum(mode="DMRG", n_gs>1)` runs a hidden DMRG sweep from each member of the manifold before measuring it, so on v2 and v3, at a split below `delta`, the excited member relaxes into the lower one in 16 of 24 runs and the call silently returns the single-state 2.0 against the two-state average 1.5, 33 per cent high

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `kondo` &middot; origin `30200a4` (the previous record's finding 12 fix)

**Status**: FIXED, with finding 11's change. Each member is installed by `set_gs` alone, which the next ground-state read hands to the session unswept, so no member is re-swept; the `hasattr(session, "set_wavefunction")` guard stays and v2 is accepted on purpose, since it behaves like v3 in every run, and the docstring and the comment at `spinchain.py:412-413` are rewritten. Pinned by `tests/test_audit_2026_09_24b_groundstate.py::test_n_gs_measures_each_member_unswept_and_restores_the_chain` on v3 and v2, six runs each, since the defect showed in 16 of 24 runs and a single run is no guard. NUMBERS CHANGE on v2 and v3 at a split below `delta`: at eps=1e-5, delta=1e-4, six runs on v3 went from 1 of 6 at the two-state average 1.5 to 6 of 6 (1.499988 to 1.500009).

**Where**: `src/dmrgpy/spinchain.py`, the `n_gs` loop of
`_get_kondo_spectrum_dmrg`, which installs each member with `set_gs` and
`session.set_wavefunction`; `set_wavefunction` drops the session's cached energy
(`pyitensor/chain.py:614`, `mpscpp3/chain_session.h:619`,
`mpscpp2/chain_session.h:194`); `dynamics.py:190`'s `set_initial_wf(self.wf0)`
sets `computed_gs=False` (`manybodychain.py:1133`), so `gs_is_current` returns
False before it reads the solver key (`groundstate.py:78`), and `kpmdmrg`'s
`get_gs()` reaches `gs_energy_single`, whose `session.gs_energy(skip_dmrg=True)`
has no cached energy and so runs `dmrg()` from the member (`groundstate.py:178`),
the result replacing `wf0` (`:185`).

The fix's own comment is the false premise, in writing: the dropped energy is
something "nothing reads again while the solver key is current". It is finding
11's round trip seen from the other side: there the cache is intact and the
injected state is silently replaced by the solved one; here `set_wavefunction` has
dropped the cache, so the same round trip runs a real sweep from the member. The
split warning is computed from the members' energies before the loop, so it cannot
see a relaxation, and the regime is exactly the one `n_gs` was added for, a split
below what the sweep resolved. The chain is the previous record's finding 12
crossing chain: an S=1 impurity with D*Sz0^2 + (D+eps)*Sz0 (D=1e-3) next to two
field-polarized S=1/2 spectators, eps below delta so the warning stays silent. It
survived because `test_n_gs_averages_the_degenerate_manifold_like_ed` runs the
exact crossing on `"python"`, where the members survive, and
`test_n_gs_1_never_touches_the_excited_states` stubs the path out.

**Expected**: each member measured as returned, the ED two-state average 1.5.

Repro, from the hunter, the shipped call (A) against the same call with the
re-solve suppressed (B), v3 at eps=1e-5, delta=1e-4:

```bash
cd <scratch>/kondo && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 05_ngs_v3_counterfactual.py
```

`kondo/05_ngs_v3_counterfactual.py`:

```python
# Counterfactual for 04: the same v3 near-crossing n_gs=2 call, (A) as
# shipped, with every groundstate.gs_energy_single call inside the n_gs loop
# logged, and (B) with that re-solve suppressed -- gs_energy_single replaced,
# for the duration of the get_kondo_spectrum call only, by one that marks the
# chain's current wf0 as its ground state with energy <wf0|H|wf0> and runs no
# sweep. If (B) gives the ED two-state average in every run while (A) does
# not, the re-solve is the cause.
import warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain, groundstate
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk import conductance
print("dmrgpy from", dmrgpy.__file__)
D, TP, eps = 1e-3, 2*np.pi, 1e-5
def chain():
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=3)
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + (D+eps)*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc
def sz0(sc, wf): return float(np.real(wf.aMb(sc.Sz[0], wf)/wf.dot(wf)))
eVs = np.array([-0.5e-3, 0.5e-3])
es = np.linspace(-30e-3, 30e-3, 3001)
kw = dict(site=0, T=0.0, order=2, mode="DMRG", submode="KPM", delta=1e-4, es=es)
sc0 = chain()
ks = KondoSpectrum(sc0, 0, T=0.0)
ks.p = np.zeros(ks.dim); ks.p[:2] = 0.5
print("ED two-state average /2pi =", np.round(conductance.second_order_dIdV(ks, eVs)/TP, 4))

real = groundstate.gs_energy_single
calls = []
def logged(self, *a, **k):
    b = sz0(self, self.wf0)
    out = real(self, *a, **k)
    calls.append((b, sz0(self, self.wf0)))
    return out
def no_sweep(self, *a, **k):
    wf = self.wf0
    self.e0 = float(np.real(wf.aMb(self.hamiltonian, wf)/wf.dot(wf)))
    self.computed_gs = True
    self._gs_solver_key = groundstate.solver_key(self)
    calls.append((sz0(self, wf), sz0(self, wf)))
    return self.e0
for label, repl in (("(A) shipped", logged), ("(B) re-solve suppressed", no_sweep)):
    vals = []
    for run in range(8):
        sc = chain()
        sc.gs_energy()           # solve first, with the real solver
        calls.clear()
        groundstate.gs_energy_single = repl
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _, d2 = sc.get_kondo_spectrum(eVs, n_gs=2, **kw)
        finally:
            groundstate.gs_energy_single = real
        vals.append(d2[1]/TP)
        print("%s run %d: n_gs=2 /2pi = %s   gs_energy_single calls in the loop (<Sz0> in -> out): %s"
              % (label, run, np.round(d2/TP, 4),
                 ", ".join("%+.3f->%+.3f" % c for c in calls)))
    vals = np.array(vals)
    print("%s: %d of 8 runs at the two-state average 1.5 (|x-1.5|<1e-3), values %s"
          % (label, np.sum(np.abs(vals-1.5) < 1e-3), np.round(vals, 4)))
```

`kondo/05_ngs_v3_counterfactual.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
ED two-state average /2pi = [1.5 1.5]
(A) shipped run 0: n_gs=2 /2pi = [2. 2.]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->-1.000, -1.000->-1.000, -1.000->-1.000
(A) shipped run 1: n_gs=2 /2pi = [2. 2.]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->-1.000, -1.000->-1.000, -1.000->-1.000
(A) shipped run 2: n_gs=2 /2pi = [2. 2.]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->-1.000, -1.000->-1.000, -1.000->-1.000
(A) shipped run 3: n_gs=2 /2pi = [2. 2.]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->-1.000, -1.000->-1.000, -1.000->-1.000
(A) shipped run 4: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(A) shipped run 5: n_gs=2 /2pi = [2. 2.]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->-1.000, -1.000->-1.000, -1.000->-1.000
(A) shipped run 6: n_gs=2 /2pi = [2. 2.]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->-1.000, -1.000->-1.000, -1.000->-1.000
(A) shipped run 7: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->-0.000, -0.000->-0.000, -0.000->-0.000
(A) shipped: 2 of 8 runs at the two-state average 1.5 (|x-1.5|<1e-3), values [2.  2.  2.  2.  1.5 2.  2.  1.5]
(B) re-solve suppressed run 0: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(B) re-solve suppressed run 1: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(B) re-solve suppressed run 2: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(B) re-solve suppressed run 3: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(B) re-solve suppressed run 4: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(B) re-solve suppressed run 5: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(B) re-solve suppressed run 6: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(B) re-solve suppressed run 7: n_gs=2 /2pi = [1.5 1.5]   gs_energy_single calls in the loop (<Sz0> in -> out): -1.000->-1.000, -1.000->-1.000, -1.000->-1.000, +0.000->+0.000, +0.000->+0.000, +0.000->+0.000
(B) re-solve suppressed: 8 of 8 runs at the two-state average 1.5 (|x-1.5|<1e-3), values [1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced on v3 and v2 (3 of 6 relaxed on v2
at eps=1e-5, `reviews/kondo_K1/02_hunter11_v2.py`); survives a pinned schedule
(nsweeps=40, no bond ramp, no noise, 3 of 8), which a caller cannot pin their way
out of, since v3 exposes no local-solver tolerance; and what decides relaxation is
the member's own weight on the lower level, not the first solve's convergence
(3.7e-10, 1.9e-9 and 1.8e-8 relaxed fully in one hidden sweep, 2.4e-11 and below
survived):

`reviews/kondo_K1/03_v3_pinned_schedule.py`:

```python
# Reviewer probe for K1: the hunter's 04 on itensor_version=3, but with an
# explicitly pinned, longer schedule (nsweeps=40, bond_ramp off, noise 0),
# and with the Sz0=-1 weight of each member logged, since the conjecture is
# that a member relaxes only when it carries a component along the lower
# member (a two-dimensional near-degenerate subspace is resolved exactly by
# one Krylov step, whatever the split). Per run: first-solve E-E0, and for
# each member entering/leaving the first correlator term: <Sz0>, P(Sz0=-1),
# E-E0. Anchor: full ED (KondoSpectrum with p=[0.5,0.5]).
import sys, warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk import conductance
import dmrgpy.kondospectrumtk.secondorder_dc as sdc
print("dmrgpy from", dmrgpy.__file__)
D, TP = 1e-3, 2*np.pi
eps = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-5
NSW = int(sys.argv[2]) if len(sys.argv) > 2 else 40
RAMP = (sys.argv[3] == "ramp") if len(sys.argv) > 3 else False
NOISE = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
print("eps=%g nsweeps=%d bond_ramp=%s noise=%g maxm=30" % (eps, NSW, RAMP, NOISE))
def chain():
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=3)
    sc.maxm, sc.nsweeps = 30, NSW
    sc.bond_ramp = RAMP
    sc.noise = NOISE
    h = D*sc.Sz[0]*sc.Sz[0] + (D+eps)*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc
def ev(sc, op, wf): return float(np.real(wf.aMb(op, wf)/wf.dot(wf)))
log = []
_so = sdc.second_order_dIdV_dc
def logging_so(chain, *a, **k):
    Pm = 0.5*chain.Sz[0]*chain.Sz[0] - 0.5*chain.Sz[0]
    b = (ev(chain, chain.Sz[0], chain.wf0), ev(chain, Pm, chain.wf0), ev(chain, chain.hamiltonian, chain.wf0))
    out = _so(chain, *a, **k)
    a2 = (ev(chain, chain.Sz[0], chain.wf0), ev(chain, Pm, chain.wf0), ev(chain, chain.hamiltonian, chain.wf0))
    log.append((b, a2))
    return out
sdc.second_order_dIdV_dc = logging_so
eVs = np.array([-0.5e-3, 0.5e-3])
es = np.linspace(-30e-3, 30e-3, 3001)
kw = dict(site=0, T=0.0, order=2, mode="DMRG", submode="KPM", delta=1e-4, es=es)
sc0 = chain()
E0 = sc0.get_ED_obj().gs_energy()
ks = KondoSpectrum(sc0, 0, T=0.0)
ks.p = np.zeros(ks.dim); ks.p[:2] = 0.5
avg = conductance.second_order_dIdV(ks, eVs)
print("ED: E1-E0=%.2e  two-state average /2pi=%s" % (ks.e[1], np.round(avg/TP, 4)))
vals = []
for run in range(8):
    sc = chain()
    eb = sc.gs_energy()
    log.clear()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _, d2 = sc.get_kondo_spectrum(eVs, n_gs=2, **kw)
    nw = len([w for w in caught if "n_gs" in str(w.message)])
    vals.append(d2[1]/TP)
    s = "run %d: first solve E-E0 %.2e | n_gs=2 /2pi=%s warn=%d |" % (run, eb-E0, np.round(d2/TP, 4), nw)
    for k, ((sb, pb, e1), (sa, pa, e2)) in enumerate(log):
        s += " m%d <Sz0> %+.4f->%+.4f P(-1) %.1e->%.1e E-E0 %.1e->%.1e;" % (k, sb, sa, pb, pa, e1-E0, e2-E0)
    print(s)
vals = np.array(vals)
print("%d of 8 runs at the two-state average 1.5 (|x-1.5|<1e-3), values %s"
      % (np.sum(np.abs(vals-1.5) < 1e-3), np.round(vals, 4)))
```

`reviews/kondo_K1/03_v3_pinned_schedule.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
eps=1e-05 nsweeps=40 bond_ramp=False noise=0 maxm=30
ED: E1-E0=1.00e-05  two-state average /2pi=[1.5 1.5]
run 0: first solve E-E0 5.04e-07 | n_gs=2 /2pi=[2. 2.] warn=0 | m0 <Sz0> -1.0000->-1.0000 P(-1) 1.0e+00->1.0e+00 E-E0 8.3e-10->0.0e+00; m1 <Sz0> +0.0000->-1.0000 P(-1) 3.7e-10->1.0e+00 E-E0 1.0e-05->0.0e+00;
run 1: first solve E-E0 8.33e-08 | n_gs=2 /2pi=[1.5 1.5] warn=0 | m0 <Sz0> -1.0000->-1.0000 P(-1) 1.0e+00->1.0e+00 E-E0 1.5e-10->1.7e-17; m1 <Sz0> +0.0000->-0.0000 P(-1) 2.4e-11->2.6e-12 E-E0 1.0e-05->1.0e-05;
run 2: first solve E-E0 2.67e-09 | n_gs=2 /2pi=[1.5 1.5] warn=0 | m0 <Sz0> -1.0000->-1.0000 P(-1) 1.0e+00->1.0e+00 E-E0 1.6e-09->0.0e+00; m1 <Sz0> -0.0000->-0.0000 P(-1) 9.8e-15->1.3e-14 E-E0 1.0e-05->1.0e-05;
run 3: first solve E-E0 5.97e-07 | n_gs=2 /2pi=[2. 2.] warn=0 | m0 <Sz0> -1.0000->-1.0000 P(-1) 1.0e+00->1.0e+00 E-E0 1.3e-09->0.0e+00; m1 <Sz0> +0.0000->-1.0000 P(-1) 1.9e-09->1.0e+00 E-E0 1.0e-05->0.0e+00;
run 4: first solve E-E0 9.57e-06 | n_gs=2 /2pi=[2. 2.] warn=0 | m0 <Sz0> -1.0000->-1.0000 P(-1) 1.0e+00->1.0e+00 E-E0 6.9e-11->2.4e-17; m1 <Sz0> +0.0000->-1.0000 P(-1) 1.8e-08->1.0e+00 E-E0 1.0e-05->0.0e+00;
run 5: first solve E-E0 0.00e+00 | n_gs=2 /2pi=[1.5 1.5] warn=0 | m0 <Sz0> -1.0000->-1.0000 P(-1) 1.0e+00->1.0e+00 E-E0 0.0e+00->0.0e+00; m1 <Sz0> +0.0000->-0.0000 P(-1) 1.8e-42->2.6e-28 E-E0 1.0e-05->1.0e-05;
run 6: first solve E-E0 0.00e+00 | n_gs=2 /2pi=[1.5 1.5] warn=0 | m0 <Sz0> -1.0000->-1.0000 P(-1) 1.0e+00->1.0e+00 E-E0 0.0e+00->0.0e+00; m1 <Sz0> +0.0000->-0.0000 P(-1) 4.9e-44->9.0e-26 E-E0 1.0e-05->1.0e-05;
run 7: first solve E-E0 3.54e-16 | n_gs=2 /2pi=[1.5 1.5] warn=0 | m0 <Sz0> -1.0000->-1.0000 P(-1) 1.0e+00->1.0e+00 E-E0 0.0e+00->0.0e+00; m1 <Sz0> +0.0000->-0.0000 P(-1) 6.5e-36->2.4e-25 E-E0 1.0e-05->1.0e-05;
5 of 8 runs at the two-state average 1.5 (|x-1.5|<1e-3), values [2.  1.5 1.5 2.  2.  1.5 1.5 1.5]
```

Not observed at eps=1e-6 (0 of 11, `04_v3_eps1e-6.py`), so the default `delta=2e-6`
admits no split where it was seen; not at the exact crossing; not on `"python"`.
Struck: "11 of 24" (the hunter's batches sum to 16 of 24). Confirmed a sub-claim
that is an error in the previous record rather than a separate defect: finding
12's Status says `n_gs>1` "refuses `itensor_version=2` ... (not run)", but the
guard tests `hasattr(session, "set_wavefunction")`, which `mpscpp2/bindings.cc:107`
exports, so v2 runs it and behaves like v3 in every run. Kept as a separate entry
from finding 11 with one shared fix: finding 11's naive fix, "have `set_gs` push to
the session", is exactly what finding 12 did by hand, and this is its outcome.

**Suggested fix**: finding 11's, at the round trip. The reviewer's in-memory patch
of line 190 gives 1.5 in 8 of 8 runs with no sweep and the chain back in its
pre-call state, which also closes finding 14; a bare deletion is not enough, since
a clone of a solved chain keeps `computed_gs=True` with an empty session and then
crashes:

`reviews/kondo_K1/05_v3_no_roundtrip.py`:

```python
# Reviewer counterfactual (C) for K1's fix shape: the same v3 near-crossing
# n_gs=2 call as the hunter's 05, with dynamics.get_dynamical_correlator's
# first line, self.set_initial_wf(self.wf0), deleted (source-patched in
# memory, the repository is not touched), and every gs_energy_single call
# during the n_gs=2 call counted. If this gives the ED average 1.5 in every
# run with zero calls, removing the round trip closes K1 on its own.
# Then, on "python" (which raises rather than aborts), what the deleted line
# was also doing: a clone of a solved chain keeps computed_gs=True and a
# wf0 while its session is fresh and empty.
import inspect, warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain, groundstate, dynamics
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk import conductance
print("dmrgpy from", dmrgpy.__file__)
LINE = "self.set_initial_wf(self.wf0) # set the initial wavefunction"
src = inspect.getsource(dynamics.get_dynamical_correlator)
assert src.count(LINE) == 1
orig = dynamics.get_dynamical_correlator
ns = dict(dynamics.__dict__)
exec(compile(src.replace(LINE, "pass # reviewer counterfactual"), dynamics.__file__, "exec"), ns)
patched = ns["get_dynamical_correlator"]
D, TP, eps = 1e-3, 2*np.pi, 1e-5
def chain(v=3):
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=v)
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + (D+eps)*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc
def sz0(sc, wf): return float(np.real(wf.aMb(sc.Sz[0], wf)/wf.dot(wf)))
eVs = np.array([-0.5e-3, 0.5e-3])
es = np.linspace(-30e-3, 30e-3, 3001)
kw = dict(site=0, T=0.0, order=2, mode="DMRG", submode="KPM", delta=1e-4, es=es)
sc0 = chain()
E0 = sc0.get_ED_obj().gs_energy()
ks = KondoSpectrum(sc0, 0, T=0.0)
ks.p = np.zeros(ks.dim); ks.p[:2] = 0.5
print("ED two-state average /2pi =", np.round(conductance.second_order_dIdV(ks, eVs)/TP, 4))
real = groundstate.gs_energy_single
ncalls = [0]
def counting(self, *a, **k):
    ncalls[0] += 1
    return real(self, *a, **k)
vals = []
for run in range(8):
    sc = chain()
    eb = sc.gs_energy()
    ncalls[0] = 0
    groundstate.gs_energy_single = counting
    dynamics.get_dynamical_correlator = patched
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, d2 = sc.get_kondo_spectrum(eVs, n_gs=2, **kw)
    finally:
        groundstate.gs_energy_single = real
        dynamics.get_dynamical_correlator = orig
    after = sc.gs_energy()
    vals.append(d2[1]/TP)
    print("(C) run %d: first solve E-E0 %.2e | n_gs=2 /2pi = %s | gs_energy_single calls during the call: %d | gs_energy()-E0 after %.2e, <Sz0> of get_gs() after %+.4f"
          % (run, eb-E0, np.round(d2/TP, 4), ncalls[0], after-E0, sz0(sc, sc.get_gs())))
vals = np.array(vals)
print("(C) round trip deleted: %d of 8 runs at 1.5 (|x-1.5|<1e-3), values %s"
      % (np.sum(np.abs(vals-1.5) < 1e-3), np.round(vals, 4)))

# what else the deleted line does: a clone of a solved chain
np.random.seed(0)
sp = chain("python")
sp.gs_energy()
for label, fn in (("shipped", orig), ("round trip deleted", patched)):
    cl = sp.clone()
    print("clone (%s): computed_gs=%s, wf0 is None: %s, session has wf0: %s"
          % (label, cl.computed_gs, cl.wf0 is None, getattr(cl._session, "wf0", "n/a") is not None))
    dynamics.get_dynamical_correlator = fn
    try:
        x, y = cl.get_dynamical_correlator(name=(cl.Sz[0], cl.Sz[0]), es=np.linspace(-0.03, 0.03, 7), delta=1e-3)
        print("   clone (%s): correlator ran, max|y| = %.4g" % (label, np.max(np.abs(y))))
    except Exception as e:
        print("   clone (%s): %s: %s" % (label, type(e).__name__, str(e)[:200]))
    finally:
        dynamics.get_dynamical_correlator = orig
```

`reviews/kondo_K1/05_v3_no_roundtrip.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
ED two-state average /2pi = [1.5 1.5]
(C) run 0: first solve E-E0 8.51e-06 | n_gs=2 /2pi = [1.5 1.5] | gs_energy_single calls during the call: 0 | gs_energy()-E0 after 8.51e-06, <Sz0> of get_gs() after -0.1497
(C) run 1: first solve E-E0 6.85e-06 | n_gs=2 /2pi = [1.5 1.5] | gs_energy_single calls during the call: 0 | gs_energy()-E0 after 6.85e-06, <Sz0> of get_gs() after -0.3148
(C) run 2: first solve E-E0 1.96e-16 | n_gs=2 /2pi = [1.5 1.5] | gs_energy_single calls during the call: 0 | gs_energy()-E0 after 1.96e-16, <Sz0> of get_gs() after -1.0000
(C) run 3: first solve E-E0 0.00e+00 | n_gs=2 /2pi = [1.5 1.5] | gs_energy_single calls during the call: 0 | gs_energy()-E0 after 0.00e+00, <Sz0> of get_gs() after -1.0000
(C) run 4: first solve E-E0 9.61e-06 | n_gs=2 /2pi = [1.5 1.5] | gs_energy_single calls during the call: 0 | gs_energy()-E0 after 9.61e-06, <Sz0> of get_gs() after -0.0396
(C) run 5: first solve E-E0 8.97e-06 | n_gs=2 /2pi = [1.5 1.5] | gs_energy_single calls during the call: 0 | gs_energy()-E0 after 8.97e-06, <Sz0> of get_gs() after -0.1030
(C) run 6: first solve E-E0 3.12e-06 | n_gs=2 /2pi = [1.5 1.5] | gs_energy_single calls during the call: 0 | gs_energy()-E0 after 3.12e-06, <Sz0> of get_gs() after -0.6891
(C) run 7: first solve E-E0 4.00e-06 | n_gs=2 /2pi = [1.5 1.5] | gs_energy_single calls during the call: 0 | gs_energy()-E0 after 4.00e-06, <Sz0> of get_gs() after -0.6001
(C) round trip deleted: 8 of 8 runs at 1.5 (|x-1.5|<1e-3), values [1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5]
clone (shipped): computed_gs=True, wf0 is None: False, session has wf0: False
   clone (shipped): correlator ran, max|y| = 692.6
clone (round trip deleted): computed_gs=True, wf0 is None: False, session has wf0: False
   clone (round trip deleted): RuntimeError: Chain.kpm_dynamical_correlator called before set_hamiltonian
```

The hunter's alternative, an energy argument to `set_wavefunction`, closes this
finding only, needs a rebuild and leaves the question of which energy the session
then reports. Accept v2 on purpose, since it behaves like v3, and correct the
previous record's Status line and the docstring. NUMBERS CHANGE only on v2/v3 at a
split below `delta`, where the relaxed runs move from 2.0 onto 1.5.

### 14. After `get_kondo_spectrum(mode="DMRG", n_gs>1)`, `gs_energy()` returns the last member's energy while `get_gs()` is the restored ground state, because the `finally` never restores `self.e0`; a consumer that shifts H by `gs_energy()` then picks up a phase (E_last - E0)*t, -14.36 against an exact -9.39 at zero bias for the two-time Kondo term

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `kondo` &middot; origin `30200a4`

**Status**: FIXED. The `n_gs` loop snapshots `(wf0, e0, computed_gs, _gs_solver_key, _gs_injected)` before the loop and restores it as a unit in the `finally`, pushing a copy of the solved state back to the session, and inside the loop `gs_energy()` is each member's own energy, the convention of the ED references. Pinned by `tests/test_audit_2026_09_24b_groundstate.py::test_n_gs_measures_each_member_unswept_and_restores_the_chain` (`gs_energy()` equal before and after, at a split below delta) and `::test_n_gs_gs_energy_is_each_members_own_inside_the_loop`. NUMBERS CHANGE only for a caller that read `gs_energy()` between an `n_gs>1` call and the next correlator: on `"python"` at eps=1e-5 `gs_energy()` after minus before goes from 1.000e-05 to exactly 0, and on v3 it is 0 in six runs.

**Where**: `src/dmrgpy/spinchain.py`, the `finally` of the `n_gs` loop, which
restores `wf0`/`computed_gs` through `set_gs` and the session through
`set_wavefunction` but not `self.e0`; inside the loop `gs_energy_single`
(`groundstate.py`, `self.e0 = out`) overwrites it on every correlator call, and
`gs_is_current` then returns the stale value while the solver key holds. The
amplifying consumer is `kondospectrumtk/dmrgtwotime.py:142`, `E0 =
chain.gs_energy()`.

The size is E_last - E0 exactly: zero at an exact crossing (the case the feature
exists for), up to `delta` in the silent regime, the full split beyond it (where
the warning fires, about the spectrum, not about the chain afterwards). The next
`get_dynamical_correlator` call heals it, and so does `get_kondo_spectrum` itself,
which together with a regression test that uses the exact crossing and checks
`vev` and a repeated KPM call but never `gs_energy()` is why nothing caught it.
The shipped comment ("gs_energy() then returns the Python-side e0") and the
previous record's Status ("the chain's own state is restored exactly") are the
false premise.

**Expected**: `gs_energy()` and `get_excited_states(n=1)[0][0]` equal their
pre-call values.

Repro, from the hunter (case (a) the silent near-crossing, case (b) the kondo
examples' 3-site chain; `02_ngs_stale_e0_v3.out` prints the same on v3):

```bash
cd <scratch>/kondo && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 02_ngs_stale_e0.py
```

`kondo/02_ngs_stale_e0.py`:

```python
# After get_kondo_spectrum(mode="DMRG", n_gs=g), is the chain's ground-state
# ENERGY restored along with its wavefunction? gs_energy() is compared with
# ED's E0 and with <H> of the chain's own get_gs(), before and after.
# (a) finding 12's crossing chain split by eps=1e-5 < delta=1e-4, so no warning
# (b) the kondo examples' 3-site Zeeman chain (non-degenerate, gap 0.77 meV)
#     at n_gs=2, where the split warning does fire
import sys, warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain
print("dmrgpy from", dmrgpy.__file__)
version = sys.argv[1] if len(sys.argv) > 1 else "python"
version = int(version) if version.isdigit() else version
print("itensor_version =", version)
G, MUB, D = 2.0, 5.7883818066e-5, 1e-3

def crossing(eps):
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=version)
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + (D+eps)*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc

def zeeman():
    sc = spinchain.Spin_Chain(["1/2"]*3, itensor_version=version)
    sc.maxm, sc.nsweeps = 30, 15
    h = G*MUB*10.0*sc.Sz[0]
    for i in range(2):
        h = h + 0.01*sc.SS(i, i+1)
    sc.set_hamiltonian(h)
    return sc

cases = [("(a) crossing, eps=1e-5, delta=1e-4", lambda: crossing(1e-5), 1e-4,
          np.linspace(-30e-3, 30e-3, 3001), np.array([-2e-3, 2e-3])),
         ("(b) Zeeman 3-site chain, delta=2e-5", zeeman, 2e-5,
          np.linspace(-20e-3, 20e-3, 5328), np.array([-2e-3, 2e-3]))]
for label, build, delta, es, eVs in cases:
    np.random.seed(1)
    sc = build()
    emu = np.sort(np.array(sc.get_ED_obj().get_diagonalized_hamiltonian()[0], dtype=float))
    print(label)
    print("   ED: E0 = %.12f  E1 = %.12f  E1-E0 = %.3e" % (emu[0], emu[1], emu[1]-emu[0]))
    e_before = sc.gs_energy()
    print("   before: gs_energy() = %.12f   <H> of get_gs() = %.12f"
          % (e_before, np.real(sc.vev(sc.hamiltonian))))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="DMRG",
                              submode="KPM", delta=delta, es=es, n_gs=2)
    print("   split warnings during the n_gs=2 call:",
          len([w for w in caught if "n_gs" in str(w.message)]))
    e_after = sc.gs_energy()
    print("   after:  gs_energy() = %.12f   <H> of get_gs() = %.12f"
          % (e_after, np.real(sc.vev(sc.hamiltonian))))
    print("   gs_energy() after - ED E0 = %.3e   (ED E1-E0 = %.3e)"
          % (e_after - emu[0], emu[1]-emu[0]))
    # a consumer that reads gs_energy() without going through the
    # dynamical-correlator dispatch: get_excited() energies and the gap
    print("   get_gap() after = %.6e   (ED %.6e)" % (sc.get_gap(), emu[1]-emu[0]))
    print("   gs_energy() once more = %.12f" % sc.gs_energy())
```

`kondo/02_ngs_stale_e0.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
itensor_version = python
(a) crossing, eps=1e-5, delta=1e-4
   ED: E0 = -0.009510000000  E1 = -0.009500000000  E1-E0 = 1.000e-05
   before: gs_energy() = -0.009510000000   <H> of get_gs() = -0.009510000000
   split warnings during the n_gs=2 call: 0
   after:  gs_energy() = -0.009500000000   <H> of get_gs() = -0.009510000000
   gs_energy() after - ED E0 = 1.000e-05   (ED E1-E0 = 1.000e-05)
   get_gap() after = 1.000000e-05   (ED 1.000000e-05)
   gs_energy() once more = -0.009500000000
(b) Zeeman 3-site chain, delta=2e-5
   ED: E0 = -0.010401002243  E1 = -0.009631357092  E1-E0 = 7.696e-04
   before: gs_energy() = -0.010401002243   <H> of get_gs() = -0.010401002243
   split warnings during the n_gs=2 call: 1
   after:  gs_energy() = -0.009631357092   <H> of get_gs() = -0.010401002243
   gs_energy() after - ED E0 = 7.696e-04   (ED E1-E0 = 7.696e-04)
   get_gap() after = 7.696452e-04   (ED 7.696452e-04)
   gs_energy() once more = -0.009631357092
```

**Reviewer (CONFIRMED, NARROWED)**: the consumer's whole error is the stale
`e0`: shifting `e0` by hand on a fresh chain reproduces the `n_gs` result to every
digit, and restoring it removes it to every digit:

`reviews/kondo_K2/04_consumer_attribution.py`:

```python
# Is the whole of the consumer's error (two_time_kondo_term_dmrg after an
# n_gs=2 call, hunter's 06) the stale self.e0, and nothing else the n_gs
# call leaves behind? Four runs of the same two-time term on finding 12's
# crossing chain (eps=1e-5), itensor_version="python", hunter's coarse grid:
#  (1) fresh chain
#  (2) fresh chain with self.e0 shifted by hand by +eps (no n_gs call)
#  (3) after get_kondo_spectrum(mode="DMRG", n_gs=2), as shipped
#  (4) after the same call, with self.e0 put back to its pre-call value
#      (what the suggested fix's `finally` would do)
# Anchor: ED two-time reference on the literal same (t2,tau) grid.
import warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.edtwotimeref import _levi_civita_coeff_G_chunk
from dmrgpy.kondospectrumtk.twotime import kondo_term_from_two_time
from dmrgpy.kondospectrumtk.dmrgtwotime import two_time_kondo_term_dmrg
print("dmrgpy from", dmrgpy.__file__)
D, eps = 1e-3, 1e-5
def chain():
    np.random.seed(0)
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version="python")
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + (D+eps)*sc.Sz[0] + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc
omega0, Gamma0 = 2e-3, 5e-6
n_t2_half, n_tau_half = 2, 3
dt2, dtau = 25./Gamma0/10, (2*np.pi/2e-5)/15
eVs = np.array([0.0, 1e-3, 2e-3])
t2_grid = dt2*np.arange(-n_t2_half, n_t2_half+1)
tau_grid = dtau*np.arange(-n_tau_half, n_tau_half+1)
sc = chain()
ks = KondoSpectrum(sc, 0, T=0.0)
ref = kondo_term_from_two_time(t2_grid, tau_grid,
        iter([(t2_grid, _levi_civita_coeff_G_chunk(ks, t2_grid, tau_grid))]),
        eVs, omega0, Gamma0)
tt = dict(omega0=omega0, Gamma0=Gamma0, dt2=dt2, n_t2_half=n_t2_half,
          dtau=dtau, n_tau_half=n_tau_half)
E_ed = sc.get_ED_obj().gs_energy()
print("t2 range: max|t2| = %.3e, so eps*max|t2| = %.2f rad" % (t2_grid.max(), eps*t2_grid.max()))
print("ED reference term                 :", np.round(ref, 4))
def report(label, val):
    print("%-34s: %s  max|dmrg-ED| = %.2e  (gs_energy()-E0 = %.3e)"
          % (label, np.round(val, 4), np.max(np.abs(val-ref)), sc.gs_energy()-E_ed))
e_before = sc.gs_energy()
report("(1) fresh chain", two_time_kondo_term_dmrg(sc, 0, eVs, **tt))
sc.e0 = e_before + eps
report("(2) fresh, e0 shifted by +eps", two_time_kondo_term_dmrg(sc, 0, eVs, **tt))
sc.e0 = e_before
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="DMRG", submode="KPM",
                          delta=1e-4, es=np.linspace(-30e-3, 30e-3, 3001), n_gs=2)
report("(3) after n_gs=2, shipped", two_time_kondo_term_dmrg(sc, 0, eVs, **tt))
sc.e0 = e_before
report("(4) after n_gs=2, e0 restored", two_time_kondo_term_dmrg(sc, 0, eVs, **tt))
```

`reviews/kondo_K2/04_consumer_attribution.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
t2 range: max|t2| = 1.000e+06, so eps*max|t2| = 10.00 rad
ED reference term                 : [-9.3885 -2.0773 -1.3109]
(1) fresh chain                   : [-9.3885 -2.0773 -1.3109]  max|dmrg-ED| = 1.87e-12  (gs_energy()-E0 = -2.351e-15)
(2) fresh, e0 shifted by +eps     : [-14.3581  -2.3156  -1.4612]  max|dmrg-ED| = 4.97e+00  (gs_energy()-E0 = 1.000e-05)
(3) after n_gs=2, shipped         : [-14.3581  -2.3156  -1.4612]  max|dmrg-ED| = 4.97e+00  (gs_energy()-E0 = 1.000e-05)
(4) after n_gs=2, e0 restored     : [-9.3885 -2.0773 -1.3109]  max|dmrg-ED| = 1.87e-12  (gs_energy()-E0 = -2.351e-15)
```

Its own finding, not a consequence of finding 13: it survives a finding-13 fix that
suppresses the sweep but still returns through `gs_energy_single` (variant A), on
both backends, and disappears only when the round trip is bypassed (variant B); on
`"python"` finding 13's relaxation never happens, yet this is there in full:

`reviews/kondo_K2/03_separability_from_K1.py`:

```python
# Is K2 (stale self.e0 after get_kondo_spectrum(mode="DMRG", n_gs=2)) a
# consequence of K1 (the hidden warm-start sweep inside the n_gs loop), or
# its own defect? Three variants on finding 12's crossing chain, eps=1e-5 <
# delta=1e-4 (no split warning):
#  (S) shipped code
#  (A) K1 "fixed" by suppressing the sweep but still returning through
#      gs_energy_single's slot: energy = <wf0|H|wf0>, no DMRG (the hunter's
#      05(B) counterfactual shape)
#  (B) K1 "fixed" by bypassing the re-solve entirely: the correlator
#      dispatch's set_initial_wf(self.wf0) no longer invalidates the stored
#      ground state, so gs_energy_single is never entered inside the loop
# For each: gs_energy() read inside the loop right after each member's
# second-order term (what an order=3 two-time term would read), and after
# the call gs_energy(), <H> of get_gs(), get_excited_states(n=1)'s energy,
# against ED.
import sys, warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain, groundstate, manybodychain
import dmrgpy.kondospectrumtk.secondorder_dc as sdc
print("dmrgpy from", dmrgpy.__file__)
version = sys.argv[1] if len(sys.argv) > 1 else "python"
version = int(version) if version.isdigit() else version
print("itensor_version =", version)
D, eps = 1e-3, 1e-5

def chain():
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=version)
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + (D+eps)*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc

def H(sc, wf): return float(np.real(wf.aMb(sc.hamiltonian, wf)/wf.dot(wf)))

real_single = groundstate.gs_energy_single
def no_sweep(self, *a, **k):
    wf = self.wf0
    self.e0 = H(self, wf)
    self.computed_gs = True
    self._gs_solver_key = groundstate.solver_key(self)
    return self.e0

real_siw = manybodychain.Many_Body_Chain.set_initial_wf
def keep_current(self, wf, reconverge=False):
    if wf is self.wf0 and self.computed_gs and not reconverge:
        return # the stored ground state stays current: no re-solve
    return real_siw(self, wf, reconverge=reconverge)

inloop = []
real_so = sdc.second_order_dIdV_dc
def logging_so(chain, *a, **k):
    out = real_so(chain, *a, **k)
    inloop.append((chain.gs_energy(), H(chain, chain.wf0)))
    return out
sdc.second_order_dIdV_dc = logging_so

eVs = np.array([-2e-3, 2e-3])
es = np.linspace(-30e-3, 30e-3, 3001)
for label in ("(S) shipped", "(A) sweep suppressed, e0=<wf0|H|wf0>",
              "(B) re-solve bypassed"):
    np.random.seed(1)
    sc = chain()
    emu = np.sort(np.array(sc.get_ED_obj().get_diagonalized_hamiltonian()[0], dtype=float))
    e_before = sc.gs_energy()
    inloop.clear()
    if label.startswith("(A)"): groundstate.gs_energy_single = no_sweep
    if label.startswith("(B)"): manybodychain.Many_Body_Chain.set_initial_wf = keep_current
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _, y = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="DMRG",
                                         submode="KPM", delta=1e-4, es=es, n_gs=2)
    finally:
        groundstate.gs_energy_single = real_single
        manybodychain.Many_Body_Chain.set_initial_wf = real_siw
    print(label)
    print("   ED E0 = %.12f  E1-E0 = %.3e ; before: gs_energy()-E0 = %.3e"
          % (emu[0], emu[1]-emu[0], e_before-emu[0]))
    print("   n_gs split warnings:", len([w for w in caught if "n_gs" in str(w.message)]),
          "  dIdV/2pi =", np.round(y/(2*np.pi), 4))
    for k, (eg, eh) in enumerate(inloop):
        print("   in loop, member %d: gs_energy()-E0 = %.3e  <H>(wf0)-E0 = %.3e"
              % (k, eg-emu[0], eh-emu[0]))
    ea = sc.gs_energy()
    print("   after: gs_energy()-E0 = %.3e   <H>(get_gs())-E0 = %.3e"
          % (ea-emu[0], H(sc, sc.get_gs())-emu[0]))
    es1, ws1 = sc.get_excited_states(n=1)
    print("   after: get_excited_states(n=1) energy-E0 = %.3e  <H> of its state-E0 = %.3e"
          % (es1[0]-emu[0], H(sc, ws1[0])-emu[0]))
```

`reviews/kondo_K2/03_separability_from_K1_python.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
itensor_version = python
(S) shipped
   ED E0 = -0.009510000000  E1-E0 = 1.000e-05 ; before: gs_energy()-E0 = -2.351e-15
   n_gs split warnings: 0   dIdV/2pi = [1.6923 1.6923]
   in loop, member 0: gs_energy()-E0 = -2.352e-15  <H>(wf0)-E0 = -2.352e-15
   in loop, member 1: gs_energy()-E0 = 1.000e-05  <H>(wf0)-E0 = 1.000e-05
   after: gs_energy()-E0 = 1.000e-05   <H>(get_gs())-E0 = -2.351e-15
   after: get_excited_states(n=1) energy-E0 = 1.000e-05  <H> of its state-E0 = -2.351e-15
(A) sweep suppressed, e0=<wf0|H|wf0>
   ED E0 = -0.009510000000  E1-E0 = 1.000e-05 ; before: gs_energy()-E0 = -2.356e-15
   n_gs split warnings: 0   dIdV/2pi = [1.6923 1.6923]
   in loop, member 0: gs_energy()-E0 = -2.352e-15  <H>(wf0)-E0 = -2.352e-15
   in loop, member 1: gs_energy()-E0 = 1.000e-05  <H>(wf0)-E0 = 1.000e-05
   after: gs_energy()-E0 = 1.000e-05   <H>(get_gs())-E0 = -2.349e-15
   after: get_excited_states(n=1) energy-E0 = 1.000e-05  <H> of its state-E0 = -2.349e-15
(B) re-solve bypassed
   ED E0 = -0.009510000000  E1-E0 = 1.000e-05 ; before: gs_energy()-E0 = -2.356e-15
   n_gs split warnings: 0   dIdV/2pi = [1.6923 1.6923]
   in loop, member 0: gs_energy()-E0 = -2.356e-15  <H>(wf0)-E0 = -2.352e-15
   in loop, member 1: gs_energy()-E0 = -2.356e-15  <H>(wf0)-E0 = 1.000e-05
   after: gs_energy()-E0 = -2.356e-15   <H>(get_gs())-E0 = -2.349e-15
   after: get_excited_states(n=1) energy-E0 = -2.356e-15  <H> of its state-E0 = -2.349e-15
```

Struck: the 53 per cent as a workflow number (the two-time term is internal, its
only direct in-tree caller builds fresh chains, and the public `order=3` path heals
itself first; the number is eps*max|t2| = 10 rad on a deliberately coarse grid), and
the v3 "direction flip" (not reproduced in four runs; it depends on finding 13
happening). Qualified: the 7.696e-4 is a non-degenerate chain where the warning
fires, not the silent regime.

**Suggested fix**: restore `self.e0` in the `finally`, better as a snapshot of
`(wf0, e0, computed_gs, _gs_solver_key)` restored as a unit, whichever way finding
13 is fixed; and inside the loop set `self.e0` to each member's own energy, the
convention the ED references use (`_levi_civita_coeff_G_chunk` puts no phase on the
bra, and the ED third-order sums measure from each initial state's energy), which
today holds only because of finding 13's sweep. Fixing finding 13 by bypassing the
round trip without this half would move the error inside the public `order=3,
n_gs>1` result. Pin it at a split below `delta` (eps=0 cannot see it). NUMBERS
CHANGE only for a caller that reads `gs_energy()` between an `n_gs>1` call and the
next correlator.

### 15. Under `submode="EX"`, `get_kondo_spectrum(mode="DMRG", n_gs>1)` returns the `n_gs=1` value bit for bit, because EX measures from the lowest vector of its own cached rediagonalized basis rather than from the chain's state, so at a degenerate ground state the answer is an arbitrary member anywhere in [1, 2] against the average 1.5, up to 33 per cent either side

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `kondo` &middot; origin `30200a4` (the `n_gs=` keyword)

**Status**: FIXED, both halves. `n_gs>1` accepts only the submodes on `spinchain._N_GS_SUBMODES`, each verified with the up/down test to read the chain's state (KPM, CVM, CVM_explicit, ROOTN, TD, TDZ, EX), and raises `NotImplementedError` for SECTOR and anything else; and `dcex` measures from the chain's own state projected onto its cached basis, d = C^dagger o, with transitions measured from that state's own energy d^dagger diag(E) d, raising `ValueError` when the state's weight in the basis differs from 1 by more than 1e-6. Its cache key now includes the solver parameters, since a loop over `maxm` reused the first basis for every value. Pinned by `tests/test_audit_2026_09_24b_groundstate.py::test_n_gs_averages_under_ex`, `::test_ex_at_a_unique_ground_state_is_unchanged`, `::test_ex_refuses_a_state_outside_its_basis` and `::test_n_gs_refuses_a_submode_that_does_not_read_the_state` (SECTOR, maxent). NUMBERS CHANGE for EX at a degenerate ground state, which now reads the chain's state: at the exact crossing `n_gs=2` goes from 1.1392, 1.8711 and 1.6200 on three seeds to 1.512780 on every seed (the 0.013 above 1.5 is EX's Lorentzian letting the S+ pole at 2e-3 leak in, as the reviewer predicted); at a unique ground state EX moves at the last digit only (0.9286496805310763 to ...719).

**Where**: `src/dmrgpy/dcex.py::dynamical_correlator` (its reference is
`wsex[0]` of H rediagonalized in an excited-state basis cached in
`_dcex_excited_cache`, keyed on `(nex, scale, gram_schmidt)`, cleared only by
`__init__` and `restart()`/`set_hamiltonian`, `manybodychain.py:753`), and the
`n_gs` loop in `spinchain.py`, which forwards `submode=` without restriction while
its docstring says "every term ... is averaged".

EX measuring from its own subspace's lowest state is implicit design (the user
guide's EX entry); the hole is that `n_gs` accepts EX. Clearing the cache does not
help: the rediagonalization re-picks inside the degenerate block whatever state
the search grew from.

Repro, from the hunter:

```bash
cd <scratch>/kondo && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 07_ngs_submode_EX.py
```

`kondo/07_ngs_submode_EX.py`:

```python
# get_kondo_spectrum(mode="DMRG", n_gs=2) documents that every term is
# averaged over the manifold, and forwards submode= to the correlator. Does
# the average happen under submode="EX" (dcex), which builds its own
# excited-state basis and caches it per chain (_dcex_excited_cache, keyed on
# (nex, scale, gram_schmidt), not on the ground state)?
# Finding 12's chain at the EXACT crossing (eps=0), itensor_version="python".
# Anchor: mode="ED" (KondoSpectrum's T->0+ equal-weight average), 1.5 at the
# plateau eV=+-0.5 meV; each member alone gives 1-<Sz_0> of that member.
import warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain
print("dmrgpy from", dmrgpy.__file__)
D, TP = 1e-3, 2*np.pi
def chain():
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version="python")
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + D*sc.Sz[0] + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc
eVs = np.array([-0.5e-3, 0.5e-3])
es = np.linspace(-30e-3, 30e-3, 3001)
sc = chain()
_, ed = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="ED")
print("mode=ED /2pi =", np.round(ed/TP, 4))
for submode, extra in (("KPM", {}), ("EX", dict(nex=10))):
    for seed in (1, 2, 3):
        np.random.seed(seed)
        sc = chain()
        kw = dict(site=0, T=0.0, order=2, mode="DMRG", submode=submode,
                  delta=1e-4, es=es, **extra)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, d1 = sc.get_kondo_spectrum(eVs, **kw)
            sz = np.real(sc.vev(sc.Sz[0]))
            _, d2 = sc.get_kondo_spectrum(eVs, n_gs=2, **kw)
        print("submode=%-3s seed %d: n_gs=1 /2pi=%s (1-<Sz_0> of the solved state %.4f)  n_gs=2 /2pi=%s  |n_gs=2 - ED| = %.3f"
              % (submode, seed, np.round(d1/TP, 4), 1-sz, np.round(d2/TP, 4), np.max(np.abs(d2-ed))/TP))
```

`kondo/07_ngs_submode_EX.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
mode=ED /2pi = [1.5 1.5]
submode=KPM seed 1: n_gs=1 /2pi=[1.8355 1.8355] (1-<Sz_0> of the solved state 1.8355)  n_gs=2 /2pi=[1.5 1.5]  |n_gs=2 - ED| = 0.000
submode=KPM seed 2: n_gs=1 /2pi=[1.2073 1.2073] (1-<Sz_0> of the solved state 1.2073)  n_gs=2 /2pi=[1.5 1.5]  |n_gs=2 - ED| = 0.000
submode=KPM seed 3: n_gs=1 /2pi=[1.3261 1.3261] (1-<Sz_0> of the solved state 1.3261)  n_gs=2 /2pi=[1.5 1.5]  |n_gs=2 - ED| = 0.000
submode=EX  seed 1: n_gs=1 /2pi=[1.6772 1.6772] (1-<Sz_0> of the solved state 1.8355)  n_gs=2 /2pi=[1.6772 1.6772]  |n_gs=2 - ED| = 0.177
submode=EX  seed 2: n_gs=1 /2pi=[1.2304 1.2304] (1-<Sz_0> of the solved state 1.2073)  n_gs=2 /2pi=[1.2304 1.2304]  |n_gs=2 - ED| = 0.270
submode=EX  seed 3: n_gs=1 /2pi=[1.319 1.319] (1-<Sz_0> of the solved state 1.3261)  n_gs=2 /2pi=[1.319 1.319]  |n_gs=2 - ED| = 0.181
```

**Reviewer (CONFIRMED, NARROWED)**: counted cache misses and ran `n_gs=2` first on
fresh chains: EX never averages, KPM does (1.499992), and the value lands on both
sides of 1.5:

`reviews/kondo_K3/01_repro_inert.py`:

```python
# Reproduce K3: under submode="EX", get_kondo_spectrum(mode="DMRG", n_gs=2)
# returns the n_gs=1 value bit for bit on finding 12's exact-crossing chain.
# Also: (a) n_gs=2 called FIRST on a fresh chain (no n_gs=1 warm-up filling
# dcex's cache beforehand), counting how many times dcex actually solves for
# excited states (cache misses) across the whole n_gs=2 call; (b) KPM on the
# same chains as the anchor that n_gs=2 does average there.
import warnings, time
import numpy as np
import dmrgpy
from dmrgpy import spinchain, dcex
print("dmrgpy from", dmrgpy.__file__)
D, TP = 1e-3, 2*np.pi
def chain():
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version="python")
    sc.maxm, sc.nsweeps = 30, 15
    h = D*sc.Sz[0]*sc.Sz[0] + D*sc.Sz[0] + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*sc.SS(1, 2)
    sc.set_hamiltonian(h)
    return sc
eVs = np.array([-0.5e-3, 0.5e-3])
es = np.linspace(-30e-3, 30e-3, 3001)
sc = chain()
_, ed = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="ED")
print("mode=ED /2pi =", np.round(ed/TP, 6))
# count dcex cache misses
orig_ges = None
misses = []
_orig = dcex.get_cached_excited_states
def counting(self, n=20, scale=10.0, **k):
    key = (n, scale, getattr(self, "excited_gram_schmidt", False))
    c = getattr(self, "_dcex_excited_cache", None)
    if c is None or c[0] != key: misses.append(1)
    return _orig(self, n=n, scale=scale, **k)
dcex.get_cached_excited_states = counting
kwEX = dict(site=0, T=0.0, order=2, mode="DMRG", submode="EX", delta=1e-4, es=es, nex=10)
kwKPM = dict(site=0, T=0.0, order=2, mode="DMRG", submode="KPM", delta=1e-4, es=es)
for seed in (11, 12, 13):
    np.random.seed(seed)
    sc = chain()
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        misses.clear()
        _, d1 = sc.get_kondo_spectrum(eVs, **kwEX)
        m1 = len(misses); misses.clear()
        _, d2 = sc.get_kondo_spectrum(eVs, n_gs=2, **kwEX)
        m2 = len(misses)
        _, k1 = sc.get_kondo_spectrum(eVs, **kwKPM)
        _, k2 = sc.get_kondo_spectrum(eVs, n_gs=2, **kwKPM)
    print("seed %d  EX: n_gs=1 %s (misses %d)  n_gs=2 %s (misses %d)  bitwise equal %s | KPM: n_gs=1 %s  n_gs=2 %s"
          % (seed, np.round(d1/TP, 6), m1, np.round(d2/TP, 6), m2, np.array_equal(d1, d2),
             np.round(k1/TP, 6), np.round(k2/TP, 6)))
# (a) n_gs=2 first, on a fresh chain
for seed in (21, 22, 23):
    np.random.seed(seed)
    sc = chain()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        misses.clear()
        _, d2 = sc.get_kondo_spectrum(eVs, n_gs=2, **kwEX)
        m2 = len(misses)
        misses.clear()
        _, d1 = sc.get_kondo_spectrum(eVs, **kwEX)
        m1 = len(misses)
    print("seed %d fresh chain, n_gs=2 FIRST: EX n_gs=2 %s (misses %d over 2 members x 3 channels)  then n_gs=1 %s (misses %d)  bitwise equal %s"
          % (seed, np.round(d2/TP, 6), m2, np.round(d1/TP, 6), m1, np.array_equal(d1, d2)))
dcex.get_cached_excited_states = _orig
```

`reviews/kondo_K3/01_repro_inert.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
mode=ED /2pi = [1.5 1.5]
seed 11  EX: n_gs=1 [1.36585 1.36585] (misses 1)  n_gs=2 [1.36585 1.36585] (misses 0)  bitwise equal True | KPM: n_gs=1 [1.065865 1.065865]  n_gs=2 [1.499992 1.499992]
seed 12  EX: n_gs=1 [1.1521 1.1521] (misses 1)  n_gs=2 [1.1521 1.1521] (misses 0)  bitwise equal True | KPM: n_gs=1 [1.873427 1.873427]  n_gs=2 [1.499992 1.499992]
seed 13  EX: n_gs=1 [1.70781 1.70781] (misses 1)  n_gs=2 [1.70781 1.70781] (misses 0)  bitwise equal True | KPM: n_gs=1 [1.759724 1.759724]  n_gs=2 [1.499992 1.499992]
seed 21 fresh chain, n_gs=2 FIRST: EX n_gs=2 [1.33414 1.33414] (misses 1 over 2 members x 3 channels)  then n_gs=1 [1.33414 1.33414] (misses 0)  bitwise equal True
seed 22 fresh chain, n_gs=2 FIRST: EX n_gs=2 [1.317326 1.317326] (misses 1 over 2 members x 3 channels)  then n_gs=1 [1.317326 1.317326] (misses 0)  bitwise equal True
seed 23 fresh chain, n_gs=2 FIRST: EX n_gs=2 [1.842623 1.842623] (misses 1 over 2 members x 3 channels)  then n_gs=1 [1.842623 1.842623] (misses 0)  bitwise equal True
```

EX's error for the state it picked is 0.2 to 0.4 per cent, so the distance from
1.5 is the pick, not EX's accuracy (`02_ex_reference_state.py`). Struck: the
hunter's "0.136 to 0.270 below 1.5" (seed luck; `np.random.seed` does not pin the
pick on `"python"`, so only ranges are quoted). Folded in: the hunter's related
lead (EX at `n_gs=1` 1.6772 against 1 - <Sz_0> = 1.8355 of the chain's state) is
the same defect seen from `n_gs=1`. SECTOR is inert for the same reason by reading
(`sectordc._sector_states` measures from a clone's own solve) and unreachable on
this chain.

**Suggested fix**: now, an allow-list of the submodes known to read the chain's
state, raising `NotImplementedError` for `n_gs>1` under any other (EX, SECTOR),
per section 4.10's rule against a deny-list. Better, in `dcex.py`: keep the one
cached basis (it depends on H only) and replace `wf0 = wsex[0]` by the chain's own
`wf0` reprojected onto the rediagonalized basis, d_n = <n|wf0>, checking
||d||^2 = 1; the prototype gives 1.512780 on every seed against KPM's 1.499992 (the
difference is a leaking S+ pole at 2e-3, predicted to 0.017), and matches the
shipped EX to 1.4e-15 at a unique ground state (`02_ex_reference_state.py` part C,
`04_lorentzian_and_prototype_control.py`). The hunter's "key the cache on the
ground state" is the wrong shape. NUMBERS CHANGE under the reprojection for EX at
`n_gs=1` at any degenerate ground state, since it then reads the chain's own state.

### 16. The public `Many_Body_Chain.get_dynamical_correlator` takes `i=`/`j=` as named parameters and hands them to `str2MO`, which ignores them next to an operator pair, so every wrapper that builds its own pair and forwards `**kwargs` drops a caller's `i=`/`j=` silently; `get_kondo_spectrum(mode="DMRG", i=1, j=1)` returns the site-0 spectrum bit for bit where the same call on `mode="ED"` raises

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `kondo` &middot; origin `1b87543` (the signature), `30200a4` (the ED/DMRG asymmetry)

**Status**: FIXED. `Many_Body_Chain.get_dynamical_correlator` defaults to `i=None, j=None`, raises `TypeError` when either is given with a non-string `name`, and passes 0 in their place on the string branch; `kpmdmrg.py` and the lower-level `_MB` route are untouched. Every inheriting wrapper, the Kondo route included, now raises on the same call. Pinned by `tests/test_audit_2026_09_24b_kpm.py::test_sites_next_to_an_operator_pair_raise` (`i=1, j=1` and `i=0, j=0` on DMRG and ED) and `::test_string_name_still_honours_the_sites`. No number changes: a pair with `i=1, j=1` returned C[Sz0,Sz0] exactly (0.306 from C[Sz1,Sz1]) on both modes and now raises, while `"ZZ", i=1, j=1` still gives C[Sz1,Sz1] with a difference of 0.

**Where**: `src/dmrgpy/manybodychain.py::get_dynamical_correlator(self,
mode="DMRG", name=None, i=0, j=0, **kwargs)`. Inheriting wrappers in `src/` (an
AST scan): the two Kondo terms (`kondospectrumtk/secondorder_dc.py:46`,
`potentialdc.py:59`), three calls in `atomtk/iets.py`,
`dynamicstk/spincorrelators.py:25`, `fermionchaintk/dynamicalcorrelator.py:15`;
by reading also `dmrgpy.atom.get_spinflip`.

`str2MO` returning a pair unchanged is documented design (`i`/`j` mean something
only next to a string name); the fault is the integer defaults, which cannot tell
"not passed" from "passed 0". It is neither the previous record's finding 10 nor
its open ROOTN item (a string name on the lower-level route, where `i`/`j` should
be honoured) nor ED TD's swallowed keywords (a `**kwargs` with no consumer). What
`30200a4` created is the asymmetry: since finding 11's fix `mode="ED"` raises on
exactly this call, and finding 11's text ("On `mode="DMRG"` the same typos raise,
only because the dynamical correlator downstream is strict") is false for `i` and
`j`.

Repro, from the hunter:

```bash
cd <scratch>/kondo && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 09_dmrg_route_swallows_i_j.py
```

`kondo/09_dmrg_route_swallows_i_j.py`:

```python
# Finding 11 made mode="ED" raise on any keyword it does not read. Does
# mode="DMRG" still swallow some? Candidates: i= and j=, which
# Many_Body_Chain.get_dynamical_correlator consumes for a string name= and
# drops when name= is already an operator pair (which is what
# second_order_dIdV_dc/third_order_potential_dIdV_dc pass).
# 3-site chain of the kondo examples (field on site 0 only, so sites 0 and
# 1 have different spectra), itensor_version="python". Anchor: mode="ED" at
# site=0 and site=1.
import numpy as np
import dmrgpy
from dmrgpy import spinchain
print("dmrgpy from", dmrgpy.__file__)
G, MUB, TP = 2.0, 5.7883818066e-5, 2*np.pi
np.random.seed(0)
sc = spinchain.Spin_Chain(["1/2"]*3, itensor_version="python")
sc.maxm, sc.nsweeps = 30, 15
h = G*MUB*10.0*sc.Sz[0]
for k in range(2): h = h + 0.01*sc.SS(k, k+1)
sc.set_hamiltonian(h)
eVs = np.array([-2e-3, 0.0, 2e-3])
es = np.linspace(-20e-3, 20e-3, 5328)
kw = dict(T=0.0, order=2, mode="DMRG", submode="KPM", delta=2e-5, es=es)
_, ed0 = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="ED")
_, ed1 = sc.get_kondo_spectrum(eVs, site=1, T=0.0, order=2, mode="ED")
print("mode=ED   site=0 /2pi:", np.round(ed0/TP, 4), "  site=1 /2pi:", np.round(ed1/TP, 4))
try:
    sc.get_kondo_spectrum(eVs, i=1, T=0.0, order=2, mode="ED")
    print("mode=ED   i=1: accepted")
except TypeError as e:
    print("mode=ED   i=1: TypeError:", str(e)[:90])
_, d0 = sc.get_kondo_spectrum(eVs, site=0, **kw)
_, d1 = sc.get_kondo_spectrum(eVs, site=1, **kw)
_, di = sc.get_kondo_spectrum(eVs, i=1, j=1, **kw)
print("mode=DMRG site=0 /2pi:", np.round(d0/TP, 4))
print("mode=DMRG site=1 /2pi:", np.round(d1/TP, 4))
print("mode=DMRG i=1,j=1 /2pi:", np.round(di/TP, 4),
      " max|i=1 - site=0| = %.3e   max|i=1 - site=1| = %.3e" % (np.max(np.abs(di-d0)), np.max(np.abs(di-d1))))
for bad in (dict(Jrho=0.05), dict(ngs=2), dict(sub_mode="CVM")):
    try:
        sc.get_kondo_spectrum(eVs, site=0, **kw, **bad)
        print("mode=DMRG %s: accepted" % bad)
    except TypeError as e:
        print("mode=DMRG %s: TypeError: %s" % (bad, str(e)[:80]))
```

`kondo/09_dmrg_route_swallows_i_j.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
mode=ED   site=0 /2pi: [0.3487 0.1286 0.3487]   site=1 /2pi: [0.0892 0.0331 0.0892]
mode=ED   i=1: TypeError: get_kondo_spectrum(mode="ED") got unexpected keyword argument(s): i. They are either mode=
mode=DMRG site=0 /2pi: [0.3487 0.1286 0.3487]
mode=DMRG site=1 /2pi: [0.0892 0.0331 0.0892]
mode=DMRG i=1,j=1 /2pi: [0.3487 0.1286 0.3487]  max|i=1 - site=0| = 0.000e+00   max|i=1 - site=1| = 1.631e+00
mode=DMRG {'Jrho': 0.05}: TypeError: Unexpected keyword argument(s) for the KPM dynamical correlator: Jrho. Solver pa
mode=DMRG {'ngs': 2}: TypeError: Unexpected keyword argument(s) for the KPM dynamical correlator: ngs. Solver par
mode=DMRG {'sub_mode': 'CVM'}: TypeError: Unexpected keyword argument(s) for the KPM dynamical correlator: sub_mode. Solve
```

**Reviewer (CONFIRMED, NARROWED)**: moved the defect one frame up, to the public
method: independent of solver and submode (DMRG KPM, DMRG CVM, ED KPM and ED ED
each return C[Sz0,Sz0] bit for bit for a pair plus `i=1, j=1`), and on the Kondo
route both correlator-based terms drop it:

`reviews/kondo_K4/02_where_ij_die.py`:

```python
# Reviewer probe for K4: where do i=/j= die, and is it the Kondo route or the
# public correlator? Same 3-site chain as the hunter (field on site 0 only, so
# sites 0 and 1 differ), itensor_version="python", sweeps pinned, seeded.
#  A. public Many_Body_Chain.get_dynamical_correlator with an explicit pair
#     (Sz0,Sz0) plus i=1,j=1: on mode="DMRG" submode KPM and CVM, and on
#     mode="ED". Expect bit-identical to (Sz0,Sz0), i.e. silently ignored.
#  B. get_kondo_spectrum(mode="DMRG", name=...): expect Python's own
#     duplicate-keyword TypeError (the hunter's fix adds name= to its list).
#  C. get_kondo_spectrum(mode="DMRG", i=1) alone, without j.
#  D. third_order_potential_dIdV_dc called directly with i=1: the other
#     correlator-based term, same mechanism.
import warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
print("dmrgpy from", dmrgpy.__file__)
G, MUB, TP = 2.0, 5.7883818066e-5, 2*np.pi
np.random.seed(0)
sc = spinchain.Spin_Chain(["1/2"]*3, itensor_version="python")
sc.maxm, sc.nsweeps = 30, 15
h = G*MUB*10.0*sc.Sz[0]
for k in range(2): h = h + 0.01*sc.SS(k, k+1)
sc.set_hamiltonian(h)

def cmp(label, y, y0, y1):
    print("%-44s max|y-C00| = %.3e   max|y-C11| = %.3e" %
          (label, np.max(np.abs(y-y0)), np.max(np.abs(y-y1))))

# A. public correlator, pair plus i/j
esA = np.linspace(-5e-3, 25e-3, 60)
P0, P1 = (sc.Sz[0], sc.Sz[0]), (sc.Sz[1], sc.Sz[1])
for mode, sub, kw in (("DMRG", "KPM", dict(delta=1e-3)),
                      ("DMRG", "CVM", dict(delta=1e-3)),
                      ("ED", "KPM", dict(delta=1e-3)),
                      ("ED", "ED", dict(delta=1e-3))):
    _, y0 = sc.get_dynamical_correlator(mode=mode, submode=sub, name=P0, es=esA, **kw)
    _, y1 = sc.get_dynamical_correlator(mode=mode, submode=sub, name=P1, es=esA, **kw)
    _, yi = sc.get_dynamical_correlator(mode=mode, submode=sub, name=P0, i=1, j=1, es=esA, **kw)
    _, ys = sc.get_dynamical_correlator(mode=mode, submode=sub, name="ZZ", i=1, j=1, es=esA, **kw)
    print("[A] mode=%s submode=%s  max|C00-C11| = %.3e" % (mode, sub, np.max(np.abs(y0-y1))))
    cmp("    name=(Sz0,Sz0), i=1, j=1", yi, y0, y1)
    cmp("    name='ZZ', i=1, j=1 (string form)", ys, y0, y1)

# B, C. the Kondo DMRG route
eVs = np.array([-2e-3, 0.0, 2e-3])
es = np.linspace(-20e-3, 20e-3, 5328)
kw = dict(T=0.0, order=2, mode="DMRG", submode="KPM", delta=2e-5, es=es)
try:
    sc.get_kondo_spectrum(eVs, site=0, name=P1, **kw)
    print("[B] kondo DMRG name=(Sz1,Sz1): accepted")
except TypeError as e:
    print("[B] kondo DMRG name=(Sz1,Sz1): TypeError:", str(e)[:100])
_, d0 = sc.get_kondo_spectrum(eVs, site=0, **kw)
_, d1 = sc.get_kondo_spectrum(eVs, site=1, **kw)
_, di = sc.get_kondo_spectrum(eVs, i=1, **kw)
_, dj = sc.get_kondo_spectrum(eVs, site=1, i=0, j=0, **kw)
print("[C] kondo DMRG site=0 /2pi:", np.round(d0/TP, 4), " site=1 /2pi:", np.round(d1/TP, 4))
cmp("[C] kondo DMRG i=1 (no j)", di, d0, d1)
print("[C] kondo DMRG site=1,i=0,j=0: max|y-site1| = %.3e" % np.max(np.abs(dj-d1)))

# D. potential term directly
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    pkw = dict(T0=1.0, mode="DMRG", submode="KPM", delta=2e-5, es=es)
    p0 = third_order_potential_dIdV_dc(sc, 0, eVs, -0.05, 0.25, **pkw)
    p1 = third_order_potential_dIdV_dc(sc, 1, eVs, -0.05, 0.25, **pkw)
    pi = third_order_potential_dIdV_dc(sc, 0, eVs, -0.05, 0.25, i=1, j=1, **pkw)
print("[D] potential site=0:", np.round(p0, 5), " site=1:", np.round(p1, 5))
cmp("[D] potential site=0, i=1, j=1", pi, p0, p1)
```

`reviews/kondo_K4/02_where_ij_die.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
[A] mode=DMRG submode=KPM  max|C00-C11| = 9.705e+01
    name=(Sz0,Sz0), i=1, j=1                 max|y-C00| = 0.000e+00   max|y-C11| = 9.705e+01
    name='ZZ', i=1, j=1 (string form)        max|y-C00| = 9.705e+01   max|y-C11| = 0.000e+00
[A] mode=DMRG submode=CVM  max|C00-C11| = 4.735e+01
    name=(Sz0,Sz0), i=1, j=1                 max|y-C00| = 0.000e+00   max|y-C11| = 4.735e+01
    name='ZZ', i=1, j=1 (string form)        max|y-C00| = 4.735e+01   max|y-C11| = 0.000e+00
[A] mode=ED submode=KPM  max|C00-C11| = 9.707e+01
    name=(Sz0,Sz0), i=1, j=1                 max|y-C00| = 0.000e+00   max|y-C11| = 9.707e+01
    name='ZZ', i=1, j=1 (string form)        max|y-C00| = 9.707e+01   max|y-C11| = 0.000e+00
[A] mode=ED submode=ED  max|C00-C11| = 4.735e+01
    name=(Sz0,Sz0), i=1, j=1                 max|y-C00| = 0.000e+00   max|y-C11| = 4.735e+01
    name='ZZ', i=1, j=1 (string form)        max|y-C00| = 4.735e+01   max|y-C11| = 0.000e+00
[B] kondo DMRG name=(Sz1,Sz1): TypeError: dmrgpy.manybodychain.Many_Body_Chain.get_dynamical_correlator() got multiple values for keyword argu
[C] kondo DMRG site=0 /2pi: [0.3487 0.1286 0.3487]  site=1 /2pi: [0.0892 0.0331 0.0892]
[C] kondo DMRG i=1 (no j)                    max|y-C00| = 0.000e+00   max|y-C11| = 1.631e+00
[C] kondo DMRG site=1,i=0,j=0: max|y-site1| = 0.000e+00
[D] potential site=0: [ 0.03912 -0.      -0.03912]  site=1: [ 0.02248 -0.      -0.02248]
[D] potential site=0, i=1, j=1               max|y-C00| = 0.000e+00   max|y-C11| = 1.664e-02
```

Kept with a caveat: "3.9 times the site-1 value" is this chain's own site contrast
(the field is on site 0), not the defect's size, which is "the `site=` spectrum
regardless of `i`". Struck from the fix: "(and `name`)", since Python already raises
on a duplicated `name`. `site=` is the only spelling anywhere (46 call sites), but
the docstring sends a reader to exactly the signature where `i`/`j` live.

**Suggested fix**: in `Many_Body_Chain.get_dynamical_correlator`, default `i=None,
j=None`, raise `TypeError` when either is given with a non-string `name`, and pass
0 in their place on the string branch. An AST scan of 292 calls across `src/`,
`tests/`, `examples/` and `benchmarks/` finds none that passes a non-literal `name=`
with `i=`/`j=` (`reviews/kondo_K4/00_ast_scan.py`). One behaviour change:
`get_kondo_spectrum(site=1, i=0, j=0)`, harmless today, raises. Correct the
previous record's finding 11 sentence. No number changes.

### 17. Both correlator-based Kondo terms assume an increasing `es`, which nothing says and which `30200a4`'s docstrings now invite breaking: on a refined block appended after a coarse grid, the second-order term is silently 3.03 off a 6.97 peak and the potential term 0.654 off a 0.671 peak, where the same points sorted are right

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `kondo` &middot; origin pre-`30200a4` in the code, `30200a4` in the docstrings

**Status**: FIXED. `_cumulative_theta0_weight` and `_convolved_F0_weight` `np.asarray` the grid and sort `x` and `S` together with a stable `np.argsort`; the docstrings now say "may be non-uniform and in any order", and the second-order docstring's "the cumulative integral starts at es[0]" reads "the lowest point of es". Pinned by `tests/test_audit_2026_09_24b_misc.py::test_kondo_terms_do_not_depend_on_the_order_of_es` (coarse then fine, fine then coarse, descending: agreement with the sorted grid to 1e-12 of the peak on both terms, and no sum-rule warning) and `::test_kondo_terms_accept_es_as_a_list`. NUMBERS CHANGE only on non-increasing grids, which were wrong: on the single S=1/2 at 10 T (`mode="ED"`, `submode="ED"`), coarse block then fine block, the potential term goes from 0.6536 to 0.0005091 off exact and the second-order term from 3.035 to 0.03183; fine then coarse was 0.6633 and 3.11, descending 1.341 and 10.99, and all three now give the sorted answer; increasing grids are `array_equal` before and after.

**Where**: `src/dmrgpy/kondospectrumtk/secondorder_dc.py::_cumulative_theta0_weight`
(`cumulative_trapezoid(S, x)` then `np.interp(eVs, x, cum)`, undefined over a
non-increasing `xp`) and `potentialdc.py::_convolved_F0_weight` (trapezoid weights
from `np.diff(x)`).

`get_dynamical_correlator` returns `x` exactly as the caller passed `es` and
evaluates `S` pointwise on every route the helpers reach, so a non-increasing `es`
arrives unchanged. No docstring says it must be increasing (the only hint is
"the cumulative integral starts at es[0]"), while `secondorder_dc.py:79` ("may be
non-uniform"), `potentialdc.py:53` and the `es` docstring ("may vary along the
grid") and `docs/user_guide.md:3266` invite a non-uniform grid, which is naturally
built by concatenation, and `get_kondo_spectrum(mode="DMRG")` computes the
second-order term on every call whatever `order` is.

Repro, from the hunter, on the exact Lehmann route (`mode="ED"`, `submode="ED"`,
the single S=1/2 at 10 T of `tests/test_audit_2026_09_24_kondo.py`):

```bash
cd <scratch>/kondo && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 10_es_ordering.py
```

`kondo/10_es_ordering.py`:

```python
# Finding 13's fix: trapezoid weights on the correlator's own grid. What
# happens on a grid that is not increasing -- descending, a refined block
# appended after the coarse one (out of order), a grid with duplicated
# points, and es passed as a list? Exact Lehmann correlator
# (mode="ED", submode="ED") on the single S=1/2 at 10 T of
# tests/test_audit_2026_09_24_kondo.py, so only the quadrature differs.
# Anchor: conductance.third_order_potential_dIdV / second_order_dIdV
# (explicit excited-state sums, no grid).
import warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV, second_order_dIdV
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
print("dmrgpy from", dmrgpy.__file__)
G, MUB = 2.0, 5.7883818066e-5
sc = spinchain.Spin_Chain(["1/2"])
sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
ks = KondoSpectrum(sc, site=0, T=0.0)
eVs = np.linspace(-1e-3, 2e-3, 21)
Jrho_s, U, delta = 0.1, 0.3, 2e-6
ref = third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0)
ref2 = second_order_dIdV(ks, eVs, T0=1.0, U=U)
uni = np.linspace(-1e-3, 3e-3, 40_000)
coarse = np.linspace(-1e-3, 3e-3, 8000)
fine = np.linspace(1.0e-3, 1.4e-3, 4000)
grids = {
    "increasing (the test's)": uni,
    "descending": uni[::-1].copy(),
    "coarse block then refined block (unsorted)": np.concatenate([coarse, fine]),
    "same points, sorted": np.sort(np.concatenate([coarse, fine])),
    "every point duplicated": np.repeat(uni, 2),
    "as a Python list": list(uni),
}
kw = dict(T0=1.0, mode="ED", submode="ED", delta=delta)
print("peaks: potential %.4f, second order %.4f" % (np.max(np.abs(ref)), np.max(np.abs(ref2))))
for label, es in grids.items():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        p = third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, es=es, **kw)
    nw = len([w for w in caught if "third_order_potential_dIdV_dc" in str(w.message)])
    s = second_order_dIdV_dc(sc, 0, eVs, U=U, es=es, **kw)
    print("%-44s potential max|dc-exact| = %.4g (warn=%d)   second order max|dc-exact| = %.4g"
          % (label, np.max(np.abs(p-ref)), nw, np.max(np.abs(s-ref2))))
```

`kondo/10_es_ordering.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
peaks: potential 0.6709, second order 6.9743
increasing (the test's)                      potential max|dc-exact| = 0.0005091 (warn=0)   second order max|dc-exact| = 0.03183
descending                                   potential max|dc-exact| = 1.341 (warn=1)   second order max|dc-exact| = 10.99
coarse block then refined block (unsorted)   potential max|dc-exact| = 0.6536 (warn=1)   second order max|dc-exact| = 3.035
same points, sorted                          potential max|dc-exact| = 0.0005091 (warn=0)   second order max|dc-exact| = 0.03183
every point duplicated                       potential max|dc-exact| = 0.0005091 (warn=0)   second order max|dc-exact| = 0.03183
as a Python list                             potential max|dc-exact| = 0.0005091 (warn=0)   second order max|dc-exact| = 0.03183
```

**Reviewer (CONFIRMED, NARROWED)**: the correlator upstream is exactly
order-preserving (0.00e+00 on ED and on DMRG KPM for reversed, permuted and
concatenated grids, `02_correlator_order_preserving.py`), so the defect lives only
in the two helpers; the same holds on `mode="DMRG"`, `submode="KPM"`, v3 (3.40
against 0.54 sorted, `04_dmrg_route_and_warning.py`). Measured against the same
points sorted, against the pre-`30200a4` potential term, and with a sort
monkeypatched into both helpers:

`reviews/kondo_K5/03_size_prefix_and_fix.py`:

```python
# Size of the defect against the SAME points sorted, the pre-30200a4
# potential code on the same grids, the sum-rule weight the warning sees,
# and the hunter's fix (sort x and S together inside both helpers)
# applied by monkeypatching. Same chain as the finding.
import sys, warnings
import numpy as np
import dmrgpy
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV, second_order_dIdV
from dmrgpy.kondospectrumtk import potentialdc, secondorder_dc
sys.path.insert(0, ".")
import old_potentialdc_30200a4parent as oldp
print("dmrgpy from", dmrgpy.__file__)
G, MUB = 2.0, 5.7883818066e-5
sc = spinchain.Spin_Chain(["1/2"])
sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
ks = KondoSpectrum(sc, site=0, T=0.0)
eVs = np.linspace(-1e-3, 2e-3, 21)
Jrho_s, U, delta = 0.1, 0.3, 2e-6
ref = third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0)
ref2 = second_order_dIdV(ks, eVs, T0=1.0, U=U)
kw = dict(T0=1.0, mode="ED", submode="ED", delta=delta)
uni = np.linspace(-1e-3, 3e-3, 40_000)
coarse = np.linspace(-1e-3, 3e-3, 8000)
fine = np.linspace(1.0e-3, 1.4e-3, 4000)
cat = np.concatenate([coarse, fine])
rng = np.random.default_rng(3)
grids = {
    "increasing uniform": uni,
    "descending uniform": uni[::-1].copy(),
    "coarse then fine (unsorted)": cat,
    "coarse then fine, sorted": np.sort(cat),
    "fine then coarse (unsorted)": np.concatenate([fine, coarse]),
    "random permutation of uniform": uni[rng.permutation(len(uni))],
}
d = np.diff(cat)
print("concat grid: coarse dx %.3g, fine dx %.3g, dx at the join %.3g (index %d)"
      % (coarse[1]-coarse[0], fine[1]-fine[0], d[len(coarse)-1], len(coarse)-1))
print("line (Zeeman) at %.5g eV; peaks: potential %.4f, second order %.4f"
      % (G*MUB*10.0, np.max(np.abs(ref)), np.max(np.abs(ref2))))

def run(es, label, tag):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        p = potentialdc.third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, es=es, **kw)
    msg = [str(w.message) for w in caught if "third_order_potential_dIdV_dc" in str(w.message)]
    held = msg[0].split("holds ")[1].split(" of")[0] if msg else "-"
    s = secondorder_dc.second_order_dIdV_dc(sc, 0, eVs, U=U, es=es, **kw)
    return np.max(np.abs(p-ref)), held, np.max(np.abs(s-ref2)), p, s

# weight = sum_k int S_kk, computed directly (no warning threshold)
def weight(es):
    tot = 0.
    for op in (sc.Sx[0], sc.Sy[0], sc.Sz[0]):
        tot += potentialdc._convolved_F0_weight(sc, op, np.array([0.]), 20e-3, 5e-6,
                                                "ED", "ED", delta, es)[1]
    return tot

res_now = {}
print("\n-- current HEAD --")
for label, es in grids.items():
    ep, held, es2, p, s = run(es, label, "now")
    res_now[label] = (p, s)
    print("%-32s potential err %.4g (warn says held=%s, sum rule weight %.4f of 0.75)   second order err %.4g"
          % (label, ep, held, weight(es), es2))

print("\n-- pre-30200a4 potential term (first spacing times plain sum) --")
for label, es in grids.items():
    p = oldp.third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, es=es, **kw)
    print("%-32s potential err %.4g" % (label, np.max(np.abs(p-ref))))

# the hunter's fix, sort x and S together after the correlator call
orig_F0, orig_th = potentialdc._convolved_F0_weight, secondorder_dc._cumulative_theta0_weight
def _sorted_dc(chain):
    real = chain.get_dynamical_correlator
    def wrapped(*a, **k):
        x, S = real(*a, **k)
        x = np.asarray(x, dtype=float); S = np.asarray(S)
        i = np.argsort(x, kind="stable")
        return x[i], S[i]
    return wrapped
class Proxy:
    def __init__(self, c): self._c = c; self.get_dynamical_correlator = _sorted_dc(c)
    def __getattr__(self, n): return getattr(self._c, n)
potentialdc._convolved_F0_weight = lambda chain, *a, **k: orig_F0(Proxy(chain), *a, **k)
secondorder_dc._cumulative_theta0_weight = lambda chain, *a, **k: orig_th(Proxy(chain), *a, **k)
print("\n-- with x,S sorted inside both helpers --")
base_p, base_s = res_now["increasing uniform"]
for label, es in grids.items():
    ep, held, es2, p, s = run(es, label, "fix")
    print("%-32s potential err %.4g (warn held=%s)   second order err %.4g   |fix - HEAD| on increasing: %s"
          % (label, ep, held, es2,
             "%.1e/%.1e" % (np.max(np.abs(p-base_p)), np.max(np.abs(s-base_s)))
             if label == "increasing uniform" else "-"))
potentialdc._convolved_F0_weight, secondorder_dc._cumulative_theta0_weight = orig_F0, orig_th

print("\n-- eVs order (HEAD, increasing es): results follow the caller's eVs --")
perm = rng.permutation(len(eVs))
p1 = potentialdc.third_order_potential_dIdV_dc(sc, 0, eVs[perm], Jrho_s, U, es=uni, **kw)
s1 = secondorder_dc.second_order_dIdV_dc(sc, 0, eVs[perm], U=U, es=uni, **kw)
print("permuted eVs: potential max|p(eVs[perm]) - p(eVs)[perm]| = %.2e, second order %.2e"
      % (np.max(np.abs(p1-base_p[perm])), np.max(np.abs(s1-base_s[perm]))))
```

`reviews/kondo_K5/03_size_prefix_and_fix.out`:

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
concat grid: coarse dx 5e-07, fine dx 1e-07, dx at the join -0.002 (index 7999)
line (Zeeman) at 0.0011577 eV; peaks: potential 0.6709, second order 6.9743

-- current HEAD --
increasing uniform               potential err 0.0005091 (warn says held=-, sum rule weight 0.7495 of 0.75)   second order err 0.03183
descending uniform               potential err 1.341 (warn says held=-0.7495, sum rule weight -0.7495 of 0.75)   second order err 10.99
coarse then fine (unsorted)      potential err 0.6536 (warn says held=1.233, sum rule weight 1.2331 of 0.75)   second order err 3.035
coarse then fine, sorted         potential err 0.0005091 (warn says held=-, sum rule weight 0.7495 of 0.75)   second order err 0.03183
fine then coarse (unsorted)      potential err 0.6633 (warn says held=1.239, sum rule weight 1.2393 of 0.75)   second order err 3.11
random permutation of uniform    potential err 206 (warn says held=-355, sum rule weight -355.0291 of 0.75)   second order err 4463

-- pre-30200a4 potential term (first spacing times plain sum) --
increasing uniform               potential err 0.0005091
descending uniform               potential err 1.341
coarse then fine (unsorted)      potential err 3.342
coarse then fine, sorted         potential err 3.342
fine then coarse (unsorted)      potential err 0.1317
random permutation of uniform    potential err 317.1

-- with x,S sorted inside both helpers --
increasing uniform               potential err 0.0005091 (warn held=-)   second order err 0.03183   |fix - HEAD| on increasing: 0.0e+00/0.0e+00
descending uniform               potential err 0.0005091 (warn held=-)   second order err 0.03183   |fix - HEAD| on increasing: -
coarse then fine (unsorted)      potential err 0.0005091 (warn held=-)   second order err 0.03183   |fix - HEAD| on increasing: -
coarse then fine, sorted         potential err 0.0005091 (warn held=-)   second order err 0.03183   |fix - HEAD| on increasing: -
fine then coarse (unsorted)      potential err 0.0005091 (warn held=-)   second order err 0.03183   |fix - HEAD| on increasing: -
random permutation of uniform    potential err 0.0005091 (warn held=-)   second order err 0.03183   |fix - HEAD| on increasing: -

-- eVs order (HEAD, increasing es): results follow the caller's eVs --
permuted eVs: potential max|p(eVs[perm]) - p(eVs)[perm]| = 0.00e+00, second order 0.00e+00
```

Struck or qualified: "pre-existing in part" becomes "pre-existing in full" (the
pre-`30200a4` potential term gives the same 1.341 on a descending grid and was
wrong on the sorted concatenation too, the previous record's finding 13), and what
is new is only the docstrings' invitation; the descending grid, an exact sign flip,
is contrived and should not carry the headline; the potential term's "widen es"
advice is wrong for an excess of weight, but not specific to ordering (a sorted,
under-resolved uniform grid trips it the same way), so it is a separate lead.

**Suggested fix**: sort `x` and `S` together with a stable `np.argsort` at the top
of both helpers (`_cumulative_theta0_weight` needs an `np.asarray` first); raising
instead would forbid the concatenation the docstrings invite. It restores the
increasing-grid answer exactly on every ordering and moves increasing grids by
0.0e+00. Change "may be non-uniform" to "may be non-uniform and in any order" in
both docstrings and in the user guide in both formats. NUMBERS CHANGE only on
non-increasing grids, which were wrong.

### 18. The padding strip that `30200a4` added at the entry of the one-site TDVP route reads the process-global `set_pad_bonds` flag rather than the state it is handed, so a state padded earlier and evolved with the flag off runs the previous record's finding 16 algorithm: 0.4929 off the unpadded trajectory from a Neel start, 1.5e-6 from an entangled one

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `pyitensor` &middot; origin `30200a4`

**Status**: FIXED, with the reviewer's `fixS`. `_strip_bond_padding` keys on the state: with the flag on the lossless suspended `position(n)`/`position(1)` sweep runs in place as before; with it off the same sweep runs on a copy (`_Chain.copy()` copies the tensor list and the centre, and nothing mutates an ITensor in place) and is adopted only if some bond shrank, so an unpadded state comes back untouched (`_strip_sweep`, `_bond_dims`). `set_pad_bonds`' docstring in `pyitensor/backend.py` now says the strip keys on the state. Pinned by `tests/test_audit_2026_09_24b_misc.py::test_padded_state_evolved_with_the_flag_off_follows_the_unpadded_run` (flag cleared and `pad_bonds_suspended()`), `::test_padded_entangled_state_evolved_with_the_flag_off_follows_the_unpadded_run` and `::test_unpadded_runs_are_bit_identical` (Neel and dimer starts at `tdvp_gse_sweeps` 0 and 3). The cost is one extra lossless SVD sweep on a copy at the entry of every unpadded `TDVP_GSE` trajectory. NUMBERS CHANGE only in the padded-then-unpadded configuration: on the Neel quench (XXZ Delta=0.7, hz=0.1, dt=0.05) at n=10, K=8, 40 steps, `tdvp_gse_sweeps=0`, the distance to the unpadded run goes from 4.929e-01 to 1.166e-15 (n=8, K=4: from 1.925e-01 to 7.772e-15; dimer start: from 7.034e-08 to 6.578e-15), and at `tdvp_gse_sweeps=3` it stays at the gauge level (2.540e-09 to 3.178e-09, the reviewer's own `fixS` value).

**Where**: `src/dmrgpy/pyitensor/chain.py::_strip_bond_padding` (returns early on
`if not _bk.pad_bonds()`), called from `Chain.evolve_and_measure_tdvp_gse`. Reached
through `timedependent.evolve_and_measure(..., tevol_method="TDVP_GSE")` with an
explicit padded `wf=`, or with the chain's cached padded `wf0` via `h=H1` and no
`set_hamiltonian` in between; the flag off either by `set_pad_bonds(None)` or by
`backend.pad_bonds_suspended()`, which give identical numbers.

With the padding intact at `tdvp_gse_sweeps=0`, `qr_split` completes the padded
zero directions into live ones, finding 16's mechanism, on a bond-K manifold. The
affected run is closer to ED (4.9e-8, where the unpadded one-site run is the frozen
product state, the correct one-site answer, 0.4929 from ED), so the defect is that
it is not the method asked for and that the same state object gives two
trajectories depending on hidden global state, not an accuracy loss. It contradicts
`set_pad_bonds`' docstring ("strip the padded zeros from the evolved state once at
trajectory entry") and the user guide ("the padding being stripped once at
trajectory entry"). It is visible only at `tdvp_gse_sweeps=0` with the true bond
dimension below K; at the default 3 the first expansion's suspended SVDs strip the
padding anyway (gauge level, 2.5e-9 and 4.7e-11). `quench_tdvp_gse` is immune
(`_apply_mpo` rebuilds A|gs> at its true rank). It survived because every
finding-16 test sets the flag before the ground state and clears it only in a
`finally`, after the evolution.

**Expected**: the unpadded trajectory.

Repro, from the hunter, XXZ (Delta=0.7) plus hz=0.1 quench on 10 sites, K=maxm=8,
40 steps of dt=0.05:

```bash
cd <scratch>/pyitensor && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 05_padded_state_unpadded_run.py
```

`pyitensor/05_padded_state_unpadded_run.py`:

```python
# _strip_bond_padding is gated on the CURRENT global flag (backend.pad_bonds()),
# not on whether the state it is handed carries padded zeros. Scenario: the
# ground state is computed under set_pad_bonds(K) (as the backend docstring
# describes: "the ground state before it ... keep[s its] padding"), padding is
# then switched off, and TDVP_GSE evolves that cached, padded wf0.
# Anchors: the fully unpadded run (the method as asked for), and numpy ED.
# Same model and start as 04: n=10 XXZ (Delta=0.7) + hz=0.1 quench from Neel,
# maxm=K=8, 40 steps of dt=0.05, <Sz_0>(t).
import sys
import numpy as np
from scipy.linalg import expm
from dmrgpy import spinchain, timedependent
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor import chain as pychain

N = 10; MAXM = 8; NT = 40; DT = 0.05; DELTA = 0.7; HZ = 0.1


def run(pad_gs, pad_evolve, sweeps):
    np.random.seed(11)
    try:
        bk.set_pad_bonds(pad_gs)
        sc = spinchain.Spin_Chain([2]*N, itensor_version="python")
        sc.maxm = MAXM; sc.nsweeps = 10; sc.cutoff = 1e-12
        h0 = 0
        for i in range(N):
            h0 = h0 + (-1)**i*sc.Sz[i]
        sc.set_hamiltonian(h0)
        wf = sc.get_gs()
        links = [wf.cpp_handle.A(i).inds for i in range(1, N+1)]
        dims = [max(ix.dim for ix in inds if ix.hastags("Link")) for inds in links]
        h = 0
        for i in range(N-1):
            h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + DELTA*sc.Sz[i]*sc.Sz[i+1]
        for i in range(N):
            h = h + HZ*sc.Sz[i]
        sc.set_hamiltonian(h)
        sc.tevol_method = "TDVP_GSE"; sc.tdvp_gse_sweeps = sweeps
        bk.set_pad_bonds(pad_evolve)
        out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=NT, dt=DT,
                                               wf=wf, mode="DMRG", return_wf=True)
        wf_end = out[2]
        ent = wf_end.get_bond_entropy(N//2 - 1) if hasattr(wf_end, "get_bond_entropy") else float("nan")
        return np.real(np.asarray(out[1])), dims, ent
    finally:
        bk.set_pad_bonds(None)


sx = np.array([[0, 1], [1, 0]])/2; sy = np.array([[0, -1j], [1j, 0]])/2; sz = np.diag([0.5, -0.5])
def op(o, i):
    m = np.array([[1.0]])
    for j in range(N):
        m = np.kron(m, o if j == i else np.eye(2))
    return m
H = sum(op(sx, i) @ op(sx, i+1) + op(sy, i) @ op(sy, i+1) + DELTA*op(sz, i) @ op(sz, i+1) for i in range(N-1))
H = H + HZ*sum(op(sz, i) for i in range(N))
idx = int("".join("1" if i % 2 == 0 else "0" for i in range(N)), 2)
psi = np.zeros(2**N, dtype=complex); psi[idx] = 1.0
U = expm(-1j*DT*H); Sz0 = op(sz, 0); ED = []
for _ in range(NT):
    ED.append(np.vdot(psi, Sz0 @ psi).real); psi = U @ psi
ED = np.asarray(ED)

for sweeps in (0, 3):
    ref, d0, e0 = run(None, None, sweeps)
    print("sweeps=%d  gs pad=None evolve pad=None : wf0 link dims %s  <Sz_0>(end)=%.6f max|DMRG-ED|=%.3e"
          % (sweeps, d0, ref[-1], np.max(np.abs(ref-ED))))
    for pg, pe in ((MAXM, MAXM), (MAXM, None)):
        y, d, e = run(pg, pe, sweeps)
        print("sweeps=%d  gs pad=%-4s evolve pad=%-4s: wf0 link dims %s  <Sz_0>(end)=%.6f "
              "max|run-unpadded|=%.3e max|DMRG-ED|=%.3e"
              % (sweeps, pg, pe, d, y[-1], np.max(np.abs(y-ref)), np.max(np.abs(y-ED))))
        sys.stdout.flush()
```

`pyitensor/05_padded_state_unpadded_run.out`:

```
sweeps=0  gs pad=None evolve pad=None : wf0 link dims [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  <Sz_0>(end)=-0.500000 max|DMRG-ED|=4.929e-01
sweeps=0  gs pad=8    evolve pad=8   : wf0 link dims [8, 8, 8, 8, 8, 8, 8, 8, 8, 8]  <Sz_0>(end)=-0.500000 max|run-unpadded|=7.772e-16 max|DMRG-ED|=4.929e-01
sweeps=0  gs pad=8    evolve pad=None: wf0 link dims [8, 8, 8, 8, 8, 8, 8, 8, 8, 8]  <Sz_0>(end)=-0.007077 max|run-unpadded|=4.929e-01 max|DMRG-ED|=4.902e-08
sweeps=3  gs pad=None evolve pad=None : wf0 link dims [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  <Sz_0>(end)=-0.007077 max|DMRG-ED|=1.419e-07
sweeps=3  gs pad=8    evolve pad=8   : wf0 link dims [8, 8, 8, 8, 8, 8, 8, 8, 8, 8]  <Sz_0>(end)=-0.007077 max|run-unpadded|=5.058e-10 max|DMRG-ED|=1.420e-07
sweeps=3  gs pad=8    evolve pad=None: wf0 link dims [8, 8, 8, 8, 8, 8, 8, 8, 8, 8]  <Sz_0>(end)=-0.007077 max|run-unpadded|=7.434e-10 max|DMRG-ED|=1.427e-07
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced on its own seed against its own ED,
through both triggers, and through the cached `wf0` (`05_cached_wf0_via_h.py`);
the entangled start (true bonds [2,4,8,5,8,5,8,4,2] at K=16) gives 1.517e-06, the
size of the unpadded run's own one-site error against ED (1.531e-06), so outside
the degenerate product start the effect sits inside the method's own error:

`reviews/pyitensor_P1/common.py`:

```python
# Shared helpers for 02/03: the two starts, own ED, and the two candidate
# replacements for _strip_bond_padding, monkeypatched in-process (the repo is
# not touched).
#   fix1  the hunter's: the same lossless position(n)/position(1) sweep with
#         the gate on the global flag removed
#   fixS  a state-keyed gate: run the sweep on a COPY with padding suspended,
#         and adopt the copy only if some bond dimension shrank; otherwise
#         return psi exactly as it came in
import numpy as np
from scipy.linalg import expm, eigh
from dmrgpy import spinchain, timedependent
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor import chain as pychain
from dmrgpy.pyitensor.mpscontainer import _link_at

ORIG = pychain._strip_bond_padding
N = 10; NT = 40; DT = 0.05; DELTA = 0.7; HZ = 0.1


def bonds(psi):
    return [_link_at(psi, i, i+1).dim for i in range(1, psi.length())]


def _sweep(psi):
    n = psi.length()
    with bk.pad_bonds_suspended():
        if psi.center is None:
            psi.center = n
        psi.position(n)
        psi.position(1)
    return psi


def fix1(psi):
    if psi.length() < 2:
        return psi
    return _sweep(psi)


def fixS(psi):
    if psi.length() < 2:
        return psi
    before = bonds(psi)
    trial = _sweep(psi.copy())
    if any(a < b for a, b in zip(bonds(trial), before)):
        return trial
    return psi


STRIPS = {"now": ORIG, "fix1": fix1, "fixS": fixS}


def build(start, pad_gs, K):
    np.random.seed(5)
    bk.set_pad_bonds(pad_gs)
    sc = spinchain.Spin_Chain([2]*N, itensor_version="python")
    sc.maxm = K
    h0 = 0
    if start == "neel":
        sc.nsweeps = 10; sc.cutoff = 1e-12
        for i in range(N):
            h0 = h0 + (-1)**i*sc.Sz[i]
    else:  # dimerized Heisenberg + 0.3 staggered field
        sc.nsweeps = 12; sc.cutoff = 1e-10
        for i in range(N-1):
            J = 1.0 if i % 2 == 0 else 0.2
            h0 = h0 + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
        for i in range(N):
            h0 = h0 + 0.3*(-1)**i*sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    h = 0
    for i in range(N-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + DELTA*sc.Sz[i]*sc.Sz[i+1]
    for i in range(N):
        h = h + HZ*sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.tevol_method = "TDVP_GSE"
    return sc, wf


def run(start, K, pad_gs, sweeps, strip="now"):
    """Ground state with padding pad_gs, then the flag CLEARED and the
    evolution run (T1). pad_gs=None is the fully unpadded run."""
    pychain._strip_bond_padding = STRIPS[strip]
    try:
        sc, wf = build(start, pad_gs, K)
        entry = bonds(wf.cpp_handle)
        sc.tdvp_gse_sweeps = sweeps
        bk.set_pad_bonds(None)
        out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=NT, dt=DT, wf=wf)
        return np.real(np.asarray(out[1])), entry
    finally:
        bk.set_pad_bonds(None)
        pychain._strip_bond_padding = ORIG


def ed(start):
    sx = np.array([[0, 1], [1, 0]])/2; sy = np.array([[0, -1j], [1j, 0]])/2; sz = np.diag([0.5, -0.5])
    def op(o, i):
        m = np.array([[1.0]])
        for j in range(N):
            m = np.kron(m, o if j == i else np.eye(2))
        return m
    H1 = sum(op(sx, i)@op(sx, i+1) + op(sy, i)@op(sy, i+1) + DELTA*op(sz, i)@op(sz, i+1) for i in range(N-1))
    H1 = H1 + HZ*sum(op(sz, i) for i in range(N))
    if start == "neel":
        bits = "".join("1" if i % 2 == 0 else "0" for i in range(N))
        psi = np.zeros(2**N, dtype=complex); psi[int(bits, 2)] = 1.0
        gap = float("nan")
    else:
        H0 = sum((1.0 if i % 2 == 0 else 0.2)*(op(sx, i)@op(sx, i+1) + op(sy, i)@op(sy, i+1) + op(sz, i)@op(sz, i+1)) for i in range(N-1))
        H0 = H0 + 0.3*sum((-1)**i*op(sz, i) for i in range(N))
        w, V = eigh(H0); psi = V[:, 0].astype(complex); gap = w[1]-w[0]
    U = expm(-1j*DT*H1); Sz0 = op(sz, 0); out = []
    for _ in range(NT):
        out.append(np.vdot(psi, Sz0 @ psi).real); psi = U @ psi
    return np.asarray(out), gap


def report(start, K):
    E, gap = ed(start)
    print("start=%s K=maxm=%d  ED H0 gap=%.4f  ED <Sz_0>(0)=%.6f" % (start, K, gap, E[0]), flush=True)
    for sweeps in (0, 3):
        ref, d0 = run(start, K, None, sweeps)
        now, d1 = run(start, K, K, sweeps)
        print(" sweeps=%d unpadded entry bonds %s max|unpadded-ED|=%.3e" % (sweeps, d0, np.max(np.abs(ref-E))))
        print("          padded   entry bonds %s T1 now : max|run-unpadded|=%.3e max|run-ED|=%.3e"
              % (d1, np.max(np.abs(now-ref)), np.max(np.abs(now-E))), flush=True)
        for s in ("fix1", "fixS"):
            y, _ = run(start, K, K, sweeps, s)
            yu, _ = run(start, K, None, sweeps, s)
            print("          %s: T1 max|run-unpadded|=%.3e | unpadded run vs now max|diff|=%.3e array_equal=%s"
                  % (s, np.max(np.abs(y-ref)), np.max(np.abs(yu-ref)), np.array_equal(yu, ref)), flush=True)
```

`reviews/pyitensor_P1/01_repro_triggers.py`:

```python
# Reviewer repro of pyitensor P1: _strip_bond_padding gated on the global
# flag, not on the state. Own seed (5), own ED (dense scipy expm, written
# here, not dmrgpy's), and three ways a padded state can reach an
# evolution that runs with the flag off:
#   T1  explicit wf= from a padded ground state, flag cleared before evolving
#   (T2, "no wf=, the chain's cached wf0", is not a trigger: set_hamiltonian
#   drops wf0, and evolve_and_measure(mode="DMRG") without wf= then raises
#   AttributeError on None.cpp_handle, see 01_repro_triggers_firstattempt_T2crash.out)
#   T3  flag left ON, but the evolution wrapped in backend.pad_bonds_suspended()
# plus the fully padded run (the finding-16 configuration, the fixed one) and
# the fully unpadded run (the method as asked for).
# n=10 XXZ Delta=0.7 + hz=0.1, Neel start (ground state of a staggered field),
# K=maxm=8, 40 steps of dt=0.05, observable <Sz_0>(t).
import numpy as np
from scipy.linalg import expm
from dmrgpy import spinchain, timedependent
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor.mpscontainer import _link_at

N = 10; K = 8; NT = 40; DT = 0.05; DELTA = 0.7; HZ = 0.1


def build(pad_gs):
    np.random.seed(5)
    bk.set_pad_bonds(pad_gs)
    sc = spinchain.Spin_Chain([2]*N, itensor_version="python")
    sc.maxm = K; sc.nsweeps = 10; sc.cutoff = 1e-12
    h0 = 0
    for i in range(N):
        h0 = h0 + (-1)**i*sc.Sz[i]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs()
    dims = [_link_at(wf.cpp_handle, i, i+1).dim for i in range(1, N)]
    h = 0
    for i in range(N-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + DELTA*sc.Sz[i]*sc.Sz[i+1]
    for i in range(N):
        h = h + HZ*sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.tevol_method = "TDVP_GSE"
    return sc, wf, dims


def run(pad_gs, trigger, sweeps):
    try:
        sc, wf, dims = build(pad_gs)
        sc.tdvp_gse_sweeps = sweeps
        if trigger == "padded":          # finding-16 configuration
            out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=NT, dt=DT, wf=wf, return_wf=True)
        elif trigger == "T1":
            bk.set_pad_bonds(None)
            out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=NT, dt=DT, wf=wf, return_wf=True)
        elif trigger == "T2":
            bk.set_pad_bonds(None)
            out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=NT, dt=DT, return_wf=True)
        elif trigger == "T3":
            with bk.pad_bonds_suspended():
                out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=NT, dt=DT, wf=wf, return_wf=True)
        else:
            raise ValueError(trigger)
        wf_end = out[2].cpp_handle
        dims_end = [_link_at(wf_end, i, i+1).dim for i in range(1, N)]
        return np.real(np.asarray(out[1])), dims, dims_end
    finally:
        bk.set_pad_bonds(None)


# independent ED
sx = np.array([[0, 1], [1, 0]])/2; sy = np.array([[0, -1j], [1j, 0]])/2; sz = np.diag([0.5, -0.5])
def op(o, i):
    m = np.array([[1.0]])
    for j in range(N):
        m = np.kron(m, o if j == i else np.eye(2))
    return m
H = sum(op(sx, i) @ op(sx, i+1) + op(sy, i) @ op(sy, i+1) + DELTA*op(sz, i) @ op(sz, i+1) for i in range(N-1))
H = H + HZ*sum(op(sz, i) for i in range(N))
# Neel with site 0 DOWN: the staggered field (-1)^i Sz_i is minimised by Sz_0=-1/2
bits = "".join("1" if i % 2 == 0 else "0" for i in range(N))   # basis 0=up, 1=down
psi = np.zeros(2**N, dtype=complex); psi[int(bits, 2)] = 1.0
U = expm(-1j*DT*H); Sz0 = op(sz, 0); ED = []
for _ in range(NT):
    ED.append(np.vdot(psi, Sz0 @ psi).real); psi = U @ psi
ED = np.asarray(ED)
print("ED <Sz_0>(0)=%.4f  <Sz_0>(end)=%.6f" % (ED[0], ED[-1]))

for sweeps in (0, 3):
    ref, d0, de0 = run(None, "unpadded", sweeps) if False else (None, None, None)
    # fully unpadded reference: gs unpadded, evolution unpadded
    try:
        sc, wf, d0 = build(None)
        sc.tdvp_gse_sweeps = sweeps
        out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=NT, dt=DT, wf=wf, return_wf=True)
        ref = np.real(np.asarray(out[1]))
        de0 = [_link_at(out[2].cpp_handle, i, i+1).dim for i in range(1, N)]
    finally:
        bk.set_pad_bonds(None)
    print("sweeps=%d unpadded          : entry bonds %s end bonds %s max|run-ED|=%.3e"
          % (sweeps, d0, de0, np.max(np.abs(ref-ED))), flush=True)
    for trig in ("padded", "T1", "T3"):
        y, d, de = run(K, trig, sweeps)
        print("sweeps=%d gs pad=%d %-7s: entry bonds %s end bonds %s max|run-unpadded|=%.3e max|run-ED|=%.3e"
              % (sweeps, K, trig, d, de, np.max(np.abs(y-ref)), np.max(np.abs(y-ED))), flush=True)
```

`reviews/pyitensor_P1/01_repro_triggers.out`:

```
ED <Sz_0>(0)=-0.5000  <Sz_0>(end)=-0.007077
sweeps=0 unpadded          : entry bonds [1, 1, 1, 1, 1, 1, 1, 1, 1] end bonds [1, 1, 1, 1, 1, 1, 1, 1, 1] max|run-ED|=4.929e-01
sweeps=0 gs pad=8 padded : entry bonds [8, 8, 8, 8, 8, 8, 8, 8, 8] end bonds [1, 1, 1, 1, 1, 1, 1, 1, 1] max|run-unpadded|=1.166e-15 max|run-ED|=4.929e-01
sweeps=0 gs pad=8 T1     : entry bonds [8, 8, 8, 8, 8, 8, 8, 8, 8] end bonds [2, 4, 8, 8, 8, 8, 8, 4, 2] max|run-unpadded|=4.929e-01 max|run-ED|=4.902e-08
sweeps=0 gs pad=8 T3     : entry bonds [8, 8, 8, 8, 8, 8, 8, 8, 8] end bonds [2, 4, 8, 8, 8, 8, 8, 4, 2] max|run-unpadded|=4.929e-01 max|run-ED|=4.902e-08
sweeps=3 unpadded          : entry bonds [1, 1, 1, 1, 1, 1, 1, 1, 1] end bonds [2, 4, 8, 8, 8, 8, 8, 4, 2] max|run-ED|=1.393e-07
sweeps=3 gs pad=8 padded : entry bonds [8, 8, 8, 8, 8, 8, 8, 8, 8] end bonds [2, 4, 8, 8, 8, 8, 8, 4, 2] max|run-unpadded|=3.178e-09 max|run-ED|=1.424e-07
sweeps=3 gs pad=8 T1     : entry bonds [8, 8, 8, 8, 8, 8, 8, 8, 8] end bonds [2, 4, 8, 8, 8, 8, 8, 4, 2] max|run-unpadded|=2.540e-09 max|run-ED|=1.418e-07
sweeps=3 gs pad=8 T3     : entry bonds [8, 8, 8, 8, 8, 8, 8, 8, 8] end bonds [2, 4, 8, 8, 8, 8, 8, 4, 2] max|run-unpadded|=2.540e-09 max|run-ED|=1.418e-07
```

Struck: `_strip_bond_padding`'s own docstring as a broken promise (it documents the
gate); the `set_pad_bonds` quotation about the ground state keeping its padding
(it describes a session padded throughout); `quench_tdvp_gse` as affected (9.4e-16,
`04_quench_immune.py`); and "unpadded runs stay `array_equal`" under the hunter's
first fix, which held only on the hunter's seed. New, but a gap in finding 16's fix
rather than its recorded residual (`Chain.tdvp_step` on the TDZ route never strips,
which is justified there). No document describes toggling the flag between the
state and the evolution, and none forbids it; LOW, below finding 16.

**Suggested fix**: the gate must key on the state, with unpadded runs bit-identical.
The hunter's first fix, an ungated strip, breaks bit-identity (7.2e-16 on the Neel
start at sweeps 0, 1.9e-9 at sweeps 3, 7.8e-15 and 3.8e-11 on the entangled one);
marking what `svd()` pads has the right constraint but misses this finding's input,
since a padded ground state is made in `dmrg.py::_apply_local_update` (a direct
`svd()` then `set_A`), so the mark would have to travel with the bond index. What
measured clean, as `fixS` in `common.py`: with the flag off, run the same suspended
`position(n)`/`position(1)` sweep on a copy and adopt it only if some bond shrank,
otherwise return psi untouched; it repairs the toggled case to 1.2e-15 (Neel) and
3.6e-14 (entangled) and is `array_equal` on unpadded runs at sweeps 0 and 3 on both
starts (`02_entangled_fixes.py`, `03_neel_fixes.py`). Its caveat: it keys on exactly
zero singular values, a wider set than "was padded", which the flag-on branch
already strips today. NUMBERS CHANGE only in the padded-then-unpadded configuration.

## Ruled out

Executed and found clean, so the next pass need not redo them.

- **operators**: `get_dagger()` on sums, products and powers, `ISy` squared, `ISy`
  next to `Sx`/`Sz` on one site, complex coefficients, spin-1 and spin-2, mixed
  C/A terms, bosons, Z3 parafermions, the native spinful and the mixed spin-fermion
  chains, checked through <w1|O^dagger|w2> = conj(<w2|O|w1>), whose right-hand side
  never goes through `get_dagger()`, agrees to 2.8e-16 on `"python"`, v3 and v2,
  and daggering twice gives the operator back exactly; `_dagger_name` plus
  `_dagger_phase` cover every non-Hermitian name on every site type except spin-2
  `S2`, the documented vendored typo; about 1000 landed proofs per alphabet over
  seven alphabets (spinless, native spinful, bosons with projectors, spin-1/2,
  spin-1, spin plus fermion), against hand-built Jordan-Wigner matrices, gave zero
  false Hermiticity, anti-Hermiticity, zero or dagger-pair proofs and zero value
  changes under `simplify()`; nothing downstream re-sorts a mixed C/A term with the
  old sign (`vev` against dense, 1.8e-15 on `"python"`, 6.4e-9 on v3); the chain's
  probe gives 72 right verdicts out of 72 on `ISy` operators; a
  Dzyaloshinskii-Moriya chain written through `ISy` is answered correctly by
  `is_hermitian`, `gs_energy()` and `exponential()`; `get_dagger(conjugate=False)`
  has not read its keyword since `79f9fef` (dead, not wrong); the TD second
  evolution receives the right adjoint pair.
- **kpm**: every DMRG route, v2 included, hands the reconstruction exactly the
  calibrated n moments, accelerated or not, n even or odd, self and cross pairs,
  and the moments are exact to 1e-12; `julia_live` computes the same n and a curve
  within 8.4e-12 of v3's; `kpm_n_scale` in (1.5, 2.0, 0, -1, True, `np.int64(2)`,
  `np.float64(2.0)`, `np.bool_(True)`) is rejected before the ground state on every
  Hermitian route, the fallbacks included, and the Kondo second-order DMRG route,
  `fermionchain.get_gr` and `kpm_finite` via `window_chain_kwargs` raise before the
  ground state too; ED `get_distribution` on a non-diagonal X integrates to
  1.000000 with an imaginary part of exactly 0 and per-line weights within 3.8e-6
  of the exact ones; nothing in `src/`, `examples/`, `tests/` or `benchmarks/`
  consumed the old n+2 count, the old ED normalization or a float `kpm_n_scale`;
  the `kpm_emax` cache is keyed on e0 and cleared by `set_conserved_sector` and
  `set_hamiltonian`; and, from finding 3's reviewer, 80 correct 12-site DMRG runs
  across truncation, `kpm_scale`, delta, acceleration and pair shape never put a
  moment above the exact Cauchy-Schwarz bound.
- **realtime**: the direct Fourier sum sits at the same 2.9e-4 trapezoid floor
  against the exact continuous transform on descending, non-uniform, 2-D and
  default grids, under `gaussian` and `parzen` damping, with negative poles and
  frequencies, odd and even nt, and with `predict=True`; TD and TDZ in a conserved
  sector whose ground state is not the global one use the sector's own state and
  origin on ED, `"python"` and v3 (1.2e-6 and 4.2e-4 against the restricted
  Lehmann sum); a misspelled TD keyword raises on v2, v3 and `"python"` on both the
  public and the lower-level route; no TDZ caller forwards a keyword it no longer
  takes.
- **kondo**: the `"python"` members survive the hidden re-solve of finding 13 in
  every run, at the Sz-protected split and at an avoided crossing; v3 at the exact
  crossing gives the right average; the chain's state (not its energy) is restored
  after `n_gs`; the second-order docstring's requirement that `es` start several
  delta below zero holds (0.753 to 0.810 lost against a predicted 0.785);
  `mode="DMRG"` raises on the other typos tried (`Jrho=`, `ngs=`, `sub_mode=`);
  duplicated `es` points and a list `es` give the sorted-grid answer;
  `secondorder_dc.py`'s hunk in `30200a4` is docstring only.
- **pyitensor**: the deflated solver on the previous record's cell agrees with
  dense to 1.0e-11 at k=0 and 3e-15 elsewhere over n=2 to 5, is reproducible run to
  run, returns vectors orthonormal to 4.4e-16, the degenerate-pair subspace and the
  pair sums of Sx and Sz equal dense's, and v3 agrees to 2.3e-12; the Haldane
  triplet, exact and split by 1e-7 to 1e-3, n cutting through it, matches dense to
  6e-15; planted spectra (exact and split pairs and triples, a quintuplet under a
  triplet, a 20-level cluster inside 1e-3, n up to the dimension) match to 5e-13 or
  better; padded against unpadded, the ground-state energy is identical, two-site
  TDVP agrees to 3.6e-12, TEBD to 1.7e-14, `TDVP_GSE` to 3e-15 at sweeps 0, KPM to
  6.6e-15, and continuing from the one-site route's returned state agrees to
  2.5e-9; `pad_bonds_suspended` restores the caller's setting on normal exit, on an
  exception and when nested; the `svd.py` and `backend.py` hunks of `30200a4` are
  comment and docstring only.

## New leads, not reviewed

Observed along the way and not handed to a reviewer of their own, so they are leads
rather than findings.

- `Many_Body_Chain.__init__(**kwargs)` drops every keyword through
  `initialize(**kwargs)`: `Spin_Chain(['1/2']*4, maxm=50).maxm` is 30, and
  `kpm_nscale=3, bogus_key=1` are accepted. The 2026-08 record noted it at its line
  1723 as a reviewer's sharpening ("`maxm=77`/`nsweeps=3` are discarded
  library-wide"); it never got a finding of its own. (2026-09-25: fixed, see `audit_2026_09_25_open_items.md`, item 1.)
- ED clips the ground-state pole just above `kpm_scale=1/2`: for `kpm_scale` in
  (0.5, 0.526], ED's |x| <= 0.95 evaluation cutoff drops the E0 weight (0.2398
  against DMRG's 0.2499 at 0.52). The DMRG counterpart: at exactly 0.5 the 0.99 in
  `kpmdmrg.dynamical_correlator_from_moments`' grid drops a pole at |x0| just below
  1 (a fermion (N_0,N_0) pair integrates to 0.248 against 0.558, the missing
  <N_0>^2 = 0.311), roughly for (0.5, 0.505].
- The elastic weight is clipped on every ground-state-anchored DMRG curve: E0 sits
  at x = -0.9875, 0.0025 inside the -0.99 mask, so about 39 per cent of any elastic
  weight is lost, a 2.2 per cent sum-rule deficit on the finding 4 chain; and the
  default `kpm_scale=0.7` with `kpm_energy_truncate=True` reaches -1.46 against a
  0.75 peak on that 4-site chain with only 1.5 per cent of the weight outside the
  window (below the known-issue file's "92.4 per cent inside means correct"; not
  converged at the default truncation settings). Both are `kpm_energy_truncate`
  territory.
- `get_dynamical_correlator_moments` has no `get_mode()`: on a 2-site v3 chain it
  aborts the interpreter ("Chain::kpm_dynamical_correlator called before
  set_hamiltonian") where `get_dynamical_correlator` falls back to ED.
- `mpsjulialive/dynamics.py::_kpm_dynamical_correlator` still takes `n=1000`,
  `deconvolve=None` and `**kwargs` and reads none of them (by reading).
- `edtk/distribution.py::get_distribution` returns `None` for a `method` other
  than "KPM" or "INV" (by reading).
- `kpm_finite` with `mode="ED"` in `window_chain_kwargs` dies with `AttributeError`
  (the bare `Many_Body_Chain` has no `get_ED_obj`), and its docstring still lists
  `deconvolve` among the forwarded keywords, which now raise.
- ED TD after `set_gs(State)` on a chain whose ED ground state was never solved
  raises `TypeError` from `e0=None` (`edtk/timedependent.py:97`), where before
  `30200a4` it re-solved and ran; ED KPM crashes the same way, and did before.
- `set_hamiltonian(H2, restart=False)` followed by `mode="ED"` TD was not tested;
  by reading ED answers with the stale EDchain's own Hamiltonian while DMRG TD
  re-solves for H2.
- `evolve_and_measure(mode="DMRG")` without `wf=` after a `set_hamiltonian` raises
  `AttributeError: 'NoneType' object has no attribute 'cpp_handle'` instead of
  computing the ground state.
- On `"python"`, `gs_energy(wf0=x)` mutates the caller's `x` in place (`<H>` of
  `x` goes from -0.278 to -1.000); v3 does not.
- `groundstate.best_gs` continues one solve n times rather than taking n
  independent tries (all n return the same `<Sz_tot>`), and `get_gs(best=True,
  **kwargs)` forwards keywords to a `best_gs(sc, n)` that takes none; it has no
  caller.
- `copy.deepcopy` of a solved chain keeps `computed_gs=True` and `wf0` with an
  empty session, which only `dynamics.py:190`'s round trip papers over (finding
  13's reviewer).
- `disentangle_manifold`'s `eig` branch also loses orthonormality on a normal but
  non-Hermitian `ma` with a repeated eigenvalue (a proven anti-Hermitian `1j*Sz0`,
  the unitary Z3 charge `Tau0*Tau1`), where a complex Schur form would give 1e-15;
  the ED twin takes `eigh` unconditionally, the mirror problem for a non-Hermitian
  operator.
- `_deflated_lanczos_run`'s `beta < 1e-12` refusal can never fire, since the
  residual test one line above has already passed whenever beta is that small;
  from a start inside an invariant subspace it returns 0.384615 against an exact
  0. Production starts are seeded Gaussians, so no returned number is affected.
- `SpinBoson_Chain` masks `A`/`Adag` with the literal int 0 at spin sites
  (`bosonchain.py:194`), which degrades under multiplication and has no
  `get_dagger()` (by reading; a crash, not a number).
- `submode="SECTOR"` raises on an empty target sector although its own message says
  the correlator is then identically zero, so a caller looping over channels has to
  special-case it.
- The potential term's sum-rule warning advises "widen es" whenever the held
  weight is off, including an excess, which missing coverage can never cause (a
  sorted, under-resolved uniform grid at dx/delta=10 holds 1.552 of 0.75).
- EX is 7 to 11 per cent low at delta=1e-4 on a non-degenerate ground state three
  delta from the crossing, where KPM is 0.24 per cent low: poles two delta below
  the threshold lose their Lorentzian tail (predicted to 4e-3); by construction,
  not a bug.
- The sorted-grid second-order error of 0.54 on a 6.97 peak on the 3-site v3 KPM
  route (finding 17's reviewer) was not examined.
- `julia_live`'s `evolve_and_measure_tdvp` triggers ITensorMPS's
  `make_inds_match` deprecation warning in `inner(psi,Aop,psi)`; the result is
  right today.

## Statements this hunt overturns

Each of these is written as established and is contradicted by an entry above;
when the corresponding fix lands, the older text should gain a pointer here.

- The previous record's finding 12 Status: "refuses `itensor_version=2` ... (not
  run)" (finding 13) and "the chain's own state is restored exactly" (finding 14);
  the comment at `spinchain.py:412-413`, that the dropped session energy is
  something "nothing reads again while the solver key is current" (findings 13
  and 14).
- The previous record's finding 11: "On `mode="DMRG"` the same typos raise, only
  because the dynamical correlator downstream is strict" (finding 16).
- The previous record's finding 7 Status, its reason for keeping the short-window
  callers on the FFT stage ("the direct stage moved the test's TD peak to the
  window edge", findings 7 and 8), and the two sentences of
  `_fourier_transform_correlator`'s docstring naming `sxt_to_skomega` as the
  switch's one caller and the two short-window callers as getting the direct
  evaluation (finding 7).
- The previous record's finding 5 Status, "the same reach as on `mode="DMRG"`",
  and `documentation.md:4664`'s placement of the check (finding 5).
- The previous record's finding 16 Status, FIXED for `evolve_and_measure_tdvp_gse`,
  `set_pad_bonds`' docstring and the user guide's "the padding being stripped once
  at trajectory entry" (finding 18).
- `docs/audit_2026_09_hole_hunt.md:212`, that `gs_energy_single` re-applies a
  `set_initial_wf()` state through `set_wavefunction`; `set_initial_wf_guess`'s
  docstring; `tests/test_bond_dimension_ramp.py`'s module docstring (finding 12).
- `gs_is_current`'s docstring, that an injected state "is returned unconditionally"
  (finding 11).
- `known_issue_kpm_energy_truncation_window.md`'s "the ED path is unaffected"
  (findings 2 and 4); `CLAUDE.md`'s O2 sentence "the ED route adopted the DMRG
  rescaling so the two share the same x at the same physical energy" and
  `dynamics.py`'s "they now share both", under the flag (finding 4).
- The user guide's KPM promise that the recursion "detects the resulting
  exponential divergence and aborts with an explicit error rather than returning a
  silently wrong spectrum" (findings 2 and 3).
- `documentation.md`'s "both agree on where the dominant spectral weight sits"
  (finding 8).
- The comment at `examples/time_evolution/time_evolution_ABA/main.py:13`
  (finding 9).

## Shared helpers

Modules imported by the scripts above, verbatim. `reviews/kondo_K5/03_size_prefix_and_fix.py`
also imports `old_potentialdc_30200a4parent`, which is
`git show 30200a4^:src/dmrgpy/kondospectrumtk/potentialdc.py` with its one
relative import made absolute, so it is not repeated here.

The three-slot wrapper every later script ran through:

`run3.sh`:

```bash
#!/bin/bash
# Runs one python script under a workflow-wide cap of three concurrent runs.
# Usage, from the folder holding the script:
#   <this file> NN_slug.py [args] 2>&1 | tee NN_slug.out
# It waits (checking every 5 s) until one of three lock slots is free, then runs
# the script with threads pinned and this worktree's src on PYTHONPATH. The
# slot is released when the script exits, including when it is killed.
LOCKDIR=<scratch>/locks
while true; do
  for i in 1 2 3; do
    exec 9>"$LOCKDIR/slot$i.lock"
    if flock -n 9; then
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

`reviews/kpm_A/reviewlib.py`:

```python
"""Shared helpers for the kpm_A review probes (not a probe itself).

exact_moments reproduces what the ED KPM route SHOULD compute from an
exact eigendecomposition: the rescaling edtk/dynamics.py uses
(x_n = (E_n-E0-W/2)/(W*kpm_scale)), the weights w_n = <vj|n><n|vi> with
vi = B|gs>, vj = A^dag|gs>, and mu_k = sum_n w_n T_k(x_n) restricted to
the poles inside [-1,1] (the only ones a Chebyshev expansion can
represent). With every weighted pole inside, these ARE the ED moments up
to roundoff."""
import numpy as np
from dmrgpy import spinchain
from dmrgpy.algebra import kpm

def build(v="python", ns=4, field=0.2):
    np.random.seed(5)
    sc = spinchain.Spin_Chain(["S=1/2"]*ns, itensor_version=v)
    h = 0
    for i in range(ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    if field != 0.0: h = h + field*sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 30
    return sc

def spectrum(sc, A, B):
    ed = sc.get_ED_obj()
    H = np.array(ed.get_hamiltonian().todense())
    e, V = np.linalg.eigh(H); gs = V[:, 0]
    Am = np.array(ed.MO2matrix(A).todense()); Bm = np.array(ed.MO2matrix(B).todense())
    vi = Bm@gs; vj = np.conjugate(Am.T)@gs
    w = np.conjugate(np.conjugate(V.T)@vj)*(np.conjugate(V.T)@vi)
    return e, w, np.linalg.norm(vi)*np.linalg.norm(vj)

def exact_moments(e, w, ks, n):
    W = e[-1]-e[0]
    x = (e-e[0]-W/2.)/(W*ks)
    inside = np.abs(x) <= 1.0
    k = np.arange(n)
    T = np.cos(np.outer(k, np.arccos(np.clip(x[inside], -1, 1))))
    return T@w[inside], x

class Spy:
    """Record the moments ED's route computes; optionally substitute the
    exact in-window moments, so the SAME reconstruction pipeline yields
    the reference curve."""
    def __init__(self):
        self.orig = kpm.get_moments_vivj; self.rec = []; self.sub = None
    def __enter__(self):
        spy = self
        def f(m0, vi, vj, n=100, use_fortran=False):
            mus = spy.orig(m0, vi, vj, n=n, use_fortran=use_fortran)
            spy.rec.append(mus)
            if spy.sub is not None: return spy.sub(n)
            return mus
        kpm.get_moments_vivj = f
        return self
    def __exit__(self, *a):
        kpm.get_moments_vivj = self.orig
```

`reviews/kpm_A2/a2lib.py`:

```python
"""Shared helpers for the kpm_A2 review probes (not a probe itself).

Everything a probe compares against comes from an exact eigendecomposition
of the chain's own Hamiltonian (the ED object's dense matrix), never from
another KPM route:

- weights w_n = <vj|n><n|vi>, vi = B|gs>, vj = A^dag|gs>;
- exact_moments_allpoles: mu_k = sum_n w_n T_k(x_n) over EVERY pole, with
  T_k continued outside [-1,1] (cosh/sign form), at the DMRG run's own
  emin/emax/scale. This is what a faithful Chebyshev recursion must return,
  so DMRG matching it says the MPS machinery is right and the window is not;
- exact_moments_inwindow: the same sum restricted to |x_n| <= 1, i.e. the
  best a Chebyshev expansion on this window can do; pushed through the same
  kpmdmrg.dynamical_correlator_from_moments it is the ideal answer for the
  run's own rescaling and moment count;
- lorentz: the exact spectrum broadened by a Lorentzian of HWHM delta, as a
  second, KPM-free anchor for the peak height.
"""
import numpy as np
from dmrgpy import spinchain, kpmdmrg


def heis(v, ns=4, field=0.2, stagger=0.0, maxm=30, nsweeps=30):
    np.random.seed(5)
    sc = spinchain.Spin_Chain(["S=1/2"]*ns, itensor_version=v)
    h = 0
    for i in range(ns-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    if field != 0.0:
        h = h + field*sc.Sz[0]
    if stagger != 0.0:
        for i in range(ns):
            h = h + stagger*(-1)**i*sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = maxm, nsweeps
    return sc


def exact_weights(sc, A, B):
    ed = sc.get_ED_obj()
    H = np.array(ed.get_hamiltonian().todense())
    e, V = np.linalg.eigh(H)
    gs = V[:, 0]
    Am = np.array(ed.MO2matrix(A).todense())
    Bm = np.array(ed.MO2matrix(B).todense())
    vi = Bm @ gs
    vj = np.conjugate(Am.T) @ gs
    w = np.conjugate(np.conjugate(V.T) @ vj)*(np.conjugate(V.T) @ vi)
    return e, w, np.linalg.norm(vi)*np.linalg.norm(vj)


def xpoles(e, emin, emax, scale):
    return (e - (emin+emax)/2.)*scale


def cheb(k, x):
    """T_k(x) for real x of any magnitude."""
    x = np.asarray(x, dtype=float)
    out = np.empty((len(k), len(x)))
    ins = np.abs(x) <= 1.0
    th = np.arccos(np.clip(x[ins], -1, 1))
    out[:, ins] = np.cos(np.outer(k, th))
    xo = x[~ins]
    if xo.size:
        ch = np.arccosh(np.abs(xo))
        out[:, ~ins] = np.cosh(np.outer(k, ch))*np.sign(xo)[None, :]**k[:, None]
    return out


def exact_moments_allpoles(e, w, emin, emax, scale, n):
    x = xpoles(e, emin, emax, scale)
    return cheb(np.arange(n), x) @ w


def exact_moments_inwindow(e, w, emin, emax, scale, n):
    x = xpoles(e, emin, emax, scale)
    ins = np.abs(x) <= 1.0
    return cheb(np.arange(n), x[ins]) @ w[ins]


def reconstruct(mus, emin, emax, scale, n, es, delta):
    return kpmdmrg.dynamical_correlator_from_moments(
        np.asarray(mus), emin, emax, scale, n, es, delta=delta)[1]


def lorentz(e, w, es, delta):
    om = e - e[0]
    return np.array([np.sum(w.real*delta/np.pi/((x-om)**2+delta**2)) for x in es])


def exact_moments_clamped(e, w, emin, emax, scale, n):
    """All poles, with any pole outside [-1,1] moved onto the nearest window
    edge: the nearest thing to the true spectrum a Chebyshev expansion on
    this window can represent (the displacement is (|x|-1)/scale in energy)."""
    x = np.clip(xpoles(e, emin, emax, scale), -1.0, 1.0)
    return cheb(np.arange(n), x) @ w


def fermions(v, ns=6, t=1.0, mu=0.3, U=0.5, maxm=40, nsweeps=30):
    from dmrgpy import fermionchain
    np.random.seed(5)
    fc = fermionchain.Fermionic_Chain(ns, itensor_version=v)
    h = 0
    for i in range(ns-1):
        h = h - t*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i])
        h = h + U*fc.N[i]*fc.N[i+1]
    for i in range(ns):
        h = h - mu*fc.N[i]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = maxm, nsweeps
    return fc
```

The in-process patches the reviewers of findings 9 and 10 used to check their
suggested fixes against the consumer tests, run outside the repository and not
applied to it:

`reviews/realtime_RB/rb_patch.py`:

```python
"""In-process patch of the hunter's fix (the repository is untouched):
edtk.timedependent.evolution_ABC with its two evolve() calls fed -Hop, so
tdtk's e^{+iHt} becomes e^{-iHt}. RB_PATCH=conj instead conjugates the
returned array (the alternative a reader might reach for). Imported by the
driver scripts, or loaded as a pytest plugin with -p rb_patch."""
import os
import numpy as np
from dmrgpy.edtk import timedependent as tded
from dmrgpy.edtk.tdtk import evolve
from dmrgpy.edtk.edchain import State

_ORIG = tded.evolution_ABC

def evolution_ABC_fixed(self, h, A=None, B=None, C=None, wf=None, nt=100, dt=0.01):
    nt = int(nt)
    Aop = self.get_operator(A); Bop = self.get_operator(B); Cop = self.get_operator(C)
    ts = np.array([dt*ii for ii in range(nt)])
    Hop = self.get_operator(h)
    if wf is None: wf = self.get_gs_array()
    elif type(wf) == State: wf = wf.v
    wfA = Aop@wf; wfC = Cop@wf
    cs = []
    for it in range(nt):
        cs.append(np.conjugate(wfC)@Bop@wfA)
        wfA = evolve(wfA, -Hop, t=dt, dt=dt)   # the only change: -Hop
        wfC = evolve(wfC, -Hop, t=dt, dt=dt)
    return ts, np.array(cs)

def evolution_ABC_conj(self, h, **kw):
    ts, cs = _ORIG(self, h, **kw)
    return ts, np.conjugate(cs)

MODE = os.environ.get("RB_PATCH", "sign")
if MODE == "sign":
    tded.evolution_ABC = evolution_ABC_fixed
elif MODE == "conj":
    tded.evolution_ABC = evolution_ABC_conj
elif MODE == "none":
    pass
print("[rb_patch] evolution_ABC patch mode:", MODE)

_CALLS = [0]
_fixed_inner = evolution_ABC_fixed
def _counting(*a, **k):
    _CALLS[0] += 1
    return _fixed_inner(*a, **k)
if MODE == "sign":
    tded.evolution_ABC = _counting

def pytest_unconfigure(config):
    print("\n[rb_patch] patched evolution_ABC was called %d times" % _CALLS[0])
```

`reviews/realtime_RB2/rb2_fixplugin.py`:

```python
"""pytest plugin: the suggested R-B2 fix applied in-process (repository
untouched). timedependent.evolve_and_measure_dmrg's return is conjugated
back, i.e. the call returns the raw session list <psi(t)|O|psi(t)>.
Counts calls and records the largest |Im| of any returned array, which is
exactly the size of what the fix flips."""
import numpy as np
from dmrgpy import timedependent as td

_ORIG = td.evolve_and_measure_dmrg
_STATS = {"calls": 0, "max_imag": 0.0}

def _fixed(self, *a, **k):
    out = _ORIG(self, *a, **k)
    cs = np.conj(np.asarray(out[1]))
    _STATS["calls"] += 1
    _STATS["max_imag"] = max(_STATS["max_imag"], float(np.max(np.abs(cs.imag))) if cs.size else 0.0)
    return (out[0], cs) + tuple(out[2:])

td.evolve_and_measure_dmrg = _fixed
print("[rb2_fixplugin] evolve_and_measure_dmrg return un-conjugated")

def pytest_unconfigure(config):
    print("\n[rb2_fixplugin] patched evolve_and_measure_dmrg called %d times, "
          "largest |Im| returned = %.2e" % (_STATS["calls"], _STATS["max_imag"]))
```
