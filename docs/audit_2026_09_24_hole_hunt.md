# Audit, 2026-09-24: five-lens hole hunt

Findings from a five-lens automated audit of the Python layer, run 2026-09-24 on
`765b537` (clean tree, both compiled extensions current, `make -n pybind`
reporting nothing to do in `mpscpp2` and `mpscpp3`), scoped to the eleven commits
that landed after the 2026-09 audit's own fixes, `1c2606d..765b537`. Each lens
hunted one class of problem; every finding below carries a repro that was
actually executed and its verbatim output, and every one was then handed to an
independent reviewer whose brief was to *refute* it. None was refuted outright,
and every one came back narrowed in at least one sub-claim, so each entry keeps
its struck sub-claims visible next to what survived.

This file is the evidence, not a task list, the same convention as
`audit_2026_08_hole_hunt.md` and `audit_2026_09_hole_hunt.md`: it records what
was observed and how to reproduce it, so a fix (or a decision that the behaviour
is intended after all) does not have to re-derive any of it. Fixed entries gain
a `**Status**` line under their classification line and are kept rather than
deleted, since the repro doubles as the regression check. Several entries
overturn a statement an earlier record or the documentation makes; those are
collected in "Statements this hunt overturns" at the end, so that the older
records can be pointed at this one when the fixes land.

## The five lenses

| Lens | Brief |
|---|---|
| `kpm-calibration` | Does `delta` mean the same broadening on every KPM route after the O2 calibration (`765b537`): same moment count, same window, same curve, and does every reader of `kpm_n_scale` go through the calibration? |
| `td-convention` | Is the two-sided TD/TDZ Lehmann density of O1 (`765b537`) right for every pair shape, every backend and every integrator, and is the operator pair read in the same order everywhere? |
| `canonical-form` | Can the canonical rewrite of `MultiOperator` terms (`593b394`, `ea9ef9b`) prove an operator Hermitian, zero or a dagger pair when it is not, or change an operator's value, and what does each consumer do on "not proven"? |
| `infinite-chain` | The lazy transfer chain of finding #20, the bond-candidate-first fixed points on `itensor_version=3` (`1fb32ec`), JAX iDMRG (`ec14a23`) and the deflated excitation solver of ROADMAP item 10. |
| `recent-misc` | The Kondo tunnelling spectrum (`b75003d`), `set_pad_bonds` no longer padding the Hamiltonian MPO (`2d86cff`), and the device-residency changes in `pyitensor/gse.py` and `tebd.py`. |

## Scope

Out of scope by construction: vendored ITensor (`mpscpp2/ITensor/`,
`mpscpp3/ITensor/`); the legacy bugs `CLAUDE.md` says are deliberately reproduced
(`evoloperator`'s z^3/6 term on `H2`, the `"moise"` key, the unreachable
`"tevol_fit_td"` branch); the open `docs/known_issue_*.md` items, in particular
the `kpm_energy_truncate` window problem; anything already in either earlier
record, including the 2026-09 record's "Demoted to notes" items; and gaps
`ROADMAP.md` marks as absent. The brief also excluded the two O2 residuals as
such (the `sqrt(1-x^2)` width profile and the Jackson line shape) and TDZ's
"contour residual"; finding 7 shows that the second exclusion rested on a
misdiagnosis, which is why it is reported anyway.

`itensor_version="julia_live"` was in scope for exactly two probes, since the
juliacall JIT cost dominates any lens that touches it: whether its KPM line
computes the same moment count and curve as v3 (it does, n=52 and a curve equal
to v3's to 4.2e-15, carrying finding 4 as v3 does), and whether it has a TD/TDZ
route left on the old convention (TDZ agrees with v3 to 4.1e-11 on an
imaginary-weight pair, and TD through `get_dynamical_correlator` raises
`NotImplementedError`, visibly). Everything else excludes it.

As a baseline, the window's own tests were run on the frame before any lens
reported: `pytest tests -k "not julia_live and (kondo or kpm or dynamical or
correlator or canonical or multioperator or infinite or vumps or idmrg or
pad_bonds or gpu or excitation or audit_2026_09)"` gave 832 passed, 2 skipped,
489 deselected in 36 minutes, so every finding below is a hole the suite passes
through rather than code already known to be red.

Every repro was run as

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLBACKEND=Agg PYTHONPATH=<repo>/src python3 <script> [args]
```

from the folder holding the script, with at most two processes at a time on a
14-core workstation. The scripts and their outputs are spliced into this file
from disk rather than retyped, because the scratch folders they ran in do not
survive the session; `<repo>` stands for the checkout and `<scratch>` for the
scratch folder, the only edit made to either. The outputs drop two kinds of
line and nothing else: the repeated notice "ITensor v3's two-site DMRG can't
handle a chain this short ..., using default ED routines" that `mode.py` prints
on every call below three sites, NH-DMRG's banner lines, the per-frequency
"CVM in E = ... CG iterations" progress lines, the missing-extension
`UserWarning` a Python-only archive of an older commit prints, and shell timing
lines. Helper modules that
several scripts import are collected once at the end, in "Shared helpers". Most
reviewers returned their reports as text, since the harness refuses report
files from subagents, so the reviewer paragraphs below are condensed from those
reports while the scripts and outputs are theirs verbatim.

## The findings at a glance

| # | Finding | Severity | Numbers change on fix | Fix cluster |
|---|---|---|---|---|
| 1 | `canonical.py` sorts mixed `A`/`C` terms with the wrong sign: false Hermiticity, zero and dagger-pair proofs | MEDIUM | only for operators mixing the two representations | operators |
| 2 | `get_dagger()` treats `ISy` as Hermitian: KPM, EX, TD, TDZ return minus the correlator | LOW | only with `ISy` in the first operator | operators |
| 3 | `get_distribution(mode="ED")` has total weight 1/scale, and an imaginary part copied from the real one | MEDIUM | yes, by exactly `scale` | kpm |
| 4 | every DMRG KPM route reconstructs from n+2 moments where ED and the calibration use n | LOW | yes, by about 2/n | kpm |
| 5 | a non-integer `kpm_n_scale` is silently rounded down on the Python routes; C++ and Python disagree at or below 0 | LOW | no (validation) | kpm |
| 6 | the documentation never says a ground-state correlator's weight sits where the calibrated width is 0.70 of the requested one | LOW | no (documentation) | kpm |
| 7 | TDZ's recorded "contour residual" is an FFT-grid interpolation error of 14 per cent, shared with TD `predict=False` | MEDIUM | yes, towards exact | real-time |
| 8 | `mode="ED"` TD on a degenerate ground state builds its two halves on two different random states | MEDIUM | only on a degenerate ground state | real-time |
| 9 | since `765b537` TDZ accepts any misspelled keyword and ignores it, and its symbolic-operator check is gone | LOW to MEDIUM | no (validation) | real-time |
| 10 | since `765b537` the lower-level correlator route drops `i=`/`j=` for TD and TDZ | LOW | no (public route unaffected) | real-time |
| 11 | `get_kondo_spectrum(mode="ED")` swallows every unknown keyword: a misspelled `Jrho_s` removes the Kondo peak | MEDIUM | no (validation) | kondo |
| 12 | `get_kondo_spectrum(mode="DMRG")` at a degenerate ground state returns one arbitrary member's spectrum, where `mode="ED"` averages | MEDIUM | only where a caller opts in to the average | kondo |
| 13 | the potential-interference term weights a non-uniform `es` with its first spacing | LOW to MEDIUM | only on a non-uniform `es` | kondo |
| 14 | the potential-term example violates its own `es` coverage requirement and blames the resulting 10 per cent on KPM | LOW | example only | kondo |
| 15 | the `"python"` excitation ansatz can drop one copy of a degenerate `H_eff(k)` level at n>=2 on the default path | MEDIUM | yes, towards dense | pyitensor |
| 16 | `set_pad_bonds` silently turns `TDVP_GSE` into a different algorithm, and `svd.py` says the padding is inert | LOW | yes, back onto unpadded | pyitensor |

The fix clusters are file-disjoint: operators (`multioperatortk/canonical.py`,
`multioperator.py`), kpm (`algebra/kpm.py`, `kpmdmrg.py`,
`edtk/distribution.py`, `manybodychain.py`'s ED branch,
`mpsjulialive/dynamics.py`), real-time (`timedependent.py`, `tdz.py`,
`edtk/timedependent.py`, the TD branch of `edtk/dynamics.py`), kondo
(`kondospectrumtk/`, `spinchain.py`, `examples/kondo/`) and pyitensor
(`pyitensor/gse.py`, `chain.py`, `backend.py`, `svd.py`,
`idmrg_excitations.py`). Every fix recommended below is Python-side, so none
needs a rebuild of either extension; three of them touch
`mpsjulialive/dynamics.py` or an entry point it calls (findings 4, 5 and 9), so
the `julia_live` tests belong to that cluster's finish.

## Findings

### 1. `canonical.py` grades the post-Jordan-Wigner names `A`/`Adag` as even against `C`/`Cdag`, so a mixed term is sorted with the wrong sign: an exactly anti-Hermitian operator is proven Hermitian, an operator of norm 4 is proven zero, and TD returns a correlator 81 per cent off its peak

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `canonical-form`

**Status**: FIXED -- `canonical.py` leaves a term spelled as written, with no reordering and no identity dropping, when it names both a pre-Jordan-Wigner `C`-type name (`_ODD`) and a bare `A`-type name (`_LOCAL_LADDER`: `A`, `Adag`, `Aup`, `Adagup`, `Adn`, `Adagdn`), the reviewer's narrow fix; the `A`-type names keep parity 0, so the boson-hopping and Jordan-Wigner-transformed proofs are kept, and `all_names_known`/`is_dagger_pair` needed no change of their own (adding the rule there would have made `is_dagger_pair(B.get_dagger(), B)` False). `P.is_hermitian()` and `fc.is_hermitian(P)` go True to False, `X.is_zero()` and `is_zero_operator(X)` True to False, `simplify()` keeps the value of a mixed term, and `is_dagger_pair` on the TD pair goes True to False while `is_dagger_pair(B.get_dagger(), B)` stays True. The cost is one-sided and safe: an exactly zero mixed operator such as `C0*A1 - A1*C0` is no longer proven zero, and no dmrgpy chain class puts `C`-type names next to boson ladder operators, so no ordinary Hamiltonian loses its proof. Pinned by `tests/test_audit_2026_09_24_operators.py::test_mixed_anti_hermitian_operator_is_not_proven_hermitian`, `::test_mixed_commutator_of_norm_four_is_not_proven_zero`, `::test_simplify_keeps_the_value_of_a_mixed_term`, `::test_td_pair_mixing_representations_is_not_a_dagger_pair`, `::test_single_representation_proofs_are_kept` (the guard against a fix that is too broad) and `tests/test_multioperator_canonical.py::test_a_term_mixing_c_and_a_names_is_left_as_written`; with the fix reverted by a scratch plugin every new test except the guard fails. NUMBERS CHANGE, only for operators that mix the two representations, from wrong to right: `<w1|simplify(Cdag3 Adag0)|w2>` on a 4-site `Fermionic_Chain` (`"python"`, `np.random.seed(3)`) goes from (-0.020637+0.145197j) to (0.020637-0.145197j), the literal value; the TD pair on the seeded 4-site complex-hopping chain (mu=1.5, delta=0.4, dt=0.1) goes from max|y-exact| 3.616e-02 with one evolution to 3.063e-05 with two, on `"python"` and v3; `gs_energy()` on hop + 1.5 P (L=4) now returns -1.246695 on both backends, the exact lowest real part, where it returned run-dependent values that are not eigenvalues.

**Where**: `src/dmrgpy/multioperatortk/canonical.py` (`_EVEN`, which lists
`A`/`Adag`/`Aup`/`Adagup`/`Adn`/`Adagdn`; the comment at lines 72 to 75 calling
them "plain local matrices by construction, hence even";
`_canonical_signature`). Consumers of the wrong proof:
`MultiOperator.simplify()`/`is_zero()`/`is_hermitian()`;
`Many_Body_Chain.is_hermitian`, reached from `gs_energy`,
`gs_energy_generalized`, `mpsalgebra.exponential`'s gate, `dynamics.py:167` and
`:250`, `excited.py:100`, `sectordc.py:118` and `degeneracy.py:9`;
`operator_norm`/`is_zero_operator`, which feed `cvm.py:70` and
`nonhermitian/dynamics.py:28`; `canonical.is_dagger_pair`, which decides
`timedependent.py`'s one-evolution shortcut and `distribution.py:59`.

The canonical form of `593b394` lists the six post-transform local names as
even. That is right against each other, against `F`, `N` and the spin names, and
wrong against a pre-transform `C`-type name, because dmrgpy writes
C_j = F_0 ... F_{j-1} A_j (`multioperatortk/jordanwigner.py::C`, and the same
string in `jordanwigner_spinful.one_fermion`): an `A`-type factor on site i
anticommutes with a `C`-type factor on site j>i and commutes with one on j<i,
meaning that the relation between the two is not a parity at all but depends on
which of them sits on the lower site. The trigger is a `C`-type factor on a
higher site written to the left of an `A`-type factor on a lower site; the
stable site sort then drops one minus sign for every such pair, and the value of
`simplify()` of the term is exactly (-1)^k times the true value, with k the
number of such pairs. A term written in one representation only (pure
`C`/`Cdag`/`N`/`F`, pure `A`/`Adag`/`N`/`F`, every Jordan-Wigner-transformed
operator, every boson operator) keeps its value.

It survived every existing test for three reasons. `mode="ED"` has no
`A`/`Adag` names on a fermionic chain (it prints "Unrecognised operator Adag"
and dies on a bare `raise`), so every ED-anchored check of the canonical form,
the commit's own 800 random products included, cannot contain them. Nothing
internal mixes the two representations inside one term: `to_terms()` and
`write()` run `jordan_wigner` and never canonicalize afterwards, and
`operatornames.recognize()` has no composite name pairing the two. And the
sympy path that `593b394` replaced never reordered anything, so it answered
"not zero", "not Hermitian" and "not a dagger pair" here, correct by accident.
The way in is public but unadvertised: `fc.A`/`fc.Adag` are attributes of
`Fermionic_Chain`, `Majorana_Chain` keeps them after deleting `C`/`Cdag`/`N`
while its documented `G` operators are `C`-built, the spinful classes inherit
them, and the native spinful chain reaches them through
`get_operator("Aup", i)`. No documentation presents them as fermionic operators;
the user guide names them only in the TEBD backend note.

**Expected**: an operator that mixes the two representations keeps its value
under `simplify()`, and a mixed term is either sorted with its
position-dependent sign or left as written and refused by the proofs, the way
parafermion names already are.

Repro, from the hunter, anchored on dense Jordan-Wigner matrices built from
scratch with numpy:

```python
"""canonical.py lists the post-Jordan-Wigner local names A/Adag (and
Aup/Adagup/Adn/Adagdn) as EVEN, i.e. commuting with a fermion on another
site. On every DMRG backend C_j = F_0...F_{j-1} A_j, and A_i anticommutes
with F_i, so A_i ANTIcommutes with C_j for j>i (and commutes for j<i).
Anchor: dense matrices built from scratch in numpy, no dmrgpy code."""
import numpy as np
from functools import reduce
from dmrgpy import fermionchain
from dmrgpy.multioperatortk import canonical

L = 4
a = np.array([[0,1],[0,0]],dtype=complex); F = np.diag([1.,-1.]).astype(complex)
I2 = np.eye(2,dtype=complex)
def site(op,i): return reduce(np.kron,[op if k==i else I2 for k in range(L)])
Am = [site(a,i) for i in range(L)]; Ad = [x.conj().T for x in Am]
Fm = [site(F,i) for i in range(L)]
Cm = [reduce(np.matmul,[Fm[k] for k in range(j)]+[Am[j]]) for j in range(L)]
Cd = [x.conj().T for x in Cm]
Hhop = sum([Cd[i]@Cm[i+1]+Cd[i+1]@Cm[i] for i in range(L-1)])
print("anchor: ||A0 C1 + C1 A0|| = %.3f   ||A0 C1 - C1 A0|| = %.3f"
      %(np.linalg.norm(Am[0]@Cm[1]+Cm[1]@Am[0]),np.linalg.norm(Am[0]@Cm[1]-Cm[1]@Am[0])))
Pm = Ad[0]@Cm[1] + Am[0]@Cd[1]
print("anchor: P = Adag0 C1 + A0 Cdag1 has ||P-P^dag|| = %.4f, ||P+P^dag|| = %.4f (anti-Hermitian)"
      %(np.linalg.norm(Pm-Pm.conj().T),np.linalg.norm(Pm+Pm.conj().T)))

for v in (3,"python"):
    print("=== itensor_version =",repr(v))
    fc = fermionchain.Fermionic_Chain(L,itensor_version=v)
    h0 = 0
    for i in range(L-1): h0 = h0 + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    fc.set_hamiltonian(h0); fc.maxm = 20; fc.nsweeps = 20
    P = fc.Adag[0]*fc.C[1] + fc.A[0]*fc.Cdag[1]
    X = fc.C[1]*fc.A[0] - fc.A[0]*fc.C[1]
    print("P.is_hermitian()   =",P.is_hermitian(),"   fc.is_hermitian(P) =",fc.is_hermitian(P))
    print("X = C1 A0 - A0 C1: X.is_zero() =",X.is_zero(),
          "  exact ||X||_F = %.3f"%np.linalg.norm(Cm[1]@Am[0]-Am[0]@Cm[1]))
    np.random.seed(3)
    w1 = fc.random_mps(); w2 = fc.random_mps()
    T = fc.Cdag[3]*fc.Adag[0]
    print("<w1|Cdag3 Adag0|w2> = %s   <w1|simplify(Cdag3 Adag0)|w2> = %s"
          %(np.round(fc.aMb(w1,T,w2),6),np.round(fc.aMb(w1,T.simplify(),w2),6)))
    print("operator_norm(X) = %.4f   operator_norm(X,simplify=False) = %.4f   is_zero_operator(X) = %s"
          %(fc.operator_norm(X,ntries=3),fc.operator_norm(X,ntries=3,simplify=False),fc.is_zero_operator(X)))
    for g in (0.8,1.5,3.0):
        ev = np.linalg.eigvals(Hhop+g*Pm); ev = ev[np.argsort(ev.real)]
        fc1 = fermionchain.Fermionic_Chain(L,itensor_version=v)
        fc1.set_hamiltonian(h0+g*P); fc1.maxm = 20; fc1.nsweeps = 20
        e_now = fc1.gs_energy()
        # the route the chain takes when the symbolic proof is refused
        orig = canonical.is_hermitian
        canonical.is_hermitian = lambda MO: False
        try:
            fc2 = fermionchain.Fermionic_Chain(L,itensor_version=v)
            fc2.set_hamiltonian(h0+g*P); fc2.maxm = 20; fc2.nsweeps = 20
            herm2 = fc2.is_hermitian(fc2.hamiltonian)
            e_probe = fc2.gs_energy()
        finally:
            canonical.is_hermitian = orig
        print("g=%.1f  gs_energy() = %-12s exact lowest-Re eigenvalue = %-22s | proof refused: is_hermitian=%s gs_energy() = %s"
              %(g,np.round(e_now,6),np.round(ev[0],6),herm2,np.round(e_probe,6)))
```

Observed (the NH-DMRG banner lines are dropped):

```
anchor: ||A0 C1 + C1 A0|| = 0.000   ||A0 C1 - C1 A0|| = 4.000
anchor: P = Adag0 C1 + A0 Cdag1 has ||P-P^dag|| = 5.6569, ||P+P^dag|| = 0.0000 (anti-Hermitian)
=== itensor_version = 3
P.is_hermitian()   = True    fc.is_hermitian(P) = True
X = C1 A0 - A0 C1: X.is_zero() = True   exact ||X||_F = 4.000
<w1|Cdag3 Adag0|w2> = (0.124905+0.036017j)   <w1|simplify(Cdag3 Adag0)|w2> = (-0.124905-0.036017j)
operator_norm(X) = 0.0000   operator_norm(X,simplify=False) = 1.1715   is_zero_operator(X) = True
g=0.8  gs_energy() = -1.831721    exact lowest-Re eigenvalue = (-1.886796-0j)         | proof refused: is_hermitian=False gs_energy() = (-1.886796+0j)
g=1.5  gs_energy() = -1.196506    exact lowest-Re eigenvalue = (-1.246695-0j)         | proof refused: is_hermitian=False gs_energy() = (-1.246695+0.896799j)
g=3.0  gs_energy() = -1.885963    exact lowest-Re eigenvalue = (-1.059767+2.668915j)  | proof refused: is_hermitian=False gs_energy() = (-1.059767+2.668915j)
=== itensor_version = 'python'
P.is_hermitian()   = True    fc.is_hermitian(P) = True
X = C1 A0 - A0 C1: X.is_zero() = True   exact ||X||_F = 4.000
<w1|Cdag3 Adag0|w2> = (0.020637-0.145197j)   <w1|simplify(Cdag3 Adag0)|w2> = (-0.020637+0.145197j)
operator_norm(X) = 0.0000   operator_norm(X,simplify=False) = 0.7896   is_zero_operator(X) = True
g=0.8  gs_energy() = -1.886796    exact lowest-Re eigenvalue = (-1.886796-0j)         | proof refused: is_hermitian=False gs_energy() = (-1.886796+0j)
g=1.5  gs_energy() = -2.675159    exact lowest-Re eigenvalue = (-1.246695-0j)         | proof refused: is_hermitian=False gs_energy() = (-1.246695+0.896799j)
g=3.0  gs_energy() = -5.759571    exact lowest-Re eigenvalue = (-1.059767+2.668915j)  | proof refused: is_hermitian=False gs_energy() = (-1.059767+0j)
```

`spectrum_check.py` builds the same matrices and prints the degenerate lowest
real parts, which is why the "proof refused" column may land on either member:

```
g=0.8 lowest six by real part: [-1.886796-0.j -1.481915-0.j -1.481915+0.j -1.077033+0.j -0.404882-0.j
 -0.404882+0.j]
       Hermitian part (=hopping) ground state: -2.236068
g=1.5 lowest six by real part: [-1.246695-0.896799j -1.246695-0.j       -1.246695+0.j
 -1.246695+0.896799j  0.      -0.896799j -0.      -0.896799j]
       Hermitian part (=hopping) ground state: -2.236068
g=3.0 lowest six by real part: [-1.059767-2.668915j -1.059767-0.j       -1.059767+0.j
 -1.059767+2.668915j -0.      -2.668915j  0.      -2.668915j]
       Hermitian part (=hopping) ground state: -2.236068
```

The TD consequence, from the `td-convention` lens: `is_dagger_pair` reads the
same grading, so it proves the pair B = C_2 C_3 + A_0 C_1,
A = Cdag_3 Cdag_2 + Adag_0 Cdag_1 to be a dagger pair although the two differ by
the sign of one term, and TD takes the one-evolution `Re F` shortcut on complex
Lehmann weights:

```python
"""TD consequence of the canonical-form lens's A/Adag parity finding:
is_dagger_pair proves a pair that is not a dagger pair, and the
one-evolution Re F shortcut then runs on complex Lehmann weights.

B = C_2 C_3 + A_0 C_1, A = Cdag_3 Cdag_2 + Adag_0 Cdag_1. The literal
dagger of B is Cdag_3 Cdag_2 + Cdag_1 Adag_0, and Adag_0 Cdag_1 =
-Cdag_1 Adag_0 because C_1 carries a Jordan-Wigner string through site 0
while A_0 is the bare local operator, so A != B^dagger."""
import numpy as np
from dmrgpy import fermionchain, timedependent
from dmrgpy.multioperatortk import canonical

n, seed, mu = 4, 3, 1.5
ES = np.linspace(-1.0, 8.0, 46)
DELTA, DT = 0.4, 0.1
rng = np.random.RandomState(seed)
t = rng.random((n, n)) + 1j * rng.random((n, n)); t = t + t.conj().T

# exact reference by hand: a_j local, F = diag(1,-1), C_j = F_0..F_{j-1} a_j
a = np.array([[0, 1], [0, 0]], dtype=complex)   # |occ> -> |emp>, basis (emp, occ)
F = np.diag([1.0, -1.0]).astype(complex); I2 = np.eye(2, dtype=complex)
def kron_list(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops: out = np.kron(out, o)
    return out
def loc(o, j): return kron_list([o if k == j else I2 for k in range(n)])
def Cj(j): return kron_list([F if k < j else (a if k == j else I2) for k in range(n)])
Cm = [Cj(j) for j in range(n)]; Cd = [c.conj().T for c in Cm]
Nm = [Cd[j] @ Cm[j] for j in range(n)]
Am = [loc(a, j) for j in range(n)]; Ad = [x.conj().T for x in Am]
H = sum(t[i, j] * Cd[i] @ Cm[j] for i in range(n) for j in range(n))
H = H + sum(0.8 * Nm[i] @ Nm[i + 1] for i in range(n - 1)) - mu * sum(Nm)
e, U = np.linalg.eigh(H)
Bm = Cm[2] @ Cm[3] + Am[0] @ Cm[1]
Aex = Cd[3] @ Cd[2] + Ad[0] @ Cd[1]
print("||A - B^dag|| (matrices) = %.3f,  ||Adag0 Cdag1 + Cdag1 Adag0|| = %.1e"
      % (np.linalg.norm(Aex - Bm.conj().T), np.linalg.norm(Ad[0] @ Cd[1] + Cd[1] @ Ad[0])))
Uh = U.conj().T
M = (Uh @ Aex @ U)[0, :] * (Uh @ Bm @ U)[:, 0]; D = e - e[0]
ref = np.array([np.sum(M * (DELTA/np.pi) / ((w - D)**2 + DELTA**2)) for w in ES])
g0 = U[:, 0]
print("exact <N> in GS = %.6f" % (g0.conj() @ sum(Nm) @ g0).real)
print("exact E0 = %.10f  max|Im M_n| = %.3e  peak |C| = %.4f"
      % (e[0], np.max(np.abs(M.imag)), np.max(np.abs(ref))))

calls = []
orig = timedependent.evolution_DC
def counted(*a_, **k): calls.append(1); return orig(*a_, **k)
timedependent.evolution_DC = counted
orig_sa = timedependent._pair_is_self_adjoint

for version in ("python", 3):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    h = 0
    for i in range(n):
        for j in range(n): h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1): h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    for i in range(n): h = h - mu * fc.N[i]
    fc.set_hamiltonian(h); fc.maxm, fc.nsweeps = 30, 12
    B = fc.C[2] * fc.C[3] + fc.A[0] * fc.C[1]
    A = fc.Cdag[3] * fc.Cdag[2] + fc.Adag[0] * fc.Cdag[1]
    print("\nitensor_version=%r gs_energy = %.10f   is_dagger_pair(A,B) = %s"
          % (version, fc.gs_energy(), canonical.is_dagger_pair(A, B)))
    for label, force in (("default (shortcut decides)", None),
                         ("shortcut disabled (two runs)", False)):
        timedependent._pair_is_self_adjoint = orig_sa if force is None else (lambda p: False)
        del calls[:]
        _x, y = fc.get_dynamical_correlator(mode="DMRG", submode="TD", name=[A, B],
                                            es=ES, delta=DELTA, dt=DT)
        y = np.asarray(y)
        print("  %-30s runs=%d  max|y-exact| = %.3e (%.0f%% of peak)  max|Im y| = %.2e"
              % (label, len(calls), np.max(np.abs(y - ref)),
                 100*np.max(np.abs(y - ref))/np.max(np.abs(ref)), np.max(np.abs(y.imag))))
    timedependent._pair_is_self_adjoint = orig_sa
```

```
||A - B^dag|| (matrices) = 4.000,  ||Adag0 Cdag1 + Cdag1 Adag0|| = 0.0e+00
exact <N> in GS = 2.000000
exact E0 = -4.0492958122  max|Im M_n| = 4.578e-02  peak |C| = 0.0447

itensor_version='python' gs_energy = -4.0492958122   is_dagger_pair(A,B) = True
  default (shortcut decides)     runs=1  max|y-exact| = 3.616e-02 (81% of peak)  max|Im y| = 0.00e+00
  shortcut disabled (two runs)   runs=2  max|y-exact| = 3.063e-05 (0% of peak)  max|Im y| = 3.59e-02

itensor_version=3 gs_energy = -4.0492958122   is_dagger_pair(A,B) = True
  default (shortcut decides)     runs=1  max|y-exact| = 3.616e-02 (81% of peak)  max|Im y| = 0.00e+00
  shortcut disabled (two runs)   runs=2  max|y-exact| = 3.063e-05 (0% of peak)  max|Im y| = 3.59e-02
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced on `itensor_version=3` and on
`"python"`, with every boolean and every `"python"` digit identical to the
hunter's (the v3 energies move from run to run, see below), and the anchor
confirmed on the backend itself without numpy: through `aMb` with
`simplify=False`, which goes `to_terms()` to `jordan_wigner` to the MPO with no
canonicalization anywhere, `A0 C1 + C1 A0` and `A1 C0 - C0 A1` are exactly zero
while `A0 C1 - C1 A0` and `A1 C0 + C0 A1` are not, and the canonical form reports
exactly the reverse for the j>i pair. The reviewer's sharper scan predicts the
sign of `simplify()` exactly, for 80 of 80 random mixed terms per backend. Four
narrowings. The hunter's prose had the trigger's ordering reversed ("an A-type
factor written left of a C-type one on a higher site"); its code and every
example it printed had it right, and the statement above is the corrected one.
The reach is as stated above and no wider: 0 of 160 single-representation terms
per backend change value. The `operator_norm` witnesses the hunter quoted (1.17
on v3, 0.79 on `"python"`) are single random-witness values that move between
runs; what the record should hold is the exact zero against a nonzero witness.
And `gs_energy()` returning numbers that are not eigenvalues is struck as a
property of this defect and kept as the demonstration that the false proof
reaches a dispatch: forcing `canonical.is_hermitian` to True on the pure-`C`
spelling Q = Cdag_0 C_1 - Cdag_1 C_0 of the very same operator gives the same
kind of garbage (v3 -1.034971, -0.889513, -1.304455; `"python"` -2.407147,
-2.265325, -2.535319, against an exact -1.246695 +- 0.896799i), so any false
Hermiticity proof does this. It also holds only where the spectrum is complex
(at g=0.8 the spectrum is real and `"python"` returns the exact -1.886796), and
"changes from run to run" holds on both backends when unseeded (the hunter's
`"python"` digits repeat only because it seeds its start). The TD consequence
reproduces to every printed digit. The reviewer's two probes:

```python
"""Reviewer probe 1: backend-only anchor, no numpy matrices.
(a) Does the backend itself make A_i anticommute with C_j for j>i and
    commute for j<i?  <w1|T|w2> with simplify=False, T going straight
    through to_terms() -> jordan_wigner -> backend.
(b) What does the canonical form say about the same operators?
(c) Sharp version of the hunter's scan: for random mixed terms, predict
    <simplify(T)> = (-1)^k <T>, k = number of (C-type written first on
    site j, A-type written later on site i<j) pairs, and check it exactly.
(d) Pure-representation scans: pure C/Cdag/N/F terms and pure A/Adag/F/N
    terms, simplify() must preserve the value."""
import numpy as np
from dmrgpy import fermionchain
from dmrgpy.multioperator import MultiOperator

L = 4
ODD = ("C","Cdag"); LOC = ("A","Adag")

def term(fs, c=1.0):
    m = MultiOperator(c=c)
    for (nm,i) in fs: m.add_operator(nm,i)
    return m

def predicted_sign(fs):
    k = 0
    for a in range(len(fs)):
        for b in range(a+1,len(fs)):
            if fs[a][0] in ODD and fs[b][0] in LOC and fs[b][1] < fs[a][1]: k += 1
    return (-1)**k

for v in (3,"python"):
    np.random.seed(11)
    fc = fermionchain.Fermionic_Chain(L,itensor_version=v)
    h0 = sum([fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] for i in range(L-1)])
    fc.set_hamiltonian(h0); fc.maxm = 16
    w1 = fc.random_mps(); w2 = fc.random_mps()
    A,C = fc.A,fc.C
    print("=== itensor_version =",repr(v))
    ref = abs(fc.aMb(w1,A[0]*C[1],w2))
    for label,T in [("A0 C1 + C1 A0 (anticommutator, j>i)",A[0]*C[1]+C[1]*A[0]),
                    ("A0 C1 - C1 A0 (commutator,     j>i)",A[0]*C[1]-C[1]*A[0]),
                    ("A1 C0 - C0 A1 (commutator,     j<i)",A[1]*C[0]-C[0]*A[1]),
                    ("A1 C0 + C0 A1 (anticommutator, j<i)",A[1]*C[0]+C[0]*A[1]),
                    ("A1 C3 + C3 A1 (anticommutator, j>i)",A[1]*C[3]+C[3]*A[1]),
                    ("A2 C0 - C0 A2 (commutator,     j<i)",A[2]*C[0]-C[0]*A[2])]:
        x = fc.aMb(w1,T,w2)
        print("  backend |<w1|%s|w2>| = %.3e   canonical: is_zero()=%s"
              %(label,abs(x),T.is_zero()))
    print("  scale: |<w1|A0 C1|w2>| = %.3e"%ref)
    rng = np.random.default_rng(5)
    def scan(names, ntr, check_pred):
        n_bad = n_pred_ok = n_nonzero = 0
        for t in range(ntr):
            nf = int(rng.integers(2,5))
            fs = [(names[rng.integers(len(names))], int(rng.integers(L))) for k in range(nf)]
            T = term(fs)
            x0 = fc.aMb(w1,T,w2); x1 = fc.aMb(w1,T.simplify(),w2)
            if abs(x0) > 1e-8: n_nonzero += 1
            if abs(x0-x1) > 1e-8: n_bad += 1
            if abs(x1-predicted_sign(fs)*x0) < 1e-8: n_pred_ok += 1
        return n_bad, n_pred_ok, n_nonzero
    for label,names in [("pure C/Cdag/N/F ",["C","Cdag","N","F"]),
                        ("pure A/Adag/N/F ",["A","Adag","N","F"]),
                        ("mixed C,A,N,F   ",["C","Cdag","A","Adag","N","F"])]:
        nb,npred,nnz = scan(names,80,True)
        print("  %s 80 terms (%d with nonzero value): simplify() changed %d; <simplify(T)> == (-1)^k <T> in %d/80"
              %(label,nnz,nb,npred))
```

```
=== itensor_version = 3
  backend |<w1|A0 C1 + C1 A0 (anticommutator, j>i)|w2>| = 0.000e+00   canonical: is_zero()=False
  backend |<w1|A0 C1 - C1 A0 (commutator,     j>i)|w2>| = 1.707e-02   canonical: is_zero()=True
  backend |<w1|A1 C0 - C0 A1 (commutator,     j<i)|w2>| = 0.000e+00   canonical: is_zero()=True
  backend |<w1|A1 C0 + C0 A1 (anticommutator, j<i)|w2>| = 1.707e-02   canonical: is_zero()=False
  backend |<w1|A1 C3 + C3 A1 (anticommutator, j>i)|w2>| = 2.554e-17   canonical: is_zero()=False
  backend |<w1|A2 C0 - C0 A2 (commutator,     j<i)|w2>| = 0.000e+00   canonical: is_zero()=True
  scale: |<w1|A0 C1|w2>| = 8.535e-03
  pure C/Cdag/N/F  80 terms (67 with nonzero value): simplify() changed 0; <simplify(T)> == (-1)^k <T> in 80/80
  pure A/Adag/N/F  80 terms (64 with nonzero value): simplify() changed 0; <simplify(T)> == (-1)^k <T> in 80/80
  mixed C,A,N,F    80 terms (66 with nonzero value): simplify() changed 4; <simplify(T)> == (-1)^k <T> in 80/80
=== itensor_version = 'python'
  backend |<w1|A0 C1 + C1 A0 (anticommutator, j>i)|w2>| = 0.000e+00   canonical: is_zero()=False
  backend |<w1|A0 C1 - C1 A0 (commutator,     j>i)|w2>| = 1.516e-01   canonical: is_zero()=True
  backend |<w1|A1 C0 - C0 A1 (commutator,     j<i)|w2>| = 0.000e+00   canonical: is_zero()=True
  backend |<w1|A1 C0 + C0 A1 (anticommutator, j<i)|w2>| = 1.516e-01   canonical: is_zero()=False
  backend |<w1|A1 C3 + C3 A1 (anticommutator, j>i)|w2>| = 2.772e-17   canonical: is_zero()=False
  backend |<w1|A2 C0 - C0 A2 (commutator,     j<i)|w2>| = 0.000e+00   canonical: is_zero()=True
  scale: |<w1|A0 C1|w2>| = 7.580e-02
  pure C/Cdag/N/F  80 terms (67 with nonzero value): simplify() changed 0; <simplify(T)> == (-1)^k <T> in 80/80
  pure A/Adag/N/F  80 terms (64 with nonzero value): simplify() changed 0; <simplify(T)> == (-1)^k <T> in 80/80
  mixed C,A,N,F    80 terms (66 with nonzero value): simplify() changed 4; <simplify(T)> == (-1)^k <T> in 80/80
```

```python
"""Reviewer probe 2.
(a) What each proposed fix does to the candidate operators, to a boson
    hopping, to a pure-C Hamiltonian and to a Jordan-Wigner-transformed one.
(b) mode="ED" on an A-type name.
(c) Is gs_energy()'s wrong number specific to the A/C defect, or what any
    false Hermiticity proof does? Force the proof on a pure-C non-Hermitian
    H and compare. Also repeat the candidate's g=1.5 run to see which
    backend moves from run to run."""
import numpy as np
from dmrgpy import fermionchain, bosonchain, multioperator
from dmrgpy.multioperatortk import canonical

A_TYPE = ("A","Adag","Aup","Adagup","Adn","Adagdn")
orig_sig = canonical._canonical_signature
orig_par = dict(canonical._PARITY)
def narrow(term):
    names = [o[0] for o in term[1:]]
    if any(n in canonical._ODD for n in names) and any(n in A_TYPE for n in names):
        return tuple((o[0],o[1]) for o in term[1:]), term[0]
    return orig_sig(term)
def use(which):
    canonical._PARITY.clear(); canonical._PARITY.update(orig_par)
    canonical._canonical_signature = orig_sig
    if which=="none":
        for n in A_TYPE: canonical._PARITY[n] = None
    if which=="narrow": canonical._canonical_signature = narrow

L = 4
fc = fermionchain.Fermionic_Chain(L,itensor_version="python")
bc = bosonchain.Bosonic_Chain(3,maxnb=[3,3,3],itensor_version="python")
A,C,Cd,Ad = fc.A,fc.C,fc.Cdag,fc.Adag
hop = sum([Cd[i]*C[i+1]+Cd[i+1]*C[i] for i in range(L-1)])
hub = hop + sum([0.7*fc.N[i]*fc.N[i+1] for i in range(L-1)]) + 0.3*(Cd[0]*Cd[2]+C[2]*C[0])
ops = [("P = Adag0 C1 + A0 Cdag1 (anti-Hermitian)   is_hermitian", lambda: (Ad[0]*C[1]+A[0]*Cd[1]).is_hermitian()),
       ("X = C1 A0 - A0 C1 (norm 4)                 is_zero     ", lambda: (C[1]*A[0]-A[0]*C[1]).is_zero()),
       ("Y = C0 A1 - A1 C0 (exactly zero)           is_zero     ", lambda: (C[0]*A[1]-A[1]*C[0]).is_zero()),
       ("pure-C hopping+NN+pairing                  is_hermitian", lambda: hub.is_hermitian()),
       ("jordan_wigner(same), A/Adag/F only         is_hermitian", lambda: multioperator.jordan_wigner(hub).is_hermitian()),
       ("boson hopping Adag0 A1 + Adag1 A0          is_hermitian", lambda: (bc.Adag[0]*bc.A[1]+bc.Adag[1]*bc.A[0]).is_hermitian()),
       ("boson Adag2 A0 - A0 Adag2                  is_zero     ", lambda: (bc.Adag[2]*bc.A[0]-bc.A[0]*bc.Adag[2]).is_zero())]
print("%-62s %-8s %-8s %-8s"%("","current","none","narrow"))
for label,f in ops:
    row = []
    for w in ("current","none","narrow"):
        use(w); row.append(str(f()))
    use("current")
    print("%-62s %-8s %-8s %-8s"%(label,*row))

# (b) ED on an A-type name
try:
    fc.vev(A[0]*Ad[0],mode="ED"); print("mode=ED accepted A/Adag")
except Exception as e:
    print("mode=ED on A0*Adag0 raises %s: %s"%(type(e).__name__,str(e)[:80]))

# (c) generic consequence: a pure-C non-Hermitian H, proof forced to True
Q = Cd[0]*C[1] - Cd[1]*C[0]   # anti-Hermitian, pure C
g = 1.5
from functools import reduce
a_ = np.array([[0,1],[0,0]],dtype=complex); F_ = np.diag([1.,-1.]).astype(complex); I2 = np.eye(2,dtype=complex)
def site(op,i): return reduce(np.kron,[op if k==i else I2 for k in range(L)])
Am = [site(a_,i) for i in range(L)]; Fm = [site(F_,i) for i in range(L)]
Cm = [reduce(np.matmul,[Fm[k] for k in range(j)]+[Am[j]]) for j in range(L)]; Cdm = [x.conj().T for x in Cm]
Hm = sum([Cdm[i]@Cm[i+1]+Cdm[i+1]@Cm[i] for i in range(L-1)]) + g*(Cdm[0]@Cm[1]-Cdm[1]@Cm[0])
ev = np.linalg.eigvals(Hm); ev = ev[np.lexsort((ev.imag,np.round(ev.real,8)))]
print("pure-C H = hop + 1.5 (Cdag0 C1 - Cdag1 C0): exact lowest-Re eigenvalues (numpy)", np.round(ev[:3],6))
orig_is_herm = canonical.is_hermitian
for v in (3,"python"):
    f1 = fermionchain.Fermionic_Chain(L,itensor_version=v); f1.set_hamiltonian(hop+g*Q)
    f1.maxm = 20; f1.nsweeps = 20
    e_ok = f1.gs_energy(); herm = f1.is_hermitian(f1.hamiltonian)
    canonical.is_hermitian = lambda MO: True
    try:
        es = []
        for rep in range(3):
            f2 = fermionchain.Fermionic_Chain(L,itensor_version=v); f2.set_hamiltonian(hop+g*Q)
            f2.maxm = 20; f2.nsweeps = 20
            es.append(np.round(f2.gs_energy(),6))
    finally:
        canonical.is_hermitian = orig_is_herm
    print("v=%-6s pure-C: proof says %s -> gs_energy()=%s | proof FORCED True, 3 runs: %s"%(v,herm,np.round(e_ok,6),es))

# the candidate's own H = hop + 1.5 P, unseeded, 3 fresh chains per backend
P = Ad[0]*C[1] + A[0]*Cd[1]
Pm = Am[0].conj().T@Cm[1] + Am[0]@Cdm[1]
Hc = sum([Cdm[i]@Cm[i+1]+Cdm[i+1]@Cm[i] for i in range(L-1)]) + g*Pm
evc = np.linalg.eigvals(Hc); evc = evc[np.lexsort((evc.imag,np.round(evc.real,8)))]
print("candidate H = hop + 1.5 P: exact lowest-Re eigenvalues (numpy)", np.round(evc[:4],6))
for v in (3,"python"):
    es = []
    for rep in range(3):
        f3 = fermionchain.Fermionic_Chain(L,itensor_version=v); f3.set_hamiltonian(hop+g*P)
        f3.maxm = 20; f3.nsweeps = 20
        es.append(np.round(f3.gs_energy(),6))
    print("v=%-6s candidate H, current code, 3 unseeded runs: %s"%(v,es))
```

```
                                                               current  none     narrow
P = Adag0 C1 + A0 Cdag1 (anti-Hermitian)   is_hermitian        True     False    False
X = C1 A0 - A0 C1 (norm 4)                 is_zero             True     False    False
Y = C0 A1 - A1 C0 (exactly zero)           is_zero             True     False    False
pure-C hopping+NN+pairing                  is_hermitian        True     True     True
jordan_wigner(same), A/Adag/F only         is_hermitian        True     False    True
boson hopping Adag0 A1 + Adag1 A0          is_hermitian        True     False    True
boson Adag2 A0 - A0 Adag2                  is_zero             True     False    True
mode=ED on A0*Adag0 raises AttributeError: 'NoneType' object has no attribute 'op'
pure-C H = hop + 1.5 (Cdag0 C1 - Cdag1 C0): exact lowest-Re eigenvalues (numpy) [-1.246695-0.896799j -1.246695-0.j       -1.246695+0.j      ]
v=3      pure-C: proof says False -> gs_energy()=(-1.246695+0.896799j) | proof FORCED True, 3 runs: [np.float64(-1.034971), np.float64(-0.889513), np.float64(-1.304455)]
v=python pure-C: proof says False -> gs_energy()=(-1.246695+0.896799j) | proof FORCED True, 3 runs: [np.float64(-2.407147), np.float64(-2.265325), np.float64(-2.535319)]
candidate H = hop + 1.5 P: exact lowest-Re eigenvalues (numpy) [-1.246695-0.896799j -1.246695-0.j       -1.246695+0.j
 -1.246695+0.896799j]
v=3      candidate H, current code, 3 unseeded runs: [np.float64(-0.991981), np.float64(-1.150282), np.float64(-0.511334)]
v=python candidate H, current code, 3 unseeded runs: [np.float64(-2.431932), np.float64(-2.472752), np.float64(-2.399882)]
```

The two candidate fixes as pytest plugins (`fix_none.py`, `fix_narrow.py`, in
"Shared helpers") against `tests/test_multioperator_canonical.py`:

```
tests/test_multioperator_canonical.py:117: AssertionError
FAILED tests/test_multioperator_canonical.py::test_heisenberg_hamiltonian_is_proven_hermitian
1 failed, 16 passed in 3.26s
```

```
17 passed in 3.27s
```

and the TD consequence rerun with the narrow fix applied
(`run_leadF_with_narrow_fix.py`, in "Shared helpers"):

```
||A - B^dag|| (matrices) = 4.000,  ||Adag0 Cdag1 + Cdag1 Adag0|| = 0.0e+00
exact <N> in GS = 2.000000
exact E0 = -4.0492958122  max|Im M_n| = 4.578e-02  peak |C| = 0.0447
itensor_version='python' gs_energy = -4.0492958122   is_dagger_pair(A,B) = False
  default (shortcut decides)     runs=2  max|y-exact| = 3.063e-05 (0% of peak)  max|Im y| = 3.59e-02
  shortcut disabled (two runs)   runs=2  max|y-exact| = 3.063e-05 (0% of peak)  max|Im y| = 3.59e-02
itensor_version=3 gs_energy = -4.0492958122   is_dagger_pair(A,B) = False
  default (shortcut decides)     runs=2  max|y-exact| = 3.063e-05 (0% of peak)  max|Im y| = 3.59e-02
  shortcut disabled (two runs)   runs=2  max|y-exact| = 3.063e-05 (0% of peak)  max|Im y| = 3.59e-02
```

**Suggested fix**: leave unreordered any term that contains both an `_ODD` name
and an `A`-type name, which makes the proofs refuse it (the hunter's second
option, measured above as `narrow`). It refuses P, X and the TD pair, keeps the
boson, pure-`C` and Jordan-Wigner-transformed proofs, passes all 17 tests of
`test_multioperator_canonical.py`, closes the TD consequence on its own
(3.616e-02 back to 3.063e-05, two evolutions) and still proves
`A = B.get_dagger()` for a mixed B, since identical spellings collect; its only
cost is that an exactly zero mixed operator such as `C0 A1 - A1 C0` is no longer
proven zero, the safe direction. The hunter's first option, parity `None` for
the six names, is wrong: `A`/`Adag` are also the ladder operators of
`Bosonic_Chain` and `SpinBoson_Chain`, where even is right, so a boson hopping
and every Jordan-Wigner-transformed operator lose their proof and
`test_heisenberg_hamiltonian_is_proven_hermitian` fails. Modelling the sign
instead would bake the finite-chain string convention (string from site 0) into
`canonical.py`, which the infinite-chain path does not use
(`to_terms(jordan_wigner_transform=False)`). With the narrow fix the "exact"
statements about the rewrite stay true; what needs correcting is the comment at
lines 72 to 75 and the docstring's list of what the module refuses, mirrored in
`CLAUDE.md` and `documentation.md` section 4.2a. Regression: `P.is_hermitian()`
False, `X.is_zero()` False, `<simplify(Cdag3 Adag0)>` equal to `<Cdag3 Adag0>`,
`is_dagger_pair` False on the TD pair, the boson hopping still proven. Numbers
change only for operators that mix the two representations, from wrong to
right.

### 2. `get_dagger()` passes the undocumented anti-Hermitian name `ISy` through unchanged, so KPM (the default), EX, TD and TDZ all return exactly minus the correlator whenever `ISy` sits in the first operator, and O1's "never a wrong number" guarantee is false

`bug` &middot; severity **LOW** (reach), correctness class &middot; CONFIRMED, NARROWED and widened &middot; lens `td-convention`

**Status**: FIXED -- `get_dagger()` carries a phase as well as a renaming, `multioperator._dagger_phase = {"ISy": -1.0}`, so ISy^dagger = -ISy; the phase multiplies the conjugated coefficient once per occurrence, daggering twice gives the operator back, and a term without `ISy` comes out byte-identical. Every other name the site types in `get_sites.h` build is Hermitian or paired in `_dagger_name` (`U+`/`U-` exist only in the never-instantiated `customspin.h`, and `XUp`/`XZ0`/`XDn` are state vectors). `ISy` stays off `canonical._PARITY`, so the proof refuses it and the chain's probe now sees the phase: `sc.is_hermitian(ISy)` is False and `sc.is_hermitian(1j*ISy)` True, where both were the other way round. The shared "adjoint for a correlator" helper the reviewer suggested as defence in depth is not added, its five call sites lying in other clusters' files, and stays open. Pinned by `tests/test_audit_2026_09_24_operators.py::test_isy_dagger_is_minus_isy`, `::test_chain_hermiticity_probe_sees_the_phase_of_isy` and `::test_isy_correlator_equals_its_1j_sy_spelling` (KPM and TD). NUMBERS CHANGE only where `ISy` sits in the first operator of a correlator or inside a hermitized operator: on the 4-site S=1/2 Heisenberg chain with hx=0.35, hz=0.15 (`"python"`, pair `(get_operator("ISy",0), Sz[1])`), KPM at delta=0.2 goes from exactly minus the `1j*Sy[0]` curve (a difference of 0.375 on a 0.187 peak) to identical, TD at delta=0.4, dt=0.1 from a difference of 0.086 on a 0.043 peak to 0.0, and the KPM first moment from +0.079470 to -0.079470 (exact -0.079685).

**Where**: `src/dmrgpy/multioperator.py` (`get_dagger` and its `_dagger_name`
map, which leave every unmapped name untouched). Consumers that turn the first
operator A into a ket through g(A) = A.get_dagger(): `timedependent.py:164` (the
TD forward run), `tdz.py:307` (the TDZ forward run),
`timedependent._adjoint_pair` (the second run `765b537` added), `kpmdmrg.py:127`
(KPM) and `dcex.py:61` (EX). The false guarantee:
`timedependent.lehmann_density_from_one_sided`'s docstring and the O1 Status
paragraph of `audit_2026_09_hole_hunt.md` (near line 2209), both saying that
refusing a pair "costs a second evolution and never a wrong number".

`is_dagger_pair` correctly refuses a pair containing a name with no known
adjoint, so TD runs twice, but both runs are built on `get_dagger()`, which
treats that name as Hermitian. The only non-Hermitian single-site name outside
the map on the reachable site types is `ISy`, which is i*Sy, anti-Hermitian, so
its true adjoint is -ISy. It is defined on the spin-1/2, spin-1 and spin-2 sites
of both vendored ITensors and mirrored in `pyitensor/sites/spin.py`, and it is
reachable only as a raw backend string through
`Many_Body_Chain.get_operator("ISy", i)`: it appears nowhere in `docs/`,
`README.md`, `examples/`, `tests/` or `operatornames.py`, and `mode="ED"` does
not know it (`pychain/build.py:66` dies on a bare `raise`), so there is no
in-library cross-check for it. The forward-run half predates `765b537` (`git log
-S` puts it at `7b71d8b`); `_adjoint_pair` and the "never a wrong number"
sentence first appear in `765b537`. Every other reachable name outside the
parity table is Hermitian (projUp/projDn, projEmp/projOcc, n, Nupdn,
FermiPhase), and Sig/Tau are inside `_dagger_name`, which is why nothing else
shows this.

**Expected**: `<GS|A^dagger(t) B|GS>`-type correlators of `get_operator("ISy",0)`
equal those of `1j*Sy[0]` pointwise.

Repro, from the hunter (it imports `anchor.py`, in "Shared helpers"):

```python
"""Lead 1: a pair naming an operator with no known adjoint ("ISy") is
refused by is_dagger_pair, so TD takes the two-run path -- but the second
run's adjoint pair is built with get_dagger(), which leaves "ISy"
untouched, i.e. treats the anti-Hermitian ISy as Hermitian."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from anchor import site_op, sx, sy, sz, heis_field, lehmann_dense, density
from dmrgpy import spinchain, timedependent
from dmrgpy.multioperatortk import canonical

n, hx, hz = 4, 0.35, 0.15
ES = np.linspace(-1.0, 6.0, 40)
DELTA, DT = 0.4, 0.1

# exact reference, by hand: A = i*Sy_0 (anti-Hermitian), B = Sz_1
H = heis_field(n, hx=hx, hz=hz)
A = 1j * site_op(sy, 0, n)
B = site_op(sz, 1, n)
D, M, e = lehmann_dense(H, A, B)
ref = density(D, M, ES, DELTA)
print("exact E0 = %.10f, max|Im M_n| = %.3e, max|Re M_n| = %.3e, peak |C| = %.4f"
      % (e[0], np.max(np.abs(M.imag)), np.max(np.abs(M.real)),
         np.max(np.abs(ref))))

calls = []
orig = timedependent.evolution_DC
def counted(*a, **k):
    calls.append(1)
    return orig(*a, **k)
timedependent.evolution_DC = counted

for version in ("python", 3):
    sc = spinchain.Spin_Chain([2] * n, itensor_version=version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n): h = h + hx * sc.Sx[i] + hz * sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 12
    print("\nitensor_version=%r  gs_energy = %.10f" % (version, sc.gs_energy()))
    ISy0 = sc.get_operator("ISy", 0)
    spellings = [("get_operator('ISy',0)", ISy0),
                 ("1j*Sy[0]", 1j * sc.Sy[0])]
    for label, Aop in spellings:
        pair = [Aop, sc.Sz[1]]
        del calls[:]
        _x, y = sc.get_dynamical_correlator(mode="DMRG", submode="TD",
                    name=pair, es=ES, delta=DELTA, dt=DT)
        y = np.asarray(y)
        err = np.max(np.abs(y - ref))
        print("  A=%-22s is_dagger_pair=%s  evolutions=%d  max|y-exact| = %.3e"
              "  (rel. to peak %.1f%%)  max|y+exact| = %.3e"
              % (label, canonical.is_dagger_pair(*pair), len(calls), err,
                 100 * err / np.max(np.abs(ref)), np.max(np.abs(y + ref))))
    # and the adjoint the two-run path actually used for the ISy spelling
    adj = timedependent._adjoint_pair([ISy0, sc.Sz[1]])
    # the forward run alone already daggers A with get_dagger()
    print("  ISy0.get_dagger().op =", ISy0.get_dagger().op)
    print("  _adjoint_pair([ISy_0,Sz_1])[1].op =", adj[1].op,
          " (the true adjoint of ISy is -ISy)")
```

Observed:

```
exact E0 = -1.6160254038, max|Im M_n| = 0.000e+00, max|Re M_n| = 6.025e-02, peak |C| = 0.0431

itensor_version='python'  gs_energy = -1.6160254038
  A=get_operator('ISy',0)  is_dagger_pair=False  evolutions=2  max|y-exact| = 8.609e-02  (rel. to peak 199.8%)  max|y+exact| = 8.482e-05
  A=1j*Sy[0]               is_dagger_pair=False  evolutions=2  max|y-exact| = 8.482e-05  (rel. to peak 0.2%)  max|y+exact| = 8.609e-02
  ISy0.get_dagger().op = [[np.float64(1.0), ['ISy', 0]]]
  _adjoint_pair([ISy_0,Sz_1])[1].op = [[np.float64(1.0), ['ISy', 0]]]  (the true adjoint of ISy is -ISy)

itensor_version=3  gs_energy = -1.6160254038
  A=get_operator('ISy',0)  is_dagger_pair=False  evolutions=2  max|y-exact| = 8.609e-02  (rel. to peak 199.8%)  max|y+exact| = 8.482e-05
  A=1j*Sy[0]               is_dagger_pair=False  evolutions=2  max|y-exact| = 8.482e-05  (rel. to peak 0.2%)  max|y+exact| = 8.609e-02
  ISy0.get_dagger().op = [[np.float64(1.0), ['ISy', 0]]]
  _adjoint_pair([ISy_0,Sz_1])[1].op = [[np.float64(1.0), ['ISy', 0]]]  (the true adjoint of ISy is -ISy)
```

**Reviewer (CONFIRMED, NARROWED and widened)**: reproduced to every printed
digit on `"python"` and v3, and also on v2, which the hunter did not run. The
claim narrows in one direction and widens in another. Narrower: the number is
wrong only when `ISy` sits in A. With `ISy` only in B every submode is right,
because g is an involution even where it is wrong, so the second TD run's
quench(g(g(B)), g(A)) undoes the wrong dagger `_adjoint_pair` took of B and both
runs depend on g only through g(A). Wider: the same exact sign flip happens
under the default `submode="KPM"`, under `"EX"` and under `"TDZ"`, so it is not a
TD defect; CVM on its conjugate-gradient path, which applies A directly at
`cvm.py:259`, is unaffected in both positions. On KPM the sign is anchored on
the first moment: exact sum_n M_n D_n = -0.079685, the `1j*Sy` spelling
-0.079470, the `ISy` spelling +0.079470. The hunter's suggested fix (raise at
the TD sites) has the wrong scope and outcome: it leaves KPM and EX returning
minus the correlator and turns a computable answer into a refusal. The
reviewer's probes (both import from their own folder):

```python
"""Reviewer probe 1: submode="TD" with an unknown anti-Hermitian name, in
BOTH positions of the pair, on v2, v3 and "python".

Prediction from the code (g = get_dagger, dagger = true adjoint):
  forward run  = quench(g(A), B)     = <GS| g(A)^dag e^{-iHt} B |GS>
  second run   = quench(g(g(B)), g(A)) = <GS| B^dag e^{-iHt} g(A) |GS>
so both runs are right iff g(A) = A^dag, and B's name never matters
(g is an involution, so the dagger _adjoint_pair takes of B is undone by
the second run's own forward dagger). Expect (ISy,Sz) = -exact and
(Sz,ISy) = +exact.

Exact reference: own numpy build, independent of dmrgpy and of the
hunter's anchor.py."""
import numpy as np
from dmrgpy import spinchain, timedependent
from dmrgpy.multioperatortk import canonical

n, hx, hz = 4, 0.35, 0.15
ES = np.linspace(-1.0, 6.0, 40)
DELTA, DT = 0.4, 0.1

s = {"x": np.array([[0, 1], [1, 0]]) / 2.0,
     "y": np.array([[0, -1j], [1j, 0]]) / 2.0,
     "z": np.array([[1, 0], [0, -1]]) / 2.0}


def op(m, i):
    out = np.eye(1)
    for k in range(n):
        out = np.kron(out, m if k == i else np.eye(2))
    return out


H = sum(op(s[a], i) @ op(s[a], i + 1) for i in range(n - 1) for a in "xyz")
H = H + sum(hx * op(s["x"], i) + hz * op(s["z"], i) for i in range(n))
e, U = np.linalg.eigh(H)
g0 = U[:, 0]


def exact(Am, Bm):
    # C_AB(w) = sum_n <0|A|n><n|B|0> L_delta(w - (E_n - E_0))
    M = (g0.conj() @ Am @ U) * (U.conj().T @ Bm @ g0)
    D = e - e[0]
    return np.array([np.sum(M * (DELTA / np.pi) / ((w - D) ** 2 + DELTA ** 2))
                     for w in ES])


ISy_m = 1j * op(s["y"], 0)       # ISy on site 0, the matrix ITensor defines
Sz_m = op(s["z"], 1)
print("check: ISy matrix anti-Hermitian: %s" % np.allclose(ISy_m.conj().T, -ISy_m))
refs = {"AB": exact(ISy_m, Sz_m), "BA": exact(Sz_m, ISy_m)}
for k, v in refs.items():
    print("exact peak |C| for order %s = %.4f" % (k, np.max(np.abs(v))))

calls = []
orig = timedependent.evolution_DC
def counted(*a, **k):
    calls.append(1)
    return orig(*a, **k)
timedependent.evolution_DC = counted

for version in ("python", 3, 2):
    sc = spinchain.Spin_Chain([2] * n, itensor_version=version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n): h = h + hx * sc.Sx[i] + hz * sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 12
    print("\nitensor_version=%r  gs_energy = %.10f (exact %.10f)"
          % (version, sc.gs_energy(), e[0]))
    ISy0 = sc.get_operator("ISy", 0)
    cases = [("(ISy_0, Sz_1)", [ISy0, sc.Sz[1]], "AB"),
             ("(1j*Sy_0, Sz_1)", [1j * sc.Sy[0], sc.Sz[1]], "AB"),
             ("(Sz_1, ISy_0)", [sc.Sz[1], ISy0], "BA"),
             ("(Sz_1, 1j*Sy_0)", [sc.Sz[1], 1j * sc.Sy[0]], "BA")]
    for label, pair, key in cases:
        ref = refs[key]
        del calls[:]
        _x, y = sc.get_dynamical_correlator(mode="DMRG", submode="TD",
                    name=list(pair), es=ES, delta=DELTA, dt=DT)
        y = np.asarray(y)
        pk = np.max(np.abs(ref))
        print("  %-17s is_dagger_pair=%-5s runs=%d  max|y-exact|=%.3e (%5.1f%% of peak)"
              "  max|y+exact|=%.3e"
              % (label, canonical.is_dagger_pair(*pair), len(calls),
                 np.max(np.abs(y - ref)), 100 * np.max(np.abs(y - ref)) / pk,
                 np.max(np.abs(y + ref))))
```

```
check: ISy matrix anti-Hermitian: True
exact peak |C| for order AB = 0.0431
exact peak |C| for order BA = 0.0431

itensor_version='python'  gs_energy = -1.6160254038 (exact -1.6160254038)
  (ISy_0, Sz_1)     is_dagger_pair=False runs=2  max|y-exact|=8.609e-02 (199.8% of peak)  max|y+exact|=8.482e-05
  (1j*Sy_0, Sz_1)   is_dagger_pair=False runs=2  max|y-exact|=8.482e-05 (  0.2% of peak)  max|y+exact|=8.609e-02
  (Sz_1, ISy_0)     is_dagger_pair=False runs=2  max|y-exact|=8.482e-05 (  0.2% of peak)  max|y+exact|=8.609e-02
  (Sz_1, 1j*Sy_0)   is_dagger_pair=False runs=2  max|y-exact|=8.482e-05 (  0.2% of peak)  max|y+exact|=8.609e-02

itensor_version=3  gs_energy = -1.6160254038 (exact -1.6160254038)
  (ISy_0, Sz_1)     is_dagger_pair=False runs=2  max|y-exact|=8.609e-02 (199.8% of peak)  max|y+exact|=8.482e-05
  (1j*Sy_0, Sz_1)   is_dagger_pair=False runs=2  max|y-exact|=8.482e-05 (  0.2% of peak)  max|y+exact|=8.609e-02
  (Sz_1, ISy_0)     is_dagger_pair=False runs=2  max|y-exact|=8.482e-05 (  0.2% of peak)  max|y+exact|=8.609e-02
  (Sz_1, 1j*Sy_0)   is_dagger_pair=False runs=2  max|y-exact|=8.482e-05 (  0.2% of peak)  max|y+exact|=8.609e-02

itensor_version=2  gs_energy = -1.6160254038 (exact -1.6160254038)
  (ISy_0, Sz_1)     is_dagger_pair=False runs=2  max|y-exact|=8.608e-02 (199.8% of peak)  max|y+exact|=1.672e-04
  (1j*Sy_0, Sz_1)   is_dagger_pair=False runs=2  max|y-exact|=1.672e-04 (  0.4% of peak)  max|y+exact|=8.608e-02
  (Sz_1, ISy_0)     is_dagger_pair=False runs=2  max|y-exact|=1.672e-04 (  0.4% of peak)  max|y+exact|=8.608e-02
  (Sz_1, 1j*Sy_0)   is_dagger_pair=False runs=2  max|y-exact|=1.672e-04 (  0.4% of peak)  max|y+exact|=8.608e-02
```

```python
"""Reviewer probe 3: which KPM spelling is the right one? Anchor on the
integrated spectral weight, which does not depend on the kernel shape:
int C_AB(w) dw over a window holding every pole = sum_n M_n = <GS|A B|GS>
(exact, by numpy). Also print where the ED path's bare raise comes from."""
import traceback
import numpy as np
from dmrgpy import spinchain

n, hx, hz = 4, 0.35, 0.15
s = {"x": np.array([[0, 1], [1, 0]]) / 2.0,
     "y": np.array([[0, -1j], [1j, 0]]) / 2.0,
     "z": np.array([[1, 0], [0, -1]]) / 2.0}
def op(m, i):
    out = np.eye(1)
    for k in range(n):
        out = np.kron(out, m if k == i else np.eye(2))
    return out
H = sum(op(s[a], i) @ op(s[a], i + 1) for i in range(n - 1) for a in "xyz")
H = H + sum(hx * op(s["x"], i) + hz * op(s["z"], i) for i in range(n))
e, U = np.linalg.eigh(H)
g0 = U[:, 0]
A = 1j * op(s["y"], 0); B = op(s["z"], 1)
M = (g0.conj() @ A @ U) * (U.conj().T @ B @ g0)
print("poles D_n in [%.3f, %.3f]" % (0.0, e[-1] - e[0]))
D = e - e[0]
print("exact sum_n M_n = <GS|ISy_0 Sz_1|GS> = %.6f%+.6fj" % (M.sum().real, M.sum().imag))
print("exact first moment sum_n M_n D_n = %.6f%+.6fj" % ((M*D).sum().real, (M*D).sum().imag))

ES = np.linspace(-3.0, 12.0, 3001)
sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
for i in range(n): h = h + hx * sc.Sx[i] + hz * sc.Sz[i]
sc.set_hamiltonian(h)
sc.maxm, sc.nsweeps = 30, 12
sc.gs_energy()
for label, Aop in (("ISy_0", sc.get_operator("ISy", 0)), ("1j*Sy_0", 1j * sc.Sy[0])):
    _x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                name=[Aop, sc.Sz[1]], es=ES, delta=0.2)
    w = np.trapezoid(np.asarray(y), ES)
    m1 = np.trapezoid(ES*np.asarray(y), ES)
    ref = np.array([np.sum(M*(0.2/np.pi)/((x-D)**2+0.2**2)) for x in ES])
    proj = np.trapezoid(np.asarray(y).real*ref.real, ES)/np.trapezoid(ref.real**2, ES)
    print("KPM A=%-8s integrated weight = %.6f%+.6fj  first moment = %.6f%+.6fj"
          "  projection on exact Lorentzian = %+.3f"
          % (label, w.real, w.imag, m1.real, m1.imag, proj))

print("\nED traceback for the name ISy:")
try:
    sc.vev(sc.get_operator("ISy", 0), mode="ED")
except Exception:
    traceback.print_exc(limit=-3)
```

```
Traceback (most recent call last):
  File "<repo>/src/dmrgpy/pychain/build.py", line 51, in get_operator
    return self.sector_restrict(multioperator.MO2matrix(name,self),name)
                                ~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^
  File "<repo>/src/dmrgpy/multioperator.py", line 371, in MO2matrix
    otmp = otmp@obj.get_operator(term[0],term[1]) # multiply
                ~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^
  File "<repo>/src/dmrgpy/pychain/build.py", line 66, in get_operator
    raise
RuntimeError: No active exception to reraise
poles D_n in [0.000, 3.128]
exact sum_n M_n = <GS|ISy_0 Sz_1|GS> = 0.000000+0.000000j
exact first moment sum_n M_n D_n = -0.079685+0.000000j
KPM A=ISy_0    integrated weight = 0.000001-0.000000j  first moment = 0.079470-0.000000j  projection on exact Lorentzian = -1.614
KPM A=1j*Sy_0  integrated weight = -0.000001+0.000000j  first moment = -0.079470+0.000000j  projection on exact Lorentzian = +1.614

ED traceback for the name ISy:
ISy
```

**Suggested fix**: give the name map a phase, so that `get_dagger` maps `ISy` to
`ISy` with -1 on the coefficient; that fixes every consumer at once, including
ones not measured here (CVM's analytic continuation, `nonhermitian/kpm.py`,
`mpsjulialive`'s KPM and TD, `correlationentropy`, the hermitizations at
`manybodychain.py:853` and `thermal.py:65`). Only after that, if wanted, add
`ISy` to `canonical._EVEN`; doing it first would let `is_hermitian` prove `ISy`
Hermitian. As defence in depth, one shared "adjoint for a correlator" helper used
at all five sites, raising unless every name is paired, carries a phase, or is in
an explicit self-adjoint set (not `all_names_known`, which would refuse Sig/Tau
and the Hermitian projectors needlessly). Correct the "never a wrong number"
sentence in the docstring and in the O1 Status paragraph. Regression: KPM and TD
on `(get_operator("ISy",0), Sz[1])` equal `(1j*Sy[0], Sz[1])` pointwise. Numbers
change only where `ISy` sits in an A slot or in a hermitized operator, from
wrong to right.

### 3. `get_distribution(mode="ED")` returns a distribution whose total weight is exactly 1/scale (0.1000 at the default `scale=10`, exact 1), and with `xs=` its imaginary part is a copy of its real part; a 2-site `itensor_version=3` chain gets this without asking for ED

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `kpm-calibration`

**Status**: FIXED -- `edtk/distribution.py::distribution_kpm` calls `kpm.dm_vivj_energy(..., x=xs)`, the `x=` branch, and multiplies by `scale/np.pi`; both interpolation lines are gone, so the `.real` copy cannot come back. `dm_vivj_energy`'s own normalization is unchanged (it feeds the O2-calibrated ED correlator), and its docstring now says that both branches return pi/scale times the density and names both callers' compensations; the old `\delta` docstring, which raised a `SyntaxWarning` at byte-compile, is replaced. Pinned by `tests/test_audit_2026_09_24_kpm.py::test_ed_distribution_integrates_to_one` (`scale` 2, 10 and `None` on Sz_0+Sz_1, which auto-scales to 2), `::test_ed_distribution_peak_weights_match_the_exact_eigendecomposition` (abs 5e-4), `::test_ed_distribution_on_requested_points_is_real_and_normalized` and `::test_two_site_v3_chain_gets_a_normalized_distribution_too`. NUMBERS CHANGE: ED `get_distribution` moves by exactly the factor `scale`; on the 4-site open Heisenberg chain, X=Sz_0, default `scale=10`, delta=0.05, the weight goes 0.1000 to 1.0000 and the peak 0.517316 to 5.173159; for X=Sz_0+Sz_1 at `scale=2` the per-peak weights are now 0.022357/0.955282/0.022357 against the exact 0.022329/0.955342/0.022329; with `xs=` max|Im y| goes from 4.96 to exactly 0. Nothing moves on any DMRG backend.

**Where**: `src/dmrgpy/edtk/distribution.py::distribution_kpm` (line 33, which
compensates only with `/np.pi`; lines 35 and 36, which both interpolate
`ys2.real`); the root in `src/dmrgpy/algebra/kpm.py::dm_vivj_energy`
(`/scale*np.pi` then `ys = ys/scale` at lines 483 and 485, one factor of `scale`
too many, which its sibling `dm_ij_energy` does not have). Reached through
`get_distribution(mode="ED")`, through `sc.mode="ED"`, through `mode.py`'s
`ns<3` fallback on `itensor_version=3`, and through an explicit v2/v3 with no
compiled extension (read in `mode.py`, not run).

`dm_vivj_energy` returns pi/scale times the physical density and has exactly two
callers: `edtk/dynamics.py:181` compensates with `*half/np.pi`,
`distribution_kpm` only with `/np.pi`, so the ED distribution integrates to
exactly 1/scale. The contract is a normalized density: user guide section 14
defines P(x) = <GS|delta(X-x)|GS>, and `distribution_kpm`'s own docstring says
the same. The shape and the relative weights are right. At `scale=None` the
auto-scale is 2*max|eig X|, so the weight is 1/(2*max|eig X|), right only when
max|eig X| = 1/2, which is exactly the Sz_0 of spin 1/2 every quick check uses.
It survived because no test reads the numbers
(`test_get_distribution_dispatches_through_get_mode` checks only that
`get_distribution_moments(mode="ED")` raises), and because the 2026-08 record's
finding #14 fixed the dispatch onto this function with a reviewer note that
`mode="ED"` "returns a correct spectrum", written without reading a value (the
only number printed there is a 1.36e-07 tail). The code predates the window (the
double `/scale` was there before `30dc282`).

**Expected**: weight 1 on any normalized state, and a real result for a
Hermitian X on its own ground state.

Repro, from the hunter:

```python
# ED get_distribution: is the integrated weight 1/scale?  Exact: 1.
import numpy as np
from dmrgpy import spinchain
sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=3)
h = 0
for i in range(3):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
for scale in (None, 1.0, 2.0, 5.0, 10.0):
    kw = {} if scale == 10.0 else {"scale": scale}
    x, y = sc.get_distribution(mode="ED", X=sc.Sz[0], delta=0.05, **kw)
    print("ED  scale=%-5s weight=%.4f  (1/scale=%s)" % (scale, np.trapezoid(np.real(y), x), None if scale is None else round(1/scale, 4)))
for scale in (None, 2.0, 10.0):
    kw = {} if scale is None else {"scale": scale}
    x, y = sc.get_distribution(mode="DMRG", X=sc.Sz[0], delta=0.05, **kw)
    print("DMRG scale=%-5s weight=%.4f" % (scale, np.trapezoid(np.real(y), x)))
# xs given: the imaginary part
xs = np.linspace(-1, 1, 5)
x, y = sc.get_distribution(mode="ED", X=sc.Sz[0], delta=0.2, xs=xs)
print("ED xs=", xs, "\n  y =", np.round(y, 4))
```

```
ED  scale=None  weight=1.0000  (1/scale=None)
ED  scale=1.0   weight=1.0000  (1/scale=1.0)
ED  scale=2.0   weight=0.5000  (1/scale=0.5)
ED  scale=5.0   weight=0.2000  (1/scale=0.2)
ED  scale=10.0  weight=0.1000  (1/scale=0.1)
DMRG scale=None  weight=0.9997
DMRG scale=2.0   weight=1.0000
DMRG scale=10.0  weight=1.0000
ED xs= [-1.  -0.5  0.   0.5  1. ]
  y = [0.0002+0.0002j 0.1298+0.1298j 0.0003+0.0003j 0.1298+0.1298j
 0.0002+0.0002j]
```

```python
# Finding B, no-user-action leg: a 2-site chain on itensor_version=3 is
# routed to ED by mode.py (v3 two-site DMRG aborts below 3 sites), so
# get_distribution() lands in edtk/distribution.py without anybody asking.
# Exact: the distribution of Sz_0 integrates to 1 on any normalized state.
import numpy as np
from dmrgpy import spinchain
for v in (2, "python", 3):
    sc = spinchain.Spin_Chain(["S=1/2"]*2, itensor_version=v)
    sc.set_hamiltonian(sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1])
    x, y = sc.get_distribution(X=sc.Sz[0], delta=0.05)
    xs = np.linspace(-1, 1, 3)
    _, y3 = sc.get_distribution(X=sc.Sz[0], delta=0.2, xs=np.array([-0.5, 0.0, 0.5]))
    print("itensor_version=%-7s get_mode()=%-4s weight=%.4f   y(xs=[-0.5,0,0.5])=%s"
          % (v, sc.get_mode(), np.trapezoid(np.real(y), x), np.round(np.asarray(y3), 4)))
```

```
itensor_version=2       get_mode()=DMRG weight=0.9997   y(xs=[-0.5,0,0.5])=[1.8696+0.j 0.0054+0.j 1.8696+0.j]
itensor_version=python  get_mode()=DMRG weight=0.9997   y(xs=[-0.5,0,0.5])=[1.8696+0.j 0.0054+0.j 1.8696-0.j]
itensor_version=3       get_mode()=ED   weight=0.1000   y(xs=[-0.5,0,0.5])=[0.1298+0.1298j 0.0003+0.0003j 0.1298+0.1298j]
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced digit for digit, against anchors
independent of every DMRG backend (the exact eigendecomposition of X in the ED
ground state, and the same file's own `method="INV"` resolvent route, which
integrates to 0.983 to 0.997, the shortfall being the Lorentzian tail leaving the
finite window). Per-peak weights match the exact eigendecomposition times
1/scale to three digits, 0.0112/0.4774/0.0112 at `scale=2` against exact
0.022329/0.955342/0.022329. No caller in `src/`, `examples/` or `tests/`
multiplies back by `scale`; the two examples that would compare against ED have
that line commented out. Struck: the hunter's "14x" peak ratio at x=+-0.5, which
bundles the 10x normalization with a roughly 1.4x width ratio from the two
routes' separate, documented-off-calibration moment counts
(`4*int(scale/delta)` on ED against `int(3*scale/delta)` on DMRG); after the fix
max|ED-DMRG| is still 2.34 on a DMRG peak of 7.58 at `scale=2`, a separate width
item that is not this finding. The reviewer's probes (its `r1_anchor.py` has one
probe bug of its own, rows labelled `scale=None` there are really the default
`scale=10`, which is why `r3` is the one quoted):

```python
# r1 had a probe bug: its "scale=None" row passed nothing, i.e. the default
# scale=10. Here scale=None is passed explicitly, so distribution_kpm's own
# auto-scale (2*max|eig X|) runs. Prediction if the weight is 1/scale:
# Sz_0 -> 1/1 = 1, Sz_0+Sz_1 -> 1/2, Sz_tot -> 1/4.
import numpy as np, dmrgpy
from dmrgpy import spinchain
print("dmrgpy from", dmrgpy.__file__)
sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=3)
h = 0
for i in range(3):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
ops = (("Sz_0", sc.Sz[0]), ("Sz_0+Sz_1", sc.Sz[0]+sc.Sz[1]),
       ("Sz_tot", sc.Sz[0]+sc.Sz[1]+sc.Sz[2]+sc.Sz[3]))
for name, X in ops:
    x, y = sc.get_distribution(mode="ED", X=X, delta=0.05, scale=None)
    xd, yd = sc.get_distribution(mode="DMRG", X=X, delta=0.05)  # DMRG default is scale=None too
    print("X=%-9s scale=None: ED weight=%.4f (grid half-width %.3f -> auto-scale %.3f)   DMRG weight=%.4f"
          % (name, np.trapezoid(np.real(y), x), x.max(), x.max()/0.95, np.trapezoid(np.real(yd), xd)))
```

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
X=Sz_0      scale=None: ED weight=1.0000 (grid half-width 0.950 -> auto-scale 1.000)   DMRG weight=0.9997
X=Sz_0+Sz_1 scale=None: ED weight=0.5000 (grid half-width 1.900 -> auto-scale 2.000)   DMRG weight=1.0000
X=Sz_tot    scale=None: ED weight=0.2500 (grid half-width 3.800 -> auto-scale 4.000)   DMRG weight=1.0000
```

```python
# Reviewer probe for candidate B: the no-user-action legs (2-site v3
# fallback, explicit sc.mode="ED"), and both candidate fixes applied as a
# local monkeypatch of distribution_kpm (repo untouched).
import numpy as np, dmrgpy
from scipy.interpolate import interp1d
from dmrgpy import spinchain
from dmrgpy.edtk import distribution as edd
from dmrgpy.algebra import kpm
print("dmrgpy from", dmrgpy.__file__)

def chain(n, v):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = 40, 30
    return sc

xs = np.array([-0.5, 0.0, 0.5])
print("--- 2-site chain, all defaults, only itensor_version differs")
for v in (2, "python", 3):
    sc = chain(2, v)
    x, y = sc.get_distribution(X=sc.Sz[0], delta=0.05)
    print("v=%-7s get_mode()=%-4s weight=%.4f" % (v, sc.get_mode(), np.trapezoid(np.real(y), x)))

print("--- 5-site chain, sc.mode='ED' set on the chain, get_distribution() called with no mode=")
sc = chain(5, 3); sc.mode = "ED"
x, y = sc.get_distribution(X=sc.Sz[2], delta=0.05)
print("weight=%.4f" % np.trapezoid(np.real(y), x))

orig = edd.distribution_kpm

def fix_hunter(wf0, X=None, scale=10.0, delta=1e-1, xs=None):
    # hunter's fix: ys2 *= scale/pi, second interpolation reads .imag
    M = X
    n = int(scale/delta)
    xs2, ys2 = kpm.dm_vivj_energy(M, wf0, wf0, scale=scale, npol=n*4, ne=n*10)
    ys2 = ys2*scale/np.pi
    if xs is None: return xs2, ys2
    ys = interp1d(xs2, ys2.real, fill_value=0., bounds_error=False)(xs)
    ys = ys + 1j*interp1d(xs2, ys2.imag, fill_value=0., bounds_error=False)(xs)
    return xs, ys

def fix_direct(wf0, X=None, scale=10.0, delta=1e-1, xs=None):
    # alternative: let dm_vivj_energy evaluate at xs itself (its x= branch,
    # added in 30dc282 for exactly this), no interpolation at all
    M = X
    n = int(scale/delta)
    xs2, ys2 = kpm.dm_vivj_energy(M, wf0, wf0, scale=scale, npol=n*4, ne=n*10, x=xs)
    return xs2, ys2*scale/np.pi

sc4 = chain(4, 3)
X = sc4.Sz[0] + sc4.Sz[1]
xg = np.linspace(-1.5, 1.5, 3001)
_, yd = sc4.get_distribution(mode="DMRG", X=X, delta=0.05, scale=2.0, xs=xg)
for label, f in (("as shipped", orig), ("hunter fix", fix_hunter), ("direct-x fix", fix_direct)):
    edd.distribution_kpm = f
    rows = []
    for scale in (2.0, 10.0):
        x, y = sc4.get_distribution(mode="ED", X=X, delta=0.05, scale=scale)
        _, yx = sc4.get_distribution(mode="ED", X=X, delta=0.05, scale=scale, xs=xg)
        rows.append("scale=%-4s weight=%.4f  xs-weight=%.4f  max|Im|(xs)=%.1e"
                    % (scale, np.trapezoid(np.real(y), x), np.trapezoid(np.real(yx), xg), np.max(np.abs(np.imag(yx)))))
    # pointwise vs DMRG at the SAME scale=2 (same rescaling on both routes)
    _, ye = sc4.get_distribution(mode="ED", X=X, delta=0.05, scale=2.0, xs=xg)
    rows.append("vs DMRG(scale=2): max|ED-DMRG|=%.3e  on DMRG peak %.3f" % (np.max(np.abs(ye.real-yd.real)), yd.real.max()))
    print("---", label); [print("   ", r) for r in rows]
edd.distribution_kpm = orig
```

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
--- 2-site chain, all defaults, only itensor_version differs
v=2       get_mode()=DMRG weight=0.9997
v=python  get_mode()=DMRG weight=0.9997
v=3       get_mode()=ED   weight=0.1000
--- 5-site chain, sc.mode='ED' set on the chain, get_distribution() called with no mode=
weight=0.1000
--- as shipped
    scale=2.0  weight=0.5000  xs-weight=0.5000  max|Im|(xs)=5.0e+00
    scale=10.0 weight=0.1000  xs-weight=0.1000  max|Im|(xs)=9.9e-01
    vs DMRG(scale=2): max|ED-DMRG|=2.696e+00  on DMRG peak 7.579
--- hunter fix
    scale=2.0  weight=1.0000  xs-weight=1.0000  max|Im|(xs)=0.0e+00
    scale=10.0 weight=1.0000  xs-weight=1.0000  max|Im|(xs)=0.0e+00
    vs DMRG(scale=2): max|ED-DMRG|=2.342e+00  on DMRG peak 7.579
--- direct-x fix
    scale=2.0  weight=1.0000  xs-weight=1.0000  max|Im|(xs)=0.0e+00
    scale=10.0 weight=1.0000  xs-weight=1.0000  max|Im|(xs)=0.0e+00
    vs DMRG(scale=2): max|ED-DMRG|=2.343e+00  on DMRG peak 7.579
```

**Suggested fix**: call `kpm.dm_vivj_energy(..., x=xs)`, the `x=` branch that
`30dc282` added for exactly this caller, multiply by `scale/np.pi`, and delete
both interpolation lines, so the `.real` copy cannot come back. The hunter's
local fix (`ys2 *= scale/np.pi`, `.imag` on line 36) measures equally right, up
to about 1e-3 of interpolation. Normalizing `dm_vivj_energy` itself to match
`dm_ij_energy` is cleaner but touches the O2-calibrated ED correlator, so the
local fix is the lower-risk one, with a comment in `dm_vivj_energy` naming its
pi/scale normalization. The regression should pin the sum rule rather than
ED-versus-DMRG agreement: weight close to 1 at `scale=2`, `scale=10` and
`scale=None` on an operator whose auto-scale is not 1, per-peak weights against
the exact eigendecomposition, and max|Im| close to 0 with `xs=`. NUMBERS CHANGE:
ED `get_distribution` moves by exactly the factor `scale` (10x at the default);
nothing moves on any DMRG backend.

### 4. Every DMRG KPM route reconstructs the spectrum from n+2 Chebyshev moments where the O2 calibration and ED use n, so a DMRG line is narrower and taller than ED's by about 2/n, which on small chains is almost the whole of the ED-versus-DMRG residual O2 reports

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `kpm-calibration`

**Status**: FIXED -- `mus = np.array(moments)[:n]` in `kpmdmrg.dynamical_correlator_moments`, after both session calls (so v2, v3, v3's `kpm_dynamical_correlator_truncated`, `"python"` and through them `infinitechain.kpm_finite`), and in `mpsjulialive/dynamics.py::_kpm_dynamical_correlator`, both before `kpm_extrapolate`; the C++, pyitensor and `kpm.jl` loops are untouched, since `general_kpm` shares them. Both docstrings are re-measured: `dynamical_correlator_moments` now says n is the calibrated count and `len(mus) == n` unless `kpm_extrapolate` resamples, and `polynomials_for_broadening`'s "4 to 6 per cent, the kernel's own asymptotics" is replaced by the numpy-anchor numbers at N=n, FWHM/(2*delta) = 0.964 at n=16, 0.978 at 26, 0.992 at 52, 1.001 at 104, 1.005 at 200 and 1.009 at 1000, with the requirement that the kernel see exactly npol moments. Pinned by `tests/test_audit_2026_09_24_kpm.py::test_dmrg_kpm_returns_exactly_the_calibrated_moment_count` (n=52 and n=17, accelerate on and off, `"python"` and v3), `::test_energy_truncated_kpm_returns_exactly_the_calibrated_moment_count`, `::test_band_centre_pole_is_the_same_curve_on_dmrg_and_ed` (at the reconstruction's own grid points, measured 1e-14), `::test_kpm_accelerate_does_not_change_the_curve_at_odd_n` and `::test_julia_live_kpm_uses_exactly_the_calibrated_moment_count`. Two residuals remain and are recorded as leads rather than defects: `dynamical_correlator_from_moments` rebuilds the curve on its own 10n-point grid and interpolates linearly onto `es`, which costs 5.8e-4 of the peak on the dimer pole at delta=0.1, and v3's band-edge estimate `emax` lands up to 2e-2 below the true E_max from run to run (25 runs on the staggered 4-site chain: -2.1e-5 to -2.0e-2), which is now the whole of the ED-versus-v3 residual. NUMBERS CHANGE: every DMRG KPM spectrum moves by roughly 2/n onto ED's; on two decoupled Heisenberg dimers, <Sz_0;Sz_0> with its pole at the band centre, the public DMRG peak goes 0.666761 to 0.620591 at delta=0.2 and 1.266336 to 1.220236 at delta=0.1, on v3 and `"python"` alike, both now equal to ED's; the O2 record's residual max|ED-DMRG| on `"python"` goes 1.32e-02 to 2.55e-04 (4-site uniform, delta=0.15), 1.26e-02 to 2.09e-04 (6-site uniform, 0.20) and 1.09e-02 to 3.50e-04 (6-site plus field 0.3, 0.25), 0.07 to 0.20 per cent of the resolvent peak (the O2 table's 4-site row is the uniform chain, not the staggered one).

**Where**: the moment loops, `pyitensor/chain.py:1923`/`:1943`,
`mpscpp2/chain_session.h:1324`/`:1354`, `mpscpp3/chain_session.h:11845`/`:11872`
and the truncated variant at `:12198`/`:12260`, `kpm.jl`'s two loops, which
return 2+n moments (`_kpm_moments_full`) or 2+2*(n//2) (the accelerated path);
`kpmdmrg.dynamical_correlator_moments` (`kpmdmrg.py:155`,
`mus = np.array(moments)`) and `mpsjulialive/dynamics.py:173`, which pass them
all on; `algebra/kpm.py::jackson_kernel` (line 813), which sets N = len(mus).
Reached by v2, v3 through both KPM entry points, `"python"`, `julia_live` and
`infinitechain.kpm_finite`. The two docstrings it contradicts:
`kpmdmrg.dynamical_correlator_moments` ("n is the number of polynomials the
backend chose") and `polynomials_for_broadening` ("reproduces the observed FWHM
to 4 to 6 per cent at npol of order 10 to 100, the residual being the kernel's
own asymptotics").

`polynomials_for_broadening` solves for the n whose Jackson reconstruction has
FWHM = 2*delta at the band centre, and the ED route
(`edtk/dynamics.py::dynamical_correlator_kpm`) feeds the kernel exactly n. The
DMRG loops, which mirror the old file-based `kpmcorrelator.h` and predate the
calibration, return two more, so `765b537` routed the calibrated n into a loop
that already returned n+2. The moments themselves are exact, which is why the sum
rule cannot see it: this is a count, not a weight. It survived because
`test_ed_and_dmrg_kpm_agree_pointwise` tolerates 10 per cent of the peak and
`test_kpm_moment_count_is_calibrated_to_the_requested_broadening` tests the
function, never `len(mus)`.

**Expected**: `len(mus) == n` on every route, and the DMRG curve equal to ED's
where the moments are exact.

Repro, from the hunter: a single isolated pole at the band centre (two decoupled
Heisenberg dimers, one pole at omega=1, x=0):

```python
# Single isolated pole at the band centre: two decoupled Heisenberg dimers,
# J=J'=1.  <Sz_0;Sz_0> has exactly one pole, omega=1, weight 1/4, and
# the many-body band is [E0,Emax]=[-1.5,0.5], so omega=1 is x=0 exactly.
# Contract under test (algebra/kpm.py::polynomials_for_broadening):
# FWHM = 2*delta at the band centre.  Measure it on each route.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kpmdmrg import dynamical_correlator_from_moments

def chain(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=v)
    h = 0
    for (i,j) in [(0,1),(2,3)]:
        h = h + sc.Sx[i]*sc.Sx[j] + sc.Sy[i]*sc.Sy[j] + sc.Sz[i]*sc.Sz[j]
    sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = 20, 20
    return sc

def fwhm(x, y):
    y = np.real(y); k = np.argmax(y); h = y[k]/2.
    l = k
    while y[l] > h: l -= 1
    r = k
    while y[r] > h: r += 1
    xl = x[l] + (h-y[l])*(x[l+1]-x[l])/(y[l+1]-y[l])
    xr = x[r-1] + (h-y[r-1])*(x[r]-x[r-1])/(y[r]-y[r-1])
    return xr-xl, x[k], y[k]

es = np.linspace(-1.0, 3.0, 8001)
print("%-6s %-8s %5s %6s %9s %9s %9s  %s" % ("delta","route","n","nmus","FWHM","FWHM/2d","peak","weight"))
for delta in (0.05, 0.1, 0.2, 0.3, 0.45):
    sc = chain(3)
    x, y = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=[sc.Sz[0], sc.Sz[0]], delta=delta, es=es)
    f, xp, yp = fwhm(x, y)
    w = np.trapezoid(np.real(y), x)
    print("%-6.2f %-8s %5s %6s %9.5f %9.4f %9.5f  %.6f" % (delta, "ED", "", "", f, f/(2*delta), yp, w))
    for v in (2, 3, "python"):
        sc = chain(v)
        mus, emin, emax, scale, n, d = sc.get_dynamical_correlator_moments(name=[sc.Sz[0], sc.Sz[0]], delta=delta)
        x, y = dynamical_correlator_from_moments(mus, emin, emax, scale, n, es)
        f, xp, yp = fwhm(x, y); w = np.trapezoid(np.real(y), x)
        print("%-6.2f %-8s %5d %6d %9.5f %9.4f %9.5f  %.6f" % (delta, str(v), n, len(mus), f, f/(2*delta), yp, w))
        if v == 3:
            x, y = dynamical_correlator_from_moments(np.array(mus)[:n], emin, emax, scale, n, es)
            f, xp, yp = fwhm(x, y); w = np.trapezoid(np.real(y), x)
            print("%-6.2f %-8s %5d %6d %9.5f %9.4f %9.5f  %.6f" % (delta, "3,mus[:n]", n, n, f, f/(2*delta), yp, w))
    # the public call on v3, for completeness
    sc = chain(3)
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=[sc.Sz[0], sc.Sz[0]], delta=delta, es=es)
    f, xp, yp = fwhm(x, y)
    print("%-6.2f %-8s %5s %6s %9.5f %9.4f %9.5f" % (delta, "3,public", "", "", f, f/(2*delta), yp))
```

```
delta  route        n   nmus      FWHM   FWHM/2d      peak  weight
0.05   ED                      0.09964    0.9964   2.41850  0.250000
0.05   2          104    106   0.09779    0.9779   2.46458  0.250000
0.05   3          104    106   0.09779    0.9779   2.46458  0.250000
0.05   3,mus[:n]   104    104   0.09965    0.9965   2.41850  0.250000
0.05   python     104    106   0.09779    0.9779   2.46458  0.250000
0.05   3,public                0.09779    0.9779   2.46458
0.10   ED                      0.19762    0.9881   1.22024  0.249996
0.10   2           52     54   0.19042    0.9521   1.26634  0.249999
0.10   3           52     54   0.19042    0.9521   1.26634  0.249999
0.10   3,mus[:n]    52     52   0.19764    0.9882   1.22024  0.249998
0.10   python      52     54   0.19042    0.9521   1.26634  0.249999
0.10   3,public                0.19042    0.9521   1.26634
0.20   ED                      0.38954    0.9738   0.62059  0.249973
0.20   2           26     28   0.36240    0.9060   0.66676  0.249989
0.20   3           26     28   0.36240    0.9060   0.66676  0.249989
0.20   3,mus[:n]    26     26   0.38956    0.9739   0.62059  0.249986
0.20   python      26     28   0.36240    0.9060   0.66676  0.249989
0.20   3,public                0.36240    0.9060   0.66676
0.30   ED                      0.58844    0.9807   0.41258  0.249894
0.30   2           17     18   0.55680    0.9280   0.43571  0.249961
0.30   3           17     18   0.55680    0.9280   0.43571  0.249961
0.30   3,mus[:n]    17     17   0.58844    0.9807   0.41258  0.249950
0.30   python      17     18   0.55680    0.9280   0.43571  0.249961
0.30   3,public                0.55680    0.9280   0.43571
0.45   ED                      0.62410    0.6934   0.38940  0.249896
0.45   2           16     18   0.55680    0.6187   0.43571  0.249960
0.45   3           16     18   0.55680    0.6187   0.43571  0.249960
0.45   3,mus[:n]    16     16   0.62410    0.6934   0.38940  0.249936
0.45   python      16     18   0.55680    0.6187   0.43571  0.249960
0.45   3,public                0.55680    0.6187   0.43571
```

and the O2 regression's own chain, with the kpm_accelerate side effect at odd n:

```python
# The O2 record's residual ED-vs-DMRG disagreement (1.32e-02 on the
# 4-site staggered Heisenberg chain at delta=0.15, A=B=Sz_0), and whether
# it is the moment count: DMRG reconstructs from len(mus)=n+2 (or n+1),
# ED from exactly n.  Also: does kpm_accelerate (a speed flag) change
# the curve when n is odd?
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kpmdmrg import dynamical_correlator_from_moments

def staggered_heisenberg(n=4, v=3):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): h = h + 0.3*(-1)**i*sc.Sz[i]
    sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = 30, 20
    return sc

es = np.linspace(0.01, 5.0, 200)
for (nsite, delta) in [(4, 0.15), (6, 0.2), (4, 0.4)]:
    sc = staggered_heisenberg(nsite)
    _, yed = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=[sc.Sz[0], sc.Sz[0]], es=es, delta=delta)
    yed = np.real(yed)
    for v in (3, "python"):
        sc = staggered_heisenberg(nsite, v)
        _, ydm = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=[sc.Sz[0], sc.Sz[0]], es=es, delta=delta)
        sc2 = staggered_heisenberg(nsite, v)
        mus, emin, emax, scale, n, d = sc2.get_dynamical_correlator_moments(name=[sc2.Sz[0], sc2.Sz[0]], delta=delta)
        _, ytr = dynamical_correlator_from_moments(np.array(mus)[:n], emin, emax, scale, n, es)
        print("n=%d delta=%.2f v=%-6s npol=%d len(mus)=%d  peak ED=%.5f  max|ED-DMRG|=%.3e  max|ED-DMRG(mus[:n])|=%.3e"
              % (nsite, delta, v, n, len(mus), yed.max(), np.max(np.abs(yed-np.real(ydm))), np.max(np.abs(yed-np.real(ytr)))))

print("-- kpm_accelerate on/off, odd n")
for delta in (0.15, 0.3):
    for acc in (True, False):
        sc = staggered_heisenberg(4, 3); sc.kpm_accelerate = acc
        mus, emin, emax, scale, n, d = sc.get_dynamical_correlator_moments(name=[sc.Sz[0], sc.Sz[0]], delta=delta)
        _, y = dynamical_correlator_from_moments(mus, emin, emax, scale, n, es)
        print("delta=%.2f kpm_accelerate=%-5s n=%d len(mus)=%d peak=%.6f" % (delta, acc, n, len(mus), np.real(y).max()))
```

```
n=4 delta=0.15 v=3      npol=47 len(mus)=48  peak ED=0.56193  max|ED-DMRG|=1.114e-02  max|ED-DMRG(mus[:n])|=1.774e-04
n=4 delta=0.15 v=python npol=47 len(mus)=48  peak ED=0.56193  max|ED-DMRG|=1.111e-02  max|ED-DMRG(mus[:n])|=4.424e-04
n=6 delta=0.20 v=3      npol=56 len(mus)=58  peak ED=0.47511  max|ED-DMRG|=1.643e-02  max|ED-DMRG(mus[:n])|=2.902e-04
n=6 delta=0.20 v=python npol=56 len(mus)=58  peak ED=0.47511  max|ED-DMRG|=1.641e-02  max|ED-DMRG(mus[:n])|=2.931e-04
n=4 delta=0.40 v=3      npol=18 len(mus)=20  peak ED=0.22347  max|ED-DMRG|=2.318e-02  max|ED-DMRG(mus[:n])|=1.975e-04
n=4 delta=0.40 v=python npol=18 len(mus)=20  peak ED=0.22347  max|ED-DMRG|=2.317e-02  max|ED-DMRG(mus[:n])|=1.996e-04
-- kpm_accelerate on/off, odd n
delta=0.15 kpm_accelerate=True  n=47 len(mus)=48 peak=0.573068
delta=0.15 kpm_accelerate=False n=47 len(mus)=49 peak=0.584880
delta=0.30 kpm_accelerate=True  n=23 len(mus)=24 peak=0.293217
delta=0.30 kpm_accelerate=False n=23 len(mus)=25 peak=0.304933
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced on v3 and `"python"`; a Jackson
reconstruction in plain numpy, importing nothing from dmrgpy, shows that ED is
the route that honours the calibration and DMRG the one two moments off
(0.9922 at N=n and 0.9560 at N=n+2 at n=52, which times the rounding factor
0.99585 gives the measured 0.9881 and 0.9520). Nothing declares n+2 deliberate,
and `JACKSON_FWHM_FACTOR = 2*sqrt(2 ln 2)*pi` is the analytic constant, not a
fitted one, so it did not absorb the extra two. Three narrowings. "Changes every
DMRG KPM peak by 2 to 9 per cent" is struck as a general statement: that is the
isolated-line figure at n=16 to 104, and on a 20-site Heisenberg chain at n=116
to 348 the change is 0.5 to 0.7 per cent; the effect is roughly 2/n, a
statement about n rather than about L. "2 to 11 per cent narrower than the
contract" holds against ED and the calibration relation for n up to about 100
to 150; at the n=16 floor both routes sit below 2*delta by design, so the 11 per
cent there is DMRG against ED. "The docstring's 4 to 6 per cent is actually
this" is narrowed to what the anchor proves, that 4 to 6 per cent cannot be the
kernel's own asymptotics at any n the floor allows (3.6 per cent at n=16, 0.8 at
52, 0.0 at 104); the O2 record's pre-fix DMRG constant of 6.9 to 7.2 does not
decide between the two routes and should not be cited for it. The hunter's
`kpm_accelerate` digits do not reproduce (0.573170/0.586233 against
0.573068/0.584880) because the band-edge estimate `emax` varies from run to run,
a separate lead recorded under "New leads"; the effect stands. The reviewer's
anchor and its size-in-practice run:

```python
# Independent anchor, numpy only (no dmrgpy import): the band-centre FWHM of
# a Jackson-damped Chebyshev reconstruction of a delta peak at x0=0, as a
# function of the number of moments N (the kernel's own N = len(mus), the
# WWAF convention g_k, k=0..N-1).  Contract under test:
#   FWHM = 7.39786 * half_width / npol   with npol = the calibrated count.
# Report FWHM(N)*n/7.39786 for N=n and N=n+2 (and N=n+1), i.e. the fraction
# of the calibrated width a route that feeds the kernel N moments delivers.
import numpy as np
C = 2.0*np.sqrt(2.0*np.log(2.0))*np.pi

def jackson_g(N):
    k = np.arange(N); q = np.pi/(N+1.)
    return ((N-k+1)*np.cos(q*k) + np.sin(q*k)/np.tan(q))/(N+1)

def line(N, x, x0=0.0):
    k = np.arange(N)
    mu = np.cos(k*np.arccos(x0))          # moments of delta(x-x0)
    c = jackson_g(N)*mu; c[1:] *= 2
    T = np.cos(np.outer(np.arccos(x), k))
    return (T@c)/(np.pi*np.sqrt(1-x**2))

def fwhm(N):
    x = np.linspace(-0.6, 0.6, 400001)
    y = line(N, x); h = y.max()/2
    above = np.where(y > h)[0]; l, r = above[0], above[-1]
    xl = x[l-1] + (h-y[l-1])*(x[l]-x[l-1])/(y[l]-y[l-1])
    xr = x[r] + (h-y[r])*(x[r+1]-x[r])/(y[r+1]-y[r])
    return xr-xl, y.max()

print("%5s %12s %12s %12s   %s" % ("n", "N=n", "N=n+1", "N=n+2", "peak(n+2)/peak(n)"))
for n in (16, 17, 18, 23, 26, 47, 52, 56, 104, 134, 259, 500):
    f0, p0 = fwhm(n); f1, p1 = fwhm(n+1); f2, p2 = fwhm(n+2)
    print("%5d %12.4f %12.4f %12.4f   %.4f" % (n, f0*n/C, f1*n/C, f2*n/C, p2/p0))
# the O2 record's pre-fix DMRG measurement: constant 6.9 to 7.2 at npol 9 to 51
print("O2 record, pre-fix DMRG constant if the kernel saw N=n+2 moments:")
for n in (9, 15, 20, 30, 40, 51):
    f2, _ = fwhm(n+2); f0, _ = fwhm(n)
    print("  npol=%2d  FWHM(n+2)*n = %.3f   FWHM(n)*n = %.3f   (analytic %.4f)" % (n, f2*n, f0*n, C))
```

```
    n          N=n        N=n+1        N=n+2   peak(n+2)/peak(n)
   16       0.9641       0.9090       0.8601   1.1189
   17       0.9659       0.9139       0.8672   1.1122
   18       0.9677       0.9182       0.8737   1.1062
   23       0.9746       0.9351       0.8987   1.0838
   26       0.9779       0.9426       0.9098   1.0744
   47       0.9905       0.9702       0.9508   1.0417
   52       0.9922       0.9738       0.9560   1.0378
   56       0.9933       0.9762       0.9596   1.0351
  104       1.0006       0.9911       0.9818   1.0191
  134       1.0026       0.9952       0.9879   1.0148
  259       1.0060       1.0021       0.9983   1.0077
  500       1.0078       1.0058       1.0038   1.0040
O2 record, pre-fix DMRG constant if the kernel saw N=n+2 moments:
  npol= 9  FWHM(n+2)*n = 5.767   FWHM(n)*n = 7.013   (analytic 7.3979)
  npol=15  FWHM(n+2)*n = 6.305   FWHM(n)*n = 7.117   (analytic 7.3979)
  npol=20  FWHM(n+2)*n = 6.547   FWHM(n)*n = 7.181   (analytic 7.3979)
  npol=30  FWHM(n+2)*n = 6.817   FWHM(n)*n = 7.260   (analytic 7.3979)
  npol=40  FWHM(n+2)*n = 6.965   FWHM(n)*n = 7.306   (analytic 7.3979)
  npol=51  FWHM(n+2)*n = 7.065   FWHM(n)*n = 7.338   (analytic 7.3979)
```

```python
# Does the size survive at the chain lengths DMRG is used for?  The extra
# two moments shift the width by about 2/n, and n = C*half_width/(2*delta)
# grows with the many-body bandwidth, i.e. with L.  Uniform S=1/2
# Heisenberg chain, v3, <Sz_{L/2};Sz_{L/2}>, same moments reconstructed
# full and cut to mus[:n].
import numpy as np, time
import dmrgpy
from dmrgpy import spinchain
from dmrgpy.kpmdmrg import dynamical_correlator_from_moments

def heis(L):
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=3)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = 40, 12
    return sc

es = np.linspace(0.0, 4.0, 2001)
for (L, delta) in ((8, 0.1), (20, 0.1), (20, 0.3)):
    t0 = time.time()
    sc = heis(L); c = L//2
    mus, emin, emax, scale, n, d = sc.get_dynamical_correlator_moments(name=[sc.Sz[c], sc.Sz[c]], delta=delta)
    mus = np.array(mus)
    _, yf = dynamical_correlator_from_moments(mus, emin, emax, scale, n, es)
    _, yc = dynamical_correlator_from_moments(mus[:n], emin, emax, scale, n, es)
    yf, yc = np.real(yf), np.real(yc)
    print("L=%2d delta=%.2f W=%.3f half=%.3f n=%d len(mus)=%d  peak full=%.5f [:n]=%.5f ratio=%.4f  max|full-cut|/peak=%.4f  2/n=%.4f  (%.0fs)"
          % (L, delta, emax-emin, 1/scale, n, len(mus), yf.max(), yc.max(), yf.max()/yc.max(), np.max(np.abs(yf-yc))/yc.max(), 2./n, time.time()-t0))
```

```
L= 8 delta=0.10 W=5.125 half=3.587 n=133 len(mus)=134  peak full=0.45096 [:n]=0.44763 ratio=1.0074  max|full-cut|/peak=0.0074  2/n=0.0150  (0s)
L=20 delta=0.10 W=13.431 half=9.402 n=348 len(mus)=350  peak full=0.25603 [:n]=0.25457 ratio=1.0057  max|full-cut|/peak=0.0057  2/n=0.0057  (14s)
L=20 delta=0.30 W=13.432 half=9.402 n=116 len(mus)=118  peak full=0.15153 [:n]=0.15077 ratio=1.0050  max|full-cut|/peak=0.0071  2/n=0.0172  (5s)
```

**Suggested fix**: cut on the Python side, `mus = np.array(moments)[:n]` at
`kpmdmrg.py:155` and at `mpsjulialive/dynamics.py:173`, both before the
`kpm_extrapolate` branch (a cut there yields fac*n, the calibrated count for the
caller's delta). That is better than editing the loops: the accelerated
recursion emits moments in pairs, so an exact odd count can only come from a
cut; the loops are shared with `general_kpm`, which `get_distribution()` reaches
with its deliberately uncalibrated count; and one cut covers every backend with
no rebuild, at the cost of one or two wasted MPO applications. Pin
`len(mus) == n` directly rather than a curve, since after the cut the ED-DMRG
residual is set by the `emax` noise. Re-measure both docstrings and the O2 Status
paragraph's "4 to 7 per cent". NUMBERS CHANGE: every DMRG KPM spectrum, by
roughly 2/n onto ED's (7.4 per cent at a band-centre pole at delta=0.2 on 4
sites; 0.5 to 0.7 per cent on a 20-site chain).

### 5. A non-integer `kpm_n_scale` is silently rounded down on `"python"`, ED and `julia_live` (1.5 gives exactly 1x, anything below 1 gives the 16-moment floor at every delta), where v2/v3 raise; and at `kpm_n_scale<=0` the C++ and Python copies of the calibration disagree, 61 moments against 16

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `kpm-calibration`

**Status**: FIXED -- `algebra/kpm.py::validate_kpm_n_scale` accepts `numbers.Integral` only (numpy integers included, `bool` excluded, > 0), raises `TypeError` or `ValueError` naming `kpm_n_scale`, and returns a plain int; it is called once before dispatch in `kpmdmrg.dynamical_correlator_moments` (ahead of `get_gs` and both session calls, which receive the validated int), in the `mode="ED"` push of `Many_Body_Chain.get_dynamical_correlator` (only when the submode is KPM, so an invalid value does not break ED's other submodes, the same reach as on `mode="DMRG"`), at the top of the `julia_live` KPM, and inside `polynomials_for_broadening`, which leaves the C++ `n_scale>0 ? : 1` branch unreachable. Behaviour change: `kpm_n_scale=2.0` now raises on `"python"` and ED, where it used to give 2x, matching v2/v3. Pinned by `tests/test_audit_2026_09_24_kpm.py::test_non_positive_integer_kpm_n_scale_raises_on_dmrg` (`"python"` and v3; 1.5, 0.5, 2.0, 0, -1, True; `match="kpm_n_scale"`, so pybind's own "incompatible function arguments" cannot pass it), `::test_non_positive_integer_kpm_n_scale_raises_on_ed`, `::test_non_positive_integer_kpm_n_scale_raises_on_julia_live`, `::test_kpm_n_scale_is_only_checked_where_it_is_read`, `::test_integer_kpm_n_scale_keeps_its_moment_count_on_dmrg`, `::test_integer_kpm_n_scale_keeps_its_moment_count_on_ed` and `::test_validator_accepts_integers_only`. No integer caller's numbers change.

**Where**: `src/dmrgpy/algebra/kpm.py:469`
(`return max(int(nmin),npol*int(n_scale))`); the C++ copies at
`mpscpp2/chain_session.h:77` and `mpscpp3/chain_session.h:325`
(`npol *= (n_scale>0 ? n_scale : 1)`) and the pybind11 signatures, which take
`int kpm_n_scale`; the readers `kpmdmrg.py:143` to `:154`,
`manybodychain.py:972` (the ED push) and `mpsjulialive/dynamics.py:158`.

`765b537` moved `kpm_n_scale`'s default from 3 to 1 and documents it as a
proportional multiplier on the calibrated count ("larger values buy a sharper
curve at proportionally higher cost"; the user guide's "a multiplier on that
calibrated count, default 1"), nowhere as an integer, and
`examples/dynamical_correlator/hodc_VS_jackson_kernel` even uses it as a real
factor in its own arithmetic. Before `765b537`, `"python"` computed a float count
and its `range(n)` raised `TypeError`, the same as C++, and ED ignored
`kpm_n_scale` altogether; so `765b537` turned a loud failure into a silent one,
and introduced the divergence at non-positive values. No test or example uses a
non-integer or non-positive value.

**Expected**: one rule on every route, either a validated positive integer or a
real multiplier.

Repro, from the hunter:

```python
# kpm_n_scale is documented as a multiplier on the calibrated count.
# What does a non-integer value do on each route?
import numpy as np
from dmrgpy import spinchain
def chain(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=v)
    h = 0
    for i in range(3):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = 20, 20
    return sc
for ns in (1, 2, 0.5, 1.5, 2.5):
    row = []
    for v in (3, "python"):
        sc = chain(v); sc.kpm_n_scale = ns
        try:
            mus, emin, emax, scale, n, d = sc.get_dynamical_correlator_moments(name=[sc.Sz[0], sc.Sz[0]], delta=0.1)
            row.append("%s: n=%d" % (v, n))
        except Exception as e:
            row.append("%s: %s: %s" % (v, type(e).__name__, str(e).splitlines()[0][:60]))
    print("kpm_n_scale=%-4s " % ns + " | ".join(row))
```

```
kpm_n_scale=1    3: n=61 | python: n=61
kpm_n_scale=2    3: n=122 | python: n=122
kpm_n_scale=0.5  3: TypeError: kpm_dynamical_correlator(): incompatible function arguments. | python: n=16
kpm_n_scale=1.5  3: TypeError: kpm_dynamical_correlator(): incompatible function arguments. | python: n=61
kpm_n_scale=2.5  3: TypeError: kpm_dynamical_correlator(): incompatible function arguments. | python: n=122
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced, with ED measured separately (its
spectrum at 1.5 and 1.999 is bit-identical to the one at 1), and the pre-commit
behaviour measured on a `git archive` of `765b537^` (Python-only, so only
`"python"` and ED run there). The reviewer added two things the hunter had not
claimed: the non-positive split (C++ clamps to 1x, Python floors to 16), and the
reverse inconsistency that `2.0` raises on v2/v3 while giving 2x on the Python
routes. The hunter's fix, raising inside `polynomials_for_broadening`, is right
in spirit and wrong in reach: it would not reach the C++ copies, so after it
`kpm_n_scale=0` would raise on `"python"`/ED and quietly give 1x on v3.

```python
# Reviewer probe: what each route does with kpm_n_scale values that are
# not positive integers, measured by the moment count actually used
# (spied on algebra.kpm.polynomials_for_broadening for the ED route, read
# off get_dynamical_correlator_moments for the DMRG routes) and, for ED
# and "python", by whether the returned spectrum differs at all.
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain
from dmrgpy.algebra import kpm as kpmmod

seen = []
_orig = kpmmod.polynomials_for_broadening
def _spy(*a, **k):
    n = _orig(*a, **k); seen.append(n); return n
kpmmod.polynomials_for_broadening = _spy

def chain(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=v)
    h = 0
    for i in range(3):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = 20, 20
    return sc

def dmrg_n(v, ns, delta=0.1):
    sc = chain(v); sc.kpm_n_scale = ns
    try:
        out = sc.get_dynamical_correlator_moments(name=[sc.Sz[0], sc.Sz[0]], delta=delta)
        return "n=%d" % out[4]
    except Exception as e:
        return type(e).__name__

es = np.linspace(-0.5, 4.0, 301)
def ed_curve(ns, delta=0.1):
    sc = chain(3); sc.kpm_n_scale = ns
    del seen[:]
    x, y = sc.get_dynamical_correlator(mode="ED", submode="KPM",
                                       name=[sc.Sz[0], sc.Sz[0]], delta=delta, es=es)
    return seen[-1], np.real(np.array(y))

def py_curve(ns, delta=0.1):
    sc = chain("python"); sc.kpm_n_scale = ns
    x, y = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                       name=[sc.Sz[0], sc.Sz[0]], delta=delta, es=es)
    return np.real(np.array(y))

print("\nA. moment count per route at delta=0.1 (4-site Heisenberg, Sz0 Sz0)")
ref_ed = {}
for ns in (1, 2, 2.0, 0.5, 1.5, 1.999, 0, -1, np.int64(2)):
    nED, yED = ed_curve(ns)
    ref_ed[repr(ns)] = yED
    print("kpm_n_scale=%-12r v2: %-10s v3: %-10s python: %-10s ED: n=%d"
          % (ns, dmrg_n(2, ns), dmrg_n(3, ns), dmrg_n("python", ns), nED))

print("\nB. ED spectrum identical to the kpm_n_scale=1 one, bit for bit?")
for k in ("1.5", "1.999", "2", "0.5", "0", "-1"):
    print("  %-6s max|y - y(1)| = %.3e" % (k, np.max(np.abs(ref_ed[k] - ref_ed["1"]))))
print("  python route: max|y(1.5)-y(1)| = %.3e, max|y(2)-y(1)| = %.3e"
      % (np.max(np.abs(py_curve(1.5) - py_curve(1))), np.max(np.abs(py_curve(2) - py_curve(1)))))

print("\nC. kpm_n_scale=0.5 on ED across delta: count used vs count 0.5x asks for")
for d in (0.2, 0.1, 0.05, 0.02):
    n05, y05 = ed_curve(0.5, delta=d)
    n1, y1 = ed_curve(1, delta=d)
    print("  delta=%-5s n(1)=%-4d n(0.5)=%-4d (0.5x would be ~%d)  peak(1)=%.4f peak(0.5)=%.4f"
          % (d, n1, n05, int(round(0.5*n1)), y1.max(), y05.max()))
```

```
dmrgpy from <repo>/src/dmrgpy/__init__.py

A. moment count per route at delta=0.1 (4-site Heisenberg, Sz0 Sz0)
kpm_n_scale=1            v2: n=61       v3: n=61       python: n=61       ED: n=61
kpm_n_scale=2            v2: n=122      v3: n=122      python: n=122      ED: n=122
kpm_n_scale=2.0          v2: TypeError  v3: TypeError  python: n=122      ED: n=122
kpm_n_scale=0.5          v2: TypeError  v3: TypeError  python: n=16       ED: n=16
kpm_n_scale=1.5          v2: TypeError  v3: TypeError  python: n=61       ED: n=61
kpm_n_scale=1.999        v2: TypeError  v3: TypeError  python: n=61       ED: n=61
kpm_n_scale=0            v2: n=61       v3: n=61       python: n=16       ED: n=16
kpm_n_scale=-1           v2: n=61       v3: n=61       python: n=16       ED: n=16
kpm_n_scale=np.int64(2)  v2: n=122      v3: n=122      python: n=122      ED: n=122

B. ED spectrum identical to the kpm_n_scale=1 one, bit for bit?
  1.5    max|y - y(1)| = 0.000e+00
  1.999  max|y - y(1)| = 0.000e+00
  2      max|y - y(1)| = 8.150e-01
  0.5    max|y - y(1)| = 6.017e-01
  0      max|y - y(1)| = 6.017e-01
  -1     max|y - y(1)| = 6.017e-01
  python route: max|y(1.5)-y(1)| = 7.772e-15, max|y(2)-y(1)| = 8.276e-01

C. kpm_n_scale=0.5 on ED across delta: count used vs count 0.5x asks for
  delta=0.2   n(1)=31   n(0.5)=16   (0.5x would be ~16)  peak(1)=0.4304 peak(0.5)=0.2318
  delta=0.1   n(1)=61   n(0.5)=16   (0.5x would be ~30)  peak(1)=0.8335 peak(0.5)=0.2318
  delta=0.05  n(1)=123  n(0.5)=16   (0.5x would be ~62)  peak(1)=1.6617 peak(0.5)=0.2318
  delta=0.02  n(1)=306  n(0.5)=16   (0.5x would be ~153)  peak(1)=4.0197 peak(0.5)=0.2318
EXIT=0
```

```python
# Reviewer probe: the same non-integer kpm_n_scale values on the tree
# BEFORE 765b537 (fd01679, git-archived, no compiled extension), on the
# two routes that tree can run: "python" DMRG and mode="ED".
import numpy as np
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
from dmrgpy import spinchain

def chain(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version=v)
    h = 0
    for i in range(3):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm, sc.nsweeps = 20, 20
    return sc

es = np.linspace(-0.5, 4.0, 301)
ed = {}
for ns in (3, 1, 0.5, 1.5, 2.5):
    sc = chain("python"); sc.kpm_n_scale = ns
    try:
        out = sc.get_dynamical_correlator_moments(name=[sc.Sz[0], sc.Sz[0]], delta=0.1)
        p = "n=%r" % (out[4],)
    except Exception as e:
        p = "%s: %s" % (type(e).__name__, str(e).splitlines()[0][:70])
    sc = chain("python"); sc.kpm_n_scale = ns
    x, y = sc.get_dynamical_correlator(mode="ED", submode="KPM",
                                       name=[sc.Sz[0], sc.Sz[0]], delta=0.1, es=es)
    ed[ns] = np.real(np.array(y))
    print("kpm_n_scale=%-4r python: %-70s ED max|y-y(3)|=%.1e"
          % (ns, p, np.max(np.abs(ed[ns] - ed[3]))))
```

```
dmrgpy from <scratch>/review-kpm-C/old/src/dmrgpy/__init__.py
kpm_n_scale=3    python: n=72                                                                   ED max|y-y(3)|=0.0e+00
kpm_n_scale=1    python: n=24                                                                   ED max|y-y(3)|=0.0e+00
kpm_n_scale=0.5  python: TypeError: 'float' object cannot be interpreted as an integer          ED max|y-y(3)|=0.0e+00
kpm_n_scale=1.5  python: TypeError: 'float' object cannot be interpreted as an integer          ED max|y-y(3)|=0.0e+00
kpm_n_scale=2.5  python: TypeError: 'float' object cannot be interpreted as an integer          ED max|y-y(3)|=0.0e+00
```

**Suggested fix**: validate the attribute once, before dispatch, where it is
read off the chain (the documentation.md section 4.10 shape): one shared helper
rejecting anything that is not a positive integer (`numbers.Integral`, not
`bool`, >0), called from `kpmdmrg.dynamical_correlator_moments` ahead of both
session calls, from the ED push at `manybodychain.py:972` and from the
`julia_live` read, with `polynomials_for_broadening` calling it too. That unifies
every backend without a rebuild and leaves the C++ `>0 ? : 1` branch
unreachable. Accepting a real multiplier (`int(round(npol*n_scale))` on both
sides) is defensible but needs both bindings to take a double and a rebuild.
Numbers do not change for any integer caller.

### 6. The documentation never says that a ground-state correlator's weight sits where the calibrated Jackson line is 0.70 of the requested width, so `delta` is typically delivered as 0.7 to 0.9 of itself on the default window and 0.3 to 0.7 on the ground-state-anchored one

`docs` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `kpm-calibration`

**Status**: FIXED -- in the documentation, as the finding asked, with no code change: `algebra/kpm.py::polynomials_for_broadening`'s docstring, `dynamics.py`'s module docstring, the user guide's KPM entry (both formats), `documentation.{md,tex}`'s "delta is one broadening" paragraph, `CLAUDE.md`'s O2 paragraph and a pointer sentence in the O2 Status paragraph of `audit_2026_09_hole_hunt.md` now say where the calibrated width holds: E0 at x0 = -1/(2*kpm_scale) on every chain, the band centre at E0 + W/2, a ground-state correlator's width tending to 0.70 of the requested one (0.157 at E0 on the anchored window), the measured 12- and 20-site numbers, and the compensation delta/sqrt(1-x(omega)^2) with x(omega) = x0 + omega/(kpm_scale*W); the 1.6x peak-ratio sentence is now about 1.5x (1.514, the kernel's own ratio at equal FWHM), with the 0.315/0.195 = 1.61 read as including the off-centre narrowing of its dominant pole at x=-0.527. No numbers change.

**Where**: `algebra/kpm.py::polynomials_for_broadening`'s docstring,
`dynamics.py`'s module docstring, `docs/user_guide.{md,tex}`'s KPM section,
`docs/documentation.{md,tex}`, `CLAUDE.md`'s O2 paragraph and the O2 Status
paragraph of `audit_2026_09_hole_hunt.md`. No code defect.

The O2 calibration makes the Jackson line FWHM = 2*delta at x=0 of the Chebyshev
window, which on the default bandwidth-centred window is the energy E0 + W/2.
The ground state itself sits at x0 = -1/(2*kpm_scale) = -0.714 on every chain,
exactly, since `shift = -(emin+emax)/2` and `scale = 1/((emax-emin)*kpm_scale)`
are the same on all five routes, and an excitation at omega sits at
x0 + omega/(kpm_scale*W). So as W grows extensively every intensive excitation of
a ground-state correlator comes out at FWHM 2*delta*sqrt(1-x^2), which tends to
2*delta*sqrt(1-1/(4*kpm_scale^2)) = 0.70 x 2*delta. The documentation states the
`sqrt(1-x^2)` profile, and in the energy-truncation section says a correlator's
weight "sits just above E0", but never puts the two together, and its phrasing
("exact at the band centre and tightens as sqrt(1-x^2) towards the band edges")
presents the off-centre case as the exception, where for a ground-state
correlator it is the rule.

**Expected**: documentation that says where the calibrated width holds and what a
ground-state correlator actually gets.

Repro, from the hunter, on a toy that slides one pole towards the lower edge:

```
centred window (default), kpm_scale=0.7, delta=0.10
    hz       W  x_pole sqrt(1-x^2)     ED/2d   v3[:n]/2d     v3/2d
  0.00    2.00   0.000       1.000     0.988       0.988     0.952
  1.45    3.45  -0.300       0.954     0.956       0.957     0.946
  4.67    6.67  -0.500       0.866     0.868       0.868     0.863
 10.50   12.50  -0.600       0.800     0.805       0.805     0.800
 20.20   22.20  -0.650       0.760     0.766       0.766     0.765
 40.00   42.00  -0.680       0.733     0.740       0.740     0.739
ground-state-anchored window (kpm_energy_truncate=True), v3
    hz    kpm_sc       W  x_pole sqrt(1-x^2)   v3[:n]/2d     v3/2d
  4.67      0.70    6.67  -0.564       0.825       0.829     0.819   (n=87, weight=0.2500)
 10.50      0.70   12.50  -0.762       0.648       0.650     0.642   (n=164, weight=0.2500)
 20.20      0.70   22.20  -0.860       0.510       0.513     0.511   (n=291, weight=0.2500)
 20.20      0.30   22.20  -0.691       0.723       0.723     0.717   (n=125, weight=0.2500)
```

**Reviewer (CONFIRMED, NARROWED)**: the profile holds to about 1 per cent at
every position, on the real chain too, and every KPM route does exactly what its
docstrings say (dmrgpy's 12-site curve equals an independent Jackson
reconstruction from the exact poles to 5.4e-15). Measured on open Heisenberg
chains, `<Sz_i;Sz_i>`: the isolated lowest pole at 0.755 x 2*delta on 12 sites
(`mode="ED"`) and 0.725 on 20 sites (v3, `mus[:n]`), against predicted 0.748 and
0.719; weight-averaged, 0.87 to 0.90 on 8 sites, 0.83 to 0.86 on 12, 0.78 to
0.80 on 20; and on 12 sites dmrgpy's curve matches the exact Lehmann sum
broadened by per-pole Gaussians of FWHM 2*delta*sqrt(1-x_n^2) to 3.2 per cent of
its peak, while it misses Gaussians of the requested width by 26.6 per cent. On
the anchored window (`kpm_energy_truncate=True`) E0 sits at -0.9875 by
construction and the factor there is 0.157; the lowest pole of the 12-site
chain comes out at 0.464 x 2*delta (predicted 0.462). That half is not the open
known issue, which is weight clipped out of the window, since only 5.4e-10 of
0.25 lies above it here. Struck: the hunter's "(roughly a 20-site Heisenberg
chain)" for the toy's 0.805; the real 20-site chain gives 0.72 on its lowest
pole. Narrowed: "no ground-state weight lies at the band centre" becomes "a
vanishing fraction", since the local weight is bounded by roughly pi J while
W/2 grows as about 0.34 L J, crossing near 20 sites. A side point in the same
passage: the documented "at equal FWHM and equal integrated weight its peak
stands about 1.6x higher than a resolvent submode's" was not measured at equal
FWHM (its dominant pole sits at x=-0.527, where the line is 0.850 of 2*delta);
the Jackson kernel's own ratio at equal FWHM is 1.514. The reviewer's probes
(they import `exactlehmann.py`, in "Shared helpers"):

```python
# Reviewer probe 2: what width does dmrgpy's default KPM actually deliver on
# a real 12-site open Heisenberg chain, <Sz_6;Sz_6>, mode="ED" (so the
# n-versus-n+2 moment question of the DMRG routes stays out)?
#  (a) pointwise: dmrgpy's curve against the exact Lehmann sum broadened by
#      a Gaussian of FWHM 2*delta (what was requested) and by one of FWHM
#      2*delta*sqrt(1-x_n^2) per pole (what the documented profile predicts
#      at the pole's own x), plus an independent Jackson reconstruction
#      from the exact poles with dmrgpy's own moment count;
#  (b) the isolated lowest pole's FWHM at delta=0.05, against 2*delta.
import sys, numpy as np
sys.path.insert(0, "<scratch>/review-kpm-H4")
from exactlehmann import sector, szop, lanczos_poles
import scipy.sparse.linalg as sla
from dmrgpy import spinchain
from dmrgpy.algebra.kpm import polynomials_for_broadening

L, site, ks = 12, 6, 0.7
st, H = sector(L)
e0, g = sla.eigsh(H, k=1, which="SA", tol=1e-13); e0 = e0[0]; g = g[:, 0]
W = 0.25*(L-1) - e0
e, wts = lanczos_poles(H, szop(st, site)*g, 160)
om = e - e0; keep = wts > 1e-14; om, wts = om[keep], wts[keep]
xn = (om - W/2)/(W*ks); fn = np.sqrt(1-xn**2)

def chain():
    sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=3)
    h = 0
    for i in range(L-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc

def jackson_exact(es, delta):
    N = polynomials_for_broadening(W*ks, delta)
    k = np.arange(N); q = np.pi/(N+1)
    gk = ((N-k+1)*np.cos(q*k) + np.sin(q*k)/np.tan(q))/(N+1)
    mu = np.array([np.sum(wts*np.cos(kk*np.arccos(xn))) for kk in k])
    x = (es - W/2)/(W*ks)
    T = np.cos(np.outer(k, np.arccos(x)))
    c = gk*mu; c[1:] *= 2
    return (c @ T)/(np.pi*np.sqrt(1-x**2))/(W*ks), N

def gauss(es, fw):
    s = fw/(2*np.sqrt(2*np.log(2)))
    return np.sum(wts[:, None]*np.exp(-(es[None, :]-om[:, None])**2/(2*s[:, None]**2))
                  /(np.sqrt(2*np.pi)*s[:, None]), axis=0)

def lorentz(es, delta):
    return np.sum(wts[:, None]*delta/np.pi/((es[None, :]-om[:, None])**2+delta**2), axis=0)

print("L=%d, site %d, E0=%.6f, W=%.4f, centred window kpm_scale=%.1f" % (L, site, e0, W, ks))
es = np.linspace(0.0, 4.0, 2001)
for delta in (0.1, 0.2):
    sc = chain()
    _, y = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=[sc.Sz[site], sc.Sz[site]],
                                       delta=delta, es=es)
    y = np.real(y)
    yj, N = jackson_exact(es, delta)
    yreq = gauss(es, 2*delta*np.ones_like(om))
    ypred = gauss(es, 2*delta*fn)
    ylor = lorentz(es, delta)
    pk = np.max(y)
    print("delta=%.2f  npol=%d  dmrgpy peak %.4f  weight %.6f" % (delta, N, pk, np.trapezoid(y, es)))
    print("   max|dmrgpy - Jackson from exact poles|         = %.2e" % np.max(np.abs(y-yj)))
    print("   max|dmrgpy - Gauss FWHM 2*delta|               = %.4f  (%.1f%% of peak)" % (np.max(np.abs(y-yreq)), 100*np.max(np.abs(y-yreq))/pk))
    print("   max|dmrgpy - Gauss FWHM 2*delta*sqrt(1-x_n^2)| = %.4f  (%.1f%% of peak)" % (np.max(np.abs(y-ypred)), 100*np.max(np.abs(y-ypred))/pk))
    print("   max|dmrgpy - Lorentz FWHM 2*delta (resolvent)| = %.4f  (%.1f%% of peak; resolvent peak %.4f)" % (np.max(np.abs(y-ylor)), 100*np.max(np.abs(y-ylor))/pk, np.max(ylor)))

def fwhm(x, y):
    k = np.argmax(y); h = y[k]/2.; l = k
    while y[l] > h: l -= 1
    r = k
    while y[r] > h: r += 1
    xl = x[l] + (h-y[l])*(x[l+1]-x[l])/(y[l+1]-y[l])
    xr = x[r-1] + (h-y[r-1])*(x[r]-x[r-1])/(y[r]-y[r-1])
    return xr-xl

print("isolated lowest pole, omega=%.4f, x=%.4f, sqrt(1-x^2)=%.4f, next weighted pole at omega=%.4f"
      % (om[0], xn[0], fn[0], om[1]))
es2 = np.linspace(0.0, 0.6, 3001)
for delta in (0.03, 0.05):
    sc = chain()
    _, y = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=[sc.Sz[site], sc.Sz[site]],
                                       delta=delta, es=es2)
    y = np.real(y)
    print("delta=%.2f  FWHM=%.5f  FWHM/(2*delta)=%.4f  predicted sqrt(1-x^2)=%.4f"
          % (delta, fwhm(es2, y), fwhm(es2, y)/(2*delta), fn[0]))
```

```
L=12, site 6, E0=-5.142091, W=7.8921, centred window kpm_scale=0.7
delta=0.10  npol=204  dmrgpy peak 0.3605  weight 0.249959
   max|dmrgpy - Jackson from exact poles|         = 5.44e-15
   max|dmrgpy - Gauss FWHM 2*delta|               = 0.0959  (26.6% of peak)
   max|dmrgpy - Gauss FWHM 2*delta*sqrt(1-x_n^2)| = 0.0116  (3.2% of peak)
   max|dmrgpy - Lorentz FWHM 2*delta (resolvent)| = 0.1726  (47.9% of peak; resolvent peak 0.2122)
delta=0.20  npol=102  dmrgpy peak 0.1813  weight 0.249631
   max|dmrgpy - Jackson from exact poles|         = 3.66e-15
   max|dmrgpy - Gauss FWHM 2*delta|               = 0.0460  (25.4% of peak)
   max|dmrgpy - Gauss FWHM 2*delta*sqrt(1-x_n^2)| = 0.0071  (3.9% of peak)
   max|dmrgpy - Lorentz FWHM 2*delta (resolvent)| = 0.0764  (42.1% of peak; resolvent peak 0.1350)
isolated lowest pole, omega=0.2809, x=-0.6634, sqrt(1-x^2)=0.7482, next weighted pole at omega=0.6288
delta=0.03  FWHM=0.04528  FWHM/(2*delta)=0.7546  predicted sqrt(1-x^2)=0.7482
delta=0.05  FWHM=0.07532  FWHM/(2*delta)=0.7532  predicted sqrt(1-x^2)=0.7482
```

```python
# Reviewer probe 4 (side check): the documented "at equal FWHM and equal
# weight the Jackson peak stands about 1.6x higher than a Lorentzian",
# measured 0.315 against 0.195 on a 6-site Heisenberg chain at delta=0.2.
# Is that ratio at equal FWHM, or does it already contain the off-centre
# narrowing this candidate is about?
import sys, numpy as np
sys.path.insert(0, "<scratch>/review-kpm-H4")
from exactlehmann import sector, szop, lanczos_poles
import scipy.sparse.linalg as sla
from dmrgpy import spinchain

L, delta, ks = 6, 0.2, 0.7
st, H = sector(L)
e0, g = sla.eigsh(H, k=1, which="SA", tol=1e-13); e0 = e0[0]; g = g[:, 0]
W = 0.25*(L-1) - e0
sc = spinchain.Spin_Chain(["S=1/2"]*L, itensor_version=3)
h = 0
for i in range(L-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
es = np.linspace(0.0, 4.0, 4001)
for site in (0,):
    e, w = lanczos_poles(H, szop(st, site)*g, 20)
    om = e - e0; k = w > 1e-12; om, w = om[k], w[k]
    x = (om - W/2)/(W*ks)
    print("L=6 site %d, W=%.4f; poles (omega, weight, x, sqrt(1-x^2)):" % (site, W))
    for a, b, c in zip(om, w, x): print("   %.4f %.5f %.4f %.4f" % (a, b, c, np.sqrt(1-c*c)))
    _, yk = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=[sc.Sz[site], sc.Sz[site]], delta=delta, es=es)
    _, yi = sc.get_dynamical_correlator(mode="ED", submode="INV", name=[sc.Sz[site], sc.Sz[site]], delta=delta, es=es)
    print("   KPM peak %.4f, INV peak %.4f, ratio %.3f" % (np.max(np.real(yk)), np.max(np.real(yi)), np.max(np.real(yk))/np.max(np.real(yi))))
# the kernel's own ratio at equal FWHM, one pole at x=0, large N
for N in (26, 100, 400):
    kk = np.arange(N); q = np.pi/(N+1)
    gk = ((N-kk+1)*np.cos(q*kk) + np.sin(q*kk)/np.tan(q))/(N+1)
    xs = np.linspace(-8.0/N, 8.0/N, 20001)
    c = gk*np.cos(kk*np.pi/2); c[1:] *= 2
    y = (c @ np.cos(np.outer(kk, np.arccos(xs))))/(np.pi*np.sqrt(1-xs**2))
    hm = y.max()/2; ab = xs[y > hm]; fw = ab[-1]-ab[0]
    print("N=%3d single pole at x=0: Jackson peak*FWHM = %.4f, Lorentzian peak*FWHM = 2/pi = %.4f, ratio at equal FWHM %.3f (Gaussian: %.3f)"
          % (N, y.max()*fw, 2/np.pi, y.max()*fw/(2/np.pi), np.sqrt(np.pi*np.log(2))))
```

```
L=6 site 0, W=3.7436; poles (omega, weight, x, sqrt(1-x^2)):
   0.4916 0.11042 -0.5267 0.8501
   1.0721 0.09910 -0.3052 0.9523
   1.4794 0.03343 -0.1497 0.9887
   1.6208 0.00444 -0.0958 0.9954
   2.0082 0.00249 0.0521 0.9986
   2.5366 0.00012 0.2537 0.9673
   2.8759 0.00000 0.3832 0.9237
   3.1504 0.00000 0.4879 0.8729
   3.4572 0.00000 0.6050 0.7962
   KPM peak 0.3148, INV peak 0.1950, ratio 1.614
N= 26 single pole at x=0: Jackson peak*FWHM = 0.9669, Lorentzian peak*FWHM = 2/pi = 0.6366, ratio at equal FWHM 1.519 (Gaussian: 1.476)
N=100 single pole at x=0: Jackson peak*FWHM = 0.9638, Lorentzian peak*FWHM = 2/pi = 0.6366, ratio at equal FWHM 1.514 (Gaussian: 1.476)
N=400 single pole at x=0: Jackson peak*FWHM = 0.9636, Lorentzian peak*FWHM = 2/pi = 0.6366, ratio at equal FWHM 1.514 (Gaussian: 1.476)
```

**Suggested fix**: documentation only, in the places listed above: on the centred
window E0 sits at x0 = -1/(2*kpm_scale) on every chain and the band centre is
E0 + W/2; a ground-state correlator's width therefore tends to 0.70 of the
requested one, and to 0.157 at E0 on the anchored window, with the measured
12- and 20-site numbers; to get FWHM 2*delta at a chosen omega, pass
delta/sqrt(1-x(omega)^2), the only compensation a user has since
`kpm_n_scale` cannot ask for fewer moments; and the peak-ratio sentence says
about 1.5x. Two code options exist and are recorded rather than recommended:
calibrating at x0 (0.70 times the moments, but everything above comes out up
to 1.43x broader, giving up the "never broader than requested" property the
`nmin` floor's docstring relies on, and nonsensical on the anchored window), or
at the correlator's own centroid x = mu_1/mu_0 (one extra Hamiltonian
application, a correlator-dependent count, and every KPM curve changing again
one release after O2).

### 7. The 3.31e-02 TDZ residual that `765b537` attributes to the contour is an FFT-grid interpolation error of 14 per cent of the peak, shared with TD at `predict=False`; the contour's own share is 6.2e-07, and the documented remedy (`alpha0`, `n_max`) does nothing

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED &middot; lens `td-convention`

**Status**: FIXED -- `_fourier_transform_correlator` evaluates the damped trapezoid sum directly at each requested `es`, chunked over `es` (`_damped_sum_at`, at most 2^21 elements per chunk), with the FFT's own time origin t_k = k*dtnew; frequencies outside the old FFT band still return 0, as the interpolation did there. At a frequency on the FFT grid the direct sum and the FFT agree to 3.9e-16 to 1.3e-15 (predict False and True, nt even and odd). `sxt_to_skomega` stays on the old stage through a private `_evaluation="fft"` keyword and returns the same bits as before, and so do the two direct callers at short windows, `tests/test_infinite_chain.py` and `examples/idmrg/td_dynamical_correlator` (delta*T of 0.375 and 0.9, where the direct stage moved the test's TD peak to the window edge), which now pass `_evaluation="fft"` explicitly until that regime is measured. `test_td_linear_prediction_sharpens_peak` had been measuring the grid: at its default nt=1200 (delta*T=6) the old stage gave FWHM 0.1275 without prediction and 0.0975 with it, the direct stage 0.0975 both ways, the exact width; prediction genuinely narrows where truncation sets the width, 0.2125 to 0.0975 at nt=200 (delta*T=1), so the test now runs there with its assertions kept, and `dynamical_correlator`'s docstring says the same. Pinned by `tests/test_audit_2026_09_24_realtime.py::test_real_time_routes_are_held_to_their_own_accuracy`, `::test_the_contour_error_on_its_own`, `::test_direct_evaluation_is_the_fft_on_its_own_grid` and `::test_sxt_to_skomega_stays_on_the_fft_stage`. NUMBERS CHANGE, on the audit chain (4-site `Fermionic_Chain`, seed 3, A=Cdag_0, B=C_2, es=linspace(-1,6,60), dt=0.1, v3, exact peak 0.2313), max|y-exact| at delta=0.4: TDZ default 3.314e-02 to 5.408e-04; TD `predict=False` 3.314e-02 to 5.398e-04 on DMRG and to 5.393e-04 on ED; TD default 2.571e-04 to 3.169e-05 on DMRG and 2.567e-04 to 3.145e-05 on ED; the returned curves move by up to 3.268e-02 (14 per cent of the peak) for TDZ and TD `predict=False` and 2.89e-04 for TD default; at delta=0.25 (peak 0.3656) TDZ goes 6.994e-02 to 6.263e-04 and TD default 5.753e-04 to 2.040e-05, a move of up to 19 per cent of the peak. The contour on its own, max|y_TDZ - y_TD(predict=False)|, is 9.85e-07.

**Where**: `src/dmrgpy/timedependent.py::_fourier_transform_correlator` (the FFT
on the grid of spacing 2*pi/(nt*dt) and the linear interpolation onto `es`),
shared by `submode="TD"` and `submode="TDZ"`. The misattribution sits in the O1
Status paragraph of `audit_2026_09_hole_hunt.md` (near line 2228), the
`765b537` commit message, `src/dmrgpy/dynamics.py:129`,
`docs/user_guide.md:2013` to `2018` and `:3773`, `docs/user_guide.tex:2299` and
`:5296`, and the comment, docstring and 6e-2 tolerance of
`test_real_time_submodes_return_the_complex_lehmann_density`
(`tests/test_audit_2026_09_correlator-conventions.py:343` to `352`).

With `predict=False`, TDZ's default and TD's documented opt-out, the damped
trapezoid sum is evaluated only on the FFT grid, whose spacing is
2*pi/(nt*dt) = 2*pi*delta/damping_periods, 1.047*delta at the default
`damping_periods=6`, and is then interpolated linearly onto `es`, meaning that a
Lorentzian of half-width delta gets about one sample per delta. On the audit's own
complex-hopping chain that alone accounts for the recorded 3.31e-02, and TD with
`predict=False`, which has no contour at all, gives the same number on DMRG and
on ED. It survived because TD defaults to `predict=True`, whose tenfold linear
prediction puts the grid at 0.105*delta, while TDZ defaults to `predict=False`;
and because the TDZ test sits at 6e-2 against exact, about 110 times its real
accuracy, so it cannot see a contour regression either. This is the second TDZ
residual the records have blamed on the reconstruction: the 2026-08 record's
"further ~2x" turned out to be `DoNormalize` (2026-09 finding #7).

**Expected**: a frequency evaluation whose error is set by the finite time
window and the method, not by how many grid points fall inside one line width.

Repro, from the hunter (its helper `leadE_fft_grid_helpers.py` is in "Shared
helpers"):

```python
"""Candidate E: with predict=False (TDZ's default, and TD's documented
opt-out) the spectrum is an FFT on a grid of spacing 2*pi/(nt*dt) =
2*pi*delta/damping_periods ~ 1.05*delta, linearly interpolated onto es,
which is a Lorentzian of HWHM delta sampled once per ~delta."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from dmrgpy import fermionchain

def complex_hopping_chain(n=4, itensor_version=3, seed=3):
    # verbatim from tests/test_audit_2026_09_correlator-conventions.py
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    rng = np.random.RandomState(seed)
    t = rng.random((n, n)) + 1j * rng.random((n, n))
    t = t + t.conj().T
    h = 0
    for i in range(n):
        for j in range(n): h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1): h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 30, 12
    return fc

def lehmann(chain, A, B):
    ed = chain.get_ED_obj()
    h = np.array(ed.get_hamiltonian().todense())
    emu, vs = np.linalg.eigh(h)
    U = np.array(vs); Uh = np.conjugate(U.T)
    Ae = Uh @ np.array(ed.MO2matrix(A).todense()) @ U
    Be = Uh @ np.array(ed.MO2matrix(B).todense()) @ U
    return emu - emu[0], Ae[0, :] * Be[:, 0]

def density(D, M, es, delta):
    return np.array([np.sum(M * (delta / np.pi) / ((w - D) ** 2 + delta ** 2))
                     for w in es])

DELTA, DT = 0.4, 0.1
nt = int(6 / DELTA / DT); dw = 2 * np.pi / (nt * DT)
fc = complex_hopping_chain()
A, B = fc.Cdag[0], fc.C[2]
D, M = lehmann(fc, A, B)
grids = [("audit es=linspace(-1,6,60)", np.linspace(-1.0, 6.0, 60)),
         ("es on the FFT grid k*%.4f" % dw, dw * np.arange(-2, 15))]
print("nt=%d  T=%.1f  FFT spacing 2pi/T = %.4f = %.2f*delta" % (nt, nt*DT, dw, dw/DELTA))
for label, ES in grids:
    ref = density(D, M, ES, DELTA)
    print("\n%s  (exact peak %.4f)" % (label, np.max(np.abs(ref))))
    for mode, sub, extra in (("DMRG", "TDZ", {}), ("DMRG", "TDZ", dict(predict=True)),
                             ("DMRG", "TD", dict(predict=False)), ("DMRG", "TD", {}),
                             ("ED", "TD", dict(predict=False))):
        _x, y = fc.get_dynamical_correlator(mode=mode, submode=sub, name=[A, B],
                    es=ES, delta=DELTA, dt=DT, **extra)
        print("  %-4s %-3s %-15s max|y-exact| = %.3e"
              % (mode, sub, str(extra) if extra else "(defaults)",
                 np.max(np.abs(np.asarray(y) - ref))))
```

```
nt=150  T=15.0  FFT spacing 2pi/T = 0.4189 = 1.05*delta

audit es=linspace(-1,6,60)  (exact peak 0.2313)
  DMRG TDZ (defaults)      max|y-exact| = 3.314e-02
  DMRG TDZ {'predict': True} max|y-exact| = 2.570e-04
  DMRG TD  {'predict': False} max|y-exact| = 3.314e-02
  DMRG TD  (defaults)      max|y-exact| = 2.570e-04
  ED   TD  {'predict': False} max|y-exact| = 3.314e-02

es on the FFT grid k*0.4189  (exact peak 0.2230)
  DMRG TDZ (defaults)      max|y-exact| = 3.457e-04
  DMRG TDZ {'predict': True} max|y-exact| = 3.147e-05
  DMRG TD  {'predict': False} max|y-exact| = 3.452e-04
  DMRG TD  (defaults)      max|y-exact| = 3.143e-05
  ED   TD  {'predict': False} max|y-exact| = 3.450e-04
```

```python
"""Candidate E, second half: the residual is set by the FFT spacing
2*pi/(nt*dt) = 2*pi*delta/damping_periods, not by dt or by the contour."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from leadE_fft_grid_helpers import complex_hopping_chain, lehmann, density
DELTA = 0.4
ES = np.linspace(-1.0, 6.0, 60)
fc = complex_hopping_chain()
A, B = fc.Cdag[0], fc.C[2]
D, M = lehmann(fc, A, B)
ref = density(D, M, ES, DELTA)
print("exact peak %.4f" % np.max(np.abs(ref)))
for dt, dp in ((0.1, 6), (0.05, 6), (0.1, 12), (0.1, 24)):
    T = int(dp/DELTA/dt)*dt
    _x, y = fc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=[A, B], es=ES,
                                        delta=DELTA, dt=dt, damping_periods=dp)
    print("TDZ dt=%.2f damping_periods=%2d  FFT spacing=%.3f*delta  max|y-exact| = %.3e"
          % (dt, dp, 2*np.pi/T/DELTA, np.max(np.abs(np.asarray(y) - ref))), flush=True)
```

```
exact peak 0.2313
TDZ dt=0.10 damping_periods= 6  FFT spacing=1.047*delta  max|y-exact| = 3.313e-02
TDZ dt=0.05 damping_periods= 6  FFT spacing=1.047*delta  max|y-exact| = 3.315e-02
TDZ dt=0.10 damping_periods=12  FFT spacing=0.524*delta  max|y-exact| = 1.424e-02
TDZ dt=0.10 damping_periods=24  FFT spacing=0.262*delta  max|y-exact| = 1.794e-03
```

**Reviewer (CONFIRMED)**: the 3.31e-02 is entirely the grid. Isolated with no
reference and no patch, max|y_TDZ - y_TD(predict=False)| on the same `es` is
2.3e-03, 1.2e-04, 9.7e-06 and 6.2e-07 for `n_max` 1 to 4 (4 is the default), and
TDZ's reconstructed C(t) lies within 5e-06 of the exact one. Evaluating the same
trapezoid sum directly at each `es` takes TDZ from 3.313e-02 to 5.4e-04, the
finite-T floor it shares with TD `predict=False` and ED, and on the FFT grid
direct evaluation and the FFT agree to 5.6e-17. The residual does not move with
dt, falls with the spacing, and is 6.99e-02 (19 per cent) at delta=0.25 on the
same chain. The user guide's remedy, "reduce it with a smaller `alpha0` or a
larger `n_max`", does nothing: `alpha0` of 0.05, 0.2 and 0.3 and `n_max=2` all
leave 3.312e-02 to 3.324e-02, and `n_max` above `_MAX_SUPPORTED_ORDER=4` raises
`NotImplementedError`. The record's Hermitian-pair 1.8e-2 has the same cause:
on its own configuration (6-site Heisenberg, A=B=Sz_0, delta=0.3, dt=0.05,
nt=400, T=20) TDZ and TD `predict=False` both give 1.830e-02 and both drop to
2.2e-04 by direct evaluation; the record had compared TDZ against TD on
`predict=True`, a different frequency pipeline. Corrections: "about 200x looser"
becomes about 110x; ED's evolution is `solve_ivp` RK45 at scipy's default
tolerances, 3.2e-08 from the exact C(t), exact at this level rather than exact
by construction. "Changes numbers towards exact for every `predict=False`
caller" is corrected twice: the fix also moves `predict=True` callers, TD's
default included (2.565e-04 to 3.1e-05 here), and "towards exact" is struck for
`sxt_to_skomega`, which is on `predict=False` but takes the caller's `nt`: at
`td_dynamical_correlator`'s defaults delta*T=1 and the grid is 6.3*delta, so the
series is cut at e^-1 and a direct evaluation would expose truncation ringing
the interpolation now smooths away (not measured; no test pins S(k,omega)
values). The reviewer's probes (`revhelpers.py` and
`probe1_complex_chain_build.py` are in "Shared helpers"):

```python
"""Reviewer probe 3: isolate the contour + Taylor-in-alpha0 error with no
reference and no patch, as max|y_TDZ - y_TD(predict=False)| on the same
es (the two share the whole FFT tail), across n_max and alpha0. Shows how
large the contour's own error is, and that a TDZ-vs-TD comparison sees it
where the pinned 6e-2 against exact cannot."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe1_complex_chain_build import build

ES = np.linspace(-1.0, 6.0, 60)
DELTA, DT = 0.4, 0.1
fc, D, M, _gap = build()
A, B = fc.Cdag[0], fc.C[2]
_x, ytd = fc.get_dynamical_correlator(mode="DMRG", submode="TD", name=[A, B],
                                      es=ES, delta=DELTA, dt=DT, predict=False)
ytd = np.asarray(ytd)
for kw in (dict(n_max=1), dict(n_max=2), dict(n_max=3), dict(n_max=4),
           dict(alpha0=0.3, n_max=2), dict(alpha0=0.3, n_max=4)):
    _x, y = fc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=[A, B],
                                        es=ES, delta=DELTA, dt=DT, **kw)
    print("TDZ %-28s max|y_TDZ - y_TD(predict=False)| = %.3e"
          % (str(kw), np.max(np.abs(np.asarray(y) - ytd))), flush=True)
```

```
TDZ {'n_max': 1}                 max|y_TDZ - y_TD(predict=False)| = 2.319e-03
TDZ {'n_max': 2}                 max|y_TDZ - y_TD(predict=False)| = 1.194e-04
TDZ {'n_max': 3}                 max|y_TDZ - y_TD(predict=False)| = 9.666e-06
TDZ {'n_max': 4}                 max|y_TDZ - y_TD(predict=False)| = 6.234e-07
TDZ {'alpha0': 0.3, 'n_max': 2}  max|y_TDZ - y_TD(predict=False)| = 2.797e-03
TDZ {'alpha0': 0.3, 'n_max': 4}  max|y_TDZ - y_TD(predict=False)| = 9.216e-05
```

```python
"""Reviewer probe 1: the audit's own complex-hopping chain (seed 3, 4 sites,
A=Cdag_0, B=C_2), exact reference by a hand-built Jordan-Wigner kron, then
(1) TDZ against TD pointwise, (2) the same calls with the FFT stage
replaced by zero-padding or by direct evaluation, (3) the time-domain
error of TDZ's reconstructed C(t) against the exact one, (4) the recorded
remedy (alpha0, n_max), (5) a second delta."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import revhelpers as rh
from dmrgpy import fermionchain

N, SEED, DT = 4, 3, 0.1
ES = np.linspace(-1.0, 6.0, 60)


def build(n=N, seed=SEED):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=3)
    rng = np.random.RandomState(seed)
    t = rng.random((n, n)) + 1j * rng.random((n, n))
    t = t + t.conj().T
    h = 0
    for i in range(n):
        for j in range(n): h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1): h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 30, 12
    # hand-built JW reference, no dmrgpy code
    a = np.array([[0, 1], [0, 0]], dtype=complex)
    Z = np.diag([1.0, -1.0]).astype(complex)
    I2 = np.eye(2, dtype=complex)
    def c(i):
        out = np.array([[1.0 + 0j]])
        for k in range(n):
            out = np.kron(out, Z if k < i else (a if k == i else I2))
        return out
    C = [c(i) for i in range(n)]
    Cd = [x.conj().T for x in C]
    Nn = [Cd[i] @ C[i] for i in range(n)]
    H = sum(t[i, j] * Cd[i] @ C[j] for i in range(n) for j in range(n))
    H = H + sum(0.8 * Nn[i] @ Nn[i + 1] for i in range(n - 1))
    e = np.linalg.eigvalsh(H)
    D, M = rh.lehmann_dense(H, Cd[0], C[2])
    return fc, D, M, e[1] - e[0]


def err(y, ref):
    return np.max(np.abs(np.asarray(y) - ref))


fc, D, M, gap = build()
A, B = fc.Cdag[0], fc.C[2]
print("gap above GS %.4f   max|Im M_n| = %.3f" % (gap, np.max(np.abs(M.imag))))
rh.install()

DELTA = 0.4
ref = rh.density(D, M, ES, DELTA)
nt = int(6 / DELTA / DT); T = nt * DT; dw = 2 * np.pi / T
print("delta=%.2f nt=%d T=%.1f FFT spacing %.4f = %.3f*delta   exact peak %.4f"
      % (DELTA, nt, T, dw, dw / DELTA, np.max(np.abs(ref))))

calls = [("DMRG", "TDZ", {}), ("DMRG", "TD", dict(predict=False)),
         ("ED", "TD", dict(predict=False)), ("DMRG", "TD", {}),
         ("DMRG", "TDZ", dict(predict=True))]
res = {}
for variant in ("orig", "pad", "direct"):
    print("\n-- FFT stage: %s%s" % (variant, " (K=16)" if variant == "pad" else ""))
    for mode, sub, extra in calls:
        rh.use(variant)
        _x, y = fc.get_dynamical_correlator(mode=mode, submode=sub, name=[A, B],
                                            es=ES, delta=DELTA, dt=DT, **extra)
        res[(variant, mode, sub, str(extra))] = np.asarray(y)
        print("  %-4s %-3s %-18s max|y-exact| = %.3e"
              % (mode, sub, str(extra) if extra else "(defaults)", err(y, ref)),
              flush=True)
        if variant == "direct" and (mode, sub, extra) == ("DMRG", "TDZ", {}):
            tails = list(rh.STATE["tail"])
            logs = list(rh.STATE["log"])
        if variant == "direct" and (mode, sub, extra) == ("ED", "TD", dict(predict=False)):
            logs_ed = list(rh.STATE["log"])
        if variant == "direct" and (mode, sub, extra) == ("DMRG", "TD", dict(predict=False)):
            logs_td = list(rh.STATE["log"])

for variant in ("orig", "direct"):
    a = res[(variant, "DMRG", "TDZ", "{}")]
    b = res[(variant, "DMRG", "TD", str(dict(predict=False)))]
    c = res[(variant, "ED", "TD", str(dict(predict=False)))]
    print("\n[%s] pointwise max|y_TDZ - y_TD(predict=False)| = %.3e   "
          "max|y_TDZ - y_ED,TD(predict=False)| = %.3e" % (variant, err(a, b), err(a, c)))

print("\ndamped tail |c(T)|e^{-dT}/max|c e^{-dt}|, and e^{-dT}:",
      ["%.2e, %.2e" % tt for tt in tails])

# time domain: what TDZ hands the FFT, against the exact correlator
print("\ntime domain (pair, then adjoint pair), against sum_n M_n e^{+-i D_n t}:")
for label, lg in (("TDZ", logs), ("TD DMRG", logs_td), ("TD ED", logs_ed)):
    for k, (ts, cs) in enumerate(lg):
        Mk = M if k == 0 else M.conj()
        ep = rh.ctime(D, Mk, ts, +1); em = rh.ctime(D, Mk, ts, -1)
        dc = cs - ep
        bound = np.sum(np.abs(dc) * np.exp(-DELTA * ts)) * DT / np.pi
        print("  %-8s run %d: max|c-exact(+)| = %.2e  max|c-exact(-)| = %.2e  "
              "|c(0)| = %.3f  freq-domain bound (dt/pi)sum|dc|e^{-dt} = %.2e"
              % (label, k, np.max(np.abs(dc)), np.max(np.abs(cs - em)),
                 np.abs(cs[0]), bound))

# on-grid identity: direct evaluation reproduces the FFT exactly at its own grid
EG = dw * np.arange(-2, 15)
refg = rh.density(D, M, EG, DELTA)
out = {}
for variant in ("orig", "direct"):
    rh.use(variant)
    _x, y = fc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=[A, B],
                                        es=EG, delta=DELTA, dt=DT)
    out[variant] = np.asarray(y)
    print("\non the FFT grid, TDZ default [%s]: max|y-exact| = %.3e"
          % (variant, err(y, refg)))
print("on the FFT grid, max|direct - orig| = %.2e" % err(out["direct"], out["orig"]))

# the recorded remedy: smaller alpha0 or larger n_max
print("\nrecorded remedy (user guide: 'reduce it with a smaller alpha0 or a larger n_max'):")
for kw in (dict(alpha0=0.1, n_max=4), dict(alpha0=0.05), dict(alpha0=0.2),
           dict(n_max=2), dict(alpha0=0.3, n_max=4)):
    line = "  TDZ %-28s" % str(kw)
    for variant in ("orig", "direct"):
        rh.use(variant)
        _x, y = fc.get_dynamical_correlator(mode="DMRG", submode="TDZ", name=[A, B],
                                            es=ES, delta=DELTA, dt=DT, **kw)
        line += "  [%s] %.3e" % (variant, err(y, ref))
    print(line, flush=True)

# a second delta on the same chain
DELTA2 = 0.25
ref2 = rh.density(D, M, ES, DELTA2)
nt2 = int(6 / DELTA2 / DT)
print("\ndelta=%.2f  nt=%d  FFT spacing = %.3f*delta  exact peak %.4f"
      % (DELTA2, nt2, 2 * np.pi / (nt2 * DT) / DELTA2, np.max(np.abs(ref2))))
for mode, sub, extra in calls[:3]:
    line = "  %-4s %-3s %-18s" % (mode, sub, str(extra) if extra else "(defaults)")
    for variant in ("orig", "direct"):
        rh.use(variant)
        _x, y = fc.get_dynamical_correlator(mode=mode, submode=sub, name=[A, B],
                                            es=ES, delta=DELTA2, dt=DT, **extra)
        line += "  [%s] %.3e" % (variant, err(y, ref2))
    print(line, flush=True)
```

```
gap above GS 0.1219   max|Im M_n| = 0.271
delta=0.40 nt=150 T=15.0 FFT spacing 0.4189 = 1.047*delta   exact peak 0.2313

-- FFT stage: orig
  DMRG TDZ (defaults)         max|y-exact| = 3.313e-02
  DMRG TD  {'predict': False} max|y-exact| = 3.313e-02
  ED   TD  {'predict': False} max|y-exact| = 3.314e-02
  DMRG TD  (defaults)         max|y-exact| = 2.565e-04
  DMRG TDZ {'predict': True}  max|y-exact| = 2.565e-04

-- FFT stage: pad (K=16)
  DMRG TDZ (defaults)         max|y-exact| = 6.689e-04
  DMRG TD  {'predict': False} max|y-exact| = 6.689e-04
  ED   TD  {'predict': False} max|y-exact| = 6.691e-04
  DMRG TD  (defaults)         max|y-exact| = 3.170e-05
  DMRG TDZ {'predict': True}  max|y-exact| = 3.253e-05

-- FFT stage: direct
  DMRG TDZ (defaults)         max|y-exact| = 5.405e-04
  DMRG TD  {'predict': False} max|y-exact| = 5.395e-04
  ED   TD  {'predict': False} max|y-exact| = 5.393e-04
  DMRG TD  (defaults)         max|y-exact| = 3.145e-05
  DMRG TDZ {'predict': True}  max|y-exact| = 3.225e-05

[orig] pointwise max|y_TDZ - y_TD(predict=False)| = 6.234e-07   max|y_TDZ - y_ED,TD(predict=False)| = 8.555e-07

[direct] pointwise max|y_TDZ - y_TD(predict=False)| = 9.854e-07   max|y_TDZ - y_ED,TD(predict=False)| = 1.174e-06

damped tail |c(T)|e^{-dT}/max|c e^{-dt}|, and e^{-dT}: ['2.58e-03, 2.58e-03', '2.58e-03, 2.58e-03']

time domain (pair, then adjoint pair), against sum_n M_n e^{+-i D_n t}:
  TDZ      run 0: max|c-exact(+)| = 5.01e-06  max|c-exact(-)| = 5.86e-01  |c(0)| = 0.293  freq-domain bound (dt/pi)sum|dc|e^{-dt} = 1.21e-06
  TDZ      run 1: max|c-exact(+)| = 5.41e-06  max|c-exact(-)| = 5.86e-01  |c(0)| = 0.293  freq-domain bound (dt/pi)sum|dc|e^{-dt} = 1.36e-06
  TD DMRG  run 0: max|c-exact(+)| = 3.49e-07  max|c-exact(-)| = 5.86e-01  |c(0)| = 0.293  freq-domain bound (dt/pi)sum|dc|e^{-dt} = 2.82e-07
  TD DMRG  run 1: max|c-exact(+)| = 3.49e-07  max|c-exact(-)| = 5.86e-01  |c(0)| = 0.293  freq-domain bound (dt/pi)sum|dc|e^{-dt} = 2.82e-07
  TD ED    run 0: max|c-exact(+)| = 3.15e-08  max|c-exact(-)| = 5.86e-01  |c(0)| = 0.293  freq-domain bound (dt/pi)sum|dc|e^{-dt} = 4.12e-09
  TD ED    run 1: max|c-exact(+)| = 3.15e-08  max|c-exact(-)| = 5.86e-01  |c(0)| = 0.293  freq-domain bound (dt/pi)sum|dc|e^{-dt} = 4.12e-09

on the FFT grid, TDZ default [orig]: max|y-exact| = 3.453e-04

on the FFT grid, TDZ default [direct]: max|y-exact| = 3.453e-04
on the FFT grid, max|direct - orig| = 5.55e-17

recorded remedy (user guide: 'reduce it with a smaller alpha0 or a larger n_max'):
  TDZ {'alpha0': 0.1, 'n_max': 4}   [orig] 3.313e-02  [direct] 5.405e-04
  TDZ {'alpha0': 0.05}              [orig] 3.313e-02  [direct] 5.395e-04
  TDZ {'alpha0': 0.2}               [orig] 3.313e-02  [direct] 5.621e-04
  TDZ {'n_max': 2}                  [orig] 3.324e-02  [direct] 6.499e-04
  TDZ {'alpha0': 0.3, 'n_max': 4}   [orig] 3.312e-02  [direct] 6.589e-04

delta=0.25  nt=240  FFT spacing = 1.047*delta  exact peak 0.3656
  DMRG TDZ (defaults)          [orig] 6.994e-02  [direct] 6.257e-04
  DMRG TD  {'predict': False}  [orig] 6.994e-02  [direct] 6.224e-04
  ED   TD  {'predict': False}  [orig] 6.994e-02  [direct] 6.221e-04
```

```python
"""Reviewer probe 2: the O1 record's own Hermitian configuration (6-site
S=1/2 Heisenberg, A=B=Sz_0, delta=0.3, es=linspace(0.01,4,30), TDZ at
dt=0.05 nt=400, TD at dt=0.025 nt=800, both T=20), which is where the
recorded "1.8e-2 on a Hermitian pair" comes from. Uniform and staggered
(0.3*(-1)^i Sz_i) chains, exact reference from a numpy kron."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import revhelpers as rh
from dmrgpy import spinchain

n = 6
ES = np.linspace(0.01, 4.0, 30)
DELTA = 0.3
sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2


def site(o, i):
    out = np.array([[1.0 + 0j]])
    for k in range(n):
        out = np.kron(out, o if k == i else np.eye(2))
    return out


rh.install()
for stag in (0.0, 0.3):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=3)
    h = 0
    H = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
        H = H + sum(site(o, i) @ site(o, i + 1) for o in (sx, sy, sz))
    for i in range(n):
        h = h + stag * (-1) ** i * sc.Sz[i]
        H = H + stag * (-1) ** i * site(sz, i)
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 30, 20
    D, M = rh.lehmann_dense(H, site(sz, 0), site(sz, 0))
    ref = rh.density(D, M, ES, DELTA)
    print("\nstaggered field %.1f: exact peak |C_AB| = %.4f  (record: 0.1421)"
          % (stag, np.max(np.abs(ref))))
    for sub, dt, nt, extra in (("TDZ", 0.05, 400, {}), ("TD", 0.05, 400, dict(predict=False)),
                               ("TD", 0.025, 800, dict(predict=False)), ("TD", 0.025, 800, {})):
        line = "  %-3s dt=%.3f nt=%d %-18s" % (sub, dt, nt, str(extra) if extra else "(defaults)")
        for variant in ("orig", "direct"):
            rh.use(variant)
            _x, y = sc.get_dynamical_correlator(mode="DMRG", submode=sub,
                                                name=[sc.Sz[0], sc.Sz[0]], es=ES,
                                                delta=DELTA, dt=dt, nt=nt, **extra)
            y = np.asarray(y)
            line += "  [%s] max|y-exact| = %.3e" % (variant, np.max(np.abs(y - ref)))
        print(line, flush=True)
```

```
staggered field 0.0: exact peak |C_AB| = 0.1421  (record: 0.1421)
  TDZ dt=0.050 nt=400 (defaults)          [orig] max|y-exact| = 1.830e-02  [direct] max|y-exact| = 2.230e-04
  TD  dt=0.050 nt=400 {'predict': False}  [orig] max|y-exact| = 1.830e-02  [direct] max|y-exact| = 2.235e-04
  TD  dt=0.025 nt=800 {'predict': False}  [orig] max|y-exact| = 1.831e-02  [direct] max|y-exact| = 2.200e-04
  TD  dt=0.025 nt=800 (defaults)          [orig] max|y-exact| = 1.526e-04  [direct] max|y-exact| = 3.351e-06

staggered field 0.3: exact peak |C_AB| = 0.1512  (record: 0.1421)
  TDZ dt=0.050 nt=400 (defaults)          [orig] max|y-exact| = 7.691e-03  [direct] max|y-exact| = 3.289e-04
  TD  dt=0.050 nt=400 {'predict': False}  [orig] max|y-exact| = 7.691e-03  [direct] max|y-exact| = 3.289e-04
  TD  dt=0.025 nt=800 {'predict': False}  [orig] max|y-exact| = 7.693e-03  [direct] max|y-exact| = 3.315e-04
  TD  dt=0.025 nt=800 (defaults)          [orig] max|y-exact| = 2.934e-04  [direct] max|y-exact| = 9.154e-06
```

**Suggested fix**: evaluate the damped trapezoid sum directly at each requested
`es` rather than zero-padding by a fixed factor (K=16 still leaves 6.7e-04
against the 5.4e-04 floor, and the K needed scales as 2*pi/(nt*dt*delta), 1.05
here and 6.3 at the infinite-chain defaults), chunked over `es` so the outer
product stays bounded (about 150 MB at defaults otherwise); if padding is kept,
derive K from delta, spacing at most delta/8. Do not fix it by switching TDZ to
`predict=True`, which leaves TD `predict=False` and `sxt_to_skomega` on the
coarse grid and ties TDZ to an extrapolation never evaluated for it
(`docs/td_dynamical_correlator_sharpening_plan.md`). Leave `sxt_to_skomega` on
the current path until its own numbers are measured. Test the contour on its
own, max|y_TDZ - y_TD(predict=False)| < ~3e-06 on the same `es`, and tighten
TDZ-versus-exact to about 1e-3 and the TD/ED rows to about 1e-4. Correct the
sentences listed under Where. NUMBERS CHANGE: TDZ at its defaults and TD at
`predict=False` by up to 14 to 19 per cent of the peak on this chain, towards
exact; TD's default by up to about 1e-3 of the peak.

### 8. `mode="ED"` TD on a degenerate ground state builds the two halves of the O1 density on two different randomly chosen ground states, so it returns the density of no state and no mixture on the manifold

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `td-convention`

**Status**: FIXED -- `edtk/timedependent.evolution_DC(..., wf0=None)` measures in `wf0` (an array or a State), defaulting to the EDchain's cached `get_gs_array()`, and `edtk/dynamics.py`'s TD branch forwards `wf0`; the energy origin is always the cached `self.e0`, as KPM does, so the shift is now `Hop - e0*I` (the sign trap: keeping `+` would have put the result 1.30e-01 off exact on a non-degenerate chain). `evolution_ABC` uses the cached state too when `wf is None`, so `evolution_ABA(mode="ED")` and `evolve_and_measure(mode="ED")` without `wf=` change as well on a degenerate ground state. After the fix ED TD on a degenerate ground state agrees with KPM, INV, CVM and ROOTN and still not with `submode="ED"`'s manifold average, the split between the two ED conventions that predates this fix. Pinned by `tests/test_audit_2026_09_24_realtime.py::test_ed_td_on_a_degenerate_ground_state_uses_the_cached_state`, `::test_ed_td_measures_the_state_it_is_given`, `::test_ed_td_on_a_non_degenerate_ground_state_is_unchanged` and `::test_evolution_aba_on_ed_uses_the_cached_state`. NUMBERS CHANGE only on a degenerate or quasi-degenerate ground state: on the degenerate 3-site Heisenberg chain with (Sp_0, Sm_2), delta=0.3, dt=0.05, es=linspace(-1,3,41), peak 1.768e-01, three identical calls used to spread by 2.03e-01 and to sit 1.42e-01 to 3.42e-01 from the cached-state density, and now are identical and 3.32e-06 from it; on the non-degenerate chain (hx=0.2, hz=0.1 on site 0) the result differs from the pre-fix construction by 2.2e-12 to 2.9e-11 (ARPACK's random start).

**Where**: `src/dmrgpy/edtk/timedependent.py::evolution_DC` (a fresh
`eigsh(-H, k=1)` with a random start on every call; `evolution_ABC` does the
same when `wf is None`); `src/dmrgpy/edtk/dynamics.py`'s TD branch, which does
not forward `wf0=`; `timedependent.lehmann_density_from_one_sided`, which since
`765b537` makes two `evolution_DC` calls for any pair not provably
self-adjoint.

`evolution_DC` ignores both states it has access to, the EDchain's cached
`get_gs_array()` (which `timedependent.dynamical_correlator` has already
computed on this route) and any `wf0=` the caller passed. That part predates
`765b537`: on a degenerate ground state it already made ED TD non-deterministic
and inconsistent with the ED submodes that use the cached state (KPM, INV, CVM,
ROOTN), even for a self-adjoint pair. What `765b537` adds is the second call: a
pair `is_dagger_pair` cannot prove self-adjoint, which includes every
(A_i, B_j) with i != j, now takes one call per half, so on a degenerate ground
state the two halves come from two different members of the manifold, and the
result is the density of the average state plus a dispersive term proportional
to Tr((rho1-rho2) X^n), the same kind of leftover O1 was closed to remove. DMRG
TD reuses the session's `wf0` for both halves and is unaffected; so are
non-degenerate ground states.

**Expected**: both halves on one state, the cached one, so that ED TD agrees
with the other ED submodes that use the cached state.

Repro, from the hunter (it instruments `edtk/timedependent`'s `eigsh`):

```python
"""Lead 2: mode="ED" submode="TD" re-solves the ground state with a
randomly started eigsh inside every evolution_DC call. On a pair that is
not provably self-adjoint the correlator now takes two such calls, so on a
degenerate ground state its two halves come from two different states."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import scipy.sparse.linalg as _slg
from anchor import heis_field
from dmrgpy import spinchain
from dmrgpy.edtk import timedependent as edtd

n = 3
ES = np.linspace(-1.0, 3.0, 41)
DELTA, DT = 0.3, 0.05

recorded = []
class RecordingSlg:
    """scipy.sparse.linalg with eigsh wrapped, for edtk/timedependent only"""
    force = None
    def __getattr__(self, k): return getattr(_slg, k)
    def eigsh(self, *a, **k):
        e, v = _slg.eigsh(*a, **k)
        if RecordingSlg.force is not None: v = RecordingSlg.force.copy()
        recorded.append(v.reshape(-1).copy())
        return e, v
edtd.slg = RecordingSlg()

sc = spinchain.Spin_Chain([2] * n)
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
ed = sc.get_ED_obj()
H = np.array(ed.get_hamiltonian().todense())
e, U = np.linalg.eigh(H)
print("ED spectrum (lowest 4):", np.round(e[:4], 10),
      " hand-built kron spectrum:", np.round(np.linalg.eigvalsh(heis_field(n))[:4], 10))

Aop, Bop = sc.Sz[1], sc.Sx[1]   # Hermitian, but A != B^dagger: two runs
Am = np.array(ed.MO2matrix(Aop).todense())
Bm = np.array(ed.MO2matrix(Bop).todense())

def exact_one_sided(g, X, Y):
    """(1/pi) sum_n <g|X|n><n|Y|g> / (delta + i(w-D_n)), the code's F"""
    Uh = U.conj().T
    M = (g.conj() @ X @ U) * (Uh @ Y @ g)
    D = e - e[0]
    return np.array([np.sum(M / (DELTA + 1j*(w - D))) for w in ES]) / np.pi

def exact_density(g):
    return 0.5*(exact_one_sided(g, Am, Bm)
                + np.conj(exact_one_sided(g, Bm.conj().T, Am.conj().T)))

for trial in range(3):
    del recorded[:]
    RecordingSlg.force = None
    _x, y = sc.get_dynamical_correlator(mode="ED", submode="TD",
                name=[Aop, Bop], es=ES, delta=DELTA, dt=DT)
    y = np.asarray(y)
    g1, g2 = recorded[0], recorded[1]
    mixed = 0.5*(exact_one_sided(g1, Am, Bm)
                 + np.conj(exact_one_sided(g2, Bm.conj().T, Am.conj().T)))
    c1, c2 = exact_density(g1), exact_density(g2)
    print("\ntrial %d: eigsh calls=%d  |<g1|g2>| = %.4f  E(g1)-E0=%.1e E(g2)-E0=%.1e"
          % (trial, len(recorded), abs(np.vdot(g1, g2)),
             (g1.conj()@H@g1).real-e[0], (g2.conj()@H@g2).real-e[0]))
    print("  peak |C_g1| = %.4f   max|C_g1 - C_g2| = %.3e"
          % (np.max(np.abs(c1)), np.max(np.abs(c1 - c2))))
    print("  returned vs exact density of g1: %.3e   of g2: %.3e   "
          "vs half-from-g1 + half-from-g2: %.3e"
          % (np.max(np.abs(y - c1)), np.max(np.abs(y - c2)),
             np.max(np.abs(y - mixed))))
    print("  returned vs density of the equal mixture of g1 and g2: %.3e"
          % np.max(np.abs(y - 0.5*(c1 + c2))))
    # the same call with both runs pinned to one ground state
    del recorded[:]
    RecordingSlg.force = g1
    _x, yp = sc.get_dynamical_correlator(mode="ED", submode="TD",
                name=[Aop, Bop], es=ES, delta=DELTA, dt=DT)
    print("  pinned to g1 for both runs: vs exact density of g1: %.3e"
          % np.max(np.abs(np.asarray(yp) - c1)))
```

```
ED spectrum (lowest 4): [-1. -1. -0. -0.]  hand-built kron spectrum: [-1. -1.  0.  0.]

trial 0: eigsh calls=2  |<g1|g2>| = 0.9348  E(g1)-E0=5.6e-16 E(g2)-E0=-2.2e-16
  peak |C_g1| = 0.0551   max|C_g1 - C_g2| = 4.962e-02
  returned vs exact density of g1: 2.491e-02   of g2: 2.477e-02   vs half-from-g1 + half-from-g2: 6.856e-05
  returned vs density of the equal mixture of g1 and g2: 1.400e-02
  pinned to g1 for both runs: vs exact density of g1: 1.248e-04

trial 1: eigsh calls=2  |<g1|g2>| = 0.3951  E(g1)-E0=-4.4e-16 E(g2)-E0=2.2e-16
  peak |C_g1| = 0.1158   max|C_g1 - C_g2| = 1.964e-01
  returned vs exact density of g1: 9.833e-02   of g2: 9.825e-02   vs half-from-g1 + half-from-g2: 1.293e-04
  returned vs density of the equal mixture of g1 and g2: 5.539e-02
  pinned to g1 for both runs: vs exact density of g1: 2.625e-04

trial 2: eigsh calls=2  |<g1|g2>| = 0.5835  E(g1)-E0=2.2e-16 E(g2)-E0=-4.4e-16
  peak |C_g1| = 0.0413   max|C_g1 - C_g2| = 1.956e-02
  returned vs exact density of g1: 9.861e-03   of g2: 9.718e-03   vs half-from-g1 + half-from-g2: 7.142e-05
  returned vs density of the equal mixture of g1 and g2: 5.517e-03
  pinned to g1 for both runs: vs exact density of g1: 9.359e-05
```

**Reviewer (CONFIRMED, NARROWED)**: the mechanism reproduces (the returned curve
matches "half from g1 plus half from g2" to 4.4e-06 to 5.0e-05 in three trials);
the reviewer's own probe fits each returned curve by least squares over every
Hermitian 2x2 rho on the ground doublet, which is if anything too generous, and
asks whether the curve is the density of any state at all. Three narrowings. The
fresh `eigsh` per call is old; at `fd01679` five identical calls on the
degenerate chain already differed by 1.36e-01, each the one-sided transform of
a single random member with the pair read reversed, so what `765b537` adds is
the mixing of two members. It does not reach every pair, only those whose
density depends on the manifold member: on an SU(2) doublet zz and xx at i != j
are the same for every member and fit to the method's own 1.4e-04 (a property of
that symmetry, not a rule). And the hunter's magnitudes (|<g1|g2>| 0.40 to 0.93,
"24 to 85 per cent of the peak") are a random variable, not sizes for the
record; for the hunter's pair (Sz_1, Sx_1) both the cached-state answer and the
manifold average are exactly zero, so (Sp_0, Sm_2) is the cleaner case, where the
best Hermitian rho misses the curve by 9.6e-04 to 8.0e-02 on a peak of 0.44.
"Every other `mode="ED"` submode uses the cached `get_gs_array()`" narrows to
KPM, INV, CVM and ROOTN, since `submode="ED"` averages the manifold instead. The
reviewer's reach table and the fix check (`probe5c_fix.py` swaps
`evolution_DC` in-process for a copy that reads the cached state, and also
measures the sign trap below):

```python
"""Reviewer probe for candidate 5's reach: which ordinary pairs expose it?
Same degenerate 3-site Heisenberg chain, public mode="ED" submode="TD",
three calls per pair. For each pair: the spread over manifold members of
the exact density (is the answer state-dependent at all?), the spread
between calls, and the best fit of each call over all Hermitian rho on
the doublet (is it the density of ANY state or mixture?)."""
import numpy as np
from dmrgpy import spinchain
from dmrgpy.timedependent import _pair_is_self_adjoint

n, DELTA, DT = 3, 0.3, 0.05
ES = np.linspace(-1.0, 3.0, 41)
sc = spinchain.Spin_Chain([2] * n)
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
ed = sc.get_ED_obj()
H = np.array(ed.get_hamiltonian().todense())
e, U = np.linalg.eigh(H)
D = e - e[0]
L = (DELTA/np.pi)/((ES[:, None] - D[None, :])**2 + DELTA**2)
basis = U[:, :2]
E = [np.array([[1, 0], [0, 0]]), np.array([[0, 0], [0, 1]]),
     np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]])]
Sp = [sc.Sx[i] + 1j*sc.Sy[i] for i in range(n)]
Sm = [sc.Sx[i] - 1j*sc.Sy[i] for i in range(n)]
pairs = {"(Sz_0,Sz_2)": (sc.Sz[0], sc.Sz[2]), "(Sx_0,Sx_2)": (sc.Sx[0], sc.Sx[2]),
         "(Sp_0,Sm_2)": (Sp[0], Sm[2]), "(Sp_1,Sm_1)": (Sp[1], Sm[1]),
         "(Sz_1,Sx_1)": (sc.Sz[1], sc.Sx[1])}
rng = np.random.default_rng(7)
for label, (A, B) in pairs.items():
    Am = np.array(ed.MO2matrix(A).todense()); Bm = np.array(ed.MO2matrix(B).todense())
    a = basis.conj().T @ Am @ U; b = U.conj().T @ Bm @ basis
    dens = lambda rho: L @ np.einsum("ji,in,nj->n", rho, a, b)
    cols = [dens(Ek.astype(complex)) for Ek in E]
    Amat = np.vstack([np.column_stack([c.real for c in cols]),
                      np.column_stack([c.imag for c in cols])])
    # spread of the exact density over 200 random pure manifold members
    cs = []
    for _ in range(200):
        c = rng.normal(size=2) + 1j*rng.normal(size=2); c /= np.linalg.norm(c)
        cs.append(dens(np.outer(c, c.conj())))
    spread = max(np.max(np.abs(x - cs[0])) for x in cs)
    ys = [np.asarray(sc.get_dynamical_correlator(mode="ED", submode="TD", name=[A, B],
          es=ES, delta=DELTA, dt=DT)[1]) for _ in range(3)]
    between = max(np.max(np.abs(ys[i] - ys[j])) for i in range(3) for j in range(i+1, 3))
    res = []
    for y in ys:
        p, *_ = np.linalg.lstsq(Amat, np.concatenate([y.real, y.imag]), rcond=None)
        res.append(np.max(np.abs(y - sum(pk*ck for pk, ck in zip(p, cols)))))
    print("%-12s self-adjoint=%-5s  peak|C| %.3e  spread over members %.2e  "
          "between calls %.2e  best-rho residual %s"
          % (label, _pair_is_self_adjoint([A, B]), max(np.max(np.abs(x)) for x in cs),
             spread, between, " ".join("%.1e" % r for r in res)))
```

```
(Sz_0,Sz_2)  self-adjoint=False  peak|C| 1.129e-01  spread over members 1.11e-16  between calls 1.19e-13  best-rho residual 1.4e-04 1.4e-04 1.4e-04
(Sx_0,Sx_2)  self-adjoint=False  peak|C| 1.129e-01  spread over members 8.33e-17  between calls 1.60e-13  best-rho residual 1.4e-04 1.4e-04 1.4e-04
(Sp_0,Sm_2)  self-adjoint=False  peak|C| 4.434e-01  spread over members 4.16e-01  between calls 1.74e-01  best-rho residual 5.4e-02 9.6e-04 8.0e-02
(Sp_1,Sm_1)  self-adjoint=True   peak|C| 7.044e-01  spread over members 3.13e-01  between calls 3.99e-01  best-rho residual 5.3e-04 1.3e-03 1.3e-03
(Sz_1,Sx_1)  self-adjoint=False  peak|C| 1.165e-01  spread over members 1.26e-01  between calls 4.52e-02  best-rho residual 5.3e-02 1.4e-02 8.9e-03
```

```python
"""Reviewer probe for candidate 5's suggested fix, applied in-process only
(the repo is untouched): edtk/timedependent.evolution_DC re-pointed at the
EDchain's cached ground state. Also the sign trap: the current code's e0
comes from eigsh(-H), i.e. it is -E_gs, while the cached self.e0 is +E_gs,
so a mechanical substitution keeping `Hop + e0*I` evolves with H + E_gs."""
import numpy as np
from scipy.sparse import identity
from dmrgpy import spinchain
from dmrgpy.edtk import timedependent as tded
from dmrgpy.edtk.tdtk import evolve

n, DELTA, DT = 3, 0.3, 0.05
ES = np.linspace(-1.0, 3.0, 41)
ORIG = tded.evolution_DC


def make_patched(sign):
    def evolution_DC(self, h=None, name=None, nt=100, dt=0.01, **kwargs):
        (A, B) = name[1], name[0]
        Hop = self.get_operator(h)
        Aop = self.get_operator(A); Bop = self.get_operator(B)
        ts = np.array([dt*ii for ii in range(nt)])
        wf0 = self.get_gs_array()            # the cached state
        e0 = self.e0                          # +E_gs
        wf = Aop@wf0.copy(); wfc = np.conjugate(wf0)
        ht = Hop + sign*e0*identity(Hop.shape[0], dtype=np.complex128)
        cs = []
        for it in range(nt):
            cs.append(wfc@Bop@wf)
            wf = evolve(wf, ht, t=dt, dt=dt).reshape(-1)
        return ts, np.array(cs)
    return evolution_DC


def chain(hx=0.0, hz0=0.0):
    sc = spinchain.Spin_Chain([2] * n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    for i in range(n):
        h = h + hx*sc.Sx[i]
    sc.set_hamiltonian(h + hz0*sc.Sz[0])
    return sc


def exact_density_of_cached(sc, A, B):
    ed = sc.get_ED_obj()
    H = np.array(ed.get_hamiltonian().todense()); e, U = np.linalg.eigh(H)
    g = ed.get_gs_array()
    Am = np.array(ed.MO2matrix(A).todense()); Bm = np.array(ed.MO2matrix(B).todense())
    M = (g.conj() @ Am @ U) * (U.conj().T @ Bm @ g)
    L = (DELTA/np.pi)/((ES[:, None] - (e - e[0])[None, :])**2 + DELTA**2)
    return L @ M


def td(sc, A, B):
    return np.asarray(sc.get_dynamical_correlator(mode="ED", submode="TD", name=[A, B],
                      es=ES, delta=DELTA, dt=DT)[1])


for label, hx, hz0 in (("degenerate", 0.0, 0.0), ("non-degenerate hx=0.2,hz0=0.1", 0.2, 0.1)):
    print("== %s" % label)
    for pname in ("(Sp_0,Sm_2)", "(Sz_1,Sx_1)", "(Sp_1,Sm_1)"):
        sc = chain(hx, hz0)
        Sp = [sc.Sx[i] + 1j*sc.Sy[i] for i in range(n)]
        Sm = [sc.Sx[i] - 1j*sc.Sy[i] for i in range(n)]
        A, B = {"(Sp_0,Sm_2)": (Sp[0], Sm[2]), "(Sz_1,Sx_1)": (sc.Sz[1], sc.Sx[1]),
                "(Sp_1,Sm_1)": (Sp[1], Sm[1])}[pname]
        ref = exact_density_of_cached(sc, A, B)
        tded.evolution_DC = ORIG
        y_head = [td(sc, A, B) for _ in range(3)]
        tded.evolution_DC = make_patched(-1.0)
        y_fix = [td(sc, A, B) for _ in range(3)]
        tded.evolution_DC = make_patched(+1.0)
        y_trap = td(sc, A, B)
        tded.evolution_DC = ORIG
        y_kpm = np.asarray(sc.get_dynamical_correlator(mode="ED", submode="KPM",
                           name=[A, B], es=ES, delta=DELTA)[1])
        sp = lambda ys: max(np.max(np.abs(ys[i]-ys[j])) for i in range(3) for j in range(i+1, 3))
        print("  %-12s peak|C_cached| %.3e | HEAD: spread %.2e, max|y-C_cached| %.2e | "
              "fixed: spread %.1e, max|y-C_cached| %.2e, max|y-KPM| %.2e | "
              "sign trap: max|y-C_cached| %.2e | fixed vs HEAD call0 %.1e"
              % (pname, np.max(np.abs(ref)), sp(y_head), max(np.max(np.abs(y-ref)) for y in y_head),
                 sp(y_fix), max(np.max(np.abs(y-ref)) for y in y_fix),
                 np.max(np.abs(y_fix[0]-y_kpm)), np.max(np.abs(y_trap-ref)),
                 np.max(np.abs(y_fix[0]-y_head[0]))))
```

```
== degenerate
  (Sp_0,Sm_2)  peak|C_cached| 1.768e-01 | HEAD: spread 1.60e-01, max|y-C_cached| 3.14e-01 | fixed: spread 0.0e+00, max|y-C_cached| 3.96e-04, max|y-KPM| 3.18e-01 | sign trap: max|y-C_cached| 1.73e-01 | fixed vs HEAD call0 3.0e-01
  (Sz_1,Sx_1)  peak|C_cached| 0.000e+00 | HEAD: spread 6.74e-02, max|y-C_cached| 6.53e-02 | fixed: spread 0.0e+00, max|y-C_cached| 0.00e+00, max|y-KPM| 0.00e+00 | sign trap: max|y-C_cached| 0.00e+00 | fixed vs HEAD call0 5.0e-02
  (Sp_1,Sm_1)  peak|C_cached| 7.074e-01 | HEAD: spread 1.10e-01, max|y-C_cached| 1.27e-01 | fixed: spread 0.0e+00, max|y-C_cached| 1.58e-03, max|y-KPM| 1.27e+00 | sign trap: max|y-C_cached| 6.92e-01 | fixed vs HEAD call0 2.1e-02
== non-degenerate hx=0.2,hz0=0.1
  (Sp_0,Sm_2)  peak|C_cached| 1.403e-01 | HEAD: spread 4.54e-12, max|y-C_cached| 1.74e-04 | fixed: spread 0.0e+00, max|y-C_cached| 1.74e-04, max|y-KPM| 1.72e-01 | sign trap: max|y-C_cached| 1.30e-01 | fixed vs HEAD call0 7.4e-12
  (Sz_1,Sx_1)  peak|C_cached| 1.434e-02 | HEAD: spread 1.74e-14, max|y-C_cached| 8.20e-05 | fixed: spread 0.0e+00, max|y-C_cached| 8.20e-05, max|y-KPM| 1.90e-02 | sign trap: max|y-C_cached| 1.47e-02 | fixed vs HEAD call0 1.4e-14
  (Sp_1,Sm_1)  peak|C_cached| 4.829e-01 | HEAD: spread 6.18e-13, max|y-C_cached| 5.01e-04 | fixed: spread 0.0e+00, max|y-C_cached| 5.01e-04, max|y-KPM| 3.50e-01 | sign trap: max|y-C_cached| 4.72e-01 | fixed vs HEAD call0 5.6e-13
```

**Suggested fix**: give `evolution_DC` a `wf0=None` parameter, the way
`evolution_ABC` has `wf`, defaulting to `self.get_gs_array()` and `self.e0`, and
forward `wf0` from `edtk/dynamics.py`'s TD branch, which also fixes the dropped
`wf0=`. There is a sign trap a test on a self-adjoint pair would not catch: the
current `e0` comes from `eigsh(-Hop)`, so `e0[0]` is -E_gs and the shift line
reads `Hop + e0[0]*I`; the cached `self.e0` is +E_gs, so the line has to become
`Hop - self.e0*I`, and keeping the `+` evolves with H + E_gs, an error of 1.3e-01
to 6.9e-01 (the `sign trap` column above). After the fix ED TD on a degenerate
ground state agrees with KPM, INV, CVM and ROOTN but still not with
`submode="ED"`'s manifold average, the pre-existing split between the two ED
conventions. Numbers change only on a degenerate or quasi-degenerate ground
state, by at most 7.4e-12 on a non-degenerate one.

### 9. Since `765b537`, `submode="TDZ"` accepts any misspelled or unknown keyword and returns the default-parameter spectrum bit for bit, and its check that the operators are symbolic is gone

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `td-convention`

**Status**: FIXED -- `tdz.dynamical_correlator_tdz` has no `**kwargs` any more and takes `i=0, j=0` as named parameters, and calls `str2MO(self, name, i=i, j=j, require_symbolic_for="submode='TDZ'")` before anything else; `alpha=`, `nmax=` and `foo=` now raise `TypeError` (Python itself suggests `alpha0`/`n_max`), and a `toMPO()` operator raises the `TypeError` naming `toMPO`. Every call site in `tests/`, `examples/`, `benchmarks/` and `mpsjulialive/dynamics.py` passes named parameters only; on `julia_live` a typo now raises through the forwarded `**kwargs` (by reading). Pinned by `tests/test_audit_2026_09_24_realtime.py::test_tdz_rejects_an_unknown_keyword` and `::test_tdz_names_a_compiled_operator`, and by `tests/test_tompo_correlator_operators.py::test_symbolic_only_submodes_name_the_problem[TDZ]`. No correctly spelled call changes. A sibling is left open and recorded as a lead: `mode="ED"` `submode="TD"` still accepts unknown keywords through `edtk/timedependent.evolution_DC`'s own `**kwargs`, which removing would make DMRG-only keywords such as `restart=` raise on ED.

**Where**: `src/dmrgpy/tdz.py::dynamical_correlator_tdz(..., **kwargs)`, whose
body no longer reads `kwargs`; `timedependent.lehmann_density_from_one_sided`,
which resolves `name` with `str2MO(self, name)` and no keywords;
`mpsjulialive/dynamics.py:106`, which forwards `**kwargs` into the same function
(by reading).

Before `765b537` the body forwarded `**kwargs` to `operatornames.str2MO`, which
raised `TypeError` for anything other than `i`, `j` and `require_symbolic_for`.
The rewrite moved name resolution into `lehmann_density_from_one_sided`, so on
the public route `get_dynamical_correlator(mode="DMRG", submode="TDZ", ...)`
`alpha` for `alpha0`, `nmax` for `n_max`, or `foo` are accepted and ignored,
while `submode="TD"` still raises (its keywords reach `evolution_dmrg_DC`'s
`str2MO`). `docs/documentation.md` section 4.10 states the policy this breaks:
"Solver parameters are chain attributes, never call arguments, so anything
unrecognized now raises `TypeError`". What is lost is exactly what the typo was
meant to change: on the hunter's chain a real `n_max=0` gives 3.58e-02 against
exact where the default gives 1.61e-02. The same rewrite dropped
`require_symbolic_for="submode='TDZ'"`, so a `toMPO()` operator now crashes
several frames deep instead of raising the `TypeError` naming `toMPO` that the
symbolic-only submodes are meant to raise; `test_symbolic_only_submodes_name_the_problem`
parametrizes only KPM and TD, which is why nothing caught it.

**Expected**: `TypeError` naming the unknown keyword, as TD raises and as TDZ
raised at `1c2606d`.

Repro, from the hunter, at HEAD and on a `git archive` of `1c2606d`
(Python-only):

```python
"""Candidate A: since 765b537 dynamical_correlator_tdz no longer passes
**kwargs anywhere, so a misspelled TDZ keyword is silently dropped."""
import numpy as np
import dmrgpy
from dmrgpy import spinchain
print("dmrgpy from", dmrgpy.__file__)
n = 4
sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm, sc.nsweeps = 20, 8
ES = np.linspace(0.0, 3.0, 13)
kw = dict(mode="DMRG", name=[sc.Sz[0], sc.Sz[0]], es=ES, delta=0.4, dt=0.1)
for sub in ("TDZ", "TD"):
    try:
        _x, y_typo = sc.get_dynamical_correlator(submode=sub, alpha=0.3, **kw)
        _x, y_def = sc.get_dynamical_correlator(submode=sub, **kw)
        print("%-3s alpha=0.3 (typo for alpha0) accepted; max|y_typo - y_default| = %.1e"
              % (sub, np.max(np.abs(np.asarray(y_typo) - np.asarray(y_def)))))
    except TypeError as ex:
        print("%-3s alpha=0.3 raised TypeError: %s" % (sub, ex))
_x, y_real = sc.get_dynamical_correlator(submode="TDZ", alpha0=0.3, **kw)
_x, y_def = sc.get_dynamical_correlator(submode="TDZ", **kw)
print("TDZ alpha0=0.3 (the real keyword) vs default: max diff = %.1e"
      % np.max(np.abs(np.asarray(y_real) - np.asarray(y_def))))
# a knob that moves the answer more: the Taylor order, misspelled
from anchor import heis_field, site_op, sz, lehmann_dense, density
D, M, _e = lehmann_dense(heis_field(n), site_op(sz, 0, n), site_op(sz, 0, n))
ref = density(D, M, ES, 0.4)
try:
    _x, y_typo = sc.get_dynamical_correlator(submode="TDZ", nmax=0, alpha0=0.3, **kw)
    _x, y_true = sc.get_dynamical_correlator(submode="TDZ", n_max=0, alpha0=0.3, **kw)
    _x, y_d = sc.get_dynamical_correlator(submode="TDZ", alpha0=0.3, **kw)
    print("TDZ nmax=0 (typo for n_max): max|y - exact| = %.3e, n_max=0 really: %.3e,"
          " default n_max=4: %.3e, typo == default to %.1e"
          % (np.max(np.abs(np.asarray(y_typo)-ref)), np.max(np.abs(np.asarray(y_true)-ref)),
             np.max(np.abs(np.asarray(y_d)-ref)), np.max(np.abs(np.asarray(y_typo)-np.asarray(y_d)))))
except TypeError as ex:
    print("TDZ nmax=0 raised TypeError: %s" % ex)
```

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
TDZ alpha=0.3 (typo for alpha0) accepted; max|y_typo - y_default| = 0.0e+00
TD  alpha=0.3 raised TypeError: str2MO() got an unexpected keyword argument 'alpha'
TDZ alpha0=0.3 (the real keyword) vs default: max diff = 4.9e-05
TDZ nmax=0 (typo for n_max): max|y - exact| = 1.605e-02, n_max=0 really: 3.581e-02, default n_max=4: 1.605e-02, typo == default to 0.0e+00
```

```
dmrgpy from <scratch>/td-convention/old_1c2606d/src/dmrgpy/__init__.py
TDZ alpha=0.3 raised TypeError: str2MO() got an unexpected keyword argument 'alpha'
TD  alpha=0.3 raised TypeError: str2MO() got an unexpected keyword argument 'alpha'
TDZ alpha0=0.3 (the real keyword) vs default: max diff = 5.4e-05
TDZ nmax=0 raised TypeError: str2MO() got an unexpected keyword argument 'nmax'
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced through the public route only, on
`"python"` and v3 at HEAD and on `"python"` at `fd01679`, the direct parent of
`765b537`. No validation happens upstream: `Many_Body_Chain.get_dynamical_correlator`
pops `mode`, `name`, `i` and `j`, `dynamics.get_dynamical_correlator` pops
`submode`, and everything else reaches `dynamical_correlator_tdz` untouched.
Struck: "a compiled operator now dies on `.get_dagger()`";
`StaticOperator.get_dagger()` exists and succeeds, and the crash comes one call
later, inside `toMPO`, as `AttributeError: 'StaticOperator' object has no
attribute 'to_terms'`. The conclusion stands, a crash rather than a silent wrong
answer, with the clear message gone. The reviewer's probe, which also carries
the sibling measurement recorded as finding 10:

```python
"""Reviewer probe for candidate 4: does submode="TDZ" swallow unknown
keywords on the public route, on more than one backend, and what does a
compiled operator do there?  Usage: python3 probe4_tdz_kwargs.py <backend>
with <backend> in {python, 3}."""
import sys, traceback
import numpy as np
import dmrgpy
from dmrgpy import spinchain

backend = sys.argv[1] if len(sys.argv) > 1 else "python"
backend = 3 if backend == "3" else backend
print("dmrgpy from", dmrgpy.__file__, " backend", repr(backend))

n = 4
sc = spinchain.Spin_Chain([2] * n, itensor_version=backend)
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm, sc.nsweeps = 20, 8          # pinned, not inherited
ES = np.linspace(0.0, 3.0, 13)
base = dict(mode="DMRG", es=ES, delta=0.4, dt=0.1)


def run(sub, name, **kw):
    try:
        return np.asarray(sc.get_dynamical_correlator(submode=sub, name=name,
                                                     **base, **kw)[1]), None
    except Exception as ex:
        return None, "%s: %s" % (type(ex).__name__, str(ex).splitlines()[0][:110])


P0 = [sc.Sz[0], sc.Sz[0]]
y_def, _ = run("TDZ", P0)
for label, kw in [("foo=1 (no such keyword)", dict(foo=1)),
                  ("alpha=0.3 (typo, alpha0)", dict(alpha=0.3)),
                  ("nmax=0 (typo, n_max)", dict(nmax=0)),
                  ("tevol='MPO' (typo, attribute)", dict(tevol="MPO"))]:
    y, err = run("TDZ", P0, **kw)
    if err: print("TDZ %-30s -> %s" % (label, err))
    else: print("TDZ %-30s -> accepted, max|y - y_default| = %.1e"
                % (label, np.max(np.abs(y - y_def))))
y_n0, _ = run("TDZ", P0, n_max=0)
y_a3, _ = run("TDZ", P0, alpha0=0.3, n_max=1)
print("TDZ real n_max=0 vs default: %.3e ; real alpha0=0.3,n_max=1 vs default: %.3e"
      % (np.max(np.abs(y_n0 - y_def)), np.max(np.abs(y_a3 - y_def))))
for label, kw in [("foo=1", dict(foo=1)), ("alpha=0.3", dict(alpha=0.3))]:
    y, err = run("TD", P0, **kw)
    print("TD  %-30s -> %s" % (label, err if err else "accepted"))

# compiled operator through the public TDZ route: which exception, where
M0 = sc.toMPO(sc.Sz[0])
try:
    sc.get_dynamical_correlator(submode="TDZ", name=[M0, M0], **base)
    print("TDZ compiled operator: accepted")
except Exception as ex:
    tb = traceback.extract_tb(ex.__traceback__)
    print("TDZ compiled operator -> %s: %s" % (type(ex).__name__, str(ex).splitlines()[0][:120]))
    for fr in tb[-3:]:
        print("     at %s:%d in %s: %s" % (fr.filename.split("dmrgpy/")[-1], fr.lineno, fr.name, fr.line))
try:
    sc.get_dynamical_correlator(submode="TD", name=[M0, M0], **base)
    print("TD  compiled operator: accepted")
except Exception as ex:
    print("TD  compiled operator -> %s: %s" % (type(ex).__name__, str(ex).splitlines()[0][:120]))

# i=/j= on the lower-level get_dynamical_correlator_MB route (string name)
for sub in ("TDZ", "TD"):
    try:
        y_str = np.asarray(sc.get_dynamical_correlator_MB(submode=sub, name="ZZ",
                           i=1, j=1, **{k: v for k, v in base.items() if k != "mode"})[1])
        y11, _ = run(sub, [sc.Sz[1], sc.Sz[1]])
        y00, _ = run(sub, [sc.Sz[0], sc.Sz[0]])
        print("%-3s _MB(name='ZZ',i=1,j=1): max|y - C[Sz1,Sz1]| = %.3e, max|y - C[Sz0,Sz0]| = %.3e"
              % (sub, np.max(np.abs(y_str - y11)), np.max(np.abs(y_str - y00))))
    except Exception as ex:
        print("%-3s _MB(name='ZZ',i=1,j=1) -> %s: %s" % (sub, type(ex).__name__, str(ex)[:100]))
```

```
dmrgpy from <repo>/src/dmrgpy/__init__.py  backend 'python'
TDZ foo=1 (no such keyword)        -> accepted, max|y - y_default| = 0.0e+00
TDZ alpha=0.3 (typo, alpha0)       -> accepted, max|y - y_default| = 0.0e+00
TDZ nmax=0 (typo, n_max)           -> accepted, max|y - y_default| = 0.0e+00
TDZ tevol='MPO' (typo, attribute)  -> accepted, max|y - y_default| = 0.0e+00
TDZ real n_max=0 vs default: 8.729e-03 ; real alpha0=0.3,n_max=1 vs default: 4.531e-03
TD  foo=1                          -> TypeError: str2MO() got an unexpected keyword argument 'foo'
TD  alpha=0.3                      -> TypeError: str2MO() got an unexpected keyword argument 'alpha'
TDZ compiled operator -> AttributeError: 'StaticOperator' object has no attribute 'to_terms'
     at manybodychain.py:832 in toMPO: return mpsalgebra.toMPO(self,H,**kwargs)
     at mpsalgebra.py:399 in toMPO: return StaticOperator(H,self)
     at multioperatortk/staticoperator.py:13 in __init__: self.cpp_handle = MBO._session.build_operator(MO.to_terms())
TD  compiled operator -> TypeError: submode='TD' rebuilds its operators inside the backend from their symbolic term list, so it needs the MultiOperator itse
TDZ _MB(name='ZZ',i=1,j=1): max|y - C[Sz1,Sz1]| = 3.737e-02, max|y - C[Sz0,Sz0]| = 0.000e+00
TD  _MB(name='ZZ',i=1,j=1): max|y - C[Sz1,Sz1]| = 4.146e-02, max|y - C[Sz0,Sz0]| = 0.000e+00
```

```
dmrgpy from <repo>/src/dmrgpy/__init__.py  backend 3
TDZ foo=1 (no such keyword)        -> accepted, max|y - y_default| = 0.0e+00
TDZ alpha=0.3 (typo, alpha0)       -> accepted, max|y - y_default| = 0.0e+00
TDZ nmax=0 (typo, n_max)           -> accepted, max|y - y_default| = 0.0e+00
TDZ tevol='MPO' (typo, attribute)  -> accepted, max|y - y_default| = 0.0e+00
TDZ real n_max=0 vs default: 8.729e-03 ; real alpha0=0.3,n_max=1 vs default: 4.531e-03
TD  foo=1                          -> TypeError: str2MO() got an unexpected keyword argument 'foo'
TD  alpha=0.3                      -> TypeError: str2MO() got an unexpected keyword argument 'alpha'
TDZ compiled operator -> AttributeError: 'StaticOperator' object has no attribute 'to_terms'
     at manybodychain.py:832 in toMPO: return mpsalgebra.toMPO(self,H,**kwargs)
     at mpsalgebra.py:399 in toMPO: return StaticOperator(H,self)
     at multioperatortk/staticoperator.py:13 in __init__: self.cpp_handle = MBO._session.build_operator(MO.to_terms())
TD  compiled operator -> TypeError: submode='TD' rebuilds its operators inside the backend from their symbolic term list, so it needs the MultiOperator itse
TDZ _MB(name='ZZ',i=1,j=1): max|y - C[Sz1,Sz1]| = 3.737e-02, max|y - C[Sz0,Sz0]| = 0.000e+00
TD  _MB(name='ZZ',i=1,j=1): max|y - C[Sz1,Sz1]| = 4.146e-02, max|y - C[Sz0,Sz0]| = 0.000e+00
```

```
dmrgpy from <scratch>/review-td-45/old_fd01679/src/dmrgpy/__init__.py  backend 'python'
TDZ foo=1 (no such keyword)        -> TypeError: str2MO() got an unexpected keyword argument 'foo'
TDZ alpha=0.3 (typo, alpha0)       -> TypeError: str2MO() got an unexpected keyword argument 'alpha'
TDZ nmax=0 (typo, n_max)           -> TypeError: str2MO() got an unexpected keyword argument 'nmax'
TDZ tevol='MPO' (typo, attribute)  -> TypeError: str2MO() got an unexpected keyword argument 'tevol'
TDZ real n_max=0 vs default: 8.743e-03 ; real alpha0=0.3,n_max=1 vs default: 4.942e-03
TD  foo=1                          -> TypeError: str2MO() got an unexpected keyword argument 'foo'
TD  alpha=0.3                      -> TypeError: str2MO() got an unexpected keyword argument 'alpha'
TDZ compiled operator -> TypeError: submode='TDZ' rebuilds its operators inside the backend from their symbolic term list, so it needs the MultiOperator its
     at tdz.py:300 in dynamical_correlator_tdz: name = operatornames.str2MO(self, name,
     at operatornames.py:181 in str2MO: require_symbolic(A,require_symbolic_for)
     at operatornames.py:118 in require_symbolic: raise TypeError(
TD  compiled operator -> TypeError: submode='TD' rebuilds its operators inside the backend from their symbolic term list, so it needs the MultiOperator itse
TDZ _MB(name='ZZ',i=1,j=1): max|y - C[Sz1,Sz1]| = 0.000e+00, max|y - C[Sz0,Sz0]| = 4.183e-02
TD  _MB(name='ZZ',i=1,j=1): max|y - C[Sz1,Sz1]| = 0.000e+00, max|y - C[Sz0,Sz0]| = 4.668e-02
```

**Suggested fix**: delete `**kwargs` from `dynamical_correlator_tdz`'s
signature and let Python name the stray key (every call site in `tests/`,
`examples/` and `mpsjulialive/dynamics.py` passes only named parameters), or
raise `TypeError` on a non-empty `kwargs`; restore
`operatornames.str2MO(self, name, require_symbolic_for="submode='TDZ'")` in
`tdz.py` before `lehmann_density_from_one_sided` is called; add "TDZ" to
`test_symbolic_only_submodes_name_the_problem` plus a bad-keyword test. Mind
finding 10 when doing it: deleting `**kwargs` turns TDZ's `i=`/`j=` on the
lower-level route into a `TypeError`, better than a silently wrong site but short
of `fd01679`, which honoured them; forwarding `i`/`j` through
`lehmann_density_from_one_sided` closes both. No returned number changes for a
correctly spelled call.

### 10. Since `765b537`, the lower-level correlator route drops `i=`/`j=` for TD and TDZ: `get_dynamical_correlator_MB(name="ZZ", i=1, j=1)` returns C[Sz_0,Sz_0] bit for bit

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `td-convention` (found by the reviewer of finding 9)

**Status**: FIXED for TD and TDZ -- `lehmann_density_from_one_sided(self, name, transform, i=0, j=0)` resolves the pair with `i`/`j`; `timedependent.dynamical_correlator` pops only `i` and `j` and passes them on, so every other key still reaches `str2MO` and TD keeps its `TypeError` on typos, and TDZ forwards its named `i`/`j`. On a 4-site Heisenberg chain (v3) `get_dynamical_correlator_MB(name="ZZ", i=1, j=1)` was 4.146e-02 (TD) and 3.737e-02 (TDZ) from C[Sz_1,Sz_1] and now matches it to 0.0, and (0,2) matches too; the public route is unchanged. ROOTN's sibling (`rootndmrg.py:54`) still drops `i`/`j` on this route and stays open. Pinned by `tests/test_audit_2026_09_24_realtime.py::test_lower_level_route_honours_the_sites` (TD and TDZ, sites (1,1) and (0,2)). No number changes on the public route.

**Where**: `src/dmrgpy/timedependent.py::lehmann_density_from_one_sided` (line
371, a bare `str2MO(self, name)`), which since `765b537` resolves the operator
names for TD and TDZ; `tdz.py:317`; reached through
`Many_Body_Chain.get_dynamical_correlator_MB` or `dynamics.get_dynamical_correlator`
called directly with a string `name=` and `i=`/`j=`. The pre-existing sibling on
the same route: `rootndmrg.py:54`, a bare `str2MO(self, name)` written when the
2026-08 record's finding #12 replaced its bare `raise`, whose `**kwargs` swallows
`i`/`j`.

At `fd01679` both `tdz.dynamical_correlator_tdz` and `evolution_dmrg_DC` called
`str2MO(..., **kwargs)` with `i`/`j` inside `kwargs`, so the lower-level route
honoured the site keywords for TD and TDZ. After `765b537` the pair is resolved
with i=j=0 before either function sees it; TD's `i`/`j` still reach
`evolution_dmrg_DC`'s `str2MO`, but by then the pair is already resolved and
`str2MO` is idempotent on a pair, so they are ignored. The public
`Many_Body_Chain.get_dynamical_correlator` resolves `name`/`i`/`j` itself and is
unaffected on every backend. Nothing that ships reaches the affected
combination: `get_dynamical_correlator_MB` appears in no test, example,
notebook, README or user-guide section, the two other functions called
`get_dynamical_correlator` that call a `dynamics` module import `edtk.dynamics`,
the `fermionchaintk/dynamicalcorrelator.py` wrappers that would pass a string
with `i`/`j` are imported nowhere, and the last real caller of `_MB`, the
parafermion override, was deleted in the 2026-09 audit.

**Expected**: C[Sz_1,Sz_1] for `name="ZZ", i=1, j=1` on every submode of the
route, as at `fd01679`.

Repro: the finder's measurement is the last two lines of
`probe4_head_python.out` under finding 9 (TDZ 3.737e-02 and TD 4.146e-02 from
C[Sz_1,Sz_1], exactly 0 from C[Sz_0,Sz_0]).

**Reviewer (CONFIRMED, NARROWED)**: reproduced on `"python"` and on 3 at HEAD
(digit for digit identical) and on `"python"` at `fd01679`, with every submode of
the route compared against the explicit pairs on the same chain. Narrowed in
framing only, nothing struck: TD and TDZ are not the only submodes that drop
`i`/`j` there, since ROOTN has done so since the 2026-08 fix (unchanged between
the two trees), while KPM, CVM, EX and SECTOR honour them and CVM_explicit
refuses them with a `TypeError`. The route's own design intent weighs against it
(the public method's comment: "Resolve name= here, once, for every submode and
both solvers"), which sets the severity, but four submodes still honour `i`/`j`
there and TD/TDZ provably did before `765b537`. The absolute |C11-C00| of KPM, TD
and TDZ changes between the trees because of O1/O2, which `CLAUDE.md` already
records as not comparable. The reviewer's probe (the per-frequency
`CVM in E = ... CG iterations` progress lines and the old tree's
missing-extension `UserWarning` are dropped):

```python
"""Refutation probe: on the lower-level route (get_dynamical_correlator_MB,
i.e. dynamics.get_dynamical_correlator called directly) with a string
name="ZZ" and i=1,j=1, which submodes honour the site keywords?

For each submode, compare the string call against the explicit pairs
[Sz1,Sz1] and [Sz0,Sz0] on the same chain (ground state cached, so a
matching pair is bit-identical). Also check the documented public route
get_dynamical_correlator(name="ZZ",i=1,j=1).
Usage: python3 probe_ij_all_submodes.py <backend> [submodes...]"""
import sys
import numpy as np
import dmrgpy
from dmrgpy import spinchain

backend = sys.argv[1] if len(sys.argv) > 1 else "python"
backend = 3 if backend == "3" else backend
subs = sys.argv[2:] or ["KPM", "CVM", "EX", "ROOTN", "TD", "TDZ",
                        "CVM_explicit", "CVMimag"]
print("dmrgpy from", dmrgpy.__file__, " backend", repr(backend))

n = 4
sc = spinchain.Spin_Chain([2] * n, itensor_version=backend)
h = 0
for i in range(n - 1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm, sc.nsweeps = 20, 10          # pinned
ES = np.linspace(0.0, 3.0, 10)
args = {"KPM": dict(es=ES, delta=0.4),
        "CVM": dict(es=ES, delta=0.4),
        "EX": dict(es=ES, delta=0.4, nex=8),
        "ROOTN": dict(es=ES, delta=0.4, N=4, nkry=10),
        "TD": dict(es=ES, delta=0.4, dt=0.1),
        "TDZ": dict(es=ES, delta=0.4, dt=0.1),
        "CVM_explicit": dict(es=ES, delta=0.4),
        "CVMimag": dict(es=ES, delta=0.4)}
P1 = [sc.Sz[1], sc.Sz[1]]
P0 = [sc.Sz[0], sc.Sz[0]]


def call(f, **kw):
    try:
        return np.asarray(f(**kw)[1]), None
    except Exception as ex:
        return None, "%s: %s" % (type(ex).__name__, str(ex).splitlines()[0][:90])


for sub in subs:
    a = args[sub]
    y11, e11 = call(sc.get_dynamical_correlator, submode=sub, name=P1, **a)
    y00, e00 = call(sc.get_dynamical_correlator, submode=sub, name=P0, **a)
    ymb, emb = call(sc.get_dynamical_correlator_MB, submode=sub, name="ZZ",
                    i=1, j=1, **a)
    ypub, epub = call(sc.get_dynamical_correlator, submode=sub, name="ZZ",
                      i=1, j=1, **a)
    if e11 or e00:
        print("%-12s explicit pair itself fails: %s" % (sub, e11 or e00))
        continue
    sep = np.max(np.abs(y11 - y00))
    if emb:
        mb = "_MB -> " + emb
    else:
        mb = "_MB: |y-C11|=%.3e |y-C00|=%.3e" % (np.max(np.abs(ymb - y11)),
                                                  np.max(np.abs(ymb - y00)))
    pub = ("public -> " + epub) if epub else \
        "public: |y-C11|=%.3e" % np.max(np.abs(ypub - y11))
    print("%-12s |C11-C00|=%.3e  %s  %s" % (sub, sep, mb, pub))
    sys.stdout.flush()
```

```
dmrgpy from <repo>/src/dmrgpy/__init__.py  backend 'python'
KPM          |C11-C00|=9.536e-02  _MB: |y-C11|=0.000e+00 |y-C00|=9.536e-02  public: |y-C11|=0.000e+00
CVM          |C11-C00|=4.336e-02  _MB: |y-C11|=0.000e+00 |y-C00|=4.336e-02  public: |y-C11|=0.000e+00
EX           |C11-C00|=4.687e-02  _MB: |y-C11|=0.000e+00 |y-C00|=4.687e-02  public: |y-C11|=0.000e+00
ROOTN        |C11-C00|=4.336e-02  _MB: |y-C11|=4.336e-02 |y-C00|=0.000e+00  public: |y-C11|=0.000e+00
TD           |C11-C00|=4.326e-02  _MB: |y-C11|=4.326e-02 |y-C00|=0.000e+00  public: |y-C11|=0.000e+00
TDZ          |C11-C00|=3.737e-02  _MB: |y-C11|=3.737e-02 |y-C00|=0.000e+00  public: |y-C11|=0.000e+00
CVM_explicit |C11-C00|=4.336e-02  _MB -> TypeError: dynamical_correlator_cvm_explicit() got an unexpected keyword argument 'i'  public: |y-C11|=0.000e+00
CVMimag      explicit pair itself fails: NotImplementedError: submode='CVMimag' needs the Pade analytic-continuation module dmrgpy.padetk, which is not
```

```
dmrgpy from <scratch>/review-td-45/old_fd01679/src/dmrgpy/__init__.py  backend 'python'
KPM          |C11-C00|=1.043e-01  _MB: |y-C11|=0.000e+00 |y-C00|=1.043e-01  public: |y-C11|=0.000e+00
EX           |C11-C00|=4.687e-02  _MB: |y-C11|=0.000e+00 |y-C00|=4.687e-02  public: |y-C11|=0.000e+00
ROOTN        |C11-C00|=4.336e-02  _MB: |y-C11|=4.336e-02 |y-C00|=0.000e+00  public: |y-C11|=0.000e+00
TD           |C11-C00|=4.650e-02  _MB: |y-C11|=0.000e+00 |y-C00|=4.650e-02  public: |y-C11|=0.000e+00
TDZ          |C11-C00|=4.280e-02  _MB: |y-C11|=0.000e+00 |y-C00|=4.280e-02  public: |y-C11|=0.000e+00
CVM_explicit |C11-C00|=4.336e-02  _MB -> TypeError: dynamical_correlator_cvm_explicit() got an unexpected keyword argument 'i'  public: |y-C11|=0.000e+00
```

**Suggested fix**: forward `i`/`j` through `lehmann_density_from_one_sided`,
with `timedependent.dynamical_correlator` popping them from `kwargs` rather than
letting them ride on into `evolution_DC`, and `tdz.dynamical_correlator_tdz`
taking `i=0, j=0` as named parameters (which fits finding 9's removal of its
`**kwargs`). The broader alternative is the shape the 2026-08 record's #12 note
prescribed: resolve `name`/`i`/`j` once in `dynamics.get_dynamical_correlator`
before submode dispatch, which would also make ROOTN honour them and let
CVM_explicit accept a string name, provided `i`/`j` are popped before dispatch
(CVM_explicit's and CVMimag's signatures take neither). No number changes on the
public route.

### 11. `get_kondo_spectrum(mode="ED")` accepts `**kwargs` and reads none of them, so a misspelled physical parameter silently falls back to its default: `Jrho=` for `Jrho_s=` removes the whole third-order Kondo peak (34 per cent of it on the paper's Fig. 7b setup)

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `recent-misc`

**Status**: FIXED -- `get_kondo_spectrum(mode="ED")` raises `TypeError` on any keyword outside its named parameters, naming them sorted and saying they are `mode="DMRG"` parameters or misspellings (no allow-list); the docstring no longer offers the dead `n=` and names what the KPM correlator consumes (`kernel`, `hodc_order`, `hodc_eta`), with the moment count set through the chain's `kpm_*` attributes, and the same stale `n` comment in `kondo_third_order_dmrg_VS_ED` is corrected. No in-tree `mode="ED"` caller (examples, tests, the user-guide snippet, notebooks) passes an extra keyword. Pinned by `tests/test_audit_2026_09_24_kondo.py::test_ed_mode_rejects_keywords_it_does_not_read` (`Jrho=`, `delta=`, `submode=`, `n_gs=`), `::test_ed_mode_names_every_unknown_keyword_sorted` and `::test_ed_mode_correct_spelling_is_unchanged` (1.1374). No number changes for a correctly spelled call.

**Where**: `src/dmrgpy/spinchain.py:148` to `150` (`get_kondo_spectrum(eV, ...,
mode="ED", **kwargs)`) and `:244` to `246`, where the default `mode="ED"` calls
`_get_kondo_spectrum_ed(eV, site, Jrho_s, U, T, T0, omega0, Gamma0, order, kB)`
and passes nothing from `kwargs`; the companion docstring defect at
`spinchain.py:215` to `218` and the matching comment in
`examples/kondo/kondo_third_order_dmrg_VS_ED`.

Every keyword that is not one of the eleven named parameters is accepted and
thrown away on the default mode, whatever its name or value, and since the
defaults are physical values rather than sentinels the result looks like a
legitimate spectrum: `Jrho_s` defaults to 0.0 and the third-order term is
proportional to it, so the typo returns exactly the `order=2` curve; at finite
field it returns a Zeeman step, nothing that looks empty. On `mode="DMRG"` the
same typos raise, only because the dynamical correlator downstream is strict.
This is `documentation.md` section 4.10's "`**kwargs` with no consumer". It came
in with the feature (`2d47de9`) and survived `b75003d`, which edited the same
docstring and branch; no test passes an unknown keyword to `get_kondo_spectrum`.

**Expected**: `TypeError` naming the unknown keyword on the ED route, as the
DMRG route already raises.

Repro, from the hunter (the paper's Fig. 3b/7b parameters, whose zero-bias peak
`tests/test_kondo_spectrum_paper_fig7.py` pins at 1.137 in the figure's units,
7.1463/2pi):

```python
# get_kondo_spectrum(mode="ED") accepts **kwargs and never reads them: a
# misspelled keyword is silently dropped. S=1/2, B=0, T=1 K, Jrho_s=-0.05,
# the paper's Fig. 3b/7b parameters (zero-bias peak 1.137 in this module's units).
import numpy as np
from dmrgpy import spinchain
sc = spinchain.Spin_Chain(["1/2"])
sc.set_hamiltonian(0*sc.Sz[0])
eVs = np.array([-4e-3, 0.0, 4e-3])
_, right = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0, omega0=20e-3)
_, typo = sc.get_kondo_spectrum(eVs, site=0, Jrho=-0.05, T=1.0, omega0=20e-3)
_, none = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=0.0, T=1.0, omega0=20e-3)
print("Jrho_s=-0.05     :", np.round(right, 4))
print("Jrho=-0.05 (typo):", np.round(typo, 4))
print("Jrho_s=0         :", np.round(none, 4))
# a DMRG-only keyword on the ED route is dropped the same way
_, ed_delta = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0, omega0=20e-3, delta=1e-3, submode="bogus")
print("ED with delta=1e-3, submode='bogus':", np.round(ed_delta, 4))
# the same typo on mode="DMRG" is not silent
try:
    sc3 = spinchain.Spin_Chain(["1/2"]*3, itensor_version=3)
    sc3.set_hamiltonian(0.01*(sc3.Sx[0]*sc3.Sx[1]+sc3.Sy[0]*sc3.Sy[1]+sc3.Sz[0]*sc3.Sz[1]))
    sc3.get_kondo_spectrum(eVs, site=0, Jrho=-0.05, T=0.0, order=2, mode="DMRG",
                           es=np.linspace(-0.03, 0.03, 600), delta=2e-4)
    print("DMRG typo: accepted")
except Exception as e:
    print("DMRG typo:", type(e).__name__, str(e)[:120])
```

```
Jrho_s=-0.05     : [5.5574 7.1463 5.5574]
Jrho=-0.05 (typo): [4.7124 4.7124 4.7124]
Jrho_s=0         : [4.7124 4.7124 4.7124]
ED with delta=1e-3, submode='bogus': [5.5574 7.1463 5.5574]
DMRG typo: TypeError Unexpected keyword argument(s) for the KPM dynamical correlator: Jrho. Solver parameters are attributes of the chain (se
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced to every printed digit, by exact
equality against the `Jrho_s=0` and `order=2` calls on a deterministic full
diagonalization, so no convergence or tolerance enters. Any other misspelled
name behaves the same way: `u=` for `U=` returns exactly the `U=0` curve (the 10
T bias asymmetry of -0.0125 becomes exactly 0.0), `t=0.3` for `T=0.3` quietly
runs at the default T=1 K, `omega_0=` quietly runs at the default cutoff, and
DMRG-route keywords and plain junk are all accepted unchanged. Struck: the
hunter's "why it survived" reason, that the ED branch discards `kwargs` by
design so a call written for the DMRG route also runs on ED; nothing records
that intent (not the commit, not PR #70, not the docstring, which says the
extra keywords "are forwarded" in its `mode="DMRG"` paragraph), and no in-tree
caller relies on it. The companion docstring defect stands in a sharper form:
`spinchain.py:215` to `218` names `n`, the number of KPM moments, as its example
of a forwarded keyword, but `n` never had an effect (at `2d47de9`
`kpmdmrg.get_dynamical_correlator` overwrote it) and since `1b87543` it raises
under both `submode="KPM"` and `"CVM"`, while a keyword the correlator really
has, `kernel=`, does go through. The reviewer's probe:

```python
# Reviewer probe 1: is get_kondo_spectrum(mode="ED") a sink for every
# unknown keyword, and does the dropped keyword leave a plausible curve?
import numpy as np
import dmrgpy
from dmrgpy import spinchain
print("dmrgpy from:", dmrgpy.__file__)

G, MUB, TP = 2.0, 5.7883818066e-5, 2*np.pi   # TP: figure units, as in tests/test_kondo_spectrum_paper_fig7.py

def chain(B):
    sc = spinchain.Spin_Chain(["1/2"])
    sc.set_hamiltonian(G*MUB*B*sc.Sz[0])
    return sc

eVs = np.array([-4e-3, -1e-3, -0.2e-3, 0.0, 0.2e-3, 1e-3, 4e-3])
base = dict(site=0, T=1.0, omega0=20e-3)

# (a) Jrho typo at B=0 is the Jrho_s=0 curve bit for bit
sc = chain(0.0)
_, right = sc.get_kondo_spectrum(eVs, Jrho_s=-0.05, **base)
_, typo  = sc.get_kondo_spectrum(eVs, Jrho=-0.05, **base)
_, none  = sc.get_kondo_spectrum(eVs, Jrho_s=0.0, **base)
_, o2    = sc.get_kondo_spectrum(eVs, order=2, **base)
print("(a) B=0  right/2pi     :", np.round(right/TP, 4))
print("    B=0  Jrho typo/2pi :", np.round(typo/TP, 4))
print("    max|typo-Jrho_s=0| =", np.max(np.abs(typo-none)),
      " max|typo-order2| =", np.max(np.abs(typo-o2)))
print("    third-order share of zero-bias peak:", round((right[3]-typo[3])/right[3], 4))

# (b) finite field: the typo returns a Zeeman step curve, not a flat or empty one
sc = chain(2.5)
_, right = sc.get_kondo_spectrum(eVs, Jrho_s=-0.05, **base)
_, typo  = sc.get_kondo_spectrum(eVs, Jrho=-0.05, **base)
print("(b) B=2.5T right/2pi   :", np.round(right/TP, 4))
print("    B=2.5T typo/2pi    :", np.round(typo/TP, 4))

# (c) a second misspelling, u= for U=, at 10 T where the odd term is nonzero
sc = chain(10.0)
_, U25  = sc.get_kondo_spectrum(eVs, Jrho_s=-0.05, U=0.25, **base)
_, utyp = sc.get_kondo_spectrum(eVs, Jrho_s=-0.05, u=0.25, **base)
_, U0   = sc.get_kondo_spectrum(eVs, Jrho_s=-0.05, U=0.0, **base)
print("(c) B=10T U=0.25/2pi   :", np.round(U25/TP, 4))
print("    B=10T u=0.25/2pi   :", np.round(utyp/TP, 4))
print("    max|u typo - U=0| =", np.max(np.abs(utyp-U0)),
      " asymmetry d(+4mV)-d(-4mV): U=0.25", round((U25[-1]-U25[0])/TP, 4),
      " u typo", round((utyp[-1]-utyp[0])/TP, 4))

# (d) misspelled temperature and cutoff fall back to the defaults T=1, omega0=20e-3
sc = chain(0.0)
_, T03  = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=0.3, omega0=20e-3)
_, t03  = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, t=0.3, omega0=20e-3)
_, T1   = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, omega0=20e-3)
print("(d) zero bias/2pi: T=0.3", round(T03[3]/TP, 4), " t=0.3 (typo)", round(t03[3]/TP, 4),
      " default T=1", round(T1[3]/TP, 4), " max|typo-default| =", np.max(np.abs(t03-T1)))
_, w10  = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0, omega0=10e-3)
_, w10t = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0, omega_0=10e-3)
print("    zero bias/2pi: omega0=10e-3", round(w10[3]/TP, 4), " omega_0=10e-3 (typo)",
      round(w10t[3]/TP, 4), " max|typo-default| =", np.max(np.abs(w10t-T1)))

# (e) nonsense keywords of any name are accepted
_, junk = sc.get_kondo_spectrum(eVs, site=0, Jrho_s=-0.05, T=1.0, omega0=20e-3,
                                delta=1e-3, submode="bogus", es="not an array",
                                dt2=-1, n=7, frobnicate=True)
print("(e) junk kwargs accepted, max|junk-right| =", np.max(np.abs(junk-T1)))
```

```
dmrgpy from: <repo>/src/dmrgpy/__init__.py
(a) B=0  right/2pi     : [0.8845 0.9804 1.1077 1.1374 1.1077 0.9804 0.8845]
    B=0  Jrho typo/2pi : [0.75 0.75 0.75 0.75 0.75 0.75 0.75]
    max|typo-Jrho_s=0| = 0.0  max|typo-order2| = 0.0
    third-order share of zero-bias peak: 0.3406
(b) B=2.5T right/2pi   : [0.8846 0.982  0.649  0.5349 0.649  0.982  0.8846]
    B=2.5T typo/2pi    : [0.75   0.7491 0.4464 0.367  0.4464 0.7491 0.75  ]
(c) B=10T U=0.25/2pi   : [1.1429 0.773  0.578  0.5732 0.5693 0.7088 1.1304]
    B=10T u=0.25/2pi   : [0.8866 0.4909 0.3237 0.3232 0.3237 0.4909 0.8866]
    max|u typo - U=0| = 0.0  asymmetry d(+4mV)-d(-4mV): U=0.25 -0.0125  u typo 0.0
(d) zero bias/2pi: T=0.3 1.2221  t=0.3 (typo) 1.1374  default T=1 1.1374  max|typo-default| = 0.0
    zero bias/2pi: omega0=10e-3 1.086  omega_0=10e-3 (typo) 1.1374  max|typo-default| = 0.0
(e) junk kwargs accepted, max|junk-right| = 0.0
```

```
versions: [3, 'python']
3 KPM {'n': 50} -> TypeError: Unexpected keyword argument(s) for the KPM dynamical correlator: n. Solver parameters are attributes of the ch
3 CVM {'n': 50} -> TypeError: dynamical_correlator() got an unexpected keyword argument 'n'
3 KPM {'Jrho': -0.05} -> TypeError: Unexpected keyword argument(s) for the KPM dynamical correlator: Jrho. Solver parameters are attributes of the
3 CVM {'Jrho': -0.05} -> TypeError: dynamical_correlator() got an unexpected keyword argument 'Jrho'
python KPM {'n': 50} -> TypeError: Unexpected keyword argument(s) for the KPM dynamical correlator: n. Solver parameters are attributes of the ch
python CVM {'n': 50} -> TypeError: dynamical_correlator() got an unexpected keyword argument 'n'
python KPM {'Jrho': -0.05} -> TypeError: Unexpected keyword argument(s) for the KPM dynamical correlator: Jrho. Solver parameters are attributes of the
python CVM {'Jrho': -0.05} -> TypeError: dynamical_correlator() got an unexpected keyword argument 'Jrho'
kernel=jackson: [0.00064 0.00025 0.00064]
kernel=lorentz: [0.16323 0.11175 0.16323]  max diff 0.1625897034207086  (0.2 s)
```

**Suggested fix**: in the `mode == "ED"` branch, raise `TypeError` on a non-empty
`kwargs`, naming the sorted keywords and saying they are `mode="DMRG"`
parameters or unknown, so a shared keyword dict fails loudly and the caller
splits it. The hunter's allow-list variant (accept the DMRG-route names
silently, raise on the rest) is not the right shape: `delta=` or
`submode="CVM"` on the default mode would stay inert, the same hole one level
down, and nothing in the tree needs ED to tolerate them. Fix the docstring at the
same time: drop the `n` example, name what the KPM correlator consumes
(`kernel`, `hodc_order`, `hodc_eta`) and say the moment count is set through the
chain attributes the kpmdmrg `TypeError` already points to. Regression: one
`pytest.raises(TypeError)` each for `Jrho=` and for `delta=` on `mode="ED"`, on a
1-site chain. No number changes for a correctly spelled call.

### 12. `get_kondo_spectrum(mode="DMRG")` at an accidental ground-state degeneracy returns the spectrum of whichever member of the manifold the random start converged to, anywhere in [1.0, 2.0] against the T->0+ value 1.5 that `mode="ED"` of the same method defines, documents and tests

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `recent-misc`

**Status**: FIXED -- `get_kondo_spectrum(mode="DMRG", n_gs=g)` averages every term (second order, third-order Kondo, potential) with equal weight over `get_excited_states(n=g, purify=True)`, each member set as the ground state in turn and the solved state restored afterwards in a `finally`; the default `n_gs=1` is the unchanged single-state path, and the docstring says which T=0 each mode takes. The reviewer's `set_gs` is only half of the switch: measured, `set_gs` alone leaves KPM on the solved state (both members returned 1.3206), because the KPM and DDMRG routes read the session's own `wf0`, so each switch also calls `session.set_wavefunction`; the session's band edges stay cached from the solved state, so every member is measured from the same E0. The third-order two-time term reads the Python-side `wf0`/`e0`, so it averages in the same pass (against a grid-consistent reference it matches the ED two-member average to 1.1e-11, where the two members differ by 6 at zero bias). At the crossing `n_gs=2` gives 1.50002 to 1.50004 in every run (`"python"` seeds 1 to 4, v3 seeds 0 to 2) where `n_gs=1` gives anywhere from 1.06 to 1.49; the chain's own state is restored exactly (a repeated default call agrees to 1.6e-13 on `"python"`, 3.2e-11 on v3); a `RuntimeWarning` fires when a member lies more than `delta` from E0; `n_gs>1` refuses `itensor_version=2` and `"julia_live"` with `NotImplementedError` (not run). Pinned by `::test_n_gs_averages_the_degenerate_manifold_like_ed`, `::test_n_gs_averages_the_two_time_kondo_term_too`, `::test_n_gs_1_never_touches_the_excited_states`, `::test_n_gs_must_be_a_positive_integer` and `::test_n_gs_needs_an_mps_session`. No number changes by default.

**Where**: `src/dmrgpy/spinchain.py::_get_kondo_spectrum_dmrg` (line 272),
which builds the second-order and potential terms from `get_dynamical_correlator`
starting from one ground-state vector; `kondospectrumtk/dmrgtwotime.py:142` and
`:146` (the third-order two-time term, from `chain.gs_energy()` and
`chain.get_gs()`); against `kondospectrumtk/edkondo.py:40` to `55`, which defines
T=0 as the equal-weight average over the degenerate ground manifold, pinned by
`tests/test_kondo_spectrum_paper_fig7.py::test_degenerate_ground_state_at_T0_is_averaged`.

`KondoSpectrum` defines T=0 as the T->0+ limit, the equal-weight average over a
degenerate ground manifold, and its test names the very case: an S=1 impurity
with D*Sz^2 + g*muB*B*Sz at the crossing g*muB*B = D. The DMRG route has no such
average and nothing says so: the docstring and the user guide's "T=0 and the
DMRG backend" say only that at T=0 "only the ground state is thermally
populated". Every Kondo test and example uses a Zeeman-split impurity or an
SU(2)-symmetric doublet, where at second order a pure state gives the ensemble
answer exactly, so only an accidental crossing shows it, which is where an STM
spectrum is taken on purpose.

**Expected**: the same T=0 on both modes of one method, either the manifold
average or a documented statement that the DMRG route uses one state.

Repro, from the hunter (an S=1 impurity at the |0>/|-1> crossing next to two
field-polarized spectators, so the DMRG backends have 3 sites):

```python
# T=0 Kondo spectrum at an accidental ground-state crossing: S=1 impurity,
# D*Sz^2 + g*muB*B*Sz with g*muB*B = D, so |0> and |-1> are degenerate.
# Two field-polarized S=1/2 spectators so the DMRG backends have 3 sites.
# mode="ED" averages over the manifold (KondoSpectrum's T->0+ limit);
# mode="DMRG" uses whatever ground state the solver converged to.
import sys
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk import conductance as C
version = sys.argv[1] if len(sys.argv) > 1 else "3"
version = int(version) if version.isdigit() else version
D = 1e-3
def chain():
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=version)
    h = D*sc.Sz[0]*sc.Sz[0] + D*sc.Sz[0]                       # crossing: |0>,|-1> degenerate
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*(sc.Sx[1]*sc.Sx[2] + sc.Sy[1]*sc.Sy[2] + sc.Sz[1]*sc.Sz[2])
    sc.set_hamiltonian(h)
    return sc
eVs = np.array([-3e-3, -1e-3, 0.0, 1e-3, 3e-3])
es = np.linspace(-30e-3, 30e-3, 6001)
sc = chain()
ks = KondoSpectrum(sc, 0, T=0.0)
print("lowest levels (meV):", np.round(ks.e[:4]*1e3, 6), " ground-manifold weights:", ks.p[:3])
_, ed = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="ED")
print("mode=ED   order=2 dI/dV/2pi:", np.round(ed/(2*np.pi), 4))
for run in range(4):
    np.random.seed(run)
    sc = chain()
    _, dm = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="DMRG",
                                  submode="KPM", delta=2e-5, es=es)
    gs = sc.get_gs()
    print("mode=DMRG run %d dI/dV/2pi:" % run, np.round(dm/(2*np.pi), 4),
          "  <Sz_0> of its ground state = %+.4f" % np.real(sc.vev(sc.Sz[0])))
```

```
lowest levels (meV): [0. 0. 2. 8.]  ground-manifold weights: [0.5 0.5 0. ]
mode=ED   order=2 dI/dV/2pi: [2.  1.5 1.5 1.5 2. ]
mode=DMRG run 0 dI/dV/2pi: [2.     1.0602 1.06   1.0602 2.    ]   <Sz_0> of its ground state = -0.0601
mode=DMRG run 1 dI/dV/2pi: [2.     1.0709 1.0707 1.0709 2.    ]   <Sz_0> of its ground state = -0.0709
mode=DMRG run 2 dI/dV/2pi: [2.     1.9242 1.9239 1.9242 2.    ]   <Sz_0> of its ground state = -0.9242
mode=DMRG run 3 dI/dV/2pi: [2.     1.6868 1.6864 1.6868 2.    ]   <Sz_0> of its ground state = -0.6867
```

```
lowest levels (meV): [0. 0. 2. 8.]  ground-manifold weights: [0.5 0.5 0. ]
mode=ED   order=2 dI/dV/2pi: [2.  1.5 1.5 1.5 2. ]
mode=DMRG run 0 dI/dV/2pi: [2.     1.8902 1.8898 1.8902 2.    ]   <Sz_0> of its ground state = -0.8901
mode=DMRG run 1 dI/dV/2pi: [2.     1.8355 1.8352 1.8355 2.    ]   <Sz_0> of its ground state = -0.8355
mode=DMRG run 2 dI/dV/2pi: [2.     1.2073 1.2071 1.2073 2.    ]   <Sz_0> of its ground state = -0.2073
mode=DMRG run 3 dI/dV/2pi: [2.     1.3262 1.3259 1.3262 2.    ]   <Sz_0> of its ground state = -0.3261
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced on both backends, and the cause is
the route, not the solver: on a chain forced onto the ED backend (`sc.mode="ED"`)
the same `mode="DMRG"` call returns 2.0 against the method's own `mode="ED"` 1.5,
deterministically, because ED hands it eigh's first vector. The DMRG state is
the exact ground state in every run (|E-E0| at most 2.3e-15, manifold weight 1
to twelve digits), and the zero-bias value equals 1 - <Sz_0> of the converged
state to 3e-4, so the possible range is exactly [1.0, 2.0], the values of the two
members themselves. Struck: "four seeds on v3", since `np.random.seed` does not
reach ITensor's `randomMPS`, so the v3 values cannot be reproduced from any seed
a caller can set (the reviewer got 1.26/1.98/1.08/1.92 and 1.49/1.95/1.38/1.97
where the hunter got 1.06/1.07/1.92/1.69; `"python"` reproduces the hunter to
every digit). "1.06 to 1.92" is a sample, not a bound. The third-order term on
the DMRG route inherits the behaviour by construction and was not measured.
Near a crossing the window is wider than the crossing itself: at a 1e-6 eV
splitting v3 at the default 15 sweeps returned a superposition in 1 of 3 runs,
and at 1e-9 eV neither backend resolves the lower level reliably, while
`KondoSpectrum`'s `degeneracy_tol` already commits ED to it. The reviewer's
probes:

```python
# Reviewer probe: is the single-state pick a property of the correlator
# ROUTE rather than of DMRG? (a) get_kondo_spectrum(mode="DMRG") on a chain
# whose sc.mode="ED" forces every correlator call onto the ED backend;
# (b) the easy-axis doublet re-measured off the |+-1>->|0> threshold
# (eV=+-1.2 meV instead of exactly +-1 meV).
import sys
import numpy as np
from dmrgpy import spinchain
version = sys.argv[1] if len(sys.argv) > 1 else "python"
version = int(version) if version.isdigit() else version
D = 1e-3
def chain(h0):
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=version)
    h = h0(sc) + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*(sc.Sx[1]*sc.Sx[2] + sc.Sy[1]*sc.Sy[2] + sc.Sz[1]*sc.Sz[2])
    sc.set_hamiltonian(h)
    return sc
es = np.linspace(-30e-3, 30e-3, 6001)
crossing = lambda sc: D*sc.Sz[0]*sc.Sz[0] + D*sc.Sz[0]
eVs = np.array([-1e-3, 0.0, 1e-3])
sc = chain(crossing)
sc.mode = "ED"
_, v = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="DMRG",
                             submode="KPM", delta=2e-5, es=es)
print("crossing, sc.mode='ED', get_kondo_spectrum(mode=DMRG) /2pi =", np.round(v/(2*np.pi), 4),
      "  ED-backend vev(Sz_0) = %+.4f" % np.real(sc.vev(sc.Sz[0], mode="ED")))
_, v = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="ED")
print("crossing, get_kondo_spectrum(mode=ED)                 /2pi =", np.round(v/(2*np.pi), 4))
doublet = lambda sc: -D*sc.Sz[0]*sc.Sz[0]
eVs = np.array([-1.2e-3, 0.0, 1.2e-3])
sc = chain(doublet)
_, ed = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="ED")
print("doublet, eV=(-1.2,0,1.2) meV, mode=ED /2pi =", np.round(ed/(2*np.pi), 4))
for run in range(3):
    np.random.seed(run)
    sc = chain(doublet)
    _, dm = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="DMRG",
                                  submode="KPM", delta=2e-5, es=es)
    print("   DMRG run %d: <Sz_0>=%+.4f  /2pi = %s" % (run, np.real(sc.vev(sc.Sz[0])), np.round(dm/(2*np.pi), 4)))
```

```
crossing, sc.mode='ED', get_kondo_spectrum(mode=DMRG) /2pi = [2.     1.9999 2.    ]   ED-backend vev(Sz_0) = -1.0000
crossing, get_kondo_spectrum(mode=ED)                 /2pi = [1.5 1.5 1.5]
doublet, eV=(-1.2,0,1.2) meV, mode=ED /2pi = [2. 1. 2.]
   DMRG run 0: <Sz_0>=-0.8747  /2pi = [2.     1.0001 2.    ]
   DMRG run 1: <Sz_0>=-0.4745  /2pi = [2.     1.0001 2.    ]
   DMRG run 2: <Sz_0>=-0.2457  /2pi = [2.     1.0001 2.    ]
```

```python
# Reviewer probe: is the DMRG ground state at the crossing exact (energy,
# variance, weight inside the ED manifold), does the DMRG dI/dV equal the
# pure-state value 1-<Sz_0> of THAT state, and what do the other routes of
# the same chain (ED vev, ED dynamical correlator submodes) do at T=0?
import sys
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
version = sys.argv[1] if len(sys.argv) > 1 else "python"
version = int(version) if version.isdigit() else version
D = 1e-3
def chain():
    sc = spinchain.Spin_Chain(["1", "1/2", "1/2"], itensor_version=version)
    h = D*sc.Sz[0]*sc.Sz[0] + D*sc.Sz[0]
    h = h + 10e-3*(sc.Sz[1] + sc.Sz[2]) + 2e-3*(sc.Sx[1]*sc.Sx[2] + sc.Sy[1]*sc.Sy[2] + sc.Sz[1]*sc.Sz[2])
    sc.set_hamiltonian(h)
    return sc
def manifold_projector(sc):
    # site 0 in {Sz=0,-1}: 1-(Sz^2+Sz)/2 ; spectators both down
    p0 = 1 - 0.5*(sc.Sz[0]*sc.Sz[0] + sc.Sz[0])
    return p0*(0.5 - sc.Sz[1])*(0.5 - sc.Sz[2])
eVs = np.array([-1e-3, 0.0, 1e-3])
es = np.linspace(-30e-3, 30e-3, 6001)
sc = chain()
edobj = sc.get_ED_obj()
emu, _ = edobj.get_diagonalized_hamiltonian()
emu = np.sort(np.array(emu, dtype=float))
E0 = emu[0]
print("ED: E0=%.15f  E1-E0=%.3e  E2-E0=%.3e (eV)" % (E0, emu[1]-E0, emu[2]-E0))
ks = KondoSpectrum(sc, 0, T=0.0)
print("KondoSpectrum T=0 weights p[:3] =", ks.p[:3])
print("ED vev(Sz_0)                      = %+.6f" % np.real(sc.vev(sc.Sz[0], mode="ED")))
_, ed = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="ED")
print("get_kondo_spectrum(mode=ED)       /2pi =", np.round(ed/(2*np.pi), 4))
for sm in ("ED", "KPM"):
    v = second_order_dIdV_dc(sc, 0, eVs, mode="ED", submode=sm, delta=2e-5, es=es)
    print("second_order_dIdV_dc(ED,%-3s)      /2pi =" % sm, np.round(v/(2*np.pi), 4))
for run in range(4):
    np.random.seed(run)
    sc = chain()
    _, dm = sc.get_kondo_spectrum(eVs, site=0, T=0.0, order=2, mode="DMRG",
                                  submode="KPM", delta=2e-5, es=es)
    e = sc.gs_energy()
    var = sc.gs_energy_fluctuation()
    w = np.real(sc.vev(manifold_projector(sc)))
    sz = np.real(sc.vev(sc.Sz[0]))
    print("DMRG run %d: E-E0=%+.2e  sqrt(var)=%.2e  manifold weight=%.12f  <Sz_0>=%+.4f"
          "  dIdV/2pi=%s  1-<Sz_0>=%.4f  max|diff|=%.1e"
          % (run, e-E0, var, w, sz, np.round(dm/(2*np.pi), 4), 1-sz,
             np.max(np.abs(dm/(2*np.pi) - (1-sz)))))
```

```
ED: E0=-0.009500000000000  E1-E0=0.000e+00  E2-E0=2.000e-03 (eV)
KondoSpectrum T=0 weights p[:3] = [0.5 0.5 0. ]
ED vev(Sz_0)                      = -1.000000
get_kondo_spectrum(mode=ED)       /2pi = [1.5 1.5 1.5]
second_order_dIdV_dc(ED,ED )      /2pi = [1.5034 1.5024 1.5034]
second_order_dIdV_dc(ED,KPM)      /2pi = [2.     1.9999 2.    ]
DMRG run 0: E-E0=-1.73e-18  sqrt(var)=0.00e+00  manifold weight=1.000000000000  <Sz_0>=-0.4949  dIdV/2pi=[1.495  1.4947 1.495 ]  1-<Sz_0>=1.4949  max|diff|=2.4e-04
DMRG run 1: E-E0=+0.00e+00  sqrt(var)=1.65e-10  manifold weight=1.000000000000  <Sz_0>=-0.9463  dIdV/2pi=[1.9464 1.946  1.9464]  1-<Sz_0>=1.9463  max|diff|=3.2e-04
DMRG run 2: E-E0=-1.73e-18  sqrt(var)=0.00e+00  manifold weight=1.000000000000  <Sz_0>=-0.3832  dIdV/2pi=[1.3833 1.383  1.3833]  1-<Sz_0>=1.3832  max|diff|=2.3e-04
DMRG run 3: E-E0=-1.73e-18  sqrt(var)=1.65e-10  manifold weight=1.000000000000  <Sz_0>=-0.9745  dIdV/2pi=[1.9745 1.9742 1.9745]  1-<Sz_0>=1.9745  max|diff|=3.2e-04
```

**Suggested fix**: every term is linear in the ground-state density matrix, so
averaging over an orthonormal basis of the manifold gives Tr[P0 X]/g whichever
basis DMRG found, exactly the ED definition. Detecting the degeneracy from
energies and raising is the wrong shape, since no tolerance separates "exactly
degenerate" from "split below what the sweep resolved" on DMRG. The better shape
is a caller-supplied manifold size (`n_gs=`, the same contract as `dex`),
averaging over `get_excited_states(n=n_gs, purify=True)` with each member set as
the ground state in turn; the default path stays byte-identical. Whichever fix
lands, one sentence in the docstring and one in the user guide should say that
the correlator route's T=0 is the single converged state, not the manifold
average. Numbers change only where a caller opts in.

### 13. The potential-interference DMRG term weights every point of `es` with its first spacing, so a non-uniform grid mis-weights the correlator by the ratio of the first spacing to the local one: 101.7 times the exact value on a smooth grid dense at the line, while the second-order sibling on the same grid is right

`bug` &middot; severity **LOW to MEDIUM** &middot; CONFIRMED &middot; lens `recent-misc`

**Status**: FIXED -- `_convolved_F0_weight` uses trapezoid weights of the actual grid, computed once, with `np.einsum('w,w,ew->e', wts, S, kernel)`. Pinned by `tests/test_audit_2026_09_24_kondo.py::test_potential_term_on_a_non_uniform_grid_matches_the_exact_sum` (within 1e-3 of the peak of the exact sum, measured 7.6e-4, and within 1e-5 of the uniform-grid result, measured 1.6e-7). NUMBERS CHANGE only on a non-uniform `es`: on the test's single S=1/2 at 10 T (delta=2e-6, Jrho_s=0.1, U=0.3) the peak at eV=1.10 meV on the 12000-point sinh-mapped grid goes from 68.2332 to 0.6704 against an exact 0.6709; on a uniform grid the result moves by 2.40e-8, the two endpoint half-bins.

**Where**: `src/dmrgpy/kondospectrumtk/potentialdc.py:39` and `:44`
(`dw = x[1] - x[0]`, `return dw*np.einsum('w,ew->e', S, kernel)`), inside
`_convolved_F0_weight`; reached by `third_order_potential_dIdV_dc` and by
`get_kondo_spectrum(mode="DMRG", order=3, U!=0)`, which hands the same `es` to
this term and to `secondorder_dc`, whose cumulative trapezoid handles a
non-uniform grid correctly.

Nothing documents a uniform grid and nothing checks for one: the private
docstring uses an undefined `dw`, and the public function defers to
`second_order_dIdV_dc`'s docstring, which asks only for "several times finer
spacing than delta". The grid reaches the quadrature unchanged (the ED Lehmann
route evaluates at every requested point and KPM interpolates onto `es`), so the
weighting is the only thing that goes wrong. It survived because the one test
and the one example that call this function both pass an `np.linspace`.

**Expected**: the same answer on any grid that resolves delta near the lines.

Repro, from the hunter (the test's own single spin, exact Lehmann correlator):

```python
# The test's own single-spin system (tests/test_kondo_spectrum_potentialdc.py),
# its own uniform es, and a non-uniform es over the same range that is
# refined around the Zeeman line; exact Lehmann correlator (submode="ED").
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV, second_order_dIdV
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
G = 2.0; MUB = 5.7883818066e-5
sc = spinchain.Spin_Chain(["1/2"])
sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
ks = KondoSpectrum(sc, site=0, T=0.0)
eVs = np.linspace(-1e-3, 2e-3, 21)
delta, Jrho_s, U = 2e-6, 0.1, 0.3
ref = third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0)
ref2 = second_order_dIdV(ks, eVs, T0=1.0, U=U)
uni = np.linspace(-1e-3, 3e-3, 40_000)                 # the test's grid, 1e-7
non = np.concatenate([np.linspace(-1e-3, 1.0e-3, 4000, endpoint=False),   # 5e-7
                      np.linspace(1.0e-3, 1.4e-3, 4000, endpoint=False),  # 1e-7 around the 1.158 meV line
                      np.linspace(1.4e-3, 3e-3, 3200)])                   # 5e-7
for label, es in (("uniform (the test's es)", uni), ("non-uniform, refined at the line", non)):
    p = third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, T0=1.0, mode="ED", submode="ED", delta=delta, es=es)
    s = second_order_dIdV_dc(sc, 0, eVs, T0=1.0, U=U, mode="ED", submode="ED", delta=delta, es=es)
    print("%-34s potential max|dc-exact| = %.4f (peak %.4f)   second order max|dc-exact| = %.4f (peak %.4f)"
          % (label, np.max(np.abs(p-ref)), np.max(np.abs(ref)), np.max(np.abs(s-ref2)), np.max(np.abs(ref2))))
```

```
uniform (the test's es)            potential max|dc-exact| = 0.0005 (peak 0.6709)   second order max|dc-exact| = 0.0318 (peak 6.9743)
non-uniform, refined at the line   potential max|dc-exact| = 2.6735 (peak 0.6709)   second order max|dc-exact| = 0.0318 (peak 6.9743)
```

**Reviewer (CONFIRMED)**: reproduced exactly. On the test's system the shipped
result at the peak is the exact one times the ratio of the first spacing to the
spacing at the 1.1577 meV line, in both directions: 4.985 on a coarse-first
grid, 0.200 on a fine-first one, and 101.7 on a smooth sinh-mapped grid dense at
the line, the one a user refining near a line would most likely write. With
`np.trapezoid` in its place all four land on the exact answer (max error 0.0005,
the same as the uniform grid), and the uniform grid moves by 2.4e-8. Struck: the
hunter's motivation that a uniform grid would need about 2e7 points; +-40 meV at
delta/4 is 1.6e5 points, measured at 1.21 s and 325 MB, so a local refinement is
a natural choice, not a forced one.

```python
# Candidate 1 probe: potentialdc._convolved_F0_weight weights every es point
# with x[1]-x[0]. Exact Lehmann correlator (mode="ED", submode="ED") so only
# the quadrature differs. For each grid: the shipped code, and the same
# function with the rectangle rule swapped for np.trapezoid (monkeypatched
# in-process only, the repo is untouched).
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import potentialdc
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV
from dmrgpy.kondospectrumtk.stepfunctions import F0

shipped = potentialdc._convolved_F0_weight

def trapezoid_weight(chain, op, eVs, omega0, Gamma0, mode, submode, delta, es, **kw):
    x, S = chain.get_dynamical_correlator(mode=mode, submode=submode,
            name=(op.get_dagger(), op), delta=delta, es=es, **kw)
    x = np.asarray(x, dtype=float); S = np.asarray(S).real
    d = eVs[:, None] - x[None, :]; s = eVs[:, None] + x[None, :]
    k = F0(d, omega0=omega0, Gamma0=Gamma0) - F0(s, omega0=omega0, Gamma0=Gamma0)
    return np.trapezoid(S[None, :]*k, x, axis=1)

def run(sc, eVs, Jrho_s, U, delta, es, fn):
    potentialdc._convolved_F0_weight = fn
    try:
        return potentialdc.third_order_potential_dIdV_dc(
            sc, 0, eVs, Jrho_s, U, T0=1.0, mode="ED", submode="ED", delta=delta, es=es)
    finally:
        potentialdc._convolved_F0_weight = shipped

G = 2.0; MUB = 5.7883818066e-5
print("=== single spin, the test's system: S=1/2, 10 T, delta=2e-6, Jrho_s=0.1, U=0.3")
sc = spinchain.Spin_Chain(["1/2"])
sc.set_hamiltonian(G*MUB*10.0*sc.Sz[0])
ks = KondoSpectrum(sc, site=0, T=0.0)
print("line at %.4f meV" % (ks.e[1]*1e3))
eVs = np.linspace(-1e-3, 2e-3, 21)
ref = third_order_potential_dIdV(ks, eVs, 0.1, 0.3, T0=1.0)
ipk = np.argmax(np.abs(ref))
c = ks.e[1]
grids = {
 "uniform (the test's es, 1e-7)": np.linspace(-1e-3, 3e-3, 40_000),
 "coarse-first: 5e-7, 1e-7 at the line, 5e-7": np.concatenate([
     np.linspace(-1e-3, 1.0e-3, 4000, endpoint=False),
     np.linspace(1.0e-3, 1.4e-3, 4000, endpoint=False),
     np.linspace(1.4e-3, 3e-3, 3200)]),
 "fine-first: 1e-7 first 1000 pts, then 5e-7": np.concatenate([
     np.linspace(-1e-3, -0.9e-3, 1000, endpoint=False),
     np.linspace(-0.9e-3, 3e-3, 7801)]),
 # smooth sinh map concentrating points at the line: spacing ~2e-8 at the
 # line, ~2.6e-6 at the ends -- resolves delta=2e-6 near the line only
 "sinh-mapped, dense at the line": c + 2e-5*np.sinh(np.linspace(np.arcsinh(-2.158e-3/2e-5), np.arcsinh(1.842e-3/2e-5), 12000)),
}
for label, es in grids.items():
    dx = np.diff(es)
    a = run(sc, eVs, 0.1, 0.3, 2e-6, es, shipped)
    b = run(sc, eVs, 0.1, 0.3, 2e-6, es, trapezoid_weight)
    print("%-45s n=%5d dx[0]=%.1e dx range [%.1e,%.1e]" % (label, len(es), dx[0], dx.min(), dx.max()))
    print("    shipped:    max|dc-exact| = %.4f  (peak %.4f, dc/exact at peak = %.3f)"
          % (np.max(np.abs(a-ref)), np.abs(ref[ipk]), a[ipk]/ref[ipk]))
    print("    trapezoid:  max|dc-exact| = %.4f  (dc/exact at peak = %.4f);  max|shipped-trapezoid| = %.3e"
          % (np.max(np.abs(b-ref)), b[ipk]/ref[ipk], np.max(np.abs(a-b))))

print("=== 3-site chain (the example's Hamiltonian), delta=2e-5, Jrho_s=0.05, U=0.2, es over +-40 meV")
sc3 = spinchain.Spin_Chain(["1/2"]*3)
h = G*MUB*10.0*sc3.Sz[0]
for i in range(2):
    h = h + 0.01*(sc3.Sx[i]*sc3.Sx[i+1] + sc3.Sy[i]*sc3.Sy[i+1] + sc3.Sz[i]*sc3.Sz[i+1])
sc3.set_hamiltonian(h)
ks3 = KondoSpectrum(sc3, site=0, T=0.0)
eVs3 = np.linspace(-2e-3, 2e-3, 41)
ref3 = third_order_potential_dIdV(ks3, eVs3, 0.05, 0.2, T0=1.0)
g3 = {
 "uniform 2e-6": np.linspace(-40e-3, 40e-3, 40001),
 "5e-6 outside +-3 meV, 2e-6 inside": np.concatenate([
     np.linspace(-40e-3, -3e-3, 7400, endpoint=False),
     np.linspace(-3e-3, 3e-3, 3000, endpoint=False),
     np.linspace(3e-3, 40e-3, 7401)]),
}
res = {}
for label, es in g3.items():
    a = run(sc3, eVs3, 0.05, 0.2, 2e-5, es, shipped)
    b = run(sc3, eVs3, 0.05, 0.2, 2e-5, es, trapezoid_weight)
    res[label] = (a, b)
    print("%-36s shipped max|dc-exact| = %.4f   trapezoid max|dc-exact| = %.4f   (peak %.4f)"
          % (label, np.max(np.abs(a-ref3)), np.max(np.abs(b-ref3)), np.max(np.abs(ref3))))
(au, bu), (an, bn) = res.values()
print("grid-only change: shipped moves by %.4f, trapezoid moves by %.2e" % (np.max(np.abs(an-au)), np.max(np.abs(bn-bu))))
```

```
=== single spin, the test's system: S=1/2, 10 T, delta=2e-6, Jrho_s=0.1, U=0.3
line at 1.1577 meV
uniform (the test's es, 1e-7)                 n=40000 dx[0]=1.0e-07 dx range [1.0e-07,1.0e-07]
    shipped:    max|dc-exact| = 0.0005  (peak 0.6709, dc/exact at peak = 0.999)
    trapezoid:  max|dc-exact| = 0.0005  (dc/exact at peak = 0.9992);  max|shipped-trapezoid| = 2.401e-08
coarse-first: 5e-7, 1e-7 at the line, 5e-7    n=11200 dx[0]=5.0e-07 dx range [1.0e-07,5.0e-07]
    shipped:    max|dc-exact| = 2.6735  (peak 0.6709, dc/exact at peak = 4.985)
    trapezoid:  max|dc-exact| = 0.0005  (dc/exact at peak = 0.9992);  max|shipped-trapezoid| = 2.674e+00
fine-first: 1e-7 first 1000 pts, then 5e-7    n= 8801 dx[0]=1.0e-07 dx range [1.0e-07,5.0e-07]
    shipped:    max|dc-exact| = 0.5369  (peak 0.6709, dc/exact at peak = 0.200)
    trapezoid:  max|dc-exact| = 0.0005  (dc/exact at peak = 0.9992);  max|shipped-trapezoid| = 5.364e-01
sinh-mapped, dense at the line                n=12000 dx[0]=1.9e-06 dx range [1.8e-08,1.9e-06]
    shipped:    max|dc-exact| = 67.5622  (peak 0.6709, dc/exact at peak = 101.698)
    trapezoid:  max|dc-exact| = 0.0005  (dc/exact at peak = 0.9992);  max|shipped-trapezoid| = 6.756e+01
=== 3-site chain (the example's Hamiltonian), delta=2e-5, Jrho_s=0.05, U=0.2, es over +-40 meV
uniform 2e-6                         shipped max|dc-exact| = 0.0067   trapezoid max|dc-exact| = 0.0067   (peak 0.1110)
5e-6 outside +-3 meV, 2e-6 inside    shipped max|dc-exact| = 0.1433   trapezoid max|dc-exact| = 0.0067   (peak 0.1110)
grid-only change: shipped moves by 0.1500, trapezoid moves by 9.59e-11
```

**Suggested fix**: trapezoid weights computed once (`wts[1:] += dx/2;
wts[:-1] += dx/2`) and `np.einsum('w,w,ew->e', wts, S, kernel)`, the same
quadrature as `secondorder_dc` without materializing a second
(2 n_eV x n_es) array. Regression on the test's own system with the sinh-mapped
grid, which fails the shipped code by a factor 100. Numbers change only on a
non-uniform `es`, where they were wrong; on a uniform grid by the two endpoint
half-bins, about 2e-8 here.

### 14. The only DMRG example of the potential term violates the `es` coverage its own library docstring requires, and blames the resulting 9.7 per cent DMRG-versus-ED gap on KPM broadening; the second-order docstring the library points to states a weaker requirement

`docs` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `recent-misc`

**Status**: FIXED -- the example's `es` now spans +-20 meV at the old spacing (5328 points), its assert is 0.05 of the peak, its sweep schedule is pinned, and its comments name the 0.77 meV lowest transition and the 10 to 16 meV lines, with the same wrong "~1.16meV Zeeman gap" comment fixed in `kondo_third_order_dmrg_VS_ED` (result unchanged there, 0.013 on 2.148); `potentialdc` and `secondorder_dc` now say which term needs which coverage, and `get_kondo_spectrum` says a shared `es` must meet the stricter one; the potential term issues a `RuntimeWarning` when sum_k int S_kk misses more than 1e-2 of S(S+1), checked on `Spin_Chain` sites where S(S+1) = (d^2-1)/4, which counts the elastic weight at omega=0, so the docstring says `es` must include omega=0. Pinned by `tests/test_audit_2026_09_24_kondo.py::test_potential_term_warns_when_es_misses_spectral_weight` and by the example's own assert. NUMBERS CHANGE in the example only: max|DMRG-ED| goes from 0.0108 to 0.0023 on a 0.1110 peak (9.7 to 2.1 per cent), measured on v3 before the kpm cluster's moment-count cut landed.

**Where**: `examples/kondo/kondo_potential_term_dmrg_VS_ED/main.py` (lines 41
and 57 to 60: `es` over +-3 meV, the comment "must cover the ~1.16meV Zeeman
gap", the residual called "KPM delta-broadening error", an assert at 0.15 of the
peak); the same wrong gap comment at
`examples/kondo/kondo_third_order_dmrg_VS_ED/main.py:62` (harmless there, U=0);
`potentialdc.py:60` to `63`, which requires `es` to "cover every eigenstate
transition energy from the ground state" and says it is needed "for the same
reason" as in `second_order_dIdV_dc`, whose docstring asks only for the
transitions "that matter for the eVs sweep"; `get_kondo_spectrum`, which hands
one `es` to both terms and points the user at "either function's docstring".

The chain's S_k-active transitions run from 0.77 to 15.98 meV, and those above 3
meV (10.36 to 15.64 meV) carry 0.40 of the total weight S(S+1)=0.75. The
potential-interference kernel F0(eV-w) - F0(eV+w) is the paper's eq. 22 closed
form, which the ED reference sums over every level, and inside the band it falls
off only as 2 eV omega0/(w(omega0+w)), so those high lines contribute 0.01078 at
eV=+-2 meV, computed analytically with no correlator call, which is the whole of
the example's 0.0108 gap on a 0.1110 peak. The second-order Theta0 kernel, by
contrast, vanishes for every level above max|eV|, which is why its docstring's
weaker requirement is right for it and wrong here.

**Expected**: an example whose `es` meets the stricter requirement, and
docstrings that say which term needs which.

Repro, from the hunter:

```python
# Is the 10% DMRG-vs-ED gap of the potential-term example KPM broadening,
# or the es grid (+-3 meV) not covering the 10-15 meV transitions?
import sys, time
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV, second_order_dIdV
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
import dmrgpy
print("dmrgpy from", dmrgpy.__file__)
version = sys.argv[1] if len(sys.argv) > 1 else "3"
version = int(version) if version.isdigit() else version
G = 2.0; MUB = 5.7883818066e-5
sc = spinchain.Spin_Chain(["1/2"]*3, itensor_version=version)
h = G*MUB*10.0*sc.Sz[0]
for i in range(2):
    h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
sc.get_gs()
ks = KondoSpectrum(sc, site=0, T=0.0)
print("transition energies from GS (meV):", np.round(ks.e*1e3, 4))
# weight of each transition for sum_k |<m|Sk|GS>|^2
Xi = np.stack([ks.Sx, ks.Sy, ks.Sz], axis=-1)
loop = np.real(np.einsum('mk,mk->m', Xi[:, 0, :], np.conj(Xi[:, 0, :])))
print("sum_k |<m|Sk|GS>|^2 per m:", np.round(loop, 4))

eVs = np.linspace(-2e-3, 2e-3, 41)
Jrho_s, U = 0.05, 0.2
ref = third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0)
# the same exact sum restricted to transitions below 3 meV
class Restricted: pass
ksr = Restricted()
keep = ks.e < 3e-3
for a in ("T", "kB", "p"): setattr(ksr, a, getattr(ks, a))
ksr.e = ks.e[keep]; ksr.dim = int(keep.sum()); ksr.p = ks.p[keep]
for a in ("Sx", "Sy", "Sz"): setattr(ksr, a, getattr(ks, a)[np.ix_(keep, keep)])
ref_lo = third_order_potential_dIdV(ksr, eVs, Jrho_s, U, T0=1.0)
print("max|ED exact| = %.4f" % np.max(np.abs(ref)))
print("max|ED exact - ED restricted to eps<3meV| = %.4f" % np.max(np.abs(ref-ref_lo)))
for label, es in (("es +-3meV, 800 pts (example)", np.linspace(-3e-3, 3e-3, 800)),
                  ("es +-40meV, 10667 pts (same spacing)", np.linspace(-40e-3, 40e-3, 10667))):
    t = time.time()
    d = third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, T0=1.0, mode="DMRG",
                                       submode="KPM", delta=2e-5, es=es)
    print("%-38s max|DMRG-ED| = %.4f  max|DMRG-EDrestricted| = %.4f  (%.1fs)"
          % (label, np.max(np.abs(d-ref)), np.max(np.abs(d-ref_lo)), time.time()-t))
```

```
dmrgpy from <repo>/src/dmrgpy/__init__.py
transition energies from GS (meV): [ 0.      0.7696 10.364  10.3716 14.8222 15.2602 15.6406 15.9798]
sum_k |<m|Sk|GS>|^2 per m: [0.1286 0.2201 0.0613 0.1716 0.0707 0.0601 0.0375 0.    ]
max|ED exact| = 0.1110
max|ED exact - ED restricted to eps<3meV| = 0.0108
es +-3meV, 800 pts (example)           max|DMRG-ED| = 0.0108  max|DMRG-EDrestricted| = 0.0023  (2.5s)
es +-40meV, 10667 pts (same spacing)   max|DMRG-ED| = 0.0023  max|DMRG-EDrestricted| = 0.0108  (2.7s)
```

**Reviewer (CONFIRMED, NARROWED)**: struck as a library defect, "the potential
term silently drops every transition outside `es`", since the requirement is
stated in full in its docstring; what is wrong is the example that violates it
and the cross-reference that makes violating it easy. Corrected: "F0 decays only
as omega0/|x|" describes the behaviour beyond omega0; the 10 to 16 meV lines sit
inside the band, where the kernel difference is about 2eV/w, a formula matching
the exact kernel to 0.5 per cent at eV=1 meV, and the tail is the model the ED
reference sums, not an artefact. Corrected: the 0.0108 equals "ED restricted to
levels below 3 meV minus full ED" in the max norm (both at the sweep edge), not
pointwise; with `es` widened to +-40 meV the remaining 0.0023 sits at eV=-0.80
meV next to the 0.77 meV line, below 5e-5 more than 0.3 meV away from it, and
falls to 0.0006 at delta=1e-5, so the example's comment is right only about that
remainder.

```python
# Candidate 2 probe: the example's 3-site chain. (1) an anchor that calls
# no dmrgpy correlator at all: the analytic per-line contribution of the
# transitions above 3 meV; (2) the exact Lehmann correlator (submode="ED")
# on the example's es and on a +-40 meV es, so KPM is out of the picture;
# (3) KPM on itensor_version=3 on both grids; (4) the sum-rule deficit
# S(S+1) - sum_k int_es S_kk on each, i.e. what the hunter's suggested
# detector would read; (5) second-order term on both grids.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum
from dmrgpy.kondospectrumtk.conductance import third_order_potential_dIdV, second_order_dIdV
from dmrgpy.kondospectrumtk.potentialdc import third_order_potential_dIdV_dc
from dmrgpy.kondospectrumtk.secondorder_dc import second_order_dIdV_dc
from dmrgpy.kondospectrumtk.stepfunctions import F0

G = 2.0; MUB = 5.7883818066e-5
omega0, Gamma0 = 20e-3, 5e-6
Jrho_s, U = 0.05, 0.2
sc = spinchain.Spin_Chain(["1/2"]*3, itensor_version=3)
h = G*MUB*10.0*sc.Sz[0]
for i in range(2):
    h = h + 0.01*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
sc.get_gs()
ks = KondoSpectrum(sc, site=0, T=0.0)
Xi = np.stack([ks.Sx, ks.Sy, ks.Sz], axis=-1)
w = np.real(np.einsum('mk,mk->m', Xi[:, 0, :], np.conj(Xi[:, 0, :])))
print("transition energies (meV):", np.round(ks.e*1e3, 4))
print("weights sum_k|<m|Sk|GS>|^2:", np.round(w, 4), " total %.6f (S(S+1)=0.75)" % w.sum())

eVs = np.linspace(-2e-3, 2e-3, 41)
ref = third_order_potential_dIdV(ks, eVs, Jrho_s, U, T0=1.0)
ref2 = second_order_dIdV(ks, eVs, T0=1.0, U=U)
pre = 4*np.pi*Jrho_s*U

# (1) analytic anchor: for eV>0, Theta0(-eV)=0, so dIdV = pre*sum_m w_m [F0(eV-e_m)-F0(eV+e_m)]
hi = ks.e > 3e-3
print("\n(1) analytic contribution of the lines above 3 meV (no correlator involved)")
for eV in (0.5e-3, 1e-3, 2e-3):
    k = F0(eV-ks.e[hi], omega0, Gamma0) - F0(eV+ks.e[hi], omega0, Gamma0)
    approx = 2*eV*omega0/(ks.e[hi]*(omega0+ks.e[hi]))
    print("  eV=%.1f meV: pre*sum w_m k_m = %.5f ; kernel k_m = %s ; 2eV*w0/(w(w0+w)) = %s"
          % (eV*1e3, pre*np.sum(w[hi]*k), np.round(k, 4), np.round(approx, 4)))
print("  max over the eVs sweep of the exact full sum minus exact sum without those lines:")
k_all = (F0(eVs[:, None]-ks.e[None, hi], omega0, Gamma0) - F0(eVs[:, None]+ks.e[None, hi], omega0, Gamma0))
contrib = pre*np.einsum('m,em->e', w[hi], k_all)  # h(eV)-h(-eV) = weighted(eV) for eV!=0 (weighted is odd)
contrib[eVs == 0] = 0.
print("  max|analytic high-line contribution| = %.5f at eV=%.2f meV (peak of exact term %.5f)"
      % (np.max(np.abs(contrib)), eVs[np.argmax(np.abs(contrib))]*1e3, np.max(np.abs(ref))))

grids = {"example es +-3 meV, 800": np.linspace(-3e-3, 3e-3, 800),
         "+-40 meV, 10667 (same spacing)": np.linspace(-40e-3, 40e-3, 10667)}
ops = (sc.Sx[0], sc.Sy[0], sc.Sz[0])
print("\n(2)-(5)")
out = {}
for label, es in grids.items():
    for mode, submode in (("ED", "ED"), ("DMRG", "KPM")):
        p = third_order_potential_dIdV_dc(sc, 0, eVs, Jrho_s, U, T0=1.0, mode=mode,
                                           submode=submode, delta=2e-5, es=es)
        s2 = second_order_dIdV_dc(sc, 0, eVs, T0=1.0, U=U, mode=mode, submode=submode,
                                   delta=2e-5, es=es)
        tot = 0.
        for op in ops:
            x, S = sc.get_dynamical_correlator(mode=mode, submode=submode,
                    name=(op.get_dagger(), op), delta=2e-5, es=es)
            tot += np.trapezoid(np.asarray(S).real, np.asarray(x, dtype=float))
        out[(label, submode)] = p
        print("  %-32s %-8s potential max|dc-exact| = %.4f (%.1f%% of peak)  pointwise max|dc-exact-analytic_high| = %.4f"
              % (label, mode+"/"+submode, np.max(np.abs(p-ref)), 100*np.max(np.abs(p-ref))/np.max(np.abs(ref)),
                 np.max(np.abs(p-ref+contrib))))
        print("  %-32s %-8s second order max|dc-exact| = %.2e (peak %.4f)   sum-rule deficit 0.75-int S = %.3e (rel %.2e)"
              % ("", "", np.max(np.abs(s2-ref2)), np.max(np.abs(ref2)), 0.75-tot, (0.75-tot)/0.75))
print("\n  KPM minus exact-Lehmann on the same grid: +-3 meV %.4f, +-40 meV %.4f"
      % tuple(np.max(np.abs(out[(l, "KPM")]-out[(l, "ED")])) for l in grids))
```

```
transition energies (meV): [ 0.      0.7696 10.364  10.3716 14.8222 15.2602 15.6406 15.9798]
weights sum_k|<m|Sk|GS>|^2: [0.1286 0.2201 0.0613 0.1716 0.0707 0.0601 0.0375 0.    ]  total 0.750000 (S(S+1)=0.75)

(1) analytic contribution of the lines above 3 meV (no correlator involved)
  eV=0.5 meV: pre*sum w_m k_m = 0.00266 ; kernel k_m = [0.0636 0.0636 0.0388 0.0372 0.0359 0.0348] ; 2eV*w0/(w(w0+w)) = [0.0636 0.0635 0.0387 0.0372 0.0359 0.0348]
  eV=1.0 meV: pre*sum w_m k_m = 0.00533 ; kernel k_m = [0.1277 0.1276 0.0777 0.0745 0.0719 0.0697] ; 2eV*w0/(w(w0+w)) = [0.1271 0.127  0.0775 0.0743 0.0718 0.0696]
  eV=2.0 meV: pre*sum w_m k_m = 0.01078 ; kernel k_m = [0.2589 0.2587 0.1565 0.1501 0.1448 0.1403] ; 2eV*w0/(w(w0+w)) = [0.2542 0.254  0.155  0.1487 0.1435 0.1391]
  max over the eVs sweep of the exact full sum minus exact sum without those lines:
  max|analytic high-line contribution| = 0.01078 at eV=-2.00 meV (peak of exact term 0.11099)

(2)-(5)
  example es +-3 meV, 800          ED/ED    potential max|dc-exact| = 0.0110 (9.9% of peak)  pointwise max|dc-exact-analytic_high| = 0.0067
                                            second order max|dc-exact| = 2.60e-01 (peak 3.1965)   sum-rule deficit 0.75-int S = 4.027e-01 (rel 5.37e-01)
  example es +-3 meV, 800          DMRG/KPM potential max|dc-exact| = 0.0108 (9.7% of peak)  pointwise max|dc-exact-analytic_high| = 0.0023
                                            second order max|dc-exact| = 1.00e-02 (peak 3.1965)   sum-rule deficit 0.75-int S = 4.013e-01 (rel 5.35e-01)
  +-40 meV, 10667 (same spacing)   ED/ED    potential max|dc-exact| = 0.0067 (6.1% of peak)  pointwise max|dc-exact-analytic_high| = 0.0108
                                            second order max|dc-exact| = 2.53e-01 (peak 3.1965)   sum-rule deficit 0.75-int S = 2.531e-04 (rel 3.37e-04)
  +-40 meV, 10667 (same spacing)   DMRG/KPM potential max|dc-exact| = 0.0023 (2.1% of peak)  pointwise max|dc-exact-analytic_high| = 0.0108
                                            second order max|dc-exact| = 1.36e-02 (peak 3.1965)   sum-rule deficit 0.75-int S = 2.873e-07 (rel 3.83e-07)

  KPM minus exact-Lehmann on the same grid: +-3 meV 0.0090, +-40 meV 0.0090
```

```
delta=2e-05  n_es=10667  max|KPM-exact| = 0.0023 at eV=-0.80 meV (line at 0.77 meV); max err with |eV-0.77meV|>0.3meV: 0.0000  (2.7s)
delta=1e-05  n_es=21334  max|KPM-exact| = 0.0006 at eV=-0.80 meV (line at 0.77 meV); max err with |eV-0.77meV|>0.3meV: 0.0000  (7.2s)
```

**Suggested fix**: widen the example's `es` to the top of the S_k spectrum plus
several delta (+-20 meV suffices), tighten its assert toward the measured 2.1
per cent (0.05 of the peak), correct its comment (the 0.77 meV lowest
transition, the 10 to 16 meV lines) and keep comparing against
`conductance.third_order_potential_dIdV` rather than `submode="ED"`, which on the
widened grid is itself 0.0067 off the exact sum where KPM is 0.0023; fix the same
comment in `kondo_third_order_dmrg_VS_ED`. In the library, say plainly that the
potential term, unlike the second-order one, needs every transition up to the top
of the S_k spectrum and that a shared `es` in `get_kondo_spectrum` must meet the
stricter requirement; giving the potential term its own grid in
`_get_kondo_spectrum_dmrg` is the structural version. A sum-rule check is sound
and free (sum_k S_k^2 = S(S+1) exactly), but it should warn rather than raise, at
about 1e-2: a flat 1e-3 would fire on the repo's own passing test grid (relative
deficit 7.1e-4 with the Lorentzian `submode="ED"`). A sharper criterion is a true
upper bound on the dI/dV error, deficit x the largest |F0(eV-w) - F0(eV+w)| at
the grid edges x 4 pi T0^2 Jrho_s U, valid when both edges of `es` lie beyond
max|eV|; for the example it is 0.072 against a 0.111 peak and fires. Numbers
change in the example only.

### 15. On `itensor_version="python"`, `excitation_energies(k, n>=2)` on the default iterative path can return n genuine eigenvalues of H_eff(k) that are not the n lowest, dropping one copy of an exactly degenerate level: 7.5 per cent off the third value on the critical 2-site Heisenberg cell at D=10, 24 of 36 calls wrong, while `ROADMAP.md` and five documentation files state that this path is not exposed

`bug` &middot; severity **MEDIUM** &middot; CONFIRMED, NARROWED &middot; lens `infinite-chain`

**Status**: FIXED -- `idmrg_excitations._lowest_iterative` keeps its single ARPACK call from the constant start for n=1, unchanged bit for bit, and for n>=2 calls `_lowest_iterative_deflated`, a port of `vx_lanczos_lowest`: one Lanczos run per value (`_deflated_lanczos_run`, full reorthogonalization, stopped on the residual beta*|s_m| < 1e-10*max(1,|val|), capped at `_ITERATIVE_EIG_MAX_ITER=300`), each from its own seeded generic start (`np.random.default_rng(run)`, which leaves the global numpy stream alone) orthogonalized against the vectors already found, deflation by a shift sized from the first run's Ritz range, each value's residual re-measured against the undeflated operator, and a non-ascending value, a non-converged run or a breakdown handing the answer to the dense path. A hand-written Lanczos rather than chained `eigsh(k=1)` because scipy refuses `which="BE"` on a complex operator, so ARPACK cannot give the Ritz range; it also needs fewer applications (203 against 243, 250 against 333 at n=3). On the n_uc=2 critical Heisenberg cell at maxm=10 (dim=300) over four momenta and n=1 to 4 the old call was wrong in 7 of 16 calls and the new path in 0 of 16 (max error 1.9e-15), at 0.42 to 0.77 s against 0.86 to 0.89 s for dense; from n~4 dense is cheaper, which the `_DENSE_EIG_MAX` comment now says. The reviewer's struck sub-claim needs a correction of its own: on this pass's VUMPS state of the Haldane cell (spin-1, maxm=12, dim=288) the old call returned n=3 at k=pi/2 off by 0.55, where the review's state had given the triplet complete, so a clean Haldane run says nothing general; the new path is 0 of 12 off there. VUMPS on the maxm=4 test cell is not reproducible under a fixed `np.random.seed` (e0 differs at 1e-8 between runs), so the tests compare the two solvers on one shared environment. The `excitation_energies`, `spectral_weights` and `dynamical_structure_factor` docstrings (in `idmrg_excitations.py` and `infinitechain.py`) now state the multiplicity and drop the stale ARPACK sentence. Pinned by `tests/test_audit_2026_09_24_pyitensor.py::test_iterative_path_returns_the_degenerate_pair_with_multiplicity` (maxm=4, dim=48, `_DENSE_EIG_MAX` 0 against 1e9, n=2 and 3 to 1e-10, the degenerate pair asserted as a premise), `::test_single_arpack_call_misses_the_copy_on_this_cell` (the pre-fix call rebuilt in-test), `::test_n1_iterative_path_is_the_single_call_bit_for_bit`, `::test_non_ascending_runs_are_refused_and_dense_answers` and `::test_spectral_weights_multiplet_sums_match_dense`. The C++ comment at `chain_session.h:8742` to `8756` is left for a C++ pass, since editing it would make `make -n pybind` report work. NUMBERS CHANGE: `excitation_energies(k, n>=2)`, `spectral_weights(n>=2)` and `dynamical_structure_factor(n>=2)` on `itensor_version="python"`, on a cell whose H_eff(k) has a degenerate level and dim > 256, move to the dense answer; at k=0.37 on the maxm=10 cell [0.289896072 0.291597913 0.313428946] becomes [0.289896072 0.289896072 0.291597913], and the Sx and Sz weight fractions go 0.058 to 0.184014 and 0.339 to 0.041282, both equal to the dense path's. n=1, `excitation_gap` and dim <= 256 are byte-identical.

**Where**: `src/dmrgpy/pyitensor/idmrg_excitations.py::_lowest_iterative`
(line 938, `v0` at `:973`, a single
`eigsh(op, k=n, which="SA", v0=ones, tol=_ITERATIVE_EIG_TOL)` at `:975`),
taken whenever dim = D*D*(d_g-1) > `_DENSE_EIG_MAX` = 256, which on an
`n_uc=2` spin-1/2 cell means `maxm>=10`; consumers `spectral_weights` and
`dynamical_structure_factor` at n>=2. The statements it contradicts:
`ROADMAP.md:386` to `391`, `docs/user_guide.md:3329`, `docs/user_guide.tex:4160`
to `4161`, `docs/documentation.md:1394`, `docs/documentation.tex:1858`,
`docs/idmrg_improvement_plan.md:100` to `101`, and the `vx_deterministic_start`
comment at `mpscpp3/chain_session.h:8742` to `8756`.

One Krylov space holds only one direction of a degenerate eigenspace in exact
arithmetic, so a single ARPACK call for n pairs can miss the second copy, return
the next distinct level in its place and shift everything after it up one slot;
each returned value is a genuine eigenpair, so the residual check passes it,
nothing falls back to dense and nothing warns. `n=1` and `excitation_gap` are
never affected, which is why every cross-backend dispersion test passes; every
v3-versus-`"python"` test runs at D<=4 (dim<=48), below the dense threshold; and
the only test that forces the iterative path
(`test_excitation_iterative_eigensolver_matches_dense`) is a D=2 TFIM with dim=4,
with no degeneracy to lose. The "not exposed" statement was measured at D=2,
dim=12, where scipy's `ncv = min(max(2n+1,20),dim)` equals dim and the Krylov
basis is the whole space, so the measurement could not discriminate. The defect
came in with `dc2d1db`; the "not exposed" claim with `765b537`.

**Expected**: the n lowest eigenvalues with multiplicity, which dense and v3
return, and which `spectral_weights`' docstring relies on when it tells the user
to group branches by multiplet size and to raise n when the lowest branch carries
no weight.

Repro, from the hunter, on the default path with nothing monkeypatched on the
`"python"` side except to build the dense reference:

```python
# The DEFAULT user path, no monkeypatching: n_uc=2 Heisenberg at D=10 gives
# dim = D*D*(d_g-1) = 300 > _DENSE_EIG_MAX = 256, so itensor_version="python"
# excitation_energies(k, n) goes through ARPACK from one constant start.
# References: the dense path on the same environment, and itensor_version=3
# (deflated Lanczos, checked against its own dense path to 1e-14 above).
import sys, time
import numpy as np
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_excitations as ie

def heis(D, v):
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=v)
    h = (ic.SxC[0]*ic.SxC[1] + ic.SyC[0]*ic.SyC[1] + ic.SzC[0]*ic.SzC[1]
         + ic.SxC[1]*ic.SxR[0] + ic.SyC[1]*ic.SyR[0] + ic.SzC[1]*ic.SzR[0])
    ic.set_hamiltonian(h); ic.gs_method = "vumps"; ic.maxm = D; ic.vumps_nrestarts = 3
    return ic

D = 10
np.random.seed(5)
icp = heis(D, "python"); t = time.time(); e0p = icp.gs_energy()
env = icp._get_excitation_environment()
print("python D=%d dim=%d e0=%.12f conv=%s (%.1fs), _DENSE_EIG_MAX=%d"
      % (D, env.D*env.D*(env.d_g-1), e0p, icp.converged, time.time()-t, ie._DENSE_EIG_MAX)); sys.stdout.flush()
ic3 = heis(D, 3); t = time.time(); e03 = ic3.gs_energy()
print("v3     D=%d e0=%.12f conv=%s (%.1fs)" % (D, e03, ic3.converged, time.time()-t)); sys.stdout.flush()
for k in (0.37, 1.0):
    n = 3
    t = time.time(); dflt = np.asarray(icp.excitation_energies(k, n=n)); t1 = time.time()-t
    ie._DENSE_EIG_MAX = 10**9
    t = time.time(); dense = np.asarray(icp.excitation_energies(k, n=n)); t2 = time.time()-t
    ie._DENSE_EIG_MAX = 256
    v3 = np.asarray(ic3.excitation_energies(k, n=n))
    print("k=%.2f n=%d python default=%s\n              python dense  =%s\n              v3 default    =%s\n   max|default-dense|=%.3e  max|default-v3|=%.3e  (default %.1fs, dense %.1fs)"
          % (k, n, np.array2string(dflt, precision=9), np.array2string(dense, precision=9),
             np.array2string(v3, precision=9), np.max(np.abs(dflt-dense)), np.max(np.abs(dflt-v3)), t1, t2))
    sys.stdout.flush()
```

```
python D=10 dim=300 e0=-0.443041225115 conv=True (22.0s), _DENSE_EIG_MAX=256
v3     D=10 e0=-0.443041225115 conv=True (12.9s)
k=0.37 n=3 python default=[0.289896072 0.291597913 0.313428946]
              python dense  =[0.289896072 0.289896072 0.291597913]
              v3 default    =[0.289896072 0.289896072 0.291597913]
   max|default-dense|=2.183e-02  max|default-v3|=2.183e-02  (default 0.3s, dense 0.9s)
k=1.00 n=3 python default=[0.757676581 0.759652536 0.766500772]
              python dense  =[0.757676581 0.757676581 0.759652536]
              v3 default    =[0.757676581 0.757676581 0.759652536]
   max|default-dense|=6.848e-03  max|default-v3|=6.848e-03  (default 0.6s, dense 0.9s)
```

and the downstream consumer (the reviewer's rerun of the hunter's script, whose
own figures were in no saved output):

```python
# Downstream consumer of the same solver: spectral_weights(n>1) on the default
# path (D=10, dim=300 > 256) against the dense path on the same environment.
import sys, time
import numpy as np
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_excitations as ie

ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
h = (ic.SxC[0]*ic.SxC[1] + ic.SyC[0]*ic.SyC[1] + ic.SzC[0]*ic.SzC[1]
     + ic.SxC[1]*ic.SxR[0] + ic.SyC[1]*ic.SyR[0] + ic.SzC[1]*ic.SzR[0])
ic.set_hamiltonian(h); ic.gs_method = "vumps"; ic.maxm = 10; ic.vumps_nrestarts = 3
np.random.seed(5)
e0 = ic.gs_energy(); print("e0=%.12f conv=%s" % (e0, ic.converged))
k = 0.37
for op in ("Sx", "Sz"):
    e1, w1, tot1 = ic.spectral_weights(op, k, p=0, n=3, return_total=True)
    ie._DENSE_EIG_MAX = 10**9
    e2, w2, tot2 = ic.spectral_weights(op, k, p=0, n=3, return_total=True)
    ie._DENSE_EIG_MAX = 256
    print("%s k=%.2f default: E=%s w=%s sum/total=%.6f" % (op, k, np.array2string(np.asarray(e1), precision=9),
          np.array2string(np.asarray(w1), precision=6), np.sum(w1)/tot1))
    print("%s k=%.2f dense  : E=%s w=%s sum/total=%.6f" % (op, k, np.array2string(np.asarray(e2), precision=9),
          np.array2string(np.asarray(w2), precision=6), np.sum(w2)/tot2))
```

```
e0=-0.443041225116 conv=True
Sx k=0.37 default: E=[0.289896072 0.291597913 0.313428946] w=[0.015624 0.00047  0.004963] sum/total=0.057844
Sx k=0.37 dense  : E=[0.289896072 0.289896072 0.291597913] w=[0.057133 0.007179 0.00047 ] sum/total=0.177950
Sz k=0.37 default: E=[0.289896072 0.291597913 0.313428946] w=[0.005665 0.010284 0.108533] sum/total=0.339299
Sz k=0.37 dense  : E=[0.289896072 0.289896072 0.291597913] w=[0.002898 0.003946 0.010284] sum/total=0.046687
```

**Reviewer (CONFIRMED, NARROWED)**: reproduced to every printed digit on the
hunter's cell and on three VUMPS seeds of the reviewer's own, all giving the same
state: n=2 and n=3 are wrong at every momentum tried (k = 0.37, 1.0, pi/2, 2.5),
by 1.7e-3 to 2.2e-2, and n=4 is right at all four. The anchor is real: dense and
v3 agree on [0.289896072 0.289896072 0.291597913] to every digit, and the pair is
the two transverse members of the lowest magnon triplet of the finite-D state
(SU(2) broken to a U(1) about an axis in the xz-plane), exactly degenerate by the
residual symmetry, with the longitudinal member split off by 1.7e-3 of finite-D
error. Struck: "whenever the iterative path runs" and "every degenerate
eigenvalue", since the Haldane chain (spin-1, `n_uc=1`, `maxm=12`, dim=288, the
default path) returns its exactly degenerate triplet complete at k = pi, 0.8pi and
pi/2 for n=2, 3, 4; the reviewer's reading, not a measurement, is that the miss
appears when the first distinct level just outside the requested set lies close
to it. Struck: "only ncv=dim recovers the copy"; at dim=48 `ncv` of 23, 24, 25,
40 and 48 recover it while 20, 26 and 30 do not, and at dim=300 `tol=0` recovers
it at `ncv=20`, both by roundoff surfacing the copy. Replaced: the hunter's Sx
fraction 0.182 -> 0.042 by the reviewer's own rerun, 0.178 -> 0.058, with the Sz
fraction moving the other way, 0.047 -> 0.339, since the substituted branch
carries most of the Sz weight; what the record should say is that the multiplet
sum cannot be formed from one member. Half right: "what matters in the C++ fix is
the deflation"; deflation with a reused start fails exactly like today's code,
and what matters is deflation plus a fresh start per run, which
`vx_lanczos_lowest`'s own comment at `chain_session.h:8875` already says; the
misleading text is the `vx_deterministic_start` comment at `:8742` to `8756`. The
reviewer's probes:

```python
# Reviewer probe: how often the DEFAULT path (no monkeypatching) misses on the
# n_uc=2 Heisenberg cell at maxm=10 (dim=300), own seeds 1,2,3, a k-scan and
# n=2,3,4.  Reference: dense path on the same environment.
import sys, time
import numpy as np
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_excitations as ie

tally = [0, 0]
for seed in (1, 2, 3):
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
    h = (ic.SxC[0]*ic.SxC[1] + ic.SyC[0]*ic.SyC[1] + ic.SzC[0]*ic.SzC[1]
         + ic.SxC[1]*ic.SxR[0] + ic.SyC[1]*ic.SyR[0] + ic.SzC[1]*ic.SzR[0])
    ic.set_hamiltonian(h); ic.gs_method = "vumps"; ic.maxm = 10; ic.vumps_nrestarts = 3
    np.random.seed(seed)
    t = time.time(); e0 = ic.gs_energy(); env = ic._get_excitation_environment()
    print("seed=%d e0=%.12f conv=%s dim=%d (%.1fs)" % (seed, e0, ic.converged,
          env.D*env.D*(env.d_g-1), time.time()-t)); sys.stdout.flush()
    for k in (0.37, 1.0, np.pi/2, 2.5):
        ie._DENSE_EIG_MAX = 10**9
        ref = np.asarray(ie.excitation_energies(env, k, n=4))
        ie._DENSE_EIG_MAX = 256
        line = "  k=%.4f dense4=%s |" % (k, np.array2string(ref, precision=9))
        for n in (2, 3, 4):
            w = np.asarray(ie.excitation_energies(env, k, n=n))
            err = np.max(np.abs(w - ref[:n]))
            tally[1] += 1
            if err > 1e-8:
                tally[0] += 1
            line += " n=%d err=%.1e" % (n, err)
        print(line); sys.stdout.flush()
print("default-path calls off by >1e-8: %d of %d" % tuple(tally))
```

```
seed=1 e0=-0.443041225115 conv=True dim=300 (19.3s)
  k=0.3700 dense4=[0.289896072 0.289896072 0.291597913 0.313428946] | n=2 err=1.7e-03 n=3 err=2.2e-02 n=4 err=2.2e-15
  k=1.0000 dense4=[0.757676581 0.757676581 0.759652536 0.766500772] | n=2 err=2.0e-03 n=3 err=6.8e-03 n=4 err=3.8e-15
  k=1.5708 dense4=[1.119681275 1.119681275 1.125552763 1.131015047] | n=2 err=5.9e-03 n=3 err=5.9e-03 n=4 err=1.3e-15
  k=2.5000 dense4=[1.503980148 1.503980148 1.512774252 1.518400439] | n=2 err=8.8e-03 n=3 err=8.8e-03 n=4 err=8.9e-16
seed=2 e0=-0.443041225115 conv=True dim=300 (17.9s)
  k=0.3700 dense4=[0.289896072 0.289896072 0.291597913 0.313428946] | n=2 err=1.7e-03 n=3 err=2.2e-02 n=4 err=4.4e-15
  k=1.0000 dense4=[0.757676581 0.757676581 0.759652536 0.766500772] | n=2 err=2.0e-03 n=3 err=6.8e-03 n=4 err=2.9e-15
  k=1.5708 dense4=[1.119681275 1.119681275 1.125552763 1.131015047] | n=2 err=5.9e-03 n=3 err=5.9e-03 n=4 err=1.8e-15
  k=2.5000 dense4=[1.503980148 1.503980148 1.512774252 1.518400439] | n=2 err=8.8e-03 n=3 err=8.8e-03 n=4 err=8.9e-16
seed=3 e0=-0.443041225115 conv=True dim=300 (17.4s)
  k=0.3700 dense4=[0.289896072 0.289896072 0.291597913 0.313428946] | n=2 err=1.7e-03 n=3 err=2.2e-02 n=4 err=5.1e-15
  k=1.0000 dense4=[0.757676581 0.757676581 0.759652536 0.766500772] | n=2 err=2.0e-03 n=3 err=6.8e-03 n=4 err=3.2e-15
  k=1.5708 dense4=[1.119681275 1.119681275 1.125552763 1.131015047] | n=2 err=5.9e-03 n=3 err=5.9e-03 n=4 err=2.2e-15
  k=2.5000 dense4=[1.503980148 1.503980148 1.512774252 1.518400439] | n=2 err=8.8e-03 n=3 err=8.8e-03 n=4 err=4.4e-16
default-path calls off by >1e-8: 24 of 36
```

```python
# Reviewer probe, independent model: spin-1 Heisenberg (Haldane) chain on a
# ONE-site cell, so no zone folding at all; the lowest branch near k=pi is the
# SU(2) magnon triplet.  D=12 -> dim = D*D*(d_g-1) = 288 > _DENSE_EIG_MAX=256,
# i.e. the default excitation_energies path is ARPACK.  Own seed (17).
# References: dense path on the SAME environment, and itensor_version=3.
import sys, time
import numpy as np
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_excitations as ie

def haldane(D, v):
    ic = infinitechain.Infinite_Spin_Chain(["1"], itensor_version=v)
    h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]
    ic.set_hamiltonian(h); ic.gs_method = "vumps"; ic.maxm = D; ic.vumps_nrestarts = 3
    return ic

D = 12
np.random.seed(17)
icp = haldane(D, "python"); t = time.time(); e0p = icp.gs_energy()
env = icp._get_excitation_environment()
dim = env.D*env.D*(env.d_g-1)
print("python D=%d (env.D=%d) dim=%d e0=%.12f conv=%s (%.1fs) default path: %s"
      % (D, env.D, dim, e0p, icp.converged, time.time()-t,
         "iterative" if dim > ie._DENSE_EIG_MAX else "dense")); sys.stdout.flush()
ic3 = haldane(D, 3); t = time.time(); e03 = ic3.gs_energy()
print("v3     D=%d e0=%.12f conv=%s (%.1fs)" % (D, e03, ic3.converged, time.time()-t)); sys.stdout.flush()
for k in (np.pi, 0.8*np.pi, 0.5*np.pi):
    for n in (2, 3, 4):
        t = time.time(); dflt = np.asarray(icp.excitation_energies(k, n=n)); t1 = time.time()-t
        ie._DENSE_EIG_MAX = 10**9
        t = time.time(); dense = np.asarray(icp.excitation_energies(k, n=n)); t2 = time.time()-t
        ie._DENSE_EIG_MAX = 256
        v3 = np.asarray(ic3.excitation_energies(k, n=n))
        print("k=%.4f n=%d default=%s\n                dense  =%s\n                v3     =%s\n   max|default-dense|=%.3e max|dense-v3|=%.3e (default %.1fs, dense %.1fs)"
              % (k, n, np.array2string(dflt, precision=9), np.array2string(dense, precision=9),
                 np.array2string(v3, precision=9), np.max(np.abs(dflt-dense)),
                 np.max(np.abs(dense-v3)), t1, t2))
        sys.stdout.flush()
# spectral weights at k=pi, n=3: the triplet sum the docstring tells the user to form
for op in ("Sx", "Sz"):
    e1, w1, tot1 = icp.spectral_weights(op, np.pi, p=0, n=3, return_total=True)
    ie._DENSE_EIG_MAX = 10**9
    e2, w2, tot2 = icp.spectral_weights(op, np.pi, p=0, n=3, return_total=True)
    ie._DENSE_EIG_MAX = 256
    print("%s k=pi n=3 default: E=%s w=%s sum/total=%.6f" % (op, np.array2string(np.asarray(e1), precision=9),
          np.array2string(np.asarray(w1), precision=6), np.sum(w1)/tot1))
    print("%s k=pi n=3 dense  : E=%s w=%s sum/total=%.6f" % (op, np.array2string(np.asarray(e2), precision=9),
          np.array2string(np.asarray(w2), precision=6), np.sum(w2)/tot2))
```

```
python D=12 (env.D=12) dim=288 e0=-1.401380643518 conv=True (5.9s) default path: iterative
v3     D=12 e0=-1.401380643519 conv=True (3.6s)
k=3.1416 n=2 default=[0.409378474 0.409378474]
                dense  =[0.409378474 0.409378474]
                v3     =[0.409378474 0.409378474]
   max|default-dense|=1.422e-12 max|dense-v3|=2.618e-12 (default 0.2s, dense 1.3s)
k=3.1416 n=3 default=[0.409378474 0.409378474 0.409378474]
                dense  =[0.409378474 0.409378474 0.409378474]
                v3     =[0.409378474 0.409378474 0.409378474]
   max|default-dense|=3.553e-15 max|dense-v3|=2.618e-12 (default 0.4s, dense 1.2s)
k=3.1416 n=4 default=[0.409378474 0.409378474 0.409378474 2.348095922]
                dense  =[0.409378474 0.409378474 0.409378474 2.348095922]
                v3     =[0.409378474 0.409378474 0.409378474 2.348095922]
   max|default-dense|=5.329e-15 max|dense-v3|=6.362e-11 (default 0.4s, dense 1.2s)
k=2.5133 n=2 default=[1.528811029 1.528811029]
                dense  =[1.528811029 1.528811029]
                v3     =[1.528811029 1.528811029]
   max|default-dense|=7.705e-14 max|dense-v3|=1.092e-13 (default 0.3s, dense 1.2s)
k=2.5133 n=3 default=[1.528811029 1.528811029 1.528811029]
                dense  =[1.528811029 1.528811029 1.528811029]
                v3     =[1.528811029 1.528811029 1.528811029]
   max|default-dense|=8.882e-16 max|dense-v3|=1.092e-13 (default 0.5s, dense 1.3s)
k=2.5133 n=4 default=[1.528811029 1.528811029 1.528811029 2.649363146]
                dense  =[1.528811029 1.528811029 1.528811029 2.649363146]
                v3     =[1.528811029 1.528811029 1.528811029 2.649363146]
   max|default-dense|=8.882e-16 max|dense-v3|=4.255e-11 (default 0.5s, dense 1.2s)
k=1.5708 n=2 default=[2.71930372 2.71930372]
                dense  =[2.71930372 2.71930372]
                v3     =[2.71930372 2.71930372]
   max|default-dense|=7.483e-13 max|dense-v3|=1.576e-12 (default 0.4s, dense 1.2s)
k=1.5708 n=3 default=[2.71930372 2.71930372 2.71930372]
                dense  =[2.71930372 2.71930372 2.71930372]
                v3     =[2.71930372 2.71930372 2.71930372]
   max|default-dense|=3.109e-15 max|dense-v3|=1.576e-12 (default 0.7s, dense 1.2s)
k=1.5708 n=4 default=[2.71930372  2.71930372  2.71930372  3.267459703]
                dense  =[2.71930372  2.71930372  2.71930372  3.267459703]
                v3     =[2.71930372  2.71930372  2.71930372  3.267459703]
   max|default-dense|=2.665e-15 max|dense-v3|=5.707e-12 (default 0.6s, dense 1.2s)
Sx k=pi n=3 default: E=[0.409378474 0.409378474 0.409378474] w=[0.058604 2.912973 0.641483] sum/total=0.974312
Sx k=pi n=3 dense  : E=[0.409378474 0.409378474 0.409378474] w=[0.05901  2.914831 0.639242] sum/total=0.974318
Sz k=pi n=3 default: E=[0.409378474 0.409378474 0.409378474] w=[0.074196 0.686521 2.852941] sum/total=0.974473
Sz k=pi n=3 dense  : E=[0.409378474 0.409378474 0.409378474] w=[0.073754 0.684224 2.855104] sum/total=0.974318
```

```python
# Reviewer probe: shape of the fix.  On the hunter's dim=300 environment
# (n_uc=2 Heisenberg, VUMPS maxm=10, np.random.seed(5)), k=0.37 and k=1.0,
# lowest 3 by sequential eigsh(k=1) on H + shift*sum_found |u><u|:
#   (A) the SAME generic v0 for every run   (B) a FRESH generic v0 per run
#   (C) the constant v0 for every run        (D) default path, n=3, tol=0 patched
# Reference: dense path on the same environment.
import sys, time
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_excitations as ie

ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version="python")
h = (ic.SxC[0]*ic.SxC[1] + ic.SyC[0]*ic.SyC[1] + ic.SzC[0]*ic.SzC[1]
     + ic.SxC[1]*ic.SxR[0] + ic.SyC[1]*ic.SyR[0] + ic.SzC[1]*ic.SzR[0])
ic.set_hamiltonian(h); ic.gs_method = "vumps"; ic.maxm = 10; ic.vumps_nrestarts = 3
np.random.seed(5)
e0 = ic.gs_energy(); env = ic._get_excitation_environment()
D, dg = env.D, env.d_g; Dx = D*(dg-1); dim = Dx*D
print("e0=%.12f dim=%d" % (e0, dim)); sys.stdout.flush()

def seq(k, starts, nev=3, shift=10.0):
    mv = lambda x: ie._h_eff_action(k, x.reshape(Dx, D), env).reshape(-1)
    found, vals, napp = [], [], [0]
    for j in range(nev):
        def mvd(x, found=list(found)):
            napp[0] += 1
            y = mv(x)
            for u in found:
                y = y + shift*u*np.vdot(u, x)
            return y
        opd = LinearOperator((dim, dim), dtype=complex, matvec=mvd)
        w1, v1 = eigsh(opd, k=1, which="SA", v0=starts(j), tol=1e-10)
        u = v1[:, 0]
        for f in found:
            u = u - f*np.vdot(f, u)
        found.append(u/np.linalg.norm(u)); vals.append(w1[0])
    return np.asarray(vals) - env.lam_AC, napp[0]

rng = np.random.default_rng(7)
fresh = [rng.standard_normal(dim) + 1j*rng.standard_normal(dim) for _ in range(3)]
ones = np.ones(dim, dtype=complex)/np.sqrt(dim)
for k in (0.37, 1.0):
    ie._DENSE_EIG_MAX = 10**9
    ref = np.asarray(ie.excitation_energies(env, k, n=3))
    ie._DENSE_EIG_MAX = 256
    print("k=%.2f dense: %s" % (k, np.array2string(ref, precision=9)))
    for tag, starts in (("(A) same generic v0 ", lambda j: fresh[0]),
                        ("(B) fresh generic v0", lambda j: fresh[j]),
                        ("(C) constant v0     ", lambda j: ones)):
        t = time.time(); w, na = seq(k, starts)
        print("  %s -> %s max|err|=%.2e (%d applications, %.1fs)"
              % (tag, np.array2string(w, precision=9), np.max(np.abs(w-ref)), na, time.time()-t))
    saved = ie._ITERATIVE_EIG_TOL
    ie._ITERATIVE_EIG_TOL = 0.0
    t = time.time(); w = np.asarray(ie.excitation_energies(env, k, n=3))
    ie._ITERATIVE_EIG_TOL = saved
    print("  (D) default path, _ITERATIVE_EIG_TOL=0 -> %s max|err|=%.2e (%.1fs)"
          % (np.array2string(w, precision=9), np.max(np.abs(w-ref)), time.time()-t))
    sys.stdout.flush()
```

```
e0=-0.443041225115 dim=300
k=0.37 dense: [0.289896072 0.289896072 0.291597913]
  (A) same generic v0  -> [0.289896072 0.291597913 0.313428946] max|err|=2.18e-02 (223 applications, 0.7s)
  (B) fresh generic v0 -> [0.289896072 0.289896072 0.291597913] max|err|=3.77e-15 (243 applications, 0.7s)
  (C) constant v0      -> [0.289896072 0.291597913 0.313428946] max|err|=2.18e-02 (223 applications, 0.7s)
  (D) default path, _ITERATIVE_EIG_TOL=0 -> [0.289896072 0.289896072 0.291597913] max|err|=1.11e-15 (0.9s)
k=1.00 dense: [0.757676581 0.757676581 0.759652536]
  (A) same generic v0  -> [0.757676581 0.759652536 0.766500772] max|err|=6.85e-03 (313 applications, 0.9s)
  (B) fresh generic v0 -> [0.757676581 0.757676581 0.759652536] max|err|=1.22e-15 (333 applications, 1.0s)
  (C) constant v0      -> [0.757676581 0.759652536 0.766500772] max|err|=6.85e-03 (313 applications, 0.9s)
  (D) default path, _ITERATIVE_EIG_TOL=0 -> [0.757676581 0.757676581 0.759652536] max|err|=5.77e-15 (1.3s)
```

**Suggested fix**: port `vx_lanczos_lowest`'s shape to `_lowest_iterative`, one
solve per eigenvalue deflated against the vectors already found, each from its
own fresh generic start orthogonalized against them, refusing a non-ascending
result, with the deflation shift sized from the first run's Ritz range rather
than hardcoded; case (B) above matches dense to 3.8e-15 and 1.2e-15 at 0.7 to 1.0
s, against 0.3 to 0.6 s for today's single call and 0.9 to 1.9 s for dense.
Setting `_ITERATIVE_EIG_TOL=0` also recovers the copy here but only by letting
roundoff surface it, and raising `ncv` is non-monotonic, so neither is a fix;
routing n>1 to dense up to a few hundred is a sound stopgap. Correct the six
documentation passages, the `vx_deterministic_start` comment (a C++ comment edit
makes `make -n pybind` report work, so either rebuild after it or leave it to a
C++ pass) and the ARPACK half of `spectral_weights`' docstring. Regression: n=2,
3 against dense at dim>20, the D=4 cell with `_DENSE_EIG_MAX=0` (dim=48,
reproduces in about 0.1 s). NUMBERS CHANGE: n>=2 on a cell with a degenerate
H_eff(k) level and dim>256, towards the dense answer; n=1 and dim<=256
byte-identical.

### 16. Under `set_pad_bonds(K)`, `tevol_method="TDVP_GSE"` on `itensor_version="python"` silently runs a different algorithm: at the recommended K=maxm the Krylov expansion adds nothing, and the bond growth comes from a QR completion of the padded zero directions instead, which `svd.py` says contribute nothing

`bug` &middot; severity **LOW** &middot; CONFIRMED, NARROWED &middot; lens `recent-misc`

**Status**: FIXED for `quench_tdvp_gse`/`evolve_and_measure_tdvp_gse`; PARTIAL by design for `tdz.py`'s `TDVP_GSE` branch at `tdvp_gse_sweeps=0`. The one-site route is exempt from padding the way the MPO is: `Chain.global_subspace_expand` and `Chain.tdvp_step(num_center=1)` run under `backend.pad_bonds_suspended()`, and `quench_tdvp_gse`/`evolve_and_measure_tdvp_gse` call a new `_strip_bond_padding` once at trajectory entry (gated on `pad_bonds()`, a lossless SVD sweep with padding suspended, on the evolved state only, so the caller's `wf` keeps its padding), after which every expansion and one-site step runs suspended (`Chain._tdvp_onesite_step`). The reviewer's fix 1 is in as belt and braces: `gse._gse_bond_step` reads m as the true rank from the spectrum, a no-op unpadded by construction that keeps a direct padded caller of `gse.py` right. The hunter's second fix is not implemented. `svd.py`'s comment and `set_pad_bonds`' docstring now say the padding is exact for the state and not inert under `qr_split`, what the exemption does and its cost under jit. The residual is deliberate: `Chain.tdvp_step` never strips, because TDZ carries the wavefunction between calls and a per-call strip would delete the expansion's zero-weight directions, so a TDZ run at `tdvp_gse_sweeps=0` under padding still starts from the padded state (0.467 against unpadded in an emulated TDZ loop, Neel start, n=8, K=4); at `sweeps=3` its first expansion strips it (to 8e-13). Measured on the reviewer's quench (XXZ Delta=0.7, hz=0.1, Neel start, 40 steps of dt=0.05), padded against unpadded: n=10, K=maxm=4, sweeps=0 from 0.4928 to 8.9e-16; n=12, K=8, sweeps=0 from 0.4929 to 1.3e-14; n=10, K=4, sweeps=3 from 4.086e-6 to 1.70e-7; n=12, K=8, sweeps=3 from 1.911e-7 to 4.7e-10; a `quench_tdvp_gse` correlator at n=8, K=4 from 0.34 to 1.2e-15; and on the Gram SVD route (K=24, above `_GRAM_MIN_DIM`) 2.9e-15 and 1.3e-10 from a Neel start and 1.2e-12 and 1.1e-12 from an entangled start. Unpadded runs are unchanged, `array_equal` at sweeps 0 and 3. Pinned by `tests/test_audit_2026_09_24_pyitensor.py::test_padded_one_site_tdvp_follows_the_unpadded_trajectory`, `::test_padded_tdvp_gse_follows_the_unpadded_trajectory`, `::test_unpadded_runs_are_unchanged_by_the_exemption`, `::test_padded_quench_tdvp_gse_follows_the_unpadded_correlator` and `::test_krylov_expansion_sees_the_true_bond_dimension_under_padding`, with the pre-fix route rebuilt in-test by monkeypatching. NUMBERS CHANGE: every padded `TDVP_GSE` run through `quench_tdvp_gse`/`evolve_and_measure_tdvp_gse` moves onto the unpadded one; note the direction at `sweeps=0`, where the padded run used to be 9.8e-5 from ED only because it was one-site TDVP on the larger D=K manifold and is now 0.4929 from ED, like the unpadded frozen product state, which is the correct one-site answer. Under JAX with `set_jit("auto")` this route again retraces once per bond dimension it grows through, the cost padding was meant to remove; it never kept frozen shapes anyway.

**Where**: `src/dmrgpy/pyitensor/gse.py::_gse_bond_step` (`m = bond_v.dim`, the
padded dimension, so `room = bond_maxdim - m = 0` at K = maxm);
`src/dmrgpy/pyitensor/tdvp.py`'s `qr_split` (reduced QR,
`k = min(ldim, rdim)`), which returns an orthonormal completion for padded zero
directions; the statements `svd.py` makes (the padded columns "contribute
nothing to any contraction") and `backend.set_pad_bonds` makes ("the represented
state is unchanged"); `docs/user_guide.md` section 19, which recommends
`set_pad_bonds(K)` "with K your bond dimension"; `tdz.py:169` to `176`, whose
`TDVP_GSE` branch calls the same session pair (by reading).

The represented state is unchanged at every instant, which is what
`set_pad_bonds` promises, but the algorithm is not. At K = maxm the expansion
adds no direction at any bond of any call (27 of 27 bond steps at n=10, K=4 and
33 of 33 at n=12, K=8, against 1 to 4 per bond unpadded), and in the first
left-to-right one-site half-sweep of the first step `qr_split` completes the
padded zero directions into live zero-weight basis vectors that one-site TDVP
then populates, so bond growth comes from an arbitrary QR completion rather than
from the Krylov subspace. You can think of the padded route as one-site TDVP on
the D=K manifold: a larger ansatz, not a wrong integrator.
`docs/gpu_cpu_performance.md` already says padding is "exact in representation,
not in trajectory", but it names a different mechanism (perturbed truncation
decisions). No test runs a one-site method under padding.

**Expected**: padded and unpadded runs of the same method agreeing to roundoff,
as they do for two-site TDVP, TEBD, energies, excited states, `vev`, correlators,
entropy, overlaps and KPM (1e-14 to 1e-15, measured).

Repro, from the hunter:

```python
# One-site TDVP with no expansion at all (tdvp_gse_sweeps=0) from a Neel
# product state: does padding alone change the trajectory, and does the
# centre-bond entanglement grow under a method that conserves bond
# dimension by construction? 8-site Heisenberg quench, <Sz_0>(t) vs ED.
import numpy as np
from dmrgpy import spinchain, timedependent
from dmrgpy.pyitensor import backend as bk

n = 8
def quench(pad, sweeps, mode="DMRG", nt=40, dt=0.05):
    bk.set_pad_bonds(pad)
    sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
    sc.tevol_method = "TDVP_GSE"
    sc.tdvp_gse_sweeps = sweeps
    sc.maxm = 16; sc.nsweeps = 8
    h0 = 0
    for i in range(n): h0 = h0 + (-1)**i*sc.Sz[i]
    h1 = 0
    for i in range(n-1):
        h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h0)
    wf = sc.get_gs(mode=mode)
    sc.set_hamiltonian(h1)
    kw = dict(return_wf=True) if mode == "DMRG" else {}
    out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=nt, dt=dt,
                                           wf=wf, mode=mode, **kw)
    sz = out[1]
    ent = out[2].get_bond_entropy(n//2 - 1) if mode == "DMRG" else float("nan")
    bk.set_pad_bonds(None)
    return np.real(np.asarray(sz)), ent

ref, _ = quench(None, 0, mode="ED")
print("ED: <Sz_0>(0) = %.4f, <Sz_0>(t=1.95) = %.4f" % (ref[0], ref[-1]))
for sweeps in (0, 3):
    for pad in (None, 16):
        np.random.seed(7)
        y, ent = quench(pad, sweeps)
        print("tdvp_gse_sweeps=%d pad=%-4s max|DMRG-ED| = %.3e   <Sz_0>(t=1.95) = %.4f   final centre-bond entropy = %.6f"
              % (sweeps, pad, np.max(np.abs(y-ref)), y[-1], ent))
```

```
ED: <Sz_0>(0) = -0.5000, <Sz_0>(t=1.95) = -0.0216
tdvp_gse_sweeps=0 pad=None max|DMRG-ED| = 4.784e-01   <Sz_0>(t=1.95) = -0.5000   final centre-bond entropy = 0.000000
tdvp_gse_sweeps=0 pad=16   max|DMRG-ED| = 2.438e-08   <Sz_0>(t=1.95) = -0.0216   final centre-bond entropy = 1.055603
tdvp_gse_sweeps=3 pad=None max|DMRG-ED| = 1.789e-07   <Sz_0>(t=1.95) = -0.0216   final centre-bond entropy = 1.055594
tdvp_gse_sweeps=3 pad=16   max|DMRG-ED| = 2.438e-08   <Sz_0>(t=1.95) = -0.0216   final centre-bond entropy = 1.055603
```

```python
# qr_split on a tensor whose right bond was padded with zero columns:
# are the extra columns of Q zero (dead) or orthonormal (live)?
import numpy as np
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor
from dmrgpy.pyitensor.svd import qr_split, svd
from dmrgpy.pyitensor import backend as bk
l, s, r = Index(4, tags="Link"), Index(2, tags="Site"), Index(3, tags="Link")
rng = np.random.default_rng(0)
T = ITensor((l, s, r), rng.standard_normal((4, 2, 3)) + 0j)
bk.set_pad_bonds(8)
U, S, V, spec = svd(T, [l, s])          # true rank 3, padded to 8
bk.set_pad_bonds(None)
Umat = U.transpose_to([l, s] + [i for i in U.inds if i not in (l, s)]).reshape(8, -1)
print("svd U (padded): column norms =", np.round(np.linalg.norm(Umat, axis=0), 3))
M = (U*S)                                  # a site tensor carrying the padded bond
Q, C, _ = qr_split(M, [l, s], orthonormal="left")
Qm = Q.transpose_to([l, s] + [i for i in Q.inds if i not in (l, s)]).reshape(8, -1)
print("qr_split Q of the padded tensor: column norms =", np.round(np.linalg.norm(Qm, axis=0), 3),
      " rank of the tensor it factorizes =", np.linalg.matrix_rank(M.transpose_to([l, s] + [i for i in M.inds if i not in (l, s)]).reshape(8, -1)))
```

```
svd U (padded): column norms = [1. 1. 1. 0. 0. 0. 0. 0.]
qr_split Q of the padded tensor: column norms = [1. 1. 1. 1. 1. 1. 1. 1.]  rank of the tensor it factorizes = 3
```

**Reviewer (CONFIRMED, NARROWED)**: the mechanism established causally on the
route itself, not only on the isolated `qr_split` call: replacing the split with
one that keeps zero-weight directions dead puts a padded `tdvp_gse_sweeps=0` run
back onto the unpadded one to 1.8e-15. What does not survive is the framing as
an accuracy defect: against an independent scipy ED, on an XXZ quench of the
reviewer's own at K below the needed rank, the padded run is as close to ED as
the unpadded one or closer in all five configurations measured, and at the
default `tdvp_gse_sweeps=3` padding moves `<Sz_0>(t)` by only 1.9e-7 to 4.1e-6,
within TDVP's own error. Struck as a defect of its own: "at pad=16<maxm GSE adds
on top of the padded 16, defeating padding's shape freezing"; it reproduces with
the same root cause, but K<maxm is not the recommended configuration and the
one-site route never had frozen shapes (`qr_split` un-pads the edge bonds on the
first step). Narrowed: the hunter's "padded to 16 it reaches the exact
trajectory" and "2.8e-6 to 4.5e-9 at n=6" were both at K at least the full
Schmidt rank, where the padded manifold is the whole Hilbert space. Narrowed:
"documented as exact" holds for the state; what is wrong is `svd.py`'s comment,
and what is missing is any statement that one-site TDVP or GSE differ. The
reviewer's probes (its `probe_d.py` execs the head of `probe_a.py`):

```python
# Reviewer probe A: set_pad_bonds on the one-site TDVP route, instrumented
# on the route itself, at K = maxm BELOW the rank the state needs (the
# hunter used K = 2^(n/2), where the padded manifold is the whole space).
# Own chain (n=10 XXZ, Delta=0.7, plus a weak uniform field), own seed,
# independent ED reference built here with numpy/scipy (not dmrgpy).
import sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
from dmrgpy import spinchain, timedependent
from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor import tdvp as tdvpmod
from dmrgpy.pyitensor import chain as chainmod
from dmrgpy.pyitensor.svd import svd as _svd_fn
import os, importlib
if os.environ.get("NOGRAM"):
    # force the exact-SVD branch everywhere, removing the Gram route's
    # shape-dependent switch (min(m,n) >= _GRAM_MIN_DIM) from the comparison
    importlib.import_module("dmrgpy.pyitensor.svd")._GRAM_MIN_DIM = 10**9
    print("NOGRAM: exact SVD everywhere")
from dmrgpy.pyitensor.mpsalgebra import _link_at, inner

n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
K = int(sys.argv[2]) if len(sys.argv) > 2 else 4
DELTA, HZ = 0.7, 0.1
NT, DT = 40, 0.05

# ---------------- independent ED reference ----------------
sx = np.array([[0, .5], [.5, 0]]); sy = np.array([[0, -.5j], [.5j, 0]]); sz = np.diag([.5, -.5])
def op(o, i):
    return sp.kron(sp.kron(sp.identity(2**i), sp.csr_matrix(o)), sp.identity(2**(n-i-1)), format="csr")
H = sum(op(sx, i) @ op(sx, i+1) + op(sy, i) @ op(sy, i+1) + DELTA*op(sz, i) @ op(sz, i+1) for i in range(n-1))
H = H + HZ*sum(op(sz, i) for i in range(n))
# Neel ground state of h0 = sum (-1)^i Sz_i: site 0 down, site 1 up, ...
bits = [1 if i % 2 == 0 else 0 for i in range(n)]   # basis index 1 = down (sz=-1/2)
idx = int("".join(str(b) for b in bits), 2)
psi0 = np.zeros(2**n, complex); psi0[idx] = 1.0
psis = expm_multiply(-1j*H, psi0, start=0.0, stop=DT*(NT-1), num=NT, endpoint=True)
SZ0 = op(sz, 0)
ref = np.array([np.vdot(p, SZ0 @ p).real for p in psis])

# ---------------- instrumentation ----------------
log = {"dims": [], "qr_calls": 0, "qr_completions": 0}
_real_step = chainmod._tdvp_step_fn
def rec_step(psi, *a, **k):
    out = _real_step(psi, *a, **k)
    log["dims"].append([_link_at(out, i, i+1).dim for i in range(1, out.length())])
    return out
chainmod._tdvp_step_fn = rec_step

_real_qr = tdvpmod.qr_split
def rec_qr(T, left_inds, tags="Link", orthonormal="left"):
    A, B, bond = _real_qr(T, left_inds, tags=tags, orthonormal=orthonormal)
    left_inds = list(left_inds)
    right = [i for i in T.inds if i not in left_inds]
    mat = np.asarray(T.transpose_to(left_inds + right)).reshape(
        int(np.prod([i.dim for i in left_inds])), -1)
    s = np.linalg.svd(mat, compute_uv=False)
    r = int(np.sum(s > 1e-13*s[0])) if s[0] > 0 else 0
    log["qr_calls"] += 1
    if bond.dim > r:
        log["qr_completions"] += 1
    return A, B, bond
tdvpmod.qr_split = rec_qr

def dead_split(T, left_inds, tags="Link", orthonormal="left"):
    """Counterfactual: split with the (padded) svd(cutoff=0) instead of QR,
    so directions with zero weight stay zero vectors (dead) rather than
    being completed to an orthonormal basis."""
    left_inds = list(left_inds)
    U, S, V, _ = _svd_fn(T, left_inds, cutoff=0.0, maxdim=None)
    bu = next(i for i in U.inds if i not in left_inds)
    bv = next(i for i in V.inds if i not in T.inds)
    if orthonormal == "left":
        return U, S*V, bu
    return U*S, V, bv


# Hunter's fix #1: in _gse_bond_step take m as phi's true rank (drop the
# padded zero rows of V1) instead of the padded bond dimension.
from dmrgpy.pyitensor import gse as gsemod
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor
from dmrgpy.pyitensor.svd import eigh_truncate
_real_gse_step = gsemod._gse_bond_step
gse_log = []
def gse_step_logged(B_phi, left_link, B_companions, right_inds, cutoff, bond_maxdim=None):
    _, S1, V1, _ = _svd_fn(B_phi, [left_link], cutoff=0.0, maxdim=None)
    bond_v = next(ind for ind in V1.inds if ind not in right_inds)
    combined = int(np.prod([ind.dim for ind in right_inds]))
    V1m = np.asarray(V1.transpose_to([bond_v] + right_inds)).reshape(bond_v.dim, combined)
    m_true = int(np.sum(np.linalg.norm(V1m, axis=1) > 0))
    out = _real_gse_step(B_phi, left_link, B_companions, right_inds, cutoff, bond_maxdim)
    new_dim = next(ind for ind in out[0].inds if ind not in right_inds).dim
    gse_log.append((bond_v.dim, m_true, new_dim - bond_v.dim))
    return out
def gse_step_truerank(B_phi, left_link, B_companions, right_inds, cutoff, bond_maxdim=None):
    """_gse_bond_step with m = number of nonzero rows of V1 (true rank)."""
    _xp = np
    combined = int(np.prod([ind.dim for ind in right_inds])) if right_inds else 1
    _, S1, V1, _ = _svd_fn(B_phi, [left_link], cutoff=0.0, maxdim=None)
    bond_v = next(ind for ind in V1.inds if ind not in right_inds)
    V1_mat = np.asarray(V1.transpose_to([bond_v] + right_inds)).reshape(bond_v.dim, combined)
    V1_mat = V1_mat[np.linalg.norm(V1_mat, axis=1) > 0]
    m = V1_mat.shape[0]
    rho2 = np.zeros((combined, combined), complex)
    Bk_mats, Bk_lefts = [], []
    for Bk in B_companions:
        left_k = next(ind for ind in Bk.inds if ind not in right_inds)
        Bk_mat = np.asarray(Bk.transpose_to([left_k] + right_inds)).reshape(left_k.dim, combined)
        Bk_mats.append(Bk_mat); Bk_lefts.append(left_k)
        rho2 += Bk_mat.conj().T @ Bk_mat
    proj = np.eye(combined) - V1_mat.conj().T @ V1_mat
    rho2_proj = proj @ rho2 @ proj
    rho2_proj = 0.5*(rho2_proj + rho2_proj.conj().T)
    U2 = np.zeros((combined, 0))
    nr = np.linalg.norm(rho2)
    if nr > 0 and np.linalg.norm(rho2_proj)/nr >= 1e-12:
        room = None if bond_maxdim is None else max(0, bond_maxdim - m)
        U2, _ = eigh_truncate(rho2_proj, cutoff, room, mindim=0)
    new_res_mat = np.concatenate([V1_mat, U2.conj().T], axis=0)
    new_dim = new_res_mat.shape[0]
    new_link = Index(new_dim, tags="Link")
    new_res_b = ITensor((new_link,) + tuple(right_inds),
                        new_res_mat.reshape((new_dim,) + tuple(ind.dim for ind in right_inds)))
    dag = new_res_mat.conj().T
    Bphi_mat = np.asarray(B_phi.transpose_to([left_link] + right_inds)).reshape(left_link.dim, combined)
    new_B_phi = ITensor((left_link, new_link), Bphi_mat @ dag)
    new_Bc = [ITensor((Bk_lefts[k], new_link), Bk_mats[k] @ dag) for k in range(len(B_companions))]
    gse_log.append((bond_v.dim, m, new_dim - m))
    return new_res_b, new_B_phi, new_Bc

def dense(psi):
    T = psi.A(1)
    for i in range(2, psi.length()+1):
        T = T * psi.A(i)
    sites = [next(ind for ind in psi.A(i).inds if ind.hastags("Site")) for i in range(1, psi.length()+1)]
    v = np.asarray(T.transpose_to(sites)).reshape(-1)
    return v/np.linalg.norm(v)

def quench(pad, method, sweeps, maxm, seed, split=None, gse_step=None):
    log["dims"].clear(); log["qr_calls"] = 0; log["qr_completions"] = 0; gse_log.clear()
    tdvpmod.qr_split = split if split is not None else rec_qr
    gsemod._gse_bond_step = gse_step if gse_step is not None else gse_step_logged
    np.random.seed(seed)
    bk.set_pad_bonds(pad)
    try:
        sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
        sc.tevol_method = method
        sc.tdvp_gse_sweeps = sweeps
        sc.maxm = maxm; sc.nsweeps = 8
        h0 = 0
        for i in range(n): h0 = h0 + (-1)**i*sc.Sz[i]
        h1 = 0
        for i in range(n-1):
            h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + DELTA*sc.Sz[i]*sc.Sz[i+1]
        for i in range(n): h1 = h1 + HZ*sc.Sz[i]
        sc.set_hamiltonian(h0)
        wf = sc.get_gs()
        wf_dims = [_link_at(wf.cpp_handle, i, i+1).dim for i in range(1, n)]
        v0 = dense(wf.cpp_handle)
        sc.set_hamiltonian(h1)
        out = timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=NT, dt=DT,
                                               wf=wf, return_wf=True)
        y = np.real(np.asarray(out[1]))
        ent = out[2].get_bond_entropy(n//2 - 1)
        return y, v0, wf_dims, ent
    finally:
        bk.set_pad_bonds(None)
        tdvpmod.qr_split = rec_qr
        gsemod._gse_bond_step = gse_step_logged

print("n=%d, K=maxm=%d; independent ED: <Sz_0>(0)=%.4f  <Sz_0>(t=%.2f)=%.6f" % (n, K, ref[0], DT*(NT-1), ref[-1]))
SEED = 11
res = {}
cases = (
    ("GSE0 unpadded               ", None, "TDVP_GSE", 0, None, None),
    ("GSE0 pad=K                  ", K,    "TDVP_GSE", 0, None, None),
    ("GSE0 pad=K dead-split       ", K,    "TDVP_GSE", 0, dead_split, None),
    ("TDVP2 maxm=K unpadded       ", None, "TDVP",     0, None, None),
    ("GSE3 unpadded               ", None, "TDVP_GSE", 3, None, None),
    ("GSE3 pad=K                  ", K,    "TDVP_GSE", 3, None, None),
    ("GSE3 pad=K dead-split       ", K,    "TDVP_GSE", 3, dead_split, None),
    ("GSE3 unpadded dead-split    ", None, "TDVP_GSE", 3, dead_split, None),
    ("GSE3 pad=K truerank-m       ", K,    "TDVP_GSE", 3, None, gse_step_truerank),
    ("GSE3 unpadded truerank-m    ", None, "TDVP_GSE", 3, None, gse_step_truerank),
)
for label, pad, method, sweeps, split, gstep in cases:
    y, v0, wf_dims, ent = quench(pad, method, sweeps, K, SEED, split, gstep)
    res[label] = (y, v0)
    print("%s max|DMRG-ED|=%.3e  <Sz0>(end)=%.6f  S_centre=%.5f" % (label, np.max(np.abs(y-ref)), y[-1], ent))
    print("     links: start %s, after step 1 %s, end %s; qr_split calls %d, completing %d"
          % (wf_dims, log["dims"][0], log["dims"][-1], log["qr_calls"], log["qr_completions"]))
    if gse_log:
        nb = n - 1
        calls = [gse_log[c*nb:(c+1)*nb] for c in range(len(gse_log)//nb)]
        for c, cl in enumerate(calls):
            print("     GSE call %d, per bond (m_used, m_true, added): %s" % (c, cl))
v_ed = psi0
for label in res:
    print("start state %s |<ED Neel|v0>| = %.15f" % (label, abs(np.vdot(v_ed, res[label][1]))))
def d(l1, l2):
    return np.max(np.abs(res[l1][0] - res[l2][0]))
L = [c[0] for c in cases]
print("max|GSE0 pad=K - GSE0 unpadded|             = %.3e" % d(L[1], L[0]))
print("max|GSE0 pad=K dead-split - GSE0 unpadded|  = %.3e" % d(L[2], L[0]))
print("max|GSE0 pad=K - TDVP2 maxm=K|              = %.3e" % d(L[1], L[3]))
print("max|GSE3 pad=K - GSE3 unpadded|             = %.3e" % d(L[5], L[4]))
print("max|GSE3 pad=K dead-split - GSE3 unpadded|  = %.3e" % d(L[6], L[4]))
print("max|GSE3 unpadded dead-split - GSE3 unpad.| = %.3e" % d(L[7], L[4]))
print("max|GSE3 pad=K truerank-m - GSE3 unpadded|  = %.3e" % d(L[8], L[4]))
print("max|GSE3 unpadded truerank-m - GSE3 unpad.| = %.3e" % d(L[9], L[4]))
```

```
n=10, K=maxm=4; independent ED: <Sz_0>(0)=-0.5000  <Sz_0>(t=1.95)=-0.007077
GSE0 unpadded                max|DMRG-ED|=4.929e-01  <Sz0>(end)=-0.500000  S_centre=0.00000
     links: start [1, 1, 1, 1, 1, 1, 1, 1, 1], after step 1 [1, 1, 1, 1, 1, 1, 1, 1, 1], end [1, 1, 1, 1, 1, 1, 1, 1, 1]; qr_split calls 720, completing 0
GSE0 pad=K                   max|DMRG-ED|=9.777e-05  <Sz0>(end)=-0.007175  S_centre=1.00099
     links: start [4, 4, 4, 4, 4, 4, 4, 4, 4], after step 1 [2, 4, 4, 4, 4, 4, 4, 4, 2], end [2, 4, 4, 4, 4, 4, 4, 4, 2]; qr_split calls 720, completing 9
GSE0 pad=K dead-split        max|DMRG-ED|=4.929e-01  <Sz0>(end)=-0.500000  S_centre=0.00000
     links: start [4, 4, 4, 4, 4, 4, 4, 4, 4], after step 1 [4, 4, 4, 4, 4, 4, 4, 4, 4], end [4, 4, 4, 4, 4, 4, 4, 4, 4]; qr_split calls 0, completing 0
TDVP2 maxm=K unpadded        max|DMRG-ED|=7.066e-05  <Sz0>(end)=-0.007148  S_centre=1.00063
     links: start [1, 1, 1, 1, 1, 1, 1, 1, 1], after step 1 [2, 4, 4, 4, 4, 4, 4, 4, 2], end [2, 4, 4, 4, 4, 4, 4, 4, 2]; qr_split calls 0, completing 0
GSE3 unpadded                max|DMRG-ED|=9.298e-05  <Sz0>(end)=-0.007170  S_centre=1.00077
     links: start [1, 1, 1, 1, 1, 1, 1, 1, 1], after step 1 [2, 4, 4, 4, 4, 4, 4, 3, 2], end [2, 4, 4, 4, 4, 4, 4, 4, 2]; qr_split calls 720, completing 1
     GSE call 0, per bond (m_used, m_true, added): [(1, 1, 1), (1, 1, 2), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3)]
     GSE call 1, per bond (m_used, m_true, added): [(2, 2, 0), (3, 3, 1), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (2, 2, 2)]
     GSE call 2, per bond (m_used, m_true, added): [(2, 2, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (2, 2, 2)]
GSE3 pad=K                   max|DMRG-ED|=8.889e-05  <Sz0>(end)=-0.007166  S_centre=1.00099
     links: start [4, 4, 4, 4, 4, 4, 4, 4, 4], after step 1 [2, 4, 4, 4, 4, 4, 4, 4, 2], end [2, 4, 4, 4, 4, 4, 4, 4, 2]; qr_split calls 720, completing 11
     GSE call 0, per bond (m_used, m_true, added): [(4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0)]
     GSE call 1, per bond (m_used, m_true, added): [(4, 2, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 2, 0)]
     GSE call 2, per bond (m_used, m_true, added): [(4, 2, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 2, 0)]
GSE3 pad=K dead-split        max|DMRG-ED|=4.929e-01  <Sz0>(end)=-0.500000  S_centre=0.00000
     links: start [4, 4, 4, 4, 4, 4, 4, 4, 4], after step 1 [4, 4, 4, 4, 4, 4, 4, 4, 4], end [4, 4, 4, 4, 4, 4, 4, 4, 4]; qr_split calls 0, completing 0
     GSE call 0, per bond (m_used, m_true, added): [(4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0)]
     GSE call 1, per bond (m_used, m_true, added): [(4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0)]
     GSE call 2, per bond (m_used, m_true, added): [(4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0), (4, 1, 0)]
GSE3 unpadded dead-split     max|DMRG-ED|=9.046e-05  <Sz0>(end)=-0.007168  S_centre=1.00097
     links: start [1, 1, 1, 1, 1, 1, 1, 1, 1], after step 1 [2, 4, 4, 4, 4, 4, 4, 3, 2], end [2, 4, 4, 4, 4, 4, 4, 4, 2]; qr_split calls 0, completing 0
     GSE call 0, per bond (m_used, m_true, added): [(1, 1, 1), (1, 1, 2), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3)]
     GSE call 1, per bond (m_used, m_true, added): [(2, 2, 0), (3, 3, 1), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (2, 2, 2)]
     GSE call 2, per bond (m_used, m_true, added): [(2, 2, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (2, 2, 2)]
GSE3 pad=K truerank-m        max|DMRG-ED|=9.278e-05  <Sz0>(end)=-0.007170  S_centre=1.00086
     links: start [4, 4, 4, 4, 4, 4, 4, 4, 4], after step 1 [2, 4, 4, 4, 4, 4, 4, 3, 2], end [2, 4, 4, 4, 4, 4, 4, 4, 2]; qr_split calls 720, completing 1
     GSE call 0, per bond (m_used, m_true, added): [(4, 1, 1), (4, 1, 2), (4, 1, 3), (4, 1, 3), (4, 1, 3), (4, 1, 3), (4, 1, 3), (4, 1, 3), (4, 1, 3)]
     GSE call 1, per bond (m_used, m_true, added): [(4, 2, 0), (4, 3, 1), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 2, 2)]
     GSE call 2, per bond (m_used, m_true, added): [(4, 2, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 2, 2)]
GSE3 unpadded truerank-m     max|DMRG-ED|=9.298e-05  <Sz0>(end)=-0.007170  S_centre=1.00077
     links: start [1, 1, 1, 1, 1, 1, 1, 1, 1], after step 1 [2, 4, 4, 4, 4, 4, 4, 3, 2], end [2, 4, 4, 4, 4, 4, 4, 4, 2]; qr_split calls 720, completing 1
     GSE call 0, per bond (m_used, m_true, added): [(1, 1, 1), (1, 1, 2), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3), (1, 1, 3)]
     GSE call 1, per bond (m_used, m_true, added): [(2, 2, 0), (3, 3, 1), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (2, 2, 2)]
     GSE call 2, per bond (m_used, m_true, added): [(2, 2, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (4, 4, 0), (2, 2, 2)]
start state GSE0 unpadded                |<ED Neel|v0>| = 1.000000000000000
start state GSE0 pad=K                   |<ED Neel|v0>| = 1.000000000000000
start state GSE0 pad=K dead-split        |<ED Neel|v0>| = 1.000000000000000
start state TDVP2 maxm=K unpadded        |<ED Neel|v0>| = 1.000000000000000
start state GSE3 unpadded                |<ED Neel|v0>| = 1.000000000000000
start state GSE3 pad=K                   |<ED Neel|v0>| = 1.000000000000000
start state GSE3 pad=K dead-split        |<ED Neel|v0>| = 1.000000000000000
start state GSE3 unpadded dead-split     |<ED Neel|v0>| = 1.000000000000000
start state GSE3 pad=K truerank-m        |<ED Neel|v0>| = 1.000000000000000
start state GSE3 unpadded truerank-m     |<ED Neel|v0>| = 1.000000000000000
max|GSE0 pad=K - GSE0 unpadded|             = 4.928e-01
max|GSE0 pad=K dead-split - GSE0 unpadded|  = 1.832e-15
max|GSE0 pad=K - TDVP2 maxm=K|              = 2.711e-05
max|GSE3 pad=K - GSE3 unpadded|             = 4.086e-06
max|GSE3 pad=K dead-split - GSE3 unpadded|  = 4.928e-01
max|GSE3 unpadded dead-split - GSE3 unpad.| = 2.518e-06
max|GSE3 pad=K truerank-m - GSE3 unpadded|  = 1.981e-07
max|GSE3 unpadded truerank-m - GSE3 unpad.| = 0.000e+00
```

```python
# Candidate alternative fix, tested by monkeypatch: exempt the one-site route
# from padding the way 2d86cff exempted the MPO (backend.pad_bonds_suspended
# around global_subspace_expand and the num_center=1 tdvp_step), and strip
# the padded zeros from the input once, only when no GSE call has touched it
# yet (GSE's own lossless position(n) strips them otherwise). qr_split and
# gse.py are left exactly as they are.
import sys
src = open("probe_a.py").read().split('print("n=%d, K=maxm=%d')[0]
exec(src)
state = {"gse_seen": False}
_real_gse_fn = chainmod._global_subspace_expand_fn
def gse_suspended(*a, **k):
    state["gse_seen"] = True
    with bk.pad_bonds_suspended():
        return _real_gse_fn(*a, **k)
def step_suspended(psi, *a, **k):
    if k.get("num_center", 2) != 1:
        return rec_step(psi, *a, **k)
    with bk.pad_bonds_suspended():
        if not state["gse_seen"]:
            nn = psi.length()
            if psi.center is None:
                psi.center = nn      # _shift_left's SVD is exact from any gauge
            psi.position(nn); psi.position(1)
            state["gse_seen"] = True
        return rec_step(psi, *a, **k)
def run(pad, sweeps, fixed):
    state["gse_seen"] = False
    if fixed:
        chainmod._global_subspace_expand_fn = gse_suspended
        chainmod._tdvp_step_fn = step_suspended
    try:
        return quench(pad, "TDVP_GSE", sweeps, K, 11)
    finally:
        chainmod._global_subspace_expand_fn = _real_gse_fn
        chainmod._tdvp_step_fn = rec_step
print("n=%d, K=maxm=%d" % (n, K))
for sweeps in (0, 3):
    y0 = run(None, sweeps, False)
    yp = run(K, sweeps, False)
    dims_p = list(log["dims"][-1])
    yf = run(K, sweeps, True)
    dims_f = list(log["dims"][-1])
    yu = run(None, sweeps, True)
    print("sweeps=%d: |pad - unpadded| = %.3e   |pad+exempt - unpadded| = %.3e   |unpadded+exempt - unpadded| = %.3e"
          % (sweeps, np.max(np.abs(yp[0]-y0[0])), np.max(np.abs(yf[0]-y0[0])), np.max(np.abs(yu[0]-y0[0]))))
    print("          errors vs ED: unpadded %.3e, pad %.3e, pad+exempt %.3e; final links pad %s, pad+exempt %s"
          % (np.max(np.abs(y0[0]-ref)), np.max(np.abs(yp[0]-ref)), np.max(np.abs(yf[0]-ref)), dims_p, dims_f))
```

```
n=10, K=maxm=4
sweeps=0: |pad - unpadded| = 4.928e-01   |pad+exempt - unpadded| = 1.055e-15   |unpadded+exempt - unpadded| = 8.327e-16
          errors vs ED: unpadded 4.929e-01, pad 9.777e-05, pad+exempt 4.929e-01; final links pad [2, 4, 4, 4, 4, 4, 4, 4, 2], pad+exempt [1, 1, 1, 1, 1, 1, 1, 1, 1]
sweeps=3: |pad - unpadded| = 4.086e-06   |pad+exempt - unpadded| = 8.333e-08   |unpadded+exempt - unpadded| = 0.000e+00
          errors vs ED: unpadded 9.298e-05, pad 8.889e-05, pad+exempt 9.306e-05; final links pad [2, 4, 4, 4, 4, 4, 4, 4, 2], pad+exempt [2, 4, 4, 4, 4, 4, 4, 4, 2]
```

```
n=12, K=maxm=8
sweeps=0: |pad - unpadded| = 4.929e-01   |pad+exempt - unpadded| = 1.299e-14   |unpadded+exempt - unpadded| = 4.663e-15
          errors vs ED: unpadded 4.929e-01, pad 4.910e-08, pad+exempt 4.929e-01; final links pad [2, 4, 8, 8, 8, 8, 8, 8, 8, 4, 2], pad+exempt [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
sweeps=3: |pad - unpadded| = 1.911e-07   |pad+exempt - unpadded| = 7.242e-10   |unpadded+exempt - unpadded| = 0.000e+00
          errors vs ED: unpadded 1.420e-07, pad 4.910e-08, pad+exempt 1.427e-07; final links pad [2, 4, 8, 8, 8, 8, 8, 8, 8, 4, 2], pad+exempt [2, 4, 8, 8, 8, 8, 8, 8, 8, 4, 2]
```

**Suggested fix**: exempt the one-site route from padding the way `2d86cff`
exempted the MPO: `backend.pad_bonds_suspended()` inside
`Chain.global_subspace_expand` and the `num_center==1` branch of
`Chain.tdvp_step` (what `tdz.py` calls), plus one lossless SVD strip of the
padded zeros at trajectory entry in `quench_tdvp_gse`/`evolve_and_measure_tdvp_gse`
when padding is on, never per call inside `tdvp_step` (TDZ advances one step per
call and carries the wavefunction between calls, so a per-call strip would delete
the previous step's GSE directions). Measured by monkeypatch above: padded back
onto unpadded to 1.1e-15 and 1.3e-14 at `sweeps=0`, to 8.3e-8 and 7.2e-10 at
`sweeps=3`, and a no-op with padding off. The hunter's first fix (m as the true
rank of S1 in `_gse_bond_step`) is correct and necessary but does not cover
`sweeps=0`; its second fix (keep padded directions dead in the one-site
`qr_split`) is wrong, since alone it freezes padded `TDVP_GSE` at bond dimension
1 (error 0.493) and it also changes unpadded numbers (2.5e-6 at n=10, K=4),
because unpadded one-site TDVP's rank-deficient QR completion is what keeps a GSE
zero-weight direction alive and `qr_split` cannot tell the two apart. At minimum,
`svd.py`'s comment and `set_pad_bonds`' docstring should say that exactness holds
for the state and not for one-site TDVP. NUMBERS CHANGE: every padded `TDVP_GSE`
run, back onto the unpadded one.

## Ruled out

Executed and found clean, so the next pass need not redo them.

- **kpm-calibration**: the moment count n is identical on v2, v3, `"python"`, ED
  and `julia_live` for the same chain, delta and `kpm_scale` (74 at 0.7, 159 at
  1.5 on a 6-site staggered chain); every live reader of `kpm_n_scale` goes
  through `polynomials_for_broadening`; every public ED KPM route receives the
  `kpm_scale`/`kpm_n_scale` push, including `mode.py`'s fallbacks,
  `dynamicstk/spincorrelators`, `atomtk/iets` and the Kondo route; ED reproduces
  an exact numpy Jackson curve on complex Lehmann weights, imaginary part
  included, to 1e-14, so KPM has no conjugation or ordering problem.
- **td-convention**: the adjoint run inherits every knob on seven
  backend/integrator combinations (v2 MPO, v3 MPO, v3 TEBD, v3 `TDVP_GSE`,
  `"python"` TDVP, MPO and TEBD), all within 5.6e-05 to 9.0e-05 of exact on a
  purely imaginary-weight pair where `Re F` alone is 100 per cent off; a
  self-adjoint pair costs one evolution and any other two, on TD and TDZ
  (parafermion pairs cost two but are right); the Kondo code's and the examples'
  `.real`/`np.abs` all sit on dagger pairs, where the density is real.
- **canonical-form**: 0 false Hermiticity proofs and 0 false zero proofs in
  16000 ED-anchored checks over spin (S=1/2, S=1), spinless, native spinful,
  boson (including projectors) and parafermion operators; `simplify()` never
  replaces a stored operator; every "not proven" consumer falls back correctly;
  nothing imports sympy any more and `pyproject.toml` no longer lists it.
- **infinite-chain**: the v3 per-momentum resolvent cache is cleared when the
  Hamiltonian changes (a re-solve on the same session equals a fresh chain to
  every digit, dense and Lanczos); v3's deflated Lanczos agrees with its dense
  path to 1e-14 to 1e-15 away from k=0; `backend.set_backend("jax")` forces
  `jax_enable_x64` (`backend.py:94`), and a seeded JAX iDMRG cell agrees with
  NumPy to 7.3e-15; `_TransferChain` agrees with the dense route to 2.7e-15 on
  `vev`, correlators and `imps_overlap`, and nothing indexes it as a list; the
  bond-candidate-first fixed point moves a converged energy by 2.6e-14; an XXZ
  chain with a Dzyaloshinskii-Moriya term reproduces the exact k-asymmetric
  one-magnon dispersion to 1e-14 on both backends, which pins the sign of k that
  every existing test (E(k) = E(-k)) could not see.
- **recent-misc**: the Kondo DMRG route was not moved by the O1/O2 convention
  change (KPM sum rule exactly 0.250000 at 0, 2 and 10 T; second order against
  ED at 0.01 to 0.39 per cent at delta=2e-5); `T>0` with `mode="DMRG"` raises as
  documented; every static and dynamical route agrees padded against unpadded to
  1e-14 to 1e-15 and the MPO exemption covers every MPO route; the pre-`765b537`
  `gse.py` against the current one on NumPy gives exactly 0.0 over every site
  tensor and five seeds, the `tebd.py` diff is docstring only, and the JAX
  `TDVP_GSE` and TEBD tests pass against NumPy.

## New leads, not reviewed

Observed along the way and not handed to a reviewer, so they are leads rather
than findings.

- The KPM band-edge estimate `emax` is unconverged and varies from run to run:
  on the 4-site staggered Heisenberg chain it comes out 1.1e-4 to 1.9e-2 below
  the exact 0.8602399 (the reviewer of finding 4), moving KPM curves at the 1e-3
  level. After finding 4's fix it sets the ED-versus-DMRG residual.
- The iDMRG 2-site cell on the critical XX chain converges (`converged=True`) to a
  dimerized state, bonds -0.3244/-0.3121 at `maxm=12`, with an energy of
  -0.3182336 above the uniform VUMPS state's -0.3182820 at the same D; v3
  reproduces it to 1e-10 and the dimerization shrinks with `maxm`. No anchor tells
  an algorithmic fixed point from a defect here.
- `second_order_dIdV_dc(mode="ED", submode="KPM")` gives 2.0 where
  `submode="ED"` gives 1.50 at a degenerate ground state, the two ED correlator
  conventions (cached state against `dex` average) disagreeing; within `dex`'s
  documented scope, recorded for completeness.
- After finding 3's fix, ED and DMRG `get_distribution` still disagree pointwise
  (2.34 on a DMRG peak of 7.58 at `scale=2`), from their separate, documented
  off-calibration moment counts.
- `sxt_to_skomega` at `td_dynamical_correlator`'s defaults has delta*T = 1, so its
  damped series is cut at e^-1 (finding 7's reviewer, by reading); no test pins
  S(k,omega) values.
- `pychain/build.py:66` ends an unknown-name lookup in `print(name); raise`, so
  `mode="ED"` dies with "No active exception to reraise" on `ISy` or `Adag`; not in
  the bare-raise lists of either earlier record.
- `SpinTwoSite`'s `S2` in the vendored headers carries the documented upstream
  (5,4) typo; dmrgpy never requests it.
- `secondorder_dc.py`'s "0.2% at every bias point" measures 0.46 per cent at
  most on its own documented setup (median 4.8e-06).
- Stale text: `dynamics.py:69` to `72` still says TD/TDZ "are the routes still off
  this convention"; `_fourier_transform_correlator` keeps the superseded
  `S_AB = -(1/pi) Im G_AB` comment; `infinitechain.py:1220` to `1221` still says
  the Hermiticity check false-rejects a Heisenberg-style term.

## Statements this hunt overturns

Each of these was written as established and is contradicted by an entry above;
when the corresponding fix lands, the older text should gain a pointer here.

- The O1 Status paragraph of `audit_2026_09_hole_hunt.md` (near line 2228), the
  `765b537` commit message, `dynamics.py:129` and the user guide: TDZ's 3.31e-02
  is not its contour, and the `alpha0`/`n_max` remedy does nothing (finding 7).
- The same O1 paragraph and `lehmann_density_from_one_sided`'s docstring: refusing
  an unknown adjoint is not "never a wrong number" (finding 2).
- The O2 Status paragraph's "4 to 7 per cent" ED-versus-DMRG agreement, and the
  docstrings of `polynomials_for_broadening` and
  `kpmdmrg.dynamical_correlator_moments` (finding 4); O2's "about 1.6x higher at
  equal FWHM" (finding 6).
- The 2026-08 record's finding #14 reviewer note that `get_distribution(mode="ED")`
  "returns a correct spectrum" (finding 3).
- `ROADMAP.md:386` to `391` and the five documentation passages listed under
  finding 15, that the `"python"` excitation solver is not exposed to the
  degenerate-copy miss, and the `vx_deterministic_start` comment's diagnosis.
- `canonical.py`'s comment that the post-transform names are "hence even"
  (finding 1), and `svd.py`'s that padded columns "contribute nothing to any
  contraction" (finding 16).

## Shared helpers

Modules imported by the scripts above, verbatim.

`td-convention/anchor.py` (findings 2 and 8):

```python
"""Exact references built by hand with numpy kron, no dmrgpy code."""
import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
I2 = np.eye(2, dtype=complex)


def site_op(o, i, n):
    out = np.array([[1.0 + 0j]])
    for k in range(n):
        out = np.kron(out, o if k == i else I2)
    return out


def heis_field(n, hx=0.0, hz=0.0, hy=0.0, J=1.0):
    H = 0
    for i in range(n - 1):
        for o in (sx, sy, sz):
            H = H + J * site_op(o, i, n) @ site_op(o, i + 1, n)
    for i in range(n):
        H = H + hx * site_op(sx, i, n) + hz * site_op(sz, i, n) \
            + hy * site_op(sy, i, n)
    return H


def lehmann_dense(H, A, B):
    e, U = np.linalg.eigh(H)
    Uh = U.conj().T
    Ae = Uh @ A @ U
    Be = Uh @ B @ U
    return e - e[0], Ae[0, :] * Be[:, 0], e


def density(D, M, es, delta):
    return np.array([np.sum(M * (delta / np.pi) / ((w - D) ** 2 + delta ** 2))
                     for w in es])
```

`td-convention/leadE_fft_grid_helpers.py` (finding 7):

```python
import numpy as np
from dmrgpy import fermionchain
def complex_hopping_chain(n=4, itensor_version=3, seed=3):
    # verbatim from tests/test_audit_2026_09_correlator-conventions.py
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    rng = np.random.RandomState(seed)
    t = rng.random((n, n)) + 1j * rng.random((n, n))
    t = t + t.conj().T
    h = 0
    for i in range(n):
        for j in range(n): h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1): h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 30, 12
    return fc

def lehmann(chain, A, B):
    ed = chain.get_ED_obj()
    h = np.array(ed.get_hamiltonian().todense())
    emu, vs = np.linalg.eigh(h)
    U = np.array(vs); Uh = np.conjugate(U.T)
    Ae = Uh @ np.array(ed.MO2matrix(A).todense()) @ U
    Be = Uh @ np.array(ed.MO2matrix(B).todense()) @ U
    return emu - emu[0], Ae[0, :] * Be[:, 0]

def density(D, M, es, delta):
    return np.array([np.sum(M * (delta / np.pi) / ((w - D) ** 2 + delta ** 2))
                     for w in es])
```

`review-td-1/revhelpers.py` and `review-td-1/probe1_complex_chain_build.py`
(finding 7):

```python
"""Reviewer's helpers: a scratchpad copy of timedependent._fourier_transform_
correlator with two variants of its last stage, and exact references built
by dense diagonalization with numpy only (no dmrgpy ED code).

variant "orig"   -> calls the repository function unchanged
variant "pad"    -> identical body, FFT zero-padded to K*len(cs)
variant "direct" -> identical body, the trapezoid sum evaluated at each es
                    directly (no FFT grid, no interpolation)
Every call logs the raw (ts, cs) it was given and the damped tail ratio."""
import numpy as np
from scipy.interpolate import interp1d
import dmrgpy.timedependent as tdm
import dmrgpy.tdz as tdz

ORIG = tdm._fourier_transform_correlator
STATE = {"variant": "orig", "K": 16, "log": [], "tail": []}


def patched(ts, cs, dt, es=None, window=[-1, 10], delta=5e-2, factor=1,
            damping="exp", predict=False, lp_order=20, lp_extend_factor=10,
            lp_fit_start_fraction=0.5, lp_max_pole_radius=1.0):
    STATE["log"].append((np.array(ts), np.array(cs)))
    kw = dict(es=es, window=window, delta=delta, factor=factor,
              damping=damping, predict=predict, lp_order=lp_order,
              lp_extend_factor=lp_extend_factor,
              lp_fit_start_fraction=lp_fit_start_fraction,
              lp_max_pole_radius=lp_max_pole_radius)
    if STATE["variant"] == "orig":
        return ORIG(ts, cs, dt, **kw)
    # ---- verbatim body of the repository function up to the FFT ----
    if predict:
        from dmrgpy.dynamicstk.linearprediction import linear_predict_extend
        ts, cs = linear_predict_extend(ts, cs, order=lp_order,
                                       extend_factor=lp_extend_factor,
                                       fit_start_fraction=lp_fit_start_fraction,
                                       max_pole_radius=lp_max_pole_radius)
    cs = cs * tdm._damping_window(ts, delta, damping=damping)
    STATE["tail"].append((np.abs(cs[-1]) / np.max(np.abs(cs)),
                          np.exp(-delta * np.max(ts))))
    ftr = interp1d(ts, cs.real, fill_value=0.0, bounds_error=False)
    fti = interp1d(ts, cs.imag, fill_value=0.0, bounds_error=False)
    tnew = np.linspace(np.min(ts), np.max(ts), len(ts) * factor)
    cs = ftr(tnew) + 1j * fti(tnew)
    ts = tnew.copy()
    dtnew = dt / factor
    cs = cs.copy()
    cs[0] = cs[0] * 0.5
    cs[-1] = cs[-1] * 0.5
    if es is None:
        es = np.linspace(window[0], window[1], 800)
    es = np.asarray(es, dtype=float)
    # ---- the only changed stage ----
    if STATE["variant"] == "pad":
        n = STATE["K"] * len(cs)
        ss = np.fft.fft(cs, n=n) * dtnew / np.pi
        ws = np.fft.fftfreq(n, d=dtnew) * 2. * np.pi
        fr = interp1d(ws, ss.real, fill_value=0.0, bounds_error=False)
        fi = interp1d(ws, ss.imag, fill_value=0.0, bounds_error=False)
        gr = fr(es) + 1j * fi(es)
    elif STATE["variant"] == "direct":
        tt = ts - ts[0]            # the FFT's own time origin
        gr = (np.exp(-1j * np.outer(es, tt)) @ cs) * dtnew / np.pi
    else:
        raise ValueError(STATE["variant"])
    return (es, gr)


def install():
    tdm._fourier_transform_correlator = patched
    tdz._fourier_transform_correlator = patched   # tdz imported the name


def use(variant, K=16):
    STATE["variant"] = variant
    STATE["K"] = K
    STATE["log"] = []
    STATE["tail"] = []


def lehmann_dense(H, A, B):
    e, U = np.linalg.eigh(H)
    Uh = U.conj().T
    Ae = Uh @ A @ U
    Be = Uh @ B @ U
    return e - e[0], Ae[0, :] * Be[:, 0]


def density(D, M, es, delta):
    return np.array([np.sum(M * (delta / np.pi) / ((w - D) ** 2 + delta ** 2))
                     for w in es])


def ctime(D, M, ts, sign):
    """sum_n M_n exp(sign*i*D_n*t)"""
    return np.exp(sign * 1j * np.outer(ts, D)) @ M
```

```python
import numpy as np
import revhelpers as rh
from dmrgpy import fermionchain
N, SEED = 4, 3

def build(n=N, seed=SEED):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=3)
    rng = np.random.RandomState(seed)
    t = rng.random((n, n)) + 1j * rng.random((n, n))
    t = t + t.conj().T
    h = 0
    for i in range(n):
        for j in range(n): h = h + t[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1): h = h + 0.8 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 30, 12
    # hand-built JW reference, no dmrgpy code
    a = np.array([[0, 1], [0, 0]], dtype=complex)
    Z = np.diag([1.0, -1.0]).astype(complex)
    I2 = np.eye(2, dtype=complex)
    def c(i):
        out = np.array([[1.0 + 0j]])
        for k in range(n):
            out = np.kron(out, Z if k < i else (a if k == i else I2))
        return out
    C = [c(i) for i in range(n)]
    Cd = [x.conj().T for x in C]
    Nn = [Cd[i] @ C[i] for i in range(n)]
    H = sum(t[i, j] * Cd[i] @ C[j] for i in range(n) for j in range(n))
    H = H + sum(0.8 * Nn[i] @ Nn[i + 1] for i in range(n - 1))
    e = np.linalg.eigvalsh(H)
    D, M = rh.lehmann_dense(H, Cd[0], C[2])
    return fc, D, M, e[1] - e[0]
```

`review-kpm-H4/exactlehmann.py` (finding 6):

```python
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla

def sector(L):
    st = np.array([s for s in range(1 << L) if bin(s).count("1") == L//2])
    idx = {int(s): i for i, s in enumerate(st)}
    r, c, v = [], [], []
    diag = np.zeros(len(st))
    for a, s in enumerate(st):
        s = int(s)
        for i in range(L-1):
            bi, bj = (s >> i) & 1, (s >> (i+1)) & 1
            if bi == bj: diag[a] += 0.25
            else:
                diag[a] -= 0.25
                r.append(idx[s ^ ((1 << i) | (1 << (i+1)))]); c.append(a); v.append(0.5)
    H = sp.csr_matrix((v, (r, c)), shape=(len(st), len(st))) + sp.diags(diag)
    return st, H

def szop(st, i): return np.where((st >> i) & 1, 0.5, -0.5)

def lanczos_poles(H, v0, m):
    nrm = np.linalg.norm(v0); V = np.zeros((m, len(v0))); a = np.zeros(m); b = np.zeros(m)
    V[0] = v0/nrm
    for k in range(m):
        w = H @ V[k]; a[k] = V[k] @ w
        w -= V[:k+1].T @ (V[:k+1] @ w); w -= V[:k+1].T @ (V[:k+1] @ w)
        if k+1 < m:
            b[k] = np.linalg.norm(w); V[k+1] = w/b[k]
    T = np.diag(a) + np.diag(b[:m-1], 1) + np.diag(b[:m-1], -1)
    e, U = np.linalg.eigh(T)
    return e, nrm**2*U[0]**2
```

`review-canonical-1/fix_none.py`, `fix_narrow.py` and
`run_leadF_with_narrow_fix.py` (finding 1):

```python
"""pytest plugin: the hunter's fix option 1, parity None for the six A-type names"""
from dmrgpy.multioperatortk import canonical
for n in ("A","Adag","Aup","Adagup","Adn","Adagdn"): canonical._PARITY[n] = None
```

```python
"""pytest plugin: the hunter's fix option 2, refuse (leave unreordered) any
term containing both a C-type (_ODD) name and an A-type name"""
from dmrgpy.multioperatortk import canonical
_A_TYPE = frozenset(("A","Adag","Aup","Adagup","Adn","Adagdn"))
_orig = canonical._canonical_signature
def _narrow(term):
    names = [o[0] for o in term[1:]]
    if any(n in canonical._ODD for n in names) and any(n in _A_TYPE for n in names):
        return tuple((o[0],o[1]) for o in term[1:]), term[0]
    return _orig(term)
canonical._canonical_signature = _narrow
```

```python
"""Run the td-convention lens's leadF script with the narrow canonical fix
(fix_narrow.py: refuse any term mixing an _ODD name with an A-type name)
applied first, to see whether that fix alone closes the TD consequence."""
import sys, runpy
sys.path.insert(0, "<scratch>/review-canonical-1")
import fix_narrow  # patches canonical._canonical_signature
runpy.run_path("<scratch>/td-convention/leadF_parity_false_positive.py", run_name="__main__")
```
