# Audit, 2026-08: cross-backend hole hunt

Findings from a five-lens automated audit of the Python layer, run 2026-08-30.
Each lens hunted one failure class -- silent backend/method dispatch fallbacks,
dropped keyword arguments, feature x feature combinations, garbage-in-garbage-
out cases that should raise, and cross-backend numerical disagreement on chains
small enough for ED to be exact. Every finding below carries a repro that was
actually executed, and every one marked CONFIRMED was then re-run by an
independent reviewer whose brief was to refute it.

This file is the evidence, not a task list: it records what was observed and
how to reproduce it, so a fix (or a decision that the behaviour is intended
after all) does not have to re-derive any of it. Fixed entries are marked as
such rather than deleted, since the repro doubles as the regression check.

Scope notes: findings against vendored ITensor (`mpscpp2/ITensor/`,
`mpscpp3/ITensor/`) were out of scope, as were the deliberately reproduced
legacy bugs CLAUDE.md lists (`evoloperator`'s z^3/6 term on `H2`, the `"moise"`
key, the unreachable `"tevol_fit_td"` branch) and gaps ROADMAP.md already marks
absent. `itensor_version="julia_live"` was not executed (juliacall JIT cost),
so no finding below is about the Julia backend.

## Reproducing anything here

Every repro in this file was run as:

```bash
cd <scratch>
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=<checkout>/src python3 repro.py
```

The explicit `PYTHONPATH` is load-bearing, not decoration: `site-
packages/dmrgpy` is typically a symlink into one specific checkout, so a bare
`python3` can silently test a different working tree and manufacture (or hide)
a finding. The BLAS pinning keeps timings and, for the segfault cases, thread
counts comparable.

## Confirmed findings

Reproduced by the finder, then reproduced again and argued against by an independent reviewer, and still standing.

### 1. Non-Hermitian Hamiltonian silently discards submode=: every submode returns the CVM_explicit resolvent instead

- **Severity**: medium  **Class**: silent-wrong
- **Code**: `src/dmrgpy/dynamics.py:16 (and src/dmrgpy/edtk/dynamics.py:50)`
- **Status**: FIXED -- dynamics.py and edtk/dynamics.py now branch per submode: KPM and CVM/CVM_explicit have real non-Hermitian implementations, EX/maxent pass through, everything else raises NotImplementedError.

**Expected**

dynamics.py checks `if not self.is_hermitian(...)` BEFORE the submode dispatch
and, for anything other than submode="KPM", returns
dynamical_correlator_non_hermitian (= nonhermitian/dynamics.py's
dynamical_correlator_cvm_explicit). So on a non-Hermitian Hamiltonian the
user's submode= is thrown away entirely: "EX", "maxent", "ROOTN", "TD", "CVM",
"CVM_explicit" all return the same explicit-resolvent number (they agree to
2.8e-8, the CG solver's own noise, on the DMRG path and are bitwise identical
on the ED path). ROADMAP.md's Sec. 4 marks "EX" and "maxent" as backend-
agnostic and explicitly "Non-Hermitian-capable", so these are documented as
supported, not absent. The same file's julia_live branch contains a comment
stating this exact pre-dispatch hijack "was blocking a working code path" and
was therefore replaced by an explicit NotImplementedError there -- the
(2,3,"python") branch and the ED branch (edtk/dynamics.py:46-50, identical
structure) still silently substitute a different algorithm instead of either
dispatching or raising.

**Observed**

```
is_hermitian: False
Non Hermitian mode in dynamical correlator
EX -> y[:2] = [0.0641056  0.14047283]
Non Hermitian mode in dynamical correlator
CVM_explicit -> y[:2] = [0.0641056  0.14047283]
Non Hermitian mode in dynamical correlator
ROOTN -> y[:2] = [0.0641056  0.14047283]
Non Hermitian mode in dynamical correlator
maxent -> y[:2] = [0.0641056  0.14047283]
Non Hermitian mode in dynamical correlator
TD -> y[:2] = [0.0641056  0.14047283]
Non Hermitian mode in dynamical correlator
CVM -> y[:2] = [0.0641056  0.14047283]
dcex ever reached for submode='EX': False
...
direct dcex.dynamical_correlator -> TypeError: mpsiram() got an unexpected keyword argument 'scale'

(r2c.py)  max|EX-CVM_explicit| = 2.78623523246313e-08
(r9_extra.py)  ED path: max|EX-CVM| = 0.0
```

**Reviewer** (independent, brief was to refute)

Confirmed, and not explainable as intended. dynamics.py:11-16 tests `not
self.is_hermitian(...)` BEFORE the submode dispatch and, for everything except
submode="KPM", returns nonhermitian/dynamics.py's
dynamical_correlator_cvm_explicit; edtk/dynamics.py:46-50 has the identical
pre-dispatch structure. So on a non-Hermitian H the user's submode= is a no-op.
I reproduced this myself on three backends: itensor_version=3
(EX/CVM_explicit/ROOTN/maxent/TD/CVM agree to ~2-8e-8, the CG solver's own
noise), itensor_version="python" (~2e-7), and mode="ED" where
EX/CVM/ED/INV/ROOTN/TD are bitwise identical (max|diff| = 0.0). The agreement
is not coincidence between converged algorithms: it is literally the same
function object each time (dcex.dynamical_correlator is never reached, verified
by monkeypatch), so DMRG convergence is irrelevant and cannot be the false-
positive cause; the ED path's exact 0.0 proves it directly. Not intended: the
2024 comment "# submode = CVM_explicit # only mode that works with non-
Hermitian" records the original intent, but the SAME function's julia_live
branch carries a later comment from the same author repudiating exactly this
structure ("This check used to run before submode dispatch entirely, which also
rejected EX/maxent even though they work fine -- confirmed directly, this was
blocking a working code path") and replaces it with a submode-aware
NotImplementedError; ROADMAP.md:99 marks EX and maxent "Non-Hermitian-capable"
for v3/pyitensor/Julia; docs/user_guide.md:1490 documents ONLY the KPM re-
routing, never a general one. The (2,3,"python") and ED branches were simply
not updated alongside the Julia one. Blast radius: all of itensor_version
2/3/"python" plus the ED backend, via the public
get_dynamical_correlator(submode=...) for EX, maxent, ROOTN, TD, CVM,
CVM_explicit, CVMimag, TDZ, and (ED only) INV and "ED" itself -- the last being
the worst case, since the exact Lehmann reference a user would cross-validate
against silently becomes the approximate resolvent and agrees bit-for-bit while
validating nothing. Two findings cap the severity below "high", both checked:
(1) the substitution is silent only for the shared kwargs -- any submode-
specific argument fails loudly (submode="EX", nex=6 -> TypeError:
dynamical_correlator_cvm_explicit() got an unexpected keyword argument 'nex';
likewise nt, n); (2) no existing test or example is currently rendered vacuous,
since every examples/non_hermitian/* correlator script uses the default KPM
route, which is correctly special-cased. The value returned is also a valid
computation of what it actually runs, not garbage. Additional latent defect the
hijack is masking, which I confirmed: dispatching EX properly would not work
today either -- dcex.py:54 forwards its own scale=10.0 through
get_excited_states -> excited.py:85 excited_states_non_hermitian -> **kwargs ->
arpacktk.excited_states -> mpsiram, which has no catch-all, giving "TypeError:
mpsiram() got an unexpected keyword argument 'scale'". So ROADMAP's "EX is Non-
Hermitian-capable" is itself inaccurate on these backends.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/r2b_nh.py:

"""Non-Hermitian H: every submode silently returns the CVM_explicit result."""
import numpy as np
from dmrgpy import spinchain, dcex, dynamics

n = 4
sc = spinchain.Spin_Chain(["S=1/2"]*n)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
h = h + 0.3j*sc.Sz[0]
sc.set_hamiltonian(h)
sc.maxm = 20; sc.nsweeps = 6
print("is_hermitian:", sc.is_hermitian(sc.hamiltonian))
es = np.linspace(0.2, 3., 20)
name = [sc.Sz[0], sc.Sz[0]]

hit = []
orig = dcex.dynamical_correlator
dynamics.dcex.dynamical_correlator = lambda *a, **k: (hit.append(1), orig(*a, **k))[1]

out = {}
for sub in ["EX", "CVM_explicit", "ROOTN", "maxent", "TD", "CVM"]:
    _, y = sc.get_dynamical_correlator(name=name, submode=sub, es=es, delta=0.1)
    out[sub] = np.asarray(y)
    print(sub, "-> y[:2] =", out[sub][:2])
print("dcex ever reached for submode='EX':", bool(hit))
for sub in out:
    print(sub, "identical to CVM_explicit:",
          np.array_equal(out[sub], out["CVM_explicit"]))
try:
    dcex.dynamical_correlator(sc, name=name, nex=4, es=es, delta=0.1)
except Exception as e:
    print("direct dcex.dynamical_correlator ->", type(e).__name__+":", e)

Command:
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r2b_nh.py

Quantifying the agreement (r2c.py):

import numpy as np
from dmrgpy import spinchain
n = 4
sc = spinchain.Spin_Chain(["S=1/2"]*n)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
h = h + 0.3j*sc.Sz[0]
sc.set_hamiltonian(h); sc.maxm=20; sc.nsweeps=6
es = np.linspace(0.2,3.,20); name=[sc.Sz[0],sc.Sz[0]]
_,a = sc.get_dynamical_correlator(name=name,submode="EX",es=es,delta=0.1)
_,b = sc.get_dynamical_correlator(name=name,submode="CVM_explicit",es=es,delta=0.1)
print("max|EX-CVM_explicit| =", np.max(np.abs(np.asarray(a)-np.asarray(b))))

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r2c.py

Same check on the ED branch (r9_extra.py, first half):

import numpy as np
from dmrgpy import spinchain
n = 4
sc = spinchain.Spin_Chain([2]*n)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h + 0.3j*sc.Sz[0]); sc.maxm=20; sc.nsweeps=6
es = np.linspace(0.2,3.,20); name=[sc.Sz[0],sc.Sz[0]]
_,a = sc.get_dynamical_correlator(mode="ED",name=name,submode="EX",es=es,delta=0.1)
_,b = sc.get_dynamical_correlator(mode="ED",name=name,submode="CVM",es=es,delta=0.1)
print("ED path: max|EX-CVM| =", np.max(np.abs(np.asarray(a)-np.asarray(b))))

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r9_extra.py
```

### 2. Any unrecognized tevol_method string silently downgrades real-time evolution to the legacy MPO-Taylor integrator (600x larger error)

- **Severity**: medium  **Class**: silent-wrong
- **Code**: `src/dmrgpy/timedependent.py:159 (and :237; same else-branch in src/dmrgpy/tdz.py:174)`
- **Status**: FIXED -- timedependent.check_tevol_method() validates the name at all three dispatch sites; backend capability still falls back as documented.

**Expected**

evolve_and_measure_dmrg/evolution_dmrg_DC dispatch with a chain of `elif
self.itensor_version in (3,"python") and self.tevol_method=="..."` tests and an
unconditional `else: self._session.evolve_and_measure(...)` (2nd-order MPO-
Taylor). tevol_method is never validated, so ANY string that is not exactly one
of "TDVP"/"TEBD"/"AUTO"/"TDVP_GSE" -- a case typo, a trailing space, a
misspelling -- silently selects a completely different integrator whose error
against ED is 600x larger here (3.5e-2 vs 5.9e-5), with no message. The
julia_live backend validates the same attribute and raises for exactly this
input (mpsjulialive/timedependent.py's _check_tevol_method: "anything outside
(TDVP,TDVP_GSE) raises there instead of silently running plain TDVP", per
ROADMAP.md Sec. 3), so the fail-loud behavior already exists on one backend and
is missing on the default one. tdz.py:174 has the same untested else-branch
(its docstring even calls it out) so submode="TDZ" is affected identically.
This is distinct from the documented itensor_version=2 -> MPO fallback: a value
the user typed and that no backend implements should raise, not be silently
reinterpreted.

**Observed**

```
tevol_method='TDVP' max|DMRG-ED| = 5.853e-05
tevol_method='MPO'  max|DMRG-ED| = 3.544e-02
tevol_method='tdvp' max|DMRG-ED| = 3.544e-02
tevol_method='TEBD ' max|DMRG-ED| = 3.544e-02
tevol_method='TVDP' max|DMRG-ED| = 3.544e-02
tevol_method='typo' max|DMRG-ED| = 3.544e-02
'TDVP' bit-identical to 'MPO': False
'MPO'  bit-identical to 'MPO': True
'tdvp' bit-identical to 'MPO': True
'TEBD ' bit-identical to 'MPO': True
'TVDP' bit-identical to 'MPO': True
'typo' bit-identical to 'MPO': True
```

**Reviewer** (independent, brief was to refute)

REPRODUCED exactly with my own script (correct PYTHONPATH, dmrgpy loaded from
/u/40/ladovj1/data/Documents/programs/dmrgpy/src): TDVP 5.853e-05 vs 3.544e-02
for each of 'MPO','tdvp','TEBD ','TVDP','typo', with all five bit-identical to
'MPO'. Repro flaws ruled out: I tightened maxm=40/nsweeps=12 and the numbers
are unchanged, and convergence is irrelevant anyway since every unrecognized
string produces output bit-identical to the explicit 'MPO' run from the *same*
initial wavefunction -- this is a dispatch difference, not a numerical
disagreement. Also confirmed identical on itensor_version="python" (5.853e-05 /
3.544e-02, 'tdvp' bit-identical to 'MPO'). julia_live was not run.  DEFECT,
tightest form: self.tevol_method is never validated -- no property, no
__setattr__, no check at any of the three dispatch sites
(timedependent.py:138/216 chains, tdz.py:166). Each dispatch tests
`itensor_version in (3,"python") and tevol_method=="<literal>"` and falls to an
unconditional else. That conflates two distinct conditions: "this method is not
applicable on this backend" (documented and intended -- itensor_version=2 has
no TDVP/TEBD, and user_guide.md says v2 "falls back to 'MPO' even with this
default"; explicit tevol_method="MPO" is also documented) with "no such method
exists anywhere" (intended nowhere). The else branch itself is correct; the
missing input validation is the bug.  INTENDED-BEHAVIOUR DEFENSE, considered
and rejected. The strongest candidate is tdz.py:153's parenthetical "the way an
unrecognized/inapplicable tevol_method otherwise would" -- it proves the author
knew the mechanism exists, but "inapplicable" is the documented v2/MPO case,
and nothing in CLAUDE.md, ROADMAP.md, docs/user_guide.md, docs/documentation.md
or the source documents typo -> MPO as a designed alias. The project states the
opposite principle four times, one of them locked in by a test:
mpsjulialive/timedependent.py::_check_tevol_method's docstring ("Raised rather
than silently falling back to plain TDVP: a caller who set tevol_method='TEBD'
asked for a specific integrator, and quietly running a different one would make
a backend-comparison script silently compare the wrong things");
user_guide.md's "Five options" list plus "'TEBD' itself deliberately stays a
hard opt-in that raises rather than silently falling back" and "the legacy
'MPO' path raises NotImplementedError there rather than silently running plain
TDVP instead, so a backend-comparison script can't quietly end up comparing
different integrators"; ROADMAP.md row 81; and
tests/test_julia_live.py::test_unsupported_tevol_method_raises. Author
awareness of a fall-through is not documentation that it is intended. This is
also distinct from the deliberately-reproduced legacy bugs (evoloperator z^3/6,
"moise"/"tevol_fit_td"), which carry explicit comments saying so.  BLAST
RADIUS: itensor_version=3 (the default) and "python", bit-identical on both.
Public API: evolve_and_measure()/evolution_ABA() via evolve_and_measure_dmrg
(timedependent.py:237), dynamical_correlator(submode="TD") via
evolution_dmrg_DC (:159), and submode="TDZ" via tdz.py:174 (same else-pattern,
verified by code inspection, not re-run). julia_live raises for exactly this
input; itensor_version=2 is documented all-MPO, so neither is affected.
SEVERITY BOUNDING: not high -- it takes a user typo to trigger (no correct
usage hits it) and the substituted integrator is a documented, dt-convergent
method rather than corrupt data (MPO error falls 3.5e-2 -> 7.0e-4 as dt goes
0.2 -> 0.025). Not low -- fully silent, on the default backend, and directly
contradicts a design principle the project states in the user guide and
enforces with a test on a sibling backend; and the TDVP:MPO error ratio
*worsens* as dt shrinks (606x at dt=0.2, 2840x at dt=0.05), so the silent
substitution is relatively worse at production step sizes.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/r3_tevol.py:

"""An unrecognized tevol_method silently downgrades to MPO-Taylor."""
import numpy as np
from dmrgpy import spinchain, timedependent

n = 4
sc = spinchain.Spin_Chain([2]*n)
h0 = 0
for i in range(n): h0 = h0 + (-1)**i*sc.Sz[i]
h1 = 0
for i in range(n-1):
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
wf   = sc.get_gs()
wfED = sc.get_gs(mode="ED")
sc.set_hamiltonian(h1)
op = sc.Sz[0]
nt, dt = 30, 0.2

_, ed = timedependent.evolve_and_measure(sc, operator=op, nt=nt, dt=dt,
                                         wf=wfED, mode="ED")
res = {}
for m in ["TDVP", "MPO", "tdvp", "TEBD ", "TVDP", "typo"]:
    sc.tevol_method = m
    _, y = timedependent.evolve_and_measure(sc, operator=op, nt=nt, dt=dt, wf=wf)
    res[m] = np.asarray(y)
    print("tevol_method=%-6r max|DMRG-ED| = %.3e" % (m, np.max(np.abs(y-ed))))
for m in res:
    print("%-6r bit-identical to 'MPO': %s" % (m, np.array_equal(res[m], res["MPO"])))

Command:
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r3_tevol.py
```

### 3. kpm_energy_truncate=True with kpm_scale below ~0.5 (its documented purpose) silently returns a spectrum whose peak is at the wrong energy, on BOTH itensor_version=3 and "python"

- **Severity**: medium  **Class**: silent-wrong
- **Code**: `src/dmrgpy/pyitensor/chain.py:1688 (_scaled_hamiltonian_gs_anchored) and the twin native method src/dmrgpy/mpscpp3/chain_session.h::kpm_dynamical_correlator_truncated, dispatched from src/dmrgpy/kpmdmrg.py:106-121`
- **Status**: PARTIALLY GUARDED -- Holzner's Eq. (41) diagnostic is now checked and the worst regime raises, on both backends. The band between that and a healthy run is still silent: see docs/known_issue_kpm_energy_truncation_window.md.

**Expected**

The exact spectral peak of <Sz_2 Sz_2>(w) on this 6-site Heisenberg chain is at
0.525 (mode="ED", submode="ED" Lehmann sum), and plain KPM at the default
kpm_scale=0.7 reproduces it (0.441, one grid step off at the 0.085-wide grid
with delta=0.3 broadening). manybodychain.py:216-221 documents
kpm_energy_truncate as existing precisely so that "kpm_scale [can] be pushed
below its current safe floor (~0.5 ...) for higher KPM spectral resolution",
and pyitensor/chain.py::_scaled_hamiltonian_gs_anchored's docstring says the
retained window runs from E0 to E0+Ws. I printed that window via
get_dynamical_correlator_moments: at kpm_scale=0.5 it is [0, 1.872] above E0
and at 0.4 it is [0, 1.497] -- the true peak at 0.525 sits comfortably INSIDE
both, so this is not "weight outside the window is expected garbage". Yet the
returned spectrum's peak jumps to 1.797 and 1.458 (essentially pinned to the
top edge of the retained window) with 3x-7x the amplitude, and no exception, no
warning, no NaN. Anchoring on peak position rather than amplitude, since plain-
KPM vs ED amplitudes already differ by kernel convention. The two independent
implementations (mpscpp3's native kpm_dynamical_correlator_truncated and
pyitensor's setter+branch) agree digit-for-digit, so this is the shared
algorithm being wrong in exactly the regime the feature was added for, not a
porting bug in one backend. Either the anchored rescaling/truncation should
produce the correct in-window spectrum, or the call should raise instead of
returning a plausible-looking but wrongly-peaked curve. Also reproduced at
kpm_truncate_dK=20, kpm_truncate_nsweeps=6, maxm=kpmmaxm=40, so it is not an
under-convergence artifact.

**Observed**

```
EXACT (Lehmann, mode=ED submode=ED): peak at e = 0.525
v=3       trunc=True kpm_scale=0.7 : retained window = [0, 2.621] above E0 ; peak at e = 0.441 ; max|y| = 0.3650
v=3       trunc=True kpm_scale=0.5 : retained window = [0, 1.872] above E0 ; peak at e = 1.797 ; max|y| = 0.6924
v=3       trunc=True kpm_scale=0.4 : retained window = [0, 1.497] above E0 ; peak at e = 1.458 ; max|y| = 1.6930
v=python  trunc=True kpm_scale=0.7 : retained window = [0, 2.621] above E0 ; peak at e = 0.441 ; max|y| = 0.3664
v=python  trunc=True kpm_scale=0.5 : retained window = [0, 1.872] above E0 ; peak at e = 1.797 ; max|y| = 0.6862
v=python  trunc=True kpm_scale=0.4 : retained window = [0, 1.497] above E0 ; peak at e = 1.458 ; max|y| = 1.6930
v=3       trunc=False kpm_scale=0.7 (default)                                   : peak at e = 0.441 ; max|y| = 0.2382
```

**Reviewer** (independent, brief was to refute)

Reproduced exactly on both backends. Tightest statement: when
kpm_energy_truncate=True and the ground-state-anchored window [E0, E0+Ws]
leaves even ~8-10% of the correlator's spectral weight above its top edge, the
per-site Krylov energy truncation does not remove that weight, it piles it at
the window edges; the Chebyshev moments then stay bounded by mu0 (measured: at
kpm_scale=0.4 they are exactly periodic 0.25,0.141,0.25,0.141,... i.e. all
weight at rescaled x=+-1), so neither _check_kpm_moment nor mpscpp3's
check_kpm_moment can fire, and get_dynamical_correlator returns a plausible-
looking curve peaked at ~0.96x the window top with 3x-23x the correct
amplitude, no warning. I ruled out the obvious refutations: (a) not under-
convergence -- identical at the repo's own recommended
kpm_truncate_dK=30/nsweeps=10, and the untruncated path at the same kpm_scale
gives the correct peak; (b) not too few moments -- the cliff between kpm_scale
0.60 (peak 0.525, correct) and 0.55 (peak 1.966 = window edge) happens at
constant n=21; (c) not the documented "the window cuts real physical weight"
limitation -- ED integration shows 92.4% of the weight inside the 0.60 window
vs 90.3% inside the 0.55 window, i.e. two percentage points of discarded weight
flip right to garbage, a hard cliff rather than gradual lineshape degradation;
(d) the mechanism itself is sound where the paper's premise holds -- on a
7-site chain with a decoupled 20*Sz[6] term (bandwidth inflated, correlator
width unchanged) truncation reproduces the exact peak 0.525 down to
kpm_scale=0.15. So the defect is the absence of any validity check on the
discarded weight, in a regime the documentation positively advertises:
docs/user_guide.md:1145 shows sc.kpm_scale = 0.3 as its worked snippet, which
on the repo's own 6-site Heisenberg example returns a peak at 0.017 with
max|y|=2.84 (23x the exact amplitude). Both implementations already compute the
paper's Eq. (40)/(41) diagnostics (avg_truncated_weight, state_change_norm) and
discard them at the call site (`.first` in mpscpp3/chain_session.h's
kpm_moments_truncated_full/_accelerated, ignored return in
pyitensor/chain.py::_maybe_energy_truncate) -- exactly the quantity that would
flag this. Reach: public get_dynamical_correlator(submode="KPM") and
get_dynamical_correlator_moments() on itensor_version=3 and "python" only, and
only with the opt-in, off-by-default kpm_energy_truncate=True;
itensor_version=2 (no-op) and ED are unaffected. Partial mitigation:
tests/test_kpm_energy_truncation_accuracy.py's module docstring already notes
"kpm_scale=0.5 already visibly distorts the lineshape here" and calls it an
expected small-system limitation -- so the degradation was known -- but it is
characterized there as lineshape distortion, not as a peak teleporting to the
window edge with several times the correct amplitude, and no user-facing
document states any floor at all.

**Repro**

```
File /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/f1.py:

import numpy as np, io, contextlib
from dmrgpy import spinchain
def build(v,**kw):
    sc = spinchain.Spin_Chain([2]*6,itensor_version=v)
    h=0
    for i in range(5):
        h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm=40; sc.nsweeps=12; sc.kpmmaxm=40
    for k,val in kw.items(): setattr(sc,k,val)
    return sc
es = np.linspace(-1,4,60)
name = lambda sc:(sc.Sz[2],sc.Sz[2])
sc = build(3)
with contextlib.redirect_stdout(io.StringIO()):
    x,y = sc.get_dynamical_correlator(mode="ED",submode="ED",name=name(sc),es=es,delta=0.3)
print("EXACT (Lehmann, mode=ED submode=ED): peak at e = %.3f"%x[np.argmax(np.abs(y))])
for v in [3,"python"]:
    for sca in [0.7,0.5,0.4]:
        sc = build(v,kpm_energy_truncate=True,kpm_scale=sca)
        with contextlib.redirect_stdout(io.StringIO()):
            mus,emin,emax,scale,n,delta = sc.get_dynamical_correlator_moments(
                    name=name(sc),delta=0.3)
        sc = build(v,kpm_energy_truncate=True,kpm_scale=sca)
        with contextlib.redirect_stdout(io.StringIO()):
            x,y = sc.get_dynamical_correlator(name=name(sc),es=es,delta=0.3)
        print("v=%-7s trunc=True kpm_scale=%.1f : retained window = [0, %.3f] above E0 ; peak at e = %.3f ; max|y| = %.4f"
              %(v,sca,emax-emin,x[np.argmax(np.abs(y))],np.max(np.abs(y))))
sc = build(3,kpm_energy_truncate=False,kpm_scale=0.7)
with contextlib.redirect_stdout(io.StringIO()):
    x,y = sc.get_dynamical_correlator(name=name(sc),es=es,delta=0.3)
print("v=3       trunc=False kpm_scale=0.7 (default)                                   : peak at e = %.3f ; max|y| = %.4f"%(x[np.argmax(np.abs(y))],np.max(np.abs(y))))

Command:
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src timeout 900 python3 f1.py 2>&1 | grep -E "EXACT|v="
```

### 4. submode="KPM" (the default dynamical correlator) silently swallows every unrecognized keyword, including the documented deconvolve= and n= that shipped examples actually pass

- **Severity**: low  **Class**: silent-wrong
- **Code**: `src/dmrgpy/kpmdmrg.py:56 (get_dynamical_correlator: `deconvolve` and `n` are in the signature but never referenced; `n` is immediately overwritten by the tuple unpack at line 78) and src/dmrgpy/kpmdmrg.py:85 (dynamical_correlator_moments's **kwargs, which absorbs and discards everything else)`
- **Status**: FIXED -- unknown keywords raise TypeError, matching CVM/CVM_explicit; the dead `deconvolve`/`n` parameters are gone from the signature.

**Expected**

`deconvolve` is advertised as a live knob: docs/user_guide.md:2533 lists it
among the kwargs forwarded to kpmdmrg.get_dynamical_correlator ("delta, kernel,
es, deconvolve, ..."), and the shipped example
examples/dynamical_correlator/dynamical_correlator_sharpen/main.py:39 calls
get_dynamical_correlator(..., deconvolve="pm", kernel="lorentz") and plots the
result under the legend "Extrapolation". It changes nothing at all -- the
visible difference in that example comes entirely from kernel="lorentz". `n`
(number of Chebyshev polynomials) is likewise in the signature and used by two
shipped examples (examples/staticcorrelators/dcf_haldane/main.py:20 `n=300`,
examples/dynamical_correlator/dynamical_correlator_entropy/main.py:18 `n=600`);
it is unpacked-over one line later and has zero effect. Beyond those two,
dynamical_correlator_moments's bare **kwargs means the default submode accepts
ANY keyword silently, so a user writing get_dynamical_correlator(...,
kpmmaxm=2) or (..., maxm=200) -- entirely plausible, given every other dmrgpy
tuning knob is set as an attribute of the same name -- gets the un-tuned result
with no diagnostic. Either the kwargs should be honored or unknown ones
rejected with a TypeError, the way submode="CVM_explicit" already does.

**Observed**

```
{'deconvolve': 'pm'}             accepted, max|y-base| = 0.000e+00
{'n': 5}                         accepted, max|y-base| = 0.000e+00
{'kpmmaxm': 2}                   accepted, max|y-base| = 0.000e+00
{'maxm': 200}                    accepted, max|y-base| = 0.000e+00
{'nsweeps': 1}                   accepted, max|y-base| = 0.000e+00
{'window': [0, 1]}               accepted, max|y-base| = 0.000e+00
{'totally_bogus_kwarg': 42}      accepted, max|y-base| = 0.000e+00
```

**Reviewer** (independent, brief was to refute)

Reproduced verbatim with my own script (v3 extension present,
cppext.available(3)=True): all seven keywords -- deconvolve="pm", n=5,
kpmmaxm=2, maxm=200, nsweeps=1, window=[0,1], totally_bogus_kwarg=42 -- were
accepted and gave max|y-base| = 0.000e+00 exactly. Not an unconverged-DMRG
artifact: the diffs are bit-identical zeros and the mechanism is read off the
source, not inferred. Defect, tightest form: kpmdmrg.get_dynamical_correlator
(the default submode="KPM") declares deconvolve= and n= but never uses them (n
is overwritten one line later by the dynamical_correlator_moments unpack at
kpmdmrg.py:78; deconvolve has no call site anywhere -- git shows its only ever
call, kpm.deconvolution, was already commented out before the file-based
backend was removed in 45d5999, and it was overwritten-by-file-read even then),
and dynamical_correlator_moments's bare **kwargs (line 85) silently discards
every other keyword. Sibling submodes are not like this: I confirmed CVM and
CVM_explicit both raise TypeError on the same bogus keyword. I tried hard to
read it as an intended legacy vestige and it does not hold:
docs/user_guide.md:2533 lists deconvolve among the kwargs "forwarded to
kpmdmrg.get_dynamical_correlator unchanged", and
examples/dynamical_correlator/dynamical_correlator_sharpen/main.py:39 plots
deconvolve="pm",kernel="lorentz" under an "Extrapolation" legend whose visible
difference I showed comes entirely from the kernel (kernel=lorentz shifts the
curve by 8.8e-02; adding deconvolve="pm" changes it by exactly 0). Nothing in
ROADMAP.md or CLAUDE.md's deliberately-reproduced-legacy-bug list (which covers
chain_session.h items only) covers this. Reach: the Hermitian KPM path on every
backend -- itensor_version 2/3/"python" through kpmdmrg.py, and "julia_live"
through the mirrored mpsjulialive/dynamics.py::_kpm_dynamical_correlator, which
has the same two dead params and the same trailing **kwargs sink (read, not
run) -- plus infinitechain.td_dynamical_correlator_window, which forwards
**kwargs straight into it. The non-Hermitian KPM route
(nonhermitian/kpm.py::dynamical_correlator_nhkpm) is unaffected: no **kwargs,
and its n= is genuinely used. Severity low rather than medium: no computed
number is ever wrong -- every returned spectrum is exactly correct for the
parameters actually in effect, the doc-endorsed dead kwargs reproduce the
default result bit-for-bit rather than a subtly perturbed one, and the real
tuning interface (attributes, delta) works. The sharpest harm is the typo sink:
kpmmaxm=2 passed as a kwarg changes nothing while setting the same-named
attribute changes the spectrum by 5.8e-01, so a plausible mis-set knob yields
an un-tuned result with no diagnostic. Weighed against low: the docs and a
shipped, mislabeled example do advertise deconvolve as live. Colour worth
noting -- kpm.deconvolution(mode="pm") returned nan when I called it by hand on
this data, so the knob is abandoned rather than dormant; wiring it up as-is
would not work either.

**Repro**

```
File /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/f2.py:

import numpy as np
from dmrgpy import spinchain
sc = spinchain.Spin_Chain([2]*6)
h = 0
for i in range(5):
    h = h + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h); sc.maxm = 20; sc.nsweeps = 8; sc.kpmmaxm = 20
name = (sc.Sz[2], sc.Sz[2]); es = np.linspace(-1,4,50)
base = sc.get_dynamical_correlator(name=name, es=es, delta=0.2)[1]
for kw in [dict(deconvolve="pm"), dict(n=5), dict(kpmmaxm=2), dict(maxm=200),
           dict(nsweeps=1), dict(window=[0,1]), dict(totally_bogus_kwarg=42)]:
    y = sc.get_dynamical_correlator(name=name, es=es, delta=0.2, **kw)[1]
    print("%-32s accepted, max|y-base| = %.3e" % (str(kw), np.max(np.abs(y-base))))

Command:
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src timeout 600 python3 f2.py 2>&1 | grep accepted
```

### 5. itensor_version=3 + set_conserved_sector: get_rdm() returns a garbage matrix (uninitialized memory), not the reduced density matrix

- **Severity**: medium  **Class**: silent-wrong
- **Code**: `src/dmrgpy/mpscpp3/bindings.cc:571 (with src/dmrgpy/mpscpp3/chain_session.h:1143)`
- **Status**: FIXED -- Chain::reduced_dm reads the elements explicitly by index instead of via rho.visit(), and bindings.cc size-checks before copying.

**Expected**

get_rdm(i=0) must return the one-site reduced density matrix [[0.5,0],[0,0.5]]
(Hermitian, trace 1) -- exactly what the same call returns on the same chain
without a sector, and what itensor_version="python" returns with the sector
set. ROADMAP.md section 2 lists "Bond entanglement entropy / reduced density
matrix" as implemented on v3, and CLAUDE.md documents conserved-sector mode as
a supported v3 feature. Mechanism: Chain::reduced_dm
(chain_session.h:1129-1145) flattens rho with rho.visit(collect), which under
QN-carrying indices only walks the *stored* (block-sparse) elements -- 2 of
them for a 2x2 diagonal rho -- while bindings.cc:571 does
std::copy(flat.begin(),flat.end(),arr.mutable_data()) into a freshly allocated
dim*dim py::array_t with no size check, so the copied values land in the wrong
positions and the remaining slots keep whatever was in the allocation. There is
no error, no warning, and no shape mismatch to notice.

**Observed**

```
Typical run:
no sector  rdm(0) = [[(0.5+0j), 0j], [0j, (0.5+0j)]]
Sz=0       rdm(0) = [[(0.5+0j), (0.5+0j)], [0j, (0.5+0j)]]
trace = (0.9999999999998932+0j)   hermitian?  False

A different run of the very same code (inside the earlier probe driver sec.py) gave instead:
rdm [[0.5+0.j 0.5+0.j]
 [0. +0.j 0. +0.j]]
i.e. trace 0.5. The (1,0) and (1,1) entries change between runs, which is the signature of uninitialized memory rather than a deterministic transposition. The pure-Python backend on the identical chain returns the correct [[0.5,0],[0,0.5]] both with and without the sector, and v3 without a sector is also correct.
```

**Reviewer** (independent, brief was to refute)

Confirmed, and the uninitialized-memory mechanism was proven directly rather
than inferred. Defect: mpscpp3's Chain::reduced_dm (chain_session.h:1129-1145)
flattens the one-site density matrix with rho.visit(collect), which under
conserved-sector mode (QN-carrying site indices) enumerates only the STORED
block-sparse elements -- for Sz-conserving spin sites, dim 1x1 blocks instead
of dim*dim entries. bindings.cc:568-571 then std::copy's that short vector into
a freshly allocated, uninitialized py::array_t({dim,dim}) with no size check,
so the values land at wrong linear offsets (packed into the first row) and the
remaining slots are never written, returning whatever was on the heap. Proof:
with S=1 sites and a 0.3*Sz[0] field the no-sector rdm is
diag(0.0687,0.26554,0.66575); under Sz=0, poisoning the heap with 9-element
complex arrays of 7.77+7.77j immediately before each get_rdm call (sector set
before any prior rdm call) returns [[0.0687,0.26554,0.66575],[7.77+7.77j
x3],[7.77+7.77j x3]], trace 15.61+15.54j -- the poison bytes appear verbatim.
The S=1/2 claim also reproduces exactly: [[0.5,0.5],[0,0.5]], non-Hermitian,
trace ~1 only by coincidence of stale freed memory. Not a repro flaw: sector
gs_energy is the exact -2.4935771338879 for n=6, the no-sector path on the same
chain is correct, itensor_version="python" is correct with and without the
sector, and the sector rdm's stored values equal the dense rdm's diagonal, so
the state is right and only the flatten/copy is broken. Not intended:
chain_session.h:4813's purpose-built reject_sector() guard covers the paths
that genuinely cannot work under a sector and deliberately excludes
correlator/AutoMPO paths; ROADMAP.md section 2 lists reduced density matrix as
implemented on v3; user_guide.md documents promote_to_dense()/promote_mps()
only for flux-violating operators, not as a prerequisite for reading an rdm;
and it is not one of CLAUDE.md's enumerated deliberately-reproduced legacy
bugs. Blast radius, narrow but silent: only Many_Body_Chain.get_rdm()
(densitymatrix.reduced_dm -> _session.reduced_dm) on itensor_version=3 with a
conserved sector set, wrong for every site dimension > 1. mpscpp2 has no sector
support at all, pyitensor and julia_live are unaffected. Entanglement entropy
is NOT affected (entropy.py uses _session.bond_entropy; pair/site entropy and
mutual information use the pure-Python reduced_dm_projective), and the other
std::copy(flat...) bindings (correlation_matrix, four_correlation_tensor*) size
their vectors explicitly and carry no analogous defect.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/rdm_repro.py:

import numpy as np
from dmrgpy import spinchain
n=6
sc=spinchain.Spin_Chain(["S=1/2"]*n)
h=0
for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
sc.setup_cpp(version=3); sc.maxm=30; sc.nsweeps=10
sc.set_hamiltonian(h)
print("no sector  rdm(0) =",np.round(np.array(sc.get_rdm(i=0)),6).tolist())
sc.set_conserved_sector(Sz=0)
r=np.array(sc.get_rdm(i=0))
print("Sz=0       rdm(0) =",np.round(r,6).tolist())
print("trace =",np.trace(r),"  hermitian? ",np.allclose(r,r.conj().T))

Command line (run it several times -- the answer is not stable):
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 rdm_repro.py
```

### 6. itensor_version=3 + set_conserved_sector: any operator product with a repeated site (e.g. the four-point correlator tensor) aborts the whole process with SIGABRT

- **Severity**: high  **Class**: crash
- **Code**: `src/dmrgpy/mpscpp3/chain_session.h:4577 (sector_terms)`
- **Status**: FIXED -- sector_terms now also composes each site's own factors: an identically-zero product drops the term (as ED answers), an impossible charge raises std::invalid_argument.

**Expected**

vev(Cdag_0 C_1 Cdag_0 C_1) must return 0 (the operator is identically zero
after anticommutation), and get_four_correlation_tensor must return the <Cdag_i
C_j Cdag_k C_l> tensor -- both are what itensor_version="python" returns in the
same sector, and what v3 returns outside a sector. ROADMAP.md section 2 marks
all three four-point ctmodes implemented on v3. The in-scope defect is the
sector guard, not vendored ITensor: CLAUDE.md states that sector mode exists
partly because "ITensor's own Error() calls abort() ... so nothing it raises is
catchable from Python -- the sector code throws
std::invalid_argument/std::runtime_error instead" and that "every terms -> MPO
path routes through mpo_from_terms/ampo_from_terms so the check covers
vev/correlators/KPM vertices too". sector_terms() (chain_session.h:4577-4602)
only validates each term's *net* charge, summed over its factors; a term whose
factors repeat a site out of order (Adag(1) A(2) Adag(1) A(2), net Nf=0) passes
the check, and AutoMPO then accumulates the partial flux left to right,
reaching QN Nf=-2 on the internal link and aborting. This is exactly the case
the task calls "the mitigation itself has a gap".

**Observed**

```
minab.py:
gs -2.2360679774997902
about to call vev(Cdag0 C1 Cdag0 C1) -- exact answer is 0
I = (dim=2|id=444|"l=2,Link") <Out>
  1: 2 QN()
Q = QN({"Nf",-2,-1})
From line 699, file index.cc

Index does not contain given QN block.

Index does not contain given QN block.
/bin/bash: line 31: 1617498 Aborted                 ... python3 minab.py
EXIT=134

four.py 3: identical abort, EXIT=134, before the first ctmode prints anything.
Controls that all succeed: `python3 four.py python` prints
explicit ok maxdiff_vs_ED=2.21e-09 / full ok maxdiff_vs_ED=2.21e-09 / sweep ok maxdiff_vs_ED=2.21e-09;
v3 with the sector cleared prints `vev = 0j`; pyitensor with the sector set prints `vev = 0j`.
```

**Reviewer** (independent, brief was to refute)

Reproduced exactly, with my own scripts. minab.py: gs -2.2360679774997902 in
sector, then "Index does not contain given QN block." with I =
(dim=2|...|"l=2,Link") / Q = QN({"Nf",-2,-1}), Aborted, EXIT=134. four.py 3:
same abort, EXIT=134, after "gs done" prints (i.e. the ground-state solve
succeeds; the abort is in the four-point call). Controls all pass: `four.py
python` prints explicit/full/sweep ok maxdiff_vs_ED=2.21e-09; v3 with no sector
prints explicit/full/sweep ok maxdiff_vs_ED=4.77e-12 and vev = 0j.  TIGHTEST
FORM (the claim as written, "any operator product with a repeated site", is
overbroad and I refuted that part): `sector_terms`
(src/dmrgpy/mpscpp3/chain_session.h:4577-4602) validates only each term's NET
charge, summed over its factors; it never checks the charge accumulated on a
single site. A net-conserving term whose same-site factors accumulate a charge
for which that site's QN-carrying index has no block passes the guard, reaches
AutoMPO, and aborts the whole process. Repeated sites per se are fine: I
verified in-sector v3 succeeds for Cdag0 C1 Cdag1 C0 and C1 Cdag0 Cdag0...
(0,1,1,0), (1,0,0,1), (0,0,1,1), (0,2,1,3), (0,1,1,2) all return the ED value;
a boson chain (maxnb=3) with Adag0 A1 Adag0 A1, a S=1 chain with S+0 S-1 S+0
S-1, and a native Hubbard pair-hop Cdagup0 Cup1 Cdagdn0 Cdn1 (0.0637694062 vs
ED 0.0637694062) all work in-sector. The crashing set is the terms that pile >1
quantum on a site whose index cannot carry it -- on JW spinless-fermion chains
(hardcore dim-2 "Adag"/"A" sites, per MultiOperator.to_terms()) that is exactly
the strings with the same operator twice at one site, e.g. (0,1,0,1) and
(0,1,2,1), which are the identically-zero strings.  That narrowness does not
reduce impact: get_four_correlation_tensor enumerates all (i,j,k,l)
(correlationentropy.py:354, Op=(Cdag[i]*C[j])*(Cdag[k]*C[l])), so it hits those
tuples unavoidably. Run per-ctmode in separate processes, ALL THREE modes kill
the process in-sector on v3: explicit EXIT=134 and full EXIT=134 with the
guard-gap signature ("Index does not contain given QN block"), and sweep
EXIT=134 with a DIFFERENT error ("Mismatched QN Index arrows") -- a distinct QN
bug in its own contraction path, which a fix confined to sector_terms would
likely not cure. The claimed locus explains explicit/full only.  Not intended,
not a documented gap, not misuse: docs/user_guide.md promises a non-conserving
operator "raises a ValueError naming it", set_conserved_sector's docstring
promises v3 and "python" give "the same answers",
CLAUDE.md/docs/documentation.md state the sector code exists precisely because
ITensor's Error() aborts uncatchably and that "every terms -> MPO path routes
through mpo_from_terms/ampo_from_terms so the check covers vev/correlators/KPM
vertices too", and ROADMAP §2 marks all three four-point ctmodes ✅ on v3.
Nothing in ROADMAP.md or the docs restricts four-point correlators in sector
mode. Nothing here is a convergence artifact (it is a crash, and the sector gs
energy is correct and matches the dense run), and PYTHONPATH pinned the working
tree.  BLAST RADIUS: itensor_version=3 only, sector mode only
(set_conserved_sector), DMRG only. Public calls affected: vev()/aMb() on any
operator with a per-site over-accumulating factor string, and
get_four_correlation_tensor() in every ctmode (explicit/full/sweep) on a JW
fermionic chain in a sector -- an uncatchable SIGABRT of the whole Python
process, no traceback, no try/except recovery, no in-sector workaround
(promote_to_dense() works only by leaving the sector). pyitensor and
v3-without-sector are unaffected and correct.

**Repro**

```
Minimal script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/minab.py:

import numpy as np
from dmrgpy import fermionchain
n=4
fc=fermionchain.Fermionic_Chain(n)
h=0
for i in range(n-1): h=h+fc.Cdag[i]*fc.C[i+1]+fc.Cdag[i+1]*fc.C[i]
fc.setup_cpp(version=3); fc.maxm=20; fc.nsweeps=8
fc.set_hamiltonian(h)
fc.set_conserved_sector(Nf=2)
print("gs",fc.gs_energy(),flush=True)
Op=fc.Cdag[0]*fc.C[1]*fc.Cdag[0]*fc.C[1]   # net Nf change = 0, exactly zero operator
print("about to call vev(Cdag0 C1 Cdag0 C1) -- exact answer is 0",flush=True)
print("vev =",fc.vev(Op),flush=True)

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src timeout 300 python3 minab.py; echo EXIT=$?

Same abort from the public four-point API, script four.py (the whole tensor, first ctmode tried):

import sys, numpy as np
from dmrgpy import fermionchain
iv=sys.argv[1]
if iv!="python": iv=int(iv)
n=6
fc=fermionchain.Fermionic_Chain(n)
h=0
for i in range(n-1): h=h+fc.Cdag[i]*fc.C[i+1]+fc.Cdag[i+1]*fc.C[i]
for i in range(n-1): h=h+0.5*fc.N[i]*fc.N[i+1]
if iv=="python": fc.setup_python()
else: fc.setup_cpp(version=iv)
fc.maxm=30; fc.nsweeps=10
fc.set_hamiltonian(h)
ref=np.array(fc.get_gs(mode="ED").get_four_correlation_tensor(ctmode="explicit"))
fc.set_conserved_sector(Nf=3)
fc.gs_energy()
wf=fc.get_gs()
for ct in ["explicit","full","sweep"]:
    try:
        t=np.array(wf.get_four_correlation_tensor(ctmode=ct))
        print(ct,"ok maxdiff_vs_ED=%.2e"%np.max(np.abs(t-ref)))
    except Exception as ex: print(ct,"RAISED",type(ex).__name__,str(ex)[:140])

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src timeout 400 python3 four.py 3; echo EXIT=$?
(and `python3 four.py python` as the control)
```

### 7. Changing maxm/nsweeps after the first gs_energy() is silently ignored: every later call returns the first, less-converged energy

- **Severity**: high  **Class**: silent-wrong
- **Code**: `src/dmrgpy/groundstate.py:117`
- **Status**: FIXED -- groundstate.gs_is_current() compares the solver parameters against the ones the stored state was computed under.

**Expected**

A fresh DMRG solve, i.e. the same numbers as the fresh-chain column. This is
not my inference about what is reasonable -- it is the contract the code and
docs state for themselves. groundstate.py's own comment above the send-cache
key says: "Every solver parameter that a re-run would pick up (maxm, nsweeps,
cutoff, noise, and the MPO bond dimension the Hamiltonian is built with) is
part of the key: the session's energy cache survives a skipped re-send, so a
user bumping any of these between bare gs_energy() calls must get a fresh
solve, not the cached energy computed under the old params."
docs/documentation.md:2024-2037 repeats it. The send-cache key really does
include maxm/nsweeps and the Hamiltonian really is re-sent -- but
gs_energy_single() sets self.skip_dmrg_gs = True at groundstate.py:117 at the
end of every call, and line 112 then passes skip_dmrg=True into
Chain::gs_energy (mpscpp3/chain_session.h:532: `if (skip_dmrg && have_wf0_ &&
have_wf0_energy_) return wf0_energy_;`, pyitensor/chain.py:488 identically), so
DMRG is never re-run and the stale wavefunction's energy comes back. Only
restart()/set_hamiltonian() clear the flag. The failure is maximally plausible-
looking: a bond-dimension convergence study -- the single most routine DMRG
workflow -- prints five identical numbers, which reads as 'converged at maxm=2'
rather than as a bug, and is wrong by 0.21 (5%) here. Every downstream
observable (vev, correlators, entropy) is then taken on the stale wf0 too.

**Observed**

```
itensor_version=3
  naive convergence sweep on ONE chain:
     maxm=  2 -> e0 = -4.0474179540
     maxm=  4 -> e0 = -4.0474179540
     maxm= 10 -> e0 = -4.0474179540
     maxm= 20 -> e0 = -4.0474179540
     maxm= 40 -> e0 = -4.0474179540
  same sweep, fresh chain each time (truth):
     maxm=  2 -> e0 = -4.0474179540
     maxm=  4 -> e0 = -4.2519330378
     maxm= 10 -> e0 = -4.2580284450
     maxm= 20 -> e0 = -4.2580352068
     maxm= 40 -> e0 = -4.2580352073
  same chain but calling restart() before each solve:
     maxm=  2 -> e0 = -4.0474179540
     maxm=  4 -> e0 = -4.2519330378
     maxm= 10 -> e0 = -4.2580253135
     maxm= 20 -> e0 = -4.2580352068
     maxm= 40 -> e0 = -4.2580352073
itensor_version='python'
  naive convergence sweep on ONE chain:
     maxm=  2 -> e0 = -4.0423998274
     maxm=  4 -> e0 = -4.0423998274
     maxm= 10 -> e0 = -4.0423998274
     maxm= 20 -> e0 = -4.0423998274
     maxm= 40 -> e0 = -4.0423998274
  same sweep, fresh chain each time (truth):
     maxm=  2 -> e0 = -4.0423998274
     maxm=  4 -> e0 = -4.2519337370
     maxm= 10 -> e0 = -4.2580284363
     maxm= 20 -> e0 = -4.2580352068
     maxm= 40 -> e0 = -4.2580352073
(itensor_version=2 behaves identically: -4.0474179540 five times on the reused chain, the correct ramp on fresh chains.)

A second probe in the same script family (scratchpad/b20.py) shows the same for nsweeps: 'iv=3: nsweeps=1 -> -4.2572184122 ; then nsweeps=20 -> -4.2572184122'.
```

**Reviewer** (independent, brief was to refute)

REPRODUCED exactly (my own copy of the script, all three backends). ED gives
-4.2580352073, matching the fresh-chain column, so the fresh-chain ramp is the
truth and the reused chain's -4.0474179540 at maxm=40 is wrong by 0.21 (4.95%),
not unconverged DMRG. Tightening nsweeps instead of maxm shows the same (iv=3:
nsweeps=1 -> -4.2580276842, then nsweeps=40 -> -4.2580276842, fresh nsweeps=40
-> -4.2580352073).  THE CLAIMED MECHANISM IS WRONG, THE BEHAVIOUR IS REAL.
Tightest statement of the defect: manybodychain.py:908, `def gs_energy(...): if
self.computed_gs: return self.e0`, returns before groundstate.gs_energy_single
is ever entered. Monkeypatching gs_energy_single proves it is entered exactly
ONCE across the whole five-point sweep. Therefore the elaborate send-cache key
in groundstate.py:87-107 -- which correctly includes maxm, nsweeps, cutoff,
noise, mpomaxm, ramp_key and sector_key -- is unreachable dead code for
repeated bare gs_energy() calls on one chain. The claimant's groundstate.py:117
/ skip_dmrg account is disproved: mpscpp3 set_hamiltonian delegates to
set_hamiltonian_mpo, which sets have_wf0_energy_=false, so a re-send forces a
real solve regardless of skip_dmrg. Confirmed directly: clearing computed_gs
alone while LEAVING skip_dmrg_gs=True yields the correct ramp (-4.047, -4.252,
-4.258, -4.2580352068, -4.2580352073); clearing skip_dmrg_gs alone while
leaving computed_gs=True changes nothing (five identical values).  NOT
INTENDED, by the codebase's own words. groundstate.py's send-cache comment ("a
user bumping any of these between bare gs_energy() calls must get a fresh
solve, not the cached energy computed under the old params"),
docs/documentation.md 4.4, and tests/test_bond_dimension_ramp.py's own
docstring ("without them, flipping bond_ramp between two bare gs_energy() calls
would hand back the previous schedule's cached energy, silently, forever") all
state the exact contract that is violated. That test has to set
sc.computed_gs=False by hand to reach the key at all -- direct evidence the
mechanism cannot be exercised through the public path it was written to
protect. The only place the computed_gs cache is documented as intended is the
narrow sector-promotion case (user_guide.md ~403), which explicitly names
restart() as the escape hatch and does not generalise to maxm/nsweeps, where
the code documents the opposite.  BLAST RADIUS. Every solver parameter in the
send-cache key, not just maxm/nsweeps: cutoff, noise, and bond_ramp too
(bond_ramp is precisely the case the existing test docstring worries about).
Backends 2, 3 and "python" verified directly; "julia_live" is gated by the same
line by construction, since manybodychain.gs_energy sits above backend
dispatch, but was not run per instructions. Downstream observables built on the
stale wf0 are wrong too: after bumping maxm 2 -> 40, vev(Sz[0]*Sz[5]) =
-0.0483633023 against the fresh chain's -0.0437166929 (11% off), and
get_gs()/correlators/entropy read the same stale wf0. Workarounds that do work
today: restart(), a fresh chain, calling set_hamiltonian again, or clearing
computed_gs. Nothing tells a user they are needed, which is why a routine bond-
dimension convergence study prints five identical numbers and reads as
"converged at maxm=2". The fix site is the computed_gs gate at
manybodychain.py:908, not groundstate.py:117.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/repro_maxm.py:

import numpy as np
from dmrgpy import spinchain
n=10
def build(iv):
    sc=spinchain.Spin_Chain(["S=1/2"]*n,itensor_version=iv)
    h=0
    for i in range(n-1):
        h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.nsweeps=10; sc.bond_ramp=False
    return sc
for iv in [3,"python",2]:
    sc=build(iv)
    print("itensor_version=%s"%repr(iv))
    print("  naive convergence sweep on ONE chain:")
    for m in [2,4,10,20,40]:
        sc.maxm=m
        print("     maxm=%3d -> e0 = %.10f"%(m,sc.gs_energy()))
    print("  same sweep, fresh chain each time (truth):")
    for m in [2,4,10,20,40]:
        s=build(iv); s.maxm=m
        print("     maxm=%3d -> e0 = %.10f"%(m,s.gs_energy()))
    print("  same chain but calling restart() before each solve:")
    sc2=build(iv)
    for m in [2,4,10,20,40]:
        sc2.maxm=m; sc2.restart()
        print("     maxm=%3d -> e0 = %.10f"%(m,sc2.gs_energy()))

Command:
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 -u repro_maxm.py
```

### 8. pyitensor silently treats an operator on an out-of-range site index as the identity, so a bad index yields a plausible wrong energy / vev=1.0 instead of an IndexError

- **Severity**: medium  **Class**: silent-wrong
- **Code**: `src/dmrgpy/pyitensor/autompo.py:87`
- **Status**: FIXED -- pyitensor/autompo.py range-checks every operator's site, like v2, v3 and ED already did.

**Expected**

An IndexError, exactly as itensor_version=2, itensor_version=3 and mode="ED"
all give for the identical input. Instead itensor_version="python" drops the
factor and keeps the coefficient: vev returns <psi|Id|psi> = 1.0 (a physically
plausible-looking number -- <Sz> is bounded by 0.5, so 1.0 is not obviously
garbage, and a user reading it as a magnetization sees nothing wrong), and a
Hamiltonian carrying 5.0*Sz on nonexistent site 10 returns exactly
-1.61602540... + 5.0 = 3.38397459..., i.e. the term became 5.0*identity, a
constant energy shift. Root cause: pyitensor/autompo.py's HTerm.resolve()
builds the per-site factor list with `for site in range(1, n + 1)` over
`by_site`, so any key outside 1..n is silently never consulted; nothing
anywhere on that path range-checks the site. Negative indices are the same trap
and are worse because they look supported: `sc.get_operator("Sz",-1)` gives
identity/1.0 while `sc.Sz[-1]` (plain Python list indexing) really is the last
site, so the two spellings of "the last site" quietly disagree. This is exactly
the failure mode a user hits when reusing an operator list built for a longer
chain, or when an index arithmetic bug (i+1 at the boundary) walks off the end.

**Observed**

```
iv='python' reference e0 = -1.6160254038
   vev(Sz[site 10]) -> (1-3.0662592690605672e-18j)
   vev(Sz[site -1]) -> (1-3.0662592690605672e-18j)
   gs_energy with 5.0*Sz[site 10] in H -> 3.38397459621556
   same, mode='ED' RAISED IndexError: list index out of range
iv=3        reference e0 = -1.6160254038
   vev(Sz[site 10]) RAISED IndexError: vector::_M_range_check: __n (which is 10) >= this->size() (which is 4)
   vev(Sz[site -1]) RAISED IndexError: vector::_M_range_check: __n (which is 18446744073709551615) >= this->size() (which is 4)
   gs_energy with 5.0*Sz[site 10] in H RAISED IndexError: vector::_M_range_check: __n (which is 10) >= this->size() (which is 4)
   same, mode='ED' RAISED IndexError: list index out of range
iv=2        reference e0 = -1.6160254038
   vev(Sz[site 10]) RAISED IndexError: ... (same as v3)
   gs_energy with 5.0*Sz[site 10] in H RAISED IndexError: ... (same as v3)
```

**Reviewer** (independent, brief was to refute)

Reproduced exactly with my own script at tighter parameters (maxm=40,
nsweeps=20): iv="python" gives e0=-1.6160254038,
vev(get_operator("Sz",10))=1.0, vev(get_operator("Sz",-1))=1.0, and gs_energy
with 5.0*Sz[site 10] added = 3.3839745962 (= e0 + 5.0 exactly).
itensor_version=2, itensor_version=3 and mode="ED" all raise IndexError on the
identical input.  Not a convergence artifact: e0 agrees across all three
backends to 1e-10, vev=1.0 is exact identity arithmetic (residual ~5e-17), and
the energy shift is exactly +5.000000000 to ~1e-13. No degeneracy or basis
dependence is involved. Import path verified to resolve to this checkout's
src/.  Not intended: I searched CLAUDE.md (including its explicit list of
deliberately-reproduced legacy bugs -- evoloperator's z^3/6, the
"moise"/"tevol_fit_td" key typos -- none of which is this), ROADMAP.md in full,
docs/documentation.md, docs/user_guide.md, and every pyitensor module
docstring. Nothing documents a site-index validation gap or any deliberate
leniency. The opposite: autompo.py's own docstring declares itself a
transcription of mpscpp3/ITensor's autompo.cc, whose real counterpart range-
checks through the sites vector (that is precisely the "vector::_M_range_check"
v3 emits); pyitensor/__init__.py's "not a general ITensor port" disclaimer
scopes API surface, not silently accepting input the mirrored code rejects; and
ROADMAP.md calls ED "the correctness reference every other backend is checked
against", which here raises.  Tightest form of the defect: HTerm.resolve()
(src/dmrgpy/pyitensor/autompo.py:87) builds its per-site factor list with `for
site in range(1, n+1)` over the `by_site` dict, so any key outside 1..n is
never consulted -- the out-of-range factor is silently dropped while the term's
coefficient survives -- and nothing on that path range-checks the site.  Reach:
itensor_version="python" only, finite chains, every term-driven public call,
since both mpobuilder.py:44 and tebd.py:96 route through resolve():
set_hamiltonian/gs_energy, vev, static correlators, KPM/TD dynamical
correlators, TEBD/TDVP evolution. (I did not test infinitechain; its L/C/R
indexing is a separate validation path, so I make no claim there.)  Two points
where my evidence is stronger than the claim's. (1) The claimant's plausibility
example is actually backwards -- <Sz>=1.0 exceeds the 0.5 bound and so is
detectable garbage. My fermion probe is the genuinely undetectable case: on a
4-site Fermionic_Chain, vev(get_operator("N",10)) returns 0.9999999999999998, a
perfectly in-bounds occupation number. (2) The dropped factor is not always a
clean identity: Jordan-Wigner runs in multioperator.to_terms() before
resolve(), so for a fermionic term the out-of-range operator and its out-of-
range F factors drop while the in-range F-string residue survives --
vev(Cdag[0]*C[10]) becomes Cdag[0] times a partial JW string, arbitrary garbage
rather than an identity substitution. I also confirmed a two-site correlator
silently collapses to a one-site vev (vev(Sz[0]*get_operator("Sz",10)) is bit-
identical to vev(Sz[0])), and the negative-index trap: get_operator("Sz",-1)
gives 1.0 while sc.Sz[-1] (plain list indexing) correctly gives the last site's
~0, so the two spellings of "the last site" silently disagree.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/repro_oob.py:

import numpy as np
from dmrgpy import spinchain
n = 4
def build(iv):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=iv)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    return sc, h
for iv in ["python", 3, 2]:
    sc, h = build(iv)
    sc.set_hamiltonian(h); sc.maxm = 20; sc.nsweeps = 8
    print("iv=%-8s reference e0 = %.10f" % (repr(iv), sc.gs_energy()))
    try:
        print("   vev(Sz[site 10]) ->", sc.vev(sc.get_operator("Sz", 10)))
    except Exception as e:
        print("   vev(Sz[site 10]) RAISED %s: %s" % (type(e).__name__, str(e)[:90]))
    try:
        print("   vev(Sz[site -1]) ->", sc.vev(sc.get_operator("Sz", -1)))
    except Exception as e:
        print("   vev(Sz[site -1]) RAISED %s: %s" % (type(e).__name__, str(e)[:90]))
    sc2, h2 = build(iv)
    h2 = h2 + 5.0*sc2.get_operator("Sz", 10)
    sc2.set_hamiltonian(h2); sc2.maxm = 20; sc2.nsweeps = 8
    try:
        print("   gs_energy with 5.0*Sz[site 10] in H ->", sc2.gs_energy())
    except Exception as e:
        print("   gs_energy with 5.0*Sz[site 10] in H RAISED %s: %s" % (type(e).__name__, str(e)[:90]))
    try:
        print("   same, mode='ED' ->", sc2.gs_energy(mode="ED"))
    except Exception as e:
        print("   same, mode='ED' RAISED %s: %s" % (type(e).__name__, str(e)[:90]))

Command:
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 -u repro_oob.py
```

## Confirmed findings, second verification pass

These came out of the same five lenses but exceeded the first run's reviewer budget; they were verified separately, under the same brief.

### 9. itensor_version=3 evolve_and_measure with tevol_method="TDVP" (the default) or "TDVP_GSE" silently renormalizes the evolving MPS to unit norm from the first step, so for a non-unit-norm start (e.g. evolution_ABA's A|gs>) C(0) is correct but every t>0 value, and the return_wf final state, is scaled by 1/||wf||^2

- **Severity**: high  **Class**: silent-wrong
- **Code**: `/u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/mpscpp3/chain_session.h:1364 (evolve_and_measure_tdvp) and :1431 (evolve_and_measure_tdvp_gse)`
- **Status**: FIXED -- both v3 methods capture the input norm and restore it after each normalizing tdvp_step, as quench_tdvp/quench_tebd always did.
- **Affects**: Only the compiled itensor_version=3 backend, and only via Chain::evolve_and_measure_tdvp (chain_session.h ~1364) and Chain::evolve_and_measure_tdvp_gse (~1431). Public API hit: timedependent.evolve_and_measure(mode="DMRG", wf=<non-unit-norm MPS>) and timedependent.evolution_ABA(mode="DMRG", A=<non-unitary operator>), with tevol_method in ("TDVP","TDVP_GSE") -- TDVP is the default and itensor_version=3 is the default -- plus tevol_method="AUTO" whenever TEBD rejects a non-nearest-neighbour H and AUTO falls back to TDVP. return_wf=True hands back a wrongly-normalized final MPS from the same calls. NOT affected: mode="ED", itensor_version=2, v3 with tevol_method="TEBD" or "MPO", itensor_version="python" (its _tdvp_step_fn does not normalize), every quench_*/dynamical-correlator path on all backends (they all do psi.normalize(); psi *= norm0), and julia_live (the identical defect is recorded as already fixed there, docs/documentation.md:3265 -- not re-run here per instructions). The wf=None default path (unit-norm wf0) is numerically unaffected.

**Expected**

<psi(t)|O|psi(t)> for the state actually handed in, i.e. proportional to
||psi||^2, which is what ED, itensor_version=2, itensor_version=3 with
tevol_method="TEBD"/"MPO", and itensor_version="python" all return (they track
the ED trajectory to 1e-4 or better and scale correctly with `scale`). The v3
TDVP/TDVP_GSE result is identical for scale=1.0, 0.5 and 3.0, i.e. the state is
silently renormalized at the first step; the answer is off by exactly
1/||psi||^2 (a factor 4 for A=Sx[0]). This is an oversight, not a convention:
the sibling methods quench_tdvp (chain_session.h ~1355), quench_tdvp_gse
(~1424) and quench_tebd (~1499) all follow their step with `psi1.normalize();
psi1 *= norm0;`, and CLAUDE.md documents a "manual norm-restoration workaround"
needed for the C++ METTS port for the same reason -- only
evolve_and_measure_tdvp and evolve_and_measure_tdvp_gse omit it. ROADMAP.md
marks TDVP and TDVP_GSE real-time evolution as implemented on v3, and both
itensor_version=3 and tevol_method="TDVP" are the defaults, so the wrong number
is what a default-configured user gets from evolution_ABA / evolve_and_measure
with any operator-perturbed start.

**Observed**

```
s20.py:
scale=1.0  ED  -0.23618   3/TDVP  -0.23618   3/TEBD  -0.23621   3/MPO  -0.23733   python/TDVP  -0.23618   2/MPO  -0.23733
scale=0.5  ED  -0.05905   3/TDVP  -0.23618   3/TEBD  -0.05905   3/MPO  -0.05933   python/TDVP  -0.05905   2/MPO  -0.05933
scale=3.0  ED  -2.12564   3/TDVP  -0.23618   3/TEBD  -2.12587   3/MPO  -2.13594   python/TDVP  -2.12564   2/MPO  -2.13594

f_tdvp_norm.py (public evolution_ABA path, A=Sx[0] so ||A|gs>||^2 = 0.25):
3        TDVP      [-0.059045 -0.235105 -0.231891 ...]  maxdev=0.176
3        TDVP_GSE  [-0.059045 -0.235105 -0.231891 ...]  maxdev=0.176
3        TEBD      [-0.059045 -0.058777 -0.057974 ...]  maxdev=1.57e-05
3        MPO       [-0.059045 -0.058775 -0.057963 ...]  maxdev=0.000252
2        MPO       [-0.059045 -0.058775 -0.057963 ...]  maxdev=0.000252
python   TDVP      [-0.059045 -0.05881  -0.058021 ...]  maxdev=0.000243
ED                 [-0.059045 -0.058775 -0.057963 ...]
```

**Reviewer** (independent, brief was to refute)

Reproduced independently with my own scripts
(/tmp/.../scratchpad/skeptic_norm.py, skeptic_scale.py), n=5 S=1/2 XXZ chain +
0.4*Sz field, maxm=20, nsweeps=14, PYTHONPATH pinned to the repo src
(cppext.available reports both v2 and v3 compiled, default 3).  evolution_ABA
with A=Sx[0], B=Sz[2] (||Sx|gs>||^2 = 0.25 exactly):   ED          [-0.059045
-0.058776 -0.057973 ...]   3/TDVP      [-0.059045 -0.235105 -0.231891 ...]
maxdev vs ED = 0.176   3/TDVP_GSE  [-0.059045 -0.235105 -0.231891 ...]  maxdev
vs ED = 0.176   3/TEBD 1.1e-05, 3/MPO 1.6e-04, 2/MPO 1.6e-04, python/TDVP
1.7e-04 Multiplying the v3 TDVP/TDVP_GSE t>0 points by 0.25 reproduces ED to
8.2e-11 / 7.8e-11. A scale sweep of evolve_and_measure with wf = scale*gs gives
ED -0.23618 / -0.05905 / -2.12564 for scale = 1.0 / 0.5 / 3.0, while v3/TDVP
returns -0.23618 for all three: the output is entirely scale-independent, i.e.
the input norm is discarded at the first step.  Mechanism:
mpscpp3/chain_session.h's tdvp_step() (~1564) calls tdvp(...) with
{"DoNormalize",true}, so every step forces unit norm. The three sibling drivers
in the same file -- quench_tdvp (~1355), quench_tdvp_gse (~1424), quench_tebd
(~1499) -- capture Cplx norm0 = sqrt(innerC(psi1,psi1)) at entry and follow
each step with psi1.normalize(); psi1 *= norm0;. evolve_and_measure_tdvp and
evolve_and_measure_tdvp_gse omit exactly that restoration;
evolve_and_measure_tebd needs none because tebd_step does not normalize.
Because the loop measures before evolving, out.correlator[0] still sees the
true norm, which is why C(0) is right and the trajectory is discontinuous after
one step.  Not intended, by the repo's own documentation. ROADMAP.md:109 says
julia's METTS "needs no manual norm-restoration workaround (unlike the C++
port), since this backend's tdvp_step already passes normalize=false" -- i.e.
the C++ port is known to require the restoration. docs/documentation.md:3265
records the identical defect being found and fixed in mpsjulialive/tdvp.jl's
evolve_and_measure_tdvp, described verbatim as "silently rescaling every result
by 1/||wf||^2 whenever the input wasn't already unit-norm (e.g.
evolution_ABA()'s wfA = A*wf, for any non-unitary A)", fixed by the norm0
pattern. docs/documentation.md:3871 records the same tdvp_step trap producing
"exactly a factor of 2 too large" results in dmrgtwotime.py. ROADMAP.md:81-82
mark v3 TDVP and TDVP_GSE real-time evolution as implemented, not restricted,
and this is not one of CLAUDE.md's deliberately reproduced legacy bugs
(evoloperator z^3/6, "moise", "tevol_fit_td"). (Minor correction to the claim:
the "manual norm-restoration workaround" note is in ROADMAP.md:109, not
CLAUDE.md.)  Repro flaws ruled out: not unconverged DMRG -- the deviation is an
exact, scale-independent multiplicative factor recovering ED to ~1e-10, which
truncation/sweep error cannot produce, and TEBD/MPO/python on the identical
chain and parameters agree with ED to 1e-4 or better; not degeneracy -- the
0.4*Sz field breaks the symmetry and every other backend reproduces ED
pointwise; not the known v3 <3-site abort -- n=5; not a shared convention -- ED
(edtk/timedependent.py's evolution_ABC applies A and never normalizes), v2, v3
TEBD/MPO and python all preserve the input norm, so only these two v3 methods
differ.  Why no test caught it: the only ED-vs-DMRG evolution_ABA comparison in
tests/test_time_evolution.py (line 116) runs setup_python() with
tevol_method="TEBD"; the v3 evolution_ABA tests either use TEBD, expect a
RuntimeError, or compare tevol_method="TDVP" against "AUTO" (which falls back
to TDVP), i.e. two identically-wrong trajectories against each other.
tests/test_julia_live.py does cover evolution_ABA norm preservation, but only
for julia_live.

**Where a fix goes**: In mpscpp3/chain_session.h, capture Cplx norm0 = sqrt(innerC(psi,psi)) at entry of both evolve_and_measure_tdvp and evolve_and_measure_tdvp_gse and add psi.normalize(); psi *= norm0; after each tdvp_step, mirroring quench_tdvp/quench_tdvp_gse and the recorded mpsjulialive fix. Fix at those two call sites only, not by flipping tdvp_step's DoNormalize, since tdz.py's complex-time path and metts_vev depend on its current semantics.

**Repro**

```
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 s20.py

--- s20.py ---
import numpy as np, warnings
from dmrgpy import spinchain, timedependent
n=5
def mk(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    if v=="python": sc.setup_python()
    elif v!="ED": sc.setup_cpp(v)
    sc.maxm=30; sc.nsweeps=12
    H=0
    for i in range(n-1):
        H = H + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+0.6*sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): H = H + 0.4*sc.Sz[i]
    sc.set_hamiltonian(H)
    return sc
nt=10; dt=0.1
for scale in [1.0,0.5,3.0]:
    sc=mk("ED"); wf=sc.get_gs(mode="ED")
    ts,ref = timedependent.evolve_and_measure(sc,mode="ED",nt=nt,dt=dt,operator=sc.Sz[2],wf=scale*wf)
    ref=np.array(ref).real
    row=["scale=%.1f  ED %9.5f"%(scale,ref[-1])]
    for ver,tm in [(3,"TDVP"),(3,"TEBD"),(3,"MPO"),("python","TDVP"),(2,"MPO")]:
        sc=mk(ver); sc.tevol_method=tm; w=sc.get_gs()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t2,y=timedependent.evolve_and_measure(sc,mode="DMRG",nt=nt,dt=dt,operator=sc.Sz[2],wf=scale*w)
        row.append("%s/%s %9.5f"%(ver,tm,np.array(y).real[-1]))
    print("   ".join(row))

A second repro through the fully public entry point evolution_ABA (same command line, file f_tdvp_norm.py):

import numpy as np, warnings
from dmrgpy import spinchain, timedependent
n=5
def mk(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    if v=="python": sc.setup_python()
    elif v!="ED": sc.setup_cpp(v)
    sc.maxm=30; sc.nsweeps=12
    H=0
    for i in range(n-1):
        H = H + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+0.6*sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): H = H + 0.4*sc.Sz[i]
    sc.set_hamiltonian(H)
    return sc
nt=10; dt=0.1
sc=mk("ED"); sc.get_gs(mode="ED")
ts,ref = timedependent.evolution_ABA(sc,A=sc.Sx[0],B=sc.Sz[2],mode="ED",nt=nt,dt=dt)
ref=np.array(ref).real
print("ED       ", np.round(ref,6))
for ver,tm in [(3,"TDVP"),(3,"TDVP_GSE"),(3,"TEBD"),(3,"MPO"),(2,"MPO"),("python","TDVP")]:
    sc=mk(ver); sc.tevol_method=tm; sc.get_gs()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t2,y = timedependent.evolution_ABA(sc,A=sc.Sx[0],B=sc.Sz[2],mode="DMRG",nt=nt,dt=dt)
    y=np.array(y).real
    print("%-8s %-9s"%(ver,tm), np.round(y,6), " maxdev=%.3g"%np.max(np.abs(y-ref)))
```

### 10. _fourier_transform_correlator omits the documented 1/pi, so submode="TD" returns pi*S_AB(omega) + Re C(0)*dt/2 instead of S_AB(omega); submode="TDZ" inherits the same pi (its further ~2x is its own reconstruction, not the FFT)

- **Severity**: high  **Class**: silent-wrong
- **Code**: `/u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/timedependent.py:428 (_fourier_transform_correlator)`
- **Status**: FIXED -- the transform applies the 1/pi convention and a trapezoidal endpoint weight. TDZ's own residual factor of 2 lives in its reconstruction and is NOT addressed here.
- **Affects**: Backend-independent: the defect is in shared Python post-processing, so it hits every path that reaches `timedependent._fourier_transform_correlator`. Measured identical (integral 1.28391/1.28392, peak 1.5208/1.5224) on `itensor_version=2`, `3`, `"python"` and on `mode="ED"` (edtk/dynamics.py routes `submode="TD"` back into `timedependent.dynamical_correlator`). Public API calls hit: `Many_Body_Chain.get_dynamical_correlator(submode="TD")` (mode="DMRG" and mode="ED"), `get_dynamical_correlator(submode="TDZ")` (tdz.py:305 reuses the same tail), and structurally — same function, no compensating factor anywhere downstream, verified by grepping tdz.py / infinitechain.py / pyitensor/idmrg_window.py for a stray `np.pi` — `timedependent.sxt_to_skomega` and `Infinite_Many_Body_Chain.td_dynamical_correlator_komega` (infinitechain.py:1406). `itensor_version="julia_live"` inherits it too (mpsjulialive/dynamics.py explicitly forwards TD/TDZ into the same top-level timedependent.py/tdz.py), but I did not execute it per instructions.

**Expected**

All submodes should return the same spectral function. docs/user_guide.md:1331
documents submode="TD" explicitly as S_AB(omega) = (1/pi) Re \int_0^T dt e^{i
omega t} C(t) w_delta(t), and timedependent.py's own dynamical_correlator
docstring says the FFT is "normalized as a Riemann sum (factor dtnew) to match
the analytic Fourier-transform convention of the other submodes". The code at
timedependent.py:428 computes `ss = np.fft.fft(cs)*dtnew` with no 1/pi at all,
so TD returns pi times the documented quantity plus a constant DC background.
Measured directly at converged nt (s16.py): td/pi reproduces the exact Lehmann
curve at the peaks (0.39557 vs 0.39159 at w=0.60; 0.25100 vs 0.24733 at w=1.40)
but carries a constant ~0.0043 background in the tails where the exact answer
is ~0.0003, which is what pushes the integral from the expected 0.25 to 1.284.
TDZ, which reuses the same _fourier_transform_correlator, is off by a further
factor (peak ratio ~8.5, not pi) that I did not diagnose. The existing tests
only check TD's peak position, never its amplitude, which is why this survives.
Peak position is also quantized to the FFT's 2*pi/T grid (0.630 vs the exact
0.660 at default nt), since `factor` oversamples time without extending T and
so cannot improve frequency resolution.

**Observed**

```
exact sum rule  int S(w) dw = <Sz0 Sz0> = 0.24999999999999997
ED    ED    integral=0.249202  peak=0.5267 at w=0.660
ED    KPM   integral=0.250000  peak=1.6881 at w=0.660
ED    CVM   integral=0.249202  peak=0.5267 at w=0.660
DMRG  KPM   integral=0.250000  peak=1.0305 at w=0.660
DMRG  EX    integral=0.249202  peak=0.5267 at w=0.660
DMRG  TD    integral=1.283968  peak=1.5208 at w=0.630
DMRG  TDZ   integral=1.571275  peak=4.4914 at w=0.630
```

**Reviewer** (independent, brief was to refute)

Confirmed on my own script (n=4 S=1/2 Heisenberg + 0.3*Sz field, maxm=20,
nsweeps=12, delta=0.1, es on [-20,20]). Exact sum rule int S(w) dw = <Sz0 Sz0>
= 0.25. ED/KPM/CVM/EX all return real arrays with integral 0.2492-0.2500 and
peak 0.5267 at w=0.660; submode="TD" returns integral 1.2839, peak 1.5208;
submode="TDZ" returns 1.5710, peak 4.4914.  The TD deviation decomposes
exactly. Measured at dt=0.1/0.05/0.025 against a hand-built Lehmann reference
(ED eigenvectors, sum_n |<n|Sz0|GS>|^2 * Lorentzian): the far-field residual y
- pi*S_exact for |w|>15 is 0.012518 / 0.006252 / 0.003124, i.e. exactly Re
C(0)*dt/2 = 0.25*dt/2 to 4 digits at every dt; after subtracting that constant
and dividing by pi, the integral is 0.24929/0.24932/0.24932 against the exact
0.24920. So the returned object is pi*S_AB(w) + Re C(0)*dt/2, i.e. i*G_AB(w)
plus an O(dt) real DC offset (the offset is the Riemann-vs-trapezoid endpoint
error: `np.fft.fft(cs)*dtnew` gives the t=0 sample full weight dt instead of
dt/2). Peak positions are unaffected apart from the 2*pi/T FFT grid
quantization (0.630 vs 0.660).  This is a bug, not a convention, on four
independent grounds: (a) docs/user_guide.md:1331 documents TD as S_AB(w) =
(1/pi) Re int_0^T dt e^{iwt} C(t) w_delta(t); (b) docs/user_guide.md:1002/1010
fix the codebase-wide convention S_AB = -(1/pi) Im G_AB and state that "All
five submodes below compute G_AB(omega) (or equivalently S_AB(omega))"; (c) the
code's own docstring (timedependent.py:289) says the Riemann normalization
exists "to match the analytic Fourier-transform convention of the other
submodes"; (d) empirically ED/CVM/KPM/EX return purely real arrays satisfying
the sum rule, while TD returns a complex array whose real part is pi times
theirs.  Ruled out as repro flaws: not DMRG convergence (identical to 5 digits
on v2/v3/python and on exact mode="ED"); not degeneracy (the 0.3*Sz field lifts
it, and ED/CVM agree to 7e-15); not finite-time/dt (the pi survives dt scan;
only the additive 0.0125 term scales with dt as predicted); not chain length
(n=4, well above v3's 3-site floor); not "quantities never meant to be equal"
(see (b)). ROADMAP.md:97 lists TD as fully implemented on all three backends,
and none of CLAUDE.md's deliberately-reproduced legacy bugs (evoloperator
z^3/6, "moise", "tevol_fit_td") touches this code.  Why it survived the suite:
tests/test_time_evolution.py and tests/test_dynamical_correlator.py only assert
TD's peak position, its FWHM ordering, on/off consistency of the realify
rewrite, and dt-independence of `y.imag` — the pi scales everything uniformly
and the endpoint offset is purely real, so no existing assertion can see
either. examples/.../dynamical_correlator_VS_ED checks only CVM and KPM against
ED; examples/.../dynamical_correlator_time_evolution plots `np.abs(y2)` for TD
against `y1.real` for KPM, which is exactly where a 3x amplitude mismatch would
have been visible and was not caught.  TDZ caveat (narrowing the lead's claim):
after removing pi, TDZ's integrated weight is 0.50016 at dt=0.1 and 0.50016 at
dt=0.05 — a further factor of exactly 2. That extra factor is NOT in the shared
FFT: TDZ's time grid is the same uniform 0,0.1,...,5.9 and its reconstructed
C(0) is 0.25, identical to TD's, but |C(t)| at the tail is ~2.2x TD's (0.301 vs
0.138), i.e. the n_max=4 perturbative real-axis reconstruction in tdz.py grows
rather than decays at the default alpha0=0.1. That is an accuracy problem of
tdz.py's reconstruction, plausibly a separate issue, and I did not diagnose it
further. The certain, shared defect is the missing 1/pi.

**Where a fix goes**: In `timedependent._fourier_transform_correlator`, scale the transform by 1/pi and half-weight the t=0 (and t=T) samples before the FFT — `ss = np.fft.fft(cs)*dtnew/np.pi` with `cs[0] *= 0.5` (or an equivalent trapezoidal correction) — so every caller (TD, TDZ, sxt_to_skomega, the infinite-chain S(k,w) reduction) inherits it in one place. Scale the complex return value rather than projecting onto its real part, since tests/test_time_evolution.py::test_td_dynamical_correlator_is_dt_independent consumes `y.imag`; TDZ's residual factor 2 belongs in tdz.py's reconstruction and should be fixed separately.

**Repro**

```
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 f_td_sumrule.py

--- f_td_sumrule.py ---
import numpy as np, warnings
from dmrgpy import spinchain
n=4
es = np.linspace(-20.,20.,4001); delta=0.1
def mk(v=3):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    if v!="ED": sc.setup_cpp(v)
    sc.maxm=20; sc.nsweeps=12
    H=0
    for i in range(n-1):
        H = H + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): H = H + 0.3*sc.Sz[i]
    sc.set_hamiltonian(H)
    return sc
sc=mk()
print("exact sum rule  int S(w) dw = <Sz0 Sz0> =", sc.vev(sc.Sz[0]*sc.Sz[0],mode="ED").real)
for m,sub,kw in [("ED","ED",{}),("ED","KPM",{}),("ED","CVM",{}),("DMRG","KPM",{}),
                 ("DMRG","EX",{}),("DMRG","TD",dict(predict=False)),("DMRG","TDZ",{})]:
    sc=mk(); nm=(sc.Sz[0],sc.Sz[0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        (x,y)=sc.get_dynamical_correlator(mode=m,submode=sub,name=nm,es=es,delta=delta,**kw)
    y=np.array(y).real
    print("%-5s %-5s integral=%.6f  peak=%.4f at w=%.3f"%(m,sub,np.trapezoid(y,es),y.max(),es[np.argmax(y)]))

(The same numbers reproduce with itensor_version=2 and "python": s10.py/s11.py in the same directory show TD integral -0.733 vs the exact -0.215 for the off-diagonal pair (Sz0,Sz1) on all three backends and on mode="ED" as well. Convergence check s15.py shows the integral stays 1.2839 for nt=600/1200/3000/6000/12000, so this is not a finite-time/resolution artifact.)
```

### 11. In conserved-sector mode, get_correlation_matrix()'s hardcoded default dmmode="fast" always raises ValueError on both sector-capable backends, even though dmmode="explicit"/"full" compute the same fully charge-conserving matrix correctly in the same call

- **Severity**: medium  **Class**: missing-combo
- **Code**: `src/dmrgpy/entropytk/correlationentropy.py:56 (dmmode="fast" default -> correlation_matrix_fast at :137)`
- **Status**: FIXED -- dmmode now defaults against the chain's state instead of being hardcoded to "fast".
- **Affects**: itensor_version=3 (compiled mpscpp3) and itensor_version="python" (pyitensor) -- i.e. both backends that implement set_conserved_sector at all. Public API hit: Many_Body_Chain.get_correlation_matrix() (manybodychain.py:754) and everything layered on it in entanglement.py -- get_correlation_entropy, get_correlation_entropy_density, get_correlated_orbitals, get_correlated_density, plus MPS.get_correlation_entropy_from_wf. Applies to any fermionic chain with any sector (Nf, Sz, or both), since a bare C_i changes every conserved charge. Not applicable to v2 or julia_live (no sector support); ED is unreachable because a sector-mode chain deliberately refuses to fall back.

**Expected**

The observable itself, <C_i^dag C_j>, conserves Nf exactly, and the same
function proves it is computable in the sector: dmmode="explicit" and
dmmode="full" reproduce the ED reference to 1e-12 (v3) / 2e-9 (pyitensor). Only
the default path breaks, because correlation_matrix_fast
(correlationentropy.py:137) materialises the intermediate states C_i|psi>, each
of which leaves the sector, so the sector guard rejects a calculation whose
result is perfectly legal. So `sc.get_correlation_matrix()` -- and everything
built on it, get_correlation_eigenvalues/get_correlation_entropy/get_correlated
_orbitals/get_correlated_density -- is unusable in a sector on BOTH sector-
capable backends unless the caller happens to know to pass an undocumented
dmmode. ROADMAP.md section 2 marks the all-pairs correlation matrix implemented
on both. Either the default should route to a sector-safe mode automatically,
or the error should name dmmode as the fix.

**Observed**

```
== itensor_version=3
default  RAISED ValueError: Conserved-sector mode (Nf=3): the term A(1) changes the conserved quantum numbers by Nf=-1 -- every operator b
explicit  ok  maxdiff=1.51e-12
full      ok  maxdiff=1.51e-12
== itensor_version=python
default  RAISED ValueError: conserved-sector mode (Nf=3): the term A(1) changes the conserved quantum numbers by Nf=-1 -- every operator b
explicit  ok  maxdiff=2.21e-09
full      ok  maxdiff=2.21e-09
```

**Reviewer** (independent, brief was to refute)

Real, and it is a dispatch/usability defect rather than a wrong-number defect.
Reproduced independently at n=6, maxm=20, nsweeps=12 on a t-V spinless chain
with set_conserved_sector(Nf=3): the default call raises ValueError ("the term
A(1) changes the conserved quantum numbers by Nf=-1") on itensor_version=3 and
on "python", while dmmode="explicit" and dmmode="full" return the matrix in the
same session, agreeing with the ED reference to 9.0e-13 (v3) and 1.0e-09
(python). All four downstream public methods raise the same error. The
observable itself, <C_i^dag C_j>, conserves every sector charge exactly; only
correlation_matrix_fast's intermediate states C_i|psi> leave the sector, so the
(correct, deliberate) sector guard in mpscpp3/chain_session.h::sector_terms and
pyitensor/sector.py rejects a calculation whose result is perfectly legal.
Repro flaws ruled out: it is an exception, not a numeric discrepancy, so
convergence and degeneracy are irrelevant; the ED reference is valid because
the global ground state lies in the requested sector (global ED e0 =
-3.222056822921109 with total <N> = 3.0000 vs sector e0 = -3.222056822921105);
n=6 is well clear of v3's <3-site abort; vev(Cdag[0]*C[1]) works fine in-
sector, so the sector solve is healthy. Intended-behaviour defenses examined
and rejected: ROADMAP.md has no sector entry at all, so nothing marks this
combination unimplemented; docs/user_guide.md's sector section enumerates
exactly what does not work under a sector (METTS, iDMRG/VUMPS, TEBD on
itensor_version=3) and states that expectation values, static correlators and
entanglement entropies do work -- the correlation-matrix family is neither
excluded nor promised; the documented "every operator built on the chain must
conserve" rule is about operators the caller supplies, and here the caller
supplies none. dmmode is documented nowhere in docs/user_guide.{md,tex} (only
in docs/documentation.md's internal notes and one example script), so the
workaround is not discoverable from the user-facing docs, and the error message
does not name it. Two existing precedents show the project's own convention is
to auto-select an available mode in exactly this situation:
_four_correlation_tensor_default_ctmode() picks the best available ctmode per
backend, and mpsjulialive/mps.py:64,69,109 hardcodes dmmode="explicit" because
julia_live's OpSum cannot build a parity-odd bare-C MPO either (documented at
docs/documentation.md:3138). Blast radius is limited by the fact that it fails
loudly and a one-kwarg workaround exists; no wrong numbers are ever returned.
One caveat for any fix: basis="Nambu" genuinely contains non-conserving <C_i
C_j> entries and must keep raising in a sector. Incidental and out of scope: on
the same sector chain, wf.get_four_correlation_tensor() (default
ctmode="sweep", v3) hard-aborts the process with ITensor's "Mismatched QN Index
arrows" (itensor.cc:912) instead of raising -- a separate, harder failure in
the same area.

**Where a fix goes**: In get_correlation_matrix_zeroT (entropytk/correlationentropy.py:56), make dmmode default to None and resolve it against the chain's state -- when a sector is set and basis!="Nambu", pick "full" (or "explicit") instead of "fast", mirroring _four_correlation_tensor_default_ctmode()'s availability-based auto-select and mpsjulialive/mps.py's forced dmmode="explicit"; an explicit dmmode= stays a hard request. Failing that, have the sector guard's ValueError name dmmode="explicit"/"full" as the remedy, and document dmmode in docs/user_guide.{md,tex}.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/f3_corrmat.py:

import sys, numpy as np
from dmrgpy import fermionchain
iv = sys.argv[1]
n = 6
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1): h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
for i in range(n-1): h = h + 0.5*fc.N[i]*fc.N[i+1]
if iv=="python": fc.setup_python()
else: fc.setup_cpp(version=3)
fc.maxm = 30; fc.nsweeps = 10
fc.set_hamiltonian(h)
ref = fc.get_correlation_matrix(mode="ED")   # exact reference, n=6
fc.set_conserved_sector(Nf=3)
fc.gs_energy()
try:
    got = np.array(fc.get_correlation_matrix())   # dmmode defaults to "fast"
    print("default  ok  maxdiff=%.2e"%np.max(np.abs(ref-got)))
except Exception as ex:
    print("default  RAISED %s: %s"%(type(ex).__name__,str(ex)[:110]))
for dm in ["explicit","full"]:
    got = np.array(fc.get_correlation_matrix(dmmode=dm))
    print("%-9s ok  maxdiff=%.2e"%(dm,np.max(np.abs(ref-got))))

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 f3_corrmat.py 3
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 f3_corrmat.py python
```

### 12. get_dynamical_correlator's documented string name= (e.g. name="ZZ", i=, j=) is accepted only for SPIN operators on submode TD/CVM/TDZ under mode="DMRG"; the default submode="KPM", submode="EX", submode="ROOTN" and every mode="ED" path reject it with an uninformative bare `raise` ("RuntimeError: No active exception to reraise"), which breaks two shipped examples outright

- **Severity**: medium  **Class**: missing-combo
- **Code**: `src/dmrgpy/kpmdmrg.py:104 (`if type(name[0])!=multioperator.MultiOperator: raise` -- a bare raise, so the caller sees "RuntimeError: No active exception to reraise"); src/dmrgpy/dcex.py:57 (name[0].get_dagger() on a str); src/dmrgpy/edtk/dynamics.py:22 (EDOperator(name[0],...)). The working paths go through src/dmrgpy/operatornames.py:85 str2MO, which does implement the string form.`
- **Status**: FIXED -- resolved once in Many_Body_Chain.get_dynamical_correlator, and str2MO covers every operator family recognize() returns.
- **Affects**: Public API: Many_Body_Chain.get_dynamical_correlator(name=<str>, i=, j=) for mode="DMRG" with submode in {"KPM" (the default), "EX", "ROOTN"}, and for mode="ED" with EVERY submode (the ED failure is upstream of submode dispatch). Also get_dynamical_correlator_moments(), which shares kpmdmrg.dynamical_correlator_moments. Working subset: submode TD/CVM/TDZ under mode="DMRG", and only for spin names (XX/YY/ZZ/+-/...); fermionic string names ("cdc", "densitydensity", ...) fail even there, because operatornames.str2MO's inner f() only handles Sx/Sy/Sz/Sp/Sm. Tested directly on itensor_version=3 (spin + fermionic chains) plus the ED path. Not run on v2 or "python", but both failing guards (kpmdmrg.py:100, rootndmrg.py:49, dcex.py:57, edtk/edchain.py:254) sit in the backend-agnostic Python layer upstream of backend selection, so v2/v3/"python" are hit identically by inspection. itensor_version="julia_live" was not run per instructions (its mpsjulialive/timedependent.py does carry its own str2MO call, so its TD path plausibly behaves like the C++ TD one, unverified).

**Expected**

docs/user_guide.md:1919-1921 states that name follows "the same (A,B)
convention as get_dynamical_correlator's own name= (a string like \"ZZ\", or an
explicit (MultiOperator,MultiOperator) tuple/list)".
operatornames.str2MO(self,name,i=0,j=0) implements exactly that, and submodes
TD/CVM/TDZ call it, so the string form demonstrably works there. The default
submode (KPM), submode="EX", and every mode="ED" route bypass str2MO and
require a MultiOperator pair, failing with an uninformative bare-raise
("RuntimeError: No active exception to reraise") or an AttributeError.
Consequences: (1) the shipped example examples/dynamical_correlator/dynamical_c
orrelator_spatial_inhomogeneous/main.py:35 calls
sc.get_dynamical_correlator(mode="DMRG", i=i, j=i, name="ZZ", ...) with the
default submode and therefore cannot run at all; (2) there is no way to cross-
check a submode="TD" string-form call against the ED oracle, since mode="ED"
rejects the same arguments the DMRG path accepts. Either every submode should
route operator resolution through str2MO, or KPM/EX/ED should reject a string
with a message naming the actual problem.

**Observed**

```
submode=KPM           mode=DMRG FAIL RuntimeError: No active exception to reraise
submode=KPM           mode=ED   FAIL RuntimeError: No active exception to reraise
submode=TD            mode=DMRG OK   max|y|=0.40015
submode=TD            mode=ED   FAIL RuntimeError: No active exception to reraise
submode=CVM           mode=DMRG OK   max|y|=0.11784
submode=CVM           mode=ED   FAIL RuntimeError: No active exception to reraise
submode=TDZ           mode=DMRG OK   max|y|=0.89619
submode=TDZ           mode=ED   FAIL RuntimeError: No active exception to reraise
submode=EX            mode=DMRG FAIL AttributeError: 'str' object has no attribute 'get_dagger'
submode=EX            mode=ED   FAIL RuntimeError: No active exception to reraise
```

**Reviewer** (independent, brief was to refute)

Reproduced exactly, with a control the lead did not have. On a 6-site
Heisenberg chain (maxm=20, nsweeps=8), name="ZZ", i=2, j=2: KPM/DMRG, KPM/ED,
TD/ED, CVM/ED, TDZ/ED, EX/DMRG, EX/ED, ROOTN/DMRG and ROOTN/ED all fail, while
TD/CVM/TDZ under mode="DMRG" succeed. The control run with
name=(sc.Sz[2],sc.Sz[2]) succeeds on all of those same paths and returns
numerically identical values to the string form where the string form works (TD
0.40015, CVM 0.11784, TDZ 0.89619 both ways), so this is purely an argument-
resolution gap, not unconverged DMRG, degeneracy, or two quantities that were
never meant to match. Tracebacks pin it to four bare `raise` statements:
kpmdmrg.py:100 and rootndmrg.py:49 (`if
type(name[0])!=multioperator.MultiOperator: raise`), edtk/edchain.py:254
(EDOperator's type check, reached from edtk/dynamics.py:23), and dcex.py:57
(`name[0].get_dagger()` on a str -> AttributeError). The working paths all
route through operatornames.str2MO (timedependent.py:132, cvm.py:22,
tdz.py:298), which implements the documented string form; the broken ones
bypass it. I could not construct an "intended behavior" defense: ROADMAP.md §4
lists every dynamical-correlator submode per backend and marks no restriction
on name= at all; this is not among CLAUDE.md's deliberately-reproduced legacy
bugs (evoloperator's z^3/6, the "moise" key, the unreachable "tevol_fit_td"
branch); docs/user_guide.md:1920 and docs/user_guide.tex:2190 explicitly state
that name= is "a string like \"ZZ\", or an explicit
(MultiOperator,MultiOperator) tuple/list" and present that as
get_dynamical_correlator's OWN convention; and two shipped examples use exactly
the broken combinations and cannot execute -- I ran clean copies of both
outside the repo with an explicit PYTHONPATH: examples/dynamical_correlator/dyn
amical_correlator_spatial_inhomogeneous/main.py (mode="DMRG", default submode,
name="ZZ") dies at kpmdmrg.py:100, and
examples/dynamical_correlator/finite_temperature_dynamical_correlator/main.py
(mode="ED", name="ZZ") dies at edchain.py:254. A shipped example that cannot
run is not a documented deliberate choice. One nuance supporting the "wrong
error message" half of the claim: the ED failure occurs at EDOperator(name[0])
before submode dispatch, so submode="TDZ", mode="ED" with a string reports the
bare-raise RuntimeError instead of the clean "NotImplementedError:
submode='TDZ' has no ED implementation" that the tuple form correctly produces.
Honest counterweight on severity: this never yields a wrong number -- every
failure is loud and immediate. The blast radius is a documented API form that
does not work on the DEFAULT submode, no ED oracle for string-form calls (so a
TD string call cannot be cross-checked against ED), two dead shipped examples,
and an error message that names no cause.

**Where a fix goes**: Normalize the operator pair once at the two dispatch points -- dynamics.get_dynamical_correlator and edtk/dynamics.get_dynamical_correlator -- by calling operatornames.str2MO(self, name, i=i, j=j) before submode dispatch, and extend str2MO's inner f() to fermionic names via name2MO so the string form is not silently spin-only. At minimum, replace the four bare `raise`s (kpmdmrg.py:100, rootndmrg.py:49, dcex.py:57, edtk/edchain.py:254) with a TypeError naming the accepted forms.

**Repro**

```
File /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/f3.py:

import numpy as np, io, contextlib
from dmrgpy import spinchain
sc = spinchain.Spin_Chain([2]*6)
h = 0
for i in range(5):
    h = h + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h); sc.maxm = 20; sc.nsweeps = 8; sc.kpmmaxm = 20
es = np.linspace(-1,4,20)
for sub in ["KPM","TD","CVM","TDZ","EX","ROOTN"]:
    for mo in ["DMRG","ED"]:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                x,y = sc.get_dynamical_correlator(mode=mo, submode=sub,
                        name="ZZ", i=2, j=2, es=es, delta=0.3)
            print("submode=%-13s mode=%-4s OK   max|y|=%.5f"%(sub,mo,np.max(np.abs(y))))
        except Exception as e:
            print("submode=%-13s mode=%-4s FAIL %s: %s"%(sub,mo,type(e).__name__,str(e)[:60]))

Command:
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src timeout 900 python3 f3.py 2>&1 | grep "^submode"
```

### 13. Thermal_Spin_Chain.__init__ drops itensor_version (the one constructor kwarg Many_Body_Chain actually honors), silently maps T<0 onto the documented T=0 branch, and raises UnboundLocalError at T==1e-5 exactly

- **Severity**: medium  **Class**: silent-wrong
- **Code**: `/u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/mpscpp3/chain_session.h:1364 (evolve_and_measure_tdvp) and :1431 (evolve_and_measure_tdvp_gse)`
- **Status**: FIXED -- kwargs are forwarded, T<0 raises, and the temperature branch is exhaustive.
- **Affects**: All of them, because the defect sits upstream of backend dispatch: `thermal.Thermal_Spin_Chain(sites, T=..., **kwargs)` and the `get_gs()` it drives. The constructor always builds `Spin_Chain(sitesT)` with the library default `itensor_version=3` (cppext.DEFAULT_ITENSOR_VERSION), so a caller asking for "python", 2, or "julia_live" silently gets v3 (or, where no C++ extension is compiled, mode.py's ED fallback). The T-branch defects (b)/(c) are backend-independent and hit the DMRG and ED paths alike.

**Expected**

<psi(t)|O|psi(t)> for the state actually handed in, i.e. proportional to
||psi||^2, which is what ED, itensor_version=2, itensor_version=3 with
tevol_method="TEBD"/"MPO", and itensor_version="python" all return (they track
the ED trajectory to 1e-4 or better and scale correctly with `scale`). The v3
TDVP/TDVP_GSE result is identical for scale=1.0, 0.5 and 3.0, i.e. the state is
silently renormalized at the first step; the answer is off by exactly
1/||psi||^2 (a factor 4 for A=Sx[0]). This is an oversight, not a convention:
the sibling methods quench_tdvp (chain_session.h ~1355), quench_tdvp_gse
(~1424) and quench_tebd (~1499) all follow their step with `psi1.normalize();
psi1 *= norm0;`, and CLAUDE.md documents a "manual norm-restoration workaround"
needed for the C++ METTS port for the same reason -- only
evolve_and_measure_tdvp and evolve_and_measure_tdvp_gse omit it. ROADMAP.md
marks TDVP and TDVP_GSE real-time evolution as implemented on v3, and both
itensor_version=3 and tevol_method="TDVP" are the defaults, so the wrong number
is what a default-configured user gets from evolution_ABA / evolve_and_measure
with any operator-perturbed start.

**Observed**

```
s20.py:
scale=1.0  ED  -0.23618   3/TDVP  -0.23618   3/TEBD  -0.23621   3/MPO  -0.23733   python/TDVP  -0.23618   2/MPO  -0.23733
scale=0.5  ED  -0.05905   3/TDVP  -0.23618   3/TEBD  -0.05905   3/MPO  -0.05933   python/TDVP  -0.05905   2/MPO  -0.05933
scale=3.0  ED  -2.12564   3/TDVP  -0.23618   3/TEBD  -2.12587   3/MPO  -2.13594   python/TDVP  -2.12564   2/MPO  -2.13594

f_tdvp_norm.py (public evolution_ABA path, A=Sx[0] so ||A|gs>||^2 = 0.25):
3        TDVP      [-0.059045 -0.235105 -0.231891 ...]  maxdev=0.176
3        TDVP_GSE  [-0.059045 -0.235105 -0.231891 ...]  maxdev=0.176
3        TEBD      [-0.059045 -0.058777 -0.057974 ...]  maxdev=1.57e-05
3        MPO       [-0.059045 -0.058775 -0.057963 ...]  maxdev=0.000252
2        MPO       [-0.059045 -0.058775 -0.057963 ...]  maxdev=0.000252
python   TDVP      [-0.059045 -0.05881  -0.058021 ...]  maxdev=0.000243
ED                 [-0.059045 -0.058775 -0.057963 ...]
```

**Reviewer** (independent, brief was to refute)

All three facets reproduce verbatim under the mandated PYTHONPATH (my own
script, /tmp/.../scratchpad/skeptic_thermal.py, confirmed importing
src/dmrgpy/__init__.py, v3 and v2 extensions both available):  (a)
`Thermal_Spin_Chain(spins, T=0.5, itensor_version="python")` ->
`tc.MBChain.itensor_version` is `3`. thermal.py:14 builds `Spin_Chain(sitesT)`
with no kwargs while the signature (thermal.py:8) advertises `**kwargs`. I
checked the "intended" defenses and they fail: fermionchain.py:225-236 carries
an explicit docstring saying **kwargs forwarding was added precisely so
`itensor_version` can be set at construction, and ROADMAP.md:106 lists
thermal.py as backend-agnostic and supported on v3/pyitensor/Julia. Sharpening
against the claimant: the parent `Many_Body_Chain.__init__` also swallows
unknown kwargs (`initialize(**kwargs)` ignores them), so `maxm=77`/`nsweeps=3`
are discarded library-wide, not by thermal.py -- "discards every constructor
kwarg" overstates it. `itensor_version` is the only kwarg genuinely lost. Blast
radius is real, not cosmetic: on this source install the swap is one correct
backend for another (skeptic_thermal2.py: E = -1.07290359 on the silently-used
v3 vs -1.07290564 after an explicit post-hoc `tc.MBChain.setup_python()`), but
on the pure-Python PyPI wheel there is no compiled extension, so the silently-
kept default of 3 routes through mode.py's ED fallback -- a user who asked for
pyitensor DMRG gets exact diagonalization of a *doubled* (physical+ancilla)
chain, which is infeasible at modest n. The documented workaround exists
(`setup_python()` after construction), so ROADMAP's coverage claim is not
falsified; the request is simply ignored without a warning.  (b) Narrower than
claimed. docs/user_guide.md:1759 states "T=0 recovers ordinary ground-state
DMRG", so the claimant's T=0.0 row is documented, correct behavior and is not
evidence of anything; "T<=0 should raise" contradicts the docs for the T=0
half. What does survive: `T=-1.0` and `T=-100.0` both return -1.6160254038,
identical to the T=0 result, because thermal.py:45's `elif self.T<1e-5` catches
negative temperatures too. An unphysical input silently aliases onto the
documented zero-temperature answer with no diagnostic, against the codebase's
own precedent at vevtk/mettsvev.py:100 (`if T <= 0: raise ValueError(...)`).
(c) The branches at thermal.py:41/45 are not exhaustive: `T == 1e-5` satisfies
neither `>1e-5` nor `<1e-5`, so `wf0` is never bound and line 48 raises
`UnboundLocalError: cannot access local variable 'wf0'`. Not a contrived input
-- `np.linspace(1e-5, 0.5, 4)[0] == 1e-5` exactly (verified,
skeptic_thermal3.py), and a temperature sweep starting at a small floor is the
natural way to drive this class; that grid point crashes.  Repro-flaw checks
done and cleared: no DMRG convergence question (facets (a)/(c) are structural,
and (b) was run with mode="ED"); no degeneracy/basis dependence (the compared
quantity is a scalar energy, and the T<0 values are bit-identical to the T=0
one, which is the point); chain is 4 sites so the known v3 <3-site abort is not
in play; nothing here depends on julia_live; and no sum rule or normalization
convention is involved. None of CLAUDE.md's deliberately reproduced legacy bugs
(evoloperator z^3/6 on H2, the "moise" key, the unreachable "tevol_fit_td"
branch) covers this file, and ROADMAP.md marks no restriction here.

**Where a fix goes**: In src/dmrgpy/thermal.py: forward the kwargs at line 14 (`Spin_Chain(sitesT, **kwargs)`), and make the temperature dispatch in `get_gs()` exhaustive by turning line 45's `elif self.T<1e-5` into a plain `else`, with an up-front `if T < 0: raise ValueError(...)` in `__init__` (and/or at the top of `get_gs`) mirroring vevtk/mettsvev.py:100.

**Repro**

```
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 s20.py

--- s20.py ---
import numpy as np, warnings
from dmrgpy import spinchain, timedependent
n=5
def mk(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    if v=="python": sc.setup_python()
    elif v!="ED": sc.setup_cpp(v)
    sc.maxm=30; sc.nsweeps=12
    H=0
    for i in range(n-1):
        H = H + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+0.6*sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): H = H + 0.4*sc.Sz[i]
    sc.set_hamiltonian(H)
    return sc
nt=10; dt=0.1
for scale in [1.0,0.5,3.0]:
    sc=mk("ED"); wf=sc.get_gs(mode="ED")
    ts,ref = timedependent.evolve_and_measure(sc,mode="ED",nt=nt,dt=dt,operator=sc.Sz[2],wf=scale*wf)
    ref=np.array(ref).real
    row=["scale=%.1f  ED %9.5f"%(scale,ref[-1])]
    for ver,tm in [(3,"TDVP"),(3,"TEBD"),(3,"MPO"),("python","TDVP"),(2,"MPO")]:
        sc=mk(ver); sc.tevol_method=tm; w=sc.get_gs()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t2,y=timedependent.evolve_and_measure(sc,mode="DMRG",nt=nt,dt=dt,operator=sc.Sz[2],wf=scale*w)
        row.append("%s/%s %9.5f"%(ver,tm,np.array(y).real[-1]))
    print("   ".join(row))

A second repro through the fully public entry point evolution_ABA (same command line, file f_tdvp_norm.py):

import numpy as np, warnings
from dmrgpy import spinchain, timedependent
n=5
def mk(v):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    if v=="python": sc.setup_python()
    elif v!="ED": sc.setup_cpp(v)
    sc.maxm=30; sc.nsweeps=12
    H=0
    for i in range(n-1):
        H = H + sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+0.6*sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): H = H + 0.4*sc.Sz[i]
    sc.set_hamiltonian(H)
    return sc
nt=10; dt=0.1
sc=mk("ED"); sc.get_gs(mode="ED")
ts,ref = timedependent.evolution_ABA(sc,A=sc.Sx[0],B=sc.Sz[2],mode="ED",nt=nt,dt=dt)
ref=np.array(ref).real
print("ED       ", np.round(ref,6))
for ver,tm in [(3,"TDVP"),(3,"TDVP_GSE"),(3,"TEBD"),(3,"MPO"),(2,"MPO"),("python","TDVP")]:
    sc=mk(ver); sc.tevol_method=tm; sc.get_gs()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t2,y = timedependent.evolution_ABA(sc,A=sc.Sx[0],B=sc.Sz[2],mode="DMRG",nt=nt,dt=dt)
    y=np.array(y).real
    print("%-8s %-9s"%(ver,tm), np.round(y,6), " maxdev=%.3g"%np.max(np.abs(y-ref)))
```

### 14. Many_Body_Chain.get_distribution / get_distribution_moments are the only public entry points that never call get_mode(), so all three routes into ED (enforced self.mode=\"ED\", the v3-with-ns<3 fallback, and the no-compiled-extension case) crash with AttributeError instead of using the ED implementation that exists and works

- **Severity**: medium  **Class**: missing-combo
- **Code**: `src/dmrgpy/manybodychain.py:795 (and :805)`
- **Status**: FIXED -- both entry points dispatch through get_mode(); the ED branch raises NotImplementedError instead of a bare raise.
- **Affects**: Backend-independent: the crash is in shared Python (`manybodychain.get_distribution`/`get_distribution_moments` -> `distribution.py` -> `kpmdmrg.general_kpm_moments` -> `general_kpm_moments_cpp_ext`, which is `self._session`-only with no version branching), so it hits `itensor_version` 2, 3 and "python" identically. Public API calls affected: `Many_Body_Chain.get_distribution(...)` and `Many_Body_Chain.get_distribution_moments(...)` called with default `mode="DMRG"`. Three trigger conditions, verified with the compiled v2+v3 extensions present in this checkout: (a) any chain with `self.mode="ED"` set; (b) any chain with `itensor_version=3` (the default) and `ns<3`; (c) any install with no compiled extension at all -- i.e. every PyPI-wheel user, on an ordinary chain of any length with all defaults. Not tested: `itensor_version="julia_live"` (per instructions, no juliacall run); it has no `_session`, so it is presumably a fourth instance, unverified.

**Expected**

Every other public entry point (vev, gs_energy, get_gs, get_excited,
get_dynamical_correlator, ...) starts with `mode = self.get_mode(mode=mode)`,
which honors an enforced self.mode="ED" and applies mode.py's automatic
ITensor-v3-with-ns<3 -> ED fallback. get_distribution/get_distribution_moments
omit that line, so (i) a chain with self.mode="ED" runs the C++ session anyway
(my spy confirms general_kpm_moments_cpp_ext is entered) and dies on the ED
State it was handed, and (ii) a plain 2-site chain on the DEFAULT backend
crashes even though the ED implementation it should have routed to exists and
returns a correct spectrum when mode="ED" is passed explicitly. The task says
the ns<3 mitigation is reportable if the mitigation itself has a gap; this is
such a gap, and the entanglement leg (get_bond_entropy on a 2-site v3 chain,
all defaults) is a second instance of the same root cause: the fallback returns
an ED State while entropy.py is self._session-only. That also contradicts
CLAUDE.md's claim that the versions= v3 exclusions in
tests/test_entanglement.py are "redundant defense-in-depth rather than required
to avoid a crash" -- on this path the exclusion is load-bearing. Secondary:
get_distribution_moments(mode="ED") hits a bare `raise` and surfaces as
"RuntimeError: No active exception to reraise" instead of a NotImplementedError
naming the problem.

**Observed**

```
r1_dist.py:
A get_mode() -> ED
A raised AttributeError 'State' object has no attribute 'cpp_handle'
ITensor v3's two-site DMRG can't handle a chain this short (n=2 < 3 sites), using default ED routines
B get_mode() -> ED  gs_energy -> -0.75
B calling get_distribution ...
Traceback (most recent call last):
  File "r1_dist.py", line 34, in <module>
    x, y = sc2.get_distribution(X=sc2.Sz[0], delta=0.4)
  File "src/dmrgpy/manybodychain.py", line 799, in get_distribution
    return distribution.get_distribution(self,**kwargs)
  File "src/dmrgpy/kpmdmrg.py", line 208, in general_kpm_moments_cpp_ext
    mus = self._session.general_kpm(X.to_terms(),wfa.cpp_handle,wfb.cpp_handle,
AttributeError: 'State' object has no attribute 'cpp_handle'
EXIT=1

r1b.py:
mode='ED' works: [-9.5        -9.48098098 -9.46196196] [1.36449025e-07+0.j 1.52108147e-07+0.j 1.60734083e-07+0.j]

r5_probe.py entropy:
  File "src/dmrgpy/entropy.py", line 35, in compute_entropy_single
    return np.abs(self._session.bond_entropy(psi.cpp_handle,b))
AttributeError: 'State' object has no attribute 'cpp_handle'
EXIT=1

r9_extra.py:
get_distribution_moments(mode='ED') -> RuntimeError: No active exception to reraise
```

**Reviewer** (independent, brief was to refute)

REAL, and the widest leg is the one the lead understates. `mode.py::get_mode()`
is the single place that resolves DMRG-vs-ED, and it handles three fallbacks:
extension-not-compiled, `itensor_version==3 and ns<3`, and an enforced
`self.mode`. `vev`, `gs_energy`, `get_gs`, `get_excited`, `get_excited_states`,
`get_dynamical_correlator` all open with `mode = self.get_mode(mode=mode)`.
`get_distribution` (manybodychain.py:795) and `get_distribution_moments` (:805)
branch on the raw `mode` argument only, so with the default `mode="DMRG"` they
enter the C++/session path no matter what `get_mode()` would have said.
Verified with my own script (`skep06.py`, printed `dmrgpy.__file__` to rule out
a stale import; v2 and v3 both compiled here): - n=5 chain, `sc.mode="ED"`:
`get_mode() -> ED`, then `get_distribution` raises `AttributeError: 'State'
object has no attribute 'cpp_handle'` (the ED `State` returned by `get_gs()`,
which *does* honor the mode, handed to session-only code). - n=2 chain, default
v3: `get_mode() -> ED`, `gs_energy() -> -0.75` fine, `get_distribution` raises
the same AttributeError. Same chain with an explicit `mode="ED"` returns a
correct spectrum, so the ED implementation (`edtk/edchain.py:105`) exists, is
wired up, and is simply unreachable. - `skep06c.py`, patching
`cppext.get_backend -> None` to reproduce exactly a pip-wheel install (the PyPI
wheel ships no C++ at all, per CLAUDE.md's packaging section): a plain 5-site
Heisenberg chain, all defaults, `get_mode() -> ED`, `_session is None`,
`gs_energy()` works, and `get_distribution(X=sc.Sz[0],delta=0.4)` raises
`AttributeError: 'NoneType' object has no attribute 'set_sweep_params'`. This
is the real blast radius: for wheel users these two methods are unusable unless
they happen to pass `mode="ED"` by hand, while every neighbouring method
silently does the right thing.  Not intended. ROADMAP.md contains no row for
distribution/`general_kpm` and marks no such restriction. The codebase's own
stated convention is the opposite: `groundstate.py:196-198` says the other
methods "all honor self.mode=\"ED\" for cross-validation" and that where an ED
implementation is genuinely absent, `self.mode=\"ED\"` is "rejected explicitly
below rather than silently ignored" (with `tests/test_dmrg_generalized.py:277`
locking that in). `get_distribution` does neither -- it neither honors nor
rejects, it ignores and then dies on a type mismatch. It is also not one of
CLAUDE.md's three deliberately reproduced legacy bugs, and it is a gap in the
ns<3 mitigation CLAUDE.md documents as a fix rather than a fresh restriction.
Repro flaws ruled out: the failure is an AttributeError, not a numeric
discrepancy, so convergence (maxm/nsweeps), ground-state degeneracy and sum-
rule/normalization conventions are all irrelevant; I ran at maxm=20/nsweeps=8
anyway. The n<3 leg is not "v3 aborting on a short chain" -- v3 never runs, the
mitigation correctly routes to ED and the crash is downstream of it. Stale-
import trap excluded by the explicit PYTHONPATH plus printing
`dmrgpy.__file__`.  Two secondary legs from the lead, which I reproduced but
which should NOT be folded into this claim: - The entanglement leg
(`get_bond_entropy` on a 2-site v3 chain -> `'State' object has no attribute
'cpp_handle'`) reproduces, but it is a different root cause and is documented:
CLAUDE.md's benchmarking section states "entropy.py is self._session-only and
has no ED path of its own". Adding `get_mode()` there would fix nothing since
there is no ED branch to route to. The lead's "this contradicts CLAUDE.md's
redundant-defense-in-depth claim" argument does not hold either:
`tests/test_entanglement.py` has no `versions=` exclusion at all, and its
docstring shows the Bell-pair test deliberately uses the chain-level
`get_site_entropy`, which handles an ED `State` fine. -
`get_distribution_moments(mode="ED")` -> `RuntimeError: No active exception to
reraise` reproduces. The *absence* of an ED moments path is consistent with the
documented design (`get_dynamical_correlator_moments`' docstring: "there is no
ED counterpart, since the ED path builds spectra by explicit summation rather
than from moments"). Only the surfacing is defective -- a bare `raise` where
the convention calls for a NotImplementedError naming the gap.

**Where a fix goes**: Add `mode = self.get_mode(mode=mode)` as the first line of both `get_distribution` and `get_distribution_moments` in `src/dmrgpy/manybodychain.py` (~795 and ~805), matching every neighbouring entry point; in `get_distribution_moments`, replace the bare `raise` in the ED branch with an explicit `NotImplementedError` naming the missing ED moments path, following the reject-explicitly convention `groundstate.py` describes.

**Repro**

```
Script A /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/r1_dist.py:

"""get_distribution/get_distribution_moments never call get_mode()."""
from dmrgpy import spinchain

def build(n, **kw):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, **kw)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 10; sc.nsweeps = 4
    return sc

sc = build(5)
sc.mode = "ED"
print("A get_mode() ->", sc.get_mode())
import dmrgpy.kpmdmrg as kpmdmrg
orig = kpmdmrg.general_kpm_moments_cpp_ext
hit = []
def spy(self, X, wfa, wfb, num_p, acc):
    hit.append(True)
    return orig(self, X, wfa, wfb, num_p, acc)
kpmdmrg.general_kpm_moments_cpp_ext = spy
try:
    x, y = sc.get_distribution(X=sc.Sz[0], delta=0.4)
    print("A get_distribution returned, DMRG C++ session used:", bool(hit))
except Exception as e:
    print("A raised", type(e).__name__, e)

sc2 = build(2)   # itensor_version=3 default
print("B get_mode() ->", sc2.get_mode(), " gs_energy ->", sc2.gs_energy())
print("B calling get_distribution ...", flush=True)
x, y = sc2.get_distribution(X=sc2.Sz[0], delta=0.4)
print("B ok", x[:2], y[:2])

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r1_dist.py; echo EXIT=$?

Script B (the ED implementation exists and works when asked for explicitly), r1b.py:

from dmrgpy import spinchain
def build(n, **kw):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, **kw)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 10; sc.nsweeps = 4
    return sc
sc2 = build(2)
x,y = sc2.get_distribution(mode="ED", X=sc2.Sz[0], delta=0.4)
print("mode='ED' works:", x[:3], y[:3])

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r1b.py

Script C (same ED-fallback gap on the entanglement path), r5_probe.py:

import sys, numpy as np
from dmrgpy import spinchain
which = sys.argv[1]
n = 2
sc = spinchain.Spin_Chain([2]*n)
h = sc.Sx[0]*sc.Sx[1]+sc.Sy[0]*sc.Sy[1]+sc.Sz[0]*sc.Sz[1]
sc.set_hamiltonian(h); sc.maxm=10; sc.nsweeps=4
print("get_mode:", sc.get_mode())
if which=="rdm":     print(sc.get_rdm(pairs=[(0,1)]))
elif which=="corrmat": print(sc.get_correlation_matrix(name="ZZ"))
elif which=="entropy": print(sc.get_bond_entropy(sc.get_gs(),0,1))
elif which=="excited": print(sc.get_excited(n=2))

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r5_probe.py entropy

Script D (secondary symptom), r9_extra.py second half:

sc2 = spinchain.Spin_Chain([2]*4); h2=0
for i in range(3): h2 = h2 + sc2.Sz[i]*sc2.Sz[i+1]
sc2.set_hamiltonian(h2)
try:
    sc2.get_distribution_moments(mode="ED", X=sc2.Sz[0])
except Exception as e:
    print("get_distribution_moments(mode='ED') ->", type(e).__name__+":", e)

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r9_extra.py
```

### 15. get_excited_states(n>=2, scale=...) on a non-Hermitian Hamiltonian raises TypeError("mpsiram() got an unexpected keyword argument 'scale'") four frames deep; the same forward makes dcex's submode="EX" correlator unusable on a non-Hermitian H

- **Severity**: medium  **Class**: crash
- **Code**: `src/dmrgpy/excited.py:85 (return excited_states_non_hermitian(self,n=n,**kwargs))`
- **Status**: FIXED -- excited_states_non_hermitian accepts and ignores scale, alongside the existing nkry_min/nkry_max compatibility slots.
- **Affects**: Backend-agnostic: the failing frames are pure Python (excited.py -> algebra/arpacktk.py) and are reached before any session/backend call. Confirmed by direct run on itensor_version=3 (compiled mpscpp3, cppext.available(3)=True) and itensor_version="python"; v2 and "julia_live" are affected identically by inspection (same code path). mode="ED" is NOT affected (edtk/edchain.get_excited_states -> algebra.lowest_states swallows scale). Public API calls hit: Many_Body_Chain.get_excited_states(n>=2, scale=...), and its wrappers get_excited(n>=2, scale=...) and get_gap(scale=...) (get_gap forwards **kwargs with n=2), whenever is_hermitian(self.hamiltonian) is False and mode resolves to "DMRG". Also hit with no user-supplied scale at all through dcex.dynamical_correlator (submode="EX"), which unconditionally forwards scale=10.0 (dcex.py:21,56); on itensor_version="julia_live" that is a fully public route (dynamics.py:107-113 dispatches submode="EX" to dcex and deliberately does not block non-Hermitian H there), while on (2,3,"python") dynamics.py:11-16 currently pre-empts every submode for non-Hermitian H, so the EX crash is only reachable by calling dcex.dynamical_correlator directly today.

**Expected**

`scale` is the documented Lagrange-multiplier strength of the overlap-penalty
excited-state method and a normal public kwarg of get_excited_states
(excited.py:34, get_excited_states_dmrg(self,n=2,noise=0.0,scale=10.0)). On a
non-Hermitian chain excited.py:85 forwards **kwargs straight into
excited_states_non_hermitian, which has no `scale` parameter and no catch-all,
so it flows down into arpacktk.mpsiram and dies several frames deep with a
TypeError naming an internal function. It should either be accepted and ignored
-- as nkry_min/nkry_max already are: excited_states_non_hermitian's own comment
at excited.py:137-141 says they exist "only for backward compatibility"
precisely because they "would otherwise hit a TypeError several calls deep",
i.e. this exact bug class was already fixed once for other kwargs and `scale`
was missed -- or rejected with a message naming `scale`. This also compounds
finding 3: dcex.dynamical_correlator always forwards scale=10.0 explicitly
(dcex.py:56), so if the silent submode substitution at dynamics.py:15 were
removed, submode="EX" on a non-Hermitian Hamiltonian -- which ROADMAP.md
section 4 advertises as "Non-Hermitian-capable" -- would crash immediately with
this TypeError.

**Observed**

```
nh3.py:
hermitian? False
no scale: [-1.56250406+1.13437038e-12j -0.9906074 +3.66633390e-11j
 -0.948398  +1.73293488e-01j]
scale=10 RAISED TypeError mpsiram() got an unexpected keyword argument 'scale'

nh2.py (traceback tail):
  File "/u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/algebra/arpacktk.py", line 709, in excited_states
    es, wfs = mpsiram(self, H, which=which, nev=1, wfskip=wfskip, **kwargs)
TypeError: mpsiram() got an unexpected keyword argument 'scale'
```

**Reviewer** (independent, brief was to refute)

Real, and not explicable as intended. Verified with my own script under
PYTHONPATH=<repo>/src (dmrgpy resolved to
/u/40/.../dmrgpy/src/dmrgpy/__init__.py, v3 extension present): a 4-site S=1/2
Heisenberg chain plus 0.4j*Sz[0] (is_hermitian -> False), maxm=20, nsweeps=6.
get_excited_states(n=3) succeeds on both backends and returns identical complex
energies [-1.562504, -0.990607, -0.948398-0.173293j]; get_excited_states(n=3,
scale=10.0) raises TypeError at algebra/arpacktk.py:709 on both. Chain:
excited.py:85 forwards **kwargs into excited_states_non_hermitian, whose
signature (excited.py:135-141) names only recursive/maxit/ncv/nkry_min/nkry_max
and has a **kwargs catch-all that it passes straight to
arpacktk.excited_states, which passes it to mpsiram (arpacktk.py:309), and
mpsiram has no catch-all. I ruled out the usual false positives: this is a
deterministic TypeError, so no convergence/degeneracy/basis question arises;
the chain has 4 sites so mpscpp3's known <3-site abort is irrelevant; both a
compiled and a pure-Python backend give the identical failure, so it is not
extension-specific; no sum rule or normalization is involved. I then tried hard
to call it intended and could not. (1) `scale` is a documented public knob of
this exact method: docs/user_guide.md:494 tells users a warned-about excited
state can be "recomputed with a larger `scale`/more sweeps", and
docs/user_guide.md:599 documents that get_excited_states routes non-Hermitian H
with n>=2 through excited_states_non_hermitian, so the documentation invites
precisely the call that crashes. (2) The behaviour is inconsistent across every
neighbouring combination I tested, and the failing one is the odd man out:
Hermitian+DMRG accepts and uses scale (n=3, scale=10 and scale=5 both return
[-1.616025,-0.957107,-0.957107]); Hermitian+ED accepts and ignores it; non-
Hermitian with n=1 accepts and ignores it (excited.py:78 routes to
gs_energy(**kwargs)); non-Hermitian with n>=2 is the only one that raises. (3)
The codebase already fixed this exact bug class once and simply missed `scale`:
excited.py:141-147's own comment says nkry_min/nkry_max are accepted "only for
backward compatibility" because otherwise "they'd flow through **kwargs into
mpsiram, which has no catch-all" -- and I confirmed get_excited_states(n=2,
nkry_max=8) works on the same non-Hermitian chain while scale=10.0 does not.
(4) It is not a ROADMAP gap: ROADMAP.md:49 lists non-Hermitian excited states
(Arnoldi/IRAM, arpacktk) as supported on all three backends and ROADMAP.md:99
lists submode="EX" as "Non-Hermitian-capable; entirely backend-agnostic". (5)
It is not one of CLAUDE.md's deliberately reproduced legacy bugs (evoloperator
z^3/6, "moise", "tevol_fit_td"). (6) It is not merely user error, because an
internal caller triggers it with no user kwarg: I ran
dcex.dynamical_correlator(sc, name=(Sz1,Sz1), es=..., delta=0.2, nex=6) on the
same non-Hermitian v3 chain and got the same TypeError, while
excited.get_excited_states(sc, n=6, purify=False) (identical call minus the
scale forward) succeeds -- so the only thing breaking the EX correlator on a
non-Hermitian Hamiltonian is this forwarded kwarg. Blast radius is therefore: a
documented public kwarg is a hard crash on one branch of a public method with
an error message naming an internal ARPACK-port function, plus submode="EX" is
dead on non-Hermitian Hamiltonians for julia_live (advertised as working by
both ROADMAP.md:99 and dynamics.py's own comment at lines 40-48) and would be
dead on (2,3,"python") too the moment the separate non-Hermitian pre-empt at
dynamics.py:11-16 were narrowed. I did not run julia_live per instructions;
that leg is established statically plus by the direct dcex call, which
exercises the entire failing code path.

**Where a fix goes**: Add `scale=None` as an explicitly named, accepted-and-ignored parameter of `excited_states_non_hermitian` in src/dmrgpy/excited.py (right beside the existing `nkry_min`/`nkry_max` backward-compatibility slots, with a one-line comment saying the overlap-penalty weight has no meaning for the IRAM solve), or equivalently strip `scale` from `kwargs` at the excited.py:85 dispatch; either way it must not reach `arpacktk.excited_states`/`mpsiram`.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/nh3.py:

import numpy as np
from dmrgpy import spinchain
n=4
sc=spinchain.Spin_Chain(["S=1/2"]*n)
h=0
for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
h=h+0.4j*sc.Sz[0]
sc.setup_cpp(version=3); sc.maxm=20; sc.nsweeps=8
sc.set_hamiltonian(h)
print("hermitian?",sc.is_hermitian(sc.hamiltonian))
try:
    print("no scale:",sc.get_excited_states(n=3)[0])
except Exception as ex: print("no scale RAISED",type(ex).__name__,str(ex)[:120])
try:
    print("scale=10:",sc.get_excited_states(n=3,scale=10.0)[0])
except Exception as ex: print("scale=10 RAISED",type(ex).__name__,str(ex)[:120])

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src timeout 300 python3 nh3.py

The same TypeError reached through the EX path, script nh2.py, which calls dcex.dynamical_correlator directly to bypass finding 3's silent substitution:

import numpy as np
from dmrgpy import spinchain, dcex
n=4
sc=spinchain.Spin_Chain(["S=1/2"]*n)
h=0
for i in range(n-1): h=h+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
h=h+0.4j*sc.Sz[0]
sc.setup_cpp(version=3); sc.maxm=20; sc.nsweeps=8
sc.set_hamiltonian(h)
es=np.linspace(0.,4.,30)
e,cpub=sc.get_dynamical_correlator(name=(sc.Sz[1],sc.Sz[1]),submode="EX",es=es,delta=0.2)
e2,cex=dcex.dynamical_correlator(sc,name=(sc.Sz[1],sc.Sz[1]),es=es,delta=0.2)
cpub=np.array(cpub); cex=np.array(cex)
print("maxdiff",np.max(np.abs(cpub-cex)))

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src timeout 400 python3 nh2.py
```

### 16. get_rdm(i=ns-1) raises IndexError on itensor_version="python" because pyitensor's reduced_dm indexes psi.A(site+1) past the end of its tensor list; itensor_version=2/3 return the correct last-site RDM

- **Severity**: medium  **Class**: crash
- **Code**: `/u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/pyitensor/chain.py:642`
- **Status**: FIXED -- pyitensor's reduced_dm handles the last site, where there is no right link to prime.
- **Affects**: Broken: itensor_version="python" (pyitensor), via the single public call Many_Body_Chain.get_rdm(i=ns-1) -> densitymatrix.reduced_dm -> Chain.reduced_dm. By inspection the same one-line pattern exists in mpsjulialive/densitymatrix.jl (ir = commonind(psi[site],psi[site+1]); ITensors.jl psi[N+1] is a BoundsError), so julia_live is very likely affected too, but per instructions it was NOT executed and is unverified. Unaffected: itensor_version=2 and 3 (correct answer), and mode="ED" chains have no get_rdm path at all. Only get_rdm is hit; the public entropy API is not (see explanation).

**Expected**

The 2x2 reduced density matrix of the last site, matching the ED reference to
~1e-13 the way sites 0..ns-2 already do on this backend and the way all ns
sites do on itensor_version=2 and 3. pyitensor/chain.py:642 unconditionally
reaches for psi.A(site+1) to find the right bond index, which does not exist
for the rightmost site. ROADMAP.md's section-2 table marks "Bond entanglement
entropy / reduced density matrix" as implemented for pyitensor with the note
that entanglement.py has no itensor_version branching "at all -- every
backend's session exposes the same primitive", so this is a row marked
implemented that breaks for one index.

**Observed**

```
site 0  exact=[0.422138 0.577862]  2 dev=2.1e-11  3 dev=4.6e-11  python dev=3.0e-13
site 1  exact=[0.533264 0.466736]  2 dev=9.8e-12  3 dev=6.9e-12  python dev=6.9e-13
site 2  exact=[0.414677 0.585323]  2 dev=2.2e-11  3 dev=2.1e-11  python dev=4.9e-13
site 3  exact=[0.533264 0.466736]  2 dev=1.4e-11  3 dev=3.2e-12  python dev=5.3e-13
site 4  exact=[0.422138 0.577862]  2 dev=3.2e-11  3 dev=1.9e-11  python RAISED IndexError(list index out of range)

Full traceback:
  File ".../dmrgpy/pyitensor/chain.py", line 642, in reduced_dm
    ir = commonIndex(psi.A(site), psi.A(site + 1))
  File ".../dmrgpy/pyitensor/mpscontainer.py", line 45, in A
    return self._tensors[i - 1]
IndexError: list index out of range
```

**Reviewer** (independent, brief was to refute)

Confirmed with my own script (5 sites, maxm=20, nsweeps=16, XXZ + 0.3*Sx +
0.11*Sz fields to break every symmetry, PYTHONPATH pinned to the checkout).
Against an exact partial trace of the ED ground-state vector: sites 0..3 agree
on all three backends (v2 ~1.4e-11, v3 ~1.0e-11, python ~5e-13); site 4 gives
v2 1.41e-11, v3 7.21e-12, and python raises IndexError at
pyitensor/chain.py:642 -> mpscontainer.py:45. So this is a genuine cross-
backend divergence on a capability ROADMAP.md marks as implemented for
pyitensor with the note that dispatch is backend-agnostic.  Why the C++ works
and the port does not: ITensor's MPS allocates A_ with N+2 slots ("idmrg may
use A_[0] and A[N+1]", mps.cc:67/81/96/110), so psi.A(N+1) is a valid default-
constructed ITensor; commonIndex against it returns an invalid Index and
prime(T,s,ir) then primes only the site index -- exactly the right thing at the
right edge. pyitensor's mpscontainer deliberately drops the trivial boundary
Link and stores a plain Python list, so A(N+1) is an out-of-range list access.
I tried hard to read this as intended and it does not hold up. chain.py:634-638
and docs/documentation.md:3247-3252 both assert the missing bounds check is "a
pre-existing, deliberately-preserved limitation shared identically by
pyitensor/chain.py and mpscpp3/chain_session.h ... not a new risk, so left as-
is". That premise is measurably false for reduced_dm: the C++ backends do not
share it, they return the correct matrix. It is true for the neighbouring
bond_entropy (b=ns aborts v3 with ITensor's "Default constructed ITensor in
product" from itensor.cc:940, because that one multiplies by the empty tensor
rather than only asking it for a common index) -- but that path is unreachable
from the public API (mps.get_entropy loops b in range(0,ns-1) and
get_bond_entropy(i,j=i+1) caps compute_entropy_single at b=ns-1; I measured b=1
and b=ns-1 agreeing to ~1e-11 between v3 and python). This is also not one of
CLAUDE.md's three deliberately reproduced legacy bugs, and the separate
documented quirk in this function (dividing by the norm squared) is a genuine
no-op, unrelated.  Repro flaws ruled out explicitly: not unconverged DMRG --
the failure is a pure Python IndexError with no numerics involved, and the same
run's sites 0..ns-2 match ED to 1e-13 on the same backend/parameters; not a
degenerate ground state -- the transverse+longitudinal fields leave a
nondegenerate ED ground state that all three backends reproduce site by site;
not a too-short chain -- 5 sites, well clear of v3's known <3-site abort; not
an apples-to-oranges comparison -- the reference is the exact partial trace of
the normalized ED vector, the identical quantity the other four sites match.
i=ns-1 is a documented, in-range argument (docs/user_guide.md section 14
documents get_rdm(i=...) with no index restriction).  Blast radius is narrow:
one index on one (probably two) backends, a loud crash rather than a wrong
number, trivially worked around. Note it is already listed unfixed as item 16
in the committed docs/audit_2026_08_hole_hunt.md index, so this is a re-
discovery, not a new regression.

**Where a fix goes**: In pyitensor/chain.py::Chain.reduced_dm, make the right bond optional at the edge: ir = commonIndex(psi.A(site), psi.A(site+1)) if site < psi.length() else None, and prime only the site index when ir is None -- mpscontainer._link_at already encodes exactly this boundary rule and could be reused. The same one-line guard belongs in mpsjulialive/densitymatrix.jl's reduced_dm, and the stale "the original never defended against this either" comment plus the matching paragraph in docs/documentation.md should be corrected.

**Repro**

```
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 f_rdm_last.py

--- f_rdm_last.py ---
import numpy as np
from dmrgpy import spinchain
n=5
def mk(v):
    sc=spinchain.Spin_Chain(["S=1/2"]*n)
    if v=="python": sc.setup_python()
    elif v!="ED": sc.setup_cpp(v)
    sc.maxm=30; sc.nsweeps=14
    H=0
    for i in range(n-1): H=H+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+0.7*sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): H=H+0.3*sc.Sx[i]+0.11*sc.Sz[i]
    sc.set_hamiltonian(H)
    return sc
sc=mk("ED"); v=np.array(sc.get_gs(mode="ED").v); v=v/np.linalg.norm(v); T=v.reshape([2]*n)
P=np.array([[0,1],[1,0]])   # basis-order flip between ED vector and session rdm
for i in range(n):
    A=np.moveaxis(T,i,0).reshape(2,-1); rho=P@(A@A.conj().T)@P
    line="site %d  exact=%s"%(i,np.round(np.diag(rho).real,6))
    for ver in [2,3,"python"]:
        s=mk(ver)
        try:
            dm=np.array(s.get_rdm(i=i))
            line+="  %s dev=%.1e"%(ver,np.max(np.abs(dm-rho)))
        except Exception as e:
            line+="  %s RAISED %s(%s)"%(ver,type(e).__name__,str(e)[:40])
    print(line)
```

### 17. get_dynamical_correlator(submode="CVMimag") and submode="maxent" are permanently dead and terminate the caller's process with exit status 0 (SystemExit, not catchable by `except Exception`), because `dmrgpy.padetk` / `dmrgpy.maxenttk` have never existed in the tree

- **Severity**: medium  **Class**: should-raise
- **Code**: `/u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/analyticcontinuation.py:10 and /u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/distribution.py:32`
- **Status**: FIXED -- both sites raise NotImplementedError naming the missing optional module instead of terminating the process.
- **Affects**: Public API: `Many_Body_Chain.get_dynamical_correlator(mode="DMRG", submode=...)`.
- `submode="CVMimag"`: dead on `itensor_version` 2, 3 and `"python"` (all three share the single `if self.itensor_version in (2,3,"python")` branch in `dynamics.py`; I ran 3 and `"python"` directly, v2 by code-read of that shared dispatch). NOT affected on `"julia_live"` — `mpsjulialive/dynamics.py` raises a clean `NotImplementedError` for that submode.
- `submode="maxent"`: dead on 2, 3, `"python"` AND `"julia_live"` — `mpsjulialive/dynamics.py:114` explicitly routes into the very same `distribution.dynamical_correlator_positive_defined`. (julia_live established by code-read only; per instruction 4 I did not run juliacall.)
- `mode="ED"` is NOT affected: `edtk/dynamics.py:82` raises a proper `NotImplementedError` for unsupported submodes.
- Secondary entry point: `distribution.get_distribution_maxent` (module-level function), reachable only by direct import — it is not bound as a method on `Many_Body_Chain`, contrary to what `docs/user_guide.md` §14 implies.

**Expected**

ROADMAP.md's section-4 table marks `CVMimag` as implemented for both
itensor_version=3 and "python" (and `maxent` for v3, pyitensor and Julia). Both
are in fact unreachable: analyticcontinuation.imag2real does `try: from .padetk
import bruteforcepade` / `except: print("Not functional yet"); exit()`, and
distribution.get_distribution_maxent does the same with `from
.maxenttk.pymaxent import reconstruct`. Neither `src/dmrgpy/padetk/` nor
`src/dmrgpy/maxenttk/` exists in the tree (`ls src/dmrgpy | grep -i
'pade\|maxent'` returns nothing), so the except branch always fires. Two
problems: (a) two ROADMAP-implemented rows are entirely dead on every backend;
(b) a library must raise (ImportError/NotImplementedError), not call `exit()`
-- SystemExit escapes ordinary `except Exception` handlers and terminates the
caller's script, so a user sweeping submodes or running a batch job loses the
whole run with only a "Not functional yet" print.

**Observed**

```
[CVMimag/3] Not functional yet|EXC CVMimag 3 SystemExit None|SCRIPT REACHED END|
[CVMimag/python] Not functional yet|EXC CVMimag python SystemExit None|SCRIPT REACHED END|
[maxent/3] Not functional yet|EXC maxent 3 SystemExit None|SCRIPT REACHED END|
[maxent/python] Not functional yet|EXC maxent python SystemExit None|SCRIPT REACHED END|
[ROOTN/3] OK ROOTN 3 integral= 0.22385485862178298|SCRIPT REACHED END|
[ROOTN/python] OK ROOTN python integral= 0.22385485862187934|SCRIPT REACHED END|

Earlier, an unguarded sweep script (s31.py) that wrapped every submode in `try/except Exception` died mid-loop right after printing "Not functional yet", because SystemExit is not an Exception subclass.
```

**Reviewer** (independent, brief was to refute)

Real and reproduced verbatim, with the blast radius worse than the lead states:
the failure exits the process with status 0.  My unguarded-sweep script
(`.../scratchpad/sk/s2.py`, a realistic user loop over
`["CVM","maxent","TDZ","EX"]` on `itensor_version=3`, n=4, maxm=20, wrapped in
the ordinary `try/except Exception`) printed `OK CVM 0.2237744790470...`, then
`Not functional yet`, then died. `TDZ` and `EX` never ran, the final `SWEEP
COMPLETED ALL FOUR` never printed, and the shell saw **exit code 0** — a
batch/CI job silently loses the remainder of the run and reports success.
`exit()` with no argument raises `SystemExit(None)` → status 0, and
`SystemExit` inherits from `BaseException`, so it slips past every `except
Exception`.  Per-submode reproduction (n=4, maxm=20, nsweeps=12, `PYTHONPATH`
pinned to this checkout; `DMRGPY FROM:
/u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/__init__.py` printed on
every run, so no stale site-packages copy): `CVMimag/3`, `CVMimag/python`,
`maxent/3`, `maxent/python` all print `Not functional yet` and hit my `except
BaseException` arm as `SystemExit None`. Control: `CVM/3` and `CVM/python`
return normally, integral 0.2237744790470507 / 0.22377447904704728 — so this is
not unconverged DMRG, not a degenerate-ground-state artifact, not a too-short
chain, and not a comparison of two things that were never equal. Nothing
numerical is involved at all: the failure is an unconditional `except:
print(...); exit()` around a missing import, hit before any physics runs
(`analyticcontinuation.py:6-12`, `distribution.py:30-33`).  Not intended, on
every test the task asks for: - The missing packages never existed: `ls
src/dmrgpy | grep -iE 'pade|maxent'` is empty and `git log --all --
src/dmrgpy/padetk src/dmrgpy/maxenttk` returns nothing. These are not a
regression, they are permanently dead branches. - Not one of CLAUDE.md's
deliberately reproduced legacy bugs (evoloperator's z³/6 on H2, the `"moise"`
key, the unreachable `"tevol_fit_td"` branch). - Not a documented backend
restriction. `ROADMAP.md` §4 marks `CVM_explicit, CVMimag` as ✅ for v3 and
pyitensor, and `maxent` as ✅ for v3, pyitensor **and** Julia ("backend-agnostic
post-processing over `power_vev` moments"). `README.md:165` advertises both
submodes. `docs/user_guide.md` §6 documents `submode="maxent"` and §15
documents the Padé continuation behind `CVMimag`, neither with any caveat. Only
`docs/documentation.md:2798-2806` admits the maxent gap ("wired up for parity
but not actually functional on *any* backend right now ... a pre-existing,
cross-backend gap"); `CVMimag`'s deadness is documented nowhere, and even that
passage does not defend `exit()`. - Not API misuse: my call is the documented
signature, and the sibling `CVM`/`TDZ`/`EX` submodes work through the identical
call. - The codebase's own intended pattern for an unimplemented submode is
right next door: `edtk/dynamics.py:80-83` and
`mpsjulialive/dynamics.py:121-125` both `raise NotImplementedError`.
`src/dmrgpy/reconstruct.py` guards the same missing `maxenttk` import with a
module-level `raise` (a legitimate ImportError), which makes the two `exit()`
sites the outliers, not the convention.  So both halves of the claim stand: (a)
two rows ROADMAP marks implemented are unreachable on every backend that routes
to them, and (b) a library calls `exit()` where it must raise. I would sharpen
only the scope — `julia_live` is not affected for `CVMimag` (it raises
correctly), and `mode="ED"` is not affected for either.

**Where a fix goes**: Replace the two `except: print("Not functional yet"); exit()` blocks (`src/dmrgpy/analyticcontinuation.py:10-12` and `src/dmrgpy/distribution.py:30-33`) with a narrow `except ImportError: raise NotImplementedError(...)` naming the missing `padetk`/`maxenttk` module, matching `edtk/dynamics.py:80-83`'s existing pattern, and correct the two ROADMAP §4 rows (plus the README/user-guide mentions) to reflect that neither submode has an implementation. Actually vendoring a Padé / max-entropy reconstruction is separate feature work.

**Repro**

```
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && for s in CVMimag maxent ROOTN; do for v in 3 python; do echo -n "[$s/$v] "; MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 f_submode_exit.py $s $v 2>&1 | grep -vE "CVM in E|iteration" | tail -3 | tr '\n' '|'; echo; done; done

--- f_submode_exit.py ---
import sys, numpy as np, warnings
from dmrgpy import spinchain
n=4
sub=sys.argv[1]; ver=sys.argv[2]
if ver not in ("python","ED"): ver=int(ver)
es=np.linspace(0.1,5.,60)
sc=spinchain.Spin_Chain(["S=1/2"]*n)
if ver=="python": sc.setup_python()
elif ver!="ED": sc.setup_cpp(ver)
sc.maxm=20; sc.nsweeps=12
H=0
for i in range(n-1): H=H+sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
for i in range(n): H=H+0.3*sc.Sz[i]
sc.set_hamiltonian(H)
nm=(sc.Sz[0],sc.Sz[0])
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x,y=sc.get_dynamical_correlator(mode="DMRG",submode=sub,name=nm,es=es,delta=0.2)
    print("OK",sub,ver,"integral=",np.trapezoid(np.array(y).real,es))
except BaseException as e:
    print("EXC",sub,ver,type(e).__name__,str(e)[:80])
print("SCRIPT REACHED END")
```

### 18. On a non-Hermitian H, any `scale=` reaching `get_excited_states(n>=2)` leaks a TypeError from `mpsiram`; because `dcex` always forwards its own `scale=10.0`, the whole submode=\"EX\" implementation is unrunnable for non-Hermitian chains

- **Severity**: medium  **Class**: crash
- **Code**: `src/dmrgpy/dcex.py:56 (via src/dmrgpy/excited.py:85 -> algebra/arpacktk.py:309 mpsiram)`
- **Status**: FIXED -- excited_states_non_hermitian accepts and ignores scale, alongside the existing nkry_min/nkry_max compatibility slots.
- **Affects**: Backend-agnostic (the failing code is in the shared Python layer, not any backend). Confirmed on itensor_version=3; the same code path is taken for 2 and "python" (excited.py's non-Hermitian branch is not itensor_version-gated) and for "julia_live". ED (edtk/edchain.py -> algebra.lowest_states) is NOT affected: it has a catch-all **kwargs that swallows scale.

Public API calls hit:
(a) Many_Body_Chain.get_excited_states(n>=2, scale=...) and get_excited(n>=2, scale=...) on a non-Hermitian Hamiltonian -- confirmed crashing directly. (For itensor_version="julia_live" even n==1 takes this route, since excited.py's n==1 shortcut is gated on itensor_version in (2,3,"python").)
(b) dcex.dynamical_correlator (the submode="EX" implementation) with ANY argument list on a non-Hermitian H -- confirmed crashing when called directly.
(b) is currently masked from the public get_dynamical_correlator on itensor_version in (2,3,"python") and on mode="ED": dynamics.py:12-16 and edtk/dynamics.py:49-51 both re-route every non-Hermitian call to nonhermitian/dynamics.py's CVM-explicit resolvent before submode dispatch, so submode="EX" never reaches dcex there (verified: the public call returned successfully and printed "Non Hermitian mode in dynamical correlator"). itensor_version="julia_live" is the one public route that does reach dcex with a non-Hermitian H, since dynamics.py's julia branch only blocks KPM/CVM/TDZ -- established by reading mpsjulialive/dynamics.py:106-113, not executed (per the no-julia_live instruction).

**Expected**

ROADMAP.md Sec. 4 lists submode="EX" as "Non-Hermitian-capable; entirely
backend-agnostic (self.get_excited_states() + generic MPS algebra)". In fact
dcex.dynamical_correlator unconditionally passes scale=10.0 (its own default)
into get_cached_excited_states -> get_excited_states; for a non-Hermitian H
that routes to excited.py:85's excited_states_non_hermitian, whose **kwargs
flow straight into arpacktk.mpsiram, which has no catch-all -- the exact
failure mode excited_states_non_hermitian's own comment says it added
nkry_min/nkry_max to prevent for other kwargs. So the "EX" correlator can never
run on a non-Hermitian chain, on any backend, with any argument list, even when
called directly. The public API never surfaces this because the previous
finding silently reroutes submode="EX" to the CVM resolvent first; the two bugs
mask each other. Either get_excited_states should accept and ignore scale on
the non-Hermitian path, or dcex should not forward it there.

**Observed**

```
is_hermitian: False
get_excited_states(n=3) OK: [-1.586925+0.j       -0.974593-0.j       -0.952043+0.129151j]
get_excited_states(n=3,scale=10.0) -> TypeError: mpsiram() got an unexpected keyword argument 'scale'
dcex.dynamical_correlator (submode='EX' impl) -> TypeError: mpsiram() got an unexpected keyword argument 'scale'
```

**Reviewer** (independent, brief was to refute)

Real, and I could not explain it away. dcex.py:56 unconditionally passes its
own default scale=10.0 into get_cached_excited_states ->
self.get_excited_states(n=nex, purify=False, scale=10.0). excited.py:74-85's
get_excited_states routes a non-Hermitian H (n>=2) to
excited_states_non_hermitian, whose signature names only
recursive/maxit/ncv/nkry_min/nkry_max; scale falls into **kwargs and is handed
to arpacktk.excited_states ->
mpsiram(self,H,which=...,nev=...,wfskip=...,scale=10.0). mpsiram
(arpacktk.py:309) has no catch-all, so it raises TypeError. This is exactly the
failure mode excited_states_non_hermitian's own comment says nkry_min/nkry_max
were whitelisted to prevent -- the fix was applied to those two names only.  My
own run (4 sites, maxm=20, itensor_version=3, XXZ + 0.3j*Sz[0]):
is_hermitian: False   get_excited_states(n=3) -> OK [-1.586925, -0.974593,
-0.952043+0.129151j]   get_excited_states(n=3, scale=10.0) -> TypeError:
mpsiram() got an unexpected keyword argument 'scale'
dcex.dynamical_correlator(..., nex=3) -> same TypeError A second run confirmed
the diagnosis is exact rather than incidental: nkry_max=8/nkry_min=3 are
accepted fine, and dropping ONLY the scale= forwarding (monkeypatching
get_cached_excited_states) lets the full dcex body run to completion on the
same non-Hermitian chain, returning finite values. So nothing else in the EX
path objects to a non-Hermitian H.  Ruled out as intended behaviour: ROADMAP.md
Sec. 4 lists EX as "Non-Hermitian-capable; entirely backend-agnostic", i.e.
this contradicts the documented intent rather than implementing a documented
restriction; it is not one of CLAUDE.md's deliberately reproduced legacy bugs
(evoloperator z^3/6, "moise", "tevol_fit_td"); and it is not a
NotImplementedError raised for an unsupported combination but an internal
TypeError leaked several frames deep. Ruled out as repro flaws: deterministic
hard exception, nothing to do with DMRG convergence, degeneracy, chain length
(4 sites, above v3's 3-site floor), or comparing quantities that were never
equal.  Blast radius is genuinely limited by the masking, which is why I kept
severity at medium rather than high: (a) is a directly public, documented knob
(docs/user_guide.md tells users to retry "with a larger scale" when an excited
state is flagged unconverged) that crashes confusingly on non-Hermitian chains;
(b) makes the EX correlator dead code for non-Hermitian input, currently
reachable only via julia_live, and would become reachable on every backend the
moment the blanket non-Hermitian re-route in dynamics.py is narrowed.

**Where a fix goes**: In src/dmrgpy/excited.py, give excited_states_non_hermitian an explicitly named, ignored `scale=None` (and `noise=None`) parameter alongside the existing nkry_min/nkry_max backward-compatibility absorbers, so Hermitian-DMRG-only knobs stop leaking into mpsiram; alternatively/additionally have dcex.get_cached_excited_states forward scale= only when self.is_hermitian(self.hamiltonian).

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/r6_nhex.py:

"""submode='EX' on a non-Hermitian H: the implementation itself is broken."""
import numpy as np
from dmrgpy import spinchain, dcex
n = 4
sc = spinchain.Spin_Chain([2]*n)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
h = h + 0.3j*sc.Sz[0]
sc.set_hamiltonian(h); sc.maxm=20; sc.nsweeps=6
print("is_hermitian:", sc.is_hermitian(sc.hamiltonian))
es_, ws = sc.get_excited_states(n=3)
print("get_excited_states(n=3) OK:", np.round(np.asarray(es_),6))
try:
    sc.get_excited_states(n=3, scale=10.0)
except Exception as e:
    print("get_excited_states(n=3,scale=10.0) ->", type(e).__name__+":", e)
try:
    dcex.dynamical_correlator(sc, name=[sc.Sz[0],sc.Sz[0]], nex=3,
                              es=np.linspace(0.2,3,20), delta=0.1)
except Exception as e:
    print("dcex.dynamical_correlator (submode='EX' impl) ->", type(e).__name__+":", e)

Command:
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 r6_nhex.py
```

### 19. MultiOperator.jordan_wigner() raises AttributeError on any two-factor same-site C[i]*C[i] term, because jordanwigner.CC(i,i) returns the Python int 0 instead of a zero MultiOperator

- **Severity**: medium  **Class**: crash
- **Code**: `/u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/multioperatortk/jordanwigner.py:39 (CC), surfacing at /u/40/ladovj1/data/Documents/programs/dmrgpy/src/dmrgpy/multioperator.py:375`
- **Status**: FIXED -- jordanwigner.CC(i,i) returns an identically-zero MultiOperator rather than the int 0.
- **Affects**: Hits the backend-agnostic operator layer, so every non-ED backend: itensor_version=2, =3 and "python". Confirmed crashing public calls: vev(C[i]*C[i]), get_pairing(pairs=...) / get_correlator_spinless(name="cc") with any diagonal (i,i) pair, and set_hamiltonian/gs_energy on a Hamiltonian containing a C[i]*C[i] term. The crash is actually raised in MultiOperator.to_terms()/write() (multioperator.py:375), reachable with no chain solve at all. mode="ED" is unaffected (MO2matrix never calls jordan_wigner) and returns 0j. itensor_version="julia_live" appears unaffected statically -- mpsjulialive serializes the *raw* term list via mpo.py's text_mpo() and never calls multioperator.jordan_wigner (grep over mpsjulialive/ shows no reference); not executed, per the no-Julia instruction.

**Expected**

0 on every backend, exactly as ED returns (C_i C_i = 0 by fermionic
anticommutation). jordanwigner.CC(i,j) returns the Python int `0` when i==j
(line 40-41), so `c*jordanwigner.two_fermions(...)` in
multioperator.jordan_wigner produces a float 0.0 whose `.op` attribute does not
exist. Note the asymmetry that hides it: Cdag[i]*Cdag[i] works, but only by
accident -- CdagCdag(i,i) calls `CC(i,i).get_dagger()` which raises inside the
bare `except:` on multioperator.py:370, falling through to the brute-force
product2jw path. Failure scenario: any user assembling a full pairing/anomalous
correlation matrix over all site pairs (the natural companion to
get_correlation_matrix, and exactly the shape
fermionchaintk/staticcorrelator.get_pairing_spinless takes `pairs=` for)
crashes on the diagonal element on every DMRG backend but works under
mode="ED".

**Observed**

```
ED <C0 C0> = 0j
2 <C0 C0> RAISED AttributeError 'float' object has no attribute 'op'
3 <C0 C0> RAISED AttributeError 'float' object has no attribute 'op'
python <C0 C0> RAISED AttributeError 'float' object has no attribute 'op'

Full traceback (from s2.py in the same directory):
  File ".../dmrgpy/multioperator.py", line 375, in jordan_wigner
    else: outterms.extend(mi.op) # flatten in the transformed terms
AttributeError: 'float' object has no attribute 'op'
```

**Reviewer** (independent, brief was to refute)

Real, and reproduced verbatim with my own script under the mandated PYTHONPATH.
jordanwigner.CC(i,j) returns the bare Python int `0` when i==j
(jordanwigner.py:39-40). In multioperator.jordan_wigner() the line `mi =
c*jordanwigner.two_fermions(...)` then yields a float 0.0; that succeeds, so
the surrounding bare `except:` brute-force fallback never fires, and the *next*
statement outside the try, `outterms.extend(mi.op)` (multioperator.py:375),
raises `AttributeError: 'float' object has no attribute 'op'`.  Evidence (n=4
spinless chain, maxm=20, nsweeps=6): - ED <C0 C0> = 0j; itensor_version 2, 3
and "python" all RAISE AttributeError. - Backend-free minimal repro:
`(fc.C[0]*fc.C[0]).to_terms()` raises;
`.to_terms(jordan_wigner_transform=False)` returns fine. - Public documented
API: with a pairing Hamiltonian, `get_pairing(pairs=[(0,1)])` gives -0.1079 on
ED, v3 and python alike (so the machinery is otherwise correct), while
`get_pairing(pairs=[(0,0)])` and the full 4x4 pair list raise on v3/python and
return the correct 0 on ED. Same for get_correlator_spinless(name="cc"). - A
Hamiltonian containing `0.3*C[0]*C[0]` raises at DMRG setup on v3 but solves
fine under ED (-2.2360679...).  Boundaries I checked, which narrow the claim:
only the n==2 branch is affected. A four-factor term containing a same-site CC
pair survives (four_fermions does `0*MO`, which goes through
MultiOperator.__rmul__ and yields a correct zero operator), and other term
shapes fall into the bare `except:` -> product2jw and work. `Cdag[i]*Cdag[i]`
works only by accident: CdagCdag(i,i) calls `CC(i,i).get_dagger()` which raises
`'int' object has no attribute 'get_dagger'` *inside* the bare except, so it
falls through to product2jw and returns the identically-zero ("Adag","Adag")
term (verified: to_terms gives [((1+0j), [('Adag',1),('Adag',1)])], vev 0j on
all backends).  Ruled out as intended/flawed: not a convergence or degeneracy
artifact (it is a deterministic Python type error, independent of physics and
of maxm/nsweeps); not the <3-site v3 abort (n=4, and ED/other pairs work); not
one of CLAUDE.md's deliberately reproduced legacy bugs (evoloperator z^3/6,
"moise", "tevol_fit_td"); no ROADMAP.md row, docs/user_guide.md text or code
comment marks the diagonal of a pairing correlator as unsupported; `return 0`
has been in the file since its creation commit (8fc91c0) with no comment, and
the adjacent line was bug-fixed later (ad7766d, `-1*CdagC(j,i)` ->
`-1*CC(j,i)`), so the file is maintained rather than frozen; the existing
regression test tests/test_spinful_fermion_chain.py::test_get_pairing_matches_t
he_bare_expectation_value only uses nearest-neighbour pairs, which is why it
never hit this. Not an API misuse either: the diagonal is the natural companion
of a full anomalous correlation matrix, exactly the shape
get_pairing(pairs=...) accepts.  Blast radius: a crash, not a silent wrong
number, and the physical value is trivially 0 -- hence medium, not high.

**Where a fix goes**: In src/dmrgpy/multioperatortk/jordanwigner.py, make CC(i,i) return the identically-zero MultiOperator `obj2MO(["A",i])*obj2MO(["A",i])` instead of the int 0 (exactly the form product2jw already produces for the Cdag case, which every backend handles and evaluates to 0j); that also repairs CdagCdag(i,i) through get_dagger() rather than through the bare except. Optionally belt-and-braces: have multioperator.jordan_wigner treat a scalar `mi` as a zero term instead of calling `mi.op`.

**Repro**

```
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 f_cc_same_site.py

--- f_cc_same_site.py ---
import numpy as np
from dmrgpy import fermionchain
n=4
for v in ["ED",2,3,"python"]:
    fc = fermionchain.Fermionic_Chain(n)
    if v=="python": fc.setup_python()
    elif v!="ED": fc.setup_cpp(v)
    fc.maxm=20; fc.nsweeps=6
    H=0
    for i in range(n-1): H = H + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    fc.set_hamiltonian(H)
    mode = "ED" if v=="ED" else "DMRG"
    try: print(v,"<C0 C0> =",fc.vev(fc.C[0]*fc.C[0],mode=mode))
    except Exception as e: print(v,"<C0 C0> RAISED",type(e).__name__,e)
```

### 20. Conserved-sector mode + get_site_entropy / get_pair_entropy / get_mutual_information raise ValueError on both sector-capable backends (v3 and "python"), for spin chains under Sz and spinless fermionic chains under Nf, because reduced_dm_projective applies non-conserving ladder projectors (S+ / Cdag) to the state; get_bond_entropy works, and promote_to_dense()+promote_mps() recovers the exact values

- **Severity**: low  **Class**: missing-combo
- **Code**: `src/dmrgpy/densitymatrix.py:27 (P01 = self.Sx[k]+1j*self.Sy[k] in reduced_dm_projective)`
- **Status**: FIXED -- reduced_dm_projective evaluates <wf|Pa^dag Pb|wf> as one sandwich instead of applying a charge-changing projector to the state.
- **Affects**: itensor_version=3 and itensor_version="python" (the only two sector-capable backends), DMRG only, only while a conserved sector is set. Public API hit: Many_Body_Chain.get_site_entropy / get_pair_entropy / get_mutual_information (entropy.py site_entropy/pair_entropy/mutual_information -> densitymatrix.reduced_dm_projective), and the identical MPS.get_site_entropy/get_pair_entropy/get_mutual_information wrappers. Confirmed for Spin_Chain S=1/2 under set_conserved_sector(Sz=0) ("the term S+(3) changes ... by Sz=2") and for Fermionic_Chain under set_conserved_sector(Nf=2) ("the term F(1)*Adag(2) changes ... by Nf=1"). Unaffected: get_bond_entropy (uses _session.bond_entropy, works and is correct in-sector), mode="ED", any non-sector run, itensor_version=2 and julia_live (no sector mode at all). S=1 sites and Spinful_Fermionic_Chain raise NotImplementedError from reduced_dm_projective with or without a sector — pre-existing and not part of this finding.

**Expected**

site_entropy(2)=0.69314718 and mutual_information(1,4)=0.01422517 -- exactly
the ED values printed in the same run, since the Sz=0 sector contains the
global ground state of the even Heisenberg chain, and since a reduced density
matrix / entanglement entropy is a sector-conserving quantity. ROADMAP.md
section 2 marks "Bond entanglement entropy / reduced density matrix"
implemented on every backend with the note that entanglement.py has no backend
branching. The same run shows get_bond_entropy works fine in the sector, so
this is not a limit of sector mode but of one implementation choice:
reduced_dm_projective (densitymatrix.py:16-61) builds its single-site basis out
of the raising operator S+ = Sx+i*Sy (Cdag for a fermionic chain), and each of
those intermediate projections leaves the sector, so the guard rejects it. The
projective RDM is reachable from the public API only through
get_site_entropy/get_pair_entropy/get_mutual_information, so all three are
unusable in a sector on both sector-capable backends.

**Observed**

```
== itensor_version=3
ED reference: site_entropy(2)=0.69314718  mutual_information(1,4)=0.01422517
sector gs_energy = -2.4935771338879276
bond_entropy(2,3) in sector: 0.7113730475802916
site_entropy RAISED ValueError: Conserved-sector mode (Sz=0): the term S+(3) changes the conserved quantum numbers by Sz=2 -- every operator b
pair_entropy RAISED ValueError: Conserved-sector mode (Sz=0): the term Sz(2)*S+(5) changes the conserved quantum numbers by Sz=2 -- every oper
mutual_information RAISED ValueError: Conserved-sector mode (Sz=0): the term S+(2) changes the conserved quantum numbers by Sz=2 -- every operator b

== itensor_version=python
(identical, same three ValueErrors, same ED reference values, and the same working bond_entropy 0.7113728863069795)
```

**Reviewer** (independent, brief was to refute)

Reproduced independently with my own script (/tmp/.../scratchpad/sk03.py,
sk03b.py, sk03c.py, all run with PYTHONPATH pinned to the working tree, MKL/OMP
threads=1, n=6 spin / n=4 fermion, maxm=20, nsweeps=8). On itensor_version=3:
sector gs_energy = -2.493577133887926 (exactly the ED value
-2.4935771338879262, so the state is right and this is not an unconverged-DMRG
artifact), get_bond_entropy(2,3) = 0.71137 works in-sector, and all three of
site_entropy / pair_entropy / mutual_information raise ValueError naming S+(3),
Sz(2)*S+(5), S+(2). Identical on itensor_version="python" (same three raises,
lowercase message from pyitensor/sector.py). Same failure on a spinless
Fermionic_Chain at Nf=2, where the offending factor is Adag/Cdag instead.  Not
a repro flaw: the quantities are well defined and computable in the sector —
running the documented escape hatch promote_to_dense() + promote_mps(wf) on the
*same* sector ground state returns site=0.6931471805599452,
pair=1.3720691940332637, mi=0.014225167086626689 on v3 (and 0.6931471805599448
/ 1.372069297597303 / 0.014225063522584458 on "python"), matching the ED
reference 0.69314718 / 1.37206919 / 0.01422517 to machine precision. So the
state, the sector and the physics are all fine; only the in-sector code path
refuses.  Not intended, on three specific grounds. (1) docs/user_guide.md's
"What works under a sector" list explicitly includes "entanglement entropies"
alongside ground state, expectation values and static correlators. (2)
mpscpp3/chain_session.h:4813's reject_sector() shows that deliberate refusals
in this design are *enumerated* by name (METTS, iDMRG/VUMPS, TEBD on v3) with a
message telling the user to clear the sector; these three methods are in no
such list — they fail incidentally through the generic sector_terms() guard,
which is a correctness net for user-built operators, not a designed refusal for
entropies. (3) The raised message names S+(3) / Adag(1), operators the caller
never wrote: it is an implementation detail of reduced_dm_projective
(densitymatrix.py:27, P01 = Sx[k]+1j*Sy[k], and Cdag[k] for fermions) leaking
out. The counter-reading that user_guide's conserving-operator rule already
covers this is weak, because that rule as written and as exemplified
(sc.vev(sc.Sx[0])) governs operators the *user* builds on the chain.  Severity:
I rank it low rather than the lead's medium. It is a loud ValueError, never a
silent wrong number; it only affects an opt-in mode; and the documented
promote_to_dense()/promote_mps() route — presented in user_guide.md as exactly
the answer for quantities the conserving-operator rule rules out — recovers the
exact values on both backends. The argument for medium is the direct
contradiction with the documented "entanglement entropies work under a sector",
plus the misleading error text; the parent can re-rank on that basis. This
independently confirms item 9 of docs/audit_2026_08_hole_hunt.md's "Unverified
leads" section.

**Where a fix goes**: In densitymatrix.py::reduced_dm_projective, stop applying the raising projector to the state when self has a sector set: compute dm[a,b] = vev(Pa.get_dagger()*Pb) as one combined MultiOperator (which conserves the charge exactly when Pa and Pb carry equal charge) and fill cross-charge entries with an exact 0 without calling the backend at all — do NOT implement this as "the off-diagonals are zero", since the pair RDM's equal-charge coherences (e.g. n_up(i)*S+(j) against S+(i)*n_up(j)) are genuinely nonzero. The cheaper alternative is to have entropy.py/reduced_dm_projective work on an internally promoted dense copy of the chain and wavefunction (promote_to_dense/promote_mps) when a sector is set.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/f4_entropy.py:

import sys, numpy as np
from dmrgpy import spinchain
iv = sys.argv[1]
n = 6
sc = spinchain.Spin_Chain(["S=1/2"]*n)
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
if iv=="python": sc.setup_python()
else: sc.setup_cpp(version=3)
sc.maxm = 30; sc.nsweeps = 10
sc.set_hamiltonian(h)
wfed = sc.get_gs(mode="ED")
print("ED reference: site_entropy(2)=%.8f  mutual_information(1,4)=%.8f"
      %(sc.get_site_entropy(wfed,2), sc.get_mutual_information(wfed,1,4)))
sc.set_conserved_sector(Sz=0)      # global GS of even Heisenberg IS in Sz=0
print("sector gs_energy =",sc.gs_energy())
wf = sc.get_gs()
print("bond_entropy(2,3) in sector:",sc.get_bond_entropy(wf,2,3))
for label,fn in [("site_entropy",lambda: sc.get_site_entropy(wf,2)),
                 ("pair_entropy",lambda: sc.get_pair_entropy(wf,1,4)),
                 ("mutual_information",lambda: sc.get_mutual_information(wf,1,4))]:
    try: print(label,"=",fn())
    except Exception as ex: print(label,"RAISED %s: %s"%(type(ex).__name__,str(ex)[:110]))

cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 f4_entropy.py 3
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 f4_entropy.py python
```

### 21. maxm <= 0 is validated nowhere: itensor_version=2/3 die with an uncatchable SIGABRT, itensor_version=\"python\" silently returns the bond-dimension-1 energy (maxm<0) or ignores the setting entirely (maxm=0)

- **Severity**: low  **Class**: silent-wrong
- **Code**: `src/dmrgpy/pyitensor/chain.py:496`
- **Status**: FIXED -- Many_Body_Chain.maxm validates on assignment.
- **Affects**: itensor_version=2 and itensor_version=3 (compiled, both available here): SIGABRT/exit 134 from ITensor's hermitian.cc, no Python exception. itensor_version="python" (pyitensor): no crash, wrong-but-plausible number for maxm<0, setting silently ignored for maxm=0. mode="ED" unaffected (maxm is not used). julia_live not tested (per instructions). Reachable from any public call that pushes sweep params into the session: gs_energy()/get_gs() (groundstate.py:49,295), vev() (vev.py:31), get_excited_states() (excited.py:48), KPM correlators (kpmdmrg.py:103,194), time evolution (timedependent.py:135,213), MPS algebra (mpsalgebra.py:159,170), NH-DMRG (nhdmrg.py:88,218), METTS (vevtk/*), and CVM via the separate cvm_maxm knob (cvm.py:109), which has the same total absence of validation.

**Expected**

A ValueError naming maxm, on every backend -- which is precisely what the
infinite-chain path already does for the same knob (`vumps_ground_state: D must
be >= 1, got 0` / `got -3`, verified in scratchpad/b11.py), so the standard
exists in this codebase. Instead the finite-chain path validates maxm nowhere:
on itensor_version="python" a negative maxm produces -1.059016994375 against
the exact -1.616025403784, a perfectly ordinary-looking variational energy (it
is a valid energy of some product-ish state, not a NaN), silently ~35% off; on
the compiled backends the same input takes down the whole interpreter through
an uncatchable ITensor abort rather than raising. Note also the inconsistency
between maxm=0 (right answer, because the bond ramp overrides it) and maxm=-5
(wrong answer) on the same backend.

**Observed**

```
--- python ---
iv='python' maxm=20 (correct) e0 = -1.616025403784
iv='python' maxm=0          e0 = -1.616025403784
iv='python' maxm=-5         e0 = -1.059016994375
--- 3 ---
m > maxdim; m = 1, maxdim = 0
From line 145, file hermitian.cc

m > maxdim

m > maxdim
(process killed, SIGABRT, no Python exception)
--- 2 ---
m > maxm; m = 1, maxm = 0
From line 108, file hermitian.cc

m > maxm

m > maxm
(process killed, SIGABRT, no Python exception)
```

**Reviewer** (independent, brief was to refute)

Confirmed with my own script (4-site open Heisenberg S=1/2, nsweeps=8, ED
reference -1.616025403784438). itensor_version="python": maxm=20 ->
-1.616025403784, maxm=0 -> -1.616025403784, maxm=-5 -> -1.059016994375.
itensor_version=3 and =2: maxm=0 and maxm=-5 both terminate the interpreter
with SIGABRT (exit 134), "m > maxdim; m = 1, maxdim = 0" from ITensor's
hermitian.cc; nothing is raised, so no caller can catch it. maxm=20 is fine on
all three, and maxm=1/2 is fine on all three (maxm=1 gives -1.059016994374
everywhere), so the failure boundary is exactly maxm <= 0 and 1 is the natural
floor -- mpsalgebra.py's is_hermitian() probe already sets self.maxm = 1
internally, so 1 is a supported value, not an edge.  The claimant's stated
cause for the maxm=0 vs maxm=-5 split on pyitensor ("the bond ramp overrides
it") is wrong; I ruled it out by rerunning with bond_ramp=False and getting
identical values (0 -> -1.616025403784, -5 -> -1.059016994374). The actual
cause is pyitensor/svd.py lines 170 and 184, `_truncate(S_host, cutoff, maxdim
if maxdim else None, max(1, mindim))`: maxdim=0 is falsy and becomes None, i.e.
"no cap at all", which is why the exact answer comes back; maxdim=-5 is truthy,
so `hi = min(maxdim, n) = -5` and `j_forced` clamps the kept rank down to
mindim=1. The -1.059016994375 returned for maxm=-5 is therefore not noise or
NaN -- it is exactly the bond-dimension-1 variational optimum, a legitimate-
looking energy, delivered with no warning that the requested setting was
nonsense.  I could not explain any of this as intended. Neither backend's
set_sweep_params validates anything (mpscpp3/chain_session.h:305-311;
pyitensor/chain.py:218-222), Many_Body_Chain never checks self.maxm
(manybodychain.py:64), ROADMAP.md has no entry restricting maxm, and this is
not one of CLAUDE.md's three deliberately reproduced legacy bugs (evoloperator
z^3/6 on H2, the "moise" key, the unreachable "tevol_fit_td" branch). The
opposite convention is well established in this codebase:
pyitensor/vumps.py:710 "vumps_ground_state: D must be >= 1, got {}" for the
same physical knob on the infinite-chain path, plus metts nsamples/nt/njobs,
infinitechain.kpm_finite n_window, idmrg_window build_window n_window,
chain.nhkpm_moments n, and algebra/kpm.py's eta. CLAUDE.md also states
explicitly that ITensor's own Error() calls abort() and "nothing it raises is
catchable from Python", which is why the sector code deliberately throws
std::invalid_argument/std::runtime_error instead -- this path is exactly the
failure mode that design note exists to prevent.  Blast radius is narrow, which
is why I graded this below the lead's "medium": no valid input yields a wrong
number, no internal code path can produce maxm <= 0 (every internal assignment
I found is maxm*2, min(old_maxm,8), 5, or 1), and the trigger is a nonsensical
user value. The realistic harm is a maxm sweep or config-driven script whose
value goes to 0 or negative, which on a compiled backend kills a long batch job
with no traceback and on pyitensor returns a quietly mean-field-quality answer.

**Where a fix goes**: Add one Python-side guard raising ValueError("maxm must be >= 1, got %r") before any session call -- best placed where sweep params are pushed (a shared helper, or a property setter for Many_Body_Chain.maxm in manybodychain.py, plus the same check for cvm_maxm/kpmmaxm), so the compiled backends are stopped before ITensor's uncatchable abort; optionally mirror it in Chain::set_sweep_params as a std::invalid_argument, matching the sector code's existing anti-abort pattern.

**Repro**

```
Script /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad/repro_negmaxm.py:

import sys
import numpy as np
from dmrgpy import spinchain
n = 4
def build(iv):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=iv)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.nsweeps = 8
    return sc
iv = sys.argv[1]
if iv != "python": iv = int(iv)
s = build(iv); s.maxm = 20
print("iv=%-8s maxm=20 (correct) e0 = %.12f" % (repr(iv), s.gs_energy()))
for m in [0, -5]:
    s = build(iv); s.maxm = m
    print("iv=%-8s maxm=%-3d        e0 = %.12f" % (repr(iv), m, s.gs_energy()))

Commands (run once per backend):
cd /tmp/claude-3291340/-u-40-ladovj1-data-Documents-programs-dmrgpy/305c3d6c-b1e7-435b-9aa3-64438e7a94d9/scratchpad && MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/u/40/ladovj1/data/Documents/programs/dmrgpy/src python3 -u repro_negmaxm.py python
... && python3 -u repro_negmaxm.py 3
... && python3 -u repro_negmaxm.py 2
```

## Refuted

Reproduced, but explained as intended behaviour or as a flawed repro. Recorded so the same claim does not get re-investigated from scratch.

- **Non-Hermitian Hamiltonian: get_dynamical_correlator silently ignores submode= entirely on itensor_version 2/3/"python" -- even a nonsense submode returns a number**

  Reproduced exactly (v3 extension available, all five submodes returned the
same curve to ~2e-10, including "THIS_IS_NOT_A_SUBMODE"), but the behaviour is
documented, deliberate, and cross-backend consistent, so the claim as framed is
refuted. (1) Not "silent":
nonhermitian/dynamics.py::dynamical_correlator_cvm_explicit prints "Non
Hermitian mode in dynamical correlator" on every call; the claimant's own repro
command pipes it through `grep -v '^Non Hermitian'`, filtering out the very
announcement of the redirect. (2) Documented intent: docs/user_guide.md §6 says
verbatim that for a non-Hermitian H, "KPM is currently the only submode with a
genuine biorthogonal implementation ...; the other submodes fall back to a
correction-vector method that assumes A†=B" — i.e. exactly the observed
CVM_explicit substitution for CVM/TD/EX/maxent. (3) Not an oversight:
src/dmrgpy/edtk/dynamics.py has the identical non-Hermitian short-circuit
placed before its own submode dispatch, and its Hermitian else-branch was
deliberately upgraded from a bare `raise` to a NotImplementedError with an
explanatory comment while the non-Hermitian short-circuit immediately above it
was left standing; ED reproduces the same behaviour (mode="ED",
submode="NONSENSE" returns the same fallback curve), so this is one design
across all backends, not a v3/pyitensor quirk. (4) The claim's premise that a
working EX path is being shadowed fails empirically: calling
dcex.dynamical_correlator directly on the same non-Hermitian v3 chain crashes
with TypeError: mpsiram() got an unexpected keyword argument 'scale' (raised
inside dcex→excited.excited_states_non_hermitian→arpacktk, from dcex's own
scale=scale, a separate latent bug), so ROADMAP §4's "EX ... Non-Hermitian-
capable" note is aspirational rather than reachable today; the fallback shields
a broken path, it does not replace a working one. (5) The cited Julia-branch
comment is about not *raising* on EX/maxent in that branch; the (2,3,"python")
branch does not raise either, it redirects per the documented fallback — no
contradiction. The only defensible remnant is that a nonexistent submode string
is not validated on the non-Hermitian path, so a typo returns the fallback
instead of raising; that is consistent with the documented "all other submodes
fall back" semantics, yields exactly what the correctly-spelled submode would
have returned, and only bites for non-Hermitian H — a cosmetic input-validation
nicety, not a defect of substance.

