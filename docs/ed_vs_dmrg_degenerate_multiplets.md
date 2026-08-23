# The ED path used to drop members of degenerate levels (fixed)

This file has been wrong twice, in opposite directions, and the corrections
are kept because the reasoning that produced each wrong version is the part
worth learning from.

## What is true now

`mode="ED"` and `mode="DMRG"` both return **one entry per eigenSTATE**. A
degenerate level appears once per member in both, they agree index by index
up to DMRG's own convergence error, and `get_gap()` correctly returns ~0
when the ground level is degenerate. Validating DMRG against ED positionally
— the natural thing to do — works.

Guarded by `tests/test_ed_degenerate_levels.py`.

## Version 1: a real observation, a wrong cause

The original report (`known_issue_excited_states_penalty_silent_failure.md`,
deleted) said `get_excited`/`get_gap` silently return duplicate copies of the
ground state on a spinful Kondo chain, reporting a gap of ~1e-6 eV where the
true gap was 3.5e-3 eV, and blamed the DMRG penalty method.

The observation was real. The cause was not. The "duplicates" were three
separately converged members of a genuine three-fold **S=1 triplet**, which a
spin-rotation invariant model must have, and `E1 - E0 = 0` was the honest
answer: the gap to a degenerate partner is zero. Settled by dense
diagonalization of the 4096-dimensional case the report used:

```
dense eigvalsh:  -1.539890 x3   -1.531907 x2   -1.530019 x2   -1.527754
dmrgpy mode=ED:  -1.539890      -1.531907      -1.530019      -1.527754
```

DMRG was right.

## Version 2: the right retraction, promoted to a wrong law

Look at that comparison again. It was read as *"ED returns distinct levels,
DMRG returns multiplet members — an API asymmetry, a property of the
library"*, and written up here as the durable finding, complete with a table
of how the two indexings drift apart and advice to compare states by `<S^2>`
instead of by position.

**That was also wrong**, and the evidence disproving it was in the same line
of output the whole time: if ED genuinely returned distinct levels, it would
be doing so at *every* size. It does not. Below `algebra.maxsize` (2000),
where ED diagonalizes densely, it has always returned every copy:

```
maxsize = 2000
N=8   dim=  256  (dense) : [-1.75 -1.75 -1.75 -1.75 -1.75 -1.75]
N=10  dim= 1024  (dense) : [-2.25 -2.25 -2.25 -2.25 -2.25 -2.25]
```

The 4096-dimensional reproduction case sat *above* that cutoff, on the ARPACK
path — and that path had a defect. There was never an asymmetry to design
around; there was one broken branch, showing up only above a size threshold
nobody had connected to the symptom.

## The actual defect, and the actual fix

`scipy.sparse.linalg.eigsh` stops once the Ritz pairs it holds have
converged, and **one copy of a degenerate eigenvalue satisfies that**. The
remaining copies are never produced. Measured on a ferromagnetic spin-1/2
Heisenberg chain (`H = -sum S_i.S_{i+1}`, ground level exactly `-(N-1)/4` and
exactly `(N+1)`-fold by SU(2)) at N=12, dim 4096:

| asked for | ground-level copies returned | should be |
|---|---|---|
| n=4  | 1  | 4  |
| n=8  | 3  | 8  |
| n=16 | 4  | 13 |

This is the ED path — the one used to decide whether a DMRG result is right —
so dropping levels here is worse than dropping them in the method under test.

**Tuning does not fix it.** Across (N=12,14) x (n=8,16), `tol`, `ncv` and
over-requesting eigenpairs each repaired some cases and broke others: `ncv=8k`
gave 8 of 8 at N=12/n=8 but 6 of 8 at N=14/n=8, and `ncv=k+40` did the
reverse. None of them is uniformly right, because none of them gives ARPACK a
*reason* to produce a second copy of an eigenvalue it has already converged.

Deflation does (`algebra._deflated_lowest_hermitian`). Each round solves for
the lowest level of

```
P h P + sigma |V><V|,        P = 1 - |V><V|
```

where `V` spans the eigenvectors found so far and `sigma` sits far above the
region of interest. A remaining copy of an already-found level is now the
*lowest* eigenvalue of that operator — the extremal case Lanczos is reliable
at — so it has to come out. Three details are load-bearing, each measured:

* accept only the lowest verified level per round (plus its exact degenerate
  partners), not every Ritz pair returned. Taking the higher ones on faith
  admits an excited state while a ground copy is still missing: that alone was
  6 of 8 instead of 8 of 8;
* accept a pair only after a real residual check against `h` — `eigsh` returns
  `k` Ritz pairs whether or not the trailing ones converged;
* re-orthogonalize against everything already found before accepting.

Verified elementwise against dense `eigvalsh` (max error ~5e-14, vectors
orthonormal, eigen-residual <1e-13) at dim 4096 for n=1,4,8,16 on a
ferromagnetic chain, an antiferromagnetic one, a random-field chain with no
degeneracy at all, and an interacting spinless-fermion chain.

Cost: one ARPACK call per level rather than one in total — n=8 takes ~0.8 s
at dim 4096 and ~4.4 s at dim 16384. Accepted deliberately; a fast wrong
spectrum is worth less than a slow right one on the reference path. The
rounds deliberately do *not* over-request eigenpairs: a round accepts only
its own lowest level, so spare Ritz pairs are pure cost once deflation is
doing the work (n=3 at dim 4096 went 1.71 s → 0.21 s on removing them, all
18 checked (model, N, n) cases still exact). Over-requesting survives only
in the non-deflated non-Hermitian branch, where it was measured to help.

**Not fixed: the non-Hermitian branch** (`slg.eigs`). The projector above
relies on eigenvectors being orthogonal, which a non-Hermitian matrix's are
not, so that branch still over-requests and truncates and can still lose
copies. It needs biorthogonal left/right deflation. Non-Hermitian ED is not
the path used to validate DMRG, which is why it was left.

## What still holds from version 2

Degenerate members converge to slightly different energies under DMRG (172
ueV apart in the case above), and that spread can exceed the splitting being
resolved, so **energies alone cannot tell a fourth level from a
poorly-converged third multiplet member**. When that matters, identify states
by a quantum number: on a spin-rotation invariant model evaluate `<S^2>` per
state (2 for a triplet member, 0 for a singlet). And ask for enough levels
that the multiplet plus the state you want both fit — with a three-fold
degenerate ground state, the first genuine excitation is `n=4`.

## Reproduction cases

The fast one, used by the test — analytic answer, no reference
diagonalization needed:

```python
from dmrgpy import spinchain
n = 12                                    # dim 4096, above maxsize
sc = spinchain.Spin_Chain(["S=1/2"]*n)
h = 0
for i in range(n-1):
    h = h - (sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
print(sc.get_excited(n=13, mode="ED"))    # must be 13 copies of -2.75
```

The original one, a 2-orbital-per-cell (c,f) chain with a **ferromagnetic**
(`J < 0`) Kondo coupling — the sign is what makes the ground state a high-spin
multiplet. Kept because it is a realistic degenerate spectrum rather than a
constructed one:

```python
import os, sys
sys.path.append(os.environ["DMRGROOT"])
import numpy as np
from dmrgpy import fermionchain

t1, t2, t3 = 0.2580984051676324, -0.05448808971925138, -0.010905376385470619
A, B = t1 + t3, t1 - t3
ta, tb, alpha = -A/2 - B/2, -A/2 + B/2, 1j*t2
U1 = 0.51336672
J, tc, mu = -U1*0.2, U1*0.2, 0.1      # note J < 0 (ferromagnetic)

def build(N, maxm):
    ni, n = 2, 2*N
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    fc.maxm, fc.nsweeps = maxm, 60
    h = 0
    for i in range(N):
        i0 = ni*i
        for (p0, p1), amp in [((0,2),ta), ((1,3),tb), ((0,3),alpha), ((1,2),alpha)]:
            ii, jj = i0+p0, i0+p1
            if jj >= n:
                continue                      # open boundary conditions
            h = h + amp*fc.Cdagup[ii]*fc.Cup[jj] + amp*fc.Cdagdn[ii]*fc.Cdn[jj]
        ic, iff = i0, i0+1
        Tup = -1j*(fc.Cdagup[iff]*fc.Cup[ic] - fc.Cdagup[ic]*fc.Cup[iff])
        Tdn = -1j*(fc.Cdagdn[iff]*fc.Cdn[ic] - fc.Cdagdn[ic]*fc.Cdn[iff])
        h = h + U1/2*(fc.Nup[iff]-0.5)*(fc.Ndn[iff]-0.5)
        h = h - tc/2*(fc.Nup[iff]-0.5)*Tdn + tc/2*(fc.Ndn[iff]-0.5)*Tup
        h = h + J*(fc.Sx[iff]*fc.Sx[ic] + fc.Sy[iff]*fc.Sy[ic] + fc.Sz[iff]*fc.Sz[ic])
    h = h + h.get_dagger()
    for i in range(N):
        h = h + mu*fc.Cdagup[ni*i]*fc.Cup[ni*i] + mu*fc.Cdagdn[ni*i]*fc.Cdn[ni*i]
    fc.set_hamiltonian(h)
    return fc

fc = build(3, 80)
es = fc.get_excited(n=8, mode="DMRG")
```

## How this propagated, which is the part worth remembering

The version-1 diagnosis survived three independent retellings. It was reported
in good faith from a real observation; a second reader verified that *no code
had changed* since it was filed and reported the conclusion onward as current;
a third repeated it as established. No one lied and no one skipped a check —
but "no code changed" and "the diagnosis is correct" are different claims, and
only the first was ever tested.

Version 2 then made the mirror-image mistake. Having correctly established
that DMRG was innocent, it explained the discrepancy by declaring the
*remaining* behaviour intended, and wrote a table, a comparison recipe and a
piece of user-facing API documentation on top of that assumption — without
ever running the one cheap check (does ED do this at a size below `maxsize`?)
that would have shown the discrepancy was a second bug. Exonerating one
component is not the same as explaining the observation.

Note also that the version-1 document already contained the correct
explanation of its own symptom, in its notes: *"the B=0 ground state may be a
degenerate high-spin multiplet, which plausibly aggravates the penalty
method."* That is not an aggravating factor — it is the whole answer, filed as
a footnote to the wrong conclusion.
