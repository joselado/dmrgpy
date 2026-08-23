# `get_excited` on a degenerate spectrum: ED returns levels, DMRG returns states

This file replaces `known_issue_excited_states_penalty_silent_failure.md`, which
recorded a **misdiagnosis**. That document claimed `get_excited`/`get_gap` silently
return duplicate copies of the ground state on a spinful Kondo chain, reporting a gap
of ~1e-6 eV where the true gap was 3.5e-3 eV. **It was not a bug.** The retraction and
what actually causes the appearance of duplication are below, because the underlying
API asymmetry is real, is a property of this library rather than of that model, and
will mislead anyone who validates DMRG against ED the natural way.

## Retraction: the "duplicates" were a genuine degenerate multiplet

Settled by dense diagonalization of the 4096-dimensional case the original report used:

```
dense eigvalsh:  -1.539890 x3   -1.531907 x2   -1.530019 x2   -1.527754   -1.521329 x2
dmrgpy mode=ED:  -1.539890      -1.531907      -1.530019      -1.527754   ...
```

The ground state is a **three-fold degenerate S=1 triplet**, which a spin-rotation
invariant model must have. DMRG's "three copies spread over 172 ueV" were three
separately converged members of that triplet, and `E1 - E0 = 0` was the honest answer:
the gap to a degenerate partner is zero. DMRG was right.

## The durable finding: the two modes return different things

```
mode="ED"     ->  DISTINCT LEVELS      (one entry per eigenvalue)
mode="DMRG"   ->  MULTIPLET MEMBERS    (one entry per eigenstate)
```

Comparing them **index by index** makes DMRG look broken on *any* degenerate
spectrum, and index-by-index comparison against ED is the natural way to validate a
DMRG result. That is what generated the original wrong report, and it will do it
again.

Concretely, on a spectrum whose ground level is three-fold degenerate:

| index | `mode="ED"` | `mode="DMRG"` |
|---|---|---|
| 0 | ground level | triplet member 1 |
| 1 | **first excited level** | triplet member 2 |
| 2 | second excited level | triplet member 3 |
| 3 | third excited level | **first excited level** |

so `get_gap()` reads `es[1] - es[0]`, which is a level spacing under ED and a
degeneracy splitting (i.e. ~0, up to convergence error) under DMRG.

### How to compare them correctly

Identify states by a **quantum number**, not by position in the list. On a
spin-rotation invariant model, evaluate `<S^2>` on each returned state: a triplet
member reads 2, a singlet reads 0, and that distinction survives even when the
energies do not. It routinely will not: degenerate members converge to slightly
different energies (172 ueV apart in the case above), which can exceed the splitting
you are trying to resolve, so *energies alone cannot tell a fourth level from a
poorly-converged third multiplet member*.

Ask for enough levels that the multiplet plus the state you want both fit: with a
three-fold degenerate ground state, the first genuine excitation is `n=4`.

## One possibly-real, unconfirmed defect: the ARPACK tolerance

`src/dmrgpy/algebra/algebra.py` sets a module-level

```python
tol = 1e-7 # tolerancy
```

used by `lowest_states` and `lowest_energies`, i.e. the **ED path** — which is exactly
where one goes to decide whether a DMRG result is right, so a tolerance that hid
degeneracy there would be worse than one that hid it in the method under test.

It was reported that at this tolerance a three-fold multiplet collapses to a single
copy even with spare Ritz vectors, while `tol=0` recovers it at every `k >= 3`, and
that for a 13-fold multiplet even `tol=0` gives only 6 of 8 (tolerance dominating for
small multiplets, Krylov headroom for large ones).

**This has not been reproduced here.** A synthetic check — a 600-dimensional matrix
with an exactly three-fold degenerate lowest eigenvalue — recovered all three copies
at `k = 3, 4, 8, 12` at *both* `tol=1e-7` and `tol=0`. So either the effect needs the
structure of a real many-body spectrum (near-degeneracies elsewhere, sparsity) rather
than a generic random matrix, or it involves the downstream handling of degenerate
eigenvectors rather than ARPACK itself. Anyone picking this up should reproduce it on
the model below before changing the tolerance.

If it is confirmed, the fix is a parameter with a tighter default rather than a
hardcoded `tol=0`: zero makes ARPACK run to machine precision on every ED call.

## Reproduction case

A model with a genuinely degenerate ground multiplet, small enough for exact
diagonalization (N=3 cells = 6 spinful sites = 4096 states) — useful for any work on
degenerate spectra, which is why it is kept.

A 2-orbital-per-cell (c,f) chain: hoppings `ta` (c-c), `tb` (f-f) and `alpha` (c-f,
purely imaginary) between adjacent cells, Hubbard `U1` on the f orbital, a
hybridization term, and a **ferromagnetic** (`J < 0`) Kondo coupling — the sign is
what makes the ground state a high-spin multiplet.

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

# Compare BY <S^2>, not by index -- see "How to compare them correctly" above.
fc = build(3, 80)
es = fc.get_excited(n=8, mode="DMRG")
```

## How this propagated, which is the part worth remembering

The wrong diagnosis survived three independent retellings. It was reported in good
faith from a real observation; a second reader verified that *no code had changed*
since it was filed and reported the conclusion onward as current; a third repeated it
as established. No one lied and no one skipped a check — but "no code changed" and
"the diagnosis is correct" are different claims, and only the first was ever tested.
What stopped it was someone who held the project history, not anyone reasoning harder.

Note that the original document already contained the correct explanation, in its own
notes: *"the B=0 ground state may be a degenerate high-spin multiplet, which plausibly
aggravates the penalty method."* That is not an aggravating factor — it is the whole
answer, filed as a footnote to the wrong conclusion.
