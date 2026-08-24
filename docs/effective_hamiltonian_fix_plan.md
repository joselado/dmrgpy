# Plan: fix `effectivehamiltonian.py`

Status: **proposed, not implemented.** Written 2026-08-24.

`src/dmrgpy/effectivehamiltonian.py` projects `self.hamiltonian` onto the
manifold spanned by the `n` lowest eigenstates and fits that projected
matrix as a linear combination of (products of) projected spin operators,
returning the fitted couplings. Public entry points:

| entry point | defined | routes to |
|---|---|---|
| `Many_Body_Chain.get_heff()` | `manybodychain.py:346` | `get_effective_hamiltonian_coefficients` |
| `Spin_Chain.get_effective_hamiltonian()` | `spinchain.py:130` | `get_effective_hamiltonian` |
| module functions, used directly by `examples/groundstate/effective_hamiltonian` | — | `get_effective_hamiltonian{,_couplings}` |

The method itself is sound and, on the case the example exercises, gives
the physically right answer: a Hubbard dimer with Peierls phase `phi`
returns an isotropic exchange of magnitude 0.0999 whose XY part is
rotated by exactly `2*phi`, which is the correct spin-twist. The problems
are around it.

## What is wrong

Verified empirically on 2- and 3-orbital Hubbard chains (ED and DMRG).

### A. Silently wrong results

**A1. A degenerate multiplet can be cut in half with no warning.** The
Hubbard dimer's low manifold is a singlet plus a 3-fold degenerate
triplet. Asking for `n=4` is correct and returns the right couplings
(~0.1). Asking for `n=2` or `n=3` cuts the triplet, so the retained
manifold is not invariant under the model's own symmetry — and the
function returns a full, plausible-looking dictionary of couplings that
is physically meaningless:

```
n=4 -> {('S^x_1','S^x_2'): 0.0808, ('S^y_1','S^y_2'): 0.0808, ('S^z_1','S^z_2'): 0.0999, ...}
n=3 -> {('S^z_1','S^z_1'): -5.5063, ('S^y_1','S^y_2'): 1.8567, ...}   # nonsense
n=2 -> {('S^x_1','S^x_2'):  2.7018, ('S^x_1','S^y_2'): 1.9863, ...}   # nonsense
```

Nothing anywhere warns. This is the most dangerous defect: the caller
chooses `n`, and a wrong choice is indistinguishable from a right one.

**A2. `get_projection_operators` hardcodes `ns = self.ns//2`.** That is
correct for `Spinful_Fermionic_Chain` (two interleaved spinless sites per
orbital) and silently wrong for every other class. On a 4-site
`Spin_Chain` it returns operators for sites 1 and 2 only — the fit is
then run against a basis that cannot represent the Hamiltonian, and
returns whatever least-squares residual minimum it lands on. Its
`mode="spin"` argument is never read.

**A3. `operators=` is ignored by `get_effective_hamiltonian_couplings`.**
Declared in the signature, never referenced; the function always calls
`get_projection_operators(self)`. A caller who supplies their own
operator basis gets a fit against a different basis, with no indication.
Same class of trap as the `apply_hamiltonian` parameter removed in
`7c81ae3`.

**A4. `tol=` is ignored by both** `get_effective_hamiltonian_couplings`
and `get_effective_hamiltonian_coefficients`. (`get_effective_hamiltonian`
does use its own `tol`, for `dict2latex`.)

### B. Broken entry points

**B1. `Spin_Chain.get_effective_hamiltonian()` raises `TypeError`
unconditionally.** `spinchain.py:132` passes `name="XX"`, which no
function downstream accepts:

```
TypeError: get_effective_hamiltonian_couplings() got an unexpected keyword argument 'name'
```

So the method is dead for the class where an effective spin Hamiltonian
is the most natural thing to ask for. (Note A2 would make it wrong even
once it runs.)

**B2. `get_heff()` returns the `NotImplemented` singleton** when
`operators is None`, rather than raising. The caller gets an object that
fails obscurely somewhere downstream instead of a clear message.

### C. Numerics

**C1. An exactly-linear least-squares problem is solved with a
random-start nonlinear optimizer.** `fit_matrix` minimizes
`mean(|h - sum_i x_i m_i|^2)` over real `x` with `scipy.optimize.minimize`
from `np.random.random(n)`. Measured against `np.linalg.lstsq` on the same
design matrix:

| | residual `||D x - y||` |
|---|---|
| `np.linalg.lstsq` | 1.4e-14 |
| current `minimize` | 3.5e-06 |

Eight orders of magnitude, on a problem with a closed-form solution. Two
successive calls on the same input differ by 2.6e-5 in the coefficients.
That interacts badly with the hard `cutoff=1e-4` on which coefficients
are kept: a coupling near the threshold can be returned by one run and
dropped by the next, changing the *key set* of the result. The RNG is
also unseeded, so nothing is reproducible run to run.

**C2. Rank deficiency is handled by a greedy pairwise-cosine filter.**
The raw candidate set is hugely underdetermined — on the dimer, 37
candidate matrices spanning a rank-16 space, i.e. 21 null directions.
`acceptable_matrix` drops any candidate whose normalized real Frobenius
overlap with an already-kept one exceeds 0.98. In the cases tested this
happened to leave exactly full column rank, but pairwise near-orthogonality
does not imply linear independence (three vectors pairwise 60° apart can
still be dependent), and when it fails nothing detects it — the solver
just returns an arbitrary point of the solution manifold. Which
candidates survive is also dict-insertion-order dependent.

**C3. The excited states are assumed exactly orthonormal.**
`get_representation` computes `h[i,j] = <psi_i|O|psi_j>` with no Gram
matrix. ED gives `|S - I| ~ 2e-16`; DMRG's overlap-penalty excited states
gave `~2e-9` on the dimer, and that degrades with system size and bond
dimension. Nothing checks it.

**C4. `method="single"` computes `P A P B P`, not `P (A B) P`.** These
differ whenever the manifold is not closed under `A`, so `"single"` and
`"full"` are genuinely different operators, not an optimization of one
another — undocumented. Products of Hermitian projected operators are
moreover not Hermitian, so a Hermitian target is being fitted over a
non-Hermitian basis with real coefficients; the anti-Hermitian parts have
to cancel among themselves rather than being excluded by construction.

### D. Output formatting

**D1. `dict2latex` emits invalid LaTeX**: `H =` and the overall prefactor
sit outside the `\[ ... \]` display block; every term is emitted with a
trailing `+`, including the last, leaving a dangling `+` before `\]`;
`cmax` can be 0 (all couplings below cutoff), giving a division by zero;
and `\[`/`\]` appear in non-raw strings.

## The plan

Five steps, each independently testable and independently useful. Steps 1
and 2 are the ones that change returned numbers.

### Step 1 — make the manifold honest (fixes A1)

In both `get_effective_hamiltonian_coefficients` and `_couplings`, request
`n+1` states instead of `n` and compare `es[n] - es[n-1]` against the
manifold's own energy spread. If the cut falls inside a degenerate
multiplet (gap below, say, `1e-6 * max(1, spread)`, with the threshold
exposed as a parameter), raise `ValueError` naming the offending
degeneracy and the nearest safe `n`. This is a hard error rather than a
warning: every returned coupling is meaningless when it trips, and the
fix is always for the caller to change `n`.

Cost: one extra excited state per call. Cheap relative to the fit.

### Step 2 — solve the fit exactly (fixes C1, C2)

Replace `fit_matrix`'s optimizer with a rank-revealing linear solve:

1. Build the real design matrix `D` (columns = `matrix2vector(m_i)`) and
   target `y = matrix2vector(h)` — `matrix2vector` already implements the
   right real inner product, `Re Tr(A^dag B)`.
2. `x, res, rank, sv = np.linalg.lstsq(D, y, rcond=None)`.
3. If `rank < D.shape[1]`, the fit is underdetermined: report it rather
   than silently returning one arbitrary solution. Two options, and this
   is the one design decision I would like confirmed — either
   (a) raise, telling the caller the basis is degenerate on this
   manifold, or (b) return the minimum-norm solution (which is what
   `lstsq` already gives) together with the null-space dimension, and
   warn. **Recommendation: (b)** — minimum-norm is a defensible canonical
   choice, and raising would break the currently-working example if
   pruning ever changes.
4. Keep `acceptable_matrix` only as a cheap pre-filter to keep `D`
   small; correctness now comes from the rank check, not from it.

This removes the RNG entirely, so results become reproducible, and drops
the residual to machine precision. `scipy.optimize` is no longer needed
here.

### Step 3 — fix the entry points and the dead parameters (fixes A2, A3, A4, B1, B2)

- `get_projection_operators`: derive the site count from the operator
  lists actually present (`len(self.Sx)`) instead of `self.ns//2`, which
  is correct for every class that has `Sx`/`Sy`/`Sz` — including
  `Spinful_Fermionic_Chain`, whose `Sx` list already has one entry per
  orbital. Drop the unused `mode` argument.
- `get_effective_hamiltonian_couplings`: honour `operators=` (use it as
  the single-operator basis when given, `get_projection_operators`
  otherwise), matching what `_coefficients` already does.
- Remove `tol` from both signatures, or thread it into the fit's
  coefficient cutoff. **Recommendation: thread it** — a coefficient
  cutoff is exactly what a caller passing `tol` would expect, and
  `fit_matrix`'s hardcoded `cutoff=1e-4` currently has no way to be set.
- `spinchain.py:132`: drop the stray `name="XX"`.
- `get_heff`/`_coefficients`: raise `ValueError` instead of returning
  `NotImplemented`.

### Step 4 — check what is being assumed (fixes C3, C4)

- Form the Gram matrix `S[i,j] = <psi_i|psi_j>` once per call. If
  `|S - I|` exceeds a tolerance, Löwdin-orthonormalize the basis
  (`S^{-1/2}`) before building any representation, so DMRG states with
  imperfect orthogonality are handled rather than silently trusted.
- Document `method="single"` vs `"full"` in the docstring as the two
  genuinely different operators they are, and state the Hermiticity
  caveat. Optionally symmetrize `"single"`'s products to
  `(AB + BA)/2` so the basis is Hermitian by construction — this changes
  results, so it should be a new opt-in `method`, not a redefinition.

### Step 5 — output and coverage (fixes D1)

- Rewrite `dict2latex` with raw strings, the prefactor inside the math
  block, `+` as a separator rather than a terminator, and a guard for
  `cmax == 0`.
- Add `tests/test_effective_hamiltonian.py`:
  - the Hubbard dimer at `phi=0` returns isotropic Heisenberg exchange,
    all three couplings equal to `4t^2/U` within tolerance;
  - at `phi != 0` the XY part is rotated by exactly `2*phi` and the total
    magnitude is unchanged (this is the physics the example demonstrates,
    currently unguarded by any test);
  - `n` cutting a degenerate multiplet raises (Step 1);
  - two successive calls return identical coefficients (Step 2);
  - `Spin_Chain.get_effective_hamiltonian()` runs (Step 3, B1) and covers
    every site (Step 3, A2);
  - `get_heff(operators=...)` on the dimer recovers the three exchange
    couplings — this path already works today and should stay working.

## Expected impact on existing results

Steps 1 and 2 change what comes back. Step 2 changes coefficients at the
~1e-5 level (the current optimizer's own error) and makes the key set
stable; Step 1 turns some currently-returned garbage into an exception.
`examples/groundstate/effective_hamiltonian` uses `n=4` on the dimer,
which is a clean cut, so it keeps working and its numbers move only in
the 5th digit. The `phi` sweep in that example should get the analytic
`2*phi` twist overlaid on the plot, so the example demonstrates a
verifiable physical result rather than just plotting whatever was fitted.
