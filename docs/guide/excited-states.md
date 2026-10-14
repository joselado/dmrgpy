# Excited states and gaps

On this page we compute the low-energy spectrum: excited states, gaps,
non-Hermitian and generalized eigenproblems, and effective Hamiltonians of
a low-energy manifold. Here `sc` is the Heisenberg chain of [Ground state](ground-state.md).

## Low-lying spectrum

The lowest eigenvalues $E_0\le E_1\le\dots$ of $H$
and the gap $\Delta=E_1-E_0$ are

```python
es = sc.get_excited(n=4)                 # [E0, E1, E2, E3]
gap = sc.get_gap()                       # E1 - E0
es, ws = sc.get_excited_states(n=4)      # energies and states
```

The list has one entry per eigenstate, so a degenerate level appears once
for each state of the multiplet. If the ground state is degenerate,
`get_gap()` is zero, and we must ask for enough states to reach the first
excited level. When DMRG splits a multiplet slightly, identify the states
by a quantum number (for example $\langle\mathbf S^2\rangle$) and not by
their position in the list.

DMRG finds each excited state as the ground state of $H$ plus a penalty
on the states already found. This search can converge to a wrong state.
When the energy variance of a state is much larger than that of the ground
state, a warning is printed, and the state should be checked with more
sweeps, a larger `maxm` or `mode="ED"` on a smaller chain. The keywords are:

- `n`: the number of states.
- `scale`: the weight of the penalty, in units of the bandwidth of $H$
  (default 10).
- `purify=True` (default, `get_excited_states`): diagonalizes $H$ again in
  the space of the states found, which separates nearly degenerate states.

For example, an open spin-1 Heisenberg chain has four low states from the
two spin-1/2 edge modes (a singlet and a triplet), and the Haldane gap is
the gap to the next state:

```python
es = s1.get_excited(n=6)    # s1: an open spin-1 Heisenberg chain
haldane_gap = es[4] - es[0]
```

![Haldane gap and edge states of the spin-1 chain](img/excited-states_haldane.png#only-light)
![Haldane gap and edge states of the spin-1 chain](img/excited-states_haldane-dark.png#only-dark)

The figure shows the same physics on long open spin-1 chains. Here the
energies come from the ground states of three
[sectors](ground-state.md#fixed-particle-number-or-fixed-magnetization),
$S^z_\mathrm{tot}=0$, 1 and 2, which is faster and safer than
`get_excited` at these lengths. The edge splitting is $E_0(1)-E_0(0)$ and
the bulk gap is $E_0(2)-E_0(1)$, one magnon on top of the edge triplet.
On the left, the gap goes as $1/n^2$, and the fit $a+b/n^2$ on
$n=48$ to 128 gives $0.4107\,J$, against the known $\Delta=0.4105\,J$.
The inset shows that the edge splitting falls as $e^{-n/\xi}$ with
$\xi\approx6$, to below $10^{-7}$ at $n=96$, so the four edge states
become degenerate. On the right is $\langle S^z_i\rangle$ in the
$S^z_\mathrm{tot}=1$ ground state at $n=48$. Each end carries a spin 1/2
(each half sums to exactly $1/2$), spread over about 6 sites, and the
bulk is not magnetized.

## Gap to another particle number

For a conserved charge $A$, for
example the particle number $N=\sum_in_i$, `get_charge_gap(d=...)` gives
the energy to change $\langle A\rangle$ by $d$ from its ground-state value
$A_0$. It solves

$$H_\lambda=H+\lambda\big(A-A_0-d\big)^2,$$

increasing $\lambda$ until the constraint holds, and returns
$\Delta(d)=E_0(A_0+d)-E_0(A_0)$.

```python
gap2 = fc.get_charge_gap(d=2)    # fc: a fermionic chain; add two particles
```

`d=2` is the cost of adding a pair, and `d=-2` of removing one.
`Fermionic_Chain`, `Spinful_Fermionic_Chain` and the native spinful chain
have this method. The [conserved sectors](ground-state.md#fixed-particle-number-or-fixed-magnetization) give the same energies
directly, $E_0(N)$ for each $N$.

## Non-Hermitian Hamiltonians

When $H\neq H^\dagger$ (gain and loss,
complex or non-reciprocal hopping), the ground state is the eigenvalue with
the smallest real part. Eigenstates come in left and right pairs,

$$H|\psi_R\rangle=E|\psi_R\rangle,\qquad
H^\dagger|\psi_L\rangle=E^*|\psi_L\rangle,\qquad
\langle\psi_L|\psi_R\rangle=1.$$

`gs_energy()` detects a non-Hermitian $H$ and uses non-Hermitian DMRG
(NH-DMRG). Here `fc` is a fermionic chain with a non-Hermitian Hamiltonian,
for example a staggered imaginary potential $i\gamma(-1)^in_i$:

```python
e, psil, psir = fc.nhdmrg()     # E, left and right eigenvectors
e = fc.gs_energy()              # the same E; stores psir as the ground state
```

The energy is not a variational bound, so the result is checked with the
residuals $\lVert H|\psi_R\rangle-E|\psi_R\rangle\rVert$ and
$\lVert H^\dagger|\psi_L\rangle-E^*|\psi_L\rangle\rVert$. The keywords are:

- `tol`: the largest relative residual accepted (default `1e-4`).
- `ntries`: how many times to restart from a random state when the residual
  is too large (default 5).

Excited states of a non-Hermitian $H$ come from `get_excited_states(n)` and
`get_excited(n)`, which use an Arnoldi method on matrix product states. It
is less accurate than NH-DMRG, so check it against `mode="ED"` on a small
chain (`examples/non_hermitian/nhdmrg_VS_ED_VS_arnoldi`). NH-DMRG works on
`itensor_version=2`, `3`, `"python"` and `"julia_live"`.

![Skin effect of the interacting Hatano-Nelson chain](img/excited-states_nonhermitian.png#only-light)
![Skin effect of the interacting Hatano-Nelson chain](img/excited-states_nonhermitian-dark.png#only-dark)

The figure shows the many-body skin effect of an interacting
Hatano-Nelson chain at half filling, with hopping $t\pm\gamma$ and a
nearest-neighbor repulsion $V=t$. On the left is the density of the right
ground state, $\langle\psi_R|n_i|\psi_R\rangle/\langle\psi_R|\psi_R\rangle$,
from `nhdmrg()` on 16 sites. It is flat at $\gamma=0$, and charge moves to
the left end as $\gamma$ grows. On the right is the charge moved into the
left half. The dashed curves are exact, since with open ends $H$ is a
similarity transform of a Hermitian chain, and ED on 10 sites agrees with
NH-DMRG to $10^{-14}$. The problem becomes ill conditioned fast with $n$
and $\gamma$; `nhdmrg` warns when it does not converge.

## Generalized eigenvalue problem

For a Hermitian positive definite
operator $A$,

$$H|\psi\rangle=\lambda A|\psi\rangle .$$

```python
A = fc.Id + 0.5*sum(fc.N)           # for example
lam = fc.gs_energy_generalized(A)   # smallest lambda
```

DMRG solves this directly, with no inverse of $A$. A non-Hermitian $H$ is
also accepted, and then $\lambda$ is the complex eigenvalue with the
smallest real part. Not available on `itensor_version=2`, it raises an
error.

After this call the chain stores the generalized eigenstate as its ground
state, and `fc.e0` is its energy
$\langle\psi|H|\psi\rangle/\langle\psi|\psi\rangle$ (with left and right
states for a non-Hermitian $H$), while $\lambda$ is in
`fc.lam_generalized`. Every later calculation (excited states,
correlators) uses this state. So call it last, or call `fc.restart()` and
`fc.gs_energy()` to get back the true ground state.

## Effective low-energy Hamiltonians

Take the $n$ lowest eigenstates and
the projector $P=\sum_{k<n}|\psi_k\rangle\langle\psi_k|$ on them. We write
the projected Hamiltonian as a combination of projected operators $O_a$,

$$PHP\simeq\sum_aJ_a\,PO_aP,$$

and fit the couplings $J_a$ by least squares. This gives, for example, the
spin exchange of a Hubbard model at strong coupling,
$J\approx4t^2/U$:

```python
# fc: a Spinful_Fermionic_Chain(2) Hubbard dimer
J = fc.get_heff(operators=[fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1],
                           fc.Sz[0]*fc.Sz[1]], n=4)    # one J per operator

from dmrgpy import effectivehamiltonian
coef = effectivehamiltonian.get_effective_hamiltonian_couplings(fc, n=4)
latex = effectivehamiltonian.get_effective_hamiltonian(fc, n=4)
```

`get_heff` fits the operators given. `get_effective_hamiltonian_couplings`
fits all products $S^\alpha_iS^\beta_j$ of the chain spin operators and
returns a dictionary, and `get_effective_hamiltonian` returns the same fit
as a LaTeX string. The keywords are:

- `n`: the number of states in the manifold. It must not cut a degenerate
  multiplet, otherwise the fit is not meaningful and it raises an error.
- `operators`: the operators of the fit. They must be linearly independent
  in the manifold, otherwise the couplings are not unique and it raises an
  error.
- `tol`: couplings smaller than this are set to zero.
- `method`: `"single"` (default) fits $PO_aP\,PO_bP$, and `"full"` fits
  $P\,O_aO_b\,P$. They differ when the manifold is not closed under $O_a$.
- `mode="ED"`: exact eigenstates, for small systems.

See `examples/groundstate/effective_hamiltonian`, which fits the twisted
exchange of a Hubbard dimer with a Peierls phase.

![Exchange of the Hubbard dimer](img/excited-states_heff.png#only-light)
![Exchange of the Hubbard dimer](img/excited-states_heff-dark.png#only-dark)

The figure shows the exchange $J$ that `get_heff` fits on a Hubbard dimer
(`n=4`, `mode="ED"`, with $20(N-2)^2$ added to $H$ to keep the four
states at half filling). The three couplings come out equal, as for an
isotropic exchange. At large $U$ the fitted $J$ joins the singlet-triplet
splitting $E_T-E_S$ and $4t^2/U$. At small $U$ it is larger, because the
singlet has doubly occupied weight and the projected spin operators are
shorter than spin 1/2; the fit is then exactly
$(E_T-E_S)/(1/4+3w/4)$ (dashed), with $w$ the singly occupied weight.
