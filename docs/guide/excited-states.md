# Excited states and gaps

On this page we compute the low-energy spectrum: excited states, gaps,
the degeneracy of the ground state, non-Hermitian and generalized
eigenproblems, and effective Hamiltonians of a low-energy manifold. Here `sc` is the Heisenberg chain of [Ground state](ground-state.md).

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

??? example "Code of this figure (docs/figures/excited-states_haldane/)"

    ```python title="compute.py"
    --8<-- "excited-states_haldane/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "excited-states_haldane/plot.py"
    ```

## Ground-state degeneracy

A degenerate ground state signals a broken symmetry or a topological
phase. Examples are the two ordered states of an Ising chain and the four
edge states of an open Haldane chain, whose splitting is shown in
[Low-lying spectrum](#low-lying-spectrum). We count the levels close to $E_0$,

$$g\approx\sum_i\exp\!\left[-\left(\frac{(E_i-E_0)^2}{\delta}\right)^2\right],$$

adding excited states $E_i$ until the count stops growing:

```python
g = sc.get_gs_degeneracy()             # close to an integer
g = sc.get_gs_degeneracy(delta=1e-4)   # a stricter window
```

A level counts when $|E_i-E_0|$ is smaller than about $\sqrt{\delta}$,
which is 0.1 at the default `delta=1e-2`. A result that is not close to
an integer means that a level sits at the edge of the window. On a finite
chain a degeneracy is often split by an amount that decays with the
length, so choose `delta` relative to that splitting. `mode="ED"` uses
exact diagonalization.

The states of the manifold are `sc.get_gs_manifold(n=2, tol=1e-3)`, the
eigenstates within `tol` of $E_0$ (more states are searched when all `n`
turn out degenerate), and `sc.vev(A, wf=w)` measures an operator in any
of them. Here `sc` is an Ising chain
$H=-\sum_iS^z_iS^z_{i+1}-0.2\sum_iS^x_i$ in its ordered phase:

```python
wfs = sc.get_gs_manifold(n=2, tol=1e-3)         # two states
M = sum(sc.Sz)
print([sc.vev(M, wf=w).real for w in wfs])      # 0 and 0
print(sc.aMb(wfs[0], M, wfs[1]).real)           # modulus 4.82 on 10 sites
up = sc.summps(wfs[0], wfs[1]).normalize()      # the ordered states
dn = sc.summps(wfs[0], -1*wfs[1]).normalize()
print(sc.vev(M, wf=up).real, sc.vev(M, wf=dn).real)   # 4.82 and -4.82, in some order
```

The eigenstates are the symmetric and antisymmetric combinations of the
two ordered states, so the magnetization vanishes in each of them and
appears in the off-diagonal element. The ordered states are their sum
and difference. The [dynamical correlator](correlator-methods.md#degenerate-ground-states)
with `mode="ED"` averages over such a manifold by itself.

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
and $\gamma$; `nhdmrg` warns when it does not converge, and again each
time the unconverged state is used.

??? example "Code of this figure (docs/figures/excited-states_nonhermitian/)"

    ```python title="compute.py"
    --8<-- "excited-states_nonhermitian/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "excited-states_nonhermitian/plot.py"
    ```

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
error. It also raises if $A$ is not Hermitian, or if $\langle\psi|A|\psi\rangle$
is not positive during the solve.

After this call the chain stores the generalized eigenstate as its ground
state, and `fc.e0` is its energy
$\langle\psi|H|\psi\rangle/\langle\psi|\psi\rangle$ (with left and right
states for a non-Hermitian $H$), while $\lambda$ is in
`fc.lam_generalized`. Every later calculation (excited states,
correlators) uses this state. So call it last, or call `fc.restart()` and
`fc.gs_energy()` to get back the true ground state.

## Effective low-energy Hamiltonians

An effective Hamiltonian acts in a small model space and has the same
low-energy spectrum as $H$. For a Hubbard chain at strong coupling the
model space is the states with one electron per site, and the effective
model is a spin model. We take the $n$ lowest eigenstates
$|\psi_k\rangle$, project them on the model space with its projector $Q$,
and orthonormalize them (Löwdin), which gives $|\tilde k\rangle$. The
effective Hamiltonian is

$$H_\mathrm{eff}=\sum_{k<n}E_k|\tilde k\rangle\langle\tilde k|
\simeq\sum_aJ_a\,O_a,$$

and the couplings $J_a$ of the operators $O_a$ are fitted by least
squares, with the $O_a$ represented on the $|\tilde k\rangle$. This is the
des Cloizeaux effective Hamiltonian. It has the eigenvalues $E_k$, so on a
Hubbard dimer the fitted exchange is the singlet-triplet splitting
$E_T-E_S$, which tends to $4t^2/U$ at strong coupling:

```python
# fc: a Spinful_Fermionic_Chain(2) Hubbard dimer, with 20*(N-2)^2 added
# to H to keep the four lowest states at half filling
from dmrgpy import effectivehamiltonian
Q = effectivehamiltonian.single_occupancy_projector(fc)
J = fc.get_heff(operators=[fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1],
                           fc.Sz[0]*fc.Sz[1]], n=4, projector=Q)

coef = effectivehamiltonian.get_effective_hamiltonian_couplings(fc, n=4,
                                                               projector=Q)
latex = effectivehamiltonian.get_effective_hamiltonian(fc, n=4, projector=Q)
```

`get_heff` fits the operators given and returns one $J$ per operator.
`get_effective_hamiltonian_couplings`
fits all products $S^\alpha_iS^\beta_j$ of the chain spin operators and
returns a dictionary, and `get_effective_hamiltonian` returns the same fit
as a LaTeX string. The keywords are:

- `n`: the number of states. It must not cut a degenerate multiplet, and
  it must not be larger than the dimension of the model space. Otherwise
  it raises an error.
- `projector`: the projector on the model space, as an operator or a list
  of commuting operators. `single_occupancy_projector(fc)` gives the states
  with one electron per site. Without it, the model space is the span of
  the $n$ eigenstates themselves. This is right for a spin chain, but for
  a Hubbard chain the spin operators then miss the doubly occupied weight,
  and $J$ is too large at small $U$.
- `operators`: the operators of the fit. They must be linearly independent
  in the model space, otherwise the couplings are not unique and it raises
  an error.
- `tol`: couplings smaller than this are set to zero.
- `method`: `"single"` (default) fits $QO_aQ\,QO_bQ$, and `"full"` fits
  $Q\,O_aO_b\,Q$. They differ when the model space is not closed under
  $O_a$.
- `mode="ED"`: exact eigenstates, for small systems.

See `examples/groundstate/effective_hamiltonian`, which fits the twisted
exchange of a Hubbard dimer with a Peierls phase.

![Exchange of the Hubbard dimer](img/excited-states_heff.png#only-light)
![Exchange of the Hubbard dimer](img/excited-states_heff-dark.png#only-dark)

The figure shows the exchange $J$ that `get_heff` fits on a Hubbard dimer
(`n=4`, `mode="ED"`, with $20(N-2)^2$ added to $H$ to keep the four
states at half filling). The three couplings come out equal, as for an
isotropic exchange. With the projector, $J$ is the splitting $E_T-E_S$
(dashed) at every $U$, and it joins $4t^2/U$ at large $U$. Without the
projector, $J$ is too large at small $U$: the singlet has doubly occupied
weight, where $\mathbf S_0\cdot\mathbf S_1=0$, and the fit gives
$(E_T-E_S)/(1/4+3w/4)$, with $w$ the singly occupied weight. This is 40%
above the splitting at $U=t$.

??? example "Code of this figure (docs/figures/excited-states_heff/)"

    ```python title="compute.py"
    --8<-- "excited-states_heff/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "excited-states_heff/plot.py"
    ```
