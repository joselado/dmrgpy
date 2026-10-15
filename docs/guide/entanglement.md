# Entanglement

On this page we measure the entanglement of a state. Here `wf` is a
ground state, for example `wf = sc.get_gs()` for the Heisenberg chain of
[Ground state](ground-state.md).

## Entanglement entropy of a cut

We cut the chain between sites $i$ and
$i+1$. The Schmidt decomposition of the state is
$|\psi\rangle=\sum_\alpha\lambda_\alpha|\alpha\rangle_L|\alpha\rangle_R$,
and the von Neumann entropy of the left part is

$$S_i=-\sum_\alpha\lambda_\alpha^2\log\lambda_\alpha^2
=-\mathrm{Tr}\,\rho_L\log\rho_L,\qquad
\rho_L=\mathrm{Tr}_R|\psi\rangle\langle\psi|.$$

```python
s = wf.get_bond_entropy(i, i+1)     # cut between sites i and i+1
s = wf.get_entropy(b)               # the same, cut between sites b-1 and b
```

In a gapped chain $S_i$ saturates with the size of the block, and in a
critical chain it grows as $\log\ell$ (see the figure in
[Central charge](#central-charge)). It also tells us the bond dimension
we need, roughly $\chi\gtrsim e^{S}$. The bond entropy needs an MPS, so it
is not available with `mode="ED"`.

## Reduced density matrix

The reduced density matrix of one site,

$$\rho_i=\mathrm{Tr}_{j\neq i}|\mathrm{GS}\rangle\langle\mathrm{GS}|,$$

contains every local observable of that site, and its eigenvalues give the
site entropy of the next section.

```python
rho = sc.get_rdm(i=0)     # d x d matrix, d the local dimension
```

It is computed from the MPS, so it is not available with `mode="ED"`,
where it raises an error.

## Entropy of one site or two sites

`wf.get_site_entropy(i)` and
`wf.get_pair_entropy(i, j)` give the entropy of the reduced density matrix
of site $i$, or of the two sites $i$ and $j$ (not necessarily neighbors),
with the rest of the chain. They are available for spin-1/2 chains and
spinless fermionic chains. With `mode="ED"`, call them on the chain:
`sc.get_site_entropy(wf, i)`, `sc.get_pair_entropy(wf, i, j)`.

## Mutual information

For two sites $i$ and $j$,

$$I(i,j)=S_i+S_j-S_{ij},$$

```python
mi = wf.get_mutual_information(i, j)
```

$I(i,j)$ measures all the correlations, classical and quantum, between the
two sites, so it shows which sites are coupled in the state even when no
term of $H$ couples them. Same availability as the site entropy.

## Central charge

For a critical open chain of length $L$, the entropy
of a cut at position $\ell$ follows the Calabrese-Cardy formula,

$$S(\ell)=\frac c6\log\!\left[\frac{2L}{\pi}\sin\frac{\pi\ell}{L}\right]+\text{const}.$$

The central charge $c$ identifies the conformal field theory, for example
$c=1/2$ for the critical Ising chain and $c=1$ for the XX and Heisenberg
chains. `get_CFT_central_charge()` fits this formula to the entropy of
every bond:

```python
c = wf.get_CFT_central_charge()
```

On an open chain the entropy also has a term that alternates between even
and odd bonds, and the fit does not remove it. So the result converges
slowly with $L$: for the spin-1/2 XX chain it gives about 0.54 at $L=20$
and 0.74 at $L=40$. Use long chains, or fit the bond entropies yourself
with only even (or only odd) $\ell$.

![Entanglement entropy and central charge](img/entanglement_central_charge.png#only-light)
![Entanglement entropy and central charge](img/entanglement_central_charge-dark.png#only-dark)

The figure shows open chains with $L=100$ sites. On the left is the
entropy of every cut. For the XX chain the DMRG entropies lie on the
exact free-fermion values (circles). The XX, Heisenberg and critical
Ising chains follow an arch. For the XX and Heisenberg chains the
even-odd oscillation of the open chain sits on top of it. The gapped spin-1 chain is flat in the bulk: its entropy
saturates. On the right, the even bonds lie on straight lines against the
Calabrese-Cardy variable. The XX and Heisenberg lines are parallel, with
$c\approx1$, and the Ising line has half their slope, $c\approx1/2$. At
$L=100$ the fits give 1.07, 1.05 and 0.52; the exact free-fermion
entropies of the XX chain give the same 1.07, so the remaining
difference is a finite-size effect.

??? example "Code of this figure (docs/figures/entanglement_central_charge/)"

    ```python title="compute.py"
    --8<-- "entanglement_central_charge/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "entanglement_central_charge/plot.py"
    ```

## Single-particle correlation matrix

For fermions, the correlation
matrix is

$$C_{ij}=\langle c_i^\dagger c_j\rangle .$$

Its eigenvalues $n_\alpha\in[0,1]$ are the occupations of the natural
orbitals. For a single Slater determinant every $n_\alpha$ is 0 or 1.
Interactions move them away from 0 and 1, and the correlation entropy
measures how far the state is from a single Slater determinant:

$$S=-\sum_\alpha\Big[n_\alpha\log n_\alpha+(1-n_\alpha)\log(1-n_\alpha)\Big].$$

Here `fc` is a fermionic chain with its Hamiltonian set.

```python
cm = fc.get_correlation_matrix()          # C_ij
na = fc.get_correlation_eigenvalues()     # n_alpha
s = fc.get_correlation_entropy()          # S
orb = fc.get_correlated_orbitals()        # |orbital|^2, most correlated first
rho = fc.get_correlated_density(n=2)      # weight of the 2 most correlated orbitals
sd = fc.get_gs().get_correlation_entropy_density()   # S resolved by site
```

The most correlated orbitals are those with $n_\alpha$ closest to $1/2$.
The four-point correlators and the two-particle density matrix are on
the [Correlators and distributions](correlators.md) page. The keywords are:

- `T`: temperature of `get_correlation_matrix` (default 0). Finite $T$ is
  computed with exact eigenstates, for small chains.
- `dmmode`: how $C_{ij}$ is computed (`"fast"`, `"full"`, `"explicit"`,
  `"simple"`). The default chooses the fastest one for the state.
  `"explicit"` works on every backend.

![Natural orbitals of the Hubbard chain](img/entanglement_natural_orbitals.png#only-light)
![Natural orbitals of the Hubbard chain](img/entanglement_natural_orbitals-dark.png#only-dark)

The figure shows the half-filled Hubbard chain with $n=12$ sites and
hopping $t=1$. On the left are the occupations $n_\alpha$ for several
values of $U$. At $U=0$ the ground state is one Slater determinant, so
half of the occupations are 1 and the other half are 0. As $U$ grows,
the step melts and all occupations move towards $1/2$. On the right is
the correlation entropy. It starts at zero and grows towards its maximum
$2n\log 2$. For a chain with $n=4$ the DMRG values (squares) lie on the
ED curve (dashed line).

??? example "Code of this figure (docs/figures/entanglement_natural_orbitals/)"

    ```python title="compute.py"
    --8<-- "entanglement_natural_orbitals/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "entanglement_natural_orbitals/plot.py"
    ```
