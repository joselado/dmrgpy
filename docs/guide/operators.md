# Hamiltonians and operators

A Hamiltonian or an observable is written as its operator expression, with
the operator lists of [Models](models.md). For the spin-1/2 Heisenberg chain,

$$H=J\sum_{i=1}^{n-1}\mathbf S_i\cdot\mathbf S_{i+1}
=J\sum_{i=1}^{n-1}\left(S_i^xS_{i+1}^x+S_i^yS_{i+1}^y+S_i^zS_{i+1}^z\right),$$

```python
from dmrgpy import spinchain
n = 30
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
```

Any sum of products of local operators is allowed, with any range and
complex coefficients. This covers single-ion anisotropy $D(S^z_i)^2$,
Dzyaloshinskii-Moriya terms $\mathbf D\cdot(\mathbf S_i\times\mathbf S_{i+1})$,
biquadratic exchange $(\mathbf S_i\cdot\mathbf S_{i+1})^2$, hopping
$t\sum_i(c_i^\dagger c_{i+1}+\mathrm{h.c.})$, the Hubbard term
$U\sum_in_{i\uparrow}n_{i\downarrow}$, and so on. [Cookbook](cookbook.md) gives complete
models. `h.get_dagger()` returns $h^\dagger$, so `h + h.get_dagger()` adds
the Hermitian conjugate. It is best to write $H$ in units where its largest
coefficient is of order 1: a few routes, such as real-time evolution, use
absolute thresholds.

## Shortcuts on fermionic chains

These methods take a function `f(i,j)`
of two sites and add the whole family of terms:

- `fc.set_hoppings(f)`: $\sum_{ij}f(i,j)\,c_i^\dagger c_j$. On a
  `Spinful_Fermionic_Chain` the indices run over the spinless sites
  (up and down alternate).
- `fc.set_hoppings_spinful(f)`: the same hopping for both spins, with `f`
  indexed by orbital. Spinful chains.
- `fc.set_hubbard(f)`: $\sum_{ij}f(i,j)\,n_in_j$, with
  $n_i=n_{i\uparrow}+n_{i\downarrow}$ on spinful chains. Note that the
  $i=j$ term is $f(i,i)\,n_i^2=f(i,i)(n_i+2n_{i\uparrow}n_{i\downarrow})$.

```python
from dmrgpy import fermionchain
fc = fermionchain.Fermionic_Chain(4)
fc.set_hoppings(lambda i,j: 1.0 if abs(i-j)==1 else 0.0)
fc.set_hubbard(lambda i,j: 2.0 if abs(i-j)==1 else 0.0)
```

## Operators

On an operator `A`, `A.simplify()` collects equal terms,
and `A.is_hermitian()`, `A.is_antihermitian()` and `A.is_zero()` answer by
reading the terms. A `True` answer is exact. A `False` answer can be wrong
when the property depends on an identity of the local space, such as
$S^xS^x=1/4$ for spin 1/2. `sc.is_hermitian(A)` does not have this
limitation, since it also checks numerically on the chain.

`sc.toMPO(A)` builds the matrix product operator of `A` once. Two such
operators can be added, subtracted and multiplied by numbers
(`A + 2*B - A`). Not available on `"julia_live"`. To combine Hamiltonians,
it is simpler to sum the operator expressions before calling
`set_hamiltonian`.

## States and operators

The chain also gives the basic operations on
states. Here `wf` is a state, for example `wf = sc.get_gs()`, and `A` an
operator:

| Call | Returns |
|---|---|
| `sc.overlap(a, b)` | $\langle a\vert b\rangle$ |
| `sc.aMb(a, A, b)` | $\langle a\vert A\vert b\rangle$ |
| `sc.applyoperator(A, wf)`, or `A*wf` | $A\vert \psi\rangle$ |
| `sc.applyinverse(A, wf)` | $A^{-1}\vert \psi\rangle$, solved iteratively |
| `sc.exponential(A, wf)` | $e^{A}\vert \psi\rangle$ |
| `sc.summps(a, b)` | $\vert a\rangle+\vert b\rangle$ |
| `sc.scale_mps(x, wf)`, or `wf*x` | $x\vert \psi\rangle$ |
| `sc.trace(A)`, `sc.inverse_trace(A)` | $\mathrm{Tr}\,A$, $\mathrm{Tr}\,A^{-1}$ |
| `sc.operator_norm(A)` | an estimate of $\lVert A\rVert$ |
| `sc.is_zero_operator(A)` | whether $A=0$ |
| `wf.dot(wf2)`, `wf.norm()` | $\langle\psi\vert \psi_2\rangle$, $\lVert\psi\rVert$ |
| `wf.normalize()` | $\vert \psi\rangle/\lVert\psi\rVert$ |
| `wf.get_conjugate()` | the complex conjugate state |

A state remembers the backend that made it (an MPS for DMRG, a vector for
`mode="ED"`), and these operations run on that backend. Give `mode="ED"`
to `overlap` and `aMb` when the states come from ED. On DMRG,
`exponential` needs a Hermitian or an anti-Hermitian `A`, so
`sc.exponential(-1j*t*h, wf)` is a real-time step. For a general `A` use
`mode="ED"` on a small chain.
