# Mean-field decoupling

For a spin Hamiltonian with exchange couplings $J_{ij}$, the
self-consistent mean-field (Weiss-field) decoupling replaces the exchange
by a field on each site,

$$H_{\rm MF}=p\,H_{\rm exch}+H_{\rm onsite}+\sum_i\mathbf h_i\cdot\mathbf S_i,
\qquad\mathbf h_i=(1-p)\sum_jJ_{ij}\langle\mathbf S_j\rangle.$$

We solve $H_{\rm MF}$, compute the new $\langle\mathbf S_i\rangle$, rebuild
$\mathbf h_i$, mix the old and new values, and repeat until the
magnetization stops changing. With $p=0$ this is pure mean-field theory.
With $0<p<1$ a fraction $p$ of the exchange is kept and treated exactly,
together with the field. With $p=1$ it is the original model. The
couplings are read from the Hamiltonian of the chain, so any one- and
two-site spin Hamiltonian works, including anisotropic and long-range
exchange. Here `sc` is the Heisenberg chain of [Ground state](ground-state.md):

```python
from dmrgpy import meanfield
n = sc.ns
m0 = [[0, 0, 0.5*(-1)**i] for i in range(n)]       # staggered start
mf = meanfield.spinchain_meanfield(sc, p=0.0, m0=m0, mix=0.5)
m = mf.get_magnetization()       # self-consistent <S_i>, a 3 x n array
e = mf.gs_energy()               # mean-field energy
```

The function returns a new chain whose Hamiltonian is the converged
$H_{\rm MF}$. This Hamiltonian includes the constant
$-(1-p)\sum_{ij}J_{ij}\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle$,
so the exchange is not counted twice. At $p=0$ its `gs_energy()` is the
expectation value of the original Hamiltonian in the mean-field state.
The keywords are:

- `p`: the fraction of the exchange kept exactly.
- `m0`: the starting magnetization, one 3-vector per site. Random by
  default. For an antiferromagnet start from a staggered pattern, since
  simple mixing can oscillate otherwise.
- `mix`: the weight of the new magnetization in each step (default 0.9).
  Lower it if the loop oscillates.
- `maxerror`, `maxite`: the convergence tolerance and the maximal number
  of iterations. After `maxite` the loop stops with a warning.
- `mixmode`: `"default"` (linear mixing) or `"broyden"`.
- `mode="ED"`: solve each step with exact diagonalization.

`meanfield.decompose_spin_hamiltonian(sc)` returns the onsite fields
`b[i,a]`, the exchange `J[i,j,a,b]` and the constant read from the
Hamiltonian. A Hamiltonian that is not of this form, for example with a
three-site term, raises an error.
