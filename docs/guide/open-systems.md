# Open systems

On this page we compute the steady state and the slowest relaxation modes
of a chain coupled to an environment, described by a Lindblad master
equation,

$$\frac{d\rho}{dt}=\mathcal L\rho=-i[H,\rho]+\sum_\mu\Big(L_\mu\rho L_\mu^\dagger
-\tfrac12\{L_\mu^\dagger L_\mu,\rho\}\Big).$$

The jump operators $L_\mu$ describe what the environment does: decay,
pumping, dephasing, dissipative hopping.

## Steady state and dissipative gap

$\mathcal L$ is linear in $\rho$, so it has eigenvalues $\lambda_k$ with
$\mathrm{Re}\,\lambda_k\le 0$. The steady state is the eigenvector with
$\lambda=0$. The slowest mode sets the relaxation time through the
dissipative gap $\Delta=-\mathrm{Re}\,\lambda_1$, meaning that a generic
state approaches the steady state as $e^{-\Delta t}$. The way we do this
is by writing $\rho$ as a matrix product state on a chain of twice the
length, one site for the ket and one for the bra index of each physical
site. Then $\mathcal L$ becomes an operator on that chain, and we take its
slowest eigenvalues from time evolution, as in
[Excited states](excited-states.md#spectrum-from-time-evolution).

As a reference, let us take an XXZ chain pumped up at the left end, pumped
down at the right end, and dephased on every site:

```python
from dmrgpy import spinchain
from dmrgpy.lindblad import Lindblad_Chain
n = 6
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version="python")
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 0.5*sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
sp = sc.Sx[0] + 1j*sc.Sy[0]          # pumps spin up at the left end
sm = sc.Sx[n-1] - 1j*sc.Sy[n-1]      # pumps spin down at the right end
jumps = [0.6*sp, 0.6*sm] + [0.4*sc.Sz[i] for i in range(n)]
lc = Lindblad_Chain(sc, jumps)
lam, res, obs = lc.get_spectrum()           # lam[0] ~ 0 is the steady state
gap = -lam[1].real                          # dissipative gap
mz = [lc.steady_state_vev(sc.Sz[i]) for i in range(n)]
```

`lam` is sorted by real part, the steady state first. `res[k]` is how
much $\lambda_k$ moves when the last snapshot is left out; only the
eigenvalues that do not move are converged. `mz` is the magnetization profile of the steady state, a
linear ramp from $+0.18$ to $-0.18$ that carries a spin current through
the chain. Here the gap is $0.0927$ against $0.0926$ from ED, and the
profile agrees with ED to $3\times10^{-4}$. The keywords are those of
`get_krylov_spectrum`, with one difference:

- `dt`: the time step (default 2.0). The evolution must run for several
  times $1/\Delta$, so `nk*dt` has to be large when the gap is small.
- `substeps`: the number of time steps per snapshot to start from
  (default 1). A strongly dissipative $\mathcal L$ makes long steps
  unstable, so a step that does not conserve the trace is redone with
  twice as many. Only the cost changes, not `dt`.
- `mode="ED"`: the dense $\mathcal L$, as a reference on a small chain.

When $\mathcal L$ conserves a quantity, for example the particle number
for number-conserving jumps, every sector has its own steady state, and
several eigenvalues come back at $\lambda\approx 0$. To stay in one sector
we start from product states of that sector:

- `start=(psi, phi)`: one local basis state per site (the occupation, for
  bosons). The evolution starts from $|\psi\rangle\langle\psi|$ and from
  the traceless $|\psi\rangle\langle\phi|$, and the sector of $\psi$ is
  kept, for example `start=([1,1,0],[0,1,1])` for two bosons on three
  sites. Modes of other sectors, which enter through small numerical
  errors, are left out of `lam`.

!!! note "Backends"
    `Lindblad_Chain` works for spin and boson chains on
    `itensor_version="python"`. Fermions are not supported.
