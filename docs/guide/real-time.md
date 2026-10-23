# Real-time evolution

On this page we evolve a state in real time,
$|\psi(t)\rangle=e^{-iHt}|\psi(0)\rangle$, and measure an observable along
the way. The typical case is a quantum quench: we prepare the ground state
of a Hamiltonian $H_0$, change the Hamiltonian to $H_1$ at $t=0$, and follow
how observables oscillate, relax or thermalize. Here `h0` and `h1` are two
Hamiltonians of the chain `sc`, for example a staggered field and the
Heisenberg chain of [Ground state](ground-state.md).

```python
from dmrgpy import timedependent
sc.set_hamiltonian(h0)
wf0 = sc.get_gs()                          # ground state of H0
sc.set_hamiltonian(h1)                     # quench to H1
(ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0],
                                            nt=200, dt=0.02, wf=wf0)
```

This gives $\langle\psi(t)|S^z_0|\psi(t)\rangle$ at the times
$t=0,\,dt,\,2dt,\dots$. The keywords are:

- `nt`, `dt`: the number of time steps (default 1000) and the time step
  (default 0.01). A smaller `dt` is more accurate.
- `wf`: the initial state. Without it, the ground state of the chain is
  used.
- `h`: the Hamiltonian of the evolution, if not the one of the chain.
- `mode="ED"`: exact evolution, the reference for small chains.
- `return_wf=True` (DMRG): also return the final state as a third output.

The accuracy of the DMRG evolution is set by `sc.maxm` and `dt`. The
entanglement grows with time after a quench, so the bond dimension needed
grows too. When the result stops agreeing with a larger `maxm`, the
simulation has reached its maximum time.

![Light cone of a domain wall in the XX chain](img/real-time_light_cone.png#only-light)
![Light cone of a domain wall in the XX chain](img/real-time_light_cone-dark.png#only-dark)

The figure shows a quench in the XX chain,
$H=\sum_i(S^x_iS^x_{i+1}+S^y_iS^y_{i+1})$, with $n=40$. The initial state
is a domain wall, all spins up on the left half and all down on the right
half. Left: $\langle S^z_x(t)\rangle$ from TDVP. The magnetization moves
only inside a cone of slope one (dashed lines), since the fastest
excitations have velocity $v=J$. Right: the profile at three times. The
model maps to free fermions, so the exact profile is known (dashed lines).
TDVP agrees with it to $10^{-4}$. At long times the profile approaches
$-\frac1\pi\arcsin(x/t)$ (dotted line).

??? example "Code of this figure (docs/figures/real-time_light_cone/)"

    ```python title="compute.py"
    --8<-- "real-time_light_cone/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "real-time_light_cone/plot.py"
    ```

## Local quench

`evolution_ABA` first applies an operator $A$ to the
state, for example a spin flip or an added electron, and then measures $B$:

$$\langle\psi|A^\dagger e^{iHt}Be^{-iHt}A|\psi\rangle.$$

Here we flip up the spin at the edge of the Heisenberg chain and follow
how the magnetization arrives at site 3:

```python
Sp0 = sc.Sx[0] + 1j*sc.Sy[0]               # S^+ on site 0
(ts, sz3) = timedependent.evolution_ABA(sc, A=Sp0, B=sc.Sz[3], nt=200, dt=0.02)
```

It takes the same keywords as `evolve_and_measure`, and starts from the
ground state of the chain unless `wf` is given. Note that $A$ is applied as
it is, so the state $A|\psi\rangle$ is not normalized.

## Evolved states

To get the states themselves, `sc.exponential(-1j*t*h,
wf)` returns $e^{-iHt}|\psi\rangle$ for one time ([Hamiltonians and operators](operators.md)), and
`timeevolution.evolve_WF(h, wf, ts=ts)` returns a list of states, one for
each time in `ts`. Note that `evolve_WF` applies $e^{+iht}$, so pass `-h`
for the physical evolution. With the states, we compute any overlap, for
example the return amplitude $\langle\psi(0)|\psi(t)\rangle$.

## Choosing the propagator

`sc.tevol_method` selects the algorithm:

- `"TDVP"` (default): two-site TDVP. It works for any Hamiltonian and
  grows the bond dimension as needed. On `itensor_version=2` it is
  replaced by `"MPO"`.
- `"TDVP_GSE"`: one-site TDVP with global subspace expansion (Yang and
  White, [PRB 102, 094315 (2020)](https://arxiv.org/abs/2005.06104)). Use
  it when the initial state has a small bond dimension, for example a
  product state. The bond dimension grows only during the first
  `sc.tdvp_gse_sweeps` steps (default 3); the other settings are
  `sc.tdvp_gse_krylov_order` (default 3) and `sc.tdvp_gse_cutoff`
  (default $10^{-8}$).
- `"TEBD"`: second-order Trotter evolution with two-site gates. Faster
  than TDVP, but only for nearest-neighbor Hamiltonians; otherwise it
  raises an error. A `Spinful_Fermionic_Chain_Native` Hubbard chain is
  nearest-neighbor, the interleaved `Spinful_Fermionic_Chain` is not.
- `"AUTO"`: `"TEBD"` when the Hamiltonian is nearest-neighbor, `"TDVP"`
  otherwise.
- `"MPO"`: a second-order Taylor expansion of $e^{-iH\,dt}$ applied as an
  MPO. The only option on `itensor_version=2`, and less accurate than TDVP
  at the same bond dimension.

```python
sc.tevol_method = "TEBD"
```

The same propagator is used by the `"TD"` and `"TDZ"` dynamical
correlators ([Choosing a correlator method](correlator-methods.md)). See `examples/time_evolution/time_evolution_quench`.
