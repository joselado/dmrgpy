# Cookbook

The first part writes out the Hamiltonian of the common models, ready for
`set_hamiltonian`. The second part gives short recipes for calculations
that combine several pages of the guide. Every chain is built as in
[Models](models.md), `n` is its number of sites, and the operators are
those of its class. On a fermionic chain, `h + h.get_dagger()` adds the
Hermitian conjugate of every hopping and pairing term, so we write each
such term once.

## Models

### XXZ chain

$$H=J\sum_i\left(S^x_iS^x_{i+1}+S^y_iS^y_{i+1}+\Delta\,S^z_iS^z_{i+1}\right)$$

```python
sc = spinchain.Spin_Chain(["S=1/2"]*n)
Delta = 0.5
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + Delta*sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
```

$\Delta=1$ is the Heisenberg chain, `sc.SS(i, i+1)`; $\Delta=0$ is the XX
chain, which maps to free fermions; $|\Delta|>1$ is gapped.

### Transverse-field Ising chain

$$H=-J\sum_iS^z_iS^z_{i+1}-h\sum_iS^x_i$$

```python
sc = spinchain.Spin_Chain(["S=1/2"]*n)
J, hx = 1.0, 0.5
h = -J*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1)) - hx*sum(sc.Sx)
```

With spin-1/2 operators the critical point is at $h=J/2$; with Pauli
matrices ($4S^zS^z$, $2S^x$) it is at $h=J$.

### J1-J2 chain

$$H=J_1\sum_i\mathbf S_i\cdot\mathbf S_{i+1}+J_2\sum_i\mathbf S_i\cdot\mathbf S_{i+2}$$

```python
sc = spinchain.Spin_Chain(["S=1/2"]*n)
J1, J2 = 1.0, 0.4
h = J1*sum(sc.SS(i, i+1) for i in range(n-1)) + J2*sum(sc.SS(i, i+2) for i in range(n-2))
```

For $J_2/J_1>0.2411$ the chain is dimerized and gapped.

### Dzyaloshinskii-Moriya interaction

$$H=\sum_i\mathbf S_i\cdot\mathbf S_{i+1}+D\sum_i\left(S^x_iS^y_{i+1}-S^y_iS^x_{i+1}\right)$$

```python
sc = spinchain.Spin_Chain(["S=1/2"]*n)
D = 0.3
h = sum(sc.SS(i, i+1) + D*(sc.Sx[i]*sc.Sy[i+1] - sc.Sy[i]*sc.Sx[i+1]) for i in range(n-1))
```

### Spin-1 chain with anisotropy and biquadratic exchange

$$H=\sum_i\left[\mathbf S_i\cdot\mathbf S_{i+1}+\beta\,(\mathbf S_i\cdot\mathbf S_{i+1})^2\right]+D\sum_i(S^z_i)^2$$

```python
s1 = spinchain.Spin_Chain(["S=1"]*n)
D, beta = 0.2, -0.3
h = sum(s1.SS(i, i+1) + beta*s1.SS(i, i+1)*s1.SS(i, i+1) for i in range(n-1))
h = h + D*sum(s1.Sz[i]*s1.Sz[i] for i in range(n))
```

$\beta=D=0$ is the Haldane chain, $\beta=-1/3$ the AKLT point, and a
large $D>0$ gives the trivial large-$D$ phase.

### Hubbard chain

$$H=-t\sum_{i\sigma}\left(c^\dagger_{i\sigma}c_{i+1\sigma}+\mathrm{h.c.}\right)
+U\sum_i\left(n_{i\uparrow}-\tfrac12\right)\left(n_{i\downarrow}-\tfrac12\right)$$

```python
fc = fermionchain.Spinful_Fermionic_Chain(n)
t, U = 1.0, 4.0
h = -t*sum(fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1] for i in range(n-1))
h = h + h.get_dagger()
h = h + U*sum((fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5) for i in range(n))
```

The shifted interaction puts the chemical potential at zero, so the
ground state is at half filling. With $U\sum_in_{i\uparrow}n_{i\downarrow}$
instead, fix the filling with a [sector](ground-state.md#fixed-particle-number-or-fixed-magnetization).

### Interacting spinless fermions

$$H=-t\sum_i\left(c^\dagger_ic_{i+1}+\mathrm{h.c.}\right)+V\sum_i\left(n_i-\tfrac12\right)\left(n_{i+1}-\tfrac12\right)$$

```python
fc = fermionchain.Fermionic_Chain(n)
t, V = 1.0, 2.0
h = -t*sum(fc.Cdag[i]*fc.C[i+1] for i in range(n-1))
h = h + h.get_dagger()
h = h + V*sum((fc.N[i]-0.5)*(fc.N[i+1]-0.5) for i in range(n-1))
```

This is the XXZ chain in fermionic language, with $t=J/2$ and
$V=J\Delta$. For $V>2t$ the chain is a charge density wave.

### Kitaev chain

$$H=\sum_i\left(-t\,c^\dagger_ic_{i+1}+\Delta\,c_ic_{i+1}+\mathrm{h.c.}\right)-\mu\sum_in_i$$

```python
fc = fermionchain.Fermionic_Chain(n)
t, Delta, mu = 1.0, 0.5, 0.3
h = sum(-t*fc.Cdag[i]*fc.C[i+1] + Delta*fc.C[i]*fc.C[i+1] for i in range(n-1))
h = h + h.get_dagger()
h = h - mu*sum(fc.N)
```

For $|\mu|<2t$ the open chain has Majorana end modes and its two lowest
states, of opposite fermion parity, are split only by $e^{-n/\xi}$. DMRG
then returns a superposition of the two, and the
[parity](topology.md#fermion-parity) of the state is not $\pm1$ (0.07
on 8 sites in our run). Take the two lowest states from
`get_excited_states(n=2)` and measure the parity of each, or use
`mode="ED"` on a short chain.

### SSH chain

$$H=\sum_i t_i\left(c^\dagger_ic_{i+1}+\mathrm{h.c.}\right),\qquad t_i=t_1,t_2,t_1,t_2,\dots$$

```python
fc = fermionchain.Fermionic_Chain(n)
t1, t2 = 1.0, 0.5
h = sum((t1 if i % 2 == 0 else t2)*fc.Cdag[i]*fc.C[i+1] for i in range(n-1))
h = h + h.get_dagger()
```

The chain is topological, with a zero mode at each end, when the weak
bond is at the edge; see the Zak phase in [Topological invariants](topology.md).

### Bose-Hubbard chain

$$H=-t\sum_i\left(a^\dagger_ia_{i+1}+\mathrm{h.c.}\right)+\frac U2\sum_in_i(n_i-1)-\mu\sum_in_i$$

```python
bc = bosonchain.Bosonic_Chain(n, maxnb=[4]*n)    # at most 3 bosons per site
t, U, mu = 1.0, 2.0, 1.0
h = -t*sum(bc.Adag[i]*bc.A[i+1] for i in range(n-1))
h = h + h.get_dagger()
h = h + 0.5*U*sum(bc.N[i]*(bc.N[i]-1) for i in range(n)) - mu*sum(bc.N)
```

`maxnb` cuts the local space; check that $\langle n_i\rangle$ and the
energy do not change when it is increased.

### $\mathbb Z_3$ clock model

$$H=-J\sum_i\left(\sigma^\dagger_i\sigma_{i+1}+\mathrm{h.c.}\right)-f\sum_i\left(\tau_i+\tau^\dagger_i\right)$$

```python
pc = parafermionchain.Parafermionic_Chain(n, Z=3)
J, f = 1.0, 0.7
h = -J*sum(pc.Sigd[i]*pc.Sig[i+1] for i in range(n-1)) - f*sum(pc.Tau)
h = h + h.get_dagger()
```

`Sig`, `Sigd`, `Tau`, `Taud` are the clock and shift operators and their
adjoints. The model is self-dual at $f=J$, where it is critical.

### Kondo lattice

$$H=-t\sum_{i\sigma}\left(c^\dagger_{i\sigma}c_{i+1\sigma}+\mathrm{h.c.}\right)
+J_K\sum_i\mathbf S_i\cdot\mathbf s_i$$

The local moments are spin sites between the fermion sites of a mixed
chain, each coupled to the two conduction sites next to it:

```python
mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F", "1/2", "F"])   # c S c S c
t, JK = 1.0, 0.8
h = t*(mc.Cdagup[0]*mc.Cup[2] + mc.Cdagdn[0]*mc.Cdn[2]
       + mc.Cdagup[2]*mc.Cup[4] + mc.Cdagdn[2]*mc.Cdn[4])
h = h + h.get_dagger()
h = h + JK*(mc.SS(0, 1) + mc.SS(1, 2) + mc.SS(2, 3) + mc.SS(3, 4))
```

`mc.SS(i, j)` is $\mathbf S_i\cdot\mathbf S_j$ with the fermion spin on a
fermion site and the local moment on a spin site. See
`examples/fermion_models/mixed_spin_fermion_chain`.

## Recipes

### Central charge of the critical Ising chain

For
$H=\sum_iS_i^zS_{i+1}^z+\tfrac12\sum_iS_i^x$, at the critical field, the
entanglement entropy of the open chain follows the Calabrese-Cardy formula
with $c=1/2$ ([Entanglement](entanglement.md)):

```python
sc.maxm = 200          # critical chains need a large bond dimension
wf = sc.get_gs()
c = wf.get_CFT_central_charge()
```

### Haldane gap of the spin-1 chain

The open chain has four low states
from the two spin-1/2 edge states. The bulk gap is the fifth state
([Excited states](excited-states.md)):

```python
s1 = spinchain.Spin_Chain(["S=1"]*n)
s1.set_hamiltonian(sum(s1.SS(i, i+1) for i in range(n-1)))
es = s1.get_excited(n=6)
gap = es[4] - es[0]
```

### Gaps of a Hubbard chain

Here `fc` is the Hubbard chain above. `get_gap()` is the gap to the first
excited state, and `get_charge_gap(d=...)` is $E_0(N+d)-E_0(N)$
([Excited states](excited-states.md)):

```python
fc.get_gap()
fc.get_charge_gap(d=1)
fc.get_charge_gap(d=2)
```

### Fidelity susceptibility across the Ising transition

See [Fidelity susceptibility](fidelity.md).

```python
import numpy as np
from dmrgpy import fidelity
h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))   # Ising coupling
h1 = 2*sum(sc.Sx[i] for i in range(n))                # transverse field
lambdas = np.linspace(0.5, 1.5, 11)
chis = [fidelity.get_fidelity(sc, h0, h1, lam) for lam in lambdas]
```

### Dynamical structure factor $S(q,\omega)$

See [Dynamical correlators](dynamical-correlators.md). We compute it from the
site-resolved correlators and a Fourier transform,
$S(q,\omega)=\frac1n\sum_{ij}e^{iq(i-j)}S_{ij}(\omega)$:

```python
es = np.linspace(-0.5, 3, 200)
q = np.pi
Sqw = 0
for i in range(n):
    for j in range(n):
        x, y = sc.get_dynamical_correlator(submode="KPM", es=es, delta=0.1,
                                           name=(sc.Sz[i], sc.Sz[j]))
        Sqw = Sqw + np.exp(1j*q*(i-j))*y/n
```

### Extrapolation to the thermodynamic limit

See [Post-processing tools](post-processing.md). Here we extrapolate the
energy per site computed on several lengths `Ls`:

```python
from dmrgpy import extrapolate
e_inf = extrapolate.size_extrapolation(Ls, e_per_site, n=1)
```
