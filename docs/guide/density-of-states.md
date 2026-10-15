# Density of states

"Density of states" can mean two different quantities, and only the second
one has a method.

## Many-body density of states

This is the density of all the
eigenstates of $H$,

$$\rho(E)=\sum_n\delta(E-E_n).$$

There is no method for it. `get_distribution` ([Density matrices and distributions](density-matrices.md)) is not a
replacement: it computes a ground-state expectation value
$\langle\mathrm{GS}|B\,\delta(X)\,A|\mathrm{GS}\rangle$, not a trace over
all the states. For a small chain, the full ED spectrum from
`get_excited(mode="ED", n=...)` ([Excited states](excited-states.md)) gives $\rho(E)$ directly.

## Local density of states

For fermions, the density of states measured
by STM ($dI/dV$) is the local single-particle spectral function

$$A_i(\omega)=\sum_n|\langle n^{N+1}|c_i^\dagger|\mathrm{GS}\rangle|^2\delta(\omega-\omega_n^+)
+\sum_m|\langle m^{N-1}|c_i|\mathrm{GS}\rangle|^2\delta(\omega-\omega_m^-),$$

with electron addition at $\omega>0$ and electron removal at $\omega<0$.
This is $A_{ii}(\omega)$ of [Dynamical correlators](dynamical-correlators.md), and we compute it with
`get_spectral_function(i=i)`, or from the two correlators
$(c_i,c_i^\dagger)$ and $(c_i^\dagger,c_i)$ as shown there. On a spinful
chain, pass `spin="up"` or `spin="dn"` and add the two.

![Local density of states of a fermion chain](img/density-of-states_ldos.png#only-light)
![Local density of states of a fermion chain](img/density-of-states_ldos-dark.png#only-dark)

The figure shows $A_i(\omega)$ of an open chain of 20 spinless fermions at
half filling, with
$H=-t\sum_i(c^\dagger_ic_{i+1}+\mathrm{h.c.})+V\sum_i(n_i-\frac12)(n_{i+1}-\frac12)$,
from `get_spectral_function` with `nex=50` and $\delta=0.1$. The left panel
is a site in the middle of the chain and the right panel the edge site. At
$V=0$ the result lies on the exact free-fermion spectrum (dashed): 20
levels inside the band $|\omega|<2t$, with a weight that depends on the
site. In the middle of the chain the weight piles up at the band edges, as
in $1/\sqrt{4t^2-\omega^2}$. At the edge it is largest at the band center.
At $V=3t$ the chain is a charge density wave, and a gap of about $t$ opens
at $\omega=\mu$. The weight moves to higher energies, around $\pm V$. Here
the 50 states hold 79 per cent of the weight at the middle site, so part of
the spectrum at high energy is missing; `info["captured"]` tells this.

??? example "Code of this figure (docs/figures/density-of-states_ldos/)"

    ```python title="compute.py"
    --8<-- "density-of-states_ldos/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "density-of-states_ldos/plot.py"
    ```
