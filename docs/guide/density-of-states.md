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
