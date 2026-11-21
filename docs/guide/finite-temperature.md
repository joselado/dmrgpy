# Finite temperature

On this page we compute thermal averages,

$$\langle O\rangle_T=\frac{\mathrm{Tr}\,[e^{-\beta H}O]}{\mathrm{Tr}\,e^{-\beta H}},
\qquad\beta=1/T.$$

There are three methods: purification, METTS and exact diagonalization.
METTS is built on imaginary-time evolution, which has its own section.

## Purification

Each physical site gets an ancilla site. The product of
singlets between each site and its ancilla is a pure state whose reduced
density matrix on the physical chain is $\rho\propto\mathbb 1$, the
$T=\infty$ state. Evolving it in imaginary time with the physical $H$,

$$|\Psi(\beta)\rangle=e^{-\beta H/2}|\Psi(0)\rangle,$$

and tracing out the ancillas gives $\rho(T)\propto e^{-\beta H}$.

```python
from dmrgpy import thermal
n = 6
tc = thermal.Thermal_Spin_Chain(["S=1/2"]*n, T=0.5)
h = sum(tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
        for i in range(n-1))
tc.set_hamiltonian(h)
tc.get_gs()                                   # builds |Psi(beta)>
E = tc.MBChain.vev(h).real                    # thermal energy
czz = tc.MBChain.vev(tc.Sz[0]*tc.Sz[1]).real  # thermal <Sz_0 Sz_1>
```

`tc.Sx`, `tc.Sy`, `tc.Sz` are the operators of the physical sites, and
`tc.MBChain` is the doubled chain that holds the purified state. Static
expectation values on `tc.MBChain` are thermal averages. A dynamical
correlator on `tc.MBChain` is not the thermal one, so use METTS or ED for
those. `T=0` gives the ordinary ground state. After a change of `tc.T` or
of a keyword below, call `tc.get_gs()` again before reading `tc.MBChain`,
which then holds the new state. The imaginary-time evolution
is accurate to about $10^{-6}$ relative at the defaults. The keywords are:

- `tc.anneal_step`: the size of each imaginary-time step, in units of the
  spectral width (default 2). Smaller is more accurate and slower.
- `tc.anneal_order`: the order of the expansion of each step (default 8).
- `tc.anneal_method`: how $e^{-\beta H/2}$ is applied. `"taylor"` is the
  expansion above, `"implicit"` the implicit steps of
  [Imaginary-time evolution](#imaginary-time-evolution), whose number grows
  as $\log\beta$ and not as $\beta W$, with $W$ the spectral width. The
  default, `None`, takes the implicit steps on `itensor_version=3` at
  $T\le0.05$ (in units of the largest coupling) when they need far fewer
  solves than the expansion needs steps, that is on long chains, and the
  expansion otherwise. Setting `anneal_step` or
  `anneal_order` keeps the expansion. The implicit steps run on
  `itensor_version=3`, `"python"` and `"julia_live"`.
- `tc.anneal_options`: the keywords of `evolve_imaginary_time` for
  `"implicit"`, except `shift`. The default is `order=8` with `nsteps=2`:
  the $T=\infty$ state holds every energy of the spectrum, and the
  trapezoid rule is not accurate enough for it.
- `tc.MBChain.maxm`: the bond dimension of the purified state. It needs
  to be larger at intermediate temperatures.

See `examples/finite_temperature/thermal_purification_VS_exact`.

![Thermal energy and specific heat of the Heisenberg chain](img/finite-temperature_specific_heat.png#only-light)
![Thermal energy and specific heat of the Heisenberg chain](img/finite-temperature_specific_heat-dark.png#only-dark)

The figure shows the thermal energy and the specific heat per site of a
spin-1/2 Heisenberg chain with $n=10$. The specific heat is the energy
variance, $C(T)=(\langle H^2\rangle_T-\langle H\rangle_T^2)/T^2$, measured
in the purified state, so no numerical derivative is needed. The reference
is the full spectrum of the chain (dashed lines). Purification agrees with
it to $10^{-5}$ at every temperature. $C(T)$ has a broad maximum near
$T\approx0.46J$, the scale of the antiferromagnetic correlations. The
METTS points are discussed in [METTS](#metts).

??? example "Code of this figure (docs/figures/finite-temperature_specific_heat/)"

    ```python title="compute.py"
    --8<-- "finite-temperature_specific_heat/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "finite-temperature_specific_heat/plot.py"
    ```

## Imaginary-time evolution

Imaginary-time evolution projects a state onto the low energies,

$$|\psi(\tau)\rangle\propto e^{-\tau H}|\psi\rangle
=\sum_n e^{-E_n\tau}|n\rangle\langle n|\psi\rangle .$$

At long $\tau$ it gives the ground state, and at $\tau=\beta/2$ it is the
step of METTS below. Each eigenstate decays as $e^{-(E_n-E_0)\tau}$, so the
energy scale that matters falls as $\tau$ grows, and the time steps can grow
with it (Zima, Stoudenmire, White, Parcollet and Kaye, arXiv:2606.02930).
The grid has panels $[0,\tau_0]$, $[\tau_0,2\tau_0]$, $[2\tau_0,4\tau_0]$,
... up to $\tau$, with a few equal steps in each, so the number of steps
grows as $\log\tau$ and not as $\tau$. A long step needs an implicit
rule, the trapezoid rule

$$\Big(1+\tfrac{\Delta\tau}{2}(H-E_0)\Big)|\psi_k\rangle=
\Big(1-\tfrac{\Delta\tau}{2}(H-E_0)\Big)|\psi_{k-1}\rangle ,$$

which is stable for any step size. Each step needs one linear solve. Here
`h` is the Hamiltonian of the chain and
`wf` the initial state, for example a Néel state:

```python
wft = sc.evolve_imaginary_time(wf, 100.)     # normalized e^{-100 H}|wf>
taus, e = sc.evolve_imaginary_time_and_measure(h, wf, 100.)
```

The second call returns the grid times and $\langle\psi(\tau)|h|\psi(\tau)\rangle$
on them. With an ED state, `mode="ED"` gives the exact evolution on the
same grid. In `examples/finite_temperature/imaginary_time_log_grid` the
Néel state of a 12-site Heisenberg chain reaches $\tau=1000$ in 72 steps,
and its energy ends $5\times10^{-12}$ above $E_0$. The error of the trapezoid
rule falls as $1/n^2$ with the number $n$ of steps per panel. It is largest
at intermediate $\tau$, where the energy is off by about $2\times10^{-3}$
(in units of the exchange) at $n=4$ in that example. The keywords are:

- `nsteps`: the steps per panel (default 4).
- `tau0`: the end of the first panel (default 0.01 over the largest
  coupling of $H$). The first steps must resolve the highest energies the
  initial state contains.
- `order`: 2 for the trapezoid rule, 4, 6, ... for the Gauss-Legendre
  rules, more accurate per step, with `order/2` complex linear solves each.
- `shift`: the energy subtracted from $H$, by default the ground-state
  energy. It must be the ground-state energy: with a lower shift the long
  steps stop projecting.

It works with `itensor_version=3` (at least 3 sites), `"python"` and
`"julia_live"`, for a Hermitian $H$.

## METTS

Minimally entangled typical thermal states (Stoudenmire and
White, arXiv:1002.1305) sample product states $|i\rangle$ with a Markov
chain. Each one is evolved to
$|\phi_i\rangle\propto e^{-\beta H/2}|i\rangle$ and then collapsed to a new
product state. The plain average of $\langle\phi_i|O|\phi_i\rangle$ over
the chain is the thermal average. The states stay weakly entangled, so
no ancillas are needed. Here `sc` is the Heisenberg chain of [Ground state](ground-state.md):

```python
T = 0.5
mean, err = sc.metts_vev(sc.Sz[0]*sc.Sz[1], T, nsamples=300, nwarmup=30)
res = sc.metts_vev([sc.Sz[0], sc.Sz[1], h], T)   # list of (mean, err)
```

Being a Monte Carlo method, it returns a mean and a statistical error.
The METTS energies in the figure of [Purification](#purification), with
300 samples per temperature, agree with the exact energy within two error
bars.
The error assumes independent samples, so it is too small when samples
are correlated. Measuring several operators in one call is much cheaper
than one call per operator, since they share the samples. The keywords
are:

- `nsamples`, `nwarmup`: the number of samples, and the number of first
  steps discarded for equilibration.
- `dbeta_half_step`: the imaginary-time step (default 0.05).
- `basis_ops`: the two bases used in turn for the collapse (default
  `("Sz","Sx")`). Alternating them reduces autocorrelation.
- `njobs`: runs independent Markov chains in parallel processes. Only
  for `itensor_version="python"`.
- `itime`: how each sample is evolved to $\beta/2$. By default the cheaper
  of two ways at the same accuracy: TDVP steps of `dbeta_half_step` at
  high temperature, and below $T\approx0.12$ (in units of the largest
  coupling) the implicit steps of
  [Imaginary-time evolution](#imaginary-time-evolution), whose number
  grows as $\log\beta$ instead of $\beta$. `itime="tdvp"` or
  `itime="implicit"` forces one, and `itime_options` takes the keywords of
  `evolve_imaginary_time`, with 8 steps per panel by default here.

METTS works with `itensor_version=3` (at least 3 sites), `"python"` and
`"julia_live"`. It
is not available on `itensor_version=2`, where it raises an error. It
also raises if the samples cannot leave one sector of a charge or a parity
that $H$ conserves. This is the case when every basis operator commutes
with a conserved particle number or $S^z_{\rm tot}$. It is also
the case for a single $S^x$ when $H$ conserves $\sum_i S^x_i$ or
$\prod_i 2S^x_i$. On a fermion chain it holds for every basis, since $H$
conserves the fermion parity and no Hermitian operator of one fermion site
changes it. So METTS is for spin chains, and a fermion model is sampled as
a spin chain through the Jordan-Wigner transformation. The default
`basis_ops` avoids all of this on a spin chain.

The same sampling gives real-time thermal correlators
$C_{AB}(t)=\langle A(t)B\rangle_T$ (arXiv:2405.18484). Each sample
$|\psi_i\rangle$ gives $\langle\psi_i|e^{iHt}Ae^{-iHt}B|\psi_i\rangle$ by
real-time evolution, and the average is $C_{AB}(t)$:

```python
ts, c, err = sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), T,
                 nt=100, dt=0.1, nsamples=200, nwarmup=30)
```

The result is in the time domain, on $t=0,\Delta t,\dots,(n_t-1)\Delta t$.
For a spectrum, apply a window and a Fourier transform yourself.

## Exact diagonalization

For small chains we sum over the eigenstates,
$\langle O\rangle_T=\sum_n p_n\langle n|O|n\rangle$ with
$p_n=e^{-E_n/T}/Z$. Here `sc` is a short chain, `n` of about 10 sites:

```python
sc.vev(sc.Sz[0]*sc.Sz[1], mode="ED", T=0.5)
sc.get_correlation_matrix(T=0.5, operators=[sc.Sz[i] for i in range(n)])
sc.get_dynamical_correlator(mode="ED", submode="ED", T=0.5,
                            name=(sc.Sz[0], sc.Sz[0]))
```

The default `vev` on a DMRG chain raises at `T>0`, so `mode="ED"` is
needed. The sum uses the full spectrum when the Hilbert space has at most about
2000 states. Above that, only the `n=` lowest states are used, and an
error is raised if the last one still has a large Boltzmann weight. The
ED dynamical correlator is the exact reference for
`metts_dynamical_correlator` (see
`examples/finite_temperature/dynamical_metts_VS_ED`).
