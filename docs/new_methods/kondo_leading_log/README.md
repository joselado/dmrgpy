# Leading-log resummation of the STM Kondo spectrum

This folder documents a method that is new in dmrgpy and not taken from
a single paper: the leading-log resummation of the third-order STM
spectrum of Ternes (NJP 17 063016, arXiv:1505.04430) for an impurity
with an arbitrary spectrum (anisotropy, field, coupled spins). It holds
the derivation, the references it rests on, the checks it passes and
what it does not contain. The code is
`src/dmrgpy/kondospectrumtk/leadinglog.py`, the tests are
`tests/test_kondo_leading_log.py`.

**Status:** see the section "Status".

## 1. Why a resummation and not a fourth order

The third order adds to the golden rule one exchange scattering with a
sample electron, and this gives the logs F(eV - eps) at each spin
excitation. Each further order adds one more power of the same log. So
the expansion parameter is rho J ln(omega0/|x|), not rho J, and it
reaches 1 at the Kondo scale T_K. A bare fourth order adds nothing
where the log is large, which is exactly where the third order fails.
What helps is to sum the leading logs, the terms (rho J)^n ln^n, to all
orders. This is what the perturbative renormalization group does, and
in a field or with anisotropy it needs couplings that depend on the
energy of the scattered electron (Rosch, Paaske, Kroha, Woelfle,
cond-mat/0202404 and cond-mat/0408506).

## 2. The model and the vertices

The model is Ternes' one. An impurity with eigenstates |a> and energies
e_a (the `KondoSpectrum` of `edkondo.py`, or any object with the same
fields) couples to the sample electrons at its site by the exchange
J psi^dag (sigma/2) psi . S, and the tip tunnels into the same site
through the vertex T0 (sigma.S/2 + U).

We write every vertex as an operator on the electron spin times an
operator on the impurity,

    V = sum_mu sigma^mu (x) A^mu,   mu = 0, x, y, z,   sigma^0 = 1,

so A^0 is potential scattering and A^x, A^y, A^z are exchange. The
product of two electron-spin matrices closes on the same set,
sigma^mu sigma^nu = sum_lambda C^{mu nu}_lambda sigma^lambda, with
C^{00}_0 = 1, C^{0a}_a = C^{a0}_a = 1, C^{ab}_0 = delta_ab and
C^{ab}_c = i eps_abc. So a vertex is four impurity matrices in the
eigenbasis.

There are two vertices:

- the sample vertex G, with bare value A^a = g S^a / 2 and A^0 = 0,
  where g is fixed by `Jrho_s` (section 5 gives the sign);
- the tip vertex A(omega), with bare value A^a = S^a / 2 and A^0 = U,
  where omega is the energy of the incoming tip electron, so omega = eV.

## 3. The flow

Let us take the band cutoff D down step by step. A second-order process
with an intermediate impurity state c and an intermediate band electron
or hole at the cutoff gives a correction to the vertex. There are two
such processes, as in the third order:

- the direct (particle) one, where the incoming electron scatters first
  and the intermediate electron is in an empty state. Its log is cut
  when D drops below |omega - (e_c - e_b)|, the distance of the
  intermediate state from its own energy shell;
- the exchange (hole) one, where a sample electron below the Fermi
  level scatters first into the outgoing state and the incoming
  electron fills the hole. Its log is cut at |omega - (e_a - e_c)|, and
  the electron-spin matrices appear in the reversed order, with a minus
  sign from the fermion exchange.

These are the two arguments of `conductance.py`'s third-order terms,
F(eV - eps_im) and F(eV - (e_f - e_m)) (the second departs from the
paper's eq. 25; see the docstring of `conductance.py` and
`.claude/rules/kondo.md`).

The approximation that makes this closed is the one of Rosch et al.
(cond-mat/0408506, eq. 6 and the text after it). The couplings on the
right-hand side are taken on shell, at the energies where the
intermediate electron sits at the Fermi level, so that the frequency
dependence enters only through the cut. This means two sets of
matrices:

- the sample vertex G, with its incoming electron at the Fermi level;
- the tip vertex element (a,b) at its own on-shell energy e_a - e_b,
  called A_on.

For the tip vertex at a bias omega this gives

    dA^lambda_ab(omega) = - sum_c sum_{mu nu} [
          C^{mu nu}_lambda G^mu_ac A_on^nu_cb K(omega - (e_c - e_b))
        - C^{nu mu}_lambda A_on^mu_ac G^nu_cb K(omega - (e_a - e_c)) ] dD,

so A(omega) is a quadrature over D once G and A_on are known. These two
flow together and close: G with the same equation at omega = 0 (cuts
|e_c - e_b| direct and |e_a - e_c| exchange), A_on with the cuts read at
its own on-shell energy (|e_a - e_c| direct and |e_c - e_b| exchange).
The overall sign and the normalization of g are fixed in section 5.
For S=1/2 in a field these are eqs. 14 and 15 of cond-mat/0408506 and
eq. 6 of cond-mat/0301106, at every omega, and the second of these
papers checks them against bare perturbation theory to order
g^3 ln^2.

Two remarks on the closure. First, G is not symmetrized: G_ab and G_ba
have different outgoing energies and different cuts, and with
anisotropy or a tilted field the flow makes G non-Hermitian and also
generates a small potential scattering. Then the tip vertex at
omega = e_f - e_i is the sample vertex G_if^*/g only up to a term
beyond leading log (5.0e-3 relative at g=-0.05 and 1.2e-3 at g=-0.025
on the S=1 model of section 5; at fixed g it grows only 2 to 4 times
while ln(omega0/Gamma0) grows by 9, much slower than a power of the
log). For S=1/2 in a field G stays
Hermitian and the identity is exact. This is a choice beyond leading
log, between G at incoming energy zero and G on shell on the
right-hand side, and the code takes the first one. Second, the first
version of this code (commit 5f845c0) kept the running tip vertex
A(omega) on the right-hand side instead of A_on. It agreed with this
one at g^3 ln^2, and it differed by a term of order g^2 with no log,
1 to 3% of the vertex at g=0.04 for S=1/2 in a field; an independent
review found it, and the code now follows the papers.

**The kernel.** The sharp cut Theta(D - |x|) is replaced by the smooth
kernel

    K(x, D) = 1/(D + r_x) - 1/(D + omega0 + |x|),   r_x = sqrt(x^2 + Gamma0^2),

whose integral over D from 0 to infinity is exactly the F0(x) of
`stepfunctions.py`. So to first order in g the flow gives back the
third order of `conductance.py` exactly, with the same omega0 and
Gamma0, and the two methods differ only from order g^2 on. Gamma0 cuts
the logs, as the decoherence rate Gamma does in cond-mat/0202404.

**The spectrum.** The renormalized tip vertex goes into the golden rule:

    dI/dV(eV) = 2 pi T0^2 sum_{i,f} p_i * 4 sum_lambda |A^lambda_fi(eV)|^2
                * Theta(eV - (e_f - e_i))  +  (s->t),

which is the second-order formula of `conductance.py` with the bare
vertex replaced by the flowed one, as in cond-mat/0301106, eq. 11. The
s->t direction is the particle-hole conjugate of the t->s one: eV goes
to -eV and U to -U, while the exchange is unchanged. This gives
`conductance.py`'s rule that the Kondo term is even in eV and the
potential term is odd.

## 4. What the method does not contain

- It is perturbative. It fails at T_K by construction, where the
  couplings reach order 1. Two checks report this. The sample flow
  stops when the sample coupling passes `pole_threshold` (1 by
  default); then the biases within 10 D_stop of any energy difference
  of the spectrum are NaN. And at each bias the running tip coupling,
  |Jrho_s| times the largest ratio of the flowed to the bare spin
  vertex, is returned as `info["coupling"]`; above 0.5 the code warns,
  at `pole_threshold` the bias is NaN. The ratio is taken over the
  transitions out of the occupied states that are open at this bias,
  the ones that enter dI/dV there, and a bare vertex below 0.1 of the
  largest one is replaced by 0.1 of the largest one (`BARE_FLOOR`);
  the end of section 7 gives the reasons and the measurements. The
  second check matters near a threshold, where the tip vertex can grow
  while the sample coupling is still small: S=3/2 with D=1 meV,
  B=0.3 meV and Jrho_s=-0.1 has a running coupling of 0.89 at eV=B and
  dI/dV/2 pi of 59 there, against 5.5 from the third order. Neither
  number is reliable there, and the warning says so.
- One loop renormalizes couplings, not energies. The steps stay at the
  bare excitation energies, so the renormalized anisotropy and g-factor
  (Delgado, Hirjibehedin, Fernandez-Rossier, arXiv:1401.7272; Oberg et
  al., arXiv:1407.1542) are not in this method.
- Gamma0 is a fixed cut. The Korringa decay of the excited states,
  which broadens the steps by about (rho J)^2 times the splitting, is
  not yet included.
- The populations are the equilibrium ones of the impurity, so there is
  no spin pumping by the tip.
- At T>0 the cutoff kernel is its average over the thermal window of the
  third order (section 8). Which kernel the on-shell flow takes, and
  whether the amplitude or the rate is averaged, are choices beyond first
  order; section 8 says which ones are taken.

## 5. Validation

Each check is in `tests/test_kondo_leading_log.py`.

1. S=1/2 at B=0, any U. Every process is elastic, and the flow has a
   closed form: the exchange part of the tip vertex goes as
   1/(1 - g F0(eV)). This checks the whole pipeline (kernel, tails,
   integration) against a formula: with this kernel the tip vertex is
   t(eV) = 1 - int_0^inf g(D)^2/g K(eV, D) dD, with
   g(D) = g/(1 + g ln((D + omega0)/(D + Gamma0))), which is
   1/(1 - g ln(omega0/|eV|)) with a sharp cut. Measured: 1.5e-8
   relative at g=-0.05, U=0.25, on six biases from -3 to 9 meV.
2. S=1/2 in a field. The flowed sample vertex against an independent
   integration of eq. 6 of cond-mat/0301106 with the same kernel. This
   is the check of the multi-level cut structure. Measured: 1e-8 or
   better in both couplings along the whole flow, at g=-0.08, B=1 meV,
   and no potential scattering is generated at omega=0, as their g_z(0)
   does not depend on the electron spin. The tip vertex at
   omega = 0.3, 0.7, -0.99 and 2.5 meV against the same equations at
   that omega (g_z+, g_z-, g_perp(omega -+ B/2)): below 1e-8. And the
   on-shell identity of the remark above, exact here, below 1e-8.
3. First order against `conductance.py` on the models of
   `tests/test_kondo_spectrum_tmatrix.py` (anisotropic S=1, coupled
   spins, tilted fields, U on). The difference from the second plus
   third order must scale as g^2. This also fixes the sign of g against
   `Jrho_s`, and we check that `Jrho_s < 0` makes the sample coupling
   grow as D decreases, which is the Kondo sign in the paper's
   convention. Measured, with U=0.25: when g is halved from 0.01 the
   difference drops by 4.10 to 4.13 on the three models, and its odd
   part in g by 8.02, so the first order is exactly the third order.
   Jrho_s<0 is the sign that grows, as it should.
4. The stop at the Kondo scale: S=1/2 at B=0 and Jrho_s=-0.2, where
   T_K is above Gamma0. The flow stops, the biases near zero are NaN
   with a warning, the far ones are finite. And the S=3/2 case of
   section 4, where only the tip coupling is large, warns.

An independent review on 2026-11-06 rederived the flow from eq. 6 of
cond-mat/0408506, checked the ln^2 coefficients of the g^3 terms at
nonzero omega against eq. 8 of cond-mat/0301106 (agreement to three
digits on four geometries), and checked invariance under a change of
basis in a degenerate ground doublet (7e-11) and under a rotation of
the spin axes (5e-11). Its findings (the closure above, the tip coupling
check, the NaN mask over all energy differences, `Gamma0 > 0`) are
fixed.

The tests run in a few seconds. They were run on d22 on 2026-11-06.

## 6. What the resummation changes

At the paper's own coupling, Jrho_s=-0.05 with omega0=20 meV and
Gamma0=5 microeV, a free S=1/2 at B=0 with U=0.25 gives dI/dV/2 pi,
in units of e^2 T0^2/h, of 2.19 at eV=10 microeV against 1.56 from the
second plus third order, and 1.44 against 1.29 at 0.4 meV. The
difference is the g^2 ln^2 and higher terms, and it is large already
here because g ln(omega0/|eV|) is about 0.38 at 10 microeV. Far from
zero bias the two agree: 1.10 against 1.09 at 9 meV.

## 7. How it is computed and what it costs

The flow of G and A_on, and the quadrature that gives the tip vertex at
every bias, are one system of equations in ln D, solved once with an
adaptive Runge-Kutta method (`leadinglog._joint_flow`). Nothing is
stored along the flow. Three facts set the cost for an impurity with d
states.

- dI/dV reads the tip vertex only on the columns of the occupied
  states, A_fi with i occupied, and the quadrature of a column does not
  need the other columns. So the tip vertex costs d^2 per bias and per
  occupied state, where the full matrix costs d^3 per bias.
- A vertex sum_mu sigma^mu (x) A^mu is a 2d x 2d matrix, and the
  contraction with C^{mu nu}_lambda is its matrix product. The hole
  channel is the same product with the electron-spin index transposed.
  So the flow of G and A_on needs four 2d x 2d products per step, and
  this is the only d^3 work that is left.
- The potential scattering of the tip does not flow. For a vertex
  proportional to the identity the direct and the exchange term are
  the same matrix, G^lambda K(omega), and they cancel at every cutoff.
  So A(omega; U) = A(omega; 0) + U exactly, the flow is solved at U=0,
  and the two tunnelling directions share one solve.

When the sample coupling reaches `pole_threshold` at D_stop the
couplings are held, and what is left of the quadrature is the integral
of K over D from D_min to D_stop, a difference of two logs.

The first implementation solved the on-shell flow with dense output
and then one flow per bias for the full d x d tip vertex. The two give
the same numbers: 1e-11 relative on chains of 3 to 5 spins 1/2, the
same NaN biases and D_stop to 1e-12 on the stopped flows, and the
values of the first one are pinned in the tests.

Measured, on a chain of spins 1/2 with exchange 1 meV, a tilted field,
Jrho_s=-0.05 and U=0.25:

- 5 spins (32 states), 61 biases, one core: 110 s before, 1.6 s now.
- 6, 7 and 8 spins, 61 biases, one core: 7 s, 35 s and 176 s.
- 9 spins (512 states), 241 biases, on d22 (i5-12600K): 392 s with the
  default threads and 658 s on one core, 1.7 GB, 578 evaluations of the
  right-hand side. With `rtol=1e-12` the curve moves by 1.3e-11. The
  first implementation was not run at this size: it needs about 350 MB
  of dense output per step and one flow of d^3 work per bias.

What is left at 9 spins is, per evaluation of the right-hand side and
on one core, 0.59 s for the four products and 0.38 s for the exchange
sum of the tip vertex, which has one kernel for every bias, final state
and intermediate state (482 x 512 x 512). The products use the threads
of the BLAS library (0.12 s on six cores). The exchange sum runs on one
core. With numba threads it took 0.04 s on six cores; this is not in
the code, since it would not follow the BLAS thread settings.

The two kernels alone do not gain on a weak GPU. On the RTX A2000 of
d22, in float64, the four products take 0.235 s, slower than six cores,
and the exchange sum 0.060 s, slower than the threaded sum. In complex64
the products take 3 ms, but the bracket is then 3e-4 off, far above the
tolerance of the flow. What does gain is the whole flow on the card,
since then the exchange sum and the integrator do not run on one core:
see "The flow on a GPU" at the end of this section.

**Real arithmetic and a table for the exchange sum (2026-11-11).** For
a real Hamiltonian every vertex is a real 2d x 2d matrix in block form.
The block form is [[A0+A3, A1-iA2], [A1+iA2, A0-A3]], and it is real
when Sx and Sz are real and Sy is imaginary in the eigenbasis, which is
the case in a real eigenbasis. The kernel is real, and the products and
the electron transpose keep a real matrix real. So `_joint_flow` tests
once that the block form of the bare vertex is real to 1e-12 of its
largest entry (`_block_is_real`). If it is, the state is kept in block
form and the flow runs in real arithmetic: the four products are dgemm
instead of zgemm, and the exchange sum of the tip vertex has 4m columns
instead of 8m. If it is not, as with a field along y, the complex code
runs as before.

The exchange sum recomputes sqrt(x^2+Gamma0^2) and omega0+|x| for every
final state, bias and intermediate state at every evaluation, and none of
them depends on D. They are now tabulated once per flow, which takes
16 n d^2 bytes for n distinct incoming energies, and the sum is one
division per element. The table is used only below `HOLE_TABLE_BYTES`
(256 MB). Above it, and for the closed-form integral below D_stop, the
sum is computed on the fly as before. The 9-spin call with 241 biases
needs 2 GB for the table, so it is above the cap, and it is unchanged.
The table is a gain for the cheap sizes only.

Measured on d22, one core, through `get_kondo_spectrum(order="LL")`, on
the chain above with 61 biases (median of 3):

- 6 spins: 2.80 s before, 0.99 s now (2.8 times).
- 7 spins: 13.52 s before, 4.04 s now (3.3 times). The largest change of
  dI/dV is 2.7e-11 relative at 7 spins and 1.2e-12 at 6 spins.
- 6 spins with a field along y of 3e-5 eV, which takes the complex
  route: 2.83 s before, 2.41 s now. Only the table acts here, and dI/dV
  changes by 9e-16.

The change of 1e-11 on the real route is not roundoff. The integrator
sees the state in another layout, so its error norm picks other steps
(41 and 38 accepted steps at 7 spins, 614 and 518 evaluations), and the
result moves inside the tolerance rtol=1e-10. Against a reference at
rtol=1e-12 the new numbers are as close as the old ones. The table alone
changes dI/dV by roundoff.

A site with no spin on any state, such as an empty orbital in a conserved
sector, has no exchange vertex. The flow does nothing and dI/dV is the
potential scattering 2 pi 4 U^2 of the second order, with the tip
coupling at zero and no warning.

**The tip coupling check on a large impurity.** Until 2026-11-11 the
check took the ratio of the flowed to the bare vertex over every final
state with a bare vertex above 1e-8 of the largest one. On a chain this
made most biases NaN from 6 spins on, 211 of 241 at 9 spins, at
Jrho_s=-0.05. Two things were wrong with it.

- It read closed transitions. At 6 spins every NaN bias was set by a
  transition to a state above the bias, which does not enter dI/dV
  there. A closed transition can also have a large bare vertex: for the
  S=1 model of section 5 at Jrho_s=-0.15 and 13 microeV the transition
  to the first excited state has a vertex 7 times its bare one, and the
  bias was NaN, while the elastic transition, the only open one, is at
  1.9 times.
- It divided by a bare vertex that can be almost zero. The flow gives a
  weak or forbidden transition a vertex of order g ln times the large
  matrix elements. This is under control, and it is many times the bare
  value of that transition: at 9 spins a transition with a bare vertex
  9e-7 of the largest one had a flowed vertex 6.5e-4 of the largest
  one, a ratio of 700. With only the open transitions and the old cut
  at 1e-8, a chain of 4 spins 1 with a weak transverse field still had
  one NaN bias and six warnings from such transitions.

So the check now takes the open transitions only, and measures a weak
transition against 0.1 of the largest bare vertex. A transition with a
bare vertex above that is judged as before. A weak or forbidden one is
flagged only when its flowed vertex is comparable to the largest bare
vertex. Measured at Jrho_s=-0.05, on bias windows that cover the whole
spectrum: chains of 6 and 7 spins 1/2 (three fields, tip on the end and
in the middle, Ising-like exchange) and of 3 and 4 spins 1 have no NaN
and no warning, and a coupling between 0.062 and 0.065, the same to
0.004 for a floor of 0.03 or 0.3 except on the 4 spins 1 (0.085 and
0.061). At 9 spins all 241 biases are finite, where 211 were NaN,
with a coupling between 0.057 and 0.066 and no warning; the 30 that
were finite before have the same values, and the curve is between 0.997
and 1.095 times the second plus third order. The check still follows
the growth of the vertex on a chain: at 6 spins the largest coupling
goes 0.064, 0.123, 0.178, 0.252 for Jrho_s = -0.05, -0.08, -0.1, -0.12,
while dI/dV over the second plus third order goes up to 1.09, 1.26,
1.48, 1.83. The checks that should fire still do: 0.892 for the S=3/2
case of section 4, and the free S=1/2 at Jrho_s=-0.12 and -0.2 as
before. One case is now stricter. The ground state of a dimer in a
field has no spin, so its elastic vertex is zero and the old check did
not see it. Below the first step it is the only open transition. At
Jrho_s=-0.15 and 8 microeV from the step its flowed vertex is 1.1 times
the largest bare vertex, and this bias is now NaN (3 of 21 biases on
the grid of the tests, none before).

What stays as it was: when the sample flow stops at D_stop, every bias
within 10 D_stop of any energy difference of the spectrum is NaN, and
on a dense spectrum this is every bias. This happens only when the
Kondo scale is above Gamma0, where the method fails anyway.

**The flow on a GPU (2026-11-12).** With
`pyitensor.backend.set_backend("jax")` before the call, `_joint_flow`
hands the flow to `leadinglog_device.py`, and it runs on a GPU when JAX
has one. This is a port, not a new method. The right-hand side is the
same, written in jax.numpy, with the same layout of the state. The
integrator is SciPy's DOP853, with its coefficients, its initial step,
its error norm and its step control, but with the state on the device.
So every array of the size of the state stays on the card, and the host
keeps only the scalars of the integrator and the tables of the thermal
kernel.

The numbers are the NumPy ones to roundoff: below 1.1e-15 relative on
JAX's CPU device, on the real and the complex route, at T>0 with both
on-shell kernels, at a stop and with chunks
(`tests/test_kondo_leading_log_device.py`), and below 8e-15 on the
A2000 of d22 against one core, on chains of 4 to 8 spins. The step sizes
agree only to about 1e-7 of themselves. The error estimate of DOP853 is
a difference of nearly equal stages, and its last digits depend on the
order of the sum. This moves the curve far below the tolerance.

Measured on d22, on the chain of section 8 (open chain of spins 1/2,
exchange 1 meV, a field of 0.3 meV along z and 0.1 meV along x on every
site), Jrho_s=-0.05, U=0.25 and 61 biases. The six-core runs shared the
host with two other jobs.

| case | one core | six cores | RTX A2000 |
|---|---|---|---|
| 8 spins, T=0 | 26 s | 15 s | 12.6 s |
| 9 spins, T=0 | 160 s | 119 s | 41 s |
| 10 spins, T=0 | | 532 s | 283 s |
| 9 spins, field along y (complex route) | 475 s | | 116 s |
| 8 spins, T=2 K (223 occupied states) | 1179 s | | 164 s |

So even a card with weak double precision gains, by 2 to 4 against six
cores at 9 and 10 spins. The products alone are not faster there. What
it gains is the rest of the right-hand side and the integrator, which
run on one core in the NumPy route. The gain falls with size on this
card (2.9 at 9 spins, 1.9 at 10), since the products, where it is
slower than six cores, grow as d^3. A card with full double precision
(V100, A100, H100, H200) is the one for large impurities.

The memory on the device is `leadinglog_device.memory_model(d, m,
real)`. It counts about 20 copies of the state, against about 24 for
SciPy's DOP853. The couplings are 8 d^2 numbers, and the tip vertex is
4 d m per energy. It also counts about ten 2d x 2d matrices of the
right-hand side and three chunks of the tip part. It gives 6.1 GB for
the T=2 K run above, where 6.15 GB were measured, and it is the model
`mem_gb` uses on this route. For 13 spins 1/2 (d=8192) on the real route
it gives about 110 GB before the first bias, so only a 141 GB card holds
it. Each energy then adds 5.2 MB per occupied state. JAX caps its
allocator at 75% of the card by default, also with
`XLA_PYTHON_CLIENT_PREALLOCATE=false`, and on the H200 the cap is close
to 110 GB. So a 13-spin run needs `XLA_PYTHON_CLIENT_MEM_FRACTION=0.95`
in its environment. On the A2000 a 10 GB array failed under the default
and was allocated at 0.95. A run at 13 spins and 0.5 K on the chain
above keeps 24 occupied states, which needs about 124 GB for 61 biases
by the model, and only fits with the larger fraction.

The tip part has two forms (`TIP_FORM`). "columns" is `_tip_bracket`
line for line: the electron index is summed before the kernel is applied,
so each chunk of occupied columns has operands of 4 d^2 numbers per
column, and it builds the exchange kernel again. "block" applies the
kernel to the vertex and does one product per bias for all the occupied
columns. It has twice the arithmetic, but the kernel is built once per
bias. On the A2000 the column form was faster in every case: 41 s
against 78 s at 9 spins, 283 s against 364 s at 10 spins, 116 s against
185 s on the complex route, and 164 s against 308 s with 223 occupied
states. So "auto" takes the column form.

What is not done. The tip vertex is integrated as part of the state, so
its 13 stages are kept like those of the couplings. It does not enter
any right-hand side, so a quadrature that accumulates its stages would
need about 6 copies instead of 20, which matters when the number of
biases times the number of occupied states is large. The flow is in
float64 only.

## 8. Finite temperature

We will now see how the flow is taken to a temperature T. One
requirement fixes the construction: its first order in Jrho_s must be
the thermal third order of `conductance.py`. That third order uses the
thermal log F(x,T), the convolution of F0 with the weight
w(u) = Theta'(u/kT)/kT of `stepfunctions.FBuilder`. So the cutoff kernel
is replaced by its thermal average,

    K_T(x,D) = int du w(u) K(x+u,D),   int_0^inf dD K_T(x,D) = F(x,T),

and the Boltzmann populations p_i and the thermal step stay in the golden
rule, as at third order:

    dI/dV = 2 pi T0^2 sum_{i,f} p_i 4 sum_lam |A^lam_fi(eV)|^2 Theta_T(eV - eps_fi) + (s -> t).

Nothing else of the flow changes. The kernel is
`thermalkernel.ThermalKernel`, and `_joint_flow` takes it as an argument;
`leading_log_dIdV` builds it whenever `ks.T > 0`. The construction was
developed on 2026-11-11 in a project on 1T-TaSe2, as a script on top of
dmrgpy, and moved here on 2026-11-12. On the same inputs the two agree to
1e-15 (four models, both choices of the on-shell kernel, up to 16
occupied states).

Three things are not fixed by the first order. We say which choice is
taken.

- The kernel of the on-shell flow. By default (`onshell="thermal"`) K_T
  is also used for G and A_on, so the sample coupling stops growing at
  D ~ kT and not at Gamma0. `onshell="T0"` keeps the T=0 kernel there,
  and then the tip vertex is exactly the T=0 one convolved in bias with
  w. The two differ at second order in Jrho_s: in the TaSe2 project a
  free-spin peak at Jrho_s=-0.05 was 5% higher with the thermal choice
  (2.68 against 2.55 at 0.3 K).
- Amplitude or rate. The amplitude is averaged, then squared and
  multiplied by the thermal step, which is what `conductance.py` does at
  third order. Averaging the rate instead would differ at second order
  without a log. This was not computed.
- Which transitions are open, for the check of the running tip coupling.
  A thermal step is never zero, so at T>0 a transition is open when
  p_i Theta is above `open_cut` = 1e-2 of the largest p. Only which biases
  are warned or NaN depends on it.

A free spin 1/2 at B=0 has a closed form. Every process is elastic and
one coupling flows, and with the thermal kernel everywhere

    dI/dV(0)/(2 pi) = (3/4)/(1 + g F(0,T))^2,   F(0,T) = ln(omega0/(c kT)),

against (3/4)(1 - 2 g F(0,T)) at third order, with c = 1.30 for
Gamma0 << kT << omega0. With the Kondo sign, g<0, the peak grows on
cooling and diverges at kT_K = (omega0/c) exp(-1/|g|).

How the kernel is computed is in the module docstring of
`thermalkernel.py`: at each cutoff the even function x -> K_T(x,D) is
tabulated once (K_T itself within 48 kT, the relative correction
K_T/K - 1 beyond) and read from cubic splines; one table takes about
2 ms. Two shortcuts fail. A fixed short quadrature in u misses the peak
of K at u=-x once D and Gamma0 are below kT. And
Gamma0 -> sqrt(Gamma0^2 + (c kT)^2) has the wrong first order: in the
TaSe2 project the odd part of the difference to the third order dropped
by 1.95 instead of 8 when g was halved.

The checks are in `tests/test_kondo_leading_log_thermal.py`, in units of
an exchange J with kB=1, omega0=20 and Gamma0=2.5e-4. The numbers below
are measured with the same models.

1. The kernel alone. The integral of K_T over D equals F(x,T), from a
   direct quadrature of w F0, to 5e-11 relative at x = 0, Gamma0, kT and
   10 kT, for kT = 0.025 and 0.3, both through the closed-form integral
   of the tables and through K_T integrated over ln D. FBuilder differs
   from the same quadrature by 8e-7 to 2e-6, which is its own accuracy.
2. The limit T -> 0. With Jrho_s=-0.05 and U=0.25 on 12 biases the
   thermal route goes to the T=0 flow as (kT/Gamma0)^2: 6.4e-7, 4.7e-9,
   1.7e-10 and 2.3e-11 at kT/Gamma0 = 0.1 to 1e-4 for a spin 1/2 in a
   field, 2.8e-6 to 3.4e-11 for the S=1 model of section 5.
3. First order. Let D(g) be the leading-log result minus the second and
   third order of `conductance.py` at the same T. Halving g from 0.01
   divides the even part of D by 4.006 and the odd part by 8.011 (spin
   1/2, B=1, kT=0.3, U=0.25), 4.006 and 8.010 for S=1. A mismatch at first
   order would give 2 in the odd part. `onshell="T0"` gives 4.007 and
   8.014, as it should, since the on-shell kernel enters at second order.
4. The free spin. At Jrho_s=-0.05 the zero-bias peak follows the closed
   form to 1e-10 at kT = 0.0015, 0.025, 0.3 and 1 (2.5496, 1.6270, 1.1672
   and 1.0173, against 1.4365, 1.2316, 1.0476 and 0.9621 at third order),
   with c = 1.41, 1.30, 1.26 and 1.18.
5. The stop at the Kondo scale. A free spin at Jrho_s=-0.12 stops at
   D = 0.0128262 at T=0. At kT = 1e-4 Gamma0 it stops at the same cutoff,
   with the same NaN biases, and the finite ones agree to 1.6e-9. At
   kT = 0.001 and 0.005 it stops higher, and at kT = 0.025 it does not
   stop, since the temperature cuts the flow first.
6. A field along y takes the complex route and gives the curve of the
   same field along z to 7e-11, and three chunks of biases (`mem_gb`) give
   the result of one flow to 1.5e-11.

What it costs. The tip part is linear in the number of occupied states,
the ones within ln(1/p_cut) kT of the ground state, and this is what T
adds. Measured on d22, one core, through `leading_log_dIdV` on an open
chain of 6 spins 1/2 with exchange 1 meV and a tilted field
(0.3 meV along z, 0.1 meV along x on every site), Jrho_s=-0.05, U=0.25 and
61 biases: 1.1 s at T=0; 2.9 s at 1 K with
`p_cut=0.5`, so with the ground state only; 15 s at 1 K with its 27
occupied states, and 35 s and 33 s at 3 K and 10 K with all 64. So the
thermal kernel alone costs 2.6 times the T=0 one here, since its rows are
read from splines on one core. The memory of one flow is
`memory_model(d, m)`; `mem_gb` cuts the biases in chunks that fit, and
each chunk solves the couplings again.

What stays as at T=0: no Korringa broadening, no renormalized energies,
equilibrium populations, and the method fails at the Kondo temperature.

## 9. References

- A. Rosch, J. Paaske, J. Kroha, P. Woelfle, PRL 90, 076804 (2003),
  cond-mat/0202404: frequency-dependent couplings, the cut by Gamma.
- A. Rosch, J. Paaske, J. Kroha, P. Woelfle, JPSJ 74, 118 (2005),
  cond-mat/0408506: the general one-loop equation (eq. 6), the
  approximations that close it, the S=1/2 equations in a field.
- A. Rosch, T. A. Costi, J. Paaske, P. Woelfle, PRB 68, 014430 (2003),
  cond-mat/0301106: the same equations in equilibrium, their solution,
  the check against bare perturbation theory to g^3 ln^2, and the
  comparison with NRG; the benchmark for a field.
- J. Paaske, A. Rosch, J. Kroha, P. Woelfle, PRB 70, 155301 (2004),
  cond-mat/0401180: the decoherence rates that cut the logs.
- M. Ternes, NJP 17, 063016 (2015), arXiv:1505.04430: the model and the
  third order this method starts from.

## Status

2026-11-06: `leadinglog.py` implemented and validated on the ED
`KondoSpectrum` at T=0 (section 5), and exposed as
`get_kondo_spectrum(order="LL")` on `mode="ED"` with `T=0` (user guide,
STM and Kondo spectra, "Leading logs"). 2026-11-12: the thermal kernel
(section 8), at any T on `mode="ED"`. Open: Korringa rates
in place of the fixed Gamma0, and the DMRG route. For the DMRG route
the TASK spectrum has the right shape, but the sum over intermediate
states is then truncated to the Krylov space, and the flow needs
products of spin operators that the space built from S_k|GS> does not
hold; a block-Krylov space from S^a S^b|GS> is the likely fix, and it
is not tested.

2026-11-11: the flow is solved as one system on the columns of the
occupied states (section 7). No number moved beyond the tolerance of
the integration. `tip_vertex` takes any mix of U now, and `cols=`.
The tip coupling check reads the open transitions only, with a floor
on the bare vertex (end of section 7); this changes which biases are
NaN or warned, not the values. The exchange sum stays on one core, by
the author's choice on this date.

2026-11-12: the flow runs on a GPU through JAX (`leadinglog_device.py`,
end of section 7), with the numbers of the NumPy route to roundoff. On
the RTX A2000 of d22 it is 2 to 4 times faster than six cores at 9 and
10 spins. The H200 timing at 11 to 13 spins is pending.
