# The thermal third-order Kondo spectrum from Krylov spaces

This folder documents `get_kondo_spectrum(mode="ED", submode="KRYLOV")`,
the third-order STM spectrum of Ternes (NJP 17 063016, arXiv:1505.04430)
at a temperature T without a full diagonalization. The code is
`src/dmrgpy/kondospectrumtk/krylovkondo.py`, the tests are
`tests/test_kondo_krylov.py`, and the example is
`examples/kondo/kondo_krylov_thermal_VS_ED`. The route was developed in a
project on 1T-TaSe2 (November 2026) and moved here.

## 1. Which states the sums need

The third-order sums of `conductance.py` run over an initial state i, a
final state f and an intermediate state m:

    sum_i p_i sum_{f,m} <i|S_a|f><f|S_b|m><m|S_c|i> ...

`mode="ED"` takes all three from the full eigenbasis, so it diagonalizes
the whole Hamiltonian. It refuses a dense diagonalization above 2000
states (`algebra.maxsize`), which is 10 spins 1/2. Two facts make most of
that basis unnecessary.

- The initial states carry the Boltzmann weight p_i. At a temperature T
  only the states within about ln(1/p_cut) kT of the ground state matter.
  We find them with a sparse eigensolver, and raise the number of states
  until the relative weight of the last one is below `p_cut` (1e-8 by
  default).
- f and m enter through <f|S_a|i> and <m|S_c|i>, with S the spin of the
  tip site. So for each initial state they lie in the smallest H-invariant
  space that holds S_x|i>, S_y|i> and S_z|i>. A Krylov space started from
  these three vectors is this space once it is exhausted, and then its
  Ritz states are exact eigenstates of H. The middle vertex <f|S_b|m> is
  taken between two such states, which is exact for them.

So for each initial state the sums run over a basis of at most
3 `nkry` + 1 states, and the result is the full ED one once every Krylov
space is exhausted, up to the initial states dropped by `p_cut`.

## 2. One block Krylov space, not three

When H conserves Sz, S+|i>, S-|i> and Sz|i> lie in three different Sz
sectors, and three separate Lanczos spaces are orthogonal to each other.
A field that breaks Sz conservation, for example a field that turns from
site to site, mixes the sectors. The three spaces then overlap, and three
separate sets of Ritz states count the shared directions two or three
times. This gave errors of order one in the TaSe2 project, and the test
`test_three_separate_spaces_double_count_with_a_field` reproduces it:
on 5 spins with a 120 degree field in the xz plane the separate spaces are
2.7 off in relative terms, while without the field they agree with ED to
5e-11.

The fix is one block Krylov space started from the three vectors together.
Each new block H B is projected twice on the whole space so far (block
Gram-Schmidt, one product with the basis per pass), and an SVD
keeps its directions with a singular value above `DROP_TOL` (1e-9) of its
largest column norm. The space stops when a block is empty, which means
the space is exhausted, or when it would pass 3 `nkry` states. Then one
Rayleigh-Ritz step in the whole space gives the Ritz states. With a field
in the xz plane H is real, the start vectors S+|i>, S-|i>, Sz|i> are real,
and everything stays real.

Nothing of the size of the Hilbert space is built twice. The coefficients
of the projections give Q^T H Q, the residual of each Ritz state comes from the
same recurrence (the part of H B outside the space for the last block),
and the spin matrices are taken in Krylov coordinates, Q^T S Q rotated by
the Ritz vectors of the small problem, so the Ritz vectors are never formed
in the Hilbert space. In a real basis the matrix of S- is the transpose of
that of S+. The matrices then come out with Sx and Sz real and Sy i times
a real matrix, so the sums of `conductance.py` run in real arithmetic.

The space of a state can also stop before 3 `nkry` states, but only when
it is asked for (`tol`, None by default). A final state f enters the sums
only while its step is open, e_f - e_i below max|eV| (plus 40 kT at T>0),
and the stop asks that every Ritz state in that window be converged: its
residual ||(H - theta) f|| below `tol` times onenorm(H). The check is a
diagonalization of the small matrix, which is complete after each block,
so it costs no product in the Hilbert space; it is made after 10 blocks
and then each time the space has grown by a factor 1.25, or by 2 while
fewer than half of the states of the window are converged. When the stop
fires the result is the one of a fixed space of that size. It fires only
where the window is narrow against the spectrum the space reaches. On 8
spins in a field with `nkry=60` and `tol=1e-8` it stopped at 121 of 181
states at 0.1 K and |eV| < 0.5 meV and at 31 at T=0 and 0.3 mV, and dI/dV
moved by 4e-12 and 1e-11 of the peak; at 1 K and 4 mV it never stopped.
On 10 spins at 0.2 K and |eV| < 1 mV it came at 373 of 1024 states (measured before the
block Gram-Schmidt), above
the 300 of the default `nkry`, so it is a way to make a large `nkry`
cheap where it is not needed, not a gain at the default.

## 3. How the sums are assembled

The initial state is projected out of the starting block, and the Krylov
space stays orthogonal to it. Since H|i> = E_i|i>, this costs nothing:
H maps the complement of i to itself. So [i] + the Ritz states is an
orthonormal basis, with i at its own energy. In this basis we build an
object with the fields that `conductance.py` reads (e, p one-hot on i,
Sx, Sy, Sz, T, kB, dim), the way `taskkondo.TaskKondoSpectrum` does, and
call `second_order_dIdV`, `third_order_kondo_dIdV` and
`third_order_potential_dIdV` unchanged. Every prefactor, the argument of
the exchange diagram and the potential-interference term (U) are then
those of `mode="ED"`. The total is the sum over initial states with their
weights, since every term is linear in the populations. The elastic
transition is in the basis as row 0, and the other members of a
degenerate level are reached through the Ritz states. A Ritz energy within
the degeneracy tolerance of E_i is set to E_i exactly, as `taskkondo.py`
does, so that the step at zero bias is read at Theta0(0) = 1/2.

At T=0 the initial states are the ground manifold with equal weights, as
in `edkondo.KondoSpectrum`. That class takes the states within 1e-9 of
the full width of the spectrum. Here the width is not known, and we use
2 onenorm(H) instead, a bound of the width.

## 4. The SU(2) reduction

When H commutes with the total spin the route takes one member of each
spin multiplet, the one at S^z_tot = 0 for an integer total spin and at
S^z_tot = +1/2 for a half-integer one (a chain of spinful fermions with no
Nf sector has both kinds), with the weight (2S+1) exp(-(E-E0)/kT). For an
unpolarized tip each term of dI/dV is the expectation value of a
rotation-invariant operator in the initial state, so all members of a
multiplet give the same result, and this is the sum over all of them.
The code checks that H commutes with the total spin on a random state,
one component at a time. Inside a degenerate level it rotates the states
to a definite S with S^2_tot, and the check fails when 2S is not an
integer. With `su2=None`, the default, a failed check means the full sum;
`su2=True` raises instead (`SU2Error`), and `su2=False` always takes the
full sum. On a 7-spin Heisenberg chain at 1.5 K the reduction takes 32
initial states instead of 110 at `p_cut=1e-14`, and the result is the same
to 3e-15 of the peak. On the 12-spin chain of section 6 at zero field and
1 K it takes 19 multiplets instead of 47 states, 19 s instead of 66 s,
and the two agree to 2.5e-12, the Krylov error at `nkry=100`. With a field
the check fails at once (9 ms on 12 spins).

The weight cut is the one of the full sum: a multiplet is kept when
exp(-(E-E0)/kT) is above `p_cut`, whatever its 2S+1, so the two routes
take the same states. A cut on (2S+1) exp(-(E-E0)/kT) would keep the
multiplets of larger spin just below the cut of the full sum: on the
12-spin chain an S=2 and an S=1 multiplet at 18.5 kT, which move dI/dV by
2.4e-8 of the peak, against 3e-12 between `nkry=100` and 200.

## 5. Checks

Each check is in `tests/test_kondo_krylov.py`, against `mode="ED"` with
`p_cut=1e-14`, on chains of 4 to 8 spins with exchange near 1 meV,
Jrho_s=-0.08, U=0.3, omega0=20 meV and Gamma0=50 microeV. The numbers are
the largest difference over the bias window divided by the largest dI/dV.

1. A Heisenberg chain at T=0, with 6 spins (a singlet) and 5 spins (a
   ground doublet, averaged as on ED): 1e-15 and 6e-16.
2. No Sz conservation. A 120 degree field in the xz plane (H real) at
   T=0 and 2 K: 1e-15. A field along y (H complex) at 1 K: 1e-12. The
   last one is the weight cut of the ED sums (`conductance.P_CUT` =
   1e-12), which keep fewer initial states than this route at
   `p_cut=1e-14`.
3. An S=1 impurity with easy-plane anisotropy in a tilted field,
   coupled to 3 spins 1/2, with order 2 and 3 at T=0 and 1.5 K: 2e-15 to
   3e-15.
4. The default `p_cut=1e-8` on the chain in a field at 2 K: 6e-9, the
   size of the cut.
5. The sparse eigensolver on an 8-spin Heisenberg chain, where every
   level is a multiplet, at 0.5 K with eight initial states: 2e-13. A
   member of a multiplet that eigsh missed would show here.
6. `su2=True` on 7 spins at T=0 and 1.5 K: 1e-15 and 4e-12; and the
   raise with a field. The default `su2=None` against `su2=False` at T>0,
   on 7 spins at 1.5 K and on a 3-site Hubbard chain with no Nf sector:
   below 3e-15 at `p_cut=1e-14` and at the defaults; with a field it takes
   the full sum. The cut of a multiplet: with `p_cut` just above the
   Boltzmann factor of a level of S >= 3/2 both routes drop it and keep
   the same states.
7. The spin of a Hubbard orbital (3 sites, U=8, a tilted field, the
   half-filled Nf sector) at kT=0.3: 1.2e-13.

The TaSe2 project checked the same route against full ED on its own
models: 6e-6 on 12 sites with the SU(2) reduction, 2e-9 on the sqrt7 and
3x3 moire tori, and with a 120 degree mean field at 0.3 K, 3.1e-11 on a
trimer and 1.9e-10 on 9 sites.

## 6. Cost

Measured on one core, with 61 biases, on an open chain of spins 1/2 with
exchange 1 meV (times 1 + 0.3 cos i), a 120 degree field of 0.2 meV in the
xz plane, T=1 K, `nkry=100` and the default `p_cut`:

- 10 spins (1024 states): 28 initial states, 24 s, against 33 s for full
  ED, and the two differ by 6.8e-9.
- 12 spins (4096 states): 49 initial states, 65 s, of which 30 s are the
  sums and the rest the Krylov spaces. Full ED refuses this size
  (`algebra.maxsize` = 2000); with the limit raised it had not finished
  after 43 minutes on one core, and it was stopped, so at this size the
  agreement with ED is not measured here. Since 2026-11-14 a Krylov space
  takes 0.60 s per initial state (0.67 s before), and 0.39 s with the
  block Gram-Schmidt, and 90% of the
  third-order term of one state is the two tables of the exchange log F,
  one per tunneling direction, each over all the final and intermediate
  states of the space. Since F is even, the table of -eV on final x final
  states is the transpose of the one at +eV, so only its columns over the
  states that are not final are evaluated; on 10 spins at 1 K, where every
  state of the space is final, this took the third-order Kondo term of one
  state from 0.83 to 0.49 s on one core of the workstation. On one core
  of Triton the whole 12-spin call now takes 36.4 s, against 52.5 s
  before the block Gram-Schmidt and this table on the same node.

Per initial state the cost is the block Krylov space (3 `nkry` sparse
products and the reorthogonalization, O(dim nkry^2)) and the sums of
`conductance.py` on a basis of 3 `nkry` + 1 states, which at T>0 means
about (3 `nkry`)^2 evaluations of F per bias. So the time is linear in the
number of initial states, which is what T sets, and close to linear in
the dimension. The memory is that of a few vectors of the Hilbert space
per Krylov state, 3 `nkry` dim numbers.

## 7. What the route does not contain

- The leading logs. The flow of `order="LL"` sums products of spin
  operators between any two states, so it needs the full eigenbasis, and
  `submode="KRYLOV"` with `order="LL"` raises. On 13 sites this is what
  blocks the leading logs in the TaSe2 project.
- A convergence report. The residual of each Ritz state is stored
  (`StateSpectrum.residual`), but no warning reads it, unlike TASK.
  Increase `nkry` until the spectrum stops changing.
- A check that eigsh found every member of a degenerate level. ARPACK
  finds the multiplicity from rounding; on the SU(2) chains of the tests
  it did. With the reduction the levels in one Sz sector are less
  degenerate, which makes this safer.
