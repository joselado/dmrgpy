# The Gaussian basis against the site basis

This folder holds the measurement behind the guide section "A basis for
impurity models" (`docs/guide/solver.md`): DMRG of the single-impurity
Anderson model in the site basis and in the Gaussian basis of
`get_gaussian_basis` (arXiv:2212.09798), for the same Hamiltonian. The
question is how much bond dimension, and how much time, the rotation saves
at a given accuracy, and what it costs.

`compare.py` runs one case and writes `data/cmp_<case>_n<n>_U<U>.json`;
`plot_compare.py` reads every file in `data/` and draws
`compare_energy.png` and `compare_entropy_mpo.png`. The runs were:

```bash
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
taskset -c 0 python3 compare.py uniform 32 4.0 8,16,32,64,128
taskset -c 0 python3 compare.py uniform 64 4.0 8,16,32,64,128
taskset -c 0 python3 compare.py uniform 128 4.0 8,16,32,64,128
taskset -c 0 python3 compare.py wilson 32 0.5 8,16,32,64,128
python3 plot_compare.py
```

## The setup

The impurity is orbital 0 of a `Spinful_Fermionic_Chain`, with
$U(n_\uparrow-\frac12)(n_\downarrow-\frac12)$, so it is particle-hole
symmetric at half filling. The bath is either a uniform chain (hopping
$t=1$, impurity coupling $V=t$, $U=4$) or a Wilson chain
(`wilsonchain.wilson_chain`, $\Lambda=2$, band $[-1,1]$, $t_D=0.2$,
$U=0.5$). The Gaussian basis is built from the ground state at $U=0$,
with the impurity fixed and no limit on the block size. Every solve is
`itensor_version=3`, 12 sweeps, a cold start, the same schedule in both
bases. The error is $E-E_{\rm best}$, with $E_{\rm best}$ the lowest
energy of all the runs of a case, which is always the Gaussian basis at
`maxm=128`.

Two checks say that this is a fair reference. More sweeps do not close the
gap of the site basis: on 32 orbitals at `maxm=64` its error is 5.34e-3
after 12 sweeps, 5.30e-3 after 30 and 5.28e-3 after 60, so what limits it
is the bond dimension. And at $U=0$, where the exact energy is twice the sum
of the filled levels, the Gaussian basis on 128 orbitals is 6.0e-7 above it
at `maxm=8`, 2.4e-7 at 16 and 7.3e-8 at 32, so the rotated Hamiltonian is
right at the largest size of the figure.

## What we see

The energy error against the bond dimension (`compare_energy.png`, top
row):

| orbitals | site, maxm=64 | site, maxm=128 | Gaussian, maxm=8 | Gaussian, maxm=16 | Gaussian, maxm=32 |
|---|---|---|---|---|---|
| 32 (uniform) | 5.5e-3 | 3.8e-4 | 4.8e-3 | 1.2e-4 | 2.0e-6 |
| 64 (uniform) | 2.4e-2 | 3.4e-3 | 5.2e-3 | 1.7e-4 | 4.5e-6 |
| 128 (uniform) | 4.8e-2 | 1.5e-2 | 5.3e-3 | 2.0e-4 | 6.3e-6 |

In the site basis the error at fixed `maxm` grows with the length of the
bath, because the Fermi sea of the bath is entangled at every cut. In the
Gaussian basis it hardly depends on the length: the bath far from the
impurity is a product state, and the entanglement sits on a few orbitals
next to it. The entanglement profile on 128 orbitals
(`compare_entropy_mpo.png`, left) shows this directly: a median of 1.9 over
the cuts in the site basis (still growing with `maxm`), and in the
Gaussian basis 1.2 at the impurity, 0.2 after its partner orbital, and
below 1e-3 from about 15 orbitals on.

The price is the Hamiltonian. The MPO of the rotated $H$ has a bond
dimension of 34, 42 and 46 on 32, 64 and 128 orbitals, against 6 in the
site basis, growing roughly as the logarithm of the size
(`compare_entropy_mpo.png`, right). Its build is also slow: on 128
orbitals the rotated $H$ has 32020 terms, and building its MPO takes 93 s
before the first sweep (63 s in the AutoMPO of ITensor v3, 28 s in the
Jordan-Wigner expansion of `MultiOperator.to_terms`). The wall times
below include this build, and were taken before the two changes that
followed this comparison: the Jordan-Wigner strings are now built in one
step (13.8 s to 3.2 s for that expansion), and on v3 a Hamiltonian with
that many string factors keeps the names of its hoppings, so that
ITensor's AutoMPO threads the strings itself
(`groundstate.native_strings`). With both, the build on 128 orbitals
takes 11 to 18 s instead of 94 s, and the Gaussian basis at `maxm=8`
took 43 s in total instead of 128 s, with the same MPO dimension. The
two routes give energies 6e-8 apart at `maxm=16`, inside the solver
noise.

At equal accuracy the Gaussian basis is faster from 64 orbitals on. On
64 orbitals the site basis at `maxm=128` took 239 s and is 3.4e-3 off;
the Gaussian basis at `maxm=16` took 23 s and is 1.7e-4 off. On 128
orbitals the site basis at `maxm=128` took 424 s and is 1.5e-2 off; the
Gaussian basis at `maxm=8` took 128 s and is 5.3e-3 off, and at `maxm=16`
158 s and 2.0e-4. On 32 orbitals the two are close in time at the
accuracy the site basis reaches.

On the Wilson chain the gain is a factor 2 in `maxm`: the Gaussian basis
at `maxm=32` is 2.9e-5 off against 9.3e-5 for the site basis at
`maxm=64`, and 1.5e-6 at `maxm=64` against 1.2e-6 at `maxm=128`. In time the two are
close: the Gaussian basis at `maxm=32` took 12.5 s against 14.1 s for
the site basis at `maxm=64`, and 49 s at `maxm=64` against 81 s at
`maxm=128`.

The timings are indicative. The machine was shared with other jobs during
the runs, the solves were pinned to one performance core with one BLAS
thread, and the 128-orbital runs overlapped with short checks on another
core.

## The reference level and the natural orbitals

`reference.py` asks which reference gives the best basis, with $U$ written
as $Un_\uparrow n_\downarrow$, as `wilsonchain.anderson_chain` does. The
bare level $\epsilon_d$ of the state without $U$, the level
$\epsilon_d+U/2$ (the Hartree level of a half-filled impurity, the Fermi
level at particle-hole symmetry), the restricted self-consistent Hartree
level $\epsilon_d+U\bar n$, and the natural orbitals of the interacting
state, which `GaussianBasis.natural_basis` builds from the correlation
matrix of the state solved at `maxm=32` in the $\epsilon_d+U/2$ basis,
at three values of `eigval_cutoff`. 32 orbitals, v3, 12 sweeps, cold
start, the error against the lowest energy of all the runs of a case
(`reference_table.py` prints the tables):

```bash
python3 reference.py uniform 32 4.0 -2.0 8,16,32    # particle-hole symmetric
python3 reference.py uniform 32 4.0 -1.0 8,16,32
python3 reference.py wilson 32 0.5 -0.25 8,16,32    # particle-hole symmetric
python3 reference.py wilson 32 0.5 -0.1 8,16,32
python3 reference.py wilson 32 0.5 -0.1 8,16,32 30  # 30 sweeps
python3 reference_table.py
```

| case | basis | MPO | maxm=8 | maxm=16 | maxm=32 |
|---|---|---|---|---|---|
| uniform, $\epsilon_d=-U/2$ | bare | 34 | 8.2e-3 | 6.9e-4 | 1.1e-5 |
| | $\epsilon_d+U/2$ (= Hartree) | 34 | 4.8e-3 | 1.2e-4 | 2.0e-6 |
| | natural, cutoff 1e-8 | 38 | 5.2e-4 | 3.1e-5 | 6.3e-7 |
| | natural, cutoff 1e-5 | 26 | 7.5e-4 | 8.9e-5 | 4.4e-6 |
| | natural, cutoff 1e-3 | 22 | 5.8e-3 | 7.8e-4 | 7.2e-5 |
| uniform, $\epsilon_d=-U/4$ | bare | 34 | 6.2e-3 | 4.6e-4 | 7.9e-6 |
| | $\epsilon_d+U/2$ | 34 | 2.2e-3 | 6.5e-5 | 1.3e-6 |
| | Hartree | 34 | 3.5e-3 | 7.2e-5 | 1.1e-6 |
| | natural, cutoff 1e-8 | 38 | 3.8e-4 | 3.0e-5 | 4.5e-7 |
| | natural, cutoff 1e-5 | 26 | 5.6e-4 | 8.2e-5 | 2.7e-6 |
| | natural, cutoff 1e-3 | 22 | 5.5e-3 | 8.0e-4 | 5.6e-5 |
| Wilson, $\epsilon_d=-U/2$ | bare | 26 | 2.2e-3 | 6.8e-4 | 7.3e-5 |
| | $\epsilon_d+U/2$ (= Hartree) | 30 | 4.0e-3 | 2.5e-4 | 2.0e-5 |
| | natural, cutoff 1e-8 | 66 | 1.5e-3 | 2.4e-4 | 2.3e-5 |
| | natural, cutoff 1e-5 | 46 | 2.8e-3 | 1.7e-4 | 1.8e-5 |
| | natural, cutoff 1e-3 | 26 | 1.5e-3 | 2.0e-4 | 4.4e-5 |
| Wilson, $\epsilon_d=-U/5$ | bare | 26 | 3.0e-3 | 9.5e-4 | 1.0e-4 |
| | $\epsilon_d+U/2$ | 26 | 1.5e-3 | 3.9e-4 | 2.2e-5 |
| | Hartree | 30 | 2.7e-3 | 1.2e-4 | 0 |
| | natural, cutoff 1e-8 | 66 | 8.2e-3 | 8.3e-4 | 7.7e-5 |
| | natural, cutoff 1e-5 | 38 | 1.0e-3 | 8.5e-5 | 2.6e-6 |
| | natural, cutoff 1e-3 | 22 | 8.0e-4 | 1.2e-4 | 1.5e-5 |

The site basis is 4e-2 to 5e-1 off at `maxm=8` and 2e-3 to 4e-2 at
`maxm=32` in every case, so every Gaussian basis is far ahead of it. On the uniform bath
$\epsilon_d+U/2$ is the best reference at `maxm=8` and 16, and within
20% of the self-consistent Hartree level at 32; the bare level is 1.7 to
7 times worse. On the Wilson chain at particle-hole symmetry
$\epsilon_d+U/2$ is the Hartree level, and the bare level is better at
`maxm=8` and 2.7 to 3.7 times worse at 16 and 32. On the Wilson chain away
from particle-hole symmetry the self-consistent Hartree level is better at
`maxm=16` and 32 (1.2e-4 and 6.1e-6 with 30 sweeps, against 3.9e-4 and
2.8e-5 for $\epsilon_d+U/2$), and $\epsilon_d+U/2$ only at `maxm=8`
(1.5e-3 against 2.7e-3). So which level is best away from particle-hole
symmetry is not settled, and a reference for two levels was not measured
at all. $\epsilon_d+U(2L-1)/2$ is the default of
`anderson_chain(basis="gaussian")` for now.

The natural orbitals depend on the cutoff, since the interacting state is
not Gaussian and the blocks of the circuit grow over its correlated part
when the cutoff is that of a pure state. With 1e-8 the uniform bath gains
6 to 9 at `maxm=8` and 2 to 4 at `maxm=16`, but the Wilson chain away
from particle-hole symmetry loses 5 at `maxm=8` (2.4 with 30 sweeps, see
below) and the MPO is 66. With
1e-5, the default of `natural_basis`, the error at `maxm=8` is smaller
than with $\epsilon_d+U/2$ in all four cases, by 1.4 to 6, and at
`maxm=16` in three of them (by 1.3 to 4.6, and 1.3 times larger on the
uniform bath away from particle-hole symmetry); at `maxm=32` it is about 2
times larger on the uniform bath and 1.1 to 8 times smaller on the Wilson
chain. A second step of natural orbitals, from the state of the first,
gives nothing more (`natural2` in the data).

Twelve sweeps are enough for the references, but not for the larger MPO
of the natural orbitals at cutoff 1e-8 on the Wilson chain. With 30 sweeps
(`python3 reference.py wilson 32 0.5 -0.1 8,16,32 30`) the references move
by less than 1e-5, while that basis goes from 8.2e-3 to 3.6e-3 at
`maxm=8` and from 8.3e-4 to 3.2e-4 at `maxm=16`. It still loses at
`maxm=8`, by 2.4 against $\epsilon_d+U/2$, and cutoff 1e-5 stays the best:
7.6e-4, 8.7e-5 and 7.8e-6 at `maxm=8`, 16 and 32, against 1.5e-3, 3.9e-4
and 2.8e-5. The measured spin-mixing
block that `natural_basis` drops was 2e-8 to 2e-4.

## A deep Wilson chain

With `anderson_chain(21, ed=-0.25, U=0.5, tD=0.2, Lambda=4.0)` (22
orbitals, Wilson sites down to $4^{-10}$ of the band), the Gaussian basis
gives no gain. On v3 with 40 sweeps, the site basis gives -1.9874176475
at `maxm=64` and -1.9874188004 at 128, and the Gaussian basis
-1.9873881894 at 32, -1.9874177381 at 64 and -1.9874187858 at 128. Its
largest bond entropy is 1.13 against 2.08, but the energy at equal `maxm`
is the same within 1e-7. At `maxm=128` the two bases agree to 1.5e-8, and
`get_kondo_screening` gives the same $M_{\rm imp}=0.843145$ and Schmidt
values $(0.70711, 0.70711)$, so the rotation loses no low-energy hopping.
Measuring the screening takes 8 to 9 s in the Gaussian basis against 1 to
2 s, since each bath $S^z_n$ becomes a sum over all pairs of rotated
orbitals.

## What is not here

The first comparison is at particle-hole symmetry only; the reference
level and the natural orbitals away from it are in the section above, on
32 orbitals only. A bath with pairing was not tried. The MPO build, 11 to 18 s on 128 orbitals now, could be cut further
by building the MPO of a quadratic Hamiltonian directly, which was not
done.
