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
below include this build.

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

## What is not here

The comparison is at particle-hole symmetry only, where the reference at
$U=0$ is unambiguous. Away from it the reference level of the impurity
(the bare level or a Hartree level) is a choice that was not measured.
Neither the natural-orbital iteration (the basis of the correlation
matrix of the interacting state, `cm=`) nor a bath with pairing was
tried. The cost of the MPO build could be cut by building the MPO of a
quadratic Hamiltonian directly, which was not done.
