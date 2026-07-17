# TDVP
ITensor implementation of the time dependent variational principle (TDVP) algorithm for finite MPS.
It also includes a global subspace expansion algorithm to enlarge the bond dimension of the MPS for TDVP according to the paper [arXiv:2005.06104](https://arxiv.org/abs/2005.06104) or [Phys. Rev. B 102, 094315 (2020)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.094315), which would be useful for reliable time evolution of two-dimensional systems and one-dimensional systems with long-range interactions, especially for the initial stage of some time evolution when the bond dimension of the MPS is small.

Requires a working version of the ITensor (C++) library on your system.

See the sample/ directory for an example of how to run it.

If you have any questions about this repository, please contact <mingruy@uci.edu>.

### TDVP

#### Function
``tdvp(MPS psi, MPO H, Cplx t, Sweeps sweeps, Args args) -> Real energy``

``tdvp(MPS psi, MPO H, Cplx t, Sweeps sweeps, DMRGObserver obs, Args args) -> Real energy``

Note there are other interfaces available for TDVP, which are similar to their [DMRG counterparts](http://itensor.org/docs.cgi?page=classes/dmrg&vers=cppv3).

#### Parameters

`psi`: the MPS to be time evolved.

`H`: the MPO of the Hamiltonian. Currently only Hermitian Hamiltonians can be treated; for non-Hermitian Hamiltonians, please merge this [pull request](https://github.com/ITensor/ITensor/pull/410) to your ITensor and set the argument `IsHermitian` to be `false` in the `tdvp` function. (In this pull request, the `applyExp` function is modified to include the Arnoldi orthogonalization so as to treat the non-Hermitian Hamiltonians. In addition, it uses a [time-step adjusting technique](https://www.maths.uq.edu.au/expokit/paper.pdf), so if your time step is too large for the `MaxIter` (controlled by `niter` of `sweeps`, see below) argument, the algorithm will automatically reduce the time step and restart the integration for multiple times until either the requested time step is reached or the `MaxRestart` is reached. This strategy is useful when you have a limited memory resource, since you can restrict your `MaxIter` to a smaller number and thus less memory is required, though same effect can be achieved by requesting a smaller time step.)

`t`: the time step of TDVP. It can be real, imaginary, or complex. The corresponding time evolution operator of a single time step will be $e^{tH}$. Therefore, to do real time evolution, `t` need to be purely imaginary; to do imaginary time evolution, `t` need to be purely real.

`sweeps`: Specify the sweep parameters of TDVP (similar to DMRG). `nsweeps` is the number of TDVP sweeps. A TDVP sweep = a sweep from left to right with half time step + a sweep from right to left with half time step. The total evolution time = `t*nsweep`. `maxdim` is the maximum bond dimension of the sweep. `cutoff` is the truncation error of the sweep (to allow truncation for the one-site TDVP, one needs to set the `Truncate` `args` to `true` (see below)). `niter` is the maximum number of lanczos iterations used when solving each local effective TDVP equations.

`obs`: is the observer one can customize to do measurement after each sweep without recalling the `tdvp` function. Similar to its use in DMRG.

`args`: `NumCenter` can either be 1 or 2, corresponding to the one-site and two-site TDVP respectively (default is `2`). `Truncate` choose whether or not truncate when doing the SVD (for `NumCenter=1`, default is `false`; for  `NumCenter=2`, default is `true`). `DoNormalize` choose whether or not normalize the MPS after a TDVP sweep (default is `true`). `Quiet` choose to whether or not print out the information of each local update. If `WriteDim` is specified, then the environment tensors PH's for the TDVP sweep will be written to disk when the bond dimension of the MPS is larger than `WriteDim`. `WriteDir` gives the directory to write those environment tensors PH's (default is the current directory).


### Global subspace expansion

#### Function
``addBasis(MPS phi, MPO H, vector<Real> truncK, Args args)``

``addBasis(MPS phi, MPO H, vector<int> maxdimK, Args args)``

#### Parameters
`phi`: the MPS will be global subspace expanded after calling the function.

`H`: the MPO of the operator to be used to construct the subspace, e.g. the Hamiltonian $H$ or the operator $1-\mathrm{i}\tau H$.

`truncK` `maxdimK`: one can choose to either specify the truncation error `truncK` or the maximum bond dimension `maxdimK` when applying the MPO `H`. Each element in the vector corresponds to each order of application. For example, `truncK={1e-8,1e-6}` means the truncation error of `H` applying to `phi` is 1e-8 and the truncation error of `H` applying to `H*phi`(obtained from previous step) is 1e-6. `maxdimK` is useful to make sure the bond dimension growth when applying the MPO `H` not out of control. When the bond dimension of the MPS `phi` is not small, a **typical strategy** is to set `maxdimK` to be the same as `maxLinkDim(phi)` and tune the `KrylovOrd` and `Cutoff` args in the `addBasis` function and the size of the time step `t` in the `tdvp` function accordingly to get optimal performance. For a Hamiltonian `H` with long-range interactions, usually at some point during the time evolution when the bond dimension becomes large enough (one needs to test what is the minimal bond dimension to turn off GSE), the global subspace expansion could be **turned off** and we can switch to the ordinary two-site TDVP without global subspace expansion.

`args`: `Cutoff` set the truncation error when diagonalizing the sum of the reduced density matrices (default is `1e-15`). `KrylovOrd` is the dimension $k$ of the Krylov subspace $\{\phi,H\phi,...,H^{k-1}\phi\}$ (default is `2`). `Method` specify which method is used when applying the MPO, it can be e.g. `DensityMatrix` or `Fit` (default is `DensityMatrix`, the `DensityMatrix` way to apply the MPO will produce more relevant basis than the `Fit` way given the same bond dimension, though at a [higher cost](http://itensor.org/docs.cgi?vers=cppv3&page=classes/mps_mpo_algs)). `DoNormalize` choose whether or not to perform normalization after applying the MPO (default is `false`). `Nsweep` specify the number of sweeps if use the `Fit` method to apply MPO (default is `2`). `Quiet` choose whether or not print out the norm and bond dimension of the Krylov vectors and the warning messages, (default is `false`). If `WriteDim` is provided, `DensityMatrix` way of applying the MPO will write the intermediate environment tensors E's to disk when the bond dimension of the MPS is bigger than `WriteDim`, thus reducing the memory usage. `WriteDir` gives the directory to write those environment tensors E's (default is the current directory).
