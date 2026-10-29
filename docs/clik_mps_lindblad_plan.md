# Lindbladian spectra and steady states with complex-time Krylov spaces

Plan, 2026-10-29. Source: P. Westhoff, M. Moroder, U. Schollwöck and
S. Paeckel, "A Tensor Network Framework for Lindbladian Spectra and Steady
States", arXiv:2509.07709 (CLIK-MPS). The plan was written from the main text
and Appendices B, F, G and H. The last part of the paper and parts of
Appendices C to E were not read.

**Status: Phase 0 fixed on "python" (2026-10-29).** `tdvp_step` takes an
Arnoldi exponentiator for a non-Hermitian MPO, and "python" now matches v3
(9.0e-7 from the exact curve on both). This is pinned by
`tests/test_nh_time_evolution.py`, recorded in `docs/numerical_changes.md`,
and checked under the JAX backend (1e-13 from NumPy). Three things are open
in Phase 0: the ED-against-DMRG norm convention (decision 1 in section 8),
the `julia_live` check, and an audit record. Phase 1 can start.

## 1. The idea

Let us start with the question the paper answers. We have a generator $G$
whose spectrum we want near its slowest-decaying end: a Lindbladian
$\mathcal{L}$, or $-H$ for a non-Hermitian Hamiltonian. Its steady state and
its first few modes decide the long-time physics: the dissipative gap
$\Delta = -\mathrm{Re}\,\lambda_2$ sets the relaxation time, and its scaling
with $L$ gives the dynamical exponent.

The way the paper does this is by building a Krylov space out of time-evolved
MPS, not out of powers of $G$. One state evolves in real time,
$\rho(t_m) = e^{G m\,\delta t}\rho_0$, and this branch resolves the steady
state. A second, traceless state evolves along a complex contour,
$z = t\,e^{-i\alpha}$. Under it, mode $k$ decays at the rate

$$\xi_k = -\mathrm{Re}\,\lambda_k \cos\alpha - \mathrm{Im}\,\lambda_k \sin\alpha ,$$

meaning that modes with $\mathrm{Im}\,\lambda_k > 0$ decay slower and stay in
the space long enough to be resolved. Modes with $\mathrm{Im}\,\lambda_k < 0$
decay faster, but a Lindbladian preserves Hermiticity, so they are the complex
conjugates of the upper ones and need no extra work. Above an angle
$\alpha_{\max}$, set by the bulk of the spectrum, some modes grow without
bound, so in practice $\alpha \ll \alpha_{\max}$.

Then we project. With the Gram matrix $M_{ij} = \langle\langle\rho_i|\rho_j\rangle\rangle$,
$M = U S U^\dagger$, $X = U S^{-1/2}$ (directions with small $S$ deflated) and
$G_{ij} = \langle\langle\rho_i|G|\rho_j\rangle\rangle$, the small matrix
$G_\mathrm{eff} = X^\dagger G X$ is diagonalized. An eigenvector is the list of
coefficients $XR$ over the stored MPS, and it is never summed into one MPS.
Every overlap and expectation value is read from matrix elements between
stored states. The left eigenvectors come the same way from a run under
$G^\dagger$, with the identity state $|1\rangle\rangle$ added, since it is the left
steady state.

Two properties make this better than what one gets from either piece alone.
The Krylov basis is spread over many moderate-$\chi$ states, so the effective
bond dimension is roughly $N_K \max_i \chi_i$. And Appendix G shows that the
steady state converges at the rate set by $\mathrm{Re}\,\lambda_3$, not by the
gap.

dmrgpy has none of the open-system part today, not even in ED. It does have
most of the building blocks: TDVP, MPS overlaps, batched bras, a doubled chain
(`thermal.py`), charge grading (`pyitensor/sector.py`) and non-Hermitian DMRG
and Arnoldi to compare against.

## 2. Phase 0: non-Hermitian real-time evolution on "python"

Everything else evolves MPS under a non-Hermitian generator, so this must be
correct first. It is not today.

**What was measured.** We take a 4-site spin-1/2 chain. The start state is the
ground state of
$h_0 = 0.3 S^z_0 - 0.2 S^z_3 + 0.4 S^x_1 + 0.5 S^z_2 + 0.1 S^x_2$
(nondegenerate). It evolves under the Heisenberg chain plus
$0.5 S^x_0 + 0.7i\,S^z_0$, with $dt=0.05$ for 40 steps, and we measure
$\langle S^z_2\rangle$. The exact reference is a dense `expm` in numpy,
normalized and raw.

```
t index           0       8       16      24      32
exact, normed   -0.4903 -0.4294 -0.2730 -0.0860  0.0649
v3 TDVP         -0.4903 -0.4294 -0.2730 -0.0860  0.0649
"python" TDVP   -0.4903 -0.4300 -0.2752 -0.0896  0.0621
exact, raw      -0.4903 -0.3253 -0.1595 -0.0407  0.0271
mode="ED"       -0.4903 -0.3253 -0.1595 -0.0407  0.0271
```

Without the $0.7i$ term, "python" and v3 agree with ED to 6e-7. So there are
two separate problems.

- **"python" is wrong by up to 4e-3.** v3 reproduces the exact normalized
  curve with the same chain, the same time step and the same TDVP, so the
  difference is in pyitensor. The suspect is `tdvp._lanczos_expm_multiply`,
  whose docstring says it is for a Hermitian operator: the three-term
  recursion is wrong for a non-Hermitian matvec, and nothing checks for that.
- **ED and DMRG measure different quantities.** Every DMRG
  `evolve_and_measure_*` loop restores `norm0` after each step (on purpose,
  see CLAUDE.md), so it gives $\langle\psi|O|\psi\rangle/\langle\psi|\psi\rangle$.
  ED's `evolution_ABC` gives the raw $\langle\psi|O|\psi\rangle$. For a
  Hermitian $H$ the two are the same. For a non-Hermitian one they differ by
  0.12 here.

**Steps.**

1. Confirm the cause. Replace the local exponential by a dense `scipy.linalg.expm`
   of the small effective matrix, on this same chain. If "python" then agrees
   with v3 to 1e-6, the Lanczos exponentiator is the whole cause.
2. Add `_arnoldi_expm_multiply` next to the Lanczos one in `pyitensor/tdvp.py`:
   full Gram-Schmidt, a Hessenberg matrix, its exponential, and the same Saad
   a posteriori estimate written with $h_{k+1,k}$ in place of $\beta_k$,
   relative to the vector and scale-free in `coeff`, as the Lanczos one is
   now. It should share `_krylov_verdict` and the sub-stepping, so the two
   paths cannot drift apart. The Hessenberg exponential is small and runs on
   the host; the matvecs and the basis go through `backend.xp`.
3. Choose between the two once per call: Lanczos when the generator is proven
   or witnessed Hermitian (`Many_Body_Chain.is_hermitian`, cached), Arnoldi
   otherwise. The Hermitian path must stay byte-identical. The existing
   `tests/test_pyitensor_time_evolution_review.py` is the check for that.
4. Allow a complex `dt` in `tdvp_step`, since Phase 1 needs
   $e^{-H\,\delta\tau\,e^{-i\alpha}}$. Check this against a dense `expm` on the
   same 4-site chain.
5. Make ED and DMRG answer the same quantity. Which one is a decision for the
   author (section 8). Whatever is chosen goes into `docs/numerical_changes.md`,
   since one of the two modes changes its numbers for every non-Hermitian $h$.
6. Record the finding with its repro in an audit record, in the shape the
   `hole-hunt` skill fixes, and pin it in a test:
   `tests/test_nh_time_evolution.py`, with the table above as the reference,
   on "python" and v3, at `abs=1e-6`.

The `julia_live` backend was not checked. It should be checked with the same
script before Phase 0 is closed, and it should raise on a non-Hermitian $h$ if
it is also wrong.

**Checkpoint, 2026-10-29.** Step 1 is done. With `tdvp._lanczos_expm_multiply`
replaced by a dense `expm` of the local matrix (built by applying `matvec` to
unit vectors), "python" gives -0.429444, -0.272961, -0.085991, 0.064909 at the
table's times, which is the exact normalized curve. So the Lanczos
exponentiator is the whole cause, and the sweep structure and the environments
are sound. The repro scripts are in `docs/clik_mps_repro/`:
`nh_tdvp2.py <gamma>` compares "python" DMRG with ED (`nh_tdvp_v3`-style runs
change `itensor_version`), `ref.py` is the dense numpy reference, and
`nh_dense.py lanczos|dense` is the confirmation above. Run them with
`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -I`.

The fix, as agreed, in this order:

1. Read `applyExp` in `mpscpp3/ITensor/itensor/iterativesolvers.h` and use it
   as the template, since v3 is the backend that is right. Check whether it
   builds a full Hessenberg, which error estimate it uses, and whether it
   sub-steps. (Done. It is a three-term Lanczos with a real alpha and no
   reorthogonalization. That keeps $A V = V T + f e_k^T$ exact with an
   oblique basis, and this is why it is right. pyitensor's full
   reorthogonalization removes components that $T$ never records, so the
   fix is Arnoldi, as below.)
2. Write `_arnoldi_expm_multiply` as a separate function in
   `pyitensor/tdvp.py`. Leave `_lanczos_expm_multiply`, its host and device
   paths and `_krylov_verdict` untouched (the locked rule that the two paths
   share one test), and keep its name, since `idmrg_window.py` imports it.
   The new function takes the full column `h = Q^H w`, reorthogonalizes once
   more, uses `h_{k+1,k} = ||w||`, takes the first column of `expm` of the
   complex k x k Hessenberg on the host, and uses Saad's estimate with
   `h_{k+1,k}` in place of `beta_k`, with the same relative exhaustion test
   and the same dyadic sub-stepping. On a device it syncs `h` once per
   iteration; the docstring says so, and it is checked under
   `backend.set_backend("jax")`. Note that the Lanczos code drops the
   non-Hermitian information in two places: `vdot(q, w).real`, and the
   projections beyond the tridiagonal.
3. Dispatch on the operator, never on `coeff` (METTS passes an imaginary `dt`
   with a Hermitian H and must stay on Lanczos). `tdvp_step(...,
   hermitian=None)`: `True`/`False` for a caller that knows, `None` for an
   exact check of the MPO once per call, threaded down to the `_evolve_*`
   helpers. No cache on the MPO object, since `_Chain.__imul__` mutates in
   place. The check is `||H - H^dag||_F / ||H||_F <= 1e-12`, with
   `Chain.hermitian_operator` (chain.py) moved to `mpsalgebra`. First verify
   whether `_Chain.position` handles an MPO's two Site legs. If it does,
   canonicalize `D = H - H^dag` and read the norm off the center, which has
   no cancellation. If not, `2 Tr(H^dag H) - 2 Re Tr(HH)` resolves the ratio
   only to about 1e-8, and the docstring must say so.
4. Do not touch `_renormalized`: v3 also renormalizes and matches the exact
   normalized curve. The ED mismatch stays decision 1 in section 8.
5. Tests: `tests/test_nh_time_evolution.py` pins the table against the dense
   reference at `abs=1e-6`, on "python" and v3. Then
   `tests/test_pyitensor_time_evolution_review.py` and the "python"
   dynamical-correlator tests show that the Hermitian path did not change.

## 3. Phase 1: the complex-contour Krylov eigensolver

This phase does not need a Lindbladian. It works for any generator, and the
first one we test it on is $G=-H$ for the non-Hermitian chains dmrgpy already
has. Then the target is the smallest $\mathrm{Re}\,E$, the same target as
NH-DMRG, and ED and NH-DMRG are both references. For a Hermitian $H$ it
reduces to a Krylov space built from imaginary-time states. That is a cheap
check of the projection step, but not a reason to use the method.

**Algorithm** (`pyitensor/clik.py`, new).

1. Probe $\alpha_{\max}$ (the paper's Algorithm 3). Take a random state with
   no component on the slowest mode (traceless for a Lindbladian; for $-H$,
   orthogonal to a first NH-DMRG estimate). Run a few short evolutions at a
   grid of angles. The first angle at which the norm grows bounds
   $\alpha_{\max}$. The default is a fixed fraction of it, `alpha_frac=0.25`.
2. Evolve the two branches with Phase 0's `tdvp_step`, real-time from
   $\rho_p$ and complex-time from $\rho_u$, without the per-step norm
   restore. Each stored MPS is normalized and its scale is kept as a separate
   number, since the scale does not change the span and would otherwise
   overflow or underflow.
3. After each new state, add its row and column of $M$ and $G_{ij}$, the
   latter as $\langle\rho_i|G|\rho_j\rangle$ through one MPO application per new
   state. The overlaps go through `mpsalgebra.BatchedBras`, the batched path
   that already runs on the device. The cost stays linear in the number of
   states per step.
4. Solve: diagonalize $M$, drop the directions below a threshold relative to
   $\max S$ (1e-14 in the paper; relative to the largest eigenvalue here,
   following the scale rule in CLAUDE.md), form $G_\mathrm{eff}$, diagonalize
   it on the host, and sort by real part.
5. Certify each eigenpair with its residual,
   $\lVert G r - \lambda r\rVert^2 = c^\dagger(W - \lambda^* G_\mathrm{eff} - \lambda G_\mathrm{eff}^\dagger + |\lambda|^2)c$
   in the orthonormal basis, with $W_{ij} = \langle G\rho_i|G\rho_j\rangle$.
   $G\rho_i$ is applied with the `maxm` cap lifted, as `vev.energy_variance`
   does, and never as an MPO of $G^\dagger G$.
6. Left eigenvectors: the same steps under $G^\dagger$ (Phase 1b, after the
   right ones work).

**Public call.** One method on `Many_Body_Chain`. It returns the eigenvalues,
their residuals and an object that gives expectation values in each
eigenvector. It is "python" only at first; any other backend raises, and
`mode="ED"` answers with a dense `eig` of $G$. The name and signature are a
decision for the author (section 8).

**Tests** (`tests/test_clik_spectrum.py`):

- On a Hermitian Heisenberg chain with 6 sites, the lowest four energies agree
  with ED to 1e-8.
- On the PT-symmetric chain and the asymmetric-hopping chain of
  `tests/test_nh_dmrg.py`, the lowest $\mathrm{Re}\,E$ agrees with ED and with
  NH-DMRG.
- A complex-conjugate pair that `alpha=0` cannot separate within a fixed number
  of states is resolved at `alpha>0`.
- The residual of step 5 agrees with the residual computed from the dense
  vector, on ED-sized chains.
- The same results under `backend.set_backend("jax")`.

## 4. Phase 2: the vectorized Lindbladian

**Physics.** We take $\dot\rho = -i[H,\rho] + \sum_\mu \left(L_\mu\rho L_\mu^\dagger - \tfrac12\{L_\mu^\dagger L_\mu,\rho\}\right)$.
With $|\rho\rangle\rangle = \sum_{ab}\rho_{ab}|a\rangle\otimes|b\rangle$, the
generator is

$$\mathcal{L} = -i\,(H\otimes 1 - 1\otimes H^T) + \sum_\mu \left( L_\mu\otimes L_\mu^* - \tfrac12 L_\mu^\dagger L_\mu\otimes 1 - \tfrac12\, 1\otimes (L_\mu^\dagger L_\mu)^T\right).$$

The physical sites and the auxiliary sites alternate (site $i$ goes to $2i$
and $2i+1$). This keeps a nearest-neighbor $H$ next-nearest-neighbor in
$\mathcal{L}$, as in the paper. The trace is
$\mathrm{Tr}\,\rho = \langle\langle 1|\rho\rangle\rangle$, and an observable is
$\mathrm{Tr}(O\rho) = \langle\langle 1|(O\otimes 1)|\rho\rangle\rangle$.

**Steps.**

1. The transpose of an operator. `multioperator.get_dagger` takes a
   `conjugate=` argument and does not use it, so dmrgpy has no transpose today.
   Add `MultiOperator.get_transpose()` with an explicit table per name:
   `A`↔`Adag`, `Sp`↔`Sm`, `Sx`, `Sz`, `N` and the projectors symmetric, `Sy`→`-Sy`.
   Coefficients are not conjugated. A name not in the table raises (CLAUDE.md
   rule against a silent `else`). Then $L^* = (L^\dagger)^T$. Either remove the
   unused `conjugate=` argument after checking its callers, or make it work.
2. `lindbladtk/vectorize.py`: from $H$ and a list of jumps on an $n$-site chain,
   build $\mathcal{L}$ as a `MultiOperator` on the $2n$-site chain of the same
   site type. Spins and bosons only. Fermions raise: on a doubled chain the
   Jordan-Wigner strings need their own treatment (the superfermion
   representation), and it is not part of this plan.
3. The identity state $|1\rangle\rangle$ as an MPS, built directly in pyitensor:
   bond dimension $d$ inside each pair, 1 between pairs.
4. ED reference, independent of step 2: build $\mathcal{L}$ as a sparse matrix
   from the physical $H$ and $L_\mu$ matrices (`MO2matrix`), restricted to a
   particle-number sector when the jumps conserve it, and solve it with sparse
   `eigs`. For the paper's check (5 sites, 4 bosons) the physical sector has
   dimension 70 and the Liouville sector 4900.
5. The two sublattice charges. Number-conserving jumps keep both $N_\mathrm{phys}$
   and $N_\mathrm{aux}$, and the paper enforces both so that truncation does
   not leak into other sectors. `pyitensor/sector.py` grades the total charge
   only, so this needs a per-sublattice grading. It comes after the dense
   version works and is measured; the leak is then something we can see.
6. A public entry point: either a new `Lindblad_Chain` class or a keyword of
   the Phase 1 method that takes jumps. This is a decision for the author
   (section 8).

**Tests**, simplest first:

- A single spin with decay: $\lambda$ and the steady state are known in closed
  form.
- An XXZ chain with 4 sites, $S^z$ dephasing and boundary driving: the full
  Liouville space is 256, and every eigenvalue we resolve agrees with ED.
- The paper's check: Bose-Hubbard with dissipation-assisted hopping,
  $L_j = \sqrt{\kappa}\,(b^\dagger_{j+1} + b^\dagger_j)(b_{j+1} - b_j)$, at 5
  sites, 4 bosons and $J=U=\kappa=1$, gives $\lambda_1 = -0.84+0.81i$ in both
  ED and Krylov.
- $\mathrm{Tr}\,\rho_\mathrm{ss} = 1$ after normalization, and
  $\langle\langle 1|\mathcal{L} = 0$ to machine precision, as a check of the
  vectorization itself.

## 5. Phase 3: steady state, observables and the error bound

- **The steady state.** It is the right eigenvector with $\lambda \approx 0$,
  normalized by its trace. Expectation values are $\mathrm{Tr}(O\rho)/\mathrm{Tr}\,\rho$
  through the matrix elements $\langle\langle 1|(O\otimes1)|\rho_i\rangle\rangle$.
  The Hermiticity projection $(\rho+\rho^\dagger)/2$ of the paper comes later,
  once a non-Hermitian artifact is actually seen.
- **The error bound (Eq. 40).** It reads
  $|\langle\langle\rho_K|\rho_\mathrm{ss}\rangle\rangle|^2 \ge 1 - \lVert\mathcal{L}\rho_K\rVert^2/\sigma_2^2$,
  with $\sigma_2$ the smallest nonzero singular value of $\mathcal{L}$. The
  numerator is Phase 1's residual. $\sigma_2$ is known only where ED reaches.
  Beyond that we report the residual itself, and the bound only where
  $\sigma_2$ is known; we do not invent an estimate of it.
- **The gap and its scaling.** $\Delta(L) = -\mathrm{Re}\,\lambda_2$ is read
  from the spectrum directly. A finite-size scaling example is the natural
  demonstration, but at the sizes of the paper (40 sites, $\chi\sim 1400$) it
  is a Triton job, not a test.

## 6. Left for later

- **The paper's warmup of the initial states (Algorithm 6).** It is made for
  the Bose-Hubbard model, and the paper reports that its effect is subleading
  next to combining the real and complex branches. Random starts first.
- **LSE-TDVP.** One-site TDVP with a local subspace expansion is about $d$
  times cheaper than two-site TDVP. `TDVP_GSE` covers the bond growth for now.
- **v3 and `julia_live`.** v3's exponentiator is already correct on a
  non-Hermitian $h$ (Phase 0's table), so a v3 port of Phase 1 is mostly
  bookkeeping. It comes after "python" is validated.
- **Fermions, the Mpemba overlaps, and infinite chains.**

## 7. Documentation and examples

When a phase lands: a user guide page (`docs-sync` skill) on the open-system
spectrum and steady state, written like the "Ground state" page;
`docs/development.md` for the new dispatch; and two examples (`new-example`
skill, 40 lines, with a plot): the Krylov spectrum of a non-Hermitian chain
against ED, and the steady-state magnetization profile of a driven XXZ chain
against ED.

## 8. Decisions for the author

1. The quantity a non-Hermitian `evolve_and_measure` returns on both modes:
   normalized, raw, or a keyword with one default.
2. The name and signature of the Phase 1 method.
3. A `Lindblad_Chain` class, or jumps as a keyword of the Phase 1 method.
