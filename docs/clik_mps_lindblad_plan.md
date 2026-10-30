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
the `julia_live` check, and an audit record.

**Phase 1 implemented (2026-10-29)** as `Many_Body_Chain.get_krylov_spectrum`
(`pyitensor/clik.py`, `krylovspectrum.py`, `tests/test_krylov_spectrum.py`,
user guide "Excited states"). Measured at the defaults (`nk=20`, `dt=0.5`,
`gram_tol=1e-12`) on 6-site chains against ED: Heisenberg ground state
1e-11, first triplet 1e-8, and the lowest four levels of the three
non-Hermitian test chains 3e-7 to 2.6e-4. The 20-site Heisenberg ground
state is 1.1e-5 from DMRG. Three findings change the plan. First,
`dt=0.1` is too small: consecutive snapshots are nearly parallel and the
lowest levels stay 1e-3 off, which is why the default is 0.5. Second, at
equal Krylov size the contour angle gave no measurable gain on these
Hamiltonian chains, so the test planned for it was not pinned; its
benefit is expected for Lindbladian modes and is checked in Phase 2.
Third, two random starts resolve at most two copies of a degenerate
level. Left eigenvectors (Phase 1b) are not done.

**Phase 2 implemented (2026-10-29)** as `Lindblad_Chain` (`lindblad.py`,
`pyitensor/vectorize.py`, `tests/test_lindblad.py`, user guide "Open
systems"). This is the plan's decision 3, a class; its doubled chain is the
same model class as the physical one. Checked as follows. The decaying spin
gives the closed form to 1e-8. `Tr(L^k)` of the doubled-chain operator
matches the `kron` superoperator to 1e-16 for spin-1/2, spin-1 and bosons
with dissipative hopping. `<<1|L|rho>>` vanishes to 1e-17. On a driven,
dephased 4-site XXZ chain the steady state is 1e-6 off, the gap mode 2e-4
and the magnetization profile 3e-6. On 6 sites at the chain's default
`maxm` the gap is 0.0924 against 0.0926 and the profile is 3e-4 off. The
Lindbladian needs much longer evolutions than a Hamiltonian, since `nk*dt`
must be several times 1/gap, so `get_spectrum` defaults to `dt=2.0`.
Again no measurable gain from the angle at these sizes. Two things are not
done. First, per-sublattice sector seeding, needed for the paper's 5-site
Bose-Hubbard check, since number conservation gives one steady state per
sector. Second, the Hermiticity projection of the steady state.

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

## Status 2026-10-30

Committed so far: Phase 0 (Arnoldi exponentiator for non-Hermitian TDVP),
Phase 1 (`get_krylov_spectrum`), Phase 2 (`Lindblad_Chain`), the sector
seeds `start=(psi, phi)` with a global subspace expansion before the first
three steps of each branch, and the paper's bookkeeping. The bookkeeping
means two things: matrix elements grow one snapshot at a time, and `res` is
the change of each eigenvalue when the last snapshot of each branch is
dropped, with the exact residual opt-in (`exact_residuals=True`). The full
test suite has not been run since Phase 0. The user asked for it once, at
the very end, so it is the first thing to do.

**What the paper does about cost.** I read the paper for this, and these
are its statements. Matrix elements are "linear in the number of time
evolved states" per step and "can be calculated in parallel".
⟨⟨ρ|L†L|ρ⟩⟩ is computed "directly inside the Krylov subspace", but only
for the steady-state bound of Appendix G, never for every pair or every
mode. Convergence is judged by tracking the spectrum, with "no automated
escape criteria". The large runs use 340 to 400 states at χ up to 2400 in
SyTen. No CLIK-MPS runtimes or MPO bond dimensions are given.

**The tracking estimate is not an error bar.** On 6-site chains against
ED it was 3 to 40 times smaller than the true error of the eigenvalue. For
example, the fourth Heisenberg level was 1.1e-5 off while it moved by only
2.6e-7. It ranks the modes by convergence correctly (`tr.py`).

**The paper's 5-site Bose-Hubbard check was not finished.** The target is
5 sites, 4 bosons, J=U=κ=1, λ1 = -0.84+0.81i, run as `bh.py 64 20 0.5 ed`.
It was stopped by the user after about 25 minutes, still in the
evolution. The seeded 3-site version (`bh_small.py`) matches ED and is
pinned in `tests/test_lindblad.py`.

**Where the time goes (profiling subagent, `tdvp_profile/`).** The 5-site
case is a doubled chain of 10 sites with local dimension 5. The MPO of L
has bond dimensions 1, 8, 22, 14, 22, 14, 22, 14, 22, 8, and an exact SVD
sweep does not shrink it. One two-site TDVP step at χ=64 takes about 56 s,
so 2 branches × 19 steps is about 35 minutes of evolution. Of the step
time, 74% is the local matvec, which already runs at the one-thread BLAS
peak (about 35 GFlop/s), and about 17% is Arnoldi bookkeeping, mostly
copying the conjugated basis every iteration. The SVDs, the environments
and the expm are negligible. Options it measured or estimated:

- One-site TDVP. Its solves cost 23 ms per iteration against 142 ms for
  two-site, so it estimated 2 to 3 times faster steps, without the cost of
  the expansion.
- A row-major Krylov basis: 42 → 16 ms per iteration in a microbenchmark,
  about 12% overall by estimate. It was not verified end to end, because
  machine load spoiled that run.
- `errgoal` 1e-8 or 1e-6: 33 → 29 or 25 matvecs, against the documented
  1e-10 rule.
- Two things did not help: selective reorthogonalization (the second pass
  fired 97% of the time) and MPO sparsity.

The Hermitian case was not profiled.

**One-site TDVP with expansion leaks out of the sector**
(`clik._step(one_site=True)`, experimental, off everywhere). On the seeded
3-site Bose-Hubbard Lindbladian, the in-sector Schmidt rank at the middle
cut is exactly 14; counting the (N_ket, N_bra) blocks across the cut gives
that number. Two-site TDVP stays at 14. With an expansion before every
one-site step, the middle bond grows to 27, and the traceless branch
drifts from the two-site one: 1 - overlap = 2.4e-3, 2.1e-2, 0.14, 0.37 and
0.62 over five steps of dt=0.5 (`os2.py`, `os4.py`). The steady state is
still right, to 1e-6, but the slow modes are off by order 1. Two things
are measured:

- One-site TDVP without the expansion, started from the full-rank in-sector
  snapshot, agrees with two-site TDVP to integrator error: 1.3e-5, 1.7e-9
  and 1.8e-13 at dt = 0.5, 0.1 and 0.02. The bond stays at 14 (`os3.py`).
  So the integrator is fine.
- Raising the expansion cutoff from 1e-10 to 1e-4 does not stop the leak
  (`os4.py`). So it does not come from tiny-weight bond directions, and its
  source is open.

The paper enforces both sublattice U(1) symmetries "to keep truncation
errors from leaking into unphysical sectors" (Appendix H). In pyitensor a
sector grading does not confine: storage is dense, link indices carry no
charges, and `sector.py` confines variational solves only with a charge
penalty. The candidates are:

- a penalty on the generator, G' = G - λ[(N_ket-N)² + (N_bra-N)²], which
  equals G on the target sector;
- link-index charges in the tensor core, which is a rewrite;
- dropping the one-site route.

Whether one-site is faster at all was never measured on the 5-site case
(`os_time.py` was stopped before its first step). On 4 sites it was
slower, 14.3 s against 13.3 s.

**Open, in order.**

1. Run the full test suite.
2. Decide whether one-site TDVP is worth pursuing. Time it first with
   `os_time.py` on a quiet machine, and only then choose among the three
   candidates above.
3. Finish the 5-site check, possibly as a Slurm job.
4. Left eigenvectors (Phase 1b), the Hermiticity projection of the steady
   state, and the Phase 0 leftovers: the ED-against-DMRG norm convention,
   the `julia_live` check and the audit record.

## Night of 2026-10-30

The task was to debug and optimize the Lindblad route, with questions put
to a second model rather than to the author. The repro scripts are in
`docs/clik_mps_repro/night_2026_10_30/`, with a README saying what each
one measured; the numbers below are what they printed.

**Where the time goes depends on the size.** On the 4-site XXZ chain the
matrix elements were 65% of the run: 2729 `inner` calls, each a few small
contractions dominated by Python overhead. On the 5-site Bose-Hubbard case
the two-site matvec dominates. Its four contractions cost about the same,
roughly 150M multiply-adds each at χ=64, and the MPO is dense in its
virtual blocks (75% of the blocks are nonzero), so there is no exact
sparsity to use.

**Committed.**

- Batched matrix elements (`BatchedBras.sandwiches`): the 3-site
  Bose-Hubbard spectrum goes from 10.9 s to 2.8 s with the same
  eigenvalues. On a device the batch is padded to a power of two, since
  JAX compiles each new shape (the JAX test went 80 s, then 125 s, then
  71 s).
- A traceless second start for random starts (`deflate=`): the XXZ gap
  mode goes from 2.2e-4 to 8.9e-5.
- No probe for a Hermitian generator, and a bisected probe otherwise: the
  3-site Bose-Hubbard probe goes from 56 evolution steps to 24.
- A row-major Arnoldi basis on the host: 10.8 s to 4.4 s of bookkeeping
  per 5-site step, about 7% of it.

**The second branch does not stay traceless.** The exact evolution
conserves the trace, but two-site TDVP at dt=0.5 does not. On the seeded
3-site Bose-Hubbard case the in-sector trace of the normalized snapshots
grows to 1.6e-3 after one step and 5.9e-2 after five, against 1e-8 to
1e-6 for a 20 times finer reference, and the last snapshots of that branch
are about 80% steady state. Projecting the sector identity $1_N$ out after
every step keeps the trace at 1e-15. At 3 sites it moves the eigenvalues
little (one mode from 1.1e-6 to 4.4e-7, another from 5.3e-4 to 9.5e-4),
so it is not committed. Projecting the full identity instead is wrong for
a seeded run: $|1\rangle\rangle$ is not in the sector, and the occupations
went from 5e-8 to 4e-6 off ED. A general version needs the sector identity
as an MPS, a charge-counting MPS of bond dimension N+1, built from the
seed's charge.

**One-site TDVP is settled: drop it.** The leak is the QR of the one-site
sweep, which completes a rank-deficient bond basis with arbitrary vectors
that mix charges; a rank-revealing split keeps the charge variance at
1e-15. Even then the one-site branch is 0.70 off the fine reference after
five steps, against 1.6e-3 for two-site, and it is not faster on 5 sites.
A hybrid, two-site until the bonds stop growing and then one-site without
expansion, is identical on XXZ (bonds at the Hilbert-space limit) and ten
times worse on Bose-Hubbard (bonds at the in-sector rank, where one-site
keeps a projection error). Neither is committed.

**The paper's 5-site check failed at dt=0.5, and why.** Run to the end
for the first time (4825 s at χ=64, one TDVP step per snapshot), it gave
a top eigenvalue of +0.485, no steady state, and occupations that are not
even mirror symmetric; ED gives $\lambda_1=-0.8409\pm0.8135i$, the
paper's value. Two causes, both in the integrator and not in CLIK itself.
First, the Arnoldi stopping test was relative to the input vector, and a
dissipative local step can shrink it to 5e-19 of the input; it is now
relative to the result. Second, and this one is structural, two-site
TDVP evolves each bond backward by $e^{-\tau G_\mathrm{eff}}$, which on
directions decaying at rate $R\approx100$ amplifies every truncation and
Krylov error by up to $e^{\tau R}\approx e^{25}$. One step of 0.5 grew the
norm by $e^{16.8}$ at the second snapshot and left the particle sector
($\langle(N-q)^2\rangle=100$); five steps of 0.1 kept the log-norm below
0.2 and the charge variance at $10^{-6}$. Three and four sites work at
dt=0.5 because the rates grow with the number of bosons. A second model
(Fable) agreed with this reading: the stopping test cannot cure it, and
the principled fix is a forward-only integrator, the basis-update and
Galerkin (BUG) scheme, built for exactly this failure of the backward
step. That is a new integrator and is left for later. What is committed
is sub-stepping: `substeps` TDVP steps per snapshot, doubled while a step
of the real branch grows the norm by more than $e$, overflows, or moves
the trace by more than $10^{-2}$ of $\lVert 1_q\rVert$. On 3 sites with
ten times the dissipation it settles on 8 and takes the occupations from
0.86 to 3e-4 off ED at the same cost.

A rejected step restarts the whole branch at the doubled count, since at
5 sites the first step at one sub-step passes the gate (the trace off by
6e-3 of $\lVert 1_q\rVert$, the state already 3e-4 out of the sector)
and only the second fails. The trace test cannot be made tighter instead:
on the guide's 6-site XXZ chain, where nothing is unstable, truncation
moves the trace by up to 1.2e-3 of $\lVert 1\rVert$, and a stable
sub-stepped 5-site step by 1.3e-3. Norm growth and overflow are checked
on both branches whenever a left steady state is known.

**The 5-site check with sub-steps.** Run with five forced TDVP steps
of 0.1 per snapshot (code at `46fb3bd`: relative Arnoldi test and
per-step projection of $1_q$, no gate yet), nk=20, dt=0.5, alpha=0.2,
χ=64, 3955 s on a loaded machine. The steady state is right: the
occupations are 0.7114, 0.8492, 0.8789, 0.8494, 0.7117 against 0.7123,
0.8486, 0.8781, 0.8486, 0.7123 from ED, mirror symmetric to 3e-4 and
about 1e-3 off. Every one of the slowest eight ED eigenvalues, the
paper's $\lambda_1=-0.8409\pm0.8135i$ included, has a Krylov eigenvalue
within 0.05 of it (the script printed the distance to one digit only).
But the list is polluted. The steady state comes out at $-0.0102$, and
above $\lambda_1$ sit modes that ED does not have, such as
$-0.107+2.32i$ and $-0.239$, with residuals of 6e-3 to 0.4. On 4 sites
nothing like this happens, and the charge variance of each Ritz vector,
read off the stored snapshots, tracks its convergence: 5e-13 to 2e-7 on
the modes that match ED to 1e-5 or better, 1e-3 to 1e-2 on the
unconverged ones. The candidate explanation for 5 sites is the charge
leak of the snapshots (about 1e-6): modes of other sectors that decay
slowly enter the Krylov space. The same diagnostic is running on 5
sites; see the end of this section.

**The 5-site check with the committed default.** Run with
`substeps=1`, so that the gate has to find the count (code at
`47fb3e2`), nk=20, dt=0.5, alpha=0.2, χ=64: the gate settled on 2
sub-steps in both branches, and the run took 2249 s, against 3955 s
with five forced sub-steps. This reproduces the paper's check:
$\lambda_1=-0.8415+0.8139i$ against $-0.8409+0.8135i$ from ED, 7e-4 off
with a residual of 2.6e-4, and its partner 8e-3 off. The occupations are
0.7116, 0.8489, 0.8784, 0.8484, 0.7114, 1e-3 off ED. The same pollution
as above is there: the steady state at $-0.0057$, and $-0.094-2.32i$,
$-0.178$, $-0.243-3.78i$ and $-0.50+2.23i$ above $\lambda_1$. Their
imaginary parts, of the order of U times a change of particle number,
point to coherences between particle sectors, which is what the
diagnostic below checks.

**The spurious 5-site modes are other sectors, and are now dropped.**
The charge variance $\langle(Q_\mathrm{ket}-q)^2+(Q_\mathrm{bra}-q)^2\rangle$
of each Ritz vector, read off the stored snapshots, is 2.4 to 8.8 for
every spurious mode listed above and 2e-5 to 2e-3 for every mode that
matches ED ($\lambda_1$: 2.1e-4). A mode of another sector has at least
1, since the charges are integers. A seeded run whose sector identity
checks out now drops the modes above 0.5. On the 5-site run the list
that remains is the steady state ($-0.0057$), $\lambda_1$ (7e-4 off),
its partner (8e-3), $-1.328$ (2e-2 off $-1.348$), $-1.453+1.450i$
(4e-3) and $-1.492-1.477i$ (5e-2): the ED spectrum in order. The steady
state eigenvalue, 5.7e-3 from zero, is the largest remaining error; the
occupations are 1e-3 off.

**A looser Krylov error goal is not a lever.** Earlier tonight
`errgoal=1e-6` looked 2.3 times faster on 5 sites; that was the broken
stopping test, and its result was garbage. With the test relative to the
result, `errgoal=1e-6` on the stiff 3-site case gives occupations 0.26 off
ED, against 3e-4 at the default 1e-10. The goal stays where it is.

**The contour angle may matter more than the probe assumes.** On the
seeded 3-site Bose-Hubbard case at dt=0.5, the probe picks 0.39 and the
slowest pair is 2.6e-4 off ED; at a fixed angle of 0.2 it is 2.9e-6 off.
One data point, recorded and not acted on.

**The Davidson commit broke five tests**, all outside the Lindblad code.
`b131b53` (2026-10-25) fails them and its parent passes them:
`test_audit_2026_09_24b_groundstate.py` (three Kondo crossing tests),
`test_audit_2026_09_24c_session.py::test_sector_makes_no_solve_on_the_callers_chain[python]`
and `test_audit_2026_09_25_scale.py::test_small_units_ground_state_still_open[python-1e-10]`,
a strict xfail that now passes.

**Open, in order, after the night of 2026-10-30.**

1. A forward-only integrator, the basis-update and Galerkin (BUG) scheme,
   as the principled cure for the stiffness; sub-stepping only works
   around it. The augmented basis must be split with an SVD, not a QR, or
   it leaks out of the sector as one-site TDVP did.
2. The steady-state eigenvalue of the 5-site check is 5.7e-3 from zero at
   χ=64. The next check is χ=96 or 128, as a Triton job, not a local run.
3. The growth gate may reject a legitimate step. `_GROWTH_TOL` is 1 while
   the norm of a normalized density matrix can grow by up to
   $\tfrac12\ln D$ as it purifies (2.1 at D=70). A fast-purifying
   Lindbladian at a long dt could then be doubled up to 64 sub-steps and
   end in "reduce dt" for nothing. No case tonight hit it.
4. The contour angle: on 3 and 4 sites a fixed 0.2 gave the slowest pair
   2.9e-6 and 9.5e-6 off ED, against 2.6e-4 and 4.5e-4 at the probed 0.39.
   Two data points, not acted on.
5. The five test failures of the Davidson commit `b131b53`, listed above.
6. The two examples of section 7, which do not exist yet.
7. The Phase 0 leftovers: the ED-against-DMRG norm convention, the
   `julia_live` check and the audit record. Left eigenvectors (Phase 1b)
   and the Hermiticity projection of the steady state.

One-site TDVP is settled: it is dropped.
