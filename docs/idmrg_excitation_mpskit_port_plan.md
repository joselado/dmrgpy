# Plan: reimplement iDMRG D>1 excitations by mirroring MPSKit.jl's own architecture

## Status and decision

After eight investigation passes (see `src/dmrgpy/pyitensor/idmrg_excitations.py`'s
own module docstring, and `ROADMAP.md` item 4 under "iDMRG") repeatedly produced
"internally self-consistent but wrong" results for a *hand-derived, patched*
mixed-gauge tangent-space excitation ansatz, the decision (2026-08-07) is to
**stop patching the existing diagram-list construction and instead write a
whole new implementation that mirrors MPSKit.jl's actual algorithm
architecture** -- not just its final equations (which is all passes 4/5/7 ever
had access to, via paper images), but its real, working source code, which is
now available locally on this machine (see "Environment" below).

This file is the handoff: everything needed to pick this up in a fresh
session without re-deriving what this session already found.

## Why this, not another patch attempt

- Eight passes across two gauge formulations (uniform, then mixed) and four
  structurally different diagram sets all failed the same qualitative way
  (suppressed dispersion, non-Hermitian H_eff at D>1) while being D=1-exact
  every time -- strong evidence of a systematic derivation/convention error
  recurring across independent hand-derivations, not bad luck on any single
  attempt.
- **New this session: an actual numerical oracle now exists and confirms the
  physics is not the problem.** MPSKit.jl's own `QuasiparticleAnsatz` gives a
  D=2 TFIM (J=1, h=2.5) dispersion matching the exact free-fermion band to
  5-6 significant figures at every momentum tested (0, 0.5, ..., 3.0) -- see
  "Oracle run" below for the exact numbers. This *proves* the tangent-space
  ansatz works correctly at D>1 in a correct implementation; the bug is
  entirely dmrgpy's own from-scratch port, not some subtlety that makes the
  method itself unreliable at D>1.
- **New this session: a concrete, checkable discrepancy was found between
  the sixth investigation pass's own hand-derived mixed-transfer fixed
  points and MPSKit's actual ground-truth ones** -- see "Fixed-point
  discrepancy" below. This is a genuinely new, high-priority lead, distinct
  from anything the first eight passes checked.

## Environment: how to reach MPSKit.jl again

Julia is already installed on this machine via `juliapkg` (the mechanism
`juliacall`/`install_julia.py` use) -- **no download needed**:

```
JULIA=/u/40/ladovj1/unix/.julia/environments/pyjuliapkg/pyjuliapkg/install/bin/julia
$JULIA --version   # 1.12.6, confirmed working
```

(This path is specific to this machine/user account -- if working from a
different checkout or account, first check whether `python3 -c "import
juliapkg; print(juliapkg.executable())"` still resolves it; dmrgpy's own
`juliacall` dependency is already what provides `juliapkg`.)

A **separate, scratch Julia project** (deliberately NOT touching dmrgpy's own
shipped `src/dmrgpy/juliapkg.json`, which controls the actual `julia` extra's
dependencies -- see `CLAUDE.md`'s packaging section for why that file must
stay minimal) was created and has `MPSKit`, `TensorKit`, and `MPSKitModels`
already resolved and precompiled:

```
/tmp/claude-3291340/.../scratchpad/mpskit_oracle/{Project.toml,Manifest.toml}
```

That exact `/tmp/...` path is this session's ephemeral scratchpad and will
not survive to a future session. To recreate it anywhere (~2-3 minutes,
mostly precompilation, confirmed this session):

```bash
mkdir -p /path/to/scratch/mpskit_oracle && cd /path/to/scratch/mpskit_oracle
$JULIA --project=. -e 'using Pkg; Pkg.add(["MPSKit", "TensorKit", "MPSKitModels"]); Pkg.precompile()'
```

Then run any `.jl` script against it with `$JULIA --project=. script.jl`.
`juliacall` (Python) could also be pointed at this same project via
`JULIA_PROJECT` if a Python-Julia round-trip becomes worthwhile instead of
file-based data exchange -- not attempted this session, see "Convention
mapping" below for why file-based (or PythonCall array conversion) is
probably the pragmatic choice regardless.

## Oracle run (already done -- exact output for reference)

`tests/test_vumps.py::_tfim_chain`-equivalent convention: `H = -4*Sx_i*Sx_{i+1}
- 2*g*Sz_i`, i.e. Pauli `H = -sigma^x sigma^x - g*sigma^z`. MPSKitModels'
`transverse_field_ising` uses a *different* axis convention (`H =
-J*(sigma^z_i sigma^z_j + g*sigma^x_i)`) but the exact free-fermion
dispersion `eps(k) = 2*sqrt(J^2+g^2-2*J*g*cos(k))` is basis-independent, so
this is a valid apples-to-apples check regardless:

```julia
using MPSKit, MPSKitModels, TensorKit
H = transverse_field_ising(; J=1.0, g=2.5)
psi0 = InfiniteMPS([ComplexSpace(2)], [ComplexSpace(2)])   # D=2
psi, envs, delta = find_groundstate(psi0, H, VUMPS(; maxiter=400, tol=1e-10))
energies, qps = excitations(H, QuasiparticleAnsatz(), [0.0,0.5,1.0,1.5,2.0,2.5,3.0], psi, envs)
# energies is an (nmomenta, num) matrix, see MPSKit's own excitations() docstring
```

Result (`e0 = -2.601031372261`, VUMPS converged to `delta=8.0e-11`):

| p | MPSKit E(p) | exact eps(p) |
|---|---|---|
| 0.0 | 3.000000 | 3.000000 |
| 0.5 | 3.383488 | 3.383541 |
| 1.0 | 4.265419 | 4.265437 |
| 1.5 | 5.252138 | 5.252167 |
| 2.0 | 6.109245 | 6.109250 |
| 2.5 | 6.709876 | 6.709909 |
| 3.0 | 6.985686 | 6.985689 |

## MPSKit's actual architecture (read directly from source, not images)

Package source, once installed as above:
`~/.julia/packages/MPSKit/<hash>/src/` (hash was `NXsUL` this session, will
differ on a fresh `Pkg.add`; find it via
`find ~/.julia/packages/MPSKit -maxdepth 1 -type d`). Key files:

- `environments/qp_envs.jl` -- `InfiniteQPEnvironments`, the channel-resolved
  GBL/GBR construction (**this is the "channel-resolved GBL/GBR recursion"
  the seventh pass attempted from paper images and found internally
  self-consistent but wrong on its rightward half** -- MPSKit's actual
  working code is right there now, not reconstructed).
- `algorithms/excitation/quasiparticleexcitation.jl` -- `QuasiparticleAnsatz`,
  `EffectiveExcitationHamiltonian`, `_effective_excitation_local_apply` (the
  actual per-site H_eff formula), `effective_excitation_renormalization_energy`
  (the energy-density subtraction, dmrgpy's own `e_cell` equivalent).
- `states/infinitemps.jl` -- `l_RL`/`r_RL`/`l_LR`/`r_LR` (the mixed-transfer
  fixed points pass six had to hand-derive; MPSKit's own ground truth is
  right there, see "Fixed-point discrepancy" below).

### Tensor convention (confirmed by direct extraction, not assumed)

```julia
AL1 = psi.AL[1]
space(AL1)      # (ℂ^2 ⊗ ℂ^2) ← ℂ^2   i.e. codomain=(left,phys), domain=(right)
convert(Array, AL1)   # size (2,2,2)
```

`convert(Array, T)` on an MPS site tensor gives array shape `(D_left, d_phys,
D_right)` **with no extra conjugation/transpose needed** -- this matches
dmrgpy's own `(D, d_g, D)` `(l,p,r)` convention directly. This was a real,
previously-open risk (index-order mismatches are exactly the kind of thing
that silently corrupts a port) and is now confirmed cheap to get right.

### The channel-resolved GBL/GBR recursion (ground truth, `qp_envs.jl`)

For an `InfiniteMPOHamiltonian` (the common case; dmrgpy's own automaton W
with S/F/pending channels is the direct analogue -- see
`idmrg._build_periodic_mpo`):

```julia
AL = exci.left_gs.AL;  AR = exci.right_gs.AR
lBs[pos+1] = lBs[pos] * TransferMatrix(AR[pos], H[pos], AL[pos]) / cis(momentum)
           + leftenv(lenvs, pos, left_gs) * TransferMatrix(exci[pos], H[pos], AL[pos]) / cis(momentum)
# (+ a `regularize!` step using l_RL/r_RL for the trivial/momentum-0 sector, then a
#   `left_excitation_transfer_system` linsolve to close the geometric sum)

rBs[pos-1] = TransferMatrix(AL[pos], H[pos], AR[pos]) * rBs[pos] * cis(momentum)
           + TransferMatrix(exci[pos], H[pos], AR[pos]) * rightenv(renvs, pos, right_gs) * cis(momentum)
# (mirror, using l_LR/r_LR)
```

Confirmed, previously-ambiguous details this resolves (pass seven tried "4
direction/index combinations" without success, per its own writeup):
- **Phase direction**: `lBs` (left-growing, i.e. dmrgpy's own "diagram 6a"
  role) divides by `cis(momentum) = e^{+ip}` each step; `rBs` (right-growing,
  "diagram 6b" role) *multiplies* by `cis(momentum)`. These are NOT simply
  `e^{+ip}` vs `e^{-ip}` mirror images of each other in the naive way passes
  six/seven assumed -- both use the *same* `cis(momentum)`, just one divides
  and one multiplies. Check this carefully against whatever convention a new
  port picks for its own momentum sign.
- **Which mixed transfer belongs to which growth direction**: `lBs` grows
  via `TransferMatrix(AR, H, AL)` (MPSKit's own "RL" transfer, ket=AR,
  bra=AL); `rBs` grows via `TransferMatrix(AL, H, AR)` ("LR" transfer,
  ket=AL, bra=AR). Naming convention confirmed from source: `TransferMatrix`'s
  own first positional arg is the ket, third is the bra -- "XY" names mean
  ket=A_X, bra=A_Y. Cross-checked against the seventh pass's own literature
  notation resolution (`E_L^R`/`E_R^L`, bra-subscript/ket-superscript): MPSKit's
  "RL" = ket-R,bra-L matches the paper's `E_L^R` exactly (same object, just
  MPSKit's own name order is ket-then-bra rather than bra-then-ket) --
  independent confirmation the seventh pass's own notation resolution was
  right, at least for *this* piece.

### The final per-site H_eff formula (ground truth, `quasiparticleexcitation.jl`)

Startlingly simple compared to dmrgpy's own diagram-1-through-6b-plus-(i)-
plus-(ii) list -- exactly THREE terms per site, matching what the fifth pass
already suspected from MPSKit's high-level structure but never got the exact
formula for:

```julia
function _effective_excitation_local_apply(site, ϕ, H::MPOHamiltonian, E, envs)
    B = ϕ[site]
    GL = leftenv(envs.leftenvs, site, ϕ.left_gs)      # ordinary background, dmrgpy's Lh/GL
    GR = rightenv(envs.rightenvs, site, ϕ.right_gs)   # ordinary background, dmrgpy's Rh/GR
    B′ = scale(B, -E)                                  # energy-density renormalization
    # B in center:
    B′ += GL * B * H[site] * GR                        # (index contraction, see @plansor in source)
    # B to the left (uses the channel-resolved GBL for this site):
    B′ += GBL[site] * AR[site] * H[site] * GR
    # B to the right (uses the channel-resolved GBR for this site):
    B′ += GL * AL[site] * H[site] * GBR[site]
    return B′
end
```

This means MPSKit's own H_eff is **matrix-free / iterative**: GBL/GBR are
resolved (via `linsolve`/GMRES) FRESH for every trial vector inside the outer
`eigsolve` (Arnoldi) call, not precomputed once and reused the way dmrgpy's
own `_right_momentum_resolvent` is (built once per momentum, reused across
every dense basis vector). A new dmrgpy port does not have to copy this --
dmrgpy's own D is small enough in every test case so far that assembling a
small dense H_eff(p) matrix (as the current module already does) and
diagonalizing it directly remains a reasonable implementation choice; the
important thing to port faithfully is the *formula* and the *fixed-point
conventions*, not necessarily the matrix-free solver architecture.

### Fixed-point discrepancy -- the single most concrete lead for the next session

`src/dmrgpy/pyitensor/idmrg_excitations.py`'s docstring (sixth investigation
pass) derived, by hand, that the (A_L-ket, A_R-bra) mixed transfer's dominant
LEFT fixed point is `conj(C)` (verified numerically at the time, "confirmed
to ~1e-14 by direct diagonalization"). **MPSKit's own ground truth
(`infinitemps.jl`, `l_LR`) is `ψ.C[loc-1]'` -- i.e. `C^dagger` (conjugate
TRANSPOSE), not plain `conj(C)`.** These are only equal when `C` is
symmetric (`C = C^T`), which a generic gauge transformation matrix is not.
Full correspondence extracted from MPSKit source (`states/infinitemps.jl`,
functions `l_RL`/`r_RL`/`l_LR`/`r_LR`):

| MPSKit name | transfer type (ket-bra) | fixed point | pass six's own claim |
|---|---|---|---|
| `l_RL` | (A_R-ket, A_L-bra) | `C` | (pass six didn't derive this side's left fp explicitly) |
| `r_RL` | (A_R-ket, A_L-bra) | `C^dagger` (`C'`) | `C^dagger` -- **matches** |
| `l_LR` | (A_L-ket, A_R-bra) | `C^dagger` (`C'`) | `conj(C)` -- **does NOT match** |
| `r_LR` | (A_L-ket, A_R-bra) | `C` | `C` -- matches |

This is a concrete, cheaply-checkable candidate explanation for (at least
part of) the historical non-Hermiticity: a `conj(C)` used where `C^dagger`
(`= conj(C)^T`) was needed would silently corrupt exactly the kind of
D>1-only, D=1-invisible (since at D=1 every matrix is its own transpose)
regularization term this investigation has been chasing for eight passes.
**First concrete thing to try in the new implementation**: use `C^dagger`
(not `conj(C)`) for the (A_L-ket,A_R-bra) transfer's own left fixed point,
and re-run the eighth pass's own D=2 TFIM Hermiticity check (see that pass's
own diagnostic pattern in git history / the module docstring) to see if this
alone measurably improves things. This was NOT tried this session (the
decision to build a fresh implementation was made before circling back to
retrofit this single fix into the old diagram-list code, which is
deliberately being retired rather than patched again).

## Recommended shape for the new implementation

1. New module (suggest `src/dmrgpy/pyitensor/idmrg_excitations_v2.py` or an
   in-place rewrite of `idmrg_excitations.py` once the old diagram-list
   approach is confirmed superseded -- the team's call at implementation
   time) built around the *channel-resolved GBL/GBR* construction above,
   not the current module's closed-form-resummed-diagram-list approach.
   dmrgpy's own automaton `W` (S/F/pending channels, `idmrg._build_periodic_mpo`)
   is already the right kind of object to play MPSKit's `H::MPOHamiltonian`
   role -- no new Hamiltonian representation should be needed.
2. Reuse VUMPS's already-validated `AL`, `AR`, `C`, `GL`, `GR` (`vumps.py`)
   as the ordinary (momentum-independent) background -- these are NOT
   suspected of being wrong (see the eighth pass's own VUMPS-bypass
   argument in `idmrg_excitations.py`'s docstring).
3. Implement the GBL/GBR recursion literally as shown above (dense
   (D^2,D^2) linear solve is fine at dmrgpy's current scale, mirroring how
   `_regularized_environments`/`_right_momentum_resolvent` already do dense
   solves rather than MPSKit's iterative GMRES) -- but get the **phase
   direction** (divide vs. multiply by `cis(momentum)`) and the **fixed-point
   convention** (the `C^dagger` vs `conj(C)` table above) exactly right from
   the start, rather than re-discovering them by trial and error.
4. Implement the exact 3-term `_effective_excitation_local_apply` formula
   (B-in-center via ordinary GL/GR, B-to-the-left via GBL, B-to-the-right via
   GBR) in place of the current module's diagrams 1-6b/(i)/(ii).
5. Validate incrementally, the same discipline all eight prior passes used:
   D=1 exactness first (must reproduce the shipped module's own exact D=1
   dispersion, e.g. via the sanity-gate pattern in the eighth pass's own
   scratch script), THEN D=2 TFIM Hermiticity, THEN the dimerized-Heisenberg
   reproducer (second pass) and the free-fermion TFIM dispersion (this
   session's own oracle numbers above are now the reference target). Do not
   trust an end-to-end dispersion match alone -- Hermiticity of H_eff(p) is
   the cheap, structural invariant that has caught every real bug so far.
6. Once working, the "decisive experiment" fable's plan originally proposed
   (feed MPSKit's own converged AL/AR/C directly into the *new* dmrgpy
   implementation and check Hermiticity) becomes a strong independent
   validation step, now that the tensor-convention risk is already resolved
   above (identical `(D,d,D)` array layout, no transpose needed).
