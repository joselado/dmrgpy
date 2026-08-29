# Infinite chains (iDMRG / VUMPS): what's worth doing next

Written 2026-08-29, from a read of `src/dmrgpy/pyitensor/idmrg.py`,
`idmrg_excitations.py`, `vumps.py`, `src/dmrgpy/infinitechain.py` and a
survey of the current literature. `ROADMAP.md` tracks what each *backend*
has; this file tracks what the *algorithms themselves* are still missing,
which is a different axis — none of the items below is a port of something
that already exists elsewhere in the tree.

Ranked by value per unit of effort. Every anchor is a real symbol in this
checkout.

| # | Item | Effort | Status |
|---|------|--------|--------|
| 1 | Iterative eigensolver for the excitation ansatz | ~1 day | **done** (pyitensor; v3 deferred, see below) |
| 2 | Spectral weights → S(k,ω) from the excitation eigenvectors | ~1 week | open |
| 3 | Two-site VUMPS (dynamic bond dimension) | ~1–2 weeks | open |
| 4 | Long-range interactions via exponential-decay channels | ~2 weeks | open |
| 5 | iTDVP on the uniform state (global quenches) | ~2–3 weeks | open |
| 6 | QN/sector conservation for infinite chains | long-term | open |

---

## 1. Iterative eigensolver for the excitation ansatz, keeping the eigenvectors

`idmrg_excitations.py::_build_H_eff_dense` assembles the **full**
`(D²(d_g−1))²` dense `H_eff(k)` by calling `_h_eff_action` once per basis
vector, then diagonalizes it with `np.linalg.eigvalsh`. Each
`_h_eff_action` call builds `GBL`/`GBR`, i.e. linear solves through
`_solve_linear_map` (GMRES above `_DENSE_SOLVE_MAX`). So a single momentum
point costs `D²(d_g−1)` linear solves, and `Infinite_Many_Body_Chain.
excitation_gap` repeats that at 41 momenta by default.

Wrapping the *same, already-written* `_h_eff_action` in a
`scipy.sparse.linalg.LinearOperator` and calling `eigsh` needs ~30–60
actions instead of ~D², i.e. roughly 30× at D=32, growing as D². It also
removes the memory wall: at D=100, d=2 the assembled matrix is 1.6 GB
today.

Second half of the same change: `excitation_energies` currently discards
the eigenvectors `X`. They are what item 2 needs, so return them (behind an
opt-in argument, keeping the existing return shape).

### What was actually done, and the correction the measurement forced

The premise above was half wrong, and the measurement is what showed it.
`_channel_resolvent` never used `_solve_linear_map`/GMRES at all: it built
a dense `(D², D²)` matrix and called `np.linalg.solve`, unconditionally.
Worse, **it was rebuilt inside every `_h_eff_action` call** even though
nothing in it depends on the excitation tensor `B` — only on `k` and the
momentum-independent environment. So each application of H_eff paid two
dense `D²`-application builds plus two fresh `O(D⁶)` solves, `dim` times
per momentum. That, not the eigensolver, was the dominant cost.

The change therefore has three parts, not one:

1. **Per-momentum resolvent cache** (`_resolvents_for`,
   `ExcitationEnvironment.resolvent_cache`), LU-factored once via
   `scipy.linalg.lu_factor` so every later solve is two triangular solves.
   Numerically exact — no returned number moves.
2. **Iterative eigensolver** (`_lowest_iterative`): `eigsh` over a
   `LinearOperator` wrapping `_h_eff_action`, above `_DENSE_EIG_MAX`, with a
   deterministic start vector, a residual check, and a dense fallback.
3. **Eigenvectors** kept and returned under `return_vectors=True`.

The dense path stays, and is not vestigial: `dim` can be `1` (a `D=1`
spin-1/2 chain) where ARPACK cannot run at all (it needs `nev < dim`), it is
genuinely faster below the threshold, and it is the fallback when the
iterative solve's residual check fails.

`_channel_resolvent` also gained its own dense-vs-iterative threshold
(`_RESOLVENT_DENSE_MAX`), set by memory rather than flops and far above
`_solve_linear_map`'s `_DENSE_SOLVE_MAX` — because a cached, reused
factorization amortizes an `O(D⁶)` cost that a one-shot solve cannot.

**v3 parity is deliberately deferred.** `Chain::vumps_build_h_eff_dense`
(`mpscpp3/chain_session.h`) is an independent C++ port with the same
rebuild-per-application pattern in `vx_regularized_solve`'s callers. It
still assembles `H_eff(k)` densely, so the two backends now differ in cost
but not in results (`tests/test_vumps_excitations_v3.py` is unchanged and
still passes). Mirroring at least the resolvent cache there is the obvious
follow-up; a C++ Lanczos is separate work again.

## 2. Spectral weights → S(k,ω) on the infinite chain (builds on 1)

With `X` in hand, `B = V_L @ X` is the quasiparticle tensor and
`⟨Ψ|O|Φ_k(B)⟩` gives the **exact δ-peak spectral weight** of each branch.
That is a momentum-resolved dynamical structure factor directly in the
thermodynamic limit — today the only infinite-chain dynamical correlators
are finite-window approximations (`infinitechain.py::kpm_finite`,
`td_dynamical_correlator`).

The machinery this needs (mixed transfer-matrix resolvents) already exists
in the same module: `_channel_resolvent`, `_mixed_fixed_points`,
`_solve_linear_map`.

References: Vanderstraeten, Haegeman & Verstraete, *Tangent-space methods
for uniform matrix product states*, [arXiv:1810.07006](https://arxiv.org/abs/1810.07006)
(SciPost Phys. Lect. Notes 7 (2019)) for the basic overlap; Osborne &
McCulloch, *Efficient and systematic calculation of arbitrary observables
for the matrix product state excitation ansatz*,
[arXiv:2408.17117](https://arxiv.org/abs/2408.17117) (Phys. Rev. Research
**7**, 023018 (2025)) for the systematic recursive scheme, which also
covers excitations of larger spatial support and multi-particle states.

## 3. Two-site VUMPS (dynamic bond dimension)

`pyitensor/vumps.py` is single-site at a fixed target `D`. Bond growth is
`_grow_initial_state` → `_random_raw_tensor`: the smaller solution embedded
in the top-left block plus random noise, propped up by a D-ramp, multiple
restarts, and a variational safety net. The module docstring is candid that
all of this is load-bearing because plain VUMPS was "genuinely unreliable
for D>1", and that a caller wanting a trustworthy D>2 result on a
less-tested model should still rerun and keep the lowest `e0`.

A two-site VUMPS grows the bond dimension by SVD instead of by noise:
solve one two-site eigenproblem plus two one-site ones, truncate at a
cutoff/`max_chi`. It can also start from a product state. Mirror: TeNPy's
[`TwoSiteVUMPSEngine`](https://tenpy.readthedocs.io/en/latest/reference/tenpy.algorithms.vumps.html);
base algorithm Zauner-Stauber *et al.*,
[arXiv:1701.07035](https://arxiv.org/abs/1701.07035) (Phys. Rev. B **97**,
045145 (2018)), and arXiv:1810.07006 Algorithm 4.

TeNPy's own documentation makes the point that matters most here: the
single-site scheme is incompatible with charge conservation *precisely
because* it grows bonds with random unitaries. So this is also the
prerequisite for item 6.

## 4. Long-range interactions via exponential-decay channels

`infinitechain.py::_canonicalize_hamiltonian` rejects any term spanning
both the previous (L) and next (R) unit cell: couplings live inside a cell
or reach the next one, and that is all.

`idmrg.py::_build_automaton` already carries an "F" accumulator channel
with an identity self-loop plus one pending channel per still-open two-site
term. Adding a channel whose self-loop is *λ* rather than the identity sums
`Σ_r λ^r` exactly, and power-law tails follow by fitting to a few
exponentials — the standard construction from Crosswhite, Doherty & Vidal,
*Applying matrix product operators to model systems with long-range
interactions*, [arXiv:0804.2504](https://arxiv.org/abs/0804.2504)
(Phys. Rev. B **78**, 035116 (2008)).

Two things make this larger than it first looks:

- **Two implementations, not one.** pyitensor's `_build_automaton` is
  shared by `idmrg.py` and `vumps.py`, but `mpscpp3/chain_session.h`'s
  `idmrg_build_row` is an independent C++ port that would need the same
  change.
- **The public API cannot state the feature.** `_canonicalize_hamiltonian`
  validates finite `MultiOperator` term lists, and an exponentially
  decaying coupling is an infinite sum no term list can express. This needs
  a new user-facing spec — something like
  `set_long_range(op_i, op_j, decay=λ)` — plus the matching validation
  changes.

A newer, surrogate-free alternative for genuinely algebraic tails (summing
the tail through the connected transfer matrix rather than fitting
exponentials) is Yang, [arXiv:2606.20522](https://arxiv.org/abs/2606.20522);
probably overkill unless the exponential-fit residual is itself the
physics under study.

## 5. iTDVP: real-time evolution of the uniform state

Infinite-system real time today exists only as the finite
infinite-boundary-condition window (`pyitensor/idmrg_window.py`, and
`mpscpp3/chain_session.h`'s `td_dynamical_correlator_window`). A
uniform-state TDVP (arXiv:1810.07006 §5, on top of Haegeman *et al.*,
*Unifying time evolution and optimization with matrix product states*,
[arXiv:1408.5056](https://arxiv.org/abs/1408.5056), Phys. Rev. B **94**,
165116 (2016)) covers **global quenches** in the thermodynamic limit —
something dmrgpy cannot do at all today.

This complements the window rather than replacing it: a dynamical
correlator is a *local* perturbation on top of a translationally invariant
background, which is exactly what the window construction is for.

## 6. QN/sector conservation for infinite chains (long-term)

Finite chains have it on two backends
(`Many_Body_Chain.set_conserved_sector`, `mpscpp3/get_sites.h` +
`pyitensor/sector.py`). Infinite chains have no quantum numbers at all —
see `infinitechain.py`'s own note that with `ConserveQNs=false` nothing
confines the search to the converged state's particle-number sector.

The largest asymptotic win in this list, and the largest job. It needs
item 3 first, for the reason given there.

---

## Deliberately not doing

- **`n_uc>=3` for `gs_method="idmrg"`.** `pyitensor/vumps_ms.py` already
  covers multi-site unit cells, and `idmrg.py`'s module docstring explains
  why the growing algorithm's `p_L`/`p_R` micro-step pairing cannot be
  extended past `n_uc=2` without a genuine redesign.
- **A standalone fix of the `idmrg_window` `U_list` gauge error** (~6e-7 on
  `S(x,t=0)`, documented in `idmrg.py`'s docstring and pinned by
  `tests/test_idmrg_window.py`). Real, but small. Revisit it alongside item
  5: Wu, [arXiv:2007.15035](https://arxiv.org/abs/2007.15035) (Phys. Rev. B
  **102**, 134306 (2020)) is TDVP for exactly this mixed setup — infinite
  uniform left and right bulk with an arbitrary finite region between — and
  is the principled version of what `idmrg_window.py` hand-rolls.
