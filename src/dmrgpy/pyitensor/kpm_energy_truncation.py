"""Energy truncation for KPM Chebyshev vectors -- Holzner, Weichselbaum,
McCulloch & von Delft, "Chebyshev matrix product state approach for
spectral functions", PRB 83, 195115 (2011), Sec. III-B.

Standalone from the Chebyshev recursion itself (chain.py's
_kpm_moments_full/_kpm_moments_accelerated): this module only knows how
to take one MPS and one (already rescaled/shifted) Hamiltonian MPO and
project high-energy components out of the MPS, one site at a time, in
dedicated truncation sweeps that move the orthogonality center losslessly
(no SVD truncation -- that would fight the energy truncation itself, see
the module-level algorithm summary in the paper's Sec. III-B).

Why this is needed at all: the KPM window Ws can be chosen much smaller
than the true many-body bandwidth W (Ws just needs to cover the
correlator's own, usually much narrower, spectral width), which is what
buys spectral resolution. But shrinking Ws makes it easy for the rescaled
Hamiltonian H' to have eigenstates with |eigenvalue| > 1 -- outside the
domain where Chebyshev polynomials are bounded. Since the recursion
|t_n> = 2 H'|t_{n-1}> - |t_{n-2}> is carried out via MPS compression, not
in H's eigenbasis, numerical noise alone is enough to leak a little
high-energy weight into |t_n>, and since |T_k(x)| grows without bound for
|x|>1, that weight blows up exponentially over subsequent recursion
steps. chain.py's _check_kpm_moment only *detects* this (raises once
moments diverge); this module is the actual fix, applied once per
Chebyshev vector before it is used for a moment or fed into the next
recursion step (mirrors the paper's own recursion ordering).

Two deliberate departures from the paper's literal Eqs. (36)-(41), noted
again at their point of use below:

- Krylov vectors are built via interleaved Gram-Schmidt (matching this
  package's own Lanczos style, dmrg.py's _lanczos_ground_state) rather
  than the paper's "compute all d_K powers of H' first, then batch
  orthogonalize". Both build the same Krylov subspace mathematically, but
  interleaving is far better conditioned once vectors start converging
  toward H's locally dominant eigendirection -- exactly the failure mode
  classical (as opposed to modified/interleaved) Gram-Schmidt is known to
  suffer from.
- The projector keeps Krylov components with |eps_alpha| < threshold
  (both signs), not just eps_alpha < threshold as in Eq. (38). The paper
  can use a one-sided cut because its rescaling (Eq. 21b) shifts by the
  exact ground-state energy E0, which pins the ground state at -W' and
  guarantees the *lower* end of the rescaled spectrum never leaks past
  -1. chain.py's own _scaled_hamiltonian shifts by the bandwidth
  *midpoint* (emin+emax)/2 instead (a cheaper, DMRG-bound-only estimate
  that doesn't require an exact E0), so once Ws is chosen smaller than W,
  leakage is possible on either side of the rescaled window.
"""
import numpy as np

from .dmrg import (_all_left_environments, _all_right_environments,
                    _extend_left, _extend_right, one_site_heff)
from .mpsalgebra import inner
from .tensor import ITensor


def _local_krylov_projection(matvec, x0, dK, threshold):
    """One site's worth of Eqs. (36)-(39): build an orthonormal Krylov
    basis of dimension <= dK for the local effective Hamiltonian
    `matvec`, starting from this site's current (flattened) tensor `x0`,
    diagonalize the dense projected matrix, and return x0 with any
    component whose Krylov-subspace energy has |eigenvalue| >= threshold
    projected out. Returns (new_x0, truncated_weight_sq) -- the latter is
    the squared norm of the removed part, i.e. one site's contribution to
    Eq. (40)."""
    dim = x0.size
    nrm = np.linalg.norm(x0)
    if nrm == 0.0:
        return x0, 0.0
    dK_eff = max(1, min(dK, dim))
    vs = [x0 / nrm]
    for _ in range(1, dK_eff):
        w = matvec(vs[-1])
        for v in vs:
            w = w - np.vdot(v, w) * v
        wn = np.linalg.norm(w)
        if wn < 1e-12:
            break  # invariant subspace found; fewer than dK vectors is fine
        vs.append(w / wn)
    k = len(vs)
    V = np.column_stack(vs)
    Hk = np.empty((k, k), dtype=complex)
    for j in range(k):
        Hk[:, j] = V.conj().T @ matvec(vs[j])
    Hk = (Hk + Hk.conj().T) / 2  # Hermitize away floating-point asymmetry
    evals, evecs = np.linalg.eigh(Hk)

    c_v = np.zeros(k, dtype=complex)
    c_v[0] = nrm  # x0 == nrm * vs[0] == nrm * V[:, 0] by construction
    c_e = evecs.conj().T @ c_v
    keep = np.abs(evals) < threshold
    truncated_weight_sq = float(np.sum(np.abs(c_e[~keep]) ** 2))
    c_e = np.where(keep, c_e, 0.0)
    new_x0 = V @ (evecs @ c_e)
    return new_x0, truncated_weight_sq


def _truncate_half_sweep(psi, H, dK, threshold, direction):
    """One directional pass over every site of psi (mutated in place),
    projecting each site's local tensor via _local_krylov_projection and
    moving the orthogonality center losslessly (cutoff=0, maxdim=None) in
    `direction` between sites -- the "standard MPS means, without any
    truncation" shift the paper calls for, so this sweep's own site-to-
    site moves never counteract the energy truncation it just did.
    Requires psi.center to already sit at the pass's starting site (site
    1 for "right", site N for "left"). Returns this pass's average
    truncated weight per site, Eq. (40)."""
    n = psi.length()
    total_truncated_sq = 0.0
    if direction == "right":
        env = _all_right_environments(H, psi)
        side_env = {0: (None, None)}
        sites = range(1, n + 1)
        extend = _extend_left
        get_lr = lambda i, side: (side[i - 1][0], side[i - 1][1], env[i + 1][0], env[i + 1][1])
    else:
        env = _all_left_environments(H, psi)
        side_env = {n + 1: (None, None)}
        sites = range(n, 0, -1)
        extend = _extend_right
        get_lr = lambda i, side: (env[i - 1][0], env[i - 1][1], side[i + 1][0], side[i + 1][1])

    for i in sites:
        L, Lbra, R, Rbra = get_lr(i, side_env)
        matvec, order_in, shape, x0 = one_site_heff(L, Lbra, H, psi, i, R, Rbra)
        new_x0, truncated_sq = _local_krylov_projection(matvec, x0, dK, threshold)
        total_truncated_sq += truncated_sq
        psi.set_A(i, ITensor(order_in, new_x0.reshape(shape)))

        is_last = (i == n) if direction == "right" else (i == 1)
        if not is_last:
            next_i = i + 1 if direction == "right" else i - 1
            psi.position(next_i)  # lossless: cutoff=0.0, maxdim=None
            if direction == "right":
                side_env[i] = extend(L, Lbra, H, psi, i)
            else:
                side_env[i] = extend(R, Rbra, H, psi, i)

    return (total_truncated_sq / n) ** 0.5


def energy_truncate(psi, H, dK=30, n_sweeps=10, threshold=1.0):
    """Project high-rescaled-energy components out of a Chebyshev vector
    `psi` (mutated in place; also returned for convenience). `H` must be
    the same rescaled/shifted Hamiltonian MPO used to build psi's own
    Chebyshev recursion (H' in the paper's notation) -- NOT the bare,
    unscaled Hamiltonian.

    Runs `n_sweeps` directional passes over the chain, alternating
    direction each pass (matching every other sweep-based algorithm in
    this package, e.g. dmrg.py's _dmrg_one_sweep) so the orthogonality
    center never needs an extra, wasted end-to-end reset between passes.

    Returns (avg_truncated_weight, state_change_norm):
      - avg_truncated_weight: this call's *last* pass's average per-site
        truncated weight, Eq. (40) -- once truncation has converged
        (residual weight stops decreasing across sweeps, as in the
        paper's own Fig. 10), this is the leftover high-energy weight.
      - state_change_norm: |||t_n>_tr - |t_n>||^2, Eq. (41), measured
        globally via a plain overlap of psi before vs. after every pass,
        not summed from the individual per-site changes above (which
        would double-count/miscount once more than one site has been
        touched in the same pass under a shifting gauge).
    """
    # raises if psi.center is unset -- every MPS this is called on
    # (mps_sum/applyMPO output) already has it set
    psi_before = psi.copy()
    psi.position(1)

    avg_truncated_weight = 0.0
    direction = "right"
    for _ in range(n_sweeps):
        avg_truncated_weight = _truncate_half_sweep(psi, H, dK, threshold, direction)
        direction = "left" if direction == "right" else "right"

    state_change_norm = (inner(psi_before, psi_before).real
                          + inner(psi, psi).real
                          - 2.0 * inner(psi_before, psi).real)
    return avg_truncated_weight, max(state_change_norm, 0.0)
