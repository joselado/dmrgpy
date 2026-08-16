"""svd(): the one tensor-factorization primitive the rest of this engine is
built on -- every truncation in the library (MPS canonicalization, sum(),
applyMPO(), nmultMPO(), the DMRG/TDVP local update) is "reshape into a
matrix across some index partition, SVD, keep the largest singular values
up to Cutoff/MaxDim, reshape back". mpscpp3/chain_session.h only calls it
directly once (bond_entropy()); every other MPS/MPO-level truncation goes
through this same function via mps.py/mpo.py in a later phase.
"""

import numpy as np

from .index import Index
from .tensor import ITensor, _find


class Spectrum:
    """The truncated singular-value spectrum of one svd() call. `probs` are
    the kept singular values squared and normalized against the *original*
    (pre-truncation) total, so sum(probs) + truncerr == 1 -- mirrors
    ITensor's Spectrum::eigs()/truncerr(), which bond_entropy() uses
    directly to compute the von Neumann entropy of a bond."""

    def __init__(self, singular_values, probs, truncerr):
        self.singular_values = singular_values
        self._probs = probs
        self.truncerr = truncerr

    def eigs(self):
        return self._probs


def _truncate(s, cutoff, maxdim, mindim):
    """s: singular values, sorted descending. Returns (keep, discarded_weight)
    -- how many of the largest to keep, and the normalized weight (sum of
    p_i over the dropped tail) that was discarded. Mirrors ITensor's own
    truncation rule: drop the smallest singular values first, stopping as
    soon as either mindim is reached, or dropping the next one would exceed
    `cutoff` of discarded weight -- except maxdim is a hard cap that's
    enforced regardless of cutoff.

    Vectorized: the cutoff floor is equivalent to keeping the smallest
    prefix of the ascending cumulative tail sum that stays <= cutoff (ties
    broken like the old one-at-a-time loop, which drops p[keep-1] only when
    the running total including it does *not* exceed cutoff), then clamped
    to [mindim, maxdim or n]."""
    n = len(s)
    p = s.astype(float) ** 2
    total = p.sum()
    if total <= 0:
        return max(1, mindim), 0.0
    p = p / total
    hi = min(maxdim, n) if maxdim is not None else n  # max keep allowed (maxdim cap)
    j_upper = max(0, n - mindim)  # can never drop past mindim, even to satisfy maxdim
    j_forced = min(max(0, n - hi), j_upper)  # maxdim-mandated drops, still capped by mindim
    if j_upper <= j_forced:
        j = j_forced
    else:
        # tail[k] = sum of the k smallest values beyond the forced drops
        # (p is descending, so the smallest are its last entries);
        # monotonically increasing in k, so the largest cutoff-satisfying
        # k is a simple prefix scan.
        p_asc = p[::-1]
        forced = float(p_asc[:j_forced].sum())
        tail = forced + np.concatenate(([0.0], np.cumsum(p_asc[j_forced:j_upper])))
        window = tail <= cutoff
        j = j_forced + (int(np.nonzero(window)[0].max()) if window.any() else 0)
    keep = n - j
    return keep, float(p[keep:].sum())


def eigh_truncate(rho, cutoff, maxdim, mindim=1):
    """Truncated eigendecomposition of a Hermitian matrix rho (typically a
    reduced/companion density matrix): diagonalize, sort eigenvalues
    descending, treat their (clipped, non-negative) square roots as an
    SVD-like singular-value spectrum, and keep the largest via _truncate()'s
    usual Cutoff/MaxDim/mindim rule. Returns (Uk, keep): Uk is rho's
    eigenvectors for the kept eigenvalues as columns (descending order),
    keep is how many were kept. Factored out because pyitensor/gse.py's
    Krylov-subspace basis enrichment and pyitensor/nhdmrg.py's fidelity
    truncation both reduce to exactly this step, once each has built its
    own (differently constructed) rho."""
    evals, evecs = np.linalg.eigh(rho)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    svals = np.sqrt(np.clip(evals, 0.0, None))
    keep, _discarded = _truncate(svals, cutoff, maxdim, mindim)
    return evecs[:, :keep], keep


# Smallest kept singular value, relative to the largest, that the
# Gram-matrix route below is still trusted for. It diagonalizes M M^dag (or
# M^dag M) instead of factorizing M, so a singular value comes back as the
# square root of an eigenvalue known only to ~eps*s_max**2 absolute -- i.e.
# with a *relative* error ~eps*(s_max/s)**2/2, which stays negligible only
# while s/s_max stays well above sqrt(eps)~1.5e-8. Anything the caller's
# Cutoff actually keeps is normally far above that (Cutoff is a weight,
# s**2/sum(s**2), so the library's own 1e-12 corresponds to s/s_max ~ 1e-6),
# but a caller keeping a genuinely rank-deficient block down to maxdim would
# not be -- hence the explicit guard and exact-SVD fallback in
# _svd_truncated() rather than an unconditional switch.
_GRAM_MIN_RELATIVE_SV = 1e-7

# Below this size the O(min(m,n)^3) eigendecomposition doesn't dominate
# anything and the extra matmuls aren't worth their own call overhead.
_GRAM_MIN_DIM = 16


def _svd_truncated(mat, cutoff, maxdim, mindim):
    """(U_keep, S_all, Vh_keep, keep, discarded) for `mat`, i.e. a thin SVD
    already truncated by _truncate()'s Cutoff/MaxDim/mindim rule.

    Prefers a Gram-matrix route -- eigendecompose the smaller of M M^dag /
    M^dag M, take singular values as square roots of its eigenvalues, and
    recover the *kept* factor on the other side by one matmul -- over
    np.linalg.svd of M itself. Same idea ITensor's own
    densityMatrixApplyMPOImpl uses (diagPosSemiDef of a companion density
    matrix rather than an SVD of the MPO-applied tensor), and for exactly
    the same reason: every MPS-algebra truncation here factorizes a
    strongly rectangular matrix -- MPS bond dim d times physical dim on one
    side, MPS bond dim times MPO bond dim on the other -- and only ever
    keeps at most maxdim of the resulting vectors. Measured on the shapes
    the KPM Chebyshev recursion actually produces at L=30, kpmmaxm=50:
    1.3-2.4x faster than np.linalg.svd (100x250: 7.7 -> 3.2 ms; 200x100:
    5.8 -> 2.8 ms), reconstructing the kept part to ~3e-14.

    Falls back to the exact SVD whenever the smallest *kept* singular value
    is too small for the squared spectrum to resolve (see
    _GRAM_MIN_RELATIVE_SV) -- so this is a speedup, never a silent
    accuracy trade."""
    m, n = mat.shape
    if min(m, n) >= _GRAM_MIN_DIM:
        left = m <= n
        gram = mat @ mat.conj().T if left else mat.conj().T @ mat
        evals, evecs = np.linalg.eigh(gram)
        evals = evals[::-1]
        evecs = evecs[:, ::-1]
        S = np.sqrt(np.clip(evals, 0.0, None))
        keep, discarded = _truncate(S, cutoff, maxdim if maxdim else None, max(1, mindim))
        if S[0] > 0.0 and S[keep - 1] >= _GRAM_MIN_RELATIVE_SV * S[0]:
            W = evecs[:, :keep]
            if left:
                U = W
                Vh = (W.conj().T @ mat) / S[:keep, None]
            else:
                Vh = W.conj().T
                U = (mat @ W) / S[None, :keep]
            return U, S, Vh, keep, discarded

    U, S, Vh = np.linalg.svd(mat, full_matrices=False)
    keep, discarded = _truncate(S, cutoff, maxdim if maxdim else None, max(1, mindim))
    return U[:, :keep], S, Vh[:keep, :], keep, discarded


def svd(T, left_inds, cutoff=0.0, maxdim=None, mindim=1, tags="Link"):
    """Split ITensor T into U * S * V (contracting U*S*V reconstructs T up
    to the requested truncation), grouping `left_inds` onto U and every
    other index of T onto V. Returns (U, S, V, spectrum) -- mirrors
    ITensor's own svd(T,U,S,V) (there an in/out-param triple plus a
    returned Spectrum; here a 4-tuple since Python has no out-params).

    U and V each get their own freshly minted bond index (both dimension
    `keep`, both tagged `tags`) rather than sharing one -- so U's columns
    and V's rows are each independently orthonormal, and callers can absorb
    S into either side (`U*S` for a right-canonical partner, `S*V` for a
    left-canonical one) depending on which way an MPS sweep is moving,
    exactly as ITensor's own three-tensor split allows.
    """
    left_inds = list(left_inds)
    for ind in left_inds:
        if not T.hasindex(ind):
            raise ValueError("svd: {} is not an index of {}".format(ind, T))
    right_inds = [ind for ind in T.inds if _find(left_inds, ind) is None]

    order = left_inds + right_inds
    arr = T.transpose_to(order)
    ldim = int(np.prod([ind.dim for ind in left_inds], dtype=int)) if left_inds else 1
    rdim = int(np.prod([ind.dim for ind in right_inds], dtype=int)) if right_inds else 1
    mat = arr.reshape(ldim, rdim)

    U, S, Vh, keep, discarded = _svd_truncated(mat, cutoff, maxdim, mindim)

    probs_full = (S.astype(float) ** 2)
    total = probs_full.sum()
    probs = probs_full[:keep] / total if total > 0 else probs_full[:keep]

    bond_u = Index(keep, tags=tags)
    bond_v = Index(keep, tags=tags)
    left_shape = tuple(ind.dim for ind in left_inds) + (keep,)
    right_shape = (keep,) + tuple(ind.dim for ind in right_inds)

    Utensor = ITensor(tuple(left_inds) + (bond_u,), U.reshape(left_shape))
    Vtensor = ITensor((bond_v,) + tuple(right_inds), Vh.reshape(right_shape))
    Stensor = ITensor((bond_u, bond_v), np.diag(S[:keep].astype(complex)))
    spectrum = Spectrum(S[:keep], probs, discarded)
    return Utensor, Stensor, Vtensor, spectrum


def qr_split(T, left_inds, tags="Link", orthonormal="left"):
    """Lossless split T = A * B via QR, for the specific case svd() is
    otherwise called with cutoff=0, maxdim=None -- i.e. no truncation is
    wanted, only an orthogonal basis for one side of the bond. QR gets
    there in about half the FLOPs of a full SVD (LAPACK's geqrf/orgqr vs
    gesdd) since it never computes the singular values or bothers
    resolving degenerate subspaces, exactly what one-site TDVP's per-site
    forward split needs (tdvp.py's _half_sweep_lr_onesite/
    _half_sweep_rl_onesite). Not a general svd() replacement: unlike SVD,
    plain QR isn't rank-revealing, so callers relying on exact-zero
    singular values being dropped (gse.py, mpsalgebra.py's randomMPS,
    kpm_energy_truncation.py) must keep using svd(cutoff=0, maxdim=None)
    instead.

    orthonormal='left' (QR): groups `left_inds` onto A with orthonormal
    columns (A^H A = I, like svd()'s U), returns (A, B) with B carrying
    every other index of T. orthonormal='right' (LQ, via QR on T^H):
    groups every index *other than* `left_inds` onto B with orthonormal
    rows (B B^H = I, like svd()'s V), returns (A, B) with A carrying
    `left_inds`. Either way A*B reconstructs T exactly, and the two
    tensors share one freshly minted bond Index (tagged `tags`) -- unlike
    svd(), which mints one for each side."""
    left_inds = list(left_inds)
    for ind in left_inds:
        if not T.hasindex(ind):
            raise ValueError("qr_split: {} is not an index of {}".format(ind, T))
    right_inds = [ind for ind in T.inds if _find(left_inds, ind) is None]

    order = left_inds + right_inds
    arr = T.transpose_to(order)
    ldim = int(np.prod([ind.dim for ind in left_inds], dtype=int)) if left_inds else 1
    rdim = int(np.prod([ind.dim for ind in right_inds], dtype=int)) if right_inds else 1
    mat = arr.reshape(ldim, rdim)

    if orthonormal == "left":
        Amat, Bmat = np.linalg.qr(mat, mode="reduced")
    elif orthonormal == "right":
        Q, R = np.linalg.qr(mat.conj().T, mode="reduced")
        Amat, Bmat = R.conj().T, Q.conj().T
    else:
        raise ValueError("qr_split: orthonormal must be 'left' or 'right', got {}".format(orthonormal))

    k = Amat.shape[1]
    bond = Index(k, tags=tags)
    Atensor = ITensor(tuple(left_inds) + (bond,), Amat.reshape(tuple(ind.dim for ind in left_inds) + (k,)))
    Btensor = ITensor((bond,) + tuple(right_inds), Bmat.reshape((k,) + tuple(ind.dim for ind in right_inds)))
    return Atensor, Btensor, bond
