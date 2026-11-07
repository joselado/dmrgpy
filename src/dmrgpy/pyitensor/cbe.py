"""Controlled bond expansion (CBE) for one-site DMRG and one-site TDVP.

References
----------
- A. Gleis, J.-W. Li, J. von Delft, "Controlled bond expansion for density
  matrix renormalization group ground state search at single-site costs",
  Phys. Rev. Lett. 130, 246402 (2023), arXiv:2207.14712 (CBE-DMRG, the
  shrewd selection of its Fig. 1, Fig. S-2 and Table I).
- A. Gleis, J.-W. Li, J. von Delft, "Projector formalism for kept and
  discarded spaces of matrix product states", Phys. Rev. B 106, 195138
  (2022), arXiv:2207.13161 (the kept/discarded projectors used below).
- J.-W. Li, A. Gleis, J. von Delft, "Time-dependent variational principle
  with controlled bond expansion for matrix product states",
  arXiv:2208.10972 (CBE-TDVP).
- FiniteMPS.jl by Qiaoyi Li (MIT license,
  https://github.com/Qiaoyi-Li/FiniteMPS.jl), src/Algorithm/CBE/ and
  src/Algorithm/DMRG.jl, TDVP/TDVP1.jl. Its `NaiveCBE` is mirrored here as
  selection="naive", and its sweep structure (expand the bond, then the
  one-site update, then a truncating SVD that moves the center) is the one
  both drivers below follow.

The idea
--------
One-site DMRG and one-site TDVP never leave the variational space set by
the current bonds. Two-site updates do, through the part of H|psi> that
lies in the two-site discarded space DD = (1-A A^dag) x (1-B B^dag) of the
bond. CBE adds to the bond only the few directions of DD on which
H|psi> has weight, so the one-site update can use them. For a bond l with
the orthogonality center on site l+1 (right-to-left sweep), psi = A_l C_{l+1}:

  A_l  ->  A_l^ex  = [A_l | A~]       (A~ orthogonal to A_l, D~ new columns)
  C_l+1 -> C^ex    = [C_l+1 ; 0]      (zero rows, so A^ex C^ex = A C)

The state is unchanged; the bond is D+D~. The one-site update of C^ex
then runs on the bigger space, and the SVD that moves the center trims
the bond back to maxdim (and the chain cutoff). The left-to-right sweep is
the mirror image, expanding B_{l+1} instead.

How A~ is chosen
----------------
selection="naive" (FiniteMPS.jl's NaiveCBE): build the two-site H|psi>,
project both legs onto the discarded spaces, SVD across the bond and keep
the D~ largest left singular vectors. Exact minimization of the cost
function C1 of the paper, but the two-site object costs O(D^3 d^2 w) and its
SVD O(D^3 d^3), so the expansion is two-site cost even though the update
is one-site.

selection="shrewd" (default; the paper's Table I): two truncations, each
at one-site cost O(D^3 d w). Preselection SVDs the right half of
H|psi> projected on (1-B^dag B), contracts its weights into the left half,
projects that on (1-A A^dag) and keeps D' = D_ex/w singular vectors of the
central MPS bond (D_ex the expanded bond dimension, w the MPO bond
dimension); regrouping the MPO bond with those D' gives the preselected
complement A^pr with at most D' w ~ D columns. Final selection projects
H|psi> onto A^pr x (1-B^dag B), with the MPO bond closed, and keeps its
D~ largest singular vectors, A~ = A^pr u~. A^pr is explicitly
re-orthogonalized against A (the paper's optional step 14), since
A^dag A^pr = 0 must hold exactly.

How many directions
-------------------
DMRG: expand to D_ex = ceil((1+delta)*maxdim) as in the paper (they use
delta=0.1), D~ = D_ex - D, then trim back to maxdim. From a product state
(D=1) this lets the bond grow as fast as the rank of the selection allows,
which is what the sweep bond-dimension ramp then controls. Singular values
are dropped only when they are numerically zero (1e-12 relative to the
norm of the object before the projection).

TDVP: the same rank cap, plus the paper's singular value thresholds,
1e-4 for the preselection and 1e-6 for the final selection. We take them
relative to the norm of the same object before the projection onto the
discarded space removes what the state already spans, so they mean the
same thing whatever the units of H.

Integrator
----------
CBE-TDVP here uses the same symmetric second-order step as the rest of
pyitensor/tdvp.py (a left-to-right half sweep by dt/2, then a
right-to-left one by dt/2), not the third-order composition of the paper,
so that it is compared like for like with TDVP and TDVP_GSE. The QR of
one-site TDVP is replaced by a truncating SVD (cutoff and maxdim), as in
the paper.

Arrays
------
Every array lives on the backend (backend.xp()); only singular values come
to the host, to count how many to keep. Nothing is written in place. The
kernels are written on plain arrays with fixed axis orders:
  L[a, w, a']      left environment: ket link, MPO link, bra link
  R[b, w, b']      right environment, same order
  W[w, s, s', w2]  MPO tensor: left MPO link, ket site, bra site, right
  M[a, s, b]       MPS tensor
with dimension-1 axes where a tensor has no such leg (chain ends).
"""

import math

import numpy as np

from . import backend as bk
from . import dmrg as _dmrg
from .dmrg import (_all_right_environments, _davidson_ground_state,
                   _extend_left, _extend_right, _noisy_start, one_site_heff,
                   zero_site_heff)
from .index import Index
from .mpsalgebra import _link_at, mpo_is_hermitian
from .svd import _svd_truncated, _truncate
from .tensor import ITensor, commonIndex

# Singular values below this fraction of the scale of the object they come
# from are numerical zeros (the paper's step 13 warns that they must be
# dropped, or A^dag A~ = 0 fails).
_NUMERICAL_ZERO = 1e-12

# The thresholds of the CBE-TDVP paper, preselection and final selection.
TDVP_EPS_PRESELECT = 1e-4
TDVP_EPS_FINAL = 1e-6

SELECTIONS = ("shrewd", "naive")


# ---------------------------------------------------------------------------
# dense kernels
# ---------------------------------------------------------------------------

def _norm(x):
    xp = bk.xp()
    return float(bk.to_host(xp.linalg.norm(x)))


def _count(s_host, scale, eps, cap):
    """How many of the descending singular values s_host to keep: those
    above max(eps, _NUMERICAL_ZERO)*scale, at most cap."""
    if scale <= 0 or len(s_host) == 0:
        return 0
    floor = max(eps, _NUMERICAL_ZERO) * scale
    keep = int(np.count_nonzero(np.asarray(s_host).real > floor))
    if cap is not None:
        keep = min(keep, int(cap))
    return max(keep, 0)


def _svd(mat):
    xp = bk.xp()
    U, s, Vh = xp.linalg.svd(mat, full_matrices=False)
    return U, s, Vh, np.asarray(bk.to_host(s)).real


def _row_basis(C):
    """Orthonormal rows spanning the row space of C (k, t, c) over (t, c):
    the right isometry B of C = Lambda B, rank-revealing, as a (r, t*c)
    matrix."""
    k = C.shape[0]
    mat = C.reshape(k, -1)
    _, _, Vh, s = _svd(mat)
    r = _count(s, float(s[0]) if len(s) else 0.0, 0.0, None)
    return Vh[:r]


def _orthonormal_complement_columns(U, A):
    """U (n, m) with the span of A's columns (n, k) removed, and made
    orthonormal again; numerically zero columns dropped. The paper's step
    14, the safety orthogonalization, applied twice."""
    xp = bk.xp()
    for _ in range(2):
        U = U - A @ (A.conj().T @ U)
    if U.shape[1] == 0:
        return U
    Q, s, _, s_host = _svd(U)
    r = _count(s_host, 1.0, 0.0, None)  # U had orthonormal columns: scale 1
    return Q[:, :r]


def expand_left_isometry(L, Wl, Wr, R, A, C, dtilde, dpre,
                         eps_pre=0.0, eps_final=0.0, selection="shrewd"):
    """New columns A~ (a, s, D~) for the left isometry A (a, s, k) of a
    bond whose center C (k, t, c) sits on its right, orthonormal and
    orthogonal to A's columns. L, Wl, Wr, R as in the module docstring.
    dtilde caps the number of new columns, dpre the preselection rank
    (shrewd only). Returns an array with at most dtilde columns, possibly
    none."""
    xp = bk.xp()
    da, ds, k = A.shape
    if dtilde <= 0:
        return xp.zeros((da, ds, 0), dtype=complex)
    Am = A.reshape(da * ds, k)
    Vb = _row_basis(C)                         # (r, t*c)
    # X[p, x, v, k] = L W_l A : the left half of H|psi>, open at the
    # central MPO bond v and the ket bond k
    X = xp.einsum("awp,wsxv,ask->pxvk", L, Wl, A, optimize=True)
    # Y[k, v, y, q] = C W_r R : the right half
    Y = xp.einsum("ktc,vtyu,cuq->kvyq", C, Wr, R, optimize=True)
    dv, dt, dc = Y.shape[1], Y.shape[2], Y.shape[3]
    Xm = X.reshape(da * ds, dv * k)

    if selection == "naive":
        # H_2s psi on (a s) x (t c), both legs projected onto the
        # discarded spaces, then the SVD across the bond (FiniteMPS.jl
        # NaiveCBE without the randomized SVD)
        Yk = Y.transpose(1, 0, 2, 3).reshape(dv * k, dt * dc)
        T = Xm @ Yk
        scale = _norm(T)
        T = T - Am @ (Am.conj().T @ T)
        T = T - (T @ Vb.conj().T) @ Vb
        U, _, _, s = _svd(T)
        keep = _count(s, scale, eps_final, dtilde)
        At = _orthonormal_complement_columns(U[:, :keep], Am)
        return At.reshape(da, ds, At.shape[1])

    if selection != "shrewd":
        raise ValueError("expand_left_isometry: unknown selection %r, "
                         "expected one of %r" % (selection, SELECTIONS))
    # (a) right half projected on (1 - B^dag B), SVD on the ket bond
    Yq = Y.reshape(k * dv, dt * dc)
    Yo = Yq - (Yq @ Vb.conj().T) @ Vb
    U1, s1, _, s1h = _svd(Yo.reshape(k, dv * dt * dc))
    r1 = _count(s1h, _norm(Y), 0.0, None)
    if r1 == 0:
        return xp.zeros((da, ds, 0), dtype=complex)
    US = U1[:, :r1] * s1[None, :r1]            # (k, r1)
    # left half with the right half's weights on the central MPS bond,
    # projected on (1 - A A^dag)
    Lt = (X.reshape(da * ds * dv, k) @ US).reshape(da * ds, dv * r1)
    scale_pre = _norm(Lt)
    Lo = Lt - Am @ (Am.conj().T @ Lt)
    # (b) truncate the central MPS bond r1 -> D'
    U2, s2, _, s2h = _svd(Lo.reshape(da * ds * dv, r1))
    r2 = _count(s2h, scale_pre, eps_pre, dpre)
    if r2 == 0:
        return xp.zeros((da, ds, 0), dtype=complex)
    us = (U2[:, :r2] * s2[None, :r2]).reshape(da * ds, dv * r2)
    # (c) regroup the MPO bond with the D' kept states: the preselected
    # complement A^pr, with numerical zeros dropped
    U3, _, _, s3h = _svd(us)
    r3 = _count(s3h, float(s3h[0]), 0.0, None)
    Apr = _orthonormal_complement_columns(U3[:, :r3], Am)
    if Apr.shape[1] == 0:
        return xp.zeros((da, ds, 0), dtype=complex)
    # (d) final selection: H|psi> on A^pr x (1 - B^dag B), MPO bond closed
    Lpr = Apr.conj().T @ Xm                     # (x, v*k)
    Ct = Lpr.reshape(-1, dv, k).transpose(0, 2, 1).reshape(-1, k * dv) @ Yq
    scale_fin = _norm(Ct)
    Co = Ct - (Ct @ Vb.conj().T) @ Vb
    U4, _, _, s4h = _svd(Co)
    r4 = _count(s4h, scale_fin, eps_final, dtilde)
    At = Apr @ U4[:, :r4]
    return At.reshape(da, ds, r4)


def expand_right_isometry(L, Wl, Wr, R, C, B, dtilde, dpre,
                          eps_pre=0.0, eps_final=0.0, selection="shrewd"):
    """Mirror of expand_left_isometry: new rows B~ (D~, t, c) for the right
    isometry B (k, t, c) of a bond whose center C (a, s, k) sits on its
    left. Computed as the left problem of the reflected chain."""
    Lr = R
    Wlr = Wr.transpose(3, 1, 2, 0)
    Wrr = Wl.transpose(3, 1, 2, 0)
    Rr = L
    Ar = B.transpose(2, 1, 0)       # (c, t, k), orthonormal columns
    Cr = C.transpose(2, 1, 0)       # (k, s, a)
    At = expand_left_isometry(Lr, Wlr, Wrr, Rr, Ar, Cr, dtilde, dpre,
                              eps_pre, eps_final, selection)
    return At.transpose(2, 1, 0)


# ---------------------------------------------------------------------------
# ITensor glue
# ---------------------------------------------------------------------------

def svd(T, left_inds, cutoff=0.0, maxdim=None):
    """svd.svd's split and truncation rule (U, S, V), without its bond
    padding, through svd.py's _svd_truncated.

    The expansion projects on the complement of U's columns, 1 - U U^dag,
    which is a projector only when U is an isometry to roundoff. The Gram
    route of _svd_truncated once recovered one factor as W^dag M / s,
    orthonormal only to eps*(s_max/s)^2; on a 40-site Heisenberg chain at
    maxm=50 that gave local eigenvalues of -30 against a ground-state
    energy of -17.54, and the sweep energy went up to -15.2. It now
    refactorizes the kept block, and both factors are isometries to
    roughly 1e-14 (tests/test_cbe.py)."""
    xp = bk.xp()
    left_inds = list(left_inds)
    right_inds = [ix for ix in T.inds if ix not in left_inds]
    ldim = int(np.prod([ix.dim for ix in left_inds])) if left_inds else 1
    rdim = int(np.prod([ix.dim for ix in right_inds])) if right_inds else 1
    mat = T.transpose_to(left_inds + right_inds).reshape(ldim, rdim)
    U, S, _, Vh, keep, _ = _svd_truncated(mat, cutoff, maxdim, 1)
    bu, bv = Index(keep, tags="Link"), Index(keep, tags="Link")
    Ut = ITensor(left_inds + [bu], U[:, :keep].reshape(
        tuple(ix.dim for ix in left_inds) + (keep,)))
    Vt = ITensor([bv] + right_inds, Vh[:keep].reshape(
        (keep,) + tuple(ix.dim for ix in right_inds)))
    St = ITensor((bu, bv), xp.diag(S[:keep].astype(complex)))
    return Ut, St, Vt, None

def _site(psi, i):
    return next(ind for ind in psi.A(i).inds if ind.hastags("Site"))


def _mps_array(T, left, s, right):
    """T as (left, s, right), with a dimension-1 axis for a missing link."""
    order = ([left] if left is not None else []) + [s] + ([right] if right is not None else [])
    arr = T.transpose_to(order)
    shape = (left.dim if left is not None else 1, s.dim,
             right.dim if right is not None else 1)
    return arr.reshape(shape)


def _mpo_array(H, i, s):
    n = H.length()
    W = H.A(i)
    wl = commonIndex(W, H.A(i - 1)) if i > 1 else None
    wr = commonIndex(W, H.A(i + 1)) if i < n else None
    order = ([wl] if wl is not None else []) + [s, s.prime(1)] + ([wr] if wr is not None else [])
    arr = W.transpose_to(order)
    shape = (wl.dim if wl is not None else 1, s.dim, s.dim,
             wr.dim if wr is not None else 1)
    return arr.reshape(shape)


def _env_array(E, Ebra, ket_link, H, i):
    """Environment E (left or right, bordering site i) as (ket, mpo, bra);
    a trivial (1,1,1) one at a chain end."""
    if E is None:
        return bk.asarray(np.ones((1, 1, 1)))
    w = commonIndex(E, H.A(i))
    return E.transpose_to([ket_link, w, Ebra])


def _target_dims(H, i, k, maxdim, delta):
    """(dtilde, dpre) for the bond between sites i and i+1 at current bond
    dimension k: expand to D_ex = ceil((1+delta)*maxdim), preselect
    D' = ceil(D_ex/w)."""
    d_ex = int(math.ceil((1.0 + delta) * maxdim))
    w = commonIndex(H.A(i), H.A(i + 1)).dim
    return max(0, d_ex - k), max(1, int(math.ceil(d_ex / float(w))))


def expand_bond_right(psi, H, i, L, Lbra, R2, R2bra, maxdim, delta,
                      eps_pre=0.0, eps_final=0.0, selection="shrewd"):
    """Expand the bond (i, i+1) of psi, center at site i and site i+1
    right-orthogonal, by rows of the right isometry (left-to-right sweep).
    L is the left environment through site i-1, R2 the right one from site
    i+2. Rewrites psi.A(i) and psi.A(i+1) with a new, larger link; psi
    itself is unchanged. Returns the number of directions added."""
    xp = bk.xp()
    n = psi.length()
    s_i, s_j = _site(psi, i), _site(psi, i + 1)
    ll = _link_at(psi, i, i - 1)
    lk = _link_at(psi, i, i + 1)
    lr = _link_at(psi, i + 1, i + 2)
    k = lk.dim
    dtilde, dpre = _target_dims(H, i, k, maxdim, delta)
    if dtilde == 0:
        return 0
    C = _mps_array(psi.A(i), ll, s_i, lk)
    B = _mps_array(psi.A(i + 1), lk, s_j, lr)
    Larr = _env_array(L, Lbra, ll, H, i)
    Rarr = _env_array(R2, R2bra, lr, H, i + 1)
    Wl, Wr = _mpo_array(H, i, s_i), _mpo_array(H, i + 1, s_j)
    Bt = expand_right_isometry(Larr, Wl, Wr, Rarr, C, B, dtilde, dpre,
                               eps_pre, eps_final, selection)
    nt = Bt.shape[0]
    if nt == 0:
        return 0
    new = Index(k + nt, tags="Link")
    Bex = xp.concatenate([B, Bt.astype(B.dtype)], axis=0)
    Cex = xp.concatenate([C, xp.zeros(C.shape[:2] + (nt,), dtype=C.dtype)], axis=2)
    left_inds = ([ll] if ll is not None else []) + [s_i, new]
    right_inds = [new, s_j] + ([lr] if lr is not None else [])
    psi.set_A(i, ITensor(left_inds, Cex.reshape(tuple(ix.dim for ix in left_inds))))
    psi.set_A(i + 1, ITensor(right_inds, Bex.reshape(tuple(ix.dim for ix in right_inds))))
    return nt


def expand_bond_left(psi, H, i, L2, L2bra, R, Rbra, maxdim, delta,
                     eps_pre=0.0, eps_final=0.0, selection="shrewd"):
    """Expand the bond (i, i+1) of psi, center at site i+1 and site i
    left-orthogonal, by columns of the left isometry (right-to-left
    sweep). L2 is the left environment through site i-1, R the right one
    from site i+2. Returns the number of directions added."""
    xp = bk.xp()
    s_i, s_j = _site(psi, i), _site(psi, i + 1)
    ll = _link_at(psi, i, i - 1)
    lk = _link_at(psi, i, i + 1)
    lr = _link_at(psi, i + 1, i + 2)
    k = lk.dim
    dtilde, dpre = _target_dims(H, i, k, maxdim, delta)
    if dtilde == 0:
        return 0
    A = _mps_array(psi.A(i), ll, s_i, lk)
    C = _mps_array(psi.A(i + 1), lk, s_j, lr)
    Larr = _env_array(L2, L2bra, ll, H, i)
    Rarr = _env_array(R, Rbra, lr, H, i + 1)
    Wl, Wr = _mpo_array(H, i, s_i), _mpo_array(H, i + 1, s_j)
    At = expand_left_isometry(Larr, Wl, Wr, Rarr, A, C, dtilde, dpre,
                              eps_pre, eps_final, selection)
    nt = At.shape[2]
    if nt == 0:
        return 0
    new = Index(k + nt, tags="Link")
    Aex = xp.concatenate([A, At.astype(A.dtype)], axis=2)
    Cex = xp.concatenate([C, xp.zeros((nt,) + C.shape[1:], dtype=C.dtype)], axis=0)
    left_inds = ([ll] if ll is not None else []) + [s_i, new]
    right_inds = [new, s_j] + ([lr] if lr is not None else [])
    psi.set_A(i, ITensor(left_inds, Aex.reshape(tuple(ix.dim for ix in left_inds))))
    psi.set_A(i + 1, ITensor(right_inds, Cex.reshape(tuple(ix.dim for ix in right_inds))))
    return nt


# ---------------------------------------------------------------------------
# CBE-DMRG
# ---------------------------------------------------------------------------

def _heff_diagonal_one_site(L, Lbra, H, i, R, Rbra, order_in):
    """Diagonal of one_site_heff at site i, flattened in order_in's axis
    order, without a matvec (the one-site counterpart of
    dmrg._heff_diagonal). Real, since Heff is Hermitian."""
    xp = bk.xp()
    n = H.length()
    s = next(ind for ind in order_in if ind.hastags("Site"))
    W = _mpo_array(H, i, s)                              # (wl, s, s', wr)
    Wd = xp.einsum("assb->asb", W)
    if L is not None:
        wl = commonIndex(L, H.A(i))
        dL = xp.einsum("lal->la", L.transpose_to([order_in[0], wl, Lbra]))
    else:
        dL = bk.asarray(np.ones((1, 1)))
    if R is not None:
        wr = commonIndex(R, H.A(i))
        dR = xp.einsum("rcr->rc", R.transpose_to([order_in[-1], wr, Rbra]))
    else:
        dR = bk.asarray(np.ones((1, 1)))
    d = xp.einsum("la,asb,rb->lsr", dL, Wd, dR, optimize=True)
    return d.reshape(-1).real


def _local_ground_state_one_site(L, Lbra, R, Rbra, H, psi, i, niter, noise,
                                 penalty=None):
    """Lowest eigenpair of the one-site effective Hamiltonian at site i,
    by the same Davidson solve two-site DMRG uses (dmrg.py), with the same
    noisy start. Returns (energy, theta). penalty, when given, is
    (P, PL, PLbra, PR, PRbra): the charge penalty MPO of a conserved sector,
    which H already contains, with its environments bordering site i. Its
    diagonal is taken out of the Davidson stop scale, as in dmrg.dmrg(),
    while the preconditioner keeps the full diagonal."""
    matvec, order_in, shape, x0 = one_site_heff(L, Lbra, H, psi, i, R, Rbra)
    x0 = _noisy_start(x0, noise)
    dim = x0.size
    if dim <= 3:
        xp = bk.xp()
        basis = bk.eye(dim)
        Hmat = xp.column_stack([matvec(basis[:, j]) for j in range(dim)])
        w, v = xp.linalg.eigh((Hmat + Hmat.conj().T) / 2)
        eval0, evec0 = w[0], v[:, 0]
    else:
        diag = _heff_diagonal_one_site(L, Lbra, H, i, R, Rbra, order_in)
        hdiag = None
        if penalty is not None:
            P, PL, PLbra, PR, PRbra = penalty
            hdiag = diag - _heff_diagonal_one_site(PL, PLbra, P, i, PR, PRbra,
                                                   order_in)
        eval0, evec0 = _davidson_ground_state(
            matvec, x0, diag, max(niter, _dmrg.DAVIDSON_NITER),
            _dmrg._DAVIDSON_RTOL, hdiag=hdiag, vtol=_dmrg._DAVIDSON_VTOL)
    theta = ITensor(tuple(order_in), evec0.reshape(shape))
    return float(np.real(bk.to_host(eval0))), theta


def _cbe_sweep(psi, H, left_env, right_env, maxdim, cutoff, niter, noise,
               delta, expand, selection, penalty=None, pleft=None,
               pright=None):
    """One left-to-right then right-to-left CBE sweep. psi enters with its
    center on site 1 and right_env valid for it; it leaves the same way,
    with right_env rebuilt. Returns the last local eigenvalue. penalty is
    the charge penalty MPO H contains (dmrg_cbe), and pleft/pright are its
    environments, kept and rebuilt exactly where left_env/right_env are.
    The bond expansion uses the full H."""
    n = psi.length()
    energy = None
    for i in range(1, n):
        L, Lbra = left_env[i - 1]
        R2, R2bra = right_env[i + 2]
        if expand:
            expand_bond_right(psi, H, i, L, Lbra, R2, R2bra, maxdim, delta,
                              selection=selection)
        R1, R1bra = _extend_right(R2, R2bra, H, psi, i + 1)
        pen = None
        if penalty is not None:
            # built after the expansion, which rewrote psi.A(i+1)
            PR1 = _extend_right(pright[i + 2][0], pright[i + 2][1], penalty,
                                psi, i + 1)
            pen = (penalty,) + pleft[i - 1] + PR1
        energy, theta = _local_ground_state_one_site(L, Lbra, R1, R1bra, H,
                                                     psi, i, niter, noise, pen)
        ll = _link_at(psi, i, i - 1)
        U, S, V, _ = svd(theta, ([ll] if ll is not None else []) + [_site(psi, i)],
                         cutoff=cutoff, maxdim=maxdim)
        psi.set_A(i, U)
        psi.set_A(i + 1, (S * V) * psi.A(i + 1))
        psi.center = i + 1
        left_env[i] = _extend_left(L, Lbra, H, psi, i)
        if penalty is not None:
            pleft[i] = _extend_left(pleft[i - 1][0], pleft[i - 1][1], penalty,
                                    psi, i)
    for i in range(n, 1, -1):
        L2, L2bra = left_env[i - 2]
        R, Rbra = right_env[i + 1]
        if expand:
            expand_bond_left(psi, H, i - 1, L2, L2bra, R, Rbra, maxdim, delta,
                             selection=selection)
        L1, L1bra = _extend_left(L2, L2bra, H, psi, i - 1)
        pen = None
        if penalty is not None:
            # built after the expansion, which rewrote psi.A(i-1)
            PL1 = _extend_left(pleft[i - 2][0], pleft[i - 2][1], penalty,
                               psi, i - 1)
            pen = (penalty,) + PL1 + pright[i + 1]
        energy, theta = _local_ground_state_one_site(L1, L1bra, R, Rbra, H,
                                                     psi, i, niter, noise, pen)
        ll = _link_at(psi, i, i - 1)
        U, S, V, _ = svd(theta, [ll], cutoff=cutoff, maxdim=maxdim)
        psi.set_A(i, V)
        psi.set_A(i - 1, psi.A(i - 1) * (U * S))
        psi.center = i - 1
        right_env[i] = _extend_right(R, Rbra, H, psi, i)
        if penalty is not None:
            pright[i] = _extend_right(pright[i + 1][0], pright[i + 1][1],
                                      penalty, psi, i)
    return energy


def dmrg_cbe(psi, H, sweeps, delta=0.1, quiet=True, noisy_start=True,
             expand=True, selection="shrewd", penalty=None):
    """Ground-state DMRG with one-site updates and controlled bond
    expansion. Same call shape and return value as dmrg.dmrg(): mutates psi
    in place and returns the last local eigenvalue. Each sweep schedule
    entry (maxdim, cutoff, noise, niter) is used as there; the noise goes
    into the first sweep of a fresh start only (dmrg._noisy_start, same
    gate as dmrg()). delta is the expansion fraction of the paper:
    bonds are expanded to ceil((1+delta)*maxdim) and trimmed back to
    maxdim. expand=False is plain one-site DMRG, which never grows a bond;
    it exists for tests. selection is "shrewd" (one-site cost, the paper)
    or "naive" (two-site cost, FiniteMPS.jl's NaiveCBE). penalty is the
    charge penalty MPO that H already contains in a conserved sector, as in
    dmrg.dmrg(): its environments are carried along and its diagonal is
    taken out of the Davidson stop scale of the one-site solves, which
    otherwise stop on the scale of the penalty. The bond expansion still
    uses the full H. None, outside a sector, changes nothing."""
    from .tdvp import _gauge_center
    if selection not in SELECTIONS:
        raise ValueError("dmrg_cbe: unknown selection %r, expected one of %r"
                         % (selection, SELECTIONS))
    _dmrg._NOISE_RNG = np.random.default_rng(0)  # same solve twice, bit for bit
    n = psi.length()
    _gauge_center(psi, 1)
    right_env = _all_right_environments(H, psi)
    left_env = {0: (None, None)}
    pleft = pright = None
    if penalty is not None:
        pright = _all_right_environments(penalty, psi)
        pleft = {0: (None, None)}
    energy = None
    for sweep_i in range(sweeps.nsweep):
        maxdim, cutoff, noise, niter = sweeps.at(sweep_i)
        energy = _cbe_sweep(psi, H, left_env, right_env, maxdim, cutoff,
                            niter, noise if noisy_start and sweep_i == 0 else 0.0,
                            delta, expand, selection, penalty, pleft, pright)
        if not quiet:
            print("CBE sweep {}: energy = {}".format(sweep_i, energy))
    psi.center = 1
    return energy


# ---------------------------------------------------------------------------
# CBE-TDVP
# ---------------------------------------------------------------------------

def _krylov(hermitian):
    from .tdvp import _krylov_expm
    return _krylov_expm(hermitian)


def _forward(L, Lbra, H, psi, i, R, Rbra, tau, niter, hermitian):
    matvec, order_in, shape, x0 = one_site_heff(L, Lbra, H, psi, i, R, Rbra)
    out = _krylov(hermitian)(matvec, x0, -1j * tau, niter=niter)
    return ITensor(tuple(order_in), out.reshape(shape))


def _backward(L, Lbra, R, Rbra, C, left_link, right_link, tau, niter, hermitian):
    matvec, order_in, shape, x0 = zero_site_heff(L, Lbra, R, Rbra, C,
                                                 left_link, right_link)
    out = _krylov(hermitian)(matvec, x0, 1j * tau, niter=niter)
    return ITensor(tuple(order_in), out.reshape(shape))


def _half_sweep_lr_cbe(psi, H, tau, cutoff, maxdim, delta, niter, hermitian,
                       selection):
    n = psi.length()
    right_env = _all_right_environments(H, psi)
    left_env = {0: (None, None)}
    for i in range(1, n + 1):
        L, Lbra = left_env[i - 1]
        if i < n:
            R2, R2bra = right_env[i + 2]
            expand_bond_right(psi, H, i, L, Lbra, R2, R2bra, maxdim, delta,
                              TDVP_EPS_PRESELECT, TDVP_EPS_FINAL, selection)
            Rn, Rnbra = _extend_right(R2, R2bra, H, psi, i + 1)
        else:
            Rn, Rnbra = right_env[n + 1]
        A_new = _forward(L, Lbra, H, psi, i, Rn, Rnbra, tau, niter, hermitian)
        if i == n:
            psi.set_A(i, A_new)
            continue
        ll = _link_at(psi, i, i - 1)
        lk = _link_at(psi, i, i + 1)
        U, S, V, _ = svd(A_new, ([ll] if ll is not None else []) + [_site(psi, i)],
                         cutoff=cutoff, maxdim=maxdim)
        new_link = commonIndex(U, S)
        C = S * V                                   # (new_link, lk)
        psi.set_A(i, U)
        nxt = psi.A(i + 1)
        psi.set_A(i + 1, C * nxt)                   # placeholder, see tdvp.py
        left_env[i] = _extend_left(L, Lbra, H, psi, i)
        Lnew, Lnewbra = left_env[i]
        C_ev = _backward(Lnew, Lnewbra, Rn, Rnbra, C, new_link, lk, tau, niter,
                         hermitian)
        psi.set_A(i + 1, C_ev * nxt)
    psi.center = n


def _half_sweep_rl_cbe(psi, H, tau, cutoff, maxdim, delta, niter, hermitian,
                       selection):
    from .dmrg import _all_left_environments
    n = psi.length()
    left_env = _all_left_environments(H, psi)
    right_env = {n + 1: (None, None)}
    for i in range(n, 0, -1):
        R, Rbra = right_env[i + 1]
        if i > 1:
            L2, L2bra = left_env[i - 2]
            expand_bond_left(psi, H, i - 1, L2, L2bra, R, Rbra, maxdim, delta,
                             TDVP_EPS_PRESELECT, TDVP_EPS_FINAL, selection)
            Ln, Lnbra = _extend_left(L2, L2bra, H, psi, i - 1)
        else:
            Ln, Lnbra = None, None
        A_new = _forward(Ln, Lnbra, H, psi, i, R, Rbra, tau, niter, hermitian)
        if i == 1:
            psi.set_A(i, A_new)
            continue
        ll = _link_at(psi, i, i - 1)
        U, S, V, _ = svd(A_new, [ll], cutoff=cutoff, maxdim=maxdim)
        new_link = commonIndex(S, V)
        C = U * S                                   # (ll, new_link)
        psi.set_A(i, V)
        prv = psi.A(i - 1)
        psi.set_A(i - 1, prv * C)                   # placeholder
        right_env[i] = _extend_right(R, Rbra, H, psi, i)
        Rnew, Rnewbra = right_env[i]
        C_ev = _backward(Ln, Lnbra, Rnew, Rnewbra, C, ll, new_link, tau, niter,
                         hermitian)
        psi.set_A(i - 1, prv * C_ev)
    psi.center = 1


def tdvp_cbe_step(psi, H, dt, cutoff, maxdim, delta=0.1, niter=50,
                  hermitian=None, selection="shrewd"):
    """One step exp(-i*dt*H) of one-site TDVP with controlled bond
    expansion before every local forward step: a left-to-right half sweep
    by dt/2, then a right-to-left one by dt/2, as tdvp.tdvp_step. dt may be
    complex (TDZ, METTS). Bonds grow up to maxdim, trimmed with cutoff.
    hermitian picks the local exponentiator exactly as in tdvp.tdvp_step.
    Mutates psi in place and returns it."""
    from .tdvp import _gauge_center
    if selection not in SELECTIONS:
        raise ValueError("tdvp_cbe_step: unknown selection %r, expected one "
                         "of %r" % (selection, SELECTIONS))
    if hermitian is None:
        hermitian = mpo_is_hermitian(H)
    tau = dt / 2.0
    _gauge_center(psi, 1)
    _half_sweep_lr_cbe(psi, H, tau, cutoff, maxdim, delta, niter, hermitian,
                       selection)
    _half_sweep_rl_cbe(psi, H, tau, cutoff, maxdim, delta, niter, hermitian,
                       selection)
    return psi
