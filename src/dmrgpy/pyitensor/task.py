"""Tangent-space Krylov (TaSK) spectral functions on top of a ground-state MPS
(Kovalska, von Delft and Gleis, arXiv:2510.07279).

The idea
--------
A Lanczos recursion on MPS (Hallberg's continued fraction, or the global
Krylov of rootndmrg.py) has to compress after every H|v>, and the
compression destroys the orthogonality of the Krylov basis. TaSK runs the
recursion inside the one-site tangent space of the ground state instead,
which is an honest vector space: a vector is

    |Psi> = sum_l A_1 ... A_{l-1} T_l B_{l+1} ... B_L

with A_l (B_l) the left (right) isometries of the ground state and T_l in
the discarded space of A_l, sum_{a,s} conj(A_l[a,s,m]) T_l[a,s,r] = 0. The
gauge makes the L terms mutually orthogonal and orthogonal to |GS>, so
overlaps, sums and scalings act on the T_l alone, exactly, at O(L D^2 d).
The spectral function

    S(w) = <GS|A delta(w - H + E0) B|GS>

is approximated by projecting H - E0 and B|GS> onto that space,
H^{1perp} = P (H - E0) P and |b> = P B|GS>, and running an ordinary Lanczos
with full reorthogonalization. The Ritz pairs (w_alpha, |alpha>) give the
discrete poles and weights; no broadening enters the solve. The one
approximation is the projection P, which is exact only if the ground-state
bond dimension saturates the Hilbert space; otherwise an excitation that
needs more than a one-site change of the ground state (a multi-spinon
state, for instance) is represented through the tangent vectors closest
to it. Its size is estimated per pole by the two-site variance
||P^{2perp} H |alpha>||^2 (Hubig et al., Gleis et al.), the part of
H|alpha> that falls into the two-site discarded space.

The elastic term n = GS is not in the tangent space (which is orthogonal
to |GS>) and is returned separately as <GS|A|GS><GS|B|GS> at w = 0.

Off-diagonal pairs
------------------
The paper treats A = B^dagger. For a general pair the weight of a pole is
<GS|A|alpha><alpha|B|GS> = conj(<alpha|a>) <alpha|b> with |a> = P A^dagger|GS>,
and <alpha|a> is evaluated by expanding |a> in the Krylov basis built from
|b>. That is one approximation more than the paper makes: the part of |a>
outside the Krylov space of |b> is dropped, which converges away with
nkry like any Krylov resolvent but is not controlled by the variance.

Everything here is dense arrays of the active pyitensor backend
(backend.xp(): NumPy, or JAX on a device) on (chi_left, d, chi_right) site
arrays and (w_left, d_out, d_in, w_right) MPO arrays; the ITensor objects
are only read at the entry point. The site and MPO arrays arrive on the
host (mpo_arrays, chain._mps_arrays_lpr, the v3 export) and real_if_exact
is the one place they are moved to the backend, so on JAX the environments,
the Krylov basis and every matvec stay on the device and only scalars and
the small Krylov-space matrices come back.

Only one array of the size of the Krylov basis, Q (one row per tangent
vector), is ever alive: H Q, the Ritz vectors U^T Q and O Q are never
formed, their overlaps with Q are accumulated column by column instead.
"""

import numpy as np

from . import backend


# ---------------------------------------------------------------- arrays

def mpo_arrays(H):
    """Every site tensor of the MPO `H` as a (w_left, d_out, d_in, w_right)
    array, d_out the primed (bra) leg and d_in the unprimed (ket) one, with
    length-1 axes at the open ends."""
    from .tensor import commonIndex
    n = H.length()
    out = []
    for i in range(1, n + 1):
        W = H.A(i)
        sout = next(ind for ind in W.inds if ind.hastags("Site") and ind.plev == 1)
        sin = next(ind for ind in W.inds if ind.hastags("Site") and ind.plev == 0)
        left = commonIndex(H.A(i - 1), W) if i > 1 else None
        right = commonIndex(W, H.A(i + 1)) if i < n else None
        order = [ind for ind in (left, sout, sin, right) if ind is not None]
        arr = np.asarray(W.transpose_to(order))
        if left is None:
            arr = arr.reshape((1,) + arr.shape)
        if right is None:
            arr = arr.reshape(arr.shape + (1,))
        out.append(np.ascontiguousarray(arr, dtype=complex))
    return out


def real_if_exact(arrays):
    """The arrays as float64 when every imaginary part is exactly zero,
    unchanged otherwise. A real Hamiltonian, a real ground state and a real
    operator then run the whole recursion in real arithmetic, about four
    times fewer flops per contraction; a complex start vector on real
    environments is still correct, numpy promotes it."""
    xp = backend.xp()
    arrays = [backend.to_host(a) for a in arrays]
    if all(not np.iscomplexobj(a) or not np.any(a.imag) for a in arrays):
        return [xp.asarray(np.ascontiguousarray(np.real(a), dtype=float))
                for a in arrays]
    return [xp.asarray(np.asarray(a, dtype=complex)) for a in arrays]


def mixed_gauge(Ms):
    """Left isometries A_l, right isometries B_l and site-centre tensors
    C_l of the normalized MPS given by the arrays `Ms`, with A_l and B_l of
    the same shape, so that

        |GS> = A_1 ... A_{l-1} C_l B_{l+1} ... B_L   for every l.

    A left QR sweep then a right LQ sweep bring every bond to at most d
    times its neighbours, after which a second left QR sweep from the
    right-canonical form gives A_l of exactly the B_l bond dimensions (a
    rank-deficient bond keeps its dimension; the A_l there is completed
    arbitrarily by the QR, which leaves both gauge relations intact)."""
    xp = backend.xp()
    Ms = list(real_if_exact(Ms))
    L = len(Ms)
    # left QR sweep (only to shrink structurally oversized bonds)
    for l in range(L - 1):
        a, d, b = Ms[l].shape
        Q, R = xp.linalg.qr(Ms[l].reshape(a * d, b))
        Ms[l] = Q.reshape(a, d, Q.shape[1])
        Ms[l + 1] = xp.tensordot(R, Ms[l + 1], axes=(1, 0))
    # right LQ sweep -> right isometries B
    Bs = [None] * L
    X = Ms[L - 1]
    for l in range(L - 1, 0, -1):
        a, d, b = X.shape
        Q, R = xp.linalg.qr(X.reshape(a, d * b).T)  # X = R^T Q^T
        Bs[l] = Q.T.reshape(Q.shape[1], d, b)
        X = xp.tensordot(Ms[l - 1], R.T, axes=(2, 0))
    nrm = float(xp.linalg.norm(X))
    if nrm == 0.0:
        raise ValueError("TaSK: the ground state has zero norm")
    Bs[0] = X / nrm  # the first site carries the norm; a (1,d,D1) unit vector
    # left QR sweep from the right-canonical form -> A, C
    As, Cs = [None] * L, [None] * L
    X = Bs[0]
    for l in range(L):
        Cs[l] = X
        if l == L - 1:
            As[l] = X  # a (D,d,1) unit vector, left-normalized
            break
        a, d, b = X.shape
        Q, R = xp.linalg.qr(X.reshape(a * d, b))
        As[l] = Q.reshape(a, d, b)
        X = xp.tensordot(R, Bs[l + 1], axes=(1, 0))
    Bs[0] = Bs[0] / xp.linalg.norm(Bs[0])
    return As, Bs, Cs


# ---------------------------------------------------------------- contractions
#
# Written as tensordot chains, not np.einsum: einsum runs complex
# contractions like these through its own loops rather than BLAS, measured
# 0.04 s per pairwise step at D=64, which made one Lanczos step on a
# 12-site chain take 18 s.

def _ext_left(E, bra, W, ket):
    """Left environment E[a,w,b] (bra, MPO, ket bond) grown by one site."""
    xp = backend.xp()
    X = xp.tensordot(E, ket, axes=(2, 0))  # a w s z
    X = xp.tensordot(X, W, axes=([1, 2], [0, 2]))  # a z t y
    X = xp.tensordot(bra.conj(), X, axes=([0, 1], [0, 2]))  # x z y
    return X.transpose(0, 2, 1)


def _ext_right(E, bra, W, ket):
    """Right environment E[x,y,z] (bra, MPO, ket bond) grown by one site."""
    xp = backend.xp()
    X = xp.tensordot(ket, E, axes=(2, 2))  # b s x y
    X = xp.tensordot(X, W, axes=([1, 3], [2, 3]))  # b x w t
    X = xp.tensordot(bra.conj(), X, axes=([1, 2], [3, 1]))  # a b w
    return X.transpose(0, 2, 1)


def _open1(Le, W, K, Re):
    """The site-l bra tensor of <left env| W |ket K> |right env>."""
    xp = backend.xp()
    X = xp.tensordot(Le, K, axes=(2, 0))  # a w s z
    X = xp.tensordot(X, W, axes=([1, 2], [0, 2]))  # a z t y
    return xp.tensordot(X, Re, axes=([1, 3], [2, 1]))  # a t x


def _open2(Le, W1, W2, K2, Re):
    """Two-site version of _open1 on the bond (l, l+1), K2[b,s,v,c]."""
    xp = backend.xp()
    X = xp.tensordot(Le, K2, axes=(2, 0))  # a w s v c
    X = xp.tensordot(X, W1, axes=([1, 2], [0, 2]))  # a v c t y
    X = xp.tensordot(X, W2, axes=([1, 4], [2, 0]))  # a c t u q
    return xp.tensordot(X, Re, axes=([1, 4], [2, 1]))  # a t u x


def _proj_left(X, A):
    """Remove from X its component in the kept space of the left isometry
    A, on X's first two legs."""
    xp = backend.xp()
    c = xp.tensordot(A.conj(), X, axes=([0, 1], [0, 1]))
    return X - xp.tensordot(A, c, axes=(2, 0))


def _proj_right(X, B):
    """Remove from X its component in the kept space of the right isometry
    B, on X's last two legs."""
    xp = backend.xp()
    c = xp.tensordot(X, B.conj(), axes=([-2, -1], [1, 2]))
    return X - xp.tensordot(c, B, axes=(-1, 0))


def _pair(T1, T2):
    return backend.xp().tensordot(T1, T2, axes=(2, 0))


def _left_envs(As, Ws):
    L = len(As)
    out = [backend.xp().ones((1, 1, 1))]
    for l in range(L):
        out.append(_ext_left(out[-1], As[l], Ws[l], As[l]))
    return out  # out[l] holds sites 0..l-1


def _right_envs(Bs, Ws):
    L = len(Bs)
    out = [None] * (L + 1)
    out[L] = backend.xp().ones((1, 1, 1))
    for l in range(L - 1, -1, -1):
        out[l] = _ext_right(out[l + 1], Bs[l], Ws[l], Bs[l])
    return out  # out[l] holds sites l..L-1


# ---------------------------------------------------------------- tangent space

class TangentSpace:
    """The one-site tangent space of a ground-state MPS, orthogonal to it,
    with the Hamiltonian H projected onto it."""

    def __init__(self, Ms, Ws):
        self.As, self.Bs, self.Cs = mixed_gauge(Ms)
        self.Ws = Ws = real_if_exact(Ws)
        self.L = len(Ms)
        self.LH = _left_envs(self.As, Ws)
        self.RH = _right_envs(self.Bs, Ws)
        self.e0 = float(backend.xp().vdot(
            self.Cs[0], _open1(self.LH[0], Ws[0], self.Cs[0], self.RH[1])).real)

    # vector-space operations on lists of T_l
    def dot(self, X, Y):
        xp = backend.xp()
        return sum(xp.vdot(x, y) for x, y in zip(X, Y))

    def norm(self, X):
        return float(np.sqrt(abs(self.dot(X, X))))

    def flat_size(self):
        return int(sum(a.size for a in self.As))

    def flat(self, X):
        """A tangent vector as one array; overlaps are then vdot."""
        return backend.xp().concatenate([x.ravel() for x in X])

    def unflat(self, v):
        out, k = [], 0
        for a in self.As:
            out.append(v[k:k + a.size].reshape(a.shape))
            k += a.size
        return out

    def gauge(self, X):
        return [_proj_left(x, a) for x, a in zip(X, self.As)]

    def dimension(self):
        return int(sum((a.shape[0] * a.shape[1] - a.shape[2]) * a.shape[2]
                       for a in self.As))

    def expectation(self, op):
        """<GS|O|GS> for an operator from operator()."""
        Wo, LO, RO = op
        y = _open1(LO[0], Wo[0], self.Cs[0], RO[1])
        return complex(backend.xp().vdot(self.Cs[0], y))

    def project_operator(self, Wo):
        """P O|GS> for the MPO arrays Wo, as a tangent vector. Exact for a
        sum of one-site operators, the projection otherwise."""
        Wo = real_if_exact(Wo)
        LO = _left_envs(self.As, Wo)
        RO = _right_envs(self.Bs, Wo)
        return self.gauge([_open1(LO[l], Wo[l], self.Cs[l], RO[l + 1])
                           for l in range(self.L)])

    def operator(self, Wo):
        """An MPO with its A-A and B-B environments, ready for apply()."""
        Wo = real_if_exact(Wo)
        return (Wo, _left_envs(self.As, Wo), _right_envs(self.Bs, Wo))

    def _mixed_envs(self, T, op=None):
        """The environments of the operator `op` (default H) whose ket
        carries the tangent tensor of T to their left (resp. right), A_l
        before it and B_l after it."""
        Ws, LE, RE = op if op is not None else (self.Ws, self.LH, self.RH)
        L, As, Bs = self.L, self.As, self.Bs
        zero = backend.xp().zeros((1, 1, 1))
        Lc = [zero]
        for l in range(L):
            Lc.append(_ext_left(Lc[-1], As[l], Ws[l], Bs[l])
                      + _ext_left(LE[l], As[l], Ws[l], T[l])
                      if l > 0 else
                      _ext_left(LE[l], As[l], Ws[l], T[l]))
        Rc = [None] * (L + 1)
        Rc[L] = zero
        for l in range(L - 1, -1, -1):
            Rc[l] = (_ext_right(RE[l + 1], Bs[l], Ws[l], T[l])
                     + (_ext_right(Rc[l + 1], Bs[l], Ws[l], As[l])
                        if l < L - 1 else 0.0))
        return Lc, Rc

    def apply(self, op, T):
        """P O P acting on the tangent vector T, for an operator from
        operator(); the part of O|T> along |GS> is dropped."""
        Ws, LE, RE = op
        L, As, Bs = self.L, self.As, self.Bs
        Lc, Rc = self._mixed_envs(T, op)
        out = []
        for l in range(L):
            y = _open1(LE[l], Ws[l], T[l], RE[l + 1])
            if l > 0:
                y = y + _open1(Lc[l], Ws[l], Bs[l], RE[l + 1])
            if l < L - 1:
                y = y + _open1(LE[l], Ws[l], As[l], Rc[l + 1])
            out.append(y)
        return self.gauge(out)

    def apply_h(self, T):
        """P (H - E0) P acting on the tangent vector T."""
        y = self.apply((self.Ws, self.LH, self.RH), T)
        return [a - self.e0 * t for a, t in zip(y, T)]

    def two_site_variance(self, T):
        """||P^{2perp} H |T>||^2, the part of H|T> in the two-site discarded
        space, which no one-site tangent vector can represent."""
        L, As, Bs, Ws = self.L, self.As, self.Bs, self.Ws
        Lc, Rc = self._mixed_envs(T)
        tot = 0.0
        for l in range(L - 1):
            Le, Re = self.LH[l], self.RH[l + 2]
            K = _pair(T[l], Bs[l + 1]) + _pair(As[l], T[l + 1])
            y = _open2(Le, Ws[l], Ws[l + 1], K, Re)
            if l > 0:
                KB = _pair(Bs[l], Bs[l + 1])
                y = y + _open2(Lc[l], Ws[l], Ws[l + 1], KB, Re)
            if l < L - 2:
                KA = _pair(As[l], As[l + 1])
                y = y + _open2(Le, Ws[l], Ws[l + 1], KA, Rc[l + 2])
            y = _proj_right(_proj_left(y, As[l]), Bs[l + 1])
            tot += float(backend.xp().vdot(y, y).real)
        return tot


# ---------------------------------------------------------------- Lanczos
#
# A Krylov basis is a preallocated store, one flattened tangent vector per
# row, so that every overlap with it is a matrix product. On NumPy the
# first m rows are a view; on a device a slice would copy the whole basis,
# so the full store is used instead, whose rows beyond m are exactly zero
# and contribute exactly zero to every overlap and every combination.

def _rows(store, m):
    return store if backend.is_device() else store[:m]


def _overlaps(Q, w):
    """Q^* w, the overlaps <q_i|w>, without forming the conjugate of Q,
    which would be a second copy of the basis."""
    return (Q @ w.conj()).conj()


def lanczos(ts, b, nkry, breakdown=1e-10):
    """Lanczos on P(H-E0)P from the tangent vector b, with full
    reorthogonalization done twice per step (the paper's step vii).

    Returns (alpha, beta, store, m): the tridiagonal elements and the
    orthonormal Krylov vectors, the first m rows of the flattened store.
    Stops early when the new vector's norm falls below `breakdown` times
    the norm of H|v_i> it came from, the Krylov space then being
    invariant."""
    xp = backend.xp()
    nb = ts.norm(b)
    if nb == 0.0:
        return np.zeros(0), np.zeros(0), None, 0
    # The basis is kept flattened, one row per vector, so that the
    # reorthogonalization is two matrix products rather than a Python sum
    # over sites per pair of vectors (a third of the time at L=32, D=64).
    v0 = ts.flat(b) / nb
    store = xp.zeros((nkry + 1, v0.size), dtype=v0.dtype)
    store = backend.setrow(store, 0, v0)
    alpha, beta = [], []
    m = 1
    for i in range(nkry + 1):
        w = ts.flat(ts.apply_h(ts.unflat(store[i])))
        if w.dtype != store.dtype:  # a complex P H P on a real start
            store = store.astype(w.dtype)
        scale = float(xp.linalg.norm(w))
        Q = _rows(store, m)
        for rep in range(2):
            c = _overlaps(Q, w)
            if rep == 0:
                alpha.append(float(c[i].real))
            w = w - c @ Q
            w = ts.flat(ts.gauge(ts.unflat(w)))
        del Q
        if i == nkry:
            break
        bnext = float(xp.linalg.norm(w))
        if bnext <= breakdown * max(scale, abs(alpha[-1]), 1e-300):
            break
        beta.append(bnext)
        store = backend.setrow(store, m, w / bnext)
        m += 1
    return np.array(alpha), np.array(beta), store, m


def poles(ts, b, nkry, a=None):
    """Discrete poles (w_alpha) and weights of <GS|A delta(w-H+E0) B|GS>
    projected to the tangent space, b = P B|GS> and a = P A^dagger|GS>
    (a = None means A = B^dagger). Also returns the Ritz vectors' Krylov
    coefficients and the basis (store, m), for the variance."""
    al, be, store, m = lanczos(ts, b, nkry)
    if m == 0:
        return np.zeros(0), np.zeros(0, dtype=complex), None, (store, m)
    Tm = np.diag(al[:m]) + np.diag(be[:m - 1], 1) + np.diag(be[:m - 1], -1)
    ws, U = np.linalg.eigh(Tm)
    bproj = U[0, :] * ts.norm(b)  # <alpha|b>
    if a is None:
        weights = np.abs(bproj) ** 2 + 0j
    else:
        ak = backend.to_host(_overlaps(store, ts.flat(a)))[:m]  # <v_i|a>
        aproj = U.conj().T @ ak  # <alpha|a>
        weights = aproj.conj() * bproj
    return ws, weights, U, (store, m)


def block_krylov(ts, starts, nkry, breakdown=1e-10):
    """Orthonormal basis of the block Krylov space of P(H-E0)P from the
    tangent vectors `starts`, nkry blocks beyond the starting one, and the
    projected Hamiltonian in it, Hq[i,j] = <q_i|P(H-E0)P|q_j>.

    The basis is kept flattened, one row per vector (ts.flat), so that
    every overlap is a matrix product. Every vector is orthogonalized twice
    against the whole basis, and a candidate is dropped when what is left
    of it is below `breakdown` times its norm before the
    orthogonalization, so the basis stops growing once the space reached
    from `starts` is invariant. The threshold is 1e-10 and not tighter:
    after two passes the remainder of a vector already in the space is
    rounding noise of relative size about 1e-13, and accepting it
    normalizes that noise, gauge violation included, up to a unit vector,
    which made P S P measurably non-Hermitian (0.4) on the basis.

    H Q is never stored: when P(H-E0)P|q_j> is computed, its overlaps with
    every vector already in the basis, which includes all q_i with i <= j,
    give Hq[:, j] above the diagonal, and Hermiticity gives the rest.
    Returns (store, count, Hq), the basis being the first `count` rows of
    the flattened store (all of it padded with zero rows on a device, see
    _rows) and Hq a host (count, count) array."""
    xp = backend.xp()
    N = ts.flat_size()
    # preallocated and filled row by row: a vstack per accepted vector
    # copies the whole basis each time
    dtype = np.result_type(*[x.dtype for w in starts for x in w],
                           ts.As[0].dtype, *[W.dtype for W in ts.Ws])
    nmax = min(N, len(starts) * (nkry + 1))
    store = xp.zeros((nmax, N), dtype=dtype)
    count = 0

    def add(block):
        nonlocal count, store
        for w in block:
            if count == nmax:
                break
            Q = _rows(store, count)
            w = ts.flat(ts.gauge(w))
            scale = float(xp.linalg.norm(w))
            for rep in range(2):
                w = w - _overlaps(Q, w) @ Q
                w = ts.flat(ts.gauge(ts.unflat(w)))
            del Q
            nw = float(xp.linalg.norm(w))
            if scale > 0.0 and nw > breakdown * scale:
                store = backend.setrow(store, count, w / nw)
                count += 1

    add(starts)
    done = 0
    Hu = np.zeros((nmax, nmax), dtype=complex)  # Hu[i, j] for i <= j
    for step in range(nkry + 1):
        n = count
        if n == done:
            break
        Hb = []
        for j in range(done, n):
            h = ts.apply_h(ts.unflat(store[j]))
            Hu[:n, j] = backend.to_host(_overlaps(_rows(store, n), ts.flat(h)))[:n]
            Hb.append(h)
        done = n
        if step == nkry:
            break
        add(Hb)
        del Hb
    upper = np.triu(Hu[:count, :count])
    Hq = upper + np.triu(upper, 1).conj().T
    Hq = 0.5 * (Hq + Hq.conj().T)  # real diagonal
    if not np.iscomplexobj(np.zeros(0, dtype=dtype)):
        Hq = Hq.real
    return store, count, Hq


def projected_spectrum(ts, ops, nkry):
    """Energies and operator matrices in the basis {|GS>} + the Ritz states
    of P(H-E0)P on the block Krylov space started from P O|GS> and
    P O^dagger|GS> for every operator in `ops`.

    ops: list of (W, Wdag) MPO-array pairs, the operator and its adjoint.
    Returns (e, mats): e[0] = 0 for the ground state and e[1:] the Ritz
    energies relative to E0, and for each operator the matrix
    mats[k][a,b] = <a|O_k|b> in that basis, in which the O_k|GS> column is
    exact and the tangent-tangent block is the projection P O_k P."""
    b = [ts.project_operator(W) for W, _ in ops]  # P O|GS>
    bd = [ts.project_operator(Wd) for _, Wd in ops]  # P O^dagger|GS>
    store, n, Hq = block_krylov(ts, b + bd, nkry)
    ws, U = np.linalg.eigh(Hq)
    Uh = U.conj().T
    e = np.concatenate([[0.0], ws])
    Q = _rows(store, n)

    def overlaps(v):  # <q_i|v> for the n basis vectors, on the host
        return backend.to_host(_overlaps(Q, v))[:n]

    mats = []
    for (W, Wd), bk, bdk in zip(ops, b, bd):
        op = ts.operator(W)
        M = np.zeros((n + 1, n + 1), dtype=complex)
        M[0, 0] = ts.expectation(op)
        M[1:, 0] = Uh @ overlaps(ts.flat(bk))
        M[0, 1:] = (Uh @ overlaps(ts.flat(bdk))).conj()
        G = np.zeros((n, n), dtype=complex)  # <q_i|O|q_j>
        for j in range(n):
            G[:, j] = overlaps(ts.flat(ts.apply(op, ts.unflat(store[j]))))
        M[1:, 1:] = Uh @ G @ U
        mats.append(M)
    return e, mats


def ritz_vector(basis, u):
    """The Ritz vector with Krylov coefficients u, as a tangent vector;
    basis is (ts, store, m)."""
    ts, store, m = basis
    xp = backend.xp()
    c = np.zeros(store.shape[0], dtype=np.result_type(u.dtype, store.dtype))
    c[:m] = u
    return ts.unflat(xp.asarray(c) @ store)


def correlator_poles(Ms, Ws, Wb, Wadag=None, nkry=100, variance=False):
    """The poles of <GS|A delta(w-H+E0) B|GS> from site arrays alone: Ms
    the ground-state MPS, Ws the Hamiltonian MPO, Wb the MPO of B and Wadag
    that of A^dagger (None when A = B^dagger). Shared by every backend that
    can export its site tensors. Returns (w, weights, variances, dim), the
    variances None unless variance=True."""
    ts = TangentSpace(Ms, Ws)
    b = ts.project_operator(Wb)
    a = None if Wadag is None else ts.project_operator(Wadag)
    ws, weights, U, (store, m) = poles(ts, b, nkry, a=a)
    var = None
    if variance and m > 0:
        var = np.array([ts.two_site_variance(ritz_vector((ts, store, m), U[:, k]))
                        for k in range(len(ws))])
    return ws, weights, var, ts.dimension()


def spectrum_from_arrays(Ms, Ws, ops, nkry):
    """projected_spectrum from site arrays: ops is a list of (MPO of O,
    MPO of O^dagger). Returns (e, mats, dim)."""
    ts = TangentSpace(Ms, Ws)
    e, mats = projected_spectrum(ts, ops, nkry)
    return e, mats, ts.dimension()
