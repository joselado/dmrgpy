"""Elementwise (Hadamard) products of MPS by recursive sketched interpolation.

Reference: Z. Meng, Y. Khoo, J. Li and E. M. Stoudenmire, "Recursive
Sketched Interpolation: Efficient Hadamard Products of Tensor Trains",
arXiv:2602.17974. The authors publish a Python version
(github.com/zmeng137/Recursive-Sketched-Interpolation) that carries no
licence, so nothing here is taken from it; this module is written from the
paper's Sec. 3.

The product G(s1..sn) = T1(s1..sn) T2(s1..sn) of two MPS of bond dimension
chi has an exact MPS of bond dimension chi^2 (the core-wise Kronecker
product), and compressing it back costs O(chi^4) or more. RSI builds a
compressed G directly, one core per iteration, at O(n d^2 chi^3):

1. Sketch. At iteration j the sites j, j+1 stay open and the tail
   j+2..n of every input is contracted with k random product states (the
   "one-cluster" sketch: one Gaussian k x d matrix per site, the same for
   every input). The tail sketches are built once, in one backward sweep,
   since the tail cores are never touched by the iterations.
2. Interpolate. The sketched inputs are multiplied elementwise and the
   result, matricized as (left pivot, s_j) x (s_{j+1}, sketch column), is
   split by a row interpolative decomposition M ~ X M[I,:]. X is the new
   core of G; the skeleton M[I,:] is discarded, except at the last bond.
3. Re-interpolate. Every input is restricted to the pivot rows I of its
   own (left pivot, s_j) two-site tensor, which aligns its left basis
   with the one G now uses. Only slicing, no decomposition.

When the tail is so short that its full dimension is at most k, the exact
tail replaces the sketch and the remaining iterations are exact
interpolations (Sec. 3.1.4).

Where this departs from the paper, and why:

- The pivots come from a full-pivoting partial LU, as in the paper, but
  the stopping tolerance is relative to the first (largest) pivot, the
  rule every threshold of this codebase follows; the paper's is absolute.
- Each sketch column and each re-interpolated input is rescaled as it is
  built, and the scale is put back at the end. Amplitudes of an n-site
  state are of order d^(-n/2), so a product of two underflows past a
  thousand sites without it. The interpolation matrix X is invariant
  under both rescalings in exact arithmetic (a column scaling of M and an
  overall scale of a row both cancel in X = M[:,J] M[I,J]^-1).
- `func` (Sec. 6, the nonlinear-map extension) is applied to the sketched
  product exactly where the paper applies it. That is a heuristic outside
  the exact-tail iterations, as the paper itself says; no rescaling is
  done when it is set, since a nonlinear map does not commute with one.

Every array goes through the backend namespace (`backend.xp()`), so the
kernel runs on a device as well; only pivot positions and pivot values
come back to the host.
"""

import numpy as np

from . import backend as bk


def _row_interpolation(M, maxdim, tol):
    """Row interpolative decomposition M ~ X M[I,:] of an (m, n) matrix.

    The pivots are those of a partially rank-revealing LU with full
    pivoting: at each step the largest entry of the Schur complement picks
    a row and a column, and the step stops after `maxdim` pivots or when
    the largest remaining entry is at most `tol` times the first pivot.
    With L the (m, r) elimination columns, L[I,:] is unit lower triangular
    and X = L L[I,:]^-1, so X[I,:] is the identity and the pivot rows are
    reproduced exactly. Returns (I, X), I a list of row positions. A zero
    matrix gives one pivot and X = 0, so the bond keeps dimension 1."""
    xp = bk.xp()
    m, n = M.shape
    S = M
    rowmask = xp.ones((m,), dtype=float)
    colmask = xp.ones((n,), dtype=float)
    rows, cols = [], []
    piv0 = None
    for _ in range(min(m, n, int(maxdim))):
        a = xp.abs(S) * rowmask[:, None] * colmask[None, :]
        flat = int(np.asarray(xp.argmax(a)))
        i, j = divmod(flat, n)
        p = bk.scalar(S[i, j])
        if piv0 is None:
            if p == 0.0:
                break
            piv0 = abs(p)
        elif abs(p) <= tol * piv0:
            break
        col = S[:, j] / p
        S = S - xp.outer(col, S[i, :])
        rows.append(i)
        cols.append(col)
        rowmask = rowmask * (xp.arange(m) != i)
        colmask = colmask * (xp.arange(n) != j)
    if not rows:
        return [0], xp.zeros((m, 1), dtype=complex)
    L = xp.stack(cols, axis=1)
    L11 = L[xp.asarray(rows), :]
    X = xp.linalg.solve(L11.T, L.T).T
    return rows, X


def _rescale_columns(K):
    """K with every column divided by its largest absolute entry (a zero
    column is left as it is)."""
    xp = bk.xp()
    s = xp.max(xp.abs(K), axis=0)
    s = xp.where(s > 0, s, 1.0)
    return K / s[None, :]


def _tail_tables(tt, omegas, k):
    """The tail matrices of one input, by site m (0-based) at which the
    tail starts: exact[m] is the full contraction of sites m..n-1 as a
    (chi_{m-1}, prod d) matrix, kept only while prod d <= k, and
    sketch[m] its one-cluster sketch, (chi_{m-1}, k), built for every
    m >= 2 by one backward Khatri-Rao sweep. exact[n] = [[1]]."""
    xp = bk.xp()
    n = len(tt)
    exact = {n: xp.ones((1, 1), dtype=complex)}
    sketch = {n: xp.ones((1, k), dtype=complex)}
    for m in range(n - 1, 1, -1):
        A = tt[m]
        if m + 1 in exact and A.shape[1] * exact[m + 1].shape[1] <= k:
            E = xp.tensordot(A, exact[m + 1], axes=([2], [0]))
            exact[m] = E.reshape(A.shape[0], -1)
        B = xp.einsum("asb,ls->alb", A, omegas[m])
        K = xp.einsum("alb,bl->al", B, sketch[m + 1])
        sketch[m] = _rescale_columns(K)
    return exact, sketch


def hadamard_arrays(tts, maxdim, tol=1e-12, oversampling=8, seed=0,
                    func=None):
    """The elementwise product of several MPS, as a list of site arrays.

    tts: a list of MPS, each a list of (chi_left, d, chi_right) arrays with
        length-1 links at the open ends, all with the same physical
        dimensions. Conjugate an input before passing it in to get, for
        example, |psi(s)|^2 = conj(psi(s)) psi(s).
    maxdim: the largest bond dimension of the result.
    tol: interpolation tolerance, relative to the largest entry of each
        interpolated matrix.
    oversampling: the sketch dimension is k = ceil(maxdim/d) + oversampling
        (the paper's Eq. (9), with d the smallest physical dimension).
    seed: seed of the Gaussian sketch matrices.
    func: an elementwise map applied to the product (Sec. 6 of the paper),
        for example lambda x: np.maximum(x.real, 0). It must act on an
        array of the backend namespace. A heuristic outside the exact-tail
        iterations.

    Returns the site arrays of the result, in interpolation gauge (not
    canonical). The cost is O(n d^2 chi^3), against O(chi^4) at least for
    the exact Kronecker product followed by compression."""
    xp = bk.xp()
    ntt = len(tts)
    if ntt < 1:
        raise ValueError("hadamard_arrays needs at least one MPS")
    n = len(tts[0])
    dims = [int(a.shape[1]) for a in tts[0]]
    for tt in tts[1:]:
        if len(tt) != n or [int(a.shape[1]) for a in tt] != dims:
            raise ValueError("hadamard_arrays: the MPS do not have the same "
                             "sites")
    tts = [[bk.asarray(a) for a in tt] for tt in tts]
    maxdim = int(maxdim)
    if maxdim < 1:
        raise ValueError("hadamard_arrays: maxdim must be at least 1")

    def combine(arrs):
        G = arrs[0]
        for g in arrs[1:]:
            G = G * g
        return func(G) if func is not None else G

    if n == 1:
        return [combine([tt[0] for tt in tts])]

    k = int(np.ceil(maxdim / min(dims))) + int(oversampling)
    rng = np.random.default_rng(seed)
    omegas = {m: bk.asarray(rng.standard_normal((k, dims[m])))
              for m in range(2, n)}
    tables = [_tail_tables(tt, omegas, k) for tt in tts]

    rescale = func is None
    logscale = 0.0
    left = [tt[0] for tt in tts]
    cores = []
    r = 1
    for j in range(n - 1):
        m = j + 2
        thetas, mats = [], []
        for i in range(ntt):
            th = xp.tensordot(left[i], tts[i][j + 1], axes=([2], [0]))
            exact, sketch = tables[i]
            T = exact[m] if m in exact else sketch[m]
            thetas.append(th)
            mats.append(xp.tensordot(th, T, axes=([3], [0])))
        G = combine(mats)
        c = G.shape[3]
        M = G.reshape(r * dims[j], dims[j + 1] * c)
        rows, X = _row_interpolation(M, maxdim, tol)
        rn = len(rows)
        cores.append(X.reshape(r, dims[j], rn))
        idx = xp.asarray(rows)
        if j == n - 2:
            cores.append(M[idx, :].reshape(rn, dims[n - 1], 1))
            break
        left = []
        for th in thetas:
            L = th.reshape(r * dims[j], dims[j + 1], th.shape[3])[idx]
            if rescale:
                s = float(np.asarray(xp.max(xp.abs(L))))
                if s > 0.0:
                    L = L / s
                    logscale += np.log(s)
            left.append(L)
        r = rn
    if rescale and logscale != 0.0:
        f = np.exp(logscale / n)
        cores = [C * f for C in cores]
    return cores


def kronecker_arrays(tts):
    """The exact elementwise product of several MPS by the core-wise
    Kronecker product: bond dimension the product of the inputs' ones.
    The reference RSI is measured against, and the route for small
    inputs where chi^2 is affordable."""
    xp = bk.xp()
    out = [bk.asarray(a) for a in tts[0]]
    for tt in tts[1:]:
        new = []
        for a, b in zip(out, tt):
            b = bk.asarray(b)
            c = xp.einsum("asb,csd->acsbd", a, b)
            new.append(c.reshape(a.shape[0] * b.shape[0], a.shape[1],
                                 a.shape[2] * b.shape[2]))
        out = new
    return out


def dense_vector(arrays):
    """The full amplitude vector of an MPS given by its site arrays, site 1
    the slowest index. Exponential in the number of sites; for tests."""
    xp = bk.xp()
    v = arrays[0]
    for A in arrays[1:]:
        v = xp.tensordot(v, A, axes=([v.ndim - 1], [0]))
    return np.asarray(v).reshape(-1)
