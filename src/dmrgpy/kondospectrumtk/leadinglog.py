"""Leading-log resummation of the STM Kondo spectrum.

One-loop poor man's scaling with frequency-dependent couplings, in the
scheme of Rosch, Paaske, Kroha and Woelfle (cond-mat/0202404,
cond-mat/0408506 eq. 6) and Rosch, Costi, Paaske and Woelfle
(cond-mat/0301106), written for an impurity with an arbitrary spectrum in
its eigenbasis. Its expansion to first order in Jrho_s is exactly the
third order of conductance.py (same omega0, same Gamma0), and it sums the
leading logs (Jrho_s)^n ln^n to all orders. The derivation, the
references, the checks and what the method leaves out are in
docs/new_methods/kondo_leading_log/README.md.

Every vertex is stored as four impurity matrices A[mu], mu=0 (potential
scattering, electron identity) and mu=1,2,3 (exchange, Pauli matrix
sigma^mu on the electron), indexed [a,b] = <a|...|b> in the impurity
eigenbasis, the convention of KondoSpectrum.

As in the papers, the couplings on the right-hand side of the flow are
taken on shell, where the intermediate electron sits at the Fermi level:
the sample vertex G with its incoming electron at the Fermi level, and
the tip vertex element (x,y) at its own on-shell energy e_x - e_y (A_on).
These two sets of matrices flow together and close. The tip vertex at a
bias omega is then a quadrature over the cutoff of G times A_on, with the
cuts of omega.

How it is computed (_joint_flow). G, A_on and the tip vertex at every
bias are one system of equations, solved once, with nothing stored along
the flow. Three facts keep it cheap for an impurity of dimension d:
  - dI/dV reads the tip vertex only on the columns of the occupied
    states, and the quadrature of a column needs no other column, so the
    tip costs d^2 per bias and occupied state, not d^3;
  - a vertex is a 2d x 2d matrix (electron spin times impurity), and the
    contraction with _C is its matrix product, so the d^3 work of the
    flow of G and A_on is four 2d x 2d products per step;
  - the potential scattering of the tip does not flow: its direct and
    exchange terms cancel at every cutoff, so U is added at the end.
For a real Hamiltonian the block form of every vertex is a real matrix,
and the flow is run in real arithmetic (_block_is_real); a complex
eigenbasis, as with a field along y, takes the complex route.

For a Hamiltonian that conserves S^z the flow runs on the blocks of total
S^z (real route, CPU; sz_blocks=None, the default, at T=0 from SZ_BLOCKS_FROM
states on, True at any size and T, False never): an index (s,a) of a vertex
carries s + M_a, which every vertex conserves, so the four products and
the two sums of the tip vertex run on blocks (_SzBlocks,
_tip_bracket_blocks). The labels M are read off the tip spin
(_sz_labels); when a degenerate level mixes S^z sectors (every multiplet
at B=0, a level crossing in a field) the flow runs in the eigenbasis of
S^z_tot from the chain and the tip vertex is rotated back (_sz_rotation).

An impurity with no spin gives the potential scattering of the tip alone,
2*pi*T0^2*4*U^2 for dI/dV, without a warning: there is no exchange vertex
to flow, and the result is the second order.

At T>0 the cutoff kernel K(x,D) is replaced by its thermal average with
the weight of stepfunctions.FBuilder (thermalkernel.ThermalKernel, passed
to _joint_flow as an argument), so the first order is the thermal third
order; the populations are the Boltzmann ones and the step is thermal, as
in conductance.py. With no kernel argument the flow is the T=0 one, line
for line. README, section 8.
"""
import sys
import warnings

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp, DOP853
from scipy.optimize import brentq

from .conductance import _theta_raw, _occupied_states, P_CUT

# C[mu,nu,lam]: sigma^mu sigma^nu = sum_lam C[mu,nu,lam] sigma^lam
_C = np.zeros((4, 4, 4), dtype=complex)
_C[0, 0, 0] = 1.
for _a in range(1, 4):
    _C[0, _a, _a] = 1.
    _C[_a, 0, _a] = 1.
    _C[_a, _a, 0] = 1.
for (_i, _j, _k) in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]:
    _C[_i, _j, _k] = 1j
    _C[_j, _i, _k] = -1j
# the hole channel takes the electron matrices in the reversed order
_CR = np.transpose(_C, (1, 0, 2))

# g = SIGN*Jrho_s is the dimensionless sample exchange of the flow; the
# sign is fixed by the first-order match against conductance.py
# (tests/test_kondo_leading_log.py), and with it Jrho_s<0, the paper's
# antiferromagnetic (Kondo) sign, is the one whose coupling grows.
SIGN = 1.0

# a running coupling above this is reported with a warning: the
# leading-log result is then only qualitative
WARN_COUPLING = 0.5

_EPS = np.finfo(float).eps

# The running tip coupling compares the flowed with the bare vertex of a
# transition. A bare vertex below this fraction of the largest one out of
# the occupied states is replaced by that fraction. The flow gives a weak
# or forbidden transition a vertex of order g ln times the large matrix
# elements, which is under control and is many times its own bare value;
# measured against the floor, it is flagged only once it is comparable to
# the largest bare vertex.
BARE_FLOOR = 0.1

# _joint_flow tabulates, once per flow, the two D-independent arrays of the
# exchange sum of the tip vertex, sqrt(x^2+Gamma0^2) and omega0+|x| for every
# final state, bias and intermediate state: 16 n d^2 bytes for n distinct
# incoming energies. Above this size the sum computes them on the fly.
HOLE_TABLE_BYTES = 256e6

# a vertex is real in block form when the bare one is, to this relative
# precision (see _block_is_real)
REAL_TOL = 1e-12

# sz_blocks=None (the default) runs the flow on the S^z blocks at T=0 from
# this many states on, when the spectrum has them: below it the blocks cost
# more than they save (30 to 40% slower up to 16 states, even at 32, 1.7x
# faster at 64, T=0 and 61 biases), and at T>0 they gain nothing (0.93 to
# 1.11x on 5 to 7 spins at 1 and 2 K). The author's picks, 2026-11-14.
SZ_BLOCKS_FROM = 32

# an entry of the tip spin above this fraction of its largest one counts
# for the S^z labels of the states (_sz_labels); what lies below it off the
# blocks is the leak of the eigenbasis, 8e-13 on the 8-spin chain of
# docs/kondo_speedup_plan.md, and the blocked flow drops it
SZ_TOL = 1e-10


def _kernel(x, D, omega0, Gamma0):
    """K(x,D) = 1/(D+r_x) - 1/(D+omega0+|x|), r_x = sqrt(x^2+Gamma0^2):
    its integral over D from 0 to infinity is exactly F0(x) of
    stepfunctions.py, so the flow is the third order at first order."""
    ax = np.abs(x)
    return 1./(D + np.sqrt(ax*ax + Gamma0*Gamma0)) - 1./(D + omega0 + ax)


def _kernel_integral(x, Dlow, Dhigh, omega0, Gamma0):
    """The integral of K(x,D) over D from Dlow to Dhigh, in closed form."""
    ax = np.abs(x)
    r = np.sqrt(ax*ax + Gamma0*Gamma0)
    return (np.log((Dhigh + r)/(Dlow + r))
            - np.log((Dhigh + omega0 + ax)/(Dlow + omega0 + ax)))


def _block(A):
    """The vertex sum_mu sigma^mu (x) A[mu] as one 2d x 2d matrix, with the
    electron spin as the outer index; A: (...,4,d,d). The product of two
    vertices is then the matrix product, which is the contraction with _C
    in one 2d x 2d product instead of sixteen d x d ones."""
    d = A.shape[-1]
    B = np.empty(A.shape[:-3] + (2*d, 2*d), dtype=complex)
    B[..., :d, :d] = A[..., 0, :, :] + A[..., 3, :, :]
    B[..., :d, d:] = A[..., 1, :, :] - 1j*A[..., 2, :, :]
    B[..., d:, :d] = A[..., 1, :, :] + 1j*A[..., 2, :, :]
    B[..., d:, d:] = A[..., 0, :, :] - A[..., 3, :, :]
    return B


def _unblock(B):
    """The inverse of _block: (...,2d,2d) -> (...,4,d,d)."""
    d = B.shape[-1]//2
    A = np.empty(B.shape[:-2] + (4, d, d), dtype=complex)
    A[..., 0, :, :] = (B[..., :d, :d] + B[..., d:, d:])/2.
    A[..., 3, :, :] = (B[..., :d, :d] - B[..., d:, d:])/2.
    A[..., 1, :, :] = (B[..., :d, d:] + B[..., d:, :d])/2.
    A[..., 2, :, :] = 0.5j*(B[..., :d, d:] - B[..., d:, :d])
    return A


def _swap(B):
    """Transpose of the electron spin index only: the two off-diagonal
    d x d blocks change places."""
    d = B.shape[-1]//2
    out = B.copy()
    out[..., :d, d:] = B[..., d:, :d]
    out[..., d:, :d] = B[..., :d, d:]
    return out


def _bracket(Xb, Yb):
    """Direct minus exchange for two vertices in block form: X Y with the
    electron matrices in this order, minus the product with the electron
    matrices in the reversed order (the hole channel), where the impurity
    matrices keep the order Y X."""
    return np.matmul(Xb, Yb) - _swap(np.matmul(_swap(Yb), _swap(Xb)))


def _sz_labels(S, tol=SZ_TOL):
    """Labels M of the impurity states such that the tip spin changes M
    the way S^z, S^+ and S^- do: <a|S^z|b> is nonzero only for M_a = M_b,
    <a|S^+|b> only for M_a = M_b + 1 (and S^- the other way), an entry
    counting as nonzero above tol of the largest entry of S. Returns
    (M, comp), integer arrays, comp the connected component of each state
    under these entries (no product of the flow joins two components, and
    M starts at 0 in each), or None when no such labeling exists: the
    Hamiltonian does not conserve S^z (a transverse field), or a
    degenerate level of the eigenbasis mixes S^z sectors (a multiplet at
    B=0). For a Hamiltonian that conserves S^z, M is S^z_tot up to a
    constant per component."""
    S = np.asarray(S)
    d = S.shape[-1]
    big = tol*max(np.max(np.abs(S)), 1e-300)
    Sp = S[0] + 1j*S[1]
    moves = [(np.abs(S[2]) > big, 0), (np.abs(Sp) > big, 1),
             (np.abs(Sp.conj().T) > big, -1)]
    M = np.zeros(d, dtype=int)
    comp = np.full(d, -1)
    for start in range(d): # a spanning tree of each component
        if comp[start] >= 0:
            continue
        comp[start] = start
        stack = [start]
        while stack:
            b = stack.pop()
            for mask, shift in moves:
                for a in np.nonzero(mask[:, b] & (comp < 0))[0]:
                    comp[a], M[a] = start, M[b] + shift
                    stack.append(a)
    for mask, shift in moves: # every entry, not only the tree
        if np.any(mask & (M[:, None] != M[None, :] + shift)):
            return None
    return M, comp


def _sz_rotation(ks, cols, tol=SZ_TOL):
    """The eigenbasis of ks turned into one of definite S^z_tot, for a
    spectrum whose degenerate levels mix S^z sectors: at B=0 every
    multiplet, and in a field a level crossing (9 spins 1/2 with J = 1 meV
    and 0.3 meV along z hold four, a magnon at 5B below the ferromagnetic
    multiplet). U is the eigenbasis of diag(e)/W + 4 S^z_tot, W the spread
    of e, so it needs no energy tolerance: sectors 1 apart are 4 apart in
    it, more than the energies. It is accepted only when U^T diag(e) U is
    diagonal to tol of W, which is the measured test that H conserves
    S^z_tot (a transverse field fails it), and when U maps the columns
    cols (the occupied states) onto as many rotated ones. Returns
    (U, e, S, cols, Wc): the rotated energies, tip spin and columns, and
    Wc = U[cols, new cols], so that a tip vertex A in the rotated basis is
    U A Wc^T in the original one. None when ks has no chain to read
    S^z_tot from (ks.eigenbasis_operator, edkondo.KondoSpectrum) or when a
    test fails."""
    reader = getattr(ks, "eigenbasis_operator", None)
    chain = getattr(ks, "chain", None)
    if reader is None or chain is None or not hasattr(chain, "Sz"):
        return None
    Mm = np.asarray(reader(sum(chain.Sz)))
    if np.iscomplexobj(Mm):
        if np.max(np.abs(Mm.imag)) > tol*max(np.max(np.abs(Mm)), 1e-300):
            return None
        Mm = Mm.real
    Mm = 0.5*(Mm + Mm.T)
    e = np.asarray(ks.e, dtype=float)
    W = max(np.ptp(e), 1e-300)
    U = np.linalg.eigh(np.diag(e/W) + 4.0*Mm)[1]
    Eu = U.T @ (e[:, None]*U)
    eu = np.diag(Eu).copy()
    if np.max(np.abs(Eu - np.diag(eu))) > tol*W:
        return None
    cols = np.asarray(cols, dtype=int)
    weight = np.sum(U[cols]**2, axis=0) # of each rotated state on the columns
    new = np.nonzero(weight > 0.5)[0]
    if len(new) != len(cols) or np.max(np.minimum(weight, 1. - weight)) > tol:
        return None
    S = _spin_ops(ks)
    Sr = np.stack([U.T @ x @ U for x in S])
    return U, np.ascontiguousarray(eu), Sr, new, U[np.ix_(cols, new)]


class _SzBlocks():
    """The block structure of the flow for an impurity that conserves S^z
    (labels from _sz_labels). In block form the index (s,a), s the
    electron spin, carries s + M_a, which every vertex of the flow
    conserves, so the vertices are block diagonal in it; _swap of a vertex
    conserves M_a - s instead. The products of _bracket are then products
    of blocks, read and written through flat indices of the 2d x 2d
    matrices: the blocks of s + M_a for X Y, and for the hole channel the
    blocks of M_a - s mapped by _swap, which is its own inverse, so one map
    reads swap(Y) and writes swap(Z). What the dense flow holds off the
    blocks is the leak of the eigenbasis, below SZ_TOL of the largest
    entry; the bracket here is zero there."""

    def __init__(self, M, comp):
        d = len(M)
        n2 = 2*d
        two_s = np.concatenate([np.ones(d, dtype=int), -np.ones(d, dtype=int)])
        L, C = np.concatenate([M, M]), np.concatenate([comp, comp])
        self.d, self.M, self.comp = d, np.asarray(M), np.asarray(comp)
        self.direct, self.exchange = [], []
        for key, store in ((2*L + two_s, self.direct), (2*L - two_s, self.exchange)):
            _, inv = np.unique(np.stack([C, key]), axis=1, return_inverse=True)
            inv = np.ravel(inv)
            for g in range(inv.max() + 1):
                I = np.nonzero(inv == g)[0]
                if store is self.direct:
                    flat = I[:, None]*n2 + I[None, :]
                else: # (r,c) -> (s_c a_r, s_r b_c), the map of _swap
                    s_, a = I//d, I % d
                    flat = (s_[None, :]*d + a[:, None])*n2 + (s_[:, None]*d + a[None, :])
                store.append((np.ascontiguousarray(flat.ravel()), len(I)))

    def work(self):
        """The arithmetic of the blocked products over the dense ones."""
        return max(sum(n**3 for _, n in g) for g in (self.direct, self.exchange))/(2.*self.d)**3

    def hole_sum(self, cols, w, e, omega0, Gamma0):
        """The exchange sum of _tip_bracket on the S^z blocks at T=0, for
        the columns cols, a _BlockHoleSum. Q[a,(i,n,j),c] = sum_k A[k,a,n,c]
        G[i,c,k,b_j] needs s_k + M_a = s_n + M_c and s_i + M_c = s_k + M_b,
        all in one component, so a final state a has M_a within 1 of M_b,
        one or two of the four (i,n) and c in two or three adjacent
        sectors. In the order of (comp, M) those c are one range, so the
        sum runs on contiguous rows, and the tables of _hole_tables are
        kept on those ranges only (compact), when they fit in
        HOLE_TABLE_BYTES; otherwise the kernel is computed on the fly."""
        d, m = self.d, len(cols)
        perm = np.lexsort((self.M, self.comp)) # sorted position -> state
        M2, comp = 2*self.M[perm], self.comp[perm]
        sp = (1, -1)
        nz = np.zeros((d, 2, 2, m, d), dtype=bool) # [a,i,n,j,c], sorted a and c
        for j, b in enumerate(cols):
            one = (comp[:, None] == self.comp[b]) & (comp[None, :] == self.comp[b])
            for i in range(2):
                for n_ in range(2):
                    ok = np.zeros((d, d), dtype=bool)
                    for k in range(2):
                        ok |= ((sp[k] + M2[:, None] == sp[n_] + M2[None, :])
                               & (sp[i] + M2[None, :] == sp[k] + 2*self.M[b]))
                    nz[:, i, n_, j, :] = ok & one
        nz = nz.reshape(d, 4*m, d)
        rows = np.nonzero(nz.any(axis=(1, 2)))[0]
        cs = [np.nonzero(nz[a].any(axis=0))[0] for a in rows]
        lo = np.array([c[0] for c in cs], dtype=np.int64)
        hi = np.array([c[-1] + 1 for c in cs], dtype=np.int64)
        ll = [np.nonzero(nz[a].any(axis=1))[0] for a in rows]
        l_ptr = np.zeros(len(rows) + 1, dtype=np.int64)
        l_ptr[1:] = np.cumsum([len(x) for x in ll])
        l_idx = np.ascontiguousarray(np.concatenate(ll) if ll else np.zeros(0), dtype=np.int64)
        w = np.ascontiguousarray(w, dtype=float)
        es = np.ascontiguousarray(np.asarray(e, dtype=float)[perm])
        off = np.zeros(len(rows) + 1, dtype=np.int64)
        off[1:] = np.cumsum(len(w)*(hi - lo))
        tables = None
        if 16*off[-1] <= HOLE_TABLE_BYTES:
            R, Wt = np.empty(off[-1]), np.empty(off[-1])
            for r, a in enumerate(rows): # the expressions of _hole_tables
                ax = np.abs((w[:, None] - es[a]) + es[None, lo[r]:hi[r]])
                R[off[r]:off[r + 1]] = np.sqrt(ax*ax + Gamma0*Gamma0).ravel()
                Wt[off[r]:off[r + 1]] = (omega0 + ax).ravel()
            tables = (R, Wt, off)
        # the direct sum: P[j,(i,a,n),c] = sum_k G[i,a,k,c] A[k,c,n,b_j] needs
        # s_k + M_c = s_n + M_b and s_i + M_a = s_k + M_c
        direct = []
        M2o = 2*self.M
        for b in cols:
            ok = np.zeros((2, d, 2, d), dtype=bool) # [i,a,n,c]
            for i in range(2):
                for n_ in range(2):
                    for k in range(2):
                        ok[i, :, n_, :] |= ((sp[k] + M2o[None, :] == sp[n_] + M2o[b])
                                            & (sp[i] + M2o[:, None] == sp[k] + M2o[None, :]))
            ok &= ((self.comp[:, None] == self.comp[b])
                   & (self.comp[None, :] == self.comp[b]))[None, :, None, :]
            ok = ok.reshape(4*d, d)
            rr = np.nonzero(ok.any(axis=1))[0]
            cc = np.nonzero(ok.any(axis=0))[0]
            i_, a_, n_ = np.unravel_index(rr, (2, d, 2))
            direct.append((rr, cc, i_[:, None], a_[:, None], n_[:, None], cc[None, :]))
        return _BlockHoleSum(np.ascontiguousarray(perm, dtype=np.int64),
                             np.ascontiguousarray(rows, dtype=np.int64), lo, hi,
                             l_ptr, l_idx, tables, w, es, omega0, Gamma0, direct)

    def bracket(self, Xb, Yb):
        """_bracket on the blocks: X Y - swap(swap(Y) swap(X))."""
        x, y = Xb.ravel(), Yb.ravel()
        out = np.zeros(x.size, dtype=np.result_type(x, y))
        for I, n in self.direct:
            out[I] = (x[I].reshape(n, n) @ y[I].reshape(n, n)).ravel()
        for J, n in self.exchange:
            out[J] -= (y[J].reshape(n, n) @ x[J].reshape(n, n)).ravel()
        return out.reshape(Xb.shape)


class _BlockHoleSum():
    """The exchange sum of _tip_bracket on the S^z blocks
    (_SzBlocks.hole_sum): the states in the order perm, the final states
    rows, for each the intermediate states lo..hi and the columns
    l_idx[l_ptr[r]:l_ptr[r+1]], and the compact tables (R, Wt, off) of
    those ranges, or None to compute the kernel on the fly. Called as
    (Q, D, out) with out zeroed."""

    def __init__(self, perm, rows, lo, hi, l_ptr, l_idx, tables, w, es, omega0,
                 Gamma0, direct):
        self.direct = direct # per column: the rows (i,a,n) and the c of P
        self.perm, self.rows, self.lo, self.hi = perm, rows, lo, hi
        self.l_ptr, self.l_idx, self.tables = l_ptr, l_idx, tables
        self.w, self.es, self.omega0, self.Gamma0 = w, es, omega0, Gamma0

    def __call__(self, Q, D, out):
        p = self.perm
        Qp = np.ascontiguousarray(Q[p][:, :, p])
        pattern = (p, self.rows, self.lo, self.hi, self.l_ptr, self.l_idx)
        if self.tables is not None:
            _hole_sum_blocks_table(Qp, *self.tables, D, out, *pattern)
        else:
            _hole_sum_blocks(Qp, self.w, self.es, self.omega0, self.Gamma0, D, out,
                             *pattern)


def _rate(G, A, Kmat):
    """The bracket of the flow of a vertex with incoming energy omega,
    R[lam] = sum_c sum_{mu nu} C[mu,nu,lam] G[mu]_ac A[nu]_cb K[c,b]
                             - C[nu,mu,lam] A[mu]_ac G[nu]_cb K[a,c],
    with Kmat[x,y] = K(omega-(e_x-e_y)): the direct cut at
    omega-(e_c-e_b), the exchange cut at omega-(e_a-e_c).
    G: (4,d,d); A: (...,4,d,d), the on-shell vertex on the right-hand
    side; Kmat: (...,d,d)."""
    return _unblock(_bracket(_block(G), _block(A*Kmat[..., None, :, :])))


def _onshell_rate(G, A, Kfix):
    """The same bracket for the element (a,b) taken at its own on-shell
    energy omega = e_a - e_b: the direct cut becomes e_a - e_c and the
    exchange cut e_c - e_b, Kfix[x,y] = K(e_x - e_y)."""
    return _unblock(_bracket(_block(G*Kfix), _block(A)))


@njit(cache=True, fastmath=True)
def _hole_sum(Q, w, e, D, Dlow, omega0, Gamma0, out):
    """out[k,a,l] = sum_c Q[a,l,c] K(w[k]-(e[a]-e[c])), for real Q. The
    kernel is K(x,D), or its integral over the cutoff from Dlow to D when
    Dlow >= 0. This is the one sum of the flow with a kernel for every
    bias, final and intermediate state, so it never stores the kernel."""
    nk, d, nl = len(w), len(e), Q.shape[1]
    G2 = Gamma0*Gamma0
    row = np.empty(d)
    for a in range(d):
        for k in range(nk):
            u = w[k] - e[a]
            if Dlow < 0.:
                # K(x,D) on one denominator: the division is most of the cost
                for c in range(d):
                    ax = abs(u + e[c])
                    r = np.sqrt(ax*ax + G2)
                    row[c] = (omega0 + ax - r)/((D + r)*(D + omega0 + ax))
            else:
                for c in range(d):
                    ax = abs(u + e[c])
                    r = np.sqrt(ax*ax + G2)
                    row[c] = (np.log((D + r)/(Dlow + r))
                              - np.log((D + omega0 + ax)/(Dlow + omega0 + ax)))
            for l in range(nl):
                acc = 0.
                for c in range(d):
                    acc += Q[a, l, c]*row[c]
                out[k, a, l] = acc


@njit(cache=True, fastmath=True)
def _hole_sum_table(Q, R, Wt, D, out):
    """The same sum as _hole_sum for Dlow < 0 with the D-independent parts
    tabulated: R[a,k,c] = sqrt(x^2+Gamma0^2), Wt[a,k,c] = omega0+|x|,
    x = w[k]-e[a]+e[c], so that K = (Wt-R)/((D+R)(D+Wt)). One division per
    element, and the same value as _hole_sum to roundoff."""
    d, nk, _ = R.shape
    nl = Q.shape[1]
    row = np.empty(d)
    for a in range(d):
        for k in range(nk):
            for c in range(d):
                r = R[a, k, c]
                wt = Wt[a, k, c]
                row[c] = (wt - r)/((D + r)*(D + wt))
            for l in range(nl):
                acc = 0.
                for c in range(d):
                    acc += Q[a, l, c]*row[c]
                out[k, a, l] = acc


@njit(cache=True, fastmath=True)
def _hole_sum_blocks_table(Q, R, Wt, off, D, out, perm, rows, lo, hi, l_ptr, l_idx):
    """_hole_sum_table on the S^z blocks (_BlockHoleSum), with Q in the
    sorted order of the states: the final states rows[r] only, for each
    the intermediate states lo[r]..hi[r] and the columns
    l_idx[l_ptr[r]:l_ptr[r+1]], and the compact tables, row r at off[r]
    as (k, c). out[k, perm[a], l] is in the original order; every other
    entry of out is left as it comes, so it comes in zeroed. The loops run
    on slices with local offsets, which numba vectorizes (the same loops
    over c0..c1 were 1.6x slower on 8 spins)."""
    nk = out.shape[0]
    for r in range(len(rows)):
        a = rows[r]
        c0 = lo[r]
        nc = hi[r] - c0
        pa = perm[a]
        row = np.empty(nc)
        for k in range(nk):
            Rk = R[off[r] + k*nc:off[r] + (k + 1)*nc]
            Wk = Wt[off[r] + k*nc:off[r] + (k + 1)*nc]
            for t in range(nc):
                rr = Rk[t]
                wt = Wk[t]
                row[t] = (wt - rr)/((D + rr)*(D + wt))
            for u in range(l_ptr[r], l_ptr[r + 1]):
                l = l_idx[u]
                q = Q[a, l, c0:c0 + nc]
                acc = 0.
                for t in range(nc):
                    acc += q[t]*row[t]
                out[k, pa, l] = acc


@njit(cache=True, fastmath=True)
def _hole_sum_blocks(Q, w, es, omega0, Gamma0, D, out, perm, rows, lo, hi, l_ptr,
                     l_idx):
    """_hole_sum_blocks_table with the kernel computed on the fly, as
    _hole_sum does, es the energies in the sorted order."""
    nk = len(w)
    G2 = Gamma0*Gamma0
    for r in range(len(rows)):
        a = rows[r]
        c0 = lo[r]
        nc = hi[r] - c0
        pa = perm[a]
        row = np.empty(nc)
        ec = es[c0:c0 + nc]
        for k in range(nk):
            u0 = w[k] - es[a]
            for t in range(nc):
                ax = abs(u0 + ec[t])
                rr = np.sqrt(ax*ax + G2)
                row[t] = (omega0 + ax - rr)/((D + rr)*(D + omega0 + ax))
            for u in range(l_ptr[r], l_ptr[r + 1]):
                l = l_idx[u]
                q = Q[a, l, c0:c0 + nc]
                acc = 0.
                for t in range(nc):
                    acc += q[t]*row[t]
                out[k, pa, l] = acc


def _hole_tables(w, e, omega0, Gamma0):
    """The tables of _hole_sum_table, or None when they would take more
    than HOLE_TABLE_BYTES (two arrays of d*n*d doubles)."""
    n, d = len(w), len(e)
    if 16*n*d*d > HOLE_TABLE_BYTES:
        return None
    ax = np.abs((w[None, :, None] - e[:, None, None]) + e[None, None, :])
    return (np.ascontiguousarray(np.sqrt(ax*ax + Gamma0*Gamma0)),
            np.ascontiguousarray(omega0 + ax))


def _hole_call(Q, w, e, D, Dlow, omega0, Gamma0, out, tables):
    if tables is not None and Dlow < 0.:
        _hole_sum_table(Q, tables[0], tables[1], D, out)
    else:
        _hole_sum(Q, w, e, D, Dlow, omega0, Gamma0, out)


def _tip_bracket(Gb, Ab, cols, w, e, D, Dlow, omega0, Gamma0, tables=None,
                 kernel=None):
    """The bracket of _rate for the tip vertex at each incoming energy w[k],
    for the columns b in cols only (the initial states) and with the
    on-shell tip vertex Ab on the right-hand side, both vertices in block
    form. Returns (n,2,d,2,m): the block form with the columns cut,
    [k, electron out, a, electron in, column]. With Dlow >= 0 the kernel
    is integrated over the cutoff from Dlow to D (see _hole_sum).

    Only matrix-vector work: d^2 per bias and column, where the full
    vertex needs d^3. Real vertices (a real block form) take the real
    route: no imaginary halves, half the columns in the exchange sum.
    tables, from _hole_tables, is used by the exchange sum when Dlow < 0.
    kernel, a thermalkernel.ThermalKernel, replaces K and its integral in
    both sums (the tables are then not used)."""
    d, m, n = len(e), len(cols), len(w)
    G4 = Gb.reshape(2, d, 2, d)
    A4 = Ab.reshape(2, d, 2, d)
    # direct: sum_c G[a,c] A[c,b] K(w-(e_c-e_b)); P[j,(i,a,n),c]
    P = np.einsum('iakc,kcnj->jianc', G4, A4[:, :, :, cols]).reshape(m, 4*d, d)
    x = w[None, None, :] - (e[None, :, None] - e[cols][:, None, None]) # [j,c,k]
    if kernel is not None:
        Kp = kernel(x, D, Dlow)

        def hole_of(Qr, out):
            kernel.hole_sum(Qr, w, e, D, Dlow, out)
    else:
        if Dlow < 0.: Kp = _kernel(x, D, omega0, Gamma0)
        else: Kp = _kernel_integral(x, Dlow, D, omega0, Gamma0)

        def hole_of(Qr, out):
            _hole_call(Qr, w, e, D, Dlow, omega0, Gamma0, out, tables)
    if not (np.iscomplexobj(Gb) or np.iscomplexobj(Ab)):
        Q = np.einsum('kanc,ickj->ainjc', A4, G4[:, :, :, cols]).reshape(d, 4*m, d)
        hole = np.empty((n, d, 4*m))
        hole_of(np.ascontiguousarray(Q), hole)
        part = np.matmul(P, Kp).reshape(m, 2, d, 2, n).transpose(4, 1, 2, 3, 0)
        return part - hole.reshape(n, d, 2, 2, m).transpose(0, 2, 1, 3, 4)
    part = (np.matmul(np.ascontiguousarray(P.real), Kp)
            + 1j*np.matmul(np.ascontiguousarray(P.imag), Kp)) # [j,(i,a,n),k]
    part = part.reshape(m, 2, d, 2, n).transpose(4, 1, 2, 3, 0)
    # exchange: sum_c A[a,c] K(w-(e_a-e_c)) G[c,b]; Q[a,(i,n,j),c]
    Q = np.einsum('kanc,ickj->ainjc', A4, G4[:, :, :, cols]).reshape(d, 4*m, d)
    Qr = np.ascontiguousarray(np.concatenate([Q.real, Q.imag], axis=1))
    hole = np.empty((n, d, 8*m))
    hole_of(Qr, hole)
    hole = (hole[:, :, :4*m] + 1j*hole[:, :, 4*m:]).reshape(n, d, 2, 2, m)
    return part - hole.transpose(0, 2, 1, 3, 4)


def _tip_bracket_blocks(Gb, Ab, cols, w, e, D, omega0, Gamma0, sums):
    """_tip_bracket on the S^z blocks, for the real route at T=0 with
    Dlow < 0: the same sums, with the direct one on the rows (i,a,n) and
    the intermediate states c that the blocks allow (sums.direct), and the
    exchange one by sums, a _BlockHoleSum. The products are written out
    over k instead of einsum (4.7 ms for the two on 8 spins)."""
    d, m, n = len(e), len(cols), len(w)
    G4 = Gb.reshape(2, d, 2, d)
    A4 = Ab.reshape(2, d, 2, d)
    Ac = A4[:, :, :, cols] # [k,c,n,j]
    x = w[None, None, :] - (e[None, :, None] - e[cols][:, None, None]) # [j,c,k]
    Kp = _kernel(x, D, omega0, Gamma0)
    part = np.zeros((m, 4*d, n))
    for j in range(m):
        rr, cc, i_, a_, n_, c_ = sums.direct[j]
        P = G4[i_, a_, 0, c_]*Ac[0, c_, n_, j] + G4[i_, a_, 1, c_]*Ac[1, c_, n_, j]
        part[j, rr] = P @ Kp[j, cc]
    part = part.reshape(m, 2, d, 2, n).transpose(4, 1, 2, 3, 0)
    # exchange: Q[a,(i,n,j),c] = sum_k A[k,a,n,c] G[i,c,k,b_j]
    Gc = G4[:, :, :, cols].transpose(0, 3, 1, 2) # [i,j,c,k]
    Q = (A4[0][:, None, :, None, :]*Gc[None, :, None, :, :, 0]
         + A4[1][:, None, :, None, :]*Gc[None, :, None, :, :, 1]).reshape(d, 4*m, d)
    hole = np.zeros((n, d, 4*m))
    sums(Q, D, hole)
    return part - hole.reshape(n, d, 2, 2, m).transpose(0, 2, 1, 3, 4)


def _spin_ops(ks):
    return np.stack([np.asarray(ks.Sx, dtype=complex),
                     np.asarray(ks.Sy, dtype=complex),
                     np.asarray(ks.Sz, dtype=complex)])


def _bare_vertex(S, coupling, U):
    d = S.shape[-1]
    A = np.zeros((4, d, d), dtype=complex)
    A[0] = U*np.eye(d)
    A[1:] = coupling*S/2.
    return A


def _energy_range(ks, omega0, eVs, Gamma0):
    """D_max and D_min of the flow: the kernel tails beyond them are
    below ~1e-10 of the logs."""
    e = np.asarray(ks.e, dtype=float)
    scale = omega0 + (np.max(e) - np.min(e)) + (np.max(np.abs(eVs)) if len(eVs) else 0.)
    return 1e10*scale, 1e-10*Gamma0


def onshell_flow(ks, Jrho_s, U=0.0, omega0=20e-3, Gamma0=5e-6, rtol=1e-10,
                 atol=1e-14, pole_threshold=1.0, Dmax=None, Dmin=None):
    """Flow of the on-shell couplings from D=Dmax down to Dmin, with dense
    output: the reference for the couplings along the flow, used by the
    tests. dI/dV does not go through it (see _joint_flow), since the dense
    output is seven copies of the state per step. The state is the sample
    vertex G and the tip vertex A_on for bare potential scattering +U and
    -U (the second one is the s->t direction). Returns (sol, D_stop):
    sol.sol(ln D) gives the state, of shape (3,4,d,d) once viewed as
    complex, in the order G, A_on(+U), A_on(-U). D_stop is the cutoff where
    the sample coupling reached pole_threshold (None if it never did);
    below it the state is held at its value there. sol.max_coupling is
    the largest running sample coupling, in units of the bare S/2."""
    if not Gamma0 > 0.:
        raise ValueError("Gamma0 must be positive, it cuts the logs")
    S = _spin_ops(ks)
    d = S.shape[-1]
    e = np.asarray(ks.e, dtype=float)
    dE = e[:, None] - e[None, :]
    if Dmax is None or Dmin is None:
        Dmax0, Dmin0 = _energy_range(ks, omega0, np.zeros(1), Gamma0)
        Dmax = Dmax0 if Dmax is None else Dmax
        Dmin = Dmin0 if Dmin is None else Dmin
    y0 = np.stack([_bare_vertex(S, SIGN*Jrho_s, 0.), _bare_vertex(S, 1., U),
                   _bare_vertex(S, 1., -U)])
    smax = np.max(np.abs(S))/2.
    ref = pole_threshold*smax

    def rhs(s, y):
        D = np.exp(s)
        X = y.view(complex).reshape(3, 4, d, d)
        Kf = _kernel(dE, D, omega0, Gamma0)
        out = np.empty_like(X)
        out[0] = -D*_rate(X[0], X[0], Kf)
        out[1:] = -D*_onshell_rate(X[0], X[1:], Kf)
        return out.reshape(-1).view(float)

    def pole(s, y):
        G = y.view(complex).reshape(3, 4, d, d)[0]
        return ref - np.max(np.abs(G[1:]))
    pole.terminal = True
    pole.direction = -1

    sol = solve_ivp(rhs, (np.log(Dmax), np.log(Dmin)), y0.reshape(-1).view(float),
                    method="DOP853", rtol=rtol, atol=atol, dense_output=True,
                    events=pole)
    if not sol.success:
        raise RuntimeError("leading-log flow failed: "+str(sol.message))
    D_stop = None
    if len(sol.t_events[0]) > 0:
        D_stop = float(np.exp(sol.t_events[0][0]))
    sol.max_coupling = float(np.max([np.max(np.abs(y.view(complex).reshape(3, 4, d, d)[0, 1:]))
                                     for y in sol.y.T])/smax)
    return sol, D_stop


def sample_flow(ks, Jrho_s, **kwargs):
    """The sample vertex part of onshell_flow: (sol, D_stop) with
    sol.sol(ln D) viewed as complex and reshaped to (3,4,d,d), element 0."""
    return onshell_flow(ks, Jrho_s, **kwargs)


def _cut_columns(Tb):
    """_unblock for a vertex with its columns cut: (n,2,d,2,m) -> (n,4,d,m)."""
    return np.stack([(Tb[:, 0, :, 0] + Tb[:, 1, :, 1])/2.,
                     (Tb[:, 0, :, 1] + Tb[:, 1, :, 0])/2.,
                     0.5j*(Tb[:, 0, :, 1] - Tb[:, 1, :, 0]),
                     (Tb[:, 0, :, 0] - Tb[:, 1, :, 1])/2.], axis=1)


def _block_is_real(S, tol=REAL_TOL):
    """True when the block form of the bare vertex, [[A0+A3, A1-iA2],
    [A1+iA2, A0-A3]] for the stacked [0,Sx,Sy,Sz], is real to tol of its
    largest entry. It is for a real Hamiltonian (Sx and Sz real, Sy
    imaginary in a real eigenbasis), and then the kernel, the products and
    _swap keep every vertex of the flow real."""
    d = S.shape[-1]
    B = _block(np.concatenate([np.zeros((1, d, d), dtype=complex), S]))
    return bool(np.max(np.abs(B.imag)) <= tol*max(np.max(np.abs(B)), 1e-300))


def _on_device():
    """True when the array backend of pyitensor is JAX, so that the flow runs
    on the device. The backend module is read only when it is imported:
    set_backend("jax") imports it, so a process that never did is on NumPy,
    and importing the pyitensor package here would cost most of a second."""
    b = sys.modules.get(__name__.rsplit(".", 2)[0] + ".pyitensor.backend")
    return b is not None and b.backend_name() == "jax"


def _integrate(rhs, y0, Dmax, Dmin, sample_coupling, ref, rtol, atol):
    """DOP853 in ln D from Dmax to Dmin, stopped where sample_coupling(y)
    crosses ref from below (located inside the step by brentq on the dense
    output). Returns (y, D_stop, max_coupling): the state at the end or at
    the stop, the cutoff of the stop (None if there is none) and the
    largest sample coupling along the flow."""
    solver = DOP853(rhs, np.log(Dmax), y0, np.log(Dmin), rtol=rtol, atol=atol)
    g_old = sample_coupling(y0)
    max_coupling = g_old
    y, D_stop = solver.y, None
    while solver.status == "running":
        message = solver.step()
        if solver.status == "failed":
            raise RuntimeError("leading-log flow failed: "+str(message))
        y = solver.y
        g_new = sample_coupling(y)
        if g_old <= ref and g_new > ref:
            # the coupling crossed pole_threshold inside this step; the
            # crossing is strict, so a vertex that is zero never stops it
            dense = solver.dense_output()
            s_stop = brentq(lambda s: ref - sample_coupling(dense(s)),
                            solver.t_old, solver.t, xtol=4*_EPS, rtol=4*_EPS)
            y = dense(s_stop)
            D_stop = float(np.exp(s_stop))
            max_coupling = max(max_coupling, sample_coupling(y))
            break
        max_coupling = max(max_coupling, g_new)
        g_old = g_new
    return y, D_stop, max_coupling


def _joint_flow(ks, omegas, cols, Jrho_s, omega0, Gamma0, rtol, atol,
                pole_threshold, kernel=None, onshell_kernel=None, sz_blocks=None):
    """One flow from Dmax to Dmin of the sample vertex G, of the on-shell
    tip vertex A_on and of the tip vertex at each incoming energy
    omegas[k], the last one for the columns cols only, all at U=0.
    Returns the tip vertex, of shape (n,4,d,m), and the info of the flow.

    The tip vertex does not enter its own right-hand side, so it is a
    quadrature over G and A_on, and here it is integrated in the same
    solve with the same error control. Nothing is stored along the flow.
    If the sample coupling reaches pole_threshold at a cutoff D_stop, the
    couplings are held there, and what is left of the quadrature is the
    integral of the kernel from Dmin to D_stop, known in closed form.

    When the block form of the bare vertex is real (a real Hamiltonian)
    the state is kept in block form and in real arithmetic; otherwise the
    state is the complex (4,d,d) matrices of the vertices. The two give
    the same equations: the numbers differ by the integrator's step
    sequence, 1e-11 of dI/dV. max_coupling is in units of the bare S/2,
    and is 0 for an impurity with no spin.

    kernel, a thermalkernel.ThermalKernel, is the cutoff kernel of the tip
    vertex at a finite temperature, and onshell_kernel that of the flow of
    G and A_on; None is the T=0 kernel K(x,D) (see leading_log_dIdV).

    sz_blocks=True runs the flow on the S^z blocks of the spectrum: the
    products of G and A_on (_SzBlocks) and, at T=0, the two sums of the
    tip vertex (_tip_bracket_blocks), on the real route and on the CPU
    only. The labels come from the tip spin (_sz_labels), or, when a
    degenerate level mixes S^z sectors, from the eigenbasis of S^z_tot
    (_sz_rotation), in which the flow then runs; the tip vertex is rotated
    back. When neither works (H does not conserve S^z_tot) it warns and
    the flow is the dense one. None (the default) does the same at T=0
    (no kernel and no onshell_kernel) from SZ_BLOCKS_FROM states on and
    falls back without a warning; False is the dense flow. The info says
    which ran (sz_blocks, sz_rotated)."""
    if not Gamma0 > 0.:
        raise ValueError("Gamma0 must be positive, it cuts the logs")
    for k in (kernel, onshell_kernel):
        if k is not None and (k.omega0 != omega0 or k.Gamma0 != Gamma0):
            raise ValueError("the thermal kernel was built for another omega0 or Gamma0")
    if sz_blocks not in (None, True, False):
        raise ValueError("sz_blocks must be None, True or False, got %r" % (sz_blocks,))
    auto = sz_blocks is None # warns on a fallback only when asked for
    if auto: # at T=0 only: at T>0 the thermal tip sum is what costs
        sz_blocks = (len(ks.e) >= SZ_BLOCKS_FROM and kernel is None
                     and onshell_kernel is None)
    if _on_device():
        if sz_blocks and not auto:
            warnings.warn("leading-log flow: sz_blocks=True is for the CPU; the "
                          "flow on the device is the dense one", RuntimeWarning,
                          stacklevel=3)
        # the same flow with the state on the device (leadinglog_device.py)
        from .leadinglog_device import joint_flow
        return joint_flow(ks, omegas, cols, Jrho_s, omega0, Gamma0, rtol, atol,
                          pole_threshold, kernel, onshell_kernel)
    S = _spin_ops(ks)
    d = S.shape[-1]
    e = np.ascontiguousarray(ks.e, dtype=float)
    dE = e[:, None] - e[None, :]
    w = np.ascontiguousarray(omegas, dtype=float)
    cols = np.asarray(cols, dtype=int)
    n, m = len(w), len(cols)
    Dmax, Dmin = _energy_range(ks, omega0, w, Gamma0)
    smax = np.max(np.abs(S))/2.
    ref = pole_threshold*smax
    real = _block_is_real(S)
    bracket, pattern, rot = _bracket, None, None
    if sz_blocks:
        labels = _sz_labels(S) if real else None
        if labels is None and real:
            # a degenerate level that mixes S^z sectors (B=0, or a crossing):
            # the flow runs in the eigenbasis of S^z_tot, and the tip
            # vertex is rotated back at the end
            rot = _sz_rotation(ks, cols)
            labels = None if rot is None else _sz_labels(rot[2])
            if labels is None:
                rot = None
            else:
                U, e, S, cols, Wc = rot
                dE = e[:, None] - e[None, :]
                smax = np.max(np.abs(S))/2.
                ref = pole_threshold*smax
        if labels is None and not auto:
            warnings.warn("leading-log flow: sz_blocks=True, but %s; the flow is "
                          "the dense one" % ("the vertex is complex (the blocks "
                                             "are on the real route only)" if not real
                                             else "H does not conserve S^z_tot in "
                                             "this eigenbasis (a transverse field, "
                                             "or a spectrum with no chain to read "
                                             "S^z_tot from)"),
                          RuntimeWarning, stacklevel=3)
        if labels is not None:
            blocks = _SzBlocks(*labels)
            bracket = blocks.bracket
            if kernel is None:
                pattern = blocks.hole_sum(cols, w, e, omega0, Gamma0)
    # the blocks keep their own tables (_SzBlocks.hole_sum)
    tables = _hole_tables(w, e, omega0, Gamma0) if kernel is None and pattern is None \
        else None

    def kon(D): # the kernel of the on-shell flow at the energy differences
        if onshell_kernel is None:
            return _kernel(dE, D, omega0, Gamma0)
        return onshell_kernel(dE, D)
    A_bare = _bare_vertex(S, 1., 0.)
    Ab0 = _block(A_bare)
    T0 = Ab0.reshape(2, d, 2, d)[:, :, :, cols]
    if real:
        nb = 4*d*d
        y0 = np.concatenate([_block(_bare_vertex(S, SIGN*Jrho_s, 0.)).ravel(),
                             Ab0.ravel(),
                             np.broadcast_to(T0, (n, 2, d, 2, m)).ravel()]).real.copy()

        def couplings(y):
            return y[:2*nb].reshape(2, 2*d, 2*d)

        def rhs(s, y):
            D = np.exp(s)
            X = couplings(y)
            Gb, Ab = X[0], X[1]
            GKb = Gb*np.tile(kon(D), (2, 2))
            out = np.empty(len(y))
            out[:nb] = bracket(Gb, GKb).ravel()
            out[nb:2*nb] = bracket(GKb, Ab).ravel()
            if pattern is None:
                out[2*nb:] = _tip_bracket(Gb, Ab, cols, w, e, D, -1., omega0, Gamma0,
                                          tables, kernel).ravel()
            else:
                out[2*nb:] = _tip_bracket_blocks(Gb, Ab, cols, w, e, D, omega0, Gamma0,
                                                 pattern).ravel()
            out *= -D
            return out

        def sample_coupling(y):
            B = couplings(y)[0]
            return max(np.max(np.abs((B[:d, d:] + B[d:, :d])/2.)),
                       np.max(np.abs((B[:d, d:] - B[d:, :d])/2.)),
                       np.max(np.abs((B[:d, :d] - B[d:, d:])/2.)))

        y, D_stop, max_coupling = _integrate(rhs, y0, Dmax, Dmin, sample_coupling, ref,
                                              rtol, atol)
        Tb = y[2*nb:].reshape(n, 2, d, 2, m)
        if D_stop is not None and D_stop > Dmin:
            X = couplings(y)
            Tb = Tb + _tip_bracket(X[0], X[1], cols, w, e, D_stop, Dmin, omega0, Gamma0,
                                   kernel=kernel)
    else:
        nx = 8*d*d # complex entries of G and A_on
        y0 = np.concatenate([_bare_vertex(S, SIGN*Jrho_s, 0.).ravel(), A_bare.ravel(),
                             np.broadcast_to(T0, (n, 2, d, 2, m)).ravel()]).view(float)

        def couplings(y):
            return y.view(complex)[:nx].reshape(2, 4, d, d)

        def rhs(s, y):
            D = np.exp(s)
            X = couplings(y)
            Kf = kon(D)
            Gb, GKb, Ab = _block(X[0]), _block(X[0]*Kf), _block(X[1])
            out = np.empty(len(y)//2, dtype=complex)
            out[:nx//2] = _unblock(_bracket(Gb, GKb)).ravel()
            out[nx//2:nx] = _unblock(_bracket(GKb, Ab)).ravel()
            out[nx:] = _tip_bracket(Gb, Ab, cols, w, e, D, -1., omega0, Gamma0,
                                    tables, kernel).ravel()
            out *= -D
            return out.view(float)

        def sample_coupling(y):
            return np.max(np.abs(couplings(y)[0, 1:]))

        y, D_stop, max_coupling = _integrate(rhs, y0, Dmax, Dmin, sample_coupling, ref,
                                              rtol, atol)
        Tb = y.view(complex)[nx:].reshape(n, 2, d, 2, m)
        if D_stop is not None and D_stop > Dmin:
            X = couplings(y)
            Tb = Tb + _tip_bracket(_block(X[0]), _block(X[1]), cols, w, e, D_stop,
                                   Dmin, omega0, Gamma0, kernel=kernel)
    A = _cut_columns(Tb)
    if rot is not None:
        A = np.matmul(np.matmul(U, A), Wc.T) # back to the eigenbasis of ks
    return A, {"D_stop": D_stop,
               "max_coupling": float(max_coupling/(smax if smax > 0. else 1.)),
               "sz_blocks": bracket is not _bracket, "sz_rotated": rot is not None}


def tip_vertex(ks, omegas, Us, Jrho_s, omega0=20e-3, Gamma0=5e-6, rtol=1e-10,
               atol=1e-14, pole_threshold=1.0, cols=None, kernel=None,
               onshell_kernel=None, sz_blocks=None):
    """The tip vertex A(omega) at D -> 0 for each incoming energy
    omegas[k], with bare potential scattering Us[k]; shape (n,4,d,d), or
    (n,4,d,m) for the columns (initial states) cols only, which is what
    dI/dV needs and costs d^2 per bias and column instead of d^3. Also
    returns the flow's info: D_stop and max_coupling (see onshell_flow).

    The potential scattering of the tip does not flow: its direct and
    exchange terms cancel exactly at every cutoff. So the flow is solved
    once at U=0, and Us[k] is added to the mu=0 vertex at the end.
    kernel, onshell_kernel and sz_blocks: see _joint_flow (None is T=0)."""
    omegas = np.asarray(omegas, dtype=float)
    Us = np.broadcast_to(np.asarray(Us, dtype=float), omegas.shape)
    d = len(ks.e)
    cols = np.arange(d) if cols is None else np.asarray(cols, dtype=int)
    # a symmetric bias grid asks for each energy twice
    w, inverse = np.unique(omegas, return_inverse=True)
    A, info = _joint_flow(ks, w, cols, Jrho_s, omega0, Gamma0, rtol, atol,
                          pole_threshold, kernel, onshell_kernel, sz_blocks)
    A = A[inverse]
    A[:, 0, cols, np.arange(len(cols))] += Us[:, None]
    return A, info


def memory_model(d, m, real=False):
    """Bytes one flow of _joint_flow holds for an impurity of d states and
    m occupied columns, as (fixed, per incoming energy): the DOP853 solver
    keeps about 24 copies of the state, the couplings are 8 d^2 numbers and
    the tip vertex 4 d m per energy, and one right-hand side of
    _tip_bracket adds three arrays of 4 m d^2 numbers. The numbers are
    complex, or real when the flow takes the real route (real=True, see
    _block_is_real), which holds half the bytes. For 9 spins 1/2 at T=0
    (d=512, m=4, 1642 energies) the complex count gives 7.0 GB, and 6.67 GB
    were measured. On one flow of 241 biases (Triton, 2026-11-14, the peak
    above that of the process before the call): 9 spins at T=0, 0.65 GB on
    the real route (model 0.65) and 1.32 GB on the complex one (model
    1.30); 8 spins at 1 K with 55 occupied states, 1.78 GB on the real
    route (model 1.99, 61 biases)."""
    item = 8 if real else 16
    return item*24*8*d*d + 12*item*m*d*d, item*4*28*m*d


def usable_memory():
    """Bytes this process can still take, or None if it cannot tell: the
    smaller of MemAvailable in /proc/meminfo and, under a cgroup with a
    memory limit (a Slurm job, systemd-run -p MemoryMax), that limit minus
    what the cgroup uses. Every level of the cgroup path is read, v2
    (memory.max, memory.current) and v1 (memory.limit_in_bytes,
    memory.usage_in_bytes), and the tightest one counts."""
    import os
    free = []
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    free.append(int(line.split()[1])*1024)
    except (OSError, ValueError, IndexError):
        pass
    def read(path):
        try:
            with open(path) as f:
                s = f.read().strip()
            return None if s == "max" else int(s)
        except (OSError, ValueError):
            return None
    try:
        with open("/proc/self/cgroup") as f:
            lines = f.read().splitlines()
    except OSError:
        lines = []
    for line in lines:
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        if parts[0] == "0" and parts[1] == "": # v2
            root, limit, usage = "/sys/fs/cgroup", "memory.max", "memory.current"
        elif "memory" in parts[1].split(","): # v1
            root, limit, usage = "/sys/fs/cgroup/memory", "memory.limit_in_bytes", "memory.usage_in_bytes"
        else:
            continue
        path = parts[2].rstrip("/")
        while True:
            lim, use = read(root + path + "/" + limit), read(root + path + "/" + usage)
            # v1 writes "no limit" as a huge number, near 2^63
            if lim is not None and use is not None and lim < 2**60:
                free.append(max(lim - use, 0))
            if not path:
                break
            path = path.rsplit("/", 1)[0]
    return min(free) if free else None


def _ll_piece(ks, eVs, Jrho_s, U, T0, omega0, Gamma0, occ, rtol, atol,
              pole_threshold, kernel, onshell_kernel, open_cut, sz_blocks=None):
    """dI/dV and the running tip coupling at the biases eVs, from one
    flow; the body of leading_log_dIdV."""
    p = ks.p[occ]
    e = np.asarray(ks.e, dtype=float)
    n = len(eVs)
    # t->s at (eV, U) and its particle-hole conjugate s->t at (-eV, -U)
    omegas = np.concatenate([eVs, -eVs])
    Us = np.concatenate([np.full(n, U), np.full(n, -U)])
    # A[k,lam,f,i], i over the occupied states only
    A, info = tip_vertex(ks, omegas, Us, Jrho_s, omega0=omega0, Gamma0=Gamma0,
                         rtol=rtol, atol=atol, pole_threshold=pole_threshold,
                         cols=occ, kernel=kernel, onshell_kernel=onshell_kernel,
                         sz_blocks=sz_blocks)
    # W[k,f,i] = 4 sum_lam |A^lam_fi(omega_k)|^2
    W = 4.*np.sum(np.abs(A)**2, axis=1)
    eps_if = e[:, None] - e[occ][None, :] # [f,i] = e_f - e_i
    Th = _theta_raw(ks, omegas[:, None, None] - eps_if[None])
    g = np.einsum('i,kfi,kfi->k', p, W, Th)
    dIdV = 2*np.pi*T0**2*(g[:n] + g[n:])
    # the running tip coupling at each bias, over the open transitions. A
    # closed one does not enter dI/dV, and with many states there is
    # always one with a bare vertex near zero and a ratio far above 1.
    S = _spin_ops(ks)
    bare = np.sqrt(np.sum(np.abs(S[:, :, occ]/2.)**2, axis=0)) # [f,i]
    bare = np.maximum(bare, BARE_FLOOR*np.max(bare))
    flow = np.sqrt(np.sum(np.abs(A[:, 1:])**2, axis=1)) # [k,f,i]
    # a thermal step is never zero, so at T>0 a transition is open when
    # its weight p_i Theta is above open_cut of the largest population
    if ks.T > 0.:
        is_open = p[None, None, :]*Th > open_cut*np.max(p)
    else:
        is_open = Th > 0.
    # a zero bare vertex (a site with no spin) has nothing to compare with
    ratio = np.max(np.where(is_open, flow/np.where(bare > 0., bare, np.inf)[None], 0.),
                   axis=(1, 2))
    coupling = abs(Jrho_s)*np.maximum(ratio[:n], ratio[n:])
    return dIdV, coupling, info


def leading_log_dIdV(ks, eVs, Jrho_s, U=0.0, T0=1.0, omega0=20e-3,
                     Gamma0=5e-6, p_cut=P_CUT, rtol=1e-10, atol=1e-14,
                     pole_threshold=1.0, stop_margin=10.0, return_info=False,
                     onshell="thermal", open_cut=1e-2, mem_gb=None,
                     sz_blocks=None):
    """dI/dV with the leading logs summed to all orders, in the units of
    conductance.py (2*pi*e^2*T0^2/hbar). At first order in Jrho_s it is
    second_order_dIdV + third_order_kondo_dIdV + third_order_potential_dIdV
    exactly, at the temperature ks.T; see
    docs/new_methods/kondo_leading_log/README.md.

    At T>0 (ks.T > 0, with ks.p the Boltzmann weights) the cutoff kernel is
    its thermal average with the weight of stepfunctions.FBuilder
    (thermalkernel.ThermalKernel), so the first order is the thermal third
    order, and the amplitude is averaged and then squared and multiplied by
    the thermal step, as conductance.py does at third order. onshell sets
    the kernel of the flow of the sample vertex and of the on-shell tip
    vertex: "thermal" (the default) uses the thermal one there too, so the
    sample coupling stops growing at D ~ kT; "T0" keeps the T=0 kernel
    there. The two differ at second order in Jrho_s. At T=0 onshell has no
    effect.

    Two checks of validity, per bias:
      - the running tip coupling, |Jrho_s| times the largest ratio of the
        flowed to the bare spin vertex over the transitions out of occupied
        states that are open at this bias, the ones that enter dI/dV
        there, with the bare vertex bounded below by BARE_FLOOR. Above
        WARN_COUPLING it warns, at pole_threshold the bias is NaN. At T>0 a
        transition is open when p_i Theta is above open_cut of the largest
        p, since a thermal step is never zero.
      - if the sample coupling reaches pole_threshold during the flow (the
        Kondo scale is above Gamma0 and kT), the flow stops at that cutoff
        D_stop, and every bias within stop_margin*D_stop of zero or of any
        energy difference e_x - e_y of the spectrum is NaN, with a warning.

    A site with no spin on any state of the impurity has no exchange
    vertex: the flow does nothing, dI/dV is the potential scattering
    2*pi*T0^2*4*U^2 of the second order, and there is no warning.

    mem_gb: None solves every bias in one flow, and on the CPU warns when
    memory_model puts that flow above usable_memory(). A number cuts the biases
    in chunks so that one flow stays below about mem_gb (memory_model, or
    leadinglog_device.memory_model on the device), at
    the cost of solving the couplings again in each chunk; it raises when
    not even 8 biases fit. On a GPU the budget is at most
    leadinglog_device.allocator_budget(), a share of the limit of JAX's
    allocator: a flow that does not fit under it is cut in chunks also
    with mem_gb=None, and a larger mem_gb warns and is lowered to it. At T>0 the occupied states are the ones within
    ln(1/p_cut) kT of the ground state, and the tip part costs d^2 per
    energy and per occupied state.

    sz_blocks=None (the default) runs the flow on the blocks of total S^z
    at T=0 from SZ_BLOCKS_FROM states on when H conserves S^z, and the
    dense flow otherwise, without a warning; True asks for the blocks at
    any size and temperature and warns when they cannot run, False is the
    dense flow. On the blocks
    (see _joint_flow, CPU and real route only) run the products of G and
    A_on, and at T=0 the two sums of the tip vertex. A degenerate level
    that mixes S^z sectors (B=0, a level crossing) is rotated to the
    eigenbasis of S^z_tot, read from ks.chain; with a transverse field
    the flow is the dense one. The blocked
    flow drops what the eigenbasis leaks off the blocks, below SZ_TOL of
    the largest entry, and a rotated one moves dI/dV at the rtol of the
    flow (below 4e-11 of the peak, measured).

    return_info=True also returns a dict with D_stop, max_coupling (the
    largest running sample coupling, Jrho_s itself at the start), coupling
    (the running tip coupling at each bias), n_occ, n_chunks, sz_blocks
    (True when every flow ran on the S^z blocks) and sz_rotated (True when
    one ran in the eigenbasis of S^z_tot)."""
    if onshell not in ("thermal", "T0"):
        raise ValueError("onshell must be \"thermal\" or \"T0\", got %r" % (onshell,))
    if not ks.T >= 0.:
        raise ValueError("T must be >= 0, got %r" % (ks.T,))
    eVs = np.asarray(eVs, dtype=float)
    occ = _occupied_states(ks, p_cut)
    e = np.asarray(ks.e, dtype=float)
    n = len(eVs)
    kernel = onshell_kernel = None
    if ks.T > 0.:
        from .thermalkernel import ThermalKernel
        span = (np.max(e) - np.min(e)) + (np.max(np.abs(eVs)) if n else 0.)
        kernel = ThermalKernel(ks.kB*ks.T, omega0, Gamma0, max(1.05*span, 10*omega0))
        if onshell == "thermal":
            onshell_kernel = kernel
    budget = None if mem_gb is None else mem_gb*1e9
    cap = None
    if _on_device():
        from . import leadinglog_device as device
        fixed, per_w = device.memory_model(len(e), len(occ), _block_is_real(_spin_ops(ks)))
        cap = device.allocator_budget()
        if cap is not None:
            if mem_gb is None and fixed + 2*n*per_w > cap:
                budget = cap
            elif mem_gb is not None and budget > cap:
                warnings.warn("leading_log_dIdV: mem_gb=%g is above the %.3g GB the "
                              "device allocator gives (XLA_PYTHON_CLIENT_MEM_FRACTION, "
                              "0.75 of the card by default); the biases are cut to fit "
                              "it" % (mem_gb, cap/1e9), stacklevel=2)
                budget = cap
    else:
        fixed, per_w = memory_model(len(e), len(occ), _block_is_real(_spin_ops(ks)))
        if budget is None:
            need, free = fixed + 2*n*per_w, usable_memory()
            if free is not None and need > free:
                warnings.warn("leading_log_dIdV: one flow of %d biases at %d occupied "
                              "states of %d needs about %.3g GB (memory_model), above "
                              "the %.3g GB this process can still take; pass mem_gb= "
                              "to cut the biases in chunks" % (n, len(occ), len(e),
                                                               need/1e9, free/1e9),
                              RuntimeWarning, stacklevel=2)
    if budget is None:
        pieces = [np.arange(n)]
    else:
        nb = int((budget - fixed)//(2*per_w))
        if nb < 8:
            if budget == cap:
                raise ValueError("leading_log_dIdV: %d occupied states of %d need %.3g GB "
                                 "for 8 biases, above the %.3g GB the device allocator "
                                 "gives; raise XLA_PYTHON_CLIENT_MEM_FRACTION before Python "
                                 "starts, or p_cut to keep fewer states"
                                 % (len(occ), len(e), (fixed + 16*per_w)/1e9, cap/1e9))
            raise ValueError("leading_log_dIdV: %d occupied states of %d need %.2g GB "
                             "for 8 biases, above mem_gb=%g; raise mem_gb, or p_cut to "
                             "keep fewer states" % (len(occ), len(e),
                                                    (fixed + 16*per_w)/1e9, mem_gb))
        pieces = [np.arange(i, min(n, i + nb)) for i in range(0, n, nb)] or [np.arange(0)]
    dIdV, coupling = np.empty(n), np.empty(n)
    D_stop, max_coupling = None, 0.
    sz_blocked, sz_rotated = bool(pieces), False
    for idx in pieces:
        dIdV[idx], coupling[idx], piece = _ll_piece(
            ks, eVs[idx], Jrho_s, U, T0, omega0, Gamma0, occ, rtol, atol,
            pole_threshold, kernel, onshell_kernel, open_cut, sz_blocks)
        max_coupling = max(max_coupling, piece["max_coupling"])
        sz_blocked = sz_blocked and piece.get("sz_blocks", False)
        sz_rotated = sz_rotated or piece.get("sz_rotated", False)
        if piece["D_stop"] is not None:
            D_stop = piece["D_stop"] if D_stop is None else max(D_stop, piece["D_stop"])
    info = {"D_stop": D_stop, "max_coupling": max_coupling, "coupling": coupling,
            "n_occ": len(occ), "n_chunks": len(pieces), "sz_blocks": sz_blocked,
            "sz_rotated": sz_rotated}
    bad = coupling >= pole_threshold
    if np.any(coupling >= WARN_COUPLING):
        warnings.warn("leading_log_dIdV: the running coupling reaches %.2g at "
                      "%d bias points; the leading-log result is not "
                      "controlled there%s" % (np.max(coupling), int(np.sum(coupling >= WARN_COUPLING)),
                                              " (NaN where it passes %g)" % pole_threshold
                                              if np.any(bad) else ""))
    if D_stop is not None:
        thr = np.unique(np.concatenate([[0.], np.abs(e[:, None] - e[None, :]).ravel()]))
        near = np.min(np.abs(np.abs(eVs)[:, None] - thr[None]), axis=1) < stop_margin*D_stop
        if np.any(near):
            warnings.warn("leading_log_dIdV: the sample coupling reached order "
                          "one at D = %.3g (the Kondo scale is above Gamma0%s); "
                          "%d bias points within %g*D of an energy difference "
                          "of the spectrum are not leading-log results and are "
                          "returned as NaN" % (D_stop, " and kT" if ks.T > 0. else "",
                                               int(np.sum(near)), stop_margin))
        bad = bad | near
    dIdV = np.where(bad, np.nan, dIdV)
    if return_info:
        return dIdV, info
    return dIdV
