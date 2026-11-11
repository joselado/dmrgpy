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

An impurity with no spin gives the potential scattering of the tip alone,
2*pi*T0^2*4*U^2 for dI/dV, without a warning: there is no exchange vertex
to flow, and the result is the second order.
"""
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


def _tip_bracket(Gb, Ab, cols, w, e, D, Dlow, omega0, Gamma0, tables=None):
    """The bracket of _rate for the tip vertex at each incoming energy w[k],
    for the columns b in cols only (the initial states) and with the
    on-shell tip vertex Ab on the right-hand side, both vertices in block
    form. Returns (n,2,d,2,m): the block form with the columns cut,
    [k, electron out, a, electron in, column]. With Dlow >= 0 the kernel
    is integrated over the cutoff from Dlow to D (see _hole_sum).

    Only matrix-vector work: d^2 per bias and column, where the full
    vertex needs d^3. Real vertices (a real block form) take the real
    route: no imaginary halves, half the columns in the exchange sum.
    tables, from _hole_tables, is used by the exchange sum when Dlow < 0."""
    d, m, n = len(e), len(cols), len(w)
    G4 = Gb.reshape(2, d, 2, d)
    A4 = Ab.reshape(2, d, 2, d)
    # direct: sum_c G[a,c] A[c,b] K(w-(e_c-e_b)); P[j,(i,a,n),c]
    P = np.einsum('iakc,kcnj->jianc', G4, A4[:, :, :, cols]).reshape(m, 4*d, d)
    x = w[None, None, :] - (e[None, :, None] - e[cols][:, None, None]) # [j,c,k]
    if Dlow < 0.: Kp = _kernel(x, D, omega0, Gamma0)
    else: Kp = _kernel_integral(x, Dlow, D, omega0, Gamma0)
    if not (np.iscomplexobj(Gb) or np.iscomplexobj(Ab)):
        Q = np.einsum('kanc,ickj->ainjc', A4, G4[:, :, :, cols]).reshape(d, 4*m, d)
        hole = np.empty((n, d, 4*m))
        _hole_call(np.ascontiguousarray(Q), w, e, D, Dlow, omega0, Gamma0, hole, tables)
        part = np.matmul(P, Kp).reshape(m, 2, d, 2, n).transpose(4, 1, 2, 3, 0)
        return part - hole.reshape(n, d, 2, 2, m).transpose(0, 2, 1, 3, 4)
    part = (np.matmul(np.ascontiguousarray(P.real), Kp)
            + 1j*np.matmul(np.ascontiguousarray(P.imag), Kp)) # [j,(i,a,n),k]
    part = part.reshape(m, 2, d, 2, n).transpose(4, 1, 2, 3, 0)
    # exchange: sum_c A[a,c] K(w-(e_a-e_c)) G[c,b]; Q[a,(i,n,j),c]
    Q = np.einsum('kanc,ickj->ainjc', A4, G4[:, :, :, cols]).reshape(d, 4*m, d)
    Qr = np.ascontiguousarray(np.concatenate([Q.real, Q.imag], axis=1))
    hole = np.empty((n, d, 8*m))
    _hole_call(Qr, w, e, D, Dlow, omega0, Gamma0, hole, tables)
    hole = (hole[:, :, :4*m] + 1j*hole[:, :, 4*m:]).reshape(n, d, 2, 2, m)
    return part - hole.transpose(0, 2, 1, 3, 4)


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
                pole_threshold):
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
    and is 0 for an impurity with no spin."""
    if not Gamma0 > 0.:
        raise ValueError("Gamma0 must be positive, it cuts the logs")
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
    tables = _hole_tables(w, e, omega0, Gamma0)
    real = _block_is_real(S)
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
            GKb = Gb*np.tile(_kernel(dE, D, omega0, Gamma0), (2, 2))
            out = np.empty(len(y))
            out[:nb] = _bracket(Gb, GKb).ravel()
            out[nb:2*nb] = _bracket(GKb, Ab).ravel()
            out[2*nb:] = _tip_bracket(Gb, Ab, cols, w, e, D, -1., omega0, Gamma0,
                                      tables).ravel()
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
            Tb = Tb + _tip_bracket(X[0], X[1], cols, w, e, D_stop, Dmin, omega0, Gamma0)
    else:
        nx = 8*d*d # complex entries of G and A_on
        y0 = np.concatenate([_bare_vertex(S, SIGN*Jrho_s, 0.).ravel(), A_bare.ravel(),
                             np.broadcast_to(T0, (n, 2, d, 2, m)).ravel()]).view(float)

        def couplings(y):
            return y.view(complex)[:nx].reshape(2, 4, d, d)

        def rhs(s, y):
            D = np.exp(s)
            X = couplings(y)
            Kf = _kernel(dE, D, omega0, Gamma0)
            Gb, GKb, Ab = _block(X[0]), _block(X[0]*Kf), _block(X[1])
            out = np.empty(len(y)//2, dtype=complex)
            out[:nx//2] = _unblock(_bracket(Gb, GKb)).ravel()
            out[nx//2:nx] = _unblock(_bracket(GKb, Ab)).ravel()
            out[nx:] = _tip_bracket(Gb, Ab, cols, w, e, D, -1., omega0, Gamma0,
                                    tables).ravel()
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
                                   Dmin, omega0, Gamma0)
    return _cut_columns(Tb), {"D_stop": D_stop,
                              "max_coupling": float(max_coupling/(smax if smax > 0. else 1.))}


def tip_vertex(ks, omegas, Us, Jrho_s, omega0=20e-3, Gamma0=5e-6, rtol=1e-10,
               atol=1e-14, pole_threshold=1.0, cols=None):
    """The tip vertex A(omega) at D -> 0 for each incoming energy
    omegas[k], with bare potential scattering Us[k]; shape (n,4,d,d), or
    (n,4,d,m) for the columns (initial states) cols only, which is what
    dI/dV needs and costs d^2 per bias and column instead of d^3. Also
    returns the flow's info: D_stop and max_coupling (see onshell_flow).

    The potential scattering of the tip does not flow: its direct and
    exchange terms cancel exactly at every cutoff. So the flow is solved
    once at U=0, and Us[k] is added to the mu=0 vertex at the end."""
    omegas = np.asarray(omegas, dtype=float)
    Us = np.broadcast_to(np.asarray(Us, dtype=float), omegas.shape)
    d = len(ks.e)
    cols = np.arange(d) if cols is None else np.asarray(cols, dtype=int)
    # a symmetric bias grid asks for each energy twice
    w, inverse = np.unique(omegas, return_inverse=True)
    A, info = _joint_flow(ks, w, cols, Jrho_s, omega0, Gamma0, rtol, atol,
                          pole_threshold)
    A = A[inverse]
    A[:, 0, cols, np.arange(len(cols))] += Us[:, None]
    return A, info


def leading_log_dIdV(ks, eVs, Jrho_s, U=0.0, T0=1.0, omega0=20e-3,
                     Gamma0=5e-6, p_cut=P_CUT, rtol=1e-10, atol=1e-14,
                     pole_threshold=1.0, stop_margin=10.0, return_info=False):
    """dI/dV with the leading logs summed to all orders, in the units of
    conductance.py (2*pi*e^2*T0^2/hbar). At first order in Jrho_s it is
    second_order_dIdV + third_order_kondo_dIdV + third_order_potential_dIdV
    exactly; see docs/new_methods/kondo_leading_log/README.md.

    T=0 only (ks.T must be 0). Two checks of validity, per bias:
      - the running tip coupling, |Jrho_s| times the largest ratio of the
        flowed to the bare spin vertex over the transitions out of occupied
        states that are open at this bias, the ones that enter dI/dV
        there, with the bare vertex bounded below by BARE_FLOOR. Above
        WARN_COUPLING it warns, at pole_threshold the bias is NaN.
      - if the sample coupling reaches pole_threshold during the flow (the
        Kondo scale is above Gamma0), the flow stops at that cutoff D_stop,
        and every bias within stop_margin*D_stop of zero or of any energy
        difference e_x - e_y of the spectrum is NaN, with a warning.

    A site with no spin on any state of the impurity has no exchange
    vertex: the flow does nothing, dI/dV is the potential scattering
    2*pi*T0^2*4*U^2 of the second order, and there is no warning.

    return_info=True also returns a dict with D_stop, max_coupling (the
    largest running sample coupling, Jrho_s itself at the start) and
    coupling (the running tip coupling at each bias)."""
    if ks.T != 0.:
        raise NotImplementedError("leading_log_dIdV is T=0 only for now; "
                                  "see docs/new_methods/kondo_leading_log")
    eVs = np.asarray(eVs, dtype=float)
    occ = _occupied_states(ks, p_cut)
    p = ks.p[occ]
    e = np.asarray(ks.e, dtype=float)
    n = len(eVs)
    # t->s at (eV, U) and its particle-hole conjugate s->t at (-eV, -U)
    omegas = np.concatenate([eVs, -eVs])
    Us = np.concatenate([np.full(n, U), np.full(n, -U)])
    # A[k,lam,f,i], i over the occupied states only
    A, info = tip_vertex(ks, omegas, Us, Jrho_s, omega0=omega0, Gamma0=Gamma0,
                         rtol=rtol, atol=atol, pole_threshold=pole_threshold,
                         cols=occ)
    D_stop = info["D_stop"]
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
    # a zero bare vertex (a site with no spin) has nothing to compare with
    ratio = np.max(np.where(Th > 0., flow/np.where(bare > 0., bare, np.inf)[None], 0.),
                   axis=(1, 2))
    coupling = abs(Jrho_s)*np.maximum(ratio[:n], ratio[n:])
    info["coupling"] = coupling
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
                          "one at D = %.3g (the Kondo scale is above Gamma0); "
                          "%d bias points within %g*D of an energy difference "
                          "of the spectrum are not leading-log results and are "
                          "returned as NaN" % (D_stop, int(np.sum(near)), stop_margin))
        bad = bad | near
    dIdV = np.where(bad, np.nan, dIdV)
    if return_info:
        return dIdV, info
    return dIdV
