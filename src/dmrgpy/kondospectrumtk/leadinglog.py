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
eigenbasis, the convention of KondoSpectrum. The sample vertex G is
needed only with its incoming electron at the Fermi level, so it is one
such set of matrices that flows alone; the tip vertex at bias omega is a
second set that flows linearly, driven by G.
"""
import warnings

import numpy as np
from scipy.integrate import solve_ivp

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


def _kernel(x, D, omega0, Gamma0):
    """K(x,D) = 1/(D+r_x) - 1/(D+omega0+|x|), r_x = sqrt(x^2+Gamma0^2):
    its integral over D from 0 to infinity is exactly F0(x) of
    stepfunctions.py, so the flow is the third order at first order."""
    ax = np.abs(x)
    return 1./(D + np.sqrt(ax*ax + Gamma0*Gamma0)) - 1./(D + omega0 + ax)


def _rate(G, A, Kmat):
    """The bracket of the flow, R[lam] = sum_c sum_{mu nu}
    C[mu,nu,lam] G[mu]_ac (K*A[nu])_cb - C[nu,mu,lam] (A[mu]*K)_ac G[nu]_cb.
    G: (4,d,d); A: (...,4,d,d); Kmat: (...,d,d), Kmat[x,y] = K(omega-(e_x-e_y)).
    Both cuts are read off the same matrix: the direct one at
    omega-(e_c-e_b) is Kmat[c,b], the exchange one at omega-(e_a-e_c) is
    Kmat[a,c]."""
    KA = A*Kmat[..., None, :, :]
    part = np.matmul(G[:, None], KA[..., None, :, :, :]) # (...,mu,nu,d,d)
    hole = np.matmul(KA[..., :, None, :, :], G[None]) # (...,mu,nu,d,d)
    return (np.einsum('mnl,...mnab->...lab', _C, part)
            - np.einsum('mnl,...mnab->...lab', _CR, hole))


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


def sample_flow(ks, Jrho_s, omega0=20e-3, Gamma0=5e-6, rtol=1e-10,
                atol=1e-14, pole_threshold=1.0, Dmax=None, Dmin=None):
    """Flow of the sample vertex G (incoming electron at the Fermi level)
    from D=Dmax down to Dmin. Returns (sol, D_stop): sol.sol(ln D) gives
    G flattened, and D_stop is the cutoff where the coupling reached
    order one (None if it never did). Below D_stop G is held at its value
    there, and results that depend on that range are not leading-log
    results any more.

    pole_threshold: the flow stops once max |G| over the exchange
    components exceeds pole_threshold times max|S|/2, i.e. once the
    running dimensionless coupling reaches about pole_threshold."""
    S = _spin_ops(ks)
    d = S.shape[-1]
    e = np.asarray(ks.e, dtype=float)
    dE = e[:, None] - e[None, :]
    if Dmax is None or Dmin is None:
        Dmax0, Dmin0 = _energy_range(ks, omega0, np.zeros(1), Gamma0)
        Dmax = Dmax0 if Dmax is None else Dmax
        Dmin = Dmin0 if Dmin is None else Dmin
    G0 = _bare_vertex(S, SIGN*Jrho_s, 0.)
    smax = np.max(np.abs(S))/2.
    ref = pole_threshold*smax

    def rhs(s, y):
        D = np.exp(s)
        G = y.view(complex).reshape(4, d, d)
        K = _kernel(dE, D, omega0, Gamma0)
        return (-D*_rate(G, G, K)).reshape(-1).view(float)

    def pole(s, y):
        G = y.view(complex).reshape(4, d, d)
        return ref - np.max(np.abs(G[1:]))
    pole.terminal = True
    pole.direction = -1

    sol = solve_ivp(rhs, (np.log(Dmax), np.log(Dmin)), G0.reshape(-1).view(float),
                    method="DOP853", rtol=rtol, atol=atol, dense_output=True,
                    events=pole)
    if not sol.success:
        raise RuntimeError("leading-log sample flow failed: "+str(sol.message))
    D_stop = None
    if len(sol.t_events[0]) > 0:
        D_stop = float(np.exp(sol.t_events[0][0]))
    # the largest running coupling reached, in units of the bare spin
    # vertex S/2: the leading-log result is controlled only while it is small
    sol.max_coupling = float(np.max([np.max(np.abs(y.view(complex).reshape(4, d, d)[1:]))
                                     for y in sol.y.T])/smax)
    return sol, D_stop


def tip_vertex(ks, omegas, Us, Jrho_s, omega0=20e-3, Gamma0=5e-6, rtol=1e-10,
               atol=1e-14, pole_threshold=1.0, chunk=None):
    """The flowed tip vertex A(omega) at D -> 0 for each incoming energy
    omegas[k], with bare potential scattering Us[k]; shape (n,4,d,d).
    Also returns the sample flow's info: D_stop (see sample_flow) and
    max_coupling, the largest running coupling it reached."""
    omegas = np.asarray(omegas, dtype=float)
    Us = np.broadcast_to(np.asarray(Us, dtype=float), omegas.shape)
    S = _spin_ops(ks)
    d = S.shape[-1]
    e = np.asarray(ks.e, dtype=float)
    dE = e[:, None] - e[None, :]
    Dmax, Dmin = _energy_range(ks, omega0, omegas, Gamma0)
    gsol, D_stop = sample_flow(ks, Jrho_s, omega0=omega0, Gamma0=Gamma0,
                               rtol=rtol, atol=atol,
                               pole_threshold=pole_threshold, Dmax=Dmax, Dmin=Dmin)
    s_stop = gsol.t[-1] # the flow ran from log(Dmax) down to here
    G_frozen = gsol.y[:, -1].view(complex).reshape(4, d, d)

    def G_at(s):
        if s <= s_stop: return G_frozen
        return gsol.sol(s).view(complex).reshape(4, d, d)

    if chunk is None: chunk = max(1, int(2e6)//(16*d*d))
    out = np.empty((len(omegas), 4, d, d), dtype=complex)
    for k0 in range(0, len(omegas), chunk):
        w = omegas[k0:k0+chunk]
        n = len(w)
        A0 = np.stack([_bare_vertex(S, 1.0, u) for u in Us[k0:k0+chunk]])
        x = w[:, None, None] - dE[None]

        def rhs(s, y):
            D = np.exp(s)
            A = y.view(complex).reshape(n, 4, d, d)
            K = _kernel(x, D, omega0, Gamma0)
            return (-D*_rate(G_at(s), A, K)).reshape(-1).view(float)

        sol = solve_ivp(rhs, (np.log(Dmax), np.log(Dmin)), A0.reshape(-1).view(float),
                        method="DOP853", rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError("leading-log tip flow failed: "+str(sol.message))
        out[k0:k0+n] = sol.y[:, -1].view(complex).reshape(n, 4, d, d)
    return out, {"D_stop": D_stop, "max_coupling": gsol.max_coupling}


def leading_log_dIdV(ks, eVs, Jrho_s, U=0.0, T0=1.0, omega0=20e-3,
                     Gamma0=5e-6, p_cut=P_CUT, rtol=1e-10, atol=1e-14,
                     pole_threshold=1.0, stop_margin=10.0, return_info=False):
    """dI/dV with the leading logs summed to all orders, in the units of
    conductance.py (2*pi*e^2*T0^2/hbar). At first order in Jrho_s it is
    second_order_dIdV + third_order_kondo_dIdV + third_order_potential_dIdV
    exactly; see docs/new_methods/kondo_leading_log/README.md.

    T=0 only (ks.T must be 0). If the sample coupling reaches order one
    during the flow (the Kondo scale is above Gamma0), the flow is stopped
    at that cutoff D_stop, and every bias within stop_margin*D_stop of a
    threshold e_f - e_i out of an occupied state, or of zero bias, is
    returned as NaN with a warning: the method is not valid there.

    return_info=True also returns a dict with D_stop and max_coupling,
    the largest running sample coupling of the flow in units of its bare
    spin vertex (Jrho_s itself at the start): the result is a controlled
    leading-log one only while it stays well below 1."""
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
    A, info = tip_vertex(ks, omegas, Us, Jrho_s, omega0=omega0,
                           Gamma0=Gamma0, rtol=rtol, atol=atol,
                           pole_threshold=pole_threshold)
    D_stop = info["D_stop"]
    # W[k,f,i] = 4 sum_lam |A^lam_fi(omega_k)|^2, i over occupied states
    W = 4.*np.sum(np.abs(A[:, :, :, occ])**2, axis=1)
    eps_if = e[:, None] - e[occ][None, :] # [f,i] = e_f - e_i
    Th = _theta_raw(ks, omegas[:, None, None] - eps_if[None])
    g = np.einsum('i,kfi,kfi->k', p, W, Th)
    dIdV = 2*np.pi*T0**2*(g[:n] + g[n:])
    if D_stop is not None:
        thr = np.concatenate([[0.], np.abs(eps_if).ravel()])
        bad = np.min(np.abs(np.abs(eVs)[:, None] - thr[None]), axis=1) < stop_margin*D_stop
        if np.any(bad):
            warnings.warn("leading_log_dIdV: the sample coupling reached order "
                          "one at D = %.3g (the Kondo scale is above Gamma0); "
                          "%d bias points within %g*D of a threshold are not "
                          "leading-log results and are returned as NaN"
                          % (D_stop, int(np.sum(bad)), stop_margin))
            dIdV = np.where(bad, np.nan, dIdV)
    if return_info:
        return dIdV, info
    return dIdV
