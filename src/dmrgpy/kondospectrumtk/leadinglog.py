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

# a running coupling above this is reported with a warning: the
# leading-log result is then only qualitative
WARN_COUPLING = 0.5


def _kernel(x, D, omega0, Gamma0):
    """K(x,D) = 1/(D+r_x) - 1/(D+omega0+|x|), r_x = sqrt(x^2+Gamma0^2):
    its integral over D from 0 to infinity is exactly F0(x) of
    stepfunctions.py, so the flow is the third order at first order."""
    ax = np.abs(x)
    return 1./(D + np.sqrt(ax*ax + Gamma0*Gamma0)) - 1./(D + omega0 + ax)


def _combine(part, hole):
    return (np.einsum('mnl,...mnab->...lab', _C, part)
            - np.einsum('mnl,...mnab->...lab', _CR, hole))


def _rate(G, A, Kmat):
    """The bracket of the flow of a vertex with incoming energy omega,
    R[lam] = sum_c sum_{mu nu} C[mu,nu,lam] G[mu]_ac A[nu]_cb K[c,b]
                             - C[nu,mu,lam] A[mu]_ac G[nu]_cb K[a,c],
    with Kmat[x,y] = K(omega-(e_x-e_y)): the direct cut at
    omega-(e_c-e_b), the exchange cut at omega-(e_a-e_c).
    G: (4,d,d); A: (...,4,d,d), the on-shell vertex on the right-hand
    side; Kmat: (...,d,d)."""
    KA = A*Kmat[..., None, :, :]
    part = np.matmul(G[:, None], KA[..., None, :, :, :])
    hole = np.matmul(KA[..., :, None, :, :], G[None])
    return _combine(part, hole)


def _onshell_rate(G, A, Kfix):
    """The same bracket for the element (a,b) taken at its own on-shell
    energy omega = e_a - e_b: the direct cut becomes e_a - e_c and the
    exchange cut e_c - e_b, Kfix[x,y] = K(e_x - e_y)."""
    GK = G*Kfix
    part = np.matmul(GK[:, None], A[..., None, :, :, :])
    hole = np.matmul(A[..., :, None, :, :], GK[None])
    return _combine(part, hole)


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
    """Flow of the on-shell couplings from D=Dmax down to Dmin: the sample
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


def tip_vertex(ks, omegas, Us, Jrho_s, omega0=20e-3, Gamma0=5e-6, rtol=1e-10,
               atol=1e-14, pole_threshold=1.0, chunk=None):
    """The tip vertex A(omega) at D -> 0 for each incoming energy
    omegas[k], with bare potential scattering Us[k], which must be +U or
    -U for one U; shape (n,4,d,d). Also returns the flow's info: D_stop
    and max_coupling (see onshell_flow)."""
    omegas = np.asarray(omegas, dtype=float)
    Us = np.broadcast_to(np.asarray(Us, dtype=float), omegas.shape)
    U = float(Us[0]) if len(Us) else 0.
    if not np.all(np.isclose(np.abs(Us), abs(U))):
        raise ValueError("tip_vertex takes one U and its opposite only")
    sel = np.where(np.isclose(Us, U), 1, 2) # which A_on drives each bias
    S = _spin_ops(ks)
    d = S.shape[-1]
    e = np.asarray(ks.e, dtype=float)
    dE = e[:, None] - e[None, :]
    Dmax, Dmin = _energy_range(ks, omega0, omegas, Gamma0)
    fsol, D_stop = onshell_flow(ks, Jrho_s, U=U, omega0=omega0, Gamma0=Gamma0,
                                rtol=rtol, atol=atol, pole_threshold=pole_threshold,
                                Dmax=Dmax, Dmin=Dmin)
    s_stop = fsol.t[-1] # the flow ran from log(Dmax) down to here
    X_frozen = fsol.y[:, -1].view(complex).reshape(3, 4, d, d)

    def X_at(s):
        if s <= s_stop: return X_frozen
        return fsol.sol(s).view(complex).reshape(3, 4, d, d)

    if chunk is None: chunk = max(1, int(2e6)//(16*d*d))
    out = np.empty((len(omegas), 4, d, d), dtype=complex)
    for k0 in range(0, len(omegas), chunk):
        w = omegas[k0:k0+chunk]
        n = len(w)
        A0 = np.stack([_bare_vertex(S, 1.0, u) for u in Us[k0:k0+chunk]])
        x = w[:, None, None] - dE[None]
        idx = sel[k0:k0+chunk]

        def rhs(s, y):
            D = np.exp(s)
            X = X_at(s)
            K = _kernel(x, D, omega0, Gamma0)
            return (-D*_rate(X[0], X[idx], K)).reshape(-1).view(float)

        sol = solve_ivp(rhs, (np.log(Dmax), np.log(Dmin)), A0.reshape(-1).view(float),
                        method="DOP853", rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError("leading-log tip vertex failed: "+str(sol.message))
        out[k0:k0+n] = sol.y[:, -1].view(complex).reshape(n, 4, d, d)
    return out, {"D_stop": D_stop, "max_coupling": fsol.max_coupling}


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
        states. Above WARN_COUPLING it warns, at pole_threshold the bias is
        NaN.
      - if the sample coupling reaches pole_threshold during the flow (the
        Kondo scale is above Gamma0), the flow stops at that cutoff D_stop,
        and every bias within stop_margin*D_stop of zero or of any energy
        difference e_x - e_y of the spectrum is NaN, with a warning.

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
    A, info = tip_vertex(ks, omegas, Us, Jrho_s, omega0=omega0, Gamma0=Gamma0,
                         rtol=rtol, atol=atol, pole_threshold=pole_threshold)
    D_stop = info["D_stop"]
    # W[k,f,i] = 4 sum_lam |A^lam_fi(omega_k)|^2, i over occupied states
    W = 4.*np.sum(np.abs(A[:, :, :, occ])**2, axis=1)
    eps_if = e[:, None] - e[occ][None, :] # [f,i] = e_f - e_i
    Th = _theta_raw(ks, omegas[:, None, None] - eps_if[None])
    g = np.einsum('i,kfi,kfi->k', p, W, Th)
    dIdV = 2*np.pi*T0**2*(g[:n] + g[n:])
    # the running tip coupling at each bias
    S = _spin_ops(ks)
    bare = np.sqrt(np.sum(np.abs(S[:, :, occ]/2.)**2, axis=0)) # [f,i]
    mask = bare > 1e-8*np.max(bare)
    flow = np.sqrt(np.sum(np.abs(A[:, 1:, :, :][:, :, :, occ])**2, axis=1)) # [k,f,i]
    ratio = np.max(np.where(mask[None], flow/np.where(mask, bare, 1.)[None], 0.), axis=(1, 2))
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
