"""Relaxation between eigenstates for the ESR-STM rate equation.

The rate from n to m is the Korringa rate of eq. 42 of arXiv:2606.27830,

    Gamma_{n->m} = g K[m, n] phi(E_n - E_m),  phi(x) = x/(1 - exp(-x/kT)),

so phi(x) = x for x >> kT (emission), x exp(-|x|/kT) for absorption, and
kT at x = 0. It obeys detailed balance at T, so the stationary populations
without drive are the Boltzmann ones. The scale g is fixed by T1, the
relaxation time a free spin 1/2 on the tip site would have at the energy
eps0 = h f_ref: for that spin K = 1/2 between its two states, and

    1/T1 = Gamma_up + Gamma_down = (g/2) eps0 coth(eps0/2kT).

Every coherence rho_nm decays at gamma_nm = (Gamma_n + Gamma_m)/2 +
1/Tphi, Gamma_n the total rate out of n and Tphi the pure dephasing time,
which the free spin fixes too: 1/T2 = 1/(2 T1) + 1/Tphi.

Rates are in 1/ns when T1 and T2 are in ns.
"""
import numpy as np


def phi(x, kT):
    """x/(1 - exp(-x/kT)), the thermal factor of a Korringa rate; max(x, 0)
    at kT = 0."""
    x = np.asarray(x, dtype=float)
    if kT == 0.:
        return np.maximum(x, 0.)
    out = np.empty_like(x)
    u = x/kT
    small = np.abs(u) < 1e-8
    out[small] = kT*(1. + 0.5*u[small])
    big = ~small
    out[big] = x[big]/(-np.expm1(-u[big]))
    return out


def korringa_scale(T1, eps0, kT):
    """g such that a free spin 1/2 with splitting eps0 relaxes in T1"""
    if not eps0 > 0.:
        raise ValueError("the reference energy h*f_ref must be above zero, "
                         "got %r" % (eps0,))
    c = 1. if kT == 0. else 1./np.tanh(eps0/(2.*kT))
    return 2./(T1*eps0*c)


def dephasing_rate(T1, T2):
    """1/Tphi = 1/T2 - 1/(2 T1), the pure dephasing of the free spin"""
    if not (T1 > 0. and T2 > 0.):
        raise ValueError("T1 and T2 must be positive, got %r and %r"
                         % (T1, T2))
    r = 1./T2 - 0.5/T1
    if r < -1e-12/T2:
        raise ValueError("T2 cannot exceed 2*T1 (T1=%g, T2=%g)" % (T1, T2))
    return max(r, 0.)


def rate_matrix(e, K, g, kT):
    """Gamma[m, n], the rate from n to m (zero diagonal)"""
    x = e[:, None] - e[None, :] # x[n, m] = E_n - E_m
    G = g*K*phi(x, kT).T # Gamma[m, n] = g K[m, n] phi(E_n - E_m)
    np.fill_diagonal(G, 0.)
    return G


def generator(Gamma):
    """R with dp/dt = R p: R[m, n] = Gamma[m, n], R[n, n] = -sum_m Gamma[m, n]"""
    R = Gamma.copy()
    R[np.diag_indices_from(R)] = -Gamma.sum(axis=0)
    return R


def stationary(e, kT):
    """The Boltzmann populations of the kept states, one-hot on the lowest at
    kT = 0 (the caller checks that it is not degenerate)"""
    if kT == 0.:
        p = np.zeros(len(e))
        p[np.argmin(e)] = 1.
        return p
    p = np.exp(-(e - np.min(e))/kT)
    return p/np.sum(p)


# Two states are connected by the relaxation when their rate is above
# BLOCK_TOL of the largest rate (the spin operators do not move charge, so
# the charge sectors of a fermion chain are separate blocks).
BLOCK_TOL = 1e-14


def blocks(Gamma):
    """Labels of the blocks of states the relaxation connects"""
    from scipy.sparse.csgraph import connected_components
    A = Gamma + Gamma.T
    big = np.max(A) if A.size else 0.
    adj = A > BLOCK_TOL*big if big > 0. else np.zeros_like(A, dtype=bool)
    return connected_components(adj, directed=False)[1]


def resolvent(R, p0, labels, e):
    """G, block diagonal over the blocks the relaxation connects: for a
    vector x with sum zero on every block, y = G x solves R y = x with sum
    zero on every block. On each block G = (R_B - pi_B 1^T)^{-1}, pi_B the
    stationary populations of the block (p0 restricted and normalized, or
    its lowest state when p0 vanishes there). The populations of a block
    stay at their weight in p0 under the drive, as long as the drive does
    not connect two blocks (the caller checks). Raises when a block has
    more than one stationary state."""
    n = len(p0)
    G = np.zeros((n, n))
    for b in np.unique(labels):
        idx = np.flatnonzero(labels == b)
        pb = p0[idx]
        if np.sum(pb) > 0.:
            pi = pb/np.sum(pb)
        else:
            pi = np.zeros(len(idx))
            pi[np.argmin(e[idx])] = 1.
        RB = R[np.ix_(idx, idx)]
        M = RB - np.outer(pi, np.ones(len(idx)))
        try:
            GB = np.linalg.inv(M)
        except np.linalg.LinAlgError:
            raise ValueError(_DISCONNECTED)
        if not np.all(np.isfinite(GB)):
            raise ValueError(_DISCONNECTED)
        if len(idx) > 1:
            x = np.zeros(len(idx))
            x[0], x[-1] = 1., -1.
            if np.linalg.norm(RB @ (GB @ x) - x) > 1e-6*np.linalg.norm(x):
                raise ValueError(_DISCONNECTED)
        G[np.ix_(idx, idx)] = GB
    return G


_DISCONNECTED = (
    "get_esr_spectrum: the relaxation does not connect every kept state to "
    "the others (the rate matrix has more than one stationary state). This "
    "happens with site_weights that leave a part of the system without "
    "relaxation, or at T=0 with a state that has no lower state reached by "
    "any S_k; give a weight to more sites or a small T.")
