"""The eigenstates an ESR-STM spectrum is read from.

Every builder returns an ESRSpectrum: the energies e of the states that can
take part in a transition inside the frequency window (relative to the
lowest one, in the units of H), the matrices of the drive V and of the
detection D in that basis, and the Korringa kernel

    K[m, n] = sum_k w_k sum_a |<m|S^a_k|n>|^2,

the matrix element of the collapse operators between eigenstates of eq. 42
of arXiv:2606.27830 (from Broekhoven et al., npj Quantum Inf. 10, 92). It
is built per site from S^z_k and S^+_k, since
|S^x_mn|^2 + |S^y_mn|^2 = (|S^+_mn|^2 + |S^+_nm|^2)/2, and for real
eigenvectors both are real matrices.

ed_spectrum diagonalizes the whole Hamiltonian up to algebra.maxsize states
and takes the lowest eigenpairs of the sparse H above it. mps_spectrum takes
the lowest eigenstates from DMRG (get_excited_states, which orthogonalizes
them and diagonalizes H in their span) and reads every matrix element as
<a|O|b> between them. Both keep the states up to emax above the lowest one,
and the MPS route raises when its states do not reach emax, since a state
missing below emax would be a missing line or a missing population.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from ..algebra import algebra
from ..kondospectrumtk import edkondo

# Above this dimension the ED route takes the lowest eigenpairs of the
# sparse H with eigsh, starting from K0 states and raising k by 1.6 until
# the last one is above emax. When eigsh would need more than a third of
# the states, a dense eigh is made up to DENSE_MAX states, and the call
# raises above it (the same rule as kondospectrumtk/krylovkondo.low_states).
K0 = 32
DENSE_MAX = 16384


class ESRSpectrum():
    """Energies (relative to the lowest state), V, D and K in the basis of
    the kept eigenstates. sigma is sqrt(<H^2> - <H>^2) of each state on the
    MPS route (None on ED, where the states are exact to rounding)."""
    def __init__(self, e, V, D, K, sigma=None, route="ED", reached=None):
        self.e = np.asarray(e, dtype=float)
        self.dim = len(self.e)
        self.V = np.asarray(V)
        self.D = np.asarray(D)
        self.K = np.asarray(K, dtype=float)
        self.sigma = None if sigma is None else np.asarray(sigma, dtype=float)
        self.route = route
        # the highest energy the builder solved for, relative to the lowest
        # state (the MPS route keeps fewer states than it solves)
        self.reached = reached


def site_operators(chain):
    """(Sz, Sp): the lists of S^z_k and S^+_k = S^x_k + i S^y_k of the chain"""
    Sz = list(chain.Sz)
    Sp = [chain.Sx[k] + 1j*chain.Sy[k] for k in range(len(chain.Sx))]
    return Sz, Sp


def _weights(chain, site_weights):
    n = len(chain.Sx)
    if site_weights is None:
        return np.ones(n)
    w = np.asarray(site_weights, dtype=float)
    if w.shape != (n,):
        raise ValueError("site_weights must have one entry per site (%d), "
                         "got shape %s" % (n, w.shape))
    if np.any(w < 0.) or not np.any(w > 0.):
        raise ValueError("site_weights must be non-negative with at least "
                         "one positive entry")
    return w


def _transform(M, U):
    """U^dagger M U for a sparse or dense M, in real arithmetic when U is
    real and M is real or purely imaginary (exact zeros, as edkondo)."""
    if not sp.issparse(M):
        M = np.asarray(M)
    if not np.iscomplexobj(U):
        if edkondo._imag_is_zero(M):
            Mr = M.real if not sp.issparse(M) else _part(M, np.real)
            return U.T @ (Mr @ U)
        if edkondo._real_part_is_zero(M):
            Mi = M.imag if not sp.issparse(M) else _part(M, np.imag)
            return 1j*(U.T @ (Mi @ U))
    if sp.issparse(M):
        M = sp.csr_matrix(M, dtype=complex)
    return np.conjugate(U.T) @ (M @ U)


def _part(M, part):
    c = M.tocoo() # a copy: M is not touched (see edkondo._values)
    return sp.csr_matrix((np.array(part(c.data), dtype=float),
                          (c.row, c.col)), shape=c.shape)


def _low_states(H, emax):
    """The lowest eigenpairs of the sparse H, sorted, enough of them that
    the last one lies more than emax above the first one."""
    dim = H.shape[0]
    k = min(K0, dim - 2)
    while k <= dim//3:
        w, v = sla.eigsh(H, k=k, which="SA", tol=0.)
        o = np.argsort(w)
        w, v = w[o], v[:, o]
        if w[-1] - w[0] > emax:
            return w, v
        k = int(1.6*k) + 1
    if dim > DENSE_MAX:
        raise ValueError(
            "get_esr_spectrum(mode=\"ED\"): more than a third of the %d "
            "states lie within %.3g of the ground state, the energy the "
            "window and the temperature need, and a dense diagonalization "
            "of that size is not done; raise p_cut, lower T or narrow the "
            "window" % (dim, emax))
    return np.linalg.eigh(H.toarray())


def ed_spectrum(chain, V, D, emax, site_weights=None):
    """ESRSpectrum from ED, keeping the states up to emax above the lowest"""
    w = _weights(chain, site_weights)
    edobj = chain.get_ED_obj()
    H = edobj.get_hamiltonian()
    dim = H.shape[0]
    if dim <= algebra.maxsize:
        emu, vs = edkondo.eigenpairs(edobj) # real for a real H
        emu = np.asarray(emu, dtype=float)
        o = np.argsort(emu)
        emu, vs = emu[o], np.asarray(vs)[:, o]
    else:
        Hs = sp.csr_matrix(H)
        if edkondo._imag_is_zero(Hs):
            Hs = _part(Hs, np.real)
        emu, vs = _low_states(Hs, emax)
    e = emu - emu[0]
    keep = e <= emax
    U = np.ascontiguousarray(vs[:, keep])
    Vm = _transform(edobj.get_operator(V), U)
    Dm = _transform(edobj.get_operator(D), U)
    K = np.zeros((U.shape[1], U.shape[1]))
    Sz, Sp = site_operators(chain)
    for k in range(len(w)):
        if w[k] == 0.: continue
        z = _transform(edobj.get_operator(Sz[k]), U)
        p = _transform(edobj.get_operator(Sp[k]), U)
        p2 = np.abs(p)**2
        K += w[k]*(np.abs(z)**2 + 0.5*(p2 + p2.T))
    return ESRSpectrum(e[keep], Vm, Dm, K, route="ED", reached=e[-1])


def mps_spectrum(chain, V, D, emax, nstates, site_weights=None):
    """ESRSpectrum from the nstates lowest DMRG eigenstates. Raises when the
    highest of them is not above emax, so that no state below emax can be
    missing from the window, as long as the excited-state solver finds the
    states in order (the penalty method can skip a state of a sector its
    start does not reach; see get_excited_states)."""
    from .. import vev
    w = _weights(chain, site_weights)
    es, ws = chain.get_excited_states(n=int(nstates))
    es = np.real(np.asarray(es))
    if len(ws) < nstates:
        raise ValueError("get_esr_spectrum(mode=\"DMRG\"): asked for %d "
                         "states, the excited-state solver returned %d"
                         % (nstates, len(ws)))
    o = np.argsort(es)
    es = es[o]
    ws = [ws[i] for i in o]
    e = es - es[0]
    if not e[-1] > emax:
        raise ValueError(
            "get_esr_spectrum(mode=\"DMRG\"): the %d lowest states reach "
            "%.4g above the ground state, but the window and the "
            "temperature need every state up to %.4g; raise nstates (or "
            "raise p_cut, lower T, narrow the window)"
            % (nstates, e[-1], emax))
    keep = np.flatnonzero(e <= emax)
    ws = [ws[i] for i in keep]
    norms = np.array([np.sqrt(abs(chain.overlap(x, x))) for x in ws])
    ws = [x*(1./nx) for x, nx in zip(ws, norms)]
    n = len(ws)
    H = chain.hamiltonian
    sigma = np.array([np.sqrt(max(np.real(vev.energy_variance(chain, H, wf=x)), 0.))
                      for x in ws])

    def hermitian(op):
        M = np.zeros((n, n), dtype=complex)
        for a in range(n):
            for b in range(a, n):
                M[a, b] = chain.aMb(ws[a], op, ws[b])
                M[b, a] = np.conjugate(M[a, b])
        return M

    def full(op):
        M = np.zeros((n, n), dtype=complex)
        for a in range(n):
            for b in range(n):
                M[a, b] = chain.aMb(ws[a], op, ws[b])
        return M

    Vm = hermitian(V)
    Dm = hermitian(D)
    K = np.zeros((n, n))
    Sz, Sp = site_operators(chain)
    for k in range(len(w)):
        if w[k] == 0.: continue
        z2 = np.abs(hermitian(Sz[k]))**2
        p2 = np.abs(full(Sp[k]))**2
        K += w[k]*(z2 + 0.5*(p2 + p2.T))
    return ESRSpectrum(e[keep], Vm, Dm, K, sigma=sigma, route="DMRG",
                       reached=e[-1])
