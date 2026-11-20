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

ed_krylov_spectrum (submode="KRYLOV") needs as eigenstates only the
populated ones, up to epop = kT ln(1/p_cut) above the lowest state. The
final states come from one block Krylov space started from V|m> for every
populated m, and H is diagonalized in the span of the populated states and
that space. This follows the pattern of kondospectrumtk/krylovkondo.py, with
one difference that is measured: the relaxation needs one common basis. A
two-level Bloch line per transition, as separate Krylov spaces would give,
was 46 to 61% of the peak off the full rate equation on 8 and 10 spins,
while the full rate equation restricted to the populated states and the
final states with |V_nm|^2 above 1e-2 of the largest was within 4e-7 to
2e-4 of it. The final states of the lines in the window must converge
(ed_krylov_spectrum); the other Ritz states are kept as relaxation
partners.
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


# A direction of a Krylov block is kept when its singular value, after the
# projection on the space so far, is above DROP_TOL of the largest column
# norm of the block (krylovkondo.DROP_TOL)
DROP_TOL = 1e-9
# The Ritz vectors are formed and tested in chunks of this many columns
CHUNK = 64
# submode="KRYLOV": a line out of a populated state m into n counts as a line
# of the window when |V_nm|^2 p_m is above LINE_TOL of the strongest
LINE_TOL = 1e-6


class ESRSpectrum():
    """Energies (relative to the lowest state), V, D and K in the basis of
    the kept eigenstates. sigma bounds the error of each energy: the square
    root of <H^2> - <H>^2 on the MPS route, the Ritz residual on the Krylov
    route (None on ED, where the states are exact to rounding)."""
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

    def element(op):
        """<a|op|b> as a function of the two states. On the session
        backends the MPO is built once (StaticOperator), not once per
        matrix element as chain.aMb does for a MultiOperator: 1.7 times
        faster per element on 10 and 12 sites at maxm=64 on v3."""
        if (getattr(chain, "_session", None) is not None
                and chain.itensor_version in (2, 3, "python")):
            from ..multioperatortk.staticoperator import StaticOperator
            S = StaticOperator(op, chain)
            return S.aMb
        return lambda wa, wb: chain.aMb(wa, op, wb)

    def hermitian(op):
        amb = element(op)
        M = np.zeros((n, n), dtype=complex)
        for a in range(n):
            for b in range(a, n):
                M[a, b] = amb(ws[a], ws[b])
                M[b, a] = np.conjugate(M[a, b])
        return M

    def full(op):
        amb = element(op)
        M = np.zeros((n, n), dtype=complex)
        for a in range(n):
            for b in range(n):
                M[a, b] = amb(ws[a], ws[b])
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


def _orth(B, Q):
    """B projected twice on the complement of the columns of Q, and an
    orthonormal basis of what is left above DROP_TOL"""
    scale = float(np.max(np.linalg.norm(B, axis=0))) if B.size else 0.
    if scale == 0.:
        return B[:, :0]
    for _ in range(2):
        B = B - Q @ (Q.conj().T @ B)
    U, sv, _ = np.linalg.svd(B, full_matrices=False)
    return U[:, sv > DROP_TOL*scale]


def ed_krylov_spectrum(chain, V, D, emax, epop, window, kT, nkry, ritz_tol,
                       site_weights=None):
    """ESRSpectrum from the populated eigenstates (up to epop above the
    lowest) and the Ritz states of one block Krylov space started from V|m>
    for every populated m, up to nkry blocks, kept up to emax above the
    lowest state. window = (elo, ehi) is the energy range of the lines that
    are read (the frequency window and its margin, times h).

    The final states of the lines in the window, those n with
    |V_nm|^2 p_m above LINE_TOL of the strongest for a populated m and
    |E_n - E_m| in the window, must be converged: their residual
    ||(H - theta)|u>|| below ritz_tol, so that an eigenvalue lies within
    it, or the call raises and asks for more nkry. The other Ritz states
    are kept as they are, as relaxation partners, whose energies and spin
    matrix elements are then approximate; this needs about half the Krylov
    space of asking every state below emax to converge (measured on 10 and
    12 spins), and their sigma is set to 0."""
    w = _weights(chain, site_weights)
    edobj = chain.get_ED_obj()
    Hs = sp.csr_matrix(edobj.get_hamiltonian())
    if edkondo._imag_is_zero(Hs):
        Hs = _part(Hs, np.real)
    dim = Hs.shape[0]
    if dim <= algebra.maxsize:
        emu, vs = np.linalg.eigh(Hs.toarray())
    else:
        emu, vs = _low_states(Hs, epop)
    pop = emu - emu[0] <= epop
    P = np.ascontiguousarray(vs[:, pop])
    Vs = sp.csr_matrix(edobj.get_operator(V))
    if not edkondo._imag_is_zero(Vs) or np.iscomplexobj(P):
        Vs = sp.csr_matrix(Vs, dtype=complex)
        P = P.astype(complex)
        Hs = sp.csr_matrix(Hs, dtype=complex)
    else:
        Vs = _part(Vs, np.real)
    # the block Krylov space of V|m>, orthogonal to the populated states
    Q = P
    B = _orth(Vs @ P, Q)
    for _ in range(int(nkry)):
        if B.shape[1] == 0: break # the space is exhausted
        Q = np.hstack([Q, B])
        B = _orth(Hs @ B, Q)
    HQ = Hs @ Q
    Hk = Q.conj().T @ HQ
    theta, y = np.linalg.eigh(0.5*(Hk + Hk.conj().T))
    cand = np.flatnonzero(theta - theta[0] <= emax)
    res = np.empty(len(cand))
    for c0 in range(0, len(cand), CHUNK):
        yc = y[:, cand[c0:c0 + CHUNK]]
        R = HQ @ yc - (Q @ yc)*theta[cand[c0:c0 + CHUNK]][None, :]
        res[c0:c0 + CHUNK] = np.linalg.norm(R, axis=0)
    del HQ
    U = Q @ y[:, cand]
    del Q
    e = theta[cand] - theta[cand[0]]
    Vm = _transform(Vs, U)
    # the lines in the window out of the populated states
    popm = e <= epop
    p = np.zeros(len(e))
    p[popm] = np.exp(-e[popm]/kT) if kT > 0. else (e[popm] == 0.)
    strength = np.abs(Vm)**2*p[None, :] # zero unless m is populated
    de = np.abs(e[:, None] - e[None, :])
    inwin = (de + res[:, None] >= window[0]) & (de - res[:, None] <= window[1])
    smax = np.max(strength*inwin) if strength.size else 0.
    line = inwin & (strength > LINE_TOL*smax) & (smax > 0.)
    need = np.any(line, axis=1)
    bad = need & (res > ritz_tol)
    if np.any(bad):
        raise ValueError(
            "get_esr_spectrum(submode=\"KRYLOV\"): %d of the %d final states "
            "of the lines in the window are not converged (largest residual "
            "%.3g, ritz_tol %.3g); raise nkry (%d)"
            % (np.sum(bad), np.sum(need), np.max(res[bad]), ritz_tol, nkry))
    Dm = _transform(edobj.get_operator(D), U)
    K = np.zeros((len(e), len(e)))
    Sz, Sp = site_operators(chain)
    for k in range(len(w)):
        if w[k] == 0.: continue
        z = _transform(edobj.get_operator(Sz[k]), U)
        pk = _transform(edobj.get_operator(Sp[k]), U)
        p2 = np.abs(pk)**2
        K += w[k]*(np.abs(z)**2 + 0.5*(p2 + p2.T))
    sigma = np.where(need | popm, res, 0.)
    return ESRSpectrum(e, Vm, Dm, K, sigma=sigma, route="KRYLOV",
                       reached=theta[-1] - theta[0])
