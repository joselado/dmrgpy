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
# submode="KRYLOV" on DMRG: directions of the common basis with an overlap
# eigenvalue below OVERLAP_CUT of the largest are dropped (two populated
# states reaching one final state give two nearly parallel Ritz vectors)
OVERLAP_CUT = 1e-10
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


def local_operator(chain, op):
    """(site, matrix) of an operator whose terms all act on one MPS site k:
    the d x d matrix (bra, ket) o with op = 1 x ... x o x ... x 1, read off
    its MPO as the partial trace over every other site, divided by their
    dimensions (one site at a time), which holds for any bond dimension of
    the MPO. None for an operator on several sites (a spinful fermion's S
    spans the two spin-orbital sites) or on a backend that does not export
    MPO arrays."""
    if chain.itensor_version not in (3, "python"):
        return None
    try:
        terms = op.to_terms()
    except Exception:
        return None
    sites = {f[1] for _, fs in terms for f in fs}
    if len(sites) != 1:
        return None
    k = sites.pop() - 1 # to_terms counts from 1
    s = chain._session
    M = s.build_operator(terms)
    if chain.itensor_version == "python":
        from ..pyitensor import task
        Ws = task.mpo_arrays(M)
    else:
        Ws = s.mpo_site_arrays(M)
    Ws = [np.asarray(W) for W in Ws]
    left = np.ones(1)
    for W in Ws[:k]:
        left = left @ (np.einsum("aiib->ab", W)/W.shape[1])
    right = np.ones(1)
    for W in Ws[:k:-1]:
        right = (np.einsum("aiib->ab", W)/W.shape[1]) @ right
    return k, np.einsum("a,aijb,b->ij", left, Ws[k], right)


def site_traces(Aa, Ab):
    """G[k][s, t] = <a| with the site-k legs open |b> (s the bra leg, t the
    ket leg) for every site k, from the site arrays of the two MPS: one
    left and one right sweep of the overlap, so that <a|O_k|b> =
    sum_st O[s, t] G[k][s, t] for every one-site operator at once."""
    L = len(Aa)
    left = [None]*(L + 1)
    left[0] = np.ones((1, 1))
    for l in range(L):
        X = np.tensordot(left[l], Ab[l], axes=(1, 0)) # (x, t, y')
        left[l + 1] = np.tensordot(Aa[l].conj(), X, axes=([0, 1], [0, 1]))
    right = np.ones((1, 1))
    G = [None]*L
    for l in range(L - 1, -1, -1):
        X = np.tensordot(Ab[l], right, axes=(2, 1)) # (y, t, x')
        Y = np.tensordot(left[l], X, axes=(1, 0)) # (x, t, x')
        G[l] = np.tensordot(Aa[l].conj(), Y, axes=([0, 2], [0, 2]))
        right = np.tensordot(Aa[l].conj(), X, axes=([1, 2], [1, 2]))
    return G


class _PairTraces():
    """site_traces for every pair a <= b of a list of states, and the
    matrices of one-site operators from them: one sweep per pair instead
    of one per pair and operator (the n sites of the Korringa kernel take
    3n/2 operators)."""
    def __init__(self, chain, ws):
        s = chain._session
        arrs = [[np.asarray(x) for x in s.mps_site_arrays(w.cpp_handle)]
                for w in ws]
        n = len(arrs)
        self.n = n
        self.G = {}
        for a in range(n):
            for b in range(a, n):
                self.G[a, b] = site_traces(arrs[a], arrs[b])

    def matrix(self, site, O):
        M = np.zeros((self.n, self.n), dtype=complex)
        for (a, b), G in self.G.items():
            M[a, b] = np.sum(O*G[site])
            if a != b:
                M[b, a] = np.sum(O*G[site].T.conj())
        return M


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
    return _states_spectrum(chain, ws, e[keep], V, D, site_weights,
                            route="DMRG", reached=e[-1])


def _states_spectrum(chain, ws, e, V, D, site_weights=None, route="DMRG",
                     reached=None):
    """ESRSpectrum of a list of MPS states with energies e (relative to the
    lowest): normalized here, sigma from vev.energy_variance, and V, D and K
    from their matrix elements."""
    from .. import vev
    w = _weights(chain, site_weights)
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

    # every one-site operator from one sweep per pair of states, on the
    # backends that export site arrays; aMb otherwise, and for a drive or a
    # detection that is not a one-site operator
    Sz, Sp = site_operators(chain)
    pt = None
    if chain.itensor_version in (3, "python"):
        pt = _PairTraces(chain, ws)

    def matrix(op, herm=True):
        loc = local_operator(chain, op) if pt is not None else None
        if loc is not None:
            return pt.matrix(*loc)
        return hermitian(op) if herm else full(op)

    Vm = matrix(V)
    Dm = matrix(D)
    K = np.zeros((n, n))
    for k in range(len(w)):
        if w[k] == 0.: continue
        z2 = np.abs(matrix(Sz[k]))**2
        p2 = np.abs(matrix(Sp[k], herm=False))**2
        K += w[k]*(z2 + 0.5*(p2 + p2.T))
    return ESRSpectrum(e, Vm, Dm, K, sigma=sigma, route=route,
                       reached=reached)


def _mpo_site_arrays(chain, M):
    """The site arrays (left link, bra, ket, right link) of an MPO of the
    chain's session (v3 or "python")"""
    if chain.itensor_version == "python":
        from ..pyitensor import task
        return task.mpo_arrays(M)
    return chain._session.mpo_site_arrays(M)


def _tangent_to_mps(ts, T):
    """The tangent vector sum_l A_1..A_{l-1} T_l B_{l+1}..B_L of the
    tangent space ts as MPS site arrays of twice the bond dimension: the
    bond carries whether T has been placed to the left (then B follows) or
    not (then A), so that site l is [[B_l, 0], [T_l, A_l]]."""
    from ..pyitensor import backend
    As = [np.asarray(backend.to_host(a)) for a in ts.As]
    Bs = [np.asarray(backend.to_host(b)) for b in ts.Bs]
    Ts = [np.asarray(backend.to_host(x)) for x in T]
    L = len(As)
    if L == 1:
        return [Ts[0]]
    out = []
    for l in range(L):
        A, B, X = As[l], Bs[l], Ts[l]
        a, d, b = A.shape
        dt = np.result_type(A, B, X)
        if l == 0:
            M = np.zeros((1, d, 2*b), dtype=dt)
            M[:, :, :b] = X
            M[:, :, b:] = A
        elif l == L - 1:
            M = np.zeros((2*a, d, 1), dtype=dt)
            M[:a] = B
            M[a:] = X
        else:
            M = np.zeros((2*a, d, 2*b), dtype=dt)
            M[:a, :, :b] = B
            M[a:, :, :b] = X
            M[a:, :, b:] = A
        out.append(M)
    return out


def mps_krylov_spectrum(chain, V, D, emax, epop, window, kT, nkry, ritz_tol,
                        nstates, site_weights=None, read=None):
    """ESRSpectrum on DMRG from the eigenstates below the window and the
    final states of the window lines found in the tangent space (TASK).

    The nstates lowest DMRG eigenstates must reach max(epop, window[0]):
    the populated states, and the states below the window, into which the
    upper states of the lines relax (without the one empty state below
    the window of a 14-spin chain, the widths were 5% off and the curve
    8%). From every populated state m, a block Krylov space of nkry steps
    of P(H - E_m)P in the one-site tangent space of m, started from P V|m>
    (pyitensor/task.py), gives Ritz states; those on a line of the window
    (|<n|V|m>|^2 above LINE_TOL of the strongest line out of m) must have a
    tangent-space residual below ritz_tol, or the call raises and asks for
    more nkry. They are formed as MPS (_tangent_to_mps, truncated to maxm
    by a sum), H is diagonalized in the span of all of them and of the low
    states (generalized eigenproblem, linearly dependent directions
    dropped), and the eigenstates up to emax are formed as MPS sums and
    read as on the eigenstate route (_states_spectrum), sigma included for
    the states below the window and the final states of the window lines
    (the other eigenstates of the basis are relaxation partners, with
    sigma set to 0, as on ED).
    On 14 spins at J = 30 GHz, T = 0.05 K it agrees with ED to 3e-7 of the
    peak and takes a fifth of the time of the 58 eigenstates."""
    from ..pyitensor import task, backend
    from ..mps import MPS
    from ..multioperatortk.staticoperator import StaticOperator
    if chain.itensor_version not in (3, "python"):
        raise NotImplementedError(
            "get_esr_spectrum(mode=\"DMRG\", submode=\"KRYLOV\") runs on "
            "itensor_version=3 and \"python\" only (got %r)"
            % (chain.itensor_version,))
    # P V|m> is exact for a sum of terms each on one MPS site (task.py's
    # project_operator); a term on two (the spin of a spinful fermion spans
    # its two spin-orbital sites) would start the space from a projection
    # and could miss final states without a word
    if any(len({f[1] for f in fs}) > 1 for _, fs in V.to_terms()):
        raise NotImplementedError(
            "get_esr_spectrum(mode=\"DMRG\", submode=\"KRYLOV\"): the drive "
            "has terms on more than one MPS site (the spin of a spinful "
            "fermion spans its two spin-orbital sites), so its start in the "
            "tangent space would only be a projection of V|m>; use "
            "submode=None (the lowest eigenstates) or mode=\"ED\"")
    s = chain._session
    ecut = max(epop, window[0])
    es, ws = chain.get_excited_states(n=int(nstates))
    es = np.real(np.asarray(es))
    o = np.argsort(es)
    es = es[o]
    ws = [ws[i] for i in o]
    e = es - es[0]
    if len(ws) < nstates or not e[-1] > ecut:
        raise ValueError(
            "get_esr_spectrum(mode=\"DMRG\", submode=\"KRYLOV\"): the %d "
            "lowest states reach %.4g above the ground state, but every "
            "state below the window and every populated state, up to %.4g, "
            "is needed; raise nstates" % (len(ws), e[-1], ecut))
    low = [w*(1./np.sqrt(abs(chain.overlap(w, w))))
           for w, x in zip(ws, e) if x <= ecut]
    elow = e[e <= ecut]
    Ws = _mpo_site_arrays(chain, s.H if chain.itensor_version == "python"
                          else s.hamiltonian_mpo())
    Wv = _mpo_site_arrays(chain, s.build_operator(V.to_terms()))
    # the tangent-space Krylov space of every populated state
    spaces = []
    for m, x in zip(low, elow):
        if x > epop: continue
        ts = task.TangentSpace(s.mps_site_arrays(m.cpp_handle), Ws)
        b = ts.project_operator(Wv)
        store, cnt, Hq, (jl, R) = task.block_krylov(ts, [b], int(nkry),
                                                    residual=True)
        th, U = np.linalg.eigh(Hq)
        res = np.zeros(cnt)
        if len(R):
            Ul = U[jl:, :]
            res = np.sqrt(np.abs(np.einsum("jf,jk,kf->f", Ul.conj(), R, Ul)))
        ov = np.asarray(backend.to_host(
            backend.row_overlaps(store, ts.flat(b), cnt)))[:cnt]
        strength = np.abs(U.conj().T @ ov)**2
        inwin = ((np.abs(th) + res >= window[0])
                 & (np.abs(th) - res <= window[1]))
        spaces.append((ts, store, cnt, U, res, strength*inwin))
    finals, nbad, nneed, worst = [], 0, 0, 0.
    for ts, store, cnt, U, res, strength in spaces:
        # relative to the strongest line of this populated state, so that
        # the final states of a weakly populated one are kept as well: they
        # are relaxation partners too (relative to the strongest line of
        # all, the curve was 1.9e-5 of the peak off on 10 spins, 2.3e-6 so)
        smax = np.max(strength) if len(strength) else 0.
        line = (strength > LINE_TOL*smax) & (smax > 0.)
        nneed += int(np.sum(line))
        bad = line & (res > ritz_tol)
        if np.any(bad):
            nbad += int(np.sum(bad))
            worst = max(worst, float(np.max(res[bad])))
            continue
        for k in np.flatnonzero(line):
            T = task.ritz_vector((ts, store, cnt), U[:, k])
            h = s.mps_from_site_arrays(_tangent_to_mps(ts, T))
            h = s.sum_mps_weighted([h, h], [0.5, 0.5]) # truncated to maxm
            finals.append(MPS(chain, cpp_handle=h))
    if nbad:
        raise ValueError(
            "get_esr_spectrum(mode=\"DMRG\", submode=\"KRYLOV\"): %d of the "
            "%d final states of the lines in the window are not converged "
            "in the tangent space (largest residual %.3g, ritz_tol %.3g); "
            "raise nkry (%d)" % (nbad, nneed, worst, ritz_tol, nkry))
    # one basis: the low states and the final states, which are nearly
    # parallel when two populated states reach the same final state
    basis = low + finals
    N = len(basis)
    pt = _PairTraces(chain, basis)
    S = np.zeros((N, N), dtype=complex)
    Hb = np.zeros((N, N), dtype=complex)
    Hop = StaticOperator(chain.hamiltonian, chain)
    for (a, c), G in pt.G.items():
        S[a, c] = np.trace(G[0])
        Hb[a, c] = Hop.aMb(basis[a], basis[c])
        S[c, a], Hb[c, a] = np.conj(S[a, c]), np.conj(Hb[a, c])
    sv, Us = np.linalg.eigh(S)
    ind = sv > OVERLAP_CUT*np.max(sv)
    X = Us[:, ind]/np.sqrt(sv[ind])
    E, C = np.linalg.eigh(X.conj().T @ Hb @ X)
    Cf = X @ C
    sel = np.flatnonzero(E - E[0] <= emax)
    states = [MPS(chain, cpp_handle=s.sum_mps_weighted(
        [x.cpp_handle for x in basis], [complex(c) for c in Cf[:, j]]))
        for j in sel]
    spec = _states_spectrum(chain, states, E[sel] - E[sel][0], V, D,
                            site_weights, route="KRYLOV", reached=e[-1])
    # The energy error sigma counts for the states below the window and the
    # final states of the lines that are read (read: the window and
    # lineshape.READ_MARGIN half widths around it, reached give or take the
    # error of the two energies), whose positions matter. The other
    # eigenstates of the basis are relaxation partners or ends of lines
    # further out, which act through their tails, and those above the top
    # of the margin, where the basis does not span every eigenstate, are
    # mixtures: on 12 spins (cutoff 1e-14) their errors were 0.1 to 2 GHz
    # between 35.9 and 39.2 GHz, the same for an overlap cut of 1e-10,
    # 1e-7 or 1e-4, and in the degeneracy tolerance they made two of them a
    # "level crossing". As on ED (ed_krylov_spectrum), theirs is set to 0.
    if read is None:
        read = window
    ee = spec.e
    popm = ee <= epop
    v2 = np.abs(spec.V)**2
    de = np.abs(ee[:, None] - ee[None, :])
    err = spec.sigma[:, None] + spec.sigma[None, :]
    inwin = (de + err >= read[0]) & (de - err <= read[1])
    line = np.zeros(len(ee), dtype=bool)
    for m in np.flatnonzero(popm):
        col = v2[:, m]*inwin[:, m]
        if np.max(col) > 0.:
            line |= col > LINE_TOL*np.max(col)
    spec.sigma = np.where(line | (ee <= ecut), spec.sigma, 0.)
    return spec


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
