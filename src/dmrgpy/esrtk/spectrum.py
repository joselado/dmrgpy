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

ed_krylov_spectrum (submode="KRYLOV") takes as eigenstates the populated
ones and every state below the window, up to max(epop, window[0]) above
the lowest state, epop = kT ln(1/p_cut): the upper state of a line relaxes
into them. The final states come from one block Krylov space started from
V|m> for every populated m. A second, closing block Krylov space, started
from S^a_k|n> and V|n> for the states below the window and the converged
final states n of the window lines, adds the relaxation partners and the
states of a second photon that V|m> and H never reach: when V conserves a
quantum number of H that S^a_k does not (a tip moment along the field),
the whole of the other sectors. H is diagonalized in the span of the
eigenstates and of both spaces. This follows the pattern of
kondospectrumtk/krylovkondo.py, with one difference that is measured: the
relaxation needs one common basis. A two-level Bloch line per transition,
as separate Krylov spaces would give, was 46 to 61% of the peak off the
full rate equation on 8 and 10 spins, while the full rate equation
restricted to the populated states and the final states with |V_nm|^2
above 1e-2 of the largest was within 4e-7 to 2e-4 of it. The final states
of the lines in the window must converge (ed_krylov_spectrum); the other
Ritz states are kept as relaxation partners.

Both Krylov routes end with a closure guard (_ed_guard, _mps_guard): for
each final state n of a window line, the part of S^a_k|n> and of V|n>
outside the kept states is resolved in energy by a few Lanczos steps, and
the call raises, naming submode=None, when a part of it below emax carries
more than GUARD_TOL of the relaxation rate out of n, or of |V|n>|^2. The
DMRG route has no closing stage: its states below the window come from
DMRG, and above the window the guard is what stands for the partners
whenever H conserves a charge that the spin operators change, the case in
which whole sectors are out of reach of its tangent spaces.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from ..algebra import algebra
from ..kondospectrumtk import edkondo
from . import rates

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
# submode="KRYLOV" on DMRG: a relaxation partner is dropped as a copy when
# more than COPY_WEIGHT of its norm squared lies in the span of the vectors
# already in the basis (the states below the window, the final states of
# the lines, the partners kept before it) within ritz_tol of its energy
# (_drop_copies). The span, not the largest single overlap, decides: two
# levels A and B closer than ritz_tol that a tangent space returns rotated,
# (A+B)/sqrt2 and (A-B)/sqrt2, each have an overlap of 0.71 with a line at
# A, and a test on single overlaps dropped both and lost B. The weights of
# the two members in the span of A add up to 1, so any COPY_WEIGHT in
# [1/2, 1) keeps B; at 1/2 itself rounding dropped both members of a pair
# at equal mixing, and every copy measured on 10 and 14 spins had a weight
# of 0.994 or more.
COPY_WEIGHT = 0.9
# the pseudo-inverse of the Gram matrix of those vectors drops directions
# below SPAN_RCOND of the largest (two populated states reaching one final
# state give two nearly parallel lines)
SPAN_RCOND = 1e-10
# submode="KRYLOV": a line out of a populated state m into n counts as a line
# of the window when |V_nm|^2 p_m is above LINE_TOL of the strongest
LINE_TOL = 1e-6
# submode="KRYLOV" on ED: the closing stage starts from S^z_k|n>, S^+_k|n>,
# S^-_k|n> and V|n> for every site k with a weight and every state n below
# the window or final state of a window line, when there are at most
# CLOSE_WIDTH of these vectors, and from CLOSE_WIDTH random combinations of
# them otherwise (generator seeded with CLOSE_SEED): a combination reaches
# every eigenstate the vectors reach, and the Krylov space resolves them. It
# grows one block at a time, and checks its Ritz states each time it has
# grown by the factor CLOSE_GROWTH.
CLOSE_WIDTH = 128
CLOSE_SEED = 20261121
CLOSE_GROWTH = 1.25
# The closure guard of both Krylov routes. For each final state n of a
# window line, the part of W|n> outside the kept states, W a witness of the
# Korringa operators, is resolved in energy by GUARD_STEPS Lanczos steps of
# H. The call raises when the part below emax takes more than GUARD_TOL of
# the relaxation rate out of n (and is more than GUARD_TOL of |W|n>|^2),
# or when the part of V|n> below emax is more than GUARD_TOL of |V|n>|^2.
# The witnesses are sqrt(w_k) S^z_k and sqrt(w_k/2) S^+-_k themselves when
# there are at most GUARD_WITNESSES of them, and otherwise GUARD_WITNESSES
# random combinations of them (seed GUARD_SEED), whose mean |<m|W|n>|^2 is
# the Korringa kernel K[m, n]. A missing T+ on the guide's dimer with the
# tip along the field moves the curve by 1.2e-4 of its peak, through the
# readout D_T+T+ = 1/2, and takes about 6e-6 of the rate out of T0 (worked
# out by hand from the dimer's eigenstates, not measured). A few Lanczos
# steps give an upper bound of the lowest missing energy, so the guard
# sees a missing state below emax for sure only when its weight brings a
# Ritz value below emax in GUARD_STEPS steps. The chains that pin it
# (tests/test_audit_2026_11_21_esr_krylov.py) exhaust their Hilbert space.
# On DMRG the guard runs only when H conserves a charge that S^+_k changes
# (_sector_blocked).
GUARD_TOL = 1e-6
GUARD_STEPS = 6
GUARD_WITNESSES = 6
GUARD_SEED = 20261121
# A Lanczos recurrence of the guard stops when its next coupling is below
# LANCZOS_BREAK of the largest entry of its tridiagonal matrix so far, and
# the ED guard runs its Lanczos recurrences GUARD_CHUNK columns at a time
LANCZOS_BREAK = 1e-10
GUARD_CHUNK = 16


class ESRSpectrum():
    """Energies (relative to the lowest state), V, D and K in the basis of
    the kept eigenstates. sigma bounds the error of each energy: the square
    root of <H^2> - <H>^2 on the MPS route, the Ritz residual on the Krylov
    route (None on ED, where the states are exact to rounding). The Krylov
    routes set it to 0 for relaxation partners, so that it takes part in the
    degeneracy rule and the position warning for the lines that are read
    only; sigma_measured keeps the measured one for every state, which
    return_poles reports."""
    def __init__(self, e, V, D, K, sigma=None, route="ED", reached=None):
        self.e = np.asarray(e, dtype=float)
        self.dim = len(self.e)
        self.V = np.asarray(V)
        self.D = np.asarray(D)
        self.K = np.asarray(K, dtype=float)
        self.sigma = None if sigma is None else np.asarray(sigma, dtype=float)
        self.sigma_measured = self.sigma
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
    site_matrices = []
    for k in range(len(w)):
        if w[k] == 0.: continue
        z = matrix(Sz[k])
        p = matrix(Sp[k], herm=False)
        site_matrices.append((k, z, p))
        p2 = np.abs(p)**2
        K += w[k]*(np.abs(z)**2 + 0.5*(p2 + p2.T))
    spec = ESRSpectrum(e, Vm, Dm, K, sigma=sigma, route=route,
                       reached=reached)
    # the normalized states and <a|S^z_k|b>, <a|S^+_k|b> between them, for
    # the closure guard of the Krylov route (_mps_guard)
    spec.states = ws
    spec.site_matrices = site_matrices
    return spec


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


def _span_weight(chain, Q, p):
    """The weight of the unit MPS p in the span of the MPS Q"""
    G = np.array([[chain.overlap(a, b) for b in Q] for a in Q])
    c = np.array([chain.overlap(a, p) for a in Q])
    w = c.conj() @ np.linalg.pinv(G, rcond=SPAN_RCOND, hermitian=True) @ c
    return float(np.real(w))


def _drop_copies(chain, low, elow, finals, isline, efin, tol):
    """The final states of the lines (finals with isline), then the
    relaxation partners that are not copies of the vectors already kept: a
    partner with more than COPY_WEIGHT of its weight in the span of the
    states below the window, the lines and the partners kept before it
    within tol of its energy (efin, elow) is dropped. The tangent spaces
    of the populated states reach the same eigenstates, so most partners
    are copies: 52 of 60 on 14 spins (J = 30 GHz, T = 0.05 K), many of them
    poorly converged (true residuals up to 0.6 GHz). They left the common
    eigenproblem ill-conditioned and moved V and D of the ground state by
    3e-6 (an overlap cut of 1e-12 or 1e-8 instead of 1e-10 gave 6.3e-7 and
    4e-5 of the peak; keeping only the 23 partners with a true residual
    below ritz_tol, mostly well converged copies, 1.6e-7). With every
    partner the curve was 1.1e-6 of the peak off ED in 306 s, without the
    copies 1.6e-7 in 216 s, and 3.3e-7 without any partner, which misses
    two of the 58 states ED keeps; on 10 spins with no tip field in H,
    1.1e-6 and 4.8e-7. On both chains every copy had a weight of 0.994 or
    more and every partner kept had no vector within tol of its energy, so
    COPY_WEIGHT decides nothing there. The low states are already
    normalized."""
    def unit(w):
        return w*(1./np.sqrt(abs(chain.overlap(w, w))))
    ref = list(zip(low, elow))
    ref += [(unit(f), x) for f, l, x in zip(finals, isline, efin) if l]
    out = [f for f, l in zip(finals, isline) if l]
    for f, l, x in zip(finals, isline, efin):
        if l: continue
        fu = unit(f)
        near = [q for q, y in ref if abs(y - x) < tol]
        if near and _span_weight(chain, near, fu) > COPY_WEIGHT:
            continue
        out.append(f)
        ref.append((fu, x))
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
    by a sum), together with the converged Ritz states between the window
    bottom and emax as relaxation partners, those that copy a vector
    already kept dropped (_drop_copies). H is diagonalized in the span of
    all of them and of the low states (generalized eigenproblem, linearly
    dependent directions dropped), and the eigenstates up to emax are
    formed as MPS sums and read as on the eigenstate route
    (_states_spectrum), sigma included for
    the states below the window and the final states of the window lines
    (the other eigenstates of the basis are relaxation partners, with
    sigma set to 0, as on ED).
    On 14 spins at J = 30 GHz, T = 0.05 K with a tip field it agrees with
    ED to 1.6e-7 of the peak, in 216 s on one core (189 s and 3.3e-7
    without the relaxation partners, which miss two of the 58 states).

    The partners above the window are only those the tangent spaces of the
    populated states reach from V|m>, so a sector that V and H never reach
    (the T+ and T- of an S-T0 line when the tip moment is along the field)
    is missing from them. There is no closing stage here, as there is on
    ED: seeding the tangent-space Krylov space with S^a_k|n> would take a
    start block of 3 vectors per site and final state, each with a store of
    nkry steps. _mps_guard checks the closure instead and raises, naming
    submode=None."""
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
    # the tangent-space Krylov space of every populated state, one at a
    # time: a store of nkry steps is about 0.5 GB on 14 spins at maxm=64,
    # so each is dropped before the next is built
    finals, isline, efin, nbad, nneed, worst = [], [], [], 0, 0, 0.
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
        strength = strength*inwin
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
        elif not nbad:
            # the converged Ritz states between the window bottom and emax
            # too, as relaxation partners: with the window lines alone a
            # chain without a tip field in H was 4.3e-5 of the peak off on
            # 10 spins (the states V never reaches), 1.1e-6 with them, for
            # 15% more time (a review of 2026-11-21); those that copy a
            # vector already in the basis are dropped after the loop
            # (_drop_copies)
            keep = line | ((res <= ritz_tol) & (x + th <= emax)
                           & (x + th > ecut))
            for k in np.flatnonzero(keep):
                T = task.ritz_vector((ts, store, cnt), U[:, k])
                h = s.mps_from_site_arrays(_tangent_to_mps(ts, T))
                h = s.sum_mps_weighted([h, h], [0.5, 0.5]) # truncated to maxm
                finals.append(MPS(chain, cpp_handle=h))
                isline.append(bool(line[k]))
                efin.append(float(x + th[k]))
        del ts, store, b
    if nbad:
        raise ValueError(
            "get_esr_spectrum(mode=\"DMRG\", submode=\"KRYLOV\"): %d of the "
            "%d final states of the lines in the window are not converged "
            "in the tangent space (largest residual %.3g, ritz_tol %.3g); "
            "raise nkry (%d)" % (nbad, nneed, worst, ritz_tol, nkry))
    finals = _drop_copies(chain, low, elow, finals, isline, efin, ritz_tol)
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
    line = _line_finals(ee, spec.V, spec.sigma, popm, read)
    finals = _line_finals(ee, spec.V, spec.sigma, popm, window)
    spec.sigma_measured = spec.sigma
    spec.sigma = np.where(line | (ee <= ecut), spec.sigma, 0.)
    _mps_guard(chain, spec, Hop, E[sel][0], finals, V,
               _weights(chain, site_weights), emax, kT)
    return spec


def _line_finals(e, V, err, popm, rng):
    """The final states of the lines out of the populated states (popm)
    whose position, give or take the error err of its two energies, lies in
    rng = (lo, hi), with |V_nm|^2 above LINE_TOL of the strongest such line
    out of the same m"""
    v2 = np.abs(V)**2
    de = np.abs(e[:, None] - e[None, :])
    er = err[:, None] + err[None, :]
    inwin = (de + er >= rng[0]) & (de - er <= rng[1])
    line = np.zeros(len(e), dtype=bool)
    for m in np.flatnonzero(popm):
        col = v2[:, m]*inwin[:, m]
        if np.max(col) > 0.:
            line |= col > LINE_TOL*np.max(col)
    return line


def _H(X):
    """The adjoint of a dense array, without a copy when it is real"""
    return X.conj().T if np.iscomplexobj(X) else X.T


def _witnesses(wk):
    """The witness operators of the closure guard for the weights wk (> 0)
    of the sites that relax: a list of (n, 3) arrays, the coefficients of
    each witness on (S^z_k, S^+_k, S^-_k) of those sites, and the factor
    that turns the sum of the measures of the witnesses into the measure of
    the Korringa kernel, sum_m K[m, n] delta(E - E_m) (see GUARD_TOL)"""
    wk = np.asarray(wk, dtype=float)
    scale = np.sqrt(np.column_stack([wk, 0.5*wk, 0.5*wk]))
    if scale.size <= GUARD_WITNESSES:
        out = []
        for i, a in np.ndindex(scale.shape):
            c = np.zeros(scale.shape)
            c[i, a] = scale[i, a]
            out.append(c)
        return out, 1.
    rng = np.random.default_rng(GUARD_SEED)
    out = [scale*rng.standard_normal(scale.shape)
           for _ in range(GUARD_WITNESSES)]
    return out, 1./GUARD_WITNESSES


def _guard_verdict(e, K, idx, kT, norm, miss_w, miss_r, nv, miss_v, label):
    """Raise when a final state of a window line sends more than GUARD_TOL
    of its relaxation rate into states below emax outside the basis, or
    when a second photon from it reaches such states; idx are the finals,
    norm, miss_w and miss_r the weight of the witnesses, and the weight and
    the rate of their parts outside the basis below emax, nv and miss_v the
    same weights for V|n>."""
    kept = np.zeros(len(idx))
    for j, n in enumerate(idx):
        others = np.arange(len(e)) != n
        kept[j] = np.sum(K[others, n]*rates.phi(e[n] - e[others], kT))
    bad = (miss_w > GUARD_TOL*norm) & (miss_r > GUARD_TOL*(kept + miss_r))
    if np.any(bad):
        share = miss_r/np.maximum(kept + miss_r, 1e-300)
        j = np.flatnonzero(bad)[np.argmax(share[bad])]
        raise ValueError(
            "get_esr_spectrum(%s): the final state of a line in the window "
            "at %.6g above the lowest state (in the units of H) sends %.3g of "
            "its relaxation rate into states below emax that the Krylov basis "
            "does not hold: the space reached from V|m> misses relaxation "
            "partners, for instance the other sectors of a quantum number of "
            "H that the drive conserves and the spin operators of the "
            "relaxation do not (a tip moment along the field). Use "
            "submode=None, the eigenstates up to the top of the window"
            % (label, e[idx[j]], share[j]))
    badv = miss_v > GUARD_TOL*nv
    if np.any(badv):
        share = miss_v/np.maximum(nv, 1e-300)
        j = np.flatnonzero(badv)[np.argmax(share[badv])]
        raise ValueError(
            "get_esr_spectrum(%s): the drive takes the final state of a line "
            "in the window at %.6g above the lowest state (in the units of H) "
            "into states below emax that the Krylov basis does not hold (%.3g "
            "of |V|n>|^2), the upper states of a second photon. Use "
            "submode=None, the eigenstates up to the top of the window"
            % (label, e[idx[j]], share[j]))


def _below_emax(th, om, e0, en, emax, kT):
    """The weight and the relaxation rate from a state at en (both relative
    to the lowest state) of the Ritz values th (absolute) and weights om of
    a Lanczos measure, counting those at most emax above the lowest state,
    e0 its absolute energy"""
    x = np.asarray(th) - e0
    below = x <= emax
    if not np.any(below):
        return 0., 0.
    om = np.asarray(om)[below]
    return float(np.sum(om)), float(np.sum(om*rates.phi(en - x[below], kT)))


def _mps_lanczos(chain, Hop, r, steps):
    """(theta, weight): the Ritz values of `steps` Lanczos steps of H (the
    StaticOperator Hop) from the MPS r, and their Gauss weights times
    ||r||^2, the quadrature of the spectral measure of r. Each new vector
    is orthogonalized against the two before it, as a weighted MPS sum
    truncated to the chain's maxm."""
    from ..mps import MPS
    s = chain._session
    n2 = float(np.real(chain.overlap(r, r)))
    if not n2 > 0.:
        return np.zeros(0), np.zeros(0)
    v = r*(1./np.sqrt(n2))
    prev, b = None, 0.
    alpha, beta = [], []
    for j in range(int(steps)):
        x = Hop*v
        a = float(np.real(chain.overlap(v, x)))
        alpha.append(a)
        if j == int(steps) - 1:
            break
        hs, cs = [x.cpp_handle, v.cpp_handle], [1.+0j, complex(-a)]
        if prev is not None:
            hs.append(prev.cpp_handle)
            cs.append(complex(-b))
        x = MPS(chain, cpp_handle=s.sum_mps_weighted(hs, cs))
        b = np.sqrt(max(float(np.real(chain.overlap(x, x))), 0.))
        top = max(np.max(np.abs(alpha)), max(beta) if beta else 0.)
        if not b > LANCZOS_BREAK*top:
            break
        beta.append(b)
        prev, v = v, x*(1./b)
    m = len(alpha)
    T = (np.diag(alpha) + np.diag(beta[:m - 1], 1)
         + np.diag(beta[:m - 1], -1))
    th, Y = np.linalg.eigh(T)
    return th, n2*np.abs(Y[0])**2


def _sector_blocked(chain, ops):
    """Whether H conserves a charge of the site types (Sz, Nf) that one of
    the operators ops changes: the tangent spaces of the populated states
    then stay in the sectors V reaches from them, and the relaxation
    partners in the sectors next to those are never in them"""
    from ..multioperatortk import charge
    names = charge.conserved_qn_names(chain, chain.hamiltonian)
    if not names:
        return False
    for op in ops:
        q = charge.multioperator_charge(chain, op, names)
        if q is None or any(x != 0 for x in q):
            return True
    return False


def _mps_guard(chain, spec, Hop, e0, finals, V, w, emax, kT):
    """The closure guard of mps_krylov_spectrum (see GUARD_TOL), with MPS:
    for each final state n of a window line, W|n> for every witness W and
    V|n>, minus their parts on the kept states (whose overlaps are the
    matrices the spectrum was read from), resolved in energy by GUARD_STEPS
    Lanczos steps. w are the site weights and e0 the absolute energy of the
    lowest kept state.

    It runs only when H conserves a charge of the site types that S^+_k
    changes (_sector_blocked), the case where whole sectors of partners
    are out of reach. Without such a charge (a transverse field anywhere in
    H) the tangent spaces reach the partners, and those above the margin
    are kept when converged and dropped otherwise, by design (1.1e-6 of the
    peak on 10 spins); the guard would count the dropped ones as missing,
    and how often that would refuse a chain of the plan is not measured.
    A conserved quantity that is not a charge of the site types (Sx_tot
    for a field and a tip along x, a mirror of the chain) is not seen."""
    from ..mps import MPS
    from ..multioperatortk.staticoperator import StaticOperator
    idx = np.flatnonzero(finals)
    if len(idx) == 0:
        return
    s = chain._session
    states = spec.states
    ks = [k for k, _, _ in spec.site_matrices]
    Sz, Sp = site_operators(chain)
    if not _sector_blocked(chain, [Sp[k] for k in ks]):
        return
    wits, factor = _witnesses(np.asarray(w)[ks])

    def outside(Op, cols, n):
        """|phi|^2, and the weight and the rate below emax of the part of
        phi = Op|n> outside the kept states, cols[i] = <i|Op|n>"""
        phi = Op*states[n]
        n2 = float(np.real(chain.overlap(phi, phi)))
        big = np.flatnonzero(np.abs(cols)**2 > 1e-14*n2)
        if n2 - np.sum(np.abs(cols[big])**2) <= 1e-3*GUARD_TOL*n2:
            return n2, 0., 0.
        r = MPS(chain, cpp_handle=s.sum_mps_weighted(
            [phi.cpp_handle] + [states[i].cpp_handle for i in big],
            [1.+0j] + [complex(-cols[i]) for i in big]))
        th, om = _mps_lanczos(chain, Hop, r, GUARD_STEPS)
        return (n2,) + _below_emax(th, om, e0, spec.e[n], emax, kT)

    norm = np.zeros(len(idx))
    miss_w = np.zeros(len(idx))
    miss_r = np.zeros(len(idx))
    for c in wits:
        O = 0
        M = 0.
        for i, (k, z, p) in enumerate(spec.site_matrices):
            Sm = chain.Sx[k] - 1j*chain.Sy[k]
            for a, (op, mat) in enumerate(((Sz[k], z), (Sp[k], p),
                                           (Sm, p.conj().T))):
                if c[i, a] != 0.:
                    O = O + float(c[i, a])*op
                    M = M + float(c[i, a])*mat
        Ow = StaticOperator(O, chain)
        for j, n in enumerate(idx):
            n2, mw, mr = outside(Ow, M[:, n], n)
            norm[j] += factor*n2
            miss_w[j] += factor*mw
            miss_r[j] += factor*mr
    nv = np.zeros(len(idx))
    miss_v = np.zeros(len(idx))
    Vop = StaticOperator(V, chain)
    for j, n in enumerate(idx):
        nv[j], miss_v[j], _ = outside(Vop, spec.V[:, n], n)
    _guard_verdict(spec.e, spec.K, idx, kT, norm, miss_w, miss_r, nv, miss_v,
                   "mode=\"DMRG\", submode=\"KRYLOV\"")


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


def _rayleigh_ritz(Hs, Q, emax, width):
    """(theta, y, cand, res): the Ritz values theta (sorted) of H in the
    span of the orthonormal columns of Q, the Ritz vectors being Q @ y, and
    of the Ritz states up to emax + width above the lowest the ones whose
    residual ||(H - theta)|u>|| reaches below emax (an eigenvalue lies
    within it), with those residuals. H Q is formed CHUNK columns at a time
    and not kept."""
    N = Q.shape[1]
    Hk = np.empty((N, N), dtype=np.result_type(Q.dtype, Hs.dtype))
    QH = _H(Q)
    for c0 in range(0, N, CHUNK):
        Hk[:, c0:c0 + CHUNK] = QH @ (Hs @ Q[:, c0:c0 + CHUNK])
    theta, y = np.linalg.eigh(0.5*(Hk + _H(Hk)))
    del Hk
    cand = np.flatnonzero(theta - theta[0] <= emax + width)
    res = np.empty(len(cand))
    for c0 in range(0, len(cand), CHUNK):
        sel = cand[c0:c0 + CHUNK]
        U = Q @ y[:, sel]
        R = Hs @ U - U*theta[sel][None, :]
        res[c0:c0 + CHUNK] = np.linalg.norm(R, axis=0)
    inside = theta[cand] - theta[0] - res <= emax
    return theta, y, cand[inside], res[inside]


def _ed_lines(e, Vm, res, epop, kT, window):
    """(popm, need): the populated states (up to epop above the lowest) and
    the final states of the lines in the window out of them, those n with
    |V_nm|^2 p_m above LINE_TOL of the strongest line and |E_n - E_m| in
    the window give or take the residual of n"""
    popm = e <= epop
    p = np.zeros(len(e))
    p[popm] = np.exp(-e[popm]/kT) if kT > 0. else (e[popm] == 0.)
    strength = np.abs(Vm)**2*p[None, :] # zero unless m is populated
    de = np.abs(e[:, None] - e[None, :])
    inwin = (de + res[:, None] >= window[0]) & (de - res[:, None] <= window[1])
    smax = np.max(strength*inwin) if strength.size else 0.
    line = inwin & (strength > LINE_TOL*smax) & (smax > 0.)
    return popm, np.any(line, axis=1)


def _check_finals(need, res, ritz_tol, nkry):
    bad = need & (res > ritz_tol)
    if np.any(bad):
        raise ValueError(
            "get_esr_spectrum(submode=\"KRYLOV\"): %d of the %d final states "
            "of the lines in the window are not converged (largest residual "
            "%.3g, ritz_tol %.3g); raise nkry (%d)"
            % (np.sum(bad), np.sum(need), np.max(res[bad]), ritz_tol, nkry))


def _sparse_operator(edobj, op):
    """The ED matrix of op as a csr matrix, real when its imaginary part is
    exactly zero"""
    M = sp.csr_matrix(edobj.get_operator(op))
    if edkondo._imag_is_zero(M):
        return _part(M, np.real)
    return sp.csr_matrix(M, dtype=complex)


def _closing_start(mats, X):
    """The start of the closing stage: M X for every matrix of mats, or
    CLOSE_WIDTH random combinations of their columns when there are more
    than CLOSE_WIDTH of them (see CLOSE_WIDTH)"""
    if len(mats)*X.shape[1] <= CLOSE_WIDTH:
        return np.hstack([M @ X for M in mats])
    rng = np.random.default_rng(CLOSE_SEED)
    out = 0.
    for M in mats:
        out = out + M @ (X @ rng.standard_normal((X.shape[1], CLOSE_WIDTH)))
    return out


def _project_out(X, Q1, blocks):
    """X projected twice on the complement of the columns of Q1 and of the
    blocks, and its coefficients on the blocks, stacked"""
    dt = np.result_type(X.dtype, Q1.dtype, *[b.dtype for b in blocks])
    X = X.astype(dt, copy=False)
    c = np.zeros((sum(b.shape[1] for b in blocks), X.shape[1]), dtype=dt)
    for _ in range(2):
        X = X - Q1 @ (_H(Q1) @ X)
        i = 0
        for b in blocks:
            d = _H(b) @ X
            X = X - b @ d
            c[i:i + b.shape[1]] += d
            i += b.shape[1]
    return X, c


def _closing_krylov(Hs, Q1, S0, nkry, top, tol):
    """An orthonormal basis Q2, orthogonal to the columns of Q1, of the
    block Krylov space of the compression of H to the complement of Q1,
    started from the columns of S0 (their part above DROP_TOL, _orth). It
    grows one block at a time, H B projected twice on the complement of Q1
    and of the space so far (the coefficients on the space are its
    projected matrix, as in krylovkondo.block_krylov), and stops when a
    block is empty (the space is exhausted), after nkry blocks, or once
    every Ritz state of the compression with an energy up to top (absolute),
    and the lowest one above it, has a residual below tol: the states the
    start reaches below emax are then converged, and more blocks would
    refine the ones above it. The residual is the norm of what the last
    block leaves outside the space, and it is checked each time the space
    has grown by CLOSE_GROWTH."""
    B = _orth(S0, Q1)
    blocks = []
    Hk = np.zeros((0, 0), dtype=np.result_type(Hs.dtype, Q1.dtype, B.dtype))
    n = checked = 0
    coupling = None
    for _ in range(int(nkry)):
        if B.shape[1] == 0: break # the space is exhausted
        k = B.shape[1]
        Hk = np.pad(Hk, ((0, k), (0, k)))
        if coupling is not None:
            Hk[n:n + k, n - coupling.shape[1]:n] = coupling
        blocks.append(B)
        n += k
        X = Hs @ B
        scale = float(np.max(np.linalg.norm(X, axis=0)))
        X, c = _project_out(X, Q1, blocks)
        Hk[:n, n - k:n] = c
        if scale > 0.:
            Ux, s, Vh = np.linalg.svd(X, full_matrices=False)
            keep = s > DROP_TOL*scale
            W = s[:, None]*Vh
        else:
            Ux, keep, W = X[:, :0], np.zeros(0, dtype=bool), np.zeros((0, k))
        B, coupling = Ux[:, keep], W[keep]
        if n >= CLOSE_GROWTH*checked:
            checked = n
            th, yk = np.linalg.eigh(0.5*(Hk + _H(Hk)))
            # the Ritz values come down from above as the space grows, so
            # the lowest one above top must be converged too: a start block
            # with most of its weight high up has none below top at first
            low = th <= top
            if not np.all(low):
                low[np.argmax(~low)] = True
            if np.all(np.linalg.norm(W @ yk[n - k:, low], axis=0) < tol):
                break
    if not blocks:
        return Q1[:, :0]
    return np.hstack(blocks)


def _lanczos_columns(Hs, R, steps):
    """For every column r of R: (theta, weight), the Ritz values of `steps`
    Lanczos steps of H from r, with full reorthogonalization, and their
    Gauss weights times ||r||^2 (the quadrature of the spectral measure of
    r); a recurrence stops at a coupling below LANCZOS_BREAK."""
    nrm = np.linalg.norm(R, axis=0)
    X = R/np.where(nrm > 0., nrm, 1.)[None, :]
    Vs = []
    alpha, beta = [], []
    for _ in range(int(steps)):
        Vs.append(X)
        W = Hs @ X
        alpha.append(np.real(np.sum(X.conj()*W, axis=0)))
        for _ in range(2):
            for Y in Vs:
                W = W - Y*np.sum(Y.conj()*W, axis=0)[None, :]
        b = np.linalg.norm(W, axis=0)
        beta.append(b)
        X = W/np.where(b > 0., b, 1.)[None, :]
    alpha, beta = np.array(alpha), np.array(beta)
    out = []
    for c in range(R.shape[1]):
        if not nrm[c] > 0.:
            out.append((np.zeros(0), np.zeros(0)))
            continue
        m = len(alpha)
        for j in range(len(alpha) - 1):
            top = max(np.max(np.abs(alpha[:j + 1, c])), np.max(beta[:j, c])
                      if j else 0.)
            if not beta[j, c] > LANCZOS_BREAK*top:
                m = j + 1
                break
        T = (np.diag(alpha[:m, c]) + np.diag(beta[:m - 1, c], 1)
             + np.diag(beta[:m - 1, c], -1))
        th, Y = np.linalg.eigh(T)
        out.append((th, nrm[c]**2*np.abs(Y[0])**2))
    return out


def _ed_guard(Hs, Vs, ops, U, e, e0, K, finals, emax, kT):
    """The closure guard of ed_krylov_spectrum (see GUARD_TOL): for each
    final state n of a window line, W|n> for every witness W and V|n>,
    minus their parts on the kept Ritz states U, resolved in energy by
    GUARD_STEPS Lanczos steps. ops are (k, w_k, S^z_k, S^+_k, S^-_k) of the
    sites that relax, and e0 the absolute energy of the lowest kept state."""
    idx = np.flatnonzero(finals)
    if len(idx) == 0:
        return
    Uf = U[:, idx]
    UH = _H(U)

    def outside(Phi):
        """per column: |phi|^2, and the weight and the rate below emax of
        the part of phi outside the kept states"""
        n2 = np.sum(np.abs(Phi)**2, axis=0)
        mw = np.zeros(len(idx))
        mr = np.zeros(len(idx))
        for c0 in range(0, len(idx), GUARD_CHUNK):
            R = Phi[:, c0:c0 + GUARD_CHUNK]
            for _ in range(2):
                R = R - U @ (UH @ R)
            r2 = np.sum(np.abs(R)**2, axis=0)
            run = np.flatnonzero(r2 > 1e-3*GUARD_TOL*n2[c0:c0 + GUARD_CHUNK])
            if len(run) == 0: continue
            for j, (th, om) in zip(run, _lanczos_columns(Hs, R[:, run],
                                                         GUARD_STEPS)):
                mw[c0 + j], mr[c0 + j] = _below_emax(th, om, e0, e[idx[c0 + j]],
                                                     emax, kT)
        return n2, mw, mr

    wits, factor = _witnesses([wk for _, wk, _, _, _ in ops])
    norm = np.zeros(len(idx))
    miss_w = np.zeros(len(idx))
    miss_r = np.zeros(len(idx))
    for c in wits:
        O = None
        for i, (_, _, z, p, m) in enumerate(ops):
            for a, M in enumerate((z, p, m)):
                if c[i, a] != 0.:
                    t = float(c[i, a])*M
                    O = t if O is None else O + t
        n2, mw, mr = outside(O @ Uf)
        norm += factor*n2
        miss_w += factor*mw
        miss_r += factor*mr
    nv, miss_v, _ = outside(Vs @ Uf)
    _guard_verdict(e, K, idx, kT, norm, miss_w, miss_r, nv, miss_v,
                   "submode=\"KRYLOV\"")


def ed_krylov_spectrum(chain, V, D, emax, epop, window, kT, nkry, ritz_tol,
                       site_weights=None):
    """ESRSpectrum from the eigenstates up to max(epop, window[0]) above the
    lowest (the populated ones, up to epop, and every state below the
    window), the Ritz states of one block Krylov space started from V|m>
    for every populated m, and those of a closing block Krylov space, kept
    up to emax above the lowest state. window = (elo, ehi) is the energy
    range of the lines that are read (the frequency window and its margin,
    times h). Both spaces have at most nkry blocks.

    The final states of the lines in the window, those n with
    |V_nm|^2 p_m above LINE_TOL of the strongest for a populated m and
    |E_n - E_m| in the window, must be converged: their residual
    ||(H - theta)|u>|| below ritz_tol, so that an eigenvalue lies within
    it, or the call raises and asks for more nkry. The other Ritz states
    are kept as they are, as relaxation partners, whose energies and spin
    matrix elements are then approximate; this needs about half the Krylov
    space of asking every state below emax to converge (measured on 10 and
    12 spins), and their sigma is set to 0.

    The closing space starts from S^z_k|n>, S^+_k|n>, S^-_k|n> (every site
    k with a weight) and V|n>, for the states below the window and the
    converged final states n, orthogonal to everything before
    (_closing_start, _closing_krylov), and H is diagonalized again in the
    span of all of it. It holds the relaxation partners that the space of
    V|m> does not reach, the T+ and T- of an S-T0 line when the tip moment
    is along the field (2 of the 4 states of the dimer were kept without
    it, and its line came out -5.0 times the eigenstate one), and the upper
    states of a second photon (the |+1> of a saturated spin 1). The guard
    _ed_guard checks the closure at the end."""
    w = _weights(chain, site_weights)
    edobj = chain.get_ED_obj()
    Hs = sp.csr_matrix(edobj.get_hamiltonian())
    if edkondo._imag_is_zero(Hs):
        Hs = _part(Hs, np.real)
    dim = Hs.shape[0]
    # every eigenstate below the window, into which the upper state of a
    # line relaxes, and every populated one (a 4-spin chain with an empty
    # Sz = -2 level below a line was 2.1e-2 of its peak off without it)
    ecut = max(epop, window[0])
    if dim <= algebra.maxsize:
        emu, vs = np.linalg.eigh(Hs.toarray())
    else:
        emu, vs = _low_states(Hs, ecut)
    low = emu - emu[0] <= ecut
    lowpop = (emu - emu[0] <= epop)[low]
    P = np.ascontiguousarray(vs[:, low])
    Vs = sp.csr_matrix(edobj.get_operator(V))
    if not edkondo._imag_is_zero(Vs) or np.iscomplexobj(P):
        Vs = sp.csr_matrix(Vs, dtype=complex)
        P = P.astype(complex)
        Hs = sp.csr_matrix(Hs, dtype=complex)
    else:
        Vs = _part(Vs, np.real)
    # the block Krylov space of V|m> for the populated m, orthogonal to the
    # eigenstates
    Q = P
    B = _orth(Vs @ P[:, lowpop], Q)
    for _ in range(int(nkry)):
        if B.shape[1] == 0: break # the space is exhausted
        Q = np.hstack([Q, B])
        B = _orth(Hs @ B, Q)
    # the Ritz states up to emax, and those above whose residual reaches
    # below it, tested up to one window width above emax
    width = window[1] - window[0]
    theta, y, cand, res = _rayleigh_ritz(Hs, Q, emax, width)
    U = Q @ y[:, cand]
    e = theta[cand] - theta[cand[0]]
    Vm = _transform(Vs, U)
    popm, need = _ed_lines(e, Vm, res, epop, kT, window)
    _check_finals(need, res, ritz_tol, nkry)
    # the closing stage, from the spin operators of every site that relaxes
    # and from V, applied to the states below the window and to the final
    # states of the window lines
    Sz, Sp = site_operators(chain)
    ops = []
    for k in range(len(w)):
        if w[k] == 0.: continue
        z = _sparse_operator(edobj, Sz[k])
        p = _sparse_operator(edobj, Sp[k])
        ops.append((k, w[k], z, p, sp.csr_matrix(_H(p))))
    mats = [M for _, _, z, p, m in ops for M in (z, p, m)] + [Vs]
    S0 = _closing_start(mats, np.hstack([P, U[:, need]]))
    Q2 = _closing_krylov(Hs, Q, S0, nkry, theta[0] + emax, ritz_tol)
    del S0
    if Q2.shape[1]:
        Q = np.hstack([Q, Q2])
        del Q2
        theta, y, cand, res = _rayleigh_ritz(Hs, Q, emax, width)
        U = Q @ y[:, cand]
        e = theta[cand] - theta[cand[0]]
        Vm = _transform(Vs, U)
        popm, need = _ed_lines(e, Vm, res, epop, kT, window)
        _check_finals(need, res, ritz_tol, nkry)
    del Q
    Dm = _transform(edobj.get_operator(D), U)
    K = np.zeros((len(e), len(e)))
    for _, wk, z, p, _ in ops:
        zk = _transform(z, U)
        pk = _transform(p, U)
        p2 = np.abs(pk)**2
        K += wk*(np.abs(zk)**2 + 0.5*(p2 + p2.T))
    _ed_guard(Hs, Vs, ops, U, e, theta[cand[0]], K, need, emax, kT)
    lowk = e <= ecut
    sigma = np.where(need | popm | lowk, res, 0.)
    out = ESRSpectrum(e, Vm, Dm, K, sigma=sigma, route="KRYLOV",
                      reached=theta[-1] - theta[0])
    out.sigma_measured = res
    return out
