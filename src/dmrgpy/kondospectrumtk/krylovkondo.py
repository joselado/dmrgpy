"""The thermal third-order Kondo spectrum from Krylov spaces of the low-lying
states, with no full diagonalization (get_kondo_spectrum(mode="ED",
submode="KRYLOV")).

The sums of conductance.py run over initial states i, final states f and
intermediate states m. The initial states carry the Boltzmann weight, so
only the lowest ones matter: they are found with a sparse eigensolver, and
k is raised until the relative weight of the last one is below p_cut. The
states f and m enter through the matrix elements <f|S_a|i> and
<m|S_c|i> of the spin of the tip site, so for each initial state they lie
in the smallest H-invariant space that holds S_x|i>, S_y|i> and S_z|i>.
A block Krylov space started from the three vectors together is that
space once it is exhausted, and the Ritz states of H in it are then exact
eigenstates. The middle vertex <f|S_b|m> is taken between two such
states, which is exact for them.

One block, not three Lanczos spaces. When H conserves Sz, S+|i>, S-|i>
and Sz|i> lie in different sectors and three separate Krylov spaces are
orthogonal. With a field that breaks Sz conservation they overlap, and
three separate sets of Ritz states count the shared directions two or
three times; this gave errors of order one in the project this route comes
from. The block space is orthonormalized as a whole (two passes of
reorthogonalization against every earlier block, then an SVD that drops
the dependent directions), so every direction is counted once.

The initial state is projected out of the starting block, and the space
stays orthogonal to it. So [i] + Ritz states is an orthonormal basis with
i at its own energy, and for each initial state the code builds an object
with the six fields that conductance.py reads (e, p one-hot on i, Sx, Sy,
Sz in that basis, T, kB, dim), the way taskkondo.TaskKondoSpectrum does,
and calls the second-order, third-order Kondo and potential-interference
sums unchanged. The total is their sum with the Boltzmann weights, since
every term is linear in the populations. The elastic transition and the
transitions to states degenerate with i are in the basis: an eigenstate
other than i is reached through the Ritz states.

The route comes from a project on 1T-TaSe2 (2026-11), where it was
checked against full ED; the checks here are in
tests/test_kondo_krylov.py and docs/new_methods/kondo_krylov_thermal/
README.md. It is third order only: the leading-log flow sums over
products of spin operators between any two states and needs the full
eigenbasis.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from . import conductance

# The relative Boltzmann weight below which an initial state is dropped.
# Every term is linear in the weights, so the dropped states change dI/dV
# by about their total weight.
P_CUT = 1e-8
# Up to this dimension the initial states come from a dense eigh of the
# whole Hamiltonian, which is exact and complete; above it from eigsh.
DENSE_DIM = 1024
# When eigsh would need more than a third of the states, the dense eigh is
# done up to this dimension, and the call raises above it.
DENSE_MAX = 16384
# eigsh starts with this many states and raises k by 1.6 until the last one
# is below p_cut (or above the ground manifold at T=0).
K0 = 32
# A direction of a Krylov block is kept when its singular value, after the
# projection on the earlier blocks, is above this fraction of the largest
# column norm of the block before the projection.
DROP_TOL = 1e-9
# Matrices with an imaginary part below this fraction of their largest
# entry are taken as real.
REAL_TOL = 1e-13


class StateSpectrum():
    """One initial state (row 0) and the Ritz states of its Krylov space,
    in the shape of edkondo.KondoSpectrum: e (energies relative to the
    initial state), p (one-hot on row 0), Sx/Sy/Sz [a,b] = <a|S|b>, T, kB,
    dim. residual is the norm ||(H - e_a)|a>|| of each Ritz state (0 for
    row 0 up to the accuracy of the eigensolver)."""
    def __init__(self, e, S, T, kB, residual):
        self.e = np.asarray(e, dtype=float)
        self.dim = len(self.e)
        self.p = np.zeros(self.dim)
        self.p[0] = 1.
        self.Sx, self.Sy, self.Sz = S
        self.T = float(T)
        self.kB = kB
        self.residual = residual


def _is_real(M):
    if sp.issparse(M):
        data = np.asarray(M.tocoo().data)
    else:
        data = np.asarray(M).ravel()
    if data.size == 0 or not np.iscomplexobj(data):
        return True
    return bool(np.max(np.abs(data.imag)) <= REAL_TOL*max(np.max(np.abs(data)), 1e-300))


def _as_real(M):
    if sp.issparse(M):
        M = M.tocsr()
        return sp.csr_matrix((np.real(M.data), M.indices, M.indptr), shape=M.shape)
    return np.real(M)


def onenorm(H):
    """The 1-norm of H (largest column sum of |H|), a bound of its
    spectral radius, so 2*onenorm(H) bounds the width of the spectrum."""
    if sp.issparse(H):
        return float(abs(H).sum(axis=0).max())
    return float(np.max(np.sum(np.abs(H), axis=0)))


def low_states(H, kT, p_cut=P_CUT, tol_degenerate=0.):
    """(E, V): the lowest eigenpairs of H, sorted, enough of them that the
    last one has a relative Boltzmann weight exp(-(E-E0)/kT) below p_cut at
    kT>0, or lies above E0 + tol_degenerate at kT=0. Dense eigh up to
    DENSE_DIM, eigsh with k raised by 1.6 otherwise; when eigsh would need
    more than a third of the states, dense eigh up to DENSE_MAX, and a
    ValueError above it. Every state below the last one is
    returned, so a degenerate level is never cut in two, as long as eigsh
    finds its multiplicity."""
    dim = H.shape[0]

    def done(w):
        if kT > 0.:
            return (w[-1] - w[0])/kT > np.log(1./p_cut)
        return w[-1] - w[0] > tol_degenerate

    if dim > DENSE_DIM:
        k = min(K0, dim - 2)
        while k <= dim//3:
            w, v = sla.eigsh(H, k=k, which="SA", tol=0.)
            o = np.argsort(w)
            w, v = w[o], v[:, o]
            if done(w):
                return w, v
            k = int(1.6*k) + 1
        if dim > DENSE_MAX:
            raise ValueError(
                "get_kondo_spectrum(submode=\"KRYLOV\"): more than a third of "
                "the %d states carry a Boltzmann weight above p_cut=%g at this "
                "temperature, and a dense diagonalization of that size is not "
                "done; raise p_cut or lower T" % (dim, p_cut))
    Hd = H.toarray() if sp.issparse(H) else np.asarray(H)
    return np.linalg.eigh(Hd)


def initial_states(H, kT, p_cut=P_CUT, degeneracy_tol=1e-9):
    """(E, V, w): the initial states and their weights, normalized to 1. At
    kT>0 the states with exp(-(E-E0)/kT) > p_cut, with Boltzmann weights;
    at kT=0 the states within degeneracy_tol*2*onenorm(H) of the lowest
    one, with equal weights (the T->0+ limit, as edkondo.KondoSpectrum,
    whose tolerance is relative to the full width of the spectrum, which
    2*onenorm(H) bounds)."""
    tol = degeneracy_tol*2.*onenorm(H)
    E, V = low_states(H, kT, p_cut, tol)
    x = E - E[0]
    if kT > 0.:
        keep = np.exp(-x/kT) > p_cut
        w = np.exp(-x[keep]/kT)
    else:
        keep = x <= tol
        w = np.ones(int(np.sum(keep)))
    return E[keep], V[:, keep], w/np.sum(w)


def block_krylov(H, starts, psi, nkry):
    """(theta, F, residual): the Ritz values and vectors of H in the block
    Krylov space started from the columns of starts, orthogonal to psi, of
    dimension at most 3*nkry. Each new block H@B is projected twice on psi
    and on every earlier block, and an SVD keeps its directions with a
    singular value above DROP_TOL of its largest column norm; the space
    stops growing when a block is empty (the space is exhausted) or would
    pass 3*nkry columns."""
    blocks = []

    def orth(V):
        scale = float(np.max(np.linalg.norm(V, axis=0))) if V.size else 0.
        if scale == 0.:
            return V[:, :0]
        for _ in range(2):
            V = V - np.outer(psi, psi.conj() @ V)
            for Q in blocks:
                V = V - Q @ (Q.conj().T @ V)
        U, s, _ = np.linalg.svd(V, full_matrices=False)
        return U[:, s > DROP_TOL*scale]

    B = orth(starts)
    n = 0
    while B.shape[1] > 0 and n + B.shape[1] <= 3*nkry:
        blocks.append(B)
        n += B.shape[1]
        B = orth(H @ B)
    if not blocks:
        return np.zeros(0), np.zeros((len(psi), 0), dtype=psi.dtype), np.zeros(0)
    Q = np.concatenate(blocks, axis=1)
    HQ = H @ Q
    Hk = Q.conj().T @ HQ
    theta, y = np.linalg.eigh((Hk + Hk.conj().T)/2.)
    F = Q @ y
    residual = np.linalg.norm(HQ @ y - F*theta[None, :], axis=0)
    return theta, F, residual


def _total_spin(chain, edobj):
    """The sparse matrices of the total spin sum_j S^a_j, a = x, y, z."""
    out = []
    for name in ("Sx", "Sy", "Sz"):
        ops = getattr(chain, name)
        tot = ops[0]
        for j in range(1, len(ops)):
            tot = tot + ops[j]
        M = edobj.get_operator(tot)
        out.append(sp.csr_matrix(M) if not sp.issparse(M) else M.tocsr())
    return out


# The commutator [H, S^a_tot] counts as zero below this fraction of
# onenorm(H)*onenorm(S^a_tot), on a random state.
SU2_TOL = 1e-10


def su2_initial_states(chain, H, kT, p_cut=P_CUT, degeneracy_tol=1e-9):
    """(E, V, w, S): one member of each spin multiplet among the initial
    states, in the sector of the smallest |S^z_tot| (0, or 1/2 for a
    half-integer total spin, where every multiplet has exactly one member),
    with the weight (2S+1) exp(-(E-E0)/kT), normalized. With an
    unpolarized tip each term of dI/dV is the expectation value of a
    rotation-invariant operator, the same for every member of a multiplet,
    so this is the sum over all of them. Raises ValueError unless H
    commutes with the total spin (checked on a random state) and S^z_tot is
    diagonal in the ED basis. Inside a degenerate level the states are
    rotated to a definite S with S^2_tot, and the call raises when 2S is
    not an integer there."""
    edobj = chain.get_ED_obj()
    Stot = _total_spin(chain, edobj)
    dim = H.shape[0]
    rng = np.random.default_rng(1)
    v = rng.standard_normal(dim) + 1j*rng.standard_normal(dim)
    v /= np.linalg.norm(v)
    nH = onenorm(H)
    for a, M in zip("xyz", Stot):
        c = np.linalg.norm(H @ (M @ v) - M @ (H @ v))
        if c > SU2_TOL*nH*max(onenorm(M), 1e-300):
            raise ValueError(
                "get_kondo_spectrum(submode=\"KRYLOV\", su2=True): H does not "
                "commute with the total spin (|[H, S%s_tot] v| = %.2g on a "
                "random state), so the members of a multiplet are not "
                "equivalent. Use su2=False." % (a, c))
    Sz = Stot[2]
    offdiag = Sz - sp.diags(Sz.diagonal())
    if offdiag.nnz and np.max(np.abs(offdiag.data)) > 1e-12:
        raise ValueError("get_kondo_spectrum(submode=\"KRYLOV\", su2=True): "
                         "S^z_tot is not diagonal in the ED basis of this chain")
    m = np.real(Sz.diagonal())
    mmin = np.min(np.abs(m))
    idx = np.flatnonzero(np.abs(m - mmin) < 1e-9)
    Hs = H[idx][:, idx]
    S2 = (Stot[0] @ Stot[0] + Stot[1] @ Stot[1] + Stot[2] @ Stot[2]).tocsr()[idx][:, idx]
    tol = degeneracy_tol*2.*nH
    # the (2S+1) of a multiplet is at most 2 max|m| + 1, so the cut on
    # exp(-x/kT) is lowered by that much to keep every state above p_cut
    pc = p_cut/(2.*np.max(np.abs(m)) + 1.)
    E, V = low_states(Hs, kT, pc, tol)
    x = E - E[0]
    keep = np.exp(-x/kT) > pc if kT > 0. else x <= tol
    E, V, x = E[keep], V[:, keep], x[keep]
    spins = np.empty(len(E))
    a = 0
    while a < len(E):
        b = a + 1
        while b < len(E) and E[b] - E[a] <= tol:
            b += 1
        Vc = V[:, a:b]
        B = Vc.conj().T @ (S2 @ Vc)
        if not np.iscomplexobj(V):
            B = np.real(B) # S^2_tot is real in a real basis
        s2, U = np.linalg.eigh(B)
        S = (-1. + np.sqrt(1. + 4.*np.clip(np.real(s2), 0., None)))/2.
        if np.any(np.abs(2*S - np.round(2*S)) > 1e-6):
            raise ValueError(
                "get_kondo_spectrum(submode=\"KRYLOV\", su2=True): a level "
                "at %.6g has S^2_tot eigenvalues %s, not of a definite spin"
                % (E[a], np.array2string(np.real(s2), precision=6)))
        V[:, a:b] = Vc @ U
        spins[a:b] = np.round(2*S)/2.
        a = b
    w = (2*spins + 1.)*(np.exp(-x/kT) if kT > 0. else 1.)
    keep = w > p_cut*np.max(w)
    full = np.zeros((dim, int(np.sum(keep))), dtype=V.dtype)
    full[idx] = V[:, keep]
    return E[keep], full, w[keep]/np.sum(w[keep]), spins[keep]


class KrylovKondoSpectrum():
    """The thermal initial states of a chain's sparse ED Hamiltonian, their
    weights, and for each one a StateSpectrum from its Krylov space (see the
    module docstring). spectra() yields (weight, StateSpectrum) one at a
    time, so only one Ritz basis is held at a time.

    chain, site, T, kB: as edkondo.KondoSpectrum. nkry: the block Krylov
    space of each initial state has at most 3*nkry states. p_cut: relative
    Boltzmann weight below which an initial state is dropped. su2: take one
    member per spin multiplet (see SU2 below); H must then commute with the
    total spin, which is checked."""
    def __init__(self, chain, site, T, kB=8.617333262e-5, nkry=100, p_cut=P_CUT,
                 su2=False, degeneracy_tol=1e-9):
        if T < 0.:
            raise ValueError("T must be >= 0, got %r" % (T,))
        if isinstance(nkry, bool) or not isinstance(nkry, (int, np.integer)) or nkry < 1:
            raise ValueError("nkry must be a positive integer, got %r" % (nkry,))
        if not 0. < p_cut < 1.:
            raise ValueError("p_cut must be between 0 and 1, got %r" % (p_cut,))
        self.chain, self.site, self.T, self.kB = chain, site, float(T), kB
        self.nkry = int(nkry)
        self.degeneracy_tol = degeneracy_tol
        edobj = chain.get_ED_obj()
        H = edobj.get_hamiltonian()
        H = H.tocsr() if sp.issparse(H) else np.asarray(H)
        S = []
        for name in ("Sx", "Sy", "Sz"):
            M = edobj.get_operator(getattr(chain, name)[site])
            S.append(M.tocsr() if sp.issparse(M) else np.asarray(M))
        self.dim = H.shape[0]
        # S+ and S- are real in a real basis (Sy is imaginary there), and so
        # is everything built from a real H and real start vectors
        Sp, Sm = S[0] + 1j*S[1], S[0] - 1j*S[1]
        self.real = all(_is_real(M) for M in (H, Sp, Sm, S[2]))
        if self.real:
            H, Sp, Sm, Sz = (_as_real(M) for M in (H, Sp, Sm, S[2]))
            self._ladder = (Sp, Sm, Sz)
        else:
            H = H.astype(complex)
            self._ladder = None
        self._S = S
        self.H = H
        kT = kB*self.T
        self.tol_degenerate = degeneracy_tol*2.*onenorm(H)
        if su2:
            self.E, self.V, self.weights, self.spins = su2_initial_states(
                chain, H, kT, p_cut, degeneracy_tol)
        else:
            self.E, self.V, self.weights = initial_states(H, kT, p_cut, degeneracy_tol)
            self.spins = None
        if self.real:
            self.V = np.real(self.V)

    def state_spectrum(self, n):
        """The StateSpectrum of the initial state n."""
        psi = self.V[:, n]
        psi = psi/np.linalg.norm(psi)
        Ei = float(np.real(psi.conj() @ (self.H @ psi)))
        if self.real:
            Sp, Sm, Sz = self._ladder
            starts = np.stack([Sp @ psi, Sm @ psi, Sz @ psi], axis=1)
        else:
            starts = np.stack([M @ psi for M in self._S], axis=1)
        theta, F, residual = block_krylov(self.H, starts, psi, self.nkry)
        basis = np.concatenate([psi[:, None], F], axis=1)
        if self.real:
            Ap, Am, Az = (basis.T @ (M @ basis) for M in self._ladder)
            mats = ((Ap + Am)/2., (Ap - Am)/2j, Az.astype(complex))
        else:
            mats = tuple(basis.conj().T @ (M @ basis) for M in self._S)
        e = np.concatenate([[0.], theta - Ei])
        # a state degenerate with i comes out at +-1e-16, and Theta0 would
        # read the zero-bias point on one side at random: snapped to 0
        # within the degeneracy tolerance, as taskkondo.py does
        e[np.abs(e) <= self.tol_degenerate] = 0.
        res0 = float(np.linalg.norm(self.H @ psi - Ei*psi))
        return StateSpectrum(e, mats, self.T, self.kB, np.concatenate([[res0], residual]))

    def spectra(self):
        for n, w in enumerate(self.weights):
            yield w, self.state_spectrum(n)


def krylov_kondo_dIdV(kks, eV, Jrho_s, U=0.0, T0=1.0, omega0=20e-3, Gamma0=5e-6,
                      order=3):
    """dI/dV of conductance.py (second order, plus the third-order Kondo
    and potential-interference terms at order=3) summed over the initial
    states of a KrylovKondoSpectrum with their weights."""
    from .stepfunctions import cached_fbuilder
    eV = np.asarray(eV, dtype=float)
    Fb = None
    if order == 3 and kks.T > 0.:
        Fb = cached_fbuilder(kks.T, omega0=omega0, Gamma0=Gamma0, kB=kks.kB)
    total = np.zeros(len(eV))
    for w, ks in kks.spectra():
        d = conductance.second_order_dIdV(ks, eV, T0=T0, U=U)
        if order == 3:
            d = d + conductance.third_order_kondo_dIdV(
                ks, eV, Jrho_s, T0=T0, omega0=omega0, Gamma0=Gamma0, Fb=Fb)
            if U != 0.0:
                d = d + conductance.third_order_potential_dIdV(
                    ks, eV, Jrho_s, U, T0=T0, omega0=omega0, Gamma0=Gamma0, Fb=Fb)
        total += w*d
    return total
