import weakref
import numpy as np
from scipy.sparse import issparse, csr_matrix
from ..algebra import algebra

# Full-spectrum ED building blocks for the third-order Kondo/STM
# perturbation theory (Ternes, arXiv:1505.04430). Unlike the ground-state
# ED already used elsewhere in dmrgpy (EDchain.get_gs_array()), the third
# order Kondo term sums over *every* eigenstate as a virtual intermediate
# state, so this always needs the full spectrum, not just the low-energy
# states DMRG or a Lanczos-type ED would normally target.
#
# A Hamiltonian whose matrix is real (no field along y, no complex
# hopping) is diagonalized in real arithmetic, and the spin operators are
# transformed in real arithmetic too: Sx and Sz come out real and Sy as i
# times a real matrix, which conductance.py reads off the arrays to build its
# Levi-Civita coefficients in real arithmetic. Inside a degenerate level the
# real eigh picks another basis than the complex one; every sum of
# conductance.py is a trace over each level, so it does not depend on it. The
# real eigenpairs are kept per ED object (_REAL_EIGH), since the complex ones
# are kept by the ED object itself and a scan at fixed H must not diagonalize
# again; the ED object's own get_diagonalized_hamiltonian is left untouched.
_REAL_EIGH = weakref.WeakKeyDictionary()


def _values(m):
    """The stored values of a sparse matrix, or the entries of a dense one.
    Read through a COO copy: m.imag and m.real of a sparse matrix share
    their data with m, and a method like count_nonzero() on them sums
    duplicates in place, which reorders the imaginary part of m itself when
    m is not in canonical form (it made the Hamiltonian of a mixed S=1 and
    S=1/2 chain non-Hermitian)."""
    if issparse(m): return np.asarray(m.tocoo().data)
    return np.asarray(m)


def _imag_is_zero(m):
    """Whether a dense or sparse matrix has no imaginary part at all (a
    stored value with a nonzero imaginary part says no, even if duplicates
    would cancel it)."""
    v = _values(m)
    return not np.iscomplexobj(v) or not np.any(v.imag)


def _real_part_is_zero(m):
    """Whether a dense or sparse matrix is purely imaginary."""
    return not np.any(np.real(_values(m)))


def eigenpairs(edobj):
    """All eigenpairs of the ED Hamiltonian: in real arithmetic when its
    matrix is real, kept per ED object and H matrix; otherwise the complex
    ones of edobj.get_diagonalized_hamiltonian(), kept by the ED object."""
    H = edobj.get_hamiltonian()
    if not _imag_is_zero(H):
        return edobj.get_diagonalized_hamiltonian()
    try: cached = _REAL_EIGH.get(edobj)
    except TypeError: cached = None # an ED object that cannot be weakly referenced
    if cached is not None and cached[0] is H:
        return cached[1], cached[2]
    Hd = np.real(np.asarray(algebra.todense(H)))
    emu, vs = algebra.eigh(np.ascontiguousarray(Hd))
    try: _REAL_EIGH[edobj] = (H, emu, vs)
    except TypeError: pass
    return emu, vs


class KondoSpectrum():
    """Holds everything derived once from a chain's full ED spectrum that
    the second- and third-order conductance terms need: eigenenergies
    (shifted so the ground state is at 0), Boltzmann occupations at
    temperature T, and the impurity spin operators Sx/Sy/Sz transformed
    into the Hamiltonian eigenbasis."""
    def __init__(self, chain, site, T, kB=8.617333262e-5, degeneracy_tol=1e-9):
        # kB in eV/K by default (so T is given in Kelvin, energies in eV,
        # matching how the rest of dmrgpy's examples quote energies); pass
        # kB=1 if you'd rather work in units where T is already an energy.
        # T=0 is a valid, exact limit (not an error): both the spin
        # system's Boltzmann occupation and the tunneling-electron thermal
        # smearing (Theta/F -> Theta0/F0) collapse to their T=0 closed
        # forms, see stepfunctions.py.
        if T < 0.: raise ValueError("KondoSpectrum requires T>=0")
        self.chain = chain
        self.site = site
        self.T = T
        self.kB = kB
        edobj = chain.get_ED_obj() # full-spectrum ED backend (pychain)
        emu, vs = eigenpairs(edobj) # all eigenpairs, real for a real H
        emu = np.array(emu, dtype=float)
        order = np.argsort(emu) # eigh output should already be sorted
        emu = emu[order]
        vs = np.array(vs)[:, order]
        self.e = emu - emu[0] # ground state energy set to 0
        self.vs = vs
        self.dim = len(self.e)
        if T == 0.:
            # Equal weights over the (numerically) degenerate ground-state
            # manifold -- the T->0+ limit of the Boltzmann distribution.
            # A one-hot p on eigh's first eigenvector is the same thing for
            # a non-degenerate ground state, and for a multiplet whose
            # members are related by a symmetry the unpolarized sums do
            # not distinguish (an SU(2) irrep, a Kramers doublet), but for
            # an accidental degeneracy -- e.g. S=1 with D*Sz^2+g*muB*B*Sz
            # at the |0>/|-1> level crossing g*muB*B=D -- it silently
            # picks whatever basis eigh returned for that eigenspace, and
            # the spectrum then depends on it. degeneracy_tol is relative
            # to the spectral bandwidth.
            width = max(self.e[-1], 1e-300)
            gs = (self.e <= degeneracy_tol*width).astype(float)
            self.p = gs/np.sum(gs)
        else:
            kT = kB*T
            p = np.exp(-self.e/kT)
            self.p = p/np.sum(p) # Boltzmann occupations, eq. "boltzmann"
        self.Sx = self._eigenbasis_spin_operator("Sx")
        self.Sy = self._eigenbasis_spin_operator("Sy")
        self.Sz = self._eigenbasis_spin_operator("Sz")
    def _eigenbasis_spin_operator(self, name):
        edobj = self.chain.get_ED_obj()
        op = getattr(self.chain, name)[self.site] # MultiOperator, this site
        # get_operator restricts to the chain's conserved sector, if any,
        # and refuses an operator that leaves it: Nf=n on a Hubbard chain is
        # exact at T=0, since S_i conserves the electron number, while an
        # Sz sector raises here rather than return a zero S+ column
        mat = edobj.get_operator(op)
        if not np.iscomplexobj(self.vs):
            # real eigenvectors: a real operator (Sx, Sz) stays real and an
            # imaginary one (Sy) is i times a real transform
            if _imag_is_zero(mat): return self._transform(mat, np.real)
            if _real_part_is_zero(mat): return 1j*self._transform(mat, np.imag)
        if issparse(mat): # O(dim^2) sparse product, one dense product left
            mat = csr_matrix(mat, dtype=complex)
            return np.conjugate(self.vs.T)@(mat@self.vs) # <a|S|b>
        mat = np.array(algebra.todense(mat), dtype=complex)
        return np.conjugate(self.vs.T)@mat@self.vs # Xi^alpha_{ab} = <a|S|b>
    def _transform(self, mat, part):
        """vs^T part(mat) vs in real arithmetic, for real eigenvectors."""
        if issparse(mat):
            c = mat.tocoo() # new arrays: mat is not touched
            m = csr_matrix((np.array(part(c.data), dtype=float), (c.row, c.col)),
                           shape=c.shape)
            return self.vs.T@(m@self.vs)
        m = np.array(part(np.asarray(algebra.todense(mat))), dtype=float)
        return self.vs.T@m@self.vs
