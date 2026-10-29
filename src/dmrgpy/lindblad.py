"""Open quantum systems: the Lindbladian of a chain, its spectrum and its
steady state (Phase 2 of docs/clik_mps_lindblad_plan.md).

The master equation

    d rho/dt = -i[H,rho] + sum_mu ( L_mu rho L_mu^dag
                                    - 1/2 {L_mu^dag L_mu, rho} )

is linear in rho. With |rho>> = sum rho_ab |a>|b> it becomes
d|rho>>/dt = L|rho>> with the generator

    L = -i (H x 1 - 1 x H^T)
        + sum_mu ( L_mu x L_mu^* - 1/2 L_mu^dag L_mu x 1
                   - 1/2 1 x (L_mu^dag L_mu)^T ),

an operator on a doubled chain whose sites alternate between the ket
(physical, 2i) and the bra (auxiliary, 2i+1) index of site i. That keeps a
nearest-neighbor H next-nearest-neighbor in L. Its eigenvalues have
Re(lambda) <= 0; the steady state is the eigenvector at lambda = 0 and the
dissipative gap is -Re(lambda_2).

Spins and bosons only: on a doubled chain the Jordan-Wigner strings of
fermions need their own treatment, which is not implemented.
"""

import re

import numpy as np

from . import multioperator


# Name of the transpose of each single-site operator name, with the sign
# it carries: Sy and ISy are antisymmetric matrices in the S^z basis, the
# ladder operators transpose into each other, and every name listed as
# symmetric is a real symmetric matrix on every spin and boson site.
_TRANSPOSE = {"Sy": (-1.0, "Sy"), "ISy": (-1.0, "ISy"),
              "Sp": (1.0, "Sm"), "Sm": (1.0, "Sp"),
              "S+": (1.0, "S-"), "S-": (1.0, "S+"),
              "A": (1.0, "Adag"), "Adag": (1.0, "A")}
_SYMMETRIC = {"Id", "Sx", "Sz", "Sx2", "Sy2", "Sz2", "S2", "N", "Nb",
              "projUp", "projDn", "projUpi", "projDni", "projZ0"}
_PROJECTOR = re.compile(r"^N\d+$")


def _transpose_name(name):
    if name in _TRANSPOSE:
        return _TRANSPOSE[name]
    if name in _SYMMETRIC or _PROJECTOR.match(name):
        return 1.0, name
    raise ValueError("Lindblad_Chain: no transpose known for the operator "
                     "name " + repr(name) + "; spin and boson operators "
                     "only")


def _rebuild(terms):
    out = multioperator.MultiOperator(term=False)
    out.op = terms
    out.i = len(terms) - 1
    return multioperator.obj2MO(out)


def transpose(mo):
    """O^T for a MultiOperator of spin or boson names: coefficients kept,
    each factor transposed, the factor order reversed."""
    terms = []
    for term in mo.op:
        c = term[0]
        factors = []
        for name, i in reversed(term[1:]):
            s, tname = _transpose_name(name)
            c = c * s
            factors.append([tname, i])
        terms.append([c] + factors)
    return _rebuild(terms)


def conjugate(mo):
    """O^* = (O^dagger)^T."""
    return transpose(mo.get_dagger())


def on_sublattice(mo, offset):
    """mo with site i moved to 2i+offset: 0 the physical, 1 the auxiliary
    sublattice of the doubled chain."""
    terms = []
    for term in mo.op:
        terms.append([term[0]] + [[name, 2 * i + offset]
                                  for name, i in term[1:]])
    return _rebuild(terms)


def vectorized_lindbladian(H, jumps):
    """The generator L of the module docstring as a MultiOperator on the
    doubled chain; H may be None or 0 for pure dissipation."""
    L = 0
    if H is not None and not (np.isscalar(H) and H == 0):
        L = -1j * (on_sublattice(H, 0) - on_sublattice(transpose(H), 1))
    for J in jumps:
        JdJ = J.get_dagger() * J
        L = L + (on_sublattice(J, 0) * on_sublattice(conjugate(J), 1)
                 - 0.5 * on_sublattice(JdJ, 0)
                 - 0.5 * on_sublattice(transpose(JdJ), 1))
    return L


def _check_sites(sites):
    for s in sites:
        if not (2 <= s <= 6 or 100 < s < 200):
            raise NotImplementedError(
                "Lindblad_Chain supports spin and boson sites only (site "
                "code " + repr(s) + "); fermions on a doubled chain need "
                "their own Jordan-Wigner treatment")


def _doubled_chain(chain):
    """A chain of the same model with every site doubled, on "python",
    so that it also has the model's ED route."""
    from .bosonchain import Bosonic_Chain
    from .spinchain import Spin_Chain
    if isinstance(chain, Bosonic_Chain):
        maxnb = [m for m in chain.maxnb for _ in (0, 1)]
        return Bosonic_Chain(len(maxnb), maxnb=maxnb, itensor_version="python")
    if all(2 <= s <= 6 for s in chain.sites):
        return Spin_Chain([s for s in chain.sites for _ in (0, 1)],
                          itensor_version="python")
    raise NotImplementedError("Lindblad_Chain: build the chain as a "
                              "Spin_Chain or a Bosonic_Chain")


class Lindblad_Chain():
    """An open chain: H and the jump operators of a chain, its Lindbladian
    on the doubled chain, and the spectrum and steady state of it."""

    def __init__(self, chain, jumps, H=None):
        _check_sites(chain.sites)
        if H is None:
            H = chain.hamiltonian
        self.chain = chain          # the physical chain
        self.H = H
        self.jumps = list(jumps)
        if not self.jumps:
            raise ValueError("Lindblad_Chain: give at least one jump "
                             "operator; without them use the chain itself")
        self.L = vectorized_lindbladian(H, self.jumps)
        self.doubled = _doubled_chain(chain)
        self.doubled.maxm = chain.maxm
        self.doubled.cutoff = chain.cutoff
        self._spectrum = None

    def get_spectrum(self, nk=20, dt=2.0, alpha=None, gram_tol=1e-12,
                     seed=0, start=None, mode="DMRG"):
        """(lam, res, obs): eigenvalues of L, largest real part first (the
        steady state, lambda ~ 0, comes first), the residual of each
        normalized eigenvector, and obs(O, k) = Tr(O r_k)/Tr(r_k) for an
        operator O of the physical chain. Only the steady state has a
        nonzero trace, so for k > 0 obs gives Tr(O r_k) with ||r_k|| = 1.
        DMRG is get_krylov_spectrum on the doubled chain
        (itensor_version="python"); mode="ED" diagonalizes L densely.

        The real branch reaches the steady state only once nk*dt is several
        times the inverse dissipative gap: on a driven 4-site XXZ chain with
        a gap of 0.14, the steady state was 1.4e-2 off at nk*dt=10 and
        9e-7 off at nk*dt=40, hence dt=2.0 here against 0.5 for a
        Hamiltonian. A conserved quantity of L (a strong symmetry) gives one
        steady state per sector, and random starts then find a degenerate
        set at lambda ~ 0. `start=(psi, phi)`, two lists with one local
        basis state per physical site (the occupation, for bosons), seeds
        the two branches with |psi><psi| and the traceless |psi><phi|
        instead of random states, as in arXiv:2509.07709. Both lie in the
        sector of psi when phi has the same charge, and the evolution
        keeps them there, so only that sector's steady state and modes
        appear."""
        if mode == "ED":
            out = self._spectrum_ed()
        elif mode == "DMRG":
            out = self._spectrum_dmrg(nk, dt, alpha, gram_tol, seed, start)
        else:
            raise ValueError("Lindblad_Chain: unrecognized mode " + repr(mode))
        self._spectrum = out
        return out

    def steady_state_vev(self, O, **kwargs):
        """Tr(O rho_ss) in the normalized steady state, the eigenvector of
        largest real part; computes the spectrum if it was not yet, with
        get_spectrum's keywords."""
        if self._spectrum is None or kwargs:
            self.get_spectrum(**kwargs)
        return self._spectrum[2](O, 0)

    def _spectrum_dmrg(self, nk, dt, alpha, gram_tol, seed, start=None):
        from .krylovspectrum import krylov_spectrum_dmrg
        from .pyitensor.mpsalgebra import inner
        from .pyitensor.vectorize import identity_mps
        dc = self.doubled
        starts = None if start is None else self._seeds(*start)
        ks = krylov_spectrum_dmrg(dc, self.L, nk=nk, dt=dt, alpha=alpha,
                                  gram_tol=gram_tol, seed=seed, starts=starts)
        one = identity_mps(dc._session.sites)
        traces = np.array([inner(one, rho) for rho in ks.states])

        def obs(O, k):
            Om = dc._session._mpo(on_sublattice(O, 0).to_terms())
            t = np.array([inner(one, Om, rho) for rho in ks.states])
            c = ks.coefficients[:, k]
            value = c @ t
            tr = c @ traces
            if k == 0:
                return value / tr
            return value / np.sqrt(abs(np.vdot(c, ks.M @ c)))
        return ks.eigenvalues, ks.residuals, obs

    def _seeds(self, psi, phi):
        """Local vectors of |psi>|psi> and |psi>|phi> on the doubled chain."""
        n = len(self.chain.sites)
        if len(psi) != n or len(phi) != n:
            raise ValueError("Lindblad_Chain: start=(psi, phi) needs one "
                             "basis state per site, %d of them" % n)
        sites = self.doubled._session.sites

        def basis(k, j):
            v = np.zeros(sites.dim(k), dtype=complex)
            v[j] = 1.0
            return v
        p = [basis(2 * i + o, psi[i]) for i in range(n) for o in (1, 2)]
        u = [basis(2 * i + 1, psi[i]) if o == 1 else basis(2 * i + 2, phi[i])
             for i in range(n) for o in (1, 2)]
        return p, u

    def _spectrum_ed(self):
        """Dense L on the physical Hilbert space, vec(rho)_{ab} = rho_ab
        row-major, so that A rho B is kron(A, B^T) vec(rho)."""
        ed = self.chain.get_ED_obj()

        def dense(op):
            m = ed.get_operator(op)
            return m.toarray() if hasattr(m, "toarray") else np.asarray(m)
        Js = [dense(J) for J in self.jumps]
        D = Js[0].shape[0]
        I = np.eye(D)
        Ld = np.zeros((D * D, D * D), dtype=complex)
        if self.H is not None and not (np.isscalar(self.H) and self.H == 0):
            Hd = dense(self.H)
            Ld += -1j * (np.kron(Hd, I) - np.kron(I, Hd.T))
        for Jd in Js:
            JdJ = Jd.conj().T @ Jd
            Ld += (np.kron(Jd, Jd.conj()) - 0.5 * np.kron(JdJ, I)
                   - 0.5 * np.kron(I, JdJ.T))
        lam, V = np.linalg.eig(Ld)
        order = np.argsort(-lam.real, kind="stable")
        lam, V = lam[order], V[:, order]
        V = V / np.linalg.norm(V, axis=0)
        res = np.linalg.norm(Ld @ V - V * lam, axis=0)

        def obs(O, k):
            R = V[:, k].reshape(D, D)
            value = np.trace(dense(O) @ R)
            return value / np.trace(R) if k == 0 else value
        return lam, res, obs
