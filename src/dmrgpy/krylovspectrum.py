"""Slowest eigenvalues of a generator from a Krylov space of time-evolved
states (the eigensolver of arXiv:2509.07709, CLIK-MPS).

The generator G defaults to -H, so the eigenvalues come back as -E, lowest
energy first; for a non-Hermitian H they are complex. The algorithm, its
parameters and the floor of its residuals are in pyitensor/clik.py.
"""

import numpy as np


def get_krylov_spectrum(self, G=None, nk=20, dt=0.5, alpha=None,
                        gram_tol=1e-12, seed=0, mode="DMRG"):
    """Return (ev, res, obs) for the generator G (default -H).

    - ev: eigenvalues of G on the Krylov space, largest real part first
      (for the default G=-H, the energies are -ev)
    - res: ||G r - lambda r|| for each normalized eigenvector r; an
      eigenvalue is converged only when its residual is small, and the
      DMRG residual cannot resolve values below about 1e-8 of the
      operator scale
    - obs(O, k): <r_k|O|r_k>/<r_k|r_k> for a MultiOperator O

    nk snapshots per branch, imaginary-time step dt, contour angle alpha of
    the second branch (None: estimated), gram_tol the relative deflation
    threshold, seed for the random starts. DMRG runs on
    itensor_version="python" only; mode="ED" diagonalizes G densely and
    returns every eigenvalue, in the same order."""
    from . import groundstate
    if G is None:
        groundstate.require_hamiltonian(self)
        G = -1 * self.hamiltonian
    mode = self.get_mode(mode=mode)
    if mode == "ED":
        return _krylov_spectrum_ed(self, G)
    if mode != "DMRG":
        raise ValueError("get_krylov_spectrum: unrecognized mode " + repr(mode))
    if self.itensor_version != "python":
        raise NotImplementedError(
            "get_krylov_spectrum runs on itensor_version=\"python\" only (got "
            + repr(self.itensor_version) + "); use mode=\"ED\" for a small "
            "chain or setup_python() to switch this one")
    if self.conserved_sector:
        raise NotImplementedError(
            "get_krylov_spectrum does not support a conserved sector yet: "
            "its random starts are not confined to it")
    self._session.set_sweep_params(self.maxm, self.nsweeps, self.cutoff,
                                   self.noise)
    self._session.set_mpomaxm(max(self.maxm, self.mpomaxm))
    ks = self._session.krylov_spectrum(G.to_terms(), nk=nk, dt=dt,
                                       alpha=alpha, gram_tol=gram_tol,
                                       seed=seed)

    def obs(O, k):
        return ks.expectation(self._session._mpo(O.to_terms()), k)
    return ks.eigenvalues, ks.residuals, obs


def _krylov_spectrum_ed(self, G):
    """The ED reference: every eigenvalue of the dense G, in the DMRG
    order, with its exact residual."""
    ed = self.get_ED_obj()
    Gd = ed.get_operator(G)
    Gd = Gd.toarray() if hasattr(Gd, "toarray") else np.asarray(Gd)
    lam, V = np.linalg.eig(Gd)
    order = np.argsort(-lam.real, kind="stable")
    lam, V = lam[order], V[:, order]
    V = V / np.linalg.norm(V, axis=0)
    res = np.linalg.norm(Gd @ V - V * lam, axis=0)

    def obs(O, k):
        Od = ed.get_operator(O)
        v = V[:, k]
        return np.vdot(v, Od @ v) / np.vdot(v, v)
    return lam, res, obs
