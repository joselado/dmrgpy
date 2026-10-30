"""Slowest eigenvalues of a generator from a Krylov space of time-evolved
states (the eigensolver of arXiv:2509.07709, CLIK-MPS).

The generator G defaults to -H, so the eigenvalues come back as -E, lowest
energy first; for a non-Hermitian H they are complex. The algorithm, its
parameters and the floor of its residuals are in pyitensor/clik.py; on itensor_version=3 the same
algorithm runs through cliksession.py.
"""

import numpy as np


def get_krylov_spectrum(self, G=None, nk=20, dt=0.5, alpha=None,
                        gram_tol=1e-12, seed=0, exact_residuals=False,
                        mode="DMRG"):
    """Return (ev, res, obs) for the generator G (default -H).

    - ev: eigenvalues of G on the Krylov space, largest real part first
      (for the default G=-H, the energies are -ev)
    - res: how far each eigenvalue moves when the last snapshot of each
      branch is dropped, the paper's tracking of the spectrum, at no extra
      cost; a converged eigenvalue does not move. It ranks the modes by
      convergence but underestimates the error, by a factor 3 to 40 on
      small chains against ED. With
      exact_residuals=True it is ||G r - lambda r|| for each normalized
      eigenvector r instead, which costs a matrix element of G^dag G for
      every pair of snapshots and cannot resolve values below about 1e-8
      of the operator scale. mode="ED" always gives the exact residual.
    - obs(O, k): <r_k|O|r_k>/<r_k|r_k> for a MultiOperator O

    nk snapshots per branch, imaginary-time step dt, contour angle alpha of
    the second branch (None: estimated), gram_tol the relative deflation
    threshold, seed for the random starts. DMRG runs on
    itensor_version="python" and 3; mode="ED" diagonalizes G densely and
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
    ks = krylov_spectrum_dmrg(self, G, nk=nk, dt=dt, alpha=alpha,
                              gram_tol=gram_tol, seed=seed,
                              exact_residuals=exact_residuals)

    def obs(O, k):
        return ks.expectation(operator_handle(self, O), k)
    return ks.eigenvalues, _res(ks), obs


def operator_handle(self, O):
    """The session's MPO of the MultiOperator O, for ks.expectation."""
    if self.itensor_version == "python":
        return self._session._mpo(O.to_terms())
    return self.toMPO(O).cpp_handle


def _res(ks):
    """The exact residuals when they were computed, else the tracking."""
    return ks.residuals if ks.residuals is not None else ks.changes


def krylov_spectrum_dmrg(self, G, nk=20, dt=0.5, alpha=None, gram_tol=1e-12,
                         seed=0, starts=None, exact_residuals=False,
                         traceless=False, sector_charge=None, substeps=1):
    """The clik.KrylovSpectrum of the MultiOperator G on this chain, with
    its stored states, for a caller that needs more than <r|O|r>
    (Lindblad_Chain reads traces off it). `starts`, `traceless`,
    `sector_charge` and `substeps`: see pyitensor Chain.krylov_spectrum
    and clik.krylov_spectrum. On itensor_version=3 it is cliksession.py's
    port, whose KrylovSpectrum reads its matrix elements on the session."""
    if self.itensor_version not in ("python", 3):
        raise NotImplementedError(
            "get_krylov_spectrum runs on itensor_version=\"python\" and 3 "
            "only (got " + repr(self.itensor_version) + "); use mode=\"ED\" "
            "for a small chain or setup_python() to switch this one")
    if self.conserved_sector:
        raise NotImplementedError(
            "get_krylov_spectrum does not support a conserved sector yet: "
            "its random starts are not confined to it")
    # decided before the sweep parameters go to the session: the witness
    # probe of is_hermitian() runs at a bond dimension of at most 8 and
    # leaves the session there, which put the whole v3 evolution at 8
    hermitian = self.is_hermitian(G) if self.itensor_version == 3 else None
    self._session.set_sweep_params(self.maxm, self.nsweeps, self.cutoff,
                                   self.noise)
    self._session.set_mpomaxm(max(self.maxm, self.mpomaxm))
    if self.itensor_version == 3:
        if not hasattr(self._session, "clik_random_starts"):
            raise NotImplementedError(
                "get_krylov_spectrum: the compiled v3 extension predates "
                "CLIK-MPS; rebuild it (make pybind in src/dmrgpy/mpscpp3)")
        from . import cliksession
        return cliksession.krylov_spectrum(
                self, G, hermitian, nk=nk, dt=dt, alpha=alpha,
                gram_tol=gram_tol, seed=seed, starts=starts,
                exact_residuals=exact_residuals, traceless=traceless,
                sector_charge=sector_charge, substeps=substeps)
    return self._session.krylov_spectrum(G.to_terms(), nk=nk, dt=dt,
                                         alpha=alpha, gram_tol=gram_tol,
                                         seed=seed, starts=starts,
                                         exact_residuals=exact_residuals,
                                         traceless=traceless,
                                         sector_charge=sector_charge,
                                         substeps=substeps)


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
