"""The levels of smallest real part of a non-Hermitian Hamiltonian by
imaginary time, on itensor_version=3 and "python".

The state e^{-tau H}|psi> of a non-Hermitian H is a sum of right
eigenvectors weighted by e^{-tau E_n}, so the levels of smallest real part
grow fastest and every other one dies out relative to them. The steps are
the implicit ones of imaginarytime.py (Zima, Stoudenmire, White, Parcollet
and Kaye, arXiv:2606.02930): x -> (1 + b(H - s))^-1 (1 - b(H - s)) x on a
logarithmic tau grid, each a two-site Galerkin sweep with local GMRES,
which takes a non-Hermitian H as it is. With the shift s at the smallest
real part every eigenvalue has Re(E - s) >= 0, and the Cayley factor
(1 - b z)/(1 + b z) of such a z lies in the unit disk for b > 0, so the
steps stay A-stable. The shift used here is the real part of the running
<H>, which needs no estimate of the answer. The state is normalized after
every step.

One evolved state does not settle when several levels share the smallest
real part, which is what PT symmetry does to a complex-conjugate pair and
what a set of purely imaginary single-particle energies does to a whole
multiplet: the state ends in their span, their phases rotate at different
rates, and its own Rayleigh quotient keeps moving. The readout is
therefore a Rayleigh-Ritz in the span of the last NSNAP states, with the
residual of every Ritz pair formed from the Gram matrices <i|j>,
<i|H|j> and <Hi|Hj>, so a Ritz value of the projection that is not an
eigenvalue (they appear below the spectrum, as in NH-DMRG's projected
problems) is recognized by its residual and dropped. Measured on the
PT-symmetric XX chain with i 0.3 (-1)^i Sz_i, maxm=64
(docs/audit_2026_11_15_hole_hunt.md, finding 3): every distinct level of
smallest real part comes out, 3 of 3 on 10 to 14 sites and 9 of 9 on 16,
real parts within 1e-5 of the exact ones, where NH-DMRG stalls from 14
sites. A level exactly degenerate with one in another symmetry sector comes
out once, as one Ritz vector in their common eigenspace.

levels() returns the accepted Ritz pairs. ground_pair() selects the one of
smallest real part, builds its Ritz vector as one weighted sum of the last
states (sum_mps_weighted), evolves a left state under H^dagger from its
conjugate, which already is the left eigenvector of a complex-symmetric H,
reads the left Ritz vector closest to conj(E) and maximally overlapping
with the right one, and biorthonormalizes the pair; gs_energy_imaginary_time
stores it as nhdmrg.gs_energy_nhdmrg stores NH-DMRG's, behind the same
two-sided certificate, and falls back to NH-DMRG when the certificate fails.
"""

import numpy as np

from . import imaginarytime as it
from . import mps

TAU = 400.0      # the end of the grid, in units of 1/(largest coefficient of H)
TAU0 = 0.01      # the first panel, in the same units, as evolve_imaginary_time's default
NSTEPS = 8       # steps per panel of the logarithmic grid
NSNAP = 24       # states of the end of the grid the Rayleigh-Ritz uses
GRAM_FLOOR = 1e-10  # Gram eigenvalues kept, relative to the largest
RITZ_TOL = 1e-3     # accepted Ritz pair: residual below this times c + |theta - e_id|


def available(self):
    """True when the imaginary-time route can run on this chain: a v3
    (3 sites or more) or "python" (2 or more) session with
    implicit_step_with, and no conserved sector, which the implicit steps
    refuse. Everything else stays on NH-DMRG."""
    if self.itensor_version not in (3, "python"): return False
    if getattr(self, "conserved_sector", None): return False
    if self.itensor_version == 3 and self.ns < 3: return False
    if self.ns < 2: return False
    return hasattr(getattr(self, "_session", None), "implicit_step_with")


def _require(self):
    if not available(self):
        raise NotImplementedError(
            "the imaginary-time levels need an itensor_version=3 (3 sites or "
            "more) or 'python' (2 or more) chain without a conserved sector "
            "(got itensor_version=%r, %d sites, sector %r)"
            % (self.itensor_version, self.ns,
               getattr(self, "conserved_sector", None)))


def _unit(H):
    """H's energy unit: its largest coefficient off the identity. A
    constant offset moves every level together and is taken out by the
    shift, so it must not set the time scale: with the offset's 1e6 as the
    unit, tau would end at 4e-4 and nothing would be projected out."""
    from .nhdmrg import _is_identity_term
    cs = [abs(c) for c, ops in H.to_terms() if not _is_identity_term(ops)]
    u = max(cs + [0.0])
    return u if u > 0. else 1.0


def _grid(H, tau, nsteps):
    """The logarithmic grid in units of H's largest coefficient off the
    identity (_unit), so that the route is the same for s*H with tau/s."""
    u = _unit(H)
    tau = TAU / u if tau is None else float(tau)
    return it.log_grid(tau, TAU0 / u, int(nsteps))


def _prepare(self):
    from .groundstate import send_hamiltonian
    s = self._session
    s.set_sweep_params(self.maxm, self.nsweeps, self.cutoff, self.noise)
    s.set_verbose(self.verbose)
    s.set_mpomaxm(max(self.maxm, self.mpomaxm))
    send_hamiltonian(self)
    return s


def _evolve(self, handle, Hop, steps, nsnap):
    """The normalized states at the last nsnap points of the grid, evolved
    from the MPS handle under the MPO Hop with the shift at the running
    Re <H>."""
    s = self._session
    snaps = []
    for k, dt in enumerate(steps):
        handle = s.implicit_step_with(Hop, handle, it.step_factors(dt, 2),
                                      0.0, True, 1, 1e-10, 10, 100)
        if k >= len(steps) - nsnap:
            snaps.append(mps.MPS(self, cpp_handle=handle).copy())
    return snaps


def _rayleigh_ritz(snaps, H):
    """[(theta, residual, c)] for every Ritz pair of H in the span of the
    states: c the coefficients of the unit-norm Ritz vector sum_i c_i
    snaps[i], residual ||H y - theta y|| for that vector, from the Gram
    matrices alone."""
    Hs = [H * b for b in snaps]
    S = np.array([[a.dot(b) for b in snaps] for a in snaps])
    M = np.array([[a.dot(hb) for hb in Hs] for a in snaps])
    Q = np.array([[ha.dot(hb) for hb in Hs] for ha in Hs])
    S = 0.5 * (S + S.conj().T)
    Q = 0.5 * (Q + Q.conj().T)
    w, U = np.linalg.eigh(S)
    keep = w > GRAM_FLOOR * w.max()
    X = U[:, keep] / np.sqrt(w[keep])
    th, Y = np.linalg.eig(X.conj().T @ M @ X)
    C = X @ Y
    out = []
    for j in range(len(th)):
        c, t = C[:, j], complex(th[j])
        nrm = float((c.conj() @ S @ c).real)
        r2 = (c.conj() @ Q @ c - t * (c.conj() @ M.conj().T @ c)
              - np.conj(t) * (c.conj() @ M @ c) + abs(t)**2 * nrm).real
        out.append((t, np.sqrt(max(r2, 0.0) / nrm), c / np.sqrt(nrm)))
    return out


def _scale(H):
    from .nhdmrg import _residual_scale
    return _residual_scale(H.to_terms())


def _accepted(pairs, H, tol):
    rscale, e_id = _scale(H)
    return sorted([p for p in pairs if p[1] < tol * (rscale + abs(p[0] - e_id))],
                  key=lambda p: (p[0].real, p[0].imag))


def _combine(self, snaps, c):
    handle = self._session.sum_mps_weighted([x.cpp_handle for x in snaps],
                                            [complex(z) for z in c])
    return mps.MPS(self, cpp_handle=handle).copy()


def _right_states(self, H, tau, nsteps, nsnap):
    _require(self)
    s = _prepare(self)
    wf = self.random_state()
    wf = wf * (1.0 / np.sqrt(wf.dot(wf).real))
    Hop = s.build_operator(H.to_terms())
    return _evolve(self, wf.cpp_handle, Hop, _grid(H, tau, nsteps), nsnap)


def levels(self, H=None, nsnap=NSNAP, tau=None, nsteps=NSTEPS, tol=RITZ_TOL,
           vectors=False):
    """(energies, residuals) of the levels of smallest real part of H
    (default self.hamiltonian), and their right eigenvectors with
    vectors=True: the Ritz pairs of the imaginary-time route (module
    docstring) whose residual is below tol times c + |E - e_id|, the
    scale of nhdmrg()'s certificate, sorted by real part. The residuals
    are relative in the same way. From one random start, so a level
    that start does not reach (a symmetry sector it has no weight in) is
    missing."""
    if H is None: H = self.hamiltonian
    if H is None: raise ValueError("set a Hamiltonian first")
    snaps = _right_states(self, H, tau, nsteps, int(nsnap))
    acc = _accepted(_rayleigh_ritz(snaps, H), H, float(tol))
    rscale, e_id = _scale(H)
    es = np.array([p[0] for p in acc], dtype=complex)
    rs = np.array([p[1] / (rscale + abs(p[0] - e_id)) for p in acc])
    if not vectors: return es, rs
    return es, rs, [_combine(self, snaps, p[2]) for p in acc]


def ground_pair(self, H=None, nsnap=NSNAP, tau=None, nsteps=NSTEPS,
                tol=RITZ_TOL):
    """(E, psil, psir): the level of smallest real part of H with its right
    and left eigenvectors, <psil|psir> = 1 (module docstring). Raises
    RuntimeError when no Ritz pair is accepted or the two vectors do not
    pair."""
    if H is None: H = self.hamiltonian
    snaps = _right_states(self, H, tau, nsteps, int(nsnap))
    acc = _accepted(_rayleigh_ritz(snaps, H), H, float(tol))
    if not acc:
        raise RuntimeError("nhlevels: no Ritz pair of the imaginary-time "
                           "states reached the residual tolerance")
    theta, _, c = acc[0]
    psir = _combine(self, snaps, c)
    psir = psir * (1.0 / np.sqrt(psir.dot(psir).real))
    # the left state: evolved under H^dagger from conj(psir), the left
    # eigenvector itself when H^T = H
    Hd = H.get_dagger()
    s = self._session
    Hdop = s.build_operator(Hd.to_terms())
    lsnaps = _evolve(self, s.conjugate(psir.cpp_handle), Hdop,
                     _grid(H, tau, nsteps), int(nsnap))
    lpairs = _rayleigh_ritz(lsnaps, Hd)
    rscale, e_id = _scale(H)
    window = 1e-4 * (rscale + abs(theta - e_id))
    near = [p for p in lpairs if abs(p[0] - np.conj(theta)) <= window]
    if near:
        # within a level shared by two sectors, the left Ritz vector that
        # pairs best with psir
        cands = [_combine(self, lsnaps, p[2]) for p in near]
        psil = max(cands, key=lambda l: abs(l.dot(psir)) / np.sqrt(l.dot(l).real))
    else:
        p = min(lpairs, key=lambda p: abs(p[0] - np.conj(theta)))
        psil = _combine(self, lsnaps, p[2])
    ov = psil.dot(psir)
    if not (abs(ov) > 1e-10 * np.sqrt(psil.dot(psil).real)):
        raise RuntimeError("nhlevels: the left and right Ritz vectors do "
                           "not pair (<psil|psir> = %r)" % (ov,))
    # <psil|psir> = 1, the factor split as nhdmrg()'s pair has it
    sq = np.sqrt(complex(ov))
    psil = psil * (1.0 / np.conj(sq))
    psir = psir * (1.0 / sq)
    E = psil.dot(H * psir) / psil.dot(psir)
    return complex(E), psil, psir


def gs_energy_imaginary_time(self, **kwargs):
    """gs_energy() for a non-Hermitian Hamiltonian by the imaginary-time
    route (nh_gs_method="imaginary_time"): ground_pair(), nhdmrg()'s
    two-sided certificate at its tol= (1e-4 by default), and the pair
    stored as NH-DMRG's is (nhdmrg.store_nh_pair). When the route raises
    or the certificate fails, it falls back to NH-DMRG with the same
    keywords, so no chain NH-DMRG solves is worse off. H= is refused as
    on the NH-DMRG route, through it."""
    from . import nhdmrg as fnh
    if "H" in kwargs or not available(self):
        return fnh.gs_energy_nhdmrg(self, **kwargs)
    tol = float(kwargs.get("tol", 1e-4))
    H = self.hamiltonian
    try:
        E, psil, psir = ground_pair(self, H)
    except RuntimeError as e:
        if self.verbose > 0:
            print("nhlevels: imaginary-time route failed (%r), using NH-DMRG" % (e,))
        return fnh.gs_energy_nhdmrg(self, **kwargs)
    Hc, Hcd, up, rscale, e_id = fnh.certificate_scales(H)
    resid = fnh._certificate(Hc, Hcd, up, rscale, e_id, E, psil, psir)
    if not (resid < tol):
        if self.verbose > 0:
            print("nhlevels: imaginary-time pair at relative residual",
                  resid, "above", tol, "- using NH-DMRG")
        return fnh.gs_energy_nhdmrg(self, **kwargs)
    self._nh_last_residual = (resid, tol)
    return fnh.store_nh_pair(self, E, psil, psir)
