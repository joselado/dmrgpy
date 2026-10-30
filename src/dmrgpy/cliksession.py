"""CLIK-MPS (arXiv:2509.07709) on a compiled session, itensor_version=3.

The algorithm, its parameters, its measured defaults and the floor of its
residuals are in pyitensor/clik.py, which is the reference: this module
runs the same two-branch Krylov space, the same adaptive sub-stepping, the
same angle probe and the same projection, with every MPS and MPO operation
going to the v3 session (two-site TDVP, global subspace expansion,
overlaps and sums) and the small projected problem done by clik.py's own
host functions. Only the orchestration is repeated here; the heavy work is
C++.

The local exponentials follow pyitensor's choice: a Hermitian G steps
with the session's tdvp_step (ITensor's Lanczos applyExp), any other G
with tdvp_arnoldi_step (mpscpp3/cbe.h), since the Lanczos applyExp does
not converge on a Lindbladian: on the 3-site Bose-Hubbard one at dt=0.5
it stalls at its iteration cap and the next SVD meets a NaN. The random
starts are Gaussian MPS from a generator seeded with `seed`
(clik_random_starts), so a run is reproducible on v3 but its starts are
not those of "python", and its unconverged levels differ by draw.
"""

import numpy as np

from .pyitensor import clik as _clik


class Session:
    """The MPS/MPO operations the driver needs, on a v3 chain."""

    def __init__(self, chain, hermitian=True):
        self.chain = chain
        self.hermitian = hermitian
        self.s = chain._session
        self.maxm = chain.maxm
        self.cutoff = chain.cutoff

    def mpo(self, op):
        return self.chain.toMPO(op).cpp_handle

    def inner(self, a, b):
        return complex(self.s.overlap(a, b))

    def sandwich(self, a, A, b):
        return complex(self.s.overlap_aMb_operator(a, A, b))

    def norm(self, psi):
        return float(np.sqrt(abs(self.inner(psi, psi))))

    def scale(self, psi, z):
        return self.s.scale_mps(psi, z)

    def normalized(self, psi):
        return self.scale(psi, 1.0 / self.norm(psi))

    def add(self, a, b):
        return self.s.sum_mps(a, b)

    def step(self, G, psi, z, niter):
        # both steps apply exp(-i dt G); dt = i z gives exp(z G), and
        # since dt is then not real neither normalizes. A non-Hermitian G
        # takes the Arnoldi step: the Lanczos applyExp inside tdvp_step
        # does not converge on a Lindbladian
        if self.hermitian:
            return self.s.tdvp_step(G, psi, 1j * z, 2, niter)
        return self.s.tdvp_arnoldi_step(G, psi, 1j * z, niter)

    def expand(self, G, psi):
        return self.s.global_subspace_expand(G, psi, _clik._GSE_ORDER,
                                             self.cutoff, self.maxm)

    def dagger(self, G):
        return self.s.hermitian_operator(G)

    def product(self, first, then):
        """The MPO of `then` after `first` (first acts first)."""
        return self.s.multiply_operators(first, then)

    def apply(self, A, psi):
        return self.s.apply_pure_operator(A, psi)


def _step(ops, psi, G, z, niter, substeps=1):
    """clik._step on the session: exp(G z)|psi>, normalized, and the log
    of its norm, as `substeps` TDVP steps."""
    lognorm = 0.0
    for _ in range(substeps):
        out = ops.step(G, psi, z / substeps, niter)
        nrm = ops.norm(out)
        if not np.isfinite(nrm) or nrm == 0:
            return out, np.inf
        psi = ops.scale(out, 1.0 / nrm)
        lognorm += np.log(nrm)
    return psi, lognorm


def _deflated(ops, psi, d):
    """clik._deflated: psi - (<d|psi>/<d|d>) d, or psi itself when that
    component is below 1e-14 of it."""
    dd = ops.inner(d, d)
    c = ops.inner(d, psi) / dd
    if abs(c) * np.sqrt(abs(dd)) <= 1e-14 * ops.norm(psi):
        return psi
    return ops.add(psi, ops.scale(d, -c))


def _step_accepted(ops, meter, old, new, lognorm):
    """clik._step_accepted on the session."""
    if not np.isfinite(lognorm) or lognorm > _clik._GROWTH_TOL:
        return False
    if meter is None:
        return True
    before = ops.inner(meter, old)
    after = ops.inner(meter, new) * np.exp(lognorm)
    scale = np.sqrt(abs(ops.inner(meter, meter)))
    return abs(after - before) <= _clik._METER_TOL * scale


def _branch_attempt(ops, G, psi0, z, nk, niter, deflate, substeps, meter,
                    gated):
    """clik._branch_attempt: the snapshots, or None at the first rejected
    step, so that the branch restarts with more sub-steps."""
    psi = ops.normalized(psi0)
    states = [psi]
    for m in range(nk - 1):
        if m < _clik._GSE_STEPS:
            psi = ops.expand(G, psi)
        try:
            new, lognorm = _step(ops, psi, G, z, niter, substeps)
        except (RuntimeError, ValueError, FloatingPointError):
            if not gated:
                raise
            return None   # overflowed inside the sweep
        if gated and not _step_accepted(ops, meter, psi, new, lognorm):
            return None
        psi = new
        if deflate is not None:
            psi = ops.normalized(_deflated(ops, psi, deflate))
        states.append(psi)
    return states


def evolve_branch(ops, G, psi0, z, nk, niter=50, deflate=None, substeps=1,
                  meter=None, bounded=False):
    """clik.evolve_branch on the session: (states, substeps)."""
    gated = meter is not None or bounded
    while True:
        states = _branch_attempt(ops, G, psi0, z, nk, niter, deflate,
                                 substeps, meter, gated)
        if states is not None:
            return states, substeps
        if substeps >= _clik._MAX_SUBSTEPS:
            raise RuntimeError(
                "krylov_spectrum: a step of |z|=%.3g is still unstable with "
                "%d TDVP sub-steps; reduce dt" % (abs(z), substeps))
        substeps *= 2


def probe_alpha_max(ops, G, psi0, lam0, dt, nsteps=8, angles=None, niter=50,
                    substeps=1):
    """clik.probe_alpha_max on the session, the same bisection."""
    if angles is None:
        angles = np.pi / 2 * np.arange(1, 8) / 8.0

    def unstable(a):
        z = dt * np.exp(-1j * a)
        psi = ops.normalized(psi0)
        increments = []
        for _ in range(nsteps):
            psi, lognorm = _step(ops, psi, G, z, niter, substeps)
            increments.append(lognorm - (lam0 * z).real)
        return np.mean(increments[nsteps // 2:]) > 0.0

    if len(angles) == 0 or unstable(angles[0]):
        return float(angles[0]) if len(angles) else np.pi / 2
    lo, hi = 1, len(angles)
    while lo < hi:
        mid = (lo + hi) // 2
        if unstable(angles[mid]):
            hi = mid
        else:
            lo = mid + 1
    return float(angles[lo]) if lo < len(angles) else np.pi / 2


class KrylovMatrices:
    """clik.KrylovMatrices on the session: M_ij = <rho_i|rho_j> and
    G_ij = <rho_i|G|rho_j>, grown one state at a time, the row of G from
    G^dag unless G is Hermitian."""

    def __init__(self, ops, G, hermitian):
        self.ops = ops
        self.G = G
        self.hermitian = hermitian
        self.Gd = None if hermitian else ops.dagger(G)
        self.states = []
        self.M = np.zeros((0, 0), dtype=complex)
        self.Gm = np.zeros((0, 0), dtype=complex)

    def add(self, psi):
        ops = self.ops
        n = len(self.states)
        self.states.append(psi)
        col_M = np.array([ops.inner(r, psi) for r in self.states])
        col_G = np.array([ops.sandwich(r, self.G, psi) for r in self.states])
        if self.hermitian:
            row_G = np.conj(col_G)
        else:
            row_G = np.conj([ops.sandwich(r, self.Gd, psi) for r in self.states])
        M = np.zeros((n + 1, n + 1), dtype=complex)
        Gm = np.zeros((n + 1, n + 1), dtype=complex)
        M[:n, :n] = self.M
        Gm[:n, :n] = self.Gm
        M[:, n] = col_M
        M[n, :n] = np.conj(col_M[:n])
        Gm[:, n] = col_G
        Gm[n, :n] = row_G[:n]
        self.M, self.Gm = M, Gm


def matrix_elements(ops, states, A):
    n = len(states)
    out = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            out[i, j] = ops.sandwich(states[i], A, states[j])
    return out


class KrylovSpectrum(_clik.KrylovSpectrum):
    """clik.KrylovSpectrum whose expectation values go through the
    session; A is an MPO handle of the session."""

    def __init__(self, ops, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ops = ops

    def expectation(self, A, k):
        c = self.coefficients[:, k]
        Am = matrix_elements(self.ops, self.states, A)
        return np.vdot(c, Am @ c) / np.vdot(c, self.M @ c)


def is_left_zero_mode(ops, G, d, rtol=1e-10):
    """pyitensor Chain's _is_left_zero_mode: ||G^dag d|| <= rtol ||G d||."""
    nx = ops.norm(ops.apply(ops.dagger(G), d))
    ny = ops.norm(ops.apply(G, d))
    return nx <= rtol * max(ny, ops.norm(d))


def krylov_spectrum(chain, Gop, hermitian, nk=20, dt=0.5, alpha=None,
                    gram_tol=1e-12, seed=0, start_maxdim=10, niter=50,
                    starts=None, exact_residuals=False, traceless=False,
                    sector_charge=None, substeps=1):
    """clik.krylov_spectrum, together with pyitensor Chain.krylov_spectrum's
    choice of the deflation vector, for the MultiOperator Gop on the v3
    chain `chain`. `hermitian` is whether Gop is Hermitian (the caller
    proves it on the MultiOperator). `starts` is a pair of lists of local
    vectors, one per site. Returns a KrylovSpectrum."""
    if nk < 1:
        raise ValueError("krylov_spectrum: nk must be at least 1, got %r" % (nk,))
    if dt <= 0:
        raise ValueError("krylov_spectrum: dt must be positive, got %r" % (dt,))
    ops = Session(chain, hermitian)
    s = chain._session
    G = ops.mpo(Gop)
    deflate = None
    sector = False
    if traceless:
        deflate = s.doubled_identity_mps()
        if sector_charge is not None:
            d = s.sector_identity_mps(int(sector_charge))
            if is_left_zero_mode(ops, G, d):
                deflate = d
                sector = True
    if starts is None:
        m0 = max(1, min(int(ops.maxm), int(start_maxdim)))
        psi_p, psi_u = s.clik_random_starts(m0, int(seed))
    else:
        psi_p, psi_u = (s.product_mps([list(np.asarray(v, dtype=complex))
                                       for v in vecs]) for vecs in starts)
    if deflate is not None:
        psi_u = _deflated(ops, psi_u, deflate)
    km = KrylovMatrices(ops, G, hermitian)
    meter = None
    if deflate is not None and abs(ops.inner(deflate, psi_p)) > 1e-8 * (
            ops.norm(deflate) * ops.norm(psi_p)):
        meter = deflate
    states, substeps = evolve_branch(ops, G, psi_p, dt, nk, niter,
                                     substeps=substeps, meter=meter,
                                     bounded=deflate is not None)
    for psi in states:
        km.add(psi)
    if alpha is None and hermitian:
        alpha = 0.25 * np.pi / 2   # see clik.krylov_spectrum
    elif alpha is None:
        lam0 = _clik.project(km.M, km.Gm, gram_tol)[0][0]
        alpha = 0.25 * probe_alpha_max(ops, G, psi_u, lam0, dt, niter=niter,
                                       substeps=substeps)
    states, _ = evolve_branch(ops, G, psi_u, dt * np.exp(-1j * alpha), nk,
                              niter, deflate, substeps=substeps,
                              bounded=deflate is not None)
    for psi in states:
        km.add(psi)
    lam, C, ndim = _clik.project(km.M, km.Gm, gram_tol)
    keep = [i for i in range(2 * nk) if i not in (nk - 1, 2 * nk - 1)]
    changes = (_clik.eigenvalue_changes(lam, km.M, km.Gm, keep, gram_tol)
               if keep else np.full(len(lam), np.inf))
    res = None
    if exact_residuals:
        W = matrix_elements(ops, km.states, ops.product(G, ops.dagger(G)))
        res = _clik.residuals(lam, C, km.M, km.Gm, W)
    out = KrylovSpectrum(ops, lam, C, res, km.states, km.M, ndim,
                         float(alpha), changes)
    out.substeps = substeps
    out.sector_deflation = sector
    return out
