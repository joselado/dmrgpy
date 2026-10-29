"""Spectrum of a generator from a Krylov space of time-evolved MPS.

Port to this engine of the eigensolver of Westhoff, Moroder, Schollwoeck
and Paeckel, arXiv:2509.07709 (CLIK-MPS), for any generator G given as an
MPO, Hermitian or not: the slowest-decaying eigenvalues of G, those of
largest real part. For a Hamiltonian, G = -H and they are the lowest
energies; for a Lindbladian, G = L and they are the steady state and the
slowest modes.

The Krylov space is not built from powers of G but from snapshots of two
time evolutions, rho(m) = exp(G z m) rho(0), m = 0..nk-1:

  * a real contour, z = dt, from one random start, which resolves the
    slowest mode, and
  * a complex contour, z = dt exp(-i alpha), from a second random start.
    Along it mode k decays at the rate
        xi_k = -Re(lambda_k) cos(alpha) - Im(lambda_k) sin(alpha),
    so modes with Im(lambda_k) > 0 decay slower and stay in the space long
    enough to be resolved. Above an angle alpha_max, set by the bulk of the
    spectrum, some mode overtakes the slowest one and the branch is
    useless; `probe_alpha_max` estimates it, and the default angle is a
    quarter of it.

Each evolution step is one tdvp_step with dt -> i z (tdvp_step applies
exp(-i dt M), and the MPO passed is G itself, so -i (i z) G = z G). The
exponentiator is decided once for G (Lanczos when G is Hermitian, Arnoldi
otherwise). Every snapshot is stored normalized: the span does not
change, and the norms of exp(G z m) would otherwise overflow or underflow.

The projection: with M_ij = <rho_i|rho_j> and G_ij = <rho_i|G|rho_j>, M is
diagonalized, directions below gram_tol times its largest eigenvalue are
dropped, X = U S^{-1/2} spans the rest orthonormally, and the small
non-Hermitian G_eff = X^dag G X is diagonalized. An eigenvector is the
coefficient list c = X R over the stored states and is never summed into
one MPS: overlaps and expectation values are read off matrix elements.
X amplifies the error of every G_ij by 1/S, so gram_tol also sets how
accurate the eigenvalues can be.

The matrix elements are grown one snapshot at a time (Section II.5 of the
paper): a new snapshot costs one overlap and two elements of G with each
stored state, linear work per step, and nothing is recomputed when the
second branch is added.

Convergence is judged as in the paper, by tracking the spectrum. The
eigenvalues are solved a second time on the space without the last
snapshot of each branch, which only restricts the matrices already
computed (Appendix G), and the change of each eigenvalue is its
convergence estimate. A converged mode does not move; an unconverged one
does. It orders the modes by convergence correctly but is not an error
bar: on 6-site chains against ED it was 3 to 40 times smaller than the
true error of the eigenvalue.

The residual of an eigenpair, ||G r - lambda r|| with ||r|| = 1, is
available on request (exact_residuals=True). It is c^dag (W - lambda G^dag
- conj(lambda) G + |lambda|^2 M) c with W_ij = <rho_i|G^dag G|rho_j> over
every pair, through the MPO of G^dag G, whose bond dimension is the square
of G's; that makes it the costliest part of a run, which is why the paper
uses <<rho|L^dag L|rho>> only for the steady-state error bound. It is also
a difference of numbers of order ||G||^2, so it resolves the residual only
down to about 1e-8 of the operator scale.

Device: the overlaps go through mpsalgebra.inner, on the active array
backend; the 2nk x 2nk algebra is host work, as in task.py.
"""

import numpy as np

from . import backend as bk
from .mpsalgebra import inner, mpo_dagger, mpo_is_hermitian, nmultMPO, randomMPS
from .gse import global_subspace_expand
from .tdvp import tdvp_step


def _scalar(x):
    return complex(bk.to_host(x)) if not isinstance(x, complex) else x


def mpo_product(A, B):
    """The MPO of the operator product A B (B acts first), exact.
    nmultMPO(X, prime(Y)) followed by mapPrime(2, 1) is Y X (X acts
    first), measured against MPS applications, hence the swap."""
    out = nmultMPO(B, A.copy().prime(), cutoff=0.0, maxdim=None)
    out.mapPrime(2, 1)
    return out


def _step(psi, G, z, cutoff, maxdim, niter, hermitian, one_site=False):
    """exp(G z)|psi>, normalized, and the log of the norm it had.
    one_site=True is one-site TDVP after a global subspace expansion of
    every step, the pairing of the paper; otherwise two-site TDVP.

    one_site=True is experimental and off everywhere: on a generator with a
    conserved charge it leaks out of the sector. On the seeded 3-site
    Bose-Hubbard Lindbladian the middle bond grows past the in-sector rank
    of 14 to 27, and the traceless branch is 0.62 off the two-site one after
    five steps, for every expansion cutoff from 1e-10 to 1e-4. One-site TDVP
    without the expansion does not leak (it agrees with two-site TDVP to
    integrator error). The paper avoids this with block-sparse U(1)
    symmetries; see docs/clik_mps_lindblad_plan.md, "Status 2026-10-30"."""
    if one_site:
        with bk.pad_bonds_suspended():
            psi = global_subspace_expand(G, psi, _GSE_ORDER, cutoff,
                                         maxdim=maxdim, bond_maxdim=maxdim)
            out = tdvp_step(psi.copy(), G, 1j * z, cutoff=cutoff,
                            maxdim=maxdim, niter=niter, num_center=1,
                            hermitian=hermitian)
    else:
        out = tdvp_step(psi.copy(), G, 1j * z, cutoff=cutoff, maxdim=maxdim,
                        niter=niter, hermitian=hermitian)
    nrm = float(np.sqrt(abs(_scalar(inner(out, out)))))
    out.normalize()
    return out, np.log(nrm)


# Steps of each branch preceded by a global subspace expansion, and its
# Krylov order. Two-site TDVP cannot move a product state under terms that
# couple next-nearest neighbors, which on the doubled chain of a
# Lindbladian is every hopping and every dissipator: the projected
# generator is then diagonal, and a seeded product start never moved
# (two eigenvalues, residual 10, on the 5-site Bose-Hubbard check).
# Expanding with G|psi> and G^2|psi> first gives the bonds the directions
# those terms need; after a few steps two-site TDVP grows them itself.
_GSE_STEPS = 3
_GSE_ORDER = 3


def evolve_branch(G, psi0, z, nk, cutoff, maxdim, niter=50, hermitian=None,
                  one_site=False):
    """nk normalized snapshots exp(G z m)|psi0>, m = 0..nk-1, with a
    global subspace expansion before each of the first _GSE_STEPS steps
    (before every step for one_site=True)."""
    if hermitian is None:
        hermitian = mpo_is_hermitian(G)
    psi = psi0.copy()
    psi.normalize()
    states = [psi]
    for m in range(nk - 1):
        if m < _GSE_STEPS and not one_site:
            with bk.pad_bonds_suspended():
                psi = global_subspace_expand(G, psi, _GSE_ORDER, cutoff,
                                             maxdim=maxdim, bond_maxdim=maxdim)
        psi, _ = _step(psi, G, z, cutoff, maxdim, niter, hermitian, one_site)
        states.append(psi)
    return states


def probe_alpha_max(G, psi0, lam0, dt, cutoff, maxdim, nsteps=8,
                    angles=None, niter=50, hermitian=None, one_site=False):
    """The first angle of `angles` at which the complex contour stops being
    stable, or pi/2 if none does (the paper's Algorithm 3).

    Along z = dt exp(-i a), the log-norm of exp(G z m)|psi0> is shifted by
    Re(lam0 z m), lam0 the slowest eigenvalue found on the real contour, so
    that the slowest mode stays at a constant shifted norm. In a stable
    run every other mode decays and the shifted log-norm settles from
    below; once a bulk mode overtakes it, the shifted log-norm grows. The
    angle is taken as unstable when the mean shifted increment over the
    second half of the nsteps steps is positive. This is a heuristic, as in
    the paper: a short window can miss a slowly overtaking mode, which is
    why the default angle is a quarter of what this returns.

    For a Hermitian generator every Im(lambda) vanishes, nothing overtakes,
    and pi/2 comes back; the angle then only scales the damping."""
    if hermitian is None:
        hermitian = mpo_is_hermitian(G)
    if angles is None:
        angles = np.pi / 2 * np.arange(1, 8) / 8.0
    for a in angles:
        z = dt * np.exp(-1j * a)
        psi = psi0.copy()
        psi.normalize()
        increments = []
        for _ in range(nsteps):
            psi, lognorm = _step(psi, G, z, cutoff, maxdim, niter, hermitian,
                                 one_site)
            increments.append(lognorm - (lam0 * z).real)
        if np.mean(increments[nsteps // 2:]) > 0.0:
            return float(a)
    return np.pi / 2


class KrylovMatrices:
    """The stored states with M_ij = <rho_i|rho_j> and G_ij =
    <rho_i|G|rho_j>, grown one state at a time."""

    def __init__(self, G):
        self.G = G
        self.states = []
        self.M = np.zeros((0, 0), dtype=complex)
        self.Gm = np.zeros((0, 0), dtype=complex)

    def add(self, psi):
        n = len(self.states)
        self.states.append(psi)
        M = np.zeros((n + 1, n + 1), dtype=complex)
        Gm = np.zeros((n + 1, n + 1), dtype=complex)
        M[:n, :n] = self.M
        Gm[:n, :n] = self.Gm
        for i, rho in enumerate(self.states):
            M[i, n] = _scalar(inner(rho, psi))
            M[n, i] = np.conj(M[i, n])
            Gm[i, n] = _scalar(inner(rho, self.G, psi))
            if i < n:
                Gm[n, i] = _scalar(inner(psi, self.G, rho))
        self.M, self.Gm = M, Gm


def matrix_elements(states, A=None):
    """<rho_i|A|rho_j> over the stored states (A=None: the overlaps)."""
    n = len(states)
    out = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            if A is None:
                if j < i:
                    out[i, j] = np.conj(out[j, i])
                    continue
                out[i, j] = _scalar(inner(states[i], states[j]))
            else:
                out[i, j] = _scalar(inner(states[i], A, states[j]))
    return out


def project(M, Gm, gram_tol):
    """Eigenvalues of G on the span, sorted by decreasing real part, and
    their coefficient vectors over the stored states (columns with
    c^dag M c = 1)."""
    M = (M + M.conj().T) / 2
    s, U = np.linalg.eigh(M)
    keep = s > gram_tol * s.max()
    X = U[:, keep] / np.sqrt(s[keep])
    Geff = X.conj().T @ Gm @ X
    lam, R = np.linalg.eig(Geff)
    order = np.argsort(-lam.real, kind="stable")
    return lam[order], X @ R[:, order], int(keep.sum())


def eigenvalue_changes(lam, M, Gm, keep, gram_tol):
    """For each eigenvalue, the distance to the nearest eigenvalue of the
    space restricted to the states `keep`: the paper's tracking of the
    spectrum, with no new matrix element."""
    sub = project(M[np.ix_(keep, keep)], Gm[np.ix_(keep, keep)], gram_tol)[0]
    return np.array([np.min(np.abs(sub - l)) for l in lam])


def residuals(lam, C, M, Gm, W):
    """||G r - lambda r|| for each normalized eigenvector r = sum c_i rho_i;
    see the module docstring for its floor."""
    out = []
    for k in range(len(lam)):
        c = C[:, k]
        l = lam[k]
        Q = W - l * Gm.conj().T - np.conj(l) * Gm + abs(l) ** 2 * M
        out.append(np.sqrt(max(np.real(np.vdot(c, Q @ c)), 0.0)))
    return np.array(out)


class KrylovSpectrum:
    """What krylov_spectrum returns: the eigenvalues, their residuals, and
    the stored states with the matrix elements, so that an expectation value
    in each eigenvector costs one more matrix of overlaps."""

    def __init__(self, lam, C, res, states, M, ndim, alpha, changes=None):
        self.eigenvalues = lam
        self.coefficients = C
        self.residuals = res      # None unless exact_residuals=True
        self.changes = changes    # the tracking estimate, always there
        self.states = states
        self.M = M
        self.ndim = ndim    # Krylov dimension kept after deflation
        self.alpha = alpha  # contour angle used

    def expectation(self, A, k):
        """<r_k|A|r_k>/<r_k|r_k> for an MPO A."""
        c = self.coefficients[:, k]
        Am = matrix_elements(self.states, A)
        return np.vdot(c, Am @ c) / np.vdot(c, self.M @ c)


def krylov_spectrum(G, sites, nk=20, dt=0.5, alpha=None, maxdim=64,
                    cutoff=1e-10, gram_tol=1e-12, seed=0, start_maxdim=10,
                    niter=50, starts=None, exact_residuals=False,
                    one_site=False):
    """The slowest eigenvalues of the MPO G (largest real part first) from
    the two-branch Krylov space described in the module docstring.

    nk snapshots per branch, time step dt, contour angle alpha of the
    second branch (None: a quarter of probe_alpha_max; 0 puts it on the
    real contour too, from its own random start), maxdim and
    cutoff for the TDVP evolution, gram_tol the relative deflation
    threshold of the Gram matrix, seed for the two random starts. `starts`,
    a pair of MPS, replaces the two random starts (real branch first): the
    way to keep the space inside one sector of a conserved quantity of G.
    exact_residuals=True adds ||G r - lambda r|| for every eigenvector (see
    the module docstring for its cost).

    The defaults are measured, on 6-site chains against ED with the
    hopping as the unit. A small dt makes consecutive imaginary-time
    snapshots nearly parallel, and deflation then removes most of them:
    at dt=0.1 the lowest levels were off by 1e-3 or more, at dt=0.6 with
    nk=24 and gram_tol=1e-12 by 2e-12 to 5e-9 (Heisenberg) and 5e-7 to
    8e-6 (asymmetric hopping). At equal Krylov size the contour angle
    changed nothing measurable on those Hamiltonian chains; its gain in
    the paper is for Lindbladian modes that oscillate strongly next to
    their decay."""
    if nk < 1:
        raise ValueError("krylov_spectrum: nk must be at least 1, got %r" % (nk,))
    if dt <= 0:
        raise ValueError("krylov_spectrum: dt must be positive, got %r" % (dt,))
    rng = np.random.default_rng(seed)
    herm = mpo_is_hermitian(G)
    m0 = max(1, min(int(maxdim), int(start_maxdim)))
    if starts is None:
        psi_p = randomMPS(sites, m0, rng=rng)
        psi_u = randomMPS(sites, m0, rng=rng)
    else:
        psi_p, psi_u = starts
    km = KrylovMatrices(G)
    for psi in evolve_branch(G, psi_p, dt, nk, cutoff, maxdim, niter, herm,
                             one_site):
        km.add(psi)
    if alpha is None:
        lam0 = project(km.M, km.Gm, gram_tol)[0][0]
        alpha = 0.25 * probe_alpha_max(G, psi_u, lam0, dt, cutoff, maxdim,
                                       niter=niter, hermitian=herm,
                                       one_site=one_site)
    for psi in evolve_branch(G, psi_u, dt * np.exp(-1j * alpha), nk,
                             cutoff, maxdim, niter, herm, one_site):
        km.add(psi)
    lam, C, ndim = project(km.M, km.Gm, gram_tol)
    keep = [i for i in range(2 * nk) if i not in (nk - 1, 2 * nk - 1)]
    changes = (eigenvalue_changes(lam, km.M, km.Gm, keep, gram_tol)
               if keep else np.full(len(lam), np.inf))
    res = None
    if exact_residuals:
        W = matrix_elements(km.states, mpo_product(mpo_dagger(G), G))
        res = residuals(lam, C, km.M, km.Gm, W)
    return KrylovSpectrum(lam, C, res, km.states, km.M, ndim, float(alpha),
                          changes)
