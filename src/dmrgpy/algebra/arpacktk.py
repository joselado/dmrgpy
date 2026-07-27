"""Matrix-free, pure-Python Implicitly Restarted Arnoldi Method (IRAM),
adapted from ARPACK (https://bitbucket.org/chaoyang2013/arpack, SRC/zn*.f
-- the double-complex non-symmetric driver znaupd/znaup2/znaitr/znapps/
zneigh/zngets/zsortc) so it runs directly against dmrgpy's own MPS/ED
wavefunction objects instead of flat numpy arrays, with no BLAS/LAPACK
Fortran dependency and no reverse-communication interface (there is no
separate process here, so Op(x) is just called directly).

Why the *complex* ARPACK routines, not the real double-precision ones
(dnaupd/...)?  dmrgpy's Hamiltonians and wavefunctions are complex-valued
already (non-Hermitian models routinely have complex spectra), so there is
nothing to gain from ARPACK's real-arithmetic double-shift trick (which
exists purely to keep dnaupd's *dense linear algebra* in real BLAS/LAPACK
for speed on flat arrays). Working in complex arithmetic throughout lets
each implicit-restart shift be a single (possibly complex) Ritz value,
applied with a single-shift QR step -- no need for dnapps's paired
double-shift bookkeeping to keep complex-conjugate Ritz pairs together.

Why numpy.linalg.qr instead of porting znapps's Givens-rotation bulge
chase? Sorensen's 1992 "Implicit Application of Polynomial Filters in a
k-Step Arnoldi Method" (the paper both dnapps and znapps implement) shows
that applying a shift mu via a single-shift QR step of the *whole* current
Hessenberg matrix H (Q,R = qr(H-mu*I); H <- Q^H H Q) produces exactly the
same orthogonal transform as ARPACK's O(ncv^2) hand-rolled bulge chase --
bulge chasing is only ARPACK's dense-linear-algebra *implementation*
detail for that identity, not part of the mathematical algorithm. Since
the Krylov subspace size ncv here is tiny (order 10-50) while a single
Op(x) is an MPO*MPS application, the O(ncv^3) cost of a plain dense QR per
shift is negligible -- there is nothing to gain from bulge chasing's extra
implementation complexity. Applying each shift's QR to the *already
transformed* H (not to a freshly-formed product polynomial prod(H-mu_i*I))
keeps the conditioning equivalent to sequential bulge chasing rather than
the numerically worse direct-polynomial formulation.

Ported pieces, and their ARPACK sources:
  arnoldi_extend  <- znaitr.f  (build/extend a length-k factorization by
                      classical Gram-Schmidt with iterative refinement,
                      the 0.717 ratio test and 2-pass retry policy taken
                      directly from znaitr's STEP5/STEP6 comments)
  hessenberg_ritz <- zneigh.f  (eigenvalues + Ritz estimates of H;
                      zneigh needs zlahqr+ztrevc to get eigenvectors of a
                      real-arithmetic-friendly Schur form -- here
                      numpy.linalg.eig gives the same normalized
                      eigenvectors of the small complex H directly)
  select_shifts   <- zngets.f/zsortc.f (sort so the "wanted" Ritz values
                      per `which` end up last; unwanted ones become
                      shifts, most-uncertain-first)
  apply_shifts    <- znapps.f (implicit restart -- see above for why QR
                      replaces the bulge chase)
  mpsiram         <- znaupd.f/znaup2.f (the outer driver: extend, check
                      convergence, adjust kev to avoid stagnation exactly
                      as znaup2 does, restart, repeat)

This is deliberately the *matrix-free, mode-1* subset of ARPACK (OP=A,
B=I; no shift-invert modes 2-4, no reverse communication, no B-inner
product) -- everything dmrgpy's own MPS/ED objects need
(``self.random_mps()``/``EDchain``, ``.dot()``, ``.normalize()``,
``.copy()``, ``+``/``-``/scalar ``*``) and nothing more.
"""
import numpy as np

from .arnolditk import random_state
from . import arnolditk


def _goodness(which, ritz):
    """Score used to sort Ritz values so the *wanted* ones (per ARPACK's
    which={LM,SM,LR,SR,LI,SI} convention) end up at the end of the array
    -- see zngets.f/zsortc.f: after sorting ascending by this score, the
    last kev entries are wanted, the first np are unwanted (shift
    candidates)."""
    if which == "LM": return np.abs(ritz)
    if which == "SM": return -np.abs(ritz)
    if which == "LR": return ritz.real
    if which == "SR": return -ritz.real
    if which == "LI": return ritz.imag
    if which == "SI": return -ritz.imag
    raise ValueError("Unknown which='%s'" % which)


def arnoldi_extend(self, Op, V, H, resid, rnorm, k, npsteps,
        wfskip=None, verbose=0):
    """Port of znaitr: extend a length-k Arnoldi factorization
    (Op*V_k - V_k*H_k = resid*e_k^T, V_k^H*resid = 0, ||resid|| = rnorm)
    by npsteps more steps, using only Op(x) plus the generic wavefunction
    operations (dot/normalize/+/-/scalar*). Every new vector costs
    exactly one Op(x) call, same as arnolditk's own build_arnoldi_chain,
    so existing Op(x)-counting instrumentation (see
    examples/groundstate/arnoldi_benchmark) applies unmodified."""
    if wfskip is None: wfskip = []
    kplusp = k + npsteps
    Hnew = np.zeros((kplusp, kplusp), dtype=complex)
    if k > 0: Hnew[:k, :k] = H
    Vout = list(V)
    beta = rnorm
    for j in range(k, kplusp):
        if beta < 1e-12:
            # invariant subspace hit mid-extension: restart from a fresh
            # random vector orthogonal to the basis built so far (mirrors
            # znaitr's dgetv0-based restart, minus the give-up-after-3-
            # tries limit -- random_state() already loops until it finds
            # a usable vector)
            if verbose > 0: print("arnoldi_extend: invariant subspace, restarting")
            resid = random_state(self, orthogonal=Vout + wfskip)
            beta = 1.0
        vj = resid * (1. / beta)
        Vout.append(vj)
        w = Op(vj)
        if wfskip:
            for s in wfskip:
                w = w - s.dot(w) * s
        coeffs = np.array([Vout[i].dot(w) for i in range(j + 1)], dtype=complex)
        r = w
        for i in range(j + 1): r = r - coeffs[i] * Vout[i]
        if j > 0: Hnew[j, j - 1] = beta
        Hnew[:j + 1, j] = coeffs
        wnorm = np.sqrt(abs(w.dot(w)))
        rnorm1 = np.sqrt(abs(r.dot(r)))
        if wnorm > 0 and rnorm1 <= 0.717 * wnorm:
            # iterative refinement (re-orthogonalization), up to 2 passes
            # -- znaitr's STEP5, same 0.717 acceptance ratio and 2-try
            # give-up policy
            converged = False
            for _try in range(2):
                if wfskip:
                    for s in wfskip:
                        r = r - s.dot(r) * s
                s = np.array([Vout[i].dot(r) for i in range(j + 1)], dtype=complex)
                for i in range(j + 1): r = r - s[i] * Vout[i]
                Hnew[:j + 1, j] += s
                rnorm2 = np.sqrt(abs(r.dot(r)))
                prev = rnorm1
                rnorm1 = rnorm2
                if rnorm2 > 0.717 * prev:
                    converged = True
                    break
            if not converged:
                # r lies numerically in span(V): treat as invariant subspace
                r = r * 0.0
                rnorm1 = 0.0
        resid = r
        beta = rnorm1
    return Vout, Hnew, resid, beta


def hessenberg_ritz(H, rnorm):
    """Port of zneigh: eigenvalues of the (small, dense) Hessenberg
    matrix H and their Ritz estimates (= rnorm * |last eigenvector
    component|, from the Arnoldi relation Op*V - V*H = resid*e_k^T).
    numpy.linalg.eig already returns eigenvectors normalized to unit
    Euclidean norm, exactly what zneigh's own post-processing does by
    hand after zlahqr+ztrevc -- nothing extra needed here."""
    ritz, Y = np.linalg.eig(H)
    bounds = rnorm * np.abs(Y[-1, :])
    return ritz, bounds, Y


def apply_shifts(V, H, resid, shifts, kev):
    """Port of znapps (see module docstring for why dense QR replaces the
    bulge chase): compress a length len(V) Arnoldi factorization down to
    length kev by applying the given shifts."""
    kplusp = H.shape[0]
    Hc = H.copy()
    Q = np.eye(kplusp, dtype=complex)
    Ident = np.eye(kplusp, dtype=complex)
    for mu in shifts:
        Qi, _ = np.linalg.qr(Hc - mu * Ident)
        Hc = Qi.conj().T @ Hc @ Qi
        Q = Q @ Qi
    Vnew = []
    for i in range(kev + 1):
        wf = Q[0, i] * V[0]
        for j in range(1, kplusp):
            wf = wf + Q[j, i] * V[j]
        Vnew.append(wf)
    betak = Hc[kev, kev - 1]
    sigmak = Q[kplusp - 1, kev - 1]
    resid_new = sigmak * resid
    if abs(betak) > 1e-13:
        resid_new = resid_new + betak * Vnew[kev]
    Hnew = Hc[:kev, :kev].copy()
    rnorm_new = np.sqrt(abs(resid_new.dot(resid_new)))
    return Vnew[:kev], Hnew, resid_new, rnorm_new


def mpsiram(self, H, which="SR", nev=1, ncv=None, tol=1e-8, maxiter=200,
        verbose=0, wfskip=None):
    """Implicitly Restarted Arnoldi Method (port of znaupd/znaup2),
    matrix-free: finds the nev eigenpairs of H selected by `which`
    ({LM,SM,LR,SR,LI,SI}, ARPACK's own convention) using only Op(x)=H*x
    applications plus generic wavefunction algebra -- no auxiliary
    spectral-radius shift is needed the way arnolditk's mode="GS" uses
    one, since which="SR" already targets the smallest-real-part
    eigenvalues directly (exact shifts, not an ad hoc global one).

    Returns (es, wfs): es are <wf|H|wf> recomputed against the original
    (unshifted) H for accuracy/parity with arnolditk's own convention,
    wfs are normalized MPS/ED-state linear combinations of the final
    Krylov basis, ordered most-wanted-first.
    """
    if wfskip is None: wfskip = []
    if ncv is None: ncv = 2 * nev + 4
    ncv = max(ncv, nev + 2)
    M = self.toMPO(H, mode=arnolditk.arnoldimode)
    Op = lambda x: M * x

    eps23 = np.finfo(float).eps ** (2. / 3.)
    kev = nev
    resid = random_state(self, orthogonal=wfskip if wfskip else None)
    rnorm = 1.0
    V, Hm = [], np.zeros((0, 0), dtype=complex)
    V, Hm, resid, rnorm = arnoldi_extend(self, Op, V, Hm, resid, rnorm,
            0, kev, wfskip=wfskip, verbose=verbose)

    for it in range(maxiter):
        npcur = ncv - kev
        if npcur > 0:
            V, Hm, resid, rnorm = arnoldi_extend(self, Op, V, Hm, resid,
                    rnorm, kev, npcur, wfskip=wfskip, verbose=verbose)

        ritz, bounds, Y = hessenberg_ritz(Hm, rnorm)
        order = np.argsort(_goodness(which, ritz))
        ritz, bounds = ritz[order], bounds[order]

        np_apply = ncv - kev
        wanted_ritz = ritz[np_apply:]
        wanted_bounds = bounds[np_apply:]
        scale = np.maximum(eps23, np.abs(wanted_ritz))
        nconv = int(np.sum(wanted_bounds <= tol * scale))

        if verbose > 0:
            print("IRAM iter", it, "kev", kev, "nconv", nconv,
                    "ritz", np.round(wanted_ritz, 6))

        if nconv >= nev or np_apply == 0:
            break

        # stagnation-avoidance: temporarily grow kev, exactly as
        # znaup2's own heuristic (nev = nev + min(nconv,np/2)), so
        # repeatedly-shifted-away near-converged directions don't stall
        # the iteration
        if nconv > 0:
            bump = min(nconv, np_apply // 2)
            if bump > 0:
                kev = min(kev + bump, ncv - 2)

        np_apply = ncv - kev
        shifts = ritz[:np_apply]
        shift_bounds = bounds[:np_apply]
        # apply the most-uncertain (largest Ritz estimate) shifts first,
        # as zngets does, to limit forward instability of the sweep
        shift_order = np.argsort(-np.abs(shift_bounds))
        shifts = shifts[shift_order]

        V, Hm, resid, rnorm = apply_shifts(V, Hm, resid, shifts, kev)

    ritz, bounds, Y = hessenberg_ritz(Hm, rnorm)
    order = np.argsort(_goodness(which, ritz))
    ritz, Y = ritz[order], Y[:, order]
    nsel = min(nev, Hm.shape[0])
    Y_sel = Y[:, -nsel:]

    wfs = []
    for i in range(nsel):
        y = Y_sel[:, i]
        wf = y[0] * V[0]
        for j in range(1, len(V)):
            wf = wf + y[j] * V[j]
        wfs.append(wf.normalize())
    es = np.array([wf.aMb(H, wf) for wf in wfs])
    # most-wanted first (Y_sel's last column was the most wanted)
    return es[::-1], wfs[::-1]


def excited_states(self, H, nwf=1, which="SR", recursive=True,
        wfskip=None, **kwargs):
    """Find nwf eigenpairs of H, either recursively -- one at a time,
    each deflated against the previously found ones via wfskip -- or as
    a single simultaneous block IRAM search (nev=nwf in one call).
    Mirrors arnolditk's own recursive_arnoldi toggle: recursive is safer
    on (near-)degenerate spectra, since each state gets its own
    independently-seeded restart cycle rather than sharing one Krylov
    subspace across all nwf targets, at the cost of nwf independent
    restarts instead of one."""
    if wfskip is None: wfskip = []
    else: wfskip = list(wfskip)
    if nwf == 1:
        return mpsiram(self, H, which=which, nev=1, wfskip=wfskip, **kwargs)
    if not recursive:
        return mpsiram(self, H, which=which, nev=nwf, wfskip=wfskip, **kwargs)
    eout, wfout = [], []
    for i in range(nwf):
        es, wfs = mpsiram(self, H, which=which, nev=1, wfskip=wfskip, **kwargs)
        eout.append(es[0])
        wfout.append(wfs[0])
        wfskip.append(wfs[0])
    return np.array(eout), wfout


def lowest_energy_non_hermitian(self, H, n=1, **kwargs):
    """Smallest-real-part eigenpairs of a (possibly non-Hermitian) H via
    IRAM -- drop-in comparison point for
    arnolditk.lowest_energy_non_hermitian (mode="GS")."""
    return mpsiram(self, H, which="SR", nev=n, **kwargs)


def lowest_energy(self, H, n=1, **kwargs):
    """Ground state(s) of a Hermitian H via IRAM -- for a Hermitian H the
    smallest-real-part eigenvalue is the ground state, same convention as
    arnolditk.lowest_energy."""
    return mpsiram(self, H, which="SR", nev=n, **kwargs)
