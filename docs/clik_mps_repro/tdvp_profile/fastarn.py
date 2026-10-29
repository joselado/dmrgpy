"""Prototype: Arnoldi exponentiator with a row-major basis (no conj copy of
the basis) and DGKS selective second Gram-Schmidt pass. Host/NumPy only."""
import numpy as np
from dmrgpy.pyitensor import tdvp
from dmrgpy.pyitensor.tdvp import (_hessenberg_exp_first_column,
        _arnoldi_converged_substep, _KRYLOV_ERRGOAL, _KRYLOV_MAX_SUBSTEPS)

ETA = 0.7071   # DGKS: reorthogonalize only if the norm fell below ETA of its pre-projection value
STATS = {"second_pass": 0, "iters": 0}
SELECTIVE = True

def _run(matvec, Q, coeff, m, n, tol, errgoal):
    Hm = np.zeros((m + 1, m), dtype=complex)
    scale = 0.0
    for it in range(m):
        w = matvec(Q[it])
        k = it + 1
        Qk = Q[:k]
        w0 = np.linalg.norm(w)
        h = (Qk @ w.conj()).conj()
        w = w - h @ Qk
        hnext = np.linalg.norm(w)
        STATS["iters"] += 1
        if (not SELECTIVE) or hnext < ETA * w0:
            h2 = (Qk @ w.conj()).conj()
            w = w - h2 @ Qk
            h = h + h2
            hnext = np.linalg.norm(w)
            STATS["second_pass"] += 1
        hnext = float(hnext)
        Hm[:k, it] = h
        Hm[k, it] = hnext
        scale = max(scale, float(np.max(np.abs(h))), hnext)
        Hk = Hm[:k, :k]
        c = _hessenberg_exp_first_column(Hk, coeff)
        if hnext <= tol * scale or k >= n:
            return k, Hk, hnext, True, c
        if abs(coeff) * hnext * abs(c[-1]) < errgoal:
            return k, Hk, hnext, True, c
        if k == m:
            return k, Hk, hnext, False, c
        Q[k] = w / hnext
    raise AssertionError

def fast_arnoldi_expm_multiply(matvec, v0, coeff, niter=30, tol=1e-12, errgoal=_KRYLOV_ERRGOAL):
    beta0 = float(np.linalg.norm(v0))
    if beta0 == 0:
        return v0.copy()
    n = v0.size; m = min(niter, n)
    v = v0; remaining = 1.0; frac = 1.0
    Q = np.empty((m, n), dtype=complex)
    for _ in range(_KRYLOV_MAX_SUBSTEPS):
        Q[0] = v / beta0
        k, Hk, hnext, converged, c = _run(matvec, Q, coeff * frac, m, n, tol, errgoal)
        if not converged:
            frac, c = _arnoldi_converged_substep(Hk, hnext, coeff, frac, errgoal)
        v = beta0 * (c @ Q[:k])
        if converged and frac >= remaining:
            return v
        remaining -= frac
        frac = min(frac, remaining)
        beta0 = float(np.linalg.norm(v))
        if beta0 == 0:
            return v
    raise RuntimeError("too many substeps")
