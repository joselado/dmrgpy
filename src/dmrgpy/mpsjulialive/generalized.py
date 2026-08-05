"""Generalized-eigenvalue DMRG on the live Julia backend, i.e. the
smallest lambda solving H|psi>=lambda*A|psi> for a Hermitian
positive-definite metric operator A.

Thin wrapper over generalized.jl's get_gs_generalized (the whole
self-consistent Lagrange-multiplier iteration runs in one Julia call, see
that file's own header for the algorithm); the non-Hermitian counterpart
lives in this package's nhdmrg.py, driven by the same .jl file.
Dispatched to from groundstate.gs_energy_generalized's
itensor_version=="julia_live" branch, which is also where the shared
argument validation (mode="ED", Hermiticity, hamiltonian-is-set) happens.
"""

from contextlib import contextmanager

from .juliasession import Main as Mainjl
from .mpo import MPO
from .mps import MPS, random_mps


# Julia-level *programming* errors: a typo, a rename, a signature change
# in ITensorMPS/ITensorNHDMRG. These must never be translated into the
# retryable RuntimeError below -- nhdmrg.py's retry loops would swallow
# them ntries times and then report a physics-flavored cause for what is
# actually a broken .jl file.
_JULIA_BUG_MARKERS = ("UndefVarError", "MethodError", "UndefKeywordError",
        "TypeError", "ArgumentError", "syntax:", "LoadError")


@contextmanager
def metric_guard():
    """Re-raise solver-level Julia failures as the RuntimeError the other
    backends raise for the same conditions, so callers -- notably
    nhdmrg.py's retry loops, which catch RuntimeError to redraw a bad
    random start -- stay backend-agnostic.

    Deliberately a *denylist*, not an allowlist. An earlier version
    translated only messages containing "collapsed to ~0" (the .jl
    layer's own two guards), which left every other transient,
    redraw-able failure raised inside ITensorNHDMRG/KrylovKit -- a
    SingularException or LAPACKException out of the biorthogonalization
    step, say -- escaping the retry loop as a raw JuliaError and aborting
    a call that the very next random draw would have completed. Those are
    exactly the failures ntries exists for, and there is no way to
    enumerate their messages up front. So the default is now "retryable",
    with genuine .jl bugs (_JULIA_BUG_MARKERS) excluded so they keep
    propagating as themselves. nhdmrg.py's own final error carries the
    last attempt's message through, so a mistranslation here is
    diagnosable rather than silent.
    """
    from juliacall import JuliaError
    try:
        yield
    except JuliaError as e:
        text = str(e)
        if any(marker in text for marker in _JULIA_BUG_MARKERS):
            raise
        raise RuntimeError(text) from e


def gs_energy_generalized(self, A, lam0=None):
    """Return (lam,wf): the smallest generalized eigenvalue of
    (self.hamiltonian, A) and its eigenvector, both as ordinary
    Python-side objects (a float and an mpsjulialive MPS).

    A is a MultiOperator (already normalized by the caller via
    multioperator.obj2MO). lam0 is an optional starting lambda estimate;
    None is passed through to Julia as NaN, generalized.jl's own "unset"
    sentinel.
    """
    H = MPO(self.hamiltonian, MBO=self)
    Aop = MPO(A, MBO=self)
    psi0 = random_mps(self) # NOT self.random_state(), see julia_random_mps note
    with metric_guard():
        lam, psi = Mainjl.get_gs_generalized(H.jlmpo, Aop.jlmpo, psi0.jlmps,
                nsweeps=self.nsweeps, cutoff=self.cutoff, maxm=self.maxm,
                noise=self.noise,
                lam0=float("nan") if lam0 is None else float(lam0))
    return float(lam), MPS(psi, MBO=self)
