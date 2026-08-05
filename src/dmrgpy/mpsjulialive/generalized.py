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
from .mps import MPS


@contextmanager
def metric_guard():
    """Re-raise the .jl layer's own "collapsed to ~0" guards (ordinary
    Julia error()s, surfaced by juliacall as a JuliaError) as the
    RuntimeError the other backends raise for the same conditions, so
    callers -- notably nhdmrg.py's retry loops, which catch RuntimeError
    to redraw a bad random start -- stay backend-agnostic. Covers both
    generalized.jl's near-null-space metric guard and nhdmrg.jl's
    not-biorthogonalizable one.

    Matched on the guard's own message rather than translating every
    JuliaError: a genuine bug in the .jl code (UndefVarError, a method
    dispatch failure, ...) must keep propagating as itself, or that same
    retry loop would silently swallow it ntries times and then report
    "every attempt hit the near-null-space guard".
    """
    from juliacall import JuliaError
    try:
        yield
    except JuliaError as e:
        if "collapsed to ~0" in str(e):
            raise RuntimeError(str(e)) from e
        raise


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
    psi0 = self.random_state()
    with metric_guard():
        lam, psi = Mainjl.get_gs_generalized(H.jlmpo, Aop.jlmpo, psi0.jlmps,
                nsweeps=self.nsweeps, cutoff=self.cutoff, maxm=self.maxm,
                lam0=float("nan") if lam0 is None else float(lam0))
    return float(lam), MPS(psi, MBO=self)
