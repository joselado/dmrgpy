"""METTS-based finite-temperature vev on the live Julia backend.

Thin wrapper over metts.jl's metts_vev (the whole nwarmup+nsamples Markov
chain runs in one Julia call, see that file's own header for the
algorithm); mirrors generalized.py's shape (build MPOs, hand them to one
Julia entry point, unpack the result back into plain Python objects).
Dispatched to from vevtk/mettsvev.py's itensor_version=="julia_live"
branch, which is also where the shared argument validation (T>0,
basis_ops non-empty, njobs==1 only, ...) happens.
"""

from .juliasession import Main as Mainjl
from .juliasession import to_julia_strvec
from .mpo import MPO


def metts_vev(self, ops, T, nsamples, nwarmup, dbeta_half_step, basis_ops,
              seed, niter):
    """Thermal <op> for every op in `ops` (a list of MultiOperators), via
    METTS sampling on the live Julia backend.

    `niter` is accepted only for signature parity with the pyitensor/v3
    backends (see vevtk/mettsvev.py) -- ITensorMPS.jl's tdvp() manages its
    own internal Krylov dimension, with no direct per-step iteration-count
    knob exposed anywhere else on this backend (evolve_and_measure_tdvp/
    quench_tdvp in timedependent.py don't take one either), so it is
    silently ignored here rather than threaded through to Julia.

    Returns (means, stderrs) -- means[k]/stderrs[k] as complex/real
    Python floats, one entry per op in `ops`, same convention as the
    other two backends' Chain.metts_vev.
    """
    H = MPO(self.hamiltonian, MBO=self)
    opmpos = [MPO(op, MBO=self).jlmpo for op in ops]
    basis_ops_jl = to_julia_strvec(list(basis_ops))
    means, stderrs = Mainjl.metts_vev(H.jlmpo, self.jlsites, opmpos, T,
            int(nsamples), int(nwarmup), dbeta_half_step, basis_ops_jl,
            int(seed), self.cutoff, self.maxm)
    means = [complex(m) for m in means]
    stderrs = [float(s) for s in stderrs]
    return means, stderrs
