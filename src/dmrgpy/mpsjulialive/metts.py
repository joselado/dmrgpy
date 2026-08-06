"""METTS-based finite-temperature vev and dynamical correlator on the live
Julia backend.

Thin wrappers over metts.jl's metts_vev/metts_dynamical_correlator (each
runs its whole nwarmup+nsamples Markov chain -- including, for the
dynamical correlator, the per-sample real-time double evolution -- in one
Julia call, see that file's own header for both algorithms); mirrors
generalized.py's shape (build MPOs, hand them to one Julia entry point,
unpack the result back into plain Python objects). Dispatched to from
vevtk/mettsvev.py's and vevtk/mettsdynamicalcorrelator.py's
itensor_version=="julia_live" branches, which are also where the shared
argument validation (T>0, basis_ops non-empty, njobs==1 only, ...)
happens.
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


def metts_dynamical_correlator(self, A_op, B_op, T, nt, dt, nsamples, nwarmup,
                                dbeta_half_step, basis_ops, seed, niter,
                                tdvp_cutoff, tdvp_maxdim, tdvp_niter):
    """Real-time C_AB(t)=<A(t)B>_T for the (A_op,B_op) MultiOperator pair,
    via dynamical METTS sampling on the live Julia backend -- thin
    wrapper over metts.jl's metts_dynamical_correlator (the whole
    nwarmup+nsamples Markov chain, including the per-sample real-time
    double evolution, runs in one Julia call, see that function's own
    header).

    `niter`/`tdvp_niter` are accepted only for signature parity with the
    pyitensor/v3 backends (see vevtk/mettsdynamicalcorrelator.py) --
    ITensorMPS.jl's tdvp() manages its own internal Krylov dimension for
    both the imaginary- and real-time evolution steps, with no per-step
    iteration-count knob exposed anywhere else on this backend, so both
    are silently ignored here rather than threaded through to Julia (same
    convention as metts_vev's own `niter` above).

    tdvp_cutoff/tdvp_maxdim: truncation controls for the *real-time*
    evolution of v(t)/w(t) specifically -- default (None) to this chain's
    own cutoff/maxm, same convention as the pyitensor backend.

    Returns (means, stderrs) -- length-nt complex/real Python lists, same
    convention as the pyitensor/v3 backends' Chain.metts_dynamical_correlator.
    """
    if tdvp_cutoff is None:
        tdvp_cutoff = self.cutoff
    if tdvp_maxdim is None:
        tdvp_maxdim = self.maxm
    H = MPO(self.hamiltonian, MBO=self)
    Ampo = MPO(A_op, MBO=self).jlmpo
    Bmpo = MPO(B_op, MBO=self).jlmpo
    basis_ops_jl = to_julia_strvec(list(basis_ops))
    means, stderrs = Mainjl.metts_dynamical_correlator(
            H.jlmpo, self.jlsites, Ampo, Bmpo, T, int(nt), dt,
            int(nsamples), int(nwarmup), dbeta_half_step, basis_ops_jl,
            int(seed), self.cutoff, self.maxm, tdvp_cutoff, tdvp_maxdim)
    means = [complex(m) for m in means]
    stderrs = [float(s) for s in stderrs]
    return means, stderrs
