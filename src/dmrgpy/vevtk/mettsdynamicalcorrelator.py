"""METTS-based finite-temperature dynamical correlator -- a
Many_Body_Chain-facing wrapper around Chain.metts_dynamical_correlator()
(pyitensor.chain.Chain for itensor_version="python", mpscpp3's
Chain::metts_dynamical_correlator in-process pybind11 method for
itensor_version=3, mpsjulialive.metts.metts_dynamical_correlator ->
metts.jl's metts_dynamical_correlator for itensor_version="julia_live").

Reference: Z. Wang, P. McClarty, D. Dankova, A. Honecker, A. Wietek,
"Spectroscopy and complex-time correlations using minimally entangled
typical thermal states", arXiv:2405.18484 -- Sec. II ("Dynamical METTS
algorithm"), Algorithm 1. Extends the static thermal_vev/metts_vev
machinery (arXiv:1002.1305, see mettsvev.py) from thermal expectation
values <A>_beta to the real-time two-time correlator
C_AB(t) = <A(t)B>_beta = <e^{iHt} A e^{-iHt} B>_beta: for every METTS
sample |psi_i> sampled by the same Markov chain metts_vev already uses,
define |v_i(0)>=B|psi_i>, |w_i(0)>=|psi_i>, evolve both in real time
under H, and average C^i(t)=<w_i(t)|A|v_i(t)> over samples -- no
importance reweighting needed, for exactly the same reason the static
algorithm's plain sample average already converges to the thermal
average (see metts.py's module docstring).

Only itensor_version in ("python", 3, "julia_live") are wired up, same
restriction as metts_vev (itensor_version=2/mpscpp2 has no TDVP module at
all, so no route to real- or imaginary-time evolution).
"""

import random

import numpy as np

from .. import operatornames


def metts_dynamical_correlator(MB, name, T, nt=200, dt=0.1, nsamples=100,
                                nwarmup=20, dbeta_half_step=0.05,
                                basis_ops=("Sz", "Sx"), seed=None, niter=30,
                                tdvp_cutoff=None, tdvp_maxdim=None,
                                tdvp_niter=50, njobs=1):
    """Real-time dynamical correlator C_AB(t) = <A(t)B>_T at temperature T
    via dynamical METTS sampling (arXiv:2405.18484, Sec. II).

    MB: a Many_Body_Chain with itensor_version in ("python", 3,
        "julia_live") and a Hamiltonian already set via
        MB.set_hamiltonian(...).
    name: the (A,B) operator pair, in the same format
        operatornames.str2MO() accepts elsewhere in this codebase -- a
        string like "ZZ" (site indices given via i=/j= kwargs), or an
        explicit (MultiOperator,MultiOperator) tuple/list. A and B are
        used exactly as resolved (no dagger applied to either), matching
        edtk/dynamics.py's dynamical_correlator_ED/dynamical_correlator_kpm
        convention -- this method is meant to be validated directly
        against MB.get_dynamical_correlator(mode="ED", submode="ED",
        T=..., name=name) (see edtk/dynamics.py's
        dynamical_correlator_finite_T), which shares that exact same
        convention.
    T: temperature (>0).
    nt, dt: number of real-time measurements and the (uniform) time step
        between them, i.e. C(t) is sampled at t=0,dt,...,(nt-1)*dt --
        same convention as Many_Body_Chain.evolution_DC's own (nt,dt).
    nsamples, nwarmup, dbeta_half_step, basis_ops, seed, niter: as in
        MB.metts_vev() -- control the shared METTS Markov chain
        (imaginary-time sampling) each retained sample is measured from.
    tdvp_cutoff, tdvp_maxdim: truncation controls for the *real-time*
        evolution of |v_i(t)>/|w_i(t)> specifically, separate from this
        chain's own cutoff/maxm (used for the imaginary-time METTS
        sampling step) -- default None, meaning "same as this chain's own
        cutoff/maxm" (see pyitensor.metts.metts_dynamical_correlator's own
        docstring on why real-time evolution generally wants a looser
        cutoff/larger bond dimension). itensor_version in ("python",
        "julia_live") only -- v3's Chain::metts_dynamical_correlator has
        no such knob (always uses this chain's own cutoff/maxm for both
        imaginary- and real-time evolution, like every other v3 TDVP
        method); passing either on itensor_version=3 raises rather than
        silently ignoring it. Julia's tdvp_step (mpsjulialive/tdvp.jl)
        takes cutoff/maxdim per call like pyitensor's does, so this is
        wired through there too rather than restricted the way v3 is.
    tdvp_niter: Krylov-iteration bound for the *real-time* TDVP evolution
        of |v_i(t)>/|w_i(t)> specifically (separate from `niter`, which
        only controls the imaginary-time sampling step) -- see
        pyitensor.metts.metts_dynamical_correlator's own docstring on why
        real- and imaginary-time evolution get independent controls here.
        Accepted but silently ignored for itensor_version="julia_live",
        same reason MB.metts_vev()'s own `niter` is: ITensorMPS.jl's
        tdvp() manages its own internal Krylov dimension with no exposed
        per-step iteration-count knob.
    njobs: run njobs independent METTS Markov chains in parallel worker
        processes and pool their statistics -- itensor_version="python"
        only, same mechanism/caveats as MB.metts_vev(njobs=...).

    Returns (ts, means, stderrs): ts is the array of nt measurement
    times, means[k] the sample-averaged C_AB(ts[k]) (complex), stderrs[k]
    its naive iid standard error over nsamples retained METTS samples --
    see MB.metts_vev()'s own docstring for the same Markov-correlation
    caveat on what that stderr does and doesn't capture.
    """
    if MB.itensor_version not in ("python", 3, "julia_live"):
        raise NotImplementedError(
            "metts_dynamical_correlator is only implemented for "
            "itensor_version='python', 3 or 'julia_live' so far (got "
            "itensor_version=%r)" % (MB.itensor_version,))
    if MB.itensor_version == 3 and MB.ns < 3:
        # Same underlying vendored-ITensor limitation as metts_vev's own
        # guard (see vevtk/mettsvev.py) -- two-site TDVP SIGABRTs for
        # chains this short, confirmed directly, and there's no ED-based
        # dynamical METTS to fall back to.
        raise NotImplementedError(
            "metts_dynamical_correlator with itensor_version=3 can't "
            "handle a chain this short (n=%d < 3 sites): ITensor v3's "
            "two-site TDVP aborts the whole process for chains below 3 "
            "sites (the same vendored-ITensor limitation as its two-site "
            "dmrg(), see CLAUDE.md). Use itensor_version='python' instead "
            "for small chains." % (MB.ns,))
    if T <= 0:
        raise ValueError("metts_dynamical_correlator: T must be > 0, got %r" % (T,))
    if not basis_ops:
        raise ValueError("metts_dynamical_correlator: basis_ops must be non-empty")
    if nsamples < 1:
        raise ValueError("metts_dynamical_correlator: nsamples must be >= 1, got %r" % (nsamples,))
    if nt < 1:
        raise ValueError("metts_dynamical_correlator: nt must be >= 1, got %r" % (nt,))
    if njobs < 1:
        raise ValueError("metts_dynamical_correlator: njobs must be >= 1, got %r" % (njobs,))
    if njobs > 1 and MB.itensor_version != "python":
        # Same reason as metts_vev's own njobs>1 guard: mpscpp3's Chain is
        # a single live in-process C++ session with no per-worker copy a
        # multiprocessing pool could hand out.
        raise NotImplementedError(
            "metts_dynamical_correlator: njobs>1 (parallel independent "
            "Markov chains) is only implemented for itensor_version="
            "'python' so far (got itensor_version=%r)" % (MB.itensor_version,))
    if MB.itensor_version not in ("python", "julia_live") and (tdvp_cutoff is not None or tdvp_maxdim is not None):
        # v3's Chain::metts_dynamical_correlator has no separate real-time
        # cutoff/maxdim knob at all -- raise rather than silently ignoring
        # a value the caller explicitly asked for.
        raise NotImplementedError(
            "metts_dynamical_correlator: tdvp_cutoff/tdvp_maxdim are only "
            "implemented for itensor_version='python' or 'julia_live' so "
            "far (got itensor_version=%r)" % (MB.itensor_version,))

    resolved = operatornames.str2MO(MB, name,
            require_symbolic_for="the METTS dynamical correlator")
    A, B = resolved[0], resolved[1]

    # Both backends' metts_dynamical_correlator want a concrete integer
    # seed, same as metts_vev -- draw one here when the caller didn't pin
    # one, rather than silently becoming deterministic.
    seed_arg = seed if seed is not None else random.getrandbits(63)

    if MB.itensor_version == "julia_live":
        # This backend has no MB._session at all (its state lives in the
        # live Julia session instead, via MB.hamiltonian/MB.jlsites) --
        # mirrors mettsvev.py's own julia_live branch.
        from ..mpsjulialive.metts import metts_dynamical_correlator as metts_dc_jl
        means, stderrs = metts_dc_jl(MB, A, B, T, nt, dt, nsamples, nwarmup,
                dbeta_half_step, basis_ops, seed_arg, niter,
                tdvp_cutoff, tdvp_maxdim, tdvp_niter)
        ts = np.array([dt * k for k in range(nt)])
        return ts, np.array(means), np.array(stderrs)

    MB._session.set_sweep_params(MB.maxm, MB.nsweeps, MB.cutoff, MB.noise)
    MB._session.set_verbose(MB.verbose)
    MB._session.set_mpomaxm(max(MB.maxm, MB.mpomaxm))
    MB._session.set_hamiltonian(MB.hamiltonian.to_terms())
    terms_a, terms_b = A.to_terms(), B.to_terms()
    if MB.itensor_version == "python":
        means, stderrs = MB._session.metts_dynamical_correlator(
            terms_a, terms_b, T, nt, dt, nsamples, nwarmup, dbeta_half_step,
            list(basis_ops), seed_arg, niter, tdvp_cutoff=tdvp_cutoff,
            tdvp_maxdim=tdvp_maxdim, tdvp_niter=tdvp_niter, njobs=njobs)
    else:
        means, stderrs = MB._session.metts_dynamical_correlator(
            terms_a, terms_b, T, nt, dt, nsamples, nwarmup, dbeta_half_step,
            list(basis_ops), seed_arg, niter, tdvp_niter)
    ts = np.array([dt * k for k in range(nt)])
    return ts, np.array(means), np.array(stderrs)
