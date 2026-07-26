"""METTS-based finite-temperature vev -- a Many_Body_Chain-facing wrapper
around Chain.metts_vev() (pyitensor.chain.Chain for itensor_version=
"python", mpscpp3's Chain::metts_vev in-process pybind11 method for
itensor_version=3). See pyitensor/metts.py for the algorithm itself
(White & Stoudenmire, "Minimally entangled typical thermal state
algorithms", New J. Phys. 12, 055026 (2010), arXiv:1002.1305) --
mpscpp3/chain_session.h's own Chain::metts_vev is a direct port of that
same algorithm onto real ITensor v3 (imaginary-time TDVP + a sequential-
sampling MPS collapse using diagHermitian() for the per-site basis
rotation), not an independent implementation.

Only itensor_version in ("python", 3) are wired up -- itensor_version=2
(mpscpp2, no TDVP module at all) doesn't have a route to imaginary-time
evolution, so isn't supported. This is a genuinely different finite-
temperature method from the two already in this codebase: thermal.py's
Thermal_Spin_Chain (ancilla purification + imaginary-time annealing of a
single, growing-entanglement wavefunction) and vevtk/thermalvev.py's
thermal_vev_ex (exact Boltzmann sum over ED excited states) -- METTS
instead samples many *unentangled* typical states and averages, trading
exactness for a bond dimension that stays small even at low temperature
(no ancilla doubling).
"""

import random


def metts_vev(MB, Op, T, nsamples=200, nwarmup=20, dbeta_half_step=0.05,
              basis_ops=("Sz", "Sx"), seed=None, niter=30):
    """Thermal <Op> at temperature T via METTS sampling.

    MB: a Many_Body_Chain with itensor_version in ("python", 3) and a
        Hamiltonian already set via MB.set_hamiltonian(...).
    Op: a MultiOperator (the observable to measure).
    Returns (mean, stderr) -- see
    pyitensor.metts.metts_thermal_average's own docstring for the
    stderr's (Markov-correlated, so optimistic) meaning.
    """
    if MB.itensor_version not in ("python", 3):
        raise NotImplementedError(
            "metts_vev is only implemented for itensor_version='python' or "
            "3 so far (got itensor_version=%r)" % (MB.itensor_version,))
    if MB.itensor_version == 3 and MB.ns < 3:
        # Same underlying vendored-ITensor limitation CLAUDE.md documents
        # for v3's two-site dmrg() (mode.py's get_mode() falls back to ED
        # for it there): two-site TDVP hits the identical "LocalOp is
        # default constructed" SIGABRT for chains this short, confirmed
        # directly. metts_vev() doesn't go through get_mode() (it isn't a
        # DMRG-vs-ED choice), so it needs its own guard here rather than
        # silently inheriting that fallback -- and unlike gs_energy/vev,
        # there's no ED-based METTS to fall back to, so this raises
        # instead of substituting a different method.
        raise NotImplementedError(
            "metts_vev with itensor_version=3 can't handle a chain this "
            "short (n=%d < 3 sites): ITensor v3's two-site TDVP aborts the "
            "whole process for chains below 3 sites (the same vendored-"
            "ITensor limitation as its two-site dmrg(), see CLAUDE.md). "
            "Use itensor_version='python' instead for small chains." % (MB.ns,))
    # Validate here rather than relying on either backend's own checks:
    # mpscpp3/chain_session.h's Chain::metts_vev does validate T/basis_ops/
    # nsamples too, but via ITensor's Error(), which -- unlike a normal C++
    # exception -- calls abort() directly (see itensor/util/error.h) and so
    # can never be caught from Python; it's only a last-resort safety net
    # against direct _session.metts_vev() misuse, not a substitute for a
    # real, catchable check here (confirmed directly: a T<=0 call used to
    # SIGABRT the whole interpreter instead of raising before this check
    # was added).
    if T <= 0:
        raise ValueError("metts_vev: T must be > 0, got %r" % (T,))
    if not basis_ops:
        raise ValueError("metts_vev: basis_ops must be non-empty")
    if nsamples < 1:
        raise ValueError("metts_vev: nsamples must be >= 1, got %r" % (nsamples,))
    MB._session.set_sweep_params(MB.maxm, MB.nsweeps, MB.cutoff, MB.noise)
    MB._session.set_verbose(MB.verbose)
    MB._session.set_mpomaxm(max(MB.maxm, MB.mpomaxm))
    MB._session.set_hamiltonian(MB.hamiltonian.to_terms())
    # Both backends' metts_vev want a concrete integer seed (the C++ side
    # has no equivalent of Python's seed=None "fresh entropy" convention),
    # so draw one here when the caller didn't pin one -- keeps the
    # non-reproducible-by-default behavior every other seed= kwarg in
    # this codebase has, rather than silently becoming deterministic.
    seed_arg = seed if seed is not None else random.getrandbits(63)
    return MB._session.metts_vev(
        Op.to_terms(), T, nsamples, nwarmup,
        dbeta_half_step, list(basis_ops), seed_arg, niter)
