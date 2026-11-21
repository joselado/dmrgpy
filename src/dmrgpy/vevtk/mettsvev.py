"""METTS-based finite-temperature vev -- a Many_Body_Chain-facing wrapper
around Chain.metts_vev() (pyitensor.chain.Chain for itensor_version=
"python", mpscpp3's Chain::metts_vev in-process pybind11 method for
itensor_version=3, mpsjulialive.metts.metts_vev -> metts.jl's metts_vev
for itensor_version="julia_live"). See pyitensor/metts.py for the
algorithm itself (White & Stoudenmire, "Minimally entangled typical
thermal state algorithms", New J. Phys. 12, 055026 (2010),
arXiv:1002.1305) -- mpscpp3/chain_session.h's and mpsjulialive/metts.jl's
own metts_vev are each a direct port of that same algorithm (imaginary-
time TDVP + a sequential-sampling MPS collapse using per-site
eigendecomposition for the basis rotation), not independent
implementations.

Only itensor_version in ("python", 3, "julia_live") are wired up --
itensor_version=2 (mpscpp2, no TDVP module at all) doesn't have a route
to imaginary-time evolution, so isn't supported. This is a genuinely
different finite-temperature method from the two already in this
codebase: thermal.py's Thermal_Spin_Chain (ancilla purification +
imaginary-time annealing of a single, growing-entanglement wavefunction)
and vevtk/thermalvev.py's thermal_vev_ex (exact Boltzmann sum over ED
excited states) -- METTS instead samples many *unentangled* typical
states and averages, trading exactness for a bond dimension that stays
small even at low temperature (no ancilla doubling).
"""

import random

import numpy as np


_MATRIX_TOL = 1e-12 # relative, for the one-site matrix tests of _check_basis


def _site_types(MB):
    """The pyitensor site type of every site of MB, from its site-type
    codes, or None when a code has no table there. _check_basis reads its
    one-site matrices from these tables on every backend: a name is the
    same operator on v3 and julia_live, and every test made on the
    matrices (Hermiticity, commutation, squares) is unchanged by a common
    change of basis and by a transpose."""
    from ..pyitensor.sites.siteset import _site_class
    try: return [_site_class(code) for code in MB.sites]
    except (ValueError, TypeError): return None


def _matrix(t, name):
    """The matrix of a one-site name on site type t, or None if t does
    not know the name."""
    try: return np.asarray(t.matrix(name), dtype=complex)
    except (ValueError, KeyError, TypeError): return None


def _commutes(A, B, sign=-1.):
    """A B + sign B A vanishes, relative to ||A|| ||B|| (sign=+1 tests
    anticommutation)"""
    return np.linalg.norm(A @ B + sign*(B @ A)) <= \
            _MATRIX_TOL*np.linalg.norm(A)*np.linalg.norm(B)


def _is_scalar(q):
    """True if every one-site matrix of q is a multiple of the identity,
    so that its sum labels no sector"""
    return all(np.linalg.norm(m - np.trace(m)/len(m)*np.eye(len(m)))
               <= _MATRIX_TOL*np.linalg.norm(m) for m in q)


def _one_site_parts(types, Q):
    """[q_k] with Q = sum_k q_k, read off the one-factor terms of Q, or
    None when a term acts on more than one site or names an operator the
    site tables do not know"""
    q = [np.zeros((t.dim, t.dim), dtype=complex) for t in types]
    for term in Q.op:
        factors = [f for f in term[1:] if f[0] != "Id"]
        if term[0] == 0 or not factors: continue # a constant
        if len(factors) != 1: return None
        name, k = factors[0]
        if not 0 <= k < len(types): return None
        m = _matrix(types[k], name)
        if m is None: return None
        q[k] = q[k] + term[0]*m
    return q


def _parity_conserved(H, types, p):
    """True if every term of H commutes with P = prod_k p_k, decided
    factor by factor: a factor that anticommutes with p_k on its own site
    flips the sign of its term under P, one that commutes keeps it. A term
    with an odd number of flips, or with a factor that does neither, and
    the parity is not proven conserved. For the fermion parity the names
    are the raw C, Cdag, ..., whose Jordan-Wigner strings commute with P,
    so the local matrix of each factor is all that matters."""
    for term in H.op:
        if term[0] == 0: continue
        flips = 0
        for (name, k) in term[1:]:
            if name == "Id": continue
            if not 0 <= k < len(types): return False
            m = _matrix(types[k], name)
            if m is None: return False
            if _commutes(m, p[k]): continue
            if not _commutes(m, p[k], sign=1.): return False
            flips += 1
        if flips % 2: return False
    return True


def _check_basis(MB, basis_ops, caller="metts_vev"):
    """Refuse a collapse basis METTS cannot sample the thermal state with,
    on every backend and for both metts_vev and metts_dynamical_correlator
    (caller names the call in the message).

    Each name must be Hermitian, since the collapse diagonalizes it, and
    each backend reads a different triangle of the matrix. Then the
    Markov chain freezes when every basis name commutes, on each site,
    with the one-site part q_k of a charge that H conserves: the collapse
    then gives eigenstates of the charge and e^{-beta H/2} keeps them in
    their sector. On a 4-site hopping chain basis_ops=("N",) returned
    <Ntot> = 3, 1, 1, 4 by seed, with stderr 0, for an exact thermal
    1.7811 (2026-10-23 hole hunt, finding 1), and ("Sx",) on a Heisenberg
    chain froze the same way although Sx is not diagonal (2026-11-21 hole
    hunt, finding 7). The charges tested are:

    - the chain's conserved quantities (get_sector_charge_operators) and
      each basis name's own sum, Q = sum_k q_k: ("Nupdn",) and ("F",) on
      a Hubbard chain and ("N2",) on a Bose-Hubbard chain gave <Ntot> 6,
      4, 2 by seed (2026-10-23, leads, finding 10). [H,Q] = 0 is proven
      by the canonical form or tested on a random witness state, as
      mpsalgebra.is_hermitian does, since the canonical form cannot prove
      it when it rests on a same-site identity such as N C = 0;
    - for a name O with (2O)^2 = 1 on every site, the parity prod_k 2O_k:
      ("Sx",) on the transverse-field Ising chain with its field along x
      gave <Px> = -1 at stderr 0 for 0.038;
    - on a chain with fermion sites, the fermion parity (-1)^N. Every
      Hermitian one-site operator of a fermion site keeps it, so on a
      chain whose H does not conserve N every basis is frozen in one
      parity: ("N",) on a Kitaev chain gave <(-1)^N> = -1 at stderr 0 for
      0.043. This one is refused, by the author's pick.

    A parity is proven conserved term by term (_parity_conserved). The
    one-site matrices come from the pyitensor site tables, so a name they
    do not know is left to the backend, which raises on it."""
    from ..multioperator import obj2MO
    types = _site_types(MB)
    mats = {}
    for name in basis_ops:
        m = None if types is None else [_matrix(t, name) for t in types]
        if m is None or any(M is None for M in m):
            # not in the tables: the backend's own probe, which raises on
            # a name the backend does not know either
            if not MB.is_hermitian(obj2MO([[name, 0]])):
                raise ValueError("%s: basis operator %r is not Hermitian, "
                        "the collapse needs a Hermitian one" % (caller, name))
            mats = None
            continue
        if any(np.linalg.norm(M - M.conj().T) > _MATRIX_TOL*np.linalg.norm(M)
               for M in m):
            raise ValueError("%s: basis operator %r is not Hermitian, the "
                    "collapse needs a Hermitian one" % (caller, name))
        if mats is not None: mats[name] = m
    if mats is None: return
    from ..groundstate import require_hamiltonian
    require_hamiltonian(MB) # refused by name, as metts_plan would
    H = MB.hamiltonian
    cmax = max([abs(t[0]) for t in H.op]+[0.])
    if cmax == 0.: return
    H = H*(1./cmax)
    def frozen_by(q): # every basis name commutes with q on every site
        return all(_commutes(m[k], q[k]) for m in mats.values()
                   for k in range(len(q)))
    charges = []
    try: sector = sorted(MB.get_sector_charge_operators().items())
    except (NotImplementedError, AttributeError): sector = []
    for label, Q in sector:
        charges.append(("the charge %s" % label, Q, _one_site_parts(types, Q)))
    for name in basis_ops:
        charges.append(("sum_i %s_i" % name,
                        sum(obj2MO([[name, i]]) for i in range(MB.ns)),
                        mats[name]))
    wf = None # the random witness, drawn only when a charge needs it
    for label, Q, q in charges:
        if q is None or _is_scalar(q) or not frozen_by(q): continue
        C = H*Q - Q*H
        if not C.is_zero(): # not proven by the canonical form
            if wf is None: wf = MB.random_state()
            c = C*wf
            hq = (H*Q)*wf
            if c is None or hq is None: continue
            # written so that a norm that is not a number does not refuse
            if not abs(c.dot(c)) <= 1e-20*abs(hq.dot(hq)): continue
        raise ValueError("%s: H conserves %s and every basis operator "
                "commutes with it on each site, so the samples never leave "
                "the charge sector of the first one and the average is not "
                "thermal. Add a basis operator that changes the charge "
                "(e.g. basis_ops=(\"Sz\",\"Sx\") on spins), or use "
                "mode=\"ED\"." % (caller, label))
    for name in basis_ops:
        p = [2.*M for M in mats[name]]
        if not all(np.linalg.norm(pk @ pk - np.eye(len(pk)))
                   <= _MATRIX_TOL*len(pk) for pk in p): continue
        if not frozen_by(p) or not _parity_conserved(H, types, p): continue
        raise ValueError("%s: H conserves the parity prod_i 2%s_i and every "
                "basis operator commutes with it on each site, so the "
                "samples never leave the parity sector of the first one and "
                "the average is not thermal. Add a basis operator that "
                "anticommutes with it (e.g. basis_ops=(\"Sz\",\"Sx\") on "
                "spins), or use mode=\"ED\"." % (caller, name))
    if not any("Nf" in t.qn_names() for t in types): return
    p = [np.diag([(-1.)**c for c in t.charges("Nf")]).astype(complex)
         if "Nf" in t.qn_names() else np.eye(t.dim, dtype=complex)
         for t in types]
    if frozen_by(p) and _parity_conserved(H, types, p):
        raise ValueError("%s: H conserves the fermion parity (-1)^N and "
                "every basis operator keeps it on each site, so the samples "
                "never leave the parity sector of the first one and the "
                "average is not thermal. METTS needs a collapse basis that "
                "mixes the parity, which the fermion sites do not offer: "
                "write the model as a spin chain through the Jordan-Wigner "
                "transformation and take basis_ops=(\"Sz\",\"Sx\"), or use "
                "mode=\"ED\"." % (caller,))


def metts_vev(MB, Op, T, nsamples=200, nwarmup=20, dbeta_half_step=0.05,
              basis_ops=("Sz", "Sx"), seed=None, niter=30, njobs=1,
              itime=None, itime_options=None):
    """Thermal <Op> at temperature T via METTS sampling.

    MB: a Many_Body_Chain with itensor_version in ("python", 3, "julia_live")
        and a Hamiltonian already set via MB.set_hamiltonian(...).
    Op: a MultiOperator (the observable to measure), or a list/tuple of
        MultiOperators to measure together on one shared METTS sample
        chain. The (nwarmup+nsamples) imaginary-time evolutions are the
        expensive part of this method, and every backend's metts_vev
        already samples every requested operator off the same chain of
        METTS states, so measuring several operators this way is far
        cheaper than calling metts_vev once per operator.
    njobs: run njobs independent METTS Markov chains in parallel worker
        processes and pool their statistics -- itensor_version="python"
        only (see pyitensor.metts.metts_thermal_average's own docstring
        for why, and for the exact pooling this uses). Measured directly
        (6-site Heisenberg+field chain, nsamples=200, nwarmup=20): wall
        time went 34s/21s/17s/14s for njobs=1/2/4/8 -- real but sublinear
        speedup (diminishing returns, not a bug), since each additional
        chain repeats the *full* nwarmup rather than sharing it, and each
        spawned worker pays its own fresh-interpreter startup cost (see
        below on why 'spawn' rather than fork). Worth it whenever
        nsamples/njobs stays comfortably above nwarmup. Unlike
        itensor_version=3 or "julia_live", each a single, live in-process
        session with no multiprocessing story (no per-worker copy to hand
        out), and unlike naively raising nsamples in a single chain,
        which buys no parallelism at all. There is no equivalent knob for
        itensor_version=3/"julia_live" here; see CLAUDE.md/kernels.py's own findings
        on why the numba/JAX contraction-kernel route (the other lever
        available for this pure-Python engine) makes METTS *slower*, not
        faster: every sample restarts from a fresh bond-dimension-1
        product state, so a run exercises many distinct contraction
        shapes rather than the same one repeatedly, and the fixed
        per-shape compile tax never gets amortized within a realistic
        sample count -- confirmed directly, ~4x-24x slower with
        kernels.USE_NUMBA=True depending on chain length, not faster.
    itime: how each sample is evolved to beta/2. "tdvp" takes steps of
        dbeta_half_step by two-site TDVP; "implicit" takes the implicit
        steps of a logarithmic grid (imaginarytime.py, arXiv:2606.02930),
        whose number grows as log(beta), with itime_options the keywords
        of evolve_imaginary_time (nsteps, default 8 here, tau0, order,
        shift, nsweeps, tol, krylov, maxiter). None (the default) takes
        whichever needs fewer steps at the accuracy of the TDVP default,
        see imaginarytime.metts_plan: TDVP at high T, in a conserved
        sector and for a non-Hermitian H, the implicit steps below T of
        about 0.12 in units of the largest coupling.
    Returns (mean, stderr) for a single Op, or a list of (mean, stderr)
    pairs -- one per operator, same order as given -- for a list/tuple of
    Ops. See pyitensor.metts.metts_thermal_average's own docstring for
    the stderr's (Markov-correlated, so optimistic) meaning.
    """
    if getattr(MB, "conserved_sector", None):
        # the collapse draws product states from every sector at once, as
        # pyitensor's Chain._reject_sector says; refused here, ahead of the
        # basis check, which on v3 did not return in 3 minutes in a sector
        # (2026-10-23 guide check)
        raise ValueError(
            "conserved-sector mode (%s): metts_vev is not available inside a "
            "sector, since the METTS collapse changes the charge; clear it "
            "with set_conserved_sector()" % (MB.conserved_sector,))
    if MB.itensor_version not in ("python", 3, "julia_live"):
        raise NotImplementedError(
            "metts_vev is only implemented for itensor_version='python', "
            "3, or 'julia_live' so far (got itensor_version=%r)" % (MB.itensor_version,))
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
    if njobs < 1:
        raise ValueError("metts_vev: njobs must be >= 1, got %r" % (njobs,))
    if njobs > 1 and MB.itensor_version != "python":
        # mpscpp3's Chain and mpsjulialive's live Julia session are each a
        # single live in-process session (see CLAUDE.md's cpp_handle
        # section for the C++ side; the Julia side has no per-worker copy
        # either, same reasoning) -- there's no per-worker copy a
        # multiprocessing pool could hand out the way pyitensor's Chain.H/
        # sites (plain Python/numpy objects, confirmed picklable) allow.
        # Rebuilding a whole separate v3/Julia session per worker process
        # is a bigger undertaking than this parameter is meant to cover,
        # so this raises rather than silently ignoring njobs or silently
        # falling back to a single sequential chain.
        raise NotImplementedError(
            "metts_vev: njobs>1 (parallel independent Markov chains) is "
            "only implemented for itensor_version='python' so far (got "
            "itensor_version=%r)" % (MB.itensor_version,))
    _check_basis(MB, basis_ops)
    single = not isinstance(Op, (list, tuple))
    ops = [Op] if single else list(Op)
    if not ops:
        raise ValueError("metts_vev: Op list/tuple must be non-empty")
    # Every backend's metts_vev wants a concrete integer seed (no
    # equivalent of Python's seed=None "fresh entropy" convention), so
    # draw one here when the caller didn't pin one -- keeps the
    # non-reproducible-by-default behavior every other seed= kwarg in
    # this codebase has, rather than silently becoming deterministic.
    seed_arg = seed if seed is not None else random.getrandbits(63)
    from ..imaginarytime import metts_plan, session_plan_kwargs
    plan = metts_plan(MB, T, itime, itime_options, dbeta_half_step)
    if MB.itensor_version == "julia_live":
        # This backend has no MB._session at all (its state lives in the
        # live Julia session instead, via MB.hamiltonian/MB.jlsites) --
        # mirrors groundstate.py's gs_energy_generalized julia_live branch.
        from ..mpsjulialive.metts import metts_vev as metts_vev_jl
        means, stderrs = metts_vev_jl(MB, ops, T, nsamples, nwarmup,
                dbeta_half_step, basis_ops, seed_arg, niter, plan=plan)
        if single:
            return means[0], stderrs[0]
        return list(zip(means, stderrs))
    MB._session.set_sweep_params(MB.maxm, MB.nsweeps, MB.cutoff, MB.noise)
    MB._session.set_verbose(MB.verbose)
    MB._session.set_mpomaxm(max(MB.maxm, MB.mpomaxm))
    MB._session.set_hamiltonian(MB.hamiltonian.to_terms())
    terms_ops = [op.to_terms() for op in ops]
    if MB.itensor_version == "python":
        means, stderrs = MB._session.metts_vev(
            terms_ops, T, nsamples, nwarmup,
            dbeta_half_step, list(basis_ops), seed_arg, niter, njobs=njobs,
            **session_plan_kwargs(MB, plan))
    else:
        means, stderrs = MB._session.metts_vev(
            terms_ops, T, nsamples, nwarmup,
            dbeta_half_step, list(basis_ops), seed_arg, niter,
            **session_plan_kwargs(MB, plan))
    if single:
        return means[0], stderrs[0]
    return list(zip(means, stderrs))
