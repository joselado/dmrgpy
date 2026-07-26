"""METTS: Minimally Entangled Typical Thermal States.

Reference: E.M. Stoudenmire and Steven R. White, "Minimally entangled
typical thermal state algorithms", New J. Phys. 12, 055026 (2010),
arXiv:1002.1305.

Algorithm (paper's Sec. 2):

  1. Start from a classical product state (CPS) |i>: an unentangled
     product of single-site basis states (not necessarily the
     computational/Sz basis -- any product of single-site vectors
     counts, e.g. |+x> at some sites).
  2. Imaginary-time evolve it by half the inverse temperature:
         |phi_i> = e^{-beta H/2} |i> / || e^{-beta H/2} |i> ||.
     Every operator expectation value <A> is estimated as the average of
     <phi_i|A|phi_i> over many samples i, because
         Tr[A e^{-beta H}] / Tr[e^{-beta H}]
           = sum_i <i|e^{-beta H/2} A e^{-beta H/2}|i>
             / sum_i <i|e^{-beta H}|i>
     for *any* complete CPS basis {|i>} (Eq. (1)-(2) of the paper), and
     the specific Markov chain built by repeatedly collapsing |phi_i>
     back down to a new CPS |i'> with probability |<i'|phi_i>|^2 has
     exactly this Boltzmann-weighted ensemble as its stationary
     distribution (Eq. (3)-(5)): each |phi_i> sampled from that chain is
     already properly weighted, so a plain (unweighted) sample average
     of <phi_i|A|phi_i> converges to the thermal average -- no importance
     reweighting needed, unlike generic Monte Carlo.
  3. Collapse |phi_i> to a new CPS by sequential ("perfect") sampling: a
     single left-to-right sweep samples one site at a time from its exact
     conditional marginal (given the sites already sampled to its left),
     using only the local MPS tensor -- never the full 2^N amplitude
     vector -- because canonicalizing to site 1 makes every site to the
     right already right-orthonormal, i.e. the partial trace over "the
     rest of the chain" needed for that marginal is already built into
     the tensor's own gauge (see collapse_to_cps() below; this is the
     same "perfect sampling from an MPS" trick as Ferris & Vidal, PRB 85,
     165146 (2012), and ITensor's own MPS::sample()).
  4. Alternate the collapse basis between two mutually unbiased choices
     (e.g. Sz and Sx for spin-1/2) from one sample to the next: collapsing
     repeatedly in the *same* basis lets the Markov chain get trapped for
     many steps in that basis's symmetry sector (the paper's Sec. 3.3),
     inflating the autocorrelation time; alternating fixes this in
     practice for the models the paper tests.

This module works with plain pyitensor MPS/MPO objects (mpscontainer.py)
and tdvp.py's TDVP stepper -- imaginary time is just tdvp_step()'s own
complex-dt generalization with dt=-1j*dtau (dt purely imaginary <->
tdvp_step's usual real-time exp(-i*dt*H); dt purely real <-> the decaying
exp(-dtau*H) imaginary-time propagator wanted here, since forward/backward
half-sweep coefficients -1j*(dt/2) and +1j*(dt/2) become -dtau/2 and
+dtau/2 respectively for dt=-1j*dtau -- a real, negative forward exponent,
i.e. genuine decay, not rotation). See chain.py's Chain.metts_vev() for how
this is wired into the rest of the engine (H and the observables' MPOs are
each built once from MultiOperator terms and reused across every sample).
"""

import multiprocessing

import numpy as np

from .index import Index, reseed_id_counter_past
from .mpscontainer import MPS
from .mpsalgebra import inner
from .tdvp import tdvp_step
from .tensor import ITensor, noPrime


def _site_operator_matrix(sites, opname, i):
    """The (dim,dim) matrix of single-site operator `opname` at (1-based)
    site i, in the convention that plain matrix-vector multiplication
    reproduces exactly what sites.op(opname,i) does when contracted into
    an MPS tensor -- obtained by actually invoking that same contraction
    on each basis vector, rather than re-deriving sites/base.py's
    unprimed-in/primed-out convention algebraically here."""
    s = sites.si(i)
    d = s.dim
    op = sites.op(opname, i)
    out = np.zeros((d, d), dtype=complex)
    for a in range(d):
        vec = np.zeros(d, dtype=complex)
        vec[a] = 1.0
        ket = ITensor((s,), vec)
        res = noPrime(op * ket, "Site")
        out[:, a] = res.array
    return out


def _eigenbasis(sites, opname, i):
    """(eigvals, eigvecs) of the named single-site Hermitian operator at
    site i, eigvecs as columns expressed in the site's native
    (computational) basis -- np.linalg.eigh, since every physical
    operator dmrgpy's site tables define (Sz, Sx, N, ...) is Hermitian."""
    M = _site_operator_matrix(sites, opname, i)
    return np.linalg.eigh(M)


def _build_eigcache(sites, basis_ops):
    """{(opname, i): (evals, evecs)} for every (opname, site) pair the
    sampling run will ever need, computed once up front. Without this,
    random_cps()/collapse_to_cps() would recompute the same operator's
    eigendecomposition from scratch at every site on every one of the
    nwarmup+nsamples Markov-chain iterations, even though a given
    (opname, i)'s eigenbasis never changes across the whole run -- with
    n sites and O(200-300) iterations that's tens of thousands of
    redundant eigh() calls for no reason."""
    n = sites.length()
    return {(opname, i): _eigenbasis(sites, opname, i)
            for opname in set(basis_ops) for i in range(1, n + 1)}


def product_state(sites, vectors):
    """An MPS product state (bond dimension 1 throughout) from a list of
    per-site complex vectors (vectors[i-1], length sites.dim(i)). Trivial
    dim-1 Link indices are still built between neighboring sites so this
    looks, to every downstream routine (mpsalgebra.inner, MPS.position,
    tdvp_step, ...), like any other MPS in this engine -- see
    mpscontainer.py's own module docstring on why interior sites always
    carry a Link-tagged index here."""
    n = sites.length()
    links = [Index(1, tags="Link,l={}".format(k + 1)) for k in range(n - 1)]
    tensors = []
    for i in range(1, n + 1):
        s = sites.si(i)
        vec = np.asarray(vectors[i - 1], dtype=complex).reshape(-1)
        inds = []
        if i > 1:
            inds.append(links[i - 2])
        inds.append(s)
        if i < n:
            inds.append(links[i - 1])
        shape = tuple(ind.dim for ind in inds)
        tensors.append(ITensor(tuple(inds), vec.reshape(shape)))
    mps = MPS(tensors)
    mps.center = 1
    return mps


def random_cps(sites, opname, rng, eigcache=None):
    """A random classical product state: at each site, uniformly pick one
    eigenstate of the named single-site operator (e.g. "Sz"). Valid seed
    for the METTS Markov chain regardless of choice -- the chain's
    stationary distribution doesn't depend on the starting CPS (paper,
    Sec. 2), only how many warmup steps are needed to reach it.

    eigcache: optional {(opname, i): (evals, evecs)} from
    _build_eigcache(), used instead of recomputing on every call."""
    n = sites.length()
    vectors = []
    for i in range(1, n + 1):
        _, evecs = eigcache[(opname, i)] if eigcache is not None else _eigenbasis(sites, opname, i)
        k = rng.integers(0, sites.dim(i))
        vectors.append(evecs[:, k].copy())
    return product_state(sites, vectors)


def collapse_to_cps(psi, sites, opname, rng, eigcache=None):
    """Sequential ("perfect") sampling collapse of MPS `psi` onto a new
    classical product state built from eigenstates of the single-site
    Hermitian operator `opname` (e.g. "Sz" or "Sx") -- the paper's
    "quantum measurement" collapse step (Sec. 2.2), generalized from a
    single qubit measurement to the whole chain via one left-to-right
    sweep.

    Canonicalizing to site 1 first (work.position(1)) makes every site
    i>1 right-orthonormal; contracting the running "collapsed-so-far"
    vector L into site i's tensor and reading off
    sum_{right link} |<eigenstate k|tensor>|^2 then gives exactly the
    marginal probability of outcome k at site i *conditioned on* the
    outcomes already sampled to its left -- no explicit partial trace is
    needed because right-orthonormality already performs it. This means
    the full N-site joint distribution is sampled exactly (up to the
    original MPS's own truncation), never requiring the full
    exponentially-sized amplitude vector.

    eigcache: optional {(opname, i): (evals, evecs)} from
    _build_eigcache(), used instead of recomputing on every call.

    Returns (new_cps_mps, outcomes): outcomes[i-1] is the sampled
    eigenvalue (float) at site i, for diagnostics only.
    """
    n = sites.length()
    work = psi.copy()
    work.position(1)
    vectors = []
    outcomes = []
    L = None  # running collapsed-prefix amplitude, an ITensor over the current left link, or None at site 1
    for i in range(1, n + 1):
        evals, evecs = eigcache[(opname, i)] if eigcache is not None else _eigenbasis(sites, opname, i)
        T = work.A(i)
        if L is not None:
            T = L * T
        s = sites.si(i)
        right_link = next((ind for ind in T.inds if ind != s), None)
        if right_link is not None:
            mat = T.transpose_to((s, right_link))
        else:
            mat = T.transpose_to((s,)).reshape(s.dim, 1)
        rot = evecs.conj().T @ mat  # amplitude of each eigenbasis outcome, per remaining right-link value
        probs_raw = np.sum(np.abs(rot) ** 2, axis=1)
        total = probs_raw.sum()
        p = probs_raw / total
        p = p / p.sum()  # re-normalize away roundoff so rng.choice never complains
        k = int(rng.choice(s.dim, p=p))
        vectors.append(evecs[:, k].copy())
        outcomes.append(float(evals[k]))
        if right_link is not None:
            L = ITensor((right_link,), rot[k, :] / np.sqrt(probs_raw[k]))
    new_cps = product_state(sites, vectors)
    return new_cps, outcomes


def imaginary_time_evolve(psi, H, dbeta_half, nsteps, cutoff, maxdim, niter=30):
    """Evolve psi by e^{-dbeta_half*H}, via nsteps of imaginary-time TDVP
    (dt=-1j*(dbeta_half/nsteps) passed to tdvp_step's own -i*dt real-time
    convention -- see this module's docstring for why that gives genuine
    decay rather than rotation), renormalizing after every step since
    imaginary-time evolution isn't norm-preserving (unlike the real-time
    case tdvp_step was originally written for)."""
    dt = -1j * (dbeta_half / nsteps)
    for _ in range(nsteps):
        psi = tdvp_step(psi, H, dt, cutoff=cutoff, maxdim=maxdim, niter=niter)
        psi.normalize()
    return psi


def _metts_single_chain(H, sites, ops, beta, nsamples, nwarmup,
                         dbeta_half_step, cutoff, maxdim, basis_ops,
                         initial_cps, seed, niter):
    """One full METTS Markov chain (nwarmup discarded + nsamples retained
    iterations), sequential -- the body metts_thermal_average() always ran
    before njobs= was added. Returns (means, stderrs, nsamples) -- the
    trailing nsamples is only needed by _pool_chain_results() below, to
    weight/combine several independent chains' statistics correctly."""
    rng = np.random.default_rng(seed)
    beta_half = beta / 2.0
    nsteps = max(1, int(np.ceil(beta_half / dbeta_half_step)))
    eigcache = _build_eigcache(sites, basis_ops)

    cps = initial_cps if initial_cps is not None else random_cps(sites, basis_ops[0], rng, eigcache)

    nbasis = len(basis_ops)
    samples = [[] for _ in ops]
    total_iters = nwarmup + nsamples
    for it in range(total_iters):
        phi = imaginary_time_evolve(cps, H, beta_half, nsteps, cutoff, maxdim, niter=niter)
        if it >= nwarmup:
            for j, A in enumerate(ops):
                samples[j].append(inner(phi, A, phi))
        basis = basis_ops[it % nbasis]
        cps, _ = collapse_to_cps(phi, sites, basis, rng, eigcache)

    means, stderrs = [], []
    for vals in samples:
        arr = np.array(vals)
        means.append(arr.mean())
        stderrs.append(arr.std(ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0.0)
    return means, stderrs, nsamples


def _max_index_id(sites, H, ops):
    """Highest Index.id reachable from sites/H/ops -- see index.py's
    reseed_id_counter_past() docstring for why _run_chain_worker() below
    needs this before it lets this (freshly spawned) process mint any
    Index of its own."""
    max_id = -1
    for i in range(1, sites.length() + 1):
        max_id = max(max_id, sites.si(i).id)
    for mpo in [H] + list(ops):
        for i in range(1, mpo.length() + 1):
            for ind in mpo.A(i).inds:
                max_id = max(max_id, ind.id)
    return max_id


def _run_chain_worker(args):
    """Top-level (picklable) target for multiprocessing.Pool -- a bound
    method or closure wouldn't survive pickling under the 'spawn' start
    method, so this just unpacks a plain tuple and calls
    _metts_single_chain() directly. Reseeds this worker's own Index id
    counter first, past every id already present in the H/sites/ops it
    just received (see index.py's reseed_id_counter_past() for why a
    freshly spawned worker's ids otherwise restart at 0 and can silently
    collide with those -- confirmed directly, this used to make njobs>1
    raise "shape-mismatch for sum" deep inside a worker's own tensor
    contraction, or worse, silently mismatch same-dimension legs)."""
    H, sites, ops = args[0], args[1], args[2]
    reseed_id_counter_past(_max_index_id(sites, H, ops))
    return _metts_single_chain(*args)


def _pool_chain_results(results, nops):
    """Combine several independent chains' (means, stderrs, nsamples)
    triples into the single (means, stderrs) pair metts_thermal_average()
    would have returned had all chains' retained samples been pooled into
    one flat array -- without ever needing the raw per-sample values
    themselves (each chain only reports its own mean/stderr/count).

    Exact, not approximate: writing v_k = stderr_k^2 * n_k for chain k's
    own (ddof=1) sample variance and M for the pooled mean, the identity
    sum_i |x_i - M|^2 = sum_k [ (n_k-1)*v_k + n_k*|m_k - M|^2 ]
    holds for any partition of a dataset into groups (the cross term
    sum_{i in k} (x_i - m_k) vanishes since m_k is group k's own mean) --
    this holds equally for complex-valued samples (METTS <Op> can have a
    small numerical-noise imaginary part) since |.|^2 is still the
    relevant squared norm and the same orthogonality argument applies."""
    means_out, stderrs_out = [], []
    ns = np.array([n for _, _, n in results], dtype=float)
    N = ns.sum()
    for j in range(nops):
        chain_means = np.array([means[j] for means, _, _ in results], dtype=complex)
        chain_stderr = np.array([stderrs[j] for _, stderrs, _ in results], dtype=float)
        chain_vars = chain_stderr ** 2 * ns
        M = np.sum(ns * chain_means) / N
        SS = np.sum((ns - 1) * chain_vars + ns * np.abs(chain_means - M) ** 2)
        var_pooled = SS / (N - 1) if N > 1 else 0.0
        means_out.append(M)
        stderrs_out.append(np.sqrt(var_pooled / N) if N > 0 else 0.0)
    return means_out, stderrs_out


def metts_thermal_average(H, sites, ops, beta, nsamples, nwarmup=20,
                           dbeta_half_step=0.05, cutoff=1e-10, maxdim=200,
                           basis_ops=("Sz",), initial_cps=None, seed=None,
                           niter=30, njobs=1):
    """Estimate Tr[A e^{-beta H}]/Tr[e^{-beta H}] for each MPO A in `ops`
    via METTS sampling (arXiv:1002.1305).

    beta: inverse temperature (1/T).
    nsamples: number of retained METTS samples, after nwarmup discarded
        equilibration steps of the Markov chain.
    dbeta_half_step: TDVP imaginary-time step size for the e^{-beta H/2}
        evolution, split into ceil((beta/2)/dbeta_half_step) equal steps.
    basis_ops: sequence of single-site operator names to cycle the
        collapse basis over, one per sample (e.g. ("Sz","Sx") for
        spin-1/2 -- see module docstring on why alternation matters).
    initial_cps: an MPS product state seeding the chain; a random CPS in
        basis_ops[0] if not given. Not supported together with njobs>1
        (each parallel chain seeds its own independent random CPS -- a
        single shared starting point wouldn't make sense split across
        several otherwise-independent Markov chains).
    njobs: run njobs independent METTS Markov chains in separate worker
        processes (multiprocessing.Pool) instead of one -- nsamples is
        split as evenly as possible across them, each with its own
        nwarmup equilibration (a fresh chain has no warmup to share) and
        an independently-seeded RNG derived from `seed` via
        numpy.random.SeedSequence.spawn(), then their statistics are
        combined exactly via _pool_chain_results() (no raw samples need
        crossing the process boundary). This trades more total warmup
        work (njobs*nwarmup steps instead of one chain's nwarmup) for
        wall-clock time on a multi-core machine -- worth it whenever
        nsamples/njobs is comfortably bigger than nwarmup, which is
        exactly the profiled bottleneck: METTS's per-sample cost here is
        dominated by this pure-Python tensor engine's own Python-level
        bookkeeping overhead (see tdvp.py/tensor.py), not shared BLAS
        work that would otherwise make multiprocessing redundant. njobs=1
        (the default) is unchanged from before this parameter existed --
        same single-process sequential chain, same seed reproducibility.

    Returns (means, stderrs): means[j] is the sample mean of ops[j],
    stderrs[j] its naive iid standard error over the nsamples retained
    samples. METTS samples are Markov-correlated, so this underestimates
    the true error unless dbeta_half_step/nwarmup are generous enough
    that consecutive samples actually decorrelate (paper, Sec. 3.3).
    """
    if beta <= 0:
        raise ValueError("metts_thermal_average: beta (1/T) must be > 0, got %r" % (beta,))
    if not basis_ops:
        raise ValueError("metts_thermal_average: basis_ops must be non-empty")
    if nsamples < 1:
        raise ValueError("metts_thermal_average: nsamples must be >= 1, got %r" % (nsamples,))
    if njobs < 1:
        raise ValueError("metts_thermal_average: njobs must be >= 1, got %r" % (njobs,))

    if njobs == 1:
        means, stderrs, _ = _metts_single_chain(
            H, sites, ops, beta, nsamples, nwarmup, dbeta_half_step,
            cutoff, maxdim, basis_ops, initial_cps, seed, niter)
        return means, stderrs

    if initial_cps is not None:
        raise ValueError("metts_thermal_average: initial_cps is not supported together with njobs>1")

    base, extra = divmod(nsamples, njobs)
    chunk_sizes = [base + (1 if k < extra else 0) for k in range(njobs)]
    chunk_sizes = [n for n in chunk_sizes if n > 0]
    if len(chunk_sizes) == 1:
        means, stderrs, _ = _metts_single_chain(
            H, sites, ops, beta, chunk_sizes[0], nwarmup, dbeta_half_step,
            cutoff, maxdim, basis_ops, None, seed, niter)
        return means, stderrs

    seeds = np.random.SeedSequence(seed).spawn(len(chunk_sizes))
    tasks = [
        (H, sites, ops, beta, n_local, nwarmup, dbeta_half_step, cutoff,
         maxdim, basis_ops, None, seeds[k], niter)
        for k, n_local in enumerate(chunk_sizes)
    ]
    # 'spawn', not the Linux default 'fork': confirmed directly this
    # matters, not just a style preference. fork()ing a Pool worker from a
    # parent process that has already run any numpy/scipy linear algebra
    # (which every caller here has -- at minimum this module's own eigh()
    # calls) duplicates whatever internal BLAS/OpenMP thread-pool state
    # that parent has already spun up; the classic fork-after-threading
    # hazard applies even though this process itself never touched the
    # threading module. Measured directly on a 6-site chain, nsamples=200:
    # the default fork-based Pool got monotonically *slower* as njobs grew
    # (34s/62s/73s/100s for njobs=1/2/4/8 -- each extra worker added
    # contention, not throughput), while 'spawn' (a fresh interpreter per
    # worker, no inherited thread-pool state) gave the expected
    # monotonic speedup on the same machine/workload (33s/21s/17s/14s).
    with multiprocessing.get_context("spawn").Pool(len(chunk_sizes)) as pool:
        results = pool.map(_run_chain_worker, tasks)
    return _pool_chain_results(results, len(ops))
