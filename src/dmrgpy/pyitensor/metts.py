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

import numpy as np

from .index import Index
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


def random_cps(sites, opname, rng):
    """A random classical product state: at each site, uniformly pick one
    eigenstate of the named single-site operator (e.g. "Sz"). Valid seed
    for the METTS Markov chain regardless of choice -- the chain's
    stationary distribution doesn't depend on the starting CPS (paper,
    Sec. 2), only how many warmup steps are needed to reach it."""
    n = sites.length()
    vectors = []
    for i in range(1, n + 1):
        _, evecs = _eigenbasis(sites, opname, i)
        k = rng.integers(0, sites.dim(i))
        vectors.append(evecs[:, k].copy())
    return product_state(sites, vectors)


def collapse_to_cps(psi, sites, opname, rng):
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
        evals, evecs = _eigenbasis(sites, opname, i)
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


def metts_thermal_average(H, sites, ops, beta, nsamples, nwarmup=20,
                           dbeta_half_step=0.05, cutoff=1e-10, maxdim=200,
                           basis_ops=("Sz",), initial_cps=None, seed=None,
                           niter=30):
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
        basis_ops[0] if not given.

    Returns (means, stderrs): means[j] is the sample mean of ops[j],
    stderrs[j] its naive iid standard error over the nsamples retained
    samples. METTS samples are Markov-correlated, so this underestimates
    the true error unless dbeta_half_step/nwarmup are generous enough
    that consecutive samples actually decorrelate (paper, Sec. 3.3).
    """
    rng = np.random.default_rng(seed)
    beta_half = beta / 2.0
    nsteps = max(1, int(np.ceil(beta_half / dbeta_half_step)))

    cps = initial_cps if initial_cps is not None else random_cps(sites, basis_ops[0], rng)

    nbasis = len(basis_ops)
    samples = [[] for _ in ops]
    total_iters = nwarmup + nsamples
    for it in range(total_iters):
        phi = imaginary_time_evolve(cps, H, beta_half, nsteps, cutoff, maxdim, niter=niter)
        if it >= nwarmup:
            for j, A in enumerate(ops):
                samples[j].append(inner(phi, A, phi))
        basis = basis_ops[it % nbasis]
        cps, _ = collapse_to_cps(phi, sites, basis, rng)

    means, stderrs = [], []
    for vals in samples:
        arr = np.array(vals)
        means.append(arr.mean())
        stderrs.append(arr.std(ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0.0)
    return means, stderrs
