"""Conserved-sector (quantum-number) mode: confining every calculation on
a chain to one total-charge sector, and leaving it again afterwards.

This is the pyitensor counterpart of mpscpp3/chain_session.h's
`set_conserved_sector`/`promote_to_dense` (and of mo_terms.h's
`expand_xy_terms`/`combine_terms`), reached from Python by exactly the
same `Many_Body_Chain.set_conserved_sector(Nf=6)` / `(Sz=0)` call. What a
sector means to a caller is identical on both backends -- "the ground
state at exactly six particles", with every operator built while it is set
checked for conserving the requested quantities, and `promote_to_dense()`
to leave the sector while keeping the state computed inside it.

How it is *implemented* differs, and the difference is worth stating
plainly rather than discovering later.

ITensor v3 confines the calculation structurally: a QN-carrying Index
sorts its basis states into per-charge blocks, every tensor built on it is
block-sparse, and an amplitude outside the sector has nowhere to be
stored. That is also where v3's speedup comes from (measured 2.0x at
n=40, 4.2x at n=60, see get_sites.h) -- and where its sharp edges come
from, since building a flux-violating operator aborts the process.

pyitensor keeps storage dense (see index.py: only site indices carry a
grading, Link indices don't, and nothing here is block-sparse), so there
is no speedup, and confinement has to be arranged rather than inherited.
Two things arrange it, and they cover different failure modes:

  * The Hamiltonian is normalized and charge-checked before it is built
    (`sector_terms`), exactly as v3 does -- so the operator actually being
    diagonalized commutes with every conserved charge, and DMRG started in
    a sector has no matrix element out of it. This is what makes the
    answer the sector's answer.

  * A charge penalty lambda*sum_k (Q_k - q_k)^2 is added to the operator
    the *variational* solves minimize (`charge_penalty_mpo`, used by
    gs_energy/excited_states only). This is the part v3 doesn't need. The
    reason it is needed here: contraction preserves exact zeros, but a
    dense LAPACK SVD does not -- Householder bidiagonalization mixes rows
    across charge blocks, so every truncation injects ~1e-16 of amplitude
    into neighboring sectors. Harmless when the target sector happens to
    contain the global ground state (a Heisenberg chain at Sz=0), and not
    harmless at all when it doesn't (a fermion chain away from the filling
    that minimizes the energy, the case this feature exists for): a
    variational sweep amplifies that dust geometrically toward the lower
    sector. The penalty is exactly zero on the target sector -- it changes
    no reported number, and `Chain.gs_energy` reports <H> under the
    unpenalized H regardless -- while turning any leaked amplitude into an
    energy cost that self-damps. It also makes the solve well-posed from
    any start state, in-sector or not.

Nothing outside the variational solves needs the penalty: time evolution,
KPM and correlators are linear maps of an in-sector state under
charge-conserving operators, with no minimization to amplify dust.
"""

import numpy as np

from .index import Index
from .mpscontainer import MPO, MPS
from .tensor import ITensor

# Relative tolerance for "this coefficient cancelled exactly", mirroring
# mo_terms.h's combine_terms(): an exact cancellation leaves rounding dust,
# a genuinely small term does not.
_COMBINE_TOL = 1e-12
_ELEMENT_TOL = 1e-12  # "this matrix element is nonzero", as in op_charge()
MAX_XY_FACTORS = 4  # 2^4 = 16 expanded terms worst case, as in mo_terms.h


# -- term-list normalization (mo_terms.h's expand_xy_terms/combine_terms) --

def _site_has_ladder_ops(sites, site):
    try:
        sites.site_type(site).matrix("S+")
        sites.site_type(site).matrix("S-")
    except Exception:
        return False
    return True


def _expansion_of(name):
    """Sx = (S+ + S-)/2, Sy = (S+ - S-)/(2i) = -i/2 S+ + i/2 S-."""
    if name == "Sx":
        return ((0.5 + 0j, "S+"), (0.5 + 0j, "S-"))
    return ((-0.5j, "S+"), (0.5j, "S-"))


def expand_xy_terms(sites, terms):
    """Expand every Sx/Sy factor into its two S+/S- pieces, term by term.

    Exact for any spin S, and a pure representation change. Sector mode
    needs it because the textbook Heisenberg Hamiltonian Sx.Sx+Sy.Sy+Sz.Sz
    conserves Sz only as a *sum*: term by term it does not, and the
    offending S+S+/S-S- strings only disappear once the Sx and Sy
    expansions are added together (see combine_terms). Factors are
    replaced in place, so the surviving factor order -- which matters for
    non-commuting same-site products -- is untouched; terms with no Sx/Sy
    factor, and factors on sites that define no S+/S- at all (the Z3/Z4
    clock operators, where "Sy" is not a spin ladder operator in the first
    place), pass through verbatim.
    """
    out = []
    for coef, factors in terms:
        factors = list(factors)
        xy_pos = [k for k, (name, site) in enumerate(factors)
                  if name in ("Sx", "Sy") and _site_has_ladder_ops(sites, site)]
        if not xy_pos:
            out.append((complex(coef), list(factors)))
            continue
        if len(xy_pos) > MAX_XY_FACTORS:
            raise ValueError(
                "conserved-sector mode: a term carries %d Sx/Sy factors, more "
                "than the %d this expansion supports -- write it with S+/S- "
                "instead" % (len(xy_pos), MAX_XY_FACTORS))
        for mask in range(1 << len(xy_pos)):
            c = complex(coef)
            new = list(factors)
            for b, k in enumerate(xy_pos):
                pc, pname = _expansion_of(factors[k][0])[(mask >> b) & 1]
                c *= pc
                new[k] = (pname, factors[k][1])
            if c == 0:
                continue
            out.append((c, new))
    return out


def combine_terms(terms):
    """Sum the coefficients of identical operator strings and drop the ones
    that cancel, keeping first-appearance order otherwise. Two strings
    count as identical only if their (name,site) factor *sequences* match
    exactly: no reordering is attempted, since that would mean tracking
    fermionic anticommutation signs, and every term list dmrgpy's
    MultiOperator produces already comes out in ascending site order (so
    the Sx.Sx/Sy.Sy pair this exists to cancel does line up). Mirrors
    mo_terms.h's combine_terms()."""
    maxabs = max([abs(complex(c)) for c, _ in terms], default=0.0)
    tol = _COMBINE_TOL * maxabs
    seen, out = {}, []
    for coef, factors in terms:
        key = tuple((name, site) for name, site in factors)
        if key in seen:
            i = seen[key]
            out[i] = (out[i][0] + complex(coef), out[i][1])
        else:
            seen[key] = len(out)
            out.append((complex(coef), list(factors)))
    return [(c, f) for c, f in out if abs(c) > tol]


# -- charges -------------------------------------------------------------

def op_charge(sites, name, site, names):
    """The charge a named single-site operator carries, as one integer per
    quantity in `names`, or None if it has no definite charge at all.

    Inferred from the operator's own matrix elements plus the charge
    labels the site set gives each local basis state: element (a,b) maps
    basis state a to basis state b (see sites/base.py's axis convention),
    so it shifts the charge by q(b)-q(a), and the operator has a definite
    charge only if every nonzero element shifts it the same way. `Sx` on
    Sz-graded spin sites is the standing counterexample -- its nonzero
    elements raise *and* lower. An identically-zero operator carries no
    charge at all, exactly as in chain_session.h's op_charge().
    """
    mat = sites.site_type(site).matrix(name)
    found = None
    d = sites.dim(site)
    for a in range(1, d + 1):
        for b in range(1, d + 1):
            if abs(mat[a - 1, b - 1]) < _ELEMENT_TOL:
                continue
            delta = tuple(x - y for x, y in zip(sites.state_charge(site, b, names),
                                                sites.state_charge(site, a, names)))
            if found is None:
                found = delta
            elif delta != found:
                return None
    return found if found is not None else tuple(0 for _ in names)


def _sector_string(sector):
    return ", ".join("%s=%d" % (k, v) for k, v in sector) or "(no sector)"


def _term_string(factors):
    return "*".join("%s(%d)" % (name, site) for name, site in factors) \
        or "(empty term)"


def _charge_string(names, charge):
    parts = ["%s=%d" % (n, c) for n, c in zip(names, charge) if c != 0]
    return ", ".join(parts) or "(no charge)"


def sector_terms(sites, sector, terms):
    """Normalize a term list for sector mode and reject what still doesn't
    conserve the sector, naming it.

    Normalization is expand_xy_terms + combine_terms (see those): it is
    what lets the textbook Sx.Sx+Sy.Sy+Sz.Sz Heisenberg Hamiltonian
    through, since the strings that break Sz individually cancel exactly
    once expanded. Whatever survives that and still carries a net charge is
    an operator the caller cannot have while a sector is set -- a bare `C`
    on an Nf-conserving chain, a transverse field on an Sz-conserving one
    -- and raises ValueError rather than quietly returning an answer from
    the wrong Hilbert space.

    `sector` is a sorted list of (name, target) pairs.
    """
    names = tuple(k for k, _ in sector)
    out = combine_terms(expand_xy_terms(sites, terms))
    for coef, factors in out:
        total = tuple(0 for _ in names)
        for name, site in factors:
            q = op_charge(sites, name, site, names)
            if q is None:
                raise ValueError(
                    "conserved-sector mode (%s): operator \"%s\" at site %d "
                    "carries no definite charge, so no operator containing it "
                    "can be built while a sector is set. Write the model with "
                    "S+/S- (or Cdag/C) instead, or clear the sector."
                    % (_sector_string(sector), name, site))
            total = tuple(a + b for a, b in zip(total, q))
        if any(total):
            raise ValueError(
                "conserved-sector mode (%s): the term %s changes the conserved "
                "quantum numbers by %s -- every operator built on a chain with "
                "a sector set must conserve them. (Terms are named after the "
                "exact Sx/Sy -> S+/S- expansion and the cancellation of strings "
                "that only conserve the sector as a sum, so this may not be "
                "spelled the way it was written.)"
                % (_sector_string(sector), _term_string(factors),
                   _charge_string(names, total)))
    return out


# -- in-sector starting states -------------------------------------------

def sector_state_plan(sites, sector):
    """One per-site basis-state assignment (1-based) whose charges add up
    to exactly the requested sector.

    Both the reachability check set_conserved_sector() runs -- so an
    impossible target errors at request time rather than mid-sweep -- and
    the arrangement the starting state is built from. A small dynamic
    program over reachable partial charges rather than a greedy fill,
    because with more than one conserved quantity (a Hubbard chain at
    fixed Nf *and* Sz) greedy choices dead-end; each layer keeps one
    representative per distinct partial charge, so its size is bounded by
    the number of reachable charge tuples, not by the number of
    assignments. Transcribed from chain_session.h's sector_state_plan().
    """
    names = tuple(k for k, _ in sector)
    target = tuple(v for _, v in sector)
    n, nq = sites.length(), len(names)

    # per-suffix reachable charge range, to prune partial sums that can no
    # longer close the gap
    smin = [[0] * nq for _ in range(n + 2)]
    smax = [[0] * nq for _ in range(n + 2)]
    for i in range(n, 0, -1):
        charges = [sites.state_charge(i, st, names)
                   for st in range(1, sites.dim(i) + 1)]
        for k in range(nq):
            col = [c[k] for c in charges]
            smin[i][k] = min(col) + smin[i + 1][k]
            smax[i][k] = max(col) + smax[i + 1][k]

    layers = [[(tuple([0] * nq), -1, 0)]]  # (partial charge, prev index, state)
    for i in range(1, n + 1):
        seen, layer = set(), []
        for pi, (charge, _, _) in enumerate(layers[i - 1]):
            for st in range(1, sites.dim(i) + 1):
                c = sites.state_charge(i, st, names)
                partial = tuple(a + b for a, b in zip(charge, c))
                if any(partial[k] + smin[i + 1][k] > target[k]
                        or partial[k] + smax[i + 1][k] < target[k]
                        for k in range(nq)):
                    continue
                if partial in seen:
                    continue  # one representative is enough
                seen.add(partial)
                layer.append((partial, pi, st))
        if not layer:
            raise ValueError(
                "conserved sector: no state of this chain has %s -- the "
                "requested sector is empty" % _sector_string(sector))
        layers.append(layer)

    plan = [1] * n
    at = 0
    for i in range(n, 0, -1):
        _, prev, st = layers[i][at]
        plan[i - 1] = st
        at = prev
    return plan


def sector_state_arrangements(sites, sector, k, seed=0):
    """`k` in-sector product states: the same multiset of local basis
    states as sector_state_plan()'s, shuffled among the sites that can host
    them. Permuting within one site-type code keeps every state legal on
    the site it lands on, and any permutation is charge-neutral by
    construction, so all of these live in the requested sector."""
    base = sector_state_plan(sites, sector)
    codes = sites.site_codes()
    positions = {}
    for i, code in enumerate(codes):
        positions.setdefault(code, []).append(i)
    rng = np.random.RandomState(1234 + int(seed))
    out = [list(base)]  # the deterministic plan is always one of them
    for _ in range(1, k):
        plan = list(base)
        for pos in positions.values():
            vals = [base[i] for i in pos]
            rng.shuffle(vals)
            for j, i in enumerate(pos):
                plan[i] = vals[j]
        out.append(plan)
    return out


def sector_product_mps(sites, states):
    """A bond-dimension-1 product MPS from a per-site basis-state
    assignment (1-based states). Boundary sites carry no trivial edge link,
    matching this package's own MPS convention (mpscontainer.py)."""
    n = sites.length()
    links = [Index(1, tags="Link,l={}".format(i + 1)) for i in range(n - 1)]
    tensors = []
    for i in range(1, n + 1):
        s = sites.si(i)
        inds = ([links[i - 2]] if i > 1 else []) + [s] \
            + ([links[i - 1]] if i < n else [])
        arr = np.zeros(tuple(ind.dim for ind in inds), dtype=complex)
        idx = tuple(states[i - 1] - 1 if ind.hastags("Site") else 0
                     for ind in inds)
        arr[idx] = 1.0
        tensors.append(ITensor(tuple(inds), arr))
    psi = MPS(tensors)
    psi.center = 1
    return psi


def sector_mps(sites, sector, m, seed=0, cutoff=0.0):
    """The starting state for a sector-mode solve: a normalized random
    superposition of up to `m` in-sector product states.

    A *single* product state would be a poor start for the same reason
    chain_session.h's default_mps() doesn't use one -- this package's DMRG
    has no noise term at all (see dmrg.py), so a start that is an exact
    eigenstate of the Hamiltonian has nothing to move it. Summing many of
    them keeps every amplitude inside the sector (a sum of in-sector states
    is in-sector) while giving the sweep a state with real bond dimension
    to work with, mirroring chain_session.h's sector_mps().

    The coefficients are random rather than all equal, which is not
    cosmetic: an *equal-weight* symmetric sum is itself an exact
    eigenstate of a symmetric Hamiltonian often enough to matter, and a
    noise-free DMRG started on one never leaves it. Measured directly --
    a 2-site hopping chain at Nf=1 has exactly two in-sector product
    states, and (|10>+|01>)/sqrt(2) is the +1 eigenvector, so the
    equal-weight start returned the *highest* state of the sector instead
    of the lowest. Random coefficients break that symmetry while staying
    exactly in the sector, since every term of the sum is.
    """
    from .mpsalgebra import sum_many
    m = max(1, int(m))
    # Draw from *more* arrangements than the bond dimension asked for and
    # let the truncating sum mix them, rather than summing exactly m of
    # them. Summing exactly m leaves a bond basis whose vectors are single
    # basis paths -- a sparse, highly structured state; oversampling and
    # truncating down to m gives a generic one of the same bond dimension,
    # which both converges faster and leaks less charge on the way. Measured
    # on a 40-site Heisenberg chain at Sz=0 (maxdim 100, the ramp starting
    # the state at 10), oversampled vs not: <(Q-q)^2> after one sweep
    # 3.5e-6 vs 5.1e-5, and the two-sweep energy -17.5182 vs -17.4115
    # against a converged -17.5414733. No measurable time cost -- the extra
    # product states are bond dimension 1 and the sum truncates immediately.
    plans = sector_state_arrangements(sites, sector, 4 * m + 8, seed=seed)
    # distinct plans only: summing the same product state twice adds bond
    # dimension carrying no new direction at all
    unique, seen = [], set()
    for plan in plans:
        key = tuple(plan)
        if key not in seen:
            seen.add(key)
            unique.append(plan)
    rng = np.random.RandomState(9876 + int(seed))
    states = []
    for plan in unique:
        psi = sector_product_mps(sites, plan)
        states.append(psi * complex(rng.randn(), rng.randn()))
    psi = states[0] if len(states) == 1 else sum_many(states, cutoff=cutoff,
                                                      maxdim=m)
    psi.position(psi.length())
    psi.normalize()
    return psi


# -- the charge penalty --------------------------------------------------

def charge_penalty_mpo(sites, sector, lam):
    """lambda * sum_k (Q_k - q_k)^2 as an MPO, where Q_k is the total charge
    operator of conserved quantity k and q_k its requested target.

    Exactly zero on the target sector, so it changes nothing about the
    answer there; positive everywhere else, so it costs any amplitude that
    leaks out of the sector (see this module's docstring for why dense
    storage needs that and ITensor v3 doesn't). Built directly as one
    bond-dimension-3 MPO per conserved quantity -- the standard
    (sum_i B_i)^2 automaton, with the constant folded into site 1's B --
    rather than by handing the O(N^2) two-site product terms to
    mpobuilder.py, which would build and compress N^2 one-term MPOs for an
    operator whose exact form is known in closed form.
    """
    from .mpsalgebra import sum_many
    names = tuple(k for k, _ in sector)
    parts = []
    for k, (name, target) in enumerate(sector):
        parts.append(_single_charge_penalty_mpo(sites, names, k, target, lam))
    if len(parts) == 1:
        return parts[0]
    return sum_many(parts, cutoff=0.0)


def _single_charge_penalty_mpo(sites, names, k, target, lam):
    n = sites.length()
    # B_i: the diagonal charge operator of quantity k at site i, with the
    # target subtracted at site 1 so that sum_i B_i = Q_k - q_k exactly.
    def B(i):
        d = sites.dim(i)
        diag = np.array([sites.state_charge(i, st, names)[k]
                          for st in range(1, d + 1)], dtype=complex)
        if i == 1:
            diag = diag - target
        return np.diag(diag)

    links = [Index(3, tags="Link,l={}".format(i + 1)) for i in range(n - 1)]
    tensors = []
    for i in range(1, n + 1):
        s = sites.si(i)
        d = s.dim
        eye = np.eye(d, dtype=complex)
        b = B(i)
        # W = [[Id, B, B^2], [0, Id, 2B], [0, 0, Id]] -- a path from the
        # left state 0 to the right state 2 places either one B^2 or two
        # Bs, which is exactly (sum_i B_i)^2 summed over all placements.
        W = np.zeros((3, d, d, 3), dtype=complex)
        W[0, :, :, 0] = eye
        W[0, :, :, 1] = b
        W[0, :, :, 2] = lam * (b @ b)
        W[1, :, :, 1] = eye
        W[1, :, :, 2] = 2.0 * lam * b
        W[2, :, :, 2] = eye
        if i == 1:
            arr, inds = W[0], (s, s.prime(1), links[0])
        elif i == n:
            arr, inds = W[:, :, :, 2], (links[n - 2], s, s.prime(1))
        else:
            arr, inds = W, (links[i - 2], s, s.prime(1), links[i - 1])
        tensors.append(ITensor(inds, arr))
    mpo = MPO(tensors)
    mpo.center = 1
    return mpo


def terms_scale(terms):
    """A charge-penalty strength of the same order as the Hamiltonian's own
    energy scale: the sum of the magnitudes of its coefficients, which
    bounds its spectral range. The penalty has to outweigh that range for
    the target sector's ground state to be the *global* minimum of
    H + penalty (a wrong sector one unit away costs lambda); much larger
    than that only inflates the effective Hamiltonian's spectral range and
    slows Lanczos down."""
    return max(1.0, float(sum(abs(complex(c)) for c, _ in terms)))


# -- leaving the sector --------------------------------------------------

def promote_tensor(T, mapping):
    """One ITensor with its graded site indices replaced by their ungraded
    counterparts (`mapping`: index id -> replacement Index), prime levels
    preserved. The array is untouched: the two indices enumerate the same
    local basis states in the same order (both come from the same site
    type -- only the grading differs), so this is an exact relabeling."""
    inds = tuple(mapping[ind.id].setprime(ind.plev) if ind.id in mapping else ind
                  for ind in T.inds)
    return ITensor(inds, T.array)


def promote_chain(chain, sites_to):
    """An MPS/MPO rebased onto `sites_to`'s site indices -- the pyitensor
    counterpart of chain_session.h's to_dense_mps() (removeQNs +
    replaceSiteInds). Exact and free: with dense storage the two states
    have identical arrays, and only the index identities (and their charge
    labels) differ.

    Each tensor's own Site-tagged index is replaced by the target set's
    index for that position, rather than by looking the old one up in the
    site set it came from -- so this still works on a state handed back
    long after the sector that produced it was cleared, which is exactly
    what Many_Body_Chain.promote_mps() is for."""
    out = chain.copy()
    for i in range(1, out.length() + 1):
        T = out.A(i)
        ids = set(ind.id for ind in T.inds if ind.hastags("Site"))
        mapping = {j: sites_to.si(i) for j in ids}
        out.set_A(i, promote_tensor(T, mapping))
    return out


def chain_site_index(psi):
    """The Site-tagged index of an MPS/MPO's first tensor -- used to tell
    whether a wavefunction belongs to the site set a Chain currently
    has."""
    for ind in psi.A(1).inds:
        if ind.hastags("Site") and ind.plev == 0:
            return ind
    return None
