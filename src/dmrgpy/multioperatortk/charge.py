"""What quantum-number charge a MultiOperator carries, computed in Python
and independently of which backend the chain runs on.

The C++ backend already knows how to do this -- `Chain::op_charge` in
mpscpp3/chain_session.h infers a single-site operator's charge from its
own dense matrix elements -- but it knows it *inside* a live session, one
operator name at a time, and only for `itensor_version=3`. The
sector-resolved dynamical correlator (sectordc.py) has to decide, before
any session is touched and for `itensor_version="python"` too, which
sector `B|gs>` lands in. That decision is pure site-type bookkeeping: the
charge tables of the local Hilbert spaces plus the operator's own matrix
elements, neither of which needs an MPS.

pyitensor's site machinery has exactly those tables and is pure Python, so
it is used here as a *library* rather than as a backend -- an ungraded
`SiteX` built from the chain's site-type codes serves any chain, whichever
backend actually runs the DMRG. `pyitensor/sector.py`'s own `op_charge`,
`expand_xy_terms` and `combine_terms` are reused verbatim, so a term list
is judged by the same rules that decide whether a sector will accept it.

The Sx/Sy expansion matters more than it looks: the textbook Heisenberg
Hamiltonian Sx.Sx+Sy.Sy+Sz.Sz has terms whose individual charge is not
even defined (`Sx` both raises and lowers), and only the *sum* conserves
Sz. Running expand+combine first is what lets `conserved_qn_names()`
answer "yes, this Hamiltonian conserves Sz" for it instead of "I cannot
tell".
"""

import numpy as np


def chain_sites(chain):
    """An ungraded pyitensor `SiteX` for this chain's site-type codes.

    Nothing here is graded by a quantum number: the charge *tables* are a
    property of the site types (`state_charge`), not of the grading, so
    this works for every backend and costs nothing to build. Cached on the
    chain, keyed by the site codes, since a correlator sweep asks for it
    once per operator.
    """
    key = tuple(chain.sites)
    cached = getattr(chain, "_charge_sites_cache", None)
    if cached is not None and cached[0] == key:
        return cached[1]
    from ..pyitensor.sites import SiteX
    sites = SiteX(list(chain.sites))
    chain._charge_sites_cache = (key, sites)
    return sites


def offered_qn_names(chain):
    """Every conserved quantity some site of this chain can carry, in a
    stable order. Empty for a chain with no sector support at all (the
    parafermion sites)."""
    from ..pyitensor.sites import site_qn_names
    out = []
    for code in chain.sites:
        for name in site_qn_names(code):
            if name not in out:
                out.append(name)
    return tuple(out)


def terms_charge(sites, terms, names):
    """The charge a term list carries, as one integer per name in `names`,
    or None if it carries no definite charge.

    `terms` is the `(coefficient, [(opname, 1-based site), ...])` form
    `MultiOperator.to_terms()` produces -- i.e. already Jordan-Wigner
    transformed, which is deliberate: the JW strings are exactly what the
    backend builds, and the string operator `F` is diagonal, so it adds
    nothing to the charge and the answer is the same either way.

    Terms are expanded (Sx/Sy -> S+/S-) and combined before being judged,
    so an operator that conserves a quantity only as a sum of terms is
    recognized as conserving it.
    """
    from ..pyitensor.sector import combine_terms, expand_xy_terms, op_charge
    names = tuple(names)
    zero = tuple(0 for _ in names)
    found = None
    for _, factors in combine_terms(expand_xy_terms(sites, terms)):
        total = zero
        for opname, site in factors:
            q = op_charge(sites, opname, site, names)
            if q is None:  # e.g. a surviving Sx: raises and lowers at once
                return None
            total = tuple(a + b for a, b in zip(total, q))
        if found is None:
            found = total
        elif total != found:
            return None  # different terms move different amounts of charge
    return zero if found is None else found


def multioperator_charge(chain, MO, names):
    """The charge of a MultiOperator, one integer per name in `names`, or
    None if it has no definite charge (`Sx`, `C+Cdag`, ...)."""
    from ..multioperator import obj2MO
    MO = obj2MO(MO)
    return terms_charge(chain_sites(chain), MO.to_terms(), names)


def conserved_qn_names(chain, operator, candidates=None):
    """Which of `candidates` (default: everything the sites offer) the
    given operator -- in practice the Hamiltonian -- actually conserves.

    This is the right default for a sector-resolved calculation: taking
    every quantity the site types offer instead is correct but fails late,
    deep inside the Hamiltonian send, on a model that breaks one of them
    (a Hubbard chain with spin-orbit coupling conserves Nf but not Sz).
    """
    if candidates is None:
        candidates = offered_qn_names(chain)
    sites = chain_sites(chain)
    from ..multioperator import obj2MO
    terms = obj2MO(operator).to_terms()
    out = []
    for name in candidates:
        q = terms_charge(sites, terms, (name,))
        if q is not None and q[0] == 0:
            out.append(name)
    return tuple(out)


def sector_dimension(chain, sector):
    """How many product basis states carry exactly the charges in
    `sector` (a dict or a list of (name,target) pairs) -- i.e. the
    dimension of that sector's Hilbert space.

    Counted by the same dynamic program over reachable partial charges
    that `pyitensor/sector.py::sector_state_plan` uses to *find* one such
    state, carrying a count per charge tuple instead of a representative.
    The number of distinct reachable charge tuples is small (it grows
    polynomially with the chain length, not exponentially), so this is
    cheap even for long chains.

    A dynamical correlator needs it to cap `nex`: asking for more states
    than the target sector holds sends the overlap-penalty search hunting
    for states that do not exist, which is slow, noisy and produces
    nothing.
    """
    items = sorted(sector.items()) if isinstance(sector, dict) else sorted(sector)
    names = tuple(k for k, _ in items)
    target = tuple(int(v) for _, v in items)
    if not names:
        return None  # no sector: the full Hilbert space, not counted here
    sites = chain_sites(chain)
    counts = {tuple(0 for _ in names): 1}
    for i in range(1, sites.length() + 1):
        new = {}
        for partial, c in counts.items():
            for st in range(1, sites.dim(i) + 1):
                q = sites.state_charge(i, st, names)
                key = tuple(a + b for a, b in zip(partial, q))
                new[key] = new.get(key, 0) + c
        counts = new
    return counts.get(target, 0)


def charge_operator(chain, name):
    """The total-charge operator of a conserved quantity as a
    MultiOperator, or None when this chain does not expose one.

    Only used to *infer* the reference sector when the caller did not say
    which one they meant (sectordc.py measures it on an unconstrained
    ground state and rounds). Note the units: `Sz` is in ITensor's integer
    2*Sz convention throughout the sector API, so the operator returned
    here is 2*sum(Sz), not sum(Sz).
    """
    if name in ("Nf", "Nb"):
        for attr in ("Ntot", "N"):
            ops = getattr(chain, attr, None)
            if ops is not None and len(ops) > 0:
                return sum(ops)
        return None
    if name == "Sz":
        ops = getattr(chain, "Sz", None)
        if ops is not None and len(ops) > 0:
            return 2 * sum(ops)
        return None
    return None


# -- decomposition into charge-homogeneous parts -------------------------
#
# The charge components of the two operators whose charge is not definite
# in the first place. `Sx` is the standing case and the one that matters:
# it is how essentially every spin correlator in dmrgpy's own examples is
# written, and on an Sz-graded chain it both raises and lowers, so a
# sector-resolved correlator cannot use it directly. But
#
#     Sx = (S+ + S-)/2,    S+ = Sx + i Sy,    S- = Sx - i Sy
#
# so Sx splits exactly into two pieces of definite charge, each of which
# is still expressible in the operator names the backends actually have
# (dmrgpy has no native S+/S-: operatornames.name2MO builds them from
# Sx/Sy too). Splitting both operators and pairing up the components whose
# charges cancel turns one impossible correlator into a sum of two
# perfectly ordinary ones, over the Sz+2 and Sz-2 sectors.

_XY_PIECES = {
    # name -> {charge sign: [(coefficient, native name), ...]}, from
    # Sx = (S+ + S-)/2 and Sy = -i(S+ - S-)/2.
    "Sx": {+1: [(0.5, "Sx"), (0.5j, "Sy")],
           -1: [(0.5, "Sx"), (-0.5j, "Sy")]},
    "Sy": {+1: [(-0.5j, "Sx"), (0.5, "Sy")],
           -1: [(0.5j, "Sx"), (0.5, "Sy")]},
}


def charge_components(chain, MO, names):
    """Split a MultiOperator into parts of definite charge:
    `{charge_tuple: MultiOperator}`, summing exactly back to the original.

    Returns None when some factor cannot be split this way -- an operator
    with no definite charge on a site that has no S+/S- to expand it into,
    which no caller can do anything with either.

    Note the pieces are labelled during construction rather than measured
    afterwards, and each is a *sum* of native-operator terms (the +1 piece
    of Sx is (Sx + i Sy)/2, two terms) whose individual factors are still
    charge-indefinite. That is fine and is exactly how sector mode's own
    Hamiltonian check works: `terms_charge` expands and combines before
    judging, so it sees S+/2 and answers +2, not "cannot tell".
    """
    from ..multioperator import obj2MO
    from ..pyitensor.sector import _site_has_ladder_ops, op_charge
    MO = obj2MO(MO)
    sites = chain_sites(chain)
    names = tuple(names)
    zero = tuple(0 for _ in names)
    groups = {}
    for term in MO.op:
        # (coefficient, factors so far, charge so far); one entry per way
        # of choosing a charge piece for each factor seen so far
        partial = [(complex(term[0]), [], zero)]
        for (opname, i) in term[1:]:
            site = i + 1
            if opname in _XY_PIECES and _site_has_ladder_ops(sites, site):
                ladder = {+1: "S+", -1: "S-"}
                options = []
                for sign, pieces in _XY_PIECES[opname].items():
                    q = op_charge(sites, ladder[sign], site, names)
                    if q is None:
                        return None
                    options.append((q, pieces))
            else:
                q = op_charge(sites, opname, site, names)
                if q is None:
                    return None  # not an Sx/Sy we know how to split
                options = [(q, [(1.0, opname)])]
            grown = []
            for coef, factors, charge in partial:
                for q, pieces in options:
                    for c, native in pieces:
                        grown.append((coef * c, factors + [(native, i)],
                                      tuple(a + b for a, b in zip(charge, q))))
            partial = grown
        for coef, factors, charge in partial:
            if coef == 0:
                continue
            groups.setdefault(charge, []).append([coef] + factors)
    out = {}
    for charge, terms in groups.items():
        piece = MO.copy()
        piece.op = terms
        out[charge] = piece
    return out
