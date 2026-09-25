"""Canonical form for MultiOperator terms.

A MultiOperator term is a coefficient times an ordered product of named
single-site operators, [c,[name,site],[name,site],...]. Two terms that
name the same product in a different factor order are the same operator,
up to a sign when the reordering exchanges two fermionic factors, but
nothing in multioperator.py knows that: terms are compared, collected and
cancelled by their literal spelling. That is what this module fixes --
it rewrites every term into one canonical spelling (factors sorted by
site, fermionic sign accounted for, identities dropped) so that equal
operators get equal signatures and can be collected in a dict.

This replaces the sympy round trip that MultiOperator.simplify() used to
go through (multioperatortk/sympymultioperator.py, removed with it, and
with sympy itself as a dependency). That one mapped each factor to a
non-commutative sympy Symbol, which carries no algebra beyond "these do
not commute", so it collected terms that were already spelled
identically and nothing else: measured on an
n=20 spin-1/2 chain, MultiOperator.is_hermitian() took 0.834 s to report
False for an ordinary Sx[i]Sx[j]+Sy[i]Sy[j]+Sz[i]Sz[j] Heisenberg
Hamiltonian, the false negative that mpsalgebra.exponential and
infinitechain.set_hamiltonian both document having hit.

What the canonical form proves, and what it does not
----------------------------------------------------
The rewrite is exact, so an empty canonical form is a proof: if the
canonical form of H-H^dagger has no surviving term, H is Hermitian, full
stop, and no numerical witness is needed. The converse does not hold. Two
spellings can be the same operator for reasons this module deliberately
does not model, so a nonempty canonical form means "not proven", not
"not Hermitian":

- same-site identities, Sx Sx = 1/4 on a spin-1/2 site, N N = N, and
  every other product rule that depends on the local Hilbert space the
  chain was built with, which a MultiOperator does not know;
- aliases between names, Sp = Sx + i Sy, or Sz = (Nup-Ndn)/2 on a
  spinful fermionic site, so a single factor can already be a sum of
  other single factors;
- the order of two factors sharing a site, which is only canonicalized
  when both are diagonal (see _DIAGONAL), since operators on one site do
  not commute in general and get_dagger() reverses their order. So
  Nup[i]*Ndn[i], the Hubbard U term, is proven, while a term resting on
  two non-diagonal same-site factors commuting is not;
- the order of a pre-Jordan-Wigner C-type factor and a bare local
  A-type one (see _LOCAL_LADDER), whose exchange is neither a
  commutation nor an anticommutation but depends on which of the two
  sits on the lower site, so a term mixing the two is not reordered at
  all: C0*A1-A1*C0 is exactly zero and is not proven zero.

Modelling the first three needs the site type, which lives on the chain
and not on the operator, and modelling the fourth needs the string
convention of the Jordan-Wigner transform, which the infinite-chain path
does not use. So every prover here is one-sided: True means proven,
False means not proven. mpsalgebra.is_hermitian consumes exactly that,
trusting a True and falling back to its random-witness probe otherwise.

Which names it understands
--------------------------
Only names in _PARITY below, which are the ones dmrgpy itself builds.
Each entry says whether that operator is fermionic (odd, picking up a
sign when it is moved past another odd factor on a different site) or
bosonic (even, commuting freely across sites), and every one of them
also has a known adjoint (either self-adjoint, or paired in
multioperator._dagger_name). A term naming anything else -- a
parafermionic Sig/Tau, which reorders with a Z_n phase rather than a
sign, or a name a caller invented -- is left exactly as it was written
and only ever collects with a term spelled the same way, which is what
the sympy version did for everything. The same goes for a term that
names both a C-type fermion and a bare A-type ladder operator, since
the sign of exchanging those two is not a grading at all. Nothing is
silently reordered on a grading this module cannot check.
"""

from .. import multioperator


# Fermionic names: odd under exchange with another odd factor on a
# different site. These are the pre-Jordan-Wigner names, the ones the ED
# backends materialize with the statistics built into the occupation
# basis (pyfermion/mbfermion.py's get_c/get_cd) and the ones
# multioperator.jordan_wigner() dresses with a string.
_ODD = ("C", "Cdag", "Cup", "Cdagup", "Cdn", "Cdagdn")

# The bare post-transform names, plain local matrices with no string:
# the fermionic a_i that jordan_wigner() writes C_i in terms of, and on
# Bosonic_Chain/SpinBoson_Chain the boson ladder operators. They commute
# with each other, with F, N and the spin names on another site, so they
# are listed as even below (giving them no parity instead would lose the
# proof of every boson hopping and of every Jordan-Wigner-transformed
# operator). What they are NOT is even against a C-type name: with
# C_j = F_0...F_{j-1} A_j, A_i anticommutes with C_j when j>i and
# commutes with it when j<i, a relation that depends on which of the two
# sits on the lower site and is not a parity at all. A term naming both
# kinds is therefore left as written, see _mixes_representations.
_LOCAL_LADDER = frozenset(("A", "Adag", "Aup", "Adagup", "Adn", "Adagdn"))

_EVEN = ("Id", "X", "Y", "Z", "Sx", "Sy", "Sz", "Sp", "Sm", "S+", "S-",
         "N", "density", "Nup", "Ndn", "Ntot",
         "A", "Adag", "Aup", "Adagup", "Adn", "Adagdn",
         "F", "Fup", "Fdn")

_PARITY = dict([(n, 1) for n in _ODD] + [(n, 0) for n in _EVEN])


def _mixes_representations(factors, ps):
    """True if a term names both a pre-Jordan-Wigner C-type fermion and
    a bare A-type ladder operator, whose exchange sign this module does
    not know (see _LOCAL_LADDER). ps are the factors' parities, and the
    odd ones are exactly the C-type names."""
    return 1 in ps and any(name in _LOCAL_LADDER for (name, i) in factors)


# Names that are diagonal in their own site's basis, and therefore
# commute with each other when they share a site. Sorting by site alone
# leaves a term like Nup[i]*Ndn[i] spelled differently from its own
# dagger, which reverses the factor order -- that is the Hubbard U term
# on a native spinful site, so without this the whole Hubbard
# Hamiltonian is not proven Hermitian even though its hopping part is.
# Only these names may be reordered within a site: Sx[i] and Sz[i] are
# both even and do NOT commute, so a general same-site sort would be
# wrong.
_DIAGONAL = frozenset(("Id", "Z", "Sz", "N", "density", "Nup", "Ndn",
                       "Ntot", "F", "Fup", "Fdn"))


def is_diagonal(name):
    """True if this operator is diagonal in its site's own basis, so
    that it commutes with every other diagonal operator on that site."""
    if name in _DIAGONAL: return True
    # boson occupation projectors, diagonal by construction
    return len(name) > 1 and name[0] == "N" and name[1:].isdigit()


def parity(name):
    """Return 1 for a fermionic name, 0 for a bosonic one, and None when
    the name is not one this module is willing to reorder."""
    p = _PARITY.get(name, None)
    if p is not None: return p
    # boson occupation projectors N0,N1,...,N{maxnb-1} (bosonchain.py's
    # self.D), a family rather than a fixed list, all of them diagonal
    if len(name) > 1 and name[0] == "N" and name[1:].isdigit(): return 0
    return None


def _canonical_signature(term):
    """Rewrite one term [c,[name,site],...] into (signature,coefficient),
    the signature being a tuple of (name,site) pairs in canonical order.

    Factors are sorted by site index with a stable sort, so factors
    sharing a site keep the order they were written in (their product
    does not commute and this module does not multiply them out), and the
    coefficient picks up a minus sign for every pair of fermionic factors
    the sort exchanges. Identity factors are dropped, unless the whole
    term is identities and there would be nothing left.

    A term naming an operator of unknown grading, or mixing a C-type
    fermion with a bare A-type ladder operator, is returned spelled
    exactly as it came in, so it still collects with an identical term
    and is never reordered.
    """
    c = term[0]
    factors = [(o[0], o[1]) for o in term[1:]]
    ps = [parity(name) for (name, i) in factors]
    if any(p is None for p in ps): # not ours to reorder
        return tuple(factors), c
    if _mixes_representations(factors, ps): # sign depends on site order
        return tuple(factors), c
    keep = [(f, p) for (f, p) in zip(factors, ps) if f[0] != "Id"]
    if len(keep) == 0: # nothing but identities, leave the term alone
        return tuple(factors), c
    # sign of the permutation, restricted to exchanges of two fermions
    n = len(keep)
    for a in range(n):
        for b in range(a+1, n):
            if keep[a][0][1] > keep[b][0][1] and keep[a][1] and keep[b][1]:
                c = -c
    keep.sort(key=lambda fp: fp[0][1]) # stable, so same-site order stands
    return tuple(_sort_diagonal_runs([f for (f, p) in keep])), c


def _sort_diagonal_runs(factors):
    """Given factors already sorted by site, put the ones sharing a site
    into a canonical order too, but only where every factor on that site
    is diagonal and they therefore commute. A site carrying anything
    else keeps the order it was written in."""
    out = []
    i = 0
    n = len(factors)
    while i < n:
        j = i
        while j < n and factors[j][1] == factors[i][1]: j += 1
        run = factors[i:j]
        if len(run) > 1 and all(is_diagonal(name) for (name, site) in run):
            run = sorted(run)
        out.extend(run)
        i = j
    return out


def canonical_dict(MO):
    """Return {signature: coefficient} for a MultiOperator, with terms
    of equal signature summed and near-zero coefficients dropped.

    "Near-zero" is relative to what went INTO each sum, not to what came
    out of it: a summed coefficient is dropped when it is at or below
    multioperator.clean_threshold (1e-12) times the larger of
      (a) the largest |coefficient| among MO's own terms, and
      (b) the sum of the |coefficients| collected into that signature.
    Taking the scale from the sums instead would keep everything in the
    one case this exists for, H - H^dagger of a Hermitian H, whose sums
    are exact zeros or rounding dust with nothing larger left to compare
    against. (b) is what that dust scales with: 2000 copies of one term
    minus their adjoints leave 2.1e-12 of the largest coefficient, above
    (a) alone, but 1.1e-15 of their own sum. (a) is the floor
    _filter_small applies at every consumption point (write(),
    to_terms(), MO2matrix), so a term every backend drops cannot keep an
    operator from being proven zero. Both are scale-free, so the answer
    does not depend on the units the operator is written in; the absolute
    1e-8 this replaced proved 1e-9*(Sz0+1j*Sx0) Hermitian and 1e-9*Sz0
    zero (2026-09-25 audit, clean-threshold)."""
    cmax = max([abs(t[0]) for t in MO.op]+[0.])
    out = dict()
    weight = dict()
    for term in MO.op:
        sig, c = _canonical_signature(term)
        out[sig] = out.get(sig, 0.0) + c
        weight[sig] = weight.get(sig, 0.0) + abs(c)
    tol = multioperator.clean_threshold
    return dict([(sig, c) for (sig, c) in out.items()
                 if abs(c) > tol*max(cmax, weight[sig])])


def canonicalize(MO):
    """Return a MultiOperator in canonical form: the same operator, with
    every term rewritten into its canonical spelling and terms that are
    equal collected into one."""
    d = canonical_dict(MO)
    out = multioperator.MultiOperator(term=False)
    out.name = MO.name
    out.op = [[complex(c)] + [list(f) for f in sig] for (sig, c) in d.items()]
    out.i = len(out.op)-1
    if len(out.op) == 0: # everything cancelled
        out.op = [[0.0]]
        out.i = 0
    return out


def is_zero(MO):
    """True if this operator is provably zero. False means not proven,
    see this module's docstring."""
    return len(canonical_dict(MO)) == 0


def all_names_known(MO):
    """True if every factor of every term names an operator listed in
    _PARITY, so that both its grading and its adjoint are known here."""
    for term in MO.op:
        for o in term[1:]:
            if parity(o[0]) is None: return False
    return True


def is_hermitian(MO):
    """True if this operator is provably Hermitian. False means not
    proven, see this module's docstring."""
    # The adjoint has to be known before the cancellation means anything:
    # multioperator.get_dagger() leaves a name it does not recognize
    # untouched, so an operator built out of such a name would cancel
    # against its own "dagger" and be proven Hermitian when it is not.
    if not all_names_known(MO): return False
    return is_zero(MO - MO.get_dagger())


def is_antihermitian(MO):
    """True if this operator is provably anti-Hermitian. False means not
    proven, see this module's docstring."""
    if not all_names_known(MO): return False
    return is_zero(MO + MO.get_dagger())


def is_dagger_pair(A, B):
    """True if A is provably B^dagger. False means not proven.

    This is is_zero(A-B.get_dagger()) with the same guard is_hermitian()
    carries, and for the same reason: get_dagger() leaves a name it does
    not recognize untouched, so an operator built out of such a name
    would cancel against its own "dagger" and be declared the adjoint of
    itself. A caller asking this question is asking about adjoints, so
    it has to refuse when the adjoint is not known, which plain
    is_zero() (a statement about the operator it is handed, not about
    anybody's dagger) has no reason to do."""
    if not (all_names_known(A) and all_names_known(B)): return False
    return is_zero(A - B.get_dagger())
