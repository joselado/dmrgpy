"""AutoMPO term compiler: HTerm/AutoMPO + Jordan-Wigner threading.

Transcribes the term-building half of
mpscpp3/ITensor/itensor/mps/autompo.cc -- isFermionic() (name-based: any
operator name starting with 'C' is fermionic), HTerm::add()'s
insertion-sort-with-sign (reorders a term's factors into site-ascending
order, picking up the correct fermion anticommutation sign along the way),
and rewriteFermionic()'s F-string insertion (mpscpp3/chain_session.h's own
comment block at the top of that file has more on why v3's inner()/innerC()
split and applyMPO()/nmultMPO() priming conventions look the way they do;
this module only concerns itself with the term-compilation side) -- plus
mo_terms.h's build_ampo(), which turns dmrgpy's own MultiOperator shape
(coef, [(opname,site), ...]) into these HTerms.

What this module deliberately does *not* port: ITensor's own automaton
MPO-compression algorithm (autompo.cc's actual SVD/state-merging MPO
builder). Reproducing that exactly is unnecessary here -- since dmrgpy
never needs bit-for-bit matching internal bond dimensions between backends
(only matching physical results, bounded by whatever MaxDim/Cutoff the
caller already passes), mpo.py builds one exact, trivial bond-dimension-1
MPO per HTerm and sum-compresses them together via ordinary MPO addition
(mps.py/mpo.py's sum()), which is simpler to get right and asymptotically
equivalent for the modest term counts dmrgpy's own Hamiltonians have.

Every single-site matrix handed around in this module is in *standard*
(physicist, output-row/input-column) convention -- i.e. the transpose of
the (in,out) storage convention documented in sites/base.py -- since that's
what ordinary matrix multiplication composes correctly (see resolve()).
"""

import numpy as np

from .sites.base import is_fermionic


class HTerm:
    """One term: a coefficient times a product of named single-site
    operators. `ops` is kept in site-ascending order, stable for factors at
    the same site, exactly mirroring HTerm::add()'s incremental insertion
    (build it by calling add() in the term's original, possibly unsorted,
    factor order -- the sign correction only comes out right that way)."""

    def __init__(self, coef=1.0):
        self.coef = complex(coef)
        self.ops = []  # [(opname, site), ...], 1-based sites

    def add(self, opname, site, x=1.0):
        idx = 0
        while idx < len(self.ops) and self.ops[idx][1] <= site:
            idx += 1
        # A fermionic operator inserted ahead of operators already fixed at
        # larger sites must anticommute past every fermionic one of them,
        # each swap contributing a sign -- equivalent to (and cheaper than)
        # actually performing the swaps one at a time.
        if idx < len(self.ops) and is_fermionic(opname):
            right_ops = self.ops[idx:]
            if sum(1 for name, _ in right_ops if is_fermionic(name)) % 2 == 1:
                self.coef *= -1
        self.coef *= x
        self.ops.insert(idx, (opname, site))
        return self

    def resolve(self, sites):
        """The full per-site list of standard-convention (dim_i,dim_i)
        matrices (one per site, 1..sites.length()) that reconstructs this
        term's action once Kronecker-multiplied (dense, for testing) or
        assembled into a bond-dimension-1 MPO (mpo.py, production path) --
        mirrors autompo.cc's rewriteFermionic(): sites this term doesn't
        touch get plain Id, except while an odd number of fermionic
        operators have been "seen" scanning left to right, when they get F
        instead, to thread the Jordan-Wigner string. A touched site
        composes its own factors (in their original relative order) and,
        exactly like rewriteFermionic, may pick up an *extra* trailing F of
        its own depending on the fermion parity carried in from its left
        versus its own -- see this module's docstring; trust this
        transcription over any "obviously it should be X" intuition, and
        check examples/pyitensor/selftest_autompo.py's cross-checks against
        dmrgpy's independent ED backend if this ever looks suspect.
        """
        by_site = {}
        for name, site in self.ops:
            by_site.setdefault(site, []).append(name)

        n = sites.length()
        out = [None] * n
        carry = False  # fermion parity carried in from sites < current
        for site in range(1, n + 1):
            site_type = sites.site_type(site)
            names = by_site.get(site)
            if names is None:
                out[site - 1] = site_type.matrix("F").T if carry else site_type.matrix("Id").T
                continue
            is_site_fermionic = sum(1 for nm in names if is_fermionic(nm)) % 2 == 1
            need_F = (carry != is_site_fermionic)
            mats = [site_type.matrix(nm).T for nm in names]
            if need_F:
                mats.append(site_type.matrix("F").T)
            combined = mats[0]
            for m in mats[1:]:
                combined = combined @ m
            out[site - 1] = combined
            carry = (carry != is_site_fermionic)
        return out


class AutoMPO:
    """Accumulates HTerms, mirroring the `ampo += coef,"Name",site,...`
    builder used throughout mpscpp3/chain_session.h. Python has no comma
    operator, so `add()` (or `+=` with an explicit tuple) takes the same
    factors as separate arguments instead."""

    def __init__(self, sites):
        self.sites = sites
        self.terms = []

    def add(self, coef, *factors):
        if len(factors) % 2 != 0:
            raise ValueError("AutoMPO.add: factors must be (name,site) pairs")
        term = HTerm(coef)
        it = iter(factors)
        for name, site in zip(it, it):
            term.add(name, site)
        self.terms.append(term)
        return self

    def __iadd__(self, value):
        self.add(*value)
        return self

    @classmethod
    def from_terms(cls, sites, terms):
        """terms: [(coef, [(opname,site), ...]), ...] -- exactly
        MultiOperator.to_terms()'s shape (and bindings.cc's
        terms_from_python's target type), matching mo_terms.h's
        build_ampo()."""
        ampo = cls(sites)
        for coef, factors in terms:
            term = HTerm(coef)
            for name, site in factors:
                term.add(name, site)
            ampo.terms.append(term)
        return ampo

    def dense_matrix(self):
        """The full dense matrix by Kronecker product over resolve()'d
        per-site factors, in ascending site order (site 1 leftmost/most
        significant -- the same convention as edtk/one2many.py's own
        np.kron loop, so this is directly comparable to dmrgpy's ED
        backend). Exponential in system size: exists purely as a
        small-system correctness reference for
        examples/pyitensor/selftest_autompo.py, not a production path --
        mpo.py's to_mpo() is what actually builds an MPO for DMRG/TDVP."""
        n = self.sites.length()
        total = 1
        for i in range(1, n + 1):
            total *= self.sites.dim(i)
        out = np.zeros((total, total), dtype=complex)
        for term in self.terms:
            mats = term.resolve(self.sites)
            full = mats[0]
            for m in mats[1:]:
                full = np.kron(full, m)
            out += term.coef * full
        return out
