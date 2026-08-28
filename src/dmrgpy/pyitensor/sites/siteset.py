"""SiteX: mirrors mpscpp3/get_sites.h's SpinX(std::vector<int> site_types)
constructor -- a SiteSet built directly from in-memory per-site type codes,
with no file I/O, using the exact same type-code convention as
manybodychain.py's callers (see get_sites.h's comment and CLAUDE.md):
2=spin-1/2, 0=spinless fermion, 1=spinful fermion (Hubbard), 3=spin-1,
4=spin-3/2, 5=spin-2, 6=spin-5/2, -2=Z3, -3=Z4. Boson is a *range* of
codes, 102..199 (code = 100+dim, dim being the requested local Hilbert
space dimension -- see bosonchain.Bosonic_Chain's "maxnb"), routed to
boson.get_boson_site(dim) below instead of a single dict entry -- the
same generalization mpscpp3/get_sites.h applies to BosonFourSite.
"""

from ..index import Index
from ..tensor import ITensor
from .boson import BosonFourSite, get_boson_site
from .fermion import ElectronSite, FermionSite
from .parafermion import Z3Site, Z4Site
from .spin import SpinFiveHalfSite, SpinHalfSite, SpinOneSite, SpinThreeHalfSite, SpinTwoSite

TYPE_CODE_TO_SITE = {
    2: SpinHalfSite,
    0: FermionSite,
    1: ElectronSite,
    3: SpinOneSite,
    4: SpinThreeHalfSite,
    5: SpinTwoSite,
    6: SpinFiveHalfSite,
    104: BosonFourSite,
    -2: Z3Site,
    -3: Z4Site,
}


def _site_class(type_code):
    if 100 < type_code < 200:
        return get_boson_site(type_code - 100)
    try:
        return TYPE_CODE_TO_SITE[type_code]
    except KeyError as e:
        raise ValueError("SiteX cannot build a site of type code {}".format(e.args[0]))


def site_qn_names(type_code):
    """The conserved quantities a site-type code can carry -- the
    pyitensor counterpart of mpscpp3/get_sites.h's site_qn_names(), reading
    the per-site-type charge tables instead of duplicating them. An empty
    tuple means no conserved-sector support for that site type (the Z3/Z4
    parafermion sites: they do carry a Z_n charge, but nothing in dmrgpy's
    Python layer can name a target for it)."""
    return _site_class(type_code).qn_names()


class SiteX:
    """The site set. `conserved`, when given, is the set of quantum-number
    names the caller wants conserved (conserved-sector mode, see
    sector.py): every site index is then built carrying a charge grading
    for exactly those of its own quantities that appear in it, and the
    whole site set gets fresh Index identities -- so a wavefunction built
    in a sector cannot be silently contracted against an operator built
    outside one. Mirrors mpscpp3/get_sites.h's SpinX(site_types,conserved),
    including its rule that a site type contributing none of the requested
    quantities is an error rather than a silently ungraded site.
    """

    def __init__(self, site_types, conserved=None):
        self._types = [_site_class(t) for t in site_types]
        self._codes = tuple(site_types)
        self._conserved = None if not conserved else tuple(sorted(conserved))
        self._indices = [Index(t.dim, tags="Site,n={}".format(i + 1),
                                qnnames=self._names_at(i + 1),
                                qn=self._charges_at(i + 1))
                          for i, t in enumerate(self._types)]

    def _names_at(self, i):
        """Which of site i's own quantum numbers this site set grades it
        by: the intersection of the request with what the site type
        offers, or None outside sector mode."""
        if self._conserved is None:
            return None
        own = self._types[i - 1].qn_names()
        names = tuple(n for n in self._conserved if n in own)
        if not names:
            raise ValueError(
                "SiteX: site %d (type code %s) conserves none of %s -- it "
                "offers %s. A site set mixing graded and ungraded sites is "
                "not something conserved-sector mode supports."
                % (i, self._codes[i - 1], ", ".join(self._conserved),
                   ", ".join(own) if own else "no quantum numbers at all"))
        return names

    def _charges_at(self, i):
        names = self._names_at(i)
        if names is None:
            return None
        t = self._types[i - 1]
        return tuple(tuple(t.charges(n)[st] for n in names)
                      for st in range(t.dim))

    def conserved(self):
        """The quantum numbers this site set is graded by (a tuple), or
        None if it is an ordinary dense site set."""
        return self._conserved

    def site_codes(self):
        return self._codes

    def qn_names(self, i):
        """The quantum numbers site i is graded by, or None."""
        return self._indices[i - 1].qnnames

    def state_charge(self, i, st, names):
        """The charge of (1-based) basis state `st` of site i, as one
        integer per name in `names` -- 0 for a quantity this site type
        doesn't carry, exactly as adding a spin site's Sz to a spinless
        fermion site's (absent) Sz would. Mirrors chain_session.h's
        state_charge()."""
        t = self._types[i - 1]
        return tuple(int(t.charges(n)[st - 1]) if n in t._QN else 0
                      for n in names)

    def length(self):
        return len(self._types)

    N = length

    def si(self, i):
        """1-based physical Index at site i."""
        return self._indices[i - 1]

    def site_type(self, i):
        return self._types[i - 1]

    def dim(self, i):
        return self._types[i - 1].dim

    def state_index(self, i, name):
        """1-based basis-state index of a named state at site i."""
        return self._types[i - 1].state(name)

    def op(self, opname, i, args=None):
        """The (unprimed=in, primed=out) ITensor for `opname` at (1-based)
        site i -- mirrors SiteSet::op()/siteset.h."""
        s = self.si(i)
        mat = self._types[i - 1].matrix(opname)
        return ITensor((s, s.prime(1)), mat)
