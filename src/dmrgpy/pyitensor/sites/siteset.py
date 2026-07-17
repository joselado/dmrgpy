"""SiteX: mirrors mpscpp3/get_sites.h's SpinX(std::vector<int> site_types)
constructor -- a SiteSet built directly from in-memory per-site type codes,
with no file I/O, using the exact same type-code convention as
manybodychain.py's callers (see get_sites.h's comment and CLAUDE.md):
2=spin-1/2, 0=spinless fermion, 1=spinful fermion (Hubbard), 3=spin-1,
4=spin-3/2, 5=spin-2, 6=spin-5/2, 104=Boson (4 levels), -2=Z3, -3=Z4.
"""

from ..index import Index
from ..tensor import ITensor
from .boson import BosonFourSite
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


class SiteX:
    def __init__(self, site_types):
        try:
            self._types = [TYPE_CODE_TO_SITE[t] for t in site_types]
        except KeyError as e:
            raise ValueError("SiteX cannot build a site of type code {}".format(e.args[0]))
        self._indices = [Index(t.dim, tags="Site,n={}".format(i + 1))
                          for i, t in enumerate(self._types)]

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
