"""SiteType: a local Hilbert space plus its table of named single-site
operators, as a plain (dim,dim) complex matrix per operator name.

Every matrix here uses the same axis convention as the ITensor v3 site
headers it was transcribed from (mpscpp3/ITensor/itensor/mps/sites/*.h and
mpscpp3/extra/*.h): each of those files builds its operators as
`Op = ITensor(dag(s), prime(s))` and then calls `Op.set(unprimed_state,
primed_state, value)` per nonzero entry -- i.e. the *unprimed* leg is the
"input"/ket index and the *primed* leg is the "output" index, matching
mpscpp3/chain_session.h's own apply_mpo() convention of always noPrime()-ing
the physical index after an MPO application to restore it to a plain ket.
So here, `matrix(name)[a, b]` is the value ITensor's `Op.set(state(a+1),
prime(state(b+1)), value)` would have set: contracting axis 0 (unprimed)
against a ket's physical index and reading off axis 1 (primed) as the new
one reproduces the same operator.
"""

import numpy as np


def build_matrix(dim, entries):
    """entries: iterable of (in_state, out_state, value), 1-based (matching
    the site headers' own IndexVal numbering)."""
    m = np.zeros((dim, dim), dtype=complex)
    for i, j, v in entries:
        m[i - 1, j - 1] = v
    return m


class SiteType:
    dim = None
    _states = {}
    _OPS = None  # set by each subclass to a dict name -> (dim,dim) ndarray

    @classmethod
    def state(cls, name):
        """1-based index of a named basis state -- mirrors each site
        header's own state(string) method."""
        return cls._states[name]

    @classmethod
    def matrix(cls, opname):
        if opname == "Id":
            return np.eye(cls.dim, dtype=complex)
        ops = cls._OPS
        if opname in ops:
            return ops[opname]
        if opname in ("F", "FermiPhase"):
            # SiteSet::op()'s own generic fallback (siteset.h): a site type
            # that never defines "F" itself (every non-fermionic site here)
            # gets a trivial identity F, so Jordan-Wigner string threading
            # through it in a later phase is a no-op, as it should be for a
            # non-fermionic local Hilbert space.
            return np.eye(cls.dim, dtype=complex)
        raise ValueError("Operator '{}' not recognized for site type {}".format(
            opname, cls.__name__))

    @staticmethod
    def is_fermionic(opname):
        """Name-based classification, mirrors ITensor's own isFermionic()
        (mpscpp3/ITensor/itensor/mps/autompo.cc): any operator name
        starting with 'C' is fermionic for Jordan-Wigner purposes, on every
        site type, independent of what that site's own operator table
        actually contains."""
        return bool(opname) and opname[0] == "C"


def is_fermionic(opname):
    """Free-function form of SiteType.is_fermionic -- this classification
    never actually depends on which site it's evaluated at (autompo.cc's
    isFermionic(SiteTerm) only ever looks at the operator name), so a
    Phase-3 AutoMPO port can call this directly instead of routing through
    a particular site instance."""
    return SiteType.is_fermionic(opname)
