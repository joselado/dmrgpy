from .. import multioperator

# Jordan-Wigner transform for Hamiltonians built directly out of
# flavor-resolved spinful-fermion operators (Cup/Cdn/Cdagup/Cdagdn) living
# on a *native* spinful site: one physical MPS/ED site per orbital, with a
# local Hilbert space of dimension 4 (Empty, Up, Down, UpDn -- ITensor v3's
# own "Electron"/"Hubbard" site, mpscpp3/get_sites.h's site-type code 1),
# instead of two separate spinless-fermion sites per orbital (site-type
# code 0, see fermionchain.Spinful_Fermionic_Chain). See
# fermionchain.Spinful_Fermionic_Chain_Native.
#
# Bare (string-free) on-site operator names used here (Aup/Adagup/Adn/
# Adagdn/Fup/Fdn/F) are exactly the ones ITensor's ElectronSite already
# implements (mpscpp3/ITensor/itensor/mps/sites/electron.h) -- this module
# only ever *emits* those names, it never talks to ITensor directly, so the
# same term list is equally valid input for the ED backend once
# pyfermion.mbfermion.MBFermion.get_operator understands them (it doesn't
# need Aup/Adagup/etc: ED operators are built straight from the un-dressed
# Cup/Cdn/Cdagup/Cdagdn names, since ED's occupation-number basis already
# encodes fermion statistics exactly, with no Jordan-Wigner string needed
# at all -- see MO2matrix()/get_ED_obj(), which never call jordan_wigner()).
#
# Mode ordering: each physical site i contributes two fermionic modes,
# "up" and "down", threaded in that fixed order -- equivalent to the flat
# mode numbering 2*i (up), 2*i+1 (down) used by the interleaved
# Spinful_Fermionic_Chain (fermionchain.py), so this reproduces bit-for-bit
# the same physical convention/sign, just packaged into half as many
# tensor-network sites. This is also exactly ITensor AutoMPO's own
# internal fermionic rewrite table for "Cup"/"Cdn"/"Cdagup"/"Cdagdn"
# (autompo.cc's fermionicTerm()): Cdagup->Adagup, Cup->Aup,
# Cdagdn->Adagdn*Fup, Cdn->Adn*Fup -- confirming this is the standard,
# already-battle-tested sign convention, not something invented here.

def obj2MO(a): return multioperator.obj2MO([a])

_BARE = {"Cup":"Aup", "Cdagup":"Adagup", "Cdn":"Adn", "Cdagdn":"Adagdn"}
_SPIN = {"Cup":0, "Cdagup":0, "Cdn":1, "Cdagdn":1} # 0 = up, 1 = down


def one_fermion(name,site):
    """Jordan-Wigner-dressed single spinful-fermion operator: thread a
    full on-site parity "F" through every site strictly before `site`,
    plus (for the "down" flavor only) a partial "Fup" for `site` itself,
    since the up mode of that same site precedes the down mode in the
    fixed per-site ordering used throughout this module."""
    spin = _SPIN[name]
    m = 1.0
    for k in range(site):
        m = obj2MO(["F",k])*m
    if spin==1: m = obj2MO(["Fup",site])*m
    return m*obj2MO([_BARE[name],site])


def product2jw(names,sites):
    """Jordan-Wigner-dress a product of operators given in a fixed order
    (matching the desired physical operator ordering). Any factor that is
    a spinful fermionic operator (Cup/Cdn/Cdagup/Cdagdn) is dressed
    independently with its own preceding string -- the standard
    order/multiplicity-safe recipe a_i = (prod_{k<i} F_k) c_i, applied
    factor by factor, which is valid regardless of what other operators
    appear in the same product. Any other operator name (Nup, Ndn, Ntot,
    Sz, Id, ...) passes through untouched: those are already diagonal and
    commute freely with Jordan-Wigner strings, needing no dressing at all.
    """
    out = 1
    for (name,site) in zip(names,sites):
        if name in _BARE: out = out*one_fermion(name,site)
        else: out = out*obj2MO([name,site])
    return out
