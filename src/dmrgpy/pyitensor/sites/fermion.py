"""Fermionic site types, transcribed operator-by-operator from
mpscpp3/ITensor/itensor/mps/sites/fermion.h (spinless) and electron.h
(spinful/Hubbard). See sites/base.py for the matrix-entry axis convention.

Both define "F" themselves (the actual fermionic sign/parity operator),
unlike every non-fermionic site type in this package -- that's what makes
SiteType.is_fermionic()'s name-based ("starts with 'C'") classification and
these two sites' operator names line up correctly for Jordan-Wigner string
threading in a later phase: is_fermionic() picks out exactly the operators
whose site type below actually carries a nontrivial F.
"""

import numpy as np

from .base import SiteType, build_matrix


class FermionSite(SiteType):
    """Spinless fermion, dim=2 (Emp=1, Occ=2). mpscpp3's SpinlessSite alias."""
    dim = 2
    _states = {"Emp": 1, "0": 1, "Occ": 2, "1": 2}
    _OPS = {
        "N": build_matrix(2, [(2, 2, 1)]),
        "n": build_matrix(2, [(2, 2, 1)]),
        "C": build_matrix(2, [(2, 1, 1)]),
        "Cdag": build_matrix(2, [(1, 2, 1)]),
        "A": build_matrix(2, [(2, 1, 1)]),
        "Adag": build_matrix(2, [(1, 2, 1)]),
        "F": build_matrix(2, [(1, 1, 1), (2, 2, -1)]),
        "FermiPhase": build_matrix(2, [(1, 1, 1), (2, 2, -1)]),
        "projEmp": build_matrix(2, [(1, 1, 1)]),
        "projOcc": build_matrix(2, [(2, 2, 1)]),
    }


SpinlessSite = FermionSite


class ElectronSite(SiteType):
    """Spinful fermion / Hubbard site, dim=4 (Emp=1, Up=2, Dn=3, UpDn=4).
    mpscpp3's HubbardSite alias.

    Cdn/Cdagdn (and Adagdn but not Adn -- transcribed exactly as electron.h
    has it) carry an extra intrinsic minus sign on the doubly-occupied
    matrix element relative to Cup/Cdagup: this is the same-site up/down
    fermion-ordering convention baked into ITensor's own ElectronSite, not
    a Jordan-Wigner string (those are threaded across *different* sites in
    a later phase) -- preserved here exactly as electron.h defines it.
    """
    dim = 4
    _states = {"0": 1, "Emp": 1, "+": 2, "Up": 2, "-": 3, "Dn": 3, "S": 4, "UpDn": 4}
    _OPS = {
        "Nup": build_matrix(4, [(2, 2, 1), (4, 4, 1)]),
        "Ndn": build_matrix(4, [(3, 3, 1), (4, 4, 1)]),
        "Nupdn": build_matrix(4, [(4, 4, 1)]),
        "Ntot": build_matrix(4, [(2, 2, 1), (3, 3, 1), (4, 4, 2)]),
        "Cup": build_matrix(4, [(2, 1, 1), (4, 3, 1)]),
        "Cdagup": build_matrix(4, [(1, 2, 1), (3, 4, 1)]),
        "Cdn": build_matrix(4, [(3, 1, 1), (4, 2, -1)]),
        "Cdagdn": build_matrix(4, [(1, 3, 1), (2, 4, -1)]),
        "Aup": build_matrix(4, [(2, 1, 1), (4, 3, 1)]),
        "Adagup": build_matrix(4, [(1, 2, 1), (3, 4, 1)]),
        "Adn": build_matrix(4, [(3, 1, 1), (4, 2, 1)]),
        "Adagdn": build_matrix(4, [(1, 3, 1), (2, 4, 1)]),
        "F": build_matrix(4, [(1, 1, 1), (2, 2, -1), (3, 3, -1), (4, 4, 1)]),
        "FermiPhase": build_matrix(4, [(1, 1, 1), (2, 2, -1), (3, 3, -1), (4, 4, 1)]),
        "Fup": build_matrix(4, [(1, 1, 1), (2, 2, -1), (3, 3, 1), (4, 4, -1)]),
        "Fdn": build_matrix(4, [(1, 1, 1), (2, 2, 1), (3, 3, -1), (4, 4, -1)]),
        "Sz": build_matrix(4, [(2, 2, 0.5), (3, 3, -0.5)]),
        "Sp": build_matrix(4, [(3, 2, 1)]),
        "S+": build_matrix(4, [(3, 2, 1)]),
        "Sm": build_matrix(4, [(2, 3, 1)]),
        "S-": build_matrix(4, [(2, 3, 1)]),
        "S2": build_matrix(4, [(2, 2, 0.75), (3, 3, 0.75)]),
    }


HubbardSite = ElectronSite
