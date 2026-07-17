"""Spin site types (S=1/2, 1, 3/2, 2, 5/2), transcribed operator-by-operator
from mpscpp3/ITensor/itensor/mps/sites/spinhalf.h, spinone.h, spintwo.h and
mpscpp3/extra/spinthreehalf.h, spinfivehalf.h. See sites/base.py for the
matrix-entry axis convention.
"""

import numpy as np

from .base import SiteType, build_matrix

ISQRT2 = 1.0 / np.sqrt(2.0)
SQRT2 = np.sqrt(2.0)


class SpinHalfSite(SiteType):
    dim = 2
    _states = {"Up": 1, "Dn": 2}
    _OPS = {
        "Sz": build_matrix(2, [(1, 1, 0.5), (2, 2, -0.5)]),
        "Sx": build_matrix(2, [(1, 2, 0.5), (2, 1, 0.5)]),
        "ISy": build_matrix(2, [(1, 2, -0.5), (2, 1, 0.5)]),
        "Sy": build_matrix(2, [(1, 2, 0.5j), (2, 1, -0.5j)]),
        "Sp": build_matrix(2, [(2, 1, 1)]),
        "S+": build_matrix(2, [(2, 1, 1)]),
        "Sm": build_matrix(2, [(1, 2, 1)]),
        "S-": build_matrix(2, [(1, 2, 1)]),
        "projUp": build_matrix(2, [(1, 1, 1)]),
        "projDn": build_matrix(2, [(2, 2, 1)]),
        "S2": build_matrix(2, [(1, 1, 0.75), (2, 2, 0.75)]),
    }


class SpinOneSite(SiteType):
    dim = 3
    _states = {"Up": 1, "+": 1, "Z0": 2, "0": 2, "Dn": 3, "-": 3}
    _OPS = {
        "Sz": build_matrix(3, [(1, 1, 1.0), (3, 3, -1.0)]),
        "Sx": build_matrix(3, [(1, 2, ISQRT2), (2, 1, ISQRT2), (2, 3, ISQRT2), (3, 2, ISQRT2)]),
        "ISy": build_matrix(3, [(1, 2, -ISQRT2), (2, 1, ISQRT2), (2, 3, -ISQRT2), (3, 2, ISQRT2)]),
        "Sy": build_matrix(3, [(1, 2, ISQRT2 * 1j), (2, 1, -ISQRT2 * 1j),
                                (2, 3, ISQRT2 * 1j), (3, 2, -ISQRT2 * 1j)]),
        "Sp": build_matrix(3, [(3, 2, SQRT2), (2, 1, SQRT2)]),
        "S+": build_matrix(3, [(3, 2, SQRT2), (2, 1, SQRT2)]),
        "Sm": build_matrix(3, [(1, 2, SQRT2), (2, 3, SQRT2)]),
        "S-": build_matrix(3, [(1, 2, SQRT2), (2, 3, SQRT2)]),
        "Sz2": build_matrix(3, [(1, 1, 1), (3, 3, 1)]),
        "Sx2": build_matrix(3, [(1, 1, 0.5), (1, 3, 0.5), (2, 2, 1.0), (3, 3, 0.5), (3, 1, 0.5)]),
        "Sy2": build_matrix(3, [(1, 1, 0.5), (1, 3, -0.5), (2, 2, 1.0), (3, 3, 0.5), (3, 1, -0.5)]),
        "projUp": build_matrix(3, [(1, 1, 1)]),
        "projZ0": build_matrix(3, [(2, 2, 1)]),
        "projDn": build_matrix(3, [(3, 3, 1)]),
        # ssp1 = 2. for dim==3 (spinone.h's dim(s)==2 ? 0.75 : 2. ternary)
        "S2": build_matrix(3, [(1, 1, 2.0), (3, 3, 2.0), (2, 2, 2.0)]),
    }


class SpinThreeHalfSite(SiteType):
    dim = 4
    _states = {"Up": 1, "3": 1, "Upi": 2, "1": 2, "Dni": 3, "-1": 3, "Dn": 4, "-3": 4}
    _VAL1 = np.sqrt(3.0) / 2.0
    _OPS = {
        "Sz": build_matrix(4, [(1, 1, 1.5), (2, 2, 0.5), (3, 3, -0.5), (4, 4, -1.5)]),
        "Sx": build_matrix(4, [(1, 2, _VAL1), (2, 1, _VAL1), (2, 3, 1.0), (3, 2, 1.0),
                                (3, 4, _VAL1), (4, 3, _VAL1)]),
        "Sy": build_matrix(4, [(1, 2, _VAL1 * 1j), (2, 1, -_VAL1 * 1j), (2, 3, 1.0j), (3, 2, -1.0j),
                                (3, 4, _VAL1 * 1j), (4, 3, -_VAL1 * 1j)]),
    }


class SpinFiveHalfSite(SiteType):
    dim = 6
    _states = {"Up": 1, "5": 1, "Upi": 2, "3": 2, "Upii": 3, "1": 3,
               "Dnii": 4, "-1": 4, "Dni": 5, "-3": 5, "Dn": 6, "-5": 6}
    _VAL1 = np.sqrt(5.0) / 2.0
    _VAL2 = np.sqrt(8.0) / 2.0
    _VAL3 = np.sqrt(9.0) / 2.0
    _OPS = {
        "Sz": build_matrix(6, [(1, 1, 2.5), (2, 2, 1.5), (3, 3, 0.5),
                                (4, 4, -0.5), (5, 5, -1.5), (6, 6, -2.5)]),
        "Sx": build_matrix(6, [(1, 2, _VAL1), (2, 1, _VAL1), (2, 3, _VAL2), (3, 2, _VAL2),
                                (3, 4, _VAL3), (4, 3, _VAL3), (4, 5, _VAL2), (5, 4, _VAL2),
                                (5, 6, _VAL1), (6, 5, _VAL1)]),
        "Sy": build_matrix(6, [(1, 2, _VAL1 * 1j), (2, 1, -_VAL1 * 1j),
                                (2, 3, _VAL2 * 1j), (3, 2, -_VAL2 * 1j),
                                (3, 4, _VAL3 * 1j), (4, 3, -_VAL3 * 1j),
                                (4, 5, _VAL2 * 1j), (5, 4, -_VAL2 * 1j),
                                (5, 6, _VAL1 * 1j), (6, 5, -_VAL1 * 1j)]),
    }


class SpinTwoSite(SiteType):
    dim = 5
    _states = {"Up": 1, "4": 1, "Upi": 2, "2": 2, "Z0": 3, "0": 3,
               "Dni": 4, "-2": 4, "Dn": 5, "-4": 5}
    _VAL1 = np.sqrt(6.0) / 2.0
    _VAL2 = np.sqrt(6.0)
    _OPS = {
        "Sz": build_matrix(5, [(1, 1, 2.0), (2, 2, 1.0), (4, 4, -1.0), (5, 5, -2.0)]),
        "Sx": build_matrix(5, [(1, 2, 1.0), (2, 1, 1.0), (2, 3, _VAL1), (3, 2, _VAL1),
                                (3, 4, _VAL1), (4, 3, _VAL1), (4, 5, 1.0), (5, 4, 1.0)]),
        "ISy": build_matrix(5, [(1, 2, -1.0), (2, 1, 1.0), (2, 3, -_VAL1), (3, 2, _VAL1),
                                 (3, 4, -_VAL1), (4, 3, _VAL1), (4, 5, -1.0), (5, 4, 1.0)]),
        "Sy": build_matrix(5, [(1, 2, 1j), (2, 1, -1j), (2, 3, _VAL1 * 1j), (3, 2, -_VAL1 * 1j),
                                (3, 4, _VAL1 * 1j), (4, 3, -_VAL1 * 1j), (4, 5, 1j), (5, 4, -1j)]),
        "Sp": build_matrix(5, [(2, 1, 2.0), (3, 2, _VAL2), (4, 3, _VAL2), (5, 4, 2.0)]),
        "S+": build_matrix(5, [(2, 1, 2.0), (3, 2, _VAL2), (4, 3, _VAL2), (5, 4, 2.0)]),
        "Sm": build_matrix(5, [(1, 2, 2.0), (2, 3, _VAL2), (3, 4, _VAL2), (4, 5, 2.0)]),
        "S-": build_matrix(5, [(1, 2, 2.0), (2, 3, _VAL2), (3, 4, _VAL2), (4, 5, 2.0)]),
        "Sz2": build_matrix(5, [(1, 1, 4), (2, 2, 1), (4, 4, 1), (5, 5, 4)]),
        "Sx2": build_matrix(5, [(1, 1, 1.0), (1, 3, _VAL1), (2, 2, 2.5), (2, 4, 1.5),
                                 (3, 1, _VAL1), (3, 3, 3.0), (3, 5, _VAL1),
                                 (4, 2, 1.5), (4, 4, 2.5), (5, 3, _VAL1), (5, 5, 1.0)]),
        "Sy2": build_matrix(5, [(1, 1, 1.0), (1, 3, -_VAL1), (2, 2, 2.5), (2, 4, -1.5),
                                 (3, 1, -_VAL1), (3, 3, 3.0), (3, 5, -_VAL1),
                                 (4, 2, -1.5), (4, 4, 2.5), (5, 3, -_VAL1), (5, 5, 1.0)]),
        "projUp": build_matrix(5, [(1, 1, 1)]),
        "projUpi": build_matrix(5, [(2, 2, 1)]),
        "projZ0": build_matrix(5, [(3, 3, 1)]),
        "projDni": build_matrix(5, [(4, 4, 1)]),
        "projDn": build_matrix(5, [(5, 5, 1)]),
        # Transcribed verbatim from itensor/mps/sites/spintwo.h, including
        # its own "Op.set(Dn,DniP,6)" entry (off-diagonal Dn->Dni, not
        # Dn->Dn) -- an apparent typo in vendored upstream ITensor v3 code,
        # not touched here per CLAUDE.md's "preserve pre-existing quirks
        # rather than silently fix them" convention for ported code. dmrgpy
        # itself never requests "S2" (only Sz/Sx/Sy/Sp/Sm reach the Python
        # API, see operatornames.py), so this is inert either way.
        "S2": build_matrix(5, [(1, 1, 6), (2, 2, 6), (3, 3, 6), (4, 4, 6), (5, 4, 6)]),
    }
