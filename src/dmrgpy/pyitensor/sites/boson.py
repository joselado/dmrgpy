"""Boson-4 site (4 occupation levels 0..3), transcribed operator-by-operator
from mpscpp3/extra/bosonfour.h. See sites/base.py for the matrix-entry axis
convention."""

import numpy as np

from .base import SiteType, build_matrix


class BosonFourSite(SiteType):
    dim = 4
    _states = {"3": 1, "Up": 1, "1": 2, "Upi": 2, "-1": 3, "Dni": 3, "-3": 4, "Dn": 4}
    _OPS = {
        "N": build_matrix(4, [(1, 1, 0.0), (2, 2, 1.0), (3, 3, 2.0), (4, 4, 3.0)]),
        "N0": build_matrix(4, [(1, 1, 1.0)]),
        "N1": build_matrix(4, [(2, 2, 1.0)]),
        "N2": build_matrix(4, [(3, 3, 1.0)]),
        "N3": build_matrix(4, [(4, 4, 1.0)]),
        "Adag": build_matrix(4, [(1, 2, 1.0), (2, 3, np.sqrt(2.0)), (3, 4, np.sqrt(3.0))]),
        "A": build_matrix(4, [(2, 1, 1.0), (3, 2, np.sqrt(2.0)), (4, 3, np.sqrt(3.0))]),
    }
