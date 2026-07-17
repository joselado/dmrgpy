"""Zn parafermion site types, transcribed operator-by-operator from the
stock mpscpp3/ITensor/itensor/mps/sites/Z3.h and dmrgpy's own
mpscpp3/extra/Z4.h. See sites/base.py for the matrix-entry axis convention.
"""

import numpy as np

from .base import SiteType, build_matrix


class Z3Site(SiteType):
    dim = 3
    _states = {"0": 1, "1": 2, "2": 3}
    _OPS = {
        "N": build_matrix(3, [(2, 2, 1), (3, 3, 2)]),
        "Sig": build_matrix(3, [(1, 3, 1), (2, 1, 1), (3, 2, 1)]),
        "SigDag": build_matrix(3, [(3, 1, 1), (1, 2, 1), (2, 3, 1)]),
        "Tau": build_matrix(3, [
            (1, 1, 1),
            (2, 2, np.cos(2 * np.pi / 3) + 1j * np.sin(2 * np.pi / 3)),
            (3, 3, np.cos(4 * np.pi / 3) + 1j * np.sin(4 * np.pi / 3)),
        ]),
        "TauDag": build_matrix(3, [
            (1, 1, 1),
            (2, 2, np.cos(2 * np.pi / 3) - 1j * np.sin(2 * np.pi / 3)),
            (3, 3, np.cos(4 * np.pi / 3) - 1j * np.sin(4 * np.pi / 3)),
        ]),
        "Proj0": build_matrix(3, [(1, 1, 1)]),
        "Proj1": build_matrix(3, [(2, 2, 1)]),
        "Proj2": build_matrix(3, [(3, 3, 1)]),
    }


class Z4Site(SiteType):
    dim = 4
    _states = {"0": 1, "1": 2, "2": 3, "3": 4}
    _OPS = {
        "N": build_matrix(4, [(2, 2, 1), (3, 3, 2), (4, 4, 3)]),
        "Sig": build_matrix(4, [(1, 4, 1), (2, 1, 1), (3, 2, 1), (4, 3, 1)]),
        "SigDag": build_matrix(4, [(4, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)]),
        "Tau": build_matrix(4, [
            (1, 1, 1),
            (2, 2, np.cos(2 * np.pi / 4) + 1j * np.sin(2 * np.pi / 4)),
            (3, 3, np.cos(4 * np.pi / 4) + 1j * np.sin(4 * np.pi / 4)),
            (4, 4, np.cos(6 * np.pi / 4) + 1j * np.sin(6 * np.pi / 4)),
        ]),
        "TauDag": build_matrix(4, [
            (1, 1, 1),
            (2, 2, np.cos(2 * np.pi / 4) - 1j * np.sin(2 * np.pi / 4)),
            (3, 3, np.cos(4 * np.pi / 4) - 1j * np.sin(4 * np.pi / 4)),
            (4, 4, np.cos(6 * np.pi / 4) - 1j * np.sin(6 * np.pi / 4)),
        ]),
        # Z4.h only defines Proj0/Proj1/Proj2 (no Proj3) -- transcribed as-is.
        "Proj0": build_matrix(4, [(1, 1, 1)]),
        "Proj1": build_matrix(4, [(2, 2, 1)]),
        "Proj2": build_matrix(4, [(3, 3, 1)]),
    }
