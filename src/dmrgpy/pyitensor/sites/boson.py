"""Boson site(s), transcribed operator-by-operator from
mpscpp3/extra/bosonfour.h -- see sites/base.py for the matrix-entry axis
convention.

Generalized (mirroring the same mpscpp3/get_sites.h generalization) from
a single fixed 4-level site to an arbitrary occupation cutoff: dimension
= maxOcc+1, matching the dmrgpy-level "maxnb" convention
(bosonchain.Bosonic_Chain). `get_boson_site(dim)` builds (and caches) the
SiteType subclass for a given dimension on demand -- see siteset.py's
TYPE_CODE_TO_SITE dispatch, which routes the same 100+dim type-code
range used by mpscpp3's get_sites.h to this factory instead of a single
dict entry. `BosonFourSite` (dim=4, type code 104) is kept as a concrete
name for backwards compatibility, and is simply `get_boson_site(4)`.
"""

import numpy as np

from .base import SiteType, build_matrix


def _boson_ops(dim):
    """N, Adag, A, and the occupation-number projectors N0..N{dim-1},
    generalizing the original hardcoded N/N0-N3/A/Adag entries."""
    maxOcc = dim - 1
    ops = {
        "N": build_matrix(dim, [(k + 1, k + 1, float(k)) for k in range(dim)]),
        "Adag": build_matrix(dim, [(k + 1, k + 2, np.sqrt(k + 1.0)) for k in range(maxOcc)]),
        "A": build_matrix(dim, [(k + 2, k + 1, np.sqrt(k + 1.0)) for k in range(maxOcc)]),
    }
    for k in range(dim):
        ops["N" + str(k)] = build_matrix(dim, [(k + 1, k + 1, 1.0)])
    return ops


_boson_site_cache = {}


def get_boson_site(dim):
    """Return (building and caching on first use) the SiteType subclass
    for a boson site of local Hilbert-space dimension `dim`. Plain
    occupation-number state names only ("0".."{dim-1}") -- unlike the
    original fixed-dim=4 site, no "Up"/"Upi"/"Dni"/"Dn" aliases and no
    Sz-label numeric convention (that scheme doesn't generalize to an
    arbitrary dimension, and nothing in this codebase currently calls
    SiteX.state_index() for a boson site: DMRG here always starts from
    randomMPS(), never a named InitState -- see chain.py/mpsalgebra.py).
    """
    if dim not in _boson_site_cache:
        _boson_site_cache[dim] = type(
            "BosonSite{}".format(dim),
            (SiteType,),
            {
                "dim": dim,
                "_states": {str(k): k + 1 for k in range(dim)},
                "_OPS": _boson_ops(dim),
            },
        )
    return _boson_site_cache[dim]


BosonFourSite = get_boson_site(4)
