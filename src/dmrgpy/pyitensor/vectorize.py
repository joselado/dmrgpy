"""Density matrices as MPS on a doubled chain.

A density matrix of an n-site chain, rho = sum rho_ab |a><b|, is the state
|rho>> = sum rho_ab |a>|b> of a 2n-site chain whose physical and auxiliary
sites alternate: site 2i-1 carries the ket index of site i and site 2i its
bra index (1-based, as in this package). The trace is the overlap with the
vectorized identity, Tr(rho) = <<1|rho>>, and Tr(O rho) = <<1|O x 1|rho>>.
"""

import numpy as np

from .index import Index
from .mpscontainer import MPS
from .tensor import ITensor


def identity_mps(sites):
    """|1>> = prod_i sum_a |a>_{2i-1}|a>_{2i}: bond dimension d inside each
    pair and 1 between pairs. Not normalized: <<1|1>> is the dimension of
    the physical Hilbert space."""
    n = sites.length()
    if n % 2:
        raise ValueError("identity_mps: a doubled chain has an even number "
                         "of sites, got %d" % n)
    tensors = []
    left = None  # dim-1 link from the previous pair
    for p in range(n // 2):
        a, b = 2 * p + 1, 2 * p + 2
        d = sites.dim(a)
        if sites.dim(b) != d:
            raise ValueError("identity_mps: sites %d and %d of a pair have "
                             "dimensions %d and %d" % (a, b, d, sites.dim(b)))
        inner_link = Index(d, tags="Link,l={}".format(a))
        right = Index(1, tags="Link,l={}".format(b)) if b < n else None
        eye = np.eye(d, dtype=complex)
        inds_a = ([left] if left else []) + [sites.si(a), inner_link]
        arr_a = eye.reshape(((1,) if left else ()) + (d, d))
        inds_b = [inner_link, sites.si(b)] + ([right] if right else [])
        arr_b = eye.reshape((d, d) + ((1,) if right else ()))
        tensors.append(ITensor(tuple(inds_a), arr_a))
        tensors.append(ITensor(tuple(inds_b), arr_b))
        left = right
    return MPS(tensors)
