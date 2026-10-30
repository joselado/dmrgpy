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


def sector_identity_mps(sites, q):
    """|1_q>> = sum over the configurations a with sum_i a_i = q of
    prod_i |a_i>_{2i-1}|a_i>_{2i}, the identity of the sector of total
    charge q, with the charge of a local basis state taken as its index
    (the occupation, for bosons). A charge counter rides on the bonds:
    dimension q+1 between pairs and d(q+1) inside a pair. For a
    Lindbladian that conserves this charge on both sides of rho (a strong
    symmetry), <<1_q| is a left steady state, and it is the trace
    functional of every state of the sector: <<1_q|rho>> = Tr(rho)."""
    n = sites.length()
    if n % 2:
        raise ValueError("sector_identity_mps: a doubled chain has an even "
                         "number of sites, got %d" % n)
    if q < 0 or int(q) != q:
        raise ValueError("sector_identity_mps: the charge must be a "
                         "non-negative integer, got %r" % (q,))
    q = int(q)
    npair = n // 2
    tensors = []
    left = None
    for p in range(npair):
        a, b = 2 * p + 1, 2 * p + 2
        d = sites.dim(a)
        nin = 1 if p == 0 else q + 1
        last = p == npair - 1
        nout = 1 if last else q + 1
        inner_link = Index(nin * d, tags="Link,l={}".format(a))
        right = None if last else Index(nout, tags="Link,l={}".format(b))
        A = np.zeros((nin, d, nin * d), dtype=complex)
        B = np.zeros((nin * d, d, nout), dtype=complex)
        for k in range(nin):
            for s in range(d):
                A[k, s, k * d + s] = 1.0
                total = k + s
                if last and total == q:
                    B[k * d + s, s, 0] = 1.0
                elif not last and total <= q:
                    B[k * d + s, s, total] = 1.0
        inds_a = ([left] if left else []) + [sites.si(a), inner_link]
        inds_b = [inner_link, sites.si(b)] + ([right] if right else [])
        tensors.append(ITensor(tuple(inds_a), A if left else A[0]))
        tensors.append(ITensor(tuple(inds_b), B if right else B[:, :, 0]))
        left = right
    return MPS(tensors)
