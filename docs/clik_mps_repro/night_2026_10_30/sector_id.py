"""Charge-counting sector identity |1_q>> on a doubled chain, and the
left-zero-mode check ||L^dag 1_q||."""
import numpy as np
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.mpscontainer import MPS
from dmrgpy.pyitensor.tensor import ITensor

def sector_identity_mps(sites, q):
    """sum over configurations a with sum_i a_i = q of |a_1 a_1 a_2 a_2 ...>>,
    the charge of a local basis state being its index."""
    n = sites.length(); npair = n // 2
    tensors = []; left = None
    for p in range(npair):
        a, b = 2 * p + 1, 2 * p + 2
        d = sites.dim(a)
        nin = 1 if p == 0 else q + 1
        inner_link = Index(nin * d, tags="Link,l={}".format(a))
        last = (p == npair - 1)
        nout = 1 if last else q + 1
        right = None if last else Index(nout, tags="Link,l={}".format(b))
        A = np.zeros((nin, d, nin * d), complex)
        for k in range(nin):
            for s in range(d):
                A[k, s, k * d + s] = 1.0
        Bt = np.zeros((nin * d, d, nout), complex)
        for k in range(nin):
            cnt = 0 if p == 0 else k
            for s in range(d):
                tot = cnt + s
                if last:
                    if tot == q: Bt[k * d + s, s, 0] = 1.0
                elif tot <= q:
                    Bt[k * d + s, s, tot] = 1.0
        inds_a = ([left] if left else []) + [sites.si(a), inner_link]
        arr_a = A if left else A[0]
        inds_b = [inner_link, sites.si(b)] + ([right] if right else [])
        arr_b = Bt if right else Bt[:, :, 0]
        tensors.append(ITensor(tuple(inds_a), arr_a)); tensors.append(ITensor(tuple(inds_b), arr_b))
        left = right
    return MPS(tensors)
