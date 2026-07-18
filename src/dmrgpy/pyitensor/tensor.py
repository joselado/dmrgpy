"""ITensor: a dense, always-complex128 tensor labeled by Index objects.

This is the numeric half of the engine (index.py is the bookkeeping half).
It implements exactly the operations mpscpp3/chain_session.h calls on plain
ITensors directly (as opposed to through MPS/MPO, see mps.py/mpo.py):
contraction (`*`), addition/subtraction with an ITensor-style "different
index structure" check, scalar multiplication, prime()/noPrime()/mapPrime()/
swapPrime() (by explicit Index or by tag), dag() (conjugation -- no arrows
to flip, since nothing here ever carries a QN flux), and commonIndex().

Every array is complex128 unconditionally. v3 itself is real/complex mixed
(and mpscpp3/chain_session.h has an extended comment on the resulting
"inner() throws on complex" trap it has to route around with innerC()
everywhere) -- since this engine is written fresh rather than mirroring
that C++ type split, always-complex sidesteps the whole issue: there is only
one ITensor type, so nothing here ever needs a plain-vs-complex distinction.
"""

import numpy as np

from .index import Index


def _find(inds, target):
    """Index of the element of `inds` that equals `target` (same id & plev),
    or None. `inds` entries are compared with Index.__eq__, i.e. identity+
    prime level, not object identity -- two different Index *objects* can
    refer to the same leg."""
    for i, ind in enumerate(inds):
        if ind == target:
            return i
    return None


def mul_plan(a_inds, b_inds):
    """The index-matching decision behind ITensor.__mul__ (a*b), factored
    out as a pure function of index *structure* alone (no array data) so
    it can be computed once and reused across many calls with the same
    tensor shapes -- see kernels.py's numba matvec chain, which calls this
    at matvec-build time (once per bond) rather than per Lanczos iteration
    (dozens of times per bond), to plan a fixed sequence of transpose+
    reshape+matmul steps ahead of time.

    Returns (a_axes, b_axes, a_free, b_free): a_axes/b_axes are the
    matched-and-therefore-contracted axis positions (same order, i.e.
    a_axes[k] pairs with b_axes[k]); a_free/b_free are the remaining axis
    positions, in their original order -- exactly what np.tensordot(a, b,
    axes=(a_axes, b_axes)) needs, with the output index order being
    a_inds[a_free] + b_inds[b_free]."""
    a_axes, b_axes, a_free = [], [], []
    used_b = set()
    for i, ind in enumerate(a_inds):
        j = _find(b_inds, ind)
        if j is not None and j not in used_b:
            a_axes.append(i)
            b_axes.append(j)
            used_b.add(j)
        else:
            a_free.append(i)
    b_free = [j for j in range(len(b_inds)) if j not in used_b]
    return a_axes, b_axes, a_free, b_free


class ITensor:
    __slots__ = ("inds", "array")

    def __init__(self, inds, array=None):
        self.inds = tuple(inds)
        shape = tuple(ind.dim for ind in self.inds)
        if array is None:
            self.array = np.zeros(shape, dtype=complex)
        else:
            arr = np.asarray(array, dtype=complex)
            if arr.shape != shape:
                raise ValueError(
                    "ITensor: array shape {} doesn't match indices {}".format(
                        arr.shape, self.inds))
            self.array = arr

    @property
    def rank(self):
        return len(self.inds)

    def hasindex(self, ind):
        return _find(self.inds, ind) is not None

    def scalar(self):
        """The single value of a rank-0 ITensor (mirrors ITensor's eltC/elt),
        used after every index has been contracted away (inner/innerC/traceC/
        bond_entropy's spectrum-free path)."""
        if self.rank != 0:
            raise ValueError("scalar() called on a rank-{} ITensor".format(self.rank))
        return complex(self.array)

    def visit(self, callback):
        """Call `callback(value)` for every element, in row-major order of
        `self.inds` -- mirrors ITensor's ITensor::visit(), used by
        reduced_dm() to flatten a two-leg density matrix into a flat list
        that bindings.cc then reshapes back into a (dim,dim) numpy array."""
        for value in self.array.reshape(-1):
            callback(complex(value))

    def transpose_to(self, order):
        """This tensor's data, reordered so its axes follow `order` (a
        sequence of Index equal to a permutation of self.inds)."""
        perm = [_find(self.inds, ind) for ind in order]
        if any(p is None for p in perm):
            raise ValueError("transpose_to: {} is not a permutation of {}".format(
                order, self.inds))
        return np.transpose(self.array, perm)

    # -- contraction / algebra -------------------------------------------

    def __mul__(self, other):
        if isinstance(other, (int, float, complex, np.number)):
            return ITensor(self.inds, self.array * other)
        if not isinstance(other, ITensor):
            return NotImplemented
        self_axes, other_axes, self_free, other_free = mul_plan(self.inds, other.inds)
        arr = np.tensordot(self.array, other.array, axes=(self_axes, other_axes))
        out_inds = tuple(self.inds[i] for i in self_free) + tuple(other.inds[j] for j in other_free)
        return ITensor(out_inds, arr)

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex, np.number)):
            return ITensor(self.inds, self.array * other)
        return NotImplemented

    def __truediv__(self, scalar):
        return ITensor(self.inds, self.array / scalar)

    def __add__(self, other):
        if not isinstance(other, ITensor):
            return NotImplemented
        if set(self.inds) != set(other.inds):
            raise ValueError(
                "ITensor addition: different index structure {} vs {}".format(
                    self.inds, other.inds))
        arr = self.array + other.transpose_to(self.inds)
        return ITensor(self.inds, arr)

    def __sub__(self, other):
        return self + (-1.0) * other

    def __neg__(self):
        return ITensor(self.inds, -self.array)

    def __repr__(self):
        return "ITensor(inds={})".format(self.inds)


# -- free functions, mirroring ITensor's own free-function API ------------
#
# Kept as free functions (not methods) on purpose: mpscpp3/chain_session.h
# is written entirely in this style (prime(T,...), dag(T), noPrime(T,...),
# commonIndex(A,B)), and later phases port that file close to verbatim --
# matching its calling convention here keeps that port mechanical.

def _prime_inds(inds, filt, inc, op):
    out = []
    for ind in inds:
        match = False
        if not filt:
            match = True
        else:
            for f in filt:
                if isinstance(f, Index):
                    match = (f.id == ind.id)
                else:
                    match = ind.hastags(f)
                if match:
                    break
        out.append(op(ind) if match else ind)
    return tuple(out)


def prime(T, *filt, inc=1):
    """prime(T) primes every index; prime(T,"Site")/prime(T,"Link") primes
    only indices carrying that tag; prime(T, ind1, ind2, ...) primes only
    the specific indices in T matching those (by id, regardless of the
    passed-in Index's own prime level) -- e.g. reduced_dm()'s
    prime(psi.A(site), sites_.si(site), ir)."""
    return ITensor(_prime_inds(T.inds, filt, inc, lambda i: i.prime(inc)), T.array)


def noPrime(T, *filt):
    return ITensor(_prime_inds(T.inds, filt, 0, lambda i: i.setprime(0)), T.array)


def mapPrime(T, a, b, *filt):
    def op(ind):
        return ind.setprime(b) if ind.plev == a else ind
    new_inds = []
    for ind in T.inds:
        selected = (not filt) or any(
            (ind.id == f.id) if isinstance(f, Index) else ind.hastags(f) for f in filt)
        new_inds.append(op(ind) if selected else ind)
    return ITensor(tuple(new_inds), T.array)


def swapPrime(T, a, b, *filt):
    def op(ind):
        if ind.plev == a:
            return ind.setprime(b)
        if ind.plev == b:
            return ind.setprime(a)
        return ind
    new_inds = []
    for ind in T.inds:
        selected = (not filt) or any(
            (ind.id == f.id) if isinstance(f, Index) else ind.hastags(f) for f in filt)
        new_inds.append(op(ind) if selected else ind)
    return ITensor(tuple(new_inds), T.array)


def dag(T):
    """Complex-conjugate every element. Real ITensor also reverses Index
    arrows; there are none here (see module docstring), so conjugation is
    the whole of it."""
    return ITensor(T.inds, np.conj(T.array))


def commonIndex(A, B, tagmatch=None):
    """The Index shared between A and B (by id & plev), optionally filtered
    by tag; None if there isn't one. Mirrors ITensor's commonIndex(), used
    by reduced_dm() to find the bond between two neighboring MPS tensors."""
    bset = {ind: ind for ind in B.inds}
    for ind in A.inds:
        if ind in bset and ind.hastags(tagmatch):
            return ind
    return None


def delta(i, j):
    """Identity/Kronecker-delta ITensor over two same-dimension indices --
    used to build trivial link tensors (e.g. edges of an MPS/MPO chain)."""
    if i.dim != j.dim:
        raise ValueError("delta: mismatched dimensions {} vs {}".format(i.dim, j.dim))
    return ITensor((i, j), np.eye(i.dim, dtype=complex))


def _pairwise_result_size(a, b):
    """Element count of a*b, without actually contracting: every Index
    shared between a and b (by identity+plev, i.e. Index equality) drops
    out; everything else survives and multiplies into the result size."""
    shared = set(a.inds) & set(b.inds)
    size = 1
    for ind in a.inds:
        if ind not in shared:
            size *= ind.dim
    for ind in b.inds:
        if ind not in shared:
            size *= ind.dim
    return size


def contract_many(tensors):
    """Contract a list of ITensors together, choosing a greedy
    size-minimizing pairwise order (repeatedly contract whichever
    remaining pair would produce the *smallest* result, per
    _pairwise_result_size) instead of Python's default left-to-right `*`
    chaining.

    This exists because that naive chaining is a real, previously-shipped
    bug, not a hypothetical one: dmrg.py's environment extension and
    mpsalgebra.py's inner() both used to write `piece = a * b * c` where
    the *last* operand was an already-accumulated environment -- meaning
    the first pairwise step (a*b) left BOTH a's and b's own link legs
    dangling simultaneously, before ever touching the environment that
    would cancel one side of them away. Measured directly on a 14-site,
    maxdim=60 DMRG bond: a 92-million-element intermediate tensor, 1.5s
    for one call, that a differently-ordered but mathematically identical
    contraction computes in 0.0006s (~2500x). Routing every multi-operand
    contraction through this function instead of hand-chaining `*` is
    meant to make that whole class of bug structurally hard to
    reintroduce, not just patch the two call sites it was found in.

    Cost of the greedy search itself is O(n^3) in the number of input
    tensors, negligible for the n ~ 3-5 operand contractions this package
    actually does (a two-site effective-Hamiltonian piece, an environment
    extension, an inner() term) -- this is not meant for contracting large
    tensor networks with many operands, where a proper contraction-order
    solver (as opposed to a greedy one) would matter.
    """
    remaining = list(tensors)
    if not remaining:
        raise ValueError("contract_many: no tensors given")
    while len(remaining) > 1:
        best = None
        for a_idx in range(len(remaining)):
            for b_idx in range(a_idx + 1, len(remaining)):
                size = _pairwise_result_size(remaining[a_idx], remaining[b_idx])
                if best is None or size < best[0]:
                    best = (size, a_idx, b_idx)
        _, a_idx, b_idx = best
        b = remaining.pop(b_idx)  # pop the larger index first so a_idx stays valid
        a = remaining.pop(a_idx)
        remaining.append(a * b)
    return remaining[0]
