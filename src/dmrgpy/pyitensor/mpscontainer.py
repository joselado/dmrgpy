"""MPS/MPO: chains of ITensor site tensors, plus orthogonality-center
canonicalization.

Boundary sites have no trivial dim-1 edge index (unlike real ITensor, which
always carries one): site 1 has no left Link, site N has no right Link.
This avoids ever having to special-case a "trivial" index and keeps every
downstream contraction (inner/traceC/...) fall out of the generic
"multiply everything, matching indices auto-contract" pattern with no
leftover rank-2-of-dim-1 tensors to squeeze away.

Canonicalization (position()) doesn't care whether a site tensor carries
one "Site"-tagged leg (MPS) or two (MPO): a bond's two sides are always
identified purely by *Link*-tag identity against the neighboring tensor
(see _shift_right/_shift_left), so MPS and MPO share one implementation
via the common _Chain base -- mirrors how mpscpp3/chain_session.h itself
never needs a separate code path for canonicalizing an MPS vs an MPO.
"""

from .index import Index
from .svd import svd
from .tensor import ITensor, commonIndex
from .tensor import mapPrime as _t_mapPrime
from .tensor import noPrime as _t_noPrime
from .tensor import prime as _t_prime


def _link_at(chain, i, j):
    """The Link index shared between chain's tensors at (1-based) sites i
    and j, or None if j is out of range (i.e. i is a boundary site on that
    side)."""
    if j < 1 or j > chain.length():
        return None
    return commonIndex(chain.A(i), chain.A(j))


class _Chain:
    def __init__(self, tensors):
        self._tensors = list(tensors)
        self.center = None  # 1-based orthogonality center, or None if unknown

    def length(self):
        return len(self._tensors)

    def A(self, i):
        return self._tensors[i - 1]

    def Aref(self, i):
        return self._tensors[i - 1]

    def set_A(self, i, tensor):
        self._tensors[i - 1] = tensor

    def copy(self):
        """Shallow copy: a new chain sharing the same ITensor objects.
        Safe because every operation in this engine returns a *new*
        ITensor rather than mutating one in place (see tensor.py) -- the
        same convention mps.py/staticoperator.py already rely on for their
        C++-handle objects' own copy()."""
        out = self.__class__.__new__(self.__class__)
        out._tensors = list(self._tensors)
        out.center = self.center
        return out

    def noPrime(self, *filt):
        for i in range(1, self.length() + 1):
            self.set_A(i, _t_noPrime(self.A(i), *filt))
        return self

    def mapPrime(self, a, b, *filt):
        for i in range(1, self.length() + 1):
            self.set_A(i, _t_mapPrime(self.A(i), a, b, *filt))
        return self

    def prime(self, *filt, inc=1):
        for i in range(1, self.length() + 1):
            self.set_A(i, _t_prime(self.A(i), *filt, inc=inc))
        return self

    # -- canonicalization --------------------------------------------

    def position(self, b, cutoff=0.0, maxdim=None):
        """Move the orthogonality center to (1-based) site b, one bond at a
        time. Requires self.center to already be set (every constructor in
        this package -- randomMPS, sum, applyMPO, ... -- sets it). Lossless
        (cutoff=0, maxdim=None) by default; a truncating sweep (used by
        sum()/applyMPO()/nmultMPO() to compress after a bond-dimension-
        growing construction) is the same procedure with cutoff/maxdim
        actually enforced -- correct regardless of whether the chain was
        already canonical anywhere to begin with, see _shift_right's
        docstring."""
        if self.center is None:
            raise RuntimeError("position(): orthogonality center is unknown")
        while self.center < b:
            self._shift_right(self.center, cutoff, maxdim)
            self.center += 1
        while self.center > b:
            self._shift_left(self.center, cutoff, maxdim)
            self.center -= 1
        return self

    def _shift_right(self, i, cutoff=0.0, maxdim=None):
        """Absorb site i into site i+1 via SVD, leaving a left-orthogonal
        tensor at site i. Doesn't assume site i was already left-orthogonal
        -- an SVD split is exact (up to cutoff/maxdim) regardless of the
        tensor's prior gauge, which is what makes a single sweep from a
        totally arbitrary (e.g. freshly summed or freshly random) chain a
        correct canonicalization, not just a canonical-to-canonical move."""
        T = self.A(i)
        nxt = self.A(i + 1)
        right_link = _link_at(self, i, i + 1)
        left_inds = [ind for ind in T.inds if ind != right_link]
        U, S, V, spec = svd(T, left_inds, cutoff=cutoff, maxdim=maxdim)
        self.set_A(i, U)
        SV = S * V
        self.set_A(i + 1, SV * nxt)

    def _shift_left(self, i, cutoff=0.0, maxdim=None):
        """Mirror of _shift_right: absorb site i into site i-1, leaving a
        right-orthogonal tensor at site i."""
        T = self.A(i)
        prev = self.A(i - 1)
        left_link = _link_at(self, i, i - 1)
        U, S, V, spec = svd(T, [left_link], cutoff=cutoff, maxdim=maxdim)
        self.set_A(i, V)
        US = U * S
        self.set_A(i - 1, prev * US)


class MPS(_Chain):
    def normalize(self):
        # Scaling *any single* tensor of the chain by a scalar scales the
        # whole represented state by that scalar, regardless of gauge --
        # unlike position()'s SVD moves, this needs no canonical form.
        from .mpsalgebra import inner
        nrm = inner(self, self).real ** 0.5
        c = self.center or 1
        self.set_A(c, self.A(c) * (1.0 / nrm))
        return self

    def __imul__(self, scalar):
        c = self.center or 1
        self.set_A(c, self.A(c) * scalar)
        return self

    def __mul__(self, scalar):
        out = self.copy()
        out *= scalar
        return out

    __rmul__ = __mul__


class MPO(_Chain):
    def __mul__(self, scalar):
        out = self.copy()
        c = out.center or 1
        out.set_A(c, out.A(c) * scalar)
        return out

    __rmul__ = __mul__
