"""Index: the identity/bookkeeping half of the tensor engine.

Mirrors the subset of ITensor v3's Index type that dmrgpy's mpscpp3 backend
actually uses (see mpscpp3/chain_session.h): a stable identity ("this leg,
regardless of prime level"), a prime level (distinguishes a "bra" copy of a
leg from its "ket" original after prime()/dag()), a dimension, and a small
set of string tags ("Site", "Link", plus free-form labels like "n=3" for
readability), plus -- only in conserved-sector mode -- a *charge grading*:
one integer charge per basis state, per conserved quantity.

That grading is the whole of the quantum-number bookkeeping here, and it
is deliberately much less than real ITensor v3's. A QN-carrying v3 Index
also reorders its basis states into contiguous per-charge blocks and makes
every tensor built on it block-sparse; this one only *labels* the states it
already had, leaving storage dense (see sector.py's module docstring for
what that buys and what it costs). Two consequences worth stating plainly:
the labels never affect Index equality or contraction (a graded Index and
its dense twin differ only by identity, exactly as a v3 QN site index does),
and only *site* indices are ever graded -- Link indices minted by svd()/
sim() stay ungraded, since nothing here exploits their block structure.

Two Index objects are the "same leg" for contraction/addition purposes iff
they compare equal, which requires both the same identity and the same
prime level -- exactly ITensor's own rule (priming turns a leg into a
distinguishable one on purpose, e.g. bra vs ket copies of a physical index).
"""

import itertools

_id_counter = itertools.count()


def _parse_tags(tags):
    if tags is None:
        return frozenset()
    if isinstance(tags, str):
        return frozenset(t for t in tags.split(",") if t)
    return frozenset(tags)


class Index:
    __slots__ = ("_id", "_dim", "_tags", "_plev", "_qnnames", "_qn")

    def __init__(self, dim, tags=(), plev=0, _id=None, qnnames=None, qn=None):
        self._id = next(_id_counter) if _id is None else _id
        self._dim = int(dim)
        self._tags = _parse_tags(tags)
        self._plev = plev
        # Charge grading, conserved-sector mode only (None when off): a
        # tuple of quantum-number names, and one integer charge per name
        # per basis state -- qn[state][k] is the charge of basis state
        # `state` (0-based) under quantum number qnnames[k].
        self._qnnames = None if qnnames is None else tuple(qnnames)
        if qn is None:
            self._qn = None
        else:
            self._qn = tuple(tuple(int(c) for c in row) for row in qn)
            if len(self._qn) != self._dim:
                raise ValueError("Index: %d charge rows for a dim-%d index"
                                 % (len(self._qn), self._dim))
            nq = len(self._qnnames or ())
            if any(len(row) != nq for row in self._qn):
                raise ValueError("Index: charge rows don't match qnnames %s"
                                 % (self._qnnames,))

    @property
    def id(self):
        return self._id

    @property
    def dim(self):
        return self._dim

    @property
    def tags(self):
        return self._tags

    @property
    def plev(self):
        return self._plev

    @property
    def qnnames(self):
        """The conserved quantities this index is graded by, or None."""
        return self._qnnames

    @property
    def qn(self):
        """Per-basis-state charges (a tuple of tuples), or None when this
        index carries no grading."""
        return self._qn

    def hasqns(self):
        return self._qn is not None

    def charge(self, state):
        """The charge tuple of (0-based) basis state `state`."""
        return self._qn[state]

    def hastags(self, tagmatch):
        """True if this Index carries every tag in tagmatch (a str or iterable
        of str). None/empty matches everything -- the "no filter" case used
        throughout chain_session.h's prime(T)/noPrime(T) with no tag arg."""
        want = _parse_tags(tagmatch)
        return want <= self._tags

    def prime(self, inc=1):
        return Index(self._dim, self._tags, self._plev + inc, _id=self._id,
                     qnnames=self._qnnames, qn=self._qn)

    def setprime(self, plev):
        return Index(self._dim, self._tags, plev, _id=self._id,
                     qnnames=self._qnnames, qn=self._qn)

    def noprime(self):
        return self.setprime(0)

    def sim(self):
        """A fresh Index with the same dim/tags but a brand-new identity --
        unrelated to this one for contraction/equality purposes. Used
        whenever an algorithm must mint a new Link index (SVD splitting,
        building a product-state MPS, ...), mirroring ITensor's sim(Index)."""
        return Index(self._dim, self._tags)

    def __eq__(self, other):
        return isinstance(other, Index) and self._id == other._id and self._plev == other._plev

    def __hash__(self):
        return hash((self._id, self._plev))

    def __repr__(self):
        tagstr = ",".join(sorted(self._tags))
        qnstr = "" if self._qn is None else ",qns='{}'".format(
            ",".join(self._qnnames))
        return "Index(dim={},tags='{}',plev={},id={}{})".format(
            self._dim, tagstr, self._plev, self._id, qnstr)


def sim(index):
    return index.sim()


def reseed_id_counter_past(max_seen_id):
    """Advance the process-local id counter so the next Index() minted in
    *this* process gets an id strictly greater than max_seen_id.

    Exists for pyitensor/metts.py's njobs>1 path: a spawned
    multiprocessing worker starts with a fresh _id_counter at 0 (a brand
    new interpreter, not a fork of the parent's memory), but the H/sites/
    ops MPS/MPO/Index objects it receives were pickled in the *parent*
    process and carry whatever ids _that_ process's counter had already
    assigned -- pickling an Index round-trips its numeric _id verbatim
    (plain __slots__ pickling, no custom __reduce__), it doesn't route
    through Index.__init__/_id_counter at all. Since Index.__eq__/__hash__
    compare _id (and _plev) alone, with no dim/tags check at all (see this
    module's docstring), a freshly-started worker mid-run naturally mints
    ids 0, 1, 2, ... internally (e.g. every new MPS bond-dimension-1 link
    Index built while sampling a classical product state) that can collide
    with an unrelated, differently-shaped Index unpickled from the parent
    -- confirmed directly: a short-chain metts_vev(..., njobs=2) call
    reliably hit "ValueError: shape-mismatch for sum" deep inside a worker
    process's tensor contraction, because two structurally different
    indices with the same colliding id compared equal. Calling this once
    per worker, before any Index() this process creates, closes that
    window for ids already present in the unpickled task -- it cannot
    protect against ids assigned *after* this call by anything running
    concurrently in the same process, but METTS workers are single-
    threaded, so there is no such concurrent user here."""
    global _id_counter
    _id_counter = itertools.count(max_seen_id + 1)
