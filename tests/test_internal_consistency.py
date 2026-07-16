"""Internal consistency check for multioperator.msum(), which has no
pre-optimization equivalent to compare against (it's a new helper) -
instead it is checked here against the "H = H + term" loop idiom it is
meant to replace, on the current implementation.
"""
from dmrgpy import multioperator
from _helpers import normalize_terms


def test_msum_matches_sequential_accumulation():
    """Summing many single-site/two-site operator terms via msum()
    must give exactly the same (cleaned) set of terms as building the
    same sum one term at a time with "+"."""
    L = 8
    ops = [multioperator.obj2MO([["Sx", i], ["Sx", j]])
           for i in range(L) for j in range(L) if i != j]

    looped = multioperator.zero()
    for op in ops:
        looped = looped + op
    looped.clean()

    bulk = multioperator.msum(ops)
    bulk.clean()

    def as_terms(op):
        # op entries are [coeff, [name,site], [name,site], ...]; adapt
        # to the (coeff, ops) shape normalize_terms() expects.
        return [(complex(t[0]), t[1:]) for t in op]

    assert normalize_terms(as_terms(bulk.op)) == normalize_terms(as_terms(looped.op))
