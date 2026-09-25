"""Regression for the `"python"` MPO builder's small-units floor, found while
taking the baseline of the 2026-09-25b fix pass (docs/audit_2026_09_25b_hole_hunt.md,
"Found during the fix pass").

pyitensor's to_mpo() built every operator at the caller's units. The
automaton's identity channel carries O(1) entries whatever the coefficients,
and the two sweeps leave roundoff of order 1e-16 in them, an absolute error:
a Heisenberg chain written at s = 1e-7 came out 1e-9 off relative to itself
and at s = 2^-40 2.5e-4 off, on numpy with OpenBLAS (MKL happened to round
those entries to exact zero, which is why test_audit_2026_09_24c_pyitensor's
1e-12 tolerance held where it was written). to_mpo now builds at unit scale,
as v2/v3's to_mpo_unit() does. The property pinned is scale covariance: the
relative error at s is the relative error at s = 1, and the scaling is by a
power of two, so s = 2^-40 is exactly representable and nothing but the
builder can move the result.
"""

import numpy as np
import pytest

from dmrgpy.pyitensor import mpobuilder as mb
from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyitensor.tensor import contract_many


def _dense(mpo, sites):
    n = mpo.length()
    T = contract_many([mpo.A(i) for i in range(1, n + 1)])
    si = [sites.si(i) for i in range(1, n + 1)]
    arr = np.asarray(T.transpose_to([i.prime(1) for i in si] + si))
    dim = int(np.prod([i.dim for i in si]))
    return arr.reshape(dim, dim)


def _heis(n, s):
    return [(s, [(op, i), (op, i + 1)]) for i in range(1, n) for op in ("Sx", "Sy", "Sz")]


def _relerr(n, terms):
    sites = SiteX([2]*n)
    a = AutoMPO.from_terms(sites, terms)
    ref = a.dense_matrix()
    m = mb.to_mpo(a, cutoff=1e-14)
    return np.linalg.norm(_dense(m, sites) - ref)/np.linalg.norm(ref), m


@pytest.mark.parametrize("n", [6, 8])
@pytest.mark.parametrize("s", [1e-7, 1e-10, 2.0**-40])
def test_small_units_error_is_the_unit_scale_error(n, s):
    err1, _ = _relerr(n, _heis(n, 1.0))
    err, m = _relerr(n, _heis(n, s))
    assert err < max(10*err1, 1e-14)
    assert max(m.A(i).inds[-1].dim for i in range(1, n)) == 5


def test_unit_scale_is_untouched():
    """At a largest coefficient of 1 or more the builder takes the unscaled
    path, so the tensors are the ones it always built."""
    assert mb._unit_scale_up(1.0) == 1.0
    assert mb._unit_scale_up(3.7) == 1.0
    assert mb._unit_scale_up(0.0) == 1.0
    for c in (0.9, 1e-7, 3e-300):
        up = mb._unit_scale_up(c)
        assert 1.0 <= c*up < 2.0
        assert np.log2(up) == int(np.log2(up))
