"""pyitensor's density-matrix MPS sum, sum_many(method="densitymatrix")
(mpsalgebra._sum_densitymatrix), against the exact route (the direct sum
and one truncating sweep) and against the sum of the dense vectors.

The density-matrix route keeps at every bond the leading eigenvectors of
the reduced density matrix of the sum, cross terms included, projected on
the basis already kept on the left. These are the left Schmidt vectors that
the truncating sweep of the exact route keeps at the same bond, so the two
methods are the same truncation: the same state, an error bounded by the
Schmidt tails of the exact sum, and a result that is an orthogonal
projection of it.

The chains have 12 sites or fewer, so every state fits in a dense vector.
That vector is the reference, and it shares no code with either method.
"""
import numpy as np
import pytest

from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor import mpobuilder as mb
from dmrgpy.pyitensor import mpsalgebra as ma
from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyitensor.tensor import commonIndex, dag, prime

# Complex weights of different size, so that the inputs differ in norm
# and in phase and a missing conjugation cannot cancel out.
COEFFS = (1.0, 0.7 - 0.4j, -0.3 + 0.9j)


def _inputs(n, chis, seed=3):
    """Random complex MPS on one spin-1/2 site set, one per bond dimension
    in `chis`, scaled by the weights in COEFFS."""
    sites = SiteX([2] * n)
    rng = np.random.default_rng(seed)
    states = [ma.randomMPS(sites, chi, rng=rng) * c
              for chi, c in zip(chis, COEFFS)]
    return sites, states


def _vec(mps, sites):
    """An MPS as a dense vector, site 1 most significant."""
    T = mps.A(1)
    for i in range(2, mps.length() + 1):
        T = T * mps.A(i)
    arr = T.transpose_to([sites.si(i) for i in range(1, mps.length() + 1)])
    return np.asarray(arr).reshape(-1)


def _op(mpo, sites):
    """An MPO as a dense matrix, out legs (primed) as rows."""
    n = mpo.length()
    T = mpo.A(1)
    for i in range(2, n + 1):
        T = T * mpo.A(i)
    order = ([sites.si(i).prime(1) for i in range(1, n + 1)]
             + [sites.si(i) for i in range(1, n + 1)])
    return np.asarray(T.transpose_to(order)).reshape(2**n, 2**n)


def _dense_sum(states, sites):
    out = 0.0
    for s in states:
        out = out + _vec(s, sites)
    return out


def _rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def _schmidt_values(vec, n):
    """The Schmidt values of a dense n-site spin-1/2 vector at every bond."""
    return [np.linalg.svd(vec.reshape(2**b, 2**(n - b)), compute_uv=False)
            for b in range(1, n)]


def _bond_dims(mps):
    return [commonIndex(mps.A(i), mps.A(i + 1)).dim
            for i in range(1, mps.length())]


def _assert_center_last(out):
    """center == n means sites 1..n-1 are left-orthogonal, so the whole
    norm sits on site n."""
    n = out.length()
    assert out.center == n
    for j in range(1, n):
        T = out.A(j)
        right = commonIndex(T, out.A(j + 1))
        G = T * dag(prime(T, right))
        g = np.asarray(G.transpose_to([right, right.prime(1)]))
        assert np.max(np.abs(g - np.eye(right.dim))) < 1e-10


CASES = {"two-8": (8, (4, 10)), "two-12": (12, (6, 9)),
         "three-10": (10, (4, 7, 10)), "three-12": (12, (5, 8, 6))}


@pytest.mark.parametrize("case", sorted(CASES))
def test_same_state_as_the_exact_sum_without_truncation(case):
    """With cutoff 0 and no maxdim both methods return the sum itself. The
    random inputs have no Schmidt value near the precision floor of the
    density-matrix route (about 1e-8 of the largest), so 1e-10 is a fair
    tolerance. The inputs are not modified."""
    n, chis = CASES[case]
    sites, states = _inputs(n, chis)
    before = [_vec(s, sites) for s in states]
    ref = _dense_sum(states, sites)
    dm = ma.sum_many(states, method="densitymatrix")
    ex = ma.sum_many(states, method="exact")
    assert _rel(_vec(dm, sites), ref) < 1e-10
    assert _rel(_vec(ex, sites), ref) < 1e-10
    assert _rel(_vec(dm, sites), _vec(ex, sites)) < 1e-10
    _assert_center_last(dm)
    for s, v in zip(states, before):
        assert _rel(_vec(s, sites), v) < 1e-14
    if len(states) == 2:
        two = ma.sum(states[0], states[1], method="densitymatrix")
        assert _rel(_vec(two, sites), ref) < 1e-10


TRUNCATIONS = {"maxdim": dict(cutoff=0.0, maxdim=5),
               "cutoff": dict(cutoff=1e-3, maxdim=None)}
TRUNC_CASES = {"two": (10, (6, 8)), "three": (10, (4, 6, 5))}


@pytest.mark.parametrize("trunc", sorted(TRUNCATIONS))
@pytest.mark.parametrize("case", sorted(TRUNC_CASES))
def test_truncation_error_is_bounded_by_the_discarded_weight(case, trunc):
    """Both methods discard, bond by bond from the left, the smallest
    Schmidt weights of a sum that is already projected on the left, so the
    squared error is at most the sum over bonds of what each bond must
    lose: the Schmidt tail of the exact sum beyond maxdim, or `cutoff` of
    the norm squared. With maxdim it is also at least the largest single
    tail, which no state of that bond dimension can beat, so the check is
    not vacuous. The result is P times the sum with P an orthogonal
    projector, so its overlap with the sum is real and equal to its norm
    squared. And since both keep the same basis at every bond, they give
    the same state, as long as no Schmidt value is degenerate at a cut."""
    kw = TRUNCATIONS[trunc]
    n, chis = TRUNC_CASES[case]
    sites, states = _inputs(n, chis, seed=11)
    ref = _dense_sum(states, sites)
    norm2 = np.vdot(ref, ref).real
    if kw["maxdim"] is None:
        upper, lower = (n - 1) * kw["cutoff"] * norm2, 0.0
    else:
        m = kw["maxdim"]
        svals = _schmidt_values(ref, n)
        tails = np.array([np.sum(s[m:]**2) for s in svals])
        upper, lower = tails.sum(), tails.max()
        assert lower > 1e-6 * norm2, "maxdim does not bind in this case"
        for s in svals:
            if len(s) > m:
                assert s[m - 1] - s[m] > 1e-6 * s[0], \
                    "near-degenerate Schmidt pair at a cut, pick another seed"
    vecs = {}
    for method in ("densitymatrix", "exact"):
        out = ma.sum_many(states, method=method, **kw)
        v = _vec(out, sites)
        err2 = np.vdot(ref - v, ref - v).real
        assert lower * (1 - 1e-8) <= err2 <= upper * (1 + 1e-8) + 1e-14 * norm2
        if kw["maxdim"] is not None:
            assert max(_bond_dims(out)) <= kw["maxdim"]
        out2 = np.vdot(v, v).real
        assert abs(np.vdot(ref, v) - out2) < 1e-10 * norm2
        assert ma.inner(out, out).real == pytest.approx(out2, rel=1e-10)
        if method == "densitymatrix":
            _assert_center_last(out)
            last = np.asarray(out.A(n).array)
            assert np.vdot(last, last).real == pytest.approx(out2, rel=1e-10)
        vecs[method] = v
    assert _rel(vecs["densitymatrix"], vecs["exact"]) < 1e-8


def test_inputs_that_share_links():
    """A state and a scaled copy of it share every link, which must not be
    contracted between them. A state minus itself is the zero state, not
    NaN, and a state plus a multiple of itself truncated to its own bond
    dimension is still exact."""
    sites, states = _inputs(8, (6,))
    psi = states[0]
    v = _vec(psi, sites)
    two = ma.sum_many([psi, psi], method="densitymatrix")
    assert _rel(_vec(two, sites), 2.0 * v) < 1e-10
    w = 0.5j
    mix = ma.sum_many([psi, psi * w], cutoff=0.0, maxdim=6,
                      method="densitymatrix")
    assert _rel(_vec(mix, sites), (1.0 + w) * v) < 1e-10
    zero = _vec(ma.sum_many([psi, psi * (-1.0)], method="densitymatrix"), sites)
    assert np.all(np.isfinite(zero))
    assert np.linalg.norm(zero) < 1e-10 * np.linalg.norm(v)


def test_one_site_chains():
    """A one-site sum has no bond to truncate: it is the sum of the site
    tensors."""
    sites, states = _inputs(1, (1, 1, 1))
    out = ma.sum_many(states, maxdim=1, method="densitymatrix")
    assert _rel(_vec(out, sites), _dense_sum(states, sites)) < 1e-12
    assert out.center == 1


def test_default_and_unknown_methods():
    """The default is the density-matrix route (the author's pick of
    2026-11-15) and equals the exact one here, and an unknown name is an error
    rather than a fallback, also for a single chain, which is returned
    before any work is done."""
    assert ma.SUM_METHOD == "densitymatrix"  # the default since 2026-11-15
    sites, states = _inputs(6, (3, 4))
    a = _vec(ma.sum_many(states, maxdim=4), sites)
    b = _vec(ma.sum_many(states, maxdim=4, method="exact"), sites)
    assert _rel(a, b) < 1e-12
    c = _vec(ma.sum(states[0], states[1], maxdim=4), sites)
    assert _rel(c, b) < 1e-12
    with pytest.raises(ValueError):
        ma.sum_many(states, method="zipup")
    with pytest.raises(ValueError):
        ma.sum(states[0], states[1], method="zipup")
    with pytest.raises(ValueError):
        ma.sum_many(states[:1], method="zipup")


def test_module_default_is_read_at_call_time(monkeypatch):
    """SUM_METHOD is read when sum_many runs, and a sum of MPO keeps the
    exact route even when the density-matrix route is asked for."""
    calls = []
    real = ma._sum_densitymatrix

    def spy(chains, cutoff=0.0, maxdim=None):
        calls.append(len(chains))
        return real(chains, cutoff=cutoff, maxdim=maxdim)

    monkeypatch.setattr(ma, "_sum_densitymatrix", spy)
    monkeypatch.setattr(ma, "SUM_METHOD", "densitymatrix")
    sites, states = _inputs(6, (3, 4))
    ref = _dense_sum(states, sites)
    assert _rel(_vec(ma.sum_many(states), sites), ref) < 1e-10
    assert _rel(_vec(ma.sum(states[0], states[1]), sites), ref) < 1e-10
    assert calls == [2, 2]

    heisenberg = [(1.0, [(op, i), (op, i + 1)])
                  for i in range(1, 6) for op in ("Sx", "Sy", "Sz")]
    field = [(0.3 + 0.1j, [("Sz", i)]) for i in range(1, 7)]
    H1 = mb.to_mpo(AutoMPO.from_terms(sites, heisenberg), cutoff=1e-14)
    H2 = mb.to_mpo(AutoMPO.from_terms(sites, field), cutoff=1e-14)
    S = ma.sum_many([H1, H2], maxdim=20)
    assert calls == [2, 2]
    target = _op(H1, sites) + _op(H2, sites)
    assert _rel(_op(S, sites), target) < 1e-12


def test_densitymatrix_route_runs_on_jax():
    """The same sums under backend.set_backend("jax"), built after the
    switch, against the dense sum and against the NumPy run from the same
    random inputs."""
    pytest.importorskip("jax", reason="the JAX backend needs jax installed")
    results = {}
    try:
        for name in ("numpy", "jax"):
            bk.set_backend(name)
            sites, states = _inputs(8, (4, 6, 5))
            ref = _dense_sum(states, sites)
            full = ma.sum_many(states, method="densitymatrix")
            cut = ma.sum_many(states, cutoff=0.0, maxdim=4,
                              method="densitymatrix")
            _assert_center_last(cut)
            results[name] = (ref, _vec(full, sites), _vec(cut, sites))
    finally:
        bk.set_backend("numpy")
        bk.set_pad_bonds(None)
        bk.set_jit("auto")
    ref, full, cut = results["jax"]
    assert _rel(full, ref) < 1e-10
    assert _rel(cut, results["numpy"][2]) < 1e-8
