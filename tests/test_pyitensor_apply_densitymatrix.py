"""pyitensor's density-matrix MPO application, applyMPO(method=
"densitymatrix") (mpsalgebra._apply_densitymatrix), against the exact
route (_apply_chain) and against the product contracted out to a dense
vector.

The density-matrix algorithm (Stoudenmire and White, New J. Phys. 12,
055026 (2010)) keeps at every bond the leading eigenvectors of the reduced
density matrix of the partially projected product. These are the right
Schmidt vectors that the exact route's truncating SVD sweep keeps at the
same bond, so the two methods are the same truncation and every check
below holds for both: the same state, an error bounded by the Schmidt
tails of the exact product, and a result that is an orthogonal projection
of it.

The chains have 10 sites or fewer, so the exact product fits in a dense
vector. That vector is the reference, and it shares no code with either
method.
"""
import numpy as np
import pytest

from dmrgpy.pyitensor import backend as bk
from dmrgpy.pyitensor import mpobuilder as mb
from dmrgpy.pyitensor import mpsalgebra as ma
from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyitensor.tensor import commonIndex, dag, prime

SPIN_HALF = 2


def heisenberg_terms(n):
    """Nearest-neighbour Heisenberg chain, MPO bond dimension 5, plus a
    field on site 1 so that a one-site chain is not the zero operator."""
    terms = [(1.0, [(op, i), (op, i + 1)])
             for i in range(1, n) for op in ("Sx", "Sy", "Sz")]
    return terms + [(0.3, [("Sx", 1)])]


def long_range_terms(n):
    """All-to-all couplings decaying as 1/r^2, with complex hoppings S+_i
    S-_j that have no Hermitian partner: the operator is complex and not
    Hermitian, so a bra copy of K with its legs swapped, or a missing
    conjugation, cannot cancel out. Its MPO bond dimension grows with n."""
    terms = []
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            r = j - i
            terms.append((1.0 / r**2, [("Sz", i), ("Sz", j)]))
            terms.append((0.5 * np.exp(0.4j * r) / r**2, [("S+", i), ("S-", j)]))
    return terms + [(0.5 + 0.25j, [("Sx", 2), ("Sy", n)])]


MODELS = {"heisenberg": heisenberg_terms, "long_range": long_range_terms}


def _case(model, n, chi, seed=7):
    sites = SiteX([SPIN_HALF] * n)
    ampo = AutoMPO.from_terms(sites, MODELS[model](n))
    K = mb.to_mpo(ampo, cutoff=1e-14)
    x = ma.randomMPS(sites, chi, rng=np.random.default_rng(seed))
    return sites, ampo, K, x


def _vec(mps, sites):
    """An MPS as a dense vector, site 1 most significant. applyMPO leaves
    the physical legs primed, so they are unprimed first."""
    m = mps.copy().noPrime("Site")
    T = m.A(1)
    for i in range(2, m.length() + 1):
        T = T * m.A(i)
    arr = T.transpose_to([sites.si(i) for i in range(1, m.length() + 1)])
    return np.asarray(arr).reshape(-1)


def _exact_product(K, x, sites):
    """K|x> contracted site by site into a dense vector, with no
    factorization at all."""
    n = K.length()
    T = K.A(1) * x.A(1)
    for i in range(2, n + 1):
        T = (T * x.A(i)) * K.A(i)
    arr = T.transpose_to([sites.si(i).prime(1) for i in range(1, n + 1)])
    return np.asarray(arr).reshape(-1)


def _rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def _schmidt_tails(vec, n, m):
    """sum_{k>m} s_k^2 at every bond of a dense n-site spin-1/2 vector,
    the weight a bond dimension m must lose there."""
    tails = []
    for b in range(1, n):
        s = np.linalg.svd(vec.reshape(2**b, 2**(n - b)), compute_uv=False)
        tails.append(float(np.sum(s[m:]**2)))
    return np.array(tails)


def _bond_dims(mps):
    return [commonIndex(mps.A(i), mps.A(i + 1)).dim
            for i in range(1, mps.length())]


def _assert_center_one(out):
    """center == 1 means sites 2..n are right-orthogonal, so the whole
    norm sits on site 1."""
    assert out.center == 1
    for j in range(2, out.length() + 1):
        T = out.A(j)
        left = commonIndex(out.A(j - 1), T)
        G = T * dag(prime(T, left))
        g = np.asarray(G.transpose_to([left, left.prime(1)]))
        assert np.max(np.abs(g - np.eye(left.dim))) < 1e-10


@pytest.mark.parametrize("model,n,chi", [("heisenberg", 10, 16),
                                         ("long_range", 8, 12),
                                         ("heisenberg", 2, 2),
                                         ("heisenberg", 1, 1)],
                         ids=["heisenberg-10", "long_range-8",
                              "heisenberg-2", "one-site"])
def test_same_state_as_the_exact_product_without_truncation(model, n, chi):
    """With cutoff 0 and a maxdim that holds any n-site state, both methods
    return K|x> itself. The random MPS has no Schmidt values near the
    density-matrix route's precision floor (about 1e-8 of the largest), so
    1e-10 is a fair tolerance for it."""
    sites, ampo, K, x = _case(model, n, chi)
    ref = _exact_product(K, x, sites)
    # the reference, against the operator's own Kronecker product
    assert _rel(ampo.dense_matrix() @ _vec(x, sites), ref) < 1e-10
    full = 2 ** (n // 2)  # the largest Schmidt rank of n spin-1/2 sites
    dm = ma.applyMPO(K, x, cutoff=0.0, maxdim=full, method="densitymatrix")
    ex = ma.applyMPO(K, x, cutoff=0.0, maxdim=full, method="exact")
    assert _rel(_vec(dm, sites), ref) < 1e-10
    assert _rel(_vec(ex, sites), ref) < 1e-10
    assert _rel(_vec(dm, sites), _vec(ex, sites)) < 1e-10
    _assert_center_one(dm)


TRUNCATIONS = {"maxdim": dict(cutoff=0.0, maxdim=8),
               "cutoff": dict(cutoff=1e-3, maxdim=None)}


@pytest.mark.parametrize("trunc", sorted(TRUNCATIONS))
@pytest.mark.parametrize("model,n,chi", [("heisenberg", 10, 16),
                                         ("long_range", 10, 12)],
                         ids=["heisenberg", "long_range"])
def test_truncation_error_is_bounded_by_the_discarded_weight(model, n, chi, trunc):
    """Both methods discard, bond by bond, the smallest Schmidt weights of a
    product that is already projected on the right, so the squared error
    is at most the sum over bonds of what each bond must lose: the Schmidt
    tail of K|x> beyond maxdim, or `cutoff` of the norm squared. With
    maxdim it is also at least the largest single tail, which no state of
    that bond dimension can beat, so the check is not vacuous. And since
    both keep the same basis at every bond, they give the same state."""
    kw = TRUNCATIONS[trunc]
    sites, _ampo, K, x = _case(model, n, chi)
    ref = _exact_product(K, x, sites)
    norm2 = np.vdot(ref, ref).real
    if kw["maxdim"] is None:
        upper, lower = (n - 1) * kw["cutoff"] * norm2, 1e-12 * norm2
    else:
        tails = _schmidt_tails(ref, n, kw["maxdim"])
        upper, lower = tails.sum(), tails.max()
    vecs = {}
    for method in ("densitymatrix", "exact"):
        out = ma.applyMPO(K, x, method=method, **kw)
        v = _vec(out, sites)
        err2 = np.vdot(ref - v, ref - v).real
        assert lower * (1 - 1e-8) <= err2 <= upper * (1 + 1e-8) + 1e-14 * norm2
        if kw["maxdim"] is not None:
            assert max(_bond_dims(out)) <= kw["maxdim"]
        vecs[method] = v
    assert _rel(vecs["densitymatrix"], vecs["exact"]) < 1e-8


@pytest.mark.parametrize("method", ["densitymatrix", "exact"])
@pytest.mark.parametrize("model", ["heisenberg", "long_range"])
def test_truncated_result_is_an_orthogonal_projection(model, method):
    """A truncated result is P K|x> with P an orthogonal projector, so its
    overlap with the exact product is real and equal to its own norm
    squared. A conjugation error in the kept basis breaks this even when
    the norm happens to look right. The norm must also be the one the rest
    of the engine reads: inner(out, out), and the norm of site 1 alone."""
    sites, _ampo, K, x = _case(model, 10, 12, seed=11)
    ref = _exact_product(K, x, sites)
    out = ma.applyMPO(K, x, cutoff=0.0, maxdim=6, method=method)
    v = _vec(out, sites)
    norm2 = np.vdot(ref, ref).real
    out2 = np.vdot(v, v).real
    assert abs(np.vdot(ref, v) - out2) < 1e-10 * norm2
    assert out2 < norm2 * (1 - 1e-6)  # something was projected out
    assert ma.inner(out, out).real == pytest.approx(out2, rel=1e-10)
    site1 = np.asarray(out.A(1).array)
    assert np.vdot(site1, site1).real == pytest.approx(out2, rel=1e-10)
    _assert_center_one(out)


def test_default_is_exact_and_unknown_methods_raise():
    """The default is the density-matrix route (the author's pick of
    2026-11-15) and equals the exact one here, an unknown name is an error
    rather than a fallback, and a zero operator gives the zero state, not
    NaN."""
    assert ma.APPLY_METHOD == "densitymatrix"  # the default since 2026-11-15
    sites, _ampo, K, x = _case("heisenberg", 6, 4)
    a = _vec(ma.applyMPO(K, x, cutoff=0.0, maxdim=4), sites)
    b = _vec(ma.applyMPO(K, x, cutoff=0.0, maxdim=4, method="exact"), sites)
    assert _rel(a, b) < 1e-12
    with pytest.raises(ValueError):
        ma.applyMPO(K, x, method="zipup")
    Z = mb.to_mpo(AutoMPO.from_terms(sites, []), cutoff=1e-14)
    z = _vec(ma.applyMPO(Z, x, maxdim=4, method="densitymatrix"), sites)
    assert np.all(np.isfinite(z)) and np.max(np.abs(z)) == 0.0


def test_densitymatrix_route_runs_on_jax():
    """The same products under backend.set_backend("jax"), built after the
    switch, against the exact product and against the NumPy run from the
    same random start."""
    pytest.importorskip("jax", reason="the JAX backend needs jax installed")
    results = {}
    try:
        for name in ("numpy", "jax"):
            bk.set_backend(name)
            sites, _ampo, K, x = _case("long_range", 8, 12)
            ref = _exact_product(K, x, sites)
            full = ma.applyMPO(K, x, cutoff=0.0, maxdim=16, method="densitymatrix")
            cut = ma.applyMPO(K, x, cutoff=0.0, maxdim=6, method="densitymatrix")
            _assert_center_one(cut)
            results[name] = (ref, _vec(full, sites), _vec(cut, sites))
    finally:
        bk.set_backend("numpy")
        bk.set_pad_bonds(None)
        bk.set_jit("auto")
    ref, full, cut = results["jax"]
    assert _rel(full, ref) < 1e-10
    assert _rel(cut, results["numpy"][2]) < 1e-8
