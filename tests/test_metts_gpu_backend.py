"""METTS (pyitensor/metts.py) on the JAX array backend: same physics, and
-- the part agreement cannot check -- the collapse step's amplitudes stay
where the rest of the engine put them.

Skipped entirely when JAX is not installed, and *not* skipped when JAX has
only a CPU device: a host-only JAX run exercises the identical code path a
device does (same immutability, same dispatch, same host round trips),
which is what these tests are about. See tests/test_pyitensor_gpu_backend.py
for the same rationale at more length, and docs/pyitensor_gpu_port_plan.md
Sec. 9 item 5 for why METTS was the port's last big open item: everything
underneath it (TDVP, dmrg.py's environments, svd.py, mpsalgebra.inner) was
already resident, so the whole surface left was collapse_to_cps().

Sizes are deliberately tiny: eager JAX is 8-14x slower than NumPy at small
bond dimension and METTS runs nwarmup+nsamples full imaginary-time
evolutions, so anything bigger costs minutes for no extra coverage.
"""

import numpy as np
import pytest

from dmrgpy import spinchain

jax = pytest.importorskip("jax", reason="the JAX backend needs jax installed")

from dmrgpy.pyitensor import backend as bk       # noqa: E402
from dmrgpy.pyitensor import metts as mettstk    # noqa: E402


@pytest.fixture
def numpy_backend_restored():
    """Process-wide state, so put it back whatever the test does."""
    yield
    bk.set_backend("numpy")
    bk.set_pad_bonds(None)
    bk.set_jit("auto")


def _chain(n=4, B=0.3):
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)],
                              itensor_version="python")
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + B * sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc


def _metts_sz(seed):
    # maxdim/cutoff come from the chain (mettsvev.py forwards MB.maxm/
    # MB.cutoff through set_sweep_params), not from metts_vev's own kwargs.
    sc = _chain()
    sc.maxm = 16
    return sc.metts_vev(sc.Sz[1], T=1.0, nsamples=4, nwarmup=1,
                        dbeta_half_step=0.25, seed=seed)


def test_metts_vev_matches_numpy(numpy_backend_restored):
    """Same seed, same Markov chain, same answer on both array libraries.

    METTS is stochastic, and the two backends agree exactly only as long
    as ~1e-15 differences in the collapse probabilities never flip an
    rng.choice draw -- if one ever did, the two runs would follow
    different trajectories and this would fail by a wide margin rather
    than subtly, which is why a 1e-10 tolerance is the right shape of
    check here even though the arithmetic is deterministic."""
    mean_np, err_np = _metts_sz(seed=7)

    bk.set_backend("jax")
    mean_jax, err_jax = _metts_sz(seed=7)

    assert abs(complex(mean_jax) - complex(mean_np)) < 1e-10
    assert abs(float(err_jax) - float(err_np)) < 1e-10


def test_collapse_transfers_only_the_probabilities(numpy_backend_restored):
    """The collapse sweep's amplitude block is (d, chi) and its
    conditional-probability vector is (d,). Before the port the whole
    block came home at every site -- via np.sum/np.abs on a device array,
    silently, with no error -- and the collapsed-prefix amplitude L went
    back out again.

    Agreement cannot see that (a round trip returns the same numbers,
    only slower), so this counts the host transfers the sweep actually
    performs and asserts on their *size*: one O(d) vector per site, and
    nothing of order chi. The MPS is handed over already centered at site
    1, so position() does no SVD and every transfer counted here is the
    collapse's own."""
    bk.set_backend("jax")
    from dmrgpy.pyitensor.mpsalgebra import randomMPS
    from dmrgpy.pyitensor.sites import SiteX

    n, d, chi = 5, 2, 8
    sites = SiteX([2] * n)   # 2 = SpinHalfSite (siteset.TYPE_CODE_TO_SITE)
    psi = randomMPS(sites, chi)
    psi.position(1)
    assert max(t.array.size for t in [psi.A(i) for i in range(1, n + 1)]) > d

    # With the eigenbasis cache a real run always builds (_metts_single_chain
    # does), the sweep's only host traffic is its own; without it,
    # _site_operator_matrix reads back d columns per site to *build* that
    # cache, which is once-per-run work rather than per-sample work.
    eigcache = mettstk._build_eigcache(sites, ["Sz"])

    sizes = []
    real_to_host = bk.to_host

    def counting_to_host(a):
        sizes.append(np.asarray(a).size)
        return real_to_host(a)

    bk.to_host = counting_to_host
    try:
        cps, outcomes = mettstk.collapse_to_cps(
            psi, sites, "Sz", np.random.default_rng(3), eigcache=eigcache)
    finally:
        bk.to_host = real_to_host

    assert len(outcomes) == n
    assert len(sizes) == n            # one transfer per site, no more
    assert set(sizes) == {d}          # ...and each is the O(d) probability vector


def test_collapsed_prefix_stays_on_the_device(numpy_backend_restored):
    """The other half of the same claim, asserted on the array type: the
    running prefix amplitude L is built from the device-resident `rot`,
    so the per-site tensor it is contracted into must itself be a device
    array. A host round trip would produce a NumPy array here and still
    give identical numbers."""
    bk.set_backend("jax")
    from dmrgpy.pyitensor.mpsalgebra import randomMPS
    from dmrgpy.pyitensor.sites import SiteX
    from dmrgpy.pyitensor.tensor import ITensor

    seen = []
    real_init = ITensor.__init__

    def recording_init(self, inds, array=None):
        real_init(self, inds, array)
        seen.append((tuple(ind.tags for ind in self.inds), array))

    n = 4
    sites = SiteX([2] * n)   # 2 = SpinHalfSite (siteset.TYPE_CODE_TO_SITE)
    psi = randomMPS(sites, 8)
    psi.position(1)

    ITensor.__init__ = recording_init
    try:
        mettstk.collapse_to_cps(psi, sites, "Sz", np.random.default_rng(5))
    finally:
        ITensor.__init__ = real_init

    prefixes = [arr for tags, arr in seen
                if len(tags) == 1 and "Link" in tags[0] and arr is not None]
    assert prefixes, "the collapse built no single-Link prefix amplitude"
    assert all(isinstance(arr, jax.Array) for arr in prefixes)


def test_njobs_is_refused_on_a_device(numpy_backend_restored):
    """multiprocessing workers are started with 'spawn', so they re-import
    pyitensor with the default NumPy backend instead of inheriting this
    process's device -- njobs>1 on a device would silently run on the host.
    backend.py's whole contract is that this fails loudly instead."""
    sc = _chain(n=3)
    bk.set_backend("jax")
    with pytest.raises(ValueError, match="njobs>1 is not supported"):
        sc.metts_vev(sc.Sz[0], T=1.0, nsamples=4, nwarmup=0,
                     dbeta_half_step=0.5, seed=1, njobs=2)
