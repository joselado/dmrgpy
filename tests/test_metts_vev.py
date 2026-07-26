"""Regression coverage for Many_Body_Chain.metts_vev() (finite-temperature
expectation values via the METTS -- Minimally Entangled Typical Thermal
States -- sampling algorithm, White & Stoudenmire, New Journal of Physics
12, 055026 (2010), arXiv:1002.1305).

METTS is implemented for itensor_version="python" (pyitensor/metts.py) and
itensor_version=3 (mpscpp3/chain_session.h's Chain::metts_vev, a direct
port of the same algorithm onto real ITensor v3 -- see its own comment)::
imaginary-time TDVP evolution of a classical product state by
e^{-beta H/2}, then a sequential/"perfect sampling" collapse back down to
a new product state, alternating the collapse basis between Sz and Sx for
ergodicity. This is a genuinely different, independent finite-temperature
method from the two already covered elsewhere (thermal.py's ancilla-
purification anneal(), tested in test_thermal_purification_VS_exact's
underlying example, and vevtk/thermalvev.py's exact ED Boltzmann sum over
excited states, tested in test_thermal_vev.py) -- so this is validated
directly against the same kind of full-spectrum ED reference those use,
not against either of them.

METTS is a Monte Carlo method: exact agreement isn't expected, only
agreement within a generous multiple of the reported (Markov-correlated,
so likely optimistic) standard error -- these seeds/sample counts were
picked empirically (see the module-level comment on parameters below) to
give several-sigma headroom without an unreasonably slow test.
"""
import pytest

from dmrgpy import cppext, spinchain
from dmrgpy.pyitensor.index import reseed_id_counter_past
from dmrgpy.pyitensor.metts import _run_chain_worker

from _helpers import setup_backend

VERSIONS = [
    "python",
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
]


def _heisenberg_field_chain(n, version, J=1.0, B=0.3):
    sc = spinchain.Spin_Chain(["S=1/2" for _ in range(n)])
    setup_backend(sc, version)
    h = 0
    for i in range(n - 1):
        h = h + J * sc.Sx[i] * sc.Sx[i + 1] + J * sc.Sy[i] * sc.Sy[i + 1] + J * sc.Sz[i] * sc.Sz[i + 1]
    for i in range(n):
        h = h + B * sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc, h


def test_metts_matches_ed_two_site_magnetization():
    """2-site Heisenberg chain + field: METTS's own reported standard
    error at these sample counts is ~0.013 (measured empirically); the
    actual deviation from the exact ED thermal average is ~0.002, so
    tol=0.05 (~4x that error) is generous headroom against an unlucky
    seed while still catching a real algorithmic regression (e.g. a sign
    error in the collapse or the imaginary-time step, which previously
    produced ~10x larger discrepancies during development). python-only:
    ITensor v3's two-site TDVP aborts the whole process for chains below
    3 sites (same vendored-ITensor limitation as its two-site dmrg(), see
    CLAUDE.md and test_metts_vev_v3_rejects_short_chain below)."""
    sc, h = _heisenberg_field_chain(2, "python")
    T = 1.0

    ed_sz = sc.vev(sc.Sz[0], mode="ED", T=T)
    metts_sz, err_sz = sc.metts_vev(sc.Sz[0], T, nsamples=600, nwarmup=50,
                                     dbeta_half_step=0.05, seed=2024)

    assert metts_sz.real == pytest.approx(ed_sz.real, abs=0.05)
    assert err_sz < 0.03  # sanity check the stochastic error itself is small


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_matches_ed_three_site_energy_and_magnetization(version):
    """3-site chain: cross-checks both the thermal energy (a stringent
    check -- a bug in the imaginary-time TDVP propagator would show up as
    a systematic energy bias, not just sampling noise) and a local <Sz>,
    each against the exact ED Boltzmann average. Parametrized over the
    pure-Python and ITensor v3 backends -- both implement the identical
    algorithm (see module docstring), so both must agree with ED."""
    sc, h = _heisenberg_field_chain(3, version, B=0.4)
    T = 0.8

    ed_e = sc.vev(h, mode="ED", T=T)
    ed_sz = sc.vev(sc.Sz[1], mode="ED", T=T)

    metts_e, err_e = sc.metts_vev(h, T, nsamples=300, nwarmup=30,
                                   dbeta_half_step=0.08, seed=77)
    metts_sz, err_sz = sc.metts_vev(sc.Sz[1], T, nsamples=300, nwarmup=30,
                                     dbeta_half_step=0.08, seed=77)

    assert metts_e.real == pytest.approx(ed_e.real, abs=0.05)
    assert metts_sz.real == pytest.approx(ed_sz.real, abs=0.05)


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_vev_multi_op_matches_separate_calls(version):
    """Op may be a list/tuple of MultiOperators, measured together on one
    shared METTS sample chain (see vevtk/mettsvev.py's docstring) instead
    of forcing one full sampling run per operator. Given the same seed,
    the shared-chain multi-op path must reproduce -- bit for bit, not just
    approximately -- the same per-operator results as calling metts_vev
    separately for each operator, since it's sampling the exact same
    Markov chain of METTS states and just measuring more operators against
    each sampled state."""
    sc, h = _heisenberg_field_chain(3, version, B=0.4)
    T = 0.8
    kwargs = dict(nsamples=40, nwarmup=10, dbeta_half_step=0.08, seed=77)

    metts_e, err_e = sc.metts_vev(h, T, **kwargs)
    metts_sz, err_sz = sc.metts_vev(sc.Sz[1], T, **kwargs)

    results = sc.metts_vev([h, sc.Sz[1]], T, **kwargs)
    assert len(results) == 2
    assert results[0] == (metts_e, err_e)
    assert results[1] == (metts_sz, err_sz)


def test_metts_vev_requires_supported_backend():
    """METTS needs imaginary-time TDVP + sequential-sampling MPS collapse;
    mpscpp2 (ITensor v2) has no TDVP module at all -- metts_vev must raise
    a clear error rather than silently dispatching to the wrong backend or
    crashing deep inside cppext."""
    sc, h = _heisenberg_field_chain(2, "python")
    sc.itensor_version = 2  # simulate a non-python/non-v3 backend without touching cppext
    with pytest.raises(NotImplementedError):
        sc.metts_vev(sc.Sz[0], 1.0, nsamples=2, nwarmup=1)


def test_metts_vev_v3_rejects_short_chain():
    """ITensor v3's two-site TDVP hits the same 'LocalOp is default
    constructed' SIGABRT as its two-site dmrg() for chains below 3 sites
    (see CLAUDE.md's documented mpscpp3 bug) -- confirmed directly during
    development (a 2-site metts_vev call crashed the whole process before
    this guard was added). metts_vev must raise a catchable Python
    exception instead, same as the itensor_version dispatch guard above.
    This only exercises the Python-side guard (fires before ever touching
    self._session), so it needs no compiled extension to run."""
    sc, h = _heisenberg_field_chain(2, "python")
    sc.itensor_version = 3
    with pytest.raises(NotImplementedError):
        sc.metts_vev(sc.Sz[0], 1.0, nsamples=2, nwarmup=1)


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_vev_rejects_non_positive_temperature(version):
    """T<=0 must raise rather than silently flip the imaginary-time step's
    sign (beta=1/T negative would make the 'decay' amplify high-energy
    components instead of suppressing them, a silently wrong thermal
    average with no error). Checked in vevtk/mettsvev.py *before*
    dispatching to either backend's own session: mpscpp3/chain_session.h's
    Chain::metts_vev does validate this too, but via ITensor's Error(),
    which calls abort() directly (see itensor/util/error.h) rather than
    throwing a catchable exception -- confirmed directly, a T<=0 call used
    to SIGABRT the whole interpreter instead of raising. That C++-side
    check is only a last-resort safety net against direct
    _session.metts_vev() misuse, not something a Python test can exercise
    without crashing the test process, so this test only ever reaches the
    Python-level guard."""
    sc, h = _heisenberg_field_chain(3, version)
    with pytest.raises(ValueError):
        sc.metts_vev(sc.Sz[0], -1.0, nsamples=2, nwarmup=1)
    with pytest.raises(ValueError):
        sc.metts_vev(sc.Sz[0], 0.0, nsamples=2, nwarmup=1)


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_vev_rejects_empty_basis_ops(version):
    """An empty basis_ops would otherwise index basis_ops[0] out of bounds
    (undefined behavior in the C++ port) or KeyError deep inside the
    Python one -- caught by vevtk/mettsvev.py's own check before either
    backend's session ever sees it (see
    test_metts_vev_rejects_non_positive_temperature's docstring on why
    the C++-side check alone isn't something a Python test can rely on)."""
    sc, h = _heisenberg_field_chain(3, version)
    with pytest.raises(ValueError):
        sc.metts_vev(sc.Sz[0], 1.0, nsamples=2, nwarmup=1, basis_ops=())


def test_metts_vev_njobs_parallel_chains_match_ed():
    """njobs>1 (pyitensor/metts.py's parallel-independent-Markov-chains
    path, python-only -- see that module's own docstring) must still
    reproduce the same physics as the sequential njobs=1 path: nsamples
    split across njobs independent chains, each with its own nwarmup and
    an independently-seeded RNG, then combined via
    _pool_chain_results()'s exact pooled-mean/pooled-variance formula.
    Same generous tolerance as the other ED cross-checks above -- this
    isn't re-testing METTS's own accuracy, only that splitting the work
    across processes doesn't silently change the sampled distribution."""
    sc, h = _heisenberg_field_chain(3, "python", B=0.4)
    T = 0.8

    ed_e = sc.vev(h, mode="ED", T=T)
    ed_sz = sc.vev(sc.Sz[1], mode="ED", T=T)

    metts_e, err_e = sc.metts_vev(h, T, nsamples=300, nwarmup=30,
                                   dbeta_half_step=0.08, seed=77, njobs=3)
    metts_sz, err_sz = sc.metts_vev(sc.Sz[1], T, nsamples=300, nwarmup=30,
                                     dbeta_half_step=0.08, seed=77, njobs=3)

    assert metts_e.real == pytest.approx(ed_e.real, abs=0.05)
    assert metts_sz.real == pytest.approx(ed_sz.real, abs=0.05)


def test_metts_vev_njobs_worker_survives_colliding_index_ids():
    """A spawned njobs>1 worker is a fresh interpreter: pyitensor/
    index.py's process-local Index _id_counter restarts at 0 there, but
    the H/sites/ops MPO/MPS objects it receives were pickled in the
    *parent* process and carry whatever ids that process's own counter
    had already assigned (pickling round-trips Index._id verbatim, it
    never goes through Index.__init__/_id_counter). Since
    Index.__eq__/__hash__ compare id (and prime level) alone, with no
    dim/tags check, a worker that mints its own new indices starting from
    0 (e.g. every sample's classical-product-state bond links) can
    collide with an unrelated, differently-shaped Index from the
    unpickled task -- confirmed directly, this made
    metts_vev(..., njobs=2) reliably raise 'ValueError: shape-mismatch
    for sum' before pyitensor/metts.py's _run_chain_worker() started
    calling index.py's reseed_id_counter_past() first. Rather than rely
    on multiprocessing.Pool actually reproducing that race (order/timing
    dependent, and confirmed to intermittently pass by luck when *other*
    tests already pushed this process's own id_counter to large values
    before this one runs), this reproduces the exact failure mode
    directly and deterministically: force this process's counter back
    down into a colliding low range (standing in for a freshly spawned
    worker's real fresh-interpreter id_counter=0), then call
    _run_chain_worker() in-process (no actual multiprocessing.Pool
    needed -- it's a plain function) and confirm it still produces a
    normal, non-crashing result."""
    sc, h = _heisenberg_field_chain(3, "python", B=0.4)
    T = 0.8
    # Populate self._session.H/sites the same way metts_vev() itself does,
    # via one cheap real call, then reach into the session for the same
    # (H, sites, ops) triple _run_chain_worker() is handed under njobs>1.
    sc.metts_vev(sc.Sz[0], T, nsamples=1, nwarmup=0, dbeta_half_step=0.5, seed=1)
    chain = sc._session
    from dmrgpy.pyitensor.autompo import AutoMPO
    from dmrgpy.pyitensor.mpobuilder import to_mpo
    op = to_mpo(AutoMPO.from_terms(chain.sites, sc.Sz[1].to_terms()),
                cutoff=1e-14, maxdim=chain.mpomaxm)

    reseed_id_counter_past(-1)  # force this process's next Index() id to 0
    means, stderrs, n = _run_chain_worker(
        (chain.H, chain.sites, [op], 1.0 / T, 5, 2, 0.2, chain.cutoff,
         chain.maxm, ("Sz", "Sx"), None, 123, 30))
    assert n == 5
    assert len(means) == 1 and len(stderrs) == 1


def test_metts_vev_njobs_rejects_v3():
    """itensor_version=3's Chain is a single live in-process C++ session
    (see CLAUDE.md's cpp_handle section) with no per-worker copy a
    multiprocessing pool could hand out, unlike pyitensor's picklable
    plain-Python/numpy H/sites -- njobs>1 must raise rather than silently
    ignoring the request or crashing deep inside a forked/pickled C++
    handle."""
    sc, h = _heisenberg_field_chain(3, 3)
    with pytest.raises(NotImplementedError):
        sc.metts_vev(sc.Sz[0], 1.0, nsamples=4, nwarmup=1, njobs=2)


def test_metts_vev_rejects_non_positive_njobs():
    """njobs<1 would otherwise silently fall through the njobs==1 fast
    path's `if` check as false and hit divmod(nsamples, 0) -- caught by
    metts_thermal_average's own guard before that, same convention as
    every other parameter here."""
    sc, h = _heisenberg_field_chain(3, "python")
    with pytest.raises(ValueError):
        sc.metts_vev(sc.Sz[0], 1.0, nsamples=4, nwarmup=1, njobs=0)


@pytest.mark.parametrize("version", VERSIONS)
def test_metts_vev_rejects_non_positive_nsamples(version):
    """nsamples<1 would otherwise leave the retained-samples list empty
    and divide by zero when averaging (silently producing NaN in the C++
    port before this guard was added) -- caught by vevtk/mettsvev.py's own
    check before either backend's session ever sees it (see
    test_metts_vev_rejects_non_positive_temperature's docstring on why the
    C++-side check alone isn't something a Python test can rely on)."""
    sc, h = _heisenberg_field_chain(3, version)
    with pytest.raises(ValueError):
        sc.metts_vev(sc.Sz[0], 1.0, nsamples=0, nwarmup=1)
