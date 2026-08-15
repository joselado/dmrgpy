"""Coverage for dmrgpy.blasthreads -- the opt-in BLAS/LAPACK thread limiter.

See the module's own docstring for the measurements motivating it. These
tests deliberately assert only *behaviour* (limits applied inside, restored
after, no-op degradation when threadpoolctl is missing) and never a timing,
which would be flaky on any shared machine.
"""
import warnings

import pytest

from dmrgpy import blasthreads


def _have_threadpoolctl():
    try:
        import threadpoolctl  # noqa: F401
        return True
    except ImportError:
        return False


def test_limit_is_a_context_manager_that_restores():
    """Inside the block every threadpool reports the requested limit; after
    it, the original counts are back. Skipped rather than xfailed without
    threadpoolctl -- there is nothing to restore in that case."""
    if not _have_threadpoolctl():
        pytest.skip("threadpoolctl not installed")
    before = {(d["internal_api"], d.get("filepath")): d["num_threads"]
              for d in blasthreads.info()}
    if not before:
        pytest.skip("no threadpool-backed library loaded")
    with blasthreads.limit(1):
        inside = [d["num_threads"] for d in blasthreads.info()]
        assert inside and all(t == 1 for t in inside)
    after = {(d["internal_api"], d.get("filepath")): d["num_threads"]
             for d in blasthreads.info()}
    assert after == before


def test_limit_restores_even_when_the_block_raises():
    if not _have_threadpoolctl():
        pytest.skip("threadpoolctl not installed")
    before = [d["num_threads"] for d in blasthreads.info()]
    if not before:
        pytest.skip("no threadpool-backed library loaded")
    with pytest.raises(ValueError):
        with blasthreads.limit(1):
            raise ValueError("boom")
    assert [d["num_threads"] for d in blasthreads.info()] == before


def test_limit_is_a_no_op_without_threadpoolctl(monkeypatch):
    """Without threadpoolctl the context manager must still run the block,
    warn once, and not raise -- it is a performance aid, never a hard
    requirement."""
    monkeypatch.setattr(blasthreads, "_controller", lambda: None)
    monkeypatch.setattr(blasthreads, "_WARNED", [False])
    ran = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with blasthreads.limit(1):
            ran.append(True)
        with blasthreads.limit(1):   # second time: no further warning
            ran.append(True)
    assert ran == [True, True]
    assert sum(1 for w in caught if issubclass(w.category, RuntimeWarning)) == 1


def test_hints_are_strings():
    assert "MKL_NUM_THREADS=1" in blasthreads.env_hint()
    assert isinstance(blasthreads.current_hint(), str)
    assert blasthreads.current_hint()


def test_gs_energy_is_unchanged_under_a_thread_limit():
    """The limiter must not change any result -- same ground-state energy to
    DMRG's own convergence, inside and outside the context."""
    from dmrgpy import spinchain
    def energy():
        sc = spinchain.Spin_Chain(["1/2"] * 8, itensor_version="python")
        h = 0
        for i in range(7):
            h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] + sc.Sz[i] * sc.Sz[i + 1]
        sc.set_hamiltonian(h)
        sc.maxm = 20
        return sc.gs_energy()
    plain = energy()
    with blasthreads.limit(1):
        limited = energy()
    assert limited == pytest.approx(plain, abs=1e-8)
