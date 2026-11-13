"""The leading-log flow on the device (leadinglog_device.py) against the
NumPy flow of leadinglog.py.

The device route is selected by pyitensor.backend.set_backend("jax"). It
writes the same right-hand side in jax.numpy and takes SciPy's DOP853 steps
with the state on the device, so it must give the NumPy numbers to
roundoff: on the real route, the complex route (a field along y), at T>0
with both on-shell kernels, at a stop of the flow, in both forms of the tip
part (TIP_FORM), and with the tip part cut into chunks. Here JAX runs on the CPU; the
numbers on a GPU are in docs/new_methods/kondo_leading_log/README.md.
"""
import warnings
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import DOP853

jax = pytest.importorskip("jax")

from dmrgpy.kondospectrumtk import leadinglog, leadinglog_device
from dmrgpy.pyitensor import backend
from test_kondo_spectrum_tmatrix import _ks, _spin

KB = 8.617333262e-5
EVS = np.linspace(-3e-3, 3e-3, 9) + 1.3e-5
# the free spin at Jrho_s=-0.2 stops at D=3.7e-4, and every bias within
# 10 D of zero is NaN, so its biases reach further out
EVS_WIDE = np.linspace(-9e-3, 9e-3, 9) + 1.3e-5


@pytest.fixture
def on_jax():
    name = backend.backend_name()
    backend.set_backend("jax")
    try:
        yield
    finally:
        backend.set_backend(name)


def _chain(n, by=0., T=0.):
    """n spins 1/2, anisotropic exchange 1 meV, a field along z that grows
    along the chain, a weak one along x and by along y; tip on spin 0"""
    S = _spin(0.5)

    def at(a, i):
        return np.kron(np.kron(np.eye(2**i), S[a]), np.eye(2**(n - i - 1)))
    H = sum(1e-3*(at(0, i)@at(0, i+1) + at(1, i)@at(1, i+1) + 0.8*at(2, i)@at(2, i+1))
            for i in range(n - 1))
    H = H + sum(1e-4*(1 + 0.1*i)*at(2, i) + 3e-5*at(0, i) + by*at(1, i) for i in range(n))
    ks = _ks(H, [at(a, 0) for a in range(3)])
    if T > 0.:
        ks.T, ks.kB = T, KB
        p = np.exp(-ks.e/(KB*T))
        ks.p = p/p.sum()
    return ks


def _free_spin():
    S = _spin(0.5)
    return SimpleNamespace(e=np.zeros(2), p=np.array([.5, .5]), dim=2, T=0., kB=KB,
                           Sx=S[0], Sy=S[1], Sz=S[2])


CASES = {
    "real": (lambda: _chain(3), {}),
    "complex": (lambda: _chain(3, by=3e-5), {}),
    "thermal": (lambda: _chain(3, T=2.), {}),
    "thermal, onshell T0": (lambda: _chain(3, T=2.), {"onshell": "T0"}),
    "stopped": (_free_spin, {"Jrho_s": -0.2, "eVs": EVS_WIDE}),
}


def _ll(ks, kw):
    kw = dict(kw)
    J = kw.pop("Jrho_s", -0.05)
    eVs = kw.pop("eVs", EVS)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return leadinglog.leading_log_dIdV(ks, eVs, J, U=0.25, return_info=True, **kw)


def _same(ref, got):
    (a, ia), (b, ib) = ref, got
    assert np.array_equal(np.isnan(a), np.isnan(b))
    fin = np.isfinite(a)
    assert np.sum(fin) >= 2
    assert np.max(np.abs(a[fin] - b[fin])) < 1e-12*np.max(np.abs(a[fin]))
    assert (ia["D_stop"] is None) == (ib["D_stop"] is None)
    if ia["D_stop"] is not None:
        assert ib["D_stop"] == pytest.approx(ia["D_stop"], rel=1e-12)
    assert ib["max_coupling"] == pytest.approx(ia["max_coupling"], rel=1e-12)
    assert np.allclose(ib["coupling"], ia["coupling"], rtol=1e-12, atol=0.)


@pytest.mark.parametrize("form", ["columns", "block"])
@pytest.mark.parametrize("name", list(CASES))
def test_the_device_flow_is_the_numpy_flow(name, form, on_jax, monkeypatch):
    """Both forms of the tip part. Measured on 2026-11-12: below 1.1e-15
    relative on every case, and the stop at the same cutoff to 2e-15."""
    monkeypatch.setattr(leadinglog_device, "TIP_FORM", form)
    make, kw = CASES[name]
    ks = make()
    backend.set_backend("numpy")
    ref = _ll(ks, kw)
    backend.set_backend("jax")
    _same(ref, _ll(ks, kw))


@pytest.mark.parametrize("form", ["columns", "block"])
def test_chunks_give_one_chunk(form, on_jax, monkeypatch):
    """At T>0 there are 8 occupied columns; a tiny CHUNK_BYTES cuts the
    biases one by one in both sums of the tip vertex, and the columns one
    by one in the column form."""
    monkeypatch.setattr(leadinglog_device, "TIP_FORM", form)
    ks = _chain(3, T=2.)
    ref = _ll(ks, {})
    monkeypatch.setattr(leadinglog_device, "CHUNK_BYTES", 2e3)
    _same(ref, _ll(ks, {}))


def test_the_integrator_takes_the_steps_of_scipy(on_jax):
    """_Dop853 against scipy.integrate.DOP853 on a damped rotation. The
    error estimate is a difference of nearly equal stages, so its last
    digits depend on the order of the sum, and the next step size moves by
    about 1e-7 of itself (measured 2.3e-7 after the first step). So the
    test asks for the same initial step, the same number of steps, step
    times within 1e-5 and the same state at each of them within the
    tolerance; and for a
    dense output that ends on the two states of the last step."""
    jnp = backend.xp()
    M = np.array([[-0.3, 2.0, 0.1], [-2.0, -0.1, 0.4], [0.05, -0.4, -1.0]])

    def f(t, y):
        return M @ y*np.cos(t)
    y0 = np.array([1., 0.5, -0.2])
    ref = DOP853(f, 0., y0, 7., rtol=1e-10, atol=1e-14)
    dev = leadinglog_device._Dop853(lambda t, y: jnp.asarray(M) @ y*np.cos(t), 0.,
                                    jnp.asarray(y0), 7., 1e-10, 1e-14)
    assert dev.h_abs == pytest.approx(ref.h_abs, rel=1e-14)
    while ref.status == "running":
        ref.step()
        dev.step()
        assert dev.t == pytest.approx(ref.t, rel=1e-5)
        # the state at the device's own time, from SciPy's interpolant
        assert np.allclose(np.asarray(dev.y), ref.dense_output()(dev.t), rtol=1e-9,
                           atol=1e-12)
    assert dev.status == "finished" and dev.t == ref.t
    at = dev.dense_output()
    assert np.allclose(np.asarray(at(dev.t_old)), np.asarray(dev.y_old), rtol=1e-14, atol=1e-16)
    assert np.allclose(np.asarray(at(dev.t)), np.asarray(dev.y), rtol=1e-14, atol=1e-16)
    s = 0.37*ref.t_old + 0.63*ref.t
    assert np.allclose(np.asarray(at(s)), ref.dense_output()(s), rtol=1e-9, atol=1e-12)


def test_mem_gb_chunks_the_biases_by_the_device_model(on_jax):
    """On the device mem_gb reads leadinglog_device.memory_model: a budget
    for 8 biases (16 energies) and a little more cuts 24 biases in 3
    chunks, and the result is that of one flow."""
    ks = _chain(3, T=2.)
    eVs = np.linspace(-3e-3, 3e-3, 24) + 1.3e-5
    fixed, per_w = leadinglog_device.memory_model(ks.dim, 8, real=True)
    one = _ll(ks, {"eVs": eVs})
    cut = _ll(ks, {"eVs": eVs, "mem_gb": (fixed + 17*per_w)/1e9})
    assert one[1]["n_chunks"] == 1 and cut[1]["n_chunks"] == 3
    assert np.max(np.abs(cut[0]/one[0] - 1)) < 1e-10


def test_the_allocator_limit_caps_the_budget(on_jax, monkeypatch):
    """On a GPU the budget is at most allocator_budget(), the share of JAX's
    allocator limit: a flow above it is cut in chunks with mem_gb=None, a
    larger mem_gb warns and is lowered to it, and when not even 8 biases fit
    under it the error names XLA_PYTHON_CLIENT_MEM_FRACTION. JAX's CPU
    device reports no limit, so nothing changes there."""
    assert leadinglog_device.allocator_budget() is None
    ks = _chain(3, T=2.)
    eVs = np.linspace(-3e-3, 3e-3, 24) + 1.3e-5
    fixed, per_w = leadinglog_device.memory_model(ks.dim, 8, real=True)
    one = _ll(ks, {"eVs": eVs})
    monkeypatch.setattr(leadinglog_device, "allocator_budget", lambda: 1e15)
    assert _ll(ks, {"eVs": eVs})[1]["n_chunks"] == 1
    monkeypatch.setattr(leadinglog_device, "allocator_budget", lambda: fixed + 17*per_w)
    cut = _ll(ks, {"eVs": eVs})
    assert cut[1]["n_chunks"] == 3
    assert np.max(np.abs(cut[0]/one[0] - 1)) < 1e-10
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        big = leadinglog.leading_log_dIdV(ks, eVs, -0.05, U=0.25, return_info=True,
                                          mem_gb=1e3)
    assert big[1]["n_chunks"] == 3
    assert any("XLA_PYTHON_CLIENT_MEM_FRACTION" in str(w.message) for w in wl)
    monkeypatch.setattr(leadinglog_device, "allocator_budget", lambda: fixed + 15*per_w)
    with pytest.raises(ValueError, match="XLA_PYTHON_CLIENT_MEM_FRACTION"):
        _ll(ks, {"eVs": eVs})
