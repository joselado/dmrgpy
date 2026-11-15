"""Elementwise (Hadamard) products of MPS by recursive sketched
interpolation (arXiv:2602.17974, pyitensor/rsi.py), the public
hadamard_product/elementwise_map/get_participation_entropy built on it, on
itensor_version=3, "python" and "julia_live"."""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.pyitensor import rsi
from dmrgpy import cppext
from _helpers import julia_live_param

BACKENDS = (["python"] + ([3] if cppext.available(3) else [])
            + [julia_live_param()])


def _dense(sc, wf):
    """The dense amplitudes of an MPS of any backend, from its site arrays"""
    if hasattr(wf, "jlmps"):
        from dmrgpy.mpsjulialive.juliasession import Main
        return rsi.dense_vector([np.asarray(a) for a in
                                 Main.mps_site_arrays(wf.jlmps)])
    return rsi.dense_vector(sc._session.mps_site_arrays(wf.cpp_handle))


def _random_tt(rng, n, d, chi):
    out, dl = [], 1
    for j in range(n):
        dr = 1 if j == n-1 else min(chi, d**(j+1), d**(n-j-1))
        out.append(rng.standard_normal((dl, d, dr))
                   + 1j*rng.standard_normal((dl, d, dr)))
        dl = dr
    return out


@pytest.mark.parametrize("seed", range(4))
def test_kernel_is_exact_at_full_bond_dimension(seed):
    """At maxdim = chi1*chi2 the sketch spans the whole product, so RSI is
    the exact product, whatever the random sketch."""
    rng = np.random.default_rng(seed)
    a, b = _random_tt(rng, 9, 2, 4), _random_tt(rng, 9, 2, 4)
    g = rsi.dense_vector(a)*rsi.dense_vector(b)
    c = rsi.hadamard_arrays([a, b], 16, oversampling=0, seed=seed)
    assert np.linalg.norm(rsi.dense_vector(c)-g) <= 1e-12*np.linalg.norm(g)
    k = rsi.kronecker_arrays([a, b])
    assert np.linalg.norm(rsi.dense_vector(k)-g) <= 1e-12*np.linalg.norm(g)


def test_kernel_three_factors_and_a_map():
    rng = np.random.default_rng(7)
    a = _random_tt(rng, 8, 3, 3)
    va = rsi.dense_vector(a)
    c = rsi.hadamard_arrays([a, a, a], 27, seed=1)
    assert np.linalg.norm(rsi.dense_vector(c)-va**3) <= 1e-11*np.linalg.norm(va**3)
    m = rsi.hadamard_arrays([a], 9, seed=1, func=lambda x: x**2)
    assert np.linalg.norm(rsi.dense_vector(m)-va**2) <= 1e-11*np.linalg.norm(va**2)


def test_kernel_does_not_underflow_on_a_long_chain():
    """Amplitudes of a 1500-site state are 2^-750, and their square
    underflows unless the scale is carried apart."""
    n = 1500
    a = [np.full((1, 2, 1), 1/np.sqrt(2)) for _ in range(n)]
    c = rsi.hadamard_arrays([a, a], 4)
    logs = sum(np.log(abs(x.reshape(-1)[0])) for x in c)
    assert all(np.isfinite(x).all() for x in c)
    assert logs == pytest.approx(n*np.log(0.5), rel=1e-10)


def _heisenberg(n, backend, maxm=40):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=backend,
                              maxm=maxm, nsweeps=12)
    h = sum(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h+0.1*sc.Sz[0])
    return sc


@pytest.mark.parametrize("backend", BACKENDS)
def test_probability_of_a_ground_state(backend):
    """|psi(s)|^2 of a DMRG ground state as an MPS: it converges to the
    dense product with the bond dimension and sums to one."""
    sc = _heisenberg(10, backend)
    wf = sc.get_gs()
    psi = _dense(sc, wf)
    pex = np.abs(psi)**2
    errs = []
    for maxm in (6, 12, 64):
        p = sc.hadamard_product(wf, wf, conj=(True, False), maxm=maxm)
        pv = _dense(sc, p)
        errs.append(np.linalg.norm(pv-pex)/np.linalg.norm(pex))
    # it converges with the bond dimension, to rounding once maxm holds p
    assert errs[0] > 10*errs[1]
    assert errs[1] > 1e3*errs[2]
    assert errs[2] < 1e-10
    assert pv.sum().real == pytest.approx(1.0, abs=1e-10)
    pd = sc.hadamard_product(wf, wf, conj=(True, False), method="direct",
                             maxm=64)
    pdv = _dense(sc, pd)
    assert np.linalg.norm(pdv-pex) < 1e-4*np.linalg.norm(pex)


@pytest.mark.parametrize("backend", BACKENDS)
def test_participation_entropy_against_ED(backend):
    sc = _heisenberg(10, backend, maxm=40)
    for q in (2, 3):
        sd = sc.get_participation_entropy(q=q, maxm=64)
        se = sc.get_participation_entropy(q=q, mode="ED")
        assert sd == pytest.approx(se, abs=1e-5)


@pytest.mark.parametrize("backend", BACKENDS)
def test_elementwise_map_and_ED_mode(backend):
    sc = _heisenberg(8, backend)
    wf = sc.get_gs()
    psi = _dense(sc, wf)
    g = sc.elementwise_map(wf, lambda x: x**2, maxm=64)
    gv = _dense(sc, g)
    assert np.linalg.norm(gv-psi**2) < 1e-10*np.linalg.norm(psi**2)
    we = sc.get_gs(mode="ED")
    pe = sc.hadamard_product(we, we, conj=(True, False))
    assert np.allclose(pe.v, np.abs(we.v)**2)


@pytest.mark.parametrize("backend", BACKENDS)
def test_arguments_are_checked(backend):
    sc = _heisenberg(6, backend)
    wf = sc.get_gs()
    with pytest.raises(ValueError):
        sc.hadamard_product(wf)
    with pytest.raises(ValueError):
        sc.hadamard_product(wf, wf, conj=(True,))
    with pytest.raises(ValueError):
        sc.hadamard_product(wf, wf, method="tci")
    with pytest.raises(ValueError):
        sc.get_participation_entropy(wf, q=1.5)
    with pytest.raises(ValueError):
        sc.hadamard_product(wf, sc.get_gs(mode="ED"))
