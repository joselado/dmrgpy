"""Second- and third-order sums of kondospectrumtk/conductance.py: the
second-order term runs over the occupied initial states, and the third-order
Kondo coefficient is blocked over them with a byte budget."""
import tracemalloc
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import conductance as cd
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum


def _old_second(ks, eVs, T0=1.0, U=0.0):
    """The full sum over every initial state, as it was before."""
    eVs = np.asarray(eVs, dtype=float)
    M2 = cd._spin_matrix_elements_squared(ks)
    eps = ks.e[:, None] - ks.e[None, :]
    weight = ks.p[None, :]*M2
    x_ts = eVs[:, None, None] - eps[None, :, :]
    x_st = -eVs[:, None, None] - eps[None, :, :]
    spin = (np.einsum('fi,efi->e', weight, cd._theta_raw(ks, x_ts))
            + np.einsum('fi,efi->e', weight, cd._theta_raw(ks, x_st)))
    return 2*np.pi*T0**2*(spin + 4*U**2*np.ones_like(eVs))


def _old_third(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6, Fb=None,
               p_cut=cd.P_CUT):
    """The single-pass coefficient tensor, as it was before."""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks, p_cut)
    p = ks.p[occ]
    coeff = np.imag(cd._triple_product_coefficients(ks, occ))/2.
    eps_if = eps_im = ks.e[None, :] - ks.e[occ, None]
    e_fm = ks.e[:, None] - ks.e[None, :]
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    dim, nocc = ks.dim, len(occ)

    def one(v):
        nv = len(v)
        Th = cd._theta_raw(ks, v[:, None, None] - eps_if[None, :, :])
        Fd = Fcall((v[:, None, None] - eps_im[None, :, :]).ravel()).reshape(nv, nocc, dim)
        Fx = Fcall((v[:, None, None] - e_fm[None, :, :]).ravel()).reshape(nv, dim, dim)
        return (np.einsum('i,ifm,eif,eim->e', p, coeff, Th, Fd, optimize=True)
                + np.einsum('i,ifm,eif,efm->e', p, coeff, Th, Fx, optimize=True))
    return 4*np.pi*T0**2*Jrho_s*(one(eVs) + one(-eVs))


def _chain(T, complex_h=False, n=3):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + 3e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                      + 0.7*sc.Sz[i]*sc.Sz[i+1])
    h = h + 0.4e-3*sc.Sz[0] + 0.3e-3*sc.Sx[1]
    if complex_h: h = h + 0.5e-3*sc.Sy[0]
    sc.set_hamiltonian(h)
    return KondoSpectrum(sc, 0, T)


EV = np.linspace(-6e-3, 6e-3, 23)


@pytest.mark.parametrize("T", [0.0, 3.0])
def test_second_order_p_cut_zero_is_the_full_sum(T):
    ks = _chain(T)
    ref = _old_second(ks, EV, U=0.2)
    got = cd.second_order_dIdV(ks, EV, U=0.2, p_cut=0)
    assert np.max(np.abs(got - ref)) < 1e-13*np.max(ref)


def test_second_order_chunks_do_not_matter(monkeypatch):
    ks = _chain(3.0)
    ref = cd.second_order_dIdV(ks, EV, U=0.2)
    monkeypatch.setattr(cd, "CHUNK_ENTRIES", 1)
    got = cd.second_order_dIdV(ks, EV, U=0.2)
    assert np.max(np.abs(got - ref)) < 1e-13*np.max(ref)


def test_second_order_default_drops_only_empty_states():
    ks = _chain(1.0)
    assert len(cd._occupied_states(ks)) < ks.dim
    ref = _old_second(ks, EV)
    got = cd.second_order_dIdV(ks, EV)
    assert np.max(np.abs(got - ref)) < 1e-10*np.max(ref)


def test_second_order_onehot_state_memory():
    """A TASK-shaped spectrum (p on one state) of dim 400 and 61 biases must
    not form the (n_eV, dim, dim) tensors (about 160 MB each)."""
    rng = np.random.default_rng(1)
    dim = 400
    ks = _chain(0.0)
    ks.dim = dim
    ks.e = np.sort(rng.random(dim))*1e-2
    ks.e -= ks.e[0]
    ks.p = np.zeros(dim); ks.p[0] = 1.0
    A = rng.random((dim, dim)) + 1j*rng.random((dim, dim))
    ks.Sx = ks.Sy = ks.Sz = A
    eV = np.linspace(-5e-3, 5e-3, 61)
    tracemalloc.start()
    cd.second_order_dIdV(ks, eV)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert peak < 30e6


@pytest.mark.parametrize("T,cplx", [(0.0, False), (3.0, False), (12.0, True)])
@pytest.mark.parametrize("coeff_bytes,exchange_bytes", [(1, 1), (16*64*3, 16*64*5)])
def test_third_order_blocked_equals_single_pass(monkeypatch, T, cplx,
                                                 coeff_bytes, exchange_bytes):
    ks = _chain(T, complex_h=cplx)
    ref = _old_third(ks, EV, -0.05, p_cut=0)
    monkeypatch.setattr(cd, "COEFF_BYTES", coeff_bytes)
    monkeypatch.setattr(cd, "EXCHANGE_BYTES", exchange_bytes)
    got = cd.third_order_kondo_dIdV(ks, EV, -0.05, p_cut=0)
    assert np.max(np.abs(got - ref)) < 1e-12*np.max(np.abs(ref))


def test_third_order_coefficient_memory_follows_the_budget(monkeypatch):
    """All 256 states occupied (n=8, high T): the single-pass coefficient
    tensor is 16*256^3 = 268 MB complex, and the old code peaks above that.
    With a 20 MB block the peak stays under 100 MB and the numbers agree."""
    ks = _chain(1000.0, n=8)
    assert len(cd._occupied_states(ks)) == ks.dim == 256
    eV = np.linspace(-5e-3, 5e-3, 3)
    monkeypatch.setattr(cd, "COEFF_BYTES", 20e6)
    monkeypatch.setattr(cd, "EXCHANGE_BYTES", 60e6)
    Fb = cd.FBuilder(ks.T, omega0=20e-3, Gamma0=5e-6, kB=ks.kB)
    tracemalloc.start()
    got = cd.third_order_kondo_dIdV(ks, eV, -0.05, p_cut=0, Fb=Fb)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert peak < 100e6
    ref = _old_third(ks, eV, -0.05, p_cut=0)
    assert np.max(np.abs(got - ref)) < 1e-12*np.max(np.abs(ref))
