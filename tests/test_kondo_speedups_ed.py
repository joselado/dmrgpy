"""The speedups of the ED Kondo route (docs/kondo_speedup_plan.md, items 1
to 6). Every test pins the numbers against the form of the code before the
change, kept here, and the saving, as a count of F or Theta entries, of
coefficient entries or of diagonalizations, so that undoing the saving fails
a test."""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.kondospectrumtk import conductance as cd
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum

_EPS3 = np.zeros((3, 3, 3))
for _i, _j, _k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    _EPS3[_i, _j, _k] = 1.
for _i, _j, _k in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
    _EPS3[_i, _j, _k] = -1.


def _old_coefficients(ks, occ):
    """Im of the Levi-Civita triple product over occ x every f x every m,
    in complex arithmetic, as built before the speedups."""
    Xi = np.stack([np.asarray(ks.Sx, dtype=complex), np.asarray(ks.Sy, dtype=complex),
                   np.asarray(ks.Sz, dtype=complex)], axis=-1)
    return np.imag(np.einsum('jkl,ifl,fmk,mij->ifm', _EPS3, Xi[occ], Xi,
                             Xi[:, occ], optimize=True))


def _old_second(ks, eVs, T0=1.0, U=0.0):
    """Second order over the occupied initial states and every final state."""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks)
    S = [np.asarray(a, dtype=complex)[:, occ] for a in (ks.Sx, ks.Sy, ks.Sz)]
    M2 = (0.5*np.abs(S[0] - 1j*S[1])**2 + 0.5*np.abs(S[0] + 1j*S[1])**2
          + np.abs(S[2])**2)
    eps = ks.e[:, None] - ks.e[occ][None, :]
    w = ks.p[occ][None, :]*M2
    spin = (np.einsum('fi,efi->e', w, cd._theta_raw(ks, eVs[:, None, None] - eps))
            + np.einsum('fi,efi->e', w, cd._theta_raw(ks, -eVs[:, None, None] - eps)))
    return 2*np.pi*T0**2*(spin + 4*U**2)


def _old_third(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6, Fb=None):
    """Third-order Kondo term in one pass over every f and m, both F tables
    of both diagrams evaluated in full, as it was before the speedups."""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks)
    p = ks.p[occ]
    coeff = _old_coefficients(ks, occ)/2.
    eps = ks.e[None, :] - ks.e[occ, None]
    e_fm = ks.e[:, None] - ks.e[None, :]
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    dim, nocc = ks.dim, len(occ)

    def one(v):
        nv = len(v)
        Th = cd._theta_raw(ks, v[:, None, None] - eps[None])
        Fd = Fcall((v[:, None, None] - eps[None]).ravel()).reshape(nv, nocc, dim)
        Fx = Fcall((v[:, None, None] - e_fm[None]).ravel()).reshape(nv, dim, dim)
        return (np.einsum('i,ifm,eif,eim->e', p, coeff, Th, Fd, optimize=True)
                + np.einsum('i,ifm,eif,efm->e', p, coeff, Th, Fx, optimize=True))
    return 4*np.pi*T0**2*Jrho_s*(one(eVs) + one(-eVs))


def _old_potential(ks, eVs, Jrho_s, U, T0=1.0, omega0=20e-3, Gamma0=5e-6,
                   Fb=None):
    """Third-order potential term with four F arrays and no bias chunks."""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks)
    p = ks.p[occ]
    Xi = np.stack([np.asarray(a, dtype=complex) for a in (ks.Sx, ks.Sy, ks.Sz)], axis=-1)
    loop = np.real(np.einsum('imk,mik->im', Xi[occ], Xi[:, occ]))
    eps = ks.e[None, :] - ks.e[occ, None]
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    nocc, dim = len(occ), ks.dim

    def one(v):
        Fim = Fcall((v[:, None, None] - eps[None]).ravel()).reshape(len(v), nocc, dim)
        Fmi = Fcall((v[:, None, None] + eps[None]).ravel()).reshape(len(v), nocc, dim)
        return np.einsum('i,im,eim->e', p, loop, Fim - Fmi)*cd._theta_raw(ks, v)
    return 4*np.pi*T0**2*Jrho_s*U*(one(eVs) - one(-eVs))


def _chain(n=5, field=True, complex_h=False):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n - 1):
        h = h + 3e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                      + 0.7*sc.Sz[i]*sc.Sz[i+1])
    if field: h = h + 0.4e-3*sc.Sz[0] + 0.3e-3*sc.Sx[2]
    if complex_h: h = h + 0.5e-3*sc.Sy[1]
    sc.set_hamiltonian(h)
    return sc


EV = np.linspace(-2e-3, 2e-3, 9)


def _rel(a, b):
    return np.max(np.abs(np.asarray(a) - np.asarray(b)))/np.max(np.abs(b))


class _CountF:
    """A callable F that counts the entries it is asked for."""
    def __init__(self, f):
        self.f, self.n = f, 0
        for name in ("odd_part",):
            if hasattr(f, name): setattr(self, name, getattr(f, name))
    def __call__(self, x):
        self.n += np.size(x)
        return self.f(x)


cd_F0 = cd.F0


def _patch_F0(monkeypatch):
    calls = []
    def counting(x, **kw):
        calls.append(np.size(x))
        return cd_F0(x, **kw)
    monkeypatch.setattr(cd, "F0", counting)
    return calls


# item 1: final states closed by the step


@pytest.mark.parametrize("T", [0.0, 0.5])
def test_final_states_cut_keeps_the_numbers(T):
    ks = KondoSpectrum(_chain(), 0, T)
    occ = cd._occupied_states(ks)
    fin = cd._final_states(ks, occ, EV)
    assert len(fin) < ks.dim # the cut is not empty on this chain
    tol = 1e-13 if T == 0. else 1e-12
    assert _rel(cd.second_order_dIdV(ks, EV, U=0.2), _old_second(ks, EV, U=0.2)) < tol
    Fb = cd._get_F(ks, 20e-3, 5e-6, None)
    assert _rel(cd.third_order_kondo_dIdV(ks, EV, -0.05, Fb=Fb if T else None),
                _old_third(ks, EV, -0.05, Fb=Fb if T else None)) < tol


def test_final_states_cut_at_zero_temperature_drops_only_closed_steps():
    """Every state left out has e_f - e_i > max|eV| for every occupied i, so
    its Theta0 is exactly 0; and the T=0 F count is that of the kept ones."""
    ks = KondoSpectrum(_chain(), 0, 0.)
    occ = cd._occupied_states(ks)
    fin = cd._final_states(ks, occ, EV)
    out = np.setdiff1d(np.arange(ks.dim), fin)
    eps = ks.e[None, out] - ks.e[occ, None]
    assert np.all(cd._theta_raw(ks, np.max(np.abs(EV)) - eps) == 0.)
    assert np.all(cd._theta_raw(ks, np.max(np.abs(EV)) - (ks.e[None, fin]
                  - ks.e[occ, None])).max(axis=0) > 0.)


def test_final_states_cut_counts_F_entries(monkeypatch):
    """The saving: n_f*dim exchange entries per bias and direction, not dim^2."""
    ks = KondoSpectrum(_chain(), 0, 0.)
    occ = cd._occupied_states(ks)
    nf, nocc, dim = len(cd._final_states(ks, occ, EV)), len(occ), ks.dim
    calls = _patch_F0(monkeypatch)
    cd.third_order_kondo_dIdV(ks, EV, -0.05)
    assert sum(calls) <= 2*len(EV)*(nocc + nf)*dim
    assert sum(calls) < 2*len(EV)*(nocc + dim)*dim


def test_final_states_cut_on_an_unsorted_spectrum():
    """The TASK and Krylov spectra need not be sorted: the cut is a mask."""
    ks = KondoSpectrum(_chain(), 0, 0.5)
    ref2 = cd.second_order_dIdV(ks, EV)
    ref3 = cd.third_order_kondo_dIdV(ks, EV, -0.05)
    perm = np.random.default_rng(3).permutation(ks.dim)
    ks.e, ks.p = ks.e[perm], ks.p[perm]
    for name in ("Sx", "Sy", "Sz"):
        setattr(ks, name, getattr(ks, name)[np.ix_(perm, perm)])
    assert _rel(cd.second_order_dIdV(ks, EV), ref2) < 1e-13
    assert _rel(cd.third_order_kondo_dIdV(ks, EV, -0.05), ref3) < 1e-12
