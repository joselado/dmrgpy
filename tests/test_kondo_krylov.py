"""get_kondo_spectrum(mode="ED", submode="KRYLOV"): the thermal Krylov route
of the third-order spectrum (kondospectrumtk/krylovkondo.py). The initial
states are the low-lying eigenstates from a sparse solver, and each one
gets one block Krylov space started from S_x|i>, S_y|i>, S_z|i>. Once every
Krylov space is exhausted and the weight cut is below the one of the ED
sums (conductance.P_CUT = 1e-12), the result must be the full ED one, at
T=0 and T>0, with and without Sz conservation, for real and complex H.
See docs/new_methods/kondo_krylov_thermal/README.md."""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain
from dmrgpy.kondospectrumtk import krylovkondo

EV = np.linspace(-4e-3, 4e-3, 21) + 1.3e-5 # off the thresholds and off zero
KW = dict(site=1, Jrho_s=-0.08, U=0.3, omega0=20e-3, Gamma0=5e-5, order=3)


def heisenberg(n, field=None):
    """A spin-1/2 chain with uneven exchange around 1 meV. field: None
    (SU(2)), "xz" (a 120 degree pattern in the xz plane: H real, Sz not
    conserved) or "y" (a field along y: H complex)."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        J = 1e-3*(1.0 + 0.3*np.cos(i))
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        if field == "xz":
            h = h + 2e-4*(np.cos(2*np.pi*i/3)*sc.Sx[i] + np.sin(2*np.pi*i/3)*sc.Sz[i])
        elif field == "y":
            h = h + 2e-4*np.cos(i)*sc.Sy[i] + 1e-4*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc


def impurity_chain(nb=3):
    """An S=1 impurity with easy-plane anisotropy in a tilted field,
    coupled to a short S=1/2 chain."""
    sc = spinchain.Spin_Chain(["S=1"] + ["S=1/2"]*nb)
    h = 1e-3*sc.Sz[0]*sc.Sz[0] + 2.3e-4*sc.Sz[0] + 0.8e-4*sc.Sx[0]
    for i in range(nb):
        J = 0.6e-3 if i == 0 else 1e-3
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    return sc


def compare(sc, T, nkry=40, p_cut=1e-14, su2=False, **kw):
    kw = {**KW, **kw}
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", T=T, **kw)
    _, kr = sc.get_kondo_spectrum(EV, mode="ED", submode="KRYLOV", T=T, nkry=nkry,
                                  p_cut=p_cut, su2=su2, **kw)
    return np.max(np.abs(kr - ed))/np.max(np.abs(ed))


@pytest.mark.parametrize("n", [6, 5])
def test_heisenberg_at_zero_temperature(n):
    """n=5 has a ground doublet, averaged with equal weights as on ED."""
    assert compare(heisenberg(n), 0.) < 1e-12 # measured 1e-15


@pytest.mark.parametrize("field,T", [("xz", 2.0), ("xz", 0.), ("y", 1.0)])
def test_without_sz_conservation(field, T):
    """A field in the xz plane keeps H real, a field along y makes it
    complex; neither conserves Sz."""
    sc = heisenberg(5, field)
    assert krylovkondo.KrylovKondoSpectrum(sc, 1, T, nkry=4).real == (field == "xz")
    # measured 1e-15 (xz) and 1e-12 (y), the size of the weight cut of the
    # ED sums (conductance.P_CUT), which keep fewer initial states
    assert compare(sc, T) < 1e-11


@pytest.mark.parametrize("order,T", [(2, 0.), (3, 0.), (3, 1.5)])
def test_anisotropic_spin_one(order, T):
    assert compare(impurity_chain(), T, site=0, order=order) < 1e-11 # measured 3e-15


def test_default_weight_cut_is_at_its_own_size():
    """p_cut=1e-8 drops the initial states below that relative weight, and
    dI/dV moves by about that much."""
    err = compare(heisenberg(5, "xz"), 2.0, p_cut=krylovkondo.P_CUT)
    assert 1e-12 < err < 1e-7 # measured 6e-9


def test_sparse_eigensolver_finds_the_multiplets(monkeypatch):
    """Above DENSE_DIM the initial states come from eigsh with k raised
    until the weight cut; on an SU(2) chain every level is a multiplet,
    and a missing member would show here."""
    monkeypatch.setattr(krylovkondo, "DENSE_DIM", 16)
    monkeypatch.setattr(krylovkondo, "K0", 6)
    sc = heisenberg(8)
    kks = krylovkondo.KrylovKondoSpectrum(sc, 1, 0.5, p_cut=1e-12)
    assert 1 < len(kks.weights) < 40 # eigsh, not the dense fallback
    assert compare(sc, 0.5, nkry=60, p_cut=1e-12) < 1e-10 # measured 2e-13


def test_su2_reduction_matches_ed_and_refuses_a_field():
    sc = heisenberg(7)
    for T in (0., 1.5):
        assert compare(sc, T, su2=True) < 1e-11 # measured 1e-15 and 4e-12
    kks = krylovkondo.KrylovKondoSpectrum(sc, 1, 1.5, p_cut=1e-12, su2=True)
    # one member per multiplet, in the sector M=1/2 of the odd chain
    assert np.all(np.mod(2*kks.spins, 2) == 1) and 1.5 in kks.spins
    assert len(kks.weights) < len(krylovkondo.KrylovKondoSpectrum(
        sc, 1, 1.5, p_cut=1e-12).weights)
    with pytest.raises(ValueError, match="does not commute with the total spin"):
        krylovkondo.KrylovKondoSpectrum(heisenberg(5, "xz"), 1, 1.0, su2=True)


def test_three_separate_spaces_double_count_with_a_field(monkeypatch):
    """The trap the block space avoids: three Lanczos spaces from S+|i>,
    S-|i> and Sz|i>, each orthonormal on its own, overlap once Sz is not
    conserved, and the shared directions are counted more than once. With
    Sz conserved they lie in different sectors and agree with ED."""
    block = krylovkondo.block_krylov

    def separate(H, starts, psi, nkry):
        out = [block(H, starts[:, [a]], psi, nkry) for a in range(starts.shape[1])]
        return tuple(np.concatenate([o[k] for o in out], axis=-1 if k == 1 else 0)
                     for k in range(3))

    assert compare(heisenberg(5, "xz"), 0.5) < 1e-11 # measured 5e-13
    monkeypatch.setattr(krylovkondo, "block_krylov", separate)
    assert compare(heisenberg(5), 0.5) < 1e-8 # measured 5e-11
    assert compare(heisenberg(5, "xz"), 0.5) > 0.1 # measured 2.7


def test_hubbard_site_in_an_electron_number_sector():
    """The spin of a Hubbard orbital, with charge fluctuations, in the Nf
    sector of the half-filled chain (S_site conserves the electron
    number)."""
    fc = fermionchain.Spinful_Fermionic_Chain_Native(3)
    h = 0
    for i in range(2):
        h = h - (fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + h.get_dagger()
    for i in range(3):
        h = h + 8.*(fc.Nup[i] - 0.5)*(fc.Ndn[i] - 0.5) + 0.1*fc.Sz[i] + 0.05*fc.Sx[i]
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=3)
    eV = np.linspace(-0.6, 0.6, 9) + 1e-3
    kw = dict(site=0, Jrho_s=-0.05, U=0.2, T=0.3, kB=1., omega0=4., Gamma0=1e-2)
    _, ed = fc.get_kondo_spectrum(eV, **kw)
    _, kr = fc.get_kondo_spectrum(eV, submode="KRYLOV", nkry=20, p_cut=1e-14, **kw)
    assert np.max(np.abs(kr - ed)) < 1e-11*np.max(np.abs(ed)) # measured 1e-13


def test_keywords():
    sc = heisenberg(4)
    kw = dict(site=0, Jrho_s=-0.05, T=0.)
    with pytest.raises(TypeError, match=r"argument\(s\): nkry\. "):
        sc.get_kondo_spectrum(EV, nkry=10, **kw) # needs submode="KRYLOV"
    with pytest.raises(TypeError, match=r"argument\(s\): delta\. "):
        sc.get_kondo_spectrum(EV, submode="KRYLOV", delta=1e-5, **kw)
    with pytest.raises(TypeError, match=r"argument\(s\): submode\. "):
        sc.get_kondo_spectrum(EV, submode="TASK", **kw) # a mode="DMRG" submode
    with pytest.raises(NotImplementedError, match="full eigenbasis"):
        sc.get_kondo_spectrum(EV, submode="KRYLOV", order="LL", **kw)
    for bad in (dict(nkry=0), dict(nkry=2.5), dict(p_cut=0.)):
        with pytest.raises(ValueError):
            sc.get_kondo_spectrum(EV, submode="KRYLOV", **bad, **kw)
