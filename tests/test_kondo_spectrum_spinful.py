"""get_kondo_spectrum on the spinful fermion chains: the tip couples to the
local electron spin c^dagger sigma c/2 of a Hubbard site, through the same
KondoSpectrumMixin as Spin_Chain."""
import numpy as np
import pytest
from dmrgpy import fermionchain, spinchain

KW = dict(site=0, Jrho_s=-0.05, T=0., kB=1., omega0=4., Gamma0=1e-2)
CLASSES = (fermionchain.Spinful_Fermionic_Chain,
           fermionchain.Spinful_Fermionic_Chain_Native)


def hubbard(cls, n, U, b=0., t=1., itensor_version="python"):
    fc = cls(n, itensor_version=itensor_version)
    h = 0
    for i in range(n-1):
        h = h - t*(fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + h.get_dagger()
    for i in range(n):
        h = h + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5) + b*fc.Sz[i]
    fc.set_hamiltonian(h)
    fc.maxm = 64 ; fc.nsweeps = 20 ; fc.cutoff = 0.
    return fc


def heisenberg(n, J, b=0.):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n): h = h + b*sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc


@pytest.mark.parametrize("cls", CLASSES)
def test_task_matches_ed_on_hubbard(cls):
    eV = np.linspace(-0.6, 0.6, 9)
    fc = hubbard(cls, 4, 10., b=0.12)
    _, ed = fc.get_kondo_spectrum(eV, **KW)
    _, tk = fc.get_kondo_spectrum(eV, mode="DMRG", submode="TASK", **KW)
    assert np.max(np.abs(tk-ed)) < 1e-6*np.max(np.abs(ed))


def test_interleaved_and_native_chains_agree():
    eV = np.linspace(-0.6, 0.6, 9)
    d = [cls_ed for cls_ed in (hubbard(c, 3, 8., b=0.1).get_kondo_spectrum(
            eV, U=0.2, **KW)[1] for c in CLASSES)]
    assert np.max(np.abs(d[0]-d[1])) < 1e-10


def test_electron_number_sector_is_exact_at_zero_temperature():
    eV = np.linspace(-0.6, 0.6, 9)
    fc = hubbard(CLASSES[1], 3, 8., b=0.1)
    _, full = fc.get_kondo_spectrum(eV, **KW)
    fc.set_conserved_sector(Nf=3)
    _, sector = fc.get_kondo_spectrum(eV, **KW)
    assert np.max(np.abs(full-sector)) < 1e-10


def test_spin_sector_raises():
    fc = hubbard(CLASSES[1], 2, 8.)
    fc.set_conserved_sector(Nf=2, Sz=0) # S+- leave it
    with pytest.raises(ValueError, match="does not conserve"):
        fc.get_kondo_spectrum(np.linspace(-1., 1., 5), **KW)


def test_hubbard_dimer_approaches_heisenberg_as_t_over_U_squared():
    """J = (sqrt(U^2+16t^2)-U)/2 is the exact singlet-triplet gap, so with
    the bias in units of J only the weights differ, by O((t/U)^2); the
    grid misses eV=J, where the step would be read by rounding"""
    x = np.linspace(-3., 3., 60)
    def spec(chain, J):
        kw = dict(KW, omega0=10*J, Gamma0=1e-2*J)
        return chain.get_kondo_spectrum(x*J, **kw)[1]
    ref = spec(heisenberg(2, 1.), 1.)
    dev = []
    for U in (20., 40.):
        J = (np.sqrt(U**2+16.)-U)/2
        dev.append(np.mean(np.abs(spec(hubbard(CLASSES[1], 2, U), J)-ref)))
    assert 0. < dev[1] < dev[0] < 0.05*np.mean(ref)
    assert dev[1]/dev[0] == pytest.approx(0.25, rel=0.1)
