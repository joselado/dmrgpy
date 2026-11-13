"""A spin chain on a normal metal with the Wilson chain on any spin
(ysr.spin_chain_on_metal): the exact reductions, the mirror symmetry
between the two end spins (a plain chain against a Kondo term that jumps
the whole spin chain), and the middle spin against ED on "python" and 3."""
import numpy as np
import pytest

from dmrgpy import fermionchain, wilsonchain, ysr
from dmrgpy.ysrtk import paritystate

_ES, _DELTA = np.linspace(-1.2, 1.2, 13), 0.06
_KW = dict(Js=0.5, JK=1.0)
_ED = {}


def _mode(v):
    return ("ED", 3) if v == "ED" else ("DMRG", v)


def _ed(site):
    """ED energy and spectrum of 3 spins on 3 Wilson sites, once per site."""
    if site not in _ED:
        fc = ysr.spin_chain_on_metal(3, 3, site=site, **_KW)
        e = fc.gs_energy(mode="ED")
        d = ysr.get_stm_spectrum(fc, _ES, delta=_DELTA, mode="ED",
                                 submode="CVM")[1]
        _ED[site] = (e, d)
    return _ED[site]


def _heisenberg_energy(ns, Js, B):
    sx = np.array([[0, 1], [1, 0]])/2
    sy = np.array([[0, -1j], [1j, 0]])/2
    sz = np.diag([0.5, -0.5])
    def op(o, i):
        m = np.eye(1)
        for k in range(ns):
            m = np.kron(m, o if k == i else np.eye(2))
        return m
    H = np.zeros((2**ns, 2**ns), complex)
    for i in range(ns - 1):
        for o in (sx, sy, sz):
            H += Js*op(o, i) @ op(o, i + 1)
    for i in range(ns):
        H += B*op(sz, i)
    return np.linalg.eigvalsh(H)[0]


def test_end_spin_is_island_chain_without_pairing():
    """The default site, the end spin, is island_chain with one island
    site taken as the first Wilson site and Delta=0, term by term."""
    ns, n, Upen = 3, 3, 5.0
    eps, t = wilsonchain.wilson_chain(n)
    a = ysr.spin_chain_on_metal(ns, n, Upen=Upen, **_KW)
    b = ysr.island_chain(ns, 1, n - 1, Delta=0.0, eps_island=eps[0],
                         t_c=t[0], wilson=(eps[1:], t[1:]), Upen=Upen, **_KW)
    assert a.tip_spin == ns - 1
    assert (a.hamiltonian - b.hamiltonian).is_zero()


def test_one_spin_is_kondo_chain():
    a = ysr.spin_chain_on_metal(1, 4, JK=0.4)
    b = ysr.kondo_chain(4, J=0.4, Delta=0.0)
    assert a.gs_energy(mode="ED") == pytest.approx(b.gs_energy(mode="ED"),
                                                   abs=1e-12)


@pytest.mark.parametrize("site,B", [(0, 0.0), (1, 0.07)])
def test_decoupled_energy_is_heisenberg_plus_free_metal(site, B):
    ns, n, Js = 3, 3, 0.5
    fc = ysr.spin_chain_on_metal(ns, n, site=site, Js=Js, JK=0.0, B=B)
    eps, t = wilsonchain.wilson_chain(n)
    e = np.linalg.eigvalsh(np.diag(eps) + np.diag(t, 1) + np.diag(t, -1))
    exact = _heisenberg_energy(ns, Js, B) + 2.0*np.sum(e[e < 0])
    assert fc.gs_energy(mode="ED") == pytest.approx(exact, abs=1e-10)


def test_the_two_end_spins_are_mirror_images():
    """site=0, whose Kondo term jumps over the whole spin chain, against
    site=ns-1, a plain chain: the uniform chain is mirror symmetric, so
    the energy and the spectrum on the coupled spin are the same."""
    (e0, d0), (e2, d2) = _ed(0), _ed(2)
    assert e0 == pytest.approx(e2, abs=1e-10)
    assert np.max(np.abs(d0 - d2)) < 1e-8
    assert abs(_ed(1)[0] - e2) > 1e-3 # the middle spin is another model


def test_middle_spin_against_the_same_junction_in_another_site_order():
    """The T junction written by hand with the Wilson chain first, then
    the middle spin, then the two outer spins: the same energy and
    spectrum on ED, so the builder's long-range term is that junction."""
    ns, n, Upen = 3, 3, 2.0*(1.0 + 1.0 + 0.5)
    eps, t = wilsonchain.wilson_chain(n)
    f, s = [2, 1, 0], [5, 3, 4] # positions of f_0..f_2 and of spins 0..2
    fc = fermionchain.Spinful_Fermionic_Chain_Native(ns + n)
    dot = lambda a, b: (fc.Sx[a]*fc.Sx[b] + fc.Sy[a]*fc.Sy[b]
                        + fc.Sz[a]*fc.Sz[b])
    h = 0.5*(dot(s[0], s[1]) + dot(s[1], s[2])) + 1.0*dot(s[1], f[0])
    for a in s:
        q = fc.Nup[a] + fc.Ndn[a] - 1
        h = h + Upen*q*q
    for k in range(n - 1):
        h = h + t[k]*ysr._hop(fc, f[k], f[k+1])
    fc.set_hamiltonian(h)
    fc.tip_spin, fc.channels = s[1], [f]
    e, d = _ed(1)
    assert fc.gs_energy(mode="ED") == pytest.approx(e, abs=1e-10)
    dh = ysr.get_stm_spectrum(fc, _ES, delta=_DELTA, mode="ED",
                              submode="CVM")[1]
    assert np.max(np.abs(dh - d)) < 1e-8


@pytest.mark.parametrize("version", ["python", 3])
@pytest.mark.parametrize("site", [0, 1])
def test_energy_against_ed(version, site):
    """6 sites saturate at maxm=64, so DMRG is exact here; the state holds
    one electron on every spin site."""
    fc = ysr.spin_chain_on_metal(3, 3, site=site, itensor_version=version,
                                 **_KW)
    fc.maxm, fc.nsweeps = 64, 20
    assert fc.gs_energy() == pytest.approx(_ed(site)[0], abs=1e-8)
    q = sum((fc.Nup[a] + fc.Ndn[a] - 1)*(fc.Nup[a] + fc.Ndn[a] - 1)
            for a in fc.impurity)
    assert abs(fc.vev(q)) < 1e-8


@pytest.mark.parametrize("version", ["python", 3])
def test_stm_spectrum_on_the_middle_spin_against_ed(version):
    """TASK against ED with the tip on the middle spin, whose bath site
    is two sites away. The 6 sites saturate at maxm=64, so the tangent
    space is the whole space: 4.4e-12 ("python") and 2.6e-12 (3) on a
    peak of 0.442, with and without pad_bonds="full"."""
    fc = ysr.spin_chain_on_metal(3, 3, site=1, itensor_version=version,
                                 **_KW)
    fc.maxm, fc.nsweeps = 64, 20
    _, d = ysr.get_stm_spectrum(fc, _ES, delta=_DELTA)
    assert np.max(np.abs(d - _ed(1)[1])) < 1e-8


def test_default_vertex_is_the_coupled_spin():
    fc = ysr.spin_chain_on_metal(3, 2, site=1, **_KW)
    kw = dict(delta=_DELTA, mode="ED", submode="CVM")
    _, a = ysr.get_stm_spectrum(fc, _ES, **kw)
    _, b = ysr.get_stm_spectrum(fc, _ES, imp=1, site=3, **kw)
    _, c = ysr.get_stm_spectrum(fc, _ES, imp=0, **kw)
    assert np.max(np.abs(a - b)) < 1e-12
    assert np.max(np.abs(a - c)) > 1e-2*np.max(a)


@pytest.mark.parametrize("site", [3, -1, 1.5])
def test_a_site_that_is_not_a_spin_raises(site):
    with pytest.raises(ValueError, match="site must be one of the spins"):
        ysr.spin_chain_on_metal(3, 3, site=site)


def test_sector_tools_refuse_a_gapless_metal():
    fc = ysr.spin_chain_on_metal(2, 2)
    with pytest.raises(ValueError, match="needs a gap"):
        ysr.subgap_energies(fc, mode="ED")
    with pytest.raises(ValueError, match="island_chain"):
        ysr.island_energies(fc, mode="ED")
    with pytest.raises(ValueError, match="needs a gap"):
        paritystate.parity_ground_state(fc, mode="ED")
