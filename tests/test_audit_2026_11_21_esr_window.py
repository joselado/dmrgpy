"""The 2026-11-21 hole hunt, cluster esr-window
(docs/audit_2026_11_21_hole_hunt.md, findings 3, 4 and 5).

What is pinned: on mode="ED" the kept states follow the drive, so the state
a second photon reaches from the upper state of a driven line is kept (a
spin 1 in a field at 0 and 50 mK equals a window that keeps all three
levels and a hand-written three-level rate equation); the band kept above
the window is counted in the power-broadened width of the strongest line,
so a saturated line just above a window at T = 0 is the Bloch line there;
the DMRG eigenstate route, which does not grow its band, warns when the
drive needs more states; and T1 and T2 are those of the tip spin whatever
its weight in site_weights, so only the ratios of the weights matter."""
import warnings

import numpy as np
import pytest

from dmrgpy import spinchain

H = 4.135667696e-6 # eV/GHz
KB = 8.617333262e-5 # eV/K
F0, TH, T1, T2 = 20.0, 0.7, 100., 50.
TIP = [np.sin(TH), 0., np.cos(TH)]


def rate_eq(f, T, rabi, nlev, d=0.):
    """The rate equation of a spin 1 in a field F0 (plus h d Sz^2) with its
    lowest nlev levels, written out by hand: Korringa eq. 42 as in
    esrtk/rates.py, the drive rates of esrtk/lineshape.py, and the DC
    readout through D = V = m.S (the repro of finding 3, verbatim)"""
    ms = np.array([-1., 0., 1.])[:nlev]
    E = H*(F0*ms + d*ms**2)
    E = E - np.min(E)
    kT = KB*T
    n = nlev
    Sp = np.zeros((n, n))
    for i in range(n - 1): # <m+1|S+|m> of a spin 1
        Sp[i + 1, i] = np.sqrt(2. - ms[i]*(ms[i] + 1.))
    Sx = 0.5*(Sp + Sp.T)
    Sz = np.diag(ms)
    K = 0.5*(Sp**2 + (Sp**2).T) # |Sx|^2 + |Sy|^2, Sz has no off-diagonal
    eps0 = H*F0
    c = 1. if kT == 0. else 1./np.tanh(eps0/(2.*kT))
    g = 2./(T1*eps0*c)
    def phi(x):
        if kT == 0.: return max(x, 0.)
        return x/(-np.expm1(-x/kT))
    Gam = np.zeros((n, n)) # Gam[m, k]: rate from k to m
    for m in range(n):
        for k in range(n):
            if m != k: Gam[m, k] = g*K[m, k]*phi(E[k] - E[m])
    out = Gam.sum(axis=0)
    gphi = 1./T2 - 0.5/T1
    V = np.sin(TH)*Sx + np.cos(TH)*Sz
    Om = 4.*np.pi*rabi
    w = 2.*np.pi*f
    R = Gam.copy()
    R[np.diag_indices(n)] = -out
    for l in range(n - 1):
        u = l + 1
        gam = 0.5*(out[u] + out[l]) + gphi
        wul = 2.*np.pi*(E[u] - E[l])/H
        W = 0.5*Om**2*V[u, l]**2*gam/((w - wul)**2 + gam**2)
        R[u, l] += W; R[l, u] += W; R[u, u] -= W; R[l, l] -= W
    A = np.vstack([R, np.ones(n)])
    b = np.zeros(n + 1); b[-1] = 1.
    p = np.linalg.lstsq(A, b, rcond=None)[0]
    if kT == 0.:
        p0 = np.zeros(n); p0[0] = 1.
    else:
        p0 = np.exp(-E/kT); p0 /= p0.sum()
    return float(np.sum((p - p0)*np.diag(V)))


def bloch(f, T, rabi, r=0., phase=0., T1=T1, T2=T2):
    """One spin 1/2, field along z, tip at TH from it, V = D = m.S: the DC
    part and the homodyne part (the repro of finding 4)"""
    OR = 2*np.pi*rabi*np.sin(TH)
    xi2 = OR**2*T1*T2
    x = 2*np.pi*(f - F0)
    y0 = -np.tanh(H*F0/(2*KB*T)) if T > 0 else -1.
    y = y0*(1 + (x*T2)**2)/(1 + (x*T2)**2 + xi2)
    dc = np.cos(TH)*(y - y0)/2
    rho = -(2*np.pi*rabi)*(np.sin(TH)/2)*y/(x + 1j/T2)
    return dc + r*np.real(np.exp(1j*phase)*rho*np.sin(TH)/2)


def free_spin():
    sc = spinchain.Spin_Chain(["S=1/2"])
    sc.set_hamiltonian(H*F0*sc.Sz[0])
    return sc


def spin1():
    """A spin 1 in a field: two rungs, both at F0"""
    sc = spinchain.Spin_Chain(["S=1"])
    sc.set_hamiltonian(H*F0*sc.Sz[0])
    return sc


FA = np.linspace(19.9, 20.1, 41) # the window of the line
FC = np.append(FA, 2*F0 + 0.5) # plus one point that keeps |+1>


@pytest.mark.parametrize("T", [0.0, 0.05])
@pytest.mark.parametrize("rabi", [1e-3, 0.05])
def test_spin1_ladder_keeps_the_second_rung(T, rabi):
    """Finding 3: below about 70 mK kT ln(1/p_cut) does not reach |+1>, 40
    GHz up, and the old band dropped the rung |0> -> |+1> (the line was
    4.3e-2 low at the default drive, 50% at rabi=0.05). The kept band
    follows the drive now: window A keeps both rungs, equals window C, which
    keeps all three levels anyway, on the rates and on the Floquet route,
    and is the three-level rate equation, not the two-level one"""
    sc = spin1()
    kw = dict(tip=TIP, T=T, T1=T1, T2=T2, f_ref=F0, rabi=rabi)
    _, a, pa = sc.get_esr_spectrum(FA, return_poles=True, **kw)
    _, c = sc.get_esr_spectrum(FC, **kw)
    c = c[:len(FA)]
    assert len(pa["frequency"]) == 2
    pk = np.max(np.abs(c))
    assert np.max(np.abs(a - c)) < 1e-12*pk
    i = int(np.argmin(np.abs(FA - F0)))
    assert a[i] == pytest.approx(rate_eq(FA[i], T, rabi, 3), rel=1e-8)
    assert abs(a[i] - rate_eq(FA[i], T, rabi, 2)) > 1e-2*abs(a[i])
    _, qa = sc.get_esr_spectrum(FA, method="floquet", nharm=3, **kw)
    _, qc = sc.get_esr_spectrum(FC, method="floquet", nharm=3, **kw)
    qc = qc[:len(FA)]
    assert np.max(np.abs(qa - qc)) < 1e-10*np.max(np.abs(qc))


@pytest.mark.parametrize("r", [0.0, 0.3])
def test_saturated_line_above_the_window_at_zero_temperature(r):
    """Finding 4: at T = 0 and xi = 14.3 (rabi = 0.05) a line 0.2 GHz above
    the window lay just past 30 free-spin widths, and the window returned 0
    at its top for 5.0e-2 of the line's peak (3.0e-2 with V_RF/V_DC = 0.3).
    The band is counted in the power-broadened width of the strongest line
    now, so the window is the Bloch line, the same as a window that reaches
    past the line"""
    sc = free_spin()
    fA = np.linspace(19.70, 19.80, 11)
    fB = np.linspace(19.70, 20.30, 61)
    kw = dict(tip=TIP, T=0., T1=T1, T2=T2, f_ref=F0, rabi=0.05,
              vrf_over_vdc=r, phase=0.)
    _, a = sc.get_esr_spectrum(fA, **kw)
    _, b = sc.get_esr_spectrum(fB, **kw)
    pk = np.max(np.abs(bloch(fB, 0., 0.05, r)))
    assert np.max(np.abs(a - bloch(fA, 0., 0.05, r))) < 1e-10*pk
    iB = int(np.argmin(np.abs(fB - fA[-1])))
    assert abs(a[-1] - b[iB]) < 1e-10*pk


def test_weak_drive_tail_left_out_is_below_the_bound():
    """At weak drive (xi = 0.03) the band is still 30 free-spin widths, and
    the DC tail of the line left out is below 1/(4 WINDOW_MARGIN^2) of its
    peak at the window top (measured 2.5e-4)"""
    from dmrgpy.esrtk import lineshape
    sc = free_spin()
    fA = np.linspace(19.70, 19.80, 11)
    _, a = sc.get_esr_spectrum(fA, tip=TIP, T=0., T1=T1, T2=T2, f_ref=F0,
                               rabi=1e-4)
    pk = abs(bloch(F0, 0., 1e-4))
    bound = 1./(4.*lineshape.WINDOW_MARGIN**2)
    assert np.max(np.abs(a - bloch(fA, 0., 1e-4))) < bound*pk


def ladder_pair():
    """The tip spin 1 at 20 GHz next to an uncoupled spin 1/2 at 7 GHz, on
    two MPS sites ("python": v3 takes ED below 3 sites): levels at 0, 7,
    20, 27, 40 and 47 GHz, two rungs of the spin 1 for each state of the
    spin 1/2, and the relaxation through the spin 1/2 connects the two
    ladders"""
    sc = spinchain.Spin_Chain(["S=1", "S=1/2"], itensor_version="python")
    sc.set_hamiltonian(H*F0*sc.Sz[0] + H*7.*sc.Sz[1])
    sc.maxm, sc.nsweeps = 8, 10
    return sc


WARN = "can take part under this drive"


def _warned(rec):
    return any(WARN in str(x.message) for x in rec)


def test_ed_keeps_the_second_rung_on_two_sites():
    """At 50 mK and the default drive the upper states of the 20 GHz rungs
    carry population, and the ED band reaches the second rungs at 40 and 47
    GHz: the window of the line equals a window that keeps every level"""
    np.random.seed(0)
    sc = ladder_pair()
    f = np.linspace(19.9, 20.1, 21)
    kw = dict(tip=TIP, T=0.05, T1=T1, T2=T2, f_ref=F0, mode="ED")
    _, a, poles = sc.get_esr_spectrum(f, return_poles=True, **kw)
    _, c = sc.get_esr_spectrum(np.append(f, 47.5), **kw)
    assert len(poles["frequency"]) == 4
    c = c[:len(f)]
    assert np.max(np.abs(a - c)) < 1e-12*np.max(np.abs(c))


def test_dmrg_route_warns_when_the_drive_needs_more_states():
    """The DMRG route keeps 30 free-spin widths above the window and the
    thermal states. It warns at 50 mK and the default drive (the rungs at
    40 and 47 GHz are dropped), and at T = 0 with a saturated line 0.2 GHz
    above the window (the band of its power-broadened width is not kept);
    with a weak drive and p_cut=1e-4 nothing is populated above the
    thermal states and it does not warn"""
    np.random.seed(0)
    sc = ladder_pair()
    kw = dict(tip=TIP, T1=T1, T2=T2, f_ref=F0, mode="DMRG", nstates=5)
    f = np.linspace(19.9, 20.1, 21)
    with pytest.warns(RuntimeWarning, match=WARN):
        sc.get_esr_spectrum(f, T=0.05, **kw)
    with pytest.warns(RuntimeWarning, match=WARN):
        sc.get_esr_spectrum(np.linspace(19.70, 19.80, 11), T=0., rabi=0.05,
                            **kw)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        sc.get_esr_spectrum(f, T=0.05, rabi=1e-7, p_cut=1e-4, **kw)
    assert not _warned(rec)


@pytest.mark.parametrize("w", [[0.5], [2.0]])
def test_site_weights_keep_the_tip_T1_and_T2(w):
    """Finding 5: on one free spin, site_weights=[0.5] and [2.0] gave a
    weak-drive width of 5.570 and 7.958 MHz for 1/(pi T2) = 6.366 MHz. The
    Korringa scale is divided by the weight of the tip site now, so the
    width is 1/(pi T2) and the saturated line is the Bloch line of the T1
    and T2 passed"""
    sc = free_spin()
    f = np.linspace(19.95, 20.05, 401)
    kw = dict(tip=TIP, T=1.2, T1=T1, T2=T2, f_ref=F0, site_weights=w)
    _, _, poles = sc.get_esr_spectrum(f, rabi=1e-4, return_poles=True, **kw)
    assert poles["width"][0] == pytest.approx(1./(np.pi*T2), rel=1e-12)
    _, dI = sc.get_esr_spectrum(f, rabi=5e-3, **kw)
    ref = bloch(f, 1.2, 5e-3)
    assert np.max(np.abs(dI - ref)) < 1e-10*np.max(np.abs(ref))


def dimer():
    """The dimer of the guide: J = 15 GHz, 25 GHz field, 4 GHz tip field on
    the first spin"""
    sd = spinchain.Spin_Chain(["S=1/2"]*2)
    hd = H*15.*(sd.Sx[0]*sd.Sx[1] + sd.Sy[0]*sd.Sy[1] + sd.Sz[0]*sd.Sz[1])
    hd = hd + H*25.*(sd.Sz[0] + sd.Sz[1]) + H*4.*sd.Sz[0]
    sd.set_hamiltonian(hd)
    return sd


def test_only_the_ratios_of_the_weights_matter():
    """Finding 5 on two spins: [0.5, 0.5] is the default (it was 1.20 of the
    default's peak off) and [2, 0.5] is [1, 0.25]. A tip site of weight 0
    keeps the scale of a site of weight 1, so weights c*w at (T1, T2) are
    w at T1/c and 1/T2' = 1/T2 - 1/(2 T1) + c/(2 T1), the identity every
    weight obeyed before"""
    sd = dimer()
    f = np.linspace(5., 45., 2001)
    kw = dict(tip=TIP, T=1.2, f_ref=25., rabi=1e-3)
    def dI(w, t1=T1, t2=T2):
        return sd.get_esr_spectrum(f, site_weights=w, T1=t1, T2=t2, **kw)[1]
    ref = dI(None)
    pk = np.max(np.abs(ref))
    assert np.max(np.abs(dI([0.5, 0.5]) - ref)) < 1e-12*pk
    b = dI([1.0, 0.25])
    assert np.max(np.abs(dI([2.0, 0.5]) - b)) < 1e-12*np.max(np.abs(b))
    c = 2.0
    t2p = 1./(1./T2 - 0.5/T1 + 0.5*c/T1)
    a0 = dI([0.0, 1.0], t1=T1/c, t2=t2p)
    assert np.max(np.abs(dI([0.0, c]) - a0)) < 1e-12*np.max(np.abs(a0))
