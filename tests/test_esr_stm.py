"""get_esr_spectrum (esrstm.py, esrtk/): the ESR-STM line shape from the
rate equation between eigenstates with Korringa relaxation, read through
the magnetoresistance of the tip site (DC part) and the RF modulation of
the conductance (homodyne part). The plan and the derivation are in
docs/esr_stm_plan.md.

What is pinned: one spin 1/2 gives the Bloch formula with power
broadening and the homodyne asymmetry, and the steady state of the driven
Lindblad equation of lindblad.py in the rotating frame; the rate equation
agrees with the Floquet solution of the driven Lindblad equation without
the secular and rotating-wave approximations; the DC amplitude of one spin
goes as sin^2(theta) cos(theta) with the angle between the field and the
tip moment; the singlet to T0 line of a dimer has no DC part at zero tip
field; the four lines of the Ti dimer of Yang et al. (PRL 119, 227206,
Supplemental Material eqs. S7 to S14) carry the matrix elements of S14 and
the form factors sin^2(xi/2) and cos^2(xi/2); a zero field raises; the
sparse ED route and the MPS route agree with full ED; the saturated
solve agrees with the null vector of the dense rate matrix; the rates do
not depend on the window once f_ref is given; and the paths that must
raise do."""
import numpy as np
import pytest

from dmrgpy import cppext, fermionchain, spinchain
from dmrgpy.esrtk import lineshape, spectrum

H = 4.135667696e-6 # eV/GHz
KB = 8.617333262e-5 # eV/K


def free_spin(f0=20.0):
    sc = spinchain.Spin_Chain(["S=1/2"])
    sc.set_hamiltonian(H*f0*sc.Sz[0])
    return sc


def bloch(f, f0, th, T, T1, T2, rabi, r=0., phase=0.):
    """dI of one spin 1/2, field along z, tip moment at th from it, V = D"""
    OR = 2*np.pi*rabi*np.sin(th) # Rabi frequency of the line, rad/ns
    xi2 = OR**2*T1*T2
    x = 2*np.pi*(f - f0)
    y0 = -np.tanh(H*f0/(2*KB*T)) if T > 0 else -1.
    y = y0*(1 + (x*T2)**2)/(1 + (x*T2)**2 + xi2)
    dc = np.cos(th)*(y - y0)/2
    rho = -(2*np.pi*rabi)*(np.sin(th)/2)*y/(x + 1j/T2)
    return dc + r*np.real(np.exp(1j*phase)*rho*np.sin(th)/2)


@pytest.mark.parametrize("rabi", [1e-4, 5e-3])
@pytest.mark.parametrize("r", [0., 0.3])
def test_single_spin_is_bloch(rabi, r):
    """Weak drive and saturation (xi^2 = 4.5 at rabi=5e-3), DC and homodyne"""
    th, f = 0.7, np.linspace(19.95, 20.05, 41)
    _, dI = free_spin().get_esr_spectrum(f, tip=[np.sin(th), 0, np.cos(th)],
                                         T=1.2, T1=100., T2=50., f_ref=20.,
                                         rabi=rabi, vrf_over_vdc=r, phase=0.4)
    ref = bloch(f, 20.0, th, 1.2, 100., 50., rabi, r, 0.4)
    assert np.max(np.abs(dI - ref)) < 1e-10*np.max(np.abs(ref))


def test_single_spin_power_broadened_width():
    """FWHM = (1/pi T2) sqrt(1 + xi^2) of the DC line"""
    th, T1, T2, rabi = np.pi/2 - 0.3, 100., 50., 5e-3
    f = np.linspace(19.9, 20.1, 40001)
    _, dI = free_spin().get_esr_spectrum(f, tip=[np.sin(th), 0, np.cos(th)],
                                         T=1.2, T1=T1, T2=T2, f_ref=20.,
                                         rabi=rabi)
    above = f[dI > 0.5*np.max(dI)]
    xi2 = (2*np.pi*rabi*np.sin(th))**2*T1*T2
    assert above[-1] - above[0] == pytest.approx(np.sqrt(1 + xi2)/(np.pi*T2),
                                                 rel=1e-3)


def test_single_spin_against_lindblad_rotating_frame():
    """The saturated steady state against lindblad.py (mode="ED"): one spin
    in the rotating frame, H = (w0 - w) S^z + Omega_R S^x, jumps
    sqrt(Gamma_down) S^-, sqrt(Gamma_up) S^+ and sqrt(2/Tphi) S^z"""
    from dmrgpy.esrtk import rates
    from dmrgpy.lindblad import Lindblad_Chain
    f0, th, T, T1, T2, rabi = 20.0, 0.9, 0.5, 100., 40., 4e-3
    f = np.array([19.99, 19.997, 20.0, 20.004, 20.02])
    _, dI = free_spin(f0).get_esr_spectrum(f, tip=[np.sin(th), 0, np.cos(th)],
                                           T=T, T1=T1, T2=T2, f_ref=f0,
                                           rabi=rabi)
    kT, eps = KB*T, H*f0
    g = rates.korringa_scale(T1, eps, kT)
    down = g*0.5*rates.phi(np.array([eps]), kT)[0]
    up = g*0.5*rates.phi(np.array([-eps]), kT)[0]
    gphi = 1./T2 - 0.5/T1
    OR = 2*np.pi*rabi*np.sin(th)
    p = spinchain.Spin_Chain(["S=1/2"])
    Sp = p.Sx[0] + 1j*p.Sy[0]
    Sm = p.Sx[0] - 1j*p.Sy[0]
    jumps = [np.sqrt(down)*Sm, np.sqrt(up)*Sp, np.sqrt(2*gphi)*p.Sz[0]]
    sz0 = -0.5*np.tanh(eps/(2*kT))
    for fi, d in zip(f, dI):
        Hrot = 2*np.pi*(f0 - fi)*p.Sz[0] + OR*p.Sx[0]
        lc = Lindblad_Chain(p, jumps, H=Hrot)
        lam, res, obs = lc.get_spectrum(mode="ED")
        sz = np.real(obs(p.Sz[0], 0))
        assert d == pytest.approx(np.cos(th)*(sz - sz0), abs=1e-9)


@pytest.mark.parametrize("rabi", [1e-4, 3e-3])
def test_rates_against_floquet_on_a_dimer(rabi):
    """Two Heisenberg spins, lines at 19, 20 (two of them, sharing T0) and
    21 GHz; the rate equation against the Floquet steady state of the
    driven Lindblad equation: the difference is the secular and
    rotating-wave error, of order gamma/omega (measured 1e-3 and 7e-4)"""
    sc = dimer(1.0, 20., 0.5)
    f = np.linspace(18.5, 21.5, 301)
    kw = dict(tip=[np.sin(0.6), 0, np.cos(0.6)], T=1.2, T1=100., T2=20.,
              f_ref=20., rabi=rabi, vrf_over_vdc=0.2, phase=0.3)
    _, a = sc.get_esr_spectrum(f, **kw)
    _, b = sc.get_esr_spectrum(f, method="floquet", **kw)
    assert np.max(np.abs(a - b)) < 3e-3*np.max(np.abs(b))


def test_dc_amplitude_against_angle():
    """V = D = m.S: the DC peak of one spin is sin^2(th) cos(th)/4 times a
    constant at weak drive (at rabi=1e-5 the saturation, xi^2 = 2e-5, is
    already the whole of the difference), zero with the tip moment normal to the field and
    largest at tan^2(th) = 2"""
    f = np.array([20.0])
    ths = np.array([0.2, 0.6, np.arctan(np.sqrt(2)), 1.2, np.pi/2])
    peaks = []
    for th in ths:
        _, dI = free_spin().get_esr_spectrum(f, tip=[np.sin(th), 0, np.cos(th)],
                                             T=1.2, f_ref=20., rabi=1e-7)
        peaks.append(dI[0])
    peaks = np.array(peaks)
    ref = np.sin(ths)**2*np.cos(ths)
    assert abs(peaks[-1]) < 1e-12*np.max(np.abs(peaks))
    assert peaks[:-1]/ref[:-1] == pytest.approx(np.full(4, peaks[0]/ref[0]), rel=1e-6)
    assert np.argmax(peaks) == 2


def dimer(J, B, Bt, D=0., n_axis=None):
    """Two spins 1/2, all energies in GHz: J S1.S2 + B (S1z + S2z) + Bt S1z,
    plus the secular dipolar term D (3 S1z S2z - S1.S2) of Yang et al. S7"""
    sc = spinchain.Spin_Chain(["S=1/2"]*2)
    h = H*J*(sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1])
    h = h + H*B*(sc.Sz[0] + sc.Sz[1]) + H*Bt*sc.Sz[0]
    if D != 0.:
        h = h + H*D*(3*sc.Sz[0]*sc.Sz[1] - (sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1]
                                          + sc.Sz[0]*sc.Sz[1]))
    sc.set_hamiltonian(h)
    return sc


def test_singlet_T0_line_has_no_dc_part_at_zero_tip_field():
    """Bae et al. (arXiv:1807.04322): the S-T0 transition leaves the
    time-averaged spin of the tip site unchanged, so only the homodyne part
    sees it; with a tip field the eigenstates become product-like and the
    DC part appears"""
    f = np.linspace(0.9, 1.3, 2001)
    kw = dict(tip=[np.sin(0.6), 0, np.cos(0.6)], T=1.2, T1=100., T2=20.,
              f_ref=20., rabi=1e-4, vrf_over_vdc=0.2, return_parts=True,
              return_poles=True)
    _, _, parts, poles = dimer(1.0, 20., 0.0).get_esr_spectrum(f, **kw)
    dc0 = np.max(np.abs(parts["dc"]))
    st = np.argmin(np.abs(poles["frequency"] - 1.0))
    assert abs(poles["form_factor"][st]) < 1e-12
    # what is left of the DC part is the relaxation that carries the driven
    # population of S and T0 on to T- and T+ (measured 1.2e-3)
    assert np.max(np.abs(parts["dc"])) < 1e-2*np.max(np.abs(parts["homodyne"]))
    _, _, parts, poles = dimer(1.0, 20., 0.5).get_esr_spectrum(f, **kw)
    st = np.argmin(np.abs(poles["frequency"] - 1.118))
    assert poles["form_factor"][st] == pytest.approx(0.369, abs=1e-3)
    # the DC part goes as Omega^2 and the homodyne one as Omega, so compare
    # the DC part with itself: it grows by 180 with the mixing
    assert np.max(np.abs(parts["dc"])) > 50*dc0


@pytest.mark.parametrize("Bt", [0.2, 1.0, 3.0])
def test_yang_dimer_lines(Bt):
    """The dimer of Yang et al. eq. S7 with the tip along the field (D = S1z)
    and the transverse drive S1x of their Section 7: the lines I to IV carry
    |<n|S1x|m>|^2 = sin^2(xi/2)/4 (I, II) and cos^2(xi/2)/4 (III, IV), as
    in eq. S14, and the form factors sin^2(xi/2) and cos^2(xi/2), with
    tan(xi) = (J - D)/(gamma hbar B_tip)"""
    J, Dd, B = 1.0, 0.2, 25.2
    sc = dimer(J, B, Bt, D=Dd)
    f = np.linspace(23., 30., 701)
    _, _, poles = sc.get_esr_spectrum(f, tip=[0, 0, 1], drive=sc.Sx[0],
                                      T=1.2, f_ref=25.2, rabi=1e-4,
                                      return_poles=True)
    xi = np.arctan2(J - Dd, Bt)
    s2, c2 = np.sin(xi/2)**2, np.cos(xi/2)**2
    # eq. S10: f_I = E- - E00, f_II = E11 - E+, f_III = E+ - E00, f_IV = E11 - E-
    r = np.sqrt((J - Dd)**2 + Bt**2)
    E11, E00 = 0.25*(J + 2*Dd) + 0.5*(2*B + Bt), 0.25*(J + 2*Dd) - 0.5*(2*B + Bt)
    Ep, Em = -0.25*(J + 2*Dd) + 0.5*r, -0.25*(J + 2*Dd) - 0.5*r
    expect = {Em - E00: s2, E11 - Ep: s2, Ep - E00: c2, E11 - Em: c2}
    driven = poles["weight"] > 1e-12
    assert np.sum(driven) == 4
    for fr, w in expect.items():
        k = np.argmin(np.abs(poles["frequency"] - fr))
        assert poles["frequency"][k] == pytest.approx(fr, abs=1e-9)
        assert poles["weight"][k] == pytest.approx(w/4, rel=1e-9)
        assert poles["form_factor"][k] == pytest.approx(w, rel=1e-9)


def test_zero_field_raises():
    """An SU(2) multiplet is degenerate and coupled by S: the populations
    inside it depend on the eigensolver's basis"""
    with pytest.raises(ValueError, match="degenerate"):
        dimer(1.0, 0.0, 0.0).get_esr_spectrum(np.linspace(0.5, 1.5, 11),
                                              tip=[0.5, 0, 0.8], f_ref=1.)


def chain(n, iv="python"):
    """An AF chain, J around 30 GHz, in a 20 GHz field, with a tip field of
    2 GHz along z and 0.5 GHz along x on site 0 (no conserved quantity)"""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=iv)
    h = 0
    for i in range(n-1):
        J = 30.*(1 + 0.2*np.cos(i))
        h = h + H*J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    h = h + H*20.*sum(sc.Sz) + H*2.0*sc.Sz[0] + H*0.5*sc.Sx[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 64, 20, 1e-14
    return sc


KW = dict(tip=[0.5, 0., 0.866], site=0, T=0.05, T1=100., T2=2., f_ref=20.,
          rabi=1e-4, vrf_over_vdc=0.1, p_cut=1e-4)
F = np.linspace(15., 30., 751)


def test_sparse_ed_route_matches_dense(monkeypatch):
    """Above algebra.maxsize the lowest eigenpairs come from eigsh"""
    _, a = chain(8).get_esr_spectrum(F, **KW)
    monkeypatch.setattr(spectrum.algebra, "maxsize", 100)
    _, b = chain(8).get_esr_spectrum(F, **KW)
    assert np.max(np.abs(a - b)) < 1e-8*np.max(np.abs(a))


BACKENDS = ["python"] + ([3] if cppext.available(3) else [])


@pytest.mark.parametrize("iv", BACKENDS)
def test_mps_route_matches_ed(iv):
    """The 24 lowest DMRG states of a 6-site chain against ED"""
    _, a = chain(6, iv).get_esr_spectrum(F, mode="ED", **KW)
    _, b, poles = chain(6, iv).get_esr_spectrum(F, mode="DMRG", nstates=24,
                                                return_poles=True, **KW)
    assert np.max(np.abs(a - b)) < 1e-6*np.max(np.abs(a))
    assert "position_bound" in poles


def test_mps_route_raises_on_an_incomplete_manifold():
    with pytest.raises(ValueError, match="raise nstates"):
        chain(6).get_esr_spectrum(F, mode="DMRG", nstates=4, **KW)


@pytest.mark.parametrize("cls", [fermionchain.Spinful_Fermionic_Chain,
                                 fermionchain.Spinful_Fermionic_Chain_Native])
def test_spinful_fermion_atom_is_a_spin_half(cls):
    """One Hubbard orbital at half filling (U = 2 meV, mu = U/2) in a field:
    S does not move charge, so the empty and doubly occupied states are
    blocks of their own in the relaxation, and their thermal weight stays
    there. The spectrum is that of a free spin 1/2 times the weight w1 of
    the singly occupied block (1 - 5.8e-5 at 1.2 K)"""
    fc = cls(1)
    U = 2e-3
    h = U*fc.Nup[0]*fc.Ndn[0] - 0.5*U*(fc.Nup[0] + fc.Ndn[0]) + H*20.*fc.Sz[0]
    fc.set_hamiltonian(h)
    th, f = 0.7, np.linspace(19.95, 20.05, 21)
    _, dI = fc.get_esr_spectrum(f, tip=[np.sin(th), 0, np.cos(th)], T=1.2,
                                f_ref=20., rabi=1e-4)
    x, u = H*20./(2*KB*1.2), U/(2*KB*1.2)
    w1 = 2*np.cosh(x)/(2*np.cosh(x) + 2*np.exp(-u))
    ref = w1*bloch(f, 20.0, th, 1.2, 100., 50., 1e-4)
    assert np.max(np.abs(dI - ref)) < 1e-8*np.max(np.abs(ref))
    # the Floquet reference keeps each block at its weight too
    _, fl = fc.get_esr_spectrum(f, tip=[np.sin(th), 0, np.cos(th)], T=1.2,
                                f_ref=20., rabi=1e-4, method="floquet")
    assert np.max(np.abs(fl - ref)) < 1e-3*np.max(np.abs(ref))


def test_bad_input_raises():
    sc = free_spin()
    f = np.linspace(19.9, 20.1, 5)
    with pytest.raises(ValueError, match="tip"):
        sc.get_esr_spectrum(f)
    with pytest.raises(ValueError, match="T2 cannot exceed"):
        sc.get_esr_spectrum(f, tip=[1, 0, 1], T1=10., T2=30., f_ref=20.)
    with pytest.raises(ValueError, match="f_ref"):
        sc.get_esr_spectrum(f, tip=[1, 0, 1])
    with pytest.raises(ValueError, match="mode"):
        sc.get_esr_spectrum(f, tip=[1, 0, 1], mode="ed")
    with pytest.raises(ValueError, match="floquet"):
        sc.get_esr_spectrum(f, tip=[1, 0, 1], mode="DMRG", method="floquet")
    with pytest.raises(ValueError, match="positive"):
        sc.get_esr_spectrum(-f, tip=[1, 0, 1])
    with pytest.raises(TypeError):
        sc.get_esr_spectrum(f, tip=[1, 0, 1], T_1=10.)


def test_f_ref_pins_the_rates():
    """T1 is defined at f_ref, so zooming onto the strongest line does not
    move it, while another f_ref does. The two windows end at the same
    frequency, so that they keep the same states (the kept states end
    WINDOW_MARGIN widths above the window, and the tails of the lines
    further out are dropped, 2.3e-4 of this line)"""
    sc = chain(6)
    kw = dict(KW)
    fa, a = sc.get_esr_spectrum(F, **kw)
    k = np.argmax(np.abs(a))
    fb = np.array([fa[k] - 0.01, fa[k], fa[k] + 0.01, F[-1]])
    _, b = sc.get_esr_spectrum(fb, **kw)
    assert b[1] == pytest.approx(a[k], rel=1e-8)
    _, c = sc.get_esr_spectrum(fb, **dict(kw, f_ref=22.5)) # another f_ref
    assert abs(c[1]/a[k] - 1.) > 0.01


def test_saturated_solve_against_the_dense_rate_equation():
    """Lines.solve (pairs above SAT_TOL exactly, the others at first order)
    against the null vector of the dense R + W(w) at strong drive"""
    sc = chain(6)
    V = 0.5*sc.Sx[0] + 0.866*sc.Sz[0]
    emax = H*(30. + 30/(np.pi*2.)) + 2*KB*0.05*np.log(1e4)
    spec = spectrum.ed_spectrum(sc, V, V, emax)
    # rabi = 1 GHz: the Korringa T1 of these states is short, and a test of
    # the linear algebra needs saturated pairs next to first-order ones
    lines = lineshape.Lines(spec, F, H, KB*0.05, 100., 2., 1.0, 1e-4, 20.)
    n = spec.dim
    for fi in (18.0, 20.76, 21.83, 24.0):
        w = 2*np.pi*fi
        y, wp = lines.solve(w)
        s = -wp*lines.Kdiag
        W = np.zeros((n, n))
        for p in range(lines.npairs):
            a = np.zeros(n)
            a[lines.u[p]], a[lines.l[p]] = 1., -1.
            W -= wp[p]*np.outer(a, a)
        A = np.vstack([lines.R + W, np.ones((1, n))])
        b = np.zeros(n + 1)
        b[-1] = 1.
        pd = np.linalg.lstsq(A, b, rcond=None)[0]
        dc_dense = np.sum((pd - lines.p0)*lines.Ddiag)
        dc = np.sum(wp*y*lines.d)
        assert np.max(s) > 1. # some pairs saturate
        assert np.min(s) < lineshape.SAT_TOL # some take the first-order branch
        assert dc == pytest.approx(dc_dense, rel=1e-8, abs=1e-14)
        # the population differences, which the homodyne part reads, carry
        # the first-order error of the far pairs, about their number times
        # SAT_TOL (5e-8 measured); the DC part is second order in it
        assert y == pytest.approx(pd[lines.u] - pd[lines.l], abs=1e-6)


def test_drive_across_charge_blocks_raises():
    """A drive that moves charge on a Hubbard atom whose charge excitations
    fall in the window: the relaxation (spin operators) does not connect
    the charge sectors, so the steady state is not defined"""
    fc = fermionchain.Spinful_Fermionic_Chain(1)
    U = 2*H*20.
    fc.set_hamiltonian(U*fc.Nup[0]*fc.Ndn[0] - 0.5*U*(fc.Nup[0] + fc.Ndn[0])
                       + H*10.*fc.Sz[0])
    with pytest.raises(ValueError, match="does not connect"):
        fc.get_esr_spectrum(np.linspace(5., 35., 31), tip=[0.6, 0, 0.8],
                            drive=fc.Cdagup[0] + fc.Cup[0], T=1.2, f_ref=10.)


def test_relaxation_with_two_stationary_states_raises():
    """At T=0 two states with no lower state and one above feeding both"""
    from dmrgpy.esrtk import rates
    G = np.zeros((3, 3))
    G[0, 2] = G[1, 2] = 1.
    R = rates.generator(G)
    p0 = np.array([1., 0., 0.])
    with pytest.raises(ValueError, match="stationary"):
        rates.resolvent(R, p0, np.zeros(3, dtype=int), np.array([0., 0.5, 1.]))


def test_non_hermitian_drive_raises():
    sc = free_spin()
    with pytest.raises(ValueError, match="Hermitian"):
        sc.get_esr_spectrum(np.linspace(19.9, 20.1, 5), tip=[0.6, 0, 0.8],
                            f_ref=20., drive=sc.Sx[0] + 1j*sc.Sy[0])


def test_zero_field_raises_on_dmrg():
    """The members of a multiplet from DMRG are split by more than rounding;
    the tolerance there includes their energy uncertainty"""
    sz = spinchain.Spin_Chain(["S=1/2"]*4, itensor_version="python")
    sz.set_hamiltonian(H*30.*sum(sz.Sx[i]*sz.Sx[i+1] + sz.Sy[i]*sz.Sy[i+1]
                                 + sz.Sz[i]*sz.Sz[i+1] for i in range(3)))
    sz.maxm, sz.nsweeps = 32, 20
    with pytest.raises(ValueError, match="degenerate"):
        sz.get_esr_spectrum(np.linspace(10., 25., 31), tip=[0.6, 0, 0.8],
                            T=0.05, T2=2., f_ref=20., p_cut=1e-4, mode="DMRG",
                            nstates=16)


def test_homodyne_taper():
    """1 up to HOMODYNE_INNER free-spin widths 2/T2 (rad/ns) from the
    resonance, 0 beyond WINDOW_MARGIN, even and monotonic in between"""
    T2 = 4.
    w = 2./T2
    x = np.linspace(0., 40.*w, 4001)
    t = lineshape.homodyne_taper(x, T2)
    assert np.all(t[x <= lineshape.HOMODYNE_INNER*w] == 1.)
    assert np.all(t[x >= lineshape.WINDOW_MARGIN*w] == 0.)
    assert np.all(np.diff(t) <= 0.)
    assert np.array_equal(lineshape.homodyne_taper(-x, T2), t)


def test_homodyne_part_does_not_depend_on_the_window():
    """The dispersive tail of a line falls only as 1/detuning, so each
    coherence is cut between 15 and 30 free-spin widths from its resonance;
    without the cut, zooming the window onto the strongest line moved the
    homodyne part by 1e-2 of the DC line on this chain, since the kept
    states change with the window (measured 1.2e-11 with it). The DC part
    moves by the DC tails of the lines past the margin (1.8e-5)"""
    sc = chain(6)
    kw = dict(KW, p_cut=1e-8, return_parts=True)
    _, _, a = sc.get_esr_spectrum(F, **kw)
    k = np.argmax(np.abs(a["dc"]))
    fb = np.array([F[k] - 0.02, F[k], F[k] + 0.02])
    _, _, b = sc.get_esr_spectrum(fb, **kw)
    pk = np.max(np.abs(a["dc"]))
    assert abs(b["homodyne"][1] - a["homodyne"][k]) < 1e-8*pk
    assert abs(b["dc"][1] - a["dc"][k]) < 1e-4*pk


def test_krylov_route_matches_ed():
    """submode="KRYLOV": the populated eigenstates and one block Krylov
    space from V|m>, against the eigenbasis; with too few blocks the final
    states of the window lines are not converged and the call raises"""
    _, a = chain(8).get_esr_spectrum(F, **KW)
    _, b = chain(8).get_esr_spectrum(F, submode="KRYLOV", nkry=60, **KW)
    assert np.max(np.abs(a - b)) < 1e-8*np.max(np.abs(a))
    with pytest.raises(ValueError, match="raise nkry"):
        chain(8).get_esr_spectrum(F, submode="KRYLOV", nkry=3, **KW)


def test_krylov_route_inputs():
    sc = chain(4)
    with pytest.raises(ValueError, match="submode"):
        sc.get_esr_spectrum(F, submode="TASK", **KW)
    with pytest.raises(ValueError, match="nkry"):
        sc.get_esr_spectrum(F, submode="KRYLOV", nkry=0, **KW)


@pytest.mark.parametrize("iv", BACKENDS)
def test_dmrg_krylov_route_matches_ed(iv):
    """submode="KRYLOV" on DMRG: the DMRG states below the window and the
    final states of the window lines from the tangent space of each
    populated state (TASK), in one common basis, against ED (measured
    1.1e-5 of the peak, rising toward the top of the window: the DC tail of
    the lines past the 30-width margin, which ED keeps among its states and
    this route does not); it raises when nkry is too small or the DMRG
    states do not reach below the window"""
    _, a = chain(8).get_esr_spectrum(F, **KW)
    sc = chain(8, iv)
    _, b = sc.get_esr_spectrum(F, mode="DMRG", submode="KRYLOV", nstates=4,
                               **KW)
    assert np.max(np.abs(a - b)) < 5e-5*np.max(np.abs(a))
    with pytest.raises(ValueError, match="raise nkry"):
        chain(8, iv).get_esr_spectrum(F, mode="DMRG", submode="KRYLOV",
                                      nstates=4, nkry=3, **KW)
    with pytest.raises(ValueError, match="raise nstates"):
        chain(8, iv).get_esr_spectrum(F, mode="DMRG", submode="KRYLOV",
                                      nstates=1, **KW)


@pytest.mark.parametrize("iv", BACKENDS)
def test_pair_traces_match_amb(iv):
    """The DMRG route reads every one-site operator from one sweep per pair
    of states (spectrum._PairTraces), with the local matrix read off the
    MPO by a partial trace (spectrum.local_operator); against aMb. An
    operator on two sites is not read that way."""
    sc = chain(6, iv)
    sc.maxm = 16
    es, ws = sc.get_excited_states(n=4)
    ws = [w*(1./np.sqrt(abs(sc.overlap(w, w)))) for w in ws]
    pt = spectrum._PairTraces(sc, ws)
    Sz, Sp = spectrum.site_operators(sc)
    for op in [Sz[2], Sp[2], Sp[5], 0.5*sc.Sx[0] + 0.866*sc.Sz[0], sc.Sy[3]]:
        M = pt.matrix(*spectrum.local_operator(sc, op))
        R = np.array([[sc.aMb(a, op, b) for b in ws] for a in ws])
        assert np.max(np.abs(M - R)) < 1e-12
    assert spectrum.local_operator(sc, sc.Sz[0]*sc.Sz[1]) is None
