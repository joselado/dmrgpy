"""Second-order (golden rule) spin-excitation dI/dV with superconducting
electrodes: Ternes' second order with a BCS density of states on the
sample, on the tip, or on both.

Ternes' second order, kondospectrumtk/conductance.py::second_order_dIdV,
is a golden rule between the eigenstates |i>, |f> of the impurity spin
alone, with a flat density of states in tip and sample, so every
transition eps_fi = e_f - e_i gives a temperature-broadened step
Theta(eV - eps_fi). With a gap in either electrode the electron has to
find a state there, and the net current is

    I(V) = 2 pi T0^2 sum_{i,f} p_i |M_fi|^2 [C_ts(eV - eps_fi) - C_st(eV + eps_fi)],
    C_ts(x) = int dw rho_t(w) f(w) rho_s(w + x) (1 - f(w + x)),
    C_st(y) = int dw rho_t(w) (1 - f(w)) rho_s(w + y) f(w + y),

with w the energy of the electron measured from the Fermi level of the
tip, f the Fermi function at the temperature of the KondoSpectrum, and
|M_fi|^2 = |S^x_fi|^2 + |S^y_fi|^2 + |S^z_fi|^2 the same spin-averaged
matrix element as second_order_dIdV (its _spin_matrix_elements_squared is
reused, not copied). The elastic potential channel is the same with
eps = 0 and 4 U^2 in place of |M|^2. rho_t and rho_s are normalized to 1
in the normal state, here the Dynes form

    rho(w) = Re[(w + i gamma)/sqrt((w + i gamma)^2 - Delta^2)],

which is the BCS density of states convolved with a Lorentzian of
half-width gamma. For rho = 1, C_ts(x) = x/(1 - exp(-x/kT)), C_st(y) is
the same at -y, and dI/dV = Theta(eV - eps) + Theta(-eV - eps), exactly
second_order_dIdV. The dI/dV is then

    dI/dV = 2 pi T0^2 sum p_i |M_fi|^2 [G_ts(eV - eps_fi) - G_st(eV + eps_fi)],

G = dC/dx. At T=0 with a normal tip and a superconducting sample,
G_ts(x) = rho_s(x) for x > 0, so every step becomes the coherence peak
of the sample at eV = Delta_s + eps_fi. With two gaps the convolution of
the two 1/sqrt edges makes C jump at Delta_s + Delta_t, so dI/dV has a
peak at eV = +-(Delta_s + Delta_t + eps_fi) whose height is set by gamma,
and at T > 0 the thermally excited quasiparticles add peaks at
eV = +-|Delta_s - Delta_t| +- eps_fi. This is the experiment of
Heinrich, Braun, Pascual and Franke, arXiv:1405.6984 (Fe-OEP-Cl, S=5/2,
on Pb with a Pb tip).

How C is computed. C = C_normal + (C - C_normal): the normal part is
closed form, and the difference is a cross-correlation on a uniform
energy grid (FFT), with the normal-state integrand on the same grid
subtracted, so the part of the quadrature error that does not depend on
the gap cancels; G of the difference is a central difference on the
grid. When both gaps are zero the difference is skipped and the result
is second_order_dIdV to roundoff. The grid step defaults to a twentieth of
the smallest of gamma, kT and Delta/20, so a gap needs gamma > 0: at gamma = 0
and T = 0 the peak is a delta function no grid resolves.

What was checked (tests/test_ysr_bcs_second_order.py): Delta = 0 gives
second_order_dIdV to roundoff, and the grid route forced at Delta = 0
gives it to the quadrature accuracy; a normal tip on a superconductor at
T=0 gives rho_s(x) theta(x); two gaps at T=0 give the current jump
pi sqrt(Delta_s Delta_t)/2 at Delta_s + Delta_t (to 0.1% to 0.5% at
gamma <= 1e-3 Delta, as a Lorentzian of half-width gamma_s + gamma_t
centred there); a superconducting tip gives the same as
ysr.tip_convolution; the current agrees with an independent golden-rule
sum over the explicit Bogoliubov quasiparticles of a discretized BCS
band in tip and sample (coherence factors u_k^2, v_k^2, quasiparticle
occupations), with a difference linear in gamma (2e-3 of max|I| at
gamma = 5e-3 Delta), since the sum puts the occupation at the level
energy and the Dynes form at the broadened one; and at J=0 the
many-body spectrum of ysr.get_stm_spectrum is this golden rule with the
local density of states of the Wilson chain, up to a factor 4 pi of
normalization.

What it does not contain: it is second order in the exchange with the
substrate, so it has no Yu-Shiba-Rusinov state (a bound state of J to all
orders, see ysr.py) and no Kondo logarithm; the spin occupations p_i are
thermal, so the spin pumping of the experiment (the second peak at
2 Delta + 4D) is not there; and the excited states have no lifetime
width. substrate_relaxation_rates gives the golden-rule rate of the
substrate alone, with the coherence factors of a spin vertex (pair
creation (1 - Delta^2/(E E')), quasiparticle scattering the
Hebel-Slichter (1 + Delta^2/(E E'))), checked against a direct BCS
quadrature at T=0 and T>0. It vanishes at T=0 when |eps_fi| < 2 Delta_s,
the "protection" of an excited spin by the gap, and starts continuously
above it.

A kernel from bcs_kernel can be reused with a new bias list; it raises
when the biases plus the transition energies leave the x range it was
tabulated on, or when the spectrum is at another temperature.
"""

import numpy as np

from ..kondospectrumtk.stepfunctions import Theta, Theta0

# largest grid the FFT correlation is allowed to build
MAX_GRID = 2**23


def dynes_dos(e, Delta, gamma):
    """Dynes density of states, normalized to 1 far from the gap, and
    exactly 1 for Delta = 0 whatever gamma is."""
    e = np.asarray(e, dtype=float)
    if Delta == 0.0:
        return np.ones_like(e)
    if gamma <= 0.0:
        raise ValueError("a gap needs a Dynes broadening gamma > 0")
    from ..ysr import bcs_dos
    return bcs_dos(e, Delta, gamma)


def dynes_pair_dos(e, Delta, gamma):
    """Anomalous (pair) density of states of the Dynes form,
    Re[Delta/sqrt((e + i gamma)^2 - Delta^2)] on the branch where
    dynes_dos is positive, so that it is sign(e) Delta/sqrt(e^2 - Delta^2)
    outside the BCS gap and 0 inside it, and 0 for Delta = 0."""
    e = np.asarray(e, dtype=float)
    if Delta == 0.0:
        return np.zeros_like(e)
    if gamma <= 0.0:
        raise ValueError("a gap needs a Dynes broadening gamma > 0")
    z = e + 1j*gamma
    s = np.sqrt(z*z - Delta*Delta)
    s = np.where(np.real(z/s) < 0.0, -s, s)
    return np.real(Delta/s)


def _fermi(w, kT):
    """Fermi function, f(0) = 1/2 also at T=0."""
    w = np.asarray(w, dtype=float)
    if kT == 0.0:
        return np.where(w < 0, 1.0, np.where(w > 0, 0.0, 0.5))
    return 0.5*(1.0 - np.tanh(w/(2.0*kT)))


def _normal_current(x, kT):
    """int dw f(w)(1 - f(w + x)) = x/(1 - exp(-x/kT)), max(x, 0) at T=0."""
    x = np.asarray(x, dtype=float)
    if kT == 0.0:
        return np.maximum(x, 0.0)
    out = np.full_like(x, kT)
    nz = x != 0.0
    out[nz] = x[nz]/(-np.expm1(-x[nz]/kT))
    return out


def _normal_conductance(x, kT):
    """d/dx of _normal_current: the step Theta of stepfunctions.py."""
    x = np.asarray(x, dtype=float)
    if kT == 0.0:
        return Theta0(x)
    return Theta(x/kT)


class BCSKernel():
    """C_ts, C_st and their derivatives G_ts, G_st of the module
    docstring, tabulated once for |x| <= xmax.

    rho_t, rho_s: callables energy -> density of states (normalized to 1
    in the normal state), or None for a flat band. de is the grid step.
    force_grid computes the gap-dependent part on the grid even for two
    flat bands (for testing the quadrature against the closed form);
    subtract_normal=False tabulates C itself instead of C - C_normal."""
    def __init__(self, xmax, kT, rho_t=None, rho_s=None, de=None,
                 force_grid=False, subtract_normal=True):
        self.kT = float(kT)
        self.xmax = float(xmax)
        self.flat = (rho_t is None and rho_s is None and not force_grid)
        self.subtract_normal = subtract_normal
        if self.flat:
            return
        if de is None or de <= 0.0:
            raise ValueError("BCSKernel needs a grid step de > 0")
        from scipy.signal import fftconvolve
        h = float(de)
        W = xmax + 40.0*self.kT + 20.0*h # the Fermi window of every x
        m = int(np.ceil(W/h))
        if 2*m + 1 > MAX_GRID:
            raise ValueError("the energy grid would have %d points; raise "
                "gamma or de, or lower the bias range" % (2*m + 1))
        w = h*np.arange(-m, m + 1) # w = 0 is a grid point
        f = _fermi(w, self.kT)
        rt = np.ones_like(w) if rho_t is None else np.asarray(rho_t(w), float)
        rs = np.ones_like(w) if rho_s is None else np.asarray(rho_s(w), float)
        kx = int(np.ceil(xmax/h)) + 2
        N = len(w)
        def corr(a, b):
            """h sum_j a_j b_{j+k}, for k = -kx..kx"""
            full = fftconvolve(b, a[::-1], mode="full")
            return h*full[N - 1 - kx:N + kx]
        Cts = corr(rt*f, rs*(1.0 - f))
        Cst = corr(rt*(1.0 - f), rs*f)
        if subtract_normal:
            Cts = Cts - corr(f, 1.0 - f)
            Cst = Cst - corr(1.0 - f, f)
        self.x = h*np.arange(-kx, kx + 1)
        self.Cts, self.Cst = Cts, Cst
        self.Gts = np.gradient(Cts, h)
        self.Gst = np.gradient(Cst, h)
    def _part(self, name, x, normal):
        """normal(x) plus the tabulated gap part. At T=0 the normal part
        has a step (a kink for the current) at x=0, a grid point, so the
        sum is tabulated and interpolated as a whole: adding the exact
        step to an interpolated table would leave a spike one grid step
        wide at x=0."""
        x = np.asarray(x, dtype=float)
        if self.flat:
            return normal(x)
        table = getattr(self, name)
        if not self.subtract_normal:
            return np.interp(x, self.x, table)
        if self.kT == 0.0:
            return np.interp(x, self.x, normal(self.x) + table)
        return normal(x) + np.interp(x, self.x, table)
    def current_ts(self, x):
        return self._part("Cts", x, lambda z: _normal_current(z, self.kT))
    def current_st(self, y):
        return self._part("Cst", y, lambda z: _normal_current(-z, self.kT))
    def conductance_ts(self, x):
        return self._part("Gts", x,
                          lambda z: _normal_conductance(z, self.kT))
    def conductance_st(self, y):
        return self._part("Gst", y,
                          lambda z: -_normal_conductance(-z, self.kT))


def _default_step(kT, Delta_s, Delta_t, gamma_s, gamma_t):
    """A twentieth of the sharpest energy scale: gamma on a gapped side, kT,
    and the gaps themselves."""
    scales = [kT] if kT > 0.0 else []
    for D, g in ((Delta_s, gamma_s), (Delta_t, gamma_t)):
        if D != 0.0:
            scales += [g, abs(D)/20.0]
    if not scales: # no scale at all: the caller has to pass de
        return 0.0
    return 0.05*min(scales)


def _dos_callables(Delta_s, Delta_t, gamma_s, gamma_t, rho_s, rho_t):
    for D, g, name in ((Delta_s, gamma_s, "sample"), (Delta_t, gamma_t, "tip")):
        if D != 0.0 and g <= 0.0:
            raise ValueError("a gap on the %s side needs a Dynes "
                "broadening gamma > 0 (at gamma = 0 and T = 0 the peaks are "
                "delta functions no grid resolves)" % name)
    if rho_s is None and Delta_s != 0.0:
        rho_s = lambda w: dynes_dos(w, Delta_s, gamma_s)
    if rho_t is None and Delta_t != 0.0:
        rho_t = lambda w: dynes_dos(w, Delta_t, gamma_t)
    return rho_s, rho_t


def _weights(ks):
    """w[i,f] = p_i |M_fi|^2 and eps[i,f] = e_f - e_i over the occupied
    initial states i (conductance.P_CUT, as in the third-order terms)."""
    from ..kondospectrumtk import conductance
    occ = conductance._occupied_states(ks)
    M2 = conductance._spin_matrix_elements_squared(ks) # [f,i]
    w = ks.p[occ, None]*M2[:, occ].T
    eps = ks.e[None, :] - ks.e[occ, None]
    return w, eps


def bcs_kernel(ks, eVs, Delta_s=0.0, Delta_t=0.0, gamma_s=0.0,
               gamma_t=0.0, rho_s=None, rho_t=None, de=None, **kwargs):
    """The BCSKernel second_order_dIdV_bcs uses for ks at the biases eVs
    (kwargs go to BCSKernel)."""
    eVs = np.asarray(eVs, dtype=float)
    kT = ks.kB*ks.T
    rho_s, rho_t = _dos_callables(Delta_s, Delta_t, gamma_s, gamma_t,
                                  rho_s, rho_t)
    w, eps = _weights(ks)
    emax = np.max(np.abs(eps[w > 0])) if np.any(w > 0) else 0.0
    xmax = np.max(np.abs(eVs)) + emax
    if de is None and (rho_s is not None or rho_t is not None
                       or kwargs.get("force_grid", False)):
        de = _default_step(kT, Delta_s, Delta_t, gamma_s, gamma_t)
        if de <= 0.0:
            raise ValueError("pass a grid step de: there is no gap, gamma "
                             "or temperature to set one")
    return BCSKernel(xmax, kT, rho_t=rho_t, rho_s=rho_s, de=de, **kwargs)


def _check_kernel(K, ks, eVs, w, eps):
    """A kernel is reused only at its own temperature and inside the x
    range it was tabulated on (np.interp would clamp silently)."""
    kT = ks.kB*ks.T
    if abs(K.kT - kT) > 1e-12*max(K.kT, kT):
        raise ValueError("the kernel was built at kT = %g, the spectrum is "
                         "at kT = %g" % (K.kT, kT))
    if K.flat:
        return
    emax = np.max(np.abs(eps[w > 0])) if np.any(w > 0) else 0.0
    need = (np.max(np.abs(eVs)) if eVs.size else 0.0) + emax
    if need > K.xmax*(1.0 + 1e-12):
        raise ValueError("the kernel covers |x| <= %g, these biases and "
                         "transitions need %g: build it with bcs_kernel for "
                         "these biases" % (K.xmax, need))


def _assemble(ks, eVs, K, T0, U, which):
    eVs = np.asarray(eVs, dtype=float)
    w, eps = _weights(ks)
    _check_kernel(K, ks, eVs, w, eps)
    ts = K.current_ts if which == "I" else K.conductance_ts
    st = K.current_st if which == "I" else K.conductance_st
    xa = eVs[:, None, None] - eps[None, :, :]
    xb = eVs[:, None, None] + eps[None, :, :]
    spin = np.einsum("if,eif->e", w, ts(xa) - st(xb))
    elastic = 4.0*U**2*(ts(eVs) - st(eVs))
    return 2*np.pi*T0**2*(spin + elastic)


def second_order_dIdV_bcs(ks, eVs, Delta_s=0.0, Delta_t=0.0, gamma_s=None,
                          gamma_t=None, gamma=0.0, T0=1.0, U=0.0,
                          rho_s=None, rho_t=None, de=None, kernel=None):
    """Second-order dI/dV(eV) of the spin spectrum ks (a KondoSpectrum,
    edkondo.py, or anything with e, p, Sx, Sy, Sz, T and kB) with Dynes
    superconductors on the sample (Delta_s, gamma_s) and on the tip
    (Delta_t, gamma_t), at the temperature ks.T; gamma is the default of
    both gammas. Units and prefactors are those of
    conductance.second_order_dIdV, which this is at Delta_s = Delta_t = 0.

    rho_s, rho_t: callables energy -> density of states (1 in the normal
    state) that replace the Dynes form on that side. de: grid step of
    the energy integral (default a twentieth of the smallest of gamma, kT
    and Delta/20).
    kernel: a BCSKernel from bcs_kernel to reuse (then the gap keywords
    are ignored)."""
    gamma_s = gamma if gamma_s is None else gamma_s
    gamma_t = gamma if gamma_t is None else gamma_t
    if kernel is None:
        kernel = bcs_kernel(ks, eVs, Delta_s, Delta_t, gamma_s, gamma_t,
                            rho_s=rho_s, rho_t=rho_t, de=de)
    return _assemble(ks, eVs, kernel, T0, U, "dIdV")


def second_order_current_bcs(ks, eVs, Delta_s=0.0, Delta_t=0.0,
                             gamma_s=None, gamma_t=None, gamma=0.0, T0=1.0,
                             U=0.0, rho_s=None, rho_t=None, de=None,
                             kernel=None):
    """The current I(eV) whose derivative second_order_dIdV_bcs is, in the
    same units times energy (same keywords)."""
    gamma_s = gamma if gamma_s is None else gamma_s
    gamma_t = gamma if gamma_t is None else gamma_t
    if kernel is None:
        kernel = bcs_kernel(ks, eVs, Delta_s, Delta_t, gamma_s, gamma_t,
                            rho_s=rho_s, rho_t=rho_t, de=de)
    return _assemble(ks, eVs, kernel, T0, U, "I")


def get_spin_excitation_spectrum(chain, eVs, site=0, T=1.0, Delta_s=0.0,
                                 Delta_t=0.0, gamma=None, gamma_s=None,
                                 gamma_t=None, T0=1.0, U=0.0,
                                 kB=8.617333262e-5, de=None):
    """dI/dV(eV) of the spin on site of chain (a Spin_Chain, or a spinful
    fermion chain) at temperature T (Kelvin with the default kB in eV/K)
    on a superconductor of gap Delta_s with a tip of gap Delta_t, from the
    full ED spectrum of the chain (edkondo.KondoSpectrum), as
    second_order_dIdV_bcs. gamma is the Dynes broadening of both sides
    (default a hundredth of the larger gap)."""
    from ..kondospectrumtk.edkondo import KondoSpectrum
    if gamma is None:
        gamma = 1e-2*max(abs(Delta_s), abs(Delta_t))
    ks = KondoSpectrum(chain, site, T, kB=kB)
    return eVs, second_order_dIdV_bcs(ks, eVs, Delta_s=Delta_s,
                                      Delta_t=Delta_t, gamma=gamma,
                                      gamma_s=gamma_s, gamma_t=gamma_t,
                                      T0=T0, U=U, de=de)


def substrate_relaxation_rates(ks, Jrho_s, Delta_s=0.0, gamma_s=0.0,
                               rho_s=None, pair_s=None, de=None):
    """Golden-rule rates R[f,i] of the transition i -> f of the spin by
    scattering one substrate electron, for J S.s_0 with s_0 the spin
    density psi^dag sigma psi/2 of the substrate under the spin and
    Jrho_s = J times the normal density of states per spin:

        R[f,i] = 2 pi (J rho)^2 (|M_fi|^2/2) int dw f(w) (1 - f(w - eps_fi))
                 [rho_s(w) rho_s(w - eps_fi) + F_s(w) F_s(w - eps_fi)],

    in units of energy/hbar (the factor 1/2 is the trace over the
    electron spin, sum_{s s'} |<s'|sigma/2|s>.S_fi|^2 = |M_fi|^2/2).
    F_s is the anomalous density of states, sign(w) Delta/sqrt(w^2 -
    Delta^2) for BCS (dynes_pair_dos). The F F term is the coherence factor
    of a time-reversal odd vertex (the spin density): with w and w - eps of
    opposite sign (creating or destroying two quasiparticles) the bracket
    is rho rho' (1 - Delta^2/(E E')), so the rate starts continuously at
    |eps| = 2 Delta; with the same sign (scattering a thermal
    quasiparticle) it is rho rho' (1 + Delta^2/(E E')), the Hebel-Slichter
    factor. In a normal metal at T=0 the decay rate of an excitation eps is
    pi (J rho)^2 |M|^2 eps; on a superconductor at T=0 it vanishes for
    eps < 2 Delta_s, up to the Dynes in-gap states, a fraction linear in
    gamma. Checked against an explicit Bogoliubov sum over a discretized
    band (tests/test_ysr_bcs_second_order.py).

    rho_s: a callable density of states (1 in the normal state) in place of
    the Dynes form; with Delta_s != 0 it needs its pair_s, the anomalous
    density of states, as well (with Delta_s = 0 there is no pairing and
    pair_s defaults to zero)."""
    from ..kondospectrumtk import conductance
    kT = ks.kB*ks.T
    if rho_s is None:
        if Delta_s != 0.0:
            if gamma_s <= 0.0:
                raise ValueError("a gap needs a Dynes broadening gamma_s > 0")
            rho_s = lambda w: dynes_dos(w, Delta_s, gamma_s)
            pair_s = lambda w: dynes_pair_dos(w, Delta_s, gamma_s)
    elif Delta_s != 0.0 and pair_s is None:
        raise ValueError("a callable rho_s on a superconductor needs its "
                         "anomalous density of states pair_s too")
    M2 = conductance._spin_matrix_elements_squared(ks) # [f,i]
    eps = ks.e[:, None] - ks.e[None, :] # [f,i]
    xmax = np.max(np.abs(eps))
    if rho_s is not None and de is None:
        de = _default_step(kT, Delta_s, Delta_s, gamma_s, gamma_s)
        if de <= 0.0:
            raise ValueError("pass a grid step de: there is no gap, gamma "
                             "or temperature to set one")
    K = BCSKernel(xmax, kT, rho_t=rho_s, rho_s=rho_s, de=de)
    # an electron at w (both energies measured from the same Fermi level)
    # leaves at w - eps: C_ts at x = -eps
    C = K.current_ts(-eps)
    if pair_s is not None:
        KF = BCSKernel(xmax, kT, rho_t=pair_s, rho_s=pair_s, de=de,
                       subtract_normal=False)
        C = C + KF.current_ts(-eps)
    rate = 2*np.pi*Jrho_s**2*0.5*M2*C
    np.fill_diagonal(rate, 0.0)
    return rate
