"""Current through one Yu-Shiba-Rusinov (YSR) state with a superconducting
tip: single-electron tunneling with quasiparticle relaxation, and resonant
Andreev reflection, after Ruby, Pientka, Peng, von Oppen, Heinrich and
Franke, arXiv:1502.05048 (PRL 115, 087001), eqs. 1 and 2 (S51, S52).

The YSR state is one fermionic level at energy eps > 0, with electron
weight u^2 and hole weight v^2 at the tip position: the local density of
states of the substrate has a pole u^2 delta(w - eps) in the electron part
and v^2 delta(w + eps) in the hole part. The tip is a BCS superconductor
with gap Delta_t and Dynes broadening gamma, rho(w) = ysr.bcs_dos(w,
Delta_t, gamma) in units of its normal density nu_t, at temperature T. The
level is filled from the tip by electrons at the rate
Gamma_e(w) = g u^2 rho(w - eV) and by holes at the rate
Gamma_h(w) = g v^2 rho(w + eV), with g = 2 pi t^2 nu_t the tip coupling
(gamma_tip), and it relaxes to the continuum: it is emptied at the rate
Gamma_1 and filled by a thermal quasiparticle at the rate Gamma_2. With
Gamma(w) = Gamma_e + Gamma_h + Gamma_1 + Gamma_2 and the tip Fermi function
n(w) at w- = w - eV and w+ = w + eV,

    I_s = (e/h) int dw { Gamma_1 [Gamma_e n(w-) - Gamma_h n(w+)]
            - Gamma_2 [Gamma_e (1 - n(w-)) - Gamma_h (1 - n(w+))] }
            / [(w - eps)^2 + Gamma(w)^2/4],
    I_a = (2e/h) int dw Gamma_e Gamma_h [n(w-) - n(w+)]
            / [(w - eps)^2 + Gamma(w)^2/4].

I_s is a tip electron that fills the level and a relaxation that empties
it (or the reverse), so it needs Gamma_1 or Gamma_2; I_a is an electron
in and a hole out, a Cooper pair into the condensate, and needs no
relaxation. The main peaks are at eV = +-(Delta_t + eps), alpha_+ ~ u^2 and
alpha_- ~ v^2, and the thermal peaks at eV = +-(Delta_t - eps), beta_- ~
Gamma_2 u^2 and beta_+ ~ Gamma_2 v^2. With thermal relaxation, detailed
balance gives Gamma_1/Gamma_2 = exp(eps/T) (S53) for one level, the
default when Gamma2 is None (detailed_balance has the degeneracy factor of
a many-body multiplet). Units: energies in any unit, the current returned as
I/G0 with G0 = 2e^2/h (so in units of energy/e), and dI/dV in units of G0.
The normal-state conductance of the junction is GN = 2 pi nu_s g at weak
coupling, nu_s the substrate density of states per spin, the only regime
where the model holds; gamma_tip_from_GN inverts by default the paper's
GN = 4x/(1+x), x = pi^2 nu^2 t^2 (S89 as printed, which is not bounded by
one channel), and offers the weak-coupling form and the one-channel
Landauer form 4x/(1+x)^2.

As in the paper, the real part of the tip self-energy (a shift of eps of
order t^2), multiple Andreev reflection (Andreev reflection inside the
tip) and the frequency dependence of Gamma_1, Gamma_2 are left out; the
last is why a Lorentzian tail far from eps carries a small current even at
V = 0 (of order Gamma_1 g/W with W the integration window), which does not
change dI/dV.

Weights from a many-body chain (rate_parameters). The poles of
zaverage.stm_poles give the electron weight W_e at +E and the hole weight
W_h at -E of the YSR transition, E = |E_ysr|. For the bare electron at
the site under the tip, rate_parameters returns u^2 = A_Lambda W_e and
v^2 = A_Lambda W_h: on the BdG chain of a classical spin (Lambda = 2 on 40
sites, Lambda = 4 on 24 sites at four z) A_Lambda times the chain weight is
the continuum residue of S32 to 0.1% to 1%, the same A_Lambda as for the
coupling. For the composite tunneling vertex X of
ysr.tunneling_operators the weights are those of the vertex, and they are
divided by the mean total weight per tip spin, sum_sigma
<{X_sigma, X_sigma^dag}>/2 (1 for the bare electron), so that u^2 and v^2
are per unit of the normal-state conductance through the same vertex
(this assumes the weight of X above the gap is spread over the band like
the bare one, which holds up to corrections of order rho J). The ratio
u^2/v^2 needs no normalization, and it is what sets the asymmetry
alpha_+/alpha_- = u^2/v^2 of the linear regime. The rate model is that of
one level: for a degenerate excited (or ground) multiplet, as the doublet
of a quantum spin 1/2, the return rates carry a factor of the degeneracy
that it does not have. The linear regime, where only the rate into the
excited multiplet counts, is not affected; the crossover conductances,
which depend on the rate back, are right only up to a factor of order 2.

Quadrature. The integrand has a Lorentzian of width Gamma at eps and
inverse square-root edges at w = +-eV +- Delta_t, so it is integrated on a
mesh graded geometrically towards each of these points (Gauss-Legendre
panels halving down to hmin), on the window [-W, W]. dI/dV is a central
difference of I with a step h below the smallest width.

What was checked (tests/test_ysr_current.py, energies in units of the tip
gap): the mesh agrees with scipy.integrate.quad to 1e-14. In the linear
regime (g = 1e-8, Gamma_1 = 4e-3, Dynes 1e-2) dI_s/dV is pi g times the
dI/dV of ysr.tip_convolution of u^2 L(w - eps) + v^2 L(w + eps), L the
Lorentzian of half-width Gamma_1/2, to 6e-5 of the peak (2.6e-3 at g = 1e-6,
the next order in g). The threshold values of the supplement are
reproduced: S68 to 1e-5, S71 to 1e-4 and the saturated current
e Gamma_1/(3 hbar) of S70 to 1e-3. At weak coupling alpha_+/alpha_- =
u^2/v^2 to 0.4%, alpha_+ beta_+/(alpha_- beta_-) = 0.9999 and
beta_-/alpha_+ = Gamma_2/Gamma_1 to 1e-4; at T = 0.1, eps = 0.3 with
detailed balance beta_-/alpha_+ = exp(-eps/T) to 8e-4 (at T = 0.15 the
thermal quasiparticles of the tip lower it by exp(-(Delta_t-eps)/T), 1%).
In the Andreev regime alpha_+/alpha_- tends to the paper's (v/u)^(10/3)
only where there is a window between Gamma_1 and the hole broadening
sqrt(omega_h^3/eps) of regime (iii) (S65, S85): with u^2/v^2 = 0.13 and
Gamma_1 = 1e-7 it is 23.5 at g = 1e-5, against 29.2 from the threshold
form S57 (the maximum lies off the threshold, S III F) and 30.0 for
(v/u)^(10/3), and it falls to 12.8 at g = 1e-3 as the hole broadening
reaches omega_e; with Gamma_1 = 1e-3 there is no such window and the ratio
reaches only 5.0 (at g = 1e-2, Gamma_2 = 0). A normal tip gives the
Breit-Wigner Andreev peak 4 u^2 v^2/(u^2+v^2)^2 in 2e^2/h to 2e-5 at both
biases. A positive bias moves tip electrons into the sample, as in
ysr.tip_convolution, so u^2 is the peak at positive bias. On a 5-site kondo_chain
(Lambda = 4, J = 0.8 A_Lambda, Delta = 0.1) the TASK subgap pole is
|E_ysr| of ED to 1e-8, the bare electron gives u^2 = v^2 and total weight
1 per spin, and the composite vertex with U = 0.3 gives u^2/v^2 = 1.557.
That chain does not reach the gap; on 8 sites, the first length that does,
eps = 0.4582 Delta (ED and the TASK pole agree to 1e-7 of Delta) and
u^2/v^2 = 1.554, against 0.4808 Delta and 1.557 on 5 sites.
"""

import numpy as np


def fermi(x, T):
    """Fermi function 1/(exp(x/T)+1), a step (1/2 at x = 0) at T = 0."""
    x = np.asarray(x, dtype=float)
    if T == 0.0:
        return np.where(x < 0.0, 1.0, np.where(x > 0.0, 0.0, 0.5))
    from scipy.special import expit
    return expit(-x/T)


def tip_dos(w, Delta_t, gamma=0.0):
    """Tip density of states in units of its normal value: ysr.bcs_dos
    for a Dynes gamma > 0, the BCS |w|/sqrt(w^2 - Delta_t^2) for
    |w| > Delta_t and exactly 0 inside the gap for gamma = 0, and 1 for a
    normal tip, Delta_t = 0. (ysr.bcs_dos at gamma = 0 keeps an imaginary
    part 1e-14, which leaves a density of order 1e-14/(2 Delta d)^(3/2) at
    a distance d inside the gap edge: harmless in a convolution, but next
    to a level of width Gamma_1 = 1e-7 it added 23% to the saturated
    current.)"""
    w = np.asarray(w, dtype=float)
    if Delta_t == 0.0:
        return np.ones(np.shape(w))
    if gamma > 0.0:
        from ..ysr import bcs_dos
        return bcs_dos(w, Delta_t, gamma)
    d2 = w*w - Delta_t*Delta_t
    return np.where(d2 > 0.0, np.abs(w)/np.sqrt(np.maximum(d2, 1e-300)), 0.0)


def detailed_balance(Gamma1, eps, T, g_ratio=1.0):
    """Gamma_2 = g_ratio Gamma_1 exp(-eps/T), 0 at T = 0. With g_ratio = 1
    this is S53, the occupation of one level. For a many-body YSR
    excitation g_ratio is the degeneracy of the excited multiplet over that
    of the ground state, with Gamma_1 the rate out of one excited state and
    Gamma_2 the total rate in: 2 for a singlet ground state and a doublet
    excitation (spin 1/2, screened side), 1/2 for a doublet ground state.
    current, didv and peak_heights use g_ratio = 1 when Gamma2 is None, so
    pass Gamma2 = detailed_balance(Gamma1, eps, T, g_ratio) for the others."""
    if T == 0.0:
        return 0.0
    return g_ratio*Gamma1*np.exp(-eps/T)


def gamma_tip_from_GN(GN, nu_s=0.5, exact=True):
    """Tip coupling g = 2 pi t^2 nu_t from the normal-state conductance GN
    in units of 2e^2/h, nu_s the substrate density of states per spin
    (1/(2D) for the flat band of the Wilson chain, 0.5 at D = 1), with
    x = pi^2 nu_t nu_s t^2, so that g = 2x/(pi nu_s).

    exact=False is the weak-coupling limit GN = 4x = 2 pi nu_s g, the only
    form that is safe: the rate model itself holds only at GN << 1.
    exact=True (default) inverts GN = 4x/(1+x), the paper's S89 as printed,
    so that its figures can be compared with on the same axis; this form is
    not bounded by one spin-degenerate channel (it gives 2 at x = 1 and
    accepts GN up to 4), so it is not a better form than the weak-coupling
    one. exact="landauer" inverts the transmission of two flat bands
    coupled by t, GN = 4x/(1+x)^2, which is at most 1 (at x = 1; the
    branch x <= 1 is taken). In x, the paper's form is GN/4 above the
    weak-coupling one and the Landauer form about GN/2 above it (0.25% and
    0.5% at GN = 1e-2), and the two exact forms differ by 2.7% at
    GN = 0.1."""
    GN = np.asarray(GN, dtype=float)
    if exact == "landauer":
        if np.any(GN < 0) or np.any(GN > 1):
            raise ValueError("GN must be in [0, 1] in units of 2e^2/h for "
                             "one channel")
        # 4x/(1+x)^2 = GN on x <= 1, written without a cancellation at small GN
        x = GN/(1.0 + np.sqrt(1.0 - GN))**2
        return 2.0*x/(np.pi*nu_s)
    if exact not in (True, False):
        raise ValueError("exact must be True, False or 'landauer'")
    if np.any(GN < 0) or np.any(GN >= 4):
        raise ValueError("GN must be in [0, 4) in units of 2e^2/h")
    if not exact:
        return GN/(2*np.pi*nu_s)
    x = GN/(4.0 - GN)
    return 2.0*x/(np.pi*nu_s)


_GL = np.polynomial.legendre.leggauss(10)


def _graded_nodes(points, lo, hi, hmin):
    """Nodes and weights on [lo, hi], with Gauss-Legendre panels that halve
    towards each of the points and towards lo and hi, down to hmin."""
    pts = np.unique(np.clip(np.asarray(points, dtype=float), lo, hi))
    pts = np.unique(np.concatenate([[lo], pts, [hi]]))
    edges = []
    for a, b in zip(pts[:-1], pts[1:]):
        L = b - a
        if L <= 0.0:
            continue
        half = 0.5*L
        s = [0.0]
        h = hmin
        while h < half:
            s.append(h)
            h = 2.0*h
        s = np.array(s + [half])
        left = a + s
        right = b - s[::-1]
        edges.append(np.concatenate([left, right[1:]]))
    e = np.unique(np.concatenate(edges))
    x0, w0 = _GL
    mid = 0.5*(e[1:] + e[:-1])
    hw = 0.5*(e[1:] - e[:-1])
    nodes = (mid[:, None] + hw[:, None]*x0[None, :]).ravel()
    weights = (hw[:, None]*w0[None, :]).ravel()
    return nodes, weights


def _scales(eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2, dynes):
    """Smallest and typical widths of the problem."""
    gmax = gamma_tip*max(u2, v2)
    omega_t = 0.5*(gmax*np.sqrt(max(Delta_t, 1e-300)))**(2.0/3.0)
    wid = max(Gamma1 + Gamma2 + gmax, dynes, omega_t)
    small = [x for x in (Gamma1 + Gamma2, dynes, gmax, omega_t) if x > 0.0]
    return (min(small) if small else 1e-12*max(Delta_t, abs(eps), 1.0)), wid


def _integrands(w, V, eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2, T,
                dynes):
    Ge = gamma_tip*u2*tip_dos(w - V, Delta_t, dynes)
    Gh = gamma_tip*v2*tip_dos(w + V, Delta_t, dynes)
    G = Ge + Gh + Gamma1 + Gamma2
    nm, np_ = fermi(w - V, T), fermi(w + V, T)
    den = (w - eps)**2 + 0.25*G*G
    s = (Gamma1*(Ge*nm - Gh*np_) - Gamma2*(Ge*(1.0 - nm) - Gh*(1.0 - np_)))/den
    a = 2.0*Ge*Gh*(nm - np_)/den
    return s, a


def current(Vs, eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2=None, T=0.0,
            dynes=0.0, window=None, hmin=None):
    """Single-particle and Andreev currents I_s(V), I_a(V) of the module
    docstring at the biases Vs (energies, eV), in units of G0 = 2e^2/h
    times energy/e. Returns (Is, Ia).

    eps > 0 is the YSR energy, u2 and v2 its electron and hole weights,
    Delta_t and dynes the gap and the Dynes broadening of the tip
    (Delta_t = 0 is a normal tip), gamma_tip = 2 pi t^2 nu_t the tip
    coupling, so that gamma_tip u2 is the paper's gamma_e, Gamma1 and
    Gamma2 the rates that empty and fill the level (Gamma2 None: detailed
    balance, Gamma1 exp(-eps/T)), T the temperature of the tip and of the
    relaxation. window is the half-width W of the frequency window
    (default 10 (Delta_t + eps + max|V|)), hmin the smallest panel of the
    graded mesh (default 1e-4 of the smallest width)."""
    eps, u2, v2 = float(eps), float(u2), float(v2)
    if eps < 0.0:
        raise ValueError("eps must be >= 0; a level at -eps with weights "
                         "(u2, v2) is the level at +eps with (v2, u2)")
    if u2 < 0.0 or v2 < 0.0 or Gamma1 < 0.0 or gamma_tip < 0.0:
        raise ValueError("u2, v2, Gamma1 and gamma_tip must be >= 0")
    if Gamma2 is None:
        Gamma2 = detailed_balance(Gamma1, eps, T)
    Vs = np.atleast_1d(np.asarray(Vs, dtype=float))
    if window is None:
        window = 10.0*(Delta_t + eps + np.max(np.abs(Vs)))
    small, _ = _scales(eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2, dynes)
    if hmin is None:
        hmin = 1e-4*small
    Is, Ia = np.zeros(len(Vs)), np.zeros(len(Vs))
    for k, V in enumerate(Vs):
        pts = [eps, V - Delta_t, V + Delta_t, -V - Delta_t, -V + Delta_t, V, -V]
        w, wt = _graded_nodes(pts, -window, window, hmin)
        s, a = _integrands(w, V, eps, u2, v2, Delta_t, gamma_tip, Gamma1,
                           Gamma2, T, dynes)
        # I = (e/h) int dw (...), and I/G0 = (1/2) int dw (...) in energy/e
        Is[k] = 0.5*np.dot(wt, s)
        Ia[k] = 0.5*np.dot(wt, a)
    return Is, Ia


def didv(Vs, eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2=None, T=0.0,
         dynes=0.0, h=None, window=None, hmin=None):
    """dI_s/dV and dI_a/dV at the biases Vs, in units of G0 = 2e^2/h, by a
    central difference of current() with step h (default 0.02 of the
    width of the peak, the largest of Gamma, the Dynes gamma and the
    tunneling scale (gamma_tip u2 sqrt(Delta_t))^(2/3)/2). Returns
    (dIs, dIa); the total is their sum."""
    if Gamma2 is None:
        Gamma2 = detailed_balance(Gamma1, eps, T)
    Vs = np.atleast_1d(np.asarray(Vs, dtype=float))
    _, wid = _scales(eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2, dynes)
    if h is None:
        h = 0.02*wid
    if window is None:
        window = 10.0*(Delta_t + eps + np.max(np.abs(Vs)) + h)
    kw = dict(T=T, dynes=dynes, window=window, hmin=hmin)
    s1, a1 = current(Vs + 0.5*h, eps, u2, v2, Delta_t, gamma_tip, Gamma1,
                     Gamma2, **kw)
    s0, a0 = current(Vs - 0.5*h, eps, u2, v2, Delta_t, gamma_tip, Gamma1,
                     Gamma2, **kw)
    return (s1 - s0)/h, (a1 - a0)/h


THRESHOLDS = {"alpha+": (1, 1), "alpha-": (-1, 1), "beta+": (1, -1),
              "beta-": (-1, -1)}


def peak_heights(eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2=None, T=0.0,
                 dynes=0.0, peaks=("alpha+", "alpha-", "beta+", "beta-"),
                 npts=61):
    """Heights of the dI/dV peaks of the paper: the main peaks alpha_+-
    at eV = +-(Delta_t + eps) and the thermal peaks beta_+- at eV =
    +-(Delta_t - eps). Each is the maximum of the total dI/dV on a local
    scan of the bias around its threshold (the maximum sits off the
    threshold by up to the width of the peak, which changes the height
    by up to about 35%, S III F), with a second scan ten times finer near
    the threshold and three refinements around the best point, since with
    a sharp tip edge the maximum is a cusp within Gamma_1/2 of it.

    Returns a dict name -> dict(V, total, s, a, threshold), with s and a
    the single-particle and Andreev parts at the maximum and threshold the
    total at the threshold bias itself."""
    if Gamma2 is None:
        Gamma2 = detailed_balance(Gamma1, eps, T)
    _, wid = _scales(eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2, dynes)
    out = {}
    kw = dict(T=T, dynes=dynes)
    args = (eps, u2, v2, Delta_t, gamma_tip, Gamma1, Gamma2)
    for name in peaks:
        sign, branch = THRESHOLDS[name]
        thr = sign*(Delta_t + branch*eps)
        span = 10.0*wid
        if eps > 0.0: # stay away from the other peak of the same sign
            span = min(span, 0.9*eps)
        x = np.linspace(-1.0, 1.0, npts)
        V = np.unique(np.concatenate([thr + span*x, thr + 0.03*span*x]))
        for _ in range(3): # the maximum can be a cusp at the threshold
            ds, da = didv(V, *args, **kw)
            k = int(np.argmax(ds + da))
            Vbest = V[k]
            dV = max(V[min(k+1, len(V)-1)] - V[k], V[k] - V[max(k-1, 0)])
            V = Vbest + np.linspace(-dV, dV, 21)
        ds, da = didv(np.array([Vbest, thr]), *args, **kw)
        out[name] = dict(V=Vbest, total=ds[0] + da[0], s=ds[0], a=da[0],
                         threshold=ds[1] + da[1])
    return out


def rate_parameters(poles, gap, Lambda=None, normalize=True, rtol=1e-3,
                    E_ysr=None):
    """eps, u2, v2 of the rate model from the discrete poles of a
    many-body chain, a dict with "omega" and "weight" as returned by
    zaverage.stm_poles (electron poles at +omega, hole poles at -omega)
    or one entry of zaverage_stm_spectrum's info["subgap"].

    The YSR pole is the subgap pole (|omega| < |gap|) of largest weight,
    eps = |omega| of it; u2 is the sum of the electron weights (omega > 0)
    within rtol*|gap| of +eps and v2 the sum of the hole weights (omega < 0)
    within rtol*|gap| of -eps. The sums matter: stm_poles has one pole per
    tip spin, and the two members of a doublet come out of TASK as two
    poles. With normalize, the weights are divided by the mean total
    weight per tunneling operator, sum(all weights)/2, so that they are per
    unit of the normal-state conductance through the same vertex (1 for the
    bare electron, module docstring). This is for the output of stm_poles
    with its two operators, one per tip spin (the default vertices, or
    [Cdagup[s], Cdagdn[s]]); a different number of operators needs
    normalize=False and a divisor of your own. It needs the full set of
    poles, not only the subgap ones (info["subgap"] of
    zaverage_stm_spectrum takes normalize=False), and raises otherwise.
    With Lambda, both are multiplied by A_Lambda (ysr.lambda_factor), which
    turns the weight of the chain into the residue of the continuum. With
    E_ysr (for instance ysr.subgap_energies(fc)["E_ysr"]), it raises when
    eps differs from |E_ysr| by more than rtol*|gap|.

    A RuntimeWarning says when the numbers are not those of one converged
    YSR state: when the gathered poles are spread by more than 1e-6 of the
    gap, or poles["elastic"] is above 1e-6 of the total weight (both mean
    that DMRG stopped at a mixture of the two parities near the
    singlet-doublet degeneracy, and eps is then off: take the poles with
    zaverage.stm_poles(..., parity=True), or raise nsweeps and compare
    with subgap_energies), and when the subgap poles left out
    carry more than 1e-3 of u2 + v2 (another subgap state, which the
    one-level model does not take, or a doublet split by more than
    rtol*|gap| by the same mixture). Measured on an 8-site chain near the
    transition (Lambda = 4, J = 0.6 A_Lambda, Delta = 0.1, U = 0.3, exact
    eps = 0.01903 Delta, u2 = 0.3338): DMRG runs that stopped at a mixture
    gave the doublet at 0.01765 and 0.01833 Delta (summed here, u2 = 0.3215,
    with the split warning) or at 0.0121 and 0.0156 Delta (beyond rtol, u2
    halved, with the warning on the weight left out); a tolerance of 1e-6
    of the gap, as before, halved u2 in both without a warning.

    Returns a dict with eps, u2, v2, the raw weights We, Wh, norm (the
    divisor), spread (max - min of |omega| over the gathered poles with
    weight above 1e-6 of the largest) and
    others, the subgap poles left out (omega, weight)."""
    import warnings
    from ..ysr import lambda_factor
    om = np.asarray(poles["omega"], dtype=float)
    wt = np.real(np.asarray(poles["weight"]))
    gap = abs(float(gap))
    sub = np.abs(om) < gap
    if not np.any(sub):
        raise ValueError("rate_parameters: no pole inside the gap")
    k = np.flatnonzero(sub)[np.argmax(wt[sub])]
    eps = abs(om[k])
    tol = rtol*gap
    if E_ysr is not None and abs(eps - abs(float(E_ysr))) > tol:
        raise ValueError("rate_parameters: the subgap pole %.6g is not "
            "|E_ysr| = %.6g; the ground state handed to TASK is not converged"
            % (eps, abs(float(E_ysr))))
    e = sub & (om > 0.0) & (np.abs(om - eps) <= tol)
    hl = sub & (om < 0.0) & (np.abs(om + eps) <= tol)
    We, Wh = float(np.sum(wt[e])), float(np.sum(wt[hl]))
    g = e | hl # the spread of the poles that carry weight
    got = np.abs(om[g][np.abs(wt[g]) > 1e-6*np.max(np.abs(wt[g]))])
    spread = float(np.max(got) - np.min(got))
    norm = 1.0
    if normalize:
        if np.all(sub):
            raise ValueError("rate_parameters: normalize needs every pole of "
                "stm_poles, the continuum included; pass normalize=False for "
                "the subgap poles alone")
        norm = 0.5*float(np.sum(wt))
    rest = sub & ~e & ~hl
    elastic = float(poles.get("elastic", 0.0)) if hasattr(poles, "get") else 0.0
    if spread > 1e-6*gap or elastic > 1e-6*max(np.sum(np.abs(wt)), 1e-300):
        warnings.warn("rate_parameters: the YSR poles are split by %.3g "
            "(elastic weight %.3g), so the ground state is a mixture of the "
            "two parities and eps = %.6g is not the YSR energy; take the "
            "poles with zaverage.stm_poles(..., parity=True), or raise "
            "nsweeps and compare with ysr.subgap_energies" % (spread, elastic, eps),
            RuntimeWarning, stacklevel=2)
    left = float(np.sum(np.abs(wt[rest])))
    if left > 1e-3*(abs(We) + abs(Wh)):
        warnings.warn("rate_parameters: subgap poles away from +-eps carry "
            "weight %.3g against %.3g at +-eps: either another subgap state, "
            "which the one-level rate model leaves out, or a YSR doublet split "
            "by a ground state of mixed parity (then take the poles with "
            "parity=True, or raise nsweeps and compare eps with "
            "ysr.subgap_energies)" % (left, We + Wh),
            RuntimeWarning, stacklevel=2)
    fac = lambda_factor(Lambda) if Lambda is not None else 1.0
    return dict(eps=eps, u2=fac*We/norm, v2=fac*Wh/norm, We=We, Wh=Wh,
                norm=norm, spread=spread,
                others=dict(omega=om[rest], weight=wt[rest]))
