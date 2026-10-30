"""Wilson chains: a conduction band mapped to a semi-infinite tight-binding
chain, and the Anderson impurity coupled to its first site.

The band, with density of states rho(e) per spin on [emin, emax], is cut at
the points +-D Lambda^-m (logarithmic discretization, Lambda > 1). Each
interval keeps one level, at the mean energy of the interval, with weight
the integral of rho over it. The Lanczos tridiagonalization of that star,
started from the normalized combination of all levels, gives the chain

    H_bath = sum_n eps_n f_n^dag f_n + sum_n t_n (f_n^dag f_{n+1} + h.c.)

with t_n ~ Lambda^-n/2: site n resolves the energy scale D Lambda^-n/2, so
a chain of n sites reaches temperatures of about D Lambda^-n/2. With rho
normalized to one state, an impurity level coupled by t_D to every band
state (hybridization Gamma = pi t_D^2 rho) couples by exactly t_D to f_0.

The tridiagonalization is done in double precision with full
reorthogonalization. For the flat band it reproduces Wilson's closed form
to 1e-15 relative, down to t_n ~ 1e-12 (n = 50 at Lambda = 3), so the
precision of the coefficients is not what limits the chain length; the
bond dimension the ground state needs is.
"""

import warnings

import numpy as np

from . import fermionchain


def flat_band_hoppings(n, Lambda=2.0, D=1.0):
    """Hoppings t_0..t_{n-2} of the Wilson chain of a flat band on [-D,D],
    in closed form (Wilson, Rev. Mod. Phys. 47, 773, eq. 2.13 without the
    A_Lambda factor). eps_n is zero by particle-hole symmetry."""
    m = np.arange(n-1, dtype=float)
    L = float(Lambda)
    return (D*(1.0+1.0/L)/2.0*(1.0-L**(-m-1))*L**(-m/2.0)
            /np.sqrt((1.0-L**(-2*m-1))*(1.0-L**(-2*m-3))))


def star_levels(Lambda=2.0, band=(-1.0, 1.0), dos=None, nint=60):
    """Levels xi and weights w of the logarithmically discretized band.

    dos is the density of states per spin, a function of energy (flat if
    None); it is normalized here, so only its shape matters. nint intervals
    on each side of the Fermi level (e=0), plus the remainder [0, e Lambda^-nint].
    A one-sided band, emin = 0 < emax, is a bosonic bath: then only the
    positive side is discretized."""
    emin, emax = float(band[0]), float(band[1])
    if not emin <= 0.0 < emax:
        raise ValueError("the band must contain the Fermi level, emin < 0 < emax "
                         "(or emin = 0 < emax for a bosonic bath)")
    if Lambda <= 1.0:
        raise ValueError("Lambda must be larger than 1")
    xs, ws = [], []
    for edge in (emax, emin):
        if edge == 0.0:
            continue # one-sided (bosonic) bath
        pts = edge*float(Lambda)**(-np.arange(nint+1, dtype=float))
        pts = np.append(pts, 0.0)
        for a, b in zip(pts[1:], pts[:-1]):
            lo, hi = min(a, b), max(a, b)
            if dos is None:
                w, e1 = hi-lo, (hi*hi-lo*lo)/2.0
            else:
                from scipy.integrate import quad
                w = quad(dos, lo, hi, epsabs=0.0, epsrel=1e-13)[0]
                e1 = quad(lambda e: e*dos(e), lo, hi, epsabs=0.0, epsrel=1e-13)[0]
            if w > 0.0:
                xs.append(e1/w)
                ws.append(w)
    w = np.array(ws)
    if np.any(w < 0.0):
        raise ValueError("the density of states must be non-negative")
    return np.array(xs), w/np.sum(w)


def wilson_chain(n, Lambda=2.0, band=(-1.0, 1.0), dos=None, nint=None):
    """Onsite energies eps (n values) and hoppings t (n-1 values) of the
    Wilson chain of n sites, for a band on [emin, emax] with density of
    states dos (flat if None, normalized here)."""
    if n < 1:
        raise ValueError("a Wilson chain needs at least one site")
    if nint is None:
        nint = n + 20 # site n probes the intervals around m ~ n/2
    xi, w = star_levels(Lambda=Lambda, band=band, dos=dos, nint=nint)
    if 2*len(xi) < n:
        raise ValueError("not enough intervals for %d sites, raise nint" % n)
    V = np.zeros((len(xi), n))
    V[:, 0] = np.sqrt(w)
    eps, t = np.zeros(n), np.zeros(max(n-1, 0))
    for k in range(n):
        u = xi*V[:, k]
        eps[k] = np.dot(V[:, k], u)
        if k == n-1:
            break
        u = u - eps[k]*V[:, k]
        if k > 0:
            u = u - t[k-1]*V[:, k-1]
        for _ in range(2): # full reorthogonalization, twice
            u = u - V[:, :k+1] @ (V[:, :k+1].T @ u)
        t[k] = np.linalg.norm(u)
        V[:, k+1] = u/t[k]
    return eps, t


def anderson_chain(n, ed=-0.5, U=1.0, tD=0.2, Lambda=2.0, D=1.0, B=0.0,
                   dos=None, kind="standard", **kwargs):
    """Anderson impurity with one or two levels, each coupled by its own
    tD to the first site of its own Wilson chain of n sites, as a
    Spinful_Fermionic_Chain.

    H = sum_l ed_l N_l + U sum_{a<b} n_a n_b + B sum_l Sz_l
        + sum_{l,sigma} tD_l (d^dag_{l sigma} f_{l,0 sigma} + h.c.)
        + sum_l H_bath,l

    ed and tD are numbers (one level) or lists of two (two levels); the
    repulsion U acts between every pair a<b of the 2L impurity spin-orbitals,
    so at ed = -U/2 the impurity holds one electron in a 2L-fold manifold.
    One level sits on site 0 and its chain on sites 1 to n. Two levels sit
    on sites n and n+1, with the chain of level 0 on sites n-1 down to 0 and
    the chain of level 1 on sites n+2 to 2n+1, so every hopping is
    first-neighbor. The chain stores the layout as fc.impurity (one site
    per level) and fc.channels (the sites of each chain, first Wilson site
    first), which get_kondo_screening reads.

    The band is on [-D,D] with density of states dos (flat if None,
    normalized), so the hybridization of level l is Gamma_l = pi tD_l^2 rho.
    The field B acts on the impurity through B Sz and shifts the bath bands
    rigidly to [-D+B, D+B] for spin up and [-D-B, D-B] for spin down, that
    is eps_n -> eps_n +- B with the zero-field t_n. The rigid shift is exact
    for the bath, but it is not the same discretization as cutting the
    shifted band around the Fermi level, which changes t_n for n above
    roughly log_Lambda(D/B).

    Take n odd. With n even the ground state has an odd number of
    electrons (2n+1 for two levels, n+1 for one level at half filling), so
    it is degenerate, and one member of it reads as an unscreened moment;
    a warning says so. Extra keywords go to the chain.

    kind chooses the chain: "standard" is the Wilson chain above, and
    "C0", "C1", "C2" are the renormalized Wilson chains of the open Wilson
    chain construction (see open_wilson_chain), built with the factor
    Lambda^1/2 per site so that site n has the energy D Lambda^-n/2 on every
    kind. For a band symmetric around the Fermi level the three are the
    same chain, since every shift vanishes."""
    eds, tDs = np.atleast_1d(ed), np.atleast_1d(tD)
    if len(eds) != len(tDs) or len(eds) not in (1, 2):
        raise ValueError("ed and tD must both be numbers or both lists of two")
    if n % 2 == 0:
        warnings.warn("n = %d chain sites is even: the ground state then has "
            "an odd number of electrons and is degenerate, and DMRG returns "
            "one member of it, which reads as an unscreened moment; take n "
            "odd" % n)
    if kind == "standard":
        eps, t = wilson_chain(n, Lambda=Lambda, band=(-D, D), dos=dos)
    else:
        eps, t = renormalized_wilson_chain(n, kind=kind, Lambda=np.sqrt(Lambda),
                                           band=(-D, D), dos=dos)
    if len(eds) == 1:
        imp, channels = [0], [list(range(1, n+1))]
    else:
        imp = [n, n+1]
        channels = [list(range(n-1, -1, -1)), list(range(n+2, 2*n+2))]
    fc = fermionchain.Spinful_Fermionic_Chain(len(imp)*(n+1), **kwargs)
    hop = lambda i, j: (fc.Cdagup[i]*fc.Cup[j] + fc.Cdagdn[i]*fc.Cdn[j]
                        + fc.Cdagup[j]*fc.Cup[i] + fc.Cdagdn[j]*fc.Cdn[i])
    modes = []
    h = 0
    for l, i in enumerate(imp):
        h = h + eds[l]*(fc.Nup[i] + fc.Ndn[i]) + B*fc.Sz[i]
        h = h + tDs[l]*hop(i, channels[l][0])
        modes += [fc.Nup[i], fc.Ndn[i]]
        c = channels[l]
        for k in range(n):
            h = h + (eps[k] + B)*fc.Nup[c[k]] + (eps[k] - B)*fc.Ndn[c[k]]
        for k in range(n-1):
            h = h + t[k]*hop(c[k], c[k+1])
    for a in range(len(modes)): # repulsion between every pair of spin-orbitals
        for b in range(a+1, len(modes)):
            h = h + U*modes[a]*modes[b]
    fc.set_hamiltonian(h)
    fc.impurity = imp
    fc.channels = channels
    return fc


# ---------------------------------------------------------------------------
# Open Wilson chains, Bruognolo et al., PRB 95, 121115(R) (2017),
# arXiv:1611.05291.
#
# The bath correlator G(w) = <<b|b^dag>>_w is written as a continued
# fraction. At step k, G_k = 1/(w - eps_k - Sigma_k), and the self-energy is
# split by energy into a slow part (w inside the window band*Lambda^-(k+1))
# and a fast part (the rest), Sigma_k = tS_k^2 G^S_k + tF_k^2 G^F_k. The slow
# part becomes site k+1, the fast part stays as the own bath F_k of site k.
# The chain eps_k, tS_k together with the baths F_k and the last slow bath
# is exact (the open Wilson chain). The renormalized Wilson chain keeps
# only the static part of every bath it drops, as the onsite shift
# Re Sigma^X_k(0) (eq. 4 of the paper).
#
# Numerically the continuum bath is replaced by a fine star (nsub levels
# per factor Lambda), and every step of the continued fraction is then exact
# algebra on its poles: the poles of Sigma_k are the zeros of G_k, one
# between each pair of neighboring poles of G_k, with weights 1/|G_k'|. The
# zeros are found by bisection relative to the pole on their left, so
# energies down to Lambda^-nint keep their relative precision. The level of
# each interval sits at its harmonic mean, w/int(rho/e), so that the star
# has exactly the static correlator of the continuum,
# Re G(0) = -int A(e)/e de, whenever that integral converges. The chain
# coefficients converge to those of the continuum as 1/nsub (about 1% at
# the default nsub = 8); the static identities hold at any nsub.
# Nothing here depends on the statistics: the same chain describes a
# fermionic band and a bosonic bath (a one-sided band, emin = 0).
# ---------------------------------------------------------------------------


def chain_green(eps, t, omega=0.0):
    """Green's function G_00(omega) = [(omega - h)^-1]_00 of the first site
    of a closed tight-binding chain with onsite energies eps and hoppings t,
    that is the free correlator of the bath operator the impurity couples
    to."""
    eps = np.asarray(eps, dtype=complex)
    t = np.asarray(t, dtype=float)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    e0 = np.zeros(len(eps))
    e0[0] = 1.0
    return np.linalg.solve(omega*np.eye(len(eps)) - h, e0)[0]


def bath_star(Lambda=2.0, band=(-1.0, 1.0), dos=None, nsub=8, nint=40):
    """Fine logarithmic star of a bath: levels xi (sorted) and weights w
    (normalized to 1). The cut points are +-Dmax Lambda^(-k/nsub), with
    Dmax = max(|emin|, emax), so both sides share the same points near the
    Fermi level, down to Dmax Lambda^-nint. Each level sits at the harmonic
    mean of its interval, which keeps int A(e)/e de exact. The innermost
    interval [0, Dmax Lambda^-nint] of each side does too when dos(0) = 0;
    otherwise that integral diverges there and the level sits at the
    mean."""
    emin, emax = float(band[0]), float(band[1])
    if not emin <= 0.0 < emax:
        raise ValueError("the band must contain the Fermi level, emin < 0 < emax "
                         "(or emin = 0 < emax for a bosonic bath)")
    if Lambda <= 1.0:
        raise ValueError("Lambda must be larger than 1")
    if int(nsub) < 1:
        raise ValueError("nsub must be a positive integer")
    from scipy.integrate import quad
    f = (lambda e: 1.0) if dos is None else dos
    gapless = (dos is None) or (dos(0.0) != 0.0)
    dmax = max(-emin, emax)
    cuts = dmax*float(Lambda)**(-np.arange(int(nint)*int(nsub)+1)/float(nsub))
    xs, ws = [], []
    for edge in (emax, emin):
        if edge == 0.0:
            continue # one-sided (bosonic) bath
        a, sg = abs(edge), np.sign(edge)
        pts = np.concatenate([[a], cuts[cuts < a], [0.0]])
        g = lambda u: f(sg*u)
        for hi, lo in zip(pts[:-1], pts[1:]):
            if dos is None:
                w = hi - lo
            else:
                w = quad(g, lo, hi, epsabs=0.0, epsrel=1e-13, limit=200)[0]
            if w < 0.0:
                raise ValueError("the density of states must be non-negative")
            if w == 0.0:
                continue
            if lo > 0.0 or not gapless:
                if dos is None:
                    h = np.log(hi/lo)
                else:
                    h = quad(lambda u: g(u)/u, lo, hi, epsabs=0.0,
                             epsrel=1e-13, limit=200)[0]
                x = w/h
            elif dos is None:
                x = hi/2.0
            else:
                x = quad(lambda u: u*g(u), lo, hi, epsabs=0.0, epsrel=1e-13,
                         limit=200)[0]/w
            xs.append(sg*x)
            ws.append(w)
    x, w = np.array(xs), np.array(ws)
    order = np.argsort(x)
    return x[order], w[order]/np.sum(w)


def _self_energy_poles(x, w, nbisect=64):
    """Poles z and weights r of the self-energy of G(z) = sum_k w_k/(z-x_k),
    with x sorted and every w_k > 0: z_j is the zero of G between x_j and
    x_{j+1}, and r_j = 1/sum_k w_k/(z_j-x_k)^2."""
    if len(x) < 2:
        return np.zeros(0), np.zeros(0)
    dx = x[:-1, None] - x[None, :] # x_j - x_k
    lo, hi = np.zeros(len(x)-1), np.diff(x)
    for _ in range(nbisect): # G goes from +inf to -inf in each interval
        mid = 0.5*(lo + hi)
        pos = np.sum(w[None, :]/(dx + mid[:, None]), axis=1) > 0.0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    d = 0.5*(lo + hi)
    r = 1.0/np.sum(w[None, :]/(dx + d[:, None])**2, axis=1)
    return x[:-1] + d, r


def _pole_green(poles, omega):
    """sum_j r_j/(omega - z_j) for poles = (z, r)."""
    z, r = poles
    return np.sum(r/(omega - z))


class OpenWilsonChain(object):
    """Open Wilson chain of n sites (arXiv:1611.05291, eq. 3).

    Site k has the onsite energy eps[k], the hopping tS[k] to site k+1 (for
    the last site, the coupling to its own slow bath), and the coupling
    tF[k] to its own fast bath F_k. The normalized spectral function of
    bath F_k is a set of poles, fast[k] = (energies, weights), and that of
    the slow bath of the last site is slow = (energies, weights). The
    static shifts are dF[k] = Re Sigma^F_k(0) and dS[k] = Re Sigma^S_k(0).
    The chain with all these baths gives exactly the bath correlator of the
    star it was built from (bath_green)."""

    def __init__(self, Lambda, band, eps, tS, tF, dS, dF, fast, slow):
        self.Lambda, self.band = Lambda, band
        self.eps, self.tS, self.tF = eps, tS, tF
        self.dS, self.dF = dS, dF
        self.fast, self.slow = fast, slow
        self.n = len(eps)

    def renormalized(self, n=None, kind="C2"):
        """Onsite energies and hoppings of the renormalized Wilson chain of
        the first n sites (all if None). "C0" has no shifts, "C1" has the
        fast shifts eps_k + dF[k], "C2" has also the slow shift dS[n-1] on
        the last site (eq. 4). C2 reproduces Re Sigma_bath(0) at any n."""
        n = self.n if n is None else int(n)
        if not 1 <= n <= self.n:
            raise ValueError("n must be between 1 and %d" % self.n)
        if kind not in ("C0", "C1", "C2"):
            raise ValueError("kind must be \"C0\", \"C1\" or \"C2\", got %r"
                             % (kind,))
        eps = np.array(self.eps[:n])
        if kind in ("C1", "C2"):
            eps = eps + self.dF[:n]
        if kind == "C2":
            eps[-1] += self.dS[n-1]
        return eps, np.array(self.tS[:n-1])

    def bath_spectrum(self, k, omega, eta, which="F"):
        """Spectrum Gamma^X_k(omega) = tX_k^2 A^X_k(omega) of the own bath of
        site k ("F"), or of the slow bath of the last site ("S"), with every
        pole broadened into a Lorentzian of half width eta."""
        if which == "F":
            (z, r), t2 = self.fast[k], self.tF[k]**2
        elif which == "S":
            if k != self.n-1:
                raise ValueError("only the last site keeps its slow bath")
            (z, r), t2 = self.slow, self.tS[k]**2
        else:
            raise ValueError("which must be \"F\" or \"S\"")
        om = np.asarray(omega, dtype=float)[..., None]
        return t2*np.sum(r*eta/np.pi/((om - z)**2 + eta**2), axis=-1)

    def bath_green(self, omega):
        """Correlator G_0(omega) of the first site, at complex omega above
        the real axis, with every fast bath and the last slow bath kept
        exactly. It equals the correlator of the star: the open chain is an
        exact rewriting of the bath."""
        g = _pole_green(self.slow, omega)
        for k in range(self.n-1, -1, -1):
            g = 1.0/(omega - self.eps[k]
                     - self.tF[k]**2*_pole_green(self.fast[k], omega)
                     - self.tS[k]**2*g)
        return g


def open_wilson_chain(n, Lambda=2.0, band=(-1.0, 1.0), dos=None, nsub=8,
                      nint=None):
    """Open Wilson chain of n sites for a bath with a spectrum of the shape
    dos on band = [emin, emax] (flat if None, normalized here). The slow
    window shrinks by Lambda at every site (the paper's "version 2",
    omega_S,k = band Lambda^-(k+1)), so site k has the energy about
    max|band| Lambda^-k, for a fermionic band as for a bosonic one (a
    standard fermionic Wilson chain needs two sites per factor Lambda).

    The continuum is represented by bath_star(Lambda, band, dos, nsub, nint),
    nint = n + 15 by default. For a spectrum that vanishes at e = 0 (a
    bosonic bath with Gamma ~ omega^s, a pseudogap) the shifts are exact.
    For a finite density of states at the Fermi level the self-energies
    have a finite imaginary part at 0, the shifts are principal values
    (regularized at the scale of the innermost level), and the renormalized
    chain is an approximation."""
    if n < 1:
        raise ValueError("an open Wilson chain needs at least one site")
    if nint is None:
        nint = n + 15
    if nint < n + 5:
        raise ValueError("nint must exceed n by at least 5, so that the slow "
                         "window of the last site still holds levels")
    x, w = bath_star(Lambda=Lambda, band=band, dos=dos, nsub=nsub, nint=nint)
    gapless = (dos is None) or (dos(0.0) != 0.0)
    eta = 3.0*max(-band[0], band[1])*float(Lambda)**(-nint) if gapless else 0.0
    shift = lambda z, r: -np.sum(r*z/(z*z + eta*eta))
    eps, tS, tF, dS, dF, fast = [], [], [], [], [], []
    for k in range(n):
        eps.append(np.dot(w, x))
        z, r = _self_energy_poles(x, w)
        lo, hi = band[0]*Lambda**(-k-1), band[1]*Lambda**(-k-1)
        s = (z >= lo) & (z <= hi)
        if not np.any(s):
            raise ValueError("no level left in the slow window of site %d, "
                             "raise nint" % k)
        tS.append(np.sqrt(np.sum(r[s])))
        tF.append(np.sqrt(np.sum(r[~s])))
        dS.append(shift(z[s], r[s]))
        dF.append(shift(z[~s], r[~s]))
        fast.append((z[~s], r[~s]/max(np.sum(r[~s]), 1e-300)))
        x, w = z[s], r[s]/np.sum(r[s])
    a = np.array
    return OpenWilsonChain(Lambda, band, a(eps), a(tS), a(tF), a(dS), a(dF),
                           fast, (x, w))


def renormalized_wilson_chain(n, kind="C2", **kwargs):
    """Onsite energies (n values) and hoppings (n-1 values) of the
    renormalized Wilson chain of type kind ("C0", "C1" or "C2") with n
    sites; the keywords go to open_wilson_chain."""
    return open_wilson_chain(n, **kwargs).renormalized(kind=kind)


def boson_bath_chain(n, alpha=0.1, s=0.4, wc=1.0, Lambda=2.0, kind="C2",
                     **kwargs):
    """Coupling t_imp and chain (eps, t) of n sites of the bosonic bath
    Gamma(omega) = t_imp^2 A(omega) = 2 alpha wc^(1-s) omega^s on [0, wc]
    (arXiv:1611.05291, eq. 6). kind is "standard" (logarithmic
    discretization at the mean energy and Lanczos, as wilson_chain) or a
    renormalized chain "C0", "C1", "C2". Extra keywords go to
    open_wilson_chain."""
    if s <= 0.0:
        raise ValueError("the bath exponent s must be positive")
    timp = np.sqrt(2.0*alpha*wc/(1.0 + s))
    dos = lambda e: e**s if e > 0.0 else 0.0
    if kind == "standard":
        eps, t = wilson_chain(n, Lambda=Lambda, band=(0.0, wc), dos=dos)
    else:
        eps, t = renormalized_wilson_chain(n, kind=kind, Lambda=Lambda,
                                           band=(0.0, wc), dos=dos, **kwargs)
    return timp, eps, t


def spin_boson_chain(n, alpha=0.1, s=0.4, Delta=0.1, eps=0.0, wc=1.0,
                     Lambda=2.0, kind="C2", maxnb=8, **kwargs):
    """Spin-boson model on a Wilson chain of n bosonic sites,

    H = -Delta/2 sx + 1/2 sz (eps + t_imp (f_0 + f_0^dag)) + H_chain,

    with the bath of boson_bath_chain (arXiv:1611.05291, eq. 7). It returns
    a Bosonic_Chain whose site 0 is the spin as a two-level site
    (sz = 2N-1, sx = A + A^dag) and whose sites 1 to n are the bath, with
    maxnb levels each. bc.x_imp is sz, the operator the bias eps couples
    to. Extra keywords go to the chain."""
    return _boson_impurity(n, alpha, s, wc, Lambda, kind, maxnb, 2,
                           lambda bc: -Delta/2.0*(bc.A[0] + bc.Adag[0]),
                           lambda bc: 2*bc.N[0] - 1, eps, kwargs)


def oscillator_chain(n, alpha=0.1, s=0.4, Omega=1.0, eps=0.0, wc=1.0,
                     Lambda=2.0, kind="C2", maxnb=8, nimp=None, **kwargs):
    """Dissipative harmonic oscillator on a Wilson chain of n bosonic sites,

    H = Omega a^dag a + 1/2 x (eps + t_imp (f_0 + f_0^dag)) + H_chain,

    with x = a + a^dag and the bath of boson_bath_chain (arXiv:1611.05291,
    eq. 5). The model is quadratic, and chi = -d<x>/d eps is
    1/(Omega + t_imp^2 G_chain(0)) on the chain (chain_green) and
    1/(Omega - 2 alpha/s) for the continuum at wc = 1. Site 0 is the
    oscillator, with nimp levels (maxnb if None); bc.x_imp is x."""
    nimp = maxnb if nimp is None else nimp
    return _boson_impurity(n, alpha, s, wc, Lambda, kind, maxnb, nimp,
                           lambda bc: Omega*bc.N[0],
                           lambda bc: bc.A[0] + bc.Adag[0], eps, kwargs)


def _boson_impurity(n, alpha, s, wc, Lambda, kind, maxnb, nimp, himp, ximp,
                    eps, kwargs):
    from . import bosonchain
    timp, e, t = boson_bath_chain(n, alpha=alpha, s=s, wc=wc, Lambda=Lambda,
                                  kind=kind)
    bc = bosonchain.Bosonic_Chain(n+1, maxnb=[nimp] + [maxnb]*n, **kwargs)
    x = ximp(bc)
    h = himp(bc) + 0.5*eps*x + 0.5*timp*x*(bc.A[1] + bc.Adag[1])
    for k in range(n):
        h = h + e[k]*bc.N[k+1]
    for k in range(n-1):
        h = h + t[k]*(bc.Adag[k+1]*bc.A[k+2] + bc.Adag[k+2]*bc.A[k+1])
    bc.set_hamiltonian(h)
    bc.x_imp = x
    bc.impurity = 0
    bc.timp = timp
    bc.chain = (e, t)
    return bc
