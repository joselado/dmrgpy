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
    on each side of the Fermi level (e=0), plus the remainder [0, e Lambda^-nint]."""
    emin, emax = float(band[0]), float(band[1])
    if not emin < 0.0 < emax:
        raise ValueError("the band must contain the Fermi level, emin < 0 < emax")
    if Lambda <= 1.0:
        raise ValueError("Lambda must be larger than 1")
    xs, ws = [], []
    for edge in (emax, emin):
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
                   dos=None, **kwargs):
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
    a warning says so. Extra keywords go to the chain."""
    eds, tDs = np.atleast_1d(ed), np.atleast_1d(tD)
    if len(eds) != len(tDs) or len(eds) not in (1, 2):
        raise ValueError("ed and tD must both be numbers or both lists of two")
    if n % 2 == 0:
        warnings.warn("n = %d chain sites is even: the ground state then has "
            "an odd number of electrons and is degenerate, and DMRG returns "
            "one member of it, which reads as an unscreened moment; take n "
            "odd" % n)
    eps, t = wilson_chain(n, Lambda=Lambda, band=(-D, D), dos=dos)
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
